# Construct streamflow using all permutations of DREAM


# Import libraries -------------------------------------------------------------
pacman::p_load(tidyverse, furrr, arrow)


# Import functions -------------------------------------------------------------
source("./Functions/utility.R")
source("./Functions/streamflow_models.R")
source("./Functions/parameter_transformations.R")
source("./Functions/catchment_data_blueprint.R")
source("./Functions/cmaes_dream_summaries.R")
source("./Functions/objective_functions.R")
source("./Functions/numerical_optimiser_setup.R")
source("./Functions/generic_functions.R")
source("./Functions/DREAM.R")
source("./Functions/objective_function_setup.R")
source("./Functions/result_set.R")
source("./Functions/boxcox_logsinh_transforms.R")

# Import data ------------------------------------------------------------------
start_stop_indexes <- readr::read_csv(
  "./Data/Tidy/start_end_index.csv",
  show_col_types = FALSE
)

data <- readr::read_csv(
  "./Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  show_col_types = FALSE
) |>
  mutate(year = as.integer(year)) |>
  # required for log-sinh. Log-sinh current formulation has asymptote of zero.
  # This means zero flows of ephemeral catchments cannot be transformed
  # add a really small value
  mutate(q_mm = q_mm + .Machine$double.eps^0.5)

best_CO2_non_CO2_per_gauge <- read_csv(
  "./Modelling/Results/CMAES/best_CO2_non_CO2_per_catchment_CMAES.csv",
  show_col_types = FALSE
)

DREAM_sequences_query <- open_dataset(
  sources = "./Modelling/Results/DREAM/Parameter_Sequences"
)


CO2_off_DREAM_sequences_query <- DREAM_sequences_query |>
  mutate(
    parameter_value = if_else(str_detect(parameter, "a3"), 0, parameter_value)
  )
# Approach:
# 1. from CO2_impact_analysis get the relevant catchments (only best models with CO2 with an evidence ratio greater than 100)
# 2. filter sequence.parquet files using relevant catchments
# 3. parquet files contain 5 columns: gauge, chain, n, parameter, parameter_value
# 4. for the relevant catchments get the corresponding streamflow model
# 5. streamflow_model inputs: catchment_data, parameter_set
# 6. combine chain and n into single column to produce unique identify for parameter sets
# 7. pull all unique parameter sets, list of vectors
# 8. map(.x = parameter_set, .f = function, catchment_data) --> produce the similar output as cmaes_streamflow_results.csv
# 9. repeat step 8. except set the a3_slope/intercept parameter to zero for CO2 off (change the tibbles directly i.e., mutate(if a3 set to zero))
# 10. repeat CO2_streamflow_impact_analysis to find impact of CO2 on streamflow, use IQR to get uncertainty




# Get relevant catchments and models -------------------------------------------
## Calculate evidence ratio for filtering ======================================
evidence_ratio_calc <- best_CO2_non_CO2_per_gauge |>
  select(gauge, contains_CO2, AIC) |>
  distinct() |>
  pivot_wider(
    names_from = contains_CO2,
    values_from = AIC
  ) |>
  mutate(
    CO2_model = `TRUE`,
    non_CO2_model = `FALSE`,
    .keep = "unused"
  ) |>
  mutate(
    AIC_difference = CO2_model - non_CO2_model # CO2 is smaller than non-CO2 then negative and CO2 is better
  ) |>
  mutate(
    evidence_ratio = case_when(
      AIC_difference < 0 ~ exp(0.5 * abs(AIC_difference)), # when CO2 model is better
      AIC_difference > 0 ~ -exp(0.5 * abs(AIC_difference)) # when non-CO2 model is better
    )
  ) |>
  select(gauge, evidence_ratio) |>
  arrange(evidence_ratio)



## Extract gauges with evidence ratio > 100 and contains CO2 ===================
best_model_per_gauge <- best_CO2_non_CO2_per_gauge |>
  slice_min(
    AIC,
    by = gauge
  ) |>
  left_join(
    evidence_ratio_calc,
    by = join_by(gauge)
  ) |>
  filter(evidence_ratio > 100) |> # I don't need to filter by CO2 model because evidence ratio implicitly does this
  select(gauge, streamflow_model) |>
  distinct()




# Create function to convert sequences_parameters to streamflow ----------------

DREAM_sequences_to_matrix <- function(DREAM_sequences) {
  DREAM_sequences |>
    # combine chain and n to produce unique identifier for each parameter set
    unite(
      col = "identifier",
      chain,
      n
    ) |>
    pivot_wider(
      names_from = parameter,
      values_from = parameter_value
    ) |>
    # remove gauge and identifier cols
    select(!c(gauge, identifier)) |>
    # convert to matrix
    as.matrix() |>
    # transpose
    t() |>
    # strip names
    unname()
  
  # returns a n x m matrix.
  # n = parameters
  # m = unique parameter sets
}




make_DREAM_streamflow_sequences <- function(gauge, observed_data, start_stop_indexes, DREAM_sequences_query, best_model_per_gauge) {
  
  # Make catchment data for streamflow model
  catchment_data <- gauge |>
    catchment_data_blueprint(
      observed_data = data,
      start_stop_indexes = start_stop_indexes
    )


  # Extract sequences for a given gauge
  DREAM_sequences <- DREAM_sequences_query |>
    filter(gauge == {{ gauge }}) |>
    collect()
  
  
  # Get the best streamflow model for a given gauge
  streamflow_model <- best_model_per_gauge |>
    filter(gauge == {{ gauge }}) |>
    pull(streamflow_model) |>
    match.fun()
  
  
  # Convert sequences to a matrix
  parameter_matrix_for_streamflow_model <- DREAM_sequences_to_matrix(DREAM_sequences)
  
  streamflow_model(
    catchment_data = catchment_data,
    parameter_set = parameter_matrix_for_streamflow_model
  ) |> 
    add_column(
      gauge = {{ gauge }},
      .before = 1
    ) |> 
    select(!c(precipitation, observed_streamflow, is_drought_year, CO2, seasonal_ratio))
  
}



# Does gauge 235205 produce NA? No, something is up
# TEST
decade_1 <- seq(from = 1990, to = 1999)
decade_2 <- seq(from = 2012, to = 2021)
  
x <- make_DREAM_streamflow_sequences(
  gauge = "235205",
  observed_data = data, 
  start_stop_indexes = start_stop_indexes, 
  DREAM_sequences_query = DREAM_sequences_query, 
  best_model_per_gauge = best_model_per_gauge
) |> 
  filter(year %in% c(decade_1, decade_2)) |> 
  pivot_longer(
    cols = !c(gauge, year),
    names_to = "permutation",
    values_to = "CO2_on_streamflow"
  ) |> 
  mutate(
    decade = case_when(
      year %in% decade_1 ~ 1,
      year %in% decade_2 ~ 2,
      .default = NA
    )
  ) |> 
  summarise(
    total_decade_CO2_on_streamflow = sum(CO2_on_streamflow),
    .by = c(decade, permutation, gauge)
  ) 


y <- make_DREAM_streamflow_sequences(
  gauge = "235205",
  observed_data = data, 
  start_stop_indexes = start_stop_indexes, 
  DREAM_sequences_query = CO2_off_DREAM_sequences_query, 
  best_model_per_gauge = best_model_per_gauge
) |> 
  filter(year %in% c(decade_1, decade_2)) |> 
  pivot_longer(
    cols = !c(gauge, year),
    names_to = "permutation",
    values_to = "CO2_off_streamflow"
  ) |> 
  mutate(
    decade = case_when(
      year %in% decade_1 ~ 1,
      year %in% decade_2 ~ 2,
      .default = NA
    )
  ) |> 
  summarise(
    total_decade_CO2_off_streamflow = sum(CO2_off_streamflow),
    .by = c(decade, permutation, gauge)
  )

combine_xy <- x |> 
  left_join(
    y,
    by = join_by(gauge, permutation, decade)
  ) |> 
  mutate(
    CO2_percent_impact = ((total_decade_CO2_on_streamflow - total_decade_CO2_off_streamflow) / total_decade_CO2_off_streamflow) * 100
  ) 
  # uncertainty 

# I feel like it would be more impactful if it were +/- the other uncertainty
combine_xy |> summarise(
  IQR_test = IQR(CO2_percent_impact),
  median_test = median(CO2_percent_impact),
  .by = decade
  )
# How does this compare to the CMAES value
# The mean and median is within a percentage point to the CMAES calibrated value
# I could do a total switch for Figure 3 only using DREAM data


# Filter parquet Sequences using best_model_per_gauge for iteration ------------
filter_gauges <- best_model_per_gauge |>
  pull(gauge) |>
  unique()

round_1_gauges <- filter_gauges[1:ceiling(length(filter_gauges)/2)]
round_2_gauges <- filter_gauges[(ceiling(length(filter_gauges)/2) + 1):length(filter_gauges)]

# Iterate of gauges and save ---------------------------------------------------
# furrr and query does not work together - cannot filter in parallel
# This produces a 27 Gb tibble -> split filter gauges in two
# Partitioning into year exceeds RAM
# Likely have to stage it 

# for now split gauges vector in two and get the relevent data for
# figure 3 (e.g., total decade streamflow for the two decades of intrest)

# Temporary solution copy paste

relevant_years <- c(seq(from = 1990, to = 1999), seq(from = 2012, to = 2021))

CO2_off_DREAM_sequences_query <- DREAM_sequences_query |>
  mutate(
    parameter_value = if_else(str_detect(parameter, "a3"), 0, parameter_value)
  )

produce_and_save_decade_streamflow_totals <- function(gauge_set, DREAM_sequences_query, filename) {
  map(
    .x = gauge_set,
    .f = make_DREAM_streamflow_sequences,
    observed_data = data,
    start_stop_indexes = start_stop_indexes,
    DREAM_sequences_query = DREAM_sequences_query, 
    best_model_per_gauge = best_model_per_gauge
  ) |> 
    list_rbind() |> 
    filter(year %in% relevant_years) |> 
    pivot_longer(
      cols = !c(gauge, year),
      names_to = "streamflow_permutation",
      values_to = "streamflow"
    ) |> 
    mutate(
      decade = if_else(year %in% seq(from = 1990, to = 1999), 1, 2)
    ) |> 
    summarise(
      total_decade_streamflow = sum(streamflow),
      .by = c(decade, gauge, streamflow_permutation)
    ) |> 
    write_parquet(
      sink = paste0("./Modelling/Results/DREAM/Streamflow_Sequences/", filename, ".parquet")
    )
}

# query does not play nice with walk/map
# Do this manually to avoid RAM limitations
stop_here() # manually do this
produce_and_save_decade_streamflow_totals(
  gauge_set = round_1_gauges,
  DREAM_sequences_query = DREAM_sequences_query, 
  filename = "streamflow_sequence_CO2_on_1"
)


produce_and_save_decade_streamflow_totals(
  gauge_set = round_2_gauges,
  DREAM_sequences_query = DREAM_sequences_query, 
  filename = "streamflow_sequence_CO2_on_2"
)

produce_and_save_decade_streamflow_totals(
  gauge_set = round_1_gauges,
  DREAM_sequences_query = CO2_off_DREAM_sequences_query, 
  filename = "streamflow_sequence_CO2_off_1"
)

produce_and_save_decade_streamflow_totals(
  gauge_set = round_2_gauges,
  DREAM_sequences_query = CO2_off_DREAM_sequences_query, 
  filename = "streamflow_sequence_CO2_off_2"
)




# Import the streamflow_sequence parquet files and join CO2 on/off -------------
## copy and paste for other sections then rbind
decade_streamflow_CO2_on_sequence_1 <- open_dataset(
  sources = "./Modelling/Results/DREAM/Streamflow_Sequences/streamflow_sequence_CO2_on_1.parquet"
) |> 
  collect() |> 
  rename(total_decade_streamflow_CO2_on = total_decade_streamflow) 

decade_streamflow_CO2_off_sequence_1 <- open_dataset(
  sources = "./Modelling/Results/DREAM/Streamflow_Sequences/streamflow_sequence_CO2_off_1.parquet"
) |> 
  collect() |> 
  rename(total_decade_streamflow_CO2_off = total_decade_streamflow)

decade_streamflow_CO2_on_sequence_2 <- open_dataset(
  sources = "./Modelling/Results/DREAM/Streamflow_Sequences/streamflow_sequence_CO2_on_2.parquet"
) |> 
  collect() |> 
  rename(total_decade_streamflow_CO2_on = total_decade_streamflow)

decade_streamflow_CO2_off_sequence_2 <- open_dataset(
  sources = "./Modelling/Results/DREAM/Streamflow_Sequences/streamflow_sequence_CO2_off_2.parquet"
) |> 
  collect() |> 
  rename(total_decade_streamflow_CO2_off = total_decade_streamflow)

decade_streamflow_CO2_on <- rbind(
  decade_streamflow_CO2_on_sequence_1,
  decade_streamflow_CO2_on_sequence_2
)

decade_streamflow_CO2_off <- rbind(
  decade_streamflow_CO2_off_sequence_1,
  decade_streamflow_CO2_off_sequence_2
)


## left_join CO2_on/off sequences
decade_streamflow_sequences <- decade_streamflow_CO2_on |> 
  left_join(
    decade_streamflow_CO2_off,
    by = join_by(decade, gauge, streamflow_permutation)
  )


# calculate CO2 on/off % difference from streamflow sequences ------------------
# sum_decade_streamflow_CO2_on - sum_decade_streamflow_CO2_off / sum_decade_streamflow_CO2_off,
CO2_impact_on_streamflow_sequences <- decade_streamflow_sequences |> 
  mutate(
    CO2_impact_on_streamflow_percent = (total_decade_streamflow_CO2_on - total_decade_streamflow_CO2_off / total_decade_streamflow_CO2_off) * 100
  ) |> 
  summarise(
    mean_impact_CO2_percentage = mean(CO2_impact_on_streamflow_percent),
    .by = gauge
  )

# This is not good. Something is not right...

