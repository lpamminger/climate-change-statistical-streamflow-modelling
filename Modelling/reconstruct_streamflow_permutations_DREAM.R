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
  sources = "./Modelling/Results/DREAM/Sequences"
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


total_streamflow_years <- function(x, period_1, period_2) {
  x |>
    filter(year %in% c(period_1, period_2)) |>
    pivot_longer(
      cols = !c(gauge, year),
      names_to = "permutation",
      values_to = "streamflow"
    ) |>
    mutate(
      decade = case_when(
        year %in% decade_1 ~ 1,
        year %in% decade_2 ~ 2,
        .default = NA
      )
    ) |>
    summarise(
      total_decade_streamflow = sum(streamflow),
      .by = c(decade, permutation, gauge)
    )
}


# Create function to summarise the impact of CO2 on streamflow uncertainty -----
CO2_impact_total_decade_streamflow_uncertainty <- function(gauge, period_1, period_2, observed_data, start_stop_indexes, DREAM_sequences_query, best_model_per_gauge) {
  
  # Convert DREAM parameter sequences to DREAM streamflow sequences  
  total_decade_streamflow_CO2_on <- make_DREAM_streamflow_sequences(
    gauge = {{ gauge }},
    observed_data = observed_data,
    start_stop_indexes = start_stop_indexes,
    DREAM_sequences_query = DREAM_sequences_query,
    best_model_per_gauge = best_model_per_gauge
  ) |>
    total_streamflow_years(
      period_1 = period_1,
      period_2 = period_2
    ) |>
    rename(
      total_decade_streamflow_CO2_on = total_decade_streamflow
    )
  
  # Turn off CO2 component in DREAM_sequences_query
  CO2_off_DREAM_sequences_query <- DREAM_sequences_query |>
    mutate(
      parameter_value = if_else(str_detect(parameter, "a3"), 0, parameter_value)
    )
  
  # Convert DREAM parameter sequences to DREAM streamflow sequences with the 
  # CO2 component turned off
  total_decade_streamflow_CO2_off <- make_DREAM_streamflow_sequences(
    gauge = {{ gauge }},
    observed_data = observed_data,
    start_stop_indexes = start_stop_indexes,
    DREAM_sequences_query = CO2_off_DREAM_sequences_query,
    best_model_per_gauge = best_model_per_gauge
  ) |>
    total_streamflow_years(
      period_1 = period_1,
      period_2 = period_2
    ) |>
    rename(
      total_decade_streamflow_CO2_off = total_decade_streamflow
    )
  
  
  # Return summary tibble
  total_decade_streamflow_CO2_on |> 
    left_join(
      total_decade_streamflow_CO2_off,
      by = join_by(decade, permutation, gauge)
    ) |> 
    mutate(
      CO2_impact_on_streamflow_percentage = ((total_decade_streamflow_CO2_on - total_decade_streamflow_CO2_off) / total_decade_streamflow_CO2_off) * 100
    ) |> 
    summarise(
      IQR_CO2_impact_on_streamflow_percentage = IQR(CO2_impact_on_streamflow_percentage),
      mean_CO2_impact_on_streamflow_percentage = mean(CO2_impact_on_streamflow_percentage),
      median_CO2_impact_on_streamflow_percentage = median(CO2_impact_on_streamflow_percentage),
      .by = decade
    ) |> 
    add_column(
      gauge = {{ gauge }},
      .before = 1
    )
}



# Filter parquet Sequences using best_model_per_gauge for iteration ------------
gauges <- best_model_per_gauge |>
  pull(gauge) |>
  unique()

# Calculate and save the impact of CO2 on streamflow uncertainty ---------------
decade_1 <- seq(from = 1990, to = 1999)
decade_2 <- seq(from = 2012, to = 2021)

map(
  .x = gauges,
  .f = CO2_impact_total_decade_streamflow_uncertainty,
  period_1 = decade_1,
  period_2 = decade_2,
  observed_data = data, 
  start_stop_indexes = start_stop_indexes, 
  DREAM_sequences_query= DREAM_sequences_query, 
  best_model_per_gauge = best_model_per_gauge
) |> 
  list_rbind() |> 
  write_csv(
    file = "./Modelling/Results/DREAM/DREAM_CO2_impact_uncertainty_on_streamflow.csv"
  )




