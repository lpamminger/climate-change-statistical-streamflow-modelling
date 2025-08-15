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



# Collect DREAM sequences ------------------------------------------------------
# Insufficent RAM to load the entire dataset
# Split gauges in two. Run script 2 times
gauges <- best_model_per_gauge |>
  pull(gauge) |>
  unique()

round_1_gauges <- gauges[1:41]
round_2_gauges <- gauges[42:81]

DREAM_sequences <- open_dataset(
  sources = "./Modelling/Results/DREAM/Sequences"
) |>
  filter(gauge == "235234") |> 
  #filter(gauge %in% round_1_gauges) |> 
  collect()



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




make_DREAM_streamflow_sequences <- function(gauge, observed_data, start_stop_indexes, DREAM_sequences, best_model_per_gauge) {
  browser()
  # Make catchment data for streamflow model
  catchment_data <- gauge |>
    catchment_data_blueprint(
      observed_data = data,
      start_stop_indexes = start_stop_indexes
    )


  # Extract sequences for a given gauge
  gauge_DREAM_sequences <- DREAM_sequences |>
    filter(gauge == {{ gauge }})


  if (is_empty_tibble(gauge_DREAM_sequences)) {
    stop("gauge not found in DREAM sequences")
  }


  # Get the best streamflow model for a given gauge
  streamflow_model <- best_model_per_gauge |>
    filter(gauge == {{ gauge }}) |>
    pull(streamflow_model) |>
    match.fun()


  # Convert sequences to a matrix
  parameter_matrix_for_streamflow_model <- DREAM_sequences_to_matrix(gauge_DREAM_sequences)

  transformed_streamflow <- streamflow_model(
    catchment_data = catchment_data,
    parameter_set = parameter_matrix_for_streamflow_model
  ) |>
    add_column(
      gauge = {{ gauge }},
      .before = 1
    ) |>
    select(!c(precipitation, observed_streamflow, is_drought_year, CO2, seasonal_ratio))

  # Convert transformed streamflow to realspace streamflow
  transformed_streamflow_permutations <- transformed_streamflow |>
    select(!c(gauge, year)) |>
    as.list()

  realspace_streamflow_permutations <- map2(
    .x = parameter_matrix_for_streamflow_model[nrow(parameter_matrix_for_streamflow_model), ],
    .y = transformed_streamflow_permutations,
    .f = inverse_log_sinh_transform,
    offset = 0
  ) |>
    # need to name columns otherwise as_tibble() throws error
    `names<-`(paste0("perm_", seq(from = 1, to = length(transformed_streamflow_permutations)))) |>
    as_tibble()

  transformed_streamflow |>
    select(gauge, year) |>
    cbind(realspace_streamflow_permutations)
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


CO2_impact_total_decade_streamflow_uncertainty <- function(gauge, period_1, period_2, observed_data, start_stop_indexes, DREAM_sequences, CO2_off_DREAM_sequences, best_model_per_gauge) {
  # Convert DREAM parameter sequences to DREAM streamflow sequences
  total_decade_streamflow_CO2_on <- make_DREAM_streamflow_sequences(
    gauge = {{ gauge }},
    observed_data = observed_data,
    start_stop_indexes = start_stop_indexes,
    DREAM_sequences = DREAM_sequences,
    best_model_per_gauge = best_model_per_gauge
  ) |>
    total_streamflow_years(
      period_1 = period_1,
      period_2 = period_2
    ) |>
    rename(
      total_decade_streamflow_CO2_on = total_decade_streamflow
    )


  # Convert DREAM parameter sequences to DREAM streamflow sequences with the
  # CO2 component turned off
  total_decade_streamflow_CO2_off <- make_DREAM_streamflow_sequences(
    gauge = {{ gauge }},
    observed_data = observed_data,
    start_stop_indexes = start_stop_indexes,
    DREAM_sequences = CO2_off_DREAM_sequences,
    best_model_per_gauge = best_model_per_gauge
  ) |>
    total_streamflow_years(
      period_1 = period_1,
      period_2 = period_2
    ) |>
    rename(
      total_decade_streamflow_CO2_off = total_decade_streamflow
    )

  browser()
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


# Turn off CO2 component in DREAM_sequences ------------------------------------
CO2_off_DREAM_sequences <- DREAM_sequences |>
  mutate(
    parameter_value = if_else(str_detect(parameter, "a3"), 0, parameter_value)
  )



# Calculate and save the impact of CO2 on streamflow uncertainty ---------------
decade_1 <- seq(from = 1990, to = 1999)
decade_2 <- seq(from = 2012, to = 2021)


CO2_impact_total_decade_streamflow_uncertainty(
  gauge = "235234",
  period_1 = decade_1,
  period_2 = decade_2,
  observed_data = data,
  start_stop_indexes = start_stop_indexes,
  DREAM_sequences = DREAM_sequences,
  CO2_off_DREAM_sequences = CO2_off_DREAM_sequences,
  best_model_per_gauge = best_model_per_gauge
)

# Infs are present when back-transforming - solution?

DREAM_CO2_impact_uncertainty_on_streamflow <- map(
  .x = gauges,
  .f = CO2_impact_total_decade_streamflow_uncertainty,
  period_1 = decade_1,
  period_2 = decade_2,
  observed_data = data,
  start_stop_indexes = start_stop_indexes,
  DREAM_sequences = DREAM_sequences,
  CO2_off_DREAM_sequences = CO2_off_DREAM_sequences,
  best_model_per_gauge = best_model_per_gauge
) |>
  list_rbind()


DREAM_CO2_impact_uncertainty_on_streamflow |>
  write_csv(
    file = "./Modelling/Results/DREAM/DREAM_CO2_impact_uncertainty_on_streamflow_round_1.csv"
  )


