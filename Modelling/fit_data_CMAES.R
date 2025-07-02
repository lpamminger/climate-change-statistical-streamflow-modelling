# Fitting streamflow models to catchments


# THERE IS SOMETHING WRONG WITH THE a parameter being equal to 150
# The models components are wildly different
# READ THE LAST PARAGRAPH OF QJ'S PAPER
# This tells you the possible range and says that a and b should use a 
# log-transform - we are using a log-transform :)
# https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2011WR010973
# In article QJ cautions the use of high values as they start to 
# cause linear shifts - I think
# instead of chasing the "perfect bounds" we stick with the recommended

# bounds are [-15, 0] in realspace exp(-15) = 3E-7 to 1

# In this paper https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/wrcr.20445
# Range of parameters is different - b is adjusted based on mean streamflow
# I could try the b_hat angle? Why?

# from this paper: https://hess.copernicus.org/articles/24/2017/2020/
# a << y << b
# see what values they used

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, cmaesr, smoof, tictoc, furrr, parallel, truncnorm, sloop)

# Import and prepare data-------------------------------------------------------
start_stop_indexes <- readr::read_csv(
  "./Data/Tidy/start_end_index.csv",
  show_col_types = FALSE
)

data <- readr::read_csv(
  "./Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  show_col_types = FALSE
) |>
  mutate(year = as.integer(year))

gauge_information <- readr::read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE
)


# Functions --------------------------------------------------------------------
source("./Functions/utility.R")
source("./Functions/boxcox_logsinh_transforms.R")
source("./Functions/streamflow_models.R")
source("./Functions/parameter_transformations.R")
source("./Functions/catchment_data_blueprint.R")
source("./Functions/cmaes_dream_summaries.R")
source("./Functions/objective_functions.R")
source("./Functions/numerical_optimiser_setup.R")
source("./Functions/generic_functions.R")
source("./Functions/CMAES.R")
source("./Functions/objective_function_setup.R")
source("./Functions/result_set.R")


# Constants --------------------------------------------------------------------
# Number of times we want to repeat each catchment-optimiser-streamflow model combinations
REPLICATES <- 10L 


# Construct catchment_data objects ---------------------------------------------
catchment_data <- purrr::map(
  .x = unique(data$gauge),
  .f = catchment_data_blueprint,
  observed_data = data,
  start_stop_indexes = start_stop_indexes
)



# Build numericial_optimiser_set object ----------------------------------------
streamflow_models <- get_all_streamflow_models()


## Produce a tibble of all combinations of catchment, optimiser, etc ===========
catchment_streamflow_model_combinations <- tidyr::expand_grid(
  catchment_data,
  streamflow_models
)

## Currently, drought models are assigned to non-drought gauges. Remove these.
## Method:
### 1. Identify gauges that should have drought streamflow models (drought_gauge)
### 2. Filter out gauges that does not have a drought and has a drought streamflow
###    model

drought_gauge <- gauge_information |>
  select(gauge, drought)

get_gauge <- function(catchment_data_object) {
  catchment_data_object$gauge_ID
}

gauge_drought_for_combinations <- map_chr(
  .x = catchment_streamflow_model_combinations |> pull(catchment_data),
  .f = get_gauge
)

is_drought_streamflow_model <- function(streamflow_model) {
  streamflow_model <- noquote(streamflow_model)
  streamflow_model()$name |>
    str_detect("drought")
}

streamflow_model_drought_for_combinations <- map_lgl(
  .x = catchment_streamflow_model_combinations |> pull(streamflow_models),
  .f = is_drought_streamflow_model
)

catchment_streamflow_model_combinations <- catchment_streamflow_model_combinations |>
  add_column(
    "gauge" = gauge_drought_for_combinations,
    .before = 1
  ) |>
  add_column(
    "drought_streamflow_model" = streamflow_model_drought_for_combinations
  ) |>
  left_join(
    drought_gauge,
    by = join_by(gauge)
  ) |>
  rename(
    observed_rainfall_drought = drought
  )

# if drought_streamflow_model = TRUE and drought = FALSE -> remove
remove_combinations <- catchment_streamflow_model_combinations |>
  filter((drought_streamflow_model & !observed_rainfall_drought))

catchment_streamflow_model_combinations <- catchment_streamflow_model_combinations |>
  anti_join(
    remove_combinations,
    by = join_by(gauge, catchment_data, streamflow_models, drought_streamflow_model, observed_rainfall_drought)
  )



## Add objective_function, streamflow_transform_method and bounds to tibble ====
objective_functions <- get_all_objective_functions() # there is only one objective function at this time

streamflow_transform_methods <- list(log_sinh_transform) # we only want to use the log-sinh transform method

extracted_catchment_data <- catchment_streamflow_model_combinations |>
  dplyr::pull(catchment_data)

bounds <- map(
  .x = extracted_catchment_data,
  .f = make_default_bounds_and_transform_methods
)

catchment_streamflow_model_combinations <- catchment_streamflow_model_combinations |>
  add_column(
    "objective_functions" = objective_functions,
    "streamflow_transform_methods" = streamflow_transform_methods,
    "bounds" = bounds
  )


## Add replicates by repeating catchment_streamflow_model_combinations =========
repeat_catchment_streamflow_model_combinations <- do.call(
  "rbind",
  replicate(n = REPLICATES, catchment_streamflow_model_combinations, simplify = FALSE)
) |>
  arrange(gauge)


## Make numerical_optimiser_setup catchment_streamflow_model_combinations  =====

# This must be in the same order as numerical_optimiser inputs
# numerical_optimiser_setup(
#  streamflow_model,
#  objective_function,
#  streamflow_transform_method,
#  catchment_data,
#  bounds_and_transform_method,
#  streamflow_transform_method_offset,
#  scale = 100,
#  minimise_likelihood = TRUE
#  )

pmap_list <- list(
  repeat_catchment_streamflow_model_combinations |> pull(streamflow_models),
  repeat_catchment_streamflow_model_combinations |> pull(objective_functions),
  repeat_catchment_streamflow_model_combinations |> pull(streamflow_transform_methods),
  repeat_catchment_streamflow_model_combinations |> pull(catchment_data),
  repeat_catchment_streamflow_model_combinations |> pull(bounds)
)



numerical_optimisers <- pmap( # This is slow
  .l = pmap_list,
  .f = numerical_optimiser_setup,
  streamflow_transform_method_offset = 0,
  scale = 100,
  minimise_likelihood = TRUE
)


## Repeat numerical optimisers based on the replicates required ================

### There must not be any gauges across chunks
### To avoid this combine gauges into chunks

gauges_from_combinations <- repeat_catchment_streamflow_model_combinations |> pull(gauge)

run_length_gauges_from_combinations <- rle(gauges_from_combinations)


# Change values in rle ($values) to 1, 2, 3 etc. to construct factor values for split
# The rle $values must be split based on GAUGES_PER_CHUNK
GAUGES_PER_CHUNK <- 24L

total_gauges <- gauges_from_combinations |>
  unique() |>
  length()

number_of_chunks <- ceiling(total_gauges / GAUGES_PER_CHUNK)

replace_rle_values <- rep(seq(from = 1, to = number_of_chunks, by = 1), each = GAUGES_PER_CHUNK)

if (GAUGES_PER_CHUNK > 1) {
  remove_excess_replace_rle_values <- length(replace_rle_values) - total_gauges

  remove_excess_index <- length(replace_rle_values) - seq(from = 1, to = remove_excess_replace_rle_values, by = 1)

  replace_rle_values <- replace_rle_values[-remove_excess_index]
}


new_rle <- run_length_gauges_from_combinations

new_rle$values <- replace_rle_values

factor_sequences <- inverse.rle(new_rle)


# Chunk numerical optimisers by gauge
numerical_optimisers_by_chunk <- split(numerical_optimisers, f = factor_sequences)







# Run cmaes --------------------------------------------------------------------

## Wrapper around cmaes to add replicates ======================================
replicate_CMAES <- function(numerical_optimiser_setup, replicate_number, cmaes_control, print_monitor) {
  cmaes_result <- CMAES(
    numerical_optimiser_setup = numerical_optimiser_setup,
    cmaes_control = cmaes_control,
    print_monitor = print_monitor
  )

  cmaes_result$replicate_number <- replicate_number
  return(cmaes_result)
}


result_set_replicates <- function(cmaes_or_dream_result) {
  result <- result_set(cmaes_or_dream_result)
  result$replicate_number <- cmaes_or_dream_result$replicate_number
  return(result)
}

parameters_summary_replicates <- function(x) {
  parameters_summary(x) |>
    add_column(
      "replicate" = x$replicate_number,
      .before = 1
    )
}

streamflow_timeseries_summary_replicates <- function(x) {
  streamflow_timeseries_summary(x) |>
    add_column(
      "replicate" = x$replicate_number,
      .before = 1
    )
}



run_and_save_cmaes_chunks <- function(numerical_optimisers, chunk_index) {
  tic()

  # Seq_along numerical optimisers to get replicate_numbers
  replicate_number <- seq_along(numerical_optimisers)

  # Run in parallel
  cmaes_results <- future_map2(
    .x = numerical_optimisers,
    .y = replicate_number,
    .f = replicate_CMAES,
    cmaes_control = list(),
    print_monitor = FALSE,
    .options = furrr_options(seed = 1L),
    .progress = TRUE
  )

  # Convert cmaes objects into default format
  default_format_cmaes_results <- map(
    .x = cmaes_results,
    .f = result_set_replicates
  )

  # Get parameter_summary of cmaes_results
  cmaes_best_parameters <- map(
    .x = default_format_cmaes_results,
    .f = parameters_summary_replicates
  ) |>
    list_rbind()


  # Get the best peforming replicate
  best_performing_replicates <- cmaes_best_parameters |>
    slice_min(
      AIC, # using AIC,
      by = c(gauge, streamflow_model)
    ) |>
    # call slice_min again but on replicate
    # if AIC and parameters are the same between replicates the previous slice
    # min will take all replicates. We only want a single replicate to be
    # saved.
    slice_min(
      replicate,
      by = c(gauge, streamflow_model)
    )


  # Get streamflow_timeseries_summary of cmaes_results
  cmaes_streamflow_timeseries <- map(
    .x = default_format_cmaes_results,
    .f = streamflow_timeseries_summary_replicates
  ) |>
    list_rbind()

  # Use replicate number to get best streamflow
  filter_replicate <- best_performing_replicates |>
    pull(replicate) |>
    unique()


  best_cmaes_streamflow_timeseries <- cmaes_streamflow_timeseries |>
    filter(replicate %in% filter_replicate)

  # Save results
  write_csv(
    best_performing_replicates,
    file = paste0("Modelling/Results/CMAES/chunk_", chunk_index, "_cmaes_results_parameter.csv")
  )

  write_csv(
    best_cmaes_streamflow_timeseries,
    file = paste0("Modelling/Results/CMAES/chunk_", chunk_index, "_cmaes_results_timeseries.csv")
  )

  cat("Chunk", chunk_index, "complete\n")
  toc()
}


## Run in all chunks ===========================================================
plan(multisession, workers = length(availableWorkers())) # parallel bit inside run_and_save_cmaes_chunks
iwalk(
  .x = numerical_optimisers_by_chunk, # after testing set to numerical_optimisers_by_chunk
  .f = run_and_save_cmaes_chunks
)


# Combine chunks into a single file --------------------------------------------
## parameter results ===========================================================
combined_cmaes_parameters_file <- list.files(
  path = "./Modelling/Results/CMAES/",
  recursive = FALSE, # I don't want it looking in other folders
  pattern = "parameter",
  full.names = TRUE
) |>
  readr::read_csv(show_col_types = FALSE)


## streamflow results ==========================================================
combined_cmaes_streamflow_file <- list.files(
  path = "./Modelling/Results/CMAES/",
  recursive = FALSE, # I don't want it looking in other folders
  pattern = "timeseries",
  full.names = TRUE
) |>
  readr::read_csv(show_col_types = FALSE)


## save big files ==============================================================
write_csv(
  x = combined_cmaes_parameters_file,
  file = "Modelling/Results/CMAES/cmaes_parameter_results.csv"
)

write_csv(
  x = combined_cmaes_streamflow_file,
  file = "Modelling/Results/CMAES/cmaes_streamflow_results.csv"
)

## remove chunks ===============================================================
### using file.remove does not put the files into recycling bin - delete perm
### delete manually
