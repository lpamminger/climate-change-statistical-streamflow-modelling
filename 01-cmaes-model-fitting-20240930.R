# Fitting streamflow models to catchments
cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, cmaesr, smoof, tictoc, furrr, parallel, truncnorm, sloop, tictoc)

# Import and prepare data-------------------------------------------------------

## Import annual streamflow, precip, CO2 and gauge data ========================
start_stop_indexes <- readr::read_csv(
  "./Data/Tidy/start_end_index.csv",
  show_col_types = FALSE
)

data <- readr::read_csv(
  "./Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  col_types = "icdddddddll", # ensuring each column is a the correct type
  show_col_types = FALSE
)

gauge_information <- readr::read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE
)



## Utility functions ===========================================================
source("./Functions/utility.R")


## Import streamflow functions =================================================
source("./Functions/streamflow_models.R")
source("./Functions/parameter_transformations.R")
source("./Functions/catchment_data_blueprint.R")
source("./Functions/cmaes_dream_summaries.R")
source("./Functions/objective_functions.R")
source("./Functions/numerical_optimiser_setup.R")
source("./Functions/generic_functions.R")
source("./Functions/my_cmaes.R")
source("./Functions/my_dream.R")
source("./Functions/objective_function_setup.R")
source("./Functions/result_set.R")




# Number of times we want to repeat each catchment-optimiser-streamflow model combinations
REPEATS <- 10L 

# Split catchments for into X chunks (due to RAM limitations).
# Must be a multiple of REPEATS to avoid duplication across chunks. There is a check in code just in case
CHUNK_SIZE <- 4500 # items per batch


# Construct catchment_data objects ---------------------------------------------
### All ########################################################################
all_catchment_data <- purrr::map(
  .x = unique(data$gauge),
  .f = catchment_data_blueprint,
  observed_data = data,
  start_stop_indexes = start_stop_indexes
)



### Drought ####################################################################
drought_gauges <- gauge_information |>
  dplyr::filter(drought == TRUE) |>
  dplyr::pull(gauge)


drought_catchment_data <- map(
  .x = drought_gauges,
  .f = catchment_data_blueprint,
  observed_data = data,
  start_stop_indexes = start_stop_indexes
)





# Build objective_functions using the optimiser_set object ---------------------
all_streamflow_models <- get_non_drought_streamflow_models() 
drought_streamflow_models <- get_drought_streamflow_models() 
all_objective_functions <- get_all_objective_functions()


## Produce a tibble of all combinations of catchment, optimiser, etc ===========
ready_for_iteration <- tidyr::expand_grid(
  all_catchment_data,
  all_streamflow_models,
  all_objective_functions
) |>
  dplyr::distinct() # make sure there are no double ups


### Drought ####################################################################
drought_ready_for_iteration <- tidyr::expand_grid(
  drought_catchment_data,
  drought_streamflow_models,
  all_objective_functions
) |>
  dplyr::distinct()




## Make numerical_optimiser_setup using ready_for_iteration tibble =============
all_numerical_optimisers_cmaes <- pmap(
  .l = list(
    ready_for_iteration |> dplyr::pull(all_streamflow_models),
    ready_for_iteration |> dplyr::pull(all_objective_functions),
    ready_for_iteration |> dplyr::pull(all_catchment_data)
  ),
  .f = numerical_optimiser_setup,
  bounds_and_transform_method = make_default_bounds_and_transform_methods(),
  minimise_likelihood = TRUE # cmaes = TRUE, dream = FALSE
)


repeat_all_numerical_optimisers_cmaes <- rep(all_numerical_optimisers_cmaes, each = REPEATS)


chunk_repeat_all_numerical_optimisers_cmaes <- split( # required for chunking in parallel
  repeat_all_numerical_optimisers_cmaes,
  ceiling(seq_along(repeat_all_numerical_optimisers_cmaes) / CHUNK_SIZE)
)



### Drought ####################################################################
drought_numerical_optimisers_cmaes <- pmap(
  .l = list(
    drought_ready_for_iteration |> dplyr::pull(drought_streamflow_models),
    drought_ready_for_iteration |> dplyr::pull(all_objective_functions),
    drought_ready_for_iteration |> dplyr::pull(drought_catchment_data)
  ),
  .f = numerical_optimiser_setup,
  bounds_and_transform_method = make_default_bounds_and_transform_methods(),
  minimise_likelihood = TRUE # cmaes = TRUE, dream = FALSE
)

repeat_drought_numerical_optimisers_cmaes <- rep(drought_numerical_optimisers_cmaes, each = REPEATS)


chunk_repeat_all_drought_numerical_optimisers_cmaes <- split( # required for chunking in parallel
  repeat_drought_numerical_optimisers_cmaes,
  ceiling(seq_along(repeat_drought_numerical_optimisers_cmaes) / CHUNK_SIZE)
)


# Run cmaes --------------------------------------------------------------------


plan(multisession, workers = length(availableWorkers())) # set once for furrr
iwalk(
  .x = chunk_repeat_all_numerical_optimisers_cmaes,
  .f = run_and_save_chunks_my_cmaes_parallel,
  is_drought = FALSE
)


gc()


### Drought ####################################################################
plan(multisession, workers = length(availableWorkers())) # set once for furrr
iwalk(
  .x = chunk_repeat_all_drought_numerical_optimisers_cmaes,
  .f = run_and_save_chunks_my_cmaes_parallel,
  is_drought = TRUE
)





# Join everything into a single file and save ----------------------------------
## Parameters ==================================================================
parameters_list_of_files <- list.files(
  path = "./Results/CMAES_results/",
  recursive = FALSE, # I don't want it looking in other folders
  pattern = "parameter",
  full.names = TRUE
)


combined_cmaes_parameters <- parameters_list_of_files |>
  readr::read_csv(show_col_types = FALSE) |>
  dplyr::slice_min( # check if REPEATS go across batches
    loglikelihood,
    by = c(gauge, streamflow_model, objective_function, optimiser)
  ) |>
  readr::write_csv(
    file = paste0("./Results/CMAES_results/CMAES_parameter_results.csv")
  )


## Streamflow ==================================================================
streamflow_list_of_files <- list.files(
  path = "./Results/CMAES_results/",
  recursive = FALSE, # I don't want it looking in other folders
  pattern = "streamflow",
  full.names = TRUE
)


combined_cmaes_streamflow <- streamflow_list_of_files |>
  readr::read_csv(show_col_types = FALSE) |>
  dplyr::slice_min( # I need this because if the REPEATS are spread over multiple chunk it will return > 1 best parameter set
    loglikelihood,
    by = c(gauge, streamflow_model, objective_function, optimiser)
  ) |>
  select(!loglikelihood) |>
  readr::write_csv(
    file = paste0("./Results/CMAES_results/CMAES_streamflow_results.csv")
  )


## Delete all files ending with .csv ===========================================
delete_files <- list.files(
  path = "./Results/CMAES_results/",
  recursive = FALSE, # I don't want it looking in other folders
  pattern = "chunk", # remove all chunk files
  full.names = TRUE
) |>
  file.remove()
