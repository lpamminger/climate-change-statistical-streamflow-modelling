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

source("./Functions/boxcox_transforms.R") # can remove

## Import streamflow functions =================================================
source("./Functions/streamflow_models.R")
source("./Functions/parameter_transformations.R")
source("./Functions/catchment_data_blueprint.R")
source("./Functions/cmaes_dream_summaries.R")
source("./Functions/objective_functions.R")
source("./Functions/numerical_optimiser_setup.R")
source("./Functions/generic_functions.R")
source("./Functions/my_cmaes.R")
source("./Functions/objective_function_setup.R")
source("./Functions/result_set.R")



# Remove testing - I think changes are good. Ready to rock and roll.
gauge <- "315450"

example <- gauge |>
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  ) |>
  numerical_optimiser_setup_vary_inputs(
    streamflow_model = streamflow_model_separate_shifted_CO2,
    objective_function = CO2_variable_objective_function,
    bounds_and_transform_method = make_default_bounds_and_transform_methods(),
    minimise_likelihood = TRUE
  ) |>
  my_cmaes(print_monitor = TRUE) |>
  result_set()

x <- parameters_summary(example)
plot(example)





# Number of times we want to repeat each catchment-optimiser-streamflow model combinations
REPEATS <- 10L

# Split catchments for into X chunks (due to RAM limitations).
# I am not convinced a larger chunk size is always better
# The 4500 seemed to run faster than the 9000? Maybe try 4000?
CHUNK_SIZE <- 4000 # items per batch


# Construct catchment_data objects ---------------------------------------------
### All ########################################################################
all_catchment_data <- purrr::map(
  .x = unique(data$gauge),
  .f = catchment_data_blueprint,
  observed_data = data,
  start_stop_indexes = start_stop_indexes
)

all_catchment_data <- all_catchment_data


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


drought_catchment_data <- drought_catchment_data


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
  .f = run_and_save_chunks_optimiser_parallel, 
  optimiser = my_cmaes,
  save_streamflow = TRUE,
  save_sequences = FALSE,
  is_drought = FALSE
)


gc()


### Drought ####################################################################
plan(multisession, workers = length(availableWorkers())) # set once for furrr
iwalk(
  .x = chunk_repeat_all_drought_numerical_optimisers_cmaes,
  .f = run_and_save_chunks_optimiser_parallel, 
  optimiser = my_cmaes,
  save_streamflow = TRUE,
  save_sequences = FALSE,
  is_drought = TRUE
)


gc()



# Join everything into a single file and save ----------------------------------
## Parameters ==================================================================
parameters_list_of_files <- list.files(
  path = "./Results/my_cmaes/",
  recursive = FALSE, # I don't want it looking in other folders
  pattern = "parameter",
  full.names = TRUE
)



# Issue: sometimes replicates produce >= 2 exact same results
# if LL is the same slice_min saves everything
# I only want a single catchment-model-objective-function combinations


# Work around:
# - add another column with the parameter number
# - for each catchment-model-objective-function group take the first X
#   rows where is the number of parameters
# This removes duplicates for parameters

# For the streamflow duplicates are removed using distinct(year)


get_parameter_number <- function(streamflow_model, objective_function) {
  streamflow_model <- noquote(streamflow_model)
  objective_function <- noquote(objective_function)
  
  streamflow_model_name <- streamflow_model()[[1]]
  objective_function_name <- objective_function()[[1]]
  parameter_number <- length(c(streamflow_model()[[2]], objective_function()[[2]]))
  
  tibble(
    "streamflow_model" = streamflow_model_name,
    "objective_function" = objective_function_name,
    "parameter_number" = parameter_number
  )
}



# THIS WILL BREAK IF ANOTHER OBJECTIVE FUNCTION IS ADDED
parameter_combinations <- map(
  .x = c(all_streamflow_models, drought_streamflow_models),
  .f = get_parameter_number,
  objective_function = all_objective_functions[[1]]
) |> 
  list_rbind()

parameter_combinations_2 <- map(
  .x = c(all_streamflow_models, drought_streamflow_models),
  .f = get_parameter_number,
  objective_function = all_objective_functions[[2]]
) |> 
  list_rbind()

all_parameter_combinations <- rbind(parameter_combinations, parameter_combinations_2)




combined_cmaes_parameters <- parameters_list_of_files |>
  readr::read_csv(show_col_types = FALSE) |> 
  left_join(
    all_parameter_combinations,
    by = join_by(streamflow_model, objective_function)
  ) |>  
  dplyr::slice_min(  
    loglikelihood,
    by = c(gauge, streamflow_model, objective_function), # only get the minimum LL for gauge, streamflow model and objective function combination
    ) |>
  slice(  
    seq(from = 1, to = parameter_number[1]),
    .by = c(gauge, streamflow_model, objective_function)
  ) |> 
  select(!parameter_number) |> 
  readr::write_csv(
    file = paste0("./Results/my_cmaes/CMAES_parameter_results_", get_date(), ".csv")
  )



## Streamflow ==================================================================
streamflow_list_of_files <- list.files(
  path = "./Results/my_cmaes/",
  recursive = FALSE, # I don't want it looking in other folders
  pattern = "streamflow",
  full.names = TRUE
)


combined_cmaes_streamflow <- streamflow_list_of_files |>
  readr::read_csv(show_col_types = FALSE) |>
  dplyr::slice_min( 
    loglikelihood,
    by = c(gauge, streamflow_model, objective_function, optimiser) # only get the minimum LL for gauge, streamflow model and objective function combination
  ) |>
  distinct( # remove log-likelihood duplicates
    year, 
    gauge, 
    streamflow_model, 
    objective_function, 
    optimiser, 
    .keep_all = TRUE
    ) |> 
  select(!loglikelihood) |> 
  readr::write_csv(
    file = paste0("./Results/my_cmaes/CMAES_streamflow_results_", get_date(), ".csv")
  )



## Delete all files ending with .csv ===========================================
delete_files <- list.files(
  path = "./Results/my_cmaes/",
  recursive = FALSE, # I don't want it looking in other folders
  pattern = "chunk", # remove all chunk files
  full.names = TRUE
) |>
  file.remove()
