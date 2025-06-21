# Fitting streamflow models to catchments

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

# Game plan
# - CMAES should be given chunks to run in parallel
# - Each chunk is a single catchment. There will be 533 chunks.
# - The function must run each replicate for a given model and select the 
#   best performing.
# - Parallel = each chunk in parallel. Each model/replicate for a given catchment
#   in series.
# - Find a way to combine drought and non-drought into a single chunk


# TODO:
# 1. generate catchment_data objects - I currently do everything all at once.
#    
# 2. generate numerical_optimiser objects
#    The repeat/replicates should just repeat the numerical_optimiser objects
#    How do I add drought_models to only drought catchments?
#    Maybe I expand_grid everything (including drought models)
#    add a temporary drought column for gauge with a drought then filter out
#    non-drought gauges with drought models
#
# 3. create function to run through each chunk, select the best run from each 
#    replicate - use parameter_summary |> rbind() to join everything up - 
#    add a replicate identify and use slice_min for likelihood to ge the best 
#    results
#    Name each chunk after a gauge
#    psuedo-code:
#    chunk_and_run_cmaes <- function(chunk_of_numerical_optimisers) {
#      map(.x = chunk_of_numerical_optimisers)
#      parameter_summary() |> rbind() all results of the map
#      slice_max() to get best *
#      * I am not sure this will work by model - likely 
#      * only return a single row
#      * Instead add replicate column - not sure of method - may be a count angle
#      Save filtered parameter_summary
#      map() streamflow_timeseries summary
#      save both as .csv's (easier eventhough .rda uses less space)
#      Avoid saving variables to minimise RAM usage - gc() at end
#      
# 
# }
# There must not be catchments across chunks
# For full utilisation of cores map() the chunk_run_cmaes and future_map
# the cmaes results
# Need to assign gauges to chunk. Name chunks by gauge? Does not work 
# with multiple gauges
# Cannot do it by gauge at a time. This will produce 533 .csv files.
# Join multiple gauges together. 
# Based on previous attempts 5000 was a good amount
# Maximum cmaes per gauge 24 models * 10 replicates = 240
# 5000 / 240 = 20 gauges per chunk - MAKE SURE ONLY ONE GAUGE PER CHUNK


# Idea for adding replicate row
# parameter_summary() will return gauge, streamflow_model, obj, parameter, parameter_value etc
# I could add_column here - make wrapper around parameter_summary
# using imap. Use index to add replicate
# Then rind()

# Constants --------------------------------------------------------------------
# Number of times we want to repeat each catchment-optimiser-streamflow model combinations
REPLICATES <- 2L

# Split catchments for into X chunks (due to RAM limitations).
# I am not convinced a larger chunk size is always better
# The 4500 seemed to run faster than the 9000? Maybe try 4000?
CHUNK_SIZE <- 5000 # items per batch

# I would rather this be gauges_per_chunk


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
#numerical_optimiser_setup(
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
GAUGES_PER_CHUNK <- 1

total_gauges <- gauges_from_combinations |> unique() |> length()

number_of_chunks <- ceiling(total_gauges / GAUGES_PER_CHUNK)

replace_rle_values <- rep(seq(from = 1, to = number_of_chunks, by = 1), each = GAUGES_PER_CHUNK)

if(GAUGES_PER_CHUNK > 1) {
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

run_and_save_cmaes_chunks <- function(numerical_optimisers, chunk_index) {
  
  # Run this in parallel - CMAES
  # This will produce cmaes objects equal to the length of numerical optimisers
  cmaes_results <- map(
    .x = numerical_optimisers,
    .f = CMAES,
    cmaes_control = list(),
    print_monitor = FALSE
  )
  
  # Convert cmaes objects into default format
  default_format_cmaes_results <- map(
    .x = cmaes_results,
    .f = result_set
  )

  # Get parameter_summary of cmaes_results
  cmaes_best_parameters <- map(
    .x = default_format_cmaes_results,
    .f = parameters_summary
  ) |> 
    list_rbind()
  
  
  # Get the best peforming replicate
  best_performing_replicates <- cmaes_best_parameters |> 
    slice_min(
      AIC, # using AIC,
      by = c(gauge, streamflow_model)
    ) |> 
    # if AIC is exactly the same between replicates slice_min will return both
    # I only want a single replicate - use distinct to only present one
    distinct()
  
  write_csv(
    best_performing_replicates,
    file = paste0("Results/CMAES/chunk_", chunk_index, "_cmaes_results_parameter.csv")
  )
}


# FROM HERE
# RUN numerical_optimsers_by_chunk in parallel
# Save results
# do the list_file code below - join all chunks into cmaes_results_parameter 
# and delete small chunk files

# I either add AIC to streamflow_timeseries_summary and repeat what I have done for parameters
# Or I can use the cmaes_results_parameter.csv to generate streamflow
# Option 2 seems the safest

x <- run_and_save_cmaes_chunks(numerical_optimisers_by_chunk[[2]], 1)


# test slice_max if same model
# what happens if AIC is exactly the same
# If AIC is exactly the same slice_min() takes both
z <- x |> 
  filter(streamflow_model == "streamflow_model_precip_only")

z_prime <- z |> 
  mutate(
    AIC = 500,
    parameter_value = 1
  )

z_same <- z |> 
  mutate(
    AIC = 500,
    parameter_value = 1
  )

z_other_model <- x |> 
  filter(streamflow_model == "streamflow_model_precip_auto")

test_zs <- rbind(z_same, z_prime) |> 
  slice_min(
    AIC,
    by = c(gauge, streamflow_model),
    with_ties = TRUE # this breaks everything if changed from TRUE
  ) |> 
  distinct(gauge, streamflow_model) # this seems to work


# Join everything into a single file and save ----------------------------------
## Parameters ==================================================================
## Empty results folder before joining just in case
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



# THIS WILL BREAK IF ANOTHER OBJECTIVE FUNCTION IS ADDED OR REMOVED
parameter_combinations <- map(
  .x = c(non_drought_streamflow_models, drought_streamflow_models),
  .f = get_parameter_number,
  objective_function = all_objective_functions[[1]]
) |> 
  list_rbind()

#parameter_combinations_2 <- map(
#  .x = c(non_drought_streamflow_models, drought_streamflow_models),
#  .f = get_parameter_number,
#  objective_function = all_objective_functions[[2]]
#) |> 
#  list_rbind()

all_parameter_combinations <- parameter_combinations # rbind(parameter_combiations1, 2 etc.)




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
