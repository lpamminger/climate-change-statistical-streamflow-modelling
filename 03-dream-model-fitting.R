## Fitting streamflow models to catchments
cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, dream, coda, tictoc, furrr, parallel, truncnorm, sloop, tictoc)
# install dream using: install.packages("dream", repos="http://R-Forge.R-project.org")
# install.packages("coda")

# floor(RAM / 3Gb) = CHUNK_SIZE


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

CMAES_results <- read_csv("./Results/my_cmaes/CMAES_parameter_results_20241108.csv",
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
source("./Functions/my_dream.R")
source("./Functions/objective_function_setup.R")
source("./Functions/result_set.R")






# Identify best model combination for DREAM gauges -----------------------------
best_model_combination_per_catchment <- CMAES_results |>
  select(!c(parameter, parameter_value, optimiser, loglikelihood, exit_message, near_bounds)) |>
  slice_min(
    AIC,
    by = gauge
  ) |> 
  distinct() |> 
  select(!AIC) |> 
  rename(
    streamflow_model_name = streamflow_model,
    objective_function_name = objective_function
  ) 


## Mannually replace these catchments ==========================================
## Non-utilised variable is performing better than models with all parameters used

# 137101A REPLACE streamflow_model_separate_shifted_CO2_seasonal_ratio CO2_variable_objective_function
# 215207 REPLACE streamflow_model_drought_precip_seasonal_ratio_auto constant_sd_objective_function
# 216004 REPLACE streamflow_model_separate_shifted_CO2_auto CO2_variable_objective_function
# 415237 REPLACE streamflow_model_precip_seasonal_ratio_auto CO2_variable_objective_function
# 138004B streamflow_model_precip_only CO2_variable_objective_function
# A5130501 streamflow_model_precip_seasonal_ratio_auto constant_sd_objective_function

hard_code_change <- tribble(
  ~gauge,      ~streamflow_model_name,                                  ~objective_function_name,
  "137101A",   "streamflow_model_separate_shifted_CO2_seasonal_ratio",  "CO2_variable_objective_function",
  "215207",    "streamflow_model_drought_precip_seasonal_ratio_auto",   "constant_sd_objective_function",
  "216004",    "streamflow_model_separate_shifted_CO2_auto",            "CO2_variable_objective_function",
  "415237",    "streamflow_model_precip_seasonal_ratio_auto",           "CO2_variable_objective_function",
  "138004B",   "streamflow_model_precip_only",                          "CO2_variable_objective_function",
  "A5130501",  "streamflow_model_precip_seasonal_ratio_auto",           "constant_sd_objective_function"
)


best_model_combination_per_catchment <- best_model_combination_per_catchment |> 
  filter(!gauge %in% hard_code_change$gauge) |> 
  rbind(hard_code_change)


## Use a join to assign each row the correct streamflow and objective function =====
get_model_name_and_function <- function(model) {
  model_name <- model()$name
  tibble(
    "model_name" = model_name,
    "model" = list(model)
  )
}


streamflow_name_and_model_for_joining <- map(
  .x = c(get_non_drought_streamflow_models(), get_drought_streamflow_models()),
  .f = get_model_name_and_function
) |>
  list_rbind() |> 
  rename(
    streamflow_model_name = model_name,
    streamflow_model = model
  )


objective_function_name_and_model_for_join <- map(
  .x = get_all_objective_functions(),
  .f = get_model_name_and_function
) |> 
  list_rbind() |> 
  rename(
    objective_function_name = model_name,
    objective_function = model
  )



## Left join ===================================================================
gauge_streamflow_model_objective_function_combinations <- best_model_combination_per_catchment |> 
  left_join(
    streamflow_name_and_model_for_joining,
    by = join_by(streamflow_model_name)
  ) |> 
  left_join(
    objective_function_name_and_model_for_join,
    by = join_by(objective_function_name)
  ) |> 
  select(!c(streamflow_model_name, objective_function_name)) |> 
  unclass() |> 
  unname() # names breaks pmap? I don't know why - remove them






# Produce numerical_optimiser_objects ready for the optimiser ------------------
numerical_optimiser_objects <- function(gauge, streamflow_model, objective_function, data, start_stop) {
  gauge |>
    catchment_data_blueprint(
      observed_data = data,
      start_stop_indexes = start_stop
    ) |>
    numerical_optimiser_setup_vary_inputs(
      streamflow_model = streamflow_model,
      objective_function = objective_function,
      bounds_and_transform_method = make_default_bounds_and_transform_methods(),
      minimise_likelihood = FALSE
    )
}



ready_for_optimisation <- pmap(
  .l = gauge_streamflow_model_objective_function_combinations,
  .f = numerical_optimiser_objects,
  data = data,
  start_stop = start_stop_indexes
)






# Split ready_for_optimisation objects into chunks to avoid exceeding RAM ------
# Subject to change. I need to work out RAM requirements
CHUNK_SIZE <- 28

chunked_ready_for_optimisation <- split( # required for chunking in parallel
  ready_for_optimisation,
  ceiling(seq_along(ready_for_optimisation) / CHUNK_SIZE)
)




# Run the optimiser (calibration) ----------------------------------------------

## Non-drought optimising ======================================================
plan(multisession, workers = length(availableWorkers())) # set once for furrr
iwalk(
  .x = chunked_ready_for_optimisation,
  .f = run_and_save_chunks_optimiser_parallel,
  optimiser = my_dream,
  save_streamflow = FALSE,
  save_sequences = TRUE,
  is_drought = FALSE
)




# Join everything into a single file and save ----------------------------------
## Parameters ==================================================================
parameters_list_of_files <- list.files(
  path = "./Results/my_dream/",
  recursive = FALSE, # I don't want it looking in other folders
  pattern = "parameter",
  full.names = TRUE
)



parameters_list_of_files |>
  readr::read_csv(show_col_types = FALSE) |>
  dplyr::slice_min(
    loglikelihood,
    by = c(gauge, streamflow_model, objective_function) # only get the minimum LL for gauge, streamflow model and objective function combination
  ) |>
  readr::write_csv(
    file = paste0("./Results/my_dream/DREAM_parameter_results.csv")
  )


## Sequences ===================================================================
sequences_list_of_files <- list.files(
  path = "./Results/my_dream/",
  recursive = FALSE, # I don't want it looking in other folders
  pattern = "sequences",
  full.names = TRUE
)


sequences_list_of_files |>
  readr::read_csv(show_col_types = FALSE) |>
  readr::write_csv(
    file = paste0("./Results/my_dream/DREAM_sequence_results.csv")
  )


## Delete all files ending with .csv ===========================================
 delete_files <- list.files(
  path = "./Results/my_dream/",
  recursive = FALSE, # I don't want it looking in other folders
  pattern = "chunk", # remove all chunk files
  full.names = TRUE
 ) |>
  file.remove()










