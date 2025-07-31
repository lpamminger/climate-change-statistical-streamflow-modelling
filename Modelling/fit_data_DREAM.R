# Fitting best performing streamflow models to catchments


# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, dream, smoof, tictoc, furrr, parallel, truncnorm, sloop)


# Import and prepare data-------------------------------------------------------
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


gauge_information <- readr::read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE
)


best_CO2_non_CO2_per_gauge_CMAES <- read_csv(
  "Modelling/Results/CMAES/best_CO2_non_CO2_per_catchment_CMAES.csv",
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
source("./Functions/DREAM.R")
source("./Functions/objective_function_setup.R")
source("./Functions/result_set.R")




# Identify the best streamflow model for a given catchment ---------------------
best_model_and_parameters_per_gauge <- best_CO2_non_CO2_per_gauge_CMAES |> 
  slice_min(AIC, by = gauge) 
  

best_model_per_gauge <- best_model_and_parameters_per_gauge |> 
  select(gauge, streamflow_model) |> 
  distinct()



## Convert the streamflow_model column from characters to a list of functions ===
### match.fun can call a function using the character string of the function name
best_model_function_per_gauge <- map(
  .x = best_model_per_gauge |> pull(streamflow_model),
  .f = match.fun
)
  


# Build catchment_data_set objects ---------------------------------------------
catchment_data <- purrr::map(
  .x = best_model_per_gauge |> pull(gauge), # gauges index will be the same as best_model_function_per_gauge
  .f = catchment_data_blueprint,
  observed_data = data,
  start_stop_indexes = start_stop_indexes
)




# Build numerical optimiser for best streamflow model and catchments -----------
DREAM_bounds_and_transform_methods <- make_DREAM_bounds_and_transform_methods( 
  # function factory - use CMAES parameters to narrow down CO2 parameters
  all_gauges_best_CMAES_parameters = best_model_and_parameters_per_gauge
  )


numerical_optimisers <- map2(
  .x = catchment_data[1:22],
  .y = best_model_function_per_gauge[1:22], 
  .f = numerical_optimiser_setup_vary_inputs,
  objective_function = constant_sd_objective_function,
  streamflow_transform_method = log_sinh_transform,
  bounds_and_transform_method = DREAM_bounds_and_transform_methods, 
  streamflow_transform_method_offset = 0,
  scale = 100,
  minimise_likelihood = FALSE
)



# Run DREAM --------------------------------------------------------------------
DREAM_controls <- list(
  check_convergence_steps = 1000,
  warm_up_per_chain = 1E3,#1E5,
  burn_in_per_chain = 1E3,#3E4, 
  iterations_after_burn_in_per_chain = 1E3,#3E4, 
  eps = 1E-6, #0.1
  steps = 300, #300
  thinning = 1
)


test_DREAM <- DREAM(
  input = numerical_optimisers[[22]],
  controls = DREAM_controls
)


get_convergence_statistics(test_DREAM)
gg_trace_plot(test_DREAM)
gg_distribution_plot(test_DREAM)
