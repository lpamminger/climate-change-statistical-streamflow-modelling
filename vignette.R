# An example to fit models to annual streamflow data --------------------------- 



# Import libraries required ----------------------------------------------------
# If this is a package this step is not required
pacman::p_load(tidyverse, cmaesr, smoof, truncnorm, sloop, ggExtra)
## ggExtra is for ggMarginal used in the examine_transform plot



# Import and prepare data-------------------------------------------------------
# If this is a package this step is not required
start_stop_indexes <- readr::read_csv(
  "./Data/Tidy/start_end_index.csv",
  show_col_types = FALSE
)

data <- readr::read_csv( # data will be in the package
  "./Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  show_col_types = FALSE
) |>
  mutate(
    year = as.integer(year)
  ) |> 
  # required for log-sinh. Log-sinh current formulation has asymptote of zero. 
  # This means zero flows of ephemeral catchments cannot be transformed
  # add a really small value
  mutate(q_mm = q_mm + .Machine$double.eps^0.5) 

gauge_information <- readr::read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE
)


## Import functions ============================================================
#If this is a package this step is not required
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



# 1. Select gauge to test from data --------------------------------------------

gauge <- "302208"   


# ideally catchment_data_blueprint should have other methods of data entry such as giving vectors individually
# See catchment_data_blue_print helper... must leave observed_data blank


example_catchment <- gauge |> 
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  ) 

# plot to see if there is anything strange happening in data
plot(example_catchment, type = "streamflow-time")
plot(example_catchment, type = "rainfall-runoff")

# 2. Prepare numerical_optimiser object ----------------------------------------

# Here you can adjust:
## - the streamflow_model used. 
###   See get_all_streamflow_models() for a list off all available models.
## - objective_function used. 
###   See get_all_objective_functions() for a list of all available objective functions.
## - streamflow_transform_method used. 
###   See get_streamflow_transform_method() for a list of all available methods
## - bounds_and_transform method used. 
###   Currently only changing the tibble located in make_default_bounds_and_transform_methods 
###   can be used to change bounds and transform methods
## - minimise_likelihood. 
###   A logical argument that can either minimise or maximise negative loglikelihood (dependent of numerical optimiser - CMAES or DREAM)
## - streamflow_transform_method_offset
###   Shifts the asymptote of streamflow_transform_method



numerical_optimiser <- example_catchment |>
  numerical_optimiser_setup_vary_inputs(
    streamflow_model = streamflow_model_slope_shifted_CO2_seasonal_ratio_auto,
    objective_function = constant_sd_objective_function, 
    streamflow_transform_method = log_sinh_transform, 
    bounds_and_transform_method = make_default_bounds_and_transform_methods, 
    minimise_likelihood = TRUE,
    streamflow_transform_method_offset = 0 
  )




# 3. Put numerical_optimiser object into a numerical optimiser -----------------

# currently only works with CMAES and DREAM functions
results <- numerical_optimiser |>
  CMAES(
    cmaes_control = list(), # can alter default cmaes controls here - see cmaesr package for options
    print_monitor = TRUE # print optimiser outputs in console while CMAES is running
  ) 


# 4. Convert results into a result_set object (standardised format) ------------
standardised_results <- results |> result_set() 


# 5. Examine results -----------------------------------------------------------
### probably should add ... in plot method so the user can alter the ggplot object directly i.e.,
### plot(standardised_results, type = "streamflow-time, ...) add theme(axis...)
plot(standardised_results, type = "streamflow-time")
plot(standardised_results, type = "rainfall-runoff")
plot(standardised_results, type = "examine_transform") 
plot(standardised_results, type = "std_errors")


# The near_bound column does not functioning correctly
# Ideally, it should scale with the parameters. It does not do this currently.
parameter_table <- parameters_summary(standardised_results)
head(parameter_table)
result_table <- streamflow_timeseries_summary(standardised_results) 

# DREAM is a bit more complex...



