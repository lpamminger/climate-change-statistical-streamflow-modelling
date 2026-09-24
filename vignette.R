# An example to fit models to annual streamflow data --------------------------- 



# Import libraries required ----------------------------------------------------
pacman::p_load(tidyverse, cmaesr, smoof, truncnorm, sloop, ggExtra, tictoc, dream) # Once this becomes a package this step is not required


# There are multiple packages called dream.
# Install the following:
# - install.packages("dream", repos="http://R-Forge.R-project.org") 
# OR
# - https://rdrr.io/rforge/dream/



# Expected run time ------------------------------------------------------------
# The expected run time for the CMAES part of the vignette is less than 30 seconds
# The expected run time of the DREAM part of the vignette is approximately 15 minutes




# Import and prepare data-------------------------------------------------------
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
source("./Functions/DREAM.R")



# 1. Select gauge to test from data --------------------------------------------

gauge <- "405245"   


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




# CMAES EXAMPLE ----------------------------------------------------------------
# 2. Prepare numerical_optimiser object ========================================

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
    streamflow_model = streamflow_model_intercept_shifted_CO2,
    objective_function = constant_sd_objective_function, 
    streamflow_transform_method = log_sinh_transform, 
    bounds_and_transform_method = make_default_bounds_and_transform_methods, 
    minimise_likelihood = TRUE,
    streamflow_transform_method_offset = 0 
  )




# 3. Put numerical_optimiser object into a numerical optimiser =================
# currently only works with CMAES and DREAM functions
results <- numerical_optimiser |>
  CMAES(
    cmaes_control = list(), # can alter default cmaes controls here - see cmaesr package for options
    print_monitor = TRUE # print optimiser outputs in console while CMAES is running
  ) 


# 4. Convert results into a result_set object (standardised format) ============
standardised_results <- results |> result_set() 


# 5. Examine results ===========================================================
### probably should add ... in plot method so the user can alter the ggplot object directly i.e.,
### plot(standardised_results, type = "streamflow-time, ...) add theme(axis...)
plot(standardised_results, type = "streamflow-time")
plot(standardised_results, type = "rainfall-runoff")
plot(standardised_results, type = "examine_transform") 
plot(standardised_results, type = "std_errors")


### Data in table form
parameter_table <- parameters_summary(standardised_results) # The near_bound column does not functioning correctly. Ideally, it should scale with the parameters. It does not do this currently.
result_table <- streamflow_timeseries_summary(standardised_results) 





# DREAM EXAMPLE (repeat steps 2-5) ---------------------------------------------
# 2. Prepare numerical_optimiser object ========================================
DREAM_numerical_optimiser <- example_catchment |>
  numerical_optimiser_setup_vary_inputs(
    streamflow_model = streamflow_model_intercept_shifted_CO2,
    objective_function = constant_sd_objective_function, 
    streamflow_transform_method = log_sinh_transform, 
    bounds_and_transform_method = make_default_bounds_and_transform_methods, 
    minimise_likelihood = FALSE, # DREAM requires maximising likelihood
    streamflow_transform_method_offset = 0 
  )


# 3. Put numerical_optimiser object into a numerical optimiser =================
DREAM_results <- DREAM_numerical_optimiser |>
  DREAM(
    controls = list() # uses parameters see DREAM_controls below. Passing empty list uses the defaults.
  ) 


# Default controls for DREAM. These were found to work best for this analysis
DREAM_controls <- list(
  check_convergence_steps = 200,
  warm_up_per_chain = 5E4, 
  burn_in_per_chain = 1E4, 
  iterations_after_burn_in_per_chain = 1E4, 
  eps = 0.05, 
  steps = 10, 
  thinning = 1
)

# For more complex controls examine the following control_DREAM_defaults()
# or `dream` package directly





# 4. Convert results into a result_set object (standardised format) ============
DREAM_standardised_results <- DREAM_results |> result_set() 


# 5. Examine results ===========================================================
plot(DREAM_standardised_results, type = "streamflow-time")
plot(DREAM_standardised_results, type = "rainfall-runoff")
plot(DREAM_standardised_results, type = "examine_transform") 
plot(DREAM_standardised_results, type = "std_errors")

### DREAM also has specific plots:
gg_trace_plot(DREAM_results) # trace plot
gg_distribution_plot(DREAM_results) # distribution plot


### Data in table form
DREAM_parameter_table <- parameters_summary(DREAM_standardised_results)
DREAM_result_table <- streamflow_timeseries_summary(DREAM_standardised_results) 










