# This no longer works with the new numerical optimisers



# The outcome of this test was a single calibration of CMAES was unable to 
# consitently find the global optima. When run 10 times and taking the lowest
# logliklihood it would get the global optima.





# Explanatory analysis of restarts ---------------------------------------------
## Reasoning in Peterson and Western 2014 - Non-linear ...
cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, cmaesr, smoof, tictoc, truncnorm, sloop, tictoc)
# dream bit more challenging on linux

# Import and prepare data-------------------------------------------------------

## Import annual streamflow, precip, CO2 and gauge data ========================
start_stop_indexes <- read_csv(
  "./Data/Tidy/start_end_index.csv", 
  show_col_types = FALSE
)

data <- read_csv(
  "./Data/Tidy/with_NA_yearly_data_CAMELS.csv", 
  col_types = "icdddddddll", # ensuring each column is a the correct type
  show_col_types = FALSE
)

gauge_information <- read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv", 
  show_col_types = FALSE
)




## Utility functions ===========================================================
source("./Functions/utility.R")
#source("./Functions/boxcox_transforms.R")

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




# Approach ---------------------------------------------------------------------
## 1. pick four catchments to look at
## 2. pick the most complex model combinations
## 3. When tolX is reached restart condition
## 4. Show that the number of restarts reduces uncertainty of the LL or AIC

# 1. Catchments ----------------------------------------------------------------
## pick two difficult catchments that always require restarts and two easy ones. 2 drought and 2 non-drought. Mix of catchments
test_catchments <- c("318076", "608002",  "401009", "104001A")

test_catchment_data <- map(
  .x = test_catchments,
  .f = catchment_data_blueprint,
  observed_data = data,
  start_stop_indexes = start_stop_indexes
)




# 2. For now ignore drought models ---------------------------------------------
# non-drought most complex = streamflow_model_separate_CO2_seasonal_ratio_auto() + CO2_variable_objective_function()
# drought most complex = streamflow_model_drought_separate_CO2_seasonal_ratio_auto() + CO2_variable_objective_function()
forged_test_catchments <- map(
  .x = test_catchment_data,
  .f = numerical_optimiser_setup_vary_inputs,
  streamflow_model = streamflow_model_separate_CO2_seasonal_ratio_auto,
  objective_function = CO2_variable_objective_function,
  bounds_and_transform_method = make_default_bounds_and_transform_methods(),
  scale = 100, # this must be 100 because the control_cmaes is hard coded as 100
  minimise_likelihood = TRUE
)




# 3. make new restart condition for cmaes --------------------------------------
## vary restarts number and conditions
vary_restarts_cmaes_control <- function(restarts_count, parameter_count) {
  
  stopifnot(is.integer(restarts_count))
  
  
  stopOnFlatFitness <- makeStoppingCondition(
    name = "flatFitness",
    message = "Stop if fitness value is flat i.e., fitn vector contains all the same values",
    stop.fun = function(envir){
      return(envir$fitn.ordered[1L] == envir$fitn.ordered[ceiling(0.7 * envir$lambda)])
    },
    code = "flatFitness"
  ) 
  
  
  myStopOnTolX = function(tol = 1E-8) {
    assertNumber(tol, na.ok = FALSE) # this bit is changed
    force(tol)
    return(makeStoppingCondition(
      name = "myTolX",
      message = sprintf("Standard deviation below tolerance in all coordinates."),
      stop.fun = function(envir = parent.frame()) {
        return(all(envir$D < tol) && all((envir$sigma * envir$p.c) < tol))
      }
    ))
    
  }
  
  
  list(
    sigma = 100/3, # hydroState recommends 1/3 of parameter space. # HARD CODED - should extract from optimiser set
    lambda = 4 * floor(3 * log(parameter_count)), # taken from hydroState
    mu = floor((4 * floor(3 * log(parameter_count))) / 2), # taken from hydroState
    restart.triggers = c(
      "flatFitness", 
      "indefCovMat", 
      "noEffectAxis", 
      "noEffectCoord",
      "myTolX"
    ), 
    max.restarts = restarts_count, 
    restart.multiplier = 2L, 
    stop.ons = c(
      list(
        stopOnFlatFitness, 
        myStopOnTolX(tol = 1E-3), 
        stopOnNoEffectAxis(),
        stopOnNoEffectCoord(),
        stopOnCondCov()
      ) 
    )
  ) 
  
}



## Make control_cmaes with different number of restarts ========================
varied_restarts <- map(
  .x = as.integer(seq(from = 0, to = 1, by = 1)), # CHANGE THIS
  .f = vary_restarts_cmaes_control,
  parameter_count = unique(lengths(lapply(forged_test_catchments, "[[", "parameter_names"))) # extracts parameter lengths from each object in forged_test_catchments
)




# Expand grid to get all unique combinations of catchment and restarts ---------
unique_combinations <- expand_grid(
                         forged_test_catchments,
                         varied_restarts
                         )



# Replicate CMAES --------------------------------------------------------------
# Change the code so it runs sequentially rather than in parallel
# The replicate function is cleaner than map with dummy variable
# This is a wrapper of my_cmaes function
repeat_cmaes <- function(numerical_optimiser_setup, cmaes_control, replicate_number) {
  replicate(
    n = replicate_number,
    expr = my_cmaes(
      numerical_optimiser_setup = numerical_optimiser_setup, # do this for all the catchments
      cmaes_control = cmaes_control # do this for all the restart counts
    ),
    simplify = FALSE
  )
}


REPLICATE_NUMBER <- 2

cmaes_restart_results <- map2(
  .x = unique_combinations |> pull(forged_test_catchments),
  .y = unique_combinations |> pull(varied_restarts),
  .f = repeat_cmaes,
  replicate_number = REPLICATE_NUMBER
)




# Turn cma_result objects into result_set objects to extract information -------
summarised_cmaes_restart_results <- map(
  .x = unlist(cmaes_restart_results, recursive = FALSE),
  .f = result_set
  ) |>
  map(.f = restarts_summary) |>
  list_rbind()

#write_csv(summarised_cmaes_restart_results, file = "summarised_cmaes_restarts.csv")
# 9132.13 sec elapsed

# Get error bars (range)
error_bar <- summarised_cmaes_restart_results |> 
  summarise(
  ymin = min(loglikelihood),
  ymax = max(loglikelihood),
  .by = c(gauge, restarts)
)

# Plotting ---------------------------------------------------------------------
## Without exit message ========================================================
plot <- summarised_cmaes_restart_results |>
  ggplot(aes(x = as.factor(restarts), y = loglikelihood)) + # colour = exit_message
  geom_errorbar(
    data = error_bar,
    aes(
      x = as.factor(restarts), 
      ymin = ymin, ymax = ymax
      ), 
    inherit.aes = FALSE,
    width = 0.4
  ) +
  geom_jitter(aes(colour = exit_message)) + #shape = 1
  theme_bw() +
  labs(
    x = "Number of restarts",
    y = "Loglikelihood"
  ) +
  facet_wrap(~gauge, scale = "free") +
  theme(legend.position = "bottom")

plot
                      
#ggsave(
#  filename = paste0("restart_count_analysis_", get_date(), ".pdf"),
#  plot = plot,
#  device = "pdf",
#  width = 297,
#  height = 210,
#  units = "mm"
#)
