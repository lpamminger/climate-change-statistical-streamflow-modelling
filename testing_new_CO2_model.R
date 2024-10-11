# Testing/Implementing new model
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
source("./Functions/objective_function_setup.R")
source("./Functions/result_set.R")


# Playing with CO2 data
sample_CO2_data <- data |> 
  filter(gauge == "102101A") |> 
  pull(CO2)


a3 <- 0.05
a4 <- 0.75 * max(sample_CO2_data) # a negative value shifts the impact of CO2 earlier in the streamflow model. Do I want this?

# a max(sample_CO2_data) effectively removes the CO2 components

base_CO2_equation <- a3 * sample_CO2_data
matrix_test <- matrix(rep(sample_CO2_data - a4, times = 5), nrow = length(sample_CO2_data))
matrix_test[matrix_test <= 0] <- 0 
new_CO2_equation <- a3 * matrix_test



plot(ts(base_CO2_equation), ylim = range(c(base_CO2_equation, new_CO2_equation)), col = "blue")
lines(ts(new_CO2_equation[, 1]), col = "red")
lines(x = c(0, 59), y = c(0, 0), lty = 2)



# Steps to implementing new model ----------------------------------------------
# a5 is the CO2 shift variable
# I have hard coded the a5 parameter in make_parameter_bounds - this should require the user to input Co2 data to work out bounds.



# Test function ----------------------------------------------------------------
test_streamflow_model_separate_CO2 <- function(catchment_data, parameter_set) { 
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "test_streamflow_model_separate_CO2",
        "parameters" = c("a0", "a1", "a3", "a5")
      )
    )
  } 
  
  # Check if parameter is a matrix - if not coerce
  if(!is.matrix(parameter_set)) {
    parameter_set <- as.matrix(parameter_set, ncol = 1)
  }
  
  # Parameters
  a0 <- parameter_set[1, ]
  a1 <- parameter_set[2, ]
  a3 <- parameter_set[3, ]
  a5 <- parameter_set[4, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  CO2 <- catchment_data$CO2
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) 
  repeat_CO2 <- matrix(CO2, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) 
  repeat_a0 <- matrix(a0, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a3 <- matrix(a3, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a5 <- matrix(a5, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Model
  repeat_shifted_CO2 <- repeat_CO2 - repeat_a5
  repeat_shifted_CO2[repeat_shifted_CO2 < 0] <- 0 # max(0, CO2 - a5)
  repeat_a0 + (repeat_a1 * repeat_precipitation) + (repeat_a3 * repeat_shifted_CO2)
  
}



# Run CMAES --------------------------------------------------------------------
gauge <- "407214"

cmaes_example <- gauge |> 
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  ) |> 
  numerical_optimiser_setup_vary_inputs(
    streamflow_model = test_streamflow_model_separate_CO2,
    objective_function = constant_sd_objective_function,
    bounds_and_transform_method = make_default_bounds_and_transform_methods(),
    minimise_likelihood = TRUE
  ) |> 
  my_cmaes(print_monitor = TRUE) |> 
  result_set() 

cmaes_parameters <- cmaes_example |> 
  parameters_summary()

cmaes_example |> 
  plot()
