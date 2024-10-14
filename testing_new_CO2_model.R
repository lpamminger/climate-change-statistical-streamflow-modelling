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

# 1.
streamflow_model_separate_shifted_CO2 <- function(catchment_data, parameter_set) { 
  
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


# 2.
streamflow_model_separate_shifted_CO2_seasonal_ratio <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_separate_shifted_CO2_seasonal_ratio",
        "parameters" = c("a0", "a1", "a3", "a4", "a5")
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
  a4 <- parameter_set[4, ]
  a5 <- parameter_set[5, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  seasonal_ratio <- catchment_data$seasonal_ratio
  CO2 <- catchment_data$CO2
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) #matrix(rep(precipitation, times = ncol(parameter_set)), ncol = ncol(parameter_set), byrow = FALSE)
  repeat_seasonal_ratio <- matrix(seasonal_ratio, ncol = ncol(parameter_set), nrow = length(seasonal_ratio), byrow = FALSE)
  repeat_CO2 <- matrix(CO2, ncol = ncol(parameter_set), nrow = length(CO2), byrow = FALSE) 
  
  repeat_a0 <- matrix(a0, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a3 <- matrix(a3, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a4 <- matrix(a4, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a5 <- matrix(a5, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Model
  repeat_shifted_CO2 <- repeat_CO2 - repeat_a5
  repeat_shifted_CO2[repeat_shifted_CO2 < 0] <- 0 # max(0, CO2 - a5)
  repeat_a0 + (repeat_a1 * repeat_precipitation) + (repeat_a3 * repeat_shifted_CO2) + (repeat_a4 * repeat_seasonal_ratio)
  
}



# 3.
streamflow_model_drought_separate_shifted_CO2_seasonal_ratio_auto <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_drought_separate_shifted_CO2_seasonal_ratio_auto",
        "parameters" = c("a0_d", "a0_n", "a1", "a2", "a3", "a4", "a5")
      )
    )
  }
  
  # Check if parameter is a matrix - if not coerce
  if(!is.matrix(parameter_set)) {
    parameter_set <- as.matrix(parameter_set, ncol = 1)
  }
  
  # Parameters
  a0_d <- parameter_set[1, ]
  a0_n <- parameter_set[2, ]
  a1 <- parameter_set[3, ]
  a2 <- parameter_set[4, ]
  a3 <- parameter_set[5, ]
  a4 <- parameter_set[6, ]
  a5 <- parameter_set[7, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  seasonal_ratio <- catchment_data$seasonal_ratio
  CO2 <- catchment_data$CO2
  drought <- catchment_data$is_drought_year
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) #matrix(rep(precipitation, times = ncol(parameter_set)), ncol = ncol(parameter_set), byrow = FALSE)
  repeat_seasonal_ratio <- matrix(seasonal_ratio, ncol = ncol(parameter_set), nrow = length(seasonal_ratio), byrow = FALSE)
  repeat_CO2 <- matrix(CO2, ncol = ncol(parameter_set), nrow = length(CO2), byrow = FALSE) 
  repeat_drought <- matrix(drought, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE)
  repeat_a0_d <- matrix(a0_d, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a0_n <- matrix(a0_n, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a2 <- matrix(a2, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a3 <- matrix(a3, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a4 <- matrix(a4, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a5 <- matrix(a5, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Autocorrelation requirements
  boxcox_streamflow_average <- colMeans((repeat_a0_d * repeat_drought) + (repeat_a0_n * !repeat_drought) + (repeat_a1 * repeat_precipitation))
  previous_boxcox_streamflow <- boxcox_streamflow_average # assume the previous flows are the average. this means the first a2 component = 0.
  
  # CO2 requirements
  repeat_shifted_CO2 <- repeat_CO2 - repeat_a5
  repeat_shifted_CO2[repeat_shifted_CO2 < 0] <- 0 # max(0, CO2 - a5)
  
  # pre-allocate
  modelled_boxcox_streamflow <- matrix(numeric(length(repeat_precipitation)), ncol = ncol(parameter_set), nrow = length(precipitation))
  
  # Model: analysis is done row by row (timestep)
  for (time_step in seq_along(precipitation)) {
    current_boxcox_streamflow <- (repeat_a0_d[time_step, ] * repeat_drought[time_step, ]) + 
      (repeat_a0_n[time_step, ] * !repeat_drought[time_step, ]) + + 
      (repeat_a1[time_step, ] * repeat_precipitation[time_step, ]) + 
      (repeat_a2[time_step, ] * (previous_boxcox_streamflow - boxcox_streamflow_average)) +
      (repeat_a3[time_step, ] * repeat_shifted_CO2[time_step, ]) + 
      (repeat_a4[time_step, ] * repeat_seasonal_ratio[time_step, ])
    
    previous_boxcox_streamflow <- current_boxcox_streamflow
    modelled_boxcox_streamflow[time_step, ] <- current_boxcox_streamflow
  }
  
  return(modelled_boxcox_streamflow)
  
}



# 4. 
streamflow_model_separate_shifted_CO2_seasonal_ratio_auto <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_separate_shifted_CO2_seasonal_ratio_auto",
        "parameters" = c("a0", "a1", "a2", "a3", "a4", "a5")
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
  a2 <- parameter_set[3, ]
  a3 <- parameter_set[4, ]
  a4 <- parameter_set[5, ]
  a5 <- parameter_set[6, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  seasonal_ratio <- catchment_data$seasonal_ratio
  CO2 <- catchment_data$CO2
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) #matrix(rep(precipitation, times = ncol(parameter_set)), ncol = ncol(parameter_set), byrow = FALSE)
  repeat_seasonal_ratio <- matrix(seasonal_ratio, ncol = ncol(parameter_set), nrow = length(seasonal_ratio), byrow = FALSE)
  repeat_CO2 <- matrix(CO2, ncol = ncol(parameter_set), nrow = length(CO2), byrow = FALSE) 
  repeat_a0 <- matrix(a0, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a2 <- matrix(a2, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a3 <- matrix(a3, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a4 <- matrix(a4, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a5 <- matrix(a5, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Autocorrelation requirements
  boxcox_streamflow_average <- colMeans(repeat_a0 + (repeat_a1 * repeat_precipitation))
  previous_boxcox_streamflow <- boxcox_streamflow_average # assume the previous flows are the average. this means the first a2 component = 0.
  
  # pre-allocate
  modelled_boxcox_streamflow <- matrix(numeric(length(repeat_precipitation)), ncol = ncol(parameter_set), nrow = length(precipitation))
  
  
  # CO2 requirements
  repeat_shifted_CO2 <- repeat_CO2 - repeat_a5
  repeat_shifted_CO2[repeat_shifted_CO2 < 0] <- 0 # max(0, CO2 - a5)
  
  # Model: analysis is done row by row (timestep)
  for (time_step in seq_along(precipitation)) {
    current_boxcox_streamflow <- repeat_a0[time_step, ] + 
      (repeat_a1[time_step, ] * repeat_precipitation[time_step, ]) + 
      (repeat_a2[time_step, ] * (previous_boxcox_streamflow - boxcox_streamflow_average)) +
      (repeat_a3[time_step, ] * repeat_shifted_CO2[time_step, ]) + 
      (repeat_a4[time_step, ] * repeat_seasonal_ratio[time_step, ])
    
    previous_boxcox_streamflow <- current_boxcox_streamflow
    modelled_boxcox_streamflow[time_step, ] <- current_boxcox_streamflow
  }
  
  return(modelled_boxcox_streamflow)
  
}





# Run CMAES --------------------------------------------------------------------
tic()

replicate_cmaes <- function(gauge, streamflow_model, objective_function) {
  cmaes_example <- gauge |>
    catchment_data_blueprint(
      observed_data = data,
      start_stop_indexes = start_stop_indexes
    ) |>
    numerical_optimiser_setup_vary_inputs(
      streamflow_model = {{ streamflow_model }},
      objective_function = {{ objective_function }},
      bounds_and_transform_method = make_default_bounds_and_transform_methods(),
      minimise_likelihood = TRUE
    ) |>
    my_cmaes(print_monitor = TRUE) |>
    result_set() |>
    parameters_summary() 
}


REPLICATES <- 10


## Site 1: =====================================================================
site_1_gauge <- "113004A"

site_1_replicates <- replicate(
  n = REPLICATES, 
  expr = replicate_cmaes(
    site_1_gauge, 
    streamflow_model_separate_shifted_CO2_seasonal_ratio_auto, 
    constant_sd_objective_function
    ),
  simplify = FALSE
) |> 
  list_rbind() 

best_site_1_replicates <- site_1_replicates |> 
  dplyr::slice_min(
    loglikelihood
  )




## Site 2: =====================================================================
site_2_gauge <- "230210"


site_2_replicates <- replicate(
  n = REPLICATES, 
  expr = replicate_cmaes(
    site_2_gauge, 
    streamflow_model_drought_separate_shifted_CO2_seasonal_ratio_auto, 
    constant_sd_objective_function
  ),
  simplify = FALSE
) |> 
  list_rbind() 

best_site_2_replicates <- site_2_replicates |> 
  dplyr::slice_min(
    loglikelihood
  )



## Site 3: =====================================================================
site_3_gauge <- "302214"

site_3_replicates <- replicate(
  n = REPLICATES, 
  expr = replicate_cmaes(
    site_3_gauge, 
    streamflow_model_separate_shifted_CO2, 
    constant_sd_objective_function
  ),
  simplify = FALSE
) |> 
  list_rbind() 

best_site_3_replicates <- site_3_replicates |> 
  dplyr::slice_min(
    loglikelihood
  )



## Site 4: =====================================================================
site_4_gauge <- "407220"

site_4_replicates <- replicate(
  n = REPLICATES, 
  expr = replicate_cmaes(
    site_1_gauge, 
    streamflow_model_separate_shifted_CO2_seasonal_ratio_auto, 
    constant_sd_objective_function
  ),
  simplify = FALSE
) |> 
  list_rbind() 

best_site_4_replicates <- site_4_replicates |> 
  dplyr::slice_min(
    loglikelihood
  )



## Site 5: =====================================================================
site_5_gauge <- "614044"

site_5_replicates <- replicate(
  n = REPLICATES, 
  expr = replicate_cmaes(
    site_1_gauge, 
    streamflow_model_separate_shifted_CO2_seasonal_ratio_auto, 
    constant_sd_objective_function
  ),
  simplify = FALSE
) |> 
  list_rbind() 

best_site_5_replicates <- site_5_replicates |> 
  dplyr::slice_min(
    loglikelihood
  )


toc()
