# Does time beat CO2 + a5 parameters? ------------------------------------------


# Import libraries required ----------------------------------------------------
pacman::p_load(tidyverse, cmaesr, smoof, truncnorm, sloop, ggExtra, tictoc, dream, furrr) # Once this becomes a package this step is not required


# There are multiple packages called dream.
# Install the following:
# - install.packages("dream", repos="http://R-Forge.R-project.org") 
# OR
# - https://rdrr.io/rforge/dream/





# Catchments for comparison ----------------------------------------------------
best_CO2_non_CO2_per_catchment <- read_csv(
  "./Modelling/Results/CMAES/best_CO2_non_CO2_per_catchment_CMAES.csv",
  show_col_types = FALSE
) 

best_model_per_gauge <- best_CO2_non_CO2_per_catchment |> 
  select(gauge, streamflow_model, AIC) |> 
  distinct() |> 
  slice_min(order_by = AIC, by = gauge) |> 
  mutate(
    contains_CO2 = str_detect(streamflow_model, "CO2")
  ) |> 
  filter(contains_CO2)


models_needed <- best_model_per_gauge |> 
  pull(streamflow_model) |> 
  unique()



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
  mutate(q_mm = q_mm + .Machine$double.eps^0.5) |> 
  group_by(gauge) |> 
  # Make CO2 = count
  mutate(
    CO2 = as.numeric(row_number())
  ) |> 
  ungroup()

gauge_information <- readr::read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE
)




# TEST STREAMFLOW MODEL --------------------------------------------------------
# SOMETHING TO NOTE:
# - changing parameter names caused unexpected behaviour
# - I am not sure where the issue is.

# Remove a5 parameter - just have time
streamflow_model_intercept_shifted_time <- function(catchment_data, parameter_set) { 
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_intercept_shifted_time",
        "parameters" = c("a0", "a1", "a3_intercept")
      )
    )
  } 
  
  # Check if parameter is a matrix - if not coerce
  if(!is.matrix(parameter_set)) {
    parameter_set <- as.matrix(parameter_set, ncol = 1)
  }
  
  # If directly inputting catchment data into model then 
  # the following function will recursively call the streamflow model
  # for each list in catchment_data (start_stop)
  if(sloop::s3_class(catchment_data)[1] == "catchment_data"){ # Check type of catchment data - make this a function call for all streamflow models
    streamflow_results <- catchment_data_directly_to_streamflow_model(
      catchment_data = catchment_data,
      parameter_set = parameter_set, 
      streamflow_model = streamflow_model_intercept_shifted_time
    )
    return(streamflow_results)
  }
  
  
  
  # Parameters
  a0 <- parameter_set[1, ]
  a1 <- parameter_set[2, ]
  a3_intercept <- parameter_set[3, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  CO2 <- catchment_data$CO2
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) 
  repeat_CO2 <- matrix(CO2, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) 
  repeat_a0 <- matrix(a0, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a3_intercept <- matrix(a3_intercept, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)

  
  # Model
  repeat_a0 + (repeat_a1 * repeat_precipitation) + (repeat_a3_intercept * repeat_CO2)
  
} 
streamflow_model_intercept_shifted_time_auto <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_intercept_shifted_time_auto",
        "parameters" = c("a0", "a1", "a2", "a3_intercept")
      )
    )
  }
  
  # Check if parameter is a matrix - if not coerce
  if(!is.matrix(parameter_set)) {
    parameter_set <- as.matrix(parameter_set, ncol = 1)
  }
  
  # If directly inputting catchment data into model then 
  # the following function will recursively call the streamflow model
  # for each list in catchment_data (start_stop)
  if(sloop::s3_class(catchment_data)[1] == "catchment_data"){ # Check type of catchment data - make this a function call for all streamflow models
    streamflow_results <- catchment_data_directly_to_streamflow_model(
      catchment_data = catchment_data,
      parameter_set = parameter_set, 
      streamflow_model = streamflow_model_intercept_shifted_time_auto
    )
    return(streamflow_results)
  }
  
  # Parameters
  a0 <- parameter_set[1, ]
  a1 <- parameter_set[2, ]
  a2 <- parameter_set[3, ]
  a3_intercept <- parameter_set[4, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  CO2 <- catchment_data$CO2
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) #matrix(rep(precipitation, times = ncol(parameter_set)), ncol = ncol(parameter_set), byrow = FALSE)
  repeat_CO2 <- matrix(CO2, ncol = ncol(parameter_set), nrow = length(CO2), byrow = FALSE) 
  repeat_a0 <- matrix(a0, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a2 <- matrix(a2, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a3_intercept <- matrix(a3_intercept, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  
  # Autocorrelation requirements
  boxcox_streamflow_average <- colMeans(repeat_a0 + (repeat_a1 * repeat_precipitation))
  previous_boxcox_streamflow <- boxcox_streamflow_average # assume the previous flows are the average. this means the first a2 component = 0.
  
  # pre-allocate
  modelled_boxcox_streamflow <- matrix(numeric(length(repeat_precipitation)), ncol = ncol(parameter_set), nrow = length(precipitation))
  
  
  # Model: analysis is done row by row (timestep)
  for (time_step in seq_along(precipitation)) {
    current_boxcox_streamflow <- repeat_a0[time_step, ] + 
      (repeat_a1[time_step, ] * repeat_precipitation[time_step, ]) + 
      (repeat_a2[time_step, ] * (previous_boxcox_streamflow - boxcox_streamflow_average)) +
      (repeat_a3_intercept[time_step, ] * repeat_CO2[time_step, ])
    
    previous_boxcox_streamflow <- current_boxcox_streamflow
    modelled_boxcox_streamflow[time_step, ] <- current_boxcox_streamflow
  }
  
  return(modelled_boxcox_streamflow)
} 
streamflow_model_intercept_shifted_time_seasonal_ratio <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_intercept_shifted_time_seasonal_ratio",
        "parameters" = c("a0", "a1", "a3_intercept", "a4")
      )
    )
  } 
  
  # Check if parameter is a matrix - if not coerce
  if(!is.matrix(parameter_set)) {
    parameter_set <- as.matrix(parameter_set, ncol = 1)
  }
  
  # If directly inputting catchment data into model then 
  # the following function will recursively call the streamflow model
  # for each list in catchment_data (start_stop)
  if(sloop::s3_class(catchment_data)[1] == "catchment_data"){ # Check type of catchment data - make this a function call for all streamflow models
    streamflow_results <- catchment_data_directly_to_streamflow_model(
      catchment_data = catchment_data,
      parameter_set = parameter_set, 
      streamflow_model = streamflow_model_intercept_shifted_time_seasonal_ratio
    )
    return(streamflow_results)
  }
  
  # Parameters
  a0 <- parameter_set[1, ]
  a1 <- parameter_set[2, ]
  a3_intercept <- parameter_set[3, ]
  a4 <- parameter_set[4, ]
  
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
  repeat_a3_intercept <- matrix(a3_intercept, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a4 <- matrix(a4, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Model
  repeat_a0 + (repeat_a1 * repeat_precipitation) + (repeat_a3_intercept * repeat_CO2) + (repeat_a4 * repeat_seasonal_ratio)
  
} 
streamflow_model_intercept_shifted_time_seasonal_ratio_auto <- function(catchment_data, parameter_set) {
  
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_intercept_shifted_time_seasonal_ratio_auto",
        "parameters" = c("a0", "a1", "a2", "a3_intercept", "a4", "a5")
      )
    )
  }
  
  
  # Check if parameter is a matrix - if not coerce
  if(!is.matrix(parameter_set)) {
    parameter_set <- as.matrix(parameter_set, ncol = 1)
  }
  
  # If directly inputting catchment data into model then 
  # the following function will recursively call the streamflow model
  # for each list in catchment_data (start_stop)
  if(sloop::s3_class(catchment_data)[1] == "catchment_data"){ # Check type of catchment data - make this a function call for all streamflow models
    streamflow_results <- catchment_data_directly_to_streamflow_model(
      catchment_data = catchment_data,
      parameter_set = parameter_set, 
      streamflow_model = streamflow_model_intercept_shifted_time_seasonal_ratio_auto
    )
    return(streamflow_results)
  }
  
  # Parameters
  a0 <- parameter_set[1, ]
  a1 <- parameter_set[2, ]
  a2 <- parameter_set[3, ]
  a3_intercept <- parameter_set[4, ]
  a4 <- parameter_set[5, ]
  
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
  repeat_a3_intercept <- matrix(a3_intercept, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a4 <- matrix(a4, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)

  
  # Autocorrelation requirements
  boxcox_streamflow_average <- colMeans(repeat_a0 + (repeat_a1 * repeat_precipitation))
  previous_boxcox_streamflow <- boxcox_streamflow_average # assume the previous flows are the average. this means the first a2 component = 0.
  
  # pre-allocate
  modelled_boxcox_streamflow <- matrix(numeric(length(repeat_precipitation)), ncol = ncol(parameter_set), nrow = length(precipitation))
  

  
  # Model: analysis is done row by row (timestep)
  for (time_step in seq_along(precipitation)) {
    current_boxcox_streamflow <- repeat_a0[time_step, ] + 
      (repeat_a1[time_step, ] * repeat_precipitation[time_step, ]) + 
      (repeat_a2[time_step, ] * (previous_boxcox_streamflow - boxcox_streamflow_average)) +
      (repeat_a3_intercept[time_step, ] * repeat_CO2[time_step, ]) + 
      (repeat_a4[time_step, ] * repeat_seasonal_ratio[time_step, ])
    
    previous_boxcox_streamflow <- current_boxcox_streamflow
    modelled_boxcox_streamflow[time_step, ] <- current_boxcox_streamflow
  }
  
  return(modelled_boxcox_streamflow)
  
} 
streamflow_model_slope_shifted_time <- function(catchment_data, parameter_set) { 
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_slope_shifted_time",
        "parameters" = c("a0", "a1", "a3_slope")
      )
    )
  } 
  
  # Check if parameter is a matrix - if not coerce
  if(!is.matrix(parameter_set)) {
    parameter_set <- as.matrix(parameter_set, ncol = 1)
  }
  
  # If directly inputting catchment data into model then 
  # the following function will recursively call the streamflow model
  # for each list in catchment_data (start_stop)
  if(sloop::s3_class(catchment_data)[1] == "catchment_data"){ # Check type of catchment data - make this a function call for all streamflow models
    streamflow_results <- catchment_data_directly_to_streamflow_model(
      catchment_data = catchment_data,
      parameter_set = parameter_set, 
      streamflow_model = streamflow_model_slope_shifted_time
    )
    return(streamflow_results)
  }
  
  # Parameters
  a0 <- parameter_set[1, ]
  a1 <- parameter_set[2, ]
  a3_slope <- parameter_set[3, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  CO2 <- catchment_data$CO2
  
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) #matrix(rep(precipitation, times = ncol(parameter_set)), ncol = ncol(parameter_set), byrow = FALSE)
  repeat_CO2 <- matrix(CO2, ncol = ncol(parameter_set), nrow = length(CO2), byrow = FALSE) 
  repeat_a0 <- matrix(a0, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a3_slope <- matrix(a3_slope, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  
  # Model
  repeat_a0 + (repeat_a1 * repeat_precipitation) + (repeat_a3_slope * repeat_precipitation * repeat_CO2)
  
} 
streamflow_model_slope_shifted_time_auto <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_slope_shifted_time_auto",
        "parameters" = c("a0", "a1", "a2", "a3_slope")
      )
    )
  }
  
  # Check if parameter is a matrix - if not coerce
  if(!is.matrix(parameter_set)) {
    parameter_set <- as.matrix(parameter_set, ncol = 1)
  }
  
  # If directly inputting catchment data into model then 
  # the following function will recursively call the streamflow model
  # for each list in catchment_data (start_stop)
  if(sloop::s3_class(catchment_data)[1] == "catchment_data"){ # Check type of catchment data - make this a function call for all streamflow models
    streamflow_results <- catchment_data_directly_to_streamflow_model(
      catchment_data = catchment_data,
      parameter_set = parameter_set, 
      streamflow_model = streamflow_model_slope_shifted_time_auto
    )
    return(streamflow_results)
  }
  
  # Parameters
  a0 <- parameter_set[1, ]
  a1 <- parameter_set[2, ]
  a2 <- parameter_set[3, ]
  a3_slope <- parameter_set[4, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  CO2 <- catchment_data$CO2
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) #matrix(rep(precipitation, times = ncol(parameter_set)), ncol = ncol(parameter_set), byrow = FALSE)
  repeat_CO2 <- matrix(CO2, ncol = ncol(parameter_set), nrow = length(CO2), byrow = FALSE) 
  repeat_a0 <- matrix(a0, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a2 <- matrix(a2, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a3_slope <- matrix(a3_slope, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Autocorrelation requirements
  boxcox_streamflow_average <- colMeans(repeat_a0 + (repeat_a1 * repeat_precipitation))
  previous_boxcox_streamflow <- boxcox_streamflow_average # assume the previous flows are the average. this means the first a2 component = 0.
  
  # pre-allocate
  modelled_boxcox_streamflow <- matrix(numeric(length(repeat_precipitation)), ncol = ncol(parameter_set), nrow = length(precipitation))
  
  
  # Model: analysis is done row by row (timestep)
  for (time_step in seq_along(precipitation)) {
    current_boxcox_streamflow <- repeat_a0[time_step, ] + 
      (repeat_a1[time_step, ] * repeat_precipitation[time_step, ]) + 
      (repeat_a2[time_step, ] * (previous_boxcox_streamflow - boxcox_streamflow_average)) +
      (repeat_a3_slope[time_step, ] * repeat_precipitation[time_step, ] * repeat_CO2[time_step, ])
    
    previous_boxcox_streamflow <- current_boxcox_streamflow
    modelled_boxcox_streamflow[time_step, ] <- current_boxcox_streamflow
  }
  
  return(modelled_boxcox_streamflow)
} 
streamflow_model_slope_shifted_time_seasonal_ratio <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_slope_shifted_time_seasonal_ratio",
        "parameters" = c("a0", "a1", "a3_slope", "a4")
      )
    )
  } 
  
  # Check if parameter is a matrix - if not coerce
  if(!is.matrix(parameter_set)) {
    parameter_set <- as.matrix(parameter_set, ncol = 1)
  }
  
  # If directly inputting catchment data into model then 
  # the following function will recursively call the streamflow model
  # for each list in catchment_data (start_stop)
  if(sloop::s3_class(catchment_data)[1] == "catchment_data"){ # Check type of catchment data - make this a function call for all streamflow models
    streamflow_results <- catchment_data_directly_to_streamflow_model(
      catchment_data = catchment_data,
      parameter_set = parameter_set, 
      streamflow_model = streamflow_model_slope_shifted_time_seasonal_ratio
    )
    return(streamflow_results)
  }
  
  # Parameters
  a0 <- parameter_set[1, ]
  a1 <- parameter_set[2, ]
  a3_slope <- parameter_set[3, ]
  a4 <- parameter_set[4, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  seasonal_ratio <- catchment_data$seasonal_ratio
  CO2 <- catchment_data$CO2
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) 
  repeat_seasonal_ratio <- matrix(seasonal_ratio, ncol = ncol(parameter_set), nrow = length(seasonal_ratio), byrow = FALSE)
  repeat_CO2 <- matrix(CO2, ncol = ncol(parameter_set), nrow = length(CO2), byrow = FALSE) 
  
  repeat_a0 <- matrix(a0, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a3_slope <- matrix(a3_slope, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a4 <- matrix(a4, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  repeat_a0 + (repeat_a1 * repeat_precipitation) + (repeat_a3_slope * repeat_precipitation * repeat_CO2) + (repeat_a4 * repeat_seasonal_ratio)
  
} 
streamflow_model_slope_shifted_time_seasonal_ratio_auto <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_slope_shifted_time_seasonal_ratio_auto",
        "parameters" = c("a0", "a1", "a2", "a3_slope", "a4")
      )
    )
  }
  
  # Check if parameter is a matrix - if not coerce
  if(!is.matrix(parameter_set)) {
    parameter_set <- as.matrix(parameter_set, ncol = 1)
  }
  
  # If directly inputting catchment data into model then 
  # the following function will recursively call the streamflow model
  # for each list in catchment_data (start_stop)
  if(sloop::s3_class(catchment_data)[1] == "catchment_data"){ # Check type of catchment data - make this a function call for all streamflow models
    streamflow_results <- catchment_data_directly_to_streamflow_model(
      catchment_data = catchment_data,
      parameter_set = parameter_set, 
      streamflow_model = streamflow_model_slope_shifted_time_seasonal_ratio_auto
    )
    return(streamflow_results)
  }
  
  # Parameters
  a0 <- parameter_set[1, ]
  a1 <- parameter_set[2, ]
  a2 <- parameter_set[3, ]
  a3_slope <- parameter_set[4, ]
  a4 <- parameter_set[5, ]
  
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
  repeat_a3_slope <- matrix(a3_slope, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a4 <- matrix(a4, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Autocorrelation requirements
  boxcox_streamflow_average <- colMeans(repeat_a0 + (repeat_a1 * repeat_precipitation))
  previous_boxcox_streamflow <- boxcox_streamflow_average # assume the previous flows are the average. this means the first a2 component = 0.
  
  # pre-allocate
  modelled_boxcox_streamflow <- matrix(numeric(length(repeat_precipitation)), ncol = ncol(parameter_set), nrow = length(precipitation))
  

  
  # Model: analysis is done row by row (timestep)
  for (time_step in seq_along(precipitation)) {
    current_boxcox_streamflow <- repeat_a0[time_step, ] + 
      (repeat_a1[time_step, ] * repeat_precipitation[time_step, ]) + 
      (repeat_a2[time_step, ] * (previous_boxcox_streamflow - boxcox_streamflow_average)) +
      (repeat_a3_slope[time_step, ] * repeat_precipitation[time_step, ] * repeat_CO2[time_step, ]) + 
      (repeat_a4[time_step, ] * repeat_seasonal_ratio[time_step, ])
    
    previous_boxcox_streamflow <- current_boxcox_streamflow
    modelled_boxcox_streamflow[time_step, ] <- current_boxcox_streamflow
  }
  
  return(modelled_boxcox_streamflow)
  
}

streamflow_model_drought_intercept_shifted_time <- function(catchment_data, parameter_set) { 
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_drought_intercept_shifted_time",
        "parameters" = c("a0_d", "a0_n", "a1", "a3_intercept")
      )
    )
  } 
  
  # Check if parameter is a matrix - if not coerce
  if(!is.matrix(parameter_set)) {
    parameter_set <- as.matrix(parameter_set, ncol = 1)
  }
  
  # If directly inputting catchment data into model then 
  # the following function will recursively call the streamflow model
  # for each list in catchment_data (start_stop)
  if(sloop::s3_class(catchment_data)[1] == "catchment_data"){ # Check type of catchment data - make this a function call for all streamflow models
    streamflow_results <- catchment_data_directly_to_streamflow_model(
      catchment_data = catchment_data,
      parameter_set = parameter_set, 
      streamflow_model = streamflow_model_drought_intercept_shifted_time
    )
    return(streamflow_results)
  }
  
  # Parameters
  a0_d <- parameter_set[1, ]
  a0_n <- parameter_set[2, ]
  a1 <- parameter_set[3, ]
  a3_intercept <- parameter_set[4, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  drought <- catchment_data$is_drought_year
  CO2 <- catchment_data$CO2
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) 
  repeat_CO2 <- matrix(CO2, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) 
  repeat_drought <- matrix(drought, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE)
  repeat_a0_d <- matrix(a0_d, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a0_n <- matrix(a0_n, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a3_intercept <- matrix(a3_intercept, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)

  
  # Model
  (repeat_a0_d * repeat_drought) + (repeat_a0_n * !repeat_drought) + (repeat_a1 * repeat_precipitation) + (repeat_a3_intercept * repeat_CO2)
  
}
streamflow_model_drought_intercept_shifted_time_auto <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_drought_intercept_shifted_time_auto",
        "parameters" = c("a0_d", "a0_n", "a1", "a2", "a3_intercept")
      )
    )
  }
  
  # Check if parameter is a matrix - if not coerce
  if(!is.matrix(parameter_set)) {
    parameter_set <- as.matrix(parameter_set, ncol = 1)
  }
  
  # If directly inputting catchment data into model then 
  # the following function will recursively call the streamflow model
  # for each list in catchment_data (start_stop)
  if(sloop::s3_class(catchment_data)[1] == "catchment_data"){ # Check type of catchment data - make this a function call for all streamflow models
    streamflow_results <- catchment_data_directly_to_streamflow_model(
      catchment_data = catchment_data,
      parameter_set = parameter_set, 
      streamflow_model = streamflow_model_drought_intercept_shifted_time_auto
    )
    return(streamflow_results)
  }
  
  # Parameters
  a0_d <- parameter_set[1, ]
  a0_n <- parameter_set[2, ]
  a1 <- parameter_set[3, ]
  a2 <- parameter_set[4, ]
  a3_intercept <- parameter_set[5, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  CO2 <- catchment_data$CO2
  drought <- catchment_data$is_drought_year
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) #matrix(rep(precipitation, times = ncol(parameter_set)), ncol = ncol(parameter_set), byrow = FALSE)
  repeat_CO2 <- matrix(CO2, ncol = ncol(parameter_set), nrow = length(CO2), byrow = FALSE) 
  repeat_drought <- matrix(drought, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE)
  repeat_a0_d <- matrix(a0_d, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a0_n <- matrix(a0_n, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a2 <- matrix(a2, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a3_intercept <- matrix(a3_intercept, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Autocorrelation requirements
  boxcox_streamflow_average <-colMeans((repeat_a0_d * repeat_drought) + (repeat_a0_n * !repeat_drought) + (repeat_a1 * repeat_precipitation))
  previous_boxcox_streamflow <- boxcox_streamflow_average # assume the previous flows are the average. this means the first a2 component = 0.
  
  # pre-allocate
  modelled_boxcox_streamflow <- matrix(numeric(length(repeat_precipitation)), ncol = ncol(parameter_set), nrow = length(precipitation))
  
  
  # Model: analysis is done row by row (timestep)
  for (time_step in seq_along(precipitation)) {
    current_boxcox_streamflow <- (repeat_a0_d[time_step, ] * repeat_drought[time_step, ]) + 
      (repeat_a0_n[time_step, ] * !repeat_drought[time_step, ]) + 
      (repeat_a1[time_step, ] * repeat_precipitation[time_step, ]) + 
      (repeat_a2[time_step, ] * (previous_boxcox_streamflow - boxcox_streamflow_average)) +
      (repeat_a3_intercept[time_step, ] * repeat_CO2[time_step, ])
    
    previous_boxcox_streamflow <- current_boxcox_streamflow
    modelled_boxcox_streamflow[time_step, ] <- current_boxcox_streamflow
  }
  
  return(modelled_boxcox_streamflow)
}
streamflow_model_drought_intercept_shifted_time_seasonal_ratio <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_drought_intercept_shifted_time_seasonal_ratio",
        "parameters" = c("a0_d", "a0_n", "a1", "a3_intercept", "a4")
      )
    )
  } 
  
  # Check if parameter is a matrix - if not coerce
  if(!is.matrix(parameter_set)) {
    parameter_set <- as.matrix(parameter_set, ncol = 1)
  }
  
  # If directly inputting catchment data into model then 
  # the following function will recursively call the streamflow model
  # for each list in catchment_data (start_stop)
  if(sloop::s3_class(catchment_data)[1] == "catchment_data"){ # Check type of catchment data - make this a function call for all streamflow models
    streamflow_results <- catchment_data_directly_to_streamflow_model(
      catchment_data = catchment_data,
      parameter_set = parameter_set, 
      streamflow_model = streamflow_model_drought_intercept_shifted_time_seasonal_ratio
    )
    return(streamflow_results)
  }
  
  # Parameters
  a0_d <- parameter_set[1, ]
  a0_n <- parameter_set[2, ]
  a1 <- parameter_set[3, ]
  a3_intercept <- parameter_set[4, ]
  a4 <- parameter_set[5, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  seasonal_ratio <- catchment_data$seasonal_ratio
  drought <- catchment_data$is_drought_year
  CO2 <- catchment_data$CO2
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) #matrix(rep(precipitation, times = ncol(parameter_set)), ncol = ncol(parameter_set), byrow = FALSE)
  repeat_seasonal_ratio <- matrix(seasonal_ratio, ncol = ncol(parameter_set), nrow = length(seasonal_ratio), byrow = FALSE)
  repeat_CO2 <- matrix(CO2, ncol = ncol(parameter_set), nrow = length(CO2), byrow = FALSE) 
  repeat_drought <- matrix(drought, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE)
  repeat_a0_d <- matrix(a0_d, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a0_n <- matrix(a0_n, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a3_intercept <- matrix(a3_intercept, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a4 <- matrix(a4, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Model
  (repeat_a0_d * repeat_drought) + (repeat_a0_n * !repeat_drought) + (repeat_a1 * repeat_precipitation) + (repeat_a3_intercept * repeat_CO2) + (repeat_a4 * repeat_seasonal_ratio)
  
}
streamflow_model_drought_intercept_shifted_time_seasonal_ratio_auto <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_drought_intercept_shifted_time_seasonal_ratio_auto",
        "parameters" = c("a0_d", "a0_n", "a1", "a2", "a3_intercept", "a4")
      )
    )
  }
  
  # Check if parameter is a matrix - if not coerce
  if(!is.matrix(parameter_set)) {
    parameter_set <- as.matrix(parameter_set, ncol = 1)
  }
  
  # If directly inputting catchment data into model then 
  # the following function will recursively call the streamflow model
  # for each list in catchment_data (start_stop)
  if(sloop::s3_class(catchment_data)[1] == "catchment_data"){ # Check type of catchment data - make this a function call for all streamflow models
    streamflow_results <- catchment_data_directly_to_streamflow_model(
      catchment_data = catchment_data,
      parameter_set = parameter_set, 
      streamflow_model = streamflow_model_drought_intercept_shifted_time_seasonal_ratio_auto
    )
    return(streamflow_results)
  }
  
  # Parameters
  a0_d <- parameter_set[1, ]
  a0_n <- parameter_set[2, ]
  a1 <- parameter_set[3, ]
  a2 <- parameter_set[4, ]
  a3_intercept <- parameter_set[5, ]
  a4 <- parameter_set[6, ]
  
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
  repeat_a3_intercept <- matrix(a3_intercept, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a4 <- matrix(a4, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Autocorrelation requirements
  boxcox_streamflow_average <- colMeans((repeat_a0_d * repeat_drought) + (repeat_a0_n * !repeat_drought) + (repeat_a1 * repeat_precipitation))
  previous_boxcox_streamflow <- boxcox_streamflow_average # assume the previous flows are the average. this means the first a2 component = 0.
  
  
  # pre-allocate
  modelled_boxcox_streamflow <- matrix(numeric(length(repeat_precipitation)), ncol = ncol(parameter_set), nrow = length(precipitation))
  
  # Model: analysis is done row by row (timestep)
  for (time_step in seq_along(precipitation)) {
    current_boxcox_streamflow <- (repeat_a0_d[time_step, ] * repeat_drought[time_step, ]) + 
      (repeat_a0_n[time_step, ] * !repeat_drought[time_step, ]) + + 
      (repeat_a1[time_step, ] * repeat_precipitation[time_step, ]) + 
      (repeat_a2[time_step, ] * (previous_boxcox_streamflow - boxcox_streamflow_average)) +
      (repeat_a3_intercept[time_step, ] * repeat_CO2[time_step, ]) + 
      (repeat_a4[time_step, ] * repeat_seasonal_ratio[time_step, ])
    
    previous_boxcox_streamflow <- current_boxcox_streamflow
    modelled_boxcox_streamflow[time_step, ] <- current_boxcox_streamflow
  }
  
  return(modelled_boxcox_streamflow)
  
}
streamflow_model_drought_slope_shifted_time <- function(catchment_data, parameter_set) { 
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_drought_slope_shifted_time",
        "parameters" = c("a0_d", "a0_n", "a1", "a3_slope")
      )
    )
  } 
  
  # Check if parameter is a matrix - if not coerce
  if(!is.matrix(parameter_set)) {
    parameter_set <- as.matrix(parameter_set, ncol = 1)
  }
  
  # If directly inputting catchment data into model then 
  # the following function will recursively call the streamflow model
  # for each list in catchment_data (start_stop)
  if(sloop::s3_class(catchment_data)[1] == "catchment_data"){ # Check type of catchment data - make this a function call for all streamflow models
    streamflow_results <- catchment_data_directly_to_streamflow_model(
      catchment_data = catchment_data,
      parameter_set = parameter_set, 
      streamflow_model = streamflow_model_drought_slope_shifted_time
    )
    return(streamflow_results)
  }
  
  # Parameters
  a0_d <- parameter_set[1, ]
  a0_n <- parameter_set[2, ]
  a1 <- parameter_set[3, ]
  a3_slope <- parameter_set[4, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  drought <- catchment_data$is_drought_year
  CO2 <- catchment_data$CO2
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) 
  repeat_CO2 <- matrix(CO2, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) 
  repeat_drought <- matrix(drought, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE)
  repeat_a0_d <- matrix(a0_d, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a0_n <- matrix(a0_n, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a3_slope <- matrix(a3_slope, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  

  
  # Model
  (repeat_a0_d * repeat_drought) + (repeat_a0_n * !repeat_drought) + (repeat_a1 * repeat_precipitation) + (repeat_a3_slope * repeat_precipitation * repeat_CO2)
  
}
streamflow_model_drought_slope_shifted_time_auto <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_drought_slope_shifted_time_auto",
        "parameters" = c("a0_d", "a0_n", "a1", "a2", "a3_slope")
      )
    )
  }
  
  # Check if parameter is a matrix - if not coerce
  if(!is.matrix(parameter_set)) {
    parameter_set <- as.matrix(parameter_set, ncol = 1)
  }
  
  # If directly inputting catchment data into model then 
  # the following function will recursively call the streamflow model
  # for each list in catchment_data (start_stop)
  if(sloop::s3_class(catchment_data)[1] == "catchment_data"){ # Check type of catchment data - make this a function call for all streamflow models
    streamflow_results <- catchment_data_directly_to_streamflow_model(
      catchment_data = catchment_data,
      parameter_set = parameter_set, 
      streamflow_model = streamflow_model_drought_slope_shifted_time_auto
    )
    return(streamflow_results)
  }
  
  # Parameters
  a0_d <- parameter_set[1, ]
  a0_n <- parameter_set[2, ]
  a1 <- parameter_set[3, ]
  a2 <- parameter_set[4, ]
  a3_slope <- parameter_set[5, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  CO2 <- catchment_data$CO2
  drought <- catchment_data$is_drought_year
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) #matrix(rep(precipitation, times = ncol(parameter_set)), ncol = ncol(parameter_set), byrow = FALSE)
  repeat_CO2 <- matrix(CO2, ncol = ncol(parameter_set), nrow = length(CO2), byrow = FALSE) 
  repeat_drought <- matrix(drought, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE)
  repeat_a0_d <- matrix(a0_d, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a0_n <- matrix(a0_n, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a2 <- matrix(a2, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a3_slope <- matrix(a3_slope, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Autocorrelation requirements
  boxcox_streamflow_average <-colMeans((repeat_a0_d * repeat_drought) + (repeat_a0_n * !repeat_drought) + (repeat_a1 * repeat_precipitation))
  previous_boxcox_streamflow <- boxcox_streamflow_average # assume the previous flows are the average. this means the first a2 component = 0.
  
  # pre-allocate
  modelled_boxcox_streamflow <- matrix(numeric(length(repeat_precipitation)), ncol = ncol(parameter_set), nrow = length(precipitation))
  
  
  # Model: analysis is done row by row (timestep)
  for (time_step in seq_along(precipitation)) {
    current_boxcox_streamflow <- (repeat_a0_d[time_step, ] * repeat_drought[time_step, ]) + 
      (repeat_a0_n[time_step, ] * !repeat_drought[time_step, ]) + 
      (repeat_a1[time_step, ] * repeat_precipitation[time_step, ]) + 
      (repeat_a2[time_step, ] * (previous_boxcox_streamflow - boxcox_streamflow_average)) +
      (repeat_a3_slope[time_step, ] * repeat_precipitation[time_step, ] * repeat_CO2[time_step, ])
    
    previous_boxcox_streamflow <- current_boxcox_streamflow
    modelled_boxcox_streamflow[time_step, ] <- current_boxcox_streamflow
  }
  
  return(modelled_boxcox_streamflow)
}
streamflow_model_drought_slope_shifted_time_seasonal_ratio <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_drought_slope_shifted_time_seasonal_ratio",
        "parameters" = c("a0_d", "a0_n", "a1", "a3_slope", "a4")
      )
    )
  } 
  
  # Check if parameter is a matrix - if not coerce
  if(!is.matrix(parameter_set)) {
    parameter_set <- as.matrix(parameter_set, ncol = 1)
  }
  
  # If directly inputting catchment data into model then 
  # the following function will recursively call the streamflow model
  # for each list in catchment_data (start_stop)
  if(sloop::s3_class(catchment_data)[1] == "catchment_data"){ # Check type of catchment data - make this a function call for all streamflow models
    streamflow_results <- catchment_data_directly_to_streamflow_model(
      catchment_data = catchment_data,
      parameter_set = parameter_set, 
      streamflow_model = streamflow_model_drought_slope_shifted_time_seasonal_ratio
    )
    return(streamflow_results)
  }
  
  # Parameters
  a0_d <- parameter_set[1, ]
  a0_n <- parameter_set[2, ]
  a1 <- parameter_set[3, ]
  a3_slope <- parameter_set[4, ]
  a4 <- parameter_set[5, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  seasonal_ratio <- catchment_data$seasonal_ratio
  drought <- catchment_data$is_drought_year
  CO2 <- catchment_data$CO2
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) #matrix(rep(precipitation, times = ncol(parameter_set)), ncol = ncol(parameter_set), byrow = FALSE)
  repeat_seasonal_ratio <- matrix(seasonal_ratio, ncol = ncol(parameter_set), nrow = length(seasonal_ratio), byrow = FALSE)
  repeat_CO2 <- matrix(CO2, ncol = ncol(parameter_set), nrow = length(CO2), byrow = FALSE) 
  repeat_drought <- matrix(drought, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE)
  repeat_a0_d <- matrix(a0_d, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a0_n <- matrix(a0_n, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a3_slope <- matrix(a3_slope, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a4 <- matrix(a4, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)

  
  # Model
  (repeat_a0_d * repeat_drought) + (repeat_a0_n * !repeat_drought) + (repeat_a1 * repeat_precipitation) + (repeat_a3_slope * repeat_precipitation * repeat_CO2) + (repeat_a4 * repeat_seasonal_ratio)
  
}
streamflow_model_drought_slope_shifted_time_seasonal_ratio_auto <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_drought_slope_shifted_time_seasonal_ratio_auto",
        "parameters" = c("a0_d", "a0_n", "a1", "a2", "a3_slope", "a4")
      )
    )
  }
  
  # Check if parameter is a matrix - if not coerce
  if(!is.matrix(parameter_set)) {
    parameter_set <- as.matrix(parameter_set, ncol = 1)
  }
  
  # If directly inputting catchment data into model then 
  # the following function will recursively call the streamflow model
  # for each list in catchment_data (start_stop)
  if(sloop::s3_class(catchment_data)[1] == "catchment_data"){ # Check type of catchment data - make this a function call for all streamflow models
    streamflow_results <- catchment_data_directly_to_streamflow_model(
      catchment_data = catchment_data,
      parameter_set = parameter_set, 
      streamflow_model = streamflow_model_drought_slope_shifted_time_seasonal_ratio_auto
    )
    return(streamflow_results)
  }
  
  # Parameters
  a0_d <- parameter_set[1, ]
  a0_n <- parameter_set[2, ]
  a1 <- parameter_set[3, ]
  a2 <- parameter_set[4, ]
  a3_slope <- parameter_set[5, ]
  a4 <- parameter_set[6, ]
  
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
  repeat_a3_slope <- matrix(a3_slope, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a4 <- matrix(a4, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Autocorrelation requirements
  boxcox_streamflow_average <- colMeans((repeat_a0_d * repeat_drought) + (repeat_a0_n * !repeat_drought) + (repeat_a1 * repeat_precipitation))
  previous_boxcox_streamflow <- boxcox_streamflow_average # assume the previous flows are the average. this means the first a2 component = 0.
  
  # pre-allocate
  modelled_boxcox_streamflow <- matrix(numeric(length(repeat_precipitation)), ncol = ncol(parameter_set), nrow = length(precipitation))
  
  # Model: analysis is done row by row (timestep)
  for (time_step in seq_along(precipitation)) {
    current_boxcox_streamflow <- (repeat_a0_d[time_step, ] * repeat_drought[time_step, ]) + 
      (repeat_a0_n[time_step, ] * !repeat_drought[time_step, ]) + + 
      (repeat_a1[time_step, ] * repeat_precipitation[time_step, ]) + 
      (repeat_a2[time_step, ] * (previous_boxcox_streamflow - boxcox_streamflow_average)) +
      (repeat_a3_slope[time_step, ] * repeat_precipitation[time_step, ] * repeat_CO2[time_step, ]) + 
      (repeat_a4[time_step, ] * repeat_seasonal_ratio[time_step, ])
    
    previous_boxcox_streamflow <- current_boxcox_streamflow
    modelled_boxcox_streamflow[time_step, ] <- current_boxcox_streamflow
  }
  
  return(modelled_boxcox_streamflow)
  
}

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



# 1. Build catchment_data object from data -------------------------------------
test_gauges <- best_model_per_gauge |>
  pull(gauge) 
  
catchment_data_list <- map(
  .x = test_gauges,
  .f = catchment_data_blueprint,
  observed_data = data,
  start_stop_indexes = start_stop_indexes
)



# 2. Prepare numerical_optimiser objects ---------------------------------------

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

# Make the character streamflow models into functions
streamflow_model_list_char <- best_model_per_gauge |> 
  pull(streamflow_model) |> 
  # replace CO2 with time to reference the adjusted models
  str_replace(pattern = "CO2", "time")

streamflow_model_list_fun <- map(
  .x = streamflow_model_list_char, # didn't use the right models
  .f = match.fun
)



numerical_optimiser_list <- map2(
  .x = catchment_data_list,
  .y = streamflow_model_list_fun,
  .f = numerical_optimiser_setup_vary_inputs,
  objective_function = constant_sd_objective_function, 
  streamflow_transform_method = log_sinh_transform, 
  bounds_and_transform_method = make_default_bounds_and_transform_methods, 
  minimise_likelihood = TRUE,
  streamflow_transform_method_offset = 0 
)



# Repeat it 10 times
numerical_optimiser_list <- rep(numerical_optimiser_list, each = 10)


# 3. Run through CMAES ---------------------------------------------------------

## Separate numerical_optimiser_list into chunks for parallel ==================
GAUGES_PER_CHUNK <- 16L
length_numerical_optimsers <- length(numerical_optimiser_list)
total_chunks <- ceiling(length_numerical_optimsers / GAUGES_PER_CHUNK)

chunk_index_for_split <- seq(from = 1, to = total_chunks) |> 
  rep(length.out = length(numerical_optimiser_list)) |> 
  sort() |> 
  as.factor()

chunked_numerical_optimiser_list <- split(numerical_optimiser_list, f = chunk_index_for_split)


## CMAES =======================================================================
run_and_save_CMAES <- function(numerical_optimiser_chunk, run, user_path) {
  
  # Run CMAES
  CMAES_results <- future_map(
    .x = numerical_optimiser_chunk,
    .f = CMAES,
    cmaes_control = list(), 
    print_monitor = FALSE,
    .options = furrr_options(seed = TRUE)
  )
  
  # Convert results to standardised format
  standardised_results <- map(
    .x = CMAES_results,
    .f = result_set
  ) |> 
  map(
    .f = parameters_summary
  ) |> 
  list_rbind(
    names_to = "run"
  ) 
  
  # Save results
  save_path = paste0(user_path, "/time_only_run_", run, ".csv")
  write_csv(standardised_results, file = save_path)
  
  # Wipe variables
  rm(CMAES_results, standardised_results)
  gc()
  
}


## Run =========================================================================
plan(multisession, workers = 16)
iwalk(
  .x = chunked_numerical_optimiser_list,
  .f = run_and_save_CMAES,
  user_path = "Modelling/Results/CMAES"
)


stop_here()
# Estimated time ~40 per chunk 
# 121 chunks
# time = 121 * 40 = 1 hour 20 minutes

# 4. Import files into a single table ------------------------------------------
combined_cmaes_parameters_file <- list.files(
  path = "./Modelling/Results/CMAES/",
  recursive = FALSE, # I don't want it looking in other folders
  pattern = "time_only_run",
  full.names = TRUE
) |>
  read_csv(show_col_types = FALSE)

best_run_cmaes_parameters_file <- combined_cmaes_parameters_file |> 
  slice_min(
    order_by = AIC,
    by = gauge
  ) |> 
  select(!run)

write_csv(best_run_cmaes_parameters_file, "Modelling/Results/CMAES/time_only_runs.csv")


# 5. Analyse the results -------------------------------------------------------

## Get time only runs
time_only_runs <- read_csv("Modelling/Results/CMAES/time_only_runs.csv") |> 
  select(gauge, AIC) |> 
  distinct() |> 
  rename(time_only_AIC = AIC)


evidence_ratio <- read_csv("Modelling/Results/CMAES/evidence_ratio_results.csv")
time_of_activation <- read_csv("Modelling/Results/all_time_of_emergence.csv") |> 
  select(gauge, year_time_of_emergence)




time_only_evidence_ratio <- best_model_per_gauge |> 
  select(!c(streamflow_model, contains_CO2)) |> 
  left_join(
    time_only_runs,
    by = join_by(gauge)
  ) |> 
  # rename for joining
  rename(
    CO2_model = AIC
  ) |> 
  # Join existing evidence ratio
  left_join(
    evidence_ratio,
    by = join_by(gauge, CO2_model)
  ) |> 
  # calculate evidence ratio of time only vs. CO2
  mutate(
    AIC_diff_time_vs_CO2 = time_only_AIC - CO2_model
  ) |> 
  mutate(
    evidence_ratio_time_vs_CO2 = exp(0.5 * AIC_diff_time_vs_CO2)
  ) |> 
  mutate(
    plus_minus_evidence_ratio_time_vs_CO2 = case_when(
      AIC_diff_time_vs_CO2 >= 0 ~ exp(0.5 * abs(AIC_diff_time_vs_CO2)), # when CO2 model is better
      AIC_diff_time_vs_CO2 < 0 ~ -exp(0.5 * abs(AIC_diff_time_vs_CO2)) # when non-CO2 model is better
    )
  ) |> 
  # add time of activation
  left_join(
    time_of_activation,
    by = join_by(gauge)
  ) |> 
  arrange(evidence_ratio_time_vs_CO2)
  




## change shape of moderately strong or greater
ab_tibble <- tribble(
  ~year_time_of_emergence, ~evidence_ratio_time_vs_CO2, ~label_name, ~high_evidence_ratio,
  1960, 1E7, "a", "Less than Moderately Strong",
  1960, 1E7, "b", "Moderately Strong or Greater"
)

## Plot
time_only_evidence_ratio |> 
  mutate(
    high_evidence_ratio = if_else(evidence_ratio > 100, "Moderately Strong or Greater", "Less than Moderately Strong")
  ) |> 
  ggplot(aes(x = year_time_of_emergence, y = evidence_ratio_time_vs_CO2)) +
  geom_point(
    fill = "grey",
    colour = "black",
    stroke = 0.2,
    size = 3,
    shape = 21
  ) +
  geom_text(
    data = ab_tibble,
    aes(label = label_name),
    fontface = "bold"
  ) +
  geom_hline(yintercept = c(0.1, 10), colour = "red", linetype = "dashed") +
  scale_y_log10(breaks = (10^seq(from = -10, to = 9, by = 1))) +
  scale_x_continuous(n.breaks = 8) +
  labs(
    x = "Time of Activation",
    y = bquote("Evidence Ratio (Linear Trend vs."~CO[2]~"Model)") # big number CO2 is better
  ) +
  theme_bw() +
  facet_wrap(~high_evidence_ratio, scales = "fixed")

start_and_end <- data |> 
  drop_na() |> 
  summarise(
    max_year = max(year),
    min_year = min(year),
    .by = gauge
  ) |> 
  right_join(
    time_of_activation,
    by = join_by(gauge)
  ) |> 
  mutate(
    near_bounds_threshold = ceiling((max_year - min_year) * 0.1) # is the ToA in final or last 10 % of t-series
  ) |> 
  mutate(
    near_start_end = near(x = max_year, y = year_time_of_emergence, tol = near_bounds_threshold) | near(x = min_year, y = year_time_of_emergence, tol = near_bounds_threshold)
  ) 



x <- time_only_evidence_ratio |> 
  left_join(
    start_and_end,
    by = join_by(gauge)
  ) #|> 
  #filter(evidence_ratio > 100) 

x |> 
  filter(evidence_ratio_time_vs_CO2 > 100)

# they roughly have the same shape distribution
wilcox.test(
  x = x |> filter(near_start_end) |> pull(evidence_ratio_time_vs_CO2),
  y = x |> filter(!near_start_end) |> pull(evidence_ratio_time_vs_CO2),
  alternative = "two.sided"
)

x |> 
  summarise(
    median_evidence_ratio_time_vs_CO2 = median(plus_minus_evidence_ratio_time_vs_CO2),
    mean_evidence_ratio_time_vs_CO2 = mean(plus_minus_evidence_ratio_time_vs_CO2),
    sd_evidence_ratio_time_vs_CO2 = sd(plus_minus_evidence_ratio_time_vs_CO2),
    .by = near_start_end
  )
  #filter(evidence_ratio_time_vs_CO2 > 10) # > 10 = 27 (33%), > 100 = 19 (23 %)


ggsave(
  filename = "time_evidence_ratio_vs_ToA.pdf",
  plot = get_last_plot(),
  device = "pdf",
  path = "Figures/Extended_Data/",
  width = 180,
  height = 140,
  units = "mm"
)




stop_here()
# SINGLE TEST CATCHMENT
gauge <- "607220"   


example_catchment <- gauge |> 
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  ) 

numerical_optimiser <- example_catchment |>
  numerical_optimiser_setup_vary_inputs(
    streamflow_model = streamflow_model_intercept_shifted_time_seasonal_ratio_auto,
    objective_function = constant_sd_objective_function, 
    streamflow_transform_method = log_sinh_transform, 
    bounds_and_transform_method = make_default_bounds_and_transform_methods, 
    minimise_likelihood = TRUE,
    streamflow_transform_method_offset = 0 
  )


results <- numerical_optimiser |>
  CMAES(
    cmaes_control = list(), # can alter default cmaes controls here - see cmaesr package for options
    print_monitor = TRUE # print optimiser outputs in console while CMAES is running
  ) 


standardised_results <- results |> result_set() 

plot(standardised_results, type = "streamflow-time")
plot(standardised_results, type = "rainfall-runoff")
plot(standardised_results, type = "examine_transform") 
plot(standardised_results, type = "std_errors")


parameters_summary(standardised_results) 
 




