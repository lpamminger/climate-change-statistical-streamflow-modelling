# Non-drought models -----------------------------------------------------------
streamflow_model_precip_only <- function(catchment_data, parameter_set) { 
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_precip_only",
        "parameters" = c("a0", "a1")
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
  
  # Get data
  precipitation <- catchment_data$precipitation
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) #matrix(rep(precipitation, times = ncol(parameter_set)), ncol = ncol(parameter_set), byrow = FALSE)
  repeat_a0 <- matrix(a0, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Model
  repeat_a0 + (repeat_a1 * repeat_precipitation)
  
}





streamflow_model_separate_CO2 <- function(catchment_data, parameter_set) { 
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_separate_CO2",
        "parameters" = c("a0", "a1", "a3")
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
  
  # Get data
  precipitation <- catchment_data$precipitation
  CO2 <- catchment_data$CO2
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) 
  repeat_CO2 <- matrix(CO2, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) 
  repeat_a0 <- matrix(a0, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a3 <- matrix(a3, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Model
  repeat_a0 + (repeat_a1 * repeat_precipitation) + (repeat_a3 * repeat_CO2)
  
}





streamflow_model_precip_seasonal_ratio <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_precip_seasonal_ratio",
        "parameters" = c("a0", "a1", "a4")
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
  a4 <- parameter_set[3, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  seasonal_ratio <- catchment_data$seasonal_ratio
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) #matrix(rep(precipitation, times = ncol(parameter_set)), ncol = ncol(parameter_set), byrow = FALSE)
  repeat_seasonal_ratio <- matrix(seasonal_ratio, ncol = ncol(parameter_set), nrow = length(seasonal_ratio), byrow = FALSE)
  repeat_a0 <- matrix(a0, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a4 <- matrix(a4, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Model
  repeat_a0 + (repeat_a1 * repeat_precipitation) + (repeat_a4 * repeat_seasonal_ratio)
  
}





streamflow_model_separate_CO2_seasonal_ratio <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_separate_CO2_seasonal_ratio",
        "parameters" = c("a0", "a1", "a3", "a4")
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
  
  # Model
  repeat_a0 + (repeat_a1 * repeat_precipitation) + (repeat_a3 * repeat_CO2) + (repeat_a4 * repeat_seasonal_ratio)
  
}




streamflow_model_precip_auto <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_precip_auto",
        "parameters" = c("a0", "a1", "a2")
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
  
  # Get data
  precipitation <- catchment_data$precipitation
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) #matrix(rep(precipitation, times = ncol(parameter_set)), ncol = ncol(parameter_set), byrow = FALSE)
  repeat_a0 <- matrix(a0, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a2 <- matrix(a2, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Autocorrelation requirements
  boxcox_streamflow_average <- colMeans(repeat_a0 + (repeat_a1 * repeat_precipitation))
  previous_boxcox_streamflow <- boxcox_streamflow_average # assume the previous flows are the average. this means the first a2 component = 0.
  
  # pre-allocate
  modelled_boxcox_streamflow <- matrix(numeric(length(repeat_precipitation)), ncol = ncol(parameter_set), nrow = length(precipitation))
  
  # Model: analysis is done row by row (timestep)
  for (time_step in seq_along(precipitation)) {
    current_boxcox_streamflow <- repeat_a0[time_step, ] + (repeat_a1[time_step, ] * repeat_precipitation[time_step, ]) + (repeat_a2[time_step, ] * (previous_boxcox_streamflow - boxcox_streamflow_average))
    previous_boxcox_streamflow <- current_boxcox_streamflow
    modelled_boxcox_streamflow[time_step, ] <- current_boxcox_streamflow
  }
  
  return(modelled_boxcox_streamflow)
}




streamflow_model_separate_CO2_auto <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_separate_CO2_auto",
        "parameters" = c("a0", "a1", "a2", "a3")
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
  
  # Get data
  precipitation <- catchment_data$precipitation
  CO2 <- catchment_data$CO2
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) #matrix(rep(precipitation, times = ncol(parameter_set)), ncol = ncol(parameter_set), byrow = FALSE)
  repeat_CO2 <- matrix(CO2, ncol = ncol(parameter_set), nrow = length(CO2), byrow = FALSE) 
  repeat_a0 <- matrix(a0, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a2 <- matrix(a2, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a3 <- matrix(a3, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
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
      (repeat_a3[time_step, ] * repeat_CO2[time_step, ])
    
    previous_boxcox_streamflow <- current_boxcox_streamflow
    modelled_boxcox_streamflow[time_step, ] <- current_boxcox_streamflow
  }
  
  return(modelled_boxcox_streamflow)
}




streamflow_model_precip_seasonal_ratio_auto <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_precip_seasonal_ratio_auto",
        "parameters" = c("a0", "a1", "a2", "a4")
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
  a4 <- parameter_set[4, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  seasonal_ratio <- catchment_data$seasonal_ratio
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) #matrix(rep(precipitation, times = ncol(parameter_set)), ncol = ncol(parameter_set), byrow = FALSE)
  repeat_seasonal_ratio <- matrix(seasonal_ratio, ncol = ncol(parameter_set), nrow = length(seasonal_ratio), byrow = FALSE)
  repeat_a0 <- matrix(a0, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a2 <- matrix(a2, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
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
      (repeat_a4[time_step, ] * repeat_seasonal_ratio[time_step, ])
    
    previous_boxcox_streamflow <- current_boxcox_streamflow
    modelled_boxcox_streamflow[time_step, ] <- current_boxcox_streamflow
  }
  
  return(modelled_boxcox_streamflow)
  
}




streamflow_model_separate_CO2_seasonal_ratio_auto <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_separate_CO2_seasonal_ratio_auto",
        "parameters" = c("a0", "a1", "a2", "a3", "a4")
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
      (repeat_a3[time_step, ] * repeat_CO2[time_step, ]) + 
      (repeat_a4[time_step, ] * repeat_seasonal_ratio[time_step, ])
    
    previous_boxcox_streamflow <- current_boxcox_streamflow
    modelled_boxcox_streamflow[time_step, ] <- current_boxcox_streamflow
  }
  
  return(modelled_boxcox_streamflow)
  
}








# Drought streamflow models ----------------------------------------------------
streamflow_model_drought_precip_only <- function(catchment_data, parameter_set) { 
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_drought_precip_only",
        "parameters" = c("a0_d", "a0_n", "a1")
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
  
  # Get data
  precipitation <- catchment_data$precipitation
  drought <- catchment_data$is_drought_year
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) #matrix(rep(precipitation, times = ncol(parameter_set)), ncol = ncol(parameter_set), byrow = FALSE)
  repeat_drought <- matrix(drought, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE)
  repeat_a0_d <- matrix(a0_d, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a0_n <- matrix(a0_n, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Model
  (repeat_a0_d * repeat_drought) + (repeat_a0_n * !repeat_drought) + (repeat_a1 * repeat_precipitation)
  
}



streamflow_model_drought_separate_CO2 <- function(catchment_data, parameter_set) { 
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_drought_separate_CO2",
        "parameters" = c("a0_d", "a0_n", "a1", "a3")
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
  a3 <- parameter_set[4, ]
  
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
  repeat_a3 <- matrix(a3, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Model
  (repeat_a0_d * repeat_drought) + (repeat_a0_n * !repeat_drought) + (repeat_a1 * repeat_precipitation) + (repeat_a3 * repeat_CO2)
  
}




streamflow_model_drought_precip_seasonal_ratio <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_drought_precip_seasonal_ratio",
        "parameters" = c("a0_d", "a0_n", "a1", "a4")
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
  a4 <- parameter_set[4, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  seasonal_ratio <- catchment_data$seasonal_ratio
  drought <- catchment_data$is_drought_year
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) #matrix(rep(precipitation, times = ncol(parameter_set)), ncol = ncol(parameter_set), byrow = FALSE)
  repeat_seasonal_ratio <- matrix(seasonal_ratio, ncol = ncol(parameter_set), nrow = length(seasonal_ratio), byrow = FALSE)
  repeat_drought <- matrix(drought, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE)
  repeat_a0_d <- matrix(a0_d, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a0_n <- matrix(a0_n, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a4 <- matrix(a4, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Model
  (repeat_a0_d * repeat_drought) + (repeat_a0_n * !repeat_drought) + (repeat_a1 * repeat_precipitation) + (repeat_a4 * repeat_seasonal_ratio)
  
}


streamflow_model_drought_separate_CO2_seasonal_ratio <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_drought_separate_CO2_seasonal_ratio",
        "parameters" = c("a0_d", "a0_n", "a1", "a3", "a4")
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
  a3 <- parameter_set[4, ]
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
  repeat_a3 <- matrix(a3, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a4 <- matrix(a4, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Model
  (repeat_a0_d * repeat_drought) + (repeat_a0_n * !repeat_drought) + (repeat_a1 * repeat_precipitation) + (repeat_a3 * repeat_CO2) + (repeat_a4 * repeat_seasonal_ratio)
  
}





streamflow_model_drought_precip_auto <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_drought_precip_auto",
        "parameters" = c("a0_d", "a0_n", "a1", "a2")
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
  
  # Get data
  precipitation <- catchment_data$precipitation
  drought <- catchment_data$is_drought_year
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) #matrix(rep(precipitation, times = ncol(parameter_set)), ncol = ncol(parameter_set), byrow = FALSE)
  repeat_drought <- matrix(drought, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE)
  repeat_a0_d <- matrix(a0_d, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a0_n <- matrix(a0_n, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a2 <- matrix(a2, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Autocorrelation requirements
  boxcox_streamflow_average <- colMeans((repeat_a0_d * repeat_drought) + (repeat_a0_n * !repeat_drought) + (repeat_a1 * repeat_precipitation))
  previous_boxcox_streamflow <- boxcox_streamflow_average # assume the previous flows are the average. this means the first a2 component = 0.
  
  # pre-allocate
  modelled_boxcox_streamflow <- matrix(numeric(length(repeat_precipitation)), ncol = ncol(parameter_set), nrow = length(precipitation))
  
  # Model: analysis is done row by row (timestep)
  for (time_step in seq_along(precipitation)) {
    current_boxcox_streamflow <- (repeat_a0_d[time_step, ] * repeat_drought[time_step, ]) + 
                                 (repeat_a0_n[time_step, ] * !repeat_drought[time_step, ]) + 
                                 (repeat_a1[time_step, ] * repeat_precipitation[time_step, ]) + 
                                 (repeat_a2[time_step, ] * (previous_boxcox_streamflow - boxcox_streamflow_average))
    
    previous_boxcox_streamflow <- current_boxcox_streamflow
    modelled_boxcox_streamflow[time_step, ] <- current_boxcox_streamflow
  }
  
  return(modelled_boxcox_streamflow)
}




streamflow_model_drought_separate_CO2_auto <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_drought_separate_CO2_auto",
        "parameters" = c("a0_d", "a0_n", "a1", "a2", "a3")
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
  repeat_a3 <- matrix(a3, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
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
                                 (repeat_a3[time_step, ] * repeat_CO2[time_step, ])
    
    previous_boxcox_streamflow <- current_boxcox_streamflow
    modelled_boxcox_streamflow[time_step, ] <- current_boxcox_streamflow
  }
  
  return(modelled_boxcox_streamflow)
}



streamflow_model_drought_precip_seasonal_ratio_auto <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_drought_precip_seasonal_ratio_auto",
        "parameters" = c("a0_d", "a0_n", "a1", "a2", "a4")
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
  a4 <- parameter_set[5, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  seasonal_ratio <- catchment_data$seasonal_ratio
  drought <- catchment_data$is_drought_year
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) #matrix(rep(precipitation, times = ncol(parameter_set)), ncol = ncol(parameter_set), byrow = FALSE)
  repeat_seasonal_ratio <- matrix(seasonal_ratio, ncol = ncol(parameter_set), nrow = length(seasonal_ratio), byrow = FALSE)
  repeat_drought <- matrix(drought, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE)
  repeat_a0_d <- matrix(a0_d, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a0_n <- matrix(a0_n, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a2 <- matrix(a2, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a4 <- matrix(a4, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Autocorrelation requirements
  boxcox_streamflow_average <- colMeans((repeat_a0_d * repeat_drought) + (repeat_a0_n * !repeat_drought) + (repeat_a1 * repeat_precipitation))
  previous_boxcox_streamflow <- boxcox_streamflow_average # assume the previous flows are the average. this means the first a2 component = 0.
  
  # pre-allocate
  modelled_boxcox_streamflow <- matrix(numeric(length(repeat_precipitation)), ncol = ncol(parameter_set), nrow = length(precipitation))
  
  # Model: analysis is done row by row (timestep)
  for (time_step in seq_along(precipitation)) {
    current_boxcox_streamflow <- (repeat_a0_d[time_step, ] * repeat_drought[time_step, ]) + 
                                 (repeat_a0_n[time_step, ] * !repeat_drought[time_step, ]) +
                                 (repeat_a1[time_step, ] * repeat_precipitation[time_step, ]) + 
                                 (repeat_a2[time_step, ] * (previous_boxcox_streamflow - boxcox_streamflow_average)) +
                                 (repeat_a4[time_step, ] * repeat_seasonal_ratio[time_step, ])
    
    previous_boxcox_streamflow <- current_boxcox_streamflow
    modelled_boxcox_streamflow[time_step, ] <- current_boxcox_streamflow
  }
  
  return(modelled_boxcox_streamflow)
  
}



streamflow_model_drought_separate_CO2_seasonal_ratio_auto <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_drought_separate_CO2_seasonal_ratio_auto",
        "parameters" = c("a0_d", "a0_n", "a1", "a2", "a3", "a4")
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
                                 (repeat_a3[time_step, ] * repeat_CO2[time_step, ]) + 
                                 (repeat_a4[time_step, ] * repeat_seasonal_ratio[time_step, ])
    
    previous_boxcox_streamflow <- current_boxcox_streamflow
    modelled_boxcox_streamflow[time_step, ] <- current_boxcox_streamflow
  }
  
  return(modelled_boxcox_streamflow)
  
}


# Get streamflow_models function -----------------------------------------------
get_drought_streamflow_models <- function() {
  c(streamflow_model_drought_precip_only,
    streamflow_model_drought_precip_auto,
    streamflow_model_drought_separate_CO2,
    streamflow_model_drought_separate_CO2_auto,
    streamflow_model_drought_separate_CO2_seasonal_ratio,
    streamflow_model_drought_precip_seasonal_ratio,
    streamflow_model_drought_precip_seasonal_ratio_auto,
    streamflow_model_drought_separate_CO2_seasonal_ratio_auto
  )
}


get_non_drought_streamflow_models <- function() {
  c(streamflow_model_precip_only,
    streamflow_model_precip_auto,
    streamflow_model_separate_CO2,
    streamflow_model_separate_CO2_auto,
    streamflow_model_separate_CO2_seasonal_ratio,
    streamflow_model_precip_seasonal_ratio,
    streamflow_model_precip_seasonal_ratio_auto,
    streamflow_model_separate_CO2_seasonal_ratio_auto
  )
}


