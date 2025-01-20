make_parameter_names <- function(streamflow_model, objective_function) {
  c(streamflow_model()[[2]], objective_function()[[2]]) # [[2]] extracts parameters from the list
  # streamflow_model() and objective_function() must return character atomic of parameters
}




# Constructor ------------------------------------------------------------------
## Builds object to be passed to either CMAES or DREAM
## I have limited the class to only have one model, objective_function and gauge

new_numerical_optimiser_setup <- function(streamflow_model, transform_parameter_methods, lower_bound, upper_bound, objective_function, catchment_data, scale, minimise_likelihood) {
  
  ## Check types ===============================================================
  stopifnot(sloop::s3_class(catchment_data)[1] == "catchment_data")
  stopifnot(is.function(streamflow_model)) #if putting list of functions in stopifnot(all(map_lgl(.x = streamflow_model, is.function)))
  stopifnot(is.function(objective_function)) # leave this for now. It would be good if I could have multiple objective_functions
  stopifnot(is.logical(minimise_likelihood))
  stopifnot(is.double(lower_bound))
  stopifnot(is.double(upper_bound))
  stopifnot(is.list(transform_parameter_methods))
  stopifnot(is.double(scale))
  stopifnot(all(purrr::map_lgl(.x = transform_parameter_methods, is.function)))
  
  
  ## Objective function ========================================================
  ready_objective_function <- objective_function_setup(
    streamflow_model = streamflow_model,
    transform_parameter_methods = transform_parameter_methods,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    objective_function = objective_function,
    catchment_data = catchment_data, 
    scale = scale,
    minimise_likelihood = minimise_likelihood
  )
  
  
  ## Make list of parameters ===================================================
  parameter_names <- make_parameter_names(streamflow_model, objective_function)
  
  ## Make class ================================================================
  structure(
    list(
      "catchment_data" = catchment_data,
      "streamflow_model" = streamflow_model,
      "objective_function" = objective_function,
      "transform_parameter_methods" = transform_parameter_methods,
      "lower_bound" = lower_bound,
      "upper_bound" = upper_bound, 
      "scale" = scale,
      "minimise_likelihood" = minimise_likelihood,
      "ready_objective_function" = ready_objective_function,
      "parameter_names" = parameter_names
    ),
    class = c("numerical_optimiser_setup", "list")
  )
  
}


# Validator --------------------------------------------------------------------
validate_optimiser_set <- function(numerical_optimiser_setup) {
  
  values <- unclass(numerical_optimiser_setup)
  
  # Checks still to be implemented:
  ##  putting non-drought catchment in a drought model should throw an error - put in optimiser class
  
  ## Check length of transform parameter method and lower/upper bound ==========
  check_length_values <- values[c("transform_parameter_methods", "lower_bound", "upper_bound", "parameter_names")]
  
  if (length(unique(lengths(check_length_values))) != 1) {
    stop(
      "transform_parameter_methods, lower_bound, upper_bound all must be the same length",
      call. = FALSE
    )
  }
  
  
  
  ## Lower bound must be less then the upper bound =============================
  if(any(values$lower_bound >= values$upper_bound)) {
    stop(
      "lower_bound must be less then upper_bound",
      call. = FALSE
    )
  }
  
  
  
  # If using log-transform cannot have bound less then zero ===================
  logarithmic_parameter_transform_index <- names(values$transform_parameter_methods) == "logarithmic_parameter_transform"
  
  if(any(values$lower_bound[logarithmic_parameter_transform_index] <= 0)) {
    stop(
      "logarithmic_parameter_transform requires a lower bound greater than zero",
      call. = FALSE
    )
  }
  
  
  # Catchments without a drought cannot be used with drought models ============
  if(any(values$streamflow_model()[[2]] %in% c("a0_d", "a0_n")) & (values$catchment_data$contains_drought == FALSE)) { #values$streamflow_model()[[2] extracts parameter from the list
    stop(
      "Cannot use a drought streamflow model for a catchment without a drought",
      call. = FALSE
    )
  } 
  
  return(numerical_optimiser_setup)
  
}




# Helper -----------------------------------------------------------------------
# give the user suggested pre-sets?
# i.e., recommended transform methods based on streamflow model and optimiser combinations, single bounds tibble

numerical_optimiser_setup <- function(streamflow_model, objective_function, catchment_data, bounds_and_transform_method, scale = 100, minimise_likelihood) { # same as constructor
  
  parameter_names <- make_parameter_names(streamflow_model, objective_function)
  
  selected_defaults <- bounds_and_transform_method |> 
    dplyr::filter(parameter %in% parameter_names)
  
  
  ## Override defaults if NULL is replaced =====================================
  ### Add later...
  
  ## Make and validate optimiser_set class =====================================
  new_numerical_optimiser_setup(
    streamflow_model = streamflow_model,
    transform_parameter_methods = selected_defaults |> dplyr::pull(transform_method),
    lower_bound = selected_defaults |> dplyr::pull(lower_bound),
    upper_bound = selected_defaults |> dplyr::pull(upper_bound),
    objective_function = objective_function,
    catchment_data = catchment_data,
    scale = scale,
    minimise_likelihood = minimise_likelihood
  ) |> 
    validate_optimiser_set()
  
}


numerical_optimiser_setup_vary_inputs <- function(catchment_data, ...) {
  numerical_optimiser_setup(catchment_data = catchment_data, ...)
}


# Make default bounds and transfer methods -------------------------------------
make_default_bounds_and_transform_methods <- function() {
  tibble::tribble(
    ~parameter, ~lower_bound,   ~upper_bound,  ~transform_method,
    "a0",        -300,           100,           linear_parameter_transform, # intercept
    "a0_d",      -300,           50,            linear_parameter_transform, # intercept - no drought
    "a0_n",      -300,           50,            linear_parameter_transform, # intercept - drought
    "a1",         1E-5 ,         1,             logarithmic_parameter_transform, # slope
    "a2",        -1,             1,             linear_parameter_transform, # autocorrelation
    "a3",        -25,            50,            linear_parameter_transform, # CO2 coefficient 
    "a4",        -250,           600,           linear_parameter_transform, # seasonal parameter
    "a5",         0,             150,           linear_parameter_transform, # 97.70 is CO2 - 280 at 2004. CO2 shift parameter - 138.53 HARD CODED - SHOULD CHANGE BASED ON CO2 INPUT. not all gauges reach max CO2. Record ends earlier for some
    "sd",         1E-8,          200,           logarithmic_parameter_transform, # constant sd objective function 
    "scale_CO2",  1E-8,          2,             logarithmic_parameter_transform # CO2 scaler for objective function
  )
}


