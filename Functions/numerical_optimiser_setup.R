make_parameter_names <- function(streamflow_model, streamflow_transform_method, objective_function) {
  c(streamflow_model()[[2]], streamflow_transform_method()[[2]], objective_function()[[2]]) # [[2]] extracts parameters from the list
  # streamflow_model() and objective_function() must return character atomic of parameters
}




# Constructor ------------------------------------------------------------------
## Builds object to be passed to either CMAES or DREAM
## I have limited the class to only have one model, objective_function and gauge

new_numerical_optimiser_setup <- function(streamflow_model, parameter_transform_method, streamflow_transform_method, lower_bound, upper_bound, objective_function, catchment_data, scale, minimise_likelihood, streamflow_transform_method_offset) {
  ## Check types ===============================================================
  stopifnot(sloop::s3_class(catchment_data)[1] == "catchment_data")
  stopifnot(is.function(streamflow_model)) # if putting list of functions in stopifnot(all(map_lgl(.x = streamflow_model, is.function)))
  stopifnot(is.function(objective_function)) # leave this for now. It would be good if I could have multiple objective_functions
  stopifnot(is.function(streamflow_transform_method))
  stopifnot(is.logical(minimise_likelihood))
  stopifnot(is.double(lower_bound))
  stopifnot(is.double(upper_bound))
  stopifnot(is.list(parameter_transform_method))
  stopifnot(is.double(scale))
  stopifnot(all(purrr::map_lgl(.x = parameter_transform_method, is.function)))
  stopifnot(is.double(streamflow_transform_method_offset))


  ## Objective function ========================================================
  ready_objective_function <- objective_function_setup(
    streamflow_model = streamflow_model,
    parameter_transform_method = parameter_transform_method,
    streamflow_transform_method = streamflow_transform_method,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    objective_function = objective_function,
    catchment_data = catchment_data,
    scale = scale,
    minimise_likelihood = minimise_likelihood,
    streamflow_transform_method_offset = streamflow_transform_method_offset
  )


  ## Make list of parameters ===================================================
  parameter_names <- make_parameter_names(streamflow_model, streamflow_transform_method, objective_function)

  ## Make class ================================================================
  structure(
    list(
      "catchment_data" = catchment_data,
      "streamflow_model" = streamflow_model,
      "objective_function" = objective_function,
      "parameter_transform_method" = parameter_transform_method,
      "streamflow_transform_method" = streamflow_transform_method,
      "lower_bound" = lower_bound,
      "upper_bound" = upper_bound,
      "scale" = scale,
      "minimise_likelihood" = minimise_likelihood,
      "ready_objective_function" = ready_objective_function,
      "parameter_names" = parameter_names,
      "streamflow_transform_method_offset" = streamflow_transform_method_offset
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
  check_length_values <- values[c("parameter_transform_method", "lower_bound", "upper_bound", "parameter_names")]

  if (length(unique(lengths(check_length_values))) != 1) {
    stop(
      "parameter_transform_method, lower_bound, upper_bound all must be the same length. Default method: make_default_bounds_and_transform_methods may not have parameter listed.",
      call. = FALSE
    )
  }



  ## Lower bound must be less then the upper bound =============================
  if (any(values$lower_bound >= values$upper_bound)) {
    stop(
      "lower_bound must be less then upper_bound",
      call. = FALSE
    )
  }



  # If using log-transform cannot have bound less then zero ===================
  logarithmic_parameter_transform_index <- names(values$parameter_transform_method) == "logarithmic_parameter_transform"

  if (any(values$lower_bound[logarithmic_parameter_transform_index] <= 0)) {
    stop(
      "logarithmic_parameter_transform requires a lower bound greater than zero",
      call. = FALSE
    )
  }


  # Catchments without a drought cannot be used with drought models ============
  if (any(values$streamflow_model()[[2]] %in% c("a0_d", "a0_n")) & (values$catchment_data$contains_drought == FALSE)) { # values$streamflow_model()[[2] extracts parameter from the list
    stop(
      "Cannot use a drought streamflow model for a catchment without a drought",
      call. = FALSE
    )
  }

  return(numerical_optimiser_setup)
}


# Use to combat zero flow values
alter_observed_streamflow <- function(observed_data_tibble, operation, value) {
  observed_data_tibble |>
    mutate(observed_streamflow = operation(observed_streamflow, value))
}


# Helper -----------------------------------------------------------------------
numerical_optimiser_setup <- function(
    streamflow_model, 
    objective_function, 
    streamflow_transform_method, 
    catchment_data, 
    bounds_and_transform_method, 
    streamflow_transform_method_offset, 
    scale = 100, 
    minimise_likelihood
    ) { # same as constructor

  ## Make a vector of parameter names ==========================================
  parameter_names <- make_parameter_names(streamflow_model, streamflow_transform_method, objective_function)


  ## log-sinh offset involves altering observed flow ===========================
  streamflow_transform_method <- noquote(streamflow_transform_method)

  if (streamflow_transform_method()$name == "log_sinh_transform") {
    full_data_replacement <- alter_observed_streamflow(
      observed_data_tibble = catchment_data$full_data_set,
      operation = `+`,
      value = streamflow_transform_method_offset
    )

    stop_start_data_replacement <- modify(
      .x = catchment_data$stop_start_data_set,
      .f = alter_observed_streamflow,
      operation = `+`,
      value = streamflow_transform_method_offset
    )


    catchment_data$full_data_set <- full_data_replacement
    catchment_data$stop_start_data_set <- stop_start_data_replacement
  }



  ## called bounds_and_transfer_method to make bounds ==========================
  bounds_and_transform_method <- noquote(bounds_and_transform_method)
  bounds_and_transform <- bounds_and_transform_method(catchment_data)

  selected_defaults <- bounds_and_transform |>
    dplyr::filter(parameter %in% parameter_names)



  ## Override defaults if NULL is replaced =====================================
  ### Add later...


  ## Make and validate optimiser_set class =====================================
  new_numerical_optimiser_setup(
    streamflow_model = streamflow_model,
    parameter_transform_method = selected_defaults |> dplyr::pull(transform_method),
    streamflow_transform_method = streamflow_transform_method,
    lower_bound = selected_defaults |> dplyr::pull(lower_bound),
    upper_bound = selected_defaults |> dplyr::pull(upper_bound),
    objective_function = objective_function,
    catchment_data = catchment_data,
    scale = scale,
    minimise_likelihood = minimise_likelihood,
    streamflow_transform_method_offset = streamflow_transform_method_offset
  ) |>
    validate_optimiser_set()
}

# Helpers helper - this is not great
numerical_optimiser_setup_vary_inputs <- function(catchment_data, ...) {
  numerical_optimiser_setup(catchment_data = catchment_data, ...)
}





# Make default bounds and transfer methods -------------------------------------
## 1. repeat this function but have catchment_data as input
## 2. find maximum CO2 value used in calibrate
## 3. replace upper bound for a5 with that value

log_sinh_asymptote <- function(a, b) {
  # Vertical asymptotes occur when log() approaches zero
  # sinh(a + bx) approaches zero
  # a + bx = 0
  # asymptote --> -a/b
  -a / b
}


## log-sinh needs special treatment to find acceptable bounds
find_acceptable_log_sinh_bounds <- function(streamflow) {
  # See log-sinh-bounds.xlsx for the method

  # rows and cols for matrix
  possible_a <- sort(c(-10^-seq(from = 1, to = 7, by = 1), -.Machine$double.eps^0.5))
  names_possible_a <- as.character(possible_a)
  length_a <- length(possible_a)
  possible_b <- sort(c(10^-seq(from = -3, to = 7, by = 1), -.Machine$double.eps^0.5))
  names_possible_b <- as.character(possible_b)
  length_b <- length(possible_b)


  matrix_of_values <- matrix(numeric(length_a * length_b), nrow = length_a, ncol = length_b) |>
    `colnames<-`(names_possible_b) |>
    `rownames<-`(names_possible_a)

  # min_streamflow <- min(streamflow)

  for (i in 1:length_b) {
    for (j in 1:length_a) {
      matrix_of_values[j, i] <- log_sinh_asymptote(
        a = possible_a[j],
        b = possible_b[i]
      )
    }
  }



  matrix_of_values[matrix_of_values > min(streamflow)] <- NA


  matrix_of_values[abs(matrix_of_values) < .Machine$double.eps^0.5] <- NA

  # range of b is always between 1E-3 and 1
  range_b <- c(1E-3, 10)

  # Use upper range_a to focus on row to get range_b
  range_a <- matrix_of_values[, as.character(range_b[1])] |>
    na.omit() |>
    names() |>
    as.numeric() |>
    range()




  list(
    "range_a" = range_a,
    "range_b" = range_b
  )
}






make_default_bounds_and_transform_methods <- function(catchment_data_set) {
  
  # Certain parameters depend on observed data ---------------------------------
  
  observed_streamflow <- catchment_data_set$stop_start_data_set |>
    list_rbind() |>
    pull(observed_streamflow)
  
  
  ## a5 (CO2 time of emergence) cannot be greater than maximum observed ========
  upper_a5_bound <- max(catchment_data_set$stop_start_data_set[[length(catchment_data_set$stop_start_data_set)]]$CO2)


  ## the intercept terms are related to the maximum observed streamflow ========
  max_observed_streamflow <- max(observed_streamflow) * 5

  ## many combination of log-sinh a and b - select acceptable ones =============
  log_sinh_bounds <- find_acceptable_log_sinh_bounds(observed_streamflow)



  tibble::tribble(
    ~parameter,    ~lower_bound,               ~upper_bound,               ~transform_method,
    "a0",           -max_observed_streamflow,   max_observed_streamflow,    linear_parameter_transform, # intercept
    "a0_d",         -max_observed_streamflow,   max_observed_streamflow,    linear_parameter_transform, # intercept - no drought
    "a0_n",         -max_observed_streamflow,   max_observed_streamflow,    linear_parameter_transform, # intercept - drought
    "a1",            1E-8,                      1,                          logarithmic_parameter_transform, # slope
    "a2",           -1,                         1,                          linear_parameter_transform, # autocorrelation
    "a3_intercept", -max_observed_streamflow,   max_observed_streamflow,    linear_parameter_transform, # CO2 coefficient for intercept
    "a3_slope",     -1,                         1,                          linear_parameter_transform, # CO2 coefficent for slope
    "a4",           -max_observed_streamflow,   max_observed_streamflow,    linear_parameter_transform, # seasonal parameter
    "a5",           0,                          upper_a5_bound,             linear_parameter_transform, # Changes depending on last CO2 value in calibration
    "a",            log_sinh_bounds$range_a[1], log_sinh_bounds$range_a[2], logarithmic_parameter_transform,
    "b",            log_sinh_bounds$range_b[1], log_sinh_bounds$range_b[2], logarithmic_parameter_transform,
    "lambda",       0,                          2,                          linear_parameter_transform, # boxcox recommended range is -2 to 2(0-2 because we want to shift it from positive skew)
    "sd",           1E-8,                       2000,                       logarithmic_parameter_transform # constant sd objective function
  )
}
