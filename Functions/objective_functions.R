# Alpha checking ---------------------------------------------------------------
# lower truncation term adjusting
# Alpha using in pnrom. If -mu/sigma is large positve number then pnorm(alpha) = 1
# This causes a divide by zero warning. To avoid warning just set to NA.
# Let NA propogate without warnings
lower_bound_correction <- function(uncorrected_mean_flow, uncorrected_uncertainty) {
  alpha <- (0 - uncorrected_mean_flow) / uncorrected_uncertainty

  # cut_off_constant <- qnorm(p = 1 - .Machine$double.eps^0.5) # this is the limit using the inverse cdf alpha = max(alpha, qnorm(p = ...))

  # alpha[alpha >= cut_off_constant] <- cut_off_constant

  near_machine_precision <- dplyr::near(pnorm(alpha), y = 1, tol = .Machine$double.eps^0.5) # check if within machine precision

  alpha[near_machine_precision] <- NA # Set values that are near 1 to NA to avoid divide by zero warnings

  return(alpha)
}




# Correcting modelled streamflow and standard deviation due to trun norm -------
# https://en.wikipedia.org/wiki/Truncated_normal_distribution
correct_mean_flow <- function(uncorrected_mean_flow, uncorrected_uncertainty) {
  alpha <- lower_bound_correction(
    uncorrected_mean_flow = uncorrected_mean_flow,
    uncorrected_uncertainty = uncorrected_uncertainty
  )

  uncorrected_mean_flow + ((uncorrected_uncertainty * dnorm(alpha)) / (1 - pnorm(alpha)))
}


correct_uncertainty_flow <- function(uncorrected_mean_flow, uncorrected_uncertainty) {
  alpha <- lower_bound_correction(
    uncorrected_mean_flow = uncorrected_mean_flow,
    uncorrected_uncertainty = uncorrected_uncertainty
  )

  sqrt(uncorrected_uncertainty^2 * (1 + ((alpha * dnorm(alpha)) / (1 - pnorm(alpha))) - (dnorm(alpha) / (1 - pnorm(alpha)))^2))
}






# Objective functions ----------------------------------------------------------
# put observed streamflow transformation in the objective function

constant_sd_no_transform_objective_function <- function(modelled_streamflow, observed_streamflow, parameter_set) {
  if (is.null(names(as.list(match.call())[-1]))) { # if no arguments provided return description
    return(
      list(
        "name" = "constant_sd_objective_function",
        "parameters" = "sd"
      )
    )
  }

  # browser()
  # dtruncnorm only works with vectors (double)
  # matrix is a subclass of double and gets coerced into a double (double atomic vector)
  constant_sd <- parameter_set[nrow(parameter_set), ]

  matrix_error_sd <- matrix(constant_sd,
    nrow = nrow(modelled_streamflow),
    ncol = ncol(modelled_streamflow),
    byrow = TRUE
  )

  # Correct modelled streamflow and uncertainty
  corrected_modelled_streamflow <- correct_mean_flow(
    uncorrected_mean_flow = modelled_streamflow,
    uncorrected_uncertainty = matrix_error_sd
  )

  corrected_uncertainty <- correct_uncertainty_flow(
    uncorrected_mean_flow = modelled_streamflow,
    uncorrected_uncertainty = matrix_error_sd
  )


  # Produce probabilities using trunnorm
  prob_boxcox_observed <- truncnorm::dtruncnorm(
    x = observed_streamflow,
    a = 0,
    b = Inf,
    mean = corrected_modelled_streamflow,
    sd = corrected_uncertainty
  )

  # Convert vector back into a matrix
  prob_boxcox_observed <- matrix(prob_boxcox_observed,
    nrow = nrow(modelled_streamflow),
    ncol = ncol(modelled_streamflow)
  )

  negative_log_likelihood <- colSums(-1 * log(prob_boxcox_observed))

  return(negative_log_likelihood)
}






constant_sd_log_sinh_objective_function <- function(modelled_streamflow, observed_streamflow, parameter_set) {
  if (is.null(names(as.list(match.call())[-1]))) { # if no arguments provided return description
    return(
      list(
        "name" = "constant_sd_log_sinh_objective_function",
        "parameters" = c("sd", "a", "b")
      )
    )
  }


  constant_sd <- parameter_set[(nrow(parameter_set) - 2), ]
  log_sinh_a <- parameter_set[(nrow(parameter_set) - 1), ]
  log_sinh_b <- parameter_set[nrow(parameter_set), ]

  matrix_error_sd <- matrix(
    constant_sd,
    nrow = nrow(modelled_streamflow),
    ncol = ncol(modelled_streamflow),
    byrow = TRUE
  )

  matrix_log_sinh_a <- matrix(
    log_sinh_a,
    nrow = nrow(modelled_streamflow),
    ncol = ncol(modelled_streamflow),
    byrow = TRUE
  )


  matrix_log_sinh_b <- matrix(
    log_sinh_b,
    nrow = nrow(modelled_streamflow),
    ncol = ncol(modelled_streamflow),
    byrow = TRUE
  )


  # Correct modelled streamflow and uncertainty
  corrected_modelled_streamflow <- correct_mean_flow(
    uncorrected_mean_flow = modelled_streamflow,
    uncorrected_uncertainty = matrix_error_sd
  )

  corrected_uncertainty <- correct_uncertainty_flow(
    uncorrected_mean_flow = modelled_streamflow,
    uncorrected_uncertainty = matrix_error_sd
  )


  # Transform observed_streamflow into log-sinh space
  transformed_observed_streamflow <- log_sinh_transform(
    a = matrix_log_sinh_a,
    b = matrix_log_sinh_b,
    y = observed_streamflow,
    offset = 300
  )


  # Produce probabilities using trunnorm
  # dtruncnorm only works with vectors (double)
  # matrix is a subclass of double and gets coerced into a double (double atomic vector)
  # when matrix is coerced into double convert it back into matrix
  prob_boxcox_observed <- truncnorm::dtruncnorm(
    x = transformed_observed_streamflow,
    a = 0,
    b = Inf,
    mean = corrected_modelled_streamflow,
    sd = corrected_uncertainty
  )

  # Convert vector back into a matrix
  prob_boxcox_observed <- matrix(
    prob_boxcox_observed,
    nrow = nrow(modelled_streamflow),
    ncol = ncol(modelled_streamflow)
  )


  negative_log_likelihood <- colSums(-1 * log(prob_boxcox_observed))

  return(negative_log_likelihood)
}



constant_sd_boxcox_objective_function <- function(modelled_streamflow, observed_streamflow, parameter_set) {
  if (is.null(names(as.list(match.call())[-1]))) { # if no arguments provided return description
    return(
      list(
        "name" = "constant_sd_boxcox_objective_function",
        "parameters" = c("sd", "boxcox_lambda")
      )
    )
  }


  constant_sd <- parameter_set[(nrow(parameter_set) - 1), ]
  boxcox_lambda <- parameter_set[nrow(parameter_set), ]

  matrix_error_sd <- matrix(constant_sd,
    nrow = nrow(modelled_streamflow),
    ncol = ncol(modelled_streamflow),
    byrow = TRUE
  )

  matrix_boxcox_lambda <- matrix(boxcox_lambda,
    nrow = nrow(modelled_streamflow),
    ncol = ncol(modelled_streamflow),
    byrow = TRUE
  )


  # Correct modelled streamflow and uncertainty
  corrected_modelled_streamflow <- correct_mean_flow(
    uncorrected_mean_flow = modelled_streamflow,
    uncorrected_uncertainty = matrix_error_sd
  )

  corrected_uncertainty <- correct_uncertainty_flow(
    uncorrected_mean_flow = modelled_streamflow,
    uncorrected_uncertainty = matrix_error_sd
  )



  transformed_observed_streamflow <- boxcox_transform(
    y = observed_streamflow,
    lambda = matrix_boxcox_lambda,
    lambda_2 = 1
  )


  # Produce probabilities using trunnorm
  # dtruncnorm only works with vectors (double)
  # matrix is a subclass of double and gets coerced into a double (double atomic vector)
  # when matrix is coerced into double convert it back into matrix
  prob_boxcox_observed <- truncnorm::dtruncnorm(
    x = transformed_observed_streamflow,
    a = 0,
    b = Inf,
    mean = corrected_modelled_streamflow,
    sd = corrected_uncertainty
  )

  # Convert vector back into a matrix
  prob_boxcox_observed <- matrix(
    prob_boxcox_observed,
    nrow = nrow(modelled_streamflow),
    ncol = ncol(modelled_streamflow)
  )

  negative_log_likelihood <- colSums(-1 * log(prob_boxcox_observed))

  return(negative_log_likelihood)
}



# Get function -----------------------------------------------------------------
## Remove CO2_variable_objective_function from objective functions tested
get_all_objective_functions <- function() {
  c(constant_sd_objective_function)
}
