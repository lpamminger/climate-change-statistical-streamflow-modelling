# Boxcox transforms ------------------------------------------------------------

boxcox_transform <- function(y, lambda = 0, lambda_2 = 0) {

  transformed_result <- ((y + lambda_2)^lambda - 1) / lambda
  
  # to make it work with log
  # get true/false vector of lambda values if near zero
  # if true replace with log(...)
  near_zero_lambda <- lambda <= .Machine$double.eps^0.5 #dplyr::near(lambda, y = 0, tol = .Machine$double.eps^0.5)
  
  transformed_result[near_zero_lambda] <- log(y[near_zero_lambda] + lambda_2)
  return(transformed_result)
 

}



boxcox_inverse_transform <- function(yt, lambda = 0, lambda_2 = 0) {
  realspace_result <- ((yt * lambda) + 1)^(1 / lambda) - lambda_2
  
  # to make it work with log
  # get true/false vector of lambda values if near zero
  # if true replace with log(...)
  near_zero_lambda <- lambda <= .Machine$double.eps^0.5#dplyr::near(lambda, y = 0, tol = .Machine$double.eps^0.5)
  
  realspace_result[near_zero_lambda] <- exp(yt[near_zero_lambda]) - lambda_2
  
  return(realspace_result)
}





# Log-sinh transform -----------------------------------------------------------
log_sinh_transform <- function(a, b, y, offset = 300) {
  (1 / b) * log(sinh(a + (b * (y + offset))))
}

asinh_exp_approximation <- function(x) {
  # Input: asinh(exp(x))
  
  #Testing tim's log-sinh - add to inverse_log_sinh then delete
  #x <- 0.1:1000
  #y <- asinh(exp(x))
  #y2 <- log(exp(x) + sqrt(exp(2 * x) + 1))
  #y3 <- x + log(1 + sqrt(exp(2 * x) + 1) / exp(x))
  #y4 <- x + log(1 + sqrt(1 + (1 / exp(2 * x))))
  #y5 <- (is.infinite(y4) * y4) + (is.infinite(y4) * (x + log(2))) - I am not sure what this is for
  
  x + log(1 + sqrt(1 + (1 / exp(2 * x))))
  
}



inverse_log_sinh_transform <- function(a, b, z, offset = 300) {
  # If any streamflow value is infinite the approximation is done for everything
  if (any(is.infinite(exp(b * z)))) {
    asinh_component <- asinh_exp_approximation(b * z)
  } else {
    asinh_component <- asinh(exp(b * z))
  }
  
  (asinh_component / b) - (a / b) - offset
}



