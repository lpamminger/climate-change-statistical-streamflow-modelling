linear_parameter_transform <- function(parameter_set, lower_bound, upper_bound, scale){
  lower_bound + ((upper_bound - lower_bound) * (parameter_set / scale)) 
}


logarithmic_parameter_transform <- function(parameter_set, lower_bound, upper_bound, scale){
  lower_bound * ((upper_bound / lower_bound)^(parameter_set / scale))
}


squared_parameter_transform <- function(parameter_set, lower_bound, upper_bound, scale){
  lower_bound + ((upper_bound - lower_bound) * (parameter_set / scale)^2)
}


no_parameter_transform <- function(parameter_set, lower_bound = NULL, upper_bound = NULL, scale = NULL) {
  parameter_set
}

# Transform parameters from 0-scale to real parameter values
transform_parameter_method <- function(method, idx, lower_bound, upper_bound, parameter_set, scale) {
  method <- noquote(method)
  
  # parameter_set is a matrix
  # cols = different combination of parameters and rows = parameter
  # only apply transform methods to rows (the same parameters)
  method(parameter_set = parameter_set[idx,], lower_bound = lower_bound, upper_bound = upper_bound, scale = scale)
}





# Converting real space into transformed space ---------------------------------
real_space_to_scaled_linear_transform <- function(real_space_parameter, lower_bound, upper_bound, scale) {
  (scale * (real_space_parameter - lower_bound)) / (upper_bound - lower_bound)
}

real_space_to_scaled_logarithmic_transform <- function(real_space_parameter, lower_bound, upper_bound, scale) {
  (scale * log(real_space_parameter / lower_bound)) / (log(upper_bound / lower_bound))
}
