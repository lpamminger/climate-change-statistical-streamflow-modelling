# Generic functions ------------------------------------------------------------


## get the best parameters from the optimiser (still in scaled form) ===========
coef.cma_result <- function(cma_result) {
  cmaes_coefs <- cma_result$best.param
  names(cmaes_coefs) <- names(cma_result$par.set$pars)
  return(cmaes_coefs)
}


## get the best loglikelihood ==================================================
get_best_fitness <- function(x) {
  UseMethod("get_best_fitness")
}

get_best_fitness.cma_result <- function(cma_result) {
  cma_result$best.fitness
}


get_best_fitness.dream <- function(dream) {
  # copied from getS3method("coef", class = "dream"). How coef.dream get the best parameter set
  sss <- window(dream)
  ppp <- window(as.mcmc(dream$hist.logp), start = start(sss), thin = thin(sss))
  return(-1 * max(ppp))
}


## get the exit message from the optimiser =====================================
get_exit_message <- function(x) {
  UseMethod("get_exit_message")
}

get_exit_message.cma_result <- function(cma_result) {
  cma_result$message
}

get_exit_message.dream <- function(dream) {
  dream$EXITMSG
}


## dream only - get the range of values tested =================================
# future does not like local generics. If I want to make this work 
# I need to turn it into a package
# Solution copy each function

transform_to_realspace <- function(x, ...) { # generic ... allows extra arguments
  UseMethod("transform_to_realspace")
}


single_col_mcmc_transform_to_realspace <- function(single_col_mcmc_result, transform_function, lower_bound, upper_bound, scale) {
  
  transform_function <- noquote(transform_function)
  
  transform_function(
    parameter_set = single_col_mcmc_result,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    scale = 100
  )
  
}


# DELETE when package is made
transform_to_realspace_mcmc <- function(mcmc_result, transform_functions, lower_bounds, upper_bounds, scale) {
  
  # Convert mcmc_result into list for pmap
  mcmc_as_list <- mcmc_result |> 
    unclass() |> 
    as_tibble() |> 
    as.list()
  
  
  # Map over each column
  real_space_as_list <- pmap(
    .l = list(mcmc_as_list, transform_functions, lower_bounds, upper_bounds),
    .f = single_col_mcmc_transform_to_realspace,
    scale = scale
  ) 
  
  # convert back to mcmc
  real_space_as_matrix <- do.call(cbind, real_space_as_list)
  colnames(real_space_as_matrix) <- paste0("par", seq_len(ncol(real_space_as_matrix)))
  
  return(as.mcmc(real_space_as_matrix))
  
}




transform_to_realspace.mcmc <- function(mcmc_result, transform_functions, lower_bounds, upper_bounds, scale) {
  
  # Convert mcmc_result into list for pmap
  mcmc_as_list <- mcmc_result |> 
    unclass() |> 
    as_tibble() |> 
    as.list()
  
  
  # Map over each column
  real_space_as_list <- pmap(
    .l = list(mcmc_as_list, transform_functions, lower_bounds, upper_bounds),
    .f = single_col_mcmc_transform_to_realspace,
    scale = scale
  ) 
  
  # convert back to mcmc
  real_space_as_matrix <- do.call(cbind, real_space_as_list)
  colnames(real_space_as_matrix) <- paste0("par", seq_len(ncol(real_space_as_matrix)))
  
  return(as.mcmc(real_space_as_matrix))
  
}

# DELETE with package
transform_to_realspace_mcmc_list <- function(mcmc_result, transform_functions, lower_bounds, upper_bounds, scale) {
  
  # repeat based on the length of mcmc.list
  real_space_as_list <- map(
    .x = mcmc_result,
    .f = transform_to_realspace_mcmc,
    transform_functions = transform_functions,
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    scale = scale
  )
  
  return(as.mcmc.list(real_space_as_list))
  
}

transform_to_realspace.mcmc.list <- function(mcmc_result, transform_functions, lower_bounds, upper_bounds, scale) {
  
  # repeat based on the length of mcmc.list
  real_space_as_list <- map(
    .x = mcmc_result,
    .f = transform_to_realspace,
    transform_functions = transform_functions,
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    scale = scale
  )
  
  return(as.mcmc.list(real_space_as_list))
  
}


transform_to_realspace_dream <- function(mcmc_result, transform_functions, lower_bounds, upper_bounds, scale) {
  
  Sequences <- transform_to_realspace_mcmc_list(
    mcmc_result = mcmc_result[["Sequences"]],
    transform_functions = transform_functions,
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    scale = scale
  )
  
  mcmc_result$Sequences <- Sequences
  return(mcmc_result)
}


transform_to_realspace.dream <- function(mcmc_result, transform_functions, lower_bounds, upper_bounds, scale) {
  
  Sequences <- transform_to_realspace(
    mcmc_result = mcmc_result[["Sequences"]],
    transform_functions = transform_functions,
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    scale = scale
  )
  
  mcmc_result$Sequences <- Sequences
  return(mcmc_result)
}







#get_sequences <- function(x, ...) { # generic ... allows extra arguments
#  UseMethod("get_sequences")
#}


#get_sequences.dream <- function(dream) {
  # This is an mcmc object
 # transform_to_realspace(
  #  mcmc_result = dream$Sequences,
   # transform_functions = dream$numerical_optimiser_setup$transform_parameter_methods, 
    #lower_bounds = dream$numerical_optimiser_setup$lower_bound, 
    #upper_bounds = dream$numerical_optimiser_setup$upper_bound, 
    #scale = dream$numerical_optimiser_setup$scale
  #)

#}



#get_sequences.default <- function(x, ...) {
 # NULL
#}


## AIC =========================================================================
get_AIC <- function(x) {
  (2 * length(coef(x))) + (2 * get_best_fitness(x))
}

## Get number of restarts - CMAES only =========================================
get_restart_count <- function(cmaes_result, ...) {
  UseMethod("get_restart_count")
}


get_restart_count.cma_result <- function(cmaes_result, ...) {
  cmaes_result$n.restarts
}

get_restart_count.default <- function(cmaes_result, ...) {
  NULL
}

## Covert parameter into real space ============================================
## uses coef
get_best_parameters_real_space <- function(cmaes_or_dream_result) {
  
  force(cmaes_or_dream_result)
  
  scaled_parameters <- coef(cmaes_or_dream_result)
  
  best_parameters <- purrr::pmap( 
    .l = list(
      cmaes_or_dream_result$numerical_optimiser_setup$transform_parameter_methods, 
      seq(from = 1, to = length(cmaes_or_dream_result$numerical_optimiser_setup$parameter_names)), 
      cmaes_or_dream_result$numerical_optimiser_setup$lower_bound, 
      cmaes_or_dream_result$numerical_optimiser_setup$upper_bound
    ),
    .f = transform_parameter_method,
    parameter_set = as.matrix(scaled_parameters, ncol = 1),
    scale = cmaes_or_dream_result$numerical_optimiser_setup$scale
  )
  
  return(unlist(best_parameters))
  
}



get_boxcox_streamflow <- function(cmaes_or_dream_result) {
  best_parameters <- get_best_parameters_real_space(cmaes_or_dream_result)
  
  cmaes_or_dream_result$numerical_optimiser_setup$streamflow_model(
    catchment_data = cmaes_or_dream_result$numerical_optimiser_setup$catchment_data$full_data_set,
    parameter_set = as.matrix(best_parameters, ncol = 1)
  )
}



plot.catchment_data <- function(x) {
  
  plot(
    x = x$full_data_set$year,
    y = x$full_data_set$observed_boxcox_streamflow,
    type = "b",
    xlab = "Year",
    ylab = "Box-Cox Streamflow"
  )
  
}

plot.result_set <- function(x) {
  
  # should only plot what I have optimised for
  # remove data not in calibration
  # easiest way is to use the full observed streamflow and an index i.e., if
  # there is a NA it was not calibrated on
  keep_data <- !is.na(x$numerical_optimiser_setup$catchment_data$full_data_set$observed_boxcox_streamflow)
  
  plot(
    x = x$numerical_optimiser_setup$catchment_data$full_data_set$year[keep_data],
    y = x$optimised_boxcox_streamflow[keep_data], 
    type = "b", 
    xlab = "Year",
    ylab = "Box-Cox Streamflow", 
    ylim = range(c(x$optimised_boxcox_streamflow[keep_data], x$numerical_optimiser_setup$catchment_data$full_data_set$observed_boxcox_streamflow[keep_data]), na.rm = TRUE))
  lines(
    x = x$numerical_optimiser_setup$catchment_data$full_data_set$year[keep_data],
    y = x$numerical_optimiser_setup$catchment_data$full_data_set$observed_boxcox_streamflow[keep_data],
    type = "b",
    col = "red"
  )
  legend(
    x = "topleft", 
    c("modelled", "observed"), 
    lty = c(1, 1), 
    col = c("black", "red"))
  
}


summary.result_set <- function(x) {
  
  #summary() # print a summary
  #should show the the exit_message, best AIC and fitted parameters. List the models used and gauge
  cat("Best fitnesss (AIC):", x$AIC_best_parameter_set)
  cat("\nBest Parameters:", paste(names(dream_example$best_parameter_set), signif(dream_example$best_parameter_set, 3), sep = ":", collapse = " "))
  
}




# is empty tibble
is_empty_tibble <- function(x) {
  if_else(nrow(x) * ncol(x) == 0, TRUE, FALSE)
}
