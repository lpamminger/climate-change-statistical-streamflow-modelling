# Generic functions ------------------------------------------------------------

## override coef.dream because it doesn't work well ============================
# the default coef.dream does not work when the mcmc.list has be transformed
# Instead, I transform in DREAM function itself

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







# get_sequences <- function(x, ...) { # generic ... allows extra arguments
#  UseMethod("get_sequences")
# }


# get_sequences.dream <- function(dream) {
# This is an mcmc object
# transform_to_realspace(
#  mcmc_result = dream$Sequences,
# transform_functions = dream$numerical_optimiser_setup$transform_parameter_methods,
# lower_bounds = dream$numerical_optimiser_setup$lower_bound,
# upper_bounds = dream$numerical_optimiser_setup$upper_bound,
# scale = dream$numerical_optimiser_setup$scale
# )

# }



# get_sequences.default <- function(x, ...) {
# NULL
# }


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
      cmaes_or_dream_result$numerical_optimiser_setup$parameter_transform_method,
      seq(from = 1, to = length(cmaes_or_dream_result$numerical_optimiser_setup$parameter_names)),
      cmaes_or_dream_result$numerical_optimiser_setup$lower_bound,
      cmaes_or_dream_result$numerical_optimiser_setup$upper_bound
    ),
    .f = transform_parameter_method,
    parameter_set = as.matrix(scaled_parameters, ncol = 1),
    scale = cmaes_or_dream_result$numerical_optimiser_setup$scale
  ) |> 
    unlist()
  
  # If streamflow_transform_method is boxcox and lambda is less than machine tol
  # then set lambda to zero. This is what boxcox_transform() does. Reflect this
  # in result
  if(cmaes_or_dream_result$numerical_optimiser_setup$streamflow_transform_method()$name == "boxcox_transform") {
    near_zero_lambda <- best_parameters[["lambda"]] <= .Machine$double.eps^0.5
    
    if(near_zero_lambda) {
      best_parameters[["lambda"]] <- 0
    }
  }
  

  return(best_parameters)
}




get_boxcox_streamflow <- function(cmaes_or_dream_result) {
  # This should only show the streamflow used in calibration

  best_parameters <- get_best_parameters_real_space(cmaes_or_dream_result)

  # Apply to start_stop_data_set (list of tibbles)
  start_stop_data_set <- cmaes_or_dream_result$numerical_optimiser_setup$catchment_data$stop_start_data_set

  # Streamflow model
  streamflow_model <- cmaes_or_dream_result$numerical_optimiser_setup$streamflow_model

  # Make streamflow using best parameters and start_stop_data_set
  streamflow_results <- map(
    .x = start_stop_data_set,
    .f = streamflow_model,
    parameter_set = best_parameters
  ) |>
    unname() |>
    unlist()

  return(streamflow_results)
}


get_transformed_observed_streamflow <- function(cmaes_or_dream_result) {
  
  # get realspace streamflow
  realspace_streamflow <- cmaes_or_dream_result$numerical_optimiser_setup$catchment_data$stop_start_data_set |>
    list_rbind() |> 
    pull(observed_streamflow)
  
  # get transform function
  # get parameters for transform
  # apply transform
  
  select_streamflow_transform_method(
    streamflow_transform_method = cmaes_or_dream_result$numerical_optimiser_setup$streamflow_transform_method,
    parameter_set = as.matrix(get_best_parameters_real_space(cmaes_or_dream_result), ncol = 1), # function relies on matrices as inputs
    timeseries = as.matrix(realspace_streamflow, ncol = 1),
    offset = cmaes_or_dream_result$numerical_optimiser_setup$streamflow_transform_method_offset
  )
  
}

get_transformed_optimised_streamflow <- function(cmaes_or_dream_result) {
  
  # This should only show the streamflow used in calibration
  best_parameters <- get_best_parameters_real_space(cmaes_or_dream_result)

  # Apply to start_stop_data_set (list of tibbles)
  start_stop_data_set <- cmaes_or_dream_result$numerical_optimiser_setup$catchment_data$stop_start_data_set

  # Streamflow model
  streamflow_model <- cmaes_or_dream_result$numerical_optimiser_setup$streamflow_model

  # Make streamflow using best parameters and start_stop_data_set
  streamflow_results <- map(
    .x = start_stop_data_set,
    .f = streamflow_model,
    parameter_set = best_parameters
  ) |>
    unname() |>
    unlist()

  return(streamflow_results)
}

get_realspace_optimised_streamflow <- function(cmaes_or_dream_result) {
  
  # Requires and inverse transform on modelled streamflow
  transformed_modelled_streamflow <- get_transformed_optimised_streamflow(cmaes_or_dream_result)
  

  # Get inverse transform - add inverse to front of function
  # this requires all streamflow_transform methods to have inverse_function_name style
  inverse_streamflow_transform_method_name <- paste0("inverse_", cmaes_or_dream_result$numerical_optimiser_setup$streamflow_transform_method()$name)
  
  realspace_modelled_streamflow <- select_streamflow_transform_method(
    streamflow_transform_method = match.fun(FUN = inverse_streamflow_transform_method_name),
    parameter_set = as.matrix(get_best_parameters_real_space(cmaes_or_dream_result), ncol = 1), # function relies on matrices as inputs
    timeseries = as.matrix(transformed_modelled_streamflow, ncol = 1),
    offset = cmaes_or_dream_result$numerical_optimiser_setup$streamflow_transform_method_offset
  )
  
  
  # realspace_modelled_streamflow cannot be less than zero
  # if less than zero set to zero
  # modify the results
  
  realspace_modelled_streamflow[realspace_modelled_streamflow < 0] <- 0


  return(realspace_modelled_streamflow)
}





# is empty tibble
is_empty_tibble <- function(x) {
  if_else(nrow(x) * ncol(x) == 0, TRUE, FALSE)
}





plot.result_set <- function(x, type) {
  
  stopifnot(type %in% c("streamflow-time", "rainfall-runoff", "examine_transform"))
  # This should only show the streamflow used in calibration

  # Get precipitation and observed streamfow from stop_start_index
  observed_data <- x$numerical_optimiser_setup$catchment_data$stop_start_data_set |>
    list_rbind()
  
  # Identify streamflow transformation method in objective function
  streamflow_transformation_method <- x$numerical_optimiser_setup$streamflow_transform_method()[[1]]
  



  # Plotting
  if (type == "streamflow-time") {
  
    
    # Create tibble for plotting
    streamflow_results <- list(
      year = observed_data |> pull(year),
      precipitation = observed_data |> pull(precipitation),
      observed_streamflow = observed_data |> pull(observed_streamflow),
      modelled_streamflow = x$optimised_modelled_streamflow_realspace
    ) |>
      as_tibble() |>
      pivot_longer(
        cols = contains("streamflow"),
        names_to = "observed_or_modelled",
        values_to = "streamflow"
      ) |>
      mutate(
        observed_or_modelled = if_else(observed_or_modelled == "modelled_streamflow", "Modelled Streamflow", "Observed Streamflow")
      )




    streamflow_results |>
      ggplot(aes(x = year, y = streamflow, colour = observed_or_modelled)) +
      geom_line() +
      geom_point() +
      labs(
        x = "Year",
        y = "Streamflow",
        colour = NULL
      ) +
      scale_colour_brewer(palette = "Set1") +
      labs(
        x = "Year",
        y = "Streamflow (mm)"
      ) +
      theme_bw() +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9),
        legend.background = element_rect(colour = "black")
      )
    
    
  } else if (type == "rainfall-runoff") {
    # Create tibble for plotting
    streamflow_results <- list(
      year = observed_data |> pull(year),
      precipitation = observed_data |> pull(precipitation),
      observed_streamflow = x$transformed_observed_streamflow,
      modelled_streamflow = x$optimised_modelled_streamflow_transformed_space
    ) |>
      as_tibble() |>
      pivot_longer(
        cols = contains("streamflow"),
        names_to = "observed_or_modelled",
        values_to = "streamflow"
      ) |>
      mutate(
        observed_or_modelled = if_else(observed_or_modelled == "modelled_streamflow", "Modelled Streamflow", "Observed Streamflow")
      )


    streamflow_results |>
      ggplot(aes(x = precipitation, y = streamflow, colour = observed_or_modelled)) +
      geom_smooth(
        formula = y ~ x,
        method = lm,
        se = FALSE
      ) +
      geom_point() +
      labs(
        x = "Precipitation (mm)",
        y = paste0("Streamflow (", streamflow_transformation_method, ")"),
        colour = NULL
      ) +
      scale_colour_brewer(palette = "Set1") +
      theme_bw() +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.1, 0.9),
        legend.background = element_rect(colour = "black")
      )
    
    
    
  } else if (type == "examine_transform") {
    # observed streamflow - x-axis and transformed streamflow on y-axis
    
    modelled_streamflow_data <- list(
      realspace_modelled_streamflow = x$optimised_modelled_streamflow_realspace,#observed_data |> pull(observed_streamflow),
      transformed_modelled_streamflow = x$optimised_modelled_streamflow_transformed_space#x$transformed_observed_streamflow
    ) |> 
      as_tibble() |> 
      arrange(realspace_modelled_streamflow)
      
    
    # This is for the plot only using model parameters
    range_of_values <- modelled_streamflow_data |> 
      pull(realspace_modelled_streamflow) |> 
      range()
    
    realspace_modelled_axis <- seq(from = range_of_values[1], to = range_of_values[2], by = 0.01)
    
    # Plot the curve using calibrated values
    transform_name <- x$numerical_optimiser_setup$streamflow_transform_method()$name
    
    if(transform_name == "log_sinh_transform") {
      a <- x$best_parameter_set[length(x$best_parameter_set) - 2]
      b <- x$best_parameter_set[length(x$best_parameter_set) - 1]
      
      transformed_modelled_axis <- x$numerical_optimiser_setup$streamflow_transform_method(a, b, realspace_modelled_axis)

    } else if(transform_name == "boxcox_transform") {
      lambda <- x$best_parameter_set[length(x$best_parameter_set) - 1]
      
      lambda_2 <- standardised_results$numerical_optimiser_setup$streamflow_transform_method_offset
      
      transformed_modelled_axis <- x$numerical_optimiser_setup$streamflow_transform_method(realspace_modelled_axis, lambda, lambda_2)
      
    } else {
      stop("Unrecognised transform")
    }
    
    # Get curve data into tibble
    curve_data <- list(
      "realspace_modelled_streamflow" = realspace_modelled_axis,
      "transformed_modelled_streamflow" = transformed_modelled_axis
    ) |> 
      as_tibble()
    
    
    
    
    modelled_streamflow_data |> 
      ggplot(aes(x = realspace_modelled_streamflow, y = transformed_modelled_streamflow)) +
      geom_line(data = curve_data, colour = "red") +
      geom_point() +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme_bw()
    
  }
}





plot.catchment_data <- function(x, type) {
  # This should be the entire dataset

  stopifnot(type %in% c("streamflow-time", "rainfall-runoff"))

  if (any(is.na(x$full_data_set$observed_streamflow))) {
    message("NA's in observed_streamflow. Removed from plot.")
  }


  if (type == "streamflow-time") {
    x$full_data_set |>
      ggplot(aes(x = year, y = observed_streamflow)) +
      geom_line(na.rm = TRUE) +
      geom_point(na.rm = TRUE) +
      labs(
        x = "Year",
        y = "Observed Streamflow (mm)"
      ) +
      theme_bw()
  } else if (type == "rainfall-runoff") {
    x$full_data_set |>
      ggplot(aes(x = precipitation, y = observed_streamflow)) +
      labs(
        x = "Precipitation (mm)",
        y = "Observed Streamflow (mm)"
      ) +
      geom_point(na.rm = TRUE) +
      theme_bw()
  }
}
