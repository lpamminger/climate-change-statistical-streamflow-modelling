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
  # This should only show the streamflow used in calibration
  best_parameters <- get_best_parameters_real_space(cmaes_or_dream_result)

  # Apply to start_stop_data_set (list of tibbles)
  start_stop_data_set <- cmaes_or_dream_result$numerical_optimiser_setup$catchment_data$stop_start_data_set |>
    list_rbind()

  observed_streamflow <- start_stop_data_set$observed_streamflow

  # Transforms change depending on objective function used
  if (cmaes_or_dream_result$numerical_optimiser_setup$objective_function()$name == "constant_sd_boxcox_objective_function") {
    best_lambda <- best_parameters[length(best_parameters)]

    transformed_observed_streamflow <- boxcox_transform(
      y = observed_streamflow,
      lambda = best_lambda,
      lambda_2 = 1
    )
  } else if (cmaes_or_dream_result$numerical_optimiser_setup$objective_function()$name == "constant_sd_log_sinh_objective_function") {
    best_a <- best_parameters[length(best_parameters) - 1]

    best_b <- best_parameters[length(best_parameters)]

    transformed_observed_streamflow <- log_sinh_transform(
      a = best_a,
      b = best_b,
      y = observed_streamflow
    )
  } else {
    stop("Name of objective function not found")
  }

  return(transformed_observed_streamflow)
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
  transformed_modelled_streamflow <- get_transformed_optimised_streamflow(cmaes_or_dream_result)

  best_parameters <- get_best_parameters_real_space(cmaes_or_dream_result)

  # this will change based on the transformation in the objective_function
  if (cmaes_or_dream_result$numerical_optimiser_setup$objective_function()$name == "constant_sd_boxcox_objective_function") {
    best_lambda <- best_parameters[length(best_parameters)]

    realspace_modelled_streamflow <- boxcox_inverse_transform(
      yt = transformed_modelled_streamflow,
      lambda = best_lambda,
      lambda_2 = 1
    )
  } else if (cmaes_or_dream_result$numerical_optimiser_setup$objective_function()$name == "constant_sd_log_sinh_objective_function") {
    best_a <- best_parameters[length(best_parameters) - 1]

    best_b <- best_parameters[length(best_parameters)]

    realspace_modelled_streamflow <- inverse_log_sinh_transform(
      a = best_a,
      b = best_b,
      z = transformed_modelled_streamflow,
      offset = 300
    )
  } else {
    stop("Name of objective function not found")
  }
  
  # realspace_modelled_streamflow cannot be less than zero
  # if less than zero set to zero
  realspace_modelled_streamflow[realspace_modelled_streamflow < 0] <- 0


  return(realspace_modelled_streamflow)
}







summary.result_set <- function(x) {
  # summary() # print a summary
  # should show the the exit_message, best AIC and fitted parameters. List the models used and gauge
  cat("Best fitnesss (AIC):", x$AIC_best_parameter_set)
  cat("\nBest Parameters:", paste(names(dream_example$best_parameter_set), signif(dream_example$best_parameter_set, 3), sep = ":", collapse = " "))
}




# is empty tibble
is_empty_tibble <- function(x) {
  if_else(nrow(x) * ncol(x) == 0, TRUE, FALSE)
}








plot_result_set_v2 <- function(x, type) {
  stopifnot(type %in% c("streamflow-time", "rainfall-runoff"))
  # This should only show the streamflow used in calibration

  # Get precipitation and observed streamfow from stop_start_index
  observed_data <- x$numerical_optimiser_setup$catchment_data$stop_start_data_set |>
    list_rbind()
  
  # Identify streamflow transformation method in objective function
  streamflow_transformation_method <- x$numerical_optimiser_setup$objective_function()[1] |> 
    str_remove("constant_sd_") |> 
    str_remove("_objective_function")
  



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
      theme_bw() +
      theme(
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9)
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
        legend.position.inside = c(0.9, 0.9)
      )
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
