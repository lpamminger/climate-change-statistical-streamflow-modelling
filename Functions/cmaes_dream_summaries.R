check_near_bounds <- function(result_set_object) {
  
  # Get information ------------------------------------------------------------
  lower_bound <- result_set_object$numerical_optimiser_setup$lower_bound
  upper_bound <- result_set_object$numerical_optimiser_setup$upper_bound
  bound_range <- upper_bound - lower_bound
  calibrated_parameters <- result_set_object$best_parameter_set
  
  # Use order of magnitude of bound range to calculate tolerance to bounds -----
  order_of_magnitude_bound_range <- get_order_magnitude(bound_range)
  
  # Take 3 orders of magnitude off for near()
  remove_orders_of_magnitude <- 4
  new_exponent <- log10(order_of_magnitude_bound_range) - remove_orders_of_magnitude
  near_tolerance_per_parameter <- 10^new_exponent
  
  
  # Check lower_bound ----------------------------------------------------------
  near_lower_bound <- near(calibrated_parameters, lower_bound, tol = near_tolerance_per_parameter)
  
  # Check upper_bound ----------------------------------------------------------
  near_upper_bound <- near(calibrated_parameters, upper_bound, tol = near_tolerance_per_parameter)
  
  # Return string
  output_message <- case_when(
    # parameters can only be near either upper or lower bound
    near_lower_bound ~ "near_lower_bound",
    near_upper_bound ~ "near_upper_bound",
    .default = NA
  )
  
  return(output_message)
}





parameters_summary <- function(x) {
  tibble::as_tibble(
    list(
      "gauge" = x$numerical_optimiser_setup$catchment_data$gauge_ID,   
      "streamflow_model" = x$numerical_optimiser_setup$streamflow_model()[[1]],
      "objective_function" = x$numerical_optimiser_setup$objective_function()[[1]],
      "optimiser" = s3_class(x$optimised_object)[1], 
      "parameter" = x$numerical_optimiser_setup$parameter_names, 
      "parameter_value" = x$best_parameter_set, 
      "loglikelihood" = x$LL_best_parameter_set,
      "AIC" = x$AIC_best_parameter_set,
      "exit_message" = x$exit_message,
      "near_bounds" = check_near_bounds(x)
    )
  )
}



sequences_summary <- function(x) { # NOT IN USE
  
  # This only take a result_object 
  
  x$sequences |> 
    tidyr::pivot_longer(
      cols = everything(),
      names_to = "parameter",
      values_to = "parameter_values" 
    ) |> 
    tibble::add_column(
      gauge = x$numerical_optimiser_setup$catchment_data$gauge_ID,
      .before = 1
    ) |> 
    tibble::add_column(
      streamflow_model = x$numerical_optimiser_setup$streamflow_model()$name,
      .before = 2
    ) |> 
    tibble::add_column(
      objective_function = x$numerical_optimiser_setup$objective_function()$name,
      .before = 3
    )
  
}







streamflow_timeseries_summary <- function(x) {
  # The modelled_boxcox_streamflow is just the streamflow optimised
  # We must use data being optimised instead of the full_data_set
  
  # Only works with result_set objects - force?
  calibrated_data <- x$numerical_optimiser_setup$catchment_data$stop_start_data_set |> 
    list_rbind()
  

  
  tibble::as_tibble(
    list(
      "gauge" = x$numerical_optimiser_setup$catchment_data$gauge_ID,
      "year" = calibrated_data |> pull(year),
      "precipitation" = calibrated_data |> pull(precipitation),
      "realspace_observed_streamflow" = calibrated_data |> pull(observed_streamflow),
      "realspace_modelled_streamflow" = x$optimised_modelled_streamflow_realspace |> as.double(), # as.double converts from matrix column to vector
      "transformed_observed_streamflow" = x$transformed_observed_streamflow |> as.double(), # as.double converts from matrix column to vector
      "transformed_modelled_streamflow" = x$optimised_modelled_streamflow_transformed_space,
      "streamflow_model" = x$numerical_optimiser_setup$streamflow_model()$name,
      "objective_function" = x$numerical_optimiser_setup$objective_function()$name,
      "streamflow_transform_method" = x$numerical_optimiser_setup$streamflow_transform_method()$name
    )
  ) 
}






# Allows running in parallel with chunking to not exceed RAM
run_and_save_chunks_optimiser_parallel <- function(chunked_numerical_optimisers, chunk_id, optimiser, save_streamflow, save_sequences, is_drought) {
  
  tictoc::tic()
  
  calibrated_results <- furrr::future_map(
    .x = chunked_numerical_optimisers,
    .f = optimiser,
    print_monitor = FALSE,
    .options = furrr_options(
      globals = TRUE,
      seed = TRUE
    )
  )
  
  
  # Purrr does not work in parallel, so I don't need plan(sequential)
  sort_results <- purrr::map(.x = calibrated_results, .f = result_set)
  
  optimiser_name <- as.character(substitute(optimiser))
  
  # I do not want to assign a variable name. Garbage collector works better like this.
  purrr::map(.x = sort_results, .f = parameters_summary) |>
    purrr::list_rbind() |>
    readr::write_csv(
      file = paste0(
        "./Results/",
        optimiser_name,
        "/",
        if_else(is_drought, "drought_", ""),
        "parameter_results_chunk_",
        chunk_id,
        "_",
        get_date(),
        ".csv"
      )
    )
  
  if (save_streamflow) {
    purrr::map(
      .x = sort_results,
      .f = modelled_streamflow_summary
    ) |>
      purrr::list_rbind() |>
      readr::write_csv(
        file = paste0(
          "./Results/",
          optimiser_name,
          "/",
          if_else(is_drought, "drought_", ""),
          "streamflow_results_chunk_",
          chunk_id,
          "_",
          get_date(),
          ".csv"
        )
      )
  }
  
  if (save_sequences) {
    purrr::map(
      .x = sort_results,
      .f = sequences_summary
    ) |>
      purrr::list_rbind() |>
      readr::write_csv(
        file = paste0(
          "./Results/",
          optimiser_name,
          "/",
          if_else(is_drought, "drought_", ""),
          "sequences_results_chunk_",
          chunk_id,
          "_",
          get_date(),
          ".csv"
        )
      )
  }
  
  
  
  cat(paste0("Chunk ", chunk_id, " complete "))
  tictoc::toc()
  cat("\n")
  
  
  # Remove objects for garbage collection
  rm(list = c("calibrated_results", "sort_results", "optimiser_name", "save_streamflow", "save_sequences", "is_drought"))
  
  # Call garbage collection
  gc()
}




