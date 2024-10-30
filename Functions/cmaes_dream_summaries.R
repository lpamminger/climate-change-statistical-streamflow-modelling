parameters_summary <- function(x) {
  tibble::as_tibble(
    list(
      "gauge" = x$numerical_optimiser_setup$catchment_data$gauge_ID,   
      "streamflow_model" = x$numerical_optimiser_setup$streamflow_model()[[1]],
      "objective_function" = x$numerical_optimiser_setup$objective_function()[[1]],
      "optimiser" = s3_class(x$optimised_object)[1], # [1] needed to get child class
      "parameter" = x$numerical_optimiser_setup$parameter_names, 
      "parameter_value" = x$best_parameter_set, 
      "loglikelihood" = x$LL_best_parameter_set,
      "AIC" = x$AIC_best_parameter_set,
      "exit_message" = x$exit_message,
      "near_bounds" = (near(x$numerical_optimiser_setup$lower_bound, x$best_parameter_set)) | (near(x$numerical_optimiser_setup$upper_bound, x$best_parameter_set))
    )
  )
}


restarts_summary <- function(x) {
  tibble::as_tibble(
    list(
      "gauge" = x$numerical_optimiser_setup$catchment_data$gauge_ID,   
      "streamflow_model" = x$numerical_optimiser_setup$streamflow_model()[[1]],
      "objective_function" = x$numerical_optimiser_setup$objective_function()[[1]],
      "optimiser" = sloop::s3_class(x$optimised_object)[1], # [1] needed to get child class
      "loglikelihood" = x$LL_best_parameter_set,
      "AIC" = x$AIC_best_parameter_set,
      "restarts" = x$restart_count,
      "exit_message" = x$exit_message
    )
  )
}



sequences_summary <- function(x) {
  
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

modelled_streamflow_summary <- function(x) {
  tibble::as_tibble(
    list(
      "year" = x$numerical_optimiser_setup$catchment_data$full_data_set$year,
      "precipitation" = x$numerical_optimiser_setup$catchment_data$full_data_set$precipitation,
      "observed_boxcox_streamflow" = x$numerical_optimiser_setup$catchment_data$full_data_set$observed_boxcox_streamflow,
      "modelled_boxcox_streamflow" = c(x$optimised_boxcox_streamflow),
      "gauge" = x$numerical_optimiser_setup$catchment_data$gauge_ID,
      "streamflow_model" = x$numerical_optimiser_setup$streamflow_model()$name,
      "objective_function" = x$numerical_optimiser_setup$objective_function()$name,
      "optimiser" = sloop::s3_class(x$optimised_object)[1],
      "loglikelihood" = x$LL_best_parameter_set
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




