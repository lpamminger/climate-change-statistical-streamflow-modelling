objective_function_setup <- function(streamflow_model, transform_parameter_methods, lower_bound, upper_bound, objective_function, catchment_data, scale, minimise_likelihood) {
  force(catchment_data)
  force(streamflow_model)
  force(objective_function)
  force(transform_parameter_methods)
  force(lower_bound)
  force(upper_bound)
  force(scale)

  # Unpack everything
  streamflow_model <- noquote(streamflow_model)
  objective_function <- noquote(objective_function)
  stop_start_data_set <- catchment_data$stop_start_data_set # I am only using the stop start (BAD)



  function(parameter_set) {
    # check if parameter is a matrix
    if (!is.matrix(parameter_set)) {
      parameter_set <- as.matrix(parameter_set, ncol = 1)
    }

    stopifnot(is.matrix(parameter_set))

    # 1. Transform parameters into actual values (then matrix it)
    # 2. Put values into streamflow model (then matrix it)
    # 3. Put values into objective_function
    # 4. Repeat 2. and 3. for each row in start_stop_index using map2_dbl
    # 5. Make optimiser_parameter_set
    # 6. sum LL from optimiser

    real_parameter_set <- purrr::pmap(
      .l = list(transform_parameter_methods, seq(from = 1, to = nrow(parameter_set)), lower_bound, upper_bound),
      .f = transform_parameter_method,
      parameter_set = parameter_set,
      scale = scale
    )

    real_parameter_set <- matrix(unlist(real_parameter_set),
      nrow = nrow(parameter_set),
      byrow = TRUE
    )


    modelled_streamflow <- purrr::map(
      .x = stop_start_data_set,
      .f = streamflow_model,
      parameter_set = real_parameter_set
    )


    observed_streamflow <- purrr::map(
      .x = stop_start_data_set,
      .f = pull_vector_from_split_data_set,
      column_name = observed_boxcox_streamflow
    )


    observed_streamflow <- purrr::map2(
      .x = observed_streamflow,
      .y = lengths(observed_streamflow), # nrow
      .f = matrix,
      ncol = ncol(parameter_set)
    )
    
    LL_scores <- purrr::pmap(
      .l = list(
        modelled_streamflow,
        observed_streamflow,
        stop_start_data_set
      ),
      .f = objective_function,
      parameter_set = real_parameter_set
    )


    LL_scores <- colSums(matrix(unlist(LL_scores), ncol = ncol(real_parameter_set), byrow = TRUE))

    if (!minimise_likelihood) {
      return(-LL_scores)
    }

    return(LL_scores)
  }
}
