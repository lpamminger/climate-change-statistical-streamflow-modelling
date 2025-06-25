objective_function_setup <- function(streamflow_model, parameter_transform_method, streamflow_transform_method, lower_bound, upper_bound, objective_function, catchment_data, scale, minimise_likelihood, streamflow_transform_method_offset) {
  force(streamflow_model)
  force(parameter_transform_method)
  force(streamflow_transform_method)
  force(lower_bound)
  force(upper_bound)
  force(catchment_data)
  force(objective_function)
  force(scale)
  force(minimise_likelihood)
  force(streamflow_transform_method_offset)

  # Unpack everything
  streamflow_model <- noquote(streamflow_model)
  objective_function <- noquote(objective_function)
  streamflow_transform_method <- noquote(streamflow_transform_method)
  stop_start_data_set <- catchment_data$stop_start_data_set # I am only using the stop start (BAD)
  
  # apply_truncnorm
  if(streamflow_transform_method()$name == "boxcox_transform") {
    apply_truncnorm <-  TRUE
  } else {
    apply_truncnorm <-  FALSE
  }
  



  function(parameter_set) {

    # check if parameter is a matrix
    if (!is.matrix(parameter_set)) {
      parameter_set <- as.matrix(parameter_set, ncol = 1)
    }

    stopifnot(is.matrix(parameter_set))


    # 1. Transform parameters into actual values (then matrix it) --------------
    real_parameter_set <- purrr::pmap(
      .l = list(parameter_transform_method, seq(from = 1, to = nrow(parameter_set)), lower_bound, upper_bound),
      .f = transform_parameter_method,
      parameter_set = parameter_set,
      scale = scale
    )

    real_parameter_set <- matrix(
      unlist(real_parameter_set), # cannot use unlist by itself. unlist coerces data into atomic vector
      nrow = nrow(parameter_set),
      byrow = TRUE
    )


    # 2. Put values into streamflow model (then matrix it) --------------------
    transformed_modelled_streamflow <- purrr::map(
      .x = stop_start_data_set,
      .f = streamflow_model,
      parameter_set = real_parameter_set
    )


    # streamflow model already produces modelled streamflow in transformed
    # space

    # 3. Transform observed_streamflow from realspace into transformed space ---
    observed_streamflow <- purrr::map(
      .x = stop_start_data_set,
      .f = pull_vector_from_split_data_set,
      column_name = observed_streamflow
    )


    observed_streamflow <- purrr::map2(
      .x = observed_streamflow,
      .y = lengths(observed_streamflow),
      .f = matrix,
      ncol = ncol(parameter_set)
    )


    # repeat - for discontinuous gauges
    repeat_real_parameter_set <- rep(list(real_parameter_set), times = length(observed_streamflow))
    
    
    transformed_observed_streamflow <- map2(
      .x = observed_streamflow,
      .y = repeat_real_parameter_set,
      .f = select_streamflow_transform_method,
      streamflow_transform_method = streamflow_transform_method,
      offset = streamflow_transform_method_offset
    )


    

    # 4. Put values into objective_function ------------------------------------
    # needs to be in a list for discontinous data
    LL_scores <- pmap(
      .l = list(
        transformed_modelled_streamflow,
        transformed_observed_streamflow
      ),
      .f = objective_function,
      parameter_set = real_parameter_set,
      apply_truncnorm = apply_truncnorm
    )


    LL_scores <- do.call("rbind", LL_scores) |> colSums()


    if (!minimise_likelihood) {
      return(-LL_scores)
    }

    return(LL_scores)
  }
}
