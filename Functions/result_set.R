# result_set constructor -------------------------------------------------------

result_set <- function(cmaes_or_dream_result) {

  stopifnot(s3_class(cmaes_or_dream_result)[1] %in% c("dream", "cma_result"))
  stopifnot(s3_class(cmaes_or_dream_result$numerical_optimiser_setup)[1] == "numerical_optimiser_setup")
  
  
  ## Make class ================================================================
  structure(
    list(
      "optimised_object" = cmaes_or_dream_result,
      "numerical_optimiser_setup" = cmaes_or_dream_result$numerical_optimiser_setup,
      "best_parameter_set" = get_best_parameters_real_space(cmaes_or_dream_result),
      "LL_best_parameter_set" = get_best_fitness(cmaes_or_dream_result),
      "AIC_best_parameter_set" = get_AIC(cmaes_or_dream_result),
      "exit_message" = get_exit_message(cmaes_or_dream_result),
      "transformed_observed_streamflow" = get_transformed_observed_streamflow(cmaes_or_dream_result),
      "optimised_modelled_streamflow_transformed_space" = get_transformed_optimised_streamflow(cmaes_or_dream_result),
      "optimised_modelled_streamflow_realspace" = get_realspace_optimised_streamflow(cmaes_or_dream_result)
    ),
    class = c("result_set", "list")
  )
  
}

# I don't think a helper is required...