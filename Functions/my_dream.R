make_default_dream_parameters <- function(PARAMETER_NUMBER) {
  list(
    ndim = PARAMETER_NUMBER, # number of parameters in model
    nseq = 4 * PARAMETER_NUMBER, # number of chains to evolve (number of parameters * factor)
    ndraw = 200000, # CHANGE TO 500,000 or 1e7 
    Rthres = 1.2, # Vrugt recommendation
    boundHandling = "rand", # Method used to handle parameter values outside of parameter bounds. One of: "reflect", "bound", "fold", "none","rand"
    thin.t = 1, # CHANGE TO 32. Save the memory. 
    burnin.length = 0 # hopefully this removes some of the burn-in errors. With a 0.9 (90 % of evaluation using burn in still error). Set to zero no burn in 
    #REPORT = 1000
  )
}




my_dream <- function(numerical_optimiser_setup, dream_control = NULL) {
  
  ## Dream maximises makes sure minimise_likelihood is FALSE
  stopifnot(s3_class(numerical_optimiser_setup)[1] == "numerical_optimiser_setup")
  stopifnot(!numerical_optimiser_setup$minimise_likelihood)
  
  # Use default parameters if dream_control is not specified
  PARAMETER_NUMBER <- length(numerical_optimiser_setup$parameter_names)
  
  if (is.null(dream_control)) {
    dream_control <- make_default_dream_parameters(PARAMETER_NUMBER)
  }
  
  ### Make pars which is a list of variable ranges #############################
  SCALE <- numerical_optimiser_setup$scale
  dream_parameters <- rep(list(c(0, SCALE)), times = PARAMETER_NUMBER)
  names(dream_parameters) <- numerical_optimiser_setup$parameter_names
  
  ## Run dream =================================================================
  dream_result <- dream(
    FUN = numerical_optimiser_setup$ready_objective_function,
    func.type = "logposterior.density",
    pars = dream_parameters,
    control = dream_control
  )
  
  ## Add numerical optimiser object to dream object ============================
  dream_result$numerical_optimiser_setup <- numerical_optimiser_setup
  return(dream_result)
  
}


