make_default_cmaes_parameters <- function(SCALE, PARAMETER_NUMBER) {
  
  ## set up custom conditions for cmaes ========================================
  ### TolX stopping condition ##################################################
  myStopOnTolX = function(tol = 1E-8) {
    checkmate::assertNumber(tol, na.ok = FALSE) # this bit is changed
    force(tol)
    return(cmaesr::makeStoppingCondition(
      name = "myTolX",
      message = sprintf("Standard deviation below tolerance in all coordinates."),
      stop.fun = function(envir = parent.frame()) {
        return(all(envir$D < tol) && all((envir$sigma * envir$p.c) < tol))
      }
    ))
    
  }
  
  ### Flat fitness stopping condition ##########################################
  stopOnFlatFitness <- cmaesr::makeStoppingCondition(
    name = "flatFitness",
    message = "Stop if fitness value is flat i.e., fitn vector contains all the same values",
    stop.fun = function(envir){
      return(envir$fitn.ordered[1L] == envir$fitn.ordered[ceiling(0.7 * envir$lambda)])
    },
    code = "flatFitness"
  ) 
  
  ## Control parameters output =================================================
  list(
    sigma = SCALE / 3, # hydroState recommends 1/3 of parameter space. # HARD CODED - should extract from optimiser set
    lambda = 4 * floor(3 * log(PARAMETER_NUMBER)), # taken from hydroState
    mu = floor((4 * floor(3 * log(PARAMETER_NUMBER))) / 2), # taken from hydroState
    restart.triggers = c(
      "flatFitness", 
      "indefCovMat", 
      "noEffectAxis", 
      "noEffectCoord",
      "myTolX"
    ), 
    max.restarts = 4L,  
    restart.multiplier = 2L, 
    stop.ons = c(
      list(
        stopOnFlatFitness, 
        myStopOnTolX(tol = 1E-5), 
        cmaesr::stopOnNoEffectAxis(),
        cmaesr::stopOnNoEffectCoord(),
        cmaesr::stopOnCondCov(),
        stopOnTimeBudget(budget = 180)
      ) 
    )
  )
}


my_cmaes <- function(numerical_optimiser_setup, cmaes_control = NULL, print_monitor = FALSE) {
  
  ## Check if numerical_optimiser_setup and its minimising the likelihood ==================
  stopifnot(sloop::s3_class(numerical_optimiser_setup)[1] == "numerical_optimiser_setup")
  stopifnot(numerical_optimiser_setup$minimise_likelihood)
  
  
  ## cmaes requires smoof objective function ===================================
  SCALE <- numerical_optimiser_setup$scale
  
  param_smoof <- map(
    .x = numerical_optimiser_setup$parameter_names,
    .f = ParamHelpers::makeNumericParam,
    lower = 0, # Scale is always between 0 and upper limit defined by user
    upper = SCALE
  )
  
  param_smoof <- ParamHelpers::makeParamSet(params = param_smoof)
  
  cmaes_objective_function <- smoof::makeSingleObjectiveFunction(
    fn = numerical_optimiser_setup$ready_objective_function, 
    par.set = param_smoof,
    vectorized = TRUE # vectorized goes brrrrrrrrr
  )
  
  
  ### Default cmaes_control ####################################################
  PARAMETER_NUMBER <- length(numerical_optimiser_setup$parameter_names)
  
  if (is.null(cmaes_control)) { 
    cmaes_control <- make_default_cmaes_parameters(SCALE, PARAMETER_NUMBER)
  }

  
  ### Make monitor #############################################################
  if(print_monitor) {
    create_monitor <- cmaesr::makeSimpleMonitor()
  } else {
    create_monitor <- NULL
  }

  
  ## Run cmaes =================================================================
  cmaes_result <- cmaesr::cmaes(
    objective.fun = cmaes_objective_function,
    start.point = rep((SCALE * 0.5), times = PARAMETER_NUMBER),
    monitor = create_monitor,
    control = cmaes_control
  )
  
  ## Add the numerical_optimiser_setup to cma_result object
  cmaes_result$numerical_optimiser_setup <- numerical_optimiser_setup
  
  # Remove anything that is not the cmaes_result
  rm(list = c("numerical_optimiser_setup", "SCALE", "param_smoof", "cmaes_objective_function", "cmaes_control"))
  gc()
  return(cmaes_result)
  
}
