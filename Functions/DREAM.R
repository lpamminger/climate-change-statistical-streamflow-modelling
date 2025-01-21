# Related dream functions ------------------------------------------------------

## function defining starting point for dream ==================================
### dream only accepts a function with pars, nseq
converged_parameters <- function(converged_DREAM){
  force(converged_DREAM)
  function(pars, nseq) {
    # pars is list of parameters
    
    # need to convert sequences back into transformed space 
    # Ways to do this:
    # 1. coef(dream) in the DREAM dream function, then add to dream$transformed_coefs (easy)
    # 2. make a function to transform and un-transform using the dream_object
    #    i.e., each parameter is transformed using table in numerical_optimiser_setup
    
    matrix(converged_DREAM$transformed_coefs, nrow = nseq, ncol = length(pars), byrow = TRUE)
  }
}








## Make list of variable ranges for dream ======================================
make_dream_pars <- function(numerical_optimiser_setup) {
  scale <- numerical_optimiser_setup$scale
  dream_pars <- rep(list(c(0, scale)), times = length(numerical_optimiser_setup$parameter_names)) 
  names(dream_pars) <- numerical_optimiser_setup$parameter_names
  return(dream_pars)
  
}


## Make DREAM default parameters ===============================================
control_DREAM_defaults <- function(numerical_optimiser_setup) {
  
  parameter_number <- length(numerical_optimiser_setup$parameter_names)
  chain_number <-  4 * parameter_number
  de_pairs <- floor((chain_number - 1) / 2) 
  parameter_number <- length(numerical_optimiser_setup$parameter_names)
  
  list(
    parameter_number = parameter_number,
    chain_number = chain_number, # is user manually changes chain_number error because DEpairs does not get updated
    warm_up_per_chain = 5E4,
    burn_in_per_chain = 1E4, 
    iterations_after_burn_in_per_chain = 1E4, 
    bound_handling = "rand",
    convergence_gelman = 1.2, 
    check_convergence_steps = 200, # checks for convergences every specified iterations 
    thinning = 1, 
    gamma = 0, 
    nCR = parameter_number, 
    eps = 0.05, 
    steps = 10,
    DEpairs = de_pairs # this is the maximum value it can be
  )
  
}


## A check to see if the users DREAM_controls are acceptable ===================
check_user_DREAM_controls <- function(DREAM_controls) {
  
  # Probably should put a type check here
  
  # Update DEpairs if chain_number was altered (TODO)
  
  # Specific parameter checks
  if (DREAM_controls$chain_number < DREAM_controls$parameter_number) {
    stop(
      "chain_number must be greater than number of parameters in numerical_optmiser_setup",
      call. = FALSE
    )
  }
  
  if(DREAM_controls$warm_up_per_chain <= 0) {
    stop(
      "warm_up_per_chain cannot be less than zero",
      call. = FALSE
    )
  }
  
  if(DREAM_controls$burn_in_per_chain <= 0) {
    stop(
      "burn_in_per_chain cannot be less than zero",
      call. = FALSE
    )
  }
  
  if(DREAM_controls$iterations_after_burn_in_per_chain <= 0) {
    stop(
      "iterations_after_burn_in_per_chain cannot be less than zero",
      call. = FALSE
    )
  }
  
  if(DREAM_controls$convergence_gelman < 1) {
    stop(
      "convergence_gelman cannot be less than one",
      call. = FALSE
    )
  }
  
  if(DREAM_controls$check_convergence_steps < 0) {
    stop(
      "check_convergence_steps cannot be less than zero",
      call. = FALSE
    )
  }
  
  if(DREAM_controls$thinning < 0) {
    stop(
      "thinning cannot be less than zero",
      call. = FALSE
    )
  }
  
  if(DREAM_controls$DEpairs <= 0) {
    stop(
      "DEpairs cannot be less than zero",
      call. = FALSE
    )
  }
  
  return(DREAM_controls)
  
}




## converting DREAM_controls to dream_controls =================================
convert_to_dream_form <- function(DREAM_controls) {
  
  force(DREAM_controls)
  
  # validate user inputs
  DREAM_controls <- check_user_DREAM_controls(DREAM_controls)
  
  # Adjust from DREAM to dream form
  total_iterations_per_chain <- DREAM_controls$burn_in_per_chain + DREAM_controls$iterations_after_burn_in_per_chain
  
  if(DREAM_controls$thinning == 1) {DREAM_controls$thinning <- NA}
  
  
  list(
    ndim = DREAM_controls$parameter_number, 
    nseq = DREAM_controls$chain_number,
    ndraw = total_iterations_per_chain * DREAM_controls$chain_number, # maximum function evaluations (this must be across all chains)
    Rthres = DREAM_controls$convergence_gelman, 
    boundHandling = DREAM_controls$bound_handling, # Method used to handle parameter values outside of parameter bounds. One of: "reflect", "bound", "fold", "none","rand"
    thin.t = DREAM_controls$thinning, 
    burnin.length = DREAM_controls$warm_up_per_chain * DREAM_controls$chain_number, # warm-up period. Not included in chains
    gamma = DREAM_controls$gamma, # (beta parameter?) jump rate - try to improve mixing (Vrugt mentions papers using 0.25 to 0.5)
    REPORT = DREAM_controls$check_convergence_steps, # How often to check Rthres is met. Lower number = more checks = more likely to be converged and high computational cost
    nCR = DREAM_controls$nCR,
    eps = DREAM_controls$eps,
    steps = DREAM_controls$steps, 
    DEpairs = DREAM_controls$DEpairs  # max number chain pair proposal 
  )
  
}


get_sequences <- function(dream_object) {
  dream_object$Sequences
}



remove_burnin_from_trace <- function(dream_object) {
  
  # dream can convergence during burn-in period. If this occurs do not 
  # shorten based on burn-in period. 
  # Take half of iterations otherwise we get the initial wiggles
  if(dream_object$in.burnin) {
    sequences <- get_sequences(dream_object)
    return(window(sequences, start = end(sequences) / 2))
  }
  
  sequences <- get_sequences(dream_object)
  burn_in <- dream_object$controls$burn_in_per_chain
  
  window(sequences, start = (burn_in + 1))
}



mcmc_list_to_tibble <- function(dream_object) {
  
  burnt_in_trace <- remove_burnin_from_trace(dream_object)
  
  map(
    .x = burnt_in_trace,
    .f = unclass
  ) |>
    map(.f = as_tibble) |>
    map2(
      .y = seq(1, length(burnt_in_trace)),
      .f = add_column,
      .before = 1,
    ) |>
    map(
      .f = `colnames<-`,
      value = c("chain", dream_object$numerical_optimiser_setup$parameter_names)
    ) |>
    list_rbind() |> # not sure to include this or not
    mutate(
      n = row_number(),
      .by = chain
    ) |>
    pivot_longer(
      cols = !c(chain, n),
      names_to = "parameter",
      values_to = "parameter_value"
    )
}



# DREAM function (wrapper around dream) ----------------------------------------

DREAM <- function(input, controls) {
  
  ## Check type ================================================================
  ### input can either be a numerical_optimiser_setup object
  ### or a previously run dream object
  check_type <- sloop::s3_class(input)[1]
  stopifnot(check_type %in% c("dream", "numerical_optimiser_setup"))
  
  
  if (check_type == "dream") {
    # get the numerical_optimiser_setup from dream object
    numerical_optimiser_setup <- input$numerical_optimiser_setup
    
    # use converged points as initial starting point
    initial_sampling_function <- converged_parameters(input)
    
  } else {
    initial_sampling_function <- LHSInit # can be in controls (add later)
    numerical_optimiser_setup <- input
  }
  
  
  ## Dream maximises LL. Makes sure minimise_likelihood is FALSE ===============
  stopifnot(!numerical_optimiser_setup$minimise_likelihood)
  
  
  
  ## Check controls ============================================================
  DREAM_controls <- modifyList( # replace defaults with user specified values
    control_DREAM_defaults(numerical_optimiser_setup), 
    val = controls
  ) 
  
  dream_controls <- convert_to_dream_form(DREAM_controls)
  
  
  
  ## Make pars which is a list of variable ranges ==============================
  dream_pars <- make_dream_pars(numerical_optimiser_setup)
  
  
  ## Run dream =================================================================
  dream_result <- dream(
    FUN = numerical_optimiser_setup$ready_objective_function,
    func.type = "logposterior.density",
    pars = dream_pars,
    control = dream_controls,
    INIT = initial_sampling_function # function or LHSInit()
  )
  
  
  ## Add numerical optimiser object to dream object ============================
  dream_result$numerical_optimiser_setup <- numerical_optimiser_setup
  
  
  ## Add controls to dream object ==============================================
  dream_result$controls <- dream_controls 
  
  ## TEMPORARY - save transformed best coefs for re-inputting into dream =======
  dream_result$transformed_coefs <- coef(dream_result, method = "sample.ml")
  
  
  ## Convert the sequences from transformed to real-space ======================
  # future does not like this function. This is because its a generic
  # IT IS BAD PRACTICE HAVING LOCAL S3 METHODS. MAKE PACKAGE.
  transform_to_realspace_dream( 
    mcmc_result = dream_result,
    transform_functions = numerical_optimiser_setup$transform_parameter_methods, 
    lower_bounds = numerical_optimiser_setup$lower_bound, 
    upper_bounds = numerical_optimiser_setup$upper_bound, 
    scale = numerical_optimiser_setup$scale
  )
  
}


# DREAM summary statistics -----------------------------------------------------
get_convergence_statistics <- function(dream_object) {
  
  # dream calculates gelman using 50 % of burn-in sequences if burn-in is complete
  
  acceptance_rate <- 1 - rejectionRate(dream_object$Sequences) 
  
  gelman_statistic <- gelman.diag(
    dream_object$Sequences,
    autoburnin = FALSE
  )
  
  list(
    gauge = dream_object$numerical_optimiser_setup$catchment_data$gauge_ID,
    model = dream_object$numerical_optimiser_setup$streamflow_model()$name,
    parameter = dream_object$numerical_optimiser_setup$parameter_names, 
    exit_message = dream_object$EXITMSG,
    acceptance_rate = acceptance_rate, 
    gelman_statistic = gelman_statistic[[1]][,1] 
  ) |> 
    as_tibble()
  
}


# DREAM plots ------------------------------------------------------------------
gg_trace_plot <- function(dream_object) {
  
  dream_object |>
    mcmc_list_to_tibble() |>
    ggplot(
      aes(
        x = n,
        y = parameter_value,
        colour = as.factor(chain)
      )
    ) +
    geom_line(show.legend = FALSE) +
    labs(
      x = "Iterations",
      y = "Parameter Value",
      title = paste0(
        "Gauge: ",
        dream_object$numerical_optimiser_setup$catchment_data$gauge_ID,
        "\n",
        "Model: ",
        dream_object$numerical_optimiser_setup$streamflow_model()$name
      ),
      subtitle = dream_object$EXITMSG
    ) +
    theme_bw() +
    facet_wrap(~parameter, scales = "free") +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  
}


gg_distribution_plot <- function(dream_object) {
  dream_object |>
    mcmc_list_to_tibble() |>
    ggplot(aes(x = parameter_value)) +
    geom_histogram(
      binwidth = binwidth_bins(30),
      fill = "darkgrey",
      colour = "black"
    ) +
    labs(
      x = "Parameter Value",
      y = "Count",
      title = paste0(
        "Gauge: ",
        dream_object$numerical_optimiser_setup$catchment_data$gauge_ID,
        "\n",
        "Model: ",
        dream_object$numerical_optimiser_setup$streamflow_model()$name
      )
    ) +
    theme_bw() +
    facet_wrap(~parameter, scales = "free") +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
}
