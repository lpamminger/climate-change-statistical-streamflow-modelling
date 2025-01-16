## Fitting streamflow models to catchments
cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(
  tidyverse,
  dream,
  coda, 
  lattice, 
  tictoc, 
  truncnorm, 
  sloop
)




# Import and prepare data-------------------------------------------------------

## Import annual streamflow, precip, CO2 and gauge data ========================
start_stop_indexes <- readr::read_csv(
  "./Data/Tidy/start_end_index.csv",
  show_col_types = FALSE
)

data <- readr::read_csv(
  "./Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  col_types = "icdddddddll", # ensuring each column is a the correct type
  show_col_types = FALSE
)

gauge_information <- readr::read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE
)

CMAES_results <- read_csv(
  "./Results/my_cmaes/CMAES_parameter_results_20241130.csv",
  show_col_types = FALSE
) |>
  filter(objective_function != "CO2_variable_objective_function") # temporary solution

## Utility functions ===========================================================
source("./Functions/utility.R")


## Import streamflow functions =================================================
source("./Functions/streamflow_models.R")
source("./Functions/parameter_transformations.R")
source("./Functions/catchment_data_blueprint.R")
source("./Functions/cmaes_dream_summaries.R")
source("./Functions/objective_functions.R")
source("./Functions/numerical_optimiser_setup.R")
source("./Functions/generic_functions.R")
source("./Functions/my_dream.R")
source("./Functions/objective_function_setup.R")
source("./Functions/result_set.R")



# What I have so far:
# 1. DREAM function 
# 2. gg_trace_plot function

# What I need:
# 1. a method for chunking i.e., run DREAM many times (save dream objects),
#    check convergence using gg_trace_plot, put converged dream objects 
#    back into DREAM function with altered controls

# I have all functions required to complete the DREAM analysis. 
# I need to integrate them in and find a way to run them many times.
# Ideas:

# 1. Run a chunk until convergence. Check trace. If good then run to get nice distributions
#### This is more involved. Maybe test a single chunk (after testing a single catchment) like this.
#### Once I am satisfied it works with the test chunk then move to idea 2.

# 2. Run a all chunks until convergence. Check all traces. If good run all to get nice distributions. 





# Q. What I want to happen?
# A. I want to build adequate distribution for each parameter

# Q. How will I do this?
# A. DREAM can be used to construct probability distributions for each parameter

# Q. What is the best way to use DREAM?
# A. I have found two DREAM runs produces the best outcome. The first
#    run involves getting the chains to converge. The second run takes the 
#    converged chains and produces enough samples to create an adequate 
#    distribution


# Q. What is the cleanest way to code this?
# A. Repeat the my_dream function 2 times.
#    Alter my_dream to take either numerical_optimiser_object or dream_object as
#    input to make feeding results from one to another easier. Both require
#    different control parameters. Separate them.


# Q. What will the function(s) do?
# A. I think it must be at least a two step process. The first step is getting
#    dream to converge. To check convergence use plot (best), gelman, acceptance rate 
#    and autocorrelation. If it were just number we could make it a one step
#    process. However, I think the trace plots it the most important indicator.
#    Second step involves building the parameter distributions by cranking 
#    iterations.


# Q. The dream controls are confusing. Can you rewrite them?
# A. Write a wrapper around dream controls with easy to understand inputs




# Get best parameters from converged dream object to be ready to be put back into dream

# Related dream functions ------------------------------------------------------
converged_parameters <- function(converged_DREAM){
  function(pars, nseq) {
    matrix(coef(converged_DREAM), nrow = nseq, ncol = length(pars), byrow = TRUE)
  }
}


make_dream_pars <- function(numerical_optimiser_setup) {
  # list of variable ranges 
  scale <- numerical_optimiser_setup$scale
  dream_pars <- rep(list(c(0, scale)), times = length(numerical_optimiser_setup$parameter_names)) 
  names(dream_pars) <- numerical_optimiser_setup$parameter_names
  return(dream_pars)
  
}

# Make DREAM default paramaters
control_DREAM_defaults <- function(numerical_optimiser_setup) {
  
  parameter_number <- length(numerical_optimiser_setup$parameter_names)
  chain_number <-  (2 * parameter_number) + 1
  
  list(
    parameter_number = parameter_number,
    chain_number = chain_number, 
    warm_up_per_chain = 5E4,
    burn_in_per_chain = 1E4, 
    iterations_after_burn_in_per_chain = 1E4, 
    bound_handling = "rand",
    convergence_gelman = 1.2, 
    check_convergence_steps = 200, # checks for convergences every specified iterations 
    thinning = 1, 
    gamma = 1, # idk
    nCR = 3, # idk
    eps = 0.05, #  idk
    steps = 10,# idk
    DEpairs = floor((chain_number - 1) / 2) 
  )
  
}


check_user_DREAM_controls <- function(DREAM_controls) {
  
  # Probably should put a type check here
  
  
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


# I think this is doing what it is supposed to be doing
# To test crank up the controls

DREAM <- function(input, controls) {
  
  ## Check type ================================================================
  ### input can either be a numerical_optimiser_setup object
  ### or a previously run dream object
  check_type <- s3_class(input)[1]
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
  dream_result$controls <- controls
  
  ## Convert the sequences from transformed to real-space ======================
  transform_to_realspace(
    mcmc_result = dream_result,
    transform_functions = numerical_optimiser_setup$transform_parameter_methods, 
    lower_bounds = numerical_optimiser_setup$lower_bound, 
    upper_bounds = numerical_optimiser_setup$upper_bound, 
    scale = numerical_optimiser_setup$scale
  )
  

}







# Summarise this into a single list then unpack

# build_parameter_distribution controls
# parameters
# parameter_number, chain_number, warmup, burnin, lengthafter burnin,
# thinning, convergence_gelman, bound_handling, gamma, report, nCR, eps, steps, DEpairs, burn_in_length


# Order of parameter checks
# 1. user inputs desired parameters
# 2. compare with defaults. User overwrites default parameters (in DREAM function)
# 3. makes sure the users parameters are feasible using check_user_dream_controls
# 4. convert into form ready for dream

# test
user_DREAM_control <- list( # this will the list in DREAM(input, control) - change required for DREAM
  warm_up_per_chain = 100, #
  burn_in_per_chain = 100,
  iterations_after_burn_in_per_chain = 100,
  convergence_gelman = 1.2,
  check_convergence_steps = 0, # check convergence every 200 steps
  thinning = 1
)


most_complex_gauge <- "218005"
example_DREAM <- most_complex_gauge |>
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  ) |>
  numerical_optimiser_setup_vary_inputs(
    streamflow_model = streamflow_model_drought_separate_shifted_CO2_seasonal_ratio_auto,
    objective_function = constant_sd_objective_function,
    bounds_and_transform_method = make_default_bounds_and_transform_methods(),
    minimise_likelihood = FALSE
  ) |> 
  DREAM(controls = user_DREAM_control)

# ERROR
# example_DREAM_again <- DREAM(example_DREAM, control = user_DREAM_control) # this causes an error. Hopefully it resolves itself when I do it properly
# "Error in if (!any(delta.tot > 0)) stop("AdaptpCR: no changes in X, would cause NaN in pCR") : missing value where TRUE/FALSE needed"

# Check convergence criteria
# Gelman should be converged because its build in
# Trace diagrams plotted
# Acceptance rate
# Autocorrelation?

# I can make custom plots using ggplot or use mcmc?

# Q. Do I convert to real-space before plotting?
# A. Yes

# Q. Should the sequences be converted to realspace in the DREAM function?
# A. Should be converted in the DREAM function.


# Q. Which plotting method do I use?
# A. ggplot gives me the most flexibility, mcmc is the easiest. 
#    Commit to ggplot.


# Q. How to extract traces from dream object?
# A. use get_sequences function in generic_functions.R. This gives me a
#    mcmc object. I want a tibble. Need to convert.


# Visually inspect trace diagrams for each catchments.
# Make a mega plot to see if traces are stationary and mixing well.
# If okay then build distributions using converged chains. Use previous code by
# adjust the control parameter of dream


# Order for plotting
# 1. get sequences
# 2. remove burn-in (window method or use burn-in from controls)
# 3. current form is a mcmc.list object. I want it as a tibble. Conversion
#    involves: unclassing, see x <- map... below


get_sequences <- function(dream_object) {
  dream_object$Sequences
}

remove_burnin_from_trace <- function(dream_object) {
  
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
      )
    ) +
    theme_bw() +
    facet_wrap(~parameter, scales = "free") +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
  
}


trace_plot <- gg_trace_plot(dream_object = example_DREAM)
trace_plot


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

dist_plot <- gg_distribution_plot(example_DREAM)
dist_plot

# use the conversion build into dream
1 - rejectionRate(y)










# Draft functions for dream ----------------------------------------------------
most_complex_gauge <- "218005"
example_DREAM <- most_complex_gauge |>
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  ) |>
  numerical_optimiser_setup_vary_inputs(
    streamflow_model = streamflow_model_drought_separate_shifted_CO2_seasonal_ratio_auto,
    objective_function = constant_sd_objective_function,
    bounds_and_transform_method = make_default_bounds_and_transform_methods(),
    minimise_likelihood = FALSE
  ) |> 
  my_dream()


## Wrapper for dream_default parameters ========================================
parameter_number <- length(example_DREAM$parameter_names) # this should a apparent from the function
chain_number <- (2 * parameter_number) + 1
warm_up_per_chain <- 5E4 # not included in sequences
burn_in_per_chain <- 1E4 # remove sequences from chain
length_after_burn_in_per_chain <- 1E4
total_iterations_per_chain <- burn_in_per_chain + length_after_burn_in_per_chain 


dream_control_parameters <- list(
  ndim = parameter_number, 
  nseq = chain_number,
  ndraw = total_iterations_per_chain * chain_number, # maximum function evaluations (this must be across all chains)
  Rthres = 1.2, 
  boundHandling = "rand", # Method used to handle parameter values outside of parameter bounds. One of: "reflect", "bound", "fold", "none","rand"
  thin.t = NA, 
  burnin.length = warm_up_per_chain * chain_number, # warm-up period. Not included in chains
  gamma = 1, # (beta parameter?) jump rate - try to improve mixing (Vrugt mentions papers using 0.25 to 0.5)
  REPORT = total_iterations_per_chain / 4, # How often to check Rthres is met. Lower number = more checks = more likely to be converged and high computational cost
  nCR = 3,
  eps = 0.05,
  steps = 10, 
  DEpairs = (chain_number - 1) / 2  # max number chain pair proposal 
)


## Run until converged =========================================================
set.seed(1)
tic <- as.numeric(Sys.time()) 
test_dream_v2 <- my_dream(
  numerical_optimiser_setup = example_DREAM,
  dream_control = dream_control_parameters
) 
toc <- as.numeric(Sys.time()) - tic

#burn_sequences <- window(test_dream_v2$Sequences, start = burn_in_per_chain) # if converged before burn-in then this will not work
#burn_sequences <- window(test_dream_v2$Sequences, start = end(high_chains_DREAM$Sequences)/2 + 1)
burn_sequences <- test_dream_v2$Sequences

1 - rejectionRate(test_dream_v2$Sequences)
gelman.diag(test_dream_v2$Sequences, autoburnin = FALSE)
gelman.plot(test_dream_v2$Sequences)

# Maybe let it converge then pass the results to init to be ran again without convergence checks?
plot(burn_sequences)


## Using converged results stick it back into DREAM ============================

### Don't use the separate function - use this - easy fixing stuff
my_dream_v2 <- function(numerical_optimiser_setup, converged_DREAM, dream_control = NULL) {
  
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
  dream_parameters <- rep(list(c(0, SCALE)), times = PARAMETER_NUMBER) # bounds
  names(dream_parameters) <- numerical_optimiser_setup$parameter_names
  
  ## Run dream =================================================================
  dream_result <- dream(
    FUN = numerical_optimiser_setup$ready_objective_function,
    func.type = "logposterior.density",
    pars = dream_parameters,
    control = dream_control,
    INIT = converged_parameters,
    INIT.pars = list("converged_DREAM" = converged_DREAM)
  )
  
  
  ## Add numerical optimiser object to dream object ============================
  dream_result$numerical_optimiser_setup <- numerical_optimiser_setup
  return(dream_result)
  
}



converged_parameters <- function(pars, nseq, converged_DREAM) {
  matrix(coef(converged_DREAM), nrow = nseq, ncol = length(pars), byrow = TRUE)
}



# run without convergence 
no_converge_dream_control_parameters <- list(
  ndim = parameter_number, 
  nseq = chain_number,
  ndraw = length_after_burn_in_per_chain * chain_number, # maximum function evaluations (this must be across all chains)
  Rthres = 1.2, 
  boundHandling = "rand", # Method used to handle parameter values outside of parameter bounds. One of: "reflect", "bound", "fold", "none","rand"
  thin.t = 10, 
  burnin.length = 0, # warm-up period. Not included in chains
  gamma = 1, # (beta parameter?) jump rate - try to improve mixing (Vrugt mentions papers using 0.25 to 0.5)
  REPORT = 0, # How often to check Rthres is met
  nCR = 3,
  eps = 0.05,
  steps = 10, 
  DEpairs = (chain_number - 1) / 2  # max number chain pair proposal 
)

test_dream_v3 <- my_dream_v2(
  numerical_optimiser_setup = example_DREAM,
  converged_DREAM = test_dream_v2,
  dream_control = no_converge_dream_control_parameters
)



real_space_DREAM_v3 <- transform_to_realspace(
  mcmc_result = test_dream_v3,
  transform_functions = example_DREAM$transform_parameter_methods, 
  lower_bounds = example_DREAM$lower_bound, 
  upper_bounds = example_DREAM$upper_bound, 
  scale = 100
)

plot(real_space_DREAM_v3$Sequences)
half_burnin_real_real_space_DREAM_v3_sequences <- window(real_space_DREAM_v3$Sequences, start = end(real_space_DREAM_v3$Sequences)/2 + 1)


x <- map(
  .x = real_space_DREAM_v3$Sequences,
  .f = unclass
) |>
  map(.f = as_tibble) |>
  map2(
    .y = seq(1, length(real_space_DREAM_v3$Sequences)),
    .f = add_column,
    .before = 1,
  ) |> 
  map(
    .f = `colnames<-`,
    value = c("Chain", real_space_DREAM_v3$numerical_optimiser_setup$parameter_names)
  ) |> 
  list_rbind()



best_parameters <- get_best_parameters_real_space(test_dream_v3)

x |>
  pivot_longer(
    cols = !Chain,
    names_to = "parameter",
    values_to = "parameter_value"
  ) |>
  ggplot(aes(x = parameter_value)) +
  geom_histogram(binwidth = binwidth_bins(30), fill = "darkgrey", colour = "black") +
  theme_bw() +
  facet_wrap(~parameter, scales = "free")


x |>
  pivot_longer(
    cols = !Chain,
    names_to = "parameter",
    values_to = "parameter_value"
  ) |>
  mutate(
    n = row_number(),
    .by = Chain
  ) |>
  ggplot(
    aes(
      x = n,
      y = parameter_value,
      colour = as.factor(Chain)
    )
  ) +
  geom_line(show.legend = FALSE) +
  theme_bw() +
  facet_wrap(~parameter, scales = "free")
# What am I learning from this. All the CO2 parameters want to be zero...  



stop_here <- 1







# What about start with CMAES best parameters? ---------------------------------
most_complex_gauge <- "218005"
CMAES_best_parameters <- CMAES_results |> 
  filter(gauge == most_complex_gauge) |> 
  filter(streamflow_model == "streamflow_model_drought_separate_shifted_CO2_seasonal_ratio_auto") |> 
  filter(objective_function == "CO2_variable_objective_function") |> 
  select(parameter, parameter_value) |> 
  left_join(
    make_default_bounds_and_transform_methods(),
    by = join_by(parameter)
  ) |> 
  add_column(
    transform_method_reverse = list(
      real_space_to_scaled_linear_transform,
      real_space_to_scaled_linear_transform,
      real_space_to_scaled_logarithmic_transform,
      real_space_to_scaled_linear_transform,
      real_space_to_scaled_linear_transform,
      real_space_to_scaled_linear_transform,
      real_space_to_scaled_linear_transform,
      real_space_to_scaled_logarithmic_transform,
      real_space_to_scaled_logarithmic_transform
    )
  )



transform_test <- function(transform_method, parameter, lower_bound, upper_bound, scale) {
  transform_method <- noquote(transform_method)
  transform_method(parameter, lower_bound, upper_bound, scale)
}


CMAES_coef <- pmap_dbl(
  .l = list(
    CMAES_best_parameters$transform_method_reverse,
    CMAES_best_parameters$parameter_value,
    CMAES_best_parameters$lower_bound,
    CMAES_best_parameters$upper_bound
  ),
  .f = transform_test,
  scale = 100
)

# Stick it into dream ==========================================================
converged_init <- function(pars, nseq, converged_coef) {
  matrix(converged_coef, nrow = nseq, ncol = length(pars), byrow = TRUE)
}



my_dream_CMAES_input <- function(numerical_optimiser_setup, converged_coef_CMAES, dream_control = NULL) {
  
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
  dream_parameters <- rep(list(c(0, SCALE)), times = PARAMETER_NUMBER) # bounds
  names(dream_parameters) <- numerical_optimiser_setup$parameter_names
  
  ## Run dream =================================================================
  dream_result <- dream(
    FUN = numerical_optimiser_setup$ready_objective_function,
    func.type = "logposterior.density",
    pars = dream_parameters,
    control = dream_control,
    INIT = converged_init,
    INIT.pars = list("converged_coef" = converged_coef_CMAES)
  )
  
  
  ## Add numerical optimiser object to dream object ============================
  dream_result$numerical_optimiser_setup <- numerical_optimiser_setup
  return(dream_result)
  
}



## Set up dream ================================================================

### objective function #########################################################
most_complex_gauge <- "218005"

example_DREAM <- most_complex_gauge |>
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  ) |>
  numerical_optimiser_setup_vary_inputs(
    streamflow_model = streamflow_model_drought_separate_shifted_CO2_seasonal_ratio_auto,
    objective_function = CO2_variable_objective_function,
    bounds_and_transform_method = make_default_bounds_and_transform_methods(),
    minimise_likelihood = FALSE
  ) 

parameter_number <- length(example_DREAM$parameter_names) # this should a apparent from the function
chain_number <- (2 * parameter_number) + 1
warm_up_per_chain <- 5E3 # not included in sequences
burn_in_per_chain <- 0 # remove sequences from chain
length_after_burn_in_per_chain <- 1E3
total_iterations_per_chain <- burn_in_per_chain + length_after_burn_in_per_chain 


no_converge_dream_control_parameters <- list(
  ndim = parameter_number, 
  nseq = chain_number,
  ndraw = length_after_burn_in_per_chain * chain_number, # maximum function evaluations (this must be across all chains)
  Rthres = 1.2, 
  boundHandling = "rand", # Method used to handle parameter values outside of parameter bounds. One of: "reflect", "bound", "fold", "none","rand"
  thin.t = 10, 
  burnin.length = 0, # warm-up period. Not included in chains
  gamma = 1, # (beta parameter?) jump rate - try to improve mixing (Vrugt mentions papers using 0.25 to 0.5)
  REPORT = 0, # How often to check Rthres is met
  nCR = 3,
  eps = 0.05,
  steps = 10, 
  DEpairs = (chain_number - 1) / 2  # max number chain pair proposal 
)

test_dream_v3 <- my_dream_CMAES_input(
  numerical_optimiser_setup = example_DREAM,
  converged_coef_CMAES = CMAES_coef,
  dream_control = no_converge_dream_control_parameters
)

plot(test_dream_v3$Sequences)
gelman.diag(test_dream_v3$Sequences, autoburnin = FALSE)
1 - rejectionRate(test_dream_v3$Sequences)

