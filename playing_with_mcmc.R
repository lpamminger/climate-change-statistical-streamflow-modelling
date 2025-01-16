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
  fmcmc,
  sloop#,
  #BayesianTools,
  #lhs
)
# install dream using: install.packages("dream", repos="http://R-Forge.R-project.org")
# install.packages("coda")



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
)

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


# Identify best model combination for DREAM gauges -----------------------------
best_model_combination_per_catchment <- CMAES_results |>
  select(!c(parameter, parameter_value, optimiser, loglikelihood, exit_message, near_bounds)) |>
  slice_min(
    AIC,
    by = gauge
  ) |> 
  distinct() |> 
  select(!AIC) |> 
  rename(
    streamflow_model_name = streamflow_model,
    objective_function_name = objective_function
  ) 


# Demonstrate how to use fmcmc -------------------------------------------------
# For the demonstration I will use select
# 212202 from the best_model_combination_per_catchment
# This is because it has the simplest streamflow_model and objective_function
# combination



## Step 1. Build log-likelihood function =======================================
gauge <- "212209"



example_setup <- gauge |> # this is using my numerical_optimiser_setup object. Not important
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  ) |>
  numerical_optimiser_setup_vary_inputs(
    streamflow_model = streamflow_model_precip_only,
    objective_function = constant_sd_objective_function,
    bounds_and_transform_method = make_default_bounds_and_transform_methods(),
    minimise_likelihood = FALSE
  ) 


example_likelihood_function <- example_setup$ready_objective_function




## Step 2. Run MCMC with default parameters and inspect the trace diagram ======
set.seed(1)
simplest_mcmc <- MCMC(
  initial = c(77, 50, 90), # initial starting point 
  fun = example_likelihood_function,
  nsteps = 1E3
)

plot(simplest_mcmc)

### Problems:
### 1. For par1 and par3 the initial part of the trace diagrams are very 
###    different to the rest of the diagram. This suggest the initial distribution
###    is very different to the target distribution. It does slowly converge.
### 2. All pars have poor mixing (flat lines). We want the distribution to 
###    sample a range of parameters when building our posterior. Stair-stepping
###    indicates getting stuck on one part of the distribution. This may indicate
###    it is jumping too far.
### 3. No bounds means the log-likelihood function can produce NA


## Step 3. Addressing issues with simplest_mcmc ================================
### Addressing problems:
### 1. introduce a burn-in period
### 2. use another kernel and/or reduce the scale parameter to avoid it jumping
###    to far. 
### 3. add bounds


improved_mcmc <- MCMC(
  initial = simplest_mcmc,
  fun = example_likelihood_function,
  nsteps = 1E3,
  burnin = 100, # problem 1
  kernel = kernel_normal_reflective(
    scale = 0.25, # problem 2
    ub = 100, # problem 3
    lb = 0
  )
)

plot(improved_mcmc)
### Inspect improved_mcmc

### Step 3.1 - New Problems: ###################################################
### 1. Better mixing, however all pars look like there is high autocorrelation 
###    in the trace diagrams. To fix increase nsteps. Increasing thinning can 
###    also decrease autocorrelation.

improved_v2_mcmc <- MCMC(
  initial = improved_mcmc,
  fun = example_likelihood_function,
  nsteps = 5E4, # new problem 1
  burnin = 1E3, 
  thin = 10, # new problem 1
  kernel = kernel_normal_reflective(
    scale = 0.25, 
    ub = 100, 
    lb = 0
  )
)


plot(improved_v2_mcmc)
### Inspect improved_v2_mcmc
### Step 3.2 - New Problems: ###################################################
### 1. par1 and par2 are still not stationary. Try a different kernel

improved_v3_mcmc <- MCMC(
  initial = improved_mcmc,
  fun = example_likelihood_function,
  kernel = kernel_adapt(
    ub = 100,
    lb = 0,
    eps = 0.5 # this makes it mix better. Found using trial and error
  ),
  nsteps = 1E4, 
  burnin = 1E3 
)

plot(improved_v3_mcmc)
### Yay its mixing well and stationary


## Step 4. Checking convergence of trace diagrams ==============================
converged_mcmc <- MCMC(
  initial = improved_mcmc,
  fun = example_likelihood_function,
  kernel = kernel_adapt(
    ub = 100,
    lb = 0,
    eps = 0.5 # this makes it mix better. Found using trial and error
  ),
  nsteps = 1E4, 
  burnin = 1E3,
  nchains = 2,
  conv_checker = convergence_gelman(freq = 1E3L, threshold = 1.2)
)

check_acceptance_rate <- 1 - rejectionRate(converged_mcmc)



### Gelman and Rubin statistic is < 1.2 indicating convergence (Vrugt)
### Acceptance rate is between 15 % to 35 % (Vrugt). This is acceptable
### Given both metrics are meet the distribution has converged


## Step 5. Build posterior distribution for sampling ===========================
### Re-run the converged_mcmc with many iterations to build adequate 
### posterior distribution. No convergence checks are required.

final_mcmc <- MCMC(
  initial = converged_mcmc,
  fun = example_likelihood_function,
  kernel = kernel_adapt(
    ub = 100,
    lb = 0,
    eps = 0.5 # this makes it mix better. Found using trial and error
  ),
  nsteps = 5E4, 
  burnin = 1E3,
  nchains = 2,
  thin = 10
)

plot(final_mcmc)

### Side-note: #################################################################
# - I have transformed each parameter between 0 and 100
# - Scaling improves search space? - CHECK WORDING  

### To convert from transformed space to real space using this functions:

single_col_mcmc_transform_to_realspace <- function(single_col_mcmc_result, transform_function, lower_bound, upper_bound, scale) {
  
  transform_function <- noquote(transform_function)
  
  transform_function(
    parameter_set = single_col_mcmc_result,
    lower_bound = lower_bound,
    upper_bound = upper_bound,
    scale = 100
  )
  
}

transform_to_realspace <- function(x, ...) { # generic ... allows extra arguments
  UseMethod("transform_to_realspace")
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

transform_to_realspace.dream <- function(mcmc_result, transform_functions, lower_bounds, upper_bounds, scale) {

  transform_to_realspace(
    mcmc_result = mcmc_result[["Sequences"]],
    transform_functions = transform_functions,
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    scale = scale
  )
  
}



## Step 6. Final result ========================================================
final_mcmc_real_space <- transform_to_realspace(
  mcmc_result = final_mcmc,
  transform_functions = c(linear_parameter_transform, logarithmic_parameter_transform, logarithmic_parameter_transform), 
  lower_bounds = example_setup$lower_bound, 
  upper_bounds = example_setup$upper_bound, 
  scale = 100
  )


plot(final_mcmc_real_space)
summary(final_mcmc_real_space)

# Compare results to CMAES
CMAES_results |> 
  filter(gauge == {{ gauge }}) |> 
  filter(streamflow_model == "streamflow_model_precip_only") |> 
  filter(objective_function == "constant_sd_objective_function")






# DREAM example ----------------------------------------------------------------

## Number of chains are very important for mixing ==============================

### Low number of chains #######################################################
low_chains_control_parameters <- list(
  ndim = 3, # number of parameters in model
  nseq = 3, # number of chains to evolve (number of parameters * 2 + 1)
  ndraw = 3.5E2, # maximum function evaluations (this must be across all chains)
  Rthres = 1.01, # Vrugt recommendation (1.2) - convergence criteria
  boundHandling = "rand", # Method used to handle parameter values outside of parameter bounds. One of: "reflect", "bound", "fold", "none","rand"
  thin.t = NA, 
  burnin.length = 5E2, # hopefully this removes some of the burn-in errors. With a 0.9 (90 % of evaluation using burn in still error). Set to zero no burn in
  REPORT = 0
)


low_chains_DREAM <- gauge |>
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  ) |>
  numerical_optimiser_setup_vary_inputs(
    streamflow_model = streamflow_model_precip_only,
    objective_function = constant_sd_objective_function,
    bounds_and_transform_method = make_default_bounds_and_transform_methods(),
    minimise_likelihood = FALSE
  ) |>
  my_dream(dream_control = low_chains_control_parameters) 


low_chains_sequence <- window(low_chains_DREAM$Sequences, start = end(low_chains_DREAM$Sequences)/2 + 1)
plot(low_chains_sequence)



### High number of chains ######################################################
high_chains_control_parameters <- list(
  ndim = 3, # number of parameters in model
  nseq = (2 * 3) + 1, # number of chains to evolve (number of parameters * 2 + 1)
  ndraw = 3.5E4, # maximum function evaluations (this must be across all chains)
  Rthres = 1.01, # Vrugt recommendation (1.2) - convergence criteria
  boundHandling = "rand", # Method used to handle parameter values outside of parameter bounds. One of: "reflect", "bound", "fold", "none","rand"
  thin.t = NA, 
  burnin.length = 5E3, # hopefully this removes some of the burn-in errors. With a 0.9 (90 % of evaluation using burn in still error). Set to zero no burn in
  REPORT = 0
)



high_chains_DREAM <- gauge |>
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  ) |>
  numerical_optimiser_setup_vary_inputs(
    streamflow_model = streamflow_model_precip_only,
    objective_function = constant_sd_objective_function,
    bounds_and_transform_method = make_default_bounds_and_transform_methods(),
    minimise_likelihood = FALSE
  ) |>
  my_dream(dream_control = high_chains_control_parameters) 


high_chains_sequence <- window(high_chains_DREAM$Sequences, start = end(high_chains_DREAM$Sequences)/2 + 1)
plot(high_chains_sequence)







## Final result for DREAM using trial and error ================================
gauge <- "212209"  
dream_control_parameters <- list(
  ndim = 3, # number of parameters in model
  nseq = (3 * 2) + 1, # number of chains to evolve (number of parameters * 2 + 1)
  ndraw = 2E5, # maximum function evaluations (this must be across all chains)
  Rthres = 1.2, # Vrugt recommendation - convergence criteria
  boundHandling = "rand", # Method used to handle parameter values outside of parameter bounds. One of: "reflect", "bound", "fold", "none","rand"
  thin.t = NA, # 10 this messes up the plots
  burnin.length = 1.5E5, # hopefully this removes some of the burn-in errors. With a 0.9 (90 % of evaluation using burn in still error). Set to zero no burn in
  gamma = 1, # (beta parameter?) jump rate - try to improve mixing (Vrugt suggests 0.25 to 0.5)
  REPORT = 1E3,
  nCR = 3,
  steps = 10 # maybe this alters steps after convergence?
  #DEpairs = 3 # this is not a parameter. You are forced to use the default
)


set.seed(1)
example_DREAM <- gauge |>
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  ) |>
  numerical_optimiser_setup_vary_inputs(
    streamflow_model = streamflow_model_precip_only,
    objective_function = constant_sd_objective_function,
    bounds_and_transform_method = make_default_bounds_and_transform_methods(),
    minimise_likelihood = FALSE
  ) |>
  my_dream(dream_control = dream_control_parameters) 

plot(example_DREAM) # steps = 10 -> 550, 100 -> 3500, 500 -> 




# diy
diy_example_DREAM <- example_DREAM

diy_example_DREAM$Sequences <- transform_to_realspace(
  mcmc_result = example_DREAM$Sequences,
  transform_functions = c(linear_parameter_transform, logarithmic_parameter_transform, logarithmic_parameter_transform), 
  lower_bounds = example_setup$lower_bound, 
  upper_bounds = example_setup$upper_bound, 
  scale = 100
)

plot(diy_example_DREAM$Sequences) # why are these different
plot(diy_example_DREAM) # why is this different to ^

# I don't know how window works
ss <- window(diy_example_DREAM$Sequences, start = end(diy_example_DREAM$Sequences)/2 + 1) # how does this represent the burn-in period?
plot(ss)

1 - rejectionRate(ss)
# ndraw = 150,000, with thin.t = 10 -> 15,000
# burn-in is 100,000
# window(...) removes the burn-in period
# it is not straight forward because ndraw and burn-in is for everything not
# for a single chain
# The number of iterations shown using ss is:
# ndraw - burnin.length / nseq


# DREAM diy --------------------------------------------------------------------
# Try most complex model combination
# Ignore convergence. Just run it with burn-in
# Longer burin and more chains are probably required.

# The DREAM R packages is not good
# a5 parameter (par7) unconverging - try boundHandleing "bound"?
# I would like a runs after convergence parameter

most_complex_gauge <- "218005"

dream_control_parameters <- list(
  ndim = 3, # number of parameters in model
  nseq = 3 * 9, # number of chains to evolve (number of parameters * 2 + 1). Works well for 2*3 + 1 -> 3 * number_parameters
  ndraw = 3E5, # maximum function evaluations (this must be across all chains)
  Rthres = 1.2, # I don't want it to converge. Set very close to 1
  boundHandling = "rand", # Method used to handle parameter values outside of parameter bounds. One of: "reflect", "bound", "fold", "none","rand"
  thin.t = 10, # 10 this messes up the plots
  burnin.length = 2E5, # hopefully this removes some of the burn-in errors. With a 0.9 (90 % of evaluation using burn in still error). Set to zero no burn in
  gamma = 1, # (beta parameter?) jump rate - try to improve mixing (Vrugt suggests 0.25 to 0.5)
  REPORT = 0,
  nCR = 3,
  steps = 10 # maybe this is after convergence. Not really
  #DEpairs = 3 # this is not a parameter. You are forced to use the default
)

tic <- as.numeric(Sys.time()) 
set.seed(1)
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
  ) |>
  my_dream(dream_control = dream_control_parameters) 

toc <- as.numeric(Sys.time()) - tic

#plot(example_DREAM)

real_space_example_DREAM <- example_DREAM

real_space_example_DREAM$Sequences <- transform_to_realspace(
  mcmc_result = real_space_example_DREAM$Sequences,
  transform_functions = c(
    linear_parameter_transform, #a0_d - par1
    linear_parameter_transform, #a0_n - par2
    logarithmic_parameter_transform, #a1 - par3
    linear_parameter_transform, #a2 - par4
    linear_parameter_transform, #a3 - par5
    linear_parameter_transform, #a4 - par6
    linear_parameter_transform, #a5 - par7
    logarithmic_parameter_transform, #sd -par8
    logarithmic_parameter_transform # CO2_scale -par9
    ), 
  lower_bounds = real_space_example_DREAM$numerical_optimiser_setup$lower_bound, 
  upper_bounds = real_space_example_DREAM$numerical_optimiser_setup$upper_bound, 
  scale = 100
)

plot(real_space_example_DREAM)
plot(real_space_example_DREAM$Sequences)
ss <- window(real_space_example_DREAM$Sequences, start = end(real_space_example_DREAM$Sequences)/2 + 1) # how does this represent the burn-in period?
plot(ss)














# Playing with Baysian Tools ---------------------------------------------------
# I am not sure how faithful this is to the OG paper...
gauge <- "212209"


example_setup <- gauge |> # this is using my numerical_optimiser_setup object. Not important
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  ) |>
  numerical_optimiser_setup_vary_inputs(
    streamflow_model = streamflow_model_precip_only,
    objective_function = constant_sd_objective_function,
    bounds_and_transform_method = make_default_bounds_and_transform_methods(),
    minimise_likelihood = FALSE
  ) 


example_likelihood_function <- example_setup$ready_objective_function


## Step 1: Create baysianSetup
bayesian_setup <- createBayesianSetup(
  likelihood = example_likelihood_function,
  lower = c(0, 0, 0),
  upper = c(100, 100, 100),
  names = example_setup$parameter_names
)

## Step 2: Specify DREAM settings

## Number of chains is tied to start value matrix (startValue)
## If not provided then it does 2 * number of parameters and randomly samples from prior

set.seed(1)
parameter_number <- 3
chain_number <- 3 * parameter_number
initial_values <- improvedLHS(n = chain_number, k = parameter_number, dup = 1) * 100 #randomLHS(n = chain_number, k = parameter_number) * 100, optimumLHS(n = chain_number, k = parameter_number, maxSweeps = 2, eps = 0.1, verbose = FALSE) * 100

burn_in_per_chain <- 5E3
length_after_burn_in <- 1E3
total_iterations_per_chain <- length_after_burn_in + burn_in_per_chain # this includes burn_in

settings_DREAM_bayesian_tools <- list(
  startValue = initial_values,
  iterations = total_iterations_per_chain * chain_number,
  nCR = 3, # number of cross-over proposals
  gamma = 1, # jump rate - doesn't do a whole lot
  eps = 0, # Ergodicity Term?
  e = 0.05, # Ergodicity Term?
  pCRupdate = TRUE, # Update cross-over probabilities
  updateInterval = 10, # is this the steps in dream?
  burnin = burn_in_per_chain * chain_number, # this counts towards iterations. Not included in chains.
  thin = 1,
  adaptation = (burn_in_per_chain * chain_number) - 1, # number or percentage of samples
  parallel = NULL,
  DEpairs = 3,
  consoleUpdates = 100, # purely for monitoring
  currentChain = 1,
  message = TRUE # purely for monitoring
)

## Step 3: Run
results_DREAM_bayesian_tools <- runMCMC(
  bayesianSetup = bayesian_setup,
  sampler = "DREAM",
  settings = settings_DREAM_bayesian_tools
)



## Step 4: Analyse Results
tracePlot(
  sampler = results_DREAM_bayesian_tools
)

gelmanDiagnostics(
  sampler = results_DREAM_bayesian_tools,
  thin = "auto",
  plot = TRUE
)


#plotDiagnostic(
 # out = results_DREAM_bayesian_tools
#)
