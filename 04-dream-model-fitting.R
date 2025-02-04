## Fitting streamflow models to catchments
cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, dream, coda, lattice, tictoc, furrr, parallel, truncnorm, sloop)
# install dream using: install.packages("dream", repos="http://R-Forge.R-project.org")
# install.packages("coda")

# TODO -> run it with all catchments

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

CMAES_results <- read_csv("./Results/my_cmaes/CMAES_parameter_results_20250122.csv",
  show_col_types = FALSE
) 

best_CO2_non_CO2_per_gauge <- read_csv(
  "./Results/my_cmaes/unmodified_best_CO2_non_CO2_per_catchment_CMAES_20250128.csv",
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
source("./Functions/DREAM.R")
source("./Functions/objective_function_setup.R")
source("./Functions/result_set.R")



# Identify best model combination for DREAM gauges -----------------------------
best_model_combination_per_catchment <- best_CO2_non_CO2_per_gauge |>
  select(gauge, streamflow_model, AIC) |> 
  slice_min(
    AIC,
    by = gauge
  ) |> 
  distinct() |> 
  select(!AIC) |> 
  rename(
    streamflow_model_name = streamflow_model
  ) |> 
  add_column(
    objective_function_name = constant_sd_objective_function()$name,
    .after = 2
    )




## Use a join to assign each row the correct streamflow and objective function =====
get_model_name_and_function <- function(model) {
  model_name <- model()$name
  tibble(
    "model_name" = model_name,
    "model" = list(model)
  )
}


streamflow_name_and_model_for_joining <- map(
  .x = c(get_non_drought_streamflow_models(), get_drought_streamflow_models()),
  .f = get_model_name_and_function
) |>
  list_rbind() |> 
  rename(
    streamflow_model_name = model_name,
    streamflow_model = model
  )


objective_function_name_and_model_for_join <- map(
  .x = get_all_objective_functions(),
  .f = get_model_name_and_function
) |> 
  list_rbind() |> 
  rename(
    objective_function_name = model_name,
    objective_function = model
  )



## Left join ===================================================================
gauge_streamflow_model_objective_function_combinations <- best_model_combination_per_catchment |> 
  left_join(
    streamflow_name_and_model_for_joining,
    by = join_by(streamflow_model_name)
  ) |> 
  left_join(
    objective_function_name_and_model_for_join,
    by = join_by(objective_function_name)
  ) |> 
  select(!c(streamflow_model_name, objective_function_name)) |> 
  unclass() |> 
  unname() # names breaks pmap? I don't know why - remove them






# Produce numerical_optimiser_objects ready for the optimiser ------------------
numerical_optimiser_objects <- function(gauge, streamflow_model, objective_function, data, start_stop) {
  
  catchment_data_set <- gauge |>
    catchment_data_blueprint(
      observed_data = data,
      start_stop_indexes = start_stop
    ) 
  
  catchment_data_set |> 
    numerical_optimiser_setup_vary_inputs(
      streamflow_model = streamflow_model,
      objective_function = objective_function,
      bounds_and_transform_method = make_default_bounds_and_transform_methods(catchment_data_set), # this must adjust for CO2 from catchment_data_set
      minimise_likelihood = FALSE
    )
}



ready_for_optimisation <- pmap(
  .l = gauge_streamflow_model_objective_function_combinations,
  .f = numerical_optimiser_objects,
  data = data,
  start_stop = start_stop_indexes
)


# Split ready_for_optimisation objects into chunks to avoid exceeding RAM ------
# Subject to change. I need to work out RAM requirements

# 16 * number_of_catchments_in_chunk * RAM_usage_of_single_catchment < 32 * 0.7 
# Assume RAM_usage_of_single_catchment is 0.5 Gb
# Solve or number_of_catchments_in_chunk:
# number_of_catchments_in_chunk < (32 * 0.7) / (16 * 0.5) ~ 2.8
# I think this is worst than running the chunks in series but catchments in
# parallel

CHUNK_SIZE <- 32

chunked_ready_for_optimisation <- split( # required for chunking in parallel
  ready_for_optimisation,
  ceiling(seq_along(ready_for_optimisation) / CHUNK_SIZE)
)




# Findings (text-only) ---------------------------------------------------------
# Everything is working as expected for gauge 003303A and 105101Aa

# Plan for scaling up:
# 1. Run DREAM with convergence for a given chunk.
# 2. Save dream_objects, trace plots (plots in list and gridExtra::marrangeGrob) and convergences stats (save .csv)
# 3. If good DREAM without convergence. Save sequences and distribution plots. See
#    if I can repurpose run_and_save_chunks_optimiser_parallel to save results

# Size information:
# - The convergence dream object is 1.5 Gb. This limits my chunk size to
#   32 / 1.5 = ~ 20 per chunk
# during convergence with 20 per chunk RAM maxed out at 5gb. 
# this suggests I could increase it?




# Learning from trace plots (all text) -----------------------------------------

## First (mk_1)
## 1. Single catchment did not converge (121003A). Failed in burn-in (burn.in flag is TRUE)
## 2. I think there are insufficient chains for some catchments. Not enough wiggle
##    Up default chain number from (2 * parameter_number) + 1 to
##    3 * parameter_number


## Second (mk_2)
## 1. Still not enough wiggle for some catchments. For example, 107001B is 
##    not exploring the sample space efficiently. However, with enough iterations
##    the problem will go away. I would like to improve the efficiency (more up/down).
##    Ways to improve efficiency are more chains, gamma?, nCR?, eps?, DEpairs (DEpairs default is maximum allowed. I can make it smaller but I think this would harm)
##    See paper for question marks


## Third (mk_3)
## 1. Looked at the DREAM (vrugt) and DREAM (r code).
## 2. steps represents the number of generations. In DREAM manual its > 1000
##    increase to 100 for testing.
## 3. nCR alters pCR  (probability of cross-over). As nCR increases pCR decreases.
##    Small pCR means smaller steps. Stick to 3.
## 4. gamma does not do anything when using logposterior.density. So we ignore it.
## 5. There are two eps values in "offde" function. This reminds me of 
##    kernel_adapt from fmcmc. fmcmc documentation: eps both sets the initial 
##    scale for the multivariate normal kernel 
##    (later replaced with actual variance-covariance of model after warm-up)
##    and ensures the variance-covariance is greater than zero. 
##    The fixed eps (not user input) looks like it ensures variance-cov 
##    is greater than zero. From testing a larger eps equals larger amplitude of
##    wiggles. Change eps from 0.05 to 0.1

## Lessons from mk 3
## 1. increasing steps improves mixing.
## 2. Bad catchments are: 
##    - 121003A (did not converge, poor acceptance rate, poor mixing, not stationary)
##    - 122004A (poor mixing - improved either with steps, chains or nCR)
## 121003A performs poorly work out why.
## 3. increasing steps does not increase iterations

## Maybe try running for all catchments - identify bad catchments, then see
## what works.
## a3 is not good if a5 > max_CO2. Maybe limit a5 to max CO2 to 
## force a3 to be meaningful? Get close enough to max. Since each 
## gauge ends at a different time it will vary.
## Based off distribution results and trace diagrams a second pass into 
## dream may not be required. If the traces are good, leave it.

## Forth (mk_4)
## 1. I want to get 121003A to converge. Traces are still not stationary.
## 2. Method to try and get 121003A to converge:
##    - increase eps, increase nCR, adjust bounds, increase iterations, check if model is correct
## 3. I am not concerned with 122004A. Cranking steps should fix the issue.

### mk_4 rapid fire testing:
## Default (check_convergence_steps = 1000, warm_up_per_chain = 5E4, burn_in_per_chain = 5E3, iterations_after_burn_in_per_chain = 5E3, nCR = 3, eps = 0.1, steps = 10,
## Test 1: increase nCR from 3 to 6 --> result: no noticeable impact
## Test 2: increase eps from 0.1 to 0.5 --> result: worse for a0_d and a5
## Test 3: decrease eps from 0.1 to 0.02 --> result: no major improvement
## Test 4: lower a5 bound to max CO2 in t-series (136.41) --> result: extreme improvement. With steps can get hairy caterpillar 
## Test 5: can I get away with lowering a5 bound to naive max (last CO2 not for the t-series 138.53) --> result: just as good. But likely to break if t-series ends well before 2022
## Test 6: full run with a5 = 136.41 (check_convergence_steps = 1000, warm_up_per_chain = 5E4, burn_in_per_chain = 1E4, iterations_after_burn_in_per_chain = 1E4, nCR = 3, eps = 0.1, steps = 100)
##         test 6 results -> not great more wiggle required, low acceptance rate and did not converge
## Test 7: increase nCR to 8 (param #), steps = 10, way better.
## Test 8: increase nCR to 16 - no meaningful change
## Test 9: nCR = 8, limit bounds, eps = c(0.02, 0.1, 0.5) does it matter? Yes, 0.5 too much and 0.02 is too little. I think stick to 0.1. Both 0.02 and 0.1 converged.
## Test 10: Vary chains. 4 * param number. Does it matter? Yes, more chains = better. Did not converge.
## Test 11: Full gas. nCR = param, limit bounds, eps = 0.1, steps = 500, double iterations for everything
##          results - Everything is good except for the CO2 models. I think
##          I should redo the CMAES with the adjusted limits or at least test them.
##          If the a5 parameter is near the end I am sceptical.
##          Used around ~ 6Gbs



## Lessons from mk_4
## 1. a5 bound must be adjusted to the maximum CO2 for the given catchment.
##    This means the a5 upper bound must change for each catchment. CHANGE MADE
## 


# Dream manual research (text only) ---------------------------------------
# Direct parameters:
# nCR = 3 works well in practice but larger nCR maybe required for high-dimensional problems (d > 50)
# nCR = alters pCR (cross-over or selection probability). Larger nCR smaller pCR (1/nCR). Smaller jumps?
# eps = 0.1 (is it c? pg. 20) random error for ergodicity (is this randomisation or ergodicity or combination of the two)
# gamma is used in CalcCbWb(gamma) equation Thiemann et al. WRR 2001 equation 20. Does not do anything if func.type = "log-posterior"
# steps is the number of times through loops. This looks the same as DREAM. I think it represents number of samples from the prior distribution. 
# gamma = kurtosis parameter baysian (is this the shaping factor in manual)
# "rand" bounds handling is good for me - fold is preferred statistically (dream manual) but if posterior is on edege of search domain problems

# Functions that take control parameters as an input (hidden to user): 
# GenCR, AdaptCR, metrop, CompDensity (gamma, Wb, Cb) and offde 
# - GenCR takes nCR
# - in CompDensity gamma has no impact on logposterior.density
# - Adapt CR similar to genCR. Only nCR
# - metrop gamma does not do anything if func.type = "logposterior.density"
# - offde = eps = 1e-6 * randn(nseq, ndim). eps is being used for noise.x
#         There are 2 eps values. This reminds me of kernel_adapt from fmcmc
#         In fmcmc eps both sets the initital scale for the multivariate normal kernel (replaced with actual variance-covariance of model after warmup)
#         and ensures the variance-covariance is greater than zero. The fixed eps looks like it ensures variance-cov is greater than zero.
# DEstrategy

# Lowering steps significantly speeds up processes. If it is too low it produces an
# error Error in AdaptpCR(CR, delta.tot, lCR, control) : AdaptpCR: no changes in X, would cause NaN in pCR
# Something to do with the 'evolution of chains'

# There is no jump-rate parameter in dream R. It is done for you.
# The larger nCR is the larger steps needs to be

# Summary:
# gamma has no impact because 
# The parameters are dependent on the posterior. Problem is I don't know the 
# posterior. Requires trial-and-error. Hopefully, I can use a single control
# parameter set. I think I can
# Small steps = speeeeeed.
# Large steps = really good mixing = slow. Scales linearly i.e., double steps = double time
# The more steps the better. For [[11]] max steps out, apply for 
# every other catchment. This will likely take a while
# Steps of 200 worked really well for catchment [[1]] 
# Increasing steps will be really good for converged catchments
# For non-converged either increase iterations or 
# lower acceptance rate means (< 15 %) indicates the posterior surface
# is difficult to traverse in pursuit of the target distribution. Can reduce 
# jump rate to increase acceptance rate (not in R). Alternatively, lower nCR
# For non-converged either increase iterations (burn-in mainly), nCR or chains
# A low acceptance rate means (< 15 %) indicates the posterior surface
# is difficult to traverse in pursuit of the target distribution. In other languages 
# you would reduce jump rate. Jump rate is fixed in DREAM R. Instead we could
# reduce eps or increase nCR. If all else fails more chains and iterations.


# Run DREAM --------------------------------------------------------------------

## Run dream on chunks and save results ========================================
# Ideally this would use run_and_save_chunks_optimiser_parallel but
# for now this is good enough
optimise_chunks <- function(chunk_ready_for_optimisation, chunk_iter, chunk_controls) {
  
  tic()
  
  ## Run dream for each catchment in chunk =====================================
  converged_dream_objects <- future_map(
    .x = chunk_ready_for_optimisation, # this must be repeated for all chunks
    .f = DREAM,
    controls = chunk_controls,
    .options = furrr_options(
      seed = TRUE
      # Future cannot resolve local S3 methods 
      # I have added work around functions in generic functions. These are ugly
    )
  )
  
  ## Create and save trace diagrams ============================================
  converged_trace_plots <- map( 
    .x = converged_dream_objects,
    .f = gg_trace_plot
  ) 
  
  ggsave(
    filename = paste0("trace_plots_chunk_", chunk_iter, "_", get_date(), ".pdf"),
    plot = gridExtra::marrangeGrob(
      converged_trace_plots, 
      nrow = 1, 
      ncol = 1,
      top = NULL
      ),
    device = "pdf",
    path = "./Graphs/DREAM_graphs",
    width = 297,
    height = 210,
    units = "mm"
  )
  
  ## Create and save convergence statistics ====================================
  converged_stats <- map( 
    .x = converged_dream_objects, 
    .f = get_convergence_statistics
  ) |> 
    list_rbind() |> 
    write_csv(
      file = paste0("./Results/my_dream/converged_stats_chunk_", chunk_iter, "_", get_date(), ".csv")
    )
  
  ## Create and save distribution plots ========================================
  converged_distributions_plots <- map(
    .x = converged_dream_objects,
    .f = gg_distribution_plot
  )
  
  ggsave(
    filename = paste0("distribution_plots_chunk_", chunk_iter, "_", get_date(), ".pdf"),
    plot = gridExtra::marrangeGrob(
      converged_distributions_plots, 
      nrow = 1, 
      ncol = 1,
      top = NULL
      ),
    # remove page numbers
    device = "pdf",
    path = "./Graphs/DREAM_graphs",
    width = 297,
    height = 210,
    units = "mm"
  )
  
  cat(paste0("Chunk ", chunk_iter, " complete "))
  toc()
  cat("\n")
}


## DREAM controls ==============================================================
chunk_controls <- list(
  check_convergence_steps = 1000,
  warm_up_per_chain = 1E5,
  burn_in_per_chain = 3E4, 
  iterations_after_burn_in_per_chain = 3E4, 
  eps = 0.1,
  steps = 500
)

## Run it all ==================================================================
plan(multisession, workers = length(availableWorkers())) # set once for furrr
iwalk(
  .x = chunked_ready_for_optimisation,
  .f = optimise_chunks,
  chunk_controls = chunk_controls
)


## Combine chunks into single file (plots and stats) and removes chunks ========

### .csv's #####################################################################

converge_stats_list_of_files <- list.files( # get
  path = "./Results/my_dream/",
  recursive = FALSE, # I don't want it looking in other folders
  pattern = "converged_stats_chunk",
  full.names = TRUE
)

converge_stats_list_of_files |> # merge and save
  read_csv(show_col_types = FALSE) |> 
  write_csv(
    paste0("./Results/my_dream/converged_stats_", get_date(), ".csv")
  )

converge_stats_list_of_files |> file.remove() # remove chunks


### .pdfs ######################################################################
trace_plots_list_of_files <- list.files( # get
  path = "./Graphs/DREAM_graphs/",
  recursive = FALSE, # I don't want it looking in other folders
  pattern = "trace_plots_chunk",
  full.names = TRUE
)


qpdf::pdf_combine( # combine
  input = trace_plots_list_of_files, 
  output = paste0("./Graphs/DREAM_graphs/trace_plots_", get_date(), ".pdf")
)


distribution_plots_list_of_files <- list.files( # get 
  path = "./Graphs/DREAM_graphs/", 
  recursive = FALSE, # I don't want it looking in other folders
  pattern = "distribution_plots_chunk",
  full.names = TRUE
)

qpdf::pdf_combine( # combine
  input = distribution_plots_list_of_files, 
  output = paste0("./Graphs/DREAM_graphs/distribution_plots_", get_date(), ".pdf")
)

c(trace_plots_list_of_files, distribution_plots_list_of_files) |> file.remove() # remove



