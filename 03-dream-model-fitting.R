## Fitting streamflow models to catchments
cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, dream, coda, lattice, tictoc, furrr, parallel, truncnorm, sloop)
# install dream using: install.packages("dream", repos="http://R-Forge.R-project.org")
# install.packages("coda")

# floor(RAM / 3Gb) = CHUNK_SIZE


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

CMAES_results <- read_csv("./Results/my_cmaes/CMAES_parameter_results_20241130.csv",
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
source("./Functions/DREAM.R")
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
  gauge |>
    catchment_data_blueprint(
      observed_data = data,
      start_stop_indexes = start_stop
    ) |>
    numerical_optimiser_setup_vary_inputs(
      streamflow_model = streamflow_model,
      objective_function = objective_function,
      bounds_and_transform_method = make_default_bounds_and_transform_methods(),
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
CHUNK_SIZE <- 16 # 28 double change size?

chunked_ready_for_optimisation <- split( # required for chunking in parallel
  ready_for_optimisation,
  ceiling(seq_along(ready_for_optimisation) / CHUNK_SIZE)
)



# Everything is working as expected for gauge 003303A and 105101A

# Plan for scaling up:
# 1. Run DREAM with convergence for a given chunk.
# 2. Save dream_objects, trace plots (plots in list and gridExtra::marrangeGrob) and convergences stats (save .csv)
# 3. Inspect if the trace plots and stats are good
# 4. If good DREAM without convergence. Save sequences and distribution plots. See
#    if I can repurpose run_and_save_chunks_optimiser_parallel to save results

# Size information:
# - The convergence dream object is 1.5 Gb. This limits my chunk size to
#   32 / 1.5 = ~ 20 per chunk
# during convergence with 20 per chunk RAM maxed out at 5gb. 
# this suggests I could increase it?

# I think it will be start with a chunk. 
# Run to converge then run to make distribution angle.
# RAM limitations.
# Manually rince and repeat
# Maximum chunk size is around 20. This means I must repeat the process
# 27 times :(


# Learning from trace plots:
## 1. Single catchment did not converge (121003A)
## 2. I think there are insufficient chains for some catchments. Not enough wiggle
##    Up default chain number from (2 * parameter_number) + 1 to
##    3 * parameter_number

## Second (mk_2)
## 1. Still not enough wiggle for some catchments. For example, 107001B is 
##    not exploring the sample space efficiently. However with enough iterations
##    the problem will go away. I would like to improve the efficiency (more up/down).
##    Ways to improve efficiency are more chains, gamma?, nCR?, eps?, DEpairs (DEpairs default is maximum allowed. I can make it smaller but I think this would harm)
##    See paper for question marks


## Third (mk_3)
## 1. Looked at the DREAM (vrugt) and DREAM (r code).
## 2. steps represents the number of generations. In DREAM manual its > 1000
##    increase to 200 for testing.
## 3. nCR alters pCR  (probability of cross-over). As nCR increases pCR decreases.
##    Small pCR means smaller steps.
## 4. gamma does not do anything when using logposterior.density. So we ignore it.
## 5. There are two eps values in "offde" function. This reminds me of 
##    kernel_adapt from fmcmc. fmcmc documentation: eps both sets the initial 
##    scale for the multivariate normal kernel 
##    (later replaced with actual variance-covariance of model after warm-up)
##    and ensures the variance-covariance is greater than zero. 
##    The fixed eps (not user input) looks like it ensures variance-cov 
##    is greater than zero. From testing a larger eps equals larger amplitude of
##    wiggles.

# Is sem number of generations (T in matlab)?

## 2. Gauge 121003A still has not converged. Ran out of iterations. 
##    Failed in burn-in (burn.in flag is TRUE)

# From dream manual
# nCR = 3 works well in practice but larger nCR maybe required for high-dimensional problems (d > 50)
# nCR = alters pCR (cross-over or selection probability). Larger nCR smaller pCR (1/nCR). Smaller jumps?
# eps = 0.1 (is it c? pg. 20) random error for ergodicity (is this randomisation or ergodicity or combination of the two)
# I don't think eps is connected to anything
# gamma is used in CalcCbWb(gamma) equation Thiemann et al. WRR 2001 equation 20. I don't think it is connected to anything? See GenCR function
# steps seems like T value from matlab (number of generations)
# maybe a side-by-side of code is required to work out what parameters do what.
# steps is the number of times through loops. This looks the same as DREAM. How did vrugt pick T in manual?
# gamma = kurtosis parameter baysian (is this the shaping factor in manual)

# functions that take control as an input
# GenCR (only nCR), AdaptCR (), metrop, CompDensity (gamma, Wb, Cb), offde, 
# in CompDensity gamma has no impact on logposterior.density
# Adapt CR similar to genCR. Only nCR
# metrop gamma has no impact on logposterior.density
# offde = eps = 1e-6 * randn(nseq, ndim). eps is being used for noise.x
#         There are 2 eps values. This reminds me of kernel_adapt from fmcmc
#         In fmcmc eps both sets the initital scale for the multivariate normal kernel (replaced with actual variance-covariance of model after warmup)
#         and ensures the variance-covariance is greater than zero. The fixed eps looks like it ensures variance-cov is greater than zero.
# DEstrategy

# Lowering steps significantly speeds up processes. If it is too low it produces an
# error Error in AdaptpCR(CR, delta.tot, lCR, control) : AdaptpCR: no changes in X, would cause NaN in pCR
# Something to do with the 'evolution of chains'

# There is no jump-rate parameter in dream R. It is done for you.
# The larger nCR is the larger steps needs to be

# From research:
# gamma has no impact because func.type = "logposterior.density"
# eps does not matter as it gets replaced with iterations after warm-up. Pick a sensible value i.e., default
# a larger eps does seem to matter for convergence. Larger value better spread around chain
# larger nCR encourages smaller jumps. If too small it cannot improve getting stuck = error. Can be fixed by increasing steps.
# steps is the number of generations (in the examples T = 5000)
# large nCR, small steps requires large eps
# The parameters are dependent on the posterior. Problem is I don't know the 
# posterior. Requires trial-and-error
# Get steps small = speeeeeed.
# steps large = really good mixing = slow. Scales linearly i.e., double steps = double time
# The more steps the better. For [[11]] max steps out, apply for 
# every other catchment. This will likely take a while
# Steps of 200 worked really well for catchment [[1]] 
# Increasing steps will be really good for converged catchments
# For non-converged either increase iterations or 

## An efficient way to see if the changes have worked is test the single catchments
#test_catchment <- chunked_ready_for_optimisation[[1]][[11]] # 11 = non-converging (700 sec steps = 200, max iter limit)

# does increasing steps increase iterations? Try using [[1]] 
# increasing steps does not increase iterations

#tic()
#set.seed(1)
#test_converged_dream <- DREAM(
#  input = test_catchment,
#  controls = list(
#    check_convergence_steps = 1000,
#    #warm_up_per_chain = 1E5,
#    #burn_in_per_chain = 2E4, 
#    #iterations_after_burn_in_per_chain = 2E4, 
#    nCR = 3,
#    eps = 0.1,
#    steps = 7,
#    bound_handling = "rand" # rand is good for me - fold is preferred statistically (dream manual) but if posterior is on edege of search domain problems
#    )
#)
#toc()

# nCR = 3, eps = 0.1 and steps = 7 = convergences for non-converging catchment, chain = 3 * param number 
#gg_trace_plot(test_converged_dream)
#get_convergence_statistics(test_converged_dream)

# lower acceptance rate means (< 15 %) indicates the posterior surface
# is difficult to traverse in pursuit of the target distribution. Can reduce 
# jump rate to increase acceptance rate (not in R). Alternatively, lower nCR
# 

# Convergence check ------------------------------------------------------------
CHUNK_ITER <- 1 # length(chunked_ready_for_optimisation)

plan(multisession, workers = length(availableWorkers())) # set once for furrr

converged_dream_objects <- future_map(
  .x = chunked_ready_for_optimisation[[CHUNK_ITER]], # this must be repeated for all chunks
  .f = DREAM,
  controls = list(
    check_convergence_steps = 1000,
    warm_up_per_chain = 1E5,
    burn_in_per_chain = 2E4, 
    iterations_after_burn_in_per_chain = 2E4, 
    nCR = 3,
    eps = 0.1,
    steps = 200,
    ),
  .options = furrr_options(
    seed = TRUE
    # Future cannot resolve local S3 methods 
    # I have added work around functions in generic functions. These are ugly
    )
  )


converged_trace_plots <- map( 
  .x = converged_dream_objects,
  .f = gg_trace_plot
) 

ggsave(
  filename = paste0("trace_plots_chunk_", CHUNK_ITER, "_", get_date(), ".pdf"),
  plot = gridExtra::marrangeGrob(converged_trace_plots, nrow = 1, ncol = 1),
  device = "pdf",
  path = "./Graphs/DREAM_graphs",
  width = 297,
  height = 210,
  units = "mm"
)


converged_stats <- map( 
  .x = converged_dream_objects, 
  .f = get_convergence_statistics
) |> 
  list_rbind() |> 
  write_csv(
    file = paste0("./Results/my_dream/converged_stats_chunk_", CHUNK_ITER, "_", get_date(), ".csv")
  )


stop_here <- 1


# put back into DREAM and run without convergences to build the distribution?
test_distribution_DREAM <- tactical_typo_DREAM(
  input = test_convergence_DREAM,
  controls = list(
    check_convergence_steps = 0, # don't check for convergence,
    thinning = 10
  )
)

distribution_plot <- gg_distribution_plot(test_distribution_DREAM)
distribution_plot



stop_here <- 1
# replace current optimisation strategy


# Run the optimiser (calibration) ----------------------------------------------

## Non-drought optimising ======================================================

# tactical typo to stop runaway code when run all is used
tactical_typo_plan(multisession, workers = length(availableWorkers())) # set once for furrr
iwalk(
  .x = chunked_ready_for_optimisation, 
  .f = run_and_save_chunks_optimiser_parallel,
  optimiser = my_dream,
  save_streamflow = FALSE,
  save_sequences = TRUE,
  is_drought = FALSE
)




# Join everything into a single file and save ----------------------------------
## Parameters ==================================================================
parameters_list_of_files <- list.files(
  path = "./Results/my_dream/",
  recursive = FALSE, # I don't want it looking in other folders
  pattern = "parameter",
  full.names = TRUE
)



parameters_list_of_files |>
  readr::read_csv(show_col_types = FALSE) |>
  dplyr::slice_min(
    loglikelihood,
    by = c(gauge, streamflow_model, objective_function) # only get the minimum LL for gauge, streamflow model and objective function combination
  ) |>
  readr::write_csv(
    file = paste0("./Results/my_dream/DREAM_parameter_results.csv")
  )


## Sequences ===================================================================
sequences_list_of_files <- list.files(
  path = "./Results/my_dream/",
  recursive = FALSE, # I don't want it looking in other folders
  pattern = "sequences",
  full.names = TRUE
)


sequences_list_of_files |>
  readr::read_csv(show_col_types = FALSE) |>
  readr::write_csv(
    file = paste0("./Results/my_dream/DREAM_sequence_results.csv")
  )


## Delete all files ending with .csv ===========================================
 delete_files <- list.files(
  path = "./Results/my_dream/",
  recursive = FALSE, # I don't want it looking in other folders
  pattern = "chunk", # remove all chunk files
  full.names = TRUE
 ) |>
  file.remove()










