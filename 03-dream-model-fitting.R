## Fitting streamflow models to catchments
cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, dream, coda, tictoc, furrr, parallel, truncnorm, sloop, tictoc)
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


# Psuedo-code for dream 
# Options
## 1. repeat 01-cmaes approach
## 2. Repeat the run dream pipe thing 

# pipe thing
# repeat_pipe <- function(gauge, streamflow_model, objective_function, data, start_stop) {
#}
# Need to separate into drought and non-drought
# separate again into batches based on ram requirements
# I could package everything (gauge, streamflow_model and objective function) into
# a single list then iterate over the list rather than using pmap?
# This would make batching/chunking easier.


# Steps.
# Separate drought catchments vs. all catchments. Drought catchment only get drought streamflow models
# Combine all gauge, streamflow model and objective function combinations in a list for iteration
# Split list up based on batches. Batches determined based on RAM requirements
# Run the repeat_pipe
# Save the results similar to 01-cmaes i.e., do not assign a variable,
# write_csv, append results and gc(). 
## Saving results - I would like the optimal parameters and sequences. Ignore the streamflow-time.
## sequences need to have |gauge|model|obj|parameter|parameter_value and rbind()



all_gauges <- unique(gauge_information$gauge)[1:2]

drought_gauges <- gauge_information |> 
  filter(drought) |> 
  pull(gauge) |> 
  unique()


non_drought_streamflow_models <- get_non_drought_streamflow_models()[1:2]

drought_streamflow_models <- get_drought_streamflow_models()[1:2]

all_objective_functions <- get_all_objective_functions()


# I need to find all the combinations of all_gauges, non_drought_streamflow_models and all_objective_functions
non_drought_combinations <- tidyr::expand_grid(
  all_gauges,
  non_drought_streamflow_models,
  all_objective_functions
) |> 
  unclass() |> 
  unname() # names breaks pmap? I don't know why

# remove a list level for the streamflow_model and objective_functions
#x <- list_flatten(non_drought_combinations$non_drought_streamflow_models)


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




ready_for_optimisation <- pmap(.l = non_drought_combinations,
          .f = numerical_optimiser_objects,
          data = data,
          start_stop = start_stop_indexes
          )








# Allows running in parallel with chunking to not exceed RAM
run_and_save_chunks_optimiser_parallel <- function(chunked_numerical_optimisers, chunk_id, optimiser, save_streamflow, save_sequences, is_drought) {
  
  tictoc::tic()
  
  calibrated_results <- furrr::future_map(
    .x = chunked_numerical_optimisers,
    .f = optimiser,
    print_monitor = FALSE,
    .options = furrr_options(
      globals = TRUE,
      seed = TRUE
    )
  )
  
  
  # Purrr does not work in parallel, so I don't need plan(sequential)
  sort_results <- purrr::map(.x = calibrated_results, .f = result_set)
  
  optimiser_name <- as.character(substitute(optimiser))
  
  # I do not want to assign a variable name. Garbage collector works better like this.
  purrr::map(.x = sort_results, .f = parameters_summary) |>
    purrr::list_rbind() |>
    readr::write_csv(
      file = paste0(
        "./Results/",
        optimiser_name,
        "/",
        if_else(is_drought, "drought_", ""),
        "parameter_results_chunk_",
        chunk_id,
        "_",
        get_date(),
        ".csv"
      )
    )
  
  if (save_streamflow) {
    purrr::map(
      .x = sort_results,
      .f = modelled_streamflow_summary
    ) |>
      purrr::list_rbind() |>
      readr::write_csv(
        file = paste0(
          "./Results/",
          optimiser_name,
          "/",
          if_else(is_drought, "drought_", ""),
          "streamflow_results_chunk_",
          chunk_id,
          "_",
          get_date(),
          ".csv"
        )
      )
  }
  
  if (save_sequences) {
    purrr::map(
      .x = sort_results,
      .f = sequences_summary
    ) |>
      purrr::list_rbind() |>
      readr::write_csv(
        file = paste0(
          "./Results/",
          optimiser_name,
          "/",
          if_else(is_drought, "drought_", ""),
          "sequences_results_chunk_",
          chunk_id,
          "_",
          get_date(),
          ".csv"
        )
      )
  }
  
  
  
  cat(paste0("Chunk ", chunk_id, " complete "))
  tictoc::toc()
  cat("\n")
  
  
  # Remove objects for garbage collection
  rm(list = c("calibrated_results", "sort_results", "optimiser_name"))
  
  # Call garbage collection
  gc()
}


x <- run_and_save_chunks_optimiser_parallel(
  chunked_numerical_optimisers = ready_for_optimisation,
  chunk_id = 1, 
  optimiser = my_dream,
  save_streamflow = TRUE,
  save_sequences = TRUE,
  is_drought = FALSE
  )




# Run DREAM --------------------------------------------------------------------
tic()
gauge <- "407214"

dream_example <- gauge |> 
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  ) |> 
  numerical_optimiser_setup_vary_inputs(
    streamflow_model = streamflow_model_precip_only,
    objective_function = CO2_variable_objective_function,
    bounds_and_transform_method = make_default_bounds_and_transform_methods(),
    minimise_likelihood = FALSE
  ) |> 
  my_dream() |> 
  result_set() 

dream_parameters <- dream_example |> 
  parameters_summary()

dream_example |> 
  plot()

toc()




test_sequences |> 
  pivot_longer(
    cols = everything(),
    names_to = "parameter",
    values_to = "parameter_values"
  ) |> 
  ggplot(aes(x = parameter_values)) +
  geom_histogram(
    binwidth = binwidth_bins(30),
    fill = "grey",
    colour = "black"
    ) +
  theme_bw() +
  scale_y_sqrt() +
  labs(
    x = "Range of Parameter Values",
    y = "Frequency"
  ) +
  facet_wrap(~parameter, scales = "free")





# for gauge = 407214, model = streamflow_model_drought_separate_CO2_seasonal_ratio_auto
# and objective function CO2_variable_objective_function expect LL of 95
# 550,000 nseq and 425 sec (~7min) gives 105 LL
# I could bump it more. Make a stepper curve
# New equation used -> round_any(((x - 2) ^ 2.5) * 1E4, 1E4, ceiling)
# Untested for 8 parameter models. 

# N_PARAMETERS = 3 -> nseq = 50,000 
# Max parameters are 8. Extrapolate? if linear it would suggest round_any((50000/3), 1E4, ceiling) = 20000 per parameter
# I don't this this will work. I think the relationship between nseq and number of parameters is non-linear
# Want to get the likelihoods to the nearest whole number of DREAM and CMAES

# For gauge 407214 maximum parameters DREAM LL = 99.xxx, CMAES LL = 95.xxx
# need to increase from 20000 * PARAMETER_NUMBER to or make it a non-linear relationship 160,000 for 8 params its too little. Try 200,000? Not enough.
# Try 300,000? Not enough.
# Playing with non-linear relationship for nseq and n_parameters:
x <- seq(from = 3, to = 8, by = 1)
y <- ((x ^ 2) - (x + 1)) * 1E4
plot(x, y)

yy <- round_any(((x - 2) ^ 2.5) * 1E4, 1E4, ceiling)
yy
# non-linear relationship works with gauge 407214 and 3, 4 parameters
# test for 5, 6, 7, 8? Go straight to 8? Good idea to test all.


# Code for dealing with negative mean box-cox streamflow values ----------------
# matrix'ify this after DREAM
# Trials for handling repression predictions with Q<0 when
# using a truncated normal likelihood

# Let the lower bound be Q=0 and upper bound Q=inf
a <- 0
b <- Inf

# Define a truncated normal st. dev
sigma <- 2

# Let's predict streamflow using annual P
P <- 200
a0 <- -12
a1 <- 0.05
mu <- a0 + a1 * P

# Now lets adjust the predicted mean flow for P to be >0
# Specifically, lets treal Q_reg as the mean of the non-truncated distribution of
# Q for the given P
alpha <- (0 - mu) / sigma
beta <- (Inf - mu) / sigma
mu_true <- mu + sigma * (dnorm(alpha) - dnorm(beta)) / (pnorm(beta) - pnorm(alpha))

sigma_true <- sqrt(sigma^2 * (1 - (-alpha * dnorm(alpha)) / (pnorm(beta) - pnorm(alpha)) -
                                ((dnorm(alpha) - dnorm(beta)) / (pnorm(beta) - pnorm(alpha)))^2))
