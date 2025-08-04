# Fitting best performing streamflow models to catchments


# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, dream, smoof, furrr, parallel, truncnorm, sloop, arrow)


# Import and prepare data-------------------------------------------------------
start_stop_indexes <- readr::read_csv(
  "./Data/Tidy/start_end_index.csv",
  show_col_types = FALSE
)

data <- readr::read_csv(
  "./Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  show_col_types = FALSE
) |>
  mutate(year = as.integer(year)) |> 
  # required for log-sinh. Log-sinh current formulation has asymptote of zero. 
  # This means zero flows of ephemeral catchments cannot be transformed
  # add a really small value
  mutate(q_mm = q_mm + .Machine$double.eps^0.5) 


gauge_information <- readr::read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE
)


best_CO2_non_CO2_per_gauge_CMAES <- read_csv(
  "Modelling/Results/CMAES/best_CO2_non_CO2_per_catchment_CMAES.csv",
  show_col_types = FALSE
) 


# Functions --------------------------------------------------------------------
source("./Functions/utility.R")
source("./Functions/boxcox_logsinh_transforms.R")
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




# Identify the best streamflow model for a given catchment ---------------------
best_model_and_parameters_per_gauge <- best_CO2_non_CO2_per_gauge_CMAES |> 
  slice_min(AIC, by = gauge) |> 
  arrange(contains_CO2) # try to get CO2 and non-CO2 models in the same chunk for calibration - makes partitioned parquet models more efficient
  

best_model_per_gauge <- best_model_and_parameters_per_gauge |> 
  select(gauge, streamflow_model) |> 
  distinct()



## Convert the streamflow_model column from characters to a list of functions ===
### match.fun can call a function using the character string of the function name
best_model_function_per_gauge <- map(
  .x = best_model_per_gauge |> pull(streamflow_model),
  .f = match.fun
)
  


# Build catchment_data_set objects ---------------------------------------------
catchment_data <- purrr::map(
  .x = best_model_per_gauge |> pull(gauge), # gauges index will be the same as best_model_function_per_gauge
  .f = catchment_data_blueprint,
  observed_data = data,
  start_stop_indexes = start_stop_indexes
)




# Build numerical optimiser for best streamflow model and catchments -----------
DREAM_bounds_and_transform_methods <- make_DREAM_bounds_and_transform_methods( 
  # function factory - use CMAES parameters to narrow down CO2 parameters
  all_gauges_best_CMAES_parameters = best_model_and_parameters_per_gauge
  )


numerical_optimisers <- map2(
  .x = catchment_data,
  .y = best_model_function_per_gauge, 
  .f = numerical_optimiser_setup_vary_inputs,
  objective_function = constant_sd_objective_function,
  streamflow_transform_method = log_sinh_transform,
  bounds_and_transform_method = DREAM_bounds_and_transform_methods, 
  streamflow_transform_method_offset = 0,
  scale = 100,
  minimise_likelihood = FALSE
)



# Chunk numerical optimisers ---------------------------------------------------
# Based off previous DREAM runs 17 gauges per chunks seems to be the limit in term of RAM
GAUGES_PER_CHUNK <- 10L
length_numerical_optimsers <- length(numerical_optimisers)
total_chunks <- ceiling(length_numerical_optimsers / GAUGES_PER_CHUNK)

chunk_index_for_split <- seq(from = 1, to = total_chunks) |> 
  rep(length.out = length(numerical_optimisers)) |> 
  sort() |> 
  as.factor()


chunked_numerical_optimisers <- split(numerical_optimisers, f = chunk_index_for_split)



# Define DREAM controls --------------------------------------------------------
DREAM_controls <- list(
  check_convergence_steps = 1000,
  warm_up_per_chain = 1E5, # low = 1E3
  burn_in_per_chain = 3E4, # low = 1E3
  iterations_after_burn_in_per_chain = 3E4, # low = 1E3
  eps = 0.1, # 1E-6
  steps = 300, 
  thinning = 1
)




# DREAM chunk function ---------------------------------------------------------
# This function:
# 1. runs numerical optimiser chunks in DREAM in parallel
# 2. save the sequences as a parquet file 
# 3. save the convergence stats as a .csv
# 4. saves distribution plots for inspection
# 5. saves trace plots for inspection


chunk_DREAM <- function(single_chunk_numerical_optimiser, chunk_iter, DREAM_controls) {
  
  # Run numerical_optimiser chunk in DREAM in parallel
  chunk_DREAM_results <- future_map(
    .x = single_chunk_numerical_optimiser,
    .f = DREAM,
    controls = DREAM_controls,
    .options = furrr_options(seed = 1L)
  )

  
  # Save sequences 
  map(
    .x = chunk_DREAM_results,
    .f = mcmc_list_to_tibble,
    add_gauge = TRUE
  ) |> 
    list_rbind() |> 
    write_parquet(sink = paste0("./Modelling/Results/DREAM/Sequences/sequence_",chunk_iter, ".parquet"))
  
  
  # Save convergence stats
  map( 
    .x = chunk_DREAM_results, 
    .f = get_convergence_statistics
  ) |> 
    list_rbind() |> 
    write_csv(
      file = paste0("./Modelling/Results/DREAM/convergence_stats_", chunk_iter, ".csv")
    )
  
  
  # Save trace plots
  converged_trace_plots <- map( 
    .x = chunk_DREAM_results,
    .f = gg_trace_plot
  ) 
  
  ggsave(
    filename = paste0("trace_plots_", chunk_iter, ".pdf"),
    plot = gridExtra::marrangeGrob(
      converged_trace_plots, 
      nrow = 1, 
      ncol = 1,
      top = NULL
    ),
    device = "pdf",
    path = "./Figures/Supplementary/",
    width = 297,
    height = 210,
    units = "mm"
  )
  
  
  # Save distribution plots
  converged_distributions_plots <- map(
    .x = chunk_DREAM_results,
    .f = gg_distribution_plot
  )
  
  ggsave(
    filename = paste0("distribution_plots_", chunk_iter, ".pdf"),
    plot = gridExtra::marrangeGrob(
      converged_distributions_plots, 
      nrow = 1, 
      ncol = 1,
      top = NULL
    ),
    # remove page numbers
    device = "pdf",
    path = "./Figures/Supplementary/",
    width = 297,
    height = 210,
    units = "mm"
  )
  
  
  rm(chunk_DREAM_results, converged_trace_plots, converged_distributions_plots)
  gc()
  #browser() 
}



# Run DREAM --------------------------------------------------------------------
plan(multisession, workers = length(availableWorkers()))
iwalk(
  .x = chunked_numerical_optimisers[35:52],
  .f = chunk_DREAM,
  DREAM_controls = DREAM_controls
)




# Combine files ----------------------------------------------------------------
## combine convergence stats
convergence_stats <- list.files(
  path = "./Modelling/Results/DREAM",
  pattern = "convergence",
  full.names = TRUE,
  recursive = FALSE
) |> 
  read_csv(show_col_types = FALSE) |> 
  write_csv("./Modelling/Results/DREAM/DREAM_convergence_stats.csv")


## combine trace plots
trace_plot_paths <- list.files(
  path = "./Figures/Supplementary",
  pattern = "trace",
  full.names = TRUE,
  recursive = FALSE
)

qpdf::pdf_combine( 
  input = trace_plot_paths, 
  output = "./Figures/Supplementary/DREAM_trace_plots.pdf"
)

## combine distribution plots
distribution_plot_paths <- list.files(
  path = "./Figures/Supplementary",
  pattern = "distribution",
  full.names = TRUE,
  recursive = FALSE
)

qpdf::pdf_combine( 
  input = distribution_plot_paths, 
  output = "./Figures/Supplementary/DREAM_distribution_plots.pdf"
)







# Single catchment test
#test_DREAM <- DREAM(
#  input = numerical_optimisers[[22]],
#  controls = DREAM_controls
#)


#get_convergence_statistics(test_DREAM)
#gg_trace_plot(test_DREAM)
#gg_distribution_plot(test_DREAM)
#save_sequences(test_DREAM, sink = "./Modelling/Results/DREAM/test.parquet")



