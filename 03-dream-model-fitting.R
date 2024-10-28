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





# Identify gauges --------------------------------------------------------------

all_gauges <- unique(gauge_information$gauge)#[1:2] # remove when done

drought_gauges <- gauge_information |> 
  filter(drought) |> 
  pull(gauge) |> 
  unique()

#drought_gauges <- drought_gauges[1:2] # remove when done



# Indentify streamflow_models and objective functions --------------------------
non_drought_streamflow_models <- get_non_drought_streamflow_models()#[1:2] # remove when done

drought_streamflow_models <- get_drought_streamflow_models()#[1:2] # remove when done

all_objective_functions <- get_all_objective_functions()



# Find all the combinations of gauges, streamflow model and objective funs -----
non_drought_combinations <- tidyr::expand_grid(
  all_gauges,
  non_drought_streamflow_models,
  all_objective_functions
) |> 
  unclass() |> 
  unname() # names breaks pmap? I don't know why - remove them


drought_combinations <- tidyr::expand_grid(
  drought_gauges,
  drought_streamflow_models,
  all_objective_functions
) |> 
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
  .l = non_drought_combinations,
  .f = numerical_optimiser_objects,
  data = data,
  start_stop = start_stop_indexes
)



drought_ready_for_optimisation <- pmap(
  .l = drought_combinations,
  .f = numerical_optimiser_objects,
  data = data,
  start_stop = start_stop_indexes
)



# Split ready_for_optimisation objects into chunks to avoid exceeding RAM ------
# Subject to change. I need to work out RAM requirements
CHUNK_SIZE <- 100

chunked_ready_for_optimisation <- split( # required for chunking in parallel
  ready_for_optimisation,
  ceiling(seq_along(ready_for_optimisation) / CHUNK_SIZE)
)

chunked_drought_ready_for_optimisation <- split( # required for chunking in parallel
  drought_ready_for_optimisation,
  ceiling(seq_along(drought_ready_for_optimisation) / CHUNK_SIZE)
)




# Run the optimiser (calibration) ----------------------------------------------

## Non-drought optimising ======================================================
plan(multisession, workers = length(availableWorkers())) # set once for furrr
iwalk(
  .x = chunked_ready_for_optimisation,
  .f = run_and_save_chunks_optimiser_parallel, 
  optimiser = my_dream,
  save_streamflow = FALSE,
  save_sequences = TRUE,
  is_drought = FALSE
)


gc()


## Drought optimising ==========================================================
plan(multisession, workers = length(availableWorkers())) # set once for furrr
iwalk(
  .x = chunked_drought_ready_for_optimisation,
  .f = run_and_save_chunks_optimiser_parallel, 
  optimiser = my_dream,
  save_streamflow = FALSE,
  save_sequences = TRUE,
  is_drought = TRUE
)


gc()





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





# DREAM - t.thin 
# do I want everything to have the same number of sequences?
# or do I want it to scale with parameter number?




# TESTING ----------------------------------------------------------------------
stop_here <- 1
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



# Current seq
PARAMETER_NUMBER <- 9
seq <- round_any(((PARAMETER_NUMBER - 2) ^ 2.5) * 1E4, 1E4, ceiling)
total_save_seq <- 1000
thin <- seq / total_save_seq




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


# Did the results produce 1000 of each parameter? ------------------------------
DREAM_sequence_results <- read_csv("Results/my_dream/DREAM_sequence_results.csv")


x <- DREAM_sequence_results |> 
  summarise(
    n = n(),
    .by = c(gauge, streamflow_model, objective_function, parameter)
  )

# I don't know how thining works. Also if DREAM converges before 
# the max ndraw is meet it will not be 1000 different combinations
