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

CMAES_results <- read_csv("./Results/my_cmaes/CMAES_parameter_results_20241031.csv",
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
CHUNK_SIZE <- 100

chunked_ready_for_optimisation <- split( # required for chunking in parallel
  ready_for_optimisation,
  ceiling(seq_along(ready_for_optimisation) / CHUNK_SIZE)
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



# analysising dream results WILL MOVE TO ANOTHER FILE --------------------------
# how close are the test catchments parameter to CMAES
# histogram
# comparison of streamflow time graphs?
sequences <- read_csv("Results/my_dream/DREAM_sequence_results.csv", show_col_types = FALSE)


### Convert a5 parameter into year #############################################
CO2 <- data |>
  filter(gauge == "102101A") |>
  pull(CO2)

year <- data |>
  filter(gauge == "102101A") |>
  pull(year)




# make this compatible with tables
single_a5_to_year <- function(shifted_CO2_parameter, CO2, year) {
  adjusted_CO2 <- if_else(CO2 - shifted_CO2_parameter < 0, 0, CO2 - shifted_CO2_parameter)

  year_where_CO2_impacts_flow <- year[adjusted_CO2 != 0][1]

  return(year_where_CO2_impacts_flow)
}


shifted_CO2_parameter_into_year_CO2_starts_impacting_flow <- function(shifted_CO2_parameter, CO2, year) {
  map_dbl(
    .x = shifted_CO2_parameter,
    .f = single_a5_to_year,
    CO2 = CO2,
    year = year
  )
}


just_a5 <- sequences |>
  filter(parameter == "a5") |>
  mutate(
    year_where_CO2_impacts_flow = shifted_CO2_parameter_into_year_CO2_starts_impacting_flow(
      shifted_CO2_parameter = parameter_values,
      CO2 = CO2,
      year = year
    )
  )


a5_to_year_plot <- just_a5 |>
  ggplot(aes(x = year_where_CO2_impacts_flow)) +
  geom_histogram(
    binwidth = binwidth_bins(30),
    fill = "grey",
    colour = "black"
  ) +
  # scale_y_sqrt() +
  labs(
    x = "Year where CO2 starts impacting streamflow",
    y = "Frequency"
  ) +
  theme_bw() +
  facet_wrap(~gauge, scales = "free")


a5_to_year_plot




# See if thin.t produces 1000 different parameter combinations -----------------
count_combinations <- sequences |>
  summarise(
    n = n(),
    .by = c(gauge, streamflow_model, objective_function, parameter)
  )

# I don't know how thining works. Also if DREAM converges before
# the max ndraw is meet it will not be 1000 different combinations



# Plot the histograms of parameter values ======================================
plot_parameter_histogram <- function(streamflow_model, objective_function, sequence_results) {
  sequence_results |>
    filter(streamflow_model == {{ streamflow_model }}) |>
    filter(objective_function == {{ objective_function }}) |>
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
      y = "Frequency",
      title = paste0("Streamflow Model: ", streamflow_model, "\nObjective Function: ", objective_function)
    ) +
    facet_grid(gauge ~ parameter, scales = "free")
}


models <- unique(pull(sequences, streamflow_model))
objfun <- unique(pull(sequences, objective_function))

histogram_plot_1 <- plot_parameter_histogram(
  streamflow_model = models[4],
  objective_function =  objfun,
  sequence_results = sequences
)

histogram_plot_2 <- plot_parameter_histogram(
  streamflow_model = models[2],
  objective_function =  objfun,
  sequence_results = sequences
)

histogram_plot_1
histogram_plot_2



# Compare streamflow timeseries ================================================
DREAM_parameter_results <- read_csv("Results/my_dream/DREAM_parameter_results.csv", show_col_types = FALSE, col_types = "cccccdddcl") |>
  drop_na()
# these will have to be inputted into the respective models
# Temporary and hardcoded

CMAES_streamflow <- read_csv("Results/my_cmaes/CMAES_streamflow_results_20241028.csv", show_col_types = FALSE)


DREAM_join <- DREAM_parameter_results |>
  select(c(gauge, streamflow_model, objective_function)) |>
  distinct()


CMAES_filtered_streamflow <- CMAES_streamflow |>
  semi_join(DREAM_join, by = join_by(gauge, streamflow_model, objective_function))


parameters <- DREAM_parameter_results |>
  filter(gauge == "408202") |>
  filter(streamflow_model == "streamflow_model_drought_separate_shifted_CO2_seasonal_ratio_auto") |>
  pull(parameter_value)


gauge <- "408202"

catchment_data <- gauge |>
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  )

catchment_data$stop_start_data_set <- NULL

catchment_data_altered <- unlist(catchment_data, recursive = FALSE)

names(catchment_data_altered) <- c(
  "gauge_ID",
  "contains_drought",
  "year",
  "precipitation",
  "observed_boxcox_streamflow",
  "is_drought_year",
  "CO2",
  "seasonal_ratio"
)

streamflow <- streamflow_model_drought_separate_shifted_CO2_seasonal_ratio_auto(
  catchment_data = catchment_data,
  parameter_set = parameters
)




# Compare parameter results ====================================================
CMAES_parameter_results <- read_csv("Results/my_cmaes/CMAES_parameter_results_20241028.csv", show_col_types = FALSE)

compare_parameter_results <- CMAES_parameter_results |>
  right_join(DREAM_parameter_results, by = join_by(gauge, streamflow_model, objective_function)) |>
  mutate(
    AIC_diff = AIC.x - AIC.y
  )



# TESTING ----------------------------------------------------------------------
# make sure everything gets a good exit message

test_suite <- best_model_combination_per_catchment |> 
  distinct(
    streamflow_model_name,
    objective_function_name,
    .keep_all = TRUE
  )

# hard code into the test below:
# - works for gauge = 105102A, model = streamflow_model_drought_precip_only, obj = constant_sd_objective_function

tic <- as.numeric(Sys.time())

gauge <- "146095A"

dream_example <- gauge |>
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
  my_dream(print_monitor = TRUE) |>
  result_set()

toc <- round(as.numeric(Sys.time()) - tic, 1)

cat(toc, "sec")

dream_parameters <- dream_example |>
  parameters_summary()

dream_example |>
  plot()






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
seq <- round_any(((PARAMETER_NUMBER - 2)^2.5) * 1E4, 1E4, ceiling)
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
y <- ((x^2) - (x + 1)) * 1E4
plot(x, y)

yy <- round_any(((x - 2)^2.5) * 1E4, 1E4, ceiling)
yy
# non-linear relationship works with gauge 407214 and 3, 4 parameters
# test for 5, 6, 7, 8? Go straight to 8? Good idea to test all.
