# DREAM analysis

cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse)


# Import functions -------------------------------------------------------------
source("./Functions/utility.R")
source("./Functions/boxcox_transforms.R")


# Import data ------------------------------------------------------------------
data <- data <- readr::read_csv(
  "Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  col_types = "icdddddddll", # ensuring each column is a the correct type
  show_col_types = FALSE
)


# Import DREAM parameters and sequences ----------------------------------------
DREAM_params <- readr::read_csv(
  "Results/my_dream/parameter_results_chunk_1_20241115.csv",
  show_col_types = FALSE
  )


DREAM_sequences <- readr::read_csv(
  "Results/my_dream/sequences_results_chunk_1_20241115.csv",
  show_col_types = FALSE
  )



# Testing only - remove non-converged catchments -------------------------------
converged_DREAM_params <- DREAM_params |> 
  filter(exit_message == "Convergence criteria reached")

#length(unique(DREAM_params$gauge))
#length(unique(converged_DREAM_params$gauge))
# Only 9 of the 28 catchments have converged using DREAM

# converged_DREAM_params |> filter(parameter == "a5")
# Only 2 of the 28 catchments have an a5 parameter


converged_DREAM_sequences <- DREAM_sequences |> 
  semi_join( # filter join
    converged_DREAM_params,
    by = join_by(gauge)
  )



# For the converged catchments look at parameter distributions -----------------




# For the converged catchments convert a5 into Time of Emergence (ToE) ---------

## Just get a5's
converged_shifted_CO2_param <- converged_DREAM_sequences |> 
  filter(parameter == "a5")


## Convert the a5's into a time of emergence
single_a5_to_time_of_emergence <- function(a5, CO2, year) {
  adjusted_CO2 <- if_else(CO2 - a5 < 0, 0, CO2 - a5)
  
  year_where_CO2_impacts_flow <- year[adjusted_CO2 != 0][1]
  
  return(year_where_CO2_impacts_flow)
}


get_CO2_and_year_per_gauge <- function(gauge, data) {
  
  data |> 
    filter(gauge == {{ gauge }}) |> 
    select(c(year, CO2))
  
}


mutate_compatible_a5_to_time_of_emergence <- function(a5, gauge, data) {
  map_dbl(
    .x = a5,
    .f = single_a5_to_time_of_emergence,
    CO2 = get_CO2_and_year_per_gauge(gauge, data)$CO2,
    year = get_CO2_and_year_per_gauge(gauge, data)$year
  )
}


converged_shifted_CO2_param <- converged_shifted_CO2_param |> 
  mutate(
    ToE = mutate_compatible_a5_to_time_of_emergence(
      a5 = parameter_values,
      gauge = gauge,
      data = data
    )
  )





## Inspect ToE's
count_na <- function(column) {
  sum(is.na(column))
}


  
ToE_range <- converged_shifted_CO2_param |> 
  summarise(
    n = n(),
    count_na = count_na(ToE),
    .by = gauge
  ) |> 
  mutate(
    plotted_ToE = n - count_na,
    percentage_NA = (count_na / n) * 100
  )


## Plot ########################################################################
converged_shifted_CO2_param |> 
  ggplot(aes(x = gauge, y = ToE)) +  # colour = state?
  geom_boxplot(na.rm = TRUE) +
  labs(
    x = "Gauge",
    y = "Time of Emergence"
  ) +
  theme_bw()




#
# !!!!!!!!!!!! OLD !!!!!!!!!!!!!!!!!!! 
#
# Moved from 03-DREAM ----------------------------------------------------------
 




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



# Checking CMAES - can delete/move
# a3 = -25 works (tick) 108003A can have an a3 less than -10 streamflow_model_separate_shifted_CO2	constant_sd_objective_function
# a3 = 50 works. 112102A	a3 > 10 streamflow_model_separate_shifted_CO2_seasonal_ratio_auto	constant_sd_objective_function
# a3 = 20 works. G8150018 a3>10	streamflow_model_separate_shifted_CO2_seasonal_ratio_auto	constant_sd_objective_function


# Results
CMAES_parameter_results <- read_csv(
  "Results/my_cmaes/CMAES_parameter_results_20241101.csv",
  show_col_types = FALSE
)


inspect_parameter_a3 <- CMAES_parameter_results |> 
  filter(parameter == "a3") |> 
  pull(parameter_value) |> 
  range()

# Changes made based on parameter results
# - catchments still hitting a3 bounds - make slightly larger
# - scale_CO2 want to be very close to zero. Cannot be zero due to log transform
# - try to reduce bounds to speed up time. Bounds to reduce:
# a0_d/a0_n, a4, sd

inspect_parameter <- CMAES_parameter_results |> # is the CO2 being shut-off? Sometimes. Hopefully the AIC deems these models as not good.
  filter(parameter == "a5") |> 
  filter(parameter_value > 110)

# Testing
random_test <- CMAES_parameter_results |> 
  filter(streamflow_model == streamflow_model_drought_separate_shifted_CO2_seasonal_ratio_auto()$name) |> 
  slice_sample(
    n = 1
  )

gauge <- "227225A"
example <- gauge |>
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  ) |>
  numerical_optimiser_setup_vary_inputs(
    streamflow_model = streamflow_model_drought_separate_shifted_CO2_seasonal_ratio_auto,
    objective_function = constant_sd_objective_function,
    bounds_and_transform_method = make_default_bounds_and_transform_methods(),
    minimise_likelihood = TRUE
  ) |>
  my_cmaes(print_monitor = TRUE) |>
  result_set() 

x <- example |> parameters_summary()
example |> plot()
# based on trial and error set a3 to -25 to 50
