# Fitting streamflow models to catchments
cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, cmaesr, smoof, tictoc, furrr, parallel, truncnorm, sloop, tictoc)

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


obs_p_and_q <- data |>
  mutate(
    obs_q_greater_than_p = q_mm > p_mm
  ) |> 
  filter(obs_q_greater_than_p) |> 
  relocate(
    q_mm,
    .after = p_mm
  ) |> 
  select(gauge, year, p_mm, q_mm)


## Utility functions ===========================================================
source("./Functions/utility.R")

source("./Functions/boxcox_transforms.R") 

## Import streamflow functions =================================================
source("./Functions/streamflow_models.R")
source("./Functions/parameter_transformations.R")
source("./Functions/catchment_data_blueprint.R")
source("./Functions/cmaes_dream_summaries.R")
source("./Functions/objective_functions.R")
source("./Functions/numerical_optimiser_setup.R")
source("./Functions/generic_functions.R")
source("./Functions/my_cmaes.R")
source("./Functions/objective_function_setup.R")
source("./Functions/result_set.R")


# Gauges to test:
# 1. Streamflow w/ CO2 off > streamflow = "207015" - works. streamflow is always less than precip. Deals with zero values better as well - in slides streamflow_model_slope_shifted_CO2_seasonal_ratio_auto
# 2. 219001 (neg intercept change) - get model - streamflow_model_intercept_shifted_CO2_seasonal_ratio_auto - in slides
# 3. Good gauge - 137101A streamflow_model_drought_precip_seasonal_ratio
# 4. Good gauge - low flow - 613146 - streamflow_model_intercept_shifted_CO2_seasonal_ratio_auto

gauge <- "207015" 

example_catchment <- gauge |>
  catchment_data_blueprint_v2(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  ) 


plot_catchment_data_v2(example_catchment, type = "streamflow-time")


cmaes_results <- example_catchment |> 
  numerical_optimiser_setup_vary_inputs_v2(
    streamflow_model = streamflow_model_slope_shifted_CO2_seasonal_ratio_auto,
    objective_function = constant_sd_objective_function_log_sinh,
    bounds_and_transform_method = make_default_bounds_and_transform_methods(example_catchment),
    minimise_likelihood = TRUE
  ) |>
  my_cmaes(
    print_monitor = TRUE
  ) 


summarised_results <- cmaes_results |> 
  result_set_v2()


parameters_summary(summarised_results)

plot_result_set_v2(summarised_results, type = "rainfall-runoff")
plot_result_set_v2(summarised_results, type = "streamflow-time") 



# What does the streamflow time plot look like when a3 is turned off? ----------
best_parameters_no_CO2 <- summarised_results$best_parameter_set
best_parameters_no_CO2[4] <- 0 # hard coded will need to change - same with a and b

x <- example_catchment |> 
  streamflow_model_slope_shifted_CO2_seasonal_ratio_auto(
    parameter_set = best_parameters_no_CO2
  ) |> 
  mutate(
    modelled_streamflow_CO2_off = inverse_log_sinh_transform(a = best_parameters_no_CO2[8], b = best_parameters_no_CO2[9], z = streamflow_results)
  ) |> 
  select(!streamflow_results) |> 
  cbind("modelled_streamflow_CO2_on" = summarised_results$optimised_modelled_streamflow_realspace)


log_sinh_plot <- x |> 
  pivot_longer(
    cols = contains("streamflow"),
    names_to = "streamflow_type",
    values_to = "streamflow"
  ) |> 
  ggplot(aes(x = year, y = streamflow, colour = streamflow_type)) +
  geom_line(linewidth = 0.8) +
  geom_point(alpha = 0.9, size = 2) +
  geom_line(aes(x = year, y = precipitation), colour = "black", linetype = "dashed", linewidth = 0.8) +
  # add secondary y axis
  scale_y_continuous(
    name = "Streamflow (mm)", # Features of the first axis
    sec.axis = sec_axis(transform = ~.*1, name = "Precipitation (mm)") # Add a second axis and specify its features
  ) +
  labs(
    x = "Year",
    y = "Streamflow (mm)",
    colour = NULL,
    title = "Log-sinh transform"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )


log_sinh_plot

