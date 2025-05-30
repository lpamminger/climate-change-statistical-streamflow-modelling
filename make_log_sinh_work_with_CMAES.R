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



# Get gauges and data from previous attempts - hard code -----------------------
gauge_streamflow_model_combinations <- tribble(
  ~gauge,    ~streamflow_model,
  "207015",  streamflow_model_slope_shifted_CO2_seasonal_ratio_auto,
  "219001",  streamflow_model_intercept_shifted_CO2_seasonal_ratio_auto
  # add 4 more gauges
  # positve slope, positive intercept, drought, very low flow
)

gauges <- gauge_streamflow_model_combinations |> pull(gauge)


# Make catchment data ----------------------------------------------------------
catchment_data_per_gauge <- map(
  .x = gauges,
  .f = catchment_data_blueprint_v2,
  observed_data = data,
  start_stop_indexes = start_stop_indexes
)


#plot_catchment_data_v2(catchment_data_per_gauge[[2]], type = "streamflow-time")


# Calibrate data ---------------------------------------------------------------

objective_function_per_gauge <- c(
  constant_sd_boxcox_objective_function,
  constant_sd_log_sinh_objective_function
)


## Add streamflow models based on the gauges ###################################


gauge_objective_function_combinations <- expand_grid(catchment_data_per_gauge, objective_function_per_gauge) |>
  add_column(
    "gauge" = rep(gauges, each = 2), # hard coded - only works if there are two objective functions
    .before = 1
  ) |> 
  left_join(
    gauge_streamflow_model_combinations,
    by = join_by(gauge)
  ) 

get_objective_function_name <- function(objective_function) {
  objective_function()$name
}

objective_function_names <- map(
  .x = gauge_objective_function_combinations |> pull(objective_function_per_gauge),
  .f = get_objective_function_name
) 

bounds_and_transform_per_gauge <- map(
  .x = gauge_objective_function_combinations |> pull(catchment_data_per_gauge),
  .f = make_default_bounds_and_transform_methods
)

gauge_objective_function_combinations <- gauge_objective_function_combinations |> 
  add_column(bounds_and_transform_per_gauge) |> 
  add_column("transform_method" = objective_function_names)


## pmap out cmaes results ######################################################
map_list <- list(
  gauge_objective_function_combinations |> pull(streamflow_model),
  gauge_objective_function_combinations |> pull(objective_function_per_gauge),
  gauge_objective_function_combinations |> pull(catchment_data_per_gauge),
  gauge_objective_function_combinations |> pull(bounds_and_transform_per_gauge)
)

numerical_optimiser_setup_combinations <- pmap(
  .l = map_list,
  .f = numerical_optimiser_setup_v2,
  minimise_likelihood = TRUE
)

# Make this parallel - later
cmaes_results <- map(
  .x = numerical_optimiser_setup_combinations,
  .f = my_cmaes,
  print_monitor = TRUE
)


## Summarise results ###########################################################
summarise_cmaes_results <- map(
  .x = cmaes_results,
  .f = result_set_v2
)


#plot_result_set_v2(summarise_cmaes_results[[4]], type = "rainfall-runoff")
#plot_result_set_v2(summarise_cmaes_results[[4]], type = "streamflow-time")




# Turn off CO2 component - add turn off info to summarise_cmaes_results --------

# need code to turn-off the CO2 term and plot the results
# - need transformed streamflow with CO2 term off
# - need realspace streamflow with CO2 term off
# Make function with result_set as input
# 1. get best parameters
# 2. find a3 and turn it off
# 3. re-run model to get transformed streamflow
# 4. turn transformed streamflow into realspace
turn_off_CO2_component <- function(result_set) {
  stopifnot(s3_class(result_set)[1] == "result_set")

  # Get best parameters
  best_parameters_CO2_on <- result_set$best_parameter_set

  # Find a3 and turn it off
  best_parameters_CO2_off <- best_parameters_CO2_on
  a3_term_position <- str_detect(names(best_parameters_CO2_off), "a3")
  best_parameters_CO2_off[a3_term_position] <- 0

  # Re-run model with CO2 off
  transformed_CO2_off_modelled_streamflow <- result_set$numerical_optimiser_setup$streamflow_model(
    catchment_data = result_set$numerical_optimiser_setup$catchment_data,
    parameter_set = best_parameters_CO2_off
  ) |>
    pull(streamflow_results)

  # Transform to realspace
  if (result_set$numerical_optimiser_setup$objective_function()$name == "constant_sd_boxcox_objective_function") {
    best_lambda <- best_parameters_CO2_off[length(best_parameters_CO2_off)]

    realspace_modelled_streamflow_CO2_off <- boxcox_inverse_transform(
      yt = transformed_CO2_off_modelled_streamflow,
      lambda = best_lambda,
      lambda_2 = 1
    ) 
    
  } else if (result_set$numerical_optimiser_setup$objective_function()$name == "constant_sd_log_sinh_objective_function") {
    best_a <- best_parameters_CO2_off[length(best_parameters_CO2_off) - 1]
    best_b <- best_parameters_CO2_off[length(best_parameters_CO2_off)]

    realspace_modelled_streamflow_CO2_off <- inverse_log_sinh_transform(
      a = best_a,
      b = best_b,
      z = transformed_CO2_off_modelled_streamflow
    ) 
    
  } else {
    stop("Name of objective function not found")
  }


  # Return a tibble of both results
  return(list("transformed_CO2_off" = transformed_CO2_off_modelled_streamflow, "realspace_CO2_off" = realspace_modelled_streamflow_CO2_off))
}


## Get turn of CO2 streamflow values ###########################################
turn_off_CO2_component_results <- map(
  .x = summarise_cmaes_results,
  .f = turn_off_CO2_component
)



## Turn into a tibble for plotting #############################################
# Tibble requirements:
# 1. year, precip, etc.
# 2. observed streamflow realspace and transformed
# 3. modelled streamflow CO2 on realspace and transformed
# 4. modelled streamflow CO2 off realspace and transformed


result_set_to_plotting_data <- function(result_set, turn_off_CO2_component) {
  result_set$numerical_optimiser_setup$catchment_data$stop_start_data_set |>
    list_rbind() |>
    cbind(
      "transformed_obs_flow" = result_set$transformed_observed_streamflow,
      "transformed_mod_flow_CO2_on" = result_set$optimised_modelled_streamflow_transformed_space,
      "realspace_mod_flow_CO2_on" = result_set$optimised_modelled_streamflow_realspace,
      "transformed_mod_flow_CO2_off" = turn_off_CO2_component(result_set)$transformed_CO2_off,
      "realspace_mod_flow_CO2_off" = turn_off_CO2_component(result_set)$realspace_CO2_off
    ) |> 
    select(-c(is_drought_year, CO2, seasonal_ratio)) |> 
    rename(
      realspace_obs_flow = observed_streamflow
    )
  
}




## Produce data for plotting ###################################################
gauge_transform_method_key <- gauge_objective_function_combinations |> 
  select(gauge, transform_method) |> 
  unite(col = "gauge_transform", sep = "-") |> 
  pull(gauge_transform)


plotting_data <- map2(
  .x = summarise_cmaes_results,
  .y = turn_off_CO2_component_results,
  .f = result_set_to_plotting_data
) |>
  `names<-`(gauge_transform_method_key) |>
  list_rbind(names_to = "gauge_transform") |> 
  separate_wider_delim(
    cols = gauge_transform,
    delim = "-",
    names = c("gauge", "transform_method")
  )



# What does the streamflow time plot look like when a3 is turned off? ----------

## Rainfall-runoff #############################################################
plotting_data |>
  select(!contains("transformed")) |> 
  pivot_longer(
    cols = contains("realspace"),
    names_to = "streamflow_type",
    values_to = "realspace_streamflow"
  ) |> 
  mutate(
    transform_method = if_else(transform_method == "constant_sd_boxcox_objective_function", "Box-Cox Transform", "Log-Sinh Transform")
  ) |> 
  mutate(
    streamflow_type = case_when(
      streamflow_type == "realspace_mod_flow_CO2_off" ~ "Modelled Streamflow CO2 Off",
      streamflow_type == "realspace_mod_flow_CO2_on" ~ "Modelled Streamflow CO2 On",
      streamflow_type == "realspace_obs_flow" ~ "Observed Streamflow",
      .default = NA
    ),
    streamflow_type = factor(streamflow_type, levels = c("Observed Streamflow", "Modelled Streamflow CO2 On", "Modelled Streamflow CO2 Off"))
  ) |> 
  ggplot(aes(x = year, y = realspace_streamflow, colour = streamflow_type)) +
  geom_line(linewidth = 0.8) +
  geom_point(alpha = 0.9, size = 2) +
  geom_line(aes(x = year, y = precipitation), colour = "black", linetype = "dashed", linewidth = 0.8) +
  labs(
    x = "Year",
    y = "Streamflow (mm)",
    colour = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  ) +
  facet_grid(gauge ~ transform_method, scales = "free_y")



## Streamflow-time #############################################################
plotting_data |> 
  select(!contains("realspace")) |> 
  pivot_longer(
    cols = contains("transformed"),
    names_to = "streamflow_type",
    values_to = "transformed_streamflow"
  ) |> 
  mutate(
    transform_method = if_else(transform_method == "constant_sd_boxcox_objective_function", "Box-Cox Transform", "Log-Sinh Transform")
  ) |> 
  mutate(
    streamflow_type = case_when(
      streamflow_type == "transformed_mod_flow_CO2_off" ~ "Modelled Streamflow CO2 Off",
      streamflow_type == "transformed_mod_flow_CO2_on" ~ "Modelled Streamflow CO2 On",
      streamflow_type == "transformed_obs_flow" ~ "Observed Streamflow",
      .default = NA
    ),
    streamflow_type = factor(streamflow_type, levels = c("Observed Streamflow", "Modelled Streamflow CO2 On", "Modelled Streamflow CO2 Off"))
  ) |>
  ggplot(aes(x = precipitation, y = transformed_streamflow, colour = streamflow_type)) +
  geom_point() +
  geom_smooth(
    method = lm,
    formula = y ~ x,
    se = FALSE,
    na.rm = FALSE
  ) +
  labs(
    x = "Precipitation",
    y = "Transformed Streamflow",
    colour = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  ) +
  facet_grid(transform_method ~ gauge, scales = "free")
  





