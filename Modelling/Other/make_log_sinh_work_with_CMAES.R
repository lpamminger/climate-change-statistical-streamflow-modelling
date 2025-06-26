# Comparing box-cox and log-sinh streamflow transformation


# Outcome: log-sinh performs better
# Figures:





# Fitting streamflow models to catchments

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
  show_col_types = FALSE
) |>
  mutate(
    year = as.integer(year)
  )

gauge_information <- readr::read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE
)





## Import functions ============================================================
source("./Functions/utility.R")
source("./Functions/boxcox_logsinh_transforms.R")
source("./Functions/streamflow_models.R")
source("./Functions/parameter_transformations.R")
source("./Functions/catchment_data_blueprint.R")
source("./Functions/cmaes_dream_summaries.R")
source("./Functions/objective_functions.R")
source("./Functions/numerical_optimiser_setup.R")
source("./Functions/generic_functions.R")
source("./Functions/CMAES.R")
source("./Functions/objective_function_setup.R")
source("./Functions/result_set.R")


# TODO:
# - make this work with the new numerical optimiser






# Get gauges and data from previous attempts - hard code -----------------------
gauge_streamflow_model_combinations <- tribble(
  ~gauge, ~streamflow_model,
  "207015", streamflow_model_slope_shifted_CO2_seasonal_ratio_auto, # large peak in streamflow when CO2 turned off
  "219001", streamflow_model_intercept_shifted_CO2_seasonal_ratio_auto, # observed streamflow > precip
  "112101B", streamflow_model_slope_shifted_CO2_seasonal_ratio_auto, # positive slope
  "225213", streamflow_model_intercept_shifted_CO2_auto, # positive intercept
  "G0050115", streamflow_model_intercept_shifted_CO2_seasonal_ratio, # very low flow
  "603005", streamflow_model_drought_slope_shifted_CO2_seasonal_ratio_auto, # drought
  "405240", streamflow_model_intercept_shifted_CO2_seasonal_ratio_auto # very high evidence ratio
)

gauges <- gauge_streamflow_model_combinations |> pull(gauge)




# Make catchment data ----------------------------------------------------------
catchment_data_per_gauge <- map(
  .x = gauges,
  .f = catchment_data_blueprint,
  observed_data = data,
  start_stop_indexes = start_stop_indexes
)


# plot(catchment_data_per_gauge[[6]], type = "streamflow-time")
# plot(catchment_data_per_gauge[[5]], type = "rainfall-runoff")





# Calibrate data ---------------------------------------------------------------

streamflow_transform_methods <- c(
  log_sinh_transform,
  boxcox_transform
)

model_combinations <- expand_grid(catchment_data_per_gauge, streamflow_transform_methods) |>
  add_column(
    "gauge" = rep(gauges, each = 2), # hard coded - only works if there are two objective functions
    .before = 1
  ) |>
  left_join(
    gauge_streamflow_model_combinations,
    by = join_by(gauge)
  )


bounds_and_transform_per_gauge <- map(
  .x = model_combinations |> pull(catchment_data_per_gauge),
  .f = make_default_bounds_and_transform_methods
)

model_combinations <- model_combinations |>
  add_column("objective_function" = list(constant_sd_objective_function)) |>
  add_column(bounds_and_transform_per_gauge) |>
  add_column("streamflow_transform_method_offset" = rep(c(0, 1), times = 7))


## pmap out cmaes results ######################################################
# order matters --> streamflow_model, objective_function, streamflow_transform_method, catchment_data, bounds_and_transform_method

map_list <- list(
  model_combinations |> pull(streamflow_model),
  model_combinations |> pull(objective_function),
  model_combinations |> pull(streamflow_transform_methods),
  model_combinations |> pull(catchment_data_per_gauge),
  model_combinations |> pull(bounds_and_transform_per_gauge),
  model_combinations |> pull(streamflow_transform_method_offset)
)

numerical_optimiser_setup_combinations <- pmap(
  .l = map_list,
  .f = numerical_optimiser_setup,
  minimise_likelihood = TRUE
)


plan(multisession, workers = length(availableWorkers())) # set once for furrr
cmaes_results <- future_map(
  .x = numerical_optimiser_setup_combinations,
  .f = CMAES,
  print_monitor = FALSE,
  .options = furrr_options(seed = 1L)
)


## Summarise results ###########################################################
summarise_cmaes_results <- map(
  .x = cmaes_results,
  .f = result_set
)

plot(summarise_cmaes_results[[1]], type = "rainfall-runoff")
plot(summarise_cmaes_results[[2]], type = "streamflow-time")

best_parameters <- map(
  .x = summarise_cmaes_results,
  .f = parameters_summary
) |>
  list_rbind() 

check_bounds <- best_parameters |> 
  filter(!is.na(near_bounds))



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

  # Transform to realspace using select_streamflow_transform_method
  inverse_streamflow_transform_method_name <- paste0("inverse_", result_set$numerical_optimiser_setup$streamflow_transform_method()$name)

  realspace_modelled_streamflow_CO2_off <- select_streamflow_transform_method(
    timeseries = as.matrix(transformed_CO2_off_modelled_streamflow, ncol = 1), # this function requires matrix input for timeseries and parameter_set
    parameter_set = as.matrix(best_parameters_CO2_off, ncol = 1),
    streamflow_transform_method = match.fun(FUN = inverse_streamflow_transform_method_name), # get the inverse from name then match the function
    offset = result_set$numerical_optimiser_setup$streamflow_transform_method_offset # get result set
  ) |>
    as.numeric() # convert back into vector

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
      "realspace_mod_flow_CO2_off" = turn_off_CO2_component(result_set)$realspace_CO2_off,
      "streamflow_transform_method" = result_set$numerical_optimiser_setup$streamflow_transform_method()$name,
      "gauge" = result_set$numerical_optimiser_setup$catchment_data$gauge_ID
    ) |>
    select(-c(is_drought_year, CO2, seasonal_ratio)) |>
    rename(
      realspace_obs_flow = observed_streamflow
    )
}




### Produce data for plotting ##################################################

plotting_data <- map2(
  .x = summarise_cmaes_results,
  .y = turn_off_CO2_component_results,
  .f = result_set_to_plotting_data
) |>
  list_rbind()



# Compare boxcox and log-sinh AIC and residuals values -------------------------
residuals_check <- plotting_data |>
  mutate(
    obs_minus_mod_residual = transformed_obs_flow - transformed_mod_flow_CO2_on,
    abs_obs_minus_mod_residual = abs(obs_minus_mod_residual)
    ) |>
  summarise(
    sum_abs_residual = sum(abs_obs_minus_mod_residual),
    .by = c(gauge, streamflow_transform_method)
  )


get_likelihood_information <- function(result_set) {
  list(
    "gauge" = result_set$numerical_optimiser_setup$catchment_data$gauge_ID,
    "loglikelihood" = result_set$LL_best_parameter_set,
    "AIC" = result_set$AIC_best_parameter_set,
    "streamflow_transform_method" = result_set$numerical_optimiser_setup$streamflow_transform_method()$name
  ) |>
    as_tibble()
}

loglikelihood_information <- map(
  .x = summarise_cmaes_results,
  .f = get_likelihood_information
) |>
  list_rbind()


numerical_comparison_boxcox_logsinh <- residuals_check |>
  left_join(
    loglikelihood_information,
    by = join_by(gauge, streamflow_transform_method)
  )



# What does the streamflow time plot look like when a3 is turned off? ----------

## Transformed streamflow time =================================================
transformed_streamflow_time_plot <- plotting_data |>
  select(!contains("realspace")) |>
  pivot_longer(
    cols = contains("transformed"),
    names_to = "streamflow_type",
    values_to = "transformed_streamflow"
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
  filter(streamflow_type != "Modelled Streamflow CO2 Off") |> # ignore for now
  ggplot(aes(x = year, y = transformed_streamflow, colour = streamflow_type)) +
  geom_line(linewidth = 0.8) +
  geom_point(alpha = 0.9, size = 2) +
  labs(
    x = "Year",
    y = "Transformed Streamflow",
    colour = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  ) +
  facet_wrap(~ gauge + streamflow_transform_method, ncol = 2, nrow = 7, scales = "free_y")

ggsave(
  filename = "boxcox_logsinh_transformed_timeseries_comparison.pdf",
  plot = transformed_streamflow_time_plot,
  device = "pdf",
  path = "Modelling/Other",
  width = 320,
  height = 420,
  units = "mm"
)

## Streamflow time =============================================================
streamflow_time_plot <- plotting_data |>
  select(!contains("transformed")) |>
  pivot_longer(
    cols = contains("realspace"),
    names_to = "streamflow_type",
    values_to = "realspace_streamflow"
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
  geom_point(alpha = 0.9, size = 1.2) +
  geom_line(aes(x = year, y = precipitation), colour = "black", linetype = "dashed", linewidth = 0.8) +
  # geom_label(
  #  aes(x = x_pos, y = y_pos, label = AIC),
  #  data = AIC_streamflow_transform_comparison,
  #  inherit.aes = FALSE
  # ) +
  labs(
    x = "Year",
    y = "Streamflow (mm)",
    colour = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  ) +
  facet_wrap(~ gauge + streamflow_transform_method, ncol = 2, nrow = 7, scales = "free_y")

ggsave(
  filename = "boxcox_logsinh_realspace_timeseries_comparison.pdf",
  plot = streamflow_time_plot,
  device = "pdf",
  path = "Modelling/Other",
  width = 320,
  height = 420,
  units = "mm"
)


## Rainfall-runoff =============================================================
rainfall_runoff_plot <- plotting_data |>
  select(!contains("realspace")) |>
  pivot_longer(
    cols = contains("transformed"),
    names_to = "streamflow_type",
    values_to = "transformed_streamflow"
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
  arrange(gauge) |>
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
  facet_wrap(~ gauge + streamflow_transform_method, nrow = 7, ncol = 2, scales = "free")

ggsave(
  filename = "boxcox_logsinh_rainfall_runoff_comparison.pdf",
  plot = rainfall_runoff_plot,
  device = "pdf",
  path = "Modelling/Other",
  width = 250,
  height = 594,
  units = "mm"
)



# Why are the AIC so different? ------------------------------------------------

# TODO:
# - get sd for each gauge and streamflow_transform method
# - join sd to plotting data
# - created lower and upper bound around the calibrated terms using +/-2 sd
# - plot

# Check error bar using sd around transformed timeseries -----------------------

only_sd <- best_parameters |>
  filter(
    parameter == "sd"
  ) |> 
  select(-c(streamflow_model, objective_function, optimiser, exit_message))

compare_sd <- numerical_comparison_boxcox_logsinh |>
  left_join(
    only_sd,
    by = join_by(gauge, AIC, loglikelihood)
  ) 

sd_plotting_data <- plotting_data |> 
  left_join(
    compare_sd,
    by = join_by(gauge, streamflow_transform_method)
  ) |> 
  select(
    c(year, 
      precipitation, 
      transformed_obs_flow, 
      transformed_mod_flow_CO2_on, 
      streamflow_transform_method, 
      gauge, 
      parameter_value
      )
    ) |> 
  rename(
    sd = parameter_value
  ) |> 
  # add lower and upper bound uncertainty for modelled streamflow
  mutate(
    lower_bound_mod_flow = transformed_mod_flow_CO2_on - (2 * sd),
    upper_bound_mod_flow = transformed_mod_flow_CO2_on + (2 * sd)
  )

  
sd_uncertainty_bars <- sd_plotting_data |> 
  pivot_longer(
    cols = starts_with("transformed"),
    names_to = "observed_or_modelled",
    values_to = "transformed_streamflow"
  ) |> 
  mutate(
    lower_bound_mod_flow = if_else(observed_or_modelled == "transformed_obs_flow", NA, lower_bound_mod_flow),
    upper_bound_mod_flow = if_else(observed_or_modelled == "transformed_obs_flow", NA, upper_bound_mod_flow)
  ) |> 
  # Make the names nice
  mutate(
    observed_or_modelled = if_else(observed_or_modelled == "transformed_obs_flow", "Observed Streamflow", "Modelled Streamflow") 
  ) |> 
  ggplot(aes(x = year, y = transformed_streamflow, colour = observed_or_modelled, fill = observed_or_modelled)) +
  geom_ribbon(
    aes(ymin = lower_bound_mod_flow, ymax = upper_bound_mod_flow),
    na.rm = FALSE, # I have mannually set obs flow to NA
    alpha = 0.25,
    colour = NA
  ) +
  geom_line() +
  scale_colour_brewer(palette = "Set1") +
  labs(
    x = "Year",
    y = "Transformed Streamflow",
    colour = NULL,
    fill = NULL
  ) +
  theme_bw() +
  facet_wrap(~gauge + streamflow_transform_method, scales = "free", ncol = 2) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  filename = "sd_uncertainty_transformed_streamflow_time.pdf",
  plot = sd_uncertainty_bars,
  device = "pdf",
  path = "Modelling/Other",
  width = 250,
  height = 594,
  units = "mm"
)


# qq plots ---------------------------------------------------------------------
residual_plotting_data <- plotting_data |>
  mutate(obs_minus_mod_residual = transformed_obs_flow - transformed_mod_flow_CO2_on) 


qqplots <- residual_plotting_data |> 
  ggplot(aes(sample = obs_minus_mod_residual)) +
  geom_qq_line() +
  geom_qq() +
  labs(
    x = "Theoretical",
    y = "Sample (obs - mod streamflow residual)"
  ) +
  theme_bw() +
  facet_wrap(~gauge + streamflow_transform_method, ncol = 2, scales = "free")
  

ggsave(
  filename = "boxcox_logsinh_qqplot.pdf",
  plot = qqplots,
  device = "pdf",
  path = "Modelling/Other",
  width = 250,
  height = 594,
  units = "mm"
)


# Look at objective_function ---------------------------------------------------

test_catchment_data_log_sinh <- plotting_data |> 
  filter(gauge == "207015") |> 
  filter(streamflow_transform_method == "log_sinh_transform")

test_catchment_data_boxcox <- plotting_data |> 
  filter(gauge == "207015") |> 
  filter(streamflow_transform_method == "boxcox_transform")

# If the gauge has chunks I must do it by chunks - test gauge without chunks first
# Gauge 207015 has a single chunk

# Sticking results in get the same calibrated results

# The log-likelihood is calculated using different observed streamflow values
# Does this mean they can be compared?

# AIC myths - https://robjhyndman.com/hyndsight/aic/

# Log-sinh
constant_sd_objective_function(
  modelled_streamflow = test_catchment_data_log_sinh |> pull(transformed_mod_flow_CO2_on) |> as.matrix(ncol = 1), 
  observed_streamflow = test_catchment_data_log_sinh |> pull(transformed_obs_flow) |> as.matrix(ncol = 1),
  parameter_set = as.matrix(summarise_cmaes_results[[1]]$best_parameter_set, ncol = 1),
  apply_truncnorm = FALSE
)

# Box-Cox
constant_sd_objective_function(
  modelled_streamflow = test_catchment_data_boxcox |> pull(transformed_mod_flow_CO2_on) |> as.matrix(ncol = 1), 
  observed_streamflow = test_catchment_data_boxcox |> pull(transformed_obs_flow) |> as.matrix(ncol = 1),
  parameter_set = as.matrix(summarise_cmaes_results[[2]]$best_parameter_set, ncol = 1),
  apply_truncnorm = TRUE
)



# 1. generate observed streamflow - realspace
set.seed(1)
observed_streamflow <- runif(n = 100, min = 25, 400)

# 2. generate modelled streamflow - realspace
modelled_streamflow <- observed_streamflow + rnorm(n = 100, mean = 0, sd = 20)

# 3. transform modelled streamflow into boxcox and log-sinh space
boxcox_parameters <- c(1, 0.1) # c(sd, lambda)
log_sinh_parameters <- c(75, 0.1, 0.5) # c(sd, a, b)

boxcox_modelled_streamflow <- boxcox_transform(y = modelled_streamflow, lambda = boxcox_parameters[2], lambda_2 = 1)
log_sinh_modelled_streamflow <- log_sinh_transform(a = log_sinh_parameters[2], b = log_sinh_parameters[3], y = modelled_streamflow, offset = 300)

# 4. calculate log-likelihoods
boxcox_loglikelihood <- constant_sd_boxcox_objective_function(
  modelled_streamflow = as.matrix(boxcox_modelled_streamflow, ncol = 1),
  observed_streamflow = as.matrix(observed_streamflow, ncol = 1),
  parameter_set = as.matrix(boxcox_parameters, ncol = 1)
)

log_sinh_loglikelihood <- constant_sd_log_sinh_objective_function(
  modelled_streamflow = as.matrix(log_sinh_modelled_streamflow, ncol = 1),
  observed_streamflow = as.matrix(observed_streamflow, ncol = 1),
  parameter_set = as.matrix(log_sinh_parameters, ncol = 1)
)

cat("boxcox LL =", boxcox_loglikelihood, "\nlog_sinh LL =", log_sinh_loglikelihood)
