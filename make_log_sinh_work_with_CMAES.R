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
source("./Functions/my_cmaes.R")
source("./Functions/objective_function_setup.R")
source("./Functions/result_set.R")


# Testing changes here ---------------------------------------------------------


# probably problem with matrices in objective_function_setup
gauge <- "112101B" 

example_catchment <- gauge |>
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  ) 


plot(example_catchment, type = "streamflow-time")

results <- example_catchment |>
  numerical_optimiser_setup_vary_inputs(
    streamflow_model = streamflow_model_slope_shifted_CO2_seasonal_ratio_auto,
    objective_function = constant_sd_objective_function,
    streamflow_transform_method = boxcox_transform, #boxcox_transform
    bounds_and_transform_method = make_default_bounds_and_transform_methods(example_catchment),
    minimise_likelihood = TRUE,
    streamflow_transform_method_offset = 300
  ) |>
  my_cmaes(
    print_monitor = TRUE
  ) 


x <- results |> result_set() 
plot(x, type = "streamflow-time")






# Get gauges and data from previous attempts - hard code -----------------------
gauge_streamflow_model_combinations <- tribble(
  ~gauge, ~streamflow_model,
   "207015",  streamflow_model_slope_shifted_CO2_seasonal_ratio_auto, # large peak in streamflow when CO2 turned off
   "219001",  streamflow_model_intercept_shifted_CO2_seasonal_ratio_auto, # observed streamflow > precip
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


#plot(catchment_data_per_gauge[[5]], type = "streamflow-time")
#plot(catchment_data_per_gauge[[5]], type = "rainfall-runoff")





# Calibrate data ---------------------------------------------------------------

## Add to gauge_streamflow_model_combinations tibble:
## - objective function
## - bounds
## - streamflow_transform_methods
## - offset?

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
  add_column("streamflow_transform_method_offset" = rep(c(300, 1), times = 7))


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
  .f = my_cmaes,
  print_monitor = FALSE,
  .options = furrr_options(seed = 1L)
)


## Summarise results ###########################################################
# cmaes_results are correct


summarise_cmaes_results <- map(
  .x = cmaes_results,
  .f = result_set
)

plot(summarise_cmaes_results[[14]], type = "rainfall-runoff")
plot(summarise_cmaes_results[[2]], type = "streamflow-time")

check_bounds <- map(
  .x = summarise_cmaes_results,
  .f = parameters_summary
) |>
  list_rbind()

stop_here()



# TODO:
# 1. fix the code below here to work with the changes





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
  ) |>
  # Manually adjust transform - realspace cannot be less than zero
  mutate(
    realspace_mod_flow_CO2_on = if_else(realspace_mod_flow_CO2_on < 0, 0, realspace_mod_flow_CO2_on),
    realspace_mod_flow_CO2_off = if_else(realspace_mod_flow_CO2_off < 0, 0, realspace_mod_flow_CO2_off)
  )


# Compare boxcox and log-sinh AIC values ---------------------------------------
get_streamflow_transform_method <- function(result_set) {
  # only works for result_set objects
  result_set$numerical_optimiser_setup$objective_function()$name[1] |>
    str_remove("constant_sd_") |>
    str_remove("_objective_function")
}

gauge_AIC_streamflow_transform <- function(result_set) {
  list(
    "gauge" = result_set$numerical_optimiser_setup$catchment_data$gauge_ID,
    "transform_method" = get_streamflow_transform_method(result_set),
    "AIC" = result_set$AIC_best_parameter_set,
    "LL" = result_set$LL_best_parameter_set
  ) |>
    as_tibble()
}





coord_for_labels <- plotting_data |>
  summarise(
    y_pos = max(transformed_obs_flow),
    .by = c(gauge, transform_method)
  ) |>
  mutate(
    transform_method = str_remove(transform_method, "constant_sd_"),
    transform_method = str_remove(transform_method, "_objective_function"),
    y_pos = ceiling(y_pos) - (ceiling(y_pos) * 0.05)
  ) |>
  add_column(
    x_pos = 1965
  )

AIC_streamflow_transform_comparison <- map(
  .x = summarise_cmaes_results,
  .f = gauge_AIC_streamflow_transform
) |>
  list_rbind() |>
  left_join(
    coord_for_labels,
    by = join_by(gauge, transform_method)
  ) |>
  mutate(
    AIC = paste0("AIC: ", round(AIC, digits = 2)),
    transform_method = case_when(
      transform_method == "boxcox" ~ "Box-Cox Transform",
      transform_method == "log_sinh" ~ "Log-Sinh Transform",
      .default = NA
    )
  )

# What does the streamflow time plot look like when a3 is turned off? ----------


## Transformed streamflow time graphs ==========================================
### Look at the results
residuals_check <- plotting_data |>
  mutate(obs_minus_mod_residual = transformed_obs_flow - transformed_mod_flow_CO2_on) |>
  summarise(
    sum_residual = sum(obs_minus_mod_residual),
    .by = c(gauge, transform_method)
  )

### Look at uncertainty parameter
uncertainty_parameter_comparison <- check_bounds |>
  filter(parameter == "sd") |>
  select(-c(streamflow_model, optimiser, exit_message, near_bounds))


transformed_streamflow_time_plot <- plotting_data |>
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
  filter(streamflow_type != "Modelled Streamflow CO2 Off") |> # ingore for now
  ggplot(aes(x = year, y = transformed_streamflow, colour = streamflow_type)) +
  geom_line(linewidth = 0.8) +
  geom_point(alpha = 0.9, size = 2) +
  geom_label(
    aes(x = x_pos, y = y_pos, label = AIC),
    data = AIC_streamflow_transform_comparison,
    inherit.aes = FALSE
  ) +
  labs(
    x = "Year",
    y = "Transformed Streamflow",
    colour = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  ) +
  facet_wrap(~ gauge + transform_method, ncol = 2, nrow = 7, scales = "free_y")

# ggsave(
#  filename = "offset_log_sinh_showing_tim_transform_issue_timeseries.pdf",
#  plot = transformed_streamflow_time_plot,
#  device = "pdf",
#  path = "./Graphs/Supplementary_Figures",
#  width = 320,
#  height = 420,
#  units = "mm"
# )


## Streamflow time =============================================================
streamflow_time_plot <- plotting_data |>
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
  facet_wrap(~ gauge + transform_method, ncol = 2, nrow = 7, scales = "free_y")

# ggsave(
#  filename = "offset_log_sinh_testing_streamflow_transform_methods_timeseries.pdf",
#  plot = streamflow_time_plot,
#  device = "pdf",
#  path = "./Graphs/Supplementary_Figures",
#  width = 320,
#  height = 420,
#  units = "mm"
# )


## Rainfall-runoff =============================================================
rainfall_runoff_plot <- plotting_data |>
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
  facet_wrap(~ gauge + transform_method, nrow = 7, ncol = 2, scales = "free")

# ggsave(
#  filename = "offset_log_sinh_testing_streamflow_transform_methods_rainfall_runoff.pdf",
#  plot = rainfall_runoff_plot,
#  device = "pdf",
#  path = "./Graphs/Supplementary_Figures",
#  width = 250,
#  height = 594,
#  units = "mm"
# )



# Why are the AIC so different? ------------------------------------------------

# Does the different transform methods produce different LL?

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


# Try Tim's suggestion of qq and *2 sd here...
