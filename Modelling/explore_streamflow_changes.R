# Purpose of script:
# 1. create streamflow timeseries for counterfactual
# 2. Plot streamflow time and rainfall-runoff plots of all gauges


# Figures produced in this R file ----------------------------------------------
# 1. Figures/Other/mega_timeseries_plot.pdf
# 2. Figures/Other/mega_rainfall_runoff_plot.pdf


# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, truncnorm, sloop)


# Import functions -------------------------------------------------------------
source("./Functions/utility.R")
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
source("./Functions/boxcox_logsinh_transforms.R")


# Import data ------------------------------------------------------------------
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

parameter_results <- read_csv(
  "./Modelling/Results/CMAES/cmaes_parameter_results.csv",
  show_col_types = FALSE
)

best_CO2_non_CO2_per_gauge <- read_csv(
  "./Modelling/Results/CMAES/best_CO2_non_CO2_per_catchment_CMAES.csv",
  show_col_types = FALSE
)

streamflow_results <- read_csv(
  "./Modelling/Results/CMAES/cmaes_streamflow_results.csv",
  show_col_types = FALSE
)


# Calculate counterfactual for all models with CO2 in them ---------------------

## Set a3 terms (slope and intercept to zero) ==================================
only_CO2_model_parameters <- best_CO2_non_CO2_per_gauge |>
  filter(str_detect(streamflow_model, "CO2"))


set_a3_zero_model_parameters <- only_CO2_model_parameters |>
  mutate(
    altered_parameter = if_else(str_detect(parameter, "a3"), 0, parameter_value)
  )


## Generate streamflow with altered_parameter ==================================
### Build catchment_dataset objects ############################################
gauges <- set_a3_zero_model_parameters |>
  pull(gauge) |>
  unique()

catchment_data <- map(
  .x = gauges,
  .f = catchment_data_blueprint,
  observed_data = data,
  start_stop_indexes = start_stop_indexes
)

### Run the streamflow models ##################################################
### streamflow_models
get_model_name_and_function <- function(model) {
  model_name <- model()$name
  tibble(
    "streamflow_model" = model_name,
    "model_function" = list(model)
  )
}


streamflow_name_and_model_for_joining <- map(
  .x = c(get_non_drought_streamflow_models(), get_drought_streamflow_models()),
  .f = get_model_name_and_function
) |>
  list_rbind()


zero_CO2_streamflow_models <- set_a3_zero_model_parameters |>
  select(gauge, streamflow_model) |>
  distinct() |>
  left_join(
    streamflow_name_and_model_for_joining,
    by = join_by(streamflow_model)
  ) |>
  pull(model_function)


### Parameter sets
zero_CO2_parameter_sets <- set_a3_zero_model_parameters |>
  select(gauge, parameter, altered_parameter)

# It might be worth generalising this...
tibble_column_to_parameter_vector <- function(gauge, data) {
  data |>
    filter(gauge == {{ gauge }}) |>
    pull(altered_parameter)
}

altered_parameter_sets <- map(
  .x = gauges,
  .f = tibble_column_to_parameter_vector,
  data = set_a3_zero_model_parameters
)


### Put it together now
generate_altered_streamflow <- function(catchment_data, parameter_set, streamflow_model) {
  streamflow_model <- noquote(streamflow_model)
  streamflow_model(catchment_data, parameter_set)
}


altered_streamflow <- pmap(
  .l = list(catchment_data, altered_parameter_sets, zero_CO2_streamflow_models),
  .f = generate_altered_streamflow
) |>
  list_rbind()


###  Summary table #############################################################
no_CO2_data_joining <- set_a3_zero_model_parameters |>
  select(gauge, streamflow_model) |>
  distinct()

only_CO2_streamflow_data <- streamflow_results |>
  semi_join(
    no_CO2_data_joining,
    by = join_by(gauge, streamflow_model)
  ) |>
  select(
    gauge,
    year,
    precipitation,
    transformed_modelled_streamflow,
    realspace_modelled_streamflow,
    transformed_observed_streamflow,
    realspace_observed_streamflow
  )


gauge_join <- data |>
  select(gauge, year, p_mm) |>
  rename(
    precipitation = p_mm
  )

streamflow_data_a3_off <- altered_streamflow |>
  left_join(
    gauge_join,
    by = join_by(year, precipitation)
  ) |>
  relocate(
    gauge,
    .before = 1
  ) |>
  # add the original streamflow on
  rename(
    a3_off_modelled_transformed_streamflow = streamflow_results
  ) |>
  left_join(
    only_CO2_streamflow_data,
    by = join_by(gauge, year, precipitation)
  )


## Convert a3_off_transformed_modelled streamflow into realspace ===============
### Pull a part streamflow_data_a3_off and put it back together

convert_a3_off_transformed_to_realspace <- function(gauge, streamflow_data_a3_off, best_model_per_gauge) {
  filtered_streamflow_data_a3_off <- streamflow_data_a3_off |>
    filter(gauge == {{ gauge }})

  a3_off_transformed_streamflow <- filtered_streamflow_data_a3_off |>
    pull(a3_off_modelled_transformed_streamflow)

  parameters <- best_model_per_gauge |>
    filter(gauge == {{ gauge }}) |>
    pull(parameter_value)

  b <- parameters[length(parameters) - 1]

  a3_off_realspace_streamflow <- inverse_log_sinh_transform(
    b = b,
    a3_off_transformed_streamflow,
    offset = 0
  )

  # add it back to streamflow_data_a3_off
  filtered_streamflow_data_a3_off <- filtered_streamflow_data_a3_off |>
    add_column(
      "a3_off_modelled_realspace_streamflow" = a3_off_realspace_streamflow,
      .before = 9
    ) |>
    # realspace streamflow must be greater than zero
    mutate(
      a3_off_modelled_realspace_streamflow = if_else(a3_off_modelled_realspace_streamflow < 0, 0, a3_off_modelled_realspace_streamflow)
    )


  return(filtered_streamflow_data_a3_off)
}


streamflow_data_with_a3_off <- map(
  .x = streamflow_data_a3_off |> pull(gauge) |> unique(),
  .f = convert_a3_off_transformed_to_realspace,
  streamflow_data_a3_off = streamflow_data_a3_off,
  best_model_per_gauge = only_CO2_model_parameters
) |>
  list_rbind()


# Compare percentage difference in a3 on vs. a3 off and save -------------------
a3_on_off_difference_data <- streamflow_data_with_a3_off |>
  select(
    gauge,
    year,
    precipitation,
    a3_off_modelled_realspace_streamflow,
    realspace_modelled_streamflow,
    realspace_observed_streamflow,
    a3_off_modelled_transformed_streamflow,
    transformed_modelled_streamflow,
    transformed_observed_streamflow
  ) |>
  rename(
    realspace_a3_off_streamflow = a3_off_modelled_realspace_streamflow,
    realspace_a3_on_streamflow = realspace_modelled_streamflow,
    transformed_a3_off_streamflow = a3_off_modelled_transformed_streamflow,
    transformed_a3_on_streamflow = transformed_modelled_streamflow
  )


# Add best non-CO2 streamflow to modelled_streamflow_with_CO2_on_off -----------
only_non_CO2_models <- best_CO2_non_CO2_per_gauge |>
  filter(!str_detect(streamflow_model, "CO2")) |>
  select(gauge, streamflow_model) |>
  distinct()


non_CO2_model_streamflow <- streamflow_results |>
  semi_join(
    only_non_CO2_models,
    by = join_by(gauge, streamflow_model)
  ) |>
  select(
    gauge,
    year,
    precipitation,
    realspace_modelled_streamflow,
    transformed_modelled_streamflow
  )


### USE THIS FOR ALL CALCULATIONS ##############################################
master_streamflow_table <- a3_on_off_difference_data |>
  left_join(
    non_CO2_model_streamflow,
    by = join_by(gauge, year, precipitation)
  ) |>
  # rename and rearrange columns
  rename(
    realspace_streamflow_CO2_model_on = realspace_a3_on_streamflow,
    realspace_streamflow_CO2_model_off = realspace_a3_off_streamflow,
    realspace_streamflow_no_CO2_model = realspace_modelled_streamflow,
    transformed_streamflow_CO2_model_on = transformed_a3_on_streamflow,
    transformed_streamflow_CO2_model_off = transformed_a3_off_streamflow,
    transformed_streamflow_no_CO2_model = transformed_modelled_streamflow,
  ) |>
  relocate(
    gauge, year, precipitation, realspace_observed_streamflow,
    realspace_streamflow_no_CO2_model, realspace_streamflow_CO2_model_on,
    realspace_streamflow_CO2_model_off, transformed_observed_streamflow,
    transformed_streamflow_no_CO2_model, transformed_streamflow_CO2_model_on,
    transformed_streamflow_CO2_model_off
  )


write_csv(
  x = master_streamflow_table,
  file = "Modelling/Results/master_streamflow_table.csv"
)


# Calculate NSE for all models -------------------------------------------------
best_model_per_gauge <- best_CO2_non_CO2_per_gauge |>
  slice_min(AIC, by = gauge) |>
  select(gauge, streamflow_model) |>
  distinct() |>
  add_column(best_model = TRUE)

all_NSE_results <- streamflow_results |>
  summarise(
    NSE_value = nash_sutcliffe_efficiency(
      observed = realspace_observed_streamflow,
      modelled = realspace_modelled_streamflow
    ),
    .by = c(gauge, streamflow_model)
  ) |>
  arrange(gauge, desc(NSE_value)) |>
  left_join(
    best_model_per_gauge,
    by = join_by(gauge, streamflow_model)
  ) |>
  mutate(
    best_model = if_else(is.na(best_model), FALSE, best_model)
  ) |>
  arrange(best_model)


all_NSE_results |>
  ggplot(aes(x = NSE_value)) +
  geom_histogram(colour = "black", fill = "grey") +
  labs(
    x = "NSE",
    y = "Frequency"
  ) +
  theme_bw()


# plot NSE vs. number of parameters in a model
# convert char into function match.fun --> then length(streamflow_model$parameters)
# match.fun(NSE_parameter_number[1]) doesn't work with vector
# this is a map angle
get_parameter_number <- function(streamflow_model_name) {
  streamflow_function <- match.fun(streamflow_model_name)
  return(length(streamflow_function()$parameter))
}

NSE_parameter_number <- all_NSE_results |>
  # cannot use match.fun with vector. Use rowwise for row-by-row operation
  rowwise() |>
  mutate(
    parameter_number = get_parameter_number(streamflow_model_name = streamflow_model)
  ) |>
  ungroup()


NSE_parameter_number |>
  ggplot(aes(x = NSE_value)) +
  geom_histogram(colour = "black", fill = "grey") +
  labs(
    x = "NSE",
    y = "Frequency"
  ) +
  theme_bw() +
  facet_wrap(~ as.factor(parameter_number))


NSE_ecdf <- NSE_parameter_number |>
  ggplot(
    aes(
      x = NSE_value,
      colour = as.factor(parameter_number)
    )
  ) +
  stat_ecdf(geom = "step", pad = TRUE) +
  labs(
    x = "Nash Sutcliffe Efficiency",
    y = "Cumulative Probability",
    colour = "Number of Parameters"
  ) +
  scale_colour_brewer(palette = "Dark2") +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.2, 0.85),
    legend.background = element_rect(colour = "black", linewidth = 0.2),
    legend.title = element_text(hjust = 0.5),
    legend.key.spacing.x = unit(0.5, units = "cm"),
  ) +
  guides(
    colour = guide_legend(
      override.aes = aes(linewidth = 1.5),
      ncol = 2
    )
  )

ggsave(
  filename = "Figures/Supplementary/NSE_ecdf.pdf",
  device = "pdf",
  plot = NSE_ecdf,
  height = 130,
  width = 140,
  units = "mm"
)


# Compare NSE values between CO2 and non-CO2 models ----------------------------
# Does the best model using AIC match NSE?
best_model_per_gauge <- best_CO2_non_CO2_per_gauge |>
  slice_min(AIC, by = gauge) |>
  select(gauge, streamflow_model) |>
  distinct() |>
  mutate(
    AIC_CO2_is_best = str_detect(streamflow_model, "CO2")
  ) |>
  select(gauge, AIC_CO2_is_best)

NSE_value_comparison <- master_streamflow_table |>
  summarise(
    NSE_CO2_model = nash_sutcliffe_efficiency(
      observed = realspace_observed_streamflow,
      modelled = realspace_streamflow_CO2_model_on
    ),
    NSE_non_CO2_model = nash_sutcliffe_efficiency(
      observed = realspace_observed_streamflow,
      modelled = realspace_streamflow_no_CO2_model
    ),
    .by = gauge
  ) |>
  left_join(
    best_model_per_gauge,
    by = join_by(gauge)
  ) |>
  mutate(
    NSE_CO2_is_best = NSE_CO2_model > NSE_non_CO2_model
  )

# For catchments where AIC says CO2 model is best NSE is always higher
# For catchment where AIC says non-CO2 model is best NSE isn't always higher.
# They are after very similar. Again, the streamflow values
# are very simlar due to a5 turning off at the last couple year of the timeseries.
# Additionally, the extra parameters maybe assisting fit.


# Timeseries of master table ---------------------------------------------------
## Simplify master master_streamflow_table for only streamflow timeseries info
## Select and pivot the master table to ignore precip and transformed streamflow
timeseries_plotting_data <- master_streamflow_table |>
  select(
    gauge, year, realspace_observed_streamflow,
    realspace_streamflow_no_CO2_model, realspace_streamflow_CO2_model_on,
    realspace_streamflow_CO2_model_off
  ) |>
  pivot_longer(
    cols = starts_with("realspace"),
    names_to = "type",
    values_to = "streamflow"
  ) |>
  mutate(
    type = case_when(
      type == "realspace_streamflow_CO2_model_on" ~ "CO2 Model",
      type == "realspace_streamflow_CO2_model_off" ~ "Counterfactual",
      type == "realspace_observed_streamflow" ~ "Observed",
      type == "realspace_streamflow_no_CO2_model" ~ "non-CO2 Model",
      .default = NA
    )
  ) |>
  # set order
  mutate(
    type = factor(type, levels = c("Observed", "CO2 Model", "Counterfactual", "non-CO2 Model"))
  )


## massive plot ================================================================
mega_timeseries_plot <- timeseries_plotting_data |>
  ggplot(aes(x = year, y = streamflow, colour = type)) +
  geom_line(linewidth = 0.25) +
  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  labs(
    x = "Year",
    y = "Streamflow (mm)",
    colour = NULL
  ) +
  theme(legend.position = "bottom") +
  facet_wrap(~gauge, scales = "free_y")

ggsave(
  file = "Figures/Other/mega_timeseries_plot.pdf",
  device = "pdf",
  plot = mega_timeseries_plot,
  width = 1189,
  height = 841,
  units = "mm"
)


# For gauges where CO2 is not the best model there is high overlap with non-CO2
# model. This is because the a5 term (time of activation term) only kicks in
# in the last few years of the time series.


# Rainfall-runoff relationship of master table ---------------------------------
rainfall_runoff_plotting_data <- master_streamflow_table |>
  select(
    gauge, precipitation, transformed_observed_streamflow,
    transformed_streamflow_no_CO2_model, transformed_streamflow_CO2_model_on,
    transformed_streamflow_CO2_model_off
  ) |>
  pivot_longer(
    cols = starts_with("transformed"),
    names_to = "type",
    values_to = "streamflow"
  ) |>
  mutate(
    type = case_when(
      type == "transformed_streamflow_CO2_model_on" ~ "CO2 Model",
      type == "transformed_streamflow_CO2_model_off" ~ "Counterfactual",
      type == "transformed_observed_streamflow" ~ "Observed",
      type == "transformed_streamflow_no_CO2_model" ~ "non-CO2 Model",
      .default = NA
    )
  ) |>
  # set order
  mutate(
    type = factor(type, levels = c("Observed", "CO2 Model", "Counterfactual", "non-CO2 Model"))
  )


## massive plot ================================================================
mega_rainfall_runoff_plot <- rainfall_runoff_plotting_data |>
  ggplot(aes(x = precipitation, y = streamflow, colour = type)) +
  geom_point(size = 0.25) +
  geom_smooth(
    formula = y ~ x,
    method = lm,
    se = FALSE,
    linewidth = 0.25
  ) +
  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  labs(
    x = "Precipitation (mm)",
    y = "Transformed Streamflow",
    colour = NULL
  ) +
  theme(legend.position = "bottom") +
  facet_wrap(~gauge, scales = "free")


# why is there a larger difference in the transformed streamflow than
# realspace streamflow?
# It is doing weird things at low flows. It attempts to get to as close to
# zero as possible. When in realspace they are all small.

ggsave(
  file = "Figures/Other/mega_rainfall_runoff_plot.pdf",
  device = "pdf",
  plot = mega_rainfall_runoff_plot,
  width = 1189,
  height = 841,
  units = "mm"
)
