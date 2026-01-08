# Illustration of rainfall-runoff changes

# Figures produced in this R file ----------------------------------------------
# - add something here


# CODE


# Import libraries -------------------------------------------------------------
pacman::p_load(tidyverse, sloop, patchwork)


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


best_CO2_and_non_CO2_model_and_params_per_gauge <- read_csv(
  "./Modelling/Results/CMAES/best_CO2_non_CO2_per_catchment_CMAES.csv",
  show_col_types = FALSE
)


start_stop_indexes <- read_csv(
  "Data/Tidy/start_end_index.csv",
  show_col_types = FALSE
)

# Calculate evidence ratio for filtering ---------------------------------------
high_evi_gauges <- best_CO2_and_non_CO2_model_and_params_per_gauge |>
  select(gauge, contains_CO2, AIC) |>
  distinct() |>
  pivot_wider(
    names_from = contains_CO2,
    values_from = AIC
  ) |>
  mutate(
    CO2_model = `TRUE`,
    non_CO2_model = `FALSE`,
    .keep = "unused"
  ) |>
  mutate(
    AIC_difference = CO2_model - non_CO2_model # CO2 is smaller than non-CO2 then negative and CO2 is better
  ) |>
  mutate(
    evidence_ratio = case_when(
      AIC_difference < 0 ~ exp(0.5 * abs(AIC_difference)), # when CO2 model is better
      AIC_difference > 0 ~ -exp(0.5 * abs(AIC_difference)) # when non-CO2 model is better
    )
  ) |>
  arrange(evidence_ratio) |>
  filter(evidence_ratio > 100) |>
  pull(gauge)


# Find gauge to use as an example ----------------------------------------------
## Get gauges with high record lengths
high_record_gauges <- gauge_information |>
  filter(record_length > 50) |>
  pull(gauge)

## Try using the simplest CO2 models first
best_CO2_model_gauge <- best_CO2_and_non_CO2_model_and_params_per_gauge |>
  filter(contains_CO2) |>
  filter(gauge %in% high_evi_gauges) |>
  filter(gauge %in% high_record_gauges) |>
  select(gauge, streamflow_model) |>
  distinct() |>
  arrange(desc(streamflow_model))


best_CO2_model_gauge |>
  filter(streamflow_model == "streamflow_model_slope_shifted_CO2")

# try 614044 for intercept shifted - this looks good

# 614044 = great
# 235219 = not as clear as 614044



rearrange_catchment_data_blueprint <- function(observed_data, gauge_ID, start_stop_indexes) {
  catchment_data_blueprint(
    gauge_ID = gauge_ID, 
    observed_data = observed_data, 
    start_stop_indexes = start_stop_indexes
  )
}

replace_CO2_in_data <- function(new_CO2, data) {
  data |> 
    mutate(
      CO2 = new_CO2
    )
}


illustration_plots <- function(gauge) {
  
  b <- best_CO2_and_non_CO2_model_and_params_per_gauge |>
    filter(contains_CO2) |>
    filter(gauge == {{ gauge }}) |>
    filter(parameter == "b") |>
    pull(parameter_value)
  
  
  modified_data <- data |>
    filter(gauge == {{ gauge }}) |>
    mutate(
      q_log_sinh = log_sinh_transform(b = b, y = q_mm, offset = 0)
    )
  
  
  # Plotting lines for the illustration
  CO2_for_selected_years <- modified_data |>
    filter(year %in% c(1960, 1990, 2010)) |>
    pull(CO2)
  
  
  different_CO2_modified_data <- map(
    .x = CO2_for_selected_years,
    .f = replace_CO2_in_data,
    data = modified_data
  )
  
  
  different_CO2_catchment_data <- map(
    .x = different_CO2_modified_data,
    .f = rearrange_catchment_data_blueprint,
    gauge = gauge,
    start_stop_indexes = start_stop_indexes
  )
  
  
  parameter_set <- best_CO2_and_non_CO2_model_and_params_per_gauge |>
    filter(contains_CO2) |>
    filter(gauge == {{ gauge }}) |> 
    # make sure we get a straight line by turning off a2, a4 etc.
    mutate(
      parameter_value = case_when(
        parameter == "a2" ~ 0,
        parameter == "a4" ~ 0,
        .default = parameter_value
      )
    ) |> 
    pull(parameter_value)
  
  
  streamflow_model <- best_CO2_model_gauge |> 
    filter(gauge == {{ gauge }}) |> 
    pull(streamflow_model) |> 
    match.fun()
  
  
  
  different_CO2_streamflow <- map(
    .x = different_CO2_catchment_data,
    .f = streamflow_model,
    parameter = parameter_set
  ) 
  
  
  streamflow_1 <- different_CO2_streamflow[[1]] |> 
    select(year, precipitation, streamflow_results) |> 
    rename(streamflow_1960 = streamflow_results)
  
  streamflow_2 <- different_CO2_streamflow[[2]] |> 
    select(year, precipitation, streamflow_results) |> 
    rename(streamflow_1990 = streamflow_results)
  
  streamflow_3 <- different_CO2_streamflow[[3]] |> 
    select(year, precipitation, streamflow_results) |> 
    rename(streamflow_2010 = streamflow_results)
  
  rainfall_runoff_relationship <- streamflow_1 |> 
    left_join(
      streamflow_2,
      by = join_by(year, precipitation)
    ) |> 
    left_join(
      streamflow_3,
      by = join_by(year, precipitation)
    ) |> 
    pivot_longer(
      cols = starts_with("streamflow"),
      names_to = "rainfall_runoff_year",
      values_to = "streamflow"
    )
  
  
  rainfall_runoff_relationship |> 
    ggplot(aes(x = precipitation, y = streamflow, colour = rainfall_runoff_year)) +
    geom_line()
  
  
  # put it together now
  plot <- modified_data |>
    ggplot(aes(x = p_mm, y = q_log_sinh, fill = year)) +
    geom_line(
      data = rainfall_runoff_relationship,
      aes(x = precipitation, y = streamflow, colour = rainfall_runoff_year),
      inherit.aes = FALSE
    ) +
    geom_point(
      colour = "black",
      shape = 21,
      size = 3
    ) +
    labs(
      x = "Annual Precipitation (mm)",
      y = "Log-sinh Annual Streamflow (mm)",
      title = gauge
    ) +
    scale_colour_manual(values = c("#440154FF", "#33638DFF", "#B8DE29FF")) +
    scale_fill_continuous(palette = "viridis") + # options: viridis, plasma
    theme_bw()
  
  return(plot)
}

# plots
illustration <- illustration_plots(gauge = "614044") | illustration_plots(gauge = "238235")
final_illustration <- illustration + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(
  filename = "illustration.pdf",
  plot = final_illustration,
  device = "pdf",
  path = "Figures/Other",
  width = 297,
  height = 180,
  units = "mm"
)
