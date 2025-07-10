pacman::p_load(tidyverse, gridExtra)


# Import functions -------------------------------------------------------------
source("./Functions/utility.R")


# Import the calibrated .csv's -------------------------------------------------
parameter_results <- read_csv(
  "./Modelling/Results/CMAES/cmaes_parameter_results.csv",
  show_col_types = FALSE
)


streamflow_results <- read_csv(
  "./Modelling/Results/CMAES/cmaes_streamflow_results.csv",
  show_col_types = FALSE
)


data <- readr::read_csv(
  "./Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  show_col_types = FALSE
) |>
  mutate(year = as.integer(year))


start_stop_indexes <- readr::read_csv(
  "./Data/Tidy/start_end_index.csv",
  show_col_types = FALSE
)


# 1. Are there any major issues with CMAES fitting? - visual inspection --------
## streamflow results only include data calibrated on
## page-per-gauge with unique models



## rainfall-runoff comparison ==================================================
rainfall_runoff_plot <- function(gauge, streamflow_results, parameter_results) {
  gauge_streamflow_results <- streamflow_results |>
    filter(gauge == {{ gauge }}) |>
    pivot_longer(
      cols = starts_with("transformed"),
      names_to = "observed_or_modelled",
      values_to = "streamflow"
    )

  # Get AIC information
  AIC_for_graphs <- parameter_results |>
    filter(gauge == {{ gauge }}) |>
    distinct(streamflow_model, .keep_all = TRUE) |>
    select(streamflow_model, AIC) |>
    add_column(
      x_pos_rainfall_runoff = gauge_streamflow_results |> pull(precipitation) |> min(), # minimum observed rainfall
      y_pos = gauge_streamflow_results |> pull(realspace_observed_streamflow) |> max(), # maximum observed streamflow
      x_pos_streamflow_time = 1959 # first value
    ) |>
    mutate(
      label = paste0("AIC = ", round(AIC, digits = 2))
    )


  # Make plot
  gauge_streamflow_results |>
    ggplot(aes(x = precipitation, y = streamflow, colour = observed_or_modelled)) +
    geom_smooth(
      method = lm,
      formula = y ~ x,
      se = FALSE
    ) +
    geom_point(size = 0.5) +
    # The aic y position is not correct - fix
    geom_label(
      aes(x = x_pos_rainfall_runoff, y = y_pos, label = label),
      data = AIC_for_graphs,
      inherit.aes = FALSE,
      size = 3,
      vjust = 1,
      hjust = -0.25
    ) +
    scale_color_brewer(palette = "Set1") +
    theme_bw() +
    labs(
      x = "Precipitation (mm)",
      y = "Log-sinh streamflow (mm)",
      colour = NULL,
      title = paste0("Gauge: ", gauge)
    ) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 5),
      plot.title = element_text(hjust = 0.5)
    ) +
    facet_wrap(~streamflow_model, scales = "free_y")
}


rainfall_runoff_plots <- map(
  .x = parameter_results |> pull(gauge) |> unique(),
  .f = rainfall_runoff_plot,
  streamflow_results = streamflow_results,
  parameter_results = parameter_results
)


ggsave(
  filename = "all_rainfall_runoff_plots.pdf",
  path = "./Figures/Supplementary/",
  # used to append all plots into a single pdf
  plot = gridExtra::marrangeGrob(rainfall_runoff_plots, nrow = 1, ncol = 1),
  device = "pdf",
  units = "mm",
  width = 297,
  height = 210
)


## Visual inspection of rainfall-runoff graphs =================================
# REDO this

# line-of-best fit not perfectly overlapping with observed
# - 234203 (check bounds)
# - 234209 - all wonky probably something related to the bounds
# - 609002
# - G8150098
# There are more than mentioned above


## streamflow-time comparison ==================================================
streamflow_time_plot <- function(gauge, streamflow_results, parameter_results) {
  gauge_streamflow_results <- streamflow_results |>
    filter(gauge == {{ gauge }}) |>
    pivot_longer(
      cols = starts_with("realspace"),
      names_to = "observed_or_modelled",
      values_to = "streamflow"
    )

  # Get AIC information
  AIC_for_graphs <- parameter_results |>
    filter(gauge == {{ gauge }}) |>
    distinct(streamflow_model, .keep_all = TRUE) |>
    select(streamflow_model, AIC) |>
    add_column(
      x_pos_rainfall_runoff = gauge_streamflow_results |> pull(precipitation) |> min(), # minimum observed rainfall
      y_pos = gauge_streamflow_results |> filter(observed_or_modelled == "realspace_observed_streamflow") |> pull(streamflow) |> max(), # maximum observed streamflow
      x_pos_streamflow_time = 1959 # first value
    ) |>
    mutate(
      label = paste0("AIC = ", round(AIC, digits = 2))
    )


  # Make plot
  gauge_streamflow_results |>
    ggplot(aes(x = year, y = streamflow, colour = observed_or_modelled)) +
    geom_line() +
    geom_point(size = 0.5) +
    geom_label(
      aes(x = x_pos_streamflow_time, y = y_pos, label = label),
      data = AIC_for_graphs,
      inherit.aes = FALSE,
      size = 3,
      vjust = 1,
      hjust = -0.25
    ) +
    scale_color_brewer(palette = "Set1") +
    theme_bw() +
    labs(
      x = "Year",
      y = "Streamflow (mm)",
      colour = NULL,
      title = paste0("Gauge: ", gauge)
    ) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 5),
      plot.title = element_text(hjust = 0.5)
    ) +
    facet_wrap(~streamflow_model, scales = "free_y")
}



streamflow_time_plots <- map(
  .x = parameter_results |> pull(gauge) |> unique(),
  .f = streamflow_time_plot,
  streamflow_results = streamflow_results,
  parameter_results = parameter_results
)


ggsave(
  filename = "all_streamflow_time_plots.pdf",
  path = "./Figures/Supplementary",
  # used to append all plots into a single pdf
  plot = gridExtra::marrangeGrob(streamflow_time_plots, nrow = 1, ncol = 1),
  device = "pdf",
  units = "mm",
  width = 297,
  height = 210
)


## Visual inspection ===========================================================
## Nothing seems obviously wrong


# Maybe transforms here?


# 2. Are the fitted parameters "acceptable"? -----------------------------------

## Utilisation of parameters ===================================================
# What parameters are being turned off that should not be turned off? ##########
check_near_zero_parameter_values <- parameter_results |>
  mutate(
    is_parameter_turned_off = near(parameter_value, 0, tol = .Machine$double.eps^0.5)
  ) |>
  filter(is_parameter_turned_off)

check_near_zero_parameter_values |> pull(parameter) |> unique()
# Only the `a` parameter want to be turned-off

check_near_zero_parameter_values |> pull(parameter) |> length()
# The `a` parameter wants to be turned off for x catchment-model combinations (out of 9684)


# Results from latest batch:
parameter_results |> filter(parameter == "a") |> pull(parameter_value) |> range(na.rm = T)

# New results in progress

# - range of a values is a = -9.99E-7 to 1E-1
# - lower range of a0 needs increasing (do the same for drought terms) maybe x3 max_Q
# - leave slope terms
# - sd and b parameters are good
# - a parameter is not good. Because set a to lower 1E-8 and upper larger?
# --> I think a is a log-transform angle set the lower to 1E-8 since most
#     catchments do not want to be shifted
# change a so lower is 1E-8 and upper is calculated - can use log-trans
# change intercept to x3


# Results
# - a3_intercept needs to be larger 
# - a3_slope and a1 are hitting the bounds for some catchment
# - 327 `a` values have hit upper bound and 204 near lower - for 9684 its fine
# What do I do with the slope terms? make them bigger? Leave them?
# Talk with Murray
# do i even need an `a` variable?
# Do a = 1E-8 or a = 0.1 correspond to catchments with zero flow? Would an offset = 1 fix this?
# A lot of them are - maybe try an offset?

# QJ cautions against large values as it turn into linear with no normalisation stabilitation
# Recommended range is -15 to 0



## Check if parameters are at the bounds =======================================
bound_issues <- parameter_results |> 
  filter(!is.na(near_bounds)) |> 
  arrange(parameter_value)


## Only issues for a and a5 = good
## a can be near zero and a5 can be near the upper/lower bound
bound_issues |> pull(parameter) |> unique()

# Look like everything hit the bounds 

## Mannually check if values are okay
manual_check_bounds <- parameter_results |>
  filter(parameter == "a") |> # change parameter here
  arrange(desc(parameter_value)) |> 
  filter(!is.na(parameter_value))

head(manual_check_bounds)
tail(manual_check_bounds)

# Are bounds being hit on zero flow catchments
gauges_with_a_hitting_bounds <- manual_check_bounds |> 
  filter(!is.na(near_bounds)) |> 
  pull(gauge) |> 
  unique()

y <- manual_check_bounds |> 
  filter(!is.na(near_bounds)) 


x <- data |> 
  filter(gauge %in% gauges_with_a_hitting_bounds) |> 
  summarise(
    min_flow = min(q_mm, na.rm = TRUE),
    .by = gauge
  ) |> 
  left_join(
    y,
    by = join_by(gauge)
  ) |> 
  distinct() |> 
  select(gauge, min_flow, streamflow_model, parameter_value) |> 
  arrange(min_flow)


# Be methodical:
# - determine what catchments are hitting the bounds 
# - determine what parameters are hitting the bounds
# - change one parameter at a time. Test in vignette
# - my initial reaction is don't change slope terms to greater than 1
# - play with intercept terms first (i.e, a0, a0_d, a0_n, a4, a3_intercept)
# - start with streamflow_model_precip_only 


precip_only_check <- parameter_results |> 
  filter(streamflow_model == "streamflow_model_precip_only") |> 
  filter(parameter %in% c("a0", "a1")) |> 
  filter(!is.na(near_bounds))


negative_transformed_streamflow <- streamflow_results |> 
  filter(transformed_modelled_streamflow < 0) |> 
  arrange(transformed_modelled_streamflow)


## Do the wobbly catchments have strange bounds? ===============================
check_parameters_wobbly_catchments <- parameter_results |> 
  filter(gauge %in% wobbly_catchments) |> 
  filter(parameter == "sd") |> 
  arrange(parameter_value)

head(check_parameters_wobbly_catchments)
tail(check_parameters_wobbly_catchments)

# The near bounds are only for a and a5. This is okay
# I do not see anything inherently wrong with parameters - they are not
# near any of the bounds





# 3. Save best CO2 and non-CO2 results per gauge -------------------------------
best_CO2_and_non_CO2_per_catchment <- parameter_results |>
  select(!c(optimiser, loglikelihood, exit_message, near_bounds, objective_function)) |>
  mutate(
    contains_CO2 = str_detect(streamflow_model, "CO2"),
    .after = 2
  ) |>
  slice_min(
    AIC,
    by = c(gauge, contains_CO2)
  )

write_csv(
  best_CO2_and_non_CO2_per_catchment,
  file = "./Modelling/Results/CMAES/best_CO2_non_CO2_per_catchment_CMAES.csv"
)


