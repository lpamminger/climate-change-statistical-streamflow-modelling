pacman::p_load(tidyverse, gridExtra, ggExtra)


# Import functions -------------------------------------------------------------
source("./Functions/utility.R")
source("./Functions/boxcox_logsinh_transforms.R")



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
  mutate(year = as.integer(year)) |> 
  # required for log-sinh. Log-sinh current formulation has asymptote of zero. 
  # This means zero flows of ephemeral catchments cannot be transformed
  # add a really small value
  mutate(q_mm = q_mm + .Machine$double.eps^0.5) 


start_stop_indexes <- readr::read_csv(
  "./Data/Tidy/start_end_index.csv",
  show_col_types = FALSE
)


gauge_information <- readr::read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
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
  path = "./Figures/Other/",
  # used to append all plots into a single pdf
  plot = gridExtra::marrangeGrob(rainfall_runoff_plots, nrow = 1, ncol = 1),
  device = "pdf",
  units = "mm",
  width = 297,
  height = 210
)


## Visual inspection of rainfall-runoff graphs =================================
## 238231, 418025, 614044, A5090503, 226222 - autocorrelation has an impact on the slope
## 613146 - auto only an issue without a CO2 term


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
  path = "./Figures/Other",
  # used to append all plots into a single pdf
  plot = gridExtra::marrangeGrob(streamflow_time_plots, nrow = 1, ncol = 1),
  device = "pdf",
  units = "mm",
  width = 297,
  height = 210
)


## Visual inspection ===========================================================
## Thoughts go here...






## Transformation curves and histograms ========================================

## This does not work - FIX

recreate_log_sinh_using_parameter_results <- function(streamflow_model, gauge_parameter_results, xaxis_range) {
  
  log_sinh_parameters <- gauge_parameter_results |>
    filter(streamflow_model == {{ streamflow_model }}) |>
    filter(parameter == "b") |>
    pull(parameter_value)
  
  xaxis <- seq(from = xaxis_range[1], to = xaxis_range[2], by = 0.01)
  
  yaxis <- log_sinh_transform(
    b = log_sinh_parameters,
    y = xaxis,
    offset = 0
  )
  
  list(
    "xaxis" = xaxis,
    "yaxis" = yaxis,
    "streamflow_model" = streamflow_model
  ) |> 
    as_tibble()
}





transformation_curves <- function(gauge, streamflow_results, parameter_results) {
  
  ## Get modelled_realspace_streamflow and modelled_transformed_streamflow =====
  curve_data <- streamflow_results |> 
    filter(gauge == {{ gauge }}) |> 
    select(realspace_modelled_streamflow, transformed_modelled_streamflow, streamflow_model)
  
  ## Recreate curve using parameters to compare with model results =============
  xaxis_range <- curve_data |> pull(realspace_modelled_streamflow) |> range()
  
  # extract parameter values for a given streamflow model
  recreate_log_sinh_results <- map(
    .x = parameter_results |> pull(streamflow_model) |> unique(),
    .f = recreate_log_sinh_using_parameter_results,
    gauge_parameter_results = parameter_results |> filter(gauge == {{ gauge }}),
    xaxis_range = xaxis_range
  ) |> 
    list_rbind()

  
  
  ## Plot ======================================================================
  curve_data |> 
    ggplot(aes(x = realspace_modelled_streamflow, y = transformed_modelled_streamflow)) +
    geom_line(
      aes(x = xaxis, y = yaxis, colour = "red"), 
      data = recreate_log_sinh_results,
      show.legend = FALSE
      ) +
    geom_point() +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
    theme_bw() +
    facet_wrap(~streamflow_model)
  
}


transformation_curves(
  gauge = "302208",
  streamflow_results = streamflow_results,
  parameter_results = parameter_results
)

#transformation_curves <- map(
#  .x = parameter_results |> pull(gauge) |> unique(),
#  .f = transformation_curves,
#  streamflow_results = streamflow_results,
#  parameter_results = parameter_results
#)


#ggsave(
#  filename = "all_transformation_curve_plots.pdf",
#  path = "./Figures/Supplementary",
  # used to append all plots into a single pdf
#  plot = gridExtra::marrangeGrob(transformation_curves, nrow = 1, ncol = 1),
#  device = "pdf",
#  units = "mm",
#  width = 297,
#  height = 210
#)



# 2. Are the fitted parameters "acceptable"? -----------------------------------

## Utilisation of parameters ===================================================
# What parameters are being turned off that should not be turned off? ##########
check_near_zero_parameter_values <- parameter_results |>
  mutate(
    is_parameter_turned_off = near(parameter_value, 0, tol = .Machine$double.eps^0.5)
  ) |>
  filter(is_parameter_turned_off)

check_near_zero_parameter_values |> pull(parameter) |> unique()




## Check if parameters are at the bounds =======================================
bound_issues <- parameter_results |> 
  filter(!is.na(near_bounds)) |> 
  arrange(parameter_value)


## a5 can be near the upper bound --> means CO2 kicks in during the last year
## a5 can be near the lower bound --> means CO2 kicks in during the first year
## b can be near the upper bound (equal to 1) --> means transform is linear


# Check all parameters
bound_issues |> pull(parameter) |> unique()
# Only b and a5 hit the bounds


# Check b parameter
bound_issues |> filter(parameter == "b") |> pull(parameter_value) |> range()

# No issues with parameters hitting bounds






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


# There are some catchments with the exact same AIC
check_duplicates <- best_CO2_and_non_CO2_per_catchment |> 
  select(gauge, streamflow_model) |> 
  distinct() |> 
  count(gauge) |> 
  filter(n > 2)


duplicate_gauges <- check_duplicates |> pull(gauge) |> unique()

# There are 48 gauges with the exact same AIC values
# It only occurs for CO2 models
# The parameters are almost exactly the same


duplicates_best_results <- best_CO2_and_non_CO2_per_catchment |> 
  filter(gauge %in% duplicate_gauges)
  

# If the non-CO2 model is better than the CO2 model, then
# selecting which CO2 does not matter 
# --> filter out duplicates where non-CO2 is better

duplicates_where_CO2_model_is_best <- duplicates_best_results |> 
  mutate(
    contains_CO2 = if_else(str_detect(streamflow_model, "CO2"), "CO2_model", "no_CO2_model")
  ) |> 
  select(gauge, contains_CO2, AIC) |>
  distinct() |> 
  pivot_wider(
    names_from = contains_CO2,
    values_from = AIC
  ) |> 
  mutate(
    is_CO2_model_better = CO2_model < no_CO2_model
  ) |> 
  filter(is_CO2_model_better) |> 
  pull(gauge)


duplicates_to_focus_on <- best_CO2_and_non_CO2_per_catchment |> 
  filter(gauge %in% duplicates_where_CO2_model_is_best)
  

# Since AIC cannot be used to differentiate the best model for some gauges,
# visual inspect will be used. The graphs used are: 
## - rainfall-runoff graphs (all look very similar - cannot be used to differentiate)
## - streamflow-time graphs (all look very similar - cannot be used to differentiate)
## - transform streamflow vs streamflow graphs (all look very similar - cannot be used to differentiate)


# Notes:
# - all four catchment have the exact same model
# - I am inclined to select the slope variant because this is easier to detect


# Adjust best_CO2_and_non_CO2_per_catchment based on findings above:
## Take the best CO2 slope model. This will not matter for cases where the 
## non-CO2 model is better. Justification = slope is easier to detect


key_removed_intercept_model_duplicates <- duplicates_best_results |> 
  select(gauge, replicate, contains_CO2, streamflow_model) |> 
  distinct() |> 
  mutate(
    contains_intercept = if_else(str_detect(streamflow_model, "intercept"), TRUE, FALSE)
  ) |> 
  mutate(
    contains_intercept_and_CO2 = contains_CO2 & contains_intercept
  ) |> 
  # remove intercept CO2 streamflow models
  filter(!contains_intercept_and_CO2) |> 
  select(gauge, replicate, streamflow_model)


removed_intercept_model_duplicates <- duplicates_best_results |> 
  semi_join(
    key_removed_intercept_model_duplicates,
    by = join_by(gauge, replicate, streamflow_model)
  )



adjusted_best_CO2_and_non_CO2_per_catchment <- best_CO2_and_non_CO2_per_catchment |> 
  # filter out duplicate gauges
  filter(!gauge %in% duplicate_gauges) |> 
  # add back the altered duplicates (removed intercept model duplicates)
  rbind(removed_intercept_model_duplicates)



write_csv(
  adjusted_best_CO2_and_non_CO2_per_catchment,
  file = "./Modelling/Results/CMAES/best_CO2_non_CO2_per_catchment_CMAES.csv"
)


