# Tim's suggestions:
# 1. Pick out 5 sample sites where CO2 is the best in detail. Pick out based on t-series and CO2-nonCO2 residual.
# 2. For each site:
##  - The best CO2, best non-CO2 and observed streamflow.
##  - The observed and CO2 and observed non-CO2 residuals
##  - The time series of CO2 minus non-CO2 streamflow
##  - Acf or partial acf of observed residuals for CO2 (a trend indicates we are missing something – possible fix may be introducing a lag-2 term or try equation on next slide)



# Clear console ----------------------------------------------------------------
cat("\014")

# Import libraries -------------------------------------------------------------
pacman::p_load(tidyverse, pdftools)


# Import functions -------------------------------------------------------------
source("./Functions/utility.R")
source("./Functions/boxcox_transforms.R")


# Import data ------------------------------------------------------------------
CMAES_results <- read_csv("./Results/CMAES_results/CMAES_parameter_results.csv",
  show_col_types = FALSE
)

streamflow_results <- read_csv("Results/CMAES_results/CMAES_streamflow_results.csv",
  show_col_types = FALSE,
  col_select = !optimiser
)

gauge_information <- read_csv("./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE,
  col_select = c(gauge, bc_lambda, state, lat, lon)
) # CAMELS is Australia wide.



# Selecting sites --------------------------------------------------------------
selected_sites <- c("614044", "113004A", "407220", "230210", "302214")

parameter_utilisation_selected_sites <- CMAES_results |>
  filter(gauge %in% selected_sites) |>
  select(!c(optimiser, loglikelihood, exit_message, near_bounds)) |>
  distinct() |>
  unite(
    col = streamflow_model_objective_function,
    c(streamflow_model, objective_function),
    sep = "_",
    remove = FALSE,
    na.rm = FALSE
  ) |>
  mutate(
    contains_CO2 = str_detect(streamflow_model_objective_function, "CO2"),
    contains_CO2 = if_else(contains_CO2, "CO2", "no_CO2"),
    .after = 2
  ) |>
  slice_min(
    AIC,
    by = c(gauge, contains_CO2)
  ) |>
  select(!streamflow_model_objective_function) |>
  filter(contains_CO2 == "no_CO2") #|> # only interested in the CO2 models
  filter(abs(parameter_value) < 1E-4)

# All parameters are non-zero


# Produce a timeseries plot for each site --------------------------------------
## Compare the observed, best CO2 and best non-CO2 (5 x 1 graph)

best_CO2_and_non_CO2_per_catchment <- CMAES_results |>
  filter(gauge %in% selected_sites) |>
  select(!c(parameter, parameter_value, optimiser, loglikelihood, exit_message, near_bounds)) |>
  distinct() |>
  unite(
    col = streamflow_model_objective_function,
    c(streamflow_model, objective_function),
    sep = "_",
    remove = FALSE,
    na.rm = FALSE
  ) |>
  mutate(
    contains_CO2 = str_detect(streamflow_model_objective_function, "CO2"),
    contains_CO2 = if_else(contains_CO2, "CO2", "no_CO2"),
    .after = 2
  ) |>
  slice_min(
    AIC,
    by = c(gauge, contains_CO2)
  )


filter_streamflow_results <- streamflow_results |>
  semi_join(
    best_CO2_and_non_CO2_per_catchment,
    by = join_by(gauge, streamflow_model, objective_function)
  )


tidy_boxcox_streamflow <- filter_streamflow_results |>
  pivot_longer(
    cols = c(observed_boxcox_streamflow, modelled_boxcox_streamflow),
    names_to = "name",
    values_to = "boxcox_streamflow"
  ) |>
  unite(
    col = streamflow_model_objective_function,
    c(streamflow_model, objective_function),
    sep = "_",
    remove = TRUE,
    na.rm = FALSE
  ) |>
  mutate(
    streamflow_model_objective_function = if_else(name == "observed_boxcox_streamflow", "observed", streamflow_model_objective_function)
  ) |>
  select(!name) |>
  mutate(
    streamflow_type = case_when(
      str_detect(streamflow_model_objective_function, "CO2") & !str_detect(streamflow_model_objective_function, "observed") ~ "CO2",
      !str_detect(streamflow_model_objective_function, "CO2") & !str_detect(streamflow_model_objective_function, "observed") ~ "non_CO2",
      .default = "observed"
    )
  ) |>
  select(!c(streamflow_model_objective_function))


tidy_streamflow <- tidy_boxcox_streamflow |>
  left_join(
    gauge_information,
    by = join_by(gauge)
  ) |>
  select(!c(state, lat, lon)) |>
  mutate(
    streamflow = boxcox_inverse_transform(yt = boxcox_streamflow, lambda = bc_lambda),
    .by = gauge
  )


plot_streamflow_timeseries <- tidy_streamflow |>
  ggplot(aes(x = year, y = streamflow, colour = streamflow_type)) +
  geom_line(na.rm = TRUE, alpha = 0.9) +
  theme_bw() +
  scale_colour_brewer(palette = "Dark2") +
  labs(
    x = "Year",
    y = "Streamflow (mm)",
    title = "Streamflow timeseries"
  ) +
  facet_wrap(~gauge, scales = "free_y", nrow = 5) +
  theme(legend.title = element_blank())


ggsave(
  filename = paste0("selected_sites_streamflow_timeseries_comparison_", get_date(), ".pdf"),
  plot = plot_streamflow_timeseries,
  device = "pdf",
  path = "./Graphs/examination_selected_sites",
  width = 297,
  height = 210,
  units = "mm"
)


# Plot the observed - CO2 and observed - nonCO2 residuals ----------------------
difference_to_observed_streamflow <- tidy_streamflow |>
  select(!c(bc_lambda, boxcox_streamflow)) |>
  distinct() |>
  pivot_wider(
    names_from = streamflow_type,
    values_from = streamflow,
  ) |>
  mutate(
    CO2_minus_non_CO2 = CO2 - non_CO2,
    observed_minus_CO2 = observed - CO2,
    observed_minus_non_CO2 = observed - non_CO2
  )


plot_difference_observed_residuals <- difference_to_observed_streamflow |>
  pivot_longer(
    cols = c(observed_minus_CO2, observed_minus_non_CO2),
    names_to = "residual_type",
    values_to = "residual_value"
  ) |>
  ggplot(aes(x = year, y = residual_value, colour = residual_type)) +
  geom_line(na.rm = TRUE) +
  theme_bw() +
  labs(
    x = "Year",
    y = "Observed streamflow minus model streamflow (mm)",
    colour = "Residual Type",
    title = "Observed minus modelled streamflow residuals"
  ) +
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~gauge, scales = "free_y", nrow = 5)


ggsave(
  filename = paste0("selected_sites_observed_residual_", get_date(), ".pdf"),
  plot = plot_difference_observed_residuals,
  device = "pdf",
  path = "./Graphs/examination_selected_sites",
  width = 297,
  height = 210,
  units = "mm"
)


# difference_CO2_minus_non_CO2 timeseries --------------------------------------
plot_difference_CO2_minus_non_CO2 <- difference_to_observed_streamflow |>
  ggplot(aes(x = year, y = CO2_minus_non_CO2)) +
  geom_line(na.rm = TRUE) +
  theme_bw() +
  labs(
    x = "Year",
    y = "Streamflow from best CO2 model minus streamflow from best non-CO2 model (mm)",
    title = "CO2 model streamflow minus non-CO2 model streamflow"
  ) +
  facet_wrap(~gauge, scales = "free_y", nrow = 5)


ggsave(
  filename = paste0("selected_sites_difference_CO2_minus_non_CO2_", get_date(), ".pdf"),
  plot = plot_difference_CO2_minus_non_CO2,
  device = "pdf",
  path = "./Graphs/examination_selected_sites",
  width = 297,
  height = 210,
  units = "mm"
)


# acf or partial acf plots of the observed residuals ---------------------------
## only interested in observed - CO2 plots


## acf using ggplot ============================================================
## https://stackoverflow.com/questions/17788859/acf-plot-with-ggplot2-setting-width-of-geom-bar

confidence_interval_acf <- function(acf_object, alpha = 0.05) {
  return(qnorm((1 + (1 - alpha)) / 2) / sqrt(acf_object$n.used))
}


tibble_confidence_interval_acf <- function(acf_object, name = NULL, alpha = 0.05) {
  upper_bound <- confidence_interval_acf(acf_object, alpha)
  lower_bound <- -upper_bound
  result_tibble <- as_tibble(cbind(name, upper_bound, lower_bound))
  result_tibble$upper_bound <- as.numeric(upper_bound)
  result_tibble$lower_bound <- as.numeric(lower_bound)
  return(result_tibble)
}


list_of_observed_minus_CO2_per_catchment <- difference_to_observed_streamflow |>
  select(year, gauge, observed_minus_CO2) |>
  pivot_wider(
    names_from = gauge,
    values_from = observed_minus_CO2
  ) |>
  select(!year) |>
  unclass()


remove_na_vector <- function(vector) {
  vector[!is.na(vector)]
}


no_NA_list_of_observed_minus_CO2_per_catchment <- map(
  .x = list_of_observed_minus_CO2_per_catchment,
  .f = remove_na_vector
)


acf_objects_per_gauge <- map(
  .x = no_NA_list_of_observed_minus_CO2_per_catchment,
  .f = acf,
  plot = FALSE
)


get_lags_from_acf <- function(acf_object, name = NULL) {
  acf_tibble <- with(acf_object, tibble(lag, acf))
  acf_tibble |>
    add_column(
      gauge = {{ name }},
      .before = 1
    )
}


lags_from_acf <- map2(
  .x = acf_objects_per_gauge,
  .y = names(acf_objects_per_gauge),
  .f = get_lags_from_acf
) |>
  list_rbind()


confidence_intervals_acf <- map2(
  .x = acf_objects_per_gauge,
  .y = names(acf_objects_per_gauge),
  .f = tibble_confidence_interval_acf
) |>
  list_rbind() |>
  rename(
    gauge = name
  )

complete_acf_per_gauge <- lags_from_acf |>
  left_join(
    confidence_intervals_acf,
    by = join_by(gauge)
  )




acf_ggplot <- complete_acf_per_gauge |>
  ggplot(aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0), linewidth = 1) +
  geom_segment(mapping = aes(xend = lag, yend = 0), linewidth = 1) +
  geom_hline(
    aes(yintercept = lower_bound),
    linetype = 2,
    linewidth = 1,
    colour = "blue"
  ) +
  geom_hline(
    aes(yintercept = upper_bound),
    linetype = 2,
    linewidth = 1,
    colour = "blue"
  ) +
  labs(
    x = "Lag",
    y = "ACF",
    title = "ACF graphs"
  ) +
  theme_bw() +
  facet_wrap(~gauge)


ggsave(
  filename = paste0("selected_sites_acf_plots_", get_date(), ".pdf"),
  plot = acf_ggplot,
  device = "pdf",
  path = "./Graphs/examination_selected_sites",
  width = 297,
  height = 210,
  units = "mm"
)

# Combine all files and save
pdf_paths_to_combine <- list.files(
  path = "./Graphs/examination_selected_sites",
  recursive = FALSE, # I don't want it looking in other folders
  full.names = TRUE
)

pdf_combine(input = pdf_paths_to_combine, output = "./Graphs/examination_selected_sites/five_site_examination.pdf")
