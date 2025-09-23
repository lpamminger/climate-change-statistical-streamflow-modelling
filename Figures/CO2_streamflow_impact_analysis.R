# CO2 streamflow impact analysis

# Figures produced in this R file ----------------------------------------------
# 1. Main --> streamflow_percentage_difference_with_timeseries.pdf
# 2. Supplementary --> streamflow_timeseries_data.pdf and captions
# 3. Supplementary --> streamflow_percentage_difference_best_CO2_model_vs_best_non_CO2_model.pdf
# 4. Other --> CO2_on_off_decade_histogram.pdf
# 5. Other --> CO2_on_off_rainfall_runoff_comparison.pdf (mega) -  only catchments where CO2 model is better than non-CO2 model (evidence ratio > 0)
# 6. Other --> streamflow_CO2_percentage_change_vs_prop_forested.pdf
# 7. Other --> CO2_model_vs_non_CO2_model_rainfall_runoff.pdf - The transformed_realspace is different depending on model used. This means 2 observed transformed streamflow is required.
# 8. Other --> complete_timeseries_plot.pdf






# CODE









# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, truncnorm, sloop, patchwork, ozmaps, sf, patchwork, metR, ggmagnify)



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

lat_lon_gauge <- gauge_information |>
  select(gauge, lat, lon, state)


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

DREAM_CO2_impact_uncertainty_on_streamflow <- read_csv(
  "Modelling/Results/DREAM/DREAM_CO2_impact_uncertainty_on_streamflow.csv",
  show_col_types = FALSE
)



# Compare best CO2 and non-CO2 models ------------------------------------------
## Summarise modelled streamflow data ==========================================
only_gauge_model_best_CO2_non_CO2_per_gauge <- best_CO2_non_CO2_per_gauge |>
  select(gauge, streamflow_model) |>
  distinct()



streamflow_data_best_CO2_non_CO2 <- streamflow_results |>
  semi_join(
    only_gauge_model_best_CO2_non_CO2_per_gauge,
    by = join_by(gauge, streamflow_model)
  )



## Are the best CO2 and non-CO2 models equivalent? ==============================
compare_equivalent_models <- only_gauge_model_best_CO2_non_CO2_per_gauge |>
  mutate(
    CO2_or_non_CO2_model = if_else(str_detect(streamflow_model, "CO2"), "CO2_model", "non_CO2_model")
  ) |>
  pivot_wider(
    id_cols = gauge,
    names_from = CO2_or_non_CO2_model,
    values_from = streamflow_model
  )



### Comparison table ### --> match equivalent models one with CO2 one without
streamflow_model_equivalent_comparison_table <- tribble(
  ~non_CO2_model,                                        ~CO2_model,
  "streamflow_model_precip_only",                        "streamflow_model_intercept_shifted_CO2",
  "streamflow_model_precip_only",                        "streamflow_model_slope_shifted_CO2",
  "streamflow_model_precip_auto",                        "streamflow_model_intercept_shifted_CO2_auto",
  "streamflow_model_precip_auto",                        "streamflow_model_slope_shifted_CO2_auto",
  "streamflow_model_precip_seasonal_ratio",              "streamflow_model_intercept_shifted_CO2_seasonal_ratio",
  "streamflow_model_precip_seasonal_ratio",              "streamflow_model_slope_shifted_CO2_seasonal_ratio",
  "streamflow_model_precip_seasonal_ratio_auto",         "streamflow_model_intercept_shifted_CO2_seasonal_ratio_auto",
  "streamflow_model_precip_seasonal_ratio_auto",         "streamflow_model_slope_shifted_CO2_seasonal_ratio_auto",
  "streamflow_model_drought_precip_only",                "streamflow_model_drought_intercept_shifted_CO2",
  "streamflow_model_drought_precip_only",                "streamflow_model_drought_slope_shifted_CO2",
  "streamflow_model_drought_precip_auto",                "streamflow_model_drought_intercept_shifted_CO2_auto",
  "streamflow_model_drought_precip_auto",                "streamflow_model_drought_slope_shifted_CO2_auto",
  "streamflow_model_drought_precip_seasonal_ratio",      "streamflow_model_drought_intercept_shifted_CO2_seasonal_ratio",
  "streamflow_model_drought_precip_seasonal_ratio",      "streamflow_model_drought_slope_shifted_CO2_seasonal_ratio",
  "streamflow_model_drought_precip_seasonal_ratio_auto", "streamflow_model_drought_intercept_shifted_CO2_seasonal_ratio_auto",
  "streamflow_model_drought_precip_seasonal_ratio_auto", "streamflow_model_drought_slope_shifted_CO2_seasonal_ratio_auto"
)



### For each row in compare_equivalent_models check if column CO2_model and non_CO2_model is in the comparison table
### Do this by using anti_join and counting row numbers
non_equivalent_models <- compare_equivalent_models |>
  anti_join(
    streamflow_model_equivalent_comparison_table,
    by = join_by(non_CO2_model, CO2_model)
  )

percent_non_equivalent_models <- (nrow(non_equivalent_models) / nrow(compare_equivalent_models)) * 100
cat("The percentage of gauges with non-equivalent models is", percent_non_equivalent_models, "%")

# Compare best CO2 with CO2 component turned off -------------------------------
## Get catchments where the CO2 model is the best ==============================
best_model_per_gauge <- best_CO2_non_CO2_per_gauge |>
  slice_min(
    AIC,
    by = gauge
  )

only_CO2_best_models <- best_model_per_gauge |>
  mutate(
    is_CO2_model = str_detect(streamflow_model, "CO2")
  ) |>
  filter(is_CO2_model) |>
  select(!is_CO2_model)


## Set a3 terms (slope and intercept to zero) ==================================
set_a3_zero_CO2_best_models <- only_CO2_best_models |>
  mutate(
    altered_parameter = if_else(str_detect(parameter, "a3"), 0, parameter_value)
  )


## Generate streamflow with altered_parameter ==================================
### Build catchment_dataset objects ############################################
CO2_gauges <- set_a3_zero_CO2_best_models |>
  pull(gauge) |>
  unique()

CO2_catchment_data <- map(
  .x = CO2_gauges,
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



zero_CO2_streamflow_models <- set_a3_zero_CO2_best_models |>
  select(gauge, streamflow_model) |>
  distinct() |>
  left_join(
    streamflow_name_and_model_for_joining,
    by = join_by(streamflow_model)
  ) |>
  pull(model_function)


### Parameter sets
zero_CO2_parameter_sets <- set_a3_zero_CO2_best_models |>
  select(gauge, parameter, altered_parameter)

# It might be worth generalising this...
tibble_column_to_parameter_vector <- function(gauge, data) {
  data |>
    filter(gauge == {{ gauge }}) |>
    pull(altered_parameter)
}

altered_parameter_sets <- map(
  .x = CO2_gauges,
  .f = tibble_column_to_parameter_vector,
  data = set_a3_zero_CO2_best_models
)


### Put it together now
generate_altered_streamflow <- function(catchment_data, parameter_set, streamflow_model) {
  streamflow_model <- noquote(streamflow_model)
  streamflow_model(catchment_data, parameter_set)
}


altered_streamflow <- pmap(
  .l = list(CO2_catchment_data, altered_parameter_sets, zero_CO2_streamflow_models),
  .f = generate_altered_streamflow
) |>
  list_rbind()

gauge_join <- data |>
  select(gauge, year, p_mm) |>
  rename(
    precipitation = p_mm
  )




###  Summary table #############################################################
no_CO2_data_joining <- set_a3_zero_CO2_best_models |>
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
# can add observed data here as well

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
    filter(contains_CO2) |>
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
  best_model_per_gauge = best_model_per_gauge
) |>
  list_rbind()




# Compare percentage difference in a3 on vs. a3 off ----------------------------
# Aim: Compare the difference in a3_on vs a3_off for a specified decade
# Method:
# 1. In realspace, sum streamflow for a given gauges with the a3 parameter on.
#    Then repeat with the a3 parameter off.
# 2. Compare the percentage difference (CO2_on - CO2_off) / number of years

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




## Percentage difference calculation ===========================================

### Method:
### - select two 10 year periods
### - sum the modelled streamflow with CO2 off (natural) and CO2 on (anthropogenic)
###   over the years during the 10 year periods
### - find the difference in the two 10 year periods
### - average by number of years during decade
### - percentage change is ((CO2_on - CO2_off) / CO2_off) * 100

decade_1 <- seq(from = 1990, to = 1999)
decade_2 <- seq(from = 2012, to = 2021)



percentage_difference_a3_on_off_data <- a3_on_off_difference_data |>
  filter(year %in% c(decade_1, decade_2)) |>
  # add decade group for summarising
  mutate(
    decade = case_when( # year - (year %% 10)
      year %in% decade_1 ~ 1,
      year %in% decade_2 ~ 2,
      .default = NA
    )
  ) |>
  filter(!is.na(decade)) |>
  # sum streamflow for each decade
  summarise(
    sum_decade_realspace_CO2_off_streamflow = sum(realspace_a3_off_streamflow),
    sum_decade_realspace_CO2_on_streamflow = sum(realspace_a3_on_streamflow),
    sum_decade_realspace_observed_streamflow = sum(realspace_observed_streamflow),
    years_of_data = n(),
    .by = c(gauge, decade)
  ) |>
  # find the absolution and percentage difference
  mutate(
    realspace_CO2_off_streamflow_per_year = sum_decade_realspace_CO2_off_streamflow / years_of_data,
    realspace_a3_on_streamflow_per_year = sum_decade_realspace_CO2_on_streamflow / years_of_data,
    CO2_impact_on_streamflow_mm_per_year = (realspace_a3_on_streamflow_per_year - realspace_CO2_off_streamflow_per_year),
    CO2_impact_on_streamflow_percent = (CO2_impact_on_streamflow_mm_per_year / realspace_CO2_off_streamflow_per_year) * 100
  ) |>
  arrange(desc(CO2_impact_on_streamflow_percent)) # Large percentage changes are not tied to years_of_data




## Percentage difference plot ==================================================

### Calculate evidence ratio for possible filtering ############################
evidence_ratio_calc <- best_CO2_non_CO2_per_gauge |>
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
  select(gauge, evidence_ratio) |>
  arrange(evidence_ratio)



plot_ready_percentage_difference_a3_on_off_data <- percentage_difference_a3_on_off_data |>
  left_join(
    lat_lon_gauge,
    by = join_by(gauge)
  ) |>
  left_join(
    evidence_ratio_calc,
    by = join_by(gauge)
  ) |>
  left_join(
    DREAM_CO2_impact_uncertainty_on_streamflow,
    by = join_by(gauge, decade)
  )


# Compare DREAM values to CMAES values
compare_CMAES_and_DREAM <- plot_ready_percentage_difference_a3_on_off_data |>
  select(gauge, decade, CO2_impact_on_streamflow_percent, IQR_CO2_impact_on_streamflow_percentage, median_CO2_impact_on_streamflow_percentage) |>
  drop_na() |>
  mutate(
    CMAES_DREAM_diff = CO2_impact_on_streamflow_percent - median_CO2_impact_on_streamflow_percentage,
    relative_CMAES_DREAM_diff = (CMAES_DREAM_diff / median_CO2_impact_on_streamflow_percentage) * 100
  ) |>
  arrange(CMAES_DREAM_diff)

compare_CMAES_and_DREAM |>
  ggplot(aes(y = relative_CMAES_DREAM_diff)) +
  geom_boxplot(
    staplewidth = 0.5
  ) +
  labs(y = "Relative percentage difference between best CMAES percentage and median DREAM percentage") +
  theme_bw()



### Calculate limits and breaks ################################################
make_limits <- function(timeseries) {
  # round up to next whole number
  limits <- timeseries |> range()
  sign_limits <- sign(limits)

  sign_limits * ceiling(abs(limits))
}


CO2_impact_on_streamflow_percent_limits <- plot_ready_percentage_difference_a3_on_off_data |>
  pull(CO2_impact_on_streamflow_percent) |>
  make_limits() |>
  as.double() |>
  round_any(accuracy = 10, f = ceiling) # round-up to nearest 10


hard_coded_breaks_CO2_impact_of_streamflow <- c(-75, -50, -25, -10, -1, 0, 1, 10, 25, 50, 75)




### Custom colour palette
big_palette <- function(x) {
  c(
    "#67001f",
    "#b2182b",
    "#d6604d",
    "#f4a582",
    "#fddbc7",
    "white",
    "white",
    "#d1e5f0",
    "#92c5de",
    "#4393c3",
    "#2166ac",
    "#053061"
  )
}


# Get shapefiles for Australia ------------------------------------------------
aus_map <- generate_aus_map_sf() 

## Breaks for DREAM uncertainty IQR
scale_size_limits <- plot_ready_percentage_difference_a3_on_off_data |>
  pull(IQR_CO2_impact_on_streamflow_percentage) |>
  range(na.rm = T) |> # can round up if I want to
  round(digits = 0)

percentage_IQR_breaks <- c(0, 2.5, 5, 10, 15, 50, 100) # custom breaks

dot_transparency <- 0.8

# Plotting function ============================================================

# I not convinced this function is required

make_CO2_streamflow_percentage_change_map <- function(data, title) {
  ## Generate Insets ===========================================================
  QLD_data <- data |>
    filter(state == "QLD")

  NSW_data <- data |>
    filter(state == "NSW")

  VIC_data <- data |>
    filter(state == "VIC")

  WA_data <- data |>
    filter(state == "WA")

  TAS_data <- data |>
    filter(state == "TAS")



  ### Generate inset plots #######################################################
  inset_dot_size <- 1.8

  inset_plot_QLD <- aus_map |>
    filter(state == "QLD") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = QLD_data,
      aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent),
      show.legend = FALSE,
      size = inset_dot_size,
      alpha = dot_transparency,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    guides(size = guide_bins(show.limits = TRUE)) +
    theme_void()


  inset_plot_NSW <- aus_map |>
    filter(state == "NSW") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = NSW_data,
      aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent),
      show.legend = FALSE,
      size = inset_dot_size,
      alpha = dot_transparency,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    guides(size = guide_bins(show.limits = TRUE)) +
    theme_void()



  inset_plot_VIC <- aus_map |>
    filter(state == "VIC") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = VIC_data,
      aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent),
      show.legend = FALSE,
      size = inset_dot_size,
      alpha = dot_transparency,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    guides(size = guide_bins(show.limits = TRUE)) +
    theme_void()



  inset_plot_WA <- aus_map |>
    filter(state == "WA") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = WA_data,
      aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent),
      show.legend = FALSE,
      size = inset_dot_size,
      alpha = dot_transparency,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    guides(size = guide_bins(show.limits = TRUE)) +
    theme_void()



  inset_plot_TAS <- aus_map |>
    filter(state == "TAS") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = TAS_data,
      aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent),
      show.legend = FALSE,
      size = inset_dot_size,
      alpha = dot_transparency,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    guides(size = guide_bins(show.limits = TRUE)) +
    theme_void()



  ## Put it together =============================================================
  single_map_aus <- aus_map |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = data,
      mapping = aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent),
      alpha = dot_transparency,
      size = inset_dot_size,
      colour = "black",
      shape = 21,
      inherit.aes = FALSE,
      stroke = 0.1
    ) +
    theme_bw() +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    # expand map
    coord_sf(xlim = c(95, 176), ylim = c(-60, 0)) +
    # magnify WA
    geom_magnify(
      from = c(114, 118, -35.5, -30),
      to = c(93, 112, -36, -10),
      shadow = FALSE,
      expand = 0,
      plot = inset_plot_WA,
      proj = "single"
    ) +
    # magnify VIC
    geom_magnify(
      # aes(from = state == "VIC"), # use aes rather than manually selecting area
      from = c(141, 149.5, -39, -34),
      to = c(95, 136, -38, -60),
      shadow = FALSE,
      plot = inset_plot_VIC,
      proj = "single"
    ) +
    # magnify QLD
    geom_magnify(
      from = c(145, 155, -29.2, -15),
      to = c(157, 178, -29.5, 1.5),
      shadow = FALSE,
      expand = 0,
      plot = inset_plot_QLD,
      proj = "single"
    ) +
    # magnify NSW
    geom_magnify(
      from = c(146.5, 154, -38, -28.1),
      to = c(157, 178, -61, -30.5),
      shadow = FALSE,
      expand = 0,
      plot = inset_plot_NSW,
      proj = "single"
    ) +
    # magnify TAS
    geom_magnify(
      from = c(144, 149, -40, -44),
      to = c(140, 155, -45, -61),
      shadow = FALSE,
      expand = 0,
      plot = inset_plot_TAS,
      proj = "single"
    ) +
    labs(
      x = NULL, # "Latitude",
      y = NULL, # "Longitude",
      fill = bquote('Average Impact of'~CO[2]~'on Streamflow (%)'), 
      size = "Percentage Impact Uncertainty (IQR)",
      title = {{ title }}
    ) +
    theme(
      legend.key = element_rect(fill = "white"),
      legend.title = element_text(hjust = 0.5),
      # legend.background = element_rect(colour = "black"), #this cuts off the negative sign
      axis.text = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.351, 0.9),
      legend.box = "horizontal", # side-by-side legends
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(margin = margin(l = 25, r = 0, t = 30, b = -30), size = 18) # push title into plot
    ) +
    guides(
      fill = guide_coloursteps(
        barwidth = unit(10, "cm"),
        show.limits = TRUE,
        even.steps = TRUE,
        title.position = "top",
        direction = "horizontal"
      ),
      size = guide_bins(
        override.aes = aes(stroke = 0.5),
        show.limits = TRUE,
        direction = "horizontal",
        title.position = "top", # warnings says its ignore these parameter - The warnings are wrong
        barwidth = unit(1, "cm")
      )
    )

  return(single_map_aus)
}






# Compare percentage changes using a histogram ---------------------------------
CO2_on_off_decade_histogram <- plot_ready_percentage_difference_a3_on_off_data |>
  mutate(
    decade = if_else(decade == 1, "1990-1999", "2012-2021")
  ) |>
  ggplot(aes(x = CO2_impact_on_streamflow_percent, fill = decade, colour = decade)) +
  geom_histogram(
    alpha = 0.25,
    position = "identity",
    bins = 40
  ) +
  labs(
    x = "Average Impact of CO2 on Streamflow (%)",
    y = "Frequency",
    colour = "Period",
    fill = "Period"
  ) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.8),
    legend.background = element_rect(colour = "black")
  )


plot_ready_percentage_difference_a3_on_off_data |>
  summarise(
    mean = mean(CO2_impact_on_streamflow_percent),
    median = median(CO2_impact_on_streamflow_percent),
    .by = decade
  )






# A time series plot of on/off of all gauges -----------------------------------
timeseries_plotting_data <- a3_on_off_difference_data |>
  pivot_longer(
    cols = starts_with("realspace"),
    names_to = "type",
    values_to = "streamflow"
  )

rainfall_runoff_plotting_data <- a3_on_off_difference_data |>
  pivot_longer(
    cols = starts_with("transformed"),
    names_to = "type",
    values_to = "streamflow"
  )

mega_timeseries_plot <- timeseries_plotting_data |>
  ggplot(aes(x = year, y = streamflow, colour = type)) +
  geom_line() +
  # uncomment for rainfall
  # geom_line(
  #  aes(x = year, y = precipitation),
  #  colour = "black",
  #  linetype = "dashed",
  #  inherit.aes = FALSE
  #  ) +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~gauge, scales = "free_y")





mega_rainfall_runoff_plot <- rainfall_runoff_plotting_data |>
  ggplot(aes(x = precipitation, y = streamflow, colour = type)) +
  geom_smooth(
    method = lm,
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.5
  ) +
  geom_point(size = 0.75) +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~gauge, scales = "free")









# Are non-forested catchments more vulnerable to change? -----------------------
# percentage change in Co2 vs prop forest
gauge_prop_forest_info <- gauge_information |>
  select(gauge, prop_forested)

all_years_percentage_difference_a3_on_off_data <- a3_on_off_difference_data |>
  filter(year %in% c(decade_1, decade_2)) |>
  # sum streamflow for each decade
  summarise(
    sum_decade_realspace_CO2_off_streamflow = sum(realspace_a3_off_streamflow),
    sum_decade_realspace_CO2_on_streamflow = sum(realspace_a3_on_streamflow),
    .by = c(gauge)
  ) |>
  # find the absolution and percentage difference
  mutate(
    CO2_impact_on_streamflow_mm_per_year = (sum_decade_realspace_CO2_on_streamflow - sum_decade_realspace_CO2_off_streamflow),
    CO2_impact_on_streamflow_percent = (CO2_impact_on_streamflow_mm_per_year / sum_decade_realspace_CO2_off_streamflow) * 100
  ) |>
  arrange(desc(CO2_impact_on_streamflow_percent)) |> # Large percentage changes are not tied to years_of_data
  # add prop forested
  left_join(
    gauge_prop_forest_info,
    by = join_by(gauge)
  )


streamflow_CO2_percentage_change_vs_prop_forested <- all_years_percentage_difference_a3_on_off_data |>
  ggplot(aes(x = prop_forested, y = CO2_impact_on_streamflow_percent)) +
  geom_point() +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
  labs(
    x = "Proportion of Catchment Forested",
    y = "Impact of CO2 on streamflow for entire gauge record (%)",
  ) +
  theme_bw()









# Compare percentage difference in best CO2 vs best non-CO2 model --------------


## Extract best CO2 and non-CO2 model streamflow ===============================
only_models_best_CO2_non_CO2_per_gauge <- best_CO2_non_CO2_per_gauge |>
  select(gauge, streamflow_model) |>
  distinct()

best_CO2_non_CO2_streamflow <- streamflow_results |>
  semi_join(
    only_models_best_CO2_non_CO2_per_gauge,
    by = join_by(gauge, streamflow_model)
  )


## Calculate streamflow difference =============================================
difference_best_CO2_non_CO2_streamflow <- best_CO2_non_CO2_streamflow |>
  select(gauge, year, realspace_modelled_streamflow, streamflow_model) |>
  mutate(
    streamflow_model = if_else(str_detect(streamflow_model, "CO2"), "CO2_model", "non_CO2_model")
  ) |>
  pivot_wider(
    names_from = streamflow_model,
    values_from = realspace_modelled_streamflow
  ) |>
  # Sum streamflow for each decade
  # We are interested only in 1990-1999 and 2012-2021
  mutate(
    decade = case_when(
      year %in% seq(from = 1990, to = 1999) ~ 1,
      year %in% seq(from = 2012, to = 2021) ~ 2,
      .default = NA
    )
  ) |>
  filter(!is.na(decade)) |>
  summarise(
    sum_decade_non_CO2_model_streamflow = sum(non_CO2_model),
    sum_decade_CO2_model_streamflow = sum(CO2_model),
    .by = c(gauge, decade)
  ) |>
  # Find percetange difference between modelled Co2 and non-CO2 streamflow
  mutate(
    CO2_impact_streamflow_mm_decade = sum_decade_CO2_model_streamflow - sum_decade_non_CO2_model_streamflow,
    CO2_impact_streamflow_percent = (CO2_impact_streamflow_mm_decade / sum_decade_non_CO2_model_streamflow) * 100
  ) |>
  arrange(desc(CO2_impact_streamflow_percent))


## Plot percentage difference graph ============================================
### add lat-lon and evidence ratio for plotting
plotting_best_CO2_non_CO2_streamflow <- difference_best_CO2_non_CO2_streamflow |>
  left_join(
    lat_lon_gauge,
    by = join_by(gauge)
  ) |>
  left_join(
    evidence_ratio_calc,
    by = join_by(gauge)
  ) |>
  # rename columns for I can use existing plotting function --> make_CO2_streamflow_percentage_change_map
  rename(
    CO2_impact_on_streamflow_percent = CO2_impact_streamflow_percent
  ) |>
  # filter evidence ratio
  filter(evidence_ratio > 100)


### redo limits
CO2_impact_on_streamflow_percent_limits <- plotting_best_CO2_non_CO2_streamflow |>
  pull(CO2_impact_on_streamflow_percent) |>
  make_limits() |>
  as.double()

hard_coded_breaks_CO2_impact_of_streamflow <- c(-75, -50, -25, -10, -1, 0, 1, 10, 25, 50, 75)

## Plots for 1990s and 2010s ===================================================

percentage_difference_CO2_model_non_CO2_model_1990s <- plotting_best_CO2_non_CO2_streamflow |>
  filter(decade == 1)

percentage_difference_CO2_model_non_CO2_model_2010s <- plotting_best_CO2_non_CO2_streamflow |>
  filter(decade == 2)

patchwork_CO2_model_and_non_CO2_model_percentage_differences <- (make_CO2_streamflow_percentage_change_map(percentage_difference_CO2_model_non_CO2_model_1990s, "1990-1999") | make_CO2_streamflow_percentage_change_map(percentage_difference_CO2_model_non_CO2_model_2010s, "2012-2021")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")


plotting_best_CO2_non_CO2_streamflow |>
  summarise(
    mean = mean(CO2_impact_on_streamflow_percent),
    median = median(CO2_impact_on_streamflow_percent),
    .by = decade
  )



# Rainfall-runoff and streamflow time graphs for best CO2 and nonCO2 models ----

# The CO2_on_off_streamflow_time_comparision only includes catchments
# where the evidence ratio is greater than zero
evidence_ratio_greater_0_gauges <- a3_on_off_difference_data |>
  pull(gauge) |>
  unique()


## Streamflow-time =============================================================
plot_streamflow_time_best_CO2_non_CO2 <- best_CO2_non_CO2_streamflow |>
  select(
    gauge,
    year,
    precipitation,
    realspace_observed_streamflow,
    realspace_modelled_streamflow,
    streamflow_model
  ) |>
  mutate(
    streamflow_model = if_else(str_detect(streamflow_model, "CO2"), "CO2_model", "non_CO2_model")
  ) |>
  pivot_wider(
    names_from = streamflow_model,
    values_from = realspace_modelled_streamflow
  ) |>
  pivot_longer(
    cols = realspace_observed_streamflow:CO2_model,
    names_to = "type",
    values_to = "streamflow"
  ) |>
  mutate(
    type = factor(type, levels = c("non_CO2_model", "CO2_model", "realspace_observed_streamflow"))
  ) |>
  left_join(
    evidence_ratio_calc,
    by = join_by(gauge)
  ) |>
  filter(gauge %in% evidence_ratio_greater_0_gauges) |> # this is done so the CO2_on_off_streamflow is an equal comparison
  ggplot(aes(x = year, y = streamflow, colour = type)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~gauge, scales = "free_y")








## Rainfall-runoff =============================================================

# The transformed_realspace is different depending on model used.
# This means 2 observed transformed streamflow is required.
# They overlap see all_rainfall_runoff plots


plot_rainfall_runoff_best_CO2_non_CO2 <- best_CO2_non_CO2_streamflow |>
  select(
    gauge,
    year,
    precipitation,
    transformed_observed_streamflow,
    transformed_modelled_streamflow,
    streamflow_model
  ) |>
  pivot_longer(
    cols = starts_with("transformed"),
    names_to = "type",
    values_to = "streamflow"
  ) |>
  # simplify streamflow names
  mutate(
    model = if_else(str_detect(streamflow_model, "CO2"), "CO2_model", "no_CO2_model")
  ) |>
  unite(
    col = "streamflow_type_and_model",
    type,
    model
  ) |>
  filter(gauge %in% evidence_ratio_greater_0_gauges) |> # this is done so the CO2_on_off_streamflow is an equal comparison
  ggplot(aes(x = precipitation, y = streamflow, colour = streamflow_type_and_model)) +
  geom_smooth(
    method = lm,
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.5
  ) +
  geom_point(size = 0.75) +
  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~gauge, scales = "free")











# Plot timeseries, observed, CO2_off, CO2_on, best_non_CO2 model ---------------
timeseries_data_obs_CO2_off_CO2_on <- timeseries_plotting_data |>
  select(gauge, year, precipitation, type, streamflow)


timeseries_data_obs_best_CO2_model_best_non_CO2_model <- best_CO2_non_CO2_streamflow |>
  select(
    gauge,
    year,
    precipitation,
    realspace_observed_streamflow,
    realspace_modelled_streamflow,
    streamflow_model
  ) |>
  mutate(
    streamflow_model = if_else(str_detect(streamflow_model, "CO2"), "CO2_model", "non_CO2_model")
  ) |>
  pivot_wider(
    names_from = streamflow_model,
    values_from = realspace_modelled_streamflow
  ) |>
  pivot_longer(
    cols = realspace_observed_streamflow:CO2_model,
    names_to = "type",
    values_to = "streamflow"
  ) |>
  mutate(
    type = factor(type, levels = c("non_CO2_model", "CO2_model", "realspace_observed_streamflow"))
  ) |>
  left_join(
    evidence_ratio_calc,
    by = join_by(gauge)
  ) |>
  filter(gauge %in% evidence_ratio_greater_0_gauges) |>
  # remove observed data and best_CO2 model to avoid duplicates while joining
  filter(type == "non_CO2_model") |>
  select(!evidence_ratio)


## Join the two datasets and make ready for plotting ===========================
all_timeseries_data <- timeseries_data_obs_CO2_off_CO2_on |>
  rbind(timeseries_data_obs_best_CO2_model_best_non_CO2_model) |>
  arrange(gauge) |>
  # give type good names
  mutate(
    type = case_when(
      type == "realspace_a3_on_streamflow" ~ "CO2 Model",
      type == "realspace_a3_off_streamflow" ~ "Counterfactual",
      type == "realspace_observed_streamflow" ~ "Observed",
      type == "non_CO2_model" ~ "non-CO2 Model",
      .default = NA
    )
  ) |>
  # set order
  mutate(
    type = factor(type, levels = c("Observed", "CO2 Model", "Counterfactual", "non-CO2 Model"))
  )


# Split up into a4 size and repeat
# Here
# It needs to be A, B, C, D labels rather than facet_labels
# A, B, C and D must be in the same relative location for all plots
# Save A, B, C and D catchments. Print out caption
# Colours are not good - need colour blind safe selection

make_facet_labels <- function(data, facet_column, x_axis_column, y_axis_column, label_type = LETTERS, hjust = 0, vjust = 0) {
  
  # The embrace operator does not work correctly in summarise i.e., max({{ y_axis_column }})
  # Link: https://forum.posit.co/t/embrace-operator-for-tidy-selection-vs-data-masking/173084
  # Possible cause: {{ y_axis_column }} isn't unquoting when it's doing the mutate 
  # Work around using rlang::ensym
  col <- rlang::ensym(y_axis_column)
  
  data |> 
    summarise(
      ylab = max(!!col),
      .by = {{ facet_column }}
    ) |> 
    # Add xlab - constant x-axis
    add_column(
      xlab = data |> pull(x_axis_column) |> min(),
      .before = 2
    ) |>  # add row numbers to tibble
    mutate(
      row_number = row_number(),
      .before = 1
    ) |>  # add label type based on row number
    mutate(
      label_name = label_type[row_number]
    ) |> 
    # apply hjust and vjust
    mutate(
      xlab = xlab + (xlab * hjust),
      ylab = ylab + (ylab * vjust)
    )

}





plot_and_save_timeseries_data <- function(plotting_data, label_data, identifier) {
  
  plot <- plotting_data |> 
    ggplot(aes(x = year, y = streamflow, colour = type, shape = type)) +
    geom_line(alpha = 0.8) +
    geom_point() +
    geom_text(
      mapping = aes(x = xlab, y = ylab, label = label_name),
      data = label_data,
      inherit.aes = FALSE,
      fontface = "bold",
      size = 10,
      size.unit = "pt"
    ) +
    scale_colour_brewer(palette = "Set1") +
    labs(
      x = "Time (Year)", 
      y = "Streamflow (mm)", 
      colour = NULL, 
      shape = NULL
    ) +
    scale_x_continuous(expand = c(0.01,0.01)) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      text = element_text(family = "sans", size = 9), # default fonts are serif, sans and mono, text size is in pt
      strip.background = element_blank(), # remove facet_strip gauge numbers
      strip.text = element_blank() # remove facet_strip gauge numbers
    ) +
    facet_wrap(~gauge, ncol = 1, scales = "free_y") 
  
  ggsave(
    filename = paste0("streamflow_timeseries_supp_plot_",identifier,".pdf"),
    path = "Figures/Supplementary",
    device = "pdf",
    plot = plot,
    width = 183,
    height = 197,
    units = "mm"
  )
}





# create txt file with figure captions for overleaf
create_caption <- function(label_data, identifier) {
  
  abc <- label_data |> pull(label_name)
  gauge <- label_data |> pull(gauge)
  gauge_abc <- paste0(gauge, " (", abc, ")")
  # concatenate everything but last value
  start_gauge_abc <- paste0(gauge_abc[1:(length(gauge_abc) - 2)], ", ", collapse = "")
  end_gauge_abc <- paste0(gauge_abc[(length(gauge_abc) - 1)], " and ", gauge_abc[length(gauge_abc)], ".")
  gauge_text <- paste(c(start_gauge_abc, end_gauge_abc), collapse = "")
  
  cat("\\begin{figure}") 
  cat("\n")
  cat("\t\\centering")
  cat("\n")
  cat(paste0("\t\\includegraphics[width=\\textwidth]{Figures/streamflow_timeseries_supp_plot_", identifier, ".pdf}"))
  cat("\t\n")
  # The line below must change
  cat(paste0("\t\\caption{\\textbf{The impact of CO$_2$ on the streamflow timeseries for gauges ", gauge_text, "} The streamflow timeseries compares observed streamflow (Observed), modelled streamflow using a model that includes CO$_2$ (CO$_2$ model), modelled streamflow using a model that includes CO$_2$ with CO$_2$ turned off (Counterfactual) and modelled streamflow using a model that does not include CO$_2$ (non-CO$_2$ Model).}"))
  cat("\n")
  # The line below must change
  cat(paste0("\t\\label{fig:supp_streamflow_timeseries_", identifier, "}")) 
  cat("\n")
  cat("\\end{figure}")
  cat("\n")
  cat("\n")

}



save_plot_and_caption_timeseries_data <- function(data_chunk, identifier) {
  
  label_data <- make_facet_labels(
    data = data_chunk,
    facet_column = "gauge",
    x_axis_column = "year",
    y_axis_column = "streamflow",
    label_type = LETTERS,
    hjust = 0.0005,
    vjust = -0.05
  )
  
  plot_and_save_timeseries_data(
    plotting_data = data_chunk, 
    label_data = label_data, 
    identifier = identifier
    )
  
  create_caption(
    label_data = label_data,
    identifier = identifier
  )
  
}

# divide all_timeseries_data into lots of 7 using split
chunk <- 7
n <- length(CO2_gauges)
split_group <- rep(rep(1:ceiling(n/chunk), each = chunk))[1:n]
split_tibble <- tibble(
  "gauge" = CO2_gauges,
  "split" = split_group
)

# left_join and split
all_timeseries_data <- all_timeseries_data |> 
  left_join(
    split_tibble,
    by = join_by(gauge)
  )

chunked_timeseries_data <- all_timeseries_data |> # converting table to list by groups https://stackoverflow.com/questions/7060272/split-up-a-dataframe-by-number-of-rows
  group_by(split) |> 
  group_map(~ .x)








# Combine percentage change map with timeseries --------------------------------
## Select handful of catchments for main paper =================================

# short list - 230210, 415226, 617003, 701002, 235234, 231211, 303203, 407246, 407253, 208004, 406214, 614044, 405230, 227210, 606195
short_list_catchments <- c("401210", "606195", "701002", "407246") # select 4



## Make label table ============================================================
map_label_table <- lat_lon_gauge |> 
  filter(gauge %in% short_list_catchments) |> 
  mutate(label_name = LETTERS[1:4])


font_size <- 3L # default size is GeomLabel$default_aes$size = 3.88


## Manually add labels in the correct spots in function
make_CO2_streamflow_percentage_change_map_uncertainty <- function(data, title) {
  
  ## Generate Insets ===========================================================
  QLD_data <- data |>
    filter(state == "QLD")
  
  NSW_data <- data |>
    filter(state == "NSW")
  
  VIC_data <- data |>
    filter(state == "VIC")
  
  WA_data <- data |>
    filter(state == "WA")
  
  TAS_data <- data |>
    filter(state == "TAS")
  
  
  
  ### Generate inset plots #######################################################
  inset_dot_size <- 1.8
  
  inset_plot_QLD <- aus_map |>
    filter(state == "QLD") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = QLD_data,
      aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent, size = IQR_CO2_impact_on_streamflow_percentage),
      show.legend = FALSE,
      alpha = dot_transparency,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    guides(size = guide_bins(show.limits = TRUE)) +
    theme_void()
  
  
  inset_plot_NSW <- aus_map |>
    filter(state == "NSW") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = NSW_data,
      aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent, size = IQR_CO2_impact_on_streamflow_percentage),
      show.legend = FALSE,
      alpha = dot_transparency,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    guides(size = guide_bins(show.limits = TRUE)) +
    theme_void()
  
  
  # Labels for VIC
  VIC_map_label_table <- map_label_table |> 
    filter(state == "VIC")
  
  inset_plot_VIC <- aus_map |>
    filter(state == "VIC") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = VIC_data,
      aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent, size = IQR_CO2_impact_on_streamflow_percentage),
      show.legend = FALSE,
      alpha = dot_transparency,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    geom_text(
      data = VIC_map_label_table,
      aes(x = lon, y = lat, label = label_name),
      nudge_x = -0.35,
      size = font_size
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    guides(size = guide_bins(show.limits = TRUE)) +
    theme_void()
  
  
  # Labels for WA
  WA_map_label_table <- map_label_table |> 
    filter(gauge == "606195") # 701002 not in subset plot so exclude
  
  inset_plot_WA <- aus_map |>
    filter(state == "WA") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = WA_data,
      aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent, size = IQR_CO2_impact_on_streamflow_percentage),
      show.legend = FALSE,
      alpha = dot_transparency,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    geom_text(
      data = WA_map_label_table,
      aes(x = lon, y = lat, label = label_name),
      nudge_y = -0.25,
      size = font_size
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    guides(size = guide_bins(show.limits = TRUE)) +
    theme_void()
  
  
  
  inset_plot_TAS <- aus_map |>
    filter(state == "TAS") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = TAS_data,
      aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent, size = IQR_CO2_impact_on_streamflow_percentage),
      show.legend = FALSE,
      alpha = dot_transparency,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    guides(size = guide_bins(show.limits = TRUE)) +
    theme_void()
  
  
  
  ## Put it together =============================================================
  
  # Big map label table
  big_map_label_table <- map_label_table |> 
    filter(gauge == "701002")
  
  single_map_aus <- aus_map |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = data,
      mapping = aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent, size = IQR_CO2_impact_on_streamflow_percentage),
      alpha = dot_transparency,
      colour = "black",
      shape = 21,
      inherit.aes = FALSE,
      stroke = 0.1
    ) +
    geom_text(
      data = big_map_label_table,
      aes(x = lon, y = lat, label = label_name),
      nudge_y = 2,
      size = font_size
    ) +
    theme_bw() +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    # expand map
    coord_sf(xlim = c(95, 176), ylim = c(-60, 0)) +
    # magnify WA
    geom_magnify(
      from = c(114, 118, -35.5, -30),
      to = c(93, 112, -36, -10),
      shadow = FALSE,
      expand = 0,
      plot = inset_plot_WA,
      proj = "single"
    ) +
    # magnify VIC
    geom_magnify(
      # aes(from = state == "VIC"), # use aes rather than manually selecting area
      from = c(141, 149.5, -39, -34),
      to = c(95, 136, -38, -60),
      shadow = FALSE,
      plot = inset_plot_VIC,
      proj = "single"
    ) +
    # magnify QLD
    geom_magnify(
      from = c(145, 155, -29.2, -15),
      to = c(157, 178, -29.5, 1.5),
      shadow = FALSE,
      expand = 0,
      plot = inset_plot_QLD,
      proj = "single"
    ) +
    # magnify NSW
    geom_magnify(
      from = c(146.5, 154, -38, -28.1),
      to = c(157, 178, -61, -30.5),
      shadow = FALSE,
      expand = 0,
      plot = inset_plot_NSW,
      proj = "single"
    ) +
    # magnify TAS
    geom_magnify(
      from = c(144, 149, -40, -44),
      to = c(140, 155, -45, -61),
      shadow = FALSE,
      expand = 0,
      plot = inset_plot_TAS,
      proj = "single"
    ) +
    labs(
      x = NULL, # "Latitude",
      y = NULL, # "Longitude",
      fill = bquote('Average Impact of'~CO[2]~'on Streamflow (%)'), 
      size = "Percentage Impact Uncertainty (IQR)",
      title = {{ title }}
    ) +
    theme(
      legend.key = element_rect(fill = "white"),
      legend.title = element_text(hjust = 0.5),
      # legend.background = element_rect(colour = "black"), #this cuts off the negative sign
      axis.text = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.351, 0.9),
      legend.box = "horizontal", # side-by-side legends
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(margin = margin(l = 25, r = 0, t = 30, b = -30), size = 18) # push title into plot
    ) +
    guides(
      fill = guide_coloursteps(
        barwidth = unit(10, "cm"),
        show.limits = TRUE,
        even.steps = TRUE,
        title.position = "top",
        direction = "horizontal"
      ),
      size = guide_bins(
        override.aes = aes(stroke = 0.5),
        show.limits = TRUE,
        direction = "horizontal",
        title.position = "top", # warnings says its ignore these parameter - The warnings are wrong
        barwidth = unit(1, "cm")
      )
    )
  
  return(single_map_aus)
}






## Combining percentage change and timeseries plots ============================

### Order percentage impact uncertainty ########################################
# force small dots to be on top and large dots on the bottom

plot_ready_percentage_difference_a3_on_off_1990s <- plot_ready_percentage_difference_a3_on_off_data |>
  filter(decade == 1) |>
  filter(evidence_ratio > 100) |>  # 100 is moderately strong
  arrange(desc(IQR_CO2_impact_on_streamflow_percentage))

plot_ready_percentage_difference_a3_on_off_2010s <- plot_ready_percentage_difference_a3_on_off_data |>
  filter(decade == 2) |>
  filter(evidence_ratio > 100) |>  # 100 is moderately strong
  arrange(desc(IQR_CO2_impact_on_streamflow_percentage))

  

### Percentage change components ###############################################
percent_change_1990 <- make_CO2_streamflow_percentage_change_map_uncertainty(
  plot_ready_percentage_difference_a3_on_off_1990s, 
  "1990-1999"
  )

percent_change_2012 <- make_CO2_streamflow_percentage_change_map_uncertainty(
  plot_ready_percentage_difference_a3_on_off_2010s, 
  "2012-2021"
  )


## timeseries components #######################################################

### annotating facets
facet_annotation <- all_timeseries_data |> 
  filter(gauge %in% short_list_catchments) |>
  summarise(
    streamflow = max(streamflow) - (max(streamflow) * 0.1),
    .by = gauge
  ) |> 
  add_column(
    year = 1960,
    label_name = LETTERS[1:4]
  )


### Shading decades 
shade_decade_1 <- all_timeseries_data |> 
  filter(gauge %in% short_list_catchments) |>
  group_by(gauge) |> 
  mutate(upper = max(streamflow) * 1.2) |> 
  filter(year %in% decade_1)

shade_decade_2 <- all_timeseries_data |> 
  filter(gauge %in% short_list_catchments) |>
  group_by(gauge) |> 
  mutate(upper = max(streamflow) * 1.2) |> 
  filter(year %in% decade_2) 
# 606195 is missing 2021
# 407246 is missing 2020 and 2021

# easiest solution is extract a year - replace year and streamflow and rbind
missing_years_606195 <- shade_decade_2 |> 
  filter(gauge == "606195") |> 
  filter(year == 2020) |> 
  mutate(
    year = 2021,
    precipitation = NA,
    streamflow = NA
  )

missing_years_407246 <- shade_decade_2 |> 
  filter(gauge == "407246") |> 
  filter(year %in% c(2018, 2019)) |> 
  mutate(
    year = if_else(year == 2018, 2020, 2021),
    precipitation = NA,
    streamflow = NA
  )

shade_decade_2 <- rbind(shade_decade_2, missing_years_606195, missing_years_407246)


### Plotting ###################################################################
bottom <- all_timeseries_data |>
  filter(gauge %in% short_list_catchments) |>
  ggplot(aes(x = year, y = streamflow, colour = type)) +
  geom_area(
    aes(x = year, y = upper),
    inherit.aes = FALSE,
    data = shade_decade_1,
    alpha = 0.08
  ) +
  geom_area(
    aes(x = year, y = upper),
    inherit.aes = FALSE,
    data = shade_decade_2,
    alpha = 0.08
  ) +
  geom_line(alpha = 0.8) +
  geom_label(
    aes(x = year, y = streamflow, label = label_name),
    data = facet_annotation,
    inherit.aes = FALSE,
    fill = NA,
    label.size = NA
    ) +
  #scale_colour_brewer(palette = "Set1") +
  scale_colour_brewer(
    labels = c(
      "Observed",
      bquote(~CO[2]~"Model"),
      "Counterfactual",
      bquote("non-"*CO[2]~"Model")
    ),
    palette = "Set1"
    #values = c("red", "green", "blue", "orange")
  ) +
  labs(x = "Year", y = "Streamflow (mm)", colour = "Streamflow Timeseries") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(), # remove strip labels from faceting
    strip.text = element_blank(), # remove strip labels from faceting
    panel.grid.minor = element_blank(), # remove minor gridlines
    axis.title = element_text(size = 10)
    ) +
  scale_y_continuous(expand = c(0, 0)) + # remove y-axis padding
  facet_wrap(~gauge, scales = "free_y", ncol = 2)



### Combining ##################################################################

# Alternative methods attempted:
## - plot_spacer and | + / operations - too much faff
## - individually plotting timeseries - gaps

# Defining layout it the best method as I have direct control

# use the area() constructor
# top, left, bottom, right bounds (t < b and l < r)
layout <- c(
  area(t = 1, l = 1, b = 3, r = 3), # 1990s percentage change
  area(t = 1, l = 4, b = 3, r = 6), # 2010s percentage change
  area(t = 4, l = 1, b = 4, r = 6) # timeseries
)

plot(layout) # check the patches are working

streamflow_percentage_difference_with_timeseries <- (percent_change_1990 + percent_change_2012 + bottom) + 
  plot_layout(design = layout, guides = "collect") & 
  theme(legend.position = "bottom") &
  guides(colour = guide_legend(title.hjust = 0.5, title.position = "top", ncol = 2))



# Instead of having timeseries spread across multiple pages - stick em together ---
plot_all_timeseries_data <- all_timeseries_data |>
  ggplot(aes(x = year, y = streamflow, colour = type)) +
  geom_line(alpha = 0.8) +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "Year", y = "Streamflow (mm)", colour = NULL) +
  theme_bw() +
  facet_wrap(~gauge, scales = "free_y")




# Saving all figures here ------------------------------------------------------
# 1. Main --> streamflow_percentage_difference_with_timeseries.pdf
# 2. Supplementary --> streamflow_timeseries_data.pdf and captions
# 3. Supplementary --> streamflow_percentage_difference_best_CO2_model_vs_best_non_CO2_model.pdf
# 4. Other --> CO2_on_off_decade_histogram.pdf
# 5. Other --> CO2_on_off_rainfall_runoff_comparison.pdf (mega) -  only catchments where CO2 model is better than non-CO2 model (evidence ratio > 0)
# 6. Other --> streamflow_CO2_percentage_change_vs_prop_forested.pdf
# 7. Other --> CO2_model_vs_non_CO2_model_rainfall_runoff.pdf - The transformed_realspace is different depending on model used. This means 2 observed transformed streamflow is required.
# 8. Other --> complete_timeseries_plot.pdf

# 1.
ggsave(
  filename = "./Figures/Main/streamflow_percentage_difference_with_timeseries.pdf",
  plot = streamflow_percentage_difference_with_timeseries,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)

# 2.
# Supplementary time series figures and captions:
sink(file = "Figures/Supplementary/streamflow_time_captions_supp.txt") # filename must change
iwalk(
  .x = chunked_timeseries_data,
  .f = save_plot_and_caption_timeseries_data
)
sink()
stop_here()

# 3.
ggsave(
  filename = "streamflow_percentage_difference_best_CO2_model_vs_best_non_CO2_model.pdf",
  path = "Figures/Supplementary",
  plot = patchwork_CO2_model_and_non_CO2_model_percentage_differences,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)

# 4.
ggsave(
  filename = "CO2_on_off_decade_histogram.pdf",
  plot = CO2_on_off_decade_histogram,
  device = "pdf",
  path = "Figures/Other",
  width = 297,
  height = 210,
  units = "mm"
)

# 5.
ggsave(
  filename = "CO2_on_off_rainfall_runoff_comparison.pdf",
  path = "Figures/Other",
  device = "pdf",
  plot = mega_rainfall_runoff_plot,
  width = 1189,
  height = 841,
  units = "mm"
)


# 6.
ggsave(
  filename = "streamflow_CO2_percentage_change_vs_prop_forested.pdf",
  path = "Figures/Other",
  device = "pdf",
  plot = streamflow_CO2_percentage_change_vs_prop_forested,
  width = 297,
  height = 210,
  units = "mm"
)


# 7.
ggsave(
  filename = "CO2_model_vs_non_CO2_model_rainfall_runoff.pdf",
  path = "Figures/Other",
  device = "pdf",
  plot = plot_rainfall_runoff_best_CO2_non_CO2,
  width = 1189,
  height = 841,
  units = "mm"
)

# 8.
ggsave(
  filename = "complete_timeseries_plot.pdf",
  path = "Figures/Other",
  device = "pdf",
  plot = plot_all_timeseries_data,
  width = 1189,
  height = 841,
  units = "mm"
  )










