# CO2 streamflow impact analysis

# Figures produced in this R file ----------------------------------------------

# 1. Main --> CO2_streamflow_on_off_decade_comparison.pdf 
# 2. Supplementary --> CO2_on_off_decade_histogram.pdf
# 3. Supplementary --> CO2_on_off_rainfall_runoff_comparison.pdf (mega)
# 4. Supplementary --> CO2_on_off_streamflow_time_comparison.pdf (mega)
# 5. Supplementary --> CO2_on_off_streamflow_time_comparison_with_rainfall.pdf (mega)
# 6. Supplementary --> streamflow_CO2_percentage_change_vs_prop_forested.pdf








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
  mutate(year = as.integer(year))

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
    names_from = CO2_or_non_CO2_model,
    values_from = streamflow_model
  )


### Comparison table ###
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
  
  a <- parameters[length(parameters) - 2]
  b <- parameters[length(parameters) - 1]
  
  a3_off_realspace_streamflow <- inverse_log_sinh_transform(
    a = a,
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

years_of_intrest <- c(seq(1990, 1999), seq(2012, 2021))

percentage_difference_a3_on_off_data <- a3_on_off_difference_data |> 
  filter(year %in% years_of_intrest) |> 
  # add decade group for summarising
  mutate(
    decade = year %/% 1000
  ) |> 
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
  )




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
  as.double()

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


# Plotting function ============================================================
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
  inset_dot_size = 1.8
  
  inset_plot_QLD <- aus_map |>
    filter(state == "QLD") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = QLD_data,
      aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent),
      show.legend = FALSE,
      size = inset_dot_size,
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
    theme_void()
  
  
  
  ## Put it together =============================================================
  single_map_aus <- aus_map |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = data,
      mapping = aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent),
      size = 2.2,
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
      fill = "Average impact of CO2 on streamflow per year (%)",
      title = {{ title }}
    ) +
    theme(
      legend.key = element_rect(fill = "grey80"),
      legend.title = element_text(hjust = 0.5),
      #legend.background = element_rect(colour = "black"), this cuts off the negative sign
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
        barwidth = unit(15, "cm"),
        show.limits = TRUE,
        even.steps = TRUE,
        title.position = "top",
        direction = "horizontal"
      )
    )
  
  return(single_map_aus)

}

# Patchwork results
plot_ready_percentage_difference_a3_on_off_1990s <- plot_ready_percentage_difference_a3_on_off_data |> 
  filter(decade == 1) |> 
  filter(evidence_ratio > 100) # 100 is moderately strong

plot_ready_percentage_difference_a3_on_off_2010s <- plot_ready_percentage_difference_a3_on_off_data |> 
  filter(decade == 2) |> 
  filter(evidence_ratio > 100) # 100 is moderately strong
  
patchwork_percentage_differences <- (make_CO2_streamflow_percentage_change_map(plot_ready_percentage_difference_a3_on_off_1990s, "1990-1999") | make_CO2_streamflow_percentage_change_map(plot_ready_percentage_difference_a3_on_off_2010s, "2012-2021")) + 
  plot_layout(guides = "collect") & theme(legend.position = "bottom")




ggsave(
  filename = "./Figures/Main/streamflow_percentage_difference_CO2_on_off.pdf",
  plot = patchwork_percentage_differences,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)



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

ggsave(
  filename = "CO2_on_off_decade_histogram.pdf",
  plot = CO2_on_off_decade_histogram,
  device = "pdf",
  path = "Figures/Supplementary",
  width = 297,
  height = 210,
  units = "mm"
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

ggsave(
  filename = "CO2_on_off_streamflow_time_comparison_with_rainfall.pdf",
  path = "Figures/Supplementary",
  device = "pdf",
  plot = mega_timeseries_plot,
  width = 1189,
  height = 841,
  units = "mm"
)



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


ggsave(
  filename = "CO2_on_off_rainfall_runoff_comparison.pdf",
  path = "Figures/Supplementary",
  device = "pdf",
  plot = mega_rainfall_runoff_plot,
  width = 1189,
  height = 841,
  units = "mm"
)






# Are non-forested catchments more vulnerable to change? -----------------------
# percentage change in Co2 vs prop forest 
gauge_prop_forest_info <- gauge_information |> 
  select(gauge, prop_forested)

all_years_percentage_difference_a3_on_off_data <- a3_on_off_difference_data |> 
  filter(year %in% years_of_intrest) |> 
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
  arrange(desc(CO2_impact_on_streamflow_percent)) |>  # Large percentage changes are not tied to years_of_data
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


ggsave(
  filename = "streamflow_CO2_percentage_change_vs_prop_forested.pdf",
  path = "Figures/Supplementary",
  device = "pdf",
  plot = streamflow_CO2_percentage_change_vs_prop_forested,
  width = 297,
  height = 210,
  units = "mm"
)




# Other ideas ------------------------------------------------------------------
# Average not by decade?







# TIM CURVES




# Ignore this bit here ------------

# Sanity check the rainfall-runoff and streamflow-time of large differences
# from decade comparison
# Input: gauge
# Output: rainfall-runoff and streamflow time
# test_gauge 207015
# get gauges from --> average_percent_diff_by_decade

tims_curves_data <- altered_realspace_streamflow_data_a3_off |>
  left_join(
    observed_streamflow,
    by = join_by(gauge, precipitation)
  ) |>
  select(gauge, year, precipitation, observed_boxcox_streamflow, a3_off_modelled_boxcox_streamflow, modelled_boxcox_streamflow, realspace_a3_on, realspace_a3_off, q_mm) |>
  rename(
    boxcox_a3_on = modelled_boxcox_streamflow,
    boxcox_a3_off = a3_off_modelled_boxcox_streamflow,
    observed_realspace_streamflow = q_mm
  ) |>
  mutate(
    a3_on_off_diff = realspace_a3_on - realspace_a3_off,
    standardised_a3_on_off_diff = a3_on_off_diff / observed_realspace_streamflow,
    standardised_precipitation = precipitation / mean(precipitation),
    .by = gauge
  )


# Do for all catchments with CO2 as best model
# Would be good to colour the dots based on type of change (neg-slope, pos-slope etc)
type_of_CO2_change <- best_model_per_gauge |>
  filter(contains_CO2) |>
  filter(str_detect(parameter, "a3")) |>
  mutate(
    a3_sign = if_else(sign(parameter_value) == -1, "negative", "positive"),
    change_type_no_sign = str_remove(parameter, "a3_")
  ) |>
  # make new column joining sign of parameter and type of parameter
  unite(
    col = change_type,
    a3_sign, change_type_no_sign
  ) |>
  select(gauge, change_type)

# Check if a3 off is larger than precipitation - mark this on the graphs
check_precip_and_a3_values <- tims_curves_data |>
  mutate(
    a3_off_bigger_than_precip = realspace_a3_off > precipitation
  ) |>
  summarise(
    streamflow_bigger_than_precip = any(a3_off_bigger_than_precip),
    .by = gauge
  )

# Join back together for plotting
tims_curves_data <- tims_curves_data |>
  left_join(
    type_of_CO2_change,
    by = join_by(gauge)
  ) |>
  left_join(
    check_precip_and_a3_values,
    by = join_by(gauge)
  ) |>
  mutate(
    change_type = factor(change_type, levels = c("negative_intercept", "positive_intercept", "negative_slope", "positive_slope"))
  )





absolute_curves <- tims_curves_data |>
  ggplot(aes(x = precipitation, y = a3_on_off_diff, colour = change_type, shape = streamflow_bigger_than_precip)) +
  geom_point() +
  labs(x = "Precipitation", y = "Q_CO2_on - Q_CO2_off") +
  theme_bw() +
  facet_wrap(~gauge, scales = "free") +
  theme(
    legend.position = "top"
  )

standardised_curves <- tims_curves_data |>
  # filter(gauge %in% c("112101B", "405228", "610001", "116011A")) |>
  ggplot(aes(x = standardised_precipitation, y = standardised_a3_on_off_diff, colour = change_type, shape = streamflow_bigger_than_precip)) + #
  geom_point() +
  labs(x = "P / P_ave", y = "(Q_CO2_on - Q_CO2_off) / Q_obs") +
  theme_bw() +
  facet_wrap(~gauge, scales = "free") +
  theme(
    legend.position = "top"
  )



ggsave(
  filename = "tim_curves_absolute.pdf",
  plot = absolute_curves,
  device = "pdf",
  path = "Graphs/Supplementary_Figures",
  width = 1189,
  height = 841,
  units = "mm"
)

ggsave(
  filename = "tim_curves_standardised.pdf",
  plot = standardised_curves,
  device = "pdf",
  path = "Graphs/Supplementary_Figures",
  width = 1189,
  height = 841,
  units = "mm"
)
















# SENS SLOPE ANALYSIS





# Rate of change CO2-on vs CO2-off using sens slope ----------------------------
rate_of_change <- altered_realspace_streamflow_data_a3_off |>
  select(gauge, year, realspace_a3_on, realspace_a3_off) |> # real space
  # select(gauge, year, a3_off_modelled_boxcox_streamflow, modelled_boxcox_streamflow) |> # boxcox
  mutate(
    relative_difference = realspace_a3_on - realspace_a3_off, # real space
    # relative_difference = modelled_boxcox_streamflow - a3_off_modelled_boxcox_streamflow
  ) |>
  # check for missing years - there are missing years - use my sen's slope to account for missing
  mutate(
    n = max(year) - year,
    .by = gauge
  ) |>
  mutate(
    n_1 = lag(n - 1),
    .by = gauge
  ) |>
  mutate(
    check = if_else(n - n_1 == 0, TRUE, FALSE)
  ) |>
  # I am not subtracting negatives
  mutate(
    is_negative_a3_on = realspace_a3_on < 0,
    is_negative_a3_off = realspace_a3_off < 0
  )

# If there is only a single occurrence then the following code errors
remove_single_entires <- rate_of_change |>
  summarise(
    n = n(),
    .by = gauge
  ) |>
  filter(n == 1)

rate_of_change <- rate_of_change |>
  anti_join(
    remove_single_entires,
    by = join_by(gauge)
  )



# Vectors into sens.slope must be:
# - in chronological order
# - be a timeseries object
# - make two checks
# - I think this is already being done - the rate of change is in gauge and year order

# Problem - deleting rows
# some things are not in order


# Prove to Tim I am correct
# - Excel
# - My own sen slope
# - Show him on powerpoint



sens_slope_estimator <- function(x) { # function(x, t)
  sens.slope(x)$estimates |> unname()
}


test_flow <- rate_of_change |>
  filter(gauge == "407253") |>
  pull(relative_difference)

test_year <- rate_of_change |>
  filter(gauge == "407253") |>
  pull(year)


sum(1:(length(test) - 1))


# My own sen slope function
my_sens_slope <- function(x, t) {
  # the length of x and t must be the same
  stopifnot(length(x) == length(t))
  
  # t must be continuous
  t_lag_1 <- lag(t, n = 1L)
  diff_t <- t - t_lag_1
  stopifnot(any(diff_t[-1] == 1))
  
  # pre-allocate array
  length_pre_allocation <- sum(1:(length(x) - 1))
  d <- numeric(length = length_pre_allocation)
  seperate_index <- 1
  
  for (j in 2:length(x)) {
    for (i in j:length(x)) {
      d[seperate_index] <- (x[i] - x[j - 1]) / (t[i] - t[j - 1])
      seperate_index <- seperate_index + 1
    }
  }
  
  # return(d)
  return(median(d, na.rm = TRUE))
}



my_sens_slope(x = test_flow, t = test_year)
sens_slope_estimator(test) # there is a non-continuous year. Does this matter?
sens.slope(test)


gauge_rate_of_change <- rate_of_change |>
  summarise(
    sens_slope = my_sens_slope(x = relative_difference, t = year),
    n = n(),
    .by = gauge
  ) |>
  arrange(sens_slope) |>
  left_join(
    lat_lon_gauge,
    by = join_by(gauge)
  )

# Put the sens_slope on a map mm/year ------------------------------------------
sens_slope_palette <- function(x) {
  c(
    "#67001f",
    "#b2182b",
    "#d6604d",
    "#f4a582",
    "#fddbc7",
    "white",
    "#d1e5f0",
    "#92c5de",
    "#4393c3",
    "#2166ac",
    "#053061"
  )
}

gauge_rate_of_change |>
  pull(sens_slope) |>
  quantile()

## Generate Insets =============================================================
### Filter data by state #######################################################

QLD_data <- gauge_rate_of_change |>
  filter(state == "QLD")

NSW_data <- gauge_rate_of_change |>
  filter(state == "NSW")

VIC_data <- gauge_rate_of_change |>
  filter(state == "VIC")

WA_data <- gauge_rate_of_change |>
  filter(state == "WA")

TAS_data <- gauge_rate_of_change |>
  filter(state == "TAS")


### Generate inset plots #######################################################

inset_plot_QLD <- aus_map |>
  filter(state == "QLD") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = QLD_data,
    aes(x = lon, y = lat, fill = sens_slope),
    show.legend = FALSE,
    size = 2.5,
    colour = "black",
    stroke = 0.1,
    shape = 21
  ) +
  binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
    aesthetics = "fill",
    palette = sens_slope_palette,
    breaks = c(-5, -2.5, -1, -0.5, -0.001, 0.001, 0.5, 1, 2.5, 5), # range()
    limits = c(-15, 10),
    show.limits = TRUE,
    guide = "colorsteps"
  ) +
  theme_void()


inset_plot_NSW <- aus_map |>
  filter(state == "NSW") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = NSW_data,
    aes(x = lon, y = lat, fill = sens_slope),
    show.legend = FALSE,
    size = 2.5,
    colour = "black",
    stroke = 0.1,
    shape = 21
  ) +
  binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
    aesthetics = "fill",
    palette = sens_slope_palette,
    breaks = c(-5, -2.5, -1, -0.5, -0.001, 0.001, 0.5, 1, 2.5, 5), # range()
    limits = c(-15, 10),
    show.limits = TRUE,
    guide = "colorsteps"
  ) +
  theme_void()




inset_plot_VIC <- aus_map |>
  filter(state == "VIC") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = VIC_data,
    aes(x = lon, y = lat, fill = sens_slope),
    show.legend = FALSE,
    size = 2.5,
    colour = "black",
    stroke = 0.1,
    shape = 21
  ) +
  binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
    aesthetics = "fill",
    palette = sens_slope_palette,
    breaks = c(-5, -2.5, -1, -0.5, -0.001, 0.001, 0.5, 1, 2.5, 5), # range()
    limits = c(-15, 10),
    show.limits = TRUE,
    guide = "colorsteps"
  ) +
  theme_void()




inset_plot_WA <- aus_map |>
  filter(state == "WA") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = WA_data,
    aes(x = lon, y = lat, fill = sens_slope),
    show.legend = FALSE,
    size = 2.5,
    colour = "black",
    stroke = 0.1,
    shape = 21
  ) +
  binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
    aesthetics = "fill",
    palette = sens_slope_palette,
    breaks = c(-5, -2.5, -1, -0.5, -0.001, 0.001, 0.5, 1, 2.5, 5), # range()
    limits = c(-15, 10),
    show.limits = TRUE,
    guide = "colorsteps"
  ) +
  theme_void()




inset_plot_TAS <- aus_map |>
  filter(state == "TAS") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = TAS_data,
    aes(x = lon, y = lat, fill = sens_slope),
    show.legend = FALSE,
    size = 2.5,
    colour = "black",
    stroke = 0.1,
    shape = 21
  ) +
  binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
    aesthetics = "fill",
    palette = sens_slope_palette,
    breaks = c(-5, -2.5, -1, -0.5, -0.001, 0.001, 0.5, 1, 2.5, 5), # range()
    limits = c(-15, 10),
    show.limits = TRUE,
    guide = "colorsteps"
  ) +
  theme_void()




## Put it together =============================================================

sens_slope_map_aus <- aus_map |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = gauge_rate_of_change,
    mapping = aes(x = lon, y = lat, fill = sens_slope),
    size = 3,
    colour = "black",
    shape = 21,
    inherit.aes = FALSE,
    stroke = 0.1
  ) +
  theme_bw() +
  binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
    aesthetics = "fill",
    palette = sens_slope_palette,
    breaks = c(-5, -2.5, -1, -0.5, -0.001, 0.001, 0.5, 1, 2.5, 5), # range()
    limits = c(-15, 10),
    show.limits = TRUE,
    guide = "colorsteps"
  ) +
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
    from = c(145, 155, -29.2, -16),
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
    fill = "Sen's Slope (mm/year)",
  ) +
  theme(
    legend.key = element_rect(fill = "grey80"),
    legend.title = element_text(hjust = 0.5),
    legend.background = element_rect(colour = "black"),
    legend.margin = margin(l = 10, r = 10, t = 5, b = 5), # add a bit more room to the left and right of the legend box
    axis.text = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.36, 0.91),
    legend.box = "horizontal", # side-by-side legends
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(margin = margin(l = 25, r = 0, t = 30, b = -30), size = 22), # push title into plot
    panel.border = element_blank()
  ) +
  guides(
    fill = guide_coloursteps(
      barwidth = unit(15, "cm"),
      show.limits = TRUE,
      even.steps = TRUE,
      title.position = "top",
      direction = "horizontal"
    )
  )


sens_slope_map_aus

ggsave(
  filename = "./Graphs/Figures/sens_slope_map_2010_2021.pdf",
  plot = sens_slope_map_aus,
  device = "pdf",
  width = 232,
  height = 200,
  units = "mm"
)














#