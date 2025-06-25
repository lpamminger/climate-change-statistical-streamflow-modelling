## Fitting streamflow models to catchments
cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, truncnorm, sloop, patchwork, ozmaps, sf, patchwork, metR, ggmagnify, trend)



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
  mutate(year = as.integer(year))

gauge_information <- readr::read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE
)

lat_lon_gauge_info <- gauge_information |> 
  select(gauge, lat, lon, state)


parameter_results <- read_csv(
  "./Results/CMAES/cmaes_parameter_results.csv", 
  show_col_types = FALSE
)

best_CO2_non_CO2_per_gauge <- read_csv(
  "./Results/CMAES/best_CO2_non_CO2_per_catchment_CMAES.csv",
  show_col_types = FALSE
) 

streamflow_results <- read_csv(
  "./Results/CMAES/cmaes_streamflow_results.csv", 
  show_col_types = FALSE
)




## Import functions ============================================================
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

(nrow(non_equivalent_models) / nrow(compare_equivalent_models)) * 100


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
# 2. Compare the difference (CO2_on - CO2_off) / number of years 

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
  ) |> 
  # add decade column
  mutate(
    decade = year - (year %% 10)
  )
  

decade_a3_on_off_difference_data <-  a3_on_off_difference_data |> 
  summarise(
    sum_realspace_a3_off_streamflow = sum(realspace_a3_off_streamflow),
    sum_realspace_a3_on_streamflow = sum(realspace_a3_on_streamflow),
    sum_precipitation = sum(precipitation),
    n = n(),
    .by = c(gauge, decade)
  ) |> 
  # average by year
  mutate(
    a3_on_off_difference = sum_realspace_a3_on_streamflow - sum_realspace_a3_off_streamflow,
  ) |> 
  mutate(
    by_year_a3_on_off_difference = a3_on_off_difference / n
  )


## Calculate evidence ratio ====================================================
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
  arrange(evidence_ratio)


## Tidy evidence ratio data for plotting =======================================
lat_lon_gauge <- gauge_information |> 
  select(gauge, lat, lon)

lat_long_evidence_ratio <- evidence_ratio_calc |>
  select(!c(AIC_difference)) |>
  left_join(
    lat_lon_gauge,
    by = join_by(gauge)
  )

### Add qualitative labels instead of using numerical evidence ratio ###########
state_gauge <- gauge_information |> 
  select(gauge, state)

binned_lat_lon_evidence_ratio <- lat_long_evidence_ratio |>
  mutate(
    binned_evidence_ratio = case_when(
      between(evidence_ratio, -1E1, 1E1) ~ "Weak",
      between(evidence_ratio, 1E1, 1E2) ~ "Moderate",
      between(evidence_ratio, 1E2, 1E3) ~ "Moderately Strong",
      between(evidence_ratio, 1E3, 1E4) ~ "Strong",
      between(evidence_ratio, 1E4, 1E6) ~ "Very Strong",
      between(evidence_ratio, 1E6, Inf) ~ "Extremely Strong",
      .default = NA
    )
  ) |>
  left_join(
    state_gauge,
    by = join_by(gauge)
  ) |>
  mutate(
    binned_evidence_ratio = factor(
      binned_evidence_ratio,
      levels = c("Weak", "Moderate", "Moderately Strong", "Strong", "Very Strong", "Extremely Strong")
    )
  )

gauges_with_evi_greater_than_moderate <- binned_lat_lon_evidence_ratio |> 
  #filter(!binned_evidence_ratio %in% c("Weak", "Moderate")) |> 
  pull(gauge) |> 
  unique()


# Get decade_a3_on_off_difference_data ready for plotting
plot_ready_decade_a3_on_off_difference_data <- decade_a3_on_off_difference_data |> 
  left_join(
    lat_lon_gauge_info,
    by = join_by(gauge)
  ) |> 
  filter(gauge %in% gauges_with_evi_greater_than_moderate) |> 
  arrange(desc(by_year_a3_on_off_difference))


# Hard code custom breaks
plot_ready_decade_a3_on_off_difference_data |> 
  pull(by_year_a3_on_off_difference) |> 
  range()


## Over entire observation period rather than decade ===========================
year_a3_on_off_difference_data <-  a3_on_off_difference_data |> 
  summarise(
    sum_realspace_a3_off_streamflow = sum(realspace_a3_off_streamflow),
    sum_realspace_a3_on_streamflow = sum(realspace_a3_on_streamflow),
    sum_precipitation = sum(precipitation),
    n = n(),
    .by = gauge
  ) |> 
  # average by year
  mutate(
    a3_on_off_difference = sum_realspace_a3_on_streamflow - sum_realspace_a3_off_streamflow,
  ) |> 
  mutate(
    by_year_a3_on_off_difference = a3_on_off_difference / n
  ) |> 
  # add location data
  left_join(
    lat_lon_gauge_info,
    by = join_by(gauge)
  ) |> 
  filter(gauge %in% gauges_with_evi_greater_than_moderate) |> 
  arrange(desc(by_year_a3_on_off_difference))


# boxplot

# add n = to boxplots
state_count <- year_a3_on_off_difference_data |> 
  summarise(
    n = n(),
    .by = state
  ) |> 
  add_column(
    y_pos = 175
  ) |> 
  mutate(
    label = paste0("n = ", n)
  )

year_a3_on_off_difference_data |> 
  ggplot(aes(x = state, y = by_year_a3_on_off_difference)) +
  geom_boxplot(staplewidth = 0.5) +
  geom_hline(yintercept = 0, colour = "red", linetype = "dashed") +
  geom_label(
    aes(x = state, y = y_pos, label = label),
    data = state_count
  ) +
  theme_bw() +
  labs(
    x = "State",
    y = "Average change in streamflow from CO2 per year (mm)"
  )

year_a3_on_off_difference_data |> 
  summarise(
    median = median(by_year_a3_on_off_difference),
    .by = state
  )


# Get shapefiles for Australia ------------------------------------------------
aus_map <- ozmaps::ozmap(x = "states") |>
  filter(!NAME %in% c("Other Territories")) |>
  rename(state = NAME) |>
  mutate(
    state = case_when(
      state == "New South Wales" ~ "NSW",
      state == "Victoria" ~ "VIC",
      state == "Queensland" ~ "QLD",
      state == "South Australia" ~ "SA",
      state == "Western Australia" ~ "WA",
      state == "Tasmania" ~ "TAS",
      state == "Northern Territory" ~ "NT",
      state == "Australian Capital Territory" ~ "ACT",
    )
  )

combine_NSW_ACT <- aus_map |>
  filter(state %in% c("NSW", "ACT")) |>
  st_union()

aus_map[1, 2] <- list(combine_NSW_ACT)

aus_map <- aus_map |>
  filter(state != "ACT")




## Plot map percentage difference ==============================================

### REPLACE WITH FULL FAT ###
# 1990 vs. 2010

# 1990
# Make final plot --------------------------------------------------------------

### Custom colour palette 
big_palette <- function(x) {
  c("#67001f",
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
    "#053061")
}


decade_comparison_CO2_impact <- function(decade) {
  ## Generate Insets =============================================================
  ### Filter data by state #######################################################
  
  QLD_data <- plot_ready_decade_a3_on_off_difference_data |>
    filter(state == "QLD") |> 
    filter(decade == {{ decade }})
  
  NSW_data <- plot_ready_decade_a3_on_off_difference_data |>
    filter(state == "NSW") |> 
    filter(decade == {{ decade }})
  
  VIC_data <- plot_ready_decade_a3_on_off_difference_data |>
    filter(state == "VIC") |> 
    filter(decade == {{ decade }})
  
  WA_data <- plot_ready_decade_a3_on_off_difference_data |>
    filter(state == "WA") |> 
    filter(decade == {{ decade }})
  
  TAS_data <- plot_ready_decade_a3_on_off_difference_data |>
    filter(state == "TAS") |> 
    filter(decade == {{ decade }})
  
  
  #plot_ready_decade_a3_on_off_difference_data |> 
  # pull(by_year_a3_on_off_difference) |> 
  # range()
  custom_limits <- c(-800, 500)
  custom_breaks <- c(-400, -100, -50, -10, -1, 0, 1, 10, 50, 100, 250)
  
  
  ### Generate inset plots #######################################################
  
  inset_plot_QLD <- aus_map |>
    filter(state == "QLD") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = QLD_data,
      aes(x = lon, y = lat, fill = by_year_a3_on_off_difference),
      show.legend = FALSE,
      size = 2.5,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = custom_breaks, # range()
      limits = custom_limits,
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
      aes(x = lon, y = lat, fill = by_year_a3_on_off_difference),
      show.legend = FALSE,
      size = 2.5,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = custom_breaks,
      limits = custom_limits,
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
      aes(x = lon, y = lat, fill = by_year_a3_on_off_difference),
      show.legend = FALSE,
      size = 2.5,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = custom_breaks, # range()
      limits = custom_limits,
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
      aes(x = lon, y = lat, fill = by_year_a3_on_off_difference),
      show.legend = FALSE,
      size = 2.5,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = custom_breaks,
      limits = custom_limits,
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
      aes(x = lon, y = lat, fill = by_year_a3_on_off_difference),
      show.legend = FALSE,
      size = 2.5,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = custom_breaks, # range()
      limits = custom_limits,
      show.limits = TRUE, 
      guide = "colorsteps"
    ) +
    theme_void()
  
  
  
  ## Put it together =============================================================
  
  single_map_aus <- aus_map |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = plot_ready_decade_a3_on_off_difference_data |> filter(decade == {{ decade }}),
      mapping = aes(x = lon, y = lat, fill = by_year_a3_on_off_difference),
      size = 3,
      colour = "black",
      shape = 21,
      inherit.aes = FALSE,
      stroke = 0.1
    ) +
    theme_bw() +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = custom_breaks,
      limits = custom_limits,
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
      #aes(from = state == "VIC"), # use aes rather than manually selecting area
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
      x = NULL,#"Latitude",
      y = NULL,#"Longitude",
      fill = "Average impact of CO2 on streamflow per year (mm/year)",
      title = paste0(decade)
    ) +
    theme(
      legend.key = element_rect(fill = "grey80"),
      legend.title = element_text(hjust = 0.5),
      #legend.background = element_rect(colour = "black"),
      axis.text = element_blank(), 
      legend.position = "inside",
      legend.position.inside = c(0.351, 0.9),
      legend.box = "horizontal", # side-by-side legends
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(margin = margin(l = 25, r = 0, t = 30, b = -30), size = 22) # push title into plot
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


# Patchwork together
patchwork_CO2_impact_decade <- (decade_comparison_CO2_impact(1990) | decade_comparison_CO2_impact(2020)) + plot_layout(guides = "collect") & theme(legend.position='bottom')
patchwork_CO2_impact_decade

ggsave(
  filename = "./Graphs/Figures/log_sinh_CO2_on_off_mm_1990_vs_2020.pdf",
  plot = patchwork_CO2_impact_decade,
  device = "pdf",
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
  #geom_line(
  #  aes(x = year, y = precipitation), 
  #  colour = "black", 
  #  linetype = "dashed",
  #  inherit.aes = FALSE
  #  ) +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~gauge, scales = "free_y")

ggsave(
  filename = "mega_timeseries_plot.pdf",
  device = "pdf",
  plot = mega_timeseries_plot,
  width = 1189,
  height = 841,
  units = 'mm'
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
  filename = "mega_rainfall_runoff_plot.pdf",
  device = "pdf",
  plot = mega_rainfall_runoff_plot,
  width = 1189,
  height = 841,
  units = 'mm'
)


# Other ideas ------------------------------------------------------------------
# Average not by decade?



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
  #filter(gauge %in% c("112101B", "405228", "610001", "116011A")) |> 
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






realspace_t_series <- altered_realspace_streamflow_data_a3_off |> 
  filter(gauge == "207015") |> 
  left_join(
    observed_streamflow,
    by = join_by(gauge, precipitation)
  ) |> 
  select(year, precipitation, realspace_a3_on, realspace_a3_off, q_mm) |>
  rename(
    observed_streamflow = q_mm,
    modelled_streamflow_CO2_on = realspace_a3_on,
    modelled_streamflow_CO2_off = realspace_a3_off
  ) |> 
  pivot_longer(
    contains("streamflow"),
    names_to = "streamflow_type",
    values_to = "streamflow"
  ) |> 
  ggplot(aes(x = year, y = streamflow, colour = streamflow_type)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2, alpha = 0.9) +
  geom_line(aes(x = year, y = precipitation), colour = "black", linewidth = 0.8, linetype = "dashed") +
  labs(
    x = "Year",
    y = "Streamflow (mm)",
    colour = NULL,
    title = "Box-cox transform"
  ) +
  scale_y_continuous(
    name = "Streamflow (mm)",
    sec.axis = sec_axis(~.*1, name = "Precipitation (mm)")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

realspace_t_series

test_plot <- boxcox_t_series / realspace_t_series

ggsave(
  filename = "test_plot.pdf",
  plot = test_plot,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)














# Rate of change CO2-on vs CO2-off using sens slope ----------------------------
rate_of_change <- altered_realspace_streamflow_data_a3_off |> 
  select(gauge, year, realspace_a3_on, realspace_a3_off) |> # real space
  #select(gauge, year, a3_off_modelled_boxcox_streamflow, modelled_boxcox_streamflow) |> # boxcox
  mutate(
    relative_difference = realspace_a3_on - realspace_a3_off, # real space
    #relative_difference = modelled_boxcox_streamflow - a3_off_modelled_boxcox_streamflow
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
  
  #return(d)
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
    lat_lon_gauge_info,
    by = join_by(gauge)
  ) 

# Put the sens_slope on a map mm/year ------------------------------------------
sens_slope_palette <- function(x) {
  c("#67001f",
    "#b2182b",
    "#d6604d",
    "#f4a582",
    "#fddbc7",
    "white",
    "#d1e5f0",
    "#92c5de",
    "#4393c3",
    "#2166ac",
    "#053061")
}

gauge_rate_of_change |> pull(sens_slope) |> quantile()

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
    breaks = c(-5, -2.5, -1, -0.5, -0.001, 0.001, 0.5, 1, 2.5, 5),# range()
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
    #aes(from = state == "VIC"), # use aes rather than manually selecting area
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
    x = NULL,#"Latitude",
    y = NULL,#"Longitude",
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














# Compare evidence ratio for best non-CO2 vs. CO2 model ------------------------
# Method:
# 1. streamflow_data_best_CO2_non_CO2
# 2. Find percentage difference for each gauge best_Co2_model - best_non_CO2_model / best_CO2_model
# 3. summarise percentage difference
# 4. Join with evidence ratio
# 5. plot

## Evidence ratio for joining ==================================================
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
  arrange(evidence_ratio)



## percentage difference calculation ===========================================
calculate_percent_diff_best_non_CO2_and_CO2 <- function(gauge, best_CO2_non_CO2_data) {
  
  # Pivot wider for mutate
  wide_streamflow_data <- best_CO2_non_CO2_data |> 
    select(!c(objective_function, optimiser)) |> 
    filter(gauge == {{ gauge }}) |> 
    pivot_wider(
      names_from = streamflow_model,
      values_from = modelled_boxcox_streamflow
    ) 
  
  # Rename columns for generic percent diff
  column_names <- names(wide_streamflow_data)
  CO2_model_name_index <- str_detect(names(wide_streamflow_data), "streamflow_model") & str_detect(names(wide_streamflow_data), "CO2")
  no_CO2_model_name_index <- str_detect(names(wide_streamflow_data), "streamflow_model") & !str_detect(names(wide_streamflow_data), "CO2")
  column_names[CO2_model_name_index] <- "CO2_model"
  column_names[no_CO2_model_name_index] <- "non_CO2_model"
  
  # Percent diff calculation - sum everything then percent diff
  calculation_streamflow_data <- wide_streamflow_data |> 
    `names<-`(column_names) |>  
    summarise(
      sum_CO2_model = sum(CO2_model),
      sum_non_CO2_model = sum(non_CO2_model),
      .by = gauge
    ) |> 
    mutate(
      percent_diff = ((sum_CO2_model - sum_non_CO2_model) / sum_CO2_model) * 100
    )
  
  return(calculation_streamflow_data)
}


percent_diff_streamflow_data_best_CO2_non_CO2 <- map(
  .x = streamflow_data_best_CO2_non_CO2 |> pull(gauge) |> unique(),
  .f = calculate_percent_diff_best_non_CO2_and_CO2,
  best_CO2_non_CO2_data = streamflow_data_best_CO2_non_CO2
) |> 
  list_rbind()


## Summarise percentage difference and join evidence ratio =====================
summarise_percent_diff_streamflow_data_best_CO2_non_CO2 <- percent_diff_streamflow_data_best_CO2_non_CO2 |> 
  left_join(
    evidence_ratio_calc,
    by = join_by(gauge)
  )
  
  

## Plot percentage_diff vs. evidence ratio =====================================
### 

summarise_percent_diff_streamflow_data_best_CO2_non_CO2 |> 
  ggplot(aes(x = evidence_ratio, y = percent_diff)) +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(
    x = "Evidence Ratio (log10)",
    y = "Mean Percentage Difference in Modelled Streamflow"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )




# Compare non-CO2 vs. CO2 and a3-on vs a3-off ----------------------------------

# Would be a good idea to get streamflow_data_a3_off and streamflow_data_best_CO2_non_CO2
# in the same format

# suggested format?

gauge_key <- "121003A"  #"A2390531" #"226410" #"403226"

## Rainfall-runoff relationship
### CO2 vs. non CO2
test_1 <- streamflow_data_best_CO2_non_CO2 |> 
  filter(gauge == {{ gauge_key }}) |> 
  select(!c(objective_function, optimiser)) |> 
  mutate(
    is_CO2_model = str_detect(streamflow_model, "CO2")
  ) |> 
  pivot_longer(
    cols = ends_with("streamflow"),
    names_to = "modelled_or_observed",
    values_to = "streamflow"
  ) |>
  mutate(
    modelled_or_observed = if_else(is_CO2_model, paste0("CO2_", modelled_or_observed), paste0("no_CO2_", modelled_or_observed))
  ) |> 
  mutate(
    modelled_or_observed = if_else(str_detect(modelled_or_observed, "observed"), "observed_boxcox_streamflow", modelled_or_observed)
  )
  # remove duplicated observed values
  #filter(!is_CO2_model & modelled_or_observed != "modelled_boxcox_streamflow")
  

graph_CO2_vs_non_CO2 <- test_1 |> 
  mutate(
    modelled_or_observed = case_when(
      modelled_or_observed == "observed_boxcox_streamflow" ~ "Observed Box-Cox Streamflow",
      modelled_or_observed == "CO2_modelled_boxcox_streamflow" ~ "Best CO2 Modelled Box-Cox Streamflow",
      modelled_or_observed == "no_CO2_modelled_boxcox_streamflow" ~ "Best non-CO2 Modelled Box-Cox Streamflow",
      .default = NA
    )
  ) |> 
  ggplot(aes(x = precipitation, y = streamflow, colour = modelled_or_observed)) +
  geom_point(size = 2) +
  geom_smooth(
    formula = y ~ x,
    method = lm,
    se = FALSE
  ) +
  labs(
    x = "Precipitation (mm/year)",
    y = "Box-Cox Streamflow (mm/year)",
    colour = NULL,
    title = "Best CO2 modelled vs. best non-CO2 model"
  ) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.35, 0.9),
    legend.background = element_rect(colour = "black"),
    plot.title = element_text(hjust = 0.5)
  )
  

timeseries_graph_CO2_vs_non_CO2 <- test_1 |> 
  ggplot(aes(x = year, y = streamflow, colour = modelled_or_observed)) +
  geom_point(size = 2) +
  geom_line() +
  labs(
    x = "Year",
    y = "Box-Cox Streamflow",
    colour = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.9, 0.9)
  )

timeseries_graph_CO2_vs_non_CO2

### a3 off vs. a3 on
test_2 <- streamflow_data_a3_off |> 
  filter(gauge == {{ gauge_key }}) |> 
  pivot_longer(
    cols = ends_with("streamflow"),
    names_to = "modelled_or_observed",
    values_to = "boxcox_streamflow"
  )

graph_CO2_off_vs_CO2_on <- test_2 |> 
  mutate(
    modelled_or_observed = case_when(
      modelled_or_observed == "observed_boxcox_streamflow" ~ "Observed Box-Cox Streamflow",
      modelled_or_observed == "modelled_boxcox_streamflow" ~ "Best CO2 Modelled Box-Cox Streamflow",
      modelled_or_observed == "a3_off_modelled_boxcox_streamflow" ~ "CO2 Component Turned Off Modelled Box-Cox Streamflow",
      .default = NA
    )
  ) |> 
  ggplot(aes(x = precipitation, y = boxcox_streamflow, colour = modelled_or_observed)) +
  geom_point(size = 2) +
  geom_smooth(
    formula = y ~ x,
    method = lm,
    se = FALSE
  ) +
  theme_bw() +
  labs(
    x = "Precipitation (mm/year)",
    y = "Box-Cox Streamflow (mm/year)",
    colour = NULL,
    title = "Turning CO2 Component Off"
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.4, 0.9),
    legend.background = element_rect(colour = "black"),
    plot.title = element_text(hjust = 0.5)
  )


save_rainfall_runoff_test_plot <- graph_CO2_vs_non_CO2 + graph_CO2_off_vs_CO2_on
ggsave(
  filename = "rainfall_runoff_comparison.pdf",
  plot = save_rainfall_runoff_test_plot,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)






# Compare 3 methods to estimate the impact of CO2 on streamflow ----------------
## Rank all the models by gauge (CMAES results)
ranked_models_per_gauge <- CMAES_results |>
  group_by(gauge) |>
  arrange(AIC) |>
  ungroup()

## Comparing the three methods
### Method 1: set a3 to zero
### Method 2: compare best CO2 to best non-CO2
### Method 3: compare best CO2 to equivalent non-CO2 model (could also be best non-CO2)

### How to compare the methods:
### 1. streamflow-time with observed
### 2. rainfall-runoff with observed
### 3. Difference between models vs. evidence ratio - not sure if this will work
###    with method 1 - I could do difference between observed - modelled vs. evi ratio

## Method 2 ====================================================================


### Difference in streamflow models ############################################
### Difference in streamflow models per year mean(mm/y) vs. evidence ratio or AIC
### Use best_CO2_non_CO2_per_gauge to filter streamflow_data
### Pivot wider, mutate model 1 - model 2 -> summarise mean .by gauge
### y axis-options:
# mean absolute difference in streamflow (mm/y)
# total difference in streamflow over the entire period
# mean percentage difference in streamflow per year




# percentage_diff <- function(x, y) {
#  (abs(x - y) / ((x + y) / 2)) * 100
# }

### Difference in streamflow for the models
difference_streamflow_per_year_best_CO2_non_CO2 <- streamflow_data_best_CO2_non_CO2 |>
  # Rename streamflow models to non-CO2 and CO2
  mutate(
    CO2_or_non_CO2 = if_else(str_detect(streamflow_model, "CO2"), "CO2_streamflow", "non_CO2_streamflow")
  ) |>
  # Remove streamflow model to make CO2 and non-CO2 on same row
  select(!c(streamflow_model, objective_function, optimiser)) |>
  pivot_wider(
    names_from = CO2_or_non_CO2,
    values_from = modelled_boxcox_streamflow
  ) |>
  mutate(
    yearly_CO2_non_CO2_difference = CO2_streamflow - non_CO2_streamflow,
    percentage_yearly_CO2_non_CO2_difference = ((CO2_streamflow - non_CO2_streamflow) / CO2_streamflow) * 100 # ,
    # alternative_percentage = percentage_diff(CO2_streamflow, non_CO2_streamflow)
  )

### Summarise yearly differences
summary_streamflow_best_CO2_non_CO2 <- difference_streamflow_per_year_best_CO2_non_CO2 |>
  summarise(
    mean_yearly_CO2_non_CO2_difference = mean(yearly_CO2_non_CO2_difference),
    sum_yearly_CO2_non_CO2_difference = sum(yearly_CO2_non_CO2_difference),
    mean_percent_yearly_CO2_non_CO2_difference = mean(percentage_yearly_CO2_non_CO2_difference),
    sum_CO2_streamflow = sum(CO2_streamflow),
    sum_non_CO2_streamflow = sum(non_CO2_streamflow),
    .by = gauge
  ) |>
  # add evidence ratio
  left_join(
    evidence_ratio_calc,
    by = join_by(gauge)
  ) |>
  # percentage diff

  arrange(sum_yearly_CO2_non_CO2_difference)

### Best CO2 vs. non-CO2 model comparison ######################################
summary_streamflow_best_CO2_non_CO2 |>
  filter(evidence_ratio > 0) |>
  # filter(abs(mean_percent_yearly_CO2_non_CO2_difference) < 1) |>
  ggplot(aes(x = evidence_ratio, y = mean_percent_yearly_CO2_non_CO2_difference)) +
  geom_point() +
  scale_x_log10() +
  theme_bw() +
  labs(
    y = "mean(CO2_t - non-CO2_t / CO2_t)"
  )

# This graph tells us
# x-axis is evidence ratio between CO2 and non-CO2 (log10 scaled)
# x-axis values < 0 removed for log10 scale. This removes catchments
# where the non-CO2 model is better
# the CO2 model on average produces -0.17 % less streamflow per year
# compared to the non-CO2 models
# There is no relationship between evidence ratio and difference in models
# We expected to see a larger difference in CO2 and non-CO2 models as
# the evidence ratio grew.

# We really need to compare it to the observed...

summary_streamflow_best_CO2_non_CO2 |>
  pull(mean_percent_yearly_CO2_non_CO2_difference) |>
  quantile()