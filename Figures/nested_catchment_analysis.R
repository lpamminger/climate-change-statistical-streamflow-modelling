# Evidence ratio analysis for nested catchments

# Figures produced in this R file ----------------------------------------------


# CODE


# Import libraries -------------------------------------------------------------
pacman::p_load(tidyverse, ozmaps, sf, ggmagnify, trend, patchwork)


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

lat_lon_gauge <- gauge_information |>
  select(gauge, lat, lon, state, status, parent_catchment)


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


start_stop_indexes <- readr::read_csv(
  "./Data/Tidy/start_end_index.csv",
  show_col_types = FALSE
)



# Nested catchment statistics --------------------------------------------------
subcatchment_information <- gauge_information |> 
  filter(status == "sub_catchment") |> 
  count(state)

total_catchments <- gauge_information |> 
  count(state) |> 
  rename(total = n)

subcatchment_stats <- subcatchment_information |> 
  left_join(
    total_catchments,
    by = join_by(state)
  ) |> 
  mutate(
    percent = (n / total) * 100
  )

subcatchment_stats |> 
  mutate(
    total_label = paste0(n, "/", total)
  ) |> 
  ggplot(aes(x = state, y = percent, label = total_label)) +
  geom_col() +
  geom_text(vjust = 0) +
  labs(
    x = "State",
    y = "Nested Catchments (%)"
  ) +
  theme_bw()


# Calculate evidence ratio -----------------------------------------------------
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

lat_long_evidence_ratio <- evidence_ratio_calc |>
  select(!c(AIC_difference)) |>
  left_join(
    lat_lon_gauge,
    by = join_by(gauge)
  )

### Add qualitative labels instead of using numerical evidence ratio ###########

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
  mutate(
    binned_evidence_ratio = factor(
      binned_evidence_ratio,
      levels = c("Weak", "Moderate", "Moderately Strong", "Strong", "Very Strong", "Extremely Strong")
    )
  )


### Add direction of change and whether the slope/intercept changed ############
best_model_per_gauge <- best_CO2_non_CO2_per_gauge |>
  slice_min(
    AIC,
    by = gauge
  ) |>
  distinct()


direction_of_a3_change <- best_model_per_gauge |>
  filter(parameter %in% c("a3_intercept", "a3_slope")) |>
  mutate(
    intercept_or_slope = if_else(str_detect(streamflow_model, "intercept"), "Intercept", "Slope")
  ) |>
  select(gauge, streamflow_model, parameter, parameter_value, intercept_or_slope) |>
  mutate(
    CO2_direction = if_else(parameter_value < 0, "Negative", "Positive")
  ) |>
  select(gauge, CO2_direction, intercept_or_slope)


a3_direction_binned_lat_lon_evidence_ratio <- binned_lat_lon_evidence_ratio |>
  left_join(
    direction_of_a3_change,
    by = join_by(gauge)
  ) |>
  replace_na(list(CO2_direction = "No CO2 Term", intercept_or_slope = "No CO2 Term")) |>
  unite(
    col = "impact_of_CO2_term",
    CO2_direction,
    intercept_or_slope,
    sep = "-"
  ) |>
  mutate(
    impact_of_CO2_term = if_else(impact_of_CO2_term == "No CO2 Term-No CO2 Term", "No CO2 Term", impact_of_CO2_term)
  ) |>
  mutate(
    impact_of_CO2_term = factor(
      impact_of_CO2_term,
      levels = c("No CO2 Term", "Negative-Intercept", "Positive-Intercept", "Negative-Slope", "Positive-Slope")
    )
  ) |>  
  # what would the results look like if subcatchment were not allowed
  filter(status != "no_nested") |> #the results roughly the same 
  group_by(parent_catchment) |> 
  mutate(group_number = cur_group_id()) |> 
  ungroup() |> 
  arrange(group_number)

remove_catchments_without_parents <- a3_direction_binned_lat_lon_evidence_ratio |> 
  count(parent_catchment) |> 
  arrange(n) |> 
  left_join(
    a3_direction_binned_lat_lon_evidence_ratio,
    by = join_by(parent_catchment)
  ) |> 
  # there are parent catchments that are not present in the data - remove
  filter(n < 2) |> 
  pull(group_number)

a3_direction_binned_lat_lon_evidence_ratio <- a3_direction_binned_lat_lon_evidence_ratio |> 
  filter(!group_number %in% remove_catchments_without_parents)

# Make final plot --------------------------------------------------------------

## Load shapefiles
aus_map <- generate_aus_map_sf()

major_streamflow_network_sf <- st_read(
  "./Data/Maps/filtered_streamflow_network_major_sf"
  )

minor_streamflow_network_sf <- st_read(
  "./Data/Mapsfiltered_streamflow_network_minor_sf"
  )

# only include minor on the inset plots - make the linewidth < major

boundary_sf <- st_read(
  "./Data/Maps/filtered_catchment_boundaries_sf"
  )

## Custom colour palette
custom_palette <- function(x) {
  rev(c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#f7f7f7"))
}


## Generate Insets =============================================================
### Filter data by state #######################################################

# this is a bad way of doing this --> use map
QLD_data <- a3_direction_binned_lat_lon_evidence_ratio |>
  filter(state == "QLD")
QLD_major_network_sf <- major_streamflow_network_sf |> 
  filter(state == "QLD")
QLD_minor_network_sf <- minor_streamflow_network_sf |> 
  filter(state == "QLD")
QLD_boundaries <- boundary_sf |> 
  filter(state == "QLD")


NSW_data <- a3_direction_binned_lat_lon_evidence_ratio |>
  filter(state == "NSW")
NSW_major_network_sf <- major_streamflow_network_sf |> 
  filter(state == "NSW")
NSW_minor_network_sf <- minor_streamflow_network_sf |> 
  filter(state == "NSW")
NSW_boundaries <- boundary_sf |> 
  filter(state == "NSW")


VIC_data <- a3_direction_binned_lat_lon_evidence_ratio |>
  filter(state == "VIC")
VIC_major_network_sf <- major_streamflow_network_sf |> 
  filter(state == "VIC")
VIC_minor_network_sf <- minor_streamflow_network_sf |> 
  filter(state == "VIC")
VIC_boundaries <- boundary_sf |> 
  filter(state == "VIC")


WA_data <- a3_direction_binned_lat_lon_evidence_ratio |>
  filter(state == "WA")
WA_major_network_sf <- major_streamflow_network_sf |> 
  filter(state == "WA")
WA_minor_network_sf <- minor_streamflow_network_sf |> 
  filter(state == "WA")
WA_boundaries <- boundary_sf |> 
  filter(state == "WA")


TAS_data <- a3_direction_binned_lat_lon_evidence_ratio |>
  filter(state == "TAS")
TAS_major_network_sf <- major_streamflow_network_sf |> 
  filter(state == "TAS")
TAS_minor_network_sf <- minor_streamflow_network_sf |> 
  filter(state == "TAS")
TAS_boundaries <- boundary_sf |> 
  filter(state == "TAS")







### Generate inset plots #######################################################

inset_plot_QLD <- aus_map |>
  filter(state == "QLD") |>
  ggplot() +
  geom_sf() +
  geom_sf(
    data = QLD_boundaries
  ) +
  geom_point(
    data = QLD_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1
  ) +
  geom_sf(
    data = QLD_major_network_sf,
    colour = "blue",
    linewidth = 0.1
  ) +
  geom_sf(
    data = QLD_minor_network_sf,
    colour = "blue",
    linewidth = 0.02
  ) +
  geom_text(
    data = QLD_data,
    aes(x = lon, y = lat, label = group_number),
    show.legend = FALSE,
    size = 1.2
  ) +
  scale_fill_manual(
    values = custom_palette(),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(21, 22, 23, 25, 24),
    drop = FALSE # Very important - every plots has the same discrete values
  ) +
  theme_void()


inset_plot_NSW <- aus_map |>
  filter(state == "NSW") |>
  ggplot() +
  geom_sf() +
  geom_sf(
    data = NSW_boundaries
  ) +
  geom_point(
    data = NSW_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1
  ) +
  geom_sf(
    data = NSW_major_network_sf,
    colour = "blue",
    linewidth = 0.1
  ) +
  geom_sf(
    data = NSW_minor_network_sf,
    colour = "blue",
    linewidth = 0.02
  ) +
  geom_text(
    data = NSW_data,
    aes(x = lon, y = lat, label = group_number),
    show.legend = FALSE,
    size = 1.2
  ) +
  scale_fill_manual(
    values = custom_palette(),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(21, 22, 23, 25, 24),
    drop = FALSE
  ) +
  theme_void()


inset_plot_VIC <- aus_map |>
  filter(state == "VIC") |>
  ggplot() +
  geom_sf() +
  geom_sf(
    data = VIC_boundaries
  ) +
  geom_point(
    data = VIC_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1
  ) +
  geom_sf(
    data = VIC_major_network_sf,
    colour = "blue",
    linewidth = 0.1
  ) +
  geom_sf(
    data = VIC_minor_network_sf,
    colour = "blue",
    linewidth = 0.02
  ) +
  geom_text(
    data = VIC_data,
    aes(x = lon, y = lat, label = group_number),
    show.legend = FALSE,
    size = 1.2
  ) +
  scale_fill_manual(
    values = custom_palette(),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(21, 22, 23, 25, 24),
    drop = FALSE
  ) +
  theme_void()


inset_plot_WA <- aus_map |>
  filter(state == "WA") |>
  ggplot() +
  geom_sf() +
  geom_sf(
    data = WA_boundaries
  ) +
  geom_point(
    data = WA_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1
  ) +
  geom_sf(
    data = WA_network_sf,
    colour = "blue",
    linewidth = 0.1
  ) +
  geom_sf(
    data = WA_minor_network_sf,
    colour = "blue",
    linewidth = 0.02
  ) +
  geom_text(
    data = WA_data,
    aes(x = lon, y = lat, label = group_number),
    show.legend = FALSE,
    size = 1.2
  ) +
  scale_fill_manual(
    values = custom_palette(),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(21, 22, 23, 25, 24),
    drop = FALSE
  ) +
  theme_void()


inset_plot_TAS <- aus_map |>
  filter(state == "TAS") |>
  ggplot() +
  geom_sf() +
  geom_sf(
    data = TAS_boundaries
  ) +
  geom_point(
    data = TAS_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1,
  ) +
  geom_sf(
    data = TAS_network_sf,
    colour = "blue",
    linewidth = 0.1
  ) +
  geom_sf(
    data = TAS_minor_network_sf,
    colour = "blue",
    linewidth = 0.02
  ) +
  geom_text(
    data = TAS_data,
    aes(x = lon, y = lat, label = group_number),
    show.legend = FALSE,
    size = 1.2
  ) +
  scale_fill_manual(
    values = custom_palette(),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(21, 22, 23, 25, 24),
    drop = FALSE
  ) +
  theme_void()


## Put it together =============================================================

single_map_aus <- aus_map |>
  ggplot() +
  geom_sf() +
  geom_sf(
    data = boundary_sf
  ) +
  geom_point(
    data = a3_direction_binned_lat_lon_evidence_ratio,
    mapping = aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
    size = 3,
    colour = "black",
    stroke = 0.1
  ) +
  geom_sf(
    data = major_streamflow_network_sf,
    colour = "blue",
    linewidth = 0.1
  ) +
  geom_text(
    data = a3_direction_binned_lat_lon_evidence_ratio,
    aes(x = lon, y = lat, label = group_number),
    show.legend = FALSE,
    size = 1.2
  ) +
  theme_bw() +
  scale_fill_manual(
    values = custom_palette(),
    drop = FALSE
  ) +
  scale_shape_manual(
    labels = c(
      bquote("No" ~ CO[2] ~ "Term"),
      "Negative-Intercept",
      "Positive-Intercept",
      "Negative-Slope",
      "Positive-Slope"
    ),
    values = c(21, 22, 23, 25, 24),
    drop = FALSE
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
    fill = "Evidence Ratio",
    shape = bquote("Impact of" ~ CO[2] ~ "Term")
  ) +
  theme(
    # legend.key = element_rect(fill = "grey80"),
    legend.title = element_text(hjust = 0.5),
    legend.background = element_rect(colour = "black"),
    axis.text = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.346, 0.9), # constants used to move the legend in the right place
    legend.box = "horizontal", # side-by-side legends
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank()
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 5, shape = 21), nrow = 3), # Wrap legend with nrow
    shape = guide_legend(override.aes = list(size = 5, fill = "grey50"), nrow = 3)
  )


ggsave(
  filename = "./Figures/Main/TESTING_subcatchment_only_evidence_ratio_aus_map.pdf",
  plot = single_map_aus,
  device = "pdf",
  width = 232,
  height = 200, # 210,
  units = "mm"
)









# Percentage changes -----------------------------------------------------------
## Summarise modelled streamflow data ==========================================
only_gauge_model_best_CO2_non_CO2_per_gauge <- best_CO2_non_CO2_per_gauge |>
  select(gauge, streamflow_model) |>
  distinct()



streamflow_data_best_CO2_non_CO2 <- streamflow_results |>
  semi_join(
    only_gauge_model_best_CO2_non_CO2_per_gauge,
    by = join_by(gauge, streamflow_model)
  )



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





### Calculate limits and breaks ################################################
make_limits <- function(timeseries) {
  # round up to next whole number
  limits <- timeseries |> range()
  sign_limits <- sign(limits)
  
  sign_limits * ceiling(abs(limits))
}







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
    geom_text(
      data = QLD_data,
      aes(x = lon, y = lat, label = group_number),
      show.legend = FALSE,
      size = 1.2
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
    geom_text(
      data = NSW_data,
      aes(x = lon, y = lat, label = group_number),
      show.legend = FALSE,
      size = 1.2
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
    geom_text(
      data = VIC_data,
      aes(x = lon, y = lat, label = group_number),
      show.legend = FALSE,
      size = 1.2
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
    geom_text(
      data = WA_data,
      aes(x = lon, y = lat, label = group_number),
      show.legend = FALSE,
      size = 1.2
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
    geom_text(
      data = TAS_data,
      aes(x = lon, y = lat, label = group_number),
      show.legend = FALSE,
      size = 1.2
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
    geom_text(
      data = data,
      aes(x = lon, y = lat, label = group_number),
      show.legend = FALSE,
      size = 1.2
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
      plot.title = element_text(margin = margin(l = 25, r = 0, t = 30, b = -30), size = 12, face = "bold") # push title into plot
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
  ) #|>
  # filter evidence ratio
  #filter(evidence_ratio > 100)




### redo limits
CO2_impact_on_streamflow_percent_limits <- plotting_best_CO2_non_CO2_streamflow |>
  drop_na() |> 
  pull(CO2_impact_on_streamflow_percent) |>
  make_limits() |>
  as.double()

hard_coded_breaks_CO2_impact_of_streamflow <- c(-75, -50, -25, -10, -1, 0, 1, 10, 25, 50, 75)



## Plots for 1990s and 2010s ===================================================

### get the group numbers from 
group_numbers <- a3_direction_binned_lat_lon_evidence_ratio |> 
  select(gauge, group_number)

percentage_difference_CO2_model_non_CO2_model_1990s <- plotting_best_CO2_non_CO2_streamflow |>
  filter(decade == 1) |> 
  ## Only include subcatchments
  filter(status != "no_nested") |>
  ## use the same subcatchment numbers for all plots
  left_join(
    group_numbers,
    by = join_by(gauge)
  ) |> 
  drop_na() # remove NA - due to gauges without parents
  
percentage_difference_CO2_model_non_CO2_model_2010s <- plotting_best_CO2_non_CO2_streamflow |>
  filter(decade == 2) |>  
  ## Only include subcatchments
  filter(status != "no_nested") |>
  ## use the same subcatchment numbers for all plots
  left_join(
    group_numbers,
    by = join_by(gauge)
  ) |> 
  drop_na() # remove NA - due to gauges without parents



patchwork_CO2_model_and_non_CO2_model_percentage_differences <- (make_CO2_streamflow_percentage_change_map(percentage_difference_CO2_model_non_CO2_model_1990s, "a") | make_CO2_streamflow_percentage_change_map(percentage_difference_CO2_model_non_CO2_model_2010s, "b")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(
  filename = "./Figures/Main/TESTING_subcatchment_streamflow_percentage_difference.pdf",
  plot = patchwork_CO2_model_and_non_CO2_model_percentage_differences,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)



