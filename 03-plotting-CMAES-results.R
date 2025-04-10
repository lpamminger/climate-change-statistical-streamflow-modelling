# Script produces:
# - a map of best CO2 vs best non-CO2 of Australia
# - streamflow timeseries comparing observed, best CO2 and best non-CO2



# TODO - transfer code/graphs from playing_with_graphs
# Key graphs to transfer are:
# 2. Time of emergence histogram
# 3. Map of time of emergence
# 4. Component map

# - account for the different types of CO2 models - slope vs. intercept
# 1. graphically compare (colour gradient evidence ratio, shape pos/neg, outline
#    or stroke colour can be slope or intercept)
# 2. some sort of slope and intercept analysis. Like proportion of pos/neg slope,
#    pos/neg intercept


cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, ozmaps, sf, ggmagnify) # need to update this linux


## Utility functions ===========================================================
source("./Functions/utility.R")


## Import streamflow functions =================================================
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
source("./Functions/boxcox_transforms.R")


# Import data ------------------------------------------------------------------
CMAES_results <- read_csv(
  "./Results/my_cmaes/CMAES_parameter_results_20250331.csv", 
  show_col_types = FALSE
) 

data <- readr::read_csv(
  "./Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  col_types = "icdddddddll", # ensuring each column is a the correct type
  show_col_types = FALSE
)

streamflow_results <- read_csv(
  "Results/my_cmaes/CMAES_streamflow_results_20250331.csv",
  show_col_types = FALSE,
  col_select = !optimiser
)

gauge_information <- read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE,
  col_select = c(gauge, bc_lambda, state, lat, lon)
) # CAMELS is Australia wide.

best_CO2_non_CO2_per_gauge <- read_csv(
  "./Results/my_cmaes/unmodified_best_CO2_non_CO2_per_catchment_CMAES_20250331.csv",
  show_col_types = FALSE
) 
  




# 1. A map of best CO2 vs best non-CO2 of Australia -----------------------------


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
### Add lat and lon ############################################################
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



# Make final plot --------------------------------------------------------------

### Custom colour palette 
custom_palette <- function(x) {
  rev(c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#f7f7f7"))
}


## Generate Insets =============================================================
### Filter data by state #######################################################

QLD_data <- a3_direction_binned_lat_lon_evidence_ratio |>
  filter(state == "QLD")

NSW_data <- a3_direction_binned_lat_lon_evidence_ratio |>
  filter(state == "NSW")

VIC_data <- a3_direction_binned_lat_lon_evidence_ratio |>
  filter(state == "VIC")

WA_data <- a3_direction_binned_lat_lon_evidence_ratio |>
  filter(state == "WA")

TAS_data <- a3_direction_binned_lat_lon_evidence_ratio |>
  filter(state == "TAS")


### Generate inset plots #######################################################

inset_plot_QLD <- aus_map |>
  filter(state == "QLD") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = QLD_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1
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


inset_plot_NSW <- aus_map |>
  filter(state == "NSW") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = NSW_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1
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
  geom_point(
    data = VIC_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1
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
  geom_point(
    data = WA_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1
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
  geom_point(
    data = TAS_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1
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
  geom_point(
    data = a3_direction_binned_lat_lon_evidence_ratio,
    mapping = aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
    size = 3,
    colour = "black",
    stroke = 0.1
  ) +
  theme_bw() +
  scale_fill_manual(
    values = custom_palette(),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(21, 22, 23, 25, 24),
    drop = FALSE
  ) +
  # expand map
  coord_sf(xlim = c(95, 180), ylim = c(-60, 0)) +
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
    from = c(145, 155, -30, -16),
    to = c(157, 178, -29.5, 1.5),
    shadow = FALSE,
    expand = 0,
    plot = inset_plot_QLD,
    proj = "single"
  ) +
  # magnify NSW
  geom_magnify(
    from = c(146.5, 154, -38, -28),
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
    fill = "Evidence Ratio",
    shape = bquote("Impact of "~CO[2]~"Term")
  ) +
  theme(
    legend.key = element_rect(fill = "grey80"),
    legend.title = element_text(hjust = 0.5),
    legend.background = element_rect(colour = "black"),
    axis.text = element_text(size = 6), 
    legend.position = "inside",
    legend.position.inside = c(0.351, 0.9),
    legend.box = "horizontal"#, # side-by-side legends
    #plot.margin = margin(20, 1, 2, 2, unit = "mm") # white area around figure
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 5, shape = 21), nrow = 3), # Wrap legend with nrow
    shape = guide_legend(override.aes = list(size = 5, fill = "grey50"), nrow = 3)
  )

#single_map_aus

ggsave(
  filename = "./Graphs/CMAES_graphs/evidence_ratio_aus_with_zoom_v2.pdf",
  plot = single_map_aus,
  device = "pdf",
  width = 237,
  height = 210,
  units = "mm"
)

stop_here <- tactical_typo()

# 3. Streamflow plots of best-CO2, best-non-CO2 and observed -------------------


## Filter based on best non-CO2 and CO2 models =================================
best_streamflow_results <- streamflow_results |>
  semi_join(
    best_CO2_non_CO2_per_gauge,
    by = join_by(gauge, streamflow_model)
  )


## Only include best streamflow that was calibrated on =========================
## join the included_in_calibration column
in_calibration <- data |> 
  select(year, gauge, included_in_calibration)

best_calibration_streamflow_results <- best_streamflow_results |> 
  left_join(
    in_calibration,
    by = join_by(year, gauge)
  ) |> 
  filter(included_in_calibration)



## Summarise results into a tidy format ========================================
tidy_boxcox_streamflow <- best_calibration_streamflow_results |>
  drop_na() |>  # only include if observed streamflow is present
  pivot_longer(
    cols = c(observed_boxcox_streamflow, modelled_boxcox_streamflow),
    names_to = "name",
    values_to = "boxcox_streamflow"
  ) |>
  mutate(
    name = if_else(name == "observed_boxcox_streamflow", "observed", "modelled")
  ) |> 
  mutate(
    contains_CO2 = str_detect(streamflow_model, "CO2")
  ) |> 
  mutate(
    legend = case_when(
      contains_CO2 & (name != "observed") ~ "modelled_CO2",
      !contains_CO2 & (name != "observed") ~ "modelled_non_CO2",
      .default = name
    )
  ) |> 
  select(
    !c(streamflow_model, objective_function, included_in_calibration, name, contains_CO2)
    )


## Convert from box-cox space to real space ====================================
### bc lambda found in gauge_information
# boxcox_inverse_transform()
tidy_streamflow <- tidy_boxcox_streamflow |>
  left_join(
    gauge_information,
    by = join_by(gauge)
  ) |>
  select(!c(state, lat, lon)) |>
  mutate(
    streamflow = boxcox_inverse_transform(yt = boxcox_streamflow, lambda = bc_lambda, lambda_2 = 1),
    .by = gauge
  )






## Plot results ================================================================
### having 534 graphs on a single page is too much = slows pc.
### Spread across multiple

# rep by gauge

# Split the tibble into X groups
# Gauges must not be across multiple groups
# Randomly assign 1,2 or 3 to each group then split? This works
# but is probably not the best way of doing it
ready_split_tidy_streamflow <- tidy_streamflow |> 
  mutate(
    split = sample(c(1, 2, 3), size = 1, replace = TRUE),
    .by = gauge,
    .before = 1
  ) 

split_tidy_streamflow <- ready_split_tidy_streamflow |> 
  split(f = ready_split_tidy_streamflow$split)


## Plot rainfall-runoff relationships ==========================================
## Need to split plot into smaller chunks. PC struggles to load it.

repeat_rainfall_runoff_plots <- function(segmented_tidy_streamflow) {
  rainfall_runoff_plots <- segmented_tidy_streamflow |> 
    ggplot(
      aes(
        x = precipitation, 
        y = boxcox_streamflow, 
        colour = legend, 
        shape = legend
        )
      ) +
    geom_point() +
    geom_line(
      stat = "smooth",
      method = "lm",
      formula = y ~ x,
      alpha = 0.4,
      linewidth = 1
    ) +
    labs(
      x = "Annual Precipitation (mm)",
      y = "Annual Boxcox Streamflow (mm)"
    ) +
    scale_colour_brewer(palette = "Set1") +
    theme_bw() +
    facet_wrap(~gauge, scales = "free") +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom"
    )
}


rainfall_runoff_graphs <- map(
  .x = split_tidy_streamflow,
  .f = repeat_rainfall_runoff_plots
)

ggsave(
  filename = paste0("rainfall_runoff_comparison_", get_date(), ".pdf"),
  plot = gridExtra::marrangeGrob(rainfall_runoff_graphs, nrow = 1, ncol = 1),
  device = "pdf",
  path = "./Graphs/CMAES_graphs",
  width = 1189,
  height = 841,
  units = "mm"
)


## Plot streamflow time ========================================================
chunk_streamflow_timeseries_plot <- function(data) {
  
  data |>
    ggplot(aes(x = year, y = streamflow, colour = legend)) +
    geom_line(na.rm = TRUE, alpha = 0.5) +
    geom_point(na.rm = TRUE, size = 0.5, alpha = 0.5) +
    theme_bw() +
    scale_colour_brewer(palette = "Set1") +
    labs(
      x = "Year",
      y = "Streamflow (mm)"
    ) +
    facet_wrap(~gauge, scales = "free_y") +
    theme(legend.title = element_blank())
  
}


plot_streamflow_timeseries <- map(
  .x = split_tidy_streamflow,
  .f = chunk_streamflow_timeseries_plot
)



ggsave(
  filename = paste0("streamflow_timeseries_comparison_", get_date(), ".pdf"),
  plot = gridExtra::marrangeGrob(plot_streamflow_timeseries, nrow = 1, ncol = 1),
  device = "pdf",
  path = "./Graphs/CMAES_graphs",
  width = 1189,
  height = 841,
  units = "mm"
)

