# Evidence ratio analysis

# Figures produced in this R file ----------------------------------------------

# 1. Main --> evidence_ratio_aus_map.pdf
# 2. Supplementary --> evidence_ratio_vs_catchment_area.pdf
# 3. Supplementary --> evidence_ratio_vs_record_length.pdf
# 4. Supplementary --> evidence_ratio_vs_prop_forested.pdf
# 5. Supplementary --> sens_slope_evidence_ratio.pdf
# 6. Other --> sens_slope_map.pdf


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
  col_select = c(
    "gauge",
    "lat",
    "lon",
    "state",
    "status",
    "parent_catchment",
    "catchment_area_sq_km",
    "record_length",
    "prop_forested",
    "river_di"
  ),
  show_col_types = FALSE
)




best_CO2_non_CO2_per_gauge <- read_csv(
  "./Modelling/Results/CMAES/best_CO2_non_CO2_per_catchment_CMAES.csv",
  show_col_types = FALSE
)


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
    gauge_information,
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
  ) #|>
  #filter(river_di < 0.05)


## save evidence ratio for other analysis
write_csv(
  a3_direction_binned_lat_lon_evidence_ratio,
  "./Modelling/Results/CMAES/evidence_ratio_results.csv",
)



# Make final plot --------------------------------------------------------------

aus_map <- generate_aus_map_sf()


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
    drop = FALSE # Very important - every plots has the same discrete values
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
  filename = "./Figures/Main/evidence_ratio_aus_map.pdf",
  plot = single_map_aus,
  device = "pdf",
  width = 232,
  height = 200, # 210,
  units = "mm"
)



# How similar are nested catchments --------------------------------------------
## Do nested catchments share the same binned ratio (e.g., Weak) ===============

a3_direction_binned_lat_lon_evidence_ratio |> 
  filter(status != "no_nested") |> 
  ggplot(aes(x = parent_catchment, fill = binned_evidence_ratio)) +
  geom_bar(position = "stack") +
  theme_bw() +
  labs(
    x = "Parent Catchment",
    y = "Catchments in Parent Catchment",
    fill = "Evidence Ratio"
  ) +
  theme(
    axis.text.x = element_text(angle = 90),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(breaks = seq(from = 0, to = 14))

# based on this graph nested catchments - the vast majority share the same
# binned ratio
# 26/85 nested catchments have different evidence ratio
# 70 % of all nested catchments share the same qualitative evidence ratio
# Of the 26 which don't:
# - 9 only have at least moderate or greater
# - 17 have a mix - tend to be more of one type >> weak
# - 50/50 weak and not weak = 5


# Do nested catchments share the same type of change if not weak ===============
# (e.g., negative slope) or at least same direction

# filter out catchments with only weak evi
multiple_parent_catchment_evi_ratio <- a3_direction_binned_lat_lon_evidence_ratio |> 
  filter(status != "no_nested") |> 
  mutate(
    binned_evidence_ratio = as.character(binned_evidence_ratio)
  ) |> 
  summarise(
    unique_evi_ratio = length(unique(binned_evidence_ratio)),
    .by = parent_catchment
  ) |> 
  filter(unique_evi_ratio != 1) |> 
  pull(parent_catchment)
  

a3_direction_binned_lat_lon_evidence_ratio |> 
  filter(parent_catchment %in% multiple_parent_catchment_evi_ratio) |> 
  ggplot(aes(x = parent_catchment, fill = impact_of_CO2_term)) +
  geom_bar(position = "stack") +
  theme_bw() +
  labs(
    x = "Parent Catchment",
    y = "Catchments in Parent Catchment",
    fill = "Impact of CO2 Term"
  ) +
  theme(
    axis.text.x = element_text(angle = 90),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(breaks = seq(from = 0, to = 14))

# From graph (24 parent catchments)
# - of parent catchments with multiple evi ratio values
# - 3 have altering signs (positive and negative change) - pretty much all
#   are in the same direction


# Compare evidence ratio between wild and non-wild catchments ------------------

## Make nice breaks for all evidence ratio plots
## use y_axis_scale_transform variable 
evi_ratio_range <- evidence_ratio_calc |> pull(evidence_ratio) |> range()
y_axis_evi_min <- round_any(x[1], accuracy = 10, f = floor)
y_axis_evi_max <- round_any(x[2], accuracy = 1E16, f = ceiling)
nice_breaks_log10 <- seq(from = log10(abs(y_axis_evi_min)), to = log10(y_axis_evi_max), by = 3)
nice_breaks_log10[1] <- -1
nice_breaks <- 10^nice_breaks_log10
y_axis_scale_transform <- scale_y_continuous(
  transform = scales::pseudo_log_trans(base = 10),
  breaks = nice_breaks
)


stat_comparison_wild_non_wild <- a3_direction_binned_lat_lon_evidence_ratio |> 
  mutate(
    wild = river_di < 0.05 # 0.05 from CAMELS-AUSv2
  )

wild_vs_non_wild_boxplot <- stat_comparison_wild_non_wild |> 
  ggplot(aes(x = wild, y = evidence_ratio)) +
  geom_boxplot(
    staplewidth = 0.5,
    outlier.alpha = 0.7,
    outlier.colour = "black",
    outlier.fill = "grey",
    outlier.shape = 21,
    whisker.linewidth = 0.2,
    box.linewidth = 0.2,
    median.linewidth = 0.4,
    staple.linewidth = 0.2,
    fill = "grey90"
    ) +
  geom_text(
    data = single_label(x_pos = FALSE, y_pos = 5E15, label_name = "e"),
    mapping = aes(x = x_pos, y = y_pos, label = label_name),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 12,
    size.unit = "pt",
    hjust = 7
  ) +
  y_axis_scale_transform +
  labs(
    x = "Is \"Near Natural\" Catchment",# (river_di < 0.05)",
    y = "Evidence Ratio"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8)
  )





wilcox.test(
  evidence_ratio ~ wild, 
  data = stat_comparison_wild_non_wild,
  alternative = "two.sided"
  ) 
# Retain the null hypothesis (no difference in evidence ratio between the two groups)




# Relationship between evidence ratio and catchment area -----------------------
single_label <- function(x_pos, y_pos, label_name) { # for adding a, b, c labels
  tribble(
    ~x_pos, ~y_pos, ~label_name,
    x_pos,  y_pos,  label_name
  )
}


## Get catchment area and record length from gauge data
gauge_area_and_record_length <- gauge_information |>
  select(gauge, catchment_area_sq_km, record_length, prop_forested)


evidence_ratio_vs_catchment_area <- a3_direction_binned_lat_lon_evidence_ratio |>
  ggplot(aes(x = catchment_area_sq_km, y = evidence_ratio)) +
  geom_point(
    fill = "grey",
    colour = "black",
    stroke = 0.1,
    shape = 21,
    size = 2,
    alpha = 0.8
  ) +
  geom_text(
    data = single_label(x_pos = 5, y_pos = 5E15, label_name = "a"),
    mapping = aes(x = x_pos, y = y_pos, label = label_name),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 12,
    size.unit = "pt"
  ) +
  y_axis_scale_transform +
  scale_x_log10() +
  labs(
    x = bquote("Catchment Area ("*km^2*")"), # bquote("X Axis Label ("*m^2*")")
    y = "Evidence Ratio"
    ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8)
  )


ggsave(
  filename = "evidence_ratio_vs_catchment_area.pdf",
  plot = evidence_ratio_vs_catchment_area,
  path = "Figures/Supplementary",
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)


# Relationship between evidence ratio and record length ------------------------
evidence_ratio_vs_record_length <- a3_direction_binned_lat_lon_evidence_ratio |>
  ggplot(aes(x = record_length, evidence_ratio)) +
  geom_jitter( # stop dots overlapping
    fill = "grey",
    colour = "black",
    stroke = 0.1,
    shape = 21,
    size = 2,
    alpha = 0.8
  ) +
  geom_text(
    data = single_label(x_pos = 30, y_pos = 5E15, label_name = "b"),
    mapping = aes(x = x_pos, y = y_pos, label = label_name),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 12,
    size.unit = "pt"
  ) +
  y_axis_scale_transform +
  labs(x = "Record Length (Years)", y = "Evidence Ratio") +
  theme_bw() +
  theme(
    axis.title.y = element_blank(), # remove double up of y axis when combining
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8)
  )


ggsave(
  filename = "evidence_ratio_vs_record_length.pdf",
  plot = evidence_ratio_vs_record_length,
  path = "Figures/Supplementary",
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)


# Relationship between evidence ratio and forested catchment -------------------
evidence_ratio_vs_prop_forested <- a3_direction_binned_lat_lon_evidence_ratio |>
  ggplot(aes(x = prop_forested, evidence_ratio)) +
  geom_point(
    fill = "grey",
    colour = "black",
    stroke = 0.1,
    shape = 21,
    size = 2,
    alpha = 0.8
  ) +
  geom_text(
    data = single_label(x_pos = 0.03, y_pos = 5E15, label_name = "c"),
    mapping = aes(x = x_pos, y = y_pos, label = label_name),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 12,
    size.unit = "pt"
  ) +
  y_axis_scale_transform +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Proportion of Catchment Forested", y = "Evidence Ratio") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8)
  )


ggsave(
  filename = "evidence_ratio_vs_prop_forested.pdf",
  plot = evidence_ratio_vs_prop_forested,
  path = "Figures/Supplementary",
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)


# Relationship between evidence ratio and sens slope of annual rainfall --------

## Calculate sen's slope of rainfall ===========================================
simple_sens_slope <- function(x) {
  result <- x |>
    as.ts() |>
    sens.slope()

  return(unname(result$estimates))
}

get_slope <- function(x, y, ...) {
  lm(y ~ x, ...)$coefficients[2] |> unname() # position of slope
}

rainfall_sens_slope_data <- data |>
  summarise(
    sens_slope = simple_sens_slope(p_mm),
    mean_annual_rainfall = mean(p_mm),
    lin_reg_slope = get_slope(x = year, y = p_mm),
    .by = gauge
  ) |>
  mutate(
    standard_sens_slope = sens_slope / mean_annual_rainfall
  )

# Test decade changes
# rainfall_sens_slope_data <- data |>
#  mutate(decade = year - year %% 10 ) |>
#  summarise(
#    mean_decade_rainfall = mean(p_mm),
#    .by = c(decade, gauge)
#  ) |>
#  summarise(
#    sens_slope = simple_sens_slope(mean_decade_rainfall),
#    mean_decade_rainfall = mean(mean_decade_rainfall),
#    .by = gauge
#  )


all_evidence_ratio_information <- additional_info_a3_direction_binned_evidence_ratio |>
  left_join(
    rainfall_sens_slope_data,
    by = join_by(gauge)
  ) |>
  arrange(standard_sens_slope)


# test decade changes
# all_evidence_ratio_information |>
#  ggplot(aes(x = sens_slope, y = evidence_ratio)) +
#  geom_point() +
#  theme_bw() +
#  scale_y_log10()

# Sen's slope/trend of rainfall is standardised
# You get the same result if you compare against the anomaly


## Plot  =======================================================================

### Custom scale for annual rainfall ###########################################
mean_annual_rainfall_limits <- all_evidence_ratio_information |>
  pull(mean_annual_rainfall) |>
  range()

mean_annual_rainfall_limits <- c(floor(mean_annual_rainfall_limits[1]), ceiling(mean_annual_rainfall_limits[2]))

mean_annual_rainfall_breaks <- all_evidence_ratio_information |>
  pull(mean_annual_rainfall) |>
  quantile(probs = seq(0, 1, length.out = 8L)) |>
  unname() |>
  round_any(accuracy = 100, f = round)

mean_annual_rainfall_breaks <- seq(from = 600, to = 1100, by = 100)

# mean_annual_rainfall_breaks <- mean_annual_rainfall_breaks[-c(1, length(mean_annual_rainfall_breaks))]

mean_annual_rainfall_palette <- function(x) {
  c("#f7fcfd", "#e0ecf4", "#bfd3e6", "#9ebcda", "#8c96c6", "#8c6bb1", "#88419d", "#810f7c", "#4d004b")
}





sens_slope_evidence_ratio_plot <- all_evidence_ratio_information |>
  ggplot(aes(x = standard_sens_slope, y = evidence_ratio, fill = mean_annual_rainfall)) +
  geom_point(
    colour = "black",
    stroke = 0.1,
    shape = 21,
    size = 2,
    alpha = 0.8
  ) +
  geom_text(
    data = single_label(x_pos = -0.004375, y_pos = 5E15, label_name = "d"),
    mapping = aes(x = x_pos, y = y_pos, label = label_name),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 12,
    size.unit = "pt"
  ) +
  y_axis_scale_transform +
  binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
    aesthetics = "fill",
    palette = mean_annual_rainfall_palette,
    breaks = mean_annual_rainfall_breaks,
    limits = mean_annual_rainfall_limits,
    show.limits = TRUE,
    guide = "colorsteps"
  ) +
  labs(
    x = "Normalised Annual Rainfall Sen's Slope",
    y = "Evidence Ratio",
    fill = "Mean Annual Rainfall (mm)"
  ) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.685, 0.87),
    legend.background = element_rect(colour = "black", fill = NULL, linewidth = 0.1),
    legend.title = element_text(hjust = 0.5, size = 7),
    legend.text = element_text(size = 5),
    legend.margin = margin(t = 1, r = 8, b = 1, l = 8, unit = "pt"),
    axis.title.y = element_blank(), # remove double up of y axis when combining
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 9),
    axis.text = element_text(size = 8)
  ) +
  guides(
    fill = guide_coloursteps(
      barwidth = unit(4, "cm"),
      barheight = unit(0.2, "cm"),
      show.limits = TRUE,
      even.steps = TRUE,
      title.position = "top",
      direction = "horizontal"
    )
  )


ggsave(
  filename = "sens_slope_evidence_ratio_plot.pdf",
  plot = sens_slope_evidence_ratio_plot,
  path = "Figures/Supplementary",
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)



# Add rainfall seasonality vs. evi ratio graph ---------------------------------
# original code in rainfall_seasonality_analysis.R
# this is just to add it to the plot
seasonal_rainfall_classification_gauge_info <- read_csv(
  "Modelling/seasonal_rainfall_data.csv",
  show_col_types = FALSE
)

evidence_ratio_filter <- -100

plotting_seasonal_rainfall_classification_gauge_info <- seasonal_rainfall_classification_gauge_info |>
  filter(evidence_ratio > evidence_ratio_filter)


count_occurences <- plotting_seasonal_rainfall_classification_gauge_info |>
  count(classification) |>
  mutate(
    y_pos = min(plotting_seasonal_rainfall_classification_gauge_info$evidence_ratio)
  ) |>
  mutate(
    n_label = paste0("n = ", n)
  )




seasonal_evidence_ratio_boxplot <- plotting_seasonal_rainfall_classification_gauge_info |>
  ggplot(aes(x = classification, y = evidence_ratio)) +
  geom_boxplot(
    # exact same as wild boxplot
    staplewidth = 0.5,
    outlier.alpha = 0.7,
    outlier.colour = "black",
    outlier.fill = "grey",
    outlier.shape = 21,
    whisker.linewidth = 0.2,
    box.linewidth = 0.2,
    median.linewidth = 0.4,
    staple.linewidth = 0.2,
    fill = "grey90"
  ) +
  #geom_text(
  #  aes(x = classification, y = y_pos, label = n_label),
  #  data = count_occurences,
  #  vjust = 1.5
  #) +
  geom_text(
    data = single_label(x_pos = "Arid", y_pos = 5E15, label_name = "f"),
    mapping = aes(x = x_pos, y = y_pos, label = label_name),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 12,
    size.unit = "pt",
    hjust = 3
  ) +
  # this must be the same as the wild catchment
  y_axis_scale_transform +
  scale_x_discrete(drop = FALSE, guide = guide_axis(n.dodge = 2)) +
  labs(
    x = "Seasonal Rainfall Zones",
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  ) 

# Combine into a 3 x 2 figure --------------------------------------------------
evidence_ratio_extended_data_3x2 <- (evidence_ratio_vs_catchment_area |
  evidence_ratio_vs_record_length) /
  (evidence_ratio_vs_prop_forested |
  sens_slope_evidence_ratio_plot) /
  (wild_vs_non_wild_boxplot |
  seasonal_evidence_ratio_boxplot)


ggsave(
  filename = "evidence_ratio_extended_data_3x2.pdf",
  plot = evidence_ratio_extended_data_3x2,
  path = "Figures/Extended_Data",
  device = "pdf",
  width = 180,
  height = 190,
  units = "mm"
)





# My own sen slope function - same as function from package
# my_sens_slope <- function(x, t) {
# the length of x and t must be the same
# stopifnot(length(x) == length(t))

# t must be continuous
# t_lag_1 <- lag(t, n = 1L)
# diff_t <- t - t_lag_1
# stopifnot(any(diff_t[-1] == 1))

# pre-allocate array
# length_pre_allocation <- sum(1:(length(x) - 1))
# d <- numeric(length = length_pre_allocation)
# seperate_index <- 1

# for (j in 2:length(x)) {
#  for (i in j:length(x)) {
#    d[seperate_index] <- (x[i] - x[j - 1]) / (t[i] - t[j - 1])
#    seperate_index <- seperate_index + 1
#  }
# }

# return(d)
# return(median(d, na.rm = TRUE))
# }

# x <- data |>
# filter(gauge == "314213")

# check sens slope values
# my_sens_slope(x = x$p_mm, t = x$year)



