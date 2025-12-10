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
pacman::p_load(tidyverse, ozmaps, sf, ggmagnify, trend)



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
  select(gauge, lat, lon)


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






# Relationship between evidence ratio and catchment area -----------------------
## Get catchment area and record length from gauge data
gauge_area_and_record_length <- gauge_information |>
  select(gauge, catchment_area_sq_km, record_length, prop_forested)


## Add gauge information to a3_direction_binned_lat_lon_evidence_ratio
additional_info_a3_direction_binned_evidence_ratio <- a3_direction_binned_lat_lon_evidence_ratio |>
  left_join(
    gauge_area_and_record_length,
    by = join_by(gauge)
  )


evidence_ratio_vs_catchment_area <- additional_info_a3_direction_binned_evidence_ratio |>
  filter(evidence_ratio > 0) |>
  ggplot(aes(x = catchment_area_sq_km, evidence_ratio)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "Catchment Area (km2)", y = "Evidence Ratio") +
  theme_bw()


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
evidence_ratio_vs_record_length <- additional_info_a3_direction_binned_evidence_ratio |>
  filter(evidence_ratio > 0) |>
  ggplot(aes(x = record_length, evidence_ratio)) +
  geom_jitter() + # stop dots overlapping
  scale_y_log10() +
  labs(x = "Record Length (Years)", y = "Evidence Ratio") +
  theme_bw()


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
evidence_ratio_vs_prop_forested <- additional_info_a3_direction_binned_evidence_ratio |>
  filter(evidence_ratio > 0) |>
  ggplot(aes(x = prop_forested, evidence_ratio)) +
  geom_point() +
  scale_y_log10() +
  labs(x = "Proportion of forested", y = "Evidence Ratio") +
  theme_bw()


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
#rainfall_sens_slope_data <- data |>
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
  filter(evidence_ratio > 0) |>
  arrange(standard_sens_slope) 
  

# test decade changes
#all_evidence_ratio_information |> 
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

#mean_annual_rainfall_breaks <- mean_annual_rainfall_breaks[-c(1, length(mean_annual_rainfall_breaks))]

mean_annual_rainfall_palette <- function(x) {
  c("#f7fcfd", "#e0ecf4", "#bfd3e6", "#9ebcda", "#8c96c6", "#8c6bb1", "#88419d", "#810f7c", "#4d004b")
}


sens_slope_evidence_ratio_plot <- all_evidence_ratio_information |>
  ggplot(aes(x = standard_sens_slope, y = evidence_ratio, fill = mean_annual_rainfall)) +
  geom_point(
    colour = "black",
    stroke = 0.1,
    shape = 21,
    size = 3
  ) +
  scale_y_log10() +
  binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
    aesthetics = "fill",
    palette = mean_annual_rainfall_palette,
    breaks = mean_annual_rainfall_breaks,
    limits = mean_annual_rainfall_limits,
    show.limits = TRUE,
    guide = "colorsteps"
  ) +
  labs(
    x = "Normalized Annual Rainfall Sen's Slope",
    y = "Evidence Ratio",
    fill = "Mean Annual Rainfall (mm)"
  ) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.75, 0.9),
    legend.background = element_rect(colour = "black", fill = NULL),
    legend.title = element_text(hjust = 0.5),
    legend.margin = margin(t = 10, r = 20, b = 10, l = 20, unit = "pt"),
  ) +
  guides(
    fill = guide_coloursteps(
      barwidth = unit(10, "cm"),
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


# This can be deleted ...
## Map =========================================================================
QLD_data <- all_evidence_ratio_information |>
  filter(state == "QLD")

NSW_data <- all_evidence_ratio_information |>
  filter(state == "NSW")

VIC_data <- all_evidence_ratio_information |>
  filter(state == "VIC")

WA_data <- all_evidence_ratio_information |>
  filter(state == "WA")

TAS_data <- all_evidence_ratio_information |>
  filter(state == "TAS")


### All colour scales must be the same #########################################
sens_range <- all_evidence_ratio_information |>
  pull(standard_sens_slope) |>
  range()

# round by itself does not do a good job - a single variable outside of range
sens_range <- c(
  round_any(sens_range[1], accuracy = 0.001, f = floor),
  round_any(sens_range[2], accuracy = 0.001, f = ceiling)
)



### Generate inset plots #######################################################

inset_plot_QLD <- aus_map |>
  filter(state == "QLD") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = QLD_data,
    aes(x = lon, y = lat, fill = standard_sens_slope),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1,
    colour = "black",
    shape = 21
  ) +
  scale_fill_distiller(
    palette = "RdBu",
    limits = sens_range,
    direction = 1
  ) +
  theme_void()


inset_plot_NSW <- aus_map |>
  filter(state == "NSW") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = NSW_data,
    aes(x = lon, y = lat, fill = standard_sens_slope),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1,
    colour = "black",
    shape = 21
  ) +
  scale_fill_distiller(
    palette = "RdBu",
    limits = sens_range,
    direction = 1
  ) +
  theme_void()


inset_plot_VIC <- aus_map |>
  filter(state == "VIC") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = VIC_data,
    aes(x = lon, y = lat, fill = standard_sens_slope),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1,
    colour = "black",
    shape = 21
  ) +
  scale_fill_distiller(
    palette = "RdBu",
    limits = sens_range,
    direction = 1
  ) +
  theme_void()



inset_plot_WA <- aus_map |>
  filter(state == "WA") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = WA_data,
    aes(x = lon, y = lat, fill = standard_sens_slope),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1,
    colour = "black",
    shape = 21
  ) +
  scale_fill_distiller(
    palette = "RdBu",
    limits = sens_range,
    direction = 1
  ) +
  theme_void()



inset_plot_TAS <- aus_map |>
  filter(state == "TAS") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = TAS_data,
    aes(x = lon, y = lat, fill = standard_sens_slope),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1,
    colour = "black",
    shape = 21
  ) +
  scale_fill_distiller(
    palette = "RdBu",
    limits = sens_range,
    direction = 1
  ) +
  theme_void()



## Put it together =============================================================

sens_slope_map <- aus_map |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = all_evidence_ratio_information,
    aes(x = lon, y = lat, fill = standard_sens_slope),
    size = 2.5,
    stroke = 0.1,
    colour = "black",
    shape = 21
  ) +
  scale_fill_distiller(
    palette = "RdBu",
    limits = sens_range,
    direction = 1
  ) +
  theme_bw() +
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
    fill = "Normalized Annual Rainfall Sen's Slope"
  ) +
  theme(
    legend.title = element_text(hjust = 0.5),
    legend.title.position = "top",
    legend.background = element_rect(colour = "black"),
    axis.text = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.346, 0.9), # constants used to move the legend in the right place
    legend.box = "horizontal", # side-by-side legends
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    legend.margin = margin(t = 5, b = 5, r = 20, l = 20, unit = "pt") # add extra padding around legend box to avoid -1.6 intersecting with line
  ) +
  guides(
    fill = guide_colourbar(
      direction = "horizontal",
      barwidth = unit(10, "cm"),
      show.limits = TRUE
    )
  )
  


ggsave(
  filename = "./Figures/Other/sens_slope_map.pdf",
  plot = sens_slope_map,
  device = "pdf",
  width = 232,
  height = 200, # 210,
  units = "mm"
)
