# Examine PET


# Import libraries required ----------------------------------------------------
pacman::p_load(tidyverse, ozmaps, sf, ggmagnify)


# Import functions -------------------------------------------------------------
source("./Functions/utility.R")

# Import data ------------------------------------------------------------------
# CAMELS-AUSv2 paper suggests using Morton Wet for APET
areal_potential_evap_SILO_daily <- read_csv(
  "Data/Raw/et_morton_wet_SILO.csv",
  show_col_types = FALSE
)

data <- readr::read_csv( # data will be in the package
  "./Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  show_col_types = FALSE
) |>
  mutate(
    year = as.integer(year)
  ) |>
  # required for log-sinh. Log-sinh current formulation has asymptote of zero.
  # This means zero flows of ephemeral catchments cannot be transformed
  # add a really small value
  mutate(q_mm = q_mm + .Machine$double.eps^0.5)


gauge_information <- readr::read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE
)

lat_lon_gauge <- gauge_information |>
  select(gauge, state, lat, lon)


best_CO2_non_CO2_per_gauge <- read_csv(
  "./Modelling/Results/CMAES/best_CO2_non_CO2_per_catchment_CMAES.csv",
  show_col_types = FALSE
)

trends_in_counterfactual_flow <- read_csv(
  "./Modelling/Other/trends_in_counterfactual_flow.csv",
  show_col_types = FALSE
  )

# Tidy data --------------------------------------------------------------------
relevant_years <- data |>
  pull(year) |>
  unique()
relevant_gauges <- data |>
  pull(gauge) |>
  unique()

evap_areal_potential_annual <- areal_potential_evap_SILO_daily |>
  pivot_longer(
    cols = !c(year, month, day),
    names_to = "gauge",
    values_to = "APET_mm"
  ) |>
  # only include years we are interested in
  filter(year %in% relevant_years) |>
  # only include gauges we are interested in
  filter(gauge %in% relevant_gauges) |>
  # sum daily PET to get annual - check for missing data
  summarise(
    annual_APET_mm = sum(APET_mm),
    n = n(),
    .by = c(year, gauge)
  ) |>
  # filter out incomplete years
  filter(n %in% c(365, 366))



# Plot PET data over time ------------------------------------------------------
PET_timeseries <- evap_areal_potential_annual |>
  ggplot(aes(x = year, y = annual_APET_mm)) +
  geom_line() +
  geom_smooth(
    method = lm,
    formula = y ~ x,
    se = FALSE,
    colour = "red",
    linetype = "dashed"
  ) +
  theme_bw() +
  labs(
    x = "Year",
    y = "PET (mm)"
  ) +
  facet_wrap(~gauge, scales = "free_y")

ggsave(
  filename = "./Figures/Other/silo_PET_timeseries.pdf",
  plot = PET_timeseries,
  device = "pdf",
  width = 1189,
  height = 841,
  units = "mm"
)



# Find the PET trend in data using linear slope --------------------------------
get_slope <- function(x, y, ...) {
  lm(y ~ x, ...)$coefficients[2] |> unname() # position of slope
}

timeseries_slopes_evap_areal_potential <- evap_areal_potential_annual |>
  #filter(year %in% c(seq(1990, 1999))) |> 
  #filter(year %in% c(seq(2012, 2021))) |> 
  right_join(
    trends_in_counterfactual_flow,
    by = join_by(gauge, year)
  ) |> 
  mutate(
    flow_partitioning_diff = realspace_a3_off_streamflow - realspace_a3_on_streamflow
  ) |> 
  drop_na() |> 
  summarise(
    trend_mm_per_y = get_slope(x = year, y = annual_APET_mm),
    trend_mm_per_y_flow = get_slope(x = year, y = flow_partitioning_diff),
    .by = gauge
  ) |>
  left_join(
    lat_lon_gauge,
    by = join_by(gauge)
  ) 

# is there no APET data for 405263? --> look into this
timeseries_slopes_evap_areal_potential |> 
  summarise(
    mean_APET_slope = mean(trend_mm_per_y, na.rm = T), 
    mean_Q_slope = mean(trend_mm_per_y_flow, na.rm = T)
    )


# Plot PET trends on a map -----------------------------------------------------
aus_map <- generate_aus_map_sf()


## Generate Insets =============================================================
### Filter data by state #######################################################

QLD_data <- timeseries_slopes_evap_areal_potential |>
  filter(state == "QLD")

NSW_data <- timeseries_slopes_evap_areal_potential |>
  filter(state == "NSW")

VIC_data <- timeseries_slopes_evap_areal_potential |>
  filter(state == "VIC")

WA_data <- timeseries_slopes_evap_areal_potential |>
  filter(state == "WA")

TAS_data <- timeseries_slopes_evap_areal_potential |>
  filter(state == "TAS")


### All colour scales must be the same #########################################
trend_range <- timeseries_slopes_evap_areal_potential |>
  pull(trend_mm_per_y) |>
  range()

# round by itself does not do a good job - a single variable outside of range
trend_range <- c(
  round_any(trend_range[1], accuracy = 0.01, f = floor),
  round_any(trend_range[2], accuracy = 0.01, f = ceiling)
)

trend_breaks <- c(trend_range[1], -0.8, 0, 0.8, trend_range[2]) # need to to show nice breaks

### Generate inset plots #######################################################

inset_plot_QLD <- aus_map |>
  filter(state == "QLD") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = QLD_data,
    aes(x = lon, y = lat, fill = trend_mm_per_y),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1,
    colour = "black",
    shape = 21
  ) +
  scale_fill_distiller(
    palette = "RdBu",
    limits = trend_range,
    breaks = trend_breaks
  ) +
  theme_void()


inset_plot_NSW <- aus_map |>
  filter(state == "NSW") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = NSW_data,
    aes(x = lon, y = lat, fill = trend_mm_per_y),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1,
    colour = "black",
    shape = 21
  ) +
  scale_fill_distiller(
    palette = "RdBu",
    limits = trend_range,
    breaks = trend_breaks
  ) +
  theme_void()


inset_plot_VIC <- aus_map |>
  filter(state == "VIC") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = VIC_data,
    aes(x = lon, y = lat, fill = trend_mm_per_y),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1,
    colour = "black",
    shape = 21
  ) +
  scale_fill_distiller(
    palette = "RdBu",
    limits = trend_range,
    breaks = trend_breaks
  ) +
  theme_void()



inset_plot_WA <- aus_map |>
  filter(state == "WA") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = WA_data,
    aes(x = lon, y = lat, fill = trend_mm_per_y),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1,
    colour = "black",
    shape = 21
  ) +
  scale_fill_distiller(
    palette = "RdBu",
    limits = trend_range,
    breaks = trend_breaks
  ) +
  theme_void()



inset_plot_TAS <- aus_map |>
  filter(state == "TAS") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = TAS_data,
    aes(x = lon, y = lat, fill = trend_mm_per_y),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1,
    colour = "black",
    shape = 21
  ) +
  scale_fill_distiller(
    palette = "RdBu",
    limits = trend_range,
    breaks = trend_breaks
  ) +
  theme_void()



## Put it together =============================================================

single_map_aus <- aus_map |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = timeseries_slopes_evap_areal_potential,
    aes(x = lon, y = lat, fill = trend_mm_per_y),
    size = 2.5,
    stroke = 0.1,
    colour = "black",
    shape = 21
  ) +
  scale_fill_distiller(
    palette = "RdBu",
    limits = trend_range,
    breaks = trend_breaks
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
    fill = "Trend in PET (mm/y)"
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
      barwidth = unit(8, "cm")
    )
  )


single_map_aus

ggsave(
  filename = "./Figures/Supplementary/PET_trends_aus_map.pdf",
  plot = single_map_aus,
  device = "pdf",
  width = 232,
  height = 200, # 210,
  units = "mm"
)



# Correlations -----------------------------------------------------------------
correlations_data <- data |>
  left_join(
    evap_areal_potential_annual,
    by = join_by(gauge, year)
  ) |> # join t_max if required
  # removing missing values
  drop_na()

## compare APET vs. Precip =====================================================
correlation_APET_vs_P <- correlations_data |>
  summarise(
    corr_P_vs_APET = cor(p_mm, annual_APET_mm, use = "complete.obs"),
    xlab = max(p_mm) * 0.95,
    ylab = max(annual_APET_mm) * 0.99,
    .by = gauge
  ) |>
  mutate(
    R2_P_vs_APET = corr_P_vs_APET^2,
    R2_label = round(R2_P_vs_APET, digits = 2)
  ) |> 
  left_join(
    lat_lon_gauge,
    by = join_by(gauge)
  )

APET_vs_P_plot <- correlations_data |>
  ggplot(aes(x = p_mm, y = annual_APET_mm)) +
  geom_smooth(
    method = lm,
    formula = y ~ x,
    se = FALSE,
    colour = "red"
  ) +
  geom_point() +
  geom_label(
    aes(x = xlab, y = ylab, label = R2_label),
    data = correlation_APET_vs_P,
    inherit.aes = FALSE
  ) +
  theme_bw() +
  labs(
    x = "Precipitation (mm)",
    y = "APET (mm)"
  ) +
  facet_wrap(~gauge, scales = "free")


ggsave(
  filename = "./Figures/Other/correlation_APET_vs_P_plot.pdf",
  plot = APET_vs_P_plot,
  device = "pdf",
  width = 1189,
  height = 841,
  units = "mm"
)


## Plot correlation on a map for the supp. =====================================
## Generate Insets =============================================================
### Filter data by state #######################################################

QLD_data <- correlation_APET_vs_P |>
  filter(state == "QLD")

NSW_data <- correlation_APET_vs_P |>
  filter(state == "NSW")

VIC_data <- correlation_APET_vs_P |>
  filter(state == "VIC")

WA_data <- correlation_APET_vs_P |>
  filter(state == "WA")

TAS_data <- correlation_APET_vs_P |>
  filter(state == "TAS")


### All colour scales must be the same #########################################
corr_range <- correlation_APET_vs_P |>
  pull(corr_P_vs_APET) |>
  range()

# round by itself does not do a good job - a single variable outside of range
corr_range <- c(
  round_any(corr_range[1], accuracy = 0.01, f = floor),
  round_any(corr_range[2], accuracy = 0.01, f = ceiling)
)

inbetween_breaks <- seq(from = corr_range[1], to = corr_range[2], length.out = 5) |> 
  round(digits = 2)

# corr_break limits cannot be rounded other risk of being NA as data is not within limits
corr_breaks <- c(corr_range[1], inbetween_breaks, corr_range[2])  # need to to show nice breaks
  

### Generate inset plots #######################################################

inset_plot_QLD <- aus_map |>
  filter(state == "QLD") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = QLD_data,
    aes(x = lon, y = lat, fill = corr_P_vs_APET),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1,
    colour = "black",
    shape = 21
  ) +
  scale_fill_distiller(
    palette = "YlOrRd",
    limits = corr_range,
    breaks = corr_breaks
  ) +
  theme_void()


inset_plot_NSW <- aus_map |>
  filter(state == "NSW") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = NSW_data,
    aes(x = lon, y = lat, fill = corr_P_vs_APET),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1,
    colour = "black",
    shape = 21
  ) +
  scale_fill_distiller(
    palette = "YlOrRd",
    limits = corr_range,
    breaks = corr_breaks
  ) +
  theme_void()


inset_plot_VIC <- aus_map |>
  filter(state == "VIC") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = VIC_data,
    aes(x = lon, y = lat, fill = corr_P_vs_APET),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1,
    colour = "black",
    shape = 21
  ) +
  scale_fill_distiller(
    palette = "YlOrRd",
    limits = corr_range,
    breaks = corr_breaks
  ) +
  theme_void()



inset_plot_WA <- aus_map |>
  filter(state == "WA") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = WA_data,
    aes(x = lon, y = lat, fill = corr_P_vs_APET),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1,
    colour = "black",
    shape = 21
  ) +
  scale_fill_distiller(
    palette = "YlOrRd",
    limits = corr_range,
    breaks = corr_breaks
  ) +
  theme_void()



inset_plot_TAS <- aus_map |>
  filter(state == "TAS") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = TAS_data,
    aes(x = lon, y = lat, fill = corr_P_vs_APET),
    show.legend = FALSE,
    size = 2.5,
    stroke = 0.1,
    colour = "black",
    shape = 21
  ) +
  scale_fill_distiller(
    palette = "YlOrRd",
    limits = corr_range,
    breaks = corr_breaks
  ) +
  theme_void()



## Put it together =============================================================

single_map_aus <- aus_map |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = correlation_APET_vs_P,
    aes(x = lon, y = lat, fill = corr_P_vs_APET),
    size = 2.5,
    stroke = 0.1,
    colour = "black",
    shape = 21
  ) +
  scale_fill_distiller(
    palette = "YlOrRd",
    limits = corr_range,
    breaks = corr_breaks
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
    fill = "Correlation between annual precipitation and PET"
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
      barwidth = unit(8, "cm")
    )
  )


single_map_aus

ggsave(
  filename = "./Figures/Supplementary/correlation_P_and_PET.pdf",
  plot = single_map_aus,
  device = "pdf",
  width = 232,
  height = 200, # 210,
  units = "mm"
)

## compare t_max vs. CO2 - probably not required
## compare t_max vs. APET - probably not required
## compare APET vs. CO2 - this is interesting







# Does the observed streamflow follow the budyko curve? ------------------------

## Filter catchment with an evidence ratio > 100 ===============================
gauges_moderately_strong_evidence_ratio <- best_CO2_non_CO2_per_gauge |>
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
  filter(evidence_ratio > 100) |>
  pull(gauge)


## Filter streamflow, precip and PET data ======================================
decade_1 <- seq(from = 1990, to = 1999)
decade_2 <- seq(from = 2012, to = 2021)

filtered_data_for_budyko <- correlations_data |>
  # We only want gauges in gauges_moderately_strong_evidence_ratio
  filter(gauge %in% gauges_moderately_strong_evidence_ratio) |>
  # Only interested in two decades
  mutate(
    decade = case_when(
      year %in% decade_1 ~ 1,
      year %in% decade_2 ~ 2,
      .default = NA
    )
  )



## AET from water balance vs. budyko ===========================================
budyko_curve <- function(P, PET) {
  sqrt(PET / P * tanh(P / PET) * (1 - exp(-PET / P)))
}


### This does not AET per changes in rainfall
AET_comparison_calc <- filtered_data_for_budyko |>
  # remove years not included in the decade_1 and 2
  drop_na() |>
  mutate(
    evapotranspiration_ratio = budyko_curve(p_mm, annual_APET_mm)
  ) |>
  summarise(
    sum_q_mm = sum(q_mm),
    sum_p_mm = sum(p_mm),
    ave_budyko_AET = mean(evapotranspiration_ratio) * mean(p_mm),
    n = n(),
    .by = c(gauge, decade)
  ) |>
  # AET
  mutate(
    sum_AET = sum_p_mm - sum_q_mm,
    ave_waterbalance_AET = sum_AET / n
  )



## Find AET differences between decades using different approaches =============
AET_comparison <- AET_comparison_calc |>
  select(gauge, decade, ave_budyko_AET, ave_waterbalance_AET) |>
  mutate(
    AET_waterbalance_minus_budyko = ave_waterbalance_AET - ave_budyko_AET
  )


ecdf_AET_comparison_function_decade_1 <- AET_comparison |>
  filter(decade == 1) |>
  pull(AET_waterbalance_minus_budyko) |>
  ecdf()

ecdf_AET_comparison_function_decade_2 <- AET_comparison |>
  filter(decade == 2) |>
  pull(AET_waterbalance_minus_budyko) |>
  ecdf()


AET_comparison <- AET_comparison |>
  mutate(
    ecdf = case_when(
      decade == 1 ~ ecdf_AET_comparison_function_decade_1(AET_waterbalance_minus_budyko),
      decade == 2 ~ ecdf_AET_comparison_function_decade_2(AET_waterbalance_minus_budyko),
      .default = NA
    )
  ) |> # change tibble to make plotting nicer
  mutate(
    Decade = if_else(decade == 1, "1990-1999", "2012-2021")
  )



AET_comparison_plot <- AET_comparison |>
  ggplot(aes(x = AET_waterbalance_minus_budyko, y = ecdf, colour = Decade)) +
  geom_step() +
  labs(
    x = bquote(Delta~"AET"), # bquote
    y = "Cumulative Probability"
  ) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.1, 0.9),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  )

ggsave(
  filename = "./Figures/Supplementary/AET_estimates_waterbalance_vs_budyko.pdf",
  plot = AET_comparison_plot,
  device = "pdf",
  width = 232,
  height = 200, # 210,
  units = "mm"
)
