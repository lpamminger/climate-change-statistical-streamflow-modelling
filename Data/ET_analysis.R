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
  summarise(
    trend_mm_per_y = get_slope(x = year, y = annual_APET_mm),
    .by = gauge
  ) |>
  left_join(
    lat_lon_gauge,
    by = join_by(gauge)
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

## compare t_max vs. CO2 - probably not required
## compare t_max vs. APET - probably not required
## compare APET vs. CO2 - this is interesting


