# Examine PET


# Import libraries required ----------------------------------------------------
pacman::p_load(tidyverse, ozmaps, sf, ggmagnify, patchwork)


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

trends_in_counterfactual_flow_data <- read_csv(
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

## Only interested in reference periods ========================================
decade_1 <- seq(from = 1990, to = 1999)
decade_2 <- seq(from = 2012, to = 2021)

get_slope <- function(x, y, ...) {
  lm(y ~ x, ...)$coefficients[2] |> unname() # position of slope
}

timeseries_APET_vs_partitioning <- evap_areal_potential_annual |>
  left_join(
    trends_in_counterfactual_flow_data,
    by = join_by(gauge, year)
  ) |>
  filter(gauge %in% gauges_moderately_strong_evidence_ratio) |>
  # add decade
  mutate(
    decade = case_when(
      year %in% decade_1 ~ "1990-1999",
      year %in% decade_2 ~ "2012-2021",
      .default = "other"
    )
  ) |>
  # remove other decade - don't remove NA's yet
  filter(decade != "other")


## Do APET trends and CO2-partitioning trends separately due to missing streamflow data years ====
APET_trends <- timeseries_APET_vs_partitioning |>
  summarise(
    APET_trend_mm_per_y = get_slope(x = year, y = annual_APET_mm),
    .by = c(gauge, decade)
  )

CO2_partitioning_trends <- timeseries_APET_vs_partitioning |>
  # remove missing years
  drop_na() |>
  summarise(
    sum_a3_off = sum(realspace_a3_off_streamflow),
    sum_a3_on = sum(realspace_a3_on_streamflow),
    n = n(),
    .by = c(gauge, decade)
  ) |>
  # find impact of CO2 partitioning
  mutate(
    CO2_partitioning_trend_mm_per_y = (sum_a3_off - sum_a3_on) / n
  )

Q_PET_ratio <- APET_trends |>
  left_join(
    CO2_partitioning_trends,
    by = join_by(gauge, decade)
  ) |>
  mutate(
    Q_PET_ratio = abs(CO2_partitioning_trend_mm_per_y) / abs(APET_trend_mm_per_y)
  ) |>
  left_join(
    lat_lon_gauge,
    by = join_by(gauge)
  )

Q_PET_ratio |>
  summarise(
    mean = mean(Q_PET_ratio),
    max = max(Q_PET_ratio),
    min = min(Q_PET_ratio),
    .by = decade
  )

#  Q_PET_ratio needs to be log-scaled
Q_PET_ratio |>
  pull(Q_PET_ratio) |>
  sort()
# c(0, 0.01, 0.1, 0, 1, 10, 100, 1000)

my_palette <- function(x) {
  rev(c(
    "#d73027",
    "#fc8d59",
    "#fee090",
    "#e0f3f8",
    "#91bfdb",
    "#4575b4"
  ))
}

# Map plotting function --------------------------------------------------------
map_plot <- function(plotting_variable, data, scale_range = NULL, scale_breaks = NULL, colour_palette, legend_title) {
  ## rename tibble columns for plotting
  ## This is more reliable than using braces {{ }}
  data <- data |>
    rename(
      plotting_variable = {{ plotting_variable }}
    )


  ## Make map template using ozmaps ============================================
  aus_map <- generate_aus_map_sf()

  ## Get inset data ============================================================
  ### Filter data by state #####################################################
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


  ## Generate same scale for all plots =========================================
  if (is.null(scale_range) | is.null(scale_breaks)) {
    scale_range <- data |>
      pull(plotting_variable) |>
      range()

    # round by itself does not do a good job - a single variable outside of range
    scale_range <- c(
      round_any(scale_range[1], accuracy = 0.01, f = floor),
      round_any(scale_range[2], accuracy = 0.01, f = ceiling)
    )

    scale_breaks <- seq(from = scale_range[1], to = scale_range[2], length.out = 5L) |> # this could be a function argument
      round(digits = 2)

    scale_breaks <- c(scale_range[1], scale_breaks[-c(1, length(scale_breaks))], scale_range[2]) # need to to show nice breaks

    # need a default palette and a check to see if the palette length matches breaks
  }




  ## Generate inset plots ======================================================
  inset_plot_QLD <- aus_map |>
    filter(state == "QLD") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = QLD_data,
      aes(x = lon, y = lat, fill = plotting_variable),
      show.legend = FALSE,
      size = 2.5,
      stroke = 0.1,
      colour = "black",
      shape = 21
    ) +
    binned_scale(
      aesthetics = "fill",
      palette = my_palette,
      breaks = scale_breaks,
      limits = scale_range,
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
      aes(x = lon, y = lat, fill = plotting_variable),
      show.legend = FALSE,
      size = 2.5,
      stroke = 0.1,
      colour = "black",
      shape = 21
    ) +
    binned_scale(
      aesthetics = "fill",
      palette = my_palette,
      breaks = scale_breaks,
      limits = scale_range,
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
      aes(x = lon, y = lat, fill = plotting_variable),
      show.legend = FALSE,
      size = 2.5,
      stroke = 0.1,
      colour = "black",
      shape = 21
    ) +
    # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
    binned_scale(
      aesthetics = "fill",
      palette = my_palette,
      breaks = scale_breaks,
      limits = scale_range,
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
      aes(x = lon, y = lat, fill = plotting_variable),
      show.legend = FALSE,
      size = 2.5,
      stroke = 0.1,
      colour = "black",
      shape = 21
    ) +
    binned_scale(
      aesthetics = "fill",
      palette = my_palette,
      breaks = scale_breaks,
      limits = scale_range,
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
      aes(x = lon, y = lat, fill = plotting_variable),
      show.legend = FALSE,
      size = 2.5,
      stroke = 0.1,
      colour = "black",
      shape = 21
    ) +
    binned_scale(
      aesthetics = "fill",
      palette = my_palette,
      breaks = scale_breaks,
      limits = scale_range,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    theme_void()


  ## The big map ===============================================================
  aus_map |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = data,
      aes(x = lon, y = lat, fill = plotting_variable),
      size = 2.5,
      stroke = 0.1,
      colour = "black",
      shape = 21
    ) +
    binned_scale(
      aesthetics = "fill",
      palette = my_palette,
      breaks = scale_breaks,
      limits = scale_range,
      show.limits = TRUE,
      guide = "colorsteps"
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
      fill = legend_title
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
      fill = guide_coloursteps(
        direction = "horizontal",
        barwidth = unit(8, "cm"),
        even.steps = TRUE
      )
    )
}


# custom range and breaks
scale_range <- Q_PET_ratio |>
  pull(Q_PET_ratio) |>
  range()

# round by itself does not do a good job - a single variable outside of range
scale_range <- c(
  round_any(scale_range[1], accuracy = 0.01, f = floor),
  round_any(scale_range[2], accuracy = 0.01, f = ceiling)
)

## Plot 1990-1999 ==============================================================
figure_label_1990 <- tribble(
  ~lon, ~lat, ~label_name,
  95,   0,   "A"
)

decade_label_1990 <- tribble(
  ~lon, ~lat, ~label_name,
  105,   0,   "1990-1999"
)

map_Q_PET_ratio_1990 <- map_plot(
  plotting_variable = Q_PET_ratio,
  data = Q_PET_ratio |> filter(decade == "1990-1999"),
  scale_range = scale_range,
  scale_breaks = c(scale_range[1], 0.01, 0.1, 1, 10, 100, scale_range[2]),
  colour_palette = my_palette,
  legend_title = bquote(Delta ~ "Q/APET")
) +
  geom_text(
    data = figure_label_1990,
    aes(x = lon, y = lat, label = label_name),
    fontface = "bold",
    size = 10,
    size.unit = "pt"
  ) +
  geom_text(
    data = decade_label_1990,
    aes(x = lon, y = lat, label = label_name),
    size = 10,
    size.unit = "pt"
  )


## Plot 2012-2021 ==============================================================
figure_label_2012 <- tribble(
  ~lon, ~lat, ~label_name,
  95,   0,   "B"
)

decade_label_2012 <- tribble(
  ~lon, ~lat, ~label_name,
  105,   0,   "2012-2021"
)

map_Q_PET_ratio_2012 <- map_plot(
  plotting_variable = Q_PET_ratio,
  data = Q_PET_ratio |> filter(decade == "2012-2021"),
  scale_range = scale_range,
  scale_breaks = c(scale_range[1], 0.01, 0.1, 1, 10, 100, scale_range[2]),
  colour_palette = my_palette,
  legend_title = bquote(Delta ~ "Q/APET")
) +
  geom_text(
    data = figure_label_2012,
    aes(x = lon, y = lat, label = label_name),
    fontface = "bold",
    size = 10,
    size.unit = "pt"
  ) +
  geom_text(
    data = decade_label_2012,
    aes(x = lon, y = lat, label = label_name),
    size = 10,
    size.unit = "pt"
  )

## patchwork together and save =================================================
final_plot_Q_PET_ratio <- (map_Q_PET_ratio_1990 | map_Q_PET_ratio_2012) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  filename = "./Figures/Supplementary/Q_PET_ratio_map.pdf",
  plot = final_plot_Q_PET_ratio,
  device = "pdf",
  width = 297,
  height = 210,
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

### get mean and sd correlations
correlation_APET_vs_P |>
  summarise(
    mean_corr = mean(corr_P_vs_APET),
    sd_corr = sd(corr_P_vs_APET),
    .by = state
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
corr_breaks <- c(corr_range[1], inbetween_breaks, corr_range[2]) # need to to show nice breaks


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


## Filter streamflow, precip and PET data ======================================
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
    x = bquote(Delta ~ "AET"), # bquote
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
