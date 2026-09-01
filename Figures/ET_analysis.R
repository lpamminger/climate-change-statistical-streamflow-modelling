# Examine PET
# Figures/Other/Q_PET_ratio_map.pdf
# Figures/Supplementary/AET_estimates_waterbalance_vs_budyko.pdf
# Figures/Supplementary/delta_AET_map.pdf
# Figures/Other/delta_Q_delta_AET_ratio_map.pdf


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


evidence_ratio_calc <- read_csv(
  "Modelling/Results/CMAES/evidence_ratio_results.csv",
  show_col_types = FALSE
)


master_streamflow_table <- read_csv(
  "Modelling/Results/master_streamflow_table.csv",
  show_col_types = FALSE
)


# Tidy APET data ---------------------------------------------------------------
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


# Define catchment aridity for later use ---------------------------------------
# From Han et al., 2020
# Arid --> Ep/P > 1.35
# Semi-arid --> 1 < Ep/P < 1.35
# Subhumid --> 0.76 < Ep/P < 1
# Humid --> Ep/P < 0.76

# Use total Ep/P over timeseries or use annual average? - Try both
# Once categorised comapre the evdidence ratio between the climate types using
# boxplot
define_aridity <- data |>
  left_join(
    evap_areal_potential_annual,
    by = join_by(gauge, year)
  ) |>
  mutate(
    aridity = annual_APET_mm / p_mm
  ) |>
  summarise(
    mean_aridity = mean(aridity, na.rm = TRUE),
    sum_precip = sum(p_mm, na.rm = TRUE),
    sum_APET = sum(annual_APET_mm, na.rm = TRUE),
    .by = gauge
  ) |>
  mutate(
    mean_aridity_from_sum = sum_APET / sum_precip
  ) |>
  # sort aridity into labels
  mutate(
    dryness_zone = case_when(
      mean_aridity_from_sum > 1.35 ~ "Arid",
      between(mean_aridity_from_sum, 1, 1.35) ~ "Semi-Arid",
      between(mean_aridity_from_sum, 0.76, 1) ~ "Sub-Humid",
      mean_aridity_from_sum < 0.76 ~ "Humid",
      .default = NA
    )
  ) |>
  # add evidence ratio
  left_join(
    evidence_ratio_calc,
    by = join_by(gauge)
  )

# see if the results make sense - looks good to me
generate_aus_map_sf() |>
  ggplot() +
  geom_sf() +
  geom_point(
    mapping = aes(x = lon, y = lat, colour = dryness_zone),
    data = define_aridity,
    inherit.aes = FALSE
  ) +
  theme_bw()

write_csv(
  define_aridity,
  file = "./Data/Tidy/aridity_information.csv"
)

count_labels <- define_aridity |>
  filter(evidence_ratio > 100) |>
  count(dryness_zone) |>
  mutate(
    count_label = paste0("n = ", n)
  ) |>
  mutate(
    y_lab = 50 # -25
  )


define_aridity |>
  # set the order for the boxplot labels
  mutate(
    dryness_zone = factor(dryness_zone, levels = c("Arid", "Semi-Arid", "Sub-Humid", "Humid"))
  ) |>
  filter(evidence_ratio > 100) |>
  ggplot(aes(x = dryness_zone, y = evidence_ratio)) +
  geom_boxplot(staplewidth = 0.5) +
  geom_text(
    aes(y = y_lab, label = count_label),
    data = count_labels
  ) +
  scale_y_continuous(
    transform = scales::pseudo_log_trans(base = 10),
    breaks = c(-10, 0, 10^seq(from = 0, to = 16, by = 1))
  ) +
  labs(
    x = NULL,
    y = "Evidence Ratio"
  ) +
  theme_bw()


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

## Only include gauges with evidence ratio > 100 ===============================
high_evidence_ratio_gauges <- evidence_ratio_calc |>
  filter(evidence_ratio > 100) |>
  pull(gauge)

high_evidence_master_streamflow_table <- master_streamflow_table |>
  filter(gauge %in% high_evidence_ratio_gauges)


## Calculate trends in counterfactual flow data ================================
get_slope <- function(x, y, ...) {
  lm(y ~ x, ...)$coefficients[2] |> unname() # position of slope
}

trends_in_counterfactual_flow_data <- high_evidence_master_streamflow_table |>
  select(gauge, year, realspace_streamflow_CO2_model_off, realspace_streamflow_CO2_model_on)


## Only interested in reference periods ========================================
decade_1 <- seq(from = 1990, to = 1999)
decade_2 <- seq(from = 2012, to = 2021)


timeseries_APET_vs_partitioning <- evap_areal_potential_annual |>
  left_join(
    trends_in_counterfactual_flow_data,
    by = join_by(gauge, year)
  ) |>
  filter(gauge %in% high_evidence_ratio_gauges) |>
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
    sum_a3_off = sum(realspace_streamflow_CO2_model_off),
    sum_a3_on = sum(realspace_streamflow_CO2_model_on),
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


my_palette <- function(x) {
  rev(c(
    "#a50026",
    "#d73027",
    "#f46d43",
    "#fdae61",
    "#fee090",
    "#e0f3f8",
    "#abd9e9",
    "#74add1",
    "#4575b4",
    "#313695"
  ))
}

# Map plotting function --------------------------------------------------------
map_plot <- function(plotting_variable, data, scale_range = NULL, scale_breaks = NULL, colour_palette, legend_title) {
  ## rename tibble columns for plotting
  ## This is more reliable than using braces {{ }}
  colour_palette <- noquote(colour_palette)

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
      palette = colour_palette,
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
      palette = colour_palette,
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
      palette = colour_palette,
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
      palette = colour_palette,
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
      palette = colour_palette,
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
      palette = colour_palette,
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
      # panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      legend.margin = margin(t = 5, b = 5, r = 20, l = 20, unit = "pt") # add extra padding around legend box to avoid -1.6 intersecting with line
    ) +
    guides(
      fill = guide_coloursteps(
        direction = "horizontal",
        barwidth = unit(12, "cm"),
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
  95, 0, "A"
)

decade_label_1990 <- tribble(
  ~lon, ~lat, ~label_name,
  105, 0, "1990-1999"
)

map_Q_PET_ratio_1990 <- map_plot(
  plotting_variable = Q_PET_ratio,
  data = Q_PET_ratio |> filter(decade == "1990-1999"),
  scale_range = scale_range,
  scale_breaks = c(scale_range[1], 0.01, 0.1, 0.2, 0.5, 1, 2, 5, 10, 100, scale_range[2]),
  colour_palette = my_palette,
  legend_title = bquote("abs(" * Delta * "Q [mm/y] / " * Delta * "APET [mm/y])")
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
  95, 0, "B"
)

decade_label_2012 <- tribble(
  ~lon, ~lat, ~label_name,
  105, 0, "2012-2021"
)

map_Q_PET_ratio_2012 <- map_plot(
  plotting_variable = Q_PET_ratio,
  data = Q_PET_ratio |> filter(decade == "2012-2021"),
  scale_range = scale_range,
  scale_breaks = c(scale_range[1], 0.01, 0.1, 0.2, 0.5, 1, 2, 5, 10, 100, scale_range[2]),
  colour_palette = my_palette,
  legend_title = bquote("abs(" * Delta * "Q [mm/y] / " * Delta * "APET [mm/y])")
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
  filename = "./Figures/Other/Q_PET_ratio_map.pdf",
  plot = final_plot_Q_PET_ratio,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)


# Correlations -----------------------------------------------------------------
data <- data |>
  left_join(
    evap_areal_potential_annual,
    by = join_by(gauge, year)
  ) |> # join t_max if required
  # removing missing values
  drop_na()

## compare APET vs. Precip =====================================================
correlation_APET_vs_P <- data |>
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


APET_vs_P_plot <- data |>
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
aus_map <- generate_aus_map_sf()
# this needs to use my plotting function
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


# Find difference in AET between water balance and budyko ----------------------
# Game plan:
# 1. select 3 periods with varying lengths
# 2. calculate average difference - don't worry about the uncertainty
# 3. produce 3 facets each with the different time spans

## What is the minimum and maximum year for each gauge
min_and_max_year_per_gauge <- data |>
  summarise(
    min_year = min(year),
    max_year = max(year),
    .by = gauge
  ) |>
  filter(gauge %in% high_evidence_ratio_gauges)
# Not all gauges have the required lengths...
# We are taking averages, however it may error if there is nothing to average

high_evidence_ratio_data <- data |>
  filter(gauge %in% high_evidence_ratio_gauges)


## Define time spans ===========================================================
timespan_comparison <- list(
  "single_decade" = list(decade_1, decade_2), # 1990-1999 and 2012-2021
  "double_decade" = list(seq(from = 1982, to = 2001), seq(from = 2002, to = 2021)), # 1980-1999 and 2002-2021
  "triple_decade" = list(seq(from = 1962, to = 1991), seq(from = 1992, to = 2021)) # 1992-2021
)


## Chunk data based on timespans ===============================================
chunk_data_timespans <- function(timespan, data) {
  data |>
    mutate(
      decade = case_when(
        year %in% timespan[[1]] ~ 1,
        year %in% timespan[[2]] ~ 2,
        .default = NA
      )
    ) |>
    drop_na()
}

timespan_data <- map(
  .x = timespan_comparison,
  .f = chunk_data_timespans,
  data = high_evidence_ratio_data
)


## AET from water balance vs. budyko ===========================================
budyko_curve <- function(P, PET) {
  sqrt(PET / P * tanh(P / PET) * (1 - exp(-PET / P)))
}


AET_comparison_calc <- function(chunked_data) {
  chunked_data |>
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
}


AET_comparison <- map(
  .x = timespan_data,
  .f = AET_comparison_calc
)

## Minimum number of year of data for each time span ===========================
fraction_of_acceptable <- 0.8
acceptable_lengths <- c(10, 20, 30) * fraction_of_acceptable

# remove gauges that do not meet this criteria
filter_acceptable_gauges <- function(AET_chunked_data, acceptable_lengths) {
  acceptable_gauges <- AET_chunked_data |>
    select(gauge, decade, n) |>
    pivot_wider(
      id_cols = gauge,
      names_from = decade,
      values_from = n
    ) |>
    left_join(
      min_and_max_year_per_gauge,
      by = join_by(gauge)
    ) |>
    filter(
      (`1` >= acceptable_lengths) & (`2` >= acceptable_lengths)
    ) |>
    pull(gauge)

  AET_chunked_data |>
    filter(gauge %in% acceptable_gauges)
}

## Calculate AET difference ====================================================
AET_difference <- map2(
  .x = AET_comparison,
  .y = acceptable_lengths,
  .f = filter_acceptable_gauges
) |>
  list_rbind(names_to = "source") |>
  mutate(
    AET_diff = ave_waterbalance_AET - ave_budyko_AET
  )

AET_difference |>
  summarise(
    n = n(),
    .by = c(source, decade)
  )

## Apply ecdf ==================================================================
make_ecdf <- function(x) {
  # function factory - returns a function
  ecdf_function <- ecdf(x)

  # return cdf
  return(ecdf_function(x))
}


AET_with_ecdf <- AET_difference |>
  group_by(source, decade) |>
  mutate(
    cdf = make_ecdf(AET_diff)
  )


## Plot ========================================================================
plot_AET_with_ecdf <- AET_with_ecdf |>
  # Hard coded from timespan comparison
  mutate(
    source = case_when(
      source == "single_decade" ~ "1990-1999 and 2012-2021",
      source == "double_decade" ~ "1982-2001 and 2002-2021",
      source == "triple_decade" ~ "1962-1991 and 1992-2021",
      .default = NA
    )
  ) |>
  mutate(
    decade = if_else(decade == 1, "First Timespan", "Second Timespan")
  ) |>
  ggplot(aes(x = AET_diff, y = cdf, colour = decade)) +
  geom_step() +
  scale_color_brewer(palette = "Set1") +
  labs(
    x = bquote(Delta * "AET [mm]"),
    y = "Cumulative Probability",
    colour = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.1, 0.9),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  ) +
  facet_wrap(~source)


ggsave(
  filename = "./Figures/Supplementary/varying_lengths_AET_difference_plot.pdf",
  plot = plot_AET_with_ecdf,
  device = "pdf",
  width = 232,
  height = 200, # 210,
  units = "mm"
)


## Peform statistical test between decades =====================================

### KS test
# assumptions of ks.test is the two samples are independent
# Are the difference in water balance between one decade independent to the next decade?
# Maybe? Storage? 
ks_test_results <- AET_with_ecdf |> 
  ungroup() |> 
  select(source, gauge, decade, AET_diff) |> 
  pivot_wider(
    id_cols = c(gauge, source),
    names_from = decade,
    values_from = AET_diff
  ) |>
  rename(
    AET_diff_decade_1 = `1`,
    AET_diff_decade_2 = `2`
  ) |> 
  summarise(
    ks_pvalue = ks.test(x = AET_diff_decade_1, y = AET_diff_decade_2)$p.value,
    .by = source
  ) |> 
  mutate(
    signif = ks_pvalue < 0.05
  )


### Mann U/wilcoxon test
wilcox_test_results <- AET_with_ecdf |> 
  ungroup() |> 
  select(source, gauge, decade, AET_diff) |> 
  pivot_wider(
    id_cols = c(gauge, source),
    names_from = decade,
    values_from = AET_diff
  ) |>
  rename(
    AET_diff_decade_1 = `1`,
    AET_diff_decade_2 = `2`
  ) |> 
  summarise(
    wilcox_pvalue = wilcox.test(x = AET_diff_decade_1, y = AET_diff_decade_2, paired = TRUE)$p.value,
    .by = source
  ) |> 
  mutate(
    signif = wilcox_pvalue < 0.05
  )



## Find AET differences between decades using different approaches =============
storage_error <- 0.03 # from Han et al, 2020 https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020WR027392


AET_comparison_single_comparison <- AET_comparison[["single_decade"]] |> # get the single decade for this analysis
  select(gauge, decade, ave_budyko_AET, ave_waterbalance_AET) |>
  mutate(
    AET_waterbalance_minus_budyko_ave = ave_waterbalance_AET - ave_budyko_AET,
    percent_AET_waterbalance_minus_budyko_ave = (ave_waterbalance_AET - ave_budyko_AET) / ave_budyko_AET
  ) |>
  # add uncertainty
  mutate(
    error_ave_waterbalance_AET = ave_waterbalance_AET * storage_error,
    error_ave_budyko_AET = ave_budyko_AET * storage_error,
  ) |>
  # calculate bounds
  mutate(
    lower_ave_waterbalance_AET = ave_waterbalance_AET - error_ave_waterbalance_AET,
    upper_ave_waterbalance_AET = ave_waterbalance_AET + error_ave_waterbalance_AET,
    lower_ave_budyko_AET = ave_budyko_AET - error_ave_budyko_AET,
    upper_ave_budyko_AET = ave_budyko_AET + error_ave_budyko_AET,
  ) |>
  mutate(
    # lower_bound = (lower_ave_waterbalance_AET - lower_ave_budyko_AET),
    lower_bound = (lower_ave_waterbalance_AET - upper_ave_budyko_AET),
    # upper_bound = (upper_ave_waterbalance_AET - upper_ave_budyko_AET),
    upper_bound = (upper_ave_waterbalance_AET - lower_ave_budyko_AET)
  ) |>
  # swap if upper_bound < lower_bound
  mutate(
    new_lower_bound = if_else(upper_bound > lower_bound, lower_bound, upper_bound),
    new_upper_bound = if_else(upper_bound > lower_bound, upper_bound, lower_bound),
    # all percentage changes relative to ave_budyko_AET
    percent_new_lower_bound = lower_bound / ave_budyko_AET,
    percent_new_upper_bound = upper_bound / ave_budyko_AET
  ) |>
  # remove calculation columns
  select(
    gauge,
    decade,
    new_lower_bound,
    AET_waterbalance_minus_budyko_ave,
    new_upper_bound,
    percent_new_lower_bound,
    percent_AET_waterbalance_minus_budyko_ave,
    percent_new_upper_bound
  )


## Make cumulative distribution functions ======================================
# I need 6 ecdf() function calls
# decade 1: median, lower, upper
# decade 2: median, lower, upper

# function factory?
# the function must take the AET_comparison data
# filter by decade
# call ecdf()
# return ecdf() column
# this could be a |> function


make_ecdf <- function(x) {
  # function factory - returns a function
  ecdf_function <- ecdf(x)

  # return cdf
  return(ecdf_function(x))
}

# # this returns the exact same thing - makes sense since they have the same rank

# arid_gauges <- define_aridity |>
# filter(dryness_zone == "Semi-Arid") |>
#  pull(gauge)

AET_comparison_single_comparison <- AET_comparison_single_comparison |>
  # filter(gauge %in% arid_gauges) |>
  # forces by decade operation
  group_by(decade) |>
  mutate(
    ecdf_AET_middle = make_ecdf(AET_waterbalance_minus_budyko_ave),
    ecdf_AET_upper = make_ecdf(new_upper_bound),
    ecdf_AET_lower = make_ecdf(new_lower_bound)
  ) |>
  mutate(
    decade = if_else(decade == 1, "1990-1999", "2012-2021")
  )


# it is like I need to rearrange then put back together
xmiddle <- AET_comparison_single_comparison |>
  select(decade, AET_waterbalance_minus_budyko_ave, ecdf_AET_middle) |>
  arrange(ecdf_AET_middle) |>
  rename(ecdf = ecdf_AET_middle)

xlower <- AET_comparison_single_comparison |>
  select(decade, new_lower_bound, ecdf_AET_lower) |>
  arrange(ecdf_AET_lower) |>
  rename(ecdf = ecdf_AET_lower)

xupper <- AET_comparison_single_comparison |>
  select(decade, new_upper_bound, ecdf_AET_upper) |>
  arrange(ecdf_AET_upper) |>
  rename(ecdf = ecdf_AET_upper)

plotting_AET_comparison <- xmiddle |>
  left_join(
    xlower,
    by = join_by(ecdf, decade)
  ) |>
  left_join(
    xupper,
    by = join_by(ecdf, decade)
  )


AET_comparison_plot <- plotting_AET_comparison |>
  ggplot(
    aes(x = ecdf, y = AET_waterbalance_minus_budyko_ave, colour = decade)
  ) +
  geom_step() +
  pammtools::geom_stepribbon(
    aes(ymin = new_lower_bound, ymax = new_upper_bound, fill = decade),
    alpha = 0.1,
    colour = NA
  ) +
  # flip the axis - for some reason this produces nicer results
  coord_flip() +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  labs(
    y = bquote(Delta * "AET [mm]"),
    x = "Cumulative Probability",
    colour = "Decade",
    fill = "Decade"
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.1, 0.9),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  )

AET_comparison_plot


ggsave(
  filename = "./Figures/Supplementary/3_percent_error_AET_estimates_waterbalance_vs_budyko.pdf",
  plot = AET_comparison_plot,
  device = "pdf",
  width = 232,
  height = 200, # 210,
  units = "mm"
)


# Map delta_AET ----------------------------------------------------------------
# custom range and breaks
AET_comparison_single_decade <- AET_comparison[["single_decade"]] |>
  left_join(
    lat_lon_gauge,
    by = join_by(gauge)
  ) |>
  mutate(
    AET_waterbalance_minus_budyko_ave = ave_waterbalance_AET - ave_budyko_AET
  )

scale_range <- AET_comparison_single_decade |>
  pull(AET_waterbalance_minus_budyko_ave) |>
  range()

# round by itself does not do a good job - a single variable outside of range
scale_range <- c(
  round_any(scale_range[1], accuracy = 1, f = floor),
  round_any(scale_range[2], accuracy = 1, f = ceiling)
)


### Plot 1990-1999 #############################################################
figure_label_1990 <- tribble(
  ~lon, ~lat, ~label_name,
  95, 0, "A"
)

decade_label_1990 <- tribble(
  ~lon, ~lat, ~label_name,
  105, 0, "1990-1999"
)

map_AET_ratio_1990 <- map_plot(
  plotting_variable = AET_waterbalance_minus_budyko_ave,
  data = AET_comparison_single_decade |> filter(decade == 1),
  scale_range = scale_range,
  scale_breaks = c(scale_range[1], -150, -100, -50, -25, 0, 25, 50, 100, 150, scale_range[2]),
  colour_palette = my_palette,
  legend_title = bquote("Mean" ~ Delta * "AET [mm/year]")
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


### Plot 2012-2021 #############################################################
figure_label_2012 <- tribble(
  ~lon, ~lat, ~label_name,
  95, 0, "B"
)

decade_label_2012 <- tribble(
  ~lon, ~lat, ~label_name,
  105, 0, "2012-2021"
)

map_AET_ratio_2012 <- map_plot(
  plotting_variable = AET_waterbalance_minus_budyko_ave,
  data = AET_comparison_single_decade |> filter(decade == 2),
  scale_range = scale_range,
  scale_breaks = c(scale_range[1], -150, -100, -50, -25, 0, 25, 50, 100, 150, scale_range[2]),
  colour_palette = my_palette,
  legend_title = bquote("Mean" ~ Delta * "AET [mm/year]")
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

### patchwork together and save ################################################
final_plot_delta_AET_ratio <- (map_AET_ratio_1990 | map_AET_ratio_2012) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  filename = "./Figures/Supplementary/delta_AET_map.pdf",
  plot = final_plot_delta_AET_ratio,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)


# Combine final_plot_delta_AET_ratio and plot_AET_with_ecdf for plot -----------

## Plot layout =================================================================
x <- final_plot_delta_AET_ratio / plot_AET_with_ecdf # --> look vaguely like this

ggsave(
  filename = "./Figures/Other/TESTING_AET_combined.pdf",
  plot = x,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)
## Use patchwork to join everything together
# Alternative methods attempted:
## - plot_spacer and | + / operations - too much faff
## - individually plotting timeseries - gaps

# Defining layout it the best method as I have direct control

# use the area() constructor
# top, left, bottom, right bounds (t < b and l < r)
layout <- c(
  area(t = 1, l = 1, b = 3, r = 3), # 1990s percentage change
  area(t = 1, l = 4, b = 3, r = 6), # 2010s percentage change
  area(t = 4, l = 1, b = 4, r = 6)#, # 10-year cdf
  #area(t = 4, l = 3, b = 4, r = 4), # 20-year cdf
  #area(t = 4, l = 5, b = 4, r = 6) # 30-year cdf
)

plot(layout) # check the patches are working
# the layout is correct - need to change cdf


## Need to split up the plot_AET_with_ecdf =====================================

### Make abc labels 
cde_labels <- tribble(
  ~label_name, ~source,
  "C", "single_decade",
  "D", "double_decade",
  "E", "triple_decade"
) |> 
  # same axis between facets = same x and y
  add_column(
    x_pos = -250,
    y_pos = 0.95
  ) |> 
  mutate(
    source = factor(source, levels = c("single_decade", "double_decade", "triple_decade"))
  ) 

### Plotting with facets
AET_with_ecdf_split <- AET_with_ecdf |>
  # Hard coded from timespan comparison
  #mutate(
  #  source = case_when(
  #    source == "single_decade" ~ "1990-1999 and 2012-2021",
  #    source == "double_decade" ~ "1982-2001 and 2002-2021",
  #    source == "triple_decade" ~ "1962-1991 and 1992-2021",
  #    .default = NA
  #  )
  #) |>
  mutate(
    decade = if_else(decade == 1, "First Period", "Second Period")
  ) |>
  mutate(
    source = factor(source, levels = c("single_decade", "double_decade", "triple_decade"))
  ) |> 
  ggplot(aes(x = AET_diff, y = cdf, colour = decade)) +
  geom_step() +
  geom_text(
    aes(x = x_pos, y = y_pos, label = label_name),
    data = cde_labels,
    inherit.aes = FALSE,
    fontface = "bold"
  ) +
  scale_color_brewer(palette = "Set1") +
  labs(
    x = bquote(Delta * "AET [mm]"),
    y = "Cumulative Probability",
    colour = "Comparison Period"
  ) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.1, 0.9),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black"),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) +
  facet_wrap(~source)




## Put it together =============================================================
combined_delta_AET_map <- (map_AET_ratio_1990 + map_AET_ratio_2012 + AET_with_ecdf_split) +
  plot_layout(design = layout, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black")
  )  &
 guides(colour = guide_legend(
   override.aes = list(linewidth = 2),
   title.hjust = 0.5, 
   title.position = "top", 
   ncol = 1)
   )

ggsave(
  filename = "./Figures/Main/combined_delta_AET_map.pdf",
  plot = combined_delta_AET_map,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)
