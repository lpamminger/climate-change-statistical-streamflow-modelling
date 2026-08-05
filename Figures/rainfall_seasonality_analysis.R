# Questions to answer:
## 1. Is there a relationship between summer and winter dominated rainfall
##    to evidence ratio?
## 2. Do catchments that exhibit strong rainfall seasonality have a
##    seasonality model component?


# What do I need?
# - evidence ratio
# - a way to determine winter/summer dominated rainfall (Koppen zones?)
#   Koppen has dry winter/dry summer for temperate catchments
#   I have this information in gauge_information.R
# - OR do I look at the precipitation data i.e., look at summer months
#   look at winter months


# Import libraries -------------------------------------------------------------
pacman::p_load(tidyverse, ozmaps, sf, ggmagnify, trend, patchwork)

# Import functions -------------------------------------------------------------
source("./Functions/utility.R")
# source("./Functions/streamflow_models.R")
# source("./Functions/parameter_transformations.R")
# source("./Functions/catchment_data_blueprint.R")
# source("./Functions/cmaes_dream_summaries.R")
# source("./Functions/objective_functions.R")
# source("./Functions/numerical_optimiser_setup.R")
# source("./Functions/generic_functions.R")
# source("./Functions/DREAM.R")
# source("./Functions/objective_function_setup.R")
# source("./Functions/result_set.R")
# source("./Functions/boxcox_logsinh_transforms.R")


# Import data ------------------------------------------------------------------
gauge_information_CAMELS <- read_csv(
  "Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE
)


daily_precipitation <- read_csv(
  "Data/Raw/precipitation_AGCD.csv",
  show_col_types = FALSE
)


evidence_ratio <- read_csv(
  "Modelling/Results/CMAES/evidence_ratio_results.csv",
  show_col_types = FALSE
) |>
  mutate(
    binned_evidence_ratio = factor(
      binned_evidence_ratio,
      levels = c("Weak", "Moderate", "Moderately Strong", "Strong", "Very Strong", "Extremely Strong")
    )
  ) |>
  mutate(
    impact_of_CO2_term = factor(
      impact_of_CO2_term,
      levels = c("No CO2 Term", "Negative-Intercept", "Positive-Intercept", "Negative-Slope", "Positive-Slope")
    )
  )


# Wrangle precipitation data into seasonal averages ----------------------------
monthly_precipitation <- daily_precipitation |>
  pivot_longer(
    cols = !c(year:day),
    names_to = "gauge",
    values_to = "p_mm"
  ) |>
  summarise(
    p_mm = sum(p_mm),
    .by = c(year, month, gauge)
  )


warm_season <- c(11, 12, 1, 2, 3, 4) # Nov to Apr
cool_season <- 5:10 # May to Oct

# calculate seasonal ratio = warm season median rainfall / cool season median rainfall (annual rainfall)
# ratio of median rainfalls

seasonal_precipitation <- monthly_precipitation |>
  mutate(
    season = case_when(
      month %in% warm_season ~ "warm_season",
      month %in% cool_season ~ "cool_season",
      .default = NA
    )
  ) |>
  summarise(
    seasonal_p_mm = sum(p_mm),
    .by = c(gauge, year, season)
  )


annual_precipitation <- seasonal_precipitation |>
  summarise(
    median_annual_p_mm = median(seasonal_p_mm),
    .by = gauge
  )


seasonal_rainfall_classification <- seasonal_precipitation |>
  pivot_wider(
    id_cols = c(gauge, year),
    names_from = season,
    values_from = seasonal_p_mm
  ) |>
  # seasonal ratio
  mutate(
    seasonal_ratio_per_year = warm_season / cool_season
  ) |>
  # take the median
  summarise(
    seasonal_ratio = median(seasonal_ratio_per_year),
    .by = gauge
  ) |>
  # add annual rainfall
  left_join(
    annual_precipitation,
    by = join_by(gauge)
  ) |>
  # classification based on BOM
  # (table 2 of https://www.bom.gov.au/climate/maps/averages/climate-classification/files/climate-classification-maps-technical-details.pdf)
  mutate(
    classification = case_when(
      (seasonal_ratio > 3) & (median_annual_p_mm > 350) ~ "Summer Dominant",
      (seasonal_ratio > 1.3) & (seasonal_ratio <= 3) & (median_annual_p_mm > 350) ~ "Summer",
      (seasonal_ratio > 1 / 1.3) & (seasonal_ratio <= 1.3) & (median_annual_p_mm > 250) ~ "Uniform",
      (seasonal_ratio >= 1 / 3) & (seasonal_ratio < 1 / 1.3) & (median_annual_p_mm > 250) ~ "Winter",
      (seasonal_ratio < 1 / 3) & (median_annual_p_mm > 250) ~ "Winter Dominant",
      (seasonal_ratio >= 1.3) & (median_annual_p_mm < 350) ~ "Arid",
      (seasonal_ratio < 1.3) & (median_annual_p_mm < 250) ~ "Arid",
      .default = NA
    )
  )


## Combine with Evidence ratio =================================================
seasonal_rainfall_classification_gauge_info <- seasonal_rainfall_classification |>
  right_join(
    evidence_ratio,
    by = join_by(gauge)
  ) |>
  select(
    gauge,
    classification,
    evidence_ratio,
    binned_evidence_ratio,
    impact_of_CO2_term,
    state,
    lat,
    lon
  ) |> 
  # add levels to classification
  mutate(
    classification = factor(classification, levels = c("Summer Dominant", "Summer", "Uniform", "Winter", "Winter Dominant", "Arid"))
  )


# Plotting ---------------------------------------------------------------------
## Checking classification is correct ==========================================
aus_map <- generate_aus_map_sf()


aus_map |>
  ggplot() +
  geom_sf() +
  geom_point(
    mapping = aes(x = lon, y = lat, fill = classification),
    data = seasonal_rainfall_classification_gauge_info,
    shape = 21
  ) +
  theme_bw()


## Box plot ====================================================================
evidence_ratio_filter <- 100

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
  geom_boxplot(staplewidth = 0.5) +
  geom_text(
    aes(x = classification, y = y_pos, label = n_label),
    data = count_occurences,
    vjust = 1.5
  ) +
  scale_y_continuous(
    transform = scales::pseudo_log_trans(base = 10),
    breaks = c(-10, 0, 10^seq(from = 0, to = 16, by = 1))
  ) +
  scale_x_discrete(drop = FALSE) +
  labs(
    x = "Seasonal Rainfall Major Zones",
    y = "Evidence Ratio"
  ) +
  theme_bw()

ggsave(
  filename = "Figures/Other/seasonal_high_evidence_ratio_boxplot.pdf",
  plot = seasonal_evidence_ratio_boxplot,
  device = "pdf",
  width = 232,
  height = 200, # 210,
  units = "mm"
)


## Box plot for streamflow percentage changes ==================================
CO2_streamflow_percentage_changes <- read_csv(
  "Modelling/Results/CO2_streamflow_percentage_changes.csv",
  show_col_types = FALSE
  )

seasonal_CO2_streamflow_percentage_changes <- seasonal_rainfall_classification_gauge_info |> 
  select(gauge, classification) |> 
  right_join(
    CO2_streamflow_percentage_changes,
    by = join_by(gauge)
  )

count_occurences <- count_occurences |> 
  mutate(
    y_pos = -75
  )

seasonal_CO2_streamflow_percentage_changes_plot <- seasonal_CO2_streamflow_percentage_changes |> 
  mutate(
    decade = if_else(decade == 1, "1990-1999", "2012-2021") 
  ) |> 
  filter(decade == "2012-2021") |> 
  ggplot(aes(x = classification, y = CO2_impact_on_streamflow_percent)) +
  geom_boxplot(staplewidth = 0.5) +
  geom_text(
    aes(x = classification, y = y_pos, label = n_label),
    data = count_occurences
  ) +
  scale_x_discrete(drop = FALSE) +
  labs(
    x = "Seasonal Rainfall Major Zones",
    y = "CO2 Impact on Streamflow (%)"
  ) +
  theme_bw() +
  facet_wrap(~decade, ncol = 2, scales = "fixed")

ggsave(
  filename = "Figures/Other/seasonal_CO2_streamflow_percentage_changes.pdf",
  plot = seasonal_CO2_streamflow_percentage_changes_plot,
  device = "pdf",
  width = 232,
  height = 200, # 210,
  units = "mm"
)


## Evidence ratio plot =========================================================
### Custom colour palette
custom_palette <- function(x) {
  rev(c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#f7f7f7"))
}


## Generate Insets =============================================================
### Filter data by state #######################################################
QLD_data <- seasonal_rainfall_classification_gauge_info |>
  filter(state == "QLD")

NSW_data <- seasonal_rainfall_classification_gauge_info |>
  filter(state == "NSW")

VIC_data <- seasonal_rainfall_classification_gauge_info |>
  filter(state == "VIC")

WA_data <- seasonal_rainfall_classification_gauge_info |>
  filter(state == "WA")

TAS_data <- seasonal_rainfall_classification_gauge_info |>
  filter(state == "TAS")


### Generate inset plots #######################################################

inset_plot_QLD <- aus_map |>
  filter(state == "QLD") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = QLD_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term, colour = classification),
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
  scale_colour_brewer(
    palette = "Set1",
    drop = FALSE
  ) +
  theme_void()


inset_plot_NSW <- aus_map |>
  filter(state == "NSW") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = NSW_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term, colour = classification),
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
  scale_colour_brewer(
    palette = "Set1",
    drop = FALSE
  ) +
  theme_void()


inset_plot_VIC <- aus_map |>
  filter(state == "VIC") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = VIC_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term, colour = classification),
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
  scale_colour_brewer(
    palette = "Set1",
    drop = FALSE
  ) +
  theme_void()


inset_plot_WA <- aus_map |>
  filter(state == "WA") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = WA_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term, colour = classification),
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
  scale_colour_brewer(
    palette = "Set1",
    drop = FALSE
  ) +
  theme_void()


inset_plot_TAS <- aus_map |>
  filter(state == "TAS") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = TAS_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term, colour = classification),
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
  scale_colour_brewer(
    palette = "Set1",
    drop = FALSE
  ) +
  theme_void()


## Put it together =============================================================

single_map_aus <- aus_map |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = seasonal_rainfall_classification_gauge_info,
    mapping = aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term, colour = classification),
    size = 3,
    stroke = 0.1
  ) +
  theme_bw() +
  scale_fill_manual(
    values = custom_palette(),
    drop = FALSE
  ) +
  scale_colour_brewer(
    palette = "Set1",
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
    shape = bquote("Impact of" ~ CO[2] ~ "Term"),
    colour = "Seasonal Rainfall"
  ) +
  theme(
    # legend.key = element_rect(fill = "grey80"),
    legend.title = element_text(hjust = 0.5, size = 8),
    legend.text = element_text(size = 6),
    legend.background = element_rect(colour = "black"),
    legend.margin = margin(2.5, 2.5, 2.5, 2.5),
    axis.text = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.37, 0.9), # constants used to move the legend in the right place
    legend.box = "horizontal", # side-by-side legends
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    legend.spacing = unit(0.1, "cm")
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 5, shape = 21), nrow = 3), # Wrap legend with nrow
    shape = guide_legend(override.aes = list(size = 5, fill = "grey50"), nrow = 3),
    colour = guide_legend(override.aes = list(size = 2, stroke = 2, shape = 21, fill = "grey50"), nrow = 3)
  )


ggsave(
  filename = "./Figures/Other/seasonal_rainfall_evidence_ratio_aus_map.pdf",
  plot = single_map_aus,
  device = "pdf",
  width = 232,
  height = 200, # 210,
  units = "mm"
)



# Rainfall seasonality vs. best model components -------------------------------

## Repeat seasonality calc using our values rather than BoMs ===================
# Our seasonality are slightly different to BoMs
# Our cool Apr - Sept --> BoM cool May - Oct 
# Our warm Oct - March --> BoM warm Nov - Apr
# use our seasonality metric because our models used them

our_warm_season <- c(10, 11, 12, 1, 2, 3) 
our_cool_season <- 4:9 

# calculate seasonal ratio = warm season median rainfall / cool season median rainfall (annual rainfall)
# ratio of median rainfalls

our_seasonal_precipitation <- monthly_precipitation |>
  mutate(
    season = case_when(
      month %in% our_warm_season ~ "warm_season",
      month %in% our_cool_season ~ "cool_season",
      .default = NA
    )
  ) |>
  summarise(
    seasonal_p_mm = sum(p_mm),
    .by = c(gauge, year, season)
  ) |> 
  pivot_wider(
    id_cols = c(gauge, year),
    names_from = season,
    values_from = seasonal_p_mm
  ) |>
  # seasonal ratio
  mutate(
    seasonal_ratio_per_year = warm_season / cool_season
  )

our_median_seasonal_ratio <- our_seasonal_precipitation |> 
  summarise(
    median_seasonal_ratio = median(seasonal_ratio_per_year),
    .by = gauge
  )


## Identify catchment that include seasonality component in best model ========= 
best_CO2_non_CO2 <- read_csv(
  "Modelling/Results/CMAES/best_CO2_non_CO2_per_catchment_CMAES.csv",
  show_col_types = FALSE
  )

contains_seasonality_term <- function(parameter_names) {
  any(parameter_names %in% "a4")
}

best_model_with_seasonality_term <- best_CO2_non_CO2 |>
  slice_min(
    AIC,
    by = gauge
  ) |> 
  summarise(
    contains_seasonality_term = contains_seasonality_term(parameter),
    .by = gauge
  ) # add the seasonality ratio (median warm to cool)
  

## Compare =====================================================================
compare_model_and_observed_seasonality <- best_model_with_seasonality_term |> 
  left_join(
    our_median_seasonal_ratio,
    by = join_by(gauge)
  )

# join the BoM labels - I acknowledge they use different months
compare_model_and_observed_seasonality <- seasonal_rainfall_classification_gauge_info |> 
  select(gauge, classification, evidence_ratio) |> 
  right_join(
    compare_model_and_observed_seasonality,
    by = join_by(gauge)
  )

# I think we count them 
# the seasonal ratio doesn't show winter dominated rainfall well 1/0.8
# count number of TRUE/FALSE per classification

compare_model_and_observed_seasonality |> 
  #filter(evidence_ratio > 100) |> 
  summarise(
    model_contains_seasonality_term = sum(contains_seasonality_term),
    n = n(),
    .by = classification
  ) |> 
  mutate(
    percentage = round(model_contains_seasonality_term / n, digits = 2) * 100
  )


compare_model_and_observed_seasonality |>
  count(contains_seasonality_term, classification) |> 
  pivot_wider(
    names_from = contains_seasonality_term,
    values_from = n
  ) |> 
  mutate(
    total = (`TRUE` + `FALSE`),
    contains_percentages = `TRUE` / (`TRUE` + `FALSE`)
  )

proportion_of_seasonality_per_model <- compare_model_and_observed_seasonality |>
  count(contains_seasonality_term, classification) |> 
  # add levels to classification
  mutate(
    classification = factor(classification, levels = c("Arid", "Summer Dominant", "Summer", "Uniform", "Winter", "Winter Dominant"))
  ) |>   
  ggplot(aes(x = classification, y = n, fill = contains_seasonality_term)) +
  geom_bar(position = "fill", stat = "identity", colour = "black") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Seasonal Rainfall Zones",
    y = "Proportion",
    fill = "Contains Seasonality Term"
  ) +
  scale_fill_manual(values = c("grey95", "grey50")) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(
  filename = "./Figures/Supplementary/proportion_of_seasonality_per_model.pdf",
  plot = proportion_of_seasonality_per_model,
  device = "pdf",
  width = 180,
  height = 140, # 210,
  units = "mm"
)










