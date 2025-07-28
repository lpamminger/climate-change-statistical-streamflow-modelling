# Model component analysis

# Figures produced in this R file ----------------------------------------------

# 1. Supplementary --> model_components.pdf













# CODE












# Import libraries -------------------------------------------------------------
pacman::p_load(tidyverse, ozmaps, sf)



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


best_CO2_non_CO2_per_gauge <- read_csv(
  "./Modelling/Results/CMAES/best_CO2_non_CO2_per_catchment_CMAES.csv",
  show_col_types = FALSE
)




# Generate shapefiles for map of Australia -------------------------------------

lat_lon_gauge <- gauge_information |>
  select(gauge, lat, lon)

aus_map <- generate_aus_map_sf()


# Plot best model componets by gauge -------------------------------------------


## Tidy data for plotting ======================================================
single_aus_map <- ozmaps::ozmap("country") |>
  uncount(5) |> # repeat the geometry 4 times
  mutate(
    simple_name = c("Drought", "Autocorrelation", "CO2 Intercept", "CO2 Slope", "Rainfall Seasonality")
  )

best_model_per_gauge <- best_CO2_non_CO2_per_gauge |>
  slice_min(
    AIC,
    by = gauge
  ) |>
  distinct()


parameters <- c("Autocorrelation", "CO2 Intercept", "CO2 Slope", "Drought", "Rainfall Seasonality")


gauges <- best_model_per_gauge |>
  pull(gauge) |>
  unique()


repeat_best_model_components <- tibble(
  gauge = gauges |> rep(each = length(parameters)),
  simple_name = parameters |> rep(times = length(gauges))
)


condensed_best_model_components <- best_model_per_gauge |>
  select(gauge, parameter) |>
  filter(parameter %in% c("a2", "a3_intercept", "a3_slope", "a0_d", "a4")) |>
  mutate(
    simple_name = case_when(
      parameter == "a2" ~ "Autocorrelation",
      parameter == "a3_intercept" ~ "CO2 Intercept",
      parameter == "a3_slope" ~ "CO2 Slope",
      parameter == "a0_d" ~ "Drought",
      parameter == "a4" ~ "Rainfall Seasonality",
      .default = NA
    )
  ) |>
  add_column(show = TRUE) |>
  select(!parameter)


plot_best_model_components <- repeat_best_model_components |>
  left_join(
    condensed_best_model_components,
    by = join_by(gauge, simple_name)
  ) |>
  mutate(
    show = coalesce(show, FALSE)
  ) |>
  left_join(
    lat_lon_gauge,
    by = join_by(gauge)
  )


## Count of components and plot labels =========================================
total_gauges <- best_model_per_gauge |>
  pull(gauge) |>
  unique() |>
  length()


count_components <- plot_best_model_components |>
  summarise(
    n = n(),
    .by = c(simple_name, show)
  ) |>
  filter(show) |>
  mutate(
    label = paste0(round((n / total_gauges) * 100, digits = 2), "%", " (n = ", n, ")")
  ) |>
  add_column(
    lon = 122
  ) |>
  add_column(
    lat = -40
  )


## plotting ====================================================================
model_components <- single_aus_map |>
  ggplot(aes(geometry = geometry)) +
  geom_sf(
    colour = "black",
    fill = "grey80"
  ) +
  coord_sf(xlim = c(111, 155), ylim = c(-44.5, -9.5)) +
  geom_point(
    mapping = aes(x = lon, y = lat, fill = show),
    data = plot_best_model_components,
    inherit.aes = FALSE,
    size = 1,
    show.legend = TRUE,
    shape = 21,
    colour = "black",
    stroke = 0.1
  ) +
  geom_text(
    mapping = aes(x = lon, y = lat, label = label),
    data = count_components,
    inherit.aes = FALSE,
    size = 10,
    size.unit = "pt"
  ) +
  scale_fill_brewer(palette = "Set1") +
  metR::scale_x_longitude(ticks = 7) +
  metR::scale_y_latitude(ticks = 7) +
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Best Streamflow Model Contains Component"
  ) +
  facet_wrap(~simple_name) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.825, 0.25),
    legend.background = element_rect(colour = "black"),
    axis.text = element_text(size = 7)
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 5, shape = 21, colour = "black"))
  )



ggsave(
  filename = "./Figures/Supplementary/model_components.pdf",
  plot = model_components,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)
