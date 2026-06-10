# Plot catchment boundary on a map



# Load packages ----------------------------------------------------------------
pacman::p_load(tidyverse, ozmaps, sf, ggmagnify, trend, patchwork)



# Import functions -------------------------------------------------------------
source("./Functions/utility.R")



# Load data --------------------------------------------------------------------
gauge_information <- readr::read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE
)

CAMELS_gauges <- gauge_information |> 
  filter(nested_status == "Not nested") |> 
  pull(gauge)

CAMELSAUS_boundary <- st_read(
  dsn = "./Data/Maps/02_location_boundary_area/shp/",
  as_tibble = TRUE
) |> 
  rename(gauge = CatchID) |> 
  filter(gauge %in% CAMELS_gauges)


evidence_ratio <- read_csv(
  "./Modelling/Results/CMAES/evidence_ratio_results.csv",
  show_col_types = FALSE
) |> 
# factor impact of CO2 term 
mutate(
  impact_of_CO2_term = factor(
    impact_of_CO2_term,
    levels = c("No CO2 Term", "Negative-Intercept", "Positive-Intercept", "Negative-Slope", "Positive-Slope")
  )
) |>
# factor qualitative evidence ratio values
mutate(
  binned_evidence_ratio = factor(
    binned_evidence_ratio,
    levels = c("Weak", "Moderate", "Moderately Strong", "Strong", "Very Strong", "Extremely Strong")
  )
)


australian_streamflow_network_sf <- st_read(
  dsn = "Data/Maps/Australian_streamflow_network",
  as_tibble = TRUE
)

aus_map <- generate_aus_map_sf()

aus_map_crs <- st_crs(aus_map)

working_crs <- aus_map_crs[[1]]



# Add catchment boundary to evidence ratio plot --------------------------------
## Add state column to CAMELS_AUS_boundaries for filtering =====================

CAMELSAUS_boundary <- evidence_ratio |> 
  select(gauge, state) |> 
  right_join(
    CAMELSAUS_boundary,
    by = join_by(gauge)
  ) |> 
  # make ACT = NSW
  mutate(
    state = if_else(state == "ACT", "NSW", state)
  ) |> 
  st_as_sf() |> 
  st_transform(crs = working_crs)


filter_by_state <- function(state, data) {
  data |> 
    filter(state == {{ state }}) |> 
    # ensure the sf_class is kept 
    st_as_sf()
}

states <- CAMELSAUS_boundary |> pull(state) |> unique()

state_CAMELSAUS_boundary <- map(
  .x = states,
  .f = filter_by_state,
  data = CAMELSAUS_boundary
) |> 
  `names<-`(states)


# Filter Australian Streamflow Network -----------------------------------------
## Steps:
## 1. Convert to working crs
## 2. Set what streamflow belongs to each gauge using st_contains
## 3. Filter major_streamflow_network using st_contains results
## 4. Only select major streams if catchment is above a certain size

transformed_streamflow_network <- australian_streamflow_network_sf |> 
  select(
    OBJECTID,
    Hierarchy,
    Perennial,
    Shape_Leng,
    geometry
  ) |> 
  st_transform(crs = working_crs)


sf_use_s2(FALSE) # makes it work - see https://github.com/r-spatial/sf/issues/1762
# This lists OBJECTID values per catchment boundary shapefile
streamflow_network_within_boundary <- st_contains(
  # y within x
  x = CAMELSAUS_boundary,
  y = transformed_streamflow_network,
) |> 
  # name each item in the list with CAMELSAUS_boundary gauge
  `names<-`(CAMELSAUS_boundary |> pull(gauge))


keep_OBJECTID <- streamflow_network_within_boundary |> 
  # use stack to assign each OBJECTID a gauge
  stack() |> 
  rename(
    OBJECTID = values,
    gauge = ind
  )


# Keep only the OBJECTIDs that below to CAMELS catchments of interest
transformed_streamflow_network <- keep_OBJECTID |>
  left_join(
    transformed_streamflow_network,
    by = join_by(OBJECTID)
  )
  


median_catchment_area <- gauge_information |> 
  pull(catchment_area_sq_km) |> 
  quantile(probs = 0.8)


# Split streamflow network in two - main and minor
states <- CAMELSAUS_boundary |> pull(state) |> unique()

## Main ========================================================================
major_streamflow_network <- gauge_information |> 
  select(gauge, state, catchment_area_sq_km) |> 
  right_join(
    transformed_streamflow_network,
    by = join_by(gauge)
  ) |> 
  filter(Hierarchy == "Major") |> 
  st_as_sf()


state_major_streamflow_network <- map(
  .x = states,
  .f = filter_by_state,
  data = major_streamflow_network
) |> 
  `names<-`(states)
 
## Minor =======================================================================
## TRY:
# - making minor really small linewidth (okay)
# - cutting minor based on catchment size (no)
# - shape length > X - only include long minor streams
allowable_stream_length <- gauge_information |> 
  select(gauge, state, catchment_area_sq_km) |> 
  right_join(
    transformed_streamflow_network,
    by = join_by(gauge)
  ) |> 
  filter(Hierarchy == "Minor") |> 
  pull(Shape_Leng) |> 
  quantile(probs = 0.95)

minor_streamflow_network <- gauge_information |> 
  select(gauge, state, catchment_area_sq_km) |> 
  right_join(
    transformed_streamflow_network,
    by = join_by(gauge)
  ) |> 
  filter(Hierarchy == "Minor") |> 
  filter(catchment_area_sq_km < median_catchment_area) |> 
  filter(Shape_Leng > allowable_stream_length) |> 
  st_as_sf()
  

state_minor_streamflow_network <- map(
  .x = states,
  .f = filter_by_state,
  data = minor_streamflow_network
) |> 
  `names<-`(states)


## Evidence ratio plot is copied from evidence_ratio_analysis.R ================

### Custom colour palette
custom_palette <- function(x) {
  rev(c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#f7f7f7"))
}


## Generate Insets =============================================================
### Filter data by state #######################################################

QLD_data <- evidence_ratio |>
  filter(state == "QLD")

NSW_data <- evidence_ratio |>
  filter(state == "NSW")

VIC_data <- evidence_ratio |>
  filter(state == "VIC")

WA_data <- evidence_ratio |>
  filter(state == "WA")

TAS_data <- evidence_ratio |>
  filter(state == "TAS")



### Plotting variables #########################################################
dot_size <- 1
catchment_border_linewidth <- 0.1
streamflow_network_linewidth <- 0.1
minor_streamflow_network_linewidth <- 1E-3

### Generate inset plots #######################################################
inset_plot_QLD <- aus_map |>
  filter(state == "QLD") |>
  ggplot() +
  geom_sf() +
  geom_sf(
    data = state_CAMELSAUS_boundary[["QLD"]],
    colour = "black",
    fill = NA,
    linewidth = catchment_border_linewidth
  ) +
  geom_sf(
    data = state_major_streamflow_network[["QLD"]],
    colour = "blue",
    linewidth = streamflow_network_linewidth
  ) +
  geom_sf(
    data = state_minor_streamflow_network[["QLD"]],
    colour = "blue",
    linewidth = minor_streamflow_network_linewidth
  ) +
  geom_point(
    data = QLD_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
    show.legend = FALSE,
    size = dot_size,
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
  geom_sf(
    data = state_CAMELSAUS_boundary[["NSW"]],
    colour = "black",
    fill = NA,
    linewidth = catchment_border_linewidth
  ) +
  geom_sf(
    data = state_major_streamflow_network[["NSW"]],
    colour = "blue",
    linewidth = streamflow_network_linewidth
  ) +
  geom_sf(
    data = state_minor_streamflow_network[["NSW"]],
    colour = "blue",
    linewidth = minor_streamflow_network_linewidth
  ) +
  geom_point(
    data = NSW_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
    show.legend = FALSE,
    size = dot_size,
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
  geom_sf(
    data = state_CAMELSAUS_boundary[["VIC"]],
    colour = "black",
    fill = NA,
    linewidth = catchment_border_linewidth
  ) +
  geom_sf(
    data = state_major_streamflow_network[["VIC"]],
    colour = "blue",
    linewidth = streamflow_network_linewidth
  ) +
  geom_sf(
    data = state_minor_streamflow_network[["VIC"]],
    colour = "blue",
    linewidth = minor_streamflow_network_linewidth
  ) +
  geom_point(
    data = VIC_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
    show.legend = FALSE,
    size = dot_size,
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
  geom_sf(
    data = state_CAMELSAUS_boundary[["WA"]],
    colour = "black",
    fill = NA,
    linewidth = catchment_border_linewidth
  ) +
  geom_sf(
    data = state_major_streamflow_network[["WA"]],
    colour = "blue",
    linewidth = streamflow_network_linewidth
  ) +
  geom_sf(
    data = state_minor_streamflow_network[["WA"]],
    colour = "blue",
    linewidth = minor_streamflow_network_linewidth
  ) +
  geom_point(
    data = WA_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
    show.legend = FALSE,
    size = dot_size,
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
  geom_sf(
    data = state_CAMELSAUS_boundary[["TAS"]],
    colour = "black",
    fill = NA,
    linewidth = catchment_border_linewidth
  ) +
  geom_sf(
    data = state_major_streamflow_network[["TAS"]],
    colour = "blue",
    linewidth = streamflow_network_linewidth
  ) +
  geom_sf(
    data = state_minor_streamflow_network[["TAS"]],
    colour = "blue",
    linewidth = minor_streamflow_network_linewidth
  ) +
  geom_point(
    data = TAS_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
    show.legend = FALSE,
    size = dot_size,
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
  geom_sf(
    data = CAMELSAUS_boundary,
    colour = "black",
    fill = NA,
    linewidth = catchment_border_linewidth
  ) +
  geom_sf(
    data = major_streamflow_network,
    colour = "blue",
    linewidth = streamflow_network_linewidth
  ) +
  geom_sf(
    data = minor_streamflow_network,
    colour = "blue",
    linewidth = minor_streamflow_network_linewidth
  ) +
  geom_point(
    data = evidence_ratio,
    mapping = aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
    size = dot_size + 0.5,
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


#single_map_aus

ggsave(
  filename = "./Figures/Other/evidence_ratio_aus_map_with_catchment_boundaries.pdf",
  plot = single_map_aus,
  device = "pdf",
  width = 232,
  height = 200, # 210,
  units = "mm"
)





# REPEAT FOR CHANGE IN STREAMFLOW ANALYSIS -------------------------------------
# TODO:

