# Filtering streamflow network shape file layer to specific catchments


# Load packages ----------------------------------------------------------------
pacman::p_load(tidyverse, ozmaps, sf, ggmagnify, trend, patchwork)


# Import gauge_info ------------------------------------------------------------
gauge_information <- read_csv(
  "Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE
  ) |> 
  # only include the parent catchment (main_catchment) for nested catchments
  # with nested catchment there will be overlap 
  # I don't think overlapping shape files will play nice with `st_contains`
  filter(status == "main_catchment")


# Import functions -------------------------------------------------------------
source("./Functions/utility.R")
source("./Functions/generic_functions.R")


# Import shape files -----------------------------------------------------------
## CAMELS boundaries
CAMELSAUS_boundary <- st_read(
  dsn = "./Data/Maps/02_location_boundary_area/shp/",
  as_tibble = TRUE
) |> 
  rename(gauge = CatchID)
  
## Stream network
australian_streamflow_network_sf <- st_read(
  dsn = "./Data/Australian_Streamflow_Network/Australian_Streamflow_Network.shp",
  as_tibble = TRUE
)

## Map of Australia
aus_map <- generate_aus_map_sf()



# Convert everything to the same coordinate systems as aus_map -----------------
aus_map_crs <- st_crs(aus_map)
working_crs <- aus_map_crs[[1]]


boundary_sf <- st_transform(CAMELSAUS_boundary, crs = working_crs) |> 
  left_join(
    gauge_information,
    by = join_by(gauge)
  ) |> 
  select(gauge, geometry, state) |> 
  drop_na()


australian_streamflow_network_sf <- st_transform(
  australian_streamflow_network_sf,
  crs = working_crs
  ) |> 
  select(
    OBJECTID,
    Hierarchy,
    geometry
  )



# Filter australian_streamflow_network_sf based on boundary_sf -----------------
## Only include streamflow network within the specified catchment boundaries
sf_use_s2(FALSE) # makes it work - see https://github.com/r-spatial/sf/issues/1762
streamflow_network_within_boundary <- st_contains(
  # y within x
  x = boundary_sf,
  y = australian_streamflow_network_sf,
) |> 
  # name the gauge
  `names<-`(boundary_sf |> pull(gauge))

# output is a list of 83 items - catchment boundaries
# each item has OBJECTID from australian_streamflow_network_sf

# streamflow_network_within_boundary contains streams within a catchment

keep_OBJECTID <- streamflow_network_within_boundary |> 
  # use stack to assign each OBJECTID a gauge
  stack() |> 
  rename(
    OBJECTID = values,
    gauge = ind
  )


filtered_australian_streamflow_network_sf <- australian_streamflow_network_sf |> 
  # try this for now to make the plot nicer - TEMP
  # filter(Hierarchy == "Major") |> 
  # right join to filter out streamflow not in keep_OBJECTID
  right_join(
    keep_OBJECTID,
    by = join_by(OBJECTID)
  ) |> 
  left_join(
    gauge_information,
    by = join_by(gauge)
  ) |> 
  select(
    OBJECTID,
    Hierarchy,
    geometry,
    gauge,
    state
  )


## Save filtered_australian_streamflow_network_sf to be used elsewhere
### current it only contains main catchments
filtered_australian_streamflow_network_sf |> 
  filter(Hierarchy == "Major") |> 
  write_sf(,
  dsn = "./Data/Maps/filtered_streamflow_network_major_sf",
  driver = "ESRI Shapefile"
  )

# this takes a long time - might need another approach i.e., not saving it read and alter only
# this option is not that much code
# this is not the play
#filtered_australian_streamflow_network_sf |> 
#  filter(Hierarchy == "Minor") |> 
#  write_sf(
#  filtered_australian_streamflow_network_sf,
#  dsn = "./Data/Maps/filtered_streamflow_network_minor_sf",
#  driver = "ESRI Shapefile"
#  )

write_sf(
  boundary_sf,
  dsn = "./Data/Maps/filtered_catchment_boundaries_sf",
  driver = "ESRI Shapefile"
)

# testing if it works
#rm(filtered_australian_streamflow_network_sf)
#filtered_australian_streamflow_network_sf <- st_read("./Data/Maps/filtered_streamflow_network_sf")


### Generate inset plots #######################################################



make_inset_maps <- function(state, aus_map, boundary_sf, stream_network_sf) {
  
  state_CAMELS_boundary <- boundary_sf |> 
    filter(state == {{ state }})
  
  state_streamflow_network <- stream_network_sf |> 
    filter(state == {{ state }})
  
  aus_map |> 
    filter(state == {{ state }}) |> 
    ggplot() +
    geom_sf() +
    geom_sf(
      data = state_CAMELS_boundary) +
    geom_sf(
      data = state_streamflow_network,
      colour = "blue",
      linewidth = 0.1
      ) +
    theme_void()
}

state_names <- aus_map |> pull(state)

map_inset <- map(
  .x = state_names,
  .f = make_inset_maps,
  aus_map = aus_map,
  boundary_sf = boundary_sf,
  stream_network_sf = filtered_australian_streamflow_network_sf
) |> 
  `names<-`(state_names)



## Put it together =============================================================

nice_map_aus <- aus_map |>
  ggplot() +
  geom_sf() +
  geom_sf(
    data = boundary_sf
  ) +
  geom_sf(
    data = filtered_australian_streamflow_network_sf,
    colour = "blue",
    linewidth = 0.1
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
    plot = map_inset[["WA"]],
    proj = "single"
  ) +
  # magnify VIC
  geom_magnify(
    # aes(from = state == "VIC"), # use aes rather than manually selecting area
    from = c(141, 149.5, -39, -34),
    to = c(95, 136, -38, -60),
    shadow = FALSE,
    plot = map_inset[["VIC"]],
    proj = "single"
  ) +
  # magnify QLD
  geom_magnify(
    from = c(145, 155, -29.2, -15),
    to = c(157, 178, -29.5, 1.5),
    shadow = FALSE,
    expand = 0,
    plot = map_inset[["QLD"]],
    proj = "single"
  ) +
  # magnify NSW
  geom_magnify(
    from = c(146.5, 154, -38, -28.1),
    to = c(157, 178, -61, -30.5),
    shadow = FALSE,
    expand = 0,
    plot = map_inset[["NSW"]],
    proj = "single"
  ) +
  # magnify TAS
  geom_magnify(
    from = c(144, 149, -40, -44),
    to = c(140, 155, -45, -61),
    shadow = FALSE,
    expand = 0,
    plot = map_inset[["TAS"]],
    proj = "single"
  ) +
  labs(
    x = NULL, # "Latitude",
    y = NULL # "Longitude",
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



nice_map_aus

#ggsave(
#  filename = "./Figures/Main/TESTING_catchment_and_stream_network.pdf",
#  plot = nice_map_aus,
#  device = "pdf",
#  width = 232,
#  height = 200, # 210,
#  units = "mm"
#)
