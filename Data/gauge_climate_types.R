# TODO:
# - Import gauge data with lat and lon
# - Use the function to get climate types 
# - Save results
# - Plot on map to see if the results make sense


# Import packages --------------------------------------------------------------
pacman::p_load(tidyverse, sf, ozmaps)


# Import gauge data ------------------------------------------------------------
gauge_data <- read_csv(
  "Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE
  ) |> 
  select(gauge, station_name, state, lat, lon)


# Assign gauge a climate type using 
# 5. Add climate type to gauge_info using kgc ----------------------------------

## `LookupCZ` function provides the climate zone based on lon and lat 
## Relies on the climatezone dataframe
climatezones <- kgc::climatezones

## To use LookupCZ the data must be in |site_ID|lon|lat| format ================
formatted_gauge_data <- gauge_data |> 
  select(gauge, lon, lat) |> # mass overwrites dplyr select 
  mutate(
    rndCoord.lon = kgc::RoundCoordinates(lon),
    rndCoord.lat = kgc::RoundCoordinates(lat)
  ) 


## Gauge information with climate type =========================================
climate_type_gauge_data <- cbind(
  formatted_gauge_data, 
  "climate_type" = kgc::LookupCZ(data = formatted_gauge_data)
) |> 
  as_tibble() |> 
  mutate(
    major_climate_type = str_sub(climate_type, start = 1L, end = 1L)
  ) |> 
  # Nice names for plotting
  mutate(
    major_climate_type = case_when(
      major_climate_type == "A" ~ "Tropical (A)",
      major_climate_type == "B" ~ "Dry (B)",
      major_climate_type == "C" ~ "Temperate (C)",
      .default = major_climate_type
    )
  ) 


# See what the climate types look like on a map --------------------------------
aus_map <- ozmaps::ozmap("country")

aus_map |> 
  ggplot(aes(geometry = geometry)) +
  geom_sf(
    colour = "black",
    fill = "grey80"
  ) +
  geom_point(
    aes(x = lon, y = lat, colour = climate_type),
    data = climate_type_gauge_data,
    inherit.aes = FALSE
  ) +
  theme_bw()

# Save results 
write_csv(climate_type_gauge_data, "Data/Tidy/gauge_information_with_climate_data_CAMELS.csv")
