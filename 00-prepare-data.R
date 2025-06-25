# Preparing CAMEL data

# Clear environment and console ------------------------------------------------
rm(list = ls())
cat("\014")

# Import libraries -------------------------------------------------------------
pacman::p_load(tidyverse, checkmate)


# Functions for tidying --------------------------------------------------------
source("./Functions/utility.R")



# Import data ------------------------------------------------------------------

daily_streamflow_raw <- readr::read_csv(
  "./Data/Raw/streamflow_mmd.csv",
  na = c("-99.99"), # NA is represented using -99.99 in CAMELS dataset. Convert to NA
  show_col_types = FALSE # quiets column specification
)

daily_precip_raw <- readr::read_csv(
  "./Data/Raw/precipitation_AGCD.csv",
  na = c("-99.99"),
  show_col_types = FALSE
)

yearly_CO2_ppm <- readr::read_csv(
  "./Data/Raw/20241125_Mauna_Loa_CO2_data.csv",
  skip = 43,
  col_select = !unc,
  show_col_types = FALSE
)

catchment_information <- readr::read_csv(
  "./Data/Raw/CAMELS_AUS_Attributes&Indices_MasterTable.csv",
  col_select = c(
    "station_id",
    "station_name",
    "state_outlet",
    "lat_outlet",
    "long_outlet"
  ),
  show_col_types = FALSE
) |> 
  rename(
    gauge = station_id,
    lat = lat_outlet,
    lon = long_outlet,
    state = state_outlet
  )




# Tidying data -----------------------------------------------------------------

## CONSTANTS ===================================================================
acceptable_missing_streamflow_days <- 10
min_record_length_years <- 30 # at least 30 years of data required. Copied from HRS
min_run_length <- 2
pre_ind_CO2_ppm <- 280


## Get ratio of warm season rainfall over cool season rainfall =================
monthly_rainfall <- daily_precip_raw |>
  pivot_longer(
    cols = !c(year:day),
    names_to = "gauge",
    values_to = "p_mm"
  ) |>
  summarise(
    p_mm = sum(p_mm),
    .by = c(year, month, gauge)
  ) |>
  arrange(gauge)


cool_seasons <- seq(from = 4, to = 9, by = 1)

warm_seasons <- c(10, 11, 12, 1, 2, 3)

### For each year sum the cool season rainfall and the warm season rainfall
### We assume cool season is between April (#4) to Sept (#9) and
### warm season is between Oct (#10) to Mar (#3)

warm_to_cool_season_rainfall_ratio <- function(monthly_rainfall, month, cool_seasons, warm_seasons) {
  # Check if a given year has all the months
  if (any(sort(month) != sort(c(cool_seasons, warm_seasons)))) {
    message("Missing months")
  }

  sum_warm_seasons_rainfall <- sum(monthly_rainfall[month %in% warm_seasons])
  sum_cool_seasons_rainfall <- sum(monthly_rainfall[month %in% cool_seasons])

  warm_to_cool_season_rainfall_ratio <- sum_warm_seasons_rainfall / sum_cool_seasons_rainfall
}


warm_season_to_annual_rainfall_ratio <- function(monthly_rainfall, month, warm_seasons) {
  sum_warm_seasons_rainfall <- sum(monthly_rainfall[month %in% warm_seasons])
  annual_rainfall <- sum(monthly_rainfall)

  warm_season_to_annual_rainfall_ratio <- sum_warm_seasons_rainfall / annual_rainfall
}





### annual_seasonal_rainfall_ratio #############################################
annual_seasonal_rainfall_ratio <- monthly_rainfall |>
  summarise(
    warm_to_cool_season_rainfall_ratio = warm_to_cool_season_rainfall_ratio(
      monthly_rainfall = p_mm,
      month = month,
      cool_seasons = cool_seasons,
      warm_seasons = warm_seasons
    ),
    warm_season_to_annual_rainfall_ratio = warm_season_to_annual_rainfall_ratio(
      monthly_rainfall = p_mm,
      month = month,
      warm_seasons = warm_seasons
    ),
    .by = c(year, gauge)
  ) |>
  mutate(
    warm_to_cool_season_rainfall_ratio = na_if(
      warm_to_cool_season_rainfall_ratio, Inf
    )
  ) |>
  mutate(
    mean_warm_to_cool_season_rainfall_ratio = mean(warm_to_cool_season_rainfall_ratio, na.rm = TRUE),
    mean_warm_season_to_annual_rainfall_ratio = mean(warm_season_to_annual_rainfall_ratio, na.rm = TRUE),
    .by = gauge
  ) |>
  arrange(gauge) |>
  mutate(
    standardised_warm_to_cool_season_rainfall_ratio = warm_to_cool_season_rainfall_ratio - mean_warm_to_cool_season_rainfall_ratio,
    standardised_warm_season_to_annual_rainfall_ratio = warm_season_to_annual_rainfall_ratio - mean_warm_season_to_annual_rainfall_ratio
  ) |>
  select(
    c(
      year,
      gauge,
      standardised_warm_to_cool_season_rainfall_ratio,
      standardised_warm_season_to_annual_rainfall_ratio
    )
  )



## Yearly data =================================================================
yearly_precip <- daily_precip_raw |> # The entire precip data set is continuous
  pivot_longer(
    cols = !c(year:day),
    names_to = "gauge",
    values_to = "p_mm"
  ) |>
  summarise(
    p_mm = sum(p_mm),
    .by = c(year, gauge)
  ) |>
  arrange(gauge, year) |>
  left_join(
    annual_seasonal_rainfall_ratio,
    by = join_by(gauge, year)
  )


yearly_streamflow <- daily_streamflow_raw |>
  pivot_longer(
    cols = !c(year:day),
    names_to = "gauge",
    values_to = "q_mm"
  ) |>
  mutate(is_na = is.na(q_mm)) |>
  summarise(
    q_mm = sum(q_mm, na.rm = TRUE), # add up everything. Ignore missing values
    q_na_count = sum(is_na),
    .by = c(year, gauge)
  ) |>
  mutate(
    q_mm = if_else(q_na_count >= acceptable_missing_streamflow_days, NA, q_mm) # if there are more than x missing values then set to NA
  ) |>
  arrange(gauge, year)


yearly_data <- yearly_precip |>
  right_join(yearly_streamflow, by = join_by(gauge, year)) |>
  left_join(yearly_CO2_ppm, by = join_by(year)) |> # right join makes one gauge = NA
  rename(CO2 = mean) |>
  filter(!is.na(CO2)) |>
  mutate(
    drought = saft_drought_algorithm(p_mm),
    .by = gauge
  )


# Check if streamflow is greater than rainfall ---------------------------------
is_streamflow_greater_than_rainfall <- yearly_data |>
  mutate(
    q_greater_than_p = q_mm > p_mm,
  ) |> 
  filter(q_greater_than_p) |> 
  # how much greater?
  mutate(
    magnitude_q_greater_than_p = q_mm - p_mm
  ) |> 
  arrange(desc(magnitude_q_greater_than_p))


count_gauge_is_streamflow_greater_than_rainfall <- is_streamflow_greater_than_rainfall |> 
  summarise(
    n = n(),
    .by = gauge
  ) |> 
  arrange(desc(n))


# Remove from data
yearly_data <- yearly_data |> 
  anti_join(
    is_streamflow_greater_than_rainfall,
    by = join_by(gauge, year, p_mm)
  )


# Tom's slice function ---------------------------------------------------------
## Find the start and end of all continous streamflow measurements.
## Return a n x m matrix where m[1] = start_index, m[2] = end_index


gauge_continous_start_end <- function(single_gauge, data, min_run_length) {
  gauge_specific_data <- data |>
    filter(single_gauge == gauge)

  start_end_matrix <- continuous_run_start_end(gauge_specific_data$q_mm) # apply it to every gauge.

  start_end_matrix |>
    as_tibble() |>
    mutate(length = end_index - start_index + 1) |>
    filter(length > min_run_length - 1) # + 1 to make the input more easy to understand
}


list_start_end_index <- map(
  .x = unique(yearly_data$gauge),
  .f = gauge_continous_start_end,
  data = yearly_data,
  min_run_length = min_run_length 
) 

names(list_start_end_index) <- paste0("gauge_", unique(yearly_data$gauge))

start_end_index <- do.call("rbind", list_start_end_index)

start_end_index <- start_end_index |>
  rownames_to_column(var = "gauge") |>
  as_tibble() |>
  mutate(gauge = str_remove(gauge, "\\.[1-9]")) |>
  mutate(gauge = str_remove(gauge, "gauge_")) # remove anything with "." followed by a number


## Add a column to show what data will be used in calibration ==================
### I need to be more strict with drought catchment -
### There is a drought when the drought column is TRUE between the start_stop_indexes.
### If there is a drought but not between the start_stop_indexes the model will
### calibrate as if there is no drought.

rle_q_mm <- rle(!is.na(yearly_data$q_mm)) # finds NA values. Include = FALSE
rle_q_mm$values[rle_q_mm$lengths < min_run_length] <- FALSE # If not na (i.e., a number) but the run is < 2 set to Include = FALSE
include <- inverse.rle(rle_q_mm)

yearly_data <- yearly_data |>
  mutate(included_in_calibration = include)


counting_chunks <- function(input_vector) {
  nrow(continuous_run_start_end(input_vector))
}


## Gauge specific data (with catchment information) ============================
gauge_data <- yearly_data |>
  summarise(
    drought = any(drought == TRUE & included_in_calibration == TRUE), # if there is a drought but we are not calibrating on it we ignore that there is a drought
    max_run_start = max_continuous_run_start_end(q_mm)[1],
    max_run_end = max_continuous_run_start_end(q_mm)[2],
    record_length = sum(included_in_calibration),
    chunks = counting_chunks(q_mm),
    .by = gauge
  ) |>
  mutate(
    max_run = (max_run_end - max_run_start) + 1
  ) |>
  filter(record_length >= min_record_length_years) |>
  left_join( # Add chunks to gauge_data i.e., continuous runs
    catchment_information, 
    by = join_by(gauge)
  ) |> 
  # remove chunk and run length stuff
  select(!c(max_run_start, max_run_end, chunks, max_run))



# Add climate type to gauge data -----------------------------------------------
# Assign gauge a climate type using kgc 

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
  ) |> 
  select(gauge, lat, lon, climate_type, major_climate_type)


## Add to gauge data ===========================================================
gauge_data <- gauge_data |> 
  left_join(
    climate_type_gauge_data,
    by = join_by(gauge, lat, lon)
  )


# gauge_data filtered out gauges that do not meet record length requirements ---
# remove these gauges from start_stop and yearly data
yearly_data <- yearly_data |> 
  filter(gauge %in% gauge_data$gauge)

start_end_index <- start_end_index |> 
  filter(gauge %in% gauge_data$gauge)


# Account for pre-industrial-CO2 -----------------------------------------------
with_NA_yearly_data <- yearly_data |>
  mutate(CO2 = CO2 - pre_ind_CO2_ppm) |> 
  # I am currently not using standardised_warm_to_cool_ratio - remove
  select(!standardised_warm_to_cool_season_rainfall_ratio)




# Save .csv  -------------------------------------------------------------------
write_csv(gauge_data, paste0("./Data/Tidy/gauge_information_CAMELS.csv"))
write_csv(start_end_index, paste0("./Data/Tidy/start_end_index.csv"))
write_csv(with_NA_yearly_data, paste0("./Data/Tidy/with_NA_yearly_data_CAMELS.csv"))






