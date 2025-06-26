# Time of emergence analysis

# Figures produced in this R file ----------------------------------------------

# 1. Main --> ToE_map_aus_uncertainty.pdf (not working)
# 2. Supplementary --> ToE_vs_uncertainty_and_evidence_ratio.pdf (get working - called ToE_vs_uncertainty in old files)
# 3. Supplementary --> ToE_vs_record_length.pdf (get working)
# 4. Testing --> ToE_cdf.pdf



# DREAM uncertainty using log-sinh incomplete




# CODE









# Figure 3. Time of emergence map ----------------------------------------------


# Time of emergence calculation ------------------------------------------------
best_model_per_gauge <- best_CO2_non_CO2_per_gauge |>
  slice_min(
    AIC,
    by = gauge
  ) |>
  distinct()


## Turn the a5 parameter (CO2 pmm) into a given year ===========================
CO2_time_of_emergence <- best_model_per_gauge |>
  filter(parameter == "a5")

CO2_data <- readr::read_csv(
  "./Data/Raw/20241125_Mauna_Loa_CO2_data.csv",
  skip = 43,
  col_select = !unc,
  show_col_types = FALSE
) |>
  mutate(
    CO2_280 = mean - 280
  )

CO2 <- CO2_data |> pull(CO2_280)
year <- CO2_data |> pull(year)


## If CO2 for a given year minus the a5 parameter is negative then, the a5
## parameter has not kicked in. Take the first positive value and get the
## year using the CO2 index
single_a5_to_time_of_emergence <- function(a5, CO2, year) {
  adjusted_CO2 <- if_else(CO2 - a5 < 0, 0, CO2 - a5)
  
  year_where_CO2_impacts_flow <- year[adjusted_CO2 != 0][1]
  
  return(year_where_CO2_impacts_flow)
}

a5_to_time_of_emergence <- function(CO2, year) {
  # Save CO2 and year vectors to the function itself
  # They do not change
  
  force(CO2)
  force(year)
  stopifnot(is.numeric(CO2))
  stopifnot(is.numeric(year))
  
  function(a5) {
    future_map_dbl(
      .x = a5,
      .f = single_a5_to_time_of_emergence,
      CO2 = CO2,
      year = year,
      .progress = TRUE
    )
  }
}

# Function factory
adjusted_a5_to_ToE <- a5_to_time_of_emergence(CO2 = CO2, year = year)

# Add state
gauge_state <- gauge_information |>
  select(gauge, state)

time_of_emergence_data <- CO2_time_of_emergence |>
  mutate(
    year_time_of_emergence = adjusted_a5_to_ToE(parameter_value)
  ) |>
  left_join(
    gauge_state,
    by = join_by(gauge)
  )


# Does the time of emergence kick in at the last couple of years 
best_gauges_time_of_emergence <- time_of_emergence_data |> 
  pull(gauge) |> 
  unique()

final_observed_year <- streamflow_results |> 
  filter(gauge %in% best_gauges_time_of_emergence) |> 
  select(gauge, year) |> 
  distinct() |> 
  slice_max(year)
# All the years end during 2021

check_late_ToE <- time_of_emergence_data |> 
  filter(year_time_of_emergence >= 2021)

# 21 out of the 253 catchments where CO2 is the best model have a CO2
# ToE in the last year



## Interquartile range from DREAM ==============================================

### Summarising parameter uncertainty ##########################################
# I must do this in two steps because my PC does not have enough RAM
# Files are too big to combine - combine in R


# Comment out when complete and use the summarised tibble

# Method:
# 1. read the sequence .csv
# 2. REMOVE --> remove the duplicate gauges (maybe use chain and n to identify duplicates?)
# 3. put load and filter by a5
# 4. run plan() --> adjusted_a5_to_ToE code to turn CO2 ppm into year
# 5. summarise to find IQR
# 6. Save results
# 7. Repeat 1-6 for the number of parts. Cannot do them all at once due to
#    RAM limiations
# 8. Combine repeated files into one --> summarised_sequences_ToE_20250429.csv

# Get all sequences
# all_sequence_files <- list.files(
#  path = "Results/my_dream/",
#  pattern = "big_chunk",
#  full.names = TRUE,
#  recursive = FALSE
# )

# parameter_uncertainty <- all_sequence_files[1] |> # we cannot do everything in one go > 32 gb of RAM
#  read_csv(
#    show_col_types = FALSE
#  ) |>
#  filter(parameter == "a5")


# - Find IQR of ToE
# - I want this in years not CO2 pmm
# - To do this I must convert the a5 values from p_mm to years
# - CO2 and year is constant for all catchments - function factory
#   to add the CO2 and year from data, then only input a5

# This works but it is really slow and RAM intensive

# plan(multisession, workers = length(availableWorkers()))
# filtered_parameter_uncertainty <- parameter_uncertainty |>
#  mutate(
#    year_ToE = adjusted_a5_to_ToE(parameter_value)
#  )


# ToE_range_uncertainty <- filtered_parameter_uncertainty |>
#  summarise(
#    DREAM_ToE_IQR = IQR(year_ToE),
#    DREAM_ToE_median = median(year_ToE),
#    .by = gauge
#  )


### Save ToE_range_uncertainty ###
### Repeat with part_1 and part_2
### Join and save the summarised information for plotting
# write_csv(
#  ToE_range_uncertainty,
#  file = "Results/my_dream/part_5_summarised_sequences_20250429.csv"
# )

# combined_summarised_ToE_data <- list.files( # get
#  path = "./Results/my_dream/",
#  recursive = FALSE, # I don't want it looking in other folders
#  pattern = "summarised_sequence",
#  full.names = TRUE
# )

# combined_summarised_ToE_data |> # merge and save
#  read_csv(show_col_types = FALSE) |>
#  write_csv(
#    paste0("./Results/my_dream/summarised_sequences_ToE_", get_date(), ".csv")
#  )

# duplicated_gauges <- test |> pull(gauge)
# duplicated_gauges[duplicated(duplicated_gauges)]


# There are duplicates - that is probably why the code broke
# The gauges are "612034" "406235" "406250" "407213"
# They all produce the same uncertainty. For now
# use distinct. This temporarily solves the problem

summarised_sequences_ToE <- read_csv(
  "Results/my_dream/summarised_sequences_ToE_20250429.csv",
  show_col_types = FALSE
) |>
  distinct()









## Create my own bins ==========================================================
min_CO2 <- data |>
  pull(CO2) |>
  min()

lat_lon_gauge <- gauge_information |>
  select(gauge, lat, lon)

custom_bins_time_of_emergence_data <- time_of_emergence_data |>
  # add uncertainty
  left_join(
    summarised_sequences_ToE,
    by = join_by(gauge)
  ) |>
  # add state
  left_join(
    state_gauge,
    by = join_by(gauge)
  ) |>
  mutate(custom_bins = year_time_of_emergence - (year_time_of_emergence %% 10)) |>
  mutate(custom_bins = as.character(custom_bins)) |>
  mutate(
    custom_bins = if_else(parameter_value < min_CO2, "Before 1959", custom_bins)
  ) |>
  # Add evidence ratio for potential filtering
  left_join(
    a3_direction_binned_lat_lon_evidence_ratio,
    by = join_by(gauge)
  ) |>
  # Only include moderately strong and above evidence ratio
  filter(!binned_evidence_ratio %in% c("Weak", "Moderate", "Moderately Strong")) |>
  # clean it up and add lat lon
  select(gauge, state, year_time_of_emergence, custom_bins, DREAM_ToE_IQR, binned_evidence_ratio) |>
  left_join(
    lat_lon_gauge,
    by = join_by(gauge)
  ) |>
  # factor custom_bins to force then into order
  mutate(
    custom_bins = factor(custom_bins, levels = c("Before 1959", "1960", "1970", "1980", "1990", "2000", "2010", "2020"))
  ) |>
  # Binning is required to DREAM_ToE_IQR
  mutate(
    binned_DREAM_ToE_IQR = case_when(
      DREAM_ToE_IQR <= 10 ~ "<=10",
      between(DREAM_ToE_IQR, 11, 20) ~ "11-20",
      between(DREAM_ToE_IQR, 21, 30) ~ "21-30",
      between(DREAM_ToE_IQR, 31, 40) ~ "31-40",
      between(DREAM_ToE_IQR, 41, 50) ~ "41-50",
      DREAM_ToE_IQR > 50 ~ ">50",
      .default = NA
    )
  ) |>
  # factor it so all plots ahve the same bins
  mutate(
    binned_DREAM_ToE_IQR = factor(binned_DREAM_ToE_IQR, levels = rev(c("<=10", "11-20", "21-30", "31-40", "41-50", "51-60", ">60")))
  ) |>
  # Points are plotted based on the order they appear in the tibble
  # Order the tibble from largest IQR to smalleest before plotting.
  arrange(desc(DREAM_ToE_IQR))






# Examine the relationship between climate type and time of emergence ----------
climate_type_ToE <- custom_bins_time_of_emergence_data |>
  left_join(
    gauge_information_with_climate,
    by = join_by(gauge, lat, lon)
  )



climate_type_aus_map <- ozmaps::ozmap("country") |>
  uncount(7) |> # repeat the geometry 4 times
  mutate(
    custom_bins = c("Before 1959", "1960", "1970", "1980", "1990", "2000", "2010")
  ) |>
  mutate(
    custom_bins = factor(custom_bins, levels = c("Before 1959", "1960", "1970", "1980", "1990", "2000", "2010"))
  )


scale_size_limits <- custom_bins_time_of_emergence_data |>
  pull(DREAM_ToE_IQR) |>
  range() # can round up if I want to

climate_type_ToE_plot <- climate_type_aus_map |>
  ggplot(aes(geometry = geometry)) +
  geom_sf(
    colour = "black",
    fill = "grey80"
  ) +
  coord_sf(xlim = c(111, 155), ylim = c(-44.5, -9.5)) +
  geom_point(
    aes(x = lon, y = lat, fill = climate_type, size = DREAM_ToE_IQR), # circle size as uncertainty?
    data = climate_type_ToE,
    inherit.aes = FALSE,
    colour = "black",
    stroke = 0.1,
    shape = 21
  ) +
  scale_size_binned(limits = scale_size_limits, breaks = c(5, 10, 15, 20, 25, 30, 35), range = c(0.1, 3)) +
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Climate type",
    size = "Time of Emergence Uncertaity (Years)"
  ) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.7, 0.15)
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 5, shape = 21), nrow = 3), # Wrap legend with nrow
    size = guide_bins(show.limits = TRUE, direction = "horizontal")
  ) +
  facet_wrap(~custom_bins)



ggsave(
  filename = "climate_type_ToE.pdf",
  plot = climate_type_ToE_plot,
  path = "Graphs/Supplementary_Figures",
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)