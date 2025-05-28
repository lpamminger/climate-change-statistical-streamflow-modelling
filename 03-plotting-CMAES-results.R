# Script produces the four key figures of the second paper


# TODO - transfer code/graphs from playing_with_graphs
# Key graphs to transfer are:
# 1. Time of emergence histogram
# - account for the different types of CO2 models - slope vs. intercept
# 2. graphically compare (colour gradient evidence ratio, shape pos/neg, outline
#    or stroke colour can be slope or intercept)
# 3. some sort of slope and intercept analysis. Like proportion of pos/neg slope,
#    pos/neg intercept



cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, ozmaps, sf, ggmagnify, furrr, parallel) 


## Utility functions ===========================================================
source("./Functions/utility.R")


## Import streamflow functions =================================================
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
source("./Functions/boxcox_transforms.R")


# Import data ------------------------------------------------------------------
CMAES_results <- read_csv(
  "./Results/my_cmaes/CMAES_parameter_results_20250331.csv", 
  show_col_types = FALSE
) 

data <- readr::read_csv(
  "./Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  col_types = "icdddddddll", # ensuring each column is a the correct type
  show_col_types = FALSE
)

streamflow_results <- read_csv(
  "Results/my_cmaes/CMAES_streamflow_results_20250331.csv",
  show_col_types = FALSE,
  col_select = !optimiser
)

gauge_information <- read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE,
  col_select = c(gauge, bc_lambda, state, lat, lon)
) # CAMELS is Australia wide.

best_CO2_non_CO2_per_gauge <- read_csv(
  "./Results/my_cmaes/unmodified_best_CO2_non_CO2_per_catchment_CMAES_20250331.csv",
  show_col_types = FALSE
) 
  

# Uncertainty important and calculations done in section "Figure 3. ToE map" 
  

# Objects used for all figures -------------------------------------------------
### Add lat and lon ############################################################
lat_lon_gauge <- gauge_information |> 
  select(gauge, lat, lon)

# Get shapefiles for Australia ------------------------------------------------
aus_map <- ozmaps::ozmap(x = "states") |>
  filter(!NAME %in% c("Other Territories")) |>
  rename(state = NAME) |>
  mutate(
    state = case_when(
      state == "New South Wales" ~ "NSW",
      state == "Victoria" ~ "VIC",
      state == "Queensland" ~ "QLD",
      state == "South Australia" ~ "SA",
      state == "Western Australia" ~ "WA",
      state == "Tasmania" ~ "TAS",
      state == "Northern Territory" ~ "NT",
      state == "Australian Capital Territory" ~ "ACT",
    )
  )

combine_NSW_ACT <- aus_map |>
  filter(state %in% c("NSW", "ACT")) |>
  st_union()

aus_map[1, 2] <- list(combine_NSW_ACT)

aus_map <- aus_map |>
  filter(state != "ACT")





# Now Supplementary - Figure S1. A of best components of streamflow models of Australia -------------
single_aus_map <- ozmaps::ozmap("country") |> 
  uncount(5) |>  # repeat the geometry 4 times
  mutate(
    simple_name = c("Drought", "Autocorrelation", "CO2 Intercept", "CO2 Slope", "Rainfall Seasonality")
  )

best_model_per_gauge <- best_CO2_non_CO2_per_gauge |>
  slice_min(
    AIC,
    by = gauge
  ) |> 
  distinct()

lat_long_gauge <- gauge_information |> 
  select(gauge, lat, lon)


# Organise data
# - facet on component
# - legend/colour/fill on included or not for gauge
# structure |gauge|component_name|show|

# 123456|Autocorrelation|FALSE
# 123456|Drought|TRUE
# etc.

# make my own tibble
parameters <- c("Autocorrelation", "CO2 Intercept", "CO2 Slope", "Drought", "Rainfall Seasonality") 
gauges <- best_model_per_gauge |> pull(gauge) |> unique() 

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
    lat_long_gauge,
    by = join_by(gauge)
  )


# Count of components and plot labels
total_gauges <- best_model_per_gauge |> pull(gauge) |> unique() |> length()
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

# Plot map
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

#model_components

ggsave(
  filename = "./Graphs/Supplementary_Figures/model_components_v2.pdf",
  plot = model_components,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)


# Figure 2. A map of best CO2 vs best non-CO2 of Australia ---------------------

## Calculate evidence ratio ====================================================
evidence_ratio_calc <- best_CO2_non_CO2_per_gauge |>
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
  arrange(evidence_ratio)


## Tidy evidence ratio data for plotting =======================================

lat_long_evidence_ratio <- evidence_ratio_calc |>
  select(!c(AIC_difference)) |>
  left_join(
    lat_lon_gauge,
    by = join_by(gauge)
  )

### Add qualitative labels instead of using numerical evidence ratio ###########
state_gauge <- gauge_information |> 
  select(gauge, state)

binned_lat_lon_evidence_ratio <- lat_long_evidence_ratio |>
  mutate(
    binned_evidence_ratio = case_when(
      between(evidence_ratio, -1E1, 1E1) ~ "Weak",
      between(evidence_ratio, 1E1, 1E2) ~ "Moderate",
      between(evidence_ratio, 1E2, 1E3) ~ "Moderately Strong",
      between(evidence_ratio, 1E3, 1E4) ~ "Strong",
      between(evidence_ratio, 1E4, 1E6) ~ "Very Strong",
      between(evidence_ratio, 1E6, Inf) ~ "Extremely Strong",
      .default = NA
    )
  ) |>
  left_join(
    state_gauge,
    by = join_by(gauge)
  ) |>
  mutate(
    binned_evidence_ratio = factor(
      binned_evidence_ratio,
      levels = c("Weak", "Moderate", "Moderately Strong", "Strong", "Very Strong", "Extremely Strong")
    )
  )



### Add direction of change and whether the slope/intercept changed ############
best_model_per_gauge <- best_CO2_non_CO2_per_gauge |>
  slice_min(
    AIC,
    by = gauge
  ) |> 
  distinct()


direction_of_a3_change <- best_model_per_gauge |> 
  filter(parameter %in% c("a3_intercept", "a3_slope")) |>
  mutate(
    intercept_or_slope = if_else(str_detect(streamflow_model, "intercept"), "Intercept", "Slope")
  ) |> 
  select(gauge, streamflow_model, parameter, parameter_value, intercept_or_slope) |> 
  mutate(
    CO2_direction = if_else(parameter_value < 0, "Negative", "Positive")
  ) |> 
  select(gauge, CO2_direction, intercept_or_slope)



a3_direction_binned_lat_lon_evidence_ratio <- binned_lat_lon_evidence_ratio |>
  left_join(
    direction_of_a3_change,
    by = join_by(gauge)
  ) |>
  replace_na(list(CO2_direction = "No CO2 Term", intercept_or_slope = "No CO2 Term")) |>
  unite(
    col = "impact_of_CO2_term",
    CO2_direction,
    intercept_or_slope,
    sep = "-"
  ) |>
  mutate(
    impact_of_CO2_term = if_else(impact_of_CO2_term == "No CO2 Term-No CO2 Term", "No CO2 Term", impact_of_CO2_term)
  ) |>
  mutate(
    impact_of_CO2_term = factor(
      impact_of_CO2_term,
      levels = c("No CO2 Term", "Negative-Intercept", "Positive-Intercept", "Negative-Slope", "Positive-Slope")
    )
  ) 




# Make final plot --------------------------------------------------------------

### Custom colour palette 
custom_palette <- function(x) {
  rev(c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#f7f7f7"))
}


## Generate Insets =============================================================
### Filter data by state #######################################################

QLD_data <- a3_direction_binned_lat_lon_evidence_ratio |>
  filter(state == "QLD")

NSW_data <- a3_direction_binned_lat_lon_evidence_ratio |>
  filter(state == "NSW")

VIC_data <- a3_direction_binned_lat_lon_evidence_ratio |>
  filter(state == "VIC")

WA_data <- a3_direction_binned_lat_lon_evidence_ratio |>
  filter(state == "WA")

TAS_data <- a3_direction_binned_lat_lon_evidence_ratio |>
  filter(state == "TAS")


### Generate inset plots #######################################################

inset_plot_QLD <- aus_map |>
  filter(state == "QLD") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = QLD_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
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
  theme_void()


inset_plot_NSW <- aus_map |>
  filter(state == "NSW") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = NSW_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
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
  theme_void()



inset_plot_VIC <- aus_map |>
  filter(state == "VIC") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = VIC_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
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
  theme_void()



inset_plot_WA <- aus_map |>
  filter(state == "WA") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = WA_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
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
  theme_void()



inset_plot_TAS <- aus_map |>
  filter(state == "TAS") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = TAS_data,
    aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
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
  theme_void()



## Put it together =============================================================

single_map_aus <- aus_map |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = a3_direction_binned_lat_lon_evidence_ratio,
    mapping = aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term),
    size = 3,
    colour = "black",
    stroke = 0.1
  ) +
  theme_bw() +
  scale_fill_manual(
    values = custom_palette(),
    drop = FALSE
  ) +
  scale_shape_manual(
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
    #aes(from = state == "VIC"), # use aes rather than manually selecting area
    from = c(141, 149.5, -39, -34),
    to = c(95, 136, -38, -60),
    shadow = FALSE,
    plot = inset_plot_VIC,
    proj = "single"
  ) +
  # magnify QLD
  geom_magnify(
    from = c(145, 155, -29.2, -16),
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
    x = NULL,#"Latitude",
    y = NULL,#"Longitude",
    fill = "Evidence Ratio",
    shape = bquote("Impact of "~CO[2]~"Term")
  ) +
  theme(
    legend.key = element_rect(fill = "grey80"),
    legend.title = element_text(hjust = 0.5),
    legend.background = element_rect(colour = "black"),
    axis.text = element_blank(), 
    legend.position = "inside",
    legend.position.inside = c(0.351, 0.9),
    legend.box = "horizontal", # side-by-side legends
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank()
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 5, shape = 21), nrow = 3), # Wrap legend with nrow
    shape = guide_legend(override.aes = list(size = 5, fill = "grey50"), nrow = 3)
  )



ggsave(
  filename = "./Graphs/Figures/evidence_ratio_aus_with_zoom_v4.pdf",
  plot = single_map_aus,
  device = "pdf",
  width = 232,
  height = 200,#210,
  units = "mm"
)



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
#all_sequence_files <- list.files(
#  path = "Results/my_dream/",
#  pattern = "big_chunk",
#  full.names = TRUE,
#  recursive = FALSE
#) 

#parameter_uncertainty <- all_sequence_files[1] |> # we cannot do everything in one go > 32 gb of RAM
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

#plan(multisession, workers = length(availableWorkers()))
#filtered_parameter_uncertainty <- parameter_uncertainty |> 
#  mutate(
#    year_ToE = adjusted_a5_to_ToE(parameter_value)
#  )


#ToE_range_uncertainty <- filtered_parameter_uncertainty |> 
#  summarise(
#    DREAM_ToE_IQR = IQR(year_ToE),
#    DREAM_ToE_median = median(year_ToE),
#    .by = gauge
#  )


### Save ToE_range_uncertainty ###
### Repeat with part_1 and part_2 
### Join and save the summarised information for plotting
#write_csv(
#  ToE_range_uncertainty,
#  file = "Results/my_dream/part_5_summarised_sequences_20250429.csv"
#)

#combined_summarised_ToE_data <- list.files( # get
#  path = "./Results/my_dream/",
#  recursive = FALSE, # I don't want it looking in other folders
#  pattern = "summarised_sequence",
#  full.names = TRUE
#)

#combined_summarised_ToE_data |> # merge and save
#  read_csv(show_col_types = FALSE) |> 
#  write_csv(
#    paste0("./Results/my_dream/summarised_sequences_ToE_", get_date(), ".csv")
#  )

#duplicated_gauges <- test |> pull(gauge) 
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


## Time of emergence plotting ==================================================

## Generate Insets =============================================================
### Filter data by state #######################################################

QLD_data <- custom_bins_time_of_emergence_data |>
  filter(state == "QLD")

NSW_data <- custom_bins_time_of_emergence_data |>
  filter(state == "NSW")

VIC_data <- custom_bins_time_of_emergence_data |>
  filter(state == "VIC")

WA_data <- custom_bins_time_of_emergence_data |>
  filter(state == "WA")

TAS_data <- custom_bins_time_of_emergence_data |>
  filter(state == "TAS")


### Generate inset plots #######################################################


# Generate inset histgrams of time of emergence for each plot
generate_ToE_histograms <- function(STATE_data) {
  STATE_data |> 
    summarise(
      count = n(),
      .by = binned_evidence_ratio
    ) |> 
    # I want to show all the possible outcomes in the histogram
    # This currently does not work
    mutate(
      binned_evidence_ratio = factor(binned_evidence_ratio, levels = c("Weak", "Moderate", "Moderately Strong", "Strong", "Very Strong", "Extremely Strong"))
    ) |> 
    ggplot(aes(x = binned_evidence_ratio, y = count)) +
    geom_col() +
    scale_x_discrete(drop = FALSE) +
    labs(
      x = "Frequency",
      y = "Evidence Ratio"
    ) +
    theme_bw()
}

state_ToE_histograms <- map(
  .x = list(QLD_data, NSW_data, VIC_data, WA_data, TAS_data),
  .f = generate_ToE_histograms
)

state_ToE_histograms
# Need to add these histograms as grobs in the inset
# Main plot -> inset map -> inset histogram
# Does patchwork do inset histograms?
# See documentation: https://patchwork.data-imaginist.com/reference/inset_element.html

# scale_size_binned()
# Using scale_size_binned() is technically between because it is a 
# does the binning for me.

#
scale_size_limits <- custom_bins_time_of_emergence_data |> 
  pull(DREAM_ToE_IQR) |> 
  range() # can round up if I want to 


transparent_dots_constant <- 0.75

inset_plot_QLD <- aus_map |>
  filter(state == "QLD") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = QLD_data,
    aes(x = lon, y = lat, fill = custom_bins, size = DREAM_ToE_IQR),
    show.legend = FALSE,
    stroke = 0.1,
    alpha = transparent_dots_constant,
    shape = 21,
    colour = "black"
  ) +
  scale_fill_brewer(palette = "BrBG", drop = FALSE) +
  scale_size_binned(limits = scale_size_limits) + # range = c(0, 2) dictates the size of the dots (important)
  guides(size = guide_bins(show.limits = TRUE)) +
  theme_void()




inset_plot_NSW <- aus_map |>
  filter(state == "NSW") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = NSW_data,
    aes(x = lon, y = lat, fill = custom_bins, size = DREAM_ToE_IQR),
    show.legend = FALSE,
    stroke = 0.1,
    alpha = transparent_dots_constant,
    shape = 21,
    colour = "black"
  ) +
  scale_fill_brewer(palette = "BrBG", drop = FALSE) +
  scale_size_binned(limits = scale_size_limits) + # range = c(0, 2) dictates the size of the dots (important)
  guides(size = guide_bins(show.limits = TRUE)) +
  theme_void()



inset_plot_VIC <- aus_map |>
  filter(state == "VIC") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = VIC_data,
    aes(x = lon, y = lat, fill = custom_bins, size = DREAM_ToE_IQR),
    show.legend = FALSE,
    alpha = transparent_dots_constant,
    stroke = 0.1,
    shape = 21,
    colour = "black"
  ) +
  scale_fill_brewer(palette = "BrBG", drop = FALSE) +
  scale_size_binned(limits = scale_size_limits) + # range = c(0, 2) dictates the size of the dots (important)
  guides(size = guide_bins(show.limits = TRUE)) +
  theme_void()



inset_plot_WA <- aus_map |>
  filter(state == "WA") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = WA_data,
    aes(x = lon, y = lat, fill = custom_bins, size = DREAM_ToE_IQR),
    show.legend = FALSE,
    stroke = 0.1,
    alpha = transparent_dots_constant,
    shape = 21,
    colour = "black"
  ) +
  scale_fill_brewer(palette = "BrBG", drop = FALSE) +
  scale_size_binned(limits = scale_size_limits) + # range = c(0, 2) dictates the size of the dots (important)
  guides(size = guide_bins(show.limits = TRUE)) +
  theme_void()



inset_plot_TAS <- aus_map |>
  filter(state == "TAS") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = TAS_data,
    aes(x = lon, y = lat, fill = custom_bins, size = DREAM_ToE_IQR),
    show.legend = FALSE,
    stroke = 0.1,
    alpha = transparent_dots_constant,
    shape = 21,
    colour = "black",
  ) +
  scale_fill_brewer(palette = "BrBG", drop = FALSE) +
  scale_size_binned(limits = scale_size_limits) + # range = c(0, 2) dictates the size of the dots (important)
  guides(size = guide_bins(show.limits = TRUE)) +
  theme_void()



## Put it together =============================================================
## This probably could be functioned...
ToE_map_aus <- aus_map |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = custom_bins_time_of_emergence_data,
    mapping = aes(x = lon, y = lat, fill = custom_bins, size = DREAM_ToE_IQR),
    stroke = 0.1,
    shape = 21,
    colour = "black"
  ) +
  scale_fill_brewer(palette = "BrBG", drop = FALSE) +
  scale_size_binned(limits = scale_size_limits) + # range = c(0, 2) dictates the size of the dots (important)
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
    #aes(from = state == "VIC"), # use aes rather than manually selecting area
    from = c(141, 149.5, -39, -34),
    to = c(95, 136, -38, -60),
    shadow = FALSE,
    plot = inset_plot_VIC,
    proj = "single"
  ) +
  # magnify QLD
  geom_magnify(
    from = c(145, 155, -30, -15),
    to = c(157, 178, -29.5, 1.5),
    shadow = FALSE,
    expand = 0,
    plot = inset_plot_QLD,
    proj = "single"
  ) +
  # magnify NSW
  geom_magnify(
    from = c(146.5, 154, -38, -28),
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
    x = NULL,#"Latitude",
    y = NULL,#"Longitude",
    fill = "Time of Emergence",
    size = "Time of Emergence Uncertainty Years (IQR)"
  ) +
  theme(
    legend.key = element_rect(fill = "grey80"),
    legend.title = element_text(hjust = 0.5),
    legend.title.position = "top",
    legend.background = element_rect(colour = "black"),
    axis.text = element_blank(), 
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.325, 0.9),
    legend.box = "horizontal"#, # side-by-side legends
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 5, shape = 21), nrow = 3), # Wrap legend with nrow
    size = guide_bins(show.limits = TRUE, direction = "horizontal")
  )


#ToE_map_aus

ggsave(
  filename = "./Graphs/Figures/ToE_map_aus_uncertainty_v7_moderate_strong_above.pdf",
  plot = ToE_map_aus,
  device = "pdf",
  width = 232,
  height = 200,
  units = "mm"
)

## Temp histogram
x <- custom_bins_time_of_emergence_data |> 
  summarise(
    n = n(),
    .by = c(custom_bins, state)
  ) |>
  arrange(state) 
  ggplot(aes(x = custom_bins, y = n)) +
  geom_col() +
  theme_bw()

## TEMP What does ToE vs. DREAM_IQR_ToE look like
ToE_against_uncertainty <- custom_bins_time_of_emergence_data |> 
  ggplot(aes(x = year_time_of_emergence, y = DREAM_ToE_IQR)) +
  geom_point() +
  #geom_smooth(method = lm, formula = y ~ x) +
  labs(
    x = "Time of Emergence",
    y = "Time of Emergence Uncertainty Years (IQR)"
  ) +
  facet_wrap(~binned_evidence_ratio, scales = "fixed") +
  theme_bw()


# TEMP. Compare evidence ratio with ToE
# evidence ratio data is lat_long_evidence_ratio
# ToE data is custom_bins_time_of_emergence_data
compare_ToE_and_evi_ratio <- custom_bins_time_of_emergence_data |> 
  left_join(
    lat_long_evidence_ratio,
    by = join_by(gauge, lat, lon)
  )

ToE_against_evi_ratio <- compare_ToE_and_evi_ratio |> 
  ggplot(aes(x = year_time_of_emergence, y = evidence_ratio)) +
  geom_point() +
  scale_y_log10() +
  labs(
    x = "Time of Emergence",
    y = "Evidence Ratio (Log Scale)"
  ) +
  theme_bw()



ggsave(
  filename = "ToE_vs_uncertainty.pdf",
  plot = gridExtra::marrangeGrob(
    list(ToE_against_uncertainty, ToE_against_evi_ratio), # add other graphs here:  
    nrow = 1, 
    ncol = 1,
    top = NULL # no page numbers
  ),
  device = "pdf",
  path = "./Graphs/Supplementary_Figures",
  width = 297,
  height = 210,
  units = "mm"
)




# Figure 4. How has CO2 impacted streamflow map --------------------------------
stop_here <- tactical_typo()



# Other figure - Streamflow plots of best-CO2, best-non-CO2 and observed -------


## Filter based on best non-CO2 and CO2 models =================================
best_streamflow_results <- streamflow_results |>
  semi_join(
    best_CO2_non_CO2_per_gauge,
    by = join_by(gauge, streamflow_model)
  )


## Only include best streamflow that was calibrated on =========================
## join the included_in_calibration column
in_calibration <- data |> 
  select(year, gauge, included_in_calibration)

best_calibration_streamflow_results <- best_streamflow_results |> 
  left_join(
    in_calibration,
    by = join_by(year, gauge)
  ) |> 
  filter(included_in_calibration)



## Summarise results into a tidy format ========================================
tidy_boxcox_streamflow <- best_calibration_streamflow_results |>
  drop_na() |>  # only include if observed streamflow is present
  pivot_longer(
    cols = c(observed_boxcox_streamflow, modelled_boxcox_streamflow),
    names_to = "name",
    values_to = "boxcox_streamflow"
  ) |>
  mutate(
    name = if_else(name == "observed_boxcox_streamflow", "observed", "modelled")
  ) |> 
  mutate(
    contains_CO2 = str_detect(streamflow_model, "CO2")
  ) |> 
  mutate(
    legend = case_when(
      contains_CO2 & (name != "observed") ~ "modelled_CO2",
      !contains_CO2 & (name != "observed") ~ "modelled_non_CO2",
      .default = name
    )
  ) |> 
  select(
    !c(streamflow_model, objective_function, included_in_calibration, name, contains_CO2)
    )


## Convert from box-cox space to real space ====================================
### bc lambda found in gauge_information
# boxcox_inverse_transform()
tidy_streamflow <- tidy_boxcox_streamflow |>
  left_join(
    gauge_information,
    by = join_by(gauge)
  ) |>
  select(!c(state, lat, lon)) |>
  mutate(
    streamflow = boxcox_inverse_transform(yt = boxcox_streamflow, lambda = bc_lambda, lambda_2 = 1),
    .by = gauge
  )






## Plot results ================================================================
### having 534 graphs on a single page is too much = slows pc.
### Spread across multiple

# rep by gauge

# Split the tibble into X groups
# Gauges must not be across multiple groups
# Randomly assign 1,2 or 3 to each group then split? This works
# but is probably not the best way of doing it
ready_split_tidy_streamflow <- tidy_streamflow |> 
  mutate(
    split = sample(c(1, 2, 3), size = 1, replace = TRUE),
    .by = gauge,
    .before = 1
  ) 

split_tidy_streamflow <- ready_split_tidy_streamflow |> 
  split(f = ready_split_tidy_streamflow$split)


## Plot rainfall-runoff relationships ==========================================
## Need to split plot into smaller chunks. PC struggles to load it.

repeat_rainfall_runoff_plots <- function(segmented_tidy_streamflow) {
  rainfall_runoff_plots <- segmented_tidy_streamflow |> 
    ggplot(
      aes(
        x = precipitation, 
        y = boxcox_streamflow, 
        colour = legend, 
        shape = legend
        )
      ) +
    geom_point() +
    geom_line(
      stat = "smooth",
      method = "lm",
      formula = y ~ x,
      alpha = 0.4,
      linewidth = 1
    ) +
    labs(
      x = "Annual Precipitation (mm)",
      y = "Annual Boxcox Streamflow (mm)"
    ) +
    scale_colour_brewer(palette = "Set1") +
    theme_bw() +
    facet_wrap(~gauge, scales = "free") +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom"
    )
}


rainfall_runoff_graphs <- map(
  .x = split_tidy_streamflow,
  .f = repeat_rainfall_runoff_plots
)

ggsave(
  filename = paste0("rainfall_runoff_comparison_", get_date(), ".pdf"),
  plot = gridExtra::marrangeGrob(rainfall_runoff_graphs, nrow = 1, ncol = 1),
  device = "pdf",
  path = "./Graphs/CMAES_graphs",
  width = 1189,
  height = 841,
  units = "mm"
)


## Plot streamflow time ========================================================
chunk_streamflow_timeseries_plot <- function(data) {
  
  data |>
    ggplot(aes(x = year, y = streamflow, colour = legend)) +
    geom_line(na.rm = TRUE, alpha = 0.5) +
    geom_point(na.rm = TRUE, size = 0.5, alpha = 0.5) +
    theme_bw() +
    scale_colour_brewer(palette = "Set1") +
    labs(
      x = "Year",
      y = "Streamflow (mm)"
    ) +
    facet_wrap(~gauge, scales = "free_y") +
    theme(legend.title = element_blank())
  
}


plot_streamflow_timeseries <- map(
  .x = split_tidy_streamflow,
  .f = chunk_streamflow_timeseries_plot
)



ggsave(
  filename = paste0("streamflow_timeseries_comparison_", get_date(), ".pdf"),
  plot = gridExtra::marrangeGrob(plot_streamflow_timeseries, nrow = 1, ncol = 1),
  device = "pdf",
  path = "./Graphs/CMAES_graphs",
  width = 1189,
  height = 841,
  units = "mm"
)


# TEMPORARY ---
# I am redoing the a3_slope
# I need to remove any a3_slope from the distribution pdf and trace pdfs
# In the pdfs search for `a3_slope`
# Get the page number
# Delete the pages

# Loading the trace plots takes a long time.
# The distribution plots have the same order. Use the order to sort the
# trace plots.

library(pdfsearch)
library(pdftools)
file <- list.files(path = "./Graphs/DREAM_graphs", full.names = TRUE)[3]

result <- keyword_search(
  file,
  keyword = c("a3_slope"),
  path = TRUE
) 

remove_pages <- result |> 
  pull(line_num)

keep_pages <- seq(from = 1, to = 534, by = 1)
keep_pages <- keep_pages[-remove_pages]


pdf_subset(input = file, pages = keep_pages, output = "./Graphs/DREAM_graphs/removed_a3_slope_trace_plots.pdf")
# Output can now be joined to a3_slope graphs