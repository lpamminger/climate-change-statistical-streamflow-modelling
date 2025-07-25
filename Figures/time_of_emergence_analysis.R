# Time of emergence analysis

# Figures produced in this R file ----------------------------------------------

# 1. Main --> ToE_map_aus_uncertainty.pdf (requires dream results - currently working without dream results)
# 2. Supplementary --> ToE_vs_uncertainty_and_evidence_ratio.pdf (requires dream results - called ToE_vs_uncertainty in old files)
# 3. Supplementary --> ToE_vs_record_length.pdf 
# 4. Supplementary --> ToE_vs_catchment_area.pdf
# 5. Testing --> ToE_cdf.pdf (not in file)
# 6. Testing --> something related to climate types (not in file)



# DREAM uncertainty using log-sinh incomplete




# CODE







# Import libraries -------------------------------------------------------------
pacman::p_load(tidyverse, ozmaps, sf, ggmagnify, furrr)



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
  mutate(year = as.integer(year))


gauge_information <- readr::read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE
)

lat_lon_gauge <- gauge_information |>
  select(gauge, lat, lon, state)


best_CO2_non_CO2_per_gauge <- read_csv(
  "./Modelling/Results/CMAES/best_CO2_non_CO2_per_catchment_CMAES.csv",
  show_col_types = FALSE
)

streamflow_results <- read_csv(
  "./Modelling/Results/CMAES/cmaes_streamflow_results.csv",
  show_col_types = FALSE
)



# Calculate time of emergence --------------------------------------------------
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

min_CO2 <- data |>
  pull(CO2) |>
  min()


time_of_emergence_data <- CO2_time_of_emergence |>
  mutate(
    year_time_of_emergence = adjusted_a5_to_ToE(parameter_value)
  ) |>
  left_join(
    lat_lon_gauge,
    by = join_by(gauge)
  ) |>
  # make decade time of emergence
  mutate(
    decade_time_of_emergence = year_time_of_emergence - year_time_of_emergence %% 10
  ) |>
  # turn it into a factor for plotting
  mutate(
    decade_time_of_emergence = as.character(decade_time_of_emergence),
    decade_time_of_emergence = if_else(
      decade_time_of_emergence == "1950",
      "Before 1960",
      decade_time_of_emergence
    ),
    decade_time_of_emergence = factor(
      decade_time_of_emergence,
      levels = c("Before 1960", as.character(seq(from = 1960, to = 2020, by = 10)))
    )
  )




# Calculate evidence ratio -----------------------------------------------------
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



## Add qualitative labels to evidence ratio ====================================


time_of_emergence_data <- time_of_emergence_data |>
  left_join(
    evidence_ratio_calc,
    by = join_by(gauge)
  ) |>
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
  # factor to force order of evidence ratio
  mutate(
    binned_evidence_ratio = factor(
      binned_evidence_ratio,
      levels = c("Weak", "Moderate", "Moderately Strong", "Strong", "Very Strong", "Extremely Strong")
    )
  ) |> 
  # For the purpose of analysis we are only interested in moderate strong or above. 
  filter(binned_evidence_ratio %in% c("Moderately Strong", "Strong", "Very Strong", "Extremely Strong")) 



# Does the time of emergence kick in at the last couple of years? -------------- 
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

cat(nrow(check_late_ToE), "catchments with a evidence ratio of moderately strong or greater have a ToE in the final year of data\n")







# Map of time of emergence------------------------------------------------------

## Get shapefiles for Australia ================================================
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


## Generate Insets =============================================================

### Filter data by state #######################################################

QLD_data <- time_of_emergence_data |>
  filter(state == "QLD")

NSW_data <- time_of_emergence_data |>
  filter(state == "NSW")

VIC_data <- time_of_emergence_data |>
  filter(state == "VIC")

WA_data <- time_of_emergence_data |>
  filter(state == "WA")

TAS_data <- time_of_emergence_data |>
  filter(state == "TAS")


### Generate inset plots #######################################################


# Need to add these histograms as grobs in the inset
# Main plot -> inset map -> inset histogram
# Does patchwork do inset histograms?
# See documentation: https://patchwork.data-imaginist.com/reference/inset_element.html

# scale_size_binned() - used for dream
# Using scale_size_binned() is technically between because it is a
# does the binning for me.


#scale_size_limits <- custom_bins_time_of_emergence_data |>
#  pull(DREAM_ToE_IQR) |>
#  range() # can round up if I want to


transparent_dots_constant <- 0.75
size <- 4 # remove when adding DREAM_ToE_IQR

inset_plot_QLD <- aus_map |>
  filter(state == "QLD") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = QLD_data,
    aes(x = lon, y = lat, fill = decade_time_of_emergence),#, size = DREAM_ToE_IQR),
    show.legend = FALSE,
    stroke = 0.1,
    alpha = transparent_dots_constant,
    shape = 21, # remove this with size
    colour = "black",
    size = size
  ) +
  scale_fill_brewer(palette = "BrBG", drop = FALSE) +
  #scale_size_binned(limits = scale_size_limits) + # range = c(0, 2) dictates the size of the dots (important)
  guides(size = guide_bins(show.limits = TRUE)) +
  theme_void() #+
# annotation_custom(
#  ggplotGrob(state_ToE_histograms[["QLD"]]),
#  xmin = 145, xmax = 150, ymin = -15, ymax = -10 # dial in the coords
# )
# inset_element(state_ToE_histograms[["QLD"]], left = 0.5, bottom = 0.5, right = 1, top = 1)
# It looks as if geom_magnify treats it as two plots instead of a single plot - it only take the 2nd element in the list
# It must be a single plot (i.e., a list of length 1) - inset_element does not work


#scale_size_binned(limits = scale_size_limits) + # range = c(0, 2) dictates the size of the dots (important)
#guides(size = guide_bins(show.limits = TRUE)) +




inset_plot_NSW <- aus_map |>
  filter(state == "NSW") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = NSW_data,
    aes(x = lon, y = lat, fill = decade_time_of_emergence),#, size = DREAM_ToE_IQR),
    show.legend = FALSE,
    stroke = 0.1,
    alpha = transparent_dots_constant,
    shape = 21,
    colour = "black",
    size = size
  ) +
  scale_fill_brewer(palette = "BrBG", drop = FALSE) +
  #scale_size_binned(limits = scale_size_limits) + # range = c(0, 2) dictates the size of the dots (important)
  #guides(size = guide_bins(show.limits = TRUE)) +
  theme_void()



inset_plot_VIC <- aus_map |>
  filter(state == "VIC") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = VIC_data,
    aes(x = lon, y = lat, fill = decade_time_of_emergence),#, size = DREAM_ToE_IQR),
    show.legend = FALSE,
    alpha = transparent_dots_constant,
    stroke = 0.1,
    shape = 21,
    colour = "black",
    size = size
  ) +
  scale_fill_brewer(palette = "BrBG", drop = FALSE) +
  #scale_size_binned(limits = scale_size_limits) + # range = c(0, 2) dictates the size of the dots (important)
  #guides(size = guide_bins(show.limits = TRUE)) +
  theme_void()



inset_plot_WA <- aus_map |>
  filter(state == "WA") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = WA_data,
    aes(x = lon, y = lat, fill = decade_time_of_emergence),#, size = DREAM_ToE_IQR),
    show.legend = FALSE,
    stroke = 0.1,
    alpha = transparent_dots_constant,
    shape = 21,
    colour = "black",
    size = size
  ) +
  scale_fill_brewer(palette = "BrBG", drop = FALSE) +
  #scale_size_binned(limits = scale_size_limits) + # range = c(0, 2) dictates the size of the dots (important)
  #guides(size = guide_bins(show.limits = TRUE)) +
  theme_void()



inset_plot_TAS <- aus_map |>
  filter(state == "TAS") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = TAS_data,
    aes(x = lon, y = lat, fill = decade_time_of_emergence),#, size = DREAM_ToE_IQR),
    show.legend = FALSE,
    stroke = 0.1,
    alpha = transparent_dots_constant,
    shape = 21,
    colour = "black",
    size = size
  ) +
  scale_fill_brewer(palette = "BrBG", drop = FALSE) +
  #scale_size_binned(limits = scale_size_limits) + # range = c(0, 2) dictates the size of the dots (important)
  #guides(size = guide_bins(show.limits = TRUE)) +
  theme_void()



## Put it together =============================================================
## This probably could be functioned...
ToE_map_aus <- aus_map |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = time_of_emergence_data,
    mapping = aes(x = lon, y = lat, fill = decade_time_of_emergence), #size = DREAM_ToE_IQR),
    stroke = 0.1,
    shape = 21,
    colour = "black",
    size = size
  ) +
  scale_fill_brewer(palette = "BrBG", drop = FALSE) +
  #scale_size_binned(limits = scale_size_limits) + # range = c(0, 2) dictates the size of the dots (important)
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
    x = NULL, # "Latitude",
    y = NULL, # "Longitude",
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
    legend.box = "horizontal" # , # side-by-side legends
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 5, shape = 21), nrow = 3)#, # Wrap legend with nrow
    #size = guide_bins(show.limits = TRUE, direction = "horizontal")
  )


ToE_map_aus

ggsave(
  filename = "./Figures/Main/ToE_map_aus_uncertainty.pdf", #"./Graphs/CMAES_graphs/log_sinh_no_uncertainty_ToE_map_aus.pdf",
  plot = ToE_map_aus,
  device = "pdf",
  width = 232,
  height = 200,
  units = "mm"
)







# Relationships with time of emergence -----------------------------------------

## Time of emergence, record length, uncertainty graph
## See ToE vs. uncertainty in RQ2\Testing-Figures\CMAES-DREAM-old-figures-boxcox\Graphs\Supplementary_Figures



## Time of emergence vs. record length =========================================
record_length_data <- gauge_information |> 
  select(gauge, record_length)

time_of_emergence_data <- time_of_emergence_data |> 
  left_join(
    record_length_data,
    by = join_by(gauge)
  )


time_of_emergence_data |> 
  ggplot(aes(x = record_length, y = year_time_of_emergence)) +
  geom_point() +
  theme_bw()




## Time of emergence vs. catchment area ========================================

catchment_area_data <- gauge_information |> 
  select(gauge, catchment_area_sq_km)


time_of_emergence_data <- time_of_emergence_data |> 
  left_join(
    catchment_area_data,
    by = join_by(gauge)
  )


time_of_emergence_data |> 
  ggplot(aes(x = catchment_area_sq_km, y = year_time_of_emergence)) +
  geom_point() +
  scale_x_log10() +
  theme_bw()














