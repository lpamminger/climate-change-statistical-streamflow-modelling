# Time of emergence analysis

# Figures produced in this R file ----------------------------------------------

# 1. Main --> ToE_map_aus_uncertainty.pdf (requires dream results - currently working without dream results)
# 2. Supplementary --> time_of_emergence_decade_histogram.pdf
# 3. Supplementary --> ToE_vs_uncertainty_and_evidence_ratio.pdf (requires dream results - called ToE_vs_uncertainty in old files)
# 4. Supplementary --> ToE_vs_record_length.pdf 
# 5. Supplementary --> ToE_vs_catchment_area.pdf
# 6. Supplementary --> time_of_emergence_and_uncertainty_plot.pdf
# 7. Supplementary --> time_of_emergence_uncertainty_vs_evidence_ratio.pdf
# 8. Supplementary --> ToE_cdf.pdf 
# X. Testing --> something related to climate types (not in file)






# CODE







# Import libraries -------------------------------------------------------------
pacman::p_load(tidyverse, ozmaps, sf, ggmagnify, furrr, arrow)



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


year_DREAM_sequence_data <- open_dataset(
  source = "./Modelling/Results/DREAM/year_DREAM_sequence_data.parquet",
  format = "parquet"
)



# Calculate time of emergence using CMAES parameters ---------------------------
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



## Add DREAM Time of emergence IQR to time of emergence data ===================
### To save load time only load in gauges that will be used in the plot
filtered_gauges <- time_of_emergence_data |> pull(gauge) |> unique()

DREAM_time_of_emergence_IQR <- year_DREAM_sequence_data |> 
  filter(gauge %in% filtered_gauges) |> 
  select(gauge, ToE) |> 
  collect() |> 
  summarise(
    ToE_IQR = IQR(ToE),
    .by = gauge
  )

time_of_emergence_data <- time_of_emergence_data |> 
  left_join(
    DREAM_time_of_emergence_IQR,
    by = join_by(gauge)
  ) |> 
  # in the plot we want big dots (large uncertainty) on the bottom and small dots (small uncertainty) on top
  # ggplot plots thing in order they appear in the tibble
  arrange(desc(ToE_IQR))




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
scale_size_limits <- time_of_emergence_data |>
  pull(ToE_IQR) |>
  range() # can round up if I want to


transparent_dots_constant <- 0.75
ToE_IQR_breaks <- c(5, 10, 15, 20, 25, 30, 35, 40)


inset_plot_QLD <- aus_map |>
  filter(state == "QLD") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = QLD_data,
    aes(x = lon, y = lat, fill = decade_time_of_emergence, size = ToE_IQR),
    show.legend = FALSE,
    stroke = 0.1,
    alpha = transparent_dots_constant,
    shape = 21, # remove this with size
    colour = "black"
  ) +
  scale_fill_brewer(palette = "BrBG", drop = FALSE) +
  scale_size_binned(limits = scale_size_limits, breaks = ToE_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
  guides(size = guide_bins(show.limits = TRUE)) +
  theme_void() 



# Need to add these histograms as grobs in the inset
# Main plot -> inset map -> inset histogram
# Does patchwork do inset histograms?
# See documentation: https://patchwork.data-imaginist.com/reference/inset_element.html
# Potential code:
# annotation_custom(
#  ggplotGrob(state_ToE_histograms[["QLD"]]),
#  xmin = 145, xmax = 150, ymin = -15, ymax = -10 # dial in the coords
# )
# inset_element(state_ToE_histograms[["QLD"]], left = 0.5, bottom = 0.5, right = 1, top = 1)
# It looks as if geom_magnify treats it as two plots instead of a single plot - it only take the 2nd element in the list
# It must be a single plot (i.e., a list of length 1) - inset_element does not work





inset_plot_NSW <- aus_map |>
  filter(state == "NSW") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = NSW_data,
    aes(x = lon, y = lat, fill = decade_time_of_emergence, size = ToE_IQR),
    show.legend = FALSE,
    stroke = 0.1,
    alpha = transparent_dots_constant,
    shape = 21,
    colour = "black"
  ) +
  scale_fill_brewer(palette = "BrBG", drop = FALSE) +
  scale_size_binned(limits = scale_size_limits, breaks = ToE_IQR_breaks) + 
  guides(size = guide_bins(show.limits = TRUE)) +
  theme_void()



inset_plot_VIC <- aus_map |>
  filter(state == "VIC") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = VIC_data,
    aes(x = lon, y = lat, fill = decade_time_of_emergence, size = ToE_IQR),
    show.legend = FALSE,
    alpha = transparent_dots_constant,
    stroke = 0.1,
    shape = 21,
    colour = "black"
  ) +
  scale_fill_brewer(palette = "BrBG", drop = FALSE) +
  scale_size_binned(limits = scale_size_limits, breaks = ToE_IQR_breaks) +
  guides(size = guide_bins(show.limits = TRUE)) +
  theme_void()



inset_plot_WA <- aus_map |>
  filter(state == "WA") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = WA_data,
    aes(x = lon, y = lat, fill = decade_time_of_emergence, size = ToE_IQR),
    show.legend = FALSE,
    stroke = 0.1,
    alpha = transparent_dots_constant,
    shape = 21,
    colour = "black"
  ) +
  scale_fill_brewer(palette = "BrBG", drop = FALSE) +
  scale_size_binned(limits = scale_size_limits, breaks = ToE_IQR_breaks) + 
  guides(size = guide_bins(show.limits = TRUE)) +
  theme_void()



inset_plot_TAS <- aus_map |>
  filter(state == "TAS") |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = TAS_data,
    aes(x = lon, y = lat, fill = decade_time_of_emergence, size = ToE_IQR),
    show.legend = FALSE,
    stroke = 0.1,
    alpha = transparent_dots_constant,
    shape = 21,
    colour = "black"
  ) +
  scale_fill_brewer(palette = "BrBG", drop = FALSE) +
  scale_size_binned(limits = scale_size_limits, breaks = ToE_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
  guides(size = guide_bins(show.limits = TRUE)) +
  theme_void()



## Put it together =============================================================
## This probably could be functioned...
ToE_map_aus <- aus_map |>
  ggplot() +
  geom_sf() +
  geom_point(
    data = time_of_emergence_data,
    mapping = aes(x = lon, y = lat, fill = decade_time_of_emergence, size = ToE_IQR),
    stroke = 0.1,
    shape = 21,
    colour = "black"
  ) +
  scale_fill_brewer(palette = "BrBG", drop = FALSE) +
  scale_size_binned(limits = scale_size_limits, breaks = ToE_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
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
    #legend.key = element_rect(fill = "grey80"),
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
    fill = guide_legend(override.aes = list(size = 5, shape = 21), nrow = 3), # Wrap legend with nrow
    size = guide_bins(show.limits = TRUE, direction = "horizontal")
  )


ToE_map_aus

ggsave(
  filename = "./Figures/Main/ToE_map_aus_uncertainty.pdf", 
  plot = ToE_map_aus,
  device = "pdf",
  width = 232,
  height = 200,
  units = "mm"
)





# Time of Emergence Empirical Cumulative Distribution --------------------------
# Construct Empirical Cumulative Distribution for each catchment using
# time of emergence data

# Method:
# 1. Pull ToE data for each gauge
# 2. Put ToE data into ecdf to produce a ecdf function
# 3. Use ecdf function to produce y-axis of cdf (xaxis remains constant)
# 4. Average y-axis (for every value of x) for each state

year_DREAM_ToE_data <- year_DREAM_sequence_data |> 
  filter(gauge %in% filtered_gauges) |> 
  select(gauge, ToE) |> 
  collect()

ecdf_function_factory <- function(gauge, data) {
  data |> 
    filter(gauge == {{ gauge }}) |>   
    pull(ToE) |> 
    ecdf()
}

gauge_ecdf_functions <- map(
  .x = filtered_gauges,
  .f = ecdf_function_factory,
  data = year_DREAM_ToE_data
)


# You cannot map over function - make a function that allows this
map_ecdf <- function(specific_ecdf, x_axis) {
  specific_ecdf(x_axis)
}

cdf_x_axis <- seq(from = 1959, to = 2010, by = 1)

gauge_cdf_results <- map(
  .x = gauge_ecdf_functions,
  .f = map_ecdf,
  x_axis = cdf_x_axis
) |> 
  `names<-`(filtered_gauges) |> 
  as_tibble() |> 
  add_column(
    "cdf_x_axis" = cdf_x_axis,
    .before = 1
  ) 


## Construct state cdf's =======================================================
state_cdf_results <- gauge_cdf_results |>
  pivot_longer(
    cols = !cdf_x_axis,
    names_to = "gauge",
    values_to = "cdf_value"
  ) |> 
  # add state
  left_join(
    lat_lon_gauge,
    by = join_by(gauge)
  ) |> 
  # summarise by state
  summarise(
    median_cdf = median(cdf_value),
    lower_bound_cdf = quantile(cdf_value, probs = 0.1),
    upper_bound_cdf = quantile(cdf_value, probs = 0.9),
    n = n(),
    .by = c(cdf_x_axis, state)
  ) 


## add n = ??? to plots
count_state_cdf <- state_cdf_results |> 
  filter(cdf_x_axis == 1959) |> 
  select(state, n) |>
  mutate(
    n_label = paste0("n = ", n)
  ) |> 
  add_column(
    "x_pos" = 1963,
    "y_pos" = 0.9
  )


## Plot ========================================================================
state_cdf_plot <- state_cdf_results |> 
  ggplot(aes(x = cdf_x_axis, y = median_cdf)) +
  geom_step() +
  geom_label(
    aes(x = x_pos, y = y_pos, label = n_label),
    data = count_state_cdf
  ) +
  geom_ribbon(aes(ymin = lower_bound_cdf, ymax = upper_bound_cdf), alpha = 0.2) +
  labs(x = "Time of Emergence (Year)", y = "Cumulative Probabiliy") +
  theme_bw() +
  facet_wrap(~state)


ggsave(
  filename = "./Figures/Supplementary/ToE_cdf.pdf", 
  plot = state_cdf_plot,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)







# Decade histogram time of emergence -------------------------------------------
time_of_emergence_histogram <- time_of_emergence_data |> 
  ggplot(aes(x = decade_time_of_emergence)) +
  geom_histogram(stat = "count") +
  labs(
    x = "Time of Emergence Decade",
    y = "Frequency"
  ) +
  theme_bw()


ggsave(
  filename = "./Figures/Supplementary/time_of_emergence_decade_histogram.pdf", 
  plot = time_of_emergence_histogram,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)







# Relationships with time of emergence -----------------------------------------

## Time of emergence vs. uncertainty ===========================================
time_of_emergence_and_uncertainty_plot <- time_of_emergence_data |> 
  ggplot(aes(x = year_time_of_emergence, y = ToE_IQR)) +
  geom_point() +
  facet_wrap(~binned_evidence_ratio) +
  labs(
    x = "Time of Emergence (Year)",
    y = "Time of Emergence Interquantile Range"
  ) +
  theme_bw()


ggsave(
  filename = "./Figures/Supplementary/time_of_emergence_and_uncertainty_plot.pdf", 
  plot = time_of_emergence_and_uncertainty_plot,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)


## Time of emergence uncertainty vs. evidence ratio ============================
time_of_emergence_uncertainty_vs_evidence_ratio <- time_of_emergence_data |> 
  ggplot(aes(x = evidence_ratio, y = ToE_IQR)) +
  geom_point() +
  scale_x_log10() +
  labs(
    x = "Evidence Ratio",
    y = "Time of Emergence Interquantile Range"
  ) +
  theme_bw()


ggsave(
  filename = "./Figures/Supplementary/time_of_emergence_uncertainty_vs_evidence_ratio.pdf", 
  plot = time_of_emergence_uncertainty_vs_evidence_ratio,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)


## Time of emergence vs. record length =========================================
record_length_data <- gauge_information |> 
  select(gauge, record_length)

time_of_emergence_data <- time_of_emergence_data |> 
  left_join(
    record_length_data,
    by = join_by(gauge)
  )


ToE_vs_record_length <- time_of_emergence_data |> 
  ggplot(aes(x = record_length, y = year_time_of_emergence)) +
  geom_point() +
  theme_bw()



ggsave(
  filename = "./Figures/Supplementary/ToE_vs_record_length.pdf", #"./Graphs/CMAES_graphs/log_sinh_no_uncertainty_ToE_map_aus.pdf",
  plot = ToE_vs_record_length,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)


## Time of emergence vs. catchment area ========================================

catchment_area_data <- gauge_information |> 
  select(gauge, catchment_area_sq_km)


time_of_emergence_data <- time_of_emergence_data |> 
  left_join(
    catchment_area_data,
    by = join_by(gauge)
  )


ToE_vs_catchment_area <- time_of_emergence_data |> 
  ggplot(aes(x = catchment_area_sq_km, y = year_time_of_emergence)) +
  geom_point() +
  scale_x_log10() +
  theme_bw()


ggsave(
  filename = "./Figures/Supplementary/ToE_vs_catchment_area.pdf", #"./Graphs/CMAES_graphs/log_sinh_no_uncertainty_ToE_map_aus.pdf",
  plot = ToE_vs_catchment_area,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)











