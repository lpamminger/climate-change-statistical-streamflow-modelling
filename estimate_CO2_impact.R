## Fitting streamflow models to catchments
cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, truncnorm, sloop, patchwork, ozmaps, sf, patchwork, metR, ggmagnify)



# Import and prepare data-------------------------------------------------------

## Import annual streamflow, precip, CO2 and gauge data ========================
start_stop_indexes <- readr::read_csv(
  "./Data/Tidy/start_end_index.csv",
  show_col_types = FALSE
)

data <- readr::read_csv(
  "./Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  col_types = "icdddddddll", # ensuring each column is a the correct type
  show_col_types = FALSE
)

gauge_information <- readr::read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE
)

lat_lon_gauge_info <- gauge_information |> 
  select(gauge, lat, lon, state)


CMAES_results <- read_csv("./Results/my_cmaes/CMAES_parameter_results_20250331.csv",
  show_col_types = FALSE
)

best_CO2_non_CO2_per_gauge <- read_csv(
  "./Results/my_cmaes/unmodified_best_CO2_non_CO2_per_catchment_CMAES_20250331.csv",
  show_col_types = FALSE
)

streamflow_data <- read_csv(
  "Results/my_cmaes/CMAES_streamflow_results_20250331.csv",
  show_col_types = FALSE
)


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



# TODO:
# Methods of comparison
# - streamflow timeseries 
# - rainfall-runoff relationship
# - average yearly percentage difference


# Compare best CO2 and non-CO2 models ------------------------------------------
## Summarise modelled streamflow data ==========================================
only_gauge_model_best_CO2_non_CO2_per_gauge <- best_CO2_non_CO2_per_gauge |>
  select(gauge, streamflow_model) |>
  distinct()



streamflow_data_best_CO2_non_CO2 <- streamflow_data |>
  semi_join(
    only_gauge_model_best_CO2_non_CO2_per_gauge,
    by = join_by(gauge, streamflow_model)
  )



## Are the best CO2 and non-CO2 models equivalent? ==============================
compare_equivalent_models <- only_gauge_model_best_CO2_non_CO2_per_gauge |> 
  mutate(
    CO2_or_non_CO2_model = if_else(str_detect(streamflow_model, "CO2"), "CO2_model", "non_CO2_model")
  ) |> 
  pivot_wider(
    names_from = CO2_or_non_CO2_model,
    values_from = streamflow_model
  )


### Comparison table ###
streamflow_model_equivalent_comparison_table <- tribble(
  ~non_CO2_model,                                        ~CO2_model, 
  "streamflow_model_precip_only",                        "streamflow_model_intercept_shifted_CO2",
  "streamflow_model_precip_only",                        "streamflow_model_slope_shifted_CO2",
  "streamflow_model_precip_auto",                        "streamflow_model_intercept_shifted_CO2_auto",
  "streamflow_model_precip_auto",                        "streamflow_model_slope_shifted_CO2_auto",
  "streamflow_model_precip_seasonal_ratio",              "streamflow_model_intercept_shifted_CO2_seasonal_ratio",
  "streamflow_model_precip_seasonal_ratio",              "streamflow_model_slope_shifted_CO2_seasonal_ratio",
  "streamflow_model_precip_seasonal_ratio_auto",         "streamflow_model_intercept_shifted_CO2_seasonal_ratio_auto",
  "streamflow_model_precip_seasonal_ratio_auto",         "streamflow_model_slope_shifted_CO2_seasonal_ratio_auto",
  "streamflow_model_drought_precip_only",                "streamflow_model_drought_intercept_shifted_CO2",
  "streamflow_model_drought_precip_only",                "streamflow_model_drought_slope_shifted_CO2",
  "streamflow_model_drought_precip_auto",                "streamflow_model_drought_intercept_shifted_CO2_auto",
  "streamflow_model_drought_precip_auto",                "streamflow_model_drought_slope_shifted_CO2_auto",
  "streamflow_model_drought_precip_seasonal_ratio",      "streamflow_model_drought_intercept_shifted_CO2_seasonal_ratio",
  "streamflow_model_drought_precip_seasonal_ratio",      "streamflow_model_drought_slope_shifted_CO2_seasonal_ratio",
  "streamflow_model_drought_precip_seasonal_ratio_auto", "streamflow_model_drought_intercept_shifted_CO2_seasonal_ratio_auto",
  "streamflow_model_drought_precip_seasonal_ratio_auto", "streamflow_model_drought_slope_shifted_CO2_seasonal_ratio_auto"
  
)



### For each row in compare_equivalent_models check if column CO2_model and non_CO2_model is in the comparison table
### Do this by using anti_join and counting row numbers
non_equivalent_models <- compare_equivalent_models |> 
  anti_join(
    streamflow_model_equivalent_comparison_table,
    by = join_by(non_CO2_model, CO2_model)
  ) 

(nrow(non_equivalent_models) / nrow(compare_equivalent_models)) * 100


# Compare best CO2 with CO2 component turned off -------------------------------
## Get catchments where the CO2 model is the best ==============================
best_model_per_gauge <- best_CO2_non_CO2_per_gauge |> 
  slice_min(
    AIC,
    by = gauge
    )

only_CO2_best_models <- best_model_per_gauge |> 
  mutate(
    is_CO2_model = str_detect(streamflow_model, "CO2")
  ) |> 
  filter(is_CO2_model) |> 
  select(!is_CO2_model)


## Set a3 terms (slope and intercept to zero) ==================================
set_a3_zero_CO2_best_models <- only_CO2_best_models |> 
  mutate(
    altered_parameter = if_else(str_detect(parameter, "a3"), 0, parameter_value)
  )
  

## Generate streamflow with altered_parameter ==================================
### Build catchment_dataset objects ############################################
CO2_gauges <- set_a3_zero_CO2_best_models |> 
  pull(gauge) |> 
  unique()

CO2_catchment_data <- map(
  .x = CO2_gauges,
  .f = catchment_data_blueprint,
  observed_data = data,
  start_stop_indexes = start_stop_indexes
)

### Run the streamflow models ##################################################
### streamflow_models 
get_model_name_and_function <- function(model) {
  model_name <- model()$name
  tibble(
    "streamflow_model" = model_name,
    "model_function" = list(model)
  )
}


streamflow_name_and_model_for_joining <- map(
  .x = c(get_non_drought_streamflow_models(), get_drought_streamflow_models()),
  .f = get_model_name_and_function
) |>
  list_rbind() 



zero_CO2_streamflow_models <- set_a3_zero_CO2_best_models |> 
  select(gauge, streamflow_model) |> 
  distinct() |> 
  left_join(
    streamflow_name_and_model_for_joining,
    by = join_by(streamflow_model)
  ) |> 
  pull(model_function)


### Parameter sets
zero_CO2_parameter_sets <- set_a3_zero_CO2_best_models |> 
  select(gauge, parameter, altered_parameter)

# It might be worth generalising this...
tibble_column_to_parameter_vector <- function(gauge, data) {
  data |> 
    filter(gauge == {{ gauge }}) |> 
    pull(altered_parameter)
} 

altered_parameter_sets <- map(
  .x = CO2_gauges,
  .f = tibble_column_to_parameter_vector,
  data = set_a3_zero_CO2_best_models
)


### Put it together now 
generate_altered_streamflow <- function(catchment_data, parameter_set, streamflow_model) {
  streamflow_model <- noquote(streamflow_model)
  streamflow_model(catchment_data, parameter_set)
}


altered_streamflow <- pmap(
  .l = list(CO2_catchment_data, altered_parameter_sets, zero_CO2_streamflow_models),
  .f = generate_altered_streamflow
) |>  
  list_rbind() 

gauge_join <- data |> 
  select(gauge, year, p_mm) |> 
  rename(
    precipitation = p_mm
  )

###  Summary table #############################################################
no_CO2_data_joining <- set_a3_zero_CO2_best_models |> 
  select(gauge, streamflow_model) |> 
  distinct()

only_CO2_streamflow_data <- streamflow_data |> 
  semi_join(
    no_CO2_data_joining,
    by = join_by(gauge, streamflow_model)
  ) |> 
  select(gauge, year, precipitation, modelled_boxcox_streamflow)

streamflow_data_a3_off <- altered_streamflow |> 
  left_join(
    gauge_join,
    by = join_by(year, precipitation)
  ) |> 
  relocate(
    gauge,
    .before = 1
  ) |> 
# add the original streamflow on
rename(
  a3_off_modelled_boxcox_streamflow = modelled_boxcox_streamflow
) |> 
  left_join(
    only_CO2_streamflow_data,
    by = join_by(gauge, year, precipitation)
  )



# Compare percentage difference in a3 on vs. a3 off ----------------------------
## Convert boxcox streamflow to actual streamflow ==============================
## this must be done so summing each decade does not result in negative flow
boxcox_lambda_gauge <- gauge_information |>
  select(gauge, bc_lambda)

realspace_streamflow_data_a3_off <- streamflow_data_a3_off |>
  left_join(
    boxcox_lambda_gauge,
    by = join_by(gauge)
  ) |>
  mutate(
    realspace_a3_on = boxcox_inverse_transform(yt = modelled_boxcox_streamflow, lambda = bc_lambda, lambda_2 = 1),
    realspace_a3_off = boxcox_inverse_transform(yt = a3_off_modelled_boxcox_streamflow, lambda = bc_lambda, lambda_2 = 1),
    .by = gauge
  )


## Check for strange streamflow values
## - There are some streamflow values that are negative or NA
observed_streamflow <- data |> 
  select(gauge, p_mm, q_mm) |> 
  rename(
    precipitation = p_mm
  )

problem_flows <- realspace_streamflow_data_a3_off |> 
  filter(is.na(realspace_a3_on) | is.na(realspace_a3_off) | (realspace_a3_on < 0) | (realspace_a3_off < 0)) |> 
  left_join(
    observed_streamflow,
    by = join_by(gauge, precipitation)
  )


## Solution to problem flows
## - All occur when annual streamflow is near zero.
## - For simplicity set to zero
## - Should there a lower bound for years in decade?

average_percent_diff_by_decade <- realspace_streamflow_data_a3_off |>
  # TEMP - set weird values to zero
  mutate(
    realspace_a3_off = case_when(
      realspace_a3_off < 0 ~ 0,
      is.na(realspace_a3_off) ~ 0,
      .default = realspace_a3_off
    ),
    realspace_a3_on = case_when(
      realspace_a3_on < 0 ~ 0,
      is.na(realspace_a3_on) ~ 0,
      .default = realspace_a3_on
    )
  ) |>
  mutate(
    decade = year - (year %% 10)
  ) |>
  summarise(
    sum_streamflow_a3_off = sum(realspace_a3_off),
    sum_streamflow_a3_on = sum(realspace_a3_on),
    n = n(),
    .by = c(gauge, decade)
  ) |>
  mutate(
    average_diff = ((sum_streamflow_a3_on - sum_streamflow_a3_off) / sum_streamflow_a3_on) * 100
  ) |>
  arrange(average_diff) |>
  left_join(
    lat_lon_gauge_info,
    by = join_by(gauge)
  ) |> 
  # only include decades with 10 years - removes some weird behaving catchment
  filter(n == 10) |> 
  # turning CO2 off makes streamflow go negative?
  filter(sum_streamflow_a3_off > 0) 


average_percent_diff_by_decade |> pull(average_diff) |> range()
average_percent_diff_by_decade |> pull(average_diff) |> mean()

## Count gauges per decade
average_percent_diff_by_decade |> 
  summarise(
    n = n(),
    .by = decade
  )



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




## Plot map percentage difference ==============================================

### REPLACE WITH FULL FAT ###
# 1990 vs. 2010

# 1990
# Make final plot --------------------------------------------------------------

### Custom colour palette 
big_palette <- function(x) {
  c("#67001f",
    "#b2182b",
    "#d6604d",
    "#f4a582",
    "#fddbc7",
    "white",
    "white",
    "#d1e5f0",
    "#92c5de",
    "#4393c3",
    "#2166ac",
    "#053061")
}


decade_comparison_CO2_impact <- function(decade) {
  ## Generate Insets =============================================================
  ### Filter data by state #######################################################
  
  QLD_data <- average_percent_diff_by_decade |>
    filter(state == "QLD") |> 
    filter(decade == {{ decade }})
  
  NSW_data <- average_percent_diff_by_decade |>
    filter(state == "NSW") |> 
    filter(decade == {{ decade }})
  
  VIC_data <- average_percent_diff_by_decade |>
    filter(state == "VIC") |> 
    filter(decade == {{ decade }})
  
  WA_data <- average_percent_diff_by_decade |>
    filter(state == "WA") |> 
    filter(decade == {{ decade }})
  
  TAS_data <- average_percent_diff_by_decade |>
    filter(state == "TAS") |> 
    filter(decade == {{ decade }})
  
  
  ### Generate inset plots #######################################################
  
  inset_plot_QLD <- aus_map |>
    filter(state == "QLD") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = QLD_data,
      aes(x = lon, y = lat, fill = average_diff),
      show.legend = FALSE,
      size = 2.5,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = c(-50, -25, -5, -1, -0.1, 0, 0.1, 1, 5, 25, 50), # range()
      limits = c(-1100, 90),
      show.limits = TRUE, 
      guide = "colorsteps"
    ) +
    theme_void()
  
  
  inset_plot_NSW <- aus_map |>
    filter(state == "NSW") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = NSW_data,
      aes(x = lon, y = lat, fill = average_diff),
      show.legend = FALSE,
      size = 2.5,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = c(-50, -25, -5, -1, -0.1, 0, 0.1, 1, 5, 25, 50), # range()
      limits = c(-1100, 90),
      show.limits = TRUE, 
      guide = "colorsteps"
    ) +
    theme_void()
  
  
  
  inset_plot_VIC <- aus_map |>
    filter(state == "VIC") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = VIC_data,
      aes(x = lon, y = lat, fill = average_diff),
      show.legend = FALSE,
      size = 2.5,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = c(-50, -25, -5, -1, -0.1, 0, 0.1, 1, 5, 25, 50), # range()
      limits = c(-1100, 90),
      show.limits = TRUE, 
      guide = "colorsteps"
    ) +
    theme_void()
  
  
  
  inset_plot_WA <- aus_map |>
    filter(state == "WA") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = WA_data,
      aes(x = lon, y = lat, fill = average_diff),
      show.legend = FALSE,
      size = 2.5,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = c(-50, -25, -5, -1, -0.1, 0, 0.1, 1, 5, 25, 50), # range()
      limits = c(-1100, 90),
      show.limits = TRUE, 
      guide = "colorsteps"
    ) +
    theme_void()
  
  
  
  inset_plot_TAS <- aus_map |>
    filter(state == "TAS") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = TAS_data,
      aes(x = lon, y = lat, fill = average_diff),
      show.legend = FALSE,
      size = 2.5,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = c(-50, -25, -5, -1, -0.1, 0, 0.1, 1, 5, 25, 50), # range()
      limits = c(-1100, 90),
      show.limits = TRUE, 
      guide = "colorsteps"
    ) +
    theme_void()
  
  
  
  ## Put it together =============================================================
  
  single_map_aus <- aus_map |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = average_percent_diff_by_decade |> filter(decade == {{ decade }}),
      mapping = aes(x = lon, y = lat, fill = average_diff),
      size = 3,
      colour = "black",
      shape = 21,
      inherit.aes = FALSE,
      stroke = 0.1
    ) +
    theme_bw() +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = c(-50, -25, -5, -1, -0.1, 0, 0.1, 1, 5, 25, 50), # range()
      limits = c(-1100, 90),
      show.limits = TRUE, 
      guide = "colorsteps"
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
      fill = "Percentage Difference",
      title = paste0(decade)
    ) +
    theme(
      legend.key = element_rect(fill = "grey80"),
      legend.title = element_text(hjust = 0.5),
      #legend.background = element_rect(colour = "black"),
      axis.text = element_blank(), 
      legend.position = "inside",
      legend.position.inside = c(0.351, 0.9),
      legend.box = "horizontal", # side-by-side legends
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(margin = margin(l = 25, r = 0, t = 30, b = -30), size = 22) # push title into plot
    ) +
    guides(
      fill = guide_coloursteps(
        barwidth = unit(15, "cm"), 
        show.limits = TRUE, 
        even.steps = TRUE,
        title.position = "top",
        direction = "horizontal"
      )
    ) 
  
  return(single_map_aus)
}


# Patchwork together
patchwork_CO2_impact_decade <- (decade_comparison_CO2_impact(1990) | decade_comparison_CO2_impact(2010)) + plot_layout(guides = "collect") & theme(legend.position='bottom')

ggsave(
  filename = "./Graphs/CMAES_graphs/CO2_on_off_v5.pdf",
  plot = patchwork_CO2_impact_decade,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)


## OLD --> REMOVE ##

### Map
single_aus_map <- ozmaps::ozmap("country") |> 
  uncount(average_percent_diff_by_decade |> pull(decade) |> unique() |> length()) |>   # repeat the geometry by number of decades in average_percent_diff_by_decade
  mutate(
    decade = average_percent_diff_by_decade |> pull(decade) |> unique()
  )


big_palette <- function(x) {
  c("#67001f",
  "#b2182b",
  "#d6604d",
  "#f4a582",
  "#fddbc7",
  "white",
  "white",
  "#d1e5f0",
  "#92c5de",
  "#4393c3",
  "#2166ac",
  "#053061")
}


plot_CO2_on_off_percent_diff <- single_aus_map |> 
  ggplot(aes(geometry = geometry)) +
  geom_sf(
    colour = "black",
    fill = "grey80"
  ) +
  geom_point(
    mapping = aes(x = lon, y = lat, fill = average_diff),
    data = average_percent_diff_by_decade,
    inherit.aes = FALSE,
    size = 1.5,
    shape = 21,
    colour = "black",
    stroke = 0.1
  ) + 
  coord_sf(xlim = c(111, 155), ylim = c(-44.5, -9.5)) +
  metR::scale_x_longitude(ticks = 10) +
  metR::scale_y_latitude(ticks = 10) +
  labs(
    x = "Longitude",
    y = "Latitude",
    colour = NULL
  ) +
  binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
    aesthetics = "fill",
    palette = big_palette,
    breaks = c(-50, -25, -5, -1, -0.1, 0, 0.1, 1, 5, 25, 50), # range()
    limits = c(-1100, 90),
    show.limits = TRUE, 
    guide = "colorsteps"
    ) +
  facet_wrap(~decade) +
  labs(
    fill = "Mean Percentage Difference (CO2_on - CO2_off)"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    #legend.position.inside = c(0.7, 0.15),
    legend.title = element_text(size = 10, hjust = 0.5)
  ) +
  guides(
    fill = guide_coloursteps(
      barwidth = unit(15, "cm"), 
      show.limits = TRUE, 
      even.steps = TRUE,
      title.position = "top",
      direction = "horizontal"
      )
    ) 

#plot_CO2_on_off_percent_diff

ggsave(
  filename = "./Graphs/CMAES_graphs/CO2_on_off_v4.pdf",
  plot = plot_CO2_on_off_percent_diff,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)

### END OF OLD ###











stop_here <- tactical_typo()




# Compare evidence ratio for best non-CO2 vs. CO2 model ------------------------
# Method:
# 1. streamflow_data_best_CO2_non_CO2
# 2. Find percentage difference for each gauge best_Co2_model - best_non_CO2_model / best_CO2_model
# 3. summarise percentage difference
# 4. Join with evidence ratio
# 5. plot

## Evidence ratio for joining ==================================================
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



## percentage difference calculation ===========================================
calculate_percent_diff_best_non_CO2_and_CO2 <- function(gauge, best_CO2_non_CO2_data) {
  
  # Pivot wider for mutate
  wide_streamflow_data <- best_CO2_non_CO2_data |> 
    select(!c(objective_function, optimiser)) |> 
    filter(gauge == {{ gauge }}) |> 
    pivot_wider(
      names_from = streamflow_model,
      values_from = modelled_boxcox_streamflow
    ) 
  
  # Rename columns for generic percent diff
  column_names <- names(wide_streamflow_data)
  CO2_model_name_index <- str_detect(names(wide_streamflow_data), "streamflow_model") & str_detect(names(wide_streamflow_data), "CO2")
  no_CO2_model_name_index <- str_detect(names(wide_streamflow_data), "streamflow_model") & !str_detect(names(wide_streamflow_data), "CO2")
  column_names[CO2_model_name_index] <- "CO2_model"
  column_names[no_CO2_model_name_index] <- "non_CO2_model"
  
  # Percent diff calculation - sum everything then percent diff
  calculation_streamflow_data <- wide_streamflow_data |> 
    `names<-`(column_names) |>  
    summarise(
      sum_CO2_model = sum(CO2_model),
      sum_non_CO2_model = sum(non_CO2_model),
      .by = gauge
    ) |> 
    mutate(
      percent_diff = ((sum_CO2_model - sum_non_CO2_model) / sum_CO2_model) * 100
    )
  
  return(calculation_streamflow_data)
}


percent_diff_streamflow_data_best_CO2_non_CO2 <- map(
  .x = streamflow_data_best_CO2_non_CO2 |> pull(gauge) |> unique(),
  .f = calculate_percent_diff_best_non_CO2_and_CO2,
  best_CO2_non_CO2_data = streamflow_data_best_CO2_non_CO2
) |> 
  list_rbind()


## Summarise percentage difference and join evidence ratio =====================
summarise_percent_diff_streamflow_data_best_CO2_non_CO2 <- percent_diff_streamflow_data_best_CO2_non_CO2 |> 
  left_join(
    evidence_ratio_calc,
    by = join_by(gauge)
  )
  
  

## Plot percentage_diff vs. evidence ratio =====================================
### 

summarise_percent_diff_streamflow_data_best_CO2_non_CO2 |> 
  ggplot(aes(x = evidence_ratio, y = percent_diff)) +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(
    x = "Evidence Ratio (log10)",
    y = "Mean Percentage Difference in Modelled Streamflow"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )




# Compare non-CO2 vs. CO2 and a3-on vs a3-off ----------------------------------

# Would be a good idea to get streamflow_data_a3_off and streamflow_data_best_CO2_non_CO2
# in the same format

# suggested format?

gauge_key <- "121003A"  #"A2390531" #"226410" #"403226"

## Rainfall-runoff relationship
### CO2 vs. non CO2
test_1 <- streamflow_data_best_CO2_non_CO2 |> 
  filter(gauge == {{ gauge_key }}) |> 
  select(!c(objective_function, optimiser)) |> 
  mutate(
    is_CO2_model = str_detect(streamflow_model, "CO2")
  ) |> 
  pivot_longer(
    cols = ends_with("streamflow"),
    names_to = "modelled_or_observed",
    values_to = "streamflow"
  ) |>
  mutate(
    modelled_or_observed = if_else(is_CO2_model, paste0("CO2_", modelled_or_observed), paste0("no_CO2_", modelled_or_observed))
  ) |> 
  mutate(
    modelled_or_observed = if_else(str_detect(modelled_or_observed, "observed"), "observed_boxcox_streamflow", modelled_or_observed)
  )
  # remove duplicated observed values
  #filter(!is_CO2_model & modelled_or_observed != "modelled_boxcox_streamflow")
  

graph_CO2_vs_non_CO2 <- test_1 |> 
  mutate(
    modelled_or_observed = case_when(
      modelled_or_observed == "observed_boxcox_streamflow" ~ "Observed Box-Cox Streamflow",
      modelled_or_observed == "CO2_modelled_boxcox_streamflow" ~ "Best CO2 Modelled Box-Cox Streamflow",
      modelled_or_observed == "no_CO2_modelled_boxcox_streamflow" ~ "Best non-CO2 Modelled Box-Cox Streamflow",
      .default = NA
    )
  ) |> 
  ggplot(aes(x = precipitation, y = streamflow, colour = modelled_or_observed)) +
  geom_point(size = 2) +
  geom_smooth(
    formula = y ~ x,
    method = lm,
    se = FALSE
  ) +
  labs(
    x = "Precipitation (mm/year)",
    y = "Box-Cox Streamflow (mm/year)",
    colour = NULL,
    title = "Best CO2 modelled vs. best non-CO2 model"
  ) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.35, 0.9),
    legend.background = element_rect(colour = "black"),
    plot.title = element_text(hjust = 0.5)
  )
  

timeseries_graph_CO2_vs_non_CO2 <- test_1 |> 
  ggplot(aes(x = year, y = streamflow, colour = modelled_or_observed)) +
  geom_point(size = 2) +
  geom_line() +
  labs(
    x = "Year",
    y = "Box-Cox Streamflow",
    colour = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.9, 0.9)
  )

timeseries_graph_CO2_vs_non_CO2

### a3 off vs. a3 on
test_2 <- streamflow_data_a3_off |> 
  filter(gauge == {{ gauge_key }}) |> 
  pivot_longer(
    cols = ends_with("streamflow"),
    names_to = "modelled_or_observed",
    values_to = "boxcox_streamflow"
  )

graph_CO2_off_vs_CO2_on <- test_2 |> 
  mutate(
    modelled_or_observed = case_when(
      modelled_or_observed == "observed_boxcox_streamflow" ~ "Observed Box-Cox Streamflow",
      modelled_or_observed == "modelled_boxcox_streamflow" ~ "Best CO2 Modelled Box-Cox Streamflow",
      modelled_or_observed == "a3_off_modelled_boxcox_streamflow" ~ "CO2 Component Turned Off Modelled Box-Cox Streamflow",
      .default = NA
    )
  ) |> 
  ggplot(aes(x = precipitation, y = boxcox_streamflow, colour = modelled_or_observed)) +
  geom_point(size = 2) +
  geom_smooth(
    formula = y ~ x,
    method = lm,
    se = FALSE
  ) +
  theme_bw() +
  labs(
    x = "Precipitation (mm/year)",
    y = "Box-Cox Streamflow (mm/year)",
    colour = NULL,
    title = "Turning CO2 Component Off"
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.4, 0.9),
    legend.background = element_rect(colour = "black"),
    plot.title = element_text(hjust = 0.5)
  )


save_rainfall_runoff_test_plot <- graph_CO2_vs_non_CO2 + graph_CO2_off_vs_CO2_on
ggsave(
  filename = "rainfall_runoff_comparison.pdf",
  plot = save_rainfall_runoff_test_plot,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)






# Compare 3 methods to estimate the impact of CO2 on streamflow ----------------
## Rank all the models by gauge (CMAES results)
ranked_models_per_gauge <- CMAES_results |>
  group_by(gauge) |>
  arrange(AIC) |>
  ungroup()

## Comparing the three methods
### Method 1: set a3 to zero
### Method 2: compare best CO2 to best non-CO2
### Method 3: compare best CO2 to equivalent non-CO2 model (could also be best non-CO2)

### How to compare the methods:
### 1. streamflow-time with observed
### 2. rainfall-runoff with observed
### 3. Difference between models vs. evidence ratio - not sure if this will work
###    with method 1 - I could do difference between observed - modelled vs. evi ratio

## Method 2 ====================================================================


### Difference in streamflow models ############################################
### Difference in streamflow models per year mean(mm/y) vs. evidence ratio or AIC
### Use best_CO2_non_CO2_per_gauge to filter streamflow_data
### Pivot wider, mutate model 1 - model 2 -> summarise mean .by gauge
### y axis-options:
# mean absolute difference in streamflow (mm/y)
# total difference in streamflow over the entire period
# mean percentage difference in streamflow per year




# percentage_diff <- function(x, y) {
#  (abs(x - y) / ((x + y) / 2)) * 100
# }

### Difference in streamflow for the models
difference_streamflow_per_year_best_CO2_non_CO2 <- streamflow_data_best_CO2_non_CO2 |>
  # Rename streamflow models to non-CO2 and CO2
  mutate(
    CO2_or_non_CO2 = if_else(str_detect(streamflow_model, "CO2"), "CO2_streamflow", "non_CO2_streamflow")
  ) |>
  # Remove streamflow model to make CO2 and non-CO2 on same row
  select(!c(streamflow_model, objective_function, optimiser)) |>
  pivot_wider(
    names_from = CO2_or_non_CO2,
    values_from = modelled_boxcox_streamflow
  ) |>
  mutate(
    yearly_CO2_non_CO2_difference = CO2_streamflow - non_CO2_streamflow,
    percentage_yearly_CO2_non_CO2_difference = ((CO2_streamflow - non_CO2_streamflow) / CO2_streamflow) * 100 # ,
    # alternative_percentage = percentage_diff(CO2_streamflow, non_CO2_streamflow)
  )

### Summarise yearly differences
summary_streamflow_best_CO2_non_CO2 <- difference_streamflow_per_year_best_CO2_non_CO2 |>
  summarise(
    mean_yearly_CO2_non_CO2_difference = mean(yearly_CO2_non_CO2_difference),
    sum_yearly_CO2_non_CO2_difference = sum(yearly_CO2_non_CO2_difference),
    mean_percent_yearly_CO2_non_CO2_difference = mean(percentage_yearly_CO2_non_CO2_difference),
    sum_CO2_streamflow = sum(CO2_streamflow),
    sum_non_CO2_streamflow = sum(non_CO2_streamflow),
    .by = gauge
  ) |>
  # add evidence ratio
  left_join(
    evidence_ratio_calc,
    by = join_by(gauge)
  ) |>
  # percentage diff

  arrange(sum_yearly_CO2_non_CO2_difference)

### Best CO2 vs. non-CO2 model comparison ######################################
summary_streamflow_best_CO2_non_CO2 |>
  filter(evidence_ratio > 0) |>
  # filter(abs(mean_percent_yearly_CO2_non_CO2_difference) < 1) |>
  ggplot(aes(x = evidence_ratio, y = mean_percent_yearly_CO2_non_CO2_difference)) +
  geom_point() +
  scale_x_log10() +
  theme_bw() +
  labs(
    y = "mean(CO2_t - non-CO2_t / CO2_t)"
  )

# This graph tells us
# x-axis is evidence ratio between CO2 and non-CO2 (log10 scaled)
# x-axis values < 0 removed for log10 scale. This removes catchments
# where the non-CO2 model is better
# the CO2 model on average produces -0.17 % less streamflow per year
# compared to the non-CO2 models
# There is no relationship between evidence ratio and difference in models
# We expected to see a larger difference in CO2 and non-CO2 models as
# the evidence ratio grew.

# We really need to compare it to the observed...

summary_streamflow_best_CO2_non_CO2 |>
  pull(mean_percent_yearly_CO2_non_CO2_difference) |>
  quantile()