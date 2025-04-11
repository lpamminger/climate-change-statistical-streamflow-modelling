## Fitting streamflow models to catchments
cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, truncnorm, sloop, patchwork, ozmaps, sf, patchwork, metR)



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
  select(gauge, lat, lon)


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
CO2_on_off_analysis <- streamflow_data_a3_off |> 
  mutate(
    # CO2 on in the model is the baseline. CO2 off is the change
    percent_diff = ((a3_off_modelled_boxcox_streamflow - modelled_boxcox_streamflow) / a3_off_modelled_boxcox_streamflow) * 100
  ) |> 
  mutate(
    decade = year - (year %% 10)
  )


average_percent_diff_by_decade <- CO2_on_off_analysis |> 
  summarise(
    average_percent_diff = mean(percent_diff),
    median_percent_diff = median(percent_diff),
    .by = c(decade, gauge)
  ) |> 
  # add lat and lon for plotting
  left_join(
    lat_lon_gauge_info,
    by = join_by(gauge)
  ) |>
  filter(decade != 1950) |> # there is only a single year during this decade
  arrange(average_percent_diff) 

## Count gauges per decade
average_percent_diff_by_decade |> 
  summarise(
    n = n(),
    .by = decade
  )


average_percent_diff_by_decade |> 
  filter(decade == 1980) |> 
  ggplot(aes(x = lon, y = lat, colour = median_percent_diff)) +
  geom_point()



## Plot map percentage difference ==============================================
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
    mapping = aes(x = lon, y = lat, colour = average_percent_diff),
    data = average_percent_diff_by_decade,
    inherit.aes = FALSE,
    size = 0.75
  ) + 
  metR::scale_x_longitude(ticks = 10) +
  metR::scale_y_latitude(ticks = 10) +
  labs(
    x = "Longitude",
    y = "Latitude",
    colour = NULL
  ) +
  binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
    aesthetics = "colour",
    palette = big_palette,
    breaks = c(-50, -25, -5, -1, 0, 1, 5, 25, 50),
    limits = c(-5000, 800),
    show.limits = TRUE, 
    guide = "colorsteps"
    ) +
  facet_wrap(~decade) +
  labs(
    colour = "Mean Percentage Difference 
    (Modelled Streamflow CO2 Component Off - Modelled Streamflow CO2 Component On / 
    Modelled Streamflow CO2 Component Off)"
  ) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.7, 0.15),
    legend.title = element_text(size = 10)
  ) +
  guides(
    colour = guide_coloursteps(
      barwidth = unit(10, "cm"), 
      show.limits = TRUE, 
      even.steps = TRUE,
      title.position = "top",
      direction = "horizontal"
      )
    ) 

#plot_CO2_on_off_percent_diff

ggsave(
  filename = "./Graphs/CMAES_graphs/CO2_on_off.pdf",
  plot = plot_CO2_on_off_percent_diff,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)













stop_here <- tactical_typo()


# Compare non-CO2 vs. CO2 and a3-on vs a3-off ----------------------------------

# Would be a good idea to get streamflow_data_a3_off and streamflow_data_best_CO2_non_CO2
# in the same format

# suggested format?

gauge_key <- "A2390531" #"226410" #"403226"

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