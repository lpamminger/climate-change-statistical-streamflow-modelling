# Map plots - work out if models with CO2 perform better than ones without
# Streamflow-time plots - what is the difference in streamflow?

cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, sf)


# Import functions -------------------------------------------------------------
source("./Functions/utility.R")
source("./Functions/boxcox_transforms.R")


# Import data ------------------------------------------------------------------
CMAES_results <- read_csv("./Results/my_cmaes/CMAES_parameter_results_20241031.csv", show_col_types = FALSE)


streamflow_results <- read_csv("Results/my_cmaes/CMAES_streamflow_results_20241031.csv",
  show_col_types = FALSE,
  col_select = !optimiser
)

gauge_information <- read_csv("./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE,
  col_select = c(gauge, bc_lambda, state, lat, lon)
) # CAMELS is Australia wide.



aus_map <- read_sf(
  dsn = "./Data/Maps",
  layer = "AUS_2016_AUST"
)

vic_map <- read_sf(
  dsn = "./Data/Maps",
  layer = "vic_map_simple"
)


# Gauges that do weird things - remove for now ---------------------------------
## Gauges visually inspected using check_model_fit in graphs
#removed_gauges <- c("G0060005", "A2390531")



# Map plots --------------------------------------------------------------------
# For a given catchment work out if a model with CO2 outperforms a model without CO2
# remove unnecessary columns
# only have unique streamflow model, gauge combinations
# combine the streamflow model and objective function using unite
# regex to see if the united columns contain "CO2"
# if yes then contains_CO2 = TRUE otherwise contains_CO2 is false
# for each catchment (gauge) get the two lowest AIC values by contains_CO2

best_CO2_and_non_CO2_per_catchment <- CMAES_results |>
  #filter(!gauge %in% removed_gauges) |>
  select(!c(parameter, parameter_value, optimiser, loglikelihood, exit_message, near_bounds)) |>
  unite(
    col = streamflow_model_objective_function,
    c(streamflow_model, objective_function),
    sep = "_",
    remove = FALSE,
    na.rm = FALSE
  ) |>
  mutate(
    contains_CO2 = str_detect(streamflow_model_objective_function, "CO2"),
    contains_CO2 = if_else(contains_CO2, "CO2", "no_CO2"),
    .after = 2
  ) |>
  slice_min(
    AIC,
    by = c(gauge, contains_CO2)
  ) |> 
  distinct()


# This will be used for selecting models and gauges for DREAM
best_model_combination_per_catchment <- best_CO2_and_non_CO2_per_catchment |> 
  slice_min(
    AIC,
    by = gauge
  )




# Evidence ratio calculation ---------------------------------------------------
# What to do (Taken from Burnham and Anderson 2002):
## - assume non-CO2 is the best. best_CO2_AIC - best_non_CO2_AIC = AIC_diff. (This assumption does not impact the results. See chapter 2.11.2)
## - Positive values mean best_non_CO2 is the better model. Negative values mean the best_CO2 model is better.
## - calculate relative likelihood exp(0.5 * AIC_diff) pg. 72
## - evidence ratio is a relative measure of how much better the model is over the other model


evidence_ratio_calc <- best_CO2_and_non_CO2_per_catchment |>
  select(!c(streamflow_model_objective_function, streamflow_model, objective_function)) |>
  pivot_wider(
    names_from = contains_CO2,
    values_from = AIC
  ) |>
  mutate(
    AIC_difference = CO2 - no_CO2 # CO2 is smaller than non-CO2 then negative and CO2 is better
  ) |>
  mutate(
    evidence_ratio = case_when(
      AIC_difference < 0 ~ -exp(0.5 * abs(AIC_difference)), # when CO2 model is better
      AIC_difference > 0 ~ exp(0.5 * abs(AIC_difference)) # when non-CO2 model is better
    )
  )





# Plot time --------------------------------------------------------------------
plot_ready_data <- evidence_ratio_calc |>
  select(!c(CO2, no_CO2, AIC_difference)) |>
  left_join(
    gauge_information,
    by = join_by(gauge)
  )


# I want the scale in these steps:
# From https://www.semanticscholar.org/paper/On-the-interpretation-of-likelihood-ratios-in-and-Martire-Kemp/3fc3690678409d8e8aa9352dce346565cf8fd0ea
# Weak or limited -> 1-10
# Moderate -> 10-100
# Moderately strong -> 100-1,000
# Strong -> 1,000-10,000
# Very strong -> 10,000-1,000,000
# Extremely strong ->  >1,000,000



## Victorian plot ==============================================================
vic_evidence_ratio_map <- ggplot() +
  geom_sf(
    data = vic_map,
    colour = "black",
    fill = "grey50"
  ) +
  geom_point(
    data = plot_ready_data |> filter(state == "VIC"),
    aes(x = lon, y = lat, colour = evidence_ratio)
  ) +
  binned_scale( # binned scale code taken from: https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
    aesthetics = "colour",
    palette = function(x) c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#f7f7f7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac", "#053061"), # length should be length(breaks + limits) - 1
    breaks = c(-1E6, -1E4, -1E3, -1E2, -1E1, 1E1, 1E2, 1E3, 1E4, 1E6),
    limits = c(-1E13, 1E13),
    show.limits = FALSE,
    guide = "coloursteps"
  ) +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude",
    colour = "Evidence Ratio"
  ) +
  theme(
    legend.key.height = unit(25, "mm")
  )

ggsave(
  filename = paste0("vic_evidence_ratio_map_", get_date(), ".pdf"),
  plot = vic_evidence_ratio_map,
  device = "pdf",
  path = "./Graphs/CMAES_graphs",
  height = 210,
  width = 297,
  units = "mm"
)

# Australian plot ==============================================================
aus_evidence_ratio_map <- ggplot() +
  geom_sf(
    data = aus_map,
    colour = "black",
    fill = "grey50"
  ) +
  coord_sf(xlim = c(110, 155)) +
  geom_point(
    data = plot_ready_data,
    aes(x = lon, y = lat, colour = evidence_ratio),
    size = 1
  ) +
  binned_scale( # binned scale code taken from: https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
    aesthetics = "colour",
    palette = function(x) c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#f7f7f7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac", "#053061"), # length should be length(breaks + limits) - 1
    breaks = c(-1E6, -1E4, -1E3, -1E2, -1E1, 1E1, 1E2, 1E3, 1E4, 1E6),
    limits = c(-1E13, 1E13),
    show.limits = FALSE,
    guide = "coloursteps"
  ) +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude",
    colour = "Evidence Ratio"
  ) +
  theme(
    legend.key.height = unit(25, "mm")
  )

ggsave(
  filename = paste0("aus_evidence_ratio_map_", get_date(), ".pdf"),
  plot = aus_evidence_ratio_map,
  device = "pdf",
  path = "./Graphs/CMAES_graphs",
  height = 210,
  width = 297,
  units = "mm"
)



# Streamflow-timeseries plots --------------------------------------------------
# the objective is to plot the best CO2, best non-CO2 and observed streamflow
# (non box-cox) on a timeseries graph and examine the difference


## Filter based on best non-CO2 and CO2 models =================================
best_streamflow_results <- streamflow_results |>
  semi_join(
    best_CO2_and_non_CO2_per_catchment,
    by = join_by(gauge, streamflow_model, objective_function)
  )



## Summarise results into a tidy format ========================================
tidy_boxcox_streamflow <- best_streamflow_results |>
  pivot_longer(
    cols = c(observed_boxcox_streamflow, modelled_boxcox_streamflow),
    names_to = "name",
    values_to = "boxcox_streamflow"
  ) |>
  unite(
    col = streamflow_model_objective_function,
    c(streamflow_model, objective_function),
    sep = "_",
    remove = TRUE,
    na.rm = FALSE
  ) |>
  mutate(
    streamflow_model_objective_function = if_else(name == "observed_boxcox_streamflow", "observed", streamflow_model_objective_function)
  ) |>
  select(!name) |>
  mutate(
    streamflow_type = case_when(
      str_detect(streamflow_model_objective_function, "CO2") & !str_detect(streamflow_model_objective_function, "observed") ~ "CO2",
      !str_detect(streamflow_model_objective_function, "CO2") & !str_detect(streamflow_model_objective_function, "observed") ~ "non_CO2",
      .default = "observed"
    )
  ) |>
  select(!c(streamflow_model_objective_function))



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
plot_streamflow_timeseries <- tidy_streamflow |>
  ggplot(aes(x = year, y = streamflow, colour = streamflow_type)) +
  geom_line(na.rm = TRUE, alpha = 0.5) +
  theme_bw() +
  scale_colour_brewer(palette = "Set1") +
  labs(
    x = "Year",
    y = "Streamflow (mm)"
  ) +
  facet_wrap(~gauge, scales = "free_y") +
  theme(legend.title = element_blank())


ggsave(
  filename = paste0("streamflow_timeseries_comparison_", get_date(), ".pdf"),
  plot = plot_streamflow_timeseries,
  device = "pdf",
  path = "./Graphs/CMAES_graphs",
  width = 1189,
  height = 841,
  units = "mm"
)



# I am not sure if the code from line 302 below works...



## Examine the difference to the observed ======================================
### i.e., observed - best_CO2 and observed - best_non_CO2

difference_to_observed_streamflow <- tidy_streamflow |>
  select(!c(bc_lambda, boxcox_streamflow)) |>
  distinct() |>
  pivot_wider(
    names_from = streamflow_type,
    values_from = streamflow,
  ) |>
  mutate(
    CO2_minus_non_CO2 = CO2 - non_CO2,
    observed_minus_CO2 = observed - CO2,
    observed_minus_non_CO2 = observed - non_CO2
  )



# If you purely compare the model with CO2 and without CO2 than...
difference_to_observed_streamflow |>
  # filter(gauge == "221207") |>
  ggplot(aes(x = year, y = CO2_minus_non_CO2)) +
  geom_line(na.rm = TRUE) +
  theme_bw() +
  facet_wrap(~gauge, scales = "free_y")

# Produces linear changes in CO2_minus_non_CO2 overtime
# What does this mean?
# A positive linear slope means as timeseries progresses the CO2 streamflow increases and non_CO2 decreases (opposite is true)



# If you compare the model (both CO2 and non-CO2) to observed than...
difference_to_observed_streamflow |>
  pivot_longer(
    cols = starts_with("observed_minus"),
    names_to = "residual_type",
    values_to = "streamflow_residual"
  ) |>
  filter(gauge == "403217") |>
  ggplot(aes(x = year, y = streamflow_residual, colour = residual_type)) +
  geom_line(na.rm = TRUE) +
  theme_bw() #+
# facet_wrap(~gauge, scales = "free_y")

x <- difference_to_observed_streamflow |>
  drop_na() |>
  filter(gauge == "403217") |>
  pull(observed_minus_non_CO2) |>
  acf()

# Totals and averages
totals_and_averages_streamflow <- difference_to_observed_streamflow |>
  summarise(
    n = n(),
    sum_observed = sum(observed, na.rm = TRUE),
    sum_CO2 = sum(CO2, na.rm = TRUE),
    sum_non_CO2 = sum(non_CO2, na.rm = TRUE),
    sum_CO2_minus_non_CO2 = sum(CO2_minus_non_CO2, na.rm = TRUE),
    sum_observed_minus_CO2 = sum(observed_minus_CO2, na.rm = TRUE),
    sum_observed_minus_non_CO2 = sum(observed_minus_non_CO2, na.rm = TRUE),
    .by = gauge
  ) |>
  mutate(
    ave_CO2_minus_non_CO2 = sum_CO2_minus_non_CO2 / n
  )


# Directly comparing model results
totals_and_averages_streamflow |>
  ggplot(aes(x = gauge, y = ave_CO2_minus_non_CO2)) +
  geom_col() +
  theme_bw() +
  labs(
    x = "Gauge",
    y = "Average annual streamflow difference between CO2 and non-CO2 models (mm/year)"
  ) +
  theme(
    axis.text.x = element_text(angle = 90),
    panel.grid.minor = element_blank()
  )





# Explore parameter combinations -----------------------------------------------
# Questions to answer:

# 1. Are the no CO2 and CO2 models selected for each catchment the same models with CO2?
# I must consider two cases:
# Case 1: when the streamflow model are the same but the objective functions are different
# Case 2: when the objective functions are the same and the streamflow models are equivalent except one has CO2 and the other does not

## Case 1 ======================================================================
# if the streamflow model is the same for the CO2 and non-CO2 and the objective function is different
# This must check if the unique streamflow models contain CO2? NO. Both streamflow models cannot have CO2

compare_objective_function <- function(streamflow_models, objective_functions) {
  if ((length(unique(streamflow_models)) == 1) & (length(unique(objective_functions)) != 1)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

same_streamflow_model_different_objective_function <- best_CO2_and_non_CO2_per_catchment |>
  summarise(
    same_model_different_objective_function = compare_objective_function(streamflow_model, objective_function),
    .by = gauge
  ) |>
  filter(same_model_different_objective_function) |>
  pull(gauge)



## Case 2 ======================================================================
# I think the best way it to write out CO2 and non_CO2 model pairs (equivalent models) in a tibble then check against CMAES_results
# This is hard coded. Would be better to use get_non_drought_streamflow_models() and get_drought_streamflow_models()
non_CO2_streamflow_model_vs_CO2_equivalent <- tribble(
  ~non_CO2_model, ~CO2_equivalent,
  "streamflow_model_precip_only", "streamflow_model_separate_CO2",
  "streamflow_model_precip_auto", "streamflow_model_separate_CO2_auto",
  "streamflow_model_precip_seasonal_ratio", "streamflow_model_separate_CO2_seasonal_ratio",
  "streamflow_model_precip_seasonal_ratio_auto", "streamflow_model_separate_CO2_seasonal_ratio_auto",
  "streamflow_model_drought_precip_only", "streamflow_model_drought_separate_CO2",
  "streamflow_model_drought_precip_auto", "streamflow_model_drought_separate_CO2_auto",
  "streamflow_model_drought_precip_seasonal_ratio", "streamflow_model_drought_separate_CO2_seasonal_ratio",
  "streamflow_model_drought_precip_seasonal_ratio_auto", "streamflow_model_drought_separate_CO2_seasonal_ratio_auto"
)



# if it contains both columns of non_CO2_streamflow_model_vs_CO2_equivalent for a given row then only difference is CO2 component in model
compare_streamflow_models <- function(streamflow_models, objective_functions, lookup_tibble) {
  filtered_lookup_tibble <- lookup_tibble |>
    filter(non_CO2_model == streamflow_models[1]) |> # This assumes the non_CO2 model is always first. I am not sure this is always true.
    unlist() |>
    unname()

  return(all(filtered_lookup_tibble == streamflow_models) & (length(unique(objective_functions)) == 1)) # ensures objective functions are the same
}


same_model_except_with_CO2_check <- best_CO2_and_non_CO2_per_catchment |>
  summarise(
    same_model_with_and_without_CO2 = compare_streamflow_models(
      streamflow_models = streamflow_model,
      objective_functions = objective_function,
      lookup_tibble = non_CO2_streamflow_model_vs_CO2_equivalent
    ),
    .by = gauge
  ) |>
  filter(same_model_with_and_without_CO2) |>
  pull(gauge)


# Answering the first question =================================================\
cat(
  length(same_model_except_with_CO2_check),
  "out of",
  length(unique(best_CO2_and_non_CO2_per_catchment$gauge)),
  "catchments have the equivalent CO2 and non CO2 models (i.e., same model except one has CO2)" # they also have the same objective function
)

cat(
  length(same_streamflow_model_different_objective_function),
  "out of",
  length(unique(best_CO2_and_non_CO2_per_catchment$gauge)),
  "catchments have the equivalent CO2 and non CO2 objective functions with the same streamflow model"
)



# 2. Are all parameters being utilised for the CO2 models? i.e., the model doesn’t set one value to zero. Check all values around zero?

# make sure same_streamflow_model_different_objective_function parameters are being fully utilised
# i.e., no parameters are really small

parameter_utilisation <- CMAES_results |>
  filter(gauge %in% c(same_streamflow_model_different_objective_function, same_model_except_with_CO2_check)) |>
  select(!c(optimiser, loglikelihood, exit_message, near_bounds)) |>
  distinct() |>
  unite(
    col = streamflow_model_objective_function,
    c(streamflow_model, objective_function),
    sep = "_",
    remove = FALSE,
    na.rm = FALSE
  ) |>
  mutate(
    contains_CO2 = str_detect(streamflow_model_objective_function, "CO2"),
    contains_CO2 = if_else(contains_CO2, "CO2", "no_CO2"),
    .after = 2
  ) |>
  slice_min(
    AIC,
    by = c(gauge, contains_CO2)
  ) |>
  select(!streamflow_model_objective_function) |>
  filter(abs(parameter_value) < 1E-4)

# What does a very small sd or scale CO2 mean?
# a very low scale CO2 means -> variable_sd <- reshape_constant_sd + (reshape_scale_CO2 * reshape_CO2)
# the optimiser is setting something to zero
# I don't know why sd want to be really small for the constant objective functions? THIS DOES NOT OCCUR
