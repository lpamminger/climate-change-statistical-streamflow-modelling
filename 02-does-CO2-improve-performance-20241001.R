# Map plots - work out if models with CO2 perform better than ones without
# Streamflow-time plots - what is the difference in streamflow?

cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, ozmaps, sf, patchwork)


# Import functions -------------------------------------------------------------
source("./Functions/utility.R")
source("./Functions/boxcox_transforms.R")


# Import data ------------------------------------------------------------------
data <- read_csv(
  "Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  show_col_types = FALSE
)

CMAES_results <- read_csv(
  "./Results/my_cmaes/CMAES_parameter_results_20241130.csv", 
  show_col_types = FALSE
  ) |> 
  filter(objective_function != "CO2_variable_objective_function") # temporary solution

data <- readr::read_csv(
  "./Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  col_types = "icdddddddll", # ensuring each column is a the correct type
  show_col_types = FALSE
)

streamflow_results <- read_csv(
  "Results/my_cmaes/CMAES_streamflow_results_20241130.csv",
  show_col_types = FALSE,
  col_select = !optimiser
)

gauge_information <- read_csv("./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE,
  col_select = c(gauge, bc_lambda, state, lat, lon)
) # CAMELS is Australia wide.



# Making map -------------------------------------------------------------------
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

aus_map[1,2] <- list(combine_NSW_ACT)
  
aus_map <- aus_map |> 
  filter(state != "ACT")

# Evidence ratio calculation ---------------------------------------------------
# What to do (Taken from Burnham and Anderson 2002):
## - assume non-CO2 is the best. best_CO2_AIC - best_non_CO2_AIC = AIC_diff. (This assumption does not impact the results. See chapter 2.11.2)
## - Positive values mean best_non_CO2 is the better model. Negative values mean the best_CO2 model is better.
## - calculate relative likelihood exp(0.5 * AIC_diff) pg. 72
## - evidence ratio is a relative measure of how much better the model is over the other model

best_CO2_and_non_CO2_per_catchment <- CMAES_results |>
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
      AIC_difference < 0 ~ exp(0.5 * abs(AIC_difference)), # when CO2 model is better
      AIC_difference > 0 ~ -exp(0.5 * abs(AIC_difference)) # when non-CO2 model is better
    )
  )









# Plot ready maps --------------------------------------------------------------
plot_ready_data <- evidence_ratio_calc |>
  select(!c(CO2, no_CO2, AIC_difference)) |>
  left_join(
    gauge_information,
    by = join_by(gauge)
  ) |> 
  arrange(evidence_ratio)


# I want the scale in these steps:
# From https://www.semanticscholar.org/paper/On-the-interpretation-of-likelihood-ratios-in-and-Martire-Kemp/3fc3690678409d8e8aa9352dce346565cf8fd0ea
# Weak or limited -> 1-10
# Moderate -> 10-100
# Moderately strong -> 100-1,000
# Strong -> 1,000-10,000
# Very strong -> 10,000-1,000,000
# Extremely strong ->  >1,000,000


## By state ====================================================================
evidence_ratio_by_state <- function(state, evidence_ratio_data, map_data) {
  
  # use state as key to extract correct polygon from map_data and evidence_ratio_data
  state_map <- map_data |> 
    filter(state == {{ state }})
  
  state_evidence_ratio <- evidence_ratio_data |> 
    filter(state == {{ state }})
  
  ggplot() +
  geom_sf(
    data = state_map,
    colour = "black",
    fill = "grey50"
  ) +
    geom_point(
      data = state_evidence_ratio,
      aes(x = lon, y = lat, colour = evidence_ratio),
      show.legend = FALSE,
      size = 0.5
    ) +
    binned_scale( # binned scale code taken from: https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "colour",
      palette = function(x) c("#f7f7f7", "#fddbc7", "#f4a582",  "#d6604d", "#b2182b", "#67001f"), # length should be length(breaks + limits) - 1
      breaks = c(1E1, 1E2, 1E3, 1E4, 1E6),
      limits = c(-1E1, 1E27),
      show.limits = TRUE,
      guide = "coloursteps"
    ) +
    theme_bw() +
    labs(
      x = "Longitude",
      y = "Latitude",
      colour = "Evidence Ratio"
    ) +
    theme(
      legend.key.height = unit(25, "mm"),
      legend.frame = element_rect(colour = "black"),
      axis.title = element_blank(),
      axis.text = element_text(size = 6)
    )
}

by_state_plots <- map(
  .x = aus_map |> pull(state) |> unique(),
  .f = evidence_ratio_by_state,
  evidence_ratio_data = plot_ready_data,
  map_data = aus_map
)

names(by_state_plots) <- aus_map |> pull(state) |> unique()

## By country ==================================================================
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
    size = 0.5
  ) +
  binned_scale( # binned scale code taken from: https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
    aesthetics = "colour",
    palette = function(x) c("#f7f7f7", "#fddbc7", "#f4a582",  "#d6604d", "#b2182b", "#67001f"), # length should be length(breaks + limits) - 1
    breaks = c(1E1, 1E2, 1E3, 1E4, 1E6),
    limits = c(-1E1, 1E27),
    show.limits = TRUE,
    guide = "coloursteps"
  ) +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude",
    colour = "Evidence Ratio"
  ) +
  theme(
    legend.key.height = unit(4, "mm"),
    legend.key.width = unit(40, "mm"),
    legend.frame = element_rect(colour = "black"),
    legend.position = "bottom"
  )


aus_evidence_ratio_map

## Combine to make a paper ready plot ==========================================
# make it look nice 
# reduce(.x = by_state_plots, .f = `+`)
#nice_plot <- (by_state_plots[[1]] | by_state_plots[[2]] | by_state_plots[[3]] | by_state_plots[[4]] | by_state_plots[[5]]) / 
 # (by_state_plots[[6]] + aus_evidence_ratio_map + by_state_plots[[7]]) + plot_layout(guides = "collect") + 
  #plot_annotation(tag_levels = "a")

top_nice_plot <- by_state_plots[["TAS"]] | by_state_plots[["QLD"]] | by_state_plots[["WA"]] | by_state_plots[["SA"]] | by_state_plots[["NT"]]
bottom_nice_plot <- by_state_plots[["VIC"]] | aus_evidence_ratio_map | by_state_plots[["NSW"]]
nice_plot <- top_nice_plot / bottom_nice_plot / guide_area() +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a")

nice_plot

ggsave(
  filename = paste0("nice_plot_", get_date(), ".pdf"),
  plot = nice_plot,
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


## Only include best streamflow that was calibrated on =========================
## join the included_in_calibration column
in_calibration <- data |> 
  select(year, gauge, included_in_calibration)

best_calibration_streamflow_results <-  best_streamflow_results |> 
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
  

chunk_streamflow_timeseries_plot <- function(data) {
  
  data |>
    ggplot(aes(x = year, y = streamflow, colour = streamflow_type)) +
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



## Examine the difference to the observed ======================================
### i.e., observed - best_CO2 and observed - best_non_CO2

# REMOVE (start)?

difference_to_observed_streamflow <- tidy_streamflow |>
  select(!c(bc_lambda, boxcox_streamflow)) |>
  distinct() |>
  pivot_wider(
    names_from = streamflow_type,
    values_from = streamflow,
  ) |>
  mutate(
    CO2_minus_non_CO2 = CO2 - non_CO2
    #observed_minus_CO2 = observed - CO2,
    #observed_minus_non_CO2 = observed - non_CO2
  )

# Totals and averages (not standardised meaning catchments with high rainfall will be larger)
totals_and_averages_streamflow <- difference_to_observed_streamflow |>
  summarise(
    n = n(),
    #sum_observed = sum(observed, na.rm = TRUE),
    #sum_CO2 = sum(CO2, na.rm = TRUE),
    #sum_non_CO2 = sum(non_CO2, na.rm = TRUE),
    sum_CO2_minus_non_CO2 = sum(CO2_minus_non_CO2, na.rm = TRUE),
    #sum_observed_minus_CO2 = sum(observed_minus_CO2, na.rm = TRUE),
    #sum_observed_minus_non_CO2 = sum(observed_minus_non_CO2, na.rm = TRUE),
    .by = gauge
  ) |>
  mutate(
    ave_CO2_minus_non_CO2 = sum_CO2_minus_non_CO2 / n
  ) |> 
  left_join(
    evidence_ratio_calc,
    by = join_by(gauge)
  ) |> 
  select(!c(no_CO2, CO2, AIC_difference)) |> #, #sum_CO2, sum_observed, sum_non_CO2))
  arrange(evidence_ratio)

### (end)



# How do the fitted parameters compared to the observed ones in RQ1? -----------
# RQ1 USED CAMELS-AUSv2 not v2. Expect to be slightly different.
RQ1_parameters <- c("a0", "a0_d", "a0_n", "a1", "a2", "sd")

RQ1_parameter_comparison <- CMAES_results |> 
  filter(parameter %in% RQ1_parameters)


stats_RQ1_parameter_comparison <- RQ1_parameter_comparison |> 
  summarise(
    median_parameters = median(parameter_value),
    q5_parameters = quantile(parameter_value, probs = 0.05),
    q95_parameters = quantile(parameter_value, probs = 0.95),
    .by = parameter
  ) |> 
  arrange(parameter)

RQ1_parameter_comparison |> 
  ggplot(aes(x = parameter_value)) +
  geom_histogram(
    colour = "black", 
    fill = "grey",
    bins = 30
    ) +
  geom_vline(
    aes(xintercept = median(parameter_value)), 
    colour = "red"
    ) +
  theme_bw() +
  facet_wrap(~parameter, scales = "free") 






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
  "streamflow_model_precip_only", "streamflow_model_separate_shifted_CO2",
  "streamflow_model_precip_auto", "streamflow_model_separate_shifted_CO2_auto",
  "streamflow_model_precip_seasonal_ratio", "streamflow_model_separate_shifted_CO2_seasonal_ratio",
  "streamflow_model_precip_seasonal_ratio_auto", "streamflow_model_separate_shifted_CO2_seasonal_ratio_auto",
  "streamflow_model_drought_precip_only", "streamflow_model_drought_separate_shifted_CO2",
  "streamflow_model_drought_precip_auto", "streamflow_model_drought_separate_shifted_CO2_auto",
  "streamflow_model_drought_precip_seasonal_ratio", "streamflow_model_drought_separate_shifted_CO2_seasonal_ratio",
  "streamflow_model_drought_precip_seasonal_ratio_auto", "streamflow_model_drought_separate_shifted_CO2_seasonal_ratio_auto"
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
# or if the a5 parameter is > 118.81 (this means the CO2 comp does nothing)

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
  select(!streamflow_model_objective_function) 

# Some catchments peform best when not all parameters are utlised. This is bad
utilisation_tolerence <- 1E-5

check_near_zero <- parameter_utilisation |>
  filter(abs(parameter_value) < utilisation_tolerence) |>  
  semi_join(
    best_model_combination_per_catchment,
    by = join_by(gauge, streamflow_model, objective_function)
  )

filter_check_near_zero <- check_near_zero |> 
  filter(parameter != "a5") 

# if sd is close to zero we have a problem


# ISSUE: Sometimes the model with CO2 not doing anything does "better" ---------
# the a model without the CO2 value. CHECK this.
# When selecting the best_model_combination_per_catchment we must 
# ensure each parameter is being utlised e.g., sd is being set to zero


assess_non_utilised_parmeters <- CMAES_results |> 
  filter(gauge %in% filter_check_near_zero$gauge) |> 
  anti_join(
    filter_check_near_zero,
    by = join_by(gauge, streamflow_model, objective_function) # remove non-utiliesd combiations
  ) |> 
  arrange(AIC)


# Using the assess_non_ultised_parameters find the next best model/obj combination
# that has all parameters utilised. Go down the list for each gauge.
# Approach: for each gauge, see if the model-objective function combination
#           has any parameter values that are extremely small. If so discard and
#           take the next best model-objective function combination.

check_unutilised_parameters <- assess_non_utilised_parmeters |> 
  mutate(
    is_not_utilised = parameter_value < utilisation_tolerence
  ) |> 
  summarise(
    sum_not_utilised = sum(is_not_utilised),
    AIC = max(AIC),
    .by = c(gauge, streamflow_model, objective_function)
  )

# REPLACE NOT UTLISED CATCHMENT-MODEL-OBJECTIVEFUNCTION COMBINATIONS -----------
# What does a very small sd or scale CO2 mean?
# a very low scale CO2 means -> variable_sd <- reshape_constant_sd + (reshape_scale_CO2 * reshape_CO2)
# the optimiser is setting something to zero
# I don't know why sd want to be really small for the constant objective functions? THIS DOES NOT OCCUR

# of the best_model_combination_per_catchment replace if parameters are not
# being utlised

new_best_model_combination_per_catchment <- best_model_combination_per_catchment |> 
  anti_join(
    check_unutilised_parameters,
    by = join_by(gauge, streamflow_model, objective_function)
  )

# If the size of new_best_model_combination_per_catchment and
# best_model_combination_per_catchment are the same them all parameters 
# in the best_model_combination_per_catchment are being utilised


# ignore from here





# Convert a5 into time of emergence --------------------------------------------
# find catchments with the best model including the a5 parameter (anything shifted)
# filter the a5 parameter out
# transform the a5 parameter into year
# make this compatible with tables
# if it is less than 39.58 then CO2 started impacting streamflow before 1959
# if the a5 is zero we cannot say it occur at industrialision because we have
# not fit the model to pre-industrialisation data.

# filter parameter_utilisation using the best_model_combination_per_catchment
just_a5_best_models <- parameter_utilisation |> 
  semi_join(
    best_model_combination_per_catchment,
    by = join_by(gauge, streamflow_model, objective_function)
  ) |> 
  filter(parameter == "a5") 



single_a5_to_year <- function(shifted_CO2_parameter, CO2, year) {
  adjusted_CO2 <- if_else(CO2 - shifted_CO2_parameter < 0, 0, CO2 - shifted_CO2_parameter)
  
  year_where_CO2_impacts_flow <- year[adjusted_CO2 != 0][1]
  
  return(year_where_CO2_impacts_flow)
}


map_single_a5_to_year <- function(gauge, just_a5_best_models, data) {
  year <- data |> 
    filter(gauge == {{ gauge }}) |> 
    pull(year)
  
  CO2 <- data |> 
    filter(gauge == {{ gauge }}) |> 
    pull(CO2)
  
  shift_CO2 <- just_a5_best_models |> 
    filter(gauge == {{ gauge }}) |> 
    pull(parameter_value)
  
  year_shift_occurred <- single_a5_to_year(shift_CO2, CO2, year)
  
  tibble(
    "gauge" = {{ gauge }},
    "a5" = shift_CO2,
    "year_shift_occurred" = year_shift_occurred
  )
}


time_of_emergence <- map(
  .x = just_a5_best_models$gauge,
  .f = map_single_a5_to_year,
  just_a5_best_models,
  data = data
) |> 
  list_rbind() 


more_data_time_of_emergence <-  time_of_emergence |> 
  left_join(
    gauge_information,
    by = join_by(gauge)
  ) |> 
  left_join(
    evidence_ratio_calc,
    by = join_by(gauge)
  ) 

# see when everything occurs  
more_data_time_of_emergence |> 
  ggplot(aes(x = gauge, y = year_shift_occurred, label = year_shift_occurred, colour = state)) +
  geom_point() +
  geom_text(nudge_y = 1) +
  theme_bw() +
  labs(
    x = "Gauge",
    y = "Year CO2 term kicked-in",
    colour = "State"
  ) +
  theme(
    axis.text.x = element_text(angle = 90)
  )


# is there a relationship between ToE and evi ratio
more_data_time_of_emergence |> 
  ggplot(aes(x = abs(evidence_ratio), y = year_shift_occurred)) +
  geom_smooth(
    method = "lm", 
    formula = y ~ x, 
    se = TRUE, 
    colour = "black",
    linewidth = 0.5
  ) +
  geom_point(
    aes(colour = state),
    size = 2
    ) +
  theme_bw() +
  scale_x_log10() +
  labs(
    x = "Evidence Ratio",
    y = "Year CO2 term kicked-in",
    colour = "State"
  )










# I am not sure if the code below here works...

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

#mb <- as.numeric(1:10 %o% 10 ^ (0:3)) # this would be useful to know earlier

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





