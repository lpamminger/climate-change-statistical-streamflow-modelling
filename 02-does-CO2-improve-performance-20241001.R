# Map plots - work out if models with CO2 perform better than ones without
# Streamflow-time plots - what is the difference in streamflow?

cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, sf)


# Import functions -------------------------------------------------------------
source("./Functions/utility.R")
source("./Functions/boxcox_transforms.R")


# Import data ------------------------------------------------------------------
CMAES_results <- read_csv("./Results/my_cmaes/CMAES_parameter_results_20241108.csv", show_col_types = FALSE)

data <- readr::read_csv(
  "./Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  col_types = "icdddddddll", # ensuring each column is a the correct type
  show_col_types = FALSE
)

streamflow_results <- read_csv("Results/my_cmaes/CMAES_streamflow_results_20241108.csv",
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


## Remove non-consectutive years of observed and modelled streamflow ===========
chopped_streamflow_results <- best_streamflow_results |> 
  drop_na() |> 
  mutate(
    lag_year = dplyr::lag(year, n = 1L),
    .by = c(gauge, streamflow_model, objective_function),
    .after = 1
  ) |> 
  mutate(
    diff_year = year - lag_year,
    .after = 2
  ) |> 
  mutate(
    diff_year = if_else(is.na(diff_year), 1, diff_year)
  ) |> 
  filter(
    diff_year == 1 # only have consecutive observed years in data
  )



## Summarise results into a tidy format ========================================
tidy_boxcox_streamflow <- chopped_streamflow_results |>
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
plot_streamflow_timeseries <- tidy_streamflow |>
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


ggsave(
  filename = paste0("streamflow_timeseries_comparison_", get_date(), ".pdf"),
  plot = plot_streamflow_timeseries,
  device = "pdf",
  path = "./Graphs/CMAES_graphs",
  width = 1189,
  height = 841,
  units = "mm"
)



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
    standardised_CO2_minus_non_CO2 = (CO2 - non_CO2) / precipitation#,
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
    sum_standardised_CO2_minus_non_CO2 = sum(standardised_CO2_minus_non_CO2, na.rm = TRUE),
    #sum_observed_minus_CO2 = sum(observed_minus_CO2, na.rm = TRUE),
    #sum_observed_minus_non_CO2 = sum(observed_minus_non_CO2, na.rm = TRUE),
    .by = gauge
  ) |>
  mutate(
    ave_CO2_minus_non_CO2 = sum_CO2_minus_non_CO2 / n,
    ave_standardised_CO2_minus_non_CO2 = sum_standardised_CO2_minus_non_CO2 / n
  ) |> 
  arrange(desc(ave_standardised_CO2_minus_non_CO2)) |> 
  left_join(
    evidence_ratio_calc,
    by = join_by(gauge)
  ) |> 
  select(!c(no_CO2, CO2, AIC_difference))#, #sum_CO2, sum_observed, sum_non_CO2))






# How do the fitted parameters compared to the observed ones in RQ1? -----------
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

# Some catchments peform best when not all parameters are utilised. This is bad
check_near_zero <- parameter_utilisation |>
  filter(abs(parameter_value) < 1E-6) |>  # sd and scale_CO2 are okay
  semi_join(
    best_model_combination_per_catchment,
    by = join_by(gauge, streamflow_model, objective_function)
  )


# Get next best model for non-utlised parameter model objective function combs.
# sd or CO2 has to have some impact. Deemed impactful if I can see the numbers in console i.e., 7.3E-5 + (30 * 3.9E-2) vs. 7.3E-8 + (30 * 3.9E-2)

x <- CMAES_results |> 
  filter(gauge %in% check_near_zero$gauge) |> 
  anti_join(
    check_near_zero,
    by = join_by(gauge, streamflow_model, objective_function) # remove non-utiliesd combiations
  ) |> 
  #filter(gauge == "A5130501") |> 
  arrange(AIC)


# Not okay = 137101A REPLACE streamflow_model_separate_shifted_CO2_seasonal_ratio CO2_variable_objective_function
#          = 215207 REPLACE streamflow_model_drought_precip_seasonal_ratio_auto constant_sd_objective_function
#          = 216004 REPLACE streamflow_model_separate_shifted_CO2_auto CO2_variable_objective_function
#          = 415237 REPLACE streamflow_model_precip_seasonal_ratio_auto CO2_variable_objective_function
#          = 138004B streamflow_model_precip_only CO2_variable_objective_function
#          = A5130501 streamflow_model_precip_seasonal_ratio_auto constant_sd_objective_function


# does near zero scale CO2 and sd catchments always perform worse? No. Sometimes they are better.
AIC_check_near_zero <- evidence_ratio_calc |> 
  filter(gauge %in% check_near_zero$gauge) |> 
  filter(AIC_difference < 0)


check_a5 <- parameter_utilisation |> 
  filter(parameter == "a5") |> 
  filter(parameter_value > 118.81) # good - nothing is greater than this, meaning the AIC successfully ignored the catchments

# What does a very small sd or scale CO2 mean?
# a very low scale CO2 means -> variable_sd <- reshape_constant_sd + (reshape_scale_CO2 * reshape_CO2)
# the optimiser is setting something to zero
# I don't know why sd want to be really small for the constant objective functions? THIS DOES NOT OCCUR




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





# Moved from 03-DREAM ----------------------------------------------------------
# DREAM - t.thin
# do I want everything to have the same number of sequences?
# or do I want it to scale with parameter number?



# analysising dream results WILL MOVE TO ANOTHER FILE --------------------------
# how close are the test catchments parameter to CMAES
# histogram
# comparison of streamflow time graphs?
sequences <- read_csv("Results/my_dream/DREAM_sequence_results.csv", show_col_types = FALSE)


### Convert a5 parameter into year #############################################
CO2 <- data |>
  filter(gauge == "102101A") |>
  pull(CO2)

year <- data |>
  filter(gauge == "102101A") |>
  pull(year)




# make this compatible with tables
single_a5_to_year <- function(shifted_CO2_parameter, CO2, year) {
  adjusted_CO2 <- if_else(CO2 - shifted_CO2_parameter < 0, 0, CO2 - shifted_CO2_parameter)
  
  year_where_CO2_impacts_flow <- year[adjusted_CO2 != 0][1]
  
  return(year_where_CO2_impacts_flow)
}


shifted_CO2_parameter_into_year_CO2_starts_impacting_flow <- function(shifted_CO2_parameter, CO2, year) {
  map_dbl(
    .x = shifted_CO2_parameter,
    .f = single_a5_to_year,
    CO2 = CO2,
    year = year
  )
}


just_a5 <- sequences |>
  filter(parameter == "a5") |>
  mutate(
    year_where_CO2_impacts_flow = shifted_CO2_parameter_into_year_CO2_starts_impacting_flow(
      shifted_CO2_parameter = parameter_values,
      CO2 = CO2,
      year = year
    )
  )


a5_to_year_plot <- just_a5 |>
  ggplot(aes(x = year_where_CO2_impacts_flow)) +
  geom_histogram(
    binwidth = binwidth_bins(30),
    fill = "grey",
    colour = "black"
  ) +
  # scale_y_sqrt() +
  labs(
    x = "Year where CO2 starts impacting streamflow",
    y = "Frequency"
  ) +
  theme_bw() +
  facet_wrap(~gauge, scales = "free")


a5_to_year_plot




# See if thin.t produces 1000 different parameter combinations -----------------
count_combinations <- sequences |>
  summarise(
    n = n(),
    .by = c(gauge, streamflow_model, objective_function, parameter)
  )

# I don't know how thining works. Also if DREAM converges before
# the max ndraw is meet it will not be 1000 different combinations



# Plot the histograms of parameter values ======================================
plot_parameter_histogram <- function(streamflow_model, objective_function, sequence_results) {
  sequence_results |>
    filter(streamflow_model == {{ streamflow_model }}) |>
    filter(objective_function == {{ objective_function }}) |>
    ggplot(aes(x = parameter_values)) +
    geom_histogram(
      binwidth = binwidth_bins(30),
      fill = "grey",
      colour = "black"
    ) +
    theme_bw() +
    scale_y_sqrt() +
    labs(
      x = "Range of Parameter Values",
      y = "Frequency",
      title = paste0("Streamflow Model: ", streamflow_model, "\nObjective Function: ", objective_function)
    ) +
    facet_grid(gauge ~ parameter, scales = "free")
}


models <- unique(pull(sequences, streamflow_model))
objfun <- unique(pull(sequences, objective_function))

histogram_plot_1 <- plot_parameter_histogram(
  streamflow_model = models[4],
  objective_function =  objfun,
  sequence_results = sequences
)

histogram_plot_2 <- plot_parameter_histogram(
  streamflow_model = models[2],
  objective_function =  objfun,
  sequence_results = sequences
)

histogram_plot_1
histogram_plot_2



# Compare streamflow timeseries ================================================
DREAM_parameter_results <- read_csv("Results/my_dream/DREAM_parameter_results.csv", show_col_types = FALSE, col_types = "cccccdddcl") |>
  drop_na()
# these will have to be inputted into the respective models
# Temporary and hardcoded

CMAES_streamflow <- read_csv("Results/my_cmaes/CMAES_streamflow_results_20241028.csv", show_col_types = FALSE)


DREAM_join <- DREAM_parameter_results |>
  select(c(gauge, streamflow_model, objective_function)) |>
  distinct()


CMAES_filtered_streamflow <- CMAES_streamflow |>
  semi_join(DREAM_join, by = join_by(gauge, streamflow_model, objective_function))


parameters <- DREAM_parameter_results |>
  filter(gauge == "408202") |>
  filter(streamflow_model == "streamflow_model_drought_separate_shifted_CO2_seasonal_ratio_auto") |>
  pull(parameter_value)


gauge <- "408202"

catchment_data <- gauge |>
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  )

catchment_data$stop_start_data_set <- NULL

catchment_data_altered <- unlist(catchment_data, recursive = FALSE)

names(catchment_data_altered) <- c(
  "gauge_ID",
  "contains_drought",
  "year",
  "precipitation",
  "observed_boxcox_streamflow",
  "is_drought_year",
  "CO2",
  "seasonal_ratio"
)

streamflow <- streamflow_model_drought_separate_shifted_CO2_seasonal_ratio_auto(
  catchment_data = catchment_data,
  parameter_set = parameters
)




# Compare parameter results ====================================================
CMAES_parameter_results <- read_csv("Results/my_cmaes/CMAES_parameter_results_20241028.csv", show_col_types = FALSE)

compare_parameter_results <- CMAES_parameter_results |>
  right_join(DREAM_parameter_results, by = join_by(gauge, streamflow_model, objective_function)) |>
  mutate(
    AIC_diff = AIC.x - AIC.y
  )



# TESTING ----------------------------------------------------------------------
# make sure everything gets a good exit message

test_suite <- best_model_combination_per_catchment |> 
  distinct(
    streamflow_model_name,
    objective_function_name,
    .keep_all = TRUE
  )

# hard code into the test below:
# - works for gauge = 105102A, model = streamflow_model_drought_precip_only, obj = constant_sd_objective_function

tic <- as.numeric(Sys.time())

gauge <- "146095A"

dream_example <- gauge |>
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  ) |>
  numerical_optimiser_setup_vary_inputs(
    streamflow_model = streamflow_model_precip_only,
    objective_function = constant_sd_objective_function,
    bounds_and_transform_method = make_default_bounds_and_transform_methods(),
    minimise_likelihood = FALSE
  ) |>
  my_dream(print_monitor = TRUE) |>
  result_set()

toc <- round(as.numeric(Sys.time()) - tic, 1)

cat(toc, "sec")

dream_parameters <- dream_example |>
  parameters_summary()

dream_example |>
  plot()






test_sequences |>
  pivot_longer(
    cols = everything(),
    names_to = "parameter",
    values_to = "parameter_values"
  ) |>
  ggplot(aes(x = parameter_values)) +
  geom_histogram(
    binwidth = binwidth_bins(30),
    fill = "grey",
    colour = "black"
  ) +
  theme_bw() +
  scale_y_sqrt() +
  labs(
    x = "Range of Parameter Values",
    y = "Frequency"
  ) +
  facet_wrap(~parameter, scales = "free")



# Current seq
PARAMETER_NUMBER <- 9
seq <- round_any(((PARAMETER_NUMBER - 2)^2.5) * 1E4, 1E4, ceiling)
total_save_seq <- 1000
thin <- seq / total_save_seq




# for gauge = 407214, model = streamflow_model_drought_separate_CO2_seasonal_ratio_auto
# and objective function CO2_variable_objective_function expect LL of 95
# 550,000 nseq and 425 sec (~7min) gives 105 LL
# I could bump it more. Make a stepper curve
# New equation used -> round_any(((x - 2) ^ 2.5) * 1E4, 1E4, ceiling)
# Untested for 8 parameter models.

# N_PARAMETERS = 3 -> nseq = 50,000
# Max parameters are 8. Extrapolate? if linear it would suggest round_any((50000/3), 1E4, ceiling) = 20000 per parameter
# I don't this this will work. I think the relationship between nseq and number of parameters is non-linear
# Want to get the likelihoods to the nearest whole number of DREAM and CMAES

# For gauge 407214 maximum parameters DREAM LL = 99.xxx, CMAES LL = 95.xxx
# need to increase from 20000 * PARAMETER_NUMBER to or make it a non-linear relationship 160,000 for 8 params its too little. Try 200,000? Not enough.
# Try 300,000? Not enough.
# Playing with non-linear relationship for nseq and n_parameters:
x <- seq(from = 3, to = 8, by = 1)
y <- ((x^2) - (x + 1)) * 1E4
plot(x, y)

yy <- round_any(((x - 2)^2.5) * 1E4, 1E4, ceiling)
yy
# non-linear relationship works with gauge 407214 and 3, 4 parameters
# test for 5, 6, 7, 8? Go straight to 8? Good idea to test all.



# Checking CMAES - can delete/move
# a3 = -25 works (tick) 108003A can have an a3 less than -10 streamflow_model_separate_shifted_CO2	constant_sd_objective_function
# a3 = 50 works. 112102A	a3 > 10 streamflow_model_separate_shifted_CO2_seasonal_ratio_auto	constant_sd_objective_function
# a3 = 20 works. G8150018 a3>10	streamflow_model_separate_shifted_CO2_seasonal_ratio_auto	constant_sd_objective_function


# Results
CMAES_parameter_results <- read_csv(
  "Results/my_cmaes/CMAES_parameter_results_20241101.csv",
  show_col_types = FALSE
)


inspect_parameter_a3 <- CMAES_parameter_results |> 
  filter(parameter == "a3") |> 
  pull(parameter_value) |> 
  range()

# Changes made based on parameter results
# - catchments still hitting a3 bounds - make slightly larger
# - scale_CO2 want to be very close to zero. Cannot be zero due to log transform
# - try to reduce bounds to speed up time. Bounds to reduce:
# a0_d/a0_n, a4, sd

inspect_parameter <- CMAES_parameter_results |> # is the CO2 being shut-off? Sometimes. Hopefully the AIC deems these models as not good.
  filter(parameter == "a5") |> 
  filter(parameter_value > 110)

# Testing
random_test <- CMAES_parameter_results |> 
  filter(streamflow_model == streamflow_model_drought_separate_shifted_CO2_seasonal_ratio_auto()$name) |> 
  slice_sample(
    n = 1
  )

gauge <- "227225A"
example <- gauge |>
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  ) |>
  numerical_optimiser_setup_vary_inputs(
    streamflow_model = streamflow_model_drought_separate_shifted_CO2_seasonal_ratio_auto,
    objective_function = constant_sd_objective_function,
    bounds_and_transform_method = make_default_bounds_and_transform_methods(),
    minimise_likelihood = TRUE
  ) |>
  my_cmaes(print_monitor = TRUE) |>
  result_set() 

x <- example |> parameters_summary()
example |> plot()
# based on trial and error set a3 to -25 to 50
