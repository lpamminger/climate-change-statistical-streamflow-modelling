# Q. Should I limit the a5 (maximum shift CO2) bounds? 
##
## - I originally allowed the a5 parameter to be "turned-off" to help identify
##   whether the CO2 model or non-CO2 model is better. I assumed a "turned-off"
##   model with perform worse than a model with all parameters fully utilised.
##
## - Using the 20241130 results (a5 can exceed max calibrated CO2), there are no
##   catchments where the turned-off a5 value model outpeforms the non-CO2 model.
##   This tells me the AIC can differentiate between turned-off parameters in models. 
##
## - The problem with allowing the a5 parameter to be "turned-off" is 
##   it produces a complex response surface. This makes it difficult for 
##   DREAM to build an adequate parameter distribution.
##
## - I could make DREAM and CMAES have different bounds. This feels wrong 
##   because I am using the best CMAES models an putting them into DREAM. If
##   I change the bounds the best CMAES models may not be the best anymore.
## 



# Q. How do the limited and non-limited a5 bounds compare?
# A. Catchment 112101B as the "near_bounds" flag for both the limited and non-limited approach  

# Options:
# 1. Leave it. DREAM will show the distributions.
# 2. Swap the best-CO2 and best-non-CO2 for that given gauge for all
# 3. Swap the best-CO2 and best-non-CO2 for a given gauge if delta_AIC < 5 (exp(-0.5 * -5) = ~ 10 starting to be meaning after this)
# 4. Try the next best CO2 model, see if it kick-in on the last year (very likely a5 kicks in around the same time) - my guess was correct

# Complete option 4. My guess was correct
# I will do option 3. Changes required when graphing best (maps and timeseries) and DREAM


# Q. How does the best_CO2_non_CO2 models change using the limiting vs. non-limiting?
#    If there are changes this supports my idea that DREAM and CMAES should not
#    have different bounds.
# A. Changing the bounds changes which models are best. When comparing
#    the limited a5 and not limited a5 all of the different models contain CO2.
#    This tells me the bounds cannot be changed between DREAM and CMAES.

# Q. What must I do if I use the limited a5 bounds?
# A.
## - I must account when a5 parameter is near the upper bound i.e., it only
##   kicks in during the final years of calibration. 
## - My current assumption is if the a5 parameter switches on during the last 
##   10 years of calibration and has a relatively large impact on the AIC (> 4.60517 - evidence ratio > 10), 
##   then the model is better than the non-CO2 model.
## - This assumption will complicate the map plotting which currently 
##   uses AIC to calculate evidence ratio to distinguish between the models.
## - There are cases where the a5 parameter switching on during the last 10 years
##   outperforms the non-CO2 model with roughly the same AIC. 
## - The next best CO2 models for a given catchments also have the a5 parameter
##   switching on during the last 10 years.
## - Currently the map plot will suggest this model is better based off the 
##   difference in AIC. I don't think this is the case. I could swap
##   these CO2 model with the worst performing CO2 model. A possible problem
##   could occur if the worst performing model with an a5 kick in during the
##   last 10 years still outperforms the non-CO2 model. CHECK THIS. There is 
##   a single catchment where this occurs (i.e., all CO2 models with and a5
##   parameter within last 10 years always outperforms non-CO2). I have not
##   replaced this with a worst performing catchment.


# Q. What must I do if I use the unlimited a5 bounds?
# A. 
## - Using the 20241130 results (a5 can exceed max calibrated CO2), there are no
#    catchments where the turned-off a5 value model outpeforms the non-CO2 model.
#    This tells me the AIC can differentiate between turned-off parameters in models. 
#    I do not have to account for this using unlimited a5 bounds.
## - DREAM is made more complex. It produces poor distribution plots and it
#    is harder to build adequate trace plots.


pacman::p_load(tidyverse, gridExtra)


# Import functions -------------------------------------------------------------
source("./Functions/utility.R")


# Import the calibrated .csv's -------------------------------------------------
parameter_results <- read_csv(
  "./Results/CMAES/cmaes_parameter_results.csv", 
  show_col_types = FALSE
) 


streamflow_results <- read_csv(
  "./Results/CMAES/cmaes_streamflow_results.csv", 
  show_col_types = FALSE
  )
 

data <- readr::read_csv(
  "./Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  show_col_types = FALSE
) |> 
  mutate(year = as.integer(year))


start_stop_indexes <- readr::read_csv(
  "./Data/Tidy/start_end_index.csv",
  show_col_types = FALSE
)


# 1. Are there any major issues with CMAES fitting? - visual inspection --------
## streamflow results only include data calibrated on 
## page-per-gauge with unique models



## rainfall-runoff comparison ==================================================
rainfall_runoff_plot <- function(gauge, streamflow_results, parameter_results) {
  
  gauge_streamflow_results <- streamflow_results |> 
    filter(gauge == {{ gauge }}) |> 
    pivot_longer(
      cols = starts_with("transformed"),
      names_to = "observed_or_modelled",
      values_to = "streamflow"
    ) 
  
  # Get AIC information
  AIC_for_graphs <- parameter_results |> 
    filter(gauge == {{ gauge }}) |> 
    distinct(streamflow_model, .keep_all = TRUE) |> 
    select(streamflow_model, AIC) |> 
    add_column(
      x_pos_rainfall_runoff = gauge_streamflow_results |> pull(precipitation) |> min(), # minimum observed rainfall
      y_pos = gauge_streamflow_results |> pull(realspace_observed_streamflow) |> max(), # maximum observed streamflow
      x_pos_streamflow_time = 1959 # first value
    ) |> 
    mutate(
      label = paste0("AIC = ", round(AIC, digits = 2))
    )
  
  
  # Make plot
  gauge_streamflow_results |> 
    ggplot(aes(x = precipitation, y = streamflow, colour = observed_or_modelled)) +
    geom_smooth(
      method = lm,
      formula = y ~ x,
      se = FALSE
    ) +
    geom_point(size = 0.5) +
    geom_label(
      aes(x = x_pos_rainfall_runoff, y = y_pos, label = label),
      data = AIC_for_graphs,
      inherit.aes = FALSE,
      size = 3,
      vjust = 1,
      hjust = -0.25
    ) +
    scale_color_brewer(palette = "Set1") +
    theme_bw() +
    labs(
      x = "Precipitation (mm)",
      y = "Log-sinh streamflow (mm)",
      colour = NULL,
      title = paste0("Gauge: ", gauge) 
    ) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 5),
      plot.title = element_text(hjust = 0.5)
    ) +
    facet_wrap(~streamflow_model, scales = "free_y")
}


rainfall_runoff_plots <- map(
  .x = parameter_results |> pull(gauge) |> unique(),
  .f = rainfall_runoff_plot,
  streamflow_results = streamflow_results,
  parameter_results = parameter_results
)


ggsave(
  filename = "all_rainfall_runoff_plots.pdf", 
  path = "./Graphs/CMAES_graphs",
  # used to append all plots into a single pdf
  plot = gridExtra::marrangeGrob(rainfall_runoff_plots, nrow = 1, ncol = 1), 
  device = "pdf",
  units = "mm",
  width = 297,
  height = 210
)


# check streamflow of 318150
# utilisation of a parameter? - remove?


## streamflow-time comparison ==================================================
streamflow_time_plot <- function(gauge, streamflow_results, parameter_results) {
  
  gauge_streamflow_results <- streamflow_results |> 
    filter(gauge == {{ gauge }}) |> 
    pivot_longer(
      cols = starts_with("realspace"),
      names_to = "observed_or_modelled",
      values_to = "streamflow"
    ) 
  
  # Get AIC information
  AIC_for_graphs <- parameter_results |> 
    filter(gauge == {{ gauge }}) |> 
    distinct(streamflow_model, .keep_all = TRUE) |> 
    select(streamflow_model, AIC) |> 
    add_column(
      x_pos_rainfall_runoff = gauge_streamflow_results |> pull(precipitation) |> min(), # minimum observed rainfall
      y_pos = gauge_streamflow_results |> filter(observed_or_modelled == "realspace_observed_streamflow") |> pull(streamflow) |> max(), # maximum observed streamflow
      x_pos_streamflow_time = 1959 # first value
    ) |> 
    mutate(
      label = paste0("AIC = ", round(AIC, digits = 2))
    )
  
  
  # Make plot
  gauge_streamflow_results |> 
    ggplot(aes(x = year, y = streamflow, colour = observed_or_modelled)) +
    geom_line() +
    geom_point(size = 0.5) +
    geom_label(
      aes(x = x_pos_streamflow_time, y = y_pos, label = label),
      data = AIC_for_graphs,
      inherit.aes = FALSE,
      size = 3,
      vjust = 1,
      hjust = -0.25
    ) +
    scale_color_brewer(palette = "Set1") +
    theme_bw() +
    labs(
      x = "Year",
      y = "Streamflow (mm)",
      colour = NULL,
      title = paste0("Gauge: ", gauge) 
    ) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 5),
      plot.title = element_text(hjust = 0.5)
    ) +
    facet_wrap(~streamflow_model, scales = "free_y")
}



streamflow_time_plots <- map(
  .x = parameter_results |> pull(gauge) |> unique(),
  .f = streamflow_time_plot,
  streamflow_results = streamflow_results,
  parameter_results = parameter_results
)


ggsave(
  filename = "all_streamflow_time_plots.pdf", 
  path = "./Graphs/CMAES_graphs",
  # used to append all plots into a single pdf
  plot = gridExtra::marrangeGrob(streamflow_time_plots, nrow = 1, ncol = 1), 
  device = "pdf",
  units = "mm",
  width = 297,
  height = 210
)



# 2. Are the fitted parameters "acceptable"? -----------------------------------

## Utilisation of parameters ===================================================
# What parameters are being turned off that should not be turned off? ##########
check_near_zero_parameter_values <- parameter_results |> 
  mutate(
    is_parameter_turned_off = near(parameter_value, 0, tol = .Machine$double.eps^0.5)
  ) |> 
  filter(is_parameter_turned_off)

# The a parameter wants to be turned off for 512 catchment-model combinations (out of 9684)
#parameter_results |> 
#  filter(parameter == "a") |> 
#  pull() |> 
#  length()


# Check the a5 parameter - how many catchment have a5 turn on in the last ######
# year of observed data

# I think this is only relevent for the best models...







  
# 3. Save best CO2 and non-CO2 results per gauge -------------------------------
best_CO2_and_non_CO2_per_catchment <- parameter_results |>
  select(!c(optimiser, loglikelihood, exit_message, near_bounds, objective_function)) |>
  mutate(
    contains_CO2 = str_detect(streamflow_model, "CO2"),
    .after = 2
  ) |>
  slice_min(
    AIC,
    by = c(gauge, contains_CO2)
  ) 

write_csv(
  best_CO2_and_non_CO2_per_catchment,
  file = "./Results/CMAES/best_CO2_non_CO2_per_catchment_CMAES.csv"
)








###### DELETE FROM HERE ##########


####### CAN REMOVE FROM HERE #########









## Check if parameters are being set to zero ===================================
### This represents one ways of parameters being "turned-off"
utilisation_tolerence <- 1E-4 

# if below utilisation_tolerence then we assume the parameter has been "turned-off"

near_zero_parameters <- best_CO2_and_non_CO2_per_catchment |> 
  filter(abs(parameter_value) <= utilisation_tolerence) |> 
  filter(!parameter %in% c("a5", "a3_slope")) 
## The a5 parameter can be near-zero. 
## The a3_slope is also very small. Ignore it from near zero test.
## It indicates the CO2 term kicks in at the start of the timeseries.
## Other parameters vary. I think this is okay. The a3 is a bit iffy



## THIS MAY BE OMITTED ###






## Check if a5 parameter is only being utilised in the last year ===============
### This represents another way of parameters being "turned-off"

### Assumption:
### - If the a5 parameter is being switched on the the last year
###   it probably is not the best model

### Testing:
### - If the a5 parameter switches on during the last 10 years of calibration
###   and has a relatively large impact on the AIC (> 5), then the model is better
###   than the non-CO2 model.
### * > 4.60... selected because when converted to evidence ratio [exp(-0.5 * -4.60...) = ~ 10]
###   10 is moderate evidence begins


### Extract indexes to get the CO2 for the last and 10 last year during calibration ###
last_row_start_stop <- start_stop_indexes |> 
  slice_tail(by = gauge) |> 
  select(c(gauge, end_index)) |> 
  mutate(
    second_last_index = end_index - 5 # ori -> 10 year before end
  )


# This is bad...
last_second_last_char <- function(years) {
  # only works with two years
  if(years[2] - years[1] >= 1) {
    return(c("second_last", "last"))
  } else {
    return(c("last", "second_last"))
  }
}


### Range of CO2 - 10 last and last per gauge ###
last_range_CO2_per_gauge <- data |> 
  select(c(gauge, year, CO2)) |> 
  left_join(last_row_start_stop, by = join_by(gauge)) |> 
  slice(c(second_last_index, end_index), .by = gauge) |> 
  distinct() |> # this is not an elegant solution
  mutate(
    char_index = last_second_last_char(year),
    .by = gauge
  ) |> 
  select(gauge, CO2, char_index) 



### Examine if the fitted a5 is within the CO2 values during the last 10 years of calibration ###
upper_a5_parameter <- best_CO2_and_non_CO2_per_catchment |> 
  filter(parameter == "a5") |> 
  left_join(last_range_CO2_per_gauge, by = join_by(gauge)) |> 
  rename(
    upper_a5_value = CO2
  ) |> 
  arrange(desc(parameter_value)) |> 
  pivot_wider(
    names_from = char_index,
    values_from = upper_a5_value
  ) |> 
  relocate(
    c(second_last, last),
    .after = parameter_value
  ) |> 
  mutate(
    between_last_and_second_last = between(parameter_value, second_last, last),
    .by = gauge,
    .before = AIC
  ) # a5 = 131.65 is an value in CO2. Seeing this repeated sketches me out


### Filter catchments where a5 parameter kick in during the last 10 years of calibration ###
catchments_where_a5_kicks_in_last_year <- upper_a5_parameter |>
  filter(between_last_and_second_last) |> 
  pull(gauge)


### Compare if the catchments_where_a5_kicks_in_last_year outperform non-CO2 catchments ### 
evidence_ratio_threshold <- 4.60517


compare_late_a5_kick_in_to_non_CO2 <- best_CO2_and_non_CO2_per_catchment |> 
  filter(gauge %in% catchments_where_a5_kicks_in_last_year) |> 
  select(gauge, contains_CO2, AIC) |> 
  distinct() |> 
  pivot_wider(
    names_from = contains_CO2,
    values_from = AIC
  ) |> 
  rename(
    contains_CO2 = `TRUE`,
    no_CO2 = `FALSE`
  ) |> 
  mutate(
    CO2_better = contains_CO2 < no_CO2,
    delta_AIC = no_CO2 - contains_CO2,
    non_CO2_better_model = delta_AIC <= evidence_ratio_threshold # corresponds to and evidence ratio > 10 (exp(-0.5 * -4.60517))
  ) |> 
  filter(CO2_better) |> 
  arrange(desc(delta_AIC))


## If the CO2 model is not significantly better (AIC > 4.6...), the
## we assume the CO2 model does not outperform the non-CO2 model because it is
## not fully utilising the parameters.

## Replace non-utilised CO2 models with another model ==========================

### Examine if any other models for a given non-utilised gauge has a
### a5 value kicking in before the final 10-years of calibration

non_utilised_a5_gauges <- compare_late_a5_kick_in_to_non_CO2 |> 
  filter(non_CO2_better_model) |> 
  pull(gauge)


best_CO2_per_catchment <- best_CO2_and_non_CO2_per_catchment |> 
  filter(contains_CO2) |> 
  select(gauge, streamflow_model) |> 
  add_column(best_model = TRUE) |> 
  distinct()

# If the next best CO2 models have a a5 value greater than the best instant exclude
exclude_a5_greater_than_best <- function(a5_values) { 
  # relies on arrange(gauge, AIC to work)
  best_value <- a5_values[1]
  return(best_value - a5_values <= 0)
}

# Then check if the new AIC would make a difference
best_non_CO2_per_catchment <- best_CO2_and_non_CO2_per_catchment |> 
  filter(!contains_CO2) |> 
  select(gauge, AIC) |> 
  rename(
    non_CO2_AIC = AIC
  ) |> 
  distinct()

alternative_CO2_models_non_utilised_a5 <- parameter_results |> 
  filter(parameter == "a5") |> 
  filter(gauge %in% non_utilised_a5_gauges) |> 
  select(gauge, streamflow_model, parameter_value, AIC) |> 
  left_join(
    best_CO2_per_catchment,
    by = join_by(gauge, streamflow_model)
    ) |> 
  mutate(
    best_model = if_else(is.na(best_model), FALSE, TRUE)
  ) |> 
  arrange(gauge, AIC) |>  # add column with current_best
  mutate( # If the next best CO2 models have a a5 value greater than or equal to the best instant exclude
    instant_exclude = exclude_a5_greater_than_best(parameter_value),
    .by = gauge
  ) |> 
  filter(!instant_exclude) 

# Q. Are there any alternative_CO2_models_non_utilised_a5 that would exceed
# the evidence_ratio_threshold?

exceed_evi_ratio_alternative_CO2_models_non_utilised_a5 <- alternative_CO2_models_non_utilised_a5 |> 
  left_join( # Then check if the new AIC would make a difference
    best_non_CO2_per_catchment,
    by = join_by(gauge) 
  ) |> 
  mutate(
    delta_AIC = non_CO2_AIC - AIC,
    change_model = delta_AIC > evidence_ratio_threshold
  ) |> 
  filter(change_model)

# A. No


## Replace models with a non-utilised a5 parameter =============================
### The gauges with a non-utilised a5 parameter will be replaced using two
### methods.

### The first method involves using the lowest AIC models from alternative_CO2_models_non_utilised_a5
replace_alternative_CO2_models_non_utilised_a5 <- alternative_CO2_models_non_utilised_a5 |> 
  slice_min(AIC, by = gauge)

### The second methods involves using the worst performing AIC CO2 models for a 
### given gauges. 
### This was done because all gauges not in the replace_alternative_CO2_models_non_utilised_a5 
### have an a5 parameter within the last 10-years of calibration. 
### We have deemed this as non-utilised. 

### Assumption: 
### Replacing the non-utilised models with the worst performing non-utilised a5 
### models should mean the non-CO2 model is now best. 

### Remove replace_alternative_CO2_models_non_utilised_a5 gauges from 
### replace_non_CO2_outperform_worst_CO2

worse_performing_CO2_models_off_late_kick_ins <- parameter_results |> 
  filter(parameter == "a5") |> 
  filter(gauge %in% non_utilised_a5_gauges) |> 
  select(gauge, streamflow_model, parameter_value, AIC) |> 
  slice_max(AIC, by = gauge) |> 
  filter(!gauge %in% replace_alternative_CO2_models_non_utilised_a5$gauge)

### Check assumption ### 
# Do the non-CO2 models outperform the worst performing CO2 models?

non_CO2_outperform_worst_CO2 <- worse_performing_CO2_models_off_late_kick_ins |> 
  left_join(
    best_non_CO2_per_catchment, 
    by = join_by(gauge)
    ) |> 
  mutate(
    non_CO2_better = non_CO2_AIC < AIC
  ) 

non_CO2_outperform_worst_CO2 |> filter(!non_CO2_better)
# Assumption is correct

replace_non_CO2_outperform_worst_CO2 <- non_CO2_outperform_worst_CO2 |> filter(non_CO2_better)


## Modify best_CO2_and_non_CO2_per_catchment ===================================
### Models from:
### - replace_alternative_CO2_models_non_utilised_a5
### - replace_non_CO2_outperform_worst_CO2

### Find and remove gauge and CO2 model combinations from best_CO2_and_non_CO2_per_catchment ###
replace_gauges <- c(
  replace_alternative_CO2_models_non_utilised_a5$gauge, 
  replace_non_CO2_outperform_worst_CO2$gauge
  )

non_utilised_CO2_models_filtered <- best_CO2_and_non_CO2_per_catchment |> 
  filter(gauge %in% replace_gauges) |> 
  filter(contains_CO2)

best_CO2_and_non_CO2_per_catchment_without_non_utilised <- best_CO2_and_non_CO2_per_catchment |> 
  anti_join(
    non_utilised_CO2_models_filtered,
    by = join_by(gauge, streamflow_model)
  )

### Extract replace_alternative_CO2_models_non_utilised_a5 from parameter results ###
extracted_replace_alternative_CO2_models_non_utilised_a5 <- parameter_results |> 
  semi_join(
    replace_alternative_CO2_models_non_utilised_a5,
    by = join_by(gauge, streamflow_model)
  ) |>
  select(!c(optimiser, loglikelihood, exit_message, near_bounds)) |>
  mutate(
    contains_CO2 = str_detect(streamflow_model, "CO2"),
    .after = 2
  ) |>
  slice_min(
    AIC,
    by = c(gauge, contains_CO2)
  ) |> 
  select(!objective_function) |> 
  distinct()

  

### Extract replace_non_CO2_outperform_worst_CO2 ###
extracted_replace_non_CO2_outperform_worst_CO2 <- parameter_results |> 
  semi_join(
    replace_non_CO2_outperform_worst_CO2,
    by = join_by(gauge, streamflow_model)
  ) |>
  select(!c(optimiser, loglikelihood, exit_message, near_bounds)) |>
  mutate(
    contains_CO2 = str_detect(streamflow_model, "CO2"),
    .after = 2
  ) |>
  slice_min(
    AIC,
    by = c(gauge, contains_CO2)
  ) |> 
  select(!objective_function) |> 
  distinct()


### Add flag to altered gauges/models ##########################################
with_flag <- rbind(
  extracted_replace_alternative_CO2_models_non_utilised_a5, 
  extracted_replace_non_CO2_outperform_worst_CO2
) |> 
  add_column(flag = "modified") |> 
  select(gauge, streamflow_model, contains_CO2, flag) |> 
  distinct()


## The modified_best_CO2_and_non_CO2_per_catchment =============================
### Join everything together
modified_best_CO2_and_non_CO2_per_catchment <- rbind(
  best_CO2_and_non_CO2_per_catchment_without_non_utilised,
  extracted_replace_alternative_CO2_models_non_utilised_a5,
  extracted_replace_non_CO2_outperform_worst_CO2
) |> 
  left_join(
    with_flag,
    by = join_by(gauge, streamflow_model, contains_CO2)
  ) |> 
  arrange(gauge)


### CHECKED - I am not missing any ### 

## Do my changes cause any big swings in AIC differences? ======================
## If so it probably requires a rethink
original_evi_ratio <- best_CO2_and_non_CO2_per_catchment |> 
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
  select(gauge, evidence_ratio) |> 
  rename(original_evidence_ratio = evidence_ratio)


flag_per_gauge <- modified_best_CO2_and_non_CO2_per_catchment |> 
  select(gauge, flag) |> 
  distinct() |> 
  filter(flag == "modified")


swings_in_evi_ratio <- modified_best_CO2_and_non_CO2_per_catchment |>
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
  # Add flag back in
  left_join(
    flag_per_gauge,
    by = join_by(gauge)
  ) |> 
  mutate(
    flag = if_else(is.na(flag), "unmodified", flag)
  ) |> 
  rename(
    modified_evi_ratio = evidence_ratio
  ) |> 
  select(gauge, modified_evi_ratio, flag) |> 
  arrange(flag, modified_evi_ratio) |> 
  left_join(
    original_evi_ratio,
    by = join_by(gauge)
  ) |> 
  relocate(
    flag,
    .after = 4
  )



# 3. Check if parameters are near bounds ---------------------------------------
## for modified_best_CO2_and_non_CO2_per_catchment
near_bounds_parameters <- parameter_results |> 
  semi_join(
    modified_best_CO2_and_non_CO2_per_catchment,
    by = join_by(gauge, streamflow_model)
  ) |> 
  filter(near_bounds) 

# There is catchment near bounds.

# Use to be a single catchment (112101B) with bounds at boundary
# Catchment has really high rainfall. Rainfall > bc_q
# Redo and replace the .csv directly saves time

# 4. Save modified_best_CO2_and_non_CO2_per_catchment ------------------------
# For the time being remove near bounds --> 112101B

write_csv(
  modified_best_CO2_and_non_CO2_per_catchment,
  file = paste0("./Results/my_cmaes/best_CO2_non_CO2_per_catchment_CMAES_", get_date(), ".csv")
)



### REMOVE WHEN READY ###
stop_here <- "stop_here"

# 4. Other stuff during testing like: 
## Comparing with limited a5 and non-limited a5 bounds. 

# Using the a5 bounds that allow a5 to be turned-off, are there any 
# cases where the a5 parameter being turned-off outperforms the non-CO2 model? A. no
# Steps:
# 1. Find catchments where the CO2 model outperformed the non-CO2 model
gauges_CO2_model_better <- best_CO2_and_non_CO2_per_catchment |> 
  select(gauge, contains_CO2, AIC) |>
  distinct() |> 
  pivot_wider(
    names_from = contains_CO2,
    values_from = AIC
  ) |> 
  rename(
    CO2_AIC = `TRUE`,
    non_CO2_AIC = `FALSE`
  ) |> 
  mutate(
    CO2_is_better = CO2_AIC < non_CO2_AIC
  ) |> 
  filter(CO2_is_better) |> 
  pull(gauge)

# 2. Find the maximum CO2 value used in calibration
last_row_start_stop <- start_stop_indexes |> 
  slice_tail(by = gauge) |> 
  select(c(gauge, end_index)) 

# extract the last CO2 values used in calibration per gauge
last_CO2_per_gauge <- data |> 
  select(c(gauge, year, CO2)) |> 
  left_join(last_row_start_stop, by = join_by(gauge)) |> 
  slice(end_index, .by = gauge) |> 
  distinct() |> # this is not an elegant solution
  select(gauge, CO2) 

# 3. Compare if the best a5 values exceeds the maximum CO2 value in calibration
x <- best_CO2_and_non_CO2_per_catchment |> 
  filter(gauge %in% gauges_CO2_model_better) |> 
  left_join(
    last_CO2_per_gauge,
    by = join_by(gauge)
  ) |> 
  rename(
    max_calibration_CO2 = CO2
  ) 

y <- x |> 
  filter(parameter == "a5") |> 
  mutate(
    a5_value_exceeds_max_calibration = parameter_value > max_calibration_CO2
  )

# 4. If yes, then the non-limited a5 is not better than the limited a5 bounds.
## The answer is no...
#unlimited_parameter_results <- read_csv(
#  "./Results/my_cmaes/CMAES_parameter_results_20241130.csv", 
#  # limited bounds 20250122
#  # turn-off bounds 20241130
#  show_col_types = FALSE
#) |> 
#  filter(objective_function == "constant_sd_objective_function")


#x <- limited_best_CO2_and_non_CO2_per_catchment |> 
#  select(gauge, contains_CO2, streamflow_model, AIC) |> 
#  distinct()

#unlimited_best_CO2_and_non_CO2_per_catchment <- unlimited_parameter_results |>
#  select(!c(optimiser, loglikelihood, exit_message, near_bounds)) |>
#  unite(
#    col = streamflow_model_objective_function,
#    c(streamflow_model, objective_function),
#    sep = "_",
#    remove = FALSE,
#    na.rm = FALSE
#  ) |>
#  mutate(
#    contains_CO2 = str_detect(streamflow_model_objective_function, "CO2"),
#    .after = 2
#  ) |>
#  slice_min(
#    AIC,
#    by = c(gauge, contains_CO2)
#  ) |> 
#  select(!objective_function) |> 
#  distinct()

#y <- unlimited_best_CO2_and_non_CO2_per_catchment |> 
#  select(gauge, contains_CO2, streamflow_model, AIC) |> 
#  distinct()

# Compare limited and unlimited best CO2 using anti-join
#z <- x |> 
#  anti_join(
#    y,
#    by = join_by(gauge, streamflow_model)
#  )

# Yeah, they are different. All of the different models contain CO2.
# This tells me the bounds cannot be changed between DREAM and CMAES.