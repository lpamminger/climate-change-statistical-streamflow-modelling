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
  "./Results/my_cmaes/CMAES_parameter_results_20250122.csv", 
  # limited bounds 20250122
  # turn-off bounds 20241130
  show_col_types = FALSE
) 


streamflow_results <- read_csv(
  "./Results/my_cmaes/CMAES_streamflow_results_20241130.csv", 
  # limited bounds 20250122
  # turn-off bounds 20241130
  show_col_types = FALSE
  )
 

data <- read_csv(
  "Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  show_col_types = FALSE
)


start_stop_indexes <- read_csv(
  "./Data/Tidy/start_end_index.csv",
  show_col_types = FALSE
)


# 1. Are there any major issues with CMAES fitting? - visual inspection --------

## Only include streamflow that was calibrated on ==============================

in_calibration <- data |> 
  select(year, gauge, included_in_calibration) # join the included_in_calibration column

only_calibration_streamflow_results <-  streamflow_results |> 
  left_join(
    in_calibration,
    by = join_by(year, gauge)
  ) |> 
  filter(included_in_calibration)
  


## Tidy data for plotting ======================================================
setup_streamflow_results_plotting <- only_calibration_streamflow_results |> 
  pivot_longer(
    cols = ends_with("streamflow"),
    names_to = "modelled_or_observed",
    values_to = "boxcox_streamflow"
  ) 



## Get unique model and objective function combinations ========================
unique_streamflow_model_objective_combinations <- setup_streamflow_results_plotting |>
  distinct(streamflow_model, objective_function) 


## Plotting function ===========================================================
check_results_plot <- function(streamflow_model, objective_function, streamflow_results) {
  streamflow_results |>
    filter(streamflow_model == {{ streamflow_model }}) |>
    filter(objective_function == {{ objective_function }}) |>
    ggplot(aes(x = year, y = boxcox_streamflow, colour = modelled_or_observed)) +
    geom_line(na.rm = TRUE) +
    geom_point(na.rm = TRUE) +
    scale_colour_brewer(palette = "Set1") +
    labs(
      x = "Year",
      y = "Boxcox Streamflow (mm)",
      title = paste0("Streamflow Model: ", streamflow_model, "\nObjective Function: ", objective_function)
    ) +
    theme_bw() +
    facet_wrap(~gauge, scales = "free_y") +
    theme(
      legend.position = "top"
    )
}





## Run and Save plots ==========================================================
many_plots <- map2(
  .x = unique_streamflow_model_objective_combinations |> dplyr::pull(streamflow_model),
  .y = unique_streamflow_model_objective_combinations |> dplyr::pull(objective_function),
  .f = check_results_plot,
  streamflow_results = setup_streamflow_results_plotting
)


ggsave(
  filename = paste0("check_model_fit_", get_date(), ".pdf"), 
  path = "./Graphs/CMAES_graphs",
  plot = gridExtra::marrangeGrob(many_plots, nrow = 1, ncol = 1), 
  # used to append all plots into a single pdf
  device = "pdf",
  units = "mm",
  width = 1189,
  height = 841
)







# 2. Are the fitted parameters "acceptable"? -----------------------------------

## Only concerned with the best CO2 and non CO2 models - extract ===============
best_CO2_and_non_CO2_per_catchment <- parameter_results |>
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



## Check if parameters are being set to zero ===================================
### This represents one ways of parameters being "turned-off"
utilisation_tolerence <- 1E-4 

# if below utilisation_tolerence then we assume the parameter has been "turned-off"

near_zero_parameters <- best_CO2_and_non_CO2_per_catchment |> 
  filter(abs(parameter_value) <= utilisation_tolerence) |> 
  filter(parameter != "a5") 
## The a5 parameter can be near-zero. 
## It indicates the CO2 term kicks in at the start of the timeseries.
## Other parameters vary. I think this is okay. The a3 is a bit iffy



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
    second_last_index = end_index - 10 # 10 year before end
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



## The modified_best_CO2_and_non_CO2_per_catchment =============================
### Join everything together
modified_best_CO2_and_non_CO2_per_catchment <- rbind(
  best_CO2_and_non_CO2_per_catchment_without_non_utilised,
  extracted_replace_alternative_CO2_models_non_utilised_a5,
  extracted_replace_non_CO2_outperform_worst_CO2
)

### CHECKED ## 




# 3. Check if parameters are near bounds ---------------------------------------
## for modified_best_CO2_and_non_CO2_per_catchment
near_bounds_parameters <- parameter_results |> 
  semi_join(
    modified_best_CO2_and_non_CO2_per_catchment,
    by = join_by(gauge, streamflow_model)
  ) |> 
  filter(near_bounds) 

# There is a single catchment (112101B) with bounds at boundary
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