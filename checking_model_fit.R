# Questions script answers:
# 1. How did the fitting perform?
# 2. Are the fitted parameters acceptable?
# 3. Can AIC accurately determine the best CO2 and non-CO2 model? 

# TODO
# Document changes in graphs 
# Make formatting nicer
# changes name of script to 02-check-model-fit.R
# re-run catchments near bounds
# save best_streamflow_model.csv (DREAM and streamflow plots) and best_CO2_and_CO2_model 
# using manually adjusted results based on parameters etc. Include flags done in this script
# install ozmap and patchwork packages
# evidence ratio is dependent on the co2 and nonCO2 AIC's. There are a few 
# catchments where the a5 value only kicks in during the last year of calibration
# with neglibile change in AIC. For some reason this has a better AIC value (smaller). 
# This does not make sense to me. Leave the map graph as is or remove these catchments.
# Maybe leave the a5 values gauge swaps for now. Tell Tim/Murray
# I have forced it to converge within the bounds.
# None of the graphs with a5 have weird upticks which is nice.

pacman::p_load(tidyverse, gridExtra)


# Import functions -------------------------------------------------------------
source("./Functions/utility.R")


# Import the calibrated .csv's -------------------------------------------------
parameter_results <- read_csv(
  "./Results/my_cmaes/CMAES_parameter_results_20250122.csv", 
  show_col_types = FALSE
)


streamflow_results <- read_csv(
  "./Results/my_cmaes/CMAES_streamflow_results_20250122.csv", 
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


# 1. Answer to "How did the fitting perform?" -------------------------------------

# Only include streamflow that was calibrated on -------------------------------
## join the included_in_calibration column
in_calibration <- data |> 
  select(year, gauge, included_in_calibration)

only_calibration_streamflow_results <-  streamflow_results |> 
  left_join(
    in_calibration,
    by = join_by(year, gauge)
  ) |> 
  filter(included_in_calibration)
  


# Tidy data for plotting -------------------------------------------------------
setup_streamflow_results_plotting <- only_calibration_streamflow_results |> 
  pivot_longer(
    cols = ends_with("streamflow"),
    names_to = "modelled_or_observed",
    values_to = "boxcox_streamflow"
  ) 



## Get unique model and objective function combinations ========================
unique_streamflow_model_objective_combinations <- setup_streamflow_results_plotting |>
  distinct(streamflow_model, objective_function) 


# Plotting function ============================================================
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

# top and tail function



# Results ----------------------------------------------------------------------
many_plots <- map2(
  .x = unique_streamflow_model_objective_combinations |> dplyr::pull(streamflow_model),
  .y = unique_streamflow_model_objective_combinations |> dplyr::pull(objective_function),
  .f = check_results_plot,
  streamflow_results = setup_streamflow_results_plotting
)




# Save and append all plots. Solutions required gridExtra function
ggsave(
  filename = paste0("check_model_fit_", get_date(), ".pdf"), 
  path = "./Graphs/CMAES_graphs",
  plot = gridExtra::marrangeGrob(many_plots, nrow = 1, ncol = 1), 
  device = "pdf",
  units = "mm",
  width = 1189,
  height = 841
)







# 2. Answer to "Are the fitted parameters acceptable?" -------------------------





# Check if parameters are near bounds ------------------------------------------
near_bounds_parameters <- parameter_results |> 
  filter(near_bounds) |> 
  filter(parameter != "a5")

# This is not good.
# Problem catchment is 112101B. Really high rainfall. Rainfall > bc_q
# Redo and replace the .csv directly


# Parameter utilisation --------------------------------------------------------
best_CO2_and_non_CO2_per_catchment <- parameter_results |>
  select(!c(optimiser, loglikelihood, exit_message, near_bounds)) |>
  unite(
    col = streamflow_model_objective_function,
    c(streamflow_model, objective_function),
    sep = "_",
    remove = FALSE,
    na.rm = FALSE
  ) |>
  mutate(
    contains_CO2 = str_detect(streamflow_model_objective_function, "CO2"),
    .after = 2
  ) |>
  slice_min(
    AIC,
    by = c(gauge, contains_CO2)
  ) |> 
  select(!objective_function) |> 
  distinct()

## Check if parameters are being set to zero ===================================
utilisation_tolerence <- 1E-3
near_zero_parameters <- best_CO2_and_non_CO2_per_catchment |> 
  filter(abs(parameter_value) <= utilisation_tolerence)

## The a5 parameter can be near-zero. It indicates the CO2 term kicks in at the start of the timeseries.
## Other parameters vary. I think this is okay. The a3 is a bit iffy
unique(near_zero_parameters$parameter)


## Check if a5 parameter is only being utilised in the last year ===============

## Assumption:
### - If the a5 parameter is being switched on the the last year
###   it probably is not the best model
### - If the a5 parameter switches on during the last 10 years of calibration
###   and has a relatively large impact on the AIC (> 5), then the model is better
###   than the non-CO2 model.
### * > 5 selected because when converted to evidence ratio [exp(-0.5 * -5) = ~ 10]
###   10 is moderate evidence begins


# extract last row per gauge from start_stop_index
last_row_start_stop <- start_stop_indexes |> 
  slice_tail(by = gauge) |> 
  select(c(gauge, end_index)) |> 
  mutate(
    second_last_index = end_index - 10 # 10 year before end
  )


last_second_last_char <- function(years) {
  # only works with two years
  if(years[2] - years[1] >= 1) {
    return(c("second_last", "last"))
  } else {
    return(c("last", "second_last"))
  }
}


# extract the last CO2 values used in calibration per gauge
last_CO2_per_gauge <- data |> 
  select(c(gauge, year, CO2)) |> 
  left_join(last_row_start_stop, by = join_by(gauge)) |> 
  slice(c(second_last_index, end_index), .by = gauge) |> 
  distinct() |> # this is not an elegant solution
  mutate(
    char_index = last_second_last_char(year),
    .by = gauge
  ) |> 
  select(gauge, CO2, char_index) 




# join and compare to fitted a5 parameter
### a5 = 131.65 is an value in CO2. Seeing this repeated sketches me out
upper_a5_parameter <- best_CO2_and_non_CO2_per_catchment |> 
  filter(parameter == "a5") |> 
  left_join(last_CO2_per_gauge, by = join_by(gauge)) |> 
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
  )

# Filter catchments where a5 kick in during the last year of calibration
catchments_where_a5_kicks_in_last_year <- upper_a5_parameter |>
  filter(between_last_and_second_last) |> 
  pull(gauge)


# Compare if the last kick in's are better than non-CO2s
compare_late_a5_kick_in_to_non_CO2 <- best_CO2_and_non_CO2_per_catchment |> 
  filter(gauge %in% catchments_where_a5_kicks_in_last_year) |> 
  select(
    gauge, 
    contains_CO2, 
    AIC
    ) |> 
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
    non_CO2_better_model = delta_AIC <= 5
  ) |> 
  filter(CO2_better) |> 
  arrange(desc(delta_AIC))




# 3. Answer to "Can AIC accurately determine the best CO2 and non-CO2 model?" ----
# No I will have to do some manual changes



# Try the next best-CO2 models
poor_performing_late_a5_gauges <- compare_late_a5_kick_in_to_non_CO2 |> 
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

compare_next_best_CO2_models <- parameter_results |> 
  filter(parameter == "a5") |> 
  filter(gauge %in% poor_performing_late_a5_gauges) |> 
  select(gauge, streamflow_model, parameter_value, AIC) |> 
  left_join(
    best_CO2_per_catchment,
    by = join_by(gauge, streamflow_model)
    ) |> 
  mutate(
    best_model = if_else(is.na(best_model), FALSE, TRUE)
  ) |> 
  arrange(gauge, AIC) |>  # add column with current_best
  mutate( # If the next best CO2 models have a a5 value greater than the best instant exclude
    instant_exclude = exclude_a5_greater_than_best(parameter_value),
    .by = gauge
  ) |> 
  filter(!instant_exclude) |> 
  left_join( # Then check if the new AIC would make a difference
    best_non_CO2_per_catchment,
    by = join_by(gauge) 
  ) |> 
  mutate(
    delta_AIC = non_CO2_AIC - AIC,
    change_model = delta_AIC > 5
  ) |> 
  filter(change_model)





# From `compare_late_a5_kick_in_to_non_CO2` there are 5 catchments
# where the a5 value kicks in in the last year of calibration that are
# better than the non-CO2 model. I find this hard to believe.


# Options:
# 1. Leave it. DREAM will show the distributions.
# 2. Swap the best-CO2 and best-non-CO2 for that given gauge for all
# 3. Swap the best-CO2 and best-non-CO2 for a given gauge if delta_AIC < 5 (exp(-0.5 * -5) = ~ 10 starting to be meaning after this)
# 4. Try the next best CO2 model, see if it kick-in on the last year (very likely a5 kicks in around the same time) - my guess was correct

# Complete option 4. My guess was correct
# I will do option 3. Changes required when graphing best (maps and timeseries) and DREAM


