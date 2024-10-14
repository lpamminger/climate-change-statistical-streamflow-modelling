# Testing/Implementing new model
cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, cmaesr, smoof, tictoc, furrr, parallel, truncnorm, sloop, tictoc)

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
  show_col_types = FALSE,
  col_select = c(gauge, bc_lambda, state, lat, lon)
)





## Utility functions ===========================================================
source("./Functions/utility.R")
source("./Functions/boxcox_transforms.R")


## Import streamflow functions =================================================
source("./Functions/streamflow_models.R")
source("./Functions/parameter_transformations.R")
source("./Functions/catchment_data_blueprint.R")
source("./Functions/cmaes_dream_summaries.R")
source("./Functions/objective_functions.R")
source("./Functions/numerical_optimiser_setup.R")
source("./Functions/generic_functions.R")
source("./Functions/my_cmaes.R")
source("./Functions/objective_function_setup.R")
source("./Functions/result_set.R")





# Steps to implementing new model ----------------------------------------------
# a5 is the CO2 shift variable
# I have hard coded the a5 parameter in make_parameter_bounds - this should require the user to input Co2 data to work out bounds.



# Test function ----------------------------------------------------------------

# 1.
streamflow_model_separate_shifted_CO2 <- function(catchment_data, parameter_set) { 
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "test_streamflow_model_separate_CO2",
        "parameters" = c("a0", "a1", "a3", "a5")
      )
    )
  } 
  
  # Check if parameter is a matrix - if not coerce
  if(!is.matrix(parameter_set)) {
    parameter_set <- as.matrix(parameter_set, ncol = 1)
  }
  
  # Parameters
  a0 <- parameter_set[1, ]
  a1 <- parameter_set[2, ]
  a3 <- parameter_set[3, ]
  a5 <- parameter_set[4, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  CO2 <- catchment_data$CO2
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) 
  repeat_CO2 <- matrix(CO2, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) 
  repeat_a0 <- matrix(a0, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a3 <- matrix(a3, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a5 <- matrix(a5, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Model
  repeat_shifted_CO2 <- repeat_CO2 - repeat_a5
  repeat_shifted_CO2[repeat_shifted_CO2 < 0] <- 0 # max(0, CO2 - a5)
  repeat_a0 + (repeat_a1 * repeat_precipitation) + (repeat_a3 * repeat_shifted_CO2)
  
}


# 2.
streamflow_model_separate_shifted_CO2_seasonal_ratio <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_separate_shifted_CO2_seasonal_ratio",
        "parameters" = c("a0", "a1", "a3", "a4", "a5")
      )
    )
  } 
  
  # Check if parameter is a matrix - if not coerce
  if(!is.matrix(parameter_set)) {
    parameter_set <- as.matrix(parameter_set, ncol = 1)
  }
  
  # Parameters
  a0 <- parameter_set[1, ]
  a1 <- parameter_set[2, ]
  a3 <- parameter_set[3, ]
  a4 <- parameter_set[4, ]
  a5 <- parameter_set[5, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  seasonal_ratio <- catchment_data$seasonal_ratio
  CO2 <- catchment_data$CO2
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) #matrix(rep(precipitation, times = ncol(parameter_set)), ncol = ncol(parameter_set), byrow = FALSE)
  repeat_seasonal_ratio <- matrix(seasonal_ratio, ncol = ncol(parameter_set), nrow = length(seasonal_ratio), byrow = FALSE)
  repeat_CO2 <- matrix(CO2, ncol = ncol(parameter_set), nrow = length(CO2), byrow = FALSE) 
  
  repeat_a0 <- matrix(a0, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a3 <- matrix(a3, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a4 <- matrix(a4, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a5 <- matrix(a5, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Model
  repeat_shifted_CO2 <- repeat_CO2 - repeat_a5
  repeat_shifted_CO2[repeat_shifted_CO2 < 0] <- 0 # max(0, CO2 - a5)
  repeat_a0 + (repeat_a1 * repeat_precipitation) + (repeat_a3 * repeat_shifted_CO2) + (repeat_a4 * repeat_seasonal_ratio)
  
}



# 3.
streamflow_model_drought_separate_shifted_CO2_seasonal_ratio_auto <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_drought_separate_shifted_CO2_seasonal_ratio_auto",
        "parameters" = c("a0_d", "a0_n", "a1", "a2", "a3", "a4", "a5")
      )
    )
  }
  
  # Check if parameter is a matrix - if not coerce
  if(!is.matrix(parameter_set)) {
    parameter_set <- as.matrix(parameter_set, ncol = 1)
  }
  
  # Parameters
  a0_d <- parameter_set[1, ]
  a0_n <- parameter_set[2, ]
  a1 <- parameter_set[3, ]
  a2 <- parameter_set[4, ]
  a3 <- parameter_set[5, ]
  a4 <- parameter_set[6, ]
  a5 <- parameter_set[7, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  seasonal_ratio <- catchment_data$seasonal_ratio
  CO2 <- catchment_data$CO2
  drought <- catchment_data$is_drought_year
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) #matrix(rep(precipitation, times = ncol(parameter_set)), ncol = ncol(parameter_set), byrow = FALSE)
  repeat_seasonal_ratio <- matrix(seasonal_ratio, ncol = ncol(parameter_set), nrow = length(seasonal_ratio), byrow = FALSE)
  repeat_CO2 <- matrix(CO2, ncol = ncol(parameter_set), nrow = length(CO2), byrow = FALSE) 
  repeat_drought <- matrix(drought, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE)
  repeat_a0_d <- matrix(a0_d, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a0_n <- matrix(a0_n, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a2 <- matrix(a2, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a3 <- matrix(a3, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a4 <- matrix(a4, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a5 <- matrix(a5, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Autocorrelation requirements
  boxcox_streamflow_average <- colMeans((repeat_a0_d * repeat_drought) + (repeat_a0_n * !repeat_drought) + (repeat_a1 * repeat_precipitation))
  previous_boxcox_streamflow <- boxcox_streamflow_average # assume the previous flows are the average. this means the first a2 component = 0.
  
  # CO2 requirements
  repeat_shifted_CO2 <- repeat_CO2 - repeat_a5
  repeat_shifted_CO2[repeat_shifted_CO2 < 0] <- 0 # max(0, CO2 - a5)
  
  # pre-allocate
  modelled_boxcox_streamflow <- matrix(numeric(length(repeat_precipitation)), ncol = ncol(parameter_set), nrow = length(precipitation))
  
  # Model: analysis is done row by row (timestep)
  for (time_step in seq_along(precipitation)) {
    current_boxcox_streamflow <- (repeat_a0_d[time_step, ] * repeat_drought[time_step, ]) + 
      (repeat_a0_n[time_step, ] * !repeat_drought[time_step, ]) + + 
      (repeat_a1[time_step, ] * repeat_precipitation[time_step, ]) + 
      (repeat_a2[time_step, ] * (previous_boxcox_streamflow - boxcox_streamflow_average)) +
      (repeat_a3[time_step, ] * repeat_shifted_CO2[time_step, ]) + 
      (repeat_a4[time_step, ] * repeat_seasonal_ratio[time_step, ])
    
    previous_boxcox_streamflow <- current_boxcox_streamflow
    modelled_boxcox_streamflow[time_step, ] <- current_boxcox_streamflow
  }
  
  return(modelled_boxcox_streamflow)
  
}



# 4. 
streamflow_model_separate_shifted_CO2_seasonal_ratio_auto <- function(catchment_data, parameter_set) {
  
  # If no parameters are given return description of model
  if(is.null(names(as.list(match.call())[-1]))) {
    return(
      list(
        "name" = "streamflow_model_separate_shifted_CO2_seasonal_ratio_auto",
        "parameters" = c("a0", "a1", "a2", "a3", "a4", "a5")
      )
    )
  }
  
  # Check if parameter is a matrix - if not coerce
  if(!is.matrix(parameter_set)) {
    parameter_set <- as.matrix(parameter_set, ncol = 1)
  }
  
  # Parameters
  a0 <- parameter_set[1, ]
  a1 <- parameter_set[2, ]
  a2 <- parameter_set[3, ]
  a3 <- parameter_set[4, ]
  a4 <- parameter_set[5, ]
  a5 <- parameter_set[6, ]
  
  # Get data
  precipitation <- catchment_data$precipitation
  seasonal_ratio <- catchment_data$seasonal_ratio
  CO2 <- catchment_data$CO2
  
  # Get into matrix form
  repeat_precipitation <- matrix(precipitation, ncol = ncol(parameter_set), nrow = length(precipitation), byrow = FALSE) #matrix(rep(precipitation, times = ncol(parameter_set)), ncol = ncol(parameter_set), byrow = FALSE)
  repeat_seasonal_ratio <- matrix(seasonal_ratio, ncol = ncol(parameter_set), nrow = length(seasonal_ratio), byrow = FALSE)
  repeat_CO2 <- matrix(CO2, ncol = ncol(parameter_set), nrow = length(CO2), byrow = FALSE) 
  repeat_a0 <- matrix(a0, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a1 <- matrix(a1, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a2 <- matrix(a2, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a3 <- matrix(a3, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a4 <- matrix(a4, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  repeat_a5 <- matrix(a5, ncol = ncol(parameter_set), nrow = nrow(repeat_precipitation), byrow = TRUE)
  
  # Autocorrelation requirements
  boxcox_streamflow_average <- colMeans(repeat_a0 + (repeat_a1 * repeat_precipitation))
  previous_boxcox_streamflow <- boxcox_streamflow_average # assume the previous flows are the average. this means the first a2 component = 0.
  
  # pre-allocate
  modelled_boxcox_streamflow <- matrix(numeric(length(repeat_precipitation)), ncol = ncol(parameter_set), nrow = length(precipitation))
  
  
  # CO2 requirements
  repeat_shifted_CO2 <- repeat_CO2 - repeat_a5
  repeat_shifted_CO2[repeat_shifted_CO2 < 0] <- 0 # max(0, CO2 - a5)
  
  # Model: analysis is done row by row (timestep)
  for (time_step in seq_along(precipitation)) {
    current_boxcox_streamflow <- repeat_a0[time_step, ] + 
      (repeat_a1[time_step, ] * repeat_precipitation[time_step, ]) + 
      (repeat_a2[time_step, ] * (previous_boxcox_streamflow - boxcox_streamflow_average)) +
      (repeat_a3[time_step, ] * repeat_shifted_CO2[time_step, ]) + 
      (repeat_a4[time_step, ] * repeat_seasonal_ratio[time_step, ])
    
    previous_boxcox_streamflow <- current_boxcox_streamflow
    modelled_boxcox_streamflow[time_step, ] <- current_boxcox_streamflow
  }
  
  return(modelled_boxcox_streamflow)
  
}





# Run CMAES --------------------------------------------------------------------
tic()

replicate_cmaes <- function(gauge, streamflow_model, objective_function) {
  cmaes_example <- gauge |>
    catchment_data_blueprint(
      observed_data = data,
      start_stop_indexes = start_stop_indexes
    ) |>
    numerical_optimiser_setup_vary_inputs(
      streamflow_model = {{ streamflow_model }},
      objective_function = {{ objective_function }},
      bounds_and_transform_method = make_default_bounds_and_transform_methods(),
      minimise_likelihood = TRUE
    ) |>
    my_cmaes(print_monitor = TRUE) |>
    result_set() 
}


REPLICATES <- 10


## Site 1: =====================================================================
site_1_gauge <- "113004A"

site_1_replicates <- replicate(
  n = REPLICATES, 
  expr = replicate_cmaes(
    site_1_gauge, 
    streamflow_model_separate_shifted_CO2_seasonal_ratio_auto, 
    constant_sd_objective_function
    ),
  simplify = FALSE
)  


best_site_1_parameters <- map(
  .x = site_1_replicates,
  .f = parameters_summary
) |>
  list_rbind() |>
  add_column(
    replicate = rep(seq(1, REPLICATES), each = length(site_1_replicates[[1]]$best_parameter_set)), # hard coded
    .before = 1
  ) |>
  dplyr::slice_min(
    loglikelihood
  )

# This is disgusting
# using the replicate number from best_site_1_parameters
site_1_streamflow <- modelled_streamflow_summary(site_1_replicates[[best_site_1_parameters$replicate[1]]])


## Site 2: =====================================================================
site_2_gauge <- "230210"


site_2_replicates <- replicate(
  n = REPLICATES, 
  expr = replicate_cmaes(
    site_2_gauge, 
    streamflow_model_drought_separate_shifted_CO2_seasonal_ratio_auto, 
    constant_sd_objective_function
  ),
  simplify = FALSE
)  


best_site_2_parameters <- map(
  .x = site_2_replicates,
  .f = parameters_summary
) |>
  list_rbind() |>
  add_column(
    replicate = rep(seq(1, REPLICATES), each = length(site_2_replicates[[1]]$best_parameter_set)), # hard coded
    .before = 1
  ) |>
  dplyr::slice_min(
    loglikelihood
  )

# This is disgusting
# using the replicate number from best_site_1_parameters
site_2_streamflow <- modelled_streamflow_summary(site_2_replicates[[best_site_2_parameters$replicate[1]]])



## Site 3: =====================================================================
site_3_gauge <- "302214"

site_3_replicates <- replicate(
  n = REPLICATES, 
  expr = replicate_cmaes(
    site_3_gauge, 
    streamflow_model_separate_shifted_CO2, 
    constant_sd_objective_function
  ),
  simplify = FALSE
) 


best_site_3_parameters <- map(
  .x = site_3_replicates,
  .f = parameters_summary
) |>
  list_rbind() |>
  add_column(
    replicate = rep(seq(1, REPLICATES), each = length(site_3_replicates[[1]]$best_parameter_set)), # hard coded
    .before = 1
  ) |>
  dplyr::slice_min(
    loglikelihood
  )

# This is disgusting
# using the replicate number from best_site_1_parameters
site_3_streamflow <- modelled_streamflow_summary(site_3_replicates[[best_site_3_parameters$replicate[1]]])



## Site 4: =====================================================================
site_4_gauge <- "407220"

site_4_replicates <- replicate(
  n = REPLICATES, 
  expr = replicate_cmaes(
    site_4_gauge, 
    streamflow_model_separate_shifted_CO2_seasonal_ratio_auto, 
    constant_sd_objective_function
  ),
  simplify = FALSE
) 


best_site_4_parameters <- map(
  .x = site_4_replicates,
  .f = parameters_summary
) |>
  list_rbind() |>
  add_column(
    replicate = rep(seq(1, REPLICATES), each = length(site_4_replicates[[1]]$best_parameter_set)), # hard coded
    .before = 1
  ) |>
  dplyr::slice_min(
    loglikelihood
  )

# This is disgusting
# using the replicate number from best_site_1_parameters
site_4_streamflow <- modelled_streamflow_summary(site_4_replicates[[best_site_4_parameters$replicate[1]]])



## Site 5: =====================================================================
site_5_gauge <- "614044"

site_5_replicates <- replicate(
  n = REPLICATES,
  expr = replicate_cmaes(
    site_5_gauge,
    streamflow_model_separate_shifted_CO2_seasonal_ratio_auto,
    constant_sd_objective_function
  ),
  simplify = FALSE
)


best_site_5_parameters <- map(
  .x = site_5_replicates,
  .f = parameters_summary
) |>
  list_rbind() |>
  add_column(
    replicate = rep(seq(1, REPLICATES), each = length(site_5_replicates[[1]]$best_parameter_set)), # hard coded
    .before = 1
  ) |>
  dplyr::slice_min(
    loglikelihood
  )

# This is disgusting
# using the replicate number from best_site_1_parameters
site_5_streamflow <- modelled_streamflow_summary(site_5_replicates[[best_site_5_parameters$replicate[1]]])

toc()



# Repeat the plotting in examine_five_sites ------------------------------------
combined_parameters <- rbind(best_site_1_parameters, best_site_2_parameters, best_site_3_parameters, best_site_4_parameters, best_site_5_parameters)
combined_streamflow <- rbind(site_1_streamflow, site_2_streamflow, site_3_streamflow, site_4_streamflow, site_5_streamflow)

write_csv(combined_streamflow, "./results/CMAES_results/testing_new_CO2_model_streamflow.csv")

streamflow_results <- combined_streamflow |> 
  left_join(
    gauge_information, 
    by = join_by(gauge)
    ) |> 
  pivot_longer(
    cols = ends_with("streamflow"),
    names_to = "modelled_or_observed",
    values_to = "boxcox_streamflow"
  ) |> 
  mutate(
    streamflow = boxcox_inverse_transform(boxcox_streamflow, lambda = bc_lambda),
    modelled_or_observed = if_else(modelled_or_observed == "observed_boxcox_streamflow", "observed", "modelled")
  )

## Plot 1
plot_streamflow_timeseries <- streamflow_results |>
  ggplot(aes(x = year, y = streamflow, colour = modelled_or_observed)) +
  geom_line(na.rm = TRUE, alpha = 0.9) +
  theme_bw() +
  scale_colour_brewer(palette = "Dark2") +
  labs(
    x = "Year",
    y = "Streamflow (mm)",
    title = "Streamflow timeseries"
  ) +
  facet_wrap(~gauge, scales = "free_y", nrow = 5) +
  theme(legend.title = element_blank())

plot_streamflow_timeseries


## Plot 2
difference_to_observed_streamflow <- streamflow_results |>
  select(!c(bc_lambda, boxcox_streamflow)) |>
  distinct() |>
  pivot_wider(
    names_from = modelled_or_observed,
    values_from = streamflow,
  ) |>
  mutate(
    observed_minus_CO2 = observed - modelled
  )


plot_difference_observed_residuals <- difference_to_observed_streamflow |>
  ggplot(aes(x = year, y = observed_minus_CO2)) +
  geom_line(na.rm = TRUE) +
  theme_bw() +
  labs(
    x = "Year",
    y = "Observed streamflow minus CO2 model streamflow (mm)",
    colour = "Residual Type",
    title = "Observed minus modelled streamflow residuals"
  ) +
  facet_wrap(~gauge, scales = "free_y", nrow = 5)


plot_difference_observed_residuals


# Plot 4
## acf using ggplot ============================================================
## https://stackoverflow.com/questions/17788859/acf-plot-with-ggplot2-setting-width-of-geom-bar

confidence_interval_acf <- function(acf_object, alpha = 0.05) {
  return(qnorm((1 + (1 - alpha)) / 2) / sqrt(acf_object$n.used))
}


tibble_confidence_interval_acf <- function(acf_object, name = NULL, alpha = 0.05) {
  upper_bound <- confidence_interval_acf(acf_object, alpha)
  lower_bound <- -upper_bound
  result_tibble <- as_tibble(cbind(name, upper_bound, lower_bound))
  result_tibble$upper_bound <- as.numeric(upper_bound)
  result_tibble$lower_bound <- as.numeric(lower_bound)
  return(result_tibble)
}


list_of_observed_minus_CO2_per_catchment <- difference_to_observed_streamflow |>
  select(year, gauge, observed_minus_CO2) |>
  pivot_wider(
    names_from = gauge,
    values_from = observed_minus_CO2
  ) |>
  select(!year) |>
  unclass()


remove_na_vector <- function(vector) {
  vector[!is.na(vector)]
}


no_NA_list_of_observed_minus_CO2_per_catchment <- map(
  .x = list_of_observed_minus_CO2_per_catchment,
  .f = remove_na_vector
)


acf_objects_per_gauge <- map(
  .x = no_NA_list_of_observed_minus_CO2_per_catchment,
  .f = acf,
  plot = FALSE
)


get_lags_from_acf <- function(acf_object, name = NULL) {
  acf_tibble <- with(acf_object, tibble(lag, acf))
  acf_tibble |>
    add_column(
      gauge = {{ name }},
      .before = 1
    )
}


lags_from_acf <- map2(
  .x = acf_objects_per_gauge,
  .y = names(acf_objects_per_gauge),
  .f = get_lags_from_acf
) |>
  list_rbind()


confidence_intervals_acf <- map2(
  .x = acf_objects_per_gauge,
  .y = names(acf_objects_per_gauge),
  .f = tibble_confidence_interval_acf
) |>
  list_rbind() |>
  rename(
    gauge = name
  )

complete_acf_per_gauge <- lags_from_acf |>
  left_join(
    confidence_intervals_acf,
    by = join_by(gauge)
  )




acf_ggplot <- complete_acf_per_gauge |>
  ggplot(aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0), linewidth = 1) +
  geom_segment(mapping = aes(xend = lag, yend = 0), linewidth = 1) +
  geom_hline(
    aes(yintercept = lower_bound),
    linetype = 2,
    linewidth = 1,
    colour = "blue"
  ) +
  geom_hline(
    aes(yintercept = upper_bound),
    linetype = 2,
    linewidth = 1,
    colour = "blue"
  ) +
  labs(
    x = "Lag",
    y = "ACF",
    title = "ACF graphs"
  ) +
  theme_bw() +
  facet_wrap(~gauge)

acf_ggplot



# CO2 contribution graphs for a5 -----------------------------------------------
combined_parameters <- read_csv("./Results/CMAES_results/testing_new_CO2_model_parameters.csv", show_col_types = FALSE)

CO2_contribution_streamflow <- function(a5, CO2) {
  
  altered_CO2_contribution <- CO2 - a5
  altered_CO2_contribution[altered_CO2_contribution < 0] <- 0
  
  return(altered_CO2_contribution)
}


shift_CO2_params <- combined_parameters |>
  filter(parameter == "a5") |> 
  select(gauge, parameter_value)


sample_CO2_data <- data |> 
  filter(gauge == "102101A") |> 
  pull(CO2)


CO2_contribution_timeseries <- map(
  .x = shift_CO2_params$parameter_value,
  .f = CO2_contribution_streamflow,
  CO2 = sample_CO2_data
)


names(CO2_contribution_timeseries) <- shift_CO2_params$gauge

tidy_CO2_contribution_timeseries <- as_tibble(CO2_contribution_timeseries) |> 
  pivot_longer(
    cols = everything(),
    names_to = "gauge",
    values_to = "CO2_contribution"
  ) |> 
  arrange(gauge) |> 
  mutate(
    year = row_number() + 1958, # quick and dirty way to get the year because everything starts at 1959
    .before = 1,
    .by = gauge
  )

# get year where CO2 starts to do something
# remove all zeroes by gauge
# slice_head by gauge
first_year_CO2_does_something <- tidy_CO2_contribution_timeseries |> 
  filter(CO2_contribution > 0) |> 
  slice_min(
    CO2_contribution,
    by = gauge
  )


tidy_CO2_contribution_timeseries |> 
  ggplot(aes(x = year, y = CO2_contribution)) +
  geom_line() +
  labs(
    x = "Year",
    y = "Contribution of CO2 - 280 to streamflow model",
    title = "When CO2 starts contributing to the model"
  ) +
  theme_bw() +
  facet_wrap(~gauge, scales = "free_y")


# Playing with CO2 data --------------------------------------------------------
sample_CO2_data <- data |> 
  filter(gauge == "102101A") |> 
  pull(CO2)

year <- data |> 
  filter(gauge == "102101A") |> 
  pull(year)


a3 <- 0.05
a4 <- 50 # a negative value shifts the impact of CO2 earlier in the streamflow model. Do I want this?

# a max(sample_CO2_data) effectively removes the CO2 components

base_CO2_equation <- a3 * sample_CO2_data 
matrix_test <- matrix(rep(sample_CO2_data - a4, times = 5), nrow = length(sample_CO2_data))
matrix_test[matrix_test <= 0] <- 0 
new_CO2_equation <- a3 * matrix_test



plot(x = year, y = base_CO2_equation, type = "l", ylim = range(c(base_CO2_equation, new_CO2_equation)), col = "blue")
lines(x = year, y = new_CO2_equation[, 1], col = "red")
lines(x = range(year), y = c(0, 0), lty = 2)

