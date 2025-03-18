# DREAM uncertainty timeseries

# OBJECTIVE:
# Recreate Vrugt streamflow timeseries plots using DREAM
# Try replicating the DREAM uncertainty plots using the rainfall-runoff relationship

# METHOD:
# - Get the DREAM results (sequences_xxxxx.csv)
# - Select one test catchment -> simplest = 141008A
# - Stick the sequence results into the streamflow_model
# - The best results are the solid line. The sequences can be the envelope

# Clear console ----------------------------------------------------------------
cat("\014") 

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, truncnorm, sloop, broom)

# Import functions -------------------------------------------------------------
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


# Import sequences -------------------------------------------------------------
sequences_DREAM <- read_csv(
  "Results/my_dream/sequences_20250304.csv",
  show_col_types = FALSE
  )

gauge <- "141008A"

single_catchment_sequence <- sequences_DREAM |> 
  filter(gauge == {{ gauge }})

# Import other converged data (DREAM and CMAES) --------------------------------
converged_stats_DREAM <- read_csv(
  "Results/my_dream/converged_stats_20250304.csv",
  show_col_types = FALSE
  ) |> 
  filter(gauge == {{ gauge }})

# The streamflow model use was streamflow_model_precip_auto


# Import CAMELS data -----------------------------------------------------------
CAMELS_data <- read_csv(
  "Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  show_col_types = FALSE
) |> 
  mutate(
    year = as.integer(year)
  )


start_end_index <- read_csv(
  "Data/Tidy/start_end_index.csv",
  show_col_types = FALSE
  )




# it is a continuous period
# I will have to account for start/stop with autocorrelation models
# if I want to generalise this


# Organising data to be used in the streamflow_model ---------------------------
# The data must be sorted row = parameters, columns = different parameter combinations

wide_single_catchment_sequence <- single_catchment_sequence |> 
  select(!gauge) |>
  unite(
    col = "identifier",
    c(chain, n),
    sep = "_" # need to be sep otherwise 1_11 = 11_1
  ) |> 
  pivot_wider(
    id_cols = parameter,
    values_from = parameter_value,
    names_from = identifier
  )


parameter_matrix <- wide_single_catchment_sequence |> 
  select(!parameter) |> 
  as.matrix() |> 
  unname()

# Apply streamflow model -------------------------------------------------------
gauge_data <- catchment_data_blueprint(
  gauge_ID = gauge,
  observed_data = CAMELS_data,
  start_stop_indexes = start_end_index
)

# You can't just put precip in. This is because catchment_data_blueprint
# has full_data_set and start_stop_data_set precipitation
# TEMPORARY solution is full dataset - does not account for autocorrelation 
# start/stop

streamflow <- streamflow_model_precip_auto(
  catchment_data = gauge_data,
  parameter_set = parameter_matrix 
) |> 
  `colnames<-`( # chapter 6 functions: advanced R 2nd
    paste0("col", seq_len(ncol(parameter_matrix)))
  ) |> 
  as_tibble() |> 
  add_column(
    year = CAMELS_data |> filter(gauge == {{ gauge }}) |> pull(year),
    .before = 1
  )

# Plotting
summary_streamflow <- streamflow |> 
  pivot_longer(
    cols = starts_with("col"),
    names_to = "chains",
    values_to = "streamflow"
  ) |> 
  summarise(
    percentile_95 = quantile(streamflow, 0.95),
    percentile_05 = quantile(streamflow, 0.05),
    upper = max(streamflow),
    lower = min(streamflow),
    .by = year
  )

# CMAES results?
CMAES_streamflow_results <- read_csv(
  "Results/my_cmaes/CMAES_streamflow_results_20250122.csv",
  show_col_types = FALSE
  ) |> 
  filter(gauge == {{ gauge }}) |>
  filter(streamflow_model == "streamflow_model_precip_auto") |> 
  select(c("year", "precipitation", "observed_boxcox_streamflow", "modelled_boxcox_streamflow"))


# 
final_summary_streamflow <- summary_streamflow |> 
  left_join(
    CMAES_streamflow_results,
    by = join_by(year)
  ) |> 
  relocate(
    c(
      "year", 
      "precipitation", 
      "modelled_boxcox_streamflow", 
      "percentile_05",
      "percentile_95",
      "observed_boxcox_streamflow"
      )
  ) |> 
  drop_na() |> 
  pivot_longer(
    cols = c(modelled_boxcox_streamflow, observed_boxcox_streamflow),
    names_to = "observed_or_modelled",
    values_to = "streamflow"
  ) |> 
  mutate(
    observed_or_modelled = if_else(str_detect(observed_or_modelled, "modelled"), "modelled streamflow using best parameter set", "observed streamflow")
  )

# transfering to real space flow?
bc_lambda <- read_csv(
  "Data/Tidy/gauge_information_CAMELS_with_climate.csv",
  show_col_types = FALSE
) |> 
  filter(gauge == {{ gauge }}) |> 
  pull(bc_lambda)

real_space_final_summary_streamflow <- final_summary_streamflow |> 
  select(!c(year, precipitation, observed_or_modelled)) |> 
  apply(
    MARGIN = c(1, 2), 
    FUN = boxcox_inverse_transform,
    lambda = bc_lambda,
    lambda_2 = 1
    ) |> 
  as_tibble() |> 
  add_column(
    observed_or_modelled = final_summary_streamflow |> pull(observed_or_modelled),
    .before = 1
  ) |> 
  add_column(
    precipitation = final_summary_streamflow |> pull(precipitation),
    .before = 1
  ) |> 
  add_column(
    year = final_summary_streamflow |> pull(year),
    .before = 1
  ) 
  

# boxcox plot
boxcox_timeseries <- final_summary_streamflow |> 
  ggplot(aes(x = year, y = streamflow, colour = observed_or_modelled)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.7,
    colour = NA,
    ) +
  geom_line() +
  labs(
    x = "Year",
    y = "Box-Cox Streamflow (mm)",
    colour = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.93),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.2)
  ) 
  

boxcox_timeseries
# realspace plot
realspace_timeseries <- real_space_final_summary_streamflow |> 
  ggplot(aes(x = year, y = streamflow, colour = observed_or_modelled)) +
  geom_ribbon(
    aes(ymin = percentile_05, ymax = percentile_95),
    alpha = 0.7,
    colour = NA,
  ) +
  geom_line() +
  labs(
    x = "Year",
    y = "Streamflow (mm)",
    colour = NULL,
    title = paste0("Parameter Uncertainty (5th and 95th percentiles) - Gauge: ", gauge)
  ) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.93),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.2),
    plot.title = element_text(hjust = 0.5)
  )

realspace_timeseries


# Repeat uncertainty using the rainfall-runoff relationship --------------------

# Uncertainty around the line of best using the parameters?
# I have the boxcox flow of all the different parameter sets
# Run lm(boxcox~rainfall) for all parameter sets -> get the fitted values along 
# the line
# Take the percentiles or min/max

precip_year <- CAMELS_data |> 
  filter(gauge == {{ gauge }}) |> 
  filter(included_in_calibration) |> 
  select(year, p_mm)

rainfall_streamflow_sequences <- streamflow |> 
  pivot_longer(
    cols = starts_with("col"),
    names_to = "identifier",
    values_to = "streamflow"
  ) |> 
  left_join(
    precip_year,
    by = join_by(year)
  ) |> 
  arrange(identifier) |> 
  rename(
    precipitation = p_mm
  ) |> 
  drop_na()

lm_fitted_values <- function(x, y) {
  lm(y ~ x)$fitted.values
}


fitted_values_sequences <- rainfall_streamflow_sequences |> 
  mutate(
    fitted_values = lm_fitted_values(x = precipitation, y = streamflow),
    .by = identifier
  )

joining_final_summary_streamflow <- final_summary_streamflow |> 
  select(year, observed_or_modelled, streamflow)

summary_fitted_values <- fitted_values_sequences |> 
  summarise(
    # year, precip and streamflow don't change anything. 
    # It means I don't have to join the values again
    year = max(year), 
    precipitation = max(precipitation),
    upper = max(fitted_values),
    lower = min(fitted_values),
    percentile_95 = quantile(fitted_values, 0.975),
    percentile_05 = quantile(fitted_values, 0.025),
    .by = year
  ) |> 
  right_join(
    joining_final_summary_streamflow,
    by = join_by(year)
  )
  
  


summary_fitted_values |> 
  mutate(
    observed_or_modelled = factor(
      observed_or_modelled, 
      levels = c("observed streamflow", "modelled streamflow using best parameter set")
      )
  ) |> 
  ggplot(aes(x = precipitation, y = streamflow, colour = observed_or_modelled, fill = observed_or_modelled)) +
  geom_ribbon( # parameter uncertainty
    aes(ymin = percentile_05, ymax = percentile_95),
    alpha = 0.8,
    colour = NA,
    fill = "black"
  ) +
  geom_smooth(
    formula = y ~ x,
    method = lm,
    se = TRUE, # standard error
    alpha = 0.2
  ) +
  geom_point() +
  labs(
    x = "Annual Precipitation (mm)",
    y = "Annual Box-Cox Streamflow (mm)",
    title = paste0("Rainfall-runoff parameter uncertainty DREAM - Gauge: ", gauge),
    colour = NULL,
    fill = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.15, 0.93),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.2),
    plot.title = element_text(hjust = 0.5)
  )


# Statistical tests (t-test on a3 parameters) ----------------------------------
# Lets try a single catchment
# extract all a3 values
# t.test() see if it is different from zero

# PROBLEM
# Not all sequences from catchments are saved...
# Does .csv have a maximum line capacity - max row limit?
# Something wrong with 04 how I assign gauges?

x <- sequences_DREAM |> 
  pull(gauge) |> 
  unique() 
# I am missing half of the gauges.

only_a3_sequences <- sequences_DREAM |>
  filter(parameter == "a3") 

# T-test greater than or less than zero function
# One sided t-test used to compare is statistically different from zero
# We have an idea of the direction of change based on the histograms

one_sided_t_test_a3 <- function(gauge, only_a3_sequences) {
  
  a3_sequences_per_gauge <- only_a3_sequences |> 
    filter(gauge == {{ gauge }}) |> 
    pull(parameter_value)
  
  # Determine less or greater test
  a3_mean_per_gauge <- mean(a3_sequences_per_gauge)
  test_type <- if_else(sign(a3_mean_per_gauge) == -1, "less", "greater")
  
  a3_t_test_per_gauge <- a3_sequences_per_gauge |> 
    t.test(
      alternative = test_type,
      mu = 0
    ) |> 
    tidy() |> 
    add_column(
      gauge = {{ gauge }},
      .before = 1
    )
  
}

a3_t_test_results <- map(
  .x = only_a3_sequences |> pull(gauge) |> unique(),
  .f = one_sided_t_test_a3,
  only_a3_sequences = only_a3_sequences
) |> 
  list_rbind()



