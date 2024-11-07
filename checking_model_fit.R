# How did the fitting perform? Use plots to answer this question
pacman::p_load(tidyverse, gridExtra)


# Import functions -------------------------------------------------------------
source("./Functions/utility.R")


# Import the calibrated .csv's ---------------------------------------------------
streamflow_results <- read_csv("./Results/my_cmaes/CMAES_streamflow_results_20241107.csv", show_col_types = FALSE)
 

# Tidy data for plotting -------------------------------------------------------
setup_streamflow_results_plotting <- streamflow_results |>
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


# TESTING
# Determining a3 parameter bounds ----------------------------------------------
lambda <- seq(from = 0, to = 1.5, by = 0.01)

swap_boxcox_transform <- function(lambda, y) {
  boxcox_transform(y = y, lambda = lambda)
}

y <- map_dbl(
  .x = lambda,
  .f = swap_boxcox_transform,
  y = 100
)

plot(lambda, y)
points(0.856, boxcox_transform(100, 0.856), col = "red")

# look at gauge 312061 because it has the max lambda
# 313061 does hit the bounds for a3. Expand a3 parameter. Yes
# catchments
# I am not concerned about the sd, a5 and scale_CO2 being near zero.

parameter_results <- read_csv("./Results/my_cmaes/CMAES_parameter_results_20241028.csv", show_col_types = FALSE)
gauge_information <- read_csv("Data/Tidy/gauge_information_CAMELS.csv", show_col_types = FALSE)


check_bounds <- parameter_results |> 
  filter(near_bounds) |> 
  filter(parameter == "a3") #|> 
  pull(gauge) |> 
  unique()

check_lambda <- gauge_information |> 
  filter(gauge %in% check_bounds)

# Its not entirely dependent on lambda them

# Trial and error to get a3 good?
gauge <- "G9030124"

example <- gauge |> 
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  ) |> 
  numerical_optimiser_setup_vary_inputs(
    streamflow_model = streamflow_model_separate_shifted_CO2_seasonal_ratio,
    objective_function = constant_sd_objective_function,
    bounds_and_transform_method = make_default_bounds_and_transform_methods(),
    minimise_likelihood = TRUE
  ) |> 
  my_cmaes(print_monitor = TRUE) |> 
  result_set() 

parameters <- example |> 
  parameters_summary()