# How did the fitting perform? Use plots to answer this question
pacman::p_load(tidyverse, gridExtra)


# Import functions -------------------------------------------------------------
source("./Functions/utility.R")


# Import the calibrated .csv's ---------------------------------------------------
streamflow_results <- read_csv("./Results/my_cmaes/CMAES_streamflow_results_20241028.csv", show_col_types = FALSE)
 

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



# There should only be 32? Not >5000 DELETE
#check_results_plot(streamflow_model = "streamflow_model_precip_only", 
#                   objective_function = "constant_sd_objective_function", 
#                   streamflow_results = setup_streamflow_results_plotting
#                   )

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
