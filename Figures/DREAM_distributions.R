# Make DREAM distributions for supplemental
# Objective:
# - each row site
# - each column parameter
# - each element is a distributions
# - To get the best results group by model type on the same page

# Unknowns:
# - how many catchments per a4 page
# - how do I deal with the white space


# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, arrow, ggh4x)


# Group catchments by model type -----------------------------------------------
## Only plot same model types
best_CO2_non_CO2_per_gauge <- read_csv(
  "./Modelling/Results/CMAES/best_CO2_non_CO2_per_catchment_CMAES.csv",
  show_col_types = FALSE
)

best_model_per_gauge <- best_CO2_non_CO2_per_gauge |>
  slice_min(
    AIC,
    by = gauge
  ) |>
  select(gauge, streamflow_model) |> 
  distinct() |> 
  arrange(streamflow_model)


# Chunk gauges based on streamflow models --------------------------------------
## Get indexes
rle_models <- rle(best_model_per_gauge |> pull(streamflow_model))
models_with_drought <- str_detect(rle_models$values, "drought")
drought_index <- 1:sum(rle_models$lengths[models_with_drought])
non_drought_index <- (sum(rle_models$lengths[models_with_drought]) + 1):sum(rle_models$lengths)


# Split drought and non_drought_model_range evenly
gauges_per_graph <- 12L
num_split_drought <- ceiling(length(drought_index) / gauges_per_graph)
num_split_non_drought <- ceiling(length(non_drought_index) / gauges_per_graph)


split_factor_drought <- rep(
  seq(from = 1, to = num_split_drought), 
  each = gauges_per_graph
  )[drought_index]


split_factor_non_drought <- rep(
  seq(from = 1, to = num_split_non_drought), 
  each = gauges_per_graph
)[seq_along(non_drought_index)]



drought_gauges <- best_model_per_gauge$gauge[drought_index]
non_drought_gauges <- best_model_per_gauge$gauge[non_drought_index]  


drought_gauges_split <- split(drought_gauges, f = split_factor_drought)
non_drought_gauges_split <- split(non_drought_gauges, f = split_factor_non_drought)

all_gauges_split <- c(drought_gauges_split, non_drought_gauges_split) |> 
  # cannot repeat names otherwise map does not map
  `names<-`(seq(from = 1, to = length(c(drought_gauges_split, non_drought_gauges_split))))


# Extract DREAM and plot sequences ---------------------------------------------
DREAM_sequence_data <- open_dataset(
  source = "./Modelling/Results/DREAM/Sequences",
  format = "parquet"
) 

extract_and_plot_sequences <- function(gauge_vector, DREAM_data, identifier) {

  
  DREAM_sequence_data <- DREAM_data |> 
    filter(
      gauge %in% {{ gauge_vector }}
    ) |> 
    collect()
  
  
  
  # Plot 
  DREAM_distribution_plot <- DREAM_sequence_data |> 
    ggplot(aes(x = parameter_value)) +
    geom_histogram(bins = 30) +
    labs(
      x = "Parameter Value",
      y = "Frequency"
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    theme_bw() +
    ggh4x::facet_grid2(
      gauge ~ parameter, 
      scales = "free",
      independent = "x"
    ) +
    theme(
      text = element_text(size = 9),
      axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 6),
      strip.text = element_text(size = 6)
    )
  
  
  
  ggsave(
    filename = paste0("DREAM_distributions_", identifier, ".pdf"),
    path = "Figures/Supplementary",
    plot = DREAM_distribution_plot,
    device = "pdf",
    width = 183,
    height = 197,
    units = "mm"
  )
}


# Make figures -----------------------------------------------------------------
iwalk(
  .x = all_gauges_split,
  .f = extract_and_plot_sequences,
  DREAM_data = DREAM_sequence_data
)


# Make captions ----------------------------------------------------------------

create_caption <- function(gauge_vector, identifier) {
  
  # Convert gauge_vector into text for caption
  length_gauge_vector <- length(gauge_vector)
  middle_text <- paste0(
    ", ", 
    gauge_vector[2:(length_gauge_vector - 1)], 
    collapse = ""
    )
  gauge_text <- paste0(gauge_vector[1], middle_text, " and ", gauge_vector[length_gauge_vector])
  
  cat("\\begin{figure}") 
  cat("\n")
  cat("\t\\centering")
  cat("\n")
  cat(paste0("\t\\includegraphics[width=\\textwidth]{Figures/DREAM_distributions_", identifier, ".pdf}"))
  cat("\t\n")
  # The line below must change
  cat(paste0("\t\\caption{\\textbf{Parameter uncertainty for gauges ", gauge_text, "}. Same as fig.~\\ref{fig:supp_DREAM_distributions_1}}"))
  cat("\n")
  # The line below must change
  cat(paste0("\t\\label{fig:supp_DREAM_distributions_", identifier, "}")) 
  cat("\n")
  cat("\\end{figure}")
  cat("\n")
  cat("\n")
  
}

sink(file = "Figures/Supplementary/DREAM_distributions_captions_supp.txt") # filename must change
iwalk(
  .x = all_gauges_split,
  .f = create_caption
)
sink()


