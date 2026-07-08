# Illustration of rainfall-runoff changes

# Figures produced in this R file ----------------------------------------------
# - Figures/Main/illustration.pdf


# CODE


# Import libraries -------------------------------------------------------------
pacman::p_load(tidyverse, sloop, patchwork)


# Import functions -------------------------------------------------------------
source("./Functions/utility.R")
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
source("./Functions/boxcox_logsinh_transforms.R")


# Import data ------------------------------------------------------------------
data <- readr::read_csv(
  "./Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  show_col_types = FALSE
) |>
  mutate(year = as.integer(year)) |>
  # required for log-sinh. Log-sinh current formulation has asymptote of zero.
  # This means zero flows of ephemeral catchments cannot be transformed
  # add a really small value
  mutate(q_mm = q_mm + .Machine$double.eps^0.5)


gauge_information <- readr::read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE
)


best_CO2_and_non_CO2_model_and_params_per_gauge <- read_csv(
  "./Modelling/Results/CMAES/best_CO2_non_CO2_per_catchment_CMAES.csv",
  show_col_types = FALSE
)


start_stop_indexes <- read_csv(
  "Data/Tidy/start_end_index.csv",
  show_col_types = FALSE
)

# Calculate evidence ratio for filtering ---------------------------------------
high_evi_gauges <- best_CO2_and_non_CO2_model_and_params_per_gauge |>
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
  arrange(evidence_ratio) |>
  filter(evidence_ratio > 100) |>
  pull(gauge)


# Find gauge to use as an example ----------------------------------------------
## Get gauges with high record lengths
high_record_gauges <- gauge_information |>
  filter(record_length > 50) |>
  pull(gauge)

## Try using the simplest CO2 models first
best_CO2_model_gauge <- best_CO2_and_non_CO2_model_and_params_per_gauge |>
  filter(contains_CO2) |> 
  filter(gauge %in% high_evi_gauges) |>
  #filter(gauge %in% high_record_gauges) |>
  select(gauge, streamflow_model) |>
  distinct() |>
  arrange(desc(streamflow_model))


best_CO2_model_gauge |>
  filter(streamflow_model == "streamflow_model_intercept_shifted_CO2")



rearrange_catchment_data_blueprint <- function(observed_data, gauge_ID, start_stop_indexes) {
  catchment_data_blueprint(
    gauge_ID = gauge_ID,
    observed_data = observed_data,
    start_stop_indexes = start_stop_indexes
  )
}

replace_CO2_in_data <- function(new_CO2, data) {
  data |>
    mutate(
      CO2 = new_CO2
    )
}


illustration_plots <- function(gauge, plot_label, plot_arrow) {
  b <- best_CO2_and_non_CO2_model_and_params_per_gauge |>
    filter(contains_CO2) |>
    filter(gauge == {{ gauge }}) |>
    filter(parameter == "b") |>
    pull(parameter_value)


  modified_data <- data |>
    filter(gauge == {{ gauge }}) |>
    mutate(
      q_log_sinh = log_sinh_transform(b = b, y = q_mm, offset = 0)
    )

  
  # Year range
  year_range <- modified_data |> 
    drop_na() |> 
    pull(year) |> 
    range()

  # Plotting lines for the illustration
  CO2_for_selected_years <- modified_data |>
    filter(year %in% c(1959, 1979, 1999, 2019)) |>
    pull(CO2)


  different_CO2_modified_data <- map(
    .x = CO2_for_selected_years,
    .f = replace_CO2_in_data,
    data = modified_data
  )


  different_CO2_catchment_data <- map(
    .x = different_CO2_modified_data,
    .f = rearrange_catchment_data_blueprint,
    gauge = gauge,
    start_stop_indexes = start_stop_indexes
  )

  # if it is a drought model turn off the drought intercept by setting it
  # to the non-drought intercept
  non_drought_intercept <- best_CO2_and_non_CO2_model_and_params_per_gauge |>
    filter(contains_CO2) |>
    filter(gauge == {{ gauge }}) |>
    filter(parameter == "a0_n") |> 
    pull(parameter_value)

  # empty vectors do not play nice
  if(is_empty(non_drought_intercept)){non_drought_intercept <- NA}
  
  parameter_set <- best_CO2_and_non_CO2_model_and_params_per_gauge |>
    filter(contains_CO2) |>
    filter(gauge == {{ gauge }}) |>
    # make sure we get a straight line by turning off a2, a4 etc.
    mutate(
      parameter_value = case_when(
        parameter == "a0_d" ~ non_drought_intercept,
        parameter == "a2" ~ 0,
        parameter == "a4" ~ 0,
        .default = parameter_value
      )
    ) |>
    pull(parameter_value)

  
  
  
  streamflow_model <- best_CO2_model_gauge |>
    filter(gauge == {{ gauge }}) |>
    pull(streamflow_model) |>
    match.fun()


  different_CO2_streamflow <- map(
    .x = different_CO2_catchment_data,
    .f = streamflow_model,
    parameter = parameter_set
  )


  streamflow_1 <- different_CO2_streamflow[[1]] |>
    select(year, precipitation, streamflow_results) |>
    rename(streamflow_1960 = streamflow_results)

  streamflow_2 <- different_CO2_streamflow[[2]] |>
    select(year, precipitation, streamflow_results) |>
    rename(streamflow_1990 = streamflow_results)

  streamflow_3 <- different_CO2_streamflow[[3]] |>
    select(year, precipitation, streamflow_results) |>
    rename(streamflow_2010 = streamflow_results)
  
  streamflow_4 <- different_CO2_streamflow[[4]] |>
    select(year, precipitation, streamflow_results) |>
    rename(streamflow_2020 = streamflow_results)
  
  # if there is no data beyond 2019 for a gauge do not plot it
  if(max(year_range) < 2019) {
    streamflow_4 <- streamflow_4 |> 
      mutate(
        streamflow_2020 = NA
      )
  }

  rainfall_runoff_relationship <- streamflow_1 |>
    left_join(
      streamflow_2,
      by = join_by(year, precipitation)
    ) |>
    left_join(
      streamflow_3,
      by = join_by(year, precipitation)
    ) |>
    left_join(
      streamflow_4,
      by = join_by(year, precipitation)
    ) |>
    pivot_longer(
      cols = starts_with("streamflow"),
      names_to = "rainfall_runoff_year",
      values_to = "streamflow"
    ) |> # rename streamflow_1960
    mutate(
      rainfall_runoff_year = case_when(
        rainfall_runoff_year == "streamflow_1960" ~ "Estimated Relationship 1959",
        rainfall_runoff_year == "streamflow_1990" ~ "Estimated Relationship 1979",
        rainfall_runoff_year == "streamflow_2010" ~ "Estimated Relationship 1999",
        rainfall_runoff_year == "streamflow_2020" ~ "Estimated Relationship 2019",
        .default = NA
      )
    )



  # Arrow aes is hard coded based on the main figure plot
  # I could rescale it using max and min log-sinh flow and precip?
  if (plot_label == "B") {
    #y_axis_exclude <- element_blank()
    arrow_aes <- aes(x = 980, y = 180, xend = 980, yend = 60)
  } else {
    #y_axis_exclude <- element_text()
    arrow_aes <- aes(x = 1340, y = 350, xend = 1340, yend = 120)
  }
  
  
  
  # put it together now
  plot <- modified_data |>
    drop_na() |> 
    ggplot(aes(x = p_mm, y = q_log_sinh, fill = year)) +
    geom_line(
      data = rainfall_runoff_relationship,
      aes(x = precipitation, y = streamflow, colour = rainfall_runoff_year),
      inherit.aes = FALSE
    ) +
    geom_point(
      colour = "black",
      shape = 21,
      size = 3
    ) +
    labs(
      x = "Annual Precipitation (mm)",
      y = "Log-sinh Annual Streamflow",
      fill = bquote(atop(CO[2]~"ppm", "(Rainfall-Runoff Year)")), #"atop("Observed ", "("*CO[2]~"ppm)")
      colour = "Modelled Rainfall-Runoff Relationship",
      title = plot_label
    ) +
    scale_colour_manual(values = c("#440154FF", "#33638DFF", "#55C667FF", "#B8DE29FF")) +
    # labels hard coded found using data |> filter(year %in% c(1960, 1980, 2000, 2020)) |> pull(CO2) |> unique()
    scale_fill_continuous(
      limits = c(1959, 2022), # hard code limits based on years (min and max of entire data)
      palette = "viridis", # options: viridis, plasma, RdYlBu
      breaks = c(1960, 1980, 2000, 2020), 
      labels = paste0(c(317, 339, 370, 414), "\n", c("(1960)", "(1980)", "(2000)", "(2020)"))
      ) + 
    theme_bw() +
    theme(
      legend.title.position = "top",
      legend.title = element_text(hjust = 0.5),
      plot.title = element_text(vjust = -8, hjust = 0.05, face = "bold"),
      #axis.title.y = y_axis_exclude,
      text = element_text(size = 9)
    ) +
    guides(
      fill = guide_colourbar(
        barwidth = unit(6, "cm")
      ),
      colour = guide_legend(
        nrow = 2,
        ncol = 2,
        override.aes = aes(linewidth = 1)
        )
    )

  
  # add arrow is true
  if(plot_arrow) {
    plot <- plot + 
      geom_segment(
        mapping = arrow_aes,
        arrow = arrow(type = "closed", length = unit(2, "mm"))
      ) 
  }
  
  return(plot)
}


# plots - pick two gauges to illustration (intercept, slope)
part_a <- illustration_plots(gauge = "606195", plot_label = "A", plot_arrow = TRUE) 
part_b <- illustration_plots(gauge = "238235", plot_label = "B", plot_arrow = TRUE) + 
  theme(axis.title.y = element_blank())
illustration <- part_a | part_b 
final_illustration <- illustration + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(
  filename = "illustration.pdf",
  plot = final_illustration,
  device = "pdf",
  path = "Figures/Main",
  width = 180,
  height = 130,
  units = "mm"
)


# add % change for the paragraph -----------------------------------------------
## mean annual rainfall ========================================================
mean_annual_rainfall <- data |> 
  filter(gauge %in% c("606195", "238235")) |> 
  summarise(
    mean_annual_rainfall = mean(p_mm),
    max_annual_rainfall = max(p_mm),
    median_annual_rainfall = median(p_mm),
    .by = gauge
  ) 

## use mean annual rainfall to find streamflow % change ========================
### modify illustration plots function to find % change
### What this code does:
### - makes log-sinh runoff vs. rainfall lines using the simplest CO2 model (just slope or intercept)
### - get the log-sinh runoff for the average rainfall
### - convert the log-sinh runoff back into realspace

percentage_change_based_on_rainfall <- function(gauge, mean_annual_rainfall_data) {
  b <- best_CO2_and_non_CO2_model_and_params_per_gauge |>
    filter(contains_CO2) |>
    filter(gauge == {{ gauge }}) |>
    filter(parameter == "b") |>
    pull(parameter_value)
  
  
  mean_annual_rainfall <- mean_annual_rainfall_data |> 
    filter(gauge == {{ gauge }}) |> 
    pull(mean_annual_rainfall)
    
  
  modified_data <- data |>
    filter(gauge == {{ gauge }}) |>
    mutate(
      q_log_sinh = log_sinh_transform(b = b, y = q_mm, offset = 0)
    ) |> 
    mutate(
      p_mm = mean_annual_rainfall
    )
  
  
  # Plotting lines for the illustration
  CO2_for_selected_years <- modified_data |>
    filter(year %in% c(1959, 1979, 1999, 2019)) |>
    pull(CO2)
  
  
  different_CO2_modified_data <- map(
    .x = CO2_for_selected_years,
    .f = replace_CO2_in_data,
    data = modified_data
  )
  
  
  different_CO2_catchment_data <- map(
    .x = different_CO2_modified_data,
    .f = rearrange_catchment_data_blueprint,
    gauge = gauge,
    start_stop_indexes = start_stop_indexes
  )
  
  
  parameter_set <- best_CO2_and_non_CO2_model_and_params_per_gauge |>
    filter(contains_CO2) |>
    filter(gauge == {{ gauge }}) |>
    # make sure we get a straight line by turning off a2, a4 etc.
    mutate(
      parameter_value = case_when(
        parameter == "a2" ~ 0,
        parameter == "a4" ~ 0,
        .default = parameter_value
      )
    ) |>
    pull(parameter_value)
  
  
  streamflow_model <- best_CO2_model_gauge |>
    filter(gauge == {{ gauge }}) |>
    pull(streamflow_model) |>
    match.fun()
  
  
  different_CO2_streamflow <- map(
    .x = different_CO2_catchment_data,
    .f = streamflow_model,
    parameter = parameter_set
  ) |> 
    # extract streamflow for each of the results
    list_rbind() |> 
    select(precipitation, CO2, streamflow_results) |> 
    distinct() |> 
    # transform streamflow results into realspace
    mutate(
      realspace_streamflow = inverse_log_sinh_transform(b = b, z = streamflow_results, offset = 0)
    )
  
}

model_fit_gauge_A <- percentage_change_based_on_rainfall(gauge = "606195", mean_annual_rainfall_data = mean_annual_rainfall)
(model_fit_gauge_A$realspace_streamflow[1] - model_fit_gauge_A$realspace_streamflow[4]) / model_fit_gauge_A$realspace_streamflow[1]

model_fit_gauge_B <- percentage_change_based_on_rainfall(gauge = "238235", mean_annual_rainfall_data = mean_annual_rainfall)
(model_fit_gauge_B$realspace_streamflow[1] - model_fit_gauge_B$realspace_streamflow[4]) / model_fit_gauge_B$realspace_streamflow[1]




# Repeat plot for all 81 high evidence ratio catchments ------------------------
## Remove gauges already in chapter 4 ==========================================
supp_illustration_gauges <- high_evi_gauges[!high_evi_gauges %in% c("606195", "238235")] 

## I want 8 per page - split supp_illustration_gauges ==========================
chunk <- 8
n <- length(supp_illustration_gauges)
split_group <- rep(rep(1:ceiling(n / chunk), each = chunk))[1:n]
split_supp_illustration_gauges <- split(supp_illustration_gauges, split_group)

## Plot groups =================================================================
### this is a wrapper for illustration_plots
grouped_illustration_plot <- function(split_gauges) {
  map2(
    .x = split_gauges,
    .y = letters[1:length(split_gauges)],
    .f = illustration_plots,
    plot_arrow = FALSE
  ) |> 
    reduce(
      .f = `+`
    ) + 
    plot_layout(
      nrow = 4, 
      guides = "collect",
      axis_titles = "collect"
    ) & 
    theme(
      legend.position = "bottom",
      plot.margin = margin(t = 0, l = 10, r = 10, b = 0, unit = "pt")
    )
}

### Make caption for plots
#### use the same style as the extended figures
#### something for gauges 123456 (a), 123456 (b) etc. Same as fig. X.
create_caption <- function(gauge_chunk, identifier) {
  abc <- letters[1:length(gauge_chunk)]
  gauge_abc <- paste0(gauge_chunk, " (", abc, ")")
  # concatenate everything but last value
  start_gauge_abc <- paste0(gauge_abc[1:(length(gauge_abc) - 2)], ", ", collapse = "")
  end_gauge_abc <- paste0(gauge_abc[(length(gauge_abc) - 1)], " and ", gauge_abc[length(gauge_abc)], ".")
  gauge_text <- paste(c(start_gauge_abc, end_gauge_abc), collapse = "")
  
  cat("\\begin{figure}")
  cat("\n")
  cat("\t\\centering")
  cat("\n")
  cat(paste0("\t\\includegraphics[width=\\textwidth]{Figures/illustration_plot_", identifier, ".pdf}"))
  cat("\t\n")
  # The line below must change
  cat(paste0("\t\\caption{\\textbf{Observed shifts in the rainfall-runoff relationship for gauges ", gauge_text, "} Same as fig X.}"))
  cat("\n")
  # The line below must change
  cat(paste0("\t\\label{fig:supp_illustration_", identifier, "}"))
  cat("\n")
  cat("\\end{figure}")
  cat("\n")
  cat("\\clearpage")
  cat("\n")
  cat("\n")
}

sink("Figures/Other/supp_illustration_gauge_caption.txt")
iwalk(
  .x = split_supp_illustration_gauges,
  .f = create_caption
)
sink()


# Issues with plot:
# - not all gauges have data at 1959 thus missing - leave
# - arrow is in a bad position (likely add argument to turn this on/off) or fix it

#illustration_plots(gauge = supp_illustration_gauges[8], plot_label = "h", plot_arrow = F)
supp_illustration_plots <- map(
  .x = split_supp_illustration_gauges,
  .f = grouped_illustration_plot
)


supp_illustration_plot_names <- paste0("illustration_plot_", 1:length(supp_illustration_plots), ".pdf")

walk2(
  .x = supp_illustration_plot_names,
  .y = supp_illustration_plots,
  .f = ggsave,
  path = "Figures/Other/",
  device = "pdf",
  width = 180,
  height = 254,
  units = "mm"
)


# Using observed data does a catchment with the same amount of rainfall receive less runoff? -----

## Aggregate high evidence ratio gauges by decade ==============================
aggregated_runoff_ratio <- data |> 
  filter(gauge %in% high_evi_gauges) |> 
  mutate(
    decade = year - (year %% 10)
  ) |> 
  summarise(
    # mean account for missing year sum(q) / number of observations not NA
    decade_mean_q = mean(q_mm, na.rm = TRUE),
    #use all available rainfall data
    decade_mean_p = mean(p_mm, na.rm = TRUE),
    .by = decade
  ) |> 
  mutate(
    runoff_ratio = decade_mean_q / decade_mean_p
  )

plot_aggregated_runoff_ratio <- aggregated_runoff_ratio |> 
  drop_na() |> 
  filter(decade != 1950) |> 
  ggplot(aes(x = decade, y = runoff_ratio,)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Decade",
    y = "Runoff Ratio"
  ) +
  theme_bw()

#ggsave()

## Aggreate by gauge ===========================================================
by_gauge_runoff_ratio <- data |> 
  filter(gauge %in% high_evi_gauges) |> 
  mutate(
    decade = year - (year %% 10)
  ) |> 
  summarise(
    # mean account for missing year sum(q) / number of observations not NA
    decade_mean_q = mean(q_mm, na.rm = TRUE),
    #use all available rainfall data
    decade_mean_p = mean(p_mm, na.rm = TRUE),
    .by = c(decade, gauge)
  ) |> 
  mutate(
    runoff_ratio = decade_mean_q / decade_mean_p
  ) |> 
  # there is only one data point for the 1950s remove
  filter(decade != 1950) 

### Mean by gauge ##############################################################
mean_decade_by_gauge <- by_gauge_runoff_ratio |> 
  summarise(mean_ratio = mean(runoff_ratio, na.rm = T), .by = decade) |> 
  ggplot(aes(x = decade, y = mean_ratio)) + 
  geom_line() +
  geom_point() +
  labs(
    x = "Decade",
    y = "Runoff Ratio"
  ) +
  theme_bw()

ggsave(
  filename = "Figures/Other/change_in_runoff_ratio.pdf",
  plot = mean_decade_by_gauge,
  device = "pdf",
  width = 140,
  height = 100,
  units = "mm"
)
 

all_gauge_runoff_ratio <- by_gauge_runoff_ratio |> 
  drop_na() |> 
  ggplot(aes(x = decade, y = runoff_ratio)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Decade",
    y = "Runoff Ratio"
  ) +
  theme_bw() +
  facet_wrap(~gauge, scale = "free_y")

ggsave(
  filename = "Figures/Other/mega_runoff_ratio_plot.pdf",
  plot = all_gauge_runoff_ratio,
  device = "pdf",
  width = 1189, 
  height = 841,
  units = "mm"
)

## Aggregate by state ==========================================================
by_state_runoff_ratio <- data |> 
  left_join(
    gauge_information,
    by = join_by(gauge)
    ) |> 
  filter(gauge %in% high_evi_gauges) |> 
  mutate(
    decade = year - (year %% 10)
  ) |> 
  summarise(
    # mean account for missing year sum(q) / number of observations not NA
    decade_mean_q = mean(q_mm, na.rm = TRUE),
    #use all available rainfall data
    decade_mean_p = mean(p_mm, na.rm = TRUE),
    .by = c(decade, state)
  ) |> 
  mutate(
    runoff_ratio = decade_mean_q / decade_mean_p
  )
  

# number of gauges by state
gauge_information |> 
  filter(gauge %in% high_evi_gauges) |> 
  count(state)

by_state_runoff_ratio |> 
  drop_na() |> 
  filter(decade != 1950) |> 
  filter(state %in% c("NSW", "TAS", "VIC", "WA")) |> # only 1 ACT, 1 QLD and 1 SA gauge
  ggplot(aes(x = decade, y = runoff_ratio, colour = state)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Decade",
    y = "Runoff Ratio"
  ) +
  theme_bw()








## Show a change in runoff ratio for Figure 4 in paper =========================
fig_4_specific_aggregated_runoff_ratio <- data |> 
  filter(gauge %in% high_evi_gauges) |> 
  mutate(
    decade = case_when(
      year %in% 1990:1999 ~ "1990-1999",
      year %in% 2012:2021 ~ "2012-2021",
      .default = NA
    )
  ) |>  
  summarise(
    # mean account for missing year sum(q) / number of observations not NA
    decade_mean_q = mean(q_mm, na.rm = TRUE),
    #use all available rainfall data
    decade_mean_p = mean(p_mm, na.rm = TRUE),
    .by = c(decade, gauge)
  ) |> 
  mutate(
    runoff_ratio = decade_mean_q / decade_mean_p
  )
