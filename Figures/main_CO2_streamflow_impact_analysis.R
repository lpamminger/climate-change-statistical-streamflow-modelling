# Analysis catchments with high evidence ratio 


# Figures produced in this R file ----------------------------------------------
# 1. streamflow_percentage_difference_with_timeseries.pdf
# 2. streamflow_timeseries_extended_data_plot.pdf (extended figure t-series)



# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, truncnorm, sloop, patchwork, ozmaps, sf, patchwork, metR, ggmagnify)


# Import functions -------------------------------------------------------------
source("./Functions/utility.R")


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


DREAM_CO2_impact_uncertainty_on_streamflow <- read_csv(
  "Modelling/Results/DREAM/DREAM_CO2_impact_uncertainty_on_streamflow.csv",
  show_col_types = FALSE
)


evidence_ratio <- read_csv(
  "./Modelling/Results/CMAES/evidence_ratio_results.csv",
  show_col_types = FALSE
)


best_CO2_non_CO2_per_gauge <- read_csv(
  "./Modelling/Results/CMAES/best_CO2_non_CO2_per_catchment_CMAES.csv",
  show_col_types = FALSE
)

master_streamflow_table <- read_csv(
  "Modelling/Results/master_streamflow_table.csv",
  show_col_types = FALSE
)




# Only include gauges with evidence ratio > 100 --------------------------------
high_evidence_ratio_gauges <- evidence_ratio |> 
  filter(evidence_ratio > 100) |> 
  pull(gauge)

high_evidence_master_streamflow_table <- master_streamflow_table |> 
  filter(gauge %in% high_evidence_ratio_gauges)


# Calculate percentage difference between streamflow CO2 on and off ------------
### Method:
### - select two 10 year periods
### - sum the modelled streamflow with CO2 off (natural) and CO2 on (anthropogenic)
###   over the years during the 10 year periods
### - find the difference in the two 10 year periods
### - average by number of years during decade
### - percentage change is ((CO2_on - CO2_off) / CO2_off) * 100

decade_1 <- seq(from = 1990, to = 1999)
decade_2 <- seq(from = 2012, to = 2021)


percentage_difference_a3_on_off_data <- high_evidence_master_streamflow_table |>
  filter(year %in% c(decade_1, decade_2)) |>
  # add decade group for summarising
  mutate(
    decade = case_when( # year - (year %% 10)
      year %in% decade_1 ~ 1,
      year %in% decade_2 ~ 2,
      .default = NA
    )
  ) |>
  filter(!is.na(decade)) |>
  # sum streamflow for each decade
  summarise(
    sum_decade_realspace_CO2_off_streamflow = sum(realspace_streamflow_CO2_model_off),
    sum_decade_realspace_CO2_on_streamflow = sum(realspace_streamflow_CO2_model_on),
    sum_decade_realspace_observed_streamflow = sum(realspace_observed_streamflow),
    years_of_data = n(),
    .by = c(gauge, decade)
  ) |>
  # find the absolution and percentage difference
  mutate(
    realspace_CO2_off_streamflow_per_year = sum_decade_realspace_CO2_off_streamflow / years_of_data,
    realspace_a3_on_streamflow_per_year = sum_decade_realspace_CO2_on_streamflow / years_of_data,
    CO2_impact_on_streamflow_mm_per_year = (realspace_a3_on_streamflow_per_year - realspace_CO2_off_streamflow_per_year),
    CO2_impact_on_streamflow_percent = (CO2_impact_on_streamflow_mm_per_year / realspace_CO2_off_streamflow_per_year) * 100
  ) |>
  arrange(desc(CO2_impact_on_streamflow_percent)) # Large percentage changes are not tied to years_of_data



## Add other data to percentage difference for plotting ==========================
plot_ready_decade_differences <- percentage_difference_a3_on_off_data |>
  left_join(
    evidence_ratio,
    by = join_by(gauge)
  ) |>
  left_join(
    DREAM_CO2_impact_uncertainty_on_streamflow,
    by = join_by(gauge, decade)
  )




# Compare DREAM values to CMAES values -----------------------------------------
compare_CMAES_and_DREAM <- plot_ready_decade_differences |>
  select(gauge, decade, CO2_impact_on_streamflow_percent, IQR_CO2_impact_on_streamflow_percentage, median_CO2_impact_on_streamflow_percentage) |>
  drop_na() |>
  mutate(
    CMAES_DREAM_diff = CO2_impact_on_streamflow_percent - median_CO2_impact_on_streamflow_percentage,
    relative_CMAES_DREAM_diff = (CMAES_DREAM_diff / median_CO2_impact_on_streamflow_percentage) * 100
  ) |>
  arrange(CMAES_DREAM_diff)

compare_CMAES_and_DREAM |>
  ggplot(aes(y = relative_CMAES_DREAM_diff)) +
  geom_boxplot(
    staplewidth = 0.5
  ) +
  labs(y = "Relative percentage difference between best CMAES percentage and median DREAM percentage") +
  theme_bw()

# They are similar



# Looking at average impacts ---------------------------------------------------
## Average percentage impact by state ==========================================
## (counterfactual)
plot_ready_decade_differences |>
  drop_na() |>
  mutate(state = if_else(state == "ACT", "NSW", state)) |>
  summarise(
    mean = mean(CO2_impact_on_streamflow_percent),
    sd = sd(CO2_impact_on_streamflow_percent),
    median = median(CO2_impact_on_streamflow_percent),
    .by = c(decade, state)
  ) |>
  arrange(state, decade)



## Average uncertainty by state using DREAM ====================================
DREAM_CO2_impact_uncertainty_on_streamflow |>
  left_join(
    gauge_information,
    by = join_by(gauge)
  ) |>
  mutate(state = if_else(state == "ACT", "NSW", state)) |>
  summarise(
    mean = mean(IQR_CO2_impact_on_streamflow_percentage),
    median = median(IQR_CO2_impact_on_streamflow_percentage),
    .by = c(state, decade)
  )












# Timeseries of streamflow data for extended figures ---------------------------
## Simplify tibble for plotting ================================================

timeseries_plotting_data <- high_evidence_master_streamflow_table |>
  select(
    gauge, year, realspace_observed_streamflow,
    realspace_streamflow_no_CO2_model, realspace_streamflow_CO2_model_on,
    realspace_streamflow_CO2_model_off
  ) |> 
  pivot_longer(
    cols = starts_with("realspace"),
    names_to = "type",
    values_to = "streamflow"
  ) |> 
  mutate(
    type = case_when(
      type == "realspace_streamflow_CO2_model_on" ~ "CO2 Model",
      type == "realspace_streamflow_CO2_model_off" ~ "Counterfactual",
      type == "realspace_observed_streamflow" ~ "Observed",
      type == "realspace_streamflow_no_CO2_model" ~ "non-CO2 Model",
      .default = NA
    )
  ) |>
  # set order
  mutate(
    type = factor(type, levels = c("Observed", "CO2 Model", "Counterfactual", "non-CO2 Model"))
  )






## Divide plots for extended figures (6 pages worth) ===========================
## divide all_timeseries_data into lots of 7 using split


chunk <- 14
n <- length(high_evidence_ratio_gauges)
split_group <- rep(rep(1:ceiling(n / chunk), each = chunk))[1:n]
split_tibble <- tibble(
  "gauge" = high_evidence_ratio_gauges,
  "split" = split_group
)


# left_join and split
timeseries_plotting_data <- timeseries_plotting_data |>
  left_join(
    split_tibble,
    by = join_by(gauge)
  )


chunked_timeseries_data <- timeseries_plotting_data |> # converting table to list by groups https://stackoverflow.com/questions/7060272/split-up-a-dataframe-by-number-of-rows
  group_by(split) |>
  group_map(~.x)



## Calculate NSE to be added to timeseries figures  ============================
# all_timeseries_data has been filtered to only include moderately_strong_or_greater_catchments
# thus, CO2 model is the best model
# the easiest way to add NSE is in the caption
# the next easiest way is to add NSE to make_facet_labels
# - could make a wrapper around make facet labels - that appends NSE to letter

# in utility file
# nash_sutcliffe_efficiency <- function(observed, modelled) {
#  1 - (sum((observed - modelled)^2) / sum((observed - mean(observed))^2))
# }


NSE_moderately_strong_or_greater_gauges <- timeseries_plotting_data |>
  select(gauge, year, type, streamflow) |>
  filter(type != "Counterfactual") |>
  pivot_wider(
    names_from = type,
    values_from = streamflow
  ) |>
  summarise(
    NSE_CO2_model = nash_sutcliffe_efficiency(observed = Observed, modelled = `CO2 Model`),
    NSE_non_CO2_model = nash_sutcliffe_efficiency(observed = Observed, modelled = `non-CO2 Model`),
    .by = gauge
  ) |>
  # signif + label
  mutate(
    lab_NSE_CO2_model = paste0("NSE = ", format(round(NSE_CO2_model, 2), nsmall = 2)),
    lab_NSE_non_CO2_model = paste0("NSE = ", format(round(NSE_non_CO2_model, 2), nsmall = 2))
  )


# Split up into a4 size and repeat
# Here
# It needs to be A, B, C, D labels rather than facet_labels
# A, B, C and D must be in the same relative location for all plots
# Save A, B, C and D catchments. Print out caption
# Colours are not good - need colour blind safe selection

make_facet_labels <- function(data, facet_column, x_axis_column, y_axis_column, label_type = LETTERS, hjust = 0, vjust = 0) {
  # The embrace operator does not work correctly in summarise i.e., max({{ y_axis_column }})
  # Link: https://forum.posit.co/t/embrace-operator-for-tidy-selection-vs-data-masking/173084
  # Possible cause: {{ y_axis_column }} isn't unquoting when it's doing the mutate
  # Work around using rlang::ensym
  col <- rlang::ensym(y_axis_column)
  
  data |>
    summarise(
      ylab = max(!!col),
      .by = {{ facet_column }}
    ) |>
    # Add xlab - constant x-axis
    add_column(
      xlab = data |> pull(x_axis_column) |> min(),
      .before = 2
    ) |>
    # add row numbers to tibble
    mutate(
      row_number = row_number(),
      .before = 1
    ) |>
    # add label type based on row number
    mutate(
      label_name = label_type[row_number]
    ) |>
    # apply hjust and vjust
    mutate(
      xlab = xlab + (xlab * hjust),
      ylab = ylab + (ylab * vjust)
    )
}


add_NSE_to_facet_labels <- function(data, NSE_data, facet_column, x_axis_column, y_axis_column, label_type = LETTERS, hjust = 0, vjust = 0) {
  # make facet labels
  facet_labels <- make_facet_labels(
    data = data,
    facet_column = {{ facet_column }},
    x_axis_column = {{ x_axis_column }},
    y_axis_column = {{ y_axis_column }},
    label_type = label_type,
    hjust = hjust,
    vjust = vjust
  )
  
  
  NSE_data |>
    select(gauge, lab_NSE_CO2_model) |>
    # this filters out gauges
    right_join(
      facet_labels,
      by = join_by(gauge)
    )
}


plot_and_save_timeseries_data <- function(plotting_data, label_data, NSE_label_data, identifier) {
  plot <- plotting_data |>
    ggplot(aes(x = year, y = streamflow, colour = type, shape = type)) +
    geom_line(alpha = 0.8, linewidth = 0.5) +
    geom_point(size = 1) +
    geom_text(
      mapping = aes(x = xlab, y = ylab, label = label_name),
      data = label_data,
      inherit.aes = FALSE,
      fontface = "bold",
      size = 8,
      size.unit = "pt"
    ) +
    geom_text(
      mapping = aes(x = xlab, y = ylab, label = lab_NSE_CO2_model),
      data = NSE_label_data,
      inherit.aes = FALSE,
      # fontface = "bold",
      size = 8,
      size.unit = "pt",
      hjust = -4.1
    ) +
    scale_colour_brewer(palette = "Set1") +
    labs(
      x = "Time (Year)",
      y = "Streamflow (mm)",
      colour = NULL,
      shape = NULL
    ) +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      legend.text = element_text(size = 8),
      strip.background = element_blank(), # remove facet_strip gauge numbers
      strip.text = element_blank(), # remove facet_strip gauge numbers
      panel.grid.minor = element_blank()
    ) +
    facet_wrap(~gauge, ncol = 2, nrow = 7, scales = "free_y")
  
  ggsave(
    filename = paste0("streamflow_timeseries_extended_data_plot_", identifier, ".pdf"),
    path = "Figures/Extended_Data",
    device = "pdf",
    plot = plot,
    width = 180,
    height = 170,
    units = "mm"
  )
}


# create txt file with figure captions for overleaf
create_caption <- function(label_data, identifier) {
  abc <- label_data |> pull(label_name)
  gauge <- label_data |> pull(gauge)
  gauge_abc <- paste0(gauge, " (", abc, ")")
  # concatenate everything but last value
  start_gauge_abc <- paste0(gauge_abc[1:(length(gauge_abc) - 2)], ", ", collapse = "")
  end_gauge_abc <- paste0(gauge_abc[(length(gauge_abc) - 1)], " and ", gauge_abc[length(gauge_abc)], ".")
  gauge_text <- paste(c(start_gauge_abc, end_gauge_abc), collapse = "")
  
  cat("\\begin{figure}")
  cat("\n")
  cat("\t\\centering")
  cat("\n")
  cat(paste0("\t\\includegraphics[width=\\textwidth]{Figures/streamflow_timeseries_supp_plot_", identifier, ".pdf}"))
  cat("\t\n")
  # The line below must change
  cat(paste0("\t\\caption{\\textbf{The impact of CO$_2$ on the streamflow timeseries for gauges ", gauge_text, "} The streamflow timeseries compares observed streamflow (Observed), modelled streamflow using a model that includes CO$_2$ (CO$_2$ model), modelled streamflow using a model that includes CO$_2$ with CO$_2$ turned off (Counterfactual) and modelled streamflow using a model that does not include CO$_2$ (non-CO$_2$ Model).}"))
  cat("\n")
  # The line below must change
  cat(paste0("\t\\label{fig:supp_streamflow_timeseries_", identifier, "}"))
  cat("\n")
  cat("\\end{figure}")
  cat("\n")
  cat("\n")
}


save_plot_and_caption_timeseries_data <- function(data_chunk, NSE_data, identifier) {
  label_data <- make_facet_labels(
    data = data_chunk,
    facet_column = "gauge",
    x_axis_column = "year",
    y_axis_column = "streamflow",
    label_type = letters,
    hjust = 0.0005,
    vjust = -0.05
  )
  
  NSE_label_data <- add_NSE_to_facet_labels(
    data = data_chunk,
    NSE_data = NSE_data,
    facet_column = "gauge",
    x_axis_column = "year",
    y_axis_column = "streamflow",
    label_type = letters,
    hjust = 0.0005,
    vjust = -0.05
  )
  
  plot_and_save_timeseries_data(
    plotting_data = data_chunk,
    label_data = label_data,
    NSE_label_data = NSE_label_data,
    identifier = identifier
  )
  
  create_caption(
    label_data = label_data,
    identifier = identifier
  )
}





## Save to file ================================================================
# Supplementary time series figures and captions:
sink(file = "Figures/Extended_Data/streamflow_time_captions_extended_data.txt") # filename must change
iwalk(
  .x = chunked_timeseries_data,
  .f = save_plot_and_caption_timeseries_data,
  NSE_data = NSE_moderately_strong_or_greater_gauges
)
sink()








# Making percentage change map plot --------------------------------------------
aus_map <- generate_aus_map_sf()

## Limits and breaks ===========================================================
make_limits <- function(timeseries) {
  # round up to next whole number
  limits <- timeseries |> range()
  sign_limits <- sign(limits)

  sign_limits * ceiling(abs(limits))
}

CO2_impact_on_streamflow_percent_limits <- plot_ready_decade_differences |>
  pull(CO2_impact_on_streamflow_percent) |>
  make_limits() |>
  as.double()

scale_size_limits <- plot_ready_decade_differences |>
  pull(IQR_CO2_impact_on_streamflow_percentage) |>
  range(na.rm = TRUE) |> # can round up if I want to
  round(digits = 0)

percentage_IQR_breaks <- c(0, 2.5, 5, 10, 15, 50, 100) # custom breaks

hard_coded_breaks_CO2_impact_of_streamflow <- c(-75, -50, -25, -10, -1, 0, 1, 10, 25, 50, 75)

dot_transparency <- 0.8


## Custom colour palette =======================================================
big_palette <- function(x) {
  c(
    "#67001f",
    "#b2182b",
    "#d6604d",
    "#f4a582",
    "#fddbc7",
    "white",
    "white",
    "#d1e5f0",
    "#92c5de",
    "#4393c3",
    "#2166ac",
    "#053061"
  )
}



## Plotting function ===========================================================

### No uncertainty map ###
make_CO2_streamflow_percentage_change_map <- function(data, title) {
  ## Generate Insets ===========================================================
  QLD_data <- data |>
    filter(state == "QLD")

  NSW_data <- data |>
    filter(state == "NSW")

  VIC_data <- data |>
    filter(state == "VIC")

  WA_data <- data |>
    filter(state == "WA")

  TAS_data <- data |>
    filter(state == "TAS")


  ### Generate inset plots #######################################################
  inset_dot_size <- 1.8

  inset_plot_QLD <- aus_map |>
    filter(state == "QLD") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = QLD_data,
      aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent),
      show.legend = FALSE,
      size = inset_dot_size,
      alpha = dot_transparency,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    guides(size = guide_bins(show.limits = TRUE)) +
    theme_void()


  inset_plot_NSW <- aus_map |>
    filter(state == "NSW") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = NSW_data,
      aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent),
      show.legend = FALSE,
      size = inset_dot_size,
      alpha = dot_transparency,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    guides(size = guide_bins(show.limits = TRUE)) +
    theme_void()


  inset_plot_VIC <- aus_map |>
    filter(state == "VIC") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = VIC_data,
      aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent),
      show.legend = FALSE,
      size = inset_dot_size,
      alpha = dot_transparency,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    guides(size = guide_bins(show.limits = TRUE)) +
    theme_void()


  inset_plot_WA <- aus_map |>
    filter(state == "WA") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = WA_data,
      aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent),
      show.legend = FALSE,
      size = inset_dot_size,
      alpha = dot_transparency,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    guides(size = guide_bins(show.limits = TRUE)) +
    theme_void()


  inset_plot_TAS <- aus_map |>
    filter(state == "TAS") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = TAS_data,
      aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent),
      show.legend = FALSE,
      size = inset_dot_size,
      alpha = dot_transparency,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    guides(size = guide_bins(show.limits = TRUE)) +
    theme_void()


  ## Put it together =============================================================
  single_map_aus <- aus_map |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = data,
      mapping = aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent),
      alpha = dot_transparency,
      size = inset_dot_size,
      colour = "black",
      shape = 21,
      inherit.aes = FALSE,
      stroke = 0.1
    ) +
    theme_bw() +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    # expand map
    coord_sf(xlim = c(95, 176), ylim = c(-60, 0)) +
    # magnify WA
    geom_magnify(
      from = c(114, 118, -35.5, -30),
      to = c(93, 112, -36, -10),
      shadow = FALSE,
      expand = 0,
      plot = inset_plot_WA,
      proj = "single"
    ) +
    # magnify VIC
    geom_magnify(
      # aes(from = state == "VIC"), # use aes rather than manually selecting area
      from = c(141, 149.5, -39, -34),
      to = c(95, 136, -38, -60),
      shadow = FALSE,
      plot = inset_plot_VIC,
      proj = "single"
    ) +
    # magnify QLD
    geom_magnify(
      from = c(145, 155, -29.2, -15),
      to = c(157, 178, -29.5, 1.5),
      shadow = FALSE,
      expand = 0,
      plot = inset_plot_QLD,
      proj = "single"
    ) +
    # magnify NSW
    geom_magnify(
      from = c(146.5, 154, -38, -28.1),
      to = c(157, 178, -61, -30.5),
      shadow = FALSE,
      expand = 0,
      plot = inset_plot_NSW,
      proj = "single"
    ) +
    # magnify TAS
    geom_magnify(
      from = c(144, 149, -40, -44),
      to = c(140, 155, -45, -61),
      shadow = FALSE,
      expand = 0,
      plot = inset_plot_TAS,
      proj = "single"
    ) +
    labs(
      x = NULL, # "Latitude",
      y = NULL, # "Longitude",
      fill = bquote("Average Impact of" ~ CO[2] ~ "on Streamflow (%)"),
      size = "Percentage Impact Uncertainty (IQR)",
      title = {{ title }}
    ) +
    theme(
      legend.key = element_rect(fill = "white"),
      legend.title = element_text(hjust = 0.5),
      # legend.background = element_rect(colour = "black"), #this cuts off the negative sign
      axis.text = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.351, 0.9),
      legend.box = "horizontal", # side-by-side legends
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(margin = margin(l = 25, r = 0, t = 30, b = -30), size = 12, face = "bold") # push title into plot
    ) +
    guides(
      fill = guide_coloursteps(
        barwidth = unit(10, "cm"),
        show.limits = TRUE,
        even.steps = TRUE,
        title.position = "top",
        direction = "horizontal"
      ),
      size = guide_bins(
        override.aes = aes(stroke = 0.5),
        show.limits = TRUE,
        direction = "horizontal",
        title.position = "top", # warnings says its ignore these parameter - The warnings are wrong
        barwidth = unit(1, "cm")
      )
    )

  return(single_map_aus)
}



# This is required in catchment boundary plots
write_csv(
  x = plot_ready_decade_differences,
  file = "Modelling/decade_streamflow_CO2_differences.csv"
)

## Map of percentage differences between decades ===============================
percentage_difference_CO2_model_non_CO2_model_1990s <- plot_ready_decade_differences |>
  filter(decade == 1)

percentage_difference_CO2_model_non_CO2_model_2010s <- plot_ready_decade_differences |>
  filter(decade == 2)

patchwork_CO2_model_and_non_CO2_model_percentage_differences <- (make_CO2_streamflow_percentage_change_map(percentage_difference_CO2_model_non_CO2_model_1990s, "a") | make_CO2_streamflow_percentage_change_map(percentage_difference_CO2_model_non_CO2_model_2010s, "b")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")









# Make plot for main figure ----------------------------------------------------
# Must include:
# 1. patchwork_CO2_model_and_non_CO2_model_percentage_differences (map)
# 2. a handful of streamflow timeseries 
# 3. Uncertainty from DREAM

## Select a handful of catchments to present in the main figure ================
short_list_catchments <- c("401210", "606195", "701002", "407246") # select 4


## Alter make_CO2_streamflow_percentage_change_map to include uncertainty and labels =========
map_label_table <- gauge_information |>
  filter(gauge %in% short_list_catchments) |>
  mutate(label_name = letters[1:4])


# THIS FUNCTION USES THE SAME LIMITS AND BREAKS AS THE PREVIOUS
# GLOBAL VARIABLES
make_CO2_streamflow_percentage_change_map_uncertainty <- function(data, title, A_or_B) {
  font_size_pt <- 9L # default size is GeomLabel$default_aes$size = 3.88
  
  ## Generate Insets ===========================================================
  QLD_data <- data |>
    filter(state == "QLD")
  
  NSW_data <- data |>
    filter(state == "NSW")
  
  VIC_data <- data |>
    filter(state == "VIC")
  
  WA_data <- data |>
    filter(state == "WA")
  
  TAS_data <- data |>
    filter(state == "TAS")
  
  
  ### Generate inset plots #######################################################
  inset_dot_size <- 1.8
  
  inset_plot_QLD <- aus_map |>
    filter(state == "QLD") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = QLD_data,
      aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent, size = IQR_CO2_impact_on_streamflow_percentage),
      show.legend = FALSE,
      alpha = dot_transparency,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    guides(size = guide_bins(show.limits = TRUE)) +
    theme_void()
  
  
  inset_plot_NSW <- aus_map |>
    filter(state == "NSW") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = NSW_data,
      aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent, size = IQR_CO2_impact_on_streamflow_percentage),
      show.legend = FALSE,
      alpha = dot_transparency,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    guides(size = guide_bins(show.limits = TRUE)) +
    theme_void()
  
  
  # Labels for VIC
  VIC_map_label_table <- map_label_table |>
    filter(state == "VIC")
  
  inset_plot_VIC <- aus_map |>
    filter(state == "VIC") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = VIC_data,
      aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent, size = IQR_CO2_impact_on_streamflow_percentage),
      show.legend = FALSE,
      alpha = dot_transparency,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    geom_text(
      data = VIC_map_label_table,
      aes(x = lon, y = lat, label = label_name),
      nudge_x = -0.35,
      size = font_size_pt,
      size.unit = "pt"
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    guides(size = guide_bins(show.limits = TRUE)) +
    theme_void()
  
  
  # Labels for WA
  WA_map_label_table <- map_label_table |>
    filter(gauge == "606195") # 701002 not in subset plot so exclude
  
  inset_plot_WA <- aus_map |>
    filter(state == "WA") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = WA_data,
      aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent, size = IQR_CO2_impact_on_streamflow_percentage),
      show.legend = FALSE,
      alpha = dot_transparency,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    geom_text(
      data = WA_map_label_table,
      aes(x = lon, y = lat, label = label_name),
      nudge_y = -0.25,
      size = font_size_pt,
      size.unit = "pt"
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    guides(size = guide_bins(show.limits = TRUE)) +
    theme_void()
  
  
  inset_plot_TAS <- aus_map |>
    filter(state == "TAS") |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = TAS_data,
      aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent, size = IQR_CO2_impact_on_streamflow_percentage),
      show.legend = FALSE,
      alpha = dot_transparency,
      colour = "black",
      stroke = 0.1,
      shape = 21
    ) +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    guides(size = guide_bins(show.limits = TRUE)) +
    theme_void()
  
  
  ## Put it together =============================================================
  
  # Big map label table
  big_map_label_table <- map_label_table |>
    filter(gauge == "701002")
  
  figure_label <- tribble(
    ~lon, ~lat, ~label_name,
    95, 0, A_or_B
  )
  
  decade_label <- tribble(
    ~lon, ~lat, ~label_name,
    105, 0, title
  )
  
  single_map_aus <- aus_map |>
    ggplot() +
    geom_sf() +
    geom_point(
      data = data,
      mapping = aes(x = lon, y = lat, fill = CO2_impact_on_streamflow_percent, size = IQR_CO2_impact_on_streamflow_percentage),
      alpha = dot_transparency,
      colour = "black",
      shape = 21,
      inherit.aes = FALSE,
      stroke = 0.1
    ) +
    geom_text(
      data = big_map_label_table,
      aes(x = lon, y = lat, label = label_name),
      nudge_y = 2,
      size = 10,
      size.unit = "pt"
    ) +
    geom_text(
      data = figure_label,
      aes(x = lon, y = lat, label = label_name),
      fontface = "bold",
      size = 10,
      size.unit = "pt"
    ) +
    geom_text(
      data = decade_label,
      aes(x = lon, y = lat, label = label_name),
      size = 10,
      size.unit = "pt"
    ) +
    theme_bw() +
    binned_scale( # https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "fill",
      palette = big_palette,
      breaks = hard_coded_breaks_CO2_impact_of_streamflow,
      limits = CO2_impact_on_streamflow_percent_limits,
      show.limits = TRUE,
      guide = "colorsteps"
    ) +
    scale_size_binned(limits = scale_size_limits, breaks = percentage_IQR_breaks) + # range = c(0, 2) dictates the size of the dots (important)
    # expand map
    coord_sf(xlim = c(95, 176), ylim = c(-60, 0)) +
    # magnify WA
    geom_magnify(
      from = c(114, 118, -35.5, -30),
      to = c(93, 112, -36, -10),
      shadow = FALSE,
      expand = 0,
      plot = inset_plot_WA,
      proj = "single"
    ) +
    # magnify VIC
    geom_magnify(
      # aes(from = state == "VIC"), # use aes rather than manually selecting area
      from = c(141, 149.5, -39, -34),
      to = c(95, 136, -38, -60),
      shadow = FALSE,
      plot = inset_plot_VIC,
      proj = "single"
    ) +
    # magnify QLD
    geom_magnify(
      from = c(145, 155, -29.2, -15),
      to = c(157, 178, -29.5, 1.5),
      shadow = FALSE,
      expand = 0,
      plot = inset_plot_QLD,
      proj = "single"
    ) +
    # magnify NSW
    geom_magnify(
      from = c(146.5, 154, -38, -28.1),
      to = c(157, 178, -61, -30.5),
      shadow = FALSE,
      expand = 0,
      plot = inset_plot_NSW,
      proj = "single"
    ) +
    # magnify TAS
    geom_magnify(
      from = c(144, 149, -40, -44),
      to = c(140, 155, -45, -61),
      shadow = FALSE,
      expand = 0,
      plot = inset_plot_TAS,
      proj = "single"
    ) +
    labs(
      x = NULL, # "Latitude",
      y = NULL, # "Longitude",
      fill = bquote("Average Impact of" ~ CO[2] ~ "on Streamflow (%)"),
      size = "Percentage Impact Uncertainty (IQR)"
    ) +
    theme(
      legend.key = element_rect(fill = "white"),
      legend.title = element_text(hjust = 0.5),
      text = element_text(family = "sans"),
      # legend.background = element_rect(colour = "black"), #this cuts off the negative sign
      axis.text = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.351, 0.9),
      legend.box = "horizontal", # side-by-side legends
      panel.grid = element_blank(),
      axis.ticks = element_blank()
    ) +
    guides(
      fill = guide_coloursteps(
        barwidth = unit(10, "cm"),
        show.limits = TRUE,
        even.steps = TRUE,
        title.position = "top",
        direction = "horizontal"
      ),
      size = guide_bins(
        override.aes = aes(stroke = 0.5),
        show.limits = TRUE,
        direction = "horizontal",
        title.position = "top", # warnings says its ignore these parameter - The warnings are wrong
        barwidth = unit(1, "cm")
      )
    )
  
  return(single_map_aus)
}



## Make map components =========================================================
percent_change_1990 <- make_CO2_streamflow_percentage_change_map_uncertainty(
  percentage_difference_CO2_model_non_CO2_model_1990s,
  "1990-1999",
  "A"
)

percent_change_2012 <- make_CO2_streamflow_percentage_change_map_uncertainty(
  percentage_difference_CO2_model_non_CO2_model_2010s,
  "2012-2021",
  "B"
)



## Make the streamflow-time plots ==============================================
facet_annotation <- timeseries_plotting_data |>
  filter(gauge %in% short_list_catchments) |>
  summarise(
    streamflow = max(streamflow) - (max(streamflow) * 0.1),
    .by = gauge
  ) |>
  add_column(
    year = 1960,
    label_name = paste0(LETTERS[3:6], "(", letters[1:4], ")")
  )


### Shading decades ############################################################
shade_decade_1 <- timeseries_plotting_data |>
  filter(gauge %in% short_list_catchments) |>
  group_by(gauge) |>
  mutate(upper = max(streamflow) * 1.2) |>
  filter(year %in% decade_1) |> 
  select(gauge, year, upper) |> 
  distinct()

shade_decade_2 <- timeseries_plotting_data |>
  filter(gauge %in% short_list_catchments) |>
  group_by(gauge) |>
  mutate(upper = max(streamflow) * 1.2) |>
  filter(year %in% decade_2)
# 606195 is missing 2021
# 407246 is missing 2020 and 2021

# easiest solution is extract a year - replace year and streamflow and rbind
missing_years_606195 <- shade_decade_2 |>
  filter(gauge == "606195") |>
  filter(year == 2020) |>
  mutate(
    year = 2021,
    precipitation = NA,
    streamflow = NA
  )

missing_years_407246 <- shade_decade_2 |>
  filter(gauge == "407246") |>
  filter(year %in% c(2018, 2019)) |>
  mutate(
    year = if_else(year == 2018, 2020, 2021),
    precipitation = NA,
    streamflow = NA
  )

shade_decade_2 <- rbind(shade_decade_2, missing_years_606195, missing_years_407246) |> 
  select(gauge, year, upper) |> 
  distinct()


### Plotting ###################################################################
main_figure_timeseries_component <- timeseries_plotting_data |>
  filter(gauge %in% short_list_catchments) |>
  ggplot(aes(x = year, y = streamflow, colour = type, shape = type)) +
  geom_area(
    aes(x = year, y = upper),
    inherit.aes = FALSE,
    data = shade_decade_1,
    alpha = 0.08
  ) +
  geom_area(
    aes(x = year, y = upper),
    inherit.aes = FALSE,
    data = shade_decade_2,
    alpha = 0.08
  ) +
  geom_line(alpha = 0.8) +
  geom_point(size = 1) +
  geom_label(
    aes(x = year, y = streamflow, label = label_name),
    data = facet_annotation,
    inherit.aes = FALSE,
    fill = NA,
    label.size = NA,
    fontface = "bold",
    size = 10,
    size.unit = "pt"
  ) +
  # scale_colour_brewer(palette = "Set1") +
  scale_colour_brewer(
    labels = c(
      "Observed",
      bquote(~ CO[2] ~ "Model"),
      "Counterfactual",
      bquote("non-" * CO[2] ~ "Model")
    ),
    palette = "Set1"
    # values = c("red", "green", "blue", "orange")
  ) +
  scale_shape_manual(
    labels = c(
      "Observed",
      bquote(~ CO[2] ~ "Model"),
      "Counterfactual",
      bquote("non-" * CO[2] ~ "Model")
    ),
    values = c(16, 17, 15, 3)
  ) +
  labs(
    x = "Year",
    y = "Streamflow (mm)",
    colour = "Streamflow Timeseries",
    shape = "Streamflow Timeseries"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(), # remove strip labels from faceting
    strip.text = element_blank(), # remove strip labels from faceting
    panel.grid.minor = element_blank(), # remove minor gridlines
    axis.title = element_text(size = 10)
  ) +
  scale_y_continuous(expand = c(0, 0)) + # remove y-axis padding
  facet_wrap(~gauge, scales = "free_y", ncol = 2)



## Use patchwork to join everything together ===================================
# Alternative methods attempted:
## - plot_spacer and | + / operations - too much faff
## - individually plotting timeseries - gaps

# Defining layout it the best method as I have direct control

# use the area() constructor
# top, left, bottom, right bounds (t < b and l < r)
layout <- c(
  area(t = 1, l = 1, b = 3, r = 3), # 1990s percentage change
  area(t = 1, l = 4, b = 3, r = 6), # 2010s percentage change
  area(t = 4, l = 1, b = 4, r = 6) # timeseries
)

# plot(layout) # check the patches are working

streamflow_percentage_difference_with_timeseries <- (percent_change_1990 + percent_change_2012 + main_figure_timeseries_component) +
  plot_layout(design = layout, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  ) &
  guides(colour = guide_legend(title.hjust = 0.5, title.position = "top", ncol = 2))



## SAVE ========================================================================
ggsave(
  filename = "./Figures/Main/streamflow_percentage_difference_with_timeseries.pdf",
  plot = streamflow_percentage_difference_with_timeseries,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)












# OTHER ANALYSIS ---------------------------------------------------------------
## Comparing percentage changes with aridity ===================================
# file created in ET_analysis.R
aridity_information <- read_csv(
  file = "Data/Tidy/aridity_information.csv",
  show_col_types = FALSE
)

aridity_labels <- aridity_information |> pull(dryness_zone) |> unique()

aridity_and_streamflow_changes <- percentage_difference_CO2_model_non_CO2_model_2010s |>
  select(
    gauge,
    CO2_impact_on_streamflow_percent,
    state
  ) |>
  left_join(
    aridity_information,
    by = join_by(gauge, state)
  ) |>
  mutate(dryness_zone = factor(dryness_zone, levels = aridity_labels))


aridity_and_streamflow_changes |>
  select(
    gauge,
    CO2_impact_on_streamflow_percent,
    state
  ) |>
  left_join(
    aridity_information,
    by = join_by(gauge)
  ) |>
  ggplot(aes(x = dryness_zone, y = CO2_impact_on_streamflow_percent)) +
  geom_point() +
  labs(
    x = NULL,
    y = "Percentage Streamflow Change (2012-2021)"
  ) +
  theme_bw()


split_index <- aridity_and_streamflow_changes |>
  pull(dryness_zone)


split_aridity_streamflow_changes <- split(aridity_and_streamflow_changes, split_index)

streamflow_aridity_maps <- map2(
  .x = split_aridity_streamflow_changes,
  .y = aridity_labels,
  .f = make_CO2_streamflow_percentage_change_map
)


streamflow_aridity_percentage_change_map_2010s <- streamflow_aridity_maps[[1]] +
  streamflow_aridity_maps[[2]] +
  streamflow_aridity_maps[[3]] +
  streamflow_aridity_maps[[4]] +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(
  filename = "streamflow_aridity_percentage_change_map_2010s.pdf",
  path = "Figures/Other",
  plot = streamflow_aridity_percentage_change_map_2010s,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)



### Are arid catchments more likely to have an inter-year autocorrelation term? ###
aridity_inter_year_auto_comparison <- best_CO2_non_CO2_per_gauge |>
  filter(str_detect(streamflow_model, "CO2")) |> 
  select(gauge, streamflow_model) |>
  distinct() |>
  mutate(
    contains_auto = str_detect(streamflow_model, "auto")
  ) |>
  left_join(
    aridity_information,
    by = join_by(gauge)
  ) |>
  filter(evidence_ratio > 100)

aridity_inter_year_auto_comparison |>
  count(contains_auto, dryness_zone) |>
  arrange(dryness_zone) |>
  pivot_wider(
    names_from = contains_auto,
    values_from = n
  ) |>
  mutate(
    percentage = `TRUE` / (`TRUE` + `FALSE`)
  )
























