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

evidence_ratio <- read_csv(
  "./Modelling/Results/CMAES/evidence_ratio_results.csv",
  show_col_types = FALSE
)


# Calculate evidence ratio for filtering ---------------------------------------
high_evi_gauges <- evidence_ratio |> 
  filter(evidence_ratio > 100) |> 
  pull(gauge)




# Make rainfall-runoff data ----------------------------------------------------
## The following function relies on catchment_data
## Must be rearraged to be used
rearrange_catchment_data_blueprint <- function(observed_data, gauge_ID, start_stop_indexes) {
  catchment_data_blueprint(
    gauge_ID = gauge_ID,
    observed_data = observed_data,
    start_stop_indexes = start_stop_indexes
  )
}


replace_col_in_data <- function(replacement, col, data_tibble) {
  data_tibble |> 
    mutate({{ col }} := !!enquo(replacement))
}



## Input: replace_precipitation = either NULL or single value
## Return: output depends on model_year and replace_precipitation if 
##         replace precpitation is NULL then output is equal to length of 
##         obs precip x model_year.
##         If replace_precip is specific then output is length(model_years)
##         Does return log-sinh streamflow
make_rainfall_runoff_results <- function(gauge, model_years, replace_precipitation, data, model_results) {
  
  
  # 1. Transform observed streaflow into log-sinh space
  b <- model_results |>
    filter(contains_CO2) |>
    filter(gauge == {{ gauge }}) |>
    filter(parameter == "b") |>
    pull(parameter_value)

  modified_data <- data |>
    filter(gauge == {{ gauge }}) |>
    mutate(
      q_log_sinh = log_sinh_transform(b = b, y = q_mm, offset = 0)
    )

  # 2. Replace data (CO2 and replace_precipitation)
  if (!is.null(replace_precipitation)) {
    stopifnot(
      is.numeric(replace_precipitation),
      length(replace_precipitation) == 1
    )

    modified_data <- modified_data |>
      mutate(
        p_mm = replace_precipitation
      )
  }

  # check models_years are not outside of data
  # This includes na values
  year_range <- data |> pull(year) |> range()
  if(any(model_years < min(year_range))) {stop("model_year less than observed year")}
  if(any(model_years > max(year_range))) {stop("model_year greater than observed year")}
  
  CO2_for_selected_years <- modified_data |>
    filter(year %in% model_years) |>
    pull(CO2)


  different_CO2_modified_data <- map(
    .x = CO2_for_selected_years,
    .f = replace_col_in_data,
    col = CO2,
    data_tibble = modified_data
  )

  # 3. Convert replace data into catchment_blueprint objects
  different_CO2_catchment_data <- map(
    .x = different_CO2_modified_data,
    .f = rearrange_catchment_data_blueprint,
    gauge = gauge,
    start_stop_indexes = start_stop_indexes
  )


  # 4. Set parameters to zero to get simplest CO2 model (i.e., auto, drought)
  ## Identify drought
  non_drought_intercept <- model_results |>
    filter(contains_CO2) |>
    filter(gauge == {{ gauge }}) |>
    filter(parameter == "a0_n") |>
    pull(parameter_value)

  # empty vectors do not play nice
  if (is_empty(non_drought_intercept)) {
    non_drought_intercept <- NA
  }

  parameter_set <- model_results |>
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


  # 5. Get the streamflow model
  streamflow_model <- model_results |>
    filter(contains_CO2) |>
    select(gauge, streamflow_model) |>
    distinct() |>
    filter(gauge == {{ gauge }}) |>
    pull(streamflow_model) |>
    match.fun()


  # 6. Put the catchment_data object, parameter_set into streamflow model
  different_CO2_streamflow <- map(
    .x = different_CO2_catchment_data,
    .f = streamflow_model,
    parameter = parameter_set
  ) |>
    `names<-`(model_years) |> # this may be changed to another name to avoid confusion
    list_rbind(names_to = "source")

  # 7. Convert log-sinh into realspace streamflow
  back_to_realspace <- different_CO2_streamflow |>
    mutate(
      realspace_streamflow_results = inverse_log_sinh_transform(
        b = b,
        z = streamflow_results,
        offset = 0
      )
    ) |> 
    select(!c(is_drought_year, seasonal_ratio)) 

  # 8. Final return 
  return(back_to_realspace)
}



# Use make_rainfall_runoff_results for percentage difference calc --------------
# wrapper around make_rainfall_runoff_results
# it could be a function factory angle
get_percentage_difference <- function(gauge, model_years, replace_precipitation, data, model_results) {
  
  # 1. replace_precipiation must be a single digit
  stopifnot(!is.null(replace_precipitation))
  
  # 2. Repeat function (return n x 2 tibble where n is model years)
  percentage_changes <- make_rainfall_runoff_results(
    gauge = {{ gauge }},
    model_years = model_years,
    replace_precipitation = replace_precipitation,
    data = data,
    model_results = model_results
  ) |> 
    select(source, realspace_streamflow_results) |> 
    distinct()
  
  # 3. Calculate between mix and max source (earliest year and latest year)
  late_streamflow <- percentage_changes |> slice_max(source, n = 1) |> pull(realspace_streamflow_results)
  early_streamflow <- percentage_changes |> slice_min(source, n = 1) |> pull(realspace_streamflow_results)
  
  (late_streamflow - early_streamflow) / early_streamflow

}






# Use make_rainfall_runoff_results for plotting --------------------------------
# wrapper around make_rainfall_runoff_results
rainfall_runoff_plots <- function(gauge, model_years, plot_label, replace_precipitation, data, model_results, plot_arrow, overwrite_legend_text) {
  
  stopifnot(is.null(replace_precipitation))
  
  # 1. Make plotting data
  plotting_data <- make_rainfall_runoff_results(
    gauge = {{ gauge }},
    model_years = model_years,
    replace_precipitation = replace_precipitation,
    data = data,
    model_results = model_results
  ) |> 
    mutate(
      source = paste0("Estimated Relationship ", source)
    )
  
  # If !is.null(overwrite_legend_text)
  # overwrite_legend_text same length as model years
  # use order to assign them
  if (!is.null(overwrite_legend_text)) {
    plotting_data <- plotting_data |>
      mutate(
        source = recode_values(
          source,
          from = unique(source),
          to = overwrite_legend_text
        )
      )
  }
    
  
  
  # 2. Plot
  plot <- plotting_data |>
    drop_na() |> 
    ggplot(aes(x = precipitation, y = observed_streamflow, fill = year)) +
    geom_line(
      aes(x = precipitation, y = streamflow_results, colour = source)
    ) +
    geom_point(
      colour = "black",
      shape = 21,
      size = 3
    ) +
    labs(
      x = "Annual Precipitation (mm)",
      y = "Log-sinh Annual Streamflow",
      fill = bquote(atop(CO[2]~"ppm", "(Rainfall-Runoff Year)")), 
      colour = "Modelled Rainfall-Runoff Relationship",
      title = plot_label
    ) +
    scale_colour_viridis_d() +
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
      text = element_text(size = 9)
    ) +
    guides(
      colour = guide_legend(
        ncol = 2,
        override.aes = aes(linewidth = 1)
      ),
      fill = guide_colourbar(
        barwidth = unit(6, "cm")
      )
    )
  
  
  # Hard coded for main figure
  # add arrow is true
  if(plot_arrow) {
    # Arrow aes is hard coded based on the main figure plot
    # I could rescale it using max and min log-sinh flow and precip?
    if (plot_label == "B") {
      #y_axis_exclude <- element_blank()
      arrow_aes <- aes(x = 980, y = 180, xend = 980, yend = 60)
    } else {
      #y_axis_exclude <- element_text()
      arrow_aes <- aes(x = 1340, y = 350, xend = 1340, yend = 120)
    }
    
    plot <- plot + 
      geom_segment(
        mapping = arrow_aes,
        arrow = arrow(type = "closed", length = unit(2, "mm"))
      ) 
  }
  
  return(plot)
  

}





# Main figure plot -------------------------------------------------------------
main_plot_gauges <- c("606195", "238235")
  
part_a <- rainfall_runoff_plots(
  gauge = main_plot_gauges[1],
  model_years = c(1959, 1979, 1999, 2019), # this needs to be the min/max in calibration of gauge - write something here
  plot_label = "A",
  replace_precipitation = NULL,
  data = data,
  model_results = best_CO2_and_non_CO2_model_and_params_per_gauge,
  plot_arrow = TRUE,
  overwrite_legend_text = NULL
)

part_b <- rainfall_runoff_plots(
  gauge = main_plot_gauges[2],
  model_years = c(1959, 1979, 1999, 2019), # this needs to be the min/max in calibration of gauge - write something here
  plot_label = "B",
  replace_precipitation = NULL,
  data = data,
  model_results = best_CO2_and_non_CO2_model_and_params_per_gauge,
  plot_arrow = TRUE,
  overwrite_legend_text = NULL
)

illustration <- part_a | part_b 
final_illustration <- illustration +
  plot_layout(
    guides = "collect",
    axis_titles = "collect"
  ) &
  theme(legend.position = "bottom") 


ggsave(
  filename = "illustration_v2.pdf",
  plot = final_illustration,
  device = "pdf",
  path = "Figures/Main",
  width = 180,
  height = 130,
  units = "mm"
)






# Percentage change calculation ------------------------------------------------
## Mean annual rainfall ========================================================
mean_annual_rainfall <- data |> 
  filter(gauge %in% high_evi_gauges) |> 
  summarise(
    mean_annual_rainfall = mean(p_mm),
    max_annual_rainfall = max(p_mm),
    median_annual_rainfall = median(p_mm),
    .by = gauge
  ) 


## Get range of years for each gauge ===========================================
### It is the range of years in included_in_calibration
range_year_in_calibration <- data |> 
  filter(gauge %in% high_evi_gauges) |> 
  filter(included_in_calibration) |> 
  summarise(
    max_year = max(year),
    min_year = min(year),
    .by = gauge
  ) 


## Convert to pmap format for percentage change calculation ====================
info_for_percentage_change_pmap <- mean_annual_rainfall |> 
  left_join(
    range_year_in_calibration,
    by = join_by(gauge)
  ) |> 
  rowwise() |> 
  mutate(
    range_years = list(c(min_year, max_year))
  )

input_pmap <- list( # gauge, model_years, replace_precipitation
  info_for_percentage_change_pmap |> pull(gauge),
  info_for_percentage_change_pmap |> pull(range_years),
  info_for_percentage_change_pmap |> pull(mean_annual_rainfall)
)


## it is a pmap angle ==========================================================
percentage_changes <- pmap_dbl(
  .l = input_pmap,
  .f = get_percentage_difference,
  data = data,
  model_results = best_CO2_and_non_CO2_model_and_params_per_gauge
)


## Make into a nice tibble =====================================================
rainfall_runoff_percent_changes <- tibble(
  "gauge" = info_for_percentage_change_pmap |> pull(gauge),
  "rainfall_runoff_change" = percentage_changes
)


### print out changes in main figure (slightly diff because I am taking the included_in_calibration)
rainfall_runoff_percent_changes |> 
  filter(gauge %in% main_plot_gauges)




# Supplementary Plot -----------------------------------------------------------

## Chunk the gauges into groups of 8 ===========================================
### Ignore the two gauges in main
supp_high_evi_gauges <- high_evi_gauges[!high_evi_gauges %in% main_plot_gauges] |> 
  sort()
chunk <- 8
n <- length(supp_high_evi_gauges)
split_group <- rep(rep(1:ceiling(n / chunk), each = chunk))[1:n]
split_supp_high_evi_gauges <- split(supp_high_evi_gauges, split_group)


## Chunk range_years using the same method as gauges ===========================
split_supp_range_years <- split(
  info_for_percentage_change_pmap |> filter(!gauge %in% c(main_plot_gauges)) |> pull(range_years),
  split_group
)

## Labels ======================================================================
split_supp_letters <- map(
  .x = lengths(split_supp_high_evi_gauges),
  .f = \(x) rep(letters[1:x])
)


## Plot function ===============================================================
### Wrapper around rainfall_runoff_plots
make_grouped_rainfall_runoff_plot <- function(gauge_group, range_years_group, letters_group, replace_precipitation, data, model_results, plot_arrow, overwrite_legend_text) {
  
  input_plot_pmap <- list(
    gauge_group,
    range_years_group,
    letters_group
  )
  
  pmap(
    .l = input_plot_pmap,
    .f = rainfall_runoff_plots,
    replace_precipitation = replace_precipitation,
    data = data,
    model_results = model_results,
    plot_arrow = plot_arrow,
    overwrite_legend_text =overwrite_legend_text
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


## Repeat plot function for all chunks =========================================
pmap_input_supp_rainfall_runoff_plots <- list(
  split_supp_high_evi_gauges,
  split_supp_range_years,
  split_supp_letters
)

supp_rainfall_runoff_plots <- pmap(
  .l = pmap_input_supp_rainfall_runoff_plots,
  .f = make_grouped_rainfall_runoff_plot,
  replace_precipitation = NULL,
  data = data,
  model_results = best_CO2_and_non_CO2_model_and_params_per_gauge,
  plot_arrow = FALSE,
  overwrite_legend_text = c("Estimated Relationship Early", "Estimated Relationship Late")
)
  
  
## Make figure captions ========================================================
create_caption <- function(chunk_gauge, chunk_letters, identifier) {
  gauge_abc <- paste0(chunk_gauge, " (", chunk_letters, ")")
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

pmap_caption_input <- list(
  split_supp_high_evi_gauges,
  split_supp_letters,
  1:length(split_supp_high_evi_gauges)
)


sink("Figures/Other/supp_illustration_gauge_caption.txt")
pwalk(
  .l = pmap_caption_input,
  .f = create_caption
)
sink()


## Save supp_rainfall_runoff_plots =============================================
supp_rainfall_runoff_plots_names <- paste0("illustration_plot_", 1:length(supp_rainfall_runoff_plots), ".pdf")

walk2(
  .x = supp_rainfall_runoff_plots_names,
  .y = supp_rainfall_runoff_plots,
  .f = ggsave,
  path = "Figures/Other/",
  device = "pdf",
  width = 180,
  height = 254,
  units = "mm"
)






# Examine the observed rainfall-runoff ratio -----------------------------------









# Testing functions
x <- make_rainfall_runoff_results(
  gauge = "606195",
  model_years = c(1959, 1979, 1975, 1999, 2019),
  replace_precipitation = NULL,
  data = data,
  model_results = best_CO2_and_non_CO2_model_and_params_per_gauge
) 


y <- get_percentage_difference(
  gauge = "606195",
  model_years = c(1964, 2020), 
  replace_precipitation = 1043,
  data = data,
  model_results = best_CO2_and_non_CO2_model_and_params_per_gauge
) 


rainfall_runoff_plots(
  gauge = "606195",
  model_years = c(1964, 2020), 
  replace_precipitation = NULL,
  data = data,
  model_results = best_CO2_and_non_CO2_model_and_params_per_gauge,
  plot_label = "a",
  plot_arrow = FALSE,
  overwrite_legend_text = c("lower", "upper")
)












### Compare with CO2 model results #############################################
CO2_percentage_change_models <- read_csv(
  file = "Modelling/decade_streamflow_CO2_differences.csv",
  show_col_types = FALSE
) |> 
  select(gauge, decade, CO2_impact_on_streamflow_percent)

comparison_percentage_changes <- CO2_percentage_change_models |> 
  filter(decade == 2) |> # latest period (2012-2021)
  select(!decade) |> 
  right_join(
    illustration_percentage_changes,
    by = join_by(gauge)
    ) |> 
  mutate(
    abs_diff = abs(percent_change_illustration - CO2_impact_on_streamflow_percent)
  )

comparison_percentage_changes |> pull(abs_diff) |> summary()
comparison_percentage_changes |> ggplot(aes(y = abs_diff)) + geom_boxplot() + scale_y_log10() + theme_bw()
# median difference is 7 %.
# I am not using the best models. Nor am I comparing like decades.
# I think this is acceptable
# I am not sure how to graphically show this - add the runoff relationship stuff










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

### Find decades with little data ##############################################
### q_mm is the limiting factor
count_years_per_decade <- data |> 
  filter(gauge %in% high_evi_gauges) |> 
  mutate(
    decade = year - (year %% 10)
  ) |> 
  select(year, gauge, q_mm, decade) |> 
  drop_na() |> 
  summarise(
    n = n(),
    .by = c(gauge, decade)
  ) 

keep_decade_and_gauges <- count_years_per_decade |> 
  arrange(n) |> 
  # to make my life easier a decade needs at least 5 year of data
  filter(n >= 5)


### Plot #######################################################################
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
  # remove gauges and decades without enough years
  semi_join(
    keep_decade_and_gauges,
    by = join_by(gauge, decade)
  )

# there should be 684 - 114 = 534 left
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
  # remove gauges and decades without enough years
  semi_join(
    keep_decade_and_gauges,
    by = join_by(gauge, decade)
  ) |> 
  left_join(
    gauge_information,
    by = join_by(gauge)
  ) |> 
  summarise(
    mean_runoff_ratio_by_state = mean(runoff_ratio),
    .by = c(state, decade)
  )

# Percentage changes using our models by state
CO2_percentage_change_models |> 
  left_join(
    gauge_information,
    by = join_by(gauge)
  ) |> 
  filter(decade == 2) |> 
  summarise(
    mean_percent_change = median(CO2_impact_on_streamflow_percent),
    .by = state
  )


# percentage change
by_state_runoff_ratio |> 
  filter(state %in% c("NSW", "TAS", "VIC", "WA")) |>
  filter(decade %in% c(1960, 2010)) |> 
  pivot_wider(
    names_from = decade,
    values_from = mean_runoff_ratio_by_state
  ) |> 
  mutate(
    percent_diff = -(`1960` - `2010`) / `1960`
  )
  
  

# number of gauges by state
gauge_information |> 
  filter(gauge %in% high_evi_gauges) |> 
  count(state)

mean_decade_runoff_ratio_by_state <- by_state_runoff_ratio |> 
  filter(state %in% c("NSW", "TAS", "VIC", "WA")) |> # only 1 ACT, 1 QLD and 1 SA gauge
  ggplot(aes(x = decade, y = mean_runoff_ratio_by_state, colour = state)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Decade",
    y = "Runoff Ratio",
    colour = "State"
  ) +
  theme_bw()

ggsave(
  filename = "Figures/Other/change_in_runoff_ratio_by_state.pdf",
  plot = mean_decade_runoff_ratio_by_state,
  device = "pdf",
  width = 140,
  height = 100,
  units = "mm"
)


# Why does TAS have such a high ratio in the 1960s?
high_evi_tas_gauges <- gauge_information |> 
  filter(state == "TAS") |> 
  filter(gauge %in% high_evi_gauges) |> 
  pull(gauge)


check_TAS <- data |> 
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
  # remove gauges and decades without enough years
  semi_join(
    keep_decade_and_gauges,
    by = join_by(gauge, decade)
  ) |> 
  filter(gauge %in% high_evi_tas_gauges) |> 
  filter(decade == 1960)



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


# Add runoff ratio to comparison_percentage_changes ----------------------------
high_evi_ratio_ratio <- by_gauge_runoff_ratio |> 
  filter(gauge %in% high_evi_gauges) |> 
  drop_na() |> 
  select(decade, gauge, runoff_ratio) |> 
  group_by(gauge)

max_values <- slice_max(high_evi_ratio_ratio, decade, n = 1) |> 
  mutate(
    decade_name = "late"
  )

min_values <- slice_min(high_evi_ratio_ratio, decade, n = 1) |> 
  mutate(
    decade_name = "early"
  )

runoff_ratio_percent_changes <- rbind(
  max_values,
  min_values
) |> 
  select(!decade) |> 
  pivot_wider(
    names_from = decade_name,
    values_from = runoff_ratio
  ) |> 
  mutate(
    percentage_change_runoff_ratio = ((late - early) / early) * 100
  ) |> 
  arrange(gauge)


## Combine everything ==========================================================
all_comparison_percentage_changes <- runoff_ratio_percent_changes |> 
  select(gauge, percentage_change_runoff_ratio) |> 
  right_join(
    comparison_percentage_changes,
    by = join_by(gauge)
  ) |> 
  select(!abs_diff)


sign_all_comparison_percentage_changes <- all_comparison_percentage_changes |> 
  pivot_longer(
    cols = !gauge,
    names_to = "type",
    values_to = "percentage_change",
  ) |> 
  mutate(
    sign_percentage_change = sign(percentage_change)
  ) |> 
  summarise(
    sum_sign = sum(sign_percentage_change)
  ) |> 
  filter(abs(sum_sign) != 3)


wonky_gauges <- c("410713", "610008", "224213", "603007", "204036", "411003")
more_wonky_gauges <- sign_all_comparison_percentage_changes |> pull(gauge)
checking_gauges <- c(wonky_gauges, more_wonky_gauges) |> unique()

looking_at_issues <- all_comparison_percentage_changes |> 
  filter(gauge %in% checking_gauges)


all_comparison_percentage_changes |> 
  pivot_longer(
    cols = !gauge,
    names_to = "type",
    values_to = "percentage_change",
  ) |> 
  mutate(
    percentage_change = percentage_change / 100
  ) |> 
  mutate(
    type = case_when(
      type == "CO2_impact_on_streamflow_percent" ~ "Turn CO2 Off",
      type == "percentage_change_runoff_ratio" ~ "Runoff Ratio",
      type == "percent_change_illustration" ~ "Rainfall-runoff lines",
      .default = NA
    )
  ) |> 
  ggplot(aes(x = type, y = percentage_change)) +
  geom_boxplot(
    staplewidth = 0.5,
    outlier.alpha = 0.7,
    outlier.colour = "black",
    outlier.fill = "grey",
    outlier.shape = 21,
    whisker.linewidth = 0.2,
    box.linewidth = 0.2,
    median.linewidth = 0.4,
    staple.linewidth = 0.2,
    fill = "grey90"
  ) +
  labs(
    x = "Calculation Method",
    y = "Streamflow Percentage Change"
  ) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw()



