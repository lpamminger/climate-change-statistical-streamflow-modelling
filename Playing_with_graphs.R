## Fitting streamflow models to catchments
cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, truncnorm, sloop, ozmaps, sf, patchwork, ggmagnify, ggfx)



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
  show_col_types = FALSE
)

gauge_information_with_climate <- readr::read_csv(
  "./Data/Tidy/gauge_information_CAMELS_with_climate.csv",
  show_col_types = FALSE
)

CMAES_results <- read_csv("./Results/my_cmaes/CMAES_parameter_results_20250331.csv",
                          show_col_types = FALSE
) 

best_CO2_non_CO2_per_gauge <- read_csv(
  "./Results/my_cmaes/unmodified_best_CO2_non_CO2_per_catchment_CMAES_20250331.csv",
  show_col_types = FALSE
) 

streamflow_data <- read_csv(
  "Results/my_cmaes/CMAES_streamflow_results_20250331.csv",
  show_col_types = FALSE
)

## Utility functions ===========================================================
source("./Functions/utility.R")


## Import streamflow functions =================================================
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



# get gauge with a5 near middle of t-series
select_gauge <- best_CO2_non_CO2_per_gauge |> 
  filter(parameter == "a5") |> 
  arrange(parameter_value)

# gauge selected is 603004 
test_gauge <- "603004"

# ToE
a5_to_time_of_emergence <- function(a5, CO2, year) {
  adjusted_CO2 <- if_else(CO2 - a5 < 0, 0, CO2 - a5)
  
  year_where_CO2_impacts_flow <- year[adjusted_CO2 != 0][1]
  
  return(year_where_CO2_impacts_flow)
}

ToE <- a5_to_time_of_emergence(
  a5 = 75.699999995,
  CO2 = data |> filter(gauge == test_gauge) |> pull(CO2),
  year = data |> filter(gauge == test_gauge) |> pull(year)
)

test_gauge_data <- best_CO2_non_CO2_per_gauge |> 
  filter(gauge == test_gauge)

best_models <- test_gauge_data |> 
  pull(streamflow_model) |> 
  unique()


observed_data <- data |> 
  filter(gauge == test_gauge) |> 
  filter(included_in_calibration) |> 
  select(year, p_mm, bc_q) |> 
  rename(
    precipitation = p_mm,
    boxcox_streamflow = bc_q
  ) |> 
  add_column(
    metric = "observed"
  )


# convert back to streamflow (no boxcox)
bc_lambda <- gauge_information |> 
  filter(gauge == test_gauge) |> 
  pull(bc_lambda)

streamflow_gauge_data <- streamflow_data |> 
  filter(gauge == test_gauge) |> 
  filter(year %in% observed_data$year) |> 
  filter(streamflow_model %in% best_models) |> 
  select(year, precipitation, modelled_boxcox_streamflow, streamflow_model) |> 
  rename(
    metric = streamflow_model,
    boxcox_streamflow = modelled_boxcox_streamflow
  ) |> 
  rbind(observed_data) |> 
  mutate(
    # back to regular streamflow
    streamflow = boxcox_inverse_transform(
      boxcox_streamflow, 
      lambda = bc_lambda, 
      lambda_2 = 1
    )
  )




name_title <- gauge_information |> 
  filter(gauge == test_gauge) |>
  select(station_name, state) |> 
  as.character()

streamflow_timeseries <- streamflow_gauge_data |> 
  ggplot(
    aes(
      x = year, 
      y = streamflow, 
      colour = metric
    )
  ) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = ToE) +
  theme_bw() +
  labs(
    title = paste0(name_title[1], " ", name_title[2])
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )

streamflow_timeseries
# This is an interesting graph
# Show the optimization worked quite well.
rainfall_runoff_plot <- streamflow_gauge_data |> 
  ggplot(
    aes(
      x = precipitation,
      y = boxcox_streamflow,
      colour = metric,
      shape = metric,
      fill = metric
    )
  ) +
  geom_point() +
  geom_line(
    stat = "smooth",
    method = "lm",
    formula = y ~ x,
    alpha = 0.4,
    linewidth = 1
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )


# Turn off a3 term - compare t-series
test_catchment_data <- catchment_data_blueprint(
  gauge_ID = test_gauge,
  observed_data = data,
  start_stop_indexes = start_stop_indexes
)

parameters <- test_gauge_data |> 
  filter(streamflow_model == "streamflow_model_intercept_shifted_CO2") |> 
  pull(parameter_value)

# streamflow model is not clever enough to pull the necessary bits
# from catchment data. Use must give it the tibble.
# Names are confusing. 
# Maybe have a function that detects if catchment data is catchment_data object
# if so then pull either the start stop or entire tibble and remove NAs

a3_off_parameters <- parameters
a3_off_parameters[3] <- 0

a3_off <- streamflow_model_intercept_shifted_CO2(
  catchment_data = test_catchment_data$stop_start_data_set$start_index,
  parameter_set = a3_off_parameters
)

# Turn off a5 term - compare t-series
a5_off_parameters <- parameters
a5_off_parameters[4] <- 1000 # larger than max CO2 value 
a5_off <- streamflow_model_separate_shifted_CO2(
  catchment_data = test_catchment_data$stop_start_data_set$start_index,
  parameter_set = a5_off_parameters
)




parameters_off <- observed_data |> 
  select(year, precipitation) |> 
  add_column("Climate change component off" = a3_off) |> 
  #add_column("a5_off" = a5_off) |> 
  pivot_longer(
    cols = ends_with("off"),
    names_to = "metric",
    values_to = "boxcox_streamflow"
  ) |> 
  mutate(
    # back to regular streamflow
    streamflow = boxcox_inverse_transform(
      boxcox_streamflow, 
      lambda = bc_lambda, 
      lambda_2 = 1
    )
  ) |> 
  rbind(streamflow_gauge_data) |> 
  filter(metric != "streamflow_model_precip_seasonal_ratio_auto") |> 
  filter(year %in% observed_data$year) 


parameter_off_plot <- parameters_off |> 
  ggplot(
    aes(
      x = year,
      y = streamflow,
      colour = metric
    )
  ) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = ToE) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )




# Difference in differences?
# separate data set into impacted by CO2 and not impacted by CO2
diff_in_diff <- streamflow_gauge_data |> # this has the NA's removed
  filter(metric == "observed") |> 
  mutate(
    CO2_impacted = year > ToE
  )

stats_diff <- diff_in_diff |> 
  summarise(
    mean_boxcox_streamflow = mean(boxcox_streamflow),
    n = n(),
    .by = CO2_impacted
  ) |> 
  mutate(
    n_label = paste0("n = ", n)
  )


diff_in_diff |> 
  ggplot(aes(x = boxcox_streamflow)) +
  geom_histogram(bins = 20) +
  facet_wrap(~CO2_impacted)

streamflow_ToE_comparison <- diff_in_diff |> 
  ggplot(aes(x = CO2_impacted, y = streamflow)) +
  geom_boxplot(
    staplewidth = 0.25
  ) +
  geom_text(
    data = stats_diff,
    aes(x = CO2_impacted, y = 0, label = n_label)
  ) +
  theme_bw()


# Time of emergence of best histograms
best_model_per_gauge <- best_CO2_non_CO2_per_gauge |>
  slice_min(
    AIC,
    by = gauge
  ) |> 
  distinct()


# time of emergence a5 parameter
time_of_emergence <- best_model_per_gauge |> 
  filter(parameter == "a5") 


repeat_a5_to_time_of_emergence <- function(gauge, a5_per_gauge, data) {
  
  # get a5 from a5_per_gauge
  a5_parameter <- a5_per_gauge |> 
    filter(gauge == {{ gauge }}) |> 
    pull(parameter_value)
  
  
  # get ToE
  a5_to_time_of_emergence(
    a5 = a5_parameter,
    CO2 = data |> filter(gauge == {{ gauge }}) |> pull(CO2),
    year = data |> filter(gauge == {{ gauge }}) |> pull(year)
  )
}

year_ToE <- map_dbl(
  .x = time_of_emergence |> pull(gauge),
  .f = repeat_a5_to_time_of_emergence,
  a5_per_gauge = time_of_emergence,
  data = data
)

gauge_state_climate <- gauge_information_with_climate |> 
  select(gauge, state, major_climate_type)

time_of_emergence_with_years <- time_of_emergence |> 
  add_column(year_ToE) |> 
  left_join(
    gauge_state_climate,
    by = join_by(gauge)
  )
  

## Create my own bins 

min_CO2 <- data |> pull(CO2) |> min() 

custom_bins_tim_of_emergence <- time_of_emergence_with_years |> 
  mutate(
    custom_bins = case_when( # I might switch to between?
      between(year_ToE, 1950, 1959) ~ "1950's", # entries start at 1959
      between(year_ToE, 1960, 1969) ~ "1960's",
      between(year_ToE, 1970, 1979) ~ "1970's",
      between(year_ToE, 1980, 1989) ~ "1980's",
      between(year_ToE, 1990, 1999) ~ "1990's",
      between(year_ToE, 2000, 2009) ~ "2000's",
      between(year_ToE, 2010, 2019) ~ "2010's",
      between(year_ToE, 2020, 2029) ~ "2020's",
      #(year_ToE >= 1959) & (year_ToE < 1964) ~ "1959 to 1964",
      #(year_ToE >= 1964) & (year_ToE < 1969) ~ "1964 to 1969",
      #(year_ToE >= 1969) & (year_ToE < 1974) ~ "1969 to 1974",
      #(year_ToE >= 1974) & (year_ToE < 1979) ~ "1974 to 1979",
      #(year_ToE >= 1979) & (year_ToE < 1984) ~ "1979 to 1984",
      #(year_ToE >= 1984) & (year_ToE < 1989) ~ "1984 to 1989",
      #(year_ToE >= 1989) & (year_ToE < 1994) ~ "1989 to 1994",
      #(year_ToE >= 1994) & (year_ToE < 1999) ~ "1994 to 1999",
      #(year_ToE >= 1999) & (year_ToE < 2004) ~ "1999 to 2004",
      #(year_ToE >= 2004) & (year_ToE < 2009) ~ "2004 to 2009",
      #(year_ToE >= 2009) & (year_ToE < 2014) ~ "2009 to 2014",
      #(year_ToE >= 2014) & (year_ToE < 2019) ~ "2014 to 2019",
      .default = NULL
    )
  ) |> 
  mutate(
    custom_bins = if_else(parameter_value < min_CO2, "Before 1959", custom_bins)
  ) 

# not used
custom_bins_by_state <- custom_bins_tim_of_emergence |> 
  summarise(
    n = n(),
    .by = c(custom_bins, state)
  ) |> 
  mutate(
    custom_bins = factor(
      custom_bins, 
      levels = c(
        "Before 1959",
        paste0(seq(from = 1950, to = 2020, by = 10), "'s")
      )
      )
  )


# by major climate type - not used
custom_bins_by_climate <- custom_bins_tim_of_emergence |> 
  summarise(
    n = n(),
    .by = c(custom_bins, major_climate_type)
  ) |> 
  mutate(
    custom_bins = factor(
      custom_bins, 
      levels = c(
        "Before 1959",
        paste0(seq(from = 1950, to = 2020, by = 10), "'s")
      )
    )
  )


# Doing it by state or climate does not really tell us much
# by state 
ToE_hist <- custom_bins_tim_of_emergence |> 
  summarise(
    n = n(),
    .by = custom_bins
  ) |> 
  mutate(
    custom_bins = factor(
      custom_bins, 
      levels = c(
        "Before 1959",
        paste0(seq(from = 1950, to = 2020, by = 10), "'s")
      )
    )
  ) |> 
  ggplot(aes(x = custom_bins, y = n)) + 
  geom_col(colour = "black", fill = "lightblue") +
  labs(
    x = "Decade",
    y = "Count",
    title = bquote("What decade does"~CO[2]~"start impacting annual rainfall-partitioning?")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20)
    )
ToE_hist
#ggsave(
#  filename = "ToE_graphs.pdf",
#  plot = ToE_hist,
#  device = "pdf",
#  path = "./Graphs",
#  width = 297,
#  height = 210,
#  units = "mm"
#)


## optimal binwidth calculation

bw <- 2 * IQR(year_ToE) / length(year_ToE)^(1/3)


hist_ToE <- time_of_emergence_with_years |> 
  ggplot(aes(x = year_ToE)) +
  geom_histogram(
    binwidth = 1, # which binwidth to choose
    fill = "lightgrey",
    colour = "black"
    ) +
  theme_bw()


# list
graphs <- list(
  streamflow_timeseries, 
  rainfall_runoff_plot, 
  parameter_off_plot, 
  streamflow_ToE_comparison,
  hist_ToE
  )


#ggsave(
#  filename = "playing_with_graphs.pdf",
#  plot = gridExtra::marrangeGrob(
#    graphs, 
#    nrow = 1, 
#    ncol = 1,
#    top = NULL
#  ),
#  # remove page numbers
#  device = "pdf",
#  path = "./Graphs",
#  width = 297,
#  height = 210,
#  units = "mm"
#)



################
### NEW WORK ###
################





# Map of binned ToE ------------------------------------------------------------
## Making map ==================================================================
aus_map <- ozmaps::ozmap(x = "states") |> 
  filter(!NAME %in% c("Other Territories")) |> 
  rename(state = NAME) |> 
  mutate(
    state = case_when(
      state == "New South Wales" ~ "NSW",
      state == "Victoria" ~ "VIC",
      state == "Queensland" ~ "QLD",
      state == "South Australia" ~ "SA",
      state == "Western Australia" ~ "WA",
      state == "Tasmania" ~ "TAS",
      state == "Northern Territory" ~ "NT",
      state == "Australian Capital Territory" ~ "ACT",
    )
  )

combine_NSW_ACT <- aus_map |> 
  filter(state %in% c("NSW", "ACT")) |> 
  st_union()

aus_map[1,2] <- list(combine_NSW_ACT)

aus_map <- aus_map |> 
  filter(state != "ACT")


## Add lon and lat to time_of_emergence_with_years
lat_long_gauge <- gauge_information |> 
  select(gauge, lat, lon)

ToE_map_info <- time_of_emergence_with_years |> 
  left_join(
    lat_long_gauge,
    by = join_by(gauge)
  )


# I require custom bins because of the 1959 things
min_CO2 <- data |> pull(CO2) |> min() 

binned_ToE_map_info <- ToE_map_info |> 
  mutate(
    custom_bins = case_when( # I might switch to between?
      between(year_ToE, 1950, 1959) ~ "1950's", # entries start at 1959
      between(year_ToE, 1960, 1969) ~ "1960's",
      between(year_ToE, 1970, 1979) ~ "1970's",
      between(year_ToE, 1980, 1989) ~ "1980's",
      between(year_ToE, 1990, 1999) ~ "1990's",
      between(year_ToE, 2000, 2009) ~ "2000's",
      between(year_ToE, 2010, 2019) ~ "2010's",
      between(year_ToE, 2020, 2029) ~ "2020's",
      .default = NULL
    )
  ) |> 
  mutate(
    custom_bins = if_else(parameter_value < min_CO2, "Before 1959", custom_bins)
  ) |> 
  mutate(
    custom_bins = factor(
      custom_bins, 
      levels = c(
        "Before 1959",
        paste0(seq(from = 1960, to = 2020, by = 10), "'s")
      )
    )
  )

## Add record length to ToE data to see if ToE is impacted by record length
record_length_per_gauge <- gauge_information |> 
  select(gauge, record_length) |> 
  mutate(
    record_length_binned = case_when(
      between(record_length, 30, 39) ~ "30-39",
      between(record_length, 40, 49) ~ "40-49",
      between(record_length, 50, 59) ~ "50-59",
      between(record_length, 60, 69) ~ "60-69",
      .default = NA
    )
  )

binned_ToE_map_info <- binned_ToE_map_info |> 
  left_join(
    record_length_per_gauge,
    by = join_by(gauge)
  ) |> 
  mutate(
    record_length_binned = factor(record_length_binned, levels = c("30-39", "40-49", "50-59", "60-69"))
  )


## Map 
ggplot() +
  geom_sf(
    data = aus_map,
    colour = "black",
    fill = "grey50"
  ) +
  geom_point(
    data = binned_ToE_map_info,
    aes(x = lon, y = lat, colour = custom_bins, shape = record_length_binned), #shape = flag),
    size = 0.75
  ) 

# Do all states - copy 03 code

gg_map_plot <- function(state, ToE_data, map_data) {
  
  state_map_data <- map_data |> 
    filter(state == {{ state }})
  
  state_ToE_data <- ToE_data |> 
    filter(state == {{ state }})
  
  ggplot() +
    geom_sf(
      data = state_map_data,
      colour = "black",
      fill = "grey50"
    ) +
    geom_point(
      data = state_ToE_data,
      aes(x = lon, y = lat, colour = custom_bins, shape = record_length_binned), 
      size = 0.75,
      show.legend = TRUE # Required to show ever factor in legend
    ) +
    scale_colour_brewer(
      palette = "RdYlBu",
      drop = FALSE # Required to show ever factor in legend
    ) +
    scale_shape_discrete(
      drop = FALSE
    ) +
    labs(
      colour = "Time of Emergence",
      x = NULL,#"Longitude",
      y = NULL,#"Latitude",
      shape = "Record Length"
    ) +
    theme_bw() +
    theme(
      legend.key = element_rect(fill = "grey50"),
      legend.background = element_rect(colour = "black"),
      legend.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 4)
    ) +
    guides(colour = guide_legend(nrow = 2)) # Wrap legend
  
  
}


by_state_plots <- map(
  .x = aus_map |> pull(state) |> unique(),
  .f = gg_map_plot,
  ToE_data = binned_ToE_map_info,
  map_data = aus_map
) |> 
  `names<-`(aus_map |> pull(state) |> unique()) #same as names(by_state_plots) <- aus_map |> pull(state) |> unique()

# The sf plots are a fixed aspect ratio
# Changing widths and heights in plot_layout will not do anything --> https://patchwork.data-imaginist.com/articles/guides/layout.html#fixed-aspect-plots

fr_top <- by_state_plots[["TAS"]] | by_state_plots[["QLD"]] | by_state_plots[["NT"]] | by_state_plots[["SA"]] 
fr_bottom <- by_state_plots[["WA"]] | by_state_plots[["VIC"]] | by_state_plots[["NSW"]]
ToE_record_length_nice_plot <- fr_top / fr_bottom / guide_area() + plot_layout(guides = "collect")
ToE_record_length_nice_plot



ToE_vs_record_length <- binned_ToE_map_info |> 
  ggplot(aes(x = record_length, y = year_ToE)) +
  geom_jitter() +
  labs(
    x = "Record Length",
    y = "Time of Emergence"
  ) +
  theme_bw()

ToE_vs_record_length

binned_ToE_map_info |> 
  ggplot(aes(x = as.factor(record_length_binned), y = as.factor(custom_bins))) +
  geom_jitter() +
  labs(
    x = "Record Length",
    y = "Time of Emergence"
  ) +
  theme_bw()


# Components of best models ----------------------------------------------------
## Histogram of models =========================================================
best_model_per_gauge |> 
  select(gauge, streamflow_model) |> 
  distinct() |> 
  ggplot(aes(x = streamflow_model)) +
  geom_bar() +
  labs(
    x = NULL,
    y = "Count"
    ) +
  theme_bw() +
  theme(
    axis.text = element_text(angle = 90)
  )


## Histogram of model components ===============================================
state_info <- gauge_information |> 
  select(gauge, state)


best_model_components <- best_model_per_gauge |> 
  left_join(
    state_info,
    by = join_by(gauge)
  ) |> 
  filter(!parameter %in% c("a0", "a1", "sd")) |> 
  ggplot(aes(x = parameter, fill = state)) +
  geom_bar() +
  theme_bw() +
  labs(
    x = "Parameter",
    y = "Count",
    fill = "State",
    title = "What is the most used parameter of the best streamflow models?",
    subtitle = "(excluding a0, a1 and sd)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )



## Map plot of best models =====================================================
streamflow_map <- best_model_per_gauge |> 
  select(gauge, streamflow_model) |> 
  distinct() |> 
  left_join(
    lat_long_gauge,
    by = join_by(gauge)
  ) |> 
  left_join(
    state_info,
    by = join_by(gauge)
  ) |> 
  mutate(
    streamflow_model = as.factor(streamflow_model)
  )


  


ggplot() +
  geom_sf(
    data = aus_map,
    colour = "black",
    fill = "grey50"
  ) +
  geom_point(
    data = streamflow_map,
    aes(x = lon, y = lat, colour = streamflow_model), 
    size = 0.75
  ) 

gg_streamflow_map_plot <- function(state, ToE_data, map_data) {
  
  state_map_data <- map_data |> 
    filter(state == {{ state }})
  
  state_ToE_data <- ToE_data |> 
    filter(state == {{ state }})
  
  ggplot() +
    geom_sf(
      data = state_map_data,
      colour = "black",
      fill = "grey50"
    ) +
    geom_point(
      data = state_ToE_data,
      aes(x = lon, y = lat, colour = streamflow_model), 
      size = 0.75,
      show.legend = TRUE # Required to show ever factor in legend
    ) +
    scale_color_discrete(
      drop = FALSE
    ) +
    labs(
      colour = "Streamflow model",
      x = NULL,#"Longitude",
      y = NULL#"Latitude"
    ) +
    theme_bw() +
    theme(
      legend.key = element_rect(fill = "grey50"),
      legend.background = element_rect(fill = "grey50", colour = "black"),
      legend.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 4),
      legend.text = element_text(size = 6)
    ) +
    guides(colour = guide_legend(nrow = 4)) # Wrap legend
  
  
}


by_state_model_plots <- map(
  .x = aus_map |> pull(state) |> unique(),
  .f = gg_streamflow_map_plot,
  ToE_data = streamflow_map,
  map_data = aus_map
) |> 
  `names<-`(aus_map |> pull(state) |> unique()) #same as names(by_state_plots) <- aus_map |> pull(state) |> unique()

model_top <- by_state_model_plots[["TAS"]] | by_state_model_plots[["QLD"]] | by_state_model_plots[["NT"]] | by_state_model_plots[["SA"]] 
model_bottom <- by_state_model_plots[["WA"]] | by_state_model_plots[["VIC"]] | by_state_model_plots[["NSW"]]
model_nice_plot <- (model_top / model_bottom / guide_area()) + plot_layout(guides = "collect")

model_nice_plot


ggsave(
  filename = "ToE_nice_plot.pdf",
  plot = fr_nice_plot,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)

#ggsave(
#  filename = "streamflow_model_nice_plot.pdf",
#  plot = model_nice_plot,
#  device = "pdf",
#  width = 297,
#  height = 210,
#  units = "mm"
#)



# Component maps ---------------------------------------------------------------
# method: make tibble with |model component|australia_shape_file| -> facet_wrap|. 
# Use the best_model tibble. 
#Remove a0_d, a5. Replace a2 with auto, replace a3 with CO2, a0_n with drought and a4 with seasonal. 
#Assign model component and lat long. 
# Add as other layer 

# Map
single_aus_map <- ozmaps::ozmap("country") |> 
  uncount(4) |>  # repeat the geometry 4 times
  mutate(
    NAME = c("Drought", "Autocorrelation", "CO2", "Rainfall Seasonality")
  )

# Best model per gauge adjustment
simplified_components_best_model_per_gauge <- best_model_per_gauge |> 
  filter(parameter %in% c("a2", "a3_intercept", "a3_slope", "a4", "a0_d")) |> 
  select(gauge, streamflow_model, parameter) |> 
  mutate(
    NAME = case_when(
      parameter == "a0_d" ~ "Drought",
      parameter == "a2" ~ "Autocorrelation",
      parameter %in% c("a3_intercept", "a3_slope") ~ "CO2",
      parameter == "a4" ~ "Rainfall Seasonality",
      .default = NA
    )
  ) |> 
  left_join(
    lat_long_gauge,
    by = join_by(gauge)
  )

# Count of components
total_gauges <- best_model_per_gauge |> pull(gauge) |> unique() |> length()
count_components <- simplified_components_best_model_per_gauge |> 
  summarise(
    n = n(),
    .by = NAME
  ) |> 
  mutate(
    label = paste0(round((n / total_gauges) * 100, digits = 2), "%", " (n = ", n, ")")
  ) |> 
  add_column(
    lon = 160
  ) |> 
  add_column(
    lat = -15
  )

# Put it together
model_components <- single_aus_map |> 
ggplot(aes(geometry = geometry)) +
  geom_sf(
    colour = "black",
    fill = "grey80"
  ) +
  geom_point(
    mapping = aes(x = lon, y = lat, colour = NAME),
    data = simplified_components_best_model_per_gauge,
    inherit.aes = FALSE,
    size = 0.75,
    show.legend = FALSE
  ) +
  geom_text(
    mapping = aes(x = lon, y = lat, label = label),
    data = count_components,
    inherit.aes = FALSE
  ) +
  metR::scale_x_longitude(ticks = 5) +
  metR::scale_y_latitude(ticks = 5) +
  labs(
    x = "Longitude",
    y = "Latitude",
    colour = NULL
  ) +
  scale_colour_brewer(palette = "Set1") +
  facet_wrap(~NAME) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  ) 

model_components
# transform axis to lat/long degrees and that?

# Bin evidence ratios for qualitative plots---------------------------------------
# Copied from 03
evidence_ratio_calc <- best_CO2_non_CO2_per_gauge |>
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
  arrange(evidence_ratio) 


### Get into plot ready form ###################################################
lat_long_evidence_ratio <- evidence_ratio_calc |>
  select(!c(AIC_difference)) |>
  left_join(
    lat_long_gauge,
    by = join_by(gauge)
  ) 


# Weak or limited -> 1-10
# Moderate -> 10-100
# Moderately strong -> 100-1,000
# Strong -> 1,000-10,000
# Very strong -> 10,000-1,000,000
# Extremely strong ->  >1,000,000

binned_lat_lon_evidence_ratio <- lat_long_evidence_ratio |> 
  mutate(
    binned_evidence_ratio = case_when(
      between(evidence_ratio, -1E1, 1E1) ~ "Weak",
      between(evidence_ratio,  1E1, 1E2) ~ "Moderate",
      between(evidence_ratio,  1E2, 1E3) ~ "Moderately Strong",
      between(evidence_ratio,  1E3, 1E4) ~ "Strong",
      between(evidence_ratio,  1E4, 1E6) ~ "Very Strong",
      between(evidence_ratio,  1E6, Inf) ~ "Extremely Strong",
      .default = NA
    )
  ) |> 
  left_join(
    state_info,
    by = join_by(gauge)
  ) |> 
  mutate(
    binned_evidence_ratio = factor(
      binned_evidence_ratio, 
      levels = c("Weak", "Moderate", "Moderately Strong", "Strong", "Very Strong", "Extremely Strong")
      )
  )


# Add the a3 direction here
direction_of_a3_change <- best_model_per_gauge |> 
  filter(parameter %in% c("a3_intercept", "a3_slope")) |>
  mutate(
    intercept_or_slope = if_else(str_detect(streamflow_model, "intercept"), "Intercept", "Slope")
  ) |> 
  select(gauge, streamflow_model, parameter, parameter_value, intercept_or_slope) |> 
  mutate(
    CO2_direction = if_else(parameter_value < 0, "Negative", "Positive")
  ) |> 
  select(gauge, CO2_direction, intercept_or_slope)


a3_direction_binned_lat_lon_evidence_ratio <- binned_lat_lon_evidence_ratio |> 
  left_join(
    direction_of_a3_change,
    by = join_by(gauge)
  ) |> 
  replace_na(list(CO2_direction = "No CO2 Term", intercept_or_slope = "No CO2 Term")) |> 
  unite(
    col = "impact_of_CO2_term", 
    CO2_direction,
    intercept_or_slope,
    sep = "-"
  ) |> 
  mutate(
    impact_of_CO2_term = if_else(impact_of_CO2_term == "No CO2 Term-No CO2 Term", "No CO2 Term", impact_of_CO2_term)
  ) |> 
  mutate(
    impact_of_CO2_term = factor(
      impact_of_CO2_term, 
      levels = c("No CO2 Term", "Negative-Intercept", "Positive-Intercept", "Negative-Slope", "Positive-Slope")
      )
  )








# Custom colour palette
custom_palette <- function(x) {
  rev(c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7","#f7f7f7"))
}
  
# function here
gg_evidence_ratio_map <- function(state = NULL, evidence_ratio_data, map_data) {
  
  # use state as key to extract correct polygon from map_data and evidence_ratio_data
  if (!is.null(state)) {
    
    map_data <- map_data |> 
      filter(state == {{ state }})
    
    evidence_ratio_data <- evidence_ratio_data |> 
      filter(state == {{ state }})
    
  }
  
  
  ggplot() +
    geom_sf(
      data = map_data,
      colour = "black",
      fill = "grey80"
    ) +
    geom_point(
      data = evidence_ratio_data,
      aes(x = lon, y = lat, fill = binned_evidence_ratio, shape = impact_of_CO2_term), 
      size = 1,
      show.legend = TRUE,
      colour = "black",
      stroke = 0.1
    ) +
    scale_fill_manual(
      values = custom_palette(),
      drop = FALSE
    ) +
    scale_shape_manual(
      values = c(21, 22, 23, 24, 25),
      drop = FALSE
    ) +
    theme_bw() +
    labs(
      x = NULL, 
      y = NULL, 
      fill = "Evidence Ratio",
      shape = "Impact of CO2 Term"
    ) +
    theme(
      legend.key = element_rect(fill = "grey80"),
      legend.title = element_text(hjust = 0.5),
      legend.background = element_rect(colour = "black"),
      axis.text = element_text(size = 6)
    ) +
    guides(
      fill = guide_legend(override.aes = list(size = 5, shape = 21), nrow = 3), # Wrap legend with nrow
      shape = guide_legend(override.aes = list(size = 5))
      ) 
  
  
}


state_evidence <- map(
  .x = aus_map |> pull(state),
  .f = gg_evidence_ratio_map,
  evidence_ratio_data = a3_direction_binned_lat_lon_evidence_ratio,
  map_data = aus_map
) |> 
  `names<-`(aus_map |> pull(state))



evi_top <- state_evidence[["TAS"]] | state_evidence[["QLD"]] | state_evidence[["NT"]] | state_evidence[["SA"]] 
evi_bottom <- state_evidence[["WA"]] | state_evidence[["VIC"]] | state_evidence[["NSW"]]
evi_nice_plot <- (evi_top / evi_bottom / guide_area()) + plot_layout(guides = "collect")
evi_nice_plot

ggsave(
  filename = "Updated_Evi_Plot.pdf",
  plot = evi_nice_plot,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)



# Other potential information - slope vs. intercept ----------------------------
# some sort of slope and intercept analysis. Like proportion of pos/neg slope,
# pos/neg intercept

y <- a3_direction_binned_lat_lon_evidence_ratio |> 
  count(intercept_or_slope, state) |> 
  pivot_wider(
    names_from = intercept_or_slope,
    values_from = n
  ) |> 
  mutate(
    total = sum(across(Intercept:`No CO2 Term`)),
    .by = state
  ) |> 
  pivot_longer(
    cols = Intercept:`No CO2 Term`,
    names_to = "intercept_or_slope",
    values_to = "n"
  ) |> 
  relocate(
    total,
    .after = 4
  ) |> 
  mutate(
    fraction = n / total
  ) |> 
  mutate(
    intercept_or_slope = factor(intercept_or_slope, levels = c("No CO2 Term", "Slope", "Intercept"))
  )

plot_label <- y |> 
  select(state, total) |> 
  distinct() |> 
  mutate(plot_label = paste0("n = ", total)) |> 
  add_column(y_pos = 1)

y |> 
  ggplot(aes(x = state, y = fraction, fill = intercept_or_slope)) +
  geom_col() +
  geom_text(
    data = plot_label,
    aes(x = state, y = y_pos, label = plot_label),
    inherit.aes = FALSE,
    nudge_y = 0.01
  ) +
  labs(
    x = "State",
    y = "Fraction",
    fill = "CO2 model component"
  ) +
  theme_bw()

## map of intercept vs. slope - shapes of colour?
gg_temp_intercept_slope_map <- function(state = NULL, evidence_ratio_data, map_data) {
  
  # use state as key to extract correct polygon from map_data and evidence_ratio_data
  if (!is.null(state)) {
    
    map_data <- map_data |> 
      filter(state == {{ state }})
    
    evidence_ratio_data <- evidence_ratio_data |> 
      filter(state == {{ state }})
    
  }
  
  
  ggplot() +
    geom_sf(
      data = map_data,
      colour = "black",
      fill = "grey50"
    ) +
    geom_point(
      data = evidence_ratio_data,
      aes(x = lon, y = lat, colour = binned_evidence_ratio, shape = intercept_or_slope), #shape = flag),
      size = 0.75,
      show.legend = TRUE
    ) +
    scale_colour_manual(
      values = custom_palette(),
      drop = FALSE
    ) +
    theme_bw() +
    labs(
      x = NULL, # commented out due for final review plot
      y = NULL, # commented out due for final review plot
      colour = "Evidence Ratio",
      shape = "Type of CO2 change"
    ) +
    theme(
      legend.key = element_rect(fill = "grey50"),
      legend.title = element_text(hjust = 0.5),
      legend.background = element_rect(colour = "black"),
      axis.text = element_text(size = 6)
    ) +
    guides(
      colour = guide_legend(override.aes = list(size = 5), nrow = 3), # Wrap legend
      shape = guide_legend(override.aes = list(size = 5))
    ) 
  
  
}


state_intercept_and_slope <- map(
  .x = aus_map |> pull(state),
  .f = gg_temp_intercept_slope_map,
  evidence_ratio_data = a3_direction_binned_lat_lon_evidence_ratio,
  map_data = aus_map
) |> 
  `names<-`(aus_map |> pull(state))



type_top <- state_intercept_and_slope[["TAS"]] | state_intercept_and_slope[["QLD"]] | state_intercept_and_slope[["NT"]] | state_intercept_and_slope[["SA"]] 
type_bottom <- state_intercept_and_slope[["WA"]] | state_intercept_and_slope[["VIC"]] | state_intercept_and_slope[["NSW"]]
type_nice_plot <- (type_top / type_bottom / guide_area()) + plot_layout(guides = "collect")
type_nice_plot


new_graphs <- list(model_components, evi_nice_plot, ToE_vs_record_length, ToE_record_length_nice_plot, type_nice_plot)


my_ggsave <- function(plot, name) {
  ggsave(
    filename = paste0(name,".pdf"),
    plot = plot,
    device = "pdf",
    path = "./Graphs",
    width = 297,
    height = 210,
    units = "mm"
  )
}

walk2(
  .x = new_graphs,
  .y = c("model_components", "evi_nice_plot", "ToE_vs_record_length", "ToE_record_length_nice_plot", "type_nice_plot"),
  .f = my_ggsave
)


# Compare 3 methods to estimate the impact of CO2 on streamflow ----------------
## Rank all the models by gauge (CMAES results)
ranked_models_per_gauge <- CMAES_results |> 
  group_by(gauge) |> 
  arrange(AIC) |> 
  ungroup()

## Comparing the three methods 
### Method 1: set a3 to zero
### Method 2: compare best CO2 to best non-CO2
### Method 3: compare best CO2 to equivalent non-CO2 model (could also be best non-CO2)

### How to compare the methods:
### 1. streamflow-time with observed
### 2. rainfall-runoff with observed
### 3. Difference between models vs. evidence ratio - not sure if this will work
###    with method 1 - I could do difference between observed - modelled vs. evi ratio

## Method 2 ====================================================================
### Evidence ratio #############################################################
evidence_ratio_calc <- best_CO2_non_CO2_per_gauge |>
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
  arrange(evidence_ratio) #|> 


### Difference in streamflow models ############################################
### Difference in streamflow models per year mean(mm/y) vs. evidence ratio or AIC
### Use best_CO2_non_CO2_per_gauge to filter streamflow_data
### Pivot wider, mutate model 1 - model 2 -> summarise mean .by gauge
### y axis-options:
# mean absolute difference in streamflow (mm/y)
# total difference in streamflow over the entire period
# mean percentage difference in streamflow per year

only_gauge_model_best_CO2_non_CO2_per_gauge <- best_CO2_non_CO2_per_gauge |> 
  select(gauge, streamflow_model) |> 
  distinct()

streamflow_data_best_CO2_non_CO2 <- streamflow_data |> 
  semi_join(
    only_gauge_model_best_CO2_non_CO2_per_gauge,
    by = join_by(gauge, streamflow_model)
  ) 


#percentage_diff <- function(x, y) {
#  (abs(x - y) / ((x + y) / 2)) * 100
#}

### Difference in streamflow for the models
difference_streamflow_per_year_best_CO2_non_CO2 <- streamflow_data_best_CO2_non_CO2 |> 
  # Rename streamflow models to non-CO2 and CO2
  mutate(
    CO2_or_non_CO2 = if_else(str_detect(streamflow_model, "CO2"), "CO2_streamflow", "non_CO2_streamflow")
  ) |> 
  # Remove streamflow model to make CO2 and non-CO2 on same row
  select(!c(streamflow_model, objective_function, optimiser)) |> 
  pivot_wider(
    names_from = CO2_or_non_CO2,
    values_from = modelled_boxcox_streamflow
  ) |> 
  mutate(
    yearly_CO2_non_CO2_difference = CO2_streamflow - non_CO2_streamflow,
    percentage_yearly_CO2_non_CO2_difference = ((CO2_streamflow - non_CO2_streamflow) / CO2_streamflow) * 100#,
    #alternative_percentage = percentage_diff(CO2_streamflow, non_CO2_streamflow)
  )
  
### Summarise yearly differences
summary_streamflow_best_CO2_non_CO2 <- difference_streamflow_per_year_best_CO2_non_CO2 |> 
  summarise(
    mean_yearly_CO2_non_CO2_difference = mean(yearly_CO2_non_CO2_difference),
    sum_yearly_CO2_non_CO2_difference = sum(yearly_CO2_non_CO2_difference),
    mean_percent_yearly_CO2_non_CO2_difference = mean(percentage_yearly_CO2_non_CO2_difference),
    sum_CO2_streamflow = sum(CO2_streamflow),
    sum_non_CO2_streamflow = sum(non_CO2_streamflow),
    .by = gauge
  ) |> 
  # add evidence ratio
  left_join(
    evidence_ratio_calc,
    by = join_by(gauge)
  ) |> 
  # percentage diff
  
  arrange(sum_yearly_CO2_non_CO2_difference)
  
### Best CO2 vs. non-CO2 model comparison ######################################
summary_streamflow_best_CO2_non_CO2 |> 
  filter(evidence_ratio > 0) |> 
  #filter(abs(mean_percent_yearly_CO2_non_CO2_difference) < 1) |> 
  ggplot(aes(x = evidence_ratio, y = mean_percent_yearly_CO2_non_CO2_difference)) +
  geom_point() +
  scale_x_log10() +
  theme_bw() +
  labs(
    y = "mean(CO2_t - non-CO2_t / CO2_t)"
  )

# This graph tells us 
# x-axis is evidence ratio between CO2 and non-CO2 (log10 scaled)
# x-axis values < 0 removed for log10 scale. This removes catchments
# where the non-CO2 model is better
# the CO2 model on average produces -0.17 % less streamflow per year
# compared to the non-CO2 models
# There is no relationship between evidence ratio and difference in models
# We expected to see a larger difference in CO2 and non-CO2 models as 
# the evidence ratio grew.

# We really need to compare it to the observed...

summary_streamflow_best_CO2_non_CO2 |> 
  pull(mean_percent_yearly_CO2_non_CO2_difference) |> 
  quantile()


# Add slope and intercept to ToE histogram -------------------------------------
specific_CO2_comp_custom_bins_tim_of_emergence <- custom_bins_tim_of_emergence |> 
  mutate(
    intercept_or_slope = if_else(
      str_detect(streamflow_model, "intercept"), "Intercept", "Slope"
      )
    )

ToE_hist_v2 <- specific_CO2_comp_custom_bins_tim_of_emergence |> 
  summarise(
    n = n(),
    .by = c(custom_bins, intercept_or_slope)
  ) |> 
  mutate(
    custom_bins = factor(
      custom_bins, 
      levels = c(
        "Before 1959",
        paste0(seq(from = 1950, to = 2020, by = 10), "'s")
      )
    )
  ) |> 
  ggplot(aes(x = custom_bins, y = n, fill = intercept_or_slope)) + 
  geom_col(colour = "black") +
  labs(
    x = "Decade",
    y = "Count",
    title = bquote("What decade does"~CO[2]~"start impacting annual rainfall-partitioning?"),
    fill = NULL
  ) +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    legend.position = "inside",
    legend.position.inside = c(0.9, 0.9),
    legend.background = element_rect(fill = "white", colour = "black"),
    legend.text = element_text(size = 18)
  )

ToE_hist_v2

# Playing with geom_magnify ----------------------------------------------------
# Linux does not like this. Try on windows machine.

# Shape files in R are not nice:
# https://github.com/hughjonesd/ggmagnify/issues/30
# Shape files mean I cannot do two geom_magnify
# Is there alternative to shapefiles
# convert a shape file into something else? Convert shape file to polygon?
# Map an x-y map
# https://ggplot2-book.org/maps.html Polygon maps is an angle

polygon_test <- tribble(
  ~name, ~X, ~Y,
   "A",   0,   0,
   "A",   2,   0,
   "A",   1,   2,
   "B",   2,   0,
   "B",   4,   0,
   "B",   3,   2
) |> 
  ggplot(aes(x = X, y = Y, fill = name)) +
  geom_polygon()
polygon_test


## Map of Aus - magnify VIC
single_map_aus <- aus_map |> 
  ggplot() +
  geom_sf() +
  theme_bw() +
  geom_magnify(
    aes(from = state == "VIC"),
    to = c(130, 140, -45, -40), 
    shadow = FALSE, 
    aspect = "fixed", 
    expand = 0
    ) 

single_map_aus

magnify_vic_with_dots <- aus_map |> 
  ggplot() +
  geom_sf() +
  geom_point(
    data = a3_direction_binned_lat_lon_evidence_ratio,
    mapping = aes(x = lon, y = lat, colour = binned_evidence_ratio)
  ) +
  theme_bw() + 
  geom_magnify(
    aes(from = state == "VIC"),
    to = c(120, 140, -45, -35),
    shadow = FALSE
    ) 


