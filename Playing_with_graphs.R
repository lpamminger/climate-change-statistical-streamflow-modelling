## Fitting streamflow models to catchments
cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, truncnorm, sloop, ozmaps, sf, patchwork)



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

CMAES_results <- read_csv("./Results/my_cmaes/CMAES_parameter_results_20250122.csv",
                          show_col_types = FALSE
) 

best_CO2_non_CO2_per_gauge <- read_csv(
  "./Results/my_cmaes/unmodified_best_CO2_non_CO2_per_catchment_CMAES_20250128.csv",
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

streamflow_data <- read_csv(
  "Results/my_cmaes/CMAES_streamflow_results_20250122.csv",
  show_col_types = FALSE
)

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
  filter(streamflow_model == "streamflow_model_separate_shifted_CO2") |> 
  pull(parameter_value)

# streamflow model is not clever enough to pull the necessary bits
# from catchment data. Use must give it the tibble.
# Names are confusing. 
# Maybe have a function that detects if catchment data is catchment_data object
# if so then pull either the start stop or entire tibble and remove NAs

a3_off_parameters <- parameters
a3_off_parameters[3] <- 0

a3_off <- streamflow_model_separate_shifted_CO2(
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
ggsave(
  filename = "ToE_graphs.pdf",
  plot = ToE_hist,
  device = "pdf",
  path = "./Graphs",
  width = 297,
  height = 210,
  units = "mm"
)


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


ggsave(
  filename = "playing_with_graphs.pdf",
  plot = gridExtra::marrangeGrob(
    graphs, 
    nrow = 1, 
    ncol = 1,
    top = NULL
  ),
  # remove page numbers
  device = "pdf",
  path = "./Graphs",
  width = 297,
  height = 210,
  units = "mm"
)



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


## Map 
ggplot() +
  geom_sf(
    data = aus_map,
    colour = "black",
    fill = "grey50"
  ) +
  geom_point(
    data = binned_ToE_map_info,
    aes(x = lon, y = lat, colour = custom_bins), #shape = flag),
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
      aes(x = lon, y = lat, colour = custom_bins), 
      size = 0.75,
      show.legend = TRUE # Required to show ever factor in legend
    ) +
    scale_colour_brewer(
      palette = "RdYlBu",
      drop = FALSE # Required to show ever factor in legend
    ) +
    labs(
      colour = "Time of Emergence",
      x = NULL,#"Longitude",
      y = NULL#"Latitude"
    ) +
    theme_bw() +
    theme(
      legend.key = element_rect(fill = "grey50"),
      legend.background = element_rect(fill = "grey50", colour = "black", ),
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
fr_nice_plot <- fr_top / fr_bottom / guide_area() + plot_layout(guides = "collect")

fr_nice_plot



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
      legend.background = element_rect(fill = "grey50", colour = "black", ),
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

ggsave(
  filename = "streamflow_model_nice_plot.pdf",
  plot = model_nice_plot,
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)


