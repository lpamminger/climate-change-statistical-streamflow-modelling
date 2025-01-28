# Script produces:
# - a map of best CO2 vs best non-CO2 of Australia
# - streamflow timeseries comparing observed, best CO2 and best non-CO2


cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, ozmaps, sf, patchwork) # need to update this linux


# Import functions -------------------------------------------------------------
source("./Functions/utility.R")
source("./Functions/boxcox_transforms.R")


# Import data ------------------------------------------------------------------
data <- read_csv(
  "Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  show_col_types = FALSE
)

CMAES_results <- read_csv(
  "./Results/my_cmaes/CMAES_parameter_results_20250122.csv", 
  show_col_types = FALSE
) 

data <- readr::read_csv(
  "./Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  col_types = "icdddddddll", # ensuring each column is a the correct type
  show_col_types = FALSE
)

streamflow_results <- read_csv(
  "Results/my_cmaes/CMAES_streamflow_results_20250122.csv",
  show_col_types = FALSE,
  col_select = !optimiser
)

gauge_information <- read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE,
  col_select = c(gauge, bc_lambda, state, lat, lon)
) # CAMELS is Australia wide.

best_CO2_non_CO2_per_gauge <- read_csv(
  "./Results/my_cmaes/best_CO2_non_CO2_per_catchment_CMAES_20250124.csv",
  show_col_types = FALSE
) 
  




# 1. A map of best CO2 vs best non-CO2 of Australia ----------------------------

## Calculate evidence ratios ===================================================

# The flag column is produces errors here
# When pivot_wider is called flag modified and unmodified are on separate rows
# I want them on the same row. The easiest way is a left join of best_CO2...
flag_per_gauge <- best_CO2_non_CO2_per_gauge |> 
  select(gauge, flag) |> 
  distinct() |> 
  filter(flag == "modified")


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
  arrange(evidence_ratio) |> 
  # Add flag back in
  left_join(
    flag_per_gauge,
    by = join_by(gauge)
  ) |> 
  mutate(
    flag = if_else(is.na(flag), "unmodified", flag)
  )


### Get into plot ready form ###################################################
plot_ready_data <- evidence_ratio_calc |>
  select(!c(AIC_difference)) |>
  left_join(
    gauge_information,
    by = join_by(gauge)
  ) 


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

## By country ==================================================================

### How evidence ratio is binned ###############################################
# I want the scale in these steps:
# From https://www.semanticscholar.org/paper/On-the-interpretation-of-likelihood-ratios-in-and-Martire-Kemp/3fc3690678409d8e8aa9352dce346565cf8fd0ea
# Weak or limited -> 1-10
# Moderate -> 10-100
# Moderately strong -> 100-1,000
# Strong -> 1,000-10,000
# Very strong -> 10,000-1,000,000
# Extremely strong ->  >1,000,000

custom_palette <- function(x) {
  rev(c(
  "#67001f",
  "#b2182b",
  "#d6604d",
  "#f4a582",
  "#fddbc7",
  "#f7f7f7",
  "#d1e5f0",
  "#4393c3",
  "#2166ac",
  "#053061"
  ))
}

aus_evidence_ratio_map <- ggplot() +
  geom_sf(
    data = aus_map,
    colour = "black",
    fill = "grey50"
  ) +
  coord_sf(xlim = c(110, 155)) +
  geom_point(
    data = plot_ready_data,
    aes(x = lon, y = lat, colour = evidence_ratio, shape = flag),
    size = 0.75
  ) +
  binned_scale( # binned scale code taken from: https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
    aesthetics = "colour",
    palette = custom_palette, # length should be length(breaks + limits) - 1
    breaks = c(-1E4, -1E3, -1E2, -1E1, 1E1, 1E2, 1E3, 1E4, 1E6),
    limits = c(-1E6, 1E20),
    show.limits = TRUE,
    guide = "coloursteps"
  ) +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude",
    colour = "Evidence Ratio",
    shape = "Flag"
  ) +
  theme(
    legend.key.height = unit(4, "mm"),
    legend.key.width = unit(40, "mm"),
    legend.frame = element_rect(colour = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 11)
  )



ggsave(
  filename = paste0("aus_evidence_ratio_", get_date(), ".pdf"),
  plot = aus_evidence_ratio_map,
  device = "pdf",
  path = "./Graphs/CMAES_graphs",
  width = 420,
  height = 297,
  units = "mm"
)



# 2. Patchwork plots by state --------------------------------------------------

evidence_ratio_by_state <- function(state, evidence_ratio_data, map_data) {
  
  # use state as key to extract correct polygon from map_data and evidence_ratio_data
  state_map <- map_data |> 
    filter(state == {{ state }})
  
  state_evidence_ratio <- evidence_ratio_data |> 
    filter(state == {{ state }})
  
  ggplot() +
    geom_sf(
      data = state_map,
      colour = "black",
      fill = "grey50"
    ) +
    geom_point(
      data = state_evidence_ratio,
      aes(x = lon, y = lat, colour = evidence_ratio),
      show.legend = FALSE,
      size = 0.5
    ) +
    binned_scale( # binned scale code taken from: https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
      aesthetics = "colour",
      palette = function(x) c("#f7f7f7", "#fddbc7", "#f4a582",  "#d6604d", "#b2182b", "#67001f"), # length should be length(breaks + limits) - 1
      breaks = c(1E1, 1E2, 1E3, 1E4, 1E6),
      limits = c(-1E1, 1E27),
      show.limits = TRUE,
      guide = "coloursteps"
    ) +
    theme_bw() +
    labs(
      x = "Longitude",
      y = "Latitude",
      colour = "Evidence Ratio"
    ) +
    theme(
      legend.key.height = unit(25, "mm"),
      legend.frame = element_rect(colour = "black"),
      axis.title = element_blank(),
      axis.text = element_text(size = 6)
    )
}

by_state_plots <- map(
  .x = aus_map |> pull(state) |> unique(),
  .f = evidence_ratio_by_state,
  evidence_ratio_data = plot_ready_data |> select(!flag),
  map_data = aus_map
)

names(by_state_plots) <- aus_map |> pull(state) |> unique()

## Combine to make a paper ready plot ==========================================
# make it look nice - take some trial and error
# reduce(.x = by_state_plots, .f = `+`) easy add


top_nice_plot <- by_state_plots[["TAS"]] | by_state_plots[["QLD"]] | by_state_plots[["WA"]] | by_state_plots[["SA"]] | by_state_plots[["NT"]]
bottom_nice_plot <- by_state_plots[["VIC"]] | aus_evidence_ratio_map | by_state_plots[["NSW"]]
nice_plot <- top_nice_plot / bottom_nice_plot / guide_area() +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a")

nice_plot

ggsave(
  filename = paste0("nice_plot_", get_date(), ".pdf"),
  plot = nice_plot,
  device = "pdf",
  path = "./Graphs/CMAES_graphs",
  height = 210,
  width = 297,
  units = "mm"
)



# 3. Streamflow plots of best-CO2, best-non-CO2 and observed -------------------


## Filter based on best non-CO2 and CO2 models =================================
best_streamflow_results <- streamflow_results |>
  semi_join(
    best_CO2_non_CO2_per_gauge,
    by = join_by(gauge, streamflow_model)
  )


## Only include best streamflow that was calibrated on =========================
## join the included_in_calibration column
in_calibration <- data |> 
  select(year, gauge, included_in_calibration)

best_calibration_streamflow_results <- best_streamflow_results |> 
  left_join(
    in_calibration,
    by = join_by(year, gauge)
  ) |> 
  filter(included_in_calibration)



## Summarise results into a tidy format ========================================
tidy_boxcox_streamflow <- best_calibration_streamflow_results |>
  drop_na() |>  # only include if observed streamflow is present
  pivot_longer(
    cols = c(observed_boxcox_streamflow, modelled_boxcox_streamflow),
    names_to = "name",
    values_to = "boxcox_streamflow"
  ) |>
  mutate(
    name = if_else(name == "observed_boxcox_streamflow", "observed", "modelled")
  ) |> 
  mutate(
    contains_CO2 = str_detect(streamflow_model, "CO2")
  ) |> 
  mutate(
    legend = case_when(
      contains_CO2 & (name != "observed") ~ "modelled_CO2",
      !contains_CO2 & (name != "observed") ~ "modelled_non_CO2",
      .default = name
    )
  ) |> 
  select(
    !c(streamflow_model, objective_function, included_in_calibration, name, contains_CO2)
    )



## Convert from box-cox space to real space ====================================
### bc lambda found in gauge_information
# boxcox_inverse_transform()
tidy_streamflow <- tidy_boxcox_streamflow |>
  left_join(
    gauge_information,
    by = join_by(gauge)
  ) |>
  select(!c(state, lat, lon)) |>
  mutate(
    streamflow = boxcox_inverse_transform(yt = boxcox_streamflow, lambda = bc_lambda, lambda_2 = 1),
    .by = gauge
  )


## Plot results ================================================================
### having 534 graphs on a single page is too much = slows pc.
### Spread across multiple

# rep by gauge

# Split the tibble into X groups
# Gauges must not be across multiple groups
# Randomly assign 1,2 or 3 to each group then split? This works
# but is probably not the best way of doing it
ready_split_tidy_streamflow <- tidy_streamflow |> 
  mutate(
    split = sample(c(1, 2, 3), size = 1, replace = TRUE),
    .by = gauge,
    .before = 1
  ) 

split_tidy_streamflow <- ready_split_tidy_streamflow |> 
  split(f = ready_split_tidy_streamflow$split)


chunk_streamflow_timeseries_plot <- function(data) {
  
  data |>
    ggplot(aes(x = year, y = streamflow, colour = legend)) +
    geom_line(na.rm = TRUE, alpha = 0.5) +
    geom_point(na.rm = TRUE, size = 0.5, alpha = 0.5) +
    theme_bw() +
    scale_colour_brewer(palette = "Set1") +
    labs(
      x = "Year",
      y = "Streamflow (mm)"
    ) +
    facet_wrap(~gauge, scales = "free_y") +
    theme(legend.title = element_blank())
  
}


plot_streamflow_timeseries <- map(
  .x = split_tidy_streamflow,
  .f = chunk_streamflow_timeseries_plot
)



ggsave(
  filename = paste0("streamflow_timeseries_comparison_", get_date(), ".pdf"),
  plot = gridExtra::marrangeGrob(plot_streamflow_timeseries, nrow = 1, ncol = 1),
  device = "pdf",
  path = "./Graphs/CMAES_graphs",
  width = 1189,
  height = 841,
  units = "mm"
)

