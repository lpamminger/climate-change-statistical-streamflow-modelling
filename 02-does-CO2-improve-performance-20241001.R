# Map plots - work out if models with CO2 perform better than ones without

cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, sf)


# Import functions ------------------------------------------------------------- 
source("./Functions/utility.R")


# Import data ------------------------------------------------------------------
CMAES_results <- read_csv("Results/CMAES_results/CMAES_parameter_results_20240930.csv", show_col_types = FALSE)

gauge_information <- read_csv("./Data/Tidy/gauge_information_CAMELS_20240807.csv",
                              show_col_types = FALSE,
                              col_select = c(gauge, state, lat, lon)
                              ) # CAMELS is Australia wide.



aus_map <- read_sf(dsn = "./Data/Maps", 
                   layer = "AUS_2016_AUST")

vic_map <- read_sf(dsn = "./Data/Maps", 
                   layer = "vic_map_simple")


# Gauges that do weird things - remove for now ---------------------------------
removed_gauges <- c("G0050115", "G0060005", "A0030501", "226220", "226407", "226222", "308145", "225020A")

# For a given catchment work out if a model with CO2 outperforms a model without CO2
# remove unnecessary columns
# only have unique streamflow model, gauge combinations
# combine the streamflow model and objective function using unite
# regex to see if the united columns contain "CO2"
# if yes then contains_CO2 = TRUE otherwise contains_CO2 is false
# for each catchment (gauge) get the two lowest AIC values by contains_CO2

edited_results <- CMAES_results |> 
  filter(! gauge %in% removed_gauges) |> 
  select(!c(parameter, parameter_value, optimiser, loglikelihood, exit_message, near_bounds)) |> 
  distinct() |> 
  unite( 
    col = streamflow_model_objective_function,
    c(streamflow_model, objective_function),
    sep = "_",
    remove = TRUE,
    na.rm = FALSE
  ) |> 
  mutate(
    contains_CO2 = str_detect(streamflow_model_objective_function, "CO2"),
    contains_CO2 = if_else(contains_CO2, "CO2", "no_CO2"),
    .after = 2
  ) |> 
  slice_min(
    AIC,
    by = c(gauge, contains_CO2)
  ) |> 
  select(!streamflow_model_objective_function) |> 
  pivot_wider(
    names_from = contains_CO2,
    values_from = AIC
  )

# Evidence ratio calculation ---------------------------------------------------
# What to do (Taken from Burnham and Anderson 2002):
## - assume non-CO2 is the best. best_CO2_AIC - best_non_CO2_AIC = AIC_diff. (This assumption does not impact the results. See chapter 2.11.2)
## - Positive values mean best_non_CO2 is the better model. Negative values mean the best_CO2 model is better.
## - calculate relative likelihood exp(0.5 * AIC_diff) pg. 72
## - evidence ratio is a relative measure of how much better the model is over the other model

evidence_ratio_calc <- edited_results |> 
  mutate(
    AIC_difference = CO2 - no_CO2 # CO2 is smaller than non-CO2 then negative and CO2 is better
  ) |> 
  mutate(
    evidence_ratio = case_when(
      AIC_difference < 0 ~ -exp(0.5 * abs(AIC_difference)), # when CO2 model is better
      AIC_difference > 0 ~ exp(0.5 * abs(AIC_difference)) # when non-CO2 model is better
    )
  )





# Plot time --------------------------------------------------------------------
plot_ready_data <- evidence_ratio_calc |> 
  select(!c(CO2, no_CO2, AIC_difference)) |> 
  left_join(
    gauge_information,
    by = join_by(gauge)
  )


# I want the scale in these steps:
# From https://www.semanticscholar.org/paper/On-the-interpretation-of-likelihood-ratios-in-and-Martire-Kemp/3fc3690678409d8e8aa9352dce346565cf8fd0ea
# Weak or limited -> 1-10
# Moderate -> 10-100 
# Moderately strong -> 100-1,000
# Strong -> 1,000-10,000
# Very strong -> 10,000-1,000,000
# Extremely strong ->  >1,000,000



## Victorian plot ==============================================================
vic_evidence_ratio_map <- ggplot() +
  geom_sf(
    data = vic_map,
    colour = "black",
    fill = "grey50"
  ) +
  geom_point(
    data = plot_ready_data |> filter(state == "VIC"),
    aes(x = lon, y = lat, colour = evidence_ratio)
  ) +
  binned_scale( # binned scale code taken from: https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
    aesthetics = "colour",
    palette = function(x) c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#f7f7f7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac", "#053061"), # length should be length(breaks + limits) - 1
    breaks = c(-1E6, -1E4, -1E3, -1E2, -1E1, 1E1, 1E2, 1E3, 1E4, 1E6),
    limits = c(-1E13, 1E13),
    show.limits = FALSE,
    guide = "coloursteps"
  ) +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude",
    colour = "Evidence Ratio"
  ) +
  theme(
    legend.key.height = unit(25, "mm")
  )

ggsave(
  filename = paste0("vic_evidence_ratio_map_", get_date(), ".pdf"),
  plot = vic_evidence_ratio_map,
  device = "pdf",
  path = "./Graphs/CMAES_graphs",
  height = 210,
  width = 297,
  units = "mm"
)

# Australian plot ==============================================================
aus_evidence_ratio_map <- ggplot() +
  geom_sf(
    data = aus_map,
    colour = "black",
    fill = "grey50"
  ) +
  coord_sf(xlim = c(110, 155)) + 
  geom_point(
    data = plot_ready_data,
    aes(x = lon, y = lat, colour = evidence_ratio),
    size = 1
  ) +
  binned_scale( # binned scale code taken from: https://stackoverflow.com/questions/65947347/r-how-to-manually-set-binned-colour-scale-in-ggplot
    aesthetics = "colour",
    palette = function(x) c("#67001f", "#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#f7f7f7", "#d1e5f0", "#92c5de", "#4393c3", "#2166ac", "#053061"), # length should be length(breaks + limits) - 1
    breaks = c(-1E6, -1E4, -1E3, -1E2, -1E1, 1E1, 1E2, 1E3, 1E4, 1E6),
    limits = c(-1E13, 1E13),
    show.limits = FALSE,
    guide = "coloursteps"
  ) +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude",
    colour = "Evidence Ratio"
  ) +
  theme(
    legend.key.height = unit(25, "mm")
  )

ggsave(
  filename = paste0("aus_evidence_ratio_map_", get_date(), ".pdf"),
  plot = aus_evidence_ratio_map,
  device = "pdf",
  path = "./Graphs/CMAES_graphs",
  height = 210,
  width = 297,
  units = "mm"
)
