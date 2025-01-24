# Script produces:
# - a map of best CO2 vs best non-CO2 of Australia
# - streamflow timeseries comparing observed, best CO2 and best non-CO2


# TODO:
# 1. Get the fancy patchwork plots working
# 2. Get the streamflow timeseries of best plots working


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
) # CAMELS is Australia wide.


# 1. A map of best CO2 vs best non-CO2 of Australia ----------------------------

## Calculate evidence ratios ===================================================
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
    aes(x = lon, y = lat, colour = evidence_ratio),
    size = 0.5
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
    colour = "Evidence Ratio"
  ) +
  theme(
    legend.key.height = unit(4, "mm"),
    legend.key.width = unit(40, "mm"),
    legend.frame = element_rect(colour = "black"),
    legend.position = "bottom"
  )


ggsave(
  filename = "test.pdf",
  plot = aus_evidence_ratio_map,
  device = "pdf",
  width = 420,
  height = 297,
  units = "mm"
)

