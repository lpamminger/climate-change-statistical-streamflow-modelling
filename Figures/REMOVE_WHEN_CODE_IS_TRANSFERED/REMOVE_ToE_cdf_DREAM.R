cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, furrr, parallel)



## Import functions ===========================================================
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
source("./Functions/boxcox_transforms.R")


# SKIP #

## Just extract ToE from DREAM - makes data easier to work with ================
## Commment out when complete
### Load .csv
big_chunk_2 <- read_csv(
  "Results/my_dream/big_chunk_2_sequences_20250429.csv",
  show_col_types = FALSE
  )

### convert to .rda (requires significantly less space)
save(big_chunk_2, file = "Results/my_dream/big_chunk_2_sequences.rda") # Save big_chunk_2 during lunch - it takes some time


### Load file and extract only a5 parameters
# load(file = "Results/my_dream/big_chunk_2_sequences.rda")
# This is way better -> use .rda files

### Join and save only_a5_big_chunks into single .rda file
only_a5_big_chunk_2 <- big_chunk_2 |> 
  filter(parameter == "a5")

save(only_a5_big_chunk_2, file = "Results/my_dream/only_a5_big_chunk_2.rda")

load(file = "Results/my_dream/only_a5_big_chunk_1.rda")
load(file = "Results/my_dream/only_a5_big_chunk_2.rda")

only_a5_sequences <- rbind(only_a5_big_chunk_1, only_a5_big_chunk_2)

save(only_a5_sequences, file = "Results/my_dream/only_a5_sequences.rda")
# delete only_a5_big chunks 1 and 2


## Data ========================================================================

only_a5_sequences <- load(file = "Results/my_dream/only_a5_sequences.rda")

best_CO2_non_CO2_per_gauge <- read_csv(
  "./Results/my_cmaes/unmodified_best_CO2_non_CO2_per_catchment_CMAES_20250331.csv",
  show_col_types = FALSE
)

gauge_information <- read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE,
  col_select = c(gauge, bc_lambda, state, lat, lon)
) # CAMELS is Australia wide.



# Calculate evidence ratios ----------------------------------------------------
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


state_gauge <- gauge_information |>
  select(gauge, state)


binned_lat_lon_evidence_ratio <- evidence_ratio_calc |>
  mutate(
    binned_evidence_ratio = case_when(
      between(evidence_ratio, -1E1, 1E1) ~ "Weak",
      between(evidence_ratio, 1E1, 1E2) ~ "Moderate",
      between(evidence_ratio, 1E2, 1E3) ~ "Moderately Strong",
      between(evidence_ratio, 1E3, 1E4) ~ "Strong",
      between(evidence_ratio, 1E4, 1E6) ~ "Very Strong",
      between(evidence_ratio, 1E6, Inf) ~ "Extremely Strong",
      .default = NA
    )
  ) |>
  left_join(
    state_gauge,
    by = join_by(gauge)
  ) |>
  mutate(
    binned_evidence_ratio = factor(
      binned_evidence_ratio,
      levels = c("Weak", "Moderate", "Moderately Strong", "Strong", "Very Strong", "Extremely Strong")
    )
  ) |> 
  select(-c(CO2_model, non_CO2_model, AIC_difference))


# Filter DREAM results using binned_evidence_ratio and state -------------------
high_evi_gauges <- binned_lat_lon_evidence_ratio |> 
  filter(binned_evidence_ratio %in% c("Moderately Strong", "Strong", "Very Strong", "Extremely Strong")) |> 
  pull(gauge)

a5_high_evi_ratio_sequences <- only_a5_sequences |> 
  filter(gauge %in% high_evi_gauges) |> 
  select(gauge, parameter_value)


## Convert a5 parameter into a year ############################################
CO2_data <- readr::read_csv(
  "./Data/Raw/20241125_Mauna_Loa_CO2_data.csv",
  skip = 43,
  col_select = !unc,
  show_col_types = FALSE
) |>
  mutate(
    CO2_280 = mean - 280
  )

CO2 <- CO2_data |> pull(CO2_280)
year <- CO2_data |> pull(year)


### If CO2 for a given year minus the a5 parameter is negative then, the a5
### parameter has not kicked in. Take the first positive value and get the
### year using the CO2 index
single_a5_to_time_of_emergence <- function(a5, CO2, year) {
  adjusted_CO2 <- if_else(CO2 - a5 < 0, 0, CO2 - a5)
  
  year_where_CO2_impacts_flow <- year[adjusted_CO2 != 0][1]
  
  return(year_where_CO2_impacts_flow)
}

a5_to_time_of_emergence <- function(CO2, year) {
  # Save CO2 and year vectors to the function itself
  # They do not change
  
  force(CO2)
  force(year)
  stopifnot(is.numeric(CO2))
  stopifnot(is.numeric(year))
  
  function(a5) {
    future_map_dbl(
      .x = a5,
      .f = single_a5_to_time_of_emergence,
      CO2 = CO2,
      year = year,
      .progress = TRUE
    )
  }
}

# Function factory
adjusted_a5_to_ToE <- a5_to_time_of_emergence(CO2 = CO2, year = year)

### Apply to a5_vic_high_evi_ratio_sequences ###################################
plan(multisession, workers = length(availableWorkers())) # set once for furrr
ToE_high_evi_ratio_sequences <- a5_high_evi_ratio_sequences |> 
  mutate(
    ToE = adjusted_a5_to_ToE(parameter_value)
  )


# Generate ecdf for each site --------------------------------------------------
# for each gauge run ecdf function on ToE
# this will return a function for each gauge
# put the years (1959 to 2021) into each function
# join ecdf values into tibble 

extract_gauge_make_ecdf <- function(gauge, data) {
  data |> 
    filter(gauge == {{ gauge }}) |> 
    pull(ToE) |> 
    ecdf()
}

gauge_ecdf <- map(
  .x = high_evi_gauges,
  .f = extract_gauge_make_ecdf,
  data = ToE_high_evi_ratio_sequences
)

# You cannot map over function - make a function that allows this
map_ecdf <- function(specific_ecdf, x) {
  specific_ecdf(x)
}

values_gauge_ecdf <- map(
  .x = gauge_ecdf,
  .f = map_ecdf,
  x = seq(from = 1959, to = 2021, by = 1) # force the year to be the same for everything
)


ecdf_ToE_value <- do.call("cbind", values_gauge_ecdf) |> 
  `colnames<-`(high_evi_gauges) |> 
  as_tibble() |> 
  add_column(
    ToE = seq(from = 1959, to = 2021, by = 1),
    .before = 1
  )

# Get median and percentiles from ecdf_ToE_value
ecdf_ToE_summary <- ecdf_ToE_value |> 
  pivot_longer(
    cols = !ToE,
    names_to = "gauge",
    values_to = "ecdf_value"
  ) |> 
  left_join(
    state_gauge,
    by = join_by(gauge)
  ) |> 
  summarise(
    median_ToE_ecdf = median(ecdf_value),
    upper_bound_ToE_ecdf = quantile(ecdf_value, probs = 0.90), 
    lower_bound_ToE_ecdf = quantile(ecdf_value, probs = 0.10),
    .by = c(ToE, state)
  ) |> 
  # remove ACT and NT because there is only a single gauge each
  filter(!state %in% c("ACT", "NT"))


## Count gauges with moderately strong or greater ==============================
gauges_per_state <- high_evi_gauges |> 
  as_tibble() |> 
  rename(gauge = value) |> 
  left_join(
    state_gauge,
    by = join_by(gauge)
  ) |> 
  count(state)|> 
  # remove ACT and NT because there is only a single gauge each
  filter(!state %in% c("ACT", "NT")) |> 
  # add x and y position for labeller
  add_column(
    x_pos = 2010,
    y_pos = 0.05
  ) |> 
  # add nice label for plotting 
  mutate(
    n_label = paste("n = ", n)
  )

# Plotting ---------------------------------------------------------------------
ecdf_ToE_graph <- ecdf_ToE_summary |> 
  ggplot(aes(x = ToE, y = median_ToE_ecdf)) +
  geom_step() +
  geom_ribbon(
    aes(ymin = lower_bound_ToE_ecdf, ymax = upper_bound_ToE_ecdf),
    alpha = 0.2
    ) +
  geom_label(
    aes(x = x_pos, y = y_pos, label = n_label),
    data = gauges_per_state,
    inherit.aes = FALSE
  ) +
  labs(
    x = "Time of Emergence",
    y = "ecdf",
    title = "Time of Emergence Emperical Cumulative Distributions",
    subtitle = "Evidence ratio moderately strong or greater (>100). Envelope p10 and p90."
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  facet_wrap(~state)

ggsave(
  filename = "ecdf_ToE.pdf",
  plot = ecdf_ToE_graph,
  path = "Graphs/Figures",
  width = 297,
  height = 210,
  units = "mm"
)




