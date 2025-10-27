# Objective: Perform Wilcox Rank Sum Test to see if the climate change
#            a3 parameters (a3_slope and a3_intercept) are statistically
#            different from zero.

# Method:
# 1. Import a3 parameters (do a3_intercept and a3_slope separately)
# 2. Filter out gauges with an evidence ratio weak or moderate (> 100)
# 3. Facet_wrap all histograms
# 4. Using facet_wrap histograms determine which gauges require transforming
#    (Wilcox Rank Sum Test requires the distributions to be symmetrical)
# 5. Use box-cox to transform values
# 6. Apply Wilcox Rank Sum Test



# Import packages --------------------------------------------------------------
pacman::p_load(tidyverse, arrow)


# Import functions -------------------------------------------------------------
source("./Functions/boxcox_logsinh_transforms.R")


# Determine which gauges have an evidence ratio > 100 --------------------------
best_CO2_non_CO2_per_gauge <- read_csv(
  "./Modelling/Results/CMAES/best_CO2_non_CO2_per_catchment_CMAES.csv",
  show_col_types = FALSE
)


## Calculate evidence ratio ====================================================
high_evi_gauges <- best_CO2_non_CO2_per_gauge |>
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



# Import DREAM distributions ---------------------------------------------------
DREAM_sequence_data <- open_dataset(
  source = "./Modelling/Results/DREAM/Sequences",
  format = "parquet"
) |> 
  filter(parameter %in% c("a3_intercept", "a3_slope")) |> 
  filter(gauge %in% high_evi_gauges) |> 
  collect()


# Plot distributions using facet_wrap ------------------------------------------


plot_histograms <- function(DREAM_sequence_data, value) {
  DREAM_sequence_data |> 
    ggplot(aes(x = {{ value }})) +
    geom_histogram() +
    theme_bw() +
    facet_wrap(~gauge, scales = "free")
  
}

untrans_a3_intercept <- DREAM_sequence_data |> 
  filter(parameter == "a3_intercept") |> 
  plot_histograms(value = parameter_value)

untrans_a3_slope <- DREAM_sequence_data |> 
  filter(parameter == "a3_slope") |> 
  plot_histograms(value = parameter_value)

ggsave(
  filename = "a3_intercept_dist.pdf",
  plot = untrans_a3_intercept,
  device = "pdf",
  path = "Figures/Other",
  width = 1189,
  height = 841,
  units = "mm"
)


ggsave(
  filename = "a3_slope_dist.pdf",
  plot = untrans_a3_slope,
  device = "pdf",
  path = "Figures/Other",
  width = 1189,
  height = 841,
  units = "mm"
)

# Identify which gauges require transform --------------------------------------
# Approximate boxcox transformation depending on the shape of the distribution
# Negative values require shift to positive values using lambda_2
# boxcox lambda range -2 to 2 (or -3 to 3)
#  using MASS to fit values does not produce the "Best symmetrical" - do it manually

gauges_a3_transform_required <- tribble(
  # cannot use lambda 0 due to test being stat test against zero (undefined)
  ~gauge,   ~lambda, 
  # intercept
  "204036", -3,
  "224213", -2,
  "226204", 2, 
  "401210", 3,
  "405274", 2,
  "407214", 3,
  "410156", 3,
  "411003", -1,
  "603007", -1,
  "701002", -1,
  # slope
  "210022", 2,
  "401217", 2,
  "405234", 2,
  "405237", 2,
  "406226", 2,
  "415223", 2,
  "610008", 3
)

adjusted_DREAM_sequence_data <- DREAM_sequence_data |> 
  left_join(
    gauges_a3_transform_required,
    by = join_by(gauge)
  ) |> 
  # Do everything by gauge
  group_by(gauge) |> 
  # Assign lambda
  mutate(
    lambda = if_else(is.na(lambda), 1, lambda) # 1 = no transform
  ) |> 
  # Assign lambda_2
  mutate(
    min_parameter_value = min(parameter_value),
    lambda_2 = if_else(
      min_parameter_value < 0, 
      abs(min_parameter_value) + .Machine$double.eps^0.5,
      0
      )
    ) |> 
  # Transform
  mutate(
    transformed_parameter_value = boxcox_transform(
      y = parameter_value,
      lambda = lambda[1],
      lambda_2 = lambda_2[1]
    )
  ) |> # Add zero transformed term
  mutate(
    transformed_zero = boxcox_transform(
      # Transformed zero produces infs - remove by adding using a 
      # very small number instead of zero
      y = .Machine$double.eps^0.5, 
      lambda = lambda[1],
      lambda_2 = lambda_2[1]
  )
)



trans_a3_intercept <- adjusted_DREAM_sequence_data |> 
  filter(parameter == "a3_intercept") |> 
  plot_histograms(value = transformed_parameter_value)

trans_a3_slope <- adjusted_DREAM_sequence_data |> 
  filter(parameter == "a3_slope") |> 
  plot_histograms(value = transformed_parameter_value)

ggsave(
  filename = "transformed_a3_intercept_dist.pdf",
  plot = trans_a3_intercept,
  device = "pdf",
  path = "Figures/Other",
  width = 1189,
  height = 841,
  units = "mm"
)


ggsave(
  filename = "transformed_a3_slope_dist.pdf",
  plot = trans_a3_slope,
  device = "pdf",
  path = "Figures/Other",
  width = 1189,
  height = 841,
  units = "mm"
)


# Apply wilcox.test ------------------------------------------------------------

## wrap wilcox.test so it outputs a single p-value =============================

wilcox_test <- function(x, mu, ...) {
  
  force(x)
  
  if(length(mu) != 1) {
    mu <- unique(mu)
  }
  
  wilcox.test(x = x, mu = mu, ...)$p.value
}



wilcox_test_results <- adjusted_DREAM_sequence_data |>
  ungroup() |> 
  # still grouped by gauge
  summarise(
    # two-sided test
    two_side_test = wilcox_test(
      x = transformed_parameter_value, 
      mu = transformed_zero,
      alternative = "two.sided",
      conf.level = 0.95
      ),
    # left sided test
    less_test = wilcox_test(
      x = transformed_parameter_value, 
      mu = transformed_zero,
      alternative = "less",
      conf.level = 0.95
    ),
    # right sided test
    greater_test = wilcox_test(
      x = transformed_parameter_value, 
      mu = transformed_zero,
      alternative = "greater",
      conf.level = 0.95
    ),
    .by = c(gauge, parameter)
  )
  

