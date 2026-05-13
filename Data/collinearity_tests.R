# Examine the collinearity between rainfall and CO2

# Import libraries required ----------------------------------------------------
pacman::p_load(tidyverse, olsrr)


# Import functions -------------------------------------------------------------
source("./Functions/utility.R")

# Import data ------------------------------------------------------------------
data <- readr::read_csv( # data will be in the package
  "./Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  show_col_types = FALSE
) |>
  mutate(
    year = as.integer(year)
  ) |>
  # required for log-sinh. Log-sinh current formulation has asymptote of zero.
  # This means zero flows of ephemeral catchments cannot be transformed
  # add a really small value
  mutate(q_mm = q_mm + .Machine$double.eps^0.5) 


gauge_information <- readr::read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE
)

lat_lon_gauge <- gauge_information |>
  select(gauge, state, lat, lon)


best_CO2_non_CO2_per_gauge <- read_csv(
  "./Modelling/Results/CMAES/best_CO2_non_CO2_per_catchment_CMAES.csv",
  show_col_types = FALSE
)


# Get gauges with moderately_strong (100) or greater ---------------------------
evidence_ratio <- best_CO2_non_CO2_per_gauge |>
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
  select(gauge, evidence_ratio)


## add to data =================================================================
data <-  data |> 
  left_join(
    evidence_ratio,
    by = join_by(gauge)
  )



# Collinearity of rainfall and CO2 ---------------------------------------------

## Calculate correlation coefficient and R2 ====================================
correlation_of_p_mm_and_CO2 <- data |> 
  summarise(
    cor = cor(x = p_mm, y = CO2),
    max_p = max(p_mm), # for plotting
    .by = gauge
  ) |> 
  # add R^2 
  mutate(
    R2 = cor^2
  ) |> 
  # add x_axis for plotting
  mutate(
    x_lab = 55 # CO2 fixed x-axis can be constant
  ) |> 
  # round for plotting
  mutate(
    R2_round = paste0("R^2==", signif(R2, digits = 2))
  ) |> 
  arrange(desc(R2))


## Summary statistic of correlation and R2 =====================================
correlation_of_p_mm_and_CO2 |>
  left_join(
    evidence_ratio,
    by = join_by(gauge)
  ) |> 
  mutate(
    moderately_strong_or_greater = evidence_ratio > 100
  ) |> 
  summarise(
    mean_R2 = mean(R2),
    sd_R2 = sd(R2)#,
   # .by = moderately_strong_or_greater
  )

## Plot CO2 against p_mm =======================================================
CO2_precip_correlation <- data |>
  # add group to distinguish plots based on evidence ratio
  mutate(
    moderately_strong_or_greater = evidence_ratio > 100
  ) |> 
  arrange(desc(moderately_strong_or_greater)) 

ordering_facets <- CO2_precip_correlation |> pull(gauge) |> unique()

CO2_precip_correlation_plot <- CO2_precip_correlation |> 
  ggplot(aes(x = CO2, y = p_mm, shape = moderately_strong_or_greater)) +
  geom_smooth(
    method = lm,
    formula = y ~ x,
    colour = "red",
    fill = "red",
    alpha = 0.1,
    linewidth = 0.8
  ) +
  geom_point(size = 1) +
  geom_label(
    data = correlation_of_p_mm_and_CO2,
    mapping = aes(x = x_lab, y = max_p, label = R2_round),
    inherit.aes = FALSE,
    parse = TRUE,
    size = 2
  ) +
  labs(
    x = "CO2 - 280 (ppm)",
    y = "Annual Precipitation (mm)"
  ) +
  theme_bw() +
  facet_wrap(~factor(gauge, levels = ordering_facets), scales = "free_y")

ggsave(
  filename = "./Figures/Other/CO2_precip_correlation_plot.pdf",
  plot = CO2_precip_correlation_plot,
  device = "pdf",
  width = 1189,
  height = 841,
  units = "mm"
)




# Calculate the Variance Inflation Factor --------------------------------------


# identify independent variables - look at models - not sure how to deal with autocorrelation and drought (add dummy?)
# Do I only do it for catchments with evidence ratio > 100 i.e., just check the catchments where CO2 is the best?
# I think we do it for everything then filter it down later - 

# pseudo code
# function(gauge, data)
# see if there is a drought term - if yes include it in the regression
# otherwise leave it
# return tibble - append gauge to leftmost column


VIF_assessment <- function(gauge, data) {
  
  # filter by gauge
  gauge_data <- data |> 
    filter(gauge == {{ gauge }})
  
  # check if any drought is TRUE
  drought_check <- gauge_data |> pull(drought) |> any()
  
  if(drought_check) {
    
    # need to convert drought to factor for lm
    drought_gauge_data <- data |> 
      mutate(drought = as.factor(drought))
      
    lm(q_mm ~ drought + p_mm + standardised_warm_season_to_annual_rainfall_ratio + CO2, data = drought_gauge_data) |> 
      ols_vif_tol() |> 
      add_column(
        gauge = {{ gauge }},
        .before = 1
      )
    
  } else {
    
    lm(q_mm ~ p_mm + standardised_warm_season_to_annual_rainfall_ratio + CO2, data = gauge_data) |> 
      ols_vif_tol() |> 
      add_column(
        gauge = {{ gauge }},
        .before = 1
      )
    
  }
}


VIF_results <- map(
  .x = data |> pull(gauge) |> unique(),
  .f = VIF_assessment,
  data = data
) |> 
  list_rbind()


VIF_plot <- VIF_results |> 
  ggplot(aes(x = Variables, y = VIF)) +
  geom_boxplot() +
  theme_bw()

ggsave(
  filename = "./Figures/Other/VIF_independent_variables_boxplot.pdf",
  plot = VIF_plot,
  device = "pdf",
  width = 200,
  height = 200,
  units = "mm"
)
