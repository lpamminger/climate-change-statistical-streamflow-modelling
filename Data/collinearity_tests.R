# Examine the collinearity between rainfall and CO2

# Import libraries required ----------------------------------------------------
pacman::p_load(tidyverse, olsrr, patchwork)


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


# add mortons ET to data -------------------------------------------------------
mortons_wet_daily <- read_csv(
  "Data/Raw/et_morton_wet_SILO.csv",
  show_col_types = FALSE
)


mortons_wet_annual <- mortons_wet_daily |>
  pivot_longer(
    cols = !c(year, month, day),
    names_to = "gauge",
    values_to = "APET"
  ) |>
  summarise(
    APET = sum(APET),
    n = n(),
    .by = c(year, gauge)
  ) |>
  # filter out incomplete years
  filter(n >= 365)


## join to data
data <- data |> 
  left_join(
    mortons_wet_annual,
    by = join_by(gauge, year)
  )


# add temperature (tmax) to the data -------------------------------------------
# use maximum daily temperature - AGCD temperature (there is also SILO) I am not sure of the difference
daily_tmax <- read_csv(
  "Data/Raw/tmax_AGCD.csv",
  show_col_types = FALSE
)


annual_tmax <- daily_tmax |>
  pivot_longer(
    cols = !c(year, month, day),
    names_to = "gauge",
    values_to = "tmax"
  ) |>
  summarise(
    tmax = mean(tmax),
    n = n(),
    .by = c(year, gauge)
  ) 
  # there are no incomplete years






## join to data
## Annual tmax is the maximum observed temperature in a given year per gauge
data <- data |> 
  left_join(
    annual_tmax,
    by = join_by(gauge, year)
  ) |> 
  select(!c(n.x, n.y))


CO2_vs_temperature_plot <- data |> 
  ggplot(aes(x = CO2, y = tmax)) +
  geom_point() +
  geom_smooth(
    formula = y ~ x,
    method = lm,
    colour = "red",
    se = FALSE
  ) +
  theme_bw() +
  facet_wrap(~gauge, scales = "free_y") 


ggsave(
  filename = "./Figures/Other/CO2_mean_max_daily_temp_correlation_plot.pdf",
  plot = CO2_vs_temperature_plot,
  device = "pdf",
  width = 1189,
  height = 841,
  units = "mm"
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



# Correlation between all independent variables --------------------------------
# Correlation matrix per catchment 

make_correlation_matrix_per_gauge <- function(gauge, data) {
  
  # if there is not drought it returns a warning
  
  data |> 
    filter(gauge == {{ gauge }}) |>
    # avoid NA's being an issue
    filter(included_in_calibration) |> 
    # only select independent variables
    select(
      p_mm, 
      standardised_warm_season_to_annual_rainfall_ratio,
      CO2,
      drought,
      APET,
      tmax
    ) |>
    as.matrix() |> 
    cor()
  
}


correlation_matrix_per_gauge <- map(
  .x = data |> pull(gauge) |> unique(),
  .f = make_correlation_matrix_per_gauge,
  data = data
) |> 
  # convert to 3D array
  simplify2array()
# warnings occur because there are catchment without drought be corr


## Find the mean, median and standard deviation for across all matrices --------
count_greater_than_abs <- function(x, value) {
  boolean_vector <- abs(x) > value
  return(sum(boolean_vector, na.rm = TRUE))
}

correlation_threshold <- 0.5

apply(
  correlation_matrix_per_gauge,
  MARGIN = c(1, 2),
  FUN = count_greater_than_abs, # change values to test
  simplify = TRUE,
  # correlations are subject - I selected 0.6
  value = correlation_threshold
  # required for the drought/non-drought catchments
  #na.rm = TRUE
  ) |> 
  round(digits = 3)

# address problematic values for 
# - p_mm, seasonality, CO2, drought
# - p_mm and drought = 2 catchment > 0.5
# - p_mm and seasonality = 8 catchment > 0.5
# - drought and CO2 = 33 catchments > 0.5
 

# do they these catchments have a high evidence ratio?
cor_p_mm_seasonality <- data |> 
  filter(included_in_calibration) |> 
  summarise(
    cor = cor(p_mm, standardised_warm_season_to_annual_rainfall_ratio),
    .by = gauge
  ) |> 
  filter(abs(cor) > correlation_threshold)


gauges_with_drought <- gauge_information |> 
  filter(drought) |> 
  pull(gauge)


cor_drought_CO2 <- data |> 
  filter(included_in_calibration) |> 
  filter(gauge %in% gauges_with_drought) |> 
  mutate(
    drought = as.numeric(drought)
  ) |> 
  summarise(
    cor = cor(CO2, drought),
    .by = gauge
  ) |> 
  filter(abs(cor) > correlation_threshold)


cor_drought_p_mm <- data |> 
  filter(included_in_calibration) |> 
  filter(gauge %in% gauges_with_drought) |> 
  summarise(
    cor = cor(drought, p_mm),
    .by = gauge
  ) |> 
  filter(abs(cor) > correlation_threshold)




# - identify catchments cor_CO2_APET |> pull(gauge)
high_corr_gauges <- rbind(cor_p_mm_seasonality, cor_drought_CO2, cor_drought_p_mm) |> 
  pull(gauge) |> 
  unique()

# - do they have a high evidence ratio? 
evi_ratio_high_corr_gauges <- data |> 
  filter(gauge %in% high_corr_gauges) |> 
  select(gauge, evidence_ratio) |> 
  distinct() |> 
  arrange(desc(evidence_ratio))



# - if so investigate streamflow time 
# I am not convinced with the drought stuff
data |> 
  filter(gauge %in% high_corr_gauges) |> 
  ggplot(aes(y = CO2, x = as.numeric(drought))) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = lm) +
  theme_bw() +
  facet_wrap(~gauge, scales = "free_y")




median_correlation_matrix <- apply(
  correlation_matrix_per_gauge,
  MARGIN = c(1, 2),
  FUN = median,
  simplify = TRUE,
  na.rm = TRUE
) |> 
  round(digits = 3)


sd_correlation_matrix <- apply(
  correlation_matrix_per_gauge,
  MARGIN = c(1, 2),
  FUN = sd,
  simplify = TRUE,
  na.rm = TRUE
) |> 
  round(digits = 3)



# Calculate the Variance Inflation Factor --------------------------------------


# identify independent variables - look at models - not sure how to deal with autocorrelation and drought (add dummy?)
# Do I only do it for catchments with evidence ratio > 100 i.e., just check the catchments where CO2 is the best?
# I think we do it for everything then filter it down later - 

# pseudo code
# function(gauge, data)
# see if there is a drought term - if yes include it in the regression
# otherwise leave it
# return tibble - append gauge to leftmost column


VIF_assessment <- function(gauge, independent_variables, data) {
  
  # filter by gauge
  gauge_data <- data |> 
    filter(gauge == {{ gauge }}) |> 
    drop_na()
  
  # check if any drought is TRUE
  drought_check <- gauge_data |> pull(drought) |> any()
  
  
  if(drought_check) {
    
    # need to convert drought to factor for lm
    drought_gauge_data <- gauge_data |> 
      mutate(drought = as.factor(drought))
    
    # add drought to independent variables
    independent_variables <- c("drought", independent_variables)
    
    # convert character vector into formula
    formula <- reformulate(termlabels = independent_variables, response = "q_mm")
      
    lm(formula = formula, data = drought_gauge_data) |> 
      ols_vif_tol() |> 
      add_column(
        gauge = {{ gauge }},
        .before = 1
      )
    
  } else {
    
    # convert character vector into formula
    formula <- reformulate(termlabels = independent_variables, response = "q_mm")
    
    lm(formula = formula, data = gauge_data) |> 
      ols_vif_tol() |> 
      add_column(
        gauge = {{ gauge }},
        .before = 1
      )
    
  }
}


## Repeat VIF assessment for different combinations of variables ===============

### VIF of independent variables used in analysis ##############################
VIF_results_1 <- map(
  .x = data |> pull(gauge) |> unique(),
  .f = VIF_assessment,
  independent_variables = c("p_mm", "standardised_warm_season_to_annual_rainfall_ratio", "CO2"), #exclude drought
  data = data
) |> 
  list_rbind() |> 
  add_column(test = "default_variables")

# only 1 gauge over 4
VIF_results_1 |> 
  filter(VIF > 4)




### VIF of independent variables used in analysis and APET #####################
VIF_results_2 <- map(
  .x = data |> pull(gauge) |> unique(),
  .f = VIF_assessment,
  independent_variables = c("p_mm", "standardised_warm_season_to_annual_rainfall_ratio", "CO2", "APET"), #exclude drought
  data = data
) |> 
  list_rbind() |> 
  add_column(test = "with_APET")

# only 2 gauge over 4
VIF_results_2 |> 
  filter(VIF > 4)




### VIF of independent variables used in analysis and temp #####################
VIF_results_3 <- map(
  .x = data |> pull(gauge) |> unique(),
  .f = VIF_assessment,
  independent_variables = c("p_mm", "standardised_warm_season_to_annual_rainfall_ratio", "CO2", "tmax"), #exclude drought
  data = data
) |> 
  list_rbind() |> 
  add_column(test = "with_tmax")

# bigggg issues
VIF_results_3 |> 
  filter(VIF > 4)





### VIF of independent variables used in analysis, APET and temp ###############
VIF_results_4 <- map(
  .x = data |> pull(gauge) |> unique(),
  .f = VIF_assessment,
  independent_variables = c("p_mm", "standardised_warm_season_to_annual_rainfall_ratio", "CO2", "tmax", "APET"), #exclude drought
  data = data
) |> 
  list_rbind() |> 
  add_column(test = "with_APET_and_tmax")

# only 1 gauge over 4
VIF_results_4 |> 
  filter(VIF > 4)




# make a 4 x 4 VIF
VIF_results <- rbind(VIF_results_1, VIF_results_2, VIF_results_3, VIF_results_4)

# add abc to facet
# strip facet labels

# tibble with x, y, test and label_name
# will be the same between plots, y will be max(VIF)
abc_labels <- VIF_results |> 
  summarise(
    y_pos = max(VIF),
    .by = test
  ) |> 
  mutate(x_pos = "APET") |> 
  mutate(label_name = c("a", "b", "c", "d")) |> 
  mutate(
    test = factor(test, levels = c("default_variables", "with_APET", "with_tmax", "with_APET_and_tmax"))
  )


# need make this nice

VIF_plot <- VIF_results |> 
  mutate(
    test = factor(test, levels = c("default_variables", "with_APET", "with_tmax", "with_APET_and_tmax"))
  ) |> 
  ggplot(aes(x = Variables, y = VIF)) +
  geom_boxplot(
    staplewidth = 0.5
  ) +
  geom_hline(yintercept = 4, colour = "red", linetype = "dashed") +
  geom_text(
    aes(x = x_pos, y = y_pos, label = label_name),
    data = abc_labels,
    fontface = "bold",
    hjust = 2
    )+
  labs(
    y = "Variance Inflation Factor (VIF)",
    x = "Independent Variables"
  ) +
  scale_x_discrete(
    labels = c("APET", bquote(CO[2]), "Drought", "Precipitation", "Rainfall Seasonality", "Temperature")
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.text = element_text(size = 8)
  ) +
  facet_wrap(~test, scales = "free_y")




ggsave(
  filename = "./Figures/Other/VIF_independent_variables_boxplot.pdf",
  plot = VIF_plot,
  device = "pdf",
  width = 280,
  height = 200,
  units = "mm"
)
