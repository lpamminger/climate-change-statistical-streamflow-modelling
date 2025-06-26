# Making DAGs

# Install packages
library(dagitty)
library(tidyverse)
library(ggdag)
library(broom)
library(patchwork)


# Get streamflow models
source("./Functions/streamflow_models.R")
get_drought_streamflow_models()
# Construct DAG ----------------------------------------------------------------

# Key variables used:
data <- readr::read_csv(
  "./Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  col_types = "icdddddddll", # ensuring each column is a the correct type
  show_col_types = FALSE
)

# https://r-causal.github.io/ggdag/reference/repel.html ?

colnames(data)
# - p_mm
# - standardised_warm_season_to_annual_rainfall_ratio
# - drought
# - CO2

# Other variables
# - autocorrelation Qbct-1 - average(Qbct-1)

# Dependent
# - both drought and warm season to annual rainfall both depend on precip
# - autocorrelation depends on drought, warm_to_..., CO2 (everything expect auto)

# The purpose of the DAGs is not to show the full equations. Rather the 
# relationships between them



## streamflow_model_precip_only_DAG ============================================
streamflow_model_precip_only_dagify <- dagify(
  boxcox_streamflow ~ CO2 + precipitation,
  precipitation ~ CO2,
  exposure = "CO2",
  outcome = "boxcox_streamflow",
  labels = c(boxcox_streamflow = "Box-Cox Streamflow", CO2 = "CO2", precipitation = "Precipitation"),
  coords = list(x = c(boxcox_streamflow = 2, CO2 = 0, precipitation = 1), y = c(boxcox_streamflow = 0, CO2 = 0, precipitation = 2))
)

streamflow_model_precip_only_dag <- ggdag_status(
  streamflow_model_precip_only_dagify, 
  use_labels = "label", 
  nudge_y = c(0.9, -0.9, -0.9),
  text = FALSE
  ) +
  theme_dag() +
  guides(fill = "none", colour = "none") # Disable the legend
streamflow_model_precip_only_dag

## streamflow_model_drought_precip_auto_DAG ====================================
# How do we do autocorrelation in a DAG?

## streamflow_model_drought_precip_seasonal_ratio_DAG ==========================
## streamflow_model_drought_precip_seasonal_ratio_auto_DAG =====================
## streamflow_model_drought_intercept_shifted_CO2_DAG ==========================
## streamflow_model_drought_intercept_shifted_CO2_auto_DAG =====================
## streamflow_model_drought_intercept_shifted_CO2_seasonal_ratio_DAG ===========
## streamflow_model_drought_intercept_shifted_CO2_seasonal_ratio_auto_DAG ======
## streamflow_model_drought_slope_shifted_CO2_DAG ==============================
## streamflow_model_drought_slope_shifted_CO2_auto_DAG =========================
## streamflow_model_drought_slope_shifted_CO2_seasonal_ratio_DAG ===============
## streamflow_model_drought_slope_shifted_CO2_seasonal_ratio_auto_DAG ==========


# Save DAGs --------------------------------------------------------------------
ggsave(
  # list of plots - grob something
  device = "pdf",
  width = 297,
  height = 210,
  units = "mm"
)

