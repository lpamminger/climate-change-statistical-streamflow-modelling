## Fitting streamflow models to catchments
cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, dream, tictoc, furrr, parallel, truncnorm, sloop, tictoc)

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
#source("./Functions/my_cmaes.R")
source("./Functions/my_dream.R")
source("./Functions/objective_function_setup.R")
source("./Functions/result_set.R")


# Run DREAM --------------------------------------------------------------------
gauge <- "407214"

dream_example <- gauge |> 
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  ) |> 
  numerical_optimiser_setup_vary_inputs(
    streamflow_model = streamflow_model_precip_only,
    objective_function = constant_sd_objective_function,
    bounds_and_transform_method = make_default_bounds_and_transform_methods(),
    minimise_likelihood = FALSE
  ) |> 
  my_dream() |> 
  result_set() 

dream_parameters <- dream_example |> 
  parameters_summary()

dream_example |> 
  plot()

# N_PARAMETERS = 3 -> nseq = 50,000 
# Max parameters are 8. Extrapolate? if linear it would suggest round_any((50000/3) * 8, 1E4, ceiling) = 140,000
# I don't this this will work. I think the relationship between nseq and number of parameters is non-linear
