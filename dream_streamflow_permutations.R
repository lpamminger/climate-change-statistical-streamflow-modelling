# Objective:
# Using the DREAM results generate streamflow for all of the 
# parameter value combinations


# Method:
# - Get sequences file
# - Transform sequences for a given gauge into a matrix for the 
#   streamflow model (columns = different parameter sets, row = parameters a0, a1 etc.)
# - Using converged_stats get the streamflow model used.
# - Build catchment_data_set object
# - Feed catchment dataset object and parameters into the streamflow model
# - Save results (lots of space required)


# Proof of concept

# install packages -------------------------------------------------------------
pacman::p_load(tidyverse, dream, furrr, parallel, truncnorm, sloop)

# get functions ----------------------------------------------------------------
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

# import data ------------------------------------------------------------------
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

converged_DREAM_results <- read_csv(
  "Results/my_dream/converged_stats_20250429.csv",
  show_col_types = FALSE
  )




sequences_DREAM <- read_csv( # maybe join together in R? Need to remove duplicate gauges
  "Results/my_dream/big_chunk_1_sequences_20250509.csv", 
  show_col_types = FALSE
  )



# Get sequences for a specific gauge (maybe do fraction for now)
single_gauge_sequences_DREAM <- sequences_DREAM |> 
  filter(gauge == "410713")

# Transform into matrix
as_matrix_single_gauge_sequences_DREAM <- single_gauge_sequences_DREAM |> 
  select(chain, n, parameter, parameter_value) |> 
  unite(col = identifier, chain, n, sep = "_") |> 
  pivot_wider(
    names_from = identifier,
    values_from = parameter_value
  ) |> 
  select(!parameter) |>  # strip first column
  as.matrix() |> 
  unname()

# Get the corresponding streamflow model 
x <- converged_DREAM_results |> 
  filter(gauge == "410713") |> 
  pull(model) |> 
  unique() 

char_to_fun_model <- eval(parse(text = x))




# Build catchment_data_set
catchment_data <- catchment_data_blueprint(
  gauge_ID = "410713", 
  observed_data = data, 
  start_stop_indexes = start_stop_indexes
  )

# Put catchment_data and parameters into model
short_matrix_parameter_set <- as_matrix_single_gauge_sequences_DREAM#[,1:50]
#c(-5.42983753, 0.02297032, 0.24775972, -0.04284651, -6.16000385, 55.5668732, 1.84033216)
# it works with a vector
# catchment_data_directly_to_streamflow_model does not work with matrices
# make it work - I think the autocorrelation bit doesn't work well
# with a matrix of values --> use recursion to force stop-start data issue here

streamflow_variations <- char_to_fun_model(
  catchment_data = catchment_data,
  parameter_set = short_matrix_parameter_set
)

x <- streamflow_variations |> 
  select(!year:seasonal_ratio) |> 
  add_column(
    gauge = "410713",
    .before = 1
  ) # these cannot be rbinded because each gauge has a different number of 
# chains

# Matrices are very fast
# it looks like I wont have to save the tibbles
# it is a save angle if I want to have the sen.slope in another file
# I will need to repeat the script twice (big_chunk_1 and big_chunk_2)

# Saving will fill up my storage space on my PC
# Maybe do sen's slope uncertainty calculations here...
# this will involves calculating sen's slope for every permuation
# get the IQR
x <- sequences_DREAM |> 
  select(gauge, chain) |>
  distinct() |> 
  summarise(
    n = n(),
    .by = gauge
  )
