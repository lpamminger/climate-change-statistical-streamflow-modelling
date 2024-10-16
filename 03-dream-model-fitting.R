## Fitting streamflow models to catchments
cat("\014") # clear console

# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, dream, coda, tictoc, furrr, parallel, truncnorm, sloop, tictoc)
# install dream using: install.packages("dream", repos="http://R-Forge.R-project.org")
# install.packages("coda")


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
source("./Functions/my_dream.R")
source("./Functions/objective_function_setup.R")
source("./Functions/result_set.R")


# Psuedo-code for dream 
# Options
## 1. repeat 01-cmaes approach
## 2. Repeat the run dream pipe thing 

# pipe thing
# repeat_pipe <- function(gauge, streamflow_model, objective_function, data, start_stop) {
#}
# Need to separate into drought and non-drought
# separate again into batches based on ram requirements
# I could package everything (gauge, streamflow_model and objective function) into
# a single list then iterate over the list rather than using pmap?
# This would make batching/chunking easier.


# Steps.
# Separate drought catchments vs. all catchments. Drought catchment only get drought streamflow models
# Combine all gauge, streamflow model and objective function combinations in a list for iteration
# Split list up based on batches. Batches determined based on RAM requirements
# Run the repeat_pipe
# Save the results similar to 01-cmaes i.e., do not assign a variable,
# write_csv, append results and gc(). 
## Saving results - I would like the optimal parameters and sequences. Ignore the streamflow-time.
## sequences need to have |gauge|model|obj|parameter|parameter_value and rbind()



all_gauges <- unique(gauge_information$gauge)

drought_gauges <- gauge_information |> 
  filter(drought) |> 
  pull(gauge) |> 
  unique()


non_drought_streamflow_models <- get_non_drought_streamflow_models()

drought_streamflow_models <- get_drought_streamflow_models()

all_objective_functions <- get_all_objective_functions()


# I need to find all the combinations of all_gauges, non_drought_streamflow_models and all_objective_functions
non_drought_combinations <- tidyr::expand_grid(
  all_gauges,
  non_drought_streamflow_models,
  all_objective_functions
) |> 
  unclass()

# This is ready for pmap?
#pwalk(.l = non_drought_combinations,
 #    .f = )


repeat_dream <- function(gauge, streamflow_model, objective_function, data, start_stop) {
  
  dream_results <- gauge |> 
    catchment_data_blueprint(
      observed_data = data,
      start_stop_indexes = start_stop
    ) |> 
    numerical_optimiser_setup_vary_inputs(
      streamflow_model = streamflow_model,
      objective_function = objective_function,
      bounds_and_transform_method = make_default_bounds_and_transform_methods(),
      minimise_likelihood = FALSE
    ) |> 
    my_dream() |> 
    result_set() 
  
  # save parameters using parameter_summary()
  # save sequences dream_results$sequences
  
}

# If would like to save everything in the same chunk/batch
# this stops > 5000 files being made

# for the cmaes 
# function()
# future_map repeat for all gauges, models and objs in a given batch
# combine and save the results sets using two save functions.
# I already have this function run_and_save_chunks_my_cmaes_parallel()  !!!!!!!!!!!! DO THIS !!!!!!!!!!
# re-use this?

#furrr::future_pwalk()


# Run DREAM --------------------------------------------------------------------
tic()
gauge <- "407214"

dream_example <- gauge |> 
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  ) |> 
  numerical_optimiser_setup_vary_inputs(
    streamflow_model = streamflow_model_precip_only,
    objective_function = CO2_variable_objective_function,
    bounds_and_transform_method = make_default_bounds_and_transform_methods(),
    minimise_likelihood = FALSE
  ) |> 
  my_dream() |> 
  result_set() 

dream_parameters <- dream_example |> 
  parameters_summary()

dream_example |> 
  plot()

toc()

test_sequences <- dream_example$sequences

test_sequences |> 
  pivot_longer(
    cols = everything(),
    names_to = "parameter",
    values_to = "parameter_values"
  ) |> 
  ggplot(aes(x = parameter_values)) +
  geom_histogram(
    binwidth = binwidth_bins(30),
    fill = "grey",
    colour = "black"
    ) +
  theme_bw() +
  scale_y_sqrt() +
  labs(
    x = "Range of Parameter Values",
    y = "Frequency"
  ) +
  facet_wrap(~parameter, scales = "free")





# for gauge = 407214, model = streamflow_model_drought_separate_CO2_seasonal_ratio_auto
# and objective function CO2_variable_objective_function expect LL of 95
# 550,000 nseq and 425 sec (~7min) gives 105 LL
# I could bump it more. Make a stepper curve
# New equation used -> round_any(((x - 2) ^ 2.5) * 1E4, 1E4, ceiling)
# Untested for 8 parameter models. 

# N_PARAMETERS = 3 -> nseq = 50,000 
# Max parameters are 8. Extrapolate? if linear it would suggest round_any((50000/3), 1E4, ceiling) = 20000 per parameter
# I don't this this will work. I think the relationship between nseq and number of parameters is non-linear
# Want to get the likelihoods to the nearest whole number of DREAM and CMAES

# For gauge 407214 maximum parameters DREAM LL = 99.xxx, CMAES LL = 95.xxx
# need to increase from 20000 * PARAMETER_NUMBER to or make it a non-linear relationship 160,000 for 8 params its too little. Try 200,000? Not enough.
# Try 300,000? Not enough.
# Playing with non-linear relationship for nseq and n_parameters:
x <- seq(from = 3, to = 8, by = 1)
y <- ((x ^ 2) - (x + 1)) * 1E4
plot(x, y)

yy <- round_any(((x - 2) ^ 2.5) * 1E4, 1E4, ceiling)
yy
# non-linear relationship works with gauge 407214 and 3, 4 parameters
# test for 5, 6, 7, 8? Go straight to 8? Good idea to test all.


