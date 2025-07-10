# An example to fit models to annual streamflow data --------------------------- 



# Import libraries required ----------------------------------------------------
# If this is a package this step is not required
pacman::p_load(tidyverse, cmaesr, smoof, truncnorm, sloop)


# Get the offset working
# Instead of altering the transformation just add + 1 to observed flow
# Then -1 after calibration


# Import and prepare data-------------------------------------------------------
# If this is a package this step is not required
start_stop_indexes <- readr::read_csv(
  "./Data/Tidy/start_end_index.csv",
  show_col_types = FALSE
)

data <- readr::read_csv( # data will be in the package
  "./Data/Tidy/with_NA_yearly_data_CAMELS.csv",
  show_col_types = FALSE
) |>
  mutate(
    year = as.integer(year)
  )

gauge_information <- readr::read_csv(
  "./Data/Tidy/gauge_information_CAMELS.csv",
  show_col_types = FALSE
)


## Import functions ============================================================
#If this is a package this step is not required
source("./Functions/utility.R")
source("./Functions/boxcox_logsinh_transforms.R")
source("./Functions/streamflow_models.R")
source("./Functions/parameter_transformations.R")
source("./Functions/catchment_data_blueprint.R")
source("./Functions/cmaes_dream_summaries.R")
source("./Functions/objective_functions.R")
source("./Functions/numerical_optimiser_setup.R")
source("./Functions/generic_functions.R")
source("./Functions/CMAES.R")
source("./Functions/objective_function_setup.R")
source("./Functions/result_set.R")



# 1. Select gauge to test from data --------------------------------------------

gauge <- "603190"


# Test - if streamflow < 1 than round to either 1 or 0
#round_data <- data |> 
#  filter(gauge == {{ gauge }}) |> 
#  mutate(
#    q_mm = if_else(q_mm < 1, round(q_mm, digits = 0), q_mm)
#  )

# ideally catchment_data_blueprint should have other methods of data entry such as giving vectors individually
# See catchment_data_blue_print helper... must leave observed_data blank
example_catchment <- gauge |> 
  catchment_data_blueprint(
    observed_data = data,
    start_stop_indexes = start_stop_indexes
  ) 

# plot to see if there is anything strange happening in data
plot(example_catchment, type = "streamflow-time")


# 2. Prepare numerical_optimiser object ----------------------------------------

# Here you can adjust:
## - the streamflow_model used. 
###   See get_all_streamflow_models() for a list off all available models.
## - objective_function used. 
###   See get_all_objective_functions() for a list of all available objective functions.
## - streamflow_transform_method used. 
###   See get_streamflow_transform_method() for a list of all available methods
## - bounds_and_transform method used. 
###   Currently only changing the tibble located in make_default_bounds_and_transform_methods 
###   can be used to change bounds and transform methods
## - minimise_likelihood. 
###   A logical argument that can either minimise or maximise negative loglikelihood (dependent of numerical optimiser - CMAES or DREAM)
## - streamflow offset used
###   For the streamflow_transform_method. Values > 0 help avoid negative 
###   streamflow values in the rainfall-runoff relationship. Only impacts
###   streamflow in the transformed space.

streamflow_transform_method_offset <-  1

numerical_optimiser <- example_catchment |>
  numerical_optimiser_setup_vary_inputs(
    streamflow_model = streamflow_model_precip_seasonal_ratio_auto,
    objective_function = constant_sd_objective_function, 
    streamflow_transform_method = log_sinh_transform, 
    bounds_and_transform_method = make_default_bounds_and_transform_methods(example_catchment, streamflow_transform_method_offset), # requires catchment_data_set object to calculate bounds 
    minimise_likelihood = TRUE,
    streamflow_transform_method_offset = streamflow_transform_method_offset
  )

# Warnings generated from making the possible bound range of log-sinh a and b


# 3. Put numerical_optimiser object into a numerical optimiser -----------------

# currently only works with CMAES and DREAM functions
results <- numerical_optimiser |>
  CMAES(
    cmaes_control = list(), # can alter default cmaes controls here - see cmaesr package for options
    print_monitor = TRUE # print optimiser outputs in console while CMAES is running
  ) 


# 4. Convert results into a result_set object (standardised format) ------------
standardised_results <- results |> result_set() 


# 5. Examine results -----------------------------------------------------------
### probably should add ... in plot method so the user can alter the ggplot object directly i.e.,
### plot(standardised_results, type = "streamflow-time, ...) add theme(axis...)
plot(standardised_results, type = "streamflow-time")
plot(standardised_results, type = "rainfall-runoff")
plot(standardised_results, type = "examine_transform") 


# Move findings to another location

# The modelled streamflow produces normal-ish distributions
#hist(standardised_results$optimised_modelled_streamflow_transformed_space, breaks = 9)
#hist(standardised_results$optimised_modelled_streamflow_realspace, breaks = 9)

#test_data <- example_catchment$stop_start_data_set |> list_rbind()

# The near_bound column does not functioning correctly
# Ideally, it should scale with the parameters. It does not do this currently.
parameter_table <- parameters_summary(standardised_results)
tail(parameter_table)
result_table <- streamflow_timeseries_summary(standardised_results) 

# DREAM is a bit more complex...
# Is the reverse transform no good? No working correctly.
#xx <- seq(from = -20, to = 150, by = 0.1)
#yy <- log_sinh_transform(alpha = 76.1, beta = 9.89, y = xx)
#plot(xx, yy, type = "l")
#z <- inverse_log_sinh_transform(alpha = 76.1, beta = 9.89, z = yy)



log_sinh_asymptote <- function(a, b) {
  # Vertical asymptotes occur when log() approaches zero
  # sinh(a + bx) approaches zero
  # a + bx = 0
  # asymptote --> -a/b
  -a / b
}

find_acceptable_log_sinh_bounds <- function(min_offsetted_streamflow) {
  
  # Create possible range of values
  a <- -10^seq(from = -8, to = 2, by = 1)
  length_a <- length(a)
  b <- 10^seq(from = -8, to = 2, by = 1)
  length_b <- length(b)
  
  # Initial matrix
  matrix_of_values <- matrix(numeric(length_a * length_b), nrow = length_a, ncol = length_b) 
  
  
  for (i in 1:length_b) {
    for (j in 1:length_a) {
      matrix_of_values[j, i] <- log_sinh_asymptote(
        a = a[j], 
        b = b[i]
      ) 
    }
  }
  
  
  # The bounds are acceptable if:
  ## 1. log_sinh_asymptote is > 0
  ## 2. log_sinh_asymptote is < minimum observed flow with offset
  
  acceptable_matrix_values <- (matrix_of_values > 0) & (matrix_of_values < min_offsetted_streamflow)
  
  # Add names to acceptable_matrix_values
  colnames(acceptable_matrix_values) <- as.character(b)
  rownames(acceptable_matrix_values) <- as.character(a)
  
  
  return(acceptable_matrix_values)
}

x <- find_acceptable_log_sinh_bounds(1)


## log-sinh needs special treatment to find acceptable bounds
find_acceptable_log_sinh_bounds <- function(streamflow, offset) {
  
  # See log-sinh-bounds.xlsx for the method
  # From excel a is limited to 1E-8 and 0.1
  browser()
  
  # rows and cols for matrix
  possible_a <- sort(10^-seq(from = 1, to = 8, by = 1)) 
  names_possible_a <- as.character(possible_a)
  length_a <- length(possible_a)
  possible_b <- sort(10^-seq(from = -3, to = 8, by = 1))
  names_possible_b <- as.character(possible_b)
  length_b <- length(possible_b)
  
  
  matrix_of_values <- matrix(numeric(length_a * length_b), nrow = length_a, ncol = length_b) |> 
    `colnames<-`(names_possible_b) |> 
    `rownames<-`(names_possible_a)
  
  min_streamflow <- min(streamflow)
  
  for (i in 1:length_b) {
    for (j in 1:length_a) {
      matrix_of_values[j, i] <- log_sinh_transform(
        a = possible_a[j], 
        b = possible_b[i], 
        y = min_streamflow, 
        offset = offset
      ) 
    }
  }
  
  # Only keep values less than min_streamflow
  matrix_of_values[matrix_of_values > min_streamflow] <- NA
  
  # Filter out Inf
  matrix_of_values[is.infinite(matrix_of_values)] <- NA
  
  # 1E-8 is always good for b
  range_a <- possible_a |> range()
  
  # Use upper range_a to focus on row to get range_b
  range_b <- matrix_of_values[as.character(range_a[2]), ] |> 
    na.omit() |> 
    names() |> 
    as.numeric() |> 
    range()
  
  
  list(
    "range_a" = c(1E-8, 1E-5),#range_a,
    "range_b" = c(1E-5, 1)#range_b
  )
}




# Is it something related to standardised results? Yes
x <- standardised_results$numerical_optimiser_setup$catchment_data$stop_start_data_set |> 
  list_rbind() |> 
  pull(observed_streamflow) # realspace

x <- standardised_results$optimised_modelled_streamflow_realspace
y <- standardised_results$optimised_modelled_streamflow_transformed_space # what the model spits out

test_data <- test_data |> 
  add_column(
    "realspace_modelled_streamflow" = x,
    "transformed_modelled_streamflow" = y
  )

params <- standardised_results$best_parameter_set
curve_values <- log_sinh_transform(params["a"], params["b"], x)


plot(x, y, xlab = "realspace modelled streamflow", ylab = "transformed modelled streamflow")
points(x, curve_values, col = "red")


x <- seq(from = 1E-8, to = 50, by = 0.00001)
params <- standardised_results$best_parameter_set
curve_values <- log_sinh_transform(params["a"], params["b"], 50)
plot(x, curve_values, type = "l")

# I am not sure why the points are not on the line
# Points below the line indicate a negative streamflow value - YES
# I don't think this is good for the transform
# The realspace values do not have negative streamflow values
# This is a problem --> modelled streamflow cannot be negative in the realspace
# Solutions:
## 1. check inverse_log_sinh if less than zero for every iterations - if true set to NA/Inf
## 2. calc the minimum value log_sinh_transform(alpha = 76.1, beta = 9.89, y = 0)
##    transformed streamflow cannot be less than this - If true set to zero
## 3. allow negative values - but set to zero flow after calibration
## 4. offset angle?

## Problem the streamflow_model is producing values that do not exist on the
## log_sinh curve - for some catchments


# Get response surface for possible range of values of a and b
a_test <- c(10^-seq(from = 0, to = 8, by = 1), 10, 100)
a_test <- sort(c(-a_test, 0, a_test))
b_test <- sort(c(10^-seq(from = 0, to = 8, by = 1), 10, 100))
#transform_result <- log_sinh_transform(a = a_test, b = b_test, y = 0)

min_obs_realspace_q_mm <- 50

ggplot_tibble <- expand_grid(a_test, b_test) |> 
  mutate(
    result = log_sinh_transform(a = a_test, b = b_test, y = min_obs_realspace_q_mm)
  ) |> 
  mutate(
    greater_than_min_obs_q = if_else(result > min_obs_realspace_q_mm, TRUE, FALSE)
  )

# Learning:
# - If the minimum observed streamflow is 0 than `a` must be > 0
# - b must be greater than zero
zz <- ggplot_tibble |> filter(a_test < 0)
z <- ggplot_tibble |> filter(!is.na(result)) |> filter(is.finite(result))

signed_log10 <- scales::trans_new(
  name = "signed_log10",
  transform = function(x) sign(x) * log10(abs(x)),
  inverse = function(x) sign(x) * 10^abs(x)
)

ggplot_tibble |> 
  ggplot(aes(x = a_test, y = b_test, fill = greater_than_min_obs_q)) +
  geom_tile() +
  scale_x_continuous(trans = signed_log10, n.breaks = 20) +
  scale_y_log10(n.breaks = 20) +
  theme_bw()

# based off this graph the a and b bounds are:
# a = 1E-8 to 0.3 log-trans
# b = 1E-8 to 1 log-trans



y <- find_acceptable_log_sinh_bounds(1000)

x <- matrix(seq(from = 1, to = 9), ncol = 3, byrow = TRUE)
colnames(x) <- c("a", "b", "c")
rownames(x) <- c("1", "2", "3")
x[1, 1] <- NA
x[2, 1] <- NA
x[1, 3] <- NA
# I want to extract 5, 6, 8, 9 matrix
y <- !is.na(x)


