# Convert a5 time of emergence parameter into a year 


# Import libraries--------------------------------------------------------------
pacman::p_load(tidyverse, furrr, parallel, truncnorm, sloop, arrow)


# Extract only ToE parameter from DREAM sequences ------------------------------
a5_DREAM_sequence_data <- open_dataset(
  source = "./Modelling/Results/DREAM/Sequences",
  format = "parquet"
) |> 
  filter(parameter == "a5") |> 
  collect()



# Get CO2 data -----------------------------------------------------------------
## Used to convert the a5 parameter into year
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



# Creating functions to convert a5 into year -----------------------------------
single_a5_to_year <- function(a5, CO2, year) {
  year[which(CO2 - a5 > 0)[1]]
}


a5_to_year_factory <- function(CO2, year) {
  # Save CO2 and year vectors to the function itself
  # They do not change
  
  force(CO2)
  force(year)
  stopifnot(is.numeric(CO2))
  stopifnot(is.numeric(year))
  
  function(a5) {
    future_map_dbl(
      .x = a5,
      .f = single_a5_to_year,
      CO2 = CO2,
      year = year,
      .progress = TRUE
    )
  }
}

# Function factory
a5_to_year <- a5_to_year_factory(CO2 = CO2, year = year)


# Convert a5 to year -----------------------------------------------------------
plan(multisession, workers = length(availableWorkers())) # set once for furrr
year_DREAM_sequence_data <- a5_DREAM_sequence_data |> 
  mutate(
    ToE = a5_to_year(parameter_value)
  )


# Save results -----------------------------------------------------------------
write_parquet(
  x = year_DREAM_sequence_data,
  sink = "./Modelling/Results/DREAM/year_DREAM_sequence_data.parquet"
)

