# Load libraries
library(kableExtra) ## for pretty LaTex tables

# Source table-making function from GitHub
devtools::source_url("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/tables/table_of_estimates.R")

# Read in simulation results from GitHub
sim_res = read.csv("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/simulations/sett3_mar_givZ.csv")

# Make table 
sim_res |> 
  table_of_estimates()