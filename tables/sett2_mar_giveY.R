# Load libraries
library(kableExtra) ## for pretty LaTex tables

# Source table-making function from GitHub
devtools::source_url("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/tables/table_of_estimates.R")

# Read in simulation results
sim_res = read.csv("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/simulations/sett2_mar_givY.csv")

# Make table 
sim_res |> 
  table_of_estimates()
