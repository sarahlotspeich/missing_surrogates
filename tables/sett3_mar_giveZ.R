# Load libraries
library(kableExtra) ## for pretty LaTex tables

# Source table-making function from GitHub
devtools::source_url("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/tables/table_of_estimates.R")

# Read in simulation results from GitHub
p = paste0("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/simulations/sett3_mar_givZ/sett3_mar_givZ_seed", 0:19, ".csv")
sim_res = do.call(dplyr::bind_rows, 
                  lapply(X = paste0(p, list.files(p)), 
                         FUN = read.csv)) 

# Make table 
sim_res |> 
  table_of_estimates()