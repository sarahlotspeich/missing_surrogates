# Load libraries
library(kableExtra) ## for pretty LaTex tables
library(dplyr) ## for table making

# Source table-making function from GitHub
devtools::source_url("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/tables/table_of_estimates.R")

# Read in simulation results from GitHub
## Gold standard + complete case + IPW 
p = paste0("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/simulations/sett1_mcar/sett1_mcar_seed", 0:19, ".csv")
sim_res = do.call(dplyr::bind_rows, 
                  lapply(X = paste0(p, list.files(p)), 
                         FUN = read.csv))
## SMLE 
p = paste0("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/simulations/sett1_mcar_smle/sett1_mcar_seed", 20:39, ".csv")
sim_res = sim_res |> 
  dplyr::bind_cols(
    do.call(dplyr::bind_rows, 
            lapply(X = paste0(p, list.files(p)), 
                         FUN = read.csv)) |> 
      dplyr::select(dplyr::starts_with("smle"))
  )

# Make table 
sim_res |> 
  table_of_estimates()
