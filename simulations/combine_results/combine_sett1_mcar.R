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

## Save combined data 
sim_res |> 
  write.csv("simulations/sett1_mcar.csv", 
            row.names = FALSE)
