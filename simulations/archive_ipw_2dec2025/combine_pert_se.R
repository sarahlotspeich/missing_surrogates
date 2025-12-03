## Gold standard + complete case + IPW 
p = paste0("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/simulations/pert_se/sett2_mar_givY_seed", 0:19, ".csv")
sim_res = do.call(dplyr::bind_rows, 
                  lapply(X = paste0(p, list.files(p)), 
                         FUN = read.csv))
# SMLE
p = paste0("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/simulations/pert_se_smle/sett2_mar_givY_seed", 20:39, ".csv")
sim_res = sim_res |>
  dplyr::bind_cols(
    do.call(dplyr::bind_rows,
            lapply(X = paste0(p, list.files(p)),
                   FUN = read.csv)) |>
      dplyr::select(dplyr::starts_with("smle"))
  )

## Save combined data 
sim_res |> 
  write.csv("simulations/pert_se.csv", 
            row.names = FALSE)
