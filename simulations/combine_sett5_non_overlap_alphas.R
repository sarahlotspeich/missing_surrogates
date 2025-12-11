## SMLE 
#p = paste0("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/simulations/sett5_non_overlap_smle/sett5_non_overlap_seed", 20:39, ".csv")
p = paste0("~/Documents/missing_surrogates/simulations/sett5_non_overlap_smle/alphas_sett5_non_overlap_seed", 20:39, ".csv")
sim_res = dplyr::bind_cols(
    do.call(dplyr::bind_rows, 
            lapply(X = paste0(p, list.files(p)), 
                   FUN = read.csv)) |> 
      dplyr::select(dplyr::starts_with("smle"))
  )

## Save combined data 
sim_res |> 
  write.csv("simulations/sett5_non_overlap_smle_alphas.csv", 
            row.names = FALSE)
