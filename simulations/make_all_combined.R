comb_pth = "~/Documents/evaluate_missing_surrogates/simulations/combine_results/"
comb_scr = paste0(comb_pth, list.files(comb_pth))
sapply(comb_scr, source)
