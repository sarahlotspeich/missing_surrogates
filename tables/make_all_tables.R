tab_pth = "~/Documents/evaluate_missing_surrogates/tables/"
tab_scr = paste0(tab_pth, list.files(tab_pth))
sapply(tab_scr[-c(1, 6)], source)
