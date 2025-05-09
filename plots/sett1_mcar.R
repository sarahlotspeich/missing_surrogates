# Read in simulation results
p = paste0("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/simulations/sett1_mcar/sett1_mcar_seed", 0:9, ".csv")
plot_dat = do.call(dplyr::bind_rows, 
                   lapply(X = paste0(p, list.files(p)), 
                          FUN = read.csv))

# Make them long 
res_long = sim_res |> 
  tidyr::pivot_longer(cols = gs_nonparam_delta:smle_param_R.s, 
                      names_to = "method_quantity", values_to = "est") |> 
  dplyr::mutate(quantity = sub(pattern = ".*_", 
                               replacement = "", 
                               x = method_quantity), 
                truth = dplyr::case_when(
                  quantity == "delta" ~ 12, 
                  quantity == "delta.s" ~ 6, 
                  .default = 0.5), 
                quantity = factor(x = quantity, 
                                  levels = c("delta", "delta.s", "R.s"), 
                                  labels = c(TeX("$\\Delta$"), TeX("$\\Delta_S$"), TeX("$R_S$"))),
                method = sub(pattern = "_delta", 
                             replacement = "", 
                             x = sub(pattern = "_delta.s", 
                                     replacement = "", 
                                     x = sub(pattern = "_R.s", 
                                             replacement = "", 
                                             x = method_quantity)))) |> 
  dplyr::select(-method_quantity) |> 
  dplyr::mutate(method = factor(x = method, 
                                levels = c("gs_nonparam", "gs_param", 
                                           "cc_nonparam", "cc_param", 
                                           "ipw_nonparam", "ipw_param", 
                                           "smle_param"), 
                                labels = c("GS (NP)", "GS (P)", 
                                           "CC (NP)", "CC (P)", 
                                           "IPW (NP)", "IPW (P)",
                                           "SMLE (P)")), 
                parametric = !grepl(pattern = "(NP)", 
                                    x = method))

# Make a boxplot 
res_long |> 
  ggplot(aes(x = method, y = est, fill = parametric)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = truth), linetype = 2, color = "white") + 
  facet_wrap(~quantity, scales = "free", ncol = 3, labeller = label_parsed) + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 5)) + 
  theme_minimal() + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "black"), 
        strip.text = element_text(color = "white")) + 
  ggtitle(label = "Boxplot of estimates under missingness completely at random (MCAR)")