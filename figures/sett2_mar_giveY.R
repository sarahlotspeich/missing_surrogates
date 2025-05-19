# Load libraries
library(ggplot2)
library(latex2exp)

# Read in simulation results
p = paste0("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/simulations/sett2_mar_givY/sett2_mar_givY_seed", 0:9, ".csv")
sim_res = do.call(dplyr::bind_rows, 
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
  dplyr::mutate(parametric = factor(x = !grepl(pattern = "nonparam", x = method), 
                                    levels = c(FALSE, TRUE), 
                                    labels = c("Nonparametric", "Parametric")), 
                method = factor(x = method, 
                                levels = c("gs_nonparam", "gs_param", 
                                           "cc_nonparam", "cc_param", 
                                           "ipw_nonparam", "ipw_param", 
                                           "smle_param", "mle_param"), 
                                labels = c("Gold Standard", "Gold Standard", 
                                           "Complete Case", "Complete Case", 
                                           "IPW", "IPW",
                                           "SMLE", "MLE")))

# Make a boxplot 
res_long |> 
  ggplot(aes(x = method, y = est, fill = parametric)) + 
  geom_boxplot() + 
  geom_hline(aes(yintercept = truth), linetype = 2, color = "white") + 
  facet_wrap(~quantity, scales = "free", ncol = 3, labeller = label_parsed) + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 5)) + 
  scale_fill_manual(name = "PTE Estimator:", values = c("#ffbd59",  "#787ff6")) + 
  theme_minimal() + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "black"), 
        strip.text = element_text(color = "white"))
ggsave(filename = "~/Documents/missing_surrogates/figures/sett2_mar_givY_boxplot.pdf", 
       device = "pdf", width = 7, height = 5)
