# Load libraries
library(ggplot2)
library(latex2exp)

# Source plot-making function from GitHub
devtools::source_url("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/figures/boxplot_of_estimates.R")

# Read in simulation results
sim_res = read.csv("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/simulations/sett2_mar_givY.csv") |> 
  dplyr::bind_cols(data.frame(seed = rep(x = 0:19, each = 50))) |> 
  dplyr::select(-seed, -dplyr::contains(c("ci", "var")))

# Make them long 
res_long = sim_res |> 
  tidyr::pivot_longer(cols = gs_nonparam_delta:smle_param_R.s, 
                      names_to = "method_quantity", 
                      values_to = "est") |> 
  dplyr::mutate(quantity = sub(pattern = ".*_", 
                               replacement = "", 
                               x = method_quantity), 
                truth = dplyr::case_when(
                  quantity == "delta" ~ 12, 
                  quantity == "delta.s" ~ 6, 
                  .default = 0.5), 
                quantity = factor(x = quantity, 
                                  levels = c("delta", "delta.s", "R.s"), 
                                  labels = c(TeX("Quantity: $\\Delta$"), TeX("Quantity: $\\Delta_S$"), TeX("Quantity: $R_S$"))),
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
                                    labels = c("PTE Estimator: Nonparametric",
                                               "PTE Estimator: Parametric")), 
                method = factor(x = method, 
                                levels = c("gs_nonparam", "cc_nonparam", "ipw_nonparam",
                                           "gs_param", "cc_param", "ipw_param", "smle_param"), 
                                labels = c("Gold Standard", "Complete Case", "IPW",
                                           "Gold Standard",  "Complete Case", "IPW",  "SMLE")))

# Make a boxplot 
res_long |> 
  boxplot_of_estimates()
ggsave(filename = "figures/sett2_mar_givY_boxplot.pdf", 
       device = "pdf", width = 7, height = 5)
