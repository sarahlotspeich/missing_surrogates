# Load libraries
library(ggplot2)
library(latex2exp)

# Source plot-making function from GitHub
devtools::source_url("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/figures/boxplot_of_estimates.R")

# Read in simulation results
sim_res = read.csv("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/simulations/sett5_non_overlap_smle_alphas.csv")

# Make them long 
res_long = sim_res |> 
  tidyr::pivot_longer(cols = smle_param_delta:smle_param_beta3, 
                      names_to = "method_quantity", 
                      values_to = "est") |> 
  dplyr::mutate(quantity = sub(pattern = ".*_", 
                               replacement = "", 
                               x = method_quantity), 
                truth = dplyr::case_when(
                  quantity == "delta" ~ 12, 
                  quantity == "delta.s" ~ 6, 
                  quantity == "R.s" ~ 0.5, 
                  quantity == "alpha0" ~ 5, 
                  quantity == "alpha1" ~ 6, 
                  quantity == "beta0" ~ 2, 
                  quantity == "beta1" ~ 1, 
                  quantity == "beta2" ~ 5, 
                  quantity == "beta3" ~ 1,
                  .default = NA), 
                quantity = factor(x = quantity, 
                                  levels = c("delta", "delta.s", "R.s", "alpha0", "alpha1", "beta0", "beta1", "beta2", "beta3"), 
                                  labels = c(TeX("$\\Delta$"), TeX("$\\Delta_S$"), TeX("$R_S$"), 
                                             TeX("$\\alpha_0$"), TeX("$\\alpha_1$"), TeX("$\\beta_0$"),
                                             TeX("$\\beta_1$"), TeX("$\\beta_2$"), TeX("$\\beta_3$")))) |> 
  dplyr::select(-method_quantity) |> 
  dplyr::mutate(parametric = "PTE Estimator: Parametric", 
                method = "SMLE")

# Make a boxplot 
cols = c("#787ff6", "#ffbd59", "#8bdddb", "#ff99ff",  "#7dd5f6") ## color palette
res_long |> 
  dplyr::filter(quantity == "alpha[0]" | quantity == "alpha[1]") |> 
  ggplot(aes(x = method, y = est, fill = overlap)) + 
  geom_hline(aes(yintercept = truth), linetype = 2, color = "black") + 
  geom_boxplot() + 
  xlab("Method") + 
  ylab("Estimate") + 
  facet_grid(cols =  vars(quantity), 
             scales = "free", 
             labeller = labeller(quantity = label_parsed)) + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 5)) + 
  scale_fill_manual(name = "Method:", values = cols) + 
  theme_minimal() + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "black"), 
        strip.text = element_text(color = "white"))  

ggsave(filename = "figures/sett5_overlap_non_overlap_boxplot.pdf", 
       device = "pdf", width = 7, height = 5)
