# Load libraries
library(ggplot2) ## for pretty plots
library(latex2exp) ## for LaTeX labels 

# Source plot-making function from GitHub
devtools::source_url("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/figures/boxplot_of_estimates.R")

# Read in simulation results from GitHub
sim_res = read.csv("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/simulations/sett4_misspec_ipw.csv")

# Make them long 
res_long = sim_res |> 
  tidyr::pivot_longer(cols = gs_nonparam:smle_param, 
                      names_to = "method", values_to = "est") |> 
  dplyr::mutate(parametric = factor(x = !grepl(pattern = "nonparam", x = method), 
                                    levels = c(FALSE, TRUE), 
                                    labels = c("PTE Estimator: Nonparametric",
                                               "PTE Estimator: Parametric")), 
                method = factor(x = method, 
                                levels = c("gs_nonparam", "cc_nonparam", "ipw_nonparam_Yonly", "ipw_nonparam_Zonly", "ipw_nonparam_YZ",
                                           "gs_param", "cc_param", "ipw_param_Yonly", "ipw_param_Zonly", "ipw_param_YZ","smle_param"), 
                                labels = c("Gold\nStandard", "Complete\nCase", "IPW\n(Y Only)", "IPW\n(Z Only)", "IPW\n(Y and Z)",
                                           "Gold\nStandard",  "Complete\nCase", "IPW\n(Y Only)", "IPW\n(Z Only)", "IPW\n(Y and Z)", "SMLE")))

# Make a boxplot 
cols = c("#787ff6", "#ffbd59", "#8bdddb", "#ff99ff",  "#7dd5f6", "#ff914d") ## color palette
res_long |> 
  ggplot(aes(x = method, y = est, fill = method)) + 
  geom_hline(yintercept = 0.5, linetype = 2, color = "black") + 
  geom_boxplot() + 
  xlab("Method") + 
  ylab("Estimate") + 
  facet_grid(cols = vars(parametric), 
             scales = "free", 
             labeller = labeller(parametric = label_value)) + 
  #scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) + 
  scale_fill_manual(name = "Method:", values = cols, guide = "none") + 
  theme_minimal() + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "black"), 
        strip.text = element_text(color = "white"))
ggsave(filename = "figures/sett4_misspec_ipw_boxplot.pdf", 
       device = "pdf", width = 10, height = 5)
