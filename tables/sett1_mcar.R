# Load libraries
library(kableExtra) ## for pretty LaTex tables

# Source plot-making function from GitHub
# devtools::source_url("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/figures/boxplot_of_estimates.R")

# Read in simulation results from GitHub
p = paste0("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/simulations/sett1_mcar/sett1_mcar_seed", 0:19, ".csv")
sim_res = do.call(dplyr::bind_rows, 
                  lapply(X = paste0(p, list.files(p)), 
                         FUN = read.csv))

# Calculate empirical bias for each method
emp_bias = sim_res |> 
  dplyr::select(dplyr::contains("R.s"), 
                -dplyr::contains(c("var", "ci"))) |> 
  tidyr::pivot_longer(cols = gs_nonparam_R.s:smle_param_R.s, 
                      names_to = "method", values_to = "est") |> 
  dplyr::group_by(method) |> 
  dplyr::summarize(bias = mean(est - 0.5)) |>
  dplyr::mutate(perc_bias = round(100 * bias / 0.5, 1)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(method = sub(pattern = "_R.s", replacement = "", x = method))

# Calculate empirical SEs for each method
emp_se = sim_res |> 
  dplyr::select(dplyr::contains("R.s"), 
                -dplyr::contains(c("var", "ci"))) |> 
  tidyr::pivot_longer(cols = gs_nonparam_R.s:smle_param_R.s, 
                      names_to = "method", values_to = "est") |> 
  dplyr::group_by(method) |> 
  dplyr::summarize(se = sd(est)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(method = sub(pattern = "_R.s", replacement = "", x = method))

# Calculate average SE estimator for each 
avg_se = sim_res |> 
  dplyr::select(dplyr::contains("var_R.s")) |> 
  tidyr::pivot_longer(cols = gs_nonparam_var_R.s:smle_param_var_R.s, 
                      names_to = "method", values_to = "var") |> 
  dplyr::group_by(method) |> 
  dplyr::summarize(see = mean(sqrt(var))) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(method = sub(pattern = "_var_R.s", replacement = "", x = method))

# Calculate coverage with normal approx CIs 
norm_cp = sim_res |> 
  dplyr::mutate(gs_nonparam_normci_cov = gs_nonparam_normci_lb_R.s <= 0.5 & 
                  0.5 <= gs_nonparam_normci_ub_R.s, 
                gs_param_normci_cov = gs_param_normci_lb_R.s <= 0.5 & 
                  0.5 <= gs_param_normci_ub_R.s, 
                cc_nonparam_normci_cov = cc_nonparam_normci_lb_R.s <= 0.5 & 
                  0.5 <= cc_nonparam_normci_ub_R.s, 
                cc_param_normci_cov = cc_param_normci_lb_R.s <= 0.5 & 
                  0.5 <= cc_param_normci_ub_R.s, 
                ipw_nonparam_normci_cov = ipw_nonparam_normci_lb_R.s <= 0.5 & 
                  0.5 <= ipw_nonparam_normci_ub_R.s, 
                ipw_param_normci_cov = ipw_param_normci_lb_R.s <= 0.5 & 
                  0.5 <= ipw_param_normci_ub_R.s, 
                smle_param_normci_cov = smle_param_normci_lb_R.s <= 0.5 & 
                  0.5 <= smle_param_normci_ub_R.s) |> 
  dplyr::select(dplyr::contains("normci_cov")) |> 
  tidyr::pivot_longer(cols = gs_nonparam_normci_cov:smle_param_normci_cov, 
                      names_to = "method", values_to = "covered") |> 
  dplyr::group_by(method) |> 
  dplyr::summarize(norm_cp = mean(covered)) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(method = sub(pattern = "_normci_cov", replacement = "", x = method))

# Calculate coverage with quantile CIs 
quant_cp = sim_res |> 
  dplyr::mutate(gs_nonparam_quantci_cov = gs_nonparam_quantci_lb_R.s <= 0.5 & 
                  0.5 <= gs_nonparam_quantci_ub_R.s, 
                gs_param_quantci_cov = gs_param_quantci_lb_R.s <= 0.5 & 
                  0.5 <= gs_param_quantci_ub_R.s, 
                cc_nonparam_quantci_cov = cc_nonparam_quantci_lb_R.s <= 0.5 & 
                  0.5 <= cc_nonparam_quantci_ub_R.s, 
                cc_param_quantci_cov = cc_param_quantci_lb_R.s <= 0.5 & 
                  0.5 <= cc_param_quantci_ub_R.s, 
                ipw_nonparam_quantci_cov = ipw_nonparam_quantci_lb_R.s <= 0.5 & 
                  0.5 <= ipw_nonparam_quantci_ub_R.s, 
                ipw_param_quantci_cov = ipw_param_quantci_lb_R.s <= 0.5 & 
                  0.5 <= ipw_param_quantci_ub_R.s, 
                smle_param_quantci_cov = smle_param_quantci_lb_R.s <= 0.5 & 
                  0.5 <= smle_param_quantci_ub_R.s) |> 
  dplyr::select(dplyr::contains("quantci_cov")) |> 
  tidyr::pivot_longer(cols = gs_nonparam_quantci_cov:smle_param_quantci_cov, 
                      names_to = "method", values_to = "covered") |> 
  dplyr::group_by(method) |> 
  dplyr::summarize(quant_cp = mean(covered))  |> 
  dplyr::ungroup() |> 
  dplyr::mutate(method = sub(pattern = "_quantci_cov", replacement = "", x = method))

# Write function to add "padded" zeros and wrap with $$ for consistency 
format_num = function(num) {
  paste0("$", format(round(num, 3), nsmall = 3), "$")
}

# Make a table (formatted for LaTex)
emp_bias |> 
  dplyr::left_join(emp_se) |> 
  dplyr::left_join(avg_se) |> 
  dplyr::left_join(norm_cp) |> 
  dplyr::left_join(quant_cp) |> 
  dplyr::mutate(
    method = factor(x = method, 
                    levels = c("gs_nonparam", "cc_nonparam", "ipw_nonparam", 
                               "gs_param", "cc_param", "ipw_param", "smle_param"), 
                    labels = c("Gold Standard (NP)", "Complete Case (NP)", "IPW (NP)", 
                               "Gold Standard (P)", "Complete Case (P)", "IPW (P)", "SMLE"))
  ) |> 
  dplyr::arrange(method) |> 
  mutate_at(.vars = c("bias", "perc_bias", "se", "see", "norm_cp", "quant_cp"), 
            .funs = format_num) |>
  kable(format = "latex", 
        booktabs = TRUE, 
        escape = FALSE, 
        align = "rrrcccc", 
        col.names = c("Method", "Bias", "\\% Bias", "Emp.", "Avg.", "Norm.", "Quant.")) |> 
  kable_styling() |> 
  add_header_above(header = c(" " = 3, 
                              "Standard Errors" = 2, 
                              "Coverage Probability" = 2), 
                   bold = TRUE) |>
  group_rows(group_label = "PTE Estimator: Nonparametric", start_row = 1, end_row = 3, italic = TRUE) |> 
  group_rows(group_label = "PTE Estimator: Parametric", start_row = 4, end_row = 7, italic = TRUE)
