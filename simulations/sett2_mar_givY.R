# Libraries and functions
## Package with PTE estimation / missing data corrections
library(missSurrogate)
## Source data generation function (surrogate distributions overlap)
devtools::source_url("https://raw.githubusercontent.com/sarahlotspeich/evaluate_missing_surrogates/refs/heads/main/gen.data.overlap.R")

# Reproducibility 
set.seed(11422) 

# Set sample sizes 
n1 = 1000
n0 = 1000

# Run simulations
REPS = 1000 ### but we ran 50 each across 20 arrays on a cluster for efficiency
## Initialize empty dataframe for results
sim_res = data.frame(
  r = 1:REPS, 
  gs_nonparam_delta = NA, gs_nonparam_delta.s = NA, gs_nonparam_R.s = NA, gs_nonparam_var_R.s = NA, gs_nonparam_normci_lb_R.s = NA, gs_nonparam_normci_ub_R.s = NA, gs_nonparam_quantci_lb_R.s = NA, gs_nonparam_quantci_ub_R.s = NA,
  gs_param_delta = NA, gs_param_delta.s = NA, gs_param_R.s = NA, gs_param_var_R.s = NA, gs_param_normci_lb_R.s = NA, gs_param_normci_ub_R.s = NA, gs_param_quantci_lb_R.s = NA, gs_param_quantci_ub_R.s = NA,
  cc_nonparam_delta = NA, cc_nonparam_delta.s = NA, cc_nonparam_R.s = NA, cc_nonparam_var_R.s = NA, cc_nonparam_normci_lb_R.s = NA, cc_nonparam_normci_ub_R.s = NA, cc_nonparam_quantci_lb_R.s = NA, cc_nonparam_quantci_ub_R.s = NA,
  cc_param_delta = NA, cc_param_delta.s = NA, cc_param_R.s = NA, cc_param_var_R.s = NA, cc_param_normci_lb_R.s = NA, cc_param_normci_ub_R.s = NA, cc_param_quantci_lb_R.s = NA, cc_param_quantci_ub_R.s = NA,
  ipw_nonparam_delta = NA, ipw_nonparam_delta.s = NA, ipw_nonparam_R.s = NA, ipw_nonparam_var_R.s = NA, ipw_nonparam_normci_lb_R.s = NA, ipw_nonparam_normci_ub_R.s = NA, ipw_nonparam_quantci_lb_R.s = NA, ipw_nonparam_quantci_ub_R.s = NA,
  ipw_param_delta = NA, ipw_param_delta.s = NA, ipw_param_R.s = NA, ipw_param_var_R.s = NA, ipw_param_normci_lb_R.s = NA, ipw_param_normci_ub_R.s = NA, ipw_param_quantci_lb_R.s = NA, ipw_param_quantci_ub_R.s = NA, 
  smle_param_delta = NA, smle_param_delta.s = NA, smle_param_R.s = NA, smle_param_var_R.s = NA, smle_param_normci_lb_R.s = NA, smle_param_normci_ub_R.s = NA, smle_param_quantci_lb_R.s = NA, smle_param_quantci_ub_R.s = NA) 

## Loop over replications
## Initialize empty dataframe for results
sim_res = data.frame(
  r = 1:REPS, 
  gs_nonparam_delta = NA, gs_nonparam_delta.s = NA, gs_nonparam_R.s = NA, gs_nonparam_var_R.s = NA, gs_nonparam_normci_lb_R.s = NA, gs_nonparam_normci_ub_R.s = NA, gs_nonparam_quantci_lb_R.s = NA, gs_nonparam_quantci_ub_R.s = NA,
  gs_param_delta = NA, gs_param_delta.s = NA, gs_param_R.s = NA, gs_param_var_R.s = NA, gs_param_normci_lb_R.s = NA, gs_param_normci_ub_R.s = NA, gs_param_quantci_lb_R.s = NA, gs_param_quantci_ub_R.s = NA,
  cc_nonparam_delta = NA, cc_nonparam_delta.s = NA, cc_nonparam_R.s = NA, cc_nonparam_var_R.s = NA, cc_nonparam_normci_lb_R.s = NA, cc_nonparam_normci_ub_R.s = NA, cc_nonparam_quantci_lb_R.s = NA, cc_nonparam_quantci_ub_R.s = NA,
  cc_param_delta = NA, cc_param_delta.s = NA, cc_param_R.s = NA, cc_param_var_R.s = NA, cc_param_normci_lb_R.s = NA, cc_param_normci_ub_R.s = NA, cc_param_quantci_lb_R.s = NA, cc_param_quantci_ub_R.s = NA,
  ipw_nonparam_delta = NA, ipw_nonparam_delta.s = NA, ipw_nonparam_R.s = NA, ipw_nonparam_var_R.s = NA, ipw_nonparam_normci_lb_R.s = NA, ipw_nonparam_normci_ub_R.s = NA, ipw_nonparam_quantci_lb_R.s = NA, ipw_nonparam_quantci_ub_R.s = NA,
  ipw_param_delta = NA, ipw_param_delta.s = NA, ipw_param_R.s = NA, ipw_param_var_R.s = NA, ipw_param_normci_lb_R.s = NA, ipw_param_normci_ub_R.s = NA, ipw_param_quantci_lb_R.s = NA, ipw_param_quantci_ub_R.s = NA) 
for (r in 1:REPS) {
  # Generate data 
  data = gen.data.overlap(n1 = n1, 
                          n0 = n0) 
  
  # Define vectors for outcomes/surrogates in untreated/treated 
  s1 = data$s1
  y1 = data$y1
  s0 = data$s0
  y0 = data$y0
  
  ##########################################################
  #Estimates with complete data ############################
  ##########################################################
  ## Estimate R with nonparametric approach (gold standard)
  Rnonparam = R.s.miss(sone=s1, 
                       szero=s0, 
                       yone = y1, 
                       yzero = y0, 
                       type = "robust", 
                       conf.int = TRUE)
  sim_res[r, c("gs_nonparam_delta", "gs_nonparam_delta.s", "gs_nonparam_R.s")] = with(Rnonparam, c(delta, delta.s, R.s))
  sim_res[r, c("gs_nonparam_var_R.s", "gs_nonparam_normci_lb_R.s", "gs_nonparam_normci_ub_R.s", 
               "gs_nonparam_quantci_lb_R.s", "gs_nonparam_quantci_ub_R.s")] = with(Rnonparam, c(R.s.var, conf.int.normal.R.s, conf.int.quantile.R.s))
  
  ## Estimate R with parametric approach (gold standard) ###
  Rparam = R.s.miss(sone=s1, 
                    szero=s0, 
                    yone = y1, 
                    yzero = y0, 
                    type = "model", 
                    conf.int = TRUE)
  sim_res[r, c("gs_param_delta", "gs_param_delta.s", "gs_param_R.s")] = with(Rparam, c(delta, delta.s, R.s))
  sim_res[r, c("gs_param_var_R.s", "gs_param_normci_lb_R.s", "gs_param_normci_ub_R.s", 
               "gs_param_quantci_lb_R.s", "gs_param_quantci_ub_R.s")] = with(Rparam, c(R.s.var, conf.int.normal.R.s, conf.int.quantile.R.s))
  
  # Simulate non-missingness indicators ###################
  ## Under MAR, probability of missingness depends on Y continuously (logistic regression)
  m1 = rbinom(n = n1, size = 1, prob = 1 / (1 + exp(- 0.015 * y1)))
  m0 = rbinom(n = n0, size = 1, prob = 1 / (1 + exp(- 0.015 * y0)))
  s0[m0==0] = NA ### make them missing
  s1[m1==0] = NA ### make them missing
  
  ##########################################################
  #Estimates with incomplete data ##########################
  ##########################################################
  ## Estimate R with nonparametric approach (complete case)
  Rnonparam_miss = R.s.miss(sone = s1[m1==1], 
                            szero = s0[m0==1], 
                            yone = y1[m1==1], 
                            yzero = y0[m0==1], 
                            type = "robust", 
                            conf.int = TRUE)
  sim_res[r, c("cc_nonparam_delta", "cc_nonparam_delta.s", "cc_nonparam_R.s")] = with(Rnonparam_miss, c(delta, delta.s, R.s))
  sim_res[r, c("cc_nonparam_var_R.s", "cc_nonparam_normci_lb_R.s", "cc_nonparam_normci_ub_R.s", 
               "cc_nonparam_quantci_lb_R.s", "cc_nonparam_quantci_ub_R.s")] = with(Rnonparam_miss, c(R.s.var, conf.int.normal.R.s, conf.int.quantile.R.s))
  
  ## Estimate R with parametric approach (complete case)
  Rparam_miss = R.s.miss(sone = s1[m1==1], 
                         szero = s0[m0==1], 
                         yone = y1[m1==1], 
                         yzero = y0[m0==1], 
                         type = "model", 
                         conf.int = TRUE)
  sim_res[r, c("cc_param_delta", "cc_param_delta.s", "cc_param_R.s")] = with(Rparam_miss, c(delta, delta.s, R.s))
  sim_res[r, c("cc_param_var_R.s", "cc_param_normci_lb_R.s", "cc_param_normci_ub_R.s", 
               "cc_param_quantci_lb_R.s", "cc_param_quantci_ub_R.s")] = with(Rparam_miss, c(R.s.var, conf.int.normal.R.s, conf.int.quantile.R.s))
  
  ## Calculate weights for IPW approaches
  m = c(m1, m0)
  y = c(y1, y0)
  ipw_fit = glm(formula = m ~ y, 
                family = "binomial")
  p1 = 1 / (1 + exp(-(ipw_fit$coefficients[1] + ipw_fit$coefficients[2] * y1)))
  w1 = 1 / p1
  p0 = 1 / (1 + exp(-(ipw_fit$coefficients[1] + ipw_fit$coefficients[2] * y0)))
  w0 = 1 / p0  
  
  ## Estimate R with nonparametric approach (IPW)
  Rnonparam_miss_ipw = R.s.miss(sone = s1, 
                                szero = s0, 
                                yone = y1, 
                                yzero = y0,
                                type = "robust",
                                wone = w1, 
                                wzero = w0, 
                                ipw.formula = m ~ y, 
                                conf.int = TRUE)
  sim_res[r, c("ipw_nonparam_delta", "ipw_nonparam_delta.s", "ipw_nonparam_R.s")] = with(Rnonparam_miss_ipw, c(delta, delta.s, R.s))
  sim_res[r, c("ipw_nonparam_var_R.s", "ipw_nonparam_normci_lb_R.s", "ipw_nonparam_normci_ub_R.s", "ipw_nonparam_quantci_lb_R.s", "ipw_nonparam_quantci_ub_R.s")] = with(Rnonparam_miss_ipw, c(R.s.var, conf.int.normal.R.s, conf.int.quantile.R.s))
  
  ## Estimate R with parametric approach (IPW)
  Rparam_miss_ipw = R.s.miss(sone = s1, 
                             szero = s0, 
                             yone = y1, 
                             yzero = y0, 
                             type = "model",
                             wone = w1, 
                             wzero = w0, 
                             ipw.formula = m ~ y, 
                             conf.int = TRUE)
  sim_res[r, c("ipw_param_delta", "ipw_param_delta.s", "ipw_param_R.s")] = with(Rparam_miss_ipw, c(delta, delta.s, R.s))
  sim_res[r, c("ipw_param_var_R.s", "ipw_param_normci_lb_R.s", "ipw_param_normci_ub_R.s", 
               "ipw_param_quantci_lb_R.s", "ipw_param_quantci_ub_R.s")] = with(Rparam_miss_ipw, c(R.s.var, conf.int.normal.R.s, conf.int.quantile.R.s))
  
  ## Estimate R with parametric approach (SMLE)
  Rparam_miss_smle = R.s.miss(sone = s1, 
                              szero = s0,
                              yone = y1,
                              yzero = y0, 
                              type = "model", 
                              conf.int = TRUE) 
  sim_res[r, c("smle_param_delta", "smle_param_delta.s", "smle_param_R.s")] = with(Rparam_miss_smle, c(delta, delta.s, R.s))
  sim_res[r, c("smle_param_var_R.s", "smle_param_normci_lb_R.s", "smle_param_normci_ub_R.s", 
               "smle_param_quantci_lb_R.s", "smle_param_quantci_ub_R.s")] = with(Rparam_miss_smle, c(R.s.var, conf.int.normal.R.s, conf.int.quantile.R.s))
  
  ## Save 
  sim_res |> 
    write.csv(paste0("sett2_mar_givY/sett2_mar_givY_seed", sim_seed, ".csv"), 
              row.names = FALSE)
}
