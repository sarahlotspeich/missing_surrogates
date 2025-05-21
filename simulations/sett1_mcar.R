# Libraries and functions
library(Rsurrogate)
source("R.s.miss.R")

# Reproducibility 
## Random seed to be used for each simulation setting
args = commandArgs(TRUE)
## When running on the cluster, give each array a unique seed by using the array ID
sim_seed = as.integer(args)
## Be reproducible queens
set.seed(sim_seed) 

# Functions to generate data 
gen.data = function(setting, n1, n0) {
  s1 = g.1(n1)
  y1 = f.cond.1(s1)
  s0 = g.0(n0)
  y0 = f.cond.0(s0)
  return(data.frame("s1" = s1, "y1" = y1, "s0" = s0, "y0" = y0))
}
f.cond.1 = function(s.vector) {
  eps1 = rnorm(length(s.vector),0,3)
  y1 = 2+5*s.vector+1 + 1*s.vector + eps1
  return(y1)		
}
f.cond.0 = function(s.vector) {
  eps0 = rnorm(length(s.vector),0,3)
  y0 = 2+5*s.vector+ eps0
  return(y0)		
}
g.1 = function(n, alpha0=5) { return(rnorm(n, alpha0 + 1,2))}
g.0 = function(n,alpha0=5) { return(rnorm(n, alpha0,1))}

# Run simulations 
## Set number of replications per array
REPS = 100
## Set sample sizes 
n1 = 1000
n0 = 1000
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
for (r in 1:REPS) {
  # Generate data 
  data = gen.data(n1=n1, n0=n0) 
  
  # Define vectors for outcomes/surrogates in untreated/treated 
  s1 = data$s1
  y1 = data$y1
  s0 = data$s0
  y0 = data$y0
  
  ##########################################################
  #Estimates with complete data ############################
  ##########################################################
  ## Estimate R with nonparametric approach (gold standard)
  Rnonparam = R.s.estimate(sone=s1, 
                           szero=s0, 
                           yone = y1, 
                           yzero = y0, 
                           type = "robust", 
                           conf.int = TRUE)
  sim_res[r, c("gs_nonparam_delta", "gs_nonparam_delta.s", "gs_nonparam_R.s")] = with(Rnonparam, c(delta, delta.s, R.s))
  sim_res[r, c("gs_nonparam_var_R.s", "gs_nonparam_normci_lb_R.s", "gs_nonparam_normci_ub_R.s", 
               "gs_nonparam_quantci_lb_R.s", "gs_nonparam_quantci_ub_R.s")] = with(Rnonparam, c(R.s.var, conf.int.normal.R.s, conf.int.quantile.R.s))
  
  ## Estimate R with parametric approach (gold standard) ###
  Rparam = R.s.estimate(sone=s1, 
                        szero=s0, 
                        yone = y1, 
                        yzero = y0, 
                        type = "model", 
                        conf.int = TRUE)
  sim_res[r, c("gs_param_delta", "gs_param_delta.s", "gs_param_R.s")] = with(Rparam, c(delta, delta.s, R.s))
  sim_res[r, c("gs_param_var_R.s", "gs_param_normci_lb_R.s", "gs_param_normci_ub_R.s", 
               "gs_param_quantci_lb_R.s", "gs_param_quantci_ub_R.s")] = with(Rparam, c(R.s.var, conf.int.normal.R.s, conf.int.quantile.R.s))
  
  # Simulate non-missingness indicators ###################
  ## Under MCAR, everybody has 35% missingness probability
  m1 = rbinom(n1, 1, 0.65) 
  m0 = rbinom(n0, 1, 0.65)
  s0[m0==0] = NA ### make them missing
  s1[m1==0] = NA ### make them missing
  
  ##########################################################
  #Estimates with incomplete data ##########################
  ##########################################################
  ## Estimate R with nonparametric approach (complete case)
  Rnonparam_miss = R.s.estimate(sone = s1[m1==1], 
                                szero = s0[m0==1], 
                                yone = y1[m1==1], 
                                yzero = y0[m0==1], 
                                type = "robust", 
                                conf.int = TRUE)
  sim_res[r, c("cc_nonparam_delta", "cc_nonparam_delta.s", "cc_nonparam_R.s")] = with(Rnonparam_miss, c(delta, delta.s, R.s))
  sim_res[r, c("cc_nonparam_var_R.s", "cc_nonparam_normci_lb_R.s", "cc_nonparam_normci_ub_R.s", 
               "cc_nonparam_quantci_lb_R.s", "cc_nonparam_quantci_ub_R.s")] = with(Rnonparam_miss, c(R.s.var, conf.int.normal.R.s, conf.int.quantile.R.s))
  
  ## Estimate R with parametric approach (complete case)
  Rparam_miss = R.s.estimate(sone = s1[m1==1], 
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
  ipw_fit = glm(formula = m ~ 1, 
                family = "binomial")
  w1 = rep(x = 1 / (1 + exp(-ipw_fit$coefficients)), 
           times = n1)
  w0 = rep(x = 1 / (1 + exp(-ipw_fit$coefficients)), 
           times = n0)
  
  ## Estimate R with nonparametric approach (IPW)
  Rnonparam_miss_ipw = R.s.miss(sone = s1, 
                                szero = s0, 
                                yone = y1, 
                                yzero = y0, 
                                type = "robust",
                                wone = w1, 
                                wzero = w0, 
                                ipw_formula = m ~ 1, 
                                conf.int = TRUE)
  sim_res[r, c("ipw_nonparam_delta", "ipw_nonparam_delta.s", "ipw_nonparam_R.s")] = with(Rnonparam_miss_ipw, c(delta, delta.s, R.s))
  sim_res[r, c("ipw_nonparam_var_R.s", "ipw_nonparam_normci_lb_R.s", "ipw_nonparam_normci_ub_R.s", 
               "ipw_nonparam_quantci_lb_R.s", "ipw_nonparam_quantci_ub_R.s")] = with(Rnonparam_miss_ipw, c(R.s.var, conf.int.normal.R.s, conf.int.quantile.R.s))
  
  ## Estimate R with parametric approach (IPW)
  Rparam_miss_ipw = R.s.miss(sone = s1, 
                             szero = s0, 
                             yone = y1, 
                             yzero = y0, 
                             type = "model",
                             wone = w1, 
                             wzero = w0, 
                             ipw_formula = m ~ 1, 
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
    write.csv(paste0("sett1_mcar/sett1_mcar_seed", sim_seed, ".csv"), 
              row.names = FALSE)
}
