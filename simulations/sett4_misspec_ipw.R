# Libraries and functions
library(Rsurrogate)

# Source IPW/SMLE function from GitHub
devtools::source_url("https://raw.githubusercontent.com/sarahlotspeich/missing_surrogates/refs/heads/main/R.s.miss.R")

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
REPS = 50
## Set sample sizes 
n1 = 1000
n0 = 1000
## Initialize empty dataframe for results
sim_res = data.frame(
  r = 1:REPS, 
  gs_nonparam = NA, gs_param = NA, 
  cc_nonparam = NA, cc_param = NA, 
  ipw_nonparam_Yonly = NA, ipw_param_Yonly = NA, 
  ipw_nonparam_Zonly = NA, ipw_param_Zonly = NA, 
  ipw_nonparam_YZ = NA, ipw_param_YZ = NA) 
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
                           conf.int = FALSE)
  sim_res[r, "gs_nonparam"] = Rnonparam$R.s
  
  ## Estimate R with parametric approach (gold standard) ###
  Rparam = R.s.estimate(sone=s1, 
                        szero=s0, 
                        yone = y1, 
                        yzero = y0, 
                        type = "model", 
                        conf.int = FALSE)
  sim_res[r, "gs_param"] = Rparam$R.s
  
  # Simulate non-missingness indicators ###################
  ## Under MAR, probability of missingness depends on Y continuously (logistic regression)
  m1 = rbinom(n = n1, size = 1, prob = 1 / (1 + exp(- 0.030 * y1)))
  m0 = rbinom(n = n0, size = 1, prob = 1 / (1 + exp(- 0.015 * y0)))
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
                                conf.int = FALSE)
  sim_res[r, "cc_nonparam"] = Rnonparam_miss$R.s
  
  ## Estimate R with parametric approach (complete case)
  Rparam_miss = R.s.estimate(sone = s1[m1==1], 
                             szero = s0[m0==1], 
                             yone = y1[m1==1], 
                             yzero = y0[m0==1], 
                             type = "model", 
                             conf.int = FALSE)
  sim_res[r, "cc_param"] = Rparam_miss$R.s
  
  ## Calculate weights for IPW approaches
  ### IPW Model 1: M ~ Y 
  m = c(m1, m0)
  y = c(y1, y0)
  ipw_fit = glm(formula = m ~ y, 
                family = "binomial")
  p1 = 1 / (1 + exp(-(ipw_fit$coefficients[1] + ipw_fit$coefficients[2] * y1)))
  w1 = 1 / p1
  p0 = 1 / (1 + exp(-(ipw_fit$coefficients[1] + ipw_fit$coefficients[2] * y0)))
  w0 = 1 / p0  

  ### Estimate R with nonparametric approach (IPW)
  Rnonparam_miss_ipw = R.s.miss(sone = s1, 
                                szero = s0, 
                                yone = y1, 
                                yzero = y0, 
                                type = "robust",
                                wone = w1, 
                                wzero = w0, 
                                conf.int = FALSE)
  sim_res[r, "ipw_nonparam_Yonly"] = Rnonparam_miss_ipw$R.s
  
  ### Estimate R with parametric approach (IPW)
  Rparam_miss_ipw = R.s.miss(sone = s1, 
                             szero = s0, 
                             yone = y1, 
                             yzero = y0, 
                             type = "model",
                             wone = w1, 
                             wzero = w0, 
                             conf.int = FALSE)
  sim_res[r, "ipw_param_Yonly"] = Rparam_miss_ipw$R.s
  
  ### IPW Model 2: M ~ Z
  m = c(m1, m0)
  z = rep(c(1, 0), each = 1000)
  ipw_fit = glm(formula = m ~ z, 
                family = "binomial")
  p1 = 1 / (1 + exp(-(ipw_fit$coefficients[1] + ipw_fit$coefficients[2])))
  w1 = 1 / rep(p1, 1000)
  p0 = 1 / (1 + exp(-(ipw_fit$coefficients[1])))
  w0 = 1 / rep(p0, 1000)  
  
  ### Estimate R with nonparametric approach (IPW)
  Rnonparam_miss_ipw = R.s.miss(sone = s1, 
                                szero = s0, 
                                yone = y1, 
                                yzero = y0, 
                                type = "robust",
                                wone = w1, 
                                wzero = w0, 
                                conf.int = FALSE)
  sim_res[r, "ipw_nonparam_Zonly"] = Rnonparam_miss_ipw$R.s
  
  ### Estimate R with parametric approach (IPW)
  Rparam_miss_ipw = R.s.miss(sone = s1, 
                             szero = s0, 
                             yone = y1, 
                             yzero = y0, 
                             type = "model",
                             wone = w1, 
                             wzero = w0, 
                             conf.int = FALSE)
  sim_res[r, "ipw_param_Zonly"] = Rparam_miss_ipw$R.s
  
  ### IPW Model 2: M ~ Z
  ipw_fit = glm(formula = m ~ y * z, 
                family = "binomial")
  p1 = 1 / (1 + exp(-(ipw_fit$coefficients[1] + ipw_fit$coefficients[2] * y1 + ipw_fit$coefficients[3] + ipw_fit$coefficients[4] * y1)))
  w1 = 1 / p1
  p0 = 1 / (1 + exp(-(ipw_fit$coefficients[1] + ipw_fit$coefficients[2] * y0)))
  w0 = 1 / p0  
  
  ### Estimate R with nonparametric approach (IPW)
  Rnonparam_miss_ipw = R.s.miss(sone = s1, 
                                szero = s0, 
                                yone = y1, 
                                yzero = y0, 
                                type = "robust",
                                wone = w1, 
                                wzero = w0, 
                                conf.int = FALSE)
  sim_res[r, "ipw_nonparam_YZ"] = Rnonparam_miss_ipw$R.s
  
  ### Estimate R with parametric approach (IPW)
  Rparam_miss_ipw = R.s.miss(sone = s1, 
                             szero = s0, 
                             yone = y1, 
                             yzero = y0, 
                             type = "model",
                             wone = w1, 
                             wzero = w0, 
                             conf.int = FALSE)
  sim_res[r, "ipw_param_YZ"] = Rparam_miss_ipw$R.s
  
  ## Save 
  sim_res |> 
    write.csv(paste0("sett4_misspec_ipw/sett4_misspec_ipw_seed", sim_seed, ".csv"), 
              row.names = FALSE)
}
