# Libraries and functions
library(missSurrogate)

# Reproducibility 
## Random seed to be used for each simulation setting
args = commandArgs(TRUE)
## When running on the cluster, give each array a unique seed by using the array ID
sim_seed = as.integer(args) + 20
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
  smle_param_delta = NA, smle_param_delta.s = NA, smle_param_R.s = NA, smle_param_var_R.s = NA, smle_param_normci_lb_R.s = NA, smle_param_normci_ub_R.s = NA, smle_param_quantci_lb_R.s = NA, smle_param_quantci_ub_R.s = NA) 
for (r in 1:REPS) {
  # Generate data 
  data = gen.data(n1=n1, n0=n0) 
  
  # Define vectors for outcomes/surrogates in untreated/treated 
  s1 = data$s1
  y1 = data$y1
  s0 = data$s0
  y0 = data$y0
  
  #########################################################
  # Simulate non-missingness indicators ###################
  ## Under MCAR, everybody has 35% missingness probability
  m1 = rbinom(n1, 1, 0.65) 
  m0 = rbinom(n0, 1, 0.65)
  s0[m0==0] = NA ### make them missing
  s1[m1==0] = NA ### make them missing
  
  ##########################################################
  #Estimates with incomplete data ##########################
  ##########################################################  
  ## Estimate R with parametric approach (SMLE)
  Rparam_miss_smle = R.s.miss(sone = s1, 
                              szero = s0,
                              yone = y1,
                              yzero = y0, 
                              type = "model", 
                              conf.int = TRUE, 
                              orig.smle = TRUE) 
  sim_res[r, c("smle_param_delta", "smle_param_delta.s", "smle_param_R.s")] = with(Rparam_miss_smle, c(delta, delta.s, R.s))
  sim_res[r, c("smle_param_var_R.s", "smle_param_normci_lb_R.s", "smle_param_normci_ub_R.s", 
               "smle_param_quantci_lb_R.s", "smle_param_quantci_ub_R.s")] = with(Rparam_miss_smle, c(R.s.var, conf.int.normal.R.s, conf.int.quantile.R.s))
  
  ## Save 
  sim_res |> 
    write.csv(paste0("sett1_mcar_smle/sett1_mcar_seed", sim_seed, ".csv"), 
              row.names = FALSE)
}
