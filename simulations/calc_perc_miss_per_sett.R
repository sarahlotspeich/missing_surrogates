# Be reproducible queens
set.seed(0) 

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
REPS = 1000
## Set sample sizes 
n1 = 1000
n0 = 1000

# Generate data 
perc_miss_sett1 = perc_miss_sett2 = perc_miss_sett3 = vector(length = REPS)
for (r in 1:REPS) {
  ## Generate data 
  data = gen.data(n1=n1, n0=n0)   
  
  # Define vectors for outcomes/surrogates in untreated/treated 
  s1 = data$s1
  y1 = data$y1
  s0 = data$s0
  y0 = data$y0
  
  ## Sett 1: Under MCAR, everybody has 30% missingness probability
  m1 = rbinom(n1, 1, 0.65) 
  m0 = rbinom(n0, 1, 0.65)
  perc_miss_sett1[r] = mean(c(m0, m1) == 0) ### save % missing values 
  
  ## Sett 2: Under MAR, probability of missingness depends on Y continuously (logistic regression)
  m1 = rbinom(n = n1, size = 1, prob = 1 / (1 + exp(- 0.015 * y1)))
  m0 = rbinom(n = n0, size = 1, prob = 1 / (1 + exp(- 0.015 * y0)))
  perc_miss_sett2[r] = mean(c(m0, m1) == 0) ### save % missing values 
  
  ## Sett 3: Under MAR, probability of missingness depends on Z (logistic regression)
  m1 = rbinom(n = n1, size = 1, prob = 1 / (1 + exp(- 0.6)))
  m0 = rbinom(n = n0, size = 1, prob = 1 / (1 + exp(- 0.4)))
  perc_miss_sett3[r] = mean(c(m0, m1) == 0) ### save % missing values 
}

# Summarize % missing values per setting 
summary(perc_miss_sett1) ## 32-39%
summary(perc_miss_sett2) ## 35-41%
summary(perc_miss_sett3) ## 34-41%