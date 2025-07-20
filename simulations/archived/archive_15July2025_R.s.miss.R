########################################
### Kernel function ####
########################################
# it just takes an x and pops back the kernel
kf=function(x, h){return(dnorm(x/h)/h)}

######################################################################
### Kernel function but setting it up to deal with matrix/vector ####
######################################################################
#but same idea as kf above
VTM<-function(vc, dm){
  #takes vc and makes it the repeated row of a matrix, repeats it dm times
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}

Kern.FUN <- function(zz,zi,bw) ## returns an (n x nz) matrix ##
{ 
  out = (VTM(zz,length(zi))- zi)/bw
  return(dnorm(out)/bw)
  
}

########################################
### My functions to estimate delta.s ####
########################################

#I make this so I can use an apply function in the next function
pred.smooth <-function(zz,zi.one, bw=NULL,y1, weight=NULL) {
  if(is.null(bw)) { bw = bw.nrd(zz)/((length(zz))^(.10))}
  if(is.null(weight)) {weight = rep(1, length(y1))}
  return(sum(Kern.FUN(zz,zi.one,bw=bw)*y1*(1/weight))/sum(Kern.FUN(zz,zi.one,bw=bw)*1/weight))
}


delta.s.single = function(sone,szero,yone,yzero, h.select = NULL, weight.1 = NULL, weight.0=NULL, n0.all=NULL) {
  #we can talk about the bandwidth later, this default should work ok
  if(is.null(h.select)) {h.select = bw.nrd(sone)*(length(sone)^(-0.25)) 
  }
  if(is.null(weight.1) & is.null(weight.0)){
    mu.1.s0 = sapply(szero,pred.smooth,zz=sone, bw=h.select, y1=yone)
    if(sum(is.na(mu.1.s0))>0){
      print(paste("Note: ", sum(is.na(mu.1.s0)), " values extrapolated."))
      c.mat = cbind(szero, mu.1.s0)
      for(o in 1:length(mu.1.s0)) {
        if(is.na(mu.1.s0[o])){
          distance = abs(s0.new - s0.new[o])
          c.temp = cbind(c.mat, distance)
          c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where mean is not na
          new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
          mu.1.s0[o] = new.est[1]   #in case there are multiple matches
        }
      }}
    delta.s = mean(mu.1.s0) - mean(yzero)
  }
  if(!is.null(weight.1) & !is.null(weight.0)){
    mu.1.s0 = sapply(szero,pred.smooth,zz=sone, bw=h.select, y1=yone, weight=weight.1)
    if(sum(is.na(mu.1.s0))>0){
      print(paste("Note: ", sum(is.na(mu.1.s0)), " values extrapolated."))
      c.mat = cbind(szero, mu.1.s0)
      for(o in 1:length(mu.1.s0)) {
        if(is.na(mu.1.s0[o])){
          distance = abs(s0.new - s0.new[o])
          c.temp = cbind(c.mat, distance)
          c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where mean is not na
          new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
          mu.1.s0[o] = new.est[1]   #in case there are multiple matches
        }
      }}
    delta.s = sum((1/weight.0)*mu.1.s0)/n0.all-sum((1/weight.0)*yzero)/n0.all
  }
  return(delta.s)
}
##########################################################
#expit function#
##########################################################
expit = function(x){
  return((exp(x))/(1+exp(x)))
}

R.s.miss = function(sone, szero, yone, yzero, type = "robust", smle = TRUE,
                    wone = NULL, wzero = NULL, max_it = 1E4, tol = 1E-3, 
                    conf.int = FALSE, boot = FALSE, ipw_formula = m ~ y * z) {
  # Estimate parameters
  est_res = R.s.miss.estimate(sone = sone, szero = szero, 
                              yone = yone, yzero = yzero, 
                              type = type, smle = smle,
                              wone = wone, wzero = wzero, 
                              max_it = max_it, tol = tol)
  
  # Estimate standard errors 
  ## Perturbation resampling
  if (conf.int & !boot) {
    se_res = R.s.miss.se.pert(num_pert = 500, conv_res = est_res,
                              sone = sone, szero = szero, yone = yone, yzero = yzero, 
                              type = type, smle = smle, wone = wone, wzero = wzero, 
                              max_it = max_it, tol = tol, ipw_formula = ipw_formula)
  } else if (conf.int) { # Bootstrap resampling
    # boot_quant = boot(data = data, statistic = boot_fn, R = 1000)
  } else {
    ### Otherwise, just return point estimates (without SEs/CIs)
    return(est_res)
  }
  
  ### Build final return list (based on R.s.estimate syntax)
  res_list = list(
    delta = as.numeric(est_res$delta), 
    delta.s = as.numeric(est_res$delta.s), 
    R.s = as.numeric(est_res$R.s), 
    delta.var = se_res$var_delta, 
    delta.s.var = se_res$var_delta.s, 
    R.s.var = se_res$var_R.s, 
    conf.int.normal.delta = se_res$norm_ci_delta, 
    conf.int.quantile.delta = se_res$quant_ci_delta, 
    conf.int.normal.delta.s = se_res$norm_ci_delta.s, 
    conf.int.quantile.delta.s = se_res$quant_ci_delta.s, 
    conf.int.normal.R.s = se_res$norm_ci_R.s, 
    conf.int.quantile.R.s = se_res$quant_ci_R.s, 
    conf.int.fieller.R.s = c(NA, NA)
  )
  
  ## And return it 
  return(res_list)
}

R.s.miss.estimate = function(sone, szero, yone, yzero, type, smle,
                             wone, wzero, max_it, tol) {
  # Define TRUE/FALSE use IPW based on non-null weights supplied
  ipw = !is.null(wone) & !is.null(wzero)
  
  # Estimate parameters
  if (type == "model") { # Wang & Taylor's approach 
    ## using IPW to handle missing data
    if (ipw) {
      est_res = R.s.miss_model_ipw(sone = sone, szero = szero, 
                                   yone = yone, yzero = yzero,
                                   wone = wone, wzero = wzero)      
    } else { ## using MLE or SMLE to handle missing data
      est_res = R.s.miss_model_smle(sone = sones, szero = szero, 
                                    yone = yzero, yzero = yzero, 
                                    nonparam = smle, conv_res = NULL, 
                                    max_it = max_it, tol = tol)  
    }
  } else if (type == "robust") {
    est_res = R.s.miss_robust_ipw(sone = sone, szero = szero, 
                                  yone = yone, yzero = yzero, 
                                  wone = wone, wzero = wzero)
  }
  
  ### Return point estimates 
  return(est_res)
}


# Function to apply to each bootstrap sample
# This function returns the slope coefficient from lm(y ~ x)
boot_R.s.miss <- function(data, indices) {
  d <- data[indices, ]  # resample data
  
  ### Simulate perturbations from Expo(1)
  Vd_one = rexp(n = length(yone), rate = 1) 
  Vd_zero = rexp(n = length(yzero), rate = 1) 
  
  ### Multiply Y and S by them 
  Yd_one = yone * Vd_one 
  Sd_one = sone * Vd_one
  Yd_zero = yzero * Vd_zero 
  Sd_zero = szero * Vd_zero
  
  ### Re-calculate weights for IPW approaches
  if (ipw) {
    #### Define non-missingness indicators for the two treatment groups
    mone = as.numeric(!is.na(Sd_one))
    mzero = as.numeric(!is.na(Sd_zero))
    
    #### Define vectors of variables for model (all of them just in case)
    m = c(mone, mzero)
    z = rep(x = c(1, 0), times = c(length(Sd_one), length(Sd_zero)))
    s = c(Sd_one, Sd_zero)
    y = c(Yd_one, Yd_zero)
    
    #### Fit the IPW model 
    ipw_fit = glm(formula = as.formula(ipw_formula), 
                  family = "binomial")
    
    #### Get estimated weights for each patient 
    w = predict(object = ipw_fit, 
                type = "response")
    
    #### Split weights into vectors for treatment/control
    wd_one = w[1:length(Sd_one)]
    wd_zero = w[-c(1:length(Sd_one))]
  }
  ### Re-estimate parameters using perturbed data
  if (type == "model") { # Wang & Taylor's approach 
    #### Using IPW to handle missing data
    if (ipw) {
      res_d = R.s.miss_model_ipw(Sd_one, Sd_zero, Yd_one, Yd_zero, wd_one, wd_zero)      
    } else { #### using SMLE to handle missing data
      res_d = R.s.miss_model_smle(Sd_one, Sd_zero, Yd_one, Yd_zero, nonparam = smle, conv_res = conv_res, max_it, tol)  
    }
  } else if (type == "robust") {
    res_d = R.s.miss_robust_ipw(Sd_one, Sd_zero, Yd_one, Yd_zero, wd_one, wd_zero)
  }
  return(cbind(it = d, with(res_d, data.frame(delta, delta.s, R.s))))
  
  
  model <- lm(y ~ x, data = d)
  return(coef(model))  # return slope
}

R.s.miss.se.boot = function(num_boot, conv_res, sone, szero, yone, yzero, 
                            type, smle, wone, wzero, max_it, tol, ipw_formula) {
  
}

R.s.miss.se.pert = function(num_pert, conv_res, sone, szero, yone, yzero, 
                            type, smle, wone, wzero, max_it, tol, ipw_formula) {
  # Define TRUE/FALSE use IPW based on non-null weights supplied
  ipw = !is.null(wone) & !is.null(wzero)
  
  # Loop over perturbation resampling D = 500 times 
  pert_quant = do.call(what = rbind, 
                       args = #parallel::mclapply(X = 1:num_pert, 
                              sapply(X = 1:num_pert, 
                                    FUN = perturb_resample, 
                                    sone = sone, szero = szero, 
                                    yone = yone, yzero = yzero, 
                                    type = type, smle = smle, 
                                    ipw = ipw, max_it = max_it, 
                                    tol = tol, conv_res = conv_res, 
                                    ipw_formula = ipw_formula, 
                                    simplify = FALSE))
  
  # Calculate two types of 95% confidence intervals
  ## Normal approximation 
  ### Variance estimates for each quantity
  var_delta = var(pert_quant$delta)
  var_delta.s = var(pert_quant$delta.s)
  var_R.s = var(pert_quant$R.s)
  
  ### Used to compute Wald-type confidence intervals
  norm_ci_delta = conv_res$delta + c(-1.96, 1.96) * sqrt(var_delta)
  norm_ci_delta.s = conv_res$delta.s + c(-1.96, 1.96) * sqrt(var_delta.s)
  norm_ci_R.s = conv_res$R.s + c(-1.96, 1.96) * sqrt(var_R.s)
  
  ## Quantile-based 
  quant_ci_delta = as.vector(quantile(x = pert_quant$delta, probs = c(0.025, 0.975)))
  quant_ci_delta.s = as.vector(quantile(x = pert_quant$delta.s, probs = c(0.025, 0.975)))
  quant_ci_R.s = as.vector(quantile(x = pert_quant$R.s, probs = c(0.025, 0.975)))
  
  ## Return 
  return(
    res_list = list(
      var_delta = var_delta, 
      var_delta.s = var_delta.s, 
      var_R.s = var_R.s, 
      norm_ci_delta = norm_ci_delta, 
      quant_ci_delta = quant_ci_delta, 
      norm_ci_delta.s = norm_ci_delta.s, 
      quant_ci_delta.s = quant_ci_delta.s, 
      norm_ci_R.s = norm_ci_R.s, 
      quant_ci_R.s = quant_ci_R.s
    )
  )
}

perturb_resample = function(d, sone, szero, yone, yzero, type = "robust", smle = TRUE, ipw,
                            max_it = 1E4, tol = 1E-3, conv_res = NULL, ipw_formula = m ~ y * z) {
  ### Simulate perturbations from Expo(1)
  Vd_one = rexp(n = length(yone), rate = 1) 
  Vd_zero = rexp(n = length(yzero), rate = 1) 
  
  ### Multiply Y and S by them 
  Yd_one = yone * Vd_one 
  Sd_one = sone * Vd_one
  Yd_zero = yzero * Vd_zero 
  Sd_zero = szero * Vd_zero
  
  ### Re-calculate weights for IPW approaches
  if (ipw) {
    #### Define non-missingness indicators for the two treatment groups
    mone = as.numeric(!is.na(Sd_one))
    mzero = as.numeric(!is.na(Sd_zero))
    
    #### Define vectors of variables for model (all of them just in case)
    m = c(mone, mzero)
    z = rep(x = c(1, 0), times = c(length(Sd_one), length(Sd_zero)))
    s = c(Sd_one, Sd_zero)
    y = c(Yd_one, Yd_zero)
    
    #### Fit the IPW model 
    ipw_fit = glm(formula = as.formula(ipw_formula),
                  data = data.frame(m, z, s, y), 
                  family = "binomial")
    
    #### Get estimated weights for each patient 
    w = predict(object = ipw_fit, 
                type = "response")
    
    #### Split weights into vectors for treatment/control
    wd_one = w[1:length(Sd_one)]
    wd_zero = w[-c(1:length(Sd_one))]
  }
  ### Re-estimate parameters using perturbed data
  if (type == "model") { # Wang & Taylor's approach 
    #### Using IPW to handle missing data
    if (ipw) {
      res_d = R.s.miss_model_ipw(Sd_one, Sd_zero, Yd_one, Yd_zero, wd_one, wd_zero)      
    } else { #### using SMLE to handle missing data
      res_d = R.s.miss_model_smle(Sd_one, Sd_zero, Yd_one, Yd_zero, nonparam = smle, conv_res = conv_res, max_it, tol)  
    }
  } else if (type == "robust") {
    res_d = R.s.miss_robust_ipw(Sd_one, Sd_zero, Yd_one, Yd_zero, wd_one, wd_zero)
  }
  return(cbind(it = d, with(res_d, data.frame(delta, delta.s, R.s))))
}

R.s.miss_robust_ipw = function(sone, szero, yone, yzero, wone, wzero) {
  # Define non-missingness indicators for the two treatment groups
  mone = as.numeric(!is.na(sone))
  mzero = as.numeric(!is.na(szero))
  
  ## Define sample sizes 
  none = length(yone) 
  nzero = length(yzero)
  
  ## Calculate PTE 
  delta = sum((1 / wone[mone==1]) * yone[mone==1]) / none - sum((1 / wzero[mzero==1]) * yzero[mzero==1]) / nzero
  delta_S = delta.s.single(sone = s1[mone==1], 
                           szero = s0[mzero==1], 
                           yone = yone[mone==1], 
                           yzero = yzero[mzero==1], 
                           weight.1 = wone[mone==1], 
                           weight.0 = wzero[mzero==1],
                           n0.all = nzero)
  R_S = 1 - delta_S / delta
  
  ## Return 
  list(delta = delta, 
       delta.s = delta_S, 
       R.s = R_S)
}

R.s.miss_model_ipw = function(sone, szero, yone, yzero, wone, wzero) {
  # Build long dataset: (Y, Z, S, W)
  long_dat = data.frame(Y = c(yone, yzero), 
                        Z = rep(x = c(1, 0), 
                                times = c(length(yone), length(yzero))), 
                        S = c(sone, szero), 
                        W = c(wone, wzero))
  
  # Get IPW estimates for linear regression of Y ~ Z + S + X x Z
  fit_beta = lm(formula = Y ~ Z * S,
                data = long_dat, 
                weights = W)
  
  ## Extract coefficient estimates 
  beta0 = fit_beta$coefficients[1]
  beta1 = fit_beta$coefficients[2]
  beta2 = fit_beta$coefficients[3]
  beta3 = fit_beta$coefficients[4]
  
  # Get IPW estimates for linear regression of S ~ Z
  fit_alpha = lm(formula = S ~ Z, 
                 data = long_dat, 
                 weights = W)
  
  ## Extract coefficient estimates 
  alpha0 = fit_alpha$coefficients[1] ### E(S|Z=0)
  alpha1 = alpha0 + fit_alpha$coefficients[2] ### E(S|Z=1)
  
  ## Construct percent treatment effect explained
  delta = beta1 + (beta2 + beta3) * alpha1 - beta2 * alpha0
  delta_S = beta1 + beta3 * alpha0
  R_S = 1 - delta_S / delta
  
  ## Return 
  list(delta = delta, 
       delta.s = delta_S, 
       R.s = R_S, 
       betas = fit_beta$coefficients, 
       alphas = c(alpha0, alpha1))
}

R.s.miss_model_smle = function(sone, szero, yone, yzero, nonparam, conv_res, max_it = 1E4, tol = 1E-3) {
  # Save useful constants 
  N0 = length(szero) ## number in control group
  N1 = length(sone) ## number in treatment group
  n = N0 + N1 ## total sample size
  
  # Create long version of *observed* data (with missingness)
  long_dat = data.frame(Y = c(yzero, yone), ## primary outcome 
                        Z = c(rep(c(0, 1), times = c(N0, N1))), ## treatment group
                        S = c(szero, sone), ## surrogate marker 
                        R = as.numeric(!is.na(c(szero, sone)))) ## (non-)missingness indicator
  long_dat = long_dat[order(long_dat$S, decreasing = TRUE), ] ## order to put non-missing first
  long_dat$ID = 1:n ## create an ID
  
  ## Split the observed data by treatment group (with missingness)
  long_dat_z0 = long_dat[long_dat$Z == 0, ] ## control group
  long_dat_z1 = long_dat[long_dat$Z == 1, ] ## treatment group
  
  # Create initial values 
  ## Linear regression of Y ~ Z + S + X x Z
  prev_beta = beta0 = rep(0, 4)
  prev_sigma = sigma0 = 0.1
  # cc_fit = lm(formula = Y ~ Z * S, 
  #             data = long_dat)
  # prev_beta = beta0 = as.numeric(cc_fit$coefficients)
  # prev_sigma = sigma0 = sigma(cc_fit)
  
  # Count and save unique non-missing values of the surrogate
  ## Control group
  S_z0 = unique(long_dat_z0$S[long_dat_z0$R == 1]) ## unique values of non-missing surrogates (ordered descendingly)
  m_z0 = length(S_z0) ## number of unique non-missing surrogates
  n0 = sum(long_dat_z0$R == 1) ## number of patients with non-missing surrogates
  
  ## Treatment group
  S_z1 = unique(long_dat_z1$S[long_dat_z1$R == 1]) ## unique values of non-missing surrogates (ordered descendingly)
  m_z1 = length(S_z1) ## number of unique non-missing surrogates
  n1 = sum(long_dat_z1$R == 1) ## number of patients with non-missing surrogates
  
  # Conditional distribution of S given Z
  if (nonparam) {
    ## Empirical probabilities of S | Z = 0
    prev_p_z0 = p0_z0 = matrix(data = 1 / m_z0, 
                               nrow = m_z0, 
                               ncol = 1)
    
    ## Empirical probabilities of S | Z = 1
    prev_p_z1 = p0_z1 = matrix(data = 1 / m_z1, 
                               nrow = m_z1, 
                               ncol = 1)
  } else {
    ## Linear regression of S ~ Z 
    prev_gamma = gamma0 = rep(0, 2)
    prev_eta = eta0 = 0.1
  }
  
  # If converged values supplied, use them as initials 
  if (!is.null(conv_res)) {
    ## Linear regression of Y ~ Z + S + X x Z
    prev_beta = beta0 = conv_res$betas
    prev_sigma = sigma0 = conv_res$sigma
    # Conditional distribution of S given Z
    if (nonparam) {
      ## Empirical probabilities of S | Z = 0
      prev_p_z0 = p0_z0 = matrix(data = conv_res$p0, 
                                 ncol = 1)
      
      ## Empirical probabilities of S | Z = 1
      prev_p_z1 = p0_z1 = matrix(data = conv_res$p1, 
                                 ncol = 1)
    } 
  } 
  
  # Create even longer version of *complete* data (without missingness)
  cd_nonmiss = long_dat[1:(n0 + n1), ] ## patients with non-missing surrogate markers, both treatment groups
  cd_nonmiss_z0 = cd_nonmiss[cd_nonmiss$Z == 0, ] ## separate control group patients
  cd_nonmiss_z1 = cd_nonmiss[cd_nonmiss$Z == 1, ] ## separate control group patients
  
  ## Complete data for Z = 0
  cd_miss_z0 = long_dat_z0[rep(x = (n0 + 1):N0, each = m_z0), ] ### create m0 copies of each patient with missing surrogate
  cd_miss_z0$S = rep(x = S_z0, times = (N0 - n0)) ## try out different surrogate values
  
  ## Complete data for Z = 1
  cd_miss_z1 = long_dat_z1[rep(x = (n1 + 1):N1, each = m_z1), ] ### create m1 copies of each patient with missing surrogate
  cd_miss_z1$S = rep(x = S_z1, times = (N1 - n1)) ### try out different surrogate values
  
  ## Combined complete dataset
  cd_miss = rbind(cd_miss_z0, cd_miss_z1) ### only those with missing surrogates
  cd = rbind(cd_nonmiss, cd_miss) ### all patients 
  
  # EM Algorithm 
  converged = FALSE ## initialize as unconverged
  it = 1 ## initialize iteration counter 
  while (!converged & it <= max_it) {
    ## E step 
    ### Update the phi_ki = P(S=sk|Zi) for patients w/ missing surrogate -------
    ### Outcome model: P(Y|S,Z) ------------------------------------------------
    #### mu = beta0 + beta1X + beta2Z + ...
    mu_beta = prev_beta[1] + prev_beta[2] * cd_miss$Z + 
      prev_beta[3] * cd_miss$S + prev_beta[4] * cd_miss$S * cd_miss$Z
    #### Calculate P(Y|S,Z) from normal distribution ---------------------------
    pYgivSZ = dnorm(x = cd_miss$Y,
                    mean = mu_beta, 
                    sd = prev_sigma)
    ############################################################################
    ### Conditional distribution of surrogate given treatment: P(S|Z) ----------
    if (nonparam) {
      #### Nonparametric
      pSgivZ = c(prev_p_z0[rep(x = 1:m_z0, times = (N0 - n0))], #### take from p_k0 if Z = 0 
                 prev_p_z1[rep(x = 1:m_z1, times = (N1 - n1))]) #### take from p_k1 if Z = 1  
    } else {
      #### mu = gamma0 + gammaZ
      mu_gamma = prev_gamma[1] + prev_gamma[2] * cd_miss$Z 
      #### Calculate P(Y|S,Z) from normal distribution
      pSgivZ = dnorm(x = cd_miss$S,
                     mean = mu_gamma, 
                     sd = prev_eta)
    }
    ############################################################################
    ## Estimate conditional expectations ---------------------------------------
    ### Update numerator -------------------------------------------------------
    #### P(Y|S,Z)P(S|Z) --------------------------------------------------------
    phi_num = pYgivSZ * pSgivZ 
    ### Update denominator -----------------------------------------------------
    #### Sum over P(Y|S,Z)P(S|Z) per patient -----------------------------------
    phi_denom = rowsum(x = phi_num, 
                       group = cd_miss$ID, 
                       reorder = FALSE)
    #### Avoid NaN resulting from dividing by 0 --------------------------------
    phi_denom[phi_denom == 0] = 1
    ### Divide them to get phi = E{I(S=s)|Y,Z} ---------------------------------
    phi = phi_num / 
      c(rep(x = phi_denom[1:(N0 - n0)], each = m_z0), #### repeat denominator for Z = 0
        rep(x = phi_denom[-c(1:(N0 - n0))], each = m_z1)) #### repeat denominator for Z = 1
    #### Add indicators for non-missing rows -----------------------------------
    phi_aug = c(rep(x = 1, times = (n0 + n1)), phi)
    
    ## M step 
    ### Re-fit the linear regression model 
    new_fit = lm(formula = Y ~ Z * S, 
                 data = cd, 
                 weights = phi_aug)
    new_beta = as.numeric(new_fit$coefficients)
    new_sigma = sigma(new_fit)
    
    ### Re-estimate distribution of S | Z 
    if (nonparam) {
      #### Re-estimate the empirical probabilities 
      sum_phi_z0 = rowsum(x = phi_aug[cd$Z == 0], #### Sum over i = 1, ..., n of phi-hats... 
                          group = cd$S[cd$Z == 0], #### for each k = 1, ..., m ...
                          reorder = FALSE) #### and keep them in the original order
      lambda_z0 = sum(sum_phi_z0) #### Sum over k = 1, ..., m for constraint
      new_p_z0 = sum_phi_z0 / lambda_z0 #### updated probabilities
      
      sum_phi_z1 = rowsum(x = phi_aug[cd$Z == 1], #### Sum over i = 1, ..., n of phi-hats... 
                          group = cd$S[cd$Z == 1], #### for each k = 1, ..., m ...
                          reorder = FALSE) #### and keep them in the original order
      lambda_z1 = sum(sum_phi_z1) #### Sum over k = 1, ..., m for constraint
      new_p_z1 = sum_phi_z1 / lambda_z1 #### updated probabilities
    } else {
      #### Re-fit the regression model
      new_fit = lm(formula = S ~ Z, 
                   data = cd, 
                   weights = phi_aug)
      new_gamma = as.numeric(new_fit$coefficients)
      new_eta = sigma(new_fit)
    }
    ## Check for convergence 
    beta_conv = !any(abs(prev_beta - new_beta) > tol)
    sigma_conv = !(abs(prev_sigma - new_sigma) > tol)
    if (nonparam) {
      p_conv = c(!any(abs(prev_p_z0 - new_p_z0) > tol), 
                 !any(abs(prev_p_z1 - new_p_z1) > tol))
      gamma_conv = TRUE
      eta_conv = TRUE
    } else {
      p_conv = TRUE
      gamma_conv = !any(abs(prev_gamma - new_gamma) > tol)
      eta_conv = !any(abs(prev_eta - new_eta) > tol)
    }
    if (mean(c(beta_conv, sigma_conv, p_conv, gamma_conv, eta_conv)) == 1) {
      converged = TRUE ### Success! 
    }
    
    ## If not converged, prepare move on to next iteration
    it = it + 1
    prev_beta = new_beta 
    prev_sigma = new_sigma
    if (nonparam) {
      prev_p_z0 = new_p_z0
      prev_p_z1 = new_p_z1  
    } else {
      prev_gamma = new_gamma
      prev_eta = new_eta
    }
  }
  
  # Return percent of treatment effect explained 
  if (converged) {
    ## Define linear regression coefficients
    beta0 = new_beta[1]
    beta1 = new_beta[2]
    beta2 = new_beta[3]
    beta3 = new_beta[4]
    
    ## Define conditional mean coefficients
    if (nonparam) {
      alpha0 = sum(S_z0 * new_p_z0) ### E(S|Z=0)
      ### Complete case mean: mean(long_dat$S[long_dat$Z == 0], na.rm = TRUE) 
      alpha1 = sum(S_z1 * new_p_z1) ### E(S|Z=1)
      ### Complete case mean: mean(long_dat$S[long_dat$Z == 1], na.rm = TRUE)   
    } else {
      alpha0 = new_gamma[1] ### E(S|Z=0)
      alpha1 = alpha0 + new_gamma[2] ### E(S|Z=1)
    }
    
    ## Construct percent treatment effect explained
    delta = beta1 + (beta2 + beta3) * alpha1 - beta2 * alpha0
    delta_S = beta1 + beta3 * alpha0
    R_S = 1 - delta_S / delta
    
    ## Return 
    list(delta = delta, 
         delta.s = delta_S, 
         R.s = R_S, 
         betas = new_beta, 
         sigma = new_sigma, 
         p0 = new_p_z0, 
         p1 = new_p_z1,
         alphas = c(alpha0, alpha1))
  } else {
    ## Return 
    list(delta = NA, 
         delta.s = NA, 
         R.s = NA, 
         betas = rep(NA, length(new_beta)), 
         sigma = NA, 
         p0 = NA, 
         p1 = NA,
         alphas = rep(NA, 2))
  }
}
