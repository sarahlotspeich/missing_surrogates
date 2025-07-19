########################################
### Kernel function ####
########################################
# it just takes an x and pops back the kernel
kf=function(x, h){return(dnorm(x/h) / h)}

######################################################################
### Kernel function but setting it up to deal with matrix/vector ####
######################################################################
#but same idea as kf above
VTM<-function(vc, dm){
  #takes vc and makes it the repeated row of a matrix, repeats it dm times
  matrix(vc, ncol=length(vc), nrow = dm, byrow = T)
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
          distance = abs(szero - szero[o])
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
          distance = abs(szero - szero[o])
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

R.s.miss = function(sone, szero, yone, yzero, wone = NULL, wzero = NULL, 
                    type = "robust", max_it = 1E4, tol = 1E-3, 
                    conf.int = FALSE, boot = FALSE, ipw_formula = NULL) {
  # Estimate parameters
  est_res = R.s.miss.estimate(sone = sone, szero = szero, ## surrogates outcomes
                              yone = yone, yzero = yzero, ## primary outcomes
                              wone = wone, wzero = wzero, ## weights (optional)
                              type = type, max_it = max_it, tol = tol) ## other arguments 
  
  # Estimate standard errors 
  ## Perturbation resampling
  if (conf.int & !boot) {
    se_res = R.s.miss.se.pert(num_pert = 500, conv_res = est_res,
                              sone = sone, szero = szero, yone = yone, yzero = yzero, 
                              type = type, max_it = max_it, tol = tol, ipw_formula = ipw_formula)
  } else if (conf.int) { # Bootstrap resampling
    se_res = R.s.miss.se.boot(num_boot = 500, conv_res = est_res, 
                              sone = sone, szero = szero, yone = yone, yzero = yzero, 
                              type = type, max_it = max_it, tol = tol, ipw_formula = ipw_formula)
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

R.s.miss.estimate = function(weight_perturb = NULL, sone, szero, yone, yzero, wone = NULL, wzero = NULL, conv_res = NULL, type, max_it, tol, ipw_formula = NULL) {
  ## Define sample sizes 
  none = length(yone) 
  nzero = length(yzero)
  
  ## Define non-missingness indicators for the two treatment groups
  mone = as.numeric(!is.na(sone))
  mzero = as.numeric(!is.na(szero))
  
  ### If weight_perturb provided, multiply surrogates and outcomes by it 
  if (!is.null(weight_perturb)) {
    szero = szero * weight_perturb[1:nzero]
    yzero = yzero * weight_perturb[1:nzero]
    sone = sone * weight_perturb[-c(1:nzero)]
    yone = yone * weight_perturb[-c(1:nzero)]
  }
  
  # Define TRUE/FALSE use IPW based on non-null weights supplied
  ipw = (!is.null(wone) & !is.null(wzero)) || ## either weights were supplied
    !is.null(ipw_formula) ## or the model formula was
  
  # Estimate parameters and SEs 
  if (ipw) {
    ## If ipw_formula is not NULL, re-calculate weights 
    if (!is.null(ipw_formula) & is.null(wone) & is.null(wzero)) {
      ### Define vectors of variables for model (all of them just in case)
      m = c(mone, mzero)
      z = rep(x = c(1, 0), times = c(length(sone), length(szero)))
      s = c(sone, szero)
      y = c(yone, yzero)
      
      ### Fit the IPW model 
      ipw_fit = glm(formula = as.formula(ipw_formula), 
                    data = data.frame(m, z, s, y),
                    family = "binomial")
      
      ### Get estimated weights (probabilities of being non-missing) for each observation 
      w = predict(object = ipw_fit, 
                  type = "response")
      
      ### Split weights into vectors for treatment/control
      wone = w[1:length(sone)]
      wzero = w[-c(1:length(sone))]
    }
    
    ## Using IPW to handle missing data
    est_res = R.s.miss_ipw(sone = sone, szero = szero, ### surrogates outcomes
                           yone = yone, yzero = yzero, ### primary outcomes
                           wone = wone, wzero = wzero, ### weights (required)
                           type = type) ### type of PTE estimator
  } else if (type == "model") { # Wang & Taylor's approach 
    est_res = R.s.miss_model_smle(sone = sone, 
                                  szero = szero, 
                                  yone = yone, 
                                  yzero = yzero, 
                                  nonparam = TRUE, 
                                  conv_res = conv_res, 
                                  max_it = max_it, 
                                  tol = tol)  
  }
  
  ### Return point estimates 
  return(est_res)
}

# Inverse probability weighting estimators -- type = "model" or "roboust"
R.s.miss_ipw = function(sone, szero, yone, yzero, wone, wzero, type) {
  ## Define non-missingness indicators for the two treatment groups
  mone = as.numeric(!is.na(sone))
  mzero = as.numeric(!is.na(szero))
  
  ## Define sample sizes 
  none = length(yone) 
  nzero = length(yzero)
  
  ## Calculate PTE 
  if (type == "robust") {
    delta = sum((1 / wone[mone == 1]) * yone[mone == 1]) / none - 
      sum((1 / wzero[mzero == 1]) * yzero[mzero == 1]) / nzero
    delta_S = delta.s.single(sone = sone[mone == 1], 
                             szero = szero[mzero == 1], 
                             yone = yone[mone == 1], 
                             yzero = yzero[mzero == 1], 
                             weight.1 = wone[mone == 1], 
                             weight.0 = wzero[mzero == 1],
                             n0.all = nzero)
    R_S = 1 - delta_S / delta
    
    ## Return 
    list(delta = delta, 
         delta.s = delta_S, 
         R.s = R_S)
  } else if (type == "model") {
    R.s.estimate(sone = wone[mone == 1] * sone[mone == 1], 
                 szero = wzero[mzero == 1] * szero[mzero == 1], 
                 yone = yone[mone == 1], 
                 yzero = yzero[mzero == 1], 
                 type = type, 
                 conf.int = FALSE)
  }
}

# Sieve maximum likelihood estimator -- Only type = "model"
R.s.miss_model_smle = function(sone, szero, yone, yzero, nonparam, conv_res, max_it = 1E4, tol = 1E-3, full_output = FALSE) {
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
    if (full_output) {
      list(delta = delta, 
           delta.s = delta_S, 
           R.s = R_S, 
           betas = new_beta, 
           sigma = new_sigma, 
           p0 = new_p_z0, 
           p1 = new_p_z1,
           alphas = c(alpha0, alpha1))
    } else {
      list(delta = delta, 
           delta.s = delta_S, 
           R.s = R_S)
    }
    
  } else {
    ## Return 
    if (full_output) {
      list(delta = NA, 
           delta.s = NA, 
           R.s = NA, 
           betas = rep(NA, length(new_beta)), 
           sigma = NA, 
           p0 = NA, 
           p1 = NA,
           alphas = rep(NA, 2))
    } else {
      list(delta = NA, 
           delta.s = NA, 
           R.s = NA)
    }
  }
}

# Perturbation resampling SEs 
R.s.miss.se.pert = function(num_pert, conv_res, sone, szero, yone, yzero, 
                            ipw_formula, type, max_it, tol) {
  # Create (n0 + n1) x num_pert matrix of perturbations
  weight_perturb = matrix(rexp(n = num_pert * (length(yone) + length(yzero)), rate = 1), 
                          ncol = num_pert)
  
  # Apply the estimation functions with it 
  pert_quant = do.call(what = rbind, 
                       args = apply(X = weight_perturb, 
                                    MARGIN = 2, ## apply across columns 
                                    FUN = R.s.miss.estimate, 
                                    sone = sone, 
                                    szero = szero, 
                                    yone = yone, 
                                    yzero = yzero, 
                                    ipw_formula = ipw_formula, 
                                    type = type, 
                                    max_it = max_it, 
                                    tol = tol)
                       )
  
  ## Create separate vectors for perturbed quantities
  pert_delta = unlist(pert_quant[, "delta"])
  pert_delta.s = unlist(pert_quant[, "delta.s"])
  pert_R.s = unlist(pert_quant[, "R.s"])
  
  # Calculate two types of 95% confidence intervals
  ## Normal approximation 
  ### Variance estimates for each quantity
  var_delta = var(pert_delta)
  var_delta.s = var(pert_delta.s)
  var_R.s = var(pert_R.s)
  
  ### Used to compute Wald-type confidence intervals
  norm_ci_delta = conv_res$delta + c(-1.96, 1.96) * sqrt(var_delta)
  norm_ci_delta.s = conv_res$delta.s + c(-1.96, 1.96) * sqrt(var_delta.s)
  norm_ci_R.s = conv_res$R.s + c(-1.96, 1.96) * sqrt(var_R.s)
  
  ## Quantile-based 
  quant_ci_delta = as.vector(quantile(x = pert_delta, probs = c(0.025, 0.975)))
  quant_ci_delta.s = as.vector(quantile(x = pert_delta.s, probs = c(0.025, 0.975)))
  quant_ci_R.s = as.vector(quantile(x = pert_R.s, probs = c(0.025, 0.975)))
  
  ## Return 
  return(
    list(
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

# Bootstrap resampling SEs
boot_R.s.miss <- function(data, indices, type, ipw_formula, conv_res, max_it, tol) {
  ## Resample data with replacement (but stratified on treatment)
  d <- data[indices, ]  
  
  ### Split into vectors 
  which_zzero = which(d$z == 0)
  which_zone = which(d$z == 1)
  Sd_zero = d$s[which_zzero]
  Sd_one = d$s[which_zone]
  Yd_zero = d$y[which_zzero]
  Yd_one = d$y[which_zone]
  
  res_d = R.s.miss.estimate(sone = Sd_one, 
                            szero = Sd_zero, 
                            yone = Yd_one, 
                            yzero = Yd_zero, 
                            wone = NULL, 
                            wzero = NULL, 
                            conv_res = conv_res, 
                            type = type, 
                            max_it = max_it, 
                            tol = tol, 
                            ipw_formula = ipw_formula)
  
  return(with(res_d, c(delta, delta.s, R.s)))
}

R.s.miss.se.boot = function(num_boot, conv_res, sone, szero, yone, yzero, 
                            type, ipw, smle, max_it, tol, ipw_formula) {
  # Build long dataset: (y, z, s, w, m)
  long_dat = data.frame(y = c(yone, yzero), 
                        z = rep(x = c(1, 0), 
                                times = c(length(yone), length(yzero))), 
                        s = c(sone, szero), 
                        m = c(as.numeric(!is.na(sone)), as.numeric(!is.na(szero))))
  
  # Initialize empty dataframe to hold bootstrapped quantities
  boot_quant = data.frame(delta = rep(NA, num_boot), 
                          delta.s = rep(NA, num_boot), 
                          R.s = rep(NA, num_boot))
  
  # Bootstrap resample from long dataset and fit estimator to it 
  for (b in 1:num_boot) {
    ## Indices for which rows to resample, preserving the treatment/control split 
    ind_b = c(sample(x = which(long_dat$z == 0), 
                     size = length(which(long_dat$z == 0)), 
                     replace = TRUE), 
              sample(x = which(long_dat$z == 1), 
                     size = length(which(long_dat$z == 1)), 
                     replace = TRUE))
    
    boot_quant[b, ] = boot_R.s.miss(data = long_dat, 
                                    indices = ind_b, 
                                    type = type,
                                    ipw_formula = ipw_formula, 
                                    conv_res = NULL, 
                                    max_it = max_it, 
                                    tol = tol)
  }
  
  # boot_quant = boot::boot(long_dat, 
  #                         boot_R.s.miss, 
  #                         R = num_boot, 
  #                         strata = long_dat$z, 
  #                         type = type,
  #                         ipw_formula = ipw_formula, 
  #                         conv_res = conv_res, 
  #                         max_it = max_it, 
  #                         tol = tol)
  
  # Calculate two types of 95% confidence intervals
  ## Normal approximation 
  ### Variance estimates for each quantity
  var_all = apply(X = boot_quant, 
                  MARGIN = 2, 
                  FUN = var)
  var_delta = var_all[1]
  var_delta.s = var_all[2]
  var_R.s = var_all[3]
  
  ### Used to compute Wald-type confidence intervals
  norm_ci_delta = conv_res$delta + c(-1.96, 1.96) * sqrt(var_delta)
  norm_ci_delta.s = conv_res$delta.s + c(-1.96, 1.96) * sqrt(var_delta.s)
  norm_ci_R.s = conv_res$R.s + c(-1.96, 1.96) * sqrt(var_R.s)
  
  ## Quantile-based 
  quant_ci_all = apply(X = boot_quant, 
                       MARGIN = 2, 
                       FUN = function(x) quantile(x = x, probs = c(0.025, 0.975)))
  quant_ci_delta = as.vector(quant_ci_all[, 1])
  quant_ci_delta.s = as.vector(quant_ci_all[, 2])
  quant_ci_R.s = as.vector(quant_ci_all[, 3])
  
  ## Return 
  return(
    list(
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
