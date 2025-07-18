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