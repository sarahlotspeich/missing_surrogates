perturb_resample = function(d, sone, szero, yone, yzero, type = "robust", smle = TRUE, ipw,
                            max_it = 1E4, tol = 1E-3, conv_res = NULL, ipw_formula) {
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

