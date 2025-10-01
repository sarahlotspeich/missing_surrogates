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
  return(sum(Kern.FUN(zz,zi.one,bw=bw)*y1*(weight))/sum(Kern.FUN(zz,zi.one,bw=bw)*weight))
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
    delta.s = sum((weight.0)*mu.1.s0)/n0.all-sum((weight.0)*yzero)/n0.all
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
    delta.var = as.numeric(se_res$var_delta), 
    delta.s.var = as.numeric(se_res$var_delta.s), 
    R.s.var = as.numeric(se_res$var_R.s), 
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
      w = 1 / predict(object = ipw_fit, 
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
    delta = sum(wone[mone == 1] * yone[mone == 1]) / none - 
      sum(wzero[mzero == 1] * yzero[mzero == 1]) / nzero
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
    # Create long version of *observed* data (with missingness)
    long_dat = data.frame(Y = c(yzero, yone), ## primary outcome 
                          Z = c(rep(c(0, 1), times = c(nzero, none))), ## treatment group
                          S = c(szero, sone), ## surrogate marker 
                          R = as.numeric(!is.na(c(szero, sone))), ## (non-)missingness indicator
                          W = c(wzero, wone)) ## weight 
    long_dat = long_dat[order(long_dat$S, decreasing = TRUE), ] ## order to put non-missing first
    
    ## Model Y ~ S * Z (saturated)
    modYgivSZ = lm(formula = Y ~ Z * S, 
                   data = long_dat, 
                   weights = W)
    
    ## Separate coefficients
    beta0 = as.numeric(modYgivSZ$coefficients["(Intercept)"])
    beta1 = as.numeric(modYgivSZ$coefficients["Z"])
    beta2 = as.numeric(modYgivSZ$coefficients["S"])
    beta3 = as.numeric(modYgivSZ$coefficients["Z:S"])
    
    ## Model S ~ Z
    modSgivZ = lm(formula = S ~ Z, 
                  data = long_dat, 
                  weights = W)
    
    ## Separate coefficients
    alpha0 = as.numeric(modSgivZ$coefficients["(Intercept)"])
    alpha1 = alpha0 + as.numeric(modSgivZ$coefficients["Z"])
    
    ## Construct percent treatment effect explained
    delta = beta1 + (beta2 + beta3) * alpha1 - beta2 * alpha0
    delta_S = beta1 + beta3 * alpha0
    R_S = 1 - delta_S / delta
    
    ## Return 
    list(delta = as.numeric(delta), 
         delta.s = as.numeric(delta_S), 
         R.s = as.numeric(R_S))
    
    # R.s.estimate(sone = wone[mone == 1] * sone[mone == 1], 
    #              szero = wzero[mzero == 1] * szero[mzero == 1], 
    #              yone = yone[mone == 1], 
    #              yzero = yzero[mzero == 1], 
    #              type = type, 
    #              conf.int = FALSE)
  }
}

# Sieve maximum likelihood estimator -- Only type = "model"
R.s.miss_model_smle <- function(sone, szero, yone, yzero,
                                nonparam = TRUE,
                                conv_res = NULL,
                                max_it = 1e4,
                                tol = 1e-3,
                                full_output = FALSE) {
  # sizes
  N0 <- length(szero)
  N1 <- length(sone)
  n  <- N0 + N1

  # Build long observed data in natural order; create ID column
  long_dat <- data.frame(
    Y = c(yzero, yone),
    Z = c(rep(0, N0), rep(1, N1)),
    S = c(szero, sone),
    R = as.numeric(!is.na(c(szero, sone))),
    ID = seq_len(n)
  )

  # Separate observed and missing
  long_nonmiss <- long_dat[long_dat$R == 1, , drop = FALSE]
  long_miss    <- long_dat[long_dat$R == 0, , drop = FALSE]

  # Support of observed S values by arm
  S_z0 <- sort(unique(na.omit(long_dat$S[long_dat$Z == 0])))
  S_z1 <- sort(unique(na.omit(long_dat$S[long_dat$Z == 1])))
  m_z0 <- length(S_z0)
  m_z1 <- length(S_z1)

  # Initialize parameters (make prev_p numeric vectors)
  prev_beta  <- rep(0, 4)   # intercept, Z, S, Z:S
  prev_sigma <- 0.1
  if (nonparam) {
    prev_p_z0 <- if (m_z0 > 0) as.numeric(rep(1 / m_z0, m_z0)) else numeric(0)
    prev_p_z1 <- if (m_z1 > 0) as.numeric(rep(1 / m_z1, m_z1)) else numeric(0)
  } else {
    prev_gamma <- c(0, 0)
    prev_eta   <- 0.1
  }

  # Use conv_res if provided
  if (!is.null(conv_res)) {
    if (!is.null(conv_res$betas)) prev_beta  <- conv_res$betas
    if (!is.null(conv_res$sigma)) prev_sigma <- conv_res$sigma
    if (nonparam && !is.null(conv_res$p0) && !is.null(conv_res$p1)) {
      prev_p_z0 <- as.numeric(conv_res$p0)
      prev_p_z1 <- as.numeric(conv_res$p1)
    } else if (!nonparam && !is.null(conv_res$gamma) && !is.null(conv_res$eta)) {
      prev_gamma <- conv_res$gamma
      prev_eta   <- conv_res$eta
    }
  }

  # helper for sigma
  calc_sigma <- function(lmobj) {
    if (is.null(lmobj)) return(NA_real_)
    res <- lmobj$residuals
    rdf <- lmobj$df.residual
    if (is.null(rdf) || rdf <= 0) return(NA_real_)
    sqrt(sum(res^2, na.rm = TRUE) / rdf)
  }

  converged <- FALSE
  it <- 1L

  while (!converged && it <= max_it) {
    # ---------- E-step ----------
    cd_nonmiss <- long_nonmiss
    # expand each missing subject to all candidate S in that arm
    expand_rows <- function(df, Svals) {
      if (nrow(df) == 0 || length(Svals) == 0) return(df[FALSE, , drop = FALSE])
      out_list <- lapply(seq_len(nrow(df)), function(i) {
        rowi <- df[i, , drop = FALSE]
        newrows <- rowi[rep(1, length(Svals)), , drop = FALSE]
        newrows$S <- Svals
        newrows
      })
      do.call(rbind, out_list)
    }
    cd_miss_z0 <- expand_rows(long_miss[long_miss$Z == 0, , drop = FALSE], S_z0)
    cd_miss_z1 <- expand_rows(long_miss[long_miss$Z == 1, , drop = FALSE], S_z1)
    cd_miss    <- if (nrow(cd_miss_z0) + nrow(cd_miss_z1) == 0) cd_miss_z0 else rbind(cd_miss_z0, cd_miss_z1)
    cd         <- if (nrow(cd_miss) == 0) cd_nonmiss else rbind(cd_nonmiss, cd_miss)

    # build phi_aug: 1's for observed rows then posterior probs for expanded rows
    phi_aug <- rep(1, nrow(cd_nonmiss))
    if (nrow(cd_miss) > 0) {
      # P(Y | S, Z) on expanded rows
      mu_beta_miss <- prev_beta[1] + prev_beta[2] * cd_miss$Z +
                      prev_beta[3] * cd_miss$S + prev_beta[4] * cd_miss$S * cd_miss$Z
      pY_miss <- dnorm(cd_miss$Y, mean = mu_beta_miss, sd = prev_sigma)

      # P(S | Z) on expanded rows -> use numeric prev_p vectors robustly
      if (nonparam) {
        pS_miss <- numeric(nrow(cd_miss))
        # fill block for Z==0 expanded rows
        if (nrow(cd_miss_z0) > 0) {
          n_miss_z0_subj <- nrow(long_miss[long_miss$Z == 0, , drop = FALSE])
          # repeat prev_p_z0 for each missing subject in z0 (order matches expand_rows)
          pS_miss[seq_len(nrow(cd_miss_z0))] <- rep(as.numeric(prev_p_z0), times = n_miss_z0_subj)
        }
        if (nrow(cd_miss_z1) > 0) {
          start1 <- if (nrow(cd_miss_z0) > 0) nrow(cd_miss_z0) + 1 else 1
          n_miss_z1_subj <- nrow(long_miss[long_miss$Z == 1, , drop = FALSE])
          if (n_miss_z1_subj > 0) {
            pS_miss[start1:(start1 + nrow(cd_miss_z1) - 1)] <- rep(as.numeric(prev_p_z1), times = n_miss_z1_subj)
          }
        }
      } else {
        mu_gamma_miss <- prev_gamma[1] + prev_gamma[2] * cd_miss$Z
        pS_miss <- dnorm(cd_miss$S, mean = mu_gamma_miss, sd = prev_eta)
      }

      phi_num <- pY_miss * pS_miss
      # group-sum phi_num by original subject ID
      denom_byID <- rowsum(phi_num, group = cd_miss$ID, reorder = FALSE)
      denom_rep <- as.numeric(denom_byID[as.character(cd_miss$ID)])
      denom_rep[is.na(denom_rep)] <- 1
      denom_rep[denom_rep == 0] <- 1
      phi_expanded <- phi_num / denom_rep
      phi_aug <- c(phi_aug, phi_expanded)
    }

    if (length(phi_aug) != nrow(cd)) stop("phi_aug length mismatch with cd rows")

    # ---------- M-step ----------
    new_fit <- tryCatch(lm(Y ~ Z * S, data = cd, weights = phi_aug), error = function(e) NULL)
    if (is.null(new_fit)) {
      if (full_output) {
        return(list(delta = NA, delta.s = NA, R.s = NA, betas = rep(NA,4), sigma = NA, p0 = NA, p1 = NA, alphas = c(NA,NA)))
      } else {
        return(list(delta = NA, delta.s = NA, R.s = NA))
      }
    }
    new_beta <- coef(new_fit)
    # ensure length 4 and order by "(Intercept)","Z","S","Z:S"
    all_names <- c("(Intercept)", "Z", "S", "Z:S")
    tmp <- setNames(rep(0, 4), all_names)
    tmp[names(new_beta)] <- new_beta
    new_beta <- as.numeric(tmp)
    new_sigma <- calc_sigma(new_fit)

    # Update p_{kz} robustly: sum phi_aug over rows with Z=z and S == s_k
    if (nonparam) {
      # Z == 0
      idx_z0 <- which(cd$Z == 0)
      if (length(S_z0) == 0) {
        new_p_z0 <- numeric(0)
      } else {
        sums_z0 <- sapply(S_z0, function(sk) {
          sum(phi_aug[idx_z0][cd$S[idx_z0] == sk], na.rm = TRUE)
        })
        if (sum(sums_z0) == 0) sums_z0 <- rep(1/length(sums_z0), length(sums_z0))
        new_p_z0 <- as.numeric(sums_z0 / sum(sums_z0))
      }

      # Z == 1
      idx_z1 <- which(cd$Z == 1)
      if (length(S_z1) == 0) {
        new_p_z1 <- numeric(0)
      } else {
        sums_z1 <- sapply(S_z1, function(sk) {
          sum(phi_aug[idx_z1][cd$S[idx_z1] == sk], na.rm = TRUE)
        })
        if (sum(sums_z1) == 0) sums_z1 <- rep(1/length(sums_z1), length(sums_z1))
        new_p_z1 <- as.numeric(sums_z1 / sum(sums_z1))
      }
    } else {
      sfit <- tryCatch(lm(S ~ Z, data = cd, weights = phi_aug), error = function(e) NULL)
      if (is.null(sfit)) {
        new_gamma <- prev_gamma; new_eta <- prev_eta
      } else {
        new_gamma <- coef(sfit); if (length(new_gamma) < 2) new_gamma <- c(new_gamma, 0)[1:2]
        new_eta <- calc_sigma(sfit)
      }
    }

    # ---------- convergence ----------
    beta_conv  <- all(!is.na(prev_beta)) && all(!is.na(new_beta)) && all(abs(prev_beta - new_beta) < tol)
    sigma_conv <- !is.na(prev_sigma) && !is.na(new_sigma) && (abs(prev_sigma - new_sigma) < tol)
    if (nonparam) {
      # check shape and non-NA
      ok0 <- length(prev_p_z0) == length(new_p_z0) && length(new_p_z0) > 0
      ok1 <- length(prev_p_z1) == length(new_p_z1) && length(new_p_z1) > 0
      if (ok0 && ok1) {
        p_conv <- all(abs(prev_p_z0 - new_p_z0) < tol) && all(abs(prev_p_z1 - new_p_z1) < tol)
      } else {
        # if no missing subjects, treat p_conv as TRUE (nothing to update)
        if (nrow(long_miss) == 0) p_conv <- TRUE else p_conv <- FALSE
      }
      gamma_conv <- TRUE; eta_conv <- TRUE
    } else {
      p_conv <- TRUE
      gamma_conv <- all(!is.na(prev_gamma)) && all(!is.na(new_gamma)) && all(abs(prev_gamma - new_gamma) < tol)
      eta_conv   <- !is.na(prev_eta) && !is.na(new_eta) && (abs(prev_eta - new_eta) < tol)
    }

    converged <- isTRUE(all(c(beta_conv, sigma_conv, p_conv, gamma_conv, eta_conv)))

    # update prev for next iter
    prev_beta  <- new_beta; prev_sigma <- new_sigma
    if (nonparam) { prev_p_z0 <- new_p_z0; prev_p_z1 <- new_p_z1 } else { prev_gamma <- new_gamma; prev_eta <- new_eta }
    it <- it + 1L
  } # end EM loop

  if (!converged) {
    if (full_output) {
      return(list(delta = NA, delta.s = NA, R.s = NA, betas = rep(NA,4), sigma = NA, p0 = NA, p1 = NA, alphas = c(NA,NA)))
    } else {
      return(list(delta = NA, delta.s = NA, R.s = NA))
    }
  }

  # compute alpha0, alpha1 and estimands
  beta0 <- prev_beta[1]; beta1 <- prev_beta[2]; beta2 <- prev_beta[3]; beta3 <- prev_beta[4]
  if (nonparam) {
    alpha0 <- if (length(S_z0) > 0 && length(prev_p_z0) > 0) sum(S_z0 * as.numeric(prev_p_z0)) else NA_real_
    alpha1 <- if (length(S_z1) > 0 && length(prev_p_z1) > 0) sum(S_z1 * as.numeric(prev_p_z1)) else NA_real_
  } else {
    alpha0 <- prev_gamma[1]; alpha1 <- alpha0 + prev_gamma[2]
  }

  delta   <- beta1 + (beta2 + beta3) * alpha1 - beta2 * alpha0
  delta.s <- beta1 + beta3 * alpha0
  R_S     <- 1 - delta.s / delta


  if (full_output) {
    return(list(delta = delta, delta.s = delta.s, R.s = R_S,
                betas = prev_beta, sigma = prev_sigma, p0 = prev_p_z0, p1 = prev_p_z1, alphas = c(alpha0, alpha1)))
  } else {
    return(list(delta = delta, delta.s = delta.s, R.s = R_S))
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
