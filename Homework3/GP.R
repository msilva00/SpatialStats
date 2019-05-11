#geoR::matern(1,2,3)
source("mcmc.R")

matern <- function(d, phi, nu) {
  #' @param
  #' d:     distance
  #' phi:   range
  #' nu:    smoothness
  #' @export
  u <- ifelse(d > 0, d / phi, .Machine$double.eps)
  logR <- -((nu-1) * log(2) + lgamma(nu)) + 
    nu * log(u) + log(besselK(u, nu))
  exp(logR)
}


gp <- function(y, X, s, 
               stepSigPsi=diag(3),
               a_tau=2, b_tau=1, 
               a_sig=2, b_sig=1, # gam2 = tau2 / sig2
               a_phi=0, b_phi=2, 
               a_z=0, b_z=3,
               B=2000, burn=1000, print_every=0) {
  
  n <- NROW(X)
  k <- NCOL(X)
  psi_dim <- 3 # gam2, phi, z
  Xt <- t(X)
  stopifnot(n == length(y) && n == NROW(s))
  stopifnot(a_z >= 0 && b_z > a_z)
  
  nu_z <- function(z) {
    stopifnot(z >= a_z && z <= b_z) # this should never happen
    floor(z) + 1/2
  }
  
  D <- as.matrix(dist(s))
  I_n <- diag(n)
  
  compute_ll <- TRUE
  
  update <- function(state) {
    out <- state
    
    # update (gam2, phi, z), beta, tau2
    ll <- function(trans_psi) {
      gam2 <- exp(trans_psi[1])
      phi <- inv_logit(trans_psi[2], a_phi, b_phi)
      z <- inv_logit(trans_psi[3], a_z, b_z)
      nu <- nu_z(z)
      
      R <- matern(D, phi, nu)
      V <- I_n + R / gam2
      
      ### See: https://stackoverflow.com/questions/39568978/how-to-calculate-variance-of-least-squares-estimator-using-qr-decomposition-in-r
      U <- chol(V)                         # U = L', and LL' = V
      z <- backsolve(U, y, transpose=TRUE) # z = L^(-1)*y
      G <- backsolve(U, X, transpose=TRUE) # G = L^(-1)*X
      #beta_hat <- qr.solve(G, z)
      qr.G <- qr(G)
      b.qty <- qr.qty(qr.G, z)
      beta_hat <- backsolve(qr.G$qr, b.qty)
      sse <- sum(qr.resid(qr.G, z) ^ 2) # S_psi: (z-G*bhat)'(z-G*bhat)
      Sig_hat <- chol2inv(qr.G$qr) # (G'G)^(-1)
      
      val <- -.5 * log_det(V) + .5 * log_det(Sig_hat) - 
        ((n-k) / 2 + a_tau) * log(sse / 2 + b_tau)
      
      list(val=val, beta_hat=beta_hat, sse=sse, Sig_hat=Sig_hat)
    }
    
    lp <- function(trans_psi) {
      # this doesn't seem to work
      #lp_log_invgamma(trans_psi[1], a_gam, b_gam) + # gam2
      # gam2: transforming parameters
      # let alpha2 = tau2, gam2 = tau2 / sig2
      # find joint of alpha2, gam2 and integrate out alpha2
      # the following is the density of log_gam2
      # which is the log density of gam2 + log(gam2)
      # = log(p_gam2(gam2)) + log(gam2)
      gam2 <- exp(log_gam2 <- trans_psi[1])
      log_gam2 + (a_sig-1) * log(gam2) - 
        (a_sig + a_tau) * log(b_sig*gam2 + b_tau) + 
        lp_logit_unif(trans_psi[2]) +  # phi
        lp_logit_unif(trans_psi[3])    # z
    }
    
    trans_curr_psi <- c(log(out$gam2), 
                        logit(out$phi, a_phi, b_phi),
                        logit(out$z, a_z, b_z))
    
    ### Metropolis
    # Doesn't help to precompute cholesky because psi_dim is small
    cand <- mvrnorm(trans_curr_psi, stepSigPsi)
    ll_cand_all <- ll(cand)
    ll_cand <- ll_cand_all$val
    ll_curr_all <- if(compute_ll) {
      compute_ll <<- FALSE
      ll_old_all <<- ll(trans_curr_psi) 
      ll_old_all
    } else ll_old_all
    
    ll_curr <- ll_curr_all$val
    
    ### Compute Acceptance Ratio
    if (ll_cand + lp(cand) - 
        ll_curr - lp(trans_curr_psi) > log(runif(1))) {
      
      out$gam2 <- exp(cand[1])
      out$phi <- inv_logit(cand[2], a_phi, b_phi)
      out$z <- inv_logit(cand[3], a_z, b_z)
      
      # update tau2
      out$tau2 <- 1 / rgamma(1, (n-k) / 2 + a_tau,
                             ll_cand_all$sse / 2 + b_tau)
      
      # update beta
      out$beta <- mvrnorm(ll_cand_all$beta_hat, 
                          ll_cand_all$Sig_hat * out$tau2)
      
      ll_old_all <<- ll_cand_all
    }
    
    
    out
  }
  
  init <- list(beta=double(k),
               tau2=b_tau, 
               gam2=b_sig, 
               phi= (a_phi + b_phi) / 2,
               z=(a_z + b_z) / 2)
  
  gibbs_out <- gibbs(init, update, B, burn, print_every)
  
  out <- matrix(NA, psi_dim + k + 1 + 2, B) 
  
  out[1:k, ] <- sapply(gibbs_out, function(x) x$beta)
  out[k+1,] <- sapply(gibbs_out, function(x) x$gam2)
  out[k+2,] <- sapply(gibbs_out, function(x) x$phi)
  out[k+3,] <- sapply(gibbs_out, function(x) x$z)
  out[k+4,] <- sapply(gibbs_out, function(x) x$tau2)
  out[k+5,] <- sapply(gibbs_out, function(x) nu_z(x$z))
  out[k+6,] <- sapply(gibbs_out, function(x) x$tau2 / x$gam2)
  
  rownames(out) <- c(paste0('beta',1:k), 
                     'gam2', 'phi', 'z', 'tau2', 'nu', 'sig2')
  t(out)
}

# post = matrix(beta, tau2, sig2, phi, nu)
gp.predict <- function(y, X, s, X_new, s_new, post) {
  n <- nrow(X)
  k <- ncol(X)
  m <- nrow(X_new)
  I_n <- diag(n)
  I_m <- diag(m)
  I_all <- diag(n+m)
  X_all <- rbind(X_new, X)
  
  stopifnot(n == length(y) && n == nrow(s))
  stopifnot(m == nrow(s_new))
  stopifnot(k == ncol(X_new))
  
  D_all <- as.matrix(dist(rbind(s_new, s)))
  
  pred <- function(state) {
    beta <- state[1:k]
    tau2 <- state['tau2']
    sig2 <- state['sig2']
    phi  <- state['phi']
    nu   <- state['nu']
    Xb <- X_all %*% beta
    
    R_all <- matern(D_all, phi, nu)
    V_all <- tau2 * I_all + sig2 * R_all
    V_new <- V_all[1:m, 1:m]
    V_old <- V_all[-c(1:m), -c(1:m)]
    V_new_old <- V_all[1:m, -c(1:m)]
    
    S <- t(solve(V_old, t(V_new_old))) # Not faster
    #S <- V_new_old %*% solve(V_old)
    EY <- Xb[1:m] + S %*% (y - Xb[-c(1:m)])
    VY <- V_new - S %*% t(V_new_old)
    
    mvrnorm(EY, VY)
  }
  
  apply(post, 1, pred)
}





