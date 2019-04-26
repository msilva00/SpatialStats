library(rootSolve)

combine <- function(x,y) {
  stopifnot(length(x) == length(y))
  n <- length(x)
  z <- rep(NA,n)
  for (i in 1:n) {
    z[2*i - 1] <- x[i]
    z[2*i] <- y[i]
  }
  z
}

# Karhunen-Loeve Representation for exponential covariance fn
KL_exp <- function(s1, s2, J=5, L=1, sig2=1, phi=1, 
                   inf=1E3, n=1E6, eps=.Machine$double.eps) {
  
  w1 <- uniroot.all(function(w) tan(w*L) - phi/w ,
                    interval=c(0,L),  n=n)
  w2 <- uniroot.all(function(w) tan(w*L) + w/phi ,
                    interval=c(0,L), n=n)
  
  w1 <- head(sort(w1[which(abs(tan(w1/L)) < inf)]), J)
  w2 <- head(sort(w2[which(abs(tan(w2/L)) < inf)]), J)
  #print(w2)
  
  if (length(w1) != length(w2)) {
    cat("length(w1): ", length(w1), "\n")
    cat("length(w2): ", length(w2), "\n")
  }
  stopifnot(length(w1) == length(w2))
  
  #lam1 <- (2 * phi) / (w1^2 + phi^2)
  #lam2 <- (2 * phi) / (w2^2 + phi^2)
  lam1 <- (2 * phi) / (w1 + phi)^2
  lam2 <- (2 * phi) / (w2 + phi)^2
  lam <- combine(lam1, lam2)
  
  # returns the vectors of length J
  psi1 <- function(s) cos(w1*s) / sqrt(L+sin(2*w1*L) / (2*w1+eps))
  psi2 <- function(s) sin(w2*s) / sqrt(L-sin(2*w2*L) / (2*w2+eps))
  
  psi <- function(s) combine(psi1(s), psi2(s))
  
  sum( lam * psi(s1) * Conj(psi(s2)) )
}


eigen_approx <- function(d, J=100, L=1) {
  j <- 1:J
  #j <- seq(-L,L, len=J)
  ### CHECK THIS ###
  f <- function(k) 1 / (1 + k^2)
  #f <- function(k) 1 / (1 + k*1i)
  
  lam <- f(j * pi / (2 * L))
  
  psi <- function(d) {
    p <- 1i * j * pi * d / (2*L)
    exp_p <- exp(p)
    #exp_p / sum(exp_p) # do I normalize like this?
    exp_p
  }
  
  # Note: 1/J == psi(0)
  #out <- sum( lam * Conj(psi(0)) * psi(d) )
  #out <- sum( lam * psi(0) * Conj(psi(d)) )
  out <- sum(lam *  psi(d)) / sum(lam * psi(0)) # normalizing?
  
  #stopifnot(Im(out) == 0)
  #Re(out)
  Re(out)
}