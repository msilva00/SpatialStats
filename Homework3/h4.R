### R
set.seed(1)

set.seed(2)

dat = read.table("alb_weightALL.csv", header = TRUE, sep = ",")
dat = dat[,-c(3,4,5)]
library(dplyr)
dat = sample_n(dat, 500)
head(dat)

y = box.cox(dat$alb_weight)

qqnorm(y)
qqline(y)

n = length(y)
lon = dat$longitude
lat = dat$latitude
loc = cbind(lon, lat)
s = cbind(lon,lat)

mod2 = lm(y~ 1 + lon + lat + lon*lat)
summary(mod2)

D = cbind(1, lon, lat, lon*lat)
k = NCOL(D)
beta = matrix(0, nburn + nmcmc, k)
dim(beta)
N = 30000
sig2 = rep(0, N)
sig2[1] = 1
tau2 = rep(0, N)
tau2[1] = 1
psi = rep(0, N)
psi[1] = 1
gamma2 = double(nburn + nmcmc)
gamma2[1] = 1

#### MCMC ####
dists = as.matrix(dist(loc[,1:2], upper = TRUE, diag= TRUE))
make.K = function(psi, gamma2, nu){
  1/gamma2*geoR::matern(dists, psi, nu) + diag(length(y))
}

#### Proposal Density ####
loglikesillrange = function(s, y, D, sig2, phi, nu, tau2){
  n = length(y)
  k = ncol(D)
  I_n = diag(n)
  K = make.K(psi, gamma2, nu)
  Kinv = solve(K)
  DVinvD = t(D) %*% Kinv %*% D
  beta.hat = solve(DKinvD) %*% (t(D) %*% Kinv %*% y)
  resid = y - D %*% beta.hat
  S2 = t(resid) %*% Kinv %*% resid
  ldS2_plus2b = determinant(S2, logarithm = T)$modulus + (2*b)
  
  ldK = determinant(K, logarithm = T)$modulus
  ldDKinvD = determinant(DKinvD, logarithm = T)$modulus
  
  # return((-n/2) * log(sig2) - (1/2) * ldV - S2/(2*sig2) - (k/2) * log(sig2) - (1/2)*ldDVinvD) + prior_psi(phi, gamma)
  p_prob = (-1/2)*ldK - (1/2)*ldDKinvD - (a + (m+k)/2)*ldS2_plus2b + prior_psi()
  return(list(K=K, Kinv=Kinv, beta.hat=beta.hat, DVinvD=DVinvD, S2=S2, p_prob = p_prob))
}





