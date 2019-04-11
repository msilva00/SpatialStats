# Covariance Functions
spherical = function(d, phi, sig2=1, nu=0) {
  ifelse (d > phi, 
          0,
          sig2 * (1 - 1.5 * d / phi + .5 * (d / phi) ^ 3))
}

pow_exp = function(d, phi, sig2, nu=1) {
  stopifnot(nu > 0 && nu <= 2)
  sig2 * exp(-abs(d / phi) ^ nu)
}

rational_quad = function(d, phi, sig2, nu=0) {
  sig2 * (1 - d^2 / (d^2 + phi^2) )
}

wave = function(d, phi, sig2, nu=0) {
  x <- d / phi
  sig2 * ifelse(d==0, 1, (sin(x) / x))
}

matern = function(d, phi, sig2, nu=0) {
  x <- sqrt(2 * nu) * d / phi
  sig2  / (2^(nu-1) * gamma(nu)) * x^nu * besselK(x, nu)
}

get_phi <- function(model, r=.05, d=1, nu=1, mn=1E-10, mx=10) {
  f <- function(phi) model(d, phi, sig2=1, nu) - r
  uniroot(f, c(mn,mx))
}



phi = rep(NA,5)


phi[1] = get_phi(spherical)$root
phi[2] = get_phi(rational_quad)$root
phi[3] = get_phi(wave)$root

phi[4] = get_phi(pow_exp)$root
phi[5] = get_phi(matern)$root

#### Covariogram Plots ####
# Spherical
xgrid = seq(0,2,length.out = 100)
plot(xgrid, spherical(xgrid,phi[1], sig2 = 1, nu = 1), type = 'l',
     ylab = "Covariogram",xlab = expression(paste("Distance ", (tau))),
     main = "Spherical", col = "blue")
legend("topright", legend=expression(paste(phi, " = 1.232")), lty=1, col="blue")

# Rational Quadratic
plot(xgrid, rational_quad(xgrid,phi[2], sig2 = 1, nu = 1), type = 'l',
     ylab = "Covariogram",xlab = expression(paste("Distance ", (tau))),
     main = "Rational Quadratic", col = "blue")
legend("topright", legend=expression(paste(phi, " = 0.229")), lty=1, col="blue")

# Wave
xgrid_wave = seq(0,5, length.out = 100)
plot(xgrid_wave, wave(xgrid_wave,phi[3], sig2 = 1, nu = 1), type = 'l',
     ylab = "Covariogram",xlab = expression(paste("Distance ", (tau))),
     main = "Wave", col = "blue", ylim = c(-0.5,1))
legend("topright", legend=expression(paste(phi, " = 0.334")), lty=1, col="blue")

# Powered Exponential
nus = seq(0.5,2, by = 0.5)
powered.exp.phi = rep(NA,4)
powered.exp.phi[1] = get_phi(pow_exp,nu = nus[1])$root
powered.exp.phi[2] = get_phi(pow_exp,nu = nus[2])$root
powered.exp.phi[3] = get_phi(pow_exp,nu = nus[3])$root
powered.exp.phi[4] = get_phi(pow_exp,nu = nus[4])$root
xgrid_powexp = seq(0,2.2, length.out = 100)
plot(xgrid_powexp, pow_exp(xgrid_powexp,powered.exp.phi[1], sig2 = 1, nu = nus[1]), type = 'l',
     ylab = "Covariogram",xlab = expression(paste("Distance ", (tau))),
     main = "Powered Exponential", col = "blue")
lines(xgrid_powexp, pow_exp(xgrid_powexp, powered.exp.phi[2], sig2 = 1, nu = nus[2]),
      type = 'l', col = "red")
lines(xgrid_powexp, pow_exp(xgrid_powexp, powered.exp.phi[3], sig2 = 1, nu = nus[3]),
      type = 'l', col = "green")
lines(xgrid_powexp, pow_exp(xgrid_powexp, powered.exp.phi[4], sig2 = 1, nu = nus[4]),
      type = 'l', col = "yellow")
legend("topright", legend=c(expression(paste(nu, " = 0.5")), expression(paste(nu, " = 1.0")), 
                            expression(paste(nu, " = 1.5")), expression(paste(nu, " = 2.0"))),
       col = c("blue", "red", "green", "yellow"), lty=1)

# Matern
matern.phi = rep(NA,4)
matern.phi[1] = get_phi(matern,nu = nus[1])$root
matern.phi[2] = get_phi(matern,nu = nus[2])$root
matern.phi[3] = get_phi(matern,nu = nus[3])$root
matern.phi[4] = get_phi(matern,nu = nus[4])$root
xgrid_matern = seq(0,3, length.out = 100)
plot(xgrid_matern, matern(xgrid_matern,matern.phi[1], sig2 = 1, nu = nus[1]), type = 'l',
     ylab = "Covariogram",xlab = expression(paste("Distance ", (tau))),
     main = "Matern", col = "blue", xlim = c(0,2))
lines(xgrid_matern, matern(xgrid_matern, matern.phi[2], sig2 = 1, nu = nus[2]),
      type = 'l', col = "red")
lines(xgrid_matern, matern(xgrid_matern, matern.phi[3], sig2 = 1, nu = nus[3]),
      type = 'l', col = "green")
lines(xgrid_matern, matern(xgrid_matern, matern.phi[4], sig2 = 1, nu = nus[4]),
      type = 'l', col = "yellow")
legend("topright", legend=c(expression(paste(nu, " = 0.5")), expression(paste(nu, " = 1.0")), 
                            expression(paste(nu, " = 1.5")), expression(paste(nu, " = 2.0"))),
       col = c("blue", "red", "green", "yellow"), lty=1)



#### Semi Variogram ####
semi_variogram = function(cov_fn, d, phi, sig2, nu, d_zero) {
  return(cov_fn(d_zero, phi, sig2, nu) - cov_fn(d, phi, sig2, nu))
}

# Spherical
plot(xgrid, semi_variogram(cov_fn = spherical,phi= phi[1], sig2 = 1, nu = 1, d_zero = 1E-10,d=xgrid),
     type = 'l',
     ylab = "Semi-Variogram",xlab = expression(paste("Distance ", (tau))),
     main = "Spherical", col = "red")
legend("bottomright", legend=expression(paste(phi, " = 1.232")), lty=1, col="red")

# Rational Quadratic
plot(xgrid, semi_variogram(cov_fn = rational_quad,phi= phi[2], sig2 = 1, nu = 1, d_zero = 1E-10,d=xgrid), 
     type = 'l',
     ylab = "Semi-variogram",xlab = expression(paste("Distance ", (tau))),
     main = "Rational Quadratic", col = "red")
legend("bottomright", legend=expression(paste(phi, " = 0.229")), lty=1, col="red")

# Wave
xgrid_wave = seq(0,5, length.out = 100)
plot(xgrid_wave, semi_variogram(cov_fn = wave,phi= phi[3], sig2 = 1, nu = 1, d_zero = 1E-10,d=xgrid_wave),
     type = 'l',
     ylab = "Semi-variogram",xlab = expression(paste("Distance ", (tau))),
     main = "Wave", col = "red", ylim = c(0,1.5))
legend("topright", legend=expression(paste(phi, " = 0.334")), lty=1, col="red")

# Powered Exponential
plot(xgrid_powexp, semi_variogram(cov_fn = pow_exp,phi= powered.exp.phi[1], sig2 = 1, nu = nus[1], d_zero = 1E-10,d=xgrid_powexp), 
     type = 'l',
     ylab = "Covariogram",xlab = expression(paste("Distance ", (tau))),
     main = "Powered Exponential", col = "blue")
lines(xgrid_powexp, semi_variogram(cov_fn = pow_exp,phi= powered.exp.phi[2], sig2 = 1, nu = nus[2], d_zero = 1E-10,d=xgrid_powexp),
      type = 'l', col = "red")
lines(xgrid_powexp, semi_variogram(cov_fn = pow_exp,phi= powered.exp.phi[3], sig2 = 1, nu = nus[3], d_zero = 1E-10,d=xgrid_powexp),
      type = 'l', col = "green")
lines(xgrid_powexp, semi_variogram(cov_fn = pow_exp,phi= powered.exp.phi[4], sig2 = 1, nu = nus[4], d_zero = 1E-10,d=xgrid_powexp),
      type = 'l', col = "yellow")
legend("bottomright", legend=c(expression(paste(nu, " = 0.5")), expression(paste(nu, " = 1.0")), 
                            expression(paste(nu, " = 1.5")), expression(paste(nu, " = 2.0"))),
       col = c("blue", "red", "green", "yellow"), lty=1)

# Matern
plot(xgrid_matern, semi_variogram(cov_fn = matern,phi= matern.phi[1], sig2 = 1, nu = nus[1], d_zero = 1E-10,d=xgrid_matern), 
     type = 'l',
     ylab = "Semi-variogram",xlab = expression(paste("Distance ", (tau))),
     main = "Matern", col = "blue")
lines(xgrid_matern, semi_variogram(cov_fn = matern,phi= matern.phi[2], sig2 = 1, nu = nus[2], d_zero = 1E-10,d=xgrid_matern),
      type = 'l', col = "red")
lines(xgrid_matern, semi_variogram(cov_fn = matern,phi= matern.phi[3], sig2 = 1, nu = nus[3], d_zero = 1E-10,d=xgrid_matern),
      type = 'l', col = "green")
lines(xgrid_matern, semi_variogram(cov_fn = matern,phi= matern.phi[4], sig2 = 1, nu = nus[4], d_zero = 1E-10,d=xgrid_matern),
      type = 'l', col = "yellow")
legend("bottomright", legend=c(expression(paste(nu, " = 0.5")), expression(paste(nu, " = 1.0")), 
                               expression(paste(nu, " = 1.5")), expression(paste(nu, " = 2.0"))),
       col = c("blue", "red", "green", "yellow"), lty=1)
