library(geoR)
# Variogram models with the same "practical" range:
#
v.f <- function(x, ...){cov.spatial(x, ...)-0.05}

sv.f = function(x,...){1*1-(cov.spatial(x,...)-0.05)}
#### Covariograms #####
# Spherical
curve(v.f(x, cov.pars = c(1, 1.232), kappa = 0.5, cov.model = "sph"), 0, 2,
     lty = 1, ylab = "Covariogram",xlab = expression(paste("Distance ", (tau))),
     main = "Spherical", col = "blue")
legend("topright", legend=expression(paste(phi, " = 1.232")), lty=1, col="blue")

# Powered Exp
curve(v.f(x, cov.pars = c(1, 0.578), kappa = 2, cov.model = "powered.exponential"),
      0, 2, lwd = 1, col = "black", xlab = expression(paste("Distance", (tau))),
      ylab = "Covariogram", main = "Powered Exponential")
curve(v.f(x, cov.pars = c(1, 0.578), kappa = 0.5, cov.model = "powered.exponential"),
      0, 2,  add = TRUE,lwd = 1, col = "green")
curve(v.f(x, cov.pars = c(1, 0.578), kappa = 1.5, cov.model = "powered.exponential"),
      0, 2, add = TRUE, lwd = 1, col = "blue")
curve(v.f(x, cov.pars = c(1, 0.578), kappa = 1, cov.model = "powered.exponential"),
      0, 2, add = TRUE, lwd = 1, col = "red")
legend("topright", legend=c(expression(paste(nu, " = 2")), expression(paste(nu, " = 1")), 
                            expression(paste(nu, " = 1.5")), expression(paste(nu, " = 0.5"))),
       col = c("black", "red", "blue", "green"), lty=1)

# Rational Quadratic / Cauchy
curve(v.f(x, cov.pars = c(1, 0.2294342), kappa = 1.0, cov.model = "cauchy"), 0, 2,
      lty = 1, ylim = c(0,1), ylab = "Covariogram",xlab = expression(paste("Distance ", (tau))),
      main = "Rational Quadratic", col = "blue")
legend("topright", legend = expression(paste(phi, " = 0.229")), lty = 1, col = "blue")

# Wave
curve(v.f(x, cov.pars = c(1, 0.3343117), kappa = 0.5, cov.model = "wave"), 0, 5,
      lty = 1, ylab = "Covariogram",xlab = expression(paste("Distance ", (tau))),
      main = "Wave", col = "blue")
legend("topright", legend = expression(paste(phi, " = 0.334")), lty = 1, col = "blue")


# Matern
curve(v.f(x, cov.pars = c(1, 0.3), kappa = 2, cov.model = "matern"),
      0, 2, lwd = 1, col = "black", xlab = expression(paste("Distance", (tau))),
      ylab = "Covariogram", main = "Matèrn")
curve(v.f(x, cov.pars = c(1, 0.3), kappa = 0.5, cov.model = "matern"),
      0, 2,  add = TRUE,lwd = 1, col = "green")
curve(v.f(x, cov.pars = c(1, 0.3), kappa = 1.5, cov.model = "matern"),
      0, 2, add = TRUE, lwd = 1, col = "blue")
curve(v.f(x, cov.pars = c(1, 0.3), kappa = 1, cov.model = "matern"),
      0, 2, add = TRUE, lwd = 1, col = "red")
legend("topright", legend=c(expression(paste(nu, " = 2")), expression(paste(nu, " = 1")), 
                            expression(paste(nu, " = 1.5")), expression(paste(nu, " = 0.5"))),
       col = c("black", "red", "blue", "green"), lty=1)


#### Semi - variograms #####
# Spherical
curve(sv.f(x, cov.pars = c(1, 1.232), kappa = 0.5, cov.model = "sph"), 0, 2,
      lty = 1, col = "red", ylab = "Semi-variogram",xlab = expression(paste("Distance ", (tau))),
      main = "Spherical")
legend("bottomright", legend=expression(paste(phi, " = 1.232")), lty=1, col = "red")

# Powered Exp
curve(sv.f(x, cov.pars = c(1, 0.578), kappa = 2, cov.model = "powered.exponential"),
      0, 2, lwd = 1, col = "black", xlab = expression(paste("Distance", (tau))),
      ylab = "Covariogram", main = "Powered Exponential")
curve(sv.f(x, cov.pars = c(1, 0.578), kappa = 0.5, cov.model = "powered.exponential"),
      0, 2,  add = TRUE,lwd = 1, col = "green")
curve(sv.f(x, cov.pars = c(1, 0.578), kappa = 1.5, cov.model = "powered.exponential"),
      0, 2, add = TRUE, lwd = 1, col = "blue")
curve(sv.f(x, cov.pars = c(1, 0.578), kappa = 1, cov.model = "powered.exponential"),
      0, 2, add = TRUE, lwd = 1, col = "red")
legend("bottomright", legend=c(expression(paste(nu, " = 2")), expression(paste(nu, " = 1")), 
                            expression(paste(nu, " = 1.5")), expression(paste(nu, " = 0.5"))),
       col = c("black", "red", "blue", "green"), lty=1)

# Rational Quadratic / Cauchy
curve(sv.f(x, cov.pars = c(1, 0.2294342), kappa = 1.0, cov.model = "cauchy"), 0, 2,
      lty = 1, ylim = c(0,1), ylab = "Covariogram",xlab = expression(paste("Distance ", (tau))),
      main = "Rational Quadratic", col = "red")
legend("bottomright", legend = expression(paste(phi, " = 0.229")), lty = 1, col = "red")

# Wave
curve(sv.f(x, cov.pars = c(1, 0.3343117), kappa = 0.5, cov.model = "wave"), 0, 5,
      lty = 1, ylab = "Covariogram",xlab = expression(paste("Distance ", (tau))),
      main = "Wave", col = "red")
legend("bottomright", legend = expression(paste(phi, " = 0.334")), lty = 1, col = "red")


# Matern
curve(sv.f(x, cov.pars = c(1, 0.3), kappa = 2, cov.model = "matern"),
      0, 2, lwd = 1, col = "black", xlab = expression(paste("Distance", (tau))),
      ylab = "Covariogram", main = "Matèrn")
curve(sv.f(x, cov.pars = c(1, 0.3), kappa = 0.5, cov.model = "matern"),
      0, 2,  add = TRUE,lwd = 1, col = "green")
curve(sv.f(x, cov.pars = c(1, 0.3), kappa = 1.5, cov.model = "matern"),
      0, 2, add = TRUE, lwd = 1, col = "blue")
curve(sv.f(x, cov.pars = c(1, 0.3), kappa = 1, cov.model = "matern"),
      0, 2, add = TRUE, lwd = 1, col = "red")
legend("bottomright", legend=c(expression(paste(nu, " = 2")), expression(paste(nu, " = 1")), 
                            expression(paste(nu, " = 1.5")), expression(paste(nu, " = 0.5"))),
       col = c("black", "red", "blue", "green"), lty=1)


