library(geoR)
set.seed(1)
phi = c(1.2324244, 0.5777683, 0.2294342, 0.3343117, 0.3338098, 0.2107771)
xgrid = seq(-3, 3, length = 100)



# Spherical Covariance Matrix
y = grf(grid = cbind(xgrid, 0), cov.model = "spherical", cov.pars = c(1, phi[1]), kappa = 0.5)
plot(xgrid, y$data, type='l', bty = 'n', main = "Spherical Covariance", main = "Spherical Covariance", xlab = "t", ylab = "X(t)")


# Powered Exponential Covariance Matrix
y = grf(grid = cbind(xgrid, 0), cov.model = "powered.exponential", cov.pars = c(1, phi[2]), kappa = 2)
plot(xgrid, y$data, type='l', bty = 'n', main = "Powered Exponential \n Covariance", col = "blue",
     ylim = c(-3,3), xlab = "t", ylab = "X(t)")


y = grf(grid = cbind(xgrid, 0), cov.model = "powered.exponential", cov.pars = c(1, phi[2]), kappa = 0.5)
lines(xgrid, y$data, type='l', bty = 'n', main = "Powered Exponential \n Covariance", col = "red")
legend("bottomleft", legend=c(expression(paste(nu, " = 2")), 
                               expression(paste(nu, " = 0.5"))), 
       col = c("red", "blue"), lty = 1)

# Rational Quadratic Covariance Matrix
y = grf(grid = cbind(xgrid, 0), cov.model = "cauchy", cov.pars = c(1, phi[3]), kappa = 1)
plot(xgrid, y$data, type='l', bty = 'n', main = "Rational Quadratic \n Covariance", col = "red",xlab = "t", ylab = "X(t)")

# Wave Covariance Matrix
y = grf(grid = cbind(xgrid, 0), cov.model = "wave", cov.pars = c(1, phi[4]), kappa = 0.5)
plot(xgrid, y$data, type='l', bty = 'n', main = "Wave Covariance", col = "red",xlab = "t", ylab = "X(t)")


# Matern Covariance Matrix
y = grf(grid = cbind(xgrid, 0), cov.model = "matern", cov.pars = c(1, phi[5]), kappa = 2)
plot(xgrid, y$data, type='l', bty = 'n', main = "Matern Covariance", col = "blue",
     ylim = c(-3,3),xlab = "t", ylab = "X(t)")
y = grf(grid = cbind(xgrid, 0), cov.model = "powered.exponential", cov.pars = c(1, phi[2]), kappa = 0.5)
lines(xgrid, y$data, type='l', bty = 'n', main = "Powered Exponential \n Covariance", col = "red")
legend("bottomleft", legend=c(expression(paste(nu, " = 2")), 
                              expression(paste(nu, " = 0.5"))), 
       col = c("red", "blue"), lty = 1)


