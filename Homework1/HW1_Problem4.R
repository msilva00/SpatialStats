library(geoR)
set.seed(1)
phi = c(1.2324244, 0.5777683, 0.2294342, 0.3343117, 0.3338098, 0.2107771)
xgrid = seq(-3, 3, length = 100)



# Spherical
y = grf(grid = cbind(xgrid, 0), cov.model = "spherical", cov.pars = c(1, phi[1]), kappa = 0.5)
plot(xgrid, y$data, type='l', bty = 'n', main = "Spherical Covariance")


# Powered Exponential
y = grf(grid = cbind(xgrid, 0), cov.model = "spherical", cov.pars = c(1, phi[2]), kappa = 2)
plot(xgrid, y$data, type='l', bty = 'n', main = "Powered Exponential \n Covariance", col = "red")
y = grf(grid = cbind(xgrid, 0), cov.model = "spherical", cov.pars = c(1, phi[2]), kappa = 0.5)
lines(xgrid, y$data, type='l', bty = 'n', main = "Powered Exponential \n Covariance", col = "blue")
legend("bottomleft", legend=c(expression(paste(nu, " = 2")), 
                               expression(paste(nu, " = 0.5"))), 
       col = c("red", "blue"), lty = 1)

