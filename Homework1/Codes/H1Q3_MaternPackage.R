library(geoR)

set.seed(9)
type = c("spherical", "powered.exponential", "cauchy", "wave", "matern", "matern", "matern","matern")
kappa = c(0.5, 2, 1, 0.5, 1, 0.5, 2, 5)
ranges = double(length(type))

# Solve for the range yielding 0.05 correlation at unit distance
for (i in 1:length(type)){
  ranges[i] = uniroot(f = function(x){
    cov.spatial(1, cov.model = type[i], cov.pars = c(1, x), kappa = kappa[i]) - 0.05},
    c(0.001, 10))$root
}

h = seq(0, 5, length = 1000)
y = cov.spatial(h, cov.model = type[5], cov.pars = c(1, ranges[5]), kappa = kappa[5])
plot(h, y, xlim = c(0,2),xlab = expression(tau),ylab = "Covariogram", type = 'l', lwd = 1, bty = 'n', col = "red")
title(main = "Mat\u{E8}rn") 

h = seq(0, 5, length = 1000)
y2 = cov.spatial(h, cov.model = type[6], cov.pars = c(1, ranges[6]), kappa = kappa[6])
lines(h, y2, type = 'l', lwd = 1, bty = 'n', col = "blue")

h = seq(0, 5, length = 1000)
y3 = cov.spatial(h, cov.model = type[7], cov.pars = c(1, ranges[7]), kappa = kappa[7])
lines(h, y3, type = 'l', lwd = 1, bty = 'n', col = "green")

h = seq(0, 5, length = 1000)
y4 = cov.spatial(h, cov.model = type[8], cov.pars = c(1, ranges[8]), kappa = kappa[8])
lines(h, y4, type = 'l', lwd = 1, bty = 'n', col = "purple")

legend("topright", legend=c(expression(paste(nu, " = 1")), expression(paste(nu, " = 1/2")), 
                            expression(paste(nu, " = 2")), expression(paste(nu, " = 5"))),
       col = c("red", "blue", "green", "purple"), lty=1)

