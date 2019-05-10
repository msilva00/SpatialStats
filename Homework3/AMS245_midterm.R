# # FULL CODE
# setwd("~/homework3-spacialstats/")
# rm(list = ls())
# libs = c('geoR','sf','fields', 'ncdf4','tidyverse')
# sapply(libs, require, character.only = TRUE)
# library(ncdf4)
# 
# temp = list.files(pattern="*.nc")
# ncin = nc_open(temp[1])
# prob_values <- ncatt_get(ncin,0,"probability_values")
# 
# #### Filter the data ####
# # dissolve world_countries_boundaries to yield landmass
# land = read_sf(
#     dsn = "UIA_World_Countries_Boundaries",
#     layer = "UIA_World_Countries_Boundaries"
#     ) %>%
#   mutate(bleh = 'bleh') %>%
#   group_by(bleh) %>%
#   summarise()
# 
# proj_intersect = function(df){
#   ## projects points dataframe subset as a feature class
#   ## then calculates the intersection with land_in_box
#   ## returns set of points that intersected
#   points = df %>%
#     na.omit %>%
#     st_as_sf(
#       coords = c('lon','lat'),
#       crs = "+proj=longlat +datum=WGS84 +no_defs"
#     )
#   ol = st_intersection(points, land_in_box)[c('geometry','bhriso', 'prob_thres')]
#   return(ol)
# }
# ?ncar_get
# get_goes_data = function(path, minlat, maxlat, minlon, maxlon){
#   # subsets satellite data to geography in question;
#   #  looking specifically at observations on land.
# 
#   # subset land data to area in question
#   # build bounding box to subset the landmass data
#   box_raw = list(
#     rbind(c(minlon, minlat), c(minlon, maxlat),
#           c(maxlon, maxlat), c(maxlon, minlat),
#           c(minlon, minlat)
#       )
#     )
#   box_poly = st_polygon(box_raw, dim = 'XY')
#   box = st_sfc(
#     list(box_poly),
#     crs = "+proj=longlat +datum=WGS84 +no_defs"
#     ) %>% st_sf
#   # intersect landmass data with bounding box
#   land_in_box = st_intersection(land, box)
# 
#   # open satellite data, pull and filter data to yield points inside bounding box
#   goes = nc_open(path)
#   data = tibble(
#     lat = as.vector(ncvar_get(goes, 'latitude')),
#     lon = as.vector(ncvar_get(goes, 'longitude')),
#     bhriso = as.vector(ncvar_get(goes, 'BHRiso')),
#     prob_thresh = as.vector(ncvar_get(goes, 'ProbabilityThreshold'))
#   ) %>% filter(
#       !is.na(bhriso) & !is.na(prob_thresh) &
#       between(lon, minlon, maxlon) &
#       between(lat, minlat, maxlat)
#   )
# 
#   # create points from coordinates; intersect with land in box
#   points = data %>%
#     st_as_sf(
#       coords = c('lon','lat'),
#       crs = "+proj=longlat +datum=WGS84 +no_defs"
#     )
#   ol = st_intersection(points, land_in_box)[c('geometry','bhriso', 'prob_thresh')]
#   tab_of_ol = cbind(st_coordinates(ol), ol$bhriso, ol$prob_thresh)
#   colnames(tab_of_ol) = c('longitude','latitude','bhriso', 'prob_thresh')
#   return(tab_of_ol)
# }
# 
# path = 'GSA_AlbedoProd_GOES_075_VIS02_2000_181.nc'
# minlon = -90; maxlon = -60; minlat = -8; maxlat = 10.5
# 
# ex = get_goes_data(path, minlat, maxlat, minlon, maxlon)
# ex %>% as_tibble %>% write_csv('prob_thresh_included.csv')
# dat = as.data.frame(ex)
# 
# ## END OF FILTERING - UNCOMMENT ABOVE TO FILTER##
# 
set.seed(2)
# READ IN DATA ##
ex = read.csv("alb_weightALL.csv", sep = ",")

head(ex)
library(dplyr)
random500 = sample_n(ex, 500)
y = random500[,6]
require(graphics)
col.vals = (y - min(y)) / (diff(range(y)))
max(col.vals)



cols = rgb(1, col.vals, 0)
?rgb()
# cols = grey(col.vals)

layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2)

layout(mat = layout.matrix,
       # heights = c(1, 6), # Heights of the two rows
       widths = c(5, 2)) # Widths of the two columns

layout.show(2)
par(mar = c(5, 4, 4, 3))
library(scales)
maps::map("world", xlim = range(random500[,1])+c(-1-1,2+1), ylim = c(-8-1,10.5+1),
          lwd = 2)
axis(side = 2, at = seq(-9,11, length.out = 6))
axis(side = 1, at = seq(-85, -60, by = 5))

points(random500, col = alpha(cols,1), pch = 16, cex = 1, lwd = 3)
title(main = "Albedo Measurement",cex.main = 1.5, outer = T, line = -2)
par(mar = c(2, 0, 2, 5))
plot(x = rep(max(random500[,1]) + 0.67, 100),
     y = seq(min(random500[,2]), max(random500[,2]), length = 100),
     pch = 15, col = rgb(rep(1, length = 100), seq(0, 1, length = 100), 0),
     cex = 2.5, axes = false)
# lines(rep(max(random500[,1]) + 0.97, 2), range(random500[,2]))
for (i in 1:6){
  lines(max(random500[,1])+ 0.97 + c(0.97, 1.07), rep(min(random500[,2])+(i-1)*diff(range(random500[,2])/5), 2))
  #   text(max(random500[,1])+1.15, min(random500[,2])+(i-1)*diff(range(random500[,2])/5),
  #       round(quantile(y, (i-1)/5), 3), pos=4)
  text(max(random500[,1])+1.12, min(random500[,2])+(i-1)*diff(range(random500[,2])/5),
       round(seq(min(y), max(y), length = 6)[i], 2), pos=4)
}
par(mfrow = c(1,2), mar = c(3,3,5,3))
hist(log(random500[,3]), freq = F, main = "")
lines(density(log(random500[,3])))

qqnorm(log(random500[,3]), main = "")
qqline(log(random500[,3]))

title("Log Transformed Albedo Measurements", outer = T, line = -2)


par(mfrow = c(1,2), mar = c(3,3,5,3))
hist((random500[,3]), freq = F, main = "", ylim = c(0,8))
# lines(density((random500[,3])))

qqnorm((random500[,3]), main = "")
qqline((random500[,3]))

title("Albedo Measurements (w/out transformation)", outer = T, line = -2)

#### Residual Analysis/Covariate Analysis ####

##### LM #####
y = log(random500$alb_weight)
lon = random500$longitude
lat = random500$latitude
loc = cbind(lon, lat)
s = cbind(lon,lat)

mod = step(lm(y ~ . + .^2 + I(loc^2), data = data.frame(loc)), scope = list("lower" = lm(y ~ 1),
                                                                            "upper" = lm(y ~ . + .^2 + I(loc^2), data = data.frame(loc))),
           direction = "both")
summary(mod)



mod2 = lm(y~ 1 + lon + lat + lon*lat)
summary(mod2)

# Residuals
yhat = predict(mod2)
res = y - yhat


#### trend detection #####
plot(lon, y, type = 'n', bty = 'n', main = "Longitude", cex.main = 1.5,
     xlab = "Longitude", ylab = "log albedo")
segments(x0 = lon, y0 = y, y1 = yhat, col = "gray80")
points(lon, y, pch = 16, col = "gray70")
points(lon, yhat, pch = 16, col = "black")

plot(lat, y, type = 'n', bty = 'n', main = "Latitude", cex.main = 1.5,
     xlab = "Latitude", ylab = "log albedo")
segments(x0 = lat, y0 = y, y1 = yhat, col = "gray80")
points(lat, y, pch = 16, col = "gray70")
points(lat, yhat, pch = 16, col = "black")

#### Varios ####
library(geoR)
plot(variog(data=res, coords=s, op='cloud'))

# par(mfrow = c(1, 2))
tmp = variog(coords = cbind(lon,lat), data = y, trend = y ~ I(lon^2) + I(lat^2))
plot(tmp$u, tmp$v, pch = 16, ylim = c(0, max(tmp$v)), type = "b",
     main = "Binned Semivariogram \n for log Albedo", cex.main = 1.5,
     xlab = "Distance", ylab = "Semivariance")
tmp = variog4(coords = cbind(lon,lat), data = y, trend = y ~ I(lon^2) + I(lat^2))
tmp = variog4(coords = cbind(lon,lat), data = y, uvec = 13)
tmp.xlim = range(sapply(tmp[-5], function(x) range(x$u)))
tmp.ylim = c(0, max(sapply(tmp[-5], function(x) max(x$v, na.rm = TRUE))))
plot(0, type='n', bty='n', xlim = tmp.xlim, ylim = tmp.ylim,
     main = "Directional Semivariogram \n for log(albedo)", cex.main = 1.5,
     xlab = "Distance", ylab = "Semivariance")
for (i in 1:4)
  lines(tmp[[2]]$u, tmp[[i]]$v, col = i, lty = i)
legend("topleft", legend = expression(0^o, 45^o, 90^o, 135^o),
       lty = 1:4, col = 1:4)
par(mfrow = c(1, 1))
dev.off()
trend = y~lon + lat + I(lon^2)*I(lat^2)

#### LSE fit to Matern correlation ####
# C(u) = sig^2 * rho(u) + kappa * (u == 0)              covariance
# gamma(u) = sig^2 * (1 - rho(u)) + kappa * (u > 0)     semi-variogram
semiv = function(u, phi, nu, sig2, kappa){
  sig2 * (1-matern(u, phi, nu)) + kappa
}
make.f = function(nu){
  f.opt = function(p){
    phi = p[1]
    sig2 = p[2]
    kappa = p[3]
    sum( (semiv(obs$u, phi, nu, sig2, kappa) - obs$v)^2)
  }
  return (f.opt)
}
objects = c("empirical_1","empirical_2","empirical_3","empirical_4")
par(mfrow = c(1, 1))
# for (nu in c(0.5,1,1.5,2)){
#   pars = optim(rep(1.1, 3), make.f(nu), lower = c(0, 0, 0), method = "L-BFGS-B")$par
#   jpeg(paste(objects[nu], ".jpg", sep=""))
#   plot(obs$u, obs$v, ylim = c(0, max(obs$v)), type = "b")
#   lines(obs$u, semiv(obs$u, pars[1], nu, pars[2], pars[3]), col = 'red')
#   text(rep(1.25, 4), seq(1, 0.6, length = 4), c(round(pars, 3), nu),
#        pos = 4, cex = 1.5)
#   text(rep(0.5, 4), seq(1, 0.6, length = 4), expression(sigma^2, phi, kappa, nu),
#        pos = 4, cex = 1.5)
# }

summary(mod)

trend = y ~ 1 + lon + lat + lon * lat
obs = variog(coords = cbind(lon,lat), data = y, trend = trend)

nu = 0.5
(pars = optim(c(0.1, 5, 1), make.f(nu), lower = c(0, 0, 0), method = "L-BFGS-B")$par)
plot(obs$u, obs$v, ylim = c(0, max(obs$v)), type = "b", xlab = "", ylab = "", 
     main = paste("Empirical Variogram, nu = ", nu))
lines(obs$u, semiv(obs$u, pars[1], nu, pars[2], pars[3]), col = 'red')
text(rep(4, 4), seq(0.6,0.4, length = 4), c(round(pars, 3), nu),
     pos = 4, cex = 1.5)
text(rep(0.5, 4), seq(0.6, 0.4, length = 4), expression(sigma^2 ~ "= ", phi ~"= ", kappa~"= ", nu ~ "= "),
     pos = 4, cex = 1.5)
# (pars = nlm(p = rep(1, 3), f = make.f(nu))$estimate)
# lines(obs$u, semiv(obs$u, pars[1], nu, pars[2], pars[3]), col = 'blue', lty = 3, lwd = 3)

nu = 1
(pars = optim(rep(1.1, 3), make.f(nu), lower = c(0, 0, 0), method = "L-BFGS-B")$par)
plot(obs$u, obs$v, ylim = c(0, max(obs$v)), type = "b", xlab = "", ylab = "", 
     main = paste("Empirical Variogram, nu = ", nu))
lines(obs$u, semiv(obs$u, pars[1], nu, pars[2], pars[3]), col = 'red')
# (pars = nlm(p = rep(1, 3), f = make.f(nu))$estimate)
# lines(obs$u, semiv(obs$u, pars[1], nu, pars[2], pars[3]), col = 'blue', lty = 3, lwd = 3)
text(rep(4, 4), seq(0.6,0.4, length = 4), c(round(pars, 3), nu),
     pos = 4, cex = 1.5)
text(rep(0.5, 4), seq(0.6, 0.4, length = 4), expression(sigma^2 ~ "= ", phi ~"= ", kappa~"= ", nu ~ "= "),
     pos = 4, cex = 1.5)
dev.off()
nu = 1.5
(pars = optim(rep(1.1, 3), make.f(nu), lower = c(0, 0, 0), method = "L-BFGS-B")$par)
plot(obs$u, obs$v, ylim = c(0, max(obs$v)), type = "b", xlab = "", ylab = "", 
     main = paste("Empirical Variogram, nu = ", nu))
lines(obs$u, semiv(obs$u, pars[1], nu, pars[2], pars[3]), col = 'red')
# (pars = nlm(p = rep(1, 3), f = make.f(nu))$estimate)
# lines(obs$u, semiv(obs$u, pars[1], nu, pars[2], pars[3]), col = 'blue', lty = 3, lwd = 3)
text(rep(4, 4), seq(0.6,0.4, length = 4), c(round(pars, 3), nu),
     pos = 4, cex = 1.5)
text(rep(0.5, 4), seq(0.6, 0.4, length = 4), expression(sigma^2 ~ "= ", phi ~"= ", kappa~"= ", nu ~ "= "),
     pos = 4, cex = 1.5)

nu = 2.5
(pars = optim(rep(1.1, 3), make.f(nu), lower = c(0, 0, 0), method = "L-BFGS-B")$par)
plot(obs$u, obs$v, ylim = c(0, max(obs$v)), type = "b", xlab = "", ylab = "", 
     main = paste("Empirical Variogram, nu = ", nu))
lines(obs$u, semiv(obs$u, pars[1], nu, pars[2], pars[3]), col = 'red')
# (pars = nlm(p = rep(1, 3), f = make.f(nu))$estimate)
# lines(obs$u, semiv(obs$u, pars[1], nu, pars[2], pars[3]), col = 'blue', lty = 3, lwd = 3)
text(rep(4, 4), seq(0.6,0.4, length = 4), c(round(pars, 3), nu),
     pos = 4, cex = 1.5)
text(rep(0.5, 4), seq(0.6, 0.4, length = 4), expression(sigma^2 ~ "= ", phi ~"= ", kappa~"= ", nu ~ "= "),
     pos = 4, cex = 1.5)

#### best params ####
library(geoR)
init <- expand.grid(seq(0,100, len=10), seq(0,1,len=10))

vario = variog(data = res, coords = loc, messages = F)
(vf1 <- variofit(vario, ini.cov.pars=init, kappa=0.5, fix.nug=F, messages = F))
(vf2 <- variofit(vario, ini.cov.pars=init, kappa=1.0, fix.nug=F,messages = F))
(vf3 <- variofit(vario, ini.cov.pars=init, kappa=1.5,  fix.nug=F,messages = F))
(vf4 <- variofit(vario, ini.cov.pars=init, kappa=2.5, fix.nug=F,messages = F))
vf <- list(vf1, vf2, vf3, vf4)



#### log like functions ####
loglikeSillRange <- function(sig2,tau2,phi,kappa,y,X,s,m) {
  n <- length(y)
  I_n <- diag(n)
  D <- as.matrix(dist(s))
  V <- sig2 * geoR::matern(D, phi=phi, kappa=kappa) + tau2 * I_n
  Vi <- solve(V)
  XViX <- t(X) %*% Vi %*% X
  b.hat <- solve(XViX) %*% (t(X) %*% Vi %*% y) 
  res <- y - X %*% b.hat
  S <- t(res) %*% Vi %*% res
  
  # ld_V <- unlist(determinant(V, log=TRUE))[1]
  # ld_XViX <- unlist(determinant(XViX, log=TRUE))[1]
  ld_V <- determinant(V, logarithm = T)$modulus
  ld_XViX <- determinant(XViX, logarithm = T)$modulus
  # (-1/2) * (ld_V + ld_XViX) - S/2
  (-1/2) * ((n+m)* log(sig2) + ld_V + ld_XViX + S^2/sig2)
}

loglikeRange <- function(phi, kappa, tau2OverSig2, y, X, s) {
  n <- length(y)
  k <- ncol(X)
  I_n <- diag(n)
  D <- as.matrix(dist(s))
  V <- geoR::matern(D, phi=phi, kappa=kappa) + tau2OverSig2 * I_n
  Vi <- solve(V)
  XViX <- t(X) %*% Vi %*% X
  b.hat <- solve(XViX) %*% (t(X) %*% Vi %*% y) 
  res <- y - X %*% b.hat
  S <- t(res) %*% Vi %*% res
  
  #ld_V <- unlist(determinant(V, log=TRUE))[1]
  #ld_XViX <- unlist(determinant(XViX, log=TRUE))[1]
  ld_V <- log(det(V))
  ld_XViX <- log(det(XViX))
  (-1/2) * (ld_V + ld_XViX) - (n-k)/2*log(S)
}

# #### Evaluate at log like sill and range ####
# # takes a long time to run
# 
# X <- as.matrix(cbind(1,mod2$model[,-1]))
# NCOL(X)
# 
# kappa_list <- as.list(c(0.5, 1, 1.5, 2.5))
# J <- 20
# 
# sig2.grid <- seq(0, 1, len=J)
# phi.grid <- seq(0, 20, len=J)
# #sig2.grid <- seq(0, 100, len=J)
# #phi.grid <- seq(0, 100, len=J)
# out <- lapply(kappa_list, function(k) {
#   out <- matrix(NA, J, J)
#   for (i in 1:J) {
#     print(paste("i:", i))
#     for (j in 1:J) {
#       out[i,j] <- loglikeSillRange(sig2=sig2.grid[j], phi=phi.grid[i],
#                                    tau2=vf.bestfit$nugget, kappa=k, y, X, s=s, m =NCOL(X))
#       # print(paste("j:", j))
#     }
#   }
#   out
# })
# # save the file because its a bitch to re-run
# write.csv(out, file = "likelihood_out.txt")
# 
# #### Kappa 0.5 ####
# (mxout = max(out[[1]], na.rm = T))
# clevels = seq(mxout-700, mxout, length.out = 10)
# 
# par(mfrow = c(1,2))
# par(mar = c(5, 4, 4, 2))
# contour(x = phi.grid[1:10], y = sig2.grid[1:20], out[[1]][1:10,1:20], 
#         xlab='phi', ylab='sig2', levels = clevels)
# 
# 
# par(mar = c(0, 1, 0, 0))
# persp( x = phi.grid, y = sig2.grid, out[[1]],
#        ylab = " \n\n\n sigma^2\n", xlab = "\n \n phi \n\n\n\n", theta = 350, phi = 20,
#        zlab = 'log likelihood', col = 'white' )
# title("nu = 0.5", outer = T,line = -2)
# 
# #### nu 1 ####
# (mxout = max(out[[2]], na.rm = T))
# clevels = seq(mxout-700, mxout, length.out = 10)
# 
# par(mfrow = c(1,2))
# par(mar = c(5, 4, 4, 2))
# contour(phi.grid[1:5], sig2.grid, out[[2]][1:5,], ylab='sig2', xlab='phi', levels = clevels)
# 
# 
# par(mar = c(0, 1, 0, 0))
# persp( phi.grid, sig2.grid, out[[2]],
#        ylab = "\n\n    sigma^2", xlab = "phi", theta = 320, phi = 5,
#        zlab = 'log likelihood', col = 'white' )
# title("nu = 1", outer = T,line = -2)
# 
# #### nu 1.5 ####
# (mxout = max(out[[3]], na.rm = T))
# clevels = seq(mxout-1000, mxout, length.out = 15)
# 
# par(mfrow = c(1,2))
# par(mar = c(5, 4, 4, 2))
# contour(phi.grid[1:3], sig2.grid, out[[3]][1:3,], ylab='sig2', xlab='phi', levels = clevels)
# 
# par(mar = c(0, 2, 0, 0))
# persp( phi.grid, sig2.grid, out[[3]],
#        ylab = "\n\n     sigma^2", xlab = "phi", theta = 350, phi = 5,
#        zlab = 'log likelihood', col = 'white' )
# title("nu = 1.5", outer = T,line = -2)
# 
# #### nu 2 ####
# (mxout = max(out[[4]], na.rm = T))
# clevels = seq(mxout-700, mxout, length.out = 10)
# 
# par(mfrow = c(1,2))
# par(mar = c(5, 4, 4, 2))
# contour(phi.grid, sig2.grid, out[[4]], ylab='sig2', xlab='phi', levels = clevels)
# 
# par(mar = c(0, 2, 0, 0))
# persp( phi.grid, sig2.grid, out[[4]],
#        ylab = "\n\n sigma^2",xlab = "phi", theta = 320, phi = 5,
#        zlab = 'log likelihood', col = 'white' )
# title("nu = 2.5", outer = T,line = -2)
# 
# #### Evaluate at marginal for optimum values of nu####
# # in this case kappa = nu = shape parameter (bad choice of variable names)
# tn = c(.212,.25,.254,.245)
# out_phi <- as.list(1:4)
# for (i in 1:4){
#   out_phi[[i]] <- sapply(phi.grid, function(phi) {
#     loglikeRange(phi=phi, 
#                  tau2OverSig2= tn[i],
#                  kappa=kappa_list[[i]], y=y, X=X, s=s)
#   })
# }
# 
# for (i in 1:4) {
#   plot(phi.grid, out_phi[[i]],type = 'b', xlab='phi', ylab='loglike', 
#        main=paste0('nu = ', kappa_list[[i]]))
#   abline(v=phi.grid[which.max(out_phi[[i]])], lty=2)
# }
# 

### MCMC
set.seed(2)
test = sort(sample(nrow(random500), 140))
train = (1:nrow(random500))[-test]

# D = cbind(1, loc[,2], loc[,3], loc[,2]^2, loc[,2] * loc[,3])
#D = cbind(1, lon, lat, lon*lat, lon^2, lat^2, ele)
D = cbind(1, lon, lat, lon*lat)
#D = matrix(1, nrow = n)
X = y
k = NCOL(D)

D.test = D[test,]
X.test = X[test]
D = D[train,]
X = X[train]

n = length(X)

dists = as.matrix(dist(loc[,1:2], upper = TRUE, diag= TRUE))
make.K = function(psi, gamma2, nu){
  1/gamma2*geoR::matern(dists, psi, nu) + diag(length(y))
}

# For tau^2 (observational error, invGamma)
prior.a = 3
prior.b = 2

# For psi (range, Gamma)
prior.c = 1
prior.d = 1

# For gamma^2 (= tau^2 / sigma^2, obs variance / spatial variance, Gamma)
prior.e = 1
prior.f = 1

nburn = 5000
nmcmc = 10000
window = 100
param.beta = matrix(1, nburn + nmcmc, k)
param.tau2 = double(nburn + nmcmc)
param.psi = double(nburn + nmcmc)
param.gamma2 = double(nburn + nmcmc)
param.nu = double(nburn + nmcmc)

nu.states = c(0.5, 1, seq(1.5, 3.5, by = 1))
param.tau2[1] = 1
param.psi[1] = 1
param.gamma2[1] = 1
param.nu[1] = nu.states[1]
cand.sig = 1e-3 * diag(2)
accept = double(nburn + nmcmc)
vec.like = double(nburn + nmcmc)

list.K = rep(list(matrix(1, length(y), length(y))), nmcmc)

# system.time({
# for (j in 1:100){
#     K = make.K(param.psi[1], param.gamma2[1], nu)
#     U = chol(K)
#     U.inv = backsolve(U, diag(n))
#     DU = t(D) %*% U.inv
#     CDU = tcrossprod(DU)
#     beta.hat = solve(CDU) %*% (DU %*% (t(U.inv) %*% X))
#     S2 = sum((t(X - D %*% beta.hat) %*% U.inv)^2)
#     Det1 = -sum(log(diag(U)))
#     Det2 = -0.5*determinant(CDU)$modulus[1]
#     post.curr = Det1 + Det2 + (-(n-k)/2-prior.a) + log(S2 + 2*prior.b) +
#         dgamma(param.psi[1], prior.c, prior.d, log = TRUE) +
#         dgamma(param.gamma2[1], prior.e, prior.f, log = TRUE)
#     }
#     })

dim(K)
pred.x = matrix(1, nmcmc, length(X.test))
K = make.K(param.psi[1], param.gamma2[1], param.nu[1])
K.inv = solve(K[train, train])
beta.hat = solve(t(D) %*% K.inv %*% D) %*% t(D) %*% K.inv %*% X
S2 = t(X - D %*% beta.hat) %*% K.inv %*% (X - D %*% beta.hat)
Det1 = -0.5*determinant(K)$modulus[1] 
Det2 = -0.5*determinant(t(D) %*% K.inv %*% D)$modulus[1]
post.curr = Det1 + Det2 + (-(n-k)/2-prior.a) * log(S2 + 2*prior.b) +
  dgamma(param.psi[1], prior.c, prior.d, log = TRUE) +
  dgamma(param.gamma2[1], prior.e, prior.f, log = TRUE)

vec.like[1] = Det1 - n/2*log(param.tau2[1]) -
  1/(2*param.tau2[1])*t(X - D %*% param.beta[1,]) %*% K.inv %*% (X - D %*% param.beta[1,])

library(MASS)

autotune = function(accept, target = 0.25, k = 2.5){
  (1+(cosh(accept-target)-1)*(k-1)/(cosh(target-
                                           ceiling(accept-target))-1))^sign(accept-target)
}

for (i in 2:(nburn + nmcmc)){
  cat(i, "/", nburn + nmcmc, "\r")
  param.beta[i,] = param.beta[i-1,]
  param.tau2[i] = param.tau2[i-1]
  param.psi[i] = param.psi[i-1]
  param.gamma2[i] = param.gamma2[i-1]
  param.nu[i] = param.nu[i-1]
  vec.like[i] = vec.like[i-1]
  
  cand = c(mvrnorm(1, c(param.psi[i-1], param.gamma2[i-1]), cand.sig),
           sample(nu.states, 1))
  if (all(cand > 0)){
    cand.K = make.K(cand[1], cand[2], cand[3])
    cand.K.inv = chol2inv(chol(cand.K[train, train]))
    cand.DKD.inv = solve(t(D) %*% cand.K.inv %*% D)
    
    cand.beta.hat = cand.DKD.inv %*% (t(D) %*% (cand.K.inv %*% X))
    cand.S2 = t(X - D %*% cand.beta.hat) %*% cand.K.inv %*% (X - D %*% cand.beta.hat)
    cand.Det1 = -0.5*determinant(cand.K)$modulus[1]
    cand.Det2 = 0.5*determinant(cand.DKD.inv)$modulus[1]
    
    post.cand = cand.Det1 + cand.Det2 + (-(n-k)/2-prior.a) * log(cand.S2 + 2*prior.b) + 
      dgamma(cand[1], prior.c, prior.d, log = TRUE) +
      dgamma(cand[2], prior.e, prior.f, log = TRUE)
    
    if (log(runif(1)) <= post.cand - post.curr){
      param.psi[i] = cand[1]
      param.gamma2[i] = cand[2]
      param.nu[i] = cand[3]
      accept[i] = 1
      
      K = cand.K
      K.inv = cand.K.inv
      DKD.inv = cand.DKD.inv
      
      beta.hat = cand.beta.hat
      S2 = cand.S2
      Det1 = cand.Det1
      Det2 = cand.Det2
      
      post.curr = post.cand
      
      param.tau2[i] = 1/rgamma(1, prior.a + (n-k)/2, prior.b + S2/2)
      param.beta[i,] = mvrnorm(1, beta.hat, param.tau2[i] * DKD.inv)
      
      vec.like[i] = Det1 - n/2*log(param.tau2[i]) -
        1/(2*param.tau2[i])*t(X - D %*% param.beta[i,]) %*% K.inv %*% (X - D %*% param.beta[i,])
    }
  }
  
  if (i > nburn){
    j = i - nburn
    list.K[[j]] = K
    ### Predictions
    mu = D.test %*% param.beta[i,] + list.K[[j]][test, train] %*% 
      (K.inv %*% (X - D %*% param.beta[i,]))
    V = list.K[[j]][test, test] - list.K[[j]][test, train] %*% 
      K.inv %*% list.K[[j]][train, test]
    pred.x[j,] = mvrnorm(1, mu, V)
  }
  
  if ((floor(i/window) == i/window) && (i <= nburn))
    cand.sig = autotune(mean(accept[(i-window+1):i]), target = 0.234, k = window / 50) *
    (cand.sig + window * var(cbind(param.psi, param.gamma2)[(i-window+1):i,]) / i)
  
}

param.beta = tail(param.beta, nmcmc)
param.tau2 = tail(param.tau2, nmcmc)
param.psi = tail(param.psi, nmcmc)
param.gamma2 = tail(param.gamma2, nmcmc)
param.nu = tail(param.nu, nmcmc)
vec.like = tail(vec.like, nmcmc)
accept = tail(accept, nmcmc)
sig2 = param.tau2 / param.gamma2


mean(accept)
par(mfrow = c(1,1))

plot(vec.like, type = 'l')
plot(param.psi, type = 'l')
plot(param.tau2, type = 'l')
plot(param.gamma2, type = 'l')
plot(sig2, type = 'l')

plot(param.beta[,c(1,2)], pch = 16)

colMeans(param.beta)
coef(mod)

plot(density(param.psi))
curve(dgamma(x, prior.c, prior.d), col = 'blue', add = TRUE)

plot(density(param.gamma2))
curve(dgamma(x, prior.e, prior.f), add = TRUE, col = 'blue')

plot(table(param.nu) / nmcmc)

(table(param.nu) / nmcmc)

plot(density(sig2))

mean(sig2)
mean(param.psi)
mean(param.tau2)

t(rbind(colMeans(param.beta),
        apply(param.beta, 2, quantile, c(0.025, 0.975))))

tmp = cbind(param.tau2, param.psi, param.gamma2)
t(rbind(apply(tmp, 2, mean),
        apply(tmp, 2, quantile, c(0.025, 0.975))))

### DIC
D.theta = -2*vec.like
var(D.theta)/2 + mean(D.theta)

### Predictions
pdf("./figs/prediction.pdf", width = 12, height = 8)
par(mfrow = c(3, 4))
for (i in 1:length(X.test)){
  plot(density(pred.x[,i]), main = i, xlab = "predicted values")
  abline(v = X.test[i], col = 'green', lwd = 2)
}
par(mfrow = c(1, 1))
dev.off()


pred.x = D %*% colMeans(param.beta)
post.K = make.K(mean(param.psi), mean(param.gamma2), nu)
mse = t(X - pred.x) %*% post.K %*% (X - pred.x)

plot(X, pred.x, pch = 16, bty = 'n')
abline(0, 1)

