# FULL CODE
setwd("~/SpatialStats/Homework3/")
rm(list = ls())
libs = c('geoR','sf','fields', 'ncdf4','tidyverse')
sapply(libs, require, character.only = TRUE)


#### Filter the data ####
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
#   ol = st_intersection(points, land_in_box)[c('geometry','bhriso')]
#   return(ol)
# }
# 
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
#     bhriso = as.vector(ncvar_get(goes, 'BHRiso'))
#   ) %>% filter(
#       !is.na(bhriso) &
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
#   ol = st_intersection(points, land_in_box)[c('geometry','bhriso')]
#   tab_of_ol = cbind(st_coordinates(ol), ol$bhriso)
#   colnames(tab_of_ol) = c('longitude','latitude','bhriso')
#   return(tab_of_ol)
# }
# 
# path = 'GSA_AlbedoProd_GOES_075_VIS02_2000_181.nc'
# minlon = -90; maxlon = -60; minlat = -8; maxlat = 10.5
# 
# ex = get_goes_data(path, minlat, maxlat, minlon, maxlon)
# ex %>% as_tibble %>% write_csv('subsetted_data.csv')
# dat = as.data.frame(ex)

## END OF FILTERING - UNVCOMMENT ABOVE TO FILTER##

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
library(VGAM)
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



mod = lm(y~ lon + lat + lon*lat)
summary(mod)

# Residuals
yhat = predict(mod)
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
for (i in 1:4){
  lines(tmp[[2]]$u, tmp[[i]]$v, col = i, lty = i)
}
legend("topleft", legend = expression(0^o, 45^o, 90^o, 135^o),
       lty = 1:4, col = 1:4)
par(mfrow = c(1, 1))
dev.off()

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
init <- expand.grid(seq(0,100, len=10), seq(0,1,len=10))

vario = variog(data = resid, coords = loc, messages = F)
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

#### Evaluate at log like sill and range ####
# takes a long time to run

X <- as.matrix(cbind(1,mod2$model[,-1]))
NCOL(X)

kappa_list <- as.list(c(0.5, 1, 1.5, 2.5))
J <- 20

sig2.grid <- seq(0, 1, len=J)
phi.grid <- seq(0, 20, len=J)
#sig2.grid <- seq(0, 100, len=J)
#phi.grid <- seq(0, 100, len=J)
out <- lapply(kappa_list, function(k) {
  out <- matrix(NA, J, J)
  for (i in 1:J) {
    print(paste("i:", i))
    for (j in 1:J) {
      out[i,j] <- loglikeSillRange(sig2=sig2.grid[j], phi=phi.grid[i],
                                   tau2=vf.bestfit$nugget, kappa=k, y, X, s=s, m =NCOL(X))
      # print(paste("j:", j))
    }
  }
  out
})
# save the file because its a bitch to re-run
write.csv(out, file = "likelihood_out.txt")

#### Kappa 0.5 ####
(mxout = max(out[[1]], na.rm = T))
clevels = seq(mxout-700, mxout, length.out = 10)

par(mfrow = c(1,2))
par(mar = c(5, 4, 4, 2))
contour(x = phi.grid[1:10], y = sig2.grid[1:20], out[[1]][1:10,1:20], 
        xlab='phi', ylab='sig2', levels = clevels)


par(mar = c(0, 1, 0, 0))
persp( x = phi.grid, y = sig2.grid, out[[1]],
       ylab = " \n\n\n sigma^2\n", xlab = "\n \n phi \n\n\n\n", theta = 350, phi = 20,
       zlab = 'log likelihood', col = 'white' )
title("nu = 0.5", outer = T,line = -2)

#### nu 1 ####
(mxout = max(out[[2]], na.rm = T))
clevels = seq(mxout-700, mxout, length.out = 10)

par(mfrow = c(1,2))
par(mar = c(5, 4, 4, 2))
contour(phi.grid[1:5], sig2.grid, out[[2]][1:5,], ylab='sig2', xlab='phi', levels = clevels)


par(mar = c(0, 1, 0, 0))
persp( phi.grid, sig2.grid, out[[2]],
       ylab = "\n\n    sigma^2", xlab = "phi", theta = 320, phi = 5,
       zlab = 'log likelihood', col = 'white' )
title("nu = 1", outer = T,line = -2)

#### nu 1.5 ####
(mxout = max(out[[3]], na.rm = T))
clevels = seq(mxout-1000, mxout, length.out = 15)

par(mfrow = c(1,2))
par(mar = c(5, 4, 4, 2))
contour(phi.grid[1:3], sig2.grid, out[[3]][1:3,], ylab='sig2', xlab='phi', levels = clevels)

par(mar = c(0, 2, 0, 0))
persp( phi.grid, sig2.grid, out[[3]],
       ylab = "\n\n     sigma^2", xlab = "phi", theta = 350, phi = 5,
       zlab = 'log likelihood', col = 'white' )
title("nu = 1.5", outer = T,line = -2)

#### nu 2 ####
(mxout = max(out[[4]], na.rm = T))
clevels = seq(mxout-700, mxout, length.out = 10)

par(mfrow = c(1,2))
par(mar = c(5, 4, 4, 2))
contour(phi.grid, sig2.grid, out[[4]], ylab='sig2', xlab='phi', levels = clevels)

par(mar = c(0, 2, 0, 0))
persp( phi.grid, sig2.grid, out[[4]],
       ylab = "\n\n sigma^2",xlab = "phi", theta = 320, phi = 5,
       zlab = 'log likelihood', col = 'white' )
title("nu = 2.5", outer = T,line = -2)

#### Evaluate at marginal for optimum values of nu####
# in this case kappa = nu = shape parameter (bad choice of variable names)
tn = c(.212,.25,.254,.245)
out_phi <- as.list(1:4)
for (i in 1:4){
  out_phi[[i]] <- sapply(phi.grid, function(phi) {
    loglikeRange(phi=phi, 
                 tau2OverSig2= tn[i],
                 kappa=kappa_list[[i]], y=y, X=X, s=s)
  })
}

for (i in 1:4) {
  plot(phi.grid, out_phi[[i]],type = 'b', xlab='phi', ylab='loglike', 
       main=paste0('nu = ', kappa_list[[i]]))
  abline(v=phi.grid[which.max(out_phi[[i]])], lty=2)
}

