library(ncdf4)
library(maptools)
library(rworldmap)
library(gtools)
library(ggmap)
library(maps)
library(geoR)
library(spBayes)

setwd("~/homework3-spacialstats/")
nc1_075 = 'GSA_AlbedoProd_GOES_075_VIS02_2000_181.nc'
nc2_135 = 'GSA_AlbedoProd_GOES_135_VIS02_2000_181.nc'

# get data
nc_075 = nc_open(nc1_075)
nc_135 = nc_open(nc2_135)

minlon = -108; maxlon = -50; minlat = 34; maxlat = 64

#### GOES 075 ####
lon = as.vector(ncvar_get(nc_075, 'longitude'))
lat = as.vector(ncvar_get(nc_075, 'latitude'))
albedo = as.vector(ncvar_get(nc_075, 'BHRiso'))
temp = data.frame(cbind(lon, lat, albedo))
temp = temp[complete.cases(temp),]
# subset my region
temp = temp[temp$lon >= minlon & 
              temp$lon <= maxlon & 
              temp$lat >= minlat & 
              temp$lat <= maxlat, ]
## Find which points fall over land
data(wrld_simpl)
pts <- SpatialPoints(temp, proj4string=CRS(proj4string(wrld_simpl)))
goes <- pts[!is.na(over(pts, wrld_simpl)$FIPS)]@coords
goes75 = goes[sample(nrow(goes), size = 1000),]
goes75 = data.frame(goes75)
head(goes75)
# plot map goes075
newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(minlon, maxlon), ylim = c(minlat, maxlat))
points(goes75[,1], goes75[,2], col = 'red', pch = 20)
setwd("~/SpatialStats/Homework5/")
write.csv(goes75, file="goes75_sample.csv")


#### GOES 135 ####
lon_135 = as.vector(ncvar_get(nc_135, 'longitude'))
lat_135 = as.vector(ncvar_get(nc_135, 'latitude'))
albedo_135 = as.vector(ncvar_get(nc_135, 'BHRiso'))
temp_135 = data.frame(cbind(lon_135, lat_135, albedo_135))
temp_135 = temp_135[complete.cases(temp_135),]
# subset my region
temp_135 = temp_135[temp_135$lon_135 >= minlon & 
                      temp_135$lon_135 <= maxlon & 
                      temp_135$lat_135 >= minlat & 
                      temp_135$lat_135 <= maxlat, ]
summary(temp_135)
## Find which points fall over land
data(wrld_simpl)
pts <- SpatialPoints(temp_135, proj4string=CRS(proj4string(wrld_simpl)))
goes <- pts[!is.na(over(pts, wrld_simpl)$FIPS)]@coords
goes135 = goes[sample(nrow(goes), size = 1000),]
goes135 = data.frame(goes135)
write.csv(goes135, file = "goes135_sample.csv")

# plot goes135
plot(newmap, xlim = c(minlon, maxlon), ylim = c(minlat, maxlat))
axis(side = 1)
abline(v = minlon)
abline(v = maxlon)
axis(side=2)
abline(h=minlat)
abline(h=maxlat)
points(goes135$lon_135, goes135$lat, col = 'green', pch = 20)
points(knots[,1], knots[,2])


# knots (where to predict)
grid <- as.matrix(expand.grid(seq(minlon, maxlon, l=30), seq(minlat, maxlat, l=30)))

# land only
pts <- SpatialPoints(grid, proj4string=CRS(proj4string(wrld_simpl)))
knots <- pts[!is.na(over(pts, wrld_simpl)$FIPS)]@coords


#---------------------------------------------------------------------------
# explore data

hist(goes75$albedo)
hist(goes135$albedo)
# skewed - try logit because albedo is in (0,1)
hist(log(goes75$albedo))
hist(log(goes135$albedo))

plot(density(log(goes75$albedo)))
plot(density(log(goes135$albedo)))

qqnorm(log(goes75$albedo))
qqline(log(goes75$albedo), col = 'red')

qqnorm(log(goes135$albedo))
qqline(log(goes135$albedo), col = 'red')

y1 = log(goes75$albedo)
y2 = log(goes135$albedo)

plot(goes75)
plot(goes135)

library(plotrix)
newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(minlon, maxlon), ylim = c(minlat, maxlat))
points(goes75$lon, goes75$lat, col = 'red', pch = 20)
points(knots[,1], knots[,2])

plot(newmap, xlim = c(minlon, maxlon), ylim = c(minlat, maxlat))
axis(side = 1)
abline(v = minlon)
abline(v = maxlon)
axis(side=2)
abline(h=minlat)
abline(h=maxlat)
points(goes135$lon_135, goes135$lat, col = 'green', pch = 20)
points(knots[,1], knots[,2])

head(goes135)
#------------------------ trend ------------------------

plot(goes75$lon, y1)
plot(goes75$lat, y1)
plot(goes135$lon, y2)
plot(goes135$lat, y2)

loc1 = goes75[,1:2]
loc2 = goes135[,1:2]

geodata1 = list(y1, loc1)
geodata2 = list(y2, loc2)
names(geodata1) = names(geodata2) = c('data', 'coords')

mod = step(lm(y1 ~ . + .^2 + I(loc1^2), data = data.frame(loc1)), scope = list("lower" = lm(y1 ~ 1), "upper" = lm(y1 ~ . + .^2 + I(loc1^2), data = data.frame(loc1))), direction = "both")
summary(mod)
AIC(mod)

mod1 = lm(y1 ~ lon + lat + I(lon^2) + I(lat^2) + lon*lat, data = loc1)
summary(mod1)
AIC(mod1) # 1348.827

mod2 = lm(y2 ~ lon + lat + I(lon^2) + I(lat^2) + lon*lat, data = loc2)
summary(mod2)
AIC(mod2) # 1193.431

mle1 = 
summary(mle)

#---------------------- predictive process ----------------
X = cbind(1, loc1[,1], loc1[,2], loc1[,1]^2, loc1[,2]^2, loc1[,1]*loc1[,2])

pri = list('beta.flat', 'sigma.sq.ig' = c(1, 0.05), 'tau.sq.ig' = c(1, 0.01), 'phi.Unif' = c(1e-6, 100), 'nu.Unif' = c(0.5, 2.5))

start = list('beta' = mod1$coefficients, 'sigma.sq' = 0.3, 'tau.sq' = 0.02, 'phi' = 1.2, 'nu' = 0.5)

tune = list("sigma.sq" = 0.01, "tau.sq" = 0.01, "phi" = 0.01, "nu" = 0.01)

pp1 = spLM(y1 ~ lon + lat + I(lon^2) + I(lat^2) + lon*lat, data = loc1, coords = as.matrix(loc1), knots = knots, n.samples = 500, cov.model = 'matern', modified.pp = FALSE, priors = pri, starting = start, tuning = tune)

pp1.samples = spPredict(pp1, pred.coords = goes75[,1:2], pred.covars = X)

pred.pp1 = apply(pp1.samples$p.y.predictive.samples, 1, mean)
#--------plot---------
plot.nice(loc1, pred.pp1, zscale = range(pp1),
          nlevels = 20, main = "Ozone - Predictive Process")
