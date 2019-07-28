library(ncdf4)
library(maptools)
library(rworldmap)
library(gtools)
library(ggmap)
library(maps)
library(geoR)
library(spBayes)
####################################################################################################################################
# Read In Data #
library(RCurl)
gitstring = "https://raw.githubusercontent.com/msilva00/SpatialStats/master/Homework5/goes135_sample.csv"
goes135 <-read.csv(text=getURL(gitstring))[,c(2:4)]
head(goes135)
gitstr2 = "https://raw.githubusercontent.com/msilva00/SpatialStats/master/Homework5/goes75_sample.csv"
goes75 = read.csv(text=getURL(gitstr2))[,c(2:4)]
head(goes75)

minlon = -108; maxlon = -50; minlat = 34; maxlat = 64

# knots (where to predict)
grid <- as.matrix(expand.grid(seq(minlon, maxlon, l=5), seq(minlat, maxlat, l=5)))

# land only
data(wrld_simpl)
pts <- SpatialPoints(grid, proj4string=CRS(proj4string(wrld_simpl)))
knots <- pts[!is.na(over(pts, wrld_simpl)$FIPS)]@coords

#### plot maps ####
newmap <- getMap(resolution = "low")
plot(newmap, xlim = c(minlon, maxlon), ylim = c(minlat, maxlat))
points(goes75$lon, goes75$lat, col = 'red', pch = 20)
points(knots[,1], knots[,2])

plot(newmap, xlim = c(minlon, maxlon), ylim = c(minlat, maxlat))
axis(side = 1)
axis(side=2)
points(goes135$lon, goes135$lat, col = 'green', pch = 20)
points(knots[,1], knots[,2])


#### EDA ####
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

#### Trend ####
plot(goes75$lon, y1)
plot(goes75$lat, y1)
plot(goes135$lon, y2)
plot(goes135$lat, y2)

loc1 = goes75[,1:2]
loc2 = goes135[,1:2]

geodata1 = list(y1, loc1)
geodata2 = list(y2, loc2)
names(geodata1) = names(geodata2) = c('data', 'coords')

# mod = step(lm(y1 ~ . + .^2 + I(loc1^2), data = data.frame(loc1)), scope = list("lower" = lm(y1 ~ 1), "upper" = lm(y1 ~ . + .^2 + I(loc1^2), data = data.frame(loc1))), direction = "both")
# summary(mod)
# AIC(mod)

mod1 = lm(y1 ~ lon + lat + I(lon^2) + I(lat^2) + lon*lat, data = loc1)
summary(mod1)
AIC(mod1) # 1495.747

mod2 = lm(y2 ~ lon_135 + lat_135 + I(lon_135^2) + I(lat_135^2) + lon_135*lat_135, data = loc2)
summary(mod2)
AIC(mod2) # 620.4768

#### predictive process ####
X = cbind(loc1$lon, loc1$lat, I(loc1$lon^2), I(loc1$lat^2), loc1$lon*loc1$lat)


priors = list("beta.flat", "sigma.sq.ig" = c(1, 0.5), "tau.sq.ig" = c(1, 0.01),
              "phi.Unif" = c(1e-6, 5), "nu.Unif" = c(0.1, 3.5))
starting = list('beta' = mod1$coefficients, "sigma.sq" = 1.317742, "tau.sq" = 0.1721688,
                "phi" = 1/0.6379521, "nu" = 2.5)
tuning = list("sigma.sq" = 0.01, "tau.sq" = 0.01, "phi" = 0.01, "nu" = 0)


pp1 = spLM(y1 ~ X, coords = X[,1:2], knots = knots, 
           n.samples = 4000, cov.model = 'matern', 
           modified.pp = FALSE, priors = priors, 
           starting = starting, tuning = tuning)

par(mfrow = c(3,1))
sigma.sq = pp1$p.theta.sample[,1]
tau.sq = pp1$p.theta.sample[,2]
phi = pp1$p.theta.sample[,3]
traceplot(sigma.sq)
traceplot(tau.sq)
plot.ts(phi)

mean(sigma.sq)
colMeans(pp1$p.theta.sample)


#### predict ####
X.test = cbind(1, test[,1], test[,2], I(test[,1]^2), I(test[,2]^2), test[,1]*test[,2])
pp1.samples = spPredict(pp1, pred.coords = test[,1:2], pred.covars = X.test)
# pp1.mod.samples = spPredict(pp1.mod, pred.coords = goes75[,1:2], pred.covars = X)

pred.pp1 = apply(pp1.samples$p.y.predictive.samples, 1, mean)
# pred.pp1.mod = apply(pp1.mod.samples$p.y.predictive.samples, 1, mean)

#### accuracy ####
pred.int = apply(pp1.samples$p.y.predictive.samples, 1, function(x) quantile(x, c(0.05, 0.95)))

# proportion of observatins that fall inside the prediction interval
1-(sum(log(test[,3]) > pred.int[2,] | log(test[,3]) < pred.int[1,])/nrow(test))

plot(density(pp1.samples$p.y.predictive.samples[1,]))
abline(v = log(test[1,3]))

plot(density(pp1.samples$p.y.predictive.samples[100,]))
abline(v = log(test[100,3]))

plot(density(pp1.samples$p.y.predictive.samples[500,]))
abline(v = log(test[500,3]))

ind=1:25
matplot(rbind(ind,ind), pred.int[,ind], type="l", lty=1,col=1, xlab="Observation",ylab="log albedo", main = '90% Prediction Intervals for Albedo')
points(ind, log(test[1:25,3]), pch=19)

#---------------------------- plot ---------------------------
par(mfrow = c(2,1))
quilt.plot(x=test[,1], y=test[,2],
           z=test[,3], main = 'Observed Albedo')
map("world", add = TRUE)

quilt.plot(x=test[,1], y=test[,2],
           z=exp(pred.pp1), main = 'Predicted Albedo')
map("world", add = TRUE)







