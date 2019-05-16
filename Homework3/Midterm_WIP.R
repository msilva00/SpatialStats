# # FULL CODE
# setwd("~/SpatialStats/Homework3/")
# rm(list = ls())
# libs = c('geoR','sf','fields', 'ncdf4','tidyverse')
# sapply(libs, require, character.only = TRUE)
# 

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
random5001 = sample_n(ex, 500)
random500 = random5001[,c(1,2,6)]
y = bcPower(random500[,3],0)
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
title(main = "Albedo Measurement \n (Box and Cox Transformation)",cex.main = 1.5, outer = T, line = -3)
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
par(mfrow = c(2,2), mar = c(3,3,5,3))
hist(log(random500[,3]), freq = F, main = "")
lines(density(log(random500[,3])))

qqnorm(log(random500[,3]), main = "")
qqline(log(random500[,3]))

title("Box and Cox Transformed Weighted Albedo Measurements", outer = T, line = -2)


hist((random500[,3]), freq = F, main = "", ylim = c(0,24), breaks = 20)
lines(density((random500[,3])))

qqnorm((random500[,3]), main = "")
qqline((random500[,3]))

title("Weighted Albedo Measurements (w/out transformation)", outer = T, line = -20)

#### Residual Analysis/Covariate Analysis ####
head(random500)
##### LM #####
y = log(random500$alb_weight)
lon = random500$longitude
lat = random500$latitude
locs = cbind(lon, lat)
s = cbind(lon,lat)

mod = step(lm(y ~ . + .^2 + I(loc^2), data = data.frame(loc)), scope = list("lower" = lm(y ~ 1),
                                                                            "upper" = lm(y ~ . + .^2 + I(loc^2), data = data.frame(loc))),
           direction = "both")
summary(mod)



mod = lm(y~ lon + lat + lon*lat)
summary(mod)


geodata = list(y, locs)
names(geodata) = c("data", "coords")

MC <- model.control(trend.d = "2nd", trend.l = "2nd",kappa = 1.5)
PC <- prior.control(phi.discrete = seq(0, 6, l = 21),
                    phi.prior = "reciprocal", tausq.rel.prior = "unif",
                    tausq.rel.discrete = seq(0, 1, l = 11))
OC <- output.control(n.post = 1000, moments = T)

?output.control
skb <- krige.bayes(geodata, model = MC, prior = PC, output = OC, locations=locs)
hist(skb$posterior$sample$beta0, xlab = expression(beta[0]))
abline(v = mean(skb$posterior$sample$beta0), col = "red", lwd = 2)
