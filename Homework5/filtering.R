# # FULL CODE
setwd("~/homework3-spacialstats/")
rm(list = ls())
libs = c('geoR','sf','fields', 'ncdf4','tidyverse')
sapply(libs, require, character.only = TRUE)
# 

#### Filter the data ####
# dissolve world_countries_boundaries to yield landmass
land = read_sf(
  dsn = "UIA_World_Countries_Boundaries",
  layer = "UIA_World_Countries_Boundaries"
) %>%
  mutate(bleh = 'bleh') %>%
  group_by(bleh) %>%
  summarise()

goes1_path = 'GSA_AlbedoProd_GOES_075_VIS02_2000_181.nc'
goes2_path = 'GSA_AlbedoProd_GOES_135_VIS02_2000_181.nc'


get_goes_data = function(minlat, maxlat, minlon, maxlon){
  # subsets satellite data to geography in question; 
  #  looking specifically at observations on land.
  
  # subset land data to area in question
  # build bounding box to subset the landmass data
  box_raw = list(
    rbind(c(minlon, minlat), c(minlon, maxlat), 
          c(maxlon, maxlat), c(maxlon, minlat), 
          c(minlon, minlat)
      )
    )
  box_poly = st_polygon(box_raw, dim = 'XY')
  box = st_sfc(
    list(box_poly), 
    crs = "+proj=longlat +datum=WGS84 +no_defs"
    ) %>% st_sf
  # intersect landmass data with bounding box
  land_in_box = st_intersection(land, box)
  
  # open satellite data, pull and filter data to yield points inside bounding box
  goes1 = nc_open(goes1_path)
  data1 = tibble(
    lat = as.vector(ncvar_get(goes1, 'latitude')),
    lon = as.vector(ncvar_get(goes1, 'longitude')),
    albedo = as.vector(ncvar_get(goes1, 'BHRiso'))
  ) %>% filter(
      !is.na(albedo) &
      between(lon, minlon, maxlon) & 
      between(lat, minlat, maxlat)
  )
  goes2 = nc_open(goes2_path)
  data2 = tibble(
    lat = as.vector(ncvar_get(goes2, 'latitude')),
    lon = as.vector(ncvar_get(goes2, 'longitude')),
    albedo = as.vector(ncvar_get(goes2, 'BHRiso'))
  ) %>% filter(
    !is.na(albedo) &
      between(lon, minlon, maxlon) & 
      between(lat, minlat, maxlat)
  )
  data = rbind(data1,data2)
  
  # create points from coordinates; intersect with land in box
  points_in_space = data %>% 
    st_as_sf(
      coords = c('lon','lat'), 
      crs = "+proj=longlat +datum=WGS84 +no_defs"
    )
  ol = st_intersection(points_in_space, land_in_box)[,c('geometry','albedo')]
  tab_of_ol = cbind(st_coordinates(ol), ol$albedo) %>% as_tibble
  colnames(tab_of_ol) = c('longitude','latitude','albedo')
  return(list(table = tab_of_ol, poly = land_in_box))
}

minlon = -108; maxlon = -50; minlat = 34; maxlat = 64

mydata = get_goes_data(minlat, maxlat, minlon, maxlon)
setwd("~/SpatialStats/Homework5/")
mydata$table %>% as_tibble %>% write_csv('subsetted_data.csv')
mydata$poly %>% st_write('land', 'land.shp', driver = 'ESRI Shapefile', delete_dsn = TRUE, delete_layer = TRUE)

# EOF