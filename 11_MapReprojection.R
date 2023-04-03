#Mapping the trap locations 

#Date: 2022-01-20
#Authors: Gonzalo + Iona 

# Empty global environment 
rm(list=ls())

#Load libraries 
library(raster)
library(rgeos)
library(rgdal) # retired
library(maptools) #retired
library(sp)

# 1. Load the spatial object ####
# Using the UK boundaries
country_boundary_uk_UTM <- readOGR(dsn = "Data/13_NUTS_Level_1_(January_2018)_Boundaries",
                                   layer = "NUTS_Level_1_(January_2018)_Boundaries")

# Look at the CRS 
proj4string(country_boundary_uk_UTM) #In UTM

# Plot map
par(mfrow=c(1,2))
plot(country_boundary_uk_UTM , main="UTM")
axis(1) ; axis(2)

# 2. Choose the right projection ####

# From:
# https://spatialreference.org/ref/epsg/4326/proj4/
# +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs 
#+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs 


# Apply the new projection 
country_boundary_uk_latlon <- spTransform(country_boundary_uk_UTM, 
                                          CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

proj4string(country_boundary_uk_latlon)

#Plot
plot(country_boundary_uk_latlon,main="LatLon")
axis(1)
axis(2)



