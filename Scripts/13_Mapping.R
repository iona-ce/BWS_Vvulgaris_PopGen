# Adding colours to the map 

# 1. Map data preparation ####

# Prepare the data: assign numerical values to different regions 
# to use later on when assigning colours to the map 

# Load libraries 
library(plyr)
library(broom)

# Getting ID names from the NUTS layers 
id <- country_boundary_uk_latlon$nuts118nm
id

# Assign values to these numbers, i.e. group by region 
value <- c(6, #North East (NE)
           6, #North West (NE)
           6, #Yorkshire and The Humber (NE)
           5, #East Midlands (EE)
           4, #West Midlands (WE)
           5, #East of England (EE)
           2, #London (SE)
           2, #South East (SE)
           2, #South West (SE)
           4, #Wales (WE)
           4, #Scotland (Scotland)
           1) #Northern Ireland

# Create data frame 
mydata <- data.frame(id, value)
str(mydata)

# Change to factor 
mydata$value <- as.factor(mydata$value)

# Tidy data 
mapdata <- tidy(country_boundary_uk_latlon, region = "nuts118nm") 

#Create dataframe 
df <- join(mapdata, mydata, by="id")

# 2. Plotting the map ####

#Load libraries 
library(ggplot2)
library(ade4)
library(adegenet)
library(ade4)
library(ggthemes) 
library(RColorBrewer)
library(scales)

# Import genind object with 2017 data (if not already imported) 
# or any dataframe with lat long data 
# bws17_genind 
load("Data/12_GenindObjects/12_BWS_50_2017_National_Genind.RData")

# Create map in black and white
# linewidth was originally size - seems to have changed under new ggplot version 
gg <-  ggplot() + 
  geom_polygon(data = df, aes(x = long, y = lat, group = group, fill = value), 
               color = "#FFFFFF", linewidth = 0) +
  scale_fill_grey(start = 0.45, end = 0.85) + 
  geom_point(data = bws17_genind@other$xy, aes(x = x, y = y), 
             col = "black", size = 1) +
  theme_clean() +
  coord_fixed(1.25)

plot(gg)

# 3. Optional: Colour map ####
# Optional and not used in manuscript:
# Some code to do the same with colours (package viridis)

# Load package 
library(viridis)

# Choose colours 
# scales::show_col(viridis(option = "plasma", n=100))
# cbPalette <- c(viridis(option = "plasma", n=81)[c(25,35,54,63,75,45)]) #Make sure that the n here is the same as above

scales::show_col(viridis(option = "magma", n=100))
cbPalette <- c(viridis(option = "magma", n=100)[c(87, 84, 90, 97, 93, 100)]) #Make sure that the n here is the same as above

show_col(cbPalette)

# Create map 
ggc <-  ggplot() + 
  geom_polygon(data = df, aes(x = long, y = lat, group = group, fill = value), 
               color = "#FFFFFF", size = 0) +
  scale_fill_manual(values=cbPalette) + 
  coord_fixed(1) + 
  theme_minimal() + 
  geom_point(data = bws17_genind@other$xy, aes(x = x, y = y), 
             col = "black", size = 1)


plot(ggc)

# 4. BWS 18 Cluster map ####

# Using the same colour scheme in the map to have an idea of the locations of the 
# 2018 traps relative to the 2017 ones 

# Load 2018 data
# bws18_genind
load("Data/12_GenindObjects/12_BWS_50_WithSibs_2018_Regional_Genind.RData")

# Check that has loaded correctly 
head(bws18_genind_regional@other$xy)

# First aim: calculate the centre of each cluster 
# This is necessary as otherwise it creates an ugly map 

# Trasnform genind information into a table 
bws18_geo_info <- as.data.frame(bws18_genind_regional@other$xy)
head(bws18_geo_info)

# Add in population (cluster) data
bws18_geo_info$pop <- bws18_genind_regional@pop

#Rename columns 

# Create a list for each cluster 
bws18_geo_info_cluster <- split(bws18_geo_info, f = bws18_geo_info$pop)  

#Calculate mean Lat and Long - probably an easier way to do this 
Lat17 <- mean(bws18_geo_info_cluster$`17`$x)
Lat22 <- mean(bws18_geo_info_cluster$`22`$x)
Lat29 <- mean(bws18_geo_info_cluster$`29`$x)
Lat62 <- mean(bws18_geo_info_cluster$`62`$x)
Lat122 <- mean(bws18_geo_info_cluster$`122`$x)
Lat193 <- mean(bws18_geo_info_cluster$`193`$x)
Lat205 <- mean(bws18_geo_info_cluster$`205`$x)
Lat286 <- mean(bws18_geo_info_cluster$`286`$x)
Lat328 <- mean(bws18_geo_info_cluster$`328`$x)

#Mean long 
Long17 <- mean(bws18_geo_info_cluster$`17`$y)
Long22 <- mean(bws18_geo_info_cluster$`22`$y)
Long29 <- mean(bws18_geo_info_cluster$`29`$y)
Long62 <- mean(bws18_geo_info_cluster$`62`$y)
Long122 <- mean(bws18_geo_info_cluster$`122`$y)
Long193 <- mean(bws18_geo_info_cluster$`193`$y)
Long205 <- mean(bws18_geo_info_cluster$`205`$y)
Long286 <- mean(bws18_geo_info_cluster$`286`$y)
Long328 <- mean(bws18_geo_info_cluster$`328`$y)

#Combine these into vectors 
Lat <- Lat <- c(Lat17, Lat22, Lat29, Lat62, Lat122, Lat193, Lat205, Lat286, Lat328)
Long <- c(Long17, Long22, Long29, Long62, Long122, Long193, Long205, Long286, Long328)

#Create data frame 
dist_table <- cbind(Long, Lat)
dist_table <- as.data.frame(dist_table)
dist_table$pop <- unique(bws18_genind$pop)

#Map these points 

gg1 <-  ggplot() + 
  geom_polygon(data = df, aes(x = long, y = lat, group = group, fill = value), 
               color = "#FFFFFF", linewidth = 0) +
  scale_fill_grey(start = 0.45, end = 0.85) + 
  geom_point(data = dist_table, aes(x = Lat, y = Long), 
             col = "black", size = 1) +
  theme_clean() +
  coord_fixed(ratio = 1.25, xlim = c(-3,1.75), ylim = c(50.5, 52)) 

plot(gg1)

gg1 <-  ggplot() + 
  geom_polygon(data = df, aes(x = long, y = lat, group = group, fill = value), 
               color = "#FFFFFF", linewidth = 0.05) +
  scale_fill_manual(values=cbPalette) + 
  geom_point(data = dist_table, aes(x = Lat, y = Long), 
             col = "black", size = 2) +
  theme_clean() +
  coord_fixed(ratio = 1.25, xlim = c(-3,1.75), ylim = c(50.5, 52)) 

plot(gg1)
