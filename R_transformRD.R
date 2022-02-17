rm(list=ls())
setwd(dir = 'C:/Users/Judith/Documents/Stage_HbR/_Scripts/R')

# Download package Spatial Data
library(sp)

# Read data, file with two columns containing RD coordinates 
r <- read.csv('cells_xy.csv', header=FALSE)
coordinates(r) = c(1, 2)  

# Assign a certain corodinate reference system to the data. 28992 is the code for RD. 
proj4string(r) <- CRS('+init=epsg:28992')

# Coordinate transform to UTM
rUTM <- spTransform(r, CRS('+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'))

# Save file
write.csv(rUTM, file = "cells_xy_UTM.csv")

# Coordinate transform to WGS84
# rWGS <- spTransform(r, CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'))

# Save file
# write.csv(rWGS, file = "adcp1_WGS.csv")


