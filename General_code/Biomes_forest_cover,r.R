############################################################
#                                                          #
#        Playing with biomes and forest cover data         #
#                                                          #
############################################################

# 12th November 2018

# Using data downloaded from https://glad.umd.edu/dataset/gfm/globaldata/global-data

# biome dataset and forest cover dataset from the above website, figuring out how to read
# them in and use with the PREDICTS data, once done, add into a main script/workflow

rm(list = ls())

# load libraries
library(sf)
library(ggplot2)
library(ggthemes)

# where are the datasets?
datadir <- "D:/BIOTA/Data"

# load in the biome dataset

# read in the biomes shapefile (this has 4 biomes in, just for forest types?)
# humid tropics, dry tropics, temperate, boreal

# setting the CRS to world lat/longs
biomes <- st_read(paste0(datadir, "/biomes/biomes.shp")) %>%
  st_set_crs(4326)

# removing the extra grid code column
biomes <- biomes[2:3]

# basic plot of the biomes
plot(biomes)



# where are the PREDICTS cropland points
pred <- read.csv(paste0(datadir, "/PREDICTS_2016/Cropland_Sites_SpRichness.csv"))


# try to extract biome for each point
# extract yield values from the predicts points

# convert pred points into an sf object?
#pred.sf <- st_multipoint(as.matrix(pred[, c("Longitude", "Latitude")]), dim = "XY")

# add in a CRS?
#pred.sf <- st_sf(pred.sf, crs = 4326)



ggplot()+
  borders("world", colour = "black", size = 0.3) +
  theme_map() +
  geom_sf(data = biomes, aes(fill = BIOME), alpha = 0.5) +
  geom_sf(data = pred.sf, col = "black") +
  theme(panel.grid.major = element_line(colour = "transparent"))



# which of the predicts points are in which biome

# had no idea how to do this, took ages to find an example that I understood:
# https://gis.stackexchange.com/questions/282750/identify-polygon-containing-point-with-r-sf-package

pred$biome <- apply(pred[ , c("Longitude", "Latitude")], 1, function(x) {  
  # transformation to palnar is required, since sf library assumes planar projection 
 
  pred.sf <- st_transform(st_sfc(st_point(x),crs = 4326), 4326)
  # st_intersects with sparse = FALSE returns a logical matrix
  # with rows corresponds to argument 1 (points) and 
  # columns to argument 2 (polygons)
  
  biomes[which(st_intersects(pred.sf, biomes, sparse = FALSE)), ]$BIOME 
})

# rename biome values from numbers
levels(biomes$BIOME)
pred$biome <- sub(1, "Boreal", pred$biome)
pred$biome <- sub(2, "Dry Tropics", pred$biome)
pred$biome <- sub(3, "Humid Tropics", pred$biome)
pred$biome <- sub(4, "Temperate", pred$biome)


ggplot()+
  borders("world", colour = "black", size = 0.3) +
  theme_map() +
  geom_sf(data = biomes, aes(fill = BIOME), alpha = 0.5) +
  geom_point(data = pred, aes(x = Longitude, y = Latitude, col = biome)) +
  theme(panel.grid.major = element_line(colour = "transparent"))


# remove sites that are not within a forest biome
nrow(pred[pred$biome == "integer(0)",])
