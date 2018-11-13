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

# where are the datasets?
datadir <- "D:/BIOTA/Data"

# load in the biome dataset
library(sf)

# read in the biomes shapefile (this has 4 biomes in, just for forest types?)
# humid tropics, dry tropics, temperate, boreal
biomes <- st_read(paste0(datadir, "/biomes/biomes.shp")) %>%
  st_set_crs(4326)

# basic plot of the biomes
plot(biomes['BIOME'])


# where are the PREDICTS cropland points
pred <- read.csv(paste0(datadir, "/PREDICTS_2016/Cropland_Sites_SpRichness.csv"))


# try to extract biome for each point
# extract yield values from the predicts points

# convert pred points into an sf object?
pred.sf <- st_multipoint(as.matrix(pred[, c("Longitude", "Latitude")]), dim = "XY")

# add in a CRS?
pred.sf <- st_sfc(pred.sf, crs = 4326)


library(ggplot2)
library(ggthemes)

ggplot()+
  borders("world", colour = "black", size = 0.3) +
  theme_map() +
  geom_sf(data = biomes, aes(fill = BIOME)) +
  geom_sf(data = pred.sf, col = "black") +
  theme(panel.grid.major = element_line(colour = "transparent"))



# which of the predicts points are in which biome

pred$biome <- st_intersects(biomes['BIOME'], pred.sf)

# remove sites that are not within a forest biome

