##%######################################################%##
#                                                          #
####          Visualising PREDICTS data points          ####
#                                                          #
##%######################################################%##

# playing around with making a map of PREDICTS data points

rm(list = ls())

# load libraries
install.packages("maps")
install.packages("mapdata")

library(maps)
library(mapdata)

# where is the data?
datadir <- "C:/Users/chaout/Dropbox/POSTDOC - biodiv & food sec/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_2016"

# load data
PRED.data <- read.csv(paste0(datadir, "/", "resource.csv")) # takes a long time to read in!


# subset to cropland sites only
crop.data <- PRED.data[PRED.data$Predominant_land_use == "Cropland", ]

# number of studies
length(unique(crop.data$Study_number)) # 21 studies

# number of sites
length(unique(crop.data$Site_name)) # 1981 sites

#number of species
length(unique(crop.data$Best_guess_binomial)) # 6947

# number of countries
length(unique(crop.data$Country)) # 44

# number of ecoregions
length(unique(crop.data$Ecoregion)) # 71

# number of biomes
length(unique(crop.data$Biome)) # 9

table(as.character(crop.data$Higher_taxon))

map("worldHires", 
    col = "grey90",
    fill = TRUE)

points(crop.data$Longitude,
       crop.data$Latitude,
       col = "blue", 
       cex = 1,
       pch = 19)



### trying ggplot mapping ###

library(ggthemes)
library(ggplot2)

ggplot(crop.data, aes(x = Longitude, y = Latitude)) +
  borders("world", colour = "gray40", fill = "gray75", size = 0.3) +
  theme_map() +
  geom_point(shape = 21, colour = "black", fill = "blue", alpha = 0.8, size = 3)
  
