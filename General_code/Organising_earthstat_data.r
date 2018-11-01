############################################################
#                                                          #
#      Playing around with Monica's EarthStat dataset      #
#                                                          #
############################################################

# started 31/10/2018

# currently using the Earthstat data processed by Monica, not sure what she did to get
# these csv files.  No units given. 

rm(list = ls())

# load libraries
library(ggplot2)
library(ggthemes)
library(reshape2)


# where is the yield data saved
datadir <- "D:/BIOTA/1_Forest_Cover_Yield/1. Data/Yield"

# list files in datadir
files <- list.files(datadir, pattern = ".csv")

# read in each yield data file (17 crops), standardise, then add to a new table
# including crop name

# space to save standardised data
all_cropdata <- NULL


# file <- files[1]

for(file in files){
  
  # read in the csv file
  cropdata <- read.csv(paste0(datadir, "/", file))
  
  # subset to desired columns
  cropdata <- cropdata[, c('lon', 'lat', 'yield')]
  
  # standardise the yield column
  cropdata$yield_stan <- scale(cropdata$yield, center = T, scale = T)
  
  # add crop name to the dataframe
  cropdata$crop <- sub("_yield.csv", "", file)
  
  # add this info to the all data table
  all_cropdata <- rbind(all_cropdata, cropdata)
  
}


# try to organise into a dif form with cast to see if any cells have more than one crop type

# just keep the standardised yield column
all_cropdata2 <- all_cropdata[, c("lon", "lat", "crop", "yield_stan")]
all_cropdata2 <- dcast(all_cropdata, lon +lat ~ crop, value.var = "yield_stan")





# load in one of the crops
cropdata <- read.csv(paste0(datadir, "/", "barley_yield.csv"))

### data is cufrrently in point format, need to create a polygon layer?



ggplot(cropdata, aes(x = lon, y = lat)) +
  borders("world", colour = "gray40", fill = "gray75", size = 0.3) +
  theme_map() +
  geom_point(shape = 21, colour = "black", fill = "blue", alpha = 0.8, size = 3)


ggplot(cropdata, aes(x = lon, y = lat))+
  borders("world", colour = "gray40", fill = "gray75", size = 0.3) +
  theme_map() +
  geom_polygon(fill = NA, colour = "black", size=0.2)+
  geom_raster(data =cropdata, aes(x = lon, y = lat, fill=yield), inherit.aes = FALSE) +
  scale_color_continuous(value = )

