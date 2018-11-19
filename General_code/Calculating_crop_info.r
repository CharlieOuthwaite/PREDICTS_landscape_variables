rm(list = ls())

# load libraries
library(raster)
library(ggplot2)
library(ggthemes)

# where are the datasets?
cropdir <- "D:/BIOTA/Data/HarvestedAreaYield175Crops_Geotiff/HarvestedAreaYield175Crops_Geotiff/HarvestedAreaYield175Crops_Geotiff"


# list the directories in this folder, one for each crop
crop.dirs <- list.dirs(cropdir, recursive = FALSE) # 175 crops


# loop through each crop folder and combine the production columns for each crop into one data frame.
crop <- crop.dirs[1]

# where to put the results
crop.prod <- NULL

for(crop in crop.dirs){
  
  # what is the name of the crop
  crop.name <- sub("_.*", "", sub(".*/", "", crop))
  
  # read in the production file
  
  crop.file <- list.files(crop, pattern = "Production.tif$")
  
  crop.tab <- raster(paste0(crop, "/", crop.file))
  
  # convert to dataframe 
  crop.tab <- as.data.frame(crop.tab, xy = TRUE)
  
  # add crop name
  crop.tab$crop <- crop.name
  
  colnames(crop.tab)[3] <- "production"
  
  # add this info to the all crops table
  crop.prod <- rbind(crop.prod, crop.tab)
  
  
}

# save the file (it takes a long time to produce)
write.csv(crop.prod, file = "D:/BIOTA/1_Forest_Cover_Yield/1. Data/Production/Crop_production_175.csv")

# read in the harvested area info





all_cropdata2 <- dcast(crop.prod, lon + lat ~ crop, value.var = "production")


ggplot()+
  borders("world", colour = "black", size = 0.3) +
  theme_map() +
  geom_raster(data = crop.test.df, aes(x = x, y = y, fill = abaca_Production)) +
  scale_fill_viridis_c() 
#  theme(panel.grid.major = element_line(colour = "transparent"))

