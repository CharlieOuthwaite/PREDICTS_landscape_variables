rm(list = ls())

# load libraries
library(raster)
library(ggplot2)
library(ggthemes)
library(reshape2)

# where are the datasets?
cropdir <- "D:/BIOTA/Data/HarvestedAreaYield175Crops_Geotiff/HarvestedAreaYield175Crops_Geotiff/HarvestedAreaYield175Crops_Geotiff"


# list the directories in this folder, one for each crop
crop.dirs <- list.dirs(cropdir, recursive = FALSE) # 175 crops


# loop through each crop folder and combine the production columns for each crop into one data frame.
#crop <- crop.dirs[1]

# where to put the results
all.crop <- NULL

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
  
  
  # reoganise
  crop.tab <- dcast(crop.tab, x + y ~ crop, value.var = "production")
  
  
  if(crop == crop.dirs[1]){all.crop <- crop.tab}else{
  
  if(identical(crop.tab[, 1:2], all.crop[, 1:2]) == TRUE){
    
    all.crop <- cbind(all.crop, crop.tab[, 3])
    
    colnames(all.crop)[grep(",", colnames(all.crop))] <- crop.name
  }else{stop}
}
  
}

# save the file (it takes a long time to produce)
write.csv(all.crop, file = "D:/BIOTA/1_Forest_Cover_Yield/1. Data/Production/Crop_production_175.csv", row.names = F)



##### read in the harvested area info  #####

# loop through each crop folder and combine the production columns for each crop into one data frame.
#crop <- crop.dirs[1]

# where to put the results
all.harv <- NULL

for(crop in crop.dirs){
  
  # what is the name of the crop
  crop.name <- sub("_.*", "", sub(".*/", "", crop))
  
  # read in the production file
  
  harv.file <- list.files(crop, pattern = "HarvestedAreaHectares.tif$")
  
  harv.tab <- raster(paste0(crop, "/", harv.file))
  
  # convert to dataframe 
  harv.tab <- as.data.frame(harv.tab, xy = TRUE)
  
  # add crop name
  harv.tab$crop <- crop.name
  
  colnames(harv.tab)[3] <- "HarvestedAreaHectares"
  
  
  # reoganise
  harv.tab <- dcast(harv.tab, x + y ~ crop, value.var = "HarvestedAreaHectares")
  
  
  if(crop == crop.dirs[1]){all.harv <- harv.tab}else{
    
    if(identical(harv.tab[, 1:2], all.harv[, 1:2]) == TRUE){
      
      all.harv <- cbind(all.harv, harv.tab[, 3])
      
      colnames(all.harv)[grep(",", colnames(all.harv))] <- crop.name
    }else{stop}
  }
  
}

# save the file (it takes a long time to produce)
write.csv(all.harv, file = "D:/BIOTA/1_Forest_Cover_Yield/1. Data/Production/Crop_HarvestArea_hectares_175.csv", row.names = F)


############################################################
#                                                          #
#             Calculate the average production             #
#               per hectare across all crops               #
#                                                          #
############################################################

# as discussed in the group meeting, use two metrics derived from yield:
  # 1. Measure of intensity (average production per hectare)
  # 2. Landscape state (total production per grid cell) - this one is already an available dataset


### calculate the production







ggplot()+
  borders("world", colour = "black", size = 0.3) +
  theme_map() +
  geom_raster(data = crop.test.df, aes(x = x, y = y, fill = abaca_Production)) +
  scale_fill_viridis_c() 
#  theme(panel.grid.major = element_line(colour = "transparent"))




