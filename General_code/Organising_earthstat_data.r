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
library(cowplot)


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
  
  # standardise the yield column (if necessary)
  #cropdata$yield_stan <- scale(cropdata$yield, center = T, scale = T)
  
  # add crop name to the dataframe
  cropdata$crop <- sub("_yield.csv", "", file)
  
  # add this info to the all data table
  all_cropdata <- rbind(all_cropdata, cropdata)
  
}


# try to organise into a dif form with cast to see if any cells have more than one crop type

# just keep the standardised yield column
#all_cropdata2 <- all_cropdata[, c("lon", "lat", "crop", "yield_stan")] # if standardising
#all_cropdata2 <- dcast(all_cropdata, lon +lat ~ crop, value.var = "yield_stan") # if standardising
all_cropdata2 <- dcast(all_cropdata, lon +lat ~ crop, value.var = "yield")


# now calculate total yield per cell
# add a column detailing total yield per cell, sum relevent columns
all_cropdata2$totals <- rowSums(all_cropdata2[, 3:ncol(all_cropdata2)], na.rm = TRUE)


# save csv file of the standardised yields and total data
#write.csv(all_cropdata2, paste0("D:/BIOTA/1_Forest_Cover_Yield/Data Exploration/All_Crop_Yield_standardised.csv"))
write.csv(all_cropdata2, paste0("D:/BIOTA/1_Forest_Cover_Yield/Data Exploration/All_Crop_Yield.csv"))

# normality checks on yield totals?
# histogram of totals
p1 <- ggplot(all_cropdata2, aes(totals)) +
  geom_histogram()

# q-q plot of totals
p2 <- ggplot(all_cropdata2, aes(sample = totals)) +
  stat_qq() + 
  stat_qq_line()

# combine into one
plot_grid(p1, p2) # organise with cowplot function

# save the plot
ggsave(filename = paste0("D:/BIOTA/1_Forest_Cover_Yield/Data Exploration/Normality_check_YieldTotal_Standardised.png"))
ggsave(filename = paste0("D:/BIOTA/1_Forest_Cover_Yield/Data Exploration/Normality_check_YieldTotal.png"))


# plot total yield on a global map

library(RColorBrewer)
cols<-brewer.pal(4, "OrRd")

ggplot(all_cropdata2, aes(x = lon, y = lat))+
  borders("world", colour = "gray40", size = 0.3) +
  theme_map() +
  geom_raster(data = all_cropdata2, aes(x = lon, y = lat, fill=totals), inherit.aes = FALSE) +
  #scale_fill_manual(values=cols) +
#  guides(fill=guide_legend(title="Yield \n(standardised)")) +
  guides(fill=guide_legend(title="Yield \n(Tons per hectare)")) +
  labs(title = paste("Total yield (17 crops), Earthstat")) +
  theme(legend.position="bottom")+
  theme_void()+ coord_equal()

# save the plot
#ggsave(filename = "D:/BIOTA/1_Forest_Cover_Yield/Data Exploration/Map_YieldTotals_standardised.png",
#        width = 8, height = 4)
ggsave(filename = "D:/BIOTA/1_Forest_Cover_Yield/Data Exploration/Map_YieldTotals.png",
       width = 8, height = 4)
