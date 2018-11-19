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
all_cropdata2 <- dcast(all_cropdata, lon + lat ~ crop, value.var = "yield")


# now calculate total yield per cell
# add a column detailing total yield per cell, sum relevent columns
all_cropdata2$totals <- rowSums(all_cropdata2[, 3:ncol(all_cropdata2)], na.rm = TRUE)


# save csv file of the standardised yields and total data
#write.csv(all_cropdata2, paste0("D:/BIOTA/1_Forest_Cover_Yield/Data Exploration/All_Crop_Yield_standardised.csv"))
write.csv(all_cropdata2, paste0("D:/BIOTA/1_Forest_Cover_Yield/Data Exploration/All_Crop_Yield.csv"), row.names = F)

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



############################################################
#                                                          #
#                  Now to try to match up                  #
#          the yield data with the PREDICTS data           #
#                                                          #
############################################################

# read in the PREDICTS data from the NHM webpage
pred.all <- read.csv("D:/BIOTA/PREDICTS_2016/resource.csv")



# subset to cropland sites only
crop.data <- pred.all[pred.all$Predominant_land_use == "Cropland", ]


### to start with, just try to look at species richness per site ###

# subset to useful columns
crop.data <- crop.data[, c(6:10, 19:21, 27, 31:34, 38:39, 65:68)]

# proportion of records that are non-zero
1-nrow(crop.data[crop.data$Measurement == 0,])/nrow(crop.data)


#### summarise to site level information on species richness ####

  # where to save the outputs
  all.res <- NULL
  
# for each SSBS get the main info and the species richness of the site, save in new table
  for(site in unique(x$SSBS)){
    
    # subset to one site
    sub.data <- x[x$SSBS == site, ]
    
    # extract information and the number of species
    result <- unlist(c(sub.data[1, c(1, 6:8, 11:15)], nrow(sub.data)))
    
    # add into table
    all.res <- rbind(all.res, result)
    
  }
  
# add column name for species richness stuff
  colnames(all.res)[10] <- "Sp.richness"
  
  # convert to dataframe
  
  all.res <- as.data.frame(all.res)

# save this dataset
  write.csv(all.res, file = "D:/BIOTA/PREDICTS_2016/Cropland_Sites_SpRichness.csv", row.names = F)



############################################################
#                                                          #
#          Try to match up points with yield data          #
#                                                          #
############################################################


# read in yield data 
yield <- read.csv("D:/BIOTA/1_Forest_Cover_Yield/Data Exploration/All_Crop_Yield.csv")

# for some reason, reading in the data means it won't plot (?!)  
# yield <- all_cropdata2

# read in PREDICTS cropland data
pred <- read.csv("D:/BIOTA/PREDICTS_2016/Cropland_Sites_SpRichness.csv")


library(RColorBrewer)
library(ggplot2)
library(ggthemes)

cols<-brewer.pal(4, "OrRd")

ggplot(data = yield, aes(x = lon, y = lat))+
  borders("world", colour = "gray40", size = 0.3) +
  theme_map() +
  geom_raster(data = yield, aes(x = lon, y = lat, fill=totals), inherit.aes = FALSE) +
  geom_point(data = pred, aes( x = Longitude, y = Latitude), col = "green") +
  #scale_fill_manual(values=cols) +
  #  guides(fill=guide_legend(title="Yield \n(standardised)")) +
  guides(fill=guide_legend(title="Yield \n(Tons per hectare)")) +
  labs(title = paste("Total yield (17 crops), Earthstat")) +
  theme(legend.position="bottom")+
  theme_void()+ coord_equal()

ggsave(file = paste0("D:/BIOTA/1_Forest_Cover_Yield/Data Exploration/Overlay_yield_PREDICTS.png"))


### figure out the yield values at the point coordinates ###

# first, convert both to spatial data objects with the same CRS

library(sf)

# remove the rows that have NAs in the lat/long columns
pred <- pred[!is.na(pred$Latitude), ]


library(raster)

# convert the lat/long/value data into a raster
yield.sp <- rasterFromXYZ(yield[, c(1,2,20)])  # add a buffer to reduce NAs?

# extract yield values from the predicts points
yield.vals <- extract(yield.sp, pred[, c("Longitude","Latitude")])

# there are some NAs, how many?
sum(is.na(yield.vals)) # 109 NAs


# reorganise data
pred$total.yield <- yield.vals


# normality checks on yield totals?
# histogram of totals
p1 <- ggplot(pred, aes(total.yield)) +
  geom_histogram()

# q-q plot of totals
p2 <- ggplot(pred, aes(sample = log(total.yield))) +
  stat_qq() + 
  stat_qq_line()

library(cowplot)
# combine into one
plot_grid(p1, p2) # organise with cowplot function

# save the plot
#ggsave(filename = paste0("D:/BIOTA/1_Forest_Cover_Yield/Data Exploration/Normality_check_YieldTotal_Standardised.png"))
ggsave(filename = paste0("D:/BIOTA/1_Forest_Cover_Yield/Data Exploration/Normality_check_YieldTotal.png"))


library(lme4)

random1 <- lmer(Sp.richness~total.yield+(1|SS), data=pred)

summary(random1)

random2 <- lmer(Sp.richness~total.yield+(1|SS)+(1|SSB), data=pred)

summary(random2)

# checking which is better
AIC(random1,random2)


### trying out some plotting ###
library(devtools)
install_github("timnewbold/StatisticalModels")

library(StatisticalModels)

PlotGLMERContinuous(model = random2, data = pred, effects = "total.yield", xlab = "Total crop yield",
                    ylab = "Species Richness", line.cols = c("#79CDCD"), 
outDir = "D:/BIOTA/1_Forest_Cover_Yield/Data Exploration")


### testing for over-dispersion ###

GLMEROverdispersion(model = random2)

# dont think there is overdispersion but testing another model anyway.

random3 <- lmer(logAbun~Use_intensity+(1|SS)+(1|SSB)+(1|SSBS), data=model.data) # doesn't work

