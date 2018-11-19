############################################################
#                                                          #
#              Workflow: Effect of yield and               #
#           proximity to forest on biodiversity            #
#                                                          #
############################################################

# 13th Novemeber 2018, starting to put together workflow for data organisation and
# analysis of the first project in BIOTA.  

rm(list = ls())


# load relevant libraries including Tim's stuff from GitHub

library(devtools)
library(ggplot2)
library(ggthemes)
#install_github("timnewbold/StatisticalModels")
#install_github("timnewbold/predicts-demo",subdir="predictsFunctions")
library(predictsFunctions)
library(StatisticalModels)
library(reshape2)
library(RColorBrewer)
library(sf)
library(raster)
library(cowplot)


# where are the various datasets saved?
datadir <- "D:/BIOTA/Data"


############################################################
#                                                          #
#              Step 1: organise PREDICTS data              #
#                                                          #
############################################################

# read in the complete PREDICTS dataset
# for some reason my download of the data doesn't like to work with the mergeSites
# function, so using the veresion Tim sent round for the PREDICTS workshop.
pred.data <- readRDS(paste0(datadir, "/PREDICTS_2016/database.rds")) # 3250404 rows

### organise using Tim's functions ###

# use correct sampling effort function (this replaces any with NA so the next function works)
pred.sites <- CorrectSamplingEffort(pred.data)

# merge sites: this combines potential subsamples within one site
pred.sites <- MergeSites(pred.sites) # 2906994 rows

# Calculate site level metrics
pred.sites.metrics <- SiteMetrics(pred.sites, extra.cols = c("Predominant_land_use", "SSB", "SSBS")) # 22678 rows


### so now have two datasets, pred.sites maintains species level metrics     ###
### for presence/absence analyses, pred.sites.metrics has site level metrics ###
### (sp richness and abundance)                                              ###


# reduce to cropland sites only for each dataset

# maintain species level dataset as well, cropland sites only
pred.sites.crop <- pred.sites[pred.sites$Predominant_land_use == "Cropland", ] # 274346 rows

# site level data cropland only
pred.sites.metrics.crop <- pred.sites.metrics[pred.sites.metrics$Predominant_land_use == "Cropland", ] # 3179 rows

# remove sites with NA in lat/long columns
pred.sites.crop <- pred.sites.crop[!is.na(pred.sites.crop$Latitude),  ] 
pred.sites.metrics.crop <- pred.sites.metrics.crop[!is.na(pred.sites.metrics.crop$Latitude),  ] 

# predicts datasets summaries

# nstudies
length(unique(pred.sites.crop$SS)) # 152 studies
# nsites
length(unique(pred.sites.crop$SSS)) # 3177 sites
# nspecies
length(unique(pred.sites.crop$Best_guess_binomial)) # 8816 species



############################################################
#                                                          #
#      Step 2: Subset data to those in forest biomes       #
#                                                          #
############################################################

# load in the biome dataset

# read in the biomes shapefile (this has 4 biomes in, just for forest types?)
# humid tropics, dry tropics, temperate, boreal
biomes <- st_read(paste0(datadir, "/biomes/biomes.shp")) %>%
  st_set_crs(4326) # setting the CRS to world lat/longs

# removing the extra grid code column
biomes <- biomes[2:3]

# basic plot of the biomes
plot(biomes)

# first get the biomes for the site level data
pred.sites.metrics.crop$biome <- apply(pred.sites.metrics.crop[ , c("Longitude", "Latitude")], 1, function(x) {  
  # transformation to palnar is required, since sf library assumes planar projection 
  
  pred.sf <- st_transform(st_sfc(st_point(x),crs = 4326), 4326)
  # st_intersects with sparse = FALSE returns a logical matrix
  # with rows corresponds to argument 1 (points) and 
  # columns to argument 2 (polygons)
  
  biomes[which(st_intersects(pred.sf, biomes, sparse = FALSE)), ]$BIOME 
})

# rename biome values from numbers
levels(biomes$BIOME)
pred.sites.metrics.crop$biome <- sub(1, "Boreal", pred.sites.metrics.crop$biome)
pred.sites.metrics.crop$biome <- sub(2, "Dry Tropics", pred.sites.metrics.crop$biome)
pred.sites.metrics.crop$biome <- sub(3, "Humid Tropics", pred.sites.metrics.crop$biome)
pred.sites.metrics.crop$biome <- sub(4, "Temperate", pred.sites.metrics.crop$biome)


ggplot()+
  borders("world", colour = "black", size = 0.3) +
  theme_map() +
  geom_sf(data = biomes, aes(fill = BIOME), alpha = 0.5) +
  geom_point(data = pred.sites.metrics.crop, aes(x = Longitude, y = Latitude, col = biome)) +
  theme(panel.grid.major = element_line(colour = "transparent"))


# remove sites that are not within a forest biome
nrow(pred.sites.metrics.crop[pred.sites.metrics.crop$biome == "integer(0)",])

pred.sites.metrics.crop <- pred.sites.metrics.crop[!pred.sites.metrics.crop$biome == "integer(0)",]




# then get the biomes for the species level data
pred.sites.crop$biome <- apply(pred.sites.crop[ , c("Longitude", "Latitude")], 1, function(x) {  
  # transformation to palnar is required, since sf library assumes planar projection 
  
  pred.sf <- st_transform(st_sfc(st_point(x),crs = 4326), 4326)
  # st_intersects with sparse = FALSE returns a logical matrix
  # with rows corresponds to argument 1 (points) and 
  # columns to argument 2 (polygons)
  
  biomes[which(st_intersects(pred.sf, biomes, sparse = FALSE)), ]$BIOME 
})

# rename biome values from numbers
levels(biomes$BIOME)
pred.sites.crop$biome <- sub(1, "Boreal", pred.sites.crop$biome)
pred.sites.crop$biome <- sub(2, "Dry Tropics", pred.sites.crop$biome)
pred.sites.crop$biome <- sub(3, "Humid Tropics", pred.sites.crop$biome)
pred.sites.crop$biome <- sub(4, "Temperate", pred.sites.crop$biome)

# remove sites that are not within a forest biome
nrow(pred.sites.crop[pred.sites.crop$biome == "integer(0)",])

pred.sites.crop <- pred.sites.crop[!pred.sites.crop$biome == "integer(0)",]

# dataset summaries

# nstudies
length(unique(pred.sites.crop$SS)) # 141 studies
# nsites
length(unique(pred.sites.crop$SSS)) # 2676 sites
# nspecies
length(unique(pred.sites.crop$Best_guess_binomial)) # 8585 species

# save the two datasets
write.csv(pred.sites.crop, file = "D:/BIOTA/1_Forest_Cover_Yield/1. Data/PREDICTS_cropland_forestBiome_sp_level.csv")
write.csv(pred.sites.metrics.crop, file = "D:/BIOTA/1_Forest_Cover_Yield/1. Data/PREDICTS_cropland_forestBiome_site_level.csv")

############################################################
#                                                          #
#     Step 3: Calculate total crop yield for each site     #
#                                                          #
############################################################

# Use the EarthStat global crop yield data to estimate total yield at each site


# where is the yield data saved
yielddir <- "D:/BIOTA/1_Forest_Cover_Yield/1. Data/Yield"

# list files in datadir
files <- list.files(yielddir, pattern = ".csv")

# read in each yield data file (17 crops) (standardise if necesssary - currently commented out)
# then add to a new table including crop name

# space to save standardised data
all_cropdata <- NULL


# file <- files[1]

for(file in files){
  
  # read in the csv file
  cropdata <- read.csv(paste0(yielddir, "/", file))
  
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


# save csv file of the total yield data (and/or standardised yields)
#write.csv(all_cropdata2, paste0("D:/BIOTA/1_Forest_Cover_Yield/Data Exploration/All_Crop_Yield_standardised.csv"))
write.csv(all_cropdata2, paste0("D:/BIOTA/1_Forest_Cover_Yield/Data Exploration/All_Crop_Yield.csv"), row.names = F)

# plot total yield on a global map

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


ggplot(data = all_cropdata2, aes(x = lon, y = lat))+
  borders("world", colour = "gray40", size = 0.3) +
  theme_map() +
  geom_raster(data = all_cropdata2, aes(x = lon, y = lat, fill=totals), inherit.aes = FALSE) +
  geom_point(data = pred.sites.metrics.crop, aes( x = Longitude, y = Latitude), col = "green") +
  #scale_fill_manual(values=cols) +
  #  guides(fill=guide_legend(title="Yield \n(standardised)")) +
  guides(fill=guide_legend(title="Yield \n(Tons per hectare)")) +
  labs(title = paste("Total yield (17 crops), Earthstat")) +
  theme(legend.position="bottom")+
  theme_void()+ coord_equal()

ggsave(file = paste0("D:/BIOTA/1_Forest_Cover_Yield/Data Exploration/Overlay_yield_PREDICTS.png"))


### figure out the yield values at the PREDICTS sites ###
# not a clue if any of this is correct or the best way to do things!

# first, convert both to spatial data objects with the same CRS
# convert the lat/long/value data into a raster
yield.sp <- rasterFromXYZ(all_cropdata2[, c(1,2,20)])  # add a buffer to reduce NAs?

# extract yield values from the predicts points for the site level and species level datasets
pred.sites.crop$yield.vals <- extract(yield.sp, pred.sites.crop[, c("Longitude","Latitude")])
pred.sites.metrics.crop$yield.vals <- extract(yield.sp, pred.sites.metrics.crop[, c("Longitude","Latitude")])



# normality checks on yield totals?
# histogram of totals
p1 <- ggplot(pred.sites.crop, aes(yield.vals)) +
  geom_histogram()

# q-q plot of totals
p2 <- ggplot(pred.sites.crop, aes(sample = log(yield.vals))) +
  stat_qq() + 
  stat_qq_line()

# combine into one
plot_grid(p1, p2) # organise with cowplot function

# save the plot
#ggsave(filename = paste0("D:/BIOTA/1_Forest_Cover_Yield/Data Exploration/Normality_check_YieldTotal_Standardised.png"))
ggsave(filename = paste0("D:/BIOTA/1_Forest_Cover_Yield/Data Exploration/Normality_check_YieldTotal.png"))


# normality checks on yield totals?
# histogram of totals
p1 <- ggplot(pred.sites.metrics.crop, aes(yield.vals)) +
  geom_histogram()

# q-q plot of totals
p2 <- ggplot(pred.sites.metrics.crop, aes(sample = log(yield.vals))) +
  stat_qq() + 
  stat_qq_line()

# combine into one
plot_grid(p1, p2) # organise with cowplot function

# save the plot
#ggsave(filename = paste0("D:/BIOTA/1_Forest_Cover_Yield/Data Exploration/Normality_check_YieldTotal_Standardised.png"))
ggsave(filename = paste0("D:/BIOTA/1_Forest_Cover_Yield/Data Exploration/Normality_check_YieldTotal.png"))



# some rows have NA in the column for total yield.  May need to remove these.


############################################################
#                                                          #
#       step 4: Calculate distance to forest habitat       #
#                                                          #
############################################################


# Monica is working on code to calculate this from the Global Land Cover dataset on forest cover.



############################################################
#                                                          #
#         Step 5: Additional covariate information         #
#                                                          #
############################################################

# include information on human population density, distance to road??
# these were included in the land use paper.
# any other things to include?



############################################################
#                                                          #
#                    Step 6: Run models                    #
#                                                          #
############################################################


# dif model options (species richness, abundance, presence/absence)

# use functions to run model selection process

# which elements need to be included as fixed and random effects?

# checks for overdispersion and spatial autocorrelation?


#### 1. Species richness models ####


r1 <- GLMER(modelData = pred.sites.metrics.crop,responseVar = "Species_richness",
            fitFamily = "poisson",fixedStruct = "yield.vals",
            randomStruct = "(1|SS)",REML = FALSE)

summary(r1$model)

# Try adding spatial block to account for spatial differences among sites within studies
r2 <- GLMER(modelData = pred.sites.metrics.crop,responseVar = "Species_richness",
            fitFamily = "poisson",fixedStruct = "yield.vals",
            randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

summary(r2$model)


# checking which is better
AIC(r1$model,r2$model)

# Species richness models are often overdispersed - is this the case?
GLMEROverdispersion(model = r2$model)

# We can control for this by fitting an observation-level random effect
r3 <- GLMER(modelData = pred.sites.metrics.crop, responseVar = "Species_richness",
            fitFamily = "poisson",fixedStruct = "yield.vals",
            randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# Is this better than the previous random-effects structure?
AIC(r2$model,r3$model)

# Has it removed the overdispersion?
GLMEROverdispersion(model = r3$model)

# Now we will run backward stepwise selection to see if land use has a significant effect
sr1 <- GLMERSelect(modelData = pred.sites.metrics.crop, responseVar = "Species_richness",
                   fitFamily = "poisson",fixedFactors = c("yield.vals", "biome"),
                   randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)")



#### 2. Abundance model ####


r1 <- GLMER(modelData = pred.sites.metrics.crop,responseVar = "Total_abundance",
            fitFamily = "gaussian",fixedStruct = "yield.vals",
            randomStruct = "(1|SS)",REML = FALSE)

summary(r1$model)

# Try adding spatial block to account for spatial differences among sites within studies
r2 <- GLMER(modelData = pred.sites.metrics.crop,responseVar = "Species_richness",
            fitFamily = "poisson",fixedStruct = "yield.vals",
            randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

summary(r2$model)


# checking which is better
AIC(r1$model,r2$model)

# Species richness models are often overdispersed - is this the case?
GLMEROverdispersion(model = r2$model)

# We can control for this by fitting an observation-level random effect
r3 <- GLMER(modelData = pred.sites.metrics.crop, responseVar = "Species_richness",
            fitFamily = "poisson",fixedStruct = "yield.vals",
            randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# Is this better than the previous random-effects structure?
AIC(r2$model,r3$model)

# Has it removed the overdispersion?
GLMEROverdispersion(model = r3$model)


abundModelSelect <- GLMERSelect(modelData = PREDICTSSites,responseVar = "LogAbund",
                                fitFamily = "gaussian",fixedFactors = "LandUse",
                                fixedTerms = list(logHPD.rs=2,logDistRd.rs=2),
                                randomStruct = "(1|SS)+(1|SSB)",verbose = TRUE)



############################################################
#                                                          #
#                 Step 7: Plots of results                 #
#                                                          #
############################################################

PlotGLMERContinuous(model = r3$model, data = r3$data, effects = "yield.vals", xlab = "Total crop yield",
                    ylab = "Species Richness", line.cols = c("#79CDCD"), 
                    outDir = "D:/BIOTA/1_Forest_Cover_Yield/Data Exploration")



############################################################
#                                                          #
#          Step 8: Checks for overdispersion etc           #
#                                                          #
############################################################


GLMEROverdispersion(model = random2)


