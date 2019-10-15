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
#install_github("timnewbold/predicts-demo",subdir="predictsFunctions", force = TRUE)
library(predictsFunctions)
library(StatisticalModels)
library(reshape2)
library(RColorBrewer)
library(sf)
library(raster)
library(cowplot)
library(rasterVis)
library(gridExtra)
library(mapview)
library(tidyverse)
library(BioStatR)


# where are the various datasets saved?
datadir <- "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/Data"


############################################################
#                                                          #
#              Step 1: organise PREDICTS data              #
#                                                          #
############################################################


### fill in info from the current rmarkdown looking at the PREDICTS data ###


# read in the complete PREDICTS dataset
# for some reason my download of the data doesn't like to work with the mergeSites
# function, so using the veresion Tim sent round for the PREDICTS workshop.
pred.data <- readRDS(paste0(datadir, "/PREDICTS_2016/database.rds")) # 3250404 rows

### organise using Tim's functions ###

# use correct sampling effort function (this replaces any with NA so the next function works)
predicts <- CorrectSamplingEffort(pred.data)

# merge sites: this combines potential subsamples within one site
predicts <- MergeSites(predicts) # 2906994 rows

# Calculate site level metrics
pred.sites.metrics <- SiteMetrics(predicts, extra.cols = c("Predominant_land_use", "SSB", "SSBS")) # 22678 rows


### so now have two datasets, pred.sites maintains species level metrics     ###
### for presence/absence analyses, pred.sites.metrics has site level metrics ###
### (sp richness and abundance)                                              ###


# subset to natural landscapes plus agriculture
pred.sub <- predicts[!predicts$Predominant_land_use %in% c("Urban", "Pasture", "Cannot decide", "Plantation forest"), ]

### remove any rows with NAs in the lat long?
pred.sub <- pred.sub[!is.na(pred.sub$Longitude), ]
# 2117928 rows


# site level data cropland only
sites.sub <- pred.sites.metrics[!pred.sites.metrics$Predominant_land_use %in% c("Urban", "Pasture", "Cannot decide", "Plantation forest"), ]

# remove sites with NA in lat/long columns
sites.sub <- sites.sub[!is.na(sites.sub$Longitude),  ] 
# 15612 rows

# predicts datasets summaries

# Nrecords
nrecs <- nrow(pred.sub)
print(nrecs)

# Nspecies
nsp <- length(unique(pred.sub$Best_guess_binomial))
print(nsp)

# Nstudies
nstudies <- length(unique(pred.sub$SS))
print(nstudies)

# Nsites
nsites <- length(unique(pred.sub$SSS))
print(nsites)

# Nsources
nsource <- length(unique(pred.sub$Source_ID))
print(nsource)

# saving just the site level data for now
write.csv(sites.sub, "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/Data/PREDICTS_2016/PREDICTS_sitelevel_agri_natural_11042019.csv")

sites.sub <- read.csv("C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/Data/PREDICTS_2016/PREDICTS_sitelevel_agri_natural_11042019.csv")

############################################################
#                                                          #
#      Step 2: Subset data to those in forest biomes       #
#                                                          #
############################################################

# This section has replaced (10/04/2019) the original biomes section.  This can still be found in the Biomes_forest_cover.r script.

# where is the dataset?
biodir <- paste0(datadir, "/WWF_Ecoregions_Biomes")

# load the WWF ecoregions shapefile
ecobio <- st_read(dsn = paste0(biodir, "/", "wwf_terr_ecos.shp"))

# there are 14 biomes, then 2 for water and ice (?)
sort(unique(ecobio['BIOME']$BIOME)) 

# what does the dataset look like.  A polygon for each region.
ecobio[ecobio$BIOME %in% c(1:14), 'BIOME']


### will need to add into workflow, just use a subset of predicts data for now to test method. ###

# convert the biome info to a factor so its not plotted as a continuous variable
ecobio$BIOME <- as.factor(ecobio$BIOME)

# plot
ggplot()+
  borders("world", colour = "black", size = 0.3) +
  theme_map() +
  geom_sf(data = ecobio[, 'BIOME'], aes(fill = BIOME), size = 0.2) +
  theme(panel.grid.major = element_line(colour = "transparent"))

# save plot
ggsave(file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/2. Figures/WWF_ALL_biomes.tiff")



# plot
ggplot()+
  borders("world", colour = "black", size = 0.3) +
  theme_map() +
  geom_sf(data = ecobio[ecobio$BIOME %in% c(1:6, 12), 'BIOME'], aes(fill = BIOME), size = 0.2) +
  theme(panel.grid.major = element_line(colour = "transparent"))

# save plot
ggsave(file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/2. Figures/WWF_Forest_biomes.tiff")



# get the lat longs only and convert to points
sites.xy <- sites.sub[, c("Longitude", "Latitude")] # 15612


# take only the biome attribute
wwfbiomes <- ecobio['BIOME']

# use this fancy mapping business in mapview package to have a look at the plotted points and biomes 

# transform the lat/longs into spatial points, testing just a few here,
# set the CRS to the same as that of the polygons
pointsT <- st_as_sf(sites.xy, coords = c("Longitude", "Latitude"))
st_crs(pointsT) <- "+init=epsg:4326 +proj=longlat +ellps=WGS84"
pointsT <- st_transform(pointsT, crs = st_crs(wwfbiomes))
print(st_crs(pointsT))

# have a check on the map again
mapview::mapview(ecobio['BIOME']) + mapview::mapview(pointsT, color = "white", col.regions = "black")


# get the point/polygon intersections

# The polygons seem to be slightly mismatched with the land area, for the NAs try the within distance function to fill in the gaps.
biomes <- sf::st_intersects(pointsT, wwfbiomes)

# replace any that are NULL with NA
biomes[sapply(biomes, function(x) length(x)==0L)] <- NA

# Get the values for the biomes
biomes <- as.character(wwfbiomes$BIOME[unlist(biomes)])

# add info to the data table
sites.sub$Forest_biome <- biomes

# check for the NAs
View(sites.sub[is.na(sites.sub$Forest_biome), ])
as.numeric(rownames(sites.sub[is.na(sites.sub$Forest_biome), ]))

# all points
mapview::mapview(ecobio['BIOME']) + mapview::mapview(pointsT, color = "white", col.regions = "black")

# NAs only
mapview::mapview(wwfbiomes) + mapview::mapview(pointsT[as.numeric(rownames(sites.sub[is.na(sites.sub$Forest_biome), ]))
                                                       , ], color = "white", col.regions = "black")

#### now try and fill in the NAs ####


# subset the points to just those with NA
pointsST_NA <- st_as_sf(sites.sub[is.na(sites.sub$Forest_biome), ], coords = c("Longitude", "Latitude"))
st_crs(pointsST_NA) <- "+init=epsg:4326 +proj=longlat +ellps=WGS84"
pointsST_NA <- st_transform(pointsST_NA, crs = st_crs(wwfbiomes))
print(st_crs(pointsST_NA))

# try and get the biomes for the NAs, using the within distance function to get the biome within 5km (? i think)
#biomes_NA <- sf::st_is_within_distance(pointsST_NA, wwfbiomes, dist = 5000, sparse = FALSE)

biomes_NA <- sf::st_join(pointsST_NA, wwfbiomes, st_is_within_distance, dist = 10000)

pointsST_NA$Forest_biome <- NA


vals <- unique(biomes_NA$X)

for(i in 1:length(vals)){
  
  x <- biomes_NA[biomes_NA$X == vals[i], ]
  nrow(x)
  
  # if the multiple values are the same, just keep one
  if(length(unique(x$BIOME)) == 1){
  
    pointsST_NA$Forest_biome[i] <- unique(x$BIOME)
  }
  else{pointsST_NA$Forest_biome[i] <- NA}
}

# check the NA ones manually on the map, what should they be.
View(pointsST_NA[is.na(pointsST_NA$Forest_biome), ])

pointsST_NA[pointsST_NA$X == 4559, "Forest_biome"] <- 14
pointsST_NA[pointsST_NA$X == 4560, "Forest_biome"] <- 14
pointsST_NA[pointsST_NA$X == 4567, "Forest_biome"] <- 14
pointsST_NA[pointsST_NA$X == 4568, "Forest_biome"] <- 14
pointsST_NA[pointsST_NA$X == 12091, "Forest_biome"] <- NA
pointsST_NA[pointsST_NA$X == 13983, "Forest_biome"] <- 1

# check the NA ones manually on the map, what should they be.
View(pointsST_NA[is.na(pointsST_NA$Forest_biome), ]) # now there is only 1 NA


# add these in to the sites.sub table
sites.sub[is.na(sites.sub$Forest_biome), "Forest_biome"] <- pointsST_NA$Forest_biome
 
nrow(sites.sub[is.na(sites.sub$Forest_biome), ]) # only 1, where biome of island was unknown



# subset to forest biomes only
# forest biomes include 1:6 and 12
table(sites.sub$Forest_biome)

sites.sub <- sites.sub[sites.sub$Forest_biome %in% c(1:6, 12), ] # 12176

table(sites.sub$Forest_biome)

nrow(sites.sub) # 12176

# rename the forest biomes
final.data.trans$Forest_biome <- sub(12, "Mediterranean Forests, Woodlands & Scrub", final.data.trans$Forest_biome)
final.data.trans$Forest_biome <- sub(1, "Tropical & Subtropical Moist Broadleaf Forests", final.data.trans$Forest_biome)
final.data.trans$Forest_biome <- sub(2, "Tropical & Subtropical Dry Broadleaf Forests", final.data.trans$Forest_biome)
final.data.trans$Forest_biome <- sub(3, "Tropical & Subtropical Coniferous Forests", final.data.trans$Forest_biome)
final.data.trans$Forest_biome <- sub(4, "Temperate Broadleaf & Mixed Forests", final.data.trans$Forest_biome)
final.data.trans$Forest_biome <- sub(5, "Temperate Conifer Forests", final.data.trans$Forest_biome)
final.data.trans$Forest_biome <- sub(6, "Boreal Forests/Taiga", final.data.trans$Forest_biome)

final.data.trans$Forest_biome <- factor(final.data.trans$Forest_biome,
                                        levels=c("Temperate Broadleaf & Mixed Forests", "Temperate Conifer Forests",
                                                 "Boreal Forests/Taiga", "Mediterranean Forests, Woodlands & Scrub",                     "Tropical & Subtropical Coniferous Forests", "Tropical & Subtropical Dry Broadleaf Forests",
                                                 "Tropical & Subtropical Moist Broadleaf Forests"))
# save the two datasets
#write.csv(pred.sub, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_sp_level.csv", row.names = F)
write.csv(sites.sub, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_site_level.csv", row.names = F)

############################################################
#                                                          #
#     Step 3: Input the produtcion and fertiliser          #
#                 intensity information                    #
#                                                          #
############################################################


# total production per grid cell and production pe rhectare per grid cell have
# been estimated globally from the EarthStat data.

# read in the relevent rasters and extract values to the predicts data tables

# the creation of the total rasters and per hectare rasters is in:
# 1. EarthStat data organisation_raster method.r script


#pred.sub <- read.csv("C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_sp_level.csv")
sites.sub <- read.csv("C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_site_level.csv")
#sites.sub <- sites.sub[, 2:ncol(sites.sub)]


total.prod <- raster("C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/Earthstat_prod_total_newmethod.tif")
prodperhec <- raster("C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/Earthstat_prod_perhec_newmethod.tif")
total.fert <- raster("C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/Earthstat_fert_total_newmethod.tif")
fertperhec <- raster("C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/Earthstat_fert_perhec_newmethod_17crops.tif")


# plots already generated in the data organisation script


# convert the PREDICTS lat/longs into spatial points
#pred.sub_xy <- pred.sub[, c("Longitude", "Latitude")]
#pred.sub_xy <- SpatialPoints(pred.sub_xy)

sites.sub_xy <- sites.sub[, c("Longitude", "Latitude")]
sites.sub_xy <- SpatialPoints(sites.sub_xy)


# get the intersection values, what are the production values for each site?
#pred.sub$prod.total <- extract(total.prod, pred.sub_xy, na.rm = FALSE)
sites.sub$prod.total <- extract(total.prod, sites.sub_xy, na.rm = FALSE)

#pred.sub$prod.per.hec <- extract(prodperhec, pred.sub_xy, na.rm = FALSE)
sites.sub$prod.per.hec <- extract(prodperhec, sites.sub_xy, na.rm = FALSE)

#pred.sub$fert.total <- extract(total.fert, pred.sub_xy, na.rm = FALSE)
sites.sub$fert.total <- extract(total.fert, sites.sub_xy, na.rm = FALSE)

#pred.sub$fert.per.hec <- extract(fertperhec, pred.sub_xy, na.rm = FALSE)
sites.sub$fert.per.hec <- extract(fertperhec, sites.sub_xy, na.rm = FALSE)

# where there are NAs, change to 0
sites.sub[is.na(sites.sub$prod.total), "prod.total"] <- 0
sites.sub[is.na(sites.sub$prod.per.hec), "prod.per.hec"] <- 0
sites.sub[is.na(sites.sub$fert.total), "fert.total"] <- 0
sites.sub[is.na(sites.sub$fert.per.hec), "fert.per.hec"] <- 0



# save this
#write.csv(pred.sub, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_sp_level.csv", row.names = FALSE)
write.csv(sites.sub, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_site_level.csv", row.names = FALSE)


############################################################
#                                                          #
#       step 4: Determine the number of crops              #
#                   in each grid square                    #  
#                                                          #
############################################################

# add in the number of crops for each grid square using the raster generated from the earthstat data

# read in dataset so far if not in the environment
#pred.sub <- read.csv("C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_sp_level.csv")
sites.sub <- read.csv("C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_site_level.csv")


# read in the ncrop raster that was made using the earthstat data
ncrops <- raster("C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/Earthstat_NCrops_newmethod.tif")


# convert the PREDICTS lat/longs into spatial points
#pred.sub_xy <- pred.sub[, c("Longitude", "Latitude")]
#pred.sub_xy <- SpatialPoints(pred.sub_xy)

#sites.sub_xy <- sites.sub[, c("Longitude", "Latitude")]
#sites.sub_xy <- SpatialPoints(sites.sub_xy)


# get the intersection values, what are the production values for each site?
#pred.sub$ncrop <- extract(ncrops, pred.sub_xy, na.rm = FALSE)
sites.sub$ncrop <- extract(ncrops, sites.sub_xy, na.rm = FALSE)

# replace NAs with 0
sites.sub[is.na(sites.sub$ncrop), "ncrop"] <- 0

#write.csv(pred.sub, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_sp_level.csv", row.names = FALSE)
write.csv(sites.sub, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_site_level.csv", row.names = FALSE)




############################################################
#                                                          #
#       step 5: Organise EarthStat data on fraction        #
#                   of area harvested                      #
#                                                          #
############################################################

# read in the data files compiled so far
#pred.sub <- read.csv(file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_sp_level.csv")
#sites.sub <- read.csv(file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_site_level.csv")


# use the raster layer creasted in the 1. Earthstat data organisation script

fract.totals <- raster("C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/Earthstat_areaharv_total_frac_newmethod.tif")



# convert the PREDICTS lat/longs into spatial points
#pred.sub_xy <- pred.sub[, c("Longitude", "Latitude")]
#pred.sub_xy <- SpatialPoints(pred.sub_xy)

#sites.sub_xy <- sites.sub[, c("Longitude", "Latitude")]
#sites.sub_xy <- SpatialPoints(sites.sub_xy)


# get the intersection values, what are the production values for each site?
#pred.sub$frac.harvested <- extract(fract.totals, pred.sub_xy, na.rm = FALSE)
sites.sub$frac.harvested <- extract(fract.totals, sites.sub_xy, na.rm = FALSE)

# replace NAs with 0
sites.sub[is.na(sites.sub$frac.harvested), "frac.harvested"] <- 0

# save
#write.csv(pred.sub, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_fracharv_sp_level.csv", row.names = FALSE)
write.csv(sites.sub, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_frac_harv_site_level.csv", row.names = FALSE)





############################################################
#                                                          #
####   step 6: Calculate distance to forest habitat     ####
#                                                          #
############################################################


# Monica has estimated the values for distance to dense forest using the Kehoe
# land cover dataset.

# second attempt by Monica of calculating distance to forest,  re-incorporate this when running final workflow.

sites.sub <- read.csv("C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_frac_harv_site_level.csv")


# read in the latest sites.sub data

datadir <- "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data"

pred<-read.csv(paste0(datadir, "/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_frac_harv_site_level.csv"))

# extract just the lat/lons
pred<-data.frame(pred$Longitude, pred$Latitude)

# change the column names
names(pred)=c("lon","lat")

# load the kehoe dataset
kehoe<-raster("C:/Users/charl/Dropbox/BIOTA/Land use datasets/Kehoe/kehoe_df_cs.tif")

plot(kehoe)
kehoe

# create a vector to store results
min_dist<-vector()

# loop through each point in the lat/lons table
for (i in 1:nrow(pred)){
  # Here Monica has created a grid around the Predicts point.
  xmin<-pred[i,]$lon -0.5
  xmax<-pred[i,]$lon +0.5
  ymin<-pred[i,]$lat -0.5
  ymax<-pred[i,]$lat +0.5
  test<-crop(kehoe, c(xmin,xmax,ymin,ymax))
  #plot(test)
  
  #single point
  pred_test=pred[i,]
  
  #subset to the cropped area
  test_pts <- rasterToPoints(test, function(x){!is.na(x)})
  
  #testing if grid does not have any forest at all, go a bit bigger
  if (dim(test_pts)[1]<=1){
    xmin<-pred[i,]$lon -1
    xmax<-pred[i,]$lon +1
    ymin<-pred[i,]$lat -1
    ymax<-pred[i,]$lat +1
    test<-crop(kehoe, c(xmin,xmax,ymin,ymax))
    #single point
    pred_test=pred[i,]
    #subset to the cropped area
    test_pts <- rasterToPoints(test, function(x){!is.na(x)})
  }
  if (dim(test_pts)[1]<=1){
    xmin<-pred[i,]$lon -10
    xmax<-pred[i,]$lon +10
    ymin<-pred[i,]$lat -10
    ymax<-pred[i,]$lat +10
    test<-crop(kehoe, c(xmin,xmax,ymin,ymax))
    #single point
    pred_test=pred[i,]
    #subset to the cropped area
    test_pts <- rasterToPoints(test, function(x){!is.na(x)})
  }
  if (dim(test_pts)[1]<=1){
    xmin<-pred[i,]$lon -15
    xmax<-pred[i,]$lon +15
    ymin<-pred[i,]$lat -15
    ymax<-pred[i,]$lat +15
    test<-crop(kehoe, c(xmin,xmax,ymin,ymax))
    #single point
    pred_test=pred[i,]
    #subset to the cropped area
    test_pts <- rasterToPoints(test, function(x){!is.na(x)})
  }
  if (dim(test_pts)[1]<=1){
    xmin<-pred[i,]$lon -50
    xmax<-pred[i,]$lon +50
    ymin<-pred[i,]$lat -10
    ymax<-pred[i,]$lat +10
    test<-crop(kehoe, c(xmin,xmax,ymin,ymax))
    #single point
    pred_test=pred[i,]
    #subset to the cropped area
    test_pts <- rasterToPoints(test, function(x){!is.na(x)})
  }
  #actual minimum distance calc - takes the minimum of the distances in the grid with non-na forest values.
  min_dist[i]<-min(pointDistance(pred_test, test_pts[,1:2], lonlat = TRUE, allpairs = TRUE)/1000) #/1000 so its in km rather than m?
  test_pts<-NULL # clear the test point
}

sites.sub <- data.frame(sites.sub, min_dist)


#write.csv(pred.sub, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_fracharv_dist_sp_level.csv", row.names = FALSE)
write.csv(sites.sub, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_frac_harv_dist_site_level.csv", row.names = FALSE)



############################################################
#                                                          #
####    step 7: Calculate number of land covers in      ####
#             in buffer around predicts site               #
#                                                          #
############################################################

# this currently uses the Kehoe land use dataset, land cover classifications have been 
# merged to ignore the livestock and suitability elements.  This was done in script:
# "Landcovers_Kehoe_data_sorting.r"

# This code is copied in from the "Calculating_n_landcovers.r" script.

#pred.sub <- read.csv("C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_fracharv_dist_sp_level.csv")
sites.sub <- read.csv("C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_frac_harv_dist_site_level.csv")



### Starting with the Kehoe landcover dataset ###

# where is the dataset
datadir <- "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/Data/Global-Land-Systems-Map_GeoTIFF"

# read in the kehoe raster
landuse <- raster(paste0(datadir, "/", "GLS_v02.bj_RECLASSIFIED.tif"))

# plot dataset for draft figures
plot(landuse)

# get the PREDICTS points, porbably read in previously if completing script all in one
sites.sub_xy <- sites.sub[, c("Longitude", "Latitude")]

# just a simple plot of the predicts points
map.world <- map_data('world')

ggplot(sites.sub_xy) +
  geom_map(data=map.world, map=map.world,
           aes(x=long, y=lat, group=group, map_id=region),
           fill="white", colour="grey", size=0.5) +
  geom_point(aes(x=Longitude, y=Latitude), size = 0.2) +
  scale_fill_distiller(palette = "Spectral", na.value = "white", name = "Kilograms") +
  ggtitle("PREDICTS sites") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme(text = element_text(size = 10))

ggsave(filename = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/2. Figures/Basic_pred_sites.tiff")


# try a selection of buffers for different scales.

# convert to spatial points
#sites.sub_xy <- sites.sub_xy[1:100, ] # for testing
sites.sub_xy <- SpatialPoints(sites.sub_xy, proj4string = CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84"))

# create the buffer options around the data points
sites.buffer.100 <- buffer(sites.sub_xy, width = 100, dissolve = FALSE)
sites.buffer.500 <- buffer(sites.sub_xy, width = 500, dissolve = FALSE)
sites.buffer.1k <- buffer(sites.sub_xy, width = 1000, dissolve = FALSE)
sites.buffer.3k <- buffer(sites.sub_xy, width = 3000, dissolve = FALSE)
sites.buffer.5k <- buffer(sites.sub_xy, width = 5000, dissolve = FALSE)


# extract the raster values for the buffer sites
landcovers_extract.100 <- extract(landuse, sites.buffer.100, df = TRUE)
landcovers_extract.500 <- extract(landuse, sites.buffer.500, df = TRUE)
landcovers_extract.1k <- extract(landuse, sites.buffer.1k, df = TRUE)
landcovers_extract.3k <- extract(landuse, sites.buffer.3k, df = TRUE)
landcovers_extract.5k <- extract(landuse, sites.buffer.5k, df = TRUE)


# write a function to look at the number of unique landcovers for each ID. 

calc_landcovers <- function(x){
  
  results <- NULL
  
  for(i in unique(x$ID)){
    
    subtab <- x[x$ID == i, ]
    
    nlandcovers <- length(unique(subtab$GLS_v02.bj))
    
    result <- c(i, nlandcovers)
    
    results <- rbind(results, result)
    
  }
  
  colnames(results) <- c("ID", "Nlandcovers")
  
  return(results)
  
}

# create the table for each set with n landcovers
res.100 <- calc_landcovers(landcovers_extract.100)
res.500 <- calc_landcovers(landcovers_extract.500)
res.1k <- calc_landcovers(landcovers_extract.1k)
res.3k <- calc_landcovers(landcovers_extract.3k)
res.5k <- calc_landcovers(landcovers_extract.5k)

# add these values to the sites dataset

sites.sub$landcovers.100 <- res.100[,2]
sites.sub$landcovers.500 <- res.500[,2]
sites.sub$landcovers.1k <- res.1k[,2]
sites.sub$landcovers.3k <- res.3k[,2]
sites.sub$landcovers.5k <- res.5k[,2]

#### need to think about altering values of kehoe dataset, one value per land cover type ####


#write.csv(pred.sub, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_fracharv_dist_sp_level.csv", row.names = FALSE)
write.csv(sites.sub, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_frac_harv_dist_landcovers_site_level.csv", row.names = FALSE)






##%######################################################%##
#                                                          #
####       Step 9: determine the homogeneity (and       ####
#        potentially other metrics) for each site          #
#                                                          #
##%######################################################%##

# This code organises a couple of heterogeneity data options from the Tuanmu 2015 paper:

# 1. Tuanmu, M. N. & Jetz, W. A global, remote sensing-based characterization of terrestrial
# habitat heterogeneity for biodiversity and ecosystem modelling. 
# Glob. Ecol. Biogeogr. 24, 1329-1339 (2015).


# oringal code in script: "Organising_heterogeneoty_data.r"

# load in the most up-to date dataset

sites.sub <- read.csv("C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_frac_harv_dist_landcovers_site_level.csv")


datadir <- "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/Data/Landscape_heterogeneity_Tuanmu2015"


# read in the required datasets
homogen <- raster(paste0(datadir, "/", "Homogeneity_01_05_5km_uint16.tif"))


# info from website:
# values of the data layers should be mulitplied by 0.0001 to obtain the actual values of the metrics.
homogen <- homogen*0.0001

# extract the PREDICTS points
sites.sub_xy <- sites.sub[, c("Longitude", "Latitude")]

sites.sub_xy <- SpatialPoints(sites.sub_xy, proj4string = CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84"))

# extract the dataset info for the PREDICTS sites
sites.sub$homogen <- extract(homogen, sites.sub_xy, na.rm = TRUE, buffer = 5000, fun = mean)

#Some NAs, have a look at these on a map
mapview::mapview(homogen) + mapview::mapview(sites.sub_xy[as.numeric(rownames(sites.sub[is.na(sites.sub$homogen), ]))
                                                       , ], color = "white", col.regions = "black")
# cannot view complete resolution as it is too big to handle, perhaps the points are slightly out of the available pixels.


# save the dataset
write.csv(sites.sub, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_frac_harv_dist_landcovers_heterog_site_level.csv", row.names = FALSE)


##%######################################################%##
#                                                          #
####                 Assess the dataset                 ####
#                                                          #
##%######################################################%##

# remove any rows that have NA in the variable columns
nrow(sites.sub[is.na(sites.sub$Forest_biome), ]) #0
nrow(sites.sub[is.na(sites.sub$prod.total), ]) #0
nrow(sites.sub[is.na(sites.sub$prod.per.hec), ]) #0
nrow(sites.sub[is.na(sites.sub$fert.total), ]) #0
nrow(sites.sub[is.na(sites.sub$fert.per.hec), ]) #0
nrow(sites.sub[is.na(sites.sub$ncrop), ]) #0
nrow(sites.sub[is.na(sites.sub$frac.harvested), ]) #0
nrow(sites.sub[is.na(sites.sub$min_dist), ]) #0
nrow(sites.sub[is.na(sites.sub$landcovers.100), ]) #0
nrow(sites.sub[is.na(sites.sub$landcovers.500), ]) #0
nrow(sites.sub[is.na(sites.sub$landcovers.1k), ]) #0
nrow(sites.sub[is.na(sites.sub$landcovers.3k), ]) #0
nrow(sites.sub[is.na(sites.sub$landcovers.5k), ]) #0
nrow(sites.sub[is.na(sites.sub$homogen), ]) # 61

sites.sub <- sites.sub[!is.na(sites.sub$homogen), ] # 12115
nrow(sites.sub)

# remove those sites that have "Cannot decide" as a use intensity
nrow(sites.sub[sites.sub$Use_intensity == "Cannot decide", ]) # 1410

sites.sub <- sites.sub[!sites.sub$Use_intensity == "Cannot decide", ] # 10705
nrow(sites.sub)

write.csv(sites.sub, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_frac_harv_dist_landcovers_heterog_site_level_CLEANED.csv", row.names = FALSE)


############################################################
#                                                          #
#               Step 8: Explore variables                  #
#                                                          #
############################################################

# make sure that the dataset with all the covariates is read in

sites.sub <- read.csv("C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_NatPlusCrop_forestBiome_Prod_Fert_ncrop_frac_harv_dist_landcovers_heterog_site_level_CLEANED_plusfieldsize_resid_RANGE_Hansen.csv")
View(sites.sub)

### The continuous variables ###

# Production total
# Production per hectare
# Fertilisation total
# Fertiliser per hectare
# number of crops
# fraction of area harvested
# distance to forest
# number of neearby land covers

### The categorical variables ###

# Use intensity (PREDICTS variable)
# Forest Biome

### The random effects ###

# SS
# SSB
# SSBS

# look at correlations of untransformed data
# This has been added to the supplementary

panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1))
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 4* abs(cor(x, y)))
}

pairs(sites.sub[!is.na(sites.sub$Hansen_mindist), c(22:27, 33:34, 45)], 
      upper.panel=panel.cor, 
      diag.panel = panel.hist, 
      main = "")




# subset to desired columns
final.data <- sites.sub[, c(1,9:11, 14:21, 24, 33, 34, 37, 45)]
#View(final.data)

# now to log some continuous variables and standardise the rest
final.data.trans <- final.data


# log transform some of the continuous variables
final.data.trans$fert.total_log <-log(final.data.trans$fert.total+1)
#final.data.trans$min.dist_log <-log(final.data.trans$min_dist)
final.data.trans$Hansen_mindist_log <-log(final.data.trans$Hansen_mindist+1) # there are 0s


# standardise all continuous variables
final.data.trans$landcovers.5k <- scale(final.data.trans$landcovers.5k)
final.data.trans$homogen <- scale(final.data.trans$homogen)
final.data.trans$fert.total_log <- scale(final.data.trans$fert.total_log)
final.data.trans$Hansen_mindist_log <-scale(final.data.trans$Hansen_mindist_log)


# get data sections for the scaling info for plotting later
Hansen_mindist_log <-final.data.trans$Hansen_mindist_log
landcovers.5k <- final.data.trans$landcovers.5k
fert.total_log <- final.data.trans$fert.total_log
homogen <- final.data.trans$homogen

table(final.data.trans$Use_intensity)

# set factor levels in the best order
final.data.trans$Use_intensity <- factor(final.data.trans$Use_intensity,
                                         levels=c("Minimal use", "Light use", "Intense use"))

table(final.data.trans$Use_intensity)


# set the forest biome to names rather than numbers
final.data.trans$Forest_biome <- sub(12, "Mediterranean Forests, Woodlands & Scrub", final.data.trans$Forest_biome)
final.data.trans$Forest_biome <- sub(1, "Tropical & Subtropical Moist Broadleaf Forests", final.data.trans$Forest_biome)
final.data.trans$Forest_biome <- sub(2, "Tropical & Subtropical Dry Broadleaf Forests", final.data.trans$Forest_biome)
final.data.trans$Forest_biome <- sub(3, "Tropical & Subtropical Coniferous Forests", final.data.trans$Forest_biome)
final.data.trans$Forest_biome <- sub(4, "Temperate Broadleaf & Mixed Forests", final.data.trans$Forest_biome)
final.data.trans$Forest_biome <- sub(5, "Temperate Conifer Forests", final.data.trans$Forest_biome)
final.data.trans$Forest_biome <- sub(6, "Boreal Forests/Taiga", final.data.trans$Forest_biome)



final.data.trans$Forest_biome <- factor(final.data.trans$Forest_biome,
                                        levels=c("Temperate Broadleaf & Mixed Forests", "Temperate Conifer Forests",
                                                 "Boreal Forests/Taiga", "Mediterranean Forests, Woodlands & Scrub",
                                                 "Tropical & Subtropical Coniferous Forests", 
                                                 "Tropical & Subtropical Dry Broadleaf Forests",
                                                 "Tropical & Subtropical Moist Broadleaf Forests"))

# set land use as character variable
final.data.trans$Predominant_land_use <- as.character(final.data.trans$Predominant_land_use)

# remove the secondary indeterminate age
#final.data.trans <- final.data.trans[!final.data.trans$Predominant_land_use == "Secondary vegetation (indeterminate age)", ] # 9490 rows

# 03/10/2019
# decided to keep secondary veg together
final.data.trans$Predominant_land_use <- sub("Mature secondary vegetation", "Secondary vegetation", final.data.trans$Predominant_land_use)
final.data.trans$Predominant_land_use <- sub("Intermediate secondary vegetation", "Secondary vegetation", final.data.trans$Predominant_land_use)
final.data.trans$Predominant_land_use <- sub("Young secondary vegetation", "Secondary vegetation", final.data.trans$Predominant_land_use)
final.data.trans[final.data.trans$Predominant_land_use == "Secondary vegetation (indeterminate age)", 'Predominant_land_use'] <- "Secondary vegetation"
# unique(final.data.trans$Predominant_land_use)

table(final.data.trans$Predominant_land_use)

final.data.trans$Predominant_land_use <- as.factor(final.data.trans$Predominant_land_use)

levels(final.data.trans$Predominant_land_use)

# set factor levels of predominant land use
final.data.trans$Predominant_land_use <- factor(final.data.trans$Predominant_land_use,
                                                levels=c("Primary vegetation","Secondary vegetation", "Cropland"))

table(final.data.trans$Predominant_land_use)

# Primary vegetation  Secondary vegetation  Cropland 
# 4799                   3935                 1971


# remove the rows with NA in Hansen info
final.data.trans <- final.data.trans[!is.na(final.data.trans$Hansen_mindist), ] # 10564



# look at correlations of rescaled data

panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1))
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 4* abs(cor(x, y)))
}

pairs(final.data.trans[, c(20, 21, 33, 34)], 
      upper.panel=panel.cor, 
      diag.panel = panel.hist, 
      main = "")


#pairs(final.data.trans[, c(20,21,33,35)], 
#      upper.panel=panel.cor, 
#      diag.panel = panel.hist, 
#      main = "")




##%######################################################%##
#                                                          #
####          Additional supplementary plots            ####
#                                                          #
##%######################################################%##



boxplot(fert.total ~ Use_intensity, data = final.data.trans, ylab = "Total fertiliser", xlab = "Use intensity", outline = FALSE)
boxplot(log(fert.total+1) ~ Use_intensity, data = final.data.trans, ylab = "log(Total fertiliser +1)", xlab = "Use intensity", outline = FALSE)

boxplot(Hansen_mindist ~ Use_intensity, data = final.data.trans, ylab = "Distance to forest", xlab = "Use intensity", outline = FALSE)
boxplot(log(Hansen_mindist+1) ~ Use_intensity, data = final.data.trans, ylab = "log(Distance to forest +1)", xlab = "Use intensity", outline = FALSE)

write.csv(final.data.trans, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/3. Outputs/Stage_5/SR_data_3PLU.csv", row.names = F)


##%######################################################%##
#                                                          #
####                  Dataset subsets                   ####
#                                                          #
##%######################################################%##

# Subsets of the data required for the abundance models and the range models

## abundance subset
final.data.abun <- final.data.trans[!is.na(final.data.trans$Total_abundance), ] # 9172 rows

# log the abundance values
final.data.abun$logAbun <- log(final.data.abun$Total_abundance+1)

write.csv(final.data.abun, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/3. Outputs/Stage_5/Abun_data_3PLU.csv", row.names = F)


## RCAR 110km 

model.data_RCAR_110km <- final.data.trans[!is.na(final.data.trans$RCAR_110km), ] # 5786 rows

write.csv(model.data_RCAR_110km, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/3. Outputs/Stage_5/RCAR_data_3PLU.csv", row.names = F)




##%######################################################%##
#                                                          #
####                 Get data summaries                 ####
#                                                          #
##%######################################################%##

# get the number of studies, sites, etc for the paper

length(unique(final.data.trans$SS)) # 484
length(unique(final.data.trans$SSB)) # 1169

length(unique(final.data.abun$SS)) # 433
length(unique(final.data.abun$SSB)) # 1112

length(unique(model.data_RCAR_110km$SS)) # 316
length(unique(model.data_RCAR_110km$SSB)) # 800






############################################################
#                                                          #
####               Step 9: Run models                   ####
#                                                          #
############################################################


# dif model options (species richness, abundance, presence/absence)

#### 1. Species richness models ####

model.data <- final.data.trans[, c(1,2,6:9,11:12,14:15, 18:19)]

system.time({sr1.2 <- GLMERSelect(modelData = model.data, 
                      responseVar = "Species_richness",
                      fitFamily = "poisson", 
                      fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity"),
                      fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1),
                      randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)", 
                      fixedInteractions = c(#"Predominant_land_use:Use_intensity", 
                                            "Predominant_land_use:fert.total_log", 
                                            "Predominant_land_use:Hansen_mindist_log",
                                            "Predominant_land_use:landcovers.5k", 
                                            "Predominant_land_use:homogen", 
                                            "Use_intensity:fert.total_log", 
                                            "Use_intensity:Hansen_mindist_log",
                                            "Use_intensity:landcovers.5k",
                                            "Use_intensity:homogen"), verbose = TRUE)})

#fixed-effect model matrix is rank deficient so dropping 4 columns / coefficients

save(sr1.2, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/3. Outputs/Stage_5/Stage5_final_models_SPECIESRICHNESS_03102019_noInt.rdata")

load(file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/3. Outputs/Stage_5/Stage5_final_models_SPECIESRICHNESS_03102019_noInt.rdata")

#load(file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/3. Outputs/Stage_5/Stage5_final_models_SPECIESRICHNESS_03102019_Int.rdata")


sr1.2$model

# rerun the final model without the polynomials, better for using the coefs in projections later

# Species_richness ~ Predominant_land_use + Forest_biome + Use_intensity +  
# poly(Hansen_mindist_log, 1) + fert.total_log +
# Predominant_land_use:fert.total_log + Use_intensity:fert.total_log  + 
# (1 | SS) + (1 | SSB) + (1 | SSBS)

# rerun the selected model, make sure using REML
# this function doesn't like unneccesary spaces
srmod <- GLMER(modelData = final.data.trans, responseVar = "Species_richness", fitFamily = "poisson",
      fixedStruct = "Predominant_land_use + Forest_biome + Use_intensity + Hansen_mindist_log + fert.total_log + Predominant_land_use:fert.total_log + Use_intensity:fert.total_log",
      randomStruct = "(1|SS) + (1|SSB) + (1|SSBS)", REML = TRUE)


srmod$model
summary(srmod$model)

# save the coefficents of the model
coefs <- fixef(srmod$model)

# save the coefficients
write.csv(coefs, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/3. Outputs/Stage_5/Stage5_final_models_SPECIESRICHNESS_nopoly_3PLU_14102019_coefs.csv")




#### 2. Abundance model ####


system.time({sr2.2 <- GLMERSelect(modelData = final.data.abun, 
                      responseVar = "logAbun",
                      fitFamily = "gaussian", 
                      fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity"),
                      fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1),
                      randomStruct = "(1|SS)+(1|SSB)", 
                      fixedInteractions = c(#"Predominant_land_use:Use_intensity", 
                                            "Predominant_land_use:fert.total_log", 
                                            "Predominant_land_use:Hansen_mindist_log", 
                                            "Predominant_land_use:landcovers.5k", 
                                            "Predominant_land_use:homogen", 
                                            "Use_intensity:fert.total_log", 
                                            "Use_intensity:Hansen_mindist_log", 
                                            "Use_intensity:landcovers.5k", 
                                            "Use_intensity:homogen"), verbose = TRUE)})

#fixed-effect model matrix is rank deficient so dropping 4 columns / coefficients

sr2.2$model
save(sr2.2, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/3. Outputs/Stage_5/Stage5_final_models_ABUNDANCE_nopoly_3PLU_14102019.rdata")



#system.time({sr2.2 <- GLMERSelect(modelData = final.data.abun, 
#                                responseVar = "logAbun",
#                                fitFamily = "gaussian", 
#                               fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity"),
#                                fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1),
#                                randomStruct = "(1|SS)+(1|SSB)", 
#                                fixedInteractions = c("Predominant_land_use:Use_intensity",
#                                                      "Predominant_land_use:poly(fert.total_log,1)", 
#                                                      "Predominant_land_use:poly(Hansen_mindist_log,1)", 
#                                                      "Predominant_land_use:poly(landcovers.5k,1)", 
#                                                      "Predominant_land_use:poly(homogen,1)", 
#                                                      "Use_intensity:poly(fert.total_log,1)", 
#                                                      "Use_intensity:poly(Hansen_mindist_log,1)", 
#                                                      "Use_intensity:poly(landcovers.5k,1)", 
#                                                      "Use_intensity:poly(homogen,1)"), verbose = TRUE)})



# logAbun ~ Predominant_land_use + Forest_biome + Use_intensity +  
# poly(fert.total_log, 1) + poly(Hansen_mindist_log, 1) + poly(landcovers.5k, 1) + poly(homogen, 1) + 
# Predominant_land_use:fert.total_log + Predominant_land_use:homogen + 
# Use_intensity:fert.total_log +  Use_intensity:homogen + 
# fert.total_log + homogen + (1 | SS) + (1 | SSB) # why are these in there twice?


abmod <- GLMER(modelData = final.data.abun, responseVar = "logAbun", fitFamily = "gaussian",
               fixedStruct = "Predominant_land_use + Forest_biome + Use_intensity + Hansen_mindist_log + fert.total_log + landcovers.5k + homogen + Predominant_land_use:fert.total_log + Predominant_land_use:homogen + Use_intensity:fert.total_log + Use_intensity:homogen",
               randomStruct = "(1|SS) + (1|SSB)", REML = TRUE)

summary(abmod$model)

# save the coefficents of the model
coefs <- fixef(abmod$model)


write.csv(coefs, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/3. Outputs/Stage_5/Stage5_final_models_ABUNDANCE_14102019_coefs.csv")




#### 3. RCAR 110km ####


system.time({sr3.1 <- GLMERSelect(modelData = model.data_RCAR_110km, responseVar = "RCAR_110km",
                   fitFamily = "gaussian", 
                   fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity"),
                   fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1),
                   randomStruct = "(1|SS)+(1|SSB)", 
                   fixedInteractions = c(#"Predominant_land_use:Use_intensity", 
                                         "Predominant_land_use:fert.total_log", 
                                         "Predominant_land_use:Hansen_mindist_log", 
                                         "Predominant_land_use:landcovers.5k", 
                                         "Predominant_land_use:homogen", 
                                         "Use_intensity:fert.total_log", 
                                         "Use_intensity:Hansen_mindist_log", 
                                         "Use_intensity:landcovers.5k",
                                         "Use_intensity:homogen"),
                   verbose = TRUE)}
)


#system.time({sr3.2 <- GLMERSelect(modelData = model.data_RCAR_110km, responseVar = "RCAR_110km",
#                                fitFamily = "gaussian", 
#                                fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity"),
#                                fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1),
#                                randomStruct = "(1|SS)+(1|SSB)", 
#                                fixedInteractions = c("Predominant_land_use:Use_intensity", 
#                                                      "Predominant_land_use:fert.total_log", 
#                                                      "Predominant_land_use:Hansen_mindist_log", 
#                                                      "Predominant_land_use:landcovers.5k", 
#                                                      "Predominant_land_use:homogen", 
#                                                      "Use_intensity:fert.total_log", 
#                                                      "Use_intensity:Hansen_mindist_log", 
#                                                      "Use_intensity:landcovers.5k",
#                                                      "Use_intensity:homogen"),
#                                verbose = TRUE)}
#)

sr3.1$model

save(sr3.1, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/3. Outputs/Stage_5/Stage5_final_models_RCAR_110km_14102019.rdata")



# RCAR_110km ~ Predominant_land_use + Forest_biome + Use_intensity +
# poly(homogen, 1) + Hansen_mindist_log + landcovers.5k + homogen +
# Predominant_land_use:homogen + Use_intensity:Hansen_mindist_log +  
# Use_intensity:landcovers.5k + Use_intensity:homogen +
# (1 | SS) + (1 | SSB)

rcarmod <- GLMER(modelData = model.data_RCAR_110km, responseVar = "RCAR_110km", fitFamily = "gaussian",
               fixedStruct = "Predominant_land_use + Forest_biome + Use_intensity + Hansen_mindist_log + homogen + landcovers.5k + Predominant_land_use:homogen + Use_intensity:Hansen_mindist_log + Use_intensity:landcovers.5k + Use_intensity:homogen",
               randomStruct = "(1|SS) + (1|SSB)", REML = TRUE)

summary(rcarmod$model)

# save the coefficents of the model
coefs <- fixef(rcarmod$model)


write.csv(coefs, file = "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/3. Outputs/Stage_5/Stage5_final_models_RCAR_nopoly_LU_14102019_coefs.csv")


##%######################################################%##
#                                                          #
####               Plots of Model Outputs               ####
#                                                          #
##%######################################################%##

# source in the edited plotting function
source("C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/BIOTA/Plotting_scripts/Tim's edited/PlotGLMERContinuous_edited.r")

plotdir <- "C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/1. ForestCover_Yield/3. Outputs/Stage_5/Figures"

################## Species Richness Model ##################



PlotGLMERFactor(model = srmod$model,data = srmod$data,
                responseVar = "Species richness",seMultiplier = 1,
                logLink = "e",catEffects = "Predominant_land_use", params = list(las = 2), xtext.srt = 15)

PlotGLMERFactor(model = srmod$model,data = srmod$data,
                responseVar = "Species richness",seMultiplier = 1,
                logLink = "e",catEffects = "Forest_biome", xtext.srt = 15)

PlotGLMERFactor(model = srmod$model,data = srmod$data,
                responseVar = "Species richness",seMultiplier = 1,
                logLink = "e",catEffects = "Use_intensity", xtext.srt = 15)


# distance only
pdf(file = paste0(plotdir, "/SR_distance_only.pdf"), width = 6, height = 4)
PlotGLMERContinuous_edit(model = srmod$model,
                         data = srmod$data, 
                         effects = "Hansen_mindist_log", 
                         otherContEffects = c("fert.total_log"),
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use"),
                         logLink = "e",
                         xlab ="Distance to forest (Km)", 
                         ylab = "Species Richness", 
                         plotRug = TRUE, 
                         line.cols = c("#458B00"),
                         seMultiplier = 1,
                         transformX = FALSE, 
                         rescaled = c(attr(Hansen_mindist_log, "scaled:scale"), attr(Hansen_mindist_log, "scaled:center")), 
                         logged = TRUE,
                         xlim = c(0, 20),
                         params = list(cex = 0.6))

dev.off()


# distance and land use
pdf(file = paste0(plotdir, "/SR_distance_LU.pdf"), width = 6, height = 4)
PlotGLMERContinuous_edit(model = srmod$model,
                         data = srmod$data, 
                         effects = "Hansen_mindist_log", 
                         otherContEffects = c("fert.total_log"),
                         byFactor = "Predominant_land_use",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use"),
                         logLink = "e",
                         xlab ="Distance to forest (Km)", 
                         ylab = "Species Richness", 
                         plotRug = TRUE, 
                         line.cols = c("#458B00", "#8B0000", "#FFA500"),
                         seMultiplier = 1,
                         transformX = FALSE, 
                         rescaled = c(attr(Hansen_mindist_log, "scaled:scale"), attr(Hansen_mindist_log, "scaled:center")), 
                         logged = TRUE,
                         xlim = c(0, 20),
                         params = list(cex = 0.6))
legend(15,18,c("Primary vegetation","Secondary vegetation", "Cropland"),
       col=c("#458B00", "#8B0000", "#FFA500"),bty="n",lty=1, cex = 1, y.intersp = 1)

dev.off()


# fertiliser only
pdf(file = paste0(plotdir, "/SR_fert_only.pdf"), width = 6, height = 4)
PlotGLMERContinuous_edit(model = srmod$model, data = srmod$data, 
                         effects = "fert.total_log", 
                         otherContEffects = c("Hansen_mindist_log"),
                         #byFactor = "Predominant_land_use",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use"),
                         logLink = "e",
                         xlab ="Total fertiliser application (kgs) /hec", 
                         ylab = "Species Richness", 
                         plotRug = TRUE, 
                         line.cols = c("#FFA500"),
                         seMultiplier = 1,
                         rescaled = c(attr(fert.total_log, "scaled:scale"), attr(fert.total_log, "scaled:center")),
                         logged = TRUE,
                         xlim = c(0, 800000),
                         params = list(cex = 0.6))

dev.off()


# fertiliser and land use
pdf(file = paste0(plotdir, "/SR_fert_LU.pdf"), width = 6, height = 4)
PlotGLMERContinuous_edit(model = srmod$model, data = srmod$data, 
                         effects = "fert.total_log", 
                         otherContEffects = c("Hansen_mindist_log"),
                         byFactor = "Predominant_land_use",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use"),
                         logLink = "e",
                         xlab ="Total fertiliser application (kgs) /hec", 
                         ylab = "Species Richness", 
                         plotRug = TRUE, 
                         line.cols = c("#458B00", "#8B0000", "#FFA500"),
                         seMultiplier = 1,
                         rescaled = c(attr(fert.total_log, "scaled:scale"), attr(fert.total_log, "scaled:center")),
                         logged = TRUE,
                         xlim = c(0, 800000),
                         params = list(cex = 0.6))
legend(600000,11,c("Primary vegetation","Secondary vegetation", "Cropland"),
       col=c("#458B00", "#8B0000", "#FFA500"),bty="n",lty=1, y.intersp = 1)

dev.off()


# fertiliser and intensity
pdf(file = paste0(plotdir, "/SR_fert_UI.pdf"), width = 6, height = 4)
PlotGLMERContinuous_edit(model = srmod$model, data = srmod$data, 
                         effects = "fert.total_log", 
                         otherContEffects = c("Hansen_mindist_log"),
                         byFactor = "Use_intensity",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use"),
                         logLink = "e",
                         xlab ="Total fertiliser application (kgs) /hec", 
                         ylab = "Species Richness", 
                         plotRug = T, 
                         line.cols = c("#458B00", "#8B0000", "#FFA500"),
                         seMultiplier = 1,
                         rescaled = c(attr(fert.total_log, "scaled:scale"), attr(fert.total_log, "scaled:center")),
                         logged = TRUE,
                         xlim = c(0, 800000),
                         params = list(cex = 0.6))
legend(650000,18,c("Minimal use", "Light use", "Intense use"),
       col=c("#458B00", "#8B0000", "#FFA500"),bty="n",lty=1, y.intersp = 1)

dev.off()




##################### Abundance Model  ##################### 

PlotGLMERFactor(model = abmod$model,data = abmod$data,
                responseVar = "Abundance",seMultiplier = 1,
                logLink = "e",catEffects = c("Predominant_land_use"), params = list(las = 2), xtext.srt = 15)

PlotGLMERFactor(model = abmod$model,data = abmod$data,
                responseVar = "Abundance",seMultiplier = 1,
                logLink = "e",catEffects = "Forest_biome", xtext.srt = 45)

PlotGLMERFactor(model = abmod$model,data = abmod$data,
                responseVar = "Abundance",seMultiplier = 1,
                logLink = "e",catEffects = "Use_intensity", xtext.srt = 45)

# Abundance, distance only
pdf(file = paste0(plotdir, "/Abun_dist_only.pdf"), width = 6, height = 4)
PlotGLMERContinuous_edit(model = abmod$model, data = abmod$data, 
                    effects = "Hansen_mindist_log", 
                    otherContEffects = c("fert.total_log", "landcovers.5k", "homogen"),
                    #byFactor = "Predominant_land_use",
                    otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Use_intensity = "Minimal use", Predominant_land_use = "Primary vegetation"),
                    logLink = "e",
                    xlab ="Distance to forest (Km)", 
                    ylab = "Abundance", 
                    plotRug = T, 
                    line.cols = c("#458B00", "#8B0000", "#FFA500"),
                    seMultiplier = 1,
                    logged = TRUE,
                    rescaled = c(attr(Hansen_mindist_log, "scaled:scale"), attr(Hansen_mindist_log, "scaled:center")),
                    params = list(cex = 0.6))
dev.off()

# Abundance, distance LU
pdf(file = paste0(plotdir, "/Abun_dist_LU.pdf"), width = 6, height = 4)
PlotGLMERContinuous_edit(model = abmod$model, data = abmod$data, 
                    effects = "Hansen_mindist_log", 
                    otherContEffects = c("fert.total_log", "landcovers.5k", "homogen"),
                    byFactor = "Predominant_land_use",
                    otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Use_intensity = "Minimal use", Predominant_land_use = "Primary vegetation"),
                    logLink = "e",
                    xlab ="Distance to forest (Km)", 
                    ylab = "Abundance", 
                    plotRug = T, 
                    line.cols = c("#458B00", "#8B0000", "#FFA500"),
                    seMultiplier = 1,
                    logged = TRUE,
                    rescaled = c(attr(Hansen_mindist_log, "scaled:scale"), attr(Hansen_mindist_log, "scaled:center")),
                    params = list(cex = 0.6))
legend(8.5,195,c("Primary vegetation", "Secondary vegetation","Cropland"),
       col=c("#458B00", "#8B0000", "#FFA500"),bty="n",lty=1)
dev.off()



# Abundance, fert only
pdf(file = paste0(plotdir, "/Abun_fert_only.pdf"), width = 6, height = 4)
PlotGLMERContinuous_edit(model = abmod$model, data = abmod$data, 
                    effects = "fert.total_log", 
                    otherContEffects = c("Hansen_mindist_log", "landcovers.5k", "homogen"),
                    #byFactor = "Predominant_land_use",
                    otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Use_intensity = "Minimal use", Predominant_land_use = "Primary vegetation"),
                    logLink = "e",
                    xlab ="Total fertiliser application (kgs) /hec", 
                    ylab = "Abundance", 
                    plotRug = T, 
                    line.cols = c("#FFA500"),
                    seMultiplier = 1,
                    rescaled = c(attr(fert.total_log, "scaled:scale"), attr(fert.total_log, "scaled:center")),
                    logged = TRUE,
                    params = list(cex = 0.6),
                    xlim = c(0, 2500000))

dev.off()


# Abundance, fert LU
pdf(file = paste0(plotdir, "/Abun_fert_LU.pdf"), width = 6, height = 4)
PlotGLMERContinuous_edit(model = abmod$model, data = abmod$data, 
                         effects = "fert.total_log", 
                         otherContEffects = c("Hansen_mindist_log", "landcovers.5k", "homogen"),
                         byFactor = "Predominant_land_use",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Use_intensity = "Minimal use", Predominant_land_use = "Primary vegetation"),
                         logLink = "e",
                         xlab ="Total fertiliser application (kgs) /hec", 
                         ylab = "Abundance", 
                         plotRug = T, 
                         line.cols = c("#458B00", "#8B0000", "#FFA500"),
                         seMultiplier = 1,
                         rescaled = c(attr(fert.total_log, "scaled:scale"), attr(fert.total_log, "scaled:center")),
                         logged = TRUE,
                         params = list(cex = 0.6 ),
                         xlim = c(0, 1000000))
legend(800000,200,c("Primary vegetation", "Secondary vegetation","Cropland"),
       col=c("#458B00", "#8B0000", "#FFA500"),bty="n",lty=1)
dev.off()



# Abundance, fert Intensity
pdf(file = paste0(plotdir, "/Abun_fert_only.pdf"), width = 6, height = 4)
PlotGLMERContinuous_edit(model = abmod$model, data = abmod$data, 
                         effects = "fert.total_log", 
                         otherContEffects = c("Hansen_mindist_log", "landcovers.5k", "homogen"),
                         byFactor = "Use_intensity",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Use_intensity = "Minimal use", Predominant_land_use = "Primary vegetation"),
                         logLink = "e",
                         xlab ="Total fertiliser application (kgs) /hec", 
                         ylab = "Abundance", 
                         plotRug = T, 
                         line.cols = c("#66CD00", "#FFA500", "#FF0000"),
                         seMultiplier = 1,
                         rescaled = c(attr(fert.total_log, "scaled:scale"), attr(fert.total_log, "scaled:center")),
                         logged = TRUE,
                         params = list(cex = 0.6),
                         xlim = c(0, 900000))
legend(750000,240,c("Minimal use", "Light use","Intense use"),
       col=c("#66CD00", "#FFA500", "#FF0000"),bty="n",lty=1)

dev.off()


# Abundance, homogen only
pdf(file = paste0(plotdir, "/Abun_homogen_only.pdf"), width = 6, height = 4)
PlotGLMERContinuous_edit(model = abmod$model, data = abmod$data, 
                         effects = "homogen", 
                         otherContEffects = c("Hansen_mindist_log", "landcovers.5k", "fert.total_log"),
                         #byFactor = "Use_intensity",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Use_intensity = "Minimal use", Predominant_land_use = "Primary vegetation"),
                         logLink = "e",
                         xlab ="Homogeneity", 
                         ylab = "Abundance", 
                         plotRug = T, 
                         line.cols = c("#66CD00", "#FFA500", "#FF0000"),
                         seMultiplier = 1,
                         rescaled = c(attr(homogen, "scaled:scale"), attr(homogen, "scaled:center")),
                         logged = F,
                         params = list(cex = 0.6),
                         xlim = c(0, 1))

dev.off()

pdf(file = paste0(plotdir, "/Abun_homogen_UI.pdf"), width = 6, height = 4)
PlotGLMERContinuous_edit(model = abmod$model, data = abmod$data, 
                         effects = "homogen", 
                         otherContEffects = c("Hansen_mindist_log", "landcovers.5k", "fert.total_log"),
                         byFactor = "Use_intensity",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Use_intensity = "Minimal use", Predominant_land_use = "Primary vegetation"),
                         logLink = "e",
                         xlab ="Homogeneity", 
                         ylab = "Abundance", 
                         plotRug = T, 
                         line.cols = c("#66CD00", "#FFA500", "#FF0000"),
                         seMultiplier = 1,
                         rescaled = c(attr(homogen, "scaled:scale"), attr(homogen, "scaled:center")),
                         logged = F,
                         params = list(cex = 0.6),
                         xlim = c(0, 1))
legend(0.8,240,c("Minimal use", "Light use","Intense use"),
       col=c("#66CD00", "#FFA500", "#FF0000"),bty="n",lty=1)

dev.off()


########################## RCAR Model ###################################
rcarmod$model


PlotGLMERFactor(model = rcarmod$model,data = rcarmod$data,
                responseVar = "RCAR_110km",seMultiplier = 1,
                logLink = "10",catEffects = c("Predominant_land_use"), params = list(las = 2), xtext.srt = 15)

PlotGLMERFactor(model = rcarmod$model,data = rcarmod$data,
                responseVar = "RCAR_110km",seMultiplier = 1,
                logLink = "10",catEffects = "Forest_biome", xtext.srt = 25)


PlotGLMERFactor(model = rcarmod$model,data = rcarmod$data,
                responseVar = "RCAR_110km",seMultiplier = 1,
                logLink = "10",catEffects = "Use_intensity", xtext.srt = 25)


# distance only
pdf(file = paste0(plotdir, "/RCAR_distance_only.pdf"), width = 6, height = 4)
PlotGLMERContinuous_edit(model = rcarmod$model, data = rcarmod$data, 
                         effects = "Hansen_mindist_log", 
                         otherContEffects = c("landcovers.5k", "homogen"),
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Use_intensity = "Minimal use", Predominant_land_use = "Primary vegetation"),
                         logLink = "10",
                         xlab ="Distance to forest (Km)", 
                         ylab = "RCAR", 
                         plotRug = TRUE, 
                         line.cols = c("#458B00"),
                         seMultiplier = 1, 
                         rescaled = c(attr(Hansen_mindist_log, "scaled:scale"), attr(Hansen_mindist_log, "scaled:center")), 
                         logged = TRUE,
                         xlim = c(0, 200), 
                         params = list(cex= 0.6, mgp=c(3.8,1,0)))

dev.off()


# distance and LU
pdf(file = paste0(plotdir, "/RCAR_distance_LU.pdf"), width = 6, height = 4)
PlotGLMERContinuous_edit(model = rcarmod$model, data = rcarmod$data, 
                    effects = "Hansen_mindist_log", 
                    otherContEffects = c("landcovers.5k", "homogen"),
                    byFactor = "Predominant_land_use",
                    otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Use_intensity = "Minimal use"),
                    logLink = "10",
                    xlab ="Distance to forest (Km)", 
                    ylab = "RCAR", 
                    plotRug = TRUE, 
                    line.cols = c("#458B00", "#8B0000", "#FFA500"),
                    seMultiplier = 1, 
                    rescaled = c(attr(Hansen_mindist_log, "scaled:scale"), attr(Hansen_mindist_log, "scaled:center")), 
                    logged = TRUE,
                    xlim = c(0, 15), 
                    params = list(cex= 0.6, mgp=c(3.8,1,0)))
legend(11,2100000,c("Primary vegetation","Secondary vegetation","Cropland"),
       col=c("#458B00", "#FFA500"),bty="n",lty=1)

dev.off()

# landcovers only
pdf(file = paste0(plotdir, "/RCAR_landcovers_only.pdf"), width = 6, height = 4)
PlotGLMERContinuous_edit(model = rcarmod$model, data = rcarmod$data, 
                         effects = "landcovers.5k", 
                         otherContEffects = c("Hansen_mindist_log", "homogen"),
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Use_intensity = "Minimal use", Predominant_land_use = "Primary vegetation"),
                         logLink = "10",
                         xlab ="Number of landcovers", 
                         ylab = "RCAR", 
                         plotRug = TRUE, 
                         line.cols = c("#1C86EE"),
                         seMultiplier = 1, 
                         rescaled = c(attr(landcovers.5k, "scaled:scale"), attr(landcovers.5k, "scaled:center")), 
                         logged = F,
                         xlim = c(0, 10), 
                         params = list(cex= 0.6, mgp=c(3.8,1,0)))

dev.off()

# homogen only
pdf(file = paste0(plotdir, "/RCAR_homogen_only.pdf"), width = 6, height = 4)
PlotGLMERContinuous_edit(model = rcarmod$model, data = rcarmod$data, 
                         effects = "homogen", 
                         otherContEffects = c("landcovers.5k", "Hansen_mindist_log"),
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Use_intensity = "Minimal use", Predominant_land_use = "Primary vegetation"),
                         logLink = "10",
                         xlab ="Homogeneity", 
                         ylab = "RCAR", 
                         plotRug = TRUE, 
                         line.cols = c("#458B00"),
                         seMultiplier = 1, 
                         rescaled = c(attr(homogen, "scaled:scale"), attr(homogen, "scaled:center")), 
                         logged = F,
                         xlim = c(0, 1), 
                         params = list(cex= 0.6, mgp=c(3.8,1,0)))

dev.off()


# homogen and LU
pdf(file = paste0(plotdir, "/RCAR_homogen_LU.pdf"), width = 6, height = 4)
PlotGLMERContinuous_edit(model = rcarmod$model, data = rcarmod$data, 
                         effects = "homogen", 
                         otherContEffects = c("landcovers.5k", "Hansen_mindist_log"),
                         byFactor = "Predominant_land_use",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Use_intensity = "Minimal use", Predominant_land_use = "Primary vegetation"),
                         logLink = "10",
                         xlab ="Homogeneity", 
                         ylab = "RCAR", 
                         plotRug = TRUE, 
                         line.cols = c("#458B00", "#8B0000", "#FFA500"),
                         seMultiplier = 1, 
                         rescaled = c(attr(homogen, "scaled:scale"), attr(homogen, "scaled:center")), 
                         logged = F,
                         xlim = c(0, 1), 
                         params = list(cex= 0.6, mgp=c(3.8,1,0)))
legend(0.7,1900000,c("Primary vegetation","Secondary vegetation","Cropland"),
       col=c("#458B00", "#FFA500"),bty="n",lty=1)

dev.off()


# distance and UI
pdf(file = paste0(plotdir, "/RCAR_distance_UI.pdf"), width = 6, height = 4)
PlotGLMERContinuous_edit(model = rcarmod$model, data = rcarmod$data, 
                         effects = "Hansen_mindist_log", 
                         otherContEffects = c("landcovers.5k", "homogen"),
                         byFactor = "Use_intensity",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Use_intensity = "Minimal use", Predominant_land_use = "Primary vegetation"),
                         logLink = "10",
                         xlab ="Distance to forest (Km)", 
                         ylab = "RCAR", 
                         plotRug = TRUE, 
                         line.cols = c("#66CD00", "#FFA500", "#FF0000"),
                         seMultiplier = 1, 
                         rescaled = c(attr(Hansen_mindist_log, "scaled:scale"), attr(Hansen_mindist_log, "scaled:center")), 
                         logged = TRUE,
                         xlim = c(0, 30), 
                         params = list(cex= 0.6, mgp=c(3.8,1,0)))
legend(22,1500000,c("Minimal use", "Light use","Intense use"),
       col=c("#66CD00", "#FFA500", "#FF0000"),bty="n",lty=1)

dev.off()

# landcovers and UI
pdf(file = paste0(plotdir, "/RCAR_landcovers_UI.pdf"), width = 6, height = 4)
PlotGLMERContinuous_edit(model = rcarmod$model, data = rcarmod$data, 
                         effects = "landcovers.5k", 
                         otherContEffects = c("homogen", "Hansen_mindist_log"),
                         byFactor = "Use_intensity",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Use_intensity = "Minimal use", Predominant_land_use = "Primary vegetation"),
                         logLink = "10",
                         xlab ="Homogeneity", 
                         ylab = "RCAR", 
                         plotRug = TRUE, 
                         line.cols = c("#66CD00", "#FFA500", "#FF0000"),
                         seMultiplier = 1, 
                         rescaled = c(attr(landcovers.5k, "scaled:scale"), attr(landcovers.5k, "scaled:center")), 
                         logged = F,
                         xlim = c(1, 10), 
                         params = list(cex= 0.6, mgp=c(3.8,1,0)))
legend(8,1700000,c("Minimal use", "Light use","Intense use"),
       col=c("#66CD00", "#FFA500", "#FF0000"),bty="n",lty=1)

dev.off()


# homogen and UI
pdf(file = paste0(plotdir, "/RCAR_homogen_UI.pdf"), width = 6, height = 4)
PlotGLMERContinuous_edit(model = rcarmod$model, data = rcarmod$data, 
                         effects = "homogen", 
                         otherContEffects = c("landcovers.5k", "Hansen_mindist_log"),
                         byFactor = "Use_intensity",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Use_intensity = "Minimal use", Predominant_land_use = "Primary vegetation"),
                         logLink = "10",
                         xlab ="Homogeneity", 
                         ylab = "RCAR", 
                         plotRug = TRUE, 
                         line.cols = c("#66CD00", "#FFA500", "#FF0000"),
                         seMultiplier = 1, 
                         rescaled = c(attr(homogen, "scaled:scale"), attr(homogen, "scaled:center")), 
                         logged = F,
                         xlim = c(0, 1), 
                         params = list(cex= 0.6, mgp=c(3.8,1,0)))
legend(0.8,1700000,c("Minimal use", "Light use","Intense use"),
       col=c("#66CD00", "#FFA500", "#FF0000"),bty="n",lty=1)

dev.off()



### Next, use the coeeficients save here and the global extent variable daata from 
### the python script to generate spatial projections.
# Projections_calc.r