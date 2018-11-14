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
install_github("timnewbold/StatisticalModels")
install_github("timnewbold/predicts-demo",subdir="predictsFunctions")
library(predictsFunctions)
library(StatisticalModels)

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

# merge sites: this combines potential subsamples within one site
pred.sites <- MergeSites(pred.data) # 2906994 rows

# use correct sampling effort function (this replaces any with NA so the next function works)
pred.sites <- CorrectSamplingEffort(pred.sites)

# Calculate site level metrics
pred.sites.metrics <- SiteMetrics(pred.sites, extra.cols = "Predominant_land_use") # 22678 rows


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



############################################################
#                                                          #
#     Step 3: Calculate total crop yield for each site     #
#                                                          #
############################################################





############################################################
#                                                          #
#       step 4: Calculate distance to forest habitat       #
#                                                          #
############################################################






############################################################
#                                                          #
#         Step 5: Additional covariate information         #
#                                                          #
############################################################

# include information on human population density, distance to road??
# this was done in the land use paper.



############################################################
#                                                          #
#                    Step 6: run models                    #
#                                                          #
############################################################


# dif model options (species richness, abundance, presence/absence)


# use functions to run model selection process

# checks for overdispersion and spatial autocorrelation?



