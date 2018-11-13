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



