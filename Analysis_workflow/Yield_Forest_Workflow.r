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


# where are the various datasets saved?
datadir <- "D:/BIOTA/Data"


############################################################
#                                                          #
#              Step 1: organise PREDICTS data              #
#                                                          #
############################################################

# read in the complete PREDICTS dataset
pred.data <- read.csv(paste0(datadir, "/PREDICTS_2016/resource.csv"))

### organise using Tim's functions ###



# merge sites

# site level metrics

# reduce to cropland sites only

# maintain species level dataset as well, cropland sites only





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




