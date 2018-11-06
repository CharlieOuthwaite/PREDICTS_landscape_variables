############################################################
#                                                          #
#                   Running basic GLMMs                    #
#                                                          #
############################################################

# 6th November 2018

# Using first attempt at sorting the data to run some basic models and look at the outputs


rm(list = ls())

# where is the data saved?
datadir <- "D:/BIOTA/1_Forest_Cover_Yield/Data Exploration"

# read in the total yields dataset
totals <- read.csv(paste0(datadir, "/", "All_Crop_Yield.csv"))

