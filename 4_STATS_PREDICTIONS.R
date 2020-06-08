##%######################################################%##
#                                                          #
####               Stats and predictions                ####
#                                                          #
##%######################################################%##

# This script uses the various model outputs to extract and/or project the 
# results presented in the text of the paper.


rm(list = ls())

setwd("C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/BIOTA")

# load libraries
library(StatisticalModels)


# directories
datadir <- '2_MODEL_SELECTION'
outdir <- '4_STATS_PROJECTIONS'  


# load model outputs - use the model selection versions as these include the stats in the model output
load(paste0(datadir, "/SPECIESRICHNESS_Model_selection.rdata"))
load(paste0(datadir, "/ABUNDANCE_Model_Selection.rdata"))
load(paste0(datadir, "/RCAR_Model_Selection.rdata"))


# check correct models
sr1$model
ab1$model
rcar1$model



##%######################################################%##
#                                                          #
####                       Stats                        ####
#                                                          #
##%######################################################%##

# Look at the model selection stats to get the significance values for each effect


# extract and save the stats tables produced as a part of the GLMERSelect function
sr1stats <- as.data.frame(sr1$stats)
ab1stats <- as.data.frame(ab1$stats)
rcar1stats <- as.data.frame(rcar1$stats)

# add a column to check significance
sr1stats$significant <- NA
ab1stats$significant <- NA
rcar1stats$significant <- NA


# function to check significance
checksig <- function(x){
  if(x <= 0.05){ 
    res <- "Yes" 
  } else { 
    res <- "No" }
  return(res)}

# add values to table
sr1stats$significant <- sapply(X = sr1stats$P, FUN = checksig)
ab1stats$significant <- sapply(X = ab1stats$P, FUN = checksig)
rcar1stats$significant <- sapply(X = rcar1stats$P, FUN = checksig)


# save the stats tables
write.csv(sr1stats, file = paste0(outdir, "/SR_Stats.csv"), row.names = FALSE)
write.csv(ab1stats, file = paste0(outdir, "/Abun_Stats.csv"), row.names = FALSE)
write.csv(rcar1stats, file = paste0(outdir, "/RCAR_Stats.csv"), row.names = FALSE)




##%######################################################%##
#                                                          #
####                    Projections                     ####
#                                                          #
##%######################################################%##


################## Distance to forest ######################


### 1. Distance to forest, SR


# read in the data files

datadir <- '1_PREDICTS_PLUS_VARIABLES'
moddir <- '2_MODEL_SELECTION'

# load in the datasets with transformed variables and the subsets 
load(paste0(datadir, "/PREDICTS_dataset_inc_variables_TRANS.rdata")) # final.data.trans
load(paste0(datadir, "/PREDICTS_abun_subset.rdata")) # final.data.abun
load(paste0(datadir, "/PREDICTS_rcar_subset.rdata")) # final.data.rcar


# load the scaling values table
scalers <- read.csv(paste0(datadir, "/Scaling_values.csv"))


# load the selected models
load(paste0(moddir, "/SRMOD_output.rdata")) # srmod
load(paste0(moddir, "/ABMOD_output.rdata")) # abmod
load(paste0(moddir, "/RCAR_output.rdata")) # rcarmod

# distance from forest to predict SR 
dist0 <- 0 # also at 0 to calc difference
dist25 <- 25 # 25km, from looking at figure. 

# rescale this as was done before modelling
#(log(x+1)-c)/s

dist25_rescaled <- (log(dist25+1)-scalers[1,3])/scalers[1,2]
dist0_rescaled <- (log(dist0+1)-scalers[1,3])/scalers[1,2]


### get median values of other variables, use the transformed variables
landcovers.5k <- median(final.data.trans$landcovers.5k)
homogen <- median(final.data.trans$homogen)
fert.total_log <- median(final.data.trans$fert.total_log)
percNH <- median(final.data.trans$percNH)

# factors required too
Forest_biome <- "Temperate Broadleaf & Mixed Forests"
Use_intensity <- "Minimal use"
Predominant_land_use <- "Primary vegetation"
Tropical <- "Temperate"


# combine into data table
pred_tab <- cbind(Predominant_land_use, Use_intensity, Forest_biome, Tropical,
                  landcovers.5k, homogen, fert.total_log, percNH)

pred_tab <- as.data.frame(rbind(pred_tab, pred_tab)) # make 2 rows for the two distances

# the distance perameter
Hansen_mindist_log <- c(dist0_rescaled, dist25_rescaled)

# add distances to table
pred_tab <- cbind(pred_tab, Hansen_mindist_log)

pred_tab[,5:8] <- t(apply(X = pred_tab[,5:8], MARGIN = 1, FUN = function(x){as.numeric(as.character(x))}))


# check levels of factor variables
levels(pred_tab$Predominant_land_use) <- levels(srmod$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(srmod$data$Use_intensity) 
levels(pred_tab$Forest_biome) <- levels(srmod$data$Forest_biome) 
levels(pred_tab$Tropical) <- levels(srmod$data$Tropical) 

# give a value of species richness
pred_tab <- cbind(pred_tab, Species_richness = rep(0,2))

# predict the SR
result <- PredictGLMER(model = srmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# back-transform the result 
result <- exp(result)

result$distance <- c(0, 25)
result$metric <- rep("SR",2)


# 1. Distance to forest, abundance

# remove the SR subset based estimates
pred_tab[, 5:8] <- NA


### get median values of other variables, use the transformed variables
pred_tab$landcovers.5k <- rep(median(final.data.abun$landcovers.5k), 2)
pred_tab$homogen <- rep(median(final.data.abun$homogen), 2)
pred_tab$fert.total_log <- rep(median(final.data.abun$fert.total_log), 2)
pred_tab$percNH <- rep(median(final.data.abun$percNH), 2)

# Add in transformed 
pred_tab$logAbun <- rep(log(0+1), 2)


# predict abundance 
result_ab <- PredictGLMER(model = abmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform result
result_ab <- exp(result_ab)

result_ab <- cbind(result_ab, distance = c(0,25), metric = "Abun")

# combine results tables  
result <- rbind(result, result_ab)



# 1. distance to forest, RCAR

# remove the abun subset based estimates
pred_tab[, 5:8] <- NA


pred_tab$landcovers.5k <- rep(median(final.data.rcar$landcovers.5k), 2)
pred_tab$homogen <- rep(median(final.data.rcar$homogen), 2)
pred_tab$fert.total_log <- rep(median(final.data.rcar$fert.total_log), 2)
pred_tab$percNH <- rep(median(final.data.rcar$percNH), 2)


# create the req RCAR data
pred_tab$RCAR_110km <- rep(0, 2)

# predict the RCAR 
result3 <- PredictGLMER(model = rcarmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform result
result3 <- 10^(result3)

result3 <- cbind(result3, distance = c(0,25), metric = "RCAR")

# combine results tables  
result <- rbind(result, result3)

# round to 2 dp
result[,1:3] <- round(result[,1:3], digits = 2)


# calculate percentage change

# function to get percentage change from starting and ending values
perc_change <- function(start, end){
  
  cng <- round(((end-start)/start)*100, 2)
  return(cng)
  
}


# add column for results
result$change <- NA

# fill in differences
result[2,6] <- perc_change(result[1,1], result[2,1])
result[4,6] <- perc_change(result[3,1], result[4,1])
result[6,6] <- perc_change(result[5,1], result[6,1])


# save result table
write.csv(result, file = paste0(outdir, "/Predictions_dist.csv"), row.names = F)



################ Percentage natural habitat ####################

# 1. SR

# percNH values
percNH0 <- 0
percNH100 <- 100

# rescale
perc0_rescaled <- (percNH0-scalers[5,3])/scalers[5,2]
perc100_rescaled<- (percNH100-scalers[5,3])/scalers[5,2]


# get median values of other values, use the transformed variables
landcovers.5k <- rep(median(final.data.trans$landcovers.5k), 2)
homogen <- rep(median(final.data.trans$homogen), 2)
fert.total_log <- rep(median(final.data.trans$fert.total_log), 2)
Hansen_mindist_log <- rep(median(final.data.trans$Hansen_mindist_log), 2)

# factors required too?
Forest_biome <- rep("Temperate Broadleaf & Mixed Forests", 2)
Use_intensity <- rep("Minimal use", 2)
Predominant_land_use <- rep("Primary vegetation", 2)
Tropical <- rep("Temperate", 2)

percNH <- c(perc0_rescaled, perc100_rescaled)

# combine into data table
pred_tab <- cbind(Predominant_land_use, Use_intensity, Forest_biome, Tropical,
                  landcovers.5k, homogen, fert.total_log,
                  Hansen_mindist_log, percNH)

pred_tab <- as.data.frame(pred_tab)


pred_tab[,5:9] <- t(apply(X = pred_tab[,5:9], MARGIN = 1, FUN = function(x){as.numeric(as.character(x))}))

# set levels of factor variables
levels(pred_tab$Predominant_land_use) <- levels(srmod$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(srmod$data$Use_intensity)
levels(pred_tab$Forest_biome) <- levels(srmod$data$Forest_biome)
levels(pred_tab$Tropical) <- levels(srmod$data$Tropical)

# doesn't matter what these values are, just need to have the columns in there
pred_tab <- cbind(pred_tab, Species_richness = rep(0, 2))

# predict the SR
result <- PredictGLMER(model = srmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# back-transform the result 
result <- exp(result)

result$percNH <- c(0, 100)
result$metric <- rep("SR", 2)


# 2. abundance

# remove the SR subset based estimates
pred_tab[, 5:8] <- NA


### get median values of other variables, use the transformed variables
pred_tab$landcovers.5k <- rep(median(final.data.abun$landcovers.5k), 2)
pred_tab$homogen <- rep(median(final.data.abun$homogen), 2)
pred_tab$fert.total_log <- rep(median(final.data.abun$fert.total_log), 2)
pred_tab$Hansen_mindist_log <- rep(median(final.data.abun$Hansen_mindist_log), 2)


# Add in transformed 
pred_tab$logAbun <- rep(log(0+1), 2)


# predict abundance 
result_ab <- PredictGLMER(model = abmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform result
result_ab <- exp(result_ab)

result_ab <- cbind(result_ab, percNH = c(0,100), metric = "Abun")

# combine results tables  
result <- rbind(result, result_ab)



# 3. RCAR


### get median values of other variables, use the transformed variables
pred_tab$landcovers.5k <- rep(median(final.data.rcar$landcovers.5k), 2)
pred_tab$homogen <- rep(median(final.data.rcar$homogen), 2)
pred_tab$fert.total_log <- rep(median(final.data.rcar$fert.total_log), 2)
pred_tab$Hansen_mindist_log <- rep(median(final.data.rcar$Hansen_mindist_log), 2)

# create the req RCAR data
pred_tab$RCAR_110km <- rep(0, 2)

# predict the RCAR 
result3 <- PredictGLMER(model = rcarmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform result
result3 <- 10^(result3)

result3 <- cbind(result3, percNH = c(0,100), metric = "RCAR")

# combine results tables  
result <- rbind(result, result3)

# round to 2 dp
result[,1:3] <- round(result[,1:3], digits = 2)


# add column for results
result$change <- NA

# fill in differences
result[2,6] <- perc_change(result[1,1], result[2,1])
result[4,6] <- perc_change(result[3,1], result[4,1])
result[6,6] <- perc_change(result[5,1], result[6,1])


# save result table
write.csv(result, file = paste0(outdir, "/Predictions_percNH.csv"), row.names = F)



##%######################################################%##
#                                                          #
####                     fertiliser                     ####
#                                                          #
##%######################################################%##


# 1. Species richness

# fert values for predicting
fert0 <- 0
fert250 <- 250

# rescale
fert0_rescaled <- (log(fert0+1)-scalers[2,3])/scalers[2,2]
fert250_rescaled <- (log(fert250+1)-scalers[2,3]+1)/scalers[2,2]


# get median values of other values, use the transformed variables
landcovers.5k <- rep(median(final.data.trans$landcovers.5k), 2)
homogen <- rep(median(final.data.trans$homogen), 2)
percNH <- rep(median(final.data.trans$percNH), 2)
Hansen_mindist_log <- rep(median(final.data.trans$Hansen_mindist_log), 2)

# factors required too?
Forest_biome <- rep("Temperate Broadleaf & Mixed Forests", 2)
Use_intensity <- rep("Minimal use", 2)
Predominant_land_use <- rep("Cropland", 2) # looking at cropland here
Tropical <- rep("Temperate", 2)

fert.total_log <- c(fert0_rescaled, fert250_rescaled)

# combine into data table
pred_tab <- cbind(Predominant_land_use, Use_intensity, Forest_biome, Tropical,
                  landcovers.5k, homogen, Hansen_mindist_log, percNH, fert.total_log)

pred_tab <- as.data.frame(pred_tab)


pred_tab[,5:9] <- t(apply(X = pred_tab[,5:9], MARGIN = 1, FUN = function(x){as.numeric(as.character(x))}))

# set levels of factor variables
levels(pred_tab$Predominant_land_use) <- levels(srmod$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(srmod$data$Use_intensity)
levels(pred_tab$Forest_biome) <- levels(srmod$data$Forest_biome)
levels(pred_tab$Tropical) <- levels(srmod$data$Tropical)

# doesn't matter what these values are, just need to have the columns in there
pred_tab <- cbind(pred_tab, Species_richness = rep(0, 2))

pred_tab$Predominant_land_use <- factor(x ="Cropland", levels = levels(srmod$data$Predominant_land_use))

# predict the SR
result <- PredictGLMER(model = srmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# back-transform the result 
result <- exp(result)

result$fert.total_log <- c(0, 250)
result$metric <- rep("SR", 2)


# 2. abundance


# remove the SR subset based estimates
pred_tab[, 5:8] <- NA


### get median values of other variables, use the transformed variables
pred_tab$landcovers.5k <- rep(median(final.data.abun$landcovers.5k), 2)
pred_tab$homogen <- rep(median(final.data.abun$homogen), 2)
pred_tab$percNH <- rep(median(final.data.abun$percNH), 2)
pred_tab$Hansen_mindist_log <- rep(median(final.data.abun$Hansen_mindist_log), 2)


# Add in transformed 
pred_tab$logAbun <- rep(log(0+1), 2)


# predict abundance 
result_ab <- PredictGLMER(model = abmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform result
result_ab <- exp(result_ab)

result_ab <- cbind(result_ab, fert.total_log = c(0,250), metric = "Abun")

# combine results tables  
result <- rbind(result, result_ab)


# round to 2 dp
result[,1:3] <- round(result[,1:3], digits = 2)


# add column for results
result$change <- NA

# fill in differences
result[2,6] <- perc_change(result[1,1], result[2,1])
result[4,6] <- perc_change(result[3,1], result[4,1])


# save result table
write.csv(result, file = paste0(outdir, "/Predictions_fert.csv"), row.names = F)

