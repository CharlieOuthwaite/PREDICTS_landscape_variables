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
datadir <- '1_PREDICTS_PLUS_VARIABLES'
moddir <- '2_MODEL_SELECTION'

# load model outputs - use the model selection versions as these include the stats in the model output
load(paste0(moddir, "/SPECIESRICHNESS_Tropical_Model_selection.rdata"))
load(paste0(moddir, "/SPECIESRICHNESS_Temperate_Model_selection.rdata"))

load(paste0(moddir, "/ABUNDANCE_Tropical_Model_Selection.rdata"))
load(paste0(moddir, "/ABUNDANCE_Temperate_Model_Selection.rdata"))


# read in the data files

# load in the datasets with transformed variables and the subsets 
load(paste0(datadir, "/PREDICTS_dataset_inc_variables_TRANS.rdata")) # final.data.trans


# Split the dataset based on realm

final.data.trans_trop <- final.data.trans[final.data.trans$Tropical == "Tropical", ]
nrow(final.data.trans_trop)
final.data.trans_trop <- droplevels(final.data.trans_trop)

final.data.trans_temp <- final.data.trans[final.data.trans$Tropical == "Temperate", ]
nrow(final.data.trans_temp)
final.data.trans_temp <- droplevels(final.data.trans_temp)


# separate out the data where abundance column is not NA
final.data.trans_trop_ABUN <- final.data.trans_trop[!is.na(final.data.trans_trop$Total_abundance), ] # 3314 rows

# separate out the data where abundance column is not NA
final.data.trans_temp_ABUN <- final.data.trans_temp[!is.na(final.data.trans_temp$Total_abundance), ] # 5740 rows


# load the scaling values table
scalers <- read.csv(paste0(datadir, "/Scaling_values.csv"))


# load the selected models
load(paste0(moddir, "/SRMOD_Tropical_output.rdata")) # srmod_trop
load(paste0(moddir, "/SRMOD_Temperate_output.rdata")) # srmod_temp

load(paste0(moddir, "/ABMOD_Tropical_output.rdata")) # abmod_temp
load(paste0(moddir, "/ABMOD_Temperate_output.rdata")) # abmod_trop




# 1. Richness model, land use and use intensity combo - TROPICAL


# basic table of median values and reference factors
pred_tab <- data.frame(landcovers.5k = median(final.data.trans_trop$landcovers.5k),
                       homogen = median(final.data.trans_trop$homogen),
                       fert.total_log = median(final.data.trans_trop$fert.total_log),
                       percNH = median(final.data.trans_trop$percNH),
                       Hansen_mindist_log =  median(final.data.trans_trop$Hansen_mindist_log),
                       Forest_biome = "Tropical & Subtropical Moist Broadleaf Forests",
                       Use_intensity = "Minimal use",
                       Predominant_land_use = "Primary vegetation",
                       #Tropical = "Temperate",
                       Species_richness = 0,
                       logAbun = 0)


levels(pred_tab$Predominant_land_use) <- levels(srmod_trop$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(srmod_trop$data$Use_intensity) 
levels(pred_tab$Forest_biome) <- levels(srmod_trop$data$Forest_biome)[c(3, 2, 1)]



# add and change factor levels of land use and intensity

pred_tab <- do.call("rbind", replicate(2, pred_tab, simplify = FALSE))


pred_tab[2, 'Predominant_land_use'] <- "Cropland"
pred_tab[2, 'Use_intensity'] <- "Intense use"

# predict the result
result <- PredictGLMER(model = srmod_trop$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result<- exp(result)

result$LU <- c("Primary", "Cropland")
result$UI <- c("Minimal", "Intense")

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

result$metric <- "SR"
result$realm <- "Tropical"

# save result table
write.csv(result, file = paste0(outdir, "/Predictions_Trop_richness.csv"), row.names = F)



# 2. Abundance model, non-tropical predictions

# basic table of median values and reference factors
pred_tab <- data.frame(landcovers.5k = median(final.data.trans_temp_ABUN$landcovers.5k),
                       homogen = median(final.data.trans_temp_ABUN$homogen),
                       fert.total_log = median(final.data.trans_temp_ABUN$fert.total_log),
                       percNH = median(final.data.trans_temp_ABUN$percNH),
                       Hansen_mindist_log =  median(final.data.trans_temp_ABUN$Hansen_mindist_log),
                       Forest_biome = "Temperate Broadleaf & Mixed Forests",
                       Use_intensity = "Minimal use",
                       Predominant_land_use = "Primary vegetation",
                       #Tropical = "Temperate",
                       Species_richness = 0,
                       logAbun = 0)


levels(pred_tab$Predominant_land_use) <- levels(abmod_temp$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(abmod_temp$data$Use_intensity) 
levels(pred_tab$Forest_biome) <- levels(abmod_temp$data$Forest_biome)


# add an extra row
pred_tab <- do.call("rbind", replicate(2, pred_tab, simplify = FALSE))

# firstly - comparing LU-UI combo
pred_tab[2, 'Predominant_land_use'] <- "Cropland"

## next looking at range of homogen
pred_tab[3:4, ] <- pred_tab[1, ]

# distance from forest to predict SR 
homogen1 <- 0.2 # rough range of data min
homogen2 <- 0.7 # rough range of data max

# rescale this as was done before modelling
#(log(x+1)-c)/s

homogen1_rescaled <- (homogen1-scalers[4,3])/scalers[4,2]
homogen2_rescaled <- (homogen2-scalers[4,3])/scalers[4,2]

pred_tab[3, 2] <- homogen1_rescaled
pred_tab[4, 2] <- homogen2_rescaled


## next looking at range of percNH
pred_tab[5:6, ] <- pred_tab[1, ]

# what is the range of interest
percNH1 <- 0 
percNH2 <- 100 

# rescale
percNH1_rescaled <- (percNH1-scalers[5,3])/scalers[5,2]
percNH2_rescaled <- (percNH2-scalers[5,3])/scalers[5,2]

# add into the predictin matrix
pred_tab[5, 4] <- percNH1_rescaled
pred_tab[6, 4] <- percNH2_rescaled

## next looking at range of fertiliser
pred_tab[7:8, ] <- pred_tab[1, ]

# what is the range of interest
fert1 <- 0 
fert2 <- 500

# rescale
fert1_rescaled <- (log(fert1+1)-scalers[2,3])/scalers[2,2]
fert2_rescaled <- (log(fert2+1)-scalers[2,3])/scalers[2,2]

# add into the predictin matrix
pred_tab[7, 3] <- fert1_rescaled
pred_tab[8, 3] <- fert2_rescaled

# change the use intensity
pred_tab[7, 7] <- "Intense use"
pred_tab[8, 7] <- "Intense use"


# predict the result
result <- PredictGLMER(model = abmod_temp$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result<- exp(result)-1

# label results table
result$test <- NA

# 1. LU UI test
result$test[1] <- "primary minimal"
result$test[2] <- "cropland minimal"

# 2.homogen test
result$test[3] <- "homogen 0.2"
result$test[4] <- "homogen 0.7"

# 3. percNH test
result$test[5] <- "percNH 0"
result$test[6] <- "percNH 100"

# 3. fert test
result$test[7] <- "fert 0"
result$test[8] <- "fert 500"

# calculate percentage change

# add column for results
result$change <- NA


# fill in differences
result[2,5] <- perc_change(result[1,1], result[2,1])
result[4,5] <- perc_change(result[3,1], result[4,1])
result[6,5] <- perc_change(result[5,1], result[6,1])
result[8,5] <- perc_change(result[7,1], result[8,1])

result$metric <- "Abun"
result$realm <- "Non_tropical"

# save result table
write.csv(result, file = paste0(outdir, "/Predictions_Temp_abundance.csv"), row.names = F)






# 3. Abundance model, tropical predictions

# basic table of median values and reference factors
pred_tab <- data.frame(landcovers.5k = median(final.data.trans_trop_ABUN$landcovers.5k),
                       homogen = median(final.data.trans_trop_ABUN$homogen),
                       fert.total_log = median(final.data.trans_trop_ABUN$fert.total_log),
                       percNH = median(final.data.trans_trop_ABUN$percNH),
                       Hansen_mindist_log =  median(final.data.trans_trop_ABUN$Hansen_mindist_log),
                       Forest_biome = "Tropical & Subtropical Moist Broadleaf Forests",
                       Use_intensity = "Minimal use",
                       Predominant_land_use = "Primary vegetation",
                       #Tropical = "Temperate",
                       Species_richness = 0,
                       logAbun = 0)


levels(pred_tab$Predominant_land_use) <- levels(abmod_trop$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(abmod_trop$data$Use_intensity) 
levels(pred_tab$Forest_biome) <- levels(abmod_trop$data$Forest_biome)[c(3, 2, 1)]




# add an extra row
pred_tab <- do.call("rbind", replicate(2, pred_tab, simplify = FALSE))

## next looking at range of homogen

# distance from forest to predict SR 
homogen1 <- 0.2 # rough range of data min
homogen2 <- 0.7 # rough range of data max

# rescale this as was done before modelling
#(log(x+1)-c)/s

homogen1_rescaled <- (homogen1-scalers[4,3])/scalers[4,2]
homogen2_rescaled <- (homogen2-scalers[4,3])/scalers[4,2]

pred_tab[1, 2] <- homogen1_rescaled
pred_tab[2, 2] <- homogen2_rescaled


## next looking at range of distance to forest in cropland
pred_tab[3:4, ] <- pred_tab[1, ]

pred_tab[3:4, 'Predominant_land_use'] <- "Cropland"

dist1 <- 0
dist2 <- 10

dist1_rescaled <- (log(dist1+1)-scalers[1,3])/scalers[1,2]
dist2_rescaled <- (log(dist2+1)-scalers[1,3])/scalers[1,2]

pred_tab[3, "Hansen_mindist_log"] <- dist1_rescaled
pred_tab[4, "Hansen_mindist_log"] <- dist2_rescaled



# predict the result
result <- PredictGLMER(model = abmod_trop$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result<- exp(result)-1

# label results table
result$test <- NA

# 1. LU UI test
result$test[1] <- "homogen 0.2"
result$test[2] <- "homogen 0.7"

# dist crop test
result$test[3] <- "dist 0 cropland"
result$test[4] <- "dist 10 cropland"



# calculate percentage change

# add column for results
result$change <- NA


# fill in differences
result[2,5] <- perc_change(result[1,1], result[2,1])
result[4,5] <- perc_change(result[3,1], result[4,1])


result$metric <- "Abun"
result$realm <- "Tropical"

# save result table
write.csv(result, file = paste0(outdir, "/Predictions_Trop_abundance.csv"), row.names = F)












##%######################################################%##
#                                                          #
####    old set of predictions before realm reruns      ####
#                                                          #
##%######################################################%##




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

