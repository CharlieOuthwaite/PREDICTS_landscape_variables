##%######################################################%##
#                                                          #
####  Assessing comparable responses across predictors  ####
#                                                          #
##%######################################################%##

# in this script, I look at the range of values for each predictor
# and predict the biodiversity change across each of them, both
# continuous and categorical to give a comparable estimate of change
# due to each variable.

rm(list = ls())


# load packages
library(StatisticalModels)
library(ggplot2)
library(raster)

# directories
datadir <- "2_MODEL_SELECTION"
outdir <- "6_Range_predictions"
#dir.create(outdir)




##%######################################################%##
#                                                          #
####     1. Compare global and predicts data ranges     ####
#                                                          #
##%######################################################%##

# Local variables: land use, use intensity
# Landscape variables: distance to forest, percentage NH, total fertiliser, 
# nlandcovers, homogeneity

# load the PREDICTS dataset
load("1_PREDICTS_PLUS_VARIABLES/PREDICTS_dataset_inc_variables_TRANS.rdata") # transformed variables

range(final.data.trans$Hansen_mindist) # 0.0000 324.4697
range(final.data.trans$percNH) # 0.992 99.988
range(final.data.trans$fert.total) # 0 2314226
range(final.data.trans$landcovers.5k) # 1 10
range(final.data.trans$homogen) # 0.1448 0.6686

# these values are presented in the table detailing the datasets.


# load in global maps
total.fert <- raster("0_DATA/Earthstat_fert_total_newmethod.tif")
hans <- read.csv("0_DATA/hans_min_dist_80.csv")
homogen <- raster("0_DATA/Homogeneity_01_05_5km_uint16.tif")
homogen <- homogen*0.0001
percNH <- raster("0_DATA/PercentNatural.tif")
percNH <- percNH/10

# no value for furthest distance to forest as just did predicts sites
cellStats(percNH, stat = 'range') # 0 100
cellStats(total.fert, stat = 'range') # 2.127861e-16 4.283350e+06
# there is a max of 10 landcovers in the dataset
cellStats(homogen, stat = 'range') # 0.065 1.000



##%######################################################%##
#                                                          #
####             predict biodiversity loss              ####
#                                                          #
##%######################################################%##

# load in the models and for each predict the abundance and richness loss across 
# the range of values, while keeping the rest at the median value

# load the scaling values table
scalers <- read.csv("1_PREDICTS_PLUS_VARIABLES/Scaling_values.csv")


# load the models
load(file = paste0(datadir, "/ABMOD_Tropical_output.rdata"))
load(file = paste0(datadir, "/ABMOD_Temperate_output.rdata"))
load(file = paste0(datadir, "/SRMOD_Tropical_output.rdata"))
load(file = paste0(datadir, "/SRMOD_Temperate_output.rdata"))


#### create the table for predicting results from each model

# interested in how this dinner in cropland, so using Cropland as the set land use

# basic table of median values and reference factors
pred_tab <- data.frame(landcovers.5kRS = median(final.data.trans$landcovers.5kRS),
                       homogenRS = median(final.data.trans$homogenRS),
                       fert.total_logRS = median(final.data.trans$fert.total_logRS),
                       percNHRS = median(final.data.trans$percNHRS),
                       Hansen_mindist_logRS =  median(final.data.trans$Hansen_mindist_logRS),
                       Forest_biome = "Tropical & Subtropical Moist Broadleaf Forests",
                       Use_intensity = "Intense use",
                       Predominant_land_use = "Cropland",
                       Species_richness = 0,
                       logAbun = 0)


pred_tab$Predominant_land_use <- factor(pred_tab$Predominant_land_use, levels = levels(abmod_trop$data$Predominant_land_use))
pred_tab$Use_intensity <- factor(pred_tab$Use_intensity, levels = levels(abmod_trop$data$Use_intensity))
pred_tab$Forest_biome <- factor(pred_tab$Forest_biome, levels = levels(abmod_trop$data$Forest_biome)[c(3, 2, 1)])


# replicate the row
pred_tab <- do.call("rbind", replicate(14, pred_tab, simplify = FALSE))

# 1: primary vs cropland
pred_tab[1, 'Predominant_land_use'] <- "Primary vegetation"

# 2: minimal vs intense
pred_tab[3, 'Use_intensity'] <- "Minimal use"

# 3: fertiliser min vs fertiliser max

# fertiliser change
# what is the range of interest
fert1 <- min(final.data.trans$fert.total) # 0
fert2 <- max(final.data.trans$fert.total) # 2314226

# rescale
fert1_rescaled <- (log(fert1+1)-scalers[2,3])/scalers[2,2]
fert2_rescaled <- (log(fert2+1)-scalers[2,3])/scalers[2,2]

# add into the predictin matrix
pred_tab[5, 3] <- fert1_rescaled
pred_tab[6, 3] <- fert2_rescaled


# 4: distance min vs distance max

dist1 <- min(final.data.trans$Hansen_mindist) # 0
dist2 <- max(final.data.trans$Hansen_mindist) # 324.4697

# rescale
dist1_rescaled <- (log(dist1+1)-scalers[1,3])/scalers[1,2]
dist2_rescaled <- (log(dist2+1)-scalers[1,3])/scalers[1,2]

# add into the predictin matrix
pred_tab[7, 5] <- dist1_rescaled
pred_tab[8, 5] <- dist2_rescaled


# 5: percNH min vs percNH max

perc1 <- min(final.data.trans$percNH) # 0.992
perc2 <- max(final.data.trans$percNH) # 99.988

# rescale
perc1_rescaled <- (perc1-scalers[5,3])/scalers[5,2]
perc2_rescaled <- (perc2-scalers[5,3])/scalers[5,2]

# add into the predictin matrix
pred_tab[9, 4] <- perc1_rescaled
pred_tab[10, 4] <- perc2_rescaled


# 6: landcovers min vs landcovers max

landc1 <- min(final.data.trans$landcovers.5k) # 1
landc2 <- max(final.data.trans$landcovers.5k) # 10

# rescale
landc1_rescaled <- (landc1-scalers[3,3])/scalers[3,2]
landc2_rescaled <- (landc2-scalers[3,3])/scalers[3,2]

# add into the predictin matrix
pred_tab[11, 1] <- landc1_rescaled
pred_tab[12, 1] <- landc2_rescaled


# 7: homogen min vs homogen max

homo1 <- min(final.data.trans$homogen) # 0.1448
homo2 <- max(final.data.trans$homogen) # 0.6686

# rescale
homo1_rescaled <- (homo1-scalers[4,3])/scalers[4,2]
homo2_rescaled <- (homo2-scalers[4,3])/scalers[4,2]

# add into the predictin matrix
pred_tab[13, 2] <- homo1_rescaled
pred_tab[14, 2] <- homo2_rescaled


# predict the result
#result_srtrop <- PredictGLMER(model = srmod_trop$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)
#result_abtrop <- PredictGLMER(model = abmod_trop$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# predict the result
result_srtrop <- PredictGLMERRandIter(model = srmod_trop$model, data = pred_tab, nIters = 10000)
result_abtrop <- PredictGLMERRandIter(model = abmod_trop$model, data = pred_tab, nIters = 10000)


# alter the table for the temperate models
pred_tab_temp <- pred_tab

pred_tab_temp$Forest_biome <- "Temperate Broadleaf & Mixed Forests"
pred_tab_temp$Forest_biome <- factor(pred_tab_temp$Forest_biome, levels = levels(srmod_temp$data$Forest_biome))

# run temperate models
#result_srtemp <- PredictGLMER(model = srmod_temp$model, data = pred_tab_temp, se.fit = TRUE, seMultiplier = 1.96)
#result_abtemp <- PredictGLMER(model = abmod_temp$model, data = pred_tab_temp, se.fit = TRUE, seMultiplier = 1.96)

# run temperate models
result_srtemp <- PredictGLMERRandIter(model = srmod_temp$model, data = pred_tab_temp, nIters = 10000)
result_abtemp <- PredictGLMERRandIter(model = abmod_temp$model, data = pred_tab_temp, nIters = 10000)



# transform the results
result_srtrop<- exp(result_srtrop)
result_srtemp<- exp(result_srtemp)
result_abtrop<- exp(result_abtrop)-1
result_abtemp<- exp(result_abtemp)-1


# calculate the percentage differences between iterations


# function to get percentage change from starting and ending values
perc_change <- function(start, end){
  
  cng <- round(((end-start)/start)*100, 2)
  return(cng)
  
}

# combine results
allresults <- as.data.frame(rbind(result_srtrop, result_srtemp, result_abtrop, result_abtemp))

dif_tab <-  NULL


# fill in differences
for(i in 1:nrow(allresults)){
  
  if(i %% 2 == 0){
    
    # determine the percentage change across all iters
    chng <- perc_change(allresults[i-1,], allresults[i,])
    
    med_change <- median(as.numeric(chng))
    LCI_change <- quantile(chng, probs = 0.025)
    UCI_change <- quantile(chng, probs = 0.975)
    
    dif_tab<- rbind(dif_tab, c(med_change, LCI_change, UCI_change))
    

    # changes_abtrop[i, 1] <- median(perc_change(result_abtrop[i-1,], result_abtrop[i,]))
    # changes_abtrop[i, 2] <- quantile(perc_change(result_abtrop[i-1,], result_abtrop[i,]), probs = 0.025)
    # changes_abtrop[i, 3] <- quantile(perc_change(result_abtrop[i-1,], result_abtrop[i,]), probs = 0.975)
    
    
    
  }else{next}
  
  
}


# organise table

dif_tab2 <- as.data.frame(dif_tab)
dif_tab2 <- do.call(rbind.data.frame, dif_tab2)


colnames(dif_tab2) <- c("median", "lower", "upper")

dif_tab2$median <- as.numeric(dif_tab2$median)
dif_tab2$lower <- as.numeric(dif_tab2$lower)
dif_tab2$upper <- as.numeric(dif_tab2$upper)

# add details of tests
tests <- c("Intense Primary to Intense Cropland", 
           "Minimal cropland to Intense cropland",
           "Fertiliser range", 
           "Distance range", 
           "Percentage NH range", 
           "Landcovers range", 
           "Homogeneity range")

dif_tab2$test <- tests

dif_tab2$model <- c(rep("srmod_trop", 7), rep("srmod_temp", 7), rep("abmod_trop", 7), rep("abmod_temp", 7))


# save the results table
write.csv(dif_tab2, paste0(outdir, "/DifferencesTable_nogroup.csv"), row.names = F)



##%######################################################%##
#                                                          #
####        Create figure to present differences        ####
#                                                          #
##%######################################################%##


#dif_tab <- read.csv(paste0(outdir, "/DifferencesTable_nogroup.csv"))


library(ggplot2)



dif_tab$trop <- sub(".*_", "", dif_tab$model)
dif_tab$metric <- sub("_.*", "", dif_tab$model)

dif_tab$trop <- sub("trop", "Tropical", dif_tab$trop)
dif_tab$trop <- sub("temp", "Non-tropical", dif_tab$trop)

dif_tab$metric <- sub("abmod", "Total Abundance", dif_tab$metric)
dif_tab$metric <- sub("srmod", "Species Richness", dif_tab$metric)

# rename info
dif_tab$test  <- sub(" range", "", dif_tab$test)
dif_tab$test  <- sub("Minimal cropland to Intense cropland", "Use Intensity", dif_tab$test)
dif_tab$test  <- sub("Intense Primary to Intense Cropland", "Land Use", dif_tab$test)
dif_tab$test  <- sub("Percentage NH", "Percentage Natural\n Habitat", dif_tab$test)
dif_tab$test  <- sub("Landcovers", "Number of\n Landcovers", dif_tab$test)
dif_tab$test  <- sub("Fertiliser", "Total Fertiliser", dif_tab$test)
dif_tab$test  <- sub("Distance", "Distance to Forest", dif_tab$test)

dif_tab$test <- factor(dif_tab$test, levels = rev(c("Land Use", "Use Intensity", "Distance to Forest",
                                                "Percentage Natural\n Habitat", "Total Fertiliser",
                                                "Number of\n Landcovers", "Homogeneity")))

dif_tab$trop <- factor(dif_tab$trop, levels = c("Tropical", "Non-tropical"))

dif_tab$median[dif_tab$median == 0] <- NA
dif_tab$lower[dif_tab$lower == 0] <- NA
dif_tab$upper[dif_tab$upper == 0] <- NA


ggplot(data = dif_tab) + 
  geom_point(aes(x = test, y = median, shape = trop, col = trop), position = position_dodge(width = 0.9), size = 2.5) +
  geom_errorbar(aes(x = test, y = median, ymin = lower, ymax = upper, col = trop), 
                position = position_dodge2(padding = 0.5),
                size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ metric) + coord_flip() +
  ylab("Percentage Change") + 
  scale_color_manual(values = c("#66CDAA", "#27408B")) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        panel.grid = element_blank(),
        legend.text = element_text(size = 10),
        aspect.ratio = 1, legend.background = element_blank(),
        text = element_text(size = 12), 
        axis.title.y = element_blank(),
        legend.position = "bottom")


ggsave(filename = paste0(outdir, "/Difference_plot.pdf"))  


