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
####             Project biodiversity loss              ####
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

theme_custom <- theme(panel.grid = element_blank(),
                      legend.position = c(0.8,0.8), legend.title = element_blank(),
                      legend.text = element_text(size = 10),
                      axis.text = element_text(size = 10),
                      axis.title = element_text(size = 10),
                      aspect.ratio = 1, legend.background = element_blank(),
                      #text = element_text(size = 8), 
                      line = element_line(size = 0.2), 
                      panel.border = element_rect(size = 0.2),
                      strip.background = element_rect(size = 0.2),
                      axis.ticks = element_line(size = 0.2),
                      axis.title.y = element_blank())


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

dif_tab$metric <- factor(dif_tab$metric, levels = c("Total Abundance", "Species Richness"))


ggplot(data = dif_tab) + 
  geom_point(aes(x = test, y = median, shape = trop, col = trop), position = position_dodge(width = 0.9), size = 2) +
  geom_errorbar(aes(x = test, y = median, ymin = lower, ymax = upper, col = trop), 
                position = position_dodge2(padding = 0.5),
                size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  facet_wrap(~ metric) + coord_flip() +
  ylab("Percentage Change") + 
  scale_color_manual(values = c("#66CDAA", "#27408B")) +
  theme_bw() +
  theme_custom +
  theme(strip.background = element_rect(fill = NA))



ggsave(filename = paste0(outdir, "/Difference_plot.pdf"), width = 6.5, height = 3.5, unit = "in")  





##%######################################################%##
#                                                          #
####       Interactions differences across range        ####
#                                                          #
##%######################################################%##

# as with the above, look at differences across the range of predictor 
# values but for interactions across land use or use intensity where the 
# interaction was signigicant in model selection process. 


#### 1. tropical results first ####

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


# 1: distance across land covers (abundance)
pred_tab1 <- do.call("rbind", replicate(6, pred_tab, simplify = FALSE))
pred_tab1$Predominant_land_use[1:4] <- c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation")
pred_tab1$Hansen_mindist_logRS[c(1,3,5)] <- min(final.data.trans$Hansen_mindist_logRS)
pred_tab1$Hansen_mindist_logRS[c(2,4,6)] <- max(final.data.trans$Hansen_mindist_logRS)

pred_tab1$metric <- "abun"
pred_tab1$test <- "Distance to forest"
pred_tab1$cat <- "Land Use"
pred_tab1$realm <- "Tropical"


# 2: Landcovers across use intensities (abundance)
pred_tab2 <- do.call("rbind", replicate(6, pred_tab, simplify = FALSE))
pred_tab2$Use_intensity[1:4] <- c("Minimal use", "Minimal use", "Light use", "Light use")
pred_tab2$landcovers.5kRS[c(1,3,5)] <- min(final.data.trans$landcovers.5kRS)
pred_tab2$landcovers.5kRS[c(2,4,6)] <- max(final.data.trans$landcovers.5kRS)

pred_tab2$metric <- "abun"
pred_tab2$test <- "Number of Landcovers"
pred_tab2$cat <- "Use Intensity"
pred_tab2$realm <- "Tropical"


# 3: percnh across use intensities and landcovers (richness)
pred_tab3 <- do.call("rbind", replicate(12, pred_tab, simplify = FALSE))
pred_tab3$Use_intensity[1:4] <- c("Minimal use", "Minimal use", "Light use", "Light use")
pred_tab3$Predominant_land_use[7:10] <- c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation")
pred_tab3$percNHRS[c(1,3,5,7,9,11)] <- min(final.data.trans$percNHRS)
pred_tab3$percNHRS[c(2,4,6,8,10,12)] <- max(final.data.trans$percNHRS)

pred_tab3$metric <- "rich"
pred_tab3$test <- "Percentage Natural Habitat"
pred_tab3$cat <- NA
pred_tab3$cat[1:6] <- "Use Intensity"
pred_tab3$cat[7:12] <- "Land Use"
pred_tab3$realm <- "Tropical"


# 4: Fertiliser across land uses (richness)
pred_tab4 <- do.call("rbind", replicate(6, pred_tab, simplify = FALSE))
pred_tab4$Predominant_land_use[1:4] <- c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation")
pred_tab4$fert.total_logRS[c(1,3,5)] <- min(final.data.trans$fert.total_logRS)
pred_tab4$fert.total_logRS[c(2,4,6)] <- max(final.data.trans$fert.total_logRS)

pred_tab4$metric <- "rich"
pred_tab4$test <- "Total Fertiliser"
pred_tab4$cat <- "Land Use"
pred_tab4$realm <- "Tropical"


# 5: Landcovers across use intensities (richness)
pred_tab5 <- do.call("rbind", replicate(6, pred_tab, simplify = FALSE))
pred_tab5$Use_intensity[1:4] <- c("Minimal use", "Minimal use", "Light use", "Light use")
pred_tab5$landcovers.5kRS[c(1,3,5)] <- min(final.data.trans$landcovers.5kRS)
pred_tab5$landcovers.5kRS[c(2,4,6)] <- max(final.data.trans$landcovers.5kRS)

pred_tab5$metric <- "rich"
pred_tab5$test <- "Number of Landcovers"
pred_tab5$cat <- "Use Intensity"
pred_tab5$realm <- "Tropical"


# 6: Homogeneity across land uses (richess)
pred_tab6 <- do.call("rbind", replicate(6, pred_tab, simplify = FALSE))
pred_tab6$Predominant_land_use[1:4] <- c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation")
pred_tab6$homogenRS[c(1,3,5)] <- min(final.data.trans$homogenRS)
pred_tab6$homogenRS[c(2,4,6)] <- max(final.data.trans$homogenRS)

pred_tab6$metric <- "rich"
pred_tab6$test <- "Homogeneity"
pred_tab6$cat <- "Land Use"
pred_tab6$realm <- "Tropical"


# bring together table for abundance  projections

pred_tab_abtrop <- rbind(pred_tab1, pred_tab2)
pred_tab_srtrop <- rbind(pred_tab3, pred_tab4, pred_tab5, pred_tab6)


#### 2. temperate results next ####



# basic table of median values and reference factors
pred_tab <- data.frame(landcovers.5kRS = median(final.data.trans$landcovers.5kRS),
                       homogenRS = median(final.data.trans$homogenRS),
                       fert.total_logRS = median(final.data.trans$fert.total_logRS),
                       percNHRS = median(final.data.trans$percNHRS),
                       Hansen_mindist_logRS =  median(final.data.trans$Hansen_mindist_logRS),
                       Forest_biome = "Temperate Broadleaf & Mixed Forests",
                       Use_intensity = "Intense use",
                       Predominant_land_use = "Cropland",
                       Species_richness = 0,
                       logAbun = 0)


pred_tab$Predominant_land_use <- factor(pred_tab$Predominant_land_use, levels = levels(abmod_trop$data$Predominant_land_use))
pred_tab$Use_intensity <- factor(pred_tab$Use_intensity, levels = levels(abmod_trop$data$Use_intensity))
pred_tab$Forest_biome <- factor(pred_tab$Forest_biome, levels = levels(srmod_temp$data$Forest_biome))



# 1: percnh across use intensities (richness)
pred_tab1 <- do.call("rbind", replicate(6, pred_tab, simplify = FALSE))
pred_tab1$Use_intensity[1:4] <- c("Minimal use", "Minimal use", "Light use", "Light use")
pred_tab1$percNHRS[c(1,3,5)] <- min(final.data.trans$percNHRS)
pred_tab1$percNHRS[c(2,4,6)] <- max(final.data.trans$percNHRS)

pred_tab1$metric <- "rich"
pred_tab1$test <- "Percentage Natural Habitat"
pred_tab1$cat<- "Use Intensity"
pred_tab1$realm <- "Non-tropical"



# 2: Fertiliser across use intensities (richness)
pred_tab2 <- do.call("rbind", replicate(6, pred_tab, simplify = FALSE))
pred_tab2$Use_intensity[1:4] <- c("Minimal use", "Minimal use", "Light use", "Light use")
pred_tab2$fert.total_logRS[c(1,3,5)] <- min(final.data.trans$fert.total_logRS)
pred_tab2$fert.total_logRS[c(2,4,6)] <- max(final.data.trans$fert.total_logRS)

pred_tab2$metric <- "rich"
pred_tab2$test <- "Total Fertiliser"
pred_tab2$cat <- "Use Intensity"
pred_tab2$realm <- "Non-tropical"



# 3: Homogeneity across land uses (richess)
pred_tab3 <- do.call("rbind", replicate(6, pred_tab, simplify = FALSE))
pred_tab3$Predominant_land_use[1:4] <- c("Primary vegetation", "Primary vegetation", "Secondary vegetation", "Secondary vegetation")
pred_tab3$homogenRS[c(1,3,5)] <- min(final.data.trans$homogenRS)
pred_tab3$homogenRS[c(2,4,6)] <- max(final.data.trans$homogenRS)

pred_tab3$metric <- "rich"
pred_tab3$test <- "Homogeneity"
pred_tab3$cat <- "Land Use"
pred_tab3$realm <- "Non-tropical"


# 4: Fertiliser across use intensities (abundance)
pred_tab4 <- do.call("rbind", replicate(6, pred_tab, simplify = FALSE))
pred_tab4$Use_intensity[1:4] <- c("Minimal use", "Minimal use", "Light use", "Light use")
pred_tab4$fert.total_logRS[c(1,3,5)] <- min(final.data.trans$fert.total_logRS)
pred_tab4$fert.total_logRS[c(2,4,6)] <- max(final.data.trans$fert.total_logRS)

pred_tab4$metric <- "abun"
pred_tab4$test <- "Total Fertiliser"
pred_tab4$cat <- "Use Intensity"
pred_tab4$realm <- "Non-tropical"



# bring together table for abundance  projections

pred_tab_abtemp<- pred_tab4
pred_tab_srtemp <- rbind(pred_tab1, pred_tab2, pred_tab3)



# predict the result
result_srtrop <- PredictGLMERRandIter(model = srmod_trop$model, data = pred_tab_srtrop, nIters = 10000)
result_abtrop <- PredictGLMERRandIter(model = abmod_trop$model, data = pred_tab_abtrop, nIters = 10000)


# run temperate models
result_srtemp <- PredictGLMERRandIter(model = srmod_temp$model, data = pred_tab_srtemp, nIters = 10000)
result_abtemp <- PredictGLMERRandIter(model = abmod_temp$model, data = pred_tab_abtemp, nIters = 10000)



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
#dif_tab2 <- do.call(rbind.data.frame, dif_tab2)


colnames(dif_tab2) <- c("median", "lower", "upper")

dif_tab2$median <- as.numeric(dif_tab2$median)
dif_tab2$lower <- as.numeric(dif_tab2$lower)
dif_tab2$upper <- as.numeric(dif_tab2$upper)


# add details of tests
tests <- c(pred_tab_srtrop$test, pred_tab_srtemp$test, pred_tab_abtrop$test, pred_tab_abtemp$test)
tests <- tests[seq(1, length(tests), 2)]
dif_tab2$test <- tests

cat <- c(pred_tab_srtrop$cat, pred_tab_srtemp$cat, pred_tab_abtrop$cat, pred_tab_abtemp$cat)
cat <- cat[seq(1, length(cat), 2)]
dif_tab2$cat <- cat

realm <- c(pred_tab_srtrop$realm, pred_tab_srtemp$realm, pred_tab_abtrop$realm, pred_tab_abtemp$realm)
realm <- realm[seq(1, length(realm), 2)]
dif_tab2$realm <- realm


dif_tab2$colour <- NA
dif_tab2$colour[dif_tab2$cat == "Use Intensity"] <- rep(c("Minimal use", " Light use", "Intense use"))
dif_tab2$colour[dif_tab2$cat == "Land Use"] <- rep(c("Primary vegetation", "Secondary vegetation", "Cropland"))

metric <- c(pred_tab_srtrop$metric, pred_tab_srtemp$metric, pred_tab_abtrop$metric, pred_tab_abtemp$metric)
metric <- metric[seq(1, length(metric), 2)]
dif_tab2$metric <- metric

# save the results table
write.csv(dif_tab2, paste0(outdir, "/DifferencesTable_Interactions.csv"), row.names = F)



##%######################################################%##
#                                                          #
####        Create figure to present differences        ####
#                                                          #
##%######################################################%##


#dif_tab2 <- read.csv(paste0(outdir, "/DifferencesTable_Interactions.csv"))

dif_tab2$test <- sub("Percentage Natural Habitat", "Percentage Natural\n Habitat", dif_tab2$test)
dif_tab2$test <- sub("Number of Landcovers", "Number of\n Landcovers", dif_tab2$test)

dif_tab2$test <- factor(dif_tab2$test, levels = rev(c("Distance to forest",
                                                    "Percentage Natural\n Habitat", "Total Fertiliser",
                                                    "Number of\n Landcovers", "Homogeneity")))

dif_tab2$realm <- factor(dif_tab2$realm, levels = c("Non-tropical", "Tropical"))

dif_tab2$colour <- factor(dif_tab2$colour, levels = rev(c("Minimal use", " Light use", "Intense use", "Primary vegetation", "Secondary vegetation", "Cropland")))

dif_tab2$metric <- sub("abun", "Total Abundance", dif_tab2$metric)
dif_tab2$metric <- sub("rich", "Species Richness", dif_tab2$metric)

dif_tab2$metric <- factor(dif_tab2$metric, levels = c("Total Abundance", "Species Richness"))

# dif_tab$median[dif_tab$median == 0] <- NA
# dif_tab$lower[dif_tab$lower == 0] <- NA
# dif_tab$upper[dif_tab$upper == 0] <- NA


dif_tab2 <- droplevels(dif_tab2)


ggplot(data = dif_tab2[dif_tab2$cat == "Land Use", ]) + 
  geom_point(aes(x = test, y = median,  col = colour, shape = colour, alpha = realm), position = position_dodge(width = 0.9), size = 1.5) +
  geom_errorbar(aes(x = test, y = median, ymin = lower, ymax = upper, col = colour, alpha = realm), 
                position = position_dodge2(padding = 0.5),
                size = 0.2,
                width = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  #facet_grid(realm~  metric) + coord_flip() +
  #facet_grid(~ realm ) 
  coord_flip() +
    ylab("Percentage Change") + 
  scale_color_manual(values = rev(c("#006400", "#8B0000", "#EEAD0E"))) +
  theme_bw() +
  theme_custom + 
  theme(legend.position = "bottom", strip.background = element_rect(fill = NA),legend.box="vertical")


ggsave(filename = paste0(outdir, "/Difference_interactions_LU_plotV2.pdf"), height = 5, width = 6, units = "in")  

#"#66CD00", "#FFB90F", "#EE0000" minimal, light, intense

#"#006400", "#8B0000", "#EEAD0E" primary, secondary, cropland


p1 <- ggplot(data = dif_tab2[dif_tab2$cat == "Use Intensity" & dif_tab2$realm == "Tropical", ]) + 
  geom_point(aes(x = test, y = median, col = colour, shape = colour), position = position_dodge(width = 0.9), size = 1.5) +
  geom_errorbar(aes(x = test, y = median, ymin = lower, ymax = upper, col = colour), 
                position = position_dodge2(padding = 0.5),
                size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  facet_grid(~  metric) + coord_flip() +
  ylab("Percentage Change") + 
  scale_color_manual(values = rev(c("#66CD00", "#FFB90F", "#EE0000"))) +
  ylim(c(-80, 220)) +
  theme_bw() +
  theme_custom + 
  theme(legend.position = "none", strip.background = element_rect(fill = NA))

p2 <- ggplot(data = dif_tab2[dif_tab2$cat == "Use Intensity" & dif_tab2$realm == "Non-tropical", ]) + 
  geom_point(aes(x = test, y = median, col = colour, shape = colour), position = position_dodge(width = 0.9), size = 1.5) +
  geom_errorbar(aes(x = test, y = median, ymin = lower, ymax = upper, col = colour), 
                position = position_dodge2(padding = 0.5),
                size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  facet_grid(~  metric) + coord_flip() +
  ylab("Percentage Change") + 
  scale_color_manual(values = rev(c("#66CD00", "#FFB90F", "#EE0000"))) +
  ylim(c(-80, 220)) +
  theme_bw() +
  theme_custom + 
  theme(legend.position = "bottom", strip.background = element_rect(fill = NA))


library(cowplot)


legend <- get_legend(p2)

cowplot::plot_grid(p1, p2 +  theme(legend.position = "none"), legend,  nrow = 3, rel_heights = c(1,1,0.2))


ggsave(filename = paste0(outdir, "/Difference_interactions_UI_plotV2.pdf"), height = 7, width = 6, units = "in")  
