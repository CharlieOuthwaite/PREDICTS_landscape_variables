##%######################################################%##
#                                                          #
####          Chao - estimated richness model           ####
#                                                          #
##%######################################################%##


rm(list = ls())

# load libraries
library(devtools)
#install_github(repo = "timnewbold/predicts-demo",subdir = "predictsFunctions")
library(predictsFunctions)
library(StatisticalModels)
library(ggplot2)
library(dplyr)
#library(yarg)


# set the data folder
datadir <- "0_DATA"

# where to save the final dataset
outdir <- "5_Additional_Tests"

# read in the complete PREDICTS dataset
pred.data <- readRDS(paste0(datadir, "/database.rds")) # 3250404 rows

# correct sampling effort 
predicts <- CorrectSamplingEffort(pred.data)

# merge sites: this combines potential subsamples within one site
predicts <- MergeSites(predicts) # 2906994 rows

# # assign sites to major groups based on Phylum/Class
# predicts$group <- NA
# 
# predicts[predicts$Class == "Aves", "group"] <- "Birds"
# predicts[predicts$Class == "Mammalia", "group"] <- "Mammals"
# predicts[predicts$Class == "Amphibia", "group"] <- "Amphibians"
# predicts[predicts$Class == "Reptilia", "group"] <- "Reptiles"
# 
# predicts[predicts$Kingdom == "Plantae", "group"] <- "Plants"
# 
# predicts[predicts$Kingdom == "Fungi", "group"] <- "Fungi_SlimeMoulds"
# 
# predicts[predicts$Kingdom == "Protozoa", "group"] <- "Protozoa"
# 
# predicts[predicts$Phylum == "Arthropoda", "group"] <- "Inverts"
# predicts[predicts$Phylum == "Annelida", "group"] <- "Inverts"
# predicts[predicts$Phylum == "Mollusca", "group"] <- "Inverts"
# predicts[predicts$Phylum == "Nematoda", "group"] <- "Inverts"
# predicts[predicts$Phylum == "Platyhelminthes", "group"] <- "Inverts"
# predicts[predicts$Phylum == "Onychophora", "group"] <- "Inverts"
# 
# 
# predicts[predicts$Phylum == "Mycetozoa", "group"] <- "Fungi_SlimeMoulds"
# predicts[predicts$Phylum == "Ascomycota", "group"] <- "Fungi_SlimeMoulds"
# predicts[predicts$Phylum == "Basidiomycota", "group"] <- "Fungi_SlimeMoulds"
# predicts[predicts$Phylum == "Glomeromycota", "group"] <- "Fungi_SlimeMoulds"
# 
# table(predicts[is.na(predicts$group), "Phylum"])
# table(predicts[is.na(predicts$group), "Class"])
# View(predicts[is.na(predicts$group), ]) # 144 rows
# 
# # check whether sites sample more than one group 
# group_count <- as.data.frame.matrix(table(predicts$SSBS, predicts$group))
# group_count[group_count >= 1] <- 1
# group_count$SSBS <- rownames(group_count)
# 
# 
# group_count <- cbind(group_count$SSBS, group_count[, 1:7] %>% mutate(sum = rowSums(.)))
# summary(group_count$sum)
# 
# # some sites look at more than one group, combine these into a "Multiple" group
# multiple <- group_count[group_count$sum > 1, 1]
# 
# predicts[predicts$SSBS %in% multiple, 'group'] <- "Multiple"
# table(predicts$group)


# Calculate site level metrics, package yarg required to get Chao richness
pred.sites.metrics <- yarg::SiteMetrics(predicts, 
                                  extra.cols = c("Predominant_land_use", "SSB", "SSBS", "Diversity_metric_is_suitable_for_Chao"),
                                  srEstimators = "Chao") # 22678 rows


# round the estimated species richness values to integers.
pred.sites.metrics$ChaoR <- round(pred.sites.metrics$ChaoR,0)

# site level data cropland only
sites.sub <- pred.sites.metrics[!pred.sites.metrics$Predominant_land_use %in% c("Urban", "Pasture", "Cannot decide", "Plantation forest"), ]

# remove sites with NA in lat/long columns
sites.sub <- sites.sub[!is.na(sites.sub$Longitude),  ] # 15612 rows

# remove sites where taxon information not available
#sites.sub <- sites.sub[!is.na(sites.sub$group),  ] # 15612 rows

# load in original data that includes all the landscape variables
load(file = "1_PREDICTS_PLUS_VARIABLES/PREDICTS_dataset_inc_variables_TRANS.rdata")

# filter by site
sites.sub <- sites.sub[sites.sub$SSBS %in% final.data.trans$SSBS, ]

sites.sub <- sites.sub[ , c("SSBS", "ChaoR")]

# merge the datasets to add the ChaoR info
final.data.trans <- merge(final.data.trans, sites.sub, by = "SSBS") # 10393 rows

# remove sites that don't have a ChaoR estimate
final.data.trans <- final.data.trans[!is.na(final.data.trans$ChaoR), ] # 5849 5019

final.data.trans$group <- factor(final.data.trans$group)


# dataset summaries

table(final.data.trans$Predominant_land_use)
# Primary vegetation Secondary vegetation             Cropland 
#               2130                 2014                 1705

# Primary vegetation Secondary vegetation             Cropland 
#               1670                 1750                 1599 

table(final.data.trans$Tropical)
# Temperate  Tropical 
#      3706      2143
#      3311      1708

table(final.data.trans$group)
# Amphibians             Birds Fungi_SlimeMoulds           Inverts           Mammals          Multiple 
#         56              1376                 1              3267               163               314 
# Plants          Reptiles 
#    654                18 

# Amphibians             Birds Fungi_SlimeMoulds           Inverts           Mammals          Multiple            Plants 
#         38              1075                48              3019               163               197               473 
# Reptiles 
#        6

# Split the dataset based on realm

final.data.trans_trop <- final.data.trans[final.data.trans$Tropical == "Tropical", ]
nrow(final.data.trans_trop)# 2143 1708
final.data.trans_trop <- droplevels(final.data.trans_trop)

final.data.trans_temp <- final.data.trans[final.data.trans$Tropical == "Temperate", ]
nrow(final.data.trans_temp) # 3706 3311
final.data.trans_temp <- droplevels(final.data.trans_temp)

#### run the model selection process ####
### 1. tropical model
system.time({ChaoR_trop <- GLMERSelect(modelData = final.data.trans_trop, 
                                    responseVar = "ChaoR",
                                    fitFamily = "poisson", 
                                    fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity"),
                                    fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1, percNH = 1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)", 
                                    fixedInteractions = c("Predominant_land_use:poly(fert.total_log,1)", 
                                                          "Predominant_land_use:poly(Hansen_mindist_log,1)",
                                                          "Predominant_land_use:poly(landcovers.5k,1)",
                                                          "Predominant_land_use:poly(homogen,1)",
                                                          "Predominant_land_use:poly(percNH,1)",
                                                          "Use_intensity:poly(fert.total_log,1)",
                                                          "Use_intensity:poly(Hansen_mindist_log,1)",
                                                          "Use_intensity:poly(landcovers.5k,1)",
                                                          "Use_intensity:poly(homogen,1)",
                                                          "Use_intensity:poly(percNH,1)",
                                                          "Predominant_land_use:Use_intensity"), verbose = T)}) 
# save the output
save(ChaoR_trop, file = paste0(outdir, "/ChaoR_Tropical_Model_selection.rdata"))
# load(paste0(outdir, "/ChaoR_Tropical_Model_selection.rdata"))

# take a look at the model output
summary(ChaoR_trop$model)

# selected model:
# ChaoR ~ Predominant_land_use + Use_intensity + 
# poly(fert.total_log, 1) + poly(homogen, 1) + 
# Predominant_land_use:Use_intensity +      
# (1 | SS) + (1 | SSB) + (1 | SSBS) + (1 | group)


table(final.data.trans_trop$Predominant_land_use)
# Primary vegetation Secondary vegetation             Cropland 
#               1033                  853                  257
table(final.data.trans_trop$Use_intensity)
# Minimal use   Light use Intense use 
#        1283         692         168 



### 2. temperate model

system.time({ChaoR_temp <- GLMERSelect(modelData = final.data.trans_temp, 
                                    responseVar = "ChaoR",
                                    fitFamily = "poisson", 
                                    fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity"),
                                    fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1, percNH = 1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)", 
                                    fixedInteractions = c("Predominant_land_use:poly(fert.total_log,1)", 
                                                          "Predominant_land_use:poly(Hansen_mindist_log,1)",
                                                          "Predominant_land_use:poly(landcovers.5k,1)",
                                                          "Predominant_land_use:poly(homogen,1)",
                                                          "Predominant_land_use:poly(percNH,1)",
                                                          "Use_intensity:poly(fert.total_log,1)",
                                                          "Use_intensity:poly(Hansen_mindist_log,1)",
                                                          "Use_intensity:poly(landcovers.5k,1)",
                                                          "Use_intensity:poly(homogen,1)",
                                                          "Use_intensity:poly(percNH,1)",
                                                          "Predominant_land_use:Use_intensity"), verbose = T)}) 
# save the output
save(ChaoR_temp, file = paste0(outdir, "/ChaoR_Temperate_Model_selection.rdata"))
#load(paste0(outdir, "/ChaoR_Temperate_Model_selection.rdata"))

# take a look at the model output
summary(ChaoR_temp$model)

# Selected model:
# ChaoR ~ Predominant_land_use + Use_intensity + 
# poly(Hansen_mindist_log, 1) + poly(landcovers.5k, 1) + poly(homogen, 1) + poly(percNH, 1) +
# Predominant_land_use:poly(Hansen_mindist_log, 1) + Predominant_land_use:poly(landcovers.5k, 1) + Predominant_land_use:poly(homogen,1) + 
# Use_intensity:poly(percNH, 1) + Predominant_land_use:Use_intensity +  
# (1 | SS) + (1 | SSB) + (1 | SSBS) + (1 | group)


table(final.data.trans_temp$Predominant_land_use)
# Primary vegetation Secondary vegetation             Cropland 
#               1097                 1161                 1448

table(final.data.trans_temp$Use_intensity)
# Minimal use   Light use Intense use 
#        1740        1007         959 



# compare the final models to those selected for richness

#load(paste0(outdir, "/SPECIESRICHNESS_Tropical_Model_selection.rdata"))
#load(paste0(outdir, "/SPECIESRICHNESS_Temperate_Model_selection.rdata"))

load(paste0("2_MODEL_SELECTION/SRMOD_Tropical_output.rdata"))
load(paste0("2_MODEL_SELECTION/SRMOD_Temperate_output.rdata"))


summary(srmod_trop$model)
summary(srmod_temp$model)


## for when issues with the function, load previously used dataset
Trop_data <- ChaoRmod_trop$data
Temp_data <- ChaoR_temp$data

##### Run the same models as was selected for species richness ####


## 1. tropical model

mod_struc1 <- "Predominant_land_use + Forest_biome + Use_intensity + poly(fert.total_log, 1) + poly(homogen, 1) + poly(percNH, 1) + poly(landcovers.5k, 1) + Predominant_land_use:poly(fert.total_log, 1) + Predominant_land_use:poly(homogen, 1) + Predominant_land_use:poly(percNH, 1) + Use_intensity:poly(landcovers.5k, 1) + Use_intensity:poly(percNH, 1)"


ChaoRmod_trop <- GLMER(modelData = Trop_data, responseVar = "ChaoR", fitFamily = "poisson",
                    fixedStruct = mod_struc1,
                    randomStruct = "(1|SS) + (1|SSB) + (1|SSBS)", REML = TRUE, maxIters = 20000)

# no warnings or errors

summary(ChaoRmod_trop$model)

# extract the coefficents of the model
coefs_trop <- fixef(ChaoRmod_trop$model)

# save the coefficients
write.csv(coefs_trop, file = paste0(outdir, "/ChaoR_Tropical_coefs_nogroup.csv"), row.names = F)

# save the model output
save(ChaoRmod_trop, file = paste0(outdir, "/ChaoR_REML_Tropical_output_nogroup.rdata"))
#load(file = paste0(outdir, "/ChaoR_REML_Tropical_output.rdata"))


## 2. temperate model


mod_struc2 <- "Predominant_land_use + Forest_biome + Use_intensity + poly(fert.total_log, 1) + poly(Hansen_mindist_log, 1) + poly(homogen, 1) + poly(percNH, 1) + Predominant_land_use:poly(homogen, 1) + Use_intensity:poly(fert.total_log, 1) + Use_intensity:poly(percNH,1) + Predominant_land_use:Use_intensity"

ChaoRmod_temp <- GLMER(modelData = Temp_data, responseVar = "ChaoR", fitFamily = "poisson",
                    fixedStruct = mod_struc2,
                    randomStruct = "(1|SS) + (1|SSB) + (1|SSBS)", REML = TRUE, maxIters = 30000)

# no warnings or errors

summary(ChaoRmod_temp$model)

# extract the coefficents of the model
coefs_temp <- fixef(ChaoRmod_temp$model)

# save the coefficients
write.csv(coefs_temp, file = paste0(outdir, "/ChaoR_Temperate_coefs_nogroup.csv"), row.names = F)

# save the model output
save(ChaoRmod_temp, file = paste0(outdir, "/ChaoR_REML_Temperate_output_nogroup.rdata"))
# load(file = paste0(outdir, "/ChaoR_REML_Temperate_output.rdata"))





##### create tables of coefficients for SR and ChaoR models #####

library(sjPlot)


tab_model(srmod_trop$model, ChaoRmod_trop$model, transform = NULL, file = paste0(outdir, "/Table_SR_ChaoR_TROP_nogroup.html"), 
          show.icc = F, show.obs = T, show.ngroups = F, show.ci = F)

tab_model(srmod_temp$model, ChaoRmod_temp$model, transform = NULL, file = paste0(outdir, "/Table_SR_ChaoR_TEMP_nogroup.html"), 
          show.icc = F, show.obs = T, show.ngroups = F, show.ci = F)



##### Look at spread of data across land uses #####

table(final.data.trans_trop$Predominant_land_use)
# Primary vegetation Secondary vegetation             Cropland 
#               1922                 1416                  381 

table(final.data.trans_temp$Predominant_land_use)
#  Primary vegetation Secondary vegetation             Cropland 
#                2728                 2390                 1556 

table(ChaoRmod_trop$data$Predominant_land_use)
# Primary vegetation Secondary vegetation             Cropland 
#               1033                  853                  257

table(ChaoRmod_temp$data$Predominant_land_use)
# Primary vegetation Secondary vegetation             Cropland 
#               1097                 1161                 1448