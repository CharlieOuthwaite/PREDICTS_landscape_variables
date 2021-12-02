##%######################################################%##
#                                                          #
####      5. Additional Tests: Taxon random effect      ####
#                                                          #
##%######################################################%##

# in this script, rerun the models with "group" as a random effect
# group was added in script 1 by looking at the taxonomic groups sampled 
# at each site. Where multiple were sampled, group = "Multiple". 

rm(list = ls())

# load libraries
library(devtools)
library(StatisticalModels)
library(roquefort)
library(cowplot)
library(gridGraphics)
library(sjPlot)
library(lme4)


# set directories
outdir <- "5_Additional_Tests"
dir.create(outdir)
datadir <- "1_PREDICTS_PLUS_VARIABLES"


#### organise data ####

# organise complete dataset with group info
load(file = paste0(datadir, "/sites_sub_inc_group.rdata"))
load(file = paste0(datadir, "/PREDICTS_dataset_inc_variables_TRANS.rdata"))

# get the group names from the sites.sub dataset, merge with organised dataset
sites.sub <- sites.sub[, c("SSBS", "group")]

final.data.trans <- merge(final.data.trans, sites.sub, by = "SSBS")

# look at spread of groups
table(final.data.trans$group)

# Amphibians      Birds      Fungi    Inverts    Mammals   Multiple     Plants   Reptiles 
#         73       1958        196       4149        903        686       2404         24

# set as a factor level
final.data.trans$group <- as.factor(final.data.trans$group)


#### run models with group level random effect ####

#### 1. Abundance models ####

# Split the dataset based on realm

final.data.trans_trop <- final.data.trans[final.data.trans$Tropical == "Tropical", ]
nrow(final.data.trans_trop) # 3719

final.data.trans_temp <- final.data.trans[final.data.trans$Tropical == "Temperate", ]
nrow(final.data.trans_temp) # 6674

# separate out the data where abundance column is not NA
final.data.trans_trop_ABUN <- final.data.trans_trop[!is.na(final.data.trans_trop$Total_abundance), ] # 3314 rows

# log the abundance values
final.data.trans_trop_ABUN$logAbun <- log(final.data.trans_trop_ABUN$Total_abundance+1)


# separate out the data where abundance column is not NA
final.data.trans_temp_ABUN <- final.data.trans_temp[!is.na(final.data.trans_temp$Total_abundance), ] # 5740 rows

# log the abundance values
final.data.trans_temp_ABUN$logAbun <- log(final.data.trans_temp_ABUN$Total_abundance+1)

final.data.trans_trop_ABUN <- droplevels(final.data.trans_trop_ABUN)
final.data.trans_temp_ABUN <- droplevels(final.data.trans_temp_ABUN)



#### 1. Tropical subset ####

# separate out the data where abundance column is not NA
mod_dat <- final.data.trans_trop_ABUN


# check for NAs in other columns
summary(is.na(mod_dat))

# run the model selection process
system.time({ab_trop <- GLMERSelect(modelData = mod_dat, 
                                    responseVar = "logAbun",
                                    fitFamily = "gaussian", 
                                    fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity"),
                                    fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1, percNH = 1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|group)", 
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
                                                          "Predominant_land_use:Use_intensity"), verbose = F)})

# take a look at the model output
summary(ab_trop$model)

# save the output
save(ab_trop, file = paste0(outdir, "/ABUNDANCE_Tropical_Model_Selection_group.rdata"))

# extract the stats produced as part of the model selection process
ab_trop_stats <- as.data.frame(ab_trop$stats)

# save these
write.csv(ab_trop_stats, file = paste0(outdir, "/ab_Trop_stats_group.csv"), row.names = F)




#### 2. Temperate subset ####

# separate out the data where abundance column is not NA
mod_dat <- final.data.trans_temp_ABUN

# run the model selection process
system.time({ab_temp <- GLMERSelect(modelData = mod_dat, 
                                    responseVar = "logAbun",
                                    fitFamily = "gaussian", 
                                    fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity"),
                                    fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1, percNH = 1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|group)", 
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
                                                          "Predominant_land_use:Use_intensity"), verbose = F)})

# take a look at the model output
summary(ab_temp$model)

# save the output
save(ab_temp, file = paste0(outdir, "/ABUNDANCE_Temperate_Model_Selection_group.rdata"))

# extract the stats produced as part of the model selection process
ab_temp_stats <- as.data.frame(ab_temp$stats)

# save these
write.csv(ab_temp_stats, file = paste0(outdir, "/ab_Temp_stats_group.csv"), row.names = F)



#### take a look at R2 values ####

abr2_trop <- R2GLMER(ab_trop$model) # [1] conditional, [2] marginal
abr2_temp <- R2GLMER(ab_temp$model)

# Look at the proportion of variance unexplained by the random effects, but
# explained by the fixed effects...  Marginal/(1 - (Conditional - Marginal))

abper_trop <- abr2_trop[[2]]/(1-(abr2_trop[[1]] - abr2_trop[[2]]))
abper_temp <- abr2_temp[[2]]/(1-(abr2_temp[[1]] - abr2_temp[[2]]))

# organise values in a table

result <- NULL

abres_trop <- unlist(c(abr2_trop[1], abr2_trop[2], abper_trop))
abres_temp  <- unlist(c(abr2_temp[1], abr2_temp[2], abper_temp))

result <- rbind(abres_trop, abres_temp)

colnames(result)[3] <- "proportion unexplained variance explained by fixed effects"

rownames(result) <- c("Abun Trop", "Abun Temp")

# save
write.csv(result, file = paste0(outdir, "/Rsquareds_Abun_group.csv"))


#### 2. Richness models ####


#### 1. Tropical subset ####

# run the model selection process

system.time({sr_trop <- GLMERSelect(modelData = final.data.trans_trop, 
                                    responseVar = "Species_richness",
                                    fitFamily = "poisson", 
                                    fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity"),
                                    fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1, percNH = 1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)+(1|group)", 
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
                                                          "Predominant_land_use:Use_intensity"), verbose = F)}) 
# save the output
save(sr_trop, file = paste0(outdir, "/SPECIESRICHNESS_Tropical_Model_selection_group.rdata"))

# take a look at the model output
summary(sr_trop$model)

# extract the stats produced as part of the model selection process
sr_trop_stats <- as.data.frame(sr_trop$stats)

# save these
write.csv(sr_trop_stats, file = paste0(outdir, "/sr_trop_stats_group.csv"), row.names = F)




#### 2. Temperate subset ####

# run the model selection process

system.time({sr_temp <- GLMERSelect(modelData = final.data.trans_temp, 
                                    responseVar = "Species_richness",
                                    fitFamily = "poisson", 
                                    fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity"),
                                    fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1, percNH = 1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)+(1|group)", 
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
                                                          "Predominant_land_use:Use_intensity"), verbose = F)}) 
# save the output
save(sr_temp, file = paste0(outdir, "/SPECIESRICHNESS_Temperate_Model_selection_group.rdata"))

# take a look at the model output
summary(sr_temp$model)

# extract the stats produced as part of the model selection process
sr_temp_stats <- as.data.frame(sr_temp$stats)

# save these
write.csv(sr_temp_stats, file = paste0(outdir, "/sr_temp_stats_group.csv"), row.names = F)


## R squareds ##


sr_trop_r2 <- R2GLMER(sr_trop$model) # [1] conditional, [2] marginal
sr_temp_r2 <- R2GLMER(sr_temp$model)

# Look at the proportion of variance unexplained by the random effects, but
# explained by the fixed effects...  Marginal/(1 - (Conditional - Marginal))


srper_trop <- sr_trop_r2[[2]]/(1-(sr_trop_r2[[1]] - sr_trop_r2[[2]]))
srper_temp <- sr_temp_r2[[2]]/(1-(sr_temp_r2[[1]] - sr_temp_r2[[2]]))


# organise values in a table

result <- NULL

srres_trop <- unlist(c(sr_trop_r2[1], sr_trop_r2[2], srper_trop))
srres_temp  <- unlist(c(sr_temp_r2[1], sr_temp_r2[2], srper_temp))

result <- rbind(srres_trop, srres_temp)

colnames(result)[3] <- "proportion unexplained variance explained by fixed effects"

rownames(result) <- c("SR Trop", "SR Temp")

# save
write.csv(result, file = paste0(outdir, "/Rsquareds_SR_groups.csv"))



# check AIC of models



load(file = paste0(outdir, "/ABUNDANCE_Tropical_Model_Selection_group.rdata"))
load(file = paste0(outdir, "/ABUNDANCE_Temperate_Model_Selection_group.rdata"))
load(file = paste0(outdir, "/SPECIESRICHNESS_Tropical_Model_selection_group.rdata"))
load(file = paste0(outdir, "/SPECIESRICHNESS_Temperate_Model_selection_group.rdata"))

ab_trop_group <- ab_trop
ab_temp_group <- ab_temp

sr_trop_group <- sr_trop
sr_temp_group <- sr_temp


# load in models without group random effect

load(file = "2_MODEL_SELECTION/ABUNDANCE_Tropical_Model_Selection.rdata")
load(file = "2_MODEL_SELECTION/ABUNDANCE_Temperate_Model_Selection.rdata")
load(file = "2_MODEL_SELECTION/SPECIESRICHNESS_Tropical_Model_selection.rdata")
load(file = "2_MODEL_SELECTION/SPECIESRICHNESS_Temperate_Model_selection.rdata")

summary(ab_trop_group$model)
# logAbun ~ Forest_biome + Use_intensity + poly(homogen, 1) + Predominant_land_use:poly(Hansen_mindist_log,  
#  1) + Use_intensity:poly(landcovers.5k, 1) + Predominant_land_use +  
#  poly(Hansen_mindist_log, 1) + poly(landcovers.5k, 1) + (1 |      SS) + (1 | SSB) + (1 | group)
summary(ab_trop$model)
# same as above 

summary(sr_trop_group$model)
# Species_richness ~ Predominant_land_use + Forest_biome + Use_intensity +  
# poly(fert.total_log, 1) + poly(homogen, 1) + Predominant_land_use:poly(fert.total_log,  
# 1) + Predominant_land_use:poly(homogen, 1) + Predominant_land_use:poly(percNH,
# 1) + Use_intensity:poly(landcovers.5k, 1) + Use_intensity:poly(percNH,  
#  1) + poly(percNH, 1) + poly(landcovers.5k, 1) + (1 | SS) +      (1 | SSB) + (1 | SSBS) + (1 | group)

summary(sr_trop$model)

# same as above 
AIC(ab_trop_group$model, ab_trop$model)
AIC(sr_trop_group$model, sr_trop$model)


AIC(ab_temp_group$model, ab_temp$model)
AIC(sr_temp_group$model, sr_temp$model)





