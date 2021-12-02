##%######################################################%##
#                                                          #
####  Alternative transformation of data including 0s   ####
#                                                          #
##%######################################################%##


# comment from reviewer:
# Perhaps conducting an inverse hyperbolic sine (or arcsinh) transformation would be more 
# appropriate/helpful in this case (https://doi.org/10.1111/obes.12325).

rm(list = ls())


# load libraries
library(devtools)
#install_github(repo = "timnewbold/StatisticalModels")
library(StatisticalModels)
library(roquefort)
library(cowplot)
library(gridGraphics)
library(sjPlot)
library(lme4)

# directories
datadir <- "1_PREDICTS_PLUS_VARIABLES"
outdir <- "2_MODEL_SELECTION"

# load in the data
load(paste0(datadir, "/PREDICTS_dataset_inc_variables_TRANS.rdata"))



# function for inverse hyperbolic sine (arcsinh)
ihs <- function(x) {
  y <- log(x + sqrt(x ^ 2 + 1))
  return(y)
}

# back transform using sinh()


# transform the fertiliser and distance data using this approach instead of log(x+1)
final.data.trans$fert.total_arc <- ihs(final.data.trans$fert.total)
final.data.trans$Hansen_mindist_arc <- ihs(final.data.trans$Hansen_mindist)

hist(final.data.trans$fert.total)
hist(final.data.trans$fert.total_arc)
hist(final.data.trans$fert.total_log)

hist(final.data.trans$Hansen_mindist)
hist(final.data.trans$Hansen_mindist_arc)
hist(final.data.trans$Hansen_mindist_log)

# rescale 
final.data.trans$fert.total_arcRS <- scale(final.data.trans$fert.total_arc)
final.data.trans$Hansen_mindist_arcRS <-scale(final.data.trans$Hansen_mindist_arc)

hist(final.data.trans$fert.total_logRS)
hist(final.data.trans$fert.total_arcRS)
hist(final.data.trans$Hansen_mindist_logRS)
hist(final.data.trans$Hansen_mindist_arcRS)

quantile(final.data.trans$fert.total_logRS)
quantile(final.data.trans$fert.total_arcRS)

quantile(final.data.trans$Hansen_mindist_logRS)
quantile(final.data.trans$Hansen_mindist_arcRS)


# split into realms and abundance datasets

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

## rerun models and compare coefficients


# 1. Tropical


mod_struc <- "Forest_biome + Use_intensity + Predominant_land_use + poly(homogenRS,1) + poly(Hansen_mindist_arcRS,1) + poly(landcovers.5kRS,1) + Predominant_land_use:poly(Hansen_mindist_arcRS,1) + Use_intensity:poly(landcovers.5kRS,1)"

mod_dat <- final.data.trans_trop_ABUN[, c(1:2, 5:8, 12, 27:39)]


abmod_trop_arc <- GLMER(modelData = mod_dat, responseVar = "logAbun", fitFamily = "gaussian",
                    fixedStruct = mod_struc,
                    randomStruct = "(1|SS) + (1|SSB)", REML = TRUE)

# take a look at the model output
summary(abmod_trop_arc$model)

# load in original model output
load("2_MODEL_SELECTION/ABMOD_Tropical_output.Rdata")
summary(abmod_trop$model)


# 2. Temperate


mod_struc <- "Predominant_land_use + Forest_biome + Use_intensity +  poly(fert.total_arcRS,1) + poly(landcovers.5kRS,1) + poly(homogenRS,1) + poly(percNHRS,1) + Use_intensity:poly(fert.total_arcRS,1) + Predominant_land_use:Use_intensity"

mod_dat <- final.data.trans_temp_ABUN[, c(1:2, 5:8, 12, 27:39)]

abmod_temp_arc <- GLMER(modelData = mod_dat, responseVar = "logAbun", fitFamily = "gaussian",
                    fixedStruct = mod_struc,
                    randomStruct = "(1|SS) + (1|SSB)", REML = TRUE)

# take a look at the model output
summary(abmod_temp_arc$model)

# load in original model output
load("2_MODEL_SELECTION/ABMOD_Temperate_output.Rdata")
summary(abmod_temp$model)


