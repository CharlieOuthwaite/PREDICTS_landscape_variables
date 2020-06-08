##%######################################################%##
#                                                          #
####               Plotting model results               ####
#                                                          #
##%######################################################%##

# This script uses the model outputs from script 2 to plot the effects.

rm(list = ls())

# load libraries
library(StatisticalModels)

# source in edited plotting function
source("PlotGLMERContinuous_edited.r")

# directories
datadir <- "2_MODEL_SELECTION" 
outdir <-  "3_PLOTTING"

# load in the final model outputs
load(paste0(datadir, "/SRMOD_output.rdata"))
load(paste0(datadir, "/ABMOD_output.rdata"))
load(paste0(datadir, "/RCAR_output.rdata"))

# read in the values used for rescaling
scalers <- read.csv("1_PREDICTS_PLUS_VARIABLES/Scaling_values.csv")


##%######################################################%##
#                                                          #
####               Species richness plots               ####
#                                                          #
##%######################################################%##


pdf(file = paste0(outdir, "/SR_all_plots.pdf"), paper = "a4", width = 7, height = 10,
    title = "Species richness model: all effects")

par(mfrow = c(4,2))


### factor level plots ###

# Land use
PlotGLMERFactor(model = srmod$model,data = srmod$data,
                responseVar = "Species richness",seMultiplier = 1.96,
                logLink = "e",catEffects = c("Predominant_land_use"), params = list(las = 2, cex = 0.6, mar = c(6,3,2,4), adj = 0.5), xtext.srt = 90)



# Use intensity
PlotGLMERFactor(model = srmod$model,data = srmod$data,
                responseVar = "Species richness",seMultiplier = 1.96,
                logLink = "e",catEffects = c( "Use_intensity"), params = list(las = 2, cex = 0.6, mar = c(6,3,2,4), adj = 0.5), xtext.srt = 90)


# forest biomes
PlotGLMERFactor(model = srmod$model,data = srmod$data,
                responseVar = "Species richness",seMultiplier = 1.96,
                logLink = "e",catEffects = "Forest_biome", params = list(las = 2, cex = 0.6, mar = c(9,5,1,4), adj = 1), xtext.srt = 50, ylim = c(-70, 100))



# Tropical
PlotGLMERFactor(model = srmod$model,data = srmod$data,
                responseVar = "Species richness",seMultiplier = 1.96,
                logLink = "e",catEffects = c("Tropical"), params = list(las = 2, cex = 0.6, mar = c(4,4,4,2), adj = 0.5), xtext.srt = 90)



### categorical ###

# distance to forest
PlotGLMERContinuous_edit(model = srmod$model,
                         data = srmod$data, 
                         effects = "Hansen_mindist_log", 
                         otherContEffects = c("fert.total_log", "percNH", "landcovers.5k", "homogen"),
                         otherFactors = list(Forest_biome = "Tropical & Subtropical Moist Broadleaf Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use", Tropical = "Temperate"),
                         logLink = "e",
                         xlab ="Distance to forest (km)", 
                         ylab = "Species Richness", 
                         plotRug = TRUE, 
                         line.cols = c("#458B00"),
                         seMultiplier = 1.96,
                         transformX = FALSE, 
                         rescaled = c(scalers[1, 2], scalers[1, 3]), 
                         logged = TRUE,
                         xlim = c(0, 350),
                         ylim = c(5, 30),
                         params = list(cex = 0.6, mgp=c(3,1,0), mar = c(6,6,3,1), adj = 0.5))


# percNH

PlotGLMERContinuous_edit(model = srmod$model,
                         data = srmod$data, 
                         effects = "percNH", 
                         otherContEffects = c("fert.total_log", "Hansen_mindist_log", "landcovers.5k", "homogen"),
                         otherFactors = list(Forest_biome = "Tropical & Subtropical Moist Broadleaf Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use", Tropical = "Temperate"),
                         logLink = "e",
                         xlab ="Percentage of Natural Habitat", 
                         ylab = "Species Richness", 
                         plotRug = TRUE, 
                         line.cols =  c("#66CD00"),
                         seMultiplier = 1.96,
                         transformX = FALSE, 
                         rescaled = c(scalers[5, 2], scalers[5, 3]), 
                         logged = FALSE,
                         xlim = c(0, 100),
                         ylim = c(5, 30),
                         params = list(cex = 0.6, mgp=c(3,1,0), mar = c(6,6,3,1), adj = 0.5))


# fertiliser

PlotGLMERContinuous_edit(model = srmod$model, data = srmod$data, 
                         effects = "fert.total_log", 
                         otherContEffects = c("Hansen_mindist_log", "percNH", "landcovers.5k", "homogen"),
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use", Tropical = "Temperate"),
                         logLink = "e",
                         xlab ="Total fertiliser application (kgs)", 
                         ylab = "Species Richness", 
                         plotRug = TRUE, 
                         line.cols = c("#8B4500"),
                         seMultiplier = 1.96,
                         rescaled = c(scalers[2, 2], scalers[2, 3]),
                         logged = TRUE,
                         xlim = c(0, 2500000),
                         ylim = c(10, 30),
                         params = list(cex = 0.6, mgp=c(3,1,0), mar = c(6,6,3,1), adj = 0.5))



# landcovers.5k

PlotGLMERContinuous_edit(model = srmod$model, data = srmod$data, 
                         effects = "landcovers.5k", 
                         otherContEffects = c("percNH", "homogen", "fert.total_log", "Hansen_mindist_log"),
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Use_intensity = "Minimal use", Predominant_land_use = "Primary vegetation", Tropical = "Temperate"),
                         logLink = "e",
                         xlab ="Number of landcovers", 
                         ylab = "Species Richness", 
                         plotRug = FALSE, 
                         line.cols = c("#8B4513"),
                         seMultiplier = 1,
                         rescaled = c(scalers[3,2], scalers[3,3]),
                         logged = F,
                         params = list(cex = 0.6, mgp=c(3,1,0), mar = c(6,6,3,1), adj = 0.5),
                         xlim = c(0, 10),
                         ylim = c(10, 30))



# homogen


PlotGLMERContinuous_edit(model = srmod$model, data = srmod$data, 
                         effects = "homogen", 
                         otherContEffects = c("landcovers.5k", "Hansen_mindist_log", "percNH", "fert.total_log"),
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Use_intensity = "Minimal use", Predominant_land_use = "Primary vegetation", Tropical = "Temperate"),
                         logLink = "e",
                         xlab ="Homogeneity", 
                         ylab = "Species Richness", 
                         plotRug = TRUE, 
                         line.cols = c("#1C86EE"),
                         seMultiplier = 1.96, 
                         rescaled = c(scalers[4, 2], scalers[4, 3]),
                         logged = F,
                         xlim = c(0, 1),
                         ylim = c(10, 30),
                         params = list(cex = 0.6, mgp=c(3,1,0), mar = c(6,6,3,1), adj = 0.5))



### selected interactions ###

# land use : Use intensity


  ### need to predict this - write a function?


PlotGLMERFactor(model = srmod$model,data = srmod$data,
                responseVar = "Species richness",seMultiplier = 1.96,
                logLink = "e",catEffects = c("Predominant_land_use", "Use_intensity"), params = list(las = 2, cex = 0.6, mar = c(6,3,2,4), adj = 1), xtext.srt = 90)
title("USING THIS TO FILL GAP", adj = 0, line = 0.2)



# Land Use:Fertiliser

PlotGLMERContinuous_edit(model = srmod$model, data = srmod$data, 
                         effects = "fert.total_log", 
                         otherContEffects = c("Hansen_mindist_log", "percNH", "landcovers.5k", "homogen"),
                         byFactor = "Predominant_land_use",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use", Tropical = "Temperate"),
                         logLink = "e",
                         xlab ="Total fertiliser application (kgs)", 
                         ylab = "Species Richness", 
                         plotRug = TRUE, 
                         line.cols = c("#458B00", "#8B0000", "#FFA500"),
                         seMultiplier = 1.96,
                         rescaled = c(scalers[2, 2], scalers[2, 3]),
                         logged = TRUE,
                         xlim = c(0, 3000),
                         ylim = c(8, 30),
                         params = list(cex = 0.6, mgp=c(3,1,0), mar = c(6,6,3,1), adj = 0.5))
legend(1500, 30,c("Primary vegetation","Secondary vegetation","Cropland"),
       col=c("#458B00", "#8B0000", "#FFA500"),bty="n",lty=1)


# Land use : homogen

PlotGLMERContinuous_edit(model = srmod$model, data = srmod$data, 
                         effects = "homogen", 
                         otherContEffects = c("Hansen_mindist_log", "percNH", "landcovers.5k", "fert.total_log"),
                         byFactor = "Predominant_land_use",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use", Tropical = "Temperate"),
                         logLink = "e",
                         xlab ="Homogeneity", 
                         ylab = "Species Richness", 
                         plotRug = TRUE, 
                         line.cols = c("#458B00", "#8B0000", "#FFA500"),
                         seMultiplier = 1.96,
                         rescaled = c(scalers[4, 2], scalers[4, 3]),
                         logged = FALSE,
                         xlim = c(0, 1),
                         ylim = c(8, 30),
                         params = list(cex = 0.6, mgp=c(3,1,0), mar = c(6,6,3,1), adj = 0.5))
legend(0.6, 30,c("Primary vegetation","Secondary vegetation","Cropland"),
       col=c("#458B00", "#8B0000", "#FFA500"),bty="n",lty=1)




# land use : landcovers

PlotGLMERContinuous_edit(model = srmod$model, data = srmod$data, 
                         effects = "landcovers.5k", 
                         otherContEffects = c("Hansen_mindist_log", "percNH", "homogen", "fert.total_log"),
                         byFactor = "Predominant_land_use",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use", Tropical = "Temperate"),
                         logLink = "e",
                         xlab ="Number of landcovers", 
                         ylab = "Species Richness", 
                         plotRug = FALSE, 
                         line.cols = c("#458B00", "#8B0000", "#FFA500"),
                         seMultiplier = 1.96,
                         rescaled = c(scalers[3, 2], scalers[3, 3]),
                         logged = FALSE,
                         xlim = c(0, 10),
                         ylim = c(8, 30),
                         params = list(cex = 0.6, mgp=c(3,1,0), mar = c(6,6,3,1), adj = 0.5))
legend(6, 30,c("Primary vegetation","Secondary vegetation","Cropland"),
       col=c("#458B00", "#8B0000", "#FFA500"),bty="n",lty=1)






# Use intensity:fertiliser

PlotGLMERContinuous_edit(model = srmod$model, data = srmod$data, 
                         effects = "fert.total_log", 
                         otherContEffects = c("Hansen_mindist_log", "percNH", "landcovers.5k", "homogen"),
                         byFactor = "Use_intensity",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use", Tropical = "Temperate"),
                         logLink = "e",
                         xlab ="Total fertiliser application (kgs)", 
                         ylab = "Species Richness", 
                         plotRug = TRUE, 
                         line.cols = c("#66CD00", "#FFB90F", "#EE0000"),
                         seMultiplier = 1.96,
                         rescaled = c(scalers[2, 2], scalers[2, 3]),
                         logged = TRUE,
                         xlim = c(0, 3000),
                         ylim = c(8, 30),
                         params = list(cex = 0.6, mgp=c(3,1,0), mar = c(6,6,3,1), adj = 0.5))
legend(1500, 30,c("Minimal use","Light use","Intense use"),
       col=c("#66CD00", "#FFB90F", "#EE0000"),bty="n",lty=1)


# Use intensity : percNH

PlotGLMERContinuous_edit(model = srmod$model, data = srmod$data, 
                         effects = "percNH", 
                         otherContEffects = c("Hansen_mindist_log", "landcovers.5k", "fert.total_log", "homogen"),
                         byFactor = "Use_intensity",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use", Tropical = "Temperate"),
                         logLink = "e",
                         xlab ="Percentage of Natural Habitat", 
                         ylab = "Species Richness", 
                         plotRug = TRUE, 
                         line.cols = c("#66CD00", "#FFB90F", "#EE0000"),
                         seMultiplier = 1.96,
                         rescaled = c(scalers[5, 2], scalers[5, 3]),
                         logged = FALSE,
                         xlim = c(0, 100),
                         ylim = c(8, 30),
                         params = list(cex = 0.6, mgp=c(3,1,0), mar = c(6,6,3,1), adj = 0.5))
legend(0, 30,c("Minimal use","Light use","Intense use"),
       col=c("#66CD00", "#FFB90F", "#EE0000"),bty="n",lty=1)



# Tropical : Distance to forest

PlotGLMERContinuous_edit(model = srmod$model,
                         data = srmod$data, 
                         effects = "Hansen_mindist_log",
                         byFactor = "Tropical",
                         otherContEffects = c("fert.total_log", "percNH", "landcovers.5k", "homogen"),
                         otherFactors = list(Forest_biome = "Tropical & Subtropical Moist Broadleaf Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use", Tropical = "Temperate"),
                         logLink = "e",
                         xlab ="Distance to forest (km)", 
                         ylab = "Species Richness", 
                         plotRug = TRUE, 
                         line.cols = c("#009ACD", "#9932CC"),
                         seMultiplier = 1.96,
                         transformX = FALSE, 
                         rescaled = c(scalers[1, 2], scalers[1, 3]), 
                         logged = TRUE,
                         xlim = c(0, 20),
                         ylim = c(5, 30),
                         params = list(cex = 0.6, mgp=c(3,1,0), mar = c(6,6,3,1), adj = 0.5))
legend(12, 30,c("Temperate","Tropical"),
       col=c("#009ACD", "#9932CC"),bty="n",lty=1)




dev.off()






##%######################################################%##
#                                                          #
####                  Abundance plots                   ####
#                                                          #
##%######################################################%##



pdf(file = paste0(outdir, "/AB_all_plots.pdf"), paper = "a4", width = 7, height = 10,
    title = "Abundance model: all effects")

par(mfrow = c(4,2))


### factor level plots ###

# Land use
PlotGLMERFactor(model = abmod$model,data = abmod$data,
                responseVar = "Total Abundance",seMultiplier = 1.96,
                logLink = "e",catEffects = c("Predominant_land_use"), params = list(las = 2, cex = 0.6, mar = c(6,3,2,4), adj = 0.5), xtext.srt = 90)



# Use intensity
PlotGLMERFactor(model = abmod$model,data = abmod$data,
                responseVar = "Total Abundance",seMultiplier = 1.96,
                logLink = "e",catEffects = c( "Use_intensity"), params = list(las = 2, cex = 0.6, mar = c(6,3,2,4), adj = 0.5), xtext.srt = 90)


# forest biomes
PlotGLMERFactor(model = abmod$model,data = abmod$data,
                responseVar = "Total Abundance",seMultiplier = 1.96,
                logLink = "e",catEffects = "Forest_biome", params = list(las = 2, cex = 0.6, mar = c(9,5,1,4), adj = 1), xtext.srt = 50, ylim = c(-200, 650))



# Tropical
PlotGLMERFactor(model = abmod$model,data = abmod$data,
                responseVar = "Total Abundance",seMultiplier = 1.96,
                logLink = "e",catEffects = c("Tropical"), params = list(las = 2, cex = 0.6, mar = c(4,4,4,2), adj = 0.5), xtext.srt = 90)



### categorical ###

# distance to forest
PlotGLMERContinuous_edit(model = abmod$model,
                         data = abmod$data, 
                         effects = "Hansen_mindist_log", 
                         otherContEffects = c("fert.total_log", "percNH", "landcovers.5k", "homogen"),
                         otherFactors = list(Forest_biome = "Tropical & Subtropical Moist Broadleaf Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use", Tropical = "Temperate"),
                         logLink = "e",
                         xlab ="Distance to forest (km)", 
                         ylab = "Total Abundance", 
                         plotRug = TRUE, 
                         line.cols = c("#458B00"),
                         seMultiplier = 1.96,
                         transformX = FALSE, 
                         rescaled = c(scalers[1, 2], scalers[1, 3]), 
                         logged = TRUE,
                         xlim = c(0, 350),
                         ylim = c(50, 900),
                         params = list(cex = 0.6, mgp=c(3,1,0), mar = c(6,6,3,1), adj = 0.5))


# percNH

PlotGLMERContinuous_edit(model = abmod$model,
                         data = abmod$data, 
                         effects = "percNH", 
                         otherContEffects = c("fert.total_log", "Hansen_mindist_log", "landcovers.5k", "homogen"),
                         otherFactors = list(Forest_biome = "Tropical & Subtropical Moist Broadleaf Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use", Tropical = "Temperate"),
                         logLink = "e",
                         xlab ="Percentage of Natural Habitat", 
                         ylab = "Total Abundance", 
                         plotRug = TRUE, 
                         line.cols =  c("#66CD00"),
                         seMultiplier = 1.96,
                         transformX = FALSE, 
                         rescaled = c(scalers[5, 2], scalers[5, 3]), 
                         logged = FALSE,
                         xlim = c(0, 100),
                         ylim = c(50, 900),
                         params = list(cex = 0.6, mgp=c(3,1,0), mar = c(6,6,3,1), adj = 0.5))


# fertiliser

PlotGLMERContinuous_edit(model = abmod$model, data = abmod$data, 
                         effects = "fert.total_log", 
                         otherContEffects = c("Hansen_mindist_log", "percNH", "landcovers.5k", "homogen"),
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use", Tropical = "Temperate"),
                         logLink = "e",
                         xlab ="Total fertiliser application (kgs)", 
                         ylab = "Total Abundance", 
                         plotRug = TRUE, 
                         line.cols = c("#8B4500"),
                         seMultiplier = 1.96,
                         rescaled = c(scalers[2, 2], scalers[2, 3]),
                         logged = TRUE,
                         xlim = c(0, 2500000),
                         ylim = c(50, 300),
                         params = list(cex = 0.6, mgp=c(3,1,0), mar = c(6,6,3,1), adj = 0.5))



# landcovers.5k

PlotGLMERContinuous_edit(model = abmod$model, data = abmod$data, 
                         effects = "landcovers.5k", 
                         otherContEffects = c("percNH", "homogen", "fert.total_log", "Hansen_mindist_log"),
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Use_intensity = "Minimal use", Predominant_land_use = "Primary vegetation", Tropical = "Temperate"),
                         logLink = "e",
                         xlab ="Number of landcovers", 
                         ylab = "Total Abundance", 
                         plotRug = FALSE, 
                         line.cols = c("#8B4513"),
                         seMultiplier = 1,
                         rescaled = c(scalers[3,2], scalers[3,3]),
                         logged = F,
                         params = list(cex = 0.6, mgp=c(3,1,0), mar = c(6,6,3,1), adj = 0.5),
                         xlim = c(0, 10),
                         ylim = c(50, 300))



# homogen


PlotGLMERContinuous_edit(model = abmod$model, data = abmod$data, 
                         effects = "homogen", 
                         otherContEffects = c("landcovers.5k", "Hansen_mindist_log", "percNH", "fert.total_log"),
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Use_intensity = "Minimal use", Predominant_land_use = "Primary vegetation", Tropical = "Temperate"),
                         logLink = "e",
                         xlab ="Homogeneity", 
                         ylab = "Total Abundance", 
                         plotRug = TRUE, 
                         line.cols = c("#1C86EE"),
                         seMultiplier = 1.96, 
                         rescaled = c(scalers[4, 2], scalers[4, 3]),
                         logged = F,
                         xlim = c(0, 1),
                         ylim = c(50, 300),
                         params = list(cex = 0.6, mgp=c(3,1,0), mar = c(6,6,3,1), adj = 0.5))



### selected interactions ###

# land use : Use intensity


### need to predict this - write a function?


PlotGLMERFactor(model = abmod$model,data = abmod$data,
                responseVar = "Total Abundance",seMultiplier = 1.96,
                logLink = "e",catEffects = c("Predominant_land_use", "Use_intensity"), params = list(las = 2, cex = 0.6, mar = c(6,3,2,4), adj = 1), xtext.srt = 90)
title("USING THIS TO FILL GAP", adj = 0, line = 0.2)



# Land Use:Fertiliser

PlotGLMERContinuous_edit(model = abmod$model, data = abmod$data, 
                         effects = "fert.total_log", 
                         otherContEffects = c("Hansen_mindist_log", "percNH", "landcovers.5k", "homogen"),
                         byFactor = "Predominant_land_use",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use", Tropical = "Temperate"),
                         logLink = "e",
                         xlab ="Total fertiliser application (kgs)", 
                         ylab = "Total Abundance", 
                         plotRug = TRUE, 
                         line.cols = c("#458B00", "#8B0000", "#FFA500"),
                         seMultiplier = 1.96,
                         rescaled = c(scalers[2, 2], scalers[2, 3]),
                         logged = TRUE,
                         xlim = c(0, 3000),
                         ylim = c(50, 300),
                         params = list(cex = 0.6, mgp=c(3,1,0), mar = c(6,6,3,1), adj = 0.5))
legend(2000, 300,c("Primary vegetation","Secondary vegetation","Cropland"),
       col=c("#458B00", "#8B0000", "#FFA500"),bty="n",lty=1)




# Land Use:percNH

PlotGLMERContinuous_edit(model = abmod$model, data = abmod$data, 
                         effects = "percNH", 
                         otherContEffects = c("Hansen_mindist_log", "fert.total_log", "landcovers.5k", "homogen"),
                         byFactor = "Predominant_land_use",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use", Tropical = "Temperate"),
                         logLink = "e",
                         xlab ="Percentage of Natural Habitat", 
                         ylab = "Total Abundance", 
                         plotRug = TRUE, 
                         line.cols = c("#458B00", "#8B0000", "#FFA500"),
                         seMultiplier = 1.96,
                         rescaled = c(scalers[5, 2], scalers[5, 3]),
                         logged = FALSE,
                         xlim = c(0, 100),
                         ylim = c(50, 300),
                         params = list(cex = 0.6, mgp=c(3,1,0), mar = c(6,6,3,1), adj = 0.5))
legend(70, 300,c("Primary vegetation","Secondary vegetation","Cropland"),
       col=c("#458B00", "#8B0000", "#FFA500"),bty="n",lty=1)



# Use intensity:fertiliser

PlotGLMERContinuous_edit(model = abmod$model, data = abmod$data, 
                         effects = "fert.total_log", 
                         otherContEffects = c("Hansen_mindist_log", "percNH", "landcovers.5k", "homogen"),
                         byFactor = "Use_intensity",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use", Tropical = "Temperate"),
                         logLink = "e",
                         xlab ="Total fertiliser application (kgs)", 
                         ylab = "Total Abundance", 
                         plotRug = TRUE, 
                         line.cols = c("#66CD00", "#FFB90F", "#EE0000"),
                         seMultiplier = 1.96,
                         rescaled = c(scalers[2, 2], scalers[2, 3]),
                         logged = TRUE,
                         xlim = c(0, 300000),
                         ylim = c(50, 300),
                         params = list(cex = 0.6, mgp=c(3,1,0), mar = c(6,6,3,1), adj = 0.5))
legend(200000, 300,c("Minimal use","Light use","Intense use"),
       col=c("#66CD00", "#FFB90F", "#EE0000"),bty="n",lty=1)


# Use intensity : percNH

PlotGLMERContinuous_edit(model = abmod$model, data = abmod$data, 
                         effects = "percNH", 
                         otherContEffects = c("Hansen_mindist_log", "landcovers.5k", "fert.total_log", "homogen"),
                         byFactor = "Use_intensity",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use", Tropical = "Temperate"),
                         logLink = "e",
                         xlab ="Percentage of Natural Habitat", 
                         ylab = "Total Abundance", 
                         plotRug = TRUE, 
                         line.cols = c("#66CD00", "#FFB90F", "#EE0000"),
                         seMultiplier = 1.96,
                         rescaled = c(scalers[5, 2], scalers[5, 3]),
                         logged = FALSE,
                         xlim = c(0, 100),
                         ylim = c(50, 300),
                         params = list(cex = 0.6, mgp=c(3,1,0), mar = c(6,6,3,1), adj = 0.5))
legend(0, 300,c("Minimal use","Light use","Intense use"),
       col=c("#66CD00", "#FFB90F", "#EE0000"),bty="n",lty=1)



# Tropical : Distance to forest

PlotGLMERContinuous_edit(model = abmod$model,
                         data = abmod$data, 
                         effects = "Hansen_mindist_log",
                         byFactor = "Tropical",
                         otherContEffects = c("fert.total_log", "percNH", "landcovers.5k", "homogen"),
                         otherFactors = list(Forest_biome = "Tropical & Subtropical Moist Broadleaf Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use", Tropical = "Temperate"),
                         logLink = "e",
                         xlab ="Distance to forest (km)", 
                         ylab = "Total Abundance", 
                         plotRug = TRUE, 
                         line.cols = c("#009ACD", "#9932CC"),
                         seMultiplier = 1.96,
                         transformX = FALSE, 
                         rescaled = c(scalers[1, 2], scalers[1, 3]), 
                         logged = TRUE,
                         xlim = c(0, 20),
                         ylim = c(50, 700),
                         params = list(cex = 0.6, mgp=c(3,1,0), mar = c(6,6,3,1), adj = 0.5))
legend(15, 700,c("Temperate","Tropical"),
       col=c("#009ACD", "#9932CC"),bty="n",lty=1)




dev.off()




##%######################################################%##
#                                                          #
####                  RCAR model plots                  ####
#                                                          #
##%######################################################%##




pdf(file = paste0(outdir, "/RCAR_all_plots.pdf"), paper = "a4", width = 7, height = 10,
    title = "RCAR model: all effects")

par(mfrow = c(4,2))


### factor level plots ###

# Land use
PlotGLMERFactor(model = rcarmod$model,data = rcarmod$data,
                responseVar = "RCAR (square kilometres)",seMultiplier = 1.96,
                logLink = "10",catEffects = c("Predominant_land_use"), params = list(las = 2, cex = 0.6, mar = c(6,3,2,4), adj = 0.5), xtext.srt = 90)



# Use intensity
PlotGLMERFactor(model = rcarmod$model,data = rcarmod$data,
                responseVar = "RCAR (square kilometres)",seMultiplier = 1.96,
                logLink = "10",catEffects = c( "Use_intensity"), params = list(las = 2, cex = 0.6, mar = c(6,3,2,4), adj = 0.5), xtext.srt = 90)


# forest biomes
PlotGLMERFactor(model = rcarmod$model,data = rcarmod$data,
                responseVar = "RCAR (square kilometres)",seMultiplier = 1.96,
                logLink = "10",catEffects = "Forest_biome", params = list(las = 2, cex = 0.6, mar = c(9,5,1,4), adj = 1), xtext.srt = 50, ylim = c(-100, 150))



### categorical ###

# distance to forest
PlotGLMERContinuous_edit(model = rcarmod$model,
                         data = rcarmod$data, 
                         effects = "Hansen_mindist_log", 
                         otherContEffects = c("percNH", "landcovers.5k", "homogen"),
                         otherFactors = list(Forest_biome = "Tropical & Subtropical Moist Broadleaf Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use"),
                         logLink = "10",
                         xlab ="Distance to forest (km)", 
                         ylab = "RCAR (square kilometres)", 
                         plotRug = TRUE, 
                         line.cols = c("#458B00"),
                         seMultiplier = 1.96,
                         transformX = FALSE, 
                         rescaled = c(scalers[1, 2], scalers[1, 3]), 
                         logged = TRUE,
                         xlim = c(0, 250),
                         ylim = c(50000, 1500000),
                         params = list(cex = 0.6, mgp=c(4,1,0), mar = c(6,6,3,1), adj = 0.5))


# percNH

PlotGLMERContinuous_edit(model = rcarmod$model,
                         data = rcarmod$data, 
                         effects = "percNH", 
                         otherContEffects = c("Hansen_mindist_log", "landcovers.5k", "homogen"),
                         otherFactors = list(Forest_biome = "Tropical & Subtropical Moist Broadleaf Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use"),
                         logLink = "10",
                         xlab ="Percentage of Natural Habitat", 
                         ylab = "RCAR (square kilometres)", 
                         plotRug = TRUE, 
                         line.cols =  c("#66CD00"),
                         seMultiplier = 1.96,
                         transformX = FALSE, 
                         rescaled = c(scalers[5, 2], scalers[5, 3]), 
                         logged = FALSE,
                         xlim = c(0, 100),
                         ylim = c(50000, 1500000),
                         params = list(cex = 0.6, mgp=c(4,1,0), mar = c(6,6,3,1), adj = 0.5))



# landcovers.5k

PlotGLMERContinuous_edit(model = rcarmod$model, data = rcarmod$data, 
                         effects = "landcovers.5k", 
                         otherContEffects = c("percNH", "homogen", "Hansen_mindist_log"),
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Use_intensity = "Minimal use", Predominant_land_use = "Primary vegetation"),
                         logLink = "10",
                         xlab ="Number of landcovers", 
                         ylab = "RCAR (square kilometres)", 
                         plotRug = FALSE, 
                         line.cols = c("#8B4513"),
                         seMultiplier = 1,
                         rescaled = c(scalers[3,2], scalers[3,3]),
                         logged = F,
                         params = list(cex = 0.6, mgp=c(4,1,0), mar = c(6,6,3,1), adj = 0.5),
                         xlim = c(0, 10),
                         ylim = c(50000, 1500000))
                         


# homogen


PlotGLMERContinuous_edit(model = rcarmod$model, data = rcarmod$data, 
                         effects = "homogen", 
                         otherContEffects = c("landcovers.5k", "Hansen_mindist_log", "percNH"),
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Use_intensity = "Minimal use", Predominant_land_use = "Primary vegetation"),
                         logLink = "10",
                         xlab ="Homogeneity", 
                         ylab = "RCAR (square kilometres)", 
                         plotRug = TRUE, 
                         line.cols = c("#1C86EE"),
                         seMultiplier = 1.96, 
                         rescaled = c(scalers[4, 2], scalers[4, 3]),
                         logged = F,
                         xlim = c(0, 1),
                         ylim = c(50000, 2000000),
                         params = list(cex = 0.6, mgp=c(4,1,0), mar = c(6,6,3,1), adj = 0.5))



### selected interactions ###

# land use : Use intensity


### need to predict this - write a function?


PlotGLMERFactor(model = rcarmod$model,data = rcarmod$data,
                responseVar = "RCAR (square kilometres)",seMultiplier = 1.96,
                logLink = "10",catEffects = c("Predominant_land_use", "Use_intensity"), params = list(las = 2, cex = 0.6, mar = c(6,3,2,4), adj = 0.5), xtext.srt = 90)
title("USING THIS TO FILL GAP", adj = 0, line = 0.2)




# Land Use:percNH

PlotGLMERContinuous_edit(model = rcarmod$model, data = rcarmod$data, 
                         effects = "percNH", 
                         otherContEffects = c("Hansen_mindist_log", "landcovers.5k", "homogen"),
                         byFactor = "Predominant_land_use",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use"),
                         logLink = "10",
                         xlab ="Percentage of Natural Habitat", 
                         ylab = "RCAR (square kilometres)", 
                         plotRug = TRUE, 
                         line.cols = c("#458B00", "#8B0000", "#FFA500"),
                         seMultiplier = 1.96,
                         rescaled = c(scalers[5, 2], scalers[5, 3]),
                         logged = FALSE,
                         xlim = c(0, 100),
                         ylim = c(50000, 2500000),
                         params = list(cex = 0.6, mgp=c(4,1,0), mar = c(6,6,3,1), adj = 0.5))
legend(70, 2500000,c("Primary vegetation","Secondary vegetation","Cropland"),
       col=c("#458B00", "#8B0000", "#FFA500"),bty="n",lty=1)


# Land Use:homogen

PlotGLMERContinuous_edit(model = rcarmod$model, data = rcarmod$data, 
                         effects = "homogen", 
                         otherContEffects = c("Hansen_mindist_log", "landcovers.5k", "percNH"),
                         byFactor = "Predominant_land_use",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use"),
                         logLink = "10",
                         xlab ="Homogeneity", 
                         ylab = "RCAR (square kilometres)", 
                         plotRug = TRUE, 
                         line.cols = c("#458B00", "#8B0000", "#FFA500"),
                         seMultiplier = 1.96,
                         rescaled = c(scalers[4, 2], scalers[4, 3]),
                         logged = FALSE,
                         xlim = c(0, 1),
                         ylim = c(50000, 2500000),
                         params = list(cex = 0.6, mgp=c(4,1,0), mar = c(6,6,3,1), adj = 0.5))
legend(0.6, 2500000,c("Primary vegetation","Secondary vegetation","Cropland"),
       col=c("#458B00", "#8B0000", "#FFA500"),bty="n",lty=1)



# Use intensity : Hansen_mindist_log

PlotGLMERContinuous_edit(model = rcarmod$model, data = rcarmod$data, 
                         effects = "Hansen_mindist_log", 
                         otherContEffects = c("percNH", "landcovers.5k", "homogen"),
                         byFactor = "Use_intensity",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use"),
                         logLink = "10",
                         xlab ="Distance to forest (Km)", 
                         ylab = "RCAR (square kilometres)", 
                         plotRug = TRUE, 
                         line.cols = c("#66CD00", "#FFB90F", "#EE0000"),
                         seMultiplier = 1.96,
                         rescaled = c(scalers[1, 2], scalers[1, 3]),
                         logged = TRUE,
                         xlim = c(0, 20),
                         ylim = c(50000, 2500000),
                         params = list(cex = 0.6, mgp=c(4,1,0), mar = c(6,6,3,1), adj = 0.5))
legend(15, 2500000,c("Minimal use","Light use","Intense use"),
       col=c("#66CD00", "#FFB90F", "#EE0000"),bty="n",lty=1)



# Use intensity : landcovers

PlotGLMERContinuous_edit(model = rcarmod$model, data = rcarmod$data, 
                         effects = "landcovers.5k", 
                         otherContEffects = c("percNH", "Hansen_mindist_log", "homogen"),
                         byFactor = "Use_intensity",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use"),
                         logLink = "10",
                         xlab ="Number of Landcovers", 
                         ylab = "RCAR (square kilometres)", 
                         plotRug = FALSE, 
                         line.cols = c("#66CD00", "#FFB90F", "#EE0000"),
                         seMultiplier = 1.96,
                         rescaled = c(scalers[3, 2], scalers[3, 3]),
                         logged = FALSE,
                         xlim = c(0, 10),
                         ylim = c(50000, 2500000),
                         params = list(cex = 0.6, mgp=c(4,1,0), mar = c(6,6,3,1), adj = 0.5))
legend(6, 2500000,c("Minimal use","Light use","Intense use"),
       col=c("#66CD00", "#FFB90F", "#EE0000"),bty="n",lty=1)


# Use intensity : homogen

PlotGLMERContinuous_edit(model = rcarmod$model, data = rcarmod$data, 
                         effects = "homogen", 
                         otherContEffects = c("percNH", "Hansen_mindist_log", "landcovers.5k"),
                         byFactor = "Use_intensity",
                         otherFactors = list(Forest_biome = "Temperate Broadleaf & Mixed Forests", Predominant_land_use = "Primary vegetation", Use_intensity = "Minimal use"),
                         logLink = "10",
                         xlab ="Homogeneity", 
                         ylab = "RCAR (square kilometres)", 
                         plotRug = TRUE, 
                         line.cols = c("#66CD00", "#FFB90F", "#EE0000"),
                         seMultiplier = 1.96,
                         rescaled = c(scalers[4, 2], scalers[4, 3]),
                         logged = FALSE,
                         xlim = c(0, 1),
                         ylim = c(50000, 2500000),
                         params = list(cex = 0.6, mgp=c(4,1,0), mar = c(6,6,3,1), adj = 0.5))
legend(0.6, 2500000,c("Minimal use","Light use","Intense use"),
       col=c("#66CD00", "#FFB90F", "#EE0000"),bty="n",lty=1)



dev.off()







##%######################################################%##
#                                                          #
####           Figure combinations for paper            ####
#                                                          #
##%######################################################%##


# Interaction between predominant land use and use intensity

