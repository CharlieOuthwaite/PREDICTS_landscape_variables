##%######################################################%##
#                                                          #
####                 2. Model selection                 ####
#                                                          #
##%######################################################%##

# 2. Richness models

# This script carries out the model selection process using the dataset generated
# in the previous script. Models are run for 3 biodiversity metrics. 

rm(list = ls())


# load libraries
library(devtools)
#install_github(repo = "timnewbold/StatisticalModels")
library(StatisticalModels)
library(roquefort)
library(cowplot)
library(gridGraphics)
library(sjPlot)


# directories
datadir <- "1_PREDICTS_PLUS_VARIABLES"
outdir <- "2_MODEL_SELECTION"


# read in the PREDICTS datasets with landscape variables
load(paste0(datadir, "/PREDICTS_dataset_inc_variables_TRANS.rdata"))
final.data.trans <- droplevels(final.data.trans)


# Split the dataset based on realm

final.data.trans_trop <- final.data.trans[final.data.trans$Tropical == "Tropical", ]
nrow(final.data.trans_trop) # 3719
final.data.trans_trop <- droplevels(final.data.trans_trop)

final.data.trans_temp <- final.data.trans[final.data.trans$Tropical == "Temperate", ]
nrow(final.data.trans_temp) # 6674
final.data.trans_temp <- droplevels(final.data.trans_temp)


##%######################################################%##
#                                                          #
####              Model selection process               ####
####         including all landscape variables          ####
#                                                          #
##%######################################################%##


#### 1. Tropical subset ####

# run the model selection process
mod_dat <- final.data.trans_trop[, c(1:2, 5:8, 11:12, 30:34)]

system.time({sr_trop <- GLMERSelect(modelData = mod_dat, 
                                responseVar = "Species_richness",
                                fitFamily = "poisson", 
                                fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity"),
                                fixedTerms = list(fert.total_logRS = 1, Hansen_mindist_logRS = 1, landcovers.5kRS = 1, homogenRS = 1, percNHRS = 1),
                                randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)", 
                                fixedInteractions = c("Predominant_land_use:poly(fert.total_logRS,1)", 
                                                      "Predominant_land_use:poly(Hansen_mindist_logRS,1)",
                                                      "Predominant_land_use:poly(landcovers.5kRS,1)", 
                                                      "Predominant_land_use:poly(homogenRS,1)", 
                                                      "Predominant_land_use:poly(percNHRS,1)",
                                                      "Use_intensity:poly(fert.total_logRS,1)", 
                                                      "Use_intensity:poly(Hansen_mindist_logRS,1)",
                                                      "Use_intensity:poly(landcovers.5kRS,1)",
                                                      "Use_intensity:poly(homogenRS,1)",
                                                      "Use_intensity:poly(percNHRS,1)",
                                                      "Predominant_land_use:Use_intensity"), verbose = F)}) 
# save the output
save(sr_trop, file = paste0(outdir, "/SPECIESRICHNESS_Tropical_Model_selection.rdata"))

# take a look at the model output
summary(sr_trop$model)

# selected model:
# Species_richness ~ Predominant_land_use + Forest_biome + Use_intensity +  
# poly(fert.total_log, 1) + poly(homogen, 1) + poly(percNH, 1) + poly(landcovers.5k, 1) +
# Predominant_land_use:poly(fert.total_log, 1) + Predominant_land_use:poly(homogen, 1) + Predominant_land_use:poly(percNH, 1) + 
# Use_intensity:poly(landcovers.5k, 1) + Use_intensity:poly(percNH, 1)  + 
# (1 | SS) +  (1 | SSB) + (1 | SSBS)

# extract the stats produced as part of the model selection process
sr_trop_stats <- as.data.frame(sr_trop$stats)

# save these
write.csv(sr_trop_stats, file = paste0(outdir, "/sr_trop_stats.csv"), row.names = F)




#### 2. Temperate subset ####

# run the model selection process

system.time({sr_temp <- GLMERSelect(modelData = final.data.trans_temp, 
                                responseVar = "Species_richness",
                                fitFamily = "poisson", 
                                fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity"),
                                fixedTerms = list(fert.total_logRS = 1, Hansen_mindist_logRS = 1, landcovers.5kRS = 1, homogenRS = 1, percNHRS = 1),
                                randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)", 
                                fixedInteractions = c("Predominant_land_use:poly(fert.total_logRS,1)", 
                                                      "Predominant_land_use:poly(Hansen_mindist_logRS,1)",
                                                      "Predominant_land_use:poly(landcovers.5kRS,1)", 
                                                      "Predominant_land_use:poly(homogenRS,1)", 
                                                      "Predominant_land_use:poly(percNHRS,1)",
                                                      "Use_intensity:poly(fert.total_logRS,1)", 
                                                      "Use_intensity:poly(Hansen_mindist_logRS,1)",
                                                      "Use_intensity:poly(landcovers.5kRS,1)",
                                                      "Use_intensity:poly(homogenRS,1)",
                                                      "Use_intensity:poly(percNHRS,1)",
                                                      "Predominant_land_use:Use_intensity"), verbose = F)}) 
# save the output
save(sr_temp, file = paste0(outdir, "/SPECIESRICHNESS_Temperate_Model_selection.rdata"))

# take a look at the model output
summary(sr_temp$model)

# selected model:
# Species_richness ~ Predominant_land_use + Forest_biome + Use_intensity +  
#  poly(fert.total_log, 1) + poly(Hansen_mindist_log, 1) + poly(homogen,1) + poly(percNH, 1) + 
# Predominant_land_use:poly(homogen, 1) + Use_intensity:poly(fert.total_log, 1) + Use_intensity:poly(percNH,1) + Predominant_land_use:Use_intensity +
# (1 | SS) + (1 | SSB) + (1 | SSBS) 

# extract the stats produced as part of the model selection process
sr_temp_stats <- as.data.frame(sr_temp$stats)

# save these
write.csv(sr_temp_stats, file = paste0(outdir, "/sr_temp_stats.csv"), row.names = F)


##%######################################################%##
#                                                          #
####       Rerun final selected models using REML       ####
#                                                          #
##%######################################################%##

# load(paste0(outdir, "/SPECIESRICHNESS_Tropical_Model_selection_nestedgrp.rdata"))
# load(paste0(outdir, "/SPECIESRICHNESS_Temperate_Model_selection_nestedgrp.rdata"))


# 1. Tropical
summary(sr_trop$model)

# selected model:
# Species_richness ~ Predominant_land_use + Forest_biome + Use_intensity +  
# poly(fert.total_log, 1) + poly(homogen, 1) + poly(percNH, 1) + poly(landcovers.5k, 1) +
# Predominant_land_use:poly(fert.total_log, 1) + Predominant_land_use:poly(homogen, 1) + Predominant_land_use:poly(percNH, 1) + 
# Use_intensity:poly(landcovers.5k, 1) + Use_intensity:poly(percNH, 1)  + 
# (1 | SS) +  (1 | SSB) + (1 | SSBS) + (1 | group:SS)

mod_struc <- "Predominant_land_use + Forest_biome + Use_intensity + poly(fert.total_logRS, 1) + poly(homogenRS, 1) + poly(percNHRS, 1) + poly(landcovers.5kRS, 1) + Predominant_land_use:poly(fert.total_logRS, 1) + Predominant_land_use:poly(homogenRS, 1) + Predominant_land_use:poly(percNHRS, 1) + Use_intensity:poly(landcovers.5kRS, 1) + Use_intensity:poly(percNHRS, 1)"


srmod_trop <- GLMER(modelData = final.data.trans_trop, responseVar = "Species_richness", fitFamily = "poisson",
                    fixedStruct = mod_struc,
                    randomStruct = "(1|SS) + (1|SSB) + (1|SSBS)", REML = TRUE,
                    maxIters = 40000)


# take a look at the model output
summary(srmod_trop$model)

# extract the coefficents of the model
coefs_trop <- fixef(srmod_trop$model)

# save the coefficients
write.csv(coefs_trop, file = paste0(outdir, "/SPECIESRICHNESS_Tropical_coefs.csv"), row.names = F)

# save the model output
save(srmod_trop, file = paste0(outdir, "/SRMOD_Tropical_output.rdata"))
#load(file = paste0(outdir, "/SRMOD_Tropical_output.rdata"))


# 2. Temperate
summary(sr_temp$model)

# selected model:
# Species_richness ~ Predominant_land_use + Forest_biome + Use_intensity +  
# poly(fert.total_log, 1) + poly(Hansen_mindist_log, 1) + poly(homogen,1) + poly(percNH, 1) + 
# Predominant_land_use:poly(homogen, 1) + Use_intensity:poly(fert.total_log, 1) + Use_intensity:poly(percNH,1) + Predominant_land_use:Use_intensity +
# (1 | SS) + (1 | SSB) + (1 | SSBS)

mod_struc <- "Predominant_land_use + Forest_biome + Use_intensity + poly(fert.total_logRS, 1) + poly(Hansen_mindist_logRS, 1) + poly(homogenRS, 1) + poly(percNHRS, 1) + Predominant_land_use:poly(homogenRS, 1) + Use_intensity:poly(fert.total_logRS, 1) + Use_intensity:poly(percNHRS,1) + Predominant_land_use:Use_intensity"

srmod_temp <- GLMER(modelData = final.data.trans_temp, responseVar = "Species_richness", fitFamily = "poisson",
                    fixedStruct = mod_struc,
                    randomStruct = "(1|SS) + (1|SSB) + (1|SSBS)", REML = TRUE,
                    maxIters = 100000)


summary(srmod_temp$model)

# extract the coefficents of the model
coefs_temp <- fixef(srmod_temp$model)

# save the coefficients
write.csv(coefs_temp, file = paste0(outdir, "/SPECIESRICHNESS_Temperate_coefs.csv"), row.names = F)

# save the model output
save(srmod_temp, file = paste0(outdir, "/SRMOD_Temperate_output.rdata"))
#load(file = paste0(outdir, "/SRMOD_Temperate_output.rdata"))


##%######################################################%##
#                                                          #
####                  R-squared values                  ####
#                                                          #
##%######################################################%##

sr_trop_r2 <- R2GLMER(srmod_trop$model) # [1] conditional, [2] marginal
sr_temp_r2 <- R2GLMER(srmod_temp$model)

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
write.csv(result, file = paste0(outdir, "/Rsquareds_SR.csv"))


# save model output tables

tab_model(srmod_trop$model, transform = NULL, file = paste0(outdir, "/SR_Trop_output_table.html"))
tab_model(srmod_temp$model, transform = NULL, file = paste0(outdir, "/SR_Temp_output_table.html"))


# combined table for paper

tab_model(srmod_trop$model, srmod_temp$model, transform = NULL, file = paste0(outdir, "/Table1_Richness.html"), 
          show.icc = F, show.obs = F, show.ngroups = F, show.ci = F)

##%######################################################%##
#                                                          #
####                    Model checks                    ####
#                                                          #
##%######################################################%##



### Tropical ###

## 1. Checking the fitted vs residuals relationship
p1 <- plot(srmod_trop$model)
# this kind of check is not informative for model with poisson error. Not included in supplementary.


## 2. Normality of Residuals

pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(srmod_trop$model), main = "")
qqline(resid(srmod_trop$model))
p2 <- recordPlot()
invisible(dev.off())



## 3. Check for spatial autocorrelation

sr_test<-roquefort::SpatialAutocorrelationTest(model=srmod_trop, all.data=final.data.trans)


summary(sr_test)

# percentage of studies that show spatial autocorrelation?
perc_auto <- (length(which(sr_test$P<0.05))/length(sr_test$P))*100

# 8.33%

sr_test_vals <- as.data.frame(sr_test$P)
sr_test_vals$`sr_test$P` <- round(sr_test_vals$`sr_test$P`, digits = 4)

label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")

p3 <- ggplot(data = sr_test_vals ) +
  geom_histogram(aes(x = sr_test_vals$`sr_test$P`)) +
  geom_vline(xintercept = 0.05, col = "red") +
  geom_text(aes(x = 0.2, y = 60, label = label1), size = 4, check_overlap = T) +
  theme_bw() +
  xlab("P-value") +
  ylab("Frequency") +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1)


# observed vs predicted
df <- cbind(exp(predict(srmod_trop$model)), final.data.trans_trop$Species_richness)
df <- as.data.frame(df)

R2 <- round(cor(df[, 1], df[, 2]), digits = 3)

r2lab <- paste("R^2 == ", R2, sep = "")


p4 <- ggplot(df, aes(x=df[, 1], y= df[,2])) +
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  annotate(geom = "text", label = r2lab, x = 55, y = 450, parse = T) +
  labs(x='Predicted Values', y='Observed Values')+
  theme_bw()

cowplot::plot_grid(p2,p3,p4,
                   labels = c("A.", "B.", "C."), nrow = 2)


ggsave(file = paste0(outdir, "/SR_Tropical_model_checks_plots.pdf"), height = 9, width = 9)



### Temperate ###

## 1. Checking the fitted vs residuals relationship
p1 <- plot(srmod_temp$model)
# this kind of check is not informative for model with poisson error. Not included in supplementary.


## 2. Normality of Residuals

pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(srmod_temp$model), main = "")
qqline(resid(srmod_temp$model))
p2 <- recordPlot()
invisible(dev.off())



## 3. Check for spatial autocorrelation

sr_test<-roquefort::SpatialAutocorrelationTest(model=srmod_temp, all.data=final.data.trans)


summary(sr_test)

# percentage of studies that show spatial autocorrelation?
perc_auto <- (length(which(sr_test$P<0.05))/length(sr_test$P))*100

# 7.02%

sr_test_vals <- as.data.frame(sr_test$P)
sr_test_vals$`sr_test$P` <- round(sr_test_vals$`sr_test$P`, digits = 4)

label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")

p3 <- ggplot(data = sr_test_vals ) +
  geom_histogram(aes(x = sr_test_vals$`sr_test$P`)) +
  geom_vline(xintercept = 0.05, col = "red") +
  geom_text(aes(x = 0.2, y = 60, label = label1), size = 4, check_overlap = T) +
  theme_bw() +
  xlab("P-value") +
  ylab("Frequency") +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1)



# observed vs predicted
df <- cbind(exp(predict(srmod_temp$model)), final.data.trans_temp$Species_richness)
df <- as.data.frame(df)

R2 <- round(cor(df[, 1], df[, 2]), digits = 3)

r2lab <- paste("R^2 == ", R2, sep = "")

p4 <- ggplot(df, aes(x=df[, 1], y= df[,2])) +
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  annotate(geom = "text", label = r2lab, x = 30, y = 255, parse = T) +
  labs(x='Predicted Values', y='Observed Values')+
  theme_bw()

cowplot::plot_grid(p2,p3,p4,
                   labels = c("A.", "B.", "C."), nrow = 2)



ggsave(file = paste0(outdir, "/SR_Temperate_model_checks_plots.pdf"), height = 9, width = 9)

