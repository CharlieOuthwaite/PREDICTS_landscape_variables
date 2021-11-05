##%######################################################%##
#                                                          #
####                 2. Model selection                 ####
#                                                          #
##%######################################################%##

# 1. Abundance models

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
library(lme4)



# directories
datadir <- "1_PREDICTS_PLUS_VARIABLES"
outdir <- "2_MODEL_SELECTION"

# read in the PREDICTS datasets with landscape variables
load(paste0(datadir, "/PREDICTS_dataset_inc_variables_TRANS.rdata"))
final.data.trans <- droplevels(final.data.trans)


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

save(final.data.trans_trop_ABUN, file = paste0(outdir, "/final.data.trans_trop_ABUN.rdata"))
save(final.data.trans_temp_ABUN, file = paste0(outdir, "/final.data.trans_temp_ABUN.rdata"))

##%######################################################%##
#                                                          #
####              Model selection process               ####
####         including all landscape variables          ####
#                                                          #
##%######################################################%##

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

# selected model

# logAbun ~ Predominant_land_use + Forest_biome + Use_intensity + poly(homogen, 1) + poly(Hansen_mindist_log, 1) + poly(landcovers.5k, 1)
# Predominant_land_use:poly(Hansen_mindist_log, 1) + Use_intensity:poly(landcovers.5k, 1) +  
#  + (1 | SS) + (1 | SSB) + (1 | group)

# save the output
save(ab_trop, file = paste0(outdir, "/ABUNDANCE_Tropical_Model_Selection_incgroup.rdata"))

# extract the stats produced as part of the model selection process
ab_trop_stats <- as.data.frame(ab_trop$stats)

# save these
write.csv(ab_trop_stats, file = paste0(outdir, "/ab_Trop_stats.csv"), row.names = F)




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

# selected model:
# logAbun ~ Predominant_land_use + Forest_biome + Use_intensity +  
# poly(fert.total_log, 1) + poly(landcovers.5k, 1) + poly(homogen, 1) + poly(percNH, 1) + poly(Hansen_mindist_log, 1) +
# Predominant_land_use:poly(Hansen_mindist_log, 1) + Use_intensity:poly(fert.total_log, 1) + Predominant_land_use:Use_intensity +  
# (1 | SS) + (1 | SSB) + (1 | group)

# save the output
save(ab_temp, file = paste0(outdir, "/ABUNDANCE_Temperate_Model_Selection_incgroup.rdata"))

# extract the stats produced as part of the model selection process
ab_temp_stats <- as.data.frame(ab_temp$stats)

# save these
write.csv(ab_temp_stats, file = paste0(outdir, "/ab_Temp_stats.csv"), row.names = F)


# load(paste0(outdir, "/ABUNDANCE_Tropical_Model_Selection.rdata"))
# load(paste0(outdir, "/ABUNDANCE_Temperate_Model_Selection.rdata"))


##%######################################################%##
#                                                          #
####        Run final selected model using REML         ####
#                                                          #
##%######################################################%##

# Take a look at the final models selected for each realm and rerun using REML

# 1. Tropical
summary(ab_trop$model)

# logAbun ~ Predominant_land_use + Forest_biome + Use_intensity + 
# poly(homogen, 1) + poly(Hansen_mindist_log, 1) + poly(landcovers.5k, 1)
# Predominant_land_use:poly(Hansen_mindist_log, 1) + Use_intensity:poly(landcovers.5k, 1) +  
#  + (1 | SS) + (1 | SSB) + (1 | group)


mod_struc <- "Forest_biome + Use_intensity + Predominant_land_use + poly(homogen,1) + poly(Hansen_mindist_log,1) + poly(landcovers.5k,1) + Predominant_land_use:poly(Hansen_mindist_log,1) + Use_intensity:poly(landcovers.5k,1)"

abmod_trop <- GLMER(modelData = final.data.trans_trop_ABUN, responseVar = "logAbun", fitFamily = "gaussian",
               fixedStruct = mod_struc,
               randomStruct = "(1|SS) + (1|SSB) + (1 | group)", REML = TRUE)

# take a look at the model output
summary(abmod_trop$model)

# extract the coefficents of the model
coefs_trop <- fixef(abmod_trop$model)

# save the coefficients
write.csv(coefs_trop, file = paste0(outdir, "/ABUNDANCE_Trop_coefs.csv"), row.names = T)

# save the model output
save(abmod_trop, file = paste0(outdir, "/ABMOD_Tropical_output.rdata"))





# 2. Temperate
summary(ab_temp$model)


# logAbun ~ Predominant_land_use + Forest_biome + Use_intensity +  
# poly(fert.total_log, 1) + poly(landcovers.5k, 1) + poly(homogen, 1) + poly(percNH, 1) + poly(Hansen_mindist_log, 1) +
# Predominant_land_use:poly(Hansen_mindist_log, 1) + Use_intensity:poly(fert.total_log, 1) + Predominant_land_use:Use_intensity +  
# (1 | SS) + (1 | SSB) + (1 | group)


mod_struc <- "Predominant_land_use + Forest_biome + Use_intensity +  poly(fert.total_log,1) + poly(landcovers.5k,1) + poly(homogen,1) + poly(percNH,1) + poly(Hansen_mindist_log, 1)+ Predominant_land_use:poly(Hansen_mindist_log, 1) + Use_intensity:poly(fert.total_log,1) + Predominant_land_use:Use_intensity"

abmod_temp <- GLMER(modelData = final.data.trans_temp_ABUN, responseVar = "logAbun", fitFamily = "gaussian",
                    fixedStruct = mod_struc,
                    randomStruct = "(1|SS) + (1|SSB) + (1 | group)", REML = TRUE)

# take a look at the model output
summary(abmod_temp$model)

# extract the coefficents of the model
coefs_temp <- fixef(abmod_temp$model)

# save the coefficients
write.csv(coefs_temp, file = paste0(outdir, "/ABUNDANCE_Temperate_coefs.csv"), row.names = F)

# save the model output
save(abmod_temp, file = paste0(outdir, "/ABMOD_Temperate_output.rdata"))



##%######################################################%##
#                                                          #
####                  R-squared values                  ####
#                                                          #
##%######################################################%##

abr2_trop <- R2GLMER(abmod_trop$model) # [1] conditional, [2] marginal
abr2_temp <- R2GLMER(abmod_temp$model)

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
write.csv(result, file = paste0(outdir, "/Rsquareds_Abun.csv"))



# save model output tables

tab_model(abmod_trop$model, transform = NULL, file = paste0(outdir, "/AB_Trop_output_table.html"))
tab_model(abmod_temp$model, transform = NULL, file = paste0(outdir, "/AB_Temp_output_table.html"))




##%######################################################%##
#                                                          #
####                    Model checks                    ####
#                                                          #
##%######################################################%##



### Temperate ###


## 1. Checking the fitted vs residuals relationship
p1 <- plot(abmod_temp$model)


## 2. Normality of Residuals

pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(abmod_temp$model), main = "")
qqline(resid(abmod_temp$model))
p2 <- recordPlot()
invisible(dev.off())



## 3. Check for spatial autocorrelation

# needed for following function

abun_dat <- final.data.trans_temp[!is.na(final.data.trans_temp$Total_abundance), ]

abmod_temp$data$SSBS <- abun_dat$SSBS
abmod_temp$data$Latitude <- abun_dat$Latitude
abmod_temp$data$Longitude <- abun_dat$Longitude


ab_test<-roquefort::SpatialAutocorrelationTest(model=abmod_temp, all.data=abun_dat)


summary(ab_test)

# percentage of studies that show spatial autocorrelation?
perc_auto <- (length(which(ab_test$P<0.05))/length(ab_test$P))*100

# 8.6%

ab_test_vals <- as.data.frame(ab_test$P)
ab_test_vals$`ab_test$P` <- round(ab_test_vals$`ab_test$P`, digits = 4)

label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")

p3 <- ggplot(data = ab_test_vals ) +
  geom_histogram(aes(x = ab_test_vals$`ab_test$P`)) +
  geom_vline(xintercept = 0.05, col = "red") +
  geom_text(aes(x = 0.2, y = 50, label = label1), size = 4, check_overlap = T) +
  theme_bw() +
  xlab("P-value") +
  ylab("Frequency") +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1)



cowplot::plot_grid(p1,p2,p3,
          labels = c("A.", "B.", "C."))


ggsave(file = paste0(outdir, "/Abun_Temperate_model_checks_plots.pdf"), height = 9, width = 9)




### Tropical ###


## 1. Checking the fitted vs residuals relationship
p1 <- plot(abmod_trop$model)


## 2. Normality of Residuals

pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(abmod_trop$model), main = "")
qqline(resid(abmod_trop$model))
p2 <- recordPlot()
invisible(dev.off())


## 3. Check for spatial autocorrelation

# needed for following function

abun_dat <- final.data.trans_trop[!is.na(final.data.trans_trop$Total_abundance), ]

abmod_trop$data$SSBS <- abun_dat$SSBS
abmod_trop$data$Latitude <- abun_dat$Latitude
abmod_trop$data$Longitude <- abun_dat$Longitude


ab_test<-roquefort::SpatialAutocorrelationTest(model=abmod_trop, all.data=abun_dat)


summary(ab_test)

# percentage of studies that show spatial autocorrelation?
perc_auto <- (length(which(ab_test$P<0.05))/length(ab_test$P))*100

# 8.5%

ab_test_vals <- as.data.frame(ab_test$P)
ab_test_vals$`ab_test$P` <- round(ab_test_vals$`ab_test$P`, digits = 4)

label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")

p3 <- ggplot(data = ab_test_vals ) +
  geom_histogram(aes(x = ab_test_vals$`ab_test$P`)) +
  geom_vline(xintercept = 0.05, col = "red") +
  geom_text(aes(x = 0.2, y = 50, label = label1), size = 4, check_overlap = T) +
  theme_bw() +
  xlab("P-value") +
  ylab("Frequency") +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1)



cowplot::plot_grid(p1,p2,p3,
                   labels = c("A.", "B.", "C."))


ggsave(file = paste0(outdir, "/Abun_Tropical_model_checks_plots.pdf"), height = 9, width = 9)

