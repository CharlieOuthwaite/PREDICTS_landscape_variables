##%######################################################%##
#                                                          #
####                 2. Model selection                 ####
#                                                          #
##%######################################################%##

# 3. RCAR models

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
library(predictsFunctions)
library(sjPlot)


# directories
datadir <- "1_PREDICTS_PLUS_VARIABLES"
outdir <- "2_MODEL_SELECTION"


# read in the PREDICTS datasets with landscape variables
load(paste0(datadir, "/PREDICTS_dataset_inc_variables_TRANS.rdata"))

# Incorporate the RCAR data
datadir <- "0_DATA"

# Set the path to your local copy of the database
predicts.path <- paste0(datadir, "/PREDICTS_Sites_Mean_RangeSizes.rds")

# Read in the PREDICTS data
predicts <- ReadPREDICTS(predicts.path)

# Use unique SSBS values to match up the two datasets.
final.data.rcar<- merge(final.data.trans, predicts, by = "SSBS", all.x= T)

# how many of these have values
sum(!is.na(final.data.rcar$RCAR_110km)) # 5742

# remove duplicated columns
final.data.rcar <- final.data.rcar[ , c(1:28, 37)]

# remove NAs
final.data.rcar <- final.data.rcar[!is.na(final.data.rcar$RCAR_110km), ] # 5742 rows

# remove the .x from column names
colnames(final.data.rcar) <- sub(".x", "", colnames(final.data.rcar))

final.data.rcar <- droplevels(final.data.rcar)


# Split the dataset based on realm

final.data.rcar_trop <- final.data.rcar[final.data.rcar$Tropical == "Tropical", ]
nrow(final.data.rcar_trop) # 1988
final.data.rcar_trop <- droplevels(final.data.rcar_trop )
  
final.data.rcar_temp <- final.data.rcar[final.data.rcar$Tropical == "Temperate", ]
nrow(final.data.rcar_temp) # 3754
final.data.rcar_temp <- droplevels(final.data.rcar_temp)
  
save(final.data.rcar_trop, file = paste0(outdir, "/final.data.trans_trop_RCAR.rdata"))
save(final.data.rcar_temp, file = paste0(outdir, "/final.data.trans_temp_RCAR.rdata"))


##%######################################################%##
#                                                          #
####              Model selection process               ####
####         including all landscape variables          ####
#                                                          #
##%######################################################%##


#### 1. Tropical subset ####

# run the model selection process


system.time({rcar_trop <- GLMERSelect(modelData = final.data.rcar_trop, responseVar = "RCAR_110km",
                                  fitFamily = "gaussian", 
                                  fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity"),
                                  fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1, percNH = 1),
                                  randomStruct = "(1|SS)+(1|SSB)", 
                                  fixedInteractions = c("Predominant_land_use:fert.total_log", 
                                                        "Predominant_land_use:Hansen_mindist_log", 
                                                        "Predominant_land_use:landcovers.5k", 
                                                        "Predominant_land_use:homogen", 
                                                        "Predominant_land_use:percNH",
                                                        "Use_intensity:fert.total_log", 
                                                        "Use_intensity:Hansen_mindist_log", 
                                                        "Use_intensity:landcovers.5k",
                                                        "Use_intensity:homogen",
                                                        "Use_intensity:percNH",
                                                        "Predominant_land_use:Use_intensity"), verbose = F)}
)

# take a look at the model output
summary(rcar_trop$model)

# save the model output
save(rcar_trop, file = paste0(outdir, "/RCAR_Tropical_Model_Selection.rdata"))

# extract the stats produced as part of the model selection process
rcar_trop_stats <- as.data.frame(rcar_trop$stats)

# save these
write.csv(rcar_trop_stats, file = paste0(outdir, "/rcar_trop_stats.csv"), row.names = F)





#### 2. Temperate subset ####

# run the model selection process

system.time({rcar_temp <- GLMERSelect(modelData = final.data.rcar_temp, responseVar = "RCAR_110km",
                                      fitFamily = "gaussian", 
                                      fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity"),
                                      fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1, percNH = 1),
                                      randomStruct = "(1|SS)+(1|SSB)", 
                                      fixedInteractions = c("Predominant_land_use:fert.total_log", 
                                                            "Predominant_land_use:Hansen_mindist_log", 
                                                            "Predominant_land_use:landcovers.5k", 
                                                            "Predominant_land_use:homogen", 
                                                            "Predominant_land_use:percNH",
                                                            "Use_intensity:fert.total_log", 
                                                            "Use_intensity:Hansen_mindist_log", 
                                                            "Use_intensity:landcovers.5k",
                                                            "Use_intensity:homogen",
                                                            "Use_intensity:percNH",
                                                            "Predominant_land_use:Use_intensity"), verbose = F)}
)

# take a look at the model output
summary(rcar_temp$model)

# save the model output
save(rcar_temp, file = paste0(outdir, "/RCAR_Temperate_Model_Selection.rdata"))

# extract the stats produced as part of the model selection process
rcar_temp_stats <- as.data.frame(rcar_temp$stats)

# save these
write.csv(rcar_temp_stats, file = paste0(outdir, "/rcar_temp_stats.csv"), row.names = F)



##%######################################################%##
#                                                          #
####          Rerun selected models using REML          ####
#                                                          #
##%######################################################%##


# load(paste0(outdir, "/RCAR_Tropical_Model_Selection.rdata"))
# load(paste0(outdir, "/RCAR_Temperate_Model_Selection.rdata"))


# 1. Tropical
summary(rcar_trop$model)

# RCAR_110km ~ Predominant_land_use + 
# homogen + percNH + 
# Predominant_land_use:homogen + Predominant_land_use:percNH +
# (1 | SS) + (1 | SSB)


# run selected model with REML

mod_struc <- "Predominant_land_use + homogen + percNH + Predominant_land_use:homogen + Predominant_land_use:percNH"

rcarmod_trop <- GLMER(modelData = final.data.rcar_trop, responseVar = "RCAR_110km", fitFamily = "gaussian",
                 fixedStruct = mod_struc,
                 randomStruct = "(1|SS) + (1|SSB)", REML = TRUE)

# take a look at the model output
summary(rcarmod_trop$model)

# extract the coefficents of the model
coefs_trop <- fixef(rcarmod_trop$model)

# save the coefficients
write.csv(coefs_trop, file = paste0(outdir, "/RCAR_Tropical_coefs.csv"), row.names = F)

# save the model output
save(rcarmod_trop, file = paste0(outdir, "/RCAR_Tropical_output.rdata"))



# 2. Temperate
summary(rcar_temp$model)


# RCAR_110km ~ Predominant_land_use + Forest_biome + Use_intensity +      
# poly(homogen, 1) + fert.total_log +      homogen + percNH +
# Predominant_land_use:fert.total_log +  Predominant_land_use:homogen + Predominant_land_use:percNH +  
# Use_intensity:fert.total_log + Use_intensity:homogen +  
# (1 | SS) + (1 | SSB)

# fixed-effect model matrix is rank deficient so dropping 1 column / coefficient

mod_struc <- "Predominant_land_use + Forest_biome + Use_intensity +  poly(homogen, 1) + fert.total_log + homogen + percNH + Predominant_land_use:fert.total_log +  Predominant_land_use:homogen + Predominant_land_use:percNH +  Use_intensity:fert.total_log + Use_intensity:homogen"
  
  
# run selected model with REML

rcarmod_temp <- GLMER(modelData = final.data.rcar_temp, responseVar = "RCAR_110km", fitFamily = "gaussian",
                      fixedStruct = mod_struc,
                      randomStruct = "(1|SS) + (1|SSB)", REML = TRUE)

# fixed-effect model matrix is rank deficient so dropping 1 column / coefficient


# take a look at the model output
summary(rcarmod_temp$model)

# extract the coefficents of the model
coefs_temp <- fixef(rcarmod_temp$model)

# save the coefficients
write.csv(coefs_temp, file = paste0(outdir, "/RCAR_Temperate_coefs.csv"), row.names = F)

# save the model output
save(rcarmod_temp, file = paste0(outdir, "/RCAR_Temperate_output.rdata"))



##%######################################################%##
#                                                          #
####                  R-squared values                  ####
#                                                          #
##%######################################################%##


rcar_trop_r2 <- R2GLMER(rcarmod_trop$model) # [1] conditional, [2] marginal
rcar_temp_r2 <- R2GLMER(rcarmod_temp$model)

# Look at the proportion of variance unexplained by the random effects, but
# explained by the fixed effects...  Marginal/(1 - (Conditional - Marginal))


rcarper_trop <- rcar_trop_r2[[2]]/(1-(rcar_trop_r2[[1]] - rcar_trop_r2[[2]]))
rcarper_temp <- rcar_temp_r2[[2]]/(1-(rcar_temp_r2[[1]] - rcar_temp_r2[[2]]))


# organise values in a table

result <- NULL

rcarres_trop <- unlist(c(rcar_trop_r2[1], rcar_trop_r2[2], rcarper_trop))
rcarres_temp  <- unlist(c(rcar_temp_r2[1], rcar_temp_r2[2], rcarper_temp))

result <- rbind(rcarres_trop, rcarres_temp)

colnames(result)[3] <- "proportion unexplained variance explained by fixed effects"

rownames(result) <- c("RCAR Trop", "RCAR Temp")

# save
write.csv(result, file = paste0(outdir, "/Rsquareds_RCAR.csv"))


# save model output tables

tab_model(rcarmod_trop$model, transform = NULL, file = paste0(outdir, "/RCAR_Trop_output_table.html"))
tab_model(rcarmod_temp$model, transform = NULL, file = paste0(outdir, "/RCAR_Temp_output_table.html"))


##%######################################################%##
#                                                          #
####                    Model checks                    ####
#                                                          #
##%######################################################%##



### Tropical ###

## 1. Checking the fitted vs residuals relationship
p1 <- plot(rcarmod_trop$model)


## 2. Normality of Residuals

pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(rcarmod_trop$model), main = "", ylim = c(-4,4))
qqline(resid(rcarmod_trop$model))
p2 <- recordPlot()
invisible(dev.off())



## 3. Check for spatial autocorrelation

rcar_dat <- final.data.rcar_trop[!is.na(final.data.rcar_trop$RCAR_110km), ]


# needed for following function
rcarmod_trop$data$SSBS <- final.data.rcar_trop$SSBS
rcarmod_trop$data$Latitude <- final.data.rcar_trop$Latitude
rcarmod_trop$data$Longitude <- final.data.rcar_trop$Longitude


rcar_test<-roquefort::SpatialAutocorrelationTest(model=rcarmod_trop, all.data=final.data.rcar_trop)


summary(rcar_test)

# percentage of studies that show spatial autocorrelation?
perc_auto <- (length(which(rcar_test$P<0.05))/length(rcar_test$P))*100

# 5.09%

rcar_test_vals <- as.data.frame(rcar_test$P)
rcar_test_vals$`ab_test$P` <- round(rcar_test_vals$`rcar_test$P`, digits = 4)

label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")

p3 <- ggplot(data = rcar_test_vals ) +
  geom_histogram(aes(x = rcar_test_vals$`ab_test$P`)) +
  geom_vline(xintercept = 0.05, col = "red") +
  geom_text(aes(x = 0.2, y = 90, label = label1), size = 4, check_overlap = T) +
  theme_bw() +
  xlab("P-value") +
  ylab("Frequency") +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1)



cowplot::plot_grid(p1,p2,p3,
          labels = c("A.", "B.", "C."))


ggsave(file = paste0(outdir, "/RCAR__Tropical_model_checks_plots.pdf"), height = 9, width = 9)




### Temperate ###

## 1. Checking the fitted vs residuals relationship
p1 <- plot(rcarmod_temp$model)


## 2. Normality of Residuals

pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(rcarmod_temp$model), main = "", ylim = c(-4,4))
qqline(resid(rcarmod_temp$model))
p2 <- recordPlot()
invisible(dev.off())



## 3. Check for spatial autocorrelation

rcar_dat <- final.data.rcar_temp[!is.na(final.data.rcar_temp$RCAR_110km), ]


# needed for following function
rcarmod_temp$data$SSBS <- final.data.rcar_temp$SSBS
rcarmod_temp$data$Latitude <- final.data.rcar_temp$Latitude
rcarmod_temp$data$Longitude <- final.data.rcar_temp$Longitude


rcar_test<-roquefort::SpatialAutocorrelationTest(model=rcarmod_temp, all.data=final.data.rcar_temp)


summary(rcar_test)

# percentage of studies that show spatial autocorrelation?
perc_auto <- (length(which(rcar_test$P<0.05))/length(rcar_test$P))*100

# 5%

rcar_test_vals <- as.data.frame(rcar_test$P)
rcar_test_vals$`ab_test$P` <- round(rcar_test_vals$`rcar_test$P`, digits = 4)

label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")

p3 <- ggplot(data = rcar_test_vals ) +
  geom_histogram(aes(x = rcar_test_vals$`ab_test$P`)) +
  geom_vline(xintercept = 0.05, col = "red") +
  geom_text(aes(x = 0.2, y = 90, label = label1), size = 4, check_overlap = T) +
  theme_bw() +
  xlab("P-value") +
  ylab("Frequency") +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1)



cowplot::plot_grid(p1,p2,p3,
          labels = c("A.", "B.", "C."))


ggsave(file = paste0(outdir, "/RCAR__Temperate_model_checks_plots.pdf"), height = 9, width = 9)

