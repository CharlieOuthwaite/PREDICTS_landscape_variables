##%######################################################%##
#                                                          #
####                 2. Model selection                 ####
#                                                          #
##%######################################################%##

# This script carries out the model selection process using the dataset generated
# in the previous script. Models are run for 3 biodiversity metrics. 

rm(list = ls())


# load libraries
library(devtools)
install_github(repo = "timnewbold/StatisticalModels")
library(StatisticalModels)
library(roquefort)
library(cowplot)
library(gridGraphics)


# directories
datadir <- "1_PREDICTS_PLUS_VARIABLES"
outdir <- "2_MODEL_SELECTION"

# read in the PREDICTS datasets with landscape variables

load(paste0(datadir, "/PREDICTS_dataset_inc_variables_TRANS.rdata"))
load(paste0(datadir, "/PREDICTS_abun_subset.rdata"))
load(paste0(datadir, "/PREDICTS_rcar_subset.rdata"))



#### Run the model selection process including all landscape variables ####


##%######################################################%##
#                                                          #
####                1. Species richness                 ####
#                                                          #
##%######################################################%##


#### 1. Species richness models ####

# run the model selection process

system.time({sr1 <- GLMERSelect(modelData = final.data.trans, 
                                  responseVar = "Species_richness",
                                  fitFamily = "poisson", 
                                  fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity", "Tropical"),
                                  fixedTerms = list(fert.total_log = 1, Hansen_mindist_log = 1, landcovers.5k = 1, homogen = 1, percNH = 1),
                                  randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)", 
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
                                                        "Predominant_land_use:Use_intensity",
                                                        "Tropical:Hansen_mindist_log"), verbose = F)}) 
# save the output
save(sr1, file = paste0(outdir, "/SPECIESRICHNESS_Model_selection.rdata"))

# take a look at the model output
summary(sr1$model)

# extract the stats produced as part of the model selection process
sr1stats <- as.data.frame(sr1$stats)

# save these
write.csv(sr1stats, file = paste0(outdir, "/sr_stats.csv"), row.names = F)


##%######################################################%##
#                                                          #
####                    2. Abundance                    ####
#                                                          #
##%######################################################%##


system.time({ab1 <- GLMERSelect(modelData = final.data.abun, 
                                  responseVar = "logAbun",
                                  fitFamily = "gaussian", 
                                  fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity", "Tropical"),
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
                                                        "Predominant_land_use:Use_intensity", 
                                                        "Tropical:Hansen_mindist_log"), verbose = F)})

# take a look at the model output
summary(ab1$model)

# save the output
save(ab1, file = paste0(outdir, "/ABUNDANCE_Model_Selection.rdata"))

# extract the stats produced as part of the model selection process
ab1stats <- as.data.frame(ab1$stats)

# save these
write.csv(ab1stats, file = paste0(outdir, "/ab_stats.csv"), row.names = F)


##%######################################################%##
#                                                          #
####                      3. RCAR                       ####
#                                                          #
##%######################################################%##


system.time({rcar1 <- GLMERSelect(modelData = final.data.rcar, responseVar = "RCAR_110km",
                                    fitFamily = "gaussian", 
                                    fixedFactors = c("Predominant_land_use", "Forest_biome", "Use_intensity", "Tropical"),
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
                                                          "Predominant_land_use:Use_intensity",
                                                          "Tropical:Hansen_mindist_log"), verbose = F)}
)

# take a look at the model output
summary(rcar1$model)

# save the model output
save(rcar1, file = paste0(outdir, "/RCAR_Model_Selection.rdata"))

# extract the stats produced as part of the model selection process
rcar1stats <- as.data.frame(rcar1$stats)

# save these
write.csv(rcar1stats, file = paste0(outdir, "/rcar_stats.csv"), row.names = F)




##%######################################################%##
#                                                          #
####        Run final selected model using REML         ####
#                                                          #
##%######################################################%##


# selected model for species richness

# Species_richness ~ Predominant_land_use + Forest_biome + Use_intensity + Tropical + 
# percNH + Hansen_mindist_log + fert.total_log + landcovers.5k + homogen + 
# Predominant_land_use:fert.total_log + Predominant_land_use:landcovers.5k + Predominant_land_use:homogen +  
# Use_intensity:fert.total_log + Use_intensity:percNH + 
# Predominant_land_use:Use_intensity + Tropical:Hansen_mindist_log +  
# (1 | SS) + (1 | SSB) + (1 | SSBS)

# rerun the selected models, using REML

srmod <- GLMER(modelData = final.data.trans, responseVar = "Species_richness", fitFamily = "poisson",
               fixedStruct = "Predominant_land_use + Forest_biome + Use_intensity + Tropical + Hansen_mindist_log + percNH + fert.total_log + landcovers.5k + homogen + Predominant_land_use:fert.total_log + Predominant_land_use:landcovers.5k + Predominant_land_use:homogen +  Use_intensity:fert.total_log + Use_intensity:percNH + Predominant_land_use:Use_intensity + Tropical:Hansen_mindist_log",
               randomStruct = "(1|SS) + (1|SSB) + (1|SSBS)", REML = TRUE)


# Warning messages:
# 1: In commonArgs(par, fn, control, environment()) :
# maxfun < 10 * length(par)^2 is not recommended.
# 2: In optwrap(optimizer, devfun, start, rho$lower, control = control,  :
# convergence code 1 from bobyqa: bobyqa -- maximum number of function evaluations exceeded
# 3: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
# Model failed to converge with max|grad| = 0.0108372 (tol = 0.001, component 1)


# take a look at the model output
summary(srmod$model)

# extract the coefficents of the model
coefs <- fixef(srmod$model)

# save the coefficients
write.csv(coefs, file = paste0(outdir, "/SPECIESRICHNESS_coefs.csv"), row.names = F)

# save the model output
save(srmod, file = paste0(outdir, "/SRMOD_output.rdata"))


# selected model for abundance

# logAbun ~ Predominant_land_use + Forest_biome + Use_intensity +  Tropical + 
# Hansen_mindist_log + percNH + fert.total_log + landcovers.5k + homogen +
# Predominant_land_use:fert.total_log + Predominant_land_use:percNH + 
# Use_intensity:fert.total_log +  Use_intensity:percNH + 
# Predominant_land_use:Use_intensity + Tropical:Hansen_mindist_log + 
# (1 | SS) + (1 | SSB)

# run selected model with REML

abmod <- GLMER(modelData = final.data.abun, responseVar = "logAbun", fitFamily = "gaussian",
               fixedStruct = "Predominant_land_use + Forest_biome + Use_intensity + Tropical + Hansen_mindist_log + percNH + fert.total_log + landcovers.5k + homogen + Predominant_land_use:fert.total_log + Predominant_land_use:percNH + Use_intensity:fert.total_log + Use_intensity:percNH + Predominant_land_use:Use_intensity + Tropical:Hansen_mindist_log",
               randomStruct = "(1|SS) + (1|SSB)", REML = TRUE)


# take a look at the model output
summary(abmod$model)

# extract the coefficents of the model
coefs <- fixef(abmod$model)

# save the coefficients
write.csv(coefs, file = paste0(outdir, "/ABUNDANCE_coefs.csv"), row.names = F)

# save the model output
save(abmod, file = paste0(outdir, "/ABMOD_output.rdata"))


# selected model for RCAR

# RCAR_110km ~ Predominant_land_use + Forest_biome + + Use_intensity +
# Hansen_mindist_log + percNH  + homogen +  landcovers.5k +
# Predominant_land_use:homogen + Predominant_land_use:percNH +  
# Use_intensity:Hansen_mindist_log + Use_intensity:landcovers.5k +
# Use_intensity:homogen + Predominant_land_use:Use_intensity +  
# (1 | SS) + (1 | SSB)

# run selected model with REML

rcarmod <- GLMER(modelData = final.data.rcar, responseVar = "RCAR_110km", fitFamily = "gaussian",
                 fixedStruct = "Predominant_land_use + Forest_biome + Use_intensity + Hansen_mindist_log + percNH + homogen + landcovers.5k + Predominant_land_use:percNH + Predominant_land_use:homogen + Use_intensity:Hansen_mindist_log + Use_intensity:landcovers.5k + Use_intensity:homogen + Predominant_land_use:Use_intensity",
                 randomStruct = "(1|SS) + (1|SSB)", REML = TRUE)

# take a look at the model output
summary(rcarmod$model)

# extract the coefficents of the model
coefs <- fixef(rcarmod$model)

# save the coefficients
write.csv(coefs, file = paste0(outdir, "/RCAR_coefs.csv"), row.names = F)

# save the model output
save(rcarmod, file = paste0(outdir, "/RCAR_output.rdata"))



##%######################################################%##
#                                                          #
####                  R-squared values                  ####
#                                                          #
##%######################################################%##

srr2 <- R2GLMER(srmod$model) # [1] conditional, [2] marginal
abr2 <- R2GLMER(abmod$model)
rcarr2 <- R2GLMER(rcarmod$model)

# Look at the percentage of variance unexplained by the random effects, but
# explained by the fixed effects...  Marginal/(1 - (Conditional - Marginal))


srper <- srr2[[2]]/(1-(srr2[[1]] - srr2[[2]]))*100
abper <- abr2[[2]]/(1-(abr2[[1]] - abr2[[2]]))*100
rcarper <- rcarr2[[2]]/(1-(rcarr2[[1]] - rcarr2[[2]]))*100


# organise values in a table

result <- NULL

srres  <- unlist(c(srr2[1], srr2[2], srper))
abres  <- unlist(c(abr2[1], abr2[2], abper))
rcarres  <- unlist(c(rcarr2[1], rcarr2[2], rcarper))

result <- rbind(srres, abres, rcarres)

colnames(result)[3] <- "% unexplained variance explained by fixed effects"

rownames(result) <- c("SR", "Abun", "RCAR")

# save
write.csv(result, file = paste0(outdir, "/Rsquareds.csv"))



# save model output tables
library(sjPlot)
library(lme4)

tab_model(srmod$model, transform = NULL, file = paste0(outdir, "/SR_output_table"))
tab_model(abmod$model, transform = NULL, file = paste0(outdir, "/AB_output_table"))
tab_model(rcarmod$model, transform = NULL, file = paste0(outdir, "/RCAR_output_table"))



##%######################################################%##
#                                                          #
####                    Model checks                    ####
#                                                          #
##%######################################################%##


### SR models ###

## 1. Checking the fitted vs residuals relationship
p1 <- plot(srmod$model)


## 2. Normality of Residuals

pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(srmod$model), main = "")
p2 <- recordPlot()
invisible(dev.off())



## 3. Check for spatial autocorrelation

sr_test<-roquefort::SpatialAutocorrelationTest(model=srmod, all.data=final.data.trans)


summary(sr_test)

# percentage of studies that show spatial autocorrelation?
perc_auto <- (length(which(sr_test$P<0.05))/length(sr_test$P))*100

# 5.29%

sr_test_vals <- as.data.frame(sr_test$P)
sr_test_vals$`sr_test$P` <- round(sr_test_vals$`sr_test$P`, digits = 4)

label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")

p3 <- ggplot(data = sr_test_vals ) +
  geom_histogram(aes(x = sr_test_vals$`sr_test$P`)) +
  geom_vline(xintercept = 0.05, col = "red") +
  geom_text(aes(x = 0.2, y = 125, label = label1), size = 4, check_overlap = T) +
  theme_bw() +
  xlab("P-value") +
  ylab("Frequency") +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1)



plot_grid(p1,p2,p3,
          labels = c("A.", "B.", "C."))


ggsave(file = paste0(outdir, "/SR_model_checks_plots.pdf"), height = 9, width = 9)



### Abundance plots ###


## 1. Checking the fitted vs residuals relationship
p1 <- plot(abmod$model)


## 2. Normality of Residuals

pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(abmod$model), main = "")
p2 <- recordPlot()
invisible(dev.off())



## 3. Check for spatial autocorrelation

# needed for following function
abmod$data$SSBS <- final.data.abun$SSBS
abmod$data$Latitude <- final.data.abun$Latitude
abmod$data$Longitude <- final.data.abun$Longitude


ab_test<-roquefort::SpatialAutocorrelationTest(model=abmod, all.data=final.data.abun)


summary(ab_test)

# percentage of studies that show spatial autocorrelation?
perc_auto <- (length(which(ab_test$P<0.05))/length(ab_test$P))*100

# 9.4%

ab_test_vals <- as.data.frame(ab_test$P)
ab_test_vals$`ab_test$P` <- round(sr_test_vals$`ab_test$P`, digits = 4)

label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")

p3 <- ggplot(data = ab_test_vals ) +
  geom_histogram(aes(x = ab_test_vals$`ab_test$P`)) +
  geom_vline(xintercept = 0.05, col = "red") +
  geom_text(aes(x = 0.2, y = 90, label = label1), size = 4, check_overlap = T) +
  theme_bw() +
  xlab("P-value") +
  ylab("Frequency") +
  theme(panel.grid = element_blank(), 
        aspect.ratio = 1)



plot_grid(p1,p2,p3,
          labels = c("A.", "B.", "C."))


ggsave(file = paste0(outdir, "/Abun_model_checks_plots.pdf"), height = 9, width = 9)



### RCAR model ###



## 1. Checking the fitted vs residuals relationship
p1 <- plot(rcarmod$model)


## 2. Normality of Residuals

pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(rcarmod$model), main = "", ylim = c(-4,4))
p2 <- recordPlot()
invisible(dev.off())



## 3. Check for spatial autocorrelation

# needed for following function
rcarmod$data$SSBS <- final.data.rcar$SSBS
rcarmod$data$Latitude <- final.data.rcar$Latitude
rcarmod$data$Longitude <- final.data.rcar$Longitude


rcar_test<-roquefort::SpatialAutocorrelationTest(model=rcarmod, all.data=final.data.rcar)


summary(rcar_test)

# percentage of studies that show spatial autocorrelation?
perc_auto <- (length(which(rcar_test$P<0.05))/length(rcar_test$P))*100

# 5.594406%

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



plot_grid(p1,p2,p3,
          labels = c("A.", "B.", "C."))


ggsave(file = paste0(outdir, "/RCAR_model_checks_plots.pdf"), height = 9, width = 9)




