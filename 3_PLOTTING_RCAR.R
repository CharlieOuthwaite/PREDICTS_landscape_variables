##%######################################################%##
#                                                          #
####          Plotting RCAR outputs, by realm           ####
#                                                          #
##%######################################################%##

rm(list = ls())

# load libraries
library(StatisticalModels)
library(gridExtra)

# directories
datadir <- '1_PREDICTS_PLUS_VARIABLES'
moddir <- '2_MODEL_SELECTION'
outdir <- '3_PLOTTING'

# load the selected models
load(paste0(moddir, "/RCAR_Tropical_output.rdata"))
load(paste0(moddir, "/RCAR_Temperate_output.rdata"))

# read in the values used for rescaling
scalers <- read.csv("1_PREDICTS_PLUS_VARIABLES/Scaling_values.csv")

# load the function for sorting data for predictions
source('functions/sort_data.r')
source('functions/rescale.r')
source('functions/unscale.r')


# load the dataset
load(paste0(moddir, "/final.data.trans_trop_RCAR.rdata"))
load(paste0(moddir, "/final.data.trans_temp_RCAR.rdata"))




##### TROPICAL PLOTS #####

summary(rcarmod_trop$model)

# RCAR_110km ~ Predominant_land_use + 
# homogen + percNH + 
# Predominant_land_use:homogen + Predominant_land_use:percNH + 
# (1 | SS) + (1 | SSB)


#### Land use ####
PlotGLMERFactor(model = rcarmod_trop$model,data = rcarmod_trop$data,
                responseVar = "RCAR (square kilometres)",seMultiplier = 1.96,
                logLink = "e",catEffects = c("Predominant_land_use"), params = list(las = 2, cex = 0.8, mar = c(6,3,2,4), adj = 0.5), xtext.srt = 90)



#### homogeneity ####

from = 0
to = 1
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'homogen'
fac <- NULL
n <- NULL
logval = F


# organise the data
# basic table of median values and reference factors
pred_tab <- data.frame(landcovers.5k = median(final.data.rcar_trop$landcovers.5k),
                       homogen = median(final.data.rcar_trop$homogen),
                       fert.total_log = median(final.data.rcar_trop$fert.total_log),
                       percNH = median(final.data.rcar_trop$percNH),
                       Hansen_mindist_log =  median(final.data.rcar_trop$Hansen_mindist_log),
                       Predominant_land_use = "Primary vegetation",
                       RCAR_110km = 0)


levels(pred_tab$Predominant_land_use) <- levels(rcarmod_trop$data$Predominant_land_use)

vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers[scalers$variable == variable, 'centre'], 
                     scale = scalers[scalers$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab2 <- do.call("rbind", replicate(length(vals), pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab2[,variable] <- vals_trans


# predict the result
result <- PredictGLMER(model = rcarmod_trop$model, data = pred_tab2, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- 10^(result)

homogen <- as.data.frame(unscale(final.data.rcar_temp$homogen, scale = scalers[4, 2], centre = scalers[4, 3], log = F))
# SR plot = full range
ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#1874CD")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#1874CD"), alpha = 0.3) +
  geom_rug(data = homogen, aes(x = V1), size = 0.1) +
  ylim(c(15000, 1000000)) +
  xlim(c(0, 1)) +
  xlab("Homogeneity") +
  ylab("RCAR (square kilometres)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) 

ggsave(filename = paste0(outdir, "/Tropical_RCAR_homogen.pdf"))



  
#### percNH ####


from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
fac <- NULL
n <- NULL
logval = F


# basic table of median values and reference factors
pred_tab <- data.frame(landcovers.5k = median(final.data.rcar_trop$landcovers.5k),
                       homogen = median(final.data.rcar_trop$homogen),
                       fert.total_log = median(final.data.rcar_trop$fert.total_log),
                       percNH = median(final.data.rcar_trop$percNH),
                       Hansen_mindist_log =  median(final.data.rcar_trop$Hansen_mindist_log),
                       Predominant_land_use = "Primary vegetation",
                       RCAR_110km = 0)


levels(pred_tab$Predominant_land_use) <- levels(rcarmod_trop$data$Predominant_land_use)

vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers[scalers$variable == variable, 'centre'], 
                     scale = scalers[scalers$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab2 <- do.call("rbind", replicate(length(vals), pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab2[,variable] <- vals_trans


# predict the result
result <- PredictGLMER(model = rcarmod_trop$model, data = pred_tab2, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- 10^(result)


# SR plot = full range
ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#66CD00")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#66CD00"), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(15000, 1000000)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("RCAR (square kilometres") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) 

ggsave(filename = paste0(outdir, "/Tropical_RCAR_percNH.pdf"))






#### Homogen, by land use ####

from = 0
to = 1
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'homogen'
fac <- 'Predominant_land_use'
n <- 3
logval = F

# basic table of median values and reference factors
pred_tab <- data.frame(landcovers.5k = median(final.data.rcar_trop$landcovers.5k),
                       homogen = median(final.data.rcar_trop$homogen),
                       fert.total_log = median(final.data.rcar_trop$fert.total_log),
                       percNH = median(final.data.rcar_trop$percNH),
                       Hansen_mindist_log =  median(final.data.rcar_trop$Hansen_mindist_log),
                       Predominant_land_use = "Primary vegetation",
                       RCAR_110km = 0)


levels(pred_tab$Predominant_land_use) <- levels(rcarmod_trop$data$Predominant_land_use)

vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers[scalers$variable == variable, 'centre'], 
                     scale = scalers[scalers$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab2 <- do.call("rbind", replicate(length(vals), pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab2[,variable] <- vals_trans


pred_tab3 <- do.call("rbind", replicate(n, pred_tab2, simplify = FALSE))
pred_tab3[, fac][(nrow(pred_tab3)/n + 1):(nrow(pred_tab3)/n + 1000)] <- levels(final.data.rcar_trop[, fac])[2]
pred_tab3[, fac][(nrow(pred_tab3)/n + 1001):(nrow(pred_tab3)/n + 2000)] <- levels(final.data.rcar_trop[, fac])[3]
  

# predict the result
result <- PredictGLMER(model = rcarmod_trop$model, data = pred_tab3, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- 10^(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab3[, fac]

homogen <- as.data.frame(unscale(final.data.rcar_trop$homogen, scale = scalers[4, 2], centre = scalers[4, 3], log = F))
homogen$Predominant_land_use <- final.data.rcar_trop$Predominant_land_use

ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = homogen, aes(x = V1, col = Predominant_land_use), size = 0.1) +
  ylim(c(15000, 5000000)) +
  xlim(c(0, 1)) +
  xlab("Homogeneity") +
  ylab("RCAR (square kilometres)") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) 

ggsave(filename = paste0(outdir, "/Tropical_RCAR_homogenLU.pdf"))


#### percNH by land use ####

from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
fac <- 'Predominant_land_use'
n <- 3
logval = F

# basic table of median values and reference factors
pred_tab <- data.frame(landcovers.5k = median(final.data.rcar_trop$landcovers.5k),
                       homogen = median(final.data.rcar_trop$homogen),
                       fert.total_log = median(final.data.rcar_trop$fert.total_log),
                       percNH = median(final.data.rcar_trop$percNH),
                       Hansen_mindist_log =  median(final.data.rcar_trop$Hansen_mindist_log),
                       Predominant_land_use = "Primary vegetation",
                       RCAR_110km = 0)


levels(pred_tab$Predominant_land_use) <- levels(rcarmod_trop$data$Predominant_land_use)

vals_trans <- sapply(vals, 
                     FUN = rescale, 
                     centre = scalers[scalers$variable == variable, 'centre'], 
                     scale = scalers[scalers$variable == variable, 'scale'],
                     logval = logval)


# add reps of pred_tab and add the new variable vals
pred_tab2 <- do.call("rbind", replicate(length(vals), pred_tab, simplify = FALSE))

# replace the variable of interest with the new ones
pred_tab2[,variable] <- vals_trans


pred_tab3 <- do.call("rbind", replicate(n, pred_tab2, simplify = FALSE))
pred_tab3[, fac][(nrow(pred_tab3)/n + 1):(nrow(pred_tab3)/n + 1000)] <- levels(final.data.rcar_trop[, fac])[2]
pred_tab3[, fac][(nrow(pred_tab3)/n + 1001):(nrow(pred_tab3)/n + 2000)] <- levels(final.data.rcar_trop[, fac])[3]


# predict the result
result <- PredictGLMER(model = rcarmod_trop$model, data = pred_tab3, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- 10^(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab3[, fac]


percNH <- as.data.frame(unscale(final.data.rcar_trop$percNH, scale = scalers[5, 2], centre = scalers[5, 3], log = F))
percNH$Predominant_land_use <- final.data.rcar_trop$Predominant_land_use

# SR plot = full range
ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1, col = Predominant_land_use), size = 0.1) +
  ylim(c(15000, 2000000)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("RCAR (square kilometres)") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.8), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1)

ggsave(filename = paste0(outdir, "/Tropical_RCAR_percNHLU.pdf"))





##### TEMPERATE PLOTS #####

summary(rcarmod_temp$model)

# RCAR_110km ~ Predominant_land_use + Forest_biome + Use_intensity +  
# poly(homogen, 1) + fert.total_log + homogen + percNH + 
# Predominant_land_use:fert.total_log +  Predominant_land_use:homogen + Predominant_land_use:percNH + 
# Use_intensity:fert.total_log + Use_intensity:homogen + 
# (1 |      SS) + (1 | SSB)



#### Land use ####
PlotGLMERFactor(model = rcarmod_temp$model,data = rcarmod_temp$data,
                responseVar = "RCAR (square kilometres)",seMultiplier = 1.96,
                logLink = "e",catEffects = c("Predominant_land_use"), params = list(las = 2, cex = 0.8, mar = c(6,3,2,4), adj = 0.5), xtext.srt = 90)



#### Use intensity ####
PlotGLMERFactor(model = rcarmod_temp$model,data = rcarmod_temp$data,
                responseVar = "RCAR (square kilometres)",seMultiplier = 1.96,
                logLink = "e",catEffects = c( "Use_intensity"), params = list(las = 2, cex = 0.8, mar = c(6,3,2,4), adj = 0.5), xtext.srt = 90)


#### forest biomes ####
PlotGLMERFactor(model = rcarmod_temp$model,data = rcarmod_temp$data,
                responseVar = "RCAR (square kilometres)",seMultiplier = 1.96,
                logLink = "e",catEffects = "Forest_biome", params = list(las = 2, cex = 0.8, mar = c(9,5,1,4), adj = 1), xtext.srt = 50)




#### land use intensity interactions ####


# basic table of median values and reference factors
pred_tab <- data.frame(landcovers.5k = median(final.data.rcar_temp$landcovers.5k),
                       homogen = median(final.data.rcar_temp$homogen),
                       fert.total_log = median(final.data.rcar_temp$fert.total_log),
                       percNH = median(final.data.rcar_temp$percNH),
                       Hansen_mindist_log =  median(final.data.rcar_temp$Hansen_mindist_log),
                       Forest_biome = "Temperate Broadleaf & Mixed Forests",
                       Use_intensity = "Minimal use",
                       Predominant_land_use = "Primary vegetation",
                       #Tropical = "Temperate",
                       Species_richness = 0,
                       logAbun = 0, 
                       RCAR_110km = 0)


# organise factor levels
# check levels of factor variables
levels(pred_tab$Predominant_land_use) <- levels(rcarmod_temp$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(rcarmod_temp$data$Use_intensity) 
levels(pred_tab$Forest_biome) <- levels(rcarmod_temp$data$Forest_biome)

# add and change factor levels of land use and intensity

pred_tab <- do.call("rbind", replicate(9, pred_tab, simplify = FALSE))


pred_tab[4:6, 'Predominant_land_use'] <- "Secondary vegetation"
pred_tab[7:9, 'Predominant_land_use'] <- "Cropland"


pred_tab[c(2,5,8), 'Use_intensity'] <- "Light use"
pred_tab[c(3,6,9), 'Use_intensity'] <- "Intense use"


### Temperate predictions ###

# predict the result
resulta <- PredictGLMERRandIter(model = rcarmod_temp$model, data = pred_tab)

# transform the results
resulta <- 10^(resulta)

resulta <- sweep(x = resulta, MARGIN = 2, STATS = resulta[1,], FUN = '/')

resulta.median <- ((apply(X = resulta, MARGIN = 1, FUN = median))*100)-100
resulta.upper <- ((apply(X = resulta, MARGIN = 1, FUN = quantile,probs = 0.975))*100)-100
resulta.lower <- ((apply(X = resulta, MARGIN = 1, FUN = quantile,probs = 0.025))*100)-100


## plot ##

errbar.cols <- c(rep("#006400",3),rep("#8B0000", 3), rep("#EEAD0E", 3))


pdf(file = paste0(outdir, "/Temperate_RCARLUUI.pdf"))

errbar(x = 1:9,y = resulta.median,yplus = resulta.upper,yminus = resulta.lower,
       col=errbar.cols,errbar.col = errbar.cols,
       ylim=c(min(resulta.lower),max(resulta.upper)),xaxt="n",
       pch =rep(c(16,17,18), 3), 
       ylab="RCAR (square kilometres)",xlab="",bty="l")

axis(side = 1,at = c(2,5,8),
     labels = c("Primary \nvegetation","Secondary\nvegetation", "Cropland"),
     padj = 0.5)

abline(h=0,col="#00000077",lty=2)

legend("topright", 
       legend = c("Minimal Use", "Light Use", "Intense Use"),
       pch = c(16,17,18), bty = "n", inset=c(0,0))

dev.off()




#### percNH ####


from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
fac <- NULL
n <- NULL
logval = F


# organise the data
pred_tab <- sort_data(modout = rcarmod_temp,
                      moddata = final.data.rcar_temp,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = rcarmod_temp$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- 10^(result)

## organise data for plotting ##
percNH <- as.data.frame(unscale(final.data.rcar_temp$percNH, scale = scalers[5, 2], centre = scalers[5, 3], log = F))

# SR plot = full range
ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#66CD00")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#66CD00"), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(15000, 3000000)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("RCAR (square kilometres") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) 

ggsave(filename = paste0(outdir, "/Temperate_RCAR_percNH.pdf"))



#### homogeneity ####

from = 0
to = 1
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'homogen'
fac <- NULL
n <- NULL
logval = F


# organise the data
pred_tab <- sort_data(modout = rcarmod_temp,
                      moddata = final.data.rcar_temp,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = rcarmod_temp$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- 10^(result)

homogen <- as.data.frame(unscale(final.data.rcar_temp$homogen, scale = scalers[4, 2], centre = scalers[4, 3], log = F))

# SR plot = full range
ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#1874CD")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#1874CD"), alpha = 0.3) +
  geom_rug(data = homogen, aes(x = V1), size = 0.1) +
  ylim(c(15000, 4000000)) +
  xlim(c(0, 1)) +
  xlab("Homogeneity") +
  ylab("RCAR (square kilometres)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) 

ggsave(filename = paste0(outdir, "/Temperate_RCAR_homogen.pdf"))



#### Fertiliser overall ####

from = 0
to = 1000000
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'fert.total_log'
fac <- NULL
n <- NULL
logval = TRUE


# organise the data
pred_tab <- sort_data(modout = rcarmod_temp,
                      moddata = final.data.rcar_temp,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = rcarmod_temp$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- 10^(result)

# plot
ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#8B6508")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#8B6508"), alpha = 0.3) +
  geom_rug(data = final.data.rcar_temp, aes(x = fert.total), size = 0.1) +
  ylim(c(15000, 3000000)) +
  xlim(c(0, 1000000)) +
  xlab("Total fertiliser application (Kgs)") +
  ylab("RCAR (square kilometres)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) 


ggsave(filename = paste0(outdir, "/Temperate_RCAR_fert.pdf"))




#### Homogen, by land use ####

from = 0
to = 1
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'homogen'
fac <- 'Predominant_land_use'
n <- 3
logval = F

# organise the data
pred_tab <- sort_data(modout = rcarmod_temp,
                      moddata = final.data.rcar_temp,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = rcarmod_temp$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- 10^(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]

homogen <- as.data.frame(unscale(final.data.rcar_temp$homogen, scale = scalers[4, 2], centre = scalers[4, 3], log = F))
homogen$Predominant_land_use <- final.data.rcar_temp$Predominant_land_use

ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = homogen, aes(x = V1, col = Predominant_land_use), size = 0.1) +
  ylim(c(15000, 4000000)) +
  xlim(c(0, 1)) +
  xlab("Homogeneity") +
  ylab("RCAR (square kilometres)") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) 

ggsave(filename = paste0(outdir, "/Temperate_RCAR_homogenLU.pdf"))


#### percNH by land use ####

from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
fac <- 'Predominant_land_use'
n <- 3
logval = F

# organise the data
pred_tab <- sort_data(modout = rcarmod_temp,
                      moddata = final.data.rcar_temp,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = rcarmod_temp$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- 10^(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]

percNH <- as.data.frame(unscale(final.data.rcar_temp$percNH, scale = scalers[5, 2], centre = scalers[5, 3], log = F))
percNH$Predominant_land_use <- final.data.rcar_temp$Predominant_land_use

# SR plot = full range
ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1, col = Predominant_land_use), size = 0.1) +
  ylim(c(15000, 4000000)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("RCAR (square kilometres)") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.8), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1)

ggsave(filename = paste0(outdir, "/Temperate_RCAR_percNHLU.pdf"))






#### fert.total by land use ####

from = 0
to = 3000
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'fert.total_log'
fac <- 'Predominant_land_use'
n <- 3
logval = T

# organise the data
pred_tab <- sort_data(modout = rcarmod_temp,
                      moddata = final.data.rcar_temp,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = rcarmod_temp$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- 10^(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]



# SR plot = full range
ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = final.data.rcar_temp, aes(x = fert.total, col = Predominant_land_use), size = 0.1) +
  ylim(c(15000, 3000000)) +
  xlim(c(0, 3000)) +
  xlab("Total fertiliser application (Kgs)") +
  ylab("RCAR (square kilometres)") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.8), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1)

ggsave(filename = paste0(outdir, "/Temperate_RCAR_fertLU.pdf"))



#### homogen by use intensity ####

from = 0
to = 1
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'homogen'
fac <- 'Use_intensity'
n <- 3
logval = F

# organise the data
pred_tab <- sort_data(modout = rcarmod_temp,
                      moddata = final.data.rcar_temp,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = rcarmod_temp$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- 10^(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]

homogen <- as.data.frame(unscale(final.data.rcar_temp$homogen, scale = scalers[4, 2], centre = scalers[4, 3], log = F))
homogen$Use_intensity <- final.data.rcar_temp$Use_intensity

# SR plot = full range
ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = homogen, aes(x = V1, col = Use_intensity), size = 0.1) +
  ylim(c(15000, 5000000)) +
  xlim(c(0, 1)) +
  xlab("Homogeneity") +
  ylab("RCAR (square kilometres)") +
  scale_colour_manual(values = c("#66CD00", "#FFB90F", "#EE0000"))+
  scale_fill_manual(values = c("#66CD00", "#FFB90F", "#EE0000")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.9), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1) 

ggsave(filename = paste0(outdir, "/Temperate_RCAR_homogenUI.pdf"))



#### Fertilier, by use intensity ####

from = 0
to = 3000
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'fert.total_log'
fac <- 'Use_intensity'
n <- 3
logval = TRUE

# organise the data
pred_tab <- sort_data(modout = rcarmod_temp,
                      moddata = final.data.rcar_temp,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = rcarmod_temp$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- 10^(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]

# plot
ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = final.data.rcar_temp, aes(x = fert.total, col = Use_intensity), size = 0.1) +
  ylim(c(15000, 3000000)) +
  xlim(c(0, 3000)) +
  xlab("Total fertiliser application (Kgs)") +
  ylab("RCAR (square kilometres)") +
  scale_colour_manual(values = c("#66CD00", "#FFB90F", "#EE0000"))+
  scale_fill_manual(values = c("#66CD00", "#FFB90F", "#EE0000")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.9), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1) 



ggsave(filename = paste0(outdir, "/Temperate_RCAR_fertUI.pdf"))




