##%######################################################%##
#                                                          #
####        Plotting Abundance outputs, by realm        ####
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
load(paste0(moddir, "/ABMOD_Tropical_output.rdata"))
load(paste0(moddir, "/ABMOD_Temperate_output.rdata"))

# read in the values used for rescaling
scalers <- read.csv("1_PREDICTS_PLUS_VARIABLES/Scaling_values.csv")

# load the function for sorting data for predictions
source('functions/sort_data.r')
source('functions/rescale.r')
source('functions/unscale.r')


# load the dataset
load(paste0(moddir, "/final.data.trans_trop_ABUN.rdata"))
load(paste0(moddir, "/final.data.trans_temp_ABUN.rdata"))


##### TROPICAL PLOTS #####

summary(abmod_trop$model)

# logAbun ~ Predominant_land_use + Forest_biome + Use_intensity +  
# poly(homogen, 1) + Hansen_mindist_log + landcovers.5k + 
# Predominant_land_use:Hansen_mindist_log + Use_intensity:landcovers.5k + 
# (1 | SS) + (1 | SSB)

#### Land use ####
PlotGLMERFactor(model = abmod_trop$model,data = abmod_trop$data,
                responseVar = "Total Abundance",seMultiplier = 1.96,
                logLink = "e",catEffects = c("Predominant_land_use"), params = list(las = 2, cex = 0.8, mar = c(6,3,2,4), adj = 0.5), xtext.srt = 90)



#### Use intensity ####
PlotGLMERFactor(model = abmod_trop$model,data = abmod_trop$data,
                responseVar = "Total Abundance",seMultiplier = 1.96,
                logLink = "e",catEffects = c( "Use_intensity"), params = list(las = 2, cex = 0.8, mar = c(6,3,2,4), adj = 0.5), xtext.srt = 90)


#### forest biomes ####
PlotGLMERFactor(model = abmod_trop$model,data = abmod_trop$data,
                responseVar = "Total Abundance",seMultiplier = 1.96,
                logLink = "e",catEffects = "Forest_biome", params = list(las = 2, cex = 0.8, mar = c(9,5,1,4), adj = 1), xtext.srt = 50)



#### land use intensity interactions ####


# basic table of median values and reference factors
pred_tab <- data.frame(landcovers.5k = median(final.data.trans_trop_ABUN$landcovers.5k),
                       homogen = median(final.data.trans_trop_ABUN$homogen),
                       fert.total_log = median(final.data.trans_trop_ABUN$fert.total_log),
                       percNH = median(final.data.trans_trop_ABUN$percNH),
                       Hansen_mindist_log =  median(final.data.trans_trop_ABUN$Hansen_mindist_log),
                       Forest_biome = "Tropical & Subtropical Moist Broadleaf Forests",
                       Use_intensity = "Minimal use",
                       Predominant_land_use = "Primary vegetation",
                       #Tropical = "Temperate",
                       Species_richness = 0,
                       logAbun = 0, 
                       RCAR_110km = 0)


# organise factor levels
# check levels of factor variables
levels(pred_tab$Predominant_land_use) <- levels(abmod_trop$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(abmod_trop$data$Use_intensity) 
levels(pred_tab$Forest_biome) <- levels(abmod_trop$data$Forest_biome)[c(3, 2, 1)]

# add and change factor levels of land use and intensity

pred_tab <- do.call("rbind", replicate(9, pred_tab, simplify = FALSE))


pred_tab[4:6, 'Predominant_land_use'] <- "Secondary vegetation"
pred_tab[7:9, 'Predominant_land_use'] <- "Cropland"



pred_tab[c(2,5,8), 'Use_intensity'] <- "Light use"
pred_tab[c(3,6,9), 'Use_intensity'] <- "Intense use"


### Tropical predictions ###

# predict the result
resulta <- PredictGLMERRandIter(model = abmod_trop$model, data = pred_tab)

# transform the results
resulta <- exp(resulta)-1

resulta <- sweep(x = resulta, MARGIN = 2, STATS = resulta[1,], FUN = '/')

resulta.median <- ((apply(X = resulta, MARGIN = 1, FUN = median))*100)-100
resulta.upper <- ((apply(X = resulta, MARGIN = 1, FUN = quantile,probs = 0.975))*100)-100
resulta.lower <- ((apply(X = resulta, MARGIN = 1, FUN = quantile,probs = 0.025))*100)-100


## plot ##

errbar.cols <- c(rep("#006400",3),rep("#8B0000", 3), rep("#EEAD0E", 3))


pdf(file = paste0(outdir, "/Tropical_Abun_LUUI.pdf"))

errbar(x = 1:9,y = resulta.median,yplus = resulta.upper,yminus = resulta.lower,
       col=errbar.cols,errbar.col = errbar.cols,
       ylim=c(-50,50),xaxt="n",
       pch =rep(c(16,17,18), 3), 
       ylab="Total Abundance (%)",xlab="",bty="l", cex.lab =1.2, cex.axis = 1.2)

axis(side = 1,at = c(2,5,8),
     labels = c("Primary \nvegetation","Secondary\nvegetation", "Cropland"),
     padj = 0.5, cex.lab =1.2)

abline(h=0,col="#00000077",lty=2)

legend("topleft", 
       legend = c("Minimal Use", "Light Use", "Intense Use"),
       pch = c(16,17,18), bty = "n", inset=c(0,0), cex = 1.2)

dev.off()



#### Distance to forest ####

from = 0
to = max(final.data.trans_trop_ABUN$Hansen_mindist)
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'Hansen_mindist_log'
fac <- NULL
n <- NULL
logval = TRUE


# organise the data
pred_tab <- sort_data(modout = abmod_trop,
                      moddata = final.data.trans_trop_ABUN,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n, 
                      log = logval)

# predict the result
result <- PredictGLMER(model = abmod_trop$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1

## organise data for plotting ##


# SR plot = full range
ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#458B00")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#458B00"), alpha = 0.3) +
  geom_rug(data = final.data.trans_trop_ABUN, aes(x = Hansen_mindist), size = 0.1) +
  ylim(c(0, 200)) +
  xlim(c(0, 100)) +
  xlab("Distance to Forest (Km)") +
  ylab("Total Abundance") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) 

ggsave(filename = paste0(outdir, "/Tropical_Abun_dist.pdf"))



#### number of landcovers ####

from = 1
to = 10
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'landcovers.5k'
fac <- NULL
n <- NULL
logval = F


# organise the data
pred_tab <- sort_data(modout = abmod_trop,
                      moddata = final.data.trans_trop_ABUN,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = abmod_trop$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1


# SR plot = full range
ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#CD950C")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#CD950C"), alpha = 0.3) +
  ylim(c(0,100)) +
  xlim(c(0, 11)) +
  xlab("Number of Landcovers") +
  ylab("Total Abundance") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  scale_x_continuous(breaks = c(seq(from = 0, to = 10, by = 2)))

ggsave(filename = paste0(outdir, "/Tropical_Abun_nlandcovers.pdf"))



#### homogeneity ####

from = 0
to = 1
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'homogen'
fac <- NULL
n <- NULL
logval = F


# organise the data
pred_tab <- sort_data(modout = abmod_trop,
                      moddata = final.data.trans_trop_ABUN,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = abmod_trop$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1

homogen <- as.data.frame(unscale(final.data.trans_trop_ABUN$homogen, scale = scalers[4, 2], centre = scalers[4, 3], log = F))

# SR plot = full range
ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#1874CD")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#1874CD"), alpha = 0.3) +
  geom_rug(data = homogen, aes(x = V1), size = 0.1) +
  ylim(c(0,200)) +
  xlim(c(0, 1)) +
  xlab("Homogeneity") +
  ylab("Total Abundance") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) 


ggsave(filename = paste0(outdir, "/Tropical_Abun_homogen.pdf"))



#### nlandcovers by UI ####

from = 0
to = 10
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'landcovers.5k'
fac <- 'Use_intensity'
n <- 3
logval = F

# organise the data
pred_tab <- sort_data(modout = abmod_trop,
                      moddata = final.data.trans_trop_ABUN,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = abmod_trop$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


# SR plot = full range
ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  #geom_rug(data = landcovers.5k, aes(x = V1)) +
  ylim(c(0, 200)) +
  xlim(c(0, 10)) +
  xlab("Number of Landcovers") +
  ylab("Total Abundance") +
  scale_colour_manual(values = c("#66CD00", "#FFB90F", "#EE0000"))+
  scale_fill_manual(values = c("#66CD00", "#FFB90F", "#EE0000")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.9), legend.title = element_blank(),
        legend.text = element_text(size = 10),
        aspect.ratio = 1, legend.background = element_blank()) 


ggsave(filename = paste0(outdir, "/Tropical_Abun_landcoversUI.pdf"))


#### distance by LU ####


from = 0
to = 80
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'Hansen_mindist_log'
fac <- 'Predominant_land_use'
n <- 3
logval = TRUE



# organise the data
pred_tab <- sort_data(modout = abmod_trop,
                      moddata = final.data.trans_trop_ABUN,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = abmod_trop$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1


## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]

# SR plot = full range
ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = final.data.trans_trop_ABUN, aes(x = Hansen_mindist, col = Predominant_land_use), size = 0.1) +
  ylim(c(0,250)) +
  xlim(c(0, 100)) +
  xlab("Distance to Forest (Km)") +
  ylab("Total Abundance") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.3,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 9),
        aspect.ratio = 1, legend.background = element_blank(),
        text = element_text(size = 14)) 

ggsave(filename = paste0(outdir, "/Tropical_Abun_distLU.pdf"))



##### TEMPERATE PLOTS #####

summary(abmod_temp$model)

# logAbun ~ Predominant_land_use + Forest_biome + Use_intensity +  
# poly(fert.total_log, 1) + poly(landcovers.5k, 1) + poly(homogen,1) + poly(percNH, 1) + fert.total_log + 
# Use_intensity:fert.total_log +  Predominant_land_use:Use_intensity
# (1 | SS) + (1 | SSB)

#### Land use ####
PlotGLMERFactor(model = abmod_temp$model,data = abmod_temp$data,
                responseVar = "Total Abundance",seMultiplier = 1.96,
                logLink = "e",catEffects = c("Predominant_land_use"), params = list(las = 2, cex = 0.8, mar = c(6,3,2,4), adj = 0.5), xtext.srt = 90)



#### Use intensity ####
PlotGLMERFactor(model = abmod_temp$model,data = abmod_temp$data,
                responseVar = "Total Abundance",seMultiplier = 1.96,
                logLink = "e",catEffects = c( "Use_intensity"), params = list(las = 2, cex = 0.8, mar = c(6,3,2,4), adj = 0.5), xtext.srt = 90)


#### forest biomes ####
PlotGLMERFactor(model = abmod_temp$model,data = abmod_temp$data,
                responseVar = "Total Abundance",seMultiplier = 1.96,
                logLink = "e",catEffects = "Forest_biome", params = list(las = 2, cex = 0.8, mar = c(9,5,1,4), adj = 1), xtext.srt = 50)




#### land use intensity interactions ####


# basic table of median values and reference factors
pred_tab <- data.frame(landcovers.5k = median(final.data.trans_temp_ABUN$landcovers.5k),
                       homogen = median(final.data.trans_temp_ABUN$homogen),
                       fert.total_log = median(final.data.trans_temp_ABUN$fert.total_log),
                       percNH = median(final.data.trans_temp_ABUN$percNH),
                       Hansen_mindist_log =  median(final.data.trans_temp_ABUN$Hansen_mindist_log),
                       Forest_biome = "Temperate Broadleaf & Mixed Forests",
                       Use_intensity = "Minimal use",
                       Predominant_land_use = "Primary vegetation",
                       #Tropical = "Temperate",
                       Species_richness = 0,
                       logAbun = 0, 
                       RCAR_110km = 0)


# organise factor levels
# check levels of factor variables
levels(pred_tab$Predominant_land_use) <- levels(abmod_temp$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(abmod_temp$data$Use_intensity) 
levels(pred_tab$Forest_biome) <- levels(abmod_temp$data$Forest_biome)

# add and change factor levels of land use and intensity

pred_tab <- do.call("rbind", replicate(9, pred_tab, simplify = FALSE))


pred_tab[4:6, 'Predominant_land_use'] <- "Secondary vegetation"
pred_tab[7:9, 'Predominant_land_use'] <- "Cropland"


pred_tab[c(2,5,8), 'Use_intensity'] <- "Light use"
pred_tab[c(3,6,9), 'Use_intensity'] <- "Intense use"


### Temperate predictions ###

# predict the result
resulta <- PredictGLMERRandIter(model = abmod_temp$model, data = pred_tab)

# transform the results
resulta <- exp(resulta)-1

resulta <- sweep(x = resulta, MARGIN = 2, STATS = resulta[1,], FUN = '/')

resulta.median <- ((apply(X = resulta, MARGIN = 1, FUN = median))*100)-100
resulta.upper <- ((apply(X = resulta, MARGIN = 1, FUN = quantile,probs = 0.975))*100)-100
resulta.lower <- ((apply(X = resulta, MARGIN = 1, FUN = quantile,probs = 0.025))*100)-100


## plot ##

errbar.cols <- c(rep("#006400",3),rep("#8B0000", 3), rep("#EEAD0E", 3))


pdf(file = paste0(outdir, "/Temperate_Abun_LUUI.pdf"))


errbar(x = 1:9,y = resulta.median,yplus = resulta.upper,yminus = resulta.lower,
       col=errbar.cols,errbar.col = errbar.cols,
       ylim=c(min(resulta.lower),max(resulta.upper)),xaxt="n",
       pch =rep(c(16,17,18), 3), 
       ylab="Total Abundance (%)",xlab="",bty="l", cex.lab = 1.2, cex.axis = 1.2)

axis(side = 1,at = c(2,5,8),
     labels = c("Primary \nvegetation","Secondary\nvegetation", "Cropland"),
     padj = 0.5, cex= 1.2)

abline(h=0,col="#00000077",lty=2)

legend("topleft", 
       legend = c("Minimal Use", "Light Use", "Intense Use"),
       pch = c(16,17,18), bty = "n", inset=c(0,0), cex = 1.2)

dev.off()


#### number of landcovers ####

from = 1
to = 10
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'landcovers.5k'
fac <- NULL
n <- NULL
logval = F


# organise the data
pred_tab <- sort_data(modout = abmod_temp,
                      moddata = final.data.trans_temp_ABUN,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = abmod_temp$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1


# SR plot = full range
ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#CD950C")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#CD950C"), alpha = 0.3) +
  ylim(c(0,250)) +
  xlim(c(0, 11)) +
  xlab("Number of Landcovers") +
  ylab("Total Abundance") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  scale_x_continuous(breaks = c(seq(from = 0, to = 10, by = 2)))

ggsave(filename = paste0(outdir, "/Temperate_Abun_nlandcovers.pdf"))



#### homogeneity ####

from = 0
to = 1
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'homogen'
fac <- NULL
n <- NULL
logval = F


# organise the data
pred_tab <- sort_data(modout = abmod_temp,
                      moddata = final.data.trans_temp_ABUN,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = abmod_temp$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1

homogen <- as.data.frame(unscale(final.data.trans_temp_ABUN$homogen, scale = scalers[4, 2], centre = scalers[4, 3], log = F))

# SR plot = full range
ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#1874CD")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#1874CD"), alpha = 0.3) +
  geom_rug(data = homogen, aes(x = V1), size = 0.1) +
  ylim(c(0,300)) +
  xlim(c(0, 1)) +
  xlab("Homogeneity") +
  ylab("Total Abundance") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) 


ggsave(filename = paste0(outdir, "/Temperate_Abun_homogen.pdf"))



#### percNH ####


from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
fac <- NULL
n <- NULL
logval = F


# organise the data
pred_tab <- sort_data(modout = abmod_temp,
                      moddata = final.data.trans_temp_ABUN,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = abmod_temp$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1

## organise data for plotting ##
percNH <- as.data.frame(unscale(final.data.trans_temp_ABUN$percNH, scale = scalers[5, 2], centre = scalers[5, 3], log = F))

# SR plot = full range
ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#66CD00")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#66CD00"), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(0,350)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Total Abundance") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) 

ggsave(filename = paste0(outdir, "/Temperate_Abun_percNH.pdf"))




#### Fertiliser overall ####

from = 0
to = 2000000
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'fert.total_log'
fac <- NULL
n <- NULL
logval = TRUE


# organise the data
pred_tab <- sort_data(modout = abmod_temp,
                      moddata = final.data.trans_temp_ABUN,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = abmod_temp$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1

# plot
ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#8B6508")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#8B6508"), alpha = 0.3) +
  geom_rug(data = final.data.trans_temp_ABUN, aes(x = fert.total), size = 0.1) +
  ylim(c(0,300)) +
  xlim(c(0, 2000000)) +
  xlab("Total fertiliser application (Kgs)") +
  ylab("Total Abundance") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) 


ggsave(filename = paste0(outdir, "/Temperate_Abun_fert.pdf"))


#### Fertilier, by use intensity ####

from = 0
to = 3000
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'fert.total_log'
fac <- 'Use_intensity'
n <- 3
logval = TRUE

# organise the data
pred_tab <- sort_data(modout = abmod_temp,
                      moddata = final.data.trans_temp_ABUN,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = abmod_temp$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)-1

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]

# plot
 ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = final.data.trans_temp_ABUN, aes(x = fert.total, col = Use_intensity), size = 0.1) +
  ylim(c(0,300)) +
  xlim(c(0, 3000)) +
  xlab("Total fertiliser application (Kgs)") +
  ylab("Total Abundance") +
  scale_colour_manual(values = c("#66CD00", "#FFB90F", "#EE0000"))+
  scale_fill_manual(values = c("#66CD00", "#FFB90F", "#EE0000")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.9), legend.title = element_blank(),
        legend.text = element_text(size = 10),
        aspect.ratio = 1, legend.background = element_blank()) 



 ggsave(filename = paste0(outdir, "/Temperate_Abun_fertUI.pdf"))
 




