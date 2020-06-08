##%######################################################%##
#                                                          #
####         Script for predicting and plotting         ####
#                                                          #
##%######################################################%##

# writing a script to predict and then plot values from models

rm(list = ls())
setwd("C:/Users/charl/Dropbox/POSTDOC - BIOTA/0. PROJECTS/BIOTA")

# load libraries
library(StatisticalModels)

# point to directories
datadir <- '1_PREDICTS_PLUS_VARIABLES'
moddir <- '2_MODEL_SELECTION'
outdir <- '3_PLOTTING'

### load the datasets and the models ###

# load the selected models
load(paste0(moddir, "/SRMOD_output.rdata")) # srmod
load(paste0(moddir, "/ABMOD_output.rdata")) # abmod
load(paste0(moddir, "/RCAR_output.rdata")) # rcarmod

# load in the datasets with transformed variables and the subsets 
load(paste0(datadir, "/PREDICTS_dataset_inc_variables_TRANS.rdata")) # final.data.trans
load(paste0(datadir, "/PREDICTS_abun_subset.rdata")) # final.data.abun
load(paste0(datadir, "/PREDICTS_rcar_subset.rdata")) # final.data.rcar

# read in the values used for rescaling
scalers <- read.csv("1_PREDICTS_PLUS_VARIABLES/Scaling_values.csv")

# load the function for sorting data for predictions
source('functions/sort_data.r')
source('functions/rescale.r')
source('functions/unscale.r')


#### 1. Species Richness results ####

### Factor level plots ###
pdf(file = paste0(outdir, "/SR_factor_plots.pdf"), paper = "a4r", 
    title = "Species richness model: all effects")

par(mfrow = c(2,2))


### factor level plots ###

# Land use
PlotGLMERFactor(model = srmod$model,data = srmod$data,
                responseVar = "Species richness",seMultiplier = 1.96,
                logLink = "e",catEffects = c("Predominant_land_use"), params = list(las = 2, cex = 0.8, mar = c(6,3,2,4), adj = 0.5), xtext.srt = 90)



# Use intensity
PlotGLMERFactor(model = srmod$model,data = srmod$data,
                responseVar = "Species richness",seMultiplier = 1.96,
                logLink = "e",catEffects = c( "Use_intensity"), params = list(las = 2, cex = 0.8, mar = c(6,3,2,4), adj = 0.5), xtext.srt = 90)


# forest biomes
PlotGLMERFactor(model = srmod$model,data = srmod$data,
                responseVar = "Species richness",seMultiplier = 1.96,
                logLink = "e",catEffects = "Forest_biome", params = list(las = 2, cex = 0.8, mar = c(9,5,1,4), adj = 1), xtext.srt = 50, ylim = c(-70, 100))



# Tropical
PlotGLMERFactor(model = srmod$model,data = srmod$data,
                responseVar = "Species richness",seMultiplier = 1.96,
                logLink = "e",catEffects = c("Tropical"), params = list(las = 2, cex = 0.8, mar = c(4,4,4,2), adj = 0.5), xtext.srt = 90)


dev.off()




### land use intensity interactions ###


# basic table of median values and reference factors
pred_tab <- data.frame(landcovers.5k = median(final.data.trans$landcovers.5k),
                       homogen = median(final.data.trans$homogen),
                       fert.total_log = median(final.data.trans$fert.total_log),
                       percNH = median(final.data.trans$percNH),
                       Hansen_mindist_log =  median(final.data.trans$Hansen_mindist_log),
                       Forest_biome = "Temperate Broadleaf & Mixed Forests",
                       Use_intensity = "Minimal use",
                       Predominant_land_use = "Primary vegetation",
                       Tropical = "Temperate",
                       Species_richness = 0,
                       logAbun = 0, 
                       RCAR_110km = 0)


# organise factor levels
# check levels of factor variables
levels(pred_tab$Predominant_land_use) <- levels(srmod$data$Predominant_land_use)
levels(pred_tab$Use_intensity) <- levels(srmod$data$Use_intensity) 
levels(pred_tab$Forest_biome) <- levels(srmod$data$Forest_biome) 
levels(pred_tab$Tropical) <- levels(srmod$data$Tropical) 

# add and change factor levels of land use and intensity

pred_tab <- do.call("rbind", replicate(9, pred_tab, simplify = FALSE))


pred_tab[4:6, 'Predominant_land_use'] <- "Secondary vegetation"
pred_tab[7:9, 'Predominant_land_use'] <- "Cropland"



pred_tab[c(2,5,8), 'Use_intensity'] <- "Light use"
pred_tab[c(3,6,9), 'Use_intensity'] <- "Intense use"


#### SR predictions ####


# predict the result
result <- PredictGLMERRandIter(model = srmod$model, data = pred_tab)

# transform the results
result <- exp(result)

result <- sweep(x = result,MARGIN = 2,STATS = result[1,],FUN = '/')

result.median <- ((apply(X = result,MARGIN = 1,FUN = median))*100)-100
result.upper <- ((apply(X = result,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
result.lower <- ((apply(X = result,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100


#### Abundance predictions ####

# predict the result
resulta <- PredictGLMERRandIter(model = abmod$model, data = pred_tab)

# transform the results
resulta <- exp(resulta)-1

resulta <- sweep(x = resulta, MARGIN = 2, STATS = resulta[1,], FUN = '/')

resulta.median <- ((apply(X = resulta, MARGIN = 1, FUN = median))*100)-100
resulta.upper <- ((apply(X = resulta, MARGIN = 1, FUN = quantile,probs = 0.975))*100)-100
resulta.lower <- ((apply(X = resulta, MARGIN = 1, FUN = quantile,probs = 0.025))*100)-100


#### RCAR predictions ####

# predict the result
resultr <- PredictGLMERRandIter(model = rcarmod$model, data = pred_tab)

# transform the results
resultr <- 10^(resultr)

resultr <- sweep(x = resultr, MARGIN = 2, STATS = resultr[1,], FUN = '/')

resultr.median <- ((apply(X = resultr, MARGIN = 1, FUN = median))*100)-100
resultr.upper <- ((apply(X = resultr, MARGIN = 1, FUN = quantile,probs = 0.975))*100)-100
resultr.lower <- ((apply(X = resultr, MARGIN = 1, FUN = quantile,probs = 0.025))*100)-100





# create the plots


pdf(file = paste0(outdir, "/LU_UI_plots.pdf"))
    

par(mfrow = c(2,2), mar = c(5, 4, 1, 2))


## SR
# colours
errbar.cols <- c(rep("#006400",3),rep("#8B0000", 3), rep("#EEAD0E", 3))


errbar(x = 1:9,y = result.median,yplus = result.upper,yminus = result.lower,
       col=errbar.cols,errbar.col = errbar.cols,
       ylim=c(min(result.lower),18),xaxt="n",
       pch =rep(c(16,17,18), 3), 
       ylab="Species Richness (%)",xlab="",bty="l")

axis(side = 1,at = c(2,5,8),
     labels = c("Primary \nvegetation","Secondary\nvegetation", "Cropland"),
     padj = 0.5)

abline(h=0,col="#00000077",lty=2)

#legend("topright", 
#       legend = c("Minimal Use", "Light Use", "Intense Use"),
#       pch = c(16,17,18), bty = "n", inset=c(0,0))


## Ab
errbar(x = 1:9,y = resulta.median,yplus = resulta.upper,yminus = resulta.lower,
       col=errbar.cols,errbar.col = errbar.cols,
       ylim=c(min(resulta.lower),38),xaxt="n",
       pch =rep(c(16,17,18), 3), 
       ylab="Total Abundance (%)",xlab="",bty="l")

axis(side = 1,at = c(2,5,8),
     labels = c("Primary \nvegetation","Secondary\nvegetation", "Cropland"),
     padj = 0.5)

abline(h=0,col="#00000077",lty=2)

#legend("topright", 
#       legend = c("Minimal Use", "Light Use", "Intense Use"),
#       pch = c(16,17,18), bty = "n", inset=c(0,0)
#)


## RCAR
errbar(x = 1:9,y = resultr.median,yplus = resultr.upper,yminus = resultr.lower,
       col=errbar.cols,errbar.col = errbar.cols,
       ylim=c(min(resultr.lower),116),xaxt="n",
       pch =rep(c(16,17,18), 3), 
       ylab="RCAR (%)",xlab="",bty="l")

axis(side = 1,at = c(2,5,8),
     labels = c("Primary \nvegetation","Secondary\nvegetation", "Cropland"),
     padj = 0.5)

abline(h=0,col="#00000077",lty=2)

legend("topleft", 
       legend = c("Minimal Use", "Light Use", "Intense Use"),
       pch = c(16,17,18), bty = "n", inset=c(0,0)
)

dev.off()



### Continous results and interactions ###

# list to save plots in
p <- list()



### Distance to forest ###

from = 0
to = 350
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'Hansen_mindist_log'
fac <- NULL
n <- NULL
logval = TRUE


# organise the data
pred_tab <- sort_data(modout = srmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n, 
                      log = logval)

# predict the result
result <- PredictGLMER(model = srmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##

# add the new vals
#result$vals <- rep(vals, n)
#result$factor <- pred_tab3[, fac]

# SR plot = full range
p[[1]] <- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#458B00")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#458B00"), alpha = 0.3) +
  geom_rug(data = final.data.trans, aes(x = Hansen_mindist), size = 0.1) +
  ylim(c(0,30)) +
  xlim(c(0, 350)) +
  xlab("Distance to Forest (Km)") +
  ylab("Species Richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("A.")
  

p[[1]] <- ggplotGrob(p[[1]])



### percNH ###


from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
fac <- NULL
n <- NULL
logval = F


# organise the data
pred_tab <- sort_data(modout = srmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = srmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##
percNH <- as.data.frame(unscale(final.data.trans$percNH, scale = scalers[5, 2], centre = scalers[5, 3], log = F))

# SR plot = full range
p[[2]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#66CD00")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#66CD00"), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(0,30)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Species Richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("A.")

p[[2]] <- ggplotGrob(p[[2]])


### Fertiliser overall ###

from = 0
to = 2500000
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'fert.total_log'
fac <- NULL
n <- NULL
logval = TRUE


# organise the data
pred_tab <- sort_data(modout = srmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = srmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

# SR plot = full range
p[[3]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#8B6508")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#8B6508"), alpha = 0.3) +
  geom_rug(data = final.data.trans, aes(x = fert.total), size = 0.1) +
  ylim(c(0,30)) +
  xlim(c(0, 2500000)) +
  xlab("Total fertiliser application (Kgs)") +
  ylab("Species Richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) 

p[[3]] <- ggplotGrob(p[[3]])


### number of landcovers ###

from = 1
to = 10
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'landcovers.5k'
fac <- NULL
n <- NULL
logval = F


# organise the data
pred_tab <- sort_data(modout = srmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = srmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)


# SR plot = full range
p[[4]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#CD950C")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#CD950C"), alpha = 0.3) +
  #geom_rug(data = final.data.trans, aes(x = landcovers.5k)) +
  ylim(c(5,30)) +
  xlim(c(0, 11)) +
  xlab("Number of Landcovers") +
  ylab("Species Richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) 

p[[4]] <- ggplotGrob(p[[4]])


### homogeneity ###

from = 0
to = 1
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'homogen'
fac <- NULL
n <- NULL
logval = F


# organise the data
pred_tab <- sort_data(modout = srmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = srmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

homogen <- as.data.frame(unscale(final.data.trans$homogen, scale = scalers[4, 2], centre = scalers[4, 3], log = F))

# SR plot = full range
p[[5]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#1874CD")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#1874CD"), alpha = 0.3) +
  geom_rug(data = homogen, aes(x = V1), size = 0.1) +
  ylim(c(5,30)) +
  xlim(c(0, 1)) +
  xlab("Homogeneity") +
  ylab("Species Richness") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) 

p[[5]] <- ggplotGrob(p[[5]])






### Fertilier, by land use ###

# in this instance, predicting fertiliser vals 
from = 0
to = 3000
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'fert.total_log'
fac <- 'Predominant_land_use'
n <- 3
logval = TRUE

# organise the data
pred_tab <- sort_data(modout = srmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = srmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]

# SR plot = full range
p[[7]] <- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = final.data.trans, aes(x = fert.total, col = Predominant_land_use), size = 0.1) +
  ylim(c(0,30)) +
  xlim(c(0, 3000)) +
  xlab("Total fertiliser application (Kgs)") +
  ylab("Species Richness") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("A.")

p[[7]] <- ggplotGrob(p[[7]])




### Homogen, by land use ###

# in this instance, predicting fertiliser vals 
from = 0
to = 1
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'homogen'
fac <- 'Predominant_land_use'
n <- 3
logval = F

# organise the data
pred_tab <- sort_data(modout = srmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = srmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]

homogen <- as.data.frame(unscale(final.data.trans$homogen, scale = scalers[4, 2], centre = scalers[4, 3], log = F))


# SR plot = full range
p[[8]] <- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = homogen, aes(x = V1), size = 0.1) +
  ylim(c(0,30)) +
  xlim(c(0, 1)) +
  xlab("Homogeneity") +
  ylab("Species Richness") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("A.")

p[[8]] <- ggplotGrob(p[[8]])





### landcovers, by land use ###

# in this instance, predicting fertiliser vals 
from = 0
to = 10
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'landcovers.5k'
fac <- 'Predominant_land_use'
n <- 3
logval = F

# organise the data
pred_tab <- sort_data(modout = srmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = srmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]

landcovers.5k <- as.data.frame(unscale(final.data.trans$landcovers.5k, scale = scalers[3, 2], centre = scalers[3, 3], logval = F))


# SR plot = full range
p[[9]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  #geom_rug(data = landcovers.5k, aes(x = V1), size = 0.1) +
  ylim(c(0,30)) +
  xlim(c(0, 10)) +
  xlab("Number of Landcovers") +
  ylab("Species Richness") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("A.") +
  scale_x_continuous(breaks = c(seq(from = 0, to = 10, by = 2)))

p[[9]] <- ggplotGrob(p[[9]])






### Fertilier, byuse intensity ###

# in this instance, predicting fertiliser vals 
from = 0
to = 3000
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'fert.total_log'
fac <- 'Use_intensity'
n <- 3
logval = TRUE

# organise the data
pred_tab <- sort_data(modout = srmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = srmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]

# SR plot = full range
p[[10]] <- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = final.data.trans, aes(x = fert.total, col = Use_intensity), size = 0.1) +
  ylim(c(0,30)) +
  xlim(c(0, 3000)) +
  xlab("Total fertiliser application (Kgs)") +
  ylab("Species Richness") +
  scale_colour_manual(values = c("#66CD00", "#FFB90F", "#EE0000"))+
  scale_fill_manual(values = c("#66CD00", "#FFB90F", "#EE0000")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.9), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1) 

p[[10]] <- ggplotGrob(p[[10]])





### percNH, byuse intensity ###

# in this instance, predicting fertiliser vals 
from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
fac <- 'Use_intensity'
n <- 3
logval = F

# organise the data
pred_tab <- sort_data(modout = srmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = srmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]

percNH <- as.data.frame(unscale(final.data.trans$percNH, scale = scalers[5, 2], centre = scalers[5, 3], log = F))


# SR plot = full range
p[[11]] <- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(0,30)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Species Richness") +
  scale_colour_manual(values = c("#66CD00", "#FFB90F", "#EE0000"))+
  scale_fill_manual(values = c("#66CD00", "#FFB90F", "#EE0000")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.9), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1) 

p[[11]] <- ggplotGrob(p[[11]])



### distance, by tropical ###

from = 0
to = 20
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'Hansen_mindist_log'
fac <- 'Tropical'
n <- 2
logval = T

# organise the data
pred_tab <- sort_data(modout = srmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = srmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


# SR plot = full range
p[[12]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = final.data.trans, aes(x = Hansen_mindist_log, col = Tropical), size = 0.1) +
  ylim(c(0,30)) +
  xlim(c(0, 20)) +
  xlab("Distance to forest (Km)") +
  ylab("Species Richness") +
  scale_colour_manual(values = c("#009ACD", "#9932CC"))+
  scale_fill_manual(values = c("#009ACD", "#9932CC")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none', legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1) 

p[[12]] <- ggplotGrob(p[[12]])


# organise species richness plots

library(gridExtra)


lay <-rbind(c(1,2),
            c(3,4),
            c(5,NA),
            c(7,8),
            c(9,10),
            c(11,12))

plots <- marrangeGrob(grobs = p[c(1:5, 7:12)], npages = 2, ncol = 2, nrow = 4,  layout_matrix =lay, top = "")

ggsave(filename = paste0(outdir, "/SR_cont_plots.pdf"), plots, height = 16, width = 8)




#### 2. Total Abundance results ####

### Factor level plots ###

pdf(file = paste0(outdir, "/AB_factor_plots.pdf"), paper = "a4r")

par(mfrow = c(2,2))


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


dev.off()





### Continous results and interactions ###


### Distance to forest ###

from = 0
to = 350
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'Hansen_mindist_log'
fac <- NULL
n <- NULL
logval = TRUE


# organise the data
pred_tab <- sort_data(modout = abmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n, 
                      log = logval)

# predict the result
result <- PredictGLMER(model = abmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##

q <-list()

q[[1]] <- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#458B00")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#458B00"), alpha = 0.3) +
  geom_rug(data = final.data.abun, aes(x = Hansen_mindist), size = 0.1) +
  ylim(c(0, 300)) +
  xlim(c(0, 350)) +
  xlab("Distance to Forest (Km)") +
  ylab("Total Abundance") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("B.")

q[[1]] <- ggplotGrob(q[[1]])



### percNH ###


from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
fac <- NULL
n <- NULL
logval = F


# organise the data
pred_tab <- sort_data(modout = abmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = abmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##
percNH <- as.data.frame(unscale(final.data.abun$percNH, scale = scalers[5, 2], centre = scalers[5, 3], log = F))

# SR plot = full range
q[[2]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#66CD00")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#66CD00"), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(0,300)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Total Abundance") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) 

q[[2]] <- ggplotGrob(q[[2]])


### Fertiliser overall ###

from = 0
to = 2500000
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'fert.total_log'
fac <- NULL
n <- NULL
logval = TRUE


# organise the data
pred_tab <- sort_data(modout = abmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = abmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

# SR plot = full range
q[[3]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#8B6508")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#8B6508"), alpha = 0.3) +
  geom_rug(data = final.data.abun, aes(x = fert.total), size = 0.1) +
  ylim(c(0,300)) +
  xlim(c(0, 2500000)) +
  xlab("Total fertiliser application (Kgs)") +
  ylab("Total Abundance") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) 

q[[3]] <- ggplotGrob(q[[3]])


### number of landcovers ###

from = 1
to = 10
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'landcovers.5k'
fac <- NULL
n <- NULL
logval = F


# organise the data
pred_tab <- sort_data(modout = abmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = abmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)


# SR plot = full range
q[[4]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#CD950C")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#CD950C"), alpha = 0.3) +
  ylim(c(0,300)) +
  xlim(c(0, 11)) +
  xlab("Number of Landcovers") +
  ylab("Total Abundance") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("B.") +
  scale_x_continuous(breaks = c(seq(from = 0, to = 10, by = 2)))


q[[4]] <- ggplotGrob(q[[4]])


### homogeneity ###

from = 0
to = 1
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'homogen'
fac <- NULL
n <- NULL
logval = F


# organise the data
pred_tab <- sort_data(modout = abmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = abmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

homogen <- as.data.frame(unscale(final.data.abun$homogen, scale = scalers[4, 2], centre = scalers[4, 3], log = F))

# SR plot = full range
q[[5]]<- ggplot(data = result) +
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
        aspect.ratio = 1) +
  ggtitle("B.")

q[[5]] <- ggplotGrob(q[[5]])




### Fertilier, by land use ###

# in this instance, predicting fertiliser vals 
from = 0
to = 3000
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'fert.total_log'
fac <- 'Predominant_land_use'
n <- 3
logval = TRUE

# organise the data
pred_tab <- sort_data(modout = abmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = abmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]

# SR plot = full range
q[[7]] <- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = final.data.abun, aes(x = fert.total, col = Predominant_land_use), size = 0.1) +
  ylim(c(0,300)) +
  xlim(c(0, 3000)) +
  xlab("Total fertiliser application (Kgs)") +
  ylab("Total Abundance") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none", legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("B.")

q[[7]] <- ggplotGrob(q[[7]])




### percNH, by land use ###

# in this instance, predicting fertiliser vals 
from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
fac <- 'Predominant_land_use'
n <- 3
logval = F

# organise the data
pred_tab <- sort_data(modout = abmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = abmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]

percNH <- as.data.frame(unscale(final.data.trans$percNH, scale = scalers[5, 2], centre = scalers[5, 3], log = F))


# SR plot = full range
q[[8]] <- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(0,300)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Total Abundance") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.8), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("B.")

q[[8]] <- ggplotGrob(q[[8]])





### Fertilier, by use intensity ###

# in this instance, predicting fertiliser vals 
from = 0
to = 3000
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'fert.total_log'
fac <- 'Use_intensity'
n <- 3
logval = TRUE

# organise the data
pred_tab <- sort_data(modout = abmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = abmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]

# SR plot = full range
q[[9]] <- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = final.data.abun, aes(x = fert.total, col = Use_intensity), size = 0.1) +
  ylim(c(0,300)) +
  xlim(c(0, 3000)) +
  xlab("Total fertiliser application (Kgs)") +
  ylab("Total Abundance") +
  scale_colour_manual(values = c("#66CD00", "#FFB90F", "#EE0000"))+
  scale_fill_manual(values = c("#66CD00", "#FFB90F", "#EE0000")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.9), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1) 

q[[9]] <- ggplotGrob(q[[9]])





### percNH, by use intensity ###

# in this instance, predicting fertiliser vals 
from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
fac <- 'Use_intensity'
n <- 3
logval = F

# organise the data
pred_tab <- sort_data(modout = abmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = abmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]

percNH <- as.data.frame(unscale(final.data.abun$percNH, scale = scalers[5, 2], centre = scalers[5, 3], log = F))


# SR plot = full range
q[[10]] <- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(0,300)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Total Abundance") +
  scale_colour_manual(values = c("#66CD00", "#FFB90F", "#EE0000"))+
  scale_fill_manual(values = c("#66CD00", "#FFB90F", "#EE0000")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.9), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1) 

q[[10]] <- ggplotGrob(q[[10]])



### distance, by tropical ###

from = 0
to = 20
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'Hansen_mindist_log'
fac <- 'Tropical'
n <- 2
logval = T

# organise the data
pred_tab <- sort_data(modout = abmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = abmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


# SR plot = full range
q[[11]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = final.data.trans, aes(x = Hansen_mindist_log, col = Tropical), size = 0.1) +
  ylim(c(0,300)) +
  xlim(c(0, 20)) +
  xlab("Distance to forest (Km)") +
  ylab("Total Abundance") +
  scale_colour_manual(values = c("#009ACD", "#9932CC"))+
  scale_fill_manual(values = c("#009ACD", "#9932CC")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.9), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) 

q[[11]] <- ggplotGrob(q[[11]])


# organise species richness plots

lay <-rbind(c(1,2),
            c(3,4),
            c(5,NA),
            c(7,8),
            c(9,10),
            c(11,NA))

plots <- marrangeGrob(grobs = q[c(1:5, 7:11)], npages = 2, ncol = 2, nrow = 4,  layout_matrix =lay, top = "")

ggsave(filename = paste0(outdir, "/AB_cont_plots.pdf"), plots, height = 16, width = 8)




### 3. RCAR Results ###

### Factor level plots ###

pdf(file = paste0(outdir, "/RCAR_factor_plots.pdf"), paper = "a4r")

par(mfrow = c(2,2))


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


dev.off()






### Continous results and interactions ###

# list to save plots in
r <- list()



### Distance to forest ###

from = 0
to = 350
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'Hansen_mindist_log'
fac <- NULL
n <- NULL
logval = TRUE


# organise the data
pred_tab <- sort_data(modout = rcarmod,
                      moddata = final.data.rcar,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n, 
                      log = logval)

# predict the result
result <- PredictGLMER(model = rcarmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- 10^(result)

## organise data for plotting ##

r[[1]] <- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#458B00")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#458B00"), alpha = 0.3) +
  geom_rug(data = final.data.rcar, aes(x = Hansen_mindist), size = 0.1) +
  ylim(c(100000, 3500000)) +
  xlim(c(0, 350)) +
  xlab("Distance to Forest (Km)") +
  ylab("RCAR (square kilometres)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) +
  ggtitle("B.")

r[[1]] <- ggplotGrob(r[[1]])



### percNH ###


from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
fac <- NULL
n <- NULL
logval = F


# organise the data
pred_tab <- sort_data(modout = rcarmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = rcarmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- 10^(result)

## organise data for plotting ##
percNH <- as.data.frame(unscale(final.data.rcar$percNH, scale = scalers[5, 2], centre = scalers[5, 3], log = F))

# SR plot = full range
r[[2]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#66CD00")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#66CD00"), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  ylim(c(15000, 4000000)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("RCAR (square kilometres") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) 

r[[2]] <- ggplotGrob(r[[2]])



### number of landcovers ###

from = 1
to = 10
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'landcovers.5k'
fac <- NULL
n <- NULL
logval = F


# organise the data
pred_tab <- sort_data(modout = rcarmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = rcarmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- 10^(result)


# SR plot = full range
r[[3]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y), col = c("#CD950C")) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus), fill = c("#CD950C"), alpha = 0.3) +
  ylim(c(15000, 4000000)) +
  xlim(c(0, 11)) +
  xlab("Number of Landcovers") +
  ylab("RCAR (square kilometres") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.8), legend.title = element_blank(),
        aspect.ratio = 1) 

r[[3]] <- ggplotGrob(r[[3]])


### homogeneity ###

from = 0
to = 1
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'homogen'
fac <- NULL
n <- NULL
logval = F


# organise the data
pred_tab <- sort_data(modout = rcarmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = rcarmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- 10^(result)

homogen <- as.data.frame(unscale(final.data.rcar$homogen, scale = scalers[4, 2], centre = scalers[4, 3], log = F))

# SR plot = full range
r[[4]]<- ggplot(data = result) +
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

r[[4]] <- ggplotGrob(r[[4]])







### percNH, by land use ###

# in this instance, predicting fertiliser vals 
from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
fac <- 'Predominant_land_use'
n <- 3
logval = F

# organise the data
pred_tab <- sort_data(modout = rcarmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = rcarmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- 10^(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]

percNH <- as.data.frame(unscale(final.data.rcar$percNH, scale = scalers[5, 2], centre = scalers[5, 3], log = F))


# SR plot = full range
r[[6]] <- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
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
        aspect.ratio = 1) +
  ggtitle("C.")

r[[6]] <- ggplotGrob(r[[6]])



### Homogen, by land use ###

# in this instance, predicting fertiliser vals 
from = 0
to = 1
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'homogen'
fac <- 'Predominant_land_use'
n <- 3
logval = F

# organise the data
pred_tab <- sort_data(modout = rcarmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = rcarmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- 10^(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]

homogen <- as.data.frame(unscale(final.data.rcar$homogen, scale = scalers[4, 2], centre = scalers[4, 3], log = F))


r[[7]] <- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = homogen, aes(x = V1), size = 0.1) +
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
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("C.")

r[[7]] <- ggplotGrob(r[[7]])






### distance, by use intensity ###

# in this instance, predicting fertiliser vals 
from = 0
to = 20
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'Hansen_mindist_log'
fac <- 'Use_intensity'
n <- 3
logval =TRUE

# organise the data
pred_tab <- sort_data(modout = rcarmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = rcarmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- 10^(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


# SR plot = full range
r[[8]] <- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = final.data.rcar, aes(x = Hansen_mindist), size = 0.1) +
  ylim(c(15000, 4000000)) +
  xlim(c(0, 20)) +
  xlab("Distance to forest (Km)") +
  ylab("RCAR (square kilometres)") +
  scale_colour_manual(values = c("#66CD00", "#FFB90F", "#EE0000"))+
  scale_fill_manual(values = c("#66CD00", "#FFB90F", "#EE0000")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.9), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1) 

r[[8]] <- ggplotGrob(r[[8]])




# land covers by Use intensity


from = 0
to = 10
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'landcovers.5k'
fac <- 'Use_intensity'
n <- 3
logval = F

# organise the data
pred_tab <- sort_data(modout = rcarmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = rcarmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- 10^(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]


# SR plot = full range
r[[9]]<- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  #geom_rug(data = landcovers.5k, aes(x = V1)) +
  ylim(c(15000, 4000000)) +
  xlim(c(0, 10)) +
  xlab("Number of Landcovers") +
  ylab("RCAR (square kilometres)") +
  scale_colour_manual(values = c("#66CD00", "#FFB90F", "#EE0000"))+
  scale_fill_manual(values = c("#66CD00", "#FFB90F", "#EE0000")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.9), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1) 

r[[9]] <- ggplotGrob(r[[9]])







# homogen by use intensity


# in this instance, predicting fertiliser vals 
from = 0
to = 1
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'homogen'
fac <- 'Use_intensity'
n <- 3
logval = F

# organise the data
pred_tab <- sort_data(modout = rcarmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = rcarmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- 10^(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]

homogen <- as.data.frame(unscale(final.data.rcar$homogen, scale = scalers[4, 2], centre = scalers[4, 3], log = F))


# SR plot = full range
r[[10]] <- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = homogen, aes(x = V1), size = 0.1) +
  ylim(c(15000, 4000000)) +
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

r[[10]] <- ggplotGrob(r[[10]])



# organise species richness plots

lay <-rbind(c(1,2),
            c(3,4),
            c(NA,6),
            c(7,8),
            c(9,10))

plots <- marrangeGrob(grobs = r[c(1:4, 6:10)], npages = 2, ncol = 2, nrow = 4,  layout_matrix =lay, top = "")

ggsave(filename = paste0(outdir, "/RCAR_cont_plots.pdf"), plots, height = 16, width = 8)






##%######################################################%##
#                                                          #
####                plots for manuscript                ####
#                                                          #
##%######################################################%##



# distance (SR and RCAR)

dist_plots <- marrangeGrob(grobs = list(p[[1]], r[[1]]), npages = 1, ncol = 2, nrow = 1, top = "")

ggsave(filename = paste0(outdir, "/Paper_dist_plots.pdf"), dist_plots, height = 4, width = 8)


# percentage natural habitat (SR, Abun interaction, RCAR interaction)

percNH_plots <- marrangeGrob(grobs = list(p[[2]], r[[6]], q[[8]]), npages = 1, ncol = 2, nrow = 2, top = "")


ggsave(filename = paste0(outdir, "/Paper_percNH_plots.pdf"), percNH_plots, height = 8, width = 8)


# homogen (SR interaction, Abun, RCAR interaction)

homogen_plots <- marrangeGrob(grobs = list(p[[8]], r[[7]], q[[5]]), npages = 1, ncol = 2, nrow = 2, top = "")

ggsave(filename = paste0(outdir, "/Paper_homogen_plots.pdf"), homogen_plots, height = 8, width = 8)


# landcovers (SR interaction, Abun)

landcovers_plots <- marrangeGrob(grobs = list(p[[9]], q[[4]]), npages = 1, ncol = 2, nrow = 1, top = "")

ggsave(filename = paste0(outdir, "/Paper_landcovers_plots.pdf"), landcovers_plots, height = 4, width = 8)


# fertiliser (SR interaction, Abun interaction)

fert_plots <- marrangeGrob(grobs = list(p[[7]], q[[7]]), npages = 1, ncol = 2, nrow = 1, top = "")

ggsave(filename = paste0(outdir, "/Paper_fert_plots.pdf"), fert_plots, height = 4, width = 8)

# tropical distance interaction

trop_plots <- marrangeGrob(grobs = list(p[[12]], q[[11]]), npages = 1, ncol = 2, nrow = 1, top = "")

ggsave(filename = paste0(outdir, "/Paper_tropical_plots.pdf"), trop_plots, height = 4, width = 8)


# fertiliser full range

### Fertilier, by land use ###

# in this instance, predicting fertiliser vals 
from = 0
to = 800000
vals <- seq(from = from, to = to, length.out = to/30)
variable <- 'fert.total_log'
fac <- 'Predominant_land_use'
n <- 3
logval = TRUE

# organise the data
pred_tab <- sort_data(modout = srmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result1 <- PredictGLMER(model = srmod$model, data = pred_tab[pred_tab$Predominant_land_use == "Primary vegetation",], se.fit = TRUE, seMultiplier = 1.96)
result2 <- PredictGLMER(model = srmod$model, data = pred_tab[pred_tab$Predominant_land_use == "Secondary vegetation",], se.fit = TRUE, seMultiplier = 1.96)
result3 <- PredictGLMER(model = srmod$model, data = pred_tab[pred_tab$Predominant_land_use == "Cropland",], se.fit = TRUE, seMultiplier = 1.96)

result <- PredictGLMER(model = srmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)


# transform the results
result <- exp(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]

# SR plot = full range
p[[13]] <- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = final.data.trans, aes(x = fert.total, col = Predominant_land_use), size = 0.1) +
  ylim(c(0,30)) +
  xlim(c(0, 800000)) +
  xlab("Total fertiliser application (Kgs)") +
  ylab("Species Richness") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.2,0.85), legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("A.")

p[[7]] <- ggplotGrob(p[[7]])



### Fertilier, by land use ###

# in this instance, predicting fertiliser vals 
from = 0
to = 3000
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'fert.total_log'
fac <- 'Predominant_land_use'
n <- 3
logval = TRUE

# organise the data
pred_tab <- sort_data(modout = abmod,
                      moddata = final.data.trans,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# predict the result
result <- PredictGLMER(model = abmod$model, data = pred_tab, se.fit = TRUE, seMultiplier = 1.96)

# transform the results
result <- exp(result)

## organise data for plotting ##

# add the new vals
result$vals <- rep(vals, n)
result$factor <- pred_tab[, fac]

# SR plot = full range
q[[7]] <- ggplot(data = result) +
  geom_line(aes(x = vals, y = y, col = factor)) +
  geom_ribbon(aes(x = vals, ymin= yminus, ymax = yplus, fill = factor), alpha = 0.3) +
  geom_rug(data = final.data.abun, aes(x = fert.total, col = Predominant_land_use), size = 0.1) +
  ylim(c(0,300)) +
  xlim(c(0, 3000)) +
  xlab("Total fertiliser application (Kgs)") +
  ylab("Total Abundance") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none", legend.title = element_blank(),
        legend.text = element_text(size = 6),
        aspect.ratio = 1, legend.background = element_blank()) +
  ggtitle("B.")

q[[7]] <- ggplotGrob(q[[7]])

