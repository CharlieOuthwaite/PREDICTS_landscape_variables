##%######################################################%##
#                                                          #
##%######################################################%##
#                                                          #
####            Plotting using GLMERRandIter            ####
####              function for predictions              ####
#                                                          #
##%######################################################%##

#                                                          #
##%######################################################%##


### NOTE: where the plot is not split by land use, the baseline LU is cropland
### NOTE: where the plot is not split by use intensity, the baseline UI is Intense use.
### NOTE: final models are models without group level random effects



rm(list = ls())

# load libraries
library(StatisticalModels)
library(gridExtra)
library(cowplot)

# directories
datadir <- '1_PREDICTS_PLUS_VARIABLES'
moddir <- '2_MODEL_SELECTION'
outdir <- '3_PLOTTING'

# create outdir if it doesn't exist already
if(!dir.exists(outdir)) dir.create(outdir)

# load the selected models
load(paste0(moddir, "/ABMOD_Tropical_output.rdata"))
load(paste0(moddir, "/ABMOD_Temperate_output.rdata"))
load(paste0(moddir, "/SRMOD_Tropical_output.rdata"))
load(paste0(moddir, "/SRMOD_Temperate_output.rdata"))

# read in the values used for rescaling
scalers <- read.csv("1_PREDICTS_PLUS_VARIABLES/Scaling_values.csv")

# load the function for sorting data for predictions
source('functions/sort_data.r')
source('functions/rescale.r')
source('functions/unscale.r')


# load the datasets
load(paste0(moddir, "/final.data.trans_trop_ABUN.rdata"))
load(paste0(moddir, "/final.data.trans_temp_ABUN.rdata"))

# load the dataset
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

# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)

#### create custom theme to be applied to all plots ####

theme_custom <- theme(panel.grid = element_blank(),
                      legend.position = c(0.8,0.8), 
                      legend.title = element_blank(),
                      legend.text = element_text(size = 8),
                      aspect.ratio = 1, legend.background = element_blank(),
                      #text = element_text(size = 8), 
                      axis.text = element_text(size = 8), 
                      line = element_line(size = 0.2), 
                      panel.border = element_rect(size = 0.2),
                      axis.title = element_text(size = 8))


##%######################################################%##
#                                                          #
####                 Distance to forest                 ####
#                                                          #
##%######################################################%##

#### Abundance, Tropical ####


range(final.data.trans$Hansen_mindist)
# 0.0000 324.4697

from = 0
to = 325
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


# reference for % difference = primary vegetation and distance closest to 0
refRow <- which((pred_tab$Predominant_land_use=="Cropland") & (pred_tab$Hansen_mindist_logRS==min(pred_tab$Hansen_mindist_logRS)))


# sort quantiles
QPV <- quantile(x = abmod_trop$data$Hansen_mindist_logRS[
  abmod_trop$data$Predominant_land_use=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = abmod_trop$data$Hansen_mindist_logRS[
  abmod_trop$data$Predominant_land_use=="Secondary vegetation"],
  probs = exclQuantiles)
QCR <- quantile(x = abmod_trop$data$Hansen_mindist_logRS[
  abmod_trop$data$Predominant_land_use=="Cropland"],
  probs = exclQuantiles)

# predict the results
result <- PredictGLMERRandIter(model = abmod_trop$model,data = pred_tab, nIters = 10000)

# back transform the abundance values
result <- exp(result)-1

# convert to relative to reference
result <- sweep(x = result,MARGIN = 2,STATS = result[refRow,],FUN = '/')

# remove anything above and below the quantiles
result[which(pred_tab$Predominant_land_use == "Primary vegetation" & pred_tab$Hansen_mindist_logRS < QPV[1]), ] <- NA
result[which(pred_tab$Predominant_land_use == "Primary vegetation" & pred_tab$Hansen_mindist_logRS > QPV[2]), ] <- NA
result[which(pred_tab$Predominant_land_use == "Secondary vegetation" & pred_tab$Hansen_mindist_logRS < QSV[1]), ] <- NA
result[which(pred_tab$Predominant_land_use == "Secondary vegetation" & pred_tab$Hansen_mindist_logRS > QSV[2]), ] <- NA
result[which(pred_tab$Predominant_land_use == "Cropland" & pred_tab$Hansen_mindist_logRS < QCR[1]), ] <- NA
result[which(pred_tab$Predominant_land_use == "Cropland" & pred_tab$Hansen_mindist_logRS > QCR[2]), ] <- NA


# Get the median, upper and lower quants for the plot
pred_tab$PredMedian <- ((apply(X = result,MARGIN = 1,
                               FUN = median,na.rm=TRUE))*100)-100
pred_tab$PredUpper <- ((apply(X = result,MARGIN = 1,
                              FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
pred_tab$PredLower <- ((apply(X = result,MARGIN = 1,
                              FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

pred_tab$dist_ori <- vals

pred_tab$realm <- "Tropical"

#pred_tab$Predominant_land_use <- sub(" vegetation", "", pred_tab$Predominant_land_use )
#pred_tab$Predominant_land_use <- factor(pred_tab$Predominant_land_use, levels = c("Primary", "Secondary", "Cropland"))


# SR plot = full range
ggplot(data = pred_tab) +
  geom_line(aes(x = dist_ori, y = PredMedian, col = Predominant_land_use)) +
  geom_ribbon(aes(x = dist_ori, ymin= PredLower, ymax = PredUpper, fill = Predominant_land_use), alpha = 0.3) +
  geom_rug(data = final.data.trans_trop_ABUN, aes(x = Hansen_mindist, col = Predominant_land_use), size = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  #facet_grid(~realm) +
  ylim(c(-80,80)) +
  xlim(c(0, 10)) +
  xlab("Distance to Forest (Km)") +
  ylab("Change in total abundance (%)") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"), labels = c("Primary", "Secondary", "Cropland"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E"), labels = c("Primary", "Secondary", "Cropland")) +
  theme_bw() +
  theme_custom
  


ggsave(filename = paste0(outdir, "/Dist_tropabunLU_refrow.pdf"), width = 3, height = 3, unit = "in")



#### Richness, Temperate ####


range(final.data.trans$Hansen_mindist)
# 0.0000 324.4697

from = 0
to = 325
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'Hansen_mindist_log'
fac <- NULL
n <- NULL
logval = TRUE



# organise the data
pred_tab <- sort_data(modout = srmod_temp,
                      moddata = final.data.trans_temp,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# reference for % difference = primary vegetation and distance closest to 0
refRow <- which((pred_tab$Predominant_land_use=="Cropland") & (pred_tab$Hansen_mindist_logRS==min(pred_tab$Hansen_mindist_logRS)))



# sort quantiles
Qdist <- quantile(x = srmod_temp$data$Hansen_mindist_logRS,
                  probs = exclQuantiles)


# predict the results
result <- PredictGLMERRandIter(model = srmod_temp$model,data = pred_tab, nIters = 10000)

# back transform the abundance values
result <- exp(result)

# convert to relative to reference
result <- sweep(x = result,MARGIN = 2,STATS = result[refRow,],FUN = '/')

# remove anything above and below the quantiles
result[which(pred_tab$Hansen_mindist_logRS < Qdist[1]), ] <- NA
result[which(pred_tab$Hansen_mindist_logRS > Qdist[2]), ] <- NA


# Get the median, upper and lower quants for the plot
pred_tab$PredMedian <- ((apply(X = result,MARGIN = 1,
                             FUN = median,na.rm=TRUE))*100)-100
pred_tab$PredUpper <- ((apply(X = result,MARGIN = 1,
                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
pred_tab$PredLower <- ((apply(X = result,MARGIN = 1,
                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

pred_tab$dist_ori <- vals

pred_tab$realm <- "Non-tropical"

# SR plot = full range
ggplot(data = pred_tab) +
  geom_line(aes(x = dist_ori, y = PredMedian), col = c("#EEAD0E")) +
  geom_ribbon(aes(x = dist_ori, ymin= PredLower, ymax = PredUpper), fill = c("#EEAD0E"), alpha = 0.3) +
  geom_rug(data = final.data.trans_trop_ABUN, aes(x = Hansen_mindist), size = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  #facet_grid(~realm) +
  ylim(c(-80,80)) +
  xlim(c(0, 20)) +
  xlab("Distance to Forest (Km)") +
  ylab("Change in species richness (%)") +
  #scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"))+
  #scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E")) +
  theme_bw() +
  theme_custom

ggsave(filename = paste0(outdir, "/Dist_richtemp_refrow.pdf"), width = 3, height = 3, uni = "in")





##%######################################################%##
#                                                          #
####                   homogeneity                      ####
#                                                          #
##%######################################################%##


#### abundance, tropical ####


range(final.data.trans$homogen)
# 0.1448 0.6686

from = 0.1
to = 0.7
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

# reference for % difference = primary vegetation and distance closest to 0
refRow <- which((pred_tab$Predominant_land_use=="Cropland") & (pred_tab$homogenRS==min(pred_tab$homogenRS)))


Qhomo <- quantile(x = abmod_trop$data$homogenRS,
  probs = exclQuantiles)


# predict the results
result <- PredictGLMERRandIter(model = abmod_trop$model,data = pred_tab, nIters = 10000)

# back transform the abundance values
result <- exp(result)-1

# convert to relative to reference
result <- sweep(x = result,MARGIN = 2,STATS = result[refRow,],FUN = '/')


# remove anything above and below the quantiles
result[which(pred_tab$homogenRS < Qhomo[1]), ] <- NA
result[which(pred_tab$homogenRS > Qhomo[2]), ] <- NA


# Get the median, upper and lower quants for the plot
pred_tab$PredMedian <- ((apply(X = result,MARGIN = 1,
                               FUN = median,na.rm=TRUE))*100)-100
pred_tab$PredUpper <- ((apply(X = result,MARGIN = 1,
                              FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
pred_tab$PredLower <- ((apply(X = result,MARGIN = 1,
                              FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

homogen_ori <- vals

ggplot(data = pred_tab) +
  geom_line(aes(x = homogen_ori, y = PredMedian), col = c("#EEAD0E")) +
  geom_ribbon(aes(x = homogen_ori, ymin= PredLower, ymax = PredUpper), fill = c("#EEAD0E"), alpha = 0.3) +
  geom_rug(data = final.data.trans_trop_ABUN, aes(x = homogen), size = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  ylim(c(-80,125)) +
  xlim(c(0, 1)) +
  xlab("Homogeneity") +
  ylab("Change in total abundance (%)") +
  theme_bw()  +
  theme_custom


ggsave(filename = paste0(outdir, "/Homogen_tropabun_refrow.pdf"), width = 3, height = 3, unit = "in")


#### Abundance, temperate ####



# organise the data
pred_tab2 <- sort_data(modout = abmod_temp,
                      moddata = final.data.trans_temp_ABUN,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# reference for % difference = primary vegetation and distance closest to 0
refRow <- which((pred_tab2$Predominant_land_use=="Cropland") & (pred_tab2$homogenRS==min(pred_tab2$homogenRS)))


# get quantiles of data
Qhomo2 <- quantile(x = abmod_temp$data$homogenRS,
                  probs = exclQuantiles)


# predict the results
result2 <- PredictGLMERRandIter(model = abmod_temp$model,data = pred_tab2, nIters = 10000)

# back transform the abundance values
result2 <- exp(result2)-1

# convert to relative to reference
result2 <- sweep(x = result2,MARGIN = 2,STATS = result2[refRow,],FUN = '/')

# remove anything above and below the quantiles
result2[which(pred_tab2$homogenRS < Qhomo[1]), ] <- NA
result2[which(pred_tab2$homogenRS > Qhomo[2]), ] <- NA


# Get the median, upper and lower quants for the plot
pred_tab2$PredMedian <- ((apply(X = result2,MARGIN = 1,
                             FUN = median,na.rm=TRUE))*100)-100
pred_tab2$PredUpper <- ((apply(X = result2,MARGIN = 1,
                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
pred_tab2$PredLower <- ((apply(X = result2,MARGIN = 1,
                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# add realm info
#pred_tab$realm <- "Tropical"
#pred_tab2$realm <- "Non-tropical"

ggplot(data = pred_tab2) +
  geom_line(aes(x = homogen_ori, y = PredMedian), col = c("#EEAD0E")) +
  geom_ribbon(aes(x = homogen_ori, ymin= PredLower, ymax = PredUpper), fill = c("#EEAD0E"), alpha = 0.3) +
  geom_rug(data = final.data.trans_temp_ABUN, aes(x = homogen), size = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(c(-80,125)) +
  xlim(c(0, 1)) +
  xlab("Homogeneity") +
  ylab("Change in total abundance (%)") +
  theme_bw()  +
  theme_custom

ggsave(filename = paste0(outdir, "/Homogen_tempabun_refrow.pdf"), width = 3, height = 3, unit = "in")


#### Richness, Tropical ####


fac <- 'Predominant_land_use'
n <- 3

# organise the data
pred_tab_3 <- sort_data(modout = srmod_trop,
                       moddata = final.data.trans_trop,
                       scalers = scalers,
                       from = from, 
                       to = to,
                       vals = vals, 
                       variable = variable,
                       fac = fac, 
                       n = n,
                       logval = logval)

# reference for % difference = primary vegetation and distance closest to 0
refRow <- which((pred_tab_3$Predominant_land_use=="Cropland") & (pred_tab_3$homogenRS==min(pred_tab_3$homogenRS)))



QPV <- quantile(x = srmod_trop$data$homogenRS[
  srmod_trop$data$Predominant_land_use=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = srmod_trop$data$homogenRS[
  srmod_trop$data$Predominant_land_use=="Secondary vegetation"],
  probs = exclQuantiles)
QCR <- quantile(x = srmod_trop$data$homogenRS[
  srmod_trop$data$Predominant_land_use=="Cropland"],
  probs = exclQuantiles)

# predict the results
result3 <- PredictGLMERRandIter(model = srmod_trop$model,data = pred_tab_3, nIters = 10000)

# back transform the abundance values
result3 <- exp(result3)

# convert to relative to reference
result3 <- sweep(x = result3,MARGIN = 2,STATS = result3[refRow,],FUN = '/')


# remove anything above and below the quantiles
result3[which(pred_tab_3$Predominant_land_use == "Primary vegetation" & pred_tab_3$homogenRS < QPV[1]), ] <- NA
result3[which(pred_tab_3$Predominant_land_use == "Primary vegetation" & pred_tab_3$homogenRS > QPV[2]), ] <- NA
result3[which(pred_tab_3$Predominant_land_use == "Secondary vegetation" & pred_tab_3$homogenRS < QSV[1]), ] <- NA
result3[which(pred_tab_3$Predominant_land_use == "Secondary vegetation" & pred_tab_3$homogenRS > QSV[2]), ] <- NA
result3[which(pred_tab_3$Predominant_land_use == "Cropland" & pred_tab_3$homogenRS < QCR[1]), ] <- NA
result3[which(pred_tab_3$Predominant_land_use == "Cropland" & pred_tab_3$homogenRS > QCR[2]), ] <- NA


# Get the median, upper and lower quants for the plot
pred_tab_3$PredMedian <- ((apply(X = result3,MARGIN = 1,
                              FUN = median,na.rm=TRUE))*100)-100
pred_tab_3$PredUpper <- ((apply(X = result3,MARGIN = 1,
                             FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
pred_tab_3$PredLower <- ((apply(X = result3,MARGIN = 1,
                             FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# add realm info
pred_tab_3$realm <- "Tropical"

pred_tab_3$homogen_ori <- vals

ggplot(data = pred_tab_3) +
  geom_line(aes(x = homogen_ori, y = PredMedian, col = Predominant_land_use)) +
  geom_ribbon(aes(x = homogen_ori, ymin= PredLower, ymax = PredUpper, fill = Predominant_land_use), alpha = 0.3) +
  geom_rug(data = final.data.trans_trop, aes(x = homogen, col = Predominant_land_use), size = 0.1) +
  #facet_grid(~ realm) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  ylim(c(-75,240)) +
  xlim(c(0, 1)) +
  xlab("Homogeneity") +
  ylab("Change in species richness (%)") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"), labels = c("Primary", "Secondary", "Cropland"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E"), labels = c("Primary", "Secondary", "Cropland")) +
  theme_bw() +
  theme_custom

ggsave(filename = paste0(outdir, "/HomogenLU_RichTrop_refrow.pdf"), width = 3, height = 3, unit = "in")



#### Richness, Non-tropical ####

# organise the data
pred_tab_4 <- sort_data(modout = srmod_temp,
                       moddata = final.data.trans_temp,
                       scalers = scalers,
                       from = from, 
                       to = to,
                       vals = vals, 
                       variable = variable,
                       fac = fac, 
                       n = n,
                       logval = logval)

# reference for % difference = primary vegetation and distance closest to 0
refRow <- which((pred_tab_4$Predominant_land_use=="Cropland") & (pred_tab_4$homogenRS==min(pred_tab_4$homogenRS)))



QPV <- quantile(x = srmod_temp$data$homogenRS[
  srmod_temp$data$Predominant_land_use=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = srmod_temp$data$homogenRS[
  srmod_temp$data$Predominant_land_use=="Secondary vegetation"],
  probs = exclQuantiles)
QCR <- quantile(x = srmod_temp$data$homogenRS[
  srmod_temp$data$Predominant_land_use=="Cropland"],
  probs = exclQuantiles)

# predict the results
result4 <- PredictGLMERRandIter(model = srmod_temp$model,data = pred_tab_4, nIters = 10000)

# back transform the abundance values
result4 <- exp(result4)

# convert to relative to reference
result4 <- sweep(x = result4,MARGIN = 2,STATS = result4[refRow,],FUN = '/')


# remove anything above and below the quantiles
result4[which(pred_tab_4$Predominant_land_use == "Primary vegetation" & pred_tab_4$homogenRS < QPV[1]), ] <- NA
result4[which(pred_tab_4$Predominant_land_use == "Primary vegetation" & pred_tab_4$homogenRS > QPV[2]), ] <- NA
result4[which(pred_tab_4$Predominant_land_use == "Secondary vegetation" & pred_tab_4$homogenRS < QSV[1]), ] <- NA
result4[which(pred_tab_4$Predominant_land_use == "Secondary vegetation" & pred_tab_4$homogenRS > QSV[2]), ] <- NA
result4[which(pred_tab_4$Predominant_land_use == "Cropland" & pred_tab_4$homogenRS < QCR[1]), ] <- NA
result4[which(pred_tab_4$Predominant_land_use == "Cropland" & pred_tab_4$homogenRS > QCR[2]), ] <- NA


# Get the median, upper and lower quants for the plot
pred_tab_4$PredMedian <- ((apply(X = result4,MARGIN = 1,
                              FUN = median,na.rm=TRUE))*100)-100
pred_tab_4$PredUpper <- ((apply(X = result4,MARGIN = 1,
                             FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
pred_tab_4$PredLower <- ((apply(X = result4,MARGIN = 1,
                             FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# add realm info
pred_tab_4$realm <- "Non-tropical"

pred_tab_4$homogen_ori <- vals

ggplot(data = pred_tab_4) +
  geom_line(aes(x = homogen_ori, y = PredMedian, col = Predominant_land_use)) +
  geom_ribbon(aes(x = homogen_ori, ymin= PredLower, ymax = PredUpper, fill = Predominant_land_use), alpha = 0.3) +
  geom_rug(data = final.data.trans_temp, aes(x = homogen, col = Predominant_land_use), size = 0.1) +
  #facet_grid(~ realm) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  ylim(c(-75,240)) +
  xlim(c(0, 1)) +
  xlab("Homogeneity") +
  ylab("Change in species richness (%)") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"), labels = c("Primary", "Secondary", "Cropland"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E"), labels = c("Primary", "Secondary", "Cropland")) +
  theme_bw() +
  theme_custom


ggsave(filename = paste0(outdir, "/Homogen_RichTemp.pdf"), width = 3, height = 3, unit = "in")




##%######################################################%##
#                                                          #
####                     Fertiliser                     ####
#                                                          #
##%######################################################%##

#### Tropical, richness ####

range(final.data.trans$fert.total)
# 0 2314226

from = 0
to = 2000
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'fert.total_log'
fac <- 'Predominant_land_use'
n <- 3
logval = TRUE

# organise the data
pred_tab <- sort_data(modout = srmod_trop,
                      moddata = final.data.trans_trop,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# reference for % difference = primary vegetation and distance closest to 0
refRow <- which((pred_tab$Predominant_land_use=="Cropland") & (pred_tab$fert.total_logRS==min(pred_tab$fert.total_logRS)))


# sort quantiles
QPV <- quantile(x = srmod_trop$data$fert.total_logRS[
  srmod_trop$data$Predominant_land_use=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = srmod_trop$data$fert.total_logRS[
  srmod_trop$data$Predominant_land_use=="Secondary vegetation"],
  probs = exclQuantiles)
QCR <- quantile(x = srmod_trop$data$fert.total_logRS[
  srmod_trop$data$Predominant_land_use=="Cropland"],
  probs = exclQuantiles)

# predict the results
result <- PredictGLMERRandIter(model = srmod_trop$model,data = pred_tab, nIters = 10000)

# back transform the abundance values
result <- exp(result)

# convert to relative to reference
result <- sweep(x = result,MARGIN = 2,STATS = result[refRow,],FUN = '/')

# remove anything above and below the quantiles
result[which(pred_tab$Predominant_land_use == "Primary vegetation" & pred_tab$fert.total_logRS < QPV[1]), ] <- NA
result[which(pred_tab$Predominant_land_use == "Primary vegetation" & pred_tab$fert.total_logRS > QPV[2]), ] <- NA
result[which(pred_tab$Predominant_land_use == "Secondary vegetation" & pred_tab$fert.total_logRS < QSV[1]), ] <- NA
result[which(pred_tab$Predominant_land_use == "Secondary vegetation" & pred_tab$fert.total_logRS > QSV[2]), ] <- NA
result[which(pred_tab$Predominant_land_use == "Cropland" & pred_tab$fert.total_logRS < QCR[1]), ] <- NA
result[which(pred_tab$Predominant_land_use == "Cropland" & pred_tab$fert.total_logRS > QCR[2]), ] <- NA

# Get the median, upper and lower quants for the plot
pred_tab$PredMedian <- ((apply(X = result,MARGIN = 1,
                               FUN = median,na.rm=TRUE))*100)-100
pred_tab$PredUpper <- ((apply(X = result,MARGIN = 1,
                              FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
pred_tab$PredLower <- ((apply(X = result,MARGIN = 1,
                              FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

pred_tab$fert_ori <- vals

pred_tab$realm <- "Tropical"

# SR plot = full range
ggplot(data = pred_tab) +
    geom_line(aes(x = fert_ori, y = PredMedian, col = Predominant_land_use)) +
    geom_ribbon(aes(x = fert_ori, ymin= PredLower, ymax = PredUpper, fill = Predominant_land_use), alpha = 0.3) +
    geom_rug(data = final.data.trans_trop_ABUN, aes(x = fert.total, col = Predominant_land_use), size = 0.1) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
    ylim(c(-100,150)) +
    xlim(c(0, 2000)) +
    xlab("Total fertiliser application (kgs)") +
    ylab("Change in species richness (%)") +
    scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"), labels = c("Primary", "Secondary", "Cropland"))+
    scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E"), labels = c("Primary", "Secondary", "Cropland")) +
    theme_bw() +
    theme_custom +
    theme(legend.position = "none")


ggsave(filename = paste0(outdir, "/FertLU_Richtrop_refrow.pdf"), width = 3, height = 3, unit = "in")


#### Richness, Temperate ####

from = 0
to = 2000
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'fert.total_log'
fac <- 'Use_intensity'
n <- 3
logval = TRUE


# organise the data
pred_tab_2 <- sort_data(modout = srmod_temp,
                        moddata = final.data.trans_temp,
                        scalers = scalers,
                        from = from, 
                        to = to,
                        vals = vals, 
                        variable = variable,
                        fac = fac, 
                        n = n,
                        logval = logval)


# reference for % difference = primary vegetation and distance closest to 0
refRow <- which((pred_tab_2$Use_intensity=="Intense use") & (pred_tab_2$fert.total_logRS==min(pred_tab_2$fert.total_logRS)))


# sort quantiles
QMU <- quantile(x = srmod_temp$data$fert.total_logRS[
  srmod_temp$data$Use_intensity=="Minimal use"],
  probs = exclQuantiles)
QLU <- quantile(x = srmod_temp$data$fert.total_logRS[
  srmod_temp$data$Use_intensity=="Light use"],
  probs = exclQuantiles)
QIU <- quantile(x = srmod_temp$data$fert.total_logRS[
  srmod_temp$data$Use_intensity=="Intense use"],
  probs = exclQuantiles)

# predict the results
result2 <- PredictGLMERRandIter(model = srmod_temp$model,data = pred_tab_2, nIters = 10000)

# back transform the abundance values
result2 <- exp(result2)

# convert to relative to reference
result2 <- sweep(x = result2,MARGIN = 2,STATS = result2[refRow,],FUN = '/')


# remove anything above and below the quantiles
result2[which(pred_tab_2$Use_intensity == "Minimal use" & pred_tab_2$fert.total_logRS < QMU[1]), ] <- NA
result2[which(pred_tab_2$Use_intensity == "Minimal use" & pred_tab_2$fert.total_logRS > QMU[2]), ] <- NA
result2[which(pred_tab_2$Use_intensity == "Light use" & pred_tab_2$fert.total_logRS < QLU[1]), ] <- NA
result2[which(pred_tab_2$Use_intensity == "Light use" & pred_tab_2$fert.total_logRS > QLU[2]), ] <- NA
result2[which(pred_tab_2$Use_intensity == "Intense use" & pred_tab_2$fert.total_logRS < QIU[1]), ] <- NA
result2[which(pred_tab_2$Use_intensity == "Intense use" & pred_tab_2$fert.total_logRS > QIU[2]), ] <- NA


# Get the median, upper and lower quants for the plot
pred_tab_2$PredMedian <- ((apply(X = result2,MARGIN = 1,
                                 FUN = median,na.rm=TRUE))*100)-100
pred_tab_2$PredUpper <- ((apply(X = result2,MARGIN = 1,
                                FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
pred_tab_2$PredLower <- ((apply(X = result2,MARGIN = 1,
                                FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

pred_tab_2$fert_ori <- vals
pred_tab_2$realm <- "Non-tropical"

# SR plot = full range
ggplot(data = pred_tab_2) +
  geom_line(aes(x = fert_ori, y = PredMedian, col = Use_intensity)) +
  geom_ribbon(aes(x = fert_ori, ymin= PredLower, ymax = PredUpper, fill = Use_intensity), alpha = 0.3) +
  geom_rug(data = final.data.trans_trop_ABUN, aes(x = fert.total, col = Use_intensity), size = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  ylim(c(-50,100)) +
  xlim(c(0, 2000)) +
  xlab("Total fertiliser application (kgs)") +
  ylab("Change in species richness (%)") +
  scale_colour_manual(values = c("#66CD00", "#FFB90F", "#EE0000"))+
  scale_fill_manual(values = c("#66CD00", "#FFB90F", "#EE0000")) +
  theme_bw() +
  theme_custom 

ggsave(filename = paste0(outdir, "/FertUI_richtemp_Refrow.pdf"), width = 3, height = 3, unit = "in")


#### Abundance, temperate ####

from = 0
to = 2000
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'fert.total_log'
fac <- 'Use_intensity'
n <- 3
logval = TRUE

# organise the data
pred_tab_3 <- sort_data(modout = abmod_temp,
                        moddata = final.data.trans_temp_ABUN,
                        scalers = scalers,
                        from = from, 
                        to = to,
                        vals = vals, 
                        variable = variable,
                        fac = fac, 
                        n = n,
                        logval = logval)

# reference for % difference = primary vegetation and distance closest to 0
refRow <- which((pred_tab_3$Use_intensity=="Intense use") & (pred_tab_3$fert.total_logRS==min(pred_tab_3$fert.total_logRS)))

# sort quantiles
QMU <- quantile(x = abmod_temp$data$fert.total_logRS[
  abmod_temp$data$Use_intensity=="Minimal use"],
  probs = exclQuantiles)
QLU <- quantile(x = abmod_temp$data$fert.total_logRS[
  abmod_temp$data$Use_intensity=="Light use"],
  probs = exclQuantiles)
QIU <- quantile(x = abmod_temp$data$fert.total_logRS[
  abmod_temp$data$Use_intensity=="Intense use"],
  probs = exclQuantiles)

# predict the results
result3 <- PredictGLMERRandIter(model = abmod_temp$model,data = pred_tab_3, nIters = 10000)

# back transform the abundance values
result3 <- exp(result3)-1

# convert to relative to reference
result3 <- sweep(x = result3,MARGIN = 2,STATS = result3[refRow,],FUN = '/')

# remove anything above and below the quantiles
result3[which(pred_tab_3$Use_intensity == "Minimal use" & pred_tab_3$fert.total_logRS < QMU[1]), ] <- NA
result3[which(pred_tab_3$Use_intensity == "Minimal use" & pred_tab_3$fert.total_logRS > QMU[2]), ] <- NA
result3[which(pred_tab_3$Use_intensity == "Light use" & pred_tab_3$fert.total_logRS < QLU[1]), ] <- NA
result3[which(pred_tab_3$Use_intensity == "Light use" & pred_tab_3$fert.total_logRS > QLU[2]), ] <- NA
result3[which(pred_tab_3$Use_intensity == "Intense use" & pred_tab_3$fert.total_logRS < QIU[1]), ] <- NA
result3[which(pred_tab_3$Use_intensity == "Intense use" & pred_tab_3$fert.total_logRS > QIU[2]), ] <- NA

# Get the median, upper and lower quants for the plot
pred_tab_3$PredMedian <- ((apply(X = result3,MARGIN = 1,
                                 FUN = median,na.rm=TRUE))*100)-100
pred_tab_3$PredUpper <- ((apply(X = result3,MARGIN = 1,
                                FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
pred_tab_3$PredLower <- ((apply(X = result3,MARGIN = 1,
                                FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

pred_tab_3$fert_ori <- vals

pred_tab_3$realm <- "Non-tropical"

# SR plot = full range
ggplot(data = pred_tab_3) +
geom_line(aes(x = fert_ori, y = PredMedian, col = Use_intensity)) +
geom_ribbon(aes(x = fert_ori, ymin= PredLower, ymax = PredUpper, fill = Use_intensity), alpha = 0.3) +
geom_rug(data = final.data.trans_trop_ABUN, aes(x = fert.total, col = Use_intensity), size = 0.1) +
geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
ylim(c(-100,100)) +
xlim(c(0, 2000)) +
xlab("Total fertiliser application (kgs)") +
ylab("Change in total abundance (%)") +
scale_colour_manual(values = c("#66CD00", "#FFB90F", "#EE0000"))+
scale_fill_manual(values = c("#66CD00", "#FFB90F", "#EE0000")) +
theme_bw() +
theme_custom + 
theme(legend.position = "none")

ggsave(filename = paste0(outdir, "/Fert_abuntemp_refrow.pdf"), width = 3, height = 3, unit = "in")




##%######################################################%##
#                                                          #
####             Percentage Natural Habitat             ####
#                                                          #
##%######################################################%##

#### Abundance, temperate ####

range(final.data.trans$percNH)
# 0.025 99.988

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

# reference for % difference = primary vegetation and distance closest to 0
refRow <- which((pred_tab$Predominant_land_use=="Cropland") & (pred_tab$percNHRS==min(pred_tab$percNHRS)))


# get quantiles to be excluded
Quant <- quantile(x = abmod_temp$data$percNHRS,
                  probs = exclQuantiles)


# predict the results
result <- PredictGLMERRandIter(model = abmod_temp$model,data = pred_tab, nIters = 10000)

# back transform the abundance values
result <- exp(result)-1

# convert to relative to reference
result <- sweep(x = result,MARGIN = 2,STATS = result[refRow,],FUN = '/')


# remove anything above and below the quantiles
result[which(pred_tab$percNHRS < Quant[1]), ] <- NA
result[which(pred_tab$percNHRS > Quant[2]), ] <- NA


# Get the median, upper and lower quants for the plot
pred_tab$PredMedian <- ((apply(X = result,MARGIN = 1,
                             FUN = median,na.rm=TRUE))*100)-100
pred_tab$PredUpper <- ((apply(X = result,MARGIN = 1,
                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
pred_tab$PredLower <- ((apply(X = result,MARGIN = 1,
                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

pred_tab$realm <- "Non-tropical"

pred_tab$percNH_ori <- vals

percNH <- as.data.frame(final.data.trans_temp_ABUN$percNH)
names(percNH) <- "V1"

# SR plot = full range
ggplot(data = pred_tab) +
  geom_line(aes(x = percNH_ori, y = PredMedian), col = c("#66CD00")) +
  geom_ribbon(aes(x = percNH_ori, ymin= PredLower, ymax = PredUpper), fill = c("#66CD00"), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1), size = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  ylim(c(-100,250)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Change in total abundance (%)") +
  theme_bw() +
  theme_custom

ggsave(filename = paste0(outdir, "/PercNH_Abuntemp_RefRow.pdf"), width = 3, height = 3, unit = "in")




#### Richness, Tropical ####


from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
fac <- 'Predominant_land_use'
n <- 3
logval = F

# organise the data
pred_tab <- sort_data(modout = srmod_trop,
                      moddata = final.data.trans_trop,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)


# reference for % difference = primary vegetation and distance closest to 0
refRow <- which((pred_tab$Predominant_land_use=="Cropland") & (pred_tab$percNHRS==min(pred_tab$percNHRS)))


# sort quantiles
QPV <- quantile(x = srmod_trop$data$percNHRS[
  srmod_trop$data$Predominant_land_use=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = srmod_trop$data$percNHRS[
  srmod_trop$data$Predominant_land_use=="Secondary vegetation"],
  probs = exclQuantiles)
QCR <- quantile(x = srmod_trop$data$percNHRS[
  srmod_trop$data$Predominant_land_use=="Cropland"],
  probs = exclQuantiles)

# predict the results
result <- PredictGLMERRandIter(model = srmod_trop$model,data = pred_tab, nIters = 10000)

# back transform
result <- exp(result)

# convert to relative to reference
result <- sweep(x = result,MARGIN = 2,STATS = result[refRow,],FUN = '/')


# remove anything above and below the quantiles
result[which(pred_tab$Predominant_land_use == "Primary vegetation" & pred_tab$percNHRS < QPV[1]), ] <- NA
result[which(pred_tab$Predominant_land_use == "Primary vegetation" & pred_tab$percNHRS > QPV[2]), ] <- NA
result[which(pred_tab$Predominant_land_use == "Secondary vegetation" & pred_tab$percNHRS < QSV[1]), ] <- NA
result[which(pred_tab$Predominant_land_use == "Secondary vegetation" & pred_tab$percNHRS > QSV[2]), ] <- NA
result[which(pred_tab$Predominant_land_use == "Cropland" & pred_tab$percNHRS < QCR[1]), ] <- NA
result[which(pred_tab$Predominant_land_use == "Cropland" & pred_tab$percNHRS > QCR[2]), ] <- NA


# Get the median, upper and lower quants for the plot
pred_tab$PredMedian <- ((apply(X = result,MARGIN = 1,
                             FUN = median,na.rm=TRUE))*100)-100
pred_tab$PredUpper <- ((apply(X = result,MARGIN = 1,
                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
pred_tab$PredLower <- ((apply(X = result,MARGIN = 1,
                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

pred_tab$percNH_ori <- vals
pred_tab$realm <- "Tropical"

percNH <- as.data.frame(final.data.trans_trop$percNH)
percNH$Predominant_land_use <- final.data.trans_trop$Predominant_land_use
names(percNH)[1] <- "V1"

# SR plot = full range
ggplot(data = pred_tab) +
  geom_line(aes(x = percNH_ori, y = PredMedian, col = Predominant_land_use)) +
  geom_ribbon(aes(x = percNH_ori, ymin= PredLower, ymax = PredUpper, fill = Predominant_land_use), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1, col = Predominant_land_use), size = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  ylim(c(-50,250)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Change in species richness (%)") +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"), labels = c("Primary", "Secondary", "Cropland"))+
  scale_fill_manual(values = c("#006400", "#8B0000", "#EEAD0E"), labels = c("Primary", "Secondary", "Cropland")) +
  theme_bw() +
  theme_custom + 
  theme(legend.position = c(0.2, 0.8))

ggsave(filename = paste0(outdir, "/PercNHLU_tropRich_refRow.pdf"), width = 3, height = 3, unit = "in")



#### Richness, Tropical ####

from = 0
to = 100
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'percNH'
fac <- 'Use_intensity'
n <- 3
logval = F

# organise the data
pred_tab <- sort_data(modout = srmod_trop,
                      moddata = final.data.trans_trop,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# reference for % difference = primary vegetation and distance closest to 0
refRow <- which((pred_tab$Use_intensity=="Intense use") & (pred_tab$percNHRS==min(pred_tab$percNHRS)))


# sort quantiles
QMU <- quantile(x = srmod_trop$data$percNHRS[
  srmod_trop$data$Use_intensity=="Minimal use"],
  probs = exclQuantiles)
QLU <- quantile(x = srmod_trop$data$percNHRS[
  srmod_trop$data$Use_intensity=="Light use"],
  probs = exclQuantiles)
QIU <- quantile(x = srmod_trop$data$percNHRS[
  srmod_trop$data$Use_intensity=="Intense use"],
  probs = exclQuantiles)

# predict the result
result <- PredictGLMERRandIter(model = srmod_trop$model,data = pred_tab, nIters = 10000)

# transform the results
result <- exp(result)

# convert to relative to reference
result <- sweep(x = result,MARGIN = 2,STATS = result[refRow,],FUN = '/')

# remove anything above and below the quantiles
result[which(pred_tab$Use_intensity == "Minimal use" & pred_tab$percNHRS < QMU[1]), ] <- NA
result[which(pred_tab$Use_intensity == "Minimal use" & pred_tab$percNHRS > QMU[2]), ] <- NA
result[which(pred_tab$Use_intensity == "Light use" & pred_tab$percNHRS < QLU[1]), ] <- NA
result[which(pred_tab$Use_intensity == "Light use" & pred_tab$percNHRS > QLU[2]), ] <- NA
result[which(pred_tab$Use_intensity == "Intense use" & pred_tab$percNHRS < QIU[1]), ] <- NA
result[which(pred_tab$Use_intensity == "Intense use" & pred_tab$percNHRS > QIU[2]), ] <- NA


# Get the median, upper and lower quants for the plot
pred_tab$PredMedian <- ((apply(X = result,MARGIN = 1,
                             FUN = median,na.rm=TRUE))*100)-100
pred_tab$PredUpper <- ((apply(X = result,MARGIN = 1,
                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
pred_tab$PredLower <- ((apply(X = result,MARGIN = 1,
                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

pred_tab$percNH_ori <- vals

pred_tab$realm <- "Tropical"

## organise data for plotting ##

# add the new vals


percNH <- as.data.frame(final.data.trans_trop$percNH)
percNH$Use_intensity <- final.data.trans_trop$Use_intensity
names(percNH)[1] <- "V1"

# SR plot = full range
ggplot(data = pred_tab) +
  geom_line(aes(x = percNH_ori, y = PredMedian, col = Use_intensity)) +
  geom_ribbon(aes(x = percNH_ori, ymin= PredLower, ymax = PredUpper, fill = Use_intensity), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1, col = Use_intensity), size = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  ylim(c(-50, 250)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Change in species richness (%)") +
  scale_colour_manual(values = c("#66CD00", "#FFB90F", "#EE0000"))+
  scale_fill_manual(values = c("#66CD00", "#FFB90F", "#EE0000")) +
  theme_bw() +
  theme_custom +
  theme(legend.position = c(0.2, 0.8))

ggsave(filename = paste0(outdir, "/Supp_Richness_percNHUI_trop_refRow.pdf"), width = 3, height = 3, unit = "in")

#### Richness,  temperate ####

# organise the data
pred_tab2 <- sort_data(modout = srmod_temp,
                      moddata = final.data.trans_temp,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# reference for % difference = primary vegetation and distance closest to 0
refRow <- which((pred_tab2$Use_intensity=="Intense use") & (pred_tab2$percNHRS==min(pred_tab2$percNHRS)))


# sort quantiles
QMU <- quantile(x = srmod_temp$data$percNHRS[
  srmod_temp$data$Use_intensity=="Minimal use"],
  probs = exclQuantiles)
QLU <- quantile(x = srmod_temp$data$percNHRS[
  srmod_temp$data$Use_intensity=="Light use"],
  probs = exclQuantiles)
QIU <- quantile(x = srmod_temp$data$percNHRS[
  srmod_temp$data$Use_intensity=="Intense use"],
  probs = exclQuantiles)

# predict the result
result2 <- PredictGLMERRandIter(model = srmod_temp$model, data = pred_tab2, nIters = 10000)

# transform the results
result2 <- exp(result2)

# convert to relative to reference
result2 <- sweep(x = result2,MARGIN = 2,STATS = result2[refRow,],FUN = '/')

# remove anything above and below the quantiles
result2[which(pred_tab2$Use_intensity == "Minimal use" & pred_tab2$percNHRS < QMU[1]), ] <- NA
result2[which(pred_tab2$Use_intensity == "Minimal use" & pred_tab2$percNHRS > QMU[2]), ] <- NA
result2[which(pred_tab2$Use_intensity == "Light use" & pred_tab2$percNHRS < QLU[1]), ] <- NA
result2[which(pred_tab2$Use_intensity == "Light use" & pred_tab2$percNHRS > QLU[2]), ] <- NA
result2[which(pred_tab2$Use_intensity == "Intense use" & pred_tab2$percNHRS < QIU[1]), ] <- NA
result2[which(pred_tab2$Use_intensity == "Intense use" & pred_tab2$percNHRS > QIU[2]), ] <- NA


# Get the median, upper and lower quants for the plot
pred_tab2$PredMedian <- ((apply(X = result2,MARGIN = 1,
                             FUN = median,na.rm=TRUE))*100)-100
pred_tab2$PredUpper <- ((apply(X = result2,MARGIN = 1,
                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
pred_tab2$PredLower <- ((apply(X = result2,MARGIN = 1,
                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

pred_tab2$percNH_ori <- vals

pred_tab2$realm <- "Non-tropical"

## organise data for plotting ##

# add the new vals


percNH <- as.data.frame(final.data.trans_temp$percNH)
percNH$Use_intensity <- final.data.trans_temp$Use_intensity
names(percNH)[1] <- "V1"

# SR plot = full range
ggplot(data = pred_tab2) +
  geom_line(aes(x = percNH_ori, y = PredMedian, col = Use_intensity)) +
  geom_ribbon(aes(x = percNH_ori, ymin= PredLower, ymax = PredUpper, fill = Use_intensity), alpha = 0.3) +
  geom_rug(data = percNH, aes(x = V1, col = Use_intensity), size = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  ylim(c(-50, 250)) +
  xlim(c(0, 100)) +
  xlab("Percentage of Natural Habitat") +
  ylab("Change in species richness (%)") +
  scale_colour_manual(values = c("#66CD00", "#FFB90F", "#EE0000"))+
  scale_fill_manual(values = c("#66CD00", "#FFB90F", "#EE0000")) +
  theme_bw() +
  theme_custom +
  theme(legend.position = "none")


ggsave(filename = paste0(outdir, "/PercNHLU_RichTemp_refRow.pdf"), width = 3, height = 3, unit = "in")




##%######################################################%##
#                                                          #
####               Number of landcovers                 ####
#                                                          #
##%######################################################%##


#### abundance, temperate ####
range(final.data.trans$landcovers.5k)

from = 1
to = 10
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'landcovers.5k'
fac <- NULL
n <- NULL
logval = F


# organise the data
pred_tab <- sort_data(modout = abmod_temp,
                      moddata = final.data.trans_temp,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# reference for % difference = primary vegetation and distance closest to 0
refRow <- which((pred_tab$Use_intensity=="Intense use") & (pred_tab$landcovers.5kRS==min(pred_tab$landcovers.5kRS)))


# sort quantiles
Qlc <- quantile(x = abmod_temp$data$landcovers.5kRS,
                  probs = exclQuantiles)


# predict the result
result <- PredictGLMERRandIter(model = abmod_temp$model, data = pred_tab, nIters = 10000)

# transform the results
result <- exp(result)-1

# convert to relative to reference
result <- sweep(x = result,MARGIN = 2,STATS = result[refRow,],FUN = '/')


# remove anything above and below the quantiles
result[which(pred_tab$landcovers.5kRS < Qlc[1]), ] <- NA
result[which(pred_tab$landcovers.5kRS > Qlc[2]), ] <- NA


# Get the median, upper and lower quants for the plot
pred_tab$PredMedian <- ((apply(X = result,MARGIN = 1,
                             FUN = median,na.rm=TRUE))*100)-100
pred_tab$PredUpper <- ((apply(X = result,MARGIN = 1,
                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
pred_tab$PredLower <- ((apply(X = result,MARGIN = 1,
                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

pred_tab$landc_ori <- vals

pred_tab$realm <- "Non-tropical"



# SR plot = full range
ggplot(data = pred_tab) +
  geom_line(aes(x = landc_ori, y = PredMedian), col = c("#EEAD0E")) +
  geom_ribbon(aes(x = landc_ori, ymin= PredLower, ymax = PredUpper), fill = c("#EEAD0E"), alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  ylim(c(-50,100)) +
  xlim(c(0, 11)) +
  xlab("Number of Landcovers") +
  ylab("Change in total abundance (%)") +
  theme_bw() +
  theme_custom

ggsave(filename = paste0(outdir, "/Landcovers_abuntemp_refRow.pdf"), width = 3, height = 3, unit = "in")




#### Richness, Tropical ####

from = 1
to = 10
vals <- seq(from = from, to = to, length.out = 1000)
variable <- 'landcovers.5k'
fac <- "Use_intensity"
n <- 3
logval = F


# organise the data
pred_tab <- sort_data(modout = srmod_trop,
                      moddata = final.data.trans_trop,
                      scalers = scalers,
                      from = from, 
                      to = to,
                      vals = vals, 
                      variable = variable,
                      fac = fac, 
                      n = n,
                      logval = logval)

# reference for % difference = primary vegetation and distance closest to 0
refRow <- which((pred_tab$Use_intensity=="Intense use") & (pred_tab$landcovers.5kRS==min(pred_tab$landcovers.5kRS)))


# sort quantiles
QMU <- quantile(x = srmod_trop$data$landcovers.5kRS[
  srmod_trop$data$Use_intensity=="Minimal use"],
  probs = exclQuantiles)
QLU <- quantile(x = srmod_trop$data$landcovers.5kRS[
  srmod_trop$data$Use_intensity=="Light use"],
  probs = exclQuantiles)
QIU <- quantile(x = srmod_trop$data$landcovers.5kRS[
  srmod_trop$data$Use_intensity=="Intense use"],
  probs = exclQuantiles)

# predict the result
result <- PredictGLMERRandIter(model = srmod_trop$model, data = pred_tab, nIters = 10000)

# transform the results
result <- exp(result)

# convert to relative to reference
result <- sweep(x = result,MARGIN = 2,STATS = result[refRow,],FUN = '/')


# remove anything above and below the quantiles
result[which(pred_tab$Use_intensity == "Minimal use" & pred_tab$landcovers.5kRS < QMU[1]), ] <- NA
result[which(pred_tab$Use_intensity == "Minimal use" & pred_tab$landcovers.5kRS > QMU[2]), ] <- NA
result[which(pred_tab$Use_intensity == "Light use" & pred_tab$landcovers.5kRS < QLU[1]), ] <- NA
result[which(pred_tab$Use_intensity == "Light use" & pred_tab$landcovers.5kRS > QLU[2]), ] <- NA
result[which(pred_tab$Use_intensity == "Intense use" & pred_tab$landcovers.5kRS < QIU[1]), ] <- NA
result[which(pred_tab$Use_intensity == "Intense use" & pred_tab$landcovers.5kRS > QIU[2]), ] <- NA


# Get the median, upper and lower quants for the plot
pred_tab$PredMedian <- ((apply(X = result,MARGIN = 1,
                             FUN = median,na.rm=TRUE))*100)-100
pred_tab$PredUpper <- ((apply(X = result,MARGIN = 1,
                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
pred_tab$PredLower <- ((apply(X = result,MARGIN = 1,
                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

pred_tab$landc_ori <- vals

pred_tab$realm <- "Tropical"

## organise data for plotting ##

# add the new vals


landcov <- as.data.frame(final.data.trans_trop$landcovers.5k)
landcov$Use_intensity <- final.data.trans_trop$Use_intensity
names(landcov)[1] <- "V1"

# SR plot = full range
ggplot(data = pred_tab) +
  geom_line(aes(x = landc_ori, y = PredMedian, col = Use_intensity)) +
  geom_ribbon(aes(x = landc_ori, ymin= PredLower, ymax = PredUpper, fill = Use_intensity), alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  #facet_grid(~realm) +
  ylim(c(-50, 100)) +
  #xlim(c(0, 10)) +
  xlab("Number of Landcovers") +
  ylab("Change in species richness (%)") +
  scale_colour_manual(values = c("#66CD00", "#FFB90F", "#EE0000"))+
  scale_fill_manual(values = c("#66CD00", "#FFB90F", "#EE0000")) +
  theme_bw() +
  theme_custom + 
  theme(legend.position = c(0.2, 0.85))

ggsave(filename = paste0(outdir, "/LandcoversUI_richtrop_refRow.pdf"), width = 3, height = 3, unit = "in")


#### Abundance, Tropical ####

# organise the data
pred_tab2 <- sort_data(modout = abmod_trop,
                       moddata = final.data.trans_trop_ABUN,
                       scalers = scalers,
                       from = from, 
                       to = to,
                       vals = vals, 
                       variable = variable,
                       fac = fac, 
                       n = n,
                       logval = logval)

# reference for % difference = primary vegetation and distance closest to 0
refRow <- which((pred_tab2$Use_intensity=="Intense use") & (pred_tab2$landcovers.5kRS==min(pred_tab2$landcovers.5kRS)))



# sort quantiles
QMU <- quantile(x = abmod_trop$data$landcovers.5kRS[
  abmod_trop$data$Use_intensity=="Minimal use"],
  probs = exclQuantiles)
QLU <- quantile(x = abmod_trop$data$landcovers.5kRS[
  abmod_trop$data$Use_intensity=="Light use"],
  probs = exclQuantiles)
QIU <- quantile(x = abmod_trop$data$landcovers.5kRS[
  abmod_trop$data$Use_intensity=="Intense use"],
  probs = exclQuantiles)

# predict the result
result2 <- PredictGLMERRandIter(model = abmod_trop$model, data = pred_tab2, nIters = 10000)

# transform the results
result2 <- exp(result2)-1

# convert to relative to reference
result2 <- sweep(x = result2,MARGIN = 2,STATS = result2[refRow,],FUN = '/')


# remove anything above and below the quantiles
result2[which(pred_tab2$Use_intensity == "Minimal use" & pred_tab2$landcovers.5kRS < QMU[1]), ] <- NA
result2[which(pred_tab2$Use_intensity == "Minimal use" & pred_tab2$landcovers.5kRS > QMU[2]), ] <- NA
result2[which(pred_tab2$Use_intensity == "Light use" & pred_tab2$landcovers.5kRS < QLU[1]), ] <- NA
result2[which(pred_tab2$Use_intensity == "Light use" & pred_tab2$landcovers.5kRS > QLU[2]), ] <- NA
result2[which(pred_tab2$Use_intensity == "Intense use" & pred_tab2$landcovers.5kRS < QIU[1]), ] <- NA
result2[which(pred_tab2$Use_intensity == "Intense use" & pred_tab2$landcovers.5kRS > QIU[2]), ] <- NA


# Get the median, upper and lower quants for the plot
pred_tab2$PredMedian <- ((apply(X = result2,MARGIN = 1,
                              FUN = median,na.rm=TRUE))*100)-100
pred_tab2$PredUpper <- ((apply(X = result2,MARGIN = 1,
                             FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
pred_tab2$PredLower <- ((apply(X = result2,MARGIN = 1,
                             FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

pred_tab2$landc_ori <- vals

pred_tab2$realm <- "Tropical"

## organise data for plotting ##

# add the new vals


landc <- as.data.frame(final.data.trans_trop_ABUN$landcovers.5k)
landc$Use_intensity <- final.data.trans_trop_ABUN$Use_intensity
names(landc)[1] <- "V1"

# SR plot = full range
ggplot(data = pred_tab2) +
  geom_line(aes(x = landc_ori, y = PredMedian, col = Use_intensity)) +
  geom_ribbon(aes(x = landc_ori, ymin= PredLower, ymax = PredUpper, fill = Use_intensity), alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  #facet_grid(~realm) +
  ylim(c(-75, 100)) +
  xlim(c(0, 10)) +
  xlab("Number of Landcovers") +
  ylab("Change in total abundance (%)") +
  scale_colour_manual(values = c("#66CD00", "#FFB90F", "#EE0000"))+
  scale_fill_manual(values = c("#66CD00", "#FFB90F", "#EE0000")) +
  theme_bw() +
  theme_custom + 
  theme(legend.position = "none")



ggsave(filename = paste0(outdir, "/LandcoversUI_Abun_trop_refRow.pdf"), width = 3, height = 3, unit = "in")




##%######################################################%##
#                                                          #
####                  Abundance LU UI                   ####
#                                                          #
##%######################################################%##


#### Abundance, Tropical ####


# basic table of median values and reference factors
pred_tab <- data.frame(landcovers.5kRS = median(final.data.trans_trop_ABUN$landcovers.5kRS),
                       homogenRS = median(final.data.trans_trop_ABUN$homogenRS),
                       fert.total_logRS = median(final.data.trans_trop_ABUN$fert.total_logRS),
                       percNHRS = median(final.data.trans_trop_ABUN$percNHRS),
                       Hansen_mindist_logRS =  median(final.data.trans_trop_ABUN$Hansen_mindist_logRS),
                       Forest_biome = "Tropical & Subtropical Moist Broadleaf Forests",
                       Use_intensity = "Minimal use",
                       Predominant_land_use = "Primary vegetation",
                       #Tropical = "Temperate",
                       Species_richness = 0,
                       logAbun = 0)


# organise factor levels
# check levels of factor variables


pred_tab$Predominant_land_use <- factor(pred_tab$Predominant_land_use, levels = levels(abmod_trop$data$Predominant_land_use))
pred_tab$Use_intensity <- factor(pred_tab$Use_intensity, levels = levels(abmod_trop$data$Use_intensity)) 
pred_tab$Forest_biome <- factor(pred_tab$Forest_biome, levels = levels(abmod_trop$data$Forest_biome)[c(3, 2, 1)]) 


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

pred_tab$median <- ((apply(X = resulta, MARGIN = 1, FUN = median))*100)-100
pred_tab$upper <- ((apply(X = resulta, MARGIN = 1, FUN = quantile,probs = 0.975))*100)-100
pred_tab$lower <- ((apply(X = resulta, MARGIN = 1, FUN = quantile,probs = 0.025))*100)-100

pred_tab$realm <- "Tropical"

# temperate


# basic table of median values and reference factors
pred_tab2 <- data.frame(landcovers.5kRS = median(final.data.trans_temp_ABUN$landcovers.5kRS),
                        homogenRS = median(final.data.trans_temp_ABUN$homogenRS),
                        fert.total_logRS = median(final.data.trans_temp_ABUN$fert.total_logRS),
                        percNHRS = median(final.data.trans_temp_ABUN$percNHRS),
                        Hansen_mindist_logRS =  median(final.data.trans_temp_ABUN$Hansen_mindist_logRS),
                        Forest_biome = "Temperate Broadleaf & Mixed Forests",
                        Use_intensity = "Minimal use",
                        Predominant_land_use = "Primary vegetation",
                        #Tropical = "Temperate",
                        Species_richness = 0,
                        logAbun = 0)


# organise factor levels
# check levels of factor variables

pred_tab2$Predominant_land_use <- factor(pred_tab2$Predominant_land_use, levels = levels(abmod_temp$data$Predominant_land_use))
pred_tab2$Use_intensity <- factor(pred_tab2$Use_intensity, levels = levels(abmod_temp$data$Use_intensity)) 
pred_tab2$Forest_biome <- factor(pred_tab2$Forest_biome, levels = levels(abmod_temp$data$Forest_biome)) 

# add and change factor levels of land use and intensity

pred_tab2 <- do.call("rbind", replicate(9, pred_tab2, simplify = FALSE))


pred_tab2[4:6, 'Predominant_land_use'] <- "Secondary vegetation"
pred_tab2[7:9, 'Predominant_land_use'] <- "Cropland"


pred_tab2[c(2,5,8), 'Use_intensity'] <- "Light use"
pred_tab2[c(3,6,9), 'Use_intensity'] <- "Intense use"


### Abundance, Temperate ####

# predict the result
resulta2 <- PredictGLMERRandIter(model = abmod_temp$model, data = pred_tab2)

# transform the results
resulta2 <- exp(resulta2)-1

resulta2 <- sweep(x = resulta2, MARGIN = 2, STATS = resulta2[1,], FUN = '/')

pred_tab2$median <- ((apply(X = resulta2, MARGIN = 1, FUN = median))*100)-100
pred_tab2$upper <- ((apply(X = resulta2, MARGIN = 1, FUN = quantile,probs = 0.975))*100)-100
pred_tab2$lower <- ((apply(X = resulta2, MARGIN = 1, FUN = quantile,probs = 0.025))*100)-100

pred_tab2$realm <- "Non-tropical"

plot_data <- rbind(pred_tab, pred_tab2)

plot_data[plot_data$upper == 0, c("upper", "lower")] <- NA

plot_data$realm <- factor(plot_data$realm, levels = c("Tropical", "Non-tropical"))

plot_data$Predominant_land_use <- sub("Primary vegetation", "Primary\nvegetation", plot_data$Predominant_land_use)
plot_data$Predominant_land_use <- sub("Secondary vegetation", "Secondary\nvegetation", plot_data$Predominant_land_use)

plot_data$Predominant_land_use <- factor(plot_data$Predominant_land_use, levels = c("Primary\nvegetation", "Secondary\nvegetation", "Cropland"))

p1 <- ggplot(data = plot_data)+
  geom_point(aes(x = Predominant_land_use, y = median, col = Predominant_land_use, shape = Use_intensity),
             position = position_dodge(width = 0.9), size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper, y = median, x = Predominant_land_use, col = Predominant_land_use),
                position = position_dodge2(padding = 0.5)) +
  facet_grid(~ realm) +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"), guide = F)+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2)+
  xlab("") +
  ylab("Total Abundance (%)") +
  theme_bw() +
  theme_custom + 
  theme(legend.position = c(0.15, 0.8), strip.background = element_rect(fill = NA, size = 0.2))
  

#ggsave(filename = paste0(outdir, "/Supp_Abun_LUUI.pdf"), width = 6, height = 3, unit = "in")



##%######################################################%##
#                                                          #
####                  Richness LU UI                    ####
#                                                          #
##%######################################################%##



# basic table of median values and reference factors
pred_tab3 <- data.frame(landcovers.5kRS = median(final.data.trans_trop$landcovers.5kRS),
                       homogenRS = median(final.data.trans_trop$homogenRS),
                       fert.total_logRS = median(final.data.trans_trop$fert.total_logRS),
                       percNHRS = median(final.data.trans_trop$percNHRS),
                       Hansen_mindist_logRS =  median(final.data.trans_trop$Hansen_mindist_logRS),
                       Forest_biome = "Tropical & Subtropical Moist Broadleaf Forests",
                       Use_intensity = "Minimal use",
                       Predominant_land_use = "Primary vegetation",
                       #Tropical = "Temperate",
                       Species_richness = 0,
                       logAbun = 0)


# organise factor levels
# check levels of factor variables
pred_tab3$Predominant_land_use <- factor(pred_tab3$Predominant_land_use, levels = levels(srmod_trop$data$Predominant_land_use))
pred_tab3$Use_intensity <- factor(pred_tab3$Use_intensity, levels = levels(srmod_trop$data$Use_intensity)) 
pred_tab3$Forest_biome <- factor(pred_tab3$Forest_biome, levels = levels(srmod_trop$data$Forest_biome)[c(3, 2, 1)]) 


# add and change factor levels of land use and intensity

pred_tab3 <- do.call("rbind", replicate(9, pred_tab3, simplify = FALSE))


pred_tab3[4:6, 'Predominant_land_use'] <- "Secondary vegetation"
pred_tab3[7:9, 'Predominant_land_use'] <- "Cropland"

pred_tab3[c(2,5,8), 'Use_intensity'] <- "Light use"
pred_tab3[c(3,6,9), 'Use_intensity'] <- "Intense use"


### Tropical predictions ###

# predict the result
resulta <- PredictGLMERRandIter(model = srmod_trop$model, data = pred_tab3)

# transform the results
resulta <- exp(resulta)

resulta <- sweep(x = resulta, MARGIN = 2, STATS = resulta[1,], FUN = '/')

pred_tab3$median <- ((apply(X = resulta, MARGIN = 1, FUN = median))*100)-100
pred_tab3$upper <- ((apply(X = resulta, MARGIN = 1, FUN = quantile,probs = 0.975))*100)-100
pred_tab3$lower <- ((apply(X = resulta, MARGIN = 1, FUN = quantile,probs = 0.025))*100)-100

pred_tab3$realm <- "Tropical"

# temperate


# basic table of median values and reference factors
pred_tab4 <- data.frame(landcovers.5kRS = median(final.data.trans_temp$landcovers.5kRS),
                        homogenRS = median(final.data.trans_temp$homogenRS),
                        fert.total_logRS = median(final.data.trans_temp$fert.total_logRS),
                        percNHRS = median(final.data.trans_temp$percNHRS),
                        Hansen_mindist_logRS =  median(final.data.trans_temp$Hansen_mindist_logRS),
                        Forest_biome = "Temperate Broadleaf & Mixed Forests",
                        Use_intensity = "Minimal use",
                        Predominant_land_use = "Primary vegetation",
                        #Tropical = "Temperate",
                        Species_richness = 0,
                        logAbun = 0)


# organise factor levels
# check levels of factor variables

pred_tab4$Predominant_land_use <- factor(pred_tab4$Predominant_land_use, levels = levels(srmod_temp$data$Predominant_land_use))
pred_tab4$Use_intensity <- factor(pred_tab4$Use_intensity, levels = levels(srmod_temp$data$Use_intensity)) 
pred_tab4$Forest_biome <- factor(pred_tab4$Forest_biome, levels = levels(srmod_temp$data$Forest_biome)) 


# add and change factor levels of land use and intensity

pred_tab4 <- do.call("rbind", replicate(9, pred_tab4, simplify = FALSE))


pred_tab4[4:6, 'Predominant_land_use'] <- "Secondary vegetation"
pred_tab4[7:9, 'Predominant_land_use'] <- "Cropland"


pred_tab4[c(2,5,8), 'Use_intensity'] <- "Light use"
pred_tab4[c(3,6,9), 'Use_intensity'] <- "Intense use"


#### Richness, Temperate ####

# predict the result
resulta2 <- PredictGLMERRandIter(model = srmod_temp$model, data = pred_tab4)

# transform the results
resulta2 <- exp(resulta2)

resulta2 <- sweep(x = resulta2, MARGIN = 2, STATS = resulta2[1,], FUN = '/')

pred_tab4$median <- ((apply(X = resulta2, MARGIN = 1, FUN = median))*100)-100
pred_tab4$upper <- ((apply(X = resulta2, MARGIN = 1, FUN = quantile,probs = 0.975))*100)-100
pred_tab4$lower <- ((apply(X = resulta2, MARGIN = 1, FUN = quantile,probs = 0.025))*100)-100

pred_tab4$realm <- "Non-tropical"

plot_data2 <- rbind(pred_tab3, pred_tab4)

plot_data2[plot_data2$upper == 0, c("upper", "lower")] <- NA

plot_data2$realm <- factor(plot_data2$realm, levels = c("Tropical", "Non-tropical"))

plot_data2$Predominant_land_use <- sub("Primary vegetation", "Primary\nvegetation", plot_data2$Predominant_land_use)
plot_data2$Predominant_land_use <- sub("Secondary vegetation", "Secondary\nvegetation", plot_data2$Predominant_land_use)

plot_data2$Predominant_land_use <- factor(plot_data2$Predominant_land_use, levels = c("Primary\nvegetation", "Secondary\nvegetation", "Cropland"))

p2 <- ggplot(data = plot_data2)+
  geom_point(aes(x = Predominant_land_use, y = median, col = Predominant_land_use, shape = Use_intensity),
             position = position_dodge(width = 0.9), size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper, y = median, x = Predominant_land_use, col = Predominant_land_use),
                position = position_dodge2(padding = 0.5)) +
  facet_grid(~ realm) +
  scale_colour_manual(values = c("#006400", "#8B0000", "#EEAD0E"), guide = F)+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2)+
  xlab("") +
  ylab("Species Richness (%)") +
  theme_bw() +
  theme_custom +
  theme(legend.position = "none", strip.background = element_rect(fill = NA, size = 0.2))

#ggsave(filename = paste0(outdir, "/Supp_Rich_LUUI.pdf"), width = 6, height = 3, uni = "in")


cowplot::plot_grid(p1, p2, nrow = 2, labels = c("(a)", "(b)"), label_size = 10)

ggsave(filename = paste0(outdir, "/Figure2_LUUI.pdf"), width = 6, height = 6, uni = "in")
