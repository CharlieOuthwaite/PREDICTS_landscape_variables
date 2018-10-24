##%######################################################%##
#                                                          #
####              Organising PREDICTS data              ####
#                                                          #
##%######################################################%##


# Start off just looking at the abundance data

# subset the dataset to include just those where metric is abundance
# need to check which column is best for this, 

rm(list = ls())

# where is the data?
datadir <- "C:/Users/chaout/Dropbox/POSTDOC - biodiv & food sec/0. PROJECTS/1. ForestCover_Yield/1. Data/PREDICTS_2016"

# load data
PRED.data <- read.csv(paste0(datadir, "/", "resource.csv")) # takes a long time to read in!


# subset to cropland sites only
crop.data <- PRED.data[PRED.data$Predominant_land_use == "Cropland", ]

# subset to abundance based metrics
crop.data.abun <- crop.data[crop.data$Diversity_metric_type == "Abundance", ]

table(crop.data.abun$Diversity_metric)

# do I need to remove any forms of diversity metric - sign density, percent cover???
crop.data.abun <- crop.data.abun[!crop.data.abun$Diversity_metric == "percent cover", ]


# subset to useful columns
crop.data.abun <-crop.data.abun[, c(6:10, 19:21, 27, 31:34, 38:39, 65:68)]

# proportion of records that are non-zero
1-nrow(crop.data.abun[crop.data.abun$Measurement == 0,])/nrow(crop.data.abun)


#### summarise to site level information ####
site_abun <- aggregate(crop.data.abun$Effort_corrected_measurement, by = crop.data.abun["SSBS"], FUN = sum)

# add in the additional info

site_abun$Predominant_land_use <- crop.data.abun[match(site_abun$SSBS, crop.data.abun$SSBS), 'Predominant_land_use']

site_abun$Use_intensity <- crop.data.abun[match(site_abun$SSBS, crop.data.abun$SSBS), 'Use_intensity']

site_abun$SS <- crop.data.abun[match(site_abun$SSBS, crop.data.abun$SSBS), 'SS']

site_abun$SSB <- crop.data.abun[match(site_abun$SSBS, crop.data.abun$SSBS), 'SSB']



# remove any NAs
model.data <- na.omit(site_abun)



#### try out a basic model of effect of use intensity on abundance in cropland ####

# remove the "cannot decide" category from land use intensity
model.data <- model.data[!model.data$Use_intensity == "Cannot decide", ] # now 1936 seperate blocks...

# reorder use intensity factors so that least impactful is first
model.data$Use_intensity <- as.factor(as.character(model.data$Use_intensity))

model.data$Use_intensity <- factor(as.character(model.data$Use_intensity), 
                                   levels = c("Minimal use", "Light use", "Intense use"))

# log transform the abundance data
colnames(model.data)[2] <- "Abundance"
model.data$logAbun <- log(model.data$Abundance )

## following example models in the tutorial.  Not including human pop density or 
## distance to road as these data are not in the PREDICTS download from NHM.

library(lme4)

# doesn't like the infs, not sure what is normally done here so just converting to 0 for now

model.data[model.data$logAbun == -Inf, 'logAbun'] <- 0

random1 <- lmer(logAbun~Use_intensity+(1|SS), data=model.data)

summary(random1)

random2 <- lmer(logAbun~Use_intensity+(1|SS)+(1|SSB), data=model.data)

summary(random2)

# checking which is better
AIC(random1,random2)


### trying out some plotting ###

library(StatisticalModels)

PlotGLMERFactor(model = random2,data = model.data,
                responseVar = "logAbun",logLink = "e",catEffects = "Use_intensity",
                xtext.srt = 45)


### testing for over-dispersion ###

GLMEROverdispersion(model = random2)

# dont think there is overdispersion but testing another model anyway.

random3 <- lmer(logAbun~Use_intensity+(1|SS)+(1|SSB)+(1|SSBS), data=model.data) # doesn't work

