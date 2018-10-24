##%######################################################%##
#                                                          #
####                 PREDICTS Tutorial                  ####
#                                                          #
##%######################################################%##

rm(list = ls())

# load libraries
install.packages("lme4")
install.packages("arm")
install.packages("emdbook")
install.packages("bbmle")
install.packages("R2WinBUGS")
# Alternatively, if you are using a Mac, install this package instead
install.packages("rjags")
# If the next line doesn't work, come to see me for a version
# of the StatisticalModels and MResModelling packages
install.packages("devtools")
library(devtools)
install.packages("survival")
install.packages("Formula")
install.packages("ggplot2")
install.packages("Hmisc")

# load data
install_github("timnewbold/StatisticalModels")
install_github("timnewbold/MResEcologicalModelling",subdir="MResModelling")


#### Exercise 1: Functional responses - Linear models and maximum likelihood estimation ####

# In this exercise, we will use data on predator functional responses 
# of East African reed frogs, from Vonesh & Bolker (2005).

library(emdbook)
data(ReedfrogFuncresp)

# First, let's plot the functional response to inspect the data:

plot(ReedfrogFuncresp$Initial,ReedfrogFuncresp$Killed)

# It looks as though there might be a linear relationship between these variables. 
# So let's try fitting a simple linear model:

m1 <- lm(Killed~Initial,data=ReedfrogFuncresp)

m1
# Coefficients:
# (Intercept)      Initial  
#       2.727        0.276

summary(m1)
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  2.72651    1.96292   1.389    0.187    
# Initial      0.27603    0.03948   6.991 6.33e-06 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 5.04 on 14 degrees of freedom
# Multiple R-squared:  0.7773,    Adjusted R-squared:  0.7614 
# F-statistic: 48.88 on 1 and 14 DF,  p-value: 6.334e-06


# It is not easy to transform the variables to make a linear relationship from this model,
# and even if we did do this the parameters would be difficult to interpret. 
# So instead we will fit the model using maximum likelihood estimation. 
# First though, as an experiment, let's refit the linear model using maximum likelihood.

# We need to load the bblme package
library(bbmle)

# First we define the likelihood function
linearNLL <- function(killed,init,m,c,sd){

# The following line calculates the expected y values
# for a given slope and intercept
y.pred <- m * init + c

# The next line calculates the likelihood for this model
# assuming that the residuals are a normal distribution
# around 0 (as in a linear regression)
# The 'dnorm' function calculates the probability density
# for a given observed value (x) compared to a normal
# distribution with a given mean and standard deviation (sd)
suppressWarnings(-sum(dnorm(x = killed,mean = y.pred,sd = sd,log = TRUE)))

}

# Now we run the maximum likelihood estimation
# This will estimate the values of three parameters:
# the slope (m) and intercept(c) of the regression,
# plus the residual error term (sd). For These# parameters,
# we specify reasonable starting values.
# We also pass two variables (our response and 
# explanatory variables): intial prey density 
# and the number killed
m2 <- mle2(minuslogl = linearNLL,start = list(m=0.2,c=2,sd=1),
data = list(init=ReedfrogFuncresp$Initial,
killed=ReedfrogFuncresp$Killed))

m2
# Coefficients:
#         m         c        sd 
# 0.2760248 2.7265615 4.7141473 


# The parameters are of course very similar to those estimated by the linear regression. 
# Before moving on to the more complex functional response model, 
# let's plot the data and the fitted linear relationship:

plot(ReedfrogFuncresp$Initial,ReedfrogFuncresp$Killed)

# Create a data frame with values of initial density to predict
preds <- data.frame(Initial=1:100)
# Predict number of prey killed
preds$Killed <- m2@coef['m']*preds$Initial+m2@coef['c']
# Plot (type="l" plots a line, and lwd=2 makes the line thick)
points(preds$Initial,preds$Killed,type="l",lwd=2,col="red")


# Now we will fit the more complex 'Type II' functional response model to the data.


# As before, we first need to define the likelihood function
binomLL <- function(killed,init,p){
  
  # The next two lines just specify that of the parameters we pass (p), 
  # the first will be attack rate (a), and the second the handling time (h)
  a = p[1]
  h = p[2]
  
  # This line defines the probability of a prey individual being killed
  # Note this is the Type II functional response equation from above, divided by
  # the intial number of prey
  pkilled <- a/(1+a*h*init)
  
  # The likelihood in this case assumes that the observed number 
  # of kills is a binomial draw where the number of trials is
  # the initial number of prey, and the probability of being killed
  # is taken from above
  -sum(dbinom(x = killed,size = init,prob = pkilled,log = TRUE))
  
}
# This line just defines the names for the parameters
parnames(binomLL) <- c("a","h")

# Now run the likelihood estimation
m3 <- mle2(minuslogl = binomLL,start = c(a=0.5,h=0.0125),
           data = list(init=ReedfrogFuncresp$Initial,
                       killed=ReedfrogFuncresp$Killed))

m3
# Coefficients:
#          a          h 
# 0.52630319 0.01664362 
# 
# Log-likelihood: -46.72

# If we compare the AIC values of these models, we can see that the Type II functional 
# response is a slightly better fit to the data (i.e. has a lower AIC value - see the 
# lecture for a reminder about AIC values):
  
  AIC(m2,m3)
#         AIC df
# 1 101.02429  3
# 2  97.44279  2
  
# Finally, we will plot the fitted Type II functional response. Note that the model 
# as defined above predicts the probability that a single prey individual is eaten,
# so we need to multiply by initial prey density to get the estimate of the total number of prey eaten.
  
  plot(ReedfrogFuncresp$Initial,ReedfrogFuncresp$Killed)
  
  preds <- data.frame(Initial=1:100)
  preds$Killed <- preds$Initial *
    m3@coef['a']/(1 + m3@coef['a'] * 
                    m3@coef['h'] * preds$Initial)
  points(preds$Initial,preds$Killed,type="l",lwd=2,col="red")


#### Exercise 2: Land use impacts on Hymenoptera in Hawaii - Generalized linear models ####    

# First, load my R package associated with these sessions:
library(MResModelling)
 
# Now, load in the Hawaii data describing the abundance of all species:
data(HawaiiHymenoptera)
  
# We will work with the data for just one species, Pimpla punicipes:
    
hh.sp <- droplevels(hh[(hh$Taxon_name_entered=="Pimpla punicipes"),])
  
hist(hh.sp$Measurement)

# For now we will analyze species presence or absence. To do so we will create a new 
# column with values of 1 where abundance is greater than zero, and values of 0 
# otherwise (i.e. where the species is absent):
hh.sp$PresAbs <- ifelse(hh.sp$Measurement>0,1,0)

# Plot the data
# The first line of code just tabulates the presences and absences by land use
presabs.counts <- table(hh.sp$Occur,hh.sp$LandUse)

# The next line divides the counts through by the column sums to calculate proportions
presabs.props <- sweep(presabs.counts,2,colSums(presabs.counts),'/')

# This line reorders the table, so that presences come before absences, 
# and the land uses are in a logical order (primary vegetation followed 
# by secondary vegetation followed by pasture)
presabs.props <- presabs.props[c(2,1),c(1,3,2)]

# Finally, we create a barplot of the proportions
# (green for presence, orange for absence)
barplot(presabs.props,col=c("#1b9e77","#d95f02"),las=1,ylab="Proportion of sites")

# Because the response variable that we are interested in (species presence or absence)
# is a binary variable, it most appropriate to use a GLM that assumes a binomial 
# distribution of errors.

# Run a model of species presence or absence as a function of land use
m1 <- glm(PresAbs~LandUse,data=hh.sp,family=binomial)

print(summary(m1))

# Also, run a null model (with just an intercept) with which to compare our land-use model
m0 <- glm(PresAbs~1,data=hh.sp,family=binomial)    

# We can plot an error bar to show the modelled result:
  
# First, create a set of data giving just the land uses we want to show
nd <- data.frame(LandUse=c("Primary","Secondary","Pasture"))
# Now predict probability of presence for these land uses, including standard error
preds <- predict(object=m1,newdata=nd,se.fit=TRUE)
# Calculate the mean predicted value for each land use, 
# back transforming the values using the inverse of the 
# link function (for a binomial GLM, the default link 
# function is logit)
y <- 1/(1+exp(-(preds$fit)))
# Now calculate the upper and lower confidence limits (multiplying by 1.96 gives 95% confidence limits)
yplus <- 1/(1+exp(-(preds$fit+1.96*preds$se.fit)))
yminus <- 1/(1+exp(-(preds$fit-1.96*preds$se.fit)))
# Finally, plot the error bar with the mean and confidence limits of the model predictions for each land use
errbar(x=nd$LandUse,y=y,yplus=yplus,yminus=yminus)

# Now do an analysis of variance to compare the model with land use to the null model:
  
anova(m0,m1,test="Chi")

# This shows that the model including land use is significantly better than the 
# null model. So land use does have a significant effect on the presence or absence 
# of this species at this particular study location in Hawaii.

# Finally, though, let's find out how much of the variation in species presence or 
# absence is explained by land use as a measure of fit of the model. Remember from the 
# lecture, that explained variation in a GLM is (null deviance - residual deviance)/null deviance.

with(m1,(null.deviance-deviance)/null.deviance)
# Note that residual deviance is labelled simply 'deviance' in the model object.

# This shows that our model explains only around 5% of the variation 
# in species presence or absence. So land use has a significant effect on the species,
# but clearly there are other factors (not considered in this model) that determine to a 
# large extent whether or not the species is present and is detected at each location in the dataset.

# Now we will work with data on the counts of all hymenopteran species at the 754 sampled sites in Hawaii:
  
data(HawaiiHymenopteraSites)

# We will work with the species richness data. 
# First, produce a boxplot of species richness as a function of land use:
  
boxplot(Species_richness~LandUse,data=hhs)

# But first, let's inspect the distribution of the species richness data:

hist(hhs$Species_richness)

# Clearly the species richness data are not normally distributed. 
# Because species richness values are counts, it would be reasonable to start 
# with a GLM with a Poisson error distribution.

# Model species richness as a function of land use
m1 <- glm(Species_richness~LandUse,data=hhs,family=poisson)

print(summary(m1))

# Also run a null, intercept-only model
m0 <- glm(Species_richness~1,data=hhs,family=poisson)

# As before, we will plot an error bar to show the modelled result:
  
# First, create a set of data with the land uses we want to show
nd <- data.frame(LandUse=c("Primary","Secondary","Pasture"))
# Now predict probability of presence for these land uses, including standard error
preds <- predict(object=m1,newdata=nd,se.fit=TRUE)
# Calculate the mean predicted value for each land use, back transforming the values using the inverse of the link function (this time, because we ran a Poisson GLM, the model used a log link function)
y <- exp(preds$fit)
# Now calculate the upper and lower confidence limits (multiplying by 1.96 gives 95% confidence limits)
yplus <- exp(preds$fit+1.96*preds$se.fit)
yminus <- exp(preds$fit-1.96*preds$se.fit)
# Finally, plot the error bar with the mean and confidence limits of the model predictions for each land use
errbar(x=nd$LandUse,y=y,yplus=yplus,yminus=yminus)

# So in this case, the model shows that species richness is, on average, 
# slightly higher in secondary vegetation than in primary vegetation, but 
# much lower in pasture than in either of the natural land-use tpyes.

# Again, do an analysis of variance to compare the land-use model with the null model:
  
anova(m0,m1,test="Chi")

# As with the model of the presence or absence of the individual species, 
# the model including land use is highly significantly better than the null model. 
# But again, let's see how much of the variation in species richness land use explains:

with(m1,(null.deviance-deviance)/null.deviance)

# So land use explains nearly 17% of the variation in species richness. 
# That's a lot for a single explanatory variable.


#### Exercise 3: Land use impacts on biodiversity - Mixed-effects models ####

# use a mixed-effects model to look at complete PREDICTS data and effect of land use,
# as in Newbold paper.

# First, load my StatisticalModels package, and read in the site-level data from Newbold et al. (2016):
  
library(StatisticalModels) 
data(PREDICTSSiteData)

# reorder land use factor levels
PREDICTSSites$LandUse <- factor(PREDICTSSites$LandUse,
                                levels=c("Primary Vegetation","Secondary Vegetation",
                                         "Plantation forest","Cropland","Pasture","Urban"))

# First, we must decide which random effects to include in the model. As a minimum, 
# we will include study identity in these models, because we expect sampled biodiversity 
# to vary considerably among studies owing to differences in sampling methods and sampling 
# effort (in this dataset, study identity designated as 'SS'):


# We need to create a modelling dataset without NAs,
# because the mixed-effects modelling functions don't handle NAs
# First, select the variables we will include in the models
# (I have included some extra terms that we will consider later)
model.data <- PREDICTSSites[,c('LogAbund','LandUse','logHPD.rs','logDistRd.rs','SS','SSB')]

# Second, remove NA values
model.data <- na.omit(model.data)

# Fixed effects are specified as in a normal linear model,
# Random effects are specified as e.g. (1|random_term)+(1|random_term2)
random1 <- lmer(LogAbund~LandUse+poly(logHPD.rs,2)+poly(logDistRd.rs,2)+(1|SS),
                data=model.data)



# here accounting for blocks within studies
random2 <- lmer(LogAbund~LandUse+poly(logHPD.rs,2)+poly(logDistRd.rs,2)+(1|SS)+(1|SSB),
                data=model.data)

AIC(random1,random2) # model 2 with spatial structuring of sites within studies is better

## once random effects structure is identified, then select the fixed effect structure.
# usually start off with all then remove, backward stepwise model selection.  Drop each
# variable and see if there is a reduction in the explained variance.  

# GLMERSelect function in Tim's package does this selection process for you
abundModelSelect <- GLMERSelect(modelData = PREDICTSSites,responseVar = "LogAbund",
                                fitFamily = "gaussian",fixedFactors = "LandUse",
                                fixedTerms = list(logHPD.rs=2,logDistRd.rs=2),
                                randomStruct = "(1|SS)+(1|SSB)",verbose = TRUE)

abundModelSelect$stats

# We can also look at the table of coefficient values from the model. From the table of 
# fixed effects, we can see that total abundance is lower in all disturbed land uses 
# compared with primary vegetation.

summary(abundModelSelect$model)

# We can plot an error bar showing the land-use effects from the final model 
# using the PlotGLMERFactor in my StatisticalModels package:
  
PlotGLMERFactor(model = abundModelSelect$model,data = abundModelSelect$data,
               responseVar = "Abundance",logLink = "e",catEffects = "LandUse",
                xtext.srt = 45)

# And we can plot the effect of human population density using PlotGLMERContinuous:
  
PlotGLMERContinuous(model = abundModelSelect$model,data = abundModelSelect$data,
                      effects = "logHPD.rs",otherFactors = list(LandUse="Primary Vegetation"),
                      xlab = "(Log) Human population density (rescaled)",ylab = "Total abundance",
                      logLink = "e")


 
# let's compare the study-only and study-plus-spatial-block random-effects structures 
# that we considered for the models of total abundance. You will probably receive warnings 
# about model convergence, but these aren't too serious (the convergence value is close 
# to the accepted threshold, and the functions in the lme4 package are very cautious about convergence).

model.data <- PREDICTSSites[,c('Species_richness','LandUse','logHPD.rs','logDistRd.rs','SS','SSB','SSBS')]
model.data <- na.omit(model.data)

# Here we will fit a <i>generalized<i> linear model, assuming a Poisson
# distribution of errors
randomR1 <- glmer(Species_richness~LandUse+poly(logHPD.rs,2)+poly(logDistRd.rs,2)+(1|SS),family="poisson",data=model.data)
randomR2 <- glmer(Species_richness~LandUse+poly(logHPD.rs,2)+poly(logDistRd.rs,2)+(1|SS)+(1|SSB),family="poisson",data=model.data)

AIC(randomR1,randomR2)
#                df      AIC
# randomR1 11 102591.5
# randomR2 12 100695.0

# As with the models of total abundance, the random-effects structure that includes 
# the effect of the spatial structure of sites is strongly favoured.

# One way to test for over-dispersion is to compare the residual deviance to the residual 
# degrees of freedom of a model. If the deviance is much larger than the degrees of freedom, 
# this is an indication of over-dispersion (there other possible reasons though). 
# There is a function in my StatisticalModels package that does this:

GLMEROverdispersion(model = randomR2)
# $residDev
# [1] 31852.27
# 
# $residDF
# [1] 15621
# 
# $ratio
# [1] 2.039068
# 
# $P.ChiSq
# [1] 0

# The large ratio of residual deviance to residual degrees of freedom indicates the 
# presence of over-dispersion (confirmed by the significant chi-square test). 
# Therefore, we will try a random-effects structure with a nested effect of site, 
# within spatial block, within study. (warning, this will take a little time to run!):

randomR3 <- glmer(Species_richness~LandUse+poly(logHPD.rs,2)+poly(logDistRd.rs,2)+(1|SS)+(1|SSB)+(1|SSBS),family="poisson",data=model.data)

AIC(randomR2,randomR3)

# The comparison of AIC values suggests that the random-effects structure including site 
# is strongly favoured, and re-running the over-dispersion test shows that including an 
# observation-level random effect has removed the over-dispersion (in fact, there is 
# now under-dispersion):
  
  GLMEROverdispersion(model = randomR3)
# $residDev
# [1] 7448.711
# 
# $residDF
# [1] 15620
# 
# $ratio
# [1] 0.4768701
# 
# $P.ChiSq
# [1] 1
  
# Now we have selected a random-effects structure for our species richness model, 
# we can perform backward stepwise model selection (again, this might take a while to run!):
    
richModelSelect <- GLMERSelect(modelData = PREDICTSSites,responseVar = "Species_richness",
                                   fitFamily = "poisson",fixedFactors = "LandUse",
                                   fixedTerms = list(logHPD.rs=2,logDistRd.rs=2),
                                   randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",verbose = TRUE)

richModelSelect$stats

summary(richModelSelect$model)

# plotting effects of land use, human population density and distance to nearest road.

PlotGLMERFactor(model = richModelSelect$model,data = richModelSelect$data,
                responseVar = "Species richness",logLink = "e",catEffects = "LandUse",
                xtext.srt = 45)

PlotGLMERContinuous(model = richModelSelect$model,data = richModelSelect$data,
                    effects = c("logHPD.rs"),otherContEffects = c("logDistRd.rs"),
                    otherFactors = list(LandUse="Primary Vegetation"),
                    xlab = "(Log) Human population density (rescaled)",
                    ylab = "Species richness",logLink = "e")

PlotGLMERContinuous(model = richModelSelect$model,data = richModelSelect$data,
                    effects = c("logDistRd.rs"),otherContEffects = c("logHPD.rs"),
                    otherFactors = list(LandUse="Primary Vegetation"),
                    xlab = "(Log) Distance to road (rescaled)",
                    ylab = "Species richness",logLink = "e")
