#this script contains all statistical analyses for the manuscript "Plant Metabarcoding as a Tool for Dietary Composition Analysis: Successes and Limitations"
#all analyses run in R version 3.6.1 Action of the Toes
#last updated: 9/25/21

setwd("<>") #set working directory

#packages
if(!require(psych)){install.packages("psych")}
if(!require(lmtest)){install.packages("lmtest")}
if(!require(rcompanion)){install.packages("rcompanion")}
if(!require(multcompView)){install.packages("multcompView")}
library(ggplot2)
library(betareg)
library(emmeans)
library(glmmTMB)
library(sjstats)
library(ggpubr)
library(DHARMa)

#data
beta.data <- read.csv("rats.only.csv") #contains sequences assigned to taxonomy and converted to relative read abundance and percent error data for all fecal samples 

#############Are the trials significantly different when on the 20% diet? 
mod.trial2 <- glmmTMB(P.CreosoteP ~ Trial, family=beta_family(link="logit"), data=beta.data.trial)
summary(mod.trial2)
joint_tests(mod.trial2) 

#############Does Juniper Correlate with the Amount of chow in the diet?
mod.chow <- glmmTMB(P.JuniperP ~ P.ChowP + (1|AID), family=beta_family(link="logit"), data=beta.data.cont4)
summary(mod.chow) 
mod.chow.res <- simulateResiduals(mod.chow)
plot(mod.chow.res)
 
#############Does Juniper Correlate With the Amount of Creosote in the Diet?
mod1 <- glmmTMB(P.JuniperP ~ P.CreosoteP + (1|AID), family=beta_family(link="logit"), data=beta.data.cont4)
summary(mod1)
mod1.res <- simulateResiduals(mod1)
plot(mod1.res)
ggplot(beta.data, aes( x=P.CreosoteP, y = P.JuniperP)) + geom_point()
#quasi r squared as most packages cannot handle beta distributed and zero inflated models.
1 - exp((2/nrow(beta.data4)) * (logLik(update(mod1, ~1))[1] - logLik(mod1)[1]))


###############Are there significant differences in the amount of creosote in each diet?
mod.creosote <- glmmTMB(P.CreosoteP ~ ActualP + (1|AID), family=beta_family(link="logit"), data=beta.data4, ziformula = ~1)
summary(mod.creosote) 
#residuals
residuals(mod.creosote, type = c("response", "pearson"))
mod.creosote.res <- simulateResiduals(mod.creosote)
plot(mod.creosote.res)
joint_tests(mod.creosote)
marginal.creosote <- emmeans(mod.creosote, ~ ActualP)
pairs(marginal.creosote, adjust ="tukey")


##########Are there significant differences in the amount of juniper per each diet? 
mod.juniper <- glmmTMB(P.JuniperP ~ ActualP + (1|AID), family=beta_family(link="logit"), data=beta.data4)
summary(mod.juniper)
mod.juniper.res <- simulateResiduals(mod.juniper.res)
plot(mod.juniper.res)
joint_tests(mod.juniper)
marginal.juniper <- emmeans(mod.juniper, ~ ActualP)
pairs(marginal.juniper, adjust ="tukey")


###########Is there a difference in percent error by diet?
#Does Juniper have difference in percent error by diet? 
error.jun <- beta.data4$Juniper.Error/100
mod.error.jun <- glmmTMB(error.jun ~ ActualP + (1|AID), family=Gamma(link ="log"), data=beta.data4)
summary(mod.error.jun)
#residuals
mod.juniper.res <- simulateResiduals(mod.juniper.res)
plot(mod.juniper.res)
joint_tests(mod.error.jun)
marginal.error.jun <- emmeans(mod.error.jun, ~ ActualP)
pairs(marginal.error.jun, adjust ="tukey") 


############Does creosote have difference in percent error by diet?
error.creo <- beta.data4$Creosote.Error/100
mod.error.creo <- glmmTMB(error.creo ~ ActualP + (1|AID), family=Gamma(link="log"), data=beta.data4)
summary(mod.error.creo)
#residuals
mod.creo.error.res <- simulateResiduals(mod.error.creo)
plot(mod.creo.error.res)
joint_tests(mod.error.creo)
marginal.error.creo <- emmeans(mod.error.creo, ~ ActualP)
pairs(marginal.error.creo, adjust ="tukey")


#########Pipeline Analysis, is pipline accuracy significantly different for our different diet components?
accuracy <- read.csv("new.accuracy.data.csv")
accuracy$AID = factor(accuracy$AID,
                      levels=unique(accuracy$AID))

accuracy$Trial = factor(accuracy$Trial,
                        levels=unique(accuracy$Trial))

#Creosote
mod.creo.accuracy <- glmmTMB(Creosote ~ Pipeline + (1|AID), data=accuracy, family=Gamma(link="log"))
summary(mod.creo.accuracy)
#residuals
mod.creo.pipeline.res <- simulateResiduals(mod.creo.accuracy)
plot(mod.creo.pipeline.res)
joint_tests(mod.creo.accuracy) 

#Juniper
mod.juniper.accuracy <- glmmTMB(Juniper ~ Pipeline + (1|AID), data=accuracy, family=Gamma(link="log"))
summary(mod.juniper.accuracy)
#residuals
mod.jun.pipeline.res <- simulateResiduals(mod.juniper.accuracy)
plot(mod.jun.pipeline.res)
joint_tests(mod.juniper.accuracy) 

#Chow
mod.chow.accuracy <- glmmTMB(Chow ~ Pipeline + (1|AID), data=accuracy, family=Gamma(link="log"))
summary(mod.chow.accuracy)
#residuals
mod.chow.pipeline.res <- simulateResiduals(mod.chow.accuracy)
plot(mod.chow.pipeline.res)
joint_tests(mod.chow.accuracy)

###############Do pooled samples have higher accuracy than non-pooled samples? 
pooled <- read.csv("pooled.vs.mean.csv")
pooled.accuracy <- glmmTMB(Percent.Error ~ Type + (1|Sample), data=pooled, family=Gamma(link="log"))
summary(pooled.accuracy)
#residuals
mod.pooled.res <- simulateResiduals(pooled.accuracy)
plot(mod.pooled.res)
joint_tests(pooled.accuracy)

########Comparing Accuracy of Diet Samples vs Fecal Samples
diet <- read.csv("new.plant.data2.csv") #table containing sequence data converted to relative read abundance for both diet and fecal samples

#change AID to factor
diet4$AID = factor(diet4$AID,
                   levels=unique(diet4$AID))
#model
fecal.accuracy <- glmmTMB(Accuracy ~ Rat + (1|Sample), data=diet4, family=Gamma(link="log"))
summary(fecal.accuracy)
#residuals
mod.fecal.res <- simulateResiduals(fecal.accuracy)
plot(mod.fecal.res)
joint_tests(fecal.accuracy) 
