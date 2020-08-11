#  ------------------------------------------------------------------------

#  Vitamin D ms
#  Fitness analyses 
#  AMS - 24/7/20

#  ------------------------------------------------------------------------

# remove all data

rm(list=ls())

# read in data - 1452 obs

VITD <- read.csv("VitD_MS_Data_240720.csv")

# libraries 

library(plyr)
library(dplyr)
library(AICcmodavg) 
library(ggplot2) 
library(lme4)
library(jtools)
library(ordinal)
library(glmmTMB)
library(bbmle)
library(sjPlot)
library(MASS)
library(cowplot)
library(data.table)
library(DHARMa)

# set factors etc

VITD$ID<-as.factor(VITD$ID)
VITD$YearF<-as.factor(VITD$Year)
VITD$Sex<-as.factor(VITD$Sex)
VITD$CoatBinF<-as.factor(VITD$CoatBin)
VITD$BirthYearF<-as.factor(VITD$BirthYear)
VITD$AgeSq<-(VITD$Age*VITD$Age)
VITD$AgeF<-as.factor(VITD$Age)
VITD$MumID<-as.factor(VITD$MumID)
VITD$AgeGroupF<-as.factor(VITD$AgeGroup)
VITD$EweFecundityO<-ordered(VITD$EweFecundity)
VITD$SurvivalF <- as.factor(VITD$Survival)


# # -----------------------------------------------------------------------


# Correlation between 25(OH)D2 and 25(OH)D3


# # -----------------------------------------------------------------------

ggplot(VITD,aes(x=X25OHD2,y=X25OHD3)) + geom_point(shape=21, alpha=0.5, colour="grey4")+
  labs(x=expression("25(OH)D"["2"]~"(nmol/l)"), y=expression("25(OH)D"["3"]~"(nmol/l)")) + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=12)) + 
  theme(axis.title.y=element_text(colour="grey6", size=12)) +
  theme(axis.line=element_line(colour="grey6")) +
  theme(text = element_text(size=12, colour="grey16")) 

# correlation and paired t-test

cor.test(x=VITD$X25OHD2,y=VITD$X25OHD3)
t.test(VITD$X25OHD2, VITD$X25OHD3, paired=TRUE)


# # -----------------------------------------------------------------------


# Phenotypic associations with vitamin D levels


# # -----------------------------------------------------------------------

## Total vitamin D

ggplot(VITD,aes(x=CoatBinF,y=Total.25D)) + geom_boxplot(outlier.colour=NA,notch=FALSE)  + xlab("Coat colour") + ylab("25(OH)D") +
  theme_classic() + theme(legend.position="none") + geom_jitter() +
  theme(axis.title.x=element_text(face="bold", colour="grey16", size=12)) + 
  theme(axis.title.y=element_text(face="bold", colour="grey16", size=12)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=14, colour="grey26")) + 
  scale_x_discrete(limit = c("1", "2"), labels = c("Dark","Light"))

# model

totalvitdmodel <- lmer(Total.25D~Sex+Age+AgeSq+CoatBinF+YearF + (1|ID) +(1|BirthYearF),
                       data=VITD)

plot(totalvitdmodel)
qqnorm(resid(totalvitdmodel))

summary(totalvitdmodel)

# sig of sex - x=3.123, df=1, p=0.07719

totalvitdmodel2<-update(totalvitdmodel, ~ . -Sex)
anova(totalvitdmodel, totalvitdmodel2)

# sig of age sq - x=219.16, df=1, p<0.001

totalvitdmodel3<-update(totalvitdmodel, ~ . -AgeSq)
anova(totalvitdmodel, totalvitdmodel3)

# sig of coat bin - x=16.526      df=1    p=4.8e-05

totalvitdmodel4<-update(totalvitdmodel, ~ . -CoatBinF)
anova(totalvitdmodel, totalvitdmodel4)

# sig of year - x=111.36      df=5  p< 2.2e-16

totalvitdmodel5<-update(totalvitdmodel, ~ . -YearF)
anova(totalvitdmodel, totalvitdmodel5)

# save model results

totalvitdmodelfixed<-as.data.frame((summary(totalvitdmodel)$coefficients))
totalvitdmodelfixed<-setDT(totalvitdmodelfixed, keep.rownames = TRUE)[]

write.csv(totalvitdmodelfixed, file="Results/TotalVitDModelCoatFixed210720.csv")

totalvitdmodelrandom<-as.data.frame(VarCorr(totalvitdmodel))

write.csv(totalvitdmodelrandom, file="Results/TotalVitDModelCoatRandom210720.csv")


## D2

ggplot(VITD,aes(x=CoatBinF,y=X25OHD2)) + geom_boxplot(outlier.colour=NA,notch=FALSE)  + xlab("Coat colour") + ylab("25(OH)D2") +
  theme_classic() + theme(legend.position="none") + geom_jitter() +
  theme(axis.title.x=element_text(face="bold", colour="grey16", size=12)) + 
  theme(axis.title.y=element_text(face="bold", colour="grey16", size=12)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=14, colour="grey26")) + 
  scale_x_discrete(limit = c("1", "2"), labels = c("Dark","Light"))

# model

d2model <- lmer(X25OHD2~Sex+Age+AgeSq+CoatBinF+YearF + (1|ID) +(1|BirthYearF),
                data=VITD)

plot(d2model)
qqnorm(resid(d2model))

summary(d2model)

# sig of sex - x=0.3594      df=1     p=0.5488

d2model2<-update(d2model, ~ . -Sex)
anova(d2model, d2model2)

# sig of age sq - x=178.17      df=1  p< 2.2e-16

d2model3<-update(d2model, ~ . -AgeSq)
anova(d2model, d2model3)

# sig of coat bin - x=4.4193      df=1    p=0.03553

d2model4<-update(d2model, ~ . -CoatBinF)
anova(d2model, d2model4)

# sig of year - x=152.61      df=5  p< 2.2e-16

d2model5<-update(d2model, ~ . -YearF)
anova(d2model, d2model5)

# save model results

d2modelfixed<-as.data.frame((summary(d2model)$coefficients))
d2modelfixed<-setDT(d2modelfixed, keep.rownames = TRUE)[]

write.csv(d2modelfixed, file="Results/D2ModelCoatFixed210720.csv")

d2modelrandom<-as.data.frame(VarCorr(d2model))

write.csv(d2modelrandom, file="Results/D2ModelCoatRandom210720.csv")


## D3

ggplot(VITD,aes(x=CoatBinF,y=X25OHD3)) + geom_boxplot(outlier.colour=NA,notch=FALSE)  + xlab("Coat colour") + ylab("25(OH)D3") +
  theme_classic() + theme(legend.position="none") + geom_jitter() +
  theme(axis.title.x=element_text(face="bold", colour="grey16", size=12)) + 
  theme(axis.title.y=element_text(face="bold", colour="grey16", size=12)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=14, colour="grey26")) + 
  scale_x_discrete(limit = c("1", "2"), labels = c("Dark","Light"))

# model

d3model <- lmer(X25OHD3~Sex+Age+AgeSq+CoatBinF+YearF + (1|ID) +(1|BirthYearF),
                data=VITD)

plot(d3model)
qqnorm(resid(d3model))

summary(d3model)

# sig of sex - x=3.6989     df=1    p=0.05445

d3model2<-update(d3model, ~ . -Sex)
anova(d3model, d3model2)

# sig of age sq - x=184.43      df=1  p< 2.2e-16

d3model3<-update(d3model, ~ . -AgeSq)
anova(d3model, d3model3)

# sig of coat bin - x=28.829      df=1  p=7.904e-08 

d3model4<-update(d3model, ~ . -CoatBinF)
anova(d3model, d3model4)

# sig of year - x=129.19      df=5  p< 2.2e-16

d3model5<-update(d3model, ~ . -YearF)
anova(d3model, d3model5)

# save model results

d3modelfixed<-as.data.frame((summary(d3model)$coefficients))
d3modelfixed<-setDT(d3modelfixed, keep.rownames = TRUE)[]

write.csv(d3modelfixed, file="Results/D3ModelCoatFixed210720.csv")

d3modelrandom<-as.data.frame(VarCorr(d3model))

write.csv(d3modelrandom, file="Results/D3ModelCoatRandom210720.csv")

## Now use AgeF to estimate increase between lambs and yearlings

totalagemodel <- lmer(Total.25D~Sex+AgeF+CoatBinF+YearF + (1|ID) +(1|BirthYearF),
                      data=VITD)
summary(totalagemodel)

d2agemodel <- lmer(X25OHD2~Sex+AgeF+CoatBinF+YearF + (1|ID) +(1|BirthYearF),
                   data=VITD)
summary(d2agemodel)

d3agemodel <- lmer(X25OHD3~Sex+AgeF+CoatBinF+YearF + (1|ID) +(1|BirthYearF),
                   data=VITD)
summary(d3agemodel)


# # -----------------------------------------------------------------------


# Survival models - lambs


# # -----------------------------------------------------------------------

# Make lamb survival model dataset - n=520

LAMBSSURV<-droplevels(subset(VITD, Age=="0" & !is.na(Survival) & !is.na(Weight))) 

# 2011 and 2012 very skewed survival

table(LAMBSSURV$YearF,LAMBSSURV$Survival) 
table(LAMBSSURV$Sex,LAMBSSURV$Survival)

# plot the raw data

ggplot(LAMBSSURV,aes(x=SurvivalF,y=X25OHD2)) + geom_boxplot(outlier.colour=NA,notch=FALSE)  + xlab("Survival") + ylab("25(OH)D2") +
  theme_classic() + theme(legend.position="none") + geom_jitter() +
  theme(axis.title.x=element_text(face="bold", colour="grey16", size=12)) + 
  theme(axis.title.y=element_text(face="bold", colour="grey16", size=12)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=14, colour="grey26")) 

ggplot(LAMBSSURV,aes(x=SurvivalF,y=X25OHD3)) + geom_boxplot(outlier.colour=NA,notch=FALSE)  + xlab("Survival") + ylab("25(OH)D3") +
  theme_classic() + theme(legend.position="none") + geom_jitter() +
  theme(axis.title.x=element_text(face="bold", colour="grey16", size=12)) + 
  theme(axis.title.y=element_text(face="bold", colour="grey16", size=12)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=14, colour="grey26")) 

ggplot(LAMBSSURV,aes(x=SurvivalF,y=Total.25D)) + geom_boxplot(outlier.colour=NA,notch=FALSE)  + xlab("Survival") + ylab("Total vitD") +
  theme_classic() + theme(legend.position="none") + geom_jitter() +
  theme(axis.title.x=element_text(face="bold", colour="grey16", size=12)) + 
  theme(axis.title.y=element_text(face="bold", colour="grey16", size=12)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=14, colour="grey26")) 

# scale continuous predictors to help glm convergence

LAMBSSURV$WeightSq <- LAMBSSURV$Weight * LAMBSSURV$Weight
LAMBSSURV$cWeight <- (LAMBSSURV$Weight - mean(LAMBSSURV$Weight, na.rm=TRUE)) / sd(LAMBSSURV$Weight, na.rm=TRUE)

LAMBSSURV$D2sq <- LAMBSSURV$X25OHD2 * LAMBSSURV$X25OHD2
LAMBSSURV$cD2 <- (LAMBSSURV$X25OHD2 - mean(LAMBSSURV$X25OHD2, na.rm=TRUE)) / sd(LAMBSSURV$X25OHD2, na.rm=TRUE)
LAMBSSURV$cD2sq <- (LAMBSSURV$D2sq - mean(LAMBSSURV$D2sq, na.rm=TRUE)) / sd(LAMBSSURV$D2sq, na.rm=TRUE)

LAMBSSURV$D3sq <- LAMBSSURV$X25OHD3 * LAMBSSURV$X25OHD3
LAMBSSURV$cD3 <- (LAMBSSURV$X25OHD3 - mean(LAMBSSURV$X25OHD3, na.rm=TRUE)) / sd(LAMBSSURV$X25OHD3, na.rm=TRUE)
LAMBSSURV$cD3sq <- (LAMBSSURV$D3sq - mean(LAMBSSURV$D3sq, na.rm=TRUE)) / sd(LAMBSSURV$D3sq, na.rm=TRUE)

LAMBSSURV$TotalDsq <- LAMBSSURV$Total.25D * LAMBSSURV$Total.25D
LAMBSSURV$cTotalD <- (LAMBSSURV$Total.25D - mean(LAMBSSURV$Total.25D, na.rm=TRUE)) / sd(LAMBSSURV$Total.25D, na.rm=TRUE)
LAMBSSURV$cTotalDsq <- (LAMBSSURV$TotalDsq - mean(LAMBSSURV$TotalDsq, na.rm=TRUE)) / sd(LAMBSSURV$TotalDsq, na.rm=TRUE)

# base model

surv.lambs.base<-glm(Survival ~ Sex + CoatBinF + cWeight + YearF, data = LAMBSSURV, family = binomial)
summary(surv.lambs.base)
drop1(surv.lambs.base,test="Chisq")

# compare with d2, d3, tot vitd, d2+d3 + year/sex interactions

surv.lambs.d2<-glm(Survival ~ Sex + CoatBinF + cWeight + YearF +  cD2, data = LAMBSSURV, family = binomial)
surv.lambs.d3<-glm(Survival ~ Sex + CoatBinF + cWeight + YearF +  cD3, data = LAMBSSURV, family = binomial)
surv.lambs.d2d3<-glm(Survival ~ Sex + CoatBinF + cWeight  + YearF +  cD2 + cD3, data = LAMBSSURV, family = binomial)
surv.lambs.dtot<-glm(Survival ~ Sex + CoatBinF + cWeight  + YearF +  cTotalD, data = LAMBSSURV, family = binomial)
surv.lambs.d2sq<-glm(Survival ~ Sex + CoatBinF + cWeight + YearF +  cD2 + cD2sq, data = LAMBSSURV, family = binomial)
surv.lambs.d3sq<-glm(Survival ~ Sex + CoatBinF + cWeight + YearF +  cD3 + cD3sq, data = LAMBSSURV, family = binomial)
surv.lambs.d2d3sq<-glm(Survival ~ Sex + CoatBinF + cWeight  + YearF +  cD2 + cD3 + cD2sq + cD3sq, data = LAMBSSURV, family = binomial)
surv.lambs.dtotsq<-glm(Survival ~ Sex + CoatBinF + cWeight  + YearF + cTotalD + cTotalDsq, data = LAMBSSURV, family = binomial)
surv.lambs.d2yr<-glm(Survival ~ Sex + CoatBinF + cWeight  + YearF + cD2 + YearF*cD2, data = LAMBSSURV, family = binomial)
surv.lambs.d3yr<-glm(Survival ~ Sex + CoatBinF + cWeight  + YearF + cD3 + YearF*cD3, data = LAMBSSURV, family = binomial)
surv.lambs.d2d3yr<-glm(Survival ~ Sex + CoatBinF + cWeight  + YearF + cD2 + cD3 + YearF*cD2 + YearF*cD3, data = LAMBSSURV, family = binomial)
surv.lambs.dtotyr<-glm(Survival ~ Sex + CoatBinF + cWeight  + YearF + cTotalD + YearF*cTotalD, data = LAMBSSURV, family = binomial)
surv.lambs.d2s<-glm(Survival ~ Sex + CoatBinF + cWeight + YearF + cD2 + Sex*cD2, data = LAMBSSURV, family = binomial)
surv.lambs.d3s<-glm(Survival ~ Sex + CoatBinF + cWeight  + YearF + cD3 + Sex*cD3, data = LAMBSSURV, family = binomial)
surv.lambs.d2d3s<-glm(Survival ~ Sex + CoatBinF + cWeight + YearF + cD2 + cD3 + Sex*cD2 + Sex*cD3, data = LAMBSSURV, family = binomial)
surv.lambs.dtots<-glm(Survival ~ Sex + CoatBinF + cWeight + YearF + cTotalD + Sex*cTotalD, data = LAMBSSURV, family = binomial)

models_list_lambsurv <- list(surv.lambs.base,surv.lambs.d2,surv.lambs.d3,surv.lambs.d2d3,surv.lambs.dtot,
                             surv.lambs.d2sq,surv.lambs.d3sq,surv.lambs.d2d3sq,surv.lambs.dtotsq,
                             surv.lambs.d2yr,surv.lambs.d3yr,surv.lambs.d2d3yr,surv.lambs.dtotyr,
                             surv.lambs.d2s, surv.lambs.d3s, surv.lambs.d2d3s, surv.lambs.dtots)
model_names_lambsurv <- c("base", "d2", "d3", "d2d3", "dtot", "d2sq", "d3sq", "d2d3sq", "dtotsq",
                          "d2yr", "d3yr", "d2d3yr", "dtotyr", "d2sx", "d3sx", "d2d3sx", "dtotsx")

# check AICc

#aictab(models_list_lambsurv, model_names_lambsurv)

#evidence(aictab(cand.set=models_list_lambsurv, modnames=model_names_lambsurv))

# check AIC 

aictab(models_list_lambsurv, model_names_lambsurv, second.ord=FALSE)

# best model within 2 AIC of base model 

summary(surv.lambs.base)

plot(surv.lambs.base)

# results are consistent without 2011 and 2012 (skewed survival years)

# # -----------------------------------------------------------------------


# Survival models - adults


# # -----------------------------------------------------------------------

# Make model dataset - 851 measures from 415 sheep

ADULTSURV<-droplevels(subset(VITD, Age>0 & !is.na(Survival) & !is.na(Weight)))

# plot raw data - blue boxes are individuals that survived, 1 are females and 2 are males 

ggplot(ADULTSURV,aes(x=Sex,y=X25OHD2, fill=SurvivalF)) + geom_boxplot(outlier.colour=NA,notch=FALSE)  + xlab("Sex") + ylab("25(OH)D2") +
  theme_classic() + geom_jitter(alpha=0.2) +
  theme(axis.title.x=element_text(face="bold", colour="grey16", size=12)) + 
  theme(axis.title.y=element_text(face="bold", colour="grey16", size=12)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=14, colour="grey26")) 

ggplot(ADULTSURV,aes(x=Sex,y=X25OHD3, fill=SurvivalF)) + geom_boxplot(outlier.colour=NA,notch=FALSE)  + xlab("Sex") + ylab("25(OH)D3") +
  theme_classic() + geom_jitter(alpha=0.2) +
  theme(axis.title.x=element_text(face="bold", colour="grey16", size=12)) + 
  theme(axis.title.y=element_text(face="bold", colour="grey16", size=12)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=14, colour="grey26")) 

ggplot(ADULTSURV,aes(x=Sex,y=Total.25D, fill=SurvivalF)) + geom_boxplot(outlier.colour=NA,notch=FALSE)  + xlab("Sex") + ylab("Total vitD") +
  theme_classic() + geom_jitter(alpha=0.2) +
  theme(axis.title.x=element_text(face="bold", colour="grey16", size=12)) + 
  theme(axis.title.y=element_text(face="bold", colour="grey16", size=12)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=14, colour="grey26")) 

# in adults - 2012 and 2016 v skewed 

table(ADULTSURV$YearF,ADULTSURV$Survival)
table(ADULTSURV$Sex,ADULTSURV$Survival)


# re-scale variables to help model convergence

ADULTSURV$cWeight <- (ADULTSURV$Weight - mean(ADULTSURV$Weight, na.rm=TRUE)) / sd(ADULTSURV$Weight, na.rm=TRUE)

ADULTSURV$D2sq <- ADULTSURV$X25OHD2 * ADULTSURV$X25OHD2
ADULTSURV$cD2 <- (ADULTSURV$X25OHD2 - mean(ADULTSURV$X25OHD2, na.rm=TRUE)) / sd(ADULTSURV$X25OHD2, na.rm=TRUE)
ADULTSURV$cD2sq <- (ADULTSURV$D2sq - mean(ADULTSURV$D2sq, na.rm=TRUE)) / sd(ADULTSURV$D2sq, na.rm=TRUE)

ADULTSURV$D3sq <- ADULTSURV$X25OHD3 * ADULTSURV$X25OHD3
ADULTSURV$cD3 <- (ADULTSURV$X25OHD3 - mean(ADULTSURV$X25OHD3, na.rm=TRUE)) / sd(ADULTSURV$X25OHD3, na.rm=TRUE)
ADULTSURV$cD3sq <- (ADULTSURV$D3sq - mean(ADULTSURV$D3sq, na.rm=TRUE)) / sd(ADULTSURV$D3sq, na.rm=TRUE)

ADULTSURV$TotalDsq <- ADULTSURV$Total.25D * ADULTSURV$Total.25D
ADULTSURV$cTotalD <- (ADULTSURV$Total.25D - mean(ADULTSURV$Total.25D, na.rm=TRUE)) / sd(ADULTSURV$Total.25D, na.rm=TRUE)
ADULTSURV$cTotalDsq <- (ADULTSURV$TotalDsq - mean(ADULTSURV$TotalDsq, na.rm=TRUE)) / sd(ADULTSURV$TotalDsq, na.rm=TRUE)

# compare base model with models including d2, d3, tot vitd, d2+d3 + year/sex/age interactions
# note there are repeat measures in this dataset, but running as a GLMM gave warnings of singularity in lme4
# as ID explained 0 variance in survival - run as glm but results are consistent as a GLMM 

surv.adults.base <- glm(Survival ~ Sex + AgeGroupF + CoatBinF + Sex*AgeGroupF + cWeight + YearF, 
                   data = ADULTSURV, family=binomial)
surv.adults.d2 <- glm(Survival ~ Sex + AgeGroupF + CoatBinF + Sex*AgeGroupF + cWeight + YearF + cD2, 
                 data = ADULTSURV, family=binomial)
surv.adults.d3 <- glm(Survival ~ Sex + AgeGroupF + CoatBinF + Sex*AgeGroupF + cWeight + YearF + cD3, 
                 data = ADULTSURV, family=binomial)
surv.adults.dtot <- glm(Survival ~ Sex + AgeGroupF + CoatBinF + Sex*AgeGroupF + cWeight + YearF + cTotalD,
                  data = ADULTSURV, family=binomial)
surv.adults.d2d3 <- glm(Survival ~ Sex + AgeGroupF + CoatBinF + Sex*AgeGroupF + cWeight + YearF + cD2 + cD3, 
                   data = ADULTSURV, family=binomial)
surv.adults.d2sq <- glm(Survival ~ Sex + AgeGroupF + CoatBinF + Sex*AgeGroupF + cWeight + YearF + cD2 + cD2sq, 
                 data = ADULTSURV, family=binomial)
surv.adults.d3sq <- glm(Survival ~ Sex + AgeGroupF + CoatBinF + Sex*AgeGroupF + cWeight + YearF + cD3 + cD3sq, 
                 data = ADULTSURV, family=binomial)
surv.adults.dtotsq <- glm(Survival ~ Sex + AgeGroupF + CoatBinF + Sex*AgeGroupF + cWeight + YearF + cTotalD + cTotalDsq,
                  data = ADULTSURV, family=binomial)
surv.adults.d2d3sq <- glm(Survival ~ Sex + AgeGroupF + CoatBinF + Sex*AgeGroupF + cWeight + YearF + cD2 + cD3 + cD2sq + cD3sq, 
                   data = ADULTSURV, family=binomial)
surv.adults.d2yr <- glm(Survival ~ Sex + AgeGroupF + CoatBinF + Sex*AgeGroupF + cWeight + YearF + cD2 + cD2*YearF,
                   data = ADULTSURV, family=binomial)
surv.adults.d3yr <- glm(Survival ~ Sex + AgeGroupF + CoatBinF + Sex*AgeGroupF + cWeight + YearF + cD3 + cD3*YearF, 
                   data = ADULTSURV, family=binomial)
surv.adults.dtotyr <- glm(Survival ~ Sex + AgeGroupF + CoatBinF + Sex*AgeGroupF + cWeight + YearF + cTotalD + cTotalD*YearF,
                  data = ADULTSURV, family=binomial)
surv.adults.d2d3yr <- glm(Survival ~ Sex + AgeGroupF + CoatBinF + Sex*AgeGroupF + cWeight + YearF + cD2 + cD3 +
                       cD2*YearF + cD3*YearF,
                     data = ADULTSURV, family=binomial)
surv.adults.d2s <- glm(Survival ~ Sex + AgeGroupF + CoatBinF + Sex*AgeGroupF + cWeight + YearF + cD2 + cD2*Sex,
                 data = ADULTSURV, family=binomial)
surv.adults.d2a <- glm(Survival ~ Sex + AgeGroupF + CoatBinF + Sex*AgeGroupF + cWeight + YearF + cD2 + cD2*AgeGroupF,
                 data = ADULTSURV, family=binomial)
surv.adults.d3s <- glm(Survival ~ Sex + AgeGroupF + CoatBinF + Sex*AgeGroupF + cWeight + YearF +
                   cD3 + cD3*Sex,
                 data = ADULTSURV, family=binomial)
surv.adults.d3a <- glm(Survival ~ Sex + AgeGroupF + CoatBinF + Sex*AgeGroupF + cWeight + YearF +
                   cD3 + cD3*AgeGroupF, 
                 data = ADULTSURV, family=binomial)
surv.adults.dtots <- glm(Survival ~ Sex + AgeGroupF + CoatBinF + Sex*AgeGroupF + cWeight + YearF +
                    cTotalD + cTotalD*Sex, 
                  data = ADULTSURV, family=binomial)
surv.adults.dtota <- glm(Survival ~ Sex + AgeGroupF + CoatBinF + Sex*AgeGroupF + cWeight + YearF +
                    cTotalD + cTotalD*AgeGroupF, 
                  data = ADULTSURV, family=binomial)
surv.adults.d2d3s <- glm(Survival ~ Sex + AgeGroupF + CoatBinF + Sex*AgeGroupF + cWeight + YearF + cD2 + cD3 +
                     cD2*Sex + cD3*Sex, 
                   data = ADULTSURV, family=binomial)
surv.adults.d2d3a <- glm(Survival ~ Sex + AgeGroupF + CoatBinF + Sex*AgeGroupF + cWeight + YearF + cD2 + cD3 +
                     cD2*AgeGroupF + cD3*AgeGroupF, 
                   data = ADULTSURV, family=binomial)

model_list_adultsurv <- list(surv.adults.base,surv.adults.d2,surv.adults.d3,surv.adults.dtot,surv.adults.d2d3,surv.adults.d2sq,
                          surv.adults.d3sq,surv.adults.dtotsq,surv.adults.d2d3sq,surv.adults.d2yr,surv.adults.d3yr,surv.adults.dtotyr,
                          surv.adults.d2d3yr,surv.adults.d2s,surv.adults.d2a,surv.adults.d3s,surv.adults.d3a,surv.adults.dtots,
                          surv.adults.dtota,surv.adults.d2d3s,surv.adults.d2d3a)
model_names_adultsurv <- c("base", "d2", "d3", "dtot", "d2d3", "d2sq", "d3sq", "dtotsq", "d2d3sq", "d2yr", "d3yr", "dtotyr", "d2d3yr",
                        "d2sx", "d2age", "d3sx", "d3age", "dtotsx", "dtotage", "d2d3sx", "d2d3age")

#aictab(model_list_adultsurv, model_names_adultsurv)

# compare AICs - d2 by sex is best fit (by greater than 2 AIC to next best fit of d2 alone) 

aictab(model_list_adultsurv, model_names_adultsurv, second.ord=FALSE) 

# best model - d2 by sex interaction

plot(surv.adults.d2s)

summary(surv.adults.d2s)

# save model results

write.csv(summary(surv.adults.d2s)$coef,"Results/AdultSurvivalD2SexAllFixed220720.csv")

# significance of fixed effects

drop1(surv.adults.d2s, test="Chisq")

####### Plot the d2 by sex model predictions - figure for ms in graph script

# model

surv.adults.d2s <- glm(Survival ~ Sex + AgeGroupF + CoatBinF + Sex*AgeGroupF + cWeight + YearF + cD2 + cD2*Sex,
               data = ADULTSURV, family=binomial)

summary(surv.adults.d2s)

# make dataframe to predict model values over

preddf2 <- data.frame(Sex="1",
                      AgeGroupF="1",
                      YearF="2011",
                      CoatBinF="1",
                      cWeight=mean(ADULTSURV$cWeight),
                      cD2=unique(sort(ADULTSURV$cD2)))

preddf3 <- data.frame(Sex="2",
                      AgeGroupF="1",
                      YearF="2011",
                      CoatBinF="1",
                      cWeight=mean(ADULTSURV$cWeight),
                      cD2=unique(sort(ADULTSURV$cD2)))


# predict and add predictions to dataframe

predd2a<- predict(surv.adults.d2s, preddf2, re.form=NA, type='response', se.fit=T)

preddf2$pred <- predd2a$fit
preddf2$upper <- predd2a$fit + predd2a$se.fit
preddf2$lower <- predd2a$fit - predd2a$se.fit

predd2b<- predict(surv.adults.d2s, preddf3, re.form=NA, type='response', se.fit=T)

preddf3$pred <- predd2b$fit
preddf3$upper <- predd2b$fit + predd2b$se.fit
preddf3$lower <- predd2b$fit - predd2b$se.fit

# plot and save graphs 

wsurvsgraph<-ggplot(ADULTSURV, aes(x=cD2, y=Survival))+ geom_point(shape=21, alpha=0.5, colour="grey40") +
  labs(x=expression("25(OH)D"["2"]~"levels"~"(standardised)"), y="Over-winter survival") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=12)) +
  theme(axis.title.y=element_text(colour="grey6", size=12)) +
  theme(axis.line=element_line(colour="grey6")) +
  theme(text = element_text(size=12, colour="grey26")) + 
  geom_line(data=preddf2, aes(x=cD2, y=pred), colour="aquamarine4", size=1) + 
  geom_line(data=preddf2, aes(x=cD2, y=upper), colour="aquamarine4", linetype="dashed", size=0.5) +
  geom_line(data=preddf2, aes(x=cD2, y=lower), colour="aquamarine4", linetype="dashed", size=0.5) + 
  geom_line(data=preddf3, aes(x=cD2, y=pred), colour="#999999", size=1) + 
  geom_line(data=preddf3, aes(x=cD2, y=upper), colour="#999999", linetype="dashed", size=0.5) +
  geom_line(data=preddf3, aes(x=cD2, y=lower), colour="#999999", linetype="dashed", size=0.5)

wsurvsgraph


# Looks like for males it might not be significant...  


# Males only model to see if d2 is significant:

ADULTSURVM <- droplevels(subset(ADULTSURV, Sex==2))

surv.adults.d2.males<-glm(Survival ~ AgeGroupF + CoatBinF + cWeight + YearF + cD2,
                          data = ADULTSURVM, family=binomial)

plot(surv.adults.d2.males)

summary(surv.adults.d2.males)

write.csv(summary(surv.adults.d2.males)$coef,"Results/AdultSurvivalD2MalesFixed220720.csv")

drop1(surv.adults.d2.males, test="Chisq") # d2 not significant

# Females only model to see if d2 significant

ADULTSURVF <- droplevels(subset(ADULTSURV, Sex==1))

surv.adults.d2.females<-glm(Survival ~ AgeGroupF + CoatBinF + cWeight + YearF + cD2,
                            data = ADULTSURVF, family=binomial)

plot(surv.adults.d2.females)

summary(surv.adults.d2.females)

write.csv(summary(surv.adults.d2.females)$coef,"Results/AdultSurvivalD2FemalesFixed220720.csv")

drop1(surv.adults.d2.females, test="Chisq") # d2 is significant for females only

# results are consistent without 2012 and 2016 (skewed survival years)

# # -----------------------------------------------------------------------


# Female lamb breeding success


# # -----------------------------------------------------------------------

# make model subset - n=112

lambffec  <- droplevels(subset(VITD, Age==0 & Sex==1 & !is.na(EweFecundity) & !is.na(Weight)))

# vit d measures by female lamb breeding success (0 didn't lamb/ 1 lambed in their first year)

ggplot(lambffec,aes(x=factor(EweFecundity), y=X25OHD2)) + geom_boxplot(outlier.colour=NA,notch=FALSE)  + xlab("Female lamb ABS") + ylab("25(OH)D2") +
  theme_classic() + geom_jitter() +
  theme(axis.title.x=element_text(face="bold", colour="grey16", size=12)) + 
  theme(axis.title.y=element_text(face="bold", colour="grey16", size=12)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=14, colour="grey26")) 

ggplot(lambffec,aes(x=factor(EweFecundity), y=X25OHD3)) + geom_boxplot(outlier.colour=NA,notch=FALSE)  + xlab("Female lamb ABS") + ylab("25(OH)D3") +
  theme_classic() + geom_jitter() +
  theme(axis.title.x=element_text(face="bold", colour="grey16", size=12)) + 
  theme(axis.title.y=element_text(face="bold", colour="grey16", size=12)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=14, colour="grey26")) 

ggplot(lambffec,aes(x=factor(EweFecundity), y=Total.25D)) + geom_boxplot(outlier.colour=NA,notch=FALSE)  + xlab("Female lamb ABS") + ylab("Total vitamin D") +
  theme_classic() + geom_jitter() +
  theme(axis.title.x=element_text(face="bold", colour="grey16", size=12)) + 
  theme(axis.title.y=element_text(face="bold", colour="grey16", size=12)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=14, colour="grey26")) 

table(lambffec$EweFecundityO)

table(lambffec$EweFecundity)

table(lambffec$EweFecundity, lambffec$YearF)
table(lambffec$Survival, lambffec$YearF)

# models with all years

# rescale variables

lambffec$cWeight <- (lambffec$Weight - mean(lambffec$Weight, na.rm=TRUE)) / sd(lambffec$Weight, na.rm=TRUE)

lambffec$D2sq <- lambffec$X25OHD2 * lambffec$X25OHD2
lambffec$cD2 <- (lambffec$X25OHD2 - mean(lambffec$X25OHD2, na.rm=TRUE)) / sd(lambffec$X25OHD2, na.rm=TRUE)
lambffec$cD2sq <- (lambffec$D2sq - mean(lambffec$D2sq, na.rm=TRUE)) / sd(lambffec$D2sq, na.rm=TRUE)

lambffec$D3sq <- lambffec$X25OHD3 * lambffec$X25OHD3
lambffec$cD3 <- (lambffec$X25OHD3 - mean(lambffec$X25OHD3, na.rm=TRUE)) / sd(lambffec$X25OHD3, na.rm=TRUE)
lambffec$cD3sq <- (lambffec$D3sq - mean(lambffec$D3sq, na.rm=TRUE)) / sd(lambffec$D3sq, na.rm=TRUE)

lambffec$TotalDsq <- lambffec$Total.25D * lambffec$Total.25D
lambffec$cTotalD <- (lambffec$Total.25D - mean(lambffec$Total.25D, na.rm=TRUE)) / sd(lambffec$Total.25D, na.rm=TRUE)
lambffec$cTotalDsq <- (lambffec$TotalDsq - mean(lambffec$TotalDsq, na.rm=TRUE)) / sd(lambffec$TotalDsq, na.rm=TRUE)

# compare base model with models including d2, d3, tot vitd, d2+d3 + year interactions

lfec.base<-glm(EweFecundity ~ CoatBinF + cWeight + YearF, data = lambffec, family = binomial)
lfec.d2<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cD2, data = lambffec, family = binomial)
lfec.d3<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cD3, data = lambffec, family = binomial)
lfec.dtot<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cTotalD, data = lambffec, family = binomial)
lfec.d2d3<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cD2 + cD3, data = lambffec, family = binomial)
lfec.d2.sq<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cD2 + cD2sq, data = lambffec, family = binomial)
lfec.d3.sq<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cD3 + cD3sq, data = lambffec, family = binomial)
lfec.dtot.sq<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cTotalD + cTotalDsq, data = lambffec, family = binomial)
lfec.d2d3.sq<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cD2 + cD3 + cD2sq + cD3sq, data = lambffec, family = binomial)
lfec.d2.yr<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cD2 + cD2*YearF, data = lambffec, family = binomial)
lfec.d3.yr<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cD3 + cD3*YearF, data = lambffec, family = binomial)
lfec.dtot.yr<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cTotalD + cTotalD*YearF, data = lambffec, family = binomial)
lfec.d2d3.yr<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cD2 + cD3 + cD2*YearF + cD2*YearF, data = lambffec, family = binomial)

# compare AIC - best fitting model within 2 AICs of base - accept base as best model

#AIC(lfec.base,lfec.d2,lfec.d3,lfec.dtot, lfec.d2d3, lfec.d2,lfec.d3,lfec.dtot, lfec.d2d3)

models_lambffec <- list(lfec.base,lfec.d2,lfec.d3,lfec.dtot, lfec.d2d3, 
                        lfec.d2.sq,lfec.d3.sq,lfec.dtot.sq, lfec.d2d3.sq,
                        lfec.d2.yr,lfec.d3.yr,lfec.dtot.yr, lfec.d2d3.yr)
model_lambffec_names <- c("base", "d2", "d3", "dtot", "d2d3",
                          "d2sq", "d3sq", "dtotsq", "d2d3sq",
                          "d2yr", "d3yr", "dtotyr", "d2d3yr")

aictab(models_lambffec, model_lambffec_names, second.ord=FALSE) 

summary(lfec.base)
plot(lfec.base)

#### Only 4 lambs with measures in 2011 - use models without 2011 in manuscript
# Model with all years in supp (results consistent)

# Sample size = 108

lambffec2  <- droplevels(subset(VITD, Age==0 & Sex==1 & Year!=2011 & !is.na(EweFecundity) & !is.na(Weight)))

# rescale variables

lambffec2$cWeight <- (lambffec2$Weight - mean(lambffec2$Weight, na.rm=TRUE)) / sd(lambffec2$Weight, na.rm=TRUE)

lambffec2$D2sq <- lambffec2$X25OHD2 * lambffec2$X25OHD2
lambffec2$cD2 <- (lambffec2$X25OHD2 - mean(lambffec2$X25OHD2, na.rm=TRUE)) / sd(lambffec2$X25OHD2, na.rm=TRUE)
lambffec2$cD2sq <- (lambffec2$D2sq - mean(lambffec2$D2sq, na.rm=TRUE)) / sd(lambffec2$D2sq, na.rm=TRUE)

lambffec2$D3sq <- lambffec2$X25OHD3 * lambffec2$X25OHD3
lambffec2$cD3 <- (lambffec2$X25OHD3 - mean(lambffec2$X25OHD3, na.rm=TRUE)) / sd(lambffec2$X25OHD3, na.rm=TRUE)
lambffec2$cD3sq <- (lambffec2$D3sq - mean(lambffec2$D3sq, na.rm=TRUE)) / sd(lambffec2$D3sq, na.rm=TRUE)

lambffec2$TotalDsq <- lambffec2$Total.25D * lambffec2$Total.25D
lambffec2$cTotalD <- (lambffec2$Total.25D - mean(lambffec2$Total.25D, na.rm=TRUE)) / sd(lambffec2$Total.25D, na.rm=TRUE)
lambffec2$cTotalDsq <- (lambffec2$TotalDsq - mean(lambffec2$TotalDsq, na.rm=TRUE)) / sd(lambffec2$TotalDsq, na.rm=TRUE)

# compare base model with models including d2, d3, tot vitd, d2+d3 + year interactions

lfec.base.2<-glm(EweFecundity ~ CoatBinF + cWeight + YearF, data = lambffec2, family = binomial)
lfec.d2.2<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cD2, data = lambffec2, family = binomial)
lfec.d3.2<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cD3, data = lambffec2, family = binomial)
lfec.dtot.2<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cTotalD, data = lambffec2, family = binomial)
lfec.d2d3.2<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cD2 + cD3, data = lambffec2, family = binomial)
lfec.d2.sq.2<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cD2 + cD2sq, data = lambffec2, family = binomial)
lfec.d3.sq.2<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cD3 + cD3sq, data = lambffec2, family = binomial)
lfec.dtot.sq.2<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cTotalD + cTotalDsq, data = lambffec2, family = binomial)
lfec.d2d3.sq.2<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cD2 + cD3 + cD2sq + cD3sq, data = lambffec2, family = binomial)
lfec.d2.yr.2<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cD2 + cD2*YearF, data = lambffec2, family = binomial)
lfec.d3.yr.2<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cD3 + cD3*YearF, data = lambffec2, family = binomial)
lfec.dtot.yr.2<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cTotalD + cTotalD*YearF, data = lambffec2, family = binomial)
lfec.d2d3.yr.2<-glm(EweFecundity ~ CoatBinF + cWeight + YearF + cD2 + cD3 + cD2*YearF + cD2*YearF, data = lambffec2, family = binomial)

#AIC(lfec.base,lfec.d2,lfec.d3,lfec.dtot, lfec.d2d3, lfec.d2,lfec.d3,lfec.dtot, lfec.d2d3)

models_lambffec2 <- list(lfec.base.2,lfec.d2.2,lfec.d3.2,lfec.dtot.2, lfec.d2d3.2, 
                        lfec.d2.sq.2,lfec.d3.sq.2,lfec.dtot.sq.2, lfec.d2d3.sq.2,
                        lfec.d2.yr.2,lfec.d3.yr.2,lfec.dtot.yr.2, lfec.d2d3.yr.2)
model_lambffec_names2 <- c("base", "d2", "d3", "dtot", "d2d3",
                          "d2sq", "d3sq", "dtotsq", "d2d3sq",
                          "d2yr", "d3yr", "dtotyr", "d2d3yr")

#aictab(models_lambffec, model_lambffec_names) 

# d3 lowest AIC, then tot d but both within 2 AIC from base

aictab(models_lambffec2, model_lambffec_names2, second.ord=FALSE) 

# With or without 2011 measures - base model is the best model

summary(lfec.base.2)
plot(lfec.base.2)


# # -----------------------------------------------------------------------


# Adult female fecundity


# # -----------------------------------------------------------------------

# make model subset - n=578 measures of 245 females

adultffec  <- droplevels(subset(VITD, Age>0 & !is.na(EweFecundityO) & !is.na(Weight))) # 578 obs

# plot ewe fecundity vs vitamin d measures

ggplot(adultffec,aes(x=factor(EweFecundity), y=X25OHD2)) + geom_boxplot(outlier.colour=NA,notch=FALSE)  + xlab("Ewe fecundity") + ylab("25(OH)D2") +
  theme_classic() + geom_jitter() +
  theme(axis.title.x=element_text(face="bold", colour="grey16", size=12)) + 
  theme(axis.title.y=element_text(face="bold", colour="grey16", size=12)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=14, colour="grey26")) 

ggplot(adultffec,aes(x=factor(EweFecundity), y=X25OHD3)) + geom_boxplot(outlier.colour=NA,notch=FALSE)  + xlab("Ewe fecundity") + ylab("25(OH)D3") +
  theme_classic() + geom_jitter() +
  theme(axis.title.x=element_text(face="bold", colour="grey16", size=12)) + 
  theme(axis.title.y=element_text(face="bold", colour="grey16", size=12)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=14, colour="grey26")) 

ggplot(adultffec,aes(x=factor(EweFecundity), y=Total.25D)) + geom_boxplot(outlier.colour=NA,notch=FALSE)  + xlab("Ewe fecundity") + ylab("Total vitamin D") +
  theme_classic() + geom_jitter() +
  theme(axis.title.x=element_text(face="bold", colour="grey16", size=12)) + 
  theme(axis.title.y=element_text(face="bold", colour="grey16", size=12)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=14, colour="grey26")) 

# re-scale

adultffec$cWeight <- (adultffec$Weight - mean(adultffec$Weight, na.rm=TRUE)) / sd(adultffec$Weight, na.rm=TRUE)

adultffec$D2sq <- adultffec$X25OHD2 * adultffec$X25OHD2
adultffec$cD2 <- (adultffec$X25OHD2 - mean(adultffec$X25OHD2, na.rm=TRUE)) / sd(adultffec$X25OHD2, na.rm=TRUE)
adultffec$cD2sq <- (adultffec$D2sq - mean(adultffec$D2sq, na.rm=TRUE)) / sd(adultffec$D2sq, na.rm=TRUE)

adultffec$D3sq <- adultffec$X25OHD3 * adultffec$X25OHD3
adultffec$cD3 <- (adultffec$X25OHD3 - mean(adultffec$X25OHD3, na.rm=TRUE)) / sd(adultffec$X25OHD3, na.rm=TRUE)
adultffec$cD3sq <- (adultffec$D3sq - mean(adultffec$D3sq, na.rm=TRUE)) / sd(adultffec$D3sq, na.rm=TRUE)

adultffec$TotalDsq <- adultffec$Total.25D * adultffec$Total.25D
adultffec$cTotalD <- (adultffec$Total.25D - mean(adultffec$Total.25D, na.rm=TRUE)) / sd(adultffec$Total.25D, na.rm=TRUE)
adultffec$cTotalDsq <- (adultffec$TotalDsq - mean(adultffec$TotalDsq, na.rm=TRUE)) / sd(adultffec$TotalDsq, na.rm=TRUE)

# check if years okay

table(adultffec$YearF, adultffec$EweFecundityO) # no need to remove any years


# compare base model with models including d2, d3, tot vitd, d2+d3 + age/year interactions

adultffec.base <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + YearF + (1|ID),
                       Hess=TRUE, nAGQ=10, data=adultffec)

adultffec.d2 <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + YearF + cD2 + (1|ID),
                     Hess=TRUE, nAGQ=10, data=adultffec)

adultffec.d3 <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + YearF +  cD3 + (1|ID),
                     Hess=TRUE, nAGQ=10, data=adultffec)

adultffec.dtot <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + YearF + cTotalD + (1|ID),
                       Hess=TRUE, nAGQ=10, data=adultffec)

adultffec.d2d3 <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + YearF + cD2 + cD3 + (1|ID),
                       Hess=TRUE, nAGQ=10, data=adultffec)

adultffec.d2.sq <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + YearF + cD2 + cD2sq + (1|ID),
                     Hess=TRUE, nAGQ=10, data=adultffec)

adultffec.d3.sq <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + YearF +  cD3 + cD3sq + (1|ID),
                     Hess=TRUE, nAGQ=10, data=adultffec)

adultffec.dtot.sq <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + YearF + cTotalD + cTotalDsq + (1|ID),
                       Hess=TRUE, nAGQ=10, data=adultffec)

adultffec.d2d3.sq <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + YearF + cD2 + cD3 + cD2sq + cD3sq + (1|ID),
                       Hess=TRUE, nAGQ=10, data=adultffec)

adultffec.d2yr <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + YearF + cD2 + cD2*YearF + (1|ID),
                       Hess=TRUE, nAGQ=10, data=adultffec)

adultffec.d3yr <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + YearF + cD3 + cD3*YearF + (1|ID),
                       Hess=TRUE, nAGQ=10, data=adultffec)

adultffec.dtotyr <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + YearF + cTotalD + cTotalD*YearF + (1|ID),
                         Hess=TRUE, nAGQ=10, data=adultffec)

adultffec.d2d3yr <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + YearF +
                           cD2 + cD3 + cD2*YearF + cD3*YearF + (1|ID),
                         Hess=TRUE, nAGQ=10, data=adultffec)

adultffec.d2age <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + YearF + cD2 + cD2*AgeGroupF + (1|ID),
                       Hess=TRUE, nAGQ=10, data=adultffec)

adultffec.d3age <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + YearF + cD3 + cD3*AgeGroupF + (1|ID),
                       Hess=TRUE, nAGQ=10, data=adultffec)

adultffec.dtotage <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + YearF + cTotalD + cTotalD*AgeGroupF + (1|ID),
                         Hess=TRUE, nAGQ=10, data=adultffec)

adultffec.d2d3age <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + YearF +
                           cD2 + cD3 + cD2*AgeGroupF + cD3*AgeGroupF + (1|ID),
                         Hess=TRUE, nAGQ=10, data=adultffec)


adultffec_model_list <- list(adultffec.base, adultffec.d2, adultffec.d3, adultffec.dtot, adultffec.d2d3,
                             adultffec.d2.sq, adultffec.d3.sq, adultffec.dtot.sq, adultffec.d2d3.sq,
                             adultffec.d2yr, adultffec.d3yr, adultffec.dtotyr, adultffec.d2d3yr,
                             adultffec.d2age, adultffec.d3age, adultffec.dtotage, adultffec.d2d3age)
adultffec_model_names <- c("base", "d2", "d3", "dtot", "d2d3", "d2sq", "d3sq", "dtotsq", "d2d3sq",
                           "d2yr", "d3yr", "dtotyr", "d2d3yr", "d2age", "d3age", "dtotage", "d2d3age")

# aictab(adultffec_model_list, adultffec_model_names)

aictab(adultffec_model_list, adultffec_model_names, second.ord=FALSE)

# Two models with the lowest AIC (delta 0.12 difference): 
# the model with a total vit D by year interaction
# and the model with a d3 by year interaction have the lowest AIC 

## model with total vit d by year interaction

adultffec.dtotyr <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + YearF + cTotalD + cTotalD*YearF + (1|ID),
                         Hess=TRUE, nAGQ=10, data=adultffec)

summary(adultffec.dtotyr)

# save model results

adultffec.dtotyr.model<-as.data.frame((summary(adultffec.dtotyr)$coefficients))
adultffec.dtotyr.model<-setDT(adultffec.dtotyr.model, keep.rownames = TRUE)[]

write.csv(adultffec.dtotyr.model, file="Results/AdultFFecTotDYearFixed220720.csv")

# check significance of fixed effects

# sig of age - x=7.881  df=2    p=0.01944 

adultffec.dtotyr.a <- clmm(EweFecundityO ~ CoatBinF + cWeight + YearF + cTotalD + cTotalD*YearF + (1|ID),
                           Hess=TRUE, nAGQ=10, data=adultffec)

anova(adultffec.dtotyr, adultffec.dtotyr.a)

# sig of coat - x=0.1142  df=1     p=0.7354

adultffec.dtotyr.c <- clmm(EweFecundityO ~ AgeGroupF + cWeight + YearF + cTotalD + cTotalD*YearF + (1|ID),
                         Hess=TRUE, nAGQ=10, data=adultffec)

anova(adultffec.dtotyr, adultffec.dtotyr.c)

# sig of weight - x=3.5271  df=1    p=0.06037

adultffec.dtotyr.w <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + YearF + cTotalD + cTotalD*YearF + (1|ID),
                         Hess=TRUE, nAGQ=10, data=adultffec)

anova(adultffec.dtotyr, adultffec.dtotyr.w)

# sig of dtotyr - x=19.304  df=5   p=0.001687

adultffec.dtotyr.d <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + YearF + cTotalD + (1|ID),
                         Hess=TRUE, nAGQ=10, data=adultffec)

anova(adultffec.dtotyr, adultffec.dtotyr.d)


## model with year by d3 interaction

adultffec.d3yr <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + YearF + cD3 + cD3*YearF + (1|ID),
                       Hess=TRUE, nAGQ=10, data=adultffec)

summary(adultffec.d3yr)

# save model results

adultffec.d3yr.model<-as.data.frame((summary(adultffec.d3yr)$coefficients))
adultffec.d3yr.model<-setDT(adultffec.d3yr.model, keep.rownames = TRUE)[]

write.csv(adultffec.d3yr.model, file="Results/AdultFFecD3YearFixed220720.csv")

# check significance of fixed effects

# sig of age - x=8.5317  df=2    p=0.01404 

adultffec.d3yr.a <- clmm(EweFecundityO ~ CoatBinF + cWeight + YearF + cD3 + cD3*YearF + (1|ID),
                       Hess=TRUE, nAGQ=10, data=adultffec)

anova(adultffec.d3yr,adultffec.d3yr.a)

# sig of coat - x=0.2149  df=1     p=0.6429

adultffec.d3yr.c <- clmm(EweFecundityO ~ AgeGroupF + cWeight + YearF + cD3 + cD3*YearF + (1|ID),
                       Hess=TRUE, nAGQ=10, data=adultffec)

anova(adultffec.d3yr,adultffec.d3yr.c)

# sig of weight - x=3.8512  df=1    p=0.04971

adultffec.d3yr.w <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + YearF + cD3 + cD3*YearF + (1|ID),
                       Hess=TRUE, nAGQ=10, data=adultffec)

anova(adultffec.d3yr,adultffec.d3yr.w)

# sig of d3*year - x=19.389  df=5   p=0.001626

adultffec.d3yr.y <- clmm(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + YearF + cD3 + (1|ID),
                       Hess=TRUE, nAGQ=10, data=adultffec)

anova(adultffec.d3yr,adultffec.d3yr.y)


## check significance of total vit d and d3 and ewe fecundity in each year

## 2011

adultffec.11 <-polr(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight,
                    data = adultffec, subset=YearF=="2011",
                    Hess = TRUE)

# sig of tot vit d - x=0.3140007 df=1 p=0.5752355

adultffec.11.dtot <-polr(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + cTotalD,
                         data = adultffec, subset=YearF=="2011",
                         Hess = TRUE)

summary(adultffec.11.dtot)

anova(adultffec.11, adultffec.11.dtot)

# sig of d3 - x=0.2602042 df=1 p=0.6099799

adultffec.11.d3 <-polr(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + cD3,
                       data = adultffec, subset=YearF=="2011",
                       Hess = TRUE)

summary(adultffec.11.d3)

anova(adultffec.11, adultffec.11.d3)

## 2012

adultffec.12 <-polr(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight,
                    data = adultffec, subset=YearF=="2012",
                    Hess = TRUE)

# sig of tot vit d - x=18.25776 df=1 p=1.929387e-05

adultffec.12.dtot <-polr(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + cTotalD,
                         data = adultffec, subset=YearF=="2012",
                         Hess = TRUE)

summary(adultffec.12.dtot)

anova(adultffec.12, adultffec.12.dtot)

# sig of d3 - x= 18.15964 df=1 p=2.031383e-05

adultffec.12.d3 <-polr(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + cD3,
                       data = adultffec, subset=YearF=="2012",
                       Hess = TRUE)

summary(adultffec.12.d3)

anova(adultffec.12, adultffec.12.d3)

## 2013

adultffec.13 <-polr(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight,
                    data = adultffec, subset=YearF=="2013",
                    Hess = TRUE)

# sig of tot vit d - x=0.4940413 df=1 p=0.4821301

adultffec.13.dtot <-polr(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + cTotalD,
                         data = adultffec, subset=YearF=="2013",
                         Hess = TRUE)

summary(adultffec.13.dtot)

anova(adultffec.13, adultffec.13.dtot)

# sig of d3 - x=1.182363 df=1 p=0.2768755

adultffec.13.d3 <-polr(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + cD3,
                       data = adultffec, subset=YearF=="2013",
                       Hess = TRUE)

summary(adultffec.13.d3)

anova(adultffec.13, adultffec.13.d3)

## 2014

adultffec.14 <-polr(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight,
                    data = adultffec, subset=YearF=="2014",
                    Hess = TRUE)

# sig of tot vit d - x=0.4162822 df=1 p=0.518798

adultffec.14.dtot <-polr(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + cTotalD,
                         data = adultffec, subset=YearF=="2014",
                         Hess = TRUE)

summary(adultffec.14.dtot)

anova(adultffec.14, adultffec.14.dtot)

# sig of d3 - x=0.4300172 df=1 p=0.5119805

adultffec.14.d3 <-polr(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + cD3,
                       data = adultffec, subset=YearF=="2014",
                       Hess = TRUE)

summary(adultffec.14.d3)

anova(adultffec.14, adultffec.14.d3)

## 2015

adultffec.15 <-polr(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight,
                    data = adultffec, subset=YearF=="2015",
                    Hess = TRUE)

# sig of tot vit d - x=2.175134 df=1 p=0.1402573

adultffec.15.dtot <-polr(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + cTotalD,
                         data = adultffec, subset=YearF=="2015",
                         Hess = TRUE)

summary(adultffec.15.dtot)

anova(adultffec.15, adultffec.15.dtot)

# sig of d3 - x=0.7561732 df=1 p=0.3845288

adultffec.15.d3 <-polr(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + cD3,
                       data = adultffec, subset=YearF=="2015",
                       Hess = TRUE)

summary(adultffec.15.d3)

anova(adultffec.15, adultffec.15.d3)

## 2016

adultffec.16 <-polr(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight,
                    data = adultffec, subset=YearF=="2016",
                    Hess = TRUE)

# sig of tot vit d - x=3.601001 df=1 p=0.05774479

adultffec.16.dtot <-polr(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + cTotalD,
                         data = adultffec, subset=YearF=="2016",
                         Hess = TRUE)

summary(adultffec.16.dtot)

anova(adultffec.16, adultffec.16.dtot)

# sig of d3 - x=3.916892 df=1 p=0.04780319

adultffec.16.d3 <-polr(EweFecundityO ~ AgeGroupF + CoatBinF + cWeight + cD3,
                       data = adultffec, subset=YearF=="2016",
                       Hess = TRUE)

summary(adultffec.16.d3)

anova(adultffec.16, adultffec.16.d3)

####### plot ewe fecundity and total vit d and d3 by year - figure for ms in graph script

adultffec$typed <- factor(ifelse(adultffec$YearF=="2012", "Highlighted", "Normal"))
adultffec$typed3 <- factor(ifelse(adultffec$YearF=="2012" | adultffec$YearF=="2016", "Highlighted", "Normal"))

adultffecd<-ggplot(adultffec,aes(x=factor(EweFecundity), y=Total.25D, fill=typed)) + geom_boxplot(outlier.colour=NA,notch=FALSE)  + 
  labs(y=expression("25(OH)D"~"(nmol/l)"), x=expression("Ewe fecundity")) + 
  theme_classic() + geom_jitter(size=1) + facet_grid(. ~ YearF) +
  theme(axis.title.x=element_text(face="bold", colour="grey16", size=12)) + 
  theme(axis.title.y=element_text(face="bold", colour="grey16", size=12)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=14, colour="grey26")) + theme(legend.position = "none") +
  scale_fill_manual(values=c("aquamarine4","#999999"))

adultffecd3<-ggplot(adultffec,aes(x=factor(EweFecundity), y=X25OHD3, fill=typed3)) + geom_boxplot(outlier.colour=NA,notch=FALSE)  +
  labs(y=expression("25(OH)D"["3"]~"(nmol/l)"), x=expression("Ewe fecundity")) +
  theme_classic() + geom_jitter(size=1) + facet_grid(. ~ YearF) + 
  theme(axis.title.x=element_text(face="bold", colour="grey16", size=12)) + 
  theme(axis.title.y=element_text(face="bold", colour="grey16", size=12)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=14, colour="grey26")) + theme(legend.position = "none") +
  scale_fill_manual(values=c("aquamarine4","#999999"))

plot_grid(adultffecd, adultffecd3, ncol = 2, nrow=1)


# # -----------------------------------------------------------------------


# Male annual breeding success - lambs


# # -----------------------------------------------------------------------

# make model subset - n=260

lambmbs  <- droplevels(subset(VITD, Age==0 & Sex==2 & !is.na(Weight) & !is.na(MaleABSBin))) 

table(lambmbs$MaleABS) # think this should be binned

table(lambmbs$MaleABSBin)

table(lambmbs$MaleABSBin, lambmbs$Survival)

table(lambmbs$MaleABSBin, lambmbs$YearF) # drop 2011??

# plot vitamin d measures by male lamb breeding success

ggplot(lambmbs,aes(x=factor(MaleABSBin), y=X25OHD2)) + geom_boxplot(outlier.colour=NA,notch=FALSE)  +
  xlab("Male lamb ABS") + ylab("25(OH)D2") +
  theme_classic() + geom_jitter() +
  theme(axis.title.x=element_text(face="bold", colour="grey16", size=12)) + 
  theme(axis.title.y=element_text(face="bold", colour="grey16", size=12)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=14, colour="grey26")) 

ggplot(lambmbs,aes(x=factor(MaleABSBin), y=X25OHD3)) + geom_boxplot(outlier.colour=NA,notch=FALSE)  +
  xlab("Male lamb ABS") + ylab("25(OH)D3") +
  theme_classic() + geom_jitter() +
  theme(axis.title.x=element_text(face="bold", colour="grey16", size=12)) + 
  theme(axis.title.y=element_text(face="bold", colour="grey16", size=12)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=14, colour="grey26")) 

ggplot(lambmbs,aes(x=factor(MaleABSBin), y=Total.25D)) + geom_boxplot(outlier.colour=NA,notch=FALSE)  +
  xlab("Male lamb ABS") + ylab("Total vitamin D") +
  theme_classic() + geom_jitter() +
  theme(axis.title.x=element_text(face="bold", colour="grey16", size=12)) + 
  theme(axis.title.y=element_text(face="bold", colour="grey16", size=12)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=14, colour="grey26")) 


# rescale variables 

lambmbs$cWeight <- (lambmbs$Weight - mean(lambmbs$Weight, na.rm=TRUE)) / sd(lambmbs$Weight, na.rm=TRUE)

lambmbs$D2sq <- lambmbs$X25OHD2 * lambmbs$X25OHD2
lambmbs$cD2 <- (lambmbs$X25OHD2 - mean(lambmbs$X25OHD2, na.rm=TRUE)) / sd(lambmbs$X25OHD2, na.rm=TRUE)
lambmbs$cD2sq <- (lambmbs$D2sq - mean(lambmbs$D2sq, na.rm=TRUE)) / sd(lambmbs$D2sq, na.rm=TRUE)

lambmbs$D3sq <- lambmbs$X25OHD3 * lambmbs$X25OHD3
lambmbs$cD3 <- (lambmbs$X25OHD3 - mean(lambmbs$X25OHD3, na.rm=TRUE)) / sd(lambmbs$X25OHD3, na.rm=TRUE)
lambmbs$cD3sq <- (lambmbs$D3sq - mean(lambmbs$D3sq, na.rm=TRUE)) / sd(lambmbs$D3sq, na.rm=TRUE)

lambmbs$TotalDsq <- lambmbs$Total.25D * lambmbs$Total.25D
lambmbs$cTotalD <- (lambmbs$Total.25D - mean(lambmbs$Total.25D, na.rm=TRUE)) / sd(lambmbs$Total.25D, na.rm=TRUE)
lambmbs$cTotalDsq <- (lambmbs$TotalDsq - mean(lambmbs$TotalDsq, na.rm=TRUE)) / sd(lambmbs$TotalDsq, na.rm=TRUE)

# compare base model with models including d2, d3, tot vitd, d2+d3 + year interactions

mlbs.base<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF, 
               data = lambmbs, family = binomial)
mlbs.d2<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cD2,
             data = lambmbs, family = binomial)
mlbs.d3<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cD3, 
             data = lambmbs, family = binomial)
mlbs.dtot<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cTotalD, 
               data = lambmbs, family = binomial)
mlbs.d2d3<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cD2 + cD3, 
               data = lambmbs, family = binomial)
mlbs.d2sq<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cD2 + cD2sq,
             data = lambmbs, family = binomial)
mlbs.d3sq<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cD3 + cD3sq, 
             data = lambmbs, family = binomial)
mlbs.dtotsq<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cTotalD + cTotalDsq, 
               data = lambmbs, family = binomial)
mlbs.d2d3sq<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cD2 + cD3 + cD2sq + cD3sq, 
               data = lambmbs, family = binomial)
mlbs.d2yr<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cD2 + cD2*YearF, 
               data = lambmbs, family = binomial)
mlbs.d3yr<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cD3 + cD3*YearF, 
               data = lambmbs, family = binomial)
mlbs.dtotyr<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cTotalD + cTotalD*YearF, 
                 data = lambmbs, family = binomial)
mlbs.d2d3yr<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cD2 + cD3 + cD2*YearF + cD3*YearF, 
                 data = lambmbs, family = binomial)

models_malelabs <- list(mlbs.base,mlbs.d2,mlbs.d3,mlbs.dtot,mlbs.d2d3,
                        mlbs.d2sq, mlbs.d3sq, mlbs.dtotsq, mlbs.d2d3sq,
                        mlbs.d2yr,mlbs.d3yr,mlbs.dtotyr,mlbs.d2d3yr)
model_malelabs_names <- c("base", "d2", "d3", "dtot", "d2d3",
                          "d2sq", "d3sq", "dtotsq", "d2d3sq",
                          "d2yr", "d3yr", "dtotyr", "d2d3yr")

# compare AICs

#aictab(models_malelabs, model_malelabs_names) 

# base lowest AIC and has lowest parameters within 2 AIC

aictab(models_malelabs, model_malelabs_names, second.ord=FALSE) 

# base model is the best model

summary(mlbs.base)
plot(mlbs.base)

# Run the models without 2011 (no males sired lambs in first year) 
# Models without 2011 are in the manuscript (all years in supp - but consistent)

# make model subset - n=213

lambmbs2  <- droplevels(subset(VITD, Age==0 & Sex==2 & !is.na(Weight) & !is.na(MaleABSBin) & Year!=2011)) # n=260

table(lambmbs2$MaleABSBin, lambmbs2$YearF)

# rescale variables

lambmbs2$cWeight <- (lambmbs2$Weight - mean(lambmbs2$Weight, na.rm=TRUE)) / sd(lambmbs2$Weight, na.rm=TRUE)

lambmbs2$D2sq <- lambmbs2$X25OHD2 * lambmbs2$X25OHD2
lambmbs2$cD2 <- (lambmbs2$X25OHD2 - mean(lambmbs2$X25OHD2, na.rm=TRUE)) / sd(lambmbs2$X25OHD2, na.rm=TRUE)
lambmbs2$cD2sq <- (lambmbs2$D2sq - mean(lambmbs2$D2sq, na.rm=TRUE)) / sd(lambmbs2$D2sq, na.rm=TRUE)

lambmbs2$D3sq <- lambmbs2$X25OHD3 * lambmbs2$X25OHD3
lambmbs2$cD3 <- (lambmbs2$X25OHD3 - mean(lambmbs2$X25OHD3, na.rm=TRUE)) / sd(lambmbs2$X25OHD3, na.rm=TRUE)
lambmbs2$cD3sq <- (lambmbs2$D3sq - mean(lambmbs2$D3sq, na.rm=TRUE)) / sd(lambmbs2$D3sq, na.rm=TRUE)

lambmbs2$TotalDsq <- lambmbs2$Total.25D * lambmbs2$Total.25D
lambmbs2$cTotalD <- (lambmbs2$Total.25D - mean(lambmbs2$Total.25D, na.rm=TRUE)) / sd(lambmbs2$Total.25D, na.rm=TRUE)
lambmbs2$cTotalDsq <- (lambmbs2$TotalDsq - mean(lambmbs2$TotalDsq, na.rm=TRUE)) / sd(lambmbs2$TotalDsq, na.rm=TRUE)

# compare base model with models including d2, d3, tot vitd, d2+d3 + year interactions

mlbs.base2<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF, 
               data = lambmbs2, family = binomial)
mlbs.d22<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cD2,
             data = lambmbs2, family = binomial)
mlbs.d32<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cD3, 
             data = lambmbs2, family = binomial)
mlbs.dtot2<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cTotalD, 
               data = lambmbs2, family = binomial)
mlbs.d2d32<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cD2 + cD3, 
               data = lambmbs2, family = binomial)
mlbs.d2sq2<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cD2 + cD2sq,
               data = lambmbs2, family = binomial)
mlbs.d3sq2<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cD3 + cD3sq, 
               data = lambmbs2, family = binomial)
mlbs.dtotsq2<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cTotalD + cTotalDsq, 
                 data = lambmbs2, family = binomial)
mlbs.d2d3sq2<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cD2 + cD3 + cD2sq + cD3sq, 
                 data = lambmbs2, family = binomial)
mlbs.d2yr2<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cD2 + cD2*YearF, 
               data = lambmbs2, family = binomial)
mlbs.d3yr2<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cD3 + cD3*YearF, 
               data = lambmbs2, family = binomial)
mlbs.dtotyr2<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cTotalD + cTotalD*YearF, 
                 data = lambmbs2, family = binomial)
mlbs.d2d3yr2<-glm(MaleABSBin ~ CoatBinF + cWeight + YearF + cD2 + cD3 + cD2*YearF + cD3*YearF, 
                 data = lambmbs2, family = binomial)

models_malelabs2 <- list(mlbs.base2,mlbs.d22,mlbs.d32,mlbs.dtot2,mlbs.d2d32,
                        mlbs.d2sq2, mlbs.d3sq2, mlbs.dtotsq2, mlbs.d2d3sq2,
                        mlbs.d2yr2,mlbs.d3yr2,mlbs.dtotyr2,mlbs.d2d3yr2)
model_malelabs_names2 <- c("base", "d2", "d3", "dtot", "d2d3",
                          "d2sq", "d3sq", "dtotsq", "d2d3sq",
                          "d2yr", "d3yr", "dtotyr", "d2d3yr")

# compare AICs

#aictab(models_malelabs2, model_malelabs_names2) 
aictab(models_malelabs2, model_malelabs_names2, second.ord=FALSE) 
 
# base model is still the best fitting model

summary(mlbs.base2)
plot(mlbs.base2)

# # -----------------------------------------------------------------------

# Adult male breeding successs

# # -----------------------------------------------------------------------

# make model subset - n=181 measures of 109 IDs

adultmbs  <- droplevels(subset(VITD, Age>0 & Sex==2 & !is.na(Weight) & !is.na(MaleABS))) 

# check

table(adultmbs$MaleABS) 
hist(adultmbs$MaleABS)
table(adultmbs$MaleABSBin, adultmbs$YearF) 

# rescale variables

adultmbs$cWeight <- (adultmbs$Weight - mean(adultmbs$Weight, na.rm=TRUE)) / sd(adultmbs$Weight, na.rm=TRUE)

adultmbs$D2sq <- adultmbs$X25OHD2 * adultmbs$X25OHD2
adultmbs$cD2 <- (adultmbs$X25OHD2 - mean(adultmbs$X25OHD2, na.rm=TRUE)) / sd(adultmbs$X25OHD2, na.rm=TRUE)
adultmbs$cD2sq <- (adultmbs$D2sq - mean(adultmbs$D2sq, na.rm=TRUE)) / sd(adultmbs$D2sq, na.rm=TRUE)

adultmbs$D3sq <- adultmbs$X25OHD3 * adultmbs$X25OHD3
adultmbs$cD3 <- (adultmbs$X25OHD3 - mean(adultmbs$X25OHD3, na.rm=TRUE)) / sd(adultmbs$X25OHD3, na.rm=TRUE)
adultmbs$cD3sq <- (adultmbs$D3sq - mean(adultmbs$D3sq, na.rm=TRUE)) / sd(adultmbs$D3sq, na.rm=TRUE)

adultmbs$TotalDsq <- adultmbs$Total.25D * adultmbs$Total.25D
adultmbs$cTotalD <- (adultmbs$Total.25D - mean(adultmbs$Total.25D, na.rm=TRUE)) / sd(adultmbs$Total.25D, na.rm=TRUE)
adultmbs$cTotalDsq <- (adultmbs$TotalDsq - mean(adultmbs$TotalDsq, na.rm=TRUE)) / sd(adultmbs$TotalDsq, na.rm=TRUE)


# looking at the best fit model from glmmTMB package - poisson, nbinom1, nbinom2 w and wo zero-inflation


mabs.zipoisson <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+(1|ID),
                          data=adultmbs,
                          ziformula=~1,
                          family=poisson)
summary(mabs.zipoisson)

mabs.zinbinom <- update(mabs.zipoisson, family=nbinom2)

mabs.zinbinom2 <- update(mabs.zipoisson, family=nbinom1)

mabs.poisson <- update(mabs.zipoisson, ziformula=~0, family=poisson)

mabs.nbinom2 <- update(mabs.zipoisson, ziformula=~0, family=nbinom2)

mabs.nbinom1 <- update(mabs.zipoisson, ziformula=~0, family=nbinom1)


AIC(mabs.zipoisson, mabs.zinbinom, mabs.zinbinom2, mabs.poisson, mabs.nbinom2, mabs.nbinom1)

# model with nbinom1 and no zero inflation has lowest AIC - use this
# check results with poisson as this is very similar (AIC diff of 0.7358)

# check these models

summary(mabs.nbinom1)
summary(mabs.poisson)

fixef(mabs.nbinom1) # parameter values > 10 are suspect - all okay
fixef(mabs.poisson) # parameter values > 10 are suspect - all okay

### nbinom family has lowest AIC

mabs.base <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+(1|ID),
                     data=adultmbs,
                     ziformula=~0,
                     family=nbinom1)

# model diagnostics - looks good

simulationOutput <- simulateResiduals(fittedModel = mabs.base, plot = T)


# # compare base model with models including d2, d3, tot vitd, d2+d3 + age/year interactions

mabs.base <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+(1|ID),
                     data=adultmbs,
                     ziformula=~0,
                     family=nbinom1)

mabs.d2<- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD2+(1|ID),
                  data=adultmbs,
                  ziformula=~0,
                  family=nbinom1)

mabs.d3 <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD3+(1|ID),
                   data=adultmbs,
                   ziformula=~0,
                   family=nbinom1)

mabs.d2d3 <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD2+cD3+(1|ID),
                     data=adultmbs,
                     ziformula=~0,
                     family=nbinom1)

mabs.dtot <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cTotalD+(1|ID),
                     data=adultmbs,
                     ziformula=~0,
                     family=nbinom1)

mabs.d2sq<- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD2+cD2sq+(1|ID),
                  data=adultmbs,
                  ziformula=~0,
                  family=nbinom1)

mabs.d3sq <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD3+cD3sq+(1|ID),
                   data=adultmbs,
                   ziformula=~0,
                   family=nbinom1)

mabs.d2d3sq <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD2+cD3+cD2sq+cD3sq+(1|ID),
                     data=adultmbs,
                     ziformula=~0,
                     family=nbinom1)

mabs.dtotsq <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cTotalD+cTotalDsq+(1|ID),
                     data=adultmbs,
                     ziformula=~0,
                     family=nbinom1)

mabs.d2yr <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD2+cD2*YearF+(1|ID),
                     data=adultmbs,
                     ziformula=~0,
                     family=nbinom1)

mabs.d3yr <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD3+cD3*YearF+(1|ID),
                     data=adultmbs,
                     ziformula=~0,
                     family=nbinom1)

mabs.d2d3yr <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD2+cD3+cD2*YearF+cD3*YearF+(1|ID),
                       data=adultmbs,
                       ziformula=~0,
                       family=nbinom1)

mabs.dtotyr <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cTotalD+cTotalD*YearF+(1|ID),
                       data=adultmbs,
                       ziformula=~0,
                       family=nbinom1)

mabs.d2age <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD2+cD2*AgeGroupF+(1|ID),
                     data=adultmbs,
                     ziformula=~0,
                     family=nbinom1)

mabs.d3age <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD3+cD3*AgeGroupF+(1|ID),
                     data=adultmbs,
                     ziformula=~0,
                     family=nbinom1)

mabs.d2d3age <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD2+cD3+cD2*AgeGroupF+cD3*AgeGroupF+(1|ID),
                       data=adultmbs,
                       ziformula=~0,
                       family=nbinom1)

mabs.dtotage <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cTotalD+cTotalD*AgeGroupF+(1|ID),
                       data=adultmbs,
                       ziformula=~0,
                       family=nbinom1)

# compare AICs of models

AIC(mabs.base, mabs.d2, mabs.d3, mabs.d2d3, mabs.dtot, mabs.d2yr, mabs.d3yr, mabs.d2d3yr, mabs.dtotyr,
    mabs.d2age, mabs.d3age, mabs.dtotage, mabs.d2d3age, mabs.d2sq, mabs.d3sq, mabs.dtotsq, mabs.d2d3sq)

mabsAIC <- data.frame(AIC(mabs.base, mabs.d2, mabs.d3, mabs.d2d3, mabs.dtot, mabs.d2yr, mabs.d3yr, mabs.d2d3yr, mabs.dtotyr,
                          mabs.d2age, mabs.d3age, mabs.dtotage, mabs.d2d3age, mabs.d2sq, mabs.d3sq, mabs.dtotsq, mabs.d2d3sq))

mabsAIC$minAIC <- min(mabsAIC$AIC)
mabsAIC$dAIC <- mabsAIC$AIC-mabsAIC$minAIC 

write.csv(mabsAIC, file="Results/MaleAdultBSVitDWOCoat210720.csv")

# base model is the best model

# check results are consistent with the poisson family (as AICs very similar)

mabs.base.2 <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+(1|ID),
                     data=adultmbs,
                     ziformula=~0,
                     family=poisson)

mabs.d2.2<- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD2+(1|ID),
                  data=adultmbs,
                  ziformula=~0,
                  family=poisson)

mabs.d3.2 <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD3+(1|ID),
                   data=adultmbs,
                   ziformula=~0,
                   family=poisson)

mabs.d2d3.2 <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD2+cD3+(1|ID),
                     data=adultmbs,
                     ziformula=~0,
                     family=poisson)

mabs.dtot.2 <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cTotalD+(1|ID),
                     data=adultmbs,
                     ziformula=~0,
                     family=poisson)

mabs.d2sq.2<- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD2+cD2sq+(1|ID),
                    data=adultmbs,
                    ziformula=~0,
                    family=poisson)

mabs.d3sq.2 <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD3+cD3sq+(1|ID),
                     data=adultmbs,
                     ziformula=~0,
                     family=poisson)

mabs.d2d3sq.2 <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD2+cD3+cD2sq+cD3sq+(1|ID),
                       data=adultmbs,
                       ziformula=~0,
                       family=poisson)

mabs.dtotsq.2 <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cTotalD+cTotalDsq+(1|ID),
                       data=adultmbs,
                       ziformula=~0,
                       family=poisson)

mabs.d2yr.2 <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD2+cD2*YearF+(1|ID),
                     data=adultmbs,
                     ziformula=~0,
                     family=poisson)

mabs.d3yr.2 <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD3+cD3*YearF+(1|ID),
                     data=adultmbs,
                     ziformula=~0,
                     family=poisson)

mabs.d2d3yr.2 <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD2+cD3+cD2*YearF+cD3*YearF+(1|ID),
                       data=adultmbs,
                       ziformula=~0,
                       family=poisson)

mabs.dtotyr.2 <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cTotalD+cTotalD*YearF+(1|ID),
                       data=adultmbs,
                       ziformula=~0,
                       family=poisson)

mabs.d2age.2 <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD2+cD2*AgeGroupF+(1|ID),
                      data=adultmbs,
                      ziformula=~0,
                      family=poisson)

mabs.d3age.2 <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD3+cD3*AgeGroupF+(1|ID),
                      data=adultmbs,
                      ziformula=~0,
                      family=poisson)

mabs.d2d3age.2 <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cD2+cD3+cD2*AgeGroupF+cD3*AgeGroupF+(1|ID),
                        data=adultmbs,
                        ziformula=~0,
                        family=poisson)

mabs.dtotage.2 <- glmmTMB(MaleABS~AgeGroupF+CoatBinF+cWeight+YearF+cTotalD+cTotalD*AgeGroupF+(1|ID),
                        data=adultmbs,
                        ziformula=~0,
                        family=poisson)

AIC(mabs.base.2, mabs.d2.2, mabs.d3.2, mabs.d2d3.2, mabs.dtot.2, mabs.d2yr.2, mabs.d3yr.2, mabs.d2d3yr.2, mabs.dtotyr.2,
    mabs.d2age.2, mabs.d3age.2, mabs.dtotage.2, mabs.d2d3age.2, mabs.d2sq.2, mabs.d3sq.2, mabs.dtotsq.2, mabs.d2d3sq.2)

mabsAIC2 <- data.frame(AIC(mabs.base.2, mabs.d2.2, mabs.d3.2, mabs.d2d3.2, mabs.dtot.2, mabs.d2yr.2, mabs.d3yr.2, mabs.d2d3yr.2, mabs.dtotyr.2,
                          mabs.d2age.2, mabs.d3age.2, mabs.dtotage.2, mabs.d2d3age.2, mabs.d2sq.2, mabs.d3sq.2, mabs.dtotsq.2, mabs.d2d3sq.2))

summary(mabsAIC2)

mabsAIC2$minAIC <- min(mabsAIC2$AIC)
mabsAIC2$dAIC <- mabsAIC2$AIC-mabsAIC2$minAIC 

# model with lowest aic within 2 aic of base model
# model results are consistent with either poisson and nbinom1 - base model is the best model

# model diagnostics are okay with both nbinom1 and poisson

simulationOutput <- simulateResiduals(fittedModel = mabs.base.2, plot = T)

