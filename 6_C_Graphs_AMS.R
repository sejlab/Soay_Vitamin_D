#  ------------------------------------------------------------------------

#  Vitamin D ms
#  Graph code
#  AMS - 26/7/20

#  ------------------------------------------------------------------------

# remove all data

rm(list=ls())

# read in data - 1452 obs

VITD <- read.csv("VitD_MS_Data_240720.csv")

# libraries 

library(ggplot2)
library(cowplot)
library(scales)

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

# Figure 3 - adult survival & d2

# # -----------------------------------------------------------------------

# model subset

ADULTSURV<-droplevels(subset(VITD, Age>0 & !is.na(Survival) & !is.na(Weight)))

# rescale variables

ADULTSURV$cWeight <- (ADULTSURV$Weight - mean(ADULTSURV$Weight, na.rm=TRUE)) / sd(ADULTSURV$Weight, na.rm=TRUE)
ADULTSURV$cD2 <- (ADULTSURV$X25OHD2 - mean(ADULTSURV$X25OHD2, na.rm=TRUE)) / sd(ADULTSURV$X25OHD2, na.rm=TRUE)
ADULTSURV$cD3 <- (ADULTSURV$X25OHD3 - mean(ADULTSURV$X25OHD3, na.rm=TRUE)) / sd(ADULTSURV$X25OHD3, na.rm=TRUE)
ADULTSURV$cTotalD <- (ADULTSURV$Total.25D - mean(ADULTSURV$Total.25D, na.rm=TRUE)) / sd(ADULTSURV$Total.25D, na.rm=TRUE)

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

plot<-wsurvsgraph

jpeg("Graphs/Figure3_AdultSurvivalD2_220720.jpeg", width=2000, height=1800, res=400)

plot

dev.off()


# # -----------------------------------------------------------------------

# Figure 4 - ewe fecundity and total vit d/d3 by year

# # -----------------------------------------------------------------------

# model subset

adultffec  <- droplevels(subset(VITD, Age>0 & !is.na(EweFecundityO) & !is.na(Weight))) # 578 obs

# graph 

adultffec$typed <- factor(ifelse(adultffec$YearF=="2012", "Highlighted", "Normal"))
adultffec$typed3 <- factor(ifelse(adultffec$YearF=="2012" | adultffec$YearF=="2016", "Highlighted", "Normal"))

adultffecd<-ggplot(adultffec,aes(x=factor(EweFecundity), y=Total.25D, fill=typed)) + geom_boxplot(outlier.colour=NA,notch=FALSE)  + 
  labs(y=expression("25(OH)D"~"(nmol/l)"), x=expression("Ewe fecundity")) + 
  theme_classic() + geom_jitter(size=1, alpha=0.3, width=0.25) + facet_grid(. ~ YearF) +
  theme(axis.title.x=element_text(face="bold", colour="grey16", size=12)) + 
  theme(axis.title.y=element_text(face="bold", colour="grey16", size=12)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=14, colour="grey26")) + theme(legend.position = "none") +
  scale_fill_manual(values=c("aquamarine4","#999999"))

adultffecd3<-ggplot(adultffec,aes(x=factor(EweFecundity), y=X25OHD3, fill=typed3)) + geom_boxplot(outlier.colour=NA,notch=FALSE)  +
  labs(y=expression("25(OH)D"["3"]~"(nmol/l)"), x=expression("Ewe fecundity")) +
  theme_classic() + geom_jitter(size=1, alpha=0.3, width=0.25) + facet_grid(. ~ YearF) + 
  theme(axis.title.x=element_text(face="bold", colour="grey16", size=12)) + 
  theme(axis.title.y=element_text(face="bold", colour="grey16", size=12)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=14, colour="grey26")) + theme(legend.position = "none") +
  scale_fill_manual(values=c("aquamarine4","#999999"))

plot_grid(adultffecd, adultffecd3, ncol = 2, nrow=1)

plot<-plot_grid(adultffecd, adultffecd3, ncol = 2, nrow=1)

jpeg("Graphs/Figure4_AdultFemaleFecundity_220720.jpeg", width=4000, height=1250, res=400)

plot

dev.off()


# # -----------------------------------------------------------------------

# Figure S4 - temporal correlation of vitamin d measures

# # -----------------------------------------------------------------------

## make sure dataframe is ordered by capt year within IDs

VITDX <- VITD[with(VITD, order(ID, Year)), ]
VITD0<-VITDX

## shift year column up by one row - next year

VITD0$NextYear <- c(VITD0$Year[-1], NA)

## do the same for vitd measures - next vitd measure

VITD0$NextD3 <- c(VITD0$X25OHD3[-1], NA)
VITD0$NextD2 <- c(VITD0$X25OHD2[-1], NA)
VITD0$NextDTOT <- c(VITD0$Total.25D[-1], NA)

## NA where next observation is of a different individual 
## only interested in comparing measures of the same individual, so NA if IDs don't match

VITD0$NextYear <- c(ifelse (VITD0$ID[-nrow(VITD0)] == VITD0$ID[-1], VITD0$NextYear, NA), NA)

VITD0$NextD3 <- c(ifelse (VITD0$ID[-nrow(VITD0)] == VITD0$ID[-1], VITD0$NextD3, NA), NA)
VITD0$NextD2 <- c(ifelse (VITD0$ID[-nrow(VITD0)] == VITD0$ID[-1], VITD0$NextD2, NA), NA)
VITD0$NextDTOT <- c(ifelse (VITD0$ID[-nrow(VITD0)] == VITD0$ID[-1], VITD0$NextDTOT, NA), NA)

# graphs of next vitamin d measure, may not be next year

#qplot(VITD0$X25OHD3, VITD0$NextD3)
#qplot(VITD0$X25OHD2, VITD0$NextD2)
#qplot(VITD0$Total.25D, VITD0$NextDTOT)

## this just looks at next observation 
## but to look at consecutive years we need observations where t+1
## so calculate diff in years - gap

VITD0$Gap <- VITD0$NextYear - VITD0$Year

#table(VITD0$Gap)

## vast majority of repeats are in successive years

## within individual change in vitd - for interest

VITD0$ChD3 <- VITD0$NextD3 - VITD0$X25OHD3
VITD0$ChD2 <- VITD0$NextD2 - VITD0$X25OHD2
VITD0$ChDTOT <- VITD0$NextDTOT - VITD0$Total.25D

## plot vitd this year againt vitd next year - but only for adults 
## (as lambs have lower measures and different distribution)

## Total 

temp2 <- subset(VITD0, !is.na(NextDTOT) & !is.na(Total.25D) & Gap==1 & Age>0)
# plot (temp2$NextDTOT ~ temp2$Total.25D)
# abline (0,1)
# hist(temp2$ChDTOT)
# mean(temp2$ChDTOT)

## D2

temp3 <- subset(VITD0, !is.na(NextD2) & !is.na(X25OHD2) & Gap==1 & Age>0)
# plot (temp3$NextD2 ~ temp3$X25OHD2)
# abline (0,1)
# hist(temp3$ChD2)
# mean(temp3$ChD2)

## D3

temp4 <- subset(VITD0, !is.na(NextD3) & !is.na(X25OHD3) & Gap==1 & Age>0)
# plot (temp4$NextD3 ~ temp4$X25OHD3)
# abline (0,1)
# hist(temp3$ChD3)
# mean(temp3$ChD3)

# graph

Reptotd<-ggplot(temp2,aes(x=Total.25D,y=NextDTOT)) + geom_point(alpha=0.2) +
  xlab("25(OH)D at age t") + ylab("25(OH)D at age t+1") + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=14)) + 
  theme(axis.title.y=element_text(colour="grey6", size=14)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=12, colour="grey26")) + geom_abline(intercept=0, slope=1, linetype="dashed", size=1)

RepD2<-ggplot(temp3,aes(x=X25OHD2,y=NextD2)) + geom_point(alpha=0.2) +
  labs(x=expression("25(OH)D"["2"]~"at age t"), y=expression("25(OH)D"["2"]~"at age t+1")) + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=14)) + 
  theme(axis.title.y=element_text(colour="grey6", size=14)) +
  theme(axis.line=element_line(colour="grey60")) +
  theme(text = element_text(size=12, colour="grey26")) + geom_abline(intercept=0, slope=1, linetype="dashed", size=1)

RepD3<-ggplot(temp4,aes(x=X25OHD3,y=NextD3)) + geom_point(alpha=0.2) +
  labs(x=expression("25(OH)D"["3"]~"at age t"), y=expression("25(OH)D"["3"]~"at age t+1")) + theme_classic() +
  theme(axis.title.x=element_text(colour="grey6", size=14)) + 
  theme(axis.title.y=element_text(colour="grey6", size=14)) +
  theme(axis.line=element_line(colour="grey60")) + 
  theme(text = element_text(size=12, colour="grey26")) + 
  geom_abline(intercept=0, slope=1, linetype="dashed", size=1) +
  scale_x_continuous(breaks=pretty_breaks(n=2)) +
  scale_y_continuous(breaks=pretty_breaks(n=2))

graph1<-plot_grid(Reptotd, RepD2, RepD3, ncol = 3, nrow=1)
graph1

cor.test(x=temp2$Total.25D,y=temp2$NextDTOT)
cor.test(x=temp3$X25OHD2,y=temp3$NextD2)
cor.test(x=temp4$X25OHD3,y=temp4$NextD3)

plot<-graph1

jpeg("Graphs/FigureS4_TemporalCorrelations_260720.jpeg", width=4000, height=1250, res=400)

plot

dev.off()
