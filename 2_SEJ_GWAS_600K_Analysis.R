#  ------------------------------------------------------------------------

#  Vitamin D ms
#  GWAS with 600K SNP Dataset
#  AMS - 14/6/19

#  ------------------------------------------------------------------------

# remove all data

rm(list=ls())

library(plyr)
library(dplyr)
library(lubridate)
library(ggplot2) 
library(lme4)
library(ggplot2)
library(GenABEL)
library(RepeatABEL)

source("r/rGLSadj.R")

load("Vit_Reg_h2.RData")

grminv2 <- grminv
grminv2[lower.tri(grminv2)] = t(grminv2)[lower.tri(grminv2)]


#~~ Load GenABEL genotype data

# load("../../../Soay Sheep Genomic Data/20200517_imputation/sheep_geno_imputed_oar31_17052020.GenABEL.RData")
# load("../../../Soay Sheep Genomic Data/20190711 Soay Sheep 50K Data/Previous_Versions/Plates_1to87_QC3.GenABEL.RData")
# 
# snp50 <- snpnames(soay87)
# snpHD <- snpnames(soayimp)
# snp50X <- snpnames(soay87[,chromosome(soay87) == 27])
# 
# #~~ Add SNPs that weren't imputed (includes the entire X chromosome)
# 
# soay87 <- soay87[,which(!snpnames(soay87) %in% snpnames(soayimp))]
# soayimp <- merge.gwaa.data(soayimp, soay87)
# rm(soay87)
# save(soayimp, file = "soayimp_genotype_data.RData")

load("soayimp_genotype_data.RData")


#~~ Get minor allele frequency data

mafinfo <- summary.snp.data(gtdata(soayimp))
mafinfo$SNP.Name <- row.names(mafinfo)
mafinfo <- subset(mafinfo, select = c(Q.2, SNP.Name))
names(mafinfo)[1] <- "MAF"

#~~ Get Map info (relative to Rambouillet)

mapdata <- data.frame(Chr = chromosome(soayimp),
                      Position = map(soayimp),
                      SNP.Name = snpnames(soayimp))
mapdata$Chr <- as.numeric(as.character(mapdata$Chr))

mapdata <- arrange(mapdata, Chr, Position)
mapdata$Chr <- as.character(mapdata$Chr)
mapdata$SNP.Name <- as.character(mapdata$SNP.Name)
mapdata$Diff <- c(0, diff(mapdata$Position))
mapdata$Diff[which(mapdata$Diff < 0)] <- 1000
mapdata$Cumu <- cumsum(mapdata$Diff)

#~~ Get Chromosome Info

chrinfo <- NULL

for(i in na.omit(unique(mapdata$Chr))){
  
  temp1 <- arrange(subset(mapdata, Chr == i), Cumu)
  
  temp2 <- data.frame(Chr = i,
                      Start = temp1[1,"Cumu"],
                      Stop = temp1[nrow(temp1),"Cumu"])
  
  chrinfo <- rbind(chrinfo, temp2)
  rm(temp1, temp2)
}

chrinfo$Mid <- chrinfo$Start + ((chrinfo$Stop - chrinfo$Start)/2)


bonf = 3.858528e-07

mapdata <- subset(mapdata, select = c(SNP.Name, Cumu))

firstRun <- FALSE

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Run the GWASs                                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

VITDSUBSET <- subset(VITD, ID %in% dimnames(grminv2)[1][[1]])
VITDSUBSET <- na.omit(subset(VITDSUBSET, select = c(ID, Total.25D, Sex, AgeF, YearF, Age, MumID, BirthYearF, Age, Age2)))

VITD2SUBSET <- subset(VITD, ID %in% dimnames(grminv2)[1][[1]])
VITD2SUBSET <- na.omit(subset(VITD2SUBSET, select = c(ID, X25OHD2, Sex, AgeF, YearF, Age, MumID, BirthYearF, Age, Age2)))

VITD3SUBSET <- subset(VITD, ID %in% dimnames(grminv2)[1][[1]])
VITD3SUBSET <- na.omit(subset(VITD3SUBSET, select = c(ID, X25OHD3, Sex, AgeF, YearF, Age, MumID, BirthYearF, Age, Age2)))



if(firstRun){
  
  #~~ ALL SHEEP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # Total D
  
  gwasvitd_prefit <- preFitModel(fixed = Total.25D ~ 1+Sex+Age+Age2+YearF,
                                 random = ~1|MumID + 1|BirthYearF, 
                                 id.name = "ID",
                                 genabel.data = soayimp,
                                 phenotype.data = VITDSUBSET, 
                                 corStruc = list(ID = list("GRM"),
                                                 MumID = list("Ind"),
                                                 BirthYearF = list("Ind")),
                                 GRM = grminv2)
  
  gwasvitd <- rGLSadj(formula.FixedEffects=Total.25D~Sex+Age+Age2+YearF,
                      genabel.data=soayimp,
                      phenotype.data=VITDSUBSET,
                      id="ID",
                      V = gwasvitd_prefit$V,
                      GRM=grminv2)
  
  gwasvitd <- process_rGLSadj_results(gwasvitd, soayimp)
  
  # Total D2
  
  gwasvitd2_prefit <- preFitModel(fixed = X25OHD2 ~ 1+Sex+Age+Age2+YearF,
                                 random = ~1|MumID + 1|BirthYearF, 
                                 id.name = "ID",
                                 genabel.data = soayimp,
                                 phenotype.data = VITD2SUBSET, 
                                 corStruc = list(ID = list("GRM"),
                                                 MumID = list("Ind"),
                                                 BirthYearF = list("Ind")),
                                 GRM = grminv2)
  
  gwasvitd2 <- rGLSadj(formula.FixedEffects=X25OHD2~Sex+Age+Age2+YearF,
                      genabel.data=soayimp,
                      phenotype.data=VITD2SUBSET,
                      id="ID",
                      V = gwasvitd2_prefit$V,
                      GRM=grminv2)
  
  gwasvitd2 <- process_rGLSadj_results(gwasvitd2, soayimp)
  
  
  # Total D3
  
  gwasvitd3_prefit <- preFitModel(fixed = X25OHD3 ~ 1+Sex+Age+Age2+YearF,
                                  random = ~1|MumID + 1|BirthYearF, 
                                  id.name = "ID",
                                  genabel.data = soayimp,
                                  phenotype.data = VITD3SUBSET, 
                                  corStruc = list(ID = list("GRM"),
                                                  MumID = list("Ind"),
                                                  BirthYearF = list("Ind")),
                                  GRM = grminv2)
  
  gwasvitd3 <- rGLSadj(formula.FixedEffects=X25OHD3~Sex+Age+Age2+YearF,
                       genabel.data=soayimp,
                       phenotype.data=VITD3SUBSET,
                       id="ID",
                       V = gwasvitd3_prefit$V,
                       GRM=grminv2)
  
  gwasvitd3 <- process_rGLSadj_results(gwasvitd3, soayimp)
  
  
  save(gwasvitd, gwasvitd2, gwasvitd3, file = "results/gwasvitd.RData")
  rm(gwasvitd, gwasvitd2, gwasvitd3, gwasvitd_prefit, gwasvitd2_prefit, gwasvitd3_prefit)
  
  gc()
  
  
}

load("results/gwasvitd.RData")

gwas.results <- rbind(cbind(gwasvitd, Age = "All", Model = "Total"),
                      cbind(gwasvitd2, Age = "All", Model = "D2"),
                      cbind(gwasvitd3, Age = "All", Model = "D3"))

rm(gwasvitd, gwasvitd2, gwasvitd3)

#~~ Get MAF information

idvec <- unique(VITD$ID) %>% as.character
idvec <- idvec[which(idvec %in% idnames(soayimp))]

snpinfo <- summary.snp.data(gtdata(soayimp[idvec,]))
head(snpinfo)

snpinfo$SNP.Name <- row.names(snpinfo)
snpinfo <- subset(snpinfo, select = c(SNP.Name, Q.2))
head(snpinfo)

gwas.results <- join(gwas.results, snpinfo)
table(gwas.results$Chromosome)
gwas.results$Chromosome <- as.numeric(as.character(gwas.results$Chromosome))

gwas.results <- join(gwas.results, mapdata)

gwas.results <- gwas.results[,c("Model", "SNP.Name", "Chromosome", "Position", "A1", "A2", "effB", "se_effB", 
                               "chi2.1df", "P1df", "Pc1df",  "ExpP",  
                               "Q.2", "Cumu")]

gwas.results$Model[which(gwas.results$Model == "Total")] <- "25(OH)D"
gwas.results$Model[which(gwas.results$Model == "D2")] <- "25(OH)D2"
gwas.results$Model[which(gwas.results$Model == "D3")] <- "25(OH)D3"

write.table(gwas.results, file = "results/2_SEJ_GWAS_Results_600K.txt", row.names = F, sep = "\t", quote = F)

rm(soayimp)
gc()



beepr::beep()
gc()
