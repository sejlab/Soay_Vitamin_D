#  ------------------------------------------------------------------------

#  Vitamin D ms
#  AnMod structure for Susie for reg h2
#  AMS - 18/5/19

#  ------------------------------------------------------------------------



library(plyr)
library(dplyr)
library(ggplot2) 
library(asreml)

source("r/ASReml4.EstEffects.R")

load("Vit_Reg_h2.RData", verbose = T)


attr(grminv, which = "INVERSE") <- TRUE

# gwas.results <- read.table("results/2_SEJ_GWAS_Results_600K.txt", header = T, stringsAsFactors = F, sep = "\t")
# gwas.results <- subset(gwas.results, Pc1df < 1.28e-6)
# 
# setwd("../../../Soay Sheep Genomic Data/20190711 Soay Sheep 50K Data/Previous_Versions/")
# 
# maptab <- read.table("Plates_1to87_QC3.bim", stringsAsFactors = F)
# 
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# # 1. Get our 50K SNPs for constructing GRM               #
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# x <- subset(maptab, V1 == 18)
# 
# xless <- which(x$V4 < gwas.results$Position[1])
# xless <- xless[length(xless)]
# xless
# 
# x[c((xless - 9):(xless+10)),]
# 
# snplist <- x[c((nrow(x)-19):nrow(x)),"V2"]
# 
# writeLines(snplist, "vitd_list.txt")
# 
# system(paste0("gcta64.exe --bfile Plates_1to87_QC3 --autosome --autosome-num 26 --extract vitd_list.txt --make-grm-gz --out vitD_reg_GRM" ))
# system(paste0("gcta64.exe --grm-gz vitD_reg_GRM --grm-adj 0 --make-grm-gz --out vitD_reg_GRM.adj"))
# 
# grm.region <- read.table("vitD_reg_GRM.adj.grm.gz")  # CONTAINS REALIZED RELATEDNESS BETWEEN ALL GENOTYPED INDIVIDUALS
# ids.region <- read.table("vitD_reg_GRM.adj.grm.id")  # CONTAINS ID LIST
# 
# grmreg <- makeGRM(grm.region, ids.region, id.vector = VITD$ID) # vector of IDs from the datasset that you use for the asreml model
# dim(grmreg)
# 
# attr(grmreg, which = "INVERSE") <- TRUE
# 
# setwd("../../../Soay Sheep Projects/Vitamin_D/20191004_Vitamin_D/")
#
# save(grmreg, file = "results/4_Regional_GRM.RData")

load("results/4_Regional_GRM.RData")
vitdped$ID2 <- vitdped$ID

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Animal Models                                       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# total vitamin D 

anmodtotvitd_old <- asreml(fixed = Total.25D~1+Sex+Age+Age2+YearF, #etc etc
                       random = ~ vm(ID, grminv) + ide(ID) + MumID +BirthYearF,   #vm(ID, ainv) is the relatedness matricx, ide(ID) is the individual identity
                       data = vitdped,
                       residual = ~idv(units)) # include this line
anmodtotvitd_old <- anmodtotvitd_old$loglik


anmodtotvitd <- asreml(fixed = Total.25D~1+Sex+Age+Age2+YearF, #etc etc
                       random = ~ vm(ID, grminv) + vm(ID2, grmreg) + ide(ID) + MumID +BirthYearF,   #vm(ID, ainv) is the relatedness matricx, ide(ID) is the individual identity
                       data = vitdped,
                       residual = ~idv(units)) # include this line

wald.asreml(anmodtotvitd)
summary.asreml(anmodtotvitd, coef = T)$coef.fixed # Fixed effects
summary.asreml(anmodtotvitd, coef = T)$varcomp    # random effects
asreml4pin(anmodtotvitd)



pchisq(2*(anmodtotvitd$loglik - anmodtotvitd_old), 1, lower.tail = F)


raneff <- summary.asreml(anmodtotvitd, coef = T)$varcomp
raneff <- raneff[-nrow(raneff),]
raneff <- cbind(raneff, asreml4pin(anmodtotvitd))

# d2

anmodd2_old <- asreml(fixed = X25OHD2~1+Sex+Age+Age2+YearF, #etc etc
                           random = ~ vm(ID, grminv) + ide(ID) + MumID +BirthYearF,   #vm(ID, ainv) is the relatedness matricx, ide(ID) is the individual identity
                           data = vitdped,
                           residual = ~idv(units)) # include this line
anmodd2_old <- anmodd2_old$loglik

anmodd2 <- asreml(fixed = X25OHD2~1+Sex+Age+Age2+YearF, #etc etc
                  random = ~ vm(ID, grminv) + vm(ID2, grmreg) + ide(ID) +MumID+BirthYearF,   #vm(ID, ainv) is the relatedness matricx, ide(ID) is the individual identity
                  data = vitdped,
                  na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units)) # include this line

wald.asreml(anmodd2)
summary.asreml(anmodd2, coef = T)$coef.fixed # Fixed effects
summary.asreml(anmodd2, coef = T)$varcomp    # random effects
asreml4pin(anmodd2)

pchisq(2*(anmodd2$loglik - anmodd2_old), 1, lower.tail = F)


raneff2 <- summary.asreml(anmodd2, coef = T)$varcomp
raneff2 <- raneff2[-nrow(raneff2),]
raneff2 <- cbind(raneff2, asreml4pin(anmodd2))

vpredict(anmodd2, h2~V4/V3+V4)


# d3


anmodd3_old <- asreml(fixed = X25OHD3~1+Sex+Age+Age2+YearF, #etc etc
                  random = ~ vm(ID, grminv) + ide(ID) +MumID+BirthYearF,   #vm(ID, ainv) is the relatedness matricx, ide(ID) is the individual identity
                  data = vitdped,
                  na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units)) # include this line
anmodd3_old <- anmodd3_old$loglik

anmodd3 <- asreml(fixed = X25OHD3~1+Sex+Age+Age2+YearF, #etc etc
                  random = ~ vm(ID, grminv) + vm(ID2, grmreg) + ide(ID) +MumID+BirthYearF,   #vm(ID, ainv) is the relatedness matricx, ide(ID) is the individual identity
                  data = vitdped,
                  na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units)) # include this line

wald.asreml(anmodd3)
summary.asreml(anmodd3, coef = T)$coef.fixed # Fixed effects
summary.asreml(anmodd3, coef = T)$varcomp    # random effects
pchisq(2*(anmodd3$loglik - anmodd3_old), 1, lower.tail = F)


raneff3 <- summary.asreml(anmodd3, coef = T)$varcomp
raneff3 <- raneff3[-nrow(raneff3),]
raneff3 <- cbind(raneff3, asreml4pin(anmodd3))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Effect sizes                                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(GenABEL)
load("soayimp_genotype_data.RData")
x <- as.character.gwaa.data(soayimp[,"oar3_OAR18_68320039"]) %>% data.frame
x$ID <- row.names(x)
head(x)
str(x)
x$ID <- as.character(x$ID)

vitdped <- join(vitdped, x)
head(vitdped)

anmodd2 <- asreml(fixed = X25OHD2~1+Sex+Age+Age2+YearF+oar3_OAR18_68320039, #etc etc
                  random = ~ vm(ID, grminv) + ide(ID) +MumID+BirthYearF,   #vm(ID, ainv) is the relatedness matricx, ide(ID) is the individual identity
                  data = vitdped,
                  na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units)) # include this line

wald.asreml(anmodd2)
summary.asreml(anmodd2, coef = T)$coef.fixed # Fixed effects
summary.asreml(anmodd2, coef = T)$varcomp    # random effects
asreml4pin(anmodd2)

