#  ------------------------------------------------------------------------

#  Vitamin D ms
#  AnMod structure for Susie for reg h2
#  AMS - 18/5/19

#  ------------------------------------------------------------------------

load("Vit_Reg_h2.RData", verbose = T)

library(plyr)
library(dplyr)
library(ggplot2) 
library(asreml)
library(asremlPlus)

source("r/ASReml4.EstEffects.R")

attr(grminv, which = "INVERSE") <- TRUE



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 5. Bivariate Models                                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

pedigree <- read.table("data/4_Updated_Pedigree_Feb2017.txt", header = T)
ainv <- ainverse(pedigree)

Vitd.D2D3 <- asreml(fixed  = cbind(X25OHD2,X25OHD3) ~ trait + trait:Sex + trait:Age + trait:Age2 + trait:YearF,
                    random = ~ us(trait):vm(ID, grminv) + us(trait):ide(ID),
                    residual   = ~ units:us(trait, init = NA),
                    data = vitdped, 
                    na.action = na.method(x = "omit", y = "omit"),
                    maxit = 20)

summary(Vitd.D2D3)
Vitd.D2D3$vparameters
vpredict(Vitd.D2D3, ra~V2/(sqrt(V1*V3)))
vpredict(Vitd.D2D3, ra~V5/(sqrt(V4*V6)))
vpredict(Vitd.D2D3, ra~V9/(sqrt(V8*V10)))


init_test <- c(0, 1, 1)
names(init_test) <- c("F", "U", "U")

Vitd.D2D3_ra0 <- asreml(fixed  = cbind(X25OHD2,X25OHD3) ~ trait + trait:Sex + trait:Age + trait:Age2 + trait:YearF,
                    random = ~ corgh(trait, init = init_test):vm(ID, grminv) + idh(trait):ide(ID),
                    residual   = ~ units:us(trait),
                    data = vitdped, 
                    na.action = na.method(x = "omit", y = "omit"),
                    maxit = 20)

summary(Vitd.D2D3_ra0)

2*(Vitd.D2D3$loglik - Vitd.D2D3_ra0$loglik)
pchisq(2*(Vitd.D2D3$loglik - Vitd.D2D3_ra0$loglik), 1, lower.tail = F)


init_test <- c(0.999, 1, 1)
names(init_test) <- c("F", "U", "U")

Vitd.D2D3_ra1 <- asreml(fixed  = cbind(X25OHD2,X25OHD3) ~ trait + trait:Sex + trait:Age + trait:Age2 + trait:YearF,
                        random = ~ corgh(trait, init = init_test):vm(ID, grminv) + idh(trait):ide(ID),
                        residual   = ~ units:us(trait),
                        data = vitdped, 
                        na.action = na.method(x = "omit", y = "omit"),
                        maxit = 20)

summary(Vitd.D2D3_ra1)

2*(Vitd.D2D3$loglik - Vitd.D2D3_ra1$loglik)
pchisq(2*(Vitd.D2D3$loglik - Vitd.D2D3_ra1$loglik), 1, lower.tail = F)
