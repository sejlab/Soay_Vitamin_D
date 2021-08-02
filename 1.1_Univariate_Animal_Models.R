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

source("r/ASReml4.EstEffects.R")

attr(grminv, which = "INVERSE") <- TRUE

vitdped$Age2 <- vitdped$Age^2
vitdped$YearF <- factor(vitdped$Year)
vitdped$BirthYearF <- factor(vitdped$BirthYear)
vitdped$Weight2 <- vitdped$Weight^2


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Animal Models                                       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# total vitamin D 

anmodtotvitd <- asreml(fixed = Total.25D~1+Sex+Age+Age2+YearF, #etc etc
                       random = ~ vm(ID, grminv) + ide(ID) + MumID +BirthYearF,   #vm(ID, ainv) is the relatedness matricx, ide(ID) is the individual identity
                       data = vitdped,
                       na.action = na.method(x = "omit"),
                       residual = ~idv(units)) # include this line

wald.asreml(anmodtotvitd)
summary.asreml(anmodtotvitd, coef = T)$coef.fixed # Fixed effects
summary.asreml(anmodtotvitd, coef = T)$varcomp    # random effects
asreml4pin(anmodtotvitd)

raneff <- summary.asreml(anmodtotvitd, coef = T)$varcomp
raneff <- raneff[-nrow(raneff),]
raneff <- cbind(raneff, asreml4pin(anmodtotvitd))

# d2

anmodd2 <- asreml(fixed = X25OHD2~1+Sex+Age+Age2+YearF, #etc etc
                  random = ~ vm(ID, grminv) + ide(ID) +MumID+BirthYearF,   #vm(ID, ainv) is the relatedness matricx, ide(ID) is the individual identity
                  data = vitdped,
                  na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units)) # include this line

wald.asreml(anmodd2)
summary.asreml(anmodd2, coef = T)$coef.fixed # Fixed effects
summary.asreml(anmodd2, coef = T)$varcomp    # random effects

raneff2 <- summary.asreml(anmodd2, coef = T)$varcomp
raneff2 <- raneff2[-nrow(raneff2),]
raneff2 <- cbind(raneff2, asreml4pin(anmodd2))


# d3

anmodd3 <- asreml(fixed = X25OHD3~1+Sex+Age+Age2+YearF, #etc etc
                  random = ~ vm(ID, grminv) + ide(ID) +MumID+BirthYearF,   #vm(ID, ainv) is the relatedness matricx, ide(ID) is the individual identity
                  data = vitdped,
                  na.action = na.method(x = "omit", y = "omit"), residual = ~idv(units)) # include this line

wald.asreml(anmodd3)
summary.asreml(anmodd3, coef = T)$coef.fixed # Fixed effects
summary.asreml(anmodd3, coef = T)$varcomp    # random effects


raneff3 <- summary.asreml(anmodd3, coef = T)$varcomp
raneff3 <- raneff3[-nrow(raneff3),]
raneff3 <- cbind(raneff3, asreml4pin(anmodd3))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Format Animal Model results                         #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

ranef.results <- rbind(cbind(raneff, Age = "All", Model = "Total"),
                       cbind(raneff2, Age = "All", Model = "D2"),
                       cbind(raneff3, Age = "All", Model = "D3"))

fixef.results <- rbind(cbind(summary(anmodtotvitd, coef = T)$coef.fixed, FixEffect = row.names(summary(anmodtotvitd, coef = T)$coef.fixed), Age = "All", Model = "Total"),
                       cbind(summary(anmodd2, coef = T)$coef.fixed, FixEffect = row.names(summary(anmodd2, coef = T)$coef.fixed), Age = "All", Model = "D2"),
                       cbind(summary(anmodd3, coef = T)$coef.fixed, FixEffect = row.names(summary(anmodd3, coef = T)$coef.fixed), Age = "All", Model = "D3"))



head(ranef.results)

ranef.results$Effect <- as.character(ranef.results$Effect)
x <- data.frame(Effect = unique(ranef.results$Effect))
x
x$Effect2 <- c("Birth Year", "Maternal", "Additive Genetic", "Permanent Environment", "Residual")
x

ranef.results <- join(ranef.results, x)
ranef.results$Effect2 <- factor(ranef.results$Effect2, levels = rev(c("Additive Genetic", "Permanent Environment", "Maternal", "Birth Year",  "Residual")))

ggplot(ranef.results, aes(Model, Estimate, fill = Effect2)) +
  geom_bar(stat = "identity", col = "grey20") +
  scale_fill_brewer(palette = "Set3") +
  labs(fill = "", x = "Model", y = "Proportion of Phenotypic Variance") +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        strip.text.y = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14))

write.table(ranef.results, "results/0_Animal_Model_Random_Effects.txt", row.names = F, sep = "\t", quote = F)
write.table(fixef.results, "results/0_Animal_Model_Fixed_Effects.txt", row.names = F, sep = "\t", quote = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~ 3. Perform LRTs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

anmodtotvitd1 <- asreml(fixed = Total.25D~1+Sex+Age+Age2+YearF, #etc etc
                        random = ~ vm(ID, grminv) + ide(ID) + MumID,   #vm(ID, ainv) is the relatedness matricx, ide(ID) is the individual identity
                        data = vitdped,
                        residual = ~idv(units)) # include this line

REMLRT(anmodtotvitd1, anmodtotvitd)

anmodtotvitd2 <- asreml(fixed = Total.25D~1+Sex+Age+Age2+YearF, #etc etc
                        random = ~ vm(ID, grminv) + ide(ID) + BirthYearF,   #vm(ID, ainv) is the relatedness matricx, ide(ID) is the individual identity
                        data = vitdped,
                        residual = ~idv(units)) # include this line

REMLRT(anmodtotvitd2, anmodtotvitd)

anmodtotvitd3 <- asreml(fixed = Total.25D~1+Sex+Age+Age2+YearF, #etc etc
                        random = ~ vm(ID, grminv) + MumID + BirthYearF,   #vm(ID, ainv) is the relatedness matricx, ide(ID) is the individual identity
                        data = vitdped,
                        residual = ~idv(units)) # include this line

REMLRT(anmodtotvitd3, anmodtotvitd)


anmodtotvitd4 <- asreml(fixed = Total.25D~1+Sex+Age+Age2+YearF, #etc etc
                        random = ~ ID + MumID + BirthYearF,   #vm(ID, ainv) is the relatedness matricx, ide(ID) is the individual identity
                        data = vitdped,
                        residual = ~idv(units)) # include this line

REMLRT(anmodtotvitd4, anmodtotvitd)

rm(anmodtotvitd1, anmodtotvitd2, anmodtotvitd3, anmodtotvitd4)
gc()

# D2


anmodd2_1 <- asreml(fixed = X25OHD2~1+Sex+Age+Age2+YearF, #etc etc
                    random = ~ vm(ID, grminv) + ide(ID) + MumID,   #vm(ID, ainv) is the relatedness matricx, ide(ID) is the individual identity
                    data = vitdped,
                    residual = ~idv(units)) # include this line

REMLRT(anmodd2_1, anmodd2)

anmodd2_2 <- asreml(fixed = X25OHD2~1+Sex+Age+Age2+YearF, #etc etc
                    random = ~ vm(ID, grminv) + ide(ID) + BirthYearF,   #vm(ID, ainv) is the relatedness matricx, ide(ID) is the individual identity
                    data = vitdped,
                    residual = ~idv(units)) # include this line

REMLRT(anmodd2_2, anmodd2)

anmodd2_3 <- asreml(fixed = X25OHD2~1+Sex+Age+Age2+YearF, #etc etc
                    random = ~ vm(ID, grminv) + MumID + BirthYearF,   #vm(ID, ainv) is the relatedness matricx, ide(ID) is the individual identity
                    data = vitdped,
                    residual = ~idv(units)) # include this line

REMLRT(anmodd2_3, anmodd2)


anmodd2_4 <- asreml(fixed = X25OHD2~1+Sex+Age+Age2+YearF, #etc etc
                    random = ~ ID + MumID + BirthYearF,   #vm(ID, ainv) is the relatedness matricx, ide(ID) is the individual identity
                    data = vitdped,
                    residual = ~idv(units)) # include this line

REMLRT(anmodd2_4, anmodd2)



# D3


anmodd3_1 <- asreml(fixed = X25OHD3~1+Sex+Age+Age2+YearF, #etc etc
                    random = ~ vm(ID, grminv) + ide(ID) + MumID,   #vm(ID, ainv) is the relatedness matricx, ide(ID) is the individual identity
                    data = vitdped,
                    residual = ~idv(units)) # include this line

REMLRT(anmodd3_1, anmodd3)

anmodd3_2 <- asreml(fixed = X25OHD3~1+Sex+Age+Age2+YearF, #etc etc
                    random = ~ vm(ID, grminv) + ide(ID) + BirthYearF,   #vm(ID, ainv) is the relatedness matricx, ide(ID) is the individual identity
                    data = vitdped,
                    residual = ~idv(units)) # include this line

REMLRT(anmodd3_2, anmodd3)

anmodd3_3 <- asreml(fixed = X25OHD3~1+Sex+Age+Age2+YearF, #etc etc
                    random = ~ vm(ID, grminv) + MumID + BirthYearF,   #vm(ID, ainv) is the relatedness matricx, ide(ID) is the individual identity
                    data = vitdped,
                    residual = ~idv(units)) # include this line

REMLRT(anmodd3_3, anmodd3)


anmodd3_4 <- asreml(fixed = X25OHD3~1+Sex+Age+Age2+YearF, #etc etc
                    random = ~ ID + MumID + BirthYearF,   #vm(ID, ainv) is the relatedness matricx, ide(ID) is the individual identity
                    data = vitdped,
                    residual = ~idv(units)) # include this line

REMLRT(anmodd3_4, anmodd3)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~ 4. Means and Variance                                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

mean(vitdped$Total.25D)
mean(vitdped$X25OHD2)
mean(vitdped$X25OHD3)
var(vitdped$Total.25D)
var(vitdped$X25OHD2)
var(vitdped$X25OHD3)


