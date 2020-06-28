#  ------------------------------------------------------------------------
#  Vitamin D ms
#  Figures and Tables
#  SEJ - 19/02/2020
#  ------------------------------------------------------------------------

library(plyr)
library(reshape2)
library(ggplot2)


load("Vit_Reg_h2.RData", verbose = T)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Format Animal Model Results                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

ranef.results <- read.table("results/0_Animal_Model_Random_Effects.txt", header = T, stringsAsFactors = F, sep = "\t")

ranef.results$Age <- factor(ranef.results$Age, levels = c("Lambs", "Adults", "All"))
ranef.results$Effect2 <- factor(ranef.results$Effect2, levels = rev(c("Additive Genetic", "Permanent Environment", "Maternal", "Birth Year",  "Residual")))
ranef.results$Model[which(ranef.results$Model == "Total")] <- "25(OH)D"
ranef.results$Model[which(ranef.results$Model == "D2")] <- "25(OH)D2"
ranef.results$Model[which(ranef.results$Model == "D3")] <- "25(OH)D3"

ggplot(ranef.results, aes(Model, Estimate, fill = Effect2)) +
  geom_bar(stat = "identity", col = "grey20") +
  scale_fill_brewer(palette = "Set3") +
  scale_x_discrete(labels = c("25(OH)D", expression("25(OH)D"[2]), expression("25(OH)D"[3]))) +
  labs(fill = "", x = "Vitamin D Measure", y = "Proportion of Phenotypic Variance") +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 12),
        strip.text.x = element_text (size = 12),
        strip.text.y = element_text (size = 12),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        legend.text = element_text(size = 12))
ggsave("figs/1_Animal_Model_All.png", width = 7, height =7 )


### Fixed effects

fixef.results <- read.table("results/0_Animal_Model_Fixed_Effects.txt", header = T, stringsAsFactors = F, sep = "\t")

head(fixef.results)

fixef.results$Model[which(fixef.results$Model == "Total")] <- "25(OH)D"
fixef.results$Model[which(fixef.results$Model == "D2")] <- "25(OH)D2"
fixef.results$Model[which(fixef.results$Model == "D3")] <- "25(OH)D3"

library(reshape2)
x <- melt(fixef.results, id.vars = c("Model", "FixEffect"))
x <- subset(x, variable != "z.ratio")
x <- dcast(x, FixEffect ~ Model + variable)

write.table(x, "results/1_Animal_Model_Fixef_Formatted.txt", row.names = F, sep = "\t", quote = F)

fixef.results <- subset(fixef.results, !FixEffect %in% c("(Intercept)", "Sex_1", "YearF_2011"))

y <- data.frame(FixEffect = unique(fixef.results$FixEffect))
y$FixEffect2 <- c("Year (2012)", "Year (2013)", "Year (2014)", "Year (2015)", "Year (2016)", "Age Sq.", "Age", "Sex (Male)")

fixef.results <- join(fixef.results, y)

ggplot(fixef.results, aes(FixEffect2, solution)) +
  geom_point() +
  geom_errorbar(aes(ymin = solution - std.error, ymax = solution + std.error), width = 0.1) +
  facet_wrap(~Model) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  labs(x = "", y = "Parameter Estimate") +
  coord_flip()
ggsave("figs/1_Animal_Model_Fixed.png", width = 7, height =5 )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Format GWAS Results                             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

gwas.results <- read.table("results/2_SEJ_GWAS_Results_600K.txt", header = T, stringsAsFactors = F)
gwas.results <- subset(gwas.results, Age == "All")


gwas.results$Model <- factor(gwas.results$Model, levels = c("25(OH)D", "25(OH)D2", "25(OH)D3"))

gwas.results <- arrange(gwas.results, Pc1df)

tophits <- NULL

for(i in unique(gwas.results$Age)){
  
  for(j in unique(gwas.results$Model)){
    
    tophits <- rbind(tophits,
                     subset(gwas.results, Age == i & Model == j)[1:30,])
    
  }
  
}

gwas.results <- arrange(gwas.results, Cumu)

chrinfo <- NULL

for(j in unique(gwas.results$Chromosome)){
  temp1 <- subset(gwas.results, Chromosome == j)
  temp2 <- data.frame(Chromosome = j,
                      Start = temp1[1,"Cumu"],
                      Stop = temp1[nrow(temp1),"Cumu"])
  chrinfo <- rbind(chrinfo, temp2)
}

chrinfo$Mid <- chrinfo$Start + ((chrinfo$Stop - chrinfo$Start)/2)

chrinfo$Chr2 <- c(0:10, "", 12, "", 14, "", 16, "", 18, "", 20, "", "", "", 24, "", "", "X")

png("figs/2_SEJ_GWAS_600K.png", width = 6, height = 12, units = "in", res = 300)

ggplot(gwas.results, aes(Cumu,-log10(Pc1df), col = factor(Chromosome %% 2))) +
  geom_point(size = 2, alpha = 0.4) +
  geom_hline(yintercept=-log10(1.28e-6),linetype=2, alpha = 0.6, size = 1) +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 16, hjust = 0),
        strip.text.y = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_rect(fill = "white")) +
  scale_x_continuous(breaks = chrinfo$Mid, labels = chrinfo$Chr2) +
  scale_colour_manual(values = c("darkgreen", "grey60")) +
  labs(x ="Chromosome", y = expression("-log"[10]*"P")) +
  facet_wrap(~Model, ncol = 1)
  
dev.off()

png("figs/2_SEJ_GWAS_600K_PP.png", width = 4, height = 12, units = "in", res = 300)

ggplot(gwas.results, aes(-log10(ExpP),-log10(Pc1df))) +
  geom_point(size = 2, alpha = 0.4) +
  geom_hline(yintercept=-log10(1.28e-6),linetype=2, alpha = 0.6, size = 1) +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 16, hjust = 0),
        strip.text.y = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        strip.background = element_rect(fill = "white")) +
  labs(x = expression("Expected -log"[10]*"P"), y = expression("Observed -log"[10]*"P")) +
  facet_wrap(~Model, ncol = 1) +
  geom_abline(slope = 1, intercept = 0)

dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Correlation between D2 & D3, Age plots          # 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

head(VITD)

ggplot(VITD, aes(X25OHD2, X25OHD3)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm") +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 16, hjust = 0),
        strip.text.y = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16)) +
  labs(x = expression("25(OH)D"[2]*" (nmol/l)"), y = expression("25(OH)D"[3]*" (nmol/l)"))
ggsave("figs/3_VitD_Corr.png", width = 6, height =6)


x <- subset(VITD, select = c(ID, Age, Sex, Total.25D, X25OHD2, X25OHD3))
x <- melt(x, id.vars = c("ID", "Age", "Sex"))
x$Sex2 <- ifelse(x$Sex == 1, "Female", ifelse(x$Sex == 2, "Male", NA))
x$variable2 <- as.character(x$variable)
x$variable2[which(x$variable == "Total.25D")] <- "25(OH)D"
x$variable2[which(x$variable == "X25OHD2")] <- "25(OH)D2"
x$variable2[which(x$variable == "X25OHD3")] <- "25(OH)D3"

ggplot(x, aes(Age, value, colour = Sex2)) +
  geom_point(alpha = 0.5) +
  stat_smooth() +
  #stat_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  theme_bw() +
  facet_wrap(~variable2, scales = "free") +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 16, hjust = 0),
        strip.text.y = element_text (size = 16),
        axis.title.y = element_text (size = 16, angle = 90),
        axis.title.x = element_text (size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = "top") +
  scale_colour_brewer(palette = "Set1") +
  labs(x = "Age", y = "Vitamin D concentration (nmol/l)", colour = "Sex")
ggsave("figs/3_VitD_Age.png", width = 10, height =4)


#~~~ t vs t+1


x <- subset*

