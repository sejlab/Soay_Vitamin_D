#  ------------------------------------------------------------------------
#  Vitamin D ms
#  Figures and Tables
#  SEJ - 19/02/2020
#  ------------------------------------------------------------------------

library(plyr)
library(reshape2)
library(ggplot2)
library(magrittr)


load("Vit_Reg_h2.RData", verbose = T)
load("soayimp_genotype_data.RData")

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

gwas.results$Model <- factor(gwas.results$Model, levels = c("25(OH)D", "25(OH)D2", "25(OH)D3"))

gwas.results <- arrange(gwas.results, Pc1df)

tophits <- NULL



for(j in unique(gwas.results$Model)){
  
  tophits <- rbind(tophits,
                   subset(gwas.results, Model == j)[1:30,])
  
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


#~~ Chromosome 18 region

chr18 <- subset(gwas.results, Model == "25(OH)D2" & Chromosome == 18 & Position > 64.3e6)

library(GenABEL)
library(reshape2)
soay18 <- soayimp[,chromosome(soayimp) == 18]
soay18 <- soay18[,map(soay18) > 64e6]
nsnps(soay18)

ld18 <- r2fast(soay18)
ld18 <- ld18 %>% melt
ld18 <- subset(ld18, value < 1.01)
ld18$Var1 <- as.character(ld18$Var1 )
ld18$Var2 <- as.character(ld18$Var2 )


topsnp <- c(which(ld18$Var1 == "oar3_OAR18_68320039"), which(ld18$Var2 == "oar3_OAR18_68320039")) %>% unique
ld18 <- ld18[topsnp,]
ld18$SNP.Name <- ifelse(ld18$Var1 == "oar3_OAR18_68320039", ld18$Var2, ld18$Var1)
ld18 <- subset(ld18, select = c(SNP.Name, value))
chr18 <- join(chr18, ld18)
chr18



ggplot(chr18, aes(Position/1e6,-log10(Pc1df), fill = value)) +
  geom_point(shape = 21, size = 3, alpha = 0.8, colour = "black") +
  geom_hline(yintercept=-log10(1.28e-6), linetype=2, alpha = 0.6, size = 1) +
  theme_bw() +
  theme(axis.text.x  = element_text (size = 12),
        axis.text.y  = element_text (size = 14),
        strip.text.x = element_text (size = 14, hjust = 0),
        strip.text.y = element_text (size = 14),
        axis.title.y = element_text (size = 14, angle = 90),
        axis.title.x = element_text (size = 14),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Chromosome 18 Position (Mb)", y = expression("Observed -log"[10]*"P"), fill = "LD (r2)")
ggsave("figs/2_Chr_18_Region.png", width = 10, height =4 )


#~~ Plot genes

chr18genes <- read.table("results/3_Top_Genes_v4.txt", header = T, stringsAsFactors = F, sep = "\t")
chr18genes <- subset(chr18genes, gene_biotype == "protein_coding")

chr18orthos <- read.table("results/3_Top_Orthologous_Regions_v4.txt", header = T, stringsAsFactors = F, sep = "\t")
chr18orthos <- subset(chr18orthos, Species == "hsapiens")
chr18orthos <- subset(chr18orthos, select = c(Gene.stable.ID, Gene.name.1, Chromosome.scaffold.name, Gene.start..bp.,Gene.end..bp.)) %>% unique
names(chr18orthos) <- c("ensembl_gene_id", "Human_Orthologue", "human_chromosome_name", "human_start_position", "human_end_position")


chr18genes <- join(chr18genes, chr18orthos)
chr18genes$gene_size <- chr18genes$end_position - chr18genes$start_position

chr18genes$colour <- "black"
chr18genes$colour[which(chr18genes$Human_Orthologue %in% c("DLK1", "PPP1R13B"))] <- "red"

chr18genes$jitter <- rep(1:2, length.out = nrow(chr18genes))
chr18genes$jitter[which(chr18genes$ensembl_gene_id == "ENSOARG00000007297")] <- 0
chr18genes <- chr18genes[-grep("IGHG", chr18genes$Human_Orthologue)[2:4],]

chr18genes_plot <- subset(chr18genes, select = c(ensembl_gene_id, colour, jitter, start_position, end_position))
chr18genes_plot <- melt(chr18genes_plot, id.vars = c("ensembl_gene_id", "colour", "jitter"))

ggplot(chr18genes_plot, aes(value, jitter, group = ensembl_gene_id, colour = colour)) +
  geom_line(size = 5) +
  scale_colour_identity() +
  theme_bw() +
  theme(panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank()) +
  coord_cartesian(ylim = c(-0.5, 2.5))
ggsave("figs/2_Chr_18_Genes.png", width = 10, height =2 )

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



  
  