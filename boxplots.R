#Title: gene signature comparison in cluster 1 vs. rest, IL1A+ vs. IL1A- and GSDME+ vs. GSDME-
#Figures: 2b & 7c,h
#Author: Laurens Lehner


library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(ggplot2)
library(EnhancedVolcano)
library(ggpubr)
library(org.Hs.eg.db)
library(data.table)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(DOSE)
library(GOSemSim)
library(dplyr)
library(pathview)
library(ggprism)

#load data as preprocessed in scanpy and converted to .h5seurat
all_cells <- LoadH5Seurat("C:\\NaCl-_new.h5seurat")
#current cell annotations
names(x = all_cells[[]])

#add average gene expression scores to seurat object

#pro-inflammatory gene set
pro_inflammatory <- list(c('CD276','IL2RA','PTPRK','CSF2','TOX2','IL12RB2','MYB','IL17F','GJB2','IFNG','IL23R','CCR7','CXCR3','IL2','IL22','IL6','IL1R1'))
all_cells <- AddModuleScore(
  object = all_cells,
  features = pro_inflammatory,
  ctrl = 5,
  name = 'pro_inflam'
)
#anti-inflammatory gene set
anti_inflammatory <- list(c('IL10','PTGDS','LRRC32','SLC7A8','CXCR6','CD180','CCR9','CTLA4','CCR5','IKZF3','P2RX7','GZMA','CCR1','NKG7','PDCD1','IL17A','HDAC4','PRNP','BATF','LGMN','LTB','ZBP1','PIK31P1','MAF','HOPX','DUSP6','CD69','TGFB1','ITGB7','PRDM1','JARID2','KLF2','BACH2','TBKBP1','NFATC2','E2F2'))
all_cells <- AddModuleScore(
  object = all_cells,
  features = anti_inflammatory,
  ctrl = 5,
  name = 'anti_inflam'
)
#T cell proliferation
T_cell_prolif <- list(c('CD24','RASGRP1','EBI3','LILRB2','BTN2A2','GPNMB','CEBPB','TNFSF13B','LILRB1','MALT1','LILRB4','RIPK3','BTN3A1','GLMN','HHLA2','CORO1A','VSIG4','TNFRSF13C','CLC','TMIGD2','CR1','RC3H1','IL23R','CTLA4','CTNNB1','CTPS1','CLECL1','CD55','DHPS','DLG1','AGER','DOCK2','EFNB1','AIF1','ELF4','EPO','ERBB2','FKBP1B','FOXJ1','TMEM131L','ICOSLG','NCSTN','SCRIB','CADM1','IL27','ABL1','FYN','IL4I1','LRRC32','PTPN22','PLA2G2D','GJA1','TNFRSF21','CCDC88B','LGALS9B','PYCARD','CD274','ANXA1','NCKAP1L','CD209','HLA-A','HLA-DMB','HLA-DPA1','HLA-DPB1','HLA-DRB1','HLA-E','HLA-G','HMGB1','HES1','CLEC4G','IGF1','IGF2','IGFBP2','IHH','IL1A','IL1B','IL2','IL2RA','IL4','IL6','IL6ST','IL10','IL12B','IL12RB1','IL15','IL18','IDO1','IRF1','JAK2','JAK3','ARG1','ARG2','LEP','LGALS3','LGALS9','LMO1','MIR181C','MIR21','MIR30B','CD46','MSN','NCK1','PNP','P2RX7','PAWR','FOXP3','ZBTB7B','IL23A','PIK3CG','IL20RB','WNT4','SASH3','RC3H2','PPP3CB','LMBR1L','PRKAR1A','PRKCQ','PRNP','CRTAM','PSMB10','TWSG1','PELI1','SH3RF1','GPAM','PTPN6','ZP4','PTPRC','BAX','SELENOK','RAC2','IL21','RPS3','RPS6','CCL5','CCL19','XCL1','SDC4','VSIR','SFTPD','PLA2G2F','SHH','MARCHF7','RASAL3','SLAMF1','BMP4','LGALS7B','SLC7A1','LGALS9C','SLC11A1','SOS1','SOS2','SPN','SPTA1','STAT5B','SYK','PRDX2','TFRC','TGFBR2','TNFRSF1B','TP53','TRAF6','TNFSF4','CCR2','TNFRSF4','TYK2','SCGB1A1','VCAM1','ZAP70','ZP3','VTCN1','ARMC5','PDCD1LG2','CD276','NDFIP1','DOCK8','CASP3','ITCH','MAD1L1','NCK2','CARD11','HAVCR2','PDE5A','CBLB','TNFSF14','TNFSF9','TNFRSF14','RIPK2','FADD','CCND3','TNFSF18','DNAJA3','CD1D','CD3E','CD6','DLG5','CD28','CD80','CD86','TNFSF8','IL27RA','MAPK8IP1','CD40LG','CD70','CD81','CD151'))
all_cells <- AddModuleScore(
  object = all_cells,
  features = T_cell_prolif,
  ctrl = 5,
  name = 'T_cell_prolif'
)
#positive regulation of T cell proliferation
positive_reg_T_cell_prolif <- list(c('CD24','HHLA2','TMIGD2','IL23R','AGER','EPO','ICOSLG','PYCARD','HMGB1','IGF1','IGF2','IGFBP2','IL2','IL2RA','IL12B','IL12RB1','IL18','JAK3','MIR21','MIR30B','IL23A','GPAM','RPS3','SLAMF1','STAT5B','TNFSF9','FADD','IL27RA'))
all_cells <- AddModuleScore(
  object = all_cells,
  features = positive_reg_T_cell_prolif,
  ctrl = 5,
  name = 'pos_reg_T_cell_prolif'
)

#scores have been added to cell annotations
names(x = all_cells[[]])

#subset data set by cluster identity
Idents(all_cells) <- "cluster_identity"
cluster1 <- subset(x = all_cells, idents = "cluster 1")
rest <- subset(x = all_cells, idents = "rest")

#save pro-inflammatory signature scores for cluster 1 and rest respectively
cluster1.pro_inflam1.values <- cluster1[["pro_inflam1"]][,1]
rest.pro_inflam1.values <- rest[["pro_inflam1"]][,1]

#save anti-inflammatory signature scores for cluster 1 and rest respectively
cluster1.anti_inflam1.values <- cluster1[["anti_inflam1"]][,1]
rest.anti_inflam1.values <- rest[["anti_inflam1"]][,1]

#save positive regulation of T cell proliferation signature scores for cluster 1 and rest respectively
cluster1.T_cell_prolif1.values <- cluster1[["T_cell_prolif1"]][,1]
rest.T_cell_prolif1.values <- rest[["T_cell_prolif1"]][,1]

#test cluster 1 against rest for positive regulation of T cell proliferation signature score with wilcoxon
wilcox_pval_clusters_T_cell_prolif <- wilcox.test(cluster1.T_cell_prolif1.values, rest.T_cell_prolif1.values, alternative = "greater")

#test cluster 1 against rest for pro-inflammatory signature score with wilcoxon
wilcox_pval_clusters_pro_inflam <- wilcox.test(cluster1.pro_inflam1.values, rest.pro_inflam1.values, alternative = "greater")

#test cluster 1 against rest for anti-inflammatory signature score with wilcoxon
wilcox_pval_clusters_anti_inflam <- wilcox.test(cluster1.anti_inflam1.values, rest.anti_inflam1.values, alternative = "less")



#subset data set by IL1A expression
Idents(all_cells) <- "contains_IL1A"
IL1A_positive <- subset(x = all_cells, idents = "IL1A+")
IL1A_negative <- subset(x = all_cells, idents = "IL1A-")

#save pro-inflammatory signature scores for IL1A+ and IL1A- respectively
IL1A_positive.pro_inflam1.values <- IL1A_positive[["pro_inflam1"]][,1]
IL1A_negative.pro_inflam1.values <- IL1A_negative[["pro_inflam1"]][,1]

#save anti-inflammatory signature scores for IL1A+ and IL1A- respectively
IL1A_positive.anti_inflam1.values <- IL1A_positive[["anti_inflam1"]][,1]
IL1A_negative.anti_inflam1.values <- IL1A_negative[["anti_inflam1"]][,1]

#save positive regulation of T cell proliferation signature scores for IL1A+ and IL1A- respectively
IL1A_positive.T_cell_prolif1.values <- IL1A_positive[["T_cell_prolif1"]][,1]
IL1A_negative.T_cell_prolif1.values <- IL1A_negative[["T_cell_prolif1"]][,1]

#test IL1A+ against IL1A- for positive regulation of T cell proliferation signature score with wilcoxon
wilcox_pval_IL1A_T_cell_prolif <- wilcox.test(IL1A_positive.T_cell_prolif1.values, IL1A_negative.T_cell_prolif1.values, alternative = "greater")

#test IL1A+ against IL1A- for pro-inflammatory signature score with wilcoxon
wilcox_pval_IL1A_pro_inflam <- wilcox.test(IL1A_positive.pro_inflam1.values, IL1A_negative.pro_inflam1.values, alternative = "greater")

#test IL1A+ against IL1A- for pro-inflammatory signature score with wilcoxon
wilcox_pval_IL1A_anti_inflam <- wilcox.test(IL1A_positive.anti_inflam1.values, IL1A_negative.anti_inflam1.values, alternative = "less")



#subset data by GSDME expression
Idents(all_cells) <- "contains_GSDME"
GSDME_positive <- subset(x = all_cells, idents = "GSDME+")
GSDME_negative <- subset(x = all_cells, idents = "GSDME-")

#save T cell proliferation signature scores for IL1A+ and IL1A- respectively
GSDME_positive.pos_reg_T_cell_prolif1.values <- GSDME_positive[["pos_reg_T_cell_prolif1"]][,1]
GSDME_negative.pos_reg_T_cell_prolif1.values <- GSDME_negative[["pos_reg_T_cell_prolif1"]][,1]

#test IL1A+ against IL1A- for pro-inflammatory signature score with wilcoxon
wilcox_pval_GSDME_pos_reg_T_cell_prolif <- wilcox.test(GSDME_positive.pos_reg_T_cell_prolif1.values, GSDME_negative.pos_reg_T_cell_prolif1.values, alternative = "greater")

#show p-values for all comparisons
wilcox_pval_clusters_pro_inflam
wilcox_pval_clusters_anti_inflam
wilcox_pval_IL1A_pro_inflam
wilcox_pval_IL1A_anti_inflam

wilcox_pval_GSDME_pos_reg_T_cell_prolif
wilcox_pval_clusters_T_cell_prolif
wilcox_pval_IL1A_T_cell_prolif

#build dataframes to plot values
df_cluster1_pro_inflam <- data.frame(cluster1[["pro_inflam1"]][,1], variable = "cluster 1")
colnames(df_cluster1_pro_inflam) <- c("score", "cluster_identity")
df_rest_pro_inflam <- data.frame(rest[["pro_inflam1"]][,1], variable = "rest clusters")
colnames(df_rest_pro_inflam) <- c("score", "cluster_identity")
df_clusters_pro_inflam <- rbind(df_cluster1_pro_inflam, df_rest_pro_inflam)

#violin plot with added boxplot showing average gene set expression for each cell by group
plot_1 <- ggplot(df_clusters_pro_inflam, aes(x=cluster_identity, y=score)) + 
  geom_violin(aes(fill=cluster_identity)) +
  geom_point(position = position_jitter(seed = 1, width = 0.25)) +
  geom_boxplot(width=0.2, outlier.color = NA, aes(fill=cluster_identity)) +
  #scale_y_continuous(limits = quantile(df_clusters_pro_inflam$score, c(0.01, 0.99))) +
  xlab("") + ylab("module score (pro-inflammatory signature)")
plot_1 + scale_fill_manual(values=c("#E25C33","#5284EC")) +
  scale_color_grey() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("wilcoxon p = 2.084391e-22") +
  theme(plot.title = element_text(hjust = 0.5))

#same as above but different gene set
df_cluster1_anti_inflam <- data.frame(cluster1[["anti_inflam1"]][,1], variable = "cluster 1")
colnames(df_cluster1_anti_inflam) <- c("score", "cluster_identity")
df_rest_anti_inflam <- data.frame(rest[["anti_inflam1"]][,1], variable = "rest clusters")
colnames(df_rest_anti_inflam) <- c("score", "cluster_identity")
df_clusters_anti_inflam <- rbind(df_cluster1_anti_inflam, df_rest_anti_inflam)

plot_2 <- ggplot(df_clusters_anti_inflam, aes(x=cluster_identity, y=score)) + 
  geom_violin(aes(fill=cluster_identity)) +
  geom_point(position = position_jitter(seed = 1, width = 0.25)) +
  geom_boxplot(width=0.2, outlier.color = NA, aes(fill=cluster_identity)) +
  #scale_y_continuous(limits = quantile(df_clusters_anti_inflam$score, c(0.01, 0.99))) +
  xlab("") + ylab("module score (anti-inflammatory signature)")
plot_2 + scale_fill_manual(values=c("#E25C33","#5284EC")) +
  scale_color_grey() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("wilcoxon p = 8.839685e-05") +
  theme(plot.title = element_text(hjust = 0.5))

#same as above but different gene set
df_IL1A_positive_pro_inflam <- data.frame(IL1A_positive[["pro_inflam1"]][,1], variable = "IL1A+")
colnames(df_IL1A_positive_pro_inflam) <- c("score", "IL1A")
df_IL1A_negative_pro_inflam <- data.frame(IL1A_negative[["pro_inflam1"]][,1], variable = "IL1A-")
colnames(df_IL1A_negative_pro_inflam) <- c("score", "IL1A")
df_IL1A_pro_inflam <- rbind(df_IL1A_positive_pro_inflam, df_IL1A_negative_pro_inflam)

plot_3 <- ggplot(df_IL1A_pro_inflam, aes(x=IL1A, y=score)) + 
  geom_violin(aes(fill=IL1A, x = factor(IL1A, level = c("IL1A+", "IL1A-")))) +
  geom_point(position = position_jitter(seed = 1, width = 0.25)) +
  geom_boxplot(width=0.2, outlier.color = NA, aes(fill=IL1A)) +
  #scale_y_continuous(limits = quantile(df_IL1A_pro_inflam$score, c(0.01, 0.99))) +
  xlab("") + ylab("module score (pro-inflammatory signature)")
plot_3 + scale_fill_manual(values=c("#5284EC","#E25C33")) +
  scale_color_grey() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("wilcoxon p = 9.783565e-18") +
  theme(plot.title = element_text(hjust = 0.5))


#same as above but different gene set
df_IL1A_positive_anti_inflam <- data.frame(IL1A_positive[["anti_inflam1"]][,1], variable = "IL1A+")
colnames(df_IL1A_positive_anti_inflam) <- c("score", "IL1A")
df_IL1A_negative_anti_inflam <- data.frame(IL1A_negative[["anti_inflam1"]][,1], variable = "IL1A-")
colnames(df_IL1A_negative_anti_inflam) <- c("score", "IL1A")
df_IL1A_anti_inflam <- rbind(df_IL1A_positive_anti_inflam, df_IL1A_negative_anti_inflam)

plot_4 <- ggplot(df_IL1A_anti_inflam, aes(x=IL1A, y=score)) + 
  geom_violin(aes(fill=IL1A, x = factor(IL1A, level = c("IL1A+", "IL1A-")))) +
  geom_point(position = position_jitter(seed = 1, width = 0.25)) +
  geom_boxplot(width=0.2, outlier.color = NA, aes(fill=IL1A)) +
  #scale_y_continuous(limits = quantile(df_IL1A_anti_inflam$score, c(0.01, 0.99))) +
  xlab("") + ylab("module score (anti-inflammatory signature)")
plot_4 + scale_fill_manual(values=c("#5284EC","#E25C33")) +
  scale_color_grey() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("wilcoxon p = 0.0001324167") +
  theme(plot.title = element_text(hjust = 0.5))


#same as above but different gene set
df_GSDME_positive_pos_reg_T_cell_prolif <- data.frame(GSDME_positive[["T_cell_prolif1"]][,1], variable = "GSDME+")
colnames(df_GSDME_positive_pos_reg_T_cell_prolif) <- c("score", "GSDME")
df_GSDME_negative_pos_reg_T_cell_prolif <- data.frame(GSDME_negative[["pos_reg_T_cell_prolif1"]][,1], variable = "GSDME-")
colnames(df_GSDME_negative_pos_reg_T_cell_prolif) <- c("score", "GSDME")
df_GSDME_pos_reg_T_cell_prolif <- rbind(df_GSDME_positive_pos_reg_T_cell_prolif, df_GSDME_negative_pos_reg_T_cell_prolif)

df_GSDME_pos_reg_T_cell_prolif$GSDME <- factor(df_GSDME_pos_reg_T_cell_prolif$GSDME,levels = c('GSDME+','GSDME-'),ordered = TRUE)
plot_5 <- ggplot(df_GSDME_pos_reg_T_cell_prolif, aes(x=GSDME, y=score)) + 
  geom_violin(aes(fill=GSDME)) +
  geom_point(position = position_jitter(seed = 1, width = 0.25)) +
  geom_boxplot(width=0.2, outlier.color = NA, aes(fill=GSDME)) +
  #scale_y_continuous(limits = quantile(df_GSDME_pos_reg_T_cell_prolif$score, c(0.01, 0.99))) +
  xlab("") + ylab("module score (positive regulation of T cell proliferation signature)")
plot_5 + scale_fill_manual(values=c("#E25C33","#5284EC")) +
  scale_color_grey() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("wilcoxon p = 1.238e-06") +
  theme(plot.title = element_text(hjust = 0.5))


#same as above but different gene set
df_cluster1_T_cell_prolif <- data.frame(cluster1[["T_cell_prolif1"]][,1], variable = "cluster 1")
colnames(df_cluster1_T_cell_prolif) <- c("score", "cluster_identity")
df_rest_T_cell_prolif <- data.frame(rest[["T_cell_prolif1"]][,1], variable = "rest clusters")
colnames(df_rest_T_cell_prolif) <- c("score", "cluster_identity")
df_clusters_T_cell_prolif <- rbind(df_cluster1_T_cell_prolif, df_rest_T_cell_prolif)

plot_6 <- ggplot(df_clusters_T_cell_prolif, aes(x=cluster_identity, y=score)) + 
  geom_violin(aes(fill=cluster_identity)) +
  geom_point(position = position_jitter(seed = 1, width = 0.25)) +
  geom_boxplot(width=0.2, outlier.color = NA, aes(fill=cluster_identity)) +
  #scale_y_continuous(limits = quantile(df_clusters_T_cell_prolif$score, c(0.01, 0.99))) +
  xlab("") + ylab("module score (T cell proliferation signature)")
plot_6 + scale_fill_manual(values=c("#E25C33","#5284EC")) +
  scale_color_grey() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("wilcoxon p = 0.0001124") +
  theme(plot.title = element_text(hjust = 0.5))


#same as above but different gene set
df_IL1A_positive_T_cell_prolif <- data.frame(IL1A_positive[["T_cell_prolif1"]][,1], variable = "IL1A+")
colnames(df_IL1A_positive_T_cell_prolif) <- c("score", "IL1A")
df_IL1A_negative_T_cell_prolif <- data.frame(IL1A_negative[["T_cell_prolif1"]][,1], variable = "IL1A-")
colnames(df_IL1A_negative_T_cell_prolif) <- c("score", "IL1A")
df_IL1A_T_cell_prolif <- rbind(df_IL1A_positive_T_cell_prolif, df_IL1A_negative_T_cell_prolif)

plot_7 <- ggplot(df_IL1A_T_cell_prolif, aes(x=IL1A, y=score)) + 
  geom_violin(aes(fill=IL1A, x = factor(IL1A, level = c("IL1A+", "IL1A-")))) +
  geom_point(position = position_jitter(seed = 1, width = 0.25)) +
  geom_boxplot(width=0.2, outlier.color = NA, aes(fill=IL1A)) +
  #scale_y_continuous(limits = quantile(df_IL1A_T_cell_prolif$score, c(0.01, 0.99))) +
  xlab("") + ylab("module score (T cell proliferation signature)")
plot_7 + scale_fill_manual(values=c("#5284EC","#E25C33")) +
  scale_color_grey() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ggtitle("wilcoxon p = < 2.2e-16") +
  theme(plot.title = element_text(hjust = 0.5))


#build data frame to convert to to excel sheet for boxplot values
#figure
figure_all <- c("Fig. 2b","Fig. 2b","Fig. 2b","Fig. 2b","Fig. 2b","Fig. 2b","Fig. 2b","Fig. 2b","Fig. 7c","Fig. 7c","Fig. 7h","Fig. 7h","Fig. 7h","Fig. 7h")
#signature
signature_all <- c("pro-inflammatory","pro-inflammatory","pro-inflammatory","pro-inflammatory","anti-inflammatory","anti-inflammatory","anti-inflammatory","anti-inflammatory","positive regulation of T cell proliferation","positive regulation of T cell proliferation","T cell proliferation", "T cell proliferation","T cell proliferation","T cell proliferation")

#group
group_all <- c("cluster 1","rest", "IL1A+","IL1A-","cluster 1","rest","IL1A+","IL1A-","GSDME+","GSDME-", "cluster 1","rest","IL1A+","IL1A-")
#median
median_all <- c(median(cluster1.pro_inflam1.values),median(rest.pro_inflam1.values),
                median(IL1A_positive.pro_inflam1.values), median(IL1A_negative.pro_inflam1.values),
                median(cluster1.anti_inflam1.values), median(rest.anti_inflam1.values),
                median(IL1A_positive.anti_inflam1.values), median(IL1A_negative.anti_inflam1.values),
                median(GSDME_positive.pos_reg_T_cell_prolif1.values), median(GSDME_negative.pos_reg_T_cell_prolif1.values),
                median(cluster1.T_cell_prolif1.values), median(rest.T_cell_prolif1.values),
                median(IL1A_positive.T_cell_prolif1.values), median(IL1A_negative.T_cell_prolif1.values))
#min
min_all <- c(min(cluster1.pro_inflam1.values),min(rest.pro_inflam1.values),
                min(IL1A_positive.pro_inflam1.values), min(IL1A_negative.pro_inflam1.values),
                min(cluster1.anti_inflam1.values), min(rest.anti_inflam1.values),
                min(IL1A_positive.anti_inflam1.values), min(IL1A_negative.anti_inflam1.values),
                min(GSDME_positive.pos_reg_T_cell_prolif1.values), min(GSDME_negative.pos_reg_T_cell_prolif1.values),
                min(cluster1.T_cell_prolif1.values), min(rest.T_cell_prolif1.values),
                min(IL1A_positive.T_cell_prolif1.values), min(IL1A_negative.T_cell_prolif1.values))
#max
max_all <- c(max(cluster1.pro_inflam1.values),max(rest.pro_inflam1.values),
                max(IL1A_positive.pro_inflam1.values), max(IL1A_negative.pro_inflam1.values),
                max(cluster1.anti_inflam1.values), max(rest.anti_inflam1.values),
                max(IL1A_positive.anti_inflam1.values), max(IL1A_negative.anti_inflam1.values),
                max(GSDME_positive.pos_reg_T_cell_prolif1.values), max(GSDME_negative.pos_reg_T_cell_prolif1.values),
                max(cluster1.T_cell_prolif1.values), max(rest.T_cell_prolif1.values),
                max(IL1A_positive.T_cell_prolif1.values), max(IL1A_negative.T_cell_prolif1.values))
#Q1
Q1_all <- c(quantile(cluster1.pro_inflam1.values, 0.25),quantile(rest.pro_inflam1.values, 0.25),
                quantile(IL1A_positive.pro_inflam1.values, 0.25), quantile(IL1A_negative.pro_inflam1.values, 0.25),
                quantile(cluster1.anti_inflam1.values, 0.25), quantile(rest.anti_inflam1.values, 0.25),
                quantile(IL1A_positive.anti_inflam1.values, 0.25), quantile(IL1A_negative.anti_inflam1.values, 0.25),
                quantile(GSDME_positive.pos_reg_T_cell_prolif1.values, 0.25), quantile(GSDME_negative.pos_reg_T_cell_prolif1.values, 0.25),
                quantile(cluster1.T_cell_prolif1.values, 0.25), quantile(rest.T_cell_prolif1.values, 0.25),
                quantile(IL1A_positive.T_cell_prolif1.values, 0.25), quantile(IL1A_negative.T_cell_prolif1.values, 0.25))

#Q3
Q3_all <- c(quantile(cluster1.pro_inflam1.values, 0.75),quantile(rest.pro_inflam1.values, 0.75),
            quantile(IL1A_positive.pro_inflam1.values, 0.75), quantile(IL1A_negative.pro_inflam1.values, 0.75),
            quantile(cluster1.anti_inflam1.values, 0.75), quantile(rest.anti_inflam1.values, 0.75),
            quantile(IL1A_positive.anti_inflam1.values, 0.75), quantile(IL1A_negative.anti_inflam1.values, 0.75),
            quantile(GSDME_positive.pos_reg_T_cell_prolif1.values, 0.75), quantile(GSDME_negative.pos_reg_T_cell_prolif1.values, 0.75),
            quantile(cluster1.T_cell_prolif1.values, 0.75), quantile(rest.T_cell_prolif1.values, 0.75),
            quantile(IL1A_positive.T_cell_prolif1.values, 0.75), quantile(IL1A_negative.T_cell_prolif1.values, 0.75))
#IQR for whiskers
IQR_cluster1_pro_inf <- (quantile(cluster1.pro_inflam1.values, 0.75) - quantile(cluster1.pro_inflam1.values, 0.25))*1.5
IQR_rest_pro_inf <- (quantile(rest.pro_inflam1.values, 0.75) - quantile(rest.pro_inflam1.values, 0.25))*1.5

IQR_IL1A_positive_pro_inf <- (quantile(IL1A_positive.pro_inflam1.values, 0.75) - quantile(IL1A_positive.pro_inflam1.values, 0.25))*1.5
IQR_IL1A_negative_pro_inf <- (quantile(IL1A_negative.pro_inflam1.values, 0.75) - quantile(IL1A_negative.pro_inflam1.values, 0.25))*1.5

IQR_cluster1_anti_inf <- (quantile(cluster1.anti_inflam1.values, 0.75) - quantile(cluster1.anti_inflam1.values, 0.25))*1.5
IQR_rest_anti_inf <- (quantile(rest.anti_inflam1.values, 0.75) - quantile(rest.anti_inflam1.values, 0.25))*1.5

IQR_IL1A_positive_anti_inf <- (quantile(IL1A_positive.anti_inflam1.values, 0.75) - quantile(IL1A_positive.anti_inflam1.values, 0.25))*1.5
IQR_IL1A_negative_anti_inf <- (quantile(IL1A_negative.anti_inflam1.values, 0.75) - quantile(IL1A_negative.anti_inflam1.values, 0.25))*1.5

IQR_cluster1_T_cell_prolif <- (quantile(cluster1.T_cell_prolif1.values, 0.75) - quantile(cluster1.T_cell_prolif1.values, 0.25))*1.5
IQR_rest_T_cell_prolif <- (quantile(rest.T_cell_prolif1.values, 0.75) - quantile(rest.T_cell_prolif1.values, 0.25))*1.5

IQR_GSDME_positive_pos_reg_T_cell_prolif <- (quantile(GSDME_positive.pos_reg_T_cell_prolif1.values, 0.75) - quantile(GSDME_positive.pos_reg_T_cell_prolif1.values, 0.25))*1.5
IQR_GSDME_negative_pos_reg_T_cell_prolif <- (quantile(GSDME_negative.pos_reg_T_cell_prolif1.values, 0.75) - quantile(GSDME_negative.pos_reg_T_cell_prolif1.values, 0.25))*1.5

IQR_IL1A_positive_T_cell_prolif <- (quantile(IL1A_positive.T_cell_prolif1.values, 0.75) - quantile(IL1A_positive.T_cell_prolif1.values, 0.25))*1.5
IQR_IL1A_negative_T_cell_prolif <- (quantile(IL1A_negative.T_cell_prolif1.values, 0.75) - quantile(IL1A_negative.T_cell_prolif1.values, 0.25))*1.5

#whisker up
whisker_up <- c(quantile(cluster1.pro_inflam1.values, 0.75) + IQR_cluster1_pro_inf ,quantile(rest.pro_inflam1.values, 0.75) + IQR_rest_pro_inf,
            quantile(IL1A_positive.pro_inflam1.values, 0.75) + IQR_IL1A_positive_pro_inf, quantile(IL1A_negative.pro_inflam1.values, 0.75) + IQR_IL1A_positive_pro_inf,
            quantile(cluster1.anti_inflam1.values, 0.75) + IQR_cluster1_anti_inf, quantile(rest.anti_inflam1.values, 0.75) + IQR_rest_anti_inf,
            quantile(IL1A_positive.anti_inflam1.values, 0.75) + IQR_IL1A_positive_anti_inf, quantile(IL1A_negative.anti_inflam1.values, 0.75) + IQR_IL1A_negative_anti_inf,
            quantile(GSDME_positive.pos_reg_T_cell_prolif1.values, 0.75) + IQR_GSDME_positive_pos_reg_T_cell_prolif, quantile(GSDME_negative.pos_reg_T_cell_prolif1.values, 0.75) + IQR_GSDME_negative_pos_reg_T_cell_prolif,
            quantile(cluster1.T_cell_prolif1.values, 0.75) + IQR_cluster1_T_cell_prolif, quantile(rest.T_cell_prolif1.values, 0.75) + IQR_rest_T_cell_prolif,
            quantile(IL1A_positive.T_cell_prolif1.values, 0.75) + IQR_IL1A_positive_T_cell_prolif, quantile(IL1A_negative.T_cell_prolif1.values, 0.75) + IQR_IL1A_negative_T_cell_prolif)
#whisker down
whisker_down <- c(quantile(cluster1.pro_inflam1.values, 0.25) - IQR_cluster1_pro_inf ,quantile(rest.pro_inflam1.values, 0.25) - IQR_rest_pro_inf,
                quantile(IL1A_positive.pro_inflam1.values, 0.25) - IQR_IL1A_positive_pro_inf, quantile(IL1A_negative.pro_inflam1.values, 0.25) - IQR_IL1A_positive_pro_inf,
                quantile(cluster1.anti_inflam1.values, 0.25) - IQR_cluster1_anti_inf, quantile(rest.anti_inflam1.values, 0.25) - IQR_rest_anti_inf,
                quantile(IL1A_positive.anti_inflam1.values, 0.25) - IQR_IL1A_positive_anti_inf, quantile(IL1A_negative.anti_inflam1.values, 0.25) - IQR_IL1A_negative_anti_inf,
                quantile(GSDME_positive.pos_reg_T_cell_prolif1.values, 0.25) - IQR_GSDME_positive_pos_reg_T_cell_prolif, quantile(GSDME_negative.pos_reg_T_cell_prolif1.values, 0.25) - IQR_GSDME_negative_pos_reg_T_cell_prolif,
                quantile(cluster1.T_cell_prolif1.values, 0.25) - IQR_cluster1_T_cell_prolif, quantile(rest.T_cell_prolif1.values, 0.25) - IQR_rest_T_cell_prolif,
                quantile(IL1A_positive.T_cell_prolif1.values, 0.25) - IQR_IL1A_positive_T_cell_prolif, quantile(IL1A_negative.T_cell_prolif1.values, 0.25) - IQR_IL1A_negative_T_cell_prolif)
#data frame
boxplot_df <- data.frame("figure" = figure_all,
                         "signature" = signature_all,
                         "group" = group_all,
                         "min" = min_all,
                         "max" = max_all,
                         "center" = median_all,
                         "Q1" = Q1_all,
                         "Q3" = Q3_all,
                         "whisker_up" = whisker_up,
                         "whisker_down" = whisker_down)
library(writexl)