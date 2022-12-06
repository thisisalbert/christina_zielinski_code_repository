#Title: clusterprofiler dotplot of enriched GO terms using upregulated genes in IL1A+ Th17 cells
#Figures: supplementary Fig. 3c
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
library(readxl)
library(writexl)
library(msigdbr)
library(forcats)
library(scales)
library(anndata)
library(tidyverse)

#gene list taken from differential expression analysis in scanpy
deg_from_excel <- read_excel("signif_upreg_genes_IL1A+.xlsx")
DEGs <- as.data.table(deg_from_excel)
DEGs[,"...1":=NULL]

#filter ribosomal genes so the enriched gene sets are no dominated by these (not always necessary)
DEGs <- dplyr::filter(DEGs, !grepl(paste(c("RPL","RPS"),collapse="|"),gene))

#create gene list for GO enrichment analysis
geneList <- DEGs$logfoldchange
names(geneList) <- as.character(DEGs[["gene"]])
geneList = sort(geneList, decreasing = TRUE)

#create gene list for KEGG enrichment analysis
entrez_ids <- bitr(as.character(DEGs[["gene"]]), fromType = "SYMBOL",
                                 toType = c("ENTREZID","SYMBOL"),
                                 OrgDb = org.Hs.eg.db)
geneList_KEGG <- DEGs$logfoldchange
names(geneList_KEGG) <- as.character(entrez_ids[["ENTREZID"]])
geneList_KEGG <- sort(geneList_KEGG, decreasing = TRUE)

#set seed to keep it reproducible
set.seed(1000)
gse_result <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              minGSSize    = 10,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              pAdjustMethod = "BH",
              keyType = "SYMBOL",
              seed = TRUE,
              verbose      = FALSE)

gse_result_KEGG <- gseKEGG(geneList     = geneList_KEGG,
               organism     = 'hsa',
               minGSSize    = 10,
               pvalueCutoff = 1,
               verbose      = FALSE)

#variable to select pathways to show on plot
selected_pathways <- gse_result$Description[1:25]

ggplot(gse_result, showCategory = selected_pathways, 
  aes(NES, fct_reorder(Description, NES),xend=0, yend = Description)) + 
  #geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE), direction = -1) +
  scale_size_continuous(range=c(2, 10)) +
  #theme_minimal() + 
  xlab("normalized enrichment score (NES)") +
  ylab(NULL) +
  ggtitle("enriched GO Terms (BP) in IL1A+ (upregulated genes)") +
  theme(plot.title = element_text(size = 10))