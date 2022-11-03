# Title: Th17 data for IL1A coexpression
# Figures: Supplementary figure 3
# Author: Mahima Arunkumar

library(dplyr)
library(readr)
library(Seurat)
library(clusterProfiler)
library(patchwork)
library(SeuratDisk)
library(SingleR)
library(SingleCellExperiment)
library(ggplot2)
library(gridExtra)
library(data.table)
library(scater)
library(loomR)
library(tibble)
library(scLink)
library(readxl)
library(writexl)
library(pheatmap)
library(corrplot)
library(igraph)
library(org.Hs.eg.db)

#Load in NaCl- h5seurat file
Th17_no_NaCl_new <- LoadH5Seurat("./input-files/Th17_noSalt/NaCl-_new.h5seurat")
exprs <- data.frame(GetAssayData(Th17_no_NaCl_new))

#Check if IL1A is significantly upregulated in cluster1 vs cluster rest (all other clusters) 
Idents(object = Th17_no_NaCl_new) <- "cluster_identity"
markers <- FindAllMarkers(Th17_no_NaCl_new)

#Find IL1A row
markers[markers$gene == "IL1A", 1:7]

#write to file
write.table(markers, file='./output-files/Task1/Cluster1_vs_rest.tsv', quote=FALSE, sep='\t', col.names = TRUE, row.names = TRUE)
write_xlsx(markers,"./output-files/Task1/Cluster1_vs_rest.xlsx")

#Check for raw counts cluster 1 vs rest
count.data <- GetAssayData(object = Th17_no_NaCl_new[["RNA"]], slot = "counts")
count.data <- as.matrix(x = count.data + 1)
new.seurat.object <- SetAssayData(
  object = Th17_no_NaCl_new,
  slot = "counts",
  new.data = count.data,
  assay = "RNA"
)

Idents(object = new.seurat.object) <- "cluster_identity"
markers <- FindAllMarkers(new.seurat.object, slot = "counts")

#Find IL1A row
markers[markers$gene == "IL1A", 1:4]


#Check IL1A+ vs IL1A-
Idents(object = Th17_no_NaCl_new) <- "contains_IL1A"
markers <- FindAllMarkers(Th17_no_NaCl_new)
#Find IL1A row
markers[markers$gene == "IL1A", 1:4]

#write to file
write.table(markers, file='./output-files/Task1/DEG_IL1A+_vs_IL1A-.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = TRUE)


#Look for IL1A coexpressing genes
Gene_expression_data_full <- as.data.frame(fread("./output-files/Task1/Gene_Count_per_Cell.tsv"))


#Find IL1A row
Gene_expression_data_full[Gene_expression_data_full$gene == "IL1A", 1:4]

# selecting columns where column is not 0 for Il1A expression
IL1A_expression <- Gene_expression_data_full[Gene_expression_data_full$gene == "IL1A", ] %>% select_if(~ any(. > 0))
IL1A_cell_names <- names(IL1A_expression) #cell names

#Subset Gene_expression_data_full to only these cells expressing IL1A
Gene_expression_data_reduced <- Gene_expression_data_full[, IL1A_cell_names]

#Remove all genes where in all columns the expression value is 0
resulting_df <- Gene_expression_data_reduced[apply(Gene_expression_data_reduced[,-1], 1, function(x) !all(x==0)),]
rownames(resulting_df) <- resulting_df$gene

#Drop the extra gene column
drop <- c("gene")
resulting_df = resulting_df[,!(names(resulting_df) %in% drop)]

#write to file
write.table(resulting_df, file='ScLink_resulting_df.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = TRUE)

#resulting dataframe resulting_df is used
resulting_df_transposed <- t(resulting_df)

#specify genes of interest; obtained primarily from https://maayanlab.cloud/archs4/gene/IL1A#correlation
genes <- c("IL1A", "PHLDA1","S100A3","S100A2","C10ORF55","ACOD1","EREG","CXCL11","RAB38","TNFSF9","CSF3","CSF2","SERPINE1","KRT6A","NIPAL4","SLAMF9","PHLDA2","KRT6B","TNFSF15","IL36G","INHBA","ELL2","MAFF","FOSL1","EHD4","CDCP1","GJB2","SPRR1B","CD274","COL17A1","IL20RB","BCL10","TXN","GPR87","PLAU","HIF1A","ICAM1","CYP27B1","IL36B","ZC3H12A","IL36RN","SMCO2","IER3","DUSP5","FGFBP1","EDN1","ANXA2","MMP1","LAMB3","CCL20","SPHK1","EMP1","LIF","PDCD1LG2","MMP10","IL1B","IL6","MMP14","KRT17","TNIP3","TMEM217","KRT16","ANXA2P2","CXCL8","PNLIPRP3","BTG3","EBI3","CXCL2","CXCL1","LAMC2","SLC7A11","AREG","CXCL3","C15ORF48","DRAM1","HRH1","IRAK2","RIPK2","SPRR2A","SPRR2D","SERPINB4","SERPINB2","RNASE7","PLAUR","F3","SERPINB8","RP11-443P15.2","SERPINB7","EPGN","CD109","B3GNT5","HBEGF","CSTB","ERO1A","CDKN1A","CST6","PCNPP3","PTGS2","ARNTL2","CCL7","UCN2")
#head(genes)

#specify genes of interest related to IL1A; obtained primarily from https://maayanlab.cloud/archs4/gene/IL1A#correlation and STRING interaction network
#genes <- c("IL1A", "IL1R2", "IL1RAP", "IL1R1", "IL10", "PHLDA1","S100A3","S100A2","C10ORF55","ACOD1","EREG","CXCL11","RAB38","TNFSF9","CSF3","CSF2","SERPINE1","KRT6A","NIPAL4","SLAMF9","PHLDA2","KRT6B","TNFSF15","IL36G","INHBA","ELL2","MAFF","FOSL1","EHD4","CDCP1","GJB2","SPRR1B","CD274","COL17A1","IL20RB","BCL10","TXN","GPR87","PLAU","HIF1A","ICAM1","CYP27B1","IL36B","ZC3H12A","IL36RN","SMCO2","IER3","DUSP5","FGFBP1","EDN1","ANXA2","MMP1","LAMB3","CCL20","SPHK1","EMP1","LIF","PDCD1LG2","MMP10","IL1B","IL6","MMP14","KRT17","TNIP3","TMEM217","KRT16","ANXA2P2","CXCL8","PNLIPRP3","BTG3","EBI3","CXCL2","CXCL1","LAMC2","SLC7A11","AREG","CXCL3","C15ORF48","DRAM1","HRH1","IRAK2","RIPK2","SPRR2A","SPRR2D","SERPINB4","SERPINB2","RNASE7","PLAUR","F3","SERPINB8","RP11-443P15.2","SERPINB7","EPGN","CD109","B3GNT5","HBEGF","CSTB","ERO1A","CDKN1A","CST6","PCNPP3","PTGS2","ARNTL2","CCL7","UCN2")
#head(genes)

#check if all gene names are found in resulting_df_transposed
func <- function(i){
  if (i %in% colnames(resulting_df_transposed) == FALSE){
    print("Meep, gene not found!\n");
  }
}
lapply(genes, func)

genes <- genes[genes %in% colnames(resulting_df_transposed)]
genes

#normalize count data
count.norm = sclink_norm(resulting_df_transposed, scale.factor = 1e6, filter.genes = FALSE, gene.names = genes)

#find correlation
corr = sclink_cor(expr = count.norm, ncores = 1)

#Plot result
corrplot(corr, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45,  col=colorRampPalette(c("blue","lightyellow","red"))(20))


#Plot heatmap
pheatmap(corr)
pheatmap(corr, cluster_cols = FALSE, cluster_rows = FALSE)

# Keep only correlations with relatively high correlation (>0.5)
corr[corr< 0.5] <- 0

#find network
networks = sclink_net(expr = count.norm, ncores = 1, lda = seq(0.5, 0.1, -0.05))
networks
networks$cor[1:3,1:3]


# Make an Igraph object from this matrix:
network <- graph_from_adjacency_matrix( corr, weighted=T, mode="undirected", diag=F)

plot(network, vertex.size=30, edge.color="black", edge.width=4, edge.arrow.size=1, edge.arrow.width=1, rescale = TRUE)
tkplot(network, vertex.size=30, vertex.color = "gold")


#enrichr for 100 IL1A coexpressing genes
i1<- c("IL1A", "PHLDA1","S100A3","S100A2","C10ORF55","ACOD1","EREG","CXCL11","RAB38","TNFSF9","CSF3","CSF2","SERPINE1","KRT6A","NIPAL4","SLAMF9","PHLDA2","KRT6B","TNFSF15","IL36G","INHBA","ELL2","MAFF","FOSL1","EHD4","CDCP1","GJB2","SPRR1B","CD274","COL17A1","IL20RB","BCL10","TXN","GPR87","PLAU","HIF1A","ICAM1","CYP27B1","IL36B","ZC3H12A","IL36RN","SMCO2","IER3","DUSP5","FGFBP1","EDN1","ANXA2","MMP1","LAMB3","CCL20","SPHK1","EMP1","LIF","PDCD1LG2","MMP10","IL1B","IL6","MMP14","KRT17","TNIP3","TMEM217","KRT16","ANXA2P2","CXCL8","PNLIPRP3","BTG3","EBI3","CXCL2","CXCL1","LAMC2","SLC7A11","AREG","CXCL3","C15ORF48","DRAM1","HRH1","IRAK2","RIPK2","SPRR2A","SPRR2D","SERPINB4","SERPINB2","RNASE7","PLAUR","F3","SERPINB8","RP11-443P15.2","SERPINB7","EPGN","CD109","B3GNT5","HBEGF","CSTB","ERO1A","CDKN1A","CST6","PCNPP3","PTGS2","ARNTL2","CCL7","UCN2")
#STRING interactions for IL1A
i2 <- c("IL1A", "IL1R1", "IL1B", "IL1R2", "IL1RAP", "IL6", "CSF2", "CSF3", "CXCL8", "IL10", "S100A13")
intersect(i1, i2)


########################################## Enrichment for IL1A co-expressed genes ###################################################

#set.seed(1000)

#Convert IL1A co-expressed genes to entrez id
IL1A_coexpressed_genes <- c("IL1A", "TNFSF9", "FOSL1", "CD109", "DRAM1", "PHLDA2", "HIF1A", "SPHK1", "MAFF", "HBEGF", "ICAM1", "PLAUR", "ZC3H12A")

genelist <- mapIds(x = org.Hs.eg.db,
                   keys = IL1A_coexpressed_genes,
                   column = "ENTREZID",
                   keytype = "SYMBOL",
                   multiVals = "first")
genelist <- as.data.frame(genelist)

ggo <- groupGO(gene     = genelist$genelist,
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",
               level    = 3,
               readable = TRUE)


ego <- enrichGO(genelist$genelist, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)

#write to file
write.table(ego, file='IL1A_coexpressing_enrichment_output.tsv', quote=FALSE, sep='\t', col.names = TRUE, row.names = TRUE)

dotplot(ego, showCategory=c("positive regulation of T cell proliferation", "cytokine production involved in inflammatory response", "regulation of cytokine production involved in inflammatory response", "production of molecular mediator involved in inflammatory response", "positive regulation of T cell proliferation", "cellular response to interleukin-1", "positive regulation of lymphocyte proliferation", "immature T cell proliferation", "mononuclear cell differentiation", "response to oxidative stress", "positive regulation of MAPK cascade"))
ggsave("IL1A_coexpression_enrichment_dotplot.pdf", width = 1883, height = 1700, units = "px")

