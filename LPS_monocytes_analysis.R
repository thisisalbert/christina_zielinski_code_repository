# Title: LPS stimulated monocytes data preprocessing and analysis
# Figures: Supplementary figure 2
# Author: Mahima Arunkumar

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(forcats)
library("readxl")

############################ Preprocessing ###############################

# Load the LPS stim. dataset
LPS.data <- Read10X_h5("./input-files/LPS_DC_scRNAseq/GSM4819712_LPS_filtered_gene_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
LPS_stim <- CreateSeuratObject(counts = LPS.data, project = "LPS_stim_monocytes", min.cells = 3, min.features = 200)
LPS_stim[["percent.mt"]] <- PercentageFeatureSet(LPS_stim, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(LPS_stim, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter 
plot1 <- FeatureScatter(LPS_stim, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(LPS_stim, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Subset, normalize and find variable features
LPS <- subset(LPS_stim, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
LPS <- NormalizeData(LPS, normalization.method = "LogNormalize", scale.factor = 10000)
LPS <- NormalizeData(LPS)
LPS <- FindVariableFeatures(LPS, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(LPS), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(LPS)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scale data
all.genes <- rownames(LPS)
LPS <- ScaleData(LPS, features = all.genes)

#Dimensionality reduction
LPS <- RunPCA(LPS, features = VariableFeatures(object = LPS))

# Examine and visualize PCA results a few different ways
print(LPS[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(LPS, reduction = "pca")

#Heatmap and elbow plot
DimHeatmap(LPS, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(LPS, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(LPS) 

#Cluster cells
LPS <- FindNeighbors(LPS, dims = 1:11)
LPS <- FindClusters(LPS, resolution = 0.5)

#UMAP
LPS <- RunUMAP(LPS, dims = 1:11)
DimPlot(LPS, reduction = "umap", label = TRUE)

#save rds
saveRDS(LPS, file = "./output-files/Task2/pbmc_LPS.rds")

LPS <- readRDS("./output-files/Task2/pbmc_LPS.rds")

############################ Analysis ###############################

#Markers for cell type annotation taken from https://panglaodb.se/index.html

#T-cells: TRBC2, CD3D, CD3G, CD3E, IL7R, GZMK, LTB, LEF1
#B-cells: PXK, CD19, MS4A1, CD74, CD79A, IGHD
#DC: FCER1A, CST3,ZBTB46, ITGAX, CX3CR1, ITGAM
#Macrophages: ID1, FAR2, IFITM1, NFIL3, NPL, OTUB2
#Monocytes: LYZ, CFP, APOBEC3A, CD7, TET2, FCGR3A, MS4A7
#NKT: IL2RB, NCAM1, CD44, IL12RB2, CXCR4
#NK: KLRD1, NKG7, GNLY, STYK1, GZMA, GZMB
#Plasma: MZB1, SSR4, IGHG1
#Platelet: PPBP

# Violinplot for Tcells
VlnPlot(LPS, features = c("TRBC2", "CD3D", "CD3G", "CD3E", "IL7R", "GZMK", "LTB", "LEF1"), slot = "counts", log = TRUE)
#VlnPlot(LPS, features = c("CD4", "CD8A", "CD8B"), slot = "counts", log = TRUE)

# Violinplot for Bcells
VlnPlot(LPS, features = c("PXK", "CD19", "MS4A1", "CD74", "CD79A", "IGHD"), slot = "counts", log = TRUE)

# Violinplot for DC
VlnPlot(LPS, features = c("ZBTB46", "ITGAX", "CX3CR1", "ITGAM", "FCER1A", "CST3"), slot = "counts", log = TRUE)

# Violinplot for Macrophages
VlnPlot(LPS, features = c("ID1", "FAR2", "IFITM1", "NFIL3", "NPL", "OTUB2"), slot = "counts", log = TRUE)

# Violinplot for Monocytes
VlnPlot(LPS, features = c("CD14", "LYZ", "CFP", "APOBEC3A", "CD7", "TET2", "FCGR3A", "MS4A7"), slot = "counts", log = TRUE)

# Violinplot for NKT
VlnPlot(LPS, features = c("IL2RB", "NCAM1", "CD44", "IL12RB2", "CXCR4"), slot = "counts", log = TRUE)

# Violinplot for NK
VlnPlot(LPS, features = c("KLRD1", "NKG7", "GNLY", "STYK1", "GZMA", "GZMB"), slot = "counts", log = TRUE)

# Violinplot for Plasma
VlnPlot(LPS, features = c("MZB1", "SSR4", "IGHG1"), slot = "counts", log = TRUE)

# Violinplot for Platelet
VlnPlot(LPS, features = c("PPBP"), slot = "counts", log = TRUE)


#Single Dotplot containing all celltype annotation markers
marker_genes <- c("TRBC2", "CD3E", "LTB", "PPBP", "IL2RB", "NCAM1", "NKG7", "GNLY", "GZMB", "CD14", "LYZ", "CFP", "MS4A7", "ITGAX", "CST3", "FAR2", "IFITM1", "OTUB2", "MS4A1", "CD19")
DotPlot(object = LPS, features = marker_genes, cols = c("lightgrey", "blue"))

#Celltype annotation and UMAP
new.cluster.ids <- c("Monocytes", "Monocytes", "Monocytes", "Monocytes", "Monocytes", "Monocytes")
names(new.cluster.ids) <- levels(LPS)
LPS <- RenameIdents(LPS, new.cluster.ids)
DimPlot(LPS, reduction = "umap", pt.size = 0.5)
ggsave("./figures/Task2/Annotated_UMAP.pdf", width = 1283, height = 978, units = "px")

#Visualize IL1A and GSDME
DefaultAssay(LPS) <- 'RNA'
FeaturePlot(LPS, c("IL1A", "GSDME")) #not found at all
ggsave("IL1A_featureplot.pdf", width = 1083, height = 978, units = "px")
# you can plot raw counts as well
VlnPlot(LPS, features = c("IL1A", "GSDME"), slot = "counts") #not found at all
ggsave("./figures/Task2/IL1A_Vlnplot.pdf", width = 1383, height = 1178, units = "px")


#Add annotation IL1A+ vs. IL1A- cells to metadata
contains_IL1A_cells <- WhichCells(LPS, expression = IL1A > 0)
LPS$contains_IL1A<- ifelse(colnames(LPS) %in% contains_IL1A_cells, "IL1A+", "IL1A-")

#DEG IL1A+ vs. IL1A-
Idents(LPS) <- "contains_IL1A"
#reorder levels from  IL1A- vs. IL1A+ to  IL1A+ vs. IL1A-
Idents(LPS) <- factor(x = Idents(LPS), levels = c("IL1A+", "IL1A-"))
markers <- FindAllMarkers(LPS)
head(markers, n = 20)

#write to file
write.table(markers, file='./output-files/Task2/DEG_IL1A+_vs_IL1A-.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = TRUE)


#check again for different IL1A+ definition
contains_no_IL1A_cells <- WhichCells(LPS, expression = IL1A == 0)
contains_zero_pos_IL1A_cells <- WhichCells(LPS, expression = IL1A > 0 & IL1A < 1)
contains_1_IL1A_cells <- WhichCells(LPS, expression = IL1A >= 1 & IL1A < 2)
contains_2_IL1A_cells <- WhichCells(LPS, expression = IL1A >= 2 & IL1A < 3)

#Add annotation IL1A+ vs. IL1A- cells to metadata
LPS$contains_IL1A<- ifelse(colnames(LPS) %in% contains_no_IL1A_cells, "IL1A-", "IL1A+")
LPS$contains_zero_pos_IL1A<- ifelse(colnames(LPS) %in% contains_zero_pos_IL1A_cells, "IL1A_0", "IL1A_no0")
LPS$contains_1_IL1A<- ifelse(colnames(LPS) %in% contains_1_IL1A_cells, "IL1A_1", "IL1A_no1")
LPS$contains_2_IL1A<- ifelse(colnames(LPS) %in% contains_2_IL1A_cells, "IL1A_2", "IL1A_no2")

#get count of how many IL1A expression level categories we have 
table(LPS$contains_IL1A) #IL1A- = 3243 and IL1A+ = 2452
table(LPS$contains_zero_pos_IL1A) #IL1A_0 = 1211 and IL1A_no0 = 4484
table(LPS$contains_1_IL1A) #IL1A_1 = 1186 and IL1A_no1 = 4509
table(LPS$contains_2_IL1A) #IL1A_2 = 55 and IL1A_no2 = 5640

#Plot barplot
# Create data
data_IL1A_freq <- data.frame(
  Category=c("No IL1A","0 < IL1A < 1", "1 <= IL1A < 2", "2 <= IL1A") ,  
  "IL1A transcript expression level"=c(3243,1211,1186,55)
)
data_IL1A_freq

# Barplot
ggplot(data_IL1A_freq, aes(x=Category, y=IL1A.transcript.expression.level)) + 
  geom_bar(stat = "identity")

#New definition: IL1A+ as cells 2 <= IL1A < 3
df1 <- subset(x = LPS, subset = IL1A >= 2)
df2 <- subset(x = LPS, subset = IL1A == 0)
LPS_new <-  merge(x = df1, y = df2)   #merge


#Add annotation IL1A+ vs. IL1A- cells to metadata
contains_threshold_IL1A_cells <- WhichCells(LPS_new, expression = IL1A >= 2 & IL1A < 3)
LPS_new$IL1A_threshold_contains <- ifelse(colnames(LPS_new) %in% contains_threshold_IL1A_cells, "IL1A+", "IL1A-")

#DEG IL1A+ vs. IL1A-
Idents(LPS_new) <- "IL1A_threshold_contains"
#reorder levels from  IL1A- vs. IL1A+ to  IL1A+ vs. IL1A-
Idents(LPS_new) <- factor(x = Idents(LPS_new), levels = c("IL1A+", "IL1A-"))
markers <- FindAllMarkers(LPS_new)
head(markers, n = 20)

#write to file
write.table(markers, file='./output-files/Task2/Threshold_DEG_IL1A+_vs_IL1A-.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = TRUE)


##################### Enrichment analysis #####################

geneList <- markers$avg_log2FC
names(geneList) <- as.character(markers[["gene"]])
geneList = sort(geneList, decreasing = TRUE)


#create gene list for KEGG enrichment analysis
entrez_ids <- bitr(as.character(markers[["gene"]]), fromType = "SYMBOL",
                   toType = c("ENTREZID","SYMBOL"),
                   OrgDb = org.Hs.eg.db)

geneList_KEGG <- markers$avg_log2FC
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

gse_result$Description
gse_result$Description[str_detect(gse_result$Description, "inflam")]

#write to files
write.table(gse_result, file='./output-files/Task2/gse_result.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = TRUE)
write.table(gse_result_KEGG, file='./output-files/Task2/gse_result_KEGG.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = TRUE)
 

selected_pathways <- c("inflammatory response", "immune response" , "response to stress" ,
                       "cell communication" , "response to cytokine", "cell motility", "cellular response to chemical stimulus" ,
                       "defense response", "signaling", "signal transduction" )

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

dev.off()

ggsave("./figures/Task2/Enrichment_dotplot_IL1A+_vs_IL1A-.pdf", width = 2083, height = 1078, units = "px")

dotplot(gse_result, x = "NES",showCategory=33, split=".sign") + facet_wrap(.~.sign, scales = "free")

#Comparison of DEG list from DC and from T cells(Th17 no Salt) where we did IL1A+ vs Il1A-
#for input in DiVenn to check for overlap

#Th17 no Salt
#upregulated
Th17_DEG_list <- read.table("IL1A+_upregulated.txt", header = TRUE)
Th17_DEG_list$number <- ifelse(Th17_DEG_list$logfoldchange > 0, "1", "2")

Th17_DEG_list_final <- Th17_DEG_list[, c("gene", "number")]
write.table(Th17_DEG_list_final, file='DiVenn_Th17_noSalt.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)

#downregulated
Th17_DEG_list <- read.table("IL1A+_upregulated.txt", header = TRUE)
Th17_DEG_list$number <- ifelse(Th17_DEG_list$logfoldchange > 0, "1", "2")

Th17_DEG_list_final <- Th17_DEG_list[, c("gene", "number")]
write.table(Th17_DEG_list_final, file='DiVenn_Th17_noSalt.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)


#DEG DC
DEG_IL1A_DC <- read.table("DEG_IL1A+_vs_IL1A-_DC.txt", header = TRUE)
DEG_IL1A_DC$number <- ifelse(DEG_IL1A_DC$avg_log2FC > 0, "1", "2")

DEG_IL1A_DC_final <- DEG_IL1A_DC[, c("gene", "number")]
write.table(DEG_IL1A_DC_final, file='DiVenn_DC.txt', quote=FALSE, sep='\t', col.names = TRUE, row.names = FALSE)




