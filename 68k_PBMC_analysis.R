# Title: 68k PBMC preprocessing and analysis
# Figures: Supplementary figure 2
# Author: Mahima Arunkumar

library(dplyr)
library(Seurat)
library(patchwork)

############################ Preprocessing ###############################

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "./input-files/PBMC_68k/filtered_matrices_mex/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc_68k <- CreateSeuratObject(counts = pbmc.data, project = "pbmc68k", min.cells = 3, min.features = 200)
#pbmc_68k
pbmc_68k[["percent.mt"]] <- PercentageFeatureSet(pbmc_68k, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc_68k, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter 
plot1 <- FeatureScatter(pbmc_68k, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc_68k, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Normalize and find variable features
pbmc <- subset(pbmc_68k, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scale data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


#Dimensionality reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(pbmc, reduction = "pca")

#Heatmap
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

#Elbow plot
ElbowPlot(pbmc) #PC 11 is what we choose

#Cluster cells
pbmc <- FindNeighbors(pbmc, dims = 1:11)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)


#UMAP
pbmc <- RunUMAP(pbmc, dims = 1:11)
DimPlot(pbmc, reduction = "umap", label = TRUE)

#save rds
saveRDS(pbmc, file = "./output-files/Task3/pbmc_68k.rds")

pbmc <- readRDS("./output-files/Task3/pbmc_68k.rds")

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
VlnPlot(pbmc, features = c("TRBC2", "CD3D", "CD3G", "CD3E", "IL7R", "GZMK", "LTB", "LEF1"), slot = "counts", log = TRUE)

# Violinplot for Bcells
VlnPlot(pbmc, features = c("PXK", "CD19", "MS4A1", "CD74", "CD79A", "IGHD"), slot = "counts", log = TRUE)

# Violinplot for DC
VlnPlot(pbmc, features = c("ZBTB46", "ITGAX", "CX3CR1", "ITGAM", "FCER1A", "CST3"), slot = "counts", log = TRUE)

# Violinplot for Macrophages
VlnPlot(pbmc, features = c("ID1", "FAR2", "IFITM1", "NFIL3", "NPL", "OTUB2"), slot = "counts", log = TRUE)

# Violinplot for Monocytes
VlnPlot(pbmc, features = c("LYZ", "CFP", "APOBEC3A", "CD7", "TET2", "FCGR3A", "MS4A7"), slot = "counts", log = TRUE)

# Violinplot for NKT
VlnPlot(pbmc, features = c("IL2RB", "NCAM1", "CD44", "IL12RB2", "CXCR4"), slot = "counts", log = TRUE)

# Violinplot for NK
VlnPlot(pbmc, features = c("KLRD1", "NKG7", "GNLY", "STYK1", "GZMA", "GZMB"), slot = "counts", log = TRUE)

# Violinplot for Plasma
VlnPlot(pbmc, features = c("MZB1", "SSR4", "IGHG1"), slot = "counts", log = TRUE)

# Violinplot for Platelet
VlnPlot(pbmc, features = c("PPBP", "ITGA2B", "GP1BA", "PF4", "TUBB1", "CCL5"), slot = "counts", log = TRUE)


#Single Dotplot containing all celltype annotation markers
marker_genes <- c("CD3E", "IL7R", "CD4", "CD8A", "CD8B", "MS4A1", "FCER1A", "CST3", "ID1", "FAR2", "IFITM1", "NFIL3", "NPL", "OTUB2",
                  "LYZ", "CFP", "FCGR3A", "MS4A7", "IL2RB", "NCAM1", "IL12RB2", "CXCR4", "KLRD1",
                  "NKG7", "GNLY", "GZMA", "GZMB", "PPBP")
DotPlot(object = pbmc, features = marker_genes)

#Celltype annotation and UMAP
new.cluster.ids <- c("T cells", "T cells", "T cells", "T cells", "NKT", "NK", "B cells", "DC", "Monocytes", "DC", "NK", "Platelet", "n.d.")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)

DefaultAssay(pbmc) <- 'RNA'

#Visualize IL1A and GSDME
FeaturePlot(pbmc, c("IL1A", "GSDME")) #not found at all
# you can plot raw counts as well
VlnPlot(pbmc, features = c("IL1A", "GSDME"), slot = "counts") #not found at all



