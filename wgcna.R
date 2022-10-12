# Title: Weighted correlation network analysis (WGCNA) network
# Author: Albert García López

pacman::p_load(tidyverse, devtools, WGCNA, Seurat, doParallel)
options(readr.show_col_types = FALSE, stringsAsFactors = FALSE)
outDir <- "path_to_dir"
dataDir <- "path_to_dir"
stopifnot(dir.exists(outDir) & dir.exists(dataDir))

# Load counts ##################################################################

# NOTE: Each row corresponds to a gene and each column to a sample or
# auxiliary information.

SeuratObj <- ReadMtx(
  mtx = file.path(dataDir, "matrix.mtx.gz"),
  cells = file.path(dataDir, "barcodes.tsv.gz"),
  features = file.path(dataDir, "features.tsv.gz")
) %>%
  CreateSeuratObject()

counts_mat <- GetAssayData(SeuratObj, slot = "counts") %>%
  as.matrix()

# Counts matrix contains NAs?
stopifnot(anyNA(counts_mat) == FALSE)

# Define metadata ##############################################################

metadata <- counts_mat %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  filter(gene_id == "IL1A") %>%
  pivot_longer(-gene_id, names_to = "cells", values_to = "values") %>%
  mutate(type = ifelse(values != 0, "IL1A_positive", "IL1A_negative")) %>%
  select(cells, values, type)

il1a_pos_cells <- metadata %>%
  filter(type == "IL1A_positive") %>%
  pull(cells)

# Filter counts for IL1A+ samples only and transpose the expression data for
# further analysis.

counts_filtered <- counts_mat %>%
  as.data.frame() %>%
  select(all_of(il1a_pos_cells)) %>%
  t() %>%
  as.data.frame()

# Data input, cleaning and pre-processing ######################################

# Checking data for excessive missing values and identification of outliers

gsg <- goodSamplesGenes(datExpr = counts_filtered)
genes_to_keep <- which(gsg$goodGenes == TRUE)
samples_to_keep <- which(gsg$goodSamples == TRUE)
counts_good <- counts_filtered[samples_to_keep, genes_to_keep]

# Automatic construction of the gene network and identification of modules #####

registerDoParallel(cores = 4)

# Pick soft-thresholding power: analysis of network topology #

sft <- pickSoftThreshold(
  data = counts_good,
  dataIsExpr = TRUE
)

# One-step network construction and module detection ###########################

net <- blockwiseModules(
  datExpr = counts_good,
  power = sft$powerEstimate,
  randomSeed = 123,
  maxBlockSize = ncol(counts_good),
  minModuleSize = 2,
  networkType = "unsigned",
  TOMType = "signed",
  saveTOMs = TRUE,
  reassignThreshold = 0,
  numericLabels = TRUE,
  verbose = 0,
  nThreads = 20
)

table(net$colors)

il1a_module <- net$colors %>%
  enframe() %>%
  filter(name == "IL1A") %>%
  pull(value)

il1a_related_genes <- net$colors %>%
  enframe() %>%
  filter(value == il1a_module) %>%
  pull(name)

# Annotation and enrichment analysis ###########################################

library(biomaRt)

il1a_related_genes_info <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id", "gene_biotype"),
  filters = "external_gene_name",
  values = il1a_related_genes,
  mart = useMart(
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = "hsapiens_gene_ensembl"
  )
) %>%
  as_tibble()

all_genes_info <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id", "gene_biotype"),
  filters = "external_gene_name",
  values = names(net$colors),
  mart = useMart(
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = "hsapiens_gene_ensembl"
  )
) %>%
  as_tibble()

unload("biomaRt")

all_genes_info_complete <- all_genes_info %>%
  drop_na(entrezgene_id)

set.seed(123)
go_enrichment <- GOenrichmentAnalysis(
  labels = net$colors[all_genes_info_complete$external_gene_name],
  entrezCodes = all_genes_info_complete$entrezgene_id,
  organism = "human",
  ontologies = c("CC", "BP", "MF"),
  evidence = "all",
  pCut = 1
)

go_enrichment_df <- go_enrichment$bestPTerms$`CC, BP, MF`$enrichment %>%
  as_tibble()

il1a_enrichment <- go_enrichment_df %>%
  filter(module == il1a_module)

write_tsv(
  x = il1a_enrichment,
  file = "file_name.tsv")
)

# Load TOM and filter for IL1A module ##########################################

load("block.RData")
TOM_mat <- as.matrix(TOM)
dimnames(TOM_mat) <- list(names(counts_good), names(counts_good))

# Get TOM for IL1A related genes

TOM_IL1A <- TOM_mat[il1a_related_genes_info$external_gene_name, il1a_related_genes_info$external_gene_name]

# Filter counts belonging to module IL1A

counts_il1a <- counts_good %>%
  as.data.frame() %>%
  dplyr::select(all_of(il1a_related_genes_info$external_gene_name)) %>%
  as.matrix()

# Make correlation matrix

il1a_cor <- cor(x = counts_il1a, method = "spearman") %>%
  as.data.frame() %>%
  rownames_to_column("Target") %>%
  pivot_longer(cols = -Target, names_to = "Source", values_to = "Corr") %>%
  filter(Target != Source) %>%
  dplyr::rename(fromAltName = Target, toAltName = Source) %>%
  mutate(`shared name` = paste(fromAltName, toAltName, sep = " (interacts with) ")) %>%
  dplyr::select(`shared name`, everything()) %>%
  as.data.frame()

write_tsv(
  x = il1a_cor,
  file = "file_name.tsv"
)

# Export to Cytoscape

cyt <- exportNetworkToCytoscape(
  adjMat = TOM_IL1A,
  edgeFile = paste0(outDir, "Cytoscape_edges_", paste(il1a_module, collapse = "_signed"), ".txt"),
  nodeFile = paste0(outDir, "Cytoscape_nodes_", paste(il1a_module, collapse = "_signed"), ".txt"),
  weighted = TRUE,
  threshold = 0,
  nodeNames = colnames(TOM_IL1A),
  altNodeNames = colnames(TOM_IL1A),
  nodeAttr = net$colors[il1a_related_genes_info$external_gene_name]
)

il1a_network_nodes_info <- il1a_cor %>%
  filter(fromAltName == "IL1A") %>%
  dplyr::select(`shared name` = toAltName, name = toAltName, Corr)

write_tsv(
  x = il1a_network_nodes_info,
  file = "file_name.tsv"
)
