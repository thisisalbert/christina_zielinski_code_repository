# Title: Differential Gene Correlation Analysis (DGCA) network
# Author: Albert García López

pacman::p_load(tidyverse, Seurat, DGCA, readxl)
outDir <- "path_to_dir"
dataDir <- "path_to_dir"
stopifnot(dir.exists(outDir) | dir.exists(dataDir))

# Load counts ##################################################################

SeuratObj <- ReadMtx(
  mtx = file.path(dataDir, "matrix.mtx.gz"),
  cells = file.path(dataDir, "barcodes.tsv.gz"),
  features = file.path(dataDir, "features.tsv.gz")
) %>%
  CreateSeuratObject()

counts_mat <- GetAssayData(SeuratObj, slot = "counts") %>%
  as.matrix()

# Sanity check

stopifnot(anyNA(counts_mat) == FALSE)

# Define design matrix #########################################################

cell_types <- counts_mat %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  filter(gene_id == "IL1A") %>%
  pivot_longer(-gene_id, names_to = "cells", values_to = "values") %>%
  mutate(type = ifelse(values != 0, "IL1A_positive", "IL1A_negative"))

il1a_pos_cells <- cell_types %>%
  filter(type == "IL1A_positive") %>%
  pull(cells)

il1a_neg_cells <- cell_types %>%
  filter(type == "IL1A_negative") %>%
  pull(cells)

design_mat <- counts_mat %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  as_tibble() %>%
  mutate(IL1A_positive = ifelse(cell_id %in% il1a_pos_cells, 1, 0)) %>%
  mutate(IL1A_negative = ifelse(cell_id %in% il1a_neg_cells, 1, 0)) %>%
  select(cell_id, IL1A_positive, IL1A_negative) %>%
  arrange(desc(IL1A_positive)) %>%
  column_to_rownames("cell_id") %>%
  as.matrix()

# Load DEGs ####################################################################

degs_all_upreg <- file.path(dataDir, "file_name.xlsx") %>%
  read_excel()

degs_all_downreg <- file.path(dataDir, "file_name.xlsx") %>%
  read_excel()

degs_cl1_upreg <- file.path(dataDir, "file_name.xlsx") %>%
  read_excel()

degs_cl1_downreg <- file.path(dataDir, "file_name.xlsx") %>%
  read_excel()

# DGCA #########################################################################

cl1 <- c(degs_cl1_upreg$gene, degs_cl1_downreg$gene)
all_cl <- c(degs_all_upreg$gene, degs_all_downreg$gene)
ref_genes <- c("^IL1A$|^CASP8$|^CASP3$|^GSDME|^RORC|^RORA$|^IL23|^NFKB|^NFAT")

registerDoParallel(cores = 5)

dgca_cl1 <- ddcorAll(
  inputMat = counts_mat,
  design = design_mat,
  compare = c("IL1A_positive", "IL1A_negative"),
  splitSet = cl1,
  nPairs = "all",
  adjust = "fdr",
  nPerms = 0,
  sigThresh = 0.05,
  corSigThresh = 0.05
)

dgca_all_cl <- ddcorAll(
  inputMat = counts_mat,
  design = design_mat,
  compare = c("IL1A_positive", "IL1A_negative"),
  splitSet = all_cl,
  nPairs = "all",
  adjust = "fdr",
  nPerms = 0,
  sigThresh = 0.05,
  corSigThresh = 0.05
)

# Significant only

signif_degs_cl1 <- dgca_cl1 %>%
  as_tibble() %>%
  filter(pValDiff_adj < 0.001) %>%
  filter(if_any(Gene1:Gene2, ~ str_detect(string = ., pattern = ref_genes)))

signif_degs_all_cl <- dgca_all_cl %>%
  as_tibble() %>%
  filter(pValDiff_adj < 0.001) %>%
  filter(if_any(Gene1:Gene2, ~ str_detect(string = ., pattern = ref_genes)))

# Export

write_tsv(
  x = signif_degs_cl1,
  file = "file_name.tsv"
)

write_tsv(
  x = signif_degs_all_cl,
  file = "file_name.tsv"
)
