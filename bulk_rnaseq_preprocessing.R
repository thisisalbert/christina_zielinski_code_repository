# Title: Bulk RNA-seq pre-processing (Geo2RNAseq pipeline)
# Author: Albert García López

# Load libraries and import functions ##########################################

library("Geo2RNAseq")
library("ggplot2")
library("ggfortify")
library("data.table")
library("gtools")
library("Rsamtools")
library("tidyverse")

# Set WD

setwd("path_to_dir")
stopifnot(dir.exists(getwd()))

# To resume previous workspaces

dir() %>% stringr::str_subset(string = ., pattern = "RData")

# Global Settings ##############################################################

# Output directories locations

indexDir <- paste0(getwd(), "/index")
fastqDir <- paste0(getwd(), "/fastq")
qualDir  <- paste0(getwd(), "/quality")
mapDir   <- paste0(getwd(), "/mapping")
countDir <- paste0(getwd(), "/counting")
tabDir   <- paste0(getwd(), "/result_tables")
plotDir  <- paste0(getwd(), "/result_plots")
degDir   <- paste0(getwd(), "/diff_exp_genes")

# Create output directories

output_directories <- c(
  indexDir, fastqDir, qualDir, mapDir, countDir, tabDir, plotDir, degDir
)

for (i in output_directories){
  if (dir.exists(i) == FALSE){
    dir.create(i, recursive = TRUE)
  } else {
    message(paste0(i, " is already there!"))
  }
}

# Sanity check

for (i in output_directories){
  if (dir.exists(i) == TRUE){
    message(paste0("Directory [", i, "] exists!"))
  } else {
    message("Something has gone wrong! Repeat previous step please!")
  }
}

# Allocate CPU's cores usage

MAX_CPUS <- 20

# Force overwritting of files

# NOTE:
# If TRUE, existing files will be overwritten without questioning.
# If FALSE, most methods will skip step for existing output files.
# In these cases, the method will return 'call[x] = "NOT USED"'

FORCE_OVERWRITE <- TRUE

# Is it single or paired-end?

paired <- TRUE

# Files ########################################################################

# NOTE: TopHat2/HISAT2 are looking for the genome file in the index directory
# move it there or create a link if the genome is located in its own directory

genome <- "path_to_fasta_reference_file"
anno   <- "path_to_gtf_reference_file"

# Sanity check

stopifnot(file.exists(genome) == TRUE)
stopifnot(file.exists(anno) == TRUE)

# Build index with make_HiSat2_index() or make_Tophat_index() ##################
# NOTE: add something like "index" after the index directory.

index <- make_HiSat2_index(
  genomeFile = genome,
  customOut = paste0(indexDir, "/index")
)

# Load Raw FASTQ files #########################################################

raw_fastq_files <- list.files(
  path = "raw_data/",
  full.names = TRUE,
  pattern = "fq.gz",
  recursive = TRUE
)

stopifnot(file.exists(raw_fastq_files) == TRUE) # QC

if (paired == TRUE){
  message(paste0("We have ", length(raw_fastq_files)/2, " samples (PAIRED)"))
} else {
  message(paste0("We have ", length(raw_fastq_files), " samples (SINGLE)"))
}

writeLines("Working with files:")
print(raw_fastq_files)

# Quality Control 1 (before Trimming) ##########################################

qualityRawDir <- file.path(qualDir, "raw")

if (length(dir(qualityRawDir)) > 0){
  warning(
    paste0("Directory \"", qualityRawDir, "\" not empty! Overwrite? ", FORCE_OVERWRITE),
    immediate. = TRUE
  )
}

if (FORCE_OVERWRITE || length(dir(qualityRawDir)) == 0){

  writeLines("FastQC - raw ...")
  fq_res <- run_FastQC(
    raw_fastq_files,
    outDir = qualityRawDir,
    cpus = MAX_CPUS,
    extend = TRUE
  )

}

# Trimming #####################################################################

trimmedDir <- fastqDir
windowsizetrimming <- 15
qualcuttrimming <- 25
phred <- "-phred33"
leading <- 3
trailing <- 3
minlen <- 30

if (length(list.files(trimmedDir, pattern = "\\.trimo(\\.pe)*\\.fastq$")) > 0){
  warning(
    paste0("Directory \"", trimmedDir, "\" not empty! Overwrite? ", FORCE_OVERWRITE),
    immediate. = TRUE
  )
}

trimming_res <- run_Trimmomatic(
  files      = raw_fastq_files,
  is.paired  = paired,
  outDir     = trimmedDir,
  cpus       = MAX_CPUS,
  windowsize = windowsizetrimming,
  qualcut    = qualcuttrimming,
  phred      = phred,
  leading    = leading,
  trailing   = trailing,
  minlen     = minlen
)

# Number of reads after trimming based on Trimmomatic output

number_raw <- trimming_res$input
number_trimmed <- trimming_res$surviving

if (paired) {
  trimmed_fastq_files <- asPairVector(trimming_res$files)
} else {
  trimmed_fastq_files <- trimming_res$files
}

# Sanity check (Can be ignored during experimental sessions)

stopifnot(length(trimmed_fastq_files) == length(raw_fastq_files))
names(number_raw) <- basename(raw_fastq_files)
names(number_trimmed) <- basename(trimmed_fastq_files)
fastq_files <- trimmed_fastq_files

# SortMeRNA (removes rRNA reads) ###############################################

filterrRNA <- TRUE
sortmeDir <- file.path(fastqDir, "sortmerna")

if (filterrRNA && length(dir(sortmeDir)) > 0){
  warning(
    paste0("Directory \"", sortmeDir, "\" not empty! Overwrite? ", FORCE_OVERWRITE),
    immediate. = TRUE
  )
}

if (filterrRNA) {

  sortmerna_res <- run_SortMeRNA(
    files     = fastq_files,
    outDir    = fastqDir,
    mode      = "fast",
    paired    = paired,
    cpus      = MAX_CPUS
  )

  non_rrna_files <- sortmerna_res$files
  number_nonrRNA <- sortmerna_res$nonrrna

}

fastq_files <- if (filterrRNA) non_rrna_files else trimmed_fastq_files

# Check whether sample/file order has changed:

stopifnot(order(trimmed_fastq_files) == order(non_rrna_files))

# Quality Control 2 (After Trimming) ###########################################

qualitySubDir <- if (filterrRNA) {
  file.path(qualDir, "non_rrna")
} else {
  file.path(qualDir, "trimmed")
}

if (length(dir(qualitySubDir)) > 0){
  warning(
    paste0("Directory \"" , qualitySubDir, "\" not empty! overwrite? ",
    FORCE_OVERWRITE)
  )
}

if (FORCE_OVERWRITE || length(dir(qualitySubDir)) == 0){

  writeLines("FastQC - trimmed ...")
  fq_res <- run_FastQC(
    fastq_files,
    outDir = qualitySubDir,
    cpus = MAX_CPUS,
    extend = FALSE
  )

}

# Mapping ######################################################################

# NOTE: Be sure that "index_map" is associated with your "index" object. Not
# doing this step as it is mentioned will end up in a HISAT2 error.

index_map <- index
anno_map  <- anno
mapper <- "hisat2"
bamDir <- file.path(mapDir, "bamfiles")
samDir <- file.path(mapDir, "samfiles")
convert_sam <- TRUE # Should SAM files be converted to BAM files if BAM files are not found?

if (mapper == "tophat2"){

  mapping_res <- run_Tophat(
    fastq_files,
    index        = index_map,
    outDir       = mapDir,
    is.paired    = paired,
    anno         = anno_map,
    addArgs      = NA,      # "-g 2 --b2-very-sensitive --no-coverage-search",
    cpus         = MAX_CPUS,
    worker       = if (paired) 5 else 10,
    use.existing = TRUE,
    overwrite = TRUE
  )
  bam_files <- mapping_res$files

} else {

  if (FORCE_OVERWRITE || length(dir(samDir)) == 0){

    mapping_res <- run_Hisat2(
      files     = fastq_files,
      index     = index_map,
      outDir    = mapDir,
      is.paired = paired,
      addArgs   = "",
      splice    = TRUE,
      cpus      = MAX_CPUS,
      overwrite = TRUE
    )

    ## NOTE: use that only if the run_Hisat2() parameter 'as.bam' is FALSE.
    # bam_files <- sam_to_bam(mapping_res$mapFiles, sort=TRUE, bamDir = bamDir)

    bam_files <- mapping_res$files

  } else {
    if (convert_sam) {
      sam_files <- list.files(samDir, pattern = "\\.sam$", full.names = TRUE)
      if (length(sam_files) != number_samples)
        stop(
        paste("Invalid number of SAM files. Expected", number_samples, "but got",
        length(sam_files))
      )
      bam_files <- sam_to_bam(mapping_res$mapFiles, sort=TRUE, bamDir = bamDir)
    } else
      bam_files <- list.files(bamDir, pattern = "\\.bam$", full.names = TRUE)

    if (length(bam_files) != number_samples)
      warning(
        paste(
          "Found",
          length(bam_files),
          "BAM files, but expected",
          number_samples,
          ". Maybe SAMtools was interrupted. Try to convert again? ",
          convert_sam),
        immediate. = TRUE)
  }
}

# Order BAM files properly

bam_files_sorted <- unlist(
  mclapply(bam_files, sortBAMs, overwrite = FORCE_OVERWRITE, mc.cores = MAX_CPUS)
)

# Counting #####################################################################

anno_count <- anno

if (length(dir(countDir)) > 0){
  warning(
    paste("Directory", countDir, "not empty! Overwrite?", FORCE_OVERWRITE),
    immediate. = TRUE
  )
}

if (FORCE_OVERWRITE || length(dir(countDir)) == 0){
  res_counting <- run_featureCounts(
    files             = bam_files_sorted,
    annotation        = anno_count,
    isGTF             = TRUE,
    IDtype            = "gene_id",
    featureType       = "exon",
    outDir            = countDir,
    isPairedEnd       = paired,
    allowMultiOverlap = FALSE,
    cpus              = MAX_CPUS
  )
} else {
  res_counting <- read.featureCounts.files(countDir)
}

# Save results

counts <- res_counting$counts
countfile <- res_counting$countFile
sumfile <- res_counting$sumFile
info <- res_counting$anno
lib_sizes <- res_counting$summary[,1]
gene_lengths <- info$Length
names(gene_lengths) <- info$GeneID

# Rename sample names for "lib_sizes" and "counts"

names(lib_sizes) <- str_replace_all(
  string = names(lib_sizes),
  pattern = "_\\d.*",
  replacement = ""
)
colnames(counts) <- str_replace_all(
  string = colnames(counts),
  pattern = "_\\d.*",
  replacement = ""
)

# Write counts to CSV file

write_count_table(
  file   = file.path(tabDir, "counts"),
  counts = counts
)

# Write counts to XLS file

write_count_table(
  file       = file.path(tabDir, "counts"),
  counts     = counts,
  as.xls     = TRUE,
  sheetNames = basename(getwd())
)

# Check for rRNA contamination #################################################

detect_high_coverage(counts, lib_sizes)

# SAMtools #####################################################################

flagstatDir <- file.path(mapDir, "flagstats")
flag_files <- make_flagstats(bam_files_sorted, flagstatDir, MAX_CPUS)

# MultiQC ######################################################################

multiqc_config_yaml <- "multiqc_config.yaml"

run_MultiQC(
  tools = c(
    "fastqc",
    "trimmomatic",
    "sortmerna",
    "tophat",
    "hisat2",
    "samtools",
    "featureCounts"
  ),
  config = multiqc_config_yaml,
  force = FORCE_OVERWRITE
)

# Mapping stats ################################################################

if (!exists("number_trimmed")) {
  warning("Variable 'number_trimmed' undefined. Setting it to NA.")
  number_trimmed <- NA
}

if (!exists("number_nonrRNA")) {
  warning("Variable 'number_nonrRNA' undefined. Setting it to NA.")
  number_nonrRNA <- NA
}

# NOTE: Precise mapping stats requires BAM sorting.
# If TRUE, Bioconductor functions are used to determine mapping stats.

precise <- TRUE

mapping_stats_df <- calc_mapping_stats(
  bamFiles    = bam_files_sorted,
  fqFiles     = raw_fastq_files,
  anno        = anno_map,       # annotation used for mapping ~ can be different from counting!
  numReads    = number_raw,     # number of raw reads - mandatory!
  numTrimmed  = number_trimmed, # number of remaining reads after trimming - or NA
  numNonrRNA  = number_nonrRNA, # number of remaining reads after SortMeRNA - or NA
  libSizes    = lib_sizes,
  paired      = paired,
  precise     = precise,
  remove.na   = FALSE,
  cpus        = MAX_CPUS,
  featureType = "exon",
  samples     = NA
)

mapping_stats_df_corrected <- mapping_stats_df %>%
  remove_rownames() %>%
  mutate(ID = str_extract(string = fastq_files, pattern = "R\\d+")) %>%
  select(ID, everything()) %>%
  mutate(ID = make.unique(names = ID, sep = "_")) %>%
  mutate(ID = str_replace_all(ID, "_1", "_2")) %>%
  mutate(ID = str_replace_all(ID, "(R\\d+$)", "\\1_1")) %>%
  arrange(ID)

# Save mapping stats

data.table::fwrite(
  x = mapping_stats_df_corrected,
  file = paste0(tabDir, "/mapping_stats_df.tsv"),
  sep = "\t",
  quote = FALSE
)

WriteXLS::WriteXLS(
  x = mapping_stats_df_corrected,
  ExcelFileName = paste0(tabDir, "/mapping_stats_df.xlsx"),
  FreezeCol = 1
)

# Sanity check

if (paired) {
  if (F %in% (asPaired(number_raw) >= lib_sizes))
    stop("Error: number of reads mapping in exons should not exceed original number of reads!")
} else {
  if (F %in% (number_raw >= lib_sizes))
    stop("Error: number of reads mapping in exons should not exceed original number of reads!")
}

# Normalized count values ######################################################

rpkm <- get_rpkm(counts, gene_lengths, lib_sizes)
tpm  <- get_tpm(counts, gene_lengths, lib_sizes)
mrn  <- get_mrn(counts)

# Write values to CSV file

write_count_table(
  file   = file.path(tabDir, "rpkm"),
  counts = rpkm
)
write_count_table(
  file   = file.path(tabDir, "tpm"),
  counts = tpm
)
write_count_table(
  file   = file.path(tabDir, "mrn"),
  counts = mrn
)

# Write to XLS file

write_count_table(
  file       = file.path(tabDir, "rpkm"),
  counts     = rpkm,
  as.xls     = TRUE,
  sheetNames = basename(getwd())
)
write_count_table(
  file       = file.path(tabDir, "tpm"),
  counts     = tpm,
  as.xls     = TRUE,
  sheetNames = basename(getwd())
)
write_count_table(
  file       = file.path(tabDir, "mrn"),
  counts     = mrn,
  as.xls     = TRUE,
  sheetNames = basename(getwd())
)

# Clustering ###################################################################

# Load sample information and design matrix

sample_info <- readxl::read_excel(path = "path_to_metadata.xlsx")

design_matrix <- sample_info %>%
  mutate(pos_vs_neg = case_when(
    Sample %in% str_subset(Sample, "pos") ~ "treatment",
    Sample %in% str_subset(Sample, "neg") ~ "control"
  )) %>%
  column_to_rownames("ID") %>%
  select(-Sample) %>%
  as.matrix()

conds <- conditions_from_design(design_matrix)
conds

# Hierarchical clustering heat map version

make_heat_clustering_plot(
  file.path(plotDir, "heat_hierarchical_clustering"),
  counts    = counts[, conds != "none"],
  conds     = conds[conds != "none"],
  overwrite = TRUE
)

# Hierarchical clustering

make_hclust_plot(
  file.path(plotDir, "hierarchical_clustering"),
  counts    = counts[, conds != "none"],
  conds     = conds[conds != "none"],
  overwrite = TRUE
)

# Pearson correlation ##########################################################

make_correlation_plots(
  dat          = mrn,
  outDir       = plotDir,
  prefix       = "corr_",
  designMatrix = design_matrix,
  overwrite    = TRUE
)

# PCA plots ####################################################################

# NOTE: if designMatrix is supplied, 'conds' is ignored.
# In that case, use 'designMatrix = NA'

mrn_to_pca <- t(mrn)[, which(apply(t(mrn), 2, var) != 0)]
vst_to_pca <- DESeq2::varianceStabilizingTransformation(object = counts) %>% t()

pca_object <- prcomp(x = vst_to_pca, center = T)
pca_var <- round(x = summary(pca_object)$importance[2,] * 100, digits = 2)

pca_df <- pca_object$x %>%
  as.data.frame() %>%
  select(PC1, PC2) %>%
  rownames_to_column("ID") %>%
  left_join(x = ., y = sample_info, by = "ID") %>%
  mutate(Group = case_when(
    Sample %in% str_subset(string = Sample, pattern = "pos") ~ "Positive",
    Sample %in% str_subset(string = Sample, pattern = "neg") ~ "Negative"
  ))

pca_plot <- ggplot(data = pca_df, mapping = aes(x = PC1, y = PC2, fill = Group)) +
  geom_point(size = 8, shape = 21, color = "black") +
  labs(
    x = paste0("PC1 (", pca_var[1], "% variance explained)"),
    y = paste0("PC2 (", pca_var[2], "% variance explained)"),
    title = "Principal Component Analysis",
    subtitle = paste0(
      "(", pca_var[1] + pca_var[2], "% cumulative variance explained)"
    )
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, size = 1),
    text = element_text(size = 15),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_text(vjust = -1),
    plot.title = element_text(face = "bold"),
    legend.key = element_blank(),
    legend.title = element_text(face = "bold", hjust = 0.5)
  )

barplot_variance <- data.frame(Variance_Explained = pca_var) %>%
  rownames_to_column("PC") %>%
  ggplot(aes(x = PC, y = Variance_Explained)) +
  geom_col(color = "black", fill = "dodgerblue") +
  geom_label(
    aes(label = Variance_Explained),
    position = position_dodge(width = 0.9),
    size = 8
  ) +
  labs(x = "Principal Components", y = "Variance Explained") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_line(colour = "grey90"),
    panel.border = element_rect(fill = NA, size = 1),
    text = element_text(size = 15),
    axis.title = element_text(face = "bold"),
    axis.title.x = element_text(vjust = -1),
    plot.title = element_text(face = "bold"),
    legend.key = element_blank(),
    legend.title = element_text(face = "bold", hjust = 0.5)
  )

cairo_pdf(
  filename = "file_name.pdf"),
  width = 10,
  height = 10,
  onefile = TRUE,
  family = "sans",
  antialias = "subpixel"
)
pca_plot
barplot_variance
dev.off()

# NOTE: row names of design matrix has to be the same order as column names of
# counts dataset:
# colnames(counts) <- samples (only in case of having a SDRF file)

if (all(colnames(counts) != rownames(design_matrix))){
  counts_ordered <- counts[, match(rownames(design_matrix), colnames(counts))]
} else {
  counts_ordered <- counts
}

# Differential Expressed Genes #################################################

pvalcut <- 0.05
logfcCut <- 1
tools <- c("DESeq", "DESeq2", "edgeR", "limma")
deg_anno <- "" # Fill in if you need additional annotation in DEG tables (optional)

deg_res <- calculate_DEGs(
  counts        = counts_ordered,
  geneLengths   = gene_lengths,
  libSizes      = lib_sizes,
  designMatrix  = design_matrix,
  pValCut       = pvalcut,
  logfcCut      = logfcCut,
  tools         = tools,
  outDir        = degDir,
  prefix        = "",
  anno          = deg_anno,
  stop.on.error = FALSE,
  cpus          = 1,
  workers       = 10
)

# Intersection of DEGs

deg_intersect <- make_deg_overview_plot(
  outDir = degDir,
  degs = deg_res$DEGs,
  tools = tools
)

print(deg_intersect)

# Export resulting tables

data.table::fwrite(
  x = as.data.frame(deg_intersect),
  file = "deg_intersect.tsv",
  quote = FALSE,
  sep = "\t"
)

WriteXLS::WriteXLS(
  x = as.data.frame(deg_intersect),
  ExcelFileName = "deg_intersect.xlsx",
  row.names = TRUE
)

# Union of DEGs

degs_union <- deg_res$DEGs$pos_vs_neg %>%
  filter(DESeq == TRUE | DESeq2 == TRUE | Limma == TRUE | EdgeR == TRUE)

paste0("Union of DEGs (n): ", degs_union$id %>% length())

# Export resulting tables

data.table::fwrite(
  x = degs_union,
  file = "degs_union.tsv",
  quote = FALSE,
  sep = "\t"
)

WriteXLS::WriteXLS(
  x = degs_union,
  ExcelFileName = "degs_union.xlsx",
  row.names = TRUE
)
