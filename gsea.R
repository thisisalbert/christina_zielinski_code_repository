# Title: GSEA of microarray data
# Author: Albert García López

pacman::p_load(tidyverse, janitor, vroom, readxl, skimr, clusterProfiler,
enrichplot, org.Hs.eg.db, DescTools, cowplot)
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("slice", "dplyr")

# Load results directory #######################################################

resDir <- "path_to_dir"

# Load processed data ##########################################################

df_processed <- vroom("diff_expression_file.csv")

# Exploratory Data Analysis ####################################################

# Data summary

skim(df_processed)

# Any duplicates?

df_processed %>%
  sapply(anyDuplicated) %>%
  enframe() %>%
  filter(value > 0)

# GSEA #########################################################################

# Prepare gene_list_object:

# 1. Numeric vector: fold change or other type of numerical variable.
# 2. Named vector: every number has a name, the corresponding gene ID.
# 3. Sorted vector: number should be sorted in decreasing order.

gene_list <- df_processed %>%
  select(entrez, logFC) %>%
  mutate(logFC = -logFC) %>%
  deframe() %>%
  sort(decreasing = TRUE)

# Gene Ontology (GO)

set.seed(2)
gsea_go <- gseGO(
  geneList = gene_list,
  ont = "ALL",
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05,
  eps = 0,
  verbose = FALSE
) %>%
  setReadable(OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")

# Ranked plots #################################################################

# Keywords

keywords <- c("proliferation", "death", "activation", "survival", "viable", "viability")

# Select pathways to plot

sel_paths <- gsea_go@result %>%
  as_tibble() %>%
  filter(str_detect(Description, paste(keywords, collapse = "|")))

sel_ids <- which(gsea_go@result$ID %in% sel_paths$ID)

# Plots

ranked_plots <- list()

for (i in 1:length(sel_ids)) {

  metrics <- gsea_go@result %>%
    as_tibble() %>%
    slice(sel_ids[i]) %>%
    select(setSize, NES, pvalue, p.adjust, qvalues)

  ranked_plots[[i]] <- gseaplot2(
    x = gsea_go,
    geneSetID = sel_ids[i],
    title = StrCap(gsea_go@result[sel_ids[i], "Description"]),
    color = "mediumseagreen",
    base_size = 25
  ) +
    annotate(
      geom = "text",
      x = 10000,
      y = -1.5,
      size = 7,
      label = paste0(
        "Gene set size = ", metrics$setSize, "\n",
        "NES = ", round(metrics$NES, 3), " / ",
        "p-value = ", formatC(x = metrics$pvalue, digits = 3, format = "e"), " / ",
        "adj. p-value = ", formatC(x = metrics$p.adjust, digits = 3, format = "e")
      )
    ) +
    annotate(
      geom = "text",
      x = 1500,
      y = 2,
      size = 9,
      fontface = 2,
      color = "#CF2E28",
      label = "IL1B"
    ) +
    annotate(
      geom = "text",
      x = 19000,
      y = 2,
      size = 9,
      fontface = 2,
      color = "#2062A6",
      label = "Th17",
    ) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.margin = unit(c(0,1,0,1), "cm")
    )

}

# Visualize and export plots ###################################################

pdf(
  file = "file_name.pdf",
  width = 35,
  height = 20,
  family = "sans",
  compress = TRUE,
  onefile = TRUE
)
plot_grid(plotlist = ranked_plots[1:6], align = "hv", axis = "tblr")
plot_grid(plotlist = ranked_plots[7:12], align = "hv", axis = "tblr")
plot_grid(plotlist = ranked_plots[13:18], align = "hv", axis = "tblr")
plot_grid(plotlist = ranked_plots[19:24], align = "hv", axis = "tblr")
plot_grid(plotlist = ranked_plots[25], align = "hv", axis = "tblr")
dev.off()

# Export metrics for the selected pathways #####################################

vroom_write(
  x = sel_paths,
  file = "file_name.tsv"
)
