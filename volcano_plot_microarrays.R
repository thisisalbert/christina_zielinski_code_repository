# Title: Volcano plot of microarray data
# Author: Albert García López

pacman::p_load(tidyverse, devtools, ggrepel)
outDir <- "path_to_dir"
stopifnot(dir.exists(outDir))

# Load processed data ##########################################################

df_processed <- read_csv("diff_expression_data.csv")

# Prepare data #################################################################

log_fc_cutoff <- 0.5
adj_pval_cutoff <- 0.05
sel_genes <- c("IL1A", "IL10", "GSDME", "GSDMD", "CASP1", "IL1B")

volc_data <- df_processed %>%
  select(genename, logFC, adj.P.Val) %>%
  mutate(DEG = case_when(
    abs(logFC) >= log_fc_cutoff & adj.P.Val < adj_pval_cutoff ~ TRUE,
    TRUE ~ FALSE
  )) %>%
  mutate(Status = case_when(
    DEG == TRUE & logFC >= log_fc_cutoff ~ "up",
    DEG == TRUE & logFC <= -log_fc_cutoff ~ "down",
    TRUE ~ "ns"
  )) %>%
  mutate(Labels = case_when(
    genename %in% sel_genes ~ genename,
    TRUE ~ NA_character_
  ))

# Plot data

font_size <- 16

volc_data_sel_genes <- volc_data %>%
  select(genename, logFC, `adj.P.Val`) %>%
  mutate(
    logFC = ifelse(genename %in% sel_genes, logFC, NA_real_),
    `adj.P.Val` = ifelse(genename %in% sel_genes, `adj.P.Val`, NA_real_)
  )

pdf(
  file = "file_name.pdf",
  width = 7,
  height = 6
)
volc_data %>%
  ggplot(aes(
    x = logFC,
    y = -log10(adj.P.Val),
    col = Status,
    label = Labels
  )) +
  geom_point(
    size = 2,
    alpha = 0.4
  ) +
  geom_point(aes(
    x = volc_data_sel_genes$logFC,
    y = -log10(volc_data_sel_genes$adj.P.Val)
  ),
  inherit.aes = FALSE,
  shape = 0,
  size = 5
  ) +
  geom_text_repel(
    size = 4,
    color = "black",
    fontface = "bold.italic",
    show.legend = FALSE,
    nudge_x = 0,
    nudge_y = 0.15,
  ) +
  geom_vline(
    xintercept = c(-log_fc_cutoff, log_fc_cutoff),
    col = "gray50",
    lty = "dashed",
    alpha = 0.5
  ) +
  geom_hline(
    yintercept = -log10(adj_pval_cutoff),
    col = "gray50",
    lty = "dashed",
    alpha = 0.5
  ) +
  scale_color_manual(
    breaks = c("up", "down", "ns"),
    values = c("#ce2627", "#4d69b1", "#989898")
  ) +
  scale_x_continuous(
    limits = c(-(ceiling(max(volc_data$logFC))), ceiling(max(volc_data$logFC))),
    n.breaks = 7
  ) +
  labs(
    x = expression(log[2]~fold~change),
    y = expression(-log[10]~adjusted~p.value)
  ) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme_bw() +
  theme(
    text = element_text(size = font_size),
    axis.text = element_text(size = font_size),
    axis.title = element_text(size = font_size),
    legend.text = element_text(size = font_size),
    legend.title = element_blank(),
    panel.border = element_rect(fill = NA, size = 1)
  )
dev.off()
