# Title: Barplots Fig 1SB
# Author: Albert García López

pacman::p_load(tidyverse, readxl, xlsx, janitor, skimr, ggpubr)
outDir <- "path_to_dir"
stopifnot(dir.exists(outDir))

# Load data ####################################################################

list_of_files <- list.files("path_to_file_dir", full.names = TRUE) %>%
  str_subset("\\.xlsx") %>%
  str_subset("Th17")

dfs_nest <- imap_dfr(
  list_of_files,
  ~ read_excel(.x) %>%
    clean_names() %>%
    mutate(DF_ID = .y, .before = everything()) %>%
    nest(data = -DF_ID)
) %>%
  mutate(DF_ID = case_when(
    DF_ID == 1 ~ "Th17_Day0_2h",
    DF_ID == 2 ~ "Th17_Day0_ut",
    DF_ID == 3 ~ "Th17_Day5_2h",
    DF_ID == 4 ~ "Th17_Day5_ut",
    DF_ID == 5 ~ "Th17-IL10+_Day5_c-MAF"
  ))

# Filter data and plot #########################################################

dfs_filtered <- dfs_nest %>%
  mutate(filtering = map(data, ~ filter(., gene_id %in% c("IL1A", "IL10")))) %>%
  select(-data) %>%
  unnest(filtering) %>%
  select(DF_ID, gene_id, matches("mean_fold_change")) %>%
  filter(DF_ID != "Th17-IL10+_Day5_c-MAF") %>%
  select(DF_ID, gene_id, mean_fold_change) %>%
  mutate(mean_fold_change = as.numeric(mean_fold_change)) %>%
  mutate(log2FC = log2(mean_fold_change) * -1) %>%
  mutate(day = str_extract(DF_ID, "Day\\d")) %>%
  mutate(gene_id = factor(gene_id, c("IL1A", "IL10")))

cairo_pdf(
  filename = "file_name.pdf",
  family = "Arial",
  onefile = TRUE,
  width = 6,
  height = 10
)
dfs_filtered %>%
  filter(day == "Day5") %>%
  ggbarplot(
    x = "gene_id",
    y = "log2FC",
    xlab = FALSE,
    add = "mean_se",
    error.plot = "errorbar"
  ) +
  stat_compare_means(
    method = "t.test",
    paired = FALSE,
    comparisons = list(c("IL1A", "IL10")),
    label = "p.signif",
    size = 10,
    bracket.size = 0.5
  ) +
  coord_fixed() +
  geom_hline(yintercept = 0, lty = 1, lwd = 0.5) +
  labs(y = "IL-10+ Th17 vs. IL-10- Th17\nlog2(fold change)", title = "T test") +
  theme(
    text = element_text(size = 18, family = "Arial"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
dfs_filtered %>%
  filter(day == "Day5") %>%
  ggbarplot(
    x = "gene_id",
    y = "log2FC",
    xlab = FALSE,
    add = "mean_se",
    error.plot = "errorbar"
  ) +
  stat_compare_means(
    method = "t.test",
    paired = FALSE,
    comparisons = list(c("IL1A", "IL10")),
    label = "p.format",
    size = 6,
    bracket.size = 0.5
  ) +
  coord_fixed() +
  geom_hline(yintercept = 0, lty = 1, lwd = 0.5) +
  labs(y = "IL-10+ Th17 vs. IL-10- Th17\nlog2(fold change)", title = "T test") +
  theme(
    text = element_text(size = 18, family = "Arial"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
dev.off()
