if (!require("pacman"))
  install.packages("pacman")
pacman::p_load("ggpubr",
               "ggplot2",
               "tidyverse",
               "ggrepel",
               "glue",
               "Seurat"
)

# export folder
result_folder = "../results/GRN/Pando/"

# Pando GRN coefs (Pando input: variable TF genes)
trt_coefs =
  read_tsv("../results/GRN/Pando/Var_TFs_of_Jaspar_CISBP/Variable_TFs-trt_Pando_GRN_coefficients.tsv")
nt_coefs =
  read_tsv("../results/GRN/Pando/Var_TFs_of_Jaspar_CISBP/Variable_TFs-nt_Pando_GRN_coefficients.tsv")

trt_tf_mean_coefs = trt_coefs %>%
  dplyr::filter(padj < 0.05) %>%
  group_by(tf) %>%
  summarize(mean_coef = mean(estimate))
trt_tf_counts = trt_coefs %>% count(tf)

nt_tf_mean_coefs = nt_coefs %>%
  dplyr::filter(padj < 0.05) %>%
  group_by(tf) %>%
  summarize(mean_coef = mean(estimate))
nt_tf_counts = nt_coefs %>% count(tf)

# load and normalize scRNA-Seq data
counts = readRDS("../data/scRNA-Seq/Nerges.counts.filter.rds")
print(paste0("Total number of cells: ", length(unique(colnames(
  counts
)))))

meta = readRDS("../data/scRNA-Seq/Nerges.meta.filter.rds")
meta = meta %>% dplyr::filter(cluster_EML != "Undef")
meta = meta %>% dplyr::filter(cellType == "EZH2i_Naive_WT" |
                                cellType == "EZH2i_Naive_D7")
counts = counts[,meta$cell]

# read into Seurat
seurat_rna = CreateSeuratObject(
  counts = counts,
  project = "scRNA_EZH2i",
  min.cells = 3,
  min.features = 200
)
seurat_rna@meta.data = cbind(seurat_rna@meta.data, meta)

seurat_rna = NormalizeData(seurat_rna,
                           normalization.method = "LogNormalize",
                           scale.factor = 10000)
norm_count = seurat_rna@assays$RNA$data

## TF activity analysis - non-treated ELCs
nt_elcs = seurat_rna@meta.data$cell[which(
  seurat_rna@meta.data$cellType == "EZH2i_Naive_WT" &
    seurat_rna@meta.data$cluster_EML == "ELC"
)]
nt_elcs_norm_count = norm_count[, nt_elcs]
nt_elcs_tf_norm_count = nt_elcs_norm_count[nt_tf_mean_coefs$tf, ]
nt_elcs_tf_means = rowMeans(nt_elcs_tf_norm_count)
nt_elcs_tf_means = tibble(tf = names(nt_elcs_tf_means),
                          nt_ELC_mean_expr = nt_elcs_tf_means)

# TF activity
nt_tf_mean_coefs = nt_tf_mean_coefs %>% inner_join(., nt_elcs_tf_means, by = "tf") %>%
  mutate(tf_activity = mean_coef * nt_ELC_mean_expr)
nt_tf_mean_coefs = nt_tf_mean_coefs %>% inner_join(., nt_tf_counts, by = "tf")

tops = nt_tf_mean_coefs %>% arrange(desc(tf_activity)) %>%
  top_n(10, wt = tf_activity) %>% pull(tf)
bottoms = nt_tf_mean_coefs %>% arrange(desc(tf_activity)) %>%
  top_n(-10, wt = tf_activity) %>% pull(tf)

# non-treated TF activity rank plot
rank_nt_order = nt_tf_mean_coefs %>%
  arrange(desc(tf_activity)) %>%
  pull(tf) %>%
  unique
rank_nt_order = factor(nt_tf_mean_coefs$tf, levels = rank_nt_order)

nt_tf_mean_coefs = nt_tf_mean_coefs %>% mutate(bar_label = ifelse(tf %in% c(tops, bottoms),
                                                                  tf, ""))

rank_nt = ggplot(nt_tf_mean_coefs,
                 aes(x = rank_nt_order, y = tf_activity, fill = nt_ELC_mean_expr)) +
  geom_point(alpha = 1,
             aes(size = n),
             shape = 21,
             color = "black") +
  geom_text_repel(aes(label = bar_label),
                  min.segment.length = 0,
                  max.overlaps = 100) +
  ylim(-60, 60) +
  scale_fill_gradient2(
    low = "white",
    mid = "#f8f3b5",
    high = "red",
    midpoint = max(nt_tf_mean_coefs$nt_ELC_mean_expr) /
      2
  ) +
  labs(
    title = "non-treated TFs in ELCs",
    subtitle = "gene expressionally variable TFs",
    size = "# of interactions",
    fill = "norm. expr.",
    x = "predicted TF",
    y = "TF activity (mean coefficient x mean expression)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 9),
    plot.title = element_text(size = 10, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black")
  ) + geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey",
    linewidth = 0.5
  )
rank_nt

## TF activity analysis - 7D EZH2i treated ELCs
trt_elcs = seurat_rna@meta.data$cell[which(
  seurat_rna@meta.data$cellType == "EZH2i_Naive_D7" &
    seurat_rna@meta.data$cluster_EML == "ELC"
)]
trt_elcs_norm_count = norm_count[,trt_elcs]
trt_elcs_tf_norm_count = trt_elcs_norm_count[trt_tf_mean_coefs$tf, ]
trt_elcs_tf_means = rowMeans(trt_elcs_tf_norm_count)
trt_elcs_tf_means = tibble(tf = names(trt_elcs_tf_means),
                           trt_ELC_mean_expr = trt_elcs_tf_means)

# TF activity
trt_tf_mean_coefs = trt_tf_mean_coefs %>% inner_join(., trt_elcs_tf_means, by = "tf") %>%
  mutate(tf_activity = mean_coef * trt_ELC_mean_expr)
trt_tf_mean_coefs = trt_tf_mean_coefs %>% inner_join(., trt_tf_counts, by = "tf")

# trt - TF activity rank
tops = trt_tf_mean_coefs %>% arrange(desc(tf_activity)) %>%
  top_n(10, wt = tf_activity) %>% pull(tf)
bottoms = trt_tf_mean_coefs %>% arrange(desc(tf_activity)) %>%
  top_n(-10, wt = tf_activity) %>% pull(tf)

# EZH2i treated TF activity rank plot
rank_trt_order =  trt_tf_mean_coefs %>%
  arrange(desc(tf_activity)) %>%
  pull(tf) %>%
  unique
rank_trt_order = factor(trt_tf_mean_coefs$tf, levels = rank_trt_order)

trt_tf_mean_coefs = trt_tf_mean_coefs %>% mutate(bar_label = ifelse(tf %in% c(tops, bottoms),
                                                                    tf, ""))

rank_trt = ggplot(trt_tf_mean_coefs,
                  aes(x = rank_trt_order, y = tf_activity, fill = trt_ELC_mean_expr)) +
  geom_point(alpha = 1,
             aes(size = n),
             shape = 21,
             color = "black") +
  geom_text_repel(aes(label = bar_label),
                  min.segment.length = 0,
                  max.overlaps = 100) +
  ylim(-3, 3) +
  scale_fill_gradient2(
    low = "white",
    mid = "#f8f3b5",
    high = "red",
    midpoint = max(trt_tf_mean_coefs$trt_ELC_mean_expr) /
      2
  ) +
  labs(
    title = "EZH2i 7D TFs in ELCs",
    subtitle = "gene expressionally variable TFs",
    size = "# of interactions",
    fill = "norm. expr.",
    x = "predicted TF",
    y = "TF activity (mean coefficient x mean expression)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 9),
    plot.title = element_text(size = 10, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black")
  ) + geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey",
    linewidth = 0.5
  )
rank_trt

ranks = ggarrange(plotlist = list(rank_nt, rank_trt))
ggsave(
  glue("{result_folder}Pando-TF_activity_rank.png"),
  plot = ranks,
  width = 12,
  height = 6,
  dpi = 300,
)

ggsave(
  glue("{result_folder}Pando-TF_activity_rank.pdf"),
  plot = ranks,
  width = 14,
  height = 6,
  device = "pdf"
)

## finding differences in TF activities
# scatter plot ranked by the highest TF activity differences
trt_tf_mean_coefs_sc = trt_tf_mean_coefs %>%
  filter(tf %in% intersect(trt_tf_mean_coefs$tf, nt_tf_mean_coefs$tf)) %>%
  mutate(sample = "treated") %>%
  rename(mean_expr = trt_ELC_mean_expr)
nt_tf_mean_coefs_sc = nt_tf_mean_coefs %>%
  filter(tf %in% intersect(nt_tf_mean_coefs$tf, nt_tf_mean_coefs$tf)) %>%
  mutate(sample = "non-treated") %>%
  rename(mean_expr = nt_ELC_mean_expr)
sc_input =  bind_rows(nt_tf_mean_coefs_sc, trt_tf_mean_coefs_sc)

tops = sc_input %>% group_by(tf) %>% summarise(max_min = diff(range(tf_activity))) %>%
  arrange(desc(max_min)) %>%
  top_n(20, wt = max_min)

sc_input = sc_input %>% filter(tf %in% tops$tf)

order = sc_input %>% group_by(tf) %>% summarise(max_min = diff(range(tf_activity))) %>%
  arrange(desc(max_min)) %>% pull(tf) %>% unique
order = factor(sc_input$tf, levels = order)

sc = ggplot(sc_input, aes(x = order, y = tf_activity, fill = sample)) +
  geom_point(shape = 21, color = "black") +
  scale_fill_manual(values = c("grey", "gold")) +
  labs(title = "",
       x = "",
       y = "",
       fill = "") +
  labs(
    title = "",
    subtitle = "",
    size = "# of interactions",
    color = "treatment",
    x = "",
    y = "TF activity (mean coefficient x mean expression)"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 9),
    plot.title = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(
      size = 12,
      color = "black",
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 12, color = "black")
  )
sc

ggsave(
  glue("{result_folder}Pando-TF_activity_differences.png"),
  plot = sc,
  width = 8,
  height = 6,
  dpi = 300,
)

ggsave(
  glue("{result_folder}Pando-TF_activity_differences.pdf"),
  plot = sc,
  width = 8,
  height = 6,
  device = "pdf"
)
