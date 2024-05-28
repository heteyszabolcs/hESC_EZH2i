suppressPackageStartupMessages({
  library("ggplot2")
  library("glue")
  library("ggrepel")
  library("tidyverse")
  library("data.table")
  library("ComplexHeatmap")
  library("circlize")
  library("viridis")
  library("ggpubr")
})

# result folder
result_folder = "../results/GRN/CellOracle/"

# network outputs (CellOracle)
trt_co_out = fread(glue("{result_folder}trt_network_score_table.tsv"))
nt_co_out = fread(glue("{result_folder}nt_network_score_table.tsv"))
tops = trt_co_out %>% group_by(cluster, V1) %>%
  summarise(max_eigenvector_centr = max(eigenvector_centrality)) %>%
  arrange(desc(max_eigenvector_centr)) %>%
  top_n(n = 5, wt = max_eigenvector_centr)

# edges
tlc_edges = fread(glue("{result_folder}trt_TLC-filtered_edges.tsv")) %>% mutate(type = "TLC")
melc_edges = fread(glue("{result_folder}trt_MeLC-filtered_edges.tsv")) %>% mutate(type = "MeLC")
hlc_edges = fread(glue("{result_folder}trt_HLC-filtered_edges.tsv")) %>% mutate(type = "HLC")
elc_edges = fread(glue("{result_folder}trt_ELC-filtered_edges.tsv")) %>% mutate(type = "ELC")
amlc_edges = fread(glue("{result_folder}trt_AMLC-filtered_edges.tsv")) %>% mutate(type = "AMLC")

edges = rbind(tlc_edges, melc_edges, hlc_edges, elc_edges, amlc_edges)

# heatmaps
create_heatmap = function(gene, cluster, save = TRUE) {
  print(gene)
  mat = edges %>% filter(source == gene) %>%
    pivot_wider(.,
                id_cols = type,
                names_from = target,
                values_from  = coef_abs) %>%
    column_to_rownames("type")
  mat[is.na(mat)] = 0
  
  if (length(mat) == 0) {
    return(print(glue("{gene} is not found")))
  }
  if (length(colnames(mat)) > 25) {
    mat = mat[, names(sort(colMeans(mat), decreasing = TRUE)[1:25])]
  }
  
  mat = as.matrix(mat)
  col_fun = colorRamp2(c(0, round(max(mat), 1) / 2, round(max(mat), 1) *
                           0.8), magma(3))
  hm = Heatmap(
    mat,
    name = "abs. coef",
    show_column_dend = FALSE,
    show_row_dend = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    col = col_fun,
    heatmap_width = unit(9, "cm"),
    heatmap_height = unit(10, "cm"),
    column_title = glue("{gene} target"),
    column_title_side = "bottom",
    column_names_gp = gpar(fontsize = 7)
  ) +
    rowAnnotation("abs coef distr" = anno_boxplot(
      mat,
      width = unit(4, "cm"),
      gp = gpar(fill = "#fc9272")
    ))
  
  if (save) {
    pdf(
      file = glue(
        "{result_folder}CellOracle_{cluster}_{gene}_abs_coefs.pdf"
      ),
      width = 6,
      height = 5
    )
    print(hm)
    dev.off()
  }
  
  print(hm)
  
}

#create_heatmap(gene = "HAND1")
lapply(tops %>% filter(cluster == "ELC") %>% pull(V1), create_heatmap, cluster = "ELC")
lapply(tops %>% filter(cluster == "TLC") %>% pull(V1), create_heatmap, cluster = "TLC")
lapply(tops %>% filter(cluster == "MeLC") %>% pull(V1), create_heatmap, cluster = "MeLC")

# scatterplots
elc_co_out = trt_co_out %>%
  filter(cluster == "ELC") %>%
  left_join(., nt_co_out %>% filter(cluster == "ELC"), by = "V1") %>%
  select(gene = V1,
         trt_degree_centrality_all = degree_centrality_all.x,
         nt_degree_centrality_all = degree_centrality_all.y) %>% 
  na.omit() %>% distinct_all(.keep_all = TRUE) %>% 
  mutate(fc = log2(trt_degree_centrality_all/nt_degree_centrality_all)) %>% 
  mutate(effect = case_when(fc > 0.9 ~ "up in EZH2i", 
                            fc < -0.9 ~ "up in non-trt",
                            .default = "unchanged")) %>% 
  mutate(effect = case_when(nt_degree_centrality_all < 0.05 &
                              trt_degree_centrality_all < 0.05 ~ "unchanged", 
                            .default = effect)) %>%
  mutate(label = ifelse(effect == "up in EZH2i" | 
                          effect == "up in non-trt", gene, ""))

labels = elc_co_out$label
elc_sc = ggplot(elc_co_out, aes(x = trt_degree_centrality_all, 
                       y = nt_degree_centrality_all,
                       color = effect)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#f0f0f0", "#fc9272", "#a6bddb")) +
  xlim(0, 0.3) +
  ylim(0., 0.3) +
  labs(
    title = "ELC",
    x = "EZH2i - degree centrality",
    y = "non-trt - degree centrality",
    color = "centrality"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")) +
  geom_text_repel(label = labels, size = 3, col = "black", 
                  max.overlaps = 30)
elc_sc   

tlc_co_out = trt_co_out %>%
  filter(cluster == "TLC") %>%
  left_join(., nt_co_out %>% filter(cluster == "TLC"), by = "V1") %>%
  select(gene = V1,
         trt_degree_centrality_all = degree_centrality_all.x,
         nt_degree_centrality_all = degree_centrality_all.y) %>% 
  na.omit() %>% distinct_all(.keep_all = TRUE) %>% 
  mutate(fc = log2(trt_degree_centrality_all/nt_degree_centrality_all)) %>% 
  mutate(effect = case_when(fc > 0.8 ~ "up in EZH2i", 
                            fc < -0.8 ~ "up in non-trt",
                            .default = "unchanged")) %>% 
  mutate(effect = case_when(nt_degree_centrality_all < 0.08 &
                              trt_degree_centrality_all < 0.08 ~ "unchanged", 
                            .default = effect)) %>%
  mutate(label = ifelse(effect == "up in EZH2i" | 
                          effect == "up in non-trt", gene, ""))

labels = tlc_co_out$label
tlc_sc = ggplot(tlc_co_out, aes(x = trt_degree_centrality_all, 
                       y = nt_degree_centrality_all,
                       color = effect)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#f0f0f0", "#fc9272", "#a6bddb")) +
  xlim(0, 0.3) +
  ylim(0., 0.3) +
  labs(
    title = "TLC",
    x = "EZH2i - degree centrality",
    y = "non-trt - degree centrality",
    color = "centrality"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")) +
  geom_text_repel(label = labels, size = 3, col = "black", 
                  max.overlaps = 30)
tlc_sc

melc_co_out = trt_co_out %>%
  filter(cluster == "MeLC") %>%
  left_join(., nt_co_out %>% filter(cluster == "MeLC"), by = "V1") %>%
  select(gene = V1,
         trt_degree_centrality_all = degree_centrality_all.x,
         nt_degree_centrality_all = degree_centrality_all.y) %>% 
  na.omit() %>% distinct_all(.keep_all = TRUE) %>% 
  mutate(fc = log2(trt_degree_centrality_all/nt_degree_centrality_all)) %>% 
  mutate(effect = case_when(fc > 0.8 ~ "up in EZH2i", 
                            fc < -0.8 ~ "up in non-trt",
                            .default = "unchanged")) %>% 
  mutate(effect = case_when(nt_degree_centrality_all < 0.08 &
                              trt_degree_centrality_all < 0.08 ~ "unchanged", 
                            .default = effect)) %>%
  mutate(label = ifelse(effect == "up in EZH2i" | 
                          effect == "up in non-trt", gene, ""))

labels = melc_co_out$label
melc_sc = ggplot(melc_co_out, aes(x = trt_degree_centrality_all, 
                       y = nt_degree_centrality_all,
                       color = effect)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#f0f0f0", "#fc9272", "#a6bddb")) +
  xlim(0, 0.3) +
  ylim(0., 0.3) +
  labs(
    title = "MeLC",
    x = "EZH2i - degree centrality",
    y = "non-trt - degree centrality",
    color = "centrality"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")) +
  geom_text_repel(label = labels, size = 3, col = "black", 
                  max.overlaps = 30)
melc_sc

scs = ggarrange(elc_sc, tlc_sc, melc_sc)
scs

ggsave(
  glue("{result_folder}CellOracle-centrality_changes.pdf"),
  plot = scs,
  width = 9,
  height = 7,
  device = "pdf"
)

elc_co_out = trt_co_out %>%
  filter(cluster == "ELC") %>%
  left_join(., nt_co_out %>% filter(cluster == "ELC"), by = "V1") %>%
  select(gene = V1,
         trt_eigenvector_centrality = eigenvector_centrality.x,
         nt_eigenvector_centrality = eigenvector_centrality.y) %>% 
  na.omit() %>% distinct_all(.keep_all = TRUE) %>% 
  mutate(fc = log2(trt_eigenvector_centrality/nt_eigenvector_centrality)) %>% 
  mutate(effect = case_when(fc > 1 ~ "up in EZH2i", 
                            fc < -1 ~ "up in non-trt",
                            .default = "unchanged")) %>% 
  mutate(effect = case_when(nt_eigenvector_centrality < 0.25 &
                              trt_eigenvector_centrality < 0.25 ~ "unchanged", 
                            .default = effect)) %>%
  mutate(label = ifelse(effect == "up in EZH2i" | 
                          effect == "up in non-trt", gene, ""))

labels = elc_co_out$label
elc_eigen_sc = ggplot(elc_co_out, aes(x = trt_eigenvector_centrality, 
                                y = nt_eigenvector_centrality,
                                color = effect)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#f0f0f0", "#fc9272", "#a6bddb")) +
  xlim(0, 1.2) +
  ylim(0., 1.2) +
  labs(
    title = "ELC",
    x = "EZH2i - eigenvector centrality",
    y = "non-trt - eigenvector centrality",
    color = "eigen centrality"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")) +
  geom_text_repel(label = labels, size = 3, col = "black", 
                  max.overlaps = 60)
elc_eigen_sc

tlc_co_out = trt_co_out %>%
  filter(cluster == "TLC") %>%
  left_join(., nt_co_out %>% filter(cluster == "TLC"), by = "V1") %>%
  select(gene = V1,
         trt_eigenvector_centrality = eigenvector_centrality.x,
         nt_eigenvector_centrality = eigenvector_centrality.y) %>% 
  na.omit() %>% distinct_all(.keep_all = TRUE) %>% 
  mutate(fc = log2(trt_eigenvector_centrality/nt_eigenvector_centrality)) %>% 
  mutate(effect = case_when(fc > 1 ~ "up in EZH2i", 
                            fc < -1 ~ "up in non-trt",
                            .default = "unchanged")) %>% 
  mutate(effect = case_when(nt_eigenvector_centrality < 0.25 &
                              trt_eigenvector_centrality < 0.25 ~ "unchanged", 
                            .default = effect)) %>%
  mutate(label = ifelse(effect == "up in EZH2i" | 
                          effect == "up in non-trt", gene, ""))

labels = tlc_co_out$label
tlc_eigen_sc = ggplot(tlc_co_out, aes(x = trt_eigenvector_centrality, 
                                      y = nt_eigenvector_centrality,
                                      color = effect)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#f0f0f0", "#fc9272", "#a6bddb")) +
  xlim(0, 1.2) +
  ylim(0., 1.2) +
  labs(
    title = "TLC",
    x = "EZH2i - eigenvector centrality",
    y = "non-trt - eigenvector centrality",
    color = "eigen centrality"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")) +
  geom_text_repel(label = labels, size = 3, col = "black", 
                  max.overlaps = 60)
tlc_eigen_sc

melc_co_out = trt_co_out %>%
  filter(cluster == "MeLC") %>%
  left_join(., nt_co_out %>% filter(cluster == "MeLC"), by = "V1") %>%
  select(gene = V1,
         trt_eigenvector_centrality = eigenvector_centrality.x,
         nt_eigenvector_centrality = eigenvector_centrality.y) %>% 
  na.omit() %>% distinct_all(.keep_all = TRUE) %>% 
  mutate(fc = log2(trt_eigenvector_centrality/nt_eigenvector_centrality)) %>% 
  mutate(effect = case_when(fc > 0.75 ~ "up in EZH2i", 
                            fc < -0.75 ~ "up in non-trt",
                            .default = "unchanged")) %>% 
  mutate(effect = case_when(nt_eigenvector_centrality < 0.25 &
                              trt_eigenvector_centrality < 0.25 ~ "unchanged", 
                            .default = effect)) %>%
  mutate(label = ifelse(effect == "up in EZH2i" | 
                          effect == "up in non-trt", gene, ""))

labels = melc_co_out$label
melc_eigen_sc = ggplot(melc_co_out, aes(x = trt_eigenvector_centrality, 
                                      y = nt_eigenvector_centrality,
                                      color = effect)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#f0f0f0", "#fc9272", "#a6bddb")) +
  xlim(0, 1.2) +
  ylim(0., 1.2) +
  labs(
    title = "MeLC",
    x = "EZH2i - eigenvector centrality",
    y = "non-trt - eigenvector centrality",
    color = "eigen centrality"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")) +
  geom_text_repel(label = labels, size = 3, col = "black", 
                  max.overlaps = 60)
melc_eigen_sc

eigen_scs = ggarrange(elc_eigen_sc, tlc_eigen_sc, melc_eigen_sc)
eigen_scs

ggsave(
  glue("{result_folder}CellOracle-eigen_centrality_changes.pdf"),
  plot = eigen_scs,
  width = 9,
  height = 7,
  device = "pdf"
)
