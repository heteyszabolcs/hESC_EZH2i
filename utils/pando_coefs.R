library("ggpubr")
library("ggplot2")
library("tidyverse")
library("ggpubr")
library("ggrepel")
library("glue")

# export folder
result_folder = "../results/GRN/Pando/"

# Pando GRN coefs
trt_coefs = 
  read_tsv("../results/GRN/Pando/scRNASeq_trt_elc_vs_nt_elc_sign_ups_cands/trt_Pando_GRN_coefficients.tsv")
nt_coefs =
  read_tsv("../results/GRN/Pando/scRNASeq_trt_elc_vs_nt_elc_sign_ups_cands/nt_Pando_GRN_coefficients.tsv")

trt_coefs = trt_coefs %>% mutate(sample = "EZH2i 7D naive ELCs")
int = intersect(trt_coefs %>% dplyr::filter(padj < 0.05) %>% pull(tf), 
          nt_coefs %>% dplyr::filter(padj < 0.05) %>% pull(tf))
trt_tfs = setdiff(trt_coefs %>% dplyr::filter(padj < 0.05) %>% pull(tf), 
          nt_coefs %>% dplyr::filter(padj < 0.05) %>% pull(tf))
nt_tfs = setdiff(nt_coefs %>% dplyr::filter(padj < 0.05) %>% pull(tf), 
        trt_coefs %>% dplyr::filter(padj < 0.05) %>% pull(tf))

# treated-only rank
trt_rank_plot = trt_coefs %>% 
  mutate(label = paste0(tf, "-", target, "-", region)) %>% 
  dplyr::filter(tf %in% trt_tfs & padj < 0.05) 

tops = trt_rank_plot %>% arrange(desc(estimate))  %>% top_n(10, wt = estimate) %>% pull(tf)
bottoms = trt_rank_plot %>% arrange(desc(estimate))  %>% top_n(-10, wt = estimate) %>% pull(tf)
rank_trt_order =  trt_rank_plot %>% arrange(desc(estimate)) %>% pull(label) %>% unique
rank_trt_order = factor(trt_rank_plot$label, levels = rank_trt_order)

trt_rank_plot = trt_rank_plot %>% mutate(bar_label = ifelse(tf %in% c(tops, bottoms) &
                                                            estimate > 10 | estimate < -10, 
                                                          tf, ""))

rank_trt = ggplot(trt_rank_plot, aes(x = rank_trt_order, y = estimate)) +
  geom_bar(stat="identity", position=position_dodge(), color = "#fc9272", 
           fill = "#fc9272") +
  geom_text_repel(aes(label = bar_label), min.segment.length = 0,
                  max.overlaps = 100, direction = "y",
                  nudge_x = 10, nudge_y = 5) +
  ylim(-100, 100) +
  labs(
    title = "EZH2i 7D only TFs in ELCs",
    subtitle = "",
    x = "",
    y = "PANDO coefficient"
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
    axis.title.y = element_text(size = 20, color = "black")
  ) 

nt_rank_plot = nt_coefs %>% 
  mutate(label = paste0(tf, "-", target, "-", region)) %>% 
  dplyr::filter(tf %in% nt_tfs & padj < 0.05) 

tops = nt_rank_plot %>% arrange(desc(estimate))  %>% top_n(10, wt = estimate) %>% pull(tf)
bottoms = nt_rank_plot %>% arrange(desc(estimate))  %>% top_n(-10, wt = estimate) %>% pull(tf)
rank_nt_order =  nt_rank_plot %>% arrange(desc(estimate)) %>% pull(label) %>% unique
rank_nt_order = factor(nt_rank_plot$label, levels = rank_nt_order)

nt_rank_plot = nt_rank_plot %>% mutate(bar_label = ifelse(tf %in% c(tops, bottoms) &
                                                              estimate > 2 | estimate < -2, 
                                                            tf, ""))

rank_nt = ggplot(nt_rank_plot, aes(x = rank_nt_order, y = estimate)) +
  geom_bar(stat="identity", position=position_dodge(), color = "black", 
           fill = "#3182bd") +
  geom_text_repel(aes(label = bar_label), min.segment.length = 0,
                    max.overlaps = 40, direction = "y",
                  nudge_x = 10, nudge_y = 5) +
  ylim(-50, 50) +
  labs(
    title = "non-treated only TFs in ELCs",
    subtitle = "",
    x = "",
    y = "PANDO coefficient"
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
    axis.title.y = element_text(size = 20, color = "black")
  ) 

int_rank_plot = nt_coefs %>% 
  mutate(label = paste0(tf, "-", target, "-", region)) %>% 
  dplyr::filter(tf %in% int & padj < 0.05) 

tops = int_rank_plot %>% arrange(desc(estimate))  %>% top_n(10, wt = estimate) %>% pull(tf)
bottoms = int_rank_plot %>% arrange(desc(estimate))  %>% top_n(-10, wt = estimate) %>% pull(tf)
int_rank_order =  int_rank_plot %>% arrange(desc(estimate)) %>% pull(label) %>% unique
int_rank_order = factor(int_rank_plot$label, levels = int_rank_order)

int_rank_plot = int_rank_plot %>% mutate(bar_label = ifelse(tf %in% c(tops, bottoms) &
                                                            estimate > 2 | estimate < -2, 
                                                          tf, ""))

rank_int_nt = ggplot(int_rank_plot, aes(x = int_rank_order, y = estimate)) +
  geom_bar(stat="identity", position=position_dodge(), color = "black", 
           fill = "#2ca25f") +
  geom_text_repel(aes(label = bar_label), min.segment.length = 0,
                  max.overlaps = 40, direction = "y",
                  nudge_x = 10, nudge_y = 5) +
  ylim(-50, 50) +
  labs(
    title = "Overlapping TFs in EZH2i 7D ELCs",
    subtitle = "",
    x = "",
    y = "PANDO coefficient"
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
    axis.title.y = element_text(size = 20, color = "black")
  ) 

int_rank_plot2= trt_coefs %>% 
  mutate(label = paste0(tf, "-", target, "-", region)) %>% 
  dplyr::filter(tf %in% int & padj < 0.05) 

tops = int_rank_plot2 %>% arrange(desc(estimate))  %>% top_n(10, wt = estimate) %>% pull(tf)
bottoms = int_rank_plot2 %>% arrange(desc(estimate))  %>% top_n(-10, wt = estimate) %>% pull(tf)
int_rank_order2 =  int_rank_plot2 %>% arrange(desc(estimate)) %>% pull(label) %>% unique
int_rank_order2 = factor(int_rank_plot2$label, levels = int_rank_order2)

int_rank_plot2 = int_rank_plot2 %>% mutate(bar_label = ifelse(tf %in% c(tops, bottoms) &
                                                              estimate > 2 | estimate < -2, 
                                                            tf, ""))

rank_int_trt = ggplot(int_rank_plot2, aes(x = int_rank_order2, y = estimate)) +
  geom_bar(stat="identity", position=position_dodge(), color = "#2ca25f", 
           fill = "#2ca25f") +
  geom_text_repel(aes(label = bar_label), min.segment.length = 0,
                  max.overlaps = 100, direction = "y",
                  nudge_x = 10, nudge_y = 5) +
  ylim(-200, 200) +
  labs(
    title = "Overlapping TFs in non-treated ELCs",
    subtitle = "",
    x = "",
    y = "PANDO coefficient"
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
    axis.title.y = element_text(size = 20, color = "black")
  ) 

ranks = ggarrange(rank_nt, rank_int_nt, rank_int_trt, rank_trt)
ggsave(
  glue("{result_folder}trtELC_vs_ntELC_upreg_genes_input-rank_plot.pdf"),
  plot = ranks,
  width = 12,
  height = 10,
  device = "pdf"
)
ggsave(
  glue("{result_folder}trtELC_vs_ntELC_upreg_genes_input-rank_plot.png"),
  plot = ranks,
  width = 12,
  height = 10,
  dpi = 500
)


nt_coefs = nt_coefs %>% mutate(sample = "non-treated naive ELCs")
coefs = rbind(trt_coefs, nt_coefs) +

trt_bp = ggplot(coefs, aes(x = sample, y = estimate, fill = sample)) +
  geom_boxplot(alpha = 1, aes(fill = sample)) +
  #ylim(-500, 500) +
  labs(title = "",
       x = "",
       y = "",
       fill = "") +
  guides(fill = "none") +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 13, color = "black", angle = 45, vjust = 1, hjust = 1))
trt_bp

top_pos_estimates = trt_coefs %>% group_by(tf) %>% summarise(mean = mean(estimate)) %>% 
  arrange(desc(mean)) %>% top_n(20, wt = mean) %>% pull(tf)
top_neg_estimates = trt_coefs %>% group_by(tf) %>% summarise(mean = mean(estimate)) %>% 
  arrange(mean) %>% top_n(-20, wt = mean) %>% pull(tf)

pos_hm_input = trt_coefs %>% dplyr::filter(tf %in% top_pos_estimates) %>% 
  dplyr::filter(padj < 0.05)

pos_hm = ggplot(pos_hm_input, aes(x = target, y = tf, fill = estimate)) +
  geom_tile(color = "black",
            lwd = 0.2,
            linetype = 1) +
  scale_fill_gradient2(
    low = "#fee0d2",
    mid = "#fc9272",
    high = "darkred",
    midpoint = max(ceiling(pos_hm_input$estimate))/2,
    limits = c(min(ceiling(pos_hm_input$estimate)-20), max(ceiling(pos_hm_input$estimate))+1)
  ) +
  xlab(label = "PANDO input") +
  ylab(label = "predicted TF") +
  labs(fill = "coefficient", title = "") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.6),
    axis.text.x = element_text(
      color = "black",
      size = 8,
      angle = 90,
      hjust = 1,
      vjust = 0.6
    ),
    axis.text.y = element_text(color = "black", size = 5),
    axis.title = element_text(size = 14)
  ) +
  coord_fixed()

neg_hm_input = trt_coefs %>% dplyr::filter(tf %in% top_neg_estimates) %>% 
  dplyr::filter(padj < 0.05)

neg_hm = ggplot(neg_hm_input, aes(x = target, y = tf, fill = estimate)) +
  geom_tile(color = "black",
            lwd = 0.2,
            linetype = 1) +
  scale_fill_gradient2(
    low = "darkblue",
    mid = "#3182bd",
    high = "#deebf7",
    midpoint = -200,
    limits = c(-400, 25)
  ) +
  xlab(label = "PANDO input") +
  ylab(label = "predicted TF") +
  labs(fill = "coefficient", title = "") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.6),
    axis.text.x = element_text(
      color = "black",
      size = 8,
      angle = 90,
      hjust = 1,
      vjust = 0.6
    ),
    axis.text.y = element_text(color = "black", size = 5),
    axis.title = element_text(size = 14)
  ) +
  coord_fixed()

ggarrange(pos_hm, neg_hm)

# CDX1 in trt network
cdx1_trt = trt_coefs %>% dplyr::filter(target == "CDX1" & pval < 0.05)
cdx1_nt = nt_coefs %>% dplyr::filter(target == "CDX1" & pval < 0.05)

order = cdx1 %>% arrange(desc(estimate)) %>% pull(tf)
order = factor(cdx1$tf, levels = order)
cdx1_bar = ggplot(cdx1, aes(x = order, y = estimate)) +
  geom_bar(stat="identity", position=position_dodge(), color = "black", fill = "lightblue") +
  #ylim(0, 0.20) +
  labs(
    title = "CDX1 downstream interactors (p-value < 0.05)",
    subtitle = "EZH2i 7D naive ELCs",
    x = "",
    y = "PANDO coefficient"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 9),
    plot.title = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 15, color = "black", angle = 0),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) 
cdx1_bar

ggsave(
  glue("{result_folder}Pando-TOBIAS_ChromVar_input-trt_coefs-CDX1_coefs.pdf"),
  plot = cdx1_bar,
  width = 8,
  height = 6,
  device = "pdf"
)

# TF specific networks
trt_grn = readRDS("../results/GRN/Pando/chromVAR_TOBIAS_candidates/trt_Pando_GRN.Rds")
nt_grn = readRDS("../results/GRN/Pando/chromVAR_TOBIAS_candidates/nt_Pando_GRN.Rds")


GetNetwork(trt_grn)
GetNetwork(nt_grn)
trt_coefs = coef(trt_grn)
nt_coefs = coef(nt_grn)

trt_modules = find_modules(
  trt_grn, 
  p_thresh = 1,
  nvar_thresh = 0, 
  min_genes_per_module = 0, 
  rsq_thresh = 0.1
)
nt_modules = find_modules(
  nt_grn, 
  p_thresh = 1,
  nvar_thresh = 0, 
  min_genes_per_module = 0, 
  rsq_thresh = 0.1
)

trt_modules_outputs = NetworkModules(trt_modules)
nt_modules_outputs = NetworkModules(nt_modules) 
trt_modules_meta = trt_modules_outputs@meta
nt_modules_meta = nt_modules_outputs@meta
trt_network_vis = get_network_graph(trt_modules, umap_method = "none")
nt_network_vis = get_network_graph(nt_modules, umap_method = "none")
get_network_graph(trt_modules, umap_method = "none")
get_network_graph(nt_modules, umap_method = "none")
plot_network_graph(trt_network_vis, layout='fr', text_size = 14)
plot_network_graph(nt_network_vis, layout='fr', text_size = 14)


trt_rsquares = trt_network_vis@grn@networks$glm_network@fit
require(tidygraph)
trt_centrality = trt_network_vis@grn@networks$glm_network@graphs$module_graph %>%
  activate(nodes) %>% 
  as_tibble() %>% 
  mutate(sample = "EZH2i 7D naive ELCs")
nt_centrality = nt_network_vis@grn@networks$glm_network@graphs$module_graph %>%
  activate(nodes) %>% 
  as_tibble() %>% 
  mutate(sample = "non-treated naive ELCs")
centralities = rbind(nt_centrality, trt_centrality)

order = centralities %>% group_by(name) %>% summarize(mean = mean(centrality)) %>% 
  arrange(desc(mean)) %>% pull(name)
order = factor(centralities$name, levels = order)
order2 = factor(centralities$sample, levels = c("non-treated naive ELCs", "EZH2i 7D naive ELCs"))

# centrality plot
ggplot(centralities, aes(x = order, y = centrality, fill = order2)) +
  geom_bar(stat="identity", position=position_dodge(), color = "black") +
  scale_fill_manual(values = c("#bdbdbd", "#fec44f")) +
  #ylim(0, 0.20) +
  labs(
    title = "GRN centralities",
    x = "candidate (ChromVar / TOBIAS outputs)",
    y = "centrality (PageRank)",
    fill = "GRN"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 9),
    plot.title = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 7, color = "black", angle = 45, vjust = 1, hjust = 0.9),
    axis.text.y = element_text(size = 7, color = "black")
  ) 

ggsave(
  glue("{result_folder}GRN_centrality_bars-ChromVar_TOBIAS_GRN.pdf"),
  plot = last_plot(),
  width = 7,
  height = 4
)


## networks
# trt
pdf(
  file = glue("{result_folder}chromVAR_TOBIAS_graph-TFCP2_network_trt.pdf"),
  width = 6,
  height = 4
)
trt_network_vis = get_tf_network(trt_network_vis, tf='TFCP2')
plot_tf_network(trt_network_vis, tf='TFCP2', circular=F)
dev.off()

pdf(
  file = glue("{result_folder}chromVAR_TOBIAS_graph-CDX1_network_trt.pdf"),
  width = 6,
  height = 4
)
trt_network_vis = get_tf_network(trt_network_vis, tf='CDX1')
plot_tf_network(trt_network_vis, tf='CDX1', circular=F)
dev.off()

# nt
pdf(
  file = glue("{result_folder}chromVAR_TOBIAS_graph-TFCP2_network_nt.pdf"),
  width = 6,
  height = 4
)
nt_network_vis = get_tf_network(nt_network_vis, tf='TFCP2')
plot_tf_network(nt_network_vis, tf='TFCP2', circular=F)
dev.off()

pdf(
  file = glue("{result_folder}chromVAR_TOBIAS_graph-CDX1_network_nt.pdf"),
  width = 6,
  height = 4
)
nt_network_vis = get_tf_network(nt_network_vis, tf='CDX1')
plot_tf_network(nt_network_vis, tf='CDX1', circular=F)
dev.off()





