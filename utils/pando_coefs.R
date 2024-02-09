library("ggpubr")
library("ggplot2")
library("tidyverse")
library("ggpubr")
library("glue")

# export folder
result_folder = "../results/GRN/Pando/"

# Pando GRN coefs
trt_coefs = 
  read_tsv("../results/GRN/Pando/chromVAR_TOBIAS_candidates/trt_Pando_GRN_coefficients.tsv")

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

ggsave(
  glue("{result_folder}Pando-TOBIAS_ChromVar_input-trt_coefs-sign_preds_hm.pdf"),
  plot = last_plot(),
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

