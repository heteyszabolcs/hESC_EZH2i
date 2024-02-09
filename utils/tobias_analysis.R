suppressPackageStartupMessages({
  library("ggplot2")
  library("ggrastr")
  library("ggrepel")
  library("data.table")
})

result_folder = "../results/motif_analysis/"

htselex_tobias = fread("../results/TOBIAS/BINDetect_HTSELEX_output/bindetect_results.txt")
htselex_tobias = htselex_tobias %>% dplyr::select(name, trt_nt_change, trt_nt_pvalue) %>% 
  separate(name, sep = "\\.", into = c("V1", "V2", "V3")) %>% 
  dplyr::select(gene_name = V3, trt_nt_change, trt_nt_pvalue)

volc_input = htselex_tobias %>% 
  mutate(group = case_when(
    trt_nt_change > 0.15 & trt_nt_pvalue < 0.05 ~ "up",
    trt_nt_change < -0.15 & trt_nt_pvalue < 0.05 ~ "down", .default = "unaltered"
  )) %>%
  mutate(sign_label = case_when(group == "up" ~ gene_name, 
                                group == "down" ~ gene_name, .default = NA))
labels = volc_input %>% pull(sign_label)

# add color, size, alpha indications
cols = c("up" = "#fc9272", "down" = "#3182bd", "unaltered" = "grey")
sizes = c("up" = 4, "down" = 4, "unaltered" = 2)
alphas = c("up" = 1, "down" = 1, "unaltered" = 0.5)

# EZH2i 7D ELC volcano
htselect_tobias_volc = volc_input %>%
  ggplot(aes(x = trt_nt_change,
             y = -log10(trt_nt_pvalue),
             fill = group,    
             size = group,
             alpha = group)) +
  ggrastr::geom_point_rast(shape = 21, # Specify shape and colour as fixed local parameters    
                           colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = 0.15,
             linetype = "dashed") +
  geom_vline(xintercept = -0.15,
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-0.5, 0.5, 0.1)),  	 
                     limits = c(-0.5, 0.5)) +
  scale_y_continuous(breaks = c(seq(0, 200, 50)),  	 
                     limits = c(0, 200)) +
  labs(
    title = "TOBIAS (HTSelex)",
    subtitle = "EZH2i 7D vs. non-treated ELCs",
    x = "fold enrichment",
    y = "-log10 p-value",
    fill = " "
  ) +
  guides(alpha = FALSE, size = FALSE, fill = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 6),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  geom_text_repel(label = labels, size = 4, max.overlaps = 110, min.segment.length = 0.1) # add labels
htselect_tobias_volc

ggsave(
  glue("{result_folder}ELC_TOBIAS_HTSelex.png"),
  plot = htselect_tobias_volc,
  width = 10,
  height = 7,
  dpi = 300,
)

ggsave(
  glue("{result_folder}ELC_TOBIAS_HTSelex.pdf"),
  plot = htselect_tobias_volc,
  width = 10,
  height = 7,
  device = "pdf"
)

# JASPAR 2020 output
jaspar_tobias = fread("../results/TOBIAS/BINDetect_Jaspar2020_output/bindetect_results.txt")
jaspar_tobias = jaspar_tobias %>% dplyr::select(name, trt_nt_change, trt_nt_pvalue) %>% 
  dplyr::select(gene_name = name, trt_nt_change, trt_nt_pvalue)

volc_input = jaspar_tobias %>% 
  mutate(group = case_when(
    trt_nt_change > 0.15 & trt_nt_pvalue < 0.05 ~ "up",
    trt_nt_change < -0.15 & trt_nt_pvalue < 0.05 ~ "down", .default = "unaltered"
  )) %>%
  mutate(sign_label = case_when(group == "up" ~ gene_name, 
                                group == "down" ~ gene_name, .default = NA))
labels = volc_input %>% pull(sign_label)

# add color, size, alpha indications
cols = c("up" = "#fc9272", "down" = "#3182bd", "unaltered" = "grey")
sizes = c("up" = 4, "down" = 4, "unaltered" = 2)
alphas = c("up" = 1, "down" = 1, "unaltered" = 0.5)

# EZH2i 7D ELC volcano
jaspar_tobias_volc = volc_input %>%
  ggplot(aes(x = trt_nt_change,
             y = -log10(trt_nt_pvalue),
             fill = group,    
             size = group,
             alpha = group)) +
  ggrastr::geom_point_rast(shape = 21, # Specify shape and colour as fixed local parameters    
                           colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = 0.15,
             linetype = "dashed") +
  geom_vline(xintercept = -0.15,
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-0.5, 0.5, 0.1)),  	 
                     limits = c(-0.5, 0.5)) +
  scale_y_continuous(breaks = c(seq(0, 200, 50)),  	 
                     limits = c(0, 200)) +
  labs(
    title = "TOBIAS (Jaspar)",
    subtitle = "EZH2i 7D vs. non-treated ELCs",
    x = "fold enrichment",
    y = "-log10 p-value",
    fill = " "
  ) +
  guides(alpha = FALSE, size = FALSE, fill = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 6),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  geom_text_repel(label = labels, size = 4, max.overlaps = 500, min.segment.length = 0.1) # add labels
jaspar_tobias_volc

ggsave(
  glue("{result_folder}ELC_TOBIAS_Jaspar.png"),
  plot = jaspar_tobias_volc,
  width = 15,
  height = 12,
  dpi = 300,
)

ggsave(
  glue("{result_folder}ELC_TOBIAS_Jaspar.pdf"),
  plot = jaspar_tobias_volc,
  width = 15,
  height = 12,
  device = "pdf"
)
