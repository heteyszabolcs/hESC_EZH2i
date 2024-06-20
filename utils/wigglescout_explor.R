if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  "tidyverse",
  "data.table",
  "wigglescout",
  "ggpubr",
  "glue",
  "RColorBrewer",
  "ComplexHeatmap",
  "viridis",
  "circlize"
)

# export folder
result_folder = "../results/Seurat_integration/scATAC_annotation/"

# input data
# bulk ATAC signals
bigwigs = list.files("../data/bulk_ATAC-Seq/KO_experiments/bigwig/", full.names = TRUE)
pgl_bigwigs = list.files("../data/bulk_ATAC-Seq/IWP2_ChIR_EZH2i/bigwig/", full.names = TRUE)
yapi_bigwigs = list.files("../data/bulk_ATAC-Seq/YAPi/bigwig/mergedLibrary/", full.names = TRUE)
# scATAC peaks
tlc_peaks = "../results/Seurat_integration/scATAC_annotation/cluster_peaks/TLC_spec_peaks_up-log2fc_2_adjp0.05_Chens_annotation.bed"
melc_peaks = "../results/Seurat_integration/scATAC_annotation/cluster_peaks/MeLC_spec_peaks_up-log2fc_2_adjp0.05_Chens_annotation.bed"
aelc_peaks = "../results/Seurat_integration/scATAC_annotation/cluster_peaks/aELC_spec_peaks_up-log2fc_1_adjp0.05_Chens_annotation.bed"

# aggregation function
aggr_signal = function(bigwigs, bed) {
  output = bw_loci(bigwigs, loci = bed)
  output = as.data.frame(output)
  return(output)
}

## visualizations of peak signals of EZH2i scATAC-Seq
tlc = aggr_signal(bigwigs = bigwigs, bed = tlc_peaks)
tlc = tlc %>% mutate(peak = paste(seqnames, as.character(start), as.character(end), sep = "-")) %>%
  dplyr::select(-width, -strand, -seqnames, -start, -end)
rownames(tlc) = tlc$peak
tlc_mat = tlc %>% dplyr::select(-peak) %>% as.matrix

melc = aggr_signal(bigwigs = bigwigs, bed = melc_peaks)
melc = melc %>% mutate(peak = paste(seqnames, as.character(start), as.character(end), sep = "-")) %>%
  dplyr::select(-width, -strand, -seqnames, -start, -end)
rownames(melc) = melc$peak
melc_mat = melc %>% dplyr::select(-peak) %>% as.matrix

aelc = aggr_signal(bigwigs = bigwigs, bed = aelc_peaks)
aelc = aelc %>% mutate(peak = paste(seqnames, as.character(start), as.character(end), sep = "-")) %>%
  dplyr::select(-width, -strand, -seqnames, -start, -end)
rownames(aelc) = aelc$peak
aelc_mat = aelc %>% dplyr::select(-peak) %>% as.matrix

melc_pgl = aggr_signal(bigwigs = pgl_bigwigs, bed = melc_peaks)
melc_pgl = melc_pgl %>% mutate(peak = paste(seqnames, as.character(start), as.character(end), sep = "-")) %>%
  dplyr::select(-width, -strand, -seqnames, -start, -end)

tlc_pgl = aggr_signal(bigwigs = pgl_bigwigs, bed = tlc_peaks)
tlc_pgl = tlc_pgl %>% mutate(peak = paste(seqnames, as.character(start), as.character(end), sep = "-")) %>%
  dplyr::select(-width, -strand, -seqnames, -start, -end)

melc_yapi = aggr_signal(bigwigs = yapi_bigwigs, bed = melc_peaks)
melc_yapi = melc_yapi %>% mutate(peak = paste(seqnames, as.character(start), as.character(end), sep = "-")) %>%
  dplyr::select(-width, -strand, -seqnames, -start, -end)

tlc_yapi = aggr_signal(bigwigs = yapi_bigwigs, bed = tlc_peaks)
tlc_yapi = tlc_yapi %>% mutate(peak = paste(seqnames, as.character(start), as.character(end), sep = "-")) %>%
  dplyr::select(-width, -strand, -seqnames, -start, -end)

# scatterplots
# bulk signals over scATAC peaks
sc_melc = melc %>% mutate(cell_type = "MeLC") %>% 
  dplyr::select(cell_type, everything()) %>% 
  pivot_longer(., names_to = "ko", values_to = "agg_sign", 
               cols = starts_with("ATAC")) %>% 
  dplyr::select(-peak) %>% 
  group_by(cell_type, ko) %>% summarise(mean = mean(agg_sign), sd = sd(agg_sign))
sc_tlc = tlc %>% mutate(cell_type = "TLC") %>% 
  dplyr::select(cell_type, everything()) %>% 
  pivot_longer(., names_to = "ko", values_to = "agg_sign", 
               cols = starts_with("ATAC")) %>% 
  dplyr::select(-peak) %>% 
  group_by(cell_type, ko) %>% summarise(mean = mean(agg_sign), sd = sd(agg_sign))
sc_aelc = aelc %>% mutate(cell_type = "aELC") %>% 
  dplyr::select(cell_type, everything()) %>% 
  pivot_longer(., names_to = "ko", values_to = "agg_sign", 
               cols = starts_with("ATAC")) %>% 
  dplyr::select(-peak) %>% 
  group_by(cell_type, ko) %>% summarise(mean = mean(agg_sign), sd = sd(agg_sign))
sc_melc_pgl = melc_pgl %>% mutate(cell_type = "MeLC") %>% 
  dplyr::select(cell_type, everything()) %>% 
  pivot_longer(., names_to = "ko", values_to = "agg_sign", 
               cols = starts_with("ATAC")) %>% 
  dplyr::select(-peak) %>% 
  group_by(cell_type, ko) %>% summarise(mean = mean(agg_sign), sd = sd(agg_sign))
sc_tlc_pgl = tlc_pgl %>% mutate(cell_type = "TLC") %>% 
  dplyr::select(cell_type, everything()) %>% 
  pivot_longer(., names_to = "ko", values_to = "agg_sign", 
               cols = starts_with("ATAC")) %>% 
  dplyr::select(-peak) %>% 
  group_by(cell_type, ko) %>% summarise(mean = mean(agg_sign), sd = sd(agg_sign))
sc_melc_yapi = melc_yapi %>% mutate(cell_type = "MeLC") %>% 
  dplyr::select(cell_type, everything()) %>% 
  pivot_longer(., names_to = "ko", values_to = "agg_sign", 
               cols = starts_with("ATAC")) %>% 
  dplyr::select(-peak) %>% 
  group_by(cell_type, ko) %>% summarise(mean = mean(agg_sign), sd = sd(agg_sign))
sc_tlc_yapi = tlc_yapi %>% mutate(cell_type = "TLC") %>% 
  dplyr::select(cell_type, everything()) %>% 
  pivot_longer(., names_to = "ko", values_to = "agg_sign", 
               cols = starts_with("ATAC")) %>% 
  dplyr::select(-peak) %>% 
  group_by(cell_type, ko) %>% summarise(mean = mean(agg_sign), sd = sd(agg_sign))

sc_input = rbind(sc_tlc, sc_melc) %>% pivot_wider(., names_from = cell_type, values_from = mean,
                                                  id_cols = ko) %>% 
  separate(ko, sep = "ATAC_", into = c("rest", "ko")) %>% 
  separate(ko, sep = ".mRp.clN.uniq", into = c("ko", "rest")) %>% 
  dplyr::select(-rest)

kos = ggplot(sc_input, aes(x = TLC, y = MeLC, color = ko)) +
  geom_point(size = 7) +
  scale_color_brewer(palette = "Set3") +
  labs(
    title = "mean bulk ATAC signals (KOs) over peak sets.",
    x = "TLC",
    y = "MeLC",
    color = " "
  ) +
  xlim(0,7) +
  ylim(0,4) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"))
ggsave(
  glue("{result_folder}mean_bulkATAC_over_peaksets.png"),
  plot = kos,
  width = 6,
  height = 6,
  dpi = 500,
)
ggsave(
  glue("{result_folder}mean_bulkATAC_over_peaksets.pdf"),
  plot = kos,
  width = 6,
  height = 6
)

sc_input_pgl = rbind(sc_tlc_pgl, sc_melc_pgl) %>%
  pivot_wider(.,
              names_from = cell_type,
              values_from = mean,
              id_cols = ko) %>%
  separate(ko, sep = "ATAC_", into = c("rest", "ko")) %>%
  separate(ko, sep = ".mRp.clN.uniq", into = c("ko", "rest")) %>%
  dplyr::select(-rest)

pgl = ggplot(sc_input_pgl, aes(x = TLC, y = MeLC, color = ko)) +
  geom_point(size = 7) +
  scale_color_brewer(palette = "Set3") +
  labs(
    title = "mean bulk ATAC signals (PGL) over peak sets.",
    x = "TLC",
    y = "MeLC",
    color = " "
  ) +
  xlim(0,7) +
  ylim(0,4) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"))
ggsave(
  glue("{result_folder}mean_bulkATAC_PGL_over_peaksets.png"),
  plot = pgl,
  width = 6,
  height = 6,
  dpi = 500,
)
ggsave(
  glue("{result_folder}mean_bulkATAC_PGL_over_peaksets.pdf"),
  plot = pgl,
  width = 6,
  height = 6
)

sc_input_yapi = rbind(sc_tlc_yapi, sc_melc_yapi) %>%
  pivot_wider(.,
              names_from = cell_type,
              values_from = mean,
              id_cols = ko) %>%
  separate(ko, sep = "ATAC_", into = c("rest", "ko")) %>%
  separate(ko, sep = "_R[0-9]_RPGC", into = c("ko", "rest")) %>%
  dplyr::select(-rest)

yapi = ggplot(sc_input_yapi, aes(x = TLC, y = MeLC, color = ko)) +
  geom_point(size = 2) +
  scale_color_brewer(palette = "Set3") +
  labs(
    title = "mean bulk ATAC signals (YAPi) over peak sets.",
    x = "TLC",
    y = "MeLC",
    color = " "
  ) +
  xlim(0,5) +
  ylim(0,5) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"))
ggsave(
  glue("{result_folder}mean_bulkATAC_YAPi_over_peaksets.png"),
  plot = yapi,
  width = 6,
  height = 6,
  dpi = 500,
)
ggsave(
  glue("{result_folder}mean_bulkATAC_YAPi_over_peaksets.pdf"),
  plot = yapi,
  width = 6,
  height = 6
)

ggarrange(kos, pgl, yapi)
ggsave(
  glue("{result_folder}mean_bulkATAC_over_peaksets-all.png"),
  plot = last_plot(),
  width = 12,
  height = 10,
  dpi = 500,
)
ggsave(
  glue("{result_folder}mean_bulkATAC_over_peaksets-all.pdf"),
  plot = last_plot(),
  width = 12,
  height = 10
)


# heatmap
# bulk ATAC signals over scATAC peaks
png(
  file = glue("{result_folder}TLC_spec_peaks_Chens-KO_exp_signals.png"),
  width = 8,
  height = 13,
  units = 'cm',
  res = 500
)
col = colorRamp2(c(0, 10, 20), c(viridis(3)[1], viridis(3)[2], viridis(3)[3]))
hm_tlc = Heatmap(
  tlc_mat,
  column_title = "",
  column_title_side = "bottom",
  row_title = "TLC peak",
  name = "aggr. signal",
  #row_km = 2,
  #column_km = 2,
  clustering_method_rows = "complete",
  col = col,
  show_column_dend = TRUE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  heatmap_width = unit(5, "cm"),
  heatmap_height = unit(12, "cm"),
  show_row_names = FALSE,
  #row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 90
)
print(hm_tlc)
dev.off()

pdf(
  file = glue("{result_folder}TLC_spec_peaks_Chens-KO_exp_signals.pdf"),
  width = 4,
  height = 5
)
print(hm_tlc)
dev.off()

png(
  file = glue("{result_folder}MeLC_spec_peaks_Chens-KO_exp_signals.png"),
  width = 8,
  height = 13,
  units = 'cm',
  res = 500
)
col = colorRamp2(c(0, 10, 20), c(viridis(3)[1], viridis(3)[2], viridis(3)[3]))
hm_melc = Heatmap(
  melc_mat,
  column_title = "",
  column_title_side = "bottom",
  row_title = "MeLC peak",
  name = "aggr. signal",
  #row_km = 2,
  #column_km = 2,
  clustering_method_rows = "complete",
  col = col,
  show_column_dend = TRUE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  heatmap_width = unit(5, "cm"),
  heatmap_height = unit(12, "cm"),
  show_row_names = FALSE,
  #row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 90
)
print(hm_melc)
dev.off()

pdf(
  file = glue("{result_folder}MeLC_spec_peaks_Chens-KO_exp_signals.pdf"),
  width = 4,
  height = 5
)
print(hm_melc)
dev.off()

png(
  file = glue("{result_folder}aELC_spec_peaks_Chens-KO_exp_signals.png"),
  width = 8,
  height = 13,
  units = 'cm',
  res = 500
)
col = colorRamp2(c(0, 10, 20), c(viridis(3)[1], viridis(3)[2], viridis(3)[3]))
hm_aelc = Heatmap(
  aelc_mat,
  column_title = "",
  column_title_side = "bottom",
  row_title = "aELC peak",
  name = "aggr. signal",
  #row_km = 2,
  #column_km = 2,
  clustering_method_rows = "complete",
  col = col,
  show_column_dend = TRUE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  heatmap_width = unit(5, "cm"),
  heatmap_height = unit(12, "cm"),
  show_row_names = FALSE,
  #row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 90
)
print(hm_aelc)
dev.off()

pdf(
  file = glue("{result_folder}aELC_spec_peaks_Chens-KO_exp_signals.pdf"),
  width = 4,
  height = 5
)
print(hm_melc)
dev.off()