if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  "tidyverse",
  "data.table",
  "ggplot2",
  "glue",
  "GenomicRanges",
  "GenomicFeatures",
  "bedscout",
  "RColorBrewer",
  "ComplexHeatmap",
  "circlize"
)

# output folder
result_folder = "../results/bulk_ATAC-Seq/"

# helper function: making GenomicRanges object
make_gr = function(path_to_bed, name) {
  bed = fread(path_to_bed)
  bed$V6 = name
  bed = GRanges(
    seqnames = bed$V1,
    ranges = IRanges(
      start = bed$V2,
      end = bed$V3,
      names = bed$V6
    )
  )
}

# create GRanges objects
CDX1 =
  make_gr(path_to_bed = "../data/bulk_ATAC-Seq/KO_experiments/MACS2/ATAC_CDX1_hESC_KO_EZH2i.mRp.clN_peaks.broadPeak", name = "CDX1_KO")
CDX2 =
  make_gr(path_to_bed = "../data/bulk_ATAC-Seq/KO_experiments/MACS2/ATAC_CDX2_hESC_KO_EZH2i.mRp.clN_peaks.broadPeak", name = "CDX2_KO")
HAND1 =
  make_gr(path_to_bed = "../data/bulk_ATAC-Seq/KO_experiments/MACS2/ATAC_HAND1_hESC_KO_EZH2i.mRp.clN_peaks.broadPeak", name = "HAND1_KO")
TBXT =
  make_gr(path_to_bed = "../data/bulk_ATAC-Seq/KO_experiments/MACS2/ATAC_TBXT_hESC_KO_EZH2i.mRp.clN_peaks.broadPeak", name = "TBXT_KO")
AP2A =
  make_gr(path_to_bed = "../data/bulk_ATAC-Seq/KO_experiments/MACS2/ATAC_AP2A_hESC_KO_EZH2i.mRp.clN_peaks.broadPeak", name = "AP2A_KO")
SOX9 =
  make_gr(path_to_bed = "../data/bulk_ATAC-Seq/KO_experiments/MACS2/ATAC_SOX9_hESC_KO_EZH2i.mRp.clN_peaks.broadPeak", name = "SOX9_KO")
MEIS2 =
  make_gr(path_to_bed = "../data/bulk_ATAC-Seq/KO_experiments/MACS2/ATAC_Meis2_hESC_KO_EZH2i.mRp.clN_peaks.broadPeak", name = "MEIS2_KO")
TBX3 =
  make_gr(path_to_bed = "../data/bulk_ATAC-Seq/KO_experiments/MACS2/ATAC_TBX3_hESC_KO_EZH2i.mRp.clN_peaks.broadPeak", name = "TBX3_KO")
WT_EZH2i =
  make_gr(path_to_bed = "../data/bulk_ATAC-Seq/KO_experiments/MACS2/ATAC_WT_hESC_EZH2i.mRp.clN_peaks.broadPeak", name = "WT_EZH2i")
WT_NT =
  make_gr(path_to_bed = "../data/bulk_ATAC-Seq/KO_experiments/MACS2/ATAC_WT_hESC_NT.mRp.clN_peaks.broadPeak", name = "WT_NT")

# Overlap analysis by bedscout
plot_euler(
  list(CDX1, CDX2, HAND1, TBXT),
  ignore.strand = TRUE,
  fills = c("#deebf7", "#9ecae1", "#3182bd"),
  names = c("CDX1 KO", "CDX2 KO", "HAND1 KO", "TBXT KO")
)
## scATAC-Seq diff GA
tlc_ga = "../results/Seurat_integration/scATAC_annotation/Diff_GA-LR_test-TLC_vs_all-Chens_annotation.tsv"
tlc_ga = fread(tlc_ga)
tlc_ga_sign = tlc_ga %>% dplyr::filter(avg_log2FC > 1 &
                                         p_val_adj < 0.05) %>% pull(name)

melc_ga = "../results/Seurat_integration/scATAC_annotation/Diff_GA-LR_test-MeLC_vs_all-Chens_annotation.tsv"
melc_ga = fread(melc_ga)
melc_ga_sign = melc_ga %>% dplyr::filter(avg_log2FC > 1 &
                                           p_val_adj < 0.05) %>% pull(name)

## scATAC-Seq peaks
# MeLC peaks
melc_peaks = "../results/Seurat_integration/scATAC_annotation/cluster_peaks/MeLC_spec_peaks_up-log2fc_2_adjp0.05_Chens_annotation.bed"
melc_peaks = rtracklayer::import(melc_peaks)
# TLC peaks
tlc_peaks = "../results/Seurat_integration/scATAC_annotation/cluster_peaks/TLC_spec_peaks_up-log2fc_2_adjp0.05_Chens_annotation.bed"
tlc_peaks = rtracklayer::import(tlc_peaks)
# aELC peaks
aelc_peaks = "../results/Seurat_integration/scATAC_annotation/cluster_peaks/aELC_spec_peaks_up-log2fc_1_adjp0.05_Chens_annotation.bed"
aelc_peaks = rtracklayer::import(aelc_peaks)
# gELC peaks
gelc_peaks = "../results/Seurat_integration/scATAC_annotation/cluster_peaks/gELC_spec_peaks_up-log2fc_1_adjp0.05_Chens_annotation.bed"
gelc_peaks = rtracklayer::import(gelc_peaks)

## bulk TF Cut&Tag peaks
# EZH2i 7D peaks
cdx1_cnt_trt = "../data/bulk_CutNTag/peaks/CDX1_EZH2i_7d_Rep1_peaks.narrowPeak"
cdx1_cnt_trt = rtracklayer::import(cdx1_cnt_trt)
cdx2_cnt_trt = "../data/bulk_CutNTag/peaks/CDX2_EZH2i_7d_Rep1_peaks.narrowPeak"
cdx2_cnt_trt = rtracklayer::import(cdx2_cnt_trt)
hand1_cnt_trt = "../data/bulk_CutNTag/peaks/HAND1_EZH2i_7d_Rep1_peaks.narrowPeak"
hand1_cnt_trt = rtracklayer::import(hand1_cnt_trt)
twist1_cnt_trt = "../data/bulk_CutNTag/peaks/Twist1_EZH2i_7d_Rep1_peaks.narrowPeak"
twist1_cnt_trt = rtracklayer::import(twist1_cnt_trt)
tbxt_cnt_trt = "../data/bulk_CutNTag/peaks/TBXT_EZH2i_7d_Rep1_peaks.narrowPeak"
tbxt_cnt_trt = rtracklayer::import(tbxt_cnt_trt)
tbx3_cnt_trt = "../data/bulk_CutNTag/peaks/TBX3_EZH2i_7d_Rep1_peaks.narrowPeak"
tbx3_cnt_trt = rtracklayer::import(tbx3_cnt_trt)
gata3_cnt_trt = "../data/bulk_CutNTag/peaks/GATA3_EZH2i_7d_Rep1_peaks.narrowPeak"
gata3_cnt_trt = rtracklayer::import(gata3_cnt_trt)
sox9_cnt_trt = "../data/bulk_CutNTag/peaks/SOX9_EZH2i_7d_Rep1_peaks.narrowPeak"
sox9_cnt_trt = rtracklayer::import(sox9_cnt_trt)
meis2_cnt_trt = "../data/bulk_CutNTag/peaks/Meis2_EZH2i_7d_Rep1_peaks.narrowPeak"
meis2_cnt_trt = rtracklayer::import(meis2_cnt_trt)
tfap2a_cnt_trt = "../data/bulk_CutNTag/peaks/TFAP2A_EZH2i_7d_Rep1_peaks.narrowPeak"
tfap2a_cnt_trt = rtracklayer::import(tfap2a_cnt_trt)
# non-treated peaks
cdx1_cnt = "../data/bulk_CutNTag/peaks/CDX1_NT_Rep1_peaks.narrowPeak"
cdx1_cnt = rtracklayer::import(cdx1_cnt)
cdx2_cnt = "../data/bulk_CutNTag/peaks/CDX2_NT_Rep1_peaks.narrowPeak"
cdx2_cnt = rtracklayer::import(cdx2_cnt)
hand1_cnt = "../data/bulk_CutNTag/peaks/HAND1_NT_Rep1_peaks.narrowPeak"
hand1_cnt = rtracklayer::import(hand1_cnt)
twist1_cnt = "../data/bulk_CutNTag/peaks/Twist1_NT_Rep1_peaks.narrowPeak"
twist1_cnt = rtracklayer::import(twist1_cnt)
tbxt_cnt = "../data/bulk_CutNTag/peaks/TBXT_NT_Rep1_peaks.narrowPeak"
tbxt_cnt = rtracklayer::import(tbxt_cnt)
tbx3_cnt = "../data/bulk_CutNTag/peaks/TBX3_NT_Rep1_peaks.narrowPeak"
tbx3_cnt = rtracklayer::import(tbx3_cnt)
gata3_cnt = "../data/bulk_CutNTag/peaks/GATA3_NT_Rep1_peaks.narrowPeak"
gata3_cnt = rtracklayer::import(gata3_cnt)
sox9_cnt = "../data/bulk_CutNTag/peaks/SOX9_NT_Rep1_peaks.narrowPeak"
sox9_cnt = rtracklayer::import(sox9_cnt)
meis2_cnt = "../data/bulk_CutNTag/peaks/Meis2_NT_Rep1_peaks.narrowPeak"
meis2_cnt = rtracklayer::import(meis2_cnt)
tfap2a_cnt = "../data/bulk_CutNTag/peaks/TFAP2A_NT_Rep1_peaks.narrowPeak"
tfap2a_cnt = rtracklayer::import(tfap2a_cnt)

## bedscout + visualizations
result_folder = "../results/Seurat_integration/scATAC_annotation/"
pdf(
  file = glue("{result_folder}scATAC_unique_peaks-Venn.pdf"),
  width = 5,
  height = 5
)
plot_euler(
  list(tlc_peaks, melc_peaks, aelc_peaks, gelc_peaks),
  ignore.strand = TRUE,
  fills = c("#deebf7", "#9ecae1", "#3182bd", "white"),
  names = c("TLC", "MeLC", "aELC", "gELC")
)
dev.off()

# RefSeq gene annotation
genes = "C:/Szabolcs/Karolinska/Data/reference_data/hs1_genes.ncbiRefSeq.bed"
genes = rtracklayer::import(genes)
tlc_ga_sign = genes[which(genes$name %in% tlc_ga_sign), ]
melc_ga_sign = genes[which(genes$name %in% melc_ga_sign), ]

# ChromHMM bed
hesc_chromhmm = "C:/Szabolcs/Karolinska/Data/reference_data/ChromHMM/ES-H1_E003_15_coreMarks.hg38.bed"
hesc_chromhmm = rtracklayer::import(hesc_chromhmm)

melc_peaks_chromhmm = annotate_nearby_features(
  melc_peaks,
  hesc_chromhmm,
  "name",
  distance_cutoff = 500,
  ignore.strand = FALSE
)
melc_peaks_chromhmm = as_tibble(melc_peaks_chromhmm) %>% mutate(cell_type = "MeLC") %>%
  group_by(cell_type, nearby_features) %>% count()
tlc_peaks_chromhmm = annotate_nearby_features(
  tlc_peaks,
  hesc_chromhmm,
  "name",
  distance_cutoff = 500,
  ignore.strand = FALSE
)
tlc_peaks_chromhmm = as_tibble(tlc_peaks_chromhmm) %>% mutate(cell_type = "TLC") %>%
  group_by(cell_type, nearby_features) %>% count()
aelc_peaks_chromhmm = annotate_nearby_features(
  aelc_peaks,
  hesc_chromhmm,
  "name",
  distance_cutoff = 500,
  ignore.strand = FALSE
)
aelc_peaks_chromhmm = as_tibble(aelc_peaks_chromhmm) %>% mutate(cell_type = "aELC") %>%
  group_by(cell_type, nearby_features) %>% count()


# ChromHMM heatmap on scATAC peaks
chromhmms = rbind(melc_peaks_chromhmm, tlc_peaks_chromhmm) %>%
  dplyr::filter(n > 10)
order = chromhmms %>% arrange(n) %>% pull(nearby_features)
order = factor(chromhmms$nearby_features, levels = unique(order))

hm1 = ggplot(chromhmms, aes(x = cell_type, y = order, fill = n)) +
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1) +
  xlab(label = "scATAC peak set") +
  ylab(label = "hESC ChromHMM") +
  labs(fill = "enrichment", subtitle = "ES-H1_E003_15_coreMarks") +
  theme(
    axis.text.x = element_text(
      color = "black",
      size = 15,
      angle = 90,
      hjust = 0.5,
      vjust = 0.5
    ),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title = element_text(size = 14),
    plot.subtitle = element_text(hjust = 1)
  ) +
  scale_fill_gradient2(
    low = "#075AFF",
    mid = "#FFFFCC",
    high = "#FF0000"
  ) +
  geom_text(aes(label = n), color = "black", size = 6) +
  coord_fixed()

hm1

ggsave(
  glue("{result_folder}scATAC_TLC_MeLC_hg38_hESC_ChromHMMs.png"),
  plot = hm1,
  width = 12,
  height = 10,
  dpi = 300,
)
ggsave(
  glue("{result_folder}scATAC_TLC_MeLC_hg38_hESC_ChromHMM.pdf"),
  plot = hm1,
  width = 12,
  height = 10
)

# bulk ATAC KO experiments
plot_pairwise_score(
  list(
    "CDX1 KO" = CDX1,
    "CDX2 KO" = CDX2,
    "HAND1 KO" = HAND1,
    "TBXT KO" = TBXT,
    "SOX9 KO" = SOX9,
    "TBX3 KO" = TBX3,
    "AP2A KO" = AP2A,
    "MEIS2 KO" = MEIS2,
    "WT EZH2i" = WT_EZH2i,
    "WT non-trt" = WT_NT
  ),
  list("MeLC GA" = melc_ga_sign, "TLC GA" = tlc_ga_sign),
  score_func = loci_overlap
)
plot_pairwise_score(
  list(
    "CDX1 KO EZH2i" = CDX1,
    "CDX2 KO EZH2i" = CDX2,
    "HAND1 KO EZH2i" = HAND1,
    "TBXT KO EZH2i" = TBXT,
    "SOX9 KO EZH2i" = SOX9,
    "TBX3 KO EZH2i" = TBX3,
    "AP2A KO EZH2i" = AP2A,
    "MEIS2 KO EZH2i" = MEIS2,
    "WT EZH2i" = WT_EZH2i,
    "WT non-trt" = WT_NT
  ),
  list(
    "MeLC" = melc_peaks,
    "TLC" = tlc_peaks,
    "gELC" = gelc_peaks,
    "aELC" = aelc_peaks
  ),
  score_func = loci_overlap
) + theme_minimal() +
  labs(y = "cluster peak set (EZH2i scATAC)", x = "MACS2 peak set (bulk ATAC)") +
  geom_text(color = "black") +
  scale_fill_gradient2(
    low = "#075AFF",
    mid = "#FFFFCC",
    high = "#FF0000"
  ) +
  ggtitle("peak overlaps") +
  theme(
    plot.title = element_text(size = 10),
    axis.text.x = element_text(
      size = 12,
      color = "black",
      angle = 90
    ),
    axis.text.y = element_text(size = 12, color = "black")
  )

ggsave(
  glue("{result_folder}scATAC-peak_overlaps.png"),
  plot = last_plot(),
  width = 6,
  height = 5,
  dpi = 300,
)
ggsave(
  glue("{result_folder}scATAC-peak_overlaps.pdf"),
  plot = last_plot(),
  width = 6,
  height = 5
)

# bulk TF Cut&Tag peak overlap analysis
result_folder = "../results/bulk_CutNTag/"

bulk_cnt_scores = pairwise_score(
  list(
    "CDX1 EZH2i" = cdx1_cnt_trt,
    "CDX2 EZH2i" = cdx2_cnt_trt,
    "HAND1 EZH2i" = hand1_cnt_trt,
    "GATA3 EZH2i" = gata3_cnt_trt,
    "TBXT EZH2i" = tbxt_cnt_trt,
    "TBX3 EZH2i" = tbx3_cnt_trt,
    "TFAP2A EZH2i" = tfap2a_cnt_trt,
    "MEIS2 EZH2i" = meis2_cnt_trt,
    "SOX9 EZH2i" = sox9_cnt_trt,
    "TWIST1 EZH2i" = twist1_cnt_trt,
    "CDX1 non-trt" = cdx1_cnt,
    "CDX2 non-trt" = cdx2_cnt,
    "HAND1 non-trt" = hand1_cnt,
    "GATA3 non-trt" = gata3_cnt,
    "TBXT non-trt" = tbxt_cnt,
    "TBX3 non-trt" = tbx3_cnt,
    "TFAP2A non-trt" = tfap2a_cnt,
    "MEIS2 non-trt" = meis2_cnt,
    "SOX9 non-trt" = sox9_cnt,
    "TWIST1 non-trt" = twist1_cnt),
    list(
      "MeLC" = melc_peaks,
      "TLC" = tlc_peaks,
      "gELC" = gelc_peaks,
      "aELC" = aelc_peaks
    ),
    score_func = loci_overlap
  )
bulk_cnt_scores = bulk_cnt_scores %>% pivot_wider(., names_from = gr1, id_cols = gr2, 
                                                  values_from = score) 
rows = bulk_cnt_scores$gr2 
bulk_cnt_scores = bulk_cnt_scores %>% dplyr::select(-gr2)
bulk_cnt_scores = as.matrix(bulk_cnt_scores)
rownames(bulk_cnt_scores) = rows

png(
  file = glue("{result_folder}bulk_CutNTag-peak_overlaps.png"),
  width = 15,
  height = 8,
  units = 'cm',
  res = 500
)
cols = colorRamp2(c(0, 5, 10), c("#9ecae1", "white", "#fc9272"))
cnt_hm = Heatmap(
  log2(bulk_cnt_scores+1),
  name = "log2(overlap)",
  # row_km = 2,
  # column_km = 1,
  clustering_method_columns = "complete",
  col = cols,
  rect_gp = gpar(col = "black", lwd = 0.5),
  show_column_dend = TRUE,
  cluster_columns = TRUE,
  cluster_rows = FALSE,
  show_row_dend = FALSE,
  heatmap_width = unit(10, "cm"),
  heatmap_height = unit(5, "cm"),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_title_side  = "bottom",
  column_names_rot = 90,
  column_title="TF Cut&Tag",
  row_title="cluster peak set (EZH2i scATAC)",
  row_title_gp = gpar(fontsize = 8)
)
print(cnt_hm)
dev.off()

pdf(
  file = glue("{result_folder}bulk_CutNTag-peak_overlaps.pdf"),
  width = 6,
  height = 3
)
print(cnt_hm)
dev.off()

bulk_cnt_only = pairwise_score(
  list(
    "CDX1 EZH2i" = cdx1_cnt_trt,
    "CDX2 EZH2i" = cdx2_cnt_trt,
    "HAND1 EZH2i" = hand1_cnt_trt,
    "GATA3 EZH2i" = gata3_cnt_trt,
    "TBXT EZH2i" = tbxt_cnt_trt,
    "TBX3 EZH2i" = tbx3_cnt_trt,
    "TFAP2A EZH2i" = tfap2a_cnt_trt,
    "MEIS2 EZH2i" = meis2_cnt_trt,
    "SOX9 EZH2i" = sox9_cnt_trt,
    "TWIST1 EZH2i" = twist1_cnt_trt),
  list(
    "CDX1 EZH2i" = cdx1_cnt_trt,
    "CDX2 EZH2i" = cdx2_cnt_trt,
    "HAND1 EZH2i" = hand1_cnt_trt,
    "GATA3 EZH2i" = gata3_cnt_trt,
    "TBXT EZH2i" = tbxt_cnt_trt,
    "TBX3 EZH2i" = tbx3_cnt_trt,
    "TFAP2A EZH2i" = tfap2a_cnt_trt,
    "MEIS2 EZH2i" = meis2_cnt_trt,
    "SOX9 EZH2i" = sox9_cnt_trt,
    "TWIST1 EZH2i" = twist1_cnt_trt),
  score_func = jaccard_index
)

bulk_cnt_only = bulk_cnt_only %>% pivot_wider(., names_from = gr1, id_cols = gr2, 
                                                  values_from = score) 
rows = bulk_cnt_only$gr2 
bulk_cnt_only = bulk_cnt_only %>% dplyr::select(-gr2)
bulk_cnt_only = as.matrix(bulk_cnt_only)
rownames(bulk_cnt_only) = rows

png(
  file = glue("{result_folder}bulk_EZH2i_CutNTag-peak_overlaps-Jaccard.png"),
  width = 12,
  height = 12,
  units = 'cm',
  res = 500
)
cols = colorRamp2(c(0, 0.2, 0.4), c("#9ecae1", "white", "#fc9272"))
cnt_hm2 = Heatmap(
  bulk_cnt_only,
  name = "Jaccard",
  # row_km = 2,
  # column_km = 1,
  clustering_method_columns = "complete",
  col = cols,
  rect_gp = gpar(col = "black", lwd = 0.5),
  show_column_dend = TRUE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  heatmap_width = unit(10, "cm"),
  heatmap_height = unit(10, "cm"),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_title_side  = "top",
  column_names_rot = 90,
  column_title="EZH2i TF Cut&Tag",
  row_title="",
  row_title_gp = gpar(fontsize = 8)
)
cnt_hm2
dev.off()

pdf(
  file = glue("{result_folder}bulk_EZH2i_CutNTag-peak_overlaps-Jaccard.pdf"),
  width = 5,
  height = 5
)
print(cnt_hm2)
dev.off()

pdf(
  file = glue("{result_folder}Venn-CDX2_CDX1_HAND1_nt.pdf"),
  width = 5,
  height = 5
)
plot_euler(
  list(cdx2_cnt, cdx1_cnt, hand1_cnt),
  ignore.strand = TRUE,
  fills = c("white", "#deebf7", "#9ecae1"),
  names = c("CDX2", "CDX1", "HAND1")
)
dev.off()

pdf(
  file = glue("{result_folder}Venn-CDX2_CDX1_HAND1_trt.pdf"),
  width = 5,
  height = 5
)
plot_euler(
  list(cdx2_cnt_trt, cdx1_cnt_trt, hand1_cnt_trt),
  ignore.strand = TRUE,
  fills = c("white", "#fee0d2", "#fc9272"),
  names = c("CDX2 EZH2i", "CDX1 EZH2i", "HAND1 EZH2i")
)
dev.off()

pdf(
  file = glue("{result_folder}Venn-TBX3_TBXT_TFAP2A_nt.pdf"),
  width = 5,
  height = 5
)
plot_euler(
  list(tbx3_cnt, tbxt_cnt, tfap2a_cnt),
  ignore.strand = TRUE,
  fills = c("white", "#deebf7", "#9ecae1"),
  names = c("TBX3", "TBXT", "TFAP2A")
)
dev.off()

pdf(
  file = glue("{result_folder}Venn-TBX3_TBXT_TFAP2A_trt.pdf"),
  width = 5,
  height = 5
)
plot_euler(
  list(tbx3_cnt_trt, tbxt_cnt_trt, tfap2a_cnt_trt),
  ignore.strand = TRUE,
  fills = c("white", "#fee0d2", "#fc9272"),
  names = c("TBX3", "TBXT", "TFAP2A")
)
dev.off()

pdf(
  file = glue("{result_folder}Venn-TBX3_TBXT_TFAP2A_HAND1_nt.pdf"),
  width = 5,
  height = 5
)
plot_euler(
  list(tbx3_cnt, tbxt_cnt, tfap2a_cnt, hand1_cnt),
  ignore.strand = TRUE,
  fills = c("white", "#deebf7", "#9ecae1", "#3182bd"),
  names = c("TBX3", "TBXT", "TFAP2A", "HAND1")
)
dev.off()

pdf(
  file = glue("{result_folder}Venn-TBX3_TBXT_TFAP2A_HAND1_trt.pdf"),
  width = 5,
  height = 5
)
plot_euler(
  list(tbx3_cnt_trt, tbxt_cnt_trt, tfap2a_cnt_trt, hand1_cnt_trt),
  ignore.strand = TRUE,
  fills = c("white", "#fee0d2", "#fc9272", "#de2d26"),
  names = c("TBX3", "TBXT", "TFAP2A", "HAND1")
)
dev.off()