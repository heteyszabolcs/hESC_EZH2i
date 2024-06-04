if (!require("pacman"))
  install.packages("pacman")
pacman::p_load("tidyverse",
               "data.table",
               "ggplot2",
               "glue",
               "GenomicRanges",
               "GenomicFeatures",
               "bedscout"
)

# output folder
result_folder = "../results/bulk_ATAC-Seq/"

# helper function
make_gr = function(path_to_bed, name) {
  bed = fread(path_to_bed)
  bed$V6 = name
  bed = GRanges(seqnames = bed$V1,
                ranges = IRanges(
                  start = bed$V2,
                  end = bed$V3,
                  names = bed$V6
                ))
}

## create GRanges objects
CDX1 = 
  make_gr(path_to_bed = "../data/bulk_ATAC-Seq/KO_experiments/MACS2/ATAC_CDX1_hESC_KO_EZH2i.mRp.clN_peaks.broadPeak",
          name = "CDX1_KO")
CDX2 = 
  make_gr(path_to_bed = "../data/bulk_ATAC-Seq/KO_experiments/MACS2/ATAC_CDX2_hESC_KO_EZH2i.mRp.clN_peaks.broadPeak",
          name = "CDX2_KO")
HAND1 = 
  make_gr(path_to_bed = "../data/bulk_ATAC-Seq/KO_experiments/MACS2/ATAC_HAND1_hESC_KO_EZH2i.mRp.clN_peaks.broadPeak",
          name = "HAND1_KO")
TBXT = 
  make_gr(path_to_bed = "../data/bulk_ATAC-Seq/KO_experiments/MACS2/ATAC_TBXT_hESC_KO_EZH2i.mRp.clN_peaks.broadPeak",
          name = "TBXT_KO")

# Overlap analysis by bedscout
plot_euler(
  list(
    CDX1,
    CDX2,
    HAND1,
    TBXT
  ),
  ignore.strand = TRUE,
  fills = c("#deebf7", "#9ecae1", "#3182bd"),
  names = c("CDX1 KO", "CDX2 KO", "HAND1 KO", "TBXT KO")
)

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


plot_euler(
  list(
    tlc_peaks,
    melc_peaks,
    aelc_peaks
  ),
  ignore.strand = TRUE,
  fills = c("#deebf7", "#9ecae1", "#3182bd"),
  names = c("TLC", "MeLC", "aELC")
)

genes = "C:/Szabolcs/Karolinska/Data/reference_data/hs1_genes.ncbiRefSeq.bed"
genes = rtracklayer::import(genes)

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


# ChromHMM heatmap
result_folder = "../results/Seurat_integration/scATAC_annotation/"
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
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") +
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
