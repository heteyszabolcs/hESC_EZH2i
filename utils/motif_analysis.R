suppressPackageStartupMessages({
  library("Signac")
  library("Seurat")
  library("BSgenome.Hsapiens.UCSC.hg38")
  library("patchwork")
  library("tidyverse")
  library("ggseqlogo")
  library("glue")
})

# for volcano
library("ggplot2")
library("ggrepel")
library("ggrastr")

result_folder = "../results/motif_analysis/"

# reading labeled scATAC objects and scATAC objects with motifs
atac_trt_seurat = readRDS("../results/Seurat_integration/scATAC_annotation/trt_scATAC-labeled.Rds")
atac_trt_motif = readRDS("../results/motif_analysis/trt_scATAC-motif.Rds")

atac_nt_seurat = readRDS("../results/Seurat_integration/scATAC_annotation/nt_scATAC-labeled.Rds")
atac_nt_motif = readRDS("../results/motif_analysis/nt_scATAC-motif.Rds")

atac_trt_seurat[["motifs"]] = atac_trt_motif
DefaultAssay(atac_trt_seurat) = "peaks"
atac_trt_seurat

atac_nt_seurat[["motifs"]] = atac_nt_motif
DefaultAssay(atac_nt_seurat) = "peaks"
atac_nt_seurat

coembed_seurat = readRDS("../results/scATAC-Seq/hESC_scATAC_nt_trt_coembedded.Rds")
coembed_motif = readRDS("../results/motif_analysis/trt_nt_coembed_scATAC-motif.Rds")
coembed_seurat[["motifs"]] = coembed_motif
DefaultAssay(coembed_seurat) = "peaks"

# differential motif analysis in EZH2i 7D treated scATAC data
# TLC vs. ELC motif analysis
elc_vs_tlc_da_trt_peaks = FindMarkers(
  object = atac_trt_seurat,
  ident.1 = "TLC",
  ident.2 = "ELC",
  group.by = "predicted.id",
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

elc_vs_tlc_top_da_trt_peak = rownames(elc_vs_tlc_da_trt_peaks[elc_vs_tlc_da_trt_peaks$p_val_adj < 0.05 & abs(elc_vs_tlc_da_trt_peaks$avg_log2FC) > 0.5, ])

atac_trt_seurat = RegionStats(
  object = atac_trt_seurat,
  assay = "motifs",
  genome = BSgenome.Hsapiens.UCSC.hg38
)

elc_vs_tlc_enriched_trt_motifs = FindMotifs(
  object = atac_trt_seurat,
  features = elc_vs_tlc_top_da_trt_peak,
  assay = "motifs"
)

elc_vs_tlc_top_enriched_trt_motifs = elc_vs_tlc_enriched_trt_motifs %>% arrange(p.adjust) %>% rownames %>% head(12)

require("ggseqlogo")
pdf(
  file = "../results/motif_analysis/trt_TLC_vs_ELC-top_motif_logos.pdf",
  width = 10,
  height = 4
)
MotifPlot(
  object = atac_trt_seurat,
  motifs = elc_vs_tlc_top_enriched_trt_motifs,
  assay = "motifs"
)
dev.off()

# MeLC vs. ELC motif analysis
elc_vs_melc_da_trt_peaks = FindMarkers(
  object = atac_trt_seurat,
  ident.1 = "MeLC",
  ident.2 = "ELC",
  group.by = "predicted.id",
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

elc_vs_melc_top_da_trt_peak = rownames(elc_vs_melc_da_trt_peaks[elc_vs_melc_da_trt_peaks$p_val_adj < 0.05 & abs(elc_vs_melc_da_trt_peaks$avg_log2FC) > 0.5, ])

elc_vs_melc_enriched_trt_motifs = FindMotifs(
  object = atac_trt_seurat,
  features = elc_vs_melc_top_da_trt_peak,
  assay = "motifs"
)

elc_vs_melc_top_enriched_trt_motifs = elc_vs_melc_enriched_trt_motifs %>% arrange(p.adjust) %>% rownames %>% head(12)

require("ggseqlogo")
pdf(
  file = "../results/motif_analysis/trt_MeLC_vs_ELC-top_motif_logos.pdf",
  width = 10,
  height = 4
)
MotifPlot(
  object = atac_trt_seurat,
  motifs = elc_vs_melc_top_enriched_trt_motifs,
  assay = "motifs"
)
dev.off()

## motif analysis on the coembedded scATAC-Seq samples (treated VS. non-treated)
trt_vs_nt_peaks = FindMarkers(
  object = coembed_seurat,
  ident.1 = "EZH2i_7D_scATAC_Seq",
  ident.2 = "non_trt_scATAC_Seq",
  group.by = "orig.ident",
  only.pos = FALSE,
  logfc.threshold = 0,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

trt_vs_nt_peaks_top_da_peak = rownames(trt_vs_nt_peaks[trt_vs_nt_peaks$p_val_adj < 0.05 & trt_vs_nt_peaks$avg_log2FC > 2, ])

coembed_seurat = RegionStats(
  object = coembed_seurat,
  assay = "motifs",
  genome = BSgenome.Hsapiens.UCSC.hg38
)

trt_vs_nt_enriched_motifs = FindMotifs(
  object = coembed_seurat,
  features = trt_vs_nt_peaks_top_da_peak,
  assay = "motifs"
)

trt_vs_nt_top_enriched_motifs = trt_vs_nt_enriched_motifs %>% arrange(p.adjust) %>% rownames %>% head(12)
trt_elc_candidates = trt_vs_nt_enriched_motifs %>% dplyr::filter(p.adjust < 0.05, fold.enrichment > 2) %>% pull(motif.name)


## visaulizations
volc_input = trt_vs_nt_enriched_motifs %>% mutate(group = case_when(
  fold.enrichment > 0.25 & p.adjust < 0.05 ~ "up", .default = "unaltered"
)) %>%
  mutate(sign_label = case_when(group == "up" ~ motif.name))
labels = volc_input %>% pull(sign_label)

# add color, size, alpha indications
cols = c("up" = "#fc9272", "unaltered" = "grey")
sizes = c("up" = 4, "unaltered" = 2)
alphas = c("up" = 1, "unaltered" = 0.5)

# EZH2i 7D motif plot
trt_elc_volc = volc_input %>%
  ggplot(aes(x = fold.enrichment,
             y = -log10(p.adjust),
             fill = group,    
             size = group,
             alpha = group)) +
  ggrastr::geom_point_rast(shape = 21, # Specify shape and colour as fixed local parameters    
                           colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = 0.25,
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(0, 5, 0.5)),  	 
                     limits = c(0, 5)) +
  scale_y_continuous(breaks = c(seq(0, 12, 2)),  	 
                     limits = c(0, 12)) +
  labs(
    title = "Enriched motifs in EZH2i 7D ELC cells",
    x = "motif fold enrichment",
    y = "-log10 adj.p-value",
    fill = " "
  ) +
  guides(alpha = FALSE, size = FALSE, fill = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  geom_text_repel(label = labels, size = 6) # add labels
trt_elc_volc

ggsave(
  glue("{result_folder}trt_ELC_motif_enrich_scatter.png"),
  plot = trt_elc_volc,
  width = 7,
  height = 7,
  dpi = 300,
)

ggsave(
  glue("{result_folder}trt_ELC_motif_enrich_scatter.pdf"),
  plot = trt_elc_volc,
  width = 7,
  height = 7,
  device = "pdf"
)

require("ggseqlogo")
pdf(
  file = "../results/motif_analysis/trt_ELC_vs_nt_ELC-top_motif_logos.pdf",
  width = 10,
  height = 4
)
MotifPlot(
  object = coembed_seurat,
  motifs = trt_vs_nt_top_enriched_motifs,
  assay = "motifs"
)
dev.off()



