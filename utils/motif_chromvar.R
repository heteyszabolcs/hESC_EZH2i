suppressPackageStartupMessages({
  library("Signac")
  library("Seurat")
  library("BSgenome.Hsapiens.UCSC.hg38")
  library("patchwork")
  library("tidyverse")
  library("glue")
})

#devtools::install_github("ge11232002/TFBSTools")

coembed_seurat = readRDS("../results/scATAC-Seq/hESC_scATAC_nt_trt_coembedded.Rds")
coembed_motif = readRDS("../results/motif_analysis/trt_nt_coembed_scATAC-motif.Rds")
coembed_seurat[["motifs"]] = coembed_motif
DefaultAssay(object = coembed_seurat) = "peaks"


valid_feats = which(seqnames(coembed_seurat[['peaks']]@ranges) %in% standardChromosomes(coembed_seurat[['peaks']]@ranges))
coembed_seurat = coembed_seurat[valid_feats, ]
coembed_seurat = RunChromVAR(
  object = coembed_seurat,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay = "motifs"
)

DefaultAssay(coembed_seurat) = 'chromvar'

trt_vs_nt_ELC_diff_activity = FindMarkers(
  object = coembed_seurat,
  ident.1 = "EZH2i_7D_scATAC_Seq",
  ident.2 = "non_trt_scATAC_Seq",
  group.by = "orig.ident",
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

trt_vs_nt_ELC_top_activities = trt_vs_nt_ELC_diff_activity %>% arrange(desc(avg_diff)) %>% rownames %>% head(12)

require("ggseqlogo")
pdf(
  file = "../results/motif_analysis/trt_ELC_vs_nt_ELC-top_chromVar_act.pdf",
  width = 10,
  height = 4
)
MotifPlot(
  object = coembed_seurat,
  motifs = trt_vs_nt_ELC_top_activities,
  assay = "motifs"
)
dev.off()


