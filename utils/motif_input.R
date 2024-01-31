suppressPackageStartupMessages({
  library("Signac")
  library("Seurat")
  library("BSgenome.Hsapiens.UCSC.hg38")
  library("patchwork")
  library("tidyverse")
  library("ggseqlogo")
})


# source: https://stuartlab.org/signac/articles/motif_vignette.html

atac_trt = readRDS("../results/Seurat_integration/scATAC_annotation/trt_scATAC-labeled.Rds")
atac_nt = readRDS("../results/Seurat_integration/scATAC_annotation/nt_scATAC-labeled.Rds")

set.seed(42)
atac_trt_ds = subset(x = atac_trt, downsample = 100)
atac_trt_ds_elc = subset(x = atac_trt_ds, subset = predicted.id == "ELC")
atac_trt_ds_elc@meta.data = atac_trt_ds_elc@meta.data %>% mutate(orig.ident = "EZH2i_7D_scATAC_Seq")
atac_nt_ds = subset(x = atac_nt, downsample = 100)
atac_nt_ds_elc = subset(x = atac_nt_ds, subset = predicted.id == "ELC")
atac_nt_ds_elc@meta.data = atac_nt_ds_elc@meta.data %>% mutate(orig.ident = "non_trt_scATAC_Seq")

coembed = merge(x = atac_nt_ds_elc, y = atac_trt_ds_elc)
saveRDS(coembed, "../results/scATAC-Seq/hESC_scATAC_nt_trt_coembedded.Rds")

library("JASPAR2020")
library("TFBSTools")
# PFMs
pfm = getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif informations
# EZH2i 7D treated
# atac_trt = AddMotifs(
#   object = atac_trt[["peaks"]],
#   genome = BSgenome.Hsapiens.UCSC.hg38,
#   pfm = pfm
# )
# 
# saveRDS(atac_trt, file = "../results/motif_analysis/trt_scATAC-motif.Rds")
# 
# # non-treated
# atac_nt = AddMotifs(
#   object = atac_nt[["peaks"]],
#   genome = BSgenome.Hsapiens.UCSC.hg38,
#   pfm = pfm
# )
# 
# saveRDS(atac_nt, file = "../results/motif_analysis/nt_scATAC-motif.Rds")

# coembed
coembed = AddMotifs(
  object = coembed[["peaks"]],
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

saveRDS(coembed, file = "../results/motif_analysis/trt_nt_coembed_scATAC-motif.Rds")
