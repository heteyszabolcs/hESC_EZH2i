# packages
suppressPackageStartupMessages({
  library("Seurat")
  library("Signac")
  library("glue")
  library("tidyverse")
  library("rtracklayer")
  library("GenomicFeatures")
  library("GenomeInfoDb")
})

result_folder = "../results/Seurat_integration/scATAC_annotation/"

# scATAC object (following multi-omic integration)
atac_trt = readRDS(glue("{result_folder}trt_scATAC-labeled.Rds"))

tlc_vs_all = FindMarkers(
  atac_trt,
  ident.1 = "TLC",
  ident.2 = NULL,
  only.pos = FALSE,
  assay = "peaks",
  group.by = "predicted.id",
  test.use = "LR", 
  latent.vars = "peak_region_fragments"
)

tlc_vs_all = tlc_vs_all %>% mutate(name = rownames(tlc_vs_all)) %>% 
  dplyr::select(name, everything())

write_tsv(tlc_vs_all,
          glue("{result_folder}Diff_peak-LR_test-TLC_vs_all.tsv"))


