# packages on Alvis
.libPaths(c(.libPaths(), "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/"))

library("Seurat", lib.loc = "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/")
library("Signac", lib.loc = "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/")
library("tidyverse", lib.loc = "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/")
library("glue", lib.loc = "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/")

result_folder = "../results/Seurat_integration/scATAC_annotation/"

# scATAC object (following multi-omic integration)
atac_trt = readRDS(glue("{result_folder}trt_scATAC-labeled.Rds"))

# one vs all diff analysis
melc_vs_all = FindMarkers(
  atac_trt,
  ident.1 = "MeLC",
  ident.2 = NULL,
  only.pos = FALSE,
  assay = "peaks",
  group.by = "predicted.id",
  test.use = "LR", 
  latent.vars = "peak_region_fragments"
)

melc_vs_all = melc_vs_all %>% mutate(name = rownames(melc_vs_all)) %>% 
  dplyr::select(name, everything())

write_tsv(melc_vs_all,
          glue("{result_folder}Diff_peak-LR_test-MeLC_vs_all.tsv"))


