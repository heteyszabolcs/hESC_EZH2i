# packages on Alvis
.libPaths(c(.libPaths(), "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/"))

library("Seurat", lib.loc = "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/")
library("Signac", lib.loc = "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/")
library("tidyverse", lib.loc = "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/")
library("glue", lib.loc = "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/")
library("data.table", lib.loc = "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/")


result_folder = "../results/Seurat_integration/scATAC_annotation/"

# scATAC object (following multi-omic integration)
atac_trt = readRDS(glue("{result_folder}trt_scATAC-labeled.Rds"))
meta = atac_trt@meta.data

# all_tlcs = fread(glue("{result_folder}scATAC-TLC_cells-union_w_Chengs.tsv"))
# all_melcs = fread(glue("{result_folder}scATAC-MeLC_cells-union_w_Chengs.tsv"))

cheng = fread("../results/Seurat_integration/scATAC_annotation/Chen_scATAC_annot-extended.tsv")
cheng = cheng %>% separate(cell, sep = ".10X.", into = c("a", "id")) %>% 
  dplyr::select(-a) %>% 
  dplyr::filter(devTime == "EZH2i")
tlc = cheng %>% dplyr::filter(str_detect(cluster_EML, "TLC")) %>% pull(id)
melc = cheng %>% dplyr::filter(str_detect(cluster_EML, "MeLC")) %>% pull(id)
aELC = cheng %>% dplyr::filter(str_detect(cluster_EML, "aELC")) %>% pull(id)
gELC = cheng %>% dplyr::filter(str_detect(cluster_EML, "gELC")) %>% pull(id)

meta = meta %>% rownames_to_column(var = "cell_id") %>% 
  mutate(predicted.id_cheng = case_when(cell_id %in% tlc ~ "TLC",
                                        cell_id %in% melc ~ "MeLC",
                                        cell_id %in% aELC ~ "aELC",
                                        cell_id %in% gELC ~ "gELC",
                                        .default = predicted.id))
rownames(meta) = meta$cell_id
meta = meta %>% dplyr::select(-cell_id)
atac_trt@meta.data = meta

# Cheng's annotation
# marker analysis
marker_analysis = function(cluster, reference) {
  differences = FindMarkers(
    atac_trt,
    ident.1 = cluster,
    ident.2 = reference,
    only.pos = FALSE,
    assay = "peaks",
    group.by = "predicted.id_cheng",
    test.use = "LR",
    latent.vars = "peak_region_fragments"
  )
  
  output = differences %>% mutate(name = rownames(differences)) %>%
    dplyr::select(name, everything())
  
  return(output)
}

marker_analysis(cluster = "gELC", reference = "aELC") %>% 
  write_tsv(., glue("{result_folder}Diff_peak-LR_test-gELC_vs_aELC-Chens_annotation.tsv"))
