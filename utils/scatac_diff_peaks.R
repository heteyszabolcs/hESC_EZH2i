# packages
suppressPackageStartupMessages({
  library("Seurat")
  library("Signac")
  library("glue")
  library("tidyverse")
  library("rtracklayer")
  library("GenomicFeatures")
  library("GenomeInfoDb")
  library("data.table")
})

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
    assay = "GA",
    group.by = "predicted.id_cheng",
    test.use = "LR",
    latent.vars = "peak_region_fragments"
  )

  output = differences %>% mutate(name = rownames(differences)) %>%
    dplyr::select(name, everything())

  return(output)
}

# diff_ga = marker_analysis(cluster = "MeLC", reference = NULL) 
# diff_ga %>% write_tsv(., glue("{result_folder}Diff_GA-LR_test-MeLC_vs_all-Chens_annotation.tsv"))

marker_analysis(cluster = "gELC", reference = "aELC") %>% write_tsv(.,
          glue("{result_folder}Diff_peak-LR_test-gELC_vs_aELC-Chens_annotation.tsv"))

## output exploration
melc = fread("../results/Seurat_integration/scATAC_annotation/Diff_peak-LR_test-MeLC_vs_all-Chens_annotation.tsv")
melc_opened = melc %>% dplyr::filter(avg_log2FC > 2 & p_val_adj < 0.05) %>% 
  separate(name, into = c("V1","V2","V3"), sep = "-") %>% 
  dplyr::select(V1, V2, V3) %>% 
  write_tsv(glue("{result_folder}MeLC_spec_peaks_up-log2fc_2_adjp0.05_Chens_annotation.bed"), col_names = FALSE)
melc_closed = melc %>% dplyr::filter(avg_log2FC < -2 & p_val_adj < 0.05) %>% 
  separate(name, into = c("V1","V2","V3"), sep = "-") %>% 
  dplyr::select(V1, V2, V3) %>% 
  write_tsv(glue("{result_folder}MeLC_spec_peaks_down-log2fc_2_adjp0.05_Chens_annotation.bed"), col_names = FALSE)

tlc = fread("../results/Seurat_integration/scATAC_annotation/Diff_peak-LR_test-TLC_vs_all-Chens_annotation.tsv")
tlc_opened = tlc %>% dplyr::filter(avg_log2FC > 2 & p_val_adj < 0.05) %>% 
  separate(name, into = c("V1","V2","V3"), sep = "-") %>% 
  dplyr::select(V1, V2, V3) %>% 
  write_tsv(glue("{result_folder}TLC_spec_peaks_up-log2fc_2_adjp0.05_Chens_annotation.bed"), col_names = FALSE)
tlc_closed = tlc %>% dplyr::filter(avg_log2FC < -2 & p_val_adj < 0.05) %>% 
  separate(name, into = c("V1","V2","V3"), sep = "-") %>% 
  dplyr::select(V1, V2, V3) %>% 
  write_tsv(glue("{result_folder}TLC_spec_peaks_down-log2fc_2_adjp0.05_Chens_annotation.bed"), col_names = FALSE)

active_elc = fread("../results/Seurat_integration/scATAC_annotation/Diff_peak-LR_test-aELC_vs_all-Chens_annotation.tsv")
active_elc_opened = active_elc %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.05) %>% 
  separate(name, into = c("V1","V2","V3"), sep = "-") %>% 
  dplyr::select(V1, V2, V3) %>% 
  write_tsv(glue("{result_folder}aELC_spec_peaks_up-log2fc_1_adjp0.05_Chens_annotation.bed"), col_names = FALSE)
active_elc_closed = active_elc %>% dplyr::filter(avg_log2FC < -1 & p_val_adj < 0.05) %>% 
  separate(name, into = c("V1","V2","V3"), sep = "-") %>% 
  dplyr::select(V1, V2, V3) %>% 
  write_tsv(glue("{result_folder}aELC_spec_peaks_down-log2fc_1_adjp0.05_Chens_annotation.bed"), col_names = FALSE)

ground_elc = fread("../results/Seurat_integration/scATAC_annotation/Diff_peak-LR_test-gELC_vs_all-Chens_annotation.tsv")
ground_elc_opened = ground_elc %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.05) %>% 
  separate(name, into = c("V1","V2","V3"), sep = "-") %>% 
  dplyr::select(V1, V2, V3) %>% 
  write_tsv(glue("{result_folder}gELC_spec_peaks_up-log2fc_1_adjp0.05_Chens_annotation.bed"), col_names = FALSE)
ground_elc_closed = ground_elc %>% dplyr::filter(avg_log2FC < -1 & p_val_adj < 0.05) %>% 
  separate(name, into = c("V1","V2","V3"), sep = "-") %>% 
  dplyr::select(V1, V2, V3) %>% 
  write_tsv(glue("{result_folder}gELC_spec_peaks_down-log2fc_1_adjp0.05_Chens_annotation.bed"), col_names = FALSE)

# genome browsing
# coverage plot function
make_coverage_plot = function(peak, label) {
  range = GRanges(
    seqnames = strsplit(peak, split = "-")[[1]][1],
    ranges = IRanges(
      start = as.numeric(strsplit(peak, split = "-")[[1]][2]),
      end = as.numeric(strsplit(peak, split = "-")[[1]][3]),
      names = "peak"
    )
  )
  
  plot = CoveragePlot(
    object = atac_trt,
    region = peak,
    ranges = range,
    ranges.title = "unique",
    assay = "peaks",
    annotation = TRUE,
    show.bulk = TRUE,
    ymax = 25,
    group.by = "predicted.id",
    peaks = TRUE,
    extend.upstream = 5000,
    extend.downstream = 5000
  )
  ggsave(
    glue("{result_folder}{label}_spec_scATAC_peak_{peak}.pdf"),
    plot = plot,
    width = 8,
    height = 6,
    device = "pdf"
  )
  
  return(print(plot))
  
}

# export plots of the top peaks
sapply(tlc %>% arrange(desc(avg_log2FC)) %>% top_n(5, wt = avg_log2FC) %>% pull(name),
       make_coverage_plot, label = "TLC")
sapply(melc %>% arrange(desc(avg_log2FC)) %>% top_n(5, wt = avg_log2FC) %>% pull(name),
       make_coverage_plot, label = "MeLC")

CoveragePlot(
  object = atac_trt,
  region = "chr10-122176886-122177269",
  ranges.title = "unique",
  assay = "peaks",
  annotation = TRUE,
  show.bulk = TRUE,
  ymax = 25,
  group.by = "predicted.id",
  peaks = TRUE,
  extend.upstream = 5000,
  extend.downstream = 5000
)

