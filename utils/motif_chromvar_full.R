suppressPackageStartupMessages({
  library("Signac")
  library("Seurat")
  library("BSgenome.Hsapiens.UCSC.hg38")
  library("patchwork")
  library("tidyverse")
  library("glue")
})

result_folder = "../results/motif_analysis/"

chromvar_analysis = function(scatac_seurat_object,
                             ident.1 = "TLC",
                             ident.2 = NULL,
                             group,
                             output_name) {
  library("JASPAR2020")
  library("TFBSTools")
  
  # PFMs (JASPAR2020)
  print("Create PFM...")
  pfm = getMatrixSet(
    x = JASPAR2020,
    opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
  )
  
  print("Motif addition...")
  seurat = AddMotifs(
    object = seurat[["peaks"]],
    genome = BSgenome.Hsapiens.UCSC.hg38,
    pfm = pfm
  )
  
  valid_feats = which(seqnames(seurat[['peaks']]@ranges) %in% 
                        standardChromosomes(seurat[['peaks']]@ranges))
  seurat = seurat[valid_feats, ]
  
  print("Run ChromVar...")
  seurat = RunChromVAR(
    object = seurat,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    assay = "motifs"
  )
  
  DefaultAssay(coembed_seurat) = 'chromvar'
  diff_activity = FindMarkers(
    object = coembed_seurat,
    ident.1 = ident.1,
    ident.2 = ident.2,
    group.by = group_name,
    only.pos = FALSE,
    mean.fxn = rowMeans,
    fc.name = "avg_diff"
  )
  
  diff_activity %>%
    mutate(gene_name = ConvertMotifID(
      seurat,
      assay = "motifs",
      id = rownames(diff_activity)
    )) %>% write_tsv(., 
                     output_name)
  
  return(diff_activity)
  
}

atac_trt = readRDS(glue("{result_folder}trt_scATAC-labeled.Rds"))
meta = atac_trt@meta.data

cheng = fread("../results/Seurat_integration/scATAC_annotation/Chen_scATAC_annot-extended.tsv")
cheng = cheng %>% separate(cell, sep = ".10X.", into = c("a", "id")) %>% 
  dplyr::select(-a) %>% 
  dplyr::filter(devTime == "EZH2i")
tlc = cheng %>% dplyr::filter(str_detect(cluster_EML, "TLC")) %>% pull(id)
melc = cheng %>% dplyr::filter(str_detect(cluster_EML, "MeLC")) %>% pull(id)

meta = meta %>% rownames_to_column(var = "cell_id") %>% 
  mutate(predicted.id_cheng = case_when(cell_id %in% tlc ~ "TLC",
                                        cell_id %in% melc ~ "MeLC",
                                        .default = predicted.id))
rownames(meta) = meta$cell_id
meta = meta %>% dplyr::select(-cell_id)
atac_trt@meta.data = meta

chromvar_analysis(scatac_seurat_object = atac_trt,
                  group = "predicted.id_cheng",
                  output_name = "../results/Seurat_integration/scATAC_annotation/scATAC_ChromVar_trt_TLC_vs_all.tsv")