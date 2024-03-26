library("utils")
# packages
suppressPackageStartupMessages({
  library("tidyverse")
  library("Pando") 
  library("Seurat") 
  library("Signac") 
  library("GenomicRanges")
  library("EnsDb.Hsapiens.v86")
  library("BSgenome.Hsapiens.UCSC.hg38")
  library("glue")
  library("data.table")
  library("doParallel")
})

# result folder
result_folder = "../results/GRN/Pando/Var_TFs_of_Jaspar_CISBP/"

# candidates
tobias_tfs = fread("../data/GRN/TOBIAS_bulkATAC-Seq_candidates.txt")
tobias_tfs = tobias_tfs %>% pull(candidate_TF) %>% unique 

chromvar_enrichments = fread("../data/GRN/ChromVar_TOBIAS_trt_ELC_vs_nt_ELC-sign_genes.txt")
chromvar_enrichments = chromvar_enrichments %>% pull(gene_name) %>% unique 

scrna_trt_elc_vs_nt_elc_sign_up = fread("../data/GRN/scRNASeq-trt_ELC_vs_nt_ELC-sign_up_genes.txt")
scrna_trt_elc_vs_nt_elc_sign_up = scrna_trt_elc_vs_nt_elc_sign_up$gene

var_tfs_of_jaspar_cisbp = fread("../data/GRN/Variable_TFs-Jaspar2024_CISBP.tsv")
var_tfs_of_jaspar_cisbp = var_tfs_of_jaspar_cisbp$TF_name

# hg38 cis-regulatory elements from SCREEN (GenomicRanges object)
screen = data('SCREEN.ccRE.UCSC.hg38')
length(SCREEN.ccRE.UCSC.hg38)

# get annotation
annotation = GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#saveRDS(annotation, "../results/EnsDB_V86_annotation.Rds")

# knn pseudocell table for mimicking coupled multiome sequencing
# nt_pseudo_x = fread("../results/Seurat_integration/nt_pseudocell_alignment.tsv")
# nt_pseudo = nt_pseudo %>% distinct(id2, .keep_all = TRUE) 

# trt_pseudo = fread("../results/Seurat_integration/trt_pseudocell_alignment.tsv")
# trt_pseudo = trt_pseudo %>% distinct(id2, .keep_all = TRUE) 

# anchor score based pseudocells - FindTransferAnchors function 
nt_pseudo = fread("../results/Seurat_integration/nt_anchor_matrix.tsv")
nt_pseudo = nt_pseudo %>% dplyr::filter(score >= 0.3) %>% 
  dplyr::select(id1 = anchor1_barcode, id2 = anchor2_barcode)  %>% 
  distinct(id2, .keep_all = TRUE) %>% 
  distinct(id1, .keep_all = TRUE)

trt_pseudo = fread("../results/Seurat_integration/trt_anchor_matrix.tsv")
trt_pseudo = trt_pseudo %>% dplyr::filter(score >= 0.3) %>% 
  dplyr::select(id1 = anchor1_barcode, id2 = anchor2_barcode)  %>% 
  distinct(id2, .keep_all = TRUE) %>% 
  distinct(id1, .keep_all = TRUE)

# hESC scRNA-Seq data
nt_counts = readRDS("../data/scRNA-Seq/Nerges.counts.filter.rds")
nt_counts = nt_counts[,nt_pseudo$id1]
colnames(nt_counts) = nt_pseudo$id2

trt_counts = readRDS("../data/scRNA-Seq/Nerges.counts.filter.rds")
trt_counts = trt_counts[,trt_pseudo$id1]
colnames(trt_counts) = trt_pseudo$id2

meta = readRDS("../data/scRNA-Seq/Nerges.meta.filter.rds")
nt_meta = meta[which(meta$cell %in% nt_pseudo$id1),]
nt_meta$SID = nt_pseudo$id2
nt_meta$cell = nt_pseudo$id2

meta = readRDS("../data/scRNA-Seq/Nerges.meta.filter.rds")
trt_meta = meta[which(meta$cell %in% trt_pseudo$id1),]
trt_meta$SID = trt_pseudo$id2
trt_meta$cell = trt_pseudo$id2

nt_rna = CreateSeuratObject(counts = nt_counts, project = "scRNA_EZH2i",
                           min.cells = 3,
                           min.features = 200)
nt_rna@meta.data = cbind(nt_rna@meta.data, nt_meta)
trt_rna = CreateSeuratObject(counts = trt_counts, project = "scRNA_EZH2i",
                            min.cells = 3,
                            min.features = 200)
trt_rna@meta.data = cbind(trt_rna@meta.data, trt_meta)

nt_ids = nt_rna@meta.data %>% dplyr::filter(cellType == "EZH2i_Naive_WT") %>% pull(SID)
trt_ids = trt_rna@meta.data %>% dplyr::filter(cellType == "EZH2i_Naive_D7") %>% pull(SID)
nt_counts = nt_counts[,nt_ids]
trt_counts = trt_counts[,trt_ids]

nt_rna = CreateSeuratObject(counts = nt_counts, project = "scRNA_non_trt",
                         min.cells = 3,
                         min.features = 200)
trt_rna = CreateSeuratObject(counts = trt_counts, project = "scRNA_EZH2i",
                            min.cells = 3,
                            min.features = 200)

## create scATAC-Seq object
# filter barcodes by pseudocell technique
nt_modified_fragments = "../data/scATAC-Seq/hESC_NT/nt_fragments_pseudocells.tsv.gz"
trt_modified_fragments = "../data/scATAC-Seq/hESC_NT/trt_fragments_pseudocells.tsv.gz"

FilterCells(
  "../data/scATAC-Seq/hESC_NT/fragments.tsv.gz",
  nt_ids,
  outfile = nt_modified_fragments,
  verbose = TRUE
)

FilterCells(
  "../data/scATAC-Seq/hESC_EZH2i_7d/fragments.tsv.gz",
  trt_ids,
  outfile = trt_modified_fragments,
  verbose = TRUE
)

# atac_meta = read.csv(
#   file = "../data/scATAC-Seq/hESC_NT/singlecell.csv",
#   header = TRUE,
#   row.names = 1
# )

nt_atac_counts = Read10X_h5(filename = "../data/scATAC-Seq/hESC_NT/filtered_peak_bc_matrix.h5")
nt_atac_counts = nt_atac_counts[,nt_ids]
saveRDS(nt_atac_counts, file = "../results/scATAC-Seq/nt_ATAC-filtered_peak_bc_matrix_h5_object.Rds")
trt_atac_counts = Read10X_h5(filename = "../data/scATAC-Seq/hESC_EZH2i_7d/filtered_peak_bc_matrix.h5")
trt_atac_counts = trt_atac_counts[,trt_ids]
saveRDS(trt_atac_counts, file = "../results/scATAC-Seq/trt_ATAC-filtered_peak_bc_matrix_h5_object.Rds")

nt_chrom_assay = CreateChromatinAssay(
  counts = nt_atac_counts,
  sep = c(":", "-"),
  fragments = nt_modified_fragments,
  min.cells = 10,
  min.features = 200,
  annotation = annotation
)

trt_chrom_assay = CreateChromatinAssay(
  counts = trt_atac_counts,
  sep = c(":", "-"),
  fragments = trt_modified_fragments,
  min.cells = 10,
  min.features = 200,
  annotation = annotation
)

### normalize scRNA count table
nt_rna = NormalizeData(nt_rna, normalization.method = "LogNormalize", scale.factor = 10000)
nt_rna = FindVariableFeatures(nt_rna, selection.method = "vst", nfeatures = 2000)

trt_rna = NormalizeData(trt_rna, normalization.method = "LogNormalize", scale.factor = 10000)
trt_rna = FindVariableFeatures(trt_rna, selection.method = "vst", nfeatures = 2000)

# add chromatin assay to Seurat object
nt_rna[["peaks"]] = nt_chrom_assay
trt_rna[["peaks"]] = trt_chrom_assay

nt_coembed = nt_rna
trt_coembed = trt_rna

### TF-IDF normalize scATAC count table
nt_coembed = RunTFIDF(nt_coembed, assay = "peaks")
nt_coembed = FindTopFeatures(nt_coembed, min.cutoff = 'q0', assay = "peaks")
nt_coembed = RunSVD(nt_coembed, assay = "peaks")

trt_coembed = RunTFIDF(trt_coembed, assay = "peaks")
trt_coembed = FindTopFeatures(trt_coembed, min.cutoff = 'q0', assay = "peaks")
trt_coembed = RunSVD(trt_coembed, assay = "peaks")

# export PANDO input object
saveRDS(nt_coembed, glue("{result_folder}nt_Pando_input.Rds"))
saveRDS(trt_coembed, glue("{result_folder}trt_Pando_input.Rds"))

#data('phastConsElements20Mammals.UCSC.hg38')

### PANDO protocol 
# source: https://quadbio.github.io/Pando/articles/getting_started.html

DefaultAssay(nt_coembed) = "peaks"
nt_coembed[['RNA']]
nt_coembed[['peaks']]

nt_atac_peaks = nt_coembed@assays$peaks@ranges
nt_screen = findOverlaps(nt_atac_peaks, 
                  SCREEN.ccRE.UCSC.hg38, type = "any", ignore.strand = FALSE)
nt_screen = nt_atac_peaks[queryHits(nt_screen)] # nt SCREEN regions

DefaultAssay(trt_coembed) = "peaks"
trt_coembed[['RNA']]
trt_coembed[['peaks']]

trt_atac_peaks = trt_coembed@assays$peaks@ranges
trt_screen = findOverlaps(trt_atac_peaks, 
                         SCREEN.ccRE.UCSC.hg38, type = "any", ignore.strand = FALSE)
trt_screen = trt_atac_peaks[queryHits(trt_screen)] # trt SCREEN regions

nt_gene_annot = Signac::Annotation(nt_coembed[["peaks"]])
trt_gene_annot = Signac::Annotation(trt_coembed[["peaks"]])

# 1.) initiating the GRN
nt_init_grn = initiate_grn(
  nt_coembed,
  rna_assay = "RNA",
  peak_assay = "peaks",
  exclude_exons = FALSE,
  regions = nt_screen
)

trt_init_grn = initiate_grn(
  trt_coembed,
  rna_assay = "RNA",
  peak_assay = "peaks",
  exclude_exons = FALSE,
  regions = trt_screen
)

nt_candidate_regions = NetworkRegions(nt_init_grn)
nt_candidate_regions = nt_candidate_regions@ranges
trt_candidate_regions = NetworkRegions(trt_init_grn)
trt_candidate_regions = trt_candidate_regions@ranges

# TF motifs from PANDO package
data('motifs')
data('motif2tf')

# 2.) finding TF motifs (long run!!!)
nt_grn_motif = find_motifs(
  object = nt_init_grn, 
  pfm = motifs, 
  genome = BSgenome.Hsapiens.UCSC.hg38,
  verbose = TRUE
)

trt_grn_motif = find_motifs(
  object = trt_init_grn, 
  pfm = motifs, 
  genome = BSgenome.Hsapiens.UCSC.hg38,
  verbose = TRUE
)

# export intermediate object
saveRDS(nt_grn_motif, glue("{result_folder}nt_Pando_findmotif.Rds"))
nt_grn_motif = readRDS( glue("{result_folder}nt_Pando_findmotif.Rds"))

saveRDS(trt_grn_motif, glue("{result_folder}trt_Pando_findmotif.Rds"))
trt_grn_motif = readRDS(glue("{result_folder}trt_Pando_findmotif.Rds"))

# 3.) infer gene regulatory network (long run!!!)
registerDoParallel(4)

# run GRN inference on the most variable scRNA-Seq genes
nt_var_genes = VariableFeatures(FindVariableFeatures(nt_coembed[["RNA"]],
                                                     nfeatures = 100))
print(nt_var_genes)
trt_var_genes = VariableFeatures(FindVariableFeatures(trt_coembed[["RNA"]], 
                                                      nfeatures = 100))
print(trt_var_genes)


nt_grn = infer_grn(
  nt_grn_motif,
  peak_to_gene_method = 'GREAT',
  genes = var_tfs_of_jaspar_cisbp,
  parallel = TRUE
)

trt_grn = infer_grn(
  trt_grn_motif,
  peak_to_gene_method = 'GREAT',
  genes = var_tfs_of_jaspar_cisbp,
  parallel = TRUE
)

# export intermediate object
saveRDS(nt_grn, glue("{result_folder}nt_Pando_GRN.Rds"))
saveRDS(trt_grn, glue("{result_folder}trt_Pando_GRN.Rds"))

GetNetwork(nt_grn)
nt_coefs = coef(nt_grn)
write_tsv(nt_coefs, glue("{result_folder}nt_Pando_GRN_coefficients.tsv"))

GetNetwork(trt_grn)
trt_coefs = coef(trt_grn)
write_tsv(trt_coefs, glue("{result_folder}trt_Pando_GRN_coefficients.tsv"))

# 3.) Module discovery
nt_modules = find_modules(
  nt_grn, 
  p_thresh = 1,
  nvar_thresh = 0, 
  min_genes_per_module = 0, 
  rsq_thresh = 0.1
)

nt_module_outputs = NetworkModules(nt_modules) 
nt_module_meta = nt_module_outputs@meta
write_tsv(nt_module_meta, glue("{result_folder}nt_Pando_GRN_modules.tsv"))

trt_modules = find_modules(
  trt_grn, 
  p_thresh = 1,
  nvar_thresh = 0, 
  min_genes_per_module = 0, 
  rsq_thresh = 0.1
)

trt_module_outputs = NetworkModules(trt_modules) 
trt_module_meta = trt_module_outputs@meta
write_tsv(trt_module_meta, glue("{result_folder}nt_Pando_GRN_modules.tsv"))

# 4.) Visualizations
# goodness-of-fit metrics
pdf(
  file = glue("{result_folder}nt_candidate_tfs-goodness_of_fits.pdf"),
  width = 4,
  height = 4
)
plot_gof(nt_modules, point_size=3)
dev.off()

pdf(
  file = glue("{result_folder}nt_candidate_tfs-module_metrics_plot.pdf"),
  width = 12,
  height = 6
)
plot_module_metrics(nt_modules)
dev.off()

pdf(
  file = glue("{result_folder}trt_candidate_tfs-goodness_of_fits.pdf"),
  width = 4,
  height = 4
)
plot_gof(trt_modules, point_size=3)
dev.off()
pdf(
  file = glue("{result_folder}trt_candidate_tfs-module_metrics_plot.pdf"),
  width = 12,
  height = 6
)
plot_module_metrics(trt_modules)
dev.off()

# network
detach("package:BSgenome.Hsapiens.UCSC.hg38", unload=TRUE)
detach("package:EnsDb.Hsapiens.v86", unload=TRUE)
nt_network_vis = get_network_graph(nt_modules, umap_method = "none")
#network_vis = RunUMAP(object = network_vis, reduction = 'lsi', dims = 2:30, assay = "peaks")
#network_vis = RunUMAP(object = network_vis, reduction = 'lsi', dims = 2:30, assay = "RNA")
pdf(
  file = glue("{result_folder}nt_candidate_tfs-network_grap-fr.pdf"),
  width = 6,
  height = 6
)
plot_network_graph(nt_network_vis, layout='fr', text_size = 14)
dev.off()

trt_network_vis = get_network_graph(trt_modules, umap_method = "none")
pdf(
  file = glue("{result_folder}trt_candidate_tfs-network_grap-fr.pdf"),
  width = 6,
  height = 6
)
plot_network_graph(trt_network_vis, layout='fr', text_size = 14)
dev.off()





