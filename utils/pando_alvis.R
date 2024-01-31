.libPaths(c(.libPaths(), "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/"))

library("cli", lib.loc = "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/")
library("Seurat", lib.loc = "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/")
library("Signac", lib.loc = "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/")
library("SeuratObject", lib.loc = "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/")

suppressPackageStartupMessages({
  library("tidyverse", lib.loc = "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/")
  library("Pando", lib.loc = "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/") # Pando requires Seurat V4!!! (e.g. Seurat_4.3.0) and SeuratObject v4
  library("GenomicRanges", lib.loc = "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/")
  library("BSgenome.Hsapiens.UCSC.hg38", lib.loc = "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/")
  library("glue", lib.loc = "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/")
  library("data.table", lib.loc = "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/")
  library("doParallel", lib.loc = "/mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils/packages/")
})

# handling v4 Seurat objects
options(Seurat.object.assay.version = "v4")

# candidates
candidate_tfs = fread("../data/GRN/candidate_tfs-Tobias_results.txt")
candidate_tfs = candidate_tfs %>% pull(candidate_TF) %>% unique 

# get annotation
annotation = readRDS("../results/EnsDB_V86_annotation.Rds")
#saveRDS(annotation, "../results/EnsDB_V86_annotation.Rds")

# knn pseudocell table for mimicking coupled multiome sequencing
nt_pseudo = fread("../results/Seurat_integration/nt_pseudocell_alignment.tsv")
nt_pseudo = nt_pseudo %>% distinct(id2, .keep_all = TRUE) # we should find out a ranking rather than distinct function!

trt_pseudo = fread("../results/Seurat_integration/trt_pseudocell_alignment.tsv")
trt_pseudo = trt_pseudo %>% distinct(id2, .keep_all = TRUE) # we should find out a ranking rather than distinct function!

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
nt_modified_fragments = "../data/scATAC-Seq/nt_fragments_pseudocells.tsv.gz"
trt_modified_fragments = "../data/scATAC-Seq/trt_fragments_pseudocells.tsv.gz"

FilterCells(
  "../data/scATAC-Seq/hESC_NT/fragments.tsv.gz",
  nt_ids,
  outfile = nt_modified_fragments,
  verbose = TRUE
)

FilterCells(
  "../data/scATAC-Seq/hESC_NT/fragments.tsv.gz",
  trt_ids,
  outfile = trt_modified_fragments,
  verbose = TRUE
)
# atac_meta = read.csv(
#   file = "../data/scATAC-Seq/hESC_NT/singlecell.csv",
#   header = TRUE,
#   row.names = 1
# )

nt_atac_counts = readRDS(file = "../results/scATAC-Seq/nt_ATAC-filtered_peak_bc_matrix_h5_object.Rds")
nt_atac_counts = nt_atac_counts[,nt_ids]
trt_atac_counts = readRDS(file = "../results/scATAC-Seq/trt_ATAC-filtered_peak_bc_matrix_h5_object.Rds")
trt_atac_counts = trt_atac_counts[,trt_ids]

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
saveRDS(nt_coembed, "../results/GRN/Pando/nt_Pando_input.Rds")
saveRDS(trt_coembed, "../results/GRN/Pando/trt_Pando_input.Rds")

#data('phastConsElements20Mammals.UCSC.hg38')

### PANDO protocol 
# source: https://quadbio.github.io/Pando/articles/getting_started.html

DefaultAssay(nt_coembed) = "peaks"
nt_coembed[['RNA']]
nt_coembed[['peaks']]

DefaultAssay(trt_coembed) = "peaks"
trt_coembed[['RNA']]
trt_coembed[['peaks']]

nt_gene_annot = Signac::Annotation(nt_coembed[["peaks"]])
trt_gene_annot = Signac::Annotation(trt_coembed[["peaks"]])

# 1.) initiating the GRN
nt_init_grn = initiate_grn(
  nt_coembed,
  rna_assay = "RNA",
  peak_assay = "peaks",
  exclude_exons = FALSE
)

trt_init_grn = initiate_grn(
  trt_coembed,
  rna_assay = "RNA",
  peak_assay = "peaks",
  exclude_exons = FALSE
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
saveRDS(nt_grn_motif, "../results/GRN/Pando/nt_Pando_findmotif.Rds")
nt_grn_motif = readRDS( "../results/GRN/Pando/nt_Pando_findmotif.Rds")

saveRDS(trt_grn_motif, "../results/GRN/Pando/trt_Pando_findmotif.Rds")
trt_grn_motif = readRDS( "../results/GRN/Pando/trt_Pando_findmotif.Rds")

# 3.) infer gene regulatory network (long run!!!)
registerDoParallel(4)

# run GRN inference on the most variable scRNA-Seq genes
nt_var_genes = FindVariableFeatures(nt_coembed[["RNA"]], nfeatures = 100)@var.features
print(nt_var_genes)
trt_var_genes = FindVariableFeatures(trt_coembed[["RNA"]], nfeatures = 100)@var.features
print(trt_var_genes)

nt_grn = infer_grn(
  nt_grn_motif,
  peak_to_gene_method = 'GREAT',
  genes = nt_var_genes,
  parallel = TRUE
)

trt_grn = infer_grn(
  trt_grn_motif,
  peak_to_gene_method = 'GREAT',
  genes = trt_var_genes,
  parallel = TRUE
)

# export intermediate object
saveRDS(nt_grn, "../results/GRN/Pando/nt_Pando_GRN.Rds")
saveRDS(trt_grn, "../results/GRN/Pando/trt_Pando_GRN.Rds")

GetNetwork(nt_grn)
nt_coefs = coef(nt_grn)
write_tsv(nt_coefs, "../results/GRN/Pando/nt_Pando_GRN_coefficients.tsv")

GetNetwork(trt_grn)
trt_coefs = coef(trt_grn)
write_tsv(trt_coefs, "../results/GRN/Pando/trt_Pando_GRN_coefficients.tsv")

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
write_tsv(nt_module_meta, "../results/GRN/Pando/nt_Pando_GRN_modules.tsv")

trt_modules = find_modules(
  trt_grn, 
  p_thresh = 1,
  nvar_thresh = 0, 
  min_genes_per_module = 0, 
  rsq_thresh = 0.1
)

trt_module_outputs = NetworkModules(trt_modules) 
trt_module_meta = trt_module_outputs@meta
write_tsv(trt_module_meta, "../results/GRN/Pando/nt_Pando_GRN_modules.tsv")

# 4.) Visualizations
# goodness-of-fit metrics
plot_gof(nt_modules, point_size=3)
pdf(
  file = "../results/GRN/Pando/nt_candidate_tfs-module_metrics_plot.pdf",
  width = 12,
  height = 6
)
plot_module_metrics(nt_modules)
dev.off()

plot_gof(trt_modules, point_size=3)
pdf(
  file = "../results/GRN/Pando/trt_candidate_tfs-module_metrics_plot.pdf",
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
  file = "../results/GRN/Pando/nt_candidate_tfs-network_grap-fr.pdf",
  width = 6,
  height = 6
)
plot_network_graph(nt_network_vis, layout='fr', text_size = 14)
dev.off()

trt_network_vis = get_network_graph(trt_modules, umap_method = "none")
pdf(
  file = "../results/GRN/Pando/trt_candidate_tfs-network_grap-fr.pdf",
  width = 6,
  height = 6
)
plot_network_graph(trt_network_vis, layout='fr', text_size = 14)
dev.off()

