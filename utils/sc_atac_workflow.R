print("Load R packages")
# packages
suppressPackageStartupMessages({
  library("Seurat")
  library("Signac")
  library("glue")
  library("tidyverse")
  library("rtracklayer")
  library("GenomicFeatures")
  library("GenomeInfoDb")
  library("ggpubr")
  library("EnsDb.Hsapiens.v86")
})

## scATAC-Seq workflow (Signac): https://stuartlab.org/signac/articles/pbmc_vignette

nt_cell_ranger = "../data/scATAC-Seq/hESC_NT/"
trt_cell_ranger = "../data/scATAC-Seq/hESC_EZH2i_7d/"

result_folder = "../results/scATAC-Seq/"

# helper function
read_seurat = function(cell_ranger_path) {
  counts = Read10X_h5(filename = glue("{cell_ranger_path}filtered_peak_bc_matrix.h5"))
  metadata = read.csv(
    file = glue("{cell_ranger_path}singlecell.csv"),
    header = TRUE,
    row.names = 1
  )
  
  chrom_assay = CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    fragments = glue("{cell_ranger_path}fragments.tsv.gz"),
    min.cells = 10,
    min.features = 200
  )
  
  seurat = CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata
  )
}

### non-treated workflow ###
# read into Seurat
nt_seurat = read_seurat(nt_cell_ranger)

# extract gene annotations from EnsDb
annotations = GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) = paste0('chr', seqlevels(annotations))
genome(annotations) = "hg38"
# add the gene information to the object
Annotation(nt_seurat) = annotations

## QC metrics
# compute nucleosome signal score per cell
nt_seurat = NucleosomeSignal(object = nt_seurat)
# compute TSS enrichment score per cell
nt_seurat = TSSEnrichment(object = nt_seurat, fast = TRUE)

nt_seurat$pct_reads_in_peaks = nt_seurat$peak_region_fragments / nt_seurat$passed_filters * 100

VlnPlot(
  object = nt_seurat,
  features = c('nCount_peaks', 'TSS.enrichment'),
  pt.size = 0.1,
  ncol = 2) 
VlnPlot(
  object = nt_seurat,
  features = c('nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 2) 

# filtering
nt_seurat = subset(
  x = nt_seurat,
  subset = nCount_peaks > 2000 &
    nCount_peaks < 30000 &
    pct_reads_in_peaks > 30 &
    nucleosome_signal < 2 &
    TSS.enrichment > 2
)

# normalization
nt_seurat = RunTFIDF(nt_seurat)
nt_seurat = FindTopFeatures(nt_seurat, min.cutoff = 'q0') # keep all features
nt_seurat = RunSVD(nt_seurat)

# UMAP
nt_seurat = RunUMAP(object = nt_seurat, reduction = 'lsi', dims = 2:30)
nt_seurat = FindNeighbors(object = nt_seurat, reduction = 'lsi', dims = 2:30)
nt_seurat = FindClusters(object = nt_seurat, verbose = FALSE, algorithm = 3, resolution = 0.2)

nt_umap = DimPlot(object = nt_seurat, label = TRUE) + 
  xlim(-8, 8) +
  ylim(-8, 8) +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("non-treated scATAC") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) 
nt_umap

# gene activities
nt_ga = GeneActivity(nt_seurat)
nt_seurat[['GA']] = CreateAssayObject(counts = nt_ga)
nt_seurat = NormalizeData(
  object = nt_seurat,
  assay = 'GA',
  normalization.method = 'LogNormalize',
  scale.factor = median(nt_seurat$nCount_GA)
)

# export Rds
saveRDS(nt_seurat, glue("{result_folder}hESC_scATAC_nt_Seurat.Rds"))

### 7d treated workflow ###
# read into Seurat
trt_seurat = read_seurat(trt_cell_ranger)
# add the gene information to the object
Annotation(trt_seurat) = annotations

## QC metrics
# compute nucleosome signal score per cell
trt_seurat = NucleosomeSignal(object = trt_seurat)
# compute TSS enrichment score per cell
trt_seurat = TSSEnrichment(object = trt_seurat, fast = TRUE)

trt_seurat$pct_reads_in_peaks = trt_seurat$peak_region_fragments / trt_seurat$passed_filters * 100

VlnPlot(
  object = trt_seurat,
  features = c('nCount_peaks', 'TSS.enrichment'),
  pt.size = 0.1,
  ncol = 2) 
VlnPlot(
  object = trt_seurat,
  features = c('nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 2) 

# filtering
trt_seurat = subset(
  x = trt_seurat,
  subset = nCount_peaks > 2000 &
    nCount_peaks < 30000 &
    pct_reads_in_peaks > 30 &
    nucleosome_signal < 2 &
    TSS.enrichment > 2
)

# normalization
trt_seurat = RunTFIDF(trt_seurat)
trt_seurat = FindTopFeatures(trt_seurat, min.cutoff = 'q0') # keep all features
trt_seurat = RunSVD(trt_seurat)

# UMAP
trt_seurat = RunUMAP(object = trt_seurat, reduction = 'lsi', dims = 2:30)
trt_seurat = FindNeighbors(object = trt_seurat, reduction = 'lsi', dims = 2:30)
trt_seurat = FindClusters(object = trt_seurat, verbose = FALSE, algorithm = 3, resolution = 0.3)

trt_umap = DimPlot(object = trt_seurat, label = TRUE) + 
  xlim(-8, 8) +
  ylim(-8, 8) +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("7d EZH2i scATAC") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) 
trt_umap


# gene activities
trt_ga = GeneActivity(trt_seurat)
trt_seurat[['GA']] = CreateAssayObject(counts = trt_ga)
trt_seurat = NormalizeData(
  object = trt_seurat,
  assay = 'GA',
  normalization.method = 'LogNormalize',
  scale.factor = median(trt_seurat$nCount_GA)
)

# export Rds
saveRDS(trt_seurat, glue("{result_folder}hESC_scATAC_7d_EZH2i_Seurat.Rds"))

umaps = ggarrange(nt_umap, trt_umap)
umaps

ggsave(
  glue("{result_folder}scATAC_Seq_UMAPs.png"),
  plot = umaps,
  width = 12,
  height = 6,
  dpi = 300,
)

ggsave(
  glue("{result_folder}scATAC_Seq_UMAPs.pdf"),
  plot = umaps,
  width = 12,
  height = 6,
  device = "pdf"
)






