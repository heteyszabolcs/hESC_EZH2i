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
  #library("DoubletFinder")
})

# export folder
result_folder = "../results/scRNA-Seq/"

# hESC scRNA-Seq data
counts = readRDS("../data/scRNA-Seq/Nerges.counts.filter.rds")

# keeping ELC cells, nontreated and D7 EZH2i treated
meta = readRDS("../data/scRNA-Seq/Nerges.meta.filter.rds")
meta = meta %>% dplyr::filter(cluster_EML != "Undef")
meta = meta %>% dplyr::filter(cellType == "EZH2i_Naive_WT" | cellType == "EZH2i_Naive_D7")
counts = counts[,meta$cell]

# read into Seurat
seurat_rna = CreateSeuratObject(counts = counts, project = "scRNA_EZH2i", min.cells = 3, min.features = 200)
seurat_rna@meta.data = cbind(seurat_rna@meta.data, meta)

## QC metrics
seurat_rna[["percent.mt"]] = PercentageFeatureSet(seurat_rna, pattern = "^MT-")
VlnPlot(seurat_rna, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# normalization (log)
seurat_rna = NormalizeData(seurat_rna, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_rna = FindVariableFeatures(seurat_rna, selection.method = "vst", nfeatures = 2000)

top10 = head(VariableFeatures(seurat_rna), 10)
highvar_plot1 = VariableFeaturePlot(seurat_rna)
highvar_plot2 = LabelPoints(plot = highvar_plot1, points = top10, repel = TRUE)
highvar_plot2 = highvar_plot2 + ggtitle("Highly variable genes")

# scaling and dim reduction
all.genes = rownames(seurat_rna)
seurat_rna = ScaleData(seurat_rna, features = all.genes)

seurat_rna = RunPCA(seurat_rna, features = VariableFeatures(object = seurat_rna))
DimHeatmap(seurat_rna, dims = 1, cells = 500, balanced = TRUE) # primary sources of heterogeneity in a dataset

ElbowPlot(seurat_rna)
seurat_rna = RunUMAP(seurat_rna, dims = 1:10)
seurat_rna = FindNeighbors(seurat_rna, dims = 1:10)
seurat_rna = FindClusters(seurat_rna, resolution = 0.4)


colnames(seurat_rna@meta.data)
dim_seurat_cls = DimPlot(object = seurat_rna, label = TRUE) + 
  xlim(-12, 12) +
  ylim(-12, 12) +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("group by Seurat clusters") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20, hjust = 0),
    plot.subtitle = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) 
dim_seurat_cls

dim_eml = DimPlot(object = seurat_rna, group.by = "cluster_EML", label = TRUE) + 
  xlim(-12, 12) +
  ylim(-12, 12) +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("group by EML") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20, hjust = 0),
    plot.subtitle = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) 
dim_eml

dim_cellt = DimPlot(object = seurat_rna, group.by = "cellType", label = TRUE) + 
  xlim(-12, 12) +
  ylim(-12, 12) +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("group by treatment") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20, hjust = 0),
    plot.subtitle = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) 
dim_cellt

umaps = ggarrange(dim_seurat_cls, dim_eml, dim_cellt)
umaps

ggsave(
  glue("{result_folder}scRNA_Seq_UMAPs-whole_set.png"),
  plot = umaps,
  width = 12,
  height = 10,
  dpi = 300,
)

ggsave(
  glue("{result_folder}scRNA_Seq_UMAPs-whole_set.pdf"),
  plot = umaps,
  width = 12,
  height = 10,
  device = "pdf"
)

# export Seurat object
saveRDS(seurat_rna, file = glue("{result_folder}hESC_EZH2i_scRNA_Seq.Rds"))

# marker analysis
tlc_vs_elc = FindMarkers(seurat_rna,
                         ident.1 = "TLC",
                         ident.2 = "ELC",
                         group.by = "cluster_EML")
tlc_vs_elc_sign = tlc_vs_elc %>% dplyr::filter(abs(avg_log2FC) > 2 & p_val_adj < 0.05)

melc_vs_elc = FindMarkers(seurat_rna,
                         ident.1 = "TLC",
                         ident.2 = "MeLC",
                         group.by = "cluster_EML")
melc_vs_elc_sign = melc_vs_elc %>% dplyr::filter(abs(avg_log2FC) > 2 & p_val_adj < 0.05)



