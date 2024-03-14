# remove package versions
# remove.packages(c("Seurat", "Signac", "SeuratObject"))
# 
# # newest versions
# install.packages("Seurat")
# install.packages("Signac")

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

result_folder = "../results/Seurat_integration/"

## scATAC
atac_nt = "../results/scATAC-Seq/hESC_scATAC_nt_Seurat.Rds"
atac_nt = readRDS(atac_nt)

atac_trt = "../results/scATAC-Seq/hESC_scATAC_7d_EZH2i_Seurat.Rds"
atac_trt = readRDS(atac_trt)

DefaultAssay(atac_nt) = "GA"
DefaultAssay(atac_trt) = "GA"

### non-treated workflow ###
# scRNA-Seq (naive ELC cells)
rna = "../results/scRNA-Seq/hESC_EZH2i_scRNA_Seq.Rds"
rna = readRDS(rna)

# naive set
naive_rna = rna@meta.data %>% dplyr::filter(str_detect(cell, "EZH2i_Naive")) %>% pull(cell)
naive_rna = subset(rna, subset = cell %in% naive_rna)
naive_nt_rna = naive_rna@meta.data %>% dplyr::filter(str_detect(SID, "Naive_WT")) %>% pull(cell)
naive_nt_rna = subset(naive_rna, subset = cell %in% naive_nt_rna)
naive_trt_rna = naive_rna@meta.data %>% dplyr::filter(str_detect(SID, "Naive_D7")) %>% pull(cell)
naive_trt_rna = subset(naive_rna, subset = cell %in% naive_trt_rna)
rm(rna)

# ELCs
naive_trt_elcs = naive_trt_rna@meta.data %>% dplyr::filter(cluster_EML == "ELC") %>% rownames
naive_nt_elcs = naive_nt_rna@meta.data %>% dplyr::filter(cluster_EML == "ELC") %>% rownames

dim_eml = DimPlot(object = naive_rna, group.by = "SID", label = TRUE) + 
  xlim(-12, 12) +
  ylim(-12, 12) +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("group by treatment (naive subset)") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 20),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) 
dim_eml

## integration by CCA
# finding anchors
nt_anchors = FindTransferAnchors(
  reference = naive_nt_rna,
  query = atac_nt,
  features = VariableFeatures(object = naive_nt_rna),
  reference.assay = "RNA",
  query.assay = "GA",
  reduction = 'cca'
)

# extract anchor matrix of AnchorSet object
nt_anchor_mat = as_tibble(nt_anchors@anchors)
nt_anchor_mat = nt_anchor_mat %>%
  mutate(anchor1_barcode = unname(naive_nt_rna$cell)[nt_anchor_mat$cell1]) %>%
  mutate(anchor2_barcode = colnames(atac_nt@assays$GA@counts)[nt_anchor_mat$cell2]) %>% 
  dplyr::filter(anchor1_barcode %in% naive_nt_elcs)
write_tsv(nt_anchor_mat, glue("{result_folder}nt_anchor_matrix.tsv"))

nt_predicted.labels = TransferData(
  anchorset = nt_anchors,
  refdata = naive_nt_rna$cluster_EML,
  weight.reduction = atac_nt[['lsi']],
  dims = 2:30
)

atac_nt = AddMetaData(object = atac_nt, metadata = nt_predicted.labels)
nt_integrated_meta = atac_nt@meta.data

# new UMAPs upon integration
umap_rna = DimPlot(
  object = naive_nt_rna,
  group.by = "cluster_EML",
  label = TRUE,
  repel = TRUE) + 
  xlim(-12, 12) +
  ylim(-12, 12) +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("naive scRNA-Seq (non-treated)") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) + NoLegend()
umap_rna

atac_nt@meta.data$predicted.id = 
  factor(atac_nt@meta.data$predicted.id, levels = c("ELC", "MeLC", "TLC", "HLC"))
umap_atac = DimPlot(
  object = atac_nt,
  group.by = "predicted.id",
  label = TRUE,
  repel = TRUE) + 
  xlim(-12, 12) +
  ylim(-12, 12) +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("non-trt scATAC-Seq") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )
umap_atac

umaps_int = ggarrange(umap_rna, umap_atac, ncol = 2)
umaps_int

ggsave(
  glue("{result_folder}CCA_integration_UMAPs-with_ntATAC.png"),
  plot = umaps_int,
  width = 12,
  height = 6,
  dpi = 300,
)

ggsave(
  glue("{result_folder}CCA_integration_UMAPs-with_ntATAC.pdf"),
  plot = umaps_int,
  width = 12,
  height = 6,
  device = "pdf"
)

# co-embedding of scRNA-Seq and scATAC-Seq
nt_var_genes = VariableFeatures(naive_nt_rna)
nt_refdata = GetAssayData(naive_nt_rna, assay = "RNA", layer = "data")[nt_var_genes, ]

nt_imputation = TransferData(anchorset = nt_anchors, refdata = nt_refdata, 
                             weight.reduction = atac_nt[["lsi"]],
                           dims = 2:30)
atac_nt[["RNA"]] = nt_imputation

nt_coembed = merge(x = naive_nt_rna, y = atac_nt)

# make scaling and dim reduction on the co-embedded dataset
nt_coembed = ScaleData(nt_coembed, features = nt_var_genes, do.scale = FALSE)
nt_coembed = RunPCA(nt_coembed, features = nt_var_genes, verbose = FALSE)
nt_coembed = RunUMAP(nt_coembed, dims = 1:30)

nt_coembed_meta = nt_coembed@meta.data
nt_coembed_meta = nt_coembed_meta %>% 
  mutate(orig.ident = str_replace(orig.ident, "SeuratProject", "nt scATAC-Seq")) %>% 
  mutate(orig.ident = str_replace(orig.ident, "EZH2i", "naive nt scRNA-Seq"))
nt_coembed@meta.data = nt_coembed_meta

# export Seurat object
saveRDS(nt_coembed, file = glue("{result_folder}nt_scATAC_scRNA_coembed.Rds"))

# UMAPs
nt_coembed_umap_ident = DimPlot(nt_coembed, group.by = "orig.ident") +
  xlim(-12, 12) +
  ylim(-12, 12) +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("coembedding - group by experiment (non-trt)") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 10, hjust = 0.1),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )  
nt_coembed_umap_ident

nt_coembed_umap_eml = DimPlot(nt_coembed, group.by = "cluster_EML") +
  xlim(-12, 12) +
  ylim(-12, 12) +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("coembedding - group by EML (non-trt)") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 10, hjust = 0.1),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )  
nt_coembed_umap_eml

nt_umaps_coembed = ggarrange(nt_coembed_umap_ident, nt_coembed_umap_eml)
nt_umaps_coembed

ggsave(
  glue("{result_folder}nt_coembed_UMAPs.png"),
  plot = nt_umaps_coembed,
  width = 12,
  height = 6,
  dpi = 300,
)

ggsave(
  glue("{result_folder}nt_coembed_UMAPs.pdf"),
  plot = nt_umaps_coembed,
  width = 12,
  height = 6,
  device = "pdf"
)

### EZH2i 7d workflow ###
rm(atac_nt)
rm(nt_coembed)
## integration by CCA
# finding anchors
trt_anchors = FindTransferAnchors(
  reference = naive_trt_rna,
  query = atac_trt,
  features = VariableFeatures(object = naive_trt_rna),
  reference.assay = "RNA",
  query.assay = "GA",
  reduction = 'cca'
)

# extract anchor matrix of AnchorSet object
trt_anchor_mat = as_tibble(trt_anchors@anchors)
trt_anchor_mat = trt_anchor_mat %>%
  mutate(anchor1_barcode = naive_trt_rna$cell[trt_anchor_mat$cell1]) %>%
  mutate(anchor2_barcode = colnames(atac_trt@assays$GA@counts)[trt_anchor_mat$cell2]) %>% 
  dplyr::filter(anchor1_barcode %in% naive_trt_elcs)
write_tsv(trt_anchor_mat, glue("{result_folder}trt_anchor_matrix.tsv"))

trt_predicted.labels = TransferData(
  anchorset = trt_anchors,
  refdata = naive_trt_rna$cluster_EML,
  weight.reduction = atac_trt[['lsi']],
  dims = 2:30
)

atac_trt = AddMetaData(object = atac_trt, metadata = trt_predicted.labels)
trt_integrated_meta = atac_trt@meta.data

# export Seurat object
#saveRDS(atac_trt, file = glue("{result_folder}trt_scATAC-labeled.Rds"))

# new UMAPs upon integration
umap_rna = DimPlot(
  object = naive_trt_rna,
  group.by = "cluster_EML",
  label = TRUE,
  repel = TRUE) + 
  xlim(-12, 12) +
  ylim(-12, 12) +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("naive scRNA-Seq (EZH2i 7D)") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  ) + NoLegend()
umap_rna

atac_trt@meta.data$predicted.id = 
  factor(atac_trt@meta.data$predicted.id, levels = c("ELC", "MeLC", "TLC", "AMLC"))
trt_umap_atac = DimPlot(
  object = atac_trt,
  group.by = "predicted.id",
  label = TRUE,
  repel = TRUE) + 
  xlim(-12, 12) +
  ylim(-12, 12) +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("EZH2i 7d scATAC-Seq") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )  
trt_umap_atac

trt_umaps_int = ggarrange(umap_rna, trt_umap_atac)
trt_umaps_int

ggsave(
  glue("{result_folder}CCA_integration_UMAPs-with_trtATAC.png"),
  plot = trt_umaps_int,
  width = 12,
  height = 6,
  dpi = 300,
)

ggsave(
  glue("{result_folder}CCA_integration_UMAPs-with_trtATAC.pdf"),
  plot = trt_umaps_int,
  width = 12,
  height = 6,
  device = "pdf"
)

# co-embedding of scRNA-Seq and scATAC-Seq
trt_var_genes = VariableFeatures(naive_trt_rna)
trt_refdata = GetAssayData(naive_trt_rna, assay = "RNA", layer = "data")[trt_var_genes, ]

trt_imputation = TransferData(anchorset = trt_anchors, refdata = trt_refdata, 
                              weight.reduction = atac_trt[["lsi"]],
                              dims = 2:30)
atac_trt[["RNA"]] = trt_imputation

trt_coembed = merge(x = naive_trt_rna, y = atac_trt)

# make scaling and dim reduction on the co-embedded dataset
trt_coembed = ScaleData(trt_coembed, features = trt_var_genes, do.scale = FALSE)
trt_coembed = RunPCA(trt_coembed, features = trt_var_genes, verbose = FALSE)
trt_coembed = RunUMAP(trt_coembed, dims = 1:30)

trt_coembed_meta = trt_coembed@meta.data
trt_coembed_meta = trt_coembed_meta %>% 
  mutate(orig.ident = str_replace(orig.ident, "EZH2i", "naive trt scRNA-Seq")) %>% 
  mutate(orig.ident = str_replace(orig.ident, "SeuratProject", "EZH2i 7D scATAC-Seq")) 
trt_coembed@meta.data = trt_coembed_meta

# export Seurat object
saveRDS(trt_coembed, file = glue("{result_folder}trt_scATAC_scRNA_coembed.Rds"))

# UMAPs
trt_coembed_umap_ident = DimPlot(trt_coembed, group.by = "orig.ident") +
  xlim(-12, 12) +
  ylim(-12, 12) +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("coembedding - group by experiment (EZH2i 7D scATAC-Seq)") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 10, hjust = 0.1),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )  
trt_coembed_umap_ident

trt_coembed_umap_eml = DimPlot(trt_coembed, group.by = "cluster_EML") +
  xlim(-12, 12) +
  ylim(-12, 12) +
  scale_colour_brewer(palette = "Set1") +
  ggtitle("coembedding - group by EML (EZH2i 7D scATAC-Seq)") +
  theme(
    text = element_text(size = 25),
    plot.title = element_text(size = 10, hjust = 0.1),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(size = 25, color = "black"),
    axis.text.y = element_text(size = 25, color = "black")
  )  
trt_coembed_umap_eml

trt_umaps_coembed = ggarrange(trt_coembed_umap_ident, trt_coembed_umap_eml)
trt_umaps_coembed

ggsave(
  glue("{result_folder}trt_coembed_UMAPs.png"),
  plot = trt_umaps_coembed,
  width = 12,
  height = 6,
  dpi = 300,
)

ggsave(
  glue("{result_folder}trt_coembed_UMAPs.pdf"),
  plot = trt_umaps_coembed,
  width = 12,
  height = 6,
  device = "pdf"
)
