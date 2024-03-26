if (!require("pacman"))
  install.packages("pacman")
pacman::p_load("Seurat",
               "SeuratData",
               "SeuratDisk",
               "Signac",
               "glue",
               "tidyverse",
               "GenomicFeatures",
               "GenomeInfoDb",
               "ggpubr",
               "EnsDb.Hsapiens.v86",
               "ggrepel"
)

# export folder
result_folder = "../results/scRNA-Seq/"

# Seurat object (e.g. scRNA-Seq data)
seurat_rna = readRDS(file = glue("{result_folder}hESC_EZH2i_scRNA_Seq.Rds"))

# subset
trt = subset(seurat_rna, subset = batch == 'Nerges_EZH2i_Naive_D7')

# converting to old Assay format
trt[["RNA3"]] = as(object = trt[["RNA"]], Class = "Assay")
write.csv(trt@assays$RNA3@counts, glue("{result_folder}hESC_EZH2i_scRNA_Seq-EZH2i_7D.csv"))
DefaultAssay(trt) = "RNA3"
trt[["RNA"]] = NULL
trt = RenameAssays(object = trt, RNA3 = 'RNA')

# save and convert (SeuratDisk package)
SaveH5Seurat(trt, filename = glue("{result_folder}hESC_EZH2i_scRNA_Seq-EZH2i_7D.h5Seurat"))
Convert(glue("{result_folder}hESC_EZH2i_scRNA_Seq-EZH2i_7D.h5Seurat"), dest = "h5ad")

# subset
nt = subset(seurat_rna, subset = batch == 'Nerges_EZH2i_Naive')

# converting to old Assay format
nt[["RNA3"]] = as(object = nt[["RNA"]], Class = "Assay")
write.csv(nt@assays$RNA3@counts, glue("{result_folder}hESC_EZH2i_scRNA_Seq-nontrt.csv"))
DefaultAssay(nt) = "RNA3"
nt[["RNA"]] = NULL
nt = RenameAssays(object = nt, RNA3 = 'RNA')

# save and convert (SeuratDisk package)
SaveH5Seurat(nt, filename = glue("{result_folder}hESC_EZH2i_scRNA_Seq-nontrt.h5Seurat"))
Convert(glue("{result_folder}hESC_EZH2i_scRNA_Seq-nontrt.h5Seurat"), dest = "h5ad")

