if (!require("pacman"))
  install.packages("pacman")
pacman::p_load("tidyverse",
               "Seurat",
               "Signac"
)

# export folder
result_folder = "../results/scRNA-Seq/"

# Cis-Bp TFs
cisbp = fread("../data/GRN/TFs/CIS-BP-TF_Information.txt", fill = TRUE)
cisbp = unique(cisbp$TF_Name)
cisbp = tibble(CIS_BP = cisbp) %>% pull(CIS_BP)

# Jaspar2024 TFs
jaspar = fread("../data/GRN/TFs/Jaspar2024_collection.txt",
               fill = TRUE,
               header = FALSE)
jaspar = jaspar$V2
jaspar = tibble(JASPAR_2024 = jaspar[which(jaspar != "")])
weirds = jaspar %>% filter(str_detect(JASPAR_2024, "\\::")) %>%
  separate(JASPAR_2024, sep = "::", into = c("V1", "V2"))

jaspar = jaspar %>% filter(!str_detect(JASPAR_2024, "\\::")) %>%
  filter(!str_detect(JASPAR_2024, "_")) %>%
  distinct_all() %>% pull(JASPAR_2024)

jaspar = unique(c(jaspar, weirds$V1, weirds$V2))
jaspar = tibble(JASPAR_2024 = jaspar) %>%
  filter(str_detect(JASPAR_2024, "^[[:upper:][:space:]]+$")) %>%
  pull(JASPAR_2024)

# union
tfs = unique(c(jaspar, cisbp))
tfs = tibble(CISBP_JASPAR2024_TFs = tfs)


# process and normalize hESC scRNA-Seq data
counts = readRDS("../data/scRNA-Seq/Nerges.counts.filter.rds")
print(paste0("Total number of cells: ", length(unique(colnames(
  counts
)))))

# keeping non-treated and D7 EZH2i treated cells
meta = readRDS("../data/scRNA-Seq/Nerges.meta.filter.rds")
meta = meta %>% dplyr::filter(cluster_EML != "Undef")
meta = meta %>% dplyr::filter(cellType == "EZH2i_Naive_WT" |
                                cellType == "EZH2i_Naive_D7")
counts = counts[,meta$cell]

seurat_rna = CreateSeuratObject(
  counts = counts,
  project = "scRNA_EZH2i",
  min.cells = 3,
  min.features = 200
)
seurat_rna@meta.data = cbind(seurat_rna@meta.data, meta)

seurat_rna = NormalizeData(seurat_rna, normalization.method = "LogNormalize", scale.factor = 10000)

# intersect with variable genes
seurat_rna = FindVariableFeatures(seurat_rna, selection.method = "vst", nfeatures = 2000)
tf_candidates = tibble(TF_name = intersect(VariableFeatures(seurat_rna), tfs$CISBP_JASPAR2024_TFs))

# write out for Pando
write_tsv(tf_candidates, "../data/GRN/Variable_TFs-Jaspar2024_CISBP.tsv")


