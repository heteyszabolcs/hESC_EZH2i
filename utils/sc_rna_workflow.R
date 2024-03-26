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

#options(Seurat.object.assay.version = "v3")

# export folder
result_folder = "../results/scRNA-Seq/"

# hESC scRNA-Seq data
counts = readRDS("../data/scRNA-Seq/Nerges.counts.filter.rds")
print(paste0("Total number of cells: ", length(unique(colnames(counts)))))

# keeping ELC cells, nontreated and D7 EZH2i treated
meta = readRDS("../data/scRNA-Seq/Nerges.meta.filter.rds")
meta = meta %>% dplyr::filter(cluster_EML != "Undef")
meta = meta %>% dplyr::filter(cellType == "EZH2i_Naive_WT" | cellType == "EZH2i_Naive_D7")
counts = counts[,meta$cell]

# explore cell annotations
nt_cell_types = meta %>% dplyr::filter(cellType == "EZH2i_Naive_WT") %>% 
  group_by(cluster_EML) %>% count() %>% mutate(sample = "non-treated naive")
trt_cell_types = meta %>% dplyr::filter(cellType == "EZH2i_Naive_D7") %>% 
  group_by(cluster_EML) %>% count()  %>% mutate(sample = "EZH2i 7D naive")
cell_types = rbind(nt_cell_types, trt_cell_types)

order = factor(cell_types$cluster_EML, levels = c("ELC", "MeLC", "TLC", "AMLC", "HLC"))
order2 = factor(cell_types$sample, levels = c("non-treated naive", "EZH2i 7D naive"))
ggplot(cell_types, aes(x = order, y = n, fill = order2)) +
  geom_bar(stat="identity", position=position_dodge(), color = "black") +
  scale_fill_manual(values = c("#bdbdbd", "#fec44f")) +
  ylim(0, 1200) +
  labs(
    title = "scRNA-Seq annotations - naive subset",
    x = "",
    y = "cell count",
    fill = ""
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 9),
    plot.title = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black")
  ) 

ggsave(
  glue("{result_folder}scRNA_Seq_naive_annotations-bar.pdf"),
  plot = last_plot(),
  width = 5,
  height = 4
)

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

# export in h5ad format (at first we must set back the Assay format)
seurat_rna[["RNA3"]] = as(object = seurat_rna[["RNA"]], Class = "Assay")
DefaultAssay(seurat_rna) = "RNA3"
seurat_rna[["RNA"]] = NULL
seurat_rna = RenameAssays(object = seurat_rna, RNA3 = 'RNA')

SaveH5Seurat(seurat_rna, filename = glue("{result_folder}hESC_EZH2i_scRNA_Seq.h5Seurat"))
Convert(glue("{result_folder}hESC_EZH2i_scRNA_Seq.h5Seurat"), dest = "h5ad")

# marker analysis
tlc_vs_elc = FindMarkers(seurat_rna,
                         ident.1 = "TLC",
                         ident.2 = "ELC",
                         group.by = "cluster_EML")
tlc_vs_elc_sign = tlc_vs_elc %>% dplyr::filter(abs(avg_log2FC) > 2 & p_val_adj < 0.05)

melc_vs_elc = FindMarkers(seurat_rna,
                         ident.1 = "MeLC",
                         ident.2 = "ELC",
                         group.by = "cluster_EML",
                        )
melc_vs_elc_sign = melc_vs_elc %>% dplyr::filter(abs(avg_log2FC) > 2 & p_val_adj < 0.05)

# marker genes between EZH2i treated ELCs and non-treated ELCs
elc_seurat_rna = subset(x = seurat_rna, subset = cluster_EML == "ELC")
trt_elc_vs_nt_elc = FindMarkers(elc_seurat_rna,
                          ident.1 = "EZH2i_Naive_D7_10X",
                          ident.2 = "EZH2i_Naive_WT_10X",
                          group.by = "SID",
                          logfc.threshold = 0.2)
trt_elc_vs_nt_elc$gene = rownames(trt_elc_vs_nt_elc)
write_tsv(trt_elc_vs_nt_elc, glue("{result_folder}scRNASeq-trt_ELC_vs_nt_ELC-fc_table.txt"))

trt_elc_vs_nt_elc_sign = trt_elc_vs_nt_elc %>% 
  dplyr::filter(abs(avg_log2FC) > 2 & p_val_adj < 0.05) 
trt_elc_vs_nt_elc_sign_up = trt_elc_vs_nt_elc %>% 
  dplyr::filter(avg_log2FC > 2 & p_val_adj < 0.05) %>% dplyr::select(gene)
write_tsv(trt_elc_vs_nt_elc_sign_up, "../data/GRN/scRNASeq-trt_ELC_vs_nt_ELC-sign_up_genes.txt")

FeaturePlot(
  elc_seurat_rna,
  features = trt_elc_vs_nt_elc_sign %>% 
    arrange(desc(avg_log2FC)) %>% 
    top_n(6, wt = avg_log2FC) %>% pull(gene),
  dims = c(1, 2), cols = c("#ece7f2", "#5e1310")) 

ggsave(
  glue("{result_folder}EZH2i_7D_ELC_vs_nt_ELC-top6_fc_featurePlots.png"),
  plot = last_plot(),
  width = 8,
  height = 8,
  dpi = 300,
)

ggsave(
  glue("{result_folder}EZH2i_7D_ELC_vs_nt_ELC-top6_fc_featurePlots.pdf"),
  plot = last_plot(),
  width = 8,
  height = 8,
  device = "pdf"
)

FeaturePlot(
  elc_seurat_rna,
  features = trt_elc_vs_nt_elc_sign %>% 
    dplyr::filter(avg_log2FC > 0) %>% 
    arrange(p_val_adj) %>% 
    top_n(-6, wt = p_val_adj) %>% pull(gene),
  dims = c(1, 2), cols = c("#ece7f2", "#de2d26"))

ggsave(
  glue("{result_folder}EZH2i_7D_ELC_vs_nt_ELC-top6_adjp_featurePlots.png"),
  plot = last_plot(),
  width = 8,
  height = 8,
  dpi = 300,
)

ggsave(
  glue("{result_folder}EZH2i_7D_ELC_vs_nt_ELC-top6_adjp_featurePlots.pdf"),
  plot = last_plot(),
  width = 8,
  height = 8,
  device = "pdf"
)

# these TFs are coming from Pando GRN protocol
# ChromVar/TOBIAS inputs
interesting_up_tfs = c("CREB3", "HMGXB4", "JUNB", "ATF7", "HIC2", "CREB1", "FOSL2", "TFAP2C")
interesting_down_tfs = c("SALL4", "ZNF589", "CREM", "TCF3", "GABPA", "MAFG", "FOSL2", "RUNX1", 
                         "ERF", "NFYA", "ZNF681")
trt_elc_vs_nt_elc %>% dplyr::filter(gene %in% interesting_up_tfs)
trt_elc_vs_nt_elc %>% dplyr::filter(gene %in% interesting_down_tfs)
# FindVariables top 100 inputs
interesting_up_tfs2 = c("CUX1", "HOXA10", "ZNF891", "KLF9", "SP5", "ZNF33B")
interesting_down_tfs2 = c("ZNF503", "SALL3", "ZNF22", "SP2", "RARB", "MEIS2")
trt_elc_vs_nt_elc %>% dplyr::filter(gene %in% interesting_up_tfs2)
trt_elc_vs_nt_elc %>% dplyr::filter(gene %in% interesting_down_tfs2)

# ChromVar / TOBIAS outputs
chromvar_tob = read_tsv("../data/GRN/ChromVar_TOBIAS_trt_ELC_vs_nt_ELC-sign_genes.txt") %>% 
  pull(gene_name)

chromvar_tob_fcs = trt_elc_vs_nt_elc %>% dplyr::filter(gene %in% chromvar_tob)

# visualizations
volc_input = trt_elc_vs_nt_elc %>% 
  mutate(group = case_when(
    avg_log2FC > 2 & p_val_adj < 1e-25 ~ "up",
    avg_log2FC < -2 & p_val_adj < 1e-25 ~ "down", .default = "unaltered"
  )) %>%
  mutate(sign_label = case_when(group == "up" ~ gene, 
                                group == "down" ~ gene, .default = NA))
labels = volc_input %>% pull(sign_label)

# add color, size, alpha indications
cols = c("up" = "#fc9272", "down" = "#3182bd", "unaltered" = "grey")
sizes = c("up" = 4, "down" = 4, "unaltered" = 2)
alphas = c("up" = 1, "down" = 1, "unaltered" = 0.5)

# EZH2i 7D ELC volcano
trt_elc_vs_nt_elc_volc = volc_input %>%
  ggplot(aes(x = avg_log2FC,
             y = -log10(p_val_adj),
             fill = group,    
             size = group,
             alpha = group)) +
  ggrastr::geom_point_rast(shape = 21, # Specify shape and colour as fixed local parameters    
                           colour = "black") +
  geom_hline(yintercept = -log10(1e-25),
             linetype = "dashed") +
  geom_vline(xintercept = 2,
             linetype = "dashed") +
  geom_vline(xintercept = -2,
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-4, 4, 1)),  	 
                     limits = c(-4, 4)) +
  scale_y_continuous(breaks = c(seq(0, 200, 50)),  	 
                     limits = c(0, 200)) +
  labs(
    title = "EZH2i 7D ELCs vs non-treated ELCs",
    subtitle = "scRNA-Seq marker analysis",
    x = "log2FoldChange",
    y = "-log10 adj. p-value",
    fill = " "
  ) +
  guides(alpha = FALSE, size = FALSE, fill = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 6),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  geom_text_repel(label = labels, size = 4, max.overlaps = 110, min.segment.length = 0.1) # add labels
trt_elc_vs_nt_elc_volc

ggsave(
  glue("{result_folder}EZH2i_7D_ELC_vs_nt_ELC-volc.png"),
  plot = trt_elc_vs_nt_elc_volc,
  width = 10,
  height = 7,
  dpi = 300,
)

ggsave(
  glue("{result_folder}EZH2i_7D_ELC_vs_nt_ELC-volc.pdf"),
  plot = trt_elc_vs_nt_elc_volc,
  width = 10,
  height = 7,
  device = "pdf"
)

volc_input = trt_elc_vs_nt_elc %>% 
  mutate(group = case_when(
    gene %in% chromvar_tob ~ "enriched in ChromVar/TOBIAS",
    avg_log2FC > 2 & p_val_adj < 1e-25 ~ "up",
    avg_log2FC < -2 & p_val_adj < 1e-25 ~ "down",
    .default = "unaltered"
  )) %>%
  mutate(sign_label = case_when(group == "enriched in ChromVar/TOBIAS" ~ gene, 
                                .default = NA))

labels = volc_input %>% pull(sign_label)

# add color, size, alpha indications
cols = c("up" = "#fc9272", "down" = "#3182bd", "unaltered" = "white", 
         "enriched in ChromVar/TOBIAS" = "yellow")
sizes = c("up" = 4, "down" = 4, "unaltered" = 2, "enriched in ChromVar/TOBIAS" = 4)
alphas = c("up" = 1, "down" = 1, "unaltered" = 0, "enriched in ChromVar/TOBIAS" = 1)

trt_elc_vs_nt_elc_volc_2 = volc_input %>%
  ggplot(aes(x = avg_log2FC,
             y = -log10(p_val_adj),
             fill = group,    
             size = group,
             alpha = group)) +
  ggrastr::geom_point_rast(shape = 21, # Specify shape and colour as fixed local parameters    
                           colour = "black") +
  geom_hline(yintercept = -log10(1e-25),
             linetype = "dashed") +
  geom_vline(xintercept = 2,
             linetype = "dashed") +
  geom_vline(xintercept = -2,
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-4, 4, 1)),  	 
                     limits = c(-4, 4)) +
  scale_y_continuous(breaks = c(seq(0, 200, 50)),  	 
                     limits = c(0, 200)) +
  labs(
    title = "EZH2i 7D ELCs vs non-treated ELCs",
    subtitle = "scRNA-Seq marker analysis",
    x = "log2FoldChange",
    y = "-log10 adj. p-value",
    fill = " "
  ) +
  guides(alpha = FALSE, size = FALSE, fill = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 6),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  geom_text_repel(label = labels, size = 4, max.overlaps = 110, min.segment.length = 0.1) # add labels
trt_elc_vs_nt_elc_volc_2

ggsave(
  glue("{result_folder}EZH2i_7D_ELC_vs_nt_ELC_TOBIAS-volc.png"),
  plot = trt_elc_vs_nt_elc_volc_2,
  width = 10,
  height = 7,
  dpi = 300,
)

ggsave(
  glue("{result_folder}EZH2i_7D_ELC_vs_nt_ELC_TOBIAS-volc.pdf"),
  plot = trt_elc_vs_nt_elc_volc_2,
  width = 10,
  height = 7,
  device = "pdf"
)

volc_input = trt_elc_vs_nt_elc %>% 
  mutate(group = case_when(
    avg_log2FC > 2 & p_val_adj < 1e-25 ~ "up",
    avg_log2FC < -2 & p_val_adj < 1e-25 ~ "down", 
    gene %in% c(interesting_up_tfs, interesting_down_tfs) ~ "TFs by Pando (ChromVar/Tobias)",
    .default = "unaltered"
  )) %>%
  mutate(sign_label = case_when(group == "TFs by Pando (ChromVar/Tobias)" ~ gene, 
                                .default = NA))

labels = volc_input %>% pull(sign_label)

# add color, size, alpha indications
cols = c("up" = "#fc9272", "down" = "#3182bd", "unaltered" = "white", 
         "TFs by Pando (ChromVar/Tobias)" = "#bcbddc")
sizes = c("up" = 4, "down" = 4, "unaltered" = 2, "TFs by Pando (ChromVar/Tobias)" = 4)
alphas = c("up" = 1, "down" = 1, "unaltered" = 0, "TFs by Pando (ChromVar/Tobias)" = 1)

trt_elc_vs_nt_elc_volc_3 = volc_input %>%
  ggplot(aes(x = avg_log2FC,
             y = -log10(p_val_adj),
             fill = group,    
             size = group,
             alpha = group)) +
  ggrastr::geom_point_rast(shape = 21, # Specify shape and colour as fixed local parameters    
                           colour = "black") +
  geom_hline(yintercept = -log10(1e-25),
             linetype = "dashed") +
  geom_vline(xintercept = 2,
             linetype = "dashed") +
  geom_vline(xintercept = -2,
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-4, 4, 1)),  	 
                     limits = c(-4, 4)) +
  scale_y_continuous(breaks = c(seq(0, 200, 50)),  	 
                     limits = c(0, 200)) +
  labs(
    title = "EZH2i 7D ELCs vs non-treated ELCs",
    subtitle = "scRNA-Seq marker analysis",
    x = "log2FoldChange",
    y = "-log10 adj. p-value",
    fill = " "
  ) +
  guides(alpha = FALSE, size = FALSE, fill = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 6),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  geom_text_repel(label = labels, size = 4, max.overlaps = 110, min.segment.length = 0.1) # add labels
trt_elc_vs_nt_elc_volc_3

ggsave(
  glue("{result_folder}EZH2i_7D_ELC_vs_nt_ELC_Pando_Chrom-Tob-volc.png"),
  plot = trt_elc_vs_nt_elc_volc_3,
  width = 10,
  height = 7,
  dpi = 300,
)

ggsave(
  glue("{result_folder}EZH2i_7D_ELC_vs_nt_ELC_Pando_Chrom-Tob-volc.pdf"),
  plot = trt_elc_vs_nt_elc_volc_3,
  width = 10,
  height = 7,
  device = "pdf"
)

volc_input = trt_elc_vs_nt_elc %>% 
  mutate(group = case_when(
    avg_log2FC > 2 & p_val_adj < 1e-25 ~ "up",
    avg_log2FC < -2 & p_val_adj < 1e-25 ~ "down", 
    gene %in% c(interesting_up_tfs2, interesting_down_tfs2) ~ "TFs by Pando (FindVars100)",
    .default = "unaltered"
  )) %>%
  mutate(sign_label = case_when(group == "TFs by Pando (FindVars100)" ~ gene, 
                                .default = NA))

labels = volc_input %>% pull(sign_label)

# add color, size, alpha indications
cols = c("up" = "#fc9272", "down" = "#3182bd", "unaltered" = "white", 
         "TFs by Pando (FindVars100)" = "pink")
sizes = c("up" = 4, "down" = 4, "unaltered" = 2, "TFs by Pando (FindVars100)" = 4)
alphas = c("up" = 1, "down" = 1, "unaltered" = 0, "TFs by Pando (FindVars100)" = 1)

trt_elc_vs_nt_elc_volc_4 = volc_input %>%
  ggplot(aes(x = avg_log2FC,
             y = -log10(p_val_adj),
             fill = group,    
             size = group,
             alpha = group)) +
  ggrastr::geom_point_rast(shape = 21, # Specify shape and colour as fixed local parameters    
                           colour = "black") +
  geom_hline(yintercept = -log10(1e-25),
             linetype = "dashed") +
  geom_vline(xintercept = 2,
             linetype = "dashed") +
  geom_vline(xintercept = -2,
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-4, 4, 1)),  	 
                     limits = c(-4, 4)) +
  scale_y_continuous(breaks = c(seq(0, 200, 50)),  	 
                     limits = c(0, 200)) +
  labs(
    title = "EZH2i 7D ELCs vs non-treated ELCs",
    subtitle = "scRNA-Seq marker analysis",
    x = "log2FoldChange",
    y = "-log10 adj. p-value",
    fill = " "
  ) +
  guides(alpha = FALSE, size = FALSE, fill = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 6),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  geom_text_repel(label = labels, size = 4, max.overlaps = 110, min.segment.length = 0.1) # add labels
trt_elc_vs_nt_elc_volc_4

ggsave(
  glue("{result_folder}EZH2i_7D_ELC_vs_nt_ELC_Pando_FindVars100-volc.pdf"),
  plot = trt_elc_vs_nt_elc_volc_4,
  width = 10,
  height = 7,
  device = "pdf"
)

ggsave(
  glue("{result_folder}EZH2i_7D_ELC_vs_nt_ELC_Pando_FindVars100-volc.png"),
  plot = trt_elc_vs_nt_elc_volc_4,
  width = 10,
  height = 7,
  dpi = 300,
)

x =seurat_rna@meta.data
