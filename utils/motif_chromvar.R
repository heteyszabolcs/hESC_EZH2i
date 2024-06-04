suppressPackageStartupMessages({
  library("Signac")
  library("Seurat")
  library("BSgenome.Hsapiens.UCSC.hg38")
  library("patchwork")
  library("tidyverse")
  library("glue")
  library("ggplot2")
  library("ggrastr")
  library("ggrepel")
})

# if returns sparse matrix error install it from source:
#devtools::install_github("ge11232002/TFBSTools")

# chromVar Nature paper: https://www.nature.com/articles/nmeth.4401

result_folder = "../results/motif_analysis/"

# read objects
coembed_seurat = readRDS("../results/scATAC-Seq/hESC_scATAC_nt_trt_coembedded.Rds")
coembed_motif = readRDS("../results/motif_analysis/trt_nt_coembed_scATAC-motif.Rds")
coembed_seurat[["motifs"]] = coembed_motif
DefaultAssay(object = coembed_seurat) = "peaks"

atac_trt_motif = readRDS(file = "../results/motif_analysis/trt_scATAC-motif.Rds")
atac_trt = readRDS(glue("../results/Seurat_integration/scATAC_annotation/trt_scATAC-labeled.Rds"))
atac_trt[["motifs"]] = atac_trt_motif
DefaultAssay(object = atac_trt) = "peaks"

# Cheng's scATAC annotation
chengs_annot = 
  read_tsv("../results/Seurat_integration/scATAC_annotation/Chen_scATAC_annot-extended.tsv")

chengs_annot_trt = chengs_annot %>% separate(cell, sep = ".10X.", into = c("a", "id")) %>% 
  dplyr::select(-a) %>% 
  dplyr::filter(devTime == "EZH2i")

# EZH2i treated TLCs, MeLCs
tlc = chengs_annot_trt %>% dplyr::filter(str_detect(cluster_EML, "TLC")) %>% pull(id)
melc = chengs_annot_trt %>% dplyr::filter(str_detect(cluster_EML, "MeLC")) %>% pull(id)

meta = atac_trt@meta.data %>% rownames_to_column(var = "cell_id") %>% 
  mutate(predicted.id_cheng = case_when(cell_id %in% tlc ~ "TLC",
                                        cell_id %in% melc ~ "MeLC",
                                        .default = predicted.id))
rownames(meta) = meta$cell_id
meta = meta %>% dplyr::select(-cell_id)
atac_trt@meta.data = meta

# ELCs
chengs_elc = chengs_annot %>% dplyr::filter(EML == "merged_ELC") %>% 
  mutate(cell = unname(sapply(cell, function(x) { strsplit(x, ".10X.")[[1]][2]})))

coembed_meta = coembed_seurat@meta.data
coembed_meta = coembed_meta %>% mutate(cell_id = 
                         unname(sapply(rownames(coembed_meta), function(x) { strsplit(x, "_")[[1]][1]})))
coembed_seurat@meta.data = coembed_meta
coembed_cheng = subset(coembed_seurat, subset = cell_id %in% chengs_elc$cell)

## ChromVar analysis
# keep only valid GRanges strings
valid_feats = which(seqnames(atac_trt[['peaks']]@ranges) %in% 
                      standardChromosomes(atac_trt[['peaks']]@ranges))
atac_trt = atac_trt[valid_feats, ]
atac_trt = RunChromVAR(
  object = atac_trt,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay = "motifs"
)

valid_feats = which(seqnames(coembed_seurat[['peaks']]@ranges) %in% standardChromosomes(coembed_seurat[['peaks']]@ranges))
coembed_seurat = coembed_seurat[valid_feats, ]
coembed_seurat = RunChromVAR(
  object = coembed_seurat,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay = "motifs"
)

# TF activity marker analysis
# ELCs
DefaultAssay(coembed_seurat) = 'chromvar'
trt_vs_nt_ELC_diff_activity = FindMarkers(
  object = coembed_seurat,
  ident.1 = "EZH2i_7D_scATAC_Seq",
  ident.2 = "non_trt_scATAC_Seq",
  group.by = "orig.ident",
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

trt_vs_nt_ELC_diff_activity %>%
  mutate(gene_name = ConvertMotifID(
    coembed_seurat,
    assay = "motifs",
    id = rownames(trt_vs_nt_ELC_diff_activity)
  )) %>% write_tsv(., 
          glue("{result_folder}ELC_trt_vs_nt-ChromVar_output.tsv"))

trt_vs_nt_ELC_top_activities = trt_vs_nt_ELC_diff_activity %>% 
  arrange(desc(avg_diff)) %>% 
  rownames %>% 
  head(10)

sign_in_trt = trt_vs_nt_ELC_diff_activity %>% 
  dplyr::filter(p_val_adj < 0.05) %>% dplyr::filter(avg_diff > 1.5) %>% rownames

# ChromVar on Cheng's ELC population
valid_feats = which(seqnames(coembed_cheng[['peaks']]@ranges) %in% standardChromosomes(coembed_cheng[['peaks']]@ranges))
coembed_cheng = coembed_cheng[valid_feats, ]
coembed_cheng = RunChromVAR(
  object = coembed_cheng,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  assay = "motifs"
)

# TF activity marker analysis
DefaultAssay(coembed_cheng) = 'chromvar'
trt_vs_nt_ELC_diff_activity_cheng = FindMarkers(
  object = coembed_cheng,
  ident.1 = "EZH2i_7D_scATAC_Seq",
  ident.2 = "non_trt_scATAC_Seq",
  group.by = "orig.ident",
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

trt_vs_nt_ELC_diff_activity_cheng = trt_vs_nt_ELC_diff_activity_cheng %>%
  mutate(gene_name = ConvertMotifID(
    coembed_seurat,
    assay = "motifs",
    id = rownames(trt_vs_nt_ELC_diff_activity_cheng)
  )) 

trt_vs_nt_ELC_diff_activity = read_tsv(glue("{result_folder}ELC_trt_vs_nt-ChromVar_output.tsv"))

length(intersect(trt_vs_nt_ELC_diff_activity_cheng %>% 
                   dplyr::filter(p_val_adj < 0.05) %>% 
                   arrange(desc(avg_diff)) %>% top_n(n = 30, wt = avg_diff) %>% 
                   pull(gene_name),
                 trt_vs_nt_ELC_diff_activity %>% 
                   dplyr::filter(p_val_adj < 0.05) %>% 
                   arrange(desc(avg_diff)) %>% top_n(n = 30, wt = avg_diff) %>% 
                   pull(gene_name)))
setdiff(
  trt_vs_nt_ELC_diff_activity_cheng %>% 
    dplyr::filter(p_val_adj < 0.05) %>% 
    arrange(desc(avg_diff)) %>% top_n(n = 30, wt = avg_diff) %>%
    pull(gene_name),
  trt_vs_nt_ELC_diff_activity %>% 
    dplyr::filter(p_val_adj < 0.05) %>% 
    arrange(desc(avg_diff)) %>% top_n(n = 30, wt = avg_diff) %>%
    pull(gene_name)
)

# EZH2i TLC and MeLC specific activities
DefaultAssay(atac_trt) = 'chromvar'
atac_trt@meta.data = meta
trt_TLC_vs_all = FindMarkers(
  object = atac_trt,
  ident.1 = "TLC",
  ident.2 = NULL,
  group.by = "predicted.id_cheng",
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

trt_TLC_vs_all %>%
  mutate(gene_name = ConvertMotifID(
    atac_trt,
    assay = "motifs",
    id = rownames(trt_TLC_vs_all)
  )) %>% write_tsv(., 
                   glue("{result_folder}trt_TLC_vs_all.tsv"))

trt_TLC_vs_all_top = trt_TLC_vs_all %>% 
  arrange(desc(avg_diff)) %>% 
  rownames %>% 
  head(10)


trt_MeLC_vs_all = FindMarkers(
  object = atac_trt,
  ident.1 = "MeLC",
  ident.2 = NULL,
  group.by = "predicted.id_cheng",
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

trt_MeLC_vs_all %>%
  mutate(gene_name = ConvertMotifID(
    atac_trt,
    assay = "motifs",
    id = rownames(trt_MeLC_vs_all)
  )) %>% write_tsv(., 
                   glue("{result_folder}trt_MeLC_vs_all.tsv"))

trt_MeLC_vs_all_top = trt_MeLC_vs_all %>% 
  arrange(desc(avg_diff)) %>% 
  rownames %>% 
  head(10)

# motif logos
require("ggseqlogo")
pdf(
  file = "../results/motif_analysis/trt_TLC_vs_all.pdf",
  width = 10,
  height = 4
)
MotifPlot(
  object = atac_trt,
  motifs = trt_TLC_vs_all_top,
  assay = "motifs"
)
dev.off()

png(
  file = "../results/motif_analysis/trt_TLC_vs_all.png",
  width = 13,
  height = 10,
  units = "cm",
  res = 300
)
MotifPlot(
  object = atac_trt,
  motifs = trt_TLC_vs_all_top,
  assay = "motifs"
)
dev.off()


pdf(
  file = "../results/motif_analysis/trt_MeLC_vs_all.pdf",
  width = 10,
  height = 4
)
MotifPlot(
  object = atac_trt,
  motifs = trt_MeLC_vs_all_top,
  assay = "motifs"
)
dev.off()

png(
  file = "../results/motif_analysis/trt_MeLC_vs_all.png",
  width = 13,
  height = 10,
  units = "cm",
  res = 300
)
MotifPlot(
  object = atac_trt,
  motifs = trt_MeLC_vs_all_top,
  assay = "motifs"
)
dev.off()



pdf(
  file = "../results/motif_analysis/trt_ELC_vs_nt_ELC-top_chromVar_act.pdf",
  width = 10,
  height = 4
)
MotifPlot(
  object = coembed_seurat,
  motifs = trt_vs_nt_ELC_top_activities,
  assay = "motifs"
)
dev.off()

pdf(
  file = "../results/motif_analysis/trt_ELC_vs_nt_ELC-1.5fc_adjp0.05_chromVar_act.pdf",
  width = 14,
  height = 14
)
MotifPlot(
  object = coembed_seurat,
  motifs = sign_in_trt,
  assay = "motifs"
)
dev.off()

MotifPlot(
  object = atac_trt,
  motifs = trt_MeLC_vs_all,
  assay = "motifs"
)

# other visualizations
# scatter
volc_input = trt_vs_nt_ELC_diff_activity %>% 
  rownames_to_column(var = "motif_id") %>% 
  mutate(group = case_when(
  avg_diff > 1.5 & p_val_adj < 0.05 ~ "up",
  avg_diff < -1.5 & p_val_adj < 0.05 ~ "down", .default = "unaltered"
)) %>%
  mutate(sign_label = case_when(group == "up" ~ motif_id, 
                                group == "down" ~ motif_id, .default = NA))
labels = volc_input %>% pull(sign_label)
labels = ConvertMotifID(coembed_seurat, assay = "motifs", id = labels)

# add color, size, alpha indications
cols = c("up" = "#fc9272", "down" = "#3182bd", "unaltered" = "grey")
sizes = c("up" = 4, "down" = 4, "unaltered" = 2)
alphas = c("up" = 1, "down" = 1, "unaltered" = 0.5)

# ChromVar motif volcano
chromvar_volc = volc_input %>%
  ggplot(aes(x = avg_diff,
             y = -log10(p_val_adj),
             fill = group,    
             size = group,
             alpha = group)) +
  ggrastr::geom_point_rast(shape = 21, # Specify shape and colour as fixed local parameters    
                           colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = 1.5,
             linetype = "dashed") +
  geom_vline(xintercept = -1.5,
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-3, 3, 1)),  	 
                     limits = c(-3, 3)) +
  scale_y_continuous(breaks = c(seq(0, 200, 25)),  	 
                     limits = c(0, 200)) +
  labs(
    title = "ChromVar motif enrichment",
    subtitle = "EZH2i 7D vs. non-treated ELCs",
    x = "motif avarage difference",
    y = "-log10 adj.p-value",
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
  geom_text_repel(label = labels, size = 4, max.overlaps = 50, min.segment.length = 0.1) # add labels
chromvar_volc

ggsave(
  glue("{result_folder}ELC_ChromVar_volc.png"),
  plot = chromvar_volc,
  width = 7,
  height = 5,
  dpi = 300,
)

ggsave(
  glue("{result_folder}ELC_ChromVar_volc.pdf"),
  plot = chromvar_volc,
  width = 7,
  height = 5,
  device = "pdf"
)

volc_input = trt_vs_nt_ELC_diff_activity_cheng %>% 
  rownames_to_column(var = "motif_id") %>% 
  mutate(group = case_when(
    avg_diff > 1.5 & p_val_adj < 0.05 ~ "up",
    avg_diff < -1.5 & p_val_adj < 0.05 ~ "down", .default = "unaltered"
  )) %>%
  mutate(sign_label = case_when(group == "up" ~ motif_id, 
                                group == "down" ~ motif_id, .default = NA))
labels = volc_input %>% pull(sign_label)
labels = ConvertMotifID(coembed_cheng, assay = "motifs", id = labels)

# add color, size, alpha indications
cols = c("up" = "#fc9272", "down" = "#3182bd", "unaltered" = "grey")
sizes = c("up" = 4, "down" = 4, "unaltered" = 2)
alphas = c("up" = 1, "down" = 1, "unaltered" = 0.5)

# ChromVar motif volcano
chromvar_volc_cheng = volc_input %>%
  ggplot(aes(x = avg_diff,
             y = -log10(p_val_adj),
             fill = group,    
             size = group,
             alpha = group)) +
  ggrastr::geom_point_rast(shape = 21, # Specify shape and colour as fixed local parameters    
                           colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = 1.5,
             linetype = "dashed") +
  geom_vline(xintercept = -1.5,
             linetype = "dashed") +
  scale_fill_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-3, 3, 1)),  	 
                     limits = c(-3, 3)) +
  scale_y_continuous(breaks = c(seq(0, 200, 25)),  	 
                     limits = c(0, 200)) +
  labs(
    title = "ChromVar motif enrichment",
    subtitle = "EZH2i 7D vs. non-treated ELCs (Cheng's annotation)",
    x = "motif avarage difference",
    y = "-log10 adj.p-value",
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
  geom_text_repel(label = labels, size = 4, max.overlaps = 50, min.segment.length = 0.1) # add labels
chromvar_volc_cheng

ggsave(
  glue("{result_folder}ELC_ChromVar_ChengELC-volc.png"),
  plot = chromvar_volc_cheng,
  width = 7,
  height = 5,
  dpi = 300,
)

ggsave(
  glue("{result_folder}ELC_ChromVar_ChengELC-volc.pdf"),
  plot = chromvar_volc_cheng,
  width = 7,
  height = 5,
  device = "pdf"
)