if (!require("pacman"))
  install.packages("pacman")
pacman::p_load("glue",
               "tidyverse",
               "data.table",
               "ggrastr",
               "ggrepel",
               "ggpubr",
               "circlize",
               "ComplexHeatmap",
               "viridis")

# result folder
result_folder = "../results/TOBIAS/KO_experiments/"

# functions
make_volcanoes = function(ko) {
  print(ko)
  res = fread(
    glue(
      "../data/bulk_ATAC-Seq/KO_experiments/TOBIAS/{ko}/bindetect_results.txt"
    )
  )
  colnames(res) = str_to_upper(colnames(res))
  res = res %>%
    select(NAME,
           ends_with("SORTED_FOOTPRINT_PVALUE"),
           ends_with("CHANGE")) %>%
    select(NAME, starts_with(glue("ATAC_{ko}_HESC_KO_EZH2I"))) %>% 
    select(NAME, contains("ATAC_WT_HESC_EZH2I"))
  colnames(res) = c(
    "gene_name",
    "KO_EZH2i_vs_WT_EZH2i_pvalue",
    "KO_EZH2i_vs_WT_EZH2i_fp_change"
  )
  
  volc_input = res %>%
    mutate(
      group = case_when(
        KO_EZH2i_vs_WT_EZH2i_fp_change > 0.10 & KO_EZH2i_vs_WT_EZH2i_pvalue < 0.05 ~ "up",
        KO_EZH2i_vs_WT_EZH2i_fp_change < -0.10 &
          KO_EZH2i_vs_WT_EZH2i_pvalue < 0.05 ~ "down",
        .default = "unaltered"
      )
    ) %>%
    mutate(sign_label = case_when(group == "up" ~ gene_name, group == "down" ~ gene_name, .default = NA))
  labels = volc_input %>% pull(sign_label)
  
  # add color, size, alpha indications
  cols = c(
    "up" = "#fc9272",
    "down" = "#3182bd",
    "unaltered" = "grey"
  )
  sizes = c("up" = 4,
            "down" = 4,
            "unaltered" = 2)
  alphas = c("up" = 1,
             "down" = 1,
             "unaltered" = 1)
  
  volc_input = volc_input %>% 
    dplyr::arrange(desc(group))
  
  labels = factor(volc_input$sign_label, levels = unique(volc_input$sign_label))
  
  # KO EZH2i vs WT EZH2i
  volc = volc_input %>%
    ggplot(aes(
      x = KO_EZH2i_vs_WT_EZH2i_fp_change,
      y = -log10(KO_EZH2i_vs_WT_EZH2i_pvalue),
      fill = group,
      size = group,
      alpha = group
    )) +
    ggrastr::geom_point_rast(shape = 21, # Specify shape and colour as fixed local parameters
                             colour = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = 0.10, linetype = "dashed") +
    geom_vline(xintercept = -0.10, linetype = "dashed") +
    scale_fill_manual(values = cols) + # Modify point colour
    scale_size_manual(values = sizes) + # Modify point size
    scale_alpha_manual(values = alphas) + # Modify point transparency
    scale_x_continuous(breaks = c(seq(-0.5, 0.5, 0.1)), limits = c(-0.5, 0.5)) +
    scale_y_continuous(breaks = c(seq(0, 200, 50)), limits = c(0, 200)) +
    labs(
      title = glue("{ko} KO EZH2i vs {ko} WT EZH2i"),
      subtitle = "",
      x = "fold enrichment",
      y = "-log10 p-value",
      fill = " "
    ) +
    guides(
      alpha = FALSE,
      size = FALSE,
      fill = guide_legend(override.aes = list(size = 5))
    ) +
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
    geom_text_repel(
      label = labels,
      size = 4,
      max.overlaps = 300,
      min.segment.length = 0.5
    ) # add labels
  volc
  
  ggsave(
    glue("{result_folder}{ko}_EZH2i_volcanoes.png"),
    plot = volc,
    width = 10,
    height = 10,
    dpi = 300,
  )
  
  ggsave(
    glue("{result_folder}{ko}_EZH2i_volcanoes.pdf"),
    plot = volc,
    width = 15,
    height = 10,
    device = "pdf"
  )
  
  # WT non-trt vs. WT EZH2i
  # res = fread(
  #   glue(
  #     "../data/bulk_ATAC-Seq/KO_experiments/TOBIAS/{ko}/bindetect_results.txt"
  #   )
  # )
  # colnames(res) = str_to_upper(colnames(res))
  # res = res %>%
  #   select(NAME,
  #          ends_with("SORTED_FOOTPRINT_PVALUE"),
  #          ends_with("CHANGE")) %>%
  #   select(NAME, starts_with(glue("ATAC_WT_HESC_NT")), ) 
  # colnames(res) = c(
  #   "gene_name",
  #   "WT_NT_vs_WT_EZH2i_pvalue",
  #   "WT_NT_vs_WT_EZH2i_fp_change"
  # )
  # 
  # volc_input = res %>%
  #   mutate(
  #     group = case_when(
  #       WT_NT_vs_WT_EZH2i_fp_change > 0.10 & WT_NT_vs_WT_EZH2i_pvalue < 0.05 ~ "up",
  #       WT_NT_vs_WT_EZH2i_fp_change < -0.10 &
  #         WT_NT_vs_WT_EZH2i_pvalue < 0.05 ~ "down",
  #       .default = "unaltered"
  #     )
  #   ) %>%
  #   mutate(sign_label = case_when(group == "up" ~ gene_name, group == "down" ~ gene_name, .default = NA))
  # 
  # volc_input = volc_input %>% 
  #   dplyr::arrange(desc(group))
  # labels = factor(volc_input$sign_label, levels = unique(volc_input$sign_label))
  # 
  # volc2 = volc_input %>%
  #   ggplot(aes(
  #     x = WT_NT_vs_WT_EZH2i_fp_change,
  #     y = -log10(WT_NT_vs_WT_EZH2i_pvalue),
  #     fill = group,
  #     size = group,
  #     alpha = group
  #   )) +
  #   ggrastr::geom_point_rast(shape = 21, # Specify shape and colour as fixed local parameters
  #                            colour = "black") +
  #   geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  #   geom_vline(xintercept = 0.10, linetype = "dashed") +
  #   geom_vline(xintercept = -0.10, linetype = "dashed") +
  #   scale_fill_manual(values = cols) + # Modify point colour
  #   scale_size_manual(values = sizes) + # Modify point size
  #   scale_alpha_manual(values = alphas) + # Modify point transparency
  #   scale_x_continuous(breaks = c(seq(-0.5, 0.5, 0.1)), limits = c(-0.5, 0.5)) +
  #   scale_y_continuous(breaks = c(seq(0, 200, 50)), limits = c(0, 200)) +
  #   labs(
  #     title = glue("{ko} WT non-trt vs {ko} WT EZH2i"),
  #     subtitle = "",
  #     x = "fold enrichment",
  #     y = "-log10 p-value",
  #     fill = " "
  #   ) +
  #   guides(
  #     alpha = FALSE,
  #     size = FALSE,
  #     fill = guide_legend(override.aes = list(size = 5))
  #   ) +
  #   theme_minimal() +
  #   theme(
  #     text = element_text(size = 16),
  #     legend.text = element_text(size = 14),
  #     plot.title = element_text(size = 12, face = "bold"),
  #     plot.subtitle = element_text(size = 6),
  #     axis.text.x = element_text(size = 14, color = "black"),
  #     axis.text.y = element_text(size = 14, color = "black"),
  #     axis.title.x = element_text(size = 20, color = "black"),
  #     axis.title.y = element_text(size = 20, color = "black")
  #   ) +
  #   geom_text_repel(
  #     label = labels,
  #     size = 4,
  #     max.overlaps = 300,
  #     min.segment.length = 0.5
  #   ) # add labels
  # volc2
  # 
  # volcs = ggarrange(volc, volc2)
  # 
  # ggsave(
  #   glue("{result_folder}{ko}_EZH2i_volcanoes.png"),
  #   plot = volcs,
  #   width = 15,
  #   height = 10,
  #   dpi = 300,
  # )
  # 

  return(volc)
  
}

# experiments
kos = c("CDX1", "CDX2", "HAND1", "MEIS2", "SOX9", "TBX3", "TBXT", "TWIST1")
lapply(kos, make_volcanoes)


get_mean_score = function(ko) {
  res = fread(
    glue(
      "../data/bulk_ATAC-Seq/KO_experiments/TOBIAS/{ko}/bindetect_results.txt"
    )
  )
  colnames(res) = str_to_upper(colnames(res))
  res = res %>%
    dplyr::select(NAME, starts_with(glue("ATAC_{ko}"))) %>%
    dplyr::select(NAME, ends_with("MEAN_SCORE"))
  colnames(res) = c("gene_name", ko)
  res = res %>% filter_at(., vars(ko), any_vars(. > 0.25))
  
  return(res)
}
mean_scores = lapply(kos, get_mean_score)
mean_scores = mean_scores %>% reduce(full_join, by = "gene_name")
rows = mean_scores$gene_name
mean_scores = as.matrix(mean_scores[,-1])
rownames(mean_scores) = rows
mean_scores = mean_scores[complete.cases(mean_scores),]

png(
  file = glue("{result_folder}mean_footprint_score_hm.png"),
  width = 15,
  height = 15,
  units = 'cm',
  res = 500
)
col_fun = colorRamp2(c(0.25, 0.375, 0.50), c(viridis(3)[1], viridis(3)[2], viridis(3)[3]))
hm = Heatmap(
  mean_scores,
  column_title = "KO",
  column_title_side = "bottom",
  row_title = "transcription factor",
  name = "mean footprint score",
  clustering_method_rows = "complete",
  col = col_fun,
  show_column_dend = TRUE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  heatmap_width = unit(5, "cm"),
  heatmap_height = unit(12, "cm"),
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 90
)
print(hm)
dev.off()

pdf(
  file = glue("{result_folder}mean_footprint_score_hm.pdf"),
  width = 4,
  height = 5
)
print(hm)
dev.off()

get_changes = function(ko) {
  res = fread(
    glue(
      "../data/bulk_ATAC-Seq/KO_experiments/TOBIAS/{ko}/bindetect_results.txt"
    )
  )
  colnames(res) = str_to_upper(colnames(res))
  res = res %>%
    dplyr::select(NAME, starts_with(glue("ATAC_{ko}"))) %>%
    dplyr::select(NAME, contains("_ATAC_WT_HESC_EZH2I")) %>%
    dplyr::select(NAME, ends_with(glue("_FOOTPRINT_CHANGE")))
  colnames(res) = c("gene_name", glue("{ko} (EZH2i)"))
  
  return(res)
}

fs_changes = lapply(kos, get_changes)
fs_changes = fs_changes %>% reduce(full_join, by = "gene_name")
rows = fs_changes$gene_name
fs_changes = as.matrix(fs_changes[,-1])
rownames(fs_changes) = rows
fs_changes = fs_changes[complete.cases(fs_changes),]

png(
  file = glue("{result_folder}footprint_score_change_hm.png"),
  width = 8,
  height = 13,
  units = 'cm',
  res = 500
)
col_fun = colorRamp2(c(-0.05, 0, 0.05), c(viridis(3)[1], viridis(3)[2], viridis(3)[3]))
hm2 = Heatmap(
  fs_changes,
  column_title = "KO",
  column_title_side = "bottom",
  row_title = "transcription factor",
  name = "fold change",
  row_km = 2,
  column_km = 1,
  #clustering_method_rows = "complete",
  col = col_fun,
  show_column_dend = TRUE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  heatmap_width = unit(5, "cm"),
  heatmap_height = unit(12, "cm"),
  row_names_gp = gpar(fontsize = 0.7),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 90
)
print(hm2)
dev.off()

pdf(
  file = glue("{result_folder}footprint_score_change_hm.pdf"),
  width = 4,
  height = 5
)
print(hm2)
dev.off()

png(
  file = glue("{result_folder}footprint_score_change_hm.png"),
  width = 10,
  height = 13,
  units = "cm",
  res = 400
)
print(hm2)
dev.off()

# cluster with higher fold changes in KO vs. nt data
higher_fcs = fread("../results/TOBIAS/KO_experiments/TOBIAS_genes_w_increased_fc_in_KO_vs_nt.txt")
higher_fcs = fs_changes[rownames(fs_changes) %in% higher_fcs$gene_name,]

png(
  file = glue("{result_folder}footprint_score_change-higher_fc_in_KOvsNT_hm.png"),
  width = 8,
  height = 13,
  units = 'cm',
  res = 500
)
col_fun = colorRamp2(c(-0.25, 0, 0.25), c(viridis(3)[1], viridis(3)[2], viridis(3)[3]))
hm3 = Heatmap(
  higher_fcs,
  column_title = "KO",
  column_title_side = "bottom",
  row_title = "transcription factor",
  name = "fold change",
  # row_km = 4,
  # column_km = 1,
  clustering_method_rows = "complete",
  col = col_fun,
  show_column_dend = TRUE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  heatmap_width = unit(5, "cm"),
  heatmap_height = unit(12, "cm"),
  row_names_gp = gpar(fontsize = 3),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 90
)
print(hm3)
dev.off()

pdf(
  file = glue("{result_folder}footprint_score_change-higher_fc_in_KOvsNT_hm.pdf"),
  width = 4,
  height = 5
)
print(hm3)
dev.off()

# cluster with higher fold changes in KO vs. nt data
lower_fcs = fread("../results/TOBIAS/KO_experiments/TOBIAS_genes_w_decreased_fc_in_KO_vs_nt.txt")
lower_fcs = fs_changes[rownames(fs_changes) %in% lower_fcs$gene_name,]

png(
  file = glue("{result_folder}footprint_score_change-lower_fc_in_KOvsNT_hm.png"),
  width = 8,
  height = 13,
  units = 'cm',
  res = 500
)
col_fun = colorRamp2(c(-0.25, 0, 0.25), c(viridis(3)[1], viridis(3)[2], viridis(3)[3]))
hm4 = Heatmap(
  lower_fcs,
  column_title = "KO",
  column_title_side = "bottom",
  row_title = "transcription factor",
  name = "fold change",
  # row_km = 4,
  # column_km = 1,
  clustering_method_rows = "complete",
  col = col_fun,
  show_column_dend = TRUE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  heatmap_width = unit(5, "cm"),
  heatmap_height = unit(12, "cm"),
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 90
)
print(hm4)
dev.off()

pdf(
  file = glue("{result_folder}footprint_score_change-lower_fc_in_KOvsNT_hm.pdf"),
  width = 4,
  height = 5
)
print(hm4)
dev.off()
