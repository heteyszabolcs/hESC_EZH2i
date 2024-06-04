# packages
suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("ggplot2")
  library("glue")
  library("ComplexHeatmap")
  library("circlize")
  library("enrichR")
  library("ggpubr")
})



annot_peaks = fread("../results/Seurat_integration/scATAC_annotation/TLC_spec_peaks_up-log2fc_2_adjp0.05_own_annot.tsv")
genes = annot_peaks$`Gene Name`

dbs = c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")
enrichr_analysis = enrichr(genes, dbs)
biol = enrichr_analysis[["GO_Biological_Process_2023"]]
biol = biol %>% dplyr::filter(0.05 > P.value) 

plot = ggplot(biol[1:15, ], aes(x = reorder(Term, -log10(P.value)), y = -log10(P.value))) +
  geom_bar(stat = 'identity', fill = "#9ecae1", color = "black") +
  coord_flip() +
  labs(title = "",
       y = "-log10(p-value)",
       x = "GO Biologocial process (2023)") +
  theme_classic() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  )
plot

ggsave(
  glue("{result_folder}enrichR_GO-TLC_spec_peaks.png"),
  plot = plot,
  width = 12,
  height = 10,
  dpi = 300,
)

ggsave(
  glue("{result_folder}enrichR_GO-TLC_spec_peaks.pdf"),
  plot = plot,
  width = 12,
  height = 10
)
