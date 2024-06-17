if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  "tidyverse",
  "data.table",
  "ggplot2",
  "glue",
  "RColorBrewer",
  "ComplexHeatmap",
  "circlize"
)

# output folder
result_folder = "../results/HOMER/"

clean_homer_result = function(dirs, label) {
  homer_results = list()
  for (dir in list.dirs(dirs)) {
    if (file.exists(glue("{dir}/knownResults.txt"))) {
      tf = str_split(dir, dirs)[[1]][2]
      tf_name = str_split(tf, "_CnT")[[1]][1]
      homer_res = glue("{dir}/knownResults.txt")
      homer_res = fread(homer_res)
      homer_res = homer_res %>% mutate(`% of Target Sequences with Motif` =
                                         as.numeric(str_remove(
                                           `% of Target Sequences with Motif`, "%"
                                         ))) %>%
        dplyr::filter(`% of Target Sequences with Motif` > 5) %>%
        dplyr::filter(`P-value` < 0.05) %>%
        top_n(n = 20, wt = `% of Target Sequences with Motif`) %>%
        separate(`Motif Name`,
                 sep = "/",
                 into = "Motif_name",
                 remove = TRUE) %>%
        dplyr::select(Motif_name,
                      `P-value`,
                      `Log P-value`,
                      pct_of_target_seq_with_motif =
                        `% of Target Sequences with Motif`) %>% 
        arrange(desc(pct_of_target_seq_with_motif)) %>% 
        mutate(sample = label) %>% 
        mutate(tf = tf_name)
    }
    homer_results[[tf]] = homer_res
  }
  return(homer_results)
}

non_trt = clean_homer_result("../results/HOMER/bulk_CnT/non_trt/", label = "NT")
trt = clean_homer_result("../results/HOMER/bulk_CnT/EZH2i_7D//", label = "EZH2i 7D")

non_trt = rbindlist(non_trt)
non_trt = non_trt %>% mutate(minuslog10p_value = -log10(as.numeric(`P-value`)))

trt = rbindlist(trt)
trt = trt %>% mutate(minuslog10p_value = -log10(as.numeric(`P-value`)))

# non trt circle heatmap plot
y_order_non_trt = non_trt %>% group_by(Motif_name) %>% 
  summarise(mean_pct = mean(pct_of_target_seq_with_motif)) %>% 
  arrange(mean_pct) %>% pull(Motif_name)
y_order_non_trt = factor(non_trt$Motif_name, levels = y_order_non_trt)

x_order_non_trt = non_trt %>% group_by(tf) %>% 
  summarise(mean_logp = mean(`Log P-value`)) %>% 
  arrange(mean_logp) %>% pull(tf)
x_order_non_trt = factor(non_trt$tf, levels = x_order_non_trt)

nt_hm = ggplot(non_trt, aes(x_order_non_trt, y_order_non_trt, fill = pct_of_target_seq_with_motif, 
                    size = minuslog10p_value)) +
  geom_point(shape = 21, stroke = 0) +
  scale_x_discrete(position = "top") +
  scale_radius(range = c(1, 5)) +
  scale_fill_gradient2(low = "#a6bddb", mid = "white", high = "#fc9272", midpoint = 25,
                       breaks = c(0, 25, 50), limits = c(0, 50)) +
  theme_minimal() +
  theme(legend.position = "right", 
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 12)) +
  guides(size = guide_legend(override.aes = list(fill = NA, color = "black", stroke = .25), 
                             label.position = "bottom",
                             title.position = "top", 
                             order = 1),
         fill = guide_colorbar(ticks.colour = NA, title.position = "top", order = 2)) +
  labs(size = "-log10(p-value)", fill = "Motif %", x = NULL, y = NULL) +
  ggtitle("HOMER motif analysis - non-treated")
nt_hm

ggsave(
  glue("{result_folder}HOMER_bubblegum_plot_nontrt.pdf"),
  width = 7,
  height = 7,
  plot = nt_hm
)

y_order_trt = trt %>% group_by(Motif_name) %>% 
  summarise(mean_pct = mean(pct_of_target_seq_with_motif)) %>% 
  arrange(mean_pct) %>% pull(Motif_name)
y_order_trt = factor(trt$Motif_name, levels = y_order_trt)

x_order_trt = trt %>% group_by(tf) %>% 
  summarise(mean_logp = mean(`Log P-value`)) %>% 
  arrange(mean_logp) %>% pull(tf)
x_order_trt = factor(trt$tf, levels = x_order_trt)

trt_hm = ggplot(trt, aes(x_order_trt, y_order_trt, fill = pct_of_target_seq_with_motif, 
                            size = minuslog10p_value)) +
  geom_point(shape = 21, stroke = 0) +
  scale_x_discrete(position = "top") +
  scale_radius(range = c(1, 5)) +
  scale_fill_gradient2(low = "#a6bddb", mid = "white", high = "#fc9272", midpoint = 25,
                       breaks = c(0, 25, 50), limits = c(0, 50)) +
  theme_minimal() +
  theme(legend.position = "right", 
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        legend.title = element_text(size = 12)) +
  guides(size = guide_legend(override.aes = list(fill = NA, color = "black", stroke = .25), 
                             label.position = "bottom",
                             title.position = "top", 
                             order = 1),
         fill = guide_colorbar(ticks.colour = NA, title.position = "top", order = 2)) +
  labs(size = "-log10(p-value)", fill = "Motif %", x = NULL, y = NULL) +
  ggtitle("HOMER motif analysis - EZH2i 7D")
trt_hm

ggsave(
  glue("{result_folder}HOMER_bubblegum_plot_trt.pdf"),
  width = 7,
  height = 7,
  plot = trt_hm
)
