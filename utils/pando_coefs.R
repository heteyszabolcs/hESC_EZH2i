library("ggpubr")
library("ggplot2")
library("tidyverse")


new_candidates = intersect(nt_coefs %>% dplyr::filter(padj < 0.05) %>% pull(tf) %>% unique(), trt_coefs %>% dplyr::filter(padj < 0.05) %>% pull(tf) %>% unique())

nt_coefs = nt_coefs %>% dplyr::filter(tf %in% new_candidates) %>% mutate(condition = "non-treated")
trt_coefs = trt_coefs %>% dplyr::filter(tf %in% new_candidates) %>% mutate(condition = "EZH2i 7D")

coefs_df = rbind(nt_coefs, trt_coefs) %>% dplyr::filter(padj < 0.05)

library("ggpubr")
x = factor(coefs_df$condition, levels = c("non-treated", "EZH2i 7D"))
ggplot(coefs_df, aes(x = x, y = estimate)) +
  geom_boxplot(color = "black", fill = "#9ecae1") +
  #scale_fill_brewer(palette = "Reds") +
  ylim(0, 8) +
  labs(
    title = "Common transcription factors",
    x = "",
    y = "PANDO coefficient"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  stat_compare_means(label.y = 6.8, label.x = 1.3, size = 10) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 50)


order = coefs_df %>% group_by(target) %>% summarize(max_coef = max(estimate)) %>% arrange(desc(max_coef)) %>% pull(target)

x = factor(coefs_df$target, levels = order)
ggplot(coefs_df, aes(x = x, y = estimate)) +
  geom_boxplot(color = "black", fill = "#9ecae1") +
  #scale_fill_brewer(palette = "Reds") +
  ylim(0, 8) +
  labs(
    title = "Coefficient distribution of candidates",
    x = "",
    y = "PANDO coefficient"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 9),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  stat_compare_means(label.y = 6.8, label.x = 1.3, size = 10) +
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", label.y = 50)
