library("ggpubr")
library("tidyverse")


# TOBIAS output
fc = read_tsv("../results/TOBIAS/BINDetect_HTSELEX_output/bindetect_results.txt")
sign_in_trt = fc %>% arrange(desc(trt_nt_change)) %>% top_n(30, wt = trt_nt_change)
sign_in_trt_genes = sign_in_trt %>% dplyr::select(output_prefix) %>% 
  separate(output_prefix, sep = "\\.", into = c("V1", "V2", "V3", "V4")) %>% 
  dplyr::select("V3", "V4") %>% mutate("V1" = paste(V3, as.character(V4), sep = ".")) %>% 
  dplyr::select(V1)

nt_beds = "../results/TOBIAS/BINDetect_HTSELEX_output/nt_bound_beds"
trt_beds = "../results/TOBIAS/BINDetect_HTSELEX_output/trt_bound_beds"

system(paste0("mkdir ", trt_beds))
system(paste0("mkdir ", nt_beds))

for(i in sign_in_trt_genes$V1) {
  system(paste0("cp ", "../results/TOBIAS/BINDetect_HTSELEX_output/*", i, "/beds/*", i, "_nt_bound.bed ", 
         nt_beds))
  system(paste0("cp ", "../results/TOBIAS/BINDetect_HTSELEX_output/*", i, "/beds/*", i, "_trt_bound.bed ", 
         trt_beds))
}

gene_symbols = sign_in_trt %>% separate(name, sep = "\\.", into = c("V1", "V2", "V3")) %>% 
  dplyr::select(output_prefix, V3)

write_tsv(gene_symbols, "../results/TOBIAS/BINDetect_HTSELEX_output/motifname2genes.tsv")
