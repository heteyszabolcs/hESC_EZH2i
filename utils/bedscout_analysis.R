if (!require("pacman"))
  install.packages("pacman")
pacman::p_load("tidyverse",
               "data.table",
               "ggplot2",
               "glue",
               "GenomicRanges",
               "GenomicFeatures",
               "bedscout"
)

# output folder
result_folder = "../results/bulk_ATAC-Seq/"

# helper function
make_gr = function(path_to_bed, name) {
  bed = fread(path_to_bed)
  bed$V6 = name
  bed = GRanges(seqnames = bed$V1,
                ranges = IRanges(
                  start = bed$V2,
                  end = bed$V3,
                  names = bed$V6
                ))
}

## create GRanges objects
CDX1 = 
  make_gr(path_to_bed = "../data/bulk_ATAC-Seq/KO_experiments/MACS2/ATAC_CDX1_hESC_KO_EZH2i.mRp.clN_peaks.broadPeak",
          name = "CDX1_KO")
CDX2 = 
  make_gr(path_to_bed = "../data/bulk_ATAC-Seq/KO_experiments/MACS2/ATAC_CDX2_hESC_KO_EZH2i.mRp.clN_peaks.broadPeak",
          name = "CDX2_KO")
HAND1 = 
  make_gr(path_to_bed = "../data/bulk_ATAC-Seq/KO_experiments/MACS2/ATAC_HAND1_hESC_KO_EZH2i.mRp.clN_peaks.broadPeak",
          name = "HAND1_KO")
TBXT = 
  make_gr(path_to_bed = "../data/bulk_ATAC-Seq/KO_experiments/MACS2/ATAC_TBXT_hESC_KO_EZH2i.mRp.clN_peaks.broadPeak",
          name = "TBXT_KO")

# Overlap analysis by bedscout
plot_euler(
  list(
    CDX1,
    CDX2,
    HAND1,
    TBXT
  ),
  ignore.strand = TRUE,
  fills = c("#deebf7", "#9ecae1", "#3182bd"),
  names = c("CDX1 KO", "CDX2 KO", "HAND1 KO", "TBXT KO")
)



