# packages
suppressPackageStartupMessages({
  library("tidyverse")
  library("glue")
  library("data.table")
  library("cicero")
  library("monocle3")
  library("Seurat")
})

# cicero needs monocle3
#devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
cicero_analysis = function(seurat_object_path,
                           label = "trt",
                           result_folder = "../results/scATAC-Seq/",
                           chroms = "../data/hg38.chrom.sizes",
                           nfeature = 50000) {
  ## Seurat scATAC object
  seurat = readRDS(seurat_object_path)
  seurat = FindVariableFeatures(seurat, nfeature = nfeature)
  var_features = VariableFeatures(seurat)
  
  DefaultAssay(seurat) = "peaks"
  counts = GetAssayData(seurat, layer = "counts")[var_features,]
  counts@x[counts@x > 0] <- 1
  
  ## create CDS
  # peak info
  peakinfo = data.frame(peak = rownames(counts)) %>%
    separate(peak, into = c("chr", "bp1", "bp2")) %>%
    mutate(site_name = paste(chr, bp1, bp2, sep = "_"))
  row.names(peakinfo) = peakinfo$site_name
  peakinfo = peakinfo %>% select(-site_name)
  
  # cell info
  cellinfo = data.frame(cell = colnames(counts))
  row.names(cellinfo) = cellinfo$cell
  
  row.names(counts) = row.names(peakinfo)
  colnames(counts) = row.names(cellinfo)
  
  input_cds = suppressWarnings(new_cell_data_set(counts,
                                                 cell_metadata = cellinfo,
                                                 gene_metadata = peakinfo))
  
  input_cds = monocle3::detect_genes(input_cds)
  #Ensure there are no peaks included with zero reads
  input_cds = input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]
  input_cds = input_cds[,Matrix::colSums(exprs(input_cds)) != 0]
  
  ## Qauality check and Filtering
  hist(Matrix::colSums(exprs(input_cds)))
  
  max_count =  15000
  min_count = 200
  input_cds = input_cds[, Matrix::colSums(exprs(input_cds)) >= min_count]
  input_cds = input_cds[, Matrix::colSums(exprs(input_cds)) <= max_count]
  
  ## Process cicero-CDS object
  set.seed(42)
  
  input_cds = detect_genes(input_cds)
  input_cds = estimate_size_factors(input_cds)
  input_cds = preprocess_cds(input_cds, method = "LSI")
  
  # UMAP
  input_cds = reduce_dimension(input_cds,
                               reduction_method = 'UMAP',
                               preprocess_method = "LSI")
  umap_coords = reducedDims(input_cds)$UMAP
  
  # make cicero object
  cicero_cds = make_cicero_cds(input_cds, reduced_coordinates = umap_coords)
  saveRDS(cicero_cds,
          glue("{result_folder}{label}_cicero_cds.Rds"))
  
  # main function: co-accessibility detection
  chroms = fread(chroms)
  conns = run_cicero(cicero_cds, chroms) # long run!
  saveRDS(conns,
          glue("{result_folder}{label}_cicero_connections.Rds"))
  
  # exports
  all_peaks = row.names(exprs(input_cds))
  write.csv(x = all_peaks,
            file = glue("{result_folder}variable_{label}_peaks.csv"))
  write.csv(
    x = conns,
    file = glue("{result_folder}{label}_cicero_connections.csv")
  )
  
}

# run
cicero_analysis(seurat_object_path="../results/scATAC-Seq/hESC_scATAC_nt_Seurat.Rds", label="nt")