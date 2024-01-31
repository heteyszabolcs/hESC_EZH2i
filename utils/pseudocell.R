library("Seurat")
library("Signac")
library("glue")

require(plyr, quietly=T)
require(tidyverse, quietly=T)
require(viridis, quietly=T)
require(gridExtra, quietly=T)
require(dbscan, quietly=T)
require(matrixStats, quietly=T)
require(igraph, quietly=T)
require(ape, quietly=T)
require(reticulate, quietly=T)
require(rdist, quietly=T)
require(ggdendro, quietly=T)
require(tidygraph, quietly=T)
require(ggraph, quietly=T)

# export folder
result_folder = "../results/Seurat_integration/"

# Seurat object (coembedded)
nt_coembed = readRDS("../results/Seurat_integration/nt_scATAC_scRNA_coembed.Rds")
trt_coembed = readRDS("../results/Seurat_integration/trt_scATAC_scRNA_coembed.Rds")

# function: force-directed k-nearest neighbour graph
# source: https://www.nature.com/articles/s41586-019-1654-9
# source: https://github.com/quadbio/organoid_regulomes - pseudocells.R file
find_pseudocells <- function(
    object,
    resolution = 0.3,
    k = 10,
    order = 1,
    reduction = 'pca',
    verbose = TRUE,
    return_object = TRUE
){
  invert_list <- function(lov){
    split(rep(names(lov), lengths(lov)), unlist(lov))
  }
  
  knn_mat <- FindNeighbors(
    object@reductions$umap@cell.embeddings,
    compute.SNN = FALSE,
    k.param = k,
    verbose = verbose
  )$nn
  
  n_cells <- ncol(knn_mat)
  mask <- runif(n_cells, 0, 1) <= resolution
  
  g <- graph_from_adjacency_matrix(knn_mat)
  tnn <- ego(g, order=order, nodes=colnames(knn_mat)[mask])
  tnn <- map(tnn, names)
  names(tnn) <- colnames(knn_mat)[mask]
  print(paste0('Found ', length(tnn), ' territories'))
  terr <- invert_list(tnn)
  nterr <- mean(map_int(terr, length))
  print(paste0('On average ', round(nterr, 1), ' territories per cell'))
  
  nsing <- 0
  terr_df <- map_dfr(colnames(knn_mat), function(x){
    if (!is.null(terr[[x]])){
      row <- c(x, sample(terr[[x]], 1))
    } else {
      row <- c(x, x)
      nsing <- nsing + 1
    }
    names(row) <- c('cell', 'pseudocell')
    return(row)
  })
  nps <- length(unique(terr_df$pseudocell))
  print(paste0('Identified ', nps, ' pseudocells',
              ' for ', ncol(knn_mat), ' cells'))
  if (return_object){
    return(AddMetaData(object, column_to_rownames(terr_df, 'cell')))
  } else {
    return(terr_df)
  }
}

nt_coembed_pseudo = find_pseudocells(object = nt_coembed)
nt_meta = nt_coembed_pseudo@meta.data
trt_coembed_pseudo = find_pseudocells(object = trt_coembed)
trt_meta = trt_coembed_pseudo@meta.data

# align cells from different data types
nt_alignment = nt_meta %>% dplyr::select(pseudocell)
nt_alignment = tibble(id1 = rownames(nt_alignment), id2 = nt_alignment$pseudocell) %>% 
  mutate(type1 = ifelse(str_detect(id1, "EZH2i"), TRUE, FALSE)) %>% 
  mutate(type2 = ifelse(grepl("EZH2i", id2), FALSE, TRUE)) %>% 
  dplyr::filter(type1 == TRUE & type2 == TRUE)

trt_alignment = trt_meta %>% dplyr::select(pseudocell)
trt_alignment = tibble(id1 = rownames(trt_alignment), id2 = trt_alignment$pseudocell) %>% 
  mutate(type1 = ifelse(str_detect(id1, "EZH2i"), TRUE, FALSE)) %>% 
  mutate(type2 = ifelse(grepl("EZH2i", id2), FALSE, TRUE)) %>% 
  dplyr::filter(type1 == TRUE & type2 == TRUE)

# write out aligned cell ids
write_tsv(nt_alignment, glue("{result_folder}nt_pseudocell_alignment.tsv"))
write_tsv(trt_alignment, glue("{result_folder}trt_pseudocell_alignment.tsv"))
