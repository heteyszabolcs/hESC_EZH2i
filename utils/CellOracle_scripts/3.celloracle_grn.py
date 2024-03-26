## packages
import numpy as np
import pandas as pd

# vis
import matplotlib.pyplot as plt

# AnnData Scanpy
import scanpy as sc
import anndata as ad

# for CellOracle
import celloracle as co

result_folder = "../results/GRN/CellOracle/"

var_tfs = "../data/GRN/Variable_TFs-Jaspar2024_CISBP.tsv"
var_tfs = pd.read_csv(var_tfs, sep='\t')
var_tfs = var_tfs["TF_name"].to_list()

counts = "../results/scRNA-Seq/hESC_EZH2i_scRNA_Seq-nontrt.csv"
counts = pd.read_csv(counts, sep=',', index_col = "Unnamed: 0")

rna = ad.read_h5ad("../results/GRN/CellOracle/hESC_EZH2i_scRNA_Seq-nt_pseudot.h5ad")
var_features = rna.var['var.features'].cat.categories.to_list()
var_features = [x for x in var_features if str(x) != 'NA']
rna = rna[:,var_features]
oracle = co.Oracle()

counts = counts.loc[var_features,:]
rna.layers["counts"] = counts.T
rna.X = rna.layers["counts"].copy()

oracle.import_anndata_as_raw_count(adata=rna,
                                   cluster_column_name="cluster_EML",
                                   embedding_name="X_umap")

# add Pando predictions and base GRN
base_grn = pd.read_parquet(result_folder + "nt_base_GRN_dataframe.parquet", engine='pyarrow')
pando_predictions = pd.read_csv("../results/GRN/Pando/Var_TFs_of_Jaspar_CISBP/nt_aggr_coef_table-adjp0.05.tsv",
                                sep = "\t")

pando_dict = {}
for TF, TGs in zip(pando_predictions.tf, pando_predictions.target_genes):
    # convert target gene to list
    TG_list = TGs.split(",")
    # store target gene list in a dictionary
    pando_dict[TF] = TG_list
pando_dict = co.utility.inverse_dictionary(pando_dict)

oracle.import_TF_data(TF_info_matrix=base_grn)
oracle.addTFinfo_dictionary(pando_dict)

# KNN imputation
oracle.perform_PCA()
plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
plt.savefig(result_folder + "nt_PCA_explained_var.png")

n_cell = oracle.adata.shape[0]
print(f"cell number is :{n_cell}")
k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}")

oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=4)

# saving CellOracle object
oracle.to_hdf5(result_folder + "nt.celloracle.oracle")

# GRN calculation
oracle = co.load_hdf5(result_folder + "nt.celloracle.oracle")
links = oracle.get_links(cluster_name_for_GRN_unit="cluster_EML", alpha=10,
                         verbose_level=10)

for cluster in list(links.links_dict.keys()):
    links.links_dict[cluster].to_csv(result_folder + f"raw_nt_GRN_for_{cluster}.csv")

# save links
links.to_hdf5(file_path=result_folder + "nt_links.celloracle.links")

# Network preprocessing
# remove uncertain edges
links.filter_links(p=0.05, weight="coef_abs", threshold_number=2000)
# Calculate network scores.
links.get_network_score()
links.merged_score.head()
links.merged_score.to_csv(result_folder + "nt_network_score_table.tsv", sep = "\t")

# Save Links object - this object will be use in the in silico perturbation
links.to_hdf5(file_path=result_folder + "nt_links.celloracle.links")

# filtered edges
for cluster in list(links.links_dict.keys()):
    df = links.filtered_links[cluster]
    df.to_csv(result_folder + "nt_" + cluster + "-filtered_edges.tsv", sep = "\t")

plt.rcParams["figure.figsize"] = [5, 5]
plt.rcParams['figure.dpi'] = 400
plt.rcParams["figure.autolayout"] = False

## visualizations
for cluster in list(links.links_dict.keys()):
    links.plot_scores_as_rank(cluster=cluster, n_gene=30,
                              save=None)
    plt.savefig(result_folder + "nt_" + cluster + "-ranked_score_figures.pdf")
    plt.close()

# Compare GRN score between two clusters
links.plot_score_comparison_2D(value="betweenness_centrality",
                               cluster1="ELC", cluster2="TLC",
                               percentile=98)
plt.savefig(result_folder + "nt_ELC_vs_TLC_networkscore_comp.pdf")
plt.close()
links.plot_score_comparison_2D(value="betweenness_centrality",
                               cluster1="ELC", cluster2="MeLC",
                               percentile=98)
plt.savefig(result_folder + "nt_ELC_vs_MeLC_networkscore_comp.pdf")
plt.close()

links.plot_score_per_cluster(goi="HAND1")
plt.savefig(result_folder + "nt_score_dynamics_in_cluster_EML_2000_HAND1.pdf")
plt.close()
links.plot_score_per_cluster(goi="CDX2")
plt.savefig(result_folder + "nt_score_dynamics_in_cluster_EML_2000_CDX2.pdf")
plt.close()
links.plot_score_per_cluster(goi="ARID3A")
plt.savefig(result_folder + "nt_score_dynamics_in_cluster_EML_2000_ARID3A.pdf")
plt.close()
links.plot_score_per_cluster(goi="NANOG")
plt.savefig(result_folder + "nt_score_dynamics_in_cluster_EML_2000_NANOG.pdf")
plt.close()

links.plot_score_discributions(values=["degree_centrality_all"],
                               method="boxplot")
plt.savefig(result_folder + "nt_degree_centrality_boxplots.pdf")
plt.close()
links.plot_score_discributions(values=["eigenvector_centrality"],
                               method="boxplot")
plt.savefig(result_folder + "nt_eigenvector_centrality_boxplots.pdf")
plt.close()
