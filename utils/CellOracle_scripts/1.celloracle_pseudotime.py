## packages
import copy
import glob
import time
import os
import shutil
import sys
import numpy as np
import pandas as pd

# AnnData Scanpy
import scanpy as sc
import anndata as ad
from tqdm.auto import tqdm

# vis
import matplotlib.pyplot as plt
import seaborn as sns
import plotly
import plotly.express as px
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# for CellOracle
import celloracle as co
from celloracle.applications import Pseudotime_calculator

# output folder
result_folder = "/proj/snic2020-6-3/SZABOLCS/hESC_EZH2i/results/scRNA-Seq/CellOracle/"

# scRNA-Seq anndata
rna = "/proj/snic2020-6-3/SZABOLCS/hESC_EZH2i/results/scRNA-Seq/hESC_EZH2i_scRNA_Seq-nontrt.h5ad"
rna = ad.read_h5ad(rna)

## Pseudotime calculation
rna.obs['cluster_EML'] = rna.obs['cluster_EML'].astype('category').values
rna.uns['cluster_EML_colors'] = ['#8dd3c7', '#bebada', '#ffffb3', '#fb8072', '#80b1d3']
pt = Pseudotime_calculator(adata=rna,
                           obsm_key="X_umap", # Dimensional reduction data name
                           cluster_column_name="cluster_EML"
                           )
pt.plot_cluster(fontsize=8)
#plt.savefig("/proj/snic2020-6-3/SZABOLCS/hESC_EZH2i/results/scRNA-Seq/pt_test.png")

ct_dictionary = {'AMLC': ['AMLC'], 'ELC': ['ELC'],
                 'HLC': ['HLC'],
                 'MeLC': ['MeLC'], 'TLC': ['TLC']}
pt.set_lineage(lineage_dictionary=ct_dictionary)
pt.plot_lineages()
#plt.savefig("/proj/snic2020-6-3/SZABOLCS/hESC_EZH2i/results/scRNA-Seq/pt_test2.png")

def plot(adata, embedding_key, cluster_column_name):
    embedding = adata.obsm[embedding_key]
    df = pd.DataFrame(embedding, columns=["x", "y"])
    df["cluster"] = adata.obs[cluster_column_name].values
    df["label"] = adata.obs.index.values
    fig = px.scatter(df, x="x", y="y", hover_name=df["label"], color="cluster")
    plotly.offline.plot(fig, filename="/proj/snic2020-6-3/SZABOLCS/hESC_EZH2i/results/scRNA-Seq/scRNA_Seq_nt_plotly.html")

#plot(adata=pt.adata,
#    embedding_key=pt.obsm_key,
#    cluster_column_name=pt.cluster_column_name)

# set root cells (based on plotly format of UMAP)
root_cells = {"ELC": "EZH2i_Naive_WT_10X_ATCACAGCACTGTGTA",
              "HLC": "EZH2i_Naive_WT_10X_GTCACTCTCCGGTTCT",
              "MeLC": "EZH2i_Naive_WT_10X_CTATCCGCACATATCG",
              "TLC": "EZH2i_Naive_WT_10X_GACCCTTGTGTATTGC"}
pt.set_root_cells(root_cells=root_cells)
pt.plot_root_cells()
#plt.savefig("/proj/snic2020-6-3/SZABOLCS/hESC_EZH2i/results/results/scRNA-Seq/pt_test3.png")

# calculate neighbors and diffusion map
sc.pp.neighbors(pt.adata, n_neighbors=30)
sc.tl.diffmap(pt.adata)

# calculate pseudotime
pt.get_pseudotime_per_each_lineage()
pt.plot_pseudotime(cmap="rainbow")
pt.adata.obs[["Pseudotime"]].head()

# export as anndata
rna.obs = pt.adata.obs
rna.write_h5ad(result_folder + "hESC_EZH2i_scRNA_Seq-nt_pseudot.h5ad")

# save multiple open images
plt.rcParams["figure.figsize"] = [5, 5]
plt.rcParams["figure.autolayout"] = True

def save_multi_image(filename):
   pp = PdfPages(filename)
   fig_nums = plt.get_fignums()
   figs = [plt.figure(n) for n in fig_nums]
   for fig in figs:
      fig.savefig(pp, format='pdf')
   pp.close()

filename = result_folder + "celloracle_nt_pseudotime_figs.pdf"
save_multi_image(filename)
