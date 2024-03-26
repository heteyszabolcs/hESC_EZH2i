## packages
import numpy as np
import pandas as pd
import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm

# vis
import matplotlib.pyplot as plt
import seaborn as sns

# for CellOracle
import celloracle as co
from celloracle import motif_analysis as ma
from celloracle.utility import save_as_pickled_object

## additional modules of Uppmax:
# module load bioinfo-tools
# module load BEDTools

# reference genome
ref_genome = "hg38"
genomes_dir = "/proj/snic2020-6-3/SZABOLCS/reference_data/genome/"

# result_folder
result_folder = "../results/GRN/CellOracle/"

# load cicero outputs
peaks = pd.read_csv("../results/scATAC-Seq/variable_trt_peaks.csv", index_col=0)
peaks = peaks.x.values
cicero_connections = pd.read_csv("../results/scATAC-Seq/trt_cicero_connections.csv", index_col=0)
cicero_connections.head()

# TSS annotation
tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome=ref_genome)
tss_annotated.tail()

# integrate TSS info with cicero outputs
integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss_annotated,
                                               cicero_connections=cicero_connections)
print(integrated.shape)
integrated.head()

# keep strong peaks
peak = integrated[integrated.coaccess >= 0.6]
peak = peak[["peak_id", "gene_short_name"]].reset_index(drop=True)
print(peak.shape)
peak.head()

# save
peak.to_csv(result_folder + "processed_trt_peak_file.csv")

## TF motif scanning for GRN construction
# call processed scATAC-Seq peaks
peaks = pd.read_csv(result_folder + "processed_trt_peak_file.csv", index_col = 0)
peaks.head()
peaks = ma.check_peak_format(peaks, ref_genome, genomes_dir=genomes_dir)

# TF motif scan (long run!)
tfi = ma.TFinfo(peak_data_frame=peaks,
                ref_genome=ref_genome,
                genomes_dir=genomes_dir)
tfi.scan(motifs=None) # using CellOracle's default database

tfi.to_hdf5(file_path=result_folder + "trt.celloracle.tfinfo")

## filtering and final base GRN construction
tfi.reset_filtering()
tfi.filter_motifs_by_score(threshold=10)
tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)

# saving
df = tfi.to_dataframe()
df.head()
df = tfi.to_dataframe()
df.to_parquet(result_folder + "trt_base_GRN_dataframe.parquet")
