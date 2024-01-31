#!/bin/bash -l
#SBATCH -A naiss2023-22-84
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH -M rackham
#SBATCH -J signac

module load R
module load R_packages

cd /crex/proj/snic2020-6-3/SZABOLCS/hESC_EZH2i/utils

Rscript motif_analysis.R