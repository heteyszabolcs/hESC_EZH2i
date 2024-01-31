#!/usr/bin/env bash
#SBATCH -A NAISS2024-22-62 -p alvis
#SBATCH -N 2 --gpus-per-node=T4:8  # We're launching 2 nodes with 8 Nvidia T4 GPUs each
#SBATCH -t 48:00:00
#SBATCH -J pando

module load R

cd /mimer/NOBACKUP/groups/naiss2024-22-62/SZABOLCS/hESC_EZH2i/utils

Rscript pando_alvis.R