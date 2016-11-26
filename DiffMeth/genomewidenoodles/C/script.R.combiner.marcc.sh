#! /bin/bash
#SBATCH
#SBATCH --time=100:0:0
#SBATCH --partition=shared

module load R/3.3.1
Rscript fisher.gw.noodles.C.R combiner 50
