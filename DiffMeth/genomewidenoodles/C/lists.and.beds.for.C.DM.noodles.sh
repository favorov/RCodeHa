#! /bin/bash
#SBATCH
#SBATCH --time=45:0:0
#SBATCH --partition=shared
#skel for marcc
module load R/3.3.1
Rscript lists.and.beds.for.C.DM.noodles.R
