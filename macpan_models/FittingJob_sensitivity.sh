#!/bin/bash
#SBATCH --account=def-earn
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=81
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-user=zhaok22@mcmaster.ca
#SBATCH --mail-type=ALL

module load r/4.5
export R_LIBS=~/local/R_libs/
Rscript fit_sensitivity.R
