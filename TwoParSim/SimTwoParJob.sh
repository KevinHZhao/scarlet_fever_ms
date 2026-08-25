#!/bin/bash
#SBATCH --account=def-earn
#SBATCH --nodes=1
#SBATCH --time=3-00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-user=zhaok22@mcmaster.ca
#SBATCH --mail-type=ALL

module load gcc/9.3.0 r/4.2.1
export R_LIBS=~/local/R_libs/
Rscript SimTwoParDataSharc.R