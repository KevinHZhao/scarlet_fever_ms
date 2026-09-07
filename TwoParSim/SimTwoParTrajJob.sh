#!/bin/bash
#SBATCH --account=def-earn
#SBATCH --nodes=1
#SBATCH --time=0-02
#SBATCH --cpus-per-task=80
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-user=zhaok22@mcmaster.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020 gcc/9.3.0 r/4.2.1
export R_LIBS=~/local/R_libs/
Rscript SimTwoParTraj.R
