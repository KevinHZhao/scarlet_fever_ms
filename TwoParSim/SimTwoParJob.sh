#!/bin/bash
#SBATCH --account=def-earn
#SBATCH --time=00:03:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=256MB
module load r/4.1.1
Rscript SimTwoParDataSharc.R