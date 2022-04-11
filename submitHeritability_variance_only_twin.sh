#!/bin/bash

#SBATCH --account=kg98
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-12:00:00
#SBATCH --mem-per-cpu=16000#SBATCH --mail-type=ALL

module load R/3.5.3-mkl
R

Rscript --vanilla S2_ACEheritabilityScript_variance_rand.R twinCovariatesDWI_only_twin.mat TwinAlignment_new.mat
