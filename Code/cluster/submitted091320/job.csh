#!/bin/bash
#SBATCH -n 32 # Number of cores

#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-01:30 # Runtime in D-HH:MM
#SBATCH -p murphy # Partition to submit to
#SBATCH --mem-per-cpu=1000 # Memory pool for all cores in MB
#SBATCH -o output_%j.out # File to which STDOUT will be written; %j inserts job ID
#SBATCH -e errors_%j.err # File to which STDERR will be written
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=xmeng@g.harvard.edu

rm errors_* slurm-*

module load R/3.6.3-fasrc01

export R_LIBS_USER=$HOME/apps/R_3.6.3:$R_LIBS_USER

R CMD BATCH '/n/holyscratch01/murphy_lab/Users/xmeng/submitted091320/2020309-Plasmode_Sim_CV-TMLE.R'

