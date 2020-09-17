#!/bin/bash
#SBATCH -n 1 # Number of cores

#SBATCH -N 1 # Ensure that all cores are on one machine

#SBATCH -t 0-00:10 # Runtime in D-HH:MM#SBATCH -p shared # Partition to submit to#SBATCH --mem=100 # Memory pool for all cores in MB#SBATCH -o output_%j.out # File to which STDOUT will be written; %j inserts job ID

#SBATCH -e errors_%j.err # File to which STDERR will be written

rm errors_* slurm-*

module load R/3.6.3-fasrc01

export R_LIBS_USER=$HOME/apps/R_3.6.3:$R_LIBS_USER

R CMD BATCH '/n/home05/xmeng/cvtmle/submitted091320/2020309-Plasmode_Sim_CV-TMLE.R'