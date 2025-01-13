#!/bin/bash -l
#SBATCH -J CarveMe
#SBATCH --mail-type=fail
#SBATCH --mail-user=ricardo.parise@uni.lu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --time=2-00:00:00
#SBATCH -p batch
#SBATCH --qos=normal
#SBATCH -o "slurm_logs/%x-%j.out"

cd /mnt/lscratch/users/rparise/Thesis
micromamba activate snakemake_env
snakemake --use-conda --conda-prefix /mnt/lscratch/users/rparise/Thesis/Envs --cores 32 --rerun-trigger mtime