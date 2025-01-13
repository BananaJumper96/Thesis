#!/bin/bash -l
#SBATCH -J metaGEMmod
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 128
#SBATCH --time=1-00:00:00
#SBATCH -p batch
#SBATCH --qos=normal
#SBATCH -o "/mnt/lscratch/users/ohickl/Thesis/slurm_logs/%x-%j.out"

cd /mnt/lscratch/users/ohickl/Thesis

# Environment variables
ENV_MANAGER="micromamba"  # or conda, mamba, depending on your environment manager
ENV_NAME="/mnt/aiongpfs/projects/imp/dev/envs/IMP_snakemake"  # Replace with your zDB conda environment name

# Get available cores from slumr env var
CORES=${SLURM_CPUS_PER_TASK}

# Initialize micromamba
echo "Initializing ${ENV_MANAGER}..."
eval "$(${ENV_MANAGER} shell hook --shell=bash)"

# Activate snakemake environment
echo "Activating snakemake environment ${ENV_NAME}..."
if ! command -v ${ENV_MANAGER} &> /dev/null; then
    echo "${ENV_MANAGER} not found. Ensure it is properly installed."
    exit 1
fi
${ENV_MANAGER} activate ${ENV_NAME} || {
    echo "Failed to activate environment ${ENV_NAME}."
    exit 1
}

snakemake --configfile /mnt/lscratch/users/ohickl/Thesis/config/config.yaml \
          --unlock

snakemake --use-conda \
          --conda-prefix /mnt/lscratch/users/ohickl/Thesis/Envs \
          --configfile /mnt/lscratch/users/ohickl/Thesis/config/config.yaml \
          --cores $CORES \
          --rerun-trigger mtime