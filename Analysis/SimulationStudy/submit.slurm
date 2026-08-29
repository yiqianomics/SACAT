#!/bin/bash
#SBATCH --job-name=dasra_formal
#SBATCH -p stat,batch
#SBATCH --prefer=sapphire
#SBATCH --array=1-100
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=7-00:00:00
#SBATCH --output=/home/zhang.16383/DORAM/dasra_formal_simulation/logs/%A_%a.out
#SBATCH --error=/home/zhang.16383/DORAM/dasra_formal_simulation/logs/%A_%a.err
#SBATCH --mail-user=zhang.16383@osu.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail

export MAMBA_EXE="$HOME/software/micromamba/bin/micromamba"
export MAMBA_ROOT_PREFIX="$HOME/micromamba"

echo "MAMBA_EXE=${MAMBA_EXE}"
echo "MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX}"

eval "$("$MAMBA_EXE" shell hook --shell bash)"
micromamba activate r451

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd /home/zhang.16383/DORAM

Rscript \
  /home/zhang.16383/DORAM/dasra_formal_hpc_simulation.R \
  replicate \
  "${SLURM_ARRAY_TASK_ID}"