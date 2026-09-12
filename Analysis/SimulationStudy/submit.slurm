#!/bin/bash
#SBATCH --job-name=sacat_simulation
#SBATCH --array=1-100
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=7-00:00:00

set -euo pipefail

: "${SACAT_PROJECT_ROOT:?Set SACAT_PROJECT_ROOT to the SACAT repository}"
: "${SACAT_SIMULATION_ROOT:?Set SACAT_SIMULATION_ROOT to the simulation output directory}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd "$SACAT_PROJECT_ROOT"
Rscript Analysis/SimulationStudy/run_simulation.R \
  replicate "${SLURM_ARRAY_TASK_ID}"
