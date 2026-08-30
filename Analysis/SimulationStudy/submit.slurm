#!/bin/bash
#SBATCH --job-name=dasra_simulation
#SBATCH --array=1-100
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=7-00:00:00

set -euo pipefail

: "${DASRA_PROJECT_ROOT:?Set DASRA_PROJECT_ROOT to the DASRA repository}"
: "${DASRA_SIMULATION_ROOT:?Set DASRA_SIMULATION_ROOT to the simulation output directory}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd "$DASRA_PROJECT_ROOT"
Rscript Analysis/SimulationStudy/run_simulation.R \
  replicate "${SLURM_ARRAY_TASK_ID}"
