#!/bin/bash
#SBATCH --job-name=sacat_joint
#SBATCH --array=1-100
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/joint_%A_%a.out
#SBATCH --error=logs/joint_%A_%a.err

set -euo pipefail
: "${SACAT_PROJECT_ROOT:?Set SACAT_PROJECT_ROOT to the SACAT repository}"
: "${SACAT_JOINT_ROOT:?Set SACAT_JOINT_ROOT to the joint simulation output directory}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd "$SACAT_PROJECT_ROOT"
Rscript Analysis/SimulationInternalValidation/code/run_joint.R replicate "${SLURM_ARRAY_TASK_ID}"
