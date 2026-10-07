#!/bin/bash
#SBATCH --job-name=sacat_meta_cdi
#SBATCH --array=1-200
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err

set -euo pipefail

# Submit from this directory, or set METASPAR_SCRIPT_DIR to its absolute path.
script_dir="${METASPAR_SCRIPT_DIR:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}}"
export METASPAR_ROOT="${METASPAR_ROOT:-${script_dir}/output}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

mode="${METASPAR_MODE:-replicate}"
if [[ "$mode" == "replicate" ]]; then
    Rscript --vanilla "${script_dir}/run_simulation.R" replicate "${SLURM_ARRAY_TASK_ID:?Supply an array task number.}"
else
    Rscript --vanilla "${script_dir}/run_simulation.R" "$mode"
fi
