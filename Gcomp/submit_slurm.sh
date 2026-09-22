#!/bin/bash
#SBATCH --job-name=gcomp_sim
#SBATCH --array=1-30
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/gcomp_%A_%a.out
#SBATCH --error=logs/gcomp_%A_%a.err

set -euo pipefail

mkdir -p logs

GROUP_ARG="${1:---group=1}"

if command -v module &> /dev/null; then
  module load r 2>/dev/null || module load R 2>/dev/null || true
fi

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MC_CORES="${SLURM_CPUS_PER_TASK:-8}"

Rscript --vanilla gcomp_worker.R "${GROUP_ARG}" --task="${SLURM_ARRAY_TASK_ID:-1}"
