#!/bin/bash
# Submit G-computation SLURM array jobs for all demographic groups
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
grid_file="${script_dir}/groups_grid.csv"
slurm_script="${script_dir}/submit_slurm.sh"

if [ ! -f "${grid_file}" ]; then
  echo "Error: groups_grid.csv not found."
  exit 1
fi

dry_run=false
if [ "${1:-}" == "--dry-run" ]; then
  dry_run=true
fi

tail -n +2 "${grid_file}" | while IFS=',' read -r group_id exposure exposure_col sex age_class age_select target_age_scaled mass_model surv_model has_repro; do
  cmd=("sbatch" "--array=1-30" "${slurm_script}" "--group=${group_id}")
  if [ "${dry_run}" = true ]; then
    echo "[DRY RUN] ${cmd[*]}"
  else
    echo "Submitting group ${group_id}: ${exposure} ${sex} ${age_class}"
    "${cmd[@]}"
  fi
done
