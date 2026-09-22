# G-Computation Pipeline

This directory contains a generic, parameterized implementation of the G-computation pipeline. Although efforts were made to make the pipeline cluster-agnostic, out-of-the-box portability cannot be guaranteed. Users will need to custom submission scripts and environment configurations to their specific computing cluster.

All 14 demographic group are centralized in `groups_grid.csv`. A single worker (`gcomp_worker.R`) handles all groups, and `gcomp_merge.R` verifies task integrity before merging.

---

## Directory Structure

```
Gcomp/
├── groups_grid.csv          # Configuration grid (14 demographic groups)
├── gcomp_worker.R           # Parameterized simulation worker
├── gcomp_merge.R            # Merge script with integrity verification
├── gcomp_collect_all.R      # Summary (slopes, 95% CI, pd, diff_abs)
├── submit_slurm.sh          # Generic SLURM job array script
├── submit_all_groups.sh     # Launcher for all 14 demographic groups
├── gcomp_global_summary.rds # Contains the simulation results reported in the article
└── README.md                # Documentation
```

---

## Workflow

### 1. Submit SLURM jobs

#### Submit all groups (30 tasks each):
```bash
# Dry run to inspect commands
bash submit_all_groups.sh --dry-run

# Submit all 14 groups
bash submit_all_groups.sh
```

#### Submit a single group:
```bash
# Submit group 1 (Duration, Female, Yearling)
sbatch --array=1-30 submit_slurm.sh --group=1
```

### 2. Merge and check integrity

Once array jobs finish, verify and merge task outputs into `df_gcomp_long.rds`:

```bash
# Merge all groups
Rscript --vanilla gcomp_merge.R --group=all

# Merge a specific group
Rscript --vanilla gcomp_merge.R --group=1
```

If any task is missing or corrupted, the script reports the failed task IDs and prints the exact `sbatch` re-run command:
```
Missing tasks: 14, 28
Rerun command: sbatch --array=14,28 submit_slurm.sh --group=1
```

### 3. Collect global summaries

Consolidate all merged group results into `gcomp_global_summary.rds` and `.csv`:
```bash
Rscript --vanilla gcomp_collect_all.R
```
Outputs are saved in `global_results/`.
