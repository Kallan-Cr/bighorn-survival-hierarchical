# Mortality in an alpine ungulate increases nonlinearly with snow depth.

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![R version](https://img.shields.io/badge/R-4.5.1-blue)](https://www.r-project.org/)
[![renv](https://img.shields.io/badge/renv-1.1.7-orange)](https://rstudio.github.io/renv/)
[![Quarto](https://img.shields.io/badge/Quarto-1.5%2B-brightgreen)](https://quarto.org/)
[![brms](https://img.shields.io/badge/brms-2.22%2B-orange)](https://paul-buerkner.github.io/brms/)

---

## Overview

This repository contains all datasets, scripts, and pipelines to reproduce the findings of this study on bighorn sheep (*Ovis canadensis*) annual survival.

We evaluated:
1. How seasonal winter snow metrics (mean snow depth, snow cover duration, median density) and threshold exceedances ($\ge 10\text{ cm}$ to $\ge 100\text{ cm}$) influence annual survival across demographic classes (lambs, prime-aged and senescent females, prime-aged and senescent males).
2. Long-term temporal trends in winter snow events and demographic survival over 45 years of continuous monitoring.
3. The decomposition of total snow effects into **direct effects** and **indirect carry-over effects** (mediated through spring and autumn body mass dynamics) using G-computation.

---

## Repository structure

```
├── README.md                        # Main project documentation
├── main.qmd                         # Master Quarto document compiling the full analysis
├── .Rprofile                        # Auto-activates renv environment on startup
├── renv.lock                        # Lockfile recording exact package versions
├── renv/                            # renv infrastructure and settings
│
├── 01_setup.qmd                     # Package loading
├── 02_custom_functions.qmd          # Functions for summary tables and extraction
├── 03_data_import.qmd               # Data import, collinearity, and variable standardization
├── 04_prelim_age_models.qmd         # Preliminary age structure comparison (linear, quadratic, splines)
├── 05_lamb_survival.qmd             # Bayesian survival models for lambs
├── 06_female_survival.qmd           # Bayesian survival models for adult females
├── 07_male_survival.qmd             # Bayesian survival models for adult males
├── 08_snow_thresholds.qmd           # Survival models across snow depth thresholds
├── 09_threshold_temporal_trends.qmd # Beta-binomial temporal trend models (90 cm & 100 cm)
├── 10_survival_temporal_trends.qmd  # Long-term temporal trends in survival across groups
├── 11_gcomp_analysis.qmd            # G-computation visualization 
├── 12_export_figures.qmd            # Multi-panel figure exports
├── 13_export_ppchecks.qmd           # Posterior predictive checks exports
│
├── Data_hashed.csv                  # Anonymized individual-level life-history and survival data
├── season_summary_thresh.csv        # Annual winter snow metrics and threshold day counts
│
├── Gcomp/                           # Parameterized G-computation cluster pipeline
    ├── README.md                    # HPC submission and documentation
    ├── groups_grid.csv              # Configuration grid (14 demographic groups)
    ├── gcomp_worker.R               # Parameterized counterfactual simulation worker
    ├── gcomp_merge.R                # Output validator and task merger
    ├── gcomp_collect_all.R          # Aggregator computing direct, indirect, and total effects
    ├── gcomp_global_summary.rds     # Consolidated G-computation effect estimates
    ├── submit_slurm.sh              # SLURM array batch job submitter
    └── submit_all_groups.sh         # Multi-group array job submitter

```

---

## Data description

### 1. `Data_hashed.csv`
Individual-level dataset (1979–2023). Each row corresponds to an individual-year observation:

| Column | Type | Description |
| :--- | :--- | :--- |
| `ID` | character | SHA-256 anonymized unique individual identifier |
| `Yr` | integer | Observation year ($1979–2023$) |
| `Age` | integer | Individual age in years ($0 = \text{lamb}$) |
| `Sex` | character | Sex (`Female` / `Male`) |
| `wtd114` | numeric | Adjusted autumn body mass in kg (standardized to September 15) |
| `pop_den` | integer | Total female (>1yr) count in year $t$ |
| `Survival` | integer | Binary overwinter survival status to year $t+1$ ($1 = \text{survived}, 0 = \text{died}$) |

### 2. `season_summary_thresh.csv`
Annual winter snow conditions derived from daily simulation with SNOWPACK ([Crémel et al. 2026](https://doi.org/10.1080/15230430.2026.2627695)):

| Column | Type | Description |
| :--- | :--- | :--- |
| `season` | character | Winter season label (e.g. `1979-1980`) |
| `Yr` | integer | Winter start year ($1979–2023$) |
| `days_count` | integer | Duration of snow cover |
| `avg_snow_depth` | numeric | Seasonal mean snow depth (cm) |
| `max_snow_depth` | numeric | Seasonal maximum snow depth (cm) |
| `median_snow_depth` | numeric | Seasonal median snow depth (cm) |
| `avg_density` | numeric | Seasonal mean snow density ($\text{kg/m}^3$) |
| `median_density` | numeric | Seasonal median snow density ($\text{kg/m}^3$) |
| `avg_SWE` | numeric | Seasonal mean Snow Water Equivalent (mm) |
| `days_gt_10` ... `days_gt_100` | integer | Number of days with snow depth exceeding $10, 20, \dots, 100\text{ cm}$ |
| `JJ_start` / `JJ_end` | integer | Day of the year of the start and end of snow season |

---

## Statistical modeling framework

All models were developed in a Bayesian framework using [`brms`](https://paul-buerkner.github.io/brms/) and Stan ([`cmdstanr`](https://mc-stan.org/cmdstanr/)):

1. **Survival Models**:
Individual annual survival is modeled using Bayesian mixed-effects logistic regression.
   
2. **Snow Threshold Trends**:
Long-term temporal trends in snow depth threshold (days exceeding 90 cm and 100 cm) are modeled using beta-binomial regression
   
3. **Model selection and validation**:
Approximate leave-one-out cross-validation (`loo` package) and Pareto k diagnostic checks ($k < 0.7$).
Posterior predictive checks with grouped bars, means, standard deviations, and zero proportions.
Prior sensitivity diagnostics (`priorsense` power-scaling).

---

## G-Computation mediation pipeline (`Gcomp/`)

To quantify how much of the total effect of winter conditions on survival operates directly versus indirectly through over-winter mass loss, we implemented a causal mediation analysis using **G-computation**.

**Prerequisite:** Before executing the G-computation pipeline, you must compute the survival models from this study as well as the body mass models from [10.17605/OSF.IO/F5M9V](https://doi.org/10.17605/OSF.IO/F5M9V), as the counterfactual simulations rely on their fitted outputs.

The pipeline evaluates **14 demographic groups** defined in `Gcomp/groups_grid.csv`:
* **Exposures**: Snow cover duration (`days_count`), Median snow density (`median_density`), and Mean snow depth (`avg_snow_depth`).
* **Sexes**: Females and Males.
* **Age classes at $t+1$**: Yearlings (age 1), Prime-aged (age 3), and Senescent (age 9).

### HPC execution workflow
Each group simulates counterfactual outcomes across 30 exposure values and 1,000 posterior draws:
1. **Launch array jobs** (e.g., on a SLURM cluster):
   ```bash
   cd Gcomp
   # Submit all groups (30 tasks per group, 8 cores and 16 GB per task):
   bash submit_all_groups.sh
   ```
2. **Validate and merge task outputs**:
   ```bash
   Rscript --vanilla gcomp_merge.R --group=all
   ```
3. **Aggregate global summaries**:
   ```bash
   Rscript --vanilla gcomp_collect_all.R
   ```

The summary table is saved directly in `Gcomp/gcomp_global_summary.rds` and read by `11_gcomp_analysis.qmd` to produce Figure 5 and summary tables.

---

## Reproducibility & environment setup

This project uses [`renv`](https://rstudio.github.io/renv/) to guarantee exact package version.

### 1. Restore the R environment
After cloning the repository, open R in the project root directory and run:

```r
install.packages("renv")
renv::restore()
```

This will automatically install the exact package versions recorded in `renv.lock`.

### 2. Install CmdStan
Models require CmdStan as backend for `brms`:

```r
cmdstanr::install_cmdstan()
```

### 3. Render the complete analysis
To run all models and compile the HTML document with all tables and figures:

```bash
quarto render main.qmd
```

---

## License & attribution

* **Code** (`.R`, `.qmd`, `.sh` scripts): [MIT License](LICENSE)
* **Data, output, and documentation** (`Data_hashed.csv`, `season_summary_thresh.csv`, `README.md`...): [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/)

See the [`LICENSE`](LICENSE) file for full terms. This project depends on third-party R
packages (see `renv.lock`) distributed under their own licenses (MIT, GPL, etc.); these
are not redistributed here and remain governed by their respective terms.

### Data use policy & citation

When reusing this material, please cite the corresponding publication.

**Note on data use:** These datasets are the result of decades of continuous individual-level field monitoring. While shared under a CC BY 4.0 license, in accordance with ethics and professional courtesy, prospective users are **strongly requested to contact the authors** prior to initiating new analyses, to prevent project overlap and to discuss potential collaboration.
