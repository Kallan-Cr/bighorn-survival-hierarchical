#!/usr/bin/env Rscript

# Collect and summarize G-computation outputs across groups
# Outputs slopes (median, 95% CI, pd) and absolute survival differences (median)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(bayestestR)
})

args <- commandArgs(trailingOnly = TRUE)
params <- list(
  grid_file = "groups_grid.csv",
  sex       = NA_character_,
  output    = "gcomp_global_summary"
)

for (a in args) {
  if (grepl("^--grid=", a)) {
    params$grid_file <- sub("^--grid=", "", a)
  } else if (grepl("^--sex=", a)) {
    params$sex <- sub("^--sex=", "", a)
  } else if (grepl("^--output=", a)) {
    params$output <- sub("^--output=", "", a)
  }
}

script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)

grid_path <- if (file.exists(params$grid_file)) {
  params$grid_file
} else {
  file.path(script_dir, params$grid_file)
}

if (!file.exists(grid_path)) {
  stop("Configuration grid not found at: ", grid_path)
}
groups_grid <- read.csv(grid_path, stringsAsFactors = FALSE)

if (!is.na(params$sex) && nzchar(params$sex)) {
  groups_grid <- groups_grid |> filter(tolower(sex) == tolower(params$sex))
}

summarize_one_group <- function(group) {
  res_dir <- file.path(script_dir, "results", group$exposure, group$sex, group$age_class)
  long_file <- file.path(res_dir, "df_gcomp_long.rds")
  
  if (!file.exists(long_file)) {
    return(NULL)
  }
  
  df_long <- readRDS(long_file)
  
  raw_col <- paste0(group$exposure_col, "_raw")
  if (!raw_col %in% names(df_long)) {
    cand_cols <- grep("_raw$", names(df_long), value = TRUE)
    if (length(cand_cols) > 0) {
      raw_col <- cand_cols[1]
    } else {
      stop("Raw exposure column missing in: ", long_file)
    }
  }
  
  df_long <- df_long |>
    rename(exposure_raw = all_of(raw_col))
  
  # 1. Slope estimation across MCMC draws
  df_slope <- df_long |>
    group_by(draw, Scenario) |>
    summarize(
      slope = cov(exposure_raw, surv_prob_marginal) / var(exposure_raw),
      .groups = "drop"
    ) |>
    group_by(Scenario) |>
    summarize(
      stats = list(bayestestR::describe_posterior(slope, central_tendency = "median", ci = 0.95)),
      .groups = "drop"
    ) |>
    unnest(stats) |>
    select(
      Scenario,
      slope = Median,
      slope_ci_lower = CI_low,
      slope_ci_upper = CI_high,
      pd
    )
  
  # 2. Absolute difference (surv_max - surv_min)
  df_diff <- df_long |>
    group_by(draw, Scenario) |>
    summarize(
      surv_min = surv_prob_marginal[which.min(exposure_raw)],
      surv_max = surv_prob_marginal[which.max(exposure_raw)],
      diff_abs = surv_max - surv_min,
      .groups = "drop"
    ) |>
    group_by(Scenario) |>
    summarize(
      diff_abs = median(diff_abs),
      .groups = "drop"
    )
  
  # 3. Combine metrics
  df_slope |>
    left_join(df_diff, by = "Scenario") |>
    mutate(
      exposure = group$exposure,
      sex      = group$sex,
      age      = group$age_class
    ) |>
    select(
      exposure, sex, age, Scenario,
      slope, slope_ci_lower, slope_ci_upper, pd,
      diff_abs
    )
}

results_list <- lapply(split(groups_grid, seq_len(nrow(groups_grid))), summarize_one_group)
results_df <- dplyr::bind_rows(results_list)

if (nrow(results_df) == 0) {
  stop("No merged group results found. Please run gcomp_merge.R first.")
}

output_dir <- file.path(script_dir, "global_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

rds_file <- file.path(output_dir, paste0(params$output, ".rds"))
csv_file <- file.path(output_dir, paste0(params$output, ".csv"))

saveRDS(results_df, rds_file)
write.csv(results_df, csv_file, row.names = FALSE)

cat(sprintf("Global summary saved to %s", rds_file))
