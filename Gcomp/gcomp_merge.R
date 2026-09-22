#!/usr/bin/env Rscript

# G-computation merge script
# Verifies all 30 tasks are present and valid before merging into df_gcomp_long.rds

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# 1. Argument parsing
args <- commandArgs(trailingOnly = TRUE)
params <- list(
  group_id    = "all",
  sex         = NA_character_,
  grid_file   = "groups_grid.csv",
  results_dir = NA_character_
)

for (a in args) {
  if (grepl("^--group=", a)) {
    params$group_id <- sub("^--group=", "", a)
  } else if (grepl("^--sex=", a)) {
    params$sex <- sub("^--sex=", "", a)
  } else if (grepl("^--grid=", a)) {
    params$grid_file <- sub("^--grid=", "", a)
  } else if (grepl("^--results_dir=", a)) {
    params$results_dir <- sub("^--results_dir=", "", a)
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

group_ids_to_merge <- if (params$group_id == "all") {
  groups_grid$group_id
} else {
  as.integer(trimws(strsplit(params$group_id, ",")[[1]]))
}

# 2. Merge function with integrity check
merge_one_group <- function(gid) {
  group <- groups_grid[groups_grid$group_id == gid, , drop = FALSE]
  if (nrow(group) == 0) {
    stop("Group ID ", gid, " not found in grid")
  }
  
  target_dir <- if (!is.na(params$results_dir) && nzchar(params$results_dir)) {
    file.path(params$results_dir, group$exposure, group$sex, group$age_class)
  } else {
    file.path(script_dir, "results", group$exposure, group$sex, group$age_class)
  }
  
  if (!dir.exists(target_dir)) {
    stop("Directory not found: ", target_dir)
  }
  
  expected_tasks <- 1:30
  expected_files <- file.path(target_dir, sprintf("gcomp_%s_%d.rds", group$exposure, expected_tasks))
  
  missing_tasks   <- integer(0)
  corrupted_tasks <- integer(0)
  loaded_dfs      <- list()
  
  for (t_id in expected_tasks) {
    f_path <- expected_files[t_id]
    
    if (!file.exists(f_path) || file.info(f_path)$size == 0) {
      missing_tasks <- c(missing_tasks, t_id)
      next
    }
    
    obj <- tryCatch(readRDS(f_path), error = function(e) NULL)
    if (is.null(obj) || !is.data.frame(obj) || nrow(obj) == 0) {
      corrupted_tasks <- c(corrupted_tasks, t_id)
      next
    }
    
    required_cols <- c("draw", "Scen_A_Total", "Scen_B_Indirect", "Scen_C_Direct")
    if (!all(required_cols %in% names(obj))) {
      corrupted_tasks <- c(corrupted_tasks, t_id)
      next
    }
    
    loaded_dfs[[length(loaded_dfs) + 1]] <- obj
  }
  
  failed_tasks <- sort(unique(c(missing_tasks, corrupted_tasks)))
  
  if (length(failed_tasks) > 0) {
    cat("Error: Group ", gid, " (", group$exposure, " ", group$sex, " ", group$age_class, 
        ") has ", length(failed_tasks), " missing/corrupted tasks.\n", sep = "")
    if (length(missing_tasks) > 0) {
      cat("Missing tasks: ", paste(missing_tasks, collapse = ", "), "\n", sep = "")
    }
    if (length(corrupted_tasks) > 0) {
      cat("Corrupted tasks: ", paste(corrupted_tasks, collapse = ", "), "\n", sep = "")
    }
    cat("Rerun command: sbatch --array=", paste(failed_tasks, collapse = ","), 
        " submit_slurm.sh --group=", gid, "\n", sep = "")
    stop("Merge aborted for group ", gid)
  }
  
  df_combined <- dplyr::bind_rows(loaded_dfs) |>
    filter(!is.na(Scen_A_Total), !is.na(Scen_B_Indirect), !is.na(Scen_C_Direct))
  
  df_long <- df_combined |>
    pivot_longer(
      cols      = starts_with("Scen_"),
      names_to  = "Scenario",
      values_to = "surv_prob_marginal"
    ) |>
    mutate(
      Scenario = recode(Scenario,
        Scen_A_Total    = "Total Effect",
        Scen_B_Indirect = "Indirect Effect (mediation)",
        Scen_C_Direct   = "Direct Effect"
      )
    )
  
  out_long_file <- file.path(target_dir, "df_gcomp_long.rds")
  saveRDS(df_long, out_long_file)
  cat("Group ", gid, " merged: ", out_long_file, " (", nrow(df_long), " rows)\n", sep = "")
  return(out_long_file)
}

# 3. Run merges
successes <- 0
for (gid in group_ids_to_merge) {
  res <- try(merge_one_group(gid), silent = FALSE)
  if (!inherits(res, "try-error")) {
    successes <- successes + 1
  }
}

cat("Merge completed: ", successes, " / ", length(group_ids_to_merge),sep = "")
