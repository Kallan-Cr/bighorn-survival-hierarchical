#!/usr/bin/env Rscript

# G-computation worker script
# Parameterized via groups_grid.csv and command-line arguments

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(brms)
  library(parallel)
  library(splines)
})

# 1. Argument 
args <- commandArgs(trailingOnly = TRUE)
params <- list(
  group_id   = NA_integer_,
  task_id    = NA_integer_,
  n_draws    = 1000L,
  data_dir   = NA_character_,
  output_dir = NA_character_
)

for (a in args) {
  if (grepl("^--group=", a)) {
    params$group_id <- as.integer(sub("^--group=", "", a))
  } else if (grepl("^--task=", a)) {
    params$task_id <- as.integer(sub("^--task=", "", a))
  } else if (grepl("^--n_draws=", a)) {
    params$n_draws <- as.integer(sub("^--n_draws=", "", a))
  } else if (grepl("^--data_dir=", a)) {
    params$data_dir <- sub("^--data_dir=", "", a)
  } else if (grepl("^--output_dir=", a)) {
    params$output_dir <- sub("^--output_dir=", "", a)
  }
}

if (is.na(params$task_id)) {
  slurm_task <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
  params$task_id <- if (nzchar(slurm_task)) as.integer(slurm_task) else 1L
}

if (is.na(params$group_id)) {
  params$group_id <- 1L
}

task_id  <- params$task_id
group_id <- params$group_id
n_draws  <- params$n_draws
n_cores  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = parallel::detectCores() - 1L))
n_cores  <- max(1L, min(n_cores, 8L))

script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)

grid_path <- file.path(script_dir, "groups_grid.csv")
if (!file.exists(grid_path)) {
  stop("Configuration grid not found at: ", grid_path)
}
groups_grid <- read.csv(grid_path, stringsAsFactors = FALSE)

if (!group_id %in% groups_grid$group_id) {
  stop("Invalid group_id: ", group_id)
}
group <- groups_grid[groups_grid$group_id == group_id, , drop = FALSE]

# 2. Data and model loading
find_file <- function(filename, custom_dir = NA) {
  candidates <- c()
  if (!is.na(custom_dir) && nzchar(custom_dir)) {
    candidates <- c(candidates, file.path(custom_dir, filename))
  }
  candidates <- c(
    candidates,
    file.path(script_dir, "data", filename),
    file.path(script_dir, filename),
    file.path(script_dir, "..", "G_comp", filename),
    file.path(script_dir, "..", "rds", filename),
    file.path(script_dir, "..", filename)
  )
  for (cand in candidates) {
    if (file.exists(cand)) return(normalizePath(cand))
  }
  stop("Required file not found: ", filename)
}

base_data_clean  <- readRDS(find_file("base_data_clean.rds", params$data_dir))
Data_females_uns <- readRDS(find_file("Data_females_uns.rds", params$data_dir))
Data_males_uns   <- readRDS(find_file("Data_males_uns.rds", params$data_dir))
Mass_loss_model  <- readRDS(find_file(group$mass_model, params$data_dir))
Survival_model   <- readRDS(find_file(group$surv_model, params$data_dir))

options(mc.cores = 1L)
Sys.setenv(MC_CORES = "1", STAN_NUM_THREADS = "1")

fem_wtd_mean <- mean(Data_females_uns$wtd114, na.rm = TRUE)
fem_wtd_sd   <- sd(Data_females_uns$wtd114,   na.rm = TRUE)
mal_wtd_mean <- mean(Data_males_uns$wtd114,   na.rm = TRUE)
mal_wtd_sd   <- sd(Data_males_uns$wtd114,     na.rm = TRUE)

# 3. Population and exposure setup
exp_col <- group$exposure_col

exposure_raw_seq <- seq(
  from = min(base_data_clean[[exp_col]], na.rm = TRUE),
  to   = max(base_data_clean[[exp_col]], na.rm = TRUE),
  length.out = 30
)

exp_mean <- mean(base_data_clean[[exp_col]], na.rm = TRUE)
exp_sd   <- sd(base_data_clean[[exp_col]],   na.rm = TRUE)
exposure_z_seq <- (exposure_raw_seq - exp_mean) / exp_sd
exposure_ref_z <- min(exposure_z_seq)

if (group$age_class == "Yearling") {
  pop_baseline <- base_data_clean |>
    filter(Age == 0, Sex == group$sex)
} else {
  pop_baseline <- base_data_clean |>
    mutate(age_select = Age) |>
    filter(age_select == group$age_select, Sex == group$sex) |>
    select(-age_select)
}

pop_baseline <- pop_baseline |>
  filter(
    !is.na(wtd114),
    !is.na(next_wtd12),
    !is.na(pop_den),
    !is.na(avg_snow_depth),
    !is.na(days_count),
    !is.na(median_density),
    !is.na(Sex),
    !is.na(wtd114_next),
    !is.na(pop_den_next)
  ) |>
  mutate(across(
    .cols = c(wtd114, pop_den, avg_snow_depth, days_count, median_density, pop_den_next),
    .fns  = ~ scale(.)[, 1]
  ))

target_raw <- exposure_raw_seq[task_id]
target_z   <- exposure_z_seq[task_id]

data_mass_target <- pop_baseline
data_mass_target[[exp_col]] <- target_z

data_mass_ref <- pop_baseline
data_mass_ref[[exp_col]] <- exposure_ref_z

# 4. Simulation across draws
set.seed(123)
tirages_mcmc <- sample(seq_len(brms::ndraws(Mass_loss_model)), n_draws)

compute_one_draw <- function(s) {
  set.seed(s)
  tryCatch({
    if (isTRUE(group$has_repro)) {
      pp_repro_tgt <- brms::posterior_predict(
        Mass_loss_model, resp = "leadrepro",
        newdata = data_mass_target,
        draw_ids = s, allow_new_levels = TRUE, cores = 1
      )
      d_tgt_repro <- data_mass_target |>
        mutate(lead_repro = as.numeric(pp_repro_tgt))
      
      pp_nxt12_tgt <- brms::posterior_predict(
        Mass_loss_model, resp = "nextwtd12",
        newdata = d_tgt_repro,
        draw_ids = s, allow_new_levels = TRUE, cores = 1
      )
      d_tgt_mass <- d_tgt_repro |>
        mutate(next_wtd12 = as.numeric(pp_nxt12_tgt))
      
      pp_w114_tgt <- brms::posterior_predict(
        Mass_loss_model, resp = "wtd114next",
        newdata = d_tgt_mass,
        draw_ids = s, allow_new_levels = TRUE, cores = 1
      )
      
      pp_repro_ref <- brms::posterior_predict(
        Mass_loss_model, resp = "leadrepro",
        newdata = data_mass_ref,
        draw_ids = s, allow_new_levels = TRUE, cores = 1
      )
      d_ref_repro <- data_mass_ref |>
        mutate(lead_repro = as.numeric(pp_repro_ref))
      
      pp_nxt12_ref <- brms::posterior_predict(
        Mass_loss_model, resp = "nextwtd12",
        newdata = d_ref_repro,
        draw_ids = s, allow_new_levels = TRUE, cores = 1
      )
      d_ref_mass <- d_ref_repro |>
        mutate(next_wtd12 = as.numeric(pp_nxt12_ref))
      
      pp_w114_ref <- brms::posterior_predict(
        Mass_loss_model, resp = "wtd114next",
        newdata = d_ref_mass,
        draw_ids = s, allow_new_levels = TRUE, cores = 1
      )
    } else {
      pp_nxt12_tgt <- brms::posterior_predict(
        Mass_loss_model, resp = "nextwtd12",
        newdata = data_mass_target,
        draw_ids = s, allow_new_levels = TRUE, cores = 1
      )
      d_tgt_mass <- data_mass_target |>
        mutate(next_wtd12 = as.numeric(pp_nxt12_tgt))
      
      pp_w114_tgt <- brms::posterior_predict(
        Mass_loss_model, resp = "wtd114next",
        newdata = d_tgt_mass,
        draw_ids = s, allow_new_levels = TRUE, cores = 1
      )
      
      pp_nxt12_ref <- brms::posterior_predict(
        Mass_loss_model, resp = "nextwtd12",
        newdata = data_mass_ref,
        draw_ids = s, allow_new_levels = TRUE, cores = 1
      )
      d_ref_mass <- data_mass_ref |>
        mutate(next_wtd12 = as.numeric(pp_nxt12_ref))
      
      pp_w114_ref <- brms::posterior_predict(
        Mass_loss_model, resp = "wtd114next",
        newdata = d_ref_mass,
        draw_ids = s, allow_new_levels = TRUE, cores = 1
      )
    }
    
    if (group$sex == "Female") {
      mass_target_scaled <- (as.numeric(pp_w114_tgt) - fem_wtd_mean) / fem_wtd_sd
      mass_ref_scaled    <- (as.numeric(pp_w114_ref) - fem_wtd_mean) / fem_wtd_sd
    } else {
      mass_target_scaled <- (as.numeric(pp_w114_tgt) - mal_wtd_mean) / mal_wtd_sd
      mass_ref_scaled    <- (as.numeric(pp_w114_ref) - mal_wtd_mean) / mal_wtd_sd
    }
    
    # Scenario A: Total Effect
    pop_A <- pop_baseline
    pop_A$wtd114     <- mass_target_scaled
    pop_A[[exp_col]] <- target_z
    pop_A$Age        <- group$target_age_scaled
    pop_A$pop_den    <- pop_baseline$pop_den_next
    p_surv_tot <- brms::posterior_epred(
      Survival_model, newdata = pop_A,
      draw_ids = s, allow_new_levels = TRUE, cores = 1
    )
    
    # Scenario B: Indirect Effect
    pop_B <- pop_baseline
    pop_B$wtd114     <- mass_target_scaled
    pop_B[[exp_col]] <- exposure_ref_z
    pop_B$Age        <- group$target_age_scaled
    pop_B$pop_den    <- pop_baseline$pop_den_next
    p_surv_ind <- brms::posterior_epred(
      Survival_model, newdata = pop_B,
      draw_ids = s, allow_new_levels = TRUE, cores = 1
    )
    
    # Scenario C: Direct Effect
    pop_C <- pop_baseline
    pop_C$wtd114     <- mass_ref_scaled
    pop_C[[exp_col]] <- target_z
    pop_C$Age        <- group$target_age_scaled
    pop_C$pop_den    <- pop_baseline$pop_den_next
    p_surv_dir <- brms::posterior_epred(
      Survival_model, newdata = pop_C,
      draw_ids = s, allow_new_levels = TRUE, cores = 1
    )
    
    res_df <- data.frame(
      raw_exposure    = target_raw,
      draw            = s,
      Scen_A_Total    = mean(p_surv_tot, na.rm = TRUE),
      Scen_B_Indirect = mean(p_surv_ind, na.rm = TRUE),
      Scen_C_Direct   = mean(p_surv_dir, na.rm = TRUE)
    )
    names(res_df)[1] <- paste0(exp_col, "_raw")
    return(res_df)
    
  }, error = function(e) {
    warning("Task ", task_id, " error in draw ", s, ": ", conditionMessage(e))
    return(NULL)
  })
}

results_list <- parallel::mclapply(
  tirages_mcmc,
  compute_one_draw,
  mc.cores   = n_cores,
  mc.cleanup = TRUE
)

df_task <- dplyr::bind_rows(results_list)

# 5. Output
out_dir <- if (!is.na(params$output_dir) && nzchar(params$output_dir)) {
  params$output_dir
} else {
  file.path(script_dir, "results", group$exposure, group$sex, group$age_class)
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_file <- file.path(out_dir, sprintf("gcomp_%s_%d.rds", group$exposure, task_id))
saveRDS(df_task, out_file)
cat("Task ", task_id, " complete (", nrow(df_task), " draws): ", out_file, "\n", sep = "")
