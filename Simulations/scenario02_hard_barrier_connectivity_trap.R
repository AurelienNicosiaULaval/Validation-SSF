# ============================================================================
# SCENARIO 2: THE HARD BARRIER (CONNECTIVITY TRAP)
# ============================================================================
# Purpose:
# Reproduces the second synthetic stress test in the manuscript.
# The observed process is reflected at x = 595, while the fitted model ignores
# the barrier and treats the landscape as smooth habitat classes.
# ============================================================================

library(amt)
library(tidyverse)
library(terra)
library(sf)
library(circular)
library(mvtnorm)

source("R/diagnose_issf.R")

set.seed(777)

# ============================================================================
# 1. Helper
# ============================================================================
fit_robust_issa <- function(trk, map) {
  df <- trk |>
    steps() |>
    random_steps(n_control = 15) |>
    extract_covariates(map)

  df_clean <- df |>
    filter(!is.na(sl_), !is.na(ta_), sl_ > 0.001) |>
    mutate(
      log_sl = log(sl_ + 0.01),
      sl_scaled = scale(sl_),
      cos_ta = cos(ta_)
    ) |>
    na.omit()

  if (var(df_clean$env) < 1e-6) {
    f_full <- case_ ~ sl_scaled + log_sl + cos_ta + strata(step_id_)
  } else {
    f_full <- case_ ~ env + sl_scaled + log_sl + cos_ta + strata(step_id_)
  }

  tryCatch(
    {
      fit_issf(df_clean, f_full)
    },
    error = function(e) {
      f_simple <- update(f_full, . ~ . - cos_ta)
      fit_issf(df_clean, f_simple)
    }
  )
}

# ============================================================================
# 2. Simulate barrier scenario
# ============================================================================
if (!dir.exists("manuscript/figures")) {
  dir.create("manuscript/figures", recursive = TRUE)
}
if (!dir.exists("results")) {
  dir.create("results")
}

r2 <- rast(
  nrows = 100,
  ncols = 100,
  xmin = 0,
  xmax = 1000,
  ymin = 0,
  ymax = 1000
)
values(r2) <- 0
r2[, 60:100] <- 1
names(r2) <- "env"

t2 <- tibble(x = 580, y = 500, t = as.POSIXct("2024-01-01"))
for (i in 1:800) {
  p <- t2[nrow(t2), ]
  dx <- rnorm(1, 5, 30)
  dy <- rnorm(1, 0, 30)
  nx <- p$x + dx
  if (nx > 595) {
    nx <- 595 - (nx - 595)
  }
  t2 <- bind_rows(t2, tibble(x = nx, y = p$y + dy, t = p$t + 3600))
}
trk2 <- make_track(t2, x, y, t, crs = 32601)

# ============================================================================
# 3. Fit model and validate
# ============================================================================
message("Running Scenario 2 (Hard barrier)...")
m2 <- fit_robust_issa(trk2, r2)
res2 <- diagnose_issf_robust(
  m2,
  trk2,
  r2,
  "env",
  n_sims = 99,
  barrier = list(axis = "x", val = 595)
)

saveRDS(res2, file.path("results", "scenario02_barrier.rds"))

ggsave(
  filename = "manuscript/figures/scenario2_barrier.png",
  plot = plot(res2),
  width = 10,
  height = 8,
  dpi = 120
)

message("Scenario 2 complete. Outputs saved in manuscript/figures and results/")
