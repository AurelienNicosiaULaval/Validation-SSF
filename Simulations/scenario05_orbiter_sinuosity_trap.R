# ============================================================================
# SCENARIO 5: THE ORBITER (SINUOSITY TRAP)
# ============================================================================
# Purpose:
# Reproduces the orbital/structural sinuosity stress test in the manuscript.
# Trajectory follows a near-constant turning pattern; the single-state model
# cannot preserve the circular structure over time.
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

if (!dir.exists("manuscript/figures")) {
  dir.create("manuscript/figures", recursive = TRUE)
}
if (!dir.exists("results")) {
  dir.create("results")
}

# ============================================================================
# 2. Simulate orbiter trajectory
# ============================================================================
r5 <- rast(
  nrows = 100,
  ncols = 100,
  xmin = 0,
  xmax = 1000,
  ymin = 0,
  ymax = 1000
)
values(r5) <- runif(ncell(r5))
names(r5) <- "env"

t5 <- tibble(x = 500, y = 500, t = as.POSIXct("2024-01-01"))
step_len <- 20
turn_ang <- 0.2
current_ang <- 0

for (i in 1:300) {
  p <- t5[nrow(t5), ]
  actual_sl <- rnorm(1, step_len, 2)
  actual_ta <- rnorm(1, turn_ang, 0.05)

  current_ang <- current_ang + actual_ta
  nx <- p$x + actual_sl * cos(current_ang)
  ny <- p$y + actual_sl * sin(current_ang)

  t5 <- bind_rows(t5, tibble(x = nx, y = ny, t = p$t + 3600))
}
trk5 <- make_track(t5, x, y, t, crs = 32601)

# ============================================================================
# 3. Fit model and validate
# ============================================================================
message("Running Scenario 5 (Orbiter)...")
m5 <- fit_robust_issa(trk5, r5)
res5 <- diagnose_issf_robust(m5, trk5, r5, "env", n_sims = 99)

saveRDS(res5, file.path("results", "scenario05_orbiter.rds"))

ggsave(
  filename = "manuscript/figures/scenario5_orbiter.png",
  plot = plot(res5),
  width = 10,
  height = 8,
  dpi = 120
)

print(res5)
message("Scenario 5 complete. Outputs saved in manuscript/figures and results/")
