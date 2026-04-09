# ============================================================================
# SCENARIO 3: THE CORRIDOR FOLLOWER (TOPOLOGY TRAP)
# ============================================================================
# Purpose:
# Reproduces the corridor stress test in the manuscript.
# The trajectory is strongly pulled toward a narrow channel in x,
# while the fitted model uses only movement and generic landscape covariates.
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
# 2. Simulate corridor trajectory
# ============================================================================
if (!dir.exists("manuscript/figures")) {
  dir.create("manuscript/figures", recursive = TRUE)
}
if (!dir.exists("results")) {
  dir.create("results")
}

r3 <- rast(
  nrows = 100,
  ncols = 100,
  xmin = 0,
  xmax = 1000,
  ymin = 0,
  ymax = 1000
)
values(r3) <- 0
r3[, 48:52] <- 1
names(r3) <- "env"

t3 <- tibble(x = 500, y = 100, t = as.POSIXct("2024-01-01"))
for (i in 1:300) {
  p <- t3[nrow(t3), ]
  dy <- rnorm(1, mean = 20, sd = 5)
  dx <- rnorm(1, 0, 5) - (p$x - 500) * 0.5
  t3 <- bind_rows(t3, tibble(x = p$x + dx, y = p$y + dy, t = p$t + 3600))
}
trk3 <- make_track(t3, x, y, t, crs = 32601)

# ============================================================================
# 3. Fit model and validate
# ============================================================================
message("Running Scenario 3 (Corridor follower)...")
m3 <- fit_robust_issa(trk3, r3)
res3 <- diagnose_issf_robust(m3, trk3, r3, "env", n_sims = 99)

saveRDS(res3, file.path("results", "scenario03_corridor.rds"))

ggsave(
  filename = "manuscript/figures/scenario3_corridor.png",
  plot = plot(res3),
  width = 10,
  height = 8,
  dpi = 120
)

message("Scenario 3 complete. Outputs saved in manuscript/figures and results/")
