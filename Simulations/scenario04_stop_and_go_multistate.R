# ============================================================================
# SCENARIO 4: THE STOP-AND-GO (MULTI-STATE TRAP)
# ============================================================================
# Purpose:
# Reproduces the multi-state switching scenario used in the manuscript.
# The data-generating process alternates between long directed and short
# tortuous steps; the fitted model assumes a single-state iSSA.
# ============================================================================

library(amt)
library(tidyverse)
library(terra)
library(sf)
library(circular)
library(mvtnorm)
library(ggplot2)

source("R/diagnose_issf.R")

set.seed(777) # Same seed as scenario definitions across the manuscript scripts.

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
# 2. Simulate stop-and-go trajectory
# ============================================================================
r4 <- rast(
  nrows = 100,
  ncols = 100,
  xmin = 0,
  xmax = 1000,
  ymin = 0,
  ymax = 1000
)
values(r4) <- runif(ncell(r4))
names(r4) <- "env"

t4 <- tibble(x = 500, y = 500, t = as.POSIXct("2024-01-01"))
state <- "transit"
n_steps <- 300

for (i in 1:n_steps) {
  p <- t4[nrow(t4), ]

  if (state == "transit") {
    if (runif(1) < 0.1) state <- "forage"
  } else {
    if (runif(1) < 0.1) state <- "transit"
  }

  if (state == "transit") {
    sl <- rgamma(1, shape = 10, scale = 5)
    ta <- as.numeric(rvonmises(1, mu = 0, kappa = 3))
  } else {
    sl <- rgamma(1, shape = 1, scale = 5)
    ta <- runif(1, -pi, pi)
  }

  prev_ang <- 0
  if (nrow(t4) > 1) {
    prev_ang <- atan2(p$y - t4$y[nrow(t4) - 1], p$x - t4$x[nrow(t4) - 1])
  }
  new_ang <- prev_ang + ta

  t4 <- bind_rows(
    t4,
    tibble(
      x = p$x + sl * cos(new_ang),
      y = p$y + sl * sin(new_ang),
      t = p$t + 3600
    )
  )
}

trk4 <- make_track(t4, x, y, t, crs = 32601)

# ============================================================================
# 3. Fit model and validate
# ============================================================================
message("Running Scenario 4 (Stop-and-go)...")
m4 <- fit_robust_issa(trk4, r4)
res4 <- diagnose_issf_robust(m4, trk4, r4, "env", n_sims = 99)
saveRDS(res4, file.path("results", "scenario04_stopgo.rds"))

# Standard diagnostic dashboard

ggsave(
  filename = "manuscript/figures/scenario4_stopgo.png",
  plot = plot(res4),
  width = 10,
  height = 8,
  dpi = 120
)

# ============================================================================
# 4. Publication figure: MSD staircase mismatch
# ============================================================================
msd_sims <- res4$stats$msd_sims_mat
msd_obs <- res4$stats$msd_obs_curve
msd_quantiles <- apply(
  msd_sims,
  1,
  quantile,
  probs = c(0.025, 0.5, 0.975),
  na.rm = TRUE
)

df_plot <- data.frame(
  Lag = 1:length(msd_obs),
  Observed = msd_obs,
  Median = msd_quantiles[2, ],
  Lower = msd_quantiles[1, ],
  Upper = msd_quantiles[3, ]
)

p_msd <- ggplot(df_plot, aes(x = Lag)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey70", alpha = 0.5) +
  geom_line(
    aes(y = Median, linetype = "Simulated (Median)"),
    color = "grey30",
    linewidth = 0.8
  ) +
  geom_line(
    aes(y = Observed, linetype = "Observed"),
    color = "black",
    linewidth = 1
  ) +
  scale_linetype_manual(
    name = "",
    values = c("Observed" = "solid", "Simulated (Median)" = "dashed")
  ) +
  labs(
    x = "Time Lag (steps)",
    y = expression(Mean ~ Squared ~ Displacement ~ (m^2)),
    title = "Failure to Reproduce Multi-State Dynamics"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = c(0.2, 0.8),
    legend.background = element_blank(),
    plot.title = element_text(face = "bold", size = 12)
  )

ggsave(
  filename = "manuscript/figures/msd_failure_stopgo.pdf",
  plot = p_msd,
  width = 6,
  height = 5,
  device = "pdf"
)

message("Scenario 4 complete. Outputs saved in manuscript/figures and results/")
