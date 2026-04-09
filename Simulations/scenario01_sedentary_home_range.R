# ============================================================================
# SCENARIO 1: THE SEDENTARY HOME RANGE (ORNSTEIN-UHLENBECK)
# ============================================================================
# Purpose:
# Reproduces the first synthetic stress test in the manuscript.
# - Homogeneous OU-like movement without an explicit home-range covariate in the fitted model.
# - Checks are made with the four diagnostic pillars (topology, MSD, sinuosity,
#   connectivity if a barrier was provided).
# ============================================================================

library(amt)
library(tidyverse)
library(terra)
library(sf)
library(circular)
library(mvtnorm)

# Source function implementing the generative validation workflow.
source("R/diagnose_issf.R")

set.seed(777)

# ============================================================================
# 1. Helper: robust iSSA fitting for synthetic trajectory
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
# 2. Simulate scenario
# ============================================================================
if (!dir.exists("manuscript/figures")) {
  dir.create("manuscript/figures", recursive = TRUE)
}

if (!dir.exists("results")) {
  dir.create("results")
}

r1 <- rast(
  nrows = 100,
  ncols = 100,
  xmin = 0,
  xmax = 1000,
  ymin = 0,
  ymax = 1000
)
values(r1) <- runif(ncell(r1))
names(r1) <- "env"

n_steps <- 1000
t1 <- tibble(x = 500, y = 500, t = as.POSIXct("2024-01-01"))
for (i in 1:n_steps) {
  p <- t1[nrow(t1), ]
  attraction_x <- -(p$x - 500) * 0.3
  attraction_y <- -(p$y - 500) * 0.3
  dx <- rnorm(1, 0, 20) + attraction_x
  dy <- rnorm(1, 0, 20) + attraction_y
  t1 <- bind_rows(t1, tibble(x = p$x + dx, y = p$y + dy, t = p$t + 3600))
}
trk1 <- make_track(t1, x, y, t, crs = 32601)

# ============================================================================
# 3. Fit model and run validation
# ============================================================================
message("Running Scenario 1 (Sedentary home range)...")
m1 <- fit_robust_issa(trk1, r1)
res1 <- diagnose_issf_robust(
  m1,
  trk1,
  r1,
  "env",
  n_sims = 99
)
saveRDS(res1, file.path("results", "scenario01_sedentary.rds"))

# Save diagnostics dashboard (same object plotted as in manuscript workflow)
ggsave(
  filename = "manuscript/figures/scenario1_sedentary.png",
  plot = plot(res1),
  width = 10,
  height = 8,
  dpi = 120
)

# ============================================================================
# 4. Publication figure: MSD failure mode for sedentary dynamics
# ============================================================================
# This is the figure used in the manuscript as Figure msd_failure_sedentary.
msd_sims <- res1$stats$msd_sims_mat
msd_obs <- res1$stats$msd_obs_curve
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
  geom_line(aes(y = Median, linetype = "Simulated (Median)"), color = "grey30", linewidth = 0.8) +
  geom_line(aes(y = Observed, linetype = "Observed"), color = "black", linewidth = 1) +
  scale_linetype_manual(
    name = "",
    values = c("Observed" = "solid", "Simulated (Median)" = "dashed")
  ) +
  labs(
    x = "Time Lag (steps)",
    y = expression(Mean ~ Squared ~ Displacement ~ (m^2)),
    title = "The Diffusion Trap (Sedentary vs Diffusive)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = c(0.2, 0.8),
    legend.background = element_blank(),
    plot.title = element_text(face = "bold", size = 12)
  )

ggsave(
  filename = "manuscript/figures/msd_failure_sedentary.pdf",
  plot = p_msd,
  width = 6,
  height = 5,
  device = "pdf"
)

message("Scenario 1 complete. Outputs saved in manuscript/figures and results/")
