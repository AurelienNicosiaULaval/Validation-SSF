# ==============================================================================
# CASE STUDY: RED DEER EMPIRICAL APPLICATION
# ==============================================================================
# Reproduces the empirical application presented in the manuscript.
# Fits a standard iSSA to a Red Deer trajectory and validates it using
# the multi-criteria generative framework.
# ==============================================================================

library(amt)
library(terra)
library(tidyverse)
library(ggplot2)

# Source the robust diagnostic function
source("R/diagnose_issf.R")

set.seed(2025)

# ==============================================================================
# 1. DATA LOADING AND PREPARATION
# ==============================================================================

data("deer")
data("sh_forest")

# A. Prepare the Environmental Map
map_deer <- terra::rast(sh_forest)
names(map_deer) <- "forest"

# B. Prepare the Observed Trajectory
# Audit data
cat("\n=== DATA AUDIT ===\n")
cat("Total fixes:", nrow(deer), "\n")
time_diffs <- as.numeric(diff(deer$t_), units = "hours")
cat("Sampling interval (hours):\n")
cat("  Median:", round(median(time_diffs), 2), "\n")
cat("  Mean:", round(mean(time_diffs), 2), "\n")
cat("  Range:", round(min(time_diffs), 2), "-", round(max(time_diffs), 2), "\n")
cat("Time span:", as.numeric(diff(range(deer$t_)), units = "days"), "days\n")

# Resample to regular 6-hour intervals
trk_real <- deer |>
  track_resample(rate = hours(6), tolerance = minutes(30)) |>
  filter_min_n_burst(min_n = 2) |>
  select(x_, y_, t_)
# Note: analyzing the full trajectory to test for non-stationarity
# head(300) removed

cat("\nAfter resampling and filtering:\n")
cat("  Total points:", nrow(trk_real), "\n")
cat("  Total steps:", nrow(trk_real) - 1, "\n")
cat(
  "  Time span:",
  round(as.numeric(diff(range(trk_real$t_)), units = "days"), 1),
  "days\n"
)
cat("==================\n\n")

# ==============================================================================
# 2. iSSA DATA PREPARATION (STEPS & CONTROL)
# ==============================================================================

# Convert points to steps, generate random available steps, and extract covariates
ssf_dat_deer <- trk_real |>
  steps() |>
  random_steps(n_control = 25) |> # 1:25 ratio as described in the paper
  extract_covariates(map_deer) |>
  mutate(
    cos_ta = cos(ta_),
    log_sl = log(sl_)
  ) |>
  na.omit()

# ==============================================================================
# 3. VISUALIZATION: SAMPLING DESIGN
# ==============================================================================
# Visualizing the "Used" vs "Available" steps to confirm the design.

# Split data into Used (Observed) and Available (Random/Control) steps
obs_steps <- ssf_dat_deer |> filter(case_ == TRUE)
rnd_steps <- ssf_dat_deer |> filter(case_ == FALSE)

p_design <- ggplot() +
  # 1. CONTROL STEPS (Background)
  geom_segment(
    data = rnd_steps,
    aes(x = x1_, y = y1_, xend = x2_, yend = y2_),
    color = "grey70",
    alpha = 0.3,
    linewidth = 0.3
  ) +
  # 2. OBSERVED TRAJECTORY (Foreground)
  geom_segment(
    data = obs_steps,
    aes(x = x1_, y = y1_, xend = x2_, yend = y2_),
    color = "firebrick",
    linewidth = 0.8
  ) +
  # 3. Endpoints
  geom_point(
    data = obs_steps,
    aes(x = x2_, y = y2_),
    color = "black",
    size = 0.8
  ) +
  # Formatting
  coord_equal() +
  theme_minimal() +
  labs(
    title = "iSSA Sampling Design",
    subtitle = "Red = Observed Trajectory | Grey = Available Control Steps (n=25)",
    x = "Easting (m)",
    y = "Northing (m)"
  )

print(p_design)
ggsave(
  "manuscript/figures/deerSSF.png",
  p_design,
  width = 8,
  height = 6,
  dpi = 300
)

# ==============================================================================
# 4. MODEL FITTING
# ==============================================================================

# Fit the conditional logistic regression
# Model structure: Selection for forest + Movement kernel (Step Length & Turn Angle)
m_deer <- fit_issf(
  ssf_dat_deer,
  case_ ~ forest + sl_ + log_sl + cos_ta + strata(step_id_)
)

# Display model coefficients
print(summary(m_deer))

# ==============================================================================
# 5. GENERATIVE VALIDATION (ROBUST)
# ==============================================================================

# Apply the multi-criteria validation framework
# We use K=99 simulations to allow for a minimum P-value of 0.01.
diag_results <- diagnose_issf_robust(
  model = m_deer,
  obs_track = trk_real,
  map = map_deer,
  env_covar_name = "forest",
  n_sims = 99,
  n_cands = 25 # Matches the n_control used in fitting
)

# A. Print Scorecard (Table)
print(diag_results)

# B. Detailed Summary
summary(diag_results)

# C. Visualization Dashboard
png(
  "manuscript/figures/result_deer.png",
  width = 12,
  height = 10,
  units = "in",
  res = 300
)
plot(diag_results)
dev.off()
