# ==============================================================================
# SCENARIO 06: RETURN-CONDITIONED SINUOSITY TRAP
# -----------------------------------------------------------------------------
# Strategy:
# - Simulate a long baseline trajectory with the same step-length/turn-angle kernel
#   used in the standard fit.
# - Select a length-N window with very low straightness and small return distance.
# - Fit the standard single-state iSSA and run diagnose_issf_robust with K=99.
# ==============================================================================

library(amt)
library(terra)
library(tidyverse)
library(circular)

source("R/diagnose_issf.R")

# Reproducibility
set.seed(20260213)

BASE_SEED <- 4000L
BASE_SHAPE <- 10
BASE_SCALE <- 2
BASE_KAPPA <- 2
MEAN_SL <- BASE_SHAPE * BASE_SCALE
DT_SEC <- 3600

ENV_MIN <- 0
ENV_MAX <- 1000
X0 <- 500
Y0 <- 500

SEARCH_GRID <- rbind(
  data.frame(N = 300, delta_factor = 1, q = 0.005, M = 5000),
  data.frame(N = 300, delta_factor = 1, q = 0.005, M = 10000),
  data.frame(N = 300, delta_factor = 2, q = 0.005, M = 5000),
  data.frame(N = 300, delta_factor = 1, q = 0.01, M = 5000),
  data.frame(N = 600, delta_factor = 1, q = 0.005, M = 5000),
  data.frame(N = 1000, delta_factor = 1, q = 0.005, M = 5000)
)

SEED_TRIES <- BASE_SEED + seq_len(40)
N_SCREEN <- 20
N_FINAL <- 99

fit_robust_issa <- function(trk, map) {
  df <- trk |
    steps() |
    random_steps(n_control = 15) |
    extract_covariates(map)

  df_clean <- df |
    filter(!is.na(sl_), !is.na(ta_), sl_ > 0.001) |
    mutate(
      log_sl = log(sl_ + 0.01),
      sl_scaled = scale(sl_),
      cos_ta = cos(ta_)
    ) |
    na.omit()

  formula_fit <- if (var(df_clean$env) < 1e-6) {
    case_ ~ sl_scaled + log_sl + cos_ta + strata(step_id_)
  } else {
    case_ ~ env + sl_scaled + log_sl + cos_ta + strata(step_id_)
  }

  tryCatch(
    fit_issf(df_clean, formula_fit),
    error = function(e) fit_issf(df_clean, update(formula_fit, . ~ . - cos_ta))
  )
}

reflect_1d <- function(v) {
  while (v < ENV_MIN || v > ENV_MAX) {
    if (v < ENV_MIN) {
      v <- ENV_MIN + (ENV_MIN - v)
    } else {
      v <- ENV_MAX - (v - ENV_MAX)
    }
  }
  v
}

calc_straightness <- function(x, y) {
  d_end <- sqrt((x[length(x)] - x[1])^2 + (y[length(y)] - y[1])^2)
  d_path <- sum(sqrt(diff(x)^2 + diff(y)^2))
  ifelse(d_path == 0, 0, d_end / d_path)
}

simulate_baseline_long <- function(M, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  df <- tibble(x = X0, y = Y0, t = as.POSIXct("2024-01-01", tz = "UTC"))
  angle <- 0

  for (i in seq_len(M)) {
    p <- df[nrow(df), ]
    sl <- rgamma(1, shape = BASE_SHAPE, scale = BASE_SCALE)
    ta <- as.numeric(rvonmises(1, mu = 0, kappa = BASE_KAPPA))
    angle <- angle + ta

    x <- p$x + sl * cos(angle)
    y <- p$y + sl * sin(angle)
    x <- reflect_1d(x)
    y <- reflect_1d(y)

    df <- bind_rows(df, tibble(x = x, y = y, t = p$t + DT_SEC))
  }

  make_track(df, x, y, t, crs = 32601)
}

select_return_window <- function(long_track, N, delta, q, min_buffer = 80) {
  xs <- long_track$x_
  ys <- long_track$y_
  starts <- seq_len(length(xs) - N)

  if (length(starts) == 0L) {
    return(NULL)
  }

  straightness <- numeric(length(starts))
  end_dist <- numeric(length(starts))
  ok_start <- rep(TRUE, length(starts))

  for (i in seq_along(starts)) {
    idx <- starts[i] + (0:N)
    straightness[i] <- calc_straightness(xs[idx], ys[idx])
    end_dist[i] <- sqrt((xs[idx[N + 1]] - xs[idx[1]])^2 + (ys[idx[N + 1]] - ys[idx[1]])^2)
    ok_start[i] <- all(xs[idx] > ENV_MIN + min_buffer, na.rm = TRUE) &&
      all(xs[idx] < ENV_MAX - min_buffer, na.rm = TRUE) &&
      all(ys[idx] > ENV_MIN + min_buffer, na.rm = TRUE) &&
      all(ys[idx] < ENV_MAX - min_buffer, na.rm = TRUE)
  }

  q_thr <- as.numeric(quantile(straightness, probs = q, na.rm = TRUE))
  cand <- which(straightness <= q_thr)

  if (length(cand) == 0L) {
    return(NULL)
  }

  if (any(ok_start)) {
    cand <- intersect(cand, which(ok_start))
    if (length(cand) == 0L) {
      cand <- which(straightness <= q_thr)
    }
  }

  cand_return <- cand[end_dist[cand] <= delta]
  if (length(cand_return) > 0L) {
    cand <- cand_return
  }

  idx_best <- cand[which.min(straightness[cand])[1]]
  idx <- idx_best + (0:N)

  out <- make_track(
    tibble(x = xs[idx], y = ys[idx], t = long_track$t_[idx]),
    x,
    y,
    t,
    crs = 32601
  )

  list(
    track = out,
    straightness = straightness[idx_best],
    end_dist = end_dist[idx_best],
    start_idx = idx_best
  )
}

run_diag <- function(obs_track, map, n_sims) {
  fit <- fit_robust_issa(obs_track, map)
  diagnose_issf_robust(
    model = fit,
    obs_track = obs_track,
    map = map,
    env_covar_name = "env",
    n_sims = n_sims,
    n_cands = 25
  )
}

out_dir <- "results/scenario06"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

map_env <- rast(
  nrows = 100,
  ncols = 100,
  xmin = ENV_MIN,
  xmax = ENV_MAX,
  ymin = ENV_MIN,
  ymax = ENV_MAX
)
values(map_env) <- 0
names(map_env) <- "env"

search_results <- list()
winner <- NULL
final_diag <- NULL
for (i in seq_len(nrow(SEARCH_GRID))) {
  cfg <- SEARCH_GRID[i, ]
  delta <- cfg$delta_factor * MEAN_SL

  if (!is.null(winner)) {
    break
  }

  for (seed_obs in SEED_TRIES) {
    seed_long <- seed_obs + 10000L
    long_track <- simulate_baseline_long(cfg$M, seed = seed_long)
    sel <- select_return_window(
      long_track,
      N = cfg$N,
      delta = delta,
      q = cfg$q
    )

    if (is.null(sel)) next

    set.seed(seed_obs)
    d <- tryCatch(
      run_diag(sel$track, map_env, n_sims = N_SCREEN),
      error = function(e) NULL
    )
    if (is.null(d)) next

    p <- d$scorecard$P_Value
    row <- tibble(
      seed_obs = seed_obs,
      seed_long = seed_long,
      N = cfg$N,
      M = cfg$M,
      q = cfg$q,
      delta_factor = cfg$delta_factor,
      delta = delta,
      p1 = p[1],
      p2 = p[2],
      p3 = p[3],
      si_obs = sel$straightness,
      end_dist = sel$end_dist,
      start_idx = sel$start_idx
    )

    search_results[[length(search_results) + 1L]] <- row

    cat(sprintf(
      "screen %4d | cfg=%d N=%d M=%d q=%0.3f df=%g | seed=%d | p=(%.3f, %.3f, %.3f) | si=%.4f | end=%.2f\n",
      length(search_results),
      i,
      cfg$N,
      cfg$M,
      cfg$q,
      cfg$delta_factor,
      seed_obs,
      p[1],
      p[2],
      p[3],
      sel$straightness,
      sel$end_dist
    ))

    if (p[1] > 0.1 && p[2] > 0.1 && p[3] < 0.08) {
      d99 <- run_diag(sel$track, map_env, n_sims = N_FINAL)
      p99 <- d99$scorecard$P_Value

      cat(sprintf(
        "final  | N=%d M=%d q=%0.3f df=%g seed=%d | p99=(%.3f, %.3f, %.3f)\n",
        cfg$N,
        cfg$M,
        cfg$q,
        cfg$delta_factor,
        p99[1],
        p99[2],
        p99[3]
      ))

      if (p99[1] > 0.1 && p99[2] > 0.1 && p99[3] < 0.05) {
        winner <- list(
          cfg = cfg,
          delta = delta,
          seed_obs = seed_obs,
          seed_long = seed_long,
          obs_track = sel$track,
          row = row,
          final_p = p99
        )
        final_diag <- d99
        break
      }
    }
  }
}

search_tbl <- bind_rows(search_results)
if (nrow(search_tbl) == 0L) {
  stop("Search produced no valid windows for the requested settings.")
}
search_tbl <- search_tbl %>% arrange(p3)
write.csv(search_tbl, file.path(out_dir, "scenario06_screening.csv"), row.names = FALSE)

if (is.null(winner) || is.null(final_diag)) {
  stop("No K=99 candidate found with p1>0.1, p2>0.1, p3<0.05. Increase search budget.")
}

final_row <- tibble(
  scenario = "07",
  seed_obs = winner$seed_obs,
  seed_long = winner$seed_long,
  N = winner$cfg$N,
  M = winner$cfg$M,
  delta = winner$delta,
  q = winner$cfg$q,
  delta_factor = winner$cfg$delta_factor,
  dt_sec = DT_SEC,
  x0 = winner$obs_track$x_[1],
  y0 = winner$obs_track$y_[1],
  step_kernel = sprintf("Gamma(shape=%.1f, scale=%.1f), vonMises(0, kappa=%.1f)", BASE_SHAPE, BASE_SCALE, BASE_KAPPA),
  p1 = winner$final_p[1],
  p2 = winner$final_p[2],
  p3 = winner$final_p[3],
  p1_screen = winner$row$p1,
  p2_screen = winner$row$p2,
  p3_screen = winner$row$p3,
  si_obs = winner$row$si_obs,
  end_dist = winner$row$end_dist
)

cat("\n=== Final Scenario 7 (K=99) candidate ===\n")
print(final_row)
cat("\n=== Scorecard ===\n")
print(final_diag$scorecard)

write.csv(final_row, file.path(out_dir, "scenario06_params_and_pvals.csv"), row.names = FALSE)

sims_df <- bind_rows(
  lapply(seq_along(final_diag$data$sims), function(i) {
    d <- final_diag$data$sims[[i]]
    d$sim_id <- i
    d
  })
)

p_ud <- ggplot() +
  geom_path(
    data = sims_df,
    aes(x = x_, y = y_, group = sim_id),
    color = "firebrick",
    alpha = 0.2,
    linewidth = 0.2
  ) +
  geom_path(
    data = final_diag$data$obs_track,
    aes(x = x_, y = y_),
    color = "black",
    linewidth = 0.8
  ) +
  coord_equal() +
  labs(
    title = "Scenario 7 - Trajectories",
    subtitle = sprintf(
      "Observed (black) vs simulated (red), N=%d, x0=(%.1f, %.1f), delta=%.1f",
      final_row$N,
      final_row$x0,
      final_row$y0,
      final_row$delta
    ),
    x = "x",
    y = "y"
  ) +
  theme_minimal()

msd_mat <- final_diag$stats$msd_sims_mat
msd_mean <- rowMeans(msd_mat, na.rm = TRUE)
msd_ci <- apply(msd_mat, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
msd_df <- data.frame(
  step = seq_along(msd_mean),
  obs = final_diag$stats$msd_obs_curve,
  mean = msd_mean,
  lo = msd_ci[1, ],
  hi = msd_ci[2, ]
)

p_msd <- ggplot(msd_df, aes(x = step)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "firebrick", alpha = 0.2) +
  geom_line(aes(y = mean), color = "firebrick", linetype = "dashed") +
  geom_line(aes(y = obs), linewidth = 0.8) +
  labs(
    title = "Scenario 7 - MSD",
    subtitle = sprintf("MSD (p=%.3f)", final_diag$scorecard$P_Value[2]),
    x = "Lag",
    y = "MSD"
  ) +
  theme_minimal()

p_sinuosity <- ggplot(data.frame(sin = final_diag$stats$sin_sims), aes(x = sin)) +
  geom_histogram(fill = "grey80", color = "white", bins = 20) +
  geom_vline(xintercept = final_diag$stats$sin_obs, color = "red", linewidth = 1.2) +
  labs(
    title = "Scenario 7 - Sinuosity",
    subtitle = sprintf("Observed SI=%.4f (p=%.3f)", final_diag$scorecard$Observed[3], final_diag$scorecard$P_Value[3]),
    x = "Straightness Index"
  ) +
  theme_minimal()


ggsave(
  filename = file.path(out_dir, "scenario06_trajectories.png"),
  plot = p_ud,
  width = 10,
  height = 8,
  dpi = 120
)
ggsave(
  filename = file.path(out_dir, "scenario06_msd.png"),
  plot = p_msd,
  width = 10,
  height = 8,
  dpi = 120
)
ggsave(
  filename = file.path(out_dir, "scenario06_sinuosity.png"),
  plot = p_sinuosity,
  width = 10,
  height = 8,
  dpi = 120
)

cat(sprintf("Results written in: %s\n", normalizePath(out_dir)))
cat(sprintf("Screening table: %s\n", file.path(out_dir, "scenario06_screening.csv")))
cat(sprintf("Params/p-values: %s\n", file.path(out_dir, "scenario06_params_and_pvals.csv")))
