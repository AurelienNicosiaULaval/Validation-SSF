#' Robust Generative Validation for Step-Selection Functions
#'
#' Performs multi-criteria generative validation of an Integrated Step-Selection Analysis (iSSA) model.
#' Compares simulated trajectories to the observed track across four dimensions:
#' Topology (Wasserstein), Dynamics (MSD), Behavior (Sinuosity), and Connectivity (Barrier crossing).
#' Propagates parameter uncertainty by sampling coefficients from the multivariate normal distribution.
#'
#' @param model A fitted iSSA model object (typically from `amt::fit_issf`).
#' @param obs_track A `track_xyt` object from the `amt` package.
#' @param map A `SpatRaster` object representing the environmental covariate.
#' @param env_covar_name Character string. Name of the environmental covariate in the model.
#' @param n_sims Integer. Number of simulations (default: 50).
#' @param max_steps Integer. Number of steps to simulate. Default: length of observed track.
#' @param n_cands Integer. Number of candidate steps per step (default: 25).
#' @param barrier Optional. List (axis, val) or sf object defining a barrier for connectivity testing.
#' @param custom_mov_params Optional list of movement parameters to override fitted ones.
#'
#' @return An object of class `issf_diag`.
#' @export
diagnose_issf_robust <- function(
  model,
  obs_track,
  map,
  env_covar_name,
  n_sims = 50,
  max_steps = NULL,
  n_cands = 25,
  barrier = NULL,
  custom_mov_params = NULL
) {
  # --- 1. SETUP & PARAMETER SAMPLING ---
  if (is.null(max_steps)) {
    max_steps <- nrow(obs_track) - 1
  }
  clean_steps <- obs_track |> amt::steps() |> stats::na.omit()

  # Fit movement distributions
  if (is.null(custom_mov_params)) {
    fit_sl <- amt::fit_distr(clean_steps$sl_, "gamma")
    fit_ta <- amt::fit_distr(clean_steps$ta_, "vonmises")

    mov_params <- list(
      shape = fit_sl$params$shape,
      scale = fit_sl$params$scale,
      kappa = fit_ta$params$kappa,
      mu = fit_ta$params$mu
    )
  } else {
    mov_params <- custom_mov_params
  }

  # --- Beta Sampling (Uncertainty Propagation) ---
  if (inherits(model, "fit_clogit")) {
    raw_model <- model$model
  } else {
    raw_model <- model
  }

  beta_means <- stats::coef(raw_model)
  beta_vcov <- stats::vcov(raw_model)

  beta_samples <- mvtnorm::rmvnorm(
    n = n_sims,
    mean = beta_means,
    sigma = beta_vcov
  )
  colnames(beta_samples) <- names(beta_means)

  # --- 2. SIMULATION ---
  x0 <- obs_track$x_[1]
  y0 <- obs_track$y_[1]
  angle0 <- atan2(obs_track$y_[2] - y0, obs_track$x_[2] - x0)

  run_sim_iter <- function(i) {
    current_betas <- beta_samples[i, ]

    # Reconstruct coefficients list dynamically
    coeffs <- list(
      env = ifelse(
        env_covar_name %in% names(current_betas),
        current_betas[env_covar_name],
        0
      ),
      sl = ifelse("sl_" %in% names(current_betas), current_betas["sl_"], 0),
      log_sl = ifelse(
        "log_sl_" %in% names(current_betas),
        current_betas["log_sl_"],
        0
      ),
      cos_ta = ifelse(
        "cos_ta_" %in% names(current_betas),
        current_betas["cos_ta_"],
        0
      )
    )

    px <- numeric(max_steps + 1)
    py <- numeric(max_steps + 1)
    px[1] <- x0
    py[1] <- y0
    curr_angle <- angle0

    for (t in 1:max_steps) {
      # Use the user-defined n_cands (default 25)

      # A. Draw from Movement Kernel
      r_sl <- stats::rgamma(
        n_cands,
        shape = mov_params$shape,
        scale = mov_params$scale
      )
      r_ta <- circular::rvonmises(
        n_cands,
        mu = mov_params$mu,
        kappa = mov_params$kappa
      )
      c_ang <- curr_angle + r_ta
      c_x <- px[t] + r_sl * cos(c_ang)
      c_y <- py[t] + r_sl * sin(c_ang)

      # B. Extract Environment
      env_vals <- terra::extract(map, cbind(c_x, c_y))[, 1]
      valid <- !is.na(env_vals)

      if (sum(valid) == 0) {
        break
      }

      # C. Calculate Selection Weights
      log_w <- (coeffs$env * env_vals) +
        (coeffs$sl * r_sl) +
        (coeffs$log_sl * log(r_sl)) +
        (coeffs$cos_ta * cos(r_ta))

      w <- exp(log_w - max(log_w, na.rm = TRUE))
      w[!valid] <- 0

      # D. Selection
      if (sum(w) > 0) {
        choice <- sample(1:n_cands, 1, prob = w)
      } else {
        choice <- sample(1:n_cands, 1)
      }

      px[t + 1] <- c_x[choice]
      py[t + 1] <- c_y[choice]
      curr_angle <- as.numeric(c_ang[choice])
    }

    return(data.frame(x_ = px[1:(t + 1)], y_ = py[1:(t + 1)]))
  }

  # Execute simulations
  sims <- lapply(1:n_sims, run_sim_iter)

  # --- 3. METRICS CALCULATION  ---

  calc_pval_right <- function(obs, sim_vec) {
    (sum(sim_vec >= obs) + 1) / (length(sim_vec) + 1)
  }
  calc_pval_left <- function(obs, sim_vec) {
    (sum(sim_vec <= obs) + 1) / (length(sim_vec) + 1)
  }

  # -- 3.1 Topology (Wasserstein) --
  calc_wass <- function(t1, t2) {
    d1 <- stats::na.omit(t1[, c("x_", "y_")])
    d2 <- stats::na.omit(t2[, c("x_", "y_")])
    if (nrow(d1) == 0 || nrow(d2) == 0) {
      return(NA)
    }
    w1 <- rep(1 / nrow(d1), nrow(d1))
    w2 <- rep(1 / nrow(d2), nrow(d2))
    transport::wasserstein(transport::wpp(d1, w1), transport::wpp(d2, w2))
  }

  null_wass <- sapply(1:n_sims, function(i) {
    other <- sample((1:n_sims)[-i], 1)
    calc_wass(sims[[i]], sims[[other]])
  })

  obs_wass_vec <- sapply(sims, function(s) calc_wass(obs_track, s))
  mean_obs_wass <- mean(obs_wass_vec, na.rm = TRUE)
  mean_null_wass <- mean(null_wass, na.rm = TRUE)
  pval_wass <- calc_pval_right(mean_obs_wass, null_wass)

  # -- 3.2 Dynamics (MSD) --
  calc_msd_vec <- function(trk) (trk$x_ - trk$x_[1])^2 + (trk$y_ - trk$y_[1])^2
  msd_obs <- calc_msd_vec(obs_track)

  max_len <- length(msd_obs)
  msd_sims_mat <- do.call(
    cbind,
    lapply(sims, function(s) {
      m <- calc_msd_vec(s)
      if (length(m) < max_len) {
        m <- c(m, rep(NA, max_len - length(m)))
      }
      if (length(m) > max_len) {
        m <- m[1:max_len]
      }
      return(m)
    })
  )

  ise_obs <- sum(
    abs(msd_obs - rowMeans(msd_sims_mat, na.rm = TRUE)),
    na.rm = TRUE
  )
  ise_null <- sapply(1:n_sims, function(i) {
    sum(
      abs(
        msd_sims_mat[, i] -
          rowMeans(msd_sims_mat[, -i, drop = FALSE], na.rm = TRUE)
      ),
      na.rm = TRUE
    )
  })
  pval_msd <- calc_pval_right(ise_obs, ise_null)

  # -- 3.3 Behavior (Sinuosity) --
  calc_sin <- function(trk) {
    d_euc <- sqrt(
      (utils::tail(trk$x_, 1) - utils::head(trk$x_, 1))^2 +
        (utils::tail(trk$y_, 1) - utils::head(trk$y_, 1))^2
    )
    dp <- sum(sqrt(diff(trk$x_)^2 + diff(trk$y_)^2))
    return(ifelse(dp == 0, 0, d_euc / dp))
  }
  sin_obs <- calc_sin(obs_track)
  sin_sims <- sapply(sims, calc_sin)
  pval_sin <- calc_pval_right(
    abs(sin_obs - mean(sin_sims, na.rm = TRUE)),
    abs(sin_sims - mean(sin_sims, na.rm = TRUE))
  )

  # -- 3.4 Connectivity (Barrier) --
  barrier_stats <- NULL
  pval_bc <- NA
  bc_obs <- NA
  if (!is.null(barrier)) {
    # Enhanced barrier crossing function - supports both simple lines and sf geometries
    count_cross <- function(trk, barrier_def) {
      # Check if barrier is an sf object (LINESTRING/MULTILINESTRING)
      if (inherits(barrier_def, "sf") || inherits(barrier_def, "sfc")) {
        # Convert to sfc if sf
        if (inherits(barrier_def, "sf")) {
          barrier_geom <- sf::st_geometry(barrier_def)
        } else {
          barrier_geom <- barrier_def
        }

        # Check geometry type
        geom_types <- unique(sf::st_geometry_type(barrier_geom))
        if (!all(geom_types %in% c("LINESTRING", "MULTILINESTRING"))) {
          stop(
            "Barrier sf object must contain only LINESTRING or MULTILINESTRING geometries"
          )
        }

        # Check CRS compatibility
        barrier_crs <- sf::st_crs(barrier_geom)
        if (is.na(barrier_crs)) {
          warning(
            "Barrier has no CRS defined. Assuming it matches trajectory coordinates."
          )
        }

        # Create trajectory segments as LINESTRING
        n <- nrow(trk)
        crossings <- 0

        for (i in 2:n) {
          # Create segment from step i-1 to i
          seg_coords <- matrix(
            c(trk$x_[i - 1], trk$y_[i - 1], trk$x_[i], trk$y_[i]),
            ncol = 2,
            byrow = TRUE
          )
          seg <- sf::st_linestring(seg_coords)
          seg_sfc <- sf::st_sfc(seg, crs = barrier_crs) # Use barrier CRS

          # Check intersection with barrier
          if (any(sf::st_intersects(seg_sfc, barrier_geom, sparse = FALSE))) {
            crossings <- crossings + 1
          }
        }
        return(crossings)
      } else if (
        is.list(barrier_def) && all(c("axis", "val") %in% names(barrier_def))
      ) {
        # Original simple axis-based crossing (backward compatible)
        coords <- if (barrier_def$axis == "x") trk$x_ else trk$y_
        return(sum(abs(diff(sign(coords - barrier_def$val)))) / 2)
      } else {
        stop(
          "Barrier must be either an sf object (LINESTRING/MULTILINESTRING) or a list with 'axis' and 'val' elements"
        )
      }
    }

    bc_obs <- count_cross(obs_track, barrier)
    bc_sims <- sapply(sims, function(s) count_cross(s, barrier))
    pval_bc <- calc_pval_left(bc_obs, bc_sims)
    barrier_stats <- list(
      obs = bc_obs,
      sims = bc_sims,
      pval = pval_bc,
      def = barrier
    )
  }

  # --- 4. OUTPUT ---
  scorecard <- data.frame(
    Pillar = c("1. Topology", "2. Dynamics", "3. Behavior"),
    Metric = c("Wasserstein", "MSD (ISE)", "Sinuosity"),
    Observed = c(mean_obs_wass, ise_obs, sin_obs),
    P_Value = c(pval_wass, pval_msd, pval_sin)
  )

  if (!is.null(barrier)) {
    scorecard <- rbind(
      scorecard,
      data.frame(
        Pillar = "4. Connectivity",
        Metric = "Barrier Crossings",
        Observed = bc_obs,
        P_Value = pval_bc
      )
    )
  }

  scorecard$Result <- ifelse(
    scorecard$P_Value < 0.05,
    "REJECTED",
    "NOT REJECTED"
  )

  output <- list(
    scorecard = scorecard,
    n_sims = n_sims,
    betas_used = beta_samples,
    stats = list(
      wass_null_dist = null_wass,
      wass_null_mean = mean_null_wass,
      wass_obs_dist = obs_wass_vec,
      wass_obs_mean = mean_obs_wass,
      wass_pval = pval_wass,
      msd_sims_mat = msd_sims_mat,
      msd_obs_curve = msd_obs,
      msd_ise_obs = ise_obs,
      msd_pval = pval_msd,
      sin_sims = sin_sims,
      sin_sim_mean = mean(sin_sims, na.rm = TRUE),
      sin_sd = stats::sd(sin_sims, na.rm = TRUE),
      sin_obs = sin_obs,
      sin_pval = pval_sin,
      barrier = barrier_stats
    ),
    data = list(sims = sims, obs_track = obs_track)
  )
  class(output) <- "issf_diag"
  return(output)
}


#' Plot Diagnostic Results
#'
#' Generates a dashboard of plots to visualize the validation results.
#'
#' @param x An object of class `issf_diag`.
#' @param ... Additional arguments (ignored).
#'
#' @export
plot.issf_diag <- function(x, ...) {
  # Unified aesthetic for all panels
  col_sim <- "#4C78A8"  # blue: simulated quantities
  col_obs <- "#D55E00"  # orange: observed quantity
  col_ref <- "#4A4A4A"  # barrier/reference guides
  theme_diag <- ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey90")
    )

  # --- PLOT 0: SPATIAL MAP (The Spaghetti Plot) ---
  sims_df <- dplyr::bind_rows(x$data$sims, .id = "sim_id")

  p_map <- ggplot2::ggplot() +
    # 1. Simulations (blue)
    ggplot2::geom_path(
      data = sims_df,
      ggplot2::aes(x = x_, y = y_, group = sim_id, color = "Simulated"),
      alpha = 0.5,
      linewidth = 0.2
    ) +
    # 2. Observation (orange)
    ggplot2::geom_path(
      data = x$data$obs_track,
      ggplot2::aes(x = x_, y = y_, color = "Observed"),
      linewidth = 0.8
    ) +
    ggplot2::labs(
      title = "0. Observed and simulated trajectories",
      subtitle = NULL,
      x = "x-coordinate",
      y = "y-coordinate"
    ) +
    ggplot2::scale_color_manual(
      values = c("Observed" = col_obs, "Simulated" = col_sim),
      name = NULL
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.spacing.x = ggplot2::unit(0.4, "cm"),
      legend.box.spacing = ggplot2::unit(0.1, "cm"),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 7),
      legend.key.size = ggplot2::unit(0.3, "cm")
    ) +
    theme_diag +
    ggplot2::coord_equal()

  # Add Barrier line if present
  if (!is.null(x$stats$barrier)) {
    b_def <- x$stats$barrier$def
    if (b_def$axis == "x") {
      p_map <- p_map +
        ggplot2::geom_vline(
          xintercept = b_def$val,
          color = col_ref,
          linetype = "dashed",
          linewidth = 1
        )
    } else {
      p_map <- p_map +
        ggplot2::geom_hline(
          yintercept = b_def$val,
          color = col_ref,
          linetype = "dashed",
          linewidth = 1
        )
    }
  }

  # --- PLOT 1: WASSERSTEIN ---
  df_wass <- dplyr::bind_rows(
    dplyr::tibble(D = x$stats$wass_null_dist, T = "Simulated"),
    dplyr::tibble(D = x$stats$wass_obs_dist, T = "Observed")
  )
  p1 <- ggplot2::ggplot(df_wass, ggplot2::aes(x = D, fill = T)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::geom_vline(
      xintercept = x$stats$wass_obs_mean,
      color = col_obs,
      linetype = "dashed"
    ) +
    ggplot2::scale_fill_manual(
      values = c("Observed" = col_obs, "Simulated" = col_sim),
      name = ""
    ) +
    ggplot2::labs(
      title = "1. Emergent utilisation distributions (Wasserstein)",
      subtitle = paste0("p = ", round(x$stats$wass_pval, 3)),
      x = "Wasserstein distance",
      y = "Density",
      fill = ""
    ) +
    theme_diag

  # --- PLOT 2: MSD ---
  msd_mean <- rowMeans(x$stats$msd_sims_mat, na.rm = TRUE)
  msd_bounds <- apply(
    x$stats$msd_sims_mat,
    1,
    stats::quantile,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )
  df_msd <- data.frame(
    T = 1:length(msd_mean),
    Obs = x$stats$msd_obs_curve,
    Mean = msd_mean,
    Lo = msd_bounds[1, ],
    Hi = msd_bounds[2, ]
  )

  p2 <- ggplot2::ggplot(df_msd, ggplot2::aes(x = T)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = Lo, ymax = Hi),
      fill = col_sim,
      alpha = 0.2
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = Mean),
      color = col_sim,
      linetype = "dashed"
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = Obs),
      color = col_obs,
      linewidth = 0.8
    ) +
    ggplot2::labs(
      title = "2. Mean squared displacement (MSD)",
      subtitle = paste0("p = ", round(x$stats$msd_pval, 3)),
      y = "Mean squared displacement",
      x = "Lag (time steps)"
    ) +
    theme_diag

  # --- PLOT 3: SINUOSITY ---
  df_sin <- data.frame(S = x$stats$sin_sims)
  p3 <- ggplot2::ggplot(df_sin, ggplot2::aes(x = S)) +
    ggplot2::geom_histogram(
      fill = col_sim,
      color = "white",
      bins = 15,
      alpha = 0.7
    ) +
    ggplot2::geom_vline(
      xintercept = x$stats$sin_obs,
      color = col_obs,
      linewidth = 1.5
    ) +
    ggplot2::labs(
      title = "3. Path sinuosity (SI)",
      subtitle = paste0("p = ", round(x$stats$sin_pval, 3)),
      x = "Straightness index",
      y = "Number of simulations"
    ) +
    theme_diag

  # --- LAYOUT ---
  if (!is.null(x$stats$barrier)) {
    df_bc <- data.frame(C = x$stats$barrier$sims)
    p4 <- ggplot2::ggplot(df_bc, ggplot2::aes(x = C)) +
      ggplot2::geom_histogram(
        fill = col_sim,
        color = "white",
        binwidth = 1,
        alpha = 0.7
      ) +
      ggplot2::geom_vline(
        xintercept = x$stats$barrier$obs,
        color = col_obs,
        linewidth = 1.5
      ) +
      ggplot2::labs(
        title = "4. Functional connectivity (Barrier crossings)",
        subtitle = paste0(
          "p = ",
          round(x$stats$barrier$pval, 3)
        ),
        x = "Number of barrier crossings",
        y = "Number of simulations"
      ) +
      theme_diag

    gridExtra::grid.arrange(
      p_map,
      gridExtra::arrangeGrob(p1, p2, p3, p4, ncol = 2),
      ncol = 2,
      widths = c(1.5, 2)
    )
  } else {
    gridExtra::grid.arrange(p_map, p1, p2, p3, ncol = 2)
  }
}

#' Summarize Validation Results
#'
#' Prints a formatted ASCII report of the multi-criteria validation.
#'
#' @param x An object of class `issf_diag`.
#' @param ... Additional arguments (ignored).
#'
#' @export
summary.issf_diag <- function(x, ...) {
  pass <- "\u2713" # Check mark
  fail <- "\u2717" # Cross mark

  cat("\n")
  cat("==============================================================\n")
  cat(sprintf(
    "   GENERATIVE iSSA VALIDATION REPORT (%d simulations)\n",
    x$n_sims
  ))
  cat("==============================================================\n\n")

  # --- PILLAR 1 ---
  cat("1. SPATIAL TOPOLOGY (Wasserstein Distance)\n")
  cat("------------------------------------------\n")
  obs_w <- x$stats$wass_obs_mean
  null_w <- x$stats$wass_null_mean
  pval_w <- x$stats$wass_pval
  ratio <- obs_w / null_w

  cat(sprintf("   Sim-Sim Mean (Null)    : %.2f\n", null_w))
  cat(sprintf("   Obs-Sim Mean (Error)   : %.2f\n", obs_w))
  cat(sprintf("   Error Ratio            : %.2f x\n", ratio))

  if (pval_w < 0.05) {
    cat(sprintf("   %s RESULT: REJECTED (p = %.3f)\n", fail, pval_w))
    cat("      -> Significant spatial deviation detected.\n")
  } else {
    cat(sprintf("   %s RESULT: NOT REJECTED (p = %.3f)\n", pass, pval_w))
    cat("      -> Spatial structure consistent with reality.\n")
  }
  cat("\n")

  # --- PILLAR 2 ---
  cat("2. TEMPORAL DYNAMICS (Mean Squared Displacement)\n")
  cat("------------------------------------------------\n")
  pval_m <- x$stats$msd_pval
  cat(sprintf("   Integrated Squared Error (Obs): %.0f\n", x$stats$msd_ise_obs))

  if (pval_m < 0.05) {
    cat(sprintf("   %s RESULT: REJECTED (p = %.3f)\n", fail, pval_m))
    cat(
      "      -> Diffusion dynamics do not match (e.g. unbounded vs home range).\n"
    )
  } else {
    cat(sprintf("   %s RESULT: NOT REJECTED (p = %.3f)\n", pass, pval_m))
    cat("      -> Diffusion rate matches the observed trajectory.\n")
  }
  cat("\n")

  # --- PILLAR 3 ---
  cat("3. BEHAVIORAL PATTERN (Sinuosity)\n")
  cat("---------------------------------\n")
  obs_s <- x$stats$sin_obs
  sim_s <- x$stats$sin_sim_mean
  pval_s <- x$stats$sin_pval
  sd_s <- x$stats$sin_sd

  cat(sprintf("   Observed Index         : %.3f\n", obs_s))
  cat(sprintf("   Simulated Mean         : %.3f (SD: %.3f)\n", sim_s, sd_s))

  if (pval_s < 0.05) {
    cat(sprintf("   %s RESULT: REJECTED (p = %.3f)\n", fail, pval_s))
    cat("      -> Behavioral phenotype (tortuosity) mismatch.\n")
  } else {
    cat(sprintf("   %s RESULT: NOT REJECTED (p = %.3f)\n", pass, pval_s))
    cat("      -> Behavioral phenotype matches observations.\n")
  }

  # --- PILLAR 4 (Optional) ---
  if (!is.null(x$stats$barrier)) {
    cat("\n")
    cat("4. CONNECTIVITY (Barrier Crossing)\n")
    cat("----------------------------------\n")
    obs_b <- x$stats$barrier$obs
    sim_b <- mean(x$stats$barrier$sims)
    pval_b <- x$stats$barrier$pval

    cat(sprintf("   Observed Crossings     : %d\n", obs_b))
    cat(sprintf("   Simulated Mean         : %.1f\n", sim_b))

    if (pval_b < 0.05) {
      cat(sprintf("   %s RESULT: REJECTED (p = %.3f)\n", fail, pval_b))
      cat(
        "      -> Permeability mismatch (Simulations cross significantly more).\n"
      )
    } else {
      cat(sprintf("   %s RESULT: NOT REJECTED (p = %.3f)\n", pass, pval_b))
      cat("      -> Barrier permeability is consistent.\n")
    }
  }

  cat("\n==============================================================\n")
  cat(" Usage: plot(x) to visualize results.\n")
  cat("==============================================================\n")
}

#' Print Scorecard
#'
#' Displays the summary table of the validation results.
#'
#' @param x An object of class `issf_diag`.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.issf_diag <- function(x, ...) {
  cat("\n")
  print(knitr::kable(
    x$scorecard,
    digits = 3,
    caption = paste0("iSSA FORMAL VALIDATION TEST (", x$n_sims, " sims)")
  ))
  cat("\n")
}
