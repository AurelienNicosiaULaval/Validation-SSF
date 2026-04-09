# ==============================================================================
# SENSITIVITY ANALYSIS: BARRIER CROSSING DIAGNOSTIC POWER
# ==============================================================================
# Description:
# Demonstrates how the statistical power of the barrier crossing diagnostic
# depends on trajectory length N and diffusivity (step length variance).
# ==============================================================================

library(amt)
library(tidyverse)
library(terra)
library(parallel)

source("R/diagnose_issf.R")
set.seed(888)

# Create figures directory
if (!dir.exists("manuscript/figures")) {
    dir.create("manuscript/figures", recursive = TRUE)
}

# ==============================================================================
# CONFIGURATION
# ==============================================================================
# Options: "SANITY" or "FULL"
SIMULATION_MODE <- "FULL"

if (SIMULATION_MODE == "SANITY") {
    message(">>> MODE: SANITY CHECK <<<")
    # Quick run to confirm power > 0 with adequate K
    n_replicates_per_grid <- 20
    n_sims_test <- 39 # K=39 -> p_min = 1/40 = 0.025 < 0.05

    # Single point that should have high power
    grid_N <- c(300)
    grid_sigma <- c(15)
} else if (SIMULATION_MODE == "FULL") {
    message(">>> MODE: FULL GRID SEARCH <<<")
    n_replicates_per_grid <- 50
    n_sims_test <- 99 # K=99 -> p_min = 0.01

    # Full grid
    grid_N <- c(100, 200, 300, 500, 800)
    grid_sigma <- c(5, 10, 15, 20, 30)
}

# Determine cores for parallel processing
# Use logical=FALSE to avoid oversubscription on hyperthreaded cores if desired,
# but usually safely using all cores minus 1 is good.
n_cores <- max(1, parallel::detectCores(logical = FALSE) - 1)
if (n_cores > 1) {
    message(sprintf("Using %d cores for parallel processing.", n_cores))
} else {
    message("Running sequentially.")
}

# ==============================================================================
# HELPERS
# ==============================================================================

# Helper: Generate trajectory with hard barrier
# Returns list with track and collision stats
simulate_barrier_trajectory <- function(
    N,
    sigma,
    start_x = 580,
    barrier_x = 595
) {
    t <- tibble(x = start_x, y = 500, t = as.POSIXct("2024-01-01"))
    collisions <- 0

    for (i in 1:N) {
        p <- t[nrow(t), ]

        # Simple biased random walk logic matching the description
        # x-bias of 5m, y-bias of 0
        dx <- rnorm(1, 5, sigma)
        dy <- rnorm(1, 0, sigma)

        nx <- p$x + dx

        # Hard barrier reflection
        if (nx > barrier_x) {
            # Collision!
            collisions <- collisions + 1
            nx <- barrier_x - (nx - barrier_x)
        }

        t <- bind_rows(t, tibble(x = nx, y = p$y + dy, t = p$t + 3600))
    }

    trk <- make_track(t, x, y, t, crs = 32601)
    return(list(trk = trk, attempted_crossings = collisions))
}

# Helper: Fit iSSA robustly (NAIVE MODEL - BLIND TO BARRIER)
fit_barrier_issa <- function(trk, r) {
    df <- trk |>
        steps() |>
        random_steps(n_control = 15) |>
        extract_covariates(r) |>
        filter(!is.na(sl_), !is.na(ta_), sl_ > 0.001) |>
        mutate(
            log_sl = log(sl_ + 0.01),
            sl_scaled = scale(sl_),
            cos_ta = cos(ta_)
        ) |>
        na.omit()

    tryCatch(
        {
            # ALWAYS fit the model WITHOUT 'env' to simulate model deficiency
            # The model assumes the landscape is homogeneous
            fit_issf(
                df,
                case_ ~ sl_scaled + log_sl + cos_ta + strata(step_id_)
            )
        },
        error = function(e) NULL
    )
}

# Helper: Run one test iteration
run_test_iteration <- function(
    scenario_name,
    start_x,
    N,
    sigma,
    n_sims_test
) {
    # Create landscape for context (barrier at 595-600)
    # We define "env" as forest/open split at 600 for the fitted model context if needed,
    # but here we fit a naive model that ignores it.

    r <- rast(
        nrows = 100,
        ncols = 100,
        xmin = 0,
        xmax = 1000,
        ymin = 0,
        ymax = 1000
    )
    values(r) <- 0
    names(r) <- "env"

    # Generate observed trajectory (Hard Barrier DGP)
    sim_res <- simulate_barrier_trajectory(
        N,
        sigma,
        start_x,
        595
    )
    trk_obs <- sim_res$trk
    attempted_crossings <- sim_res$attempted_crossings

    # Fit Naive iSSA (Blind to barrier)
    m <- fit_barrier_issa(trk_obs, r)
    if (is.null(m)) {
        return(NULL)
    }

    # Run diagnostic - Pillar 4 (Connectivity / Barrier)
    # Checks for crossings at x=595
    bstats <- tryCatch(
        {
            res <- diagnose_issf_robust(
                m,
                trk_obs,
                r,
                "env",
                n_sims = n_sims_test,
                barrier = list(axis = "x", val = 595)
            )
            res$stats$barrier
        },
        error = function(e) NULL
    )

    if (is.null(bstats)) {
        return(NULL)
    }

    # Return key metrics WITHOUT context columns (N, sigma) to avoid name clash
    return(tibble(
        attempted_crossings_dgp = attempted_crossings, # How many times it hit wall in DGP
        obs_crossings_diag = bstats$obs, # Should be 0 given reflection
        sim_cross_mean = mean(bstats$sims), # Mean crossings in simulated replicates
        sim_prop_any = mean(bstats$sims > 0), # Prop sims with >= 1 crossing
        p_value = bstats$pval
    ))
}


# ==============================================================================
# MAIN SIMULATION LOOP
# ==============================================================================
message(sprintf("Starting simulation in %s mode...", SIMULATION_MODE))
message(sprintf("N Replicates: %d", n_replicates_per_grid))
message(sprintf("K Sims/Test: %d", n_sims_test))

start_time <- Sys.time()

# Define Scenarios - User requested keeping Baseline
scenarios <- tribble(
    ~name      , ~start_x ,
    "Baseline" ,      580
)

# Build Experiment Grid
grid <- expand_grid(
    scenarios,
    N = grid_N,
    sigma = grid_sigma,
    replicate = 1:n_replicates_per_grid
)

message(sprintf("Total runs to process: %d", nrow(grid)))

# Function to process a single row index
process_row <- function(i) {
    row <- grid[i, ]

    # Run test
    res <- run_test_iteration(
        row$name,
        row$start_x,
        row$N,
        row$sigma,
        n_sims_test
    )

    if (is.null(res)) {
        return(NULL)
    }

    # Bind results with context
    bind_cols(row, res)
}

# Run using mclapply (Mac/Linux compatible)
# Windows users would need parLapply, but user is on Mac
results_list <- mclapply(
    1:nrow(grid),
    process_row,
    mc.cores = n_cores,
    mc.preschedule = TRUE # Better load balancing for uniform tasks? TRUE is default.
)

# Combine results, removing NULLs
results <- bind_rows(results_list)

end_time <- Sys.time()
runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
message(sprintf("Simulation completed in %.1f seconds", runtime))

# ==============================================================================
# ANALYSIS & OUTPUT
# ==============================================================================

if (nrow(results) == 0) {
    stop("No valid results collected!")
}

# 1. Summary
summary_stats <- results |>
    group_by(N, sigma) |>
    summarize(
        n_valid = n(),
        # Power: prop p < 0.05
        power_hat = mean(p_value < 0.05, na.rm = TRUE),
        # SE of Power
        power_se = sqrt((power_hat * (1 - power_hat)) / n_valid),

        mean_crossings_naive = mean(sim_cross_mean, na.rm = TRUE),
        prop_any_cross_naive = mean(sim_prop_any, na.rm = TRUE),
        .groups = "drop"
    )

print(summary_stats)

# 2. Save Results
if (!dir.exists("results")) {
    dir.create("results")
}
write_csv(summary_stats, "results/barrier_sensitivity_summary.csv")
write_csv(results, "results/barrier_sensitivity_details.csv")

# 3. Quick Plotting (if Full mode or just to check)
if (nrow(summary_stats) > 0) {
    # Power vs N
    p1 <- ggplot(
        summary_stats,
        aes(x = N, y = power_hat, color = factor(sigma))
    ) +
        geom_line() +
        geom_point() +
        geom_errorbar(
            aes(
                ymin = power_hat - 1.96 * power_se,
                ymax = power_hat + 1.96 * power_se
            ),
            width = 0.1
        ) +
        geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray") +
        geom_hline(yintercept = 0.80, linetype = "dashed", color = "green") +
        labs(
            title = "Power vs Trajectory Length (N)",
            y = "Power (Rate p < 0.05)",
            color = "Sigma (m)"
        ) +
        theme_minimal()

    ggsave(
        "manuscript/figures/barrier_sensitivity_N.png",
        p1,
        width = 6,
        height = 4,
        bg = "white"
    )

    # Power vs Sigma
    p2 <- ggplot(
        summary_stats,
        aes(x = sigma, y = power_hat, color = factor(N))
    ) +
        geom_line() +
        geom_point() +
        labs(
            title = "Power vs Diffusivity (Sigma)",
            y = "Power (Rate p < 0.05)",
            color = "N (steps)"
        ) +
        theme_minimal()

    ggsave(
        "manuscript/figures/barrier_sensitivity_sigma.png",
        p2,
        width = 6,
        height = 4,
        bg = "white"
    )
}

message("Done.")
