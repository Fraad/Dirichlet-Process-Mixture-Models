# ============================================================================
# Simulation study: KDE (Sheather-Jones) vs DPM (Stan stick-breaking) vs
# CRP (NIMBLE collapsed Gibbs) density estimation
#
# Densities:
#   (1) Normal:     N(0, 1)
#   (2) Bimodal:    0.5 N(-0.5, 0.25^2) + 0.5 N(0.5, 0.25^2)
#   (3) Multimodal: 4 peaks (small, large, small, medium)
#
# Sample sizes:   n in {25, 100, 200}
# Replications:   R = 20
# Eval grid:      G = 401 points on [-6, 6]  (odd -> Simpson's 1/3 rule)
#
# Priors:
#   sigma^2 (NIMBLE) / sigma (Stan) ~ InvGamma(a0 = 2, b0 = 1.0)
#   mu | sigma, tau ~ N(mu_0 = 0, sd = sigma * tau)        [as in existing files]
#   tau              ~ InvGamma(tau_a = 5, tau_b = 12)
#   alpha            ~ Gamma(shape = 2, rate = alpha_b)
#
# Alpha calibration: for each (density, n), alpha_b is set so that
#   E[alpha] * log(1 + n / E[alpha]) = K_target
# with K_target = (1, 2, 4) for (Normal, Bimodal, Multimodal). This gives
# 9 distinct (density, n) prior settings.
#
# Precompilation:
#   * Stan model compiled ONCE globally (alpha_a, alpha_b are data inputs).
#   * NIMBLE: 9 models compiled upfront (one per (density, n) cell), each
#     reused across the 20 replicates of its cell via cModel$setData(...).
#
# Outputs (in /mnt/user-data/outputs):
#   sim_density_results.rds   - full raw output (pointwise estimates + metrics)
#   sim_density_metrics.csv   - summary metrics table
# ============================================================================

suppressPackageStartupMessages({
  library(rstan)
  library(nimble)
})

set.seed(2026)
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

# ----------------------------------------------------------------------------
# 1. Settings
# ----------------------------------------------------------------------------

SMOKE_TEST <- FALSE   # set TRUE for a quick sanity check (tiny MCMC, 2 reps)

DENSITIES    <- c("Normal", "Bimodal", "Multimodal")
SAMPLE_SIZES <- c(25, 100, 200)
R_REPS       <- 15

G        <- 401                                    # odd -> Simpson works
GRID_LO  <- -6
GRID_HI  <-  6
eval_grid <- seq(GRID_LO, GRID_HI, length.out = G)
dx        <- (GRID_HI - GRID_LO) / (G - 1)

# Target prior expected number of clusters
K_TARGET <- c(Normal = 1, Bimodal = 2, Multimodal = 4)

# NIG hyperparameters (fixed across all densities, sample sizes, replications)
M_0     <- 0
A_0     <- 2
B_0     <- 1.0
TAU_A   <- 5
TAU_B   <- 12

TAU <- 3

ALPHA_A <- 2

K_TRUNC <- 15                                       # stick-breaking truncation

# MCMC settings
STAN_CHAINS <- 1
STAN_WARMUP <- 500
STAN_ITER   <- 1500        # 1000 post-warmup per chain -> 4000 total draws
NIM_CHAINS  <- 1
NIM_NITER   <- 1500
NIM_NBURN   <- 500

if (SMOKE_TEST) {
  R_REPS      <- 2
  SAMPLE_SIZES <- c(50)
  STAN_WARMUP <- 200; STAN_ITER <- 500
  NIM_NBURN   <- 200; NIM_NITER <- 500
}

STAN_FILE <- "Finite mixture.stan"
OUT_DIR   <- "Sim_outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# 2. Test densities
# ----------------------------------------------------------------------------

f_normal <- function(x) dnorm(x, 0, 1)
s_normal <- function(n) rnorm(n, 0, 1)

f_bimodal <- function(x) 0.5 * dnorm(x, -0.5, 0.25) + 0.5 * dnorm(x, 0.5, 0.25)
s_bimodal <- function(n) {
  z <- rbinom(n, 1, 0.5)
  ifelse(z == 1, rnorm(n, -0.5, 0.25), rnorm(n, 0.5, 0.25))
}

# 4 peaks (left -> right: small, large, small, medium) with equal sigma,
# so peak heights are proportional to weights.
mm_w  <- c(0.10, 0.45, 0.10, 0.35)
mm_mu <- c(-3, -1,  1,  3) 
mm_sd <- c(0.4, 0.4, 0.4, 0.4)
f_multimodal <- function(x) {
  out <- numeric(length(x))
  for (j in seq_along(mm_w)) {
    out <- out + mm_w[j] * dnorm(x, mm_mu[j], mm_sd[j])
  }
  out
}
s_multimodal <- function(n) {
  comp <- sample(seq_along(mm_w), n, replace = TRUE, prob = mm_w)
  rnorm(n, mm_mu[comp], mm_sd[comp])
}


dens_list <- list(
  Normal     = list(f = f_normal,     s = s_normal),
  Bimodal    = list(f = f_bimodal,    s = s_bimodal),
  Multimodal = list(f = f_multimodal, s = s_multimodal)
)

f_true_grid <- lapply(dens_list, function(d) d$f(eval_grid))

# ----------------------------------------------------------------------------
# 3. Composite Simpson's 1/3 rule
# ----------------------------------------------------------------------------
# For G grid points (G odd) on [a, b] with spacing dx = (b-a)/(G-1),
#   integral ~ (dx/3) * [f(x_1) + 4 f(x_2) + 2 f(x_3) + 4 f(x_4) + ... +
#                       2 f(x_{G-2}) + 4 f(x_{G-1}) + f(x_G)]

simpson <- function(fvals, dx) {
  Gloc <- length(fvals)
  if (Gloc %% 2 != 1L)
    stop("Simpson's rule requires an odd number of grid points")
  w <- rep(2, Gloc)
  w[c(1L, Gloc)] <- 1
  w[seq(2L, Gloc - 1L, by = 2L)] <- 4
  (dx / 3) * sum(w * fvals)
}

# ----------------------------------------------------------------------------
# 4. Alpha-prior calibration
# ----------------------------------------------------------------------------
# Solve  alpha_star * log(1 + n / alpha_star) = K_target   for alpha_star.
# Then alpha ~ Gamma(shape = 2, rate = alpha_b) with alpha_b = 2 / alpha_star
# gives E[alpha] = alpha_star.

solve_alpha_star <- function(n, K_target) {
  fn <- function(a) a * log1p(n / a) - K_target
  uniroot(fn, interval = c(1e-6, 1e6), tol = 1e-8)$root
}

alpha_grid <- expand.grid(
  density = DENSITIES, n = SAMPLE_SIZES,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
alpha_grid$K_target   <- K_TARGET[alpha_grid$density]
alpha_grid$alpha_star <- mapply(solve_alpha_star, alpha_grid$n, alpha_grid$K_target)
alpha_grid$alpha_b    <- ALPHA_A / alpha_grid$alpha_star

cat("\n--- Alpha-prior calibration  (Gamma(2, alpha_b), E[alpha] = alpha_star) ---\n")
print(alpha_grid, row.names = FALSE)

# Helper for cell lookup
get_alpha_b <- function(d, n) {
  alpha_grid$alpha_b[alpha_grid$density == d & alpha_grid$n == n]
}

# ----------------------------------------------------------------------------
# 5. Stan model: compile ONCE
# ----------------------------------------------------------------------------

cat("\n--- Compiling Stan model ---\n")
stan_code_text <- paste(readLines(STAN_FILE), collapse = "\n")
t_b <- Sys.time()
stan_model_obj <- stan_model(model_code = stan_code_text,
                             model_name = "dpm_stickbreaking")
stan_compile_time <- as.numeric(difftime(Sys.time(), t_b, units = "secs"))
cat(sprintf("Stan compile time: %.1f sec\n", stan_compile_time))

# ----------------------------------------------------------------------------
# 6. NIMBLE: precompile 9 models (3 densities x 3 sample sizes)
# ----------------------------------------------------------------------------

crp_code <- nimbleCode({
  for (i in 1:n) {
    y[i]   ~ dnorm(mu[i], sd = s[i] )
    mu[i] <- muTilde[xi[i]]
    s[i] <- sTilde[xi[i]]
  }
  xi[1:n] ~ dCRP(alpha, size = n)
  for (i in 1:n) {
    muTilde[i] ~ dnorm(mu_0, sd =  sqrt(tau) * sTilde[i] )
    sTilde[i] ~ dinvgamma(a_0, b_0)
  }
  #tau   ~ dinvgamma(tau_a, tau_b)
  alpha ~ dgamma(shape = alpha_a, rate = alpha_b)
})

build_nimble_cell <- function(nn, alpha_b_val, dens_fn) {
  y0 <- dens_fn(nn)                                # placeholder data
  consts <- list(
    n = nn,
    mu_0 = M_0, a_0 = A_0, b_0 = B_0,
    #tau_a = TAU_A, tau_b = TAU_B,
    tau = TAU,
    alpha_a = ALPHA_A, alpha_b = alpha_b_val
  )
  inits <- list(
    xi      = sample.int(min(nn, 5), size = nn, replace = TRUE),
    muTilde = rnorm(nn, M_0, 1),
    sTilde = nimble::rinvgamma(nn, A_0, B_0),
    #tau     = max(nimble::rinvgamma(1, TAU_A, TAU_B), 0.1),
    alpha   = max(ALPHA_A / alpha_b_val, 0.1)
  )
  rModel <- nimbleModel(crp_code,
                        data = list(y = y0),
                        constants = consts,
                        inits = inits,
                        check = FALSE, calculate = FALSE)
  cModel <- compileNimble(rModel, showCompilerOutput = FALSE)
  conf   <- configureMCMC(rModel,
                          monitors = c("xi", "muTilde", "sTilde", "alpha", "tau"),
                          print = FALSE)
  mcmcR  <- buildMCMC(conf)
  cMCMC  <- compileNimble(mcmcR, project = rModel,
                          resetFunctions = TRUE, showCompilerOutput = FALSE)
  list(cModel = cModel, cMCMC = cMCMC)
}

cat("\n--- Precompiling NIMBLE (9 models) ---\n")
t_nim <- Sys.time()
nimble_models <- list()
for (d in DENSITIES) {
  nimble_models[[d]] <- list()
  for (nn in SAMPLE_SIZES) {
    ab <- get_alpha_b(d, nn)
    t_b <- Sys.time()
    cat(sprintf("  density=%-11s n=%3d  alpha_b=%7.4f  alpha_star=%7.4f  ... ",
                d, nn, ab, ALPHA_A / ab))
    nimble_models[[d]][[as.character(nn)]] <-
      build_nimble_cell(nn, ab, dens_list[[d]]$s)
    cat(sprintf("done (%.1f s)\n",
                as.numeric(difftime(Sys.time(), t_b, units = "secs"))))
  }
}
nimble_compile_time <- as.numeric(difftime(Sys.time(), t_nim, units = "secs"))
cat(sprintf("Total NIMBLE compile time: %.1f sec\n", nimble_compile_time))

# ----------------------------------------------------------------------------
# 7. Method fits: KDE, DPM (Stan), CRP (NIMBLE)
# ----------------------------------------------------------------------------

# (a) KDE with Sheather-Jones bandwidth, refit per replication.
fit_kde <- function(y, grid) {
  bw <- tryCatch(bw.SJ(y), error = function(e) bw.nrd0(y))
  d  <- density(y, bw = bw, from = GRID_LO, to = GRID_HI, n = 1024)
  approx(d$x, d$y, xout = grid, rule = 2)$y
}

# (b) Stick-breaking DPM. Posterior mean of predictive density
#       f_hat(x) = sum_k pi_k * N(x | mu_k, sigma_k)
fit_dpm <- function(y, alpha_b_val, grid) {
  stan_data <- list(
    N = length(y),
    K = K_TRUNC,
    y = y,
    sigma_a = A_0, sigma_b = B_0,
    mu_0    = M_0,
    #tau_a   = TAU_A, tau_b = TAU_B,
    tau = TAU,
    alpha_a = ALPHA_A, alpha_b = alpha_b_val
  )
  fit <- 
    sampling(stan_model_obj,
             data    = stan_data,
             chains  = STAN_CHAINS,
             warmup  = STAN_WARMUP,
             iter    = STAN_ITER,
             refresh = 0,
             control = list(adapt_delta = 0.9, max_treedepth = 12)
  )
  ex      <- rstan::extract(fit, pars = c("mu", "sigma", "pi_var"))
  mu_s    <- ex$mu        # [draws, K]
  sigma_s <- ex$sigma     # [draws, K]
  pi_s    <- ex$pi_var    # [draws, K]
  n_draws <- nrow(mu_s)
  f_acc <- numeric(length(grid))
  for (r in seq_len(n_draws)) {
    for (k in seq_len(K_TRUNC)) {
      f_acc <- f_acc + pi_s[r, k] * dnorm(grid, mu_s[r, k], sigma_s[r, k])
    }
  }
  f_acc / n_draws
}

# (c) CRP collapsed Gibbs. Posterior mean of predictive density
#       f_hat(x | theta) = [ sum_{k in occ} m_k * N(x | mu*_k, sigma*_k)
#                          + alpha * N(x | mu_new, sigma_new) ] / (alpha + n)
fit_crp <- function(y, dens_label, nn, grid) {
  nm <- nimble_models[[dens_label]][[as.character(nn)]]
  nm$cModel$setData(list(y = y))
  inits_fn <- function() {
    list(
      xi      = sample.int(min(nn, 5), size = nn, replace = TRUE),
      muTilde = rnorm(nn, M_0, 1),
      sTilde = nimble::rinvgamma(nn, A_0, B_0),
      #tau     = max(nimble::rinvgamma(1, TAU_A, TAU_B), 0.1),
      alpha   = abs(rnorm(1, 1, 0.5)) + 0.1
    )
  }
  samples_list <- runMCMC(
    nm$cMCMC,
    niter   = NIM_NITER,
    nburnin = NIM_NBURN,
    nchains = NIM_CHAINS,
    inits   = inits_fn,
    samplesAsCodaMCMC = FALSE,
    progressBar = FALSE,
    setSeed = FALSE
  )
  samp_mat <- if (NIM_CHAINS > 1) do.call(rbind, samples_list) else samples_list

  alpha_col   <- samp_mat[, "alpha"]
  #tau_col     <- samp_mat[, "tau"]
  muT_cols    <- samp_mat[, grep("^muTilde\\[",  colnames(samp_mat))]
  sT_cols    <- samp_mat[, grep("^sTilde\\[",  colnames(samp_mat))]
  xi_cols     <- samp_mat[, grep("^xi\\[",       colnames(samp_mat))]

  n_iter <- nrow(samp_mat)
  f_acc  <- numeric(length(grid))
  for (r in seq_len(n_iter)) {
    xi_r       <- xi_cols[r, ]
    k_unique   <- unique(xi_r)
    cnt_full   <- tabulate(xi_r, nbins = nn)
    m_k        <- cnt_full[k_unique]
    unused     <- setdiff(seq_len(nn), k_unique)

    if (length(unused) > 0L) {
      kNew   <- unused[1L]
      mu_new <- muT_cols[r, kNew]
      s_new <- sT_cols[r, kNew]
    } else {
      # All n slots occupied (rare): draw a fresh sample from the base measure
      s_new <- nimble::rinvgamma(1, A_0, B_0)
      mu_new <- rnorm(1, M_0,  sqrt(TAU) * s_new) # OLD: tau_col[r] * sqrt(s2_new)
    }

    occ_term <- numeric(length(grid))
    for (k in seq_along(k_unique)) {
      occ_term <- occ_term +
        m_k[k] * dnorm(grid, muT_cols[r, k_unique[k]],
                             sqrt(sT_cols[r, k_unique[k]]))
    }
    new_term <- alpha_col[r] * dnorm(grid, mu_new, s_new)
    f_acc <- f_acc + (occ_term + new_term) / (alpha_col[r] + nn)
  }
  f_acc / n_iter
}

# ----------------------------------------------------------------------------
# 8. Main simulation loop
# ----------------------------------------------------------------------------

results    <- list()
total_runs <- length(DENSITIES) * length(SAMPLE_SIZES) * R_REPS
run_count  <- 0L
t_sim_start <- Sys.time()

set.seed(2026)

cat("\n--- Running simulation ", total_runs, " total replicates ---\n", sep = "")

for (d in DENSITIES) {
  results[[d]] <- list()
  for (nn in SAMPLE_SIZES) {
    ab <- get_alpha_b(d, nn)

    f_kde_mat <- matrix(NA_real_, R_REPS, G)
    f_dpm_mat <- matrix(NA_real_, R_REPS, G)
    f_crp_mat <- matrix(NA_real_, R_REPS, G)
    rt_kde <- numeric(R_REPS)
    rt_dpm <- numeric(R_REPS)
    rt_crp <- numeric(R_REPS)

    for (r in seq_len(R_REPS)) {
      run_count <- run_count + 1L
      y <- dens_list[[d]]$s(nn)

      t0 <- Sys.time()
      f_kde_mat[r, ] <- fit_kde(y, eval_grid)
      rt_kde[r] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

      t0 <- Sys.time()
      f_dpm_mat[r, ] <- fit_dpm(y, ab, eval_grid)
      rt_dpm[r] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

      t0 <- Sys.time()
      f_crp_mat[r, ] <- fit_crp(y, d, nn, eval_grid)
      rt_crp[r] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

      cat(sprintf("[%3d/%d] %-11s n=%3d rep=%2d/%d | KDE %.2fs  DPM %.1fs  CRP %.1fs\n",
                  run_count, total_runs, d, nn, r, R_REPS,
                  rt_kde[r], rt_dpm[r], rt_crp[r]))
    }

    results[[d]][[as.character(nn)]] <- list(
      KDE = list(f = f_kde_mat, runtime = rt_kde),
      DPM = list(f = f_dpm_mat, runtime = rt_dpm),
      CRP = list(f = f_crp_mat, runtime = rt_crp)
    )
  }
}
sim_total_time_min <- as.numeric(difftime(Sys.time(), t_sim_start, units = "mins"))
cat(sprintf("\nTotal simulation runtime: %.1f min\n", sim_total_time_min))

# ----------------------------------------------------------------------------
# 9. Metrics: MISE = IBias^2 + IVar  (1/R variance for exact decomposition),
#             MAE = average over replicates of integrated |f_hat - f|,
#             plus average runtime per replication.
# ----------------------------------------------------------------------------

compute_metrics <- function(f_mat, f_grid_true, dx) {
  Rloc <- nrow(f_mat)
  Gloc <- ncol(f_mat)
  f_bar  <- colMeans(f_mat)
  bias_sq_grid <- (f_bar - f_grid_true)^2
  # 1/R denominator -> bias_sq + var = mse exactly at each grid point
  var_grid <- colMeans(sweep(f_mat, 2, f_bar, FUN = "-")^2)
  mse_grid <- colMeans(sweep(f_mat, 2, f_grid_true, FUN = "-")^2)
  l1_per_rep <- apply(abs(sweep(f_mat, 2, f_grid_true, FUN = "-")),
                      1, simpson, dx = dx)
  list(
    MISE   = simpson(mse_grid,     dx),
    IBias2 = simpson(bias_sq_grid, dx),
    IVar   = simpson(var_grid,     dx),
    "%Ivar" = simpson(var_grid, dx) /simpson(mse_grid,     dx),
    MAE    = mean(l1_per_rep)
  )
}

rows <- list()
for (d in DENSITIES) {
  f_grid <- f_true_grid[[d]]
  for (nn in SAMPLE_SIZES) {
    cell <- results[[d]][[as.character(nn)]]
    for (m in c("KDE", "DPM", "CRP")) {
      met <- compute_metrics(cell[[m]]$f, f_grid, dx)
      rows[[length(rows) + 1L]] <- data.frame(
        density          = d,
        n                = nn,
        method           = m,
        MISE             = met$MISE,
        IBias2           = met$IBias2,
        IVar             = met$IVar,
        "%Ivar"          = met$`%Ivar`,
        MAE              = met$MAE,
        avg_runtime_sec  = mean(cell[[m]]$runtime),
        stringsAsFactors = FALSE
      )
    }
  }
}
results_table <- do.call(rbind, rows)
cat("\n--- Summary table ---\n")
print(results_table, row.names = FALSE, digits = 2)



#9b. Replicate-curve plots: 3x3 (method x sample size) panels, one per density


library(ggplot2)

library(patchwork)
# For each (density, method, n) cell we plot all R_REPS density estimates as

# thin grey curves with the true density overlaid in bold black -- in the

# spirit of the reference plot showing how the estimator varies over the

# truth across replicates.



METHODS <- c("KDE", "DPM", "CRP")



# Per-density x-axis windows: restrict to the support that contains the

# action without clipping the curves themselves (coord_cartesian, not xlim).

plot_xlim <- list(
  
  Normal     = c(-3.5,  3.5),
  
  Bimodal    = c(-2.0,  2.0),
  
  Multimodal = c(-5.0,  5.0)
  
)



build_long_df <- function(d) {
  
  out <- list()
  
  for (m in METHODS) {
    
    for (nn in SAMPLE_SIZES) {
      
      f_mat <- results[[d]][[as.character(nn)]][[m]]$f   # R_REPS x G
      
      df <- data.frame(
        
        x      = rep(eval_grid, times = nrow(f_mat)),
        
        f      = as.vector(t(f_mat)),
        
        rep    = rep(seq_len(nrow(f_mat)), each = length(eval_grid)),
        
        method = m,
        
        n      = nn
        
      )
      
      out[[length(out) + 1L]] <- df
      
    }
    
  }
  
  do.call(rbind, out)
  
}



build_truth_df <- function(d) {
  
  expand.grid(method = METHODS, n = SAMPLE_SIZES,
              
              stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE) |>
    
    (\(g) {
      
      do.call(rbind, lapply(seq_len(nrow(g)), function(i) {
        
        data.frame(x = eval_grid, f_true = f_true_grid[[d]],
                   
                   method = g$method[i], n = g$n[i])
        
      }))
      
    })()
  
}



# Per-cell mean estimate (average over the R_REPS replicate curves).

build_mean_df <- function(d) {
  
  out <- list()
  
  for (m in METHODS) {
    
    for (nn in SAMPLE_SIZES) {
      
      f_mat <- results[[d]][[as.character(nn)]][[m]]$f       # R_REPS x G
      
      f_bar <- colMeans(f_mat)
      
      out[[length(out) + 1L]] <- data.frame(
        
        x = eval_grid, f_mean = f_bar, method = m, n = nn
        
      )
      
    }
    
  }
  
  do.call(rbind, out)
  
}



plot_density_grid <- function(d) {
  
  long_df  <- build_long_df(d)
  
  truth_df <- build_truth_df(d)
  
  mean_df  <- build_mean_df(d)
  
  
  
  # Order factors so rows go KDE -> DPM -> CRP and columns 25 -> 100 -> 200
  
  long_df$method  <- factor(long_df$method,  levels = METHODS)
  
  long_df$n       <- factor(long_df$n,       levels = SAMPLE_SIZES,
                            
                            labels = paste0("n = ", SAMPLE_SIZES))
  
  truth_df$method <- factor(truth_df$method, levels = METHODS)
  
  truth_df$n      <- factor(truth_df$n,      levels = SAMPLE_SIZES,
                            
                            labels = paste0("n = ", SAMPLE_SIZES))
  
  mean_df$method  <- factor(mean_df$method,  levels = METHODS)
  
  mean_df$n       <- factor(mean_df$n,       levels = SAMPLE_SIZES,
                            
                            labels = paste0("n = ", SAMPLE_SIZES))
  
  
  
  y_top <- max(c(long_df$f, truth_df$f_true, mean_df$f_mean),
               
               na.rm = TRUE) * 1.05
  
  
  
  ggplot() +
    
    geom_line(data = long_df,
              
              aes(x = x, y = f, group = rep),
              
              colour = "grey55", alpha = 0.45, linewidth = 0.3) +
    
    geom_line(data = truth_df,
              
              aes(x = x, y = f_true),
              
              colour = "black", linewidth = 0.9) +
    
    geom_line(data = mean_df,
              
              aes(x = x, y = f_mean),
              
              colour = "#1f6feb", linewidth = 0.8) +
    
    facet_grid(method ~ n, switch = "y") +
    
    coord_cartesian(xlim = plot_xlim[[d]], ylim = c(0, y_top),
                    
                    expand = FALSE) +
    
    labs(
      
      title    = paste0("Density estimates across replicates - ", d),
      
      subtitle = paste0("Grey: ", R_REPS,
                        
                        " replicate estimates.  Blue: mean estimate.  ",
                        
                        "Black: true density."),
      
      x = "x", y = "density"
      
    ) +
    
    theme_bw(base_size = 11) +
    
    theme(
      
      strip.background = element_rect(fill = "grey90", colour = NA),
      
      strip.text       = element_text(face = "bold"),
      
      panel.grid.minor = element_blank(),
      
      plot.title       = element_text(face = "bold")
      
    )
  
}



  cat("\n--- Building replicate-curve plots ---\n")
  
  for (d in DENSITIES) {
    
    p <- plot_density_grid(d)
    
    out_pdf <- file.path(OUT_DIR, sprintf("replicate_curves_%s.pdf", d))
    
    out_png <- file.path(OUT_DIR, sprintf("replicate_curves_%s.png", d))
    
    ggsave(out_pdf, p, width = 9, height = 8, units = "in")
    
    ggsave(out_png, p, width = 9, height = 8, units = "in", dpi = 200)
    
    cat(sprintf("  Wrote %s\n", out_pdf))
    
  }
  



# ----------------------------------------------------------------------------
# 10. Persist
# ----------------------------------------------------------------------------

saveRDS(
  list(
    results             = results,
    results_table       = results_table,
    alpha_grid          = alpha_grid,
    eval_grid           = eval_grid,
    dx                  = dx,
    f_true_grid         = f_true_grid,
    stan_compile_time   = stan_compile_time,
    nimble_compile_time = nimble_compile_time,
    sim_total_time_min  = sim_total_time_min,
    settings = list(
      DENSITIES = DENSITIES, SAMPLE_SIZES = SAMPLE_SIZES, R_REPS = R_REPS,
      G = G, GRID_LO = GRID_LO, GRID_HI = GRID_HI,
      M_0 = M_0, A_0 = A_0, B_0 = B_0,
      #TAU_A = TAU_A, TAU_B = TAU_B, 
      TAU = TAU,
      ALPHA_A = ALPHA_A, K_TRUNC = K_TRUNC,
      STAN_CHAINS = STAN_CHAINS, STAN_WARMUP = STAN_WARMUP, STAN_ITER = STAN_ITER,
      NIM_CHAINS  = NIM_CHAINS,  NIM_NITER   = NIM_NITER,  NIM_NBURN  = NIM_NBURN
    )
  ),
  file.path(OUT_DIR, "sim_density_results.rds")
)
write.csv(results_table,
          file.path(OUT_DIR, "sim_density_metrics.csv"),
          row.names = FALSE)
cat(sprintf("\nWrote:\n  %s\n  %s\n",
            file.path(OUT_DIR, "sim_density_results.rds"),
            file.path(OUT_DIR, "sim_density_metrics.csv")))





