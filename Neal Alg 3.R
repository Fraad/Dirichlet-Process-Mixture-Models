




#Neal Algorithm 3 (NIG + kappa) ---------------------------------------

# NIG adaptation: predictive marginals are Student-t (not Gaussian);
# at end of sweep, instantiate (mu_c, s2_c) for each occupied cluster, then
# update kappa via Gibbs and alpha via Escobar-West.
neal_alg3_nig <- function(y, n_iter, n_burn, mu0,
                          a_sigma, b_sigma, a_kappa, b_kappa,
                          alpha_shape, alpha_rate) {
  n <- length(y)
  
  kappa   <- 1 / rgamma(1, a_kappa, b_kappa)
  alpha   <- rgamma(1, alpha_shape, alpha_rate)
  c_alloc <- rep(1L, n)
  n_k     <- as.integer(table(c_alloc))
  Nclust  <- length(n_k)
  
  # Log marginal predictive m(y_new | y_c, kappa)
  # NIG -> Student-t: df = 2*a_n, location = mu_n,
  #                   scale^2 = (b_n / a_n) * (1 + 1/lambda_n)
  log_marg <- function(y_new, y_c, kappa_cur) {
    if (length(y_c) == 0L) {
      df     <- 2 * a_sigma
      loc    <- mu0
      scale2 <- (b_sigma / a_sigma) * (1 + kappa_cur)
    } else {
      nk_c     <- length(y_c)
      ybar     <- mean(y_c)
      Sk       <- sum((y_c - ybar)^2)
      lambda_n <- 1 / kappa_cur + nk_c
      mu_n     <- (mu0 / kappa_cur + nk_c * ybar) / lambda_n
      a_n      <- a_sigma + nk_c / 2
      b_n      <- b_sigma + 0.5 * Sk +
        nk_c * (ybar - mu0)^2 / (2 * (1 + nk_c * kappa_cur))
      df     <- 2 * a_n
      loc    <- mu_n
      scale2 <- (b_n / a_n) * (1 + 1 / lambda_n)
    }
    z <- (y_new - loc) / sqrt(scale2)
    dt(z, df, log = TRUE) - 0.5 * log(scale2)
  }
  
  draws <- vector("list", n_iter)
  
  for (it in seq_len(n_iter)) {
    
    #Allocation step ---
    for (i in seq_len(n)) {
      c_i <- c_alloc[i]
      n_k[c_i] <- n_k[c_i] - 1L
      
      if (n_k[c_i] == 0L) {
        # Slot-recycling trick
        n_k[c_i] <- n_k[Nclust]
        c_alloc[c_alloc == Nclust] <- c_i
        n_k <- n_k[-Nclust]
        Nclust <- Nclust - 1L
      }
      c_alloc[i] <- -1L
      
      logp <- numeric(Nclust + 1L)
      for (cc in seq_len(Nclust)) {
        y_c <- y[c_alloc == cc]
        logp[cc] <- log(n_k[cc]) + log_marg(y[i], y_c, kappa)
      }
      logp[Nclust + 1L] <- log(alpha) + log_marg(y[i], numeric(0), kappa)
      
      logp  <- logp - max(logp)
      probs <- exp(logp); probs <- probs / sum(probs)
      newz  <- sample.int(Nclust + 1L, 1L, prob = probs)
      
      if (newz == Nclust + 1L) {
        n_k    <- c(n_k, 0L)
        Nclust <- Nclust + 1L
      }
      c_alloc[i]  <- newz
      n_k[newz]   <- n_k[newz] + 1L
    }
    
    #End : simulate (mu_c, s2_c) for occupied clusters ---
    mu_c <- numeric(Nclust)
    s2_c <- numeric(Nclust)
    for (cc in seq_len(Nclust)) {
      y_c      <- y[c_alloc == cc]
      nk_c     <- length(y_c)
      ybar     <- mean(y_c)
      Sk       <- sum((y_c - ybar)^2)
      lambda_n <- 1 / kappa + nk_c
      mu_n     <- (mu0 / kappa + nk_c * ybar) / lambda_n
      a_n      <- a_sigma + nk_c / 2
      b_n      <- b_sigma + 0.5 * Sk +
        nk_c * (ybar - mu0)^2 / (2 * (1 + nk_c * kappa))
      s2_c[cc] <- 1 / rgamma(1, shape = a_n, rate = b_n)
      mu_c[cc] <- rnorm(1, mu_n, sqrt(s2_c[cc] / lambda_n))
    }
    
    #Kappa Gibbs update (occupied only) 
    kappa <- update_kappa(mu_c, s2_c, mu0, a_kappa, b_kappa,
                          occupied = seq_len(Nclust))
    
    # Alpha update 
    eta    <- rbeta(1, alpha + 1, n)
    A      <- alpha_shape + Nclust - 1
    B      <- alpha_rate - log(eta)
    pi_eta <- A / (A + n * B)
    alpha  <- if (runif(1) < pi_eta)
      rgamma(1, alpha_shape + Nclust,     B)
    else
      rgamma(1, alpha_shape + Nclust - 1, B)
    
    # Draw a "new" cluster from G_0(kappa) for posterior predictive
    s2_new <- 1 / rgamma(1, a_sigma, b_sigma)
    mu_new <- rnorm(1, mu0, sqrt(kappa * s2_new))
    
    weights <- c(n_k / (n + alpha), alpha / (n + alpha))
    means   <- c(mu_c, mu_new)
    sds     <- c(sqrt(s2_c), sqrt(s2_new))
    
    draws[[it]] <- list(weights = weights, means = means, sds = sds,
                        alpha = alpha, kappa = kappa)
  }
  
  draws[(n_burn + 1):n_iter]
}
