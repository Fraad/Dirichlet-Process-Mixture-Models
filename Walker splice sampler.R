

# Implementation of Walker splice sampler

# Using a standard conjugate Normal-InvGamma conjugate base measure with Normal likelihood. 

# Walker (2007) slice sampler — DP mixture of Normals
# Base measure (NIG):  sigma^2_k ~ IG(sigma_a, sigma_b)
#                      mu_k | sigma^2_k ~ N(mu0, kappa * sigma^2_k)
# Concentration:       alpha ~ Gamma(a_alpha, b_alpha)  

update_kappa <- function(mu, s2, mu0, a_kappa, b_kappa, occupied) {
  mu_occ  <- mu[occupied]
  s2_occ  <- s2[occupied]
  K_occ   <- length(mu_occ)
  S_tau   <- sum((mu_occ - mu0)^2 / s2_occ)
  shape_n <- a_kappa + K_occ / 2
  rate_n  <- b_kappa + 0.5 * S_tau
  1 / rgamma(1, shape = shape_n, rate = rate_n)
}


walker_slice_nig <- function(y, n_iter = 2000, n_burn, 
                             mu0 = 0,
                             sigma_a = 2, sigma_b = 1,
                             a_kappa = 2, b_kappa = 1,
                             alpha_shape = alpha_shape , alpha_rate = alpha_rate) {
  n <- length(y)
  
  # draw (mu, sigma^2) from the NIG prior given current kappa
  rprior_NIG <- function(kappa_cur) {
    s2 <- 1 / rgamma(1, shape = sigma_a, rate = sigma_b)
    list(mu = rnorm(1, mu0, sqrt(kappa_cur * s2)), s2 = s2)
  }
  
  # --- initialise ---
  kappa <- 1 / rgamma(1, a_kappa, b_kappa)
  alpha <- rgamma(1, alpha_shape, alpha_rate)
  d     <- rep(1L, n)
  v     <- rbeta(1, 1, alpha)
  w     <- v
  init  <- rprior_NIG(kappa)
  mu    <- init$mu;  s2 <- init$s2
  K     <- 1L
  
  draws <- vector("list", n_iter)
  
  for (it in seq_len(n_iter)) {
    
    ## 1. slice variables
    u     <- runif(n, 0, w[d])
    u_min <- min(u)
    
    ## 2. extend sticks until remaining mass < u_min
    while (1 - sum(w) > u_min) {
      v_new <- rbeta(1, 1, alpha)
      w     <- c(w, v_new * (1 - sum(w)))
      v     <- c(v, v_new)
      new   <- rprior_NIG(kappa)
      mu    <- c(mu, new$mu);  s2 <- c(s2, new$s2)
      K     <- K + 1L
    }
    
    ## 3. reallocate
    for (i in seq_len(n)) {
      elig <- which(w > u[i])
      lp   <- dnorm(y[i], mu[elig], sqrt(s2[elig]), log = TRUE)
      p    <- exp(lp - max(lp))
      d[i] <- if (length(elig) == 1L) elig else sample(elig, 1, prob = p)
    }
    
    ## 3b. Trim unoccupied tail sticks: only retain j = 1, ..., k_star
    k_star_trim <- max(d)
    if (K > k_star_trim) {
      v  <- v[seq_len(k_star_trim)]
      w  <- w[seq_len(k_star_trim)]
      mu <- mu[seq_len(k_star_trim)]
      s2 <- s2[seq_len(k_star_trim)]
      K  <- k_star_trim
    }
    
    
    ## 4. Escobar–West update for alpha | d
    K_star <- length(unique(d))
    eta    <- rbeta(1, alpha + 1, n)
    A      <- alpha_shape + K_star - 1
    B      <- alpha_rate - log(eta)
    pi_eta <- A / (A + n * B)
    alpha  <- if (runif(1) < pi_eta)
      rgamma(1, alpha_shape + K_star,     B)
    else
      rgamma(1, alpha_shape + K_star - 1, B)
    
    ## 5. Walker step C: update v_j from truncated Beta(1, alpha) on (L_j, U_j)
    k_star  <- max(d)
    log_1mv <- log1p(-v)
    cs      <- cumsum(log_1mv)
    log_u   <- log(u)
    
    for (j in seq_len(k_star)) {
      
      log_prod_lt_j <- if (j == 1L) 0 else cs[j - 1L]
      idx_eq <- which(d == j)
      L_j <- if (length(idx_eq))
        exp(max(log_u[idx_eq]) - log_prod_lt_j) else 0
      
      idx_gt <- which(d > j)
      if (length(idx_gt)) {
        ki  <- d[idx_gt]
        log_prod_lt_ki <- cs[ki - 1L]
        log_ratio <- log_u[idx_gt] - log(v[ki]) - (log_prod_lt_ki - log_1mv[j])
        U_j <- 1 - exp(max(log_ratio))
      } else {
        U_j <- 1
      }
      
      L_j <- max(0, min(1, L_j))
      U_j <- max(0, min(1, U_j))
      if (U_j <= L_j) next
      
      A_cdf <- (1 - L_j)^alpha
      B_cdf <- (1 - U_j)^alpha
      xi    <- runif(1)
      v[j]  <- 1 - (A_cdf - xi * (A_cdf - B_cdf))^(1 / alpha)
      
      log_1mv[j] <- log1p(-v[j])
      cs <- cumsum(log_1mv)
    }
    
    ## Rebuild stick-breaking weights
    w <- v * c(1, cumprod(1 - v[-K]))
    
    ## 6. kappa | occupied (mu_k, sigma_k^2)  -- IG conjugate
    occupied <- sort(unique(d))
    K_occ    <- length(occupied)
    ss_kappa <- sum((mu[occupied] - mu0)^2 / s2[occupied])
    kappa    <- 1 / rgamma(1,
                           shape = a_kappa + K_occ / 2,
                           rate  = b_kappa + 0.5 * ss_kappa)
    
    
    
    ## 7. (mu_k, sigma_k^2) | d, y, kappa  -- NIG conjugate (uses updated kappa)
    for (k in seq_len(K)) {
      idx  <- which(d == k)
      nk_k <- length(idx)
      if (nk_k > 0L) {
        ybar     <- mean(y[idx])
        Sk       <- sum((y[idx] - ybar)^2)
        lambda_n <- 1/kappa + nk_k
        mu_n     <- (mu0/kappa + nk_k * ybar) / lambda_n
        alpha_n  <- sigma_a + nk_k / 2
        beta_n   <- sigma_b + 0.5 * Sk +
          nk_k * (ybar - mu0)^2 / (2 * (1 + nk_k * kappa))
        s2[k]    <- 1 / rgamma(1, shape = alpha_n, rate = beta_n)
        mu[k]    <- rnorm(1, mu_n, sqrt(s2[k] / lambda_n))
      } else {
        new   <- rprior_NIG(kappa)
        mu[k] <- new$mu;  s2[k] <- new$s2
      }
    }
    
    draws[[it]] <- list(d = d, w = w, mu = mu, s2 = s2,
                        K_active = K_star, alpha = alpha, kappa = kappa)
  }
  #draws
  draws[(n_burn + 1):n_iter] 
}

