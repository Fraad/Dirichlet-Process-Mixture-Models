//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=1> K;
  int<lower=1> N;
  int<lower=1> N_grid;
  int<lower=1> N_anchor;
  array[N] real y;
  vector[N_grid] grid_pts;
  vector[N_anchor] anchor_pts;
  real<lower=0> a_sigma;
  real<lower=0> b_sigma;
  real mu_0;
  real<lower=0> a_kappa;
  real<lower=0> b_kappa;
  real<lower=0> alpha_a;
  real<lower=0> alpha_b;
}
parameters {
  vector[K] z_mu;
  vector<lower=0>[K] sigma_sq;
  vector<lower=0, upper=1>[K-1] v;
  real<lower=0> alpha;
  real<lower=0> kappa;
}
transformed parameters {
  simplex[K] pi_var;
  vector[K] mu;
  {
    real rem = 1.0;
    for (k in 1:(K-1)) {
      pi_var[k] = v[k] * rem;
      rem -= pi_var[k];
    }
    pi_var[K] = rem;
  }
  for (k in 1:K)
    mu[k] = mu_0 + sqrt(kappa * sigma_sq[k]) * z_mu[k];
}
model {
  kappa    ~ inv_gamma(a_kappa, b_kappa);
  alpha    ~ gamma(alpha_a, alpha_b);
  v        ~ beta(1, alpha);
  sigma_sq ~ inv_gamma(a_sigma, b_sigma);
  z_mu     ~ std_normal();
  vector[K] lpi = log(pi_var);
  for (n in 1:N) {
    vector[K] lps = lpi;
    for (k in 1:K)
      lps[k] += normal_lpdf(y[n] | mu[k], sqrt(sigma_sq[k]));
    target += log_sum_exp(lps);
  }
}
generated quantities {
  array[N] int<lower=1, upper=K> z_hard;
  int<lower=1, upper=K> K_occ;
  vector[N_grid]   f_pred;
  vector[N_anchor] log_f_anchor;

  // Hard cluster assignments via categorical_logit on posterior responsibilities
  for (n in 1:N) {
    vector[K] lp = log(pi_var);
    for (k in 1:K)
      lp[k] += normal_lpdf(y[n] | mu[k], sqrt(sigma_sq[k]));
    z_hard[n] = categorical_logit_rng(lp);
  }
  {
    array[K] int occ_flag = rep_array(0, K);
    for (n in 1:N) occ_flag[z_hard[n]] = 1;
    K_occ = sum(occ_flag);
  }
  
  for (g in 1:N_grid) {
    vector[K] lp = log(pi_var);
    for (k in 1:K)
      lp[k] += normal_lpdf(grid_pts[g] | mu[k], sqrt(sigma_sq[k]));
    f_pred[g] = exp(log_sum_exp(lp));
  }
  
  for (a in 1:N_anchor) {
    vector[K] lp = log(pi_var);
    for (k in 1:K)
      lp[k] += normal_lpdf(anchor_pts[a] | mu[k], sqrt(sigma_sq[k]));
    log_f_anchor[a] = log_sum_exp(lp);
  }
  
}