
functions {

}

data {
  // Constants
  int<lower=0> T; // Number of brood years
  // Observations
  real recruits[T]; // Array of recruit abundance
  real spawners[T]; // Array of spawner abundance
  // Prior parameters
  real mu_alpha;
  real sd_alpha;
  real mu_beta[4];
  real sd_beta[4];
}

parameters {
  real alpha; // Per-capita population growth rate at low abundance
  real<lower=0> beta[4]; // Density dependence parameters
  real<lower=0> sigma; // Process error standard deviation
}

model {
  real log_recruits_hat;
  for (t in 4:T) {
    // Compute log predicted recruits
    log_recruits_hat = log(spawners[t])
    + alpha
    - beta[1] * spawners[t]
    - beta[2] * spawners[t - 1]
    - beta[3] * spawners[t - 2]
    - beta[4] * spawners[t - 3];
    // Define sampling statement
    recruits[t] ~ lognormal(log_recruits_hat, sigma);
  }
  // Define priors
  alpha ~ normal(mu_alpha, sd_alpha);
  beta ~ normal(mu_beta, sd_beta);
}

generated quantities {

}
