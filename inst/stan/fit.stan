data {
  // Constants
  int<lower=1> B; // Number of beta parameters
  int<lower=0> G; // Number of environmental covariates
  int<lower=5> N; // Number of brood years
  // Observations
  real environs[(G > 0) ? N : 0, G];
  real<lower=0> recruits[N]; // Recruitment abundance
  real<lower=0> spawners[N]; // Spawner abundance
  // Prior means
  real mu_alpha;
  real mu_beta[B];
  real mu_gamma[G];
  real mu_sigma;
  // Prior standard deviations
  real<lower=0> sd_alpha;
  real<lower=0> sd_beta[B];
  real<lower=0> sd_gamma[G];
  real<lower=0> sd_sigma;
  // Fudge factors
  real<lower=0> fudge;
}

transformed data {
  // Response variable: log recruits per spawner
  real y[N];
  for (n in 1:N) {
    y[n] = log((recruits[n] + fudge) / (spawners[n] + fudge));
  }
}

parameters {
  // Demographic rates
  real alpha; // Per-captica population growth rate at low abundance
  real beta[B]; // Density dependent mortality parameters
  real gamma[G]; // Environmental covariate parameter
  // Process error
  real<lower=0> sigma;
}

model {
  // Initialize predicted response
  real y_hat[N] = rep_array(alpha, N);
  // Populate predicted response
  for (n in 1:N) {
    // Add intrinsic effects
    for (b in 1:((n < B) ? n : B)) {
      y_hat[n] += beta[b] * spawners[n - b + 1];
    }
    // Add extrinsic effects
    if (G > 0) {
      // Iterate over extrinsic
      for (g in 1:G) {
        y_hat[n] += gamma[g] * environs[n, g];
      }
    }
  }
  // Priors
  alpha ~ normal(mu_alpha, sd_alpha);
  beta ~ normal(mu_beta, sd_beta);
  if (G > 0) {
    gamma ~ normal(mu_gamma, sd_gamma);
  }
  sigma ~ normal(mu_sigma, sd_sigma);
  // Sampling statement
  y ~ normal(y_hat, sigma);
}
