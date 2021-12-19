data {
  // Constants
  int<lower=1> B; // Number of beta parameters
  int<lower=0> G; // Number of environmental covariates
  int<lower=6> N; // Number of brood years
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
  // Response index limit
  int<lower=5> I = N - 1;
  // Response variable: log recruits per spawner
  real y[I];
  for (i in 1:I) {
    y[i] = log((recruits[i] + fudge) / (spawners[i] + fudge));
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
  real y_hat[I] = rep_array(alpha, I);
  // Populate predicted response
  for (i in 1:I) {
    // Add intrinsic effects
    for (b in 1:((i < B) ? i : B)) {
      y_hat[i] += beta[b] * spawners[i - b + 1];
    }
    // Add extrinsic effects
    if (G > 0) {
      // Iterate over extrinsic
      for (g in 1:G) {
        y_hat[i] += gamma[g] * environs[i, g];
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

generated quantities {
  // Define
  int<lower=0> time = N;
  real<lower=0> observed = recruits[N];
  // Initialize
  real mean_y_forecast = alpha;
  real y_forecast;
  real<lower=0> forecast;
  // Compute mean y forecast
  for (b in 1:B) {
    mean_y_forecast += beta[b] * spawners[N - b + 1];
  }
  if (G > 0) {
    for (g in 1:G) {
      mean_y_forecast += gamma[g] * environs[N, g];
    }
  }
  // Generate y forecast
  y_forecast = normal_rng(mean_y_forecast, sigma);
  // Define forecast
  forecast = spawners[N] * exp(y_forecast);
}
