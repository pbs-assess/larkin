data {
  // Constants
  int<lower=6> N; // Number of brood years
  int<lower=1> B; // Number of beta parameters
  int<lower=0> G; // Number of environmental covariates
  // Observations
  vector[N]<lower=0> recruits; // Recruitment abundance
  vector[N]<lower=0> spawners; // Spawning abundance
  array[G ? N : 0, G] real environs; // Environmental indices
  // Prior means
  real prior_mean_alpha;
  vector[B] prior_mean_beta;
  vector[G] prior_mean_gamma;
  real prior_mean_sigma;
  // Prior standard deviations
  real<lower=0> prior_sd_alpha;
  vector[B]<lower=0> prior_sd_beta;
  vector[G]<lower=0> prior_sd_gamma;
  real<lower=0> prior_sd_sigma;
  real<lower=0> prior_sd_phi;
  // Fudge factors
  real<lower=0> fudge;
}

transformed data {
  // Response index limit
  int<lower=5> I = N - 1;
  // Response: log recruits per spawner
  vector[I] y = log(recruits[1:I] + fudge) - log(spawners[1:I] + fudge);
  // Timevary indicator
  int<lower=0, upper=1> timevary;
  if (prior_sd_phi == 0) {
    timevary = 0;
  } else {
    timevary = 1;
  }
}

parameters {
  // Demographic rates
  vector[timevary ? I : 1] aleph; // Components of alpha
  vector[B]<upper=0> beta; // Density dependent mortality parameters
  vector[G] gamma; // Environmental covariate parameters
  // Standard deviation parameters
  real<lower=0> sigma; // Process error
  real<lower=0> phi; // Random walk on alpha
}

transformed parameters {
  // Per-capita productivity at low abundance
  vector[I] alpha;
  if (timevary) {
    alpha = aleph;
  } else {
    alpha = rep_vector(aleph[1], I);
  }
}

model {
  // Initialize expected response
  vector[I] y_hat = alpha;
  // Populate expected response
  for (i in 1:I) {
    // Add intrinsic effects
    for (b in 1:((i < B) ? i : B)) {
      y_hat[i] += beta[b] * spawners[i - b + 1];
    }
    // Add extrinsic effects
    if (G > 0) {
      for (g in 1:G) {
        y_hat[i] += gamma[g] * environs[i, g];
      }
    }
  }
  // Random walk
  // TODO


  // Priors
  alpha ~ normal(mu_alpha, sd_alpha);
  // TODO


  beta ~ normal(prior_mean_beta, prior_sd_beta);
  if (G > 0) {
    gamma ~ normal(prior_mean_gamma, prior_sd_gamma);
  }
  sigma ~ normal(prior_mean_sigma, prior_sd_sigma);
  if (timevary) {
    phi ~ prior_sd_phi * std_normal();
  }
  // Sampling statement
  y ~ normal(y_hat, sigma);
}

generated quantities {
  // TODO wrap in curly braces so not reporting parts don't need?


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
