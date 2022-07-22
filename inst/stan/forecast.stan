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
  real<lower=0> prior_mean_sigma;
  real<lower=0> prior_mean_omega;
  // Prior standard deviations
  real<lower=0> prior_sd_alpha;
  vector[B]<lower=0> prior_sd_beta;
  vector[G]<lower=0> prior_sd_gamma;
  real<lower=0> prior_sd_sigma;
  real<lower=0> prior_sd_omega;
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
  if (prior_mean_omega || prior_sd_omega) {
    timevary = 1;
  } else {
    timevary = 0;
  }
}

parameters {
  // Demographic rates
  vector[timevary ? I : 1] aleph; // Components of alpha
  vector[B]<upper=0> beta; // Density dependent mortality parameters
  vector[G] gamma; // Environmental covariate parameters
  // Standard deviation parameters
  real<lower=0> sigma; // Process error
  vector[timevary ? 1 : 0]<lower=0> omega; // Random walk on alpha
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
  // Priors
  alpha[1] ~ normal(prior_mean_alpha, prior_sd_alpha);
  beta ~ normal(prior_mean_beta, prior_sd_beta);
  if (G > 0) {
    gamma ~ normal(prior_mean_gamma, prior_sd_gamma);
  }
  sigma ~ normal(prior_mean_sigma, prior_sd_sigma);
  // Non-centred random walk
  if (timevary) {
    // Random walk
    (alpha[2:I] - alpha[1:(I-1)]) ./ omega[1] ~ std_normal();
    // Random walk standard deviation prior
    omega[1] ~ normal(prior_mean_omega, prior_sd_omega);
  }
  // Sampling statement
  y ~ normal(y_hat, sigma);
}

generated quantities {
  // Declare output
  int<lower=0> time = N;
  real<lower=0> observed = recruits[N];
  real<lower=0> forecast;
  // Declare not output
  {
    // Components of forecast
    real mean_y_forecast;
    real y_forecast
  }
  // Initialize
  if (timevary) {
    mean_y_forecast = normal_rng(alpha[I], omega);
  } else {
    mean_y_forecast = alpha[1];
  }
  // Add intrinsic effects
  for (b in 1:B) {
    mean_y_forecast += beta[b] * spawners[N - b + 1];
  }
  // Add extrinsic effects
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
