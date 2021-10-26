data {
  // Constants
  int<lower=6> T; // Number of brood years
  // Observations
  real<lower=0> recruits[T]; // Recruitment abundance
  real<lower=0> spawners[T]; // Spawner abundance
  // Prior means
  real mu_alpha;
  real mu_beta[4];
  real mu_sigma;
  // Prior standard deviations
  real<lower=0> sd_alpha;
  real<lower=0> sd_beta;
  real<lower=0> sd_sigma;
  // Fudge factors
  real<lower=0> fudge;
}

transformed data {
  // Response index limit
  int<lower=5> N = T - 1;
  // Initialize response variable
  real y[N];
  // Response variable: log recruits per spawner
  for (n in 1:N) {
    y[n] = log((recruits[n] + fudge) / (spawners[n] + fudge));
  }
}

parameters {
  // Demographic rates
  real alpha; // Per-captica population growth rate at low abundance
  real<lower=0> beta[4]; // Density dependent mortality parameters
  // Process error
  real<lower=0> sigma;
}

model {
  // Predicted response
  real y_hat[N];
  y_hat[1] = alpha - beta[1] * spawners[1];
  y_hat[2] = alpha - beta[1] * spawners[2] - beta[2] * spawners[1];
  y_hat[3] = alpha - beta[1] * spawners[3] - beta[2] * spawners[2] - beta[3] * spawners[1];
  for (n in 4:N) {
    y_hat[n] = alpha
    - beta[1] * spawners[n]
    - beta[2] * spawners[n - 1]
    - beta[3] * spawners[n - 2]
    - beta[4] * spawners[n - 3];
  }
  // Priors
  alpha ~ normal(mu_alpha, sd_alpha);
  beta ~ normal(mu_beta, sd_beta);
  sigma ~ normal(mu_sigma, sd_sigma);
  // Sampling statement
  y ~ normal(y_hat, sigma);
}

generated quantities {
  // Define
  int<lower=0> time = T;
  real<lower=0> observed = recruits[T];
  real y_forecast = normal_rng(
    alpha
    - beta[1] * spawners[T]
    - beta[2] * spawners[T - 1]
    - beta[3] * spawners[T - 2]
    - beta[4] * spawners[T - 3],
    sigma
  );
  real<lower=0> forecast = spawners[T] * exp(y_forecast);
}
