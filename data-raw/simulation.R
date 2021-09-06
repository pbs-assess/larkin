# Simulate population dynamics from a Larkin model
s1 <- sim(
  a = 7, # Per-capita population growth rate at low spawner abundance
  b = c(6, 2, 1, 1), # Density dependence parameters for spawners and lags
  init = rep(1, 8), # Initial spawner abundances
  h = 0.2, # Constant harvest rate
  p_bar = c(0.003, 0.917, 0.08), # Age at maturity mean proportions
  # omega = 0.6, # Age at maturity proportion deviation parameter
  # rho = 0.5, # Process error autocorrelation parameter
  sigma = 0.8, # Process error standard deviation parameter
  # phi_1 = 0.1, # Initial process error deviation
  burn = 100, # Number of burn-in steps before the simulation
  span = 1000, # Number of steps in the simulation
  seed = 109 # Seed value for set.seed()
)
simulation <- s1

# Write to data/
usethis::use_data(simulation, overwrite = TRUE)
