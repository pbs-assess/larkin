#' Simulate Sockeye Dynamics From A Larkin Stock-Recruitment Model
#'
#' @param a [numeric()] per-capita population growth rate at low abundance
#' @param b [numeric()] [vector()] density dependence parameters
#' @param init [numeric()] [vector()] initial spawner abundances
#' @param h [numeric()] constant harvest rate
#' @param p_bar [numeric()] [vector()] mean age at maturity proportions
#' @param omega [numeric()] age at maturity proportion deviation parameter
#' @param rho [numeric()] process error autocorrelation parameter
#' @param sigma [numeric()] process error standard deviation parameter
#' @param phi_1 [numeric()] initial process error deviation
#' @param burn [integer()] number of burn-in steps before the simulation
#' @param span [integer()] number of steps in the simulation
#' @param seed [integer()] passed internally to [set.seed()]
#'
#' @return [data.frame()] simulated spawner and recruitment at age abundances
#'
#' @export
#'
#' @examples
#' s1 <- sim()
#'
sim <- function(
  a = 7,
  b = c(6e-4, 2e-4, 1e-4, 1e-4),
  init = rep(1, 8),
  h = 0.2,
  p_bar = c(0.003, 0.917, 0.08),
  omega = 0.6,
  rho = 0.5,
  sigma = 0.8,
  phi_1 = 0.1,
  burn = 100,
  span = 1000,
  seed = NULL
) {

  # Set seed -------------------------------------------------------------------

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Initialize indexes ---------------------------------------------------------

  num_init <- 8L
  num_iter <- num_init + burn + span
  ind_init <- seq_len(num_init)
  ind_data <- (num_init + burn + 1):(num_iter)

  # Instantiate state variables ------------------------------------------------

  S <- rep(NA, num_iter)
  R <- matrix(NA, nrow = num_iter, ncol = 3)
  colnames(R) <- c("R_3", "R_4", "R_5")

  # Initialize state variables -------------------------------------------------

  S[ind_init] <- init

  # Instantiate simulation parameters ------------------------------------------

  p_t <- matrix(NA, nrow = num_iter, ncol = 3)
  phi_t <- rep(NA, num_iter)

  # Simulate state variables ---------------------------------------------------

  for (i in seq_len(num_iter)) {
    # Define age at maturity proportions
    p_numer <- exp(omega * stats::rnorm(3, 0, 1)) * p_bar
    p_denom <- sum(p_numer)
    p_t[i, ] <- p_numer / p_denom
    # Define autocorrelated process error (see Holt and Peterman 2008)
    if (i == 1) {
      phi_t[i] <- phi_1
    } else {
      phi_t[i] <- rho * phi_t[i - 1] + stats::rnorm(1, -sigma^2/2, sigma)
    }
    # Update state variables
    if (i >= 4) {
      R[i, ] <- p_t[i, ] * a * S[i] * exp(sum(-b * S[i:(i - 3)]) + phi_t[i])
    }
    if (i >= 8) {
      S[i + 1] <- (1 - h) * (R[i - 4, 3] + R[i - 3, 2] + R[i - 2, 1])
    }
    if (S[i + 1] < 2e-06) {
      S[i + 1] <- 0
    }
  }

  # Return simulated data ------------------------------------------------------

  N <- data.frame(
    Index = seq_along(ind_data),
    R_3 = R[ind_data, 1],
    R_4 = R[ind_data, 2],
    R_5 = R[ind_data, 3],
    S = S[ind_data]
  )
  return(N)
}
