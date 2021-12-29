#' Fit A Larkin or Ricker Model Via CmdStanR
#'
#' @param data [list()]
#' @param chains [integer()] number of chains
#' @param step_size [integer()] initial step size
#' @param iter_warmup [integer()] number of warmup iterations
#' @param iter_sampling [integer()] number of sampling iterations
#' @param ... additional arguments to pass to \code{$sample()} method
#'
#' @return [list()]
#' @export
#'
fit <- function (data,
                 chains = 1,
                 step_size = 0.01,
                 iter_warmup = 250,
                 iter_sampling = 750,
                 ...) {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_list(data, c("double", "integer", "array", "matrix"))

  # Create model object --------------------------------------------------------

  mod <- cmdstanr::cmdstan_model(
    system.file("stan", paste0("fit", ".stan"), package = "larkin"),
    include_path = system.file("stan", package = "larkin")
  )

  # Fit the model --------------------------------------------------------------

  fit <- mod$sample(
    data = data,
    chains = chains,
    step_size = step_size,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    ...
  )

  # Assemble summaries ---------------------------------------------------------

  # summaries <- summarise_posterior_draws(fit, data)

  # Return values --------------------------------------------------------------

  structure(list(
    data = data,
    draws = fit$draws(format = "df"),
    fit = fit,
    summary = fit$summary()),
    class = "fit")
}

#' Forecast From A Ricker Or Larkin Model Using Previous Values Only
#'
#' @param data [list()]
#' @param index [integer()]
#' @param buffer [integer()] number of forecasts prior to \code{index}
#' @param metric [character()]
#' @param beyond [logical()] TBD
#' @param cores [numeric()] number of cores for parallel processing
#' @param chains [integer()] number of chains
#' @param step_size [integer()] initial step size
#' @param iter_warmup [integer()] number of warmup iterations
#' @param iter_sampling [integer()] number of sampling iterations
#' @param ... additional arguments to pass to \code{$sample()} method
#'
#' @importFrom rlang .data
#' @importFrom rlang :=
#'
#' @return [tibble::tibble()]
#' @export
#'
forecast <- function (data,
                      index,
                      buffer = 10,
                      metric = "mamse",
                      beyond = FALSE,
                      cores = NULL,
                      chains = 1,
                      step_size = 0.01,
                      iter_warmup = 250,
                      iter_sampling = 750,
                      ...) {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_list(data, c("double", "integer"))
  checkmate::assert_integerish(data$N, lower = 6, len = 1, any.missing = FALSE)
  checkmate::assert_integerish(index, lower = 6, len = 1, any.missing = FALSE)
  # Temporary
  checkmate::assert_integerish(index, upper = data$N)

  # Define indexes -------------------------------------------------------------

  indexes <- c(index:data$N)

  # Apply forecasting ----------------------------------------------------------

  if (is.null(cores)) {
    # Apply in sequence
    forecasts <- lapply(
      unique(c(indexes - buffer, indexes)),
      FUN = forecast_single_value,
      data = data,
      chains = chains,
      step_size = step_size,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      ...
    )
  } else {
    if (.Platform$OS.type == "unix") {
      forecasts <- parallel::mclapply(
        unique(c(indexes - buffer, indexes)),
        FUN = forecast_single_value,
        data = data,
        chains = chains,
        step_size = step_size,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        # ...,
        mc.cores = cores
      )
    } else {
      stop("parallel forecast via sockets not yet implemented R/utils.R")
    }
  }

  # Return forecasts -----------------------------------------------------------

  # Bind rows
  forecasts <- forecasts %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(.data$time)
  # Define observed-forecast matrix
  m <- matrix(c(forecasts$observed, forecasts$forecast), ncol = 2L)
  # Return forecasts
  forecasts %>%
    dplyr::mutate(
      mre = runner::runner(m, f = larkin::matric, fun = larkin::mre),
      !!metric := runner::runner(m, f = larkin::matric, fun = get(metric))
    ) %>%
    dplyr::filter(.data$time %in% indexes)
}

#' Simulate Sockeye Dynamics From A Larkin Stock-Recruitment Model
#'
#' @param alpha [numeric()] per-capita log population growth rate at low
#'   spawner abundance
#' @param beta [numeric()] [vector()] density dependence parameters
#' @param init [numeric()] [vector()] initial spawner abundances
#' @param p_bar [numeric()] [vector()] mean age at maturity proportions
#' @param phi [numeric()] process error autocorrelation parameter
#' @param sigma [numeric()] process error standard deviation parameter
#' @param burn [integer()] number of burn-in steps before the simulation
#' @param span [integer()] number of steps in the simulation
#' @param harvest [numeric()] harvest rate
#' @param gamma [numeric()] environmental variable influence parameter
#' @param environs [tibble::tibble()] environmental variables
#' @param extirp [numeric()] extirpation threshold
#'
#' @return [data.frame()] simulated spawner and recruitment at age abundances
#'
#' @export
#'
#' @examples
#' s1 <- sim()
#'
#' environs <- simulate_environmental_covariates(100, rep(0, 3), 0.5, 1)
#' s2 <- sim(gamma = rep(0.5, 3), environs = environs)
#'
sim <- function (alpha = 2,
                 beta = c(-8, -6, -4, -2),
                 init = rep(0.1, 8),
                 p_bar = c(0.003, 0.917, 0.08),
                 phi = 0,
                 sigma = 0.1,
                 burn = 100,
                 span = 100,
                 harvest = 0.2,
                 gamma = NULL,
                 environs = NULL,
                 extirp = 1e-6) {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_number(alpha, finite = TRUE)
  checkmate::assert_numeric(beta, finite = TRUE, any.missing = FALSE, len = 4)
  checkmate::assert_numeric(init, finite = TRUE, any.missing = FALSE, len = 8)
  checkmate::assert_numeric(p_bar, finite = TRUE, any.missing = FALSE, len = 3)
  checkmate::assert_number(phi, finite = TRUE)
  checkmate::assert_number(sigma, lower = 0, finite = TRUE)
  checkmate::assert_integerish(burn, lower = 0, any.missing = FALSE, len = 1)
  checkmate::assert_integerish(span, lower = 0, any.missing = FALSE, len = 1)
  checkmate::assert_number(harvest, lower = 0, finite = TRUE)
  checkmate::assert_numeric(gamma, finite = TRUE, null.ok = TRUE)
  checkmate::assert_data_frame(environs, null.ok = TRUE)
  checkmate::assert_number(extirp)

  # Check more arguments -------------------------------------------------------

  if (!is.null(environs)) {
    checkmate::assert_true(nrow(environs) == span)
  }
  if (!is.null(gamma) & ! is.null(environs)) {
    checkmate::assert_true(length(gamma) == ncol(environs))
  }

  # Define index limits --------------------------------------------------------

  num_init <- 8L
  num_iter <- num_init + burn + span

  # Define indexes -------------------------------------------------------------

  ind_init <- seq_len(num_init)
  ind_span <- (num_init + burn + 1):(num_iter)

  # Buffer environs ------------------------------------------------------------

  if (!is.null(environs)) {
    means <- apply(environs, 2, mean)
    sds <- apply(environs, 2, stats::sd)
    buffer <- matrix(0, nrow = num_init + burn, ncol = ncol(environs))
    colnames(buffer) <- colnames(environs)
    for (j in seq_len(ncol(environs))) {
      buffer[, j] <- stats::rnorm(num_init + burn, means[j], sds[j])
    }
    buffer <- tibble::as_tibble(buffer)
    environs <- dplyr::bind_rows(buffer, environs)
  }

  # Initialize state variables -------------------------------------------------

  returns <- rep(NA_real_, num_iter)
  spawners <- c(init, rep(NA_real_, num_iter - num_init))
  r_3 <- rep(NA_real_, num_iter)
  r_4 <- rep(NA_real_, num_iter)
  r_5 <- rep(NA_real_, num_iter)
  recruits <- rep(NA_real_, num_iter)

  # Initialize proportions -----------------------------------------------------

  h_t <- rep(NA_real_, num_iter)
  p_t <- matrix(NA_real_, nrow = num_iter, ncol = 3L)

  # Initialize vectors ---------------------------------------------------------

  growth <- rep(NA_real_, num_iter)
  ddmort <- rep(NA_real_, num_iter)
  environ <- rep(NA_real_, num_iter)
  epsilon <- rep(NA_real_, num_iter)

  # Simulate -------------------------------------------------------------------

  for (i in seq_len(num_iter)) {
    # Define proportions
    p_t[i, ] <- p_bar
    # Define process error
    if (i == 1) {
      epsilon[i] <- stats::rnorm(1, 0, sigma)
    } else {
      epsilon[i] <- phi * epsilon[i - 1] + stats::rnorm(1, 0, sigma)
    }
    # Define environmental influence
    if (is.null(gamma) | is.null(environs)) {
      environ[i] <- 0
    } else {
      environ[i] <- sum(gamma * environs[i, ])
    }
    # Recruitment
    if (i >= 4) {
      # Realized growth rate
      ddmort[i] <- sum(beta * spawners[i:(i - 3)])
      growth[i] <- exp(alpha + ddmort[i] + environ[i] + epsilon[i])
      # Recruitment
      r_3[i] <- p_t[i, 1] * spawners[i] * growth[i]
      r_4[i] <- p_t[i, 2] * spawners[i] * growth[i]
      r_5[i] <- p_t[i, 3] * spawners[i] * growth[i]
      recruits[i] <- r_3[i] + r_4[i] + r_5[i]
    }
    # Spawners
    if (i >= 8 & i < num_iter) {
      # Harvest
      h_t[i + 1] <- harvest
      # Returns
      returns[i + 1] <- sum(r_5[i - 4], r_4[i - 3], r_3[i - 2])
      # Spawners
      spawners[i + 1] <- (1 - h_t[i + 1]) * returns[i + 1]
      # Extirpation
      if (spawners[i + 1] < extirp) {
        spawners[i + 1] <- 0
      }
    }
  }

  # Return simulated data ------------------------------------------------------

  tibble::tibble(
    time = seq_along(ind_span),
    if (is.null(environs)) NULL else environs[ind_span, ],
    returns = returns[ind_span],
    spawners = spawners[ind_span],
    r_3 = r_3[ind_span],
    r_4 = r_4[ind_span],
    r_5 = r_5[ind_span],
    recruits = recruits[ind_span]
  )
}

#' Simulate Environmental Covariates
#'
#' @param n [integer()] number of time steps
#' @param mu [numeric()][vector()] means
#' @param phi [numeric()][vector()]
#' @param sigma [numeric()][vector()]
#' @param var_names [character()]
#'
#' @return [tibble::tibble()]
#' @export
#'
#' @examples
#'
#' simulate_environmental_covariates(15, rep(0, 5), 0.5, 1)
#'
simulate_environmental_covariates <- function (n,
                                               mu = 0,
                                               phi = 0,
                                               sigma = 1,
                                               var_names = NULL) {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_integerish(n, lower = 1, any.missing = FALSE, len = 1)
  checkmate::assert_double(mu, finite = T, any.missing = FALSE, min.len = 1)
  checkmate::assert_double(phi, finite = T, any.missing = FALSE, min.len = 1)
  checkmate::assert_double(sigma, finite = T, any.missing = FALSE, min.len = 1)
  checkmate::assert_character(var_names, any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_true(length(phi) %in% c(1, length(mu)))
  checkmate::assert_true(length(sigma) %in% c(1, length(mu)))
  checkmate::assert_true(length(var_names) %in% c(0, length(mu)))

  # Convert argument lengths ---------------------------------------------------

  if (length(mu) > 1) {
    if (length(phi) == 1) {
      phi <- rep(phi, length(mu))
    }
    if (length(sigma) == 1) {
      sigma <- rep(sigma, length(mu))
    }
  }

  # Simulate -------------------------------------------------------------------

  m <- matrix(NA, nrow = n, ncol = length(mu))
  for (i in seq_len(n)) {
    for (j in seq_along(mu)) {
      if (i == 1) {
        m[i, j] <- stats::rnorm(n = 1, mean = mu, sd = sigma)
      } else {
        m[i, j] <- phi[j] * m[i - 1, j] + stats::rnorm(n = 1, mean = mu, sd = sigma)
      }
    }
  }
  if (is.null(var_names)) {
    colnames(m) <- paste0("covar_", seq_along(mu))
  } else {
    colnames(m) <- var_names
  }

  # Return ---------------------------------------------------------------------

  return(tibble::as_tibble(m, .name_repair = "universal"))
}
