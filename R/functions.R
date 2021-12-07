#' Fit A Larkin or Ricker Model Via CmdStanR
#'
#' @param data [list()]
#' @param model [character()] either \code{"larkin"} or \code{"ricker"}
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
                 model = "larkin",
                 chains = 1,
                 step_size = 0.01,
                 iter_warmup = 250,
                 iter_sampling = 750,
                 ...) {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_list(data, c("double", "integer"))
  checkmate::assert_choice(model, c("larkin", "ricker"))

  # Create model object --------------------------------------------------------

  mod <- cmdstanr::cmdstan_model(
    system.file("stan", paste0(model, ".stan"), package = "larkin"),
    include_path = system.file("stan", package = "larkin")
  )

  # Fit the model --------------------------------------------------------------

  cmdfit <- mod$sample(
    data = data,
    chains = chains,
    step_size = step_size,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    ...
  )

  # Assemble parameters --------------------------------------------------------

  parameters <- NULL

  # Return values --------------------------------------------------------------

  structure(list(
    data = data,
    output = cmdfit$output_files(),
    summary = cmdfit$summary(),
    parameters = parameters),
    class = model)
}

#' Forecast From A Ricker Or Larkin Model Using Previous Values Only
#'
#' @param inds [numeric()] [vector()] indexes of values to forecast
#' @param data [list()]
#' @param model [character()] either \code{"larkin"} or \code{"ricker"}
#' @param window [integer()] length of the buffer between training and forecast
#' @param cores [numeric()] number of cores for parallel processing
#' @param chains [integer()] number of chains
#' @param step_size [integer()] initial step size
#' @param iter_warmup [integer()] number of warmup iterations
#' @param iter_sampling [integer()] number of sampling iterations
#' @param ... additional arguments to pass to \code{$sample()} method
#'
#' @importFrom rlang .data
#'
#' @return [tibble::tibble()]
#' @export
#'
forecast <- function (inds,
                      data,
                      model = "larkin",
                      window = 10,
                      cores = NULL,
                      chains = 1,
                      step_size = 0.01,
                      iter_warmup = 250,
                      iter_sampling = 750,
                      ...) {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_list(data, c("double", "integer"))
  checkmate::assert_choice(model, c("larkin", "ricker"))

  # Define indexes -------------------------------------------------------------

  inds_forecast <- unique(c(inds - window, inds))
  inds_return <- inds

  # Apply forecasting ----------------------------------------------------------

  if (is.null(cores)) {
    # Apply in sequence
    output <- lapply(
      inds_forecast,
      FUN = forecast_single_value,
      data = data,
      model = model,
      chains = chains,
      step_size = step_size,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      ...
    )
  } else {
    if (.Platform$OS.type == "unix") {
      output <- parallel::mclapply(
        inds_forecast,
        FUN = forecast_single_value,
        data = data,
        model = model,
        chains = chains,
        step_size = step_size,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        ...,
        mc.cores = cores
      )
    } else {
      stop("parallel forecast via sockets not yet implemented R/utils.R")
    }
  }

  # Return forecasts -----------------------------------------------------------

  dplyr::bind_rows(output) %>%
    dplyr::arrange(.data$time) %>%
    dplyr::mutate(
      # Metrics
      re = re(.data$observed, .data$forecast),
      pe = pe(.data$observed, .data$forecast),
      ape = ape(.data$observed, .data$forecast),
      aape = aape(.data$observed, .data$forecast),
      # Running metrics
      cmre = runner::runner(.data$re, f = mean_or_na),
      cmae = runner::runner(abs(.data$re), f = mean_or_na),
      cmpe = runner::runner(.data$pe, f = mean_or_na),
      cmape = runner::runner(.data$ape, f = mean_or_na),
      cmaape = runner::runner(.data$aape, f = mean_or_na),
      crmse = runner::runner(.data$re, f = rms_or_na),
      crho = runner::runner(
        matrix(c(.data$observed, .data$forecast), ncol = 2),
        f = rho_or_na
      )
    ) %>%
    dplyr::filter(
      .data$time %in% inds_return
    )
}

#' Simulate Sockeye Dynamics From A Larkin Stock-Recruitment Model
#'
#' @param alpha [numeric()] per-capita log population growth rate at low
#'   spawner abundance
#' @param beta [numeric()] [vector()] density dependence parameters
#' @param init [numeric()] [vector()] initial spawner abundances
#' @param p_bar [numeric()] [vector()] mean age at maturity proportions
#' @param sigma [numeric()] process error standard deviation parameter
#' @param burn [integer()] number of burn-in steps before the simulation
#' @param span [integer()] number of steps in the simulation
#' @param harvest [numeric()] harvest rate
#' @param gamma [numeric()] environmental variable influence parameter
#' @param x_phi [numeric()] environmental variable autocorrelation parameter
#' @param extirp [numeric()] extirpation threshold
#' @param seed [integer()] passed internally to [set.seed()]
#'
#' @return [data.frame()] simulated spawner and recruitment at age abundances
#'
#' @export
#'
#' @examples
#' s1 <- sim()
#'
sim <- function (alpha = 2,
                 beta = c(8, 6, 4, 2),
                 init = rep(0.1, 8),
                 p_bar = c(0.003, 0.917, 0.08),
                 sigma = 0.1,
                 burn = 100,
                 span = 100,
                 harvest = 0.2,
                 gamma = 0,
                 x_phi = 0,
                 extirp = 1e-6,
                 seed = NULL) {

  # Check arguments ------------------------------------------------------------

  # Set seed -------------------------------------------------------------------

  set.seed(seed)

  # Define index limits --------------------------------------------------------

  num_init <- 8L
  num_iter <- num_init + burn + span

  # Define indexes -------------------------------------------------------------

  ind_init <- seq_len(num_init)
  ind_data <- (num_init + burn + 1):(num_iter)

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

  # Initialize environmental variable ------------------------------------------

  x <- rep(NA_real_, num_iter)

  # Initialize process error ---------------------------------------------------

  epsilon <- rep(NA_real_, num_iter)

  # Simulate -------------------------------------------------------------------

  for (i in seq_len(num_iter)) {
    # Define proportions
    p_t[i, ] <- p_bar
    # Define environmental influence
    if (i == 1) {
      x[i] <- stats::rnorm(1, 0, 1)
    } else {
      x[i] <- x_phi * x[i - 1] + stats::rnorm(1, 0, 1)
    }
    # Define process error
    epsilon <- stats::rnorm(1, 0, sigma)
    # Recruitment
    if (i >= 4) {
      # Realized growth rate
      ddmort <- sum(-beta * spawners[i:(i - 3)])
      envir <- gamma * x[i]
      growth <- exp(alpha + ddmort + envir + epsilon)
      # Recruitment
      r_3[i] <- p_t[i, 1] * spawners[i] * growth
      r_4[i] <- p_t[i, 2] * spawners[i] * growth
      r_5[i] <- p_t[i, 3] * spawners[i] * growth
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
    time = seq_along(ind_data),
    x = x[ind_data],
    returns = returns[ind_data],
    spawners = spawners[ind_data],
    r_3 = r_3[ind_data],
    r_4 = r_4[ind_data],
    r_5 = r_5[ind_data],
    recruits = recruits[ind_data]
  )
}
