#' Forecast From A Ricker Or Larkin Model
#'
#' @param index [integer()] [vector()]
#' @param recruits [numeric()] [vector()]
#' @param spawners [numeric()] [vector()]
#' @param environs [numeric()] [vector()]
#' @param timevary [logical()]
#' @param mu_alpha [numeric()]
#' @param mu_beta [numeric()] [vector()]
#' @param mu_gamma [numeric()] [vector()]
#' @param mu_sigma [numeric()]
#' @param sd_alpha [numeric()]
#' @param sd_beta [numeric()]
#' @param sd_gamma [numeric()] [vector()]
#' @param sd_sigma [numeric()]
#' @param sd_phi [numeric()]
#' @param id_cols [list()] of colname-value pairs
#' @param buffer [logical()]
#' @param cores [numeric()]
#' @param chains [numeric()]
#' @param step_size [numeric()]
#' @param adapt_delta [numeric()]
#' @param iter_warmup [integer()]
#' @param iter_sampling [integer()]
#' @param ... additional arguments to pass to [cmdstanr::sample()]
#'
#' @return [tibble::tibble()]
#' @export
#'
forecast <- function (index,
                      recruits,
                      spawners,
                      environs = matrix(0, 0, 0),
                      timevary = FALSE,
                      mu_alpha,
                      mu_beta,
                      mu_gamma = numeric(0),
                      mu_sigma,
                      sd_alpha,
                      sd_beta,
                      sd_gamma = numeric(0),
                      sd_sigma,
                      sd_phi = numeric(0),
                      id_cols = NULL,
                      buffer = TRUE,
                      cores = NULL,
                      chains = 3,
                      step_size = 0.01,
                      adapt_delta = 0.9,
                      iter_warmup = 250,
                      iter_sampling = 750,
                      ...) {

  # Check arguments ------------------------------------------------------------

  # max(index) not more than length(recruits) + 1


  # Define data list -----------------------------------------------------------

  data <- list(
    N = length(recruits),
    B = length(mu_beta),
    G = ncol(as.matrix(environs)),
    recruits = recruits,
    spawners = spawners,
    environs = as.matrix(environs),
    timevary = as.numeric(timevary),
    mu_alpha = mu_alpha,
    mu_beta = as.array(mu_beta),
    mu_gamma = as.array(mu_gamma),
    mu_sigma = mu_sigma,
    sd_alpha = sd_alpha,
    sd_beta = as.array(sd_beta),
    sd_gamma = as.array(sd_gamma),
    sd_sigma = sd_sigma,
    sd_phi = sd_phi,
    fudge = 1e-12
  )

  # Generate forecasts ---------------------------------------------------------

  if (length(index) > 1) {
    # Create forecasts
    if (.Platform$OS.type == "unix") {
      forecasts <- parallel::mclapply(
        X = index,
        FUN = forecast,
        recruits = recruits,
        spawners = spawners,
        environs = environs,
        timevary = timevary,
        mu_alpha = mu_alpha,
        mu_beta = mu_beta,
        mu_gamma = mu_gamma,
        mu_sigma = mu_sigma,
        sd_alpha = sd_alpha,
        sd_beta = sd_beta,
        sd_gamma = sd_gamma,
        sd_sigma = sd_sigma,
        sd_phi = sd_phi,
        id_cols = id_cols,
        buffer = FALSE,
        chains = chains,
        step_size = step_size,
        adapt_delta = adapt_delta,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        ...,
        mc.cores = cores
      )
    } else {
      forecasts <- lapply(
        X = index,
        FUN = forecast,
        recruits = recruits,
        spawners = spawners,
        environs = environs,
        timevary = timevary,
        mu_alpha = mu_alpha,
        mu_beta = mu_beta,
        mu_gamma = mu_gamma,
        mu_sigma = mu_sigma,
        sd_alpha = sd_alpha,
        sd_beta = sd_beta,
        sd_gamma = sd_gamma,
        sd_sigma = sd_sigma,
        sd_phi = sd_phi,
        id_cols = id_cols,
        buffer = FALSE,
        chains = chains,
        step_size = step_size,
        adapt_delta = adapt_delta,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        ...
      )
    }
    # Assemble forecasts
    forecasts <- dplyr::bind_rows(forecasts) %>%
      dplyr::ungroup()
  } else if (length(index) == 1) {
    # Define data list
    data <- list(
      N = index,
      B = length(mu_beta),
      G = ncol(as.matrix(environs)),
      recruits = recruits[seq_len(min(nrow(recruits), index))],
      spawners = spawners[seq_len(min(nrow(spawners), index))],
      environs = as.matrix(environs[seq_len(min(nrow(environs), index)), ]),
      timevary = as.numeric(timevary),
      mu_alpha = mu_alpha,
      mu_beta = as.array(mu_beta),
      mu_gamma = as.array(mu_gamma),
      mu_sigma = mu_sigma,
      sd_alpha = sd_alpha,
      sd_beta = as.array(sd_beta),
      sd_gamma = as.array(sd_gamma),
      sd_sigma = sd_sigma,
      sd_phi = sd_phi,
      fudge = 1e-12
    )
    # Create model object
    mod <- cmdstanr::cmdstan_model(
      stan_file = here::here("stan", "forecast.stan"),
      include_path = here::here("stan")
    )
    # Create forecasts
    forecasts <- mod$sample(
      data = data,
      chains = chains,
      step_size = step_size,
      adapt_delta = adapt_delta,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      ...
    ) %>%
      tidybayes::spread_draws(forecast) %>%
      tidybayes::summarise_draws() %>%
      dplyr::ungroup() %>%
      # dplyr::select(.data$mean:.data$ess_tail) %>%
      dplyr::mutate(forecast = .data$mean) %>%
      dplyr::relocate(.data$forecast, .before = 1) %>%
      dplyr::mutate(index = index) %>%
      dplyr::relocate(.data$index, .before = 1) %>%
      dplyr::bind_cols(id_cols) %>%
      dplyr::ungroup()
    # TODO: Make a list including forecasts and parameter summaries?
  } else {
    stop("index must have length >= 1")
  }

  # Format forecasts -----------------------------------------------------------

  # Buffer forecasts
  if (buffer) {
    forecasts <- forecasts %>%
      tibble::add_row(index = seq_len(min(index) - 1), .before = 1)
  }
  # Arrange forecasts
  forecasts <- dplyr::arrange(forecasts, .data$index)

  # Return forecasts -----------------------------------------------------------

  return(forecasts)
}











# New above here ---------------------------------------------------------------






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
