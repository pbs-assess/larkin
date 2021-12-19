#' Forecast Single Value Using Previous Values Only
#'
#' @param index [numeric()]
#' @param data [list()]
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
forecast_single_value <- function (index,
                                   data,
                                   chains = 1,
                                   step_size = 0.01,
                                   iter_warmup = 250,
                                   iter_sampling = 750,
                                   ...) {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_list(data, c("double", "integer"))
  # checkmate::assert_choice(model, c("larkin", "ricker"))

  # Truncate data --------------------------------------------------------------

  data$N <- index
  data$recruits <- data$recruits[seq_len(index)]
  data$spawners <- data$spawners[seq_len(index)]
  data$environs <- data$environs[seq_len(min(nrow(data$environs), index)), ]


  # Create model object --------------------------------------------------------

  mod <- cmdstanr::cmdstan_model(
    system.file("stan", paste0("forecast", ".stan"), package = "larkin"),
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

  # Extract forecast -----------------------------------------------------------

  # Placate R-CMD-check
  forecast <- NULL

  # Return forecast
  fit %>%
    tidybayes::spread_draws(forecast) %>%
    tidybayes::summarise_draws() %>%
    dplyr::ungroup() %>%
    dplyr::mutate(time = as.integer(index)) %>%
    dplyr::mutate(observed = data$recruits[index]) %>%
    dplyr::mutate(forecast = .data$mean) %>%
    dplyr::select(
      .data$time,
      .data$observed,
      .data$forecast,
      .data$median:.data$ess_tail
    )
}

#' Summarise Posterior Draws
#'
#' @param x [cmdstanr::sample()] model fit object
#' @param data [list()] model data
#'
#' @return [list()]
#' @export
#'
summarise_posterior_draws <- function (x, data) {
  # Placate R-CMD-check
  index <- NULL
  # lp__
  lp__ <- tidybayes::spread_draws(x, lp__) %>%
    tidybayes::summarise_draws() %>%
    dplyr::ungroup()
  # Productivity
  alpha <- tidybayes::spread_draws(x, alpha) %>%
    tidybayes::summarise_draws() %>%
    dplyr::ungroup()
  # Denisty dependence
  beta <- tidybayes::spread_draws(x, beta[index]) %>%
    tidybayes::summarise_draws() %>%
    dplyr::ungroup()
  # Environmental
  if (data$G == 0) {
    gamma <- tibble::tibble()
  } else {
    gamma <- tidybayes::spread_draws(x, gamma[index]) %>%
      tidybayes::summarise_draws() %>%
      dplyr::ungroup()
  }
  # Standard deviation
  sigma <- tidybayes::spread_draws(x, sigma) %>%
    tidybayes::summarise_draws() %>%
    dplyr::ungroup()
  # Return list
  list(
    lp__ = lp__,
    alpha = alpha,
    beta = beta,
    gamma = gamma,
    sigma = sigma
  )
}

