#' Forecast From A Ricker Or Larkin Model
#'
#' @param index [integer()] [vector()]
#' @param recruits [numeric()] [vector()]
#' @param spawners [numeric()] [vector()]
#' @param environs [numeric()] [vector()]
#' @param timevary [logical()]
#' @param prior_mean_alpha [numeric()]
#' @param prior_mean_beta [numeric()] [vector()]
#' @param prior_mean_gamma [numeric()] [vector()]
#' @param prior_mean_sigma [numeric()]
#' @param prior_mean_omega [numeric()]
#' @param prior_sd_alpha [numeric()]
#' @param prior_sd_beta [numeric()]
#' @param prior_sd_gamma [numeric()] [vector()]
#' @param prior_sd_sigma [numeric()]
#' @param prior_sd_omega [numeric()]
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
#' @importFrom rlang .data
#' @importFrom rlang :=
#'
#' @return [tibble::tibble()]
#' @export
#'
forecast <- function (index,
                      recruits,
                      spawners,
                      environs = matrix(0, 0, 0),
                      timevary = FALSE,
                      prior_mean_alpha,
                      prior_mean_beta,
                      prior_mean_gamma = numeric(0),
                      prior_mean_sigma,
                      prior_mean_omega = 0,
                      prior_sd_alpha,
                      prior_sd_beta,
                      prior_sd_gamma = numeric(0),
                      prior_sd_sigma,
                      prior_sd_omega = 0,
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
  # use timevary as a check on prior_mean_omega and prior_sd_omega

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
        prior_mean_alpha = prior_mean_alpha,
        prior_mean_beta = prior_mean_beta,
        prior_mean_gamma = prior_mean_gamma,
        prior_mean_sigma = prior_mean_sigma,
        prior_mean_omega = prior_mean_omega,
        prior_sd_alpha = prior_sd_alpha,
        prior_sd_beta = prior_sd_beta,
        prior_sd_gamma = prior_sd_gamma,
        prior_sd_sigma = prior_sd_sigma,
        prior_sd_omega = prior_sd_omega,
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
        prior_mean_alpha = prior_mean_alpha,
        prior_mean_beta = prior_mean_beta,
        prior_mean_gamma = prior_mean_gamma,
        prior_mean_sigma = prior_mean_sigma,
        prior_mean_omega = prior_mean_omega,
        prior_sd_alpha = prior_sd_alpha,
        prior_sd_beta = prior_sd_beta,
        prior_sd_gamma = prior_sd_gamma,
        prior_sd_sigma = prior_sd_sigma,
        prior_sd_omega = prior_sd_omega,
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
      B = length(prior_mean_beta),
      G = ncol(as.matrix(environs)),
      recruits = recruits[seq_len(min(nrow(recruits), index))],
      spawners = spawners[seq_len(min(nrow(spawners), index))],
      environs = as.matrix(environs[seq_len(min(nrow(environs), index)), ]),
      prior_mean_alpha = prior_mean_alpha,
      prior_mean_beta = as.array(prior_mean_beta),
      prior_mean_gamma = as.array(prior_mean_gamma),
      prior_mean_sigma = prior_mean_sigma,
      prior_mean_omega = prior_mean_omega,
      prior_sd_alpha = prior_sd_alpha,
      prior_sd_beta = as.array(prior_sd_beta),
      prior_sd_gamma = as.array(prior_sd_gamma),
      prior_sd_sigma = prior_sd_sigma,
      prior_sd_omega = prior_sd_omega,
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
    # TODO: For now, min ESS and max R-hat
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
