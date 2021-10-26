#' Forecast Single Value Using Previous Values Only
#'
#' @param index [numeric()]
#' @param data [list()]
#' @param model [character()] either \code{"larkin"} or \code{"ricker"}
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
                                   model = "larkin",
                                   chains = 1,
                                   step_size = 0.01,
                                   iter_warmup = 250,
                                   iter_sampling = 750,
                                   ...) {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_list(data, c("double", "integer"))
  checkmate::assert_choice(model, c("larkin", "ricker"))

  # Truncate data --------------------------------------------------------------

  data$T <- index
  data$recruits <- data$recruits[seq_len(index)]
  data$spawners <- data$spawners[seq_len(index)]

  # Create model object --------------------------------------------------------

  mod <- cmdstanr::cmdstan_model(
    system.file("stan", paste0(model, "_forecast.stan"), package = "larkin"),
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

  # Extract forecast -----------------------------------------------------------

  cmdfit$summary() %>%
    dplyr::filter(.data$variable %in% c("forecast")) %>%
    dplyr::mutate(
      time = as.integer(index),
      observed = data$recruits[index],
      forecast = .data$mean
    ) %>%
    dplyr::select(
      .data$time,
      .data$observed,
      .data$forecast,
      .data$sd,
      .data$rhat,
      .data$ess_bulk,
      .data$ess_tail
    )
}
