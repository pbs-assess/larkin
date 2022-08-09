#' Forecast From A Ricker Or Larkin Model
#'
#' @param data [data.frame()]
#' @param index [integer()] [vector()]
#' @param recruits [character()]
#' @param spawners [character()]
#' @param environs [character()] [vector()]
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
#' @param id_cols [character()] colnames in \code{data}
#' @param id_vals [list()] of name-value pairs
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
forecast <- function (data,
                      index,
                      recruits = "recruits",
                      spawners = "spawners",
                      environs = character(0),
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
                      id_vals = NULL,
                      buffer = TRUE,
                      cores = 1,
                      chains = 3,
                      step_size = 0.01,
                      adapt_delta = 0.9,
                      iter_warmup = 250,
                      iter_sampling = 750,
                      ...) {

  # Check arguments ------------------------------------------------------------

  # max(index) not more than length(recruits) + 1
  # use timevary as a check on prior_mean_omega and prior_sd_omega

  # Define id columns and values -----------------------------------------------

  # Define id columns
  id_columns <- NULL
  if (length(id_cols) > 0) {
    id_columns <- data[1L, id_cols]
  }
  # Define id values
  id_values <- NULL
  if (length(id_vals) > 0) {
    id_values <- tibble::as_tibble(id_vals)
  }

  # Generate forecasts ---------------------------------------------------------

  if (length(index) > 1) {
    # Create forecasts
    if (.Platform$OS.type == "unix") {
      outputs <- parallel::mclapply(
        X = index,
        FUN = larkin::forecast,
        data = data,
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
        id_vals = id_vals,
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
      outputs <- lapply(
        X = index,
        FUN = larkin::forecast,
        data = data,
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
        id_vals = id_vals,
        buffer = FALSE,
        chains = chains,
        step_size = step_size,
        adapt_delta = adapt_delta,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        ...
      )
    }
    # Transpose outputs
    outputs <- purrr::transpose(outputs)
    # Define forecasts
    forecasts <- dplyr::bind_rows(outputs$forecasts) %>%
      # dplyr::arrange(.data$index) %>%
      dplyr::ungroup()
    # Define lp__
    lp__ <- dplyr::bind_rows(outputs$lp__) %>%
      # dplyr::arrange(.data$index) %>%
      dplyr::ungroup()
    # Define alpha
    alpha <- dplyr::bind_rows(outputs$alpha) %>%
      # dplyr::arrange(.data$index, .data$n) %>%
      dplyr::ungroup()
    # Define beta
    beta <- dplyr::bind_rows(outputs$beta) %>%
      # dplyr::arrange(.data$index, .data$b) %>%
      dplyr::ungroup()
    # Define gamma
    gamma <- dplyr::bind_rows(outputs$gamma) %>%
      dplyr::ungroup()
    # Define sigma
    sigma <- dplyr::bind_rows(outputs$sigma) %>%
      # dplyr::arrange(.data$index) %>%
      dplyr::ungroup()
    # Define omega
    omega <- dplyr::bind_rows(outputs$omega) %>%
      dplyr::ungroup()
    # Define output
    output <- list(
      forecasts = forecasts,
      lp__ = lp__,
      alpha = alpha,
      beta = beta,
      gamma = gamma,
      sigma = sigma,
      omega = omega
    )
  } else if (length(index) == 1) {
    # Define observations
    recruits <- dplyr::pull(data, recruits)
    spawners <- dplyr::pull(data, spawners)
    if (length(which(environs %in% colnames(data))) > 0) {
      environs <- as.matrix(dplyr::select(data, environs))
    } else {
      environs <- matrix(0, 0, 0)
    }
    # Define Stan data
    stan_data <- list(
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
      stan_file = system.file(
        "stan", "forecast.stan",
        package = "larkin",
        mustWork = TRUE
      ),
      include_path = system.file(
        "stan",
        package = "larkin",
        mustWork = TRUE
      )
    )
    # Create samples
    samples <- mod$sample(
      data = stan_data,
      chains = chains,
      step_size = step_size,
      adapt_delta = adapt_delta,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      ...
    )
    # Placate R-CMD-check
    n <- NULL
    b <- NULL
    g <- NULL
    o <- NULL
    # Define summaries
    summaries <- samples %>%
      tidybayes::summarise_draws() %>%
      dplyr::ungroup()
    # Define maximum rhat
    max_rhat <- max(summaries$rhat, na.rm = TRUE)
    # Define minumum ess
    min_ess_bulk <- min(summaries$ess_bulk, na.rm = TRUE)
    min_ess_tail <- min(summaries$ess_tail, na.rm = TRUE)
    # Define forecasts
    forecasts <- samples %>%
      tidybayes::spread_draws(forecast) %>%
      tidybayes::summarise_draws() %>%
      dplyr::ungroup() %>%
      dplyr::select(.data$mean:.data$q95) %>%
      # dplyr::mutate(forecast = .data$mean) %>%
      dplyr::mutate(observed = recruits[index]) %>%
      dplyr::mutate(index = index) %>%
      dplyr::mutate(max_rhat = max_rhat) %>%
      dplyr::mutate(min_ess_bulk = min_ess_bulk) %>%
      dplyr::mutate(min_ess_tail = min_ess_tail) %>%
      # dplyr::relocate(.data$forecast, .before = 1) %>%
      dplyr::relocate(.data$observed, .before = 1) %>%
      dplyr::relocate(.data$index, .before = 1) %>%
      dplyr::bind_cols(id_values) %>%
      dplyr::relocate(colnames(id_values), .before = 1) %>%
      dplyr::bind_cols(id_columns) %>%
      dplyr::relocate(colnames(id_columns), .before = 1) %>%
      dplyr::ungroup()
    # Define lp__
    lp__ <- samples %>%
      tidybayes::spread_draws(lp__) %>%
      tidybayes::summarise_draws() %>%
      dplyr::mutate(index = index) %>%
      dplyr::relocate(.data$index, .before = 1) %>%
      dplyr::bind_cols(id_values) %>%
      dplyr::relocate(colnames(id_values), .before = 1) %>%
      dplyr::bind_cols(id_columns) %>%
      dplyr::relocate(colnames(id_columns), .before = 1) %>%
      dplyr::ungroup()
    # Define alpha
    alpha <- samples %>%
      tidybayes::spread_draws(alpha[n]) %>%
      tidybayes::summarise_draws() %>%
      dplyr::mutate(index = index) %>%
      dplyr::relocate(.data$index, .before = 1) %>%
      dplyr::mutate(index = index) %>%
      dplyr::relocate(.data$index, .before = 1) %>%
      dplyr::bind_cols(id_values) %>%
      dplyr::relocate(colnames(id_values), .before = 1) %>%
      dplyr::bind_cols(id_columns) %>%
      dplyr::relocate(colnames(id_columns), .before = 1) %>%
      dplyr::ungroup()
    # Define beta
    beta <- samples %>%
      tidybayes::spread_draws(beta[b]) %>%
      tidybayes::summarise_draws() %>%
      dplyr::mutate(index = index) %>%
      dplyr::relocate(.data$index, .before = 1) %>%
      dplyr::bind_cols(id_values) %>%
      dplyr::relocate(colnames(id_values), .before = 1) %>%
      dplyr::bind_cols(id_columns) %>%
      dplyr::relocate(colnames(id_columns), .before = 1) %>%
      dplyr::ungroup()
    # Define gamma
    if (length(prior_mean_gamma) > 0) {
      gamma <- samples %>%
        tidybayes::spread_draws(gamma[g]) %>%
        tidybayes::summarise_draws() %>%
        dplyr::mutate(index = index) %>%
        dplyr::relocate(.data$index, .before = 1) %>%
        dplyr::bind_cols(id_values) %>%
        dplyr::relocate(colnames(id_values), .before = 1) %>%
        dplyr::bind_cols(id_columns) %>%
        dplyr::relocate(colnames(id_columns), .before = 1) %>%
        dplyr::ungroup()
    } else {
      gamma <- tibble::tibble()
    }
    # Define sigma
    sigma <- samples %>%
      tidybayes::spread_draws(sigma) %>%
      tidybayes::summarise_draws() %>%
      dplyr::mutate(index = index) %>%
      dplyr::relocate(.data$index, .before = 1) %>%
      dplyr::bind_cols(id_values) %>%
      dplyr::relocate(colnames(id_values), .before = 1) %>%
      dplyr::bind_cols(id_columns) %>%
      dplyr::relocate(colnames(id_columns), .before = 1) %>%
      dplyr::ungroup()
    # Define omega
    if (prior_mean_omega > 0 | prior_sd_omega > 0) {
      omega <- samples %>%
        tidybayes::spread_draws(omega[o]) %>%
        tidybayes::summarise_draws() %>%
        dplyr::mutate(index = index) %>%
        dplyr::relocate(.data$index, .before = 1) %>%
        dplyr::bind_cols(id_values) %>%
        dplyr::relocate(colnames(id_values), .before = 1) %>%
        dplyr::bind_cols(id_columns) %>%
        dplyr::relocate(colnames(id_columns), .before = 1) %>%
        dplyr::ungroup() %>%
        dplyr::select(!o)
    } else {
      omega <- tibble::tibble()
    }
    # Define output
    output <- list(
      forecasts = forecasts,
      lp__ = lp__,
      alpha = alpha,
      beta = beta,
      gamma = gamma,
      sigma = sigma,
      omega = omega
    )
  } else {
    stop("index must have length >= 1")
  }

  # Buffer forecasts -----------------------------------------------------------

  if (buffer) {
    output$forecasts <- output$forecasts %>%
      tibble::add_row(
        index = seq_len(min(index) - 1),
        id_columns,
        id_values,
        .before = 1
      )
  }

  # Return output --------------------------------------------------------------

  return(output)
}
