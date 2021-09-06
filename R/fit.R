#' Fit A Larkin Model Via CmdStanR
#'
#' @param data [list()]
#' @param chains [integer()] number of chains
#' @param step_size [integer()] initial step size
#' @param iter_warmup [integer()] number of warmup iterations
#' @param iter_sampling [integer()] number of sampling iterations
#' @param ... additional arguments to pass to \code{$sample()} method
#'
#' @return [larkin()]
#' @export
#'
larkin <- function (data,
                    chains = 1,
                    step_size = 0.1,
                    iter_warmup = 250,
                    iter_sampling = 750,
                    ...) {

  # Check arguments ------------------------------------------------------------


  # Fit the model via cmdstanr -------------------------------------------------

  # Create model object
  mod <- cmdstanr::cmdstan_model(
    system.file("stan", "larkin.stan", package = "larkin"),
    include_path = system.file("stan", package = "larkin")
    # cpp_options = list(stan_threads = TRUE)
  )
  # Fit the model
  cmdfit <- mod$sample(
    data = data,
    chains = chains,
    step_size = step_size,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    ...
  )

  # Return larkin object -------------------------------------------------------

  structure(list(
    data = data,
    output = cmdfit$output_files(),
    summary = cmdfit$summary()),
    class = "larkin")
}

#' Fit A Ricker Model Via CmdStanR
#'
#' @param data [list()]
#' @param chains [integer()] number of chains
#' @param step_size [integer()] initial step size
#' @param iter_warmup [integer()] number of warmup iterations
#' @param iter_sampling [integer()] number of sampling iterations
#' @param ... additional arguments to pass to \code{$sample()} method
#'
#' @return [ricker()]
#' @export
#'
ricker <- function (data,
                    chains = 1,
                    step_size = 0.1,
                    iter_warmup = 250,
                    iter_sampling = 750,
                    ...) {

  # Check arguments ------------------------------------------------------------


  # Fit the model via cmdstanr -------------------------------------------------

  # Create model object
  mod <- cmdstanr::cmdstan_model(
    system.file("stan", "ricker.stan", package = "larkin"),
    include_path = system.file("stan", package = "larkin")
    # cpp_options = list(stan_threads = TRUE)
  )
  # Fit the model
  cmdfit <- mod$sample(
    data = data,
    chains = chains,
    step_size = step_size,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    ...
  )

  # Return ricker object -------------------------------------------------------

  structure(list(
    data = data,
    output = cmdfit$output_files(),
    summary = cmdfit$summary()),
    class = "ricker")
}
