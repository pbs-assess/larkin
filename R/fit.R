#' Fit A Larkin Model Via CmdStanR
#'
#' @param data [list()]
#' @param chains [integer()] number of chains
#' @param step_size [integer()] initial step size
#' @param iter_warmup [integer()] number of warmup iterations
#' @param iter_sampling [integer()] number of sampling iterations
#' @param use_reduce_sum [logical()] use within chain parallel threading
#' @param threads_per_chain [integer()] number of threads per chain
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
                    use_reduce_sum = FALSE,
                    threads_per_chain = 8,
                    ...) {

  # Check arguments ------------------------------------------------------------


  # Fit the model via cmdstanr -------------------------------------------------

  # Create model object
  mod <- cmdstanr::cmdstan_model(
    system.file("stan", "larkin.stan", package = "larkin"),
    include_path = system.file("stan", package = "larkin"),
    cpp_options = list(stan_threads = TRUE))
  # Fit the model
  cmdfit <- mod$sample(
    data = data,
    chains = chains,
    step_size = step_size,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    threads_per_chain = threads_per_chain,
    ...
  )
  # Convert to stanfit
  samples <- rstan::read_stan_csv(cmdfit$output_files())

  # Return mmmfit object -------------------------------------------------------

  structure(list(
    data = data,
    samples = samples,
    cmdfit = cmdfit),
    class = "larkin")
}
