% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{forecast}
\alias{forecast}
\title{Forecast From A Ricker Or Larkin Model}
\usage{
forecast(
  data,
  index = NULL,
  recruits = "recruits",
  spawners = "spawners",
  p_prime = c(0.003, 0.917, 0.08),
  environs = character(0),
  run_stan = TRUE,
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
  cores = 1,
  chains = 3,
  step_size = 0.01,
  adapt_delta = 0.9,
  iter_warmup = 250,
  iter_sampling = 750,
  ...
)
}
\arguments{
\item{data}{\code{\link[=data.frame]{data.frame()}} input data}

\item{index}{\code{\link[=integer]{integer()}} \code{\link[=vector]{vector()}}}

\item{recruits}{\code{\link[=character]{character()}}}

\item{spawners}{\code{\link[=character]{character()}}}

\item{p_prime}{\code{\link[=character]{character()}} \code{\link[=vector]{vector()}}}

\item{environs}{\code{\link[=character]{character()}} \code{\link[=vector]{vector()}}}

\item{run_stan}{\code{\link[=logical]{logical()}}}

\item{prior_mean_alpha}{\code{\link[=numeric]{numeric()}}}

\item{prior_mean_beta}{\code{\link[=numeric]{numeric()}} \code{\link[=vector]{vector()}}}

\item{prior_mean_gamma}{\code{\link[=numeric]{numeric()}} \code{\link[=vector]{vector()}}}

\item{prior_mean_sigma}{\code{\link[=numeric]{numeric()}}}

\item{prior_mean_omega}{\code{\link[=numeric]{numeric()}}}

\item{prior_sd_alpha}{\code{\link[=numeric]{numeric()}}}

\item{prior_sd_beta}{\code{\link[=numeric]{numeric()}}}

\item{prior_sd_gamma}{\code{\link[=numeric]{numeric()}} \code{\link[=vector]{vector()}}}

\item{prior_sd_sigma}{\code{\link[=numeric]{numeric()}}}

\item{prior_sd_omega}{\code{\link[=numeric]{numeric()}}}

\item{id_cols}{\code{\link[=character]{character()}} colnames in \code{data}}

\item{id_vals}{\code{\link[=list]{list()}} of name-value pairs}

\item{cores}{\code{\link[=numeric]{numeric()}}}

\item{chains}{\code{\link[=numeric]{numeric()}}}

\item{step_size}{\code{\link[=numeric]{numeric()}}}

\item{adapt_delta}{\code{\link[=numeric]{numeric()}}}

\item{iter_warmup}{\code{\link[=integer]{integer()}}}

\item{iter_sampling}{\code{\link[=integer]{integer()}}}

\item{...}{additional arguments to pass to \code{\link[cmdstanr:model-method-sample]{cmdstanr::sample()}}}
}
\value{
\code{\link[tibble:tibble]{tibble::tibble()}}
}
\description{
Forecast From A Ricker Or Larkin Model
}
