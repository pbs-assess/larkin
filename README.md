# larkin
Forecast from Larkin and Ricker stock-recruitment models

<!-- badges: start -->
[![R-CMD-check](https://github.com/luke-a-rogers/larkin/workflows/R-CMD-check/badge.svg)](https://github.com/luke-a-rogers/larkin/actions)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
<!-- badges: end -->

## Installation

1. Install the R package `cmdstanr` (see <https://mc-stan.org/cmdstanr/index.html>).

``` r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
```

2. Install CmdStan (see <https://mc-stan.org/cmdstanr/articles/cmdstanr.html>).

``` r
cmdstanr::check_cmdstan_toolchain()
cmdstanr::install_cmdstan(cores = parallel::detectCores())
```

3. Install `larkin`.

``` r
remotes::install_github("pbs-assess/larkin")
```

