#' Absolute Lag-1-Scaled Error
#'
#' @param x [numeric()] [vector()] actual values
#' @param y [numeric()] [vector()] forecast values
#'
#' @return [numeric()] [vector()] absolute scaled errors
#' @export
#'
#' @examples
#' x <- seq(0, 10, 0.1)
#' y <- seq(10, 0, -0.1)
#' m <- abs_lag1_scaled_error(x, y)
#' plot(m, type = "l")
#'
abs_lag1_scaled_error <- function (x, y) {
  # Check arguments
  checkmate::assert_numeric(x, finite = TRUE, min.len = 1)
  checkmate::assert_numeric(y, finite = TRUE, min.len = 1)
  checkmate::assert_vector(x, strict = TRUE)
  checkmate::assert_vector(y, strict = TRUE)
  checkmate::assert_true(length(x) == length(y))
  # Return value
  return(abs(y - x) / mean(abs(x - dplyr::lag(x)), na.rm = TRUE))
}

#' Absolute Mean-Scaled Error
#'
#' @param x [numeric()] [vector()] observed values
#' @param y [numeric()] [vector()] forecast values
#'
#' @return [numeric()] [vector()] absolute scaled errors
#' @export
#'
#' @examples
#' x <- seq(0, 10, 0.1)
#' y <- seq(10, 0, -0.1)
#' m <- abs_mean_scaled_error(x, y)
#' plot(m, type = "l")
#'
abs_mean_scaled_error <- function (x, y) {
  # Check arguments
  checkmate::assert_numeric(x, finite = TRUE, min.len = 1)
  checkmate::assert_numeric(y, finite = TRUE, min.len = 1)
  checkmate::assert_vector(x, strict = TRUE)
  checkmate::assert_vector(y, strict = TRUE)
  checkmate::assert_true(length(x) == length(y))
  # Return value
  return(abs(y - x) / mean(abs(x - mean(x, na.rm = TRUE)), na.rm = TRUE))
}

#' Absolute Percentage Error
#'
#' @param x [numeric()] [vector()] actual values
#' @param y [numeric()] [vector()] forecast values
#'
#' @return [numeric()] [vector()] absolute percentage errors
#' @export
#'
#' @examples
#' x <- seq(0, 10, 0.1)
#' y <- seq(10, 0, -0.1)
#' m <- abs_pct_error(x, y)
#' plot(m, type = "l")
#'
abs_pct_error <- function (x, y) {
  # Check arguments
  checkmate::assert_numeric(x, finite = TRUE, min.len = 1)
  checkmate::assert_numeric(y, finite = TRUE, min.len = 1)
  checkmate::assert_vector(x, strict = TRUE)
  checkmate::assert_vector(y, strict = TRUE)
  checkmate::assert_true(length(x) == length(y))
  # Return value
  return(abs((y - x) / x))
}

#' Arctangent Absolute Percentage Error
#'
#' @param x [numeric()] [vector()] actual values
#' @param y [numeric()] [vector()] forecast values
#'
#' @return [numeric()] [vector()] arctangent absolute percentage errors
#' @export
#'
#' @references Kim and Kim 2016 doi.org/10.1016/j.ijforecast.2015.12.003
#'
#' @examples
#' x <- seq(0, 10, 0.1)
#' y <- seq(10, 0, -0.1)
#' m <- arctan_abs_pct_error(x, y)
#' plot(m, type = "l")
#'
arctan_abs_pct_error <- function (x, y) {
  # Check arguments
  checkmate::assert_numeric(x, finite = TRUE, min.len = 1)
  checkmate::assert_numeric(y, finite = TRUE, min.len = 1)
  checkmate::assert_vector(x, strict = TRUE)
  checkmate::assert_vector(y, strict = TRUE)
  checkmate::assert_true(length(x) == length(y))
  # Return
  return(atan(abs_pct_error(x, y)))
}

#' Mean Arctangent Absolute Percentage Error
#'
#' @param x [numeric()] [vector()] actual values
#' @param y [numeric()] [vector()] forecast values
#'
#' @return [numeric()] mean arctangent absolute percentage error
#' @export
#'
#' @references Kim and Kim 2016 doi.org/10.1016/j.ijforecast.2015.12.003
#'
#' @examples
#' maape(1:10, 10:1)
#'
maape <- function (x, y) {
  # Check arguments
  checkmate::assert_numeric(x, finite = TRUE, min.len = 1)
  checkmate::assert_numeric(y, finite = TRUE, min.len = 1)
  checkmate::assert_vector(x, strict = TRUE)
  checkmate::assert_vector(y, strict = TRUE)
  checkmate::assert_true(length(x) == length(y))
  # Return value
  return(mean(arctan_abs_pct_error(x, y), na.rm = TRUE))
}

#' Mean Absolute Error
#'
#' @param x [numeric()] [vector()] actual values
#' @param y [numeric()] [vector()] forecast values
#'
#' @return [numeric()] mean absolute error
#' @export
#'
#' @examples
#' mae(1:10, 10:1)
#'
mae <- function (x, y) {
  # Check arguments
  checkmate::assert_numeric(x, finite = TRUE, min.len = 1)
  checkmate::assert_numeric(y, finite = TRUE, min.len = 1)
  checkmate::assert_vector(x, strict = TRUE)
  checkmate::assert_vector(y, strict = TRUE)
  checkmate::assert_true(length(x) == length(y))
  # Return value
  return(mean(abs(y - x), na.rm = TRUE))
}

#' Mean Absolute Mean-Scaled Error
#'
#' @param x [numeric()][vector()] observed values
#' @param y [numeric()][vector()] forecast values
#'
#' @return [numeric()] scalar
#' @export
#'
#' @examples
#' x <- rnorm(100, 0, 1)
#' y1 <- x
#' mamse(x, y1)
#'
#' y2 <- rep(mean(x), length(x))
#' mamse(x, y2)
#'
#' mamse(c(x, NA), c(y1, NA))
#'
mamse <- function(x, y) {
  # Check arguments
  checkmate::assert_numeric(x, finite = TRUE, min.len = 1)
  checkmate::assert_numeric(y, finite = TRUE, min.len = 1)
  checkmate::assert_vector(x, strict = TRUE)
  checkmate::assert_vector(y, strict = TRUE)
  checkmate::assert_true(length(x) == length(y))
  # Return value
  return(mean(abs_mean_scaled_error(x, y), na.rm = TRUE))
}

#' Mean Absolute Percentage Error
#'
#' @param x [numeric()][vector()] observed values
#' @param y [numeric()][vector()] forecast values
#'
#' @return [numeric()] mean absolute percentage error
#' @export
#'
#' @examples
#' mape(1:10, 10:1)
#'
mape <- function (x, y) {
  # Check arguments
  checkmate::assert_numeric(x, finite = TRUE, min.len = 1)
  checkmate::assert_numeric(y, finite = TRUE, min.len = 1)
  checkmate::assert_vector(x, strict = TRUE)
  checkmate::assert_vector(y, strict = TRUE)
  checkmate::assert_true(length(x) == length(y))
  # Return value
  return(mean(abs_pct_error(x, y), na.rm = TRUE))
}

#' Mean Absolute Lag-1-Scaled Error
#'
#' @param x [numeric()][vector()] observed values
#' @param y [numeric()][vector()] forecast values
#'
#' @return [numeric()] scalar
#' @export
#'
#' @examples
#' x <- rnorm(100, 0, 1)
#' y1 <- x
#' mase(x, y1)
#'
#' y2 <- dplyr::lag(x)
#' mase(x, y2)
#'
#' mase(c(x, NA), c(y1, NA))
#'
mase <- function (x, y) {
  # Check arguments
  checkmate::assert_numeric(x, finite = TRUE, min.len = 1)
  checkmate::assert_numeric(y, finite = TRUE, min.len = 1)
  checkmate::assert_vector(x, strict = TRUE)
  checkmate::assert_vector(y, strict = TRUE)
  checkmate::assert_true(length(x) == length(y))
  # Return value
  return(mean(abs_lag1_scaled_error(x, y), na.rm = TRUE))
}

#' Metric Wrapper
#'
#' @param x [numeric()][matrix()] of two columns
#' @param fun [function()] taking two [numeric()][vector()] arguments and
#'   returning one scalar [numeric()] value
#' @param ... additional arguments to pass to \code{fun()}
#'
#' @return [numeric()] scalar
#' @export
#'
#' @examples
#' x <- rnorm(100, 0, 1)
#' x1 <- matrix(c(x, x), ncol = 2)
#' matric(x1, mase)
#'
#' x2 <- rbind(x1, matrix(c(NA_real_, NA_real_), ncol = 2))
#' matric(x2, mase)
#'
matric <- function (x, fun, ...) {
  # Check arguments
  checkmate::assert_matrix(x, mode = "numeric", ncols = 2L)
  checkmate::assert_function(fun, null.ok = FALSE)
  # Compute value
  if (is.na(sum(x[nrow(x),]))) {
    value <- NA_real_
  } else {
    value <- fun(x[, 1], x[, 2], ...)
  }
  # Check value
  checkmate::assert_number(value, na.ok = TRUE)
  # Return value
  return(value)
}

#' Mean Percent Error
#'
#' @param x [numeric()] [vector()] actual values
#' @param y [numeric()] [vector()] forecast values
#'
#' @return [numeric()] mean percent error
#' @export
#'
#' @examples
#' mpe(1:10, 10:1)
#'
mpe <- function (x, y) {
  # Check arguments
  checkmate::assert_numeric(x, finite = TRUE, min.len = 1)
  checkmate::assert_numeric(y, finite = TRUE, min.len = 1)
  checkmate::assert_vector(x, strict = TRUE)
  checkmate::assert_vector(y, strict = TRUE)
  checkmate::assert_true(length(x) == length(y))
  # Return value
  return(mean(pct_error(x, y), na.rm = TRUE))
}

#' Mean Raw Error
#'
#' @param x [numeric()][vector()] observed values
#' @param y [numeric()][vector()] forecast values
#'
#' @return [numeric()] scalar
#' @export
#'
mre <- function (x, y) {
  # Check arguments
  checkmate::assert_numeric(x, finite = TRUE, min.len = 1)
  checkmate::assert_numeric(y, finite = TRUE, min.len = 1)
  checkmate::assert_vector(x, strict = TRUE)
  checkmate::assert_vector(y, strict = TRUE)
  checkmate::assert_true(length(x) == length(y))
  # Return value
  return(mean(raw_error(x, y), na.rm = TRUE))
}

#' Percent Error
#'
#' @param x [numeric()][vector()] observed values
#' @param y [numeric()][vector()] forecast values
#'
#' @return [numeric()][vector()] percent error
#' @export
#'
#' @examples
#' pct_error(0:10, 10:0)
#'
pct_error <- function (x, y) {
  # Check arguments
  checkmate::assert_numeric(x, finite = TRUE, min.len = 1)
  checkmate::assert_numeric(y, finite = TRUE, min.len = 1)
  checkmate::assert_vector(x, strict = TRUE)
  checkmate::assert_vector(y, strict = TRUE)
  checkmate::assert_true(length(x) == length(y))
  # Return value
  return((y - x) / x)
}

#' Raw Error
#'
#' @param x [numeric()][vector()] observed values
#' @param y [numeric()][vector()] forecast values
#'
#' @return [numeric()][vector()] raw error
#' @export
#'
raw_error <- function (x, y) {
  # Check arguments
  checkmate::assert_numeric(x, finite = TRUE, min.len = 1)
  checkmate::assert_numeric(y, finite = TRUE, min.len = 1)
  checkmate::assert_vector(x, strict = TRUE)
  checkmate::assert_vector(y, strict = TRUE)
  checkmate::assert_true(length(x) == length(y))
  # Return value
  return(y - x)
}

#' Root Mean Squared Error
#'
#' @param x [numeric()][vector()] observed values
#' @param y [numeric()][vector()] forecast values
#'
#' @return [numeric()] root mean squared error
#' @export
#'
#' @examples
#' rmse(1:10, 10:1)
#'
rmse <- function (x, y) {
  # Check arguments
  checkmate::assert_numeric(x, finite = TRUE, min.len = 1)
  checkmate::assert_numeric(y, finite = TRUE, min.len = 1)
  checkmate::assert_vector(x, strict = TRUE)
  checkmate::assert_vector(y, strict = TRUE)
  checkmate::assert_true(length(x) == length(y))
  # Return value
  return(sqrt(mean((y - x)^2, na.rm = TRUE)))
}

#' Scaled Mean Arctangent Absolute Percentage Error
#'
#' @param x [numeric()] [vector()] actual values
#' @param y [numeric()] [vector()] forecast values
#'
#' @return [numeric()] scaled mean arctangent absolute percentage error
#' @export
#'
#' @references Kim and Kim 2016 doi.org/10.1016/j.ijforecast.2015.12.003
#'
#' @examples
#' x <- rnorm(100, 0, 1)
#' y1 <- x
#' smaape(x, y1)
#'
#' y2 <- rep(mean(x), length(x))
#' smaape(x, y2)
#'
#' smaape(c(x, NA), c(y1, NA))
#'
smaape <- function (x, y) {
  # Check arguments
  checkmate::assert_numeric(x, finite = TRUE, min.len = 1)
  checkmate::assert_numeric(y, finite = TRUE, min.len = 1)
  checkmate::assert_vector(x, strict = TRUE)
  checkmate::assert_vector(y, strict = TRUE)
  checkmate::assert_true(length(x) == length(y))
  # Return value
  return(maape(x, y) / maape(x, rep(mean(x, na.rm = TRUE), length(x))))
}

#' Standardised Root Mean Squared Error
#'
#' @param x [numeric()][vector()] observed values
#' @param y [numeric()][vector()] forecast values
#'
#' @return [numeric()] root mean squared error
#' @export
#'
#' @examples
#' x <- rnorm(100, 0, 1)
#' y1 <- x
#' srmse(x, y1)
#'
#' y2 <- rep(mean(x), length(x))
#' srmse(x, y2)
#'
srmse <- function (x, y) {
  # Check arguments
  checkmate::assert_numeric(x, finite = TRUE, min.len = 1)
  checkmate::assert_numeric(y, finite = TRUE, min.len = 1)
  checkmate::assert_vector(x, strict = TRUE)
  checkmate::assert_vector(y, strict = TRUE)
  checkmate::assert_true(length(x) == length(y))
  # Return value
  return(rmse(x, y) / stats::sd(x))
}
