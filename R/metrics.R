#' Arctangent Absolute Percentage Error
#'
#' @param a [numeric()] [vector()] of actual values
#' @param f [numeric()] [vector()] of forecast values
#'
#' @return [numeric()] [vector()] of arctangent absolute percentage errors
#' @export
#'
#' @references Kim and Kim 2016 doi.org/10.1016/j.ijforecast.2015.12.003
#'
#' @examples
#' a <- seq(0, 10, 0.1)
#' f <- seq(10, 0, -0.1)
#' m <- aape(a, f)
#' plot(m, type = "l")
#'
aape <- function (a, f) {
  # Check arguments

  # Return arctangent absolute percentage error
  return(atan(ape(a, f)))
}

#' Absolute Percentage Error
#'
#' @param a [numeric()] [vector()] of actual values
#' @param f [numeric()] [vector()] of forecast values
#'
#' @return [numeric()] [vector()] of absolute percentage errors
#' @export
#'
#' @examples
#' a <- seq(0, 10, 0.1)
#' f <- seq(10, 0, -0.1)
#' m <- ape(a, f)
#' plot(m, type = "l")
#'
ape <- function (a, f) {
  # Check arguments

  # Return absolute percentage error
  return(abs((a - f) / a))
}

#' Mean Arctangent Absolute Percentage Error
#'
#' @param a [numeric()] [vector()] of actual values
#' @param f [numeric()] [vector()] of forecast values
#'
#' @return [numeric()] mean arctangent absolute percentage error
#' @export
#'
#' @references Kim and Kim 2016 doi.org/10.1016/j.ijforecast.2015.12.003
#'
#' @examples
#' maape(1:10, 10:1)
#'
maape <- function (a, f) {
  # Check arguments

  # Return maape
  return(mean(aape(a, f), na.rm = TRUE))
}

#' Mean Absolute Error
#'
#' @param a [numeric()] [vector()] of actual values
#' @param f [numeric()] [vector()] of forecast values
#'
#' @return [numeric()] mean absolute error
#' @export
#'
#' @examples
#' mae(1:10, 10:1)
#'
mae <- function (a, f) {
  # Check arguments

  # Return mae
  return(mean(abs(a - f), na.rm = TRUE))
}

#' Mean Absolute Percentage Error
#'
#' @param a [numeric()] [vector()] of actual values
#' @param f [numeric()] [vector()] of forecast values
#'
#' @return [numeric()] mean absolute percentage error
#' @export
#'
#' @examples
#' mape(1:10, 10:1)
#'
mape <- function (a, f) {
  # Check arguments

  # Return mape
  return(mean(ape(a, f), na.rm = TRUE))
}

#' Arithmetic Mean Or NA
#'
#' @param x [numeric()] [vector()]
#' @param max_na [integer()] scalar
#'
#' @return [numeric()] scalar or \code{NA_real_}
#' @export
#'
#' @examples
#' # Ends in NA
#' x <- c(1:9, NA)
#' mean_or_na(x)
#'
#' # Too many NA values
#' x <- c(1:3, NA, NA, NA, 7:10)
#' mean_or_na(x, 2)
#'
#' # Mean
#' x <- 1:10
#' mean_or_na(x)
#'
mean_or_na <- function (x, max_na = Inf) {
  if (is.na(x[length(x)])) {
    val <- NA_real_
  } else if (length(which(is.na(x))) > max_na) {
    val <- NA_real_
  } else {
    val <- base::mean(x = x, na.rm = TRUE)
  }
  return(val)
}

#' Mean Percentage Error
#'
#' @param a [numeric()] [vector()] of actual values
#' @param f [numeric()] [vector()] of forecast values
#'
#' @return [numeric()] mean percentage error
#' @export
#'
#' @examples
#' mpe(1:10, 10:1)
#'
mpe <- function (a, f) {
  # Check arguments

  # Return mpe
  return(mean((f - a) / a, na.rm = TRUE))
}

#' Mean Raw Error
#'
#' @param a [numeric()] [vector()] of actual values
#' @param f [numeric()] [vector()] of forecast values
#'
#' @return [numeric()] mean raw error
#' @export
#'
#' @examples
#' mre(1:10, 10:1)
#'
mre <- function (a, f) {
  # Check arguments

  # Return mre
  return(mean(f - a, na.rm = TRUE))
}

#' Percentage Error
#'
#' @param a [numeric()] [vector()] of actual values
#' @param f [numeric()] [vector()] of forecast values
#'
#' @return [numeric()] symmetric percentage error
#' @export
#'
#' @examples
#' pe(0:10, 10:0)
#'
pe <- function (a, f) {
  # Check arguments

  # Return percent error
  return((f - a) / a)
}

#' Raw Error
#'
#' @param a [numeric()] [vector()] of actual values
#' @param f [numeric()] [vector()] of forecast values
#'
#' @return [numeric()] raw error
#' @export
#'
#' @examples
#' re(1:10, 10:1)
#'
re <- function (a, f) {
  # Check arguments

  # Return raw error
  return(f - a)
}

#' Pearson Correlation Or NA
#'
#' @param x [numeric()] [matrix()] with \code{ncol = 2}
#' @param max_na [integer()] scalar
#'
#' @return [numeric()] scalar or \code{NA_real_}
#' @export
#'
#' @examples
#' # Ends in NA
#' a <- c(1:9, NA)
#' b <- c(1:10)
#' x <- matrix(c(a, b), ncol = 2)
#' rho_or_na(x)
#'
#' # Too many NA values
#' a <- c(1:3, NA, NA, NA, 7:10)
#' b <- c(1:10)
#' x <- matrix(c(a, b), ncol = 2)
#' rho_or_na(x, max_na = 2)
#'
#' # Correlation coefficient rho
#' a <- 1:10
#' b <- c(1:10)
#' x <- matrix(c(a, b), ncol = 2)
#' rho_or_na(x)
#'
rho_or_na <- function (x, max_na = Inf) {
  # Check arguments
  checkmate::assert_matrix(x, mode = "numeric", min.rows = 1, ncols = 2)
  checkmate::assert_number(max_na, lower = 0)
  # Compute value
  if (is.na(sum(x[nrow(x), ]))) {
    val <- NA_real_
  } else if (length(which(is.na(rowSums(x)))) > max_na) {
    val <- NA_real_
  } else if (length(which(!is.na(rowSums(x)))) < 3) {
    val <- NA_real_
  } else {
    val <- stats::cor(x = x[, 1], y = x[, 2], use = "pairwise.complete.obs")
  }
  return(val)
}

#' Root Mean Square Or NA
#'
#' @param x [numeric()] [vector()]
#' @param max_na [integer()] scalar
#'
#' @return [numeric()] scalar or \code{NA_real_}
#' @export
#'
#' @examples
#' # Ends in NA
#' x <- c(1:9, NA)
#' rms_or_na(x)
#'
#' # Too many NA values
#' x <- c(1:3, NA, NA, NA, 7:10)
#' rms_or_na(x, 2)
#'
#' # Root mean squared
#' x <- 1:10
#' rms_or_na(x)
#'
rms_or_na <- function (x, max_na = Inf) {
  if (is.na(x[length(x)])) {
    val <- NA_real_
  } else if (length(which(is.na(x))) > max_na) {
    val <- NA_real_
  } else {
    val <- base::sqrt(mean((x^2), na.rm = TRUE))
  }
  return(val)
}

#' Root Mean Squared Error
#'
#' @param a [numeric()] [vector()] of actual values
#' @param f [numeric()] [vector()] of forecast values
#'
#' @return [numeric()] root mean squared error
#' @export
#'
#' @examples
#' rmse(1:10, 10:1)
#'
rmse <- function (a, f) {
  # Check arguments

  # return rmse
  return(sqrt(mean((a - f)^2, na.rm = TRUE)))
}
