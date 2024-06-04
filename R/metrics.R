#' Running Or Total Mean Raw Error
#'
#' @param x [matrix()] with two columns or [numeric()][vector()]
#' @param y [numeric()][vector()]
#' @param running [logical()]
#'
#' @return [numeric()]
#' @export
#'
#' @examples
#'
#' # Without NA
#' mre(1:10, 0:9, FALSE)
#' mre(1:10, 0:9, TRUE)
#' mre(matrix(c(1:10, 0:9), ncol = 2), running = FALSE)
#' mre(matrix(c(1:10, 0:9), ncol = 2), running = TRUE)
#'
#' # With NA
#' mre(matrix(c(1:10, NA, 1:9), ncol = 2), running = FALSE)
#' mre(matrix(c(1:10, NA, 1:9), ncol = 2), running = TRUE)
#'
mre <- function (x, y = NULL, running = FALSE) {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_numeric(y, finite = TRUE, null.ok = TRUE)
  checkmate::assert_logical(running, len = 1L, any.missing = FALSE)

  # Vectors to matrix ----------------------------------------------------------

  if (!is.matrix(x)) {
    if (length(x) == length(y)) {
      x <- matrix(c(x, y), ncol = 2L)
    } else {
      stop("rmse: length(x) must equal length(y")
    }
  }

  # Running or total -----------------------------------------------------------

  if (running) {
    value <- runner::runner(x, f = mre, y = NULL, running = FALSE)
  } else {
    value <- mean((x[, 2] - x[, 1]), na.rm = TRUE)
  }

  # Replace NaN by NA ----------------------------------------------------------

  value <- dplyr::na_if(value, NaN)

  # Return value ---------------------------------------------------------------

  return(value)
}

#' Running Or Total Square Root Of The Mean Squared Error
#'
#' @param x [matrix()] with two columns or [numeric()][vector()]
#' @param y [numeric()][vector()]
#' @param running [logical()]
#'
#' @return [numeric()]
#' @export
#'
#' @examples
#'
#' # Without NA
#' rmse(1:10, 0:9, FALSE)
#' rmse(1:10, 0:9, TRUE)
#' rmse(matrix(c(1:10, 0:9), ncol = 2), running = FALSE)
#' rmse(matrix(c(1:10, 0:9), ncol = 2), running = TRUE)
#'
#' # With NA
#' rmse(matrix(c(1:10, NA, 1:9), ncol = 2), running = FALSE)
#' rmse(matrix(c(1:10, NA, 1:9), ncol = 2), running = TRUE)
#'
rmse <- function (x, y = NULL, running = FALSE) {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_numeric(y, finite = TRUE, null.ok = TRUE)
  checkmate::assert_logical(running, len = 1L, any.missing = FALSE)

  # Vectors to matrix ----------------------------------------------------------

  if (!is.matrix(x)) {
    if (length(x) == length(y)) {
      x <- matrix(c(x, y), ncol = 2L)
    } else {
      stop("rmse: length(x) must equal length(y")
    }
  }

  # Running or total -----------------------------------------------------------

  if (running) {
    value <- runner::runner(x, f = rmse, y = NULL, running = FALSE)
  } else {
    value <- sqrt(mean((x[, 2] - x[, 1])^2, na.rm = TRUE))
  }

  # Replace NaN by NA ----------------------------------------------------------

  value <- dplyr::na_if(value, NaN)

  # Return value ---------------------------------------------------------------

  return(value)
}

#' Scaled Square Root of the Mean Squared Error (SRMSE)
#'
#' @param x [numeric()] [vector()] Observed values
#' @param y [numeric()] [vector()] Forecast values
#' @param s [numeric()] [vector()] Forecast values used to scale
#' @param running [logical()]
#'
#' @return [numeric()]
#' @export
#'
#' @examples
#' # Perfect
#' srmse(0:10, 0:10, 0:10) # NaN because 0/0
#'
#' # Correlated
#' srmse(0:10, 1:11, 1:11) # 1 because y == s
#'
#' # Running
#' srmse(0:10, 1:11, 1:11, TRUE) # rep(1, 11) because y == s
#'
srmse <- function (x, y = NULL, s, running = FALSE) {
  rmse(x, y, running) / rmse(x, s, running)
}
