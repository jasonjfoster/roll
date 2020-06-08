##' Rolling Any
##'
##' A function for computing the rolling and expanding any of time-series data.
##'
##' @param x logical vector or matrix. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @return An object of the same class and dimension as \code{x} with the rolling and expanding
##' any.
##' @examples
##' n <- 15
##' x <- rnorm(n)
##' 
##' # rolling any with complete windows
##' roll_any(x < 0, width = 5)
##' 
##' # rolling any with partial windows
##' roll_any(x < 0, width = 5)
##' 
##' # expanding any with partial windows
##' roll_any(x < 0, width = n)
##' @export
roll_any <- function(x, width, min_obs = width,
                     complete_obs = FALSE, na_restore = FALSE,
                     online = TRUE) {
  return(.Call(`_roll_roll_any`,
               x,
               as.integer(width),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.logical(online)
  ))
}

##' Rolling All
##'
##' A function for computing the rolling and expanding all of time-series data.
##'
##' @param x logical vector or matrix. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @return An object of the same class and dimension as \code{x} with the rolling and expanding
##' all.
##' @examples
##' n <- 15
##' x <- rnorm(n)
##' 
##' # rolling all with complete windows
##' roll_all(x < 0, width = 5)
##' 
##' # rolling all with partial windows
##' roll_all(x < 0, width = 5)
##' 
##' # expanding all with partial windows
##' roll_all(x < 0, width = n)
##' @export
roll_all <- function(x, width, min_obs = width,
                     complete_obs = FALSE, na_restore = FALSE,
                     online = TRUE) {
  return(.Call(`_roll_roll_all`,
               x,
               as.integer(width),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.logical(online)
  ))
}

##' Rolling Sums
##'
##' A function for computing the rolling and expanding sums of time-series data.
##'
##' @param x vector or matrix. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @return An object of the same class and dimension as \code{x} with the rolling and expanding
##' sums.
##' @examples
##' n <- 15
##' x <- rnorm(n)
##' weights <- 0.9 ^ (n:1)
##' 
##' # rolling sums with complete windows
##' roll_sum(x, width = 5)
##' 
##' # rolling sums with partial windows
##' roll_sum(x, width = 5, min_obs = 1)
##' 
##' # expanding sums with partial windows
##' roll_sum(x, width = n, min_obs = 1)
##' 
##' # expanding sums with partial windows and weights
##' roll_sum(x, width = n, min_obs = 1, weights = weights)
##' @export
roll_sum <- function(x, width, weights = rep(1, width),
                     min_obs = width, complete_obs = FALSE, na_restore = FALSE,
                     online = TRUE) {
  return(.Call(`_roll_roll_sum`,
               x,
               as.integer(width),
               as.numeric(weights),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.logical(online)
  ))
}

##' Rolling Products
##'
##' A function for computing the rolling and expanding products of time-series data.
##'
##' @param x vector or matrix. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @return An object of the same class and dimension as \code{x} with the rolling and expanding
##' products.
##' @examples
##' n <- 15
##' x <- rnorm(n)
##' weights <- 0.9 ^ (n:1)
##' 
##' # rolling products with complete windows
##' roll_prod(x, width = 5)
##' 
##' # rolling products with partial windows
##' roll_prod(x, width = 5, min_obs = 1)
##' 
##' # expanding products with partial windows
##' roll_prod(x, width = n, min_obs = 1)
##' 
##' # expanding products with partial windows and weights
##' roll_prod(x, width = n, min_obs = 1, weights = weights)
##' @export
roll_prod <- function(x, width, weights = rep(1, width),
                      min_obs = width, complete_obs = FALSE, na_restore = FALSE,
                      online = TRUE) {
  return(.Call(`_roll_roll_prod`,
               x,
               as.integer(width),
               as.numeric(weights),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.logical(online)
  ))
}

##' Rolling Means
##'
##' A function for computing the rolling and expanding means of time-series data.
##'
##' @param x vector or matrix. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @return An object of the same class and dimension as \code{x} with the rolling and expanding
##' means.
##' @examples
##' n <- 15
##' x <- rnorm(n)
##' weights <- 0.9 ^ (n:1)
##' 
##' # rolling means with complete windows
##' roll_mean(x, width = 5)
##' 
##' # rolling means with partial windows
##' roll_mean(x, width = 5, min_obs = 1)
##' 
##' # expanding means with partial windows
##' roll_mean(x, width = n, min_obs = 1)
##' 
##' # expanding means with partial windows and weights
##' roll_mean(x, width = n, min_obs = 1, weights = weights)
##' @export
roll_mean <- function(x, width, weights = rep(1, width),
                      min_obs = width, complete_obs = FALSE, na_restore = FALSE,
                      online = TRUE) {
  return(.Call(`_roll_roll_mean`,
               x,
               as.integer(width),
               as.numeric(weights),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.logical(online)
  ))
}

##' Rolling Minimums
##'
##' A function for computing the rolling and expanding minimums of time-series data.
##'
##' @param x vector or matrix. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @return An object of the same class and dimension as \code{x} with the rolling and expanding
##' minimums.
##' @examples
##' n <- 15
##' x <- rnorm(n)
##' weights <- 0.9 ^ (n:1)
##' 
##' # rolling minimums with complete windows
##' roll_min(x, width = 5)
##' 
##' # rolling minimums with partial windows
##' roll_min(x, width = 5, min_obs = 1)
##' 
##' # expanding minimums with partial windows
##' roll_min(x, width = n, min_obs = 1)
##' 
##' # expanding minimums with partial windows and weights
##' roll_min(x, width = n, min_obs = 1, weights = weights)
##' @export
roll_min <- function(x, width, weights = rep(1, width),
                     min_obs = width, complete_obs = FALSE, na_restore = FALSE,
                     online = TRUE) {
  return(.Call(`_roll_roll_quantile`,
               x,
               as.integer(width),
               as.numeric(weights),
               as.numeric(0),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.logical(online)
  ))
}

##' Rolling Maximums
##'
##' A function for computing the rolling and expanding maximums of time-series data.
##'
##' @param x vector or matrix. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @return An object of the same class and dimension as \code{x} with the rolling and expanding
##' maximums.
##' @examples
##' n <- 15
##' x <- rnorm(n)
##' weights <- 0.9 ^ (n:1)
##' 
##' # rolling maximums with complete windows
##' roll_max(x, width = 5)
##' 
##' # rolling maximums with partial windows
##' roll_max(x, width = 5, min_obs = 1)
##' 
##' # expanding maximums with partial windows
##' roll_max(x, width = n, min_obs = 1)
##' 
##' # expanding maximums with partial windows and weights
##' roll_max(x, width = n, min_obs = 1, weights = weights)
##' @export
roll_max <- function(x, width, weights = rep(1, width),
                     min_obs = width, complete_obs = FALSE, na_restore = FALSE,
                     online = TRUE) {
  return(.Call(`_roll_roll_quantile`,
               x,
               as.integer(width),
               as.numeric(weights),
               as.numeric(1),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.logical(online)
  ))
}

##' Rolling Index of Minimums
##'
##' A function for computing the rolling and expanding index of minimums of time-series data.
##'
##' @param x vector or matrix. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @return An object of the same class and dimension as \code{x} with the rolling and expanding
##' index of minimums.
##' @examples
##' n <- 15
##' x <- rnorm(n)
##' weights <- 0.9 ^ (n:1)
##' 
##' # rolling index of minimums with complete windows
##' roll_idxmin(x, width = 5)
##' 
##' # rolling index of minimums with partial windows
##' roll_idxmin(x, width = 5, min_obs = 1)
##' 
##' # expanding index of minimums with partial windows
##' roll_idxmin(x, width = n, min_obs = 1)
##' 
##' # expanding index of minimums with partial windows and weights
##' roll_idxmin(x, width = n, min_obs = 1, weights = weights)
##' @export
roll_idxmin <- function(x, width, weights = rep(1, width),
                        min_obs = width, complete_obs = FALSE, na_restore = FALSE,
                        online = TRUE) {
  return(.Call(`_roll_roll_idxquantile`,
               x,
               as.integer(width),
               as.numeric(weights),
               as.numeric(0),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.logical(online)
  ))
}

##' Rolling Index of Maximums
##'
##' A function for computing the rolling and expanding index of maximums of time-series data.
##'
##' @param x vector or matrix. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @return An object of the same class and dimension as \code{x} with the rolling and expanding
##' index of maximums.
##' @examples
##' n <- 15
##' x <- rnorm(n)
##' weights <- 0.9 ^ (n:1)
##' 
##' # rolling index of maximums with complete windows
##' roll_idxmax(x, width = 5)
##' 
##' # rolling index of maximums with partial windows
##' roll_idxmax(x, width = 5, min_obs = 1)
##' 
##' # expanding index of maximums with partial windows
##' roll_idxmax(x, width = n, min_obs = 1)
##' 
##' # expanding index of maximums with partial windows and weights
##' roll_idxmax(x, width = n, min_obs = 1, weights = weights)
##' @export
roll_idxmax <- function(x, width, weights = rep(1, width),
                        min_obs = width, complete_obs = FALSE, na_restore = FALSE,
                        online = TRUE) {
  return(.Call(`_roll_roll_idxquantile`,
               x,
               as.integer(width),
               as.numeric(weights),
               as.numeric(1),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.logical(online)
  ))
}

##' Rolling Medians
##'
##' A function for computing the rolling and expanding medians of time-series data.
##'
##' @param x vector or matrix. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @return An object of the same class and dimension as \code{x} with the rolling and expanding
##' medians.
##' @examples
##' n <- 15
##' x <- rnorm(n)
##' weights <- 0.9 ^ (n:1)
##' 
##' # rolling medians with complete windows
##' roll_median(x, width = 5)
##' 
##' # rolling medians with partial windows
##' roll_median(x, width = 5, min_obs = 1)
##' 
##' # expanding medians with partial windows
##' roll_median(x, width = n, min_obs = 1)
##' 
##' # expanding medians with partial windows and weights
##' roll_median(x, width = n, min_obs = 1, weights = weights)
##' @export
roll_median <- function(x, width, weights = rep(1, width),
                        min_obs = width, complete_obs = FALSE, na_restore = FALSE,
                        online = FALSE) {
  return(.Call(`_roll_roll_quantile`,
               x,
               as.integer(width),
               as.numeric(weights),
               as.numeric(0.5),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.logical(online)
  ))
}

##' Rolling Quantiles
##'
##' A function for computing the rolling and expanding quantiles of time-series data.
##'
##' @param x vector or matrix. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param p numeric. Probability between zero and one.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @details The methodology for computing the quantiles is based on the inverse of the empirical
##' distribution function with averaging at discontinuities (see "Definition 2" in Hyndman and Fan, 1996). 
##' @return An object of the same class and dimension as \code{x} with the rolling and expanding
##' quantiles.
##' @references Hyndman, R.J. and Fan, Y. (1996). "Sample quantiles in statistical packages."
##' \emph{American Statistician}, 50(4), 361-365.
##' @examples
##' n <- 15
##' x <- rnorm(n)
##' weights <- 0.9 ^ (n:1)
##' 
##' # rolling quantiles with complete windows
##' roll_quantile(x, width = 5)
##' 
##' # rolling quantiles with partial windows
##' roll_quantile(x, width = 5, min_obs = 1)
##' 
##' # expanding quantiles with partial windows
##' roll_quantile(x, width = n, min_obs = 1)
##' 
##' # expanding quantiles with partial windows and weights
##' roll_quantile(x, width = n, min_obs = 1, weights = weights)
##' @export
roll_quantile <- function(x, width, weights = rep(1, width), p = 0.5,
                          min_obs = width, complete_obs = FALSE, na_restore = FALSE,
                          online = FALSE) {
  return(.Call(`_roll_roll_quantile`,
               x,
               as.integer(width),
               as.numeric(weights),
               as.numeric(p),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.logical(online)
  ))
}

##' Rolling Variances
##'
##' A function for computing the rolling and expanding variances of time-series data.
##'
##' @param x vector or matrix. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param center logical. If \code{TRUE} then the weighted mean of each variable is used,
##' if \code{FALSE} then zero is used.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @details The denominator used gives an unbiased estimate of the variance,
##' so if the weights are the default then the divisor \code{n - 1} is obtained.
##' @return An object of the same class and dimension as \code{x} with the rolling and expanding
##' variances.
##' @examples
##' n <- 15
##' x <- rnorm(n)
##' weights <- 0.9 ^ (n:1)
##' 
##' # rolling variances with complete windows
##' roll_var(x, width = 5)
##' 
##' # rolling variances with partial windows
##' roll_var(x, width = 5, min_obs = 1)
##' 
##' # expanding variances with partial windows
##' roll_var(x, width = n, min_obs = 1)
##' 
##' # expanding variances with partial windows and weights
##' roll_var(x, width = n, min_obs = 1, weights = weights)
##' @export
roll_var <- function(x, width, weights = rep(1, width), center = TRUE,
                     min_obs = width, complete_obs = FALSE, na_restore = FALSE,
                     online = TRUE) {
  return(.Call(`_roll_roll_var`,
               x,
               as.integer(width),
               as.numeric(weights),
               as.logical(center),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.logical(online)
  ))
}

##' Rolling Standard Deviations
##'
##' A function for computing the rolling and expanding standard deviations of time-series data.
##'
##' @param x vector or matrix. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param center logical. If \code{TRUE} then the weighted mean of each variable is used,
##' if \code{FALSE} then zero is used.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @details The denominator used gives an unbiased estimate of the standard deviation,
##' so if the weights are the default then the divisor \code{n - 1} is obtained.
##' @return An object of the same class and dimension as \code{x} with the rolling and expanding
##' standard deviations.
##' @examples
##' n <- 15
##' x <- rnorm(n)
##' weights <- 0.9 ^ (n:1)
##' 
##' # rolling standard deviations with complete windows
##' roll_sd(x, width = 5)
##' 
##' # rolling standard deviations with partial windows
##' roll_sd(x, width = 5, min_obs = 1)
##' 
##' # expanding standard deviations with partial windows
##' roll_sd(x, width = n, min_obs = 1)
##' 
##' # expanding standard deviations with partial windows and weights
##' roll_sd(x, width = n, min_obs = 1, weights = weights)
##' @export
roll_sd <- function(x, width, weights = rep(1, width), center = TRUE,
                    min_obs = width, complete_obs = FALSE, na_restore = FALSE,
                    online = TRUE) {
  return(.Call(`_roll_roll_sd`,
               x,
               as.integer(width),
               as.numeric(weights),
               as.logical(center),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.logical(online)
  ))
}

##' Rolling Scaling and Centering
##'
##' A function for computing the rolling and expanding scaling and centering of time-series data.
##'
##' @param x vector or matrix. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param center logical. If \code{TRUE} then the weighted mean of each variable is used,
##' if \code{FALSE} then zero is used.
##' @param scale logical. If \code{TRUE} then the weighted standard deviation of each variable is used,
##' if \code{FALSE} then no scaling is done.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @details If \code{center} is \code{TRUE} then centering is done by subtracting the weighted mean from 
##' each variable, if \code{FALSE} then zero is used. After centering, if \code{scale} is \code{TRUE} then 
##' scaling is done by dividing by the weighted standard deviation for each variable if \code{center} is 
##' \code{TRUE}, and the root mean square otherwise. If \code{scale} is \code{FALSE} then no scaling is 
##' done.
##' 
##' The denominator used gives an unbiased estimate of the standard deviation,
##' so if the weights are the default then the divisor \code{n - 1} is obtained.
##' @return An object of the same class and dimension as \code{x} with the rolling and expanding
##' scaling and centering.
##' @examples
##' n <- 15
##' x <- rnorm(n)
##' weights <- 0.9 ^ (n:1)
##' 
##' # rolling z-scores with complete windows
##' roll_scale(x, width = 5)
##' 
##' # rolling z-scores with partial windows
##' roll_scale(x, width = 5, min_obs = 1)
##' 
##' # expanding z-scores with partial windows
##' roll_scale(x, width = n, min_obs = 1)
##' 
##' # expanding z-scores with partial windows and weights
##' roll_scale(x, width = n, min_obs = 1, weights = weights)
##' @export
roll_scale <- function(x, width, weights = rep(1, width), center = TRUE, scale = TRUE,
                       min_obs = width, complete_obs = FALSE, na_restore = FALSE,
                       online = TRUE) {
  return(.Call(`_roll_roll_scale`,
               x,
               as.integer(width),
               as.numeric(weights),
               as.logical(center),
               as.logical(scale),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.logical(online)
  ))
}

##' Rolling Covariances
##'
##' A function for computing the rolling and expanding covariances of time-series data.
##' 
##' @param x vector or matrix. Rows are observations and columns are variables.
##' @param y vector or matrix. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param center logical. If \code{TRUE} then the weighted mean of each variable is used,
##' if \code{FALSE} then zero is used.
##' @param scale logical. If \code{TRUE} then the weighted standard deviation of each variable is used,
##' if \code{FALSE} then no scaling is done.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then pairwise is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @details The denominator used gives an unbiased estimate of the covariance,
##' so if the weights are the default then the divisor \code{n - 1} is obtained.
##' @return A cube with each slice the rolling and expanding covariances.
##' @examples
##' n <- 15
##' x <- rnorm(n)
##' y <- rnorm(n)
##' weights <- 0.9 ^ (n:1)
##' 
##' # rolling covariances with complete windows
##' roll_cov(x, y, width = 5)
##' 
##' # rolling covariances with partial windows
##' roll_cov(x, y, width = 5, min_obs = 1)
##' 
##' # expanding covariances with partial windows
##' roll_cov(x, y, width = n, min_obs = 1)
##' 
##' # expanding covariances with partial windows and weights
##' roll_cov(x, y, width = n, min_obs = 1, weights = weights)
##' @export
roll_cov <- function(x, y = NULL, width, weights = rep(1, width), center = TRUE, scale = FALSE,
                     min_obs = width, complete_obs = TRUE, na_restore = FALSE,
                     online = TRUE) {
  return(.Call(`_roll_roll_cov`,
               x, y,
               as.integer(width),
               as.numeric(weights),
               as.logical(center),
               as.logical(scale),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.logical(online)
  ))
}

##' Rolling Correlations
##'
##' A function for computing the rolling and expanding correlations of time-series data.
##' 
##' @param x vector or matrix. Rows are observations and columns are variables.
##' @param y vector or matrix. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param center logical. If \code{TRUE} then the weighted mean of each variable is used,
##' if \code{FALSE} then zero is used.
##' @param scale logical. If \code{TRUE} then the weighted standard deviation of each variable is used,
##' if \code{FALSE} then no scaling is done.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then pairwise is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @details The denominator used gives an unbiased estimate of the covariance,  
##' so if the weights are the default then the divisor \code{n - 1} is obtained.
##' @return A cube with each slice the rolling and expanding correlations.
##' @examples
##' n <- 15
##' x <- rnorm(n)
##' y <- rnorm(n)
##' weights <- 0.9 ^ (n:1)
##' 
##' # rolling correlations with complete windows
##' roll_cor(x, y, width = 5)
##' 
##' # rolling correlations with partial windows
##' roll_cor(x, y, width = 5, min_obs = 1)
##' 
##' # expanding correlations with partial windows
##' roll_cor(x, y, width = n, min_obs = 1)
##' 
##' # expanding correlations with partial windows and weights
##' roll_cor(x, y, width = n, min_obs = 1, weights = weights)
##' @export
roll_cor <- function(x, y = NULL, width, weights = rep(1, width), center = TRUE, scale = TRUE,
                     min_obs = width, complete_obs = TRUE, na_restore = FALSE,
                     online = TRUE) {
  return(.Call(`_roll_roll_cov`,
               x, y,
               as.integer(width),
               as.numeric(weights),
               as.logical(center),
               as.logical(scale),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.logical(online)
  ))
}

##' Rolling Linear Models
##'
##' A function for computing the rolling and expanding linear models of time-series data.
##' 
##' @param x vector or matrix. Rows are observations and columns are the independent variables.
##' @param y vector or matrix. Rows are observations and columns are the dependent variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param intercept logical. Either \code{TRUE} to include or \code{FALSE} to remove the intercept.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then pairwise is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @return A list containing the following components:
##' \item{coefficients}{A list of objects with the rolling and expanding coefficients for each \code{y}.
##' An object is the same class and dimension (with an added column for the intercept) as \code{x}.}
##' \item{r.squared}{A list of objects with the rolling and expanding r-squareds for each \code{y}.
##' An object is the same class as \code{x}.}
##' \item{std.error}{A list of objects with the rolling and expanding standard errors for each \code{y}.
##' An object is the same class and dimension (with an added column for the intercept) as \code{x}.}
##' @examples
##' n <- 15
##' x <- rnorm(n)
##' y <- rnorm(n)
##' weights <- 0.9 ^ (n:1)
##' 
##' # rolling regressions with complete windows
##' roll_lm(x, y, width = 5)
##' 
##' # rolling regressions with partial windows
##' roll_lm(x, y, width = 5, min_obs = 1)
##' 
##' # expanding regressions with partial windows
##' roll_lm(x, y, width = n, min_obs = 1)
##' 
##' # expanding regressions with partial windows and weights
##' roll_lm(x, y, width = n, min_obs = 1, weights = weights)
##' @export
roll_lm <- function(x, y, width, weights = rep(1, width), intercept = TRUE,
                    min_obs = width, complete_obs = TRUE, na_restore = FALSE,
                    online = TRUE) {
  return(.Call(`_roll_roll_lm`,
               x, y,
               as.integer(width),
               as.numeric(weights),
               as.logical(intercept),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.logical(online)
  ))
}