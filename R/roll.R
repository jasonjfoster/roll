##' Rolling Any
##'
##' A function for computing rolling any of time-series data.
##'
##' @param x logical matrix or xts object. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @return An object of the same class and dimension as \code{x} with the rolling any.
##' @examples
##' n_vars <- 3
##' n_obs <- 15
##' x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##' 
##' # rolling any
##' result <- roll_any(x < 0, 5)
##' 
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
##' A function for computing rolling all of time-series data.
##'
##' @param x logical matrix or xts object. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @return An object of the same class and dimension as \code{x} with the rolling all.
##' @examples
##' n_vars <- 3
##' n_obs <- 15
##' x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##' 
##' # rolling all
##' result <- roll_all(x < 0, 5)
##' 
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
##' A function for computing rolling sums of time-series data.
##'
##' @param x matrix or xts object. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @return An object of the same class and dimension as \code{x} with the rolling sums.
##' @examples
##' n_vars <- 3
##' n_obs <- 15
##' x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##' 
##' # rolling sums
##' result <- roll_sum(x, 5)
##' 
##' # rolling sums with exponential decay
##' weights <- 0.9 ^ (5:1)
##' result <- roll_sum(x, 5, weights)
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
##' A function for computing rolling products of time-series data.
##'
##' @param x matrix or xts object. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @return An object of the same class and dimension as \code{x} with the rolling products.
##' @examples
##' n_vars <- 3
##' n_obs <- 15
##' x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##' 
##' # rolling products
##' result <- roll_prod(x, 5)
##' 
##' # rolling products with exponential decay
##' weights <- 0.9 ^ (5:1)
##' result <- roll_prod(x, 5, weights)
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
##' A function for computing rolling means of time-series data.
##'
##' @param x matrix or xts object. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @return An object of the same class and dimension as \code{x} with the rolling means.
##' @examples
##' n_vars <- 3
##' n_obs <- 15
##' x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##' 
##' # rolling means
##' result <- roll_mean(x, 5)
##' 
##' # rolling means with exponential decay
##' weights <- 0.9 ^ (5:1)
##' result <- roll_mean(x, 5, weights)
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

##' Rolling Variances
##'
##' A function for computing rolling variances of time-series data.
##'
##' @param x matrix or xts object. Rows are observations and columns are variables.
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
##' @details The denominator used gives an unbiased estimate of the variance, so if the weights are the 
##' default then the divisor \code{n - 1} is obtained.
##' @return An object of the same class and dimension as \code{x} with the rolling variances.
##' @examples
##' n_vars <- 3
##' n_obs <- 15
##' x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##' 
##' # rolling variances
##' result <- roll_var(x, 5)
##' 
##' # rolling variances with exponential decay
##' weights <- 0.9 ^ (5:1)
##' result <- roll_var(x, 5, weights)
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
##' A function for computing rolling standard deviations of time-series data.
##'
##' @param x matrix or xts object. Rows are observations and columns are variables.
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
##' @details The denominator used gives an unbiased estimate of the standard deviation, so if the weights are the 
##' default then the divisor \code{n - 1} is obtained.
##' @return An object of the same class and dimension as \code{x} with the rolling standard deviations.
##' @examples
##' n_vars <- 3
##' n_obs <- 15
##' x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##' 
##' # rolling standard deviations
##' result <- roll_sd(x, 5)
##' 
##' # rolling standard deviations with exponential decay
##' weights <- 0.9 ^ (5:1)
##' result <- roll_sd(x, 5, weights)
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
##' A function for computing rolling scaling and centering of time-series data.
##'
##' @param x matrix or xts object. Rows are observations and columns are variables.
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
##' The denominator used gives an unbiased estimate of the standard deviation, so if the weights are the 
##' default then the divisor \code{n - 1} is obtained.
##' @return An object of the same class and dimension as \code{x} with the rolling scaling and centering.
##' @examples
##' n_vars <- 3
##' n_obs <- 15
##' x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##' 
##' # rolling z-scores
##' result <- roll_scale(x, 5)
##'
##' # rolling z-scores with exponential decay
##' weights <- 0.9 ^ (5:1)
##' result <- roll_scale(x, 5, weights)
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

##' Rolling Covariance Matrices
##'
##' A function for computing rolling covariance matrices of time-series data.
##' 
##' @param x matrix or xts object. Rows are observations and columns are variables.
##' @param y matrix or xts object. Rows are observations and columns are variables.
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
##' @details The denominator used gives an unbiased estimate of the covariance, so if the weights are the 
##' default then the divisor \code{n - 1} is obtained.
##' @return A cube with each slice the rolling covariance matrix.
##' @examples
##' n_vars <- 3
##' n_obs <- 15
##' x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##' 
##' # rolling covariance matrices
##' result <- roll_cov(x, width = 5)
##' 
##' # rolling covariance matrices with exponential decay
##' weights <- 0.9 ^ (5:1)
##' result <- roll_cov(x, width = 5, weights = weights)
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

##' Rolling Correlation Matrices
##'
##' A function for computing rolling correlation matrices of time-series data.
##' 
##' @param x matrix or xts object. Rows are observations and columns are variables.
##' @param y matrix or xts object. Rows are observations and columns are variables.
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
##' @details The denominator used gives an unbiased estimate of the covariance, so if the weights are the 
##' default then the divisor \code{n - 1} is obtained.
##' @return A cube with each slice the rolling correlation matrix.
##' @examples
##' n_vars <- 3
##' n_obs <- 15
##' x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##' 
##' # rolling correlation matrices
##' result <- roll_cor(x, width = 5)
##' 
##' # rolling correlation matrices with exponential decay
##' weights <- 0.9 ^ (5:1)
##' result <- roll_cor(x, width = 5, weights = weights)
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
##' A function for computing rolling linear models of time-series data.
##' 
##' @param x matrix or xts object. Rows are observations and columns are the independent variables.
##' @param y matrix or xts object. Rows are observations and columns are the dependent variables.
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
##' \item{coefficients}{A list of objects with the rolling coefficients for each \code{y}.
##' An object is the same class and dimension (with an added column for the intercept) as \code{x}.}
##' \item{r.squared}{A list of objects with the rolling r-squareds for each \code{y}.
##' An object is the same class as \code{x}.}
##' \item{std.error}{A list of objects with the rolling standard errors for each \code{y}.
##' An object is the same class and dimension (with an added column for the intercept) as \code{x}.}
##' @examples
##' n_vars <- 3
##' n_obs <- 15
##' x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##' y <- matrix(rnorm(n_obs), nrow = n_obs, ncol = 1)
##' 
##' # rolling regressions
##' result <- roll_lm(x, y, 5)
##' 
##' # rolling regressions with exponential decay
##' weights <- 0.9 ^ (5:1)
##' result <- roll_lm(x, y, 5, weights)
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