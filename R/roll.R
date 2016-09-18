##' Rolling Means
##'
##' A parallel function for computing rolling means of time-series data.
##'
##' @param data matrix or xts object. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is NA.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param parallel_for character. Executes a "for" loop in which iterations run in parallel by
##' \code{rows} or \code{cols}.
##' @return An object of the same class and dimension as \code{data} with the rolling means.
##' @seealso \code{\link[RcppParallel]{setThreadOptions}} for thread options via RcppParallel.
##' @examples
##' n_vars <- 10
##' n_obs <- 1000
##' data <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##' 
##' # 252-day rolling mean
##' result <- roll_mean(data, 252)
##' 
##' # Equivalent to 'na.rm = TRUE'
##' result <- roll_mean(data, 252, min_obs = 1)
##' 
##' # Expanding window
##' result <- roll_mean(data, n_obs, min_obs = 1)
##' 
##' # Exponential decay
##' weights <- 0.9 ^ (251:0)
##' result <- roll_mean(data, 252, weights, min_obs = 1)
##' @export
roll_mean <- function(data, width, weights = rep(1, width),
                      min_obs = width, complete_obs = FALSE, na_restore = FALSE,
                      parallel_for = c("rows", "cols")) {
  return(.Call('roll_roll_mean', PACKAGE = 'roll',
               data,
               as.integer(width),
               as.numeric(weights),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.character(match.arg(parallel_for))
  ))
}

##' Rolling Variances
##'
##' A parallel function for computing rolling variances of time-series data.
##'
##' @param data matrix or xts object. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param center logical. If \code{TRUE} then the weighted mean of each variable is used,
##' if \code{FALSE} then zero is used.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is NA.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param parallel_for character. Executes a "for" loop in which iterations run in parallel by
##' \code{rows} or \code{cols}.
##' @return An object of the same class and dimension as \code{data} with the rolling variances.
##' @seealso \code{\link[RcppParallel]{setThreadOptions}} for thread options via RcppParallel.
##' @examples
##' n_vars <- 10
##' n_obs <- 1000
##' data <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##' 
##' # 252-day rolling variance
##' result <- roll_var(data, 252)
##' 
##' # Equivalent to 'na.rm = TRUE'
##' result <- roll_var(data, 252, min_obs = 1)
##' 
##' # Expanding window
##' result <- roll_var(data, n_obs, min_obs = 1)
##' 
##' # Exponential decay
##' weights <- 0.9 ^ (251:0)
##' result <- roll_var(data, 252, weights, min_obs = 1)
##' @export
roll_var <- function(data, width, weights = rep(1, width), center = TRUE,
                     min_obs = width, complete_obs = FALSE, na_restore = FALSE,
                     parallel_for = c("rows", "cols")) {
  return(.Call('roll_roll_var', PACKAGE = 'roll',
               data,
               as.integer(width),
               as.numeric(weights),
               as.logical(center),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.character(match.arg(parallel_for))
  ))
}

##' Rolling Standard Deviations
##'
##' A parallel function for computing rolling standard deviations of time-series data.
##'
##' @param data matrix or xts object. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param center logical. If \code{TRUE} then the weighted mean of each variable is used,
##' if \code{FALSE} then zero is used.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is NA.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param parallel_for character. Executes a "for" loop in which iterations run in parallel by
##' \code{rows} or \code{cols}.
##' @return An object of the same class and dimension as \code{data} with the rolling standard deviations.
##' @seealso \code{\link[RcppParallel]{setThreadOptions}} for thread options via RcppParallel.
##' @examples
##' n_vars <- 10
##' n_obs <- 1000
##' data <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##' 
##' # 252-day rolling standard deviation
##' result <- roll_sd(data, 252)
##' 
##' # Equivalent to 'na.rm = TRUE'
##' result <- roll_sd(data, 252, min_obs = 1)
##' 
##' # Expanding window
##' result <- roll_sd(data, n_obs, min_obs = 1)
##' 
##' # Exponential decay
##' weights <- 0.9 ^ (251:0)
##' result <- roll_sd(data, 252, weights, min_obs = 1)
##' @export
roll_sd <- function(data, width, weights = rep(1, width), center = TRUE,
                    min_obs = width, complete_obs = FALSE, na_restore = FALSE,
                    parallel_for = c("rows", "cols")) {
  return(.Call('roll_roll_sd', PACKAGE = 'roll',
               data,
               as.integer(width),
               as.numeric(weights),
               as.logical(center),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.character(match.arg(parallel_for))
  ))
}

##' Rolling Scaling and Centering
##'
##' A parallel function for computing rolling scaling and centering of time-series data.
##'
##' @param data matrix or xts object. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param center logical. If \code{TRUE} then the weighted mean of each variable is used,
##' if \code{FALSE} then zero is used.
##' @param scale logical. If \code{TRUE} then the weighted standard deviation of each variable is used,
##' if \code{FALSE} then no scaling is done.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is NA.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then each value is used.
##' @param na_restore logical. Should missing values be restored?
##' @param parallel_for character. Executes a "for" loop in which iterations run in parallel by
##' \code{rows} or \code{cols}.
##' @return An object of the same class and dimension as \code{data} with the rolling scaling and centering.
##' @seealso \code{\link[RcppParallel]{setThreadOptions}} for thread options via RcppParallel.
##' @examples
##' n_vars <- 10
##' n_obs <- 1000
##' data <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##' 
##' # 252-day rolling z-score
##' result <- roll_scale(data, 252)
##' 
##' # Equivalent to 'na.rm = TRUE'
##' result <- roll_scale(data, 252, min_obs = 1)
##' 
##' # Expanding window
##' result <- roll_scale(data, n_obs, min_obs = 1)
##' 
##' # Exponential decay
##' weights <- 0.9 ^ (251:0)
##' result <- roll_scale(data, 252, weights, min_obs = 1)
##' @export
roll_scale <- function(data, width, weights = rep(1, width), center = TRUE, scale = TRUE,
                       min_obs = width, complete_obs = FALSE, na_restore = FALSE,
                       parallel_for = c("rows", "cols")) {
  return(.Call('roll_roll_scale', PACKAGE = 'roll',
               data,
               as.integer(width),
               as.numeric(weights),
               as.logical(center),
               as.logical(scale),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.character(match.arg(parallel_for))
  ))
}

##' Rolling Covariance Matrices
##'
##' A parallel function for computing rolling covariance matrices of time-series data.
##' 
##' @param data matrix or xts object. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param center logical. If \code{TRUE} then the weighted mean of each variable is used,
##' if \code{FALSE} then zero is used.
##' @param scale logical. If \code{TRUE} then the weighted standard deviation of each variable is used,
##' if \code{FALSE} then no scaling is done.
##' @param min_obs integer. Minimum number of observations required to have a value within a window, 
##' otherwise result is NA.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then pairwise is used.
##' @param na_restore logical. Should missing values be restored?
##' @param parallel_for character. Executes a "for" loop in which iterations run in parallel by
##' \code{rows} or \code{cols}.
##' @return A cube with each slice the rolling covariance matrix.
##' @seealso \code{\link[RcppParallel]{setThreadOptions}} for thread options via RcppParallel.
##' @examples
##' n_vars <- 10
##' n_obs <- 1000
##' data <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##' 
##' # 252-day rolling covariance matrix
##' result <- roll_cov(data, 252)
##' 
##' # Equivalent to 'na.rm = TRUE'
##' result <- roll_cov(data, 252, min_obs = 1)
##' 
##' # Expanding window
##' result <- roll_cov(data, n_obs, min_obs = 1)
##' 
##' # Exponential decay
##' weights <- 0.9 ^ (251:0)
##' result <- roll_cov(data, 252, weights, min_obs = 1)
##' @export
roll_cov <- function(data, width, weights = rep(1, width), center = TRUE, scale = FALSE,
                     min_obs = width, complete_obs = TRUE, na_restore = FALSE,
                     parallel_for = c("rows", "cols")) {
  return(.Call('roll_roll_cov', PACKAGE = 'roll',
               data,
               as.integer(width),
               as.numeric(weights),
               as.logical(center),
               as.logical(scale),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.character(match.arg(parallel_for))
  ))
}

##' Rolling Correlation Matrices
##'
##' A parallel function for computing rolling correlation matrices of time-series data.
##' 
##' @param data matrix or xts object. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param center logical. If \code{TRUE} then the weighted mean of each variable is used,
##' if \code{FALSE} then zero is used.
##' @param scale logical. If \code{TRUE} then the weighted standard deviation of each variable is used,
##' if \code{FALSE} then no scaling is done.
##' @param min_obs integer. Minimum number of observations required to have a value within a window, 
##' otherwise result is NA.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then pairwise is used.
##' @param na_restore logical. Should missing values be restored?
##' @param parallel_for character. Executes a "for" loop in which iterations run in parallel by
##' \code{rows} or \code{cols}.
##' @return A cube with each slice the rolling correlation matrix.
##' @seealso \code{\link[RcppParallel]{setThreadOptions}} for thread options via RcppParallel.
##' @examples
##' n_vars <- 10
##' n_obs <- 1000
##' data <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##' 
##' # 252-day rolling correlation matrix
##' result <- roll_cor(data, 252)
##' 
##' # Equivalent to 'na.rm = TRUE'
##' result <- roll_cor(data, 252, min_obs = 1)
##' 
##' # Expanding window
##' result <- roll_cor(data, n_obs, min_obs = 1)
##' 
##' # Exponential decay
##' weights <- 0.9 ^ (251:0)
##' result <- roll_cor(data, 252, weights, min_obs = 1)
##' @export
roll_cor <- function(data, width, weights = rep(1, width), center = TRUE, scale = TRUE,
                     min_obs = width, complete_obs = TRUE, na_restore = FALSE,
                     parallel_for = c("rows", "cols")) {
  return(.Call('roll_roll_cov', PACKAGE = 'roll',
               data,
               as.integer(width),
               as.numeric(weights),
               as.logical(center),
               as.logical(scale),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.character(match.arg(parallel_for))
  ))
}

##' Rolling Linear Models
##'
##' A parallel function for computing rolling linear models of time-series data.
##' 
##' @param x matrix or xts object. Rows are observations and columns are the independent variables.
##' @param y matrix or xts object. Rows are observations and columns are the dependent variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param intercept logical. Either \code{TRUE} to include or \code{FALSE} to remove the intercept.
##' @param center logical. \code{center = z} is shorthand for \code{center_x = z} and
##' \code{center_y = z}, where \code{z} is either \code{TRUE} or \code{FALSE}.
##' @param center_x logical. If \code{TRUE} then the weighted mean of each \code{x} variable is used,
##' if \code{FALSE} then zero is used.
##' @param center_y logical. Analogous to \code{center_x}.
##' @param scale logical. \code{scale = z} is shorthand for \code{scale_x = z} and
##' \code{scale_y = z}, where \code{z} is either \code{TRUE} or \code{FALSE}.
##' @param scale_x logical. If \code{TRUE} then the weighted standard deviation of each \code{x} 
##' variable is used, if \code{FALSE} then no scaling is done.
##' @param scale_y logical. Analogous to \code{scale_x}.
##' @param min_obs integer. Minimum number of observations required to have a value within a window, 
##' otherwise result is NA.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then pairwise is used.
##' @param na_restore logical. Should missing values be restored?
##' @param parallel_for character. Executes a "for" loop in which iterations run in parallel by
##' \code{rows} or \code{cols}.
##' @return A list containing the following components:
##' \item{coefficients}{A list of objects with the rolling coefficients for each \code{y}.
##' An object is the same class and dimension (with an added column for the intercept) as \code{x}.}
##' \item{r.squared}{A list of objects with the rolling r-squareds for each \code{y}.
##' An object is the same class as \code{x}.}
##' @note If users are already taking advantage of parallelism using multithreaded BLAS/LAPACK
##' libraries, then limit the number of cores in the RcppParallel package to one with the
##' \code{\link[RcppParallel]{setThreadOptions}} function.
##' @seealso \code{\link[RcppParallel]{setThreadOptions}} for thread options via RcppParallel.
##' @examples
##' n_vars <- 10
##' n_obs <- 1000
##' x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##' y <- matrix(rnorm(n_obs), nrow = n_obs, ncol = 1)
##' 
##' # 252-day rolling regression
##' result <- roll_lm(x, y, 252)
##' 
##' # Equivalent to 'na.rm = TRUE'
##' result <- roll_lm(x, y, 252, min_obs = 1)
##' 
##' # Expanding window
##' result <- roll_lm(x, y, n_obs, min_obs = 1)
##' 
##' # Exponential decay
##' weights <- 0.9 ^ (251:0)
##' result <- roll_lm(x, y, 252, weights, min_obs = 1)
##' @export
roll_lm <- function(x, y, width, weights = rep(1, width), intercept = TRUE, 
                    center = FALSE, center_x = center, center_y = center,
                    scale = FALSE, scale_x = scale, scale_y = scale,
                    min_obs = width, complete_obs = TRUE,
                    na_restore = FALSE, parallel_for = c("rows", "cols")) {
  return(.Call('roll_roll_lm', PACKAGE = 'roll',
               x,
               y,
               as.integer(width),
               as.numeric(weights),
               as.logical(intercept),
               as.logical(center_x),
               as.logical(center_y),
               as.logical(scale_x),
               as.logical(scale_y),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.character(match.arg(parallel_for))
  ))
}

##' Rolling Eigenvalues and Eigenvectors
##'
##' A parallel function for computing rolling eigenvalues and eigenvectors of time-series data.
##' 
##' @param data matrix or xts object. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param center logical. If \code{TRUE} then the weighted mean of each variable is used,
##' if \code{FALSE} then zero is used.
##' @param scale logical. If \code{TRUE} then the weighted standard deviation of each variable is used,
##' if \code{FALSE} then no scaling is done.
##' @param min_obs integer. Minimum number of observations required to have a value within a window, 
##' otherwise result is NA.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then pairwise is used.
##' @param na_restore logical. Should missing values be restored?
##' @param parallel_for character. Executes a "for" loop in which iterations run in parallel by
##' \code{rows} or \code{cols}.
##' @return A list containing the following components:
##' \item{values}{An object of the same class and dimension as \code{data} with the rolling eigenvalues.}
##' \item{vectors}{A cube with each slice the rolling eigenvectors.}
##' @note If users are already taking advantage of parallelism using multithreaded BLAS/LAPACK
##' libraries, then limit the number of cores in the RcppParallel package to one with the
##' \code{\link[RcppParallel]{setThreadOptions}} function.
##' @seealso \code{\link[RcppParallel]{setThreadOptions}} for thread options via RcppParallel.
##' @examples
##' n_vars <- 10
##' n_obs <- 1000
##' data <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##' 
##' # 252-day rolling eigenvalues and eigenvectors
##' result <- roll_eigen(data, 252)
##' 
##' # Equivalent to 'na.rm = TRUE'
##' result <- roll_eigen(data, 252, min_obs = 1)
##' 
##' # Expanding window
##' result <- roll_eigen(data, n_obs, min_obs = 1)
##' 
##' # Exponential decay
##' weights <- 0.9 ^ (251:0)
##' result <- roll_eigen(data, 252, weights, min_obs = 1)
##' @export
roll_eigen <- function(data, width, weights = rep(1, width), center = TRUE, scale = FALSE,
                       min_obs = width, complete_obs = TRUE, na_restore = FALSE,
                       parallel_for = c("rows", "cols")) {
  return(.Call('roll_roll_eigen', PACKAGE = 'roll',
               data,
               as.integer(width),
               as.numeric(weights),
               as.logical(center),
               as.logical(scale),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.character(match.arg(parallel_for))
  ))
}

##' Rolling Principal Component Regressions 
##'
##' A parallel function for computing rolling principal component regressions of time-series data.
##' 
##' @param x matrix or xts object. Rows are observations and columns are the independent variables.
##' @param y matrix or xts object. Rows are observations and columns are the dependent variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param comps integer vector. Select a subset of principal components.
##' @param intercept logical. Either \code{TRUE} to include or \code{FALSE} to remove the intercept.
##' @param center logical. \code{center = z} is shorthand for \code{center_x = z} and
##' \code{center_y = z}, where \code{z} is either \code{TRUE} or \code{FALSE}.
##' @param center_x logical. If \code{TRUE} then the weighted mean of each \code{x} variable is used,
##' if \code{FALSE} then zero is used.
##' @param center_y logical. Analogous to \code{center_x}.
##' @param scale logical. \code{scale = z} is shorthand for \code{scale_x = z} and
##' \code{scale_y = z}, where \code{z} is either \code{TRUE} or \code{FALSE}.
##' @param scale_x logical. If \code{TRUE} then the weighted standard deviation of each \code{x} 
##' variable is used, if \code{FALSE} then no scaling is done.
##' @param scale_y logical. Analogous to \code{scale_x}.
##' @param min_obs integer. Minimum number of observations required to have a value within a window, 
##' otherwise result is NA.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then pairwise is used.
##' @param na_restore logical. Should missing values be restored?
##' @param parallel_for character. Executes a "for" loop in which iterations run in parallel by
##' \code{rows} or \code{cols}.
##' @return A list containing the following components:
##' \item{coefficients}{A list of objects with the rolling coefficients for each \code{y}.
##' An object is the same class and dimension (with an added column for the intercept) as \code{x}.}
##' \item{r.squared}{A list of objects with the rolling r-squareds for each \code{y}.
##' An object is the same class as \code{x}.}
##' @note If users are already taking advantage of parallelism using multithreaded BLAS/LAPACK
##' libraries, then limit the number of cores in the RcppParallel package to one with the
##' \code{\link[RcppParallel]{setThreadOptions}} function.
##' @seealso \code{\link[RcppParallel]{setThreadOptions}} for thread options via RcppParallel.
##' @examples
##' n_vars <- 10
##' n_obs <- 1000
##' x <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##' y <- matrix(rnorm(n_obs), nrow = n_obs, ncol = 1)
##' 
##' # 252-day rolling principal component regression
##' result <- roll_pcr(x, y, 252, comps = 1)
##' 
##' # Equivalent to 'na.rm = TRUE'
##' result <- roll_pcr(x, y, 252, comps = 1, min_obs = 1)
##' 
##' # Expanding window
##' result <- roll_pcr(x, y, n_obs, comps = 1, min_obs = 1)
##' 
##' # Exponential decay
##' weights <- 0.9 ^ (251:0)
##' result <- roll_pcr(x, y, 252, comps = 1, weights, min_obs = 1)
##' @export
roll_pcr <- function(x, y, width, comps = 1:ncol(x), weights = rep(1, width), intercept = TRUE, 
                     center = FALSE, center_x = center, center_y = center,
                     scale = FALSE, scale_x = scale, scale_y = scale,
                     min_obs = width, complete_obs = TRUE,
                     na_restore = FALSE, parallel_for = c("rows", "cols")) {
  return(.Call('roll_roll_pcr', PACKAGE = 'roll',
               x,
               y,
               as.integer(width),
               as.numeric(comps),
               as.numeric(weights),
               as.logical(intercept),
               as.logical(center_x),
               as.logical(center_y),
               as.logical(scale_x),
               as.logical(scale_y),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.character(match.arg(parallel_for))
  ))
}

##' Rolling Variance Inflation Factors
##'
##' A parallel function for computing rolling variance inflation factors of time-series data.
##' 
##' @param data matrix or xts object. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param center logical. If \code{TRUE} then the weighted mean of each variable is used,
##' if \code{FALSE} then zero is used.
##' @param scale logical. If \code{TRUE} then the weighted standard deviation of each variable is used,
##' if \code{FALSE} then no scaling is done.
##' @param min_obs integer. Minimum number of observations required to have a value within a window, 
##' otherwise result is NA.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then pairwise is used.
##' @param na_restore logical. Should missing values be restored?
##' @param parallel_for character. Executes a "for" loop in which iterations run in parallel by
##' \code{rows} or \code{cols}.
##' @return An object of the same class and dimension as \code{data} with the rolling variance
##' inflation factors.
##' @note If users are already taking advantage of parallelism using multithreaded BLAS/LAPACK
##' libraries, then limit the number of cores in the RcppParallel package to one with the
##' \code{\link[RcppParallel]{setThreadOptions}} function.
##' @seealso \code{\link[RcppParallel]{setThreadOptions}} for thread options via RcppParallel.
##' @examples
##' n_vars <- 10
##' n_obs <- 1000
##' data <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)
##' 
##' # 252-day rolling variance inflation factors
##' result <- roll_vif(data, 252)
##' 
##' # Equivalent to 'na.rm = TRUE'
##' result <- roll_vif(data, 252, min_obs = 1)
##' 
##' # Expanding window
##' result <- roll_vif(data, n_obs, min_obs = 1)
##' 
##' # Exponential decay
##' weights <- 0.9 ^ (251:0)
##' result <- roll_vif(data, 252, weights, min_obs = 1)
##' @export
roll_vif <- function(data, width, weights = rep(1, width), center = FALSE, scale = FALSE,
                     min_obs = width, complete_obs = TRUE, na_restore = FALSE,
                     parallel_for = c("rows", "cols")) {
  return(.Call('roll_roll_vif', PACKAGE = 'roll',
               data,
               as.integer(width),
               as.numeric(weights),
               as.logical(center),
               as.logical(scale),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.character(match.arg(parallel_for))
  ))
}



