n_vars <- 4
n_obs <- 20
lambda <- 0.9

# test data
set.seed(5640)
dates <- seq(Sys.Date(), length.out = n_obs, by = "-1 day")
test_data <- list(zoo::zoo(matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars), dates),
                  zoo::zoo(matrix(rev(rep(1:n_vars, times = n_vars, each = n_obs / n_vars)) / 1000,
                                  nrow = n_obs, ncol = n_vars), dates))
test_data <- lapply(test_data, setNames, paste0("x", rep(1:n_vars)))
test_data[[3]] <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)

set.seed(5640)
idx <- sample(1:(n_obs * n_vars), n_obs / 4)
test_data[[3]][idx] <- NA

# test arguments
test_data_x <- lapply(test_data, function(x){x[ , 1:2]})
test_data_y <- lapply(test_data, function(x){x[ , 3:4, drop = FALSE]})
test_data_null <- c(test_data_x, list(NULL))
test_width <- c(1, 2, 10, 20)
test_intercept <- c(TRUE, FALSE)
test_center <- c(TRUE, FALSE)
test_scale <- c(TRUE, FALSE)
test_min_obs <- c(1, 2, 10, 20)
test_complete_obs <- c(TRUE, FALSE)
test_na_restore <- c(TRUE, FALSE)
test_online <- c(TRUE, FALSE)

# test functions
scale_z <- function(x, center = TRUE, scale = TRUE) {
  
  n_rows_x <- length(x)
  
  result <- scale(x, center = center, scale = scale)
  
  if ((scale && (n_rows_x > 1)) || !scale) {
    result <- result[length(result)]
  } else {
    result <- NA
  }
  
  return(result)
  
}

rollapplyr_cube <- function(f, x, y, width) {
  
  n_rows_xy <- nrow(x)
  n_cols_x <- ncol(x)
  r_cube <- array(NA, c(n_cols_x, n_cols_x, n_rows_xy))
  
  for (i in 1:n_rows_xy) {
    
    result <- f(x[max(1, i - width + 1):i, , drop = FALSE],
                y[max(1, i - width + 1):i, , drop = FALSE])
    
    if (!anyNA(result)) {
      r_cube[ , , i] <- result
    }
    
  }
  
  return(r_cube)
  
}

rollapplyr_lm <- function(x, y, width, intercept) {
  
  n_rows_xy <- nrow(x)
  n_cols_x <- ncol(x)
  
  if (intercept) {
    n_cols_x <- n_cols_x + 1
  }
  
  result <- list("coefficients" = matrix(NA, n_rows_xy, n_cols_x),
                 "r.squared" = matrix(NA, n_rows_xy, 1),
                 "std.error" = matrix(NA, n_rows_xy, n_cols_x))
  
  if (zoo::is.zoo(x)) {
    
    x_attr <- attributes(x)
    x_attr[["dim"]] <- NULL
    x_attr[["dimnames"]] <- NULL
    
  }
  
  for (i in 1:n_rows_xy) {
    
    x_subset <- x[max(1, i - width + 1):i, , drop = FALSE]
    y_subset <- y[max(1, i - width + 1):i, , drop = FALSE]
    data <- as.data.frame(cbind(y_subset, x_subset))
    
    if (intercept) {
      fit <- lm(reformulate(termlabels = ".", response = names(data)[1]), data = data)
    } else {
      fit <- lm(reformulate(termlabels = ".-1", response = names(data)[1]), data = data)
    }
    
    summary_fit <- summary(fit)
    
    if ((nrow(coef(summary_fit)) == n_cols_x)) {
      
      result[["coefficients"]][i, ] <- coef(summary_fit)[ , "Estimate"]
      result[["r.squared"]][i, ] <- summary_fit$r.squared
      
      if (!is.na(summary_fit$r.squared)) {
        result[["std.error"]][i, ] <- coef(summary_fit)[ , "Std. Error"]
      }
      
    }
    
  }
  
  if (exists("x_attr")) {
    
    attributes(result[["coefficients"]]) <- x_attr
    attributes(result[["r.squared"]]) <- x_attr
    attributes(result[["std.error"]]) <- x_attr
    
    attr(result[["coefficients"]], "dim") <- c(n_rows_xy, n_cols_x)
    attr(result[["r.squared"]], "dim") <- c(n_rows_xy, 1)
    attr(result[["std.error"]], "dim") <- c(n_rows_xy, n_cols_x)
    
  }
  
  return(result)
  
}