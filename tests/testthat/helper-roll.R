n_vars <- 4
n_obs <- 20
lambda <- 0.9

# test data
set.seed(5640)
dates <- seq(Sys.Date(), length.out = n_obs, by = "-1 day")

if (requireNamespace("zoo", quietly = TRUE)) {
  
  test_data <- list(zoo::zoo(matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars), dates),
                    zoo::zoo(matrix(rev(rep(1:n_vars, times = n_vars, each = n_obs / n_vars)) / 1000,
                                    nrow = n_obs, ncol = n_vars), dates))
  
} else {
  
  test_data <- list(matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars),
                    matrix(rev(rep(1:n_vars, times = n_vars, each = n_obs / n_vars)) / 1000,
                           nrow = n_obs, ncol = n_vars))
  
}

test_data <- lapply(test_data, setNames, paste0("x", rep(1:n_vars)))
test_data[[3]] <- matrix(rnorm(n_obs * n_vars), nrow = n_obs, ncol = n_vars)

set.seed(5640)
idx <- sample(1:(n_obs * n_vars), n_obs / 4)
test_data[[3]][idx] <- as.numeric(NA)

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
  result <- as.numeric(NA)
  
  temp <- scale(x, center = center, scale = scale)
  
  if ((scale && (n_rows_x > 1)) || !scale) {
    result <- temp[length(temp)]
  }
  
  return(result)
  
}

rollapplyr_cube <- function(f, x, y, width) {
  
  if (is.matrix(x) || is.matrix(y)) {
    
    x <- as.matrix(x)
    y <- as.matrix(y)
    
    n_rows_xy <- nrow(x)
    n_cols_x <- ncol(x)
    n_cols_y <- ncol(y)
    result <- array(as.numeric(NA), c(n_cols_x, n_cols_y, n_rows_xy))
    
    for (i in 1:n_rows_xy) {
      
      temp <- f(x[max(1, i - width + 1):i, , drop = FALSE],
                y[max(1, i - width + 1):i, , drop = FALSE])
      
      if (!anyNA(temp)) {
        result[ , , i] <- temp
      }
      
    }
    
  } else {
    
    n_rows_xy <- length(x)
    result <- array(NA, n_rows_xy)
    
    for (i in 1:n_rows_xy) {
      
      temp <- f(x[max(1, i - width + 1):i],
                y[max(1, i - width + 1):i])
      
      if (!anyNA(temp)) {
        result[i] <- temp
      }
      
    }
    
  }
  
  return(result)
  
}

rollapplyr_lm <- function(x, y, width, intercept) {
  
  n_rows_xy <- nrow(x)
  n_cols_x <- ncol(x)
  
  if (intercept) {
    n_cols_x <- n_cols_x + 1
  }
  
  result <- list("coefficients" = matrix(as.numeric(NA), n_rows_xy, n_cols_x),
                 "r.squared" = matrix(as.numeric(NA), n_rows_xy, 1),
                 "std.error" = matrix(as.numeric(NA), n_rows_xy, n_cols_x))
  
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
    summary_fit_coef <- coef(summary_fit)[ , "Estimate"]
    
    if ((nrow(coef(summary_fit)) == n_cols_x)) {
      
      result[["coefficients"]][i, ] <- summary_fit_coef
      
      # "In summary.lm(fit) : essentially perfect fit: summary may be unreliable"
      if (!(isTRUE(all.equal(as.numeric(rep(summary_fit_coef[1], length(y_subset))), as.numeric(y_subset))) &&
            isTRUE(all.equal(as.numeric(summary_fit_coef[-1]), rep(0, length(summary_fit_coef[-1])))))) {
        
        result[["r.squared"]][i, ] <- summary_fit$r.squared
        result[["std.error"]][i, ] <- coef(summary_fit)[ , "Std. Error"]
        
      } else {
        
        result[["r.squared"]][i, ] <- as.numeric(NA)
        result[["std.error"]][i, ] <- as.numeric(NA)
        
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