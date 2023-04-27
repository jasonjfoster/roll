dimnames_lm_x <- function(dimnames, n_cols_x, intercept) {
  
  if (intercept && (length(dimnames) > 1)) {
    
    result <- list(dimnames[[1]], c("(Intercept)", dimnames[[2]]))
    
  } else if (!intercept && (length(dimnames) > 1)) {
    
    result <- list(dimnames[[1]], dimnames[[2]])
    
  } else if (intercept) {
    
    result <- list(NULL, c("(Intercept)", paste0("x", rep(1:(n_cols_x - 1)))))
    
  } else {
    
    result <- list(NULL, paste0("x", rep(1:n_cols_x)))
    
  }
  
  return(result)
  
}

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
    
    if (!is.matrix(x)) {
      
      temp_attr <- attributes(x)
      x <- as.matrix(zoo::coredata(x))
      attr(x, "dimnames") <- NULL
      attr(x, "index") <- temp_attr[["index"]]
      attr(x, "class") <- temp_attr[["class"]]
      
    }
    
    if (!is.matrix(y)) {
      
      temp_attr <- attributes(y)
      y <- as.matrix(zoo::coredata(y))
      attr(y, "dimnames") <- NULL
      attr(y, "index") <- temp_attr[["index"]]
      attr(y, "class") <- temp_attr[["class"]]
      
    }
    
    n_rows_xy <- nrow(x)
    n_cols_x <- ncol(x)
    n_cols_y <- ncol(y)
    result <- array(as.numeric(NA), c(n_cols_x, n_cols_y, n_rows_xy))
    
    for (i in 1:n_rows_xy) {
      
      result[ , , i] <- f(x[max(1, i - width + 1):i, , drop = FALSE],
                          y[max(1, i - width + 1):i, , drop = FALSE])
      
    }
    
    attr(result, "dim") <- c(n_cols_x, n_cols_y, n_rows_xy)
    
    x_dimnames <- dimnames(x)
    y_dimnames <- dimnames(y)
    attr(result, "dimnames") <- list(x_dimnames[[2]], y_dimnames[[2]], NULL)
    
  } else {
    
    n_rows_xy <- length(x)
    result <- array(as.numeric(NA), n_rows_xy)
    
    for (i in 1:n_rows_xy) {
      
      result[i] <- f(x[max(1, i - width + 1):i],
                     y[max(1, i - width + 1):i])
      
    }
    
    result <- as.numeric(result)
    
  }
  
  return(result)
  
}

rollapplyr_lm <- function(x, y, width, intercept) {
  
  if (!requireNamespace("zoo", quietly = TRUE)) {
    stop("zoo package required for this function")
  }
  
  if (is.matrix(x) || is.matrix(y) || intercept) {
    
    if (!is.matrix(x)) {
      
      temp_attr <- attributes(x)
      x <- as.matrix(zoo::coredata(x))
      attr(x, "dimnames") <- NULL
      attr(x, "index") <- temp_attr[["index"]]
      attr(x, "class") <- temp_attr[["class"]]
      
    }
    
    if (!is.matrix(y)) {
      
      temp_attr <- attributes(y)
      y <- as.matrix(zoo::coredata(y))
      attr(y, "dimnames") <- NULL
      attr(y, "index") <- temp_attr[["index"]]
      attr(y, "class") <- temp_attr[["class"]]
      
    }
    
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
      
    } else if (zoo::is.zoo(y)) {
      
      x_attr <- attributes(y)
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
      
      if (nrow(coef(summary_fit)) == n_cols_x) {
        
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
      
    }
    
    attr(result[["coefficients"]], "dim") <- c(n_rows_xy, n_cols_x)
    attr(result[["r.squared"]], "dim") <- c(n_rows_xy, 1)
    attr(result[["std.error"]], "dim") <- c(n_rows_xy, n_cols_x)
    
    x_dimnames <- dimnames(x)
    y_dimnames <- dimnames(y)
    
    attr(result[["coefficients"]], "dimnames") <- dimnames_lm_x(x_dimnames, n_cols_x, intercept)
    if (length(x_dimnames) > 1) {
      attr(result[["r.squared"]], "dimnames") <- list(x_dimnames[[1]], "R-squared")
    } else {
      attr(result[["r.squared"]], "dimnames") <- list(NULL, "R-squared")
    }
    attr(result[["std.error"]], "dimnames") <- dimnames_lm_x(x_dimnames, n_cols_x, intercept)
    
  } else {
    
    n_rows_xy <- length(x)
    n_cols_x <- 1
    
    result <- list("coefficients" = rep(as.numeric(NA), n_rows_xy),
                   "r.squared" = rep(as.numeric(NA), n_rows_xy),
                   "std.error" = rep(as.numeric(NA), n_rows_xy))
    
    if (zoo::is.zoo(x)) {

      x_attr <- attributes(x)
      x_attr[["dim"]] <- NULL
      x_attr[["dimnames"]] <- NULL

    } else if (zoo::is.zoo(y)) {

      x_attr <- attributes(y)
      x_attr[["dim"]] <- NULL
      x_attr[["dimnames"]] <- NULL

    }
    
    for (i in 1:n_rows_xy) {
      
      x_subset <- x[max(1, i - width + 1):i]
      y_subset <- y[max(1, i - width + 1):i]
      data <- as.data.frame(cbind(y_subset, x_subset))
      
      if (intercept) {
        fit <- lm(reformulate(termlabels = ".", response = names(data)[1]), data = data)
      } else {
        fit <- lm(reformulate(termlabels = ".-1", response = names(data)[1]), data = data)
      }
      
      summary_fit <- summary(fit)
      summary_fit_coef <- coef(summary_fit)[ , "Estimate"]
      
      if (nrow(coef(summary_fit)) == n_cols_x) {
        
        result[["coefficients"]][i] <- summary_fit_coef
        
        # "In summary.lm(fit) : essentially perfect fit: summary may be unreliable"
        if (!(isTRUE(all.equal(as.numeric(rep(summary_fit_coef[1], length(y_subset))), as.numeric(y_subset))) &&
              isTRUE(all.equal(as.numeric(summary_fit_coef[-1]), rep(0, length(summary_fit_coef[-1])))))) {
          
          result[["r.squared"]][i] <- summary_fit$r.squared
          result[["std.error"]][i] <- coef(summary_fit)[ , "Std. Error"]
          
        } else {
          
          result[["r.squared"]][i] <- as.numeric(NA)
          result[["std.error"]][i] <- as.numeric(NA)
          
        }
        
      }
      
    }
    
    if (exists("x_attr")) {
      
      attributes(result[["coefficients"]]) <- x_attr
      attributes(result[["r.squared"]]) <- x_attr
      attributes(result[["std.error"]]) <- x_attr
      
    }
    
  }
  
  return(result)
  
}
