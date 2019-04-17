context("rolling statistics")

n_vars_x <- 4
n_vars_y <- 1
n_obs <- 15
lambda <- 0.9

set.seed(5640)

# test data
dates <- seq(Sys.Date(), length.out = n_obs, by = "-1 day")
test_data <- list(zoo::zoo(matrix(rnorm(n_obs * n_vars_x), nrow = n_obs, ncol = n_vars_x), dates),
                  zoo::zoo(matrix(rev(rep(1:n_vars_x, times = n_vars_x, each = n_obs / n_vars_x)) / 1000,
                                  nrow = n_obs, ncol = n_vars_x), dates))
test_data <- lapply(test_data, setNames, paste0("x", rep(1:n_vars_x)))
test_data[[3]] <- matrix(rnorm(n_obs * n_vars_x), nrow = n_obs, ncol = n_vars_x)

idx <- sample(1:(n_obs * n_vars_x), n_obs / 4)
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

test_that("equal to online algorithm", {
  
  # skip_on_cran()
  
  for (ax in 1:length(test_data_x)) {
    for (b in 1:length(test_width)) {
      
      width <- test_width[b]
      test_weights <- list(lambda ^ ((2 * width):1))
      # test_weights <- list(rep(1, width), lambda ^ (width:1), 1:width,
      #                      rep(1, 2 * width), lambda ^ ((2 * width):1), 1:(width * 2))
      
      for (c in 1:length(test_min_obs)) {
        for (d in 1:length(test_complete_obs)) {
          for (e in 1:length(test_na_restore)) {
            
            expect_equal(roll_any(test_data_x[[ax]] < 0, width,
                                  test_min_obs[c], test_complete_obs[d],
                                  test_na_restore[e], test_online[1]),
                         roll_any(test_data_x[[ax]] < 0, width,
                                  test_min_obs[c], test_complete_obs[d],
                                  test_na_restore[e], test_online[2]))
            
            expect_equal(roll_all(test_data_x[[ax]] < 0, width,
                                  test_min_obs[c], test_complete_obs[d],
                                  test_na_restore[e], test_online[1]),
                         roll_all(test_data_x[[ax]] < 0, width,
                                  test_min_obs[c], test_complete_obs[d],
                                  test_na_restore[e], test_online[2]))
            
            
            for (f in 1:length(test_weights)) {
              
              # 'online = TRUE' is not supported
              expect_equal(roll_min(test_data_x[[ax]], width,
                                    test_weights[[f]], test_min_obs[c],
                                    test_complete_obs[d], test_na_restore[e],
                                    test_online[1]),
                           roll_min(test_data_x[[ax]], width,
                                    test_weights[[f]], test_min_obs[c],
                                    test_complete_obs[d], test_na_restore[e],
                                    test_online[2]))
              
              expect_equal(roll_max(test_data_x[[ax]], width,
                                    test_weights[[f]], test_min_obs[c],
                                    test_complete_obs[d], test_na_restore[e],
                                    test_online[1]),
                           roll_max(test_data_x[[ax]], width,
                                    test_weights[[f]], test_min_obs[c],
                                    test_complete_obs[d], test_na_restore[e],
                                    test_online[2]))
              
              expect_equal(roll_median(test_data_x[[ax]], width,
                                       test_weights[[f]], test_min_obs[c],
                                       test_complete_obs[d], test_na_restore[e],
                                       test_online[1]),
                           roll_median(test_data_x[[ax]], width,
                                       test_weights[[f]], test_min_obs[c],
                                       test_complete_obs[d], test_na_restore[e],
                                       test_online[2]))
              
              expect_equal(roll_sum(test_data_x[[ax]], width,
                                    test_weights[[f]], test_min_obs[c],
                                    test_complete_obs[d], test_na_restore[e],
                                    test_online[1]),
                           roll_sum(test_data_x[[ax]], width,
                                    test_weights[[f]], test_min_obs[c],
                                    test_complete_obs[d], test_na_restore[e],
                                    test_online[2]))
              
              expect_equal(roll_prod(test_data_x[[ax]], width,
                                     test_weights[[f]], test_min_obs[c],
                                     test_complete_obs[d], test_na_restore[e],
                                     test_online[1]),
                           roll_prod(test_data_x[[ax]], width,
                                     test_weights[[f]], test_min_obs[c],
                                     test_complete_obs[d], test_na_restore[e],
                                     test_online[2]))
              
              expect_equal(roll_mean(test_data_x[[ax]], width,
                                     test_weights[[f]], test_min_obs[c],
                                     test_complete_obs[d], test_na_restore[e],
                                     test_online[1]),
                           roll_mean(test_data_x[[ax]], width,
                                     test_weights[[f]], test_min_obs[c],
                                     test_complete_obs[d], test_na_restore[e],
                                     test_online[2]))
              
              for (g in 1:length(test_center)) {
                
                expect_equal(roll_var(test_data_x[[ax]], width,
                                      test_weights[[f]], test_center[g],
                                      test_min_obs[c], test_complete_obs[d],
                                      test_na_restore[e], test_online[1]),
                             roll_var(test_data_x[[ax]], test_width[b],
                                      test_weights[[f]], test_center[g],
                                      test_min_obs[c], test_complete_obs[d],
                                      test_na_restore[e], test_online[2]))
                
                expect_equal(roll_sd(test_data_x[[ax]], width,
                                     test_weights[[f]], test_center[g],
                                     test_min_obs[c], test_complete_obs[d],
                                     test_na_restore[e], test_online[1]),
                             roll_sd(test_data_x[[ax]], test_width[b],
                                     test_weights[[f]], test_center[g],
                                     test_min_obs[c], test_complete_obs[d],
                                     test_na_restore[e], test_online[2]))
                
                for (h in 1:length(test_scale)) {
                  
                  expect_equal(roll_scale(test_data_x[[ax]], width,
                                          test_weights[[f]], test_center[g],
                                          test_scale[h], test_min_obs[c],
                                          test_complete_obs[d], test_na_restore[e],
                                          test_online[1]),
                               roll_scale(test_data_x[[ax]], width,
                                          test_weights[[f]], test_center[g],
                                          test_scale[h], test_min_obs[c],
                                          test_complete_obs[d], test_na_restore[e],
                                          test_online[2]))
                  
                  for (ay in 1:length(test_data_null)) {
                    
                    expect_equal(roll_cov(test_data_x[[ax]], test_data_null[[ay]],
                                          width, test_weights[[f]],
                                          test_center[g], test_scale[h],
                                          test_min_obs[c], test_complete_obs[d],
                                          test_na_restore[e], test_online[1]),
                                 roll_cov(test_data_x[[ax]], test_data_null[[ay]],
                                          width, test_weights[[f]],
                                          test_center[g], test_scale[h],
                                          test_min_obs[c], test_complete_obs[d],
                                          test_na_restore[e], test_online[2]))
                    
                    expect_equal(roll_cor(test_data_x[[ax]], test_data_null[[ay]],
                                          width, test_weights[[f]],
                                          test_center[g], test_scale[h],
                                          test_min_obs[c], test_complete_obs[d],
                                          test_na_restore[e], test_online[1]),
                                 roll_cor(test_data_x[[ax]], test_data_null[[ay]],
                                          width, test_weights[[f]],
                                          test_center[g], test_scale[h],
                                          test_min_obs[c], test_complete_obs[d],
                                          test_na_restore[e], test_online[2]))
                    
                  }
                  
                }
                
              }
              
              for (ay in 1:length(test_data_y)) {
                for (i in 1:length(test_intercept)) {
                  
                  # 'complete_obs = FALSE' is not supported
                  expect_equal(roll_lm(test_data_x[[ax]], test_data_y[[ay]],
                                       test_width[b], test_weights[[f]],
                                       test_intercept[i], test_min_obs[c],
                                       test_complete_obs[d], test_na_restore[e],
                                       test_online[1]),
                               roll_lm(test_data_x[[ax]], test_data_y[[ay]],
                                       test_width[b], test_weights[[f]],
                                       test_intercept[i], test_min_obs[c],
                                       test_complete_obs[d], test_na_restore[e],
                                       test_online[2]))
                  
                }
              }
              
            }
            
          }
        }
      }
      
    }
  }
  
})

test_that("equivalent to zoo::rollapply", {
  
  for (ax in 1:(length(test_data_x) - 1)) {
    for (b in 1:length(test_width)) {
      
      width <- test_width[b]
      test_weights <- list(rep(1, width), lambda ^ (width:1), 1:width,
                           rep(1, 2 * width), lambda ^ ((2 * width):1), 1:(width * 2))
      
      for (i in 1:length(test_online)) {
        
        expect_equivalent(roll_any(test_data_x[[ax]] < 0, width,
                                   test_min_obs[1], test_complete_obs[2],
                                   test_na_restore[2], test_online[i]),
                          zoo::rollapplyr(test_data_x[[ax]] < 0, width = width,
                                          any, partial = TRUE))
        
        expect_equivalent(roll_all(test_data_x[[ax]] < 0, width,
                                   test_min_obs[1], test_complete_obs[2],
                                   test_na_restore[2], test_online[i]),
                          zoo::rollapplyr(test_data_x[[ax]] < 0, width = width,
                                          all, partial = TRUE))
        
        expect_equivalent(roll_min(test_data_x[[ax]], width,
                                   test_weights[[1]], test_min_obs[1],
                                   test_complete_obs[2], test_na_restore[2],
                                   test_online[i]),
                          zoo::rollapplyr(test_data_x[[ax]], width = width,
                                          min, partial = TRUE))
        
        expect_equivalent(roll_max(test_data_x[[ax]], width,
                                   test_weights[[1]], test_min_obs[1],
                                   test_complete_obs[2], test_na_restore[2],
                                   test_online[i]),
                          zoo::rollapplyr(test_data_x[[ax]], width = width,
                                          max, partial = TRUE))
        
        expect_equivalent(roll_median(test_data_x[[ax]], width,
                                      test_weights[[1]], test_min_obs[1],
                                      test_complete_obs[2], test_na_restore[2],
                                      test_online[i]),
                          zoo::rollapplyr(test_data_x[[ax]], width = width,
                                          median, partial = TRUE))
        
        expect_equivalent(roll_sum(test_data_x[[ax]], width,
                                   test_weights[[1]], test_min_obs[1],
                                   test_complete_obs[2], test_na_restore[2],
                                   test_online[i]),
                          zoo::rollapplyr(test_data_x[[ax]], width = width,
                                          sum, partial = TRUE))
        
        expect_equivalent(roll_prod(test_data_x[[ax]], width,
                                    test_weights[[1]], test_min_obs[1],
                                    test_complete_obs[2], test_na_restore[2],
                                    test_online[i]),
                          zoo::rollapplyr(test_data_x[[ax]], width = width,
                                          prod, partial = TRUE))
        
        expect_equivalent(roll_mean(test_data_x[[ax]], width,
                                    test_weights[[1]], test_min_obs[1],
                                    test_complete_obs[2], test_na_restore[2],
                                    test_online[i]),
                          zoo::rollapplyr(test_data_x[[ax]], width = width,
                                          mean, partial = TRUE))
        
        expect_equivalent(roll_var(test_data_x[[ax]], width,
                                   test_weights[[1]], test_center[1],
                                   test_min_obs[1], test_complete_obs[2],
                                   test_na_restore[2], test_online[i]),
                          zoo::rollapplyr(test_data_x[[ax]], width = width,
                                          var, partial = TRUE))
        
        expect_equivalent(roll_sd(test_data_x[[ax]], width,
                                  test_weights[[1]], test_center[1],
                                  test_min_obs[1], test_complete_obs[2],
                                  test_na_restore[2], test_online[i]),
                          zoo::rollapplyr(test_data_x[[ax]], width = width,
                                          sd, partial = TRUE))
        
        for (g in 1:length(test_center)) {
          for (h in 1:length(test_scale)) {
            
            expect_equivalent(roll_scale(test_data_x[[ax]], width,
                                         test_weights[[1]], test_center[g],
                                         test_scale[h], test_min_obs[1],
                                         test_complete_obs[2], test_na_restore[2],
                                         test_online[i]),
                              zoo::rollapplyr(test_data_x[[ax]], width = width,
                                              scale_z, center = test_center[g],
                                              scale = test_scale[h], partial = TRUE))
            
          }
        }
        
        for (ay in 1:(length(test_data_y) - 1)) {
          
          expect_equivalent(roll_cov(test_data_x[[ax]], test_data_y[[ay]],
                                     width, test_weights[[1]],
                                     test_center[1], test_scale[2],
                                     test_min_obs[1], test_complete_obs[2],
                                     test_na_restore[2], test_online[i]),
                            rollapplyr_cube(cov, test_data_x[[ax]], test_data_y[[ay]],
                                            width))
          
          expect_equivalent(roll_cor(test_data_x[[ax]], test_data_y[[ay]],
                                     width, test_weights[[1]],
                                     test_center[1], test_scale[1],
                                     test_min_obs[1], test_complete_obs[2],
                                     test_na_restore[2], test_online[i]),
                            rollapplyr_cube(cor, test_data_x[[ax]], test_data_y[[ay]],
                                            width))
          
          for (i in 1:length(test_intercept)) {
            
            expect_equivalent(roll_lm(test_data_x[[ax]], test_data_y[[ay]][ , 1, drop = FALSE],
                                      width, test_weights[[1]],
                                      test_intercept[i], test_min_obs[1],
                                      test_complete_obs[2], test_na_restore[2],
                                      test_online[i]),
                              rollapplyr_lm(test_data_x[[ax]], test_data_y[[ay]][ , 1, drop = FALSE],
                                            width, test_intercept[i]))
            
          }
          
        }
        
      }
      
    }
  }
  
})