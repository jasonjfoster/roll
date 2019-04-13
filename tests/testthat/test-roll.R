context("rolling means")

n_vars_x <- 3
n_vars_y <- 1
n_obs <- 15
lambda <- 0.9

# test data
set.seed(5640)
test_data_x <- matrix(rnorm(n_obs * n_vars_x), nrow = n_obs, ncol = n_vars_x)
test_data_y <- matrix(rnorm(n_obs), nrow = n_obs, ncol = n_vars_y)
colnames(test_data_x) <- paste0("x", rep(1:n_vars_x))
colnames(test_data_y) <- paste0("y", rep(1:n_vars_y))
test_width <- c(1, 2, 10, 100)
test_online <- c(TRUE, FALSE)

last <- function(x) {
  
  return(x[length(x)])

}

rollapplyr_cube <- function(f, x, width) {
  
  n_rows_x <- nrow(x)
  n_cols_x <- ncol(x)
  r_cube <- array(NA, c(n_cols_x, n_cols_x, n_rows_x))
  
  for (i in 1:n_rows_x) {
    r_cube[ , , i] <- f(x[max(1, i - width + 1):i, , drop = FALSE])
  }
  
  return(r_cube)
  
}

rollapplyr_lm <- function(x, y, width) {
  
  n_rows_xy <- nrow(x)
  n_cols_x <- ncol(x)
  
  result <- list("coefficients" = matrix(NA, n_rows_xy, n_cols_x + 1),
                 "r.squared" = matrix(NA, n_rows_xy, 1),
                 "std.error" = matrix(NA, n_rows_xy, n_cols_x + 1))
  
  for (i in 1:n_rows_xy) {
    
    x_subset <- x[max(1, i - width + 1):i, , drop = FALSE]
    y_subset <- y[max(1, i - width + 1):i, , drop = FALSE]
    data <- as.data.frame(cbind(y_subset, x_subset))
    
    fit <- lm(reformulate(termlabels = ".", response = names(data)[1]), data = data)
    summary_fit <- summary(fit)
    
    if (nrow(coef(summary_fit)) == n_cols_x + 1) {
      
      result[["coefficients"]][i, ] <- coef(summary_fit)[ , "Estimate"]
      result[["r.squared"]][i, ] <- summary_fit$r.squared
      result[["std.error"]][i, ] <- coef(summary_fit)[ , "Std. Error"]
      
    } 
    
  }
  
  return(result)
  
}

test_that("equal to online algorithm", {
  
  
  for (b in 1:length(test_width)) {
    
    width <- test_width[b]
    test_weights <- list(rep(1, width), lambda ^ (width:1), 1:width,
                         rep(1, 2 * width), lambda ^ ((2 * width):1), 1:(width * 2))
    
    for (c in 1:length(test_weights)) {
      
      # roll functions
      expect_equal(roll_any(test_data_x < 0, width,
                            min_obs = 1, online = test_online[1]),
                   roll_any(test_data_x < 0, width,
                            min_obs = 1, online = test_online[2]))
      
      expect_equal(roll_all(test_data_x < 0, width,
                            min_obs = 1, online = test_online[1]),
                   roll_all(test_data_x < 0, width,
                            min_obs = 1, online = test_online[2]))
      
      expect_equal(roll_min(test_data_x, width, test_weights[[c]],
                            min_obs = 1, online = test_online[1]),
                   roll_min(test_data_x, width, test_weights[[c]],
                            min_obs = 1, online = test_online[2]))
      
      expect_equal(roll_max(test_data_x, width, test_weights[[c]],
                            min_obs = 1, online = test_online[1]),
                   roll_max(test_data_x, width, test_weights[[c]],
                            min_obs = 1, online = test_online[2]))
      
      expect_equal(roll_median(test_data_x, width, test_weights[[c]],
                               min_obs = 1, online = test_online[1]),
                   roll_median(test_data_x, width, test_weights[[c]],
                               min_obs = 1, online = test_online[2]))
      
      expect_equal(roll_sum(test_data_x, width, test_weights[[c]],
                            min_obs = 1, online = test_online[1]),
                   roll_sum(test_data_x, width, test_weights[[c]],
                            min_obs = 1, online = test_online[2]))
      
      expect_equal(roll_prod(test_data_x, width, test_weights[[c]],
                             min_obs = 1, online = test_online[1]),
                   roll_prod(test_data_x, width, test_weights[[c]],
                             min_obs = 1, online = test_online[2]))
      
      expect_equal(roll_mean(test_data_x, width, test_weights[[c]],
                             min_obs = 1, online = test_online[1]),
                   roll_mean(test_data_x, width, test_weights[[c]],
                             min_obs = 1, online = test_online[2]))
      
      expect_equal(roll_var(test_data_x, width, test_weights[[c]],
                            min_obs = 1, online = test_online[1]),
                   roll_var(test_data_x, width, test_weights[[c]],
                            min_obs = 1, online = test_online[2]))
      
      expect_equal(roll_sd(test_data_x, width, test_weights[[c]],
                           min_obs = 1, online = test_online[1]),
                   roll_sd(test_data_x, width, test_weights[[c]],
                           min_obs = 1, online = test_online[2]))
      
      expect_equal(roll_scale(test_data_x, width, test_weights[[c]],
                              min_obs = 1, online = test_online[1]),
                   roll_scale(test_data_x, width, test_weights[[c]],
                              min_obs = 1, online = test_online[2]))
      
      expect_equal(roll_cov(test_data_x, width = width, weights = test_weights[[c]],
                            min_obs = 1, online = test_online[1]),
                   roll_cov(test_data_x, width = width, weights = test_weights[[c]],
                            min_obs = 1, online = test_online[2]))
      
      expect_equal(roll_cor(test_data_x, width = width, weights = test_weights[[c]],
                            min_obs = 1, online = test_online[1]),
                   roll_cor(test_data_x, width = width, weights = test_weights[[c]],
                            min_obs = 1, online = test_online[2]))
      
      expect_equal(roll_lm(test_data_x, test_data_y, width = width, weights = test_weights[[c]],
                           min_obs = 1, online = test_online[1]),
                   roll_lm(test_data_x, test_data_y, width = width, weights = test_weights[[c]],
                           min_obs = 1, online = test_online[2]))
      
    }
    
  }
  
})

test_that("equivalent to zoo::rollapply", {
  
  for (b in 1:length(test_width)) {
    
    width <- test_width[b]
    test_weights <- list(rep(1, width), lambda ^ (width:1), 1:width,
                         rep(1, 2 * width), lambda ^ ((2 * width):1), 1:(width * 2))
    
    for (c in 1:length(test_weights)) {
      for (i in 1:length(test_online)) {
        
        # base functions
        expect_equivalent(roll_any(test_data_x < 0, width, min_obs = 1, online = test_online[i]),
                          zoo::rollapplyr(test_data_x < 0, width = width, any, partial = TRUE))
        
        expect_equivalent(roll_all(test_data_x < 0, width, min_obs = 1, online = test_online[i]),
                          zoo::rollapplyr(test_data_x < 0, width = width, all, partial = TRUE))
        
        expect_equivalent(roll_min(test_data_x, width, min_obs = 1, online = test_online[i]),
                          zoo::rollapplyr(test_data_x, width = width, min, partial = TRUE))
        
        expect_equivalent(roll_max(test_data_x, width, min_obs = 1, online = test_online[i]),
                          zoo::rollapplyr(test_data_x, width = width, max, partial = TRUE))
        
        expect_equivalent(roll_median(test_data_x, width, min_obs = 1, online = test_online[i]),
                          zoo::rollapplyr(test_data_x, width = width, median, partial = TRUE))
        
        expect_equivalent(roll_sum(test_data_x, width, min_obs = 1, online = test_online[i]),
                          zoo::rollapplyr(test_data_x, width = width, sum, partial = TRUE))
        
        expect_equivalent(roll_prod(test_data_x, width, min_obs = 1, online = test_online[i]),
                          zoo::rollapplyr(test_data_x, width = width, prod, partial = TRUE))
        
        expect_equivalent(roll_mean(test_data_x, width, min_obs = 1, online = test_online[i]),
                          zoo::rollapplyr(test_data_x, width = width, mean, partial = TRUE))
        
        expect_equivalent(roll_var(test_data_x, width, min_obs = 1, online = test_online[i]),
                          zoo::rollapplyr(test_data_x, width = width, var, partial = TRUE))
        
        expect_equivalent(roll_sd(test_data_x, width, min_obs = 1, online = test_online[i]),
                          zoo::rollapplyr(test_data_x, width = width, sd, partial = TRUE))
        
        expect_equivalent(roll_scale(test_data_x, width, min_obs = 1, online = test_online[i]),
                          zoo::rollapplyr(test_data_x, width = width, function(z) last(scale(z)), partial = TRUE))
        
        expect_equivalent(roll_cov(test_data_x, width = width, min_obs = 1, online = test_online[i]),
                          rollapplyr_cube(cov, test_data_x, width))
        
        expect_equivalent(roll_cor(test_data_x, width = width, min_obs = 1, online = test_online[i]),
                          rollapplyr_cube(cor, test_data_x, width))
        
        expect_equivalent(roll_lm(test_data_x, test_data_y, width = width, min_obs = 1, online = test_online[i]),
                          rollapplyr_lm(test_data_x, test_data_y, width))
        
        
      }
    }
    
  }
  
})