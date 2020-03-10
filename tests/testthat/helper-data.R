set.seed(5640)
n_vars <- 5
n_obs <- 15 # 15 # 1500
n_size <- n_obs * n_vars
lambda <- 0.9 # 0.9 # 1 / 0.9
dates <- rev(seq(Sys.Date(), length.out = n_obs, by = "-1 day"))

# test data
test_ls <- list("deterministic matrix with 0's" =
                  matrix(rep(0:(n_vars - 1) / 1000, each = n_obs / n_vars),
                         nrow = n_obs, ncol = n_vars),
                
                "random matrix with 0's" =
                  matrix(sample(c(0, rnorm(n_size - 1)), n_size, replace = TRUE,
                                prob = rep(1 / n_size, n_size)),
                         nrow = n_obs, ncol = n_vars),
                
                "random matrix with 0's and NA's" =
                  matrix(sample(c(NA, 0, rnorm(n_size - 2)), n_size, replace = TRUE,
                                prob = c(1 / n_vars, rep(1 / n_size, n_size - 1))),
                         nrow = n_obs, ncol = n_vars))

if (requireNamespace("zoo", quietly = TRUE)) {
  test_ls[1:2] <- lapply(test_ls[1:2], function(x) zoo::zoo(x, dates))
}
test_ls[1:2] <- lapply(test_ls[1:2], setNames, paste0("x", rep(1:n_vars)))

# test-roll
test_roll_x <- c(lapply(test_ls, function(x){x[ , 1:3]}),
                 list("random vector with 0's and NA's" = test_ls[[3]][ , 1]))
test_roll_y <- c(lapply(test_ls, function(x){x[ , 4:5]}),
                 list("random vector with 0's and NA's" = test_ls[[3]][ , 4]))
test_roll_null <- c(test_roll_x, "null object" = list(NULL))
names(test_roll_x[[4]]) <- zoo::index(test_roll_x[[1]])
names(test_roll_y[[4]]) <- zoo::index(test_roll_y[[1]])

# test-zoo
test_zoo_x <- c(lapply(test_ls, function(x){x[ , 1:3]})[-3],
                list("random vector with 0's" = test_ls[[2]][ , 1]))
test_zoo_y <- c(lapply(test_ls, function(x){x[ , 4:5]})[-3],
                list("random vector with 0's" = test_ls[[2]][ , 4]))
test_zoo_yy <- c(lapply(test_ls, function(x){x[ , 4, drop = FALSE]})[-3], # univariate 'y' for base::lm
                 list("random vector with 0's" = test_ls[[2]][ , 4]))

# test arguments
test_width <- c(1, 5, 10, 15)
test_intercept <- c(TRUE, FALSE)
test_center <- c(TRUE, FALSE)
test_scale <- c(TRUE, FALSE)
test_min_obs <- c(1, 5, 10, 15)
test_complete_obs <- c(TRUE, FALSE)
test_na_restore <- c(TRUE, FALSE)
test_online <- c(TRUE, FALSE)