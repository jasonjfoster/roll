test_that("input checks", {
  
  # test data
  width <- test_width[2]
  
  for (ax in 1:length(test_check_x)) {
    
    # check_rows_equal
    expect_error(roll_cov(test_check_x[[ax]], test_check_y[-1], width),
                 "number of rows in 'x' must equal the number of rows in 'y'")
    expect_error(roll_crossprod(test_check_x[[ax]], test_check_y[-1], width),
                 "number of rows in 'x' must equal the number of rows in 'y'")
    expect_error(roll_lm(test_check_x[[ax]], test_check_y[-1], width),
                 "number of rows in 'x' must equal the number of rows in 'y'")
    
    # check_pos_int: 'width'
    expect_error(roll_any(test_check_x[[ax]], width = 0),
                 "value of 'width' must be greater than or equal to one")
    expect_error(roll_all(test_check_x[[ax]], width = 0),
                 "value of 'width' must be greater than or equal to one")
    expect_error(roll_sum(test_check_x[[ax]], width = 0),
                 "value of 'width' must be greater than or equal to one")
    expect_error(roll_prod(test_check_x[[ax]], width = 0),
                 "value of 'width' must be greater than or equal to one")
    expect_error(roll_mean(test_check_x[[ax]], width = 0),
                 "value of 'width' must be greater than or equal to one")
    expect_error(roll_idxmin(test_check_x[[ax]], width = 0),
                 "value of 'width' must be greater than or equal to one")
    expect_error(roll_quantile(test_check_x[[ax]], width = 0),
                 "value of 'width' must be greater than or equal to one")
    expect_error(roll_var(test_check_x[[ax]], width = 0),
                 "value of 'width' must be greater than or equal to one")
    expect_error(roll_sd(test_check_x[[ax]], width = 0),
                 "value of 'width' must be greater than or equal to one")
    expect_error(roll_scale(test_check_x[[ax]], width = 0),
                 "value of 'width' must be greater than or equal to one")
    expect_error(roll_cov(test_check_x[[ax]], test_check_y, width = 0),
                 "value of 'width' must be greater than or equal to one")
    expect_error(roll_crossprod(test_check_x[[ax]], test_check_y, width = 0),
                 "value of 'width' must be greater than or equal to one")
    expect_error(roll_lm(test_check_x[[ax]], test_check_y, width = 0),
                 "value of 'width' must be greater than or equal to one")
    
    # check_pos_int: 'min_obs'
    expect_error(roll_any(test_check_x[[ax]], width, min_obs = 0),
                 "value of 'min_obs' must be greater than or equal to one")
    expect_error(roll_all(test_check_x[[ax]], width, min_obs = 0),
                 "value of 'min_obs' must be greater than or equal to one")
    expect_error(roll_sum(test_check_x[[ax]], width, min_obs = 0),
                 "value of 'min_obs' must be greater than or equal to one")
    expect_error(roll_prod(test_check_x[[ax]], width, min_obs = 0),
                 "value of 'min_obs' must be greater than or equal to one")
    expect_error(roll_mean(test_check_x[[ax]], width, min_obs = 0),
                 "value of 'min_obs' must be greater than or equal to one")
    expect_error(roll_idxmin(test_check_x[[ax]], width, min_obs = 0),
                 "value of 'min_obs' must be greater than or equal to one")
    expect_error(roll_quantile(test_check_x[[ax]], width, min_obs = 0),
                 "value of 'min_obs' must be greater than or equal to one")
    expect_error(roll_var(test_check_x[[ax]], width, min_obs = 0),
                 "value of 'min_obs' must be greater than or equal to one")
    expect_error(roll_sd(test_check_x[[ax]], width, min_obs = 0),
                 "value of 'min_obs' must be greater than or equal to one")
    expect_error(roll_scale(test_check_x[[ax]], width, min_obs = 0),
                 "value of 'min_obs' must be greater than or equal to one")
    expect_error(roll_cov(test_check_x[[ax]], test_check_y, width, min_obs = 0),
                 "value of 'min_obs' must be greater than or equal to one")
    expect_error(roll_crossprod(test_check_x[[ax]], test_check_y, width, min_obs = 0),
                 "value of 'min_obs' must be greater than or equal to one")
    expect_error(roll_lm(test_check_x[[ax]], test_check_y, width, min_obs = 0),
                 "value of 'min_obs' must be greater than or equal to one")
    
    # check_weights
    expect_error(roll_sum(test_check_x[[ax]], width, rep(1, width - 1)),
                 "length of 'weights' must be greater than or equal to the number of rows in 'x' or 'width'")
    expect_error(roll_prod(test_check_x[[ax]], width, rep(1, width - 1)),
                 "length of 'weights' must be greater than or equal to the number of rows in 'x' or 'width'")
    expect_error(roll_mean(test_check_x[[ax]], width, rep(1, width - 1)),
                 "length of 'weights' must be greater than or equal to the number of rows in 'x' or 'width'")
    expect_error(roll_idxmin(test_check_x[[ax]], width, rep(1, width - 1)),
                 "length of 'weights' must be greater than or equal to the number of rows in 'x' or 'width'")
    expect_error(roll_quantile(test_check_x[[ax]], width, rep(1, width - 1)),
                 "length of 'weights' must be greater than or equal to the number of rows in 'x' or 'width'")
    expect_error(roll_var(test_check_x[[ax]], width, rep(1, width - 1)),
                 "length of 'weights' must be greater than or equal to the number of rows in 'x' or 'width'")
    expect_error(roll_sd(test_check_x[[ax]], width, rep(1, width - 1)),
                 "length of 'weights' must be greater than or equal to the number of rows in 'x' or 'width'")
    expect_error(roll_scale(test_check_x[[ax]], width, rep(1, width - 1)),
                 "length of 'weights' must be greater than or equal to the number of rows in 'x' or 'width'")
    expect_error(roll_cov(test_check_x[[ax]], test_check_y, width, rep(1, width - 1)),
                 "length of 'weights' must be greater than or equal to the number of rows in 'x' (and 'y', if applicable) or 'width'",
                 fixed = TRUE)
    expect_error(roll_crossprod(test_check_x[[ax]], test_check_y, width, rep(1, width - 1)),
                 "length of 'weights' must be greater than or equal to the number of rows in 'x' (and 'y', if applicable) or 'width'",
                 fixed = TRUE)
    expect_error(roll_lm(test_check_x[[ax]], test_check_y, width, rep(1, width - 1)),
                 "length of 'weights' must be greater than or equal to the number of rows in 'x' (and 'y') or 'width'",
                 fixed = TRUE)
    
    # check_bounds_double: 'p'
    expect_error(roll_quantile(test_check_x[[ax]], width, p = -1),
                 "value of 'p' must be between")
    
    # 'online' (roll_quantile only)
    expect_error(roll_quantile(test_check_x[[ax]], width, 1:width),
                 "'online' is only supported for equal 'weights'")
    
  }
  
})

test_that("input warnings", {
  
  # test data
  width <- test_width[2]
  
  for (ax in 1:length(test_check_x)) {
    
    # check_lambda
    expect_warning(roll_sum(test_check_x[[ax]], width, 1:width),
                   "'online' is only supported for equal or exponential decay 'weights'")
    expect_silent(roll_sum(test_check_x[[ax]], width, 1:width, online = FALSE))
    expect_warning(roll_prod(test_check_x[[ax]], width, 1:width),
                   "'online' is only supported for equal or exponential decay 'weights'")
    expect_silent(roll_prod(test_check_x[[ax]], width, 1:width, online = FALSE))
    expect_warning(roll_mean(test_check_x[[ax]], width, 1:width),
                   "'online' is only supported for equal or exponential decay 'weights'")
    expect_silent(roll_mean(test_check_x[[ax]], width, 1:width, online = FALSE))
    expect_warning(roll_var(test_check_x[[ax]], width, 1:width),
                   "'online' is only supported for equal or exponential decay 'weights'")
    expect_silent(roll_var(test_check_x[[ax]], width, 1:width, online = FALSE))
    expect_warning(roll_sd(test_check_x[[ax]], width, 1:width),
                   "'online' is only supported for equal or exponential decay 'weights'")
    expect_silent(roll_sd(test_check_x[[ax]], width, 1:width, online = FALSE))
    expect_warning(roll_scale(test_check_x[[ax]], width, 1:width),
                   "'online' is only supported for equal or exponential decay 'weights'")
    expect_silent(roll_scale(test_check_x[[ax]], width, 1:width, online = FALSE))
    expect_warning(roll_cov(test_check_x[[ax]], test_check_y, width, 1:width),
                   "'online' is only supported for equal or exponential decay 'weights'")
    expect_silent(roll_cov(test_check_x[[ax]], test_check_y, width, 1:width, online = FALSE))
    expect_warning(roll_crossprod(test_check_x[[ax]], test_check_y, width, 1:width),
                   "'online' is only supported for equal or exponential decay 'weights'")
    expect_silent(roll_crossprod(test_check_x[[ax]], test_check_y, width, 1:width, online = FALSE))
    expect_warning(roll_lm(test_check_x[[ax]], test_check_y, width, 1:width),
                   "'online' is only supported for equal or exponential decay 'weights'")
    expect_silent(roll_lm(test_check_x[[ax]], test_check_y, width, 1:width, online = FALSE))
    
    # check_lambda (not called by roll_idxmin or roll_quantile)
    expect_silent(roll_idxmin(test_check_x[[ax]], width, 1:width))
    expect_silent(roll_quantile(test_check_x[[ax]], width, 1:width, online = FALSE))
    
    # 'complete_obs = FALSE' (roll_lm only)
    expect_warning(roll_lm(test_check_x[[ax]], test_check_y, width, complete_obs = FALSE),
                   "'complete_obs = FALSE' is not supported")
    
  }
  
})
