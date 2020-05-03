test_that("equal to online algorithm", {
  
  # skip("long-running test")
  
  # test data
  test_roll_x <- c(lapply(test_ls, function(x){x[ , 1:3]}),
                   list("random vector with 0's and NA's" = test_ls[[3]][ , 1]))
  test_roll_y <- c(lapply(test_ls, function(x){x[ , 4:5]}),
                   list("random vector with 0's and NA's" = test_ls[[3]][ , 4]))
  test_roll_null <- c(test_roll_x, "null object" = list(NULL))
  
  if (requireNamespace("zoo", quietly = TRUE)) {
    names(test_roll_x[[4]]) <- zoo::index(test_roll_x[[1]])
    names(test_roll_y[[4]]) <- zoo::index(test_roll_y[[1]])
  }
  
  for (ax in 1:length(test_roll_x)) {
    for (b in 1:length(test_width)) {
      
      width <- test_width[b]
      test_weights <- list(lambda ^ ((2 * width):1))
      # test_weights <- list(rep(1, width), lambda ^ (width:1), 1:width,
      #                      rep(1, 2 * width), lambda ^ ((2 * width):1), 1:(width * 2))
      # test_weights <- lapply(test_weights, "-")
      # test_weights <- list(rep(0, width))
      
      for (c in 1:length(test_min_obs)) {
        for (d in 1:length(test_complete_obs)) {
          for (e in 1:length(test_na_restore)) {
            
            expect_equal(roll_any(test_roll_x[[ax]] < 0, width,
                                  test_min_obs[c], test_complete_obs[d],
                                  test_na_restore[e], test_online[1]),
                         roll_any(test_roll_x[[ax]] < 0, width,
                                  test_min_obs[c], test_complete_obs[d],
                                  test_na_restore[e], test_online[2]))
            
            expect_equal(roll_all(test_roll_x[[ax]] < 0, width,
                                  test_min_obs[c], test_complete_obs[d],
                                  test_na_restore[e], test_online[1]),
                         roll_all(test_roll_x[[ax]] < 0, width,
                                  test_min_obs[c], test_complete_obs[d],
                                  test_na_restore[e], test_online[2]))
            
            for (f in 1:length(test_weights)) {
              
              expect_equal(roll_sum(test_roll_x[[ax]], width,
                                    test_weights[[f]], test_min_obs[c],
                                    test_complete_obs[d], test_na_restore[e],
                                    test_online[1]),
                           roll_sum(test_roll_x[[ax]], width,
                                    test_weights[[f]], test_min_obs[c],
                                    test_complete_obs[d], test_na_restore[e],
                                    test_online[2]))

              expect_equal(roll_prod(test_roll_x[[ax]], width,
                                     test_weights[[f]], test_min_obs[c],
                                     test_complete_obs[d], test_na_restore[e],
                                     test_online[1]),
                           roll_prod(test_roll_x[[ax]], width,
                                     test_weights[[f]], test_min_obs[c],
                                     test_complete_obs[d], test_na_restore[e],
                                     test_online[2]))
              
              expect_equal(roll_mean(test_roll_x[[ax]], width,
                                     test_weights[[f]], test_min_obs[c],
                                     test_complete_obs[d], test_na_restore[e],
                                     test_online[1]),
                           roll_mean(test_roll_x[[ax]], width,
                                     test_weights[[f]], test_min_obs[c],
                                     test_complete_obs[d], test_na_restore[e],
                                     test_online[2]))
              
              expect_equal(roll_min(test_roll_x[[ax]], width,
                                    test_weights[[f]], test_min_obs[c],
                                    test_complete_obs[d], test_na_restore[e],
                                    test_online[1]),
                           roll_min(test_roll_x[[ax]], width,
                                    test_weights[[f]], test_min_obs[c],
                                    test_complete_obs[d], test_na_restore[e],
                                    test_online[2]))

              expect_equal(roll_max(test_roll_x[[ax]], width,
                                    test_weights[[f]], test_min_obs[c],
                                    test_complete_obs[d], test_na_restore[e],
                                    test_online[1]),
                           roll_max(test_roll_x[[ax]], width,
                                    test_weights[[f]], test_min_obs[c],
                                    test_complete_obs[d], test_na_restore[e],
                                    test_online[2]))

              expect_equal(roll_idxmin(test_roll_x[[ax]], width,
                                       test_weights[[f]], test_min_obs[c],
                                       test_complete_obs[d], test_na_restore[e],
                                       test_online[1]),
                           roll_idxmin(test_roll_x[[ax]], width,
                                       test_weights[[f]], test_min_obs[c],
                                       test_complete_obs[d], test_na_restore[e],
                                       test_online[2]))

              expect_equal(roll_idxmax(test_roll_x[[ax]], width,
                                       test_weights[[f]], test_min_obs[c],
                                       test_complete_obs[d], test_na_restore[e],
                                       test_online[1]),
                           roll_idxmax(test_roll_x[[ax]], width,
                                       test_weights[[f]], test_min_obs[c],
                                       test_complete_obs[d], test_na_restore[e],
                                       test_online[2]))

              # "'online' is not supported"
              expect_equal(roll_median(test_roll_x[[ax]], width,
                                       test_weights[[f]], test_min_obs[c],
                                       test_complete_obs[d], test_na_restore[e],
                                       test_online[1]),
                           roll_median(test_roll_x[[ax]], width,
                                       test_weights[[f]], test_min_obs[c],
                                       test_complete_obs[d], test_na_restore[e],
                                       test_online[2]))
              
              for (g in 1:length(test_center)) {
                
                expect_equal(roll_var(test_roll_x[[ax]], width,
                                      test_weights[[f]], test_center[g],
                                      test_min_obs[c], test_complete_obs[d],
                                      test_na_restore[e], test_online[1]),
                             roll_var(test_roll_x[[ax]], width,
                                      test_weights[[f]], test_center[g],
                                      test_min_obs[c], test_complete_obs[d],
                                      test_na_restore[e], test_online[2]))
                
                expect_equal(roll_sd(test_roll_x[[ax]], width,
                                     test_weights[[f]], test_center[g],
                                     test_min_obs[c], test_complete_obs[d],
                                     test_na_restore[e], test_online[1]),
                             roll_sd(test_roll_x[[ax]], width,
                                     test_weights[[f]], test_center[g],
                                     test_min_obs[c], test_complete_obs[d],
                                     test_na_restore[e], test_online[2]))
                
                for (h in 1:length(test_scale)) {
                  
                  expect_equal(roll_scale(test_roll_x[[ax]], width,
                                          test_weights[[f]], test_center[g],
                                          test_scale[h], test_min_obs[c],
                                          test_complete_obs[d], test_na_restore[e],
                                          test_online[1]),
                               roll_scale(test_roll_x[[ax]], width,
                                          test_weights[[f]], test_center[g],
                                          test_scale[h], test_min_obs[c],
                                          test_complete_obs[d], test_na_restore[e],
                                          test_online[2]))
                  
                  for (ay in 1:length(test_roll_null)) {
                    
                    expect_equal(roll_cov(test_roll_x[[ax]], test_roll_null[[ay]],
                                          width, test_weights[[f]],
                                          test_center[g], test_scale[h],
                                          test_min_obs[c], test_complete_obs[d],
                                          test_na_restore[e], test_online[1]),
                                 roll_cov(test_roll_x[[ax]], test_roll_null[[ay]],
                                          width, test_weights[[f]],
                                          test_center[g], test_scale[h],
                                          test_min_obs[c], test_complete_obs[d],
                                          test_na_restore[e], test_online[2]))
                    
                    expect_equal(roll_cor(test_roll_x[[ax]], test_roll_null[[ay]],
                                          width, test_weights[[f]],
                                          test_center[g], test_scale[h],
                                          test_min_obs[c], test_complete_obs[d],
                                          test_na_restore[e], test_online[1]),
                                 roll_cor(test_roll_x[[ax]], test_roll_null[[ay]],
                                          width, test_weights[[f]],
                                          test_center[g], test_scale[h],
                                          test_min_obs[c], test_complete_obs[d],
                                          test_na_restore[e], test_online[2]))
                    
                  }
                  
                }
                
              }
              
              for (ay in 1:length(test_roll_y)) {
                for (i in 1:length(test_intercept)) {
                  
                  # "'complete_obs = FALSE' is not supported"
                  expect_equal(roll_lm(test_roll_x[[ax]], test_roll_y[[ay]],
                                       test_width[b], test_weights[[f]],
                                       test_intercept[i], test_min_obs[c],
                                       test_complete_obs[d], test_na_restore[e],
                                       test_online[1]),
                               roll_lm(test_roll_x[[ax]], test_roll_y[[ay]],
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