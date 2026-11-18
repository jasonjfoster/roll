test_that("equivalent to zoo::rollapply", {
  
  skip_on_cran()
  
  packages <- c("zoo")
  
  status <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  
  if (!all(status)) {
    
    skip(paste0("'", paste0(names(status[!status]), collapse = "', '"),
                "' package(s) required for this test"))
    
  }  
  
  # test data
  test_zoo_x <- c(lapply(test_ls[-3], function(x){x[ , 1:3]}),
                  list("random vector with 0's" = test_ls[[2]][ , 1]))
  test_zoo_y <- c(lapply(test_ls[-3], function(x){x[ , 4:5]}),
                  list("random vector with 0's" = test_ls[[2]][ , 4]))
  test_zoo_yy <- c(lapply(test_ls[-3], function(x){x[ , 4, drop = FALSE]}), # univariate 'y' for stats::lm
                   list("random vector with 0's" = test_ls[[2]][ , 4]))
  
  for (ax in 1:length(test_zoo_x)) {
    for (b in 1:length(test_width)) {
      
      width <- test_width[b]     
      test_weights <- list(rep(1, width), lambda ^ (width:1))
      
      for (i in 1:length(test_online)) {
        
        expect_equal(roll_any(test_zoo_x[[ax]] < 0, width,
                              test_min_obs[1], test_complete_obs[2],
                              test_na_restore[2], test_online[i]),
                     zoo::rollapplyr(test_zoo_x[[ax]] < 0, width = width,
                                     any, partial = TRUE))

        expect_equal(roll_all(test_zoo_x[[ax]] < 0, width,
                              test_min_obs[1], test_complete_obs[2],
                              test_na_restore[2], test_online[i]),
                     zoo::rollapplyr(test_zoo_x[[ax]] < 0, width = width,
                                     all, partial = TRUE))

        for (f in 1:length(test_weights)) {
          
          expect_equal(roll_sum(test_zoo_x[[ax]], width,
                                test_weights[[f]], test_min_obs[1],
                                test_complete_obs[2], test_na_restore[2],
                                test_online[i]),
                       zoo::rollapplyr(test_zoo_x[[ax]], width = width,
                                       sum_weights, weights = test_weights[[f]],
                                       partial = TRUE))
          
          expect_equal(roll_prod(test_zoo_x[[ax]], width,
                                 test_weights[[f]], test_min_obs[1],
                                 test_complete_obs[2], test_na_restore[2],
                                 test_online[i]),
                       zoo::rollapplyr(test_zoo_x[[ax]], width = width,
                                       prod_weights, weights = test_weights[[f]],
                                       partial = TRUE))
          
          expect_equal(roll_mean(test_zoo_x[[ax]], width,
                                 test_weights[[f]], test_min_obs[1],
                                 test_complete_obs[2], test_na_restore[2],
                                 test_online[i]),
                       zoo::rollapplyr(test_zoo_x[[ax]], width = width,
                                       mean_weights, weights = test_weights[[f]],
                                       partial = TRUE))
          
          expect_equal(roll_min(test_zoo_x[[ax]], width,
                                test_weights[[f]], test_min_obs[1],
                                test_complete_obs[2], test_na_restore[2],
                                test_online[i]),
                       zoo::rollapplyr(test_zoo_x[[ax]], width = width,
                                       min, partial = TRUE))
          
          expect_equal(roll_max(test_zoo_x[[ax]], width,
                                test_weights[[f]], test_min_obs[1],
                                test_complete_obs[2], test_na_restore[2],
                                test_online[i]),
                       zoo::rollapplyr(test_zoo_x[[ax]], width = width,
                                       max, partial = TRUE))
          
          expect_equal(roll_idxmin(test_zoo_x[[ax]], width,
                                   test_weights[[f]], test_min_obs[1],
                                   test_complete_obs[2], test_na_restore[2],
                                   test_online[i]),
                       zoo::rollapplyr(test_zoo_x[[ax]], width = width,
                                       which.min, partial = TRUE))
          
          expect_equal(roll_idxmax(test_zoo_x[[ax]], width,
                                   test_weights[[f]], test_min_obs[1],
                                   test_complete_obs[2], test_na_restore[2],
                                   test_online[i]),
                       zoo::rollapplyr(test_zoo_x[[ax]], width = width,
                                       which.max, partial = TRUE))
          
          # "'online' is only supported for equal 'weights'"
          expect_equal(roll_median(test_zoo_x[[ax]], width,
                                   rep(1, width), test_min_obs[1],
                                   test_complete_obs[2], test_na_restore[2],
                                   test_online[i]),
                       zoo::rollapplyr(test_zoo_x[[ax]], width = width,
                                       median, partial = TRUE))
          
          for (g in 1:length(test_p)) {
            
            # "'online' is only supported for equal 'weights'"
            expect_equal(roll_quantile(test_zoo_x[[ax]], width,
                                       rep(1, width), test_p[g],
                                       test_min_obs[1], test_complete_obs[2],
                                       test_na_restore[2], test_online[i]),
                         zoo::rollapplyr(test_zoo_x[[ax]], width = width,
                                         quantile, probs = test_p[g],
                                         type = 2, names = FALSE,
                                         partial = TRUE))
            
          }
          
          expect_equal(roll_var(test_zoo_x[[ax]], width,
                                test_weights[[f]], test_center[1],
                                test_min_obs[1], test_complete_obs[2],
                                test_na_restore[2], test_online[i]),
                       zoo::rollapplyr(test_zoo_x[[ax]], width = width,
                                       var_weights, weights = test_weights[[f]],
                                       partial = TRUE))
          
          expect_equal(roll_sd(test_zoo_x[[ax]], width,
                               test_weights[[f]], test_center[1],
                               test_min_obs[1], test_complete_obs[2],
                               test_na_restore[2], test_online[i]),
                       zoo::rollapplyr(test_zoo_x[[ax]], width = width,
                                       sd_weights, weights = test_weights[[f]],
                                       partial = TRUE))
          
          for (g in 1:length(test_center)) {
            for (h in 1:length(test_scale)) {
              
              expect_equal(roll_scale(test_zoo_x[[ax]], width,
                                      test_weights[[f]], test_center[g],
                                      test_scale[h], test_min_obs[1],
                                      test_complete_obs[2], test_na_restore[2],
                                      test_online[i]),
                           zoo::rollapplyr(test_zoo_x[[ax]], width = width,
                                           scale_weights, weights = test_weights[[f]],
                                           center = test_center[g], scale = test_scale[h],
                                           partial = TRUE))
              
            }
          }
          
          for (ay in 1:length(test_zoo_y)) {
            
            expect_equal(roll_cov(test_zoo_x[[ax]], test_zoo_y[[ay]],
                                  width, test_weights[[f]],
                                  test_center[1], test_scale[2],
                                  test_min_obs[1], test_complete_obs[2],
                                  test_na_restore[2], test_online[i]),
                         rollapplyr_cube(cov_weights, test_zoo_x[[ax]],
                                         test_zoo_y[[ay]], width,
                                         test_weights[[f]]))
            
            expect_equal(roll_cor(test_zoo_x[[ax]], test_zoo_y[[ay]],
                                  width, test_weights[[f]],
                                  test_center[1], test_scale[1],
                                  test_min_obs[1], test_complete_obs[2],
                                  test_na_restore[2], test_online[i]),
                         rollapplyr_cube(cor_weights, test_zoo_x[[ax]],
                                         test_zoo_y[[ay]], width,
                                         test_weights[[f]]))
            
            expect_equal(roll_crossprod(test_zoo_x[[ax]], test_zoo_y[[ay]],
                                        width, test_weights[[f]],
                                        test_center[2], test_scale[2],
                                        test_min_obs[1], test_complete_obs[2],
                                        test_na_restore[2], test_online[i]),
                         rollapplyr_cube(crossprod_weights, test_zoo_x[[ax]],
                                         test_zoo_y[[ay]], width,
                                         test_weights[[f]]))
            
          }
          
          for (ay in 1:length(test_zoo_yy)) {
            for (g in 1:length(test_intercept)) {
              
              # "essentially perfect fit: summary may be unreliable"
              # "'complete_obs = FALSE' is not supported"
              expect_equal(suppressWarnings(roll_lm(test_zoo_x[[ax]], test_zoo_yy[[ay]],
                                                    width, test_weights[[f]],
                                                    test_intercept[g], test_min_obs[1],
                                                    test_complete_obs[2], test_na_restore[2],
                                                    test_online[i])),
                           suppressWarnings(rollapplyr_lm(test_zoo_x[[ax]], test_zoo_yy[[ay]],
                                                          width, test_weights[[f]],
                                                          test_intercept[g])))
              
            }
          }
          
        }
        
      }
      
    }
  }
  
})
