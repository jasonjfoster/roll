context("rolling statistics")

test_that("equal to online algorithm", {
  
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
              
              # "'online' is not supported"
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
                  
                  # "'complete_obs' is not supported"
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
  
  if (!requireNamespace("zoo", quietly = TRUE))
    skip("zoo package required for this test")
  
  for (ax in 1:(length(test_data_x) - 1)) {
    for (b in 1:length(test_width)) {
      
      width <- test_width[b]
      test_weights <- list(rep(1, width), lambda ^ (width:1), 1:width,
                           rep(1, 2 * width), lambda ^ ((2 * width):1), 1:(width * 2))
      
      for (j in 1:length(test_online)) {
        
        expect_equivalent(roll_any(test_data_x[[ax]] < 0, width,
                                   test_min_obs[1], test_complete_obs[2],
                                   test_na_restore[2], test_online[j]),
                          zoo::rollapplyr(test_data_x[[ax]] < 0, width = width,
                                          any, partial = TRUE))
        
        expect_equivalent(roll_all(test_data_x[[ax]] < 0, width,
                                   test_min_obs[1], test_complete_obs[2],
                                   test_na_restore[2], test_online[j]),
                          zoo::rollapplyr(test_data_x[[ax]] < 0, width = width,
                                          all, partial = TRUE))
        
        expect_equivalent(roll_sum(test_data_x[[ax]], width,
                                   test_weights[[1]], test_min_obs[1],
                                   test_complete_obs[2], test_na_restore[2],
                                   test_online[j]),
                          zoo::rollapplyr(test_data_x[[ax]], width = width,
                                          sum, partial = TRUE))
        
        expect_equivalent(roll_prod(test_data_x[[ax]], width,
                                    test_weights[[1]], test_min_obs[1],
                                    test_complete_obs[2], test_na_restore[2],
                                    test_online[j]),
                          zoo::rollapplyr(test_data_x[[ax]], width = width,
                                          prod, partial = TRUE))
        
        expect_equivalent(roll_mean(test_data_x[[ax]], width,
                                    test_weights[[1]], test_min_obs[1],
                                    test_complete_obs[2], test_na_restore[2],
                                    test_online[j]),
                          zoo::rollapplyr(test_data_x[[ax]], width = width,
                                          mean, partial = TRUE))
        
        expect_equivalent(roll_min(test_data_x[[ax]], width,
                                   test_weights[[1]], test_min_obs[1],
                                   test_complete_obs[2], test_na_restore[2],
                                   test_online[j]),
                          zoo::rollapplyr(test_data_x[[ax]], width = width,
                                          min, partial = TRUE))
        
        expect_equivalent(roll_max(test_data_x[[ax]], width,
                                   test_weights[[1]], test_min_obs[1],
                                   test_complete_obs[2], test_na_restore[2],
                                   test_online[j]),
                          zoo::rollapplyr(test_data_x[[ax]], width = width,
                                          max, partial = TRUE))
        
        expect_equivalent(roll_median(test_data_x[[ax]], width,
                                      test_weights[[1]], test_min_obs[1],
                                      test_complete_obs[2], test_na_restore[2],
                                      test_online[j]),
                          zoo::rollapplyr(test_data_x[[ax]], width = width,
                                          median, partial = TRUE))
        
        expect_equivalent(roll_var(test_data_x[[ax]], width,
                                   test_weights[[1]], test_center[1],
                                   test_min_obs[1], test_complete_obs[2],
                                   test_na_restore[2], test_online[j]),
                          zoo::rollapplyr(test_data_x[[ax]], width = width,
                                          var, partial = TRUE))
        
        expect_equivalent(roll_sd(test_data_x[[ax]], width,
                                  test_weights[[1]], test_center[1],
                                  test_min_obs[1], test_complete_obs[2],
                                  test_na_restore[2], test_online[j]),
                          zoo::rollapplyr(test_data_x[[ax]], width = width,
                                          sd, partial = TRUE))
        
        for (g in 1:length(test_center)) {
          for (h in 1:length(test_scale)) {
            
            expect_equivalent(roll_scale(test_data_x[[ax]], width,
                                         test_weights[[1]], test_center[g],
                                         test_scale[h], test_min_obs[1],
                                         test_complete_obs[2], test_na_restore[2],
                                         test_online[j]),
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
                                     test_na_restore[2], test_online[j]),
                            rollapplyr_cube(cov, test_data_x[[ax]], test_data_y[[ay]],
                                            width))
          
          expect_equivalent(roll_cor(test_data_x[[ax]], test_data_y[[ay]],
                                     width, test_weights[[1]],
                                     test_center[1], test_scale[1],
                                     test_min_obs[1], test_complete_obs[2],
                                     test_na_restore[2], test_online[j]),
                            rollapplyr_cube(cor, test_data_x[[ax]], test_data_y[[ay]],
                                            width))
          
          for (i in 1:length(test_intercept)) {
            
            expect_equivalent(roll_lm(test_data_x[[ax]], test_data_y[[ay]][ , 1, drop = FALSE],
                                      width, test_weights[[1]],
                                      test_intercept[i], test_min_obs[1],
                                      test_complete_obs[2], test_na_restore[2],
                                      test_online[j]),
                              rollapplyr_lm(test_data_x[[ax]], test_data_y[[ay]][ , 1, drop = FALSE],
                                            width, test_intercept[i]))
            
          }
          
        }
        
      }
      
    }
  }
  
})