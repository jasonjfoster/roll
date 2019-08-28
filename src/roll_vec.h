#ifndef ROLL_VEC_H
#define ROLL_VEC_H

#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

// 'Worker' function for computing rolling means using an online algorithm
struct RollMeanOnlineVec : public Worker {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_mean;         // destination (pass by reference)
  
  // initialize with source and destination
  RollMeanOnlineVec(const NumericVector x, const int n,
                    const int n_rows_x,
                    const int width, const arma::vec arma_weights,
                    const int min_obs,
                    const bool na_restore, arma::vec& arma_mean)
    : x(x), n(n),
      n_rows_x(n_rows_x),
      width(width), arma_weights(arma_weights),
      min_obs(min_obs),
      na_restore(na_restore), arma_mean(arma_mean) { }
  
  // function call operator that iterates by column
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      
      int n_obs = 0;
      long double lambda = 0;
      long double w_new = 0;
      long double w_old = 0;
      long double x_new = 0;
      long double x_old = 0;
      long double sum_w = 0;
      long double sum_x = 0;
      
      if (width > 1) {
        lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed!
      } else {
        lambda = arma_weights[n - 1];
      }
      
      for (int i = 0; i < n_rows_x; i++) {
        
        if (std::isnan(x[i])) {
          
          w_new = 0;
          x_new = 0;
          
        } else {
          
          w_new = arma_weights[n - 1];
          x_new = x[i];
          
        }
        
        // expanding window
        if (i < width) {
          
          // don't include if missing value
          if (!std::isnan(x[i])) {
            n_obs += 1;
          }
          
          if (width > 1) {
            
            sum_w = lambda * sum_w + w_new;
            sum_x = lambda * sum_x + w_new * x_new;
            
          } else {
            
            sum_w = w_new;
            sum_x = w_new * x_new;
            
          }
          
        }
        
        // rolling window
        if (i >= width) {
          
          // don't include if missing value
          if (!std::isnan(x[i]) && std::isnan(x[i - width])) {
            
            n_obs += 1;
            
          } else if (std::isnan(x[i]) && !std::isnan(x[i - width])) {
            
            n_obs -= 1;
            
          }
          
          if (std::isnan(x[i - width])) {
            
            w_old = 0;
            x_old = 0;
            
          } else {
            
            w_old = arma_weights[n - width];
            x_old = x[i - width];
            
          }
          
          if (width > 1) {
            
            sum_w = lambda * sum_w + w_new - lambda * w_old;
            sum_x = lambda * sum_x + w_new * x_new - lambda * w_old * x_old;
            
          } else {
            
            sum_w = w_new;
            sum_x = w_new * x_new;
            
          }
          
        }
        
        // don't compute if missing value and 'na_restore' argument is TRUE
        if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
          
          // compute the mean
          if (n_obs >= min_obs) {
            arma_mean[i] = sum_x / sum_w;
          } else {
            arma_mean[i] = NA_REAL;
          }
          
        } else {
          
          // can be either NA or NaN
          arma_mean[i] = x[i];
          
        }
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling means using a standard algorithm
struct RollMeanBatchVec : public Worker {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_mean;         // destination (pass by reference)
  
  // initialize with source and destination
  RollMeanBatchVec(const NumericVector x, const int n,
                   const int n_rows_x,
                   const int width, const arma::vec arma_weights,
                   const int min_obs,
                   const bool na_restore, arma::vec& arma_mean)
    : x(x), n(n),
      n_rows_x(n_rows_x),
      width(width), arma_weights(arma_weights),
      min_obs(min_obs),
      na_restore(na_restore), arma_mean(arma_mean) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      int i = z;
      
      int count = 0;
      int n_obs = 0;
      long double sum_w = 0;
      long double sum_x = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value
          if (!std::isnan(x[i - count])) {
            
            sum_w += arma_weights[n - count - 1];
            sum_x += arma_weights[n - count - 1] * x[i - count];
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        // compute the mean
        if (n_obs >= min_obs) {
          arma_mean[i] = sum_x / sum_w;
        } else {
          arma_mean[i] = NA_REAL;
        }
        
        
      } else {
        
        // can be either NA or NaN
        arma_mean[i] = x[i];
        
      }
      
    }
  }
  
};

#endif