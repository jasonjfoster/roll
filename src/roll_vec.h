#ifndef ROLL_VEC_H
#define ROLL_VEC_H

#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

arma::ivec stl_sort_index(arma::vec& x);

// 'Worker' function for computing rolling means using an online algorithm
struct RollMeanOnlineVec {
  
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
  
  // function call operator that iterates by index
  void operator()() {
    
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

// 'Worker' function for computing rolling variances using an online algorithm
struct RollVarOnlineVec {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_var;          // destination (pass by reference)
  
  // initialize with source and destination
  RollVarOnlineVec(const NumericVector x, const int n,
                   const int n_rows_x,
                   const int width, const arma::vec arma_weights,
                   const bool center, const int min_obs,
                   const bool na_restore,
                   arma::vec& arma_var)
    : x(x), n(n),
      n_rows_x(n_rows_x),
      width(width), arma_weights(arma_weights),
      center(center), min_obs(min_obs),
      na_restore(na_restore),
      arma_var(arma_var) { }
  
  // function call operator that iterates by index
  void operator()() {
    
    int n_obs = 0;
    long double lambda = 0;
    long double w_new = 0;
    long double w_old = 0; 
    long double x_new = 0;
    long double x_old = 0;
    long double sum_w = 0;
    long double sum_x = 0;
    long double sumsq_w = 0;
    long double sumsq_x = 0;
    long double mean_prev_x = 0;
    long double mean_x = 0;
    
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
        
        sum_w = lambda * sum_w + w_new;
        sum_x = lambda * sum_x + w_new * x_new;
        sumsq_w = pow(lambda, (long double)2.0) * sumsq_w + pow(w_new, (long double)2.0);
        
        if (center && (n_obs > 0)) {
          
          // compute the mean
          mean_prev_x = mean_x;
          mean_x = sum_x / sum_w;
          
        }
        
        // compute the sum of squares
        if (!std::isnan(x[i]) && (n_obs > 1)) {
          
          sumsq_x = lambda * sumsq_x +
            w_new * (x_new - mean_x) * (x_new - mean_prev_x);
          
        } else if (std::isnan(x[i])) {
          
          sumsq_x = lambda * sumsq_x;
          
        } else if (!std::isnan(x[i]) && (n_obs == 1) && !center) {
          
          sumsq_x = w_new * pow(x_new, (long double)2.0);
          
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
        
        sum_w = lambda * sum_w + w_new - lambda * w_old;
        sum_x = lambda * sum_x + w_new * x_new - lambda * w_old * x_old;
        sumsq_w = pow(lambda, (long double)2.0) * sumsq_w +
          pow(w_new, (long double)2.0) - pow(lambda * w_old, (long double)2.0);
        
        if (center && (n_obs > 0)) {
          
          // compute the mean
          mean_prev_x = mean_x;
          mean_x = sum_x / sum_w;
          
        }
        
        // compute the sum of squares
        if (!std::isnan(x[i]) && !std::isnan(x[i - width])) {
          
          sumsq_x = lambda * sumsq_x +
            w_new * (x_new - mean_x) * (x_new - mean_prev_x) -
            lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
          
        } else if (!std::isnan(x[i]) && std::isnan(x[i - width])) {
          
          sumsq_x = lambda * sumsq_x +
            w_new * (x_new - mean_x) * (x_new - mean_prev_x);
          
        } else if (std::isnan(x[i]) && !std::isnan(x[i - width])) {
          
          sumsq_x = lambda * sumsq_x -
            lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
          
        } else if (std::isnan(x[i]) || std::isnan(x[i - width])) {
          
          sumsq_x = lambda * sumsq_x;
          
        }
        
      }
      
      // don't compute if missing value
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        // compute the unbiased estimate of variance
        if ((n_obs > 1) && (n_obs >= min_obs)) {
          
          if (std::abs(sumsq_x) <= sqrt(arma::datum::eps)) {
            arma_var[i] = 0;
          } else {
            arma_var[i] = sumsq_x / (sum_w - sumsq_w / sum_w);
          }
          
        } else {
          arma_var[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_var[i] = x[i];
        
      }
      
    }
    
  }
  
};

// 'Worker' function for computing rolling variances using a standard algorithm
struct RollVarBatchVec : public Worker {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_var;          // destination (pass by reference)
  
  // initialize with source and destination
  RollVarBatchVec(const NumericVector x, const int n,
                  const int n_rows_x,
                  const int width, const arma::vec arma_weights,
                  const bool center, const int min_obs,
                  const bool na_restore,
                  arma::vec& arma_var)
    : x(x), n(n),
      n_rows_x(n_rows_x),
      width(width), arma_weights(arma_weights),
      center(center), min_obs(min_obs),
      na_restore(na_restore),
      arma_var(arma_var) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z;
      
      long double mean_x = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        if (center) {
          
          int count = 0;
          long double sum_w = 0;
          long double sum_x = 0;
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {
            
            // don't include if missing value
            if (!std::isnan(x[i - count])) {
              
              // compute the rolling sum
              sum_w += arma_weights[n - count - 1];
              sum_x += arma_weights[n - count - 1] * x[i - count];
              
            }
            
            count += 1;
            
          }
          
          // compute the mean
          mean_x = sum_x / sum_w;
          
        }
        
        int count = 0;
        int n_obs = 0;
        long double sum_w = 0;
        long double sumsq_w = 0;
        long double sumsq_x = 0;
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value
          if (!std::isnan(x[i - count])) {
            
            sum_w += arma_weights[n - count - 1];
            sumsq_w += pow(arma_weights[n - count - 1], 2.0);
            
            // compute the rolling sum of squares with 'center' argument
            if (center) {
              sumsq_x += arma_weights[n - count - 1] *
                pow(x[i - count] - mean_x, (long double)2.0);
            } else if (!center) {
              sumsq_x += arma_weights[n - count - 1] *
                pow(x[i - count], 2.0);
            }
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        // compute the unbiased estimate of variance
        if ((n_obs > 1) && (n_obs >= min_obs)) {
          
          if (std::abs(sumsq_x) <= sqrt(arma::datum::eps)) {
            arma_var[i] = 0;
          } else {
            arma_var[i] = sumsq_x / (sum_w - sumsq_w / sum_w);
          }
          
        } else {
          arma_var[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_var[i] = x[i];
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling standard deviations using an online algorithm
struct RollSdOnlineVec {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_sd;          // destination (pass by reference)
  
  // initialize with source and destination
  RollSdOnlineVec(const NumericVector x, const int n,
                  const int n_rows_x,
                  const int width, const arma::vec arma_weights,
                  const bool center, const int min_obs,
                  const bool na_restore,
                  arma::vec& arma_sd)
    : x(x), n(n),
      n_rows_x(n_rows_x),
      width(width), arma_weights(arma_weights),
      center(center), min_obs(min_obs),
      na_restore(na_restore),
      arma_sd(arma_sd) { }
  
  // function call operator that iterates by index
  void operator()() {
    
    int n_obs = 0;
    long double lambda = 0;
    long double w_new = 0;
    long double w_old = 0; 
    long double x_new = 0;
    long double x_old = 0;
    long double sum_w = 0;
    long double sum_x = 0;
    long double sumsq_w = 0;
    long double sumsq_x = 0;
    long double mean_prev_x = 0;
    long double mean_x = 0;
    
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
        
        sum_w = lambda * sum_w + w_new;
        sum_x = lambda * sum_x + w_new * x_new;
        sumsq_w = pow(lambda, (long double)2.0) * sumsq_w + pow(w_new, (long double)2.0);
        
        if (center && (n_obs > 0)) {
          
          // compute the mean
          mean_prev_x = mean_x;
          mean_x = sum_x / sum_w;
          
        }
        
        // compute the sum of squares
        if (!std::isnan(x[i]) && (n_obs > 1)) {
          
          sumsq_x = lambda * sumsq_x +
            w_new * (x_new - mean_x) * (x_new - mean_prev_x);
          
        } else if (std::isnan(x[i])) {
          
          sumsq_x = lambda * sumsq_x;
          
        } else if (!std::isnan(x[i]) && (n_obs == 1) && !center) {
          
          sumsq_x = w_new * pow(x_new, (long double)2.0);
          
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
        
        sum_w = lambda * sum_w + w_new - lambda * w_old;
        sum_x = lambda * sum_x + w_new * x_new - lambda * w_old * x_old;
        sumsq_w = pow(lambda, (long double)2.0) * sumsq_w +
          pow(w_new, (long double)2.0) - pow(lambda * w_old, (long double)2.0);
        
        if (center && (n_obs > 0)) {
          
          // compute the mean
          mean_prev_x = mean_x;
          mean_x = sum_x / sum_w;
          
        }
        
        // compute the sum of squares
        if (!std::isnan(x[i]) && !std::isnan(x[i - width])) {
          
          sumsq_x = lambda * sumsq_x +
            w_new * (x_new - mean_x) * (x_new - mean_prev_x) -
            lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
          
        } else if (!std::isnan(x[i]) && std::isnan(x[i - width])) {
          
          sumsq_x = lambda * sumsq_x +
            w_new * (x_new - mean_x) * (x_new - mean_prev_x);
          
        } else if (std::isnan(x[i]) && !std::isnan(x[i - width])) {
          
          sumsq_x = lambda * sumsq_x -
            lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
          
        } else if (std::isnan(x[i]) || std::isnan(x[i - width])) {
          
          sumsq_x = lambda * sumsq_x;
          
        }
        
      }
      
      // don't compute if missing value
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        // compute the unbiased estimate of standard deviation
        if ((n_obs > 1) && (n_obs >= min_obs)) {
          
          if (std::abs(sumsq_x) <= sqrt(arma::datum::eps)) {
            arma_sd[i] = 0;
          } else {
            arma_sd[i] = sqrt(sumsq_x / (sum_w - sumsq_w / sum_w));
          }
          
        } else {
          arma_sd[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_sd[i] = x[i];
        
      }
      
    }
    
  }
  
};

// 'Worker' function for computing rolling standard deviations using a standard algorithm
struct RollSdBatchVec : public Worker {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_sd;          // destination (pass by reference)
  
  // initialize with source and destination
  RollSdBatchVec(const NumericVector x, const int n,
                 const int n_rows_x,
                 const int width, const arma::vec arma_weights,
                 const bool center, const int min_obs,
                 const bool na_restore,
                 arma::vec& arma_sd)
    : x(x), n(n),
      n_rows_x(n_rows_x),
      width(width), arma_weights(arma_weights),
      center(center), min_obs(min_obs),
      na_restore(na_restore),
      arma_sd(arma_sd) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z;
      
      long double mean_x = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        if (center) {
          
          int count = 0;
          long double sum_w = 0;
          long double sum_x = 0;
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {
            
            // don't include if missing value
            if (!std::isnan(x[i - count])) {
              
              // compute the rolling sum
              sum_w += arma_weights[n - count - 1];
              sum_x += arma_weights[n - count - 1] * x[i - count];
              
            }
            
            count += 1;
            
          }
          
          // compute the mean
          mean_x = sum_x / sum_w;
          
        }
        
        int count = 0;
        int n_obs = 0;
        long double sum_w = 0;
        long double sumsq_w = 0;
        long double sumsq_x = 0;
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value
          if (!std::isnan(x[i - count])) {
            
            sum_w += arma_weights[n - count - 1];
            sumsq_w += pow(arma_weights[n - count - 1], 2.0);
            
            // compute the rolling sum of squares with 'center' argument
            if (center) {
              sumsq_x += arma_weights[n - count - 1] *
                pow(x[i - count] - mean_x, (long double)2.0);
            } else if (!center) {
              sumsq_x += arma_weights[n - count - 1] *
                pow(x[i - count], 2.0);
            }
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        // compute the unbiased estimate of standard deviation
        if ((n_obs > 1) && (n_obs >= min_obs)) {
          
          if (std::abs(sumsq_x) <= sqrt(arma::datum::eps)) {
            arma_sd[i] = 0;
          } else {
            arma_sd[i] = sqrt(sumsq_x / (sum_w - sumsq_w / sum_w));
          }
          
        } else {
          arma_sd[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_sd[i] = x[i];
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling centering and scaling using an online algorithm
struct RollScaleOnlineVec {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const bool scale;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_scale;        // destination (pass by reference)
  
  // initialize with source and destination
  RollScaleOnlineVec(const NumericVector x, const int n,
                     const int n_rows_x,
                     const int width, const arma::vec arma_weights,
                     const bool center, const bool scale,
                     const int min_obs,
                     const bool na_restore, arma::vec& arma_scale)
    : x(x), n(n),
      n_rows_x(n_rows_x),
      width(width), arma_weights(arma_weights),
      center(center), scale(scale),
      min_obs(min_obs),
      na_restore(na_restore), arma_scale(arma_scale) { }
  
  // function call operator that iterates by column
  void operator()() {
    
    int n_obs = 0;
    long double lambda = 0;
    long double w_new = 0;
    long double w_old = 0; 
    long double x_new = 0;
    long double x_old = 0;
    long double sum_w = 0;
    long double sum_x = 0;
    long double sumsq_w = 0;
    long double sumsq_x = 0;
    long double mean_prev_x = 0;
    long double mean_x = 0;
    long double var_x = 0;
    long double x_ij = 0;
    
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
        x_ij = x[i];
        
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
          sumsq_w = pow(lambda, (long double)2.0) * sumsq_w + pow(w_new, (long double)2.0);
          
        } else {
          
          sum_w = w_new;
          sum_x = w_new * x_new;
          sumsq_w = pow(w_new, (long double)2.0);
          
        }
        
        if (center && (n_obs > 0)) {
          
          // compute the mean
          mean_prev_x = mean_x;
          mean_x = sum_x / sum_w;
          
        }
        
        if (scale) {
          
          // compute the sum of squares
          if (!std::isnan(x[i]) && (n_obs > 1)) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x);
            
          } else if (std::isnan(x[i])) {
            
            sumsq_x = lambda * sumsq_x;
            
          } else if (!std::isnan(x[i]) && (n_obs == 1) && !center) {
            
            sumsq_x = w_new * pow(x_new, (long double)2.0);
            
          }
          
          var_x = sumsq_x / (sum_w - sumsq_w / sum_w);
          
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
          sumsq_w = pow(lambda, (long double)2.0) * sumsq_w +
            pow(w_new, (long double)2.0) - pow(lambda * w_old, (long double)2.0);
          
        } else {
          
          sum_w = w_new;
          sum_x = w_new * x_new;
          sumsq_w = pow(w_new, (long double)2.0);
          
        }
        
        if (center && (n_obs > 0)) {
          
          // compute the mean
          mean_prev_x = mean_x;
          mean_x = sum_x / sum_w;
          
        }
        
        if (scale) {
          
          // compute the sum of squares
          if (!std::isnan(x[i]) && !std::isnan(x[i - width])) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x) -
              lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
            
          } else if (!std::isnan(x[i]) && std::isnan(x[i - width])) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x);
            
          } else if (std::isnan(x[i]) && !std::isnan(x[i - width])) {
            
            sumsq_x = lambda * sumsq_x -
              lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
            
          } else if (std::isnan(x[i]) || std::isnan(x[i - width])) {
            
            sumsq_x = lambda * sumsq_x;
            
          }
          
          var_x = sumsq_x / (sum_w - sumsq_w / sum_w);
          
        }
        
      }
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        // compute the unbiased estimate of centering and scaling
        if (n_obs >= min_obs) {
          
          if (scale && ((n_obs <= 1) || (var_x < 0) ||
              (sqrt(var_x) <= sqrt(arma::datum::eps)))) {
            arma_scale[i] = NA_REAL;
          } else if (center && scale) {
            arma_scale[i] = (x_ij - mean_x) / sqrt(var_x);
          } else if (!center && scale) {
            arma_scale[i] = x_ij / sqrt(var_x);
          } else if (center && !scale) {
            arma_scale[i] = x_ij - mean_x;
          } else if (!center && !scale) {
            arma_scale[i] = x_ij;
          }
          
        } else {
          arma_scale[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_scale[i] = x[i];
        
      }
      
    }
    
  }
  
};

// 'Worker' function for computing rolling centering and scaling using a standard algorithm
struct RollScaleBatchVec : public Worker {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_x;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const bool scale;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_scale;        // destination (pass by reference)
  
  // initialize with source and destination
  RollScaleBatchVec(const NumericVector x, const int n,
                    const int n_rows_x,
                    const int width, const arma::vec arma_weights,
                    const bool center, const bool scale,
                    const int min_obs,
                    const bool na_restore, arma::vec& arma_scale)
    : x(x), n(n),
      n_rows_x(n_rows_x),
      width(width), arma_weights(arma_weights),
      center(center), scale(scale),
      min_obs(min_obs),
      na_restore(na_restore), arma_scale(arma_scale) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z;
      
      long double mean_x = 0;
      long double var_x = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
        
        if (center) {
          
          int count = 0;
          long double sum_w = 0;
          long double sum_x = 0;
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {
            
            // don't include if missing value
            if (!std::isnan(x[i - count])) {
              
              // compute the rolling sum
              sum_w += arma_weights[n - count - 1];
              sum_x += arma_weights[n - count - 1] * x[i - count];
              
            }
            
            count += 1;
            
          }
          
          // compute the mean
          mean_x = sum_x / sum_w;
          
        }
        
        if (scale) {
          
          int count = 0;
          int n_obs = 0;
          long double sum_w = 0;
          long double sumsq_w = 0;
          long double sumsq_x = 0;
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {
            
            // don't include if missing value
            if (!std::isnan(x[i - count])) {
              
              sum_w += arma_weights[n - count - 1];
              sumsq_w += pow(arma_weights[n - count - 1], 2.0);
              
              // compute the rolling sum of squares with 'center' argument
              if (center) {
                sumsq_x += arma_weights[n - count - 1] *
                  pow(x[i - count] - mean_x, (long double)2.0);
              } else if (!center) {
                sumsq_x += arma_weights[n - count - 1] *
                  pow(x[i - count], 2.0);
              }
              
              n_obs += 1;
              
            }
            
            count += 1;
            
          }
          
          // compute the unbiased estimate of variance
          var_x = sumsq_x / (sum_w - sumsq_w / sum_w);
          
        }
        
        int count = 0;
        int n_obs = 0;
        bool any_na = false;
        long double x_ij = 0;
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value
          if (!std::isnan(x[i - count])) {
            
            // keep first non-missing value
            if (!any_na) {
              x_ij = x[i - count];
            }
            
            any_na = true;
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        // compute the unbiased estimate of centering and scaling
        if (n_obs >= min_obs) {
          
          if (scale && ((n_obs <= 1) || (var_x < 0) ||
              (sqrt(var_x) <= sqrt(arma::datum::eps)))) {
            arma_scale[i] = NA_REAL;
          } else if (center && scale) {
            arma_scale[i] = (x_ij - mean_x) / sqrt(var_x);
          } else if (!center && scale) {
            arma_scale[i] = x_ij / sqrt(var_x);
          } else if (center && !scale) {
            arma_scale[i] = x_ij - mean_x;
          } else if (!center && !scale) {
            arma_scale[i] = x_ij;
          }
          
        } else {
          arma_scale[i] = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_scale[i] = x[i];
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling covariances using an online algorithm
struct RollCovOnlineVecXX {
  
  const RVector<double> x;      // source
  const int n;
  const int n_rows_xy;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const bool scale;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_cov;          // destination (pass by reference)
  
  // initialize with source and destination
  RollCovOnlineVecXX(const NumericVector x, const int n,
                     const int n_rows_xy,
                     const int width, const arma::vec arma_weights,
                     const bool center, const bool scale,
                     const int min_obs,
                     const bool na_restore, arma::vec& arma_cov)
    : x(x), n(n),
      n_rows_xy(n_rows_xy),
      width(width), arma_weights(arma_weights),
      center(center), scale(scale),
      min_obs(min_obs),
      na_restore(na_restore), arma_cov(arma_cov) { }
  
  // function call operator that iterates by column
  void operator()() {
    
    int n_obs = 0;
    long double lambda = 0;
    long double w_new = 0;
    long double w_old = 0;      
    long double x_new = 0;
    long double x_old = 0;
    long double sum_w = 0;
    long double sum_x = 0;
    long double sumsq_w = 0;
    long double sumsq_x = 0;
    long double sumsq_xy = 0;
    long double mean_prev_x = 0;
    long double mean_x = 0;
    
    if (width > 1) {
      lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed!
    } else {
      lambda = arma_weights[n - 1];
    }
    
    for (int i = 0; i < n_rows_xy; i++) {
      
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
        
        sum_w = lambda * sum_w + w_new;
        sum_x = lambda * sum_x + w_new * x_new;
        sumsq_w = pow(lambda, (long double)2.0) * sumsq_w + pow(w_new, (long double)2.0);
        
        if (center && (n_obs > 0)) {
          
          // compute the mean
          mean_prev_x = mean_x;
          mean_x = sum_x / sum_w;
          
        }
        
        if (scale) {
          
          // compute the sum of squares
          if (!std::isnan(x[i]) && (n_obs > 1)) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x);
            
          } else if (std::isnan(x[i])) {
            
            sumsq_x = lambda * sumsq_x;
            
          } else if (!std::isnan(x[i]) && (n_obs == 1) && !center) {
            
            sumsq_x = w_new * pow(x_new, (long double)2.0);
            
          }
          
        }
        
        // compute the sum of squares
        if (!std::isnan(x[i]) && (n_obs > 1)) {

          sumsq_xy = lambda * sumsq_xy +
            w_new * (x_new - mean_x) * (x_new - mean_prev_x);

        } else if (std::isnan(x[i])) {

          sumsq_xy = lambda * sumsq_xy;

        } else if (!std::isnan(x[i]) && (n_obs == 1) && !center) {

          sumsq_xy = w_new * pow(x_new, (long double)2.0);

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
        
        sum_w = lambda * sum_w + w_new - lambda * w_old;
        sum_x = lambda * sum_x + w_new * x_new - lambda * w_old * x_old;
        sumsq_w = pow(lambda, (long double)2.0) * sumsq_w +
          pow(w_new, (long double)2.0) - pow(lambda * w_old, (long double)2.0);
        
        if (center && (n_obs > 0)) {
          
          // compute the mean
          mean_prev_x = mean_x;
          mean_x = sum_x / sum_w;
          
        }
        
        if (scale) {
          
          // compute the sum of squares
          if (!std::isnan(x[i]) && !std::isnan(x[i - width])) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x) -
              lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
            
          } else if (!std::isnan(x[i]) && std::isnan(x[i - width])) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x);
            
          } else if (std::isnan(x[i]) && !std::isnan(x[i - width])) {
            
            sumsq_x = lambda * sumsq_x -
              lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
            
          } else if (std::isnan(x[i]) || std::isnan(x[i - width])) {
            
            sumsq_x = lambda * sumsq_x;
            
          }
          
        }
        
        // compute the sum of squares
        if (!std::isnan(x[i]) && !std::isnan(x[i - width])) {

          sumsq_xy = lambda * sumsq_xy +
            w_new * (x_new - mean_x) * (x_new - mean_prev_x) -
            lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);

        } else if (!std::isnan(x[i]) && std::isnan(x[i - width])) {

          sumsq_xy = lambda * sumsq_xy +
            w_new * (x_new - mean_x) * (x_new - mean_prev_x);

        } else if (std::isnan(x[i]) && !std::isnan(x[i - width])) {

          sumsq_xy = lambda * sumsq_xy -
            lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);

        } else if (std::isnan(x[i]) || std::isnan(x[i - width])) {

          sumsq_xy = lambda * sumsq_xy;

        }
        
      }
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {
          
          // compute the unbiased estimate of variance
          if ((n_obs > 1) && (n_obs >= min_obs)) {
            
            if (scale) {
              
              // don't compute if the standard deviation is zero
              if ((sumsq_x < 0) || (sqrt(sumsq_x) <= sqrt(arma::datum::eps)))  {
                arma_cov[i] = NA_REAL;
              } else {
                
                if (std::abs(sumsq_xy) <= sqrt(arma::datum::eps)) {
                  arma_cov[i] = 0;
                } else {
                  arma_cov[i] = sumsq_xy / (sqrt(sumsq_x) * sqrt(sumsq_x));
                }
                
              }
              
            } else if (!scale) {
              
              if (std::abs(sumsq_xy) <= sqrt(arma::datum::eps)) {
                arma_cov[i] = 0;
              } else {
                arma_cov[i] = sumsq_xy / (sum_w - sumsq_w / sum_w);
              }
              
            }
            
          } else {
            arma_cov[i] = NA_REAL;
          }
          
      } else {
        
        // can be either NA or NaN
        arma_cov[i] = x[i];
        
      }
      
    }
    
  }
  
};

// 'Worker' function for computing rolling covariances using an online algorithm
struct RollCovOnlineVecXY {

  const RVector<double> x;      // source
  const RVector<double> y;      // source
  const int n;
  const int n_rows_xy;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const bool scale;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_cov;          // destination (pass by reference)
  
  // initialize with source and destination
  RollCovOnlineVecXY(const NumericVector x, const NumericVector y,
                     const int n, const int n_rows_xy,
                     const int width, const arma::vec arma_weights,
                     const bool center, const bool scale,
                     const int min_obs,
                     const bool na_restore, arma::vec& arma_cov)
    : x(x), y(y),
      n(n), n_rows_xy(n_rows_xy),
      width(width), arma_weights(arma_weights),
      center(center), scale(scale),
      min_obs(min_obs),
      na_restore(na_restore), arma_cov(arma_cov) { }
  
  // function call operator that iterates by column
  void operator()() {
    
    int n_obs = 0;
    long double lambda = 0;
    long double w_new = 0;
    long double w_old = 0;
    long double x_new = 0;
    long double x_old = 0;
    long double y_new = 0;
    long double y_old = 0;
    long double sum_w = 0;
    long double sum_x = 0;
    long double sum_y = 0;
    long double sumsq_w = 0;
    long double sumsq_x = 0;
    long double sumsq_y = 0;
    long double sumsq_xy = 0;
    long double mean_prev_x = 0;
    long double mean_prev_y = 0;
    long double mean_x = 0;
    long double mean_y = 0;
    
    if (width > 1) {
      lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed!
    } else {
      lambda = arma_weights[n - 1];
    }
    
    for (int i = 0; i < n_rows_xy; i++) {
      
      if (std::isnan(x[i]) || std::isnan(y[i])) {
        
        w_new = 0;
        x_new = 0;
        y_new = 0;
        
      } else {
        
        w_new = arma_weights[n - 1];
        x_new = x[i];
        y_new = y[i];
        
      }
      
      // expanding window
      if (i < width) {
        
        // don't include if missing value
        if (!std::isnan(x[i]) && !std::isnan(y[i])) {
          n_obs += 1;
        }
        
        sum_w = lambda * sum_w + w_new;
        sum_x = lambda * sum_x + w_new * x_new;
        sum_y = lambda * sum_y + w_new * y_new;
        sumsq_w = pow(lambda, (long double)2.0) * sumsq_w + pow(w_new, (long double)2.0);
        
        if (center && (n_obs > 0)) {
          
          // compute the mean
          mean_prev_x = mean_x;
          mean_prev_y = mean_y;
          mean_x = sum_x / sum_w;
          mean_y = sum_y / sum_w;
          
        }
        
        if (scale) {
          
          // compute the sum of squares
          if (!std::isnan(x[i]) && !std::isnan(y[i]) && (n_obs > 1)) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x);
            sumsq_y = lambda * sumsq_y +
              w_new * (y_new - mean_y) * (y_new - mean_prev_y);
            
          } else if (std::isnan(x[i]) || std::isnan(y[i])) {
            
            sumsq_x = lambda * sumsq_x;
            sumsq_y = lambda * sumsq_y;
            
          } else if (!std::isnan(x[i]) && !std::isnan(y[i]) && (n_obs == 1) && !center) {
            
            sumsq_x = w_new * pow(x_new, (long double)2.0);
            sumsq_y = w_new * pow(y_new, (long double)2.0);
            
          }
          
        }
        
        // compute the sum of squares
        if (!std::isnan(x[i]) && !std::isnan(y[i]) && (n_obs > 1)) {
          
          sumsq_xy = lambda * sumsq_xy +
            w_new * (x_new - mean_x) * (y_new - mean_prev_y);
          
        } else if (std::isnan(x[i]) || std::isnan(y[i])) {
          
          sumsq_xy = lambda * sumsq_xy;
          
        } else if (!std::isnan(x[i]) && !std::isnan(y[i]) && (n_obs == 1) && !center) {
          
          sumsq_xy = w_new * x_new * y_new;
          
        }
        
      }
      
      // rolling window
      if (i >= width) {
        
        // don't include if missing value
        if (!std::isnan(x[i]) && !std::isnan(y[i]) &&
            (std::isnan(x[i - width]) || std::isnan(y[i - width]))) {
          
          n_obs += 1;
          
        } else if ((std::isnan(x[i]) || std::isnan(y[i])) &&
          (!std::isnan(x[i - width]) && !std::isnan(y[i - width]))) {
          
          n_obs -= 1;
          
        }
        
        if (std::isnan(x[i - width]) || std::isnan(y[i - width])) {
          
          w_old = 0;
          x_old = 0;
          y_old = 0;
          
        } else {
          
          w_old = arma_weights[n - width];
          x_old = x[i - width];
          y_old = y[i - width];
          
        }
        
        sum_w = lambda * sum_w + w_new - lambda * w_old;
        sum_x = lambda * sum_x + w_new * x_new - lambda * w_old * x_old;
        sum_y = lambda * sum_y + w_new * y_new - lambda * w_old * y_old;
        sumsq_w = pow(lambda, (long double)2.0) * sumsq_w +
          pow(w_new, (long double)2.0) - pow(lambda * w_old, (long double)2.0);
        
        if (center && (n_obs > 0)) {
          
          // compute the mean
          mean_prev_x = mean_x;
          mean_prev_y = mean_y;
          mean_x = sum_x / sum_w;
          mean_y = sum_y / sum_w;
          
        }
        
        if (scale) {
          
          // compute the sum of squares
          if (!std::isnan(x[i]) && !std::isnan(y[i]) &&
              !std::isnan(x[i - width]) && !std::isnan(y[i - width])) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x) -
              lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
            sumsq_y = lambda * sumsq_y +
              w_new * (y_new - mean_y) * (y_new - mean_prev_y) -
              lambda * w_old * (y_old - mean_y) * (y_old - mean_prev_y);
            
          } else if (!std::isnan(x[i]) && !std::isnan(y[i]) &&
            (std::isnan(x[i - width]) || std::isnan(y[i - width]))) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x);
            sumsq_y = lambda * sumsq_y +
              w_new * (y_new - mean_y) * (y_new - mean_prev_y);
            
          } else if ((std::isnan(x[i]) || std::isnan(y[i])) &&
            !std::isnan(x[i - width]) && !std::isnan(y[i - width])) {
            
            sumsq_x = lambda * sumsq_x -
              lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
            sumsq_y = lambda * sumsq_y -
              lambda * w_old * (y_old - mean_y) * (y_old - mean_prev_y);
            
          } else if (std::isnan(x[i]) || std::isnan(y[i]) ||
            std::isnan(x[i - width]) || std::isnan(y[i - width])) {
            
            sumsq_x = lambda * sumsq_x;
            sumsq_y = lambda * sumsq_y;
            
          }
          
        }
        
        // compute the sum of squares
        if (!std::isnan(x[i]) && !std::isnan(y[i]) &&
            !std::isnan(x[i - width]) && !std::isnan(y[i - width])) {
          
          sumsq_xy = lambda * sumsq_xy +
            w_new * (x_new - mean_x) * (y_new - mean_prev_y) -
            lambda * w_old * (x_old - mean_x) * (y_old - mean_prev_y);
          
        } else if (!std::isnan(x[i]) && !std::isnan(y[i]) &&
          (std::isnan(x[i - width]) || std::isnan(y[i - width]))) {
          
          sumsq_xy = lambda * sumsq_xy +
            w_new * (x_new - mean_x) * (y_new - mean_prev_y);
          
        } else if ((std::isnan(x[i]) || std::isnan(y[i])) &&
          !std::isnan(x[i - width]) && !std::isnan(y[i - width])) {
          
          sumsq_xy = lambda * sumsq_xy -
            lambda * w_old * (x_old - mean_x) * (y_old - mean_prev_y);
          
        } else if (std::isnan(x[i]) || std::isnan(y[i]) ||
          std::isnan(x[i - width]) || std::isnan(y[i - width])) {
          
          sumsq_xy = lambda * sumsq_xy;
          
        }
        
      }
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]) &&
          !std::isnan(y[i]))) {
          
          // compute the unbiased estimate of variance
          if ((n_obs > 1) && (n_obs >= min_obs)) {
            
            if (scale) {
              
              // don't compute if the standard deviation is zero
              if ((sumsq_x < 0) || (sumsq_y < 0) ||
                  (sqrt(sumsq_x) <= sqrt(arma::datum::eps)) || (sqrt(sumsq_y) <= sqrt(arma::datum::eps))) {
                arma_cov[i] = NA_REAL;
              } else {
                
                if (std::abs(sumsq_xy) <= sqrt(arma::datum::eps)) {
                  arma_cov[i] = 0;
                } else {
                  arma_cov[i] = sumsq_xy / (sqrt(sumsq_x) * sqrt(sumsq_y));
                }
                
              }
              
            } else if (!scale) {
              
              if (std::abs(sumsq_xy) <= sqrt(arma::datum::eps)) {
                arma_cov[i] = 0;
              } else {
                arma_cov[i] = sumsq_xy / (sum_w - sumsq_w / sum_w);
              }
              
            }
            
          } else {
            arma_cov[i] = NA_REAL;
          }
          
      } else {
        
        // can be either NA or NaN
        if (std::isnan(x[i])) {
          arma_cov[i] = x[i];
        } else {
          arma_cov[i] = y[i];
        }
        
      }
      
    }
    
  }

};

// 'Worker' function for computing rolling covariances using a standard algorithm
struct RollCovBatchVecXX : public Worker {

  const RVector<double> x;       // source
  const int n;
  const int n_rows_xy;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const bool scale;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_cov;          // destination (pass by reference)

  // initialize with source and destination
  RollCovBatchVecXX(const NumericVector x, const int n,
                    const int n_rows_xy,
                    const int width, const arma::vec arma_weights,
                    const bool center, const bool scale,
                    const int min_obs,
                    const bool na_restore, arma::vec& arma_cov)
    : x(x), n(n),
      n_rows_xy(n_rows_xy),
      width(width), arma_weights(arma_weights),
      center(center), scale(scale),
      min_obs(min_obs),
      na_restore(na_restore), arma_cov(arma_cov) { }

  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {

      // from 1D to 3D array (lower triangle)
      int i = z;

      long double mean_x = 0;
      long double var_x = 0;

      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]))) {

          if (center) {

            int count = 0;
            long double sum_w = 0;
            long double sum_x = 0;

            // number of observations is either the window size or,
            // for partial results, the number of the current row
            while ((width > count) && (i >= count)) {

              // don't include if missing value
              if (!std::isnan(x[i - count])) {

                  // compute the rolling sum
                  sum_w += arma_weights[n - count - 1];
                sum_x += arma_weights[n - count - 1] * x[i - count];

              }

              count += 1;

            }

            // compute the mean
            mean_x = sum_x / sum_w;

          }

          if (scale) {

            int count = 0;
            long double sumsq_x = 0;

            // number of observations is either the window size or,
            // for partial results, the number of the current row
            while ((width > count) && (i >= count)) {

              // don't include if missing value
              if (!std::isnan(x[i - count])) {

                  // compute the rolling sum of squares with 'center' argument
                  if (center) {

                    sumsq_x += arma_weights[n - count - 1] *
                      pow(x[i - count] - mean_x, (long double)2.0);

                  } else if (!center) {

                    sumsq_x += arma_weights[n - count - 1] *
                      pow(x[i - count], 2.0);

                  }

              }

              count += 1;

            }

            // compute the unbiased estimate of variance
            var_x = sumsq_x;

          }

          int count = 0;
          int n_obs = 0;
          long double sum_w = 0;
          long double sumsq_w = 0;
          long double sumsq_xy = 0;

          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {

            // don't include if missing value
            if (!std::isnan(x[i - count])) {

                sum_w += arma_weights[n - count - 1];
              sumsq_w += pow(arma_weights[n - count - 1], 2.0);

              // compute the rolling sum of squares with 'center' argument
              if (center) {
                sumsq_xy += arma_weights[n - count - 1] *
                  pow(x[i - count] - mean_x, (long double)2.0);
              } else if (!center) {
                sumsq_xy += arma_weights[n - count - 1] *
                  pow(x[i - count], 2.0);
              }

              n_obs += 1;

            }

            count += 1;

          }

          // compute the unbiased estimate of covariance
          if ((n_obs > 1) && (n_obs >= min_obs)) {

            if (scale) {

              // don't compute if the standard deviation is zero
              if ((var_x < 0) || (sqrt(var_x) <= sqrt(arma::datum::eps))) {
                arma_cov[i] = NA_REAL;
              } else {

                if (std::abs(sumsq_xy) <= sqrt(arma::datum::eps)) {
                  arma_cov[i] = 0;
                } else {
                  arma_cov[i] = sumsq_xy / (sqrt(var_x) * sqrt(var_x));
                }

              }

            } else if (!scale) {

              if (std::abs(sumsq_xy) <= sqrt(arma::datum::eps)) {
                arma_cov[i] = 0;
              } else {
                arma_cov[i] = sumsq_xy / (sum_w - sumsq_w / sum_w);
              }

            }

          } else {
            arma_cov[i] = NA_REAL;
          }

      } else {

        // can be either NA or NaN
        arma_cov[i] = x[i];

      }

    }
  }

};

// 'Worker' function for computing rolling covariances using a standard algorithm
struct RollCovBatchVecXY : public Worker {

  const RVector<double> x;       // source
  const RVector<double> y;       // source
  const int n;
  const int n_rows_xy;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const bool scale;
  const int min_obs;
  const bool na_restore;
  arma::vec& arma_cov;          // destination (pass by reference)

  // initialize with source and destination
  RollCovBatchVecXY(const NumericVector x, const NumericVector y,
                    const int n, const int n_rows_xy,
                    const int width, const arma::vec arma_weights,
                    const bool center, const bool scale,
                    const int min_obs,
                    const bool na_restore, arma::vec& arma_cov)
    : x(x), y(y),
      n(n), n_rows_xy(n_rows_xy),
      width(width), arma_weights(arma_weights),
      center(center), scale(scale),
      min_obs(min_obs),
      na_restore(na_restore), arma_cov(arma_cov) { }

  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {

      // from 1D to 3D array
      int i = z;

      long double mean_x = 0;
      long double mean_y = 0;
      long double var_x = 0;
      long double var_y = 0;

      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x[i]) &&
          !std::isnan(y[i]))) {

          if (center) {

            int count = 0;
            long double sum_w = 0;
            long double sum_x = 0;
            long double sum_y = 0;

            // number of observations is either the window size or,
            // for partial results, the number of the current row
            while ((width > count) && (i >= count)) {

              // don't include if missing value
              if (!std::isnan(x[i - count]) && !std::isnan(y[i - count])) {

                  // compute the rolling sum
                  sum_w += arma_weights[n - count - 1];
                sum_x += arma_weights[n - count - 1] * x[i - count];
                sum_y += arma_weights[n - count - 1] * y[i - count];

              }

              count += 1;

            }

            // compute the mean
            mean_x = sum_x / sum_w;
            mean_y = sum_y / sum_w;

          }

          if (scale) {

            int count = 0;
            long double sumsq_x = 0;
            long double sumsq_y = 0;

            // number of observations is either the window size or,
            // for partial results, the number of the current row
            while ((width > count) && (i >= count)) {

              // don't include if missing value
              if (!std::isnan(x[i - count]) && !std::isnan(y[i - count])) {

                  // compute the rolling sum of squares with 'center' argument
                  if (center) {

                    sumsq_x += arma_weights[n - count - 1] *
                      pow(x[i - count] - mean_x, (long double)2.0);
                    sumsq_y += arma_weights[n - count - 1] *
                      pow(y[i - count] - mean_y, (long double)2.0);

                  } else if (!center) {

                    sumsq_x += arma_weights[n - count - 1] *
                      pow(x[i - count], 2.0);
                    sumsq_y += arma_weights[n - count - 1] *
                      pow(y[i - count], 2.0);

                  }

              }

              count += 1;

            }

            // compute the unbiased estimate of variance
            var_x = sumsq_x;
            var_y = sumsq_y;

          }

          int count = 0;
          int n_obs = 0;
          long double sum_w = 0;
          long double sumsq_w = 0;
          long double sumsq_xy = 0;

          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {

            // don't include if missing value
            if (!std::isnan(x[i - count]) && !std::isnan(y[i - count])) {

                sum_w += arma_weights[n - count - 1];
              sumsq_w += pow(arma_weights[n - count - 1], 2.0);

              // compute the rolling sum of squares with 'center' argument
              if (center) {
                sumsq_xy += arma_weights[n - count - 1] *
                  (x[i - count] - mean_x) * (y[i - count] - mean_y);
              } else if (!center) {
                sumsq_xy += arma_weights[n - count - 1] *
                  x[i - count] * y[i - count];
              }

              n_obs += 1;

            }

            count += 1;

          }

          // compute the unbiased estimate of covariance
          if ((n_obs > 1) && (n_obs >= min_obs)) {

            if (scale) {

              // don't compute if the standard deviation is zero
              if ((var_x < 0) || (var_y < 0) ||
                  (sqrt(var_x) <= sqrt(arma::datum::eps)) || (sqrt(var_y) <= sqrt(arma::datum::eps))) {
                arma_cov[i] = NA_REAL;
              } else {

                if (std::abs(sumsq_xy) <= sqrt(arma::datum::eps)) {
                  arma_cov[i] = 0;
                } else {
                  arma_cov[i] = sumsq_xy / (sqrt(var_x) * sqrt(var_y));
                }

              }

            } else if (!scale) {

              if (std::abs(sumsq_xy) <= sqrt(arma::datum::eps)) {
                arma_cov[i] = 0;
              } else {
                arma_cov[i] = sumsq_xy / (sum_w - sumsq_w / sum_w);
              }

            }

          } else {
            arma_cov[i] = NA_REAL;
          }

      } else {

        // can be either NA or NaN
        if (std::isnan(x[i])) {
          arma_cov[i] = x[i];
        } else {
          arma_cov[i] = y[i];
        }

      }

    }
  }

};

#endif