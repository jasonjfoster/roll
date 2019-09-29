#ifndef ROLL_MAT_H
#define ROLL_MAT_H

#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

arma::ivec stl_sort_min(arma::vec& x);
arma::ivec stl_sort_max(arma::vec& x);

// 'Worker' function for computing rolling any using an online algorithm
struct RollAnyOnlineMat : public Worker {
  
  const RMatrix<int> x;         // source
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const int min_obs;
  const RVector<int> rcpp_any_na;
  const bool na_restore;
  RMatrix<int> rcpp_any;        // destination (pass by reference)
  
  // initialize with source and destination
  RollAnyOnlineMat(const IntegerMatrix x, const int n_rows_x,
                   const int n_cols_x, const int width,
                   const int min_obs, const IntegerVector rcpp_any_na,
                   const bool na_restore, IntegerMatrix rcpp_any)
    : x(x), n_rows_x(n_rows_x),
      n_cols_x(n_cols_x), width(width),
      min_obs(min_obs), rcpp_any_na(rcpp_any_na),
      na_restore(na_restore), rcpp_any(rcpp_any) { }
  
  // function call operator that iterates by column
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      
      int count = 0;
      int n_obs = 0;
      int x_new = 0;
      int x_old = 0;
      int sum_x = 0;
      
      for (int i = 0; i < n_rows_x; i++) {
        
        if ((rcpp_any_na[i] != 0) || (x(i, j) == NA_INTEGER) || (x(i, j) == 0)) {
          
          x_new = 0;
          
        } else {
          
          x_new = 1;
          
        }
        
        // expanding window
        if (i < width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((rcpp_any_na[i] == 0) && (x(i, j) != NA_INTEGER)) {
            n_obs += 1;
          }
          
          sum_x = sum_x + x_new;
          
          count += 1;
          
        }
        
        // rolling window
        if (i >= width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((rcpp_any_na[i] == 0) && (x(i, j) != NA_INTEGER) &&
              ((rcpp_any_na[i - width] != 0) || (x(i - width, j) == NA_INTEGER))) {
            
            n_obs += 1;
            
          } else if (((rcpp_any_na[i] != 0) || (x(i, j) == NA_INTEGER)) &&
            (rcpp_any_na[i - width] == 0) && (x(i - width, j) != NA_INTEGER)) {
            
            n_obs -= 1;
            
          }
          
          if ((rcpp_any_na[i - width] != 0) || (x(i - width, j) == NA_INTEGER) ||
              (x(i - width, j) == 0)) {
            
            x_old = 0;
            
          } else {
            
            x_old = 1;
            
          }
          
          sum_x = sum_x + x_new - x_old;
          
        }
        
        // don't compute if missing value and 'na_restore' argument is TRUE
        if ((!na_restore) || (na_restore && (x(i, j) != NA_INTEGER))) {
          
          // compute any
          if (n_obs >= min_obs) {
            
            if (sum_x > 0) {
              rcpp_any(i, j) = 1;
            } else if (n_obs == count) {
              rcpp_any(i, j) = 0;
            } else {
              rcpp_any(i, j) = NA_INTEGER;
            }
            
          } else {
            rcpp_any(i, j) = NA_INTEGER;
          }
          
        } else {
          
          // can be either NA or NaN
          rcpp_any(i, j) = x(i, j);
          
        }
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling any using a standard algorithm
struct RollAnyBatchMat : public Worker {
  
  const RMatrix<int> x;         // source
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const int min_obs;
  const RVector<int> rcpp_any_na;
  const bool na_restore;
  RMatrix<int> rcpp_any;        // destination (pass by reference)
  
  // initialize with source and destination
  RollAnyBatchMat(const IntegerMatrix x, const int n_rows_x,
                  const int n_cols_x, const int width,
                  const int min_obs, const IntegerVector rcpp_any_na,
                  const bool na_restore, IntegerMatrix rcpp_any)
    : x(x), n_rows_x(n_rows_x),
      n_cols_x(n_cols_x), width(width),
      min_obs(min_obs), rcpp_any_na(rcpp_any_na),
      na_restore(na_restore), rcpp_any(rcpp_any) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z / n_cols_x;
      int j = z % n_cols_x;
      
      int count = 0;
      int n_obs = 0;
      int sum_x = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && (x(i, j) != NA_INTEGER))) {
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((rcpp_any_na[i - count] == 0) && (x(i - count, j) != NA_INTEGER)) {
            
            // compute the sum
            if (x(i - count, j) == 1) {
              sum_x = 1;
            }
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        // compute any
        if (n_obs >= min_obs) {
          
          if (sum_x > 0) {
            rcpp_any(i, j) = 1;
          } else if (n_obs == count) {
            rcpp_any(i, j) = 0;
          } else {
            rcpp_any(i, j) = NA_INTEGER;
          }
          
        } else {
          rcpp_any(i, j) = NA_INTEGER;
        }
        
      } else {
        
        // can be either NA or NaN
        rcpp_any(i, j) = x(i, j);
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling all using an online algorithm
struct RollAllOnlineMat : public Worker {
  
  const RMatrix<int> x;         // source
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const int min_obs;
  const RVector<int> rcpp_any_na;
  const bool na_restore;
  RMatrix<int> rcpp_all;        // destination (pass by reference)
  
  // initialize with source and destination
  RollAllOnlineMat(const IntegerMatrix x, const int n_rows_x,
                   const int n_cols_x, const int width,
                   const int min_obs, const IntegerVector rcpp_any_na,
                   const bool na_restore, IntegerMatrix rcpp_all)
    : x(x), n_rows_x(n_rows_x),
      n_cols_x(n_cols_x), width(width),
      min_obs(min_obs), rcpp_any_na(rcpp_any_na),
      na_restore(na_restore), rcpp_all(rcpp_all) { }
  
  // function call operator that iterates by column
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      
      int count = 0;
      int n_obs = 0;
      int x_new = 0;
      int x_old = 0;
      int sum_x = 0;
      
      for (int i = 0; i < n_rows_x; i++) {
        
        if ((rcpp_any_na[i] != 0) || (x(i, j) == NA_INTEGER) || (x(i, j) != 0)) {
          
          x_new = 0;
          
        } else {
          
          x_new = 1;
          
        }
        
        // expanding window
        if (i < width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((rcpp_any_na[i] == 0) && (x(i, j) != NA_INTEGER)) {
            n_obs += 1;
          }
          
          sum_x = sum_x + x_new;
          
          count += 1;
          
        }
        
        // rolling window
        if (i >= width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((rcpp_any_na[i] == 0) && (x(i, j) != NA_INTEGER) &&
              ((rcpp_any_na[i - width] != 0) || (x(i - width, j) == NA_INTEGER))) {
            
            n_obs += 1;
            
          } else if (((rcpp_any_na[i] != 0) || (x(i, j) == NA_INTEGER)) &&
            (rcpp_any_na[i - width] == 0) && (x(i - width, j) != NA_INTEGER)) {
            
            n_obs -= 1;
            
          }
          
          if ((rcpp_any_na[i - width] != 0) || (x(i - width, j) == NA_INTEGER) ||
              (x(i - width, j) != 0)) {
            
            x_old = 0;
            
          } else {
            
            x_old = 1;
            
          }
          
          sum_x = sum_x + x_new - x_old;
          
        }
        
        // don't compute if missing value and 'na_restore' argument is TRUE
        if ((!na_restore) || (na_restore && (x(i, j) != NA_INTEGER))) {
          
          // compute all
          if (n_obs >= min_obs) {
            
            if (sum_x > 0) {
              rcpp_all(i, j) = 0;
            } else if (n_obs == count) {
              rcpp_all(i, j) = 1;
            } else {
              rcpp_all(i, j) = NA_INTEGER;
            }
            
          } else {
            rcpp_all(i, j) = NA_INTEGER;
          }
          
        } else {
          
          // can be either NA or NaN
          rcpp_all(i, j) = x(i, j);
          
        }
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling all using a standard algorithm
struct RollAllBatchMat : public Worker {
  
  const RMatrix<int> x;         // source
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const int min_obs;
  const RVector<int> rcpp_any_na;
  const bool na_restore;
  RMatrix<int> rcpp_all;        // destination (pass by reference)
  
  // initialize with source and destination
  RollAllBatchMat(const IntegerMatrix x, const int n_rows_x,
                  const int n_cols_x, const int width,
                  const int min_obs, const IntegerVector rcpp_any_na,
                  const bool na_restore, IntegerMatrix rcpp_all)
    : x(x), n_rows_x(n_rows_x),
      n_cols_x(n_cols_x), width(width),
      min_obs(min_obs), rcpp_any_na(rcpp_any_na),
      na_restore(na_restore), rcpp_all(rcpp_all) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z / n_cols_x;
      int j = z % n_cols_x;
      
      int count = 0;
      int n_obs = 0;
      int sum_x = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && (x(i, j) != NA_INTEGER))) {
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((rcpp_any_na[i - count] == 0) && (x(i - count, j) != NA_INTEGER)) {
            
            // compute the sum
            if (x(i - count, j) == 0) {
              sum_x = 1;
            }
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        // compute all
        if (n_obs >= min_obs) {
          
          if (sum_x > 0) {
            rcpp_all(i, j) = 0;
          } else if (n_obs == count) {
            rcpp_all(i, j) = 1;
          } else {
            rcpp_all(i, j) = NA_INTEGER;
          }
          
        } else {
          rcpp_all(i, j) = NA_INTEGER;
        }
        
      } else {
        
        // can be either NA or NaN
        rcpp_all(i, j) = x(i, j);
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling sums using an online algorithm
struct RollSumOnlineMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_sum;          // destination (pass by reference)
  
  // initialize with source and destination
  RollSumOnlineMat(const NumericMatrix x, const int n,
                   const int n_rows_x, const int n_cols_x,
                   const int width, const arma::vec arma_weights,
                   const int min_obs, const arma::uvec arma_any_na,
                   const bool na_restore, arma::mat& arma_sum)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), arma_sum(arma_sum) { }
  
  // function call operator that iterates by column
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      
      int n_obs = 0;
      long double lambda = 0;
      long double w_new = 0;
      long double w_old = 0;
      long double x_new = 0;
      long double x_old = 0;
      long double sum_x = 0;
      
      if (width > 1) {
        lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed!
      } else {
        lambda = arma_weights[n - 1];
      }
      
      for (int i = 0; i < n_rows_x; i++) {
        
        if ((arma_any_na[i] != 0) || std::isnan(x(i, j))) {
          
          w_new = 0;
          x_new = 0;
          
        } else {
          
          w_new = arma_weights[n - 1];
          x_new = x(i, j);
          
        }
        
        // expanding window
        if (i < width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j))) {
            n_obs += 1;
          }
          
          if (width > 1) {
            sum_x = lambda * sum_x + w_new * x_new;
          } else {
            sum_x = w_new * x_new;
          }
          
        }
        
        // rolling window
        if (i >= width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) &&
              ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)))) {
            
            n_obs += 1;
            
          } else if (((arma_any_na[i] != 0) || std::isnan(x(i, j))) &&
            (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j))) {
            
            n_obs -= 1;
            
          }
          
          if ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j))) {
            
            w_old = 0;
            x_old = 0;
            
          } else {
            
            w_old = arma_weights[n - width];
            x_old = x(i - width, j);
            
          }
          
          if (width > 1) {
            sum_x = lambda * sum_x + w_new * x_new - lambda * w_old * x_old;
          } else {
            sum_x = w_new * x_new;
          }
          
        }
        
        // don't compute if missing value and 'na_restore' argument is TRUE
        if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
          
          // compute the sum
          if (n_obs >= min_obs) {
            arma_sum(i, j) = sum_x;
          } else {
            arma_sum(i, j) = NA_REAL;
          }
          
        } else {
          
          // can be either NA or NaN
          arma_sum(i, j) = x(i, j);
          
        }
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling sums using a standard algorithm
struct RollSumBatchMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_sum;          // destination (pass by reference)
  
  // initialize with source and destination
  RollSumBatchMat(const NumericMatrix x, const int n,
                  const int n_rows_x, const int n_cols_x,
                  const int width, const arma::vec arma_weights,
                  const int min_obs, const arma::uvec arma_any_na,
                  const bool na_restore, arma::mat& arma_sum)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), arma_sum(arma_sum) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z / n_cols_x;
      int j = z % n_cols_x;
      
      int count = 0;
      int n_obs = 0;
      long double sum_x = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j))) {
            
            // compute the rolling sum
            sum_x += arma_weights[n - count - 1] * x(i - count, j);
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        // compute the sum
        if (n_obs >= min_obs) {
          arma_sum(i, j) = sum_x;
        } else {
          arma_sum(i, j) = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_sum(i, j) = x(i, j);
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling products using an online algorithm
struct RollProdOnlineMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_prod;         // destination (pass by reference)
  
  // initialize with source and destination
  RollProdOnlineMat(const NumericMatrix x, const int n,
                    const int n_rows_x, const int n_cols_x,
                    const int width, const arma::vec arma_weights,
                    const int min_obs, const arma::uvec arma_any_na,
                    const bool na_restore, arma::mat& arma_prod)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), arma_prod(arma_prod) { }
  
  // function call operator that iterates by column
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      
      int n_obs = 0;
      long double lambda = 0;
      long double n_new = 0;
      long double n_old = 0;
      long double n_exp = 0;
      long double w_new = 0;
      long double w_old = 0;
      long double x_new = 0;
      long double x_old = 0;
      long double prod_w = 1;
      long double prod_x = 1;
      
      if (width > 1) {
        lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed!
      } else {
        lambda = arma_weights[n - 1];
      }
      
      for (int i = 0; i < n_rows_x; i++) {
        
        // expanding window
        if (i < width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j))) {
            n_obs += 1;
          }
          
          if ((arma_any_na[i] != 0) || std::isnan(x(i, j))) {
            
            n_new = n_obs;
            w_new = 1;
            x_new = 1;
            
          } else {
            
            n_new = n_obs - 1;
            w_new = arma_weights[n - 1];
            x_new = x(i, j);
            
          }
          
          if (n_new == 0) {
            n_exp = 1;
          } else if (n_new > n_old) {
            n_exp = n_exp * lambda;
          } else if (n_new < n_old) {
            n_exp = n_exp / lambda;
          }
          
          n_old = n_new;
          prod_w *= w_new * n_exp;
          prod_x *= x_new;
          
        }
        
        // rolling window
        if (i >= width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) &&
              ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)))) {
            
            n_obs += 1;
            
          } else if (((arma_any_na[i] != 0) || std::isnan(x(i, j))) &&
            (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j))) {
            
            n_obs -= 1;
            
          }
          
          if ((arma_any_na[i] != 0) || std::isnan(x(i, j))) {
            
            n_new = n_obs;
            w_new = 1;
            x_new = 1;
            
          } else {
            
            n_new = n_obs - 1;
            w_new = arma_weights[n - 1];
            x_new = x(i, j);
            
          }
          
          if ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j))) {
            
            w_old = 1;
            x_old = 1;
            
          } else {
            
            w_old = arma_weights[n - width];
            x_old = x(i - width, j);
            
          }
          
          if (n_new == 0) {
            n_exp = 1;
          } else if (n_new > n_old) {
            n_exp = n_exp * lambda;
          } else if (n_new < n_old) {
            n_exp = n_exp / lambda;
          }
          
          n_old = n_new;
          prod_w *= w_new * n_exp / w_old;
          prod_x *= x_new / x_old;
          
        }
        
        // don't compute if missing value and 'na_restore' argument is TRUE
        if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
          
          // compute the product
          if (n_obs >= min_obs) {
            arma_prod(i, j) = prod_w * prod_x;
          } else {
            arma_prod(i, j) = NA_REAL;
          }
          
        } else {
          
          // can be either NA or NaN
          arma_prod(i, j) = x(i, j);
          
        }
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling products using a standard algorithm
struct RollProdBatchMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_prod;         // destination (pass by reference)
  
  // initialize with source and destination
  RollProdBatchMat(const NumericMatrix x, const int n,
                   const int n_rows_x, const int n_cols_x,
                   const int width, const arma::vec arma_weights,
                   const int min_obs, const arma::uvec arma_any_na,
                   const bool na_restore, arma::mat& arma_prod)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), arma_prod(arma_prod) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z / n_cols_x;
      int j = z % n_cols_x;
      
      int count = 0;
      int n_obs = 0;
      long double prod_x = 1;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j))) {
            
            // compute the rolling product
            prod_x *= arma_weights[n - count - 1] * x(i - count, j);
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        // compute the product
        if (n_obs >= min_obs) {
          arma_prod(i, j) = prod_x;
        } else {
          arma_prod(i, j) = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_prod(i, j) = x(i, j);
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling means using an online algorithm
struct RollMeanOnlineMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_mean;         // destination (pass by reference)
  
  // initialize with source and destination
  RollMeanOnlineMat(const NumericMatrix x, const int n,
                    const int n_rows_x, const int n_cols_x,
                    const int width, const arma::vec arma_weights,
                    const int min_obs, const arma::uvec arma_any_na,
                    const bool na_restore, arma::mat& arma_mean)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      min_obs(min_obs), arma_any_na(arma_any_na),
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
        
        if ((arma_any_na[i] != 0) || std::isnan(x(i, j))) {
          
          w_new = 0;
          x_new = 0;
          
        } else {
          
          w_new = arma_weights[n - 1];
          x_new = x(i, j);
          
        }
        
        // expanding window
        if (i < width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j))) {
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
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) &&
              ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)))) {
            
            n_obs += 1;
            
          } else if (((arma_any_na[i] != 0) || std::isnan(x(i, j))) &&
            (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j))) {
            
            n_obs -= 1;
            
          }
          
          if ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j))) {
            
            w_old = 0;
            x_old = 0;
            
          } else {
            
            w_old = arma_weights[n - width];
            x_old = x(i - width, j);
            
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
        if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
          
          // compute the mean
          if (n_obs >= min_obs) {
            arma_mean(i, j) = sum_x / sum_w;
          } else {
            arma_mean(i, j) = NA_REAL;
          }
          
        } else {
          
          // can be either NA or NaN
          arma_mean(i, j) = x(i, j);
          
        }
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling means using a standard algorithm
struct RollMeanBatchMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_mean;         // destination (pass by reference)
  
  // initialize with source and destination
  RollMeanBatchMat(const NumericMatrix x, const int n,
                   const int n_rows_x, const int n_cols_x,
                   const int width, const arma::vec arma_weights,
                   const int min_obs, const arma::uvec arma_any_na,
                   const bool na_restore, arma::mat& arma_mean)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), arma_mean(arma_mean) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z / n_cols_x;
      int j = z % n_cols_x;
      
      int count = 0;
      int n_obs = 0;
      long double sum_w = 0;
      long double sum_x = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j))) {
            
            // compute the rolling sum
            sum_w += arma_weights[n - count - 1];
            sum_x += arma_weights[n - count - 1] * x(i - count, j);
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        // compute the mean
        if (n_obs >= min_obs) {
          arma_mean(i, j) = sum_x / sum_w;
        } else {
          arma_mean(i, j) = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_mean(i, j) = x(i, j);
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling minimums using an online algorithm
struct RollMinOnlineMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_min;          // destination (pass by reference)
  
  // initialize with source and destination
  RollMinOnlineMat(const NumericMatrix x, const int n,
                   const int n_rows_x, const int n_cols_x,
                   const int width, const arma::vec arma_weights,
                   const int min_obs, const arma::uvec arma_any_na,
                   const bool na_restore, arma::mat& arma_min)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), arma_min(arma_min) { }
  
  // function call operator that iterates by column
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      
      int n_obs = 0;
      long double min_x = 0;
      std::deque<int> deck(width);
      
      for (int i = 0; i < n_rows_x; i++) {
        
        // expanding window
        if (i < width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j))) {
            n_obs += 1;
          }
          
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j))) {
            
            while (!deck.empty() && ((arma_any_na[deck.back()] != 0) ||
                   std::isnan(x(deck.back(), j)) || (x(i, j) <= x(deck.back(), j)))) {
              
              deck.pop_back();
              
            }
            
            deck.push_back(i);
            
          }
          
          if (width > 1) {
            min_x = x(deck.front(), j);
          } else {
            min_x = x(i, j);
          }
          
        }
        
        // rolling window
        if (i >= width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) &&
              ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)))) {
            
            n_obs += 1;
            
          } else if (((arma_any_na[i] != 0) || std::isnan(x(i, j))) &&
            (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j))) {
            
            n_obs -= 1;
            
          }
          
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j))) {
            
            while (!deck.empty() && ((arma_any_na[deck.back()] != 0) ||
                   std::isnan(x(deck.back(), j)) || (x(i, j) <= x(deck.back(), j)))) {
              
              deck.pop_back();
              
            }
            
            deck.push_back(i);
            
          }
          
          while (!deck.empty() && (n_obs > 0) && (deck.front() <= i - width)) {
            deck.pop_front();
          }
          
          if (width > 1) {
            min_x = x(deck.front(), j);
          } else {
            min_x = x(i, j);
          }
          
        }
        
        // don't compute if missing value and 'na_restore' argument is TRUE
        if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
          
          // compute the minimum
          if (n_obs >= min_obs) {
            arma_min(i, j) = min_x;
          } else {
            arma_min(i, j) = NA_REAL;
          }
          
        } else {
          
          // can be either NA or NaN
          arma_min(i, j) = x(i, j);
          
        }
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling minimums using a standard algorithm
struct RollMinBatchMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_min;          // destination (pass by reference)
  
  // initialize with source and destination
  RollMinBatchMat(const NumericMatrix x, const int n,
                  const int n_rows_x, const int n_cols_x,
                  const int width, const arma::vec arma_weights,
                  const int min_obs, const arma::uvec arma_any_na,
                  const bool na_restore, arma::mat& arma_min)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), arma_min(arma_min) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z / n_cols_x;
      int j = z % n_cols_x;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
        
        int offset = std::max(0, i - width + 1);
        int n_size_x = i - offset + 1;
        arma::vec x_subset(n_size_x);
        arma::uvec arma_any_na_subset(n_size_x);
        
        std::copy(x.begin() + n_rows_x * j + offset, x.begin() + n_rows_x * j + i + 1,
                  x_subset.begin());
        std::copy(arma_any_na.begin() + offset, arma_any_na.begin() + i + 1,
                  arma_any_na_subset.begin());
        
        // similar to R's sort with 'index.return = TRUE'
        arma::ivec sort_ix = stl_sort_min(x_subset);
        
        int k = 0;
        int count = 0;
        int n_obs = 0;
        long double min_x = 0;
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (n_size_x - 1 >= count)) {
          
          k = sort_ix[n_size_x - count - 1];
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na_subset[k] == 0) && !std::isnan(x_subset[k])) {
            
            // last element of sorted array
            // note: 'weights' must be greater than 0
            min_x = x_subset[k];
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        // compute the minimum
        if ((n_obs >= min_obs)) {
          arma_min(i, j) = min_x;
        } else {
          arma_min(i, j) = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_min(i, j) = x(i, j);
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling maximums using an online algorithm
struct RollMaxOnlineMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_max;          // destination (pass by reference)
  
  // initialize with source and destination
  RollMaxOnlineMat(const NumericMatrix x, const int n,
                   const int n_rows_x, const int n_cols_x,
                   const int width, const arma::vec arma_weights,
                   const int min_obs, const arma::uvec arma_any_na,
                   const bool na_restore, arma::mat& arma_max)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), arma_max(arma_max) { }
  
  // function call operator that iterates by column
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      
      int n_obs = 0;
      long double max_x = 0;
      std::deque<int> deck(width);
      
      for (int i = 0; i < n_rows_x; i++) {
        
        // expanding window
        if (i < width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j))) {
            n_obs += 1;
          }
          
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j))) {
            
            while (!deck.empty() && ((arma_any_na[deck.back()] != 0) || 
                   std::isnan(x(deck.back(), j)) || (x(i, j) >= x(deck.back(), j)))) {
              
              deck.pop_back();
              
            }
            
            deck.push_back(i);
            
          }
          
          if (width > 1) {
            max_x = x(deck.front(), j);
          } else {
            max_x = x(i, j);
          }
          
        }
        
        // rolling window
        if (i >= width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) &&
              ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)))) {
            
            n_obs += 1;
            
          } else if (((arma_any_na[i] != 0) || std::isnan(x(i, j))) &&
            (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j))) {
            
            n_obs -= 1;
            
          }
          
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j))) {
            
            while (!deck.empty() && ((arma_any_na[deck.back()] != 0) ||
                   std::isnan(x(deck.back(), j)) || (x(i, j) >= x(deck.back(), j)))) {
              
              deck.pop_back();
              
            }
            
            deck.push_back(i);
            
          }
          
          while (!deck.empty() && (n_obs > 0) && (deck.front() <= i - width)) {
            deck.pop_front();
          }
          
          if (width > 1) {
            max_x = x(deck.front(), j);
          } else {
            max_x = x(i, j);
          }
          
        }
        
        // don't compute if missing value and 'na_restore' argument is TRUE
        if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
          
          // compute the maximum
          if (n_obs >= min_obs) {
            arma_max(i, j) = max_x;
          } else {
            arma_max(i, j) = NA_REAL;
          }
          
        } else {
          
          // can be either NA or NaN
          arma_max(i, j) = x(i, j);
          
        }
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling maximums using a standard algorithm
struct RollMaxBatchMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_max;          // destination (pass by reference)
  
  // initialize with source and destination
  RollMaxBatchMat(const NumericMatrix x, const int n,
                  const int n_rows_x, const int n_cols_x,
                  const int width, const arma::vec arma_weights,
                  const int min_obs, const arma::uvec arma_any_na,
                  const bool na_restore, arma::mat& arma_max)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), arma_max(arma_max) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z / n_cols_x;
      int j = z % n_cols_x;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
        
        int offset = std::max(0, i - width + 1);
        int n_size_x = i - offset + 1;
        arma::vec x_subset(n_size_x);
        arma::uvec arma_any_na_subset(n_size_x);
        
        std::copy(x.begin() + n_rows_x * j + offset, x.begin() + n_rows_x * j + i + 1,
                  x_subset.begin());
        std::copy(arma_any_na.begin() + offset, arma_any_na.begin() + i + 1,
                  arma_any_na_subset.begin());
        
        // similar to R's sort with 'index.return = TRUE'
        arma::ivec sort_ix = stl_sort_max(x_subset);
        
        int k = 0;
        int count = 0;
        int n_obs = 0;
        long double max_x = 0;
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (n_size_x - 1 >= count)) {
          
          k = sort_ix[n_size_x - count - 1];
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na_subset[k] == 0) && !std::isnan(x_subset[k])) {
            
            // first element of sorted array
            // note: 'weights' must be greater than 0
            max_x = x_subset[k];
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        // compute the maximum
        if ((n_obs >= min_obs)) {
          arma_max(i, j) = max_x;
        } else {
          arma_max(i, j) = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_max(i, j) = x(i, j);
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling index of minimums using an online algorithm
struct RollIdxMinOnlineMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const RVector<int> rcpp_any_na;
  const bool na_restore;
  RMatrix<int> rcpp_idxmin;     // destination (pass by reference)
  
  // initialize with source and destination
  RollIdxMinOnlineMat(const NumericMatrix x, const int n,
                      const int n_rows_x, const int n_cols_x,
                      const int width, const arma::vec arma_weights,
                      const int min_obs, const IntegerVector rcpp_any_na,
                      const bool na_restore, IntegerMatrix rcpp_idxmin)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      min_obs(min_obs), rcpp_any_na(rcpp_any_na),
      na_restore(na_restore), rcpp_idxmin(rcpp_idxmin) { }
  
  // function call operator that iterates by column
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      
      int n_obs = 0;
      int idxmin_x = 0;
      std::deque<int> deck(width);
      
      for (int i = 0; i < n_rows_x; i++) {
        
        // expanding window
        if (i < width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((rcpp_any_na[i] == 0) && !std::isnan(x(i, j))) {
            n_obs += 1;
          }
          
          if ((rcpp_any_na[i] == 0) && !std::isnan(x(i, j))) {
            
            while (!deck.empty() && ((rcpp_any_na[deck.back()] != 0) ||
                   std::isnan(x(deck.back(), j)) || (x(i, j) < x(deck.back(), j)))) {
              
              deck.pop_back();
              
            }
            
            deck.push_back(i);
            
          }
          
          if (width > 1) {
            idxmin_x = deck.front() + 1;
          } else {
            idxmin_x = 1;
          }
          
        }
        
        // rolling window
        if (i >= width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((rcpp_any_na[i] == 0) && !std::isnan(x(i, j)) &&
              ((rcpp_any_na[i - width] != 0) || std::isnan(x(i - width, j)))) {
            
            n_obs += 1;
            
          } else if (((rcpp_any_na[i] != 0) || std::isnan(x(i, j))) &&
            (rcpp_any_na[i - width] == 0) && !std::isnan(x(i - width, j))) {
            
            n_obs -= 1;
            
          }
          
          if ((rcpp_any_na[i] == 0) && !std::isnan(x(i, j))) {
            
            while (!deck.empty() && ((rcpp_any_na[deck.back()] != 0) ||
                   std::isnan(x(deck.back(), j)) || (x(i, j) < x(deck.back(), j)))) {
              
              deck.pop_back();
              
            }
            
            deck.push_back(i);
            
          }
          
          while (!deck.empty() && (n_obs > 0) && (deck.front() <= i - width)) {
            deck.pop_front();
          }
          
          if (width > 1) {
            idxmin_x = width - (i - deck.front());
          } else {
            idxmin_x = 1;
          }
          
        }
        
        // don't compute if missing value and 'na_restore' argument is TRUE
        if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
          
          // compute the index of minimum
          if (n_obs >= min_obs) {
            rcpp_idxmin(i, j) = idxmin_x;
          } else {
            rcpp_idxmin(i, j) = NA_REAL;
          }
          
        } else {
          
          // can be either NA or NaN
          rcpp_idxmin(i, j) = x(i, j);
          
        }
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling index of minimums using a standard algorithm
struct RollIdxMinBatchMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const RVector<int> rcpp_any_na;
  const bool na_restore;
  RMatrix<int> rcpp_idxmin;     // destination (pass by reference)
  
  // initialize with source and destination
  RollIdxMinBatchMat(const NumericMatrix x, const int n,
                     const int n_rows_x, const int n_cols_x,
                     const int width, const arma::vec arma_weights,
                     const int min_obs, const IntegerVector rcpp_any_na,
                     const bool na_restore, IntegerMatrix rcpp_idxmin)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      min_obs(min_obs), rcpp_any_na(rcpp_any_na),
      na_restore(na_restore), rcpp_idxmin(rcpp_idxmin) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z / n_cols_x;
      int j = z % n_cols_x;
      
      int count = 0;
      int n_obs = 0;
      int idxmin_x = i;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((rcpp_any_na[i - count] == 0) && !std::isnan(x(i - count, j))) {
            
            // last element of sorted array
            // note: 'weights' must be greater than 0
            if ((rcpp_any_na[idxmin_x] != 0) || std::isnan(x(idxmin_x, j)) ||
                (x(i - count, j) <= x(idxmin_x, j))) {
              
              idxmin_x = i - count;
              
            }
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        // compute the index of minimum
        if ((n_obs >= min_obs)) {
          
          if (i < width) {
            rcpp_idxmin(i, j) = idxmin_x + 1;
          } else if (i >= width) {
            rcpp_idxmin(i, j) = width - (i - idxmin_x);
          }
          
        } else {
          rcpp_idxmin(i, j) = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        rcpp_idxmin(i, j) = x(i, j);
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling index of maximums using an online algorithm
struct RollIdxMaxOnlineMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const RVector<int> rcpp_any_na;
  const bool na_restore;
  RMatrix<int> rcpp_idxmax;     // destination (pass by reference)
  
  // initialize with source and destination
  RollIdxMaxOnlineMat(const NumericMatrix x, const int n,
                      const int n_rows_x, const int n_cols_x,
                      const int width, const arma::vec arma_weights,
                      const int min_obs, const IntegerVector rcpp_any_na,
                      const bool na_restore, IntegerMatrix rcpp_idxmax)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      min_obs(min_obs), rcpp_any_na(rcpp_any_na),
      na_restore(na_restore), rcpp_idxmax(rcpp_idxmax) { }
  
  // function call operator that iterates by column
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      
      int n_obs = 0;
      int idxmax_x = 0;
      std::deque<int> deck(width);
      
      for (int i = 0; i < n_rows_x; i++) {
        
        // expanding window
        if (i < width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((rcpp_any_na[i] == 0) && !std::isnan(x(i, j))) {
            n_obs += 1;
          }
          
          if ((rcpp_any_na[i] == 0) && !std::isnan(x(i, j))) {
            
            while (!deck.empty() && ((rcpp_any_na[deck.back()] != 0) ||
                   std::isnan(x(deck.back(), j)) || (x(i, j) > x(deck.back(), j)))) {
              
              deck.pop_back();
              
            }
            
            deck.push_back(i);
            
          }
          
          if (width > 1) {
            idxmax_x = deck.front() + 1;
          } else {
            idxmax_x = 1;
          }
          
        }
        
        // rolling window
        if (i >= width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((rcpp_any_na[i] == 0) && !std::isnan(x(i, j)) &&
              ((rcpp_any_na[i - width] != 0) || std::isnan(x(i - width, j)))) {
            
            n_obs += 1;
            
          } else if (((rcpp_any_na[i] != 0) || std::isnan(x(i, j))) &&
            (rcpp_any_na[i - width] == 0) && !std::isnan(x(i - width, j))) {
            
            n_obs -= 1;
            
          }
          
          if ((rcpp_any_na[i] == 0) && !std::isnan(x(i, j))) {
            
            while (!deck.empty() && ((rcpp_any_na[deck.back()] != 0) ||
                   std::isnan(x(deck.back(), j)) || (x(i, j) > x(deck.back(), j)))) {
              
              deck.pop_back();
              
            }
            
            deck.push_back(i);
            
          }
          
          while (!deck.empty() && (n_obs > 0) && (deck.front() <= i - width)) {
            deck.pop_front();
          }
          
          if (width > 1) {
            idxmax_x = width - (i - deck.front());
          } else {
            idxmax_x = 1;
          }
          
        }
        
        // don't compute if missing value and 'na_restore' argument is TRUE
        if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
          
          // compute the index of maximum
          if (n_obs >= min_obs) {
            rcpp_idxmax(i, j) = idxmax_x;
          } else {
            rcpp_idxmax(i, j) = NA_REAL;
          }
          
        } else {
          
          // can be either NA or NaN
          rcpp_idxmax(i, j) = x(i, j);
          
        }
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling index of maximums using a standard algorithm
struct RollIdxMaxBatchMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const RVector<int> rcpp_any_na;
  const bool na_restore;
  RMatrix<int> rcpp_idxmax;     // destination (pass by reference)
  
  // initialize with source and destination
  RollIdxMaxBatchMat(const NumericMatrix x, const int n,
                     const int n_rows_x, const int n_cols_x,
                     const int width, const arma::vec arma_weights,
                     const int min_obs, const IntegerVector rcpp_any_na,
                     const bool na_restore, IntegerMatrix rcpp_idxmax)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      min_obs(min_obs), rcpp_any_na(rcpp_any_na),
      na_restore(na_restore), rcpp_idxmax(rcpp_idxmax) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z / n_cols_x;
      int j = z % n_cols_x;
      
      int count = 0;
      int n_obs = 0;
      int idxmax_x = i;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (i >= count)) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((rcpp_any_na[i - count] == 0) && !std::isnan(x(i - count, j))) {
            
            // first element of sorted array
            // note: 'weights' must be greater than 0
            if ((rcpp_any_na[idxmax_x] != 0) || std::isnan(x(idxmax_x, j)) ||
                (x(i - count, j) >= x(idxmax_x, j))) {
              
              idxmax_x = i - count;
              
            }
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        // compute the index of maximum
        if ((n_obs >= min_obs)) {
          
          if (i < width) {
            rcpp_idxmax(i, j) = idxmax_x + 1;
          } else if (i >= width) {
            rcpp_idxmax(i, j) = width - (i - idxmax_x);
          }
          
        } else {
          rcpp_idxmax(i, j) = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        rcpp_idxmax(i, j) = x(i, j);
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling medians using a standard algorithm
struct RollMedianBatchMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_median;       // destination (pass by reference)
  
  // initialize with source and destination
  RollMedianBatchMat(const NumericMatrix x, const int n,
                     const int n_rows_x, const int n_cols_x,
                     const int width, const arma::vec arma_weights,
                     const int min_obs, const arma::uvec arma_any_na,
                     const bool na_restore, arma::mat& arma_median)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), arma_median(arma_median) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z / n_cols_x;
      int j = z % n_cols_x;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
        
        int k = 0;
        int count = 0;
        long double sum_w = 0;
        
        int offset = std::max(0, i - width + 1);
        int n_size_x = i - offset + 1;
        arma::vec x_subset(n_size_x);
        arma::vec arma_weights_subset(n_size_x);
        arma::uvec arma_any_na_subset(n_size_x);
        
        std::copy(x.begin() + n_rows_x * j + offset, x.begin() + n_rows_x * j + i + 1,
                  x_subset.begin());
        std::copy(arma_weights.begin() + n - n_size_x, arma_weights.begin() + n,
                  arma_weights_subset.begin());
        std::copy(arma_any_na.begin() + offset, arma_any_na.begin() + i + 1,
                  arma_any_na_subset.begin());
        
        // similar to R's sort with 'index.return = TRUE'
        arma::ivec sort_ix = stl_sort_min(x_subset);
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (n_size_x - 1 >= count)) {
          
          k = sort_ix[n_size_x - count - 1];
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na_subset[k] == 0) && !std::isnan(x_subset[k])) {
            
            // compute the rolling sum
            sum_w += arma_weights_subset[k];
            
          }
          
          count += 1;
          
        }
        
        count = 0;
        int n_obs = 0;
        int temp_ix = 0;
        long double sum_upper_w = 0;
        long double sum_upper_x = 0;
        long double sum_upper_w_temp = 0;
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        while ((width > count) && (n_size_x - 1 >= count)) {
          
          k = sort_ix[n_size_x - count - 1];
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na_subset[k] == 0) && !std::isnan(x_subset[k])) {
            
            // last element of sorted array that is half of 'weights'
            // note: 'weights' must be greater than 0
            if (sum_upper_w / sum_w <= 0.5) {
              
              temp_ix = n_size_x - count - 1;
              sum_upper_w_temp = sum_upper_w;
              
              // compute the rolling sum
              sum_upper_w += arma_weights_subset[k];
              sum_upper_x = x_subset[k];
              
            }
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        // compute the median
        if ((n_obs >= min_obs)) {
          
          // average if upper and lower weight is equal
          if (std::abs(sum_upper_w_temp / sum_w - 0.5) <= sqrt(arma::datum::eps)) {
            
            k = sort_ix[temp_ix + 1];
            arma_median(i, j) = (x_subset[k] + sum_upper_x) / 2;
            
          } else {
            arma_median(i, j) = sum_upper_x;
          }
          
        } else {
          arma_median(i, j) = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_median(i, j) = x(i, j);
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling variances using an online algorithm
struct RollVarOnlineMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_var;          // destination (pass by reference)
  
  // initialize with source and destination
  RollVarOnlineMat(const NumericMatrix x, const int n,
                   const int n_rows_x, const int n_cols_x,
                   const int width, const arma::vec arma_weights,
                   const bool center, const int min_obs,
                   const arma::uvec arma_any_na, const bool na_restore,
                   arma::mat& arma_var)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      center(center), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      arma_var(arma_var) { }
  
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
        
        if ((arma_any_na[i] != 0) || std::isnan(x(i, j))) {
          
          w_new = 0;
          x_new = 0;
          
        } else {
          
          w_new = arma_weights[n - 1];
          x_new = x(i, j);
          
        }
        
        // expanding window
        if (i < width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j))) {
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
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && (n_obs > 1)) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x);
            
          } else if ((arma_any_na[i] != 0) || std::isnan(x(i, j))) {
            
            sumsq_x = lambda * sumsq_x;
            
          } else if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) &&
            (n_obs == 1) && !center) {
            
            sumsq_x = w_new * pow(x_new, (long double)2.0);
            
          }
          
        }
        
        // rolling window
        if (i >= width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) &&
              ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)))) {
            
            n_obs += 1;
            
          } else if (((arma_any_na[i] != 0) || std::isnan(x(i, j))) &&
            (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j))) {
            
            n_obs -= 1;
            
          }
          
          if ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j))) {
            
            w_old = 0;
            x_old = 0;
            
          } else {
            
            w_old = arma_weights[n - width];
            x_old = x(i - width, j);
            
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
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) &&
              (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j))) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x) -
              lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
            
          } else if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) &&
            ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)))) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x);
            
          } else if (((arma_any_na[i] != 0) || std::isnan(x(i, j))) &&
            (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j))) {
            
            sumsq_x = lambda * sumsq_x -
              lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
            
          } else if ((arma_any_na[i] != 0) || std::isnan(x(i, j)) ||
            (arma_any_na[i - width] != 0) || std::isnan(x(i - width, j))) {
            
            sumsq_x = lambda * sumsq_x;
            
          }
          
        }
        
        // don't compute if missing value and 'na_restore' argument is TRUE
        if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
          
          // compute the unbiased estimate of variance
          if ((n_obs > 1) && (n_obs >= min_obs)) {
            
            if (std::abs(sumsq_x) <= sqrt(arma::datum::eps)) {
              arma_var(i, j) = 0;
            } else {
              arma_var(i, j) = sumsq_x / (sum_w - sumsq_w / sum_w);
            }
            
          } else {
            arma_var(i, j) = NA_REAL;
          }
          
        } else {
          
          // can be either NA or NaN
          arma_var(i, j) = x(i, j);
          
        }
        
      }
    }
  }
  
};

// 'Worker' function for computing rolling variances using a standard algorithm
struct RollVarBatchMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_var;          // destination (pass by reference)
  
  // initialize with source and destination
  RollVarBatchMat(const NumericMatrix x, const int n,
                  const int n_rows_x, const int n_cols_x,
                  const int width, const arma::vec arma_weights,
                  const bool center, const int min_obs,
                  const arma::uvec arma_any_na, const bool na_restore,
                  arma::mat& arma_var)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      center(center), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      arma_var(arma_var) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z / n_cols_x;
      int j = z % n_cols_x;
      
      long double mean_x = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
        
        if (center) {
          
          int count = 0;
          long double sum_w = 0;
          long double sum_x = 0;
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
            if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j))) {
              
              // compute the rolling sum
              sum_w += arma_weights[n - count - 1];
              sum_x += arma_weights[n - count - 1] * x(i - count, j);
              
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
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j))) {
            
            sum_w += arma_weights[n - count - 1];
            sumsq_w += pow(arma_weights[n - count - 1], 2.0);
            
            // compute the rolling sum of squares with 'center' argument
            if (center) {
              sumsq_x += arma_weights[n - count - 1] *
                pow(x(i - count, j) - mean_x, (long double)2.0);
            } else if (!center) {
              sumsq_x += arma_weights[n - count - 1] *
                pow(x(i - count, j), 2.0);
            }
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        // compute the unbiased estimate of variance
        if ((n_obs > 1) && (n_obs >= min_obs)) {
          
          if (std::abs(sumsq_x) <= sqrt(arma::datum::eps)) {
            arma_var(i, j) = 0;
          } else {
            arma_var(i, j) = sumsq_x / (sum_w - sumsq_w / sum_w);
          }
          
        } else {
          arma_var(i, j) = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_var(i, j) = x(i, j);
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling standard deviations using an online algorithm
struct RollSdOnlineMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_sd;         // destination (pass by reference)
  
  // initialize with source and destination
  RollSdOnlineMat(const NumericMatrix x, const int n,
                  const int n_rows_x, const int n_cols_x,
                  const int width, const arma::vec arma_weights,
                  const bool center, const int min_obs,
                  const arma::uvec arma_any_na, const bool na_restore,
                  arma::mat& arma_sd)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      center(center), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      arma_sd(arma_sd) { }
  
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
        
        if ((arma_any_na[i] != 0) || std::isnan(x(i, j))) {
          
          w_new = 0;
          x_new = 0;
          
        } else {
          
          w_new = arma_weights[n - 1];
          x_new = x(i, j);
          
        }
        
        // expanding window
        if (i < width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j))) {
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
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && (n_obs > 1)) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x);
            
          } else if ((arma_any_na[i] != 0) || std::isnan(x(i, j))) {
            
            sumsq_x = lambda * sumsq_x;
            
          } else if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) &&
            (n_obs == 1) && !center) {
            
            sumsq_x = w_new * pow(x_new, (long double)2.0);
            
          }
          
        }
        
        // rolling window
        if (i >= width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) &&
              ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)))) {
            
            n_obs += 1;
            
          } else if (((arma_any_na[i] != 0) || std::isnan(x(i, j))) &&
            (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j))) {
            
            n_obs -= 1;
            
          }
          
          if ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j))) {
            
            w_old = 0;
            x_old = 0;
            
          } else {
            
            w_old = arma_weights[n - width];
            x_old = x(i - width, j);
            
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
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) &&
              (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j))) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x) -
              lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
            
          } else if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) &&
            ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)))) {
            
            sumsq_x = lambda * sumsq_x +
              w_new * (x_new - mean_x) * (x_new - mean_prev_x);
            
          } else if (((arma_any_na[i] != 0) || std::isnan(x(i, j))) &&
            (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j))) {
            
            sumsq_x = lambda * sumsq_x -
              lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
            
          } else if ((arma_any_na[i] != 0) || std::isnan(x(i, j)) ||
            (arma_any_na[i - width] != 0) || std::isnan(x(i - width, j))) {
            
            sumsq_x = lambda * sumsq_x;
            
          }
          
        }
        
        // don't compute if missing value and 'na_restore' argument is TRUE
        if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
          
          // compute the unbiased estimate of standard deviation
          if ((n_obs > 1) && (n_obs >= min_obs)) {
            
            if (std::abs(sumsq_x) <= sqrt(arma::datum::eps)) {
              arma_sd(i, j) = 0;
            } else {
              arma_sd(i, j) = sqrt(sumsq_x / (sum_w - sumsq_w / sum_w));
            }
            
          } else {
            arma_sd(i, j) = NA_REAL;
          }
          
        } else {
          
          // can be either NA or NaN
          arma_sd(i, j) = x(i, j);
          
        }
        
      }
    }
  }
  
};

// 'Worker' function for computing rolling standard deviations using a standard algorithm
struct RollSdBatchMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_sd;           // destination (pass by reference)
  
  // initialize with source and destination
  RollSdBatchMat(const NumericMatrix x, const int n,
                 const int n_rows_x, const int n_cols_x,
                 const int width, const arma::vec arma_weights,
                 const bool center, const int min_obs,
                 const arma::uvec arma_any_na, const bool na_restore,
                 arma::mat& arma_sd)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      center(center), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      arma_sd(arma_sd) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z / n_cols_x;
      int j = z % n_cols_x;
      
      long double mean_x = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
        
        if (center) {
          
          int count = 0;
          long double sum_w = 0;
          long double sum_x = 0;
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
            if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j))) {
              
              // compute the rolling sum
              sum_w += arma_weights[n - count - 1];
              sum_x += arma_weights[n - count - 1] * x(i - count, j);
              
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
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j))) {
            
            sum_w += arma_weights[n - count - 1];
            sumsq_w += pow(arma_weights[n - count - 1], 2.0);
            
            // compute the rolling sum of squares with 'center' argument
            if (center) {
              sumsq_x += arma_weights[n - count - 1] *
                pow(x(i - count, j) - mean_x, (long double)2.0);
            } else if (!center) {
              sumsq_x += arma_weights[n - count - 1] *
                pow(x(i - count, j), 2.0);
            }
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        // compute the unbiased estimate of standard deviation
        if ((n_obs > 1) && (n_obs >= min_obs)) {
          
          if (std::abs(sumsq_x) <= sqrt(arma::datum::eps)) {
            arma_sd(i, j) = 0;
          } else {
            arma_sd(i, j) = sqrt(sumsq_x / (sum_w - sumsq_w / sum_w));
          }
          
        } else {
          arma_sd(i, j) = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_sd(i, j) = x(i, j);
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling centering and scaling using an online algorithm
struct RollScaleOnlineMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const bool scale;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_scale;        // destination (pass by reference)
  
  // initialize with source and destination
  RollScaleOnlineMat(const NumericMatrix x, const int n,
                     const int n_rows_x, const int n_cols_x,
                     const int width, const arma::vec arma_weights,
                     const bool center, const bool scale,
                     const int min_obs, const arma::uvec arma_any_na,
                     const bool na_restore, arma::mat& arma_scale)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      center(center), scale(scale),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), arma_scale(arma_scale) { }
  
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
        
        if ((arma_any_na[i] != 0) || std::isnan(x(i, j))) {
          
          w_new = 0;
          x_new = 0;
          
        } else {
          
          w_new = arma_weights[n - 1];
          x_new = x(i, j);
          x_ij = x(i, j);
          
        }
        
        // expanding window
        if (i < width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j))) {
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
            if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && (n_obs > 1)) {
              
              sumsq_x = lambda * sumsq_x +
                w_new * (x_new - mean_x) * (x_new - mean_prev_x);
              
            } else if ((arma_any_na[i] != 0) || std::isnan(x(i, j))) {
              
              sumsq_x = lambda * sumsq_x;
              
            } else if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) &&
              (n_obs == 1) && !center) {
              
              sumsq_x = w_new * pow(x_new, (long double)2.0);
              
            }
            
            var_x = sumsq_x / (sum_w - sumsq_w / sum_w);
            
          }
          
        }
        
        // rolling window
        if (i >= width) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) &&
              ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)))) {
            
            n_obs += 1;
            
          } else if (((arma_any_na[i] != 0) || std::isnan(x(i, j))) &&
            (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j))) {
            
            n_obs -= 1;
            
          }
          
          if ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j))) {
            
            w_old = 0;
            x_old = 0;
            
          } else {
            
            w_old = arma_weights[n - width];
            x_old = x(i - width, j);
            
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
            if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) &&
                (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j))) {
              
              sumsq_x = lambda * sumsq_x +
                w_new * (x_new - mean_x) * (x_new - mean_prev_x) -
                lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
              
            } else if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) &&
              ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)))) {
              
              sumsq_x = lambda * sumsq_x +
                w_new * (x_new - mean_x) * (x_new - mean_prev_x);
              
            } else if (((arma_any_na[i] != 0) || std::isnan(x(i, j))) &&
              (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j))) {
              
              sumsq_x = lambda * sumsq_x -
                lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
              
            } else if ((arma_any_na[i] != 0) || std::isnan(x(i, j)) ||
              (arma_any_na[i - width] != 0) || std::isnan(x(i - width, j))) {
              
              sumsq_x = lambda * sumsq_x;
              
            }
            
            var_x = sumsq_x / (sum_w - sumsq_w / sum_w);
            
          }
          
        }
        
        // don't compute if missing value and 'na_restore' argument is TRUE
        if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
          
          // compute the unbiased estimate of centering and scaling
          if (n_obs >= min_obs) {
            
            if (scale && ((n_obs <= 1) || (var_x < 0) ||
                (sqrt(var_x) <= sqrt(arma::datum::eps)))) {
              arma_scale(i, j) = NA_REAL;
            } else if (center && scale) {
              arma_scale(i, j) = (x_ij - mean_x) / sqrt(var_x);
            } else if (!center && scale) {
              arma_scale(i, j) = x_ij / sqrt(var_x);
            } else if (center && !scale) {
              arma_scale(i, j) = x_ij - mean_x;
            } else if (!center && !scale) {
              arma_scale(i, j) = x_ij;
            }
            
          } else {
            arma_scale(i, j) = NA_REAL;
          }
          
        } else {
          
          // can be either NA or NaN
          arma_scale(i, j) = x(i, j);
          
        }
        
      }
    }
  }
  
};

// 'Worker' function for computing rolling centering and scaling using a standard algorithm
struct RollScaleBatchMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const bool scale;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_scale;        // destination (pass by reference)
  
  // initialize with source and destination
  RollScaleBatchMat(const NumericMatrix x, const int n,
                    const int n_rows_x, const int n_cols_x,
                    const int width, const arma::vec arma_weights,
                    const bool center, const bool scale,
                    const int min_obs, const arma::uvec arma_any_na,
                    const bool na_restore, arma::mat& arma_scale)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      center(center), scale(scale),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), arma_scale(arma_scale) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z / n_cols_x;
      int j = z % n_cols_x;
      
      long double mean_x = 0;
      long double var_x = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
        
        if (center) {
          
          int count = 0;
          long double sum_w = 0;
          long double sum_x = 0;
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
            if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j))) {
              
              // compute the rolling sum
              sum_w += arma_weights[n - count - 1];
              sum_x += arma_weights[n - count - 1] * x(i - count, j);
              
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
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
            if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j))) {
              
              sum_w += arma_weights[n - count - 1];
              sumsq_w += pow(arma_weights[n - count - 1], 2.0);
              
              // compute the rolling sum of squares with 'center' argument
              if (center) {
                sumsq_x += arma_weights[n - count - 1] *
                  pow(x(i - count, j) - mean_x, (long double)2.0);
              } else if (!center) {
                sumsq_x += arma_weights[n - count - 1] *
                  pow(x(i - count, j), 2.0);
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
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j))) {
            
            // keep first non-missing value
            if (!any_na) {
              x_ij = x(i - count, j);
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
            arma_scale(i, j) = NA_REAL;
          } else if (center && scale) {
            arma_scale(i, j) = (x_ij - mean_x) / sqrt(var_x);
          } else if (!center && scale) {
            arma_scale(i, j) = x_ij / sqrt(var_x);
          } else if (center && !scale) {
            arma_scale(i, j) = x_ij - mean_x;
          } else if (!center && !scale) {
            arma_scale(i, j) = x_ij;
          }
          
        } else {
          arma_scale(i, j) = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        arma_scale(i, j) = x(i, j);
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling covariances using an online algorithm
struct RollCovOnlineMatXX : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_xy;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const bool scale;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::cube& arma_cov;         // destination (pass by reference)
  
  // initialize with source and destination
  RollCovOnlineMatXX(const NumericMatrix x, const int n,
                     const int n_rows_xy, const int n_cols_x,
                     const int width, const arma::vec arma_weights,
                     const bool center, const bool scale,
                     const int min_obs, const arma::uvec arma_any_na,
                     const bool na_restore, arma::cube& arma_cov)
    : x(x), n(n),
      n_rows_xy(n_rows_xy), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      center(center), scale(scale),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), arma_cov(arma_cov) { }
  
  // function call operator that iterates by column
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      for (std::size_t k = 0; k <= j; k++) {
        
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
          
          if ((arma_any_na[i] != 0) || std::isnan(x(i, j)) || std::isnan(x(i, k))) {
            
            w_new = 0;
            x_new = 0;
            y_new = 0;
            
          } else {
            
            w_new = arma_weights[n - 1];
            x_new = x(i, j);
            y_new = x(i, k);
            
          }
          
          // expanding window
          if (i < width) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
            if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(x(i, k))) {
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
              if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(x(i, k)) &&
                  (n_obs > 1)) {
                
                sumsq_x = lambda * sumsq_x +
                  w_new * (x_new - mean_x) * (x_new - mean_prev_x);
                sumsq_y = lambda * sumsq_y +
                  w_new * (y_new - mean_y) * (y_new - mean_prev_y);
                
              } else if ((arma_any_na[i] != 0) || std::isnan(x(i, j)) || std::isnan(x(i, k))) {
                
                sumsq_x = lambda * sumsq_x;
                sumsq_y = lambda * sumsq_y;
                
              } else if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(x(i, k)) &&
                (n_obs == 1) && !center) {
                
                sumsq_x = w_new * pow(x_new, (long double)2.0);
                sumsq_y = w_new * pow(y_new, (long double)2.0);
                
              }
              
            }
            
            // compute the sum of squares
            if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(x(i, k)) &&
                (n_obs > 1)) {
              
              sumsq_xy = lambda * sumsq_xy +
                w_new * (x_new - mean_x) * (y_new - mean_prev_y);
              
            } else if ((arma_any_na[i] != 0) || std::isnan(x(i, j)) || std::isnan(x(i, k))) {
              
              sumsq_xy = lambda * sumsq_xy;
              
            } else if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(x(i, k)) &&
              (n_obs == 1) && !center) {
              
              sumsq_xy = w_new * x_new * y_new;
              
            }
            
          }
          
          // rolling window
          if (i >= width) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
            if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(x(i, k)) &&
                ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)) || std::isnan(x(i - width, k)))) {
              
              n_obs += 1;
              
            } else if (((arma_any_na[i] != 0) || std::isnan(x(i, j)) || std::isnan(x(i, k))) &&
              (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j)) && !std::isnan(x(i - width, k))) {
              
              n_obs -= 1;
              
            }
            
            if ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)) || std::isnan(x(i - width, k))) {
              
              w_old = 0;
              x_old = 0;
              y_old = 0;
              
            } else {
              
              w_old = arma_weights[n - width];
              x_old = x(i - width, j);
              y_old = x(i - width, k);
              
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
              if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(x(i, k)) &&
                  (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j)) && !std::isnan(x(i - width, k))) {
                
                sumsq_x = lambda * sumsq_x +
                  w_new * (x_new - mean_x) * (x_new - mean_prev_x) -
                  lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
                sumsq_y = lambda * sumsq_y +
                  w_new * (y_new - mean_y) * (y_new - mean_prev_y) -
                  lambda * w_old * (y_old - mean_y) * (y_old - mean_prev_y);
                
              } else if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(x(i, k)) &&
                ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)) || std::isnan(x(i - width, k)))) {
                
                sumsq_x = lambda * sumsq_x +
                  w_new * (x_new - mean_x) * (x_new - mean_prev_x);
                sumsq_y = lambda * sumsq_y +
                  w_new * (y_new - mean_y) * (y_new - mean_prev_y);
                
              } else if (((arma_any_na[i] != 0) || std::isnan(x(i, j)) || std::isnan(x(i, k))) &&
                (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j)) && !std::isnan(x(i - width, k))) {
                
                sumsq_x = lambda * sumsq_x -
                  lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
                sumsq_y = lambda * sumsq_y -
                  lambda * w_old * (y_old - mean_y) * (y_old - mean_prev_y);
                
              } else if ((arma_any_na[i] != 0) || std::isnan(x(i, j)) || std::isnan(x(i, k)) ||
                (arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)) || std::isnan(x(i - width, k))) {
                
                sumsq_x = lambda * sumsq_x;
                sumsq_y = lambda * sumsq_y;
                
              }
              
            }
            
            // compute the sum of squares
            if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(x(i, k)) &&
                (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j)) && !std::isnan(x(i - width, k))) {
              
              sumsq_xy = lambda * sumsq_xy +
                w_new * (x_new - mean_x) * (y_new - mean_prev_y) -
                lambda * w_old * (x_old - mean_x) * (y_old - mean_prev_y);
              
            } else if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(x(i, k)) &&
              ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)) || std::isnan(x(i - width, k)))) {
              
              sumsq_xy = lambda * sumsq_xy +
                w_new * (x_new - mean_x) * (y_new - mean_prev_y);
              
            } else if (((arma_any_na[i] != 0) || std::isnan(x(i, j)) || std::isnan(x(i, k))) &&
              (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j)) && !std::isnan(x(i - width, k))) {
              
              sumsq_xy = lambda * sumsq_xy -
                lambda * w_old * (x_old - mean_x) * (y_old - mean_prev_y);
              
            } else if ((arma_any_na[i] != 0) || std::isnan(x(i, j)) || std::isnan(x(i, k)) ||
              (arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)) || std::isnan(x(i - width, k))) {
              
              sumsq_xy = lambda * sumsq_xy;
              
            }
            
          }
          
          // don't compute if missing value and 'na_restore' argument is TRUE
          if ((!na_restore) || (na_restore && !std::isnan(x(i, j)) &&
              !std::isnan(x(i, k)))) {
              
              // compute the unbiased estimate of variance
              if ((n_obs > 1) && (n_obs >= min_obs)) {
                
                if (scale) {
                  
                  // don't compute if the standard deviation is zero
                  if ((sumsq_x < 0) || (sumsq_y < 0) ||
                      (sqrt(sumsq_x) <= sqrt(arma::datum::eps)) || (sqrt(sumsq_y) <= sqrt(arma::datum::eps))) {
                    arma_cov(j, k, i) = NA_REAL;
                  } else {
                    
                    if (std::abs(sumsq_xy) <= sqrt(arma::datum::eps)) {
                      arma_cov(j, k, i) = 0;
                    } else {
                      arma_cov(j, k, i) = sumsq_xy / (sqrt(sumsq_x) * sqrt(sumsq_y));
                    }
                    
                  }
                  
                } else if (!scale) {
                  
                  if (std::abs(sumsq_xy) <= sqrt(arma::datum::eps)) {
                    arma_cov(j, k, i) = 0;
                  } else {
                    arma_cov(j, k, i) = sumsq_xy / (sum_w - sumsq_w / sum_w);
                  }
                  
                }
                
              } else {
                arma_cov(j, k, i) = NA_REAL;
              }
              
          } else {
            
            // can be either NA or NaN
            if (std::isnan(x(i, j))) {
              arma_cov(j, k, i) = x(i, j);
            } else {
              arma_cov(j, k, i) = x(i, k);
            }
            
          }
          
          // covariance matrix is symmetric
          arma_cov(k, j, i) = arma_cov(j, k, i);
          
        }
        
      }
    }
  }
  
};

// 'Worker' function for computing rolling covariances using an online algorithm
struct RollCovOnlineMatXY : public Worker {
  
  const RMatrix<double> x;      // source
  const RMatrix<double> y;      // source
  const int n;
  const int n_rows_xy;
  const int n_cols_x;
  const int n_cols_y;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const bool scale;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::cube& arma_cov;         // destination (pass by reference)
  
  // initialize with source and destination
  RollCovOnlineMatXY(const NumericMatrix x, const NumericMatrix y,
                     const int n, const int n_rows_xy,
                     const int n_cols_x, const int n_cols_y,
                     const int width, const arma::vec arma_weights,
                     const bool center, const bool scale,
                     const int min_obs, const arma::uvec arma_any_na,
                     const bool na_restore, arma::cube& arma_cov)
    : x(x), y(y),
      n(n), n_rows_xy(n_rows_xy),
      n_cols_x(n_cols_x), n_cols_y(n_cols_y),
      width(width), arma_weights(arma_weights),
      center(center), scale(scale),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), arma_cov(arma_cov) { }
  
  // function call operator that iterates by column
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      for (int k = 0; k <= n_cols_y - 1; k++) {
        
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
          
          if ((arma_any_na[i] != 0) || std::isnan(x(i, j)) || std::isnan(y(i, k))) {
            
            w_new = 0;
            x_new = 0;
            y_new = 0;
            
          } else {
            
            w_new = arma_weights[n - 1];
            x_new = x(i, j);
            y_new = y(i, k);
            
          }
          
          // expanding window
          if (i < width) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
            if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(y(i, k))) {
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
              if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(y(i, k)) &&
                  (n_obs > 1)) {
                
                sumsq_x = lambda * sumsq_x +
                  w_new * (x_new - mean_x) * (x_new - mean_prev_x);
                sumsq_y = lambda * sumsq_y +
                  w_new * (y_new - mean_y) * (y_new - mean_prev_y);
                
              } else if ((arma_any_na[i] != 0) || std::isnan(x(i, j)) || std::isnan(y(i, k))) {
                
                sumsq_x = lambda * sumsq_x;
                sumsq_y = lambda * sumsq_y;
                
              } else if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(y(i, k)) &&
                (n_obs == 1) && !center) {
                
                sumsq_x = w_new * pow(x_new, (long double)2.0);
                sumsq_y = w_new * pow(y_new, (long double)2.0);
                
              }
              
            }
            
            // compute the sum of squares
            if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(y(i, k)) &&
                (n_obs > 1)) {
              
              sumsq_xy = lambda * sumsq_xy +
                w_new * (x_new - mean_x) * (y_new - mean_prev_y);
              
            } else if ((arma_any_na[i] != 0) || std::isnan(x(i, j)) || std::isnan(y(i, k))) {
              
              sumsq_xy = lambda * sumsq_xy;
              
            } else if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(y(i, k)) &&
              (n_obs == 1) && !center) {
              
              sumsq_xy = w_new * x_new * y_new;
              
            }
            
          }
          
          // rolling window
          if (i >= width) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
            if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(y(i, k)) &&
                ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)) || std::isnan(y(i - width, k)))) {
              
              n_obs += 1;
              
            } else if (((arma_any_na[i] != 0) || std::isnan(x(i, j)) || std::isnan(y(i, k))) &&
              (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j)) && !std::isnan(y(i - width, k))) {
              
              n_obs -= 1;
              
            }
            
            if ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)) || std::isnan(y(i - width, k))) {
              
              w_old = 0;
              x_old = 0;
              y_old = 0;
              
            } else {
              
              w_old = arma_weights[n - width];
              x_old = x(i - width, j);
              y_old = y(i - width, k);
              
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
              if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(y(i, k)) &&
                  (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j)) && !std::isnan(y(i - width, k))) {
                
                sumsq_x = lambda * sumsq_x +
                  w_new * (x_new - mean_x) * (x_new - mean_prev_x) -
                  lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
                sumsq_y = lambda * sumsq_y +
                  w_new * (y_new - mean_y) * (y_new - mean_prev_y) -
                  lambda * w_old * (y_old - mean_y) * (y_old - mean_prev_y);
                
              } else if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(y(i, k)) &&
                ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)) || std::isnan(y(i - width, k)))) {
                
                sumsq_x = lambda * sumsq_x +
                  w_new * (x_new - mean_x) * (x_new - mean_prev_x);
                sumsq_y = lambda * sumsq_y +
                  w_new * (y_new - mean_y) * (y_new - mean_prev_y);
                
              } else if (((arma_any_na[i] != 0) || std::isnan(x(i, j)) || std::isnan(y(i, k))) &&
                (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j)) && !std::isnan(y(i - width, k))) {
                
                sumsq_x = lambda * sumsq_x -
                  lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
                sumsq_y = lambda * sumsq_y -
                  lambda * w_old * (y_old - mean_y) * (y_old - mean_prev_y);
                
              } else if ((arma_any_na[i] != 0) || std::isnan(x(i, j)) || std::isnan(y(i, k)) ||
                (arma_any_na[i - width] == 0) || std::isnan(x(i - width, j)) || std::isnan(y(i - width, k))) {
                
                sumsq_x = lambda * sumsq_x;
                sumsq_y = lambda * sumsq_y;
                
              }
              
            }
            
            // compute the sum of squares
            if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(y(i, k)) &&
                (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j)) && !std::isnan(y(i - width, k))) {
              
              sumsq_xy = lambda * sumsq_xy +
                w_new * (x_new - mean_x) * (y_new - mean_prev_y) -
                lambda * w_old * (x_old - mean_x) * (y_old - mean_prev_y);
              
            } else if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(y(i, k)) &&
              ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)) || std::isnan(y(i - width, k)))) {
              
              sumsq_xy = lambda * sumsq_xy +
                w_new * (x_new - mean_x) * (y_new - mean_prev_y);
              
            } else if (((arma_any_na[i] != 0) || std::isnan(x(i, j)) || std::isnan(y(i, k))) &&
              (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j)) && !std::isnan(y(i - width, k))) {
              
              sumsq_xy = lambda * sumsq_xy -
                lambda * w_old * (x_old - mean_x) * (y_old - mean_prev_y);
              
            } else if ((arma_any_na[i] != 0) || std::isnan(x(i, j)) || std::isnan(y(i, k)) ||
              (arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)) || std::isnan(y(i - width, k))) {
              
              sumsq_xy = lambda * sumsq_xy;
              
            }
            
          }
          
          // don't compute if missing value and 'na_restore' argument is TRUE
          if ((!na_restore) || (na_restore && !std::isnan(x(i, j)) &&
              !std::isnan(y(i, k)))) {
              
              // compute the unbiased estimate of variance
              if ((n_obs > 1) && (n_obs >= min_obs)) {
                
                if (scale) {
                  
                  // don't compute if the standard deviation is zero
                  if ((sumsq_x < 0) || (sumsq_y < 0) ||
                      (sqrt(sumsq_x) <= sqrt(arma::datum::eps)) || (sqrt(sumsq_y) <= sqrt(arma::datum::eps))) {
                    arma_cov(j, k, i) = NA_REAL;
                  } else {
                    
                    if (std::abs(sumsq_xy) <= sqrt(arma::datum::eps)) {
                      arma_cov(j, k, i) = 0;
                    } else {
                      arma_cov(j, k, i) = sumsq_xy / (sqrt(sumsq_x) * sqrt(sumsq_y));
                    }
                    
                  }
                  
                } else if (!scale) {
                  
                  if (std::abs(sumsq_xy) <= sqrt(arma::datum::eps)) {
                    arma_cov(j, k, i) = 0;
                  } else {
                    arma_cov(j, k, i) = sumsq_xy / (sum_w - sumsq_w / sum_w);
                  }
                  
                }
                
              } else {
                arma_cov(j, k, i) = NA_REAL;
              }
              
          } else {
            
            // can be either NA or NaN
            if (std::isnan(x(i, j))) {
              arma_cov(j, k, i) = x(i, j);
            } else {
              arma_cov(j, k, i) = y(i, k);
            }
            
          }
          
        }
        
      }
    }
  }
  
};

// 'Worker' function for computing rolling covariances using a standard algorithm
struct RollCovBatchMatXX : public Worker {
  
  const RMatrix<double> x;       // source
  const int n;
  const int n_rows_xy;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const bool scale;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::cube& arma_cov;          // destination (pass by reference)
  
  // initialize with source and destination
  RollCovBatchMatXX(const NumericMatrix x, const int n,
                    const int n_rows_xy, const int n_cols_x,
                    const int width, const arma::vec arma_weights,
                    const bool center, const bool scale, 
                    const int min_obs, const arma::uvec arma_any_na,
                    const bool na_restore, arma::cube& arma_cov)
    : x(x), n(n),
      n_rows_xy(n_rows_xy), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      center(center), scale(scale),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), arma_cov(arma_cov) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 3D array (lower triangle)
      int n_unique = n_cols_x * (n_cols_x + 1) / 2;
      int i = z / n_unique;
      int z_unique = z % n_unique;
      int k = n_cols_x -
        floor((sqrt((long double)(4 * n_cols_x * (n_cols_x + 1) - (7 + 8 * z_unique))) - 1) / 2) - 1;
      int j = z_unique - n_cols_x * k + k * (k + 1) / 2;
      
      long double mean_x = 0;
      long double mean_y = 0;
      long double var_x = 0;
      long double var_y = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x(i, j)) &&
          !std::isnan(x(i, k)))) {
          
          if (center) {
            
            int count = 0;
            long double sum_w = 0;
            long double sum_x = 0;
            long double sum_y = 0;
            
            // number of observations is either the window size or,
            // for partial results, the number of the current row
            while ((width > count) && (i >= count)) {
              
              // don't include if missing value and 'any_na' argument is 1
              // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
              if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j)) &&
                  !std::isnan(x(i - count, k))) {
                  
                  // compute the rolling sum
                  sum_w += arma_weights[n - count - 1];
                sum_x += arma_weights[n - count - 1] * x(i - count, j);
                sum_y += arma_weights[n - count - 1] * x(i - count, k);
                
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
              
              // don't include if missing value and 'any_na' argument is 1
              // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
              if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j)) &&
                  !std::isnan(x(i - count, k))) {
                  
                  // compute the rolling sum of squares with 'center' argument
                  if (center) {
                    
                    sumsq_x += arma_weights[n - count - 1] *
                      pow(x(i - count, j) - mean_x, (long double)2.0);
                    sumsq_y += arma_weights[n - count - 1] *
                      pow(x(i - count, k) - mean_y, (long double)2.0);
                    
                  } else if (!center) {
                    
                    sumsq_x += arma_weights[n - count - 1] *
                      pow(x(i - count, j), 2.0);
                    sumsq_y += arma_weights[n - count - 1] *
                      pow(x(i - count, k), 2.0);
                    
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
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
            if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j)) &&
                !std::isnan(x(i - count, k))) {
                
                sum_w += arma_weights[n - count - 1];
              sumsq_w += pow(arma_weights[n - count - 1], 2.0);
              
              // compute the rolling sum of squares with 'center' argument
              if (center) {
                sumsq_xy += arma_weights[n - count - 1] * 
                  (x(i - count, j) - mean_x) * (x(i - count, k) - mean_y);
              } else if (!center) {
                sumsq_xy += arma_weights[n - count - 1] * 
                  x(i - count, j) * x(i - count, k);
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
                arma_cov(j, k, i) = NA_REAL;
              } else {
                
                if (std::abs(sumsq_xy) <= sqrt(arma::datum::eps)) {
                  arma_cov(j, k, i) = 0;
                } else {
                  arma_cov(j, k, i) = sumsq_xy / (sqrt(var_x) * sqrt(var_y));
                }
                
              }
              
            } else if (!scale) {
              
              if (std::abs(sumsq_xy) <= sqrt(arma::datum::eps)) {
                arma_cov(j, k, i) = 0;
              } else {
                arma_cov(j, k, i) = sumsq_xy / (sum_w - sumsq_w / sum_w);
              }
              
            }
            
          } else {
            arma_cov(j, k, i) = NA_REAL;
          }
          
      } else {
        
        // can be either NA or NaN
        if (std::isnan(x(i, j))) {
          arma_cov(j, k, i) = x(i, j);
        } else {
          arma_cov(j, k, i) = x(i, k);
        }
        
      }
      
      // covariance matrix is symmetric
      arma_cov(k, j, i) = arma_cov(j, k, i);
      
    }
  }
  
};

// 'Worker' function for computing rolling covariances using a standard algorithm
struct RollCovBatchMatXY : public Worker {
  
  const RMatrix<double> x;       // source
  const RMatrix<double> y;       // source
  const int n;
  const int n_rows_xy;
  const int n_cols_x;
  const int n_cols_y;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const bool scale;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::cube& arma_cov;          // destination (pass by reference)
  
  // initialize with source and destination
  RollCovBatchMatXY(const NumericMatrix x, const NumericMatrix y,
                    const int n, const int n_rows_xy,
                    const int n_cols_x, const int n_cols_y,
                    const int width, const arma::vec arma_weights,
                    const bool center, const bool scale,
                    const int min_obs, const arma::uvec arma_any_na,
                    const bool na_restore, arma::cube& arma_cov)
    : x(x), y(y),
      n(n), n_rows_xy(n_rows_xy),
      n_cols_x(n_cols_x), n_cols_y(n_cols_y),
      width(width), arma_weights(arma_weights),
      center(center), scale(scale),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), arma_cov(arma_cov) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 3D array
      int i = z % n_rows_xy;
      int j = z / (n_cols_y * n_rows_xy);
      int k = (z / n_rows_xy) % n_cols_y;
      
      long double mean_x = 0;
      long double mean_y = 0;
      long double var_x = 0;
      long double var_y = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x(i, j)) &&
          !std::isnan(y(i, k)))) {
          
          if (center) {
            
            int count = 0;
            long double sum_w = 0;
            long double sum_x = 0;
            long double sum_y = 0;
            
            // number of observations is either the window size or,
            // for partial results, the number of the current row
            while ((width > count) && (i >= count)) {
              
              // don't include if missing value and 'any_na' argument is 1
              // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
              if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j)) &&
                  !std::isnan(y(i - count, k))) {
                  
                  // compute the rolling sum
                  sum_w += arma_weights[n - count - 1];
                sum_x += arma_weights[n - count - 1] * x(i - count, j);
                sum_y += arma_weights[n - count - 1] * y(i - count, k);
                
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
              
              // don't include if missing value and 'any_na' argument is 1
              // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
              if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j)) &&
                  !std::isnan(y(i - count, k))) {
                  
                  // compute the rolling sum of squares with 'center' argument
                  if (center) {
                    
                    sumsq_x += arma_weights[n - count - 1] *
                      pow(x(i - count, j) - mean_x, (long double)2.0);
                    sumsq_y += arma_weights[n - count - 1] *
                      pow(y(i - count, k) - mean_y, (long double)2.0);
                    
                  } else if (!center) {
                    
                    sumsq_x += arma_weights[n - count - 1] *
                      pow(x(i - count, j), 2.0);
                    sumsq_y += arma_weights[n - count - 1] *
                      pow(y(i - count, k), 2.0);
                    
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
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
            if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j)) &&
                !std::isnan(y(i - count, k))) {
                
                sum_w += arma_weights[n - count - 1];
              sumsq_w += pow(arma_weights[n - count - 1], 2.0);
              
              // compute the rolling sum of squares with 'center' argument
              if (center) {
                sumsq_xy += arma_weights[n - count - 1] * 
                  (x(i - count, j) - mean_x) * (y(i - count, k) - mean_y);
              } else if (!center) {
                sumsq_xy += arma_weights[n - count - 1] * 
                  x(i - count, j) * y(i - count, k);
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
                arma_cov(j, k, i) = NA_REAL;
              } else {
                
                if (std::abs(sumsq_xy) <= sqrt(arma::datum::eps)) {
                  arma_cov(j, k, i) = 0;
                } else {
                  arma_cov(j, k, i) = sumsq_xy / (sqrt(var_x) * sqrt(var_y));
                }
                
              }
              
            } else if (!scale) {
              
              if (std::abs(sumsq_xy) <= sqrt(arma::datum::eps)) {
                arma_cov(j, k, i) = 0;
              } else {
                arma_cov(j, k, i) = sumsq_xy / (sum_w - sumsq_w / sum_w);
              }
              
            }
            
          } else {
            arma_cov(j, k, i) = NA_REAL;
          }
          
      } else {
        
        // can be either NA or NaN
        if (std::isnan(x(i, j))) {
          arma_cov(j, k, i) = x(i, j);
        } else {
          arma_cov(j, k, i) = y(i, k);
        }
        
      }
      
    }
  }
  
};

// 'Worker' function for computing rolling covariances using an online algorithm
struct RollCovOnlineMatLm : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_xy;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const bool intercept;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::vec& arma_n_obs;        // destination (pass by reference)
  arma::vec& arma_sum_w;
  arma::mat& arma_mean;
  arma::cube& arma_cov;
  
  // initialize with source and destination
  RollCovOnlineMatLm(const NumericMatrix x, const int n,
                     const int n_rows_xy, const int n_cols_x,
                     const int width, const arma::vec arma_weights,
                     const bool intercept, const int min_obs,
                     const arma::uvec arma_any_na, const bool na_restore,
                     arma::vec& arma_n_obs, arma::vec& arma_sum_w,
                     arma::mat& arma_mean, arma::cube& arma_cov)
    : x(x), n(n),
      n_rows_xy(n_rows_xy), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      intercept(intercept), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      arma_n_obs(arma_n_obs), arma_sum_w(arma_sum_w),
      arma_mean(arma_mean), arma_cov(arma_cov) { }
  
  // function call operator that iterates by column
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      for (std::size_t k = 0; k <= j; k++) {
        
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
        long double sumsq_xy = 0;
        // long double mean_prev_x = 0;
        long double mean_prev_y = 0;
        long double mean_x = 0;
        long double mean_y = 0;
        
        if (width > 1) {
          lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed!
        } else {
          lambda = arma_weights[n - 1];
        }
        
        for (int i = 0; i < n_rows_xy; i++) {
          
          if ((arma_any_na[i] != 0) || std::isnan(x(i, j)) || std::isnan(x(i, k))) {
            
            w_new = 0;
            x_new = 0;
            y_new = 0;
            
          } else {
            
            w_new = arma_weights[n - 1];
            x_new = x(i, j);
            y_new = x(i, k);
            
          }
          
          // expanding window
          if (i < width) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
            if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(x(i, k))) {
              n_obs += 1;
            }
            
            if (width > 1) {
              
              sum_w = lambda * sum_w + w_new;
              sum_x = lambda * sum_x + w_new * x_new;
              sum_y = lambda * sum_y + w_new * y_new;
              sumsq_w = pow(lambda, (long double)2.0) * sumsq_w + pow(w_new, (long double)2.0);
              
            } else {
              
              sum_w = w_new;
              sum_x = w_new * x_new;
              sum_y = w_new * y_new;
              sumsq_w = pow(w_new, (long double)2.0);
              
            }
            
            if (intercept && (n_obs > 0)) {
              
              // compute the mean
              // mean_prev_x = mean_x;
              mean_prev_y = mean_y;
              mean_x = sum_x / sum_w;
              mean_y = sum_y / sum_w;
              
            }
            
            // compute the sum of squares
            if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(x(i, k)) && (n_obs > 1)) {
              
              if (width > 1) {
                sumsq_xy = lambda * sumsq_xy +
                  w_new * (x_new - mean_x) * (y_new - mean_prev_y);
              } else {
                sumsq_xy = w_new * (x_new - mean_x) * (y_new - mean_prev_y);
              }
              
            } else if ((arma_any_na[i] != 0) || std::isnan(x(i, j)) || std::isnan(x(i, k))) {
              
              sumsq_xy = lambda * sumsq_xy;
              
            } else if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(x(i, k)) &&
              (n_obs == 1) && !intercept) {
              
              sumsq_xy = w_new * x_new * y_new;
              
            }
            
          }
          
          // rolling window
          if (i >= width) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
            if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(x(i, k)) &&
                ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)) || std::isnan(x(i - width, k)))) {
              
              n_obs += 1;
              
            } else if (((arma_any_na[i] != 0) || std::isnan(x(i, j)) || std::isnan(x(i, k))) &&
              (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j)) && !std::isnan(x(i - width, k))) {
              
              n_obs -= 1;
              
            }
            
            if ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)) || std::isnan(x(i - width, k))) {
              
              w_old = 0;
              x_old = 0;
              y_old = 0;
              
            } else {
              
              w_old = arma_weights[n - width];
              x_old = x(i - width, j);
              y_old = x(i - width, k);
              
            }
            
            if (width > 1) {
              
              sum_w = lambda * sum_w + w_new - lambda * w_old;
              sum_x = lambda * sum_x + w_new * x_new - lambda * w_old * x_old;
              sum_y = lambda * sum_y + w_new * y_new - lambda * w_old * y_old;
              sumsq_w = pow(lambda, (long double)2.0) * sumsq_w +
                pow(w_new, (long double)2.0) - pow(lambda * w_old, (long double)2.0);
              
            } else {
              
              sum_w = w_new;
              sum_x = w_new * x_new;
              sum_y = w_new * y_new;
              sumsq_w = pow(w_new, (long double)2.0);
              
            }
            
            if (intercept && (n_obs > 0)) {
              
              // compute the mean
              // mean_prev_x = mean_x;
              mean_prev_y = mean_y;
              mean_x = sum_x / sum_w;
              mean_y = sum_y / sum_w;
              
            }
            
            // compute the sum of squares
            if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(x(i, k)) &&
                (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j)) && !std::isnan(x(i - width, k))) {
              
              if (width > 1) {
                
                sumsq_xy = lambda * sumsq_xy +
                  w_new * (x_new - mean_x) * (y_new - mean_prev_y) -
                  lambda * w_old * (x_old - mean_x) * (y_old - mean_prev_y);
                
              } else {
                
                sumsq_xy = w_new * (x_new - mean_x) * (y_new - mean_prev_y);
                
              }
              
            } else if ((arma_any_na[i] == 0) && !std::isnan(x(i, j)) && !std::isnan(x(i, k)) &&
              ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)) || std::isnan(x(i - width, k)))) {
              
              if (width > 1) {
                
                sumsq_xy = lambda * sumsq_xy +
                  w_new * (x_new - mean_x) * (y_new - mean_prev_y);
                
              } else {
                sumsq_xy = w_new * (x_new - mean_x) * (y_new - mean_prev_y);
              }
              
            } else if (((arma_any_na[i] != 0) || std::isnan(x(i, j)) || std::isnan(x(i, k))) &&
              (arma_any_na[i - width] == 0) && !std::isnan(x(i - width, j)) && !std::isnan(x(i - width, k))) {
              
              sumsq_xy = lambda * sumsq_xy -
                lambda * w_old * (x_old - mean_x) * (y_old - mean_prev_y);
              
            } else if ((arma_any_na[i] != 0) || std::isnan(x(i, j)) || std::isnan(x(i, k)) ||
              (arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)) || std::isnan(x(i - width, k))) {
              
              sumsq_xy = lambda * sumsq_xy;
              
            }
            
          }
          
          // degrees of freedom and intercept std.error
          if (((int)j == n_cols_x - 1) && ((int)k == n_cols_x - 1)) {
            
            arma_n_obs[i] = n_obs;
            arma_sum_w[i] = sum_w;
            
          }
          
          // intercept
          if (j == k) {
            arma_mean(i, j) = mean_x;
          }
          
          // don't compute if missing value and 'na_restore' argument is TRUE
          if ((!na_restore) || (na_restore && !std::isnan(x(i, j)) &&
              !std::isnan(x(i, k)))) {
              
              // compute the unbiased estimate of variance
              // if ((n_obs > 1) && (n_obs >= min_obs)) {
              if (n_obs >= min_obs) {
                
                if (std::abs(sumsq_xy) <= sqrt(arma::datum::eps)) {
                  arma_cov(j, k, i) = 0;
                } else {
                  arma_cov(j, k, i) = sumsq_xy;
                }
                
              } else {
                arma_cov(j, k, i) = NA_REAL;
              }
              
          } else {
            
            // can be either NA or NaN
            if (std::isnan(x(i, j))) {
              arma_cov(j, k, i) = x(i, j);
            } else {
              arma_cov(j, k, i) = x(i, k);
            }
            
          }
          
          // covariance matrix is symmetric
          arma_cov(k, j, i) = arma_cov(j, k, i);
          
        }
        
      }
    }
  }
  
};

// 'Worker' function for computing rolling covariances using a standard algorithm
struct RollCovBatchMatLm : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_xy;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const bool intercept;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::vec& arma_n_obs;        // destination (pass by reference)
  arma::vec& arma_sum_w;
  arma::mat& arma_mean;
  arma::cube& arma_cov;
  
  // initialize with source and destination
  RollCovBatchMatLm(const NumericMatrix x, const int n,
                    const int n_rows_xy, const int n_cols_x,
                    const int width, const arma::vec arma_weights,
                    const bool intercept, const int min_obs,
                    const arma::uvec arma_any_na, const bool na_restore,
                    arma::vec& arma_n_obs, arma::vec& arma_sum_w,
                    arma::mat& arma_mean, arma::cube& arma_cov)
    : x(x), n(n),
      n_rows_xy(n_rows_xy), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      intercept(intercept), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      arma_n_obs(arma_n_obs), arma_sum_w(arma_sum_w),
      arma_mean(arma_mean), arma_cov(arma_cov) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 3D array (lower triangle)
      int n_unique = n_cols_x * (n_cols_x + 1) / 2;
      int i = z / n_unique;
      int z_unique = z % n_unique;
      int k = n_cols_x -
        floor((sqrt((long double)(4 * n_cols_x * (n_cols_x + 1) - (7 + 8 * z_unique))) - 1) / 2) - 1;
      int j = z_unique - n_cols_x * k + k * (k + 1) / 2;
      
      long double mean_x = 0;
      long double mean_y = 0;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if ((!na_restore) || (na_restore && !std::isnan(x(i, j)) &&
          !std::isnan(x(i, k)))) {
          
          if (intercept) {
            
            int count = 0;
            long double sum_w = 0;
            long double sum_x = 0;
            long double sum_y = 0;
            
            // number of observations is either the window size or,
            // for partial results, the number of the current row
            while ((width > count) && (i >= count)) {
              
              // don't include if missing value and 'any_na' argument is 1
              // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
              if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j)) &&
                  !std::isnan(x(i - count, k))) {
                  
                  // compute the rolling sum
                  sum_w += arma_weights[n - count - 1];
                sum_x += arma_weights[n - count - 1] * x(i - count, j);
                sum_y += arma_weights[n - count - 1] * x(i - count, k);
                
              }
              
              count += 1;
              
            }
            
            // compute the mean
            mean_x = sum_x / sum_w;
            mean_y = sum_y / sum_w;
            
          }
          
          int count = 0;
          int n_obs = 0;
          long double sum_w = 0;
          long double sumsq_w = 0;
          long double sumsq_xy = 0;
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
            if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j)) &&
                !std::isnan(x(i - count, k))) {
                
                sum_w += arma_weights[n - count - 1];
              sumsq_w += pow(arma_weights[n - count - 1], 2.0);
              
              // compute the rolling sum of squares with 'center' argument
              if (intercept) {
                sumsq_xy += arma_weights[n - count - 1] * 
                  (x(i - count, j) - mean_x) * (x(i - count, k) - mean_y);
              } else if (!intercept) {
                sumsq_xy += arma_weights[n - count - 1] * 
                  x(i - count, j) * x(i - count, k);
              }
              
              n_obs += 1;
              
            }
            
            count += 1;
            
          }
          
          // degrees of freedom and intercept std.error
          if ((j == n_cols_x - 1) && (k == n_cols_x - 1)) {
            
            arma_n_obs[i] = n_obs;
            arma_sum_w[i] = sum_w;
            
          }
          
          // intercept
          if (j == k) {
            arma_mean(i, j) = mean_x;
          }
          
          // compute the unbiased estimate of covariance
          // if ((n_obs > 1) && (n_obs >= min_obs)) {
          if (n_obs >= min_obs) {
            
            if (std::abs(sumsq_xy) <= sqrt(arma::datum::eps)) {
              arma_cov(j, k, i) = 0;
            } else {
              arma_cov(j, k, i) = sumsq_xy;
            }
            
          } else {
            arma_cov(j, k, i) = NA_REAL;
          }
          
      } else {
        
        // can be either NA or NaN
        if (std::isnan(x(i, j))) {
          arma_cov(j, k, i) = x(i, j);
        } else {
          arma_cov(j, k, i) = x(i, k);
        }
        
      }
      
      // covariance matrix is symmetric
      arma_cov(k, j, i) = arma_cov(j, k, i);
      
    }
  }
  
};

// 'Worker' function for rolling linear models
struct RollLmMatInterceptTRUE : public Worker {
  
  const arma::cube arma_cov;    // source
  const int n;
  const int n_rows_xy;
  const int n_cols_x;
  const int width;
  const arma::vec arma_n_obs;
  const arma::vec arma_sum_w;
  const arma::mat arma_mean;
  arma::mat& arma_coef;         // destination (pass by reference)
  arma::vec& arma_rsq;
  arma::mat& arma_se;
  
  // initialize with source and destination
  RollLmMatInterceptTRUE(const arma::cube arma_cov, const int n,
                         const int n_rows_xy, const int n_cols_x,
                         const int width, const arma::vec arma_n_obs,
                         const arma::vec arma_sum_w, const arma::mat arma_mean,
                         arma::mat& arma_coef, arma::vec& arma_rsq,
                         arma::mat& arma_se)
    : arma_cov(arma_cov), n(n),
      n_rows_xy(n_rows_xy), n_cols_x(n_cols_x),
      width(width), arma_n_obs(arma_n_obs),
      arma_sum_w(arma_sum_w), arma_mean(arma_mean),
      arma_coef(arma_coef), arma_rsq(arma_rsq),
      arma_se(arma_se) { }
  
  // function call operator that iterates by slice
  void operator()(std::size_t begin_slice, std::size_t end_slice) {
    for (std::size_t i = begin_slice; i < end_slice; i++) {
      
      arma::mat sigma = arma_cov.slice(i);
      arma::mat A = sigma.submat(0, 0, n_cols_x - 2, n_cols_x - 2);
      arma::mat b = sigma.submat(0, n_cols_x - 1, n_cols_x - 2, n_cols_x - 1);
      arma::vec coef(n_cols_x - 1);
      
      // check if missing value is present
      bool any_na = sigma.has_nan();
      
      // don't compute if missing value 
      if (!any_na) {
        
        // check if solution is found 
        bool status_solve = arma::solve(coef, A, b, arma::solve_opts::no_approx);
        int df_fit = n_cols_x;
        
        // don't find approximate solution for rank deficient system,
        // and the width and current row must be greater than the
        // number of variables
        if (status_solve && (arma_n_obs[i] >= df_fit)) {
          
          // intercept
          arma::mat mean_x = arma_mean.submat(i, 0, i, n_cols_x - 2);
          arma_coef(i, 0) = arma_mean(i, n_cols_x - 1) -
            as_scalar(mean_x * coef);
          
          // coefficients
          arma::mat trans_coef = trans(coef);
          arma_coef.submat(i, 1, i, n_cols_x - 1) = trans_coef;
          
          // r-squared
          long double var_y = sigma(n_cols_x - 1, n_cols_x - 1);
          if ((var_y < 0) || (sqrt(var_y) <= sqrt(arma::datum::eps))) {
            arma_rsq[i] = NA_REAL;
          } else {
            arma_rsq[i] = as_scalar(trans_coef * A * coef) / var_y;
          }
          
          // check if matrix is singular
          arma::mat A_inv(n_cols_x, n_cols_x);
          bool status_inv = arma::inv(A_inv, A);
          int df_resid = arma_n_obs[i] - n_cols_x;
          
          if (status_inv && (df_resid > 0)) {
            
            // residual variance
            long double var_resid = (1 - arma_rsq[i]) * var_y / df_resid;
            
            // standard errors
            if ((var_resid < 0) || (sqrt(var_resid) <= sqrt(arma::datum::eps))) {
              var_resid = 0;
            }
            
            arma_se(i, 0) = sqrt(var_resid * (1 / arma_sum_w[i] +
              as_scalar(mean_x * A_inv * trans(mean_x))));
            arma_se.submat(i, 1, i, n_cols_x - 1) = sqrt(var_resid * trans(diagvec(A_inv)));
            
          } else {
            
            arma::rowvec no_solution(n_cols_x);
            no_solution.fill(NA_REAL);
            
            arma_se.row(i) = no_solution;
            
          }
          
        } else {
          
          arma::rowvec no_solution(n_cols_x);
          no_solution.fill(NA_REAL);
          
          arma_coef.row(i) = no_solution;
          arma_rsq[i] = NA_REAL;
          arma_se.row(i) = no_solution;
          
        }
        
      } else {
        
        arma::rowvec no_solution(n_cols_x);
        no_solution.fill(NA_REAL);
        
        arma_coef.row(i) = no_solution;
        arma_rsq[i] = NA_REAL;
        arma_se.row(i) = no_solution;
        
      }
      
    }
  }
  
};

// 'Worker' function for rolling linear models
struct RollLmMatInterceptFALSE : public Worker {
  
  const arma::cube arma_cov;    // source
  const int n;
  const int n_rows_xy;
  const int n_cols_x;
  const int width;
  const arma::vec arma_n_obs;
  const arma::vec arma_sum_w;
  arma::mat& arma_coef;         // destination (pass by reference)
  arma::vec& arma_rsq;
  arma::mat& arma_se;
  
  // initialize with source and destination
  RollLmMatInterceptFALSE(const arma::cube arma_cov, const int n,
                          const int n_rows_xy, const int n_cols_x,
                          const int width, const arma::vec arma_n_obs,
                          const arma::vec arma_sum_w, arma::mat& arma_coef,
                          arma::vec& arma_rsq, arma::mat& arma_se)
    : arma_cov(arma_cov), n(n),
      n_rows_xy(n_rows_xy), n_cols_x(n_cols_x),
      width(width), arma_n_obs(arma_n_obs),
      arma_sum_w(arma_sum_w), arma_coef(arma_coef),
      arma_rsq(arma_rsq), arma_se(arma_se) { }
  
  // function call operator that iterates by slice
  void operator()(std::size_t begin_slice, std::size_t end_slice) {
    for (std::size_t i = begin_slice; i < end_slice; i++) {
      
      arma::mat sigma = arma_cov.slice(i);
      arma::mat A = sigma.submat(0, 0, n_cols_x - 2, n_cols_x - 2);
      arma::mat b = sigma.submat(0, n_cols_x - 1, n_cols_x - 2, n_cols_x - 1);
      arma::vec coef(n_cols_x - 2);
      
      // check if missing value is present
      bool any_na = sigma.has_nan();
      
      // don't compute if missing value 
      if (!any_na) {
        
        // check if solution is found      
        bool status_solve = arma::solve(coef, A, b, arma::solve_opts::no_approx);
        int df_fit = n_cols_x - 1;
        
        // don't find approximate solution for rank deficient system,
        // and the width and current row must be greater than the
        // number of variables
        if (status_solve && (arma_n_obs[i] >= df_fit)) {
          
          // coefficients
          arma::mat trans_coef = trans(coef);
          arma_coef.row(i) = trans_coef;
          
          // r-squared
          long double var_y = sigma(n_cols_x - 1, n_cols_x - 1);
          if ((var_y < 0) || (sqrt(var_y) <= sqrt(arma::datum::eps))) {
            arma_rsq[i] = NA_REAL;
          } else {
            arma_rsq[i] = as_scalar(trans_coef * A * coef) / var_y;
          }
          
          // check if matrix is singular
          arma::mat A_inv(n_cols_x, n_cols_x);
          bool status_inv = arma::inv(A_inv, A);
          int df_resid = arma_n_obs[i] - n_cols_x + 1;
          
          if (status_inv && (df_resid > 0)) {
            
            // residual variance
            long double var_resid = (1 - arma_rsq[i]) * var_y / df_resid;
            
            // standard errors
            if ((var_resid < 0) || (sqrt(var_resid) <= sqrt(arma::datum::eps))) {
              var_resid = 0;
            }
            
            arma_se.row(i) = sqrt(var_resid * trans(diagvec(A_inv)));
            
          } else {
            
            arma::vec no_solution(n_cols_x - 1);
            no_solution.fill(NA_REAL);
            
            arma_se.row(i) = trans(no_solution);
            
          }
          
        } else {
          
          arma::vec no_solution(n_cols_x - 1);
          no_solution.fill(NA_REAL);
          
          arma_coef.row(i) = trans(no_solution);
          arma_rsq[i] = NA_REAL;
          arma_se.row(i) = trans(no_solution);
          
        }
        
      } else {
        
        arma::vec no_solution(n_cols_x - 1);
        no_solution.fill(NA_REAL);
        
        arma_coef.row(i) = trans(no_solution);
        arma_rsq[i] = NA_REAL;
        arma_se.row(i) = trans(no_solution);
        
      }
      
    }
  }
  
};

#endif