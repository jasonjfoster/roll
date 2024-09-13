#ifndef ROLL_MAT_H
#define ROLL_MAT_H

#define ARMA_WARN_LEVEL 0

#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

void update_n_obs(int& n_obs, const  bool& is_na,
                  const bool& is_na_old, const int& i,
                  const int& width);

arma::ivec stl_sort_index(arma::vec& x);

namespace roll {

// 'Worker' function for computing the rolling statistic using an online algorithm
struct RollAnyOnlineMat : public Worker {
  
  const RMatrix<int> x;         // source
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  RMatrix<int> rcpp_any;        // destination (pass by reference)
  
  // initialize with source and destination
  RollAnyOnlineMat(const IntegerMatrix x, const int n_rows_x,
                   const int n_cols_x, const int width,
                   const int min_obs, const arma::uvec arma_any_na,
                   const bool na_restore, IntegerMatrix rcpp_any)
    : x(x), n_rows_x(n_rows_x),
      n_cols_x(n_cols_x), width(width),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), rcpp_any(rcpp_any) { }
  
  // function call operator that iterates by column
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      
      bool is_na = false;
      bool is_na_old = false;
      int n_obs = 0;
      int sum_x = 0;
      
      for (int i = 0; i < n_rows_x; i++) {
        
        is_na = (arma_any_na[i] != 0) || (x(i, j) == NA_INTEGER);
        
        if (i >= width) {
          is_na_old = (arma_any_na[i - width] != 0) || (x(i - width, j) == NA_INTEGER);
        }
        
        update_n_obs(n_obs, is_na, is_na_old, i, width);
        
        if (!is_na) {
          
          // compute the sum
          if (x(i, j) != 0) {
            sum_x += 1;
          }
          
        }
        
        // rolling window
        if (i >= width) {
          if (!is_na_old) {
            if (x(i - width, j) != 0) {
              sum_x -= 1;
            }
          }
        }
        
        // don't compute if missing value and 'na_restore' argument is TRUE
        if (!na_restore || (x(i, j) != NA_INTEGER)) {
          
          if (n_obs >= min_obs) {
            
            if (sum_x > 0) {
              rcpp_any(i, j) = 1;
            } else if ((i >= width && n_obs == width) || (i < width && n_obs == i + 1)) {
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

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollAnyOfflineMat : public Worker {
  
  const RMatrix<int> x;         // source
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  RMatrix<int> rcpp_any;        // destination (pass by reference)
  
  // initialize with source and destination
  RollAnyOfflineMat(const IntegerMatrix x, const int n_rows_x,
                    const int n_cols_x, const int width,
                    const int min_obs, const arma::uvec arma_any_na,
                    const bool na_restore, IntegerMatrix rcpp_any)
    : x(x), n_rows_x(n_rows_x),
      n_cols_x(n_cols_x), width(width),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), rcpp_any(rcpp_any) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z / n_cols_x;
      int j = z % n_cols_x;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if (!na_restore || (x(i, j) != NA_INTEGER)) {
        
        int n_obs = 0;
        int sum_x = 0;
        bool is_na = false;
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        for (int count = 0; width > count && i >= count; count++) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          is_na = (arma_any_na[i - count] != 0) || (x(i - count, j) == NA_INTEGER);
          
          if (!is_na) {
            
            // compute the sum
            if (x(i - count, j) == 1) {
              sum_x += 1;
            }
            
            n_obs += 1;
            
          }
          
        }
        
        if (n_obs >= min_obs) {
          
          if (sum_x > 0) {
            rcpp_any(i, j) = 1;
          } else if ((i >= width && n_obs == width) || (i < width && n_obs == i + 1)) {
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

// 'Worker' function for computing the rolling statistic using an online algorithm
struct RollAllOnlineMat : public Worker {
  
  const RMatrix<int> x;         // source
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  RMatrix<int> rcpp_all;        // destination (pass by reference)
  
  // initialize with source and destination
  RollAllOnlineMat(const IntegerMatrix x, const int n_rows_x,
                   const int n_cols_x, const int width,
                   const int min_obs, const arma::uvec arma_any_na,
                   const bool na_restore, IntegerMatrix rcpp_all)
    : x(x), n_rows_x(n_rows_x),
      n_cols_x(n_cols_x), width(width),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), rcpp_all(rcpp_all) { }
  
  // function call operator that iterates by column
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      
      bool is_na = false;
      bool is_na_old = false;
      int n_obs = 0;
      int sum_x = 0;
      
      for (int i = 0; i < n_rows_x; i++) {
        
        is_na = (arma_any_na[i] != 0) || (x(i, j) == NA_INTEGER);
        
        if (i >= width) {
          is_na_old = (arma_any_na[i - width] != 0) || (x(i - width, j) == NA_INTEGER);
        }
        
        update_n_obs(n_obs, is_na, is_na_old, i, width);
        
        if (!is_na) {
          
          // compute the sum
          if (x(i, j) == 0) {
            sum_x += 1;
          }
          
        }
        
        // rolling window
        if (i >= width) {
          if (!is_na_old) {
            if (x(i - width, j) == 0) {
              sum_x -= 1;
            }
          }
        }
        
        // don't compute if missing value and 'na_restore' argument is TRUE
        if (!na_restore || (x(i, j) != NA_INTEGER)) {
          
          if (n_obs >= min_obs) {
            
            if (sum_x > 0) {
              rcpp_all(i, j) = 0;
            } else if ((i >= width && n_obs == width) || (i < width && n_obs == i + 1)) {
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

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollAllOfflineMat : public Worker {
  
  const RMatrix<int> x;         // source
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  RMatrix<int> rcpp_all;        // destination (pass by reference)
  
  // initialize with source and destination
  RollAllOfflineMat(const IntegerMatrix x, const int n_rows_x,
                    const int n_cols_x, const int width,
                    const int min_obs, const arma::uvec arma_any_na,
                    const bool na_restore, IntegerMatrix rcpp_all)
    : x(x), n_rows_x(n_rows_x),
      n_cols_x(n_cols_x), width(width),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), rcpp_all(rcpp_all) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z / n_cols_x;
      int j = z % n_cols_x;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if (!na_restore || (x(i, j) != NA_INTEGER)) {
        
        int n_obs = 0;
        int sum_x = 0;
        int is_na = false;
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        for (int count = 0; width > count && i >= count; count++) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          is_na = (arma_any_na[i - count] != 0) || (x(i - count, j) == NA_INTEGER);
          
          if (!is_na) {
            
            // compute the sum
            if (x(i - count, j) == 0) {
              sum_x += 1;
            }
            
            n_obs += 1;
            
          }
          
        }
        
        if (n_obs >= min_obs) {
          
          if (sum_x > 0) {
            rcpp_all(i, j) = 0;
          } else if ((i >= width && n_obs == width) || (i < width && n_obs == i + 1)) {
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

// 'Worker' function for computing the rolling statistic using an online algorithm
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
      
      if (arma_weights[n - 1] == 0) {
        lambda = 1;
      } else if (width > 1) {
        lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed
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

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollSumOfflineMat : public Worker {
  
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
  RollSumOfflineMat(const NumericMatrix x, const int n,
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
            
            // compute the sum
            sum_x += arma_weights[n - count - 1] * x(i - count, j);
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
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

// 'Worker' function for computing the rolling statistic using an online algorithm
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
      int n_zero = 0;
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
      
      if (arma_weights[n - 1] == 0) {
        lambda = 1;
      } else if (width > 1) {
        lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed
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
            
            if (x(i, j) == 0) {
              
              x_new = 1;
              n_zero += 1;
              
            } else {
              x_new = x(i, j);
            }
            
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
            
            if (x(i, j) == 0) {
              
              x_new = 1;
              n_zero += 1;
              
            } else {
              x_new = x(i, j);
            }
            
          }
          
          if ((arma_any_na[i - width] != 0) || std::isnan(x(i - width, j))) {
            
            w_old = 1;
            x_old = 1;
            
          } else {
            
            if (x(i - width, j) == 0) {
              
              x_old = 1;
              n_zero -= 1;
              
            } else {
              x_old = x(i - width, j);
            }
            
            if (arma_weights[n - width] == 0) {
              w_old = 1;
            } else {
              w_old = arma_weights[n - width];
            }
            
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
          
          if (n_obs >= min_obs) {
            
            if (n_zero == 0) {
              arma_prod(i, j) = prod_w * prod_x;
            } else {
              arma_prod(i, j) = 0;
            }
            
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

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollProdOfflineMat : public Worker {
  
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
  RollProdOfflineMat(const NumericMatrix x, const int n,
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
            
            // compute the product
            prod_x *= arma_weights[n - count - 1] * x(i - count, j);
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
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

// 'Worker' function for computing the rolling statistic using an online algorithm
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
        lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed
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

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollMeanOfflineMat : public Worker {
  
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
  RollMeanOfflineMat(const NumericMatrix x, const int n,
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
            
            // compute the sum
            sum_w += arma_weights[n - count - 1];
            sum_x += arma_weights[n - count - 1] * x(i - count, j);
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
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

// 'Worker' function for computing the rolling statistic using an online algorithm
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
  RMatrix<double> rcpp_min;     // destination (pass by reference)
  
  // initialize with source and destination
  RollMinOnlineMat(const NumericMatrix x, const int n,
                   const int n_rows_x, const int n_cols_x,
                   const int width, const arma::vec arma_weights,
                   const int min_obs, const arma::uvec arma_any_na,
                   const bool na_restore, NumericMatrix rcpp_min)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), rcpp_min(rcpp_min) { }
  
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
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j))) {
            n_obs += 1;
          }
          
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j))) {
            
            while (!deck.empty() && ((arma_any_na[deck.back()] != 0) ||
                   std::isnan(x(deck.back(), j)) || (x(i, j) < x(deck.back(), j)))) {
              
              deck.pop_back();
              
            }
            
            deck.push_back(i);
            
          }
          
          if (width > 1) {
            idxmin_x = deck.front();
          } else {
            idxmin_x = i;
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
                   std::isnan(x(deck.back(), j)) || (x(i, j) < x(deck.back(), j)))) {
              
              deck.pop_back();
              
            }
            
            deck.push_back(i);
            
          }
          
          while (!deck.empty() && (deck.front() <= i - width)) {
            deck.pop_front();
          }
          
          if (width > 1) {
            idxmin_x = deck.front();
          } else {
            idxmin_x = i;
          }
          
        }
        
        // don't compute if missing value and 'na_restore' argument is TRUE
        if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
          
          if (n_obs >= min_obs) {
            rcpp_min(i, j) = x(idxmin_x, j);
          } else {
            rcpp_min(i, j) = NA_REAL;
          }
          
        } else {
          
          // can be either NA or NaN
          rcpp_min(i, j) = x(i, j);
          
        }
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollMinOfflineMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  RMatrix<double> rcpp_min;     // destination (pass by reference)
  
  // initialize with source and destination
  RollMinOfflineMat(const NumericMatrix x, const int n,
                    const int n_rows_x, const int n_cols_x,
                    const int width, const arma::vec arma_weights,
                    const int min_obs, const arma::uvec arma_any_na,
                    const bool na_restore, NumericMatrix rcpp_min)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), rcpp_min(rcpp_min) { }
  
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
          if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j))) {
            
            // last element of sorted array
            // note: 'weights' must be greater than 0
            if ((arma_any_na[idxmin_x] != 0) || std::isnan(x(idxmin_x, j)) ||
                (x(i - count, j) <= x(idxmin_x, j))) {
              
              idxmin_x = i - count;
              
            }
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        if ((n_obs >= min_obs)) {
          rcpp_min(i, j) = x(idxmin_x, j);
        } else {
          rcpp_min(i, j) = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        rcpp_min(i, j) = x(i, j);
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an online algorithm
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
  RMatrix<double> rcpp_max;     // destination (pass by reference)
  
  // initialize with source and destination
  RollMaxOnlineMat(const NumericMatrix x, const int n,
                   const int n_rows_x, const int n_cols_x,
                   const int width, const arma::vec arma_weights,
                   const int min_obs, const arma::uvec arma_any_na,
                   const bool na_restore, NumericMatrix rcpp_max)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), rcpp_max(rcpp_max) { }
  
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
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j))) {
            n_obs += 1;
          }
          
          if ((arma_any_na[i] == 0) && !std::isnan(x(i, j))) {
            
            while (!deck.empty() && ((arma_any_na[deck.back()] != 0) || 
                   std::isnan(x(deck.back(), j)) || (x(i, j) > x(deck.back(), j)))) {
              
              deck.pop_back();
              
            }
            
            deck.push_back(i);
            
          }
          
          if (width > 1) {
            idxmax_x = deck.front();
          } else {
            idxmax_x = i;
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
                   std::isnan(x(deck.back(), j)) || (x(i, j) > x(deck.back(), j)))) {
              
              deck.pop_back();
              
            }
            
            deck.push_back(i);
            
          }
          
          while (!deck.empty() && (deck.front() <= i - width)) {
            deck.pop_front();
          }
          
          if (width > 1) {
            idxmax_x = deck.front();
          } else {
            idxmax_x = i;
          }
          
        }
        
        // don't compute if missing value and 'na_restore' argument is TRUE
        if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
          
          if (n_obs >= min_obs) {
            rcpp_max(i, j) = x(idxmax_x, j);
          } else {
            rcpp_max(i, j) = NA_REAL;
          }
          
        } else {
          
          // can be either NA or NaN
          rcpp_max(i, j) = x(i, j);
          
        }
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollMaxOfflineMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  RMatrix<double> rcpp_max;     // destination (pass by reference)
  
  // initialize with source and destination
  RollMaxOfflineMat(const NumericMatrix x, const int n,
                    const int n_rows_x, const int n_cols_x,
                    const int width, const arma::vec arma_weights,
                    const int min_obs, const arma::uvec arma_any_na,
                    const bool na_restore, NumericMatrix rcpp_max)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), rcpp_max(rcpp_max) { }
  
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
          if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j))) {
            
            // first element of sorted array
            // note: 'weights' must be greater than 0
            if ((arma_any_na[idxmax_x] != 0) || std::isnan(x(idxmax_x, j)) ||
                (x(i - count, j) >= x(idxmax_x, j))) {
              
              idxmax_x = i - count;
              
            }
            
            n_obs += 1;
            
          }
          
          count += 1;
          
        }
        
        if ((n_obs >= min_obs)) {
          rcpp_max(i, j) = x(idxmax_x, j);
        } else {
          rcpp_max(i, j) = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        rcpp_max(i, j) = x(i, j);
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an online algorithm
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
          
          while (!deck.empty() && (deck.front() <= i - width)) {
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
          
          if (n_obs >= min_obs) {
            rcpp_idxmin(i, j) = idxmin_x;
          } else {
            rcpp_idxmin(i, j) = NA_INTEGER;
          }
          
        } else {
          
          // can be either NA or NaN
          rcpp_idxmin(i, j) = (int)x(i, j);
          
        }
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollIdxMinOfflineMat : public Worker {
  
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
  RollIdxMinOfflineMat(const NumericMatrix x, const int n,
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
        
        if ((n_obs >= min_obs)) {
          
          if (i < width) {
            rcpp_idxmin(i, j) = idxmin_x + 1;
          } else if (i >= width) {
            rcpp_idxmin(i, j) = width - (i - idxmin_x);
          }
          
        } else {
          rcpp_idxmin(i, j) = NA_INTEGER;
        }
        
      } else {
        
        // can be either NA or NaN
        rcpp_idxmin(i, j) = (int)x(i, j);
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an online algorithm
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
          
          while (!deck.empty() && (deck.front() <= i - width)) {
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
          
          if (n_obs >= min_obs) {
            rcpp_idxmax(i, j) = idxmax_x;
          } else {
            rcpp_idxmax(i, j) = NA_INTEGER;
          }
          
        } else {
          
          // can be either NA or NaN
          rcpp_idxmax(i, j) = (int)x(i, j);
          
        }
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollIdxMaxOfflineMat : public Worker {
  
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
  RollIdxMaxOfflineMat(const NumericMatrix x, const int n,
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
        
        if ((n_obs >= min_obs)) {
          
          if (i < width) {
            rcpp_idxmax(i, j) = idxmax_x + 1;
          } else if (i >= width) {
            rcpp_idxmax(i, j) = width - (i - idxmax_x);
          }
          
        } else {
          rcpp_idxmax(i, j) = NA_INTEGER;
        }
        
      } else {
        
        // can be either NA or NaN
        rcpp_idxmax(i, j) = (int)x(i, j);
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an online algorithm
struct RollQuantileOnlineMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const double p;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  RMatrix<double> rcpp_quantile;// destination (pass by reference)
  
  // initialize with source and destination
  RollQuantileOnlineMat(const NumericMatrix x, const int n,
                        const int n_rows_x, const int n_cols_x,
                        const int width, const arma::vec arma_weights,
                        const double p, const int min_obs,
                        const arma::uvec arma_any_na, const bool na_restore,
                        NumericMatrix rcpp_quantile)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      p(p), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      rcpp_quantile(rcpp_quantile) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      
      int n_obs = 0;
      bool is_na = false;
      bool is_na_old = false;
      long double lambda = 0;
      long double w_new = 0;
      long double w_old = 0;
      long double x_new = 0;
      long double x_old = 0;
      long double sum_w = 0;
      long double sum_lower_w = 0;
      long double sum_upper_w = 0;
      std::multiset<std::tuple<long double, int>> mset_lower;
      std::multiset<std::tuple<long double, int>> mset_upper;
      
      int offset = 0;
      int n_size_x = 0;
      bool status = false;
      
      for (int i = 0; i < n_rows_x; i++) {
        
        // // uncomment if exponential-weights
        // sum_w = 0;
        // 
        // // number of observations is either the window size or,
        // // for partial results, the number of the current row
        // for (int count = 0; width > count && i >= count; count++) {
        //   
        //   // don't include if missing value and 'any_na' argument is 1
        //   // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
        //   is_na = (arma_any_na[i - count] != 0) || std::isnan(x(i - count, j));
        //   
        //   if (!is_na) {
        //     
        //     // compute the sum
        //     sum_w += arma_weights[n - count - 1];
        //     
        //   }
        //   
        // }
        
        is_na = (arma_any_na[i] != 0) || std::isnan(x(i, j));
        
        if (i >= width) {
          is_na_old = (arma_any_na[i - width] != 0) || (std::isnan(x(i - width, j)));
        }
        
        update_n_obs(n_obs, is_na, is_na_old, i, width);
        
        if (width > 1) {
          
          lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed
          sum_lower_w = lambda * sum_lower_w;
          sum_upper_w = lambda * sum_upper_w;
          
          if (!is_na) {
            
            w_new = arma_weights[n - 1];
            x_new = x(i, j);
            
            std::tuple<long double, int> tpl_new = std::make_tuple(x_new, i);
            
            sum_w += lambda * w_new; // comment if exponential-weights
            
            if (mset_lower.empty() || (x_new <= std::get<0>(*mset_lower.rbegin()))) {
              
              mset_lower.insert(tpl_new);
              sum_lower_w += w_new;
              
            } else {
              
              mset_upper.insert(tpl_new);
              sum_upper_w += w_new;
              
            }
            
            status = false;
            
            while (!status && (sum_lower_w / sum_w > p)) {
              
              offset = std::get<1>(*mset_lower.rbegin());
              n_size_x = i - offset;
              
              if (n - n_size_x > 0) {
                
                mset_upper.insert(*mset_lower.rbegin());
                sum_upper_w += arma_weights[std::max(0, n - n_size_x - 1)];
                sum_lower_w -= arma_weights[std::max(0, n - n_size_x - 1)];
                mset_lower.erase(*mset_lower.rbegin());
                
              } else {
                status = true;
              }
              
            }
            
            status = false;
            
            while (!status && !mset_upper.empty() && (sum_lower_w / sum_w < p)) {
              
              offset = std::get<1>(*mset_upper.begin());
              n_size_x = i - offset;
              
              if (n - n_size_x > 0) {
                
                mset_lower.insert(*mset_upper.begin());
                sum_lower_w += arma_weights[std::max(0, n - n_size_x - 1)];
                sum_upper_w -= arma_weights[std::max(0, n - n_size_x - 1)];
                mset_upper.erase(*mset_upper.begin());
                
              } else {
                status = true;
              }
              
            }
            
          }
          
          // rolling window
          if (i >= width) {
            
            if (!is_na_old) {
              
              w_old = arma_weights[n - width];
              x_old = x(i - width, j);
              
              std::tuple<long double, int> tpl_old = std::make_tuple(x_old, i - width);
              
              sum_w -= lambda * w_old; // comment if exponential-weights
              
              if (!mset_lower.empty() && (x_old <= std::get<0>(*mset_lower.rbegin()))) {
                
                mset_lower.erase(mset_lower.find(tpl_old));
                sum_lower_w -= lambda * w_old;
                
              } else {
                
                mset_upper.erase(mset_upper.find(tpl_old));
                sum_upper_w -= lambda * w_old;
                
              }
              
              status = false;
              
              while (!status && (sum_lower_w / sum_w > p)) {
                
                offset = std::get<1>(*mset_lower.rbegin());
                n_size_x = i - offset;
                
                if (n - n_size_x > 0) {
                  
                  mset_upper.insert(*mset_lower.rbegin());
                  sum_upper_w += arma_weights[std::max(0, n - n_size_x - 1)];
                  sum_lower_w -= arma_weights[std::max(0, n - n_size_x - 1)];
                  mset_lower.erase(*mset_lower.rbegin());
                  
                } else {
                  status = true;
                }
                
              }
              
              status = false;
              
              while (!status && !mset_upper.empty() && (sum_lower_w / sum_w < p)) {
                
                offset = std::get<1>(*mset_upper.begin());
                n_size_x = i - offset;
                
                if (n - n_size_x > 0) {
                  
                  mset_lower.insert(*mset_upper.begin());
                  sum_lower_w += arma_weights[std::max(0, n - n_size_x - 1)];
                  sum_upper_w -= arma_weights[std::max(0, n - n_size_x - 1)];
                  mset_upper.erase(*mset_upper.begin());
                  
                } else {
                  status = true;
                }
                
              }
              
            }
            
          }
          
          // don't compute if missing value and 'na_restore' argument is TRUE
          if (!na_restore || !std::isnan(x(i, j))) {
            
            if (n_obs >= min_obs) {
              
              if (p == 0) {
                rcpp_quantile(i, j) = std::get<0>(*mset_upper.begin());
              } else if (p == 1) {
                rcpp_quantile(i, j) = std::get<0>(*mset_lower.rbegin());
              } else if (std::fabs(sum_lower_w / sum_w - p) <= sqrt(arma::datum::eps)) {
                rcpp_quantile(i, j) = (std::get<0>(*mset_lower.rbegin()) + std::get<0>(*mset_upper.begin())) / 2;
              } else {
                rcpp_quantile(i, j) = std::get<0>(*mset_lower.rbegin());
              }
              
            } else {
              rcpp_quantile(i, j) = NA_REAL;
            }
            
          } else {
            
            // can be either NA or NaN
            rcpp_quantile(i, j) = x(i, j);
            
          }
          
        } else if (n_obs >= min_obs) {
          rcpp_quantile(i, j) = x(i, j);
        } else {
          rcpp_quantile(i, j) = NA_REAL;
        }
        
      }
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollQuantileOfflineMat : public Worker {
  
  const RMatrix<double> x;      // source
  const int n;
  const int n_rows_x;
  const int n_cols_x;
  const int width;
  const arma::vec arma_weights;
  const double p;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  RMatrix<double> rcpp_quantile;// destination (pass by reference)
  
  // initialize with source and destination
  RollQuantileOfflineMat(const NumericMatrix x, const int n,
                         const int n_rows_x, const int n_cols_x,
                         const int width, const arma::vec arma_weights,
                         const double p, const int min_obs,
                         const arma::uvec arma_any_na, const bool na_restore,
                         NumericMatrix rcpp_quantile)
    : x(x), n(n),
      n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      p(p), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      rcpp_quantile(rcpp_quantile) { }
  
  // function call operator that iterates by index
  void operator()(std::size_t begin_index, std::size_t end_index) {
    for (std::size_t z = begin_index; z < end_index; z++) {
      
      // from 1D to 2D array
      int i = z / n_cols_x;
      int j = z % n_cols_x;
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if (!na_restore || !std::isnan(x(i, j))) {
        
        int k = 0;
        bool is_na = false;
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
        arma::ivec sort_ix = stl_sort_index(x_subset);
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        for (int count = 0; count < n_size_x; count++) {
          
          k = sort_ix[n_size_x - count - 1];
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          is_na = arma_any_na_subset[k] == 1 || std::isnan(x_subset[k]);
          
          if (!is_na) {
            
            // compute the sum
            sum_w += arma_weights_subset[k];
            
          }
          
        }
        
        int n_obs = 0;
        int k_lower = 0;
        int idxquantile1_x = 0;
        int idxquantile2_x = 0;
        bool status1 = false;
        bool status2 = false;
        long double sum_upper_w = 0;
        long double sum_upper_w_temp = 0;
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        for (int count = 0; count < n_size_x; count++) {
          
          k = sort_ix[n_size_x - count - 1];
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          is_na = arma_any_na_subset[k] == 1 || std::isnan(x_subset[k]);
          
          if (!is_na) {
            
            // compute the sum
            sum_upper_w += arma_weights_subset[k];
            
            // last element of sorted array that is 'p' of 'weights'
            // note: 'weights' must be greater than 0
            if (!status1 && (sum_upper_w / sum_w >= p)) {
              
              status1 = true;
              idxquantile1_x = n_size_x - count - 1;
              sum_upper_w_temp = sum_upper_w;
              idxquantile2_x = idxquantile1_x;
              
            }
            
            n_obs += 1;
            
          }
          
          k_lower = sort_ix[std::max(0, n_size_x - count - 2)];
          
          is_na = arma_any_na_subset[k_lower] != 0 || std::isnan(x_subset[k_lower]);
          
          if (!is_na) {
            if (status1 && !status2) {
              
              status2 = true;
              idxquantile2_x = n_size_x - count - 1;
              
            }
          }
          
        }
        
        if (n_obs >= min_obs) {
          
          k = sort_ix[idxquantile1_x];
          
          // average if upper and lower weight is equal
          if (std::fabs(sum_upper_w_temp / sum_w - p) <= sqrt(arma::datum::eps)) {
            
            k_lower = sort_ix[std::max(0, idxquantile2_x - 1)];
            
            if ((arma_any_na_subset[k_lower] == 0) && !std::isnan(x_subset[k_lower])) {
              rcpp_quantile(i, j) = (x_subset[k] + x_subset[k_lower]) / 2;
            } else {
              rcpp_quantile(i, j) = x_subset[k];
            }
            
          } else {
            rcpp_quantile(i, j) = x_subset[k];
          }
          
        } else {
          rcpp_quantile(i, j) = NA_REAL;
        }
        
      } else {
        
        // can be either NA or NaN
        rcpp_quantile(i, j) = x(i, j);
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an online algorithm
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
      bool is_na = false;
      bool is_na_old = false;
      long double lambda = 0;
      long double lambda_sq = 0;
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
      
      for (int i = 0; i < n_rows_x; i++) {
        
        is_na = (arma_any_na[i] != 0) || std::isnan(x(i, j));
        
        if (i >= width) {
          is_na_old = (arma_any_na[i - width] != 0) || (std::isnan(x(i - width, j)));
        }
        
        update_n_obs(n_obs, is_na, is_na_old, i, width);
        
        if (width > 1) {
          
          lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed
          lambda_sq = lambda * lambda;
          
          if (!is_na) {
            
            w_new = arma_weights[n - 1];
            x_new = x(i, j);
            
          } else {
            
            w_new = 0;
            x_new = 0;
            
          }
          
          // expanding window
          if (i < width) {
            
            sum_w = lambda * sum_w + w_new;
            sum_x = lambda * sum_x + w_new * x_new;
            sumsq_w = lambda_sq * sumsq_w + w_new * w_new;
            
            if (center && (n_obs > 0)) {
              
              // compute the mean
              mean_prev_x = mean_x;
              mean_x = sum_x / sum_w;
              
            }
            
            // compute the sum of squares
            if (!is_na) {
              
              sumsq_x = lambda * sumsq_x +
                w_new * (x_new - mean_x) * (x_new - mean_prev_x);
              
            } else {
              sumsq_x = lambda * sumsq_x;
            }
            
            // rolling window
          } else {
            
            if (!is_na_old) {
              
              w_old = arma_weights[n - width];
              x_old = x(i - width, j);
              
            } else {
              
              w_old = 0;
              x_old = 0;
              
            }
            
            sum_w = lambda * sum_w + w_new - lambda * w_old;
            sum_x = lambda * sum_x + w_new * x_new - lambda * w_old * x_old;
            sumsq_w = lambda_sq * sumsq_w + w_new * w_new - lambda_sq * w_old * w_old;
            
            if (center) {
              
              // compute the mean
              mean_prev_x = mean_x;
              mean_x = sum_x / sum_w;
              
            }
            
            // compute the sum of squares
            if (!is_na && !is_na_old) {
              
              sumsq_x = lambda * sumsq_x +
                w_new * (x_new - mean_x) * (x_new - mean_prev_x) -
                lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
              
            } else if (!is_na && is_na_old) {
              
              sumsq_x = lambda * sumsq_x +
                w_new * (x_new - mean_x) * (x_new - mean_prev_x);
              
            } else if (is_na && !is_na_old) {
              
              sumsq_x = lambda * sumsq_x -
                lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
              
            } else if (is_na && is_na_old) {
              sumsq_x = lambda * sumsq_x;
            }
            
          }
          
          // don't compute if missing value and 'na_restore' argument is TRUE
          if (!na_restore || !std::isnan(x(i, j))) {
            
            if ((n_obs > 1) && (n_obs >= min_obs)) {
              arma_var(i, j) = sumsq_x / (sum_w - sumsq_w / sum_w);
            } else {
              arma_var(i, j) = NA_REAL;
            }
            
          } else {
            
            // can be either NA or NaN
            arma_var(i, j) = x(i, j);
            
          }
          
        } else {
          arma_var(i, j) = NA_REAL;
        }
        
      }
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollVarOfflineMat : public Worker {
  
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
  RollVarOfflineMat(const NumericMatrix x, const int n,
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
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if (!na_restore || !std::isnan(x(i, j))) {
        
        bool is_na = false;
        long double mean_x = 0;
        
        if (center) {
          
          long double sum_w = 0;
          long double sum_x = 0;
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          for (int count = 0; width > count && i >= count; count++) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
            is_na = (arma_any_na[i - count] != 0) || std::isnan(x(i - count, j));
            
            if (!is_na) {
              
              // compute the sum
              sum_w += arma_weights[n - count - 1];
              sum_x += arma_weights[n - count - 1] * x(i - count, j);
              
            }
            
          }
          
          // compute the mean
          mean_x = sum_x / sum_w;
          
        }
        
        int n_obs = 0;
        long double sum_w = 0;
        long double sumsq_w = 0;
        long double sumsq_x = 0;
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        for (int count = 0; width > count && i >= count; count++) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          is_na = (arma_any_na[i - count] != 0) || std::isnan(x(i - count, j));
          
          if (!is_na) {
            
            sum_w += arma_weights[n - count - 1];
            sumsq_w += arma_weights[n - count - 1] * arma_weights[n - count - 1];
            
            // compute the sum of squares with 'center' argument
            if (center) {
              
              sumsq_x += arma_weights[n - count - 1] *
                (x(i - count, j) - mean_x) * (x(i - count, j) - mean_x);
              
            } else if (!center) {
              
              sumsq_x += arma_weights[n - count - 1] *
                x(i - count, j) * x(i - count, j);
              
            }
            
            n_obs += 1;
            
          }
          
        }
        
        if ((n_obs > 1) && (n_obs >= min_obs)) {
          arma_var(i, j) = sumsq_x / (sum_w - sumsq_w / sum_w);
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

// 'Worker' function for computing the rolling statistic using an online algorithm
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
      long double var_x = 0;
      
      if (width > 1) {
        lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed
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
        
        // don't compute if missing value and 'na_restore' argument is TRUE
        if ((!na_restore) || (na_restore && !std::isnan(x(i, j)))) {
          
          if ((n_obs > 1) && (n_obs >= min_obs)) {
            
            if ((var_x < 0) || (sqrt(var_x) <= sqrt(arma::datum::eps))) {
              arma_sd(i, j) = 0;
            } else {
              arma_sd(i, j) = sqrt(var_x);
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

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollSdOfflineMat : public Worker {
  
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
  RollSdOfflineMat(const NumericMatrix x, const int n,
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
              
              // compute the sum
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
            
            // compute the sum of squares with 'center' argument
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
        
        if ((n_obs > 1) && (n_obs >= min_obs)) {
          arma_sd(i, j) = sqrt(sumsq_x / (sum_w - sumsq_w / sum_w));
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

// 'Worker' function for computing the rolling statistic using an online algorithm
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
        lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed
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

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollScaleOfflineMat : public Worker {
  
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
  RollScaleOfflineMat(const NumericMatrix x, const int n,
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
              
              // compute the sum
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
              
              // compute the sum of squares with 'center' argument
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

// 'Worker' function for computing the rolling statistic using an online algorithm
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
          lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed
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
              
              if ((n_obs > 1) && (n_obs >= min_obs)) {
                
                if (scale) {
                  
                  // don't compute if negative or sqrt is zero
                  if ((sumsq_x < 0) || (sumsq_y < 0)) {
                    arma_cov(j, k, i) = NA_REAL;
                  } else if ((sqrt(sumsq_x) <= sqrt(arma::datum::eps)) || (sqrt(sumsq_y) <= sqrt(arma::datum::eps))) {
                    arma_cov(j, k, i) = NA_REAL;
                  } else {
                    arma_cov(j, k, i) = sumsq_xy / (sqrt(sumsq_x) * sqrt(sumsq_y));
                  }
                  
                } else if (!scale) {
                  arma_cov(j, k, i) = sumsq_xy / (sum_w - sumsq_w / sum_w);
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

// 'Worker' function for computing the rolling statistic using an online algorithm
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
          lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed
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
              
              if ((n_obs > 1) && (n_obs >= min_obs)) {
                
                if (scale) {
                  
                  // don't compute if negative or sqrt is zero
                  if ((sumsq_x < 0) || (sumsq_y < 0)) {
                    arma_cov(j, k, i) = NA_REAL;
                  } else if ((sqrt(sumsq_x) <= sqrt(arma::datum::eps)) || (sqrt(sumsq_y) <= sqrt(arma::datum::eps))) {
                    arma_cov(j, k, i) = NA_REAL;
                  } else {
                    arma_cov(j, k, i) = sumsq_xy / (sqrt(sumsq_x) * sqrt(sumsq_y));
                  }
                  
                } else if (!scale) {
                  arma_cov(j, k, i) = sumsq_xy / (sum_w - sumsq_w / sum_w);
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

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollCovOfflineMatXX : public Worker {
  
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
  RollCovOfflineMatXX(const NumericMatrix x, const int n,
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
      
      long double sumsq_x = 0;
      long double sumsq_y = 0;
      long double mean_x = 0;
      long double mean_y = 0;
      
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
                  
                  // compute the sum
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
            
            // number of observations is either the window size or,
            // for partial results, the number of the current row
            while ((width > count) && (i >= count)) {
              
              // don't include if missing value and 'any_na' argument is 1
              // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
              if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j)) &&
                  !std::isnan(x(i - count, k))) {
                  
                  // compute the sum of squares with 'center' argument
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
              
              // compute the sum of squares with 'center' argument
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
          
          if ((n_obs > 1) && (n_obs >= min_obs)) {
            
            if (scale) {
              
              // don't compute if negative or sqrt is zero
              if ((sumsq_x < 0) || (sumsq_y < 0)) {
                arma_cov(j, k, i) = NA_REAL;
              } else if ((sqrt(sumsq_x) <= sqrt(arma::datum::eps)) || (sqrt(sumsq_y) <= sqrt(arma::datum::eps))) {
                arma_cov(j, k, i) = NA_REAL;
              } else {
                arma_cov(j, k, i) = sumsq_xy / (sqrt(sumsq_x) * sqrt(sumsq_y));
              }
              
            } else if (!scale) {
              arma_cov(j, k, i) = sumsq_xy / (sum_w - sumsq_w / sum_w);
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

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollCovOfflineMatXY : public Worker {
  
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
  RollCovOfflineMatXY(const NumericMatrix x, const NumericMatrix y,
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
      
      long double sumsq_x = 0;
      long double sumsq_y = 0;
      long double mean_x = 0;
      long double mean_y = 0;
      
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
                  
                  // compute the sum
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
            
            // number of observations is either the window size or,
            // for partial results, the number of the current row
            while ((width > count) && (i >= count)) {
              
              // don't include if missing value and 'any_na' argument is 1
              // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
              if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j)) &&
                  !std::isnan(y(i - count, k))) {
                  
                  // compute the sum of squares with 'center' argument
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
              
              // compute the sum of squares with 'center' argument
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
          
          if ((n_obs > 1) && (n_obs >= min_obs)) {
            
            if (scale) {
              
              // don't compute if negative or sqrt is zero
              if ((sumsq_x < 0) || (sumsq_y < 0)) {
                arma_cov(j, k, i) = NA_REAL;
              } else if ((sqrt(sumsq_x) <= sqrt(arma::datum::eps)) || (sqrt(sumsq_y) <= sqrt(arma::datum::eps))) {
                arma_cov(j, k, i) = NA_REAL;
              } else {
                arma_cov(j, k, i) = sumsq_xy / (sqrt(sumsq_x) * sqrt(sumsq_y));
              }
              
            } else if (!scale) {
              arma_cov(j, k, i) = sumsq_xy / (sum_w - sumsq_w / sum_w);
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

// 'Worker' function for computing the rolling statistic using an online algorithm
struct RollCrossProdOnlineMatXX : public Worker {
  
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
  arma::vec& arma_n_obs;        // destination (pass by reference)
  arma::vec& arma_sum_w;
  arma::mat& arma_mean;
  arma::cube& arma_cov;
  
  // initialize with source and destination
  RollCrossProdOnlineMatXX(const NumericMatrix x, const int n,
                           const int n_rows_xy, const int n_cols_x,
                           const int width, const arma::vec arma_weights,
                           const bool center, const bool scale,
                           const int min_obs, const arma::uvec arma_any_na,
                           const bool na_restore, arma::vec& arma_n_obs,
                           arma::vec& arma_sum_w, arma::mat& arma_mean,
                           arma::cube& arma_cov)
    : x(x), n(n),
      n_rows_xy(n_rows_xy), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      center(center), scale(scale),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), arma_n_obs(arma_n_obs),
      arma_sum_w(arma_sum_w), arma_mean(arma_mean),
      arma_cov(arma_cov) { }
  
  // function call operator that iterates by column
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      for (std::size_t k = 0; k <= j; k++) {
        
        int n_obs = 0;
        bool is_na = false;
        bool is_na_old = false;
        long double lambda = 0;
        // long double lambda_sq = 0;
        long double w_new = 0;
        long double w_old = 0;      
        long double x_new = 0;
        long double x_old = 0;
        long double y_new = 0;
        long double y_old = 0;
        long double sum_w = 0;
        long double sum_x = 0;
        long double sum_y = 0;
        // long double sumsq_w = 0;
        long double sumsq_x = 0;
        long double sumsq_y = 0;
        long double sumsq_xy = 0;
        long double mean_prev_x = 0;
        long double mean_prev_y = 0;
        long double mean_x = 0;
        long double mean_y = 0;
        
        for (int i = 0; i < n_rows_xy; i++) {
          
          is_na = (arma_any_na[i] != 0) || std::isnan(x(i, j)) || std::isnan(x(i, k));
          
          if (i >= width) {
            
            is_na_old = (arma_any_na[i - width] != 0) || std::isnan(x(i - width, j)) ||
              std::isnan(x(i - width, k));
            
          }
          
          update_n_obs(n_obs, is_na, is_na_old, i, width);
          
          if (width > 1) {
            
            lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed
            // lambda_sq = lambda * lambda;
            
            if (!is_na) {
              
              w_new = arma_weights[n - 1];
              x_new = x(i, j);
              y_new = x(i, k);
              
            } else {
              
              w_new = 0;
              x_new = 0;
              y_new = 0;
              
            }
            
            // expanding window
            if (i < width) {
              
              sum_w = lambda * sum_w + w_new;
              sum_x = lambda * sum_x + w_new * x_new;
              sum_y = lambda * sum_y + w_new * y_new;
              // sumsq_w = lambda_sq * sumsq_w + w_new * w_new;
              
              if (center && (n_obs > 0)) {
                
                // compute the mean
                mean_prev_x = mean_x;
                mean_prev_y = mean_y;
                mean_x = sum_x / sum_w;
                mean_y = sum_y / sum_w;
                
              }
              
              if (scale) {
                
                // compute the sum of squares
                if (!is_na) {
                  
                  sumsq_x = lambda * sumsq_x +
                    w_new * (x_new - mean_x) * (x_new - mean_prev_x);
                  sumsq_y = lambda * sumsq_y +
                    w_new * (y_new - mean_y) * (y_new - mean_prev_y);
                  
                } else {
                  
                  sumsq_x = lambda * sumsq_x;
                  sumsq_y = lambda * sumsq_y;
                  
                }
                
              }
              
              // compute the sum of squares
              if (!is_na) {
                
                sumsq_xy = lambda * sumsq_xy +
                  w_new * (x_new - mean_x) * (y_new - mean_prev_y);
                
              } else {
                sumsq_xy = lambda * sumsq_xy;
              }
              
              // rolling window
            } else {
              
              if (!is_na_old) {
                
                w_old = arma_weights[n - width];
                x_old = x(i - width, j);
                y_old = x(i - width, k);
                
              } else {
                
                w_old = 0;
                x_old = 0;
                y_old = 0;
                
              }
              
              sum_w = lambda * sum_w + w_new - lambda * w_old;
              sum_x = lambda * sum_x + w_new * x_new - lambda * w_old * x_old;
              sum_y = lambda * sum_y + w_new * y_new - lambda * w_old * y_old;
              // sumsq_w = lambda_sq * sumsq_w + w_new * w_new - lambda_sq * w_old * w_old;
              
              if (center && (n_obs > 0)) {
                
                // compute the mean
                mean_prev_x = mean_x;
                mean_prev_y = mean_y;
                mean_x = sum_x / sum_w;
                mean_y = sum_y / sum_w;
                
              }
              
              if (scale) {
                
                // compute the sum of squares
                if (!is_na && !is_na_old) {
                  
                  sumsq_x = lambda * sumsq_x +
                    w_new * (x_new - mean_x) * (x_new - mean_prev_x) -
                    lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
                  sumsq_y = lambda * sumsq_y +
                    w_new * (y_new - mean_y) * (y_new - mean_prev_y) -
                    lambda * w_old * (y_old - mean_y) * (y_old - mean_prev_y);
                  
                } else if (!is_na && is_na_old) {
                  
                  sumsq_x = lambda * sumsq_x +
                    w_new * (x_new - mean_x) * (x_new - mean_prev_x);
                  sumsq_y = lambda * sumsq_y +
                    w_new * (y_new - mean_y) * (y_new - mean_prev_y);
                  
                } else if (is_na && !is_na_old) {
                  
                  sumsq_x = lambda * sumsq_x -
                    lambda * w_old * (x_old - mean_x) * (x_old - mean_prev_x);
                  sumsq_y = lambda * sumsq_y -
                    lambda * w_old * (y_old - mean_y) * (y_old - mean_prev_y);
                  
                } else if (is_na && is_na_old) {
                  
                  sumsq_x = lambda * sumsq_x;
                  sumsq_y = lambda * sumsq_y;
                  
                }
                
              }
              
              // compute the sum of squares
              if (!is_na && !is_na_old) {
                
                sumsq_xy = lambda * sumsq_xy +
                  w_new * (x_new - mean_x) * (y_new - mean_prev_y) -
                  lambda * w_old * (x_old - mean_x) * (y_old - mean_prev_y);
                
              } else if (!is_na && is_na_old) {
                
                sumsq_xy = lambda * sumsq_xy +
                  w_new * (x_new - mean_x) * (y_new - mean_prev_y);
                
              } else if (is_na && !is_na_old) {
                
                sumsq_xy = lambda * sumsq_xy -
                  lambda * w_old * (x_old - mean_x) * (y_old - mean_prev_y);
                
              } else if (is_na && is_na_old) {
                sumsq_xy = lambda * sumsq_xy;
              }
              
            }
            
          } else {
            
            lambda = arma_weights[n - 1];
            
            if (!is_na) {
              
              w_new = arma_weights[n - 1];
              x_new = x(i, j);
              y_new = x(i, k);
              
            } else {
              
              w_new = 0;
              x_new = 0;
              y_new = 0;
              
            }
            
            sum_w = w_new;
            sum_x = w_new * x_new;
            sum_y = w_new * y_new;
            // sumsq_w = w_new * w_new;
            
            if (center && (n_obs > 0)) {
              
              // compute the mean
              mean_x = sum_x / sum_w;
              mean_y = sum_y / sum_w;
              
            }
            
            if (scale) {
              
              // compute the sum of squares
              if (!is_na) {
                
                sumsq_x = w_new * (x_new - mean_x) * (x_new - mean_x);
                sumsq_y = w_new * (y_new - mean_y) * (y_new - mean_y);
                
              } else {
                
                sumsq_x = 0;
                sumsq_y = 0;
                
              }
              
            }
            
            // compute the sum of squares
            if (!is_na) {
              sumsq_xy = w_new * (x_new - mean_x) * (y_new - mean_y);
            } else {
              sumsq_xy = 0;
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
          if (!na_restore || (!std::isnan(x(i, j)) && !std::isnan(x(i, k)))) {
            
            // if ((n_obs > 1) && (n_obs >= min_obs)) {
            if (n_obs >= min_obs) {
              
              if (scale) {
                
                // don't compute if negative or sqrt is zero
                if ((sumsq_x < 0) || (sumsq_y < 0)) {
                  arma_cov(j, k, i) = NA_REAL;
                } else if ((sqrt(sumsq_x) <= sqrt(arma::datum::eps)) || (sqrt(sumsq_y) <= sqrt(arma::datum::eps))) {
                  arma_cov(j, k, i) = NA_REAL;
                } else {
                  arma_cov(j, k, i) = sumsq_xy / (sqrt(sumsq_x) * sqrt(sumsq_y));
                }
                
              } else if (!scale) {
                arma_cov(j, k, i) = sumsq_xy; // / (sum_w - sumsq_w / sum_w);
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

// 'Worker' function for computing the rolling statistic using an online algorithm
struct RollCrossProdOnlineMatXY : public Worker {
  
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
  RollCrossProdOnlineMatXY(const NumericMatrix x, const NumericMatrix y,
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
        // long double sumsq_w = 0;
        long double sumsq_x = 0;
        long double sumsq_y = 0;
        long double sumsq_xy = 0;
        long double mean_prev_x = 0;
        long double mean_prev_y = 0;
        long double mean_x = 0;
        long double mean_y = 0;
        
        if (width > 1) {
          lambda = arma_weights[n - 2] / arma_weights[n - 1]; // check already passed
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
            
            if (width > 1) {
              
              sum_w = lambda * sum_w + w_new;
              sum_x = lambda * sum_x + w_new * x_new;
              sum_y = lambda * sum_y + w_new * y_new;
              // sumsq_w = pow(lambda, (long double)2.0) * sumsq_w + pow(w_new, (long double)2.0);
              
            } else {
              
              sum_w = w_new;
              sum_x = w_new * x_new;
              sum_y = w_new * y_new;
              // sumsq_w = pow(w_new, (long double)2.0);
              
            }
            
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
            
            if (width > 1) {
              
              sum_w = lambda * sum_w + w_new - lambda * w_old;
              sum_x = lambda * sum_x + w_new * x_new - lambda * w_old * x_old;
              sum_y = lambda * sum_y + w_new * y_new - lambda * w_old * y_old;
              // sumsq_w = pow(lambda, (long double)2.0) * sumsq_w +
              //   pow(w_new, (long double)2.0) - pow(lambda * w_old, (long double)2.0);
              
            } else {
              
              sum_w = w_new;
              sum_x = w_new * x_new;
              sum_y = w_new * y_new;
              // sumsq_w = pow(w_new, (long double)2.0);
              
            }
            
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
              
              // if ((n_obs > 1) && (n_obs >= min_obs)) {
              if (n_obs >= min_obs) {
                
                if (scale) {
                  
                  // don't compute if negative or sqrt is zero
                  if ((sumsq_x < 0) || (sumsq_y < 0)) {
                    arma_cov(j, k, i) = NA_REAL;
                  } else if ((sqrt(sumsq_x) <= sqrt(arma::datum::eps)) || (sqrt(sumsq_y) <= sqrt(arma::datum::eps))) {
                    arma_cov(j, k, i) = NA_REAL;
                  } else {
                    arma_cov(j, k, i) = sumsq_xy / (sqrt(sumsq_x) * sqrt(sumsq_y));
                  }
                  
                } else if (!scale) {
                  arma_cov(j, k, i) = sumsq_xy; // / (sum_w - sumsq_w / sum_w);
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

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollCrossProdOfflineMatXX : public Worker {
  
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
  arma::vec& arma_n_obs;        // destination (pass by reference)
  arma::vec& arma_sum_w;
  arma::mat& arma_mean;
  arma::cube& arma_cov;
  
  // initialize with source and destination
  RollCrossProdOfflineMatXX(const NumericMatrix x, const int n,
                            const int n_rows_xy, const int n_cols_x,
                            const int width, const arma::vec arma_weights,
                            const bool center, const bool scale, 
                            const int min_obs, const arma::uvec arma_any_na,
                            const bool na_restore, arma::vec& arma_n_obs,
                            arma::vec& arma_sum_w, arma::mat& arma_mean,
                            arma::cube& arma_cov)
    : x(x), n(n),
      n_rows_xy(n_rows_xy), n_cols_x(n_cols_x),
      width(width), arma_weights(arma_weights),
      center(center), scale(scale),
      min_obs(min_obs), arma_any_na(arma_any_na),
      na_restore(na_restore), arma_n_obs(arma_n_obs),
      arma_sum_w(arma_sum_w), arma_mean(arma_mean),
      arma_cov(arma_cov) { }
  
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
      
      // don't compute if missing value and 'na_restore' argument is TRUE
      if (!na_restore || (!std::isnan(x(i, j)) && !std::isnan(x(i, k)))) {
        
        bool is_na = false;
        long double sumsq_x = 0;
        long double sumsq_y = 0;
        long double mean_x = 0;
        long double mean_y = 0;
        
        if (center) {
          
          long double sum_w = 0;
          long double sum_x = 0;
          long double sum_y = 0;
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          for (int count = 0; width > count && i >= count; count++) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
            is_na = (arma_any_na[i - count] != 0) || std::isnan(x(i - count, j)) ||
              std::isnan(x(i - count, k));
            
            if (!is_na) {
              
              // compute the sum
              sum_w += arma_weights[n - count - 1];
              sum_x += arma_weights[n - count - 1] * x(i - count, j);
              sum_y += arma_weights[n - count - 1] * x(i - count, k);
              
            }
            
          }
          
          // compute the mean
          mean_x = sum_x / sum_w;
          mean_y = sum_y / sum_w;
          
        }
        
        if (scale) {
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          for (int count = 0; width > count && i >= count; count++) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
            is_na = (arma_any_na[i - count] != 0) || std::isnan(x(i - count, j)) ||
              std::isnan(x(i - count, k));
            
            if (!is_na) {
              
              // compute the sum of squares with 'center' argument
              if (center) {
                
                sumsq_x += arma_weights[n - count - 1] *
                  (x(i - count, j) - mean_x) * (x(i - count, j) - mean_x);
                sumsq_y += arma_weights[n - count - 1] *
                  (x(i - count, k) - mean_y) * (x(i - count, k) - mean_y);
                
              } else if (!center) {
                
                sumsq_x += arma_weights[n - count - 1] *
                  x(i - count, j) * x(i - count, j);
                sumsq_y += arma_weights[n - count - 1] *
                  x(i - count, k) * x(i - count, k);
                
              }
              
            }
            
          }
          
        }
        
        int n_obs = 0;
        long double sum_w = 0;
        // long double sumsq_w = 0;
        long double sumsq_xy = 0;
        
        // number of observations is either the window size or,
        // for partial results, the number of the current row
        for (int count = 0; width > count && i >= count; count++) {
          
          // don't include if missing value and 'any_na' argument is 1
          // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
          is_na = (arma_any_na[i - count] != 0) || std::isnan(x(i - count, j)) ||
            std::isnan(x(i - count, k));
          
          if (!is_na) {
            
            sum_w += arma_weights[n - count - 1];
            // sumsq_w += arma_weights[n - count - 1] * arma_weights[n - count - 1];
            
            // compute the sum of squares with 'center' argument
            if (center) {
              
              sumsq_xy += arma_weights[n - count - 1] * 
                (x(i - count, j) - mean_x) * (x(i - count, k) - mean_y);
              
            } else if (!center) {
              
              sumsq_xy += arma_weights[n - count - 1] * 
                x(i - count, j) * x(i - count, k);
              
            }
            
            n_obs += 1;
            
          }
          
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
        
        // if ((n_obs > 1) && (n_obs >= min_obs)) {
        if (n_obs >= min_obs) {
          
          if (scale) {
            
            // don't compute if negative or sqrt is zero
            if ((sumsq_x < 0) || (sumsq_y < 0)) {
              arma_cov(j, k, i) = NA_REAL;
            } else if ((sqrt(sumsq_x) <= sqrt(arma::datum::eps)) || (sqrt(sumsq_y) <= sqrt(arma::datum::eps))) {
              arma_cov(j, k, i) = NA_REAL;
            } else {
              arma_cov(j, k, i) = sumsq_xy / (sqrt(sumsq_x) * sqrt(sumsq_y));
            }
            
          } else if (!scale) {
            arma_cov(j, k, i) = sumsq_xy; // / (sum_w - sumsq_w / sum_w);
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

// 'Worker' function for computing the rolling statistic using an offline algorithm
struct RollCrossProdOfflineMatXY : public Worker {
  
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
  RollCrossProdOfflineMatXY(const NumericMatrix x, const NumericMatrix y,
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
      
      long double sumsq_x = 0;
      long double sumsq_y = 0;
      long double mean_x = 0;
      long double mean_y = 0;
      
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
                  
                  // compute the sum
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
            
            // number of observations is either the window size or,
            // for partial results, the number of the current row
            while ((width > count) && (i >= count)) {
              
              // don't include if missing value and 'any_na' argument is 1
              // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
              if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j)) &&
                  !std::isnan(y(i - count, k))) {
                  
                  // compute the sum of squares with 'center' argument
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
            
          }
          
          int count = 0;
          int n_obs = 0;
          long double sum_w = 0;
          // long double sumsq_w = 0;
          long double sumsq_xy = 0;
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is FALSE
            if ((arma_any_na[i - count] == 0) && !std::isnan(x(i - count, j)) &&
                !std::isnan(y(i - count, k))) {
                
                sum_w += arma_weights[n - count - 1];
              // sumsq_w += pow(arma_weights[n - count - 1], 2.0);
              
              // compute the sum of squares with 'center' argument
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
          
          // if ((n_obs > 1) && (n_obs >= min_obs)) {
          if (n_obs >= min_obs) {
            
            if (scale) {
              
              // don't compute if negative or sqrt is zero
              if ((sumsq_x < 0) || (sumsq_y < 0)) {
                arma_cov(j, k, i) = NA_REAL;
              } else if ((sqrt(sumsq_x) <= sqrt(arma::datum::eps)) || (sqrt(sumsq_y) <= sqrt(arma::datum::eps))) {
                arma_cov(j, k, i) = NA_REAL;
              } else {
                arma_cov(j, k, i) = sumsq_xy / (sqrt(sumsq_x) * sqrt(sumsq_y));
              }
              
            } else if (!scale) {
              arma_cov(j, k, i) = sumsq_xy; // / (sum_w - sumsq_w / sum_w);
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

// 'Worker' function for computing the rolling statistic using a standard algorithm
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
            
            // standard errors
            long double var_resid = (1 - arma_rsq[i]) * var_y / df_resid;
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

// 'Worker' function for computing the rolling statistic using a standard algorithm
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
            
            // standard errors
            long double var_resid = (1 - arma_rsq[i]) * var_y / df_resid;
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

}

#endif