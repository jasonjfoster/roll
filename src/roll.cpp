// Google C++ Style Guide: https://google.github.io/styleguide/cppguide.html

#define ARMA_DONT_PRINT_ERRORS

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

List dimnames_ols(const List& input, const int& n_cols) {
  
  if (input.size() > 1) {
    
    CharacterVector dimnames_cols = input[1];
    int n = dimnames_cols.size() + 1;
    CharacterVector result(n);
    
    result(0) = "(Intercept)";
    std::copy(dimnames_cols.begin(), dimnames_cols.end(), result.begin() + 1);
    
    return(List::create(input[0], result));
    
  } else {
    
    std::string x = "x";
    CharacterVector result(n_cols + 1);
    result(0) = "(Intercept)";
    
    for (int i = 1; i < n_cols + 1; i++) {
      
      char n[sizeof(i)];
      
      sprintf(n, "%d", i);
      
      result[i] = x + n;
      
    }
    
    return(List::create(R_NilValue, result));
    
  }
  
}

CharacterVector dimnames_pc(const int& n_cols) {
  
  std::string x = "PC";
  CharacterVector result(n_cols);
  
  for (int i = 0; i < n_cols; i++) {
    
    char n[sizeof(i + 1)];
    
    sprintf(n, "%d", i + 1);
    
    result[i] = x + n;
    
  }
  
  return result;
  
}

struct AnyNaRows : public Worker {
  
  const RMatrix<double> data;
  const int n_rows;
  const int n_cols;
  arma::uvec& result;
  
  AnyNaRows(const NumericMatrix data, const int n_rows,
            const int n_cols, arma::uvec& result)
    : data(data), n_rows(n_rows),
      n_cols(n_cols), result(result) { }
  
  void operator()(std::size_t begin_row, std::size_t end_row) {
    for (std::size_t i = begin_row; i < end_row; i++) {
      int any_na = 0;
      int j = 0;
      while ((any_na == 0) && (j < n_cols)) {
        if (std::isnan(data(i, j))) {
          any_na = 1;
        }
        j += 1;
      }
      result[i] = any_na;
    }
  }
  
};

arma::uvec any_na(const NumericMatrix& data) {
  
  int n_rows = data.nrow();
  int n_cols = data.ncol();
  arma::uvec result(n_rows);
  
  AnyNaRows any_na_rows(data, n_rows, n_cols, result);
  parallelFor(0, n_rows, any_na_rows); 
  
  return result;
  
}

void check_ols(const int& x_n_rows, const int& y_n_rows) {
  
  if (x_n_rows != y_n_rows) {
    stop("number of rows in 'x' must equal the number of rows in 'y'");
  }
  
}

void check_vif(const int& x_n_rows) {
  
  if (x_n_rows == 1) {
    stop("number of columns in 'x' must be greater than one");
  }
  
}

void check_width(const int& width, const int& n_rows) {
  
  if ((width < 1) || (width > n_rows)) {
    stop("value of 'width' must be between one and number of rows in 'data'");
  }
  
}

void check_comps(const arma::uvec& comps, const unsigned int& n_cols) {
  
  if (comps.max() > n_cols) {
    stop("maximum value of 'comps' must be less than or equal to number of columns in 'x'");
  }
  
  if (comps.min() < 1) {
    stop("minimum value of 'comps' must be greater than or equal to one");
  }
  
  if (comps.size() > n_cols) {
    stop("length of 'comps' must be less than or equal to number of columns in 'x'");
  }
  
}

// this is the same as check_comps except uses < instead of <=
void check_comps_vif(const arma::uvec& comps, const unsigned int& n_cols) {
  
  if (comps.max() > n_cols) {
    stop("maximum value of 'comps' must be less than number of columns in 'x'");
  }
  
  if (comps.min() < 1) {
    stop("minimum value of 'comps' must be greater than or equal to one");
  }
  
  if (comps.size() > n_cols) {
    stop("length of 'comps' must be less than number of columns in 'x'");
  }
  
}

void check_weights(const arma::vec& weights, const unsigned int& width) {
  
  if (weights.size() != width) {
    stop("length of 'weights' must be equal to 'width'");
  }
  
}

void check_min_obs(const int& min_obs, const int& width) {
  
  if ((min_obs < 1) || (min_obs > width)) {
    stop("value of 'min_obs' must be between one and 'width'");
  }
  
}

arma::uvec seq(const int& size) {
  
  arma::uvec result(size);
  
  for (int i = 0; i < size; i++) {
    result[i] = i;
  }
  
  return(result);
  
}

// 'Worker' function for computing rolling means
struct RollMeanRows : public Worker {
  
  const RMatrix<double> data;   // source
  const int n_rows;
  const int n_cols;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_center;       // destination (pass by reference)
  
  // initialize with source and destination
  RollMeanRows(const NumericMatrix data, const int n_rows,
               const int n_cols, const int width,
               const arma::vec arma_weights, const int min_obs,
               const arma::uvec arma_any_na, const bool na_restore,
               arma::mat& arma_center)
    : data(data), n_rows(n_rows),
      n_cols(n_cols), width(width),
      arma_weights(arma_weights), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      arma_center(arma_center) { }
  
  // function call operator that iterates by row 
  void operator()(std::size_t begin_row, std::size_t end_row) {
    for (std::size_t i = begin_row; i < end_row; i++) {
      for (int j = 0; j < n_cols; j++) {
        
        int count = 0;
        int n_obs = 0;
        double sum_data = 0;
        double sum_weights = 0;
        
        // don't compute if missing value and 'na_restore' argument is true
        if ((!na_restore) || (na_restore && !std::isnan(data(i, j)))) {
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= (unsigned)count)) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is false
            if ((arma_any_na[i - count] == 0) && !std::isnan(data(i - count, j))) {
              sum_data += data(i - count, j) * arma_weights[width - count - 1];
              sum_weights += arma_weights[width - count - 1];
              n_obs += 1;
            }
            
            count += 1;
            
          }
          
          // compute the mean
          if (n_obs >= min_obs) {
            arma_center(i, j) = sum_data / sum_weights;
          } else {
            arma_center(i, j) = NA_REAL;
          }
          
        } else {
          arma_center(i, j) = NA_REAL;
        }
        
      }
    }
  }
  
};

// 'Worker' function for computing rolling means
struct RollMeanCols : public Worker {
  
  const RMatrix<double> data;   // source
  const int n_rows;
  const int n_cols;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_center;       // destination (pass by reference)
  
  // initialize with source and destination
  RollMeanCols(const NumericMatrix data, const int n_rows,
               const int n_cols, const int width,
               const arma::vec arma_weights, const int min_obs,
               const arma::uvec arma_any_na, const bool na_restore,
               arma::mat& arma_center)
    : data(data), n_rows(n_rows),
      n_cols(n_cols), width(width),
      arma_weights(arma_weights), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      arma_center(arma_center) { }
  
  // function call operator that iterates by column 
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      for (int i = 0; i < n_rows; i++) {
        
        int count = 0;
        int n_obs = 0;
        double sum_data = 0;
        double sum_weights = 0;
        
        // don't compute if missing value and 'na_restore' argument is true
        if ((!na_restore) || (na_restore && !std::isnan(data(i, j)))) {
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is false
            if ((arma_any_na[i - count] == 0) && !std::isnan(data(i - count, j))) {
              
              // compute the rolling sum
              sum_data += data(i - count, j) * arma_weights[width - count - 1];
              sum_weights += arma_weights[width - count - 1];
              n_obs += 1;
              
            }
            
            count += 1;
            
          }
          
          // compute the mean
          if (n_obs >= min_obs) {
            arma_center(i, j) = sum_data / sum_weights;
          } else {
            arma_center(i, j) = NA_REAL;
          }
          
        } else {
          arma_center(i, j) = NA_REAL;
        }
        
      }
    }
  }
  
};

// [[Rcpp::export]]
NumericMatrix roll_mean(const NumericMatrix& data, const int& width,
                        const arma::vec& weights, const int& min_obs,
                        const bool& complete_obs, const bool& na_restore,
                        const std::string& parallel_for) {
  
  int n_rows = data.nrow();
  int n_cols = data.ncol();
  arma::uvec arma_any_na(n_rows);
  arma::mat arma_center(n_rows, n_cols);
  
  // check 'width' argument for errors
  check_width(width, n_rows);
  
  // default 'weights' argument is equal-weighted,
  // otherwise check argument for errors
  check_weights(weights, width);
  
  // default 'min_obs' argument is 'width' (equivalent to 'na.rm = FALSE'),
  // otherwise check argument for errors
  check_min_obs(min_obs, width);
  
  // default 'complete_obs' argument is 'false' (equivalent to 'pairwise'),
  // otherwise check argument for errors
  if (complete_obs) {
    arma_any_na = any_na(data);
  } else {
    arma_any_na.fill(0);
  }
  
  // compute rolling means
  if (parallel_for == "rows") {
    RollMeanRows roll_mean_rows(data, n_rows, n_cols, width, weights,
                                min_obs, arma_any_na, na_restore, arma_center);
    parallelFor(0, n_rows, roll_mean_rows);
  } else if (parallel_for == "cols") {
    RollMeanCols roll_mean_cols(data, n_rows, n_cols, width, weights,
                                min_obs, arma_any_na, na_restore, arma_center);
    parallelFor(0, n_cols, roll_mean_cols);
  }
  
  // create and return a matrix or xts object
  NumericMatrix result(wrap(arma_center));
  List dimnames = data.attr("dimnames");
  result.attr("dimnames") = dimnames;
  result.attr("index") = data.attr("index");
  result.attr(".indexCLASS") = data.attr(".indexCLASS");
  result.attr(".indexTZ") = data.attr(".indexTZ");
  result.attr("tclass") = data.attr("tclass");
  result.attr("tzone") = data.attr("tzone");
  result.attr("class") = data.attr("class");
  
  return result;
  
}

// 'Worker' function for computing rolling variances
struct RollVarRows : public Worker {
  
  const RMatrix<double> data;   // source
  const int n_rows;
  const int n_cols;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const arma::mat arma_center;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_scale;        // destination (pass by reference)
  
  // initialize with source and destination
  RollVarRows(const NumericMatrix data, const int n_rows,
              const int n_cols, const int width,
              const arma::vec arma_weights, const bool center,
              const arma::mat arma_center, const int min_obs,
              const arma::uvec arma_any_na, const bool na_restore,
              arma::mat& arma_scale)
    : data(data), n_rows(n_rows),
      n_cols(n_cols), width(width),
      arma_weights(arma_weights), center(center),
      arma_center(arma_center), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      arma_scale(arma_scale) { }
  
  // function call operator that iterates by row
  void operator()(std::size_t begin_row, std::size_t end_row) {
    for (std::size_t i = begin_row; i < end_row; i++) {
      for (int j = 0; j < n_cols; j++) {
        
        int count = 0;
        int n_obs = 0;
        double sum_data = 0;
        double sum_weights = 0;
        double sum_weights_sq = 0;
        
        // don't compute if missing value and 'na_restore' argument is true
        if ((!na_restore) || (na_restore && !std::isnan(data(i, j)))) {
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= (unsigned)count)) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is false
            if ((arma_any_na[i - count] == 0) && !std::isnan(data(i - count, j))) {
              
              // compute the rolling sum of squares with 'center' argument
              if (center) {
                sum_data +=
                  pow(data(i - count, j) - arma_center(i, j), 2.0) *
                  arma_weights[width - count - 1];
              } else if (!center) {
                sum_data += pow(data(i - count, j), 2.0) *
                  arma_weights[width - count - 1];
              }
              
              sum_weights += arma_weights[width - count - 1];
              sum_weights_sq += pow(arma_weights[width - count - 1], 2.0);
              n_obs += 1;
              
            }
            
            count += 1;
            
          }
          
          // compute the unbiased estimate of variance
          if (n_obs >= min_obs) {
            arma_scale(i, j) = ((sum_data / sum_weights) /
                                  (1 - (sum_weights_sq / pow(sum_weights, 2.0))));
          } else {
            arma_scale(i, j) = NA_REAL;
          }
          
        } else {
          arma_scale(i, j) = NA_REAL;
        }
        
      }
    }
  }
  
};

// 'Worker' function for computing rolling variances
struct RollVarCols : public Worker {
  
  const RMatrix<double> data;   // source
  const int n_rows;
  const int n_cols;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const arma::mat arma_center;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_scale;        // destination (pass by reference)
  
  // initialize with source and destination
  RollVarCols(const NumericMatrix data, const int n_rows,
              const int n_cols, const int width,
              const arma::vec arma_weights, const bool center,
              const arma::mat arma_center, const int min_obs,
              const arma::uvec arma_any_na, const bool na_restore,
              arma::mat& arma_scale)
    : data(data), n_rows(n_rows),
      n_cols(n_cols), width(width),
      arma_weights(arma_weights), center(center),
      arma_center(arma_center), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      arma_scale(arma_scale) { }
  
  // function call operator that iterates by row
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      for (int i = 0; i < n_rows; i++) {
        
        int count = 0;
        int n_obs = 0;
        double sum_data = 0;
        double sum_weights = 0;
        double sum_weights_sq = 0;
        
        // don't compute if missing value and 'na_restore' argument is true
        if ((!na_restore) || (na_restore && !std::isnan(data(i, j)))) {
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is false
            if ((arma_any_na[i - count] == 0) && !std::isnan(data(i - count, j))) {
              
              // compute the rolling sum of squares with 'center' argument
              if (center) {
                sum_data +=
                  pow(data(i - count, j) - arma_center(i, j), 2.0) *
                  arma_weights[width - count - 1];
              } else if (!center) {
                sum_data +=
                  pow(data(i - count, j), 2.0) *
                  arma_weights[width - count - 1];
              }
              
              sum_weights += arma_weights[width - count - 1];
              sum_weights_sq += pow(arma_weights[width - count - 1], 2.0);
              n_obs += 1;
              
            }
            
            count += 1;
            
          }
          
          // compute the unbiased estimate of variance
          if (n_obs >= min_obs) {
            arma_scale(i, j) = ((sum_data / sum_weights) /
                                  (1 - (sum_weights_sq / pow(sum_weights, 2.0))));
          } else {
            arma_scale(i, j) = NA_REAL;
          }
          
        } else {
          arma_scale(i, j) = NA_REAL;
        }
        
      }
    }
  }
  
};

// [[Rcpp::export]]
NumericMatrix roll_var(const NumericMatrix& data, const int& width,
                       const arma::vec& weights, const bool& center,
                       const int& min_obs, const bool& complete_obs,
                       const bool& na_restore, const std::string& parallel_for) {
  
  int n_rows = data.nrow();
  int n_cols = data.ncol();
  arma::uvec arma_any_na(n_rows);
  arma::mat arma_center(n_rows, n_cols);
  arma::mat arma_scale(n_rows, n_cols);
  
  // check 'width' argument for errors
  check_width(width, n_rows);
  
  // default 'weights' argument is equal-weighted,
  // otherwise check argument for errors
  check_weights(weights, width);
  
  // default 'min_obs' argument is 'width' (equivalent to 'na.rm = FALSE'),
  // otherwise check argument for errors
  check_min_obs(min_obs, width);
  
  // default 'complete_obs' argument is 'false' (equivalent to 'pairwise'),
  // otherwise check argument for errors
  if (complete_obs) {
    arma_any_na = any_na(data);
  } else {
    arma_any_na.fill(0);
  }
  
  // default 'center' argument subtracts mean of each variable,
  // otherwise no scaling is done
  if (center) {
    if (parallel_for == "rows") {
      RollMeanRows roll_mean_rows(data, n_rows, n_cols, width, weights,
                                  min_obs, arma_any_na, na_restore,
                                  arma_center);
      parallelFor(0, n_rows, roll_mean_rows);
    } else if (parallel_for == "cols") {
      RollMeanCols roll_mean_cols(data, n_rows, n_cols, width, weights,
                                  min_obs, arma_any_na, na_restore,
                                  arma_center);
      parallelFor(0, n_cols, roll_mean_cols);
    }
  }
  
  // compute rolling variances
  if (parallel_for == "rows") {
    RollVarRows roll_var_rows(data, n_rows, n_cols, width, weights,
                              center, arma_center,
                              min_obs, arma_any_na, na_restore,
                              arma_scale);
    parallelFor(0, n_rows, roll_var_rows); 
  } else if (parallel_for == "cols") {
    RollVarCols roll_var_cols(data, n_rows, n_cols, width, weights,
                              center, arma_center,
                              min_obs, arma_any_na, na_restore,
                              arma_scale);
    parallelFor(0, n_cols, roll_var_cols);   
  }
  
  // create and return a matrix or xts object
  NumericMatrix result(wrap(arma_scale));
  List dimnames = data.attr("dimnames");
  result.attr("dimnames") = dimnames;
  result.attr("index") = data.attr("index");
  result.attr(".indexCLASS") = data.attr(".indexCLASS");
  result.attr(".indexTZ") = data.attr(".indexTZ");
  result.attr("tclass") = data.attr("tclass");
  result.attr("tzone") = data.attr("tzone");
  result.attr("class") = data.attr("class");
  
  return result;
  
}

// 'Worker' function for computing rolling standard deviations
struct RollSdRows : public Worker {
  
  const RMatrix<double> data;   // source
  const int n_rows;
  const int n_cols;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const arma::mat arma_center;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_scale;        // destination (pass by reference)
  
  // initialize with source and destination
  RollSdRows(const NumericMatrix data, const int n_rows,
             const int n_cols, const int width,
             const arma::vec arma_weights, const bool center,
             const arma::mat arma_center, const int min_obs,
             const arma::uvec arma_any_na, const bool na_restore,
             arma::mat& arma_scale)
    : data(data), n_rows(n_rows),
      n_cols(n_cols), width(width),
      arma_weights(arma_weights), center(center),
      arma_center(arma_center), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      arma_scale(arma_scale) { }
  
  // function call operator that iterates by row
  void operator()(std::size_t begin_row, std::size_t end_row) {
    for (std::size_t i = begin_row; i < end_row; i++) {
      for (int j = 0; j < n_cols; j++) {
        
        int count = 0;
        int n_obs = 0;
        double sum_data = 0;
        double sum_weights = 0;
        double sum_weights_sq = 0;
        
        // don't compute if missing value and 'na_restore' argument is true
        if ((!na_restore) || (na_restore && !std::isnan(data(i, j)))) {
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= (unsigned)count)) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is false
            if ((arma_any_na[i - count] == 0) && !std::isnan(data(i - count, j))) {
              
              // compute the rolling sum of squares with 'center' argument
              if (center) {
                sum_data +=
                  pow(data(i - count, j) - arma_center(i, j), 2.0) *
                  arma_weights[width - count - 1];
              } else if (!center) {
                sum_data += pow(data(i - count, j), 2.0) *
                  arma_weights[width - count - 1];
              }
              
              sum_weights += arma_weights[width - count - 1];
              sum_weights_sq += pow(arma_weights[width - count - 1], 2.0);
              n_obs += 1;
              
            }
            
            count += 1;
            
          }
          
          // compute the unbiased estimate of standard deviation
          if (n_obs >= min_obs) {
            arma_scale(i, j) = sqrt((sum_data / sum_weights) /
                                      (1 - (sum_weights_sq / pow(sum_weights, 2.0))));
          } else {
            arma_scale(i, j) = NA_REAL;
          }
          
        } else {
          arma_scale(i, j) = NA_REAL;
        }
        
      }
    }
  }
  
};

// 'Worker' function for computing rolling standard deviations
struct RollSdCols : public Worker {
  
  const RMatrix<double> data;   // source
  const int n_rows;
  const int n_cols;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const arma::mat arma_center;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_scale;        // destination (pass by reference)
  
  // initialize with source and destination
  RollSdCols(const NumericMatrix data, const int n_rows,
             const int n_cols, const int width,
             const arma::vec arma_weights, const bool center,
             const arma::mat arma_center, const int min_obs,
             const arma::uvec arma_any_na, const bool na_restore,
             arma::mat& arma_scale)
    : data(data), n_rows(n_rows),
      n_cols(n_cols), width(width),
      arma_weights(arma_weights), center(center),
      arma_center(arma_center), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      arma_scale(arma_scale) { }
  
  // function call operator that iterates by row
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      for (int i = 0; i < n_rows; i++) {
        
        int count = 0;
        int n_obs = 0;
        double sum_data = 0;
        double sum_weights = 0;
        double sum_weights_sq = 0;
        
        // don't compute if missing value and 'na_restore' argument is true
        if ((!na_restore) || (na_restore && !std::isnan(data(i, j)))) {
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is false
            if ((arma_any_na[i - count] == 0) && !std::isnan(data(i - count, j))) {
              
              // compute the rolling sum of squares with 'center' argument
              if (center) {
                sum_data +=
                  pow(data(i - count, j) - arma_center(i, j), 2.0) *
                  arma_weights[width - count - 1];
              } else if (!center) {
                sum_data +=
                  pow(data(i - count, j), 2.0) *
                  arma_weights[width - count - 1];
              }
              
              sum_weights += arma_weights[width - count - 1];
              sum_weights_sq += pow(arma_weights[width - count - 1], 2.0);
              n_obs += 1;
              
            }
            
            count += 1;
            
          }
          
          // compute the unbiased estimate of standard deviation
          if (n_obs >= min_obs) {
            arma_scale(i, j) = sqrt((sum_data / sum_weights) /
                                      (1 - (sum_weights_sq / pow(sum_weights, 2.0))));
          } else {
            arma_scale(i, j) = NA_REAL;
          }
          
        } else {
          arma_scale(i, j) = NA_REAL;
        }
        
      }
    }
  }
  
};

// [[Rcpp::export]]
NumericMatrix roll_sd(const NumericMatrix& data, const int& width,
                      const arma::vec& weights, const bool& center,
                      const int& min_obs, const bool& complete_obs,
                      const bool& na_restore, const std::string& parallel_for) {
  
  int n_rows = data.nrow();
  int n_cols = data.ncol();
  arma::uvec arma_any_na(n_rows);
  arma::mat arma_center(n_rows, n_cols);
  arma::mat arma_scale(n_rows, n_cols);
  
  // check 'width' argument for errors
  check_width(width, n_rows);
  
  // default 'weights' argument is equal-weighted,
  // otherwise check argument for errors
  check_weights(weights, width);
  
  // default 'min_obs' argument is 'width' (equivalent to 'na.rm = FALSE'),
  // otherwise check argument for errors
  check_min_obs(min_obs, width);
  
  // default 'complete_obs' argument is 'false' (equivalent to 'pairwise'),
  // otherwise check argument for errors
  if (complete_obs) {
    arma_any_na = any_na(data);
  } else {
    arma_any_na.fill(0);
  }
  
  // default 'center' argument subtracts mean of each variable,
  // otherwise no scaling is done
  if (center) {
    if (parallel_for == "rows") {
      RollMeanRows roll_mean_rows(data, n_rows, n_cols, width, weights,
                                  min_obs, arma_any_na, na_restore,
                                  arma_center);
      parallelFor(0, n_rows, roll_mean_rows);
    } else if (parallel_for == "cols") {
      RollMeanCols roll_mean_cols(data, n_rows, n_cols, width, weights,
                                  min_obs, arma_any_na, na_restore,
                                  arma_center);
      parallelFor(0, n_cols, roll_mean_cols);
    }
  }
  
  // compute rolling standard deviations
  if (parallel_for == "rows") {
    RollSdRows roll_sd_rows(data, n_rows, n_cols, width, weights,
                            center, arma_center,
                            min_obs, arma_any_na, na_restore,
                            arma_scale);
    parallelFor(0, n_rows, roll_sd_rows); 
  } else if (parallel_for == "cols") {
    RollSdCols roll_sd_cols(data, n_rows, n_cols, width, weights,
                            center, arma_center,
                            min_obs, arma_any_na, na_restore,
                            arma_scale);
    parallelFor(0, n_cols, roll_sd_cols);   
  }
  
  // create and return a matrix or xts object
  NumericMatrix result(wrap(arma_scale));
  List dimnames = data.attr("dimnames");
  result.attr("dimnames") = dimnames;
  result.attr("index") = data.attr("index");
  result.attr(".indexCLASS") = data.attr(".indexCLASS");
  result.attr(".indexTZ") = data.attr(".indexTZ");
  result.attr("tclass") = data.attr("tclass");
  result.attr("tzone") = data.attr("tzone");
  result.attr("class") = data.attr("class");
  
  return result;
  
}

// 'Worker' function for computing rolling scaling and centering
struct RollScaleCenterRows : public Worker {
  
  const RMatrix<double> data;   // source
  const int n_rows;
  const int n_cols;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const arma::mat arma_center;
  const bool scale;
  arma::mat& arma_scale;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_cov;          // destination (pass by reference)
  
  // initialize with source and destination
  RollScaleCenterRows(const NumericMatrix data, const int n_rows,
                      const int n_cols, const int width,
                      const arma::vec arma_weights, const bool center,
                      const arma::mat arma_center, const bool scale,
                      arma::mat& arma_scale, const int min_obs,
                      const arma::uvec arma_any_na, const bool na_restore,
                      arma::mat& arma_cov)
    : data(data), n_rows(n_rows),
      n_cols(n_cols), width(width),
      arma_weights(arma_weights), center(center),
      arma_center(arma_center), scale(scale),
      arma_scale(arma_scale), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      arma_cov(arma_cov) { }
  
  // function call operator that iterates by row
  void operator()(std::size_t begin_row, std::size_t end_row) {
    for (std::size_t i = begin_row; i < end_row; i++) {
      for (int j = 0; j < n_cols; j++) {
        
        int count = 0;
        int n_obs = 0;
        bool any_data = false;
        double sum_data = 0;
        
        // don't compute if missing value and 'na_restore' argument is true
        if ((!na_restore) || (na_restore && !std::isnan(data(i, j)))) {
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= (unsigned)count)) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is false
            if ((arma_any_na[i - count] == 0) && !std::isnan(data(i - count, j))) {
              
              // keep first non-missing value
              if (!any_data) {
                sum_data = data(i - count, j);
              }
              
              any_data = true;
              n_obs += 1;
              
            }
            
            count += 1;
            
          }
          
          // compute the scaling and centering
          if (n_obs >= min_obs) {
            
            // compute with 'center' and 'scale' arguments
            if (center && scale) {
              arma_cov(i, j) = (sum_data - arma_center(i, j)) / sqrt(arma_scale(i, j));
            } else if (!center && scale) {
              arma_cov(i, j) = (sum_data) / sqrt(arma_scale(i, j));
            } else if (center && !scale) {
              arma_cov(i, j) = (sum_data - arma_center(i, j));
            } else if (!center && !scale) {
              arma_cov(i, j) = (sum_data);
            }
            
          } else {
            arma_cov(i, j) = NA_REAL;
          }
          
        } else {
          arma_cov(i, j) = NA_REAL;
        }
        
      }
    }
  }
  
};

// 'Worker' function for computing rolling scaling and centering
struct RollScaleCenterCols : public Worker {
  
  const RMatrix<double> data;   // source
  const int n_rows;
  const int n_cols;
  const int width;
  const arma::vec arma_weights;
  const bool center;
  const arma::mat arma_center;
  const bool scale;
  arma::mat& arma_scale;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_cov;          // destination (pass by reference)
  
  // initialize with source and destination
  RollScaleCenterCols(const NumericMatrix data, const int n_rows,
                      const int n_cols, const int width,
                      const arma::vec arma_weights, const bool center,
                      const arma::mat arma_center, const bool scale,
                      arma::mat& arma_scale, const int min_obs,
                      const arma::uvec arma_any_na, const bool na_restore,
                      arma::mat& arma_cov)
    : data(data), n_rows(n_rows),
      n_cols(n_cols), width(width),
      arma_weights(arma_weights), center(center),
      arma_center(arma_center), scale(scale),
      arma_scale(arma_scale), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      arma_cov(arma_cov) { }
  
  // function call operator that iterates by row
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      for (int i = 0; i < n_rows; i++) {
        
        int count = 0;
        int n_obs = 0;
        bool any_data = false;
        double sum_data = 0;
        
        // don't compute if missing value and 'na_restore' argument is true
        if ((!na_restore) || (na_restore && !std::isnan(data(i, j)))) {
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is false
            if ((arma_any_na[i - count] == 0) && !std::isnan(data(i - count, j))) {
              
              // keep first non-missing value
              if (!any_data) {
                sum_data = data(i - count, j);
              }
              
              any_data = true;
              n_obs += 1;
              
            }
            
            count += 1;
            
          }
          
          // compute the scaling and centering
          if (n_obs >= min_obs) {
            
            // compute with 'center' and 'scale' arguments
            if (center && scale) {
              arma_cov(i, j) = (sum_data - arma_center(i, j)) / sqrt(arma_scale(i, j));
            } else if (!center && scale) {
              arma_cov(i, j) = (sum_data) / sqrt(arma_scale(i, j));
            } else if (center && !scale) {
              arma_cov(i, j) = (sum_data - arma_center(i, j));
            } else if (!center && !scale) {
              arma_cov(i, j) = (sum_data);
            }
            
          } else {
            arma_cov(i, j) = NA_REAL;
          }
          
        } else {
          arma_cov(i, j) = NA_REAL;
        }
        
      }
    }
  }
  
};

// [[Rcpp::export]]
NumericMatrix roll_scale(const NumericMatrix& data, const int& width,
                         const arma::vec& weights, const bool& center,
                         const bool& scale, const int& min_obs,
                         const bool& complete_obs, const bool& na_restore,
                         const std::string& parallel_for) {
  
  int n_rows = data.nrow();
  int n_cols = data.ncol();
  arma::uvec arma_any_na(n_rows);
  arma::mat arma_center(n_rows, n_cols);
  arma::mat arma_scale(n_rows, n_cols);
  arma::mat arma_cov(n_rows, n_cols);
  
  // check 'width' argument for errors
  check_width(width, n_rows);
  
  // default 'weights' argument is equal-weighted,
  // otherwise check argument for errors
  check_weights(weights, width);
  
  // default 'min_obs' argument is 'width' (equivalent to 'na.rm = FALSE'),
  // otherwise check argument for errors
  check_min_obs(min_obs, width);
  
  // default 'complete_obs' argument is 'false' (equivalent to 'pairwise'),
  // otherwise check argument for errors
  if (complete_obs) {
    arma_any_na = any_na(data);
  } else {
    arma_any_na.fill(0);
  }
  
  // default 'center' argument subtracts mean of each variable,
  // otherwise zero is used
  if (center) {
    if (parallel_for == "rows") {
      RollMeanRows roll_mean_rows(data, n_rows, n_cols, width, weights,
                                  min_obs, arma_any_na, na_restore,
                                  arma_center);
      parallelFor(0, n_rows, roll_mean_rows);
    } else if (parallel_for == "cols") {
      RollMeanCols roll_mean_cols(data, n_rows, n_cols, width, weights,
                                  min_obs, arma_any_na, na_restore,
                                  arma_center);
      parallelFor(0, n_cols, roll_mean_cols);
    }
  }
  
  // default 'scale' argument divides by the standard deviation of each variable
  // otherwise no scaling is done
  if (scale) {
    if (parallel_for == "rows") {
      RollVarRows roll_var_rows(data, n_rows, n_cols, width, weights,
                                center, arma_center,
                                min_obs, arma_any_na, na_restore,
                                arma_scale);
      parallelFor(0, n_rows, roll_var_rows); 
    } else if (parallel_for == "cols") {
      RollVarCols roll_var_cols(data, n_rows, n_cols, width, weights,
                                center, arma_center,
                                min_obs, arma_any_na, na_restore,
                                arma_scale);
      parallelFor(0, n_cols, roll_var_cols);   
    }
  }
  
  // compute rolling scaling and centering
  if (parallel_for == "rows") {
    RollScaleCenterRows roll_scale_center_rows(data, n_rows, n_cols, width, weights,
                                               center, arma_center, scale, arma_scale,
                                               min_obs, arma_any_na, na_restore,
                                               arma_cov);
    parallelFor(0, n_rows, roll_scale_center_rows); 
  } else if (parallel_for == "cols") {
    RollScaleCenterCols roll_scale_center_cols(data, n_rows, n_cols, width, weights,
                                               center, arma_center, scale, arma_scale,
                                               min_obs, arma_any_na, na_restore,
                                               arma_cov);
    parallelFor(0, n_cols, roll_scale_center_cols); 
  }
  
  // create and return a matrix or xts object
  NumericMatrix result(wrap(arma_cov));
  List dimnames = data.attr("dimnames");
  result.attr("dimnames") = dimnames;
  result.attr("index") = data.attr("index");
  result.attr(".indexCLASS") = data.attr(".indexCLASS");
  result.attr(".indexTZ") = data.attr(".indexTZ");
  result.attr("tclass") = data.attr("tclass");
  result.attr("tzone") = data.attr("tzone");
  result.attr("class") = data.attr("class");
  
  return result;
  
}

// 'Worker' function for computing rolling means
struct RollMeanRowsCube : public Worker {
  
  const RMatrix<double> data;     // source
  const int n_rows;
  const int n_cols;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::cube& arma_center_j;     // destination (pass by reference)
  arma::cube& arma_center_k;     // destination (pass by reference)
  
  // initialize with source and destination
  RollMeanRowsCube(const NumericMatrix data, const int n_rows,
                   const int n_cols, const int width,
                   const arma::vec arma_weights, const int min_obs,
                   const arma::uvec arma_any_na, const bool na_restore,
                   arma::cube& arma_center_j, arma::cube& arma_center_k)
    : data(data), n_rows(n_rows),
      n_cols(n_cols), width(width),
      arma_weights(arma_weights), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      arma_center_j(arma_center_j), arma_center_k(arma_center_k) { }
  
  // function call operator that iterates by row
  void operator()(std::size_t begin_row, std::size_t end_row) {
    for (std::size_t i = begin_row; i < end_row; i++) {
      for (int j = 0; j < n_cols; j++) {
        for (int k = 0; k <= j; k++) {
          
          int count = 0;
          int n_obs = 0;
          double sum_j = 0;
          double sum_k = 0;
          double sum_weights = 0;
          
          // don't compute if missing value and 'na_restore' argument is true
          if ((!na_restore) ||
              (na_restore && !std::isnan(data(i, j)) && !std::isnan(data(i, k)))) {
            
            // number of observations is either the window size or,
            // for partial results, the number of the current row
            while ((width > count) && (i >= (unsigned)count)) {
              
              // don't include if missing value and 'any_na' argument is 1
              // note: 'any_na' is set to 0 if 'complete_obs' argument is false
              if ((arma_any_na[i - count] == 0) &&
                  !std::isnan(data(i - count, j)) &&
                  !std::isnan(data(i - count, k))) {
                  
                  // compute the rolling sum
                  sum_j += data(i - count, j) * arma_weights[width - count - 1];
                sum_k += data(i - count, k) * arma_weights[width - count - 1];
                sum_weights += arma_weights[width - count - 1];
                n_obs += 1;
                
              }
              
              count += 1;
              
            }
            
            // compute the mean
            if (n_obs >= min_obs) {
              arma_center_j(k, j, i) = sum_j / sum_weights;
              arma_center_k(k, j, i) = sum_k / sum_weights;
            } else {
              arma_center_j(k, j, i) = NA_REAL;
              arma_center_k(k, j, i) = NA_REAL;
            }
            
          } else {
            arma_center_j(k, j, i) = NA_REAL;
            arma_center_k(k, j, i) = NA_REAL;
          }
          
          // covariance matrix is symmetric
          arma_center_j(j, k, i) = arma_center_j(k, j, i);
          arma_center_k(j, k, i) = arma_center_k(k, j, i);
          
        }
        
      }
    }
  }
  
};

// 'Worker' function for computing rolling means
struct RollMeanColsCube : public Worker {
  
  const RMatrix<double> data;     // source
  const int n_rows;
  const int n_cols;
  const int width;
  const arma::vec arma_weights;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::cube& arma_center_j;     // destination (pass by reference)
  arma::cube& arma_center_k;     // destination (pass by reference)
  
  // initialize with source and destination
  RollMeanColsCube(const NumericMatrix data, const int n_rows,
                   const int n_cols, const int width,
                   const arma::vec arma_weights, const int min_obs,
                   const arma::uvec arma_any_na, const bool na_restore,
                   arma::cube& arma_center_j, arma::cube& arma_center_k)
    : data(data), n_rows(n_rows),
      n_cols(n_cols), width(width),
      arma_weights(arma_weights), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      arma_center_j(arma_center_j), arma_center_k(arma_center_k) { }
  
  // function call operator that iterates by column
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      for (std::size_t k = 0; k <= j; k++) {  
        for (int i = 0; i < n_rows; i++) {
          
          int count = 0;
          int n_obs = 0;
          double sum_j = 0;
          double sum_k = 0;
          double sum_weights = 0;
          
          // don't compute if missing value and 'na_restore' argument is true
          if ((!na_restore) ||
              (na_restore && !std::isnan(data(i, j)) && !std::isnan(data(i, k)))) {
            
            // number of observations is either the window size or,
            // for partial results, the number of the current row
            while ((width > count) && (i >= count)) {
              
              // don't include if missing value and 'any_na' argument is 1
              // note: 'any_na' is set to 0 if 'complete_obs' argument is false
              if ((arma_any_na[i - count] == 0) &&
                  !std::isnan(data(i - count, j)) &&
                  !std::isnan(data(i - count, k))) {
                  
                  // compute the rolling sum
                  sum_j += data(i - count, j) * arma_weights[width - count - 1];
                sum_k += data(i - count, k) * arma_weights[width - count - 1];
                sum_weights += arma_weights[width - count - 1];
                n_obs += 1;
                
              }
              
              count += 1;
              
            }
            
            // compute the mean
            if (n_obs >= min_obs) {
              arma_center_j(k, j, i) = sum_j / sum_weights;
              arma_center_k(k, j, i) = sum_k / sum_weights;
            } else {
              arma_center_j(k, j, i) = NA_REAL;
              arma_center_k(k, j, i) = NA_REAL;
            }
            
          } else {
            arma_center_j(k, j, i) = NA_REAL;
            arma_center_k(k, j, i) = NA_REAL;
          }
          
          // covariance matrix is symmetric
          arma_center_j(j, k, i) = arma_center_j(k, j, i);
          arma_center_k(j, k, i) = arma_center_k(k, j, i);
          
        }
        
      }
    }
  }
  
};

// 'Worker' function for computing rolling variances
struct RollVarRowsCube : public Worker {
  
  const RMatrix<double> data;     // source
  const int n_rows;
  const int n_cols;
  const int width;
  const arma::vec arma_weights;
  const bool center_x;
  const bool center_y;
  const arma::cube arma_center_j;
  const arma::cube arma_center_k;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::cube& arma_scale_j;       // destination (pass by reference)
  arma::cube& arma_scale_k;       // destination (pass by reference)
  
  // initialize with source and destination
  RollVarRowsCube(const NumericMatrix data, const int n_rows,
                  const int n_cols, const int width,
                  const arma::vec arma_weights, const bool center_x,
                  const bool center_y, const arma::cube arma_center_j,
                  const arma::cube arma_center_k, const int min_obs,
                  const arma::uvec arma_any_na, const bool na_restore,
                  arma::cube& arma_scale_j, arma::cube& arma_scale_k)
    : data(data), n_rows(n_rows),
      n_cols(n_cols), width(width),
      arma_weights(arma_weights), center_x(center_x),
      center_y(center_y),  arma_center_j(arma_center_j),
      arma_center_k(arma_center_k), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      arma_scale_j(arma_scale_j), arma_scale_k(arma_scale_k) { }
  
  // function call operator that iterates by row
  void operator()(std::size_t begin_row, std::size_t end_row) {
    for (std::size_t i = begin_row; i < end_row; i++) {
      for (int j = 0; j < n_cols; j++) {
        for (int k = 0; k <= j; k++) {
          
          int count = 0;
          int n_obs = 0;
          double sum_j = 0;
          double sum_k = 0;
          double sum_weights = 0;
          double sum_weights_sq = 0;
          
          // don't compute if missing value and 'na_restore' argument is true
          if ((!na_restore) ||
              (na_restore && !std::isnan(data(i, j)) && !std::isnan(data(i, k)))) {
            
            // number of observations is either the window size or,
            // for partial results, the number of the current row
            while ((width > count) && (i >= (unsigned)count)) {
              
              // don't include if missing value and 'any_na' argument is 1
              // note: 'any_na' is set to 0 if 'complete_obs' argument is false
              if ((arma_any_na[i - count] == 0) &&
                  !std::isnan(data(i - count, j))
                  && !std::isnan(data(i - count, k))) {
                  
                  // compute the rolling sum of squares with 'center' argument
                  if (center_x) {
                    if (j != n_cols - 1) {
                      sum_j += pow(data(i - count, j) - arma_center_j(j, k, i), 2.0) *
                        arma_weights[width - count - 1];
                    }
                    if (k != n_cols - 1) {
                      sum_k += pow(data(i - count, k) - arma_center_k(j, k, i), 2.0) *
                        arma_weights[width - count - 1];
                    }
                  } else if (!center_x) {
                    if (j != n_cols - 1) {
                      sum_j += pow(data(i - count, j), 2.0) * arma_weights[width - count - 1];
                    }
                    if (k != n_cols - 1) {
                      sum_k += pow(data(i - count, k), 2.0) * arma_weights[width - count - 1];
                    }
                  }
                  
                  if (center_y) {
                    if (j == n_cols - 1) {
                      sum_j += pow(data(i - count, j) - arma_center_j(j, k, i), 2.0) *
                        arma_weights[width - count - 1];
                    }
                    if (k == n_cols - 1){
                      sum_k += pow(data(i - count, k) - arma_center_k(j, k, i), 2.0) *
                        arma_weights[width - count - 1];
                    }
                  } else if (!center_y) {
                    if (j == n_cols - 1) {
                      sum_j += pow(data(i - count, j), 2.0) * arma_weights[width - count - 1];
                    }
                    if (k == n_cols - 1) {
                      sum_k += pow(data(i - count, k), 2.0) * arma_weights[width - count - 1];
                    }
                  }
                  
                  sum_weights += arma_weights[width - count - 1];
                  sum_weights_sq += pow(arma_weights[width - count - 1], 2.0);
                  n_obs += 1;
                  
              }
              
              count += 1;
              
            }
            
            // compute the mean
            if (n_obs >= min_obs) {
              arma_scale_j(k, j, i) = ((sum_j / sum_weights) /
                                         (1 - (sum_weights_sq / pow(sum_weights, 2.0))));
              arma_scale_k(k, j, i) = ((sum_k / sum_weights) /
                                         (1 - (sum_weights_sq / pow(sum_weights, 2.0))));
            } else {
              arma_scale_j(k, j, i) = NA_REAL;
              arma_scale_k(k, j, i) = NA_REAL;
            }
            
          } else {
            arma_scale_j(k, j, i) = NA_REAL;
            arma_scale_k(k, j, i) = NA_REAL;
          }
          
          // covariance matrix is symmetric
          arma_scale_j(j, k, i) = arma_scale_j(k, j, i);
          arma_scale_k(j, k, i) = arma_scale_k(k, j, i);
          
        }
        
      }
    }
  }
  
};

// 'Worker' function for computing rolling variances
struct RollVarColsCube : public Worker {
  
  const RMatrix<double> data;     // source
  const int n_rows;
  const int n_cols;
  const int width;
  const arma::vec arma_weights;
  const bool center_x;
  const bool center_y;
  const arma::cube arma_center_j;
  const arma::cube arma_center_k;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::cube& arma_scale_j;     // destination (pass by reference)
  arma::cube& arma_scale_k;     // destination (pass by reference)
  
  // initialize with source and destination
  RollVarColsCube(const NumericMatrix data, const int n_rows,
                  const int n_cols, const int width,
                  const arma::vec arma_weights, const bool center_x, 
                  const bool center_y, const arma::cube arma_center_j,
                  const arma::cube arma_center_k, const int min_obs,
                  const arma::uvec arma_any_na, const bool na_restore,
                  arma::cube& arma_scale_j, arma::cube& arma_scale_k)
    : data(data), n_rows(n_rows),
      n_cols(n_cols), width(width),
      arma_weights(arma_weights), center_x(center_x), 
      center_y(center_y), arma_center_j(arma_center_j),
      arma_center_k(arma_center_k), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      arma_scale_j(arma_scale_j), arma_scale_k(arma_scale_k) { }
  
  // function call operator that iterates by column
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      for (std::size_t k = 0; k <= j; k++) {  
        for (int i = 0; i < n_rows; i++) {
          
          int count = 0;
          int n_obs = 0;
          double sum_j = 0;
          double sum_k = 0;
          double sum_weights = 0;
          double sum_weights_sq = 0;
          
          // don't compute if missing value and 'na_restore' argument is true
          if ((!na_restore) ||
              (na_restore && !std::isnan(data(i, j)) && !std::isnan(data(i, k)))) {
            
            // number of observations is either the window size or,
            // for partial results, the number of the current row
            while ((width > count) && (i >= count)) {
              
              // don't include if missing value and 'any_na' argument is 1
              // note: 'any_na' is set to 0 if 'complete_obs' argument is false
              if ((arma_any_na[i - count] == 0) &&
                  !std::isnan(data(i - count, j)) &&
                  !std::isnan(data(i - count, k))) {
                  
                  // compute the rolling sum of squares with 'center' argument
                  if (center_x) {
                    if ((int)j != n_cols - 1) {
                      sum_j += pow(data(i - count, j) - arma_center_j(j, k, i), 2.0) *
                        arma_weights[width - count - 1];
                    }
                    if ((int)k != n_cols - 1) {
                      sum_k += pow(data(i - count, k) - arma_center_k(j, k, i), 2.0) *
                        arma_weights[width - count - 1];
                    }
                  } else if (!center_x) {
                    if ((int)j != n_cols - 1) {
                      sum_j += pow(data(i - count, j), 2.0) * arma_weights[width - count - 1];
                    }
                    if ((int)k != n_cols - 1) {
                      sum_k += pow(data(i - count, k), 2.0) * arma_weights[width - count - 1];
                    }
                  }
                  
                  if (center_y) {
                    if ((int)j == n_cols - 1) {
                      sum_j += pow(data(i - count, j) - arma_center_j(j, k, i), 2.0) *
                        arma_weights[width - count - 1];
                    }
                    if ((int)k == n_cols - 1){
                      sum_k += pow(data(i - count, k) - arma_center_k(j, k, i), 2.0) *
                        arma_weights[width - count - 1];
                    }
                  } else if (!center_y) {
                    if ((int)j == n_cols - 1) {
                      sum_j += pow(data(i - count, j), 2.0) * arma_weights[width - count - 1];
                    }
                    if ((int)k == n_cols - 1) {
                      sum_k += pow(data(i - count, k), 2.0) * arma_weights[width - count - 1];
                    }
                  }
                  
                  sum_weights += arma_weights[width - count - 1];
                  sum_weights_sq += pow(arma_weights[width - count - 1], 2.0);
                  n_obs += 1;
                  
              }
              
              count += 1;
              
            }
            
            // compute the mean
            if (n_obs >= min_obs) {
              arma_scale_j(k, j, i) = ((sum_j / sum_weights) /
                                         (1 - (sum_weights_sq / pow(sum_weights, 2.0))));
              arma_scale_k(k, j, i) = ((sum_k / sum_weights) /
                                         (1 - (sum_weights_sq / pow(sum_weights, 2.0))));
            } else {
              arma_scale_j(k, j, i) = NA_REAL;
              arma_scale_k(k, j, i) = NA_REAL;
            }
            
          } else {
            arma_scale_j(k, j, i) = NA_REAL;
            arma_scale_k(k, j, i) = NA_REAL;
          }
          
          // covariance matrix is symmetric
          arma_scale_j(j, k, i) = arma_scale_j(k, j, i);
          arma_scale_k(j, k, i) = arma_scale_k(k, j, i);
          
        }
        
      }
    }
  }
  
};

// 'Worker' function for computing rolling covariance matrices
struct RollCovRows : public Worker {
  
  const RMatrix<double> data;     // source
  const int n_rows;
  const int n_cols;
  const int width;
  const arma::vec arma_weights;
  const bool center_x;
  const bool center_y;
  const arma::cube arma_center_j;
  const arma::cube arma_center_k;
  const bool scale_x;
  const bool scale_y;
  const arma::cube arma_scale_j;
  const arma::cube arma_scale_k;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::cube& arma_cov;           // destination (pass by reference)
  
  // initialize with source and destination
  RollCovRows(const NumericMatrix data, const int n_rows,
              const int n_cols, const int width,
              const arma::vec arma_weights, const bool center_x, 
              const bool center_y, const arma::cube arma_center_j,
              const arma::cube arma_center_k, const bool scale_x,
              const bool scale_y, const arma::cube arma_scale_j,
              const arma::cube arma_scale_k, const int min_obs,
              const arma::uvec arma_any_na, const bool na_restore,
              arma::cube& arma_cov)
    : data(data), n_rows(n_rows),
      n_cols(n_cols), width(width),
      arma_weights(arma_weights), center_x(center_x),
      center_y(center_y), arma_center_j(arma_center_j),
      arma_center_k(arma_center_k), scale_x(scale_x),
      scale_y(scale_y), arma_scale_j(arma_scale_j),
      arma_scale_k(arma_scale_k), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      arma_cov(arma_cov) { }
  
  // function call operator that iterates by row
  void operator()(std::size_t begin_row, std::size_t end_row) {
    for (std::size_t i = begin_row; i < end_row; i++) {
      for (int j = 0; j < n_cols; j++) {
        for (int k = 0; k <= j; k++) {
          
          int count = 0;
          int n_obs = 0;
          double sum_x = 0;
          double sum_y = 0;
          double sum_data = 0;
          double sum_weights = 0;
          double sum_weights_sq = 0;
          
          // don't compute if missing value and 'na_restore' argument is true
          if ((!na_restore)
                || (na_restore && !std::isnan(data(i, j)) && !std::isnan(data(i, k)))) {
            
            // number of observations is either the window size or,
            // for partial results, the number of the current row
            while ((width > count) && (i >= (unsigned)count)) {
              
              // don't include if missing value and 'any_na' argument is 1
              // note: 'any_na' is set to 0 if 'complete_obs' argument is false
              if ((arma_any_na[i - count] == 0) &&
                  !std::isnan(data(i - count, j)) &&
                  !std::isnan(data(i - count, k))) {
                  
                  // compute the rolling sum of squares with 'center' and 'scale' arguments
                  if (center_x && center_y && scale_x && scale_y) {
                    sum_data +=
                      ((data(i - count, j) - arma_center_j(k, j, i)) / sqrt(arma_scale_j(k, j, i))) *
                      ((data(i - count, k) - arma_center_k(k, j, i)) / sqrt(arma_scale_k(k, j, i))) *
                      arma_weights[width - count - 1];
                  } else if (!center_x && !center_y && scale_x && scale_y) {
                    sum_data +=
                      ((data(i - count, j)) / sqrt(arma_scale_j(k, j, i))) *
                      ((data(i - count, k)) / sqrt(arma_scale_k(k, j, i))) *
                      arma_weights[width - count - 1];
                  } else if (center_x && center_y && !scale_x && !scale_y) {
                    sum_data +=
                      ((data(i - count, j) - arma_center_j(k, j, i))) *
                      ((data(i - count, k) - arma_center_k(k, j, i))) *
                      arma_weights[width - count - 1];
                  } else if (!center_x && !center_y && !scale_x && !scale_y) {
                    sum_data +=
                      ((data(i - count, j))) *
                      ((data(i - count, k))) *
                      arma_weights[width - count - 1];
                  } else {
                    
                    if (center_x && scale_x) {
                      if (j != n_cols - 1) {
                        sum_x = (data(i - count, j) - arma_center_j(k, j, i)) / sqrt(arma_scale_j(k, j, i));
                      }
                      if (k != n_cols - 1) {
                        sum_y = (data(i - count, k) - arma_center_k(k, j, i)) / sqrt(arma_scale_k(k, j, i));
                      }
                    } else if (!center_x && scale_x) {
                      if (j != n_cols - 1) {
                        sum_x = (data(i - count, j)) / sqrt(arma_scale_j(k, j, i));
                      }
                      if (k != n_cols - 1) {
                        sum_y = (data(i - count, k)) / sqrt(arma_scale_k(k, j, i));
                      }
                    } else if (center_x && !scale_x) {
                      if (j != n_cols - 1) {
                        sum_x = data(i - count, j) - arma_center_j(k, j, i);
                      }
                      if (k != n_cols - 1) {
                        sum_y = data(i - count, k) - arma_center_k(k, j, i);
                      }
                    } else if (!center_x && !scale_x) {
                      if (j != n_cols - 1) {
                        sum_x = data(i - count, j);
                      }
                      if (k != n_cols - 1) {
                        sum_y = data(i - count, k);
                      }
                    } 
                    
                    if (center_y && scale_y) {
                      if (j == n_cols - 1) {
                        sum_x = (data(i - count, j) - arma_center_j(k, j, i)) / sqrt(arma_scale_j(k, j, i));
                      }
                      if (k == n_cols - 1) {
                        sum_y = (data(i - count, k) - arma_center_k(k, j, i)) / sqrt(arma_scale_k(k, j, i));
                      }
                    } else if (!center_y && scale_y) {
                      if (j == n_cols - 1) {
                        sum_x = (data(i - count, j)) / sqrt(arma_scale_j(k, j, i));
                      }
                      if (k == n_cols - 1) {
                        sum_y = (data(i - count, k)) / sqrt(arma_scale_k(k, j, i));
                      }
                    } else if (center_y && !scale_y) {
                      if (j == n_cols - 1) {
                        sum_x = data(i - count, j) - arma_center_j(k, j, i);
                      }
                      if (k == n_cols - 1) {
                        sum_y = data(i - count, k) - arma_center_k(k, j, i);
                      }
                    } else if (!center_y && !scale_y) {
                      if (j == n_cols - 1) {
                        sum_x = data(i - count, j);
                      }
                      if (k == n_cols - 1) {
                        sum_y = data(i - count, k);
                      }
                    } 
                    
                    sum_data += sum_x * sum_y * arma_weights[width - count - 1];
                    
                  }
                  
                  sum_weights += arma_weights[width - count - 1];
                  sum_weights_sq += pow(arma_weights[width - count - 1], 2.0);
                  n_obs += 1;
                  
              }
              
              count += 1;
              
            }
            
            // compute the unbiased estimate of covariance
            if (n_obs >= min_obs) {
              arma_cov(k, j, i) = ((sum_data / sum_weights) /
                                     (1 - (sum_weights_sq / pow(sum_weights, 2.0))));
            } else {
              arma_cov(k, j, i) = NA_REAL;
            }
            
          } else {
            arma_cov(k, j, i) = NA_REAL;
          }
          
          // covariance matrix is symmetric
          arma_cov(j, k, i) = arma_cov(k, j, i);
          
        }
        
      }
    }
  }
  
};

// 'Worker' function for computing rolling covariance matrices
struct RollCovCols : public Worker {
  
  const RMatrix<double> data;     // source
  const int n_rows;
  const int n_cols;
  const int width;
  const arma::vec arma_weights;
  const bool center_x;
  const bool center_y;
  const arma::cube arma_center_j;
  const arma::cube arma_center_k;
  const bool scale_x;
  const bool scale_y;
  const arma::cube arma_scale_j;
  const arma::cube arma_scale_k;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::cube& arma_cov;           // destination (pass by reference)
  
  // initialize with source and destination
  RollCovCols(const NumericMatrix data, const int n_rows,
              const int n_cols, const int width,
              const arma::vec arma_weights, const bool center_x, 
              const bool center_y, const arma::cube arma_center_j,
              const arma::cube arma_center_k, const bool scale_x,
              const bool scale_y, const arma::cube arma_scale_j,
              const arma::cube arma_scale_k, const int min_obs,
              const arma::uvec arma_any_na, const bool na_restore,
              arma::cube& arma_cov)
    : data(data), n_rows(n_rows),
      n_cols(n_cols), width(width),
      arma_weights(arma_weights), center_x(center_x),
      center_y(center_y), arma_center_j(arma_center_j),
      arma_center_k(arma_center_k), scale_x(scale_x),
      scale_y(scale_y), arma_scale_j(arma_scale_j),
      arma_scale_k(arma_scale_k), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      arma_cov(arma_cov) { }
  
  // function call operator that iterates by column
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      for (std::size_t k = 0; k <= j; k++) {  
        for (int i = 0; i < n_rows; i++) {
          
          int count = 0;
          int n_obs = 0;
          double sum_x = 0;
          double sum_y = 0;
          double sum_data = 0;
          double sum_weights = 0;
          double sum_weights_sq = 0;
          
          // don't compute if missing value and 'na_restore' argument is true
          if ((!na_restore) ||
              (na_restore && !std::isnan(data(i, j)) && !std::isnan(data(i, k)))) {
            
            // number of observations is either the window size or,
            // for partial results, the number of the current row
            while ((width > count) && (i >= count)) {
              
              // don't include if missing value and 'any_na' argument is 1
              // note: 'any_na' is set to 0 if 'complete_obs' argument is false
              if ((arma_any_na[i - count] == 0) &&
                  !std::isnan(data(i - count, j)) &&
                  !std::isnan(data(i - count, k))) {
                  
                  // compute the rolling sum of squares with 'center' and 'scale' arguments
                  if (center_x && center_y && scale_x && scale_y) {
                    sum_data +=
                      ((data(i - count, j) - arma_center_j(k, j, i)) / sqrt(arma_scale_j(k, j, i))) *
                      ((data(i - count, k) - arma_center_k(k, j, i)) / sqrt(arma_scale_k(k, j, i))) *
                      arma_weights[width - count - 1];
                  } else if (!center_x && !center_y && scale_x && scale_y) {
                    sum_data +=
                      ((data(i - count, j)) / sqrt(arma_scale_j(k, j, i))) *
                      ((data(i - count, k)) / sqrt(arma_scale_k(k, j, i))) *
                      arma_weights[width - count - 1];
                  } else if (center_x && center_y && !scale_x && !scale_y) {
                    sum_data +=
                      ((data(i - count, j) - arma_center_j(k, j, i))) *
                      ((data(i - count, k) - arma_center_k(k, j, i))) *
                      arma_weights[width - count - 1];
                  } else if (!center_x && !center_y && !scale_x && !scale_y) {
                    sum_data +=
                      ((data(i - count, j))) *
                      ((data(i - count, k))) *
                      arma_weights[width - count - 1];
                  } else {
                    
                    if (center_x && scale_x) {
                      if ((int)j != n_cols - 1) {
                        sum_x = (data(i - count, j) - arma_center_j(k, j, i)) / sqrt(arma_scale_j(k, j, i));
                      }
                      if ((int)k != n_cols - 1) {
                        sum_y = (data(i - count, k) - arma_center_k(k, j, i)) / sqrt(arma_scale_k(k, j, i));
                      }
                    } else if (!center_x && scale_x) {
                      if ((int)j != n_cols - 1) {
                        sum_x = (data(i - count, j)) / sqrt(arma_scale_j(k, j, i));
                      }
                      if ((int)k != n_cols - 1) {
                        sum_y = (data(i - count, k)) / sqrt(arma_scale_k(k, j, i));
                      }
                    } else if (center_x && !scale_x) {
                      if ((int)j != n_cols - 1) {
                        sum_x = data(i - count, j) - arma_center_j(k, j, i);
                      }
                      if ((int)k != n_cols - 1) {
                        sum_y = data(i - count, k) - arma_center_k(k, j, i);
                      }
                    } else if (!center_x && !scale_x) {
                      if ((int)j != n_cols - 1) {
                        sum_x = data(i - count, j);
                      }
                      if ((int)k != n_cols - 1) {
                        sum_y = data(i - count, k);
                      }
                    }
                    
                    if (center_y && scale_y) {
                      if ((int)j == n_cols - 1) {
                        sum_x = (data(i - count, j) - arma_center_j(k, j, i)) / sqrt(arma_scale_j(k, j, i));
                      }
                      if ((int)k == n_cols - 1) {
                        sum_y = (data(i - count, k) - arma_center_k(k, j, i)) / sqrt(arma_scale_k(k, j, i));
                      }
                    } else if (!center_y && scale_y) {
                      if ((int)j == n_cols - 1) {
                        sum_x = (data(i - count, j)) / sqrt(arma_scale_j(k, j, i));
                      }
                      if ((int)k == n_cols - 1) {
                        sum_y = (data(i - count, k)) / sqrt(arma_scale_k(k, j, i));
                      }
                    } else if (center_y && !scale_y) {
                      if ((int)j == n_cols - 1) {
                        sum_x = data(i - count, j) - arma_center_j(k, j, i);
                      }
                      if ((int)k == n_cols - 1) {
                        sum_y = data(i - count, k) - arma_center_k(k, j, i);
                      }
                    } else if (!center_y && !scale_y) {
                      if ((int)j == n_cols - 1) {
                        sum_x = data(i - count, j);
                      }
                      if ((int)k == n_cols - 1) {
                        sum_y = data(i - count, k);
                      }
                    }
                    
                    sum_data += sum_x * sum_y * arma_weights[width - count - 1];
                    
                  }
                  
                  sum_weights += arma_weights[width - count - 1];
                  sum_weights_sq += pow(arma_weights[width - count - 1], 2.0);
                  n_obs += 1;
                  
              }
              
              count += 1;
              
            }
            
            // compute the unbiased estimate of covariance
            if (n_obs >= min_obs) {
              arma_cov(k, j, i) = ((sum_data / sum_weights) / 
                                     (1 - (sum_weights_sq / pow(sum_weights, 2.0))));
            } else {
              arma_cov(k, j, i) = NA_REAL;
            }
            
          } else {
            arma_cov(k, j, i) = NA_REAL;
          }
          
          // covariance matrix is symmetric
          arma_cov(j, k, i) = arma_cov(k, j, i);
          
        }
        
      }
    }
  }
  
};

// [[Rcpp::export]]
NumericVector roll_cov(const NumericMatrix& data, const int& width,
                       const arma::vec& weights, const bool& center,
                       const bool& scale, const int& min_obs,
                       const bool& complete_obs, const bool& na_restore,
                       const std::string& parallel_for) {
  
  int n_rows = data.nrow();
  int n_cols = data.ncol();
  bool center_x = center;
  bool center_y = center;
  bool scale_x = scale;
  bool scale_y = scale;
  arma::uvec arma_any_na(n_rows);
  arma::cube arma_center_j(n_cols, n_cols, n_rows);
  arma::cube arma_center_k(n_cols, n_cols, n_rows);
  arma::cube arma_scale_j(n_cols, n_cols, n_rows);
  arma::cube arma_scale_k(n_cols, n_cols, n_rows);
  arma::cube arma_scale_intercept(n_cols, n_cols, n_rows);
  arma::cube arma_cov(n_cols, n_cols, n_rows);
  
  // check 'width' argument for errors
  check_width(width, n_rows);
  
  // default 'weights' argument is equal-weighted,
  // otherwise check argument for errors
  check_weights(weights, width);
  
  // default 'min_obs' argument is 'width' (equivalent to 'na.rm = FALSE'),
  // otherwise check argument for errors
  check_min_obs(min_obs, width);
  
  // default 'complete_obs' argument is 'true' (equivalent to 'complete'),
  // otherwise check argument for errors
  if (complete_obs) {
    arma_any_na = any_na(data);
  } else {
    arma_any_na.fill(0);
  }
  
  // default 'center' argument subtracts mean of each variable,
  // otherwise no scaling is done
  if (center) {
    if (parallel_for == "rows") {
      RollMeanRowsCube roll_mean_rows(data, n_rows, n_cols, width, weights,
                                      min_obs, arma_any_na, na_restore,
                                      arma_center_j, arma_center_k);
      parallelFor(0, n_rows, roll_mean_rows);
    } else if (parallel_for == "cols") {
      RollMeanColsCube roll_mean_cols(data, n_rows, n_cols, width, weights,
                                      min_obs, arma_any_na, na_restore,
                                      arma_center_j, arma_center_k);
      parallelFor(0, n_cols, roll_mean_cols);
    }
  }
  
  // default 'scale' argument is none,
  // otherwise divide by the standard deviation of each variable
  if (scale) {
    if (parallel_for == "rows") {
      RollVarRowsCube roll_var_rows(data, n_rows, n_cols, width, weights,
                                    center_x, center_y,
                                    arma_center_j, arma_center_k,
                                    min_obs, arma_any_na, na_restore,
                                    arma_scale_j, arma_scale_k);
      parallelFor(0, n_rows, roll_var_rows);
    } else if (parallel_for == "cols") {
      RollVarColsCube roll_var_cols(data, n_rows, n_cols, width, weights,
                                    center_x, center_y,
                                    arma_center_j, arma_center_k,
                                    min_obs, arma_any_na, na_restore,
                                    arma_scale_j, arma_scale_k);
      parallelFor(0, n_cols, roll_var_cols);
    }
  }
  
  // compute rolling covariance matrices
  if (parallel_for == "rows") {
    RollCovRows roll_cov_rows(data, n_rows, n_cols, width, weights,
                              center_x, center_y,
                              arma_center_j, arma_center_k,
                              scale_x, scale_y,
                              arma_scale_j, arma_scale_k,
                              min_obs, arma_any_na, na_restore,
                              arma_cov);
    parallelFor(0, n_rows, roll_cov_rows); 
  } else if (parallel_for == "cols") {
    RollCovCols roll_cov_cols(data, n_rows, n_cols, width, weights,
                              center_x, center_y,
                              arma_center_j, arma_center_k,
                              scale_x, scale_y,
                              arma_scale_j, arma_scale_k,
                              min_obs, arma_any_na, na_restore,
                              arma_cov);
    parallelFor(0, n_cols, roll_cov_cols);   
  }
  
  // create and return a matrix
  NumericVector result(wrap(arma_cov));
  result.attr("dim") = IntegerVector::create(n_cols, n_cols, n_rows);
  List dimnames = data.attr("dimnames");
  if (dimnames.size() > 1) {
    result.attr("dimnames") = List::create(dimnames[1], dimnames[1]);
  }
  
  return result;
  
}

// 'Worker' function for computing rolling means (scaled)
struct RollMeanScaleRows : public Worker {
  
  const RMatrix<double> data;   // source
  const int n_rows;
  const int n_cols;
  const int width;
  const arma::vec arma_weights;
  const arma::cube arma_scale_j;
  const arma::cube arma_scale_k;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_center;       // destination (pass by reference)
  
  // initialize with source and destination
  RollMeanScaleRows(const NumericMatrix data, const int n_rows,
                    const int n_cols, const int width,
                    const arma::vec arma_weights, const arma::cube arma_scale_j,
                    const arma::cube arma_scale_k, const int min_obs,
                    const arma::uvec arma_any_na, const bool na_restore,
                    arma::mat& arma_center)
    : data(data), n_rows(n_rows),
      n_cols(n_cols), width(width),
      arma_weights(arma_weights), arma_scale_j(arma_scale_j),
      arma_scale_k(arma_scale_k), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      arma_center(arma_center) { }
  
  // function call operator that iterates by row 
  void operator()(std::size_t begin_row, std::size_t end_row) {
    for (std::size_t i = begin_row; i < end_row; i++) {
      for (int j = 0; j < n_cols; j++) {
        
        int count = 0;
        int n_obs = 0;
        double sum_data = 0;
        double sum_weights = 0;
        
        // don't compute if missing value and 'na_restore' argument is true
        if ((!na_restore) || (na_restore && !std::isnan(data(i, j)))) {
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= (unsigned)count)) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is false
            if ((arma_any_na[i - count] == 0) && !std::isnan(data(i - count, j))) {
              
              // compute the rolling sum             
              sum_data += (data(i - count, j) / sqrt(arma_scale_j(j, j, i))) *
                arma_weights[width - count - 1];
              sum_weights += arma_weights[width - count - 1];
              n_obs += 1;
              
            }
            
            count += 1;
            
          }
          
          // compute the mean
          if (n_obs >= min_obs) {
            arma_center(i, j) = sum_data / sum_weights;
          } else {
            arma_center(i, j) = NA_REAL;
          }
          
        } else {
          arma_center(i, j) = NA_REAL;
        }
        
      }
    }
  }
  
};

// 'Worker' function for computing rolling means (scaled)
struct RollMeanScaleCols : public Worker {
  
  const RMatrix<double> data;   // source
  const int n_rows;
  const int n_cols;
  const int width;
  const arma::vec arma_weights;
  const arma::cube arma_scale_j;
  const arma::cube arma_scale_k;
  const int min_obs;
  const arma::uvec arma_any_na;
  const bool na_restore;
  arma::mat& arma_center;       // destination (pass by reference)
  
  // initialize with source and destination
  RollMeanScaleCols(const NumericMatrix data, const int n_rows,
                    const int n_cols, const int width,
                    const arma::vec arma_weights, const arma::cube arma_scale_j,
                    const arma::cube arma_scale_k, const int min_obs,
                    const arma::uvec arma_any_na, const bool na_restore,
                    arma::mat& arma_center)
    : data(data), n_rows(n_rows),
      n_cols(n_cols), width(width),
      arma_weights(arma_weights), arma_scale_j(arma_scale_j),
      arma_scale_k(arma_scale_k), min_obs(min_obs),
      arma_any_na(arma_any_na), na_restore(na_restore),
      arma_center(arma_center) { }
  
  // function call operator that iterates by column 
  void operator()(std::size_t begin_col, std::size_t end_col) {
    for (std::size_t j = begin_col; j < end_col; j++) {
      for (int i = 0; i < n_rows; i++) {
        
        int count = 0;
        int n_obs = 0;
        double sum_data = 0;
        double sum_weights = 0;
        
        // don't compute if missing value and 'na_restore' argument is true
        if ((!na_restore) || (na_restore && !std::isnan(data(i, j)))) {
          
          // number of observations is either the window size or,
          // for partial results, the number of the current row
          while ((width > count) && (i >= count)) {
            
            // don't include if missing value and 'any_na' argument is 1
            // note: 'any_na' is set to 0 if 'complete_obs' argument is false
            if ((arma_any_na[i - count] == 0) && !std::isnan(data(i - count, j))) {
              
              // compute the rolling sum
              sum_data += (data(i - count, j) / sqrt(arma_scale_j(j, j, i))) *
                arma_weights[width - count - 1];
              sum_weights += arma_weights[width - count - 1];
              n_obs += 1;
              
            }
            
            count += 1;
            
          }
          
          // compute the mean
          if (n_obs >= min_obs) {
            arma_center(i, j) = sum_data / sum_weights;
          } else {
            arma_center(i, j) = NA_REAL;
          }
          
        } else {
          arma_center(i, j) = NA_REAL;
        }
        
      }
    }
  }
  
};

// 'Worker' function for rolling ordinary least squares regressions
struct RollLmSlices : public Worker {
  
  const arma::cube arma_cov;      // source
  const int n_rows;
  const int n_cols;
  const bool center_x;
  const bool center_y;
  const arma::cube arma_center_intercept;
  const bool scale_x;
  const bool scale_y;
  const arma::mat arma_scale_intercept;
  arma::mat& arma_coef;           // destination (pass by reference)
  arma::mat& arma_rsq;            // destination (pass by reference)
  
  // initialize with source and destination
  RollLmSlices(const arma::cube arma_cov, const int n_rows,
               const int n_cols, const bool center_x,
               const bool center_y, const arma::cube arma_center_intercept,
               const bool scale_x, const bool scale_y,
               const arma::mat arma_scale_intercept, arma::mat& arma_coef,
               arma::mat& arma_rsq)
    : arma_cov(arma_cov), n_rows(n_rows),
      n_cols(n_cols), center_x(center_x),
      center_y(center_y), arma_center_intercept(arma_center_intercept),
      scale_x(scale_x), scale_y(scale_y),
      arma_scale_intercept(arma_scale_intercept), arma_coef(arma_coef),
      arma_rsq(arma_rsq) { }
  
  // function call operator that iterates by slice
  void operator()(std::size_t begin_slice, std::size_t end_slice) {
    for (std::size_t i = begin_slice; i < end_slice; i++) {
      
      arma::mat x = arma_cov.slice(i);
      arma::mat A = x.submat(0, 0, n_cols - 2, n_cols - 2);
      arma::mat b = x.submat(0, n_cols - 1, n_cols - 2, n_cols - 1);
      arma::vec coef(n_cols - 1);
      
      int j = 0;
      int k = 0;
      bool any_na = false;
      
      // check if missing value is present
      while ((!any_na) && (j * k < (n_cols - 1) * (n_cols - 1))) {
        for (j = 0; j < n_cols; j++) {
          for (k = 0; k < n_cols; k++) {
            if (std::isnan(x(j, k)))
              any_na = true;
          }
        }
      }
      
      // don't compute if missing value 
      if (!any_na) {
        
        // check if solution is found      
        bool status = arma::solve(coef, A, b);
        
        // don't find approximate solution for rank deficient system
        if (status) {
          
          // get diagonal from 'center' cube
          arma::mat center_intercept = arma_center_intercept.slice(i);
          arma::vec center_intercept_diag = center_intercept.diag();
          
          // get row from 'scale' matrix
          arma::vec scale_intercept_row = trans(arma_scale_intercept.row(i));
          
          // intercept
          double intercept_x = 0;
          double intercept_y = 0;
          if (center_x && center_y) {
            if (scale_x) {
              intercept_x = as_scalar(trans(coef) * scale_intercept_row.subvec(0, n_cols - 2));
            } else if (!scale_x) {
              intercept_x = as_scalar(trans(coef) * center_intercept_diag.subvec(0, n_cols - 2));
            }
            if (scale_y) {
              intercept_y = scale_intercept_row[n_cols - 1];
            } else if (!scale_y) {
              intercept_y = center_intercept_diag[n_cols - 1];
            }
          } else if (!center_x && center_y) {
            if (scale_y) {
              intercept_y = scale_intercept_row[n_cols - 1];
            } else if (!scale_y) {
              intercept_y = center_intercept_diag[n_cols - 1];
            }
          } else if (center_x && !center_y) {
            if (scale_x) {
              intercept_x = as_scalar(trans(coef) * scale_intercept_row.subvec(0, n_cols - 2));
            } else if (!scale_x) {
              intercept_x = as_scalar(trans(coef) * center_intercept_diag.subvec(0, n_cols - 2));
            }
          }
          arma_coef(i, 0) = intercept_y - intercept_x;
          
          // coefficients
          arma_coef.submat(i, 1, i, n_cols - 1) = trans(coef);
          
          // r-squared
          arma_rsq(i, 0) = as_scalar((trans(coef) * A * coef) /
                                       x.submat(n_cols - 1, n_cols - 1, n_cols - 1, n_cols - 1));
          
        } else if (!status) {
          
          arma::vec no_solution(n_cols);
          no_solution.fill(NA_REAL);
          
          arma_coef.row(i) = trans(no_solution);
          arma_rsq(i, 0) = NA_REAL;
          
        }
        
      } else {
        
        arma::vec no_solution(n_cols);
        no_solution.fill(NA_REAL);
        
        arma_coef.row(i) = trans(no_solution);
        arma_rsq(i, 0) = NA_REAL;
        
      }
      
    }
  }
  
};

List roll_lm_z(const NumericMatrix& x, const NumericVector& y,
               const int& width, const arma::vec& weights,
               const bool& center_x, const bool& center_y,
               const bool& scale_x, const bool& scale_y,
               const int& min_obs, const bool& complete_obs,
               const bool& na_restore, const std::string& parallel_for) {
  
  int n_rows = x.nrow();  
  int n_cols = x.ncol() + 1;
  arma::uvec arma_any_na(n_rows);
  arma::cube arma_center_j(n_cols, n_cols, n_rows);
  arma::cube arma_center_k(n_cols, n_cols, n_rows);
  arma::cube arma_scale_j(n_cols, n_cols, n_rows);
  arma::cube arma_scale_k(n_cols, n_cols, n_rows);
  arma::mat arma_scale_intercept(n_rows, n_cols);
  arma::cube arma_cov(n_cols, n_cols, n_rows);
  arma::mat arma_coef(n_rows, n_cols);
  arma::mat arma_rsq(n_rows, 1);
  
  // check 'x' and 'y' arguments for errors
  check_ols(n_rows, y.size());
  
  // cbind x and y variables
  NumericMatrix data(n_rows, n_cols);
  std::copy(x.begin(), x.end(), data.begin());
  std::copy(y.begin(), y.end(), data.begin() + n_rows * (n_cols - 1));
  
  // check 'width' argument for errors
  check_width(width, n_rows);
  
  // default 'weights' argument is equal-weighted,
  // otherwise check argument for errors
  check_weights(weights, width);
  
  // default 'min_obs' argument is 'width' (equivalent to 'na.rm = FALSE'),
  // otherwise check argument for errors
  check_min_obs(min_obs, width);
  
  // default 'complete_obs' argument is 'true' (equivalent to 'complete'),
  // otherwise check argument for errors
  if (complete_obs) {
    arma_any_na = any_na(data);
  } else {
    arma_any_na.fill(0);
  }
  
  // default 'center' argument subtracts mean of each variable,
  // otherwise zero is used
  if (center_x || center_y) {
    if (parallel_for == "rows") {
      RollMeanRowsCube roll_mean_rows(data, n_rows, n_cols, width, weights,
                                      min_obs, arma_any_na, na_restore,
                                      arma_center_j, arma_center_k);
      parallelFor(0, n_rows, roll_mean_rows);
    } else if (parallel_for == "cols") {
      RollMeanColsCube roll_mean_cols(data, n_rows, n_cols, width, weights,
                                      min_obs, arma_any_na, na_restore,
                                      arma_center_j, arma_center_k);
      parallelFor(0, n_cols, roll_mean_cols);
    }
  }
  
  // default 'scale' argument is none,
  // otherwise divide by the standard deviation of each variable
  if (scale_x || scale_y) {
    if (parallel_for == "rows") {
      RollVarRowsCube roll_var_rows(data, n_rows, n_cols, width, weights,
                                    center_x, center_y, 
                                    arma_center_j, arma_center_k,
                                    min_obs, arma_any_na, na_restore,
                                    arma_scale_j, arma_scale_k);
      parallelFor(0, n_rows, roll_var_rows);
    } else if (parallel_for == "cols") {
      RollVarColsCube roll_var_cols(data, n_rows, n_cols, width, weights,
                                    center_x, center_y,
                                    arma_center_j, arma_center_k,
                                    min_obs, arma_any_na, na_restore,
                                    arma_scale_j, arma_scale_k);
      parallelFor(0, n_cols, roll_var_cols);
    }
  }
  
  // compute rolling covariance matrices
  if (parallel_for == "rows") {
    RollCovRows roll_cov_rows(data, n_rows, n_cols, width, weights,
                              center_x, center_y,
                              arma_center_j, arma_center_k,
                              scale_x, scale_y,
                              arma_scale_j, arma_scale_k,
                              min_obs, arma_any_na, na_restore,
                              arma_cov);
    parallelFor(0, n_rows, roll_cov_rows); 
  } else if (parallel_for == "cols") {
    RollCovCols roll_cov_cols(data, n_rows, n_cols, width, weights,
                              center_x, center_y,
                              arma_center_j, arma_center_k,
                              scale_x, scale_y,
                              arma_scale_j, arma_scale_k,
                              min_obs, arma_any_na, na_restore,
                              arma_cov);
    parallelFor(0, n_cols, roll_cov_cols);   
  }
  
  // compute rolling mean of each variable (scaled)
  if ((center_x || center_y) && (scale_x || scale_y)) {
    if (parallel_for == "rows") {
      RollMeanScaleRows roll_mean_scale_rows(data, n_rows, n_cols, width, weights,
                                             arma_scale_j, arma_scale_k, 
                                             min_obs, arma_any_na, na_restore,
                                             arma_scale_intercept);
      parallelFor(0, n_rows, roll_mean_scale_rows);
    } else if (parallel_for == "cols") {
      RollMeanScaleCols roll_mean_scale_cols(data, n_rows, n_cols, width, weights,
                                             arma_scale_j, arma_scale_k,
                                             min_obs, arma_any_na, na_restore,
                                             arma_scale_intercept);
      parallelFor(0, n_cols, roll_mean_scale_cols);
    }
  }
  
  // compute rolling ordinary least squares regressions
  RollLmSlices roll_lm_slices(arma_cov, n_rows, n_cols,
                              center_x, center_y, arma_center_j,
                              scale_x, scale_y, arma_scale_intercept,
                              arma_coef, arma_rsq);
  parallelFor(0, n_rows, roll_lm_slices);
  
  // create and return a list
  List result = List::create(Named("coefficients") = arma_coef,
                             Named("r.squared") = arma_rsq);
  
  return result;
  
}

// [[Rcpp::export]]
List roll_lm(const NumericMatrix& x, const NumericMatrix& y,
             const int& width, const arma::vec& weights,
             const bool& center_x, const bool& center_y,
             const bool& scale_x, const bool& scale_y,
             const int& min_obs, const bool& complete_obs,
             const bool& na_restore, const std::string& parallel_for) {
  
  int n_rows = x.nrow();  
  int n_cols = x.ncol() + 1;
  int y_n_cols = y.ncol();
  List result_z(2);
  List result_coef(y_n_cols);
  List result_rsq(y_n_cols);
  List result(2);
  
  // create a list of matrices,
  // otherwise a list of lists
  if (y_n_cols == 1) {
    
    result_z = roll_lm_z(x, y(_, 0),
                         width, weights,
                         center_x, center_y,
                         scale_x, scale_y,
                         min_obs, complete_obs,
                         na_restore, parallel_for);
    
    arma::mat arma_coef_z = result_z[0];
    arma::mat arma_rsq_z = result_z[1];
    
    // create and return a matrix or xts object for coefficients
    NumericVector coef(wrap(arma_coef_z));
    coef.attr("dim") = IntegerVector::create(n_rows, n_cols);
    List x_dimnames = x.attr("dimnames");
    coef.attr("dimnames") = dimnames_ols(x_dimnames, n_cols - 1);
    coef.attr("index") = x.attr("index");
    coef.attr(".indexCLASS") = x.attr(".indexCLASS");
    coef.attr(".indexTZ") = x.attr(".indexTZ");
    coef.attr("tclass") = x.attr("tclass");
    coef.attr("tzone") = x.attr("tzone");
    coef.attr("class") = x.attr("class");
    
    // create and return a matrix or xts object for r-squareds
    NumericVector rsq(wrap(arma_rsq_z));
    rsq.attr("dim") = IntegerVector::create(n_rows, 1);
    if (x_dimnames.size() > 1) {
      rsq.attr("dimnames") = List::create(x_dimnames[0], "R-squared");
    } else {
      rsq.attr("dimnames") = List::create(R_NilValue, "R-squared");
    }
    rsq.attr("index") = x.attr("index");
    rsq.attr(".indexCLASS") = x.attr(".indexCLASS");
    rsq.attr(".indexTZ") = x.attr(".indexTZ");
    rsq.attr("tclass") = x.attr("tclass");
    rsq.attr("tzone") = x.attr("tzone");
    rsq.attr("class") = x.attr("class");
    
    List y_dimnames = y.attr("dimnames");
    if (y_dimnames.size() > 1) {
      result_coef.attr("names") = y_dimnames[1];
      result_rsq.attr("names") = y_dimnames[1];
    }
    
    // create and return a list
    result = List::create(Named("coefficients") = coef,
                          Named("r.squared") = rsq);
    
  } else {
    
    for (int z = 0; z < y_n_cols; z++) {
      
      result_z = roll_lm_z(x, y(_, z),
                           width, weights,
                           center_x, center_y,
                           scale_x, scale_y,
                           min_obs, complete_obs,
                           na_restore, parallel_for);
      
      arma::mat arma_coef_z = result_z[0];
      arma::mat arma_rsq_z = result_z[1];
      
      // create and return a matrix or xts object for coefficients
      NumericVector coef(wrap(arma_coef_z));
      coef.attr("dim") = IntegerVector::create(n_rows, n_cols);
      List x_dimnames = x.attr("dimnames");
      coef.attr("dimnames") = dimnames_ols(x_dimnames, n_cols - 1);
      coef.attr("index") = x.attr("index");
      coef.attr(".indexCLASS") = x.attr(".indexCLASS");
      coef.attr(".indexTZ") = x.attr(".indexTZ");
      coef.attr("tclass") = x.attr("tclass");
      coef.attr("tzone") = x.attr("tzone");
      coef.attr("class") = x.attr("class");
      
      // create and return a matrix or xts object for r-squareds
      NumericVector rsq(wrap(arma_rsq_z));
      rsq.attr("dim") = IntegerVector::create(n_rows, 1);
      if (x_dimnames.size() > 1) {
        rsq.attr("dimnames") = List::create(x_dimnames[0], "R-squared");
      } else {
        rsq.attr("dimnames") = List::create(R_NilValue, "R-squared");
      }
      rsq.attr("index") = x.attr("index");
      rsq.attr(".indexCLASS") = x.attr(".indexCLASS");
      rsq.attr(".indexTZ") = x.attr(".indexTZ");
      rsq.attr("tclass") = x.attr("tclass");
      rsq.attr("tzone") = x.attr("tzone");
      rsq.attr("class") = x.attr("class");
      
      result_coef(z) = coef;
      result_rsq(z) = rsq;
      
    }
    
    // add names to each list
    List y_dimnames = y.attr("dimnames");
    if (y_dimnames.size() > 1) {
      result_coef.attr("names") = y_dimnames[1];
      result_rsq.attr("names") = y_dimnames[1];
    }
    
    // create and return a list
    result = List::create(Named("coefficients") = result_coef,
                          Named("r.squared") = result_rsq);
    
  }
  
  return result;
  
}

// 'Worker' function for rolling eigenvalues and eigenvectors
struct RollEigenSlices : public Worker {
  
  const arma::cube arma_cov;            // source
  const int n_rows;
  const int n_cols;
  arma::mat& arma_eigen_values;         // destination (pass by reference)
  arma::cube& arma_eigen_vectors;       // destination (pass by reference)
  
  // initialize with source and destination
  RollEigenSlices(const arma::cube arma_cov, const int n_rows,
                  const int n_cols, arma::mat& arma_eigen_values,
                  arma::cube& arma_eigen_vectors)
    : arma_cov(arma_cov), n_rows(n_rows),
      n_cols(n_cols), arma_eigen_values(arma_eigen_values),
      arma_eigen_vectors(arma_eigen_vectors) { }
  
  // function call operator that iterates by slice
  void operator()(std::size_t begin_slice, std::size_t end_slice) {
    for (std::size_t i = begin_slice; i < end_slice; i++) {
      
      arma::mat x = arma_cov.slice(i);
      arma::mat A = x.submat(0, 0, n_cols - 1, n_cols - 1);
      arma::vec eigen_values(n_cols);
      arma::mat eigen_vectors(n_cols, n_cols);
      
      int j = 0;
      int k = 0;
      bool any_na = false;
      
      // check if missing value is present
      while ((!any_na) && (j * k < (n_cols - 1) * (n_cols - 1))) {
        for (j = 0; j < n_cols; j++) {
          for (k = 0; k < n_cols; k++) {
            if (std::isnan(x(j, k)))
              any_na = true;
          }
        }
      }
      
      // don't compute if missing value 
      if (!any_na) {
        
        // check if solution is found
        bool status = arma::eig_sym(eigen_values, eigen_vectors, A);
        
        // don't find approximate solution for rank deficient system
        if (status) {
          
          // reverse order for consistency with R's eigen
          std::reverse(eigen_values.begin(), eigen_values.end());
          eigen_vectors = arma::fliplr(eigen_vectors);
          
          arma_eigen_values.row(i) = trans(eigen_values);
          arma_eigen_vectors.slice(i) = eigen_vectors;
          
        } else if (!status) {
          
          arma::vec no_solution_row(n_cols);
          no_solution_row.fill(NA_REAL);
          
          arma::mat no_solution_slice(n_cols, n_cols);
          no_solution_slice.fill(NA_REAL);
          
          arma_eigen_values.row(i) = trans(no_solution_row);
          arma_eigen_vectors.slice(i) = no_solution_slice;
          
        }
        
      } else {
        
        arma::vec no_solution_row(n_cols);
        no_solution_row.fill(NA_REAL);
        
        arma::mat no_solution_slice(n_cols, n_cols);
        no_solution_slice.fill(NA_REAL);
        
        arma_eigen_values.row(i) = trans(no_solution_row);
        arma_eigen_vectors.slice(i) = no_solution_slice;
        
      }
      
    }
  }
  
};

// [[Rcpp::export]]
List roll_eigen(const NumericMatrix& data, const int& width,
                const arma::vec& weights, const bool& center,
                const bool& scale, const int& min_obs,
                const bool& complete_obs, const bool& na_restore,
                const std::string& parallel_for) {
  
  int n_rows = data.nrow();  
  int n_cols = data.ncol();
  bool center_x = center;
  bool center_y = center;
  bool scale_x = scale;
  bool scale_y = scale;
  arma::uvec arma_any_na(n_rows);
  arma::cube arma_center_j(n_cols, n_cols, n_rows);
  arma::cube arma_center_k(n_cols, n_cols, n_rows);
  arma::cube arma_scale_j(n_cols, n_cols, n_rows);
  arma::cube arma_scale_k(n_cols, n_cols, n_rows);
  arma::cube arma_scale_intercept(n_cols, n_cols, n_rows);
  arma::cube arma_cov(n_cols, n_cols, n_rows);
  arma::mat arma_eigen_values(n_rows, n_cols);
  arma::cube arma_eigen_vectors(n_cols, n_cols, n_rows);
  
  // check 'width' argument for errors
  check_width(width, n_rows);
  
  // default 'weights' argument is equal-weighted,
  // otherwise check argument for errors
  check_weights(weights, width);
  
  // default 'min_obs' argument is 'width' (equivalent to 'na.rm = FALSE'),
  // otherwise check argument for errors
  check_min_obs(min_obs, width);
  
  // default 'complete_obs' argument is 'true' (equivalent to 'complete'),
  // otherwise check argument for errors
  if (complete_obs) {
    arma_any_na = any_na(data);
  } else {
    arma_any_na.fill(0);
  }
  
  // default 'center' argument subtracts mean of each variable,
  // otherwise zero is used
  if (center) {
    if (parallel_for == "rows") {
      RollMeanRowsCube roll_mean_rows(data, n_rows, n_cols, width, weights,
                                      min_obs, arma_any_na, na_restore,
                                      arma_center_j, arma_center_k);
      parallelFor(0, n_rows, roll_mean_rows);
    } else if (parallel_for == "cols") {
      RollMeanColsCube roll_mean_cols(data, n_rows, n_cols, width, weights,
                                      min_obs, arma_any_na, na_restore,
                                      arma_center_j, arma_center_k);
      parallelFor(0, n_cols, roll_mean_cols);
    }
  }
  
  // default 'scale' argument is none,
  // otherwise divide by the standard deviation of each variable
  if (scale) {
    if (parallel_for == "rows") {
      RollVarRowsCube roll_var_rows(data, n_rows, n_cols, width, weights,
                                    center_x, center_y, arma_center_j, arma_center_k,
                                    min_obs, arma_any_na, na_restore,
                                    arma_scale_j, arma_scale_k);
      parallelFor(0, n_rows, roll_var_rows);
    } else if (parallel_for == "cols") {
      RollVarColsCube roll_var_cols(data, n_rows, n_cols, width, weights,
                                    center_x, center_y, arma_center_j, arma_center_k,
                                    min_obs, arma_any_na, na_restore,
                                    arma_scale_j, arma_scale_k);
      parallelFor(0, n_cols, roll_var_cols);
    }
  }
  
  // compute rolling covariance matrices
  if (parallel_for == "rows") {
    RollCovRows roll_cov_rows(data, n_rows, n_cols, width, weights,
                              center_x, center_y,
                              arma_center_j, arma_center_k,
                              scale_x, scale_y,
                              arma_scale_j, arma_scale_k,
                              min_obs, arma_any_na, na_restore,
                              arma_cov);
    parallelFor(0, n_rows, roll_cov_rows); 
  } else if (parallel_for == "cols") {
    RollCovCols roll_cov_cols(data, n_rows, n_cols, width, weights,
                              center_x, center_y,
                              arma_center_j, arma_center_k,
                              scale_x, scale_y,
                              arma_scale_j, arma_scale_k,
                              min_obs, arma_any_na, na_restore,
                              arma_cov);
    parallelFor(0, n_cols, roll_cov_cols);   
  }
  
  // compute rolling eigenvalues and eigenvectors
  RollEigenSlices roll_eigen_slices(arma_cov, n_rows, n_cols,
                                    arma_eigen_values, arma_eigen_vectors);
  parallelFor(0, n_rows, roll_eigen_slices);
  
  // create and return a matrix or xts object for eigenvalues
  NumericMatrix eigen_values(wrap(arma_eigen_values));
  List dimnames = data.attr("dimnames");
  if (dimnames.size() > 1) {
    eigen_values.attr("dimnames") = List::create(dimnames[0], dimnames_pc(n_cols));
  } else {
    eigen_values.attr("dimnames") = List::create(R_NilValue, dimnames_pc(n_cols));
  }
  eigen_values.attr("index") = data.attr("index");
  eigen_values.attr(".indexCLASS") = data.attr(".indexCLASS");
  eigen_values.attr(".indexTZ") = data.attr(".indexTZ");
  eigen_values.attr("tclass") = data.attr("tclass");
  eigen_values.attr("tzone") = data.attr("tzone");
  eigen_values.attr("class") = data.attr("class");
  
  // create and return a cube for eigenvectors
  NumericVector eigen_vectors(wrap(arma_eigen_vectors));
  eigen_vectors.attr("dim") = IntegerVector::create(n_cols, n_cols, n_rows);
  if (dimnames.size() > 1) {
    eigen_vectors.attr("dimnames") = List::create(dimnames[1], dimnames_pc(n_cols));
  } else {
    eigen_vectors.attr("dimnames") = List::create(R_NilValue, dimnames_pc(n_cols));
  }
  
  // create and return a list
  List result = List::create(Named("values") = eigen_values,
                             Named("vectors") = eigen_vectors);
  
  return result;
  
}

// 'Worker' function for rolling OLS regressions
struct RollPcrSlices : public Worker {
  
  const arma::cube arma_cov;            // source
  const int n_rows;
  const int n_cols;
  const arma::uvec arma_cols;
  const arma::uvec arma_comps;
  const bool center_x;
  const bool center_y;
  const arma::cube arma_center_intercept;
  const bool scale_x;
  const bool scale_y;
  const arma::mat arma_scale_intercept;
  const arma::cube arma_eigen_vectors;
  arma::mat& arma_coef;                 // destination (pass by reference)
  arma::mat& arma_rsq;
  
  // initialize with source and destination
  RollPcrSlices(const arma::cube arma_cov, const int n_rows,
                const int n_cols, const arma::uvec arma_cols,
                const arma::uvec arma_comps, const bool center_x,
                const bool center_y, const arma::cube arma_center_intercept,
                const bool scale_x, const bool scale_y,
                const arma::mat arma_scale_intercept, const arma::cube arma_eigen_vectors,
                arma::mat& arma_coef, arma::mat& arma_rsq)
    : arma_cov(arma_cov), n_rows(n_rows),
      n_cols(n_cols), arma_cols(arma_cols),
      arma_comps(arma_comps), center_x(center_x),
      center_y(center_y), arma_center_intercept(arma_center_intercept),
      scale_x(scale_x), scale_y(scale_y),
      arma_scale_intercept(arma_scale_intercept), arma_eigen_vectors(arma_eigen_vectors),
      arma_coef(arma_coef), arma_rsq(arma_rsq) { }
  
  // function call operator that iterates by slice
  void operator()(std::size_t begin_slice, std::size_t end_slice) {
    for (std::size_t i = begin_slice; i < end_slice; i++) {
      
      arma::mat x = arma_cov.slice(i);
      arma::mat A = x.submat(0, 0, n_cols - 2, n_cols - 2);
      arma::mat b = x.submat(0, n_cols - 1, n_cols - 2, n_cols - 1);
      arma::vec gamma(n_cols - 1);
      arma::mat eigen_vectors = arma_eigen_vectors.slice(i);
      
      int j = 0;
      int k = 0;
      bool any_na = false;
      
      // check if missing value is present
      while ((!any_na) && (j * k < (n_cols - 1) * (n_cols - 1))) {
        for (j = 0; j < n_cols; j++) {
          for (k = 0; k < n_cols; k++) {
            if (std::isnan(x(j, k)))
              any_na = true;
          }
        }
      }
      
      // don't compute if missing value 
      if (!any_na) {
        
        // check if solution is found
        bool status = arma::solve(gamma, A * eigen_vectors, b);
        
        // don't find approximate solution for rank deficient system
        if (status) {
          
          // get diagonal from 'center' cube
          arma::mat center_intercept = arma_center_intercept.slice(i);
          arma::vec center_intercept_diag = center_intercept.diag();
          
          // get row from 'scale' matrix
          arma::vec scale_intercept_row = trans(arma_scale_intercept.row(i));
          
          arma::vec gamma_subset = gamma(arma_comps - 1);
          arma::vec coef = eigen_vectors.submat(arma_cols, arma_comps - 1) * gamma_subset;
          
          // intercept
          double intercept_x = 0;
          double intercept_y = 0;
          if (center_x && center_y) {
            if (scale_x) {
              intercept_x = as_scalar(trans(coef) * scale_intercept_row.subvec(0, n_cols - 2));
            } else if (!scale_x) {
              intercept_x = as_scalar(trans(coef) * center_intercept_diag.subvec(0, n_cols - 2));
            }
            if (scale_y) {
              intercept_y = scale_intercept_row[n_cols - 1];
            } else if (!scale_y) {
              intercept_y = center_intercept_diag[n_cols - 1];
            }
          } else if (!center_x && center_y) {
            if (scale_y) {
              intercept_y = scale_intercept_row[n_cols - 1];
            } else if (!scale_y) {
              intercept_y = center_intercept_diag[n_cols - 1];
            }
          } else if (center_x && !center_y) {
            if (scale_x) {
              intercept_x = as_scalar(trans(coef) * scale_intercept_row.subvec(0, n_cols - 2));
            } else if (!scale_x) {
              intercept_x = as_scalar(trans(coef) * center_intercept_diag.subvec(0, n_cols - 2));
            }
          }
          arma_coef(i, 0) = intercept_y - intercept_x;
          
          // coefficients
          arma_coef.submat(i, 1, i, n_cols - 1) = trans(coef);
          
          // r-squared
          arma_rsq(i, 0) = as_scalar((trans(coef) * A * coef) /
                                       x.submat(n_cols - 1, n_cols - 1, n_cols - 1, n_cols - 1));
          
        } else if (!status) {
          
          arma::vec no_solution(n_cols);
          no_solution.fill(NA_REAL);
          
          arma_coef.row(i) = trans(no_solution);
          arma_rsq[i] = NA_REAL;
          
        }
        
      } else {
        
        arma::vec no_solution(n_cols);
        no_solution.fill(NA_REAL);
        
        arma_coef.row(i) = trans(no_solution);
        arma_rsq[i] = NA_REAL;
        
      }
      
    }
  }
  
};

List roll_pcr_z(const NumericMatrix& x, const NumericVector& y,
                const int& width, const arma::uvec& comps,
                const arma::vec& weights, const bool& center_x,
                const bool& center_y, const bool& scale_x,
                const bool& scale_y, const int& min_obs,
                const bool& complete_obs, const bool& na_restore,
                const std::string& parallel_for) {
  
  int n_rows = x.nrow();
  int n_cols = x.ncol() + 1;
  arma::uvec arma_any_na(n_rows);
  arma::cube arma_center_j(n_cols, n_cols, n_rows);
  arma::cube arma_center_k(n_cols, n_cols, n_rows);
  arma::cube arma_scale_j(n_cols, n_cols, n_rows);
  arma::cube arma_scale_k(n_cols, n_cols, n_rows);
  arma::mat arma_scale_intercept(n_rows, n_cols);
  arma::cube arma_cov(n_cols, n_cols, n_rows);
  arma::mat arma_eigen_values(n_rows, n_cols - 1);
  arma::cube arma_eigen_vectors(n_cols - 1, n_cols - 1, n_rows);
  arma::mat arma_coef(n_rows, n_cols);
  arma::mat arma_rsq(n_rows, 1);
  
  // check 'x' and 'y' arguments for errors
  check_ols(n_rows, y.size());
  
  // cbind x and y variables
  NumericMatrix data(n_rows, n_cols);
  std::copy(x.begin(), x.end(), data.begin());
  std::copy(y.begin(), y.end(), data.begin() + n_rows * (n_cols - 1));
  
  // check 'width' argument for errors
  check_width(width, n_rows);
  
  // default 'comps' argument is all components
  check_comps(comps, n_cols - 1);
  
  // default 'weights' argument is equal-weighted,
  // otherwise check argument for errors
  check_weights(weights, width);
  
  // default 'min_obs' argument is 'width' (equivalent to 'na.rm = FALSE'),
  // otherwise check argument for errors
  check_min_obs(min_obs, width);
  
  // default 'complete_obs' argument is 'true' (equivalent to 'complete'),
  // otherwise check argument for errors
  if (complete_obs) {
    arma_any_na = any_na(data);
  } else {
    arma_any_na.fill(0);
  }
  
  // default 'center' argument subtracts mean of each variable,
  // otherwise zero is used
  if (center_x || center_y) {
    if (parallel_for == "rows") {
      RollMeanRowsCube roll_mean_rows(data, n_rows, n_cols, width, weights,
                                      min_obs, arma_any_na, na_restore,
                                      arma_center_j, arma_center_k);
      parallelFor(0, n_rows, roll_mean_rows);
    } else if (parallel_for == "cols") {
      RollMeanColsCube roll_mean_cols(data, n_rows, n_cols, width, weights,
                                      min_obs, arma_any_na, na_restore,
                                      arma_center_j, arma_center_k);
      parallelFor(0, n_cols, roll_mean_cols);
    }
  }
  
  // default 'scale' argument is none,
  // otherwise divide by the standard deviation of each variable
  if (scale_x || scale_y) {
    if (parallel_for == "rows") {
      RollVarRowsCube roll_var_rows(data, n_rows, n_cols, width, weights,
                                    center_x, center_y, arma_center_j, arma_center_k,
                                    min_obs, arma_any_na, na_restore,
                                    arma_scale_j, arma_scale_k);
      parallelFor(0, n_rows, roll_var_rows);
    } else if (parallel_for == "cols") {
      RollVarColsCube roll_var_cols(data, n_rows, n_cols, width, weights,
                                    center_x, center_y, arma_center_j, arma_center_k,
                                    min_obs, arma_any_na, na_restore,
                                    arma_scale_j, arma_scale_k);
      parallelFor(0, n_cols, roll_var_cols);
    }
  }
  
  // compute rolling mean of each variable (scaled)
  if ((center_x || center_y) && (scale_x || scale_y)) {
    if (parallel_for == "rows") {
      RollMeanScaleRows roll_mean_scale_rows(data, n_rows, n_cols, width, weights,
                                             arma_scale_j, arma_scale_k, 
                                             min_obs, arma_any_na, na_restore,
                                             arma_scale_intercept);
      parallelFor(0, n_rows, roll_mean_scale_rows);
    } else if (parallel_for == "cols") {
      RollMeanScaleCols roll_mean_scale_cols(data, n_rows, n_cols, width, weights,
                                             arma_scale_j, arma_scale_k,
                                             min_obs, arma_any_na, na_restore,
                                             arma_scale_intercept);
      parallelFor(0, n_cols, roll_mean_scale_cols);
    }
  }
  
  // compute rolling covariance matrices
  if (parallel_for == "rows") {
    RollCovRows roll_cov_rows(data, n_rows, n_cols, width, weights,
                              center_x, center_y,
                              arma_center_j, arma_center_k,
                              scale_x, scale_y,
                              arma_scale_j, arma_scale_k,
                              min_obs, arma_any_na, na_restore,
                              arma_cov);
    parallelFor(0, n_rows, roll_cov_rows); 
  } else if (parallel_for == "cols") {
    RollCovCols roll_cov_cols(data, n_rows, n_cols, width, weights,
                              center_x, center_y,
                              arma_center_j, arma_center_k,
                              scale_x, scale_y,
                              arma_scale_j, arma_scale_k,
                              min_obs, arma_any_na, na_restore,
                              arma_cov);
    parallelFor(0, n_cols, roll_cov_cols);   
  }
  
  // compute rolling eigenvalues and eigenvectors
  RollEigenSlices roll_eigen_slices(arma_cov, n_rows, n_cols - 1,
                                    arma_eigen_values, arma_eigen_vectors);
  parallelFor(0, n_rows, roll_eigen_slices);
  
  // compute rolling principal component regressions
  arma::uvec arma_cols = seq(n_cols - 1);
  RollPcrSlices roll_pcr_slices(arma_cov, n_rows, n_cols, arma_cols, comps,
                                center_x, center_y, arma_center_j,
                                scale_x, scale_y, arma_scale_intercept,
                                arma_eigen_vectors, arma_coef,
                                arma_rsq);
  parallelFor(0, n_rows, roll_pcr_slices);
  
  // create and return a list
  List result = List::create(Named("coefficients") = arma_coef,
                             Named("r.squared") = arma_rsq);
  
  return result;
  
}

// [[Rcpp::export]]
List roll_pcr(const NumericMatrix& x, const NumericMatrix& y,
              const int& width, const arma::uvec& comps,
              const arma::vec& weights, const bool& center_x,
              const bool& center_y, const bool& scale_x,
              const bool& scale_y, const int& min_obs,
              const bool& complete_obs, const bool& na_restore,
              const std::string& parallel_for) {
  
  int n_rows = x.nrow();  
  int n_cols = x.ncol() + 1;
  int y_n_cols = y.ncol();
  List result_z(2);
  List result_coef(y_n_cols);
  List result_rsq(y_n_cols);
  List result(2);
  
  // create a list of matrices,
  // otherwise a list of lists
  if (y_n_cols == 1) {
    
    result_z = roll_pcr_z(x, y(_, 0),
                          width, comps, 
                          weights, center_x,
                          center_y, scale_x,
                          scale_y, min_obs,
                          complete_obs, na_restore, 
                          parallel_for);
    
    arma::mat arma_coef_z = result_z[0];
    arma::mat arma_rsq_z = result_z[1];
    
    // create and return a matrix or xts object for coefficients
    NumericVector coef(wrap(arma_coef_z));
    coef.attr("dim") = IntegerVector::create(n_rows, n_cols);
    List x_dimnames = x.attr("dimnames");
    coef.attr("dimnames") = dimnames_ols(x_dimnames, n_cols - 1);
    coef.attr("index") = x.attr("index");
    coef.attr(".indexCLASS") = x.attr(".indexCLASS");
    coef.attr(".indexTZ") = x.attr(".indexTZ");
    coef.attr("tclass") = x.attr("tclass");
    coef.attr("tzone") = x.attr("tzone");
    coef.attr("class") = x.attr("class");
    
    // create and return a matrix or xts object for r-squareds
    NumericVector rsq(wrap(arma_rsq_z));
    rsq.attr("dim") = IntegerVector::create(n_rows, 1);
    if (x_dimnames.size() > 1) {
      rsq.attr("dimnames") = List::create(x_dimnames[0], "R-squared");
    } else {
      rsq.attr("dimnames") = List::create(R_NilValue, "R-squared");
    }
    rsq.attr("index") = x.attr("index");
    rsq.attr(".indexCLASS") = x.attr(".indexCLASS");
    rsq.attr(".indexTZ") = x.attr(".indexTZ");
    rsq.attr("tclass") = x.attr("tclass");
    rsq.attr("tzone") = x.attr("tzone");
    rsq.attr("class") = x.attr("class");
    
    List y_dimnames = y.attr("dimnames");
    if (y_dimnames.size() > 1) {
      result_coef.attr("names") = y_dimnames[1];
      result_rsq.attr("names") = y_dimnames[1];
    }
    
    // create and return a list
    result = List::create(Named("coefficients") = coef,
                          Named("r.squared") = rsq);
    
  } else {
    
    for (int z = 0; z < y_n_cols; z++) {
      
      result_z = roll_pcr_z(x, y(_, z),
                            width, comps, 
                            weights, center_x,
                            center_y, scale_x,
                            scale_y, min_obs,
                            complete_obs, na_restore, 
                            parallel_for);
      
      arma::mat arma_coef_z = result_z[0];
      arma::mat arma_rsq_z = result_z[1];
      
      // create and return a matrix or xts object for coefficients
      NumericVector coef(wrap(arma_coef_z));
      coef.attr("dim") = IntegerVector::create(n_rows, n_cols);
      List x_dimnames = x.attr("dimnames");
      coef.attr("dimnames") = dimnames_ols(x_dimnames, n_cols - 1);
      coef.attr("index") = x.attr("index");
      coef.attr(".indexCLASS") = x.attr(".indexCLASS");
      coef.attr(".indexTZ") = x.attr(".indexTZ");
      coef.attr("tclass") = x.attr("tclass");
      coef.attr("tzone") = x.attr("tzone");
      coef.attr("class") = x.attr("class");
      
      // create and return a matrix or xts object for r-squareds
      NumericVector rsq(wrap(arma_rsq_z));
      rsq.attr("dim") = IntegerVector::create(n_rows, 1);
      if (x_dimnames.size() > 1) {
        rsq.attr("dimnames") = List::create(x_dimnames[0], "R-squared");
      } else {
        rsq.attr("dimnames") = List::create(R_NilValue, "R-squared");
      }
      rsq.attr("index") = x.attr("index");
      rsq.attr(".indexCLASS") = x.attr(".indexCLASS");
      rsq.attr(".indexTZ") = x.attr(".indexTZ");
      rsq.attr("tclass") = x.attr("tclass");
      rsq.attr("tzone") = x.attr("tzone");
      rsq.attr("class") = x.attr("class");
      
      result_coef(z) = coef;
      result_rsq(z) = rsq;
      
    }
    
    // add names to each list
    List y_dimnames = y.attr("dimnames");
    if (y_dimnames.size() > 1) {
      result_coef.attr("names") = y_dimnames[1];
      result_rsq.attr("names") = y_dimnames[1];
    }
    
    // create and return a list
    result = List::create(Named("coefficients") = result_coef,
                          Named("r.squared") = result_rsq);
    
  }
  
  return result;
  
}

// 'Worker' function for rolling OLS regressions
struct RollLmVifSlices : public Worker {
  
  const arma::cube arma_cov;  // source
  const int n_rows;
  const int n_cols;
  const arma::uvec arma_cols;
  arma::mat& arma_vif;        // destination (pass by reference)
  
  // initialize with source and destination
  RollLmVifSlices(const arma::cube arma_cov, const int n_rows,
                  const int n_cols, const arma::uvec arma_cols,
                  arma::mat& arma_vif)
    : arma_cov(arma_cov), n_rows(n_rows),
      n_cols(n_cols), arma_cols(arma_cols),
      arma_vif(arma_vif) { }
  
  // function call operator that iterates by slice
  void operator()(std::size_t begin_slice, std::size_t end_slice) {
    for (std::size_t i = begin_slice; i < end_slice; i++) {
      
      arma::mat x = arma_cov.slice(i);
      
      int j = 0;
      int k = 0;
      bool any_na = false;
      
      // check if missing value is present
      while ((!any_na) && (j * k < (n_cols - 1) * (n_cols - 1))) {
        for (j = 0; j < n_cols; j++) {
          for (k = 0; k < n_cols; k++) {
            if (std::isnan(x(j, k)))
              any_na = true;
          }
        }
      }
      
      // don't compute if missing value 
      if (!any_na) {
        
        for (int j = 0; j < n_cols; j++) {
          
          arma::mat x = arma_cov.slice(i);
          x.swap_cols(j, n_cols - 1);
          x.swap_rows(j, n_cols - 1);
          
          arma::mat A = x.submat(0, 0, n_cols - 2, n_cols - 2);
          arma::mat b = x.submat(0, n_cols - 1, n_cols - 2, n_cols - 1);
          arma::vec coef(n_cols - 1);
          
          // check if solution is found
          bool status = arma::solve(coef, A, b);
          
          // don't find approximate solution for rank deficient system
          if (status) {
            
            // r-squared
            double rsq = as_scalar((trans(coef) * A * coef) /
                                     x.submat(n_cols - 1, n_cols - 1, n_cols - 1, n_cols - 1));
            
            arma_vif(i, j) = 1 / (1 - rsq);
            
            
          } else if (!status) {
            arma_vif(i, j) = NA_REAL;
          }
          
        }
        
      } else {
        
        arma::vec no_solution(n_cols);
        no_solution.fill(NA_REAL);
        
        arma_vif.row(i) = trans(no_solution);
        
      }
      
    }
  }
  
};

// 'Worker' function for rolling OLS regressions
struct RollPcrVifSlices : public Worker {
  
  const arma::cube arma_cov;    // source
  const int n_rows;
  const int n_cols;
  const arma::uvec arma_cols;
  const arma::uvec arma_comps;  
  arma::mat& arma_vif;          // destination (pass by reference)
  
  // initialize with source and destination
  RollPcrVifSlices(const arma::cube arma_cov, const int n_rows,
                   const int n_cols, const arma::uvec arma_cols,
                   const arma::uvec arma_comps, arma::mat& arma_vif)
    : arma_cov(arma_cov), n_rows(n_rows),
      n_cols(n_cols), arma_cols(arma_cols),
      arma_comps(arma_comps), arma_vif(arma_vif) { }
  
  // function call operator that iterates by slice
  void operator()(std::size_t begin_slice, std::size_t end_slice) {
    for (std::size_t i = begin_slice; i < end_slice; i++) {
      
      arma::mat x = arma_cov.slice(i);
      
      int j = 0;
      int k = 0;
      bool any_na = false;
      
      // check if missing value is present
      while ((!any_na) && (j * k < (n_cols - 1) * (n_cols - 1))) {
        for (j = 0; j < n_cols; j++) {
          for (k = 0; k < n_cols; k++) {
            if (std::isnan(x(j, k)))
              any_na = true;
          }
        }
      }
      
      // don't compute if missing value 
      if (!any_na) {
        
        for (int j = 0; j < n_cols; j++) {
          
          arma::mat x = arma_cov.slice(i);
          x.swap_cols(j, n_cols - 1);
          x.swap_rows(j, n_cols - 1);
          
          arma::mat A = x.submat(0, 0, n_cols - 2, n_cols - 2);
          arma::mat b = x.submat(0, n_cols - 1, n_cols - 2, n_cols - 1);
          arma::vec eigen_values(n_cols - 1);
          arma::mat eigen_vectors(n_cols - 1, n_cols - 1);
          arma::vec gamma(n_cols - 1);
          
          // check if eig_sym solution is found
          bool status1 = arma::eig_sym(eigen_values, eigen_vectors, A);
          
          // don't find approximate solution for rank deficient system
          if (status1) {
            
            // reverse order for consistency with R's eigen
            std::reverse(eigen_values.begin(), eigen_values.end());
            eigen_vectors = arma::fliplr(eigen_vectors);
            
            // check if solution is found
            bool status2 = arma::solve(gamma, A * eigen_vectors, b);
            
            // don't find approximate solution for rank deficient system
            if (status2) {
              
              arma::vec gamma_subset = gamma(arma_comps - 1);
              arma::vec coef = eigen_vectors.submat(arma_cols, arma_comps - 1) * gamma_subset;
              
              // r-squared
              double rsq = as_scalar((trans(coef) * A * coef) /
                                       x.submat(n_cols - 1, n_cols - 1, n_cols - 1, n_cols - 1));
              
              arma_vif(i, j) = 1 / (1 - rsq);
              
            } else if (!status2) {
              arma_vif(i, j) = NA_REAL;
            }
            
          } else if (!status1) {
            arma_vif(i, j) = NA_REAL;
          }
          
        }
        
      } else {
        
        arma::vec no_solution(n_cols);
        no_solution.fill(NA_REAL);
        
        arma_vif.row(i) = trans(no_solution);
        
      }
      
    }
  }
  
};

// [[Rcpp::export]]
NumericMatrix roll_vif(const NumericMatrix& data, const int& width,
                       const arma::uvec& comps, const arma::vec& weights,
                       const bool& center, const bool& scale,
                       const int& min_obs, const bool& complete_obs,
                       const bool& na_restore, const std::string& parallel_for) {
  
  int n_rows = data.nrow();
  int n_cols = data.ncol();
  bool center_x = center;
  bool center_y = center;
  bool scale_x = scale;
  bool scale_y = scale;
  arma::uvec arma_any_na(n_rows);
  arma::cube arma_center_j(n_cols, n_cols, n_rows);
  arma::cube arma_center_k(n_cols, n_cols, n_rows);
  arma::cube arma_scale_j(n_cols, n_cols, n_rows);
  arma::cube arma_scale_k(n_cols, n_cols, n_rows);
  arma::mat arma_scale_intercept(n_rows, n_cols);
  arma::cube arma_cov(n_cols, n_cols, n_rows);
  arma::mat arma_vif(n_rows, n_cols);
  
  // check 'x' argument for errors
  check_vif(n_cols);
  
  // check 'width' argument for errors
  check_width(width, n_rows);
  
  // default 'comps' argument is all components
  check_comps_vif(comps, n_cols - 1);
  
  // default 'weights' argument is equal-weighted,
  // otherwise check argument for errors
  check_weights(weights, width);
  
  // default 'min_obs' argument is 'width' (equivalent to 'na.rm = FALSE'),
  // otherwise check argument for errors
  check_min_obs(min_obs, width);
  
  // default 'complete_obs' argument is 'true' (equivalent to 'complete'),
  // otherwise check argument for errors
  if (complete_obs) {
    arma_any_na = any_na(data);
  } else {
    arma_any_na.fill(0);
  }
  
  // default 'center' argument subtracts mean of each variable,
  // otherwise zero is used
  if (center) {
    if (parallel_for == "rows") {
      RollMeanRowsCube roll_mean_rows(data, n_rows, n_cols, width, weights,
                                      min_obs, arma_any_na, na_restore,
                                      arma_center_j, arma_center_k);
      parallelFor(0, n_rows, roll_mean_rows);
    } else if (parallel_for == "cols") {
      RollMeanColsCube roll_mean_cols(data, n_rows, n_cols, width, weights,
                                      min_obs, arma_any_na, na_restore,
                                      arma_center_j, arma_center_k);
      parallelFor(0, n_cols, roll_mean_cols);
    }
  }
  
  // default 'scale' argument is none,
  // otherwise divide by the standard deviation of each variable
  if (scale) {
    if (parallel_for == "rows") {
      RollVarRowsCube roll_var_rows(data, n_rows, n_cols, width, weights,
                                    center_x, center_y, arma_center_j, arma_center_k,
                                    min_obs, arma_any_na, na_restore,
                                    arma_scale_j, arma_scale_k);
      parallelFor(0, n_rows, roll_var_rows);
    } else if (parallel_for == "cols") {
      RollVarColsCube roll_var_cols(data, n_rows, n_cols, width, weights,
                                    center_x, center_y, arma_center_j, arma_center_k,
                                    min_obs, arma_any_na, na_restore,
                                    arma_scale_j, arma_scale_k);
      parallelFor(0, n_cols, roll_var_cols);
    }
  }
  
  // compute rolling covariance matrices
  if (parallel_for == "rows") {
    RollCovRows roll_cov_rows(data, n_rows, n_cols, width, weights,
                              center_x, center_y,
                              arma_center_j, arma_center_k,
                              scale_x, scale_y,
                              arma_scale_j, arma_scale_k,
                              min_obs, arma_any_na, na_restore,
                              arma_cov);
    parallelFor(0, n_rows, roll_cov_rows); 
  } else if (parallel_for == "cols") {
    RollCovCols roll_cov_cols(data, n_rows, n_cols, width, weights,
                              center_x, center_y,
                              arma_center_j, arma_center_k,
                              scale_x, scale_y,
                              arma_scale_j, arma_scale_k,
                              min_obs, arma_any_na, na_restore,
                              arma_cov);
    parallelFor(0, n_cols, roll_cov_cols);   
  }
  
  // compute rolling variance inflation factors
  arma::uvec arma_cols = seq(n_cols - 1);
  if (comps.size() == arma_cols.size()) {
    RollLmVifSlices roll_lm_vif_slices(arma_cov, n_rows, n_cols, arma_cols,
                                       arma_vif);
    parallelFor(0, n_rows, roll_lm_vif_slices);
  } else {
    RollPcrVifSlices roll_pcr_vif_slices(arma_cov, n_rows, n_cols, arma_cols,
                                         comps, arma_vif);
    parallelFor(0, n_rows, roll_pcr_vif_slices);
  }
  
  // create and return a matrix or xts object for variance inflation factors 
  NumericMatrix result(wrap(arma_vif));
  List dimnames = data.attr("dimnames");
  result.attr("dimnames") = dimnames;
  result.attr("index") = data.attr("index");
  result.attr(".indexCLASS") = data.attr(".indexCLASS");
  result.attr(".indexTZ") = data.attr(".indexTZ");
  result.attr("tclass") = data.attr("tclass");
  result.attr("tzone") = data.attr("tzone");
  result.attr("class") = data.attr("class");
  
  return result;
  
}


