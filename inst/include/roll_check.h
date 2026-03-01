#ifndef ROLL_CHECK_H
#define ROLL_CHECK_H

#include <RcppArmadillo.h>
using namespace Rcpp;

namespace roll {

// scalar checks: positive integer (>= 1)
inline void check_pos_int(const int& value, const char* name) {

  if (value < 1) {
    stop("value of '%s' must be greater than or equal to one", name);
  }

}

// scalar checks: bounded double [lower, upper]
inline void check_bounds_double(const double& value, const double& lower,
                                const double& upper, const char* name) {

  if ((value < lower) || (value > upper)) {
    stop("value of '%s' must be between %f and %f", name, lower, upper);
  }

}

// dimension checks: equality
inline void check_rows_equal(const int& n_rows_a, const int& n_rows_b,
                             const char* name_a, const char* name_b) {

  if (n_rows_a != n_rows_b) {
    stop("number of rows in '%s' must equal the number of rows in '%s'", name_a, name_b);
  }

}

// consolidated weights check
inline void check_weights(const int& n_rows, const int& width,
                          const arma::vec& weights, const char* context) {

  if ((int)weights.size() < std::min(width, n_rows)) {
    stop("length of 'weights' must be greater than or equal to the number of rows in %s or 'width'",
          context);
  }

}

// lambda check for online algorithm (returns bool)
inline bool check_lambda(const arma::vec& weights, const int& n_rows_x,
                         const int& width, const bool& online) {
  
  // check if equal-weights
  bool status_eq = all(weights == weights[0]);
  bool status_exp = true;
  
  // check if exponential-weights
  if (!status_eq) {
    
    int n = weights.size();
    long double lambda = 0;
    long double lambda_prev = 0;
    
    // check if constant ratio
    for (int i = 0; (i < n - 1) && status_exp; i++) {
      
      // ratio of weights
      lambda_prev = lambda;
      lambda = weights[n - i - 2] / weights[n - i - 1];
      
      // tolerance for consistency with R's all.equal
      if (((i > 0) && (std::abs(lambda - lambda_prev) > sqrt(arma::datum::eps))) ||
          ((weights[n - i - 2] > weights[n - i - 1]) && (width < n_rows_x)) ||
          (std::isnan(lambda) || (std::isinf(lambda)))) {
        
        status_exp = false;
        
      }
      
    }
    
  }
  
  if (!status_exp && online) {
    warning("'online' is only supported for equal or exponential decay 'weights'");
  }
  
  return status_exp;
  
}

}

#endif