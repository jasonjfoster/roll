#include "roll.h"

void check_width(const int& width) {
  
  if (width < 1) {
    stop("value of 'width' must be greater than zero");
  }
  
}

bool check_weights(const int& n_rows_xy, const int& width,
                   const arma::vec& weights, const bool& online) {
  
  if ((int)weights.size() < std::min(width, n_rows_xy)) {
    stop("length of 'weights' must be greater than or equal to either the number of rows or 'width'");
  }
  
  // check if equal-weights
  bool status_eq = all(weights == weights[0]);
  bool status_exp = true;
  
  // check if exponential-weights
  if (!status_eq) {
    
    int i = 0;
    int n = weights.size();
    long double lambda = 0;
    long double lambda_prev = 0;
    
    // check if constant ratio
    while (status_exp && (i <= (n - 2))) {
      
      // ratio of weights
      lambda_prev = lambda;
      lambda = weights[n - i - 2] / weights[n - i - 1];
      
      // tolerance for consistency with R's all.equal
      if ((i > 0) && (std::abs(lambda - lambda_prev) > sqrt(arma::datum::eps))) {
        status_exp = false;
      }
      
      i += 1;
      
    }
    
  }
  
  if (!status_exp && online) {
    warning("'online' is only supported for equal or exponential 'weights'");
  }
  
  return status_exp;
  
}

void check_min_obs(const int& min_obs) {
  
  if (min_obs < 1) {
    stop("value of 'min_obs' must be greater than zero");
  }
  
}

void check_lm(const int& n_rows_x, const int& n_rows_y) {
  
  if (n_rows_x != n_rows_y) {
    stop("number of rows in 'x' must equal the number of rows in 'y'");
  }
  
}

List x_dimnames_lm(const List& input, const int& n_cols_x, const bool& intercept) {
  
  if (intercept && (input.size() > 1)) {
    
    CharacterVector dimnames_cols = input[1];
    CharacterVector result(n_cols_x);
    result(0) = "(Intercept)";
    
    std::copy(dimnames_cols.begin(), dimnames_cols.end(), result.begin() + 1);
    
    return List::create(input[0], result);
    
  } else if (!intercept && (input.size() > 1)) {
    
    return List::create(input[0], input[1]);
    
  } else if (intercept) {
    
    CharacterVector result(n_cols_x);
    result(0) = "(Intercept)";
    
    for (int i = 1; i < n_cols_x; i++) {
      
      result[i] = "x";
      result[i] += i;
      
    }
    
    return List::create(R_NilValue, result);
    
  } else {
    
    CharacterVector result(n_cols_x);
    
    for (int i = 0; i < n_cols_x; i++) {
      
      result[i] = "x";
      result[i] += i + 1;
      
    }
    
    return List::create(R_NilValue, result);
    
  }
  
}

CharacterVector y_dimnames_lm(const List& input, const int& n_cols_y) {
  
  if (input.size() > 1) {
    
    return input[1];
    
  } else {
    
    CharacterVector result(n_cols_y);
    
    for (int i = 0; i < n_cols_y; i++) {
      
      result[i] = "y";
      result[i] += i + 1;
      
    }
    
    return result;
    
  }
  
}

arma::uvec any_na_x(const NumericMatrix& x) {
  
  int n_rows_x = x.nrow();
  int n_cols_x = x.ncol();
  arma::uvec result(n_rows_x);
  
  for (int i = 0; i < n_rows_x; i++) {
    
    int any_na = 0;
    int j = 0;
    
    while ((any_na == 0) && (j < n_cols_x)) {
      if (std::isnan(x(i, j))) {
        any_na = 1;
      }
      j += 1;
    }
    
    result[i] = any_na;
    
  }
  
  return result;
  
}

arma::uvec any_na_xy(const NumericMatrix& x, const NumericMatrix& y) {
  
  int n_rows_xy = x.nrow(); // check that x and y are the same!
  int n_cols_x = x.ncol();
  int n_cols_y = y.ncol();
  arma::uvec result(n_rows_xy);
  
  for (int i = 0; i < n_rows_xy; i++) {
    
    int any_na = 0;
    int j = 0;
    int k = 0;
    
    while ((any_na == 0) && (j < n_cols_x)) {
      if (std::isnan(x(i, j))) {
        any_na = 1;
      }
      j += 1;
    }
    
    while ((any_na == 0) && (k < n_cols_y)) {
      if (std::isnan(y(i, k))) {
        any_na = 1;
      }
      k += 1;
    }
    
    result[i] = any_na;
    
  }
  
  return result;
  
}

// [[Rcpp::export(.roll_sum)]]
NumericMatrix roll_sum(const NumericMatrix& x, const int& width,
                       const arma::vec& weights, const int& min_obs,
                       const bool& complete_obs, const bool& na_restore,
                       const bool& online) {
  
  int n = weights.size();
  int n_rows_x = x.nrow();
  int n_cols_x = x.ncol();
  arma::uvec arma_any_na(n_rows_x);
  arma::mat arma_sum(n_rows_x, n_cols_x);
  
  // check 'width' argument for errors
  check_width(width);
  
  // default 'weights' argument is equal-weighted,
  // otherwise check argument for errors
  bool status = check_weights(n_rows_x, width, weights, online);
  
  // default 'min_obs' argument is 'width',
  // otherwise check argument for errors
  check_min_obs(min_obs);
  
  // default 'complete_obs' argument is 'false',
  // otherwise check argument for errors
  if (complete_obs) {
    arma_any_na = any_na_x(x);
  } else {
    arma_any_na.fill(0);
  }
  
  // compute rolling sums
  if (status && online) {
    
    RollSumOnline roll_sum_online(x, n, n_rows_x, n_cols_x, width,
                                  weights, min_obs,
                                  arma_any_na, na_restore,
                                  arma_sum);
    parallelFor(0, n_cols_x, roll_sum_online);
    
  } else {
    
    RollSumParallel roll_sum_parallel(x, n, n_rows_x, n_cols_x, width,
                                      weights, min_obs,
                                      arma_any_na, na_restore,
                                      arma_sum);
    parallelFor(0, n_rows_x * n_cols_x, roll_sum_parallel);
    
  }
  
  // create and return a matrix or xts object
  NumericMatrix result(wrap(arma_sum));
  List dimnames = x.attr("dimnames");
  result.attr("dimnames") = dimnames;
  result.attr("index") = x.attr("index");
  result.attr(".indexCLASS") = x.attr(".indexCLASS");
  result.attr(".indexTZ") = x.attr(".indexTZ");
  result.attr("tclass") = x.attr("tclass");
  result.attr("tzone") = x.attr("tzone");
  result.attr("class") = x.attr("class");
  
  return result;
  
}

// [[Rcpp::export(.roll_prod)]]
NumericMatrix roll_prod(const NumericMatrix& x, const int& width,
                        const arma::vec& weights, const int& min_obs,
                        const bool& complete_obs, const bool& na_restore,
                        const bool& online) {
  
  int n = weights.size();
  int n_rows_x = x.nrow();
  int n_cols_x = x.ncol();
  arma::uvec arma_any_na(n_rows_x);
  arma::mat arma_prod(n_rows_x, n_cols_x);
  
  // check 'width' argument for errors
  check_width(width);
  
  // default 'weights' argument is equal-weighted,
  // otherwise check argument for errors
  bool status = check_weights(n_rows_x, width, weights, online);
  
  // default 'min_obs' argument is 'width',
  // otherwise check argument for errors
  check_min_obs(min_obs);
  
  // default 'complete_obs' argument is 'false',
  // otherwise check argument for errors
  if (complete_obs) {
    arma_any_na = any_na_x(x);
  } else {
    arma_any_na.fill(0);
  }
  
  // compute rolling products
  if (status && online) {
    
    RollProdOnline roll_prod_online(x, n, n_rows_x, n_cols_x, width,
                                    weights, min_obs,
                                    arma_any_na, na_restore, 
                                    arma_prod);
    parallelFor(0, n_cols_x, roll_prod_online);
    
  } else {
    
    RollProdParallel roll_prod_parallel(x, n, n_rows_x, n_cols_x, width,
                                        weights, min_obs,
                                        arma_any_na, na_restore,
                                        arma_prod);
    parallelFor(0, n_rows_x * n_cols_x, roll_prod_parallel);
    
  }
  
  // create and return a matrix or xts object
  NumericMatrix result(wrap(arma_prod));
  List dimnames = x.attr("dimnames");
  result.attr("dimnames") = dimnames;
  result.attr("index") = x.attr("index");
  result.attr(".indexCLASS") = x.attr(".indexCLASS");
  result.attr(".indexTZ") = x.attr(".indexTZ");
  result.attr("tclass") = x.attr("tclass");
  result.attr("tzone") = x.attr("tzone");
  result.attr("class") = x.attr("class");
  
  return result;
  
}

// [[Rcpp::export(.roll_mean)]]
NumericMatrix roll_mean(const NumericMatrix& x, const int& width,
                        const arma::vec& weights, const int& min_obs,
                        const bool& complete_obs, const bool& na_restore,
                        const bool& online) {
  
  int n = weights.size();
  int n_rows_x = x.nrow();
  int n_cols_x = x.ncol();
  arma::uvec arma_any_na(n_rows_x);
  arma::mat arma_mean(n_rows_x, n_cols_x);
  
  // check 'width' argument for errors
  check_width(width);
  
  // default 'weights' argument is equal-weighted,
  // otherwise check argument for errors
  bool status = check_weights(n_rows_x, width, weights, online);
  
  // default 'min_obs' argument is 'width',
  // otherwise check argument for errors
  check_min_obs(min_obs);
  
  // default 'complete_obs' argument is 'false',
  // otherwise check argument for errors
  if (complete_obs) {
    arma_any_na = any_na_x(x);
  } else {
    arma_any_na.fill(0);
  }
  
  // compute rolling means
  if (status && online) {
    
    RollMeanOnline roll_mean_online(x, n, n_rows_x, n_cols_x, width,
                                    weights, min_obs,
                                    arma_any_na, na_restore,
                                    arma_mean);
    parallelFor(0, n_cols_x, roll_mean_online);
    
  } else {
    
    RollMeanParallel roll_mean_parallel(x, n, n_rows_x, n_cols_x, width,
                                        weights, min_obs,
                                        arma_any_na, na_restore,
                                        arma_mean);
    parallelFor(0, n_rows_x * n_cols_x, roll_mean_parallel);
    
  }
  
  // create and return a matrix or xts object
  NumericMatrix result(wrap(arma_mean));
  List dimnames = x.attr("dimnames");
  result.attr("dimnames") = dimnames;
  result.attr("index") = x.attr("index");
  result.attr(".indexCLASS") = x.attr(".indexCLASS");
  result.attr(".indexTZ") = x.attr(".indexTZ");
  result.attr("tclass") = x.attr("tclass");
  result.attr("tzone") = x.attr("tzone");
  result.attr("class") = x.attr("class");
  
  return result;
  
}

// [[Rcpp::export(.roll_var)]]
NumericMatrix roll_var(const NumericMatrix& x, const int& width,
                       const arma::vec& weights, const bool& center,
                       const int& min_obs, const bool& complete_obs,
                       const bool& na_restore, const bool& online) {
  
  int n = weights.size();
  int n_rows_x = x.nrow();
  int n_cols_x = x.ncol();
  arma::uvec arma_any_na(n_rows_x);
  arma::mat arma_var(n_rows_x, n_cols_x);
  
  // check 'width' argument for errors
  check_width(width);
  
  // default 'weights' argument is equal-weighted,
  // otherwise check argument for errors
  bool status = check_weights(n_rows_x, width, weights, online);
  
  // default 'min_obs' argument is 'width',
  // otherwise check argument for errors
  check_min_obs(min_obs);
  
  // default 'complete_obs' argument is 'false',
  // otherwise check argument for errors
  if (complete_obs) {
    arma_any_na = any_na_x(x);
  } else {
    arma_any_na.fill(0);
  }
  
  // compute rolling variances
  if (status && online) {
    
    RollVarOnline roll_var_online(x, n, n_rows_x, n_cols_x, width,
                                  weights, center, min_obs,
                                  arma_any_na, na_restore,
                                  arma_var);
    parallelFor(0, n_cols_x, roll_var_online);
    
  } else {
    
    RollVarParallel roll_var_parallel(x, n, n_rows_x, n_cols_x, width,
                                      weights, center, min_obs,
                                      arma_any_na, na_restore,
                                      arma_var);
    parallelFor(0, n_rows_x * n_cols_x, roll_var_parallel);
    
  }
  
  // create and return a matrix or xts object
  NumericMatrix result(wrap(arma_var));
  List dimnames = x.attr("dimnames");
  result.attr("dimnames") = dimnames;
  result.attr("index") = x.attr("index");
  result.attr(".indexCLASS") = x.attr(".indexCLASS");
  result.attr(".indexTZ") = x.attr(".indexTZ");
  result.attr("tclass") = x.attr("tclass");
  result.attr("tzone") = x.attr("tzone");
  result.attr("class") = x.attr("class");
  
  return result;
  
}

// [[Rcpp::export(.roll_sd)]]
NumericMatrix roll_sd(const NumericMatrix& x, const int& width,
                      const arma::vec& weights, const bool& center,
                      const int& min_obs, const bool& complete_obs,
                      const bool& na_restore, const bool& online) {
  
  int n = weights.size();
  int n_rows_x = x.nrow();
  int n_cols_x = x.ncol();
  arma::uvec arma_any_na(n_rows_x);
  arma::mat arma_sd(n_rows_x, n_cols_x);
  
  // check 'width' argument for errors
  check_width(width);
  
  // default 'weights' argument is equal-weighted,
  // otherwise check argument for errors
  bool status = check_weights(n_rows_x, width, weights, online);
  
  // default 'min_obs' argument is 'width',
  // otherwise check argument for errors
  check_min_obs(min_obs);
  
  // default 'complete_obs' argument is 'false',
  // otherwise check argument for errors
  if (complete_obs) {
    arma_any_na = any_na_x(x);
  } else {
    arma_any_na.fill(0);
  }
  
  // compute rolling standard deviations
  if (status && online) {
    
    RollVarOnline roll_sd_online(x, n, n_rows_x, n_cols_x, width,
                                 weights, center, min_obs,
                                 arma_any_na, na_restore,
                                 arma_sd);
    parallelFor(0, n_cols_x, roll_sd_online);
    
  } else {
    
    RollVarParallel roll_sd_parallel(x, n, n_rows_x, n_cols_x, width,
                                     weights, center, min_obs,
                                     arma_any_na, na_restore,
                                     arma_sd);
    parallelFor(0, n_rows_x * n_cols_x, roll_sd_parallel);
    
  }
  
  // create and return a matrix or xts object
  NumericMatrix result(wrap(arma_sd));
  List dimnames = x.attr("dimnames");
  result.attr("dimnames") = dimnames;
  result.attr("index") = x.attr("index");
  result.attr(".indexCLASS") = x.attr(".indexCLASS");
  result.attr(".indexTZ") = x.attr(".indexTZ");
  result.attr("tclass") = x.attr("tclass");
  result.attr("tzone") = x.attr("tzone");
  result.attr("class") = x.attr("class");
  
  return result;
  
}

// [[Rcpp::export(.roll_scale)]]
NumericMatrix roll_scale(const NumericMatrix& x, const int& width,
                         const arma::vec& weights, const bool& center,
                         const bool& scale, const int& min_obs,
                         const bool& complete_obs, const bool& na_restore,
                         const bool& online) {
  
  int n = weights.size();
  int n_rows_x = x.nrow();
  int n_cols_x = x.ncol();
  arma::uvec arma_any_na(n_rows_x);
  arma::mat arma_scale(n_rows_x, n_cols_x);
  
  // check 'width' argument for errors
  check_width(width);
  
  // default 'weights' argument is equal-weighted,
  // otherwise check argument for errors
  bool status = check_weights(n_rows_x, width, weights, online);
  
  // default 'min_obs' argument is 'width',
  // otherwise check argument for errors
  check_min_obs(min_obs);
  
  // default 'complete_obs' argument is 'false',
  // otherwise check argument for errors
  if (complete_obs) {
    arma_any_na = any_na_x(x);
  } else {
    arma_any_na.fill(0);
  }
  
  // compute rolling centering and scaling
  if (status && online) {
    
    RollScaleOnline roll_scale_online(x, n, n_rows_x, n_cols_x, width,
                                      weights, center, scale, min_obs,
                                      arma_any_na, na_restore,
                                      arma_scale);
    parallelFor(0, n_cols_x, roll_scale_online);
    
  } else {
    
    RollScaleParallel roll_scale_parallel(x, n, n_rows_x, n_cols_x, width,
                                          weights, center, scale, min_obs,
                                          arma_any_na, na_restore,
                                          arma_scale);
    parallelFor(0, n_rows_x * n_cols_x, roll_scale_parallel);
    
  }
  
  // create and return a matrix or xts object
  NumericMatrix result(wrap(arma_scale));
  List dimnames = x.attr("dimnames");
  result.attr("dimnames") = dimnames;
  result.attr("index") = x.attr("index");
  result.attr(".indexCLASS") = x.attr(".indexCLASS");
  result.attr(".indexTZ") = x.attr(".indexTZ");
  result.attr("tclass") = x.attr("tclass");
  result.attr("tzone") = x.attr("tzone");
  result.attr("class") = x.attr("class");
  
  return result;
  
}

NumericVector roll_cov_z(const NumericMatrix& x, const NumericMatrix& y,
                         const int& width, const arma::vec& weights,
                         const bool& center, const bool& scale,
                         const int& min_obs, const bool& complete_obs,
                         const bool& na_restore, const bool& online,
                         const bool& symmetric) {
  
  int n = weights.size();
  int n_rows_xy = x.nrow();
  int n_cols_x = x.ncol();
  int n_cols_y = y.ncol();
  arma::uvec arma_any_na(n_rows_xy);
  arma::cube arma_cov(n_cols_x, n_cols_y, n_rows_xy);
  
  // check 'x' and 'y' arguments for errors
  check_lm(n_rows_xy, y.nrow());
  
  // check 'width' argument for errors
  check_width(width);
  
  // default 'weights' argument is equal-weighted,
  // otherwise check argument for errors
  bool status = check_weights(n_rows_xy, width, weights, online);
  
  // default 'min_obs' argument is 'width',
  // otherwise check argument for errors
  check_min_obs(min_obs);
  
  // default 'complete_obs' argument is 'true',
  // otherwise check argument for errors
  if (complete_obs && symmetric) {
    arma_any_na = any_na_x(x);
  } else if (complete_obs && !symmetric) {
    arma_any_na = any_na_xy(x, y);
  } else {
    arma_any_na.fill(0);
  }
  
  // compute rolling covariances
  if (status && online) {
    
    if (symmetric) {
      
      // y is null
      RollCovOnlineXX roll_cov_online(x, n, n_rows_xy, n_cols_x, width,
                                      weights, center, scale, min_obs,
                                      arma_any_na, na_restore,
                                      arma_cov);
      parallelFor(0, n_cols_x, roll_cov_online);
      
    } else if (!symmetric) {
      
      // y is not null
      RollCovOnlineXY roll_cov_online(x, y, n, n_rows_xy, n_cols_x, n_cols_y, width,
                                      weights, center, scale, min_obs,
                                      arma_any_na, na_restore,
                                      arma_cov);
      parallelFor(0, n_cols_x, roll_cov_online);
      
    }
    
  } else {
    
    if (symmetric) {
      
      // y is null
      RollCovParallelXX roll_cov_parallel(x, n, n_rows_xy, n_cols_x, width,
                                          weights, center, scale, min_obs,
                                          arma_any_na, na_restore,
                                          arma_cov);
      parallelFor(0, n_rows_xy * n_cols_x * (n_cols_x + 1) / 2, roll_cov_parallel);
      
    } else if (!symmetric) {
      
      // y is not null
      RollCovParallelXY roll_cov_parallel(x, y, n, n_rows_xy, n_cols_x, n_cols_y, width,
                                          weights, center, scale, min_obs,
                                          arma_any_na, na_restore,
                                          arma_cov);
      parallelFor(0, n_rows_xy * n_cols_x * n_cols_y, roll_cov_parallel);
      
    }
    
  }
  
  // create and return a matrix
  NumericVector result(wrap(arma_cov));
  result.attr("dim") = IntegerVector::create(n_cols_x, n_cols_y, n_rows_xy);
  List x_dimnames = x.attr("dimnames");
  List y_dimnames = y.attr("dimnames");
  if ((x_dimnames.size() > 1) && (y_dimnames.size() > 1)) {
    result.attr("dimnames") = List::create(x_dimnames[1], y_dimnames[1]);
  } else if (x_dimnames.size() > 1) {
    result.attr("dimnames") = List::create(x_dimnames[1], R_NilValue);
  } else if (y_dimnames.size() > 1) {
    result.attr("dimnames") = List::create(R_NilValue, y_dimnames[1]);
  }
  
  return result;
  
}

// [[Rcpp::export(.roll_cov)]]
NumericVector roll_cov(const NumericMatrix& x, const Nullable<NumericMatrix>& y,
                       const int& width, const arma::vec& weights,
                       const bool& center, const bool& scale,
                       const int& min_obs, const bool& complete_obs,
                       const bool& na_restore, const bool& online) {
  
  if (y.isNotNull()) {
    
    NumericMatrix yy(y);
    return roll_cov_z(x, yy, width, weights, center, scale, min_obs, complete_obs, 
                      na_restore, online, false);
    
  }
  
  return roll_cov_z(x, x, width, weights, center, scale, min_obs, complete_obs,
                    na_restore, online, true);
  
}

List roll_lm_z(const NumericMatrix& x, const NumericVector& y,
               const int& width, const arma::vec& weights,
               const bool& intercept, const int& min_obs,
               const bool& complete_obs, const bool& na_restore,
               const bool& online) {
  
  int n = weights.size();
  int n_rows_xy = x.nrow();
  int n_cols_x = x.ncol() + 1;
  arma::uvec arma_any_na(n_rows_xy);
  arma::vec arma_n_obs(n_rows_xy);
  arma::vec arma_sum_w(n_rows_xy);
  arma::mat arma_mean(n_rows_xy, n_cols_x);
  arma::cube arma_cov(n_cols_x, n_cols_x, n_rows_xy);
  arma::vec arma_rsq(n_rows_xy);
  List result(3);
  
  // check 'x' and 'y' arguments for errors
  check_lm(n_rows_xy, y.size());
  
  // check 'width' argument for errors
  check_width(width);
  
  // default 'weights' argument is equal-weighted,
  // otherwise check argument for errors
  bool status = check_weights(n_rows_xy, width, weights, online);
  
  // default 'min_obs' argument is 'width',
  // otherwise check argument for errors
  check_min_obs(min_obs);
  
  // cbind x and y variables
  NumericMatrix data(n_rows_xy, n_cols_x);
  std::copy(x.begin(), x.end(), data.begin());
  std::copy(y.begin(), y.end(), data.begin() + n_rows_xy * (n_cols_x - 1));
  
  // default 'complete_obs' argument is 'true',
  // otherwise check argument for errors
  if (complete_obs) {
    arma_any_na = any_na_x(data);
  } else {
    
    warning("pairwise regression is not currently supported");
    arma_any_na = any_na_x(data);
    
  }
  
  // compute rolling covariances
  if (status && online) {
    
    RollCovOnlineLm roll_cov_online(data, n, n_rows_xy, n_cols_x, width,
                                    weights, intercept,
                                    min_obs, arma_any_na, na_restore,
                                    arma_n_obs, arma_sum_w, arma_mean,
                                    arma_cov);
    parallelFor(0, n_cols_x, roll_cov_online);
    
  } else {
    
    RollCovParallelLm roll_cov_parallel(data, n, n_rows_xy, n_cols_x, width,
                                        weights, intercept,
                                        min_obs, arma_any_na, na_restore,
                                        arma_n_obs, arma_sum_w, arma_mean,
                                        arma_cov);
    parallelFor(0, n_rows_xy * n_cols_x * (n_cols_x + 1) / 2, roll_cov_parallel);
    
  }
  
  // compute rolling linear models
  if (intercept) {
    
    arma::mat arma_coef(n_rows_xy, n_cols_x);
    arma::mat arma_se(n_rows_xy, n_cols_x);
    RollLmInterceptTRUE roll_lm_slices(arma_cov, n, n_rows_xy, n_cols_x, width,
                                       arma_n_obs, arma_sum_w, arma_mean,
                                       arma_coef, arma_rsq, arma_se);
    parallelFor(0, n_rows_xy, roll_lm_slices);
    
    result = List::create(Named("coefficients") = arma_coef,
                          Named("r.squared") = arma_rsq,
                          Named("std.error") = arma_se);
    
  } else if (!intercept) {
    
    arma::mat arma_coef(n_rows_xy, n_cols_x - 1);
    arma::mat arma_se(n_rows_xy, n_cols_x - 1);
    RollLmInterceptFALSE roll_lm_slices(arma_cov, n, n_rows_xy, n_cols_x, width,
                                        arma_n_obs, arma_sum_w,
                                        arma_coef, arma_rsq, arma_se);
    parallelFor(0, n_rows_xy, roll_lm_slices);
    
    // create and return a list
    result = List::create(Named("coefficients") = arma_coef,
                          Named("r.squared") = arma_rsq,
                          Named("std.error") = arma_se);
    
  }
  
  return result;
  
}

// [[Rcpp::export(.roll_lm)]]
List roll_lm(const NumericMatrix& x, const NumericMatrix& y,
             const int& width, const arma::vec& weights,
             const bool& intercept, const int& min_obs,
             const bool& complete_obs, const bool& na_restore,
             const bool& online) {
  
  int n_rows_xy = x.nrow();
  int n_cols_x = x.ncol();
  int n_cols_y = y.ncol();
  List result_coef(n_cols_y);
  List result_rsq(n_cols_y);
  List result_se(n_cols_y);
  List result_z(3);
  List result(3);
  
  if (intercept) {
    n_cols_x += 1;
  }
  
  // create a list of matrices,
  // otherwise a list of lists
  if (n_cols_y == 1) {
    
    result_z = roll_lm_z(x, y(_, 0),
                         width, weights,
                         intercept, min_obs,
                         complete_obs, na_restore,
                         online);
    
    arma::mat arma_coef_z = result_z[0];
    arma::mat arma_rsq_z = result_z[1];
    arma::mat arma_se_z = result_z[2];
    
    // create and return a matrix or xts object for coefficients
    NumericVector coef(wrap(arma_coef_z));
    coef.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
    List x_dimnames = x.attr("dimnames");
    coef.attr("dimnames") = x_dimnames_lm(x_dimnames, n_cols_x, intercept);
    coef.attr("index") = x.attr("index");
    coef.attr(".indexCLASS") = x.attr(".indexCLASS");
    coef.attr(".indexTZ") = x.attr(".indexTZ");
    coef.attr("tclass") = x.attr("tclass");
    coef.attr("tzone") = x.attr("tzone");
    coef.attr("class") = x.attr("class");
    
    // create and return a matrix or xts object for r-squareds
    NumericVector rsq(wrap(arma_rsq_z));
    rsq.attr("dim") = IntegerVector::create(n_rows_xy, 1);
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
    
    // create and return a matrix or xts object for standard errors
    NumericVector se(wrap(arma_se_z));
    se.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
    se.attr("dimnames") = coef.attr("dimnames");
    se.attr("index") = x.attr("index");
    se.attr(".indexCLASS") = x.attr(".indexCLASS");
    se.attr(".indexTZ") = x.attr(".indexTZ");
    se.attr("tclass") = x.attr("tclass");
    se.attr("tzone") = x.attr("tzone");
    se.attr("class") = x.attr("class");
    
    // create and return a list
    result = List::create(Named("coefficients") = coef,
                          Named("r.squared") = rsq,
                          Named("std.error") = se);
    
  } else {
    
    for (int z = 0; z < n_cols_y; z++) {
      
      result_z = roll_lm_z(x, y(_, z),
                           width, weights,
                           intercept, min_obs,
                           complete_obs, na_restore,
                           online);
      
      arma::mat arma_coef_z = result_z[0];
      arma::mat arma_rsq_z = result_z[1];
      arma::mat arma_se_z = result_z[2];
      
      // create and return a matrix or xts object for coefficients
      NumericVector coef(wrap(arma_coef_z));
      coef.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
      List x_dimnames = x.attr("dimnames");
      coef.attr("dimnames") = x_dimnames_lm(x_dimnames, n_cols_x, intercept);
      coef.attr("index") = x.attr("index");
      coef.attr(".indexCLASS") = x.attr(".indexCLASS");
      coef.attr(".indexTZ") = x.attr(".indexTZ");
      coef.attr("tclass") = x.attr("tclass");
      coef.attr("tzone") = x.attr("tzone");
      coef.attr("class") = x.attr("class");
      
      // create and return a matrix or xts object for r-squareds
      NumericVector rsq(wrap(arma_rsq_z));
      rsq.attr("dim") = IntegerVector::create(n_rows_xy, 1);
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
      
      // create and return a matrix or xts object for standard errors
      NumericVector se(wrap(arma_se_z));
      se.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
      se.attr("dimnames") = coef.attr("dimnames");
      se.attr("index") = x.attr("index");
      se.attr(".indexCLASS") = x.attr(".indexCLASS");
      se.attr(".indexTZ") = x.attr(".indexTZ");
      se.attr("tclass") = x.attr("tclass");
      se.attr("tzone") = x.attr("tzone");
      se.attr("class") = x.attr("class");
      
      result_coef(z) = coef;
      result_rsq(z) = rsq;
      result_se(z) = se;
      
    }
    
    // add names to each list
    List y_dimnames = y.attr("dimnames");
    result_coef.attr("names") = y_dimnames_lm(y_dimnames, n_cols_y);
    result_rsq.attr("names") = result_coef.attr("names");
    result_se.attr("names") = result_coef.attr("names");
    
    // create and return a list
    result = List::create(Named("coefficients") = result_coef,
                          Named("r.squared") = result_rsq,
                          Named("std.error") = result_se);
    
    
  }
  
  return result;
  
}
