#include "roll.h"

List dimnames_lm_x(const List& input, const int& n_cols_x,
                   const bool& intercept) {
  
  CharacterVector result(n_cols_x);
  
  if (input.size() > 1) {
    
    CharacterVector dimnames_cols = input[1];
    
    if (intercept) {

      result(0) = "(Intercept)";
      
      std::copy(dimnames_cols.begin(), dimnames_cols.end(), result.begin() + 1);
      
      return List::create(input[0], result);
      
    } else {
      return List::create(input[0], dimnames_cols);
    }
    
  } else {
    
    if (intercept) {
      
      result(0) = "(Intercept)";
      
      for (int i = 1; i < n_cols_x; i++) {
        
        result[i] = "x";
        result[i] += i;
        
      }
      
      return List::create(R_NilValue, result);
      
    } else {
      
      for (int i = 0; i < n_cols_x; i++) {
        
        result[i] = "x";
        result[i] += i + 1;
        
      }
      
      return List::create(R_NilValue, result);
      
    }
    
  }
  
}

CharacterVector dimnames_lm_y(const List& input, const int& n_cols_y) {
  
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

arma::uvec any_na_i(const IntegerMatrix& x) {
  
  int n_rows_x = x.nrow();
  int n_cols_x = x.ncol();
  arma::uvec result(n_rows_x);
  
  for (int i = 0; i < n_rows_x; i++) {
    
    int any_na = 0;
    
    for (int j = 0; (j < n_cols_x) && (any_na == 0); j++) {
      if (x(i, j) == NA_INTEGER) {
        any_na = 1;
      }
    }
    
    result[i] = any_na;
    
  }
  
  return result;
  
}

arma::uvec any_na_x(const NumericMatrix& x) {
  
  int n_rows_x = x.nrow();
  int n_cols_x = x.ncol();
  arma::uvec result(n_rows_x);
  
  for (int i = 0; i < n_rows_x; i++) {
    
    int any_na = 0;
    
    for (int j = 0; (j < n_cols_x) && (any_na == 0); j++) {
      if (std::isnan(x(i, j))) {
        any_na = 1;
      }
    }
    
    result[i] = any_na;
    
  }
  
  return result;
  
}

arma::uvec any_na_xy(const NumericMatrix& x, const NumericMatrix& y) {
  
  int n_rows_xy = x.nrow();
  int n_cols_x = x.ncol();
  int n_cols_y = y.ncol();
  arma::uvec result(n_rows_xy);
  
  for (int i = 0; i < n_rows_xy; i++) {
    
    int any_na = 0;
    
    for (int j = 0; (j < n_cols_x) && (any_na == 0); j++) {
      if (std::isnan(x(i, j))) {
        any_na = 1;
      }
    }
    
    for (int k = 0; (k < n_cols_y) && (any_na == 0); k++) {
      if (std::isnan(y(i, k))) {
        any_na = 1;
      }
    }
    
    result[i] = any_na;
    
  }
  
  return result;
  
}

// [[Rcpp::export(.roll_any)]]
SEXP roll_any(const SEXP& x, const int& width,
              const int& min_obs, const bool& complete_obs,
              const bool& na_restore, const bool& online) {
  
  if (Rf_isMatrix(x)) {
    
    IntegerMatrix xx(x); // RMatrix<bool> is not supported
    int n_rows_x = xx.nrow();
    int n_cols_x = xx.ncol();
    arma::uvec arma_any_na(n_rows_x);
    IntegerMatrix rcpp_any(n_rows_x, n_cols_x);
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // default 'complete_obs' argument is 'false',
    // otherwise check argument for errors
    if (complete_obs) {
      arma_any_na = any_na_i(x);
    } else {
      arma_any_na.fill(0);
    }
    
    // compute rolling any
    if (online) {
      
      roll::RollAnyOnlineMat roll_any_online(xx, n_rows_x, n_cols_x, width,
                                             min_obs, arma_any_na, na_restore,
                                             rcpp_any);
      parallelFor(0, n_cols_x, roll_any_online);
      
    } else {
      
      roll::RollAnyOfflineMat roll_any_offline(xx, n_rows_x, n_cols_x, width,
                                               min_obs, arma_any_na, na_restore,
                                               rcpp_any);
      parallelFor(0, n_rows_x * n_cols_x, roll_any_offline);
      
    }
    
    // create and return a matrix or xts object
    LogicalMatrix result(rcpp_any);
    roll::xts_attr(result, xx);
    
    return result;
    
  } else {
    
    LogicalVector xx(x);
    int n_rows_x = xx.size();
    IntegerVector rcpp_x(xx);
    IntegerVector rcpp_any(n_rows_x);
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // compute rolling any
    if (online) {
      
      roll::RollAnyOnlineVec roll_any_online(rcpp_x, n_rows_x, width,
                                             min_obs, na_restore,
                                             rcpp_any);
      roll_any_online();
      
    } else {
      
      roll::RollAnyOfflineVec roll_any_offline(rcpp_x, n_rows_x, width,
                                               min_obs, na_restore,
                                               rcpp_any);
      parallelFor(0, n_rows_x, roll_any_offline);
      
    }
    
    // create and return a vector object
    LogicalVector result(wrap(rcpp_any));
    roll::vec_attr(result, xx);
    
    return result;
    
  }
  
}

// [[Rcpp::export(.roll_all)]]
SEXP roll_all(const SEXP& x, const int& width,
              const int& min_obs, const bool& complete_obs,
              const bool& na_restore, const bool& online) {
  
  if (Rf_isMatrix(x)) {
    
    IntegerMatrix xx(x); // RMatrix<bool> is not supported
    int n_rows_x = xx.nrow();
    int n_cols_x = xx.ncol();
    arma::uvec arma_any_na(n_rows_x);
    IntegerMatrix rcpp_all(n_rows_x, n_cols_x);
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // default 'complete_obs' argument is 'false',
    // otherwise check argument for errors
    if (complete_obs) {
      arma_any_na = any_na_i(x);
    } else {
      arma_any_na.fill(0);
    }
    
    // compute rolling all
    if (online) {
      
      roll::RollAllOnlineMat roll_all_online(xx, n_rows_x, n_cols_x, width,
                                             min_obs, arma_any_na, na_restore,
                                             rcpp_all);
      parallelFor(0, n_cols_x, roll_all_online);
      
    } else {
      
      roll::RollAllOfflineMat roll_all_offline(xx, n_rows_x, n_cols_x, width,
                                               min_obs, arma_any_na, na_restore,
                                               rcpp_all);
      parallelFor(0, n_rows_x * n_cols_x, roll_all_offline);
      
    }
    
    // create and return a matrix or xts object
    LogicalMatrix result(rcpp_all);
    roll::xts_attr(result, xx);
    
    return result;
    
  } else {
    
    LogicalVector xx(x);
    int n_rows_x = xx.size();
    IntegerVector rcpp_x(xx);
    IntegerVector rcpp_all(n_rows_x);
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // compute rolling all
    if (online) {
      
      roll::RollAllOnlineVec roll_all_online(rcpp_x, n_rows_x, width,
                                             min_obs, na_restore,
                                             rcpp_all);
      roll_all_online();
      
    } else {
      
      roll::RollAllOfflineVec roll_all_offline(rcpp_x, n_rows_x, width,
                                               min_obs, na_restore,
                                               rcpp_all);
      parallelFor(0, n_rows_x, roll_all_offline);
      
    }
    
    // create and return a vector object
    LogicalVector result(wrap(rcpp_all));
    roll::vec_attr(result, xx);
    
    return result;
    
  }
  
}

// [[Rcpp::export(.roll_sum)]]
SEXP roll_sum(const SEXP& x, const int& width,
              const arma::vec& weights, const int& min_obs,
              const bool& complete_obs, const bool& na_restore,
              const bool& online) {
  
  if (Rf_isMatrix(x)) {
    
    NumericMatrix xx(x);
    int n = weights.size();
    int n_rows_x = xx.nrow();
    int n_cols_x = xx.ncol();
    arma::uvec arma_any_na(n_rows_x);
    arma::mat arma_sum(n_rows_x, n_cols_x);
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_x, width, weights, "'x'");
    bool status = roll::check_lambda(weights, n_rows_x, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // default 'complete_obs' argument is 'false',
    // otherwise check argument for errors
    if (complete_obs) {
      arma_any_na = any_na_x(xx);
    } else {
      arma_any_na.fill(0);
    }
    
    // compute rolling sums
    if (status && online) {
      
      roll::RollSumOnlineMat roll_sum_online(xx, n, n_rows_x, n_cols_x, width,
                                             weights, min_obs,
                                             arma_any_na, na_restore,
                                             arma_sum);
      parallelFor(0, n_cols_x, roll_sum_online);
      
    } else {
      
      roll::RollSumOfflineMat roll_sum_offline(xx, n, n_rows_x, n_cols_x, width,
                                               weights, min_obs,
                                               arma_any_na, na_restore,
                                               arma_sum);
      parallelFor(0, n_rows_x * n_cols_x, roll_sum_offline);
      
    }
    
    // create and return a matrix or xts object
    NumericMatrix result(wrap(arma_sum));
    roll::xts_attr(result, xx);
    
    return result;
    
  } else {
    
    NumericVector xx(x);
    int n = weights.size();
    int n_rows_x = xx.size();
    arma::vec arma_sum(n_rows_x);
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_x, width, weights, "'x'");
    bool status = roll::check_lambda(weights, n_rows_x, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // compute rolling sums
    if (status && online) {
      
      roll::RollSumOnlineVec roll_sum_online(xx, n, n_rows_x, width,
                                             weights, min_obs,
                                             na_restore,
                                             arma_sum);
      roll_sum_online();
      
    } else {
      
      roll::RollSumOfflineVec roll_sum_offline(xx, n, n_rows_x, width,
                                               weights, min_obs,
                                               na_restore,
                                               arma_sum);
      parallelFor(0, n_rows_x, roll_sum_offline);
      
    }
    
    // create and return a vector object
    NumericVector result(wrap(arma_sum));
    roll::vec_attr(result, xx);
    
    return result;
    
  }
  
}

// [[Rcpp::export(.roll_prod)]]
SEXP roll_prod(const SEXP& x, const int& width,
               const arma::vec& weights, const int& min_obs,
               const bool& complete_obs, const bool& na_restore,
               const bool& online) {
  
  if (Rf_isMatrix(x)) {
    
    NumericMatrix xx(x);
    int n = weights.size();
    int n_rows_x = xx.nrow();
    int n_cols_x = xx.ncol();
    arma::uvec arma_any_na(n_rows_x);
    arma::mat arma_prod(n_rows_x, n_cols_x);
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_x, width, weights, "'x'");
    bool status = roll::check_lambda(weights, n_rows_x, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // default 'complete_obs' argument is 'false',
    // otherwise check argument for errors
    if (complete_obs) {
      arma_any_na = any_na_x(xx);
    } else {
      arma_any_na.fill(0);
    }
    
    // compute rolling products
    if (status && online) {
      
      roll::RollProdOnlineMat roll_prod_online(xx, n, n_rows_x, n_cols_x, width,
                                               weights, min_obs,
                                               arma_any_na, na_restore, 
                                               arma_prod);
      parallelFor(0, n_cols_x, roll_prod_online);
      
    } else {
      
      roll::RollProdOfflineMat roll_prod_offline(xx, n, n_rows_x, n_cols_x, width,
                                                 weights, min_obs,
                                                 arma_any_na, na_restore,
                                                 arma_prod);
      parallelFor(0, n_rows_x * n_cols_x, roll_prod_offline);
      
    }
    
    // create and return a matrix or xts object
    NumericMatrix result(wrap(arma_prod));
    roll::xts_attr(result, xx);
    
    return result;
    
  } else {
    
    NumericVector xx(x);
    int n = weights.size();
    int n_rows_x = xx.size();
    arma::vec arma_prod(n_rows_x);
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_x, width, weights, "'x'");
    bool status = roll::check_lambda(weights, n_rows_x, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // compute rolling products
    if (status && online) {
      
      roll::RollProdOnlineVec roll_prod_online(xx, n, n_rows_x, width,
                                               weights, min_obs,
                                               na_restore, 
                                               arma_prod);
      roll_prod_online();
      
    } else {
      
      roll::RollProdOfflineVec roll_prod_offline(xx, n, n_rows_x, width,
                                                 weights, min_obs,
                                                 na_restore,
                                                 arma_prod);
      parallelFor(0, n_rows_x, roll_prod_offline);
      
    }
    
    // create and return a vector object
    NumericVector result(wrap(arma_prod));
    roll::vec_attr(result, xx);
    
    return result;
    
  }
  
}

// [[Rcpp::export(.roll_mean)]]
SEXP roll_mean(const SEXP& x, const int& width,
               const arma::vec& weights, const int& min_obs,
               const bool& complete_obs, const bool& na_restore,
               const bool& online) {
  
  if (Rf_isMatrix(x)) {
    
    NumericMatrix xx(x);
    int n = weights.size();
    int n_rows_x = xx.nrow();
    int n_cols_x = xx.ncol();
    arma::uvec arma_any_na(n_rows_x);
    arma::mat arma_mean(n_rows_x, n_cols_x);
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_x, width, weights, "'x'");
    bool status = roll::check_lambda(weights, n_rows_x, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // default 'complete_obs' argument is 'false',
    // otherwise check argument for errors
    if (complete_obs) {
      arma_any_na = any_na_x(xx);
    } else {
      arma_any_na.fill(0);
    }
    
    // compute rolling means
    if (status && online) {
      
      roll::RollMeanOnlineMat roll_mean_online(xx, n, n_rows_x, n_cols_x, width,
                                               weights, min_obs,
                                               arma_any_na, na_restore,
                                               arma_mean);
      parallelFor(0, n_cols_x, roll_mean_online);
      
    } else {
      
      roll::RollMeanOfflineMat roll_mean_offline(xx, n, n_rows_x, n_cols_x, width,
                                                 weights, min_obs,
                                                 arma_any_na, na_restore,
                                                 arma_mean);
      parallelFor(0, n_rows_x * n_cols_x, roll_mean_offline);
      
    }
    
    // create and return a matrix or xts object
    NumericMatrix result(wrap(arma_mean));
    roll::xts_attr(result, xx);
    
    return result;
    
  } else {
    
    NumericVector xx(x);
    int n = weights.size();
    int n_rows_x = xx.size();
    arma::vec arma_mean(n_rows_x);
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_x, width, weights, "'x'");
    bool status = roll::check_lambda(weights, n_rows_x, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // compute rolling means
    if (status && online) {
      
      roll::RollMeanOnlineVec roll_mean_online(xx, n, n_rows_x, width,
                                               weights, min_obs,
                                               na_restore,
                                               arma_mean);
      roll_mean_online();
      
    } else {
      
      roll::RollMeanOfflineVec roll_mean_offline(xx, n, n_rows_x, width,
                                                 weights, min_obs,
                                                 na_restore,
                                                 arma_mean);
      parallelFor(0, n_rows_x, roll_mean_offline);
      
    }
    
    // create and return a vector object
    NumericVector result(wrap(arma_mean));
    roll::vec_attr(result, xx);
    
    return result;
    
  }
  
}

// [[Rcpp::export(.roll_idxquantile)]]
SEXP roll_idxquantile(const SEXP& x, const int& width,
                      const arma::vec& weights, const double& p,
                      const int& min_obs, const bool& complete_obs,
                      const bool& na_restore, const bool& online) {
  
  if (Rf_isMatrix(x)) {
    
    NumericMatrix xx(x);
    int n = weights.size();
    int n_rows_x = xx.nrow();
    int n_cols_x = xx.ncol();
    arma::uvec arma_any_na(n_rows_x);
    arma::imat arma_idxquantile(n_rows_x, n_cols_x); // unsigned int coerce NA to int
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_x, width, weights, "'x'");
    
    // check 'p' argument for errors
    roll::check_bounds_double(p, 0, 1, "p");
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // default 'complete_obs' argument is 'false',
    // otherwise check argument for errors
    if (complete_obs) {
      arma_any_na = any_na_x(xx);
    } else {
      arma_any_na.fill(0);
    }
    
    // compute rolling index of quantiles
    if (online) {
      
      if (p == 0) {
        
        roll::RollIdxMinOnlineMat roll_idxmin_online(xx, n, n_rows_x, n_cols_x, width,
                                                     weights, min_obs,
                                                     arma_any_na, na_restore,
                                                     arma_idxquantile);
        parallelFor(0, n_cols_x, roll_idxmin_online);
        
      } else if (p == 1) {
        
        roll::RollIdxMaxOnlineMat roll_idxmax_online(xx, n, n_rows_x, n_cols_x, width,
                                                     weights, min_obs,
                                                     arma_any_na, na_restore,
                                                     arma_idxquantile);
        parallelFor(0, n_cols_x, roll_idxmax_online);
        
      }
      
    } else {
      
      if (p == 0) {
        
        roll::RollIdxMinOfflineMat roll_idxmin_offline(xx, n, n_rows_x, n_cols_x, width,
                                                       weights, min_obs,
                                                       arma_any_na, na_restore,
                                                       arma_idxquantile);
        parallelFor(0, n_rows_x * n_cols_x, roll_idxmin_offline);
        
      } else if (p == 1) {
        
        roll::RollIdxMaxOfflineMat roll_idxmax_offline(xx, n, n_rows_x, n_cols_x, width,
                                                       weights, min_obs,
                                                       arma_any_na, na_restore,
                                                       arma_idxquantile);
        parallelFor(0, n_rows_x * n_cols_x, roll_idxmax_offline);
        
      }
      
    }
    
    // create and return a matrix or xts object
    IntegerMatrix result(wrap(arma_idxquantile));
    roll::xts_attr(result, xx);
    
    return result;
    
  } else {
    
    NumericVector xx(x);
    int n = weights.size();
    int n_rows_x = xx.size();
    arma::ivec arma_idxquantile(n_rows_x); // unsigned int coerce NA to int
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_x, width, weights, "'x'");
    
    // check 'p' argument for errors
    roll::check_bounds_double(p, 0, 1, "p");
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // compute rolling index of quantiles
    if (online) {
      
      if (p == 0) {
        
        roll::RollIdxMinOnlineVec roll_idxmin_online(xx, n, n_rows_x, width,
                                                     weights, min_obs,
                                                     na_restore,
                                                     arma_idxquantile);
        roll_idxmin_online();
        
      } else if (p == 1) {
        
        roll::RollIdxMaxOnlineVec roll_idxmax_online(xx, n, n_rows_x, width,
                                                     weights, min_obs,
                                                     na_restore,
                                                     arma_idxquantile);
        roll_idxmax_online();
        
      }
      
    } else {
      
      if (p == 0) {
        
        roll::RollIdxMinOfflineVec roll_idxmin_offline(xx, n, n_rows_x, width,
                                                       weights, min_obs,
                                                       na_restore,
                                                       arma_idxquantile);
        parallelFor(0, n_rows_x, roll_idxmin_offline);
        
      } else if (p == 1) {
        
        roll::RollIdxMaxOfflineVec roll_idxmax_offline(xx, n, n_rows_x, width,
                                                       weights, min_obs,
                                                       na_restore,
                                                       arma_idxquantile);
        parallelFor(0, n_rows_x, roll_idxmax_offline);
        
      }
      
    }
    
    // create and return a vector object
    IntegerVector result(wrap(arma_idxquantile));
    roll::vec_attr(result, xx);
    
    return result;
    
  }
  
}

// [[Rcpp::export(.roll_quantile)]]
SEXP roll_quantile(const SEXP& x, const int& width,
                   const arma::vec& weights, const double& p,
                   const int& min_obs, const bool& complete_obs,
                   const bool& na_restore, const bool& online) {
  
  if (Rf_isMatrix(x)) {
    
    NumericMatrix xx(x);
    int n = weights.size();
    int n_rows_x = xx.nrow();
    int n_cols_x = xx.ncol();
    arma::uvec arma_any_na(n_rows_x);
    arma::mat arma_quantile(n_rows_x, n_cols_x);
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_x, width, weights, "'x'");
    
    // check 'p' argument for errors
    roll::check_bounds_double(p, 0, 1, "p");
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // default 'complete_obs' argument is 'false',
    // otherwise check argument for errors
    if (complete_obs) {
      arma_any_na = any_na_x(xx);
    } else {
      arma_any_na.fill(0);
    }
    
    // compute rolling quantiles
    if (online) {
      
      if (p == 0) {
        
        roll::RollMinOnlineMat roll_min_online(xx, n, n_rows_x, n_cols_x, width,
                                               weights, min_obs,
                                               arma_any_na, na_restore,
                                               arma_quantile);
        parallelFor(0, n_cols_x, roll_min_online);
        
      } else if (p == 1) {
        
        roll::RollMaxOnlineMat roll_max_online(xx, n, n_rows_x, n_cols_x, width,
                                               weights, min_obs,
                                               arma_any_na, na_restore,
                                               arma_quantile);
        parallelFor(0, n_cols_x, roll_max_online);
        
      } else {
        
        if (any(weights != weights[0])) {
          stop("'online' is only supported for equal 'weights'");
        }
        
        roll::RollQuantileOnlineMat roll_quantile_online(xx, n, n_rows_x, n_cols_x, width,
                                                         weights, p, min_obs,
                                                         arma_any_na, na_restore,
                                                         arma_quantile);
        parallelFor(0, n_cols_x, roll_quantile_online);
        
      }
      
    } else {
      
      if (p == 0) {
        
        roll::RollMinOfflineMat roll_min_offline(xx, n, n_rows_x, n_cols_x, width,
                                                 weights, min_obs,
                                                 arma_any_na, na_restore,
                                                 arma_quantile);
        parallelFor(0, n_rows_x * n_cols_x, roll_min_offline);
        
      } else if (p == 1) {
        
        roll::RollMaxOfflineMat roll_max_offline(xx, n, n_rows_x, n_cols_x, width,
                                                 weights, min_obs,
                                                 arma_any_na, na_restore,
                                                 arma_quantile);
        parallelFor(0, n_rows_x * n_cols_x, roll_max_offline);
        
      } else {
        
        roll::RollQuantileOfflineMat roll_quantile_offline(xx, n, n_rows_x, n_cols_x, width,
                                                           weights, 1 - p, min_obs,
                                                           arma_any_na, na_restore,
                                                           arma_quantile);
        parallelFor(0, n_rows_x * n_cols_x, roll_quantile_offline);
        
      }
      
    }
    
    // create and return a matrix or xts object
    NumericMatrix result(wrap(arma_quantile));
    roll::xts_attr(result, xx);
    
    return result;
    
  } else {
    
    NumericVector xx(x);
    int n = weights.size();
    int n_rows_x = xx.size();
    arma::vec arma_quantile(n_rows_x);
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_x, width, weights, "'x'");
    
    // check 'p' argument for errors
    roll::check_bounds_double(p, 0, 1, "p");
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // compute rolling quantiles
    if (online) {
      
      if (p == 0) {
        
        roll::RollMinOnlineVec roll_min_online(xx, n, n_rows_x, width,
                                               weights, min_obs,
                                               na_restore,
                                               arma_quantile);
        roll_min_online();
        
      } else if (p == 1) {
        
        roll::RollMaxOnlineVec roll_max_online(xx, n, n_rows_x, width,
                                               weights, min_obs,
                                               na_restore,
                                               arma_quantile);
        roll_max_online();
        
      } else {
        
        if (any(weights != weights[0])) {
          stop("'online' is only supported for equal 'weights'");
        }
        
        roll::RollQuantileOnlineVec roll_quantile_online(xx, n, n_rows_x, width,
                                                         weights, p, min_obs,
                                                         na_restore,
                                                         arma_quantile);
        roll_quantile_online();
        
      }
      
    } else {
      
      if (p == 0) {
        
        roll::RollMinOfflineVec roll_min_offline(xx, n, n_rows_x, width,
                                                 weights, min_obs,
                                                 na_restore,
                                                 arma_quantile);
        parallelFor(0, n_rows_x, roll_min_offline);
        
      } else if (p == 1) {
        
        roll::RollMaxOfflineVec roll_max_offline(xx, n, n_rows_x, width,
                                                 weights, min_obs,
                                                 na_restore,
                                                 arma_quantile);
        parallelFor(0, n_rows_x, roll_max_offline);
        
      } else {
        
        roll::RollQuantileOfflineVec roll_quantile_offline(xx, n, n_rows_x, width,
                                                           weights, 1 - p, min_obs,
                                                           na_restore,
                                                           arma_quantile);
        parallelFor(0, n_rows_x, roll_quantile_offline);
        
      }
      
    }
      
    NumericVector result(wrap(arma_quantile));
    roll::vec_attr(result, xx);
    
    return result;
    
  }
  
}

// [[Rcpp::export(.roll_var)]]
SEXP roll_var(const SEXP& x, const int& width,
              const arma::vec& weights, const bool& center,
              const int& min_obs, const bool& complete_obs,
              const bool& na_restore, const bool& online) {
  
  if (Rf_isMatrix(x)) {
    
    NumericMatrix xx(x);
    int n = weights.size();
    int n_rows_x = xx.nrow();
    int n_cols_x = xx.ncol();
    arma::uvec arma_any_na(n_rows_x);
    arma::mat arma_var(n_rows_x, n_cols_x);
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_x, width, weights, "'x'");
    bool status = roll::check_lambda(weights, n_rows_x, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // default 'complete_obs' argument is 'false',
    // otherwise check argument for errors
    if (complete_obs) {
      arma_any_na = any_na_x(x);
    } else {
      arma_any_na.fill(0);
    }
    
    // compute rolling variances
    if (status && online) {
      
      roll::RollVarOnlineMat roll_var_online(xx, n, n_rows_x, n_cols_x, width,
                                             weights, center, min_obs,
                                             arma_any_na, na_restore,
                                             arma_var);
      parallelFor(0, n_cols_x, roll_var_online);
      
    } else {
      
      roll::RollVarOfflineMat roll_var_offline(xx, n, n_rows_x, n_cols_x, width,
                                               weights, center, min_obs,
                                               arma_any_na, na_restore,
                                               arma_var);
      parallelFor(0, n_rows_x * n_cols_x, roll_var_offline);
      
    }
    
    // create and return a matrix or xts object
    NumericMatrix result(wrap(arma_var));
    roll::xts_attr(result, xx);
    
    return result;
    
  } else {
    
    NumericVector xx(x);
    int n = weights.size();
    int n_rows_x = xx.size();
    arma::vec arma_var(n_rows_x);
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_x, width, weights, "'x'");
    bool status = roll::check_lambda(weights, n_rows_x, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // compute rolling variances
    if (status && online) {
      
      roll::RollVarOnlineVec roll_var_online(xx, n, n_rows_x, width,
                                             weights, center, min_obs,
                                             na_restore,
                                             arma_var);
      roll_var_online();
      
    } else {
      
      roll::RollVarOfflineVec roll_var_offline(xx, n, n_rows_x, width,
                                               weights, center, min_obs,
                                               na_restore,
                                               arma_var);
      parallelFor(0, n_rows_x, roll_var_offline);
      
    }
    
    // create and return a vector object
    NumericVector result(wrap(arma_var));
    roll::vec_attr(result, xx);
    
    return result;
    
  }
  
}

// [[Rcpp::export(.roll_sd)]]
SEXP roll_sd(const SEXP& x, const int& width,
             const arma::vec& weights, const bool& center,
             const int& min_obs, const bool& complete_obs,
             const bool& na_restore, const bool& online) {
  
  if (Rf_isMatrix(x)) {
    
    NumericMatrix xx(x);
    int n = weights.size();
    int n_rows_x = xx.nrow();
    int n_cols_x = xx.ncol();
    arma::uvec arma_any_na(n_rows_x);
    arma::mat arma_sd(n_rows_x, n_cols_x);
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_x, width, weights, "'x'");
    bool status = roll::check_lambda(weights, n_rows_x, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // default 'complete_obs' argument is 'false',
    // otherwise check argument for errors
    if (complete_obs) {
      arma_any_na = any_na_x(xx);
    } else {
      arma_any_na.fill(0);
    }
    
    // compute rolling standard deviations
    if (status && online) {
      
      roll::RollSdOnlineMat roll_sd_online(xx, n, n_rows_x, n_cols_x, width,
                                           weights, center, min_obs,
                                           arma_any_na, na_restore,
                                           arma_sd);
      parallelFor(0, n_cols_x, roll_sd_online);
      
    } else {
      
      roll::RollSdOfflineMat roll_sd_offline(xx, n, n_rows_x, n_cols_x, width,
                                             weights, center, min_obs,
                                             arma_any_na, na_restore,
                                             arma_sd);
      parallelFor(0, n_rows_x * n_cols_x, roll_sd_offline);
      
    }
    
    // create and return a matrix or xts object
    NumericMatrix result(wrap(arma_sd));
    roll::xts_attr(result, xx);
    
    return result;
    
  } else {
    
    NumericVector xx(x);
    int n = weights.size();
    int n_rows_x = xx.size();
    arma::vec arma_sd(n_rows_x);
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_x, width, weights, "'x'");
    bool status = roll::check_lambda(weights, n_rows_x, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // compute rolling standard deviations
    if (status && online) {
      
      roll::RollSdOnlineVec roll_sd_online(xx, n, n_rows_x, width,
                                           weights, center, min_obs,
                                           na_restore,
                                           arma_sd);
      roll_sd_online();
      
    } else {
      
      roll::RollSdOfflineVec roll_sd_offline(xx, n, n_rows_x, width,
                                             weights, center, min_obs,
                                             na_restore,
                                             arma_sd);
      parallelFor(0, n_rows_x, roll_sd_offline);
      
    }
    
    // create and return a vector object
    NumericVector result(wrap(arma_sd));
    roll::vec_attr(result, xx);
    
    return result;
    
  }
  
}

// [[Rcpp::export(.roll_scale)]]
SEXP roll_scale(const SEXP& x, const int& width,
                const arma::vec& weights, const bool& center,
                const bool& scale, const int& min_obs,
                const bool& complete_obs, const bool& na_restore,
                const bool& online) {
  
  if (Rf_isMatrix(x)) {
    
    NumericMatrix xx(x);
    int n = weights.size();
    int n_rows_x = xx.nrow();
    int n_cols_x = xx.ncol();
    arma::uvec arma_any_na(n_rows_x);
    arma::mat arma_scale(n_rows_x, n_cols_x);
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_x, width, weights, "'x'");
    bool status = roll::check_lambda(weights, n_rows_x, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // default 'complete_obs' argument is 'false',
    // otherwise check argument for errors
    if (complete_obs) {
      arma_any_na = any_na_x(xx);
    } else {
      arma_any_na.fill(0);
    }
    
    // compute rolling centering and scaling
    if (status && online) {
      
      roll::RollScaleOnlineMat roll_scale_online(xx, n, n_rows_x, n_cols_x, width,
                                                 weights, center, scale, min_obs,
                                                 arma_any_na, na_restore,
                                                 arma_scale);
      parallelFor(0, n_cols_x, roll_scale_online);
      
    } else {
      
      roll::RollScaleOfflineMat roll_scale_offline(xx, n, n_rows_x, n_cols_x, width,
                                                   weights, center, scale, min_obs,
                                                   arma_any_na, na_restore,
                                                   arma_scale);
      parallelFor(0, n_rows_x * n_cols_x, roll_scale_offline);
      
    }
    
    // create and return a matrix or xts object
    NumericMatrix result(wrap(arma_scale));
    roll::xts_attr(result, xx);
    
    return result;
    
  } else {
    
    NumericVector xx(x);
    int n = weights.size();
    int n_rows_x = xx.size();
    arma::vec arma_scale(n_rows_x);
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_x, width, weights, "'x'");
    bool status = roll::check_lambda(weights, n_rows_x, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // compute rolling centering and scaling
    if (status && online) {
      
      roll::RollScaleOnlineVec roll_scale_online(xx, n, n_rows_x, width,
                                                 weights, center, scale, min_obs,
                                                 na_restore,
                                                 arma_scale);
      roll_scale_online();
      
    } else {
      
      roll::RollScaleOfflineVec roll_scale_offline(xx, n, n_rows_x, width,
                                                   weights, center, scale, min_obs,
                                                   na_restore,
                                                   arma_scale);
      parallelFor(0, n_rows_x, roll_scale_offline);
      
    }
    
    // create and return a vector object
    NumericVector result(wrap(arma_scale));
    roll::vec_attr(result, xx);
    
    return result;
    
  }
  
}

SEXP roll_cov_z(const SEXP& x, const SEXP& y,
                const int& width, const arma::vec& weights,
                const bool& center, const bool& scale,
                const int& min_obs, const bool& complete_obs,
                const bool& na_restore, const bool& online,
                const bool& symmetric) {
  
  if (Rf_isMatrix(x) && Rf_isMatrix(y)) {
    
    NumericMatrix xx(x);
    NumericMatrix yy(y);
    int n = weights.size();
    int n_rows_xy = xx.nrow();
    int n_cols_x = xx.ncol();
    int n_cols_y = yy.ncol();
    arma::uvec arma_any_na(n_rows_xy);
    arma::cube arma_cov(n_cols_x, n_cols_y, n_rows_xy);
    
    // check 'x' and 'y' arguments for errors
    roll::check_rows_equal(n_rows_xy, yy.nrow(), "x", "y");
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_xy, width, weights, "'x' (and 'y', if applicable)");
    bool status = roll::check_lambda(weights, n_rows_xy, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // default 'complete_obs' argument is 'true',
    // otherwise check argument for errors
    if (complete_obs && symmetric) {
      arma_any_na = any_na_x(xx);
    } else if (complete_obs && !symmetric) {
      arma_any_na = any_na_xy(xx, yy);
    } else {
      arma_any_na.fill(0);
    }
    
    // compute rolling covariances
    if (status && online) {
      
      if (symmetric) {
        
        // y is null
        roll::RollCovOnlineMatXX roll_cov_online(xx, n, n_rows_xy, n_cols_x, width,
                                                 weights, center, scale, min_obs,
                                                 arma_any_na, na_restore,
                                                 arma_cov);
        parallelFor(0, n_cols_x, roll_cov_online);
        
      } else if (!symmetric) {
        
        // y is not null
        roll::RollCovOnlineMatXY roll_cov_online(xx, yy, n, n_rows_xy, n_cols_x, n_cols_y, width,
                                                 weights, center, scale, min_obs,
                                                 arma_any_na, na_restore,
                                                 arma_cov);
        parallelFor(0, n_cols_x, roll_cov_online);
        
      }
      
    } else {
      
      if (symmetric) {
        
        // y is null
        roll::RollCovOfflineMatXX roll_cov_offline(xx, n, n_rows_xy, n_cols_x, width,
                                                   weights, center, scale, min_obs,
                                                   arma_any_na, na_restore,
                                                   arma_cov);
        parallelFor(0, n_rows_xy * n_cols_x * (n_cols_x + 1) / 2, roll_cov_offline);
        
      } else if (!symmetric) {
        
        // y is not null
        roll::RollCovOfflineMatXY roll_cov_offline(xx, yy, n, n_rows_xy, n_cols_x, n_cols_y, width,
                                                   weights, center, scale, min_obs,
                                                   arma_any_na, na_restore,
                                                   arma_cov);
        parallelFor(0, n_rows_xy * n_cols_x * n_cols_y, roll_cov_offline);
        
      }
      
    }
    
    // create and return a matrix
    NumericVector result(wrap(arma_cov));
    List dimnames_x = xx.attr("dimnames");
    List dimnames_y = yy.attr("dimnames");
    roll::cube_attr(result, n_cols_x, n_cols_y, n_rows_xy, dimnames_x, dimnames_y);
    
    return result;
    
  } else if (Rf_isMatrix(x)) {
    
    NumericMatrix xx(x);
    NumericVector yy(y);
    // yy.attr("dim") = IntegerVector::create(yy.size(), 1);
    // NumericMatrix yyy(wrap(yy));
    NumericMatrix yyy(yy.size(), 1, yy.begin());
    
    int n = weights.size();
    int n_rows_xy = xx.nrow();
    int n_cols_x = xx.ncol();
    int n_cols_y = yyy.ncol();
    arma::uvec arma_any_na(n_rows_xy);
    arma::cube arma_cov(n_cols_x, n_cols_y, n_rows_xy);
    
    // check 'x' and 'y' arguments for errors
    roll::check_rows_equal(n_rows_xy, yyy.nrow(), "x", "y");
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_xy, width, weights, "'x' (and 'y', if applicable)");
    bool status = roll::check_lambda(weights, n_rows_xy, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // default 'complete_obs' argument is 'true',
    // otherwise check argument for errors
    if (complete_obs) {
      arma_any_na = any_na_xy(xx, yyy);
    } else {
      arma_any_na.fill(0);
    }
    
    // compute rolling covariances
    if (status && online) {
      
      roll::RollCovOnlineMatXY roll_cov_online(xx, yyy, n, n_rows_xy, n_cols_x, n_cols_y, width,
                                               weights, center, scale, min_obs,
                                               arma_any_na, na_restore,
                                               arma_cov);
      parallelFor(0, n_cols_x, roll_cov_online);
      
    } else {
      
      roll::RollCovOfflineMatXY roll_cov_offline(xx, yyy, n, n_rows_xy, n_cols_x, n_cols_y, width,
                                                 weights, center, scale, min_obs,
                                                 arma_any_na, na_restore,
                                                 arma_cov);
      parallelFor(0, n_rows_xy * n_cols_x * n_cols_y, roll_cov_offline);
      
    }
    
    // create and return a matrix
    NumericVector result(wrap(arma_cov));
    List dimnames_x = xx.attr("dimnames");
    List dimnames_y = yyy.attr("dimnames");
    roll::cube_attr(result, n_cols_x, n_cols_y, n_rows_xy, dimnames_x, dimnames_y);
    
    return result;
    
  } else if (Rf_isMatrix(y)) {
    
    NumericVector xx(x);
    NumericMatrix yy(y);
    // xx.attr("dim") = IntegerVector::create(xx.size(), 1);
    // NumericMatrix xxx(wrap(xx));
    NumericMatrix xxx(xx.size(), 1, xx.begin());
    
    int n = weights.size();
    int n_rows_xy = xxx.nrow();
    int n_cols_x = xxx.ncol();
    int n_cols_y = yy.ncol();
    arma::uvec arma_any_na(n_rows_xy);
    arma::cube arma_cov(n_cols_x, n_cols_y, n_rows_xy);
    
    // check 'x' and 'y' arguments for errors
    roll::check_rows_equal(n_rows_xy, yy.nrow(), "x", "y");
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_xy, width, weights, "'x' (and 'y', if applicable)");
    bool status = roll::check_lambda(weights, n_rows_xy, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // default 'complete_obs' argument is 'true',
    // otherwise check argument for errors
    if (complete_obs) {
      arma_any_na = any_na_xy(xxx, yy);
    } else {
      arma_any_na.fill(0);
    }
    
    // compute rolling covariances
    if (status && online) {
      
      roll::RollCovOnlineMatXY roll_cov_online(xxx, yy, n, n_rows_xy, n_cols_x, n_cols_y, width,
                                               weights, center, scale, min_obs,
                                               arma_any_na, na_restore,
                                               arma_cov);
      parallelFor(0, n_cols_x, roll_cov_online);
      
    } else {
      
      // y is not null
      roll::RollCovOfflineMatXY roll_cov_offline(xxx, yy, n, n_rows_xy, n_cols_x, n_cols_y, width,
                                                 weights, center, scale, min_obs,
                                                 arma_any_na, na_restore,
                                                 arma_cov);
      parallelFor(0, n_rows_xy * n_cols_x * n_cols_y, roll_cov_offline);
      
      
    }
    
    // create and return a matrix
    NumericVector result(wrap(arma_cov));
    List dimnames_x = xxx.attr("dimnames");
    List dimnames_y = yy.attr("dimnames");
    roll::cube_attr(result, n_cols_x, n_cols_y, n_rows_xy, dimnames_x, dimnames_y);
    
    return result;
    
  } else {
    
    NumericVector xx(x);
    NumericVector yy(y);
    int n = weights.size();
    int n_rows_xy = xx.size();
    arma::vec arma_cov(n_rows_xy);
    
    // check 'x' and 'y' arguments for errors
    roll::check_rows_equal(n_rows_xy, yy.size(), "x", "y");
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_xy, width, weights, "'x' (and 'y', if applicable)");
    bool status = roll::check_lambda(weights, n_rows_xy, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // compute rolling covariances
    if (status && online) {
      
      if (symmetric) {
        
        // y is null
        roll::RollCovOnlineVecXX roll_cov_online(xx, n, n_rows_xy, width,
                                                 weights, center, scale, min_obs,
                                                 na_restore,
                                                 arma_cov);
        roll_cov_online();
        
      } else if (!symmetric) {
        
        // y is not null
        roll::RollCovOnlineVecXY roll_cov_online(xx, yy, n, n_rows_xy, width,
                                                 weights, center, scale, min_obs,
                                                 na_restore,
                                                 arma_cov);
        roll_cov_online();
        
      }
      
    } else {
      
      if (symmetric) {
        
        // y is null
        roll::RollCovOfflineVecXX roll_cov_offline(xx, n, n_rows_xy, width,
                                                   weights, center, scale, min_obs,
                                                   na_restore,
                                                   arma_cov);
        parallelFor(0, n_rows_xy, roll_cov_offline);
        
      } else if (!symmetric) {
        
        // y is not null
        roll::RollCovOfflineVecXY roll_cov_offline(xx, yy, n, n_rows_xy, width,
                                                   weights, center, scale, min_obs,
                                                   na_restore,
                                                   arma_cov);
        parallelFor(0, n_rows_xy, roll_cov_offline);
        
      }
      
    }
    
    // create and return a vector object
    NumericVector result(wrap(arma_cov));
    roll::strip_dim(result);
    
    return result;
    
  }
  
}

// [[Rcpp::export(.roll_cov)]]
SEXP roll_cov(const SEXP& x, const SEXP& y,
              const int& width, const arma::vec& weights,
              const bool& center, const bool& scale,
              const int& min_obs, const bool& complete_obs,
              const bool& na_restore, const bool& online) {
  
  if (Rf_isNull(y)) {
    
    return roll_cov_z(x, x, width, weights, center, scale, min_obs, complete_obs, 
                      na_restore, online, true);
    
  } else {
    
    return roll_cov_z(x, y, width, weights, center, scale, min_obs, complete_obs, 
                      na_restore, online, false);
    
  }
  
}

SEXP roll_crossprod_z(const SEXP& x, const SEXP& y,
                      const int& width, const arma::vec& weights,
                      const bool& center, const bool& scale,
                      const int& min_obs, const bool& complete_obs,
                      const bool& na_restore, const bool& online,
                      const bool& symmetric) {
  
  if (Rf_isMatrix(x) && Rf_isMatrix(y)) {
    
    NumericMatrix xx(x);
    NumericMatrix yy(y);
    int n = weights.size();
    int n_rows_xy = xx.nrow();
    int n_cols_x = xx.ncol();
    int n_cols_y = yy.ncol();
    arma::uvec arma_any_na(n_rows_xy);
    arma::vec arma_n_obs(n_rows_xy);
    arma::vec arma_sum_w(n_rows_xy);
    arma::mat arma_mean(n_rows_xy, n_cols_x);
    arma::cube arma_crossprod(n_cols_x, n_cols_y, n_rows_xy);
    
    // check 'x' and 'y' arguments for errors
    roll::check_rows_equal(n_rows_xy, yy.nrow(), "x", "y");
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_xy, width, weights, "'x' (and 'y', if applicable)");
    bool status = roll::check_lambda(weights, n_rows_xy, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // default 'complete_obs' argument is 'true',
    // otherwise check argument for errors
    if (complete_obs && symmetric) {
      arma_any_na = any_na_x(xx);
    } else if (complete_obs && !symmetric) {
      arma_any_na = any_na_xy(xx, yy);
    } else {
      arma_any_na.fill(0);
    }
    
    // compute rolling crossproducts
    if (status && online) {
      
      if (symmetric) {
        
        // y is null
        roll::RollCrossProdOnlineMatXX roll_crossprod_online(xx, n, n_rows_xy, n_cols_x, width,
                                                             weights, center, scale, min_obs,
                                                             arma_any_na, na_restore,
                                                             arma_n_obs, arma_sum_w, arma_mean,
                                                             arma_crossprod);
        parallelFor(0, n_cols_x, roll_crossprod_online);
        
      } else if (!symmetric) {
        
        // y is not null
        roll::RollCrossProdOnlineMatXY roll_crossprod_online(xx, yy, n, n_rows_xy, n_cols_x, n_cols_y, width,
                                                             weights, center, scale, min_obs,
                                                             arma_any_na, na_restore,
                                                             arma_crossprod);
        parallelFor(0, n_cols_x, roll_crossprod_online);
        
      }
      
    } else {
      
      if (symmetric) {
        
        // y is null
        roll::RollCrossProdOfflineMatXX roll_crossprod_offline(xx, n, n_rows_xy, n_cols_x, width,
                                                               weights, center, scale, min_obs,
                                                               arma_any_na, na_restore,
                                                               arma_n_obs, arma_sum_w, arma_mean,
                                                               arma_crossprod);
        parallelFor(0, n_rows_xy * n_cols_x * (n_cols_x + 1) / 2, roll_crossprod_offline);
        
      } else if (!symmetric) {
        
        // y is not null
        roll::RollCrossProdOfflineMatXY roll_crossprod_offline(xx, yy, n, n_rows_xy, n_cols_x, n_cols_y, width,
                                                               weights, center, scale, min_obs,
                                                               arma_any_na, na_restore,
                                                               arma_crossprod);
        parallelFor(0, n_rows_xy * n_cols_x * n_cols_y, roll_crossprod_offline);
        
      }
      
    }
    
    // create and return a matrix
    NumericVector result(wrap(arma_crossprod));
    List dimnames_x = xx.attr("dimnames");
    List dimnames_y = yy.attr("dimnames");
    roll::cube_attr(result, n_cols_x, n_cols_y, n_rows_xy, dimnames_x, dimnames_y);
    
    return result;
    
  } else if (Rf_isMatrix(x)) {
    
    NumericMatrix xx(x);
    NumericVector yy(y);
    // yy.attr("dim") = IntegerVector::create(yy.size(), 1);
    // NumericMatrix yyy(wrap(yy));
    NumericMatrix yyy(yy.size(), 1, yy.begin());
    
    int n = weights.size();
    int n_rows_xy = xx.nrow();
    int n_cols_x = xx.ncol();
    int n_cols_y = yyy.ncol();
    arma::uvec arma_any_na(n_rows_xy);
    arma::cube arma_crossprod(n_cols_x, n_cols_y, n_rows_xy);
    
    // check 'x' and 'y' arguments for errors
    roll::check_rows_equal(n_rows_xy, yyy.nrow(), "x", "y");
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_xy, width, weights, "'x' (and 'y', if applicable)");
    bool status = roll::check_lambda(weights, n_rows_xy, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // default 'complete_obs' argument is 'true',
    // otherwise check argument for errors
    if (complete_obs) {
      arma_any_na = any_na_xy(xx, yyy);
    } else {
      arma_any_na.fill(0);
    }
    
    // compute rolling crossproducts
    if (status && online) {
      
      roll::RollCrossProdOnlineMatXY roll_crossprod_online(xx, yyy, n, n_rows_xy, n_cols_x, n_cols_y, width,
                                                           weights, center, scale, min_obs,
                                                           arma_any_na, na_restore,
                                                           arma_crossprod);
      parallelFor(0, n_cols_x, roll_crossprod_online);
      
    } else {
      
      roll::RollCrossProdOfflineMatXY roll_crossprod_offline(xx, yyy, n, n_rows_xy, n_cols_x, n_cols_y, width,
                                                             weights, center, scale, min_obs,
                                                             arma_any_na, na_restore,
                                                             arma_crossprod);
      parallelFor(0, n_rows_xy * n_cols_x * n_cols_y, roll_crossprod_offline);
      
    }
    
    // create and return a matrix
    NumericVector result(wrap(arma_crossprod));
    List dimnames_x = xx.attr("dimnames");
    List dimnames_y = yyy.attr("dimnames");
    roll::cube_attr(result, n_cols_x, n_cols_y, n_rows_xy, dimnames_x, dimnames_y);
    
    return result;
    
  } else if (Rf_isMatrix(y)) {
    
    NumericVector xx(x);
    NumericMatrix yy(y);
    // xx.attr("dim") = IntegerVector::create(xx.size(), 1);
    // NumericMatrix xxx(wrap(xx));
    NumericMatrix xxx(xx.size(), 1, xx.begin());
    
    int n = weights.size();
    int n_rows_xy = xxx.nrow();
    int n_cols_x = xxx.ncol();
    int n_cols_y = yy.ncol();
    arma::uvec arma_any_na(n_rows_xy);
    arma::cube arma_crossprod(n_cols_x, n_cols_y, n_rows_xy);
    
    // check 'x' and 'y' arguments for errors
    roll::check_rows_equal(n_rows_xy, yy.nrow(), "x", "y");
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_xy, width, weights, "'x' (and 'y', if applicable)");
    bool status = roll::check_lambda(weights, n_rows_xy, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // default 'complete_obs' argument is 'true',
    // otherwise check argument for errors
    if (complete_obs) {
      arma_any_na = any_na_xy(xxx, yy);
    } else {
      arma_any_na.fill(0);
    }
    
    // compute rolling crossproducts
    if (status && online) {
      
      roll::RollCrossProdOnlineMatXY roll_crossprod_online(xxx, yy, n, n_rows_xy, n_cols_x, n_cols_y, width,
                                                           weights, center, scale, min_obs,
                                                           arma_any_na, na_restore,
                                                           arma_crossprod);
      parallelFor(0, n_cols_x, roll_crossprod_online);
      
    } else {
      
      // y is not null
      roll::RollCrossProdOfflineMatXY roll_crossprod_offline(xxx, yy, n, n_rows_xy, n_cols_x, n_cols_y, width,
                                                             weights, center, scale, min_obs,
                                                             arma_any_na, na_restore,
                                                             arma_crossprod);
      parallelFor(0, n_rows_xy * n_cols_x * n_cols_y, roll_crossprod_offline);
      
      
    }
    
    // create and return a matrix
    NumericVector result(wrap(arma_crossprod));
    List dimnames_x = xxx.attr("dimnames");
    List dimnames_y = yy.attr("dimnames");
    roll::cube_attr(result, n_cols_x, n_cols_y, n_rows_xy, dimnames_x, dimnames_y);
    
    return result;
    
  } else {
    
    NumericVector xx(x);
    NumericVector yy(y);
    int n = weights.size();
    int n_rows_xy = xx.size();
    arma::vec arma_crossprod(n_rows_xy);
    
    // check 'x' and 'y' arguments for errors
    roll::check_rows_equal(n_rows_xy, yy.size(), "x", "y");
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_xy, width, weights, "'x' (and 'y', if applicable)");
    bool status = roll::check_lambda(weights, n_rows_xy, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // compute rolling crossproducts
    if (status && online) {
      
      if (symmetric) {
        
        // y is null
        roll::RollCrossProdOnlineVecXX roll_crossprod_online(xx, n, n_rows_xy, width,
                                                             weights, center, scale, min_obs,
                                                             na_restore,
                                                             arma_crossprod);
        roll_crossprod_online();
        
      } else if (!symmetric) {
        
        // y is not null
        roll::RollCrossProdOnlineVecXY roll_crossprod_online(xx, yy, n, n_rows_xy, width,
                                                             weights, center, scale, min_obs,
                                                             na_restore,
                                                             arma_crossprod);
        roll_crossprod_online();
        
      }
      
    } else {
      
      if (symmetric) {
        
        // y is null
        roll::RollCrossProdOfflineVecXX roll_crossprod_offline(xx, n, n_rows_xy, width,
                                                               weights, center, scale, min_obs,
                                                               na_restore,
                                                               arma_crossprod);
        parallelFor(0, n_rows_xy, roll_crossprod_offline);
        
      } else if (!symmetric) {
        
        // y is not null
        roll::RollCrossProdOfflineVecXY roll_crossprod_offline(xx, yy, n, n_rows_xy, width,
                                                               weights, center, scale, min_obs,
                                                               na_restore,
                                                               arma_crossprod);
        parallelFor(0, n_rows_xy, roll_crossprod_offline);
        
      }
      
    }
    
    // create and return a vector object
    NumericVector result(wrap(arma_crossprod));
    roll::strip_dim(result);
    
    return result;
    
  }
  
}

// [[Rcpp::export(.roll_crossprod)]]
SEXP roll_crossprod(const SEXP& x, const SEXP& y,
                    const int& width, const arma::vec& weights,
                    const bool& center, const bool& scale,
                    const int& min_obs, const bool& complete_obs,
                    const bool& na_restore, const bool& online) {
  
  if (Rf_isNull(y)) {
    
    return roll_crossprod_z(x, x, width, weights, center, scale, min_obs, complete_obs, 
                            na_restore, online, true);
    
  } else {
    
    return roll_crossprod_z(x, y, width, weights, center, scale, min_obs, complete_obs, 
                            na_restore, online, false);
    
  }
  
}

List roll_lm_z(const SEXP& x, const NumericVector& y,
               const int& width, const arma::vec& weights,
               const bool& intercept, const int& min_obs,
               const bool& complete_obs, const bool& na_restore,
               const bool& online) {
  
  if (Rf_isMatrix(x)) {
    
    NumericMatrix xx(x);
    int n = weights.size();
    int n_rows_xy = xx.nrow();
    int n_cols_x = xx.ncol() + 1;
    arma::uvec arma_any_na(n_rows_xy);
    arma::vec arma_n_obs(n_rows_xy);
    arma::vec arma_sum_w(n_rows_xy);
    arma::mat arma_mean(n_rows_xy, n_cols_x);
    arma::cube arma_cov(n_cols_x, n_cols_x, n_rows_xy);
    arma::vec arma_rsq(n_rows_xy);
    List result(3);
    
    // check 'x' and 'y' arguments for errors
    roll::check_rows_equal(n_rows_xy, y.size(), "x", "y");
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_xy, width, weights, "'x' (and 'y')");
    bool status = roll::check_lambda(weights, n_rows_xy, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // cbind x and y variables
    NumericMatrix data(n_rows_xy, n_cols_x);
    std::copy(xx.begin(), xx.end(), data.begin());
    std::copy(y.begin(), y.end(), data.begin() + n_rows_xy * (n_cols_x - 1));
    
    // default 'complete_obs' argument is 'true',
    // otherwise check argument for errors
    if (complete_obs) {
      arma_any_na = any_na_x(data);
    } else {
      
      warning("'complete_obs = FALSE' is not supported");
      arma_any_na = any_na_x(data);
      
    }
    
    // compute rolling crossproducts
    if (status && online) {
      
      roll::RollCrossProdOnlineMatXX roll_cov_online(data, n, n_rows_xy, n_cols_x, width,
                                                     weights, intercept, false, min_obs,
                                                     arma_any_na, na_restore,
                                                     arma_n_obs, arma_sum_w, arma_mean,
                                                     arma_cov);
      parallelFor(0, n_cols_x, roll_cov_online);
      
    } else {
      
      roll::RollCrossProdOfflineMatXX roll_cov_offline(data, n, n_rows_xy, n_cols_x, width,
                                                       weights, intercept, false, min_obs,
                                                       arma_any_na, na_restore,
                                                       arma_n_obs, arma_sum_w, arma_mean,
                                                       arma_cov);
      parallelFor(0, n_rows_xy * n_cols_x * (n_cols_x + 1) / 2, roll_cov_offline);
      
    }
    
    // compute rolling linear models
    if (intercept) {
      
      arma::mat arma_coef(n_rows_xy, n_cols_x);
      arma::mat arma_se(n_rows_xy, n_cols_x);
      roll::RollLmMatInterceptTRUE roll_lm_slices(arma_cov, n, n_rows_xy, n_cols_x, width,
                                                  arma_n_obs, arma_sum_w, arma_mean,
                                                  arma_coef, arma_rsq, arma_se);
      parallelFor(0, n_rows_xy, roll_lm_slices);
      
      result = List::create(Named("coefficients") = arma_coef,
                            Named("r.squared") = arma_rsq,
                            Named("std.error") = arma_se);
      
    } else if (!intercept) {
      
      arma::mat arma_coef(n_rows_xy, n_cols_x - 1);
      arma::mat arma_se(n_rows_xy, n_cols_x - 1);
      roll::RollLmMatInterceptFALSE roll_lm_slices(arma_cov, n, n_rows_xy, n_cols_x, width,
                                                   arma_n_obs, arma_sum_w,
                                                   arma_coef, arma_rsq, arma_se);
      parallelFor(0, n_rows_xy, roll_lm_slices);
      
      // create and return a list
      result = List::create(Named("coefficients") = arma_coef,
                            Named("r.squared") = arma_rsq,
                            Named("std.error") = arma_se);
      
    }
    
    return result;
    
  } else {
    
    NumericVector xx(x);
    int n = weights.size();
    int n_rows_xy = xx.size();
    int n_cols_x = 1 + 1;
    arma::uvec arma_any_na(n_rows_xy);
    arma::vec arma_n_obs(n_rows_xy);
    arma::vec arma_sum_w(n_rows_xy);
    arma::mat arma_mean(n_rows_xy, n_cols_x);
    arma::cube arma_cov(n_cols_x, n_cols_x, n_rows_xy);
    arma::vec arma_rsq(n_rows_xy);
    List result(3);
    
    // check 'x' and 'y' arguments for errors
    roll::check_rows_equal(n_rows_xy, y.size(), "x", "y");
    
    // check 'width' argument for errors
    roll::check_pos_int(width, "width");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    roll::check_weights(n_rows_xy, width, weights, "'x' (and 'y')");
    bool status = roll::check_lambda(weights, n_rows_xy, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    roll::check_pos_int(min_obs, "min_obs");
    
    // cbind x and y variables
    NumericMatrix data(n_rows_xy, n_cols_x);
    std::copy(xx.begin(), xx.end(), data.begin());
    std::copy(y.begin(), y.end(), data.begin() + n_rows_xy * (n_cols_x - 1));
    
    // default 'complete_obs' argument is 'true',
    // otherwise check argument for errors
    if (complete_obs) {
      arma_any_na = any_na_x(data);
    } else {
      
      warning("'complete_obs = FALSE' is not supported");
      arma_any_na = any_na_x(data);
      
    }
    
    // compute rolling crossproducts
    if (status && online) {
      
      roll::RollCrossProdOnlineMatXX roll_cov_online(data, n, n_rows_xy, n_cols_x, width,
                                                     weights, intercept, false, min_obs,
                                                     arma_any_na, na_restore,
                                                     arma_n_obs, arma_sum_w, arma_mean,
                                                     arma_cov);
      parallelFor(0, n_cols_x, roll_cov_online);
      
    } else {
      
      roll::RollCrossProdOfflineMatXX roll_cov_offline(data, n, n_rows_xy, n_cols_x, width,
                                                       weights, intercept, false, min_obs,
                                                       arma_any_na, na_restore,
                                                       arma_n_obs, arma_sum_w, arma_mean,
                                                       arma_cov);
      parallelFor(0, n_rows_xy * n_cols_x * (n_cols_x + 1) / 2, roll_cov_offline);
      
    }
    
    // compute rolling linear models
    if (intercept) {
      
      arma::mat arma_coef(n_rows_xy, n_cols_x);
      arma::mat arma_se(n_rows_xy, n_cols_x);
      roll::RollLmMatInterceptTRUE roll_lm_slices(arma_cov, n, n_rows_xy, n_cols_x, width,
                                                  arma_n_obs, arma_sum_w, arma_mean,
                                                  arma_coef, arma_rsq, arma_se);
      parallelFor(0, n_rows_xy, roll_lm_slices);
      
      result = List::create(Named("coefficients") = arma_coef,
                            Named("r.squared") = arma_rsq,
                            Named("std.error") = arma_se);
      
    } else if (!intercept) {
      
      arma::vec arma_coef(n_rows_xy);
      arma::vec arma_se(n_rows_xy);
      roll::RollLmVecInterceptFALSE roll_lm_slices(arma_cov, n, n_rows_xy, width,
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
  
}

// [[Rcpp::export(.roll_lm)]]
List roll_lm(const SEXP& x, const SEXP& y,
             const int& width, const arma::vec& weights,
             const bool& intercept, const int& min_obs,
             const bool& complete_obs, const bool& na_restore,
             const bool& online) {
  
  if (Rf_isMatrix(x) && Rf_isMatrix(y)) {
    
    NumericMatrix xx(x);
    NumericMatrix yy(y);
    int n_rows_xy = xx.nrow();
    int n_cols_x = xx.ncol();
    int n_cols_y = yy.ncol();
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
      
      result_z = roll_lm_z(xx, yy(_, 0), width,
                           weights, intercept, min_obs,
                           complete_obs, na_restore,
                           online);
      
      arma::mat arma_coef_z = result_z[0];
      arma::mat arma_rsq_z = result_z[1];
      arma::mat arma_se_z = result_z[2];
      
      // create and return a matrix or xts object for coefficients
      NumericVector coef(wrap(arma_coef_z));
      coef.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
      List dimnames_x = xx.attr("dimnames");
      roll::xts_attr(coef, xx, dimnames_lm_x(dimnames_x, n_cols_x, intercept));
      
      // create and return a matrix or xts object for r-squareds
      NumericVector rsq(wrap(arma_rsq_z));
      roll::rsq_attr(rsq, n_rows_xy, xx, dimnames_x);
      
      // create and return a matrix or xts object for standard errors
      NumericVector se(wrap(arma_se_z));
      se.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
      roll::xts_attr(se, xx, coef.attr("dimnames"));
      
      // create and return a list
      result = List::create(Named("coefficients") = coef,
                            Named("r.squared") = rsq,
                            Named("std.error") = se);
      
    } else {
      
      for (int z = 0; z < n_cols_y; z++) {
        
        result_z = roll_lm_z(xx, yy(_, z), width,
                             weights, intercept, min_obs,
                             complete_obs, na_restore,
                             online);
        
        arma::mat arma_coef_z = result_z[0];
        arma::mat arma_rsq_z = result_z[1];
        arma::mat arma_se_z = result_z[2];
        
        // create and return a matrix or xts object for coefficients
        NumericVector coef(wrap(arma_coef_z));
        coef.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
        List dimnames_x = xx.attr("dimnames");
        roll::xts_attr(coef, xx, dimnames_lm_x(dimnames_x, n_cols_x, intercept));
        
        // create and return a matrix or xts object for r-squareds
        NumericVector rsq(wrap(arma_rsq_z));
        roll::rsq_attr(rsq, n_rows_xy, xx, dimnames_x);
        
        // create and return a matrix or xts object for standard errors
        NumericVector se(wrap(arma_se_z));
        se.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
        roll::xts_attr(se, xx, coef.attr("dimnames"));
        
        result_coef(z) = coef;
        result_rsq(z) = rsq;
        result_se(z) = se;
        
      }
      
      // add names to each list
      List dimnames_y = yy.attr("dimnames");
      result_coef.attr("names") = dimnames_lm_y(dimnames_y, n_cols_y);
      result_rsq.attr("names") = result_coef.attr("names");
      result_se.attr("names") = result_coef.attr("names");
      
      // create and return a list
      result = List::create(Named("coefficients") = result_coef,
                            Named("r.squared") = result_rsq,
                            Named("std.error") = result_se);
      
    }
    
    return result;
    
  } else if (Rf_isMatrix(x)) {
    
    NumericMatrix xx(x);
    NumericVector yy(y);
    
    int n_rows_xy = xx.nrow();
    int n_cols_x = xx.ncol();
    List result_coef(1);
    List result_rsq(1);
    List result_se(1);
    List result_z(3);
    List result(3);
    
    if (intercept) {
      n_cols_x += 1;
    }
    
    // create a list of matrices
    result_z = roll_lm_z(xx, yy, width,
                         weights, intercept, min_obs,
                         complete_obs, na_restore,
                         online);
    
    arma::mat arma_coef_z = result_z[0];
    arma::mat arma_rsq_z = result_z[1];
    arma::mat arma_se_z = result_z[2];
    
    // create and return a matrix or xts object for coefficients
    NumericVector coef(wrap(arma_coef_z));
    coef.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
    List dimnames_x = xx.attr("dimnames");
    roll::xts_attr(coef, xx, dimnames_lm_x(dimnames_x, n_cols_x, intercept));
    
    // create and return a matrix or xts object for r-squareds
    NumericVector rsq(wrap(arma_rsq_z));
    roll::rsq_attr(rsq, n_rows_xy, xx, dimnames_x);
    
    // create and return a matrix or xts object for standard errors
    NumericVector se(wrap(arma_se_z));
    se.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
    roll::xts_attr(se, xx, coef.attr("dimnames"));
    
    // create and return a list
    result = List::create(Named("coefficients") = coef,
                          Named("r.squared") = rsq,
                          Named("std.error") = se);
    
    return result;
    
  } else if (Rf_isMatrix(y)) {
    
    NumericVector xx(x);
    NumericMatrix yy(y);
    // xx.attr("dim") = IntegerVector::create(xx.size(), 1);
    // NumericMatrix xxx(wrap(xx));
    NumericMatrix xxx(xx.size(), 1, xx.begin());
    
    int n_rows_xy = xxx.nrow();
    int n_cols_x = xxx.ncol();
    int n_cols_y = yy.ncol();
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
      
      result_z = roll_lm_z(xxx, yy(_, 0), width,
                           weights, intercept, min_obs,
                           complete_obs, na_restore,
                           online);
      
      arma::mat arma_coef_z = result_z[0];
      arma::mat arma_rsq_z = result_z[1];
      arma::mat arma_se_z = result_z[2];
      
      // create and return a matrix or xts object for coefficients
      NumericVector coef(wrap(arma_coef_z));
      coef.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
      List dimnames_x = xx.attr("dimnames");
      roll::xts_attr(coef, yy, dimnames_lm_x(dimnames_x, n_cols_x, intercept));
      
      // create and return a matrix or xts object for r-squareds
      NumericVector rsq(wrap(arma_rsq_z));
      roll::rsq_attr(rsq, n_rows_xy, yy, dimnames_x);
      
      // create and return a matrix or xts object for standard errors
      NumericVector se(wrap(arma_se_z));
      se.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
      roll::xts_attr(se, yy, coef.attr("dimnames"));
      
      // create and return a list
      result = List::create(Named("coefficients") = coef,
                            Named("r.squared") = rsq,
                            Named("std.error") = se);
      
    } else {
      
      for (int z = 0; z < n_cols_y; z++) {
        
        result_z = roll_lm_z(xxx, yy(_, z), width,
                             weights, intercept, min_obs,
                             complete_obs, na_restore,
                             online);
        
        arma::mat arma_coef_z = result_z[0];
        arma::mat arma_rsq_z = result_z[1];
        arma::mat arma_se_z = result_z[2];
        
        // create and return a matrix or xts object for coefficients
        NumericVector coef(wrap(arma_coef_z));
        coef.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
        List dimnames_x = xxx.attr("dimnames");
        roll::xts_attr(coef, yy, dimnames_lm_x(dimnames_x, n_cols_x, intercept));
        
        // create and return a matrix or xts object for r-squareds
        NumericVector rsq(wrap(arma_rsq_z));
        roll::rsq_attr(rsq, n_rows_xy, yy, dimnames_x);
        
        // create and return a matrix or xts object for standard errors
        NumericVector se(wrap(arma_se_z));
        se.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
        roll::xts_attr(se, yy, coef.attr("dimnames"));
        
        result_coef(z) = coef;
        result_rsq(z) = rsq;
        result_se(z) = se;
        
      }
      
      // add names to each list
      List dimnames_y = yy.attr("dimnames");
      result_coef.attr("names") = dimnames_lm_y(dimnames_y, n_cols_y);
      result_rsq.attr("names") = result_coef.attr("names");
      result_se.attr("names") = result_coef.attr("names");
      
      // create and return a list
      result = List::create(Named("coefficients") = result_coef,
                            Named("r.squared") = result_rsq,
                            Named("std.error") = result_se);
      
    }
    
    return result;
    
  } else {
    
    NumericVector xx(x);
    NumericVector yy(y);
    
    int n_rows_xy = xx.size();
    int n_cols_x = 1;
    int n_cols_y = 1;
    List result_coef(n_cols_y);
    List result_rsq(n_cols_y);
    List result_se(n_cols_y);
    List result_z(3);
    List result(3);
    
    if (intercept) {
      n_cols_x += 1;
    }
    
    // create a list of matrices
    result_z = roll_lm_z(xx, yy, width,
                         weights, intercept, min_obs,
                         complete_obs, na_restore,
                         online);
    
    arma::mat arma_coef_z = result_z[0];
    arma::mat arma_rsq_z = result_z[1];
    arma::mat arma_se_z = result_z[2];
    
    if (intercept) {
      
      // create and return a matrix or xts object for coefficients
      NumericVector coef(wrap(arma_coef_z));
      coef.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
      List dimnames_x = xx.attr("dimnames");
      roll::xts_attr(coef, xx, dimnames_lm_x(dimnames_x, n_cols_x, intercept));
      
      // create and return a matrix or xts object for r-squareds
      NumericVector rsq(wrap(arma_rsq_z));
      roll::rsq_attr(rsq, n_rows_xy, xx, dimnames_x);
      
      // create and return a matrix or xts object for standard errors
      NumericVector se(wrap(arma_se_z));
      se.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
      roll::xts_attr(se, xx, coef.attr("dimnames"));
      
      // create and return a list
      result = List::create(Named("coefficients") = coef,
                            Named("r.squared") = rsq,
                            Named("std.error") = se);
      
      return result;
      
    } else {
      
      // create and return a vector object for coefficients
      NumericVector coef(wrap(arma_coef_z));
      roll::vec_attr(coef, xx);
      
      // create and return a vector object for r-squareds
      NumericVector rsq(wrap(arma_rsq_z));
      roll::vec_attr(rsq, xx);
      
      // create and return a vector object for standard errors
      NumericVector se(wrap(arma_se_z));
      roll::vec_attr(se, xx);
      
      // create and return a list
      result = List::create(Named("coefficients") = coef,
                            Named("r.squared") = rsq,
                            Named("std.error") = se);
      
      return result;
      
    }
    
  }
  
}
