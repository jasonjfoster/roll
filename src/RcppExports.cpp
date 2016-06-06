// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// roll_mean
NumericMatrix roll_mean(const NumericMatrix& data, const int& width, const arma::vec& weights, const int& min_obs, const bool& complete_obs, const bool& na_restore, const std::string& parallel_for);
RcppExport SEXP roll_roll_mean(SEXP dataSEXP, SEXP widthSEXP, SEXP weightsSEXP, SEXP min_obsSEXP, SEXP complete_obsSEXP, SEXP na_restoreSEXP, SEXP parallel_forSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const int& >::type width(widthSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const int& >::type min_obs(min_obsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type complete_obs(complete_obsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type na_restore(na_restoreSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type parallel_for(parallel_forSEXP);
    __result = Rcpp::wrap(roll_mean(data, width, weights, min_obs, complete_obs, na_restore, parallel_for));
    return __result;
END_RCPP
}
// roll_var
NumericMatrix roll_var(const NumericMatrix& data, const int& width, const arma::vec& weights, const bool& center, const int& min_obs, const bool& complete_obs, const bool& na_restore, const std::string& parallel_for);
RcppExport SEXP roll_roll_var(SEXP dataSEXP, SEXP widthSEXP, SEXP weightsSEXP, SEXP centerSEXP, SEXP min_obsSEXP, SEXP complete_obsSEXP, SEXP na_restoreSEXP, SEXP parallel_forSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const int& >::type width(widthSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type center(centerSEXP);
    Rcpp::traits::input_parameter< const int& >::type min_obs(min_obsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type complete_obs(complete_obsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type na_restore(na_restoreSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type parallel_for(parallel_forSEXP);
    __result = Rcpp::wrap(roll_var(data, width, weights, center, min_obs, complete_obs, na_restore, parallel_for));
    return __result;
END_RCPP
}
// roll_sd
NumericMatrix roll_sd(const NumericMatrix& data, const int& width, const arma::vec& weights, const bool& center, const int& min_obs, const bool& complete_obs, const bool& na_restore, const std::string& parallel_for);
RcppExport SEXP roll_roll_sd(SEXP dataSEXP, SEXP widthSEXP, SEXP weightsSEXP, SEXP centerSEXP, SEXP min_obsSEXP, SEXP complete_obsSEXP, SEXP na_restoreSEXP, SEXP parallel_forSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const int& >::type width(widthSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type center(centerSEXP);
    Rcpp::traits::input_parameter< const int& >::type min_obs(min_obsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type complete_obs(complete_obsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type na_restore(na_restoreSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type parallel_for(parallel_forSEXP);
    __result = Rcpp::wrap(roll_sd(data, width, weights, center, min_obs, complete_obs, na_restore, parallel_for));
    return __result;
END_RCPP
}
// roll_cov
NumericVector roll_cov(const NumericMatrix& data, const int& width, const arma::vec& weights, const bool& center, const bool& scale, const int& min_obs, const bool& complete_obs, const bool& na_restore, const std::string& parallel_for);
RcppExport SEXP roll_roll_cov(SEXP dataSEXP, SEXP widthSEXP, SEXP weightsSEXP, SEXP centerSEXP, SEXP scaleSEXP, SEXP min_obsSEXP, SEXP complete_obsSEXP, SEXP na_restoreSEXP, SEXP parallel_forSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const int& >::type width(widthSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type center(centerSEXP);
    Rcpp::traits::input_parameter< const bool& >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const int& >::type min_obs(min_obsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type complete_obs(complete_obsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type na_restore(na_restoreSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type parallel_for(parallel_forSEXP);
    __result = Rcpp::wrap(roll_cov(data, width, weights, center, scale, min_obs, complete_obs, na_restore, parallel_for));
    return __result;
END_RCPP
}
// roll_lm
List roll_lm(const NumericMatrix& x, const NumericMatrix& y, const int& width, const arma::vec& weights, const bool& center_x, const bool& center_y, const bool& scale_x, const bool& scale_y, const int& min_obs, const bool& complete_obs, const bool& na_restore, const std::string& parallel_for);
RcppExport SEXP roll_roll_lm(SEXP xSEXP, SEXP ySEXP, SEXP widthSEXP, SEXP weightsSEXP, SEXP center_xSEXP, SEXP center_ySEXP, SEXP scale_xSEXP, SEXP scale_ySEXP, SEXP min_obsSEXP, SEXP complete_obsSEXP, SEXP na_restoreSEXP, SEXP parallel_forSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const int& >::type width(widthSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type center_x(center_xSEXP);
    Rcpp::traits::input_parameter< const bool& >::type center_y(center_ySEXP);
    Rcpp::traits::input_parameter< const bool& >::type scale_x(scale_xSEXP);
    Rcpp::traits::input_parameter< const bool& >::type scale_y(scale_ySEXP);
    Rcpp::traits::input_parameter< const int& >::type min_obs(min_obsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type complete_obs(complete_obsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type na_restore(na_restoreSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type parallel_for(parallel_forSEXP);
    __result = Rcpp::wrap(roll_lm(x, y, width, weights, center_x, center_y, scale_x, scale_y, min_obs, complete_obs, na_restore, parallel_for));
    return __result;
END_RCPP
}
// roll_eigen
List roll_eigen(const NumericMatrix& data, const int& width, const arma::vec& weights, const bool& center, const bool& scale, const int& min_obs, const bool& complete_obs, const bool& na_restore, const std::string& parallel_for);
RcppExport SEXP roll_roll_eigen(SEXP dataSEXP, SEXP widthSEXP, SEXP weightsSEXP, SEXP centerSEXP, SEXP scaleSEXP, SEXP min_obsSEXP, SEXP complete_obsSEXP, SEXP na_restoreSEXP, SEXP parallel_forSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const int& >::type width(widthSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type center(centerSEXP);
    Rcpp::traits::input_parameter< const bool& >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const int& >::type min_obs(min_obsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type complete_obs(complete_obsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type na_restore(na_restoreSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type parallel_for(parallel_forSEXP);
    __result = Rcpp::wrap(roll_eigen(data, width, weights, center, scale, min_obs, complete_obs, na_restore, parallel_for));
    return __result;
END_RCPP
}
// roll_pcr
List roll_pcr(const NumericMatrix& x, const NumericMatrix& y, const int& width, const arma::uvec& comps, const arma::vec& weights, const bool& center_x, const bool& center_y, const bool& scale_x, const bool& scale_y, const int& min_obs, const bool& complete_obs, const bool& na_restore, const std::string& parallel_for);
RcppExport SEXP roll_roll_pcr(SEXP xSEXP, SEXP ySEXP, SEXP widthSEXP, SEXP compsSEXP, SEXP weightsSEXP, SEXP center_xSEXP, SEXP center_ySEXP, SEXP scale_xSEXP, SEXP scale_ySEXP, SEXP min_obsSEXP, SEXP complete_obsSEXP, SEXP na_restoreSEXP, SEXP parallel_forSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const int& >::type width(widthSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type comps(compsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type center_x(center_xSEXP);
    Rcpp::traits::input_parameter< const bool& >::type center_y(center_ySEXP);
    Rcpp::traits::input_parameter< const bool& >::type scale_x(scale_xSEXP);
    Rcpp::traits::input_parameter< const bool& >::type scale_y(scale_ySEXP);
    Rcpp::traits::input_parameter< const int& >::type min_obs(min_obsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type complete_obs(complete_obsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type na_restore(na_restoreSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type parallel_for(parallel_forSEXP);
    __result = Rcpp::wrap(roll_pcr(x, y, width, comps, weights, center_x, center_y, scale_x, scale_y, min_obs, complete_obs, na_restore, parallel_for));
    return __result;
END_RCPP
}
// roll_vif
NumericMatrix roll_vif(const NumericMatrix& data, const int& width, const arma::uvec& comps, const arma::vec& weights, const bool& center, const bool& scale, const int& min_obs, const bool& complete_obs, const bool& na_restore, const std::string& parallel_for);
RcppExport SEXP roll_roll_vif(SEXP dataSEXP, SEXP widthSEXP, SEXP compsSEXP, SEXP weightsSEXP, SEXP centerSEXP, SEXP scaleSEXP, SEXP min_obsSEXP, SEXP complete_obsSEXP, SEXP na_restoreSEXP, SEXP parallel_forSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const int& >::type width(widthSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type comps(compsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type center(centerSEXP);
    Rcpp::traits::input_parameter< const bool& >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const int& >::type min_obs(min_obsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type complete_obs(complete_obsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type na_restore(na_restoreSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type parallel_for(parallel_forSEXP);
    __result = Rcpp::wrap(roll_vif(data, width, comps, weights, center, scale, min_obs, complete_obs, na_restore, parallel_for));
    return __result;
END_RCPP
}
