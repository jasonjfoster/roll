#ifndef ROLL_ATTR_H
#define ROLL_ATTR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

namespace roll {

// xts attributes: copy dimnames, index, .indexCLASS, .indexTZ, tclass, tzone, class
template <typename T, typename S>
inline void xts_attr(T& target, const S& source) {

  target.attr("dimnames") = source.attr("dimnames");
  target.attr("index") = source.attr("index");
  target.attr(".indexCLASS") = source.attr(".indexCLASS");
  target.attr(".indexTZ") = source.attr(".indexTZ");
  target.attr("tclass") = source.attr("tclass");
  target.attr("tzone") = source.attr("tzone");
  target.attr("class") = source.attr("class");

}

// xts attributes: copy xts attrs from source with supplied dimnames override
template <typename T, typename S>
inline void xts_attr(T& target, const S& source,
                     SEXP dimnames) {

  target.attr("dimnames") = dimnames;
  target.attr("index") = source.attr("index");
  target.attr(".indexCLASS") = source.attr(".indexCLASS");
  target.attr(".indexTZ") = source.attr(".indexTZ");
  target.attr("tclass") = source.attr("tclass");
  target.attr("tzone") = source.attr("tzone");
  target.attr("class") = source.attr("class");

}

// vector attributes: strip dim, conditionally copy names, copy index and class
template <typename T, typename S>
inline void vec_attr(T& target, const S& source) {

  target.attr("dim") = R_NilValue;
  List names = source.attr("names");

  if (names.size() > 0) {
    target.attr("names") = names;
  }

  target.attr("index") = source.attr("index");
  target.attr("class") = source.attr("class");

}

// cube attributes: set dim to (n_d1, n_d2, n_d3), conditionally attach dimnames
template <typename T>
inline void cube_attr(T& target, const int& n_d1,
                      const int& n_d2, const int& n_d3,
                      const List& dimnames_x, const List& dimnames_y) {

  target.attr("dim") = IntegerVector::create(n_d1, n_d2, n_d3);

  if ((dimnames_x.size() > 1) && (dimnames_y.size() > 1)) {
    target.attr("dimnames") = List::create(dimnames_x[1], dimnames_y[1]);
  } else if (dimnames_x.size() > 1) {
    target.attr("dimnames") = List::create(dimnames_x[1], R_NilValue);
  } else if (dimnames_y.size() > 1) {
    target.attr("dimnames") = List::create(R_NilValue, dimnames_y[1]);
  }

}

// rsq attributes: set dim to (n_rows, 1), attach "R-squared" column name and xts attributes
template <typename T, typename S>
inline void rsq_attr(T& target, const int& n_rows,
                     const S& source, SEXP dimnames_x) {

  target.attr("dim") = IntegerVector::create(n_rows, 1);
  List dimnames(dimnames_x);

  if (dimnames.size() > 1) {
    xts_attr(target, source, List::create(dimnames[0], "R-squared"));
  } else {
    xts_attr(target, source, List::create(R_NilValue, "R-squared"));
  }

}

// strip dim: remove the dim attribute
template <typename T>
inline void strip_dim(T& target) {

  target.attr("dim") = R_NilValue;

}

}

#endif