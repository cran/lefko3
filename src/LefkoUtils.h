#ifndef LEFKOUTILS_H
#define LEFKOUTILS_H

#include <Rcpp.h>
#include <RcppArmadillo.h>

namespace LefkoUtils {
  Rcpp::NumericVector concat_dbl(Rcpp::NumericVector, Rcpp::NumericVector);
  Rcpp::IntegerVector concat_int(Rcpp::IntegerVector, Rcpp::IntegerVector);
  Rcpp::StringVector concat_str(Rcpp::StringVector, Rcpp::StringVector);
  
  bool stringcompare_hard(std::string, std::string);
  Rcpp::List stringcompare_soft(std::string, std::string);
  bool stringcompare_simple(std::string, std::string, bool);
  bool stringcompare_x(std::string, std::string, std::string);
  
  Rcpp::CharacterVector stringsort(Rcpp::CharacterVector);
  Rcpp::IntegerVector int_sort(Rcpp::IntegerVector);
  Rcpp::IntegerMatrix refsort_num(Rcpp::NumericMatrix, Rcpp::NumericVector);
  Rcpp::IntegerVector refsort_str(Rcpp::CharacterVector, Rcpp::CharacterVector);
  
  arma::vec flagrantcrap(arma::mat, arma::uvec);
  arma::vec moreflagrantcrap(arma::mat);
  
  arma::sp_mat spmat_log(arma::sp_mat);
  
  Rcpp::IntegerVector shrink(Rcpp::IntegerVector);
}


#endif
