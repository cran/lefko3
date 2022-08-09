#include <RcppArmadillo.h>
#include "LefkoUtils.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Append NumericVector to the End of Another NumericVector
//' 
//' This function appends one NumericVector fully to another.
//' 
//' @name concat_dbl
//' 
//' @param A Any NumericVector.
//' @param B Any other NumericVector.
//' 
//' @return Returns a new NumericVector with elements of vector A followed by
//' elements of vector B.
//'
//' @keywords internal
//' @noRd
Rcpp::NumericVector LefkoUtils::concat_dbl(Rcpp::NumericVector x, Rcpp::NumericVector y) {
  
  std::vector<double> xconv = as<std::vector<double> >(x);
  std::vector<double> yconv = as<std::vector<double> >(y);
  std::vector<double> z(xconv.size() + yconv.size());
  
  std::copy(xconv.begin(), xconv.end(), z.begin());
  std::copy(yconv.begin(), yconv.end(), z.begin() + xconv.size());
  
  Rcpp::NumericVector zconv(z.begin(), z.end());
  
  return(zconv);
}

//' Append IntegerVector to the End of Another IntegerVector
//' 
//' Returns a new IntegerVector with elements of vector A followed by
//' elements of vector B.
//' 
//' @name concat_int
//' 
//' @param A Any IntegerVector.
//' @param B Any other IntegerVector.
//' 
//' @return Returns a new IntegerVector with elements of vector A followed by
//' elements of vector B.
//' 
//' @keywords internal
//' @noRd
Rcpp::IntegerVector LefkoUtils::concat_int(Rcpp::IntegerVector x, Rcpp::IntegerVector y) {
  
  std::vector<long long int> xconv = as<std::vector<long long int> >(x);
  std::vector<long long int> yconv = as<std::vector<long long int> >(y);
  std::vector<long long int> z(x.size() + y.size());
  
  std::copy(xconv.begin(), xconv.end(), z.begin());
  std::copy(yconv.begin(), yconv.end(), z.begin() + xconv.size());
  
  Rcpp::IntegerVector zconv(z.begin(), z.end());
  
  return(zconv);
}

//' Append StringVector to the End of Another StringVector
//' 
//' Returns a new StringVector with elements of vector A followed by
//' elements of vector B.
//' 
//' @name concat_str
//' 
//' @param A Any StringVector.
//' @param B Any other StringVector.
//' 
//' @return Returns a new StringVector with elements of vector A followed by
//' elements of vector B.
//' 
//' @keywords internal
//' @noRd
Rcpp::StringVector LefkoUtils::concat_str(Rcpp::StringVector x, Rcpp::StringVector y) {
  
  std::vector<std::string> xconv = as<std::vector<std::string> >(x);
  std::vector<std::string> yconv = as<std::vector<std::string> >(y);
  std::vector<std::string> z(x.size() + y.size());
  
  std::copy(x.begin(), x.end(), z.begin());
  std::copy(y.begin(), y.end(), z.begin() + x.size());
  
  Rcpp::StringVector zconv(z.begin(), z.end());
  
  return(zconv);
}

//' Compares Two Strings Literally
//' 
//' This function compares two strings element by element. Returns \code{FALSE}
//' in case of any differences whatsoever.
//' 
//' @name stringcompare_hard
//' 
//' @param str1 The first string
//' @param str2 The second string
//' 
//' @return A logical value. In case of any difference at all, it will return
//' \code{FALSE}.
//' 
//' @keywords internal
//' @noRd
bool LefkoUtils::stringcompare_hard(std::string str1, std::string str2) {
  int str1_length = str1.size();
  int str2_length = str2.size();
  bool same = true;
  
  if (str1_length == str2_length && str1_length > 0) {
    for (int i = 0; i < str1_length; i++) {
      if (str1[i] != str2[i]) {
        same = false;
      }
    }
  } else if (str1_length != str2_length) {
    same = false;
  }
  
  return same;
}

//' Compares Two Strings, Assessing Inclusion
//' 
//' This function compares two strings, and will assess whether \code{str2} is
//' contained within \code{str1}.
//' 
//' @name stringcompare_soft
//' 
//' @param str1 The first string
//' @param str2 The second string
//' 
//' @return A list of two values. The first is a logical value indicating
//' whether \code{str2} occurs within \code{str1}. The second element is an
//' integer indicating at what element of \code{str1} \code{str2} begins.
//' \code{FALSE}.
//' 
//' @keywords internal
//' @noRd
Rcpp::List LefkoUtils::stringcompare_soft(std::string str1, std::string str2) {
  int str1_length = str1.size();
  int str2_length = str2.size();
  int rem_check {0};
  bool same = false;
  unsigned int start_index {0};
  
  if (str1_length >= str2_length && str2_length > 0) {
    for (int i = 0; i < str1_length; i++) {
      if (str1[i] != str2[rem_check]) {
        rem_check = 0;
        // same = false;
      } else {
        if (rem_check == 0) start_index = i;
        rem_check += 1;
        if (rem_check >= str2_length) break;
      }
    }
    
    if (rem_check == str2_length) {
      same = true;
    }
  }
  
  Rcpp::List output = Rcpp::List::create(_["contains"] = same, _["start_index"] = start_index);
  
  return output;
}

//' Compares Two Strings, Assessing Inclusion
//' 
//' This function compares two strings, and will assess whether \code{str2} is
//' contained within \code{str1}. It is a simpler version of 
//' \code{stringcompare_soft()} that yields only the logical result.
//' 
//' @name stringcompare_simple
//' 
//' @param str1 The first string
//' @param str2 The second string
//' @param lower A logical value indicating whether to change all inputs to
//' lower case before checking.
//' 
//' @return A logical value indicating whether \code{str2} occurs within
//' \code{str1}.
//' 
//' @keywords internal
//' @noRd
bool LefkoUtils::stringcompare_simple(std::string str1, std::string str2, bool lower = false) {
  int str1_length = str1.size();
  int str2_length = str2.size();
  int rem_check {0};
  bool same = false;
  
  if (str1_length >= str2_length && str2_length > 0) {
    for (int i = 0; i < str1_length; i++) {
      if (!lower) {
        if (str1[i] != str2[rem_check]) {
          rem_check = 0;
        } else {
          rem_check += 1;
          if (rem_check >= str2_length) break;
        }
      } else {
        if (tolower(str1[i]) != tolower(str2[rem_check])) {
          rem_check = 0;
        } else {
          rem_check += 1;
          if (rem_check >= str2_length) break;
        }
      }
    }
    
    if (rem_check == str2_length) {
      same = true;
    }
  }
  
  return same;
}

//' Compares Three Strings for Interaction Notation
//' 
//' This function compares checks to see if one string is composed of the other
//' two strings in R's interaction notation.
//' 
//' @name stringcompare_x
//' 
//' @param str1 The first string. Used for comparison.
//' @param str2 The second string. Will be incorporated into interaction format.
//' @param str3 The third string. Will be incorporated into interaction format.
//' 
//' @return A logical value. In case of any difference at all, it will return
//' \code{FALSE}.
//' 
//' @keywords internal
//' @noRd
bool LefkoUtils::stringcompare_x(std::string str1, std::string str2, std::string str3) {
  int str1_length = str1.size();
  int str2_length = str2.size();
  int str3_length = str3.size();
  int combined_length = str2_length + str3_length + 1;
  bool same = false;
  bool same1 = true;
  bool same2 = true;
  
  if (str1_length == combined_length && str1_length > 0) {
    std::string x1 = str2;
    x1 += ":";
    x1 += str3;
    
    std::string x2 = str3;
    x2 += ":";
    x2 += str2;
    
    for (int i = 0; i < str1_length; i++) {
      if (str1[i] != x1[i]) {
        same1 = false;
      }
      if (str1[i] != x2[i]) {
        same2 = false;
      }
    }
  } else {
    same1 = false;
    same2 = false;
  }
  
  if (same1 || same2) same = true;
  
  return same;
}

//' Sort String Elements
//' 
//' This function is based on code obtained from R Bloggers
//' (see https://www.r-bloggers.com/2013/01/handling-strings-with-rcpp/). It
//' sorts the elements of a string vector in alphabetical order.
//' 
//' @name stringsort
//' 
//' @param string_input A string vector.
//' 
//' @return The sorted string vector.
//' 
//' @keywords internal
//' @noRd
Rcpp::CharacterVector LefkoUtils::stringsort(Rcpp::CharacterVector string_input ) {
  int len = string_input.size();
  
  std::vector<std::string> converted(len);
  for (int i=0; i < len; i++) converted[i] = as<std::string>(string_input(i));
  std::sort( converted.begin(), converted.end() );
  
  Rcpp::CharacterVector new_converted(len);
  new_converted = converted;
  
  return new_converted;
}

//' Sort Integer Elements
//' 
//' This function is based on code obtained from the Rcpp Gallery by Ross
//' Bennett (see https://gallery.rcpp.org/articles/sorting/). It sorts the
//' elements of an integer vector.
//' 
//' @name int_sort
//' 
//' @param int_input An integer vector.
//' 
//' @return The sorted integer vector.
//' 
//' @keywords internal
//' @noRd
Rcpp::IntegerVector LefkoUtils::int_sort(Rcpp::IntegerVector x) {
   Rcpp::IntegerVector y = clone(x);
   std::sort(y.begin(), y.end());
   
   return y;
}

//' Function to Index a Numeric Vector According to a Reference Vector
//' 
//' Function \code{refsort_num()} takes a numeric matrix and replaces it with an
//' integer vector showing the position of each element in the input vector
//' within the reference vector.
//' 
//' @name refsort_num
//' 
//' @param vec The matrix to index
//' @param ref The vector to use as a reference
//' 
//' @return An integer vector with integers referring to elements in vector
//' \code{ref}.
//' 
//' @keywords internal
//' @noRd
Rcpp::IntegerMatrix LefkoUtils::refsort_num(Rcpp::NumericMatrix vec, Rcpp::NumericVector ref) {
  int vec_length = vec.length();
  int ref_length = ref.length();
  
  Rcpp::IntegerMatrix output(vec.nrow(), vec.ncol());
  
  for (int i = 0; i < vec_length; i++) {
    for (int j = 0; j < ref_length; j++) {
      if (vec[i] == ref[j]) output[i] = j + 1;
    }
  }
  
  return output;
}

//' Function to Index a Numeric Vector According to a Reference Vector
//' 
//' Function \code{refsort_str()} takes a string vector and replaces it with an
//' integer vector showing the position of each element in the input vector
//' within the reference vector.
//' 
//' @name refsort_str
//' 
//' @param vec The vector to index
//' @param ref The vector to use as a reference
//' 
//' @return An integer vector with integers referring to elements in vector
//' \code{ref}.
//' 
//' @keywords internal
//' @noRd
Rcpp::IntegerVector LefkoUtils::refsort_str(Rcpp::CharacterVector vec, Rcpp::CharacterVector ref) {
  int vec_length = vec.length();
  int ref_length = ref.length();
  
  Rcpp::IntegerVector output(vec_length);
  
  for (int i = 0; i < vec_length; i++) {
    for (int j = 0; j < ref_length; j++) {
      if (LefkoUtils::stringcompare_hard(as<std::string>(vec[i]), as<std::string>(ref[j]))) output[i] = j + 1;
    }
  }
  
  return output;
}

//' Vectorize Matrix for Historical Mean Matrix Estimation
//' 
//' Function \code{flagrantcrap()} vectorizes core indices of matrices
//' input as list elements.
//' 
//' @name flagrantcrap
//' 
//' @param Xmat A matrix originally a part of a list object.
//' @param allindices A vector of indices to remove from the matrix
//' 
//' @return A column vector of specifically called elements from the input
//' matrix.
//' 
//' @keywords internal
//' @noRd
arma::vec LefkoUtils::flagrantcrap(arma::mat Xmat, arma::uvec allindices) {
  
  arma::vec newcol = Xmat.elem(allindices);
  
  return newcol;
}

//' Vectorize Matrix for Ahistorical Mean Matrix Estimation
//' 
//' Function \code{moreflagrantcrap()} vectorizes matrices input as list
//' elements.
//' 
//' @name moreflagrantcrap
//' 
//' @param Xmat A matrix originally a part of a list object.
//' 
//' @return A column vector of the input matrix.
//' 
//' @keywords internal
//' @noRd
arma::vec LefkoUtils::moreflagrantcrap(arma::mat Xmat) {
  
  arma::vec newcol = arma::vectorise(Xmat);
  
  return newcol;
}

//' Calculate Logarithms of Non-Zero Elements of Sparse Matrix
//' 
//' Function \code{spmat_log} finds the non-zero elements in a sparse matrix,
//' calculates their logs, and inserts them back into the matrix and returns it.
//' Based on code developed by Coatless Professor and posted by him on
//' StackOverflow.
//' 
//' @name spmat_log
//' 
//' @param B A sparse matrix. Note that this is assumed to be a population
//' projection matrix, meaning that all values are either 0 or positive.
//' 
//' @return A sparse matrix with non-zero values as logs of the elements in the
//' input matrix.
//' 
//' @keywords internal
//' @noRd
arma::sp_mat LefkoUtils::spmat_log(arma::sp_mat coremat)
{
  arma::sp_mat::const_iterator start = coremat.begin();
  arma::sp_mat::const_iterator end   = coremat.end();
  arma::sp_mat::const_iterator it = start; 
  
  int n = std::distance(start, end);
  
  if (n > 0) {
    arma::umat locs(2, n);
    arma::uvec temp(2);
    arma::vec vals(n);
    arma::vec logvals(n);
    locs.zeros();
    temp.zeros();
    vals.zeros();
    logvals.zeros();
    
    for(int i = 0; i < n; ++i) {
      temp(0) = it.row();
      temp(1) = it.col();
      locs.col(i) = temp;
      vals(i) = coremat(temp(0), temp(1));
      logvals(i) = log(vals(i));
      ++it; // increment
    }
    
    coremat = arma::sp_mat(locs, logvals, coremat.n_rows, coremat.n_cols);
  }
  
  return coremat;
}

//' Resize an IntegerVector
//' 
//' This function resizes an IntegerVector. It is based on code provided by Dirk
//' Eddelbuettel on StackExchange
//' (https://stackoverflow.com/questions/13782943/how-to-resize-a-numericvector).
//' 
//' @name shrink
//' 
//' @param x_int An IntegerVector.
//' 
//' @return A new IntegerVector exactly one element smaller than \code{x_int}.
//' 
//' @keywords internal
//' @noRd
Rcpp::IntegerVector LefkoUtils::shrink(Rcpp::IntegerVector x_int) {
  arma::ivec x = as<arma::ivec>(x_int);
  
  arma::ivec y = x;
  y.resize(x.size() - 1);
  
  Rcpp::IntegerVector y_out = as<IntegerVector>(wrap(y));
  
  return y_out;
}

