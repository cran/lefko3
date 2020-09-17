#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Estimates Mean Population Projection Matrix, Using Summed U and F Matrices
//' 
//' This function estimates the mean population projection matrices, treating the
//' mean as arithmetic across space but either arithmetic or geometric across time.
//' It differs from \code{\link{turbogeodiesel}()} in that it estimates the \code{A} matrix
//' as a sum of the associated \code{U} and \code{F} matrices. Used to power the
//' \code{\link{lmean}()} function.
//' 
//' @param loy A data frame denoting the population, patch, and time step designation
//' of each matrix. Includes a total of 9 variables.
//' @param Umats A matrix with all U matrices turned into columns.
//' @param Fmats A matrix with all F matrices turned into columns.
//' @param geom Should the mean across time be geometric (1) or arithmetic (0)?
//' @param sparse Should 0s be ignored when some matrices include non-zero entries
//' in common elements?
//' @param numofpops Number of populations to be analyzed.
//' @param numofpatches Number of patches to be analyzed, where this number should
//' include a patch total across all populations.
//' @param numofyears Number of time steps to be analyzed.
//' 
//' @return A matrix with 3n columns, where n is the sum of the number of patches and
//' populations. Each pop/patch has its own set of three columns denoting survival,
//' fecundity, and the overall sum of the previous two columns.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat geodiesel(DataFrame loy, arma::mat Umats, arma::mat Fmats, int geom, int sparse, 
                    int numofpops, int numofpatches, int numofyears) {
  
  int elements = Umats.n_rows;
  int matrices = Umats.n_cols - 1;
  
  arma::vec poppatchc = loy["poppatchc"];
  arma::vec yearsinpatch = loy["yearsinpatch"];
  
  arma::vec indexmatU = Umats.col(matrices);
  arma::vec indexmatF = Fmats.col(matrices);
  
  arma::mat out(elements, ((numofpatches*3) + (numofpops * 3)));
  out.zeros();
  
  
  if (geom == 1) {
    if (sparse == 0) {
      // Geometric without sparse
      
      for (int i = 0; i < elements; i++) {
        
        if (indexmatU(i) != 0 || indexmatF(i) != 0) {
          
          for (int j = 0; j < matrices; j++) {
            //This section creates the sum of logs for each non-zero positive element
            
            if (Umats(i, j) != 0 && out(i, poppatchc(j)) != 1000) {
              out(i, poppatchc(j)) = out(i, poppatchc(j)) + (log(Umats(i, j)) / yearsinpatch(j));
            } else {
              out(i, poppatchc(j)) = 1000;
            }
            
            if (Fmats(i, j) != 0 && out(i, (numofpatches + poppatchc(j))) != 1000) {
              out(i, (numofpatches + poppatchc(j))) = out(i, (numofpatches + poppatchc(j))) + (log(Fmats(i, j)) / yearsinpatch(j));
            } else {
              out(i, (numofpatches + poppatchc(j))) = 1000;
            }
          }
          
          for (int k = 0; k < numofpatches; k++) {
            // This section takes e^logtotal for the U and F elements and creates the A elements,
            // and then creates the grand means
            
            // U elements
            if (out(i, k) != 1000) {
              out(i, k) = exp(out(i, k));
            } else {
              out(i, k) = 0;
            }
            
            // F elements
            if (out(i, (k + numofpatches)) != 1000) {
              out(i, (k + numofpatches)) = exp(out(i, (k + numofpatches)));
            } else {
              out(i, (k + numofpatches)) = 0;
            }
            
            // A elements
            out(i, (k + 2*numofpatches)) = out(i, k) + out(i, (k + numofpatches));
            
            // Grand mean elements
            out(i, (3*numofpatches)) = out(i, (3*numofpatches)) + (out(i, k) / numofpatches);
            out(i, (3*numofpatches + 1)) = out(i, (3*numofpatches + 1)) + (out(i, (numofpatches + k)) / numofpatches);
            out(i, (3*numofpatches + 2)) = out(i, (3*numofpatches)) + out(i, (3*numofpatches + 1));
          }
        }
      }
    } else {
      // Geometric with sparse
      arma::uvec totalvec(numofpatches*2);
      
      for (int i = 0; i < elements; i++) {
        
        if (indexmatU(i) != 0 || indexmatF(i) != 0) {
          
          for (int l = 0; l < numofpatches; l++) {
            //This section creates a vector to keep track of the dividend in mean calculations
            
            int properindex = numofyears * l;
            totalvec(l) = yearsinpatch(properindex);
            totalvec(l + numofpatches) = yearsinpatch(properindex);
          }
          
          for (int j = 0; j < matrices; j++) {
            //This section creates the sum of logs for each non-zero positive element
            
            if (Umats(i, j) != 0) {
              out(i, poppatchc(j)) = out(i, poppatchc(j)) + log(Umats(i, j));
            } else {
              totalvec(poppatchc(j)) = totalvec(poppatchc(j)) - 1;
            }
            
            if (Fmats(i, j) != 0) {
              out(i, (numofpatches + poppatchc(j))) = out(i, (numofpatches + poppatchc(j))) + log(Fmats(i, j));
            } else {
              totalvec(numofpatches + poppatchc(j)) = totalvec(numofpatches + poppatchc(j)) - 1;
            }
          }
          
          for (int k = 0; k < numofpatches; k++) {
            // This section takes e^logtotal for the U and F elements and creates the A elements,
            // and then creates the grand means
            
            // U elements
            if (totalvec(k) != 0) {
              out(i, k) = exp((out(i, k) / totalvec(k)));
            } else {
              out(i, k) = 0;
            }
            
            // F elements
            if (totalvec(k + numofpatches) != 0) {
              out(i, (k + numofpatches)) = exp((out(i, (k + numofpatches)) / totalvec(k + numofpatches)));
            } else {
              out(i, (k + numofpatches)) = 0;
            }
            
            // A elements
            out(i, (k + 2*numofpatches)) = out(i, k) + out(i, (k + numofpatches));
            
            // Grand mean elements
            out(i, (3*numofpatches)) = out(i, (3*numofpatches)) + (out(i, k) / numofpatches);
            out(i, (3*numofpatches + 1)) = out(i, (3*numofpatches + 1)) + (out(i, (numofpatches + k)) / numofpatches);
            out(i, (3*numofpatches + 2)) = out(i, (3*numofpatches)) + out(i, (3*numofpatches + 1));
          }
        }
      }
    }
  } else {
    if (sparse == 0) {
      // Arithmetic without sparse
      
      for (int i = 0; i < elements; i++) {
        
        if (indexmatU(i) != 0 || indexmatF(i) != 0) {
          
          for (int j = 0; j < matrices; j++) {
            //This section creates the means for U and F elements
            
            out(i, poppatchc(j)) = out(i, poppatchc(j)) + ((Umats(i, j)) / yearsinpatch(j));
            
            out(i, (numofpatches + poppatchc(j))) = out(i, (numofpatches + poppatchc(j))) + ((Fmats(i, j)) / yearsinpatch(j));
          }
          
          for (int k = 0; k < numofpatches; k++) {
            // This section creates the A elements,
            // and then creates the grand means
            
            // A elements
            out(i, (k + 2*numofpatches)) = out(i, k) + out(i, (k + numofpatches));
            
            // Grand mean elements
            out(i, (3*numofpatches)) = out(i, (3*numofpatches)) + (out(i, k) / numofpatches);
            out(i, (3*numofpatches + 1)) = out(i, (3*numofpatches + 1)) + (out(i, (numofpatches + k)) / numofpatches);
            out(i, (3*numofpatches + 2)) = out(i, (3*numofpatches)) + out(i, (3*numofpatches + 1));
          }
        }
      }
    } else {
      // Arithmetic with sparse
      
      arma::uvec totalvec(numofpatches*2);
      
      for (int i = 0; i < elements; i++) {
        
        if (indexmatU(i) != 0 || indexmatF(i) != 0) {
          
          for (int l = 0; l < numofpatches; l++) {
            //This section creates a vector to keep track of the dividend in mean calculations
            
            int properindex = numofyears * l;
            totalvec(l) = yearsinpatch(properindex);
            totalvec(l + numofpatches) = yearsinpatch(properindex);
          }
          
          for (int j = 0; j < matrices; j++) {
            //This section creates the sum of logs for each non-zero positive element
            
            if (Umats(i, j) != 0) {
              out(i, poppatchc(j)) = out(i, poppatchc(j)) + (Umats(i, j));
            } else {
              totalvec(poppatchc(j)) = totalvec(poppatchc(j)) - 1;
            }
            
            if (Fmats(i, j) != 0) {
              out(i, (numofpatches + poppatchc(j))) = out(i, (numofpatches + poppatchc(j))) + (Fmats(i, j));
            } else {
              totalvec(numofpatches + poppatchc(j)) = totalvec(numofpatches + poppatchc(j)) - 1;
            }
          }
          
          for (int k = 0; k < numofpatches; k++) {
            // This section takes sum for the U and F elements and creates the A elements,
            // and then creates the grand means
            
            // U elements
            if (totalvec(k) != 0) {
              out(i, k) = (out(i, k) / totalvec(k));
            } else {
              out(i, k) = 0;
            }
            
            // F elements
            if (totalvec(k + numofpatches) != 0) {
              out(i, (k + numofpatches)) = (out(i, (k + numofpatches)) / totalvec(k + numofpatches));
            } else {
              out(i, (k + numofpatches)) = 0;
            }
            
            // A elements
            out(i, (k + 2*numofpatches)) = out(i, k) + out(i, (k + numofpatches));
            
            // Grand mean elements
            out(i, (3*numofpatches)) = out(i, (3*numofpatches)) + (out(i, k) / numofpatches);
            out(i, (3*numofpatches + 1)) = out(i, (3*numofpatches + 1)) + (out(i, (numofpatches + k)) / numofpatches);
            out(i, (3*numofpatches + 2)) = out(i, (3*numofpatches)) + out(i, (3*numofpatches + 1));
          }
        }
      }
    }
  }
  
  return out;
}

//' Estimates Mean Population Projection Matrix
//' 
//' This function estimates the mean population projection matrices, treating the
//' mean as arithmetic across space but either arithmetic or geometric across time.
//' It differs from \code{\link{geodiesel}()} in that it estimates the \code{A} matrix
//' indepenently of the associated \code{U} and \code{F} matrices. Used to power the
//' \code{\link{lmean}()} function.
//' 
//' @param loy2c A data frame denoting the population, patch, and time step designation
//' of each matrix. Includes 9 variables.
//' @param Umats A matrix with all U matrices turned into columns.
//' @param Fmats A matrix with all F matrices turned into columns.
//' @param Amats A matrix with all A matrices turned into columns.
//' @param geom Should the mean across time be geometric (1) or arithmetic (0)?
//' @param sparse Should 0s be ignored when some matrices include non-zero entries
//' in common elements?
//' @param numofpops Number of populations to be analyzed.
//' @param numofpatches Number of patches to be analyzed, where this number should
//' include a patch total across all populations.
//' @param numofyears Number of time steps to be analyzed.
//' 
//' @return A matrix with 3n columns, where n is the sum of the number of patches and
//' populations. Each pop/patch has its own set of three columns denoting survival,
//' fecundity, and the overall projection.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat turbogeodiesel(DataFrame loy, arma::mat Umats, arma::mat Fmats, arma::mat Amats, int geom, 
                         int sparse, int numofpops, int numofpatches, int numofyears) {
  
  int elements = Umats.n_rows;
  int matrices = Umats.n_cols - 1;
  
  arma::vec poppatchc = loy["poppatchc"];
  arma::vec yearsinpatch = loy["yearsinpatch"];
  
  arma::vec indexmatU = Umats.col(matrices);
  arma::vec indexmatF = Fmats.col(matrices);
  
  arma::mat out(elements, ((numofpatches*3) + (numofpops * 3)));
  out.zeros();
  
  
  if (geom == 1) {
    if (sparse == 0) {
      // Geometric without sparse
      
      for (int i = 0; i < elements; i++) {
        
        if (indexmatU(i) != 0 || indexmatF(i) != 0) {
          
          for (int j = 0; j < matrices; j++) {
            //This section creates the sum of logs for each non-zero positive element
            
            if (Umats(i, j) != 0 && out(i, poppatchc(j)) != 1000) {
              out(i, poppatchc(j)) = out(i, poppatchc(j)) + (log(Umats(i, j)) / yearsinpatch(j));
            } else {
              out(i, poppatchc(j)) = 1000;
            }
            
            if (Fmats(i, j) != 0 && out(i, (numofpatches + poppatchc(j))) != 1000) {
              out(i, (numofpatches + poppatchc(j))) = out(i, (numofpatches + poppatchc(j))) + (log(Fmats(i, j)) / yearsinpatch(j));
            } else {
              out(i, (numofpatches + poppatchc(j))) = 1000;
            }
            
            if (Amats(i, j) != 0 && out(i, (2 * numofpatches + poppatchc(j))) != 1000) {
              out(i, (2 * numofpatches + poppatchc(j))) = out(i, (2 * numofpatches + poppatchc(j))) + (log(Amats(i, j)) / yearsinpatch(j));
            } else {
              out(i, (2 * numofpatches + poppatchc(j))) = 1000;
            }
          }
          
          for (int k = 0; k < numofpatches; k++) {
            // This section takes e^logtotal for the U and F elements and creates the A elements,
            // and then creates the grand means
            
            // U elements
            if (out(i, k) != 1000) {
              out(i, k) = exp(out(i, k));
            } else {
              out(i, k) = 0;
            }
            
            // F elements
            if (out(i, (k + numofpatches)) != 1000) {
              out(i, (k + numofpatches)) = exp(out(i, (k + numofpatches)));
            } else {
              out(i, (k + numofpatches)) = 0;
            }
            
            // A elements
            if (out(i, (k + 2 * numofpatches)) != 1000) {
              out(i, (k + 2 * numofpatches)) = exp(out(i, (k + 2 * numofpatches)));
            } else {
              out(i, (k + 2 * numofpatches)) = 0;
            }
            
            // Grand mean elements
            out(i, (3*numofpatches)) = out(i, (3*numofpatches)) + (out(i, k) / numofpatches);
            out(i, (3*numofpatches + 1)) = out(i, (3*numofpatches + 1)) + (out(i, (numofpatches + k)) / numofpatches);
            out(i, (3*numofpatches + 2)) = out(i, (3*numofpatches)) + out(i, (3*numofpatches + 1));
          }
        }
      }
    } else {
      // Geometric with sparse
      arma::uvec totalvec(numofpatches*3);
      
      for (int i = 0; i < elements; i++) {
        
        if (indexmatU(i) != 0 || indexmatF(i) != 0) {
          
          for (int l = 0; l < numofpatches; l++) {
            //This section creates a vector to keep track of the dividend in mean calculations
            
            int properindex = numofyears * l;
            totalvec(l) = yearsinpatch(properindex);
            totalvec(l + numofpatches) = yearsinpatch(properindex);
            totalvec(l + 2 * numofpatches) = yearsinpatch(properindex);
          }
          
          for (int j = 0; j < matrices; j++) {
            //This section creates the sum of logs for each non-zero positive element
            
            if (Umats(i, j) != 0) {
              out(i, poppatchc(j)) = out(i, poppatchc(j)) + log(Umats(i, j));
            } else {
              totalvec(poppatchc(j)) = totalvec(poppatchc(j)) - 1;
            }
            
            if (Fmats(i, j) != 0) {
              out(i, (numofpatches + poppatchc(j))) = out(i, (numofpatches + poppatchc(j))) + log(Fmats(i, j));
            } else {
              totalvec(numofpatches + poppatchc(j)) = totalvec(numofpatches + poppatchc(j)) - 1;
            }
            
            if (Amats(i, j) != 0) {
              out(i, (2 * numofpatches + poppatchc(j))) = out(i, (2 * numofpatches + poppatchc(j))) + log(Amats(i, j));
            } else {
              totalvec(2 * numofpatches + poppatchc(j)) = totalvec(2 * numofpatches + poppatchc(j)) - 1;
            }
          }
          
          for (int k = 0; k < numofpatches; k++) {
            // This section takes e^logtotal for the U and F elements and creates the A elements,
            // and then creates the grand means
            
            // U elements
            if (totalvec(k) != 0) {
              out(i, k) = exp((out(i, k) / totalvec(k)));
            } else {
              out(i, k) = 0;
            }
            
            // F elements
            if (totalvec(k + numofpatches) != 0) {
              out(i, (k + numofpatches)) = exp((out(i, (k + numofpatches)) / totalvec(k + numofpatches)));
            } else {
              out(i, (k + numofpatches)) = 0;
            }
            
            // A elements
            if (totalvec(k + 2 * numofpatches) != 0) {
              out(i, (k + 2 * numofpatches)) = exp((out(i, (k + 2 * numofpatches)) / totalvec(k + 2 * numofpatches)));
            } else {
              out(i, (k + 2 * numofpatches)) = 0;
            }
            
            // Grand mean elements
            out(i, (3*numofpatches)) = out(i, (3*numofpatches)) + (out(i, k) / numofpatches);
            out(i, (3*numofpatches + 1)) = out(i, (3*numofpatches + 1)) + (out(i, (numofpatches + k)) / numofpatches);
            out(i, (3*numofpatches + 2)) = out(i, (3*numofpatches)) + out(i, (3*numofpatches + 1));
          }
        }
      }
    }
  } else {
    if (sparse == 0) {
      // Arithmetic without sparse
      
      for (int i = 0; i < elements; i++) {
        
        if (indexmatU(i) != 0 || indexmatF(i) != 0) {
          
          for (int j = 0; j < matrices; j++) {
            //This section creates the means for U, F, and A elements
            
            out(i, poppatchc(j)) = out(i, poppatchc(j)) + ((Umats(i, j)) / yearsinpatch(j));
            
            out(i, (numofpatches + poppatchc(j))) = out(i, (numofpatches + poppatchc(j))) + ((Fmats(i, j)) / yearsinpatch(j));
            
            out(i, (2 * numofpatches + poppatchc(j))) = out(i, (2 * numofpatches + poppatchc(j))) + ((Amats(i, j)) / yearsinpatch(j));
          }
          
          for (int k = 0; k < numofpatches; k++) {
            // This section creates the grand means
            
            out(i, (3*numofpatches)) = out(i, (3*numofpatches)) + (out(i, k) / numofpatches);
            out(i, (3*numofpatches + 1)) = out(i, (3*numofpatches + 1)) + (out(i, (numofpatches + k)) / numofpatches);
            out(i, (3*numofpatches + 2)) = out(i, (3*numofpatches)) + out(i, (3*numofpatches + 1));
          }
        }
      }
    } else {
      // Arithmetic with sparse
      
      arma::uvec totalvec(numofpatches*3);
      
      for (int i = 0; i < elements; i++) {
        
        if (indexmatU(i) != 0 || indexmatF(i) != 0) {
          
          for (int l = 0; l < numofpatches; l++) {
            //This section creates a vector to keep track of the dividend in mean calculations
            
            int properindex = numofyears * l;
            totalvec(l) = yearsinpatch(properindex);
            totalvec(l + numofpatches) = yearsinpatch(properindex);
            totalvec(l + 2 * numofpatches) = yearsinpatch(properindex);
          }
          
          for (int j = 0; j < matrices; j++) {
            //This section creates the sum of logs for each non-zero positive element
            
            if (Umats(i, j) != 0) {
              out(i, poppatchc(j)) = out(i, poppatchc(j)) + (Umats(i, j));
            } else {
              totalvec(poppatchc(j)) = totalvec(poppatchc(j)) - 1;
            }
            
            if (Fmats(i, j) != 0) {
              out(i, (numofpatches + poppatchc(j))) = out(i, (numofpatches + poppatchc(j))) + (Fmats(i, j));
            } else {
              totalvec(numofpatches + poppatchc(j)) = totalvec(numofpatches + poppatchc(j)) - 1;
            }
            
            if (Amats(i, j) != 0) {
              out(i, (2 * numofpatches + poppatchc(j))) = out(i, (2 * numofpatches + poppatchc(j))) + (Amats(i, j));
            } else {
              totalvec(2 * numofpatches + poppatchc(j)) = totalvec(2 * numofpatches + poppatchc(j)) - 1;
            }
          }
          
          for (int k = 0; k < numofpatches; k++) {
            // This section takes sum for the U and F elements and creates the A elements,
            // and then creates the grand means
            
            // U elements
            if (totalvec(k) != 0) {
              out(i, k) = (out(i, k) / totalvec(k));
            } else {
              out(i, k) = 0;
            }
            
            // F elements
            if (totalvec(k + numofpatches) != 0) {
              out(i, (k + numofpatches)) = (out(i, (k + numofpatches)) / totalvec(k + numofpatches));
            } else {
              out(i, (k + numofpatches)) = 0;
            }
            
            // A elements
            if (totalvec(k + 2 * numofpatches) != 0) {
              out(i, (k + 2 * numofpatches)) = (out(i, (k + 2 * numofpatches)) / totalvec(k + 2 * numofpatches));
            } else {
              out(i, (k + 2 * numofpatches)) = 0;
            }
            
            // Grand mean elements
            out(i, (3*numofpatches)) = out(i, (3*numofpatches)) + (out(i, k) / numofpatches);
            out(i, (3*numofpatches + 1)) = out(i, (3*numofpatches + 1)) + (out(i, (numofpatches + k)) / numofpatches);
            out(i, (3*numofpatches + 2)) = out(i, (3*numofpatches)) + out(i, (3*numofpatches + 1));
          }
        }
      }
    }
  }
  
  return out;
}


