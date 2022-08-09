#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "LefkoUtils.h"

using namespace Rcpp;
using namespace arma;
using namespace LefkoUtils;


//' Function \code{lambda3()} is a generic function that returns the dominant
//' eigenvalue of a matrix, set of dominant eigenvalues of a set of matrices,
//' or set of dominant eigenvalues for a \code{lefkoMat} object. It can handle
//' large and sparse matrices supplied as \code{lefkoMat} objects or as
//' individual matrices, and can be used with large historical matrices, IPMs, 
//' age x stage matrices, as well as smaller ahistorical matrices.
//' 
//' @param mpm A lefkoMat object, a list of projection matrices, or a single
//' projection matrix.
//' @param sparse A string set to \code{"auto"} (the default), \code{"yes"}, or
//' \code{"no"}. If set to \code{"auto"}, then will determine whether to use
//' sparse matrix methods automatically.
//' 
//' @return The value returned depends on the class of the \code{mats} argument.
//' If a \code{lefkoMat} object is provided, then this function will return the
//' \code{labels} data frame with a new column named \code{lambda} showing the
//' dominant eigenvalues for each matrix. If a list of matrices is provided,
//' then this function will produce a numeric vector with the dominant
//' eigenvalues provided in order of matrix. If a single matrix is provided,
//' then this function will return the dominant eigenvalue of that matrix. Only
//' the largest real parts of the eigenvalues are returned.
//' 
//' @section Notes:
//' If \code{sparse = "auto"} (the default), then R will use sparse matrix
//' eigenanalysis if the matrices are both sparse (i.e, percentage of matrix
//' elements that are non-zero <= 50%) and have more than 100 rows.
//' 
//' @seealso \code{\link{slambda3}()}
//' 
//' @examples
//' # Lathyrus example
//' data(lathyrus)
//' 
//' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
//' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
//' repvector <- c(0, 0, 0, 0, 0, 1, 0)
//' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
//' matvector <- c(0, 0, 1, 1, 1, 1, 1)
//' immvector <- c(1, 1, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0)
//' indataset <- c(0, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
//' 
//' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector,
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
//'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
//'   propstatus = propvector)
//' 
//' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
//'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9,
//'   juvcol = "Seedling1988", sizeacol = "Volume88", repstracol = "FCODE88",
//'   fecacol = "Intactseed88", deadacol = "Dead1988",
//'   nonobsacol = "Dormant1988", stageassign = lathframe, stagesize = "sizea",
//'   censorcol = "Missing1988", censorkeep = NA, censor = TRUE)
//' 
//' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl", "mat"),
//'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep", "Sdl"),
//'   stage1 = c("Sd", "rep", "Sd", "rep", "npr", "npr", "Sd"),
//'   eststage3 = c(NA, NA, NA, NA, NA, NA, "mat"),
//'   eststage2 = c(NA, NA, NA, NA, NA, NA, "Sdl"),
//'   eststage1 = c(NA, NA, NA, NA, NA, NA, "NotAlive"),
//'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054, NA),
//'   type = c(1, 1, 1, 1, 3, 3, 1), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
//'   stageframe = lathframe, historical = TRUE)
//' 
//' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = "all", 
//'   stages = c("stage3", "stage2", "stage1"), supplement = lathsupp3,
//'   yearcol = "year2", indivcol = "individ")
//' 
//' ehrlen3mean <- lmean(ehrlen3)
//' lambda3(ehrlen3mean)
//' 
//' # Cypripedium example
//' data(cypdata)
//' 
//' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
//' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
//'   "XLg")
//' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
//' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
//' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
//' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
//' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
//' 
//' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
//'   propstatus = propvector, immstatus = immvector, indataset = indataset,
//'   binhalfwidth = binvec)
//' 
//' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
//'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
//'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
//'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
//'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
//'   NRasRep = TRUE)
//' 
//' # Here we use supplemental() to provide overwrite and reproductive info
//' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
//'     "XSm", "Sm", "SD", "P1"),
//'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
//'     "rep"),
//'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
//'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
//'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
//'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
//'   stageframe = cypframe_raw, historical = FALSE)
//' 
//' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
//'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
//'   size = c("size3added", "size2added"), supplement = cypsupp2r,
//'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
//' 
//' lambda3(cypmatrix2r)
//' 
//' @export lambda3
// [[Rcpp::export(lambda3)]]
RObject lambda3(RObject mpm, String sparse = "auto") {
  
  RObject output;
  
  int sparse_check = 0;
  
  if (stringcompare_simple(sparse, "aut", false)) {
    sparse_check = 2; // Auto decision
    
  } else if (stringcompare_simple(sparse, "y", false) || stringcompare_simple(sparse, "t", false)) {
    sparse_check = 1; // Forced sparse matrix
    
  } else if (stringcompare_simple(sparse, "n", false) || stringcompare_simple(sparse, "f", false)) {
    sparse_check = 0; // Forced dense matrix
  }
  
  if (is<List>(mpm)) {
    List mpm_ = as<List>(mpm);
    
    CharacterVector mpm_names = mpm_.attr("names");
    int no_mpm_names = mpm_names.length();
    
    bool A_check = false;
    for (int i = 0; i < no_mpm_names; i++) {
      if (stringcompare_simple(as<std::string>(mpm_names(i)), "A", false)) A_check = true;
    }
    
    bool labels_check = false;
    for (int i = 0; i < no_mpm_names; i++) {
      if (stringcompare_simple(as<std::string>(mpm_names(i)), "labels", false)) labels_check = true;
    }
    
    if (!A_check || !labels_check) {
      if (!Rf_isMatrix(mpm_[0])) {
        throw Rcpp::exception("Object mpm list structure is not recognized.", false);
      }
      
      int no_matrices = mpm_.length();
      
      if (sparse_check == 2) {
        arma::mat a1 = mpm_[0];
        int mat_rows = a1.n_rows;
        int mat_cols = a1.n_cols;
        int total_elems = mat_rows * mat_cols;
        
        arma::uvec nonzeros = find(a1);
        int no_nonzeros = nonzeros.n_elem;
        
        double density = static_cast<double>(no_nonzeros) / static_cast<double>(total_elems);
        
        if (density <= 0.5 && total_elems > 400) {
          sparse_check = 1;
        } else {
          sparse_check = 0;
        }
      }
      
      NumericVector lambda_prog (no_matrices);
      
      for (int i = 0; i < no_matrices; i++) {
        if (sparse_check == 0) {
          arma::mat Amat = as<arma::mat>(mpm_(i));
          
          arma::cx_vec Aeigval;
          arma::cx_mat Aeigvecl;
          arma::cx_mat Aeigvecr;
          
          eig_gen(Aeigval, Aeigvecl, Aeigvecr, Amat);
          
          arma::vec all_eigenvalues = real(Aeigval);
          double maxval = max(all_eigenvalues);
          
          arma::uvec max_elems = find(all_eigenvalues == maxval);
          if (max_elems.n_elem == 0) {
            throw Rcpp::exception("Eigenanalysis failed.", false);
          }
          
          arma::vec chosen_eigvals = all_eigenvalues.elem(max_elems);
          arma::uvec pos_max_eigvals = find(chosen_eigvals > 0);
          
          lambda_prog(i) = chosen_eigvals(pos_max_eigvals(0));
          
        } else {
          arma::sp_mat spAmat(as<arma::mat>(mpm_(i)));
          
          arma::cx_vec Aeigval;
          arma::cx_mat Aeigvecr;
          
          eigs_gen(Aeigval, Aeigvecr, spAmat, 1, "lr");
          
          arma::vec all_eigenvalues = real(Aeigval);
          double maxval = max(all_eigenvalues);
          
          arma::uvec max_elems = find(all_eigenvalues == maxval);
          if (max_elems.n_elem == 0) {
            throw Rcpp::exception("Eigenanalysis failed.", false);
          }
          arma::vec chosen_eigvals = all_eigenvalues.elem(max_elems);
          arma::uvec pos_max_eigvals = find(chosen_eigvals > 0);
          
          if (pos_max_eigvals.n_elem > 0) lambda_prog(i) = chosen_eigvals(pos_max_eigvals(0));
        }
      }
      output = lambda_prog;
      
    } else {
      List A_list = mpm_["A"];
      DataFrame labels = as<DataFrame>(mpm_["labels"]);
      
      if (!Rf_isMatrix(A_list[0])) throw Rcpp::exception("Object mpm does not appear to contain matrices", false);
      int no_matrices = A_list.length();
      
      if (sparse_check == 2) {
        arma::mat a1 = A_list[0];
        int mat_rows = a1.n_rows;
        int mat_cols = a1.n_cols;
        int total_elems = mat_rows * mat_cols;
        
        arma::uvec nonzeros = find(a1);
        int no_nonzeros = nonzeros.n_elem;
        
        double density = static_cast<double>(no_nonzeros) / static_cast<double>(total_elems);
        
        if (density <= 0.5 && total_elems > 400) {
          sparse_check = 1;
        } else {
          sparse_check = 0;
        }
      }
      
      NumericVector lambda_prog (no_matrices);
      
      for (int i = 0; i < no_matrices; i++) {
        if (sparse_check == 0) {
          arma::mat Amat = as<arma::mat>(A_list(i));
          
          arma::cx_vec Aeigval;
          arma::cx_mat Aeigvecl;
          arma::cx_mat Aeigvecr;
          
          eig_gen(Aeigval, Aeigvecl, Aeigvecr, Amat);
          
          arma::vec all_eigenvalues = real(Aeigval);
          double maxval = max(all_eigenvalues);
          
          arma::uvec max_elems = find(all_eigenvalues == maxval);
          if (max_elems.n_elem == 0) {
            throw Rcpp::exception("Eigenanalysis failed.", false);
          }
          
          arma::vec chosen_eigvals = all_eigenvalues.elem(max_elems);
          arma::uvec pos_max_eigvals = find(chosen_eigvals > 0);
          
          lambda_prog(i) = chosen_eigvals(pos_max_eigvals(0));
        } else {
          arma::sp_mat spAmat(as<arma::mat>(A_list(i)));
          
          arma::cx_vec Aeigval;
          arma::cx_mat Aeigvecr;
          
          eigs_gen(Aeigval, Aeigvecr, spAmat, 1, "lr");
          
          arma::vec all_eigenvalues = real(Aeigval);
          double maxval = max(all_eigenvalues);
          
          arma::uvec max_elems = find(all_eigenvalues == maxval);
          if (max_elems.n_elem == 0) {
            throw Rcpp::exception("Eigenanalysis failed.", false);
          }
          arma::vec chosen_eigvals = all_eigenvalues.elem(max_elems);
          arma::uvec pos_max_eigvals = find(chosen_eigvals > 0);
          
          if (pos_max_eigvals.n_elem > 0) lambda_prog(i) = chosen_eigvals(pos_max_eigvals(0));
        }
      }
      DataFrame new_out;
      
      CharacterVector l_pop = labels["pop"];
      CharacterVector l_patch = labels["patch"];
      
      int l_length = labels.length();
      
      if (l_length == 3) {
        NumericVector l_year2 = labels["year2"];
        
        new_out = DataFrame::create(_["pop"] = l_pop, _["patch"] = l_patch,
          _["year2"] = l_year2, _["lambda"] = lambda_prog);
      } else {
        new_out = DataFrame::create(_["pop"] = l_pop, _["patch"] = l_patch,
          _["lambda"] = lambda_prog);
      }
      output = new_out;
      
    }
    
  } else if(Rf_isMatrix(mpm)) {
    arma::mat mpm_ = as<arma::mat>(mpm);
    
    if (sparse_check == 2) {
      int mat_rows = mpm_.n_rows;
      int mat_cols = mpm_.n_cols;
      int total_elems = mat_rows * mat_cols;
      
      arma::uvec nonzeros = find(mpm_);
      int no_nonzeros = nonzeros.n_elem;
      
      double density = static_cast<double>(no_nonzeros) / static_cast<double>(total_elems);
      
      if (density <= 0.5 && total_elems > 400) {
        sparse_check = 1;
      } else {
        sparse_check = 0;
      }
    }
    
    NumericVector lambda_prog (1);
    
    if (sparse_check == 0) {
      arma::cx_vec Aeigval;
      arma::cx_mat Aeigvecl;
      arma::cx_mat Aeigvecr;
      
      eig_gen(Aeigval, Aeigvecl, Aeigvecr, mpm_);
      
      arma::vec all_eigenvalues = real(Aeigval);
      
      double maxval = max(all_eigenvalues);
      
      arma::uvec max_elems = find(all_eigenvalues == maxval);
      if (max_elems.n_elem == 0) {
        throw Rcpp::exception("Eigenanalysis failed.", false);
      }
      
      arma::vec chosen_eigvals = all_eigenvalues.elem(max_elems);
      arma::uvec pos_max_eigvals = find(chosen_eigvals > 0);
      
      lambda_prog(0) = chosen_eigvals(pos_max_eigvals(0));
    } else {
      arma::sp_mat spAmat(mpm_);
      
      arma::cx_vec Aeigval;
      arma::cx_mat Aeigvecr;
      
      eigs_gen(Aeigval, Aeigvecr, spAmat, 1, "lr");
      
      arma::vec all_eigenvalues = real(Aeigval);
      
      double maxval = max(all_eigenvalues);
      
      arma::uvec max_elems = find(all_eigenvalues == maxval);
      if (max_elems.n_elem == 0) {
        throw Rcpp::exception("Eigenanalysis failed.", false);
      }
      
      arma::vec chosen_eigvals = all_eigenvalues.elem(max_elems);
      arma::uvec pos_max_eigvals = find(chosen_eigvals > 0);
      
      if (pos_max_eigvals.n_elem > 0) lambda_prog(0) = chosen_eigvals(pos_max_eigvals(0));
    }
    
    output = lambda_prog;
  } else {
    throw Rcpp::exception("Object mpm does not appear to be an appropriate MPM.", false);
  }
  
  return output;
}
