#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Vectorize Matrix for Historical Mean Matrix Estimation
//' 
//' Function \code{flagrantcrap()} vectorizes core indices of matrices
//' input as list elements.
//' 
//' @param Xmat A matrix originally a part of a list object.
//' @param allindices A vector of indices to remove from the matrix
//' 
//' @return A column vector of certain elements from the input matrix.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec flagrantcrap(arma::mat Xmat, arma::uvec allindices) {
  
  arma::vec newcol = Xmat.elem(allindices);
  
  return newcol;
}

//' Vectorize Matrix for Ahistorical Mean Matrix Estimation
//' 
//' Function \code{moreflagrantcrap()} vectorizes matrices input as list
//' elements.
//' 
//' @param Xmat A matrix originally a part of a list object.
//' 
//' @return A column vector of the input matrix.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec moreflagrantcrap(arma::mat Xmat) {
  
  arma::vec newcol = arma::vectorise(Xmat);
  
  return newcol;
}

//' Estimates Mean LefkoMat Object for Historical MPM
//' 
//' Function \code{turbogeodiesel()} estimates mean historical population
//' projection matrices, treating the mean as element-wise arithmetic.
//' 
//' @param loy A data frame denoting the population, patch, and time step
//' designation of each matrix. Includes a total of 9 variables.
//' @param Umats A matrix with all U matrices turned into columns.
//' @param Fmats A matrix with all F matrices turned into columns.
//' @param stages This is the core stageframe held by \code{mats}, equivalent to
//' \code{ahstages}.
//' @param hstages This is the \code{hstages} object held by \code{mats}.
//' @param patchmats A logical value stating whether to estimate patch-level
//' means.
//' @param popmats A logical value stating whether to estimate population-level
//' means.
//' 
//' @return A list using the basic blueprint of a lefkoMat object.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List turbogeodiesel(DataFrame loy, List Umats, List Fmats, DataFrame stages, 
  DataFrame hstages, bool patchmats, bool popmats) {
  
  StringVector pops = loy["pop"];
  arma::uvec pop_num = loy["popc"];
  StringVector patches = loy["patch"];
  arma::uvec year2 = loy["year2"];
  arma::uvec poppatchc = loy["poppatchc"];
  arma::uvec patchesinpop = loy["patchesinpop"];
  arma::uvec yearsinpatch = loy["yearsinpatch"];
  arma::uvec uniquepops = unique(pop_num);
  arma::uvec uniquepoppatches = unique(poppatchc);
  int loydim = pops.length();
  int numofpops = uniquepops.n_elem;
  int numofpatches = uniquepoppatches.n_elem;
  
  if (numofpatches == 1) popmats = 0;
  
  StringVector poporderlong(loydim);
  arma::uvec poporderlong_num(loydim);
  StringVector patchorderlong(loydim);
  arma::uvec annmatriceslong(loydim);
  arma::uvec meanassign(loydim);
  poporderlong_num.zeros();
  annmatriceslong.zeros();
  meanassign.zeros();
  
  pop_num = pop_num + 1;
  poppatchc = poppatchc + 1;
  
  poporderlong(0) = pops(0);
  poporderlong_num(0) = pop_num(0);
  patchorderlong(0) = patches(0);
  annmatriceslong(0) = 1;
  meanassign(0) = 1;
  
  int counter {0};
  
  StringVector uniquepops_str(numofpops);
  uniquepops_str(0) = pops(0);
  int popcounter {0};
  
  // In this beginning bit, we assess how many mean matrices we will need, and the overall order of means
  if (loydim > 1) {
    for (int i = 1; i < loydim; i++) {
      
      if (poppatchc(i) != poppatchc(i-1)) {
        
        counter++;
        poporderlong(counter) = pops(i);
        poporderlong_num(counter) = pop_num(i);
        patchorderlong(counter) = patches(i);
        annmatriceslong(counter) = 1;
        meanassign(i) = meanassign(i-1) + 1;
        
        if (pop_num(i) != pop_num(i-1)) {
          popcounter += 1;
          uniquepops_str(popcounter) = pops(i);
        }
        
      } else {
        
        annmatriceslong(counter) = annmatriceslong(counter) + 1;
        meanassign(i) = meanassign(i-1);
        
      }
    }
  }
  
  arma::uvec toestimate = find(poporderlong_num);
  int popcount = toestimate.n_elem;
  
  int totalmatrices = toestimate.n_elem + numofpops;
  
  if (patchmats == 1 && popmats == 0) {
    totalmatrices = toestimate.n_elem;
  } else if (patchmats == 0 && popmats == 1) {
    totalmatrices = numofpops;
  }
  
  arma::uvec poporder = poporderlong_num.elem(toestimate);
  arma::uvec patchorder = poppatchc.elem(toestimate);
  arma::uvec annmatrices = annmatriceslong.elem(toestimate);
  
  StringVector poporder_str(popcount);
  StringVector patchorder_str(popcount);
  
  for (int i = 0; i < popcount; i++) {
    poporder_str(i) = pops(toestimate(i));
    patchorder_str(i) = patches(toestimate(i));
  }
  
  // This next chunk predicts which elements will be targeted for arithmetic mean estimation
  arma::uvec astages = stages["stage_id"];
  int numstages = astages.n_elem;
  
  arma::uvec hstage3in = hstages["stage_id_2"];
  arma::uvec hstage2nin = hstages["stage_id_1"];
  int numhstages = hstage3in.n_elem;
  
  int predictedsize = numstages * numstages * numstages;
  
  arma::uvec hsindexl(predictedsize);
  hsindexl.zeros();
  
  counter = 0;
  
  for (int i1 = 0; i1 < numhstages; i1++) {
    for (int i2 = 0; i2 < numhstages; i2++) {
      if (hstage3in(i1) == hstage2nin(i2)) {
        
        hsindexl(counter) = (i1 * numhstages) + i2;
        
        counter++;
      }
    }
  }
  
  arma::uvec hsgood = find(hsindexl);
  arma::uvec hsindex = hsindexl.elem(hsgood);
  arma::uvec zerovec(1);
  zerovec.zeros();
  arma::uvec allindices = join_cols(zerovec, hsindex);
  
  // Now we build a U and F matrices of element-wise arithmetic means, where
  // each column corresponds to the predicted non-zero elements of each mean
  // matrix, and each matrix is presented as a column vector within the 
  // overall matrix. The A matrix is the sum of U and F.
  
  int core_elem = counter;
  
  arma::mat umatvec(core_elem, totalmatrices);
  arma::mat fmatvec(core_elem, totalmatrices);
  umatvec.zeros();
  fmatvec.zeros();
  
  int patchchoice {0};
  int popchoice {0};
  
  pop_num = pop_num - 1;
  poppatchc = poppatchc - 1;
  
  for (int i = 0; i < loydim; i++) {
    
    if (patchmats == 1) {
      
      patchchoice = poppatchc(i);
      
      umatvec.col(patchchoice) = umatvec.col(patchchoice) + (flagrantcrap(Umats[i], allindices) / yearsinpatch(i));
      fmatvec.col(patchchoice) = fmatvec.col(patchchoice) + (flagrantcrap(Fmats[i], allindices) / yearsinpatch(i));
      
    }
    
    if (popmats == 1) {
      if (patchmats == 1) {
        
        popchoice = numofpatches + pop_num(i);
        
      } else {
        
        popchoice = pop_num(i);
        
      }
      
      umatvec.col(popchoice) = umatvec.col(popchoice) + (flagrantcrap(Umats[i], allindices) / (yearsinpatch(i) * patchesinpop(i)));
      fmatvec.col(popchoice) = fmatvec.col(popchoice) + (flagrantcrap(Fmats[i], allindices) / (yearsinpatch(i) * patchesinpop(i)));
      
    }
  }
  
  arma::mat amatvec = umatvec + fmatvec;
  
  // Here we create the cheat sheet algorithm
  int cheatsheetlength {1};
  if (numofpatches > 1) cheatsheetlength = numofpops + numofpatches;
  StringVector poporder_redone(cheatsheetlength);
  StringVector patchorder_redone(cheatsheetlength);
  
  if (numofpatches > 1) {
    for (int i = 0; i < numofpatches; i++) {
      poporder_redone(i) = poporderlong(i);
      patchorder_redone(i) = patchorderlong(i);
    }
    
    for (int i = 0; i < numofpops; i++) {
      poporder_redone(i+numofpatches) = uniquepops_str(i);
      patchorder_redone(i+numofpatches) = "0";
    }
  } else {
    poporder_redone(0) = poporderlong(0);
    patchorder_redone(0) = patchorderlong(0);
  }
  
  DataFrame cheatsheet = DataFrame::create(Named("pop") = poporder_redone, _["patch"] = patchorder_redone);
  
  // Now we will create the main list objects holding the matrices
  arma::mat umat_base(numhstages, numhstages);
  arma::mat fmat_base(numhstages, numhstages);
  arma::mat amat_base(numhstages, numhstages);
  
  umat_base.zeros();
  fmat_base.zeros();
  amat_base.zeros();
  
  umat_base.elem(allindices) = umatvec.col(0);
  fmat_base.elem(allindices) = fmatvec.col(0);
  amat_base.elem(allindices) = amatvec.col(0);
  
  List U = List::create(umat_base);
  List F = List::create(fmat_base);
  List A = List::create(amat_base);
  
  if (totalmatrices > 1) {
    for (int i = 1; i < totalmatrices; i++) {
      umat_base.zeros();
      fmat_base.zeros();
      amat_base.zeros();
      
      umat_base.elem(allindices) = umatvec.col(i);
      fmat_base.elem(allindices) = fmatvec.col(i);
      amat_base.elem(allindices) = amatvec.col(i);
      
      U.push_back(umat_base);
      F.push_back(fmat_base);
      A.push_back(amat_base);
    }
  }
  
  // Matrix QC output
  
  arma::uvec utrans = find(umatvec);
  arma::uvec ftrans = find(fmatvec);
  int totalutrans = utrans.n_elem;
  int totalftrans = ftrans.n_elem;
  
  arma::vec matrixqc(3);
  matrixqc(0) = totalutrans; // summed number of u transitions
  matrixqc(1) = totalftrans; // summed number of f transitions
  matrixqc(2) = totalmatrices;
  
  
  // Final output
  
  List output = List::create(Named("A") = A, _["U"] = U, _["F"] = F,
    _["hstages"] = hstages, _["ahstages"] = stages, _["labels"] = cheatsheet,
    _["matrixqc"] = matrixqc);
  
  return output;
}

//' Estimates Mean LefkoMat Object for Ahistorical MPM
//' 
//' Function \code{geodiesel()} estimates mean ahistorical population projection
//' matrices, treating the mean as element-wise arithmetic.
//' 
//' @param loy A data frame denoting the population, patch, and time step
//' designation of each matrix. Includes a total of 9 variables.
//' @param Umats A matrix with all U matrices turned into columns.
//' @param Fmats A matrix with all F matrices turned into columns.
//' @param stages This is the core stageframe held by \code{mats}, equivalent to
//' \code{ahstages}.
//' @param patchmats A logical value stating whether to estimate patch-level
//' means.
//' @param popmats A logical value stating whether to estimate population-level
//' means.
//' 
//' @return A list using the basic blueprint of a LefkoMat object.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List geodiesel(DataFrame loy, List Umats, List Fmats, DataFrame stages,
  bool patchmats, bool popmats) {
  
  StringVector pops = loy["pop"];
  arma::uvec pop_num = loy["popc"];
  StringVector patches = loy["patch"];
  arma::uvec year2 = loy["year2"];
  arma::uvec poppatchc = loy["poppatchc"];
  arma::uvec patchesinpop = loy["patchesinpop"];
  arma::uvec yearsinpatch = loy["yearsinpatch"];
  arma::uvec uniquepops = unique(pop_num);
  arma::uvec uniquepoppatches = unique(poppatchc);
  int loydim = pops.length();
  int numofpops = uniquepops.n_elem;
  int numofpatches = uniquepoppatches.n_elem;
  
  if (numofpatches == 1) popmats = 0;
  
  StringVector poporderlong(loydim);
  arma::uvec poporderlong_num(loydim);
  StringVector patchorderlong(loydim);
  arma::uvec annmatriceslong(loydim);
  arma::uvec meanassign(loydim);
  poporderlong_num.zeros();
  annmatriceslong.zeros();
  meanassign.zeros();
  
  pop_num = pop_num + 1;
  poppatchc = poppatchc + 1;
  
  poporderlong(0) = pops(0);
  poporderlong_num(0) = pop_num(0);
  patchorderlong(0) = patches(0);
  annmatriceslong(0) = 1;
  meanassign(0) = 1;
  
  int counter {0};
  
  StringVector uniquepops_str(numofpops);
  uniquepops_str(0) = pops(0);
  int popcounter {0};
  
  // In this beginning bit, we assess how many mean matrices we will need, and the overall order of means
  if (loydim > 1) {
    for (int i = 1; i < loydim; i++) {
      
      if (poppatchc(i) != poppatchc(i-1)) {
        
        counter++;
        poporderlong(counter) = pops(i);
        poporderlong_num(counter) = pop_num(i);
        patchorderlong(counter) = patches(i);
        annmatriceslong(counter) = 1;
        meanassign(i) = meanassign(i-1) + 1;
        
        if (pop_num(i) != pop_num(i-1)) {
          popcounter += 1;
          uniquepops_str(popcounter) = pops(i);
        }
        
      } else {
        
        annmatriceslong(counter) = annmatriceslong(counter) + 1;
        meanassign(i) = meanassign(i-1);
        
      }
    }
  }
  
  arma::uvec toestimate = find(poporderlong_num);
  int popcount = toestimate.n_elem;
  
  int totalmatrices = toestimate.n_elem + numofpops;
  
  if (patchmats == 1 && popmats == 0) {
    totalmatrices = toestimate.n_elem;
  } else if (patchmats == 0 && popmats == 1) {
    totalmatrices = numofpops;
  }
  
  arma::uvec poporder = poporderlong_num.elem(toestimate);
  arma::uvec patchorder = poppatchc.elem(toestimate);
  arma::uvec annmatrices = annmatriceslong.elem(toestimate);
  
  StringVector poporder_str(popcount);
  StringVector patchorder_str(popcount);
  
  for (int i = 0; i < popcount; i++) {
    poporder_str(i) = pops(toestimate(i));
    patchorder_str(i) = patches(toestimate(i));
  }
  
  // This next chunk predicts which elements will be targeted for arithmetic mean estimation
  arma::uvec astages = stages["stage_id"];
  int numstages = astages.n_elem;
  
  // Now we build a U and F matrices of element-wise arithmetic means, where
  // each column corresponds to the predicted non-zero elements of each mean
  // matrix, and each matrix is presented as a column vector within the 
  // overall matrix. The A matrix is the sum of U and F.
  
  int core_elem = numstages * numstages;
  
  arma::mat umatvec(core_elem, totalmatrices);
  arma::mat fmatvec(core_elem, totalmatrices);
  umatvec.zeros();
  fmatvec.zeros();
  
  int patchchoice {0};
  int popchoice {0};
  
  pop_num = pop_num - 1;
  poppatchc = poppatchc - 1;
  
  for (int i = 0; i < loydim; i ++) {
    
    if (patchmats == 1) {
      
      patchchoice = poppatchc(i);
      
      umatvec.col(patchchoice) = umatvec.col(patchchoice) + (moreflagrantcrap(Umats[i]) / yearsinpatch(i));
      fmatvec.col(patchchoice) = fmatvec.col(patchchoice) + (moreflagrantcrap(Fmats[i]) / yearsinpatch(i));
      
    }
    
    if (popmats == 1) {
      if (patchmats == 1) {
        
        popchoice = numofpatches + pop_num(i);
        
      } else {
        
        popchoice = pop_num(i);
        
      }
      
      umatvec.col(popchoice) = umatvec.col(popchoice) + (moreflagrantcrap(Umats[i]) / (yearsinpatch(i) * patchesinpop(i)));
      fmatvec.col(popchoice) = fmatvec.col(popchoice) + (moreflagrantcrap(Fmats[i]) / (yearsinpatch(i) * patchesinpop(i)));
      
    }
  }
  
  arma::mat amatvec = umatvec + fmatvec;
  
  // Here we create the cheat sheet algorithm
  int cheatsheetlength {1};
  if (numofpatches > 1) cheatsheetlength = numofpops + numofpatches;
  StringVector poporder_redone(cheatsheetlength);
  StringVector patchorder_redone(cheatsheetlength);
  
  if (numofpatches > 1) {
    for (int i = 0; i < numofpatches; i++) {
      poporder_redone(i) = poporderlong(i);
      patchorder_redone(i) = patchorderlong(i);
    }
    
    for (int i = 0; i < numofpops; i++) {
      poporder_redone(i+numofpatches) = uniquepops_str(i);
      patchorder_redone(i+numofpatches) = "0";
    }
  } else {
    poporder_redone(0) = poporderlong(0);
    patchorder_redone(0) = patchorderlong(0);
  }
  
  DataFrame cheatsheet = DataFrame::create(Named("pop") = poporder_redone, 
    _["patch"] = patchorder_redone);
  
  // Now we will create the main list objects holding the matrices
  
  arma::mat umat_base = umatvec.col(0);
  umat_base.reshape(numstages, numstages);
  
  arma::mat fmat_base = fmatvec.col(0);
  fmat_base.reshape(numstages, numstages);
  
  arma::mat amat_base = amatvec.col(0);
  amat_base.reshape(numstages, numstages);
  
  List U = List::create(umat_base);
  List F = List::create(fmat_base);
  List A = List::create(amat_base);
  
  if (totalmatrices > 1) {
    for (int i = 1; i < totalmatrices; i++) {
      umat_base.zeros();
      fmat_base.zeros();
      amat_base.zeros();
      
      umat_base = umatvec.col(i);
      fmat_base = fmatvec.col(i);
      amat_base = amatvec.col(i);
      
      umat_base.reshape(numstages, numstages);
      fmat_base.reshape(numstages, numstages);
      amat_base.reshape(numstages, numstages);
      
      U.push_back(umat_base);
      F.push_back(fmat_base);
      A.push_back(amat_base);
    }
  }
  
  // Matrix QC output
  
  arma::uvec utrans = find(umatvec);
  arma::uvec ftrans = find(fmatvec);
  int totalutrans = utrans.n_elem;
  int totalftrans = ftrans.n_elem;
  
  arma::vec matrixqc(3);
  matrixqc(0) = totalutrans; // summed number of u transitions
  matrixqc(1) = totalftrans; // summed number of f transitions
  matrixqc(2) = totalmatrices;
  
  
  // Final output
  
  List output = List::create(Named("A") = A, _["U"] = U, _["F"] = F, 
    _["hstages"] = NULL, _["ahstages"] = stages, _["labels"] = cheatsheet, 
    _["matrixqc"] = matrixqc);
  
  return output;
}

//' Complete Full Eigen Analysis of a Single Dense Matrix
//' 
//' \code{decomp3()} returns all eigenvalues, right eigenvectors, and left
//' eigenvectors estimated for a matrix by the \code{eig_gen}() function
//' in the C++ Armadillo library. Works with dense matrices.
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//'
//' @return This function returns all estimated eigenvalues, right
//' eigenvectors, and left eigenvectors of a single matrix. This output is
//' provided as a list with three parts, named appropriately.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List decomp3(arma::mat Amat) {
  
  arma::cx_vec Aeigval;
  arma::cx_mat Aeigvecl;
  arma::cx_mat Aeigvecr;
  
  eig_gen(Aeigval, Aeigvecl, Aeigvecr, Amat);
  
  List output = List::create(Named("eigenvalues") = Aeigval, _["left_eigenvectors"] = Aeigvecl,
                             _["right_eigenvectors"] = Aeigvecr);
  
  return output;
}

//' Complete Full Eigen Analysis of a Single Sparse Matrix
//' 
//' \code{decomp3sp()} returns all eigenvalues, right eigenvectors, and left
//' eigenvectors estimated for a matrix by the \code{eigs_gen}() function
//' in the C++ Armadillo library. Works with sparse matrices.
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//'
//' @return This function returns all estimated eigenvalues, right
//' eigenvectors, and left eigenvectors of a single matrix. This output is
//' provided as a list with three parts, named appropriately.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List decomp3sp(arma::mat Amat) {
  
  arma::sp_mat spAmat(Amat);
  arma::sp_mat t_spAmat = spAmat.t();
  
  arma::cx_vec Aeigval;
  arma::cx_vec Aeigvall;
  arma::cx_mat Aeigvecl;
  arma::cx_mat Aeigvecr;
  
  eigs_gen(Aeigval, Aeigvecr, spAmat, 1);
  
  eigs_gen(Aeigvall, Aeigvecl, t_spAmat, 1);
  
  List output = List::create(Named("eigenvalues") = Aeigval, _["left_eigenvectors"] = Aeigvecl,
                             _["right_eigenvectors"] = Aeigvecr);
  
  return output;
}

//' Estimate Deterministic Population Growth Rate of a Dense Matrix
//' 
//' \code{lambda3matrix()} returns the dominant eigenvalue of a single
//' dense projection matrix.
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//'
//' @return This function returns the dominant eigenvalue of the matrix. This
//' is given as the largest real part of all eigenvalues estimated via the 
//' \code{eig_gen}() function in the C++ Armadillo library.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double lambda3matrix(arma::mat Amat) {
  
  List eigenstuff = decomp3(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  double lambda = max(realeigenvals);
  
  return lambda;
}

//' Estimate Deterministic Population Growth Rate of a Sparse Matrix
//' 
//' \code{lambda3matrixsp()} returns the dominant eigenvalue of a single
//' sparse projection matrix. This function can handle large and sparse 
//' matrices, and so can be used with large historical matrices, IPMs, 
//' age x stage matrices, as well as smaller ahistorical matrices.
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//'
//' @return This function returns the dominant eigenvalue of the matrix. This
//' is given as the largest real part of all eigenvalues estimated via the 
//' \code{eigs_gen}() function in the C++ Armadillo library.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
double lambda3matrixsp(arma::mat Amat) {
  
  List eigenstuff = decomp3sp(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  double lambda = max(realeigenvals);
  
  return lambda;
}

//' Estimate Stable Stage Distribution for a Dense Population Matrix
//' 
//' \code{ss3matrix()} returns the stable stage distribution for a 
//' dense population matrix.
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' 
//' @return This function returns the stable stage distribution corresponding to
//' the input matrix. The stable stage distribution is given as the right 
//' eigenvector associated with largest real part of the eigenvalues estimated 
//' for the matrix via the \code{eig_gen}() function in the C++ Armadillo 
//' library, divided by the sum of the associated right eigenvector. 
//' 
//' @seealso \code{\link{stablestage3}()}
//' @seealso \code{\link{stablestage3.lefkoMat}()}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec ss3matrix(arma::mat Amat) {
  
  List eigenstuff = decomp3(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  int lambda1 = realeigenvals.index_max();
  
  cx_mat rightvec = eigenstuff["right_eigenvectors"];
  arma::vec realrightvec = real(rightvec.col(lambda1));
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  
  arma::vec wcorr (rvel);
  
  for (int i = 0; i < rvel; i++) {
    wcorr(i) = realrightvec(i) / rvsum;
  }
  
  return wcorr;
}

//' Estimate Stable Stage Distribution for a Sparse Population Matrix
//' 
//' \code{ss3matrixsp()} returns the stable stage distribution for a 
//' sparse population matrix. This function can handle large and sparse 
//' matrices, and so can be used with large historical matrices, IPMs, 
//' age x stage matrices, as well as smaller ahistorical matrices.
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' 
//' @return This function returns the stable stage distribution corresponding to
//' the input matrix. The stable stage distribution is given as the right
//' eigenvector associated with largest real part of the eigenvalues estimated
//' for the matrix via the \code{eigs_gen}() function in the C++ Armadillo
//' library, divided by the sum of the associated right eigenvector. 
//' 
//' @seealso \code{\link{stablestage3}()}
//' @seealso \code{\link{stablestage3.lefkoMat}()}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec ss3matrixsp(arma::mat Amat) {
  
  List eigenstuff = decomp3sp(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  int lambda1 = realeigenvals.index_max();
  
  cx_mat rightvec = eigenstuff["right_eigenvectors"];
  arma::vec realrightvec = real(rightvec.col(lambda1));
  realrightvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-14 with 0
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  
  arma::vec wcorr (rvel);
  
  for (int i = 0; i < rvel; i++) {
    wcorr(i) = realrightvec(i) / rvsum;
  }
  
  return wcorr;
}

//' Estimate Reproductive Value for a Dense Population Matrix
//' 
//' \code{rv3matrix()} returns the reproductive values for stages in a
//' dense population matrix. The function provides standard reproductive
//' values, meaning that the overall reproductive values of basic life
//' history stages in a historical matrix are not provided (the
//' \code{\link{repvalue3.lefkoMat}()} function estimates these on the basis
//' of stage description information provided in the \code{lefkoMat} object
//' used as input in that function).
//' 
//' @param Amat A population projection matrix.
//' 
//' @return This function returns a vector characterizing the
//' reproductive values for stages of a population projection matrix. This is
//' given as the left eigenvector associated with largest real part of the
//' dominant eigenvalue estimated via the \code{eig_gen}() function in the C++
//' Armadillo library, divided by the first non-zero element of the left
//' eigenvector. 
//' 
//' @seealso \code{\link{repvalue3}()}
//' @seealso \code{\link{repvalue3.lefkoMat}()}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec rv3matrix(arma::mat Amat) {
  
  List eigenstuff = decomp3(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  int lambda1 = realeigenvals.index_max();
  
  cx_mat leftvec = eigenstuff["left_eigenvectors"];
  arma::vec realleftvec = real(leftvec.col(lambda1));
  arma::vec rlvabs = abs(realleftvec);
  
  // This section identifies the first non-zero element of the reproductive value vector
  arma::uvec rlvabsalt = find(rlvabs);
  int rlvminelem = rlvabsalt(0);
  double rlvmin = realleftvec(rlvminelem);
  
  int lvel = realleftvec.n_elem;
  
  arma::vec vcorr (lvel);
  
  for (int i = 0; i < lvel; i++) {
    vcorr(i) = realleftvec(i) / rlvmin;
  }
  
  return vcorr;
}

//' Estimate Reproductive Value for a Sparse Population Matrix
//' 
//' \code{rv3matrixsp()} returns the reproductive values for stages in a 
//' sparse population matrix. The function provides standard reproductive 
//' values, meaning that the overall reproductive values of basic life 
//' history stages in a historical matrix are not provided (the 
//' \code{\link{repvalue3.lefkoMat}()} function estimates these on the basis 
//' of stage description information provided in the \code{lefkoMat} object 
//' used as input in that function). This function can handle large and 
//' sparse matrices, and so can be used with large historical matrices, IPMs, 
//' age x stage matrices, as well as smaller ahistorical 
//' matrices.
//' 
//' @param Amat A population projection matrix.
//' 
//' @return This function returns a vector characterizing the
//' reproductive values for stages of a population projection matrix. This is
//' given as the left eigenvector associated with largest real part of the
//' dominant eigenvalue estimated via the \code{eigs_gen}() function in the C++
//' Armadillo library, divided by the first non-zero element of the left
//' eigenvector. 
//' 
//' @seealso \code{\link{repvalue3}()}
//' @seealso \code{\link{repvalue3.lefkoMat}()}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec rv3matrixsp(arma::mat Amat) {
  
  List eigenstuff = decomp3sp(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  int lambda1 = realeigenvals.index_max();
  
  cx_mat leftvec = eigenstuff["left_eigenvectors"];
  arma::vec realleftvec = real(leftvec.col(lambda1));
  realleftvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-14 with 0
  arma::vec rlvabs = abs(realleftvec);
  
  // This section identifies the first non-zero element of the reproductive value vector
  arma::uvec rlvabsalt = find(rlvabs);
  int rlvminelem = rlvabsalt(0);
  double rlvmin = realleftvec(rlvminelem);
  
  int lvel = realleftvec.n_elem;
  
  arma::vec vcorr (lvel);
  
  for (int i = 0; i < lvel; i++) {
    vcorr(i) = realleftvec(i) / rlvmin;
  }
  
  return vcorr;
}

//' Estimate Sensitivities for a Dense Population Matrix
//' 
//' \code{sens3matrix()} returns the sensitivity of lambda with respect
//' to each element in a dense matrix. This is accomplished via the
//' \code{eig_gen}() function in the C++ Armadillo library.
//' 
//' @param Amat A population projection matrix.
//' 
//' @return This function returns a matrix of sensitivities. 
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat sens3matrix(arma::mat Amat) {
  
  List eigenstuff = decomp3(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  int lambda1 = realeigenvals.index_max();
  
  // This is the w vector
  cx_mat rightvec = eigenstuff["right_eigenvectors"];
  arma::vec realrightvec = real(rightvec.col(lambda1));
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  
  // This is the v vector
  cx_mat leftvec = eigenstuff["left_eigenvectors"];
  arma::vec realleftvec = real(leftvec.col(lambda1));
  arma::vec rlvabs = abs(realleftvec);
  
  arma::uvec rlvabsalt = find(rlvabs);
  int rlvminelem = rlvabsalt(0);
  double rlvmin = realleftvec(rlvminelem);
  
  arma::vec wcorr (rvel);
  arma::vec vcorr (rvel);
  arma::vec vwprod (rvel);
  arma::mat smat (rvel, rvel);
  smat.zeros();
  
  // This loop and the following line create the scalar product vw
  for (int i = 0; i < rvel; i++) {
    wcorr(i) = realrightvec(i) / rvsum;
    vcorr(i) = realleftvec(i) / rlvmin;
    
    vwprod(i) = wcorr(i) * vcorr(i);
  }
  
  double vwscalar = sum(vwprod);
  
  // This loop populates the sensitivity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      
      smat(i, j) = vcorr(i) * wcorr(j) / vwscalar;
    }
  }
  
  return smat;
}

//' Estimate Sensitivities for a Sparse Population Matrix
//' 
//' \code{sens3matrixsp()} returns the sensitivity of lambda with respect
//' to each element in a sparse matrix. This is accomplished via the
//' \code{eig_gen}() function in the C++ Armadillo library.
//' 
//' @param Amat A population projection matrix.
//' 
//' @return This function returns a matrix of sensitivities. 
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat sens3matrixsp(arma::mat Amat) {
  
  List eigenstuff = decomp3sp(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  int lambda1 = realeigenvals.index_max();
  
  // This is the w vector
  cx_mat rightvec = eigenstuff["right_eigenvectors"];
  arma::vec realrightvec = real(rightvec.col(lambda1));
  realrightvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-14 with 0
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  
  // This is the v vector
  cx_mat leftvec = eigenstuff["left_eigenvectors"];
  arma::vec realleftvec = real(leftvec.col(lambda1));
  realleftvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-14 with 0
  arma::vec rlvabs = abs(realleftvec);
  
  arma::uvec rlvabsalt = find(rlvabs);
  int rlvminelem = rlvabsalt(0);
  double rlvmin = realleftvec(rlvminelem);
  
  arma::vec wcorr (rvel);
  arma::vec vcorr (rvel);
  arma::vec vwprod (rvel);
  arma::mat smat (rvel, rvel);
  smat.zeros();
  
  // This loop and the following line create the scalar product vw
  for (int i = 0; i < rvel; i++) {
    wcorr(i) = realrightvec(i) / rvsum;
    vcorr(i) = realleftvec(i) / rlvmin;
    
    vwprod(i) = wcorr(i) * vcorr(i);
  }
  
  double vwscalar = sum(vwprod);
  
  // This loop populates the sensitivity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      
      smat(i, j) = vcorr(i) * wcorr(j) / vwscalar;
    }
  }
  
  return smat;
}

//' Estimate Sensitivities for a Historical LefkoMat Object
//' 
//' \code{sens3hlefko()} returns the sensitivity of lambda with respect
//' to each historical stage-pair in the matrix, and the associated
//' sensitivity for each life history stage. This is accomplished via the 
//' \code{eigs_gen}() function in the C++ Armadillo library.
//' 
//' @param Amat A population projection matrix.
//' @param ahstages An integar vector of unique ahistorical stages.
//' @param hstages An integar vector of unique historical stage pairs.
//' 
//' @return This function returns a list with two sensitivity matrices:
//' \item{h_smat}{Matrix of sensitivities corresponding to the historical
//' matrix.}
//' \item{ah_smat}{Matrix of sensitivities corresponding to the ahistorical
//' matrix.}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List sens3hlefko(arma::mat Amat, DataFrame ahstages, DataFrame hstages) {
  
  arma::uvec stage_id = ahstages["stage_id"];
  arma::uvec h_stage_2 = hstages["stage_id_2"];
  arma::uvec h_stage_1 = hstages["stage_id_1"];
  
  List eigenstuff = decomp3sp(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  //double lambda = max(realeigenvals);
  int lambda1 = realeigenvals.index_max();
  
  // This is the w vector
  cx_mat rightvec = eigenstuff["right_eigenvectors"];
  arma::vec realrightvec = real(rightvec.col(lambda1));
  realrightvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-14 with 0
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  
  // This is the v vector
  cx_mat leftvec = eigenstuff["left_eigenvectors"];
  arma::vec realleftvec = real(leftvec.col(lambda1));
  realleftvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-14 with 0
  arma::vec rlvabs = abs(realleftvec);
  
  arma::uvec rlvabsalt = find(rlvabs);
  int rlvminelem = rlvabsalt(0);
  double rlvmin = realleftvec(rlvminelem);
  
  int ahstagelength = stage_id.n_elem;
  
  arma::vec wcorr (rvel);
  arma::vec vcorr (rvel);
  arma::vec vwprod (rvel);
  arma::mat smat (rvel, rvel);
  smat.zeros();
  
  arma::vec wcorrah (ahstagelength);
  wcorrah.zeros();
  arma::vec vcorrah (ahstagelength);
  vcorrah.zeros();
  arma::vec vwprodah (ahstagelength);
  vwprodah.zeros();
  arma::mat ahsens(ahstagelength, ahstagelength);
  ahsens.zeros();
  
  int ahrows {0};
  //int ahcols {0};
  
  // This loop and the following line create the scalar product vw and the ahistorical stable stage distribution w
  for (int i = 0; i < rvel; i++) {
    wcorr(i) = realrightvec(i) / rvsum;
    vcorr(i) = realleftvec(i) / rlvmin;
    
    vwprod(i) = wcorr(i) * vcorr(i);
    
    ahrows = h_stage_2(i) - 1;
    
    wcorrah(ahrows) = wcorrah(ahrows) + wcorr(i);
  }
  
  double vwscalar = sum(vwprod);
  
  // This loop and the following line create the scalar product vw and the ahistorical stable stage distribution w
  for (int i = 0; i < rvel; i++) {
    ahrows = h_stage_2(i) - 1;
    
    vcorrah(ahrows) = vcorr(i) * wcorr(i) / wcorrah(ahrows);
    
  }
  
  // This loop populates the historical and ahistorical elasticity matrices
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      
      smat(i, j) = vcorr(i) * wcorr(j) / vwscalar;
      
    }
  }
  
  for (int i = 0; i < ahstagelength; i++) {
    vwprodah(i) = wcorrah(i) * vcorrah(i);
  }
  
  double vwscalarah = sum(vwprodah);
  
  for (int i = 0; i < ahstagelength; i++) {
    for (int j = 0; j < ahstagelength; j++) {
      ahsens(i, j) = vcorrah(i) * wcorrah(j) / vwscalarah;
    }
  }
  
  List output = List::create(Named("h_smat") = smat, _["ah_smat"] = ahsens);
  
  return output;
}

//' Estimate Elasticities for a Dense Population Matrix
//' 
//' \code{elas3matrix()} returns the elasticity of lambda with respect
//' to each element in a dense matrix. This is accomplished via the
//' \code{eig_gen}() function in the C++ Armadillo library.
//' 
//' @param Amat A population projection matrix.
//' 
//' @return This function returns a matrix of elasticities. 
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat elas3matrix(arma::mat Amat) {
  
  List eigenstuff = decomp3(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  double lambda = max(realeigenvals);
  int lambda1 = realeigenvals.index_max();
  
  // This is the w vector
  cx_mat rightvec = eigenstuff["right_eigenvectors"];
  arma::vec realrightvec = real(rightvec.col(lambda1));
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  
  // This is the v vector
  cx_mat leftvec = eigenstuff["left_eigenvectors"];
  arma::vec realleftvec = real(leftvec.col(lambda1));
  arma::vec rlvabs = abs(realleftvec);
  
  arma::uvec rlvabsalt = find(rlvabs);
  int rlvminelem = rlvabsalt(0);
  double rlvmin = realleftvec(rlvminelem);
  
  arma::vec wcorr (rvel);
  arma::vec vcorr (rvel);
  arma::vec vwprod (rvel);
  arma::mat emat (rvel, rvel);
  emat.zeros();
  
  // This loop and the following line create the scalar product vw
  for (int i = 0; i < rvel; i++) {
    wcorr(i) = realrightvec(i) / rvsum;
    vcorr(i) = realleftvec(i) / rlvmin;
    
    vwprod(i) = wcorr(i) * vcorr(i);
  }
  
  double vwscalar = sum(vwprod);
  
  // This loop populates the elasticity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      
      emat(i, j) = (vcorr(i) * wcorr(j) * Amat(i, j)) / (vwscalar * lambda);
    }
  }
  
  return emat;
}

//' Estimate Elasticities for a Sparse Population Matrix
//' 
//' \code{elas3matrixsp()} returns the elasticity of lambda with respect
//' to each element in a sparse matrix. This is accomplished via the
//' \code{eigs_gen}() function in the C++ Armadillo library.
//' 
//' @param Amat A population projection matrix.
//' 
//' @return This function returns a matrix of elasticities. 
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat elas3matrixsp(arma::mat Amat) {
  
  List eigenstuff = decomp3sp(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  double lambda = max(realeigenvals);
  int lambda1 = realeigenvals.index_max();
  
  // This is the w vector
  cx_mat rightvec = eigenstuff["right_eigenvectors"];
  arma::vec realrightvec = real(rightvec.col(lambda1));
  realrightvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-14 with 0
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  
  // This is the v vector
  cx_mat leftvec = eigenstuff["left_eigenvectors"];
  arma::vec realleftvec = real(leftvec.col(lambda1));
  realleftvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-14 with 0
  arma::vec rlvabs = abs(realleftvec);
  
  arma::uvec rlvabsalt = find(rlvabs);
  int rlvminelem = rlvabsalt(0);
  double rlvmin = realleftvec(rlvminelem);
  
  arma::vec wcorr (rvel);
  arma::vec vcorr (rvel);
  arma::vec vwprod (rvel);
  arma::mat emat (rvel, rvel);
  emat.zeros();
  
  // This loop and the following line create the scalar product vw
  for (int i = 0; i < rvel; i++) {
    wcorr(i) = realrightvec(i) / rvsum;
    vcorr(i) = realleftvec(i) / rlvmin;
    
    vwprod(i) = wcorr(i) * vcorr(i);
  }
  
  double vwscalar = sum(vwprod);
  
  // This loop populates the elasticity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      
      emat(i, j) = (vcorr(i) * wcorr(j) * Amat(i, j)) / (vwscalar * lambda);
    }
  }
  
  return emat;
}

//' Estimate Elasticities for a Historical LefkoMod Object
//' 
//' \code{elas3hlefko()} returns the elasticity of lambda with respect
//' to each historical stage-pair in the matrix, and the summed elasticities
//' for each life history stage. This is accomplished via the \code{eigs_gen}()
//' function in the C++ Armadillo library.
//' 
//' @param Amat A population projection matrix.
//' @param ahstages An integar vector of unique ahistorical stages.
//' @param hstages An integar vector of unique historical stage pairs.
//' 
//' @return This function returns a list with two elasticity matrices:
//' \item{h_emat}{Matrix of elasticities corresponding to the historical matrix.}
//' \item{ah_emat}{Matrix of elasticities corresponding to the ahistorical
//' matrix, but using summed historical elasticities as the basis of estimation.}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List elas3hlefko(arma::mat Amat, DataFrame ahstages, DataFrame hstages) {
  
  arma::uvec stage_id = ahstages["stage_id"];
  arma::uvec h_stage_2 = hstages["stage_id_2"];
  arma::uvec h_stage_1 = hstages["stage_id_1"];
  
  List eigenstuff = decomp3sp(Amat);
  
  cx_vec Eigenvals = eigenstuff["eigenvalues"];
  arma::vec realeigenvals = real(Eigenvals);
  
  double lambda = max(realeigenvals);
  int lambda1 = realeigenvals.index_max();
  
  // This is the w vector
  cx_mat rightvec = eigenstuff["right_eigenvectors"];
  arma::vec realrightvec = real(rightvec.col(lambda1));
  realrightvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-14 with 0
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  
  // This is the v vector
  cx_mat leftvec = eigenstuff["left_eigenvectors"];
  arma::vec realleftvec = real(leftvec.col(lambda1));
  realleftvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-14 with 0
  arma::vec rlvabs = abs(realleftvec);
  
  arma::uvec rlvabsalt = find(rlvabs);
  int rlvminelem = rlvabsalt(0);
  double rlvmin = realleftvec(rlvminelem);
  
  arma::vec wcorr (rvel);
  arma::vec vcorr (rvel);
  arma::vec vwprod (rvel);
  arma::mat emat (rvel, rvel);
  emat.zeros();
  
  // This loop and the following line create the scalar product vw
  for (int i = 0; i < rvel; i++) {
    wcorr(i) = realrightvec(i) / rvsum;
    vcorr(i) = realleftvec(i) / rlvmin;
    
    vwprod(i) = wcorr(i) * vcorr(i);
  }
  
  double vwscalar = sum(vwprod);
  
  // The next few lines set up the empty ahistorical matrix
  int ahstagelength = stage_id.n_elem;
  
  int ahrows {0};
  int ahcols {0};
  
  arma::mat ahelas(ahstagelength, ahstagelength);
  ahelas.zeros();
  
  // This loop populates the historical and ahistorical elasticity matrices
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      
      emat(i, j) = (vcorr(i) * wcorr(j) * Amat(i, j)) / (vwscalar * lambda);
      
      ahrows = h_stage_2(i) - 1;
      ahcols = h_stage_1(i) - 1;
      
      ahelas(ahrows, ahcols) = ahelas(ahrows, ahcols) + emat(i, j);
    }
  }
  
  List output = List::create(Named("h_emat") = emat, _["ah_emat"] = ahelas);
  
  return output;
}

//' Core projection function
//' 
//' Function \code{proj3()} runs the matrix projections used in other functions
//' in package \code{lefko3}.
//' 
//' @param start_mat The starting matrix for the projection.
//' @param start_vec The starting population vector for the projection.
//' @param core_list A list of full projection matrices, corresponding to the 
//' \code{$A} list within a \code{lefkoMat} object.
//' @param mat_order A vector giving the order of matrices to use at each time.
//' @param standardize A logical value stating whether to standardize population
//' size vector to sum to 1 at each estimated time.
//' 
//' @return A matrix in which each row is the population vector at each 
//' projected time.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat proj3(arma::vec start_vec, List core_list, arma::uvec mat_order,
  bool standardize) {
  
  int nostages = start_vec.n_elem;
  int theclairvoyant = mat_order.n_elem;
  arma::vec thechosenone;
  arma::vec theseventhson;
  arma::mat theprophecy;
  arma::sp_mat theprophecy_sp;
  
  arma::mat popproj(nostages, (theclairvoyant + 1)); // This is the population vector
  popproj.zeros();
  
  popproj.col(0) = start_vec;
  thechosenone = start_vec;
  
  theprophecy = as<arma::mat>(core_list[0]);
  
  for (int i = 0; i < theclairvoyant; i++) {
    theseventhson = theprophecy * thechosenone;
    popproj.col(i+1) = theseventhson;
    
    if (standardize) {
      thechosenone = theseventhson / sum(theseventhson);
    } else {
      thechosenone = theseventhson;
    }
    
    theprophecy = as<arma::mat>(core_list[(mat_order(i))]);
  }
  
  return popproj;
}

//' Stochastic population growth rate estimation
//' 
//' Function \code{slambda3()} estimates the stochastic population growth rate,
//' \eqn{a}, defined as the long-term arithmetic mean of the log population 
//' growth estimated per simulated time (as given in equation 2 in Tuljapurkar,
//' Horvitz, and Pascarella 2003). This term is estimated via projection of 
//' randomly sampled matrices, similarly to the procedure outlined in Box 7.4 of
//' Morris and Doak (2002). Can handle both lefkoMat objects and lists of full A
//' matrices. 
//' 
//' @param mpm A matrix projection model of class \code{lefkoMat}, or a list of
//' full matrix projection matrices.
//' @param times Number of iterations to random samples. Defaults to 10,000.
//' @param tweights Numeric vector denoting the probabilistic weightings of
//' annual matrices. Defaults to equal weighting among times.
//' 
//' @return A data frame with the following variables:
//' 
//' \item{pop}{The identity of the population.}
//' \item{patch}{The identity of the patch.}
//' \item{a}{Estimate of stochastic growth rate, estimated as the arithmetic
//' mean of the log population growth rate across simulated times.}
//' \item{var}{The estimated variance of a.}
//' \item{sd}{The standard deviation of a.}
//' \item{se}{The standard error of a.}
//'
//' Stochastic growth rate is estimated both at the patch level and at the
//' population level. Population level estimates will be noted at the end of the
//' data frame with 0 entries for patch designation.
//' 
//' @examples
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
//' rep_cyp_raw <- matrix(0, 11, 11)
//' rep_cyp_raw[1:2,7:11] <- 0.5
//' 
//' cypover3r <- overwrite(stage3 = c("SD", "SD", "P1", "P1", "P2", "P3", "SL", 
//'     "SL", "SL", "D", "XSm", "Sm", "D", "XSm", "Sm"), 
//'   stage2 = c("SD", "SD", "SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", 
//'     "SL", "SL", "SL", "SL", "SL"),
//'   stage1 = c("SD", "rep", "SD", "rep", "SD", "P1", "P2", "P3", "SL", "P3", 
//'     "P3", "P3", "SL", "SL", "SL"),
//'   eststage3 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", "D", 
//'     "XSm", "Sm"), 
//'   eststage2 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", 
//'     "XSm", "XSm", "XSm"), 
//'   eststage1 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", 
//'     "XSm", "XSm", "XSm"), 
//'   givenrate = c(0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.25, 0.4, 0.4, NA, NA, NA, 
//'     NA, NA, NA), 
//'   type = c("S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", 
//'     "S", "S"))
//' 
//' cypmatrix3r <- rlefko3(data = cypraw_v1, stageframe = cypframe_raw, 
//'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
//'   size = c("size3added", "size2added", "size1added"), 
//'   repmatrix = rep_cyp_raw, overwrite = cypover3r, yearcol = "year2", 
//'   patchcol = "patchid", indivcol = "individ")
//' 
//' cypstoch <- slambda3(cypmatrix3r)
//' 
//' @export slambda3
// [[Rcpp::export]]
DataFrame slambda3(List mpm, int times = 10000, 
  Nullable<NumericVector> tweights = R_NilValue) {
  
  int theclairvoyant {0};
  
  theclairvoyant = times;
  
  if (theclairvoyant < 1) {
    stop("Option must equal a positive integer.");
  }
  
  if (!Rf_isMatrix(mpm[0])) {
    List amats = mpm["A"];
    List umats = mpm["U"];
    List fmats = mpm["F"];
    DataFrame stageframe = as<DataFrame>(mpm["ahstages"]);
    DataFrame hstages = as<DataFrame>(mpm["hstages"]);
    DataFrame labels = as<DataFrame>(mpm["labels"]);
    
    bool historical;
    
    if (hstages.length() > 1) {
      historical = true;
    } else {
      historical = false;
    }
    
    if (labels.length() < 3) {
      stop("Function 'slambda3' requires annual matrices. This lefkoMat object appears to be a set of mean matrices, and lacks annual matrices.");
    }
    
    StringVector poporder = labels["pop"];
    StringVector patchorder = labels["patch"];
    IntegerVector yearorder = labels["year2"];
    
    int loysize = poporder.length();
    StringVector poppatch = clone(poporder);
    
    for (int i = 0; i < loysize; i++) {
      poppatch(i) += " ";
      poppatch(i) += patchorder(i);
    }
    
    StringVector uniquepops = sort_unique(poporder);
    StringVector uniquepoppatches = sort_unique(poppatch);
    IntegerVector uniqueyears = sort_unique(yearorder);
    IntegerVector popc = match(poporder, uniquepops) - 1;
    IntegerVector poppatchc = match(poppatch, uniquepoppatches) - 1;
    IntegerVector year2c = match(yearorder, uniqueyears) - 1;
    int yl = uniqueyears.length();
    
    arma::vec twinput;
    
    if (tweights.isNotNull()) {
      if (as<NumericVector>(tweights).length() != yl) {
        stop("Time weight vector must be the same length as the number of times represented in the lefkoMat object used as input.");
      }
      twinput = as<arma::vec>(tweights);
    } else {
      twinput.resize(yl);
      twinput.ones();
    } // At the end of this, we have an arma::vec holding with whatever the user supplied, or a 1 for each time
    
    arma::uvec armapopc = as<arma::uvec>(popc);
    arma::uvec armapoppatchc = as<arma::uvec>(poppatchc);
    arma::uvec armayear2c = as<arma::uvec>(year2c);
    arma::uvec patchesinpop(loysize);
    arma::uvec yearsinpatch(loysize);
    patchesinpop.zeros();
    yearsinpatch.zeros();
    
    for (int i = 0; i < loysize; i++) {
      arma::uvec animalhouse = find(armapopc == armapopc(i));
      arma::uvec christmasvacation = armapoppatchc.elem(animalhouse);
      arma::uvec summervacation = unique(christmasvacation);
      int togaparty = summervacation.n_elem;
      
      patchesinpop(i) = togaparty;
      
      arma::uvec ninebelowzero = find(armapoppatchc == armapoppatchc(i));
      arma::uvec thedamned = armayear2c.elem(ninebelowzero);
      arma::uvec motorhead = unique(thedamned);
      int dexysmidnightrunners = motorhead.n_elem;
      
      yearsinpatch(i) = dexysmidnightrunners;
    }
    
    DataFrame listofyears = DataFrame::create(Named("pop") = poporder,
      _["patch"] = patchorder, _["year2"] = yearorder, _["poppatch"] = poppatch,
      _["popc"] = popc, _["poppatchc"] = poppatchc, _["year2c"] = year2c, 
      _["patchesinpop"] = patchesinpop, _["yearsinpatch"] = yearsinpatch);
    
    List mean_lefkomat;
    
    if (hstages.length() == 1) {
      mean_lefkomat = geodiesel(listofyears, umats, fmats, stageframe, 1, 1);
    } else {
      mean_lefkomat = turbogeodiesel(listofyears, umats, fmats, stageframe, hstages, 1, 1);
    }
    
    // Here we take the matrices corresponding to each individual patch, run the simulation, and
    // estimate all descriptive metrics
    List meanamats = mean_lefkomat["A"];
    List mmlabels = mean_lefkomat["labels"];
    StringVector mmpops = mmlabels["pop"];
    StringVector mmpatches = mmlabels["patch"];
    
    arma::mat thechosenone = as<arma::mat>(meanamats[0]);
    int meanmatsize = thechosenone.n_elem;
    int meanmatrows = thechosenone.n_rows;
    arma::vec startvec;
    int trials = meanamats.length();
    
    arma::uvec ppcindex = as<arma::uvec>(poppatchc);
    arma::uvec allppcs = as<arma::uvec>(sort_unique(poppatchc));
    int allppcsnem = allppcs.n_elem;
    
    arma::mat slmat(theclairvoyant, trials);
    slmat.zeros();
    
    arma::vec sl_mean(trials);
    arma::vec sl_var(trials);
    arma::vec sl_sd(trials);
    arma::vec sl_se(trials);
    sl_mean.zeros();
    sl_var.zeros();
    sl_sd.zeros();
    sl_se.zeros();
    
    arma::vec tweights_corr = twinput / sum(twinput);
    
    for (int i= 0; i < allppcsnem; i++) {
      thechosenone = as<arma::mat>(meanamats[i]);
      
      arma::uvec thenumbersofthebeast = find(ppcindex == allppcs(i));
      arma::uvec theprophecy = Rcpp::RcppArmadillo::sample(thenumbersofthebeast, theclairvoyant, true, tweights_corr);
      
      if (historical == true) {
        startvec = ss3matrixsp(thechosenone);
      } else {
        startvec = ss3matrix(thechosenone);
      }
      
      arma::mat projection = proj3(startvec, amats, theprophecy, 1);
      
      for (int j = 0; j < theclairvoyant; j++) {
        //slmat(j,i) = sum(projection.col(j+1));
        double madness = sum(projection.col(j+1));
        slmat(j,i) = log(madness);
      }
      
      sl_mean(i) = mean(slmat.col(i));
      sl_var(i) = var(slmat.col(i));
      sl_sd(i) = stddev(slmat.col(i));
      sl_se(i) = sl_sd(i) / sqrt(static_cast<double>(theclairvoyant));
    }
    
    int pop_est {1};
    StringVector allpops = unique(poporder);
    arma::uvec popmatch(loysize);
    arma::uvec yearmatch(loysize);
    popmatch.zeros();
    yearmatch.zeros();
    List meanmatyearlist(uniqueyears.length());
    
    if (allppcsnem > 1) { // This checks if there are any pop-mean matrices separate from the the patch means
      pop_est = trials - allppcsnem;
      
      for (int i = 0; i < pop_est; i++) { // This loop goes through each population
        thechosenone = as<arma::mat>(meanamats[allppcsnem + i]);
        
        if (historical == true) {
          startvec = ss3matrixsp(thechosenone);
        } else {
          startvec = ss3matrix(thechosenone);
        }
        
        for (int j = 0; j < loysize; j++) { // This checks which A matrices match the current population in the loop
          if (poporder(j) == allpops(i)) {
            popmatch(j) = 1;
          } else {
            popmatch(j) = 0;
          }
        }
        
        arma::uvec neededmatspop = find(popmatch);
        
        for (int j = 0; j < yl; j++) { // This loop checks for each year and develops a matrix mean across patches
          for (int k = 0; k < loysize; k++) { // This inner loop develops a vector to find all matrices corresponding to the current year
            
            if (yearorder(k) == uniqueyears(j)) {
              yearmatch(k) = 1;
            } else {
              yearmatch(k) = 0;
            }
          }
          
          arma::uvec neededmatsyear = find(yearmatch);
          arma::uvec crankybanky = intersect(neededmatsyear, neededmatspop);
          // This vector catches the indices of matrices that match the current year and population
          int crankybankynem = crankybanky.n_elem;
          
          arma::mat crossmat(meanmatsize, crankybankynem);
          crossmat.zeros();
          
          for (int j = 0; j < crankybankynem; j++) {
            crossmat.col(j) = as<arma::vec>(amats(crankybanky(j)));
          }
          
          arma::vec happymedium(meanmatsize);
          happymedium.zeros();
          
          for (int j = 0; j < meanmatsize; j++) {
            for (int k = 0; k < crankybankynem; k++) {
              happymedium(j) = happymedium(j) + crossmat(j, k) / (crankybankynem);
            }
          }
          
          arma::mat finalyearmat = happymedium;
          finalyearmat.reshape(meanmatrows, meanmatrows);
          
          meanmatyearlist(j) = finalyearmat;
        }
        
        int numyearsused = meanmatyearlist.length();
        arma::uvec choicevec = linspace<arma::uvec>(0, (numyearsused - 1), numyearsused);
        arma::uvec theprophecy = Rcpp::RcppArmadillo::sample(choicevec, theclairvoyant, true, tweights_corr);
        
        arma::mat projection = proj3(startvec, meanmatyearlist, theprophecy, 1);
        
        for (int j = 0; j < theclairvoyant; j++) {
          double madness = sum(projection.col(j+1));
          slmat(j,(allppcsnem +i)) = log(madness);
        }
        
        sl_mean((allppcsnem +i)) = mean(slmat.col((allppcsnem +i)));
        sl_var((allppcsnem +i)) = var(slmat.col((allppcsnem +i)));
        sl_sd((allppcsnem +i)) = stddev(slmat.col((allppcsnem +i)));
        sl_se((allppcsnem +i)) = sl_sd((allppcsnem +i)) / sqrt(static_cast<double>(theclairvoyant));
      }
    }
    return DataFrame::create(_["pop"] = mmpops, _["patch"] = mmpatches,
      _["a"] = sl_mean, _["var"] = sl_var, _["sd"] = sl_sd, _["se"] = sl_se);
  } else {
    
    List amats = mpm;
    
    int yl = amats.length();
    arma::mat firstmat = as<arma::mat>(amats[0]);
    int matrows = firstmat.n_rows;
    int matcols = firstmat.n_cols;
    
    bool historical;
    
    if (matrows > 400) {
      historical = true;
    } else {
      historical = false;
    }
    
    arma::uvec uniqueyears(yl);
    for (int i = 0; i < yl; i++) {
      uniqueyears(i) = i;
    }
    
    arma::vec twinput;
    
    if (matrows != matcols) {
      stop("Supplied matrices must be square. Please check matrix dimensions and fix.");
    }
    
    if (tweights.isNotNull()) {
      if (as<NumericVector>(tweights).length() != yl) {
        stop("Time weight vector must be the same length as the number of times represented in the lefkoMat object used as input.");
      }
      twinput = as<arma::vec>(tweights);
    } else {
      twinput.resize(yl);
      twinput.ones();
    } // At the end of this, we have an arma::vec holding with whatever the user supplied, or a 1 for each time
    
    // Now we create the mean matrix
    arma::mat thechosenone(matrows, matcols);
    thechosenone.zeros();
    
    for (int i = 0; i < yl; i++) {
      arma::mat columnified = as<arma::mat>(amats[i]);
      thechosenone = thechosenone + (columnified / yl);
    }
    
    // Here we take the matrices corresponding to each individual patch, run the simulation, and
    // estimate all descriptive metrics
    arma::vec startvec;
    int trials {1};
    
    arma::mat slmat(theclairvoyant, trials);
    slmat.zeros();
    
    arma::vec sl_mean(trials);
    arma::vec sl_var(trials);
    arma::vec sl_sd(trials);
    arma::vec sl_se(trials);
    sl_mean.zeros();
    sl_var.zeros();
    sl_sd.zeros();
    sl_se.zeros();
    
    arma::vec tweights_corr = twinput / sum(twinput);
    
    arma::uvec thenumbersofthebeast = uniqueyears;
    arma::uvec theprophecy = Rcpp::RcppArmadillo::sample(thenumbersofthebeast, theclairvoyant, true, tweights_corr);
    
    if (historical == true) {
      startvec = ss3matrixsp(thechosenone);
    } else {
      startvec = ss3matrix(thechosenone);
    }
    
    arma::mat projection = proj3(startvec, amats, theprophecy, 1);
    
    for (int j = 0; j < theclairvoyant; j++) {
      slmat(j,0) = sum(projection.col(j+1));
    }
    
    sl_mean(0) = mean(slmat.col(0));
    sl_var(0) = var(slmat.col(0));
    sl_sd(0) = stddev(slmat.col(0));
    sl_se(0) = sl_sd(0) / sqrt(static_cast<double>(theclairvoyant));
    
    CharacterVector mmpops(1);
    CharacterVector mmpatches(1);
    mmpops(0) = "1";
    mmpatches(0) = "0";
    
    return DataFrame::create(_["pop"] = mmpops, _["patch"] = mmpatches,
      _["a"] = sl_mean, _["var"] = sl_var, _["sd"] = sl_sd, _["se"] = sl_se);
  }
}

//' Creates Size Index for Elasticity Summaries of hMPMs
//' 
//' Function \code{bambi3()} creates an index of estimable elements in
//' historical matrices, and details the kind of transition that it is.
//' 
//' @param stages This is the core stageframe held by \code{mats}, equivalent to
//' \code{ahstages}.
//' @param hstages This is the \code{hstages} object held by \code{mats}.
//' 
//' @return A data frame with the following elements:
//' \item{index}{Vector index of matrix element in C++ terms.}
//' \item{transition}{Category of transition.}
//' \item{size3}{Size in time \emph{t}+1.}
//' \item{repstatus3}{Reproductive status in time \emph{t}+1.}
//' \item{entrystatus3}{Entry status in time \emph{t}+1.}
//' \item{size2}{Size in time \emph{t}.}
//' \item{repstatus2}{Reproductive status in time \emph{t}.}
//' \item{entrystatus2}{Entry status in time \emph{t}.}
//' \item{size1}{Size in time \emph{t}-1.}
//' \item{repstatus1}{Reproductive status in time \emph{t}11.}
//' \item{entrystatus1}{Entry status in time \emph{t}-1.}
//'
//' The kind of transitions conforms to the following code: \code{10}: full
//' stasis, \code{11}: stasis to growth, \code{12}: full growth, \code{13}:
//' growth to stasis, \code{14}: stasis to shrinkage, \code{15}: full shrinkage,
//' \code{16}: shrinkage to stasis, \code{17}: growth to shrinkage, \code{18}:
//' shrinkage to growth, \code{20}: stasis to fecundity, \code{21}: growth to
//' fecundity, \code{22}: shrinkage to fecundity, \code{23}: fecundity to
//' stasis, \code{24}: fecundity to growth, \code{25}: fecundity to shrinkage,
//' \code{26}: fecundity to fecundity.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
DataFrame bambi3(DataFrame stages, DataFrame hstages) {
  
  StringVector stagenames = stages["stage"];
  arma::uvec astages = stages["stage_id"];
  arma::vec sizes = stages["original_size"];
  arma::uvec repstatus = stages["repstatus"];
  arma::uvec entrystage = stages["entrystage"];
  int numstages = astages.n_elem;
  
  arma::uvec hstage3in = hstages["stage_id_2"];
  arma::uvec hstage2nin = hstages["stage_id_1"];
  //StringVector hstagenames3 = hstages["stage_2"];
  //StringVector hstagenames2 = hstages["stage_1"];
  int numhstages = hstage3in.n_elem;
  
  arma::uvec stages_order = astages - 1;
  arma::uvec hstage3_order = hstage3in - 1;
  arma::uvec hstage2_order = hstage2nin - 1;
  
  int predictedsize = numstages * numstages * numstages;
  
  arma::ivec hsindexl(predictedsize);
  hsindexl.fill(-1);
  
  arma::uvec transition_type(predictedsize);
  transition_type.zeros();
  
  arma::vec size1(predictedsize);
  arma::vec size2(predictedsize);
  arma::vec size3(predictedsize);
  size1.fill(-1);
  size2.fill(-1);
  size3.fill(-1);
  
  arma::uvec repstatus1(predictedsize);
  arma::uvec repstatus2(predictedsize);
  arma::uvec repstatus3(predictedsize);
  repstatus1.zeros();
  repstatus2.zeros();
  repstatus3.zeros();
  
  arma::uvec entrystatus1(predictedsize);
  arma::uvec entrystatus2(predictedsize);
  arma::uvec entrystatus3(predictedsize);
  entrystatus1.zeros();
  entrystatus2.zeros();
  entrystatus3.zeros();
  
  StringVector longnames3(predictedsize);
  StringVector longnames2(predictedsize);
  StringVector longnames1(predictedsize);
  
  int counter = 0;
  
  for (int i1 = 0; i1 < numhstages; i1++) {
    for (int i2 = 0; i2 < numhstages; i2++) {
      if (hstage3in(i1) == hstage2nin(i2)) {
        
        hsindexl(counter) = (i1 * numhstages) + i2;
        
        int stage1 = hstage2_order(i2);
        longnames1(counter) = stagenames(stage1);
        size1(counter) = sizes(stage1);
        repstatus1(counter) = repstatus(stage1);
        entrystatus1(counter) = entrystage(stage1);
        
        int stage2 = hstage2_order(i1);
        longnames2(counter) = stagenames(stage2);
        size2(counter) = sizes(stage2);
        repstatus2(counter) = repstatus(stage2);
        entrystatus2(counter) = entrystage(stage2);
        
        int stage3 = hstage3_order(i2);
        longnames3(counter) = stagenames(stage3);
        size3(counter) = sizes(stage3);
        repstatus3(counter) = repstatus(stage3);
        entrystatus3(counter) = entrystage(stage3);
        
        if (entrystatus3(counter) == 1 && repstatus2(counter) == 1) {
          if (entrystatus2(counter) == 1 && repstatus1(counter) == 1) {
            transition_type(counter) = 26; // Fecundity to fecundity
          } else if (size2(counter) == size1(counter)) {
            transition_type(counter) = 20; // Stasis to fecundity
          } else if (size2(counter) > size1(counter)) {
            transition_type(counter) = 21; // Growth to fecundity
          } else if (size2(counter) < size1(counter)) {
            transition_type(counter) = 22; // Shrinkage to fecundity
          }
        } else if (entrystatus2(counter) == 1 && repstatus1(counter) == 1) {
          if (size3(counter) == size2(counter)) {
            transition_type(counter) = 23; // Fecundity to stasis
          } else if (size3(counter) > size2(counter)) {
            transition_type(counter) = 24; // Fecundity to growth
          } else if (size3(counter) < size2(counter)) {
            transition_type(counter) = 25; // Fecundity to shrinkage
          }
        } else if (size3(counter) == size2(counter) && size2(counter) == size1(counter)) {
          transition_type(counter) = 10; // Full stasis
        } else if (size3(counter) > size2(counter) && size2(counter) == size1(counter)) {
          transition_type(counter) = 11; // Stasis to growth
        } else if (size3(counter) > size2(counter) && size2(counter) > size1(counter)) {
          transition_type(counter) = 12; // Full growth
        } else if (size3(counter) == size2(counter) && size2(counter) > size1(counter)) {
          transition_type(counter) = 13; // Growth to stasis
        } else if (size3(counter) < size2(counter) && size2(counter) == size1(counter)) {
          transition_type(counter) = 14; // Stasis to shrinkage
        } else if (size3(counter) < size2(counter) && size2(counter) < size1(counter)) {
          transition_type(counter) = 15; // Full shrinkage
        } else if (size3(counter) == size2(counter) && size2(counter) < size1(counter)) {
          transition_type(counter) = 16; // Shrinkage to stasis
        } else if (size3(counter) < size2(counter) && size2(counter) > size1(counter)) {
          transition_type(counter) = 17; // Growth to shrinkage
        } else if (size3(counter) > size2(counter) && size2(counter) < size1(counter)) {
          transition_type(counter) = 18; // Shrinkage to growth
        }
        
        counter++;
      }
    }
  }
  
  StringVector names3(counter);
  StringVector names2(counter);
  StringVector names1(counter);
  
  for (int i = 0; i < counter; i++) {
    names3(i) = longnames3(i);
    names2(i) = longnames2(i);
    names1(i) = longnames1(i);
  }
  
  arma::uvec targetindices = find(hsindexl > -1);
  arma::ivec hsindex = hsindexl.elem(targetindices);
  arma::uvec t_type = transition_type.elem(targetindices);
  arma::vec size3c = size3.elem(targetindices);
  arma::vec size2c = size2.elem(targetindices);
  arma::vec size1c = size1.elem(targetindices);
  arma::uvec r_status3 = repstatus3.elem(targetindices);
  arma::uvec r_status2 = repstatus2.elem(targetindices);
  arma::uvec r_status1 = repstatus1.elem(targetindices);
  arma::uvec e_status3 = entrystatus3.elem(targetindices);
  arma::uvec e_status2 = entrystatus2.elem(targetindices);
  arma::uvec e_status1 = entrystatus1.elem(targetindices);
  
  DataFrame output = DataFrame::create(Named("index") = hsindex, _["transition"] = t_type,
    _["stage3"] = names3, _["size3"] = size3c, _["repstatus3"] = r_status3, _["entrystatus3"] = e_status3,
    _["stage2"] = names2, _["size2"] = size2c, _["repstatus2"] = r_status2, _["entrystatus2"] = e_status2,
    _["stage1"] = names1, _["size1"] = size1c, _["repstatus1"] = r_status1, _["entrystatus1"] = e_status1);
  
  return output;
}

//' Creates Size Index for Elasticity Summaries of ahMPMs
//' 
//' Function \code{bambi2()} creates an index of estimable elements in
//' ahistorical matrices, and details the kind of transition that it is.
//' 
//' @param stages This is the core stageframe held by \code{mats}, equivalent to
//' \code{ahstages}.
//' 
//' @return A data frame with the following elements:
//' \item{index}{Vector index of matrix element in C++ terms.}
//' \item{transition}{Category of transition.}
//' \item{stage3}{Stage in time \emph{t}+1.}
//' \item{size3}{Size in time \emph{t}+1.}
//' \item{repstatus3}{Reproductive status in time \emph{t}+1.}
//' \item{entrystatus3}{Entry status in time \emph{t}+1.}
//' \item{stage2}{Stage in time \emph{t}.}
//' \item{size2}{Size in time \emph{t}.}
//' \item{repstatus2}{Reproductive status in time \emph{t}.}
//' \item{entrystatus2}{Entry status in time \emph{t}.}
//'
//' The kind of transitions conforms to the following code: \code{1}: stasis, 
//' \code{2}: growth, \code{3}: shrinkage, \code{4}: fecundity.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
DataFrame bambi2(DataFrame stages) {
  
  StringVector stagenames = stages["stage"];
  arma::uvec astages = stages["stage_id"];
  arma::vec sizes = stages["original_size"];
  arma::uvec repstatus = stages["repstatus"];
  arma::uvec entrystage = stages["entrystage"];
  int numstages = astages.n_elem;
  
  arma::uvec stages_order = astages - 1;
  
  int predictedsize = numstages * numstages;
  
  arma::ivec ahsindexl(predictedsize);
  ahsindexl.fill(-1);
  
  arma::uvec transition_type(predictedsize);
  transition_type.zeros();
  
  StringVector longstages3(predictedsize);
  StringVector longstages2(predictedsize);
  
  arma::vec size2(predictedsize);
  arma::vec size3(predictedsize);
  size2.fill(-1);
  size3.fill(-1);
  
  arma::uvec repstatus2(predictedsize);
  arma::uvec repstatus3(predictedsize);
  repstatus2.zeros();
  repstatus3.zeros();
  
  arma::uvec entrystatus2(predictedsize);
  arma::uvec entrystatus3(predictedsize);
  entrystatus2.zeros();
  entrystatus3.zeros();
  
  int counter = 0;
  
  for (int i1 = 0; i1 < numstages; i1++) {
    for (int i2 = 0; i2 < numstages; i2++) {
      
      ahsindexl(counter) = (i1 * numstages) + i2;
      
      int stage2 = stages_order(i1);
      longstages2(counter) = stagenames(stage2);
      size2(counter) = sizes(stage2);
      repstatus2(counter) = repstatus(stage2);
      entrystatus2(counter) = entrystage(stage2);
      
      int stage3 = stages_order(i2);
      longstages3(counter) = stagenames(stage3);
      size3(counter) = sizes(stage3);
      repstatus3(counter) = repstatus(stage3);
      entrystatus3(counter) = entrystage(stage3);
      
      if (entrystatus3(counter) == 1 && repstatus2(counter) == 1) {
        transition_type(counter) = 4; // Fecundity
      } else if (size3(counter) == size2(counter)) {
        transition_type(counter) = 1; // Stasis
      } else if (size3(counter) > size2(counter)) {
        transition_type(counter) = 2; // Growth
      } else if (size3(counter) < size2(counter)) {
        transition_type(counter) = 3; // Shrinkage
      }
      
      counter++;
    }
  }
  
  arma::uvec targetindices = find(ahsindexl > -1);
  arma::ivec ahsindex = ahsindexl.elem(targetindices);
  
  DataFrame output = DataFrame::create(Named("index") = ahsindex, _["transition"] = transition_type,
    _["stage3"] = longstages3, _["size3"] = size3, _["repstatus3"] = repstatus3, _["entrystatus3"] = entrystatus3,
    _["stage2"] = longstages2, _["size2"] = size2, _["repstatus2"] = repstatus2, _["entrystatus2"] = entrystatus2);
  
  return output;
}

//' Creates Summary Data for Elasticity Matrix Inputs
//' 
//' Function \code{demolition3()} sums elasticity values from elasticity
//' matrices according to the categories developed by functions \code{bambi2()}
//' and \code{bambi3()}.
//' 
//' @param e_amat A single elasticity matrix.
//' @param amat The A matrix corresponding to \code{e_amat}.
//' @param fmat The F matrix corresponding to \code{e_amat}.
//' @param bambesque This is the output from \code{bambi2()} or \code{bambi3()}
//' corresponding to the current lefkoMat object.
//' 
//' @return A list with two data frames, one showing the summed elasticities for
//' the historical matrix supplied (if supplied), and the other showing the
//' ahistorical summary of the historical matrix or the summed elasticities of
//' a supplied ahistorical elasticity matrix. Note that the elasticity of
//' fecundity transitions will be split with any co-occurring survival
//' transition proportionately.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List demolition3(arma::mat e_amat, arma::mat amat, arma::mat fmat, 
  DataFrame bambesque) {
  
  arma::uvec eindices = bambesque["index"];
  arma::uvec categories = bambesque["transition"];
  
  int e_amatsize = e_amat.n_elem;
  int maxelem = static_cast<int>(eindices.max());
  int minindex = static_cast<int>(categories.min());
  
  if (maxelem > e_amatsize) {
    stop("Supplied info does not seem to correspond to current matrix inputs.");
  }
  
  arma::mat corr_mat = amat;
  arma::uvec z_indices = find(corr_mat == 0);
  int z_indicesnem = z_indices.n_elem;
  
  for (int i = 0; i < z_indicesnem; i++) {
    corr_mat(z_indices(i)) = 1.0;
  }
  
  arma::mat fec_fraction = fmat / corr_mat;
  
  DataFrame histout;
  DataFrame ahistout;
  
  StringVector histcats {"Full stasis", "Stasis to growth", "Full growth",
    "Growth to stasis", "Stasis to shrinkage", "Full shrinkage",
    "Shrinkage to stasis", "Growth to shrinkage", "Shrinkage to growth",
    "Stasis to fecundity", "Growth to fecundity", "Shrinkage to fecundity",
    "Fecundity to stasis", "Fecundity to growth", "Fecundity to shrinkage",
    "Fecundity to fecundity"};
  arma::uvec histcatnums {10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21, 22, 23,
    24, 25, 26};
  arma::vec histsums(16);
  histsums.zeros();
  arma::vec hc_ahistsums(4);
  hc_ahistsums.zeros();
  
  StringVector ahistcats {"Stasis", "Growth", "Shrinkage", "Fecundity"};
  arma::uvec ahistcatnums {1, 2, 3, 4};
  arma::vec ahistsums(4);
  ahistsums.zeros();
  
  //double summerator {0};
  
  if (minindex > 9) {
    
    arma::vec size3 = bambesque["size3"];
    arma::vec size2 = bambesque["size2"];
    arma::vec size1 = bambesque["size1"];
    
    arma::uvec repstatus3 = bambesque["repstatus3"];
    arma::uvec repstatus2 = bambesque["repstatus2"];
    arma::uvec repstatus1 = bambesque["repstatus1"];
    
    arma::uvec entrystatus3 = bambesque["entrystatus3"];
    arma::uvec entrystatus2 = bambesque["entrystatus2"];
    arma::uvec entrystatus1 = bambesque["entrystatus1"];
    
    for (int i = 0; i < 16; i++) {
      arma::uvec currentguys = find(categories == histcatnums(i));
      int currentguysnem = currentguys.n_elem;
      
      if (histcatnums(i) == 20 || histcatnums(i) == 21 || histcatnums(i) == 22 || histcatnums(i) == 26) { // Fecundity transitions

        // Some transitions may be combinations of fecundity and survival, requiring
        // the associated elasticities to be split. This next section does that
        for (int j = 0; j < currentguysnem; j++) {
          double this_guy = eindices(currentguys(j));
          
          if (fec_fraction(this_guy) == 1) {
            hc_ahistsums(3) += (e_amat(this_guy));
            histsums(i) += (e_amat(this_guy));
          } else {
            hc_ahistsums(3) += (e_amat(this_guy) * fec_fraction(this_guy));
            histsums(i) += (e_amat(this_guy) * fec_fraction(this_guy));
            
            arma::uvec counter = find(eindices == this_guy);
            
            if (entrystatus2(counter(0)) == 1 && repstatus1(counter(0)) == 1) {
              if (size3(counter(0)) == size2(counter(0))) {
                hc_ahistsums(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(12) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              } else if (size3(counter(0)) > size2(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(13) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              } else if (size3(counter(0)) < size2(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(14) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              }
            } else if (size3(counter(0)) == size2(counter(0)) && size2(counter(0)) == size1(counter(0))) {
              hc_ahistsums(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              histsums(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
            } else if (size3(counter(0)) > size2(counter(0)) && size2(counter(0)) == size1(counter(0))) {
              hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              histsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
            } else if (size3(counter(0)) > size2(counter(0)) && size2(counter(0)) > size1(counter(0))) {
              hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              histsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
            } else if (size3(counter(0)) == size2(counter(0)) && size2(counter(0)) > size1(counter(0))) {
              hc_ahistsums(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              histsums(3) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
            } else if (size3(counter(0)) < size2(counter(0)) && size2(counter(0)) == size1(counter(0))) {
              hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              histsums(4) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
            } else if (size3(counter(0)) < size2(counter(0)) && size2(counter(0)) < size1(counter(0))) {
              hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              histsums(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
            } else if (size3(counter(0)) == size2(counter(0)) && size2(counter(0)) < size1(counter(0))) {
              hc_ahistsums(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              histsums(6) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
            } else if (size3(counter(0)) < size2(counter(0)) && size2(counter(0)) > size1(counter(0))) {
              hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              histsums(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
            } else if (size3(counter(0)) > size2(counter(0)) && size2(counter(0)) < size1(counter(0))) {
              hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              histsums(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
            }
          }
        }
      } else if (histcatnums(i) == 14 || histcatnums(i) == 15 || histcatnums(i) == 17 || histcatnums(i) == 25) { // Shrinkage transitions
        
        double getoutofdodge = sum(e_amat.elem(eindices(currentguys)));
        histsums(i) += getoutofdodge;
        hc_ahistsums(2) += getoutofdodge;
        
      } else if (histcatnums(i) == 10 || histcatnums(i) == 13 || histcatnums(i) == 16 || histcatnums(i) == 23) { // Stasis transitions
        
        double getoutofdodge = sum(e_amat.elem(eindices(currentguys)));
        histsums(i) += getoutofdodge;
        hc_ahistsums(0) += getoutofdodge;
        
      } else if (histcatnums(i) == 11 || histcatnums(i) == 12 || histcatnums(i) == 18 || histcatnums(i) == 24) { // Growth transitions
        
        double getoutofdodge = sum(e_amat.elem(eindices(currentguys)));
        histsums(i) += getoutofdodge;
        hc_ahistsums(1) += getoutofdodge;
      }
    }
    
    histout = DataFrame::create(Named("category") = histcats, _["elas"] = histsums);
    
    ahistout = DataFrame::create(Named("category") = ahistcats, _["elas"] = hc_ahistsums);
    
  } else {
    histout = R_NilValue;
    
    for (int i = 0; i < 4; i++) {
      arma::uvec currentguys = find(categories == ahistcatnums(i));
      double getoutofdodge = sum(e_amat.elem(eindices(currentguys)));
      ahistsums(i) += getoutofdodge;
    }
    
    ahistout = DataFrame::create(Named("category") = ahistcats, _["elas"] = ahistsums);
  }
  
  List output = List::create(Named("hist") = histout, _["ahist"] = ahistout);
  
  return output;
}
