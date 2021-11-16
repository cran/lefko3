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
//' @return A column vector of specifically called elements from the input
//' matrix.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.flagrantcrap)]]
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
// [[Rcpp::export(.moreflagrantcrap)]]
arma::vec moreflagrantcrap(arma::mat Xmat) {
  
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
//' @param B A sparse matrix. Note that this is assumed to be a population
//' projection matrix, meaning that all values are either 0 or positive.
//' 
//' @return A sparse matrix with non-zero values as logs of the elements in the
//' input matrix.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.spmat_log)]]
arma::sp_mat spmat_log(arma::sp_mat coremat)
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

//' Estimates Mean LefkoMat Object for Historical MPM
//' 
//' Function \code{turbogeodiesel()} estimates mean historical population
//' projection matrices, treating the mean as element-wise arithmetic.
//' 
//' @param loy A data frame denoting the population, patch, and occasion
//' designation for each matrix. Includes a total of 9 variables.
//' @param Umats A matrix with all U matrices turned into columns.
//' @param Fmats A matrix with all F matrices turned into columns.
//' @param hstages This is the \code{hstages} object held by \code{mats}.
//' @param agestages This is the \code{agestages} object held by \code{mats}.
//' @param stages This is the core stageframe held by \code{mats}, equivalent to
//' \code{ahstages}.
//' @param patchmats A logical value stating whether to estimate patch-level
//' means.
//' @param popmats A logical value stating whether to estimate population-level
//' means.
//' 
//' @return A list using the structure of a lefkoMat object.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.turbogeodiesel)]]
List turbogeodiesel(DataFrame loy, List Umats, List Fmats, DataFrame hstages, 
  DataFrame agestages, DataFrame stages, bool patchmats, bool popmats) {
  
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
  int format_int {0};
  arma::uvec astages = stages["stage_id"];
  StringVector stagenames = stages["stage"];
  int numstages = astages.n_elem;
  
  if (stagenames(numstages - 1) == "AlmostBorn") format_int = 1;
  
  arma::uvec hstage3in = hstages["stage_id_2"];
  arma::uvec hstage2nin = hstages["stage_id_1"];
  int numhstages = hstage3in.n_elem;
  
  int predictedsize = 2 * numstages * numstages * numstages;
  
  arma::uvec hsindexl(predictedsize);
  hsindexl.zeros();
  
  counter = 0;
  
  if (format_int == 0) {
    // This bit handles Ehrlen format
    for (int i1 = 0; i1 < numhstages; i1++) {
      for (int i2 = 0; i2 < numhstages; i2++) {
        if (hstage3in(i1) == hstage2nin(i2)) {
          
          hsindexl(counter) = (i1 * numhstages) + i2;
          
          counter++;
        }
      }
    }
  } else {
    // This bit handles deVries format
    for (int i1 = 0; i1 < numhstages; i1++) {
      for (int i2 = 0; i2 < numhstages; i2++) {
        if (hstage3in(i1) == hstage2nin(i2)) {
          
          hsindexl(counter) = (i1 * numhstages) + i2;
          
          counter++;
        } else if (hstage2nin(i2) == numstages || hstage3in(i1) == numstages) {
          
          hsindexl(counter) = (i1 * numhstages) + i2;
          
          counter++;
        }
      }
    }
  }
  
  arma::uvec hsgood = find(hsindexl);
  arma::uvec hsindex = hsindexl.elem(hsgood);
  arma::uvec zerovec(1);
  zerovec.zeros();
  arma::uvec allindices = join_cols(zerovec, hsindex);
  
  // Now we build U and F matrices of element-wise arithmetic means, where
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
    _["hstages"] = hstages, _["agestages"] = agestages, _["ahstages"] = stages,
    _["labels"] = cheatsheet, _["matrixqc"] = matrixqc);
  
  return output;
}

//' Estimates Mean LefkoMat Object for Ahistorical MPM
//' 
//' Function \code{geodiesel()} estimates mean ahistorical population
//' projection matrices, treating the mean as element-wise arithmetic. The
//' function can handle both normal ahistorical MPMs and age x stage ahistorical
//' MPMs.
//' 
//' @param loy A data frame denoting the population, patch, and occasion
//' designation of each matrix. Includes a total of 9 variables.
//' @param Umats A matrix with all U matrices turned into columns.
//' @param Fmats A matrix with all F matrices turned into columns.
//' @param agestages This is the \code{agestages} object held by \code{mats}.
//' @param stages This is the core stageframe held by \code{mats}, equivalent to
//' \code{ahstages}.
//' @param patchmats A logical value stating whether to estimate patch-level
//' means.
//' @param popmats A logical value stating whether to estimate population-level
//' means.
//' 
//' @return A list using the structure of a LefkoMat object.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.geodiesel)]]
List geodiesel(DataFrame loy, List Umats, List Fmats, DataFrame agestages,
  DataFrame stages, bool patchmats, bool popmats) {
  
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
  int initialstages = astages.n_elem;
  
  // Now we will tet for the presence of ages, and determine the matrix dimensions required
  arma::mat initUmat = Umats(0);
  int colsused = initUmat.n_cols;
  int agemultiplier = colsused / initialstages;
  
  int numstages = astages.n_elem * agemultiplier;
  
  // Now we build U and F matrices of element-wise arithmetic means, where
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
  
  // Now we will create the main list objects to hold the matrices
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
  matrixqc(0) = totalutrans; // summed number of U transitions
  matrixqc(1) = totalftrans; // summed number of F transitions
  matrixqc(2) = totalmatrices;
  
  // Final output
  
  List output = List::create(Named("A") = A, _["U"] = U, _["F"] = F, 
    _["hstages"] = NULL, _["agestages"] = agestages, _["ahstages"] = stages,
    _["labels"] = cheatsheet, _["matrixqc"] = matrixqc);
  
  return output;
}

//' Full Eigen Analysis of a Single Dense Matrix
//' 
//' Function \code{decomp3()} returns all eigenvalues, right eigenvectors, and
//' left eigenvectors estimated for a matrix by the \code{eig_gen}() function
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
// [[Rcpp::export(.decomp3)]]
List decomp3(arma::mat Amat) {
  
  arma::cx_vec Aeigval;
  arma::cx_mat Aeigvecl;
  arma::cx_mat Aeigvecr;
  
  eig_gen(Aeigval, Aeigvecl, Aeigvecr, Amat);
  
  List output = List::create(Named("eigenvalues") = Aeigval,
    _["left_eigenvectors"] = Aeigvecl, _["right_eigenvectors"] = Aeigvecr);
  
  return output;
}

//' Full Eigen Analysis of a Single Sparse Matrix
//' 
//' Function \code{decomp3sp()} returns all eigenvalues, right eigenvectors, and
//' left eigenvectors estimated for a matrix by the \code{eigs_gen}() function
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
// [[Rcpp::export(.decomp3sp)]]
List decomp3sp(arma::mat Amat) {
  
  arma::sp_mat spAmat(Amat);
  arma::sp_mat t_spAmat = spAmat.t();
  
  arma::cx_vec Aeigval;
  arma::cx_vec Aeigvall;
  arma::cx_mat Aeigvecl;
  arma::cx_mat Aeigvecr;
  
  eigs_gen(Aeigval, Aeigvecr, spAmat, 1); // Identified as problem with valgrind
  
  eigs_gen(Aeigvall, Aeigvecl, t_spAmat, 1);
  
  List output = List::create(Named("eigenvalues") = Aeigval,
    _["left_eigenvectors"] = Aeigvecl, _["right_eigenvectors"] = Aeigvecr);
  
  return output;
}

//' Full Eigen Analysis of a Single Sparse Matrix, with Sparse Input
//' 
//' \code{decomp3sp_inp()} returns all eigenvalues, right eigenvectors, and left
//' eigenvectors estimated for a matrix by the \code{eigs_gen}() function
//' in the C++ Armadillo library. Works with sparse matrices.
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//'
//' @return This function returns all estimated eigenvalues, right
//' eigenvectors, and left eigenvectors of a single matrix. This output is
//' provided as a list with three parts, named appropriately.
//' 
//' @section Notes:
//' This function works slightly differently from function \code{decomp3sp()} in
//' that the latter function requires a sparse matrix provided in dense format,
//' while this function requires a sparse matrix in sparse format.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.decomp3sp_inp)]]
List decomp3sp_inp(arma::sp_mat spAmat) {
  
  arma::sp_mat t_spAmat = spAmat.t();
  
  arma::cx_vec Aeigval;
  arma::cx_vec Aeigvall;
  arma::cx_mat Aeigvecl;
  arma::cx_mat Aeigvecr;
  
  eigs_gen(Aeigval, Aeigvecr, spAmat, 1);
  
  eigs_gen(Aeigvall, Aeigvecl, t_spAmat, 1);
  
  List output = List::create(Named("eigenvalues") = Aeigval,
    _["left_eigenvectors"] = Aeigvecl, _["right_eigenvectors"] = Aeigvecr);
  
  return output;
}

//' Estimate Deterministic Population Growth Rate of Any Matrix
//' 
//' \code{lambda3matrix()} returns the dominant eigenvalue of a single
//' dense or sparse projection matrix, provided in dense matrix format.
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' @param sparse A logical value indicating whether to use sparse matrix
//' format.
//'
//' @return This function returns the dominant eigenvalue of the matrix. This
//' is given as the largest real part of all eigenvalues estimated via the 
//' \code{eig_gen}() and \code{eigs_gen}() functions in the C++ Armadillo
//' library.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.lambda3matrix)]]
double lambda3matrix(arma::mat Amat, bool sparse) {
  
  double lambda {0};
  
  if (!sparse) {
    List eigenstuff = decomp3(Amat);
    
    arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
    
    lambda = max(realeigenvals);
  } else {
    List eigenstuff = decomp3sp(Amat);
    
    arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
    
    lambda = max(realeigenvals);
  }
  
  return lambda;
}

//' Estimate Stable Stage Distribution of Any Population Matrix
//' 
//' \code{ss3matrix()} returns the stable stage distribution for a 
//' dense or sparse population matrix.
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' @param sparse A logical value indicating whether to use sparse or dense
//' format in matrix calculations.
//' 
//' @return This function returns the stable stage distribution corresponding to
//' the input matrix.
//' 
//' @seealso \code{\link{stablestage3}()}
//' @seealso \code{\link{stablestage3.lefkoMat}()}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.ss3matrix)]]
arma::vec ss3matrix(arma::mat Amat, bool sparse) {
  List eigenstuff;
  
  if (sparse) {
    eigenstuff = decomp3sp(Amat);
  } else {
    eigenstuff = decomp3(Amat);
  }
  
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  
  int lambda1 = realeigenvals.index_max();
  
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.0000000001); // This line replaces all numbers lower than 1 x 10-10 with 0
  
  double rvsum = sum(realrightvec);
  realrightvec = realrightvec / rvsum;
  
  return realrightvec;
}

//' Estimate Reproductive Value of Any Population Matrix
//' 
//' \code{rv3matrix()} returns the reproductive values for stages in a
//' dense or sparse population matrix (both provided in dense matrix format).
//' The function provides standard reproductive values, meaning that the overall
//' reproductive values of basic life history stages in a historical matrix are
//' not provided (the \code{\link{repvalue3.lefkoMat}()} function estimates
//' these on the basis of stage description information provided in the
//' \code{lefkoMat} object used as input in that function).
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' @param sparse A logical value indicating whether to use sparse or dense
//' format in matrix calculations.
//' 
//' @return This function returns a vector characterizing the reproductive
//' values for stages of a population projection matrix.
//' 
//' @seealso \code{\link{repvalue3}()}
//' @seealso \code{\link{repvalue3.lefkoMat}()}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.rv3matrix)]]
arma::vec rv3matrix(arma::mat Amat, bool sparse) {
  List eigenstuff;
  
  if (sparse) {
    eigenstuff = decomp3sp(Amat);
  } else {
    eigenstuff = decomp3(Amat);
  }
  
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  
  int lambda1 = realeigenvals.index_max();
  
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.0000000001); // This line replaces all numbers lower than 1 x 10-10 with 0

  // This section identifies the first non-zero element of the reproductive value vector
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  
  realleftvec = realleftvec / rlvmin;
  
  return realleftvec;
}

//' Estimate Deterministic Sensitivities of Any Population Matrix
//' 
//' \code{sens3matrix()} returns the sensitivity of lambda with respect
//' to each element in a dense or sparse matrix (provided in dense matrix
//' format). This is accomplished via the \code{eig_gen}() and \code{eigs_gen}()
//' functions in the C++ Armadillo library.
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' @param sparse A logical value indicating whether to use sparse or dense
//' format in matrix calculations.
//' 
//' @return This function returns a matrix of deterministic sensitivities. 
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.sens3matrix)]]
arma::mat sens3matrix(arma::mat Amat, bool sparse) {
  List eigenstuff;
  
  if (sparse) {
    eigenstuff = decomp3sp(Amat);
  } else {
    eigenstuff = decomp3(Amat);
  }
  
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  
  int lambda1 = realeigenvals.index_max();
  
  // This is the w vector
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-10 with 0
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  realrightvec = realrightvec / rvsum;
  
  // This is the v vector
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-10 with 0

  // This section identifies the first non-zero element of the reproductive value vector
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  realleftvec = realleftvec / rlvmin;
  
  arma::vec vwprod (rvel);
  arma::mat smat (rvel, rvel);
  smat.zeros();
  
  // This loop and the following line create the scalar product vw
  for (int i = 0; i < rvel; i++) {
    vwprod(i) = realrightvec(i) * realleftvec(i);
  }
  double vwscalar = sum(vwprod);
  
  // This loop populates the sensitivity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      smat(i, j) = realleftvec(i) * realrightvec(j) / vwscalar;
    }
  }
  
  return smat;
}

//' Estimate Deterministic Sensitivities of A Spars Matrixe
//' 
//' \code{sens3sp_matrix()} returns the sensitivity of lambda with respect
//' to each element in a sparse matrix, provided in sparse matrix format. This
//' is accomplished via the \code{eigs_gen}() function in the C++ Armadillo
//' library.
//' 
//' @param Aspmat A population projection matrix in sparse matrix format.
//' @param refmat A sparse matrix used for reference to create associated 0s in
//' the sensitivity matrix.
//' 
//' @return This function returns a sparse matrix of deterministic
//' sensitivities. Zeroes are derived from the reference matrix, and replace
//' non-zero entries that will be zeroed out in the following math. Currently
//' used in LTRE estimation.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.sens3sp_matrix)]]
arma::sp_mat sens3sp_matrix(arma::sp_mat Aspmat, arma::sp_mat refmat) {
  
  List eigenstuff = decomp3sp_inp(Aspmat);

  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  
  int lambda1 = realeigenvals.index_max();
  
  // This is the w vector
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-10 with 0
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  realrightvec = realrightvec / rvsum;
  
  // This is the v vector
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.00000000000001); // This line replaces all numbers lower than 1 x 10-10 with 0

  // This section identifies the first non-zero element of the reproductive value vector
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  realleftvec = realleftvec / rlvmin;
  
  arma::vec vwprod (rvel);
  arma::sp_mat smat (rvel, rvel);
  smat.zeros();
  
  // This loop and the following line create the scalar product vw
  for (int i = 0; i < rvel; i++) {
    vwprod(i) = realrightvec(i) * realleftvec(i);
  }
  double vwscalar = sum(vwprod);
  
  // This loop populates the sensitivity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      if (refmat(i, j) != 0) {
        smat(i, j) = realleftvec(i) * realrightvec(j) / vwscalar;
      }
    }
  }
  
  return smat;
}

//' Estimate Deterministic Sensitivities of a Historical LefkoMat Object
//' 
//' \code{sens3hlefko()} returns the sensitivity of lambda with respect
//' to each historical stage-pair in the matrix, and the associated
//' sensitivity for each life history stage. This is accomplished via the 
//' \code{eigs_gen}() function in the C++ Armadillo library.
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' @param ahstages An integar vector of unique ahistorical stages.
//' @param hstages An integar vector of unique historical stage pairs.
//' 
//' @return This function returns a list with two deterministic sensitivity
//' matrices:
//' \item{h_smat}{Matrix of sensitivities corresponding to the historical
//' matrix.}
//' \item{ah_smat}{Matrix of sensitivities corresponding to the ahistorical
//' matrix.}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.sens3hlefko)]]
List sens3hlefko(arma::mat Amat, DataFrame ahstages, DataFrame hstages) {
  
  arma::uvec stage_id = ahstages["stage_id"];
  arma::uvec h_stage_2 = hstages["stage_id_2"];
  arma::uvec h_stage_1 = hstages["stage_id_1"];
  
  List eigenstuff = decomp3sp(Amat);
  
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  
  int lambda1 = realeigenvals.index_max();
  
  // This is the w vector
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.00000000000001); // Using a lower threshold than previously here
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  realrightvec = realrightvec / rvsum;
  
  // This is the v vector
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.00000000000001); // Using a lower threshold than previously here

  // This section identifies the first non-zero element of the reproductive value vector
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  realleftvec = realleftvec / rlvmin;
  
  int ahstagelength = stage_id.n_elem;
  
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
  
  // This loop and the following line create the scalar product vw and the ahistorical stable stage distribution w
  for (int i = 0; i < rvel; i++) {
    vwprod(i) = realrightvec(i) * realleftvec(i);
    
    ahrows = h_stage_2(i) - 1;
    
    wcorrah(ahrows) = wcorrah(ahrows) + realrightvec(i);
  }
  double vwscalar = sum(vwprod);
  
  // This loop creates a corrected reproductive value vector
  for (int i = 0; i < rvel; i++) {
    ahrows = h_stage_2(i) - 1;
    
    if (wcorrah(ahrows) != 0) {
      vcorrah(ahrows) = vwprod(i) / wcorrah(ahrows) + vcorrah(ahrows);
    } else {
      // This line deals with the fact that some stages are associated with expected corrected stable stage proportions of 0
      vcorrah(ahrows) = 0 + vcorrah(ahrows);
    }
  }
  
  // This loop populates the sensitivity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      smat(i, j) = realleftvec(i) * realrightvec(j) / vwscalar;
    }
  }
  
  // This next section creates the ahistorical sensitivity matrix
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

//' Estimate Deterministic Elasticities of Any Population Matrix
//' 
//' \code{elas3matrix()} returns the elasticity of lambda with respect
//' to each element in a dense or sparse matrix, both provided in dense matrix
//' format. This is accomplished via the \code{eig_gen}() and \code{eigs_gen}()
//' functions in the C++ Armadillo library.
//' 
//' @param Amat A population projection matrix of class \code{matrix}.
//' @param sparse A logical value indicating whether to use sparse or dense
//' format in matrix calculations.
//' 
//' @return This function returns a matrix of deterministic elasticities. 
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.elas3matrix)]]
arma::mat elas3matrix(arma::mat Amat, bool sparse) {
  List eigenstuff;
  
  if (sparse) {
    eigenstuff = decomp3sp(Amat);
  } else {
    eigenstuff = decomp3(Amat);
  }
  
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  
  int lambda1 = realeigenvals.index_max();
  double lambda = max(realeigenvals);

  // This is the w vector
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.00000000000001); // Lower threshold than used in w and v
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  realrightvec = realrightvec / rvsum;

  // This is the v vector
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.00000000000001); // Lower threshold than used in w and v

  // This section identifies the first non-zero element of the reproductive value vector
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  realleftvec = realleftvec / rlvmin;

  arma::vec vwprod (rvel);
  arma::mat emat (rvel, rvel);
  emat.zeros();
  
  // This loop and the following line create the scalar product vw
  for (int i = 0; i < rvel; i++) {
    vwprod(i) = realrightvec(i) * realleftvec(i);
  }
  double vwscalar = sum(vwprod);
  
  // This loop populates the elasticity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      emat(i, j) = (realleftvec(i) * realrightvec(j) * Amat(i, j)) / (vwscalar * lambda);
    }
  }
  
  return emat;
}

//' Estimate Deterministic Elasticities of a Historical LefkoMat Object
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
//' @return This function returns a list with two deterministic elasticity
//' matrices:
//' \item{h_emat}{Matrix of elasticities corresponding to the historical matrix.}
//' \item{ah_emat}{Matrix of elasticities corresponding to the ahistorical
//' matrix, but using summed historical elasticities as the basis of estimation.}
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.elas3hlefko)]]
List elas3hlefko(arma::mat Amat, DataFrame ahstages, DataFrame hstages) {
  
  arma::uvec stage_id = ahstages["stage_id"];
  arma::uvec h_stage_2 = hstages["stage_id_2"];
  arma::uvec h_stage_1 = hstages["stage_id_1"];
  
  List eigenstuff = decomp3sp(Amat);
  
  arma::vec realeigenvals = real(as<arma::cx_vec>(eigenstuff["eigenvalues"]));
  int lambda1 = realeigenvals.index_max();
  double lambda = max(realeigenvals);
  
  // This is the w vector
  arma::vec realrightvec = real(as<arma::cx_mat>(eigenstuff["right_eigenvectors"]).col(lambda1));
  realrightvec.clean(0.00000000000001);
  
  double rvsum = sum(realrightvec);
  int rvel = realrightvec.n_elem;
  realrightvec = realrightvec / rvsum;

  // This is the v vector
  arma::vec realleftvec = real(as<arma::cx_mat>(eigenstuff["left_eigenvectors"]).col(lambda1));
  realleftvec.clean(0.00000000000001);

  // This section identifies the first non-zero element of the reproductive value vector
  arma::uvec rlvabsalt = find(realleftvec);
  double rlvmin = realleftvec(static_cast<unsigned long long>(rlvabsalt(0)));
  realleftvec = realleftvec / rlvmin;

  arma::vec vwprod (rvel);
  arma::mat emat (rvel, rvel);
  emat.zeros();
  
  // This loop and the following line create the scalar product vw
  for (int i = 0; i < rvel; i++) {
    vwprod(i) = realrightvec(i) * realleftvec(i);
  }
  double vwscalar = sum(vwprod);
  
  // The next few lines set up the empty ahistorical matrix
  int ahstagelength = stage_id.n_elem;

  arma::mat ahelas(ahstagelength, ahstagelength);
  ahelas.zeros();
  
  // This loop populates the elasticity matrix
  for (int i = 0; i < rvel; i++) {
    for (int j = 0; j < rvel; j++) {
      
      emat(i, j) = (realleftvec(i) * realrightvec(j) * Amat(i, j)) / (vwscalar * lambda);
      
      ahelas((h_stage_2(i) - 1), (h_stage_1(i) - 1)) = ahelas((h_stage_2(i) - 1), (h_stage_1(i) - 1)) + emat(i, j);
    }
  }
  
  List output = List::create(Named("h_emat") = emat, _["ah_emat"] = ahelas);
  
  return output;
}

//' Core Time-based Population Matrix Projection Function
//' 
//' Function \code{proj3()} runs the matrix projections used in other functions
//' in package \code{lefko3}.
//' 
//' @param start_vec The starting population vector for the projection.
//' @param core_list A list of full projection matrices, corresponding to the 
//' \code{$A} list within a \code{lefkoMat} object.
//' @param mat_order A vector giving the order of matrices to use at each occasion.
//' @param standardize A logical value stating whether to standardize population
//' size vector to sum to 1 at each estimated occasion.
//' @param growthonly A logical value stating whether to output only a matrix
//' showing the change in population size from one year to the next for use in
//' stochastic population growth rate estimation (TRUE), or a larger matrix also
//' containing the w and v projections for stochastic perturbation analysis,
//' stage distribution estimation, and reproductive value estimation.
//' @param integeronly A logical value indicating whether to round all projected
//' numbers of individuals to the nearest integer.
//' 
//' @return A matrix in which, if \code{growthonly = TRUE}, each row is the
//' population vector at each projected occasion, and if \code{growthonly =
//' FALSE}, the top third of the matrix is the actual number of individuals in
//' each stage across time, the second third is the w projection (stage
//' distribution), and the bottom third is the v projection (reproductive
//' values) for use in estimation of stochastic sensitivities and elasticities
//' (in addition, a further row is appended to the bottom, corresponding to the
//' \emph{R} vector, which is the sum of the unstandardized \emph{w} vector
//' resulting from each occasion's projection).
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.proj3)]]
arma::mat proj3(arma::vec start_vec, List core_list, arma::uvec mat_order,
  bool standardize, bool growthonly, bool integeronly) {
  
  int sparse_switch {0};
  
  int nostages = start_vec.n_elem;
  int theclairvoyant = mat_order.n_elem;
  arma::vec theseventhson;
  arma::rowvec theseventhgrandson;
  arma::mat theprophecy;
  arma::mat thesecondprophecy;
  
  arma::mat popproj(nostages, (theclairvoyant + 1)); // This is the population vector
  arma::mat wpopproj(nostages, (theclairvoyant + 1)); // This is the population w vector
  arma::mat vpopproj(nostages, (theclairvoyant + 1)); // This is the population v vector
  arma::mat Rvecmat(1, (theclairvoyant+1));
  popproj.zeros();
  wpopproj.zeros();
  vpopproj.zeros();
  Rvecmat.zeros();
  
  theseventhson = start_vec;
  theseventhgrandson = start_vec.as_row();
  
  arma::mat finaloutput;
  
  // Here we will check if the matrix is large and sparse
  int test_elems = as<arma::mat>(core_list(0)).n_elem;
  arma::uvec nonzero_elems = find(as<arma::mat>(core_list(0)));
  int all_nonzeros = nonzero_elems.n_elem;
  double sparse_check = static_cast<double>(all_nonzeros) / static_cast<double>(test_elems);
  if (sparse_check <= 0.5) {
    sparse_switch = 1;
  } else sparse_switch = 0;
  
  // Now the projection
  popproj.col(0) = start_vec;
  if (!growthonly) {
    wpopproj.col(0) = start_vec / sum(start_vec);
    vpopproj.col(theclairvoyant) = start_vec / sum(start_vec);
    Rvecmat(0) = sum(start_vec);
  }
  
  if (sparse_switch == 0) {
    // Dense matrix projection
    
    for (int i = 0; i < theclairvoyant; i++) {
      theprophecy = as<arma::mat>(core_list[(mat_order(i))]);
      
      theseventhson = theprophecy * theseventhson; // thechosenone
      if (integeronly) {
        theseventhson = floor(theseventhson);
      }
      popproj.col(i+1) = theseventhson;
      
      Rvecmat(i+1) = sum(theseventhson);
      
      if (standardize) {
        theseventhson = theseventhson / sum(theseventhson);
      }
      
      if (!growthonly) {
        wpopproj.col(i+1) = theseventhson / Rvecmat(i+1);
        
        thesecondprophecy = as<arma::mat>(core_list[(mat_order(theclairvoyant - (i+1)))]);
        theseventhgrandson = theseventhgrandson * thesecondprophecy;
        
        double seventhgrandsum = sum(theseventhgrandson);
        arma::vec midwife = theseventhgrandson.as_col() / seventhgrandsum;
        
        theseventhgrandson = theseventhgrandson / seventhgrandsum;
        
        vpopproj.col(theclairvoyant - (i+1)) = midwife;
      }
    }
  } else {
    // Sparse matrix projection
    
    arma::sp_mat sparse_seventhson = arma::sp_mat(theseventhson);
    
    int matlist_length = core_list.size();
    arma::mat first_mat = core_list(0);
    arma::sp_mat new_sparse = arma::sp_mat(first_mat);
    Rcpp::List sparse_list = List::create(_["1"] = new_sparse);
    if(matlist_length > 1) {
      for (int i = 1; i < matlist_length; i++) {
        first_mat = as<arma::mat>(core_list(i));
        new_sparse = arma::sp_mat(first_mat);
        sparse_list.push_back(new_sparse);
      }
    }
    arma::sp_mat sparse_prophecy;
    arma::sp_mat sparse_secondprophecy;
    
    for (int i = 0; i < theclairvoyant; i++) {
      sparse_prophecy = as<arma::sp_mat>(sparse_list[(mat_order(i))]);
      
      sparse_seventhson = sparse_prophecy * sparse_seventhson;
      if (integeronly) {
        sparse_seventhson = floor(sparse_seventhson);
      }
      popproj.col(i+1) = arma::vec(arma::mat(sparse_seventhson));
      
      Rvecmat(i+1) = sum(popproj.col(i+1));
      
      if (standardize) {
        sparse_seventhson = sparse_seventhson / sum(popproj.col(i+1));
      }
      
      if (!growthonly) {
        wpopproj.col(i+1) = arma::vec(arma::mat(sparse_seventhson)) / Rvecmat(i+1);
        
        sparse_secondprophecy = as<arma::sp_mat>(sparse_list[(mat_order(theclairvoyant - (i+1)))]);
        theseventhgrandson = theseventhgrandson * sparse_secondprophecy;
        
        double seventhgrandsum = sum(theseventhgrandson);
        arma::vec midwife = theseventhgrandson.as_col() / seventhgrandsum;
        
        theseventhgrandson = theseventhgrandson / seventhgrandsum;
        
        vpopproj.col(theclairvoyant - (i+1)) = midwife;
      }
    }
  }
  
  if (growthonly) {
    return popproj;
  } else {
    arma::mat revised_vproj = join_cols(vpopproj, Rvecmat);
    arma::mat expanded_proj = join_cols(wpopproj, revised_vproj);
    
    return join_cols(popproj, expanded_proj);
  }
}

//' Slimmed-down Time-based Population Sparse Matrix Projection Function
//' 
//' Function \code{proj3sp()} runs the matrix projections used in some other
//' functions in package \code{lefko3}, but only when the input is sparse. This
//' is a slimmed down version of function \code{proj3()}
//' 
//' @param start_vec The starting population vector for the projection.
//' @param core_list A list of full projection matrices, corresponding to
//' the \code{$A} list within a \code{lefkoMat} object. Matrices must be in
//' \code{arma::sp_mat} format.
//' @param mat_order A vector giving the order of matrices to use at each occasion.
//' @param standardize A logical value stating whether to standardize population
//' size vector to sum to 1 at each estimated occasion.
//' @param growthonly A logical value stating whether to output only a matrix
//' showing the change in population size from one year to the next for use in
//' stochastic population growth rate estimation (TRUE), or a larger matrix also
//' containing the w and v projections for stochastic perturbation analysis,
//' stage distribution estimation, and reproductive value estimation.
//' @param integeronly A logical value indicating whether to round all projected
//' numbers of individuals to the nearest integer.
//' 
//' @return A matrix in which, if \code{growthonly = TRUE}, each row is the
//' population vector at each projected occasion, and if \code{growthonly =
//' FALSE}, the top third of the matrix is the actual number of individuals in
//' each stage across time, the second third is the w projection (stage
//' distribution), and the bottom third is the v projection (reproductive
//' values) for use in estimation of stochastic sensitivities and elasticities
//' (in addition, a further row is appended to the bottom, corresponding to the
//' \emph{R} vector, which is the sum of the unstandardized \emph{w} vector
//' resulting from each occasion's projection).
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.proj3sp)]]
arma::mat proj3sp(arma::vec start_vec, List core_list, arma::uvec mat_order,
  bool standardize, bool growthonly, bool integeronly) {
  
  int nostages = start_vec.n_elem;
  int theclairvoyant = mat_order.n_elem;
  arma::vec theseventhson;
  arma::rowvec theseventhgrandson;
  arma::mat theprophecy;
  arma::mat thesecondprophecy;
  
  arma::mat popproj(nostages, (theclairvoyant + 1)); // This is the population vector
  arma::mat wpopproj(nostages, (theclairvoyant + 1)); // This is the population w vector
  arma::mat vpopproj(nostages, (theclairvoyant + 1)); // This is the population v vector
  arma::mat Rvecmat(1, (theclairvoyant+1));
  popproj.zeros();
  wpopproj.zeros();
  vpopproj.zeros();
  Rvecmat.zeros();
  
  theseventhson = start_vec;
  theseventhgrandson = start_vec.as_row();
  
  arma::sp_mat sparse_seventhson = arma::sp_mat(theseventhson);
    
  arma::mat finaloutput;
  
  // Now the projection
  popproj.col(0) = start_vec;
  if (!growthonly) {
    wpopproj.col(0) = start_vec / sum(start_vec);
    vpopproj.col(theclairvoyant) = start_vec / sum(start_vec);
    Rvecmat(0) = sum(start_vec);
  }
  
  // Sparse matrix projection
  arma::sp_mat sparse_prophecy;
  arma::sp_mat sparse_secondprophecy;
    
  for (int i = 0; i < theclairvoyant; i++) {
    sparse_prophecy = as<arma::sp_mat>(core_list[(mat_order(i))]);
      
    sparse_seventhson = sparse_prophecy * sparse_seventhson;
    if (integeronly) {
      sparse_seventhson = floor(sparse_seventhson);
    }
    popproj.col(i+1) = arma::vec(arma::mat(sparse_seventhson));
    
    Rvecmat(i+1) = sum(popproj.col(i+1));
    
    if (standardize) {
      sparse_seventhson = sparse_seventhson / sum(popproj.col(i+1));
    }
    
    if (!growthonly) {
      wpopproj.col(i+1) = arma::vec(arma::mat(sparse_seventhson)) / Rvecmat(i+1);
      
      sparse_secondprophecy = as<arma::sp_mat>(core_list[(mat_order(theclairvoyant - (i+1)))]);
      theseventhgrandson = theseventhgrandson * sparse_secondprophecy;
      
      double seventhgrandsum = sum(theseventhgrandson);
      arma::vec midwife = theseventhgrandson.as_col() / seventhgrandsum;
      
      theseventhgrandson = theseventhgrandson / seventhgrandsum;
      
      vpopproj.col(theclairvoyant - (i+1)) = midwife;
    }
  }
  
  if (growthonly) {
    return popproj;
  } else {
    arma::mat revised_vproj = join_cols(vpopproj, Rvecmat);
    arma::mat expanded_proj = join_cols(wpopproj, revised_vproj);
    
    return join_cols(popproj, expanded_proj);
  }
}

//' Core Time-based Density-Dependent Population Matrix Projection Function
//' 
//' Function \code{proj3dens()} runs density-dependent matrix projections.
//' 
//' @param start_vec The starting population vector for the projection.
//' @param core_list A list of full projection matrices, corresponding to the 
//' \code{A} list within a \code{lefkoMat} object.
//' @param mat_order A vector giving the order of matrices to use at each occasion.
//' @param growthonly A logical value stating whether to output only a matrix
//' showing the change in population size from one year to the next for use in
//' stochastic population growth rate estimation (TRUE), or a larger matrix also
//' containing the w and v projections for stochastic perturbation analysis,
//' stage distribution estimation, and reproductive value estimation.
//' @param integeronly A logical value indicating whether to round all projected
//' numbers of individuals to the nearest integer.
//' @param substoch An integer value indicating whether to force survival-
//' transition matrices to be substochastic in density dependent simulations.
//' Defaults to \code{0}, which does not force substochasticity. Alternatively,
//' \code{1} forces all survival-transition elements to range from 0.0 to 1.0,
//' and \code{2} forces all column rows to total no more than 1.0.
//' @param dens_input The original \code{lefkoDens} data frame supplied through
//' the \code{\link{density_input}()} function.
//' @param dens_index A list giving the indices of elements in object
//' \code{dens_input}.
//' 
//' @return A matrix in which, if \code{growthonly = TRUE}, each row is the
//' population vector at each projected occasion, and if \code{growthonly =
//' FALSE}, the top third of the matrix is the actual number of individuals in
//' each stage across time, the second third is the w projection (stage
//' distribution), and the bottom third is the v projection (reproductive
//' values) for use in estimation of stochastic sensitivities and elasticities
//' (in addition, a further row is appended to the bottom, corresponding to the
//' \emph{R} vector, which is the sum of the unstandardized \emph{w} vector
//' resulting from each occasion's projection).
//' 
//' @section Notes:
//' There is no option to standardize population vectors here, because density
//' dependence requires the full population size to be tracked.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.proj3dens)]]
arma::mat proj3dens(arma::vec start_vec, List core_list, arma::uvec mat_order,
  bool growthonly, bool integeronly, int substoch, Rcpp::DataFrame dens_input,
  Rcpp::List dens_index) {
  
  int sparse_switch {0};
  int time_delay {1};
  double pop_size {0};
  
  int nostages = start_vec.n_elem;
  int theclairvoyant = mat_order.n_elem;
  arma::vec theseventhson;
  arma::rowvec theseventhgrandson;
  arma::mat theprophecy;
  arma::mat thesecondprophecy;
  
  // Density dependence arguments
  arma::uvec dyn_index321 = dens_index["index321"];
  arma::uvec dyn_index_col = dens_index[1];
  arma::uvec dyn_style = dens_input["style"];
  arma::vec dyn_alpha = dens_input["alpha"];
  arma::vec dyn_beta = dens_input["beta"];
  arma::uvec dyn_delay = dens_input["time_delay"];
  arma::uvec dyn_type = dens_input["type"];
  int n_dyn_elems = dyn_index321.n_elem;
  
  // Matrices and vectors for projection results
  arma::mat popproj(nostages, (theclairvoyant + 1)); // This is the population vector
  arma::mat wpopproj(nostages, (theclairvoyant + 1)); // This is the population w vector
  arma::mat vpopproj(nostages, (theclairvoyant + 1)); // This is the population v vector
  arma::mat Rvecmat(1, (theclairvoyant+1));
  popproj.zeros();
  wpopproj.zeros();
  vpopproj.zeros();
  Rvecmat.zeros();
  
  theseventhson = start_vec;
  theseventhgrandson = start_vec.as_row();
  
  arma::mat finaloutput;
  
  // Here we will check if the matrix is large and sparse
  int test_elems = as<arma::mat>(core_list(0)).n_elem;
  arma::uvec nonzero_elems = find(as<arma::mat>(core_list(0)));
  int all_nonzeros = nonzero_elems.n_elem;
  double sparse_check = static_cast<double>(all_nonzeros) / static_cast<double>(test_elems);
  if (sparse_check <= 0.5) {
    sparse_switch = 1;
  } else sparse_switch = 0;
  
  // Now the projection
  popproj.col(0) = start_vec;
  if (!growthonly) {
    wpopproj.col(0) = start_vec / sum(start_vec);
    vpopproj.col(theclairvoyant) = start_vec / sum(start_vec);
    Rvecmat(0) = sum(start_vec);
  }
  
  double changing_element {0.0};
  double changing_colsum {0.0};
  
  if (sparse_switch == 0) {
    // Dense matrix projection
    
    for (int i = 0; i < theclairvoyant; i++) {
      theprophecy = as<arma::mat>(core_list[(mat_order(i))]);
      
      // Now we modify the matrix with density dependence
      for (int j = 0; j < n_dyn_elems; j++) {
        time_delay = dyn_delay(j);
        if (time_delay > 0) time_delay = time_delay - 1;
        
        if (i >= time_delay) {
          pop_size = sum(popproj.col(i - time_delay));
          
          if (dyn_style(j) == 1) { // Ricker
            changing_element = theprophecy(dyn_index321(j)) * 
              dyn_alpha(j) * exp((-1*dyn_beta(j)) * pop_size); // Fi*ALPHA*exp(-BETA*n)
            
            if (substoch == 0 || dyn_type(j) == 2) {
              theprophecy(dyn_index321(j)) = changing_element;
            } else if (substoch == 1 && dyn_type(j) == 1) {
              if (changing_element > 1.0) {
                changing_element = 1.0;
              } else if (changing_element < 0.0) {
                changing_element = 0.0;
              }
              theprophecy(dyn_index321(j)) = changing_element;
            } else if (substoch == 2 && dyn_type(j) == 1) {
              changing_colsum = sum(theprophecy.col(dyn_index_col(j))) - theprophecy(dyn_index321(j));
              
              if (changing_element > (1.0 - changing_colsum)) {
                changing_element = (1.0 - changing_colsum);
              }
              theprophecy(dyn_index321(j)) = changing_element;
            }
          } else if (dyn_style(j) == 2) { // Beverton-Holt
            changing_element = theprophecy(dyn_index321(j)) * 
              dyn_alpha(j) / (1 + dyn_beta(j) * pop_size); // Fi*ALPHA/(1+BETA*n)
            
            if (substoch == 0 || dyn_type(j) == 2) {
              theprophecy(dyn_index321(j)) = changing_element;
            } else if (substoch == 1 && dyn_type(j) == 1) {
              if (changing_element > 1.0) {
                changing_element = 1.0;
              } else if (changing_element < 0.0) {
                changing_element = 0.0;
              }
              theprophecy(dyn_index321(j)) = changing_element;
            } else if (substoch == 2 && dyn_type(j) == 1) {
              changing_colsum = sum(theprophecy.col(dyn_index_col(j))) - theprophecy(dyn_index321(j));
              
              if (changing_element > (1.0 - changing_colsum)) {
                changing_element = (1.0 - changing_colsum);
              }
              theprophecy(dyn_index321(j)) = changing_element;
            }
          } else if (dyn_style(j) == 3) { // Usher function
            changing_element = theprophecy(dyn_index321(j)) * 
              (1 / (1 + exp(dyn_alpha(j) * pop_size + dyn_beta(j)))); // Fi*(1 / (1 + exp(alpha*N+b)))
            
            if (substoch == 0 || dyn_type(j) == 2) {
              theprophecy(dyn_index321(j)) = changing_element;
            } else if (substoch == 1 && dyn_type(j) == 1) {
              if (changing_element > 1.0) {
                changing_element = 1.0;
              } else if (changing_element < 0.0) {
                changing_element = 0.0;
              }
              theprophecy(dyn_index321(j)) = changing_element;
            } else if (substoch == 2 && dyn_type(j) == 1) {
              changing_colsum = sum(theprophecy.col(dyn_index_col(j))) - theprophecy(dyn_index321(j));
              
              if (changing_element > (1.0 - changing_colsum)) {
                changing_element = (1.0 - changing_colsum);
              }
              theprophecy(dyn_index321(j)) = changing_element;
            }
          } else if (dyn_style(j) == 4) { // Logistic function
            changing_element = theprophecy(dyn_index321(j)) * 
              (1 - pop_size / dyn_alpha(j)); // Fi*(1 - ALPHA/n)
            
            if (substoch == 0 || dyn_type(j) == 2) {
              theprophecy(dyn_index321(j)) = changing_element;
            } else if (substoch == 1 && dyn_type(j) == 1) {
              if (changing_element > 1.0) {
                changing_element = 1.0;
              } else if (changing_element < 0.0) {
                changing_element = 0.0;
              }
              theprophecy(dyn_index321(j)) = changing_element;
            } else if (substoch == 2 && dyn_type(j) == 1) {
              changing_colsum = sum(theprophecy.col(dyn_index_col(j))) - theprophecy(dyn_index321(j));
              
              if (changing_element > (1.0 - changing_colsum)) {
                changing_element = (1.0 - changing_colsum);
              }
              theprophecy(dyn_index321(j)) = changing_element;
            }
          }
        }
      }
      
      theseventhson = theprophecy * theseventhson; // thechosenone
      if (integeronly) {
        theseventhson = floor(theseventhson);
      }
      popproj.col(i+1) = theseventhson;
      
      Rvecmat(i+1) = sum(theseventhson);
      
      if (!growthonly) {
        wpopproj.col(i+1) = theseventhson / Rvecmat(i+1);
        
        thesecondprophecy = as<arma::mat>(core_list[(mat_order(theclairvoyant - (i+1)))]);
        theseventhgrandson = theseventhgrandson * thesecondprophecy;
        
        double seventhgrandsum = sum(theseventhgrandson);
        arma::vec midwife = theseventhgrandson.as_col() / seventhgrandsum;
        
        theseventhgrandson = theseventhgrandson / seventhgrandsum;
        
        vpopproj.col(theclairvoyant - (i+1)) = midwife;
      }
    }
  } else {
    // Sparse matrix projection
    
    arma::sp_mat sparse_seventhson = arma::sp_mat(theseventhson);
    
    int matlist_length = core_list.size();
    arma::mat first_mat = core_list(0);
    arma::sp_mat new_sparse = arma::sp_mat(first_mat);
    Rcpp::List sparse_list = List::create(_["1"] = new_sparse);
    if(matlist_length > 1) {
      for (int i = 1; i < matlist_length; i++) {
        first_mat = as<arma::mat>(core_list(i));
        new_sparse = arma::sp_mat(first_mat);
        sparse_list.push_back(new_sparse);
      }
    }
    arma::sp_mat sparse_prophecy;
    arma::sp_mat sparse_secondprophecy;
    
    for (int i = 0; i < theclairvoyant; i++) {
      sparse_prophecy = as<arma::sp_mat>(sparse_list[(mat_order(i))]);
      
      // Now we modify the matrix with density dependence
      for (int j = 0; j < n_dyn_elems; j++) {
        time_delay = dyn_delay(j);
        if (time_delay > 0) time_delay = time_delay - 1;
        
        if (i >= time_delay) {
          pop_size = sum(popproj.col(i - time_delay));
          
          if (dyn_style(j) == 1) { // Ricker
            changing_element = sparse_prophecy(dyn_index321(j)) * 
              dyn_alpha(j) * exp((-1*dyn_beta(j)) * pop_size); // Fi*ALPHA*exp(-BETA*n)
            
            if (substoch == 0 || dyn_type(j) == 2) {
              sparse_prophecy(dyn_index321(j)) = changing_element;
            } else if (substoch == 1 && dyn_type(j) == 1) {
              if (changing_element > 1.0) {
                changing_element = 1.0;
              } else if (changing_element < 0.0) {
                changing_element = 0.0;
              }
              sparse_prophecy(dyn_index321(j)) = changing_element;
            } else if (substoch == 2 && dyn_type(j) == 1) {
              changing_colsum = sum(sparse_prophecy.col(dyn_index_col(j))) - sparse_prophecy(dyn_index321(j));
              
              if (changing_element > (1.0 - changing_colsum)) {
                changing_element = (1.0 - changing_colsum);
              }
              sparse_prophecy(dyn_index321(j)) = changing_element;
            }
          } else if (dyn_style(j) == 2) { // Beverton-Holt
            changing_element = sparse_prophecy(dyn_index321(j)) * 
              dyn_alpha(j) / (1 + dyn_beta(j) * pop_size); // Fi*ALPHA/(1+BETA*n)
            
            if (substoch == 0 || dyn_type(j) == 2) {
              sparse_prophecy(dyn_index321(j)) = changing_element;
            } else if (substoch == 1 && dyn_type(j) == 1) {
              if (changing_element > 1.0) {
                changing_element = 1.0;
              } else if (changing_element < 0.0) {
                changing_element = 0.0;
              }
              sparse_prophecy(dyn_index321(j)) = changing_element;
            } else if (substoch == 2 && dyn_type(j) == 1) {
              changing_colsum = sum(sparse_prophecy.col(dyn_index_col(j))) - sparse_prophecy(dyn_index321(j));
              
              if (changing_element > (1.0 - changing_colsum)) {
                changing_element = (1.0 - changing_colsum);
              }
              sparse_prophecy(dyn_index321(j)) = changing_element;
            }
          } else if (dyn_style(j) == 3) { // Usher function
            changing_element = sparse_prophecy(dyn_index321(j)) * 
              (1 / (1 + exp(dyn_alpha(j) * pop_size + dyn_beta(j)))); // Fi*(1 / (1 + exp(alpha*N+b)))
            
            if (substoch == 0 || dyn_type(j) == 2) {
              sparse_prophecy(dyn_index321(j)) = changing_element;
            } else if (substoch == 1 && dyn_type(j) == 1) {
              if (changing_element > 1.0) {
                changing_element = 1.0;
              } else if (changing_element < 0.0) {
                changing_element = 0.0;
              }
              sparse_prophecy(dyn_index321(j)) = changing_element;
            } else if (substoch == 2 && dyn_type(j) == 1) {
              changing_colsum = sum(sparse_prophecy.col(dyn_index_col(j))) - sparse_prophecy(dyn_index321(j));
              
              if (changing_element > (1.0 - changing_colsum)) {
                changing_element = (1.0 - changing_colsum);
              }
              sparse_prophecy(dyn_index321(j)) = changing_element;
            }
          } else if (dyn_style(j) == 4) { // Logistic function
            changing_element = sparse_prophecy(dyn_index321(j)) * 
              (1 - pop_size / dyn_alpha(j)); // Fi*(1 - ALPHA/n)
            
            if (substoch == 0 || dyn_type(j) == 2) {
              sparse_prophecy(dyn_index321(j)) = changing_element;
            } else if (substoch == 1 && dyn_type(j) == 1) {
              if (changing_element > 1.0) {
                changing_element = 1.0;
              } else if (changing_element < 0.0) {
                changing_element = 0.0;
              }
              sparse_prophecy(dyn_index321(j)) = changing_element;
            } else if (substoch == 2 && dyn_type(j) == 1) {
              changing_colsum = sum(sparse_prophecy.col(dyn_index_col(j))) - sparse_prophecy(dyn_index321(j));
              
              if (changing_element > (1.0 - changing_colsum)) {
                changing_element = (1.0 - changing_colsum);
              }
              sparse_prophecy(dyn_index321(j)) = changing_element;
            }
          }
        }
      }
      
      sparse_seventhson = sparse_prophecy * sparse_seventhson;
      if (integeronly) {
        sparse_seventhson = floor(sparse_seventhson);
      }
      popproj.col(i+1) = arma::vec(arma::mat(sparse_seventhson));
      
      Rvecmat(i+1) = sum(popproj.col(i+1));
      
      if (!growthonly) {
        wpopproj.col(i+1) = arma::vec(arma::mat(sparse_seventhson)) / Rvecmat(i+1);
        
        sparse_secondprophecy = as<arma::sp_mat>(sparse_list[(mat_order(theclairvoyant - (i+1)))]);
        theseventhgrandson = theseventhgrandson * sparse_secondprophecy;
        
        double seventhgrandsum = sum(theseventhgrandson);
        arma::vec midwife = theseventhgrandson.as_col() / seventhgrandsum;
        
        theseventhgrandson = theseventhgrandson / seventhgrandsum;
        
        vpopproj.col(theclairvoyant - (i+1)) = midwife;
      }
    }
  }
  
  if (growthonly) {
    return popproj;
  } else {
    arma::mat revised_vproj = join_cols(vpopproj, Rvecmat);
    arma::mat expanded_proj = join_cols(wpopproj, revised_vproj);
    
    return join_cols(popproj, expanded_proj);
  }
}

//' Conduct Population Projection Simulations
//' 
//' Function \code{projection3()} runs projection simulations. It projects the
//' population an patches forward in time by a user-defined number of occasions.
//' Projections may be deterministic or stochastic, and may be density
//' dependent either way. If deterministic, then projections will be cyclical if
//' matrices exist covering multiple occasions for each population or patch. If
//' stochastic, then annual matrices will be shuffled within patches and
//' populations. Replicates may also be requested.
//' 
//' @param mpm A matrix projection model of class \code{lefkoMat}, or a list of
//' full matrix projection matrices.
//' @param nreps The number of replicate projections.
//' @param times Number of occasions to iterate per replicate. Defaults to
//' 10,000.
//' @param stochastic A logical value denoting whether to conduct a stochastic
//' projection or a deterministic / cyclical projection.
//' @param standardize A logical value denoting whether to re-standardize the
//' population size to 1.0 at each occasion. Defaults to FALSE.
//' @param growthonly A logical value indicating whether to produce only the
//' projected population size at each occasion, or a vector showing the stage
//' distribution followed by the reproductive value vector followed by the full
//' population size at each occasion. Defaults to TRUE.
//' @param integeronly A logical value indicating whether to round the number of
//' individuals projected in each stage at each occasion to the nearest
//' integer. Defaults to FALSE.
//' @param substoch An integer value indicating whether to force survival-
//' transition matrices to be substochastic in density dependent simulations.
//' Defaults to \code{0}, which does not force substochasticity. Alternatively,
//' \code{1} forces all survival-transition elements to range from 0.0 to 1.0,
//' and \code{2} forces all column rows to total no more than 1.0.
//' @param start_vec An optional numeric vector denoting the starting stage
//' distribution for the projection. Defaults to a single individual of each
//' stage.
//' @param start_frame An optional data frame characterizing stages, age-stages,
//' or stage-pairs that should be set to non-zero values in the starting vector,
//' and what those values should be. Can only be used with \code{lefkoMat}
//' objects.
//' @param tweights An optional numeric vector denoting the probabilistic
//' weightings of annual matrices. Defaults to equal weighting among occasions.
//' @param density An optional data frame describing the matrix elements that
//' will be subject to density dependence, and the exact kind of density
//' dependence that they will be subject to. The data frame used should be an
//' object of class \code{lefkoDens}, which is the output from function
//' \code{\link{density_input}()}.
//' 
//' @return A list of class \code{lefkoProj}, which always includes the first
//' three elements of the following, and also includes the remaining elements
//' below when a \code{lefkoMat} object is used as input:
//' \item{projection}{A list of lists of matrices showing the total number of
//' individuals per stage per occasion. The first list corresponds to each
//' pop-patch followed by each population. The inner list corresponds to
//' replicates within each pop-patch or population.}
//' \item{stage_dist}{A list of lists of the actual stage distribution in each
//' occasion in each replicate in each pop-patch or population. The list order
//' is the same as in \code{projection}.}
//' \item{rep_value}{A list of lists of the actual reproductive value in each
//' occasion in each replicate in each pop-patch or population. The list order
//' is the same as in \code{projection}.}
//' \item{pop_size}{A list of data frames showing the total population size in
//' each occasion per replicate (row within data frame) per pop-patch or
//' population (list element).}
//' \item{labels}{A data frame showing the order of populations and patches in
//' item \code{projection}.}
//' \item{control}{A short vector indicating the number of replicates and the
//' number of occasions projected per replicate.}
//' \item{ahstages}{The original stageframe used in the study.}
//' \item{hstages}{A data frame showing the order of historical stage pairs.}
//' \item{agestages}{A data frame showing the order of age-stage pairs.}
//' 
//' @section Notes:
//' Projections are run both at the patch level and at the population level.
//' Population level estimates will be noted at the end of the
//' data frame with 0 entries for patch designation.
//' 
//' Weightings given in \code{tweights} do not need to sum to 1. Final
//' weightings used will be based on the proportion per element of the sum of
//' elements in the user-supplied vector.
//' 
//' Starting vectors can be input in one of two ways: 1) as \code{start_vec}
//' input, which is a vector of numbers of the numbers of individuals in each
//' stage, stage pair, or age-stage, with the length of the vector necessarily
//' as long as there are rows in the matrices of the MPM; or 2) as
//' \code{start_frame} input, which is a data frame showing only those stages,
//' stage pairs, or age-stages that should begin with more than 0 individuals,
//' and the numbers of individuals that those stages should start with (this
//' object is created using the \code{\link{start_input}()} function). If both
//' are provided, then \code{start_frame} takes precedence and \code{start_vec}
//' is ignored. If neither is provided, then \code{projection3()} automatically
//' assumes that each stage, stage pair, or age-stage begins with a single
//' individual. Importantly, if a \code{lefkoMat} object is not used, and a list
//' of matrices is provided instead, then \code{start_frame} cannot be utilized
//' and a full \code{start_vec} must be provided to conduct a simulation with
//' starting numbers of individuals other than 1 per stage.
//' 
//' The resulting data frames in element \code{projection} are separated by
//' pop-patch according to the order provided in element \code{labels}, but the
//' matrices for each element of \code{projection} have the result of each
//' replicate stacked in order on top of one another without any break or
//' indication. Results for each replicate must be separated using the
//' information provided in elements \code{control} and the 3 stage
//' descriptor elements.
//' 
//' Density dependent projections are automatically set up if object
//' \code{density} is input. If this object is not included, then density
//' independent projections will be set up. Note that currently, density
//' dependent projections can only be performed with \code{lefkoMat} objects.
//' 
//' The stage distributions and reproductive values produced are not the
//' asymptotic values as would be given by the standardized right and left
//' eigenvectors associated with the dominant eigenvalue of a matrix, but are
//' vectors describing these values at the specific points in time projected.
//' See equations 14.86 and 14.88 and section 14.4 on Sensitivity and Elasticity
//' Analysis under Environmental Stochasticity in Caswell (2001, Matrix
//' Population Models, Sinauer Associates) for more details.
//' 
//' @seealso \code{\link{start_input}()}
//' @seealso \code{\link{density_input}()}
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
//' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl"), 
//'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep"),
//'   stage1 = c("Sd", "rep", "Sd", "rep", "all", "all"), 
//'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054),
//'   type = c(1, 1, 1, 1, 3, 3), type_t12 = c(1, 2, 1, 2, 1, 1),
//'   stageframe = lathframe, historical = TRUE)
//' 
//' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe,
//'   year = c(1989, 1990), stages = c("stage3", "stage2", "stage1"),
//'   supplement = lathsupp3, yearcol = "year2", indivcol = "individ")
//' 
//' lathproj <- projection3(ehrlen3, nreps = 5, stochastic = TRUE)
//' 
//' # Cypripedium example
//' rm(list = ls(all=TRUE))
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
//' cypsupp3r <- supplemental(stage3 = c("SD", "SD", "P1", "P1", "P2", "P3", "SL",
//'     "D", "XSm", "Sm", "D", "XSm", "Sm", "SD", "P1"),
//'   stage2 = c("SD", "SD", "SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL",
//'     "SL", "SL", "rep", "rep"),
//'   stage1 = c("SD", "rep", "SD", "rep", "SD", "P1", "P2", "P3", "P3", "P3",
//'     "SL", "SL", "SL", "mat", "mat"),
//'   eststage3 = c(NA, NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", "D", "XSm", "Sm",
//'     NA, NA),
//'   eststage2 = c(NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm", "XSm",
//'     "XSm", NA, NA),
//'   eststage1 = c(NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm", "XSm",
//'     "XSm", NA, NA),
//'   givenrate = c(0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.25, NA, NA, NA, NA, NA, NA,
//'     NA, NA),
//'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
//'   type = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
//'   type_t12 = c(1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
//'   stageframe = cypframe_raw, historical = TRUE)
//' 
//' cypmatrix3r <- rlefko3(data = cypraw_v1, stageframe = cypframe_raw, 
//'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
//'   size = c("size3added", "size2added", "size1added"), 
//'   supplement = cypsupp3r, yearcol = "year2", 
//'   patchcol = "patchid", indivcol = "individ")
//' 
//' cypstoch <- projection3(cypmatrix3r, nreps = 5, stochastic = TRUE)
//' 
//' @export projection3
// [[Rcpp::export]]
Rcpp::List projection3(List mpm, int nreps = 1, int times = 10000,
  bool stochastic = false, bool standardize = false, bool growthonly = true,
  bool integeronly = false, int substoch = 0,
  Nullable<NumericVector> start_vec = R_NilValue,
  Nullable<DataFrame> start_frame = R_NilValue,
  Nullable<NumericVector> tweights = R_NilValue,
  Nullable<DataFrame> density = R_NilValue) {
  
  Rcpp::List dens_index;
  Rcpp::DataFrame start_thru;
  Rcpp::DataFrame dens_input;
  
  int theclairvoyant {0};
  theclairvoyant = times;
  
  int dens_switch {0};
  int used_matsize {0};
  int total_projrows {0};
  
  if (theclairvoyant < 1) {
    throw Rcpp::exception("Option times must equal a positive integer.", false);
  }
  if (nreps < 1) {
    throw Rcpp::exception("Option nreps must be a positive integer.", false);
  }
  
  if (substoch != 0 && substoch != 1 && substoch != 2) {
    throw Rcpp::exception("Option substoch must be set to 0, 1, or 2.", false);
  }
  
  arma::uvec theprophecy(theclairvoyant);
  theprophecy.zeros();
  
  arma::vec startvec;
  arma::mat projection;
  
  if (!Rf_isMatrix(mpm[0])) {
    List amats = mpm["A"];
    List umats = mpm["U"];
    List fmats = mpm["F"];
    DataFrame stageframe = as<DataFrame>(mpm["ahstages"]);
    DataFrame hstages = as<DataFrame>(mpm["hstages"]);
    DataFrame labels = as<DataFrame>(mpm["labels"]);
    DataFrame agestages = as<DataFrame>(mpm["agestages"]);
    
    bool historical = false;
    bool agebystage = false;
    
    if (hstages.length() > 1) {
      historical = true;
    }
    if (agestages.length() > 1) {
      agebystage = true;
    }
    
    if (density.isNotNull()) { 
      Rcpp::DataFrame dens_thru(density);
      dens_input = dens_thru;
      dens_switch = 1;
      
      Rcpp::StringVector di_stage3 = dens_input["stage3"];
      Rcpp::StringVector di_stage2 = dens_input["stage2"];
      Rcpp::StringVector di_stage1 = dens_input["stage1"];
      int di_size = di_stage3.length();
      
      if (historical) {
        StringVector stage3 = hstages["stage_2"];
        StringVector stage2r = hstages["stage_1"];
        StringVector stage2c = hstages["stage_2"];
        StringVector stage1 = hstages["stage_1"];
        int hst_size = stage3.length();
        
        arma::uvec hst_3(hst_size);
        arma::uvec hst_2r(hst_size);
        arma::uvec hst_2c(hst_size);
        arma::uvec hst_1(hst_size);
        hst_3.zeros();
        hst_2r.zeros();
        hst_2c.zeros();
        hst_1.zeros();
        
        arma::uvec di_stage32_id(di_size);
        arma::uvec di_stage21_id(di_size);
        arma::uvec di_index(di_size);
        di_stage32_id.zeros();
        di_stage21_id.zeros();
        di_index.zeros();
        
        for (int i = 0; i < di_size; i++) { // This loop runs through each density_input line
          for (int j = 0; j < hst_size; j++) {
            if (di_stage3(i) == stage3(j)) {
              hst_3(j) = 1;
            } else {
              hst_3(j) = 0;
            }
          }
          
          for (int j = 0; j < hst_size; j++) {
            if (di_stage2(i) == stage2r(j)) {
              hst_2r(j) = 1;
            } else {
              hst_2r(j) = 0;
            }
          }
          
          for (int j = 0; j < hst_size; j++) {
            if (di_stage2(i) == stage2c(j)) {
              hst_2c(j) = 1;
            } else {
              hst_2c(j) = 0;
            }
          }
          
          for (int j = 0; j < hst_size; j++) {
            if (di_stage1(i) == stage1(j)) {
              hst_1(j) = 1;
            } else {
              hst_1(j) = 0;
            }
          }
          
          arma::uvec find_hst3 = find(hst_3);
          arma::uvec find_hst2r = find(hst_2r);
          arma::uvec find_hst2c = find(hst_2c);
          arma::uvec find_hst1 = find(hst_1);
          
          arma::uvec pop_32 = intersect(find_hst3, find_hst2r);
          arma::uvec pop_21 = intersect(find_hst2c, find_hst1);
          
          di_stage32_id(i) = pop_32(0);
          di_stage21_id(i) = pop_21(0);
          di_index(i) = pop_32(0) + (pop_21(0) * hst_size);
          
          // Eventually we need to zero the indexers out again
          hst_3.zeros();
          hst_2r.zeros();
          hst_2c.zeros();
          hst_1.zeros();
        }
        
        dens_index = Rcpp::List::create(_["index32"] = di_stage32_id,
          _["index21"] = di_stage21_id, _["index321"] = di_index);
      } else {
        StringVector stage3 = stageframe["stage"];
        StringVector stage2 = stageframe["stage"];
        int ahst_size = stage3.length();
        
        arma::uvec ahst_3(ahst_size);
        arma::uvec ahst_2(ahst_size);
        ahst_3.zeros();
        ahst_2.zeros();

        arma::uvec di_stage32_id(di_size);
        arma::uvec di_stage21_id(di_size);
        arma::uvec di_index(di_size);
        di_stage32_id.zeros();
        di_stage21_id.zeros();
        di_index.zeros();
        
        for (int i = 0; i < di_size; i++) { // This loop runs through each density_input line
          for (int j = 0; j < ahst_size; j++) {
            if (di_stage3(i) == stage3(j)) {
              ahst_3(j) = 1;
            } else {
              ahst_3(j) = 0;
            }
          }
          
          for (int j = 0; j < ahst_size; j++) {
            if (di_stage2(i) == stage2(j)) {
              ahst_2(j) = 1;
            } else {
              ahst_2(j) = 0;
            }
          }
          
          arma::uvec find_ahst3 = find(ahst_3);
          arma::uvec find_ahst2 = find(ahst_2);
          
          di_stage32_id(i) = find_ahst3(0);
          di_stage21_id(i) = find_ahst2(0);

          di_index(i) = find_ahst3(0) + (find_ahst2(0) * ahst_size);
          
          // Eventually we need to zero the indexers out again
          ahst_3.zeros();
          ahst_2.zeros();
        }
        
        dens_index = Rcpp::List::create(_["index3"] = di_stage32_id,
          _["index2"] = di_stage21_id, _["index321"] = di_index);
      }
    }
    
    if (labels.length() < 3) {
      throw Rcpp::exception("Function 'projection3' requires annual matrices. This lefkoMat object appears to be a set of mean matrices, and lacks annual matrices.", false);
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
        throw Rcpp::exception("Time weight vector must be the same length as the number of occasions represented in the lefkoMat object used as input.", false);
      }
      twinput = as<arma::vec>(tweights);
    } else {
      twinput.resize(yl);
      twinput.ones();
    } // At the end of this, we have an arma::vec holding with whatever the user supplied, or a 1 for each occasion
    
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
      mean_lefkomat = geodiesel(listofyears, umats, fmats, agestages,
        stageframe, 1, 1);
    } else {
      mean_lefkomat = turbogeodiesel(listofyears, umats, fmats, hstages,
        agestages, stageframe, 1, 1);
    }
    
    // Here we take the matrices corresponding to each individual patch, run the
    // simulation, and estimate all descriptive metrics
    List meanamats = mean_lefkomat["A"];
    List mmlabels = mean_lefkomat["labels"];
    StringVector mmpops = mmlabels["pop"];
    StringVector mmpatches = mmlabels["patch"];
    
    arma::mat thechosenone = as<arma::mat>(meanamats[0]);
    int meanmatsize = thechosenone.n_elem;
    int meanmatrows = thechosenone.n_rows;
    arma::vec startvec;
    int trials = meanamats.length();
    used_matsize = meanmatrows;
    
    arma::uvec ppcindex = as<arma::uvec>(poppatchc);
    arma::uvec allppcs = as<arma::uvec>(sort_unique(poppatchc));
    int allppcsnem = allppcs.n_elem;
    List plist_hold(allppcsnem);
    int pop_est {1};
    if (allppcsnem > 1) {
      pop_est = trials - allppcsnem;    
    }
    List projection_list(trials);
    
    if(start_frame.isNotNull()) {
      Rcpp::DataFrame start_thru(start_frame);
      
      startvec.set_size(meanmatrows);
      startvec.zeros();
      
      arma::uvec start_elems = start_thru["row_num"];
      start_elems = start_elems - 1;
      arma::vec start_values = start_thru["value"];
      
      if (start_elems.max() > (meanmatrows - 1)) {
        throw Rcpp::exception("Start vector input frame includes element indices too high for this MPM", false);
      }
      
      for (int i = 0; i < start_elems.n_elem; i++) {
        startvec(start_elems(i)) = start_values(i);
      }
    } else if (start_vec.isNotNull()) {
      if (as<NumericVector>(start_vec).length() != meanmatrows) {
        throw Rcpp::exception("Start vector must be the same length as the number of rows in each matrix.", false);
      }
      startvec = as<arma::vec>(start_vec);
      
    } else {
      startvec.set_size(meanmatrows);
      startvec.ones();
    } // At the end of this, we have an arma::vec holding with whatever the user supplied, or a 1 for each stage
    
    twinput = twinput / sum(twinput);
    
    for (int i= 0; i < allppcsnem; i++) {
      thechosenone = as<arma::mat>(meanamats[i]);
      
      arma::uvec thenumbersofthebeast = find(ppcindex == allppcs(i));
      int chosen_yl = thenumbersofthebeast.n_elem;
      
      // This loop takes care of multiple replicates, creating the final data frame
      // of results for each pop-patch
      for (int rep = 0; rep < nreps; rep++) {
        if (stochastic) {
          theprophecy = Rcpp::RcppArmadillo::sample(thenumbersofthebeast, theclairvoyant, true, twinput);
        } else {
          theprophecy.set_size(theclairvoyant);
          theprophecy.zeros();
          
          for (int j = 0; j < theclairvoyant; j++) {
            theprophecy(j) = thenumbersofthebeast(j % chosen_yl);
          }
        }
        
        if (dens_switch) {
          if (rep == 0) {
            projection = proj3dens(startvec, amats, theprophecy, growthonly,
              integeronly, substoch, dens_input, dens_index);
          } else {
            arma::mat nextproj = proj3dens(startvec, amats, theprophecy,
              growthonly, integeronly, substoch, dens_input, dens_index);
            projection = arma::join_cols(projection, nextproj);
          }
        } else {
          if (rep == 0) {
            projection = proj3(startvec, amats, theprophecy, standardize, growthonly,
              integeronly);
          } else {
            arma::mat nextproj = proj3(startvec, amats, theprophecy, standardize,
              growthonly, integeronly);
            projection = arma::join_cols(projection, nextproj);
          }
        }
      }
      
      projection_list(i) = projection;
    }
    
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
        int chosen_yl = choicevec.n_elem;
      
        // This loop takes care of multiple replicates, creating the final data frame
        // of results for the pop mean(s)
        for (int rep = 0; rep < nreps; rep++) {
          if (stochastic) {
            theprophecy = Rcpp::RcppArmadillo::sample(choicevec, theclairvoyant, true, twinput);
          } else {
            theprophecy.zeros();
            
            for (int j = 0; j < theclairvoyant; j++) {
              theprophecy(j) = choicevec(j % chosen_yl);
            }
          }
          
          if (dens_switch) {
            if (rep == 0) {
              projection = proj3dens(startvec, meanmatyearlist, theprophecy,
                growthonly, integeronly, substoch, dens_input, dens_index);
            } else {
              arma::mat nextproj = proj3dens(startvec, meanmatyearlist, theprophecy,
                growthonly, integeronly, substoch, dens_input, dens_index);
              projection = arma::join_cols(projection, nextproj);
            }
          } else {
            if (rep == 0) {
              projection = proj3(startvec, meanmatyearlist, theprophecy,
                standardize, growthonly, integeronly);
            } else {
              arma::mat nextproj = proj3(startvec, meanmatyearlist, theprophecy,
                standardize, growthonly, integeronly);
              projection = arma::join_cols(projection, nextproj);
            }
          }
        }
        projection_list(allppcsnem + i) = projection;
      }
    }
    
    // The final output will have a projection list with # elements = nreps, nested
    // within a list with # elements = # pop-patches
    List projection_set(nreps);
    List ss_set(nreps);
    List rv_set(nreps);
    arma::mat total_sizes_set(nreps, (times+1), fill::zeros);
    
    int length_ppy = projection_list.length();
    
    List final_projection(length_ppy);
    List final_ss(length_ppy);
    List final_rv(length_ppy);
    List final_ns(length_ppy);
    
    List output(9);
    
    if (!growthonly) {
      
      arma::mat list_proj(total_projrows, (times+1), fill::zeros);
      arma::mat extracted_proj(used_matsize, used_matsize, fill::zeros);
      
      int diversion = used_matsize * 3 + 1;
      
      for (int j = 0; j < length_ppy; j++) {
        list_proj = as<arma::mat>(projection_list[j]);
        
        for (int i = 0; i < nreps; i++) {
          extracted_proj = list_proj.rows((diversion * i), (diversion * i + used_matsize - 1));
          projection_set(i) = extracted_proj;
          
          extracted_proj = list_proj.rows((diversion * i + used_matsize), (diversion * i + (2 * used_matsize) - 1));
          ss_set(i) = extracted_proj;
          
          extracted_proj = list_proj.rows((diversion * i + (2 * used_matsize)), (diversion * i + (3 * used_matsize) - 1));
          rv_set(i) = extracted_proj;
          
          total_sizes_set.row(i) = list_proj.row(diversion * (i+1) - 1);
        }
        final_projection(j) = clone(projection_set);
        final_ss(j) = clone(ss_set);
        final_rv(j) = clone(rv_set);
        final_ns(j) = total_sizes_set;
      }
      
    } else {
      
      arma::mat list_proj(total_projrows, (times+1), fill::zeros);
      arma::mat extracted_proj(used_matsize, used_matsize, fill::zeros);
      
      int diversion = used_matsize;
      
      for (int j = 0; j < length_ppy; j++) {
        list_proj = as<arma::mat>(projection_list[j]);
        
        for (int i = 0; i < nreps; i++) {
          extracted_proj = list_proj.rows((diversion * i), (diversion * i + used_matsize - 1));
          projection_set(i) = extracted_proj;
          
          total_sizes_set.row(i) = sum(extracted_proj, 0);
        }
        
        final_projection(j) = clone(projection_set);
        final_ss(j) = NULL;
        final_rv(j) = NULL;
        final_ns(j) = total_sizes_set;
      }
      
    }
    
    DataFrame newlabels = DataFrame::create(_["pop"] = mmpops,
      _["patch"] = mmpatches);
    Rcpp::IntegerVector control = {nreps, times};
    
    output(0) = final_projection;
    output(1) = final_ss;
    output(2) = final_rv;
    output(3) = final_ns;
    output(4) = newlabels;
    output(5) = stageframe;
    output(6) = hstages;
    output(7) = agestages;
    output(8) = control;
    
    CharacterVector namevec = {"projection", "stage_dist", "rep_value", "pop_size",
      "labels", "ahstages", "hstages", "agestages", "control"};
    output.attr("names") = namevec;
    
    if (dens_switch) {
      output.push_back(dens_index);
    }
    output.attr("class") = "lefkoProj";
    
    return output;
    
  } else {
    
    List projection_list (1);
    List amats = mpm;
    
    int yl = amats.length();
    arma::mat firstmat = as<arma::mat>(amats[0]);
    int matrows = firstmat.n_rows;
    int matcols = firstmat.n_cols;
    used_matsize = matrows;
    
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
      throw Rcpp::exception("Supplied matrices must be square. Please check matrix dimensions and fix.", false);
    }
    
    if (tweights.isNotNull()) {
      if (as<NumericVector>(tweights).length() != yl) {
        throw Rcpp::exception("Time weight vector must be the same length as the number of occasions represented in the lefkoMat object used as input.", false);
      }
      twinput = as<arma::vec>(tweights);
    } else {
      twinput.resize(yl);
      twinput.ones();
    } // At the end of this, we have an arma::vec holding with whatever the user supplied, or a 1 for each occasion
    
    if (start_vec.isNotNull()) {
      if (as<NumericVector>(start_vec).length() != matrows) {
        throw Rcpp::exception("Start vector must be the same length as the number of rows in each matrix.", false);
      }
      startvec = as<arma::vec>(start_vec);
    } else {
      startvec.set_size(matrows);
      startvec.ones();
    } // At the end of this, we have an arma::vec holding with whatever the user supplied, or a 1 for each stage
    
    // Now we create the mean matrix
    arma::mat thechosenone(matrows, matcols);
    thechosenone.zeros();
    
    for (int i = 0; i < yl; i++) {
      arma::mat columnified = as<arma::mat>(amats[i]);
      thechosenone = thechosenone + (columnified / yl);
    }
    
    // Here we take the matrices corresponding to each individual patch, run the simulation, and
    // estimate all descriptive metrics
    twinput = twinput / sum(twinput);
    
    arma::uvec thenumbersofthebeast = uniqueyears;
    
    // This loop takes care of multiple replicates, creating the final data frame
    // of results
    for (int rep = 0; rep < nreps; rep++) {
      if (stochastic) {
        theprophecy = Rcpp::RcppArmadillo::sample(thenumbersofthebeast, theclairvoyant, true, twinput);
      } else {
        theprophecy.zeros();
        
        for (int i = 0; i < theclairvoyant; i++) {
          theprophecy(i) = thenumbersofthebeast(i % yl);
        }
      }
      if (rep == 0) {
        projection = proj3(startvec, amats, theprophecy, standardize, growthonly, integeronly);
      } else {
        arma::mat nextproj = proj3(startvec, amats, theprophecy, standardize, growthonly, integeronly);
        projection = arma::join_cols(projection, nextproj);
      }
    }
    
    projection_list(0) = projection;
    
    DataFrame newlabels = DataFrame::create(_["pop"] = 1,
      _["patch"] = 1);
    
    Rcpp::IntegerVector control = {nreps, times};
    
    Rcpp::List output = List::create(_["projection"] = projection_list, _["labels"] = newlabels,
      _["control"] = control);
    output.attr("class") = "lefkoProj";
    
    return output;
  }
}

//' Estimate Stochastic Population Growth Rate
//' 
//' Function \code{slambda3()} estimates the stochastic population growth rate,
//' \eqn{a}, defined as the long-term arithmetic mean of the log population 
//' growth rate estimated per simulated occasion (as given in equation 2 in
//' Tuljapurkar, Horvitz, and Pascarella 2003). This term is estimated via
//' projection of randomly sampled matrices, similarly to the procedure outlined
//' in Box 7.4 of Morris and Doak (2002). Can handle both lefkoMat objects and
//' lists of full A matrices. 
//' 
//' @param mpm A matrix projection model of class \code{lefkoMat}, or a list of
//' full matrix projection matrices.
//' @param times Number of occasions to iterate. Defaults to 10,000.
//' @param dense_only A logical value indicating whether to force matrices to be
//' run in dense format. Defaults to \code{FALSE}, and should only be used if
//' errors occur when running under default conditions.
//' @param tweights Numeric vector denoting the probabilistic weightings of
//' annual matrices. Defaults to equal weighting among occasions.
//' 
//' @return A data frame with the following variables:
//' 
//' \item{pop}{The identity of the population.}
//' \item{patch}{The identity of the patch.}
//' \item{a}{Estimate of stochastic growth rate, estimated as the arithmetic
//' mean of the log population growth rate across simulated occasions.}
//' \item{var}{The estimated variance of a.}
//' \item{sd}{The standard deviation of a.}
//' \item{se}{The standard error of a.}
//'
//' @section Notes:
//' Stochastic growth rate is estimated both at the patch level and at the
//' population level. Population level estimates will be noted at the end of the
//' data frame with 0 entries for patch designation.
//' 
//' Weightings given in \code{tweights} do not need to sum to 1. Final
//' weightings used will be based on the proportion per element of the sum of
//' elements in the user-supplied vector.
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
//' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl"), 
//'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep"),
//'   stage1 = c("Sd", "rep", "Sd", "rep", "all", "all"), 
//'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054),
//'   type = c(1, 1, 1, 1, 3, 3), type_t12 = c(1, 2, 1, 2, 1, 1),
//'   stageframe = lathframe, historical = TRUE)
//' 
//' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe,
//'   year = c(1989, 1990), stages = c("stage3", "stage2", "stage1"),
//'   supplement = lathsupp3, yearcol = "year2", indivcol = "individ")
//' 
//' slambda3(ehrlen3)
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
//' cypsupp3r <- supplemental(stage3 = c("SD", "SD", "P1", "P1", "P2", "P3", "SL",
//'     "D", "XSm", "Sm", "D", "XSm", "Sm", "SD", "P1"),
//'   stage2 = c("SD", "SD", "SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL",
//'     "SL", "SL", "rep", "rep"),
//'   stage1 = c("SD", "rep", "SD", "rep", "SD", "P1", "P2", "P3", "P3", "P3",
//'     "SL", "SL", "SL", "mat", "mat"),
//'   eststage3 = c(NA, NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", "D", "XSm", "Sm",
//'     NA, NA),
//'   eststage2 = c(NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm", "XSm",
//'     "XSm", NA, NA),
//'   eststage1 = c(NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm", "XSm",
//'     "XSm", NA, NA),
//'   givenrate = c(0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.25, NA, NA, NA, NA, NA, NA,
//'     NA, NA),
//'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
//'   type = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
//'   type_t12 = c(1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
//'   stageframe = cypframe_raw, historical = TRUE)
//' 
//' cypmatrix3r <- rlefko3(data = cypraw_v1, stageframe = cypframe_raw, 
//'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
//'   size = c("size3added", "size2added", "size1added"), 
//'   supplement = cypsupp3r, yearcol = "year2", 
//'   patchcol = "patchid", indivcol = "individ")
//' 
//' cypstoch <- slambda3(cypmatrix3r, dense_only = TRUE)
//' cypstoch
//' 
//' @export slambda3
// [[Rcpp::export]]
DataFrame slambda3(List mpm, int times = 10000, bool dense_only = false,
  Nullable<NumericVector> tweights = R_NilValue) {
  
  int theclairvoyant {0};
  int sparse_switch {0};
  
  theclairvoyant = times;
  
  if (theclairvoyant < 1) {
    throw Rcpp::exception("Option must equal a positive integer.", false);
  }
  
  if (!Rf_isMatrix(mpm[0])) {
    List amats = mpm["A"];
    List umats = mpm["U"];
    List fmats = mpm["F"];
    DataFrame stageframe = as<DataFrame>(mpm["ahstages"]);
    DataFrame hstages = as<DataFrame>(mpm["hstages"]);
    DataFrame labels = as<DataFrame>(mpm["labels"]);
    DataFrame agestages = as<DataFrame>(mpm["agestages"]);
    
    bool historical;
    
    if (hstages.length() > 1) {
      historical = true;
    } else {
      historical = false;
    }
    
    // Here we will check if the matrix is large and sparse
    int test_elems = as<arma::mat>(amats(0)).n_elem;
    arma::uvec nonzero_elems = find(as<arma::mat>(amats(0)));
    int all_nonzeros = nonzero_elems.n_elem;
    double sparse_check = static_cast<double>(all_nonzeros) / static_cast<double>(test_elems);
    if (sparse_check <= 0.5) {
      sparse_switch = 1;
    } else sparse_switch = 0;
    
    if (dense_only) sparse_switch = 0;
    
    if (labels.length() < 3) {
      throw Rcpp::exception("Function 'slambda3' requires annual matrices. This lefkoMat object appears to be a set of mean matrices, and lacks annual matrices.", false);
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
        throw Rcpp::exception("Time weight vector must be the same length as the number of occasions represented in the lefkoMat object used as input.", false);
      }
      twinput = as<arma::vec>(tweights);
    } else {
      twinput.resize(yl);
      twinput.ones();
    } // At the end of this, we have an arma::vec holding with whatever the user supplied, or a 1 for each occasion
    
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
      mean_lefkomat = geodiesel(listofyears, umats, fmats, agestages,
        stageframe, 1, 1);
    } else {
      mean_lefkomat = turbogeodiesel(listofyears, umats, fmats, hstages,
        agestages, stageframe, 1, 1);
    }
    
    // Here we take the matrices corresponding to each individual patch, run the simulation, and
    // estimate all descriptive metrics
    List meanamats = mean_lefkomat["A"];
    List mmlabels = mean_lefkomat["labels"];
    StringVector mmpops = mmlabels["pop"];
    StringVector mmpatches = mmlabels["patch"];
    
    int meanmatsize = as<arma::mat>(meanamats[0]).n_elem; // thechosenone
    int meanmatrows = as<arma::mat>(meanamats[0]).n_rows; // thechosenone
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
    
    twinput = twinput / sum(twinput);
    
    for (int i= 0; i < allppcsnem; i++) {
      arma::uvec thenumbersofthebeast = find(ppcindex == allppcs(i));
      arma::uvec theprophecy = Rcpp::RcppArmadillo::sample(thenumbersofthebeast, theclairvoyant, true, twinput);
      
      startvec = ss3matrix(as<arma::mat>(meanamats[i]), sparse_switch); // thechosenone
      
      arma::mat projection = proj3(startvec, amats, theprophecy, 1, 1, 0);
      
      for (int j = 0; j < theclairvoyant; j++) {
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
        startvec = ss3matrix(as<arma::mat>(meanamats[allppcsnem + i]), sparse_switch); // thechosenone
        
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
        arma::uvec theprophecy = Rcpp::RcppArmadillo::sample(choicevec, theclairvoyant, true, twinput);
        
        arma::mat projection = proj3(startvec, meanmatyearlist, theprophecy, 1, 1, 0);
        
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
    
    // Here we will check if the matrix is large and sparse
    int test_elems = as<arma::mat>(amats(0)).n_elem;
    arma::uvec nonzero_elems = find(as<arma::mat>(amats(0)));
    int all_nonzeros = nonzero_elems.n_elem;
    double sparse_check = static_cast<double>(all_nonzeros) / static_cast<double>(test_elems);
    if (sparse_check <= 0.5) {
      sparse_switch = 1;
    } else sparse_switch = 0;
    
    if (dense_only) sparse_switch = 0;
    
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
      throw Rcpp::exception("Supplied matrices must be square. Please check matrix dimensions and fix.", false);
    }
    
    if (tweights.isNotNull()) {
      if (as<NumericVector>(tweights).length() != yl) {
        throw Rcpp::exception("Time weight vector must be the same length as the number of occasions represented in the lefkoMat object used as input.", false);
      }
      twinput = as<arma::vec>(tweights);
    } else {
      twinput.resize(yl);
      twinput.ones();
    } // At the end of this, we have an arma::vec holding whatever the user supplied, or a 1 for each occasion
    
    // Now we create the mean matrix
    arma::mat thechosenone(matrows, matcols);
    thechosenone.zeros();
    
    for (int i = 0; i < yl; i++) {
      thechosenone = thechosenone + (as<arma::mat>(amats[i]) / yl);
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
    
    twinput = twinput / sum(twinput);
    
    arma::uvec thenumbersofthebeast = uniqueyears;
    arma::uvec theprophecy = Rcpp::RcppArmadillo::sample(thenumbersofthebeast, theclairvoyant, true, twinput);
    
    startvec = ss3matrix(thechosenone, sparse_switch);
    
    arma::mat projection = proj3(startvec, amats, theprophecy, 1, 1, 0);
    
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

//' Estimate Stochastic Sensitivity or Elasticity of Matrix Set
//' 
//' Function \code{stoch_senselas()} estimates the sensitivity and elasticity to
//' matrix elements of \eqn{a}, defined as the long-term arithmetic mean of the
//' log population growth estimated per simulated occasion (as given in equation 2
//' in Tuljapurkar, Horvitz, and Pascarella 2003). 
//' 
//' @param mpm A matrix projection model of class \code{lefkoMat}, or a list of
//' full matrix projection matrices.
//' @param times Number of occasions to iterate. Defaults to 10,000.
//' @param style An integer designating whether to estimate sensitivity matrices
//' (\code{1}) or elasticity matrices (\code{2}). Defaults to 1.
//' @param tweights Numeric vector denoting the probabilistic weightings of
//' annual matrices. Defaults to equal weighting among occasions.
//' 
//' @return A list of one or two cubes (3d array) where each slice corresponds
//' to a sensitivity or elasticity matrix for a specific pop-patch, followed by
//' the sensitivity or elasticity matrices of all populations (only if multiple
//' pop-patches occur in the input). Two such cubes are only provided when a
//' historical lefkoMat object is used as input, in which case the first
//' element is the historical sensitivity/elasticity matrix, and the second is
//' the ahistorical sensitivity/elasticity matrix.
//' 
//' @section Notes:
//' Weightings given in \code{tweights} do not need to sum to 1. Final
//' weightings used will be based on the proportion per element of the sum of
//' elements in the user-supplied vector.
//' 
//' This function currently requires all patches to have the same occasions, if
//' a \code{lefkoMat} object is used as input. Asymmetry in the number of
//' occasions across patches and/or populations will likely cause errors.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export(.stoch_senselas)]]
Rcpp::List stoch_senselas(List mpm, int times = 10000, int style = 1,
  Nullable<NumericVector> tweights = R_NilValue) {
  
  int theclairvoyant = times;
  
  if (theclairvoyant < 1) {
    throw Rcpp::exception("Option must equal a positive integer.", false);
  }
  
  List projlist = List::create(Named("delete") = 1);
  
  if (!Rf_isMatrix(mpm[0])) {
    List amats = mpm["A"];
    List umats = mpm["U"];
    List fmats = mpm["F"];
    DataFrame stageframe = as<DataFrame>(mpm["ahstages"]);
    DataFrame hstages = as<DataFrame>(mpm["hstages"]);
    DataFrame labels = as<DataFrame>(mpm["labels"]);
    DataFrame agestages = as<DataFrame>(mpm["agestages"]);
    
    // The next lines are necessary to assess ahistorical versions of historical
    // sensitivities and potentially elasticities
    arma::uvec ahstages_id = stageframe["stage_id"];
    StringVector ahstages_name = stageframe["stage"];
    int ahstages_num = ahstages_id.n_elem;
    
    arma::uvec hstages_id2(ahstages_num * ahstages_num);
    hstages_id2.zeros();
    int hstages_num {0};
    
    bool historical;
    
    if (hstages.length() > 1) {
      historical = true;
      
      arma::uvec hstages_id = hstages["stage_id_2"];
      
      hstages_num = hstages_id.n_elem;
      
      for (int i = 0; i < hstages_num; i++) {
        hstages_id2(i) = hstages_id(i);
      }
    } else {
      historical = false;
    }
    
    if (labels.length() < 3) {
      Rf_warningcall(R_NilValue, "This function requires annual matrices as input. This lefkoMat object appears to be a set of mean matrices, and may lack annual matrices.");
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
    arma::uvec uniqueyears_arma = as<arma::uvec>(uniqueyears);
    
    IntegerVector popc = match(poporder, uniquepops) - 1;
    IntegerVector poppatchc = match(poppatch, uniquepoppatches) - 1;
    IntegerVector year2c = match(yearorder, uniqueyears) - 1;
    int yl = uniqueyears.length();
    
    arma::vec twinput;
    
    if (tweights.isNotNull()) {
      if (as<NumericVector>(tweights).length() != yl) {
        throw Rcpp::exception("Time weight vector must be the same length as the number of occasions represented in the lefkoMat object used as input.", false);
      }
      twinput = as<arma::vec>(tweights);
    } else {
      twinput.resize(yl);
      twinput.ones();
    } // At the end of this, we have an arma::vec holding with whatever the user supplied, or a 1 for each occasion
    
    arma::vec twinput_corr = twinput / sum(twinput);
    arma::uvec theprophecy_allyears = Rcpp::RcppArmadillo::sample(uniqueyears_arma, theclairvoyant, true, twinput_corr);
    
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
    
    // Now we will create a set of means for patches and populations
    List mean_lefkomat;
    
    if (hstages.length() == 1) {
      mean_lefkomat = geodiesel(listofyears, umats, fmats, agestages,
        stageframe, 1, 1);
    } else {
      mean_lefkomat = turbogeodiesel(listofyears, umats, fmats, hstages,
        agestages, stageframe, 1, 1);
    }
    
    // Now we will set up the preliminaries for the stochastic simulations
    List meanamats = mean_lefkomat["A"];
    List mmlabels = mean_lefkomat["labels"];
    StringVector mmpops = mmlabels["pop"];
    StringVector mmpatches = mmlabels["patch"];
    
    int meanmatsize = as<arma::mat>(meanamats[0]).n_elem; // thechosenone
    int meanmatrows = as<arma::mat>(meanamats[0]).n_rows; // thechosenone
    arma::vec startvec(meanmatrows);
    startvec.ones();
    startvec = startvec / meanmatrows; // This is the start vector for w and v calculations
    int trials = meanamats.length();
    
    // Here we initialize two cubes to hold sensitivity or elasticity matrices, the
    // first for general use while the second is specifically for ahistorical versions
    // of historical matrices
    arma::cube senscube(meanmatrows, meanmatrows, trials);
    arma::cube senscube_ah(ahstages_num, ahstages_num, trials);
    
    senscube.zeros();
    senscube_ah.zeros();
    
    // This next matrix will hold the year values for each run
    arma::umat yearspulled(trials, theclairvoyant);
    yearspulled.zeros();
    
    arma::uvec ppcindex = as<arma::uvec>(poppatchc);
    arma::uvec allppcs = as<arma::uvec>(sort_unique(poppatchc));
    int allppcsnem = allppcs.n_elem;
    
    arma::uvec year2arma = as<arma::uvec>(yearorder);
    
    // These matrices and vectors will hold R values
    arma::mat Rvecmat(trials, theclairvoyant);
    Rvecmat.zeros();
    
    for (int i= 0; i < allppcsnem; i++) { // This loop goes through each pop-patch
      arma::uvec theprophecy = theprophecy_allyears;
      theprophecy.zeros();
      
      arma::uvec tnotb_patch = find(ppcindex == allppcs(i));
      
      for (int j = 0; j < yl; j++) { // This creates the main index that will mark the correct matrices to use; needs to be modified for situations in which patches do not have the same years
        arma::uvec tnotb_years = find(year2arma == uniqueyears(j));
        arma::uvec thenumbersofthebeast = intersect(tnotb_patch, tnotb_years);
        
        if (thenumbersofthebeast.n_elem > 0) {
          arma::uvec prophetic_yearindices = find(theprophecy_allyears == uniqueyears(j));
          
          if (prophetic_yearindices.n_elem > 0) {
            int replacement = thenumbersofthebeast(0);
            theprophecy.elem(prophetic_yearindices).fill(replacement);
          }
        }
      }
      yearspulled.row(i) = theprophecy.t();
      
      // The next section creates stable stage and rep value vectors arranged in
      // matrix format. The first two are general for whatever has been input,
      // whether historical or ahistorical, while the next two are specifically
      // for ahistorical versions of historical inputs
      arma::mat wprojection(startvec.n_elem, (theclairvoyant + 1));
      arma::mat vprojection(startvec.n_elem, (theclairvoyant + 1));
      arma::vec wprojection_ah(ahstages_num);
      arma::vec vprojection_ah(ahstages_num);
      wprojection.zeros();
      vprojection.zeros();
      wprojection_ah.zeros();
      vprojection_ah.zeros();
      
      // Here we run the control loop to create the w and v values we need
      arma::mat crazy_prophet = proj3(startvec, amats, theprophecy, 1, 0, 0);
      wprojection = crazy_prophet.submat(startvec.n_elem, 0, ((startvec.n_elem * 2) - 1), theclairvoyant);
      vprojection = crazy_prophet.submat((startvec.n_elem * 2), 0, ((startvec.n_elem * 3) - 1), theclairvoyant);
      Rvecmat.row(i) = crazy_prophet.submat((startvec.n_elem * 3), 1, (startvec.n_elem * 3), theclairvoyant); // Rvec
      
      // All references should go to senscube, which is a 3d array designed to hold the sensitivity matrices
      for (int j = 0; j < theclairvoyant; j++) { // This is the main time loop for the sensitivity matrices, 
                                                 // adding each occasion to the respective matrix for each pop-patch
        arma::vec vtplus1 = vprojection.col(j+1);
        arma::vec wtplus1 = wprojection.col(j+1);
        arma::vec wt = wprojection.col(j);
        
        arma::mat currentsens_num = vtplus1 * wt.as_row(); // This is the numerator of the key matrix equation
        arma::mat currentsens_den = (Rvecmat(i, j) * vtplus1.as_row() * wtplus1); // This is the denominator of the key matrix equation
        double cd_double = currentsens_den(0,0);
        arma::mat currentsens = currentsens_num / (cd_double * theclairvoyant);
        
        // This creates the sensitivity matrices
        if (style == 1) {
          senscube.slice(i) += currentsens;
          
          if (historical) {
            wprojection_ah.zeros();
            vprojection_ah.zeros();
            
            // This loop creates the ahistorical stable stage distribution for projected
            // occasion j+1
            for (int k1 = 0; k1 < hstages_num; k1++) {
              int current_stage2 = hstages_id2(k1);
              wprojection_ah(current_stage2 - 1) = wprojection_ah(current_stage2 - 1)  +
                wtplus1(k1);
            } // k1 loop
            
            // Now the ahistorical reproductive value vector for occasion j+1
            for (int k2 = 0; k2 < hstages_num; k2++) {
              int current_stage2 = hstages_id2(k2);
              
              if (wprojection_ah(current_stage2 - 1) > 0) {
                vprojection_ah(current_stage2 - 1) = vprojection_ah(current_stage2 - 1) +
                  (vtplus1(k2) * wtplus1(k2) / wprojection_ah(current_stage2 - 1));
              }
            } // k2 loop
            
            // Now to propagate the projection sensitivity matrix, and add it to
            // the main sensitivity matrix
            arma::rowvec wtah_tpose = wprojection_ah.as_row();
            arma::rowvec vtah_tpose = vprojection_ah.as_row();
            arma::mat csah_num = vprojection_ah * wtah_tpose;
            arma::mat csah_den = (Rvecmat(i, j) * vtah_tpose * wprojection_ah);
            double cdah_double = csah_den(0,0);
            arma::mat csah = csah_num / (cdah_double * theclairvoyant);
            senscube_ah.slice(i) += csah;
          } // if historical statement
        } else {
          // This creates the elasticity matrices
          senscube.slice(i) += currentsens % as<arma::mat>(amats[(theprophecy(j))]);
        }
      }
    }
    
    // This section works on the pop means
    int pop_est {1};
    StringVector allpops = unique(poporder);
    arma::uvec popmatch(loysize);
    arma::uvec yearmatch(loysize);
    popmatch.zeros();
    yearmatch.zeros();
    List meanmatyearlist(yl);
    
    IntegerVector tnotb_all = seq(0, (yl - 1));
    arma::uvec theprophecy = theprophecy_allyears;
    theprophecy.zeros();
    
    for (int j = 0; j < yl; j++) { // This creates the main index that will mark the correct matrices to use; needs to be modified for situations in which patches do not have the same years
      arma::uvec prophetic_yearindices = find(theprophecy_allyears == uniqueyears(j));
      
      if (prophetic_yearindices.n_elem > 0) {
        theprophecy.elem(prophetic_yearindices).fill(j);
      }
    }
    
    if (allppcsnem > 1) { // This checks for pop-mean matrices, which would only occur if there are more than 1 pop-patches
      pop_est = trials - allppcsnem;
      
      for (int i = 0; i < pop_est; i++) { // This loop goes through each population
        for (int j = 0; j < loysize; j++) { // This checks which A matrices match the current population in the loop
          if (poporder(j) == allpops(i)) {
            popmatch(j) = 1;
          } else {
            popmatch(j) = 0;
          }
        }
        
        arma::uvec neededmatspop = find(popmatch == 1);
        
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
        yearspulled.row(allppcsnem + i) = theprophecy.t();
        
        // Now we need to use meanmatyearlist in place of amats
        arma::mat wprojection(startvec.n_elem, (theclairvoyant + 1));
        arma::mat vprojection(startvec.n_elem, (theclairvoyant + 1));
        arma::vec wprojection_ah(ahstages_num);
        arma::vec vprojection_ah(ahstages_num);
        wprojection.zeros();
        vprojection.zeros();
        wprojection_ah.zeros();
        vprojection_ah.zeros();
        
        // Here we run the control loop to create the w and v values we need
        arma::mat crazy_prophet = proj3(startvec, meanmatyearlist, theprophecy, 1, 0, 0);
        wprojection = crazy_prophet.submat(startvec.n_elem, 0, ((startvec.n_elem * 2) - 1), theclairvoyant);
        vprojection = crazy_prophet.submat((startvec.n_elem * 2), 0, ((startvec.n_elem * 3) - 1), theclairvoyant);
        Rvecmat.row(allppcsnem + i) = crazy_prophet.submat((startvec.n_elem * 3), 1, (startvec.n_elem * 3), theclairvoyant); // Rvec
        
        // All references should go to senscube, which is a 3d array designed to
        // hold the sensitivity matrices
        
        // Next is the main time loop for the sensitivity matrices, adding each
        // occasion to the respective matrix for each pop-patch
        for (int j = 0; j < theclairvoyant; j++) {  
          
          arma::vec vtplus1 = vprojection.col(j+1);
          arma::vec wtplus1 = wprojection.col(j+1);
          arma::vec wt = wprojection.col(j);
          
          arma::mat currentsens_num = vtplus1 * wt.as_row(); // This is the numerator of the key matrix equation
          arma::mat currentsens_den = (Rvecmat((allppcsnem + i), j) * vtplus1.as_row() * wtplus1); // This is the denominator of the key matrix equation
          double cd_double = currentsens_den(0,0);
          arma::mat currentsens = currentsens_num / (cd_double * theclairvoyant);
          
          if (style == 1) {
            // This is the sensitivity matrix
            senscube.slice(allppcsnem + i) += currentsens; 
            
            if (historical) {
              wprojection_ah.zeros();
              vprojection_ah.zeros();
              
              // This loop creates the ahistorical stable stage distribution for projected
              // occasion j+1
              for (int k1 = 0; k1 < hstages_num; k1++) {
                int current_stage2 = hstages_id2(k1);
                wprojection_ah(current_stage2 - 1) = wprojection_ah(current_stage2 - 1)  +
                  wtplus1(k1);
              } // k1 loop
              
              // Now the ahistorical reproductive value vector for occasion j+1
              for (int k2 = 0; k2 < hstages_num; k2++) {
                int current_stage2 = hstages_id2(k2);
                
                if (wprojection_ah(current_stage2 - 1) > 0) {
                  vprojection_ah(current_stage2 - 1) = vprojection_ah(current_stage2 - 1) +
                    (vtplus1(k2) * wtplus1(k2) / wprojection_ah(current_stage2 - 1));
                }
              } // k2 loop
              
              // Now to propagate the projection sensitivity matrix, and add it to
              // the main sensitivity matrix
              arma::rowvec wtah_tpose = wprojection_ah.as_row();
              arma::rowvec vtah_tpose = vprojection_ah.as_row();
              arma::mat csah_num = vprojection_ah * wtah_tpose;
              arma::mat csah_den = (Rvecmat(i, j) * vtah_tpose * wprojection_ah);
              double cdah_double = csah_den(0,0);
              arma::mat csah = csah_num / (cdah_double * theclairvoyant);
              senscube_ah.slice(i) += csah;
              
            } // if historical statement
          } else {
            // This is the elasticity matrix
            senscube.slice(allppcsnem + i) += currentsens % as<arma::mat>(meanmatyearlist[(theprophecy(j))]); 
          }
        }
      } // for loop i, for populations
    } // if statement, checking if more than one patch and thus determining if population means need to be dealt with
    
    if (historical && style == 2) {
      for (int k = 0; k < trials; k++) {
        arma::uvec hstages_id1 = hstages["stage_id_1"];
        arma::uvec hstages_id2 = hstages["stage_id_2"];
        arma::mat elasah(ahstages_num, ahstages_num);
        elasah.zeros();
        
        arma::mat hslice = senscube.slice(k);
        
        for (int i = 0; i < hstages_num; i++) {
          for (int j = 0; j < hstages_num; j++) {
            elasah((hstages_id2(j) - 1), (hstages_id1(i) - 1)) = elasah((hstages_id2(j) - 1), (hstages_id1(i) - 1)) + hslice(j, i);
          }
        }
        senscube_ah.slice(k) = elasah;
      }
    }

    return Rcpp::List::create(_["maincube"] = senscube, _["ahcube"] = senscube_ah);
    
  } else { 
    // This chunk focuses on the core population in cases where no patch and year
    // data is given, but a single list of A matrices is provided
    
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
    
    IntegerVector uniqueyears = seq(0, (yl - 1));
    arma::uvec uniqueyears_arma = as<arma::uvec>(uniqueyears);
    arma::vec twinput;
    
    if (matrows != matcols) {
      throw Rcpp::exception("Supplied matrices must be square. Please check matrix dimensions and fix.", false);
    }
    
    if (tweights.isNotNull()) {
      if (as<NumericVector>(tweights).length() != yl) {
        throw Rcpp::exception("Time weight vector must be the same length as the number of occasions represented in the lefkoMat object used as input.", false);
      }
      twinput = as<arma::vec>(tweights);
    } else {
      twinput.resize(yl);
      twinput.ones();
    } // At the end of this, we have an arma::vec holding with whatever the user supplied, or a 1 for each occasion
    
    // Here we set up the vector of chosen occasions, sampled from all possible occasions
    arma::vec twinput_corr = twinput / sum(twinput);
    arma::uvec theprophecy = Rcpp::RcppArmadillo::sample(uniqueyears_arma, theclairvoyant, true, twinput_corr);
    
    // Here we initialize a core empty matrix and start vector for w and v calculations.
    // The matrix will be changed at each occasion.
    arma::vec startvec(matrows);
    startvec.ones();
    startvec = startvec / matrows; // The is the start vector for w and v calculations
    int trials {1};
    
    // Here we initialize a flat cube to hold the sensitivity or elasticity matrix
    arma::cube senscube(matrows, matrows, trials);
    senscube.zeros();
    
    // These matrices and vectors will hold R values
    arma::mat Rvecmat(trials, theclairvoyant);
    Rvecmat.zeros();
    
    arma::mat wprojection(startvec.n_elem, (theclairvoyant + 1));
    arma::mat vprojection(startvec.n_elem, (theclairvoyant + 1));
    wprojection.zeros();
    vprojection.zeros();
    
    // Here we run the control loop to create the w and v values we need
    arma::mat crazy_prophet = proj3(startvec, amats, theprophecy, 1, 0, 0);
    wprojection = crazy_prophet.submat(startvec.n_elem, 0, ((startvec.n_elem * 2) - 1), theclairvoyant);
    vprojection = crazy_prophet.submat((startvec.n_elem * 2), 0, ((startvec.n_elem * 3) - 1), theclairvoyant);
    Rvecmat.row(0) = crazy_prophet.submat((startvec.n_elem * 3), 1, (startvec.n_elem * 3), theclairvoyant); // Rvec
    
    // All references should go to senscube, which is a 3d array designed to hold the sensitivity matrices
    for (int j = 0; j < theclairvoyant; j++) { // This is the main time loop for the sensitivity matrices, 
                                               // adding each occasion to the respective matrix for each pop-patch
      arma::vec vtplus1 = vprojection.col(j+1); // used to be j+1
      arma::vec wt = wprojection.col(j);
      
      arma::mat currentsens_num = vtplus1 * wt.as_row(); // This is the numerator of the key matrix equation
      arma::mat currentsens_den = (Rvecmat(0, j) * vtplus1.as_row() * wprojection.col(j+1)); // This is the denominator of the key matrix equation
      double cd_double = currentsens_den(0,0);
      arma::mat currentsens = currentsens_num / (cd_double * theclairvoyant);
      
      if (style == 1) {
        senscube.slice(0) += currentsens; // This is the sensitivity matrix
      } else {
        senscube.slice(0) += currentsens % as<arma::mat>(amats[(theprophecy(j))]); // This is the elasticity matrix
      }
    }
    
  return Rcpp::List::create(_["maincube"] = senscube);    
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
//' \item{size3}{Size in occasion \emph{t}+1.}
//' \item{repstatus3}{Reproductive status in occasion \emph{t}+1.}
//' \item{entrystatus3}{Entry status in occasion \emph{t}+1.}
//' \item{size2}{Size in occasion \emph{t}.}
//' \item{repstatus2}{Reproductive status in occasion \emph{t}.}
//' \item{entrystatus2}{Entry status in occasion \emph{t}.}
//' \item{size1}{Size in occasion \emph{t}-1.}
//' \item{repstatus1}{Reproductive status in occasion \emph{t}11.}
//' \item{entrystatus1}{Entry status in occasion \emph{t}-1.}
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
// [[Rcpp::export(.bambi3)]]
DataFrame bambi3(DataFrame stages, DataFrame hstages) {
  
  StringVector stagenames = stages["stage"];
  arma::uvec astages = stages["stage_id"];
  arma::vec sizes = stages["original_size"];
  arma::uvec repstatus = stages["repstatus"];
  arma::uvec entrystage = stages["entrystage"];
  int numstages = astages.n_elem;
  
  arma::uvec hstage3in = hstages["stage_id_2"];
  arma::uvec hstage2nin = hstages["stage_id_1"];
  int numhstages = hstage3in.n_elem;
  
  hstage3in = hstage3in - 1;
  hstage2nin = hstage2nin - 1;
  
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
      if (hstage3in(i1) == (hstage2nin(i2))) { // Originally this was: if (hstage3in(i1) == (hstage2nin(i2) + 1))
        
        hsindexl(counter) = (i1 * numhstages) + i2;
        
        int stage1 = hstage2nin(i2);
        longnames1(counter) = stagenames(stage1);
        size1(counter) = sizes(stage1);
        repstatus1(counter) = repstatus(stage1);
        entrystatus1(counter) = entrystage(stage1);
        
        int stage2 = hstage2nin(i1);
        longnames2(counter) = stagenames(stage2);
        size2(counter) = sizes(stage2);
        repstatus2(counter) = repstatus(stage2);
        entrystatus2(counter) = entrystage(stage2);
        
        int stage3 = hstage3in(i2);
        longnames3(counter) = stagenames(stage3);
        size3(counter) = sizes(stage3);
        repstatus3(counter) = repstatus(stage3);
        entrystatus3(counter) = entrystage(stage3);
        
        if (entrystatus3(counter) == 1 && repstatus2(counter) == 1) {
          if (entrystatus2(counter) == 1 && repstatus1(counter) == 1) {
            transition_type(counter) = 26; // Fecundity to fecundity
          } else if (size2(counter) == size1(counter)) {
            if (repstatus2(counter) > repstatus1(counter) || entrystatus3(counter) < entrystatus2(counter)) {
              transition_type(counter) = 21; // Growth to fecundity
            } else if (repstatus2(counter) < repstatus1(counter) || entrystatus2(counter) > entrystatus1(counter)) {
              transition_type(counter) = 22; // Shrinkage to fecundity
            } else {
              transition_type(counter) = 20; // Stasis to fecundity
            }
          } else if (size2(counter) > size1(counter)) {
            transition_type(counter) = 21; // Growth to fecundity
          } else if (size2(counter) < size1(counter)) {
            transition_type(counter) = 22; // Shrinkage to fecundity
          }
        } else if (entrystatus2(counter) == 1 && repstatus1(counter) == 1) {
          if (size3(counter) == size2(counter)) {
            if (repstatus3(counter) > repstatus2(counter) || entrystatus3(counter) < entrystatus2(counter)) {
              transition_type(counter) = 24; // Fecundity to growth
            } else if (repstatus3(counter) < repstatus2(counter) || entrystatus3(counter) > entrystatus2(counter)) {
              transition_type(counter) = 25; // Fecundity to shrinkage
            } else {
              transition_type(counter) = 23; // Fecundity to stasis
            }
          } else if (size3(counter) > size2(counter)) {
            transition_type(counter) = 24; // Fecundity to growth
          } else if (size3(counter) < size2(counter)) {
            transition_type(counter) = 25; // Fecundity to shrinkage
          }
        } else if (size3(counter) == size2(counter) && size2(counter) == size1(counter)) {
          if (repstatus2(counter) > repstatus1(counter) || entrystatus2(counter) < entrystatus1(counter)) {
            if (repstatus3(counter) > repstatus2(counter) || entrystatus3(counter) < entrystatus2(counter)) {
              transition_type(counter) = 12; // Full growth
            } else if (repstatus3(counter) < repstatus2(counter) || entrystatus3(counter) > entrystatus2(counter)) {
              transition_type(counter) = 17; // Growth to shrinkage
            } else {
              transition_type(counter) = 13; // Growth to stasis
            }
          } else if (repstatus2(counter) < repstatus1(counter)) {
            if (repstatus3(counter) > repstatus2(counter) || entrystatus3(counter) < entrystatus2(counter)) {
              transition_type(counter) = 18; // Shrinkage to growth
            } else if (repstatus3(counter) < repstatus2(counter) || entrystatus3(counter) > entrystatus2(counter)) {
              transition_type(counter) = 15; // Full shrinkage
            } else {
              transition_type(counter) = 16; // Shrinkage to stasis
            }
          } else {
            if (repstatus3(counter) > repstatus2(counter) || entrystatus3(counter) < entrystatus2(counter)) {
              transition_type(counter) = 11; // Stasis to growth
            } else if (repstatus3(counter) < repstatus2(counter) || entrystatus3(counter) > entrystatus2(counter)) {
              transition_type(counter) = 14; // Stasis to shrinkage
            } else {
              transition_type(counter) = 10; // Full stasis
            }
          }
        } else if (size3(counter) > size2(counter) && size2(counter) == size1(counter)) {
          if (repstatus2(counter) > repstatus1(counter) || entrystatus2(counter) < entrystatus1(counter)) {
            transition_type(counter) = 12; // Full growth
          } else if (repstatus2(counter) < repstatus1(counter) || entrystatus2(counter) > entrystatus1(counter)) {
            transition_type(counter) = 18; // Shrinkage to growth
          } else {
            transition_type(counter) = 11; // Stasis to growth
          }
        } else if (size3(counter) > size2(counter) && size2(counter) > size1(counter)) {
          transition_type(counter) = 12; // Full growth
        } else if (size3(counter) == size2(counter) && size2(counter) > size1(counter)) {
          if (repstatus3(counter) > repstatus2(counter) || entrystatus3(counter) < entrystatus2(counter)) {
            transition_type(counter) = 12; // Full growth
          } else if (repstatus3(counter) < repstatus2(counter) || entrystatus3(counter) > entrystatus2(counter)) {
            transition_type(counter) = 17; // Growth to shrinkage
          } else {
            transition_type(counter) = 13; // Growth to stasis
          }
        } else if (size3(counter) < size2(counter) && size2(counter) == size1(counter)) {
          if (repstatus2(counter) > repstatus1(counter) || entrystatus2(counter) < entrystatus1(counter)) {
            transition_type(counter) = 17; // Growth to shrinkage
          } else if (repstatus2(counter) < repstatus1(counter) || entrystatus2(counter) > entrystatus1(counter)) {
            transition_type(counter) = 15; // Full shrinkage
          } else {
            transition_type(counter) = 14; // Stasis to shrinkage
          }
        } else if (size3(counter) < size2(counter) && size2(counter) < size1(counter)) {
          transition_type(counter) = 15; // Full shrinkage
        } else if (size3(counter) == size2(counter) && size2(counter) < size1(counter)) {
          if (repstatus3(counter) > repstatus2(counter) || entrystatus3(counter) < entrystatus2(counter)) {
            transition_type(counter) = 18; // Shrinkage to growth
          } else if (repstatus3(counter) < repstatus2(counter) || entrystatus3(counter) > entrystatus2(counter)) {
            transition_type(counter) = 15; // Full shrinkage
          } else {
            transition_type(counter) = 16; // Shrinkage to stasis
          }
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
//' \item{stage3}{Stage in occasion \emph{t}+1.}
//' \item{size3}{Size in occasion \emph{t}+1.}
//' \item{repstatus3}{Reproductive status in occasion \emph{t}+1.}
//' \item{entrystatus3}{Entry status in occasion \emph{t}+1.}
//' \item{stage2}{Stage in occasion \emph{t}.}
//' \item{size2}{Size in occasion \emph{t}.}
//' \item{repstatus2}{Reproductive status in occasion \emph{t}.}
//' \item{entrystatus2}{Entry status in occasion \emph{t}.}
//'
//' The kind of transitions conforms to the following code: \code{1}: stasis, 
//' \code{2}: growth, \code{3}: shrinkage, \code{4}: fecundity.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.bambi2)]]
DataFrame bambi2(DataFrame stages) {
  
  StringVector stagenames = stages["stage"];
  arma::uvec astages = stages["stage_id"];
  arma::vec sizes = stages["original_size"];
  arma::uvec repstatus = stages["repstatus"];
  arma::uvec entrystage = stages["entrystage"];
  int numstages = astages.n_elem;
  
  astages = astages - 1;
  
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
      
      int stage2 =  astages(i1);
      longstages2(counter) = stagenames(stage2);
      size2(counter) = sizes(stage2);
      repstatus2(counter) = repstatus(stage2);
      entrystatus2(counter) = entrystage(stage2);
      
      int stage3 = astages(i2);
      longstages3(counter) = stagenames(stage3);
      size3(counter) = sizes(stage3);
      repstatus3(counter) = repstatus(stage3);
      entrystatus3(counter) = entrystage(stage3);
      
      if (entrystatus3(counter) == 1 && repstatus2(counter) == 1) {
        transition_type(counter) = 4; // Fecundity
      } else if (size3(counter) == size2(counter)) {
        if (repstatus3(counter) > repstatus2(counter) || entrystatus3(counter) < entrystatus2(counter)) {
          transition_type(counter) = 2; // Growth
        } else if (repstatus3(counter) < repstatus2(counter) || entrystatus3(counter) > entrystatus2(counter)) {
          transition_type(counter) = 3; // Shrinkage
        } else {
          transition_type(counter) = 1; // Stasis
        }
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
//' matrices, and LTRE contributions from LTRE and sLTRE matrices, according to
//' the categories developed by functions \code{bambi2()} and \code{bambi3()}.
//' 
//' @param e_amat A single elasticity, LTRE, or sLTRE matrix.
//' @param bambesque This is the output from \code{bambi2()} or \code{bambi3()}
//' corresponding to the current lefkoMat object. The format is a data frame
//' giving the indices and characteristics of all predicted potential non-zero
//' elements in the supplied matrix.
//' @param amat_ The A matrix corresponding to \code{e_amat}. If not supplied,
//' then only \code{bambesque} is used to determine transition categories. If
//' provided, then fecundity transitions may be split between fecundity and
//' survival portions.
//' @param fmat_ The F matrix corresponding to \code{e_amat}. If not supplied,
//' then only \code{bambesque} is used to determine transition categories. If
//' provided, then fecundity transitions may be split between fecundity and
//' survival portions.
//' 
//' @return A list with two data frames, one showing the summed elasticities for
//' the historical matrix supplied (if supplied), and the other showing the
//' ahistorical summary of the historical matrix or the summed elasticities of
//' a supplied ahistorical elasticity matrix. Also includes sums of only the
//' positive elements and only the negative elements, in all cases.
//' 
//' @section Notes:
//' If the original matrices are provided, then this function was made to split
//' co-occurring survival-fecundity elasticities according to the ratio of the
//' fecundity portion of the element to the survival portion of that element.
//' However, this transition splitting capability developed using the original
//' matrices does not currently work properly, and so it is better to use this
//' function without objects \code{amat_} and \code{fmat_}, forcing co-occurring
//' survival-fecundity transitions to be treated as fecundity only.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.demolition3)]]
List demolition3(arma::mat e_amat, DataFrame bambesque,
  Nullable<Rcpp::NumericVector> amat_ = R_NilValue,
  Nullable<Rcpp::NumericVector> fmat_ = R_NilValue) {
  
  arma::uvec eindices = bambesque["index"];
  arma::uvec categories = bambesque["transition"];
  
  int e_amatsize = e_amat.n_elem;
  int e_amatrows = e_amat.n_rows;
  int maxelem = static_cast<int>(eindices.max());
  int minindex = static_cast<int>(categories.min());
  
  arma::mat amat(e_amatrows, e_amatrows);
  arma::mat fmat(e_amatrows, e_amatrows);
  
  if (maxelem > e_amatsize) {
    throw Rcpp::exception("Supplied info does not seem to correspond to current matrix inputs.", false);
  }
  
  if (amat_.isNotNull() && fmat_.isNotNull()) {
    amat = Rcpp::as<arma::mat>(amat_);
    fmat = Rcpp::as<arma::mat>(fmat_);
  } else {
    amat.ones();
    fmat.zeros();
    
    arma::uvec fec_trans = find(categories == 4);
    if (fec_trans.n_elem > 0) {
      for (int i = 0; i < fec_trans.n_elem; i ++) {
        fmat(eindices(fec_trans(i))) = 1;
      }
    }
    fec_trans = find(categories == 20);
    if (fec_trans.n_elem > 0) {
      for (int i = 0; i < fec_trans.n_elem; i ++) {
        fmat(eindices(fec_trans(i))) = 1;
      }
    }
    fec_trans = find(categories == 21);
    if (fec_trans.n_elem > 0) {
      for (int i = 0; i < fec_trans.n_elem; i ++) {
        fmat(eindices(fec_trans(i))) = 1;
      }
    }
    fec_trans = find(categories == 22);
    if (fec_trans.n_elem > 0) {
      for (int i = 0; i < fec_trans.n_elem; i ++) {
        fmat(eindices(fec_trans(i))) = 1;
      }
    }
    fec_trans = find(categories == 26);
    if (fec_trans.n_elem > 0) {
      for (int i = 0; i < fec_trans.n_elem; i ++) {
        fmat(eindices(fec_trans(i))) = 1;
      }
    }
  }
  
  // The next portion allows fecundity transitions to be split, if they include survival portions
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
  arma::vec histpos(16);
  arma::vec histneg(16);
  arma::vec hc_ahistsums(4);
  arma::vec hc_ahistpos(4);
  arma::vec hc_ahistneg(4);
  histsums.zeros();
  histpos.zeros();
  histneg.zeros();
  hc_ahistsums.zeros();
  hc_ahistpos.zeros();
  hc_ahistneg.zeros();
  
  StringVector ahistcats {"Stasis", "Growth", "Shrinkage", "Fecundity"};
  arma::uvec ahistcatnums {1, 2, 3, 4};
  arma::vec ahistsums(4);
  arma::vec ahistpos(4);
  arma::vec ahistneg(4);
  ahistsums.zeros();
  ahistpos.zeros();
  ahistneg.zeros();
  
  if (minindex > 9) {
    // Object minindex will only be above 9 if the supplied object is historical
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
          int this_guy = eindices(currentguys(j));
          
          if (fec_fraction(this_guy) == 1) {
            hc_ahistsums(3) += (e_amat(this_guy));
            histsums(i) += (e_amat(this_guy));
            
            if (e_amat(this_guy) > 0) {
              hc_ahistpos(3) += (e_amat(this_guy));
              histpos(i) += (e_amat(this_guy));
            } else if (e_amat(this_guy) < 0) {
              hc_ahistneg(3) += (e_amat(this_guy));
              histneg(i) += (e_amat(this_guy));
            }
          } else {
            hc_ahistsums(3) += (e_amat(this_guy) * fec_fraction(this_guy));
            histsums(i) += (e_amat(this_guy) * fec_fraction(this_guy));
            
            if (e_amat(this_guy) > 0) {
              hc_ahistpos(3) += (e_amat(this_guy) * (fec_fraction(this_guy)));
              histpos(i) += (e_amat(this_guy) * (fec_fraction(this_guy)));
            } else if (e_amat(this_guy) < 0) {
              hc_ahistneg(3) += (e_amat(this_guy) * (fec_fraction(this_guy)));
              histneg(i) += (e_amat(this_guy) * (fec_fraction(this_guy)));
            }
            
            arma::uvec counter = find(eindices == this_guy);
            
            if (entrystatus2(counter(0)) == 1 && repstatus1(counter(0)) == 1) {
              if (size3(counter(0)) == size2(counter(0)) && repstatus3(counter(0)) == repstatus2(counter(0)) &&
                  entrystatus3(counter(0)) == entrystatus2(counter(0))) {
                hc_ahistsums(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(12) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(12) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(12) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              } else if (size3(counter(0)) > size2(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(13) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(13) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(13) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              } else if (size3(counter(0)) < size2(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(14) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(14) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(14) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              } else if (repstatus3(counter(0)) > repstatus2(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(13) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(13) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(13) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              } else if (repstatus3(counter(0)) < repstatus2(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(14) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(14) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(14) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              } else if (entrystatus3(counter(0)) < entrystatus2(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(13) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(13) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(13) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              } else if (entrystatus3(counter(0)) > entrystatus2(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(14) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(14) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(14) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              }
            } else if (size3(counter(0)) == size2(counter(0)) && size2(counter(0)) == size1(counter(0))) {
              if (repstatus3(counter(0)) > repstatus2(counter(0))) {
                if (repstatus2(counter(0)) > repstatus1(counter(0))) {
                  hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                } else if (repstatus2(counter(0)) < repstatus1(counter(0))) {
                  hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                } else if (entrystatus2(counter(0)) < entrystatus1(counter(0))) {
                  hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                } else if (entrystatus2(counter(0)) > entrystatus1(counter(0))) {
                  hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                } else {
                  hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                }
              } else if (repstatus3(counter(0)) < repstatus2(counter(0))) {
                if (repstatus2(counter(0)) > repstatus1(counter(0))) {
                  hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                } else if (repstatus2(counter(0)) < repstatus1(counter(0))) {
                  hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                } else if (entrystatus2(counter(0)) < entrystatus1(counter(0))) {
                  hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                } else if (entrystatus2(counter(0)) > entrystatus1(counter(0))) {
                  hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                } else {
                  hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(4) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(4) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(4) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                }
              } else if (entrystatus3(counter(0)) < entrystatus2(counter(0))) {
                if (repstatus2(counter(0)) > repstatus1(counter(0))) {
                  hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                } else if (repstatus2(counter(0)) < repstatus1(counter(0))) {
                  hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                } else if (entrystatus2(counter(0)) < entrystatus1(counter(0))) {
                  hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                } else if (entrystatus2(counter(0)) > entrystatus1(counter(0))) {
                  hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                } else {
                  hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                }
              } else if (entrystatus3(counter(0)) > entrystatus2(counter(0))) {
                if (repstatus2(counter(0)) > repstatus1(counter(0))) {
                  hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                } else if (repstatus2(counter(0)) < repstatus1(counter(0))) {
                  hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                } else if (entrystatus2(counter(0)) < entrystatus1(counter(0))) {
                  hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                } else if (entrystatus2(counter(0)) > entrystatus1(counter(0))) {
                  hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                } else {
                  hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histsums(4) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                  if (e_amat(this_guy) > 0) {
                    hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histpos(4) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  } else if (e_amat(this_guy) < 0) {
                    hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                    histneg(4) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  }
                }
              } else {
                hc_ahistsums(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              }
            } else if (size3(counter(0)) > size2(counter(0)) && size2(counter(0)) == size1(counter(0))) {
              if (repstatus2(counter(0)) > repstatus1(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              } else if (repstatus2(counter(0)) < repstatus1(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              } else if (entrystatus2(counter(0)) < entrystatus1(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              } else if (entrystatus2(counter(0)) > entrystatus1(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              } else {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              }
            } else if (size3(counter(0)) > size2(counter(0)) && size2(counter(0)) > size1(counter(0))) {
              hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              histsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              
              if (e_amat(this_guy) > 0) {
                hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              } else if (e_amat(this_guy) < 0) {
                hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              }
            } else if (size3(counter(0)) == size2(counter(0)) && size2(counter(0)) > size1(counter(0))) {
              if (repstatus3(counter(0)) > repstatus2(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              } else if (repstatus3(counter(0)) < repstatus2(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              } else if (entrystatus3(counter(0)) < entrystatus2(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              } else if (entrystatus3(counter(0)) > entrystatus2(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              } else {
                hc_ahistsums(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              }
            } else if (size3(counter(0)) < size2(counter(0)) && size2(counter(0)) == size1(counter(0))) {
              if (repstatus2(counter(0)) > repstatus1(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              } else if (repstatus2(counter(0)) < repstatus1(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              } else if (entrystatus2(counter(0)) < entrystatus1(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              } else if (entrystatus2(counter(0)) > entrystatus1(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              } else {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(4) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(4) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(4) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              }
            } else if (size3(counter(0)) < size2(counter(0)) && size2(counter(0)) < size1(counter(0))) {
              hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              histsums(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              
              if (e_amat(this_guy) > 0) {
                hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histpos(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              } else if (e_amat(this_guy) < 0) {
                hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histneg(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              }
            } else if (size3(counter(0)) == size2(counter(0)) && size2(counter(0)) < size1(counter(0))) {
              if (repstatus3(counter(0)) > repstatus2(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              } else if (repstatus3(counter(0)) < repstatus2(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              } else if (entrystatus3(counter(0)) < entrystatus2(counter(0))) {
                hc_ahistsums(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(1) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(8) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              } else if (entrystatus3(counter(0)) > entrystatus2(counter(0))) {
                hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              } else {
                hc_ahistsums(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histsums(6) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                
                if (e_amat(this_guy) > 0) {
                  hc_ahistpos(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histpos(86) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                } else if (e_amat(this_guy) < 0) {
                  hc_ahistneg(0) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                  histneg(6) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                }
              }
            } else if (size3(counter(0)) < size2(counter(0)) && size2(counter(0)) > size1(counter(0))) {
              hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              histsums(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              
              if (e_amat(this_guy) > 0) {
                hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histpos(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              } else if (e_amat(this_guy) < 0) {
                hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histneg(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              }
            } else if (size3(counter(0)) > size2(counter(0)) && size2(counter(0)) < size1(counter(0))) {
              hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              histsums(7) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
            }
          }
        }
      } else if (histcatnums(i) == 14 || histcatnums(i) == 15 || histcatnums(i) == 17 || histcatnums(i) == 25) { // Shrinkage transitions
        
        arma::vec all_es = e_amat.elem(eindices(currentguys));
        arma::uvec all_es_pos = find(all_es > 0);
        arma::uvec all_es_neg = find(all_es < 0);
        
        int all_es_pos_num = all_es_pos.n_elem;
        int all_es_neg_num = all_es_neg.n_elem;
        
        double getoutofdodge = sum(all_es);
        histsums(i) += getoutofdodge;
        hc_ahistsums(2) += getoutofdodge;
        
        if (all_es_pos_num > 0) {
          histpos(i) += sum(all_es.elem(all_es_pos));
          hc_ahistpos(2) += sum(all_es.elem(all_es_pos));
        }
        if (all_es_neg_num > 0) {
          histneg(i) += sum(all_es.elem(all_es_neg));
          hc_ahistneg(2) += sum(all_es.elem(all_es_neg));
        }
        
      } else if (histcatnums(i) == 10 || histcatnums(i) == 13 || histcatnums(i) == 16 || histcatnums(i) == 23) { // Stasis transitions
        
        arma::vec all_es = e_amat.elem(eindices(currentguys));
        arma::uvec all_es_pos = find(all_es > 0);
        arma::uvec all_es_neg = find(all_es < 0);
        
        int all_es_pos_num = all_es_pos.n_elem;
        int all_es_neg_num = all_es_neg.n_elem;
        
        double getoutofdodge = sum(all_es);
        histsums(i) += getoutofdodge;
        hc_ahistsums(0) += getoutofdodge;
        
        if (all_es_pos_num > 0) {
          histpos(i) += sum(all_es.elem(all_es_pos));
          hc_ahistpos(0) += sum(all_es.elem(all_es_pos));
        }
        if (all_es_neg_num > 0) {
          histneg(i) += sum(all_es.elem(all_es_neg));
          hc_ahistneg(0) += sum(all_es.elem(all_es_neg));
        }
        
      } else if (histcatnums(i) == 11 || histcatnums(i) == 12 || histcatnums(i) == 18 || histcatnums(i) == 24) { // Growth transitions
        
        arma::vec all_es = e_amat.elem(eindices(currentguys));
        arma::uvec all_es_pos = find(all_es > 0);
        arma::uvec all_es_neg = find(all_es < 0);
        
        int all_es_pos_num = all_es_pos.n_elem;
        int all_es_neg_num = all_es_neg.n_elem;
        
        double getoutofdodge = sum(all_es);
        histsums(i) += getoutofdodge;
        hc_ahistsums(1) += getoutofdodge;
        
        if (all_es_pos_num > 0) {
          histpos(i) += sum(all_es.elem(all_es_pos));
          hc_ahistpos(1) += sum(all_es.elem(all_es_pos));
        }
        if (all_es_neg_num > 0) {
          histneg(i) += sum(all_es.elem(all_es_neg));
          hc_ahistneg(1) += sum(all_es.elem(all_es_neg));
        }
      }
    }
    
    histout = DataFrame::create(Named("category") = histcats, _["elas"] = histsums,
      _["elas_pos"] = histpos, _["elas_neg"] = histneg);
    ahistout = DataFrame::create(Named("category") = ahistcats, _["elas"] = hc_ahistsums,
      _["elas_pos"] = hc_ahistpos, _["elas_neg"] = hc_ahistneg);
  } else {
    histout = R_NilValue;
    
    for (int i = 0; i < 4; i++) {
      arma::uvec currentguys = find(categories == ahistcatnums(i));
      
      arma::vec all_es = e_amat.elem(eindices(currentguys));
      arma::uvec all_es_pos = find(all_es > 0);
      arma::uvec all_es_neg = find(all_es < 0);
      
      int all_es_pos_num = all_es_pos.n_elem;
      int all_es_neg_num = all_es_neg.n_elem;
      
      double getoutofdodge = sum(all_es);
      ahistsums(i) += getoutofdodge;
      
      if (all_es_pos_num > 0) {
        ahistpos(i) += sum(all_es.elem(all_es_pos));
      }
      if (all_es_neg_num > 0) {
        ahistneg(i) += sum(all_es.elem(all_es_neg));
      }
    }
    
    ahistout = DataFrame::create(Named("category") = ahistcats, _["elas"] = ahistsums,
      _["elas_pos"] = ahistpos, _["elas_neg"] = ahistneg);
  }
  
  List output = List::create(Named("hist") = histout, _["ahist"] = ahistout);
  
  return output;
}

//' Estimate LTRE of Any Population Matrix
//' 
//' \code{ltre3matrix()} returns the one-way fixed deterministic LTRE matrix of
//' a dense or sparse set of input matrices.
//' 
//' @param Amats A list of population projection matrices (not an entire
//' \code{lefkoMat} object.
//' @param refnum An integer vector giving the numbers of the matrices to use as
//' reference from_ \code{refmats}.
//' @param refmats_ A list of reference population projection matrices.
//' @param mean A logical value indicating whether to use the element-wise mean
//' matrix as the reference.
//' @param sparse A logical value indicating whether to use sparse or dense
//' format in matrix calculations.
//' 
//' @return This function returns a cube of LTRE contributions, with each slice
//' corresponding to each input matrix in Amats. 
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.ltre3matrix)]]
arma::cube ltre3matrix(List Amats, Rcpp::IntegerVector refnum,
  Nullable<Rcpp::List> refmats_ = R_NilValue, bool mean = true,
  bool sparse = false) {
  
  int Amatnum = Amats.length();
  int matdim = as<arma::mat>(Amats(0)).n_rows;
  
  List eigenstuff;
  List mean_set;
  
  arma::cube senscube(matdim, matdim, Amatnum);
  senscube.zeros();
  
  if (!sparse) {
    // Dense matrix analysis
    arma::mat finalrefmat(matdim, matdim);
    finalrefmat.zeros();
      
    int refmatnum = refnum.length();
    
    if (refmats_.isNotNull()) {
      Rcpp::List refmats(refmats_);
      
      if (mean && refmatnum > 1) {
        for (int i = 0; i < refmatnum; i++) {
          finalrefmat = finalrefmat + (as<arma::mat>(refmats(refnum(i))) / static_cast<double>(refmatnum));
        }
      } else {
        finalrefmat = as<arma::mat>(refmats(refnum(0)));
      }
    } else {
      if (mean && refmatnum > 1) {
        for (int i = 0; i < Amatnum; i++) {
          finalrefmat = finalrefmat + (as<arma::mat>(Amats(refnum(i))) / static_cast<double>(Amatnum));
        }
      } else {
        finalrefmat = as<arma::mat>(Amats(refnum(0)));
      }
    }
    
    // Now we create the halfway matrices and run the sensitivities
    arma::mat halfmat(matdim, matdim);
    arma::mat diffmat(matdim, matdim);
    halfmat.zeros();
    diffmat.zeros();
    
    for (int i = 0; i < Amatnum; i++) {
      halfmat = (finalrefmat + (as<arma::mat>(Amats(i)))) / static_cast<double>(2.0);
      diffmat = (as<arma::mat>(Amats(i))) - finalrefmat;
      
      senscube.slice(i) = diffmat % sens3matrix(halfmat, 0);
    }
    
    return senscube;
    
  } else {
    // Sparse matrix analysis
    arma::sp_mat finalrefmat(matdim, matdim);
    finalrefmat.zeros();
      
    int refmatnum = refnum.length();
    
    if (refmats_.isNotNull()) {
      Rcpp::List refmats(refmats_);
      
      if (mean && refmatnum > 1) {
        for (int i = 0; i < refmatnum; i++) {
          finalrefmat = finalrefmat + (arma::sp_mat(as<arma::mat>(refmats(refnum(i)))) / static_cast<double>(refmatnum));
        }
      } else {
        finalrefmat = (arma::sp_mat(as<arma::mat>(refmats(refnum(0)))));
      }
    } else {
      if (mean && refmatnum > 1) {
        for (int i = 0; i < Amatnum; i++) {
          finalrefmat = finalrefmat + (arma::sp_mat(as<arma::mat>(Amats(refnum(i)))) / static_cast<double>(Amatnum));
        }
      } else {
        finalrefmat = (arma::sp_mat(as<arma::mat>(Amats(refnum(0)))));
      }
    }
    
    // Now we create the halfway matrices and run the sensitivities
    arma::sp_mat halfmat(matdim, matdim);
    arma::sp_mat diffmat(matdim, matdim);
    halfmat.zeros();
    diffmat.zeros();
    
    for (int i = 0; i < Amatnum; i++) {
      halfmat = (finalrefmat + (arma::sp_mat(as<arma::mat>(Amats(i))))) / static_cast<double>(2.0);
      diffmat = (arma::sp_mat(as<arma::mat>(Amats(i)))) - finalrefmat;
      
      senscube.slice(i) = arma::mat(diffmat % sens3sp_matrix(halfmat, diffmat));
    }
    
    return senscube;
  }
}

//' Estimate sLTRE of Any Population Matrix
//' 
//' \code{sltre3matrix()} returns the one-way stochastic LTRE matrix of
//' a dense or sparse set of input matrices.
//' 
//' @param Amats A list of population projection matrices (not an entire
//' \code{lefkoMat} object).
//' @param loy A data frame showing the order of populations, patches, and
//' occasions of the matrices provided in object \code{Amats}.
//' @param refnum An integer vector giving the numbers of the matrices to use as
//' reference from \code{refmats}.
//' @param refmats_ A list of reference population projection matrices.
//' @param tweights_ Numeric vector denoting the probabilistic weightings of
//' annual matrices. Defaults to equal weighting among occasions.
//' @param steps The number of occasions to project the stochastic simulation
//' forward, if performing an sLTRE. Defaults to 10,000. Note that the total
//' number of occasions projected equals this number plus the number given in
//' object \code{burnin}.
//' @param burnin The number of initial occasions to project the population
//' without calculating population metrics. Defaults to 3000.
//' @param sparse A logical value indicating whether to use sparse or dense
//' format in matrix calculations.
//' 
//' @return This function returns a list of two lists of matrices. The first,
//' \code{cont_mean}, holds the sLTRE contributions of shifts in mean elements.
//' The second, \code{cont_sd}, holds the sLTRE contributions of shifts in
//' temporal standard deviations of matrix elements.
//' 
//' @section Notes:
//' This function uses the simulation approach developed in Davison et al.
//' (2010), which is a good approximation though not an analytical solution.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.sltre3matrix)]]
Rcpp::List sltre3matrix(List Amats, DataFrame loy, Rcpp::IntegerVector refnum,
  Nullable<Rcpp::List> refmats_ = R_NilValue,
  Nullable<arma::vec> tweights_ = R_NilValue, int steps = 10000,
  int burnin = 3000, bool sparse = false) {
  
  int theclairvoyant = steps + burnin;
  arma::vec tweights;
  
  int Amatnum = Amats.length();
  int matdim = as<arma::mat>(Amats(0)).n_rows;
  int matlength = as<arma::mat>(Amats(0)).n_elem;
  int refmatnum = refnum.length();
  
  List eigenstuff;
  List mean_set;
  
  List poppatch_meanmat;
  List poppatch_sdmat;
  List ref_byyear;
  List elas_mean;
  List elas_sd;
  List diff_meanmat;
  List diff_sdmat;
  List cont_meanmat;
  List cont_sdmat;
  
  // This section creates the order vectors for pop-patch-year
  StringVector pops = loy["pop"];
  arma::uvec pop_num = loy["popc"];
  StringVector patches = loy["patch"];
  arma::uvec year2c = loy["year2c"];
  arma::uvec poppatchc = loy["poppatchc"];
  arma::uvec patchesinpop = loy["patchesinpop"];
  arma::uvec yearsinpatch = loy["yearsinpatch"];
  arma::uvec uniquepops = unique(pop_num);
  arma::uvec uniquepoppatches = unique(poppatchc);
  arma::uvec uniqueyears = unique(year2c);
  
  int numpoppatches = uniquepoppatches.n_elem;
  int numyears = uniqueyears.n_elem;
  
  arma::vec twinput_corr;
  arma::uvec theprophecy;
  
  // The main loop
  if (!sparse) {
    // Dense matrix analysis
    arma::mat mat_mean(matdim, matdim);
    arma::mat mat_sd(matdim, matdim);
    arma::mat ref_matmean(matdim, matdim);
    arma::mat ref_matsd(matdim, matdim);
    mat_mean.zeros();
    mat_sd.zeros();
    ref_matmean.zeros();
    ref_matsd.zeros();
    
    // First the pop/patch means and sds
    for (int i = 0; i < numpoppatches; i++) {
      arma::uvec poppatch_chosen = find(poppatchc == uniquepoppatches(i));
      int numpoppatch_chosen = poppatch_chosen.n_elem;
      
      for (int j = 0; j < numpoppatch_chosen; j++) {
        mat_mean = mat_mean + (as<arma::mat>(Amats(poppatch_chosen(j))) / static_cast<double>(numpoppatch_chosen));
      }
      if (i == 0) {
        poppatch_meanmat = List::create(mat_mean);
      } else {
        poppatch_meanmat.push_back(mat_mean);
      }
      mat_mean.zeros();
    }
    
    for (int i = 0; i < numpoppatches; i++) {
      arma::uvec poppatch_chosen = find(poppatchc == uniquepoppatches(i));
      int numpoppatch_chosen = poppatch_chosen.n_elem;
      arma::vec mat_elems(numpoppatch_chosen);
      
      for (int j = 0; j < matlength; j++) {
        mat_elems.zeros();
        
        for (int k = 0; k < numpoppatch_chosen; k++) {
          mat_elems(k) = as<arma::mat>(Amats(poppatch_chosen(k)))(j);
        }
        
        if (sum(mat_elems) != 0) {
          mat_sd(j) = stddev(mat_elems, 1);
        }
      }
      if (i == 0) {
        poppatch_sdmat = List::create(mat_sd);
      } else {
        poppatch_sdmat.push_back(mat_sd);
      }
      mat_sd.zeros();
    }
    
    // Here we create the reference matrix sets
    if (refmats_.isNotNull()) {
      Rcpp::List refmats(refmats_);
      
      mat_mean.zeros();
      mat_sd.zeros();
      
      ref_byyear = List::create(refmats);
      
      for (int i = 0; i < refmatnum; i++) {
        mat_mean = mat_mean + (as<arma::mat>(refmats(refnum(i))) / static_cast<double>(refmatnum));
      }
      
      for (int i = 0; i < matlength; i++) {
        arma::vec mat_elems(refmatnum);
        mat_elems.zeros();
        
        for (int j = 0; j < refmatnum; j++) {
          mat_elems(j) = as<arma::mat>(refmats(refnum(j)))(i);
        }
        
        if (sum(mat_elems) != 0) {
          mat_sd(i) = stddev(mat_elems, 1);
        }
      }
      ref_matmean = mat_mean;
      ref_matsd = mat_sd;
      
      mat_mean.zeros();
      mat_sd.zeros();
      
    } else {
      if (refmatnum == Amatnum) {
        
        // Reference by year
        for (int i = 0; i < numyears; i++) {
          arma::uvec year_chosen = find(year2c == uniqueyears(i));
          int numyear_chosen = year_chosen.n_elem;
          
          for (int j = 0; j < numyear_chosen; j++) {
            mat_mean = mat_mean + (as<arma::mat>(Amats(year_chosen(j))) / static_cast<double>(numyear_chosen));
          }
          if (i == 0) {
            ref_byyear = List::create(mat_mean);
          } else {
            ref_byyear.push_back(mat_mean);
          }
          mat_mean.zeros();
        }
        
        // Reference mean and sd
        for (int i = 0; i < numyears; i++) {
          mat_mean = mat_mean + (as<arma::mat>(ref_byyear(i)) / static_cast<double>(numyears));
        }
        for (int i = 0; i < matlength; i++) {
          arma::vec mat_elems(refmatnum);
          mat_elems.zeros();
          
          for (int j = 0; j < numyears; j++) {
            mat_elems(j) = as<arma::mat>(ref_byyear(j))(i);
          }
          
          if (sum(mat_elems) != 0) {
            mat_sd(i) = stddev(mat_elems, 1);
          }
        }
        ref_matmean = mat_mean;
        ref_matsd = mat_sd;
        
        mat_mean.zeros();
        mat_sd.zeros();
      } else {
        
        // Reference by year
        for (int i = 0; i < refmatnum; i++) {
          if (i == 0) {
            ref_byyear = List::create(as<arma::mat>(Amats(i)));
          } else {
            ref_byyear.push_back(as<arma::mat>(Amats(i)));
          }
          mat_mean.zeros();
          mat_sd.zeros();
        }
        
        // Reference mean and sd
        for (int i = 0; i < refmatnum; i++) {
          mat_mean = mat_mean + (as<arma::mat>(ref_byyear(i)) / static_cast<double>(refmatnum));
        }
        for (int i = 0; i < matlength; i++) {
          arma::vec mat_elems(refmatnum);
          mat_elems.zeros();
          
          for (int j = 0; j < refmatnum; j++) {
            mat_elems(j) = as<arma::mat>(ref_byyear(j))(i);
          }
          
          if (sum(mat_elems) != 0) {
            mat_sd(i) = stddev(mat_elems, 1);
          }
        }
        ref_matmean = mat_mean;
        ref_matsd = mat_sd;
        
        mat_mean.zeros();
        mat_sd.zeros();
      }
    }
    
    // Time weights
    if (tweights_.isNotNull()) {
      if (as<NumericVector>(tweights_).length() != numyears) {
        throw Rcpp::exception("Time weight vector must be the same length as the number of occasions in the reference matrix set.", false);
      }
      tweights = as<arma::vec>(tweights_);
    } else {
      tweights.resize(numyears);
      tweights.ones();
    } // At the end of this, we have an arma::vec holding with whatever the user supplied, or a 1 for each occasion
    
    // Here we set up the vector of chosen occasions, sampled from all possible occasions
    tweights = tweights / sum(tweights);
    theprophecy = Rcpp::RcppArmadillo::sample(uniqueyears, theclairvoyant, true, tweights);
    
    arma::vec startvec(matdim);
    startvec.ones();
    startvec = startvec / matdim; // This is the start vector for w and v calculations
    
    // Stochastic elasticities
    
    // The next section creates stable stage and rep value vectors arranged in
    // matrix format. The first two are general for whatever has been input,
    // whether historical or ahistorical, while the next two are specifically
    // for ahistorical versions of historical inputs
    arma::mat wprojection(startvec.n_elem, (theclairvoyant + 1));
    arma::mat vprojection(startvec.n_elem, (theclairvoyant + 1));
    wprojection.zeros();
    vprojection.zeros();
    
    arma::rowvec Rvecmat(theclairvoyant);
    Rvecmat.zeros();
    
    // Here we run the control loop to create the w and v values we need from the reference annual matrices
    arma::mat crazy_prophet = proj3(startvec, ref_byyear, theprophecy, 1, 0, 0);
    wprojection = crazy_prophet.submat(startvec.n_elem, 0, ((startvec.n_elem * 2) - 1), theclairvoyant);
    vprojection = crazy_prophet.submat((startvec.n_elem * 2), 0, ((startvec.n_elem * 3) - 1), theclairvoyant);
    Rvecmat = crazy_prophet.submat((startvec.n_elem * 3), 1, (startvec.n_elem * 3), theclairvoyant); // Rvec
    
    arma::mat sensmat(matdim, matdim);
    arma::mat elasmean1(matdim, matdim);
    arma::mat elassd1(matdim, matdim);
    sensmat.zeros();
    elasmean1.zeros();
    elassd1.zeros();
    
    // All references should go to elasmat
    for (int j = burnin; j < theclairvoyant; j++) { // This is the main time loop for the sensitivity matrices, 
      arma::vec vtplus1 = vprojection.col(j+1);
      arma::vec wtplus1 = wprojection.col(j+1);
      arma::vec wt = wprojection.col(j);
      
      arma::mat currentsens_num = vtplus1 * wt.as_row(); // This is the numerator of the key matrix equation
      arma::mat currentsens_den = (Rvecmat(j) * vtplus1.as_row() * wtplus1); // This is the denominator of the key matrix equation
      double cd_double = currentsens_den(0,0);
      sensmat = currentsens_num / (cd_double);
      
      elasmean1 = elasmean1 + ((sensmat % ref_matmean) / static_cast<double>(theclairvoyant - burnin));
      elassd1 = elassd1 + ((sensmat % ((as<arma::mat>(ref_byyear(theprophecy(j)))) - ref_matmean)) / static_cast<double>(theclairvoyant - burnin));
    }
    elas_mean = List::create(elasmean1);
    elas_sd = List::create(elassd1);
    
    // Now the difference matrix estimation
    arma::mat diffmean1(matdim, matdim);
    arma::mat diffsd1(matdim, matdim);
    diffmean1.zeros();
    diffsd1.zeros();
    
    // This creates indices of 0 elements for use in the next control loop
    diffmean1 = log(ref_matmean);
    diffsd1 = log(ref_matsd);
    diffmean1.elem(find_nonfinite(diffmean1)).zeros();
    diffsd1.elem(find_nonfinite(diffsd1)).zeros();
    
    ref_matmean = diffmean1;
    ref_matsd = diffsd1;
    
    for (int i = 0; i < numpoppatches; i++) {
      diffmean1.zeros();
      diffsd1.zeros();
      
      diffmean1 = log(as<arma::mat>(poppatch_meanmat(i)));
      diffsd1 = log(as<arma::mat>(poppatch_sdmat(i)));
      diffmean1.elem(find_nonfinite(diffmean1)).zeros();
      diffsd1.elem(find_nonfinite(diffsd1)).zeros();
      
      poppatch_meanmat(i) = diffmean1;
      poppatch_sdmat(i) = diffsd1;
      
      diffmean1.zeros();
      diffsd1.zeros();
      
      diffmean1 = as<arma::mat>(poppatch_meanmat(i)) - ref_matmean;
      diffsd1 = as<arma::mat>(poppatch_sdmat(i)) - ref_matsd;
      
      // In Davison's original code, all elements equal to 0 in the reference
      // matrices must also be 0s in the difference matrices
      
      // Now the contributions
      diffmean1 = diffmean1 % as<arma::mat>(elas_mean(0));
      diffsd1 = diffsd1 % as<arma::mat>(elas_sd(0));
      
      if (i == 0) {
        cont_meanmat = List::create(diffmean1);
        cont_sdmat = List::create(diffsd1);
      } else {
        cont_meanmat.push_back(diffmean1);
        cont_sdmat.push_back(diffsd1);
      }
    }
  } else {
    // Sparse matrix analysis
    arma::sp_mat mat_mean(matdim, matdim);
    arma::sp_mat mat_sd(matdim, matdim);
    arma::sp_mat ref_matmean(matdim, matdim);
    arma::sp_mat ref_matsd(matdim, matdim);
    mat_mean.zeros();
    mat_sd.zeros();
    ref_matmean.zeros();
    ref_matsd.zeros();
  
    // First the pop/patch means and sds
    for (int i = 0; i < numpoppatches; i++) {
      arma::uvec poppatch_chosen = find(poppatchc == uniquepoppatches(i));
      int numpoppatch_chosen = poppatch_chosen.n_elem;
      
      for (int j = 0; j < numpoppatch_chosen; j++) {
        mat_mean = mat_mean + (arma::sp_mat(as<arma::mat>(Amats(poppatch_chosen(j)))) / static_cast<double>(numpoppatch_chosen));
      }
      if (i == 0) {
        poppatch_meanmat = List::create(mat_mean);
      } else {
        poppatch_meanmat.push_back(mat_mean);
      }
      mat_mean.zeros();
    }
    
    for (int i = 0; i < numpoppatches; i++) {
      arma::uvec poppatch_chosen = find(poppatchc == uniquepoppatches(i));
      int numpoppatch_chosen = poppatch_chosen.n_elem;
      arma::vec mat_elems(numpoppatch_chosen);
      
      for (int j = 0; j < matlength; j++) {
        mat_elems.zeros();
        
        for (int k = 0; k < numpoppatch_chosen; k++) {
          mat_elems(k) = as<arma::mat>(Amats(poppatch_chosen(k)))(j);
        }
        
        if (sum(mat_elems) != 0) {
          mat_sd(j) = stddev(mat_elems, 1);
        }
      }
      if (i == 0) {
        poppatch_sdmat = List::create(mat_sd);
      } else {
        poppatch_sdmat.push_back(mat_sd);
      }
      mat_sd.zeros();
    }
    
    // Here we create the reference matrix sets
    if (refmats_.isNotNull()) {
      Rcpp::List refmats(refmats_);
      
      mat_mean.zeros();
      mat_sd.zeros();
      
      ref_byyear = List::create(refmats);
      
      for (int i = 0; i < refmatnum; i++) {
        mat_mean = mat_mean + (arma::sp_mat(as<arma::mat>(refmats(refnum(i)))) / static_cast<double>(refmatnum));
      }
      
      for (int i = 0; i < matlength; i++) {
        arma::vec mat_elems(refmatnum);
        mat_elems.zeros();
        
        for (int j = 0; j < refmatnum; j++) {
          mat_elems(j) = arma::sp_mat(as<arma::mat>(refmats(refnum(j))))(i);
        }
        
        if (sum(mat_elems) != 0) {
          mat_sd(i) = stddev(mat_elems, 1);
        }
      }
      
      ref_matmean = mat_mean;
      ref_matsd = mat_sd;
      
    } else {
      
      if (refmatnum == Amatnum) {
        
        // Reference by year
        for (int i = 0; i < numyears; i++) {
          arma::uvec year_chosen = find(year2c == uniqueyears(i));
          int numyear_chosen = year_chosen.n_elem;
          
        for (int j = 0; j < numyear_chosen; j++) {
          mat_mean = mat_mean + (arma::sp_mat(as<arma::mat>(Amats(year_chosen(j)))) / static_cast<double>(numyear_chosen));
        }
        if (i == 0) {
          ref_byyear = List::create(mat_mean);
        } else {
          ref_byyear.push_back(mat_mean);
        }
        mat_mean.zeros();
        }
        
        // Reference mean and sd
        for (int i = 0; i < numyears; i++) {
          mat_mean = mat_mean + (as<arma::sp_mat>(ref_byyear(i)) / static_cast<double>(numyears));
        }
        for (int i = 0; i < matlength; i++) {
          arma::vec mat_elems(refmatnum);
          mat_elems.zeros();
          
          for (int j = 0; j < numyears; j++) {
            mat_elems(j) = as<arma::sp_mat>(ref_byyear(j))(i);
          }
          
          if (sum(mat_elems) != 0) {
            mat_sd(i) = stddev(mat_elems, 1);
          }
        }
        ref_matmean = mat_mean;
        ref_matsd = mat_sd;
        
        mat_mean.zeros();
        mat_sd.zeros();
      } else {
        
        // Reference by year
        for (int i = 0; i < refmatnum; i++) {
          if (i == 0) {
            ref_byyear = List::create(arma::sp_mat(as<arma::mat>(Amats(i))));
          } else {
            ref_byyear.push_back(arma::sp_mat(as<arma::mat>(Amats(i))));
          }
          mat_mean.zeros();
        }
        
        // Reference mean and sd
        for (int i = 0; i < refmatnum; i++) {
          mat_mean = mat_mean + (as<arma::sp_mat>(ref_byyear(i)) / static_cast<double>(refmatnum));
        }
        for (int i = 0; i < matlength; i++) {
          arma::vec mat_elems(refmatnum);
          mat_elems.zeros();
          
          for (int j = 0; j < refmatnum; j++) {
            mat_elems(j) = as<arma::sp_mat>(ref_byyear(j))(i);
          }
          
          if (sum(mat_elems) != 0) {
            mat_sd(i) = stddev(mat_elems, 1);
          }
        }
        ref_matmean = mat_mean;
        ref_matsd = mat_sd;
        
        mat_mean.zeros();
        mat_sd.zeros();
      }
    }
    
    // Time weights
    if (tweights_.isNotNull()) {
      if (as<NumericVector>(tweights_).length() != numyears) {
        throw Rcpp::exception("Time weight vector must be the same length as the number of occasions in the reference matrix set.", false);
      }
      tweights = as<arma::vec>(tweights_);
    } else {
      tweights.resize(numyears);
      tweights.ones();
    } // At the end of this, we have an arma::vec holding with whatever the user supplied, or a 1 for each occasion
    
    // Here we set up the vector of chosen occasions, sampled from all possible occasions
    tweights = tweights / sum(tweights);
    theprophecy = Rcpp::RcppArmadillo::sample(uniqueyears, theclairvoyant, true, tweights);
    
    arma::vec startvec(matdim);
    startvec.ones();
    startvec = startvec / matdim; // This is the start vector for w and v calculations
    
    // Stochastic elasticities
    
    // The next section creates stable stage and rep value vectors arranged in
    // matrix format. The first two are general for whatever has been input,
    // whether historical or ahistorical, while the next two are specifically
    // for ahistorical versions of historical inputs
    arma::mat wprojection(startvec.n_elem, (theclairvoyant + 1));
    arma::mat vprojection(startvec.n_elem, (theclairvoyant + 1));
    wprojection.zeros();
    vprojection.zeros();
    
    arma::rowvec Rvecmat(theclairvoyant);
    Rvecmat.zeros();
    
    // Here we run the control loop to create the w and v values we need from the reference annual matrices
    arma::mat crazy_prophet = proj3sp(startvec, ref_byyear, theprophecy, 1, 0, 0);
    wprojection = crazy_prophet.submat(startvec.n_elem, 0, ((startvec.n_elem * 2) - 1), theclairvoyant);
    vprojection = crazy_prophet.submat((startvec.n_elem * 2), 0, ((startvec.n_elem * 3) - 1), theclairvoyant);
    Rvecmat = crazy_prophet.submat((startvec.n_elem * 3), 1, (startvec.n_elem * 3), theclairvoyant); // Rvec
    
    arma::sp_mat sensmat(matdim, matdim);
    arma::sp_mat elasmean1(matdim, matdim);
    arma::sp_mat elassd1(matdim, matdim);
    sensmat.zeros();
    elasmean1.zeros();
    elassd1.zeros();
    
    // All references should go to elasmat
    for (int j = burnin; j < theclairvoyant; j++) { // This is the main time loop for the sensitivity matrices, 
      arma::vec vtplus1 = vprojection.col(j+1);
      arma::vec wtplus1 = wprojection.col(j+1);
      arma::vec wt = wprojection.col(j);
      
      arma::mat currentsens_num = vtplus1 * wt.as_row(); // This is the numerator of the key matrix equation
      arma::mat currentsens_den = (Rvecmat(j) * vtplus1.as_row() * wtplus1); // This is the denominator of the key matrix equation
      double cd_double = currentsens_den(0,0);
      sensmat = arma::sp_mat(currentsens_num / (cd_double));
      
      elasmean1 = elasmean1 + ((sensmat % ref_matmean) / (static_cast<double>(theclairvoyant - burnin)));
      elassd1 = elassd1 + ((sensmat % ((as<arma::sp_mat>(ref_byyear(theprophecy(j)))) - ref_matmean)) / (static_cast<double>(theclairvoyant - burnin)));
    }
    elas_mean = List::create(elasmean1);
    elas_sd = List::create(elassd1);
    
    // Now the difference and contribution matrix estimation
    ref_matmean = spmat_log(ref_matmean);
    ref_matsd = spmat_log(ref_matsd);
    
    arma::sp_mat diffmean1(matdim, matdim);
    arma::sp_mat diffsd1(matdim, matdim);
    
    for (int i = 0; i < numpoppatches; i++) {
      diffmean1.zeros();
      diffsd1.zeros();
      
      poppatch_meanmat(i) = spmat_log(as<arma::sp_mat>(poppatch_meanmat(i)));
      poppatch_sdmat(i) = spmat_log(as<arma::sp_mat>(poppatch_sdmat(i)));
      
      diffmean1 = as<arma::sp_mat>(poppatch_meanmat(i)) - ref_matmean;
      diffsd1 = as<arma::sp_mat>(poppatch_sdmat(i)) - ref_matsd;
      
      // Now the contributions
      diffmean1 = diffmean1 % as<arma::sp_mat>(elas_mean(0));
      diffsd1 = diffsd1 % as<arma::sp_mat>(elas_sd(0));
      
      if (i == 0) {
        cont_meanmat = List::create(arma::mat(diffmean1));
        cont_sdmat = List::create(arma::mat(diffsd1));
      } else {
        cont_meanmat.push_back(arma::mat(diffmean1));
        cont_sdmat.push_back(arma::mat(diffsd1));
      }
    }
  }
  return List::create(Named("cont_mean") = cont_meanmat,
    _["cont_sd"] = cont_sdmat);
}

