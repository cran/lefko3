#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Core Engine for cond_hmpm()
//' 
//' Creates a list of conditional ahistorical matrices in the style noted in
//' deVries and Caswell (2018).
//'
//' @param mainmat Historical matrix.
//' @param indices Data frame including the stages at times t-1, t, and t+1, as
//' well as indices corresponding to elements in the main historical matrix and
//' the conditional matrices to be produced.
//' @param ahstages The number of stages in the stageframe.
//' @param stageframe The original stageframe for the input matrices.
//'
//' @return A list of ahistorical matrices.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List hoffmannofstuttgart(arma::mat mainmat, DataFrame indices, int ahstages,
  StringVector stagenames) {
  arma::uvec stage1 = indices["stage1"];
  arma::uvec stage2 = indices["stage2"];
  arma::uvec stage3 = indices["stage3"];
  arma::ivec main_index = indices["main_index"];
  arma::ivec new_index = indices["new_index"];
  
  arma::mat newmatrix(ahstages, ahstages);
  newmatrix.zeros();
  
  arma::uvec condidx = find(stage1 == 1);
  int condlength = condidx.n_elem;
  
  for (int i = 0; i < condlength; i++) {
    if (mainmat(main_index(condidx(i))) > 0) newmatrix(new_index(condidx(i))) = 
      newmatrix(new_index(condidx(i))) + mainmat(main_index(condidx(i)));
  }
  
  Rcpp::List condlist = List::create(Named("1") = newmatrix);
  
  for (int i = 1; i < ahstages; i++) {
    newmatrix.zeros();
    
    condidx = find(stage1 == (i+1));
    condlength = condidx.n_elem;
    
    for (int j = 0; j < condlength; j++) {
      if (mainmat(main_index(condidx(j))) > 0) newmatrix(new_index(condidx(j))) = 
        newmatrix(new_index(condidx(j))) + mainmat(main_index(condidx(j)));
    }
    
    condlist.push_back(newmatrix);
  }
  condlist.names() = stagenames;
  
  return (condlist);
}

//' Extract Conditional Ahistorical Matrices from Historical MPM
//' 
//' Function \code{cond_hmpm()} takes historical MPMs and decomposes them into 
//' ahistorical matrices conditional upon stage in time \emph{t}-1. In effect,
//' the function takes each historical matrix within a lefkoMat object, and
//' forms one ahistorical matrix for each stage in time \emph{t}-1.
//' 
//' @param hmpm A historical matrix projection model of class \code{lefkoMat}.
//' @param matchoice A character denoting whether to use A, U, or F matrices.
//' Defaults to A matrices.
//' @param err_check A logical value denoting whether to include a data frame
//' of element equivalence from the conditional matrices to the original
//' matrices. Used only for debugging purposes. Defaults to FALSE.
//' 
//' @return A \code{lefkoCondMat} object, with the following elements:
//' 
//' \item{Acond}{A multi-level list holding the conditional A matrices derived
//' from the input \code{lefkoMat} object. The top level of the list corresponds
//' to each historical matrix in turn, and the lower level corresponds to each
//' stage in time \emph{t}-1, with individual conditional matrices named for the
//' latter.}
//' \item{hstages}{A data frame matrix showing the pairing of ahistorical stages
//' used to create historical stage pairs.}
//' \item{ahstages}{A data frame detailing the characteristics of associated
//' ahistorical stages.}
//' \item{labels}{A data frame showing the patch and year of each input full A 
//' matrix in order.}
//' 
//' @examples
//' data(cypdata)
//'  
//' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
//' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md",
//'   "Lg", "XLg")
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
//' cypsupp3r <- supplemental(stage3 = c("SD", "SD", "P1", "P1", "P2", "P3",
//'     "SL", "SL", "SL", "D", "XSm", "Sm", "D", "XSm", "Sm", "SD", "P1"),
//'   stage2 = c("SD", "SD", "SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL",
//'     "SL", "SL", "SL", "SL", "SL", "rep", "rep"),
//'   stage1 = c("SD", "rep", "SD", "rep", "SD", "P1", "P2", "P3", "SL", "P3",
//'     "P3", "P3", "SL", "SL", "SL", "all", "all"),
//'   eststage3 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", "D",
//'     "XSm", "Sm", NA, NA),
//'   eststage2 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm",
//'     "XSm", "XSm", "XSm", NA, NA),
//'   eststage1 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm",
//'     "XSm", "XSm", "XSm", NA, NA),
//'   givenrate = c(0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.25, 0.4, 0.4, NA, NA, NA, NA,
//'     NA, NA, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
//'     0.5, 0.5),
//'   type = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
//'   type_t12 = c(1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
//'   stageframe = cypframe_raw, historical = TRUE)
//' 
//' cypmatrix3r <- rlefko3(data = cypraw_v1, stageframe = cypframe_raw,
//'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
//'   size = c("size3added", "size2added", "size1added"),
//'   supplement = cypsupp3r, yearcol = "year2", patchcol = "patchid",
//'   indivcol = "individ")
//' 
//' cypcondmats <- cond_hmpm(cypmatrix3r)
//' summary(cypcondmats)
//' 
//' @export cond_hmpm
// [[Rcpp::export]]
Rcpp::List cond_hmpm(List hmpm, Nullable<CharacterVector> matchoice = R_NilValue,
  Nullable<LogicalVector> err_check = R_NilValue) {
  
  CharacterVector usedmats("a");
  int iusedmats {1};
  LogicalVector err_v = {false};
  bool err_c = false;
  
  if (matchoice.isNotNull()) {
    usedmats = matchoice;
    
    if (usedmats(0) == "A" || usedmats(0) == "a") {
      iusedmats = 1;
    } else if (usedmats(0) == "U" || usedmats(0) == "u") {
      iusedmats = 2;
    } else if (usedmats(0) == "F" || usedmats(0) == "f") {
      iusedmats = 3;
    } else {
      Rcpp::Rcout << "Choice of matrix not recognized. Using A matrices.";
      iusedmats = 1;
    }
  }
  
  if (err_check.isNotNull()) {
    err_v = err_check;
    if (err_v(0)) {
      err_c = true;
    }
  }
  
  List amats = hmpm["A"];
  List umats = hmpm["U"];
  List fmats = hmpm["F"];
  List hstages = hmpm["hstages"];
  List stageframe = hmpm["ahstages"];
  List labels = hmpm["labels"];
  int numofmats = amats.length();
  
  arma::uvec hstage1 = hstages["stage_id_1"];
  arma::uvec hstage2 = hstages["stage_id_2"];
  arma::uvec ahstages = stageframe["stage_id"];
  StringVector stagenames = stageframe["stage"];
  StringVector matnames(numofmats);
  
  for (int i = 0; i < numofmats; i++) {
    matnames(i) = i + 1;
  }
  
  int ahmpm_rows = ahstages.n_elem;
  int hmpm_rows = hstage1.n_elem;
  int hmpm_elems = 2 * ahmpm_rows * ahmpm_rows * ahmpm_rows;
  
  int format_int {0};
  if (stagenames(ahmpm_rows - 1) == "AlmostBorn") format_int = 1;
  
  arma::uvec stage1(hmpm_elems);
  arma::uvec stage2(hmpm_elems);
  arma::uvec stage3(hmpm_elems);
  arma::ivec main_index(hmpm_elems);
  arma::ivec new_index(hmpm_elems);
  stage1.zeros();
  stage2.zeros();
  stage3.zeros();
  main_index.fill(-1);
  new_index.fill(-1);
  
  int stage1_proxy {0};
  int stage2o_proxy {0};
  int stage2n_proxy {0};
  int stage3_proxy {0};
  
  int counter = 0;
  
  // This next section creates an index of matrix elements in the original matrices
  // and the new conditional matrices. It must be run on the hstages inputs in
  // order to handle reduced matrices properly.
  if (format_int == 0) {
    // This portion handles Ehrlen-format hMPMs
    for (int prior = 0; prior < hmpm_rows; prior++) {
      for (int post = 0; post < hmpm_rows; post++) {
        stage1_proxy = hstage1[prior];
        stage2o_proxy = hstage2[prior];
        stage2n_proxy = hstage1[post];
        stage3_proxy = hstage2[post];
        
        if (stage2o_proxy == stage2n_proxy) {
          stage1[counter] = stage1_proxy;
          stage2[counter] = stage2n_proxy;
          stage3[counter] = stage3_proxy;
              
          main_index[counter] = post + (prior * hmpm_rows);
          new_index[counter] = (stage3[counter] - 1) + ((stage2[counter] - 1) * ahmpm_rows);
          
          counter++;
        }
      }
    }
  } else {
    // This section deals with deVries-format hMPMs
    for (int prior = 0; prior < hmpm_rows; prior++) {
      for (int post = 0; post < hmpm_rows; post++) {
        stage1_proxy = hstage1[prior];
        stage2o_proxy = hstage2[prior];
        stage2n_proxy = hstage1[post];
        stage3_proxy = hstage2[post];
        
        if (stage2o_proxy == stage2n_proxy) {
          stage1[counter] = stage1_proxy;
          stage2[counter] = stage2o_proxy;
          stage3[counter] = stage3_proxy;
              
          main_index[counter] = post + (prior * hmpm_rows);
          new_index[counter] = (stage3[counter] - 1) + ((stage2[counter] - 1) * ahmpm_rows);
          
          counter++;
        } else if (stage2n_proxy == ahmpm_rows && stage1_proxy < ahmpm_rows) {
          stage1[counter] = stage1_proxy;
          stage2[counter] = stage2o_proxy;
          stage3[counter] = stage3_proxy;
              
          main_index[counter] = post + (prior * hmpm_rows);
          new_index[counter] = (stage3[counter] - 1) + ((stage2[counter] - 1) * ahmpm_rows);
              
          counter++;
        }
      }
    }
  }
  
  arma::uvec elements_to_use = find(stage1 > 0);
  stage1 = stage1.elem(elements_to_use);
  stage2 = stage2.elem(elements_to_use);
  stage3 = stage3.elem(elements_to_use);
  main_index = main_index.elem(elements_to_use);
  new_index = new_index.elem(elements_to_use);
  
  DataFrame loveontherocks = DataFrame::create(Named("stage1") = stage1,
    _["stage2"] = stage2, _["stage3"] = stage3, _["main_index"] = main_index,
    _["new_index"] = new_index);
  
  // This next line creates the conditional matrices in list form
  List mats1 = hoffmannofstuttgart(amats(0), loveontherocks, ahmpm_rows, stagenames);
  
  if (iusedmats == 3) {
    mats1 = hoffmannofstuttgart(fmats(0), loveontherocks, ahmpm_rows, stagenames);
  } else if (iusedmats == 2) {
    mats1 = hoffmannofstuttgart(umats(0), loveontherocks, ahmpm_rows, stagenames);
  } 
  
  List allout = List::create(Named("1") = mats1);
  
  if (iusedmats == 3) {
    for (int i = 1; i < numofmats; i++) {
      mats1 = hoffmannofstuttgart(fmats(i), loveontherocks, ahmpm_rows, stagenames);
      
      allout.push_back(mats1);
    }
  } else if (iusedmats == 2) {
    for (int i = 1; i < numofmats; i++) {
      mats1 = hoffmannofstuttgart(umats(i), loveontherocks, ahmpm_rows, stagenames);
      
      allout.push_back(mats1);
    }
  } else {
    for (int i = 1; i < numofmats; i++) {
      mats1 = hoffmannofstuttgart(amats(i), loveontherocks, ahmpm_rows, stagenames);
      
      allout.push_back(mats1);
    }
  }
  allout.names() = matnames;
  
  List panama = List::create(Named("Acond") = allout, _["hstages"] = hstages,
    _["ahstages"] = stageframe, _["labels"] = labels);
  
  if (err_c) {
    Rcpp::List morestuff = List::create(_["err_check"] = loveontherocks);
    panama.push_back(morestuff);
  }
  panama.attr("class") = "lefkoCondMat";
  
  return panama;
}
