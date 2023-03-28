#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <LefkoUtils.h>

using namespace Rcpp;
using namespace arma;
using namespace LefkoUtils;
using namespace LefkoMats;


//' Creates Size Index for Elasticity Summaries of hMPMs
//' 
//' Function \code{bambi3()} creates an index of estimable elements in
//' historical matrices, and details the kind of transition that it is.
//' 
//' @name .bambi3
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
//' \item{repstatus1}{Reproductive status in occasion \emph{t}-1.}
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
DataFrame bambi3(const DataFrame& stages, const DataFrame& hstages) {
  
  StringVector stagenames = as<StringVector>(stages["stage"]);
  arma::uvec astages = as<arma::uvec>(stages["stage_id"]);
  arma::vec sizes = as<arma::vec>(stages["original_size"]);
  arma::uvec repstatus = as<arma::uvec>(stages["repstatus"]);
  arma::uvec entrystage = as<arma::uvec>(stages["entrystage"]);
  int numstages = static_cast<int>(astages.n_elem);
  
  arma::uvec hstage3in = as<arma::uvec>(hstages["stage_id_2"]);
  arma::uvec hstage2nin = as<arma::uvec>(hstages["stage_id_1"]);
  int numhstages = static_cast<int>(hstage3in.n_elem);
  
  hstage3in = hstage3in - 1;
  hstage2nin = hstage2nin - 1;
  int predictedsize = numstages * numstages * numstages;
  
  arma::ivec hsindexl(predictedsize);
  arma::uvec transition_type(predictedsize);
  hsindexl.fill(-1);
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
      if (hstage3in(i1) == (hstage2nin(i2))) {
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
            if (repstatus2(counter) > repstatus1(counter) ||
                entrystatus3(counter) < entrystatus2(counter)) {
              transition_type(counter) = 21; // Growth to fecundity
            } else if (repstatus2(counter) < repstatus1(counter) ||
              entrystatus2(counter) > entrystatus1(counter)) {
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
            if (repstatus3(counter) > repstatus2(counter) ||
              entrystatus3(counter) < entrystatus2(counter)) {
              transition_type(counter) = 24; // Fecundity to growth
            } else if (repstatus3(counter) < repstatus2(counter) ||
              entrystatus3(counter) > entrystatus2(counter)) {
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
          if (repstatus2(counter) > repstatus1(counter) ||
            entrystatus2(counter) < entrystatus1(counter)) {
            if (repstatus3(counter) > repstatus2(counter) ||
              entrystatus3(counter) < entrystatus2(counter)) {
              transition_type(counter) = 12; // Full growth
            } else if (repstatus3(counter) < repstatus2(counter) ||
              entrystatus3(counter) > entrystatus2(counter)) {
              transition_type(counter) = 17; // Growth to shrinkage
            } else {
              transition_type(counter) = 13; // Growth to stasis
            }
          } else if (repstatus2(counter) < repstatus1(counter)) {
            if (repstatus3(counter) > repstatus2(counter) ||
              entrystatus3(counter) < entrystatus2(counter)) {
              transition_type(counter) = 18; // Shrinkage to growth
            } else if (repstatus3(counter) < repstatus2(counter) ||
              entrystatus3(counter) > entrystatus2(counter)) {
              transition_type(counter) = 15; // Full shrinkage
            } else {
              transition_type(counter) = 16; // Shrinkage to stasis
            }
          } else {
            if (repstatus3(counter) > repstatus2(counter) ||
              entrystatus3(counter) < entrystatus2(counter)) {
              transition_type(counter) = 11; // Stasis to growth
            } else if (repstatus3(counter) < repstatus2(counter) ||
              entrystatus3(counter) > entrystatus2(counter)) {
              transition_type(counter) = 14; // Stasis to shrinkage
            } else {
              transition_type(counter) = 10; // Full stasis
            }
          }
        } else if (size3(counter) > size2(counter) && size2(counter) == size1(counter)) {
          if (repstatus2(counter) > repstatus1(counter) ||
            entrystatus2(counter) < entrystatus1(counter)) {
            transition_type(counter) = 12; // Full growth
          } else if (repstatus2(counter) < repstatus1(counter) ||
            entrystatus2(counter) > entrystatus1(counter)) {
            transition_type(counter) = 18; // Shrinkage to growth
          } else {
            transition_type(counter) = 11; // Stasis to growth
          }
        } else if (size3(counter) > size2(counter) && size2(counter) > size1(counter)) {
          transition_type(counter) = 12; // Full growth
        } else if (size3(counter) == size2(counter) && size2(counter) > size1(counter)) {
          if (repstatus3(counter) > repstatus2(counter) ||
            entrystatus3(counter) < entrystatus2(counter)) {
            transition_type(counter) = 12; // Full growth
          } else if (repstatus3(counter) < repstatus2(counter) ||
            entrystatus3(counter) > entrystatus2(counter)) {
            transition_type(counter) = 17; // Growth to shrinkage
          } else {
            transition_type(counter) = 13; // Growth to stasis
          }
        } else if (size3(counter) < size2(counter) && size2(counter) == size1(counter)) {
          if (repstatus2(counter) > repstatus1(counter) ||
            entrystatus2(counter) < entrystatus1(counter)) {
            transition_type(counter) = 17; // Growth to shrinkage
          } else if (repstatus2(counter) < repstatus1(counter) ||
            entrystatus2(counter) > entrystatus1(counter)) {
            transition_type(counter) = 15; // Full shrinkage
          } else {
            transition_type(counter) = 14; // Stasis to shrinkage
          }
        } else if (size3(counter) < size2(counter) && size2(counter) < size1(counter)) {
          transition_type(counter) = 15; // Full shrinkage
        } else if (size3(counter) == size2(counter) && size2(counter) < size1(counter)) {
          if (repstatus3(counter) > repstatus2(counter) ||
            entrystatus3(counter) < entrystatus2(counter)) {
            transition_type(counter) = 18; // Shrinkage to growth
          } else if (repstatus3(counter) < repstatus2(counter) ||
            entrystatus3(counter) > entrystatus2(counter)) {
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
  
  DataFrame output = DataFrame::create(Named("index") = hsindex,
    _["transition"] = t_type, _["stage3"] = names3, _["size3"] = size3c,
    _["repstatus3"] = r_status3, _["entrystatus3"] = e_status3,
    _["stage2"] = names2, _["size2"] = size2c, _["repstatus2"] = r_status2,
    _["entrystatus2"] = e_status2, _["stage1"] = names1, _["size1"] = size1c,
    _["repstatus1"] = r_status1, _["entrystatus1"] = e_status1);
  
  return output;
}

//' Creates Size Index for Elasticity Summaries of ahMPMs
//' 
//' Function \code{bambi2()} creates an index of estimable elements in
//' ahistorical matrices, and details the kind of transition that it is.
//' 
//' @name .bambi2
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
DataFrame bambi2(const DataFrame& stages) {
  
  StringVector stagenames = as<StringVector>(stages["stage"]);
  arma::uvec astages = as<arma::uvec>(stages["stage_id"]);
  arma::vec sizes = as<arma::vec>(stages["original_size"]);
  arma::uvec repstatus = as<arma::uvec>(stages["repstatus"]);
  arma::uvec entrystage = as<arma::uvec>(stages["entrystage"]);
  int numstages = static_cast<int>(astages.n_elem);
  astages = astages - 1;
  int predictedsize = numstages * numstages;
  
  arma::ivec ahsindexl(predictedsize);
  arma::uvec transition_type(predictedsize);
  ahsindexl.fill(-1);
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
        if (repstatus3(counter) > repstatus2(counter) ||
          entrystatus3(counter) < entrystatus2(counter)) {
          transition_type(counter) = 2; // Growth
        } else if (repstatus3(counter) < repstatus2(counter) ||
          entrystatus3(counter) > entrystatus2(counter)) {
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
  
  DataFrame output = DataFrame::create(Named("index") = ahsindex,
    _["transition"] = transition_type, _["stage3"] = longstages3,
    _["size3"] = size3, _["repstatus3"] = repstatus3,
    _["entrystatus3"] = entrystatus3, _["stage2"] = longstages2,
    _["size2"] = size2, _["repstatus2"] = repstatus2,
    _["entrystatus2"] = entrystatus2);
  
  return output;
}

//' Sum Positive and Negative LTRE or Elasticity Contributions
//' 
//' @name .demolition4
//' 
//' Function \code{demolition4()} takes \code{lefkoElas} and \code{lefkoLTRE}
//' inputs and returns a data frame summarizing their positive and negative
//' contributions.
//' 
//' @param cmats Any \code{lefkoElas} or \code{lefkoLTRE} object.
//' 
//' @return A data frame with three columns per contribution type, corresponding
//' to summed positive, negative, and total contributions, and rows
//' corresponding to each input matrix.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.demolition4)]]
Rcpp::List demolition4 (List cmats) {
  bool ltre_input {false};
  bool sparse_input {false};
  
  StringVector cmats_class = as<StringVector>(cmats.attr("class"));
  for (int i = 0; i < cmats_class.length(); i++) {
    if (stringcompare_hard(as<std::string>(cmats_class(i)), "lefkoLTRE")) {
      ltre_input = true;
    } else if (stringcompare_hard(as<std::string>(cmats_class(i)), "lefkoElas")) {
      ltre_input = false;
    } else {
      throw Rcpp::exception("Unrecognized input in summary.lefkoLTRE().", false);
    }
  }
  
  List output_df;
  StringVector varnames_;
  IntegerVector rowlabels;
  
  if (ltre_input) { 
    List cont_0 = as<List>(cmats(0));
    int conts_num = cmats.length() - 4;
    if (conts_num > 4) conts_num -= 2;
    int mats_num = cont_0.length();
    
    if (is<S4>(cont_0(0))) sparse_input = true;
    
    List temp_df ((conts_num * 4) + 1);
    rowlabels = seq(1, mats_num);
    temp_df(0) = rowlabels;
    StringVector varnames ((conts_num * 4) + 1);
    varnames(0) = "matrix";
    
    for (int i = 0; i < conts_num; i++) {
      List current_cont = as<List>(cmats(i));
      
      NumericVector pos_cont (mats_num);
      NumericVector neg_cont (mats_num);
      NumericVector abs_cont (mats_num);
      NumericVector tot_cont (mats_num);
      
      for (int j = 0; j < mats_num; j++) {
        if (!sparse_input) { 
          arma::mat current_mat = as<arma::mat>(current_cont(j));
          arma::uvec pos_vec = find(current_mat > 0.0);
          arma::uvec neg_vec = find(current_mat < 0.0);
          
          double pos_sum = accu(current_mat.elem(pos_vec));
          double neg_sum = accu(current_mat.elem(neg_vec));
          
          pos_cont(j) = pos_sum;
          neg_cont(j) = neg_sum;
          abs_cont(j) = pos_sum + (-1.0 * neg_sum);
          tot_cont(j) = pos_sum + neg_sum;
          
        } else {
          arma::sp_mat current_mat = as<arma::sp_mat>(current_cont(j));
          
          double pos_sum {0.0};
          double neg_sum {0.0};
          
          for (int k = 0; k < static_cast<int>(current_mat.n_elem); k++) {
            if (current_mat(k) > 0.0) {
              pos_sum += static_cast<double>(current_mat(k));
            } else if (current_mat(k) < 0.0) {
              neg_sum += static_cast<double>(current_mat(k));
            }
          }
          
          pos_cont(j) = pos_sum;
          neg_cont(j) = neg_sum;
          abs_cont(j) = pos_sum + (-1.0 * neg_sum);
          tot_cont(j) = pos_sum + neg_sum;
          
        }
      }
      
      if (i == 0) {
        varnames(1) = "means_positive";
        varnames(2) = "means_negative";
        varnames(3) = "means_abs_sum";
        varnames(4) = "means_total";
      }
      
      if (i == 1) {
        if (conts_num == 2) {
          varnames(5) = "sd_positive";
          varnames(6) = "sd_negative";
          varnames(7) = "sd_abs_sum";
          varnames(8) = "sd_total";
        } else {
          varnames(5) = "elas_positive";
          varnames(6) = "elas_negative";
          varnames(7) = "elas_abs_sum";
          varnames(8) = "elas_total";
        }
      }
      
      if (i == 2) {
        varnames(9) = "cv_positive";
        varnames(10) = "cv_negative";
        varnames(11) = "cv_abs_sum";
        varnames(12) = "cv_total";
      }
      
      if (i == 3) {
        varnames(13) = "cor_positive";
        varnames(14) = "cor_negative";
        varnames(15) = "cor_abs_sum";
        varnames(16) = "cor_total";
      }
      
      temp_df(1 + (4 * i) + 0) = pos_cont;
      temp_df(1 + (4 * i) + 1) = neg_cont;
      temp_df(1 + (4 * i) + 2) = abs_cont;
      temp_df(1 + (4 * i) + 3) = tot_cont;
    }
    output_df = temp_df;
    varnames_ = varnames;
  }
  
  StringVector rowlabels_str = as<StringVector>(rowlabels);
  StringVector out_class = {"data.frame"};
  output_df.attr("class") = out_class;
  output_df.attr("names") = varnames_;
  output_df.attr("row.names") = rowlabels_str;
  
  return output_df;
}

//' Creates Summary Data for Elasticity Matrix Inputs
//' 
//' Function \code{demolition3()} sums elasticity values from elasticity
//' matrices, and LTRE contributions from LTRE and sLTRE matrices, according to
//' the categories developed by functions \code{bambi2()} and \code{bambi3()}.
//' It requires \code{matrix} class inputs.
//' 
//' @name .demolition3
//' 
//' @param e_amat A single elasticity, LTRE, or sLTRE matrix of class
//' \code{matrix}.
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
List demolition3(const arma::mat& e_amat, const DataFrame& bambesque,
  Nullable<Rcpp::NumericMatrix> amat_ = R_NilValue,
  Nullable<Rcpp::NumericMatrix> fmat_ = R_NilValue) {
  
  arma::uvec eindices = as<arma::uvec>(bambesque["index"]);
  arma::uvec categories = as<arma::uvec>(bambesque["transition"]);
  
  int e_amatsize = static_cast<int>(e_amat.n_elem);
  int e_amatrows = static_cast<int>(e_amat.n_rows);
  int maxelem = static_cast<int>(eindices.max());
  int minindex = static_cast<int>(categories.min());
  
  arma::mat amat;
  arma::mat fmat;
  
  if (maxelem > e_amatsize) {
    throw Rcpp::exception("Supplied info does not correspond to input matrices.",
      false);
  }
  
  if (amat_.isNotNull() && fmat_.isNotNull()) {
    amat = Rcpp::as<arma::mat>(amat_);
    fmat = Rcpp::as<arma::mat>(fmat_);
    
  } else {
    arma::mat amat1(e_amatrows, e_amatrows, fill::ones);
    arma::mat fmat1(e_amatrows, e_amatrows, fill::zeros);
    amat = amat1;
    fmat = fmat1;
    
    arma::uvec fec_trans = find(categories == 4);
    if (fec_trans.n_elem > 0) {
      for (int i = 0; i < static_cast<int>(fec_trans.n_elem); i ++) {
        fmat(eindices(fec_trans(i))) = 1;
      }
    }
    fec_trans = find(categories == 20);
    if (fec_trans.n_elem > 0) {
      for (int i = 0; i < static_cast<int>(fec_trans.n_elem); i ++) {
        fmat(eindices(fec_trans(i))) = 1;
      }
    }
    fec_trans = find(categories == 21);
    if (fec_trans.n_elem > 0) {
      for (int i = 0; i < static_cast<int>(fec_trans.n_elem); i ++) {
        fmat(eindices(fec_trans(i))) = 1;
      }
    }
    fec_trans = find(categories == 22);
    if (fec_trans.n_elem > 0) {
      for (int i = 0; i < static_cast<int>(fec_trans.n_elem); i ++) {
        fmat(eindices(fec_trans(i))) = 1;
      }
    }
    fec_trans = find(categories == 26);
    if (fec_trans.n_elem > 0) {
      for (int i = 0; i < static_cast<int>(fec_trans.n_elem); i ++) {
        fmat(eindices(fec_trans(i))) = 1;
      }
    }
  }
  
  // Splits fecundity transitions if they include survival portions
  arma::mat corr_mat = amat;
  arma::uvec z_indices = find(corr_mat == 0.0);
  int z_indicesnem = static_cast<int>(z_indices.n_elem);
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
    // Object minindex will only be > 9 if historical
    arma::vec size3 = as<arma::vec>(bambesque["size3"]);
    arma::vec size2 = as<arma::vec>(bambesque["size2"]);
    arma::vec size1 = as<arma::vec>(bambesque["size1"]);
    
    arma::uvec repstatus3 = as<arma::uvec>(bambesque["repstatus3"]);
    arma::uvec repstatus2 = as<arma::uvec>(bambesque["repstatus2"]);
    arma::uvec repstatus1 = as<arma::uvec>(bambesque["repstatus1"]);
    
    arma::uvec entrystatus3 = as<arma::uvec>(bambesque["entrystatus3"]);
    arma::uvec entrystatus2 = as<arma::uvec>(bambesque["entrystatus2"]);
    arma::uvec entrystatus1 = as<arma::uvec>(bambesque["entrystatus1"]);
    
    for (int i = 0; i < 16; i++) {
      arma::uvec currentguys = find(categories == histcatnums(i));
      int currentguysnem = static_cast<int>(currentguys.n_elem);
      
      if (histcatnums(i) == 20 || histcatnums(i) == 21 || histcatnums(i) == 22 ||
        histcatnums(i) == 26) { // Fecundity transitions
      
        // Splits transitions that are actually combos of fecundity and survival
        for (int j = 0; j < currentguysnem; j++) {
          int this_guy = eindices(currentguys(j));
          
          if (fec_fraction(this_guy) == 1.0) {
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
              if (size3(counter(0)) == size2(counter(0)) &&
                repstatus3(counter(0)) == repstatus2(counter(0)) &&
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
              
            } else if (size3(counter(0)) == size2(counter(0)) &&
              size2(counter(0)) == size1(counter(0))) {
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
              
            } else if (size3(counter(0)) > size2(counter(0)) &&
                size2(counter(0)) == size1(counter(0))) {
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
              
            } else if (size3(counter(0)) > size2(counter(0)) &&
              size2(counter(0)) > size1(counter(0))) {
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
              
            } else if (size3(counter(0)) < size2(counter(0)) &&
              size2(counter(0)) == size1(counter(0))) {
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
              
            } else if (size3(counter(0)) < size2(counter(0)) &&
              size2(counter(0)) < size1(counter(0))) {
              hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              histsums(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              
              if (e_amat(this_guy) > 0) {
                hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histpos(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              } else if (e_amat(this_guy) < 0) {
                hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histneg(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              }
              
            } else if (size3(counter(0)) == size2(counter(0)) &&
              size2(counter(0)) < size1(counter(0))) {
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
              
            } else if (size3(counter(0)) < size2(counter(0)) &&
              size2(counter(0)) > size1(counter(0))) {
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
      } else if (histcatnums(i) == 14 || histcatnums(i) == 15 ||
          histcatnums(i) == 17 || histcatnums(i) == 25) { // Shrinkage transitions
        
        arma::vec all_es = e_amat.elem(eindices(currentguys));
        arma::uvec all_es_pos = find(all_es > 0.0);
        arma::uvec all_es_neg = find(all_es < 0.0);
        int all_es_pos_num = static_cast<int>(all_es_pos.n_elem);
        int all_es_neg_num = static_cast<int>(all_es_neg.n_elem);
        
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
        
      } else if (histcatnums(i) == 10 || histcatnums(i) == 13 ||
          histcatnums(i) == 16 || histcatnums(i) == 23) { // Stasis transitions
        
        arma::vec all_es = e_amat.elem(eindices(currentguys));
        arma::uvec all_es_pos = find(all_es > 0.0);
        arma::uvec all_es_neg = find(all_es < 0.0);
        int all_es_pos_num = static_cast<int>(all_es_pos.n_elem);
        int all_es_neg_num = static_cast<int>(all_es_neg.n_elem);
        
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
        
      } else if (histcatnums(i) == 11 || histcatnums(i) == 12 ||
          histcatnums(i) == 18 || histcatnums(i) == 24) { // Growth transitions
        
        arma::vec all_es = e_amat.elem(eindices(currentguys));
        arma::uvec all_es_pos = find(all_es > 0.0);
        arma::uvec all_es_neg = find(all_es < 0.0);
        int all_es_pos_num = static_cast<int>(all_es_pos.n_elem);
        int all_es_neg_num = static_cast<int>(all_es_neg.n_elem);
        
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
    
    histout = DataFrame::create(Named("category") = histcats,
      _["elas"] = histsums, _["elas_pos"] = histpos, _["elas_neg"] = histneg);
    ahistout = DataFrame::create(Named("category") = ahistcats,
      _["elas"] = hc_ahistsums, _["elas_pos"] = hc_ahistpos,
      _["elas_neg"] = hc_ahistneg);
  } else {
    histout = R_NilValue;
    
    for (int i = 0; i < 4; i++) {
      arma::uvec currentguys = find(categories == ahistcatnums(i));
      arma::vec all_es = e_amat.elem(eindices(currentguys));
      arma::uvec all_es_pos = find(all_es > 0.0);
      arma::uvec all_es_neg = find(all_es < 0.0);
      int all_es_pos_num = static_cast<int>(all_es_pos.n_elem);
      int all_es_neg_num = static_cast<int>(all_es_neg.n_elem);
      
      double getoutofdodge = sum(all_es);
      ahistsums(i) += getoutofdodge;
      if (all_es_pos_num > 0) {
        ahistpos(i) += sum(all_es.elem(all_es_pos));
      }
      if (all_es_neg_num > 0) {
        ahistneg(i) += sum(all_es.elem(all_es_neg));
      }
    }
    
    ahistout = DataFrame::create(Named("category") = ahistcats,
      _["elas"] = ahistsums, _["elas_pos"] = ahistpos,
      _["elas_neg"] = ahistneg);
  }
  
  List output = List::create(Named("hist") = histout, _["ahist"] = ahistout);
  
  return output;
}

//' Creates Summary Data for Elasticity Matrix Inputs
//' 
//' Function \code{demolition3sp()} sums elasticity values from elasticity
//' matrices, and LTRE contributions from LTRE and sLTRE matrices, according to
//' the categories developed by functions \code{bambi2()} and \code{bambi3()}.
//' It requires \code{dgCMatrix} class inputs.
//' 
//' @name .demolition3sp
//' 
//' @param e_amat A single elasticity, LTRE, or sLTRE matrix of class
//' \code{dgCMatrix}.
//' @param bambesque This is the output from \code{bambi2()} or \code{bambi3()}
//' corresponding to the current lefkoMat object. The format is a data frame
//' giving the indices and characteristics of all predicted potential non-zero
//' elements in the supplied matrix.
//' @param amat_ The A matrix corresponding to \code{e_amat}. If not supplied,
//' then only \code{bambesque} is used to determine transition categories. If
//' provided, then fecundity transitions may be split between fecundity and
//' survival portions. Must also be of class \code{dgCMatrix}.
//' @param fmat_ The F matrix corresponding to \code{e_amat}. If not supplied,
//' then only \code{bambesque} is used to determine transition categories. If
//' provided, then fecundity transitions may be split between fecundity and
//' survival portions. Must also be of class \code{dgCMatrix}.
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
// [[Rcpp::export(.demolition3sp)]]
List demolition3sp(const arma::sp_mat& e_amat, const DataFrame& bambesque,
  Nullable<arma::sp_mat> amat_ = R_NilValue,
  Nullable<arma::sp_mat> fmat_ = R_NilValue) {
  
  arma::uvec eindices = as<arma::uvec>(bambesque["index"]);
  arma::uvec categories = as<arma::uvec>(bambesque["transition"]);
  
  int e_amatsize = static_cast<int>(e_amat.n_elem);
  int e_amatrows = static_cast<int>(e_amat.n_rows);
  int maxelem = static_cast<int>(eindices.max());
  int minindex = static_cast<int>(categories.min());
  
  arma::sp_mat amat;
  arma::sp_mat fmat;
  
  if (maxelem > e_amatsize) {
    throw Rcpp::exception("Supplied info does not correspond to input matrices.",
      false);
  }
  
  if (amat_.isNotNull() && fmat_.isNotNull()) {
    arma::sp_mat amat_temp = as<arma::sp_mat>(amat_);
    arma::sp_mat fmat_temp = as<arma::sp_mat>(fmat_);
    
    amat = amat_temp;
    fmat = fmat_temp;
    
  } else {
    arma::sp_mat fmat1(e_amatrows, e_amatrows);
    fmat = fmat1;
    
    arma::uvec fec_trans = find(categories == 4);
    if (fec_trans.n_elem > 0) {
      for (int i = 0; i < static_cast<int>(fec_trans.n_elem); i ++) {
        fmat(eindices(fec_trans(i))) = 1;
      }
    }
    fec_trans = find(categories == 20);
    if (fec_trans.n_elem > 0) {
      for (int i = 0; i < static_cast<int>(fec_trans.n_elem); i ++) {
        fmat(eindices(fec_trans(i))) = 1;
      }
    }
    fec_trans = find(categories == 21);
    if (fec_trans.n_elem > 0) {
      for (int i = 0; i < static_cast<int>(fec_trans.n_elem); i ++) {
        fmat(eindices(fec_trans(i))) = 1;
      }
    }
    fec_trans = find(categories == 22);
    if (fec_trans.n_elem > 0) {
      for (int i = 0; i < static_cast<int>(fec_trans.n_elem); i ++) {
        fmat(eindices(fec_trans(i))) = 1;
      }
    }
    fec_trans = find(categories == 26);
    if (fec_trans.n_elem > 0) {
      for (int i = 0; i < static_cast<int>(fec_trans.n_elem); i ++) {
        fmat(eindices(fec_trans(i))) = 1;
      }
    }
  }
  
  // Splits fecundity transitions if they include survival portions
  arma::sp_mat fec_fraction(e_amatrows, e_amatrows);
  arma::sp_mat corr_mat = amat;
  
  arma::uvec z_indices = find(corr_mat);
  int z_indicesnem = static_cast<int>(z_indices.n_elem);
  
  for (int i = 0; i < z_indicesnem; i++) {
    fec_fraction(i) = fmat(i) / amat(i);
  }
  
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
    // Object minindex will only be > 9 if historical
    arma::vec size3 = as<arma::vec>(bambesque["size3"]);
    arma::vec size2 = as<arma::vec>(bambesque["size2"]);
    arma::vec size1 = as<arma::vec>(bambesque["size1"]);
    
    arma::uvec repstatus3 = as<arma::uvec>(bambesque["repstatus3"]);
    arma::uvec repstatus2 = as<arma::uvec>(bambesque["repstatus2"]);
    arma::uvec repstatus1 = as<arma::uvec>(bambesque["repstatus1"]);
    
    arma::uvec entrystatus3 = as<arma::uvec>(bambesque["entrystatus3"]);
    arma::uvec entrystatus2 = as<arma::uvec>(bambesque["entrystatus2"]);
    arma::uvec entrystatus1 = as<arma::uvec>(bambesque["entrystatus1"]);
    
    for (int i = 0; i < 16; i++) {
      arma::uvec currentguys = find(categories == histcatnums(i));
      int currentguysnem = static_cast<int>(currentguys.n_elem);
      
      if (histcatnums(i) == 20 || histcatnums(i) == 21 || histcatnums(i) == 22 ||
        histcatnums(i) == 26) { // Fecundity transitions
      
        // Splits transitions that are actually combos of fecundity and survival
        for (int j = 0; j < currentguysnem; j++) {
          int this_guy = eindices(currentguys(j));
          
          if (fec_fraction(this_guy) == 1.0) {
            hc_ahistsums(3) += (e_amat(this_guy));
            histsums(i) += (e_amat(this_guy));
            
            if (e_amat(this_guy) > 0) {
              hc_ahistpos(3) += (e_amat(this_guy));
              histpos(i) += (e_amat(this_guy));
            } else if (e_amat(this_guy) < 0) {
              hc_ahistneg(3) += (e_amat(this_guy));
              histneg(i) += (e_amat(this_guy));
            }
            
          } else if (fec_fraction(this_guy) > 0.0) {
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
              if (size3(counter(0)) == size2(counter(0)) &&
                repstatus3(counter(0)) == repstatus2(counter(0)) &&
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
              
            } else if (size3(counter(0)) == size2(counter(0)) &&
              size2(counter(0)) == size1(counter(0))) {
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
              
            } else if (size3(counter(0)) > size2(counter(0)) &&
                size2(counter(0)) == size1(counter(0))) {
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
              
            } else if (size3(counter(0)) > size2(counter(0)) &&
              size2(counter(0)) > size1(counter(0))) {
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
              
            } else if (size3(counter(0)) < size2(counter(0)) &&
              size2(counter(0)) == size1(counter(0))) {
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
              
            } else if (size3(counter(0)) < size2(counter(0)) &&
              size2(counter(0)) < size1(counter(0))) {
              hc_ahistsums(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              histsums(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              
              if (e_amat(this_guy) > 0) {
                hc_ahistpos(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histpos(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              } else if (e_amat(this_guy) < 0) {
                hc_ahistneg(2) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
                histneg(5) += (e_amat(this_guy) * (1 - fec_fraction(this_guy)));
              }
              
            } else if (size3(counter(0)) == size2(counter(0)) &&
              size2(counter(0)) < size1(counter(0))) {
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
              
            } else if (size3(counter(0)) < size2(counter(0)) &&
              size2(counter(0)) > size1(counter(0))) {
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
      } else if (histcatnums(i) == 14 || histcatnums(i) == 15 ||
          histcatnums(i) == 17 || histcatnums(i) == 25) { // Shrinkage transitions
        
        arma::vec all_es (static_cast<int>(currentguys.n_elem), fill::zeros);
        for (int m = 0; m < static_cast<int>(currentguys.n_elem); m++) {
          all_es(m) = e_amat(eindices(currentguys(m)));
        }
        
        arma::uvec all_es_pos = find(all_es > 0.0);
        arma::uvec all_es_neg = find(all_es < 0.0);
        int all_es_pos_num = static_cast<int>(all_es_pos.n_elem);
        int all_es_neg_num = static_cast<int>(all_es_neg.n_elem);
        
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
        
      } else if (histcatnums(i) == 10 || histcatnums(i) == 13 ||
          histcatnums(i) == 16 || histcatnums(i) == 23) { // Stasis transitions
        
        arma::vec all_es (static_cast<int>(currentguys.n_elem), fill::zeros);
        for (int m = 0; m < static_cast<int>(currentguys.n_elem); m++) {
          all_es(m) = e_amat(eindices(currentguys(m)));
        }
        
        arma::uvec all_es_pos = find(all_es > 0.0);
        arma::uvec all_es_neg = find(all_es < 0.0);
        int all_es_pos_num = static_cast<int>(all_es_pos.n_elem);
        int all_es_neg_num = static_cast<int>(all_es_neg.n_elem);
        
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
        
      } else if (histcatnums(i) == 11 || histcatnums(i) == 12 ||
          histcatnums(i) == 18 || histcatnums(i) == 24) { // Growth transitions
        
        arma::vec all_es (static_cast<int>(currentguys.n_elem), fill::zeros);
        for (int m = 0; m < static_cast<int>(currentguys.n_elem); m++) {
          all_es(m) = e_amat(eindices(currentguys(m)));
        }
        
        arma::uvec all_es_pos = find(all_es > 0.0);
        arma::uvec all_es_neg = find(all_es < 0.0);
        int all_es_pos_num = static_cast<int>(all_es_pos.n_elem);
        int all_es_neg_num = static_cast<int>(all_es_neg.n_elem);
        
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
    
    histout = DataFrame::create(Named("category") = histcats,
      _["elas"] = histsums, _["elas_pos"] = histpos, _["elas_neg"] = histneg);
    ahistout = DataFrame::create(Named("category") = ahistcats,
      _["elas"] = hc_ahistsums, _["elas_pos"] = hc_ahistpos,
      _["elas_neg"] = hc_ahistneg);
  } else {
    histout = R_NilValue;
    
    for (int i = 0; i < 4; i++) {
      arma::uvec currentguys = find(categories == ahistcatnums(i));
      
      arma::vec all_es (static_cast<int>(currentguys.n_elem), fill::zeros);
      for (int m = 0; m < static_cast<int>(currentguys.n_elem); m++) {
        all_es(m) = e_amat(eindices(currentguys(m)));
      }
        
      arma::uvec all_es_pos = find(all_es > 0.0);
      arma::uvec all_es_neg = find(all_es < 0.0);
      int all_es_pos_num = static_cast<int>(all_es_pos.n_elem);
      int all_es_neg_num = static_cast<int>(all_es_neg.n_elem);
      
      double getoutofdodge = sum(all_es);
      ahistsums(i) += getoutofdodge;
      if (all_es_pos_num > 0) {
        ahistpos(i) += sum(all_es.elem(all_es_pos));
      }
      if (all_es_neg_num > 0) {
        ahistneg(i) += sum(all_es.elem(all_es_neg));
      }
    }
    
    ahistout = DataFrame::create(Named("category") = ahistcats,
      _["elas"] = ahistsums, _["elas_pos"] = ahistpos,
      _["elas_neg"] = ahistneg);
  }
  
  List output = List::create(Named("hist") = histout, _["ahist"] = ahistout);
  
  return output;
}

//' Estimate Deterministic Population Growth Rate As Dominant Eigenvalue
//' 
//' Function \code{lambda3()} is a generic function that returns the dominant
//' eigenvalue of a matrix, set of dominant eigenvalues of a set of matrices,
//' or set of dominant eigenvalues for a \code{lefkoMat} object. It can handle
//' large and sparse matrices supplied as \code{lefkoMat} objects or as
//' individual matrices, and can be used with large historical matrices, IPMs, 
//' age x stage matrices, as well as smaller ahistorical matrices.
//' 
//' @param mpm A lefkoMat object, a list of projection matrices, or a single
//' projection matrix.
//' @param force_sparse A logical value or string detailing whether to force
//' sparse matrix encoding for simple matrix input. Defaults to \code{"auto"},
//' which only forces sparse matrix coding if simple matrices are input that are
//' both sparse (i.e, percentage of matrix elements that are non-zero <= 50%)
//' and have more than 20 rows. Can also be set to \code{"yes"}, \code{"no"},
//' \code{TRUE}, or \code{FALSE}. Note that sparse matrix coding is always used
//' for \code{lefkoMat} objects with matrices in sparse format (class
//' \code{dgCMatrix}).
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
RObject lambda3(RObject& mpm, Nullable<RObject> force_sparse = R_NilValue) {
  
  RObject output;
  
  int sparse_check {2};
  bool matrix_class_input {false};
  
  if (force_sparse.isNotNull()) {
    if (is<LogicalVector>(force_sparse)) {
      LogicalVector sparse_ = as<LogicalVector>(force_sparse);
      
      if (static_cast<bool>(sparse_(0))) {
        sparse_check = 1; // Forced sparse matrix
      } else {
        sparse_check = 0; // Not forced sparse matrix
      }
    } else if (is<StringVector>(force_sparse)) {
      StringVector sparse_ = as<StringVector>(force_sparse);
      
      if (stringcompare_simple(String(sparse_(0)), "au", true)) {
        sparse_check = 2; // Auto decision
        
      } else if (stringcompare_simple(String(sparse_(0)), "y", true) ||
          stringcompare_simple(String(sparse_(0)), "t", true)) {
        sparse_check = 1; // Forced sparse matrix
        
      } else if (stringcompare_simple(String(sparse_(0)), "n", true) ||
          stringcompare_simple(String(sparse_(0)), "f", true)) {
        sparse_check = 0; // Not forced sparse matrix
      } else {  
        throw Rcpp::exception("Value entered for argument sparse not understood.",
          false);
      }
    } else {
      throw Rcpp::exception("Value entered for argument sparse not understood.",
        false);
    }
  }
  
  if (is<List>(mpm)) {
    List mpm_ = as<List>(mpm);
    
    CharacterVector mpm_names;
    if (mpm_.hasAttribute("names")) mpm_names = mpm_.attr("names");
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
      // List of matrices input
      
      if (is<NumericMatrix>(mpm_(0))) { 
        matrix_class_input = true;
      } else if (is<S4>(mpm_(0))) { 
        matrix_class_input = false;
      } else {
        throw Rcpp::exception("Object mpm list structure is not recognized.", false);
      }
      
      int no_matrices = mpm_.length();
      
      if (matrix_class_input && sparse_check == 2) {
        arma::mat a1 = mpm_[0];
        int mat_rows = a1.n_rows;
        int mat_cols = a1.n_cols;
        int total_elems = mat_rows * mat_cols;
        
        arma::uvec nonzeros = find(a1);
        int no_nonzeros = static_cast<int>(nonzeros.n_elem);
        
        double density = static_cast<double>(no_nonzeros) / static_cast<double>(total_elems);
        
        if (density <= 0.5 && total_elems > 399) {
          sparse_check = 1;
        } else {
          sparse_check = 0;
        }
      }
      
      NumericVector lambda_prog (no_matrices);
      
      for (int i = 0; i < no_matrices; i++) {
        if (sparse_check == 0 && matrix_class_input) {
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
          arma::sp_mat spAmat;
          
          if (matrix_class_input) { 
            arma::sp_mat spAmat_(as<arma::mat>(mpm_(i)));
            spAmat = spAmat_;
          } else { 
            spAmat = as<arma::sp_mat>(mpm_(i));
          }
          
          arma::cx_vec Aeigval;
          arma::cx_mat Aeigvecr;
          
          eigs_gen(Aeigval, Aeigvecr, spAmat, 1, "lr");
          
          arma::vec all_eigenvalues = real(Aeigval);
          double maxval = max(all_eigenvalues);
          
          arma::uvec max_elems = find(all_eigenvalues == maxval);
          if (max_elems.n_elem == 0) {
            throw Rcpp::exception("Eigen analysis failed.", false);
          }
          arma::vec chosen_eigvals = all_eigenvalues.elem(max_elems);
          arma::uvec pos_max_eigvals = find(chosen_eigvals > 0);
          
          if (pos_max_eigvals.n_elem > 0) lambda_prog(i) = chosen_eigvals(pos_max_eigvals(0));
          
        }
      }
      output = lambda_prog;
      
    } else {
      // lefkoMat input
      
      List A_list = mpm_["A"];
      DataFrame labels = as<DataFrame>(mpm_["labels"]);
      
      if (is<NumericMatrix>(A_list(0))) { 
        matrix_class_input = true;
      } else if (is<S4>(A_list(0))) { 
        matrix_class_input = false;
      } else {
        throw Rcpp::exception("Object mpm does not appear to contain matrices.", false);
      }
      
      int no_matrices = A_list.length();
      
      if (matrix_class_input && sparse_check == 2) {
        arma::mat a1 = A_list[0];
        int mat_rows = a1.n_rows;
        int mat_cols = a1.n_cols;
        int total_elems = mat_rows * mat_cols;
        
        arma::uvec nonzeros = find(a1);
        int no_nonzeros = static_cast<int>(nonzeros.n_elem);
        
        double density = static_cast<double>(no_nonzeros) /
          static_cast<double>(total_elems);
        
        if (density <= 0.5 && total_elems > 399) {
          sparse_check = 1;
        } else {
          sparse_check = 0;
        }
      }
      
      NumericVector lambda_prog (no_matrices);
      
      for (int i = 0; i < no_matrices; i++) {
        if (sparse_check == 0 && matrix_class_input) {
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
          arma::sp_mat spAmat;
          
          if (matrix_class_input) {
            arma::sp_mat spAmat_(as<arma::mat>(A_list(i)));
            spAmat = spAmat_;
          } else { 
            spAmat = as<arma::sp_mat>(A_list(i));
          }
          
          arma::cx_vec Aeigval;
          arma::cx_mat Aeigvecr;
          
          eigs_gen(Aeigval, Aeigvecr, spAmat, 1, "lr");
          
          arma::vec all_eigenvalues = real(Aeigval);
          double maxval = max(all_eigenvalues);
          
          arma::uvec max_elems = find(all_eigenvalues == maxval);
          if (max_elems.n_elem == 0) {
            throw Rcpp::exception("Eigen analysis failed.", false);
          }
          arma::vec chosen_eigvals = all_eigenvalues.elem(max_elems);
          arma::uvec pos_max_eigvals = find(chosen_eigvals > 0);
          
          if (pos_max_eigvals.n_elem > 0) {
            lambda_prog(i) = chosen_eigvals(pos_max_eigvals(0));
          }
        }
      }
      
      DataFrame new_out;
      CharacterVector l_pop = labels["pop"];
      CharacterVector l_patch = labels["patch"];
      
      int l_length = labels.length();
      
      if (l_length == 3) {
        CharacterVector l_year2 = labels["year2"];
        
        new_out = DataFrame::create(_["pop"] = l_pop, _["patch"] = l_patch,
          _["year2"] = l_year2, _["lambda"] = lambda_prog);
      } else {
        new_out = DataFrame::create(_["pop"] = l_pop, _["patch"] = l_patch,
          _["lambda"] = lambda_prog);
      }
      output = new_out;
    }
    
  } else if(is<NumericMatrix>(mpm)) {
    // Single matrix input
    
    arma::mat mpm_ = as<arma::mat>(mpm);
    
    if (sparse_check == 2) {
      int mat_rows = mpm_.n_rows;
      int mat_cols = mpm_.n_cols;
      int total_elems = mat_rows * mat_cols;
      
      arma::uvec nonzeros = find(mpm_);
      int no_nonzeros = static_cast<int>(nonzeros.n_elem);
      
      double density = static_cast<double>(no_nonzeros) /
        static_cast<double>(total_elems);
      
      if (density <= 0.5 && total_elems > 399) {
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
        throw Rcpp::exception("Eigen analysis failed.", false);
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
        throw Rcpp::exception("Eigen analysis failed.", false);
      }
      
      arma::vec chosen_eigvals = all_eigenvalues.elem(max_elems);
      arma::uvec pos_max_eigvals = find(chosen_eigvals > 0);
      
      if (pos_max_eigvals.n_elem > 0) lambda_prog(0) = chosen_eigvals(pos_max_eigvals(0));
    }
    
    output = lambda_prog;
    
  } else if(is<S4>(mpm)) {
    // Single sparse matrix input
    
    arma::sp_mat spAmat = as<arma::sp_mat>(mpm);
    
    NumericVector lambda_prog (1);
    arma::cx_vec Aeigval;
    arma::cx_mat Aeigvecr;
    
    eigs_gen(Aeigval, Aeigvecr, spAmat, 1, "lr");
    
    arma::vec all_eigenvalues = real(Aeigval);
    
    double maxval = max(all_eigenvalues);
    
    arma::uvec max_elems = find(all_eigenvalues == maxval);
    if (max_elems.n_elem == 0) {
      throw Rcpp::exception("Eigen analysis failed.", false);
    }
    
    arma::vec chosen_eigvals = all_eigenvalues.elem(max_elems);
    arma::uvec pos_max_eigvals = find(chosen_eigvals > 0);
    
    if (pos_max_eigvals.n_elem > 0) lambda_prog(0) = chosen_eigvals(pos_max_eigvals(0));
    
    output = lambda_prog;
    
  } else {
    throw Rcpp::exception("Object mpm does not appear to be an appropriate MPM.",
      false);
  }
  
  return output;
}
