#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include "LefkoUtils.h"

using namespace Rcpp;
using namespace arma;
using namespace LefkoUtils;


//' Calculate Actual Stage, Age, Stage-Pair, or Age-Stage Distributions
//' 
//' Function \code{actualstage3()} shows the frequencies and proportions of
//' each stage, stage pair, age-stage, or age in each year.
//' 
//' @name actualstage3
//' 
//' @param data A demographic dataset in hfv format.
//' @param check_stage A logical value indicating whether to assess frequencies
//' and proportions of stages. Defaults to \code{TRUE}.
//' @param check_age A logical value indicating whether to assess frequencies and
//' proportions of ages. Defaults to \code{FALSE}.
//' @param historical A logical value indicating whether the stage structure
//' should be ahistorical (\code{FALSE}) or historical (\code{TRUE}). Defaults to
//' \code{FALSE}.
//' @param year2 A string value indicating the name of the variable coding for
//' monitoring occasion at time \emph{t}. Defaults to \code{"year2"}.
//' @param indices A vector of three strings, indicating the stage indices for
//' times \emph{t}+1, \emph{t}, and \emph{t}-1, respectively, in \code{data}.
//' Defaults to \code{c("stage3index", "stage2index", "stage1index")}.
//' @param stagecol A vector of three strings, indicating the stage name columns
//' for times \emph{t}+1, \emph{t}, and \emph{t}-1, respectively, in \code{data}.
//' Defaults to \code{stagecol = c("stage3", "stage2", "stage1")}.
//' @param agecol A single string indicating the age of individuals in time
//' \emph{t}. Defaults to \code{"obsage"}.
//' @param remove_stage A string vector indicating the names of stages to remove
//' from consideration. Defaults to \code{"NotAlive"}.
//' @param t1_allow A string vector indicating which stages to be removed should
//' be allowed in the stage at time \emph{t}-1 portion of historical stage
//' pairs, if \code{historical = TRUE}. Defaults to \code{"NotAlive"}. Can also
//' be set to \code{"none"}.
//' 
//' @return A data frame with the following variables:
//' \item{rowid}{A string identifier term, equal to the monitoring occasion in
//' time \emph{t} and the stage index.}
//' \item{stageindex}{The stageframe index of the stage. Only output if
//' \code{check_stage = TRUE}.}
//' \item{stage}{The name of each stage, or \code{NA}. Only output if
//' \code{check_stage = TRUE}.}
//' \item{stage2}{The name of the stage in time \emph{t}. Only output if
//' \code{check_stage = TRUE}.}
//' \item{stage1}{The name of the stage in time \emph{t}-1, or \code{NA}. Only
//' output if \code{check_stage = TRUE}.}
//' \item{age}{The age at time \emph{t}. Only output if \code{check_age = TRUE}.}
//' \item{year2}{Monitoring occasion in time \emph{t}.}
//' \item{frequency}{The number of individuals in the respective stage and time.}
//' \item{actual_prop}{The proportion of individuals alive in time \emph{t} in
//' the respective stage.}
//' 
//' @section Notes:
//' This function produces frequencies and proportions of stages in hfv formatted
//' data using stage index variables rather than stage name variables, and so
//' requires the former. The latter is only required if the user wants to know
//' the associated stage names.
//' 
//' Frequencies and proportions will be calculated for all times, including the
//' last time, which is generally found in the \code{stage3} columns of the last
//' \code{year2} entry in object \code{data}. The default is to treat the
//' \code{year2} entry for that time as \code{max(year2) + 1}.
//' 
//' If \code{check_stage = TRUE} and \code{check_age = FALSE}, then this function
//' will assess frequencies and proportions of stages or historical stage-pairs.
//' If both \code{check_stage = TRUE} and \code{check_age = TRUE}, then this
//' function will assess frequencies and proportions of age-stages. If
//' \code{check_stage = FALSE} and \code{check_age = TRUE}, then the frequencies
//' and proportions of ages only will be assessed.
//' 
//' Note that no stageframe is required for this function to operate. Stage
//' names and their order are inferred directly from the object \code{data}.
//' 
//' @examples
//' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 3, 6, 11, 19.5)
//' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
//'   "XLg")
//' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
//' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
//' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
//' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
//' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1.5, 1.5, 3.5, 5)
//' comments <- c("Dormant seed", "1st yr protocorm", "2nd yr protocorm",
//'   "3rd yr protocorm", "Seedling", "Dormant adult",
//'   "Extra small adult (1 shoot)", "Small adult (2-4 shoots)",
//'   "Medium adult (5-7 shoots)", "Large adult (8-14 shoots)",
//'   "Extra large adult (>14 shoots)")
//' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector, 
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
//'   propstatus = propvector, immstatus = immvector, indataset = indataset, 
//'   binhalfwidth = binvec, comments = comments)
//' 
//' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004, 
//'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
//'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
//'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
//'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
//'   NRasRep = TRUE, age_offset = 4)
//' 
//' all_stage_props_ah <- actualstage3(cypraw_v1)
//' all_stage_props_h <- actualstage3(cypraw_v1, historical = TRUE)
//' all_stage_props_h_NANotAllow <- actualstage3(cypraw_v1, historical = TRUE,
//'   t1_allow = "none")
//' all_stage_props_as <- actualstage3(cypraw_v1, check_age = TRUE)
//' all_age_props <- actualstage3(cypraw_v1, check_stage = FALSE,
//'   check_age = TRUE)
//' 
//' @export actualstage3
// [[Rcpp::export(actualstage3)]]
List actualstage3(RObject data, bool check_stage = true, bool check_age = false,
  bool historical = false, Nullable<RObject> year2 = R_NilValue,
  Nullable<RObject> indices = R_NilValue, Nullable<RObject> stagecol = R_NilValue,
  Nullable<RObject> agecol= R_NilValue, Nullable<RObject> remove_stage = R_NilValue,
  Nullable<RObject> t1_allow = R_NilValue) {
  
  DataFrame data_;
  
  String year2_;
  String agecol_;
  StringVector remove_stage_;
  StringVector t1_allow_;
  
  int year2_num_;
  int agecol_num_;
  IntegerVector remove_stage_num_;
  IntegerVector remove_stage_num_t1a;
  IntegerVector t1_allow_num_;
  
  bool year2_num_yn = false;
  bool agecol_num_yn = false;
  bool remove_stage_num_yn = false;
  bool t1_allow_num_yn = false;
  
  StringVector indices_;
  StringVector stagecol_;
  
  IntegerVector indices_num_;
  IntegerVector stagecol_num_;
  
  bool indices_num_yn = false;
  bool stagecol_num_yn = false;
  
  int stage3index_;
  int stage2index_;
  int stage1index_;
  int stage3_;
  int stage2_;
  int stage1_;
  
  bool stages3_supplied = false;
  bool stages2_supplied = false;
  bool stages1_supplied = false;
  bool indices3_supplied = false;
  bool indices2_supplied = false;
  bool indices1_supplied = false;
  bool ages_supplied = false;
  bool years_supplied = false;
  
  if (!check_age && !check_stage) {
    throw Rcpp::exception("Options check_age and check_stage cannot both be FALSE.", false);
  }
  
  if (year2.isNotNull()) {
    if (is<StringVector>(year2)) {
      StringVector year2_long = as<StringVector>(year2);
      year2_ = year2_long(0);
      
    } else if (is<IntegerVector>(year2)) {
      IntegerVector year2_long = as<IntegerVector>(year2);
      year2_num_ = year2_long(0);
      year2_num_yn = true;
      
    } else {
      throw Rcpp::exception("Object year2 must be either a variable name or a column number from object data.", 
        false);
    }
    
  } else year2_ = "year2";
  
  if (agecol.isNotNull()) {
    if (is<StringVector>(agecol)) {
      StringVector agecol_long = as<StringVector>(agecol);
      agecol_ = agecol_long(0);
      
    } else if (is<IntegerVector>(agecol)) {
      IntegerVector agecol_long = as<IntegerVector>(agecol);
      agecol_num_ = agecol_long(0);
      agecol_num_yn = true;
      
    } else {
      throw Rcpp::exception("Object agecol must be either a variable name or a column number from object data.", 
        false);
    }
    
  } else agecol_ = "obsage";
  
  if (remove_stage.isNotNull()) {
    if (is<StringVector>(remove_stage)) {
      remove_stage_ = as<StringVector>(remove_stage);
      
      if (remove_stage_.length() == 1) {
        if (remove_stage_(0) == "none") {
          remove_stage_ = {};
          remove_stage_num_ = {};
          
        }
      }
      
    } else if (is<IntegerVector>(remove_stage)) {
      remove_stage_num_ = as<IntegerVector>(remove_stage);
      remove_stage_num_yn = true;
      
    } else {
      throw Rcpp::exception("Object remove_stage must be either a string of variable names or a vector of column number from object data.", 
        false);
    }
    
  } else remove_stage_ = {"NotAlive"};
  
  if (t1_allow.isNotNull()) {
    if (is<StringVector>(t1_allow)) {
      t1_allow_ = as<StringVector>(t1_allow);
      
      if (t1_allow_.length() == 1) {
        if (t1_allow_(0) == "none") {
          t1_allow_ = {};
          t1_allow_num_ = {};
          
        }
      }
      
    } else if (is<IntegerVector>(t1_allow)) {
      t1_allow_num_ = as<IntegerVector>(t1_allow);
      t1_allow_num_yn = true;
      
    } else {
      throw Rcpp::exception("Object t1_allow must be either a string of variable names or a vector of column number from object data.", false);
    }
  } else t1_allow_ = {"NotAlive"};
  
  if (indices.isNotNull()) {
    if (is<StringVector>(indices)) {
      indices_ = as<StringVector>(indices);
      
    } else if (is<IntegerVector>(indices)) {
      indices_num_ = as<IntegerVector>(indices);
      indices_num_yn = true;
      
    } else {
      throw Rcpp::exception("Object indices must be a vector of variable names or column numbers from object data.", 
        false);
    }
    
  } else indices_ = {"stage3index", "stage2index", "stage1index"};
  
  if (stagecol.isNotNull()) {
    if (is<StringVector>(stagecol)) {
      stagecol_ = as<StringVector>(stagecol);
      
    } else if (is<IntegerVector>(stagecol)) {
      stagecol_num_ = as<IntegerVector>(stagecol);
      stagecol_num_yn = true;
      
    } else {
      throw Rcpp::exception("Object stagecol must be a vector of variable names or column numbers from object data.", 
        false);
    }
    
  } else stagecol_ = {"stage3", "stage2", "stage1"};
  
  if (is<DataFrame>(data)) {
    data_ = as<DataFrame>(data);
    
  } else {
    throw Rcpp::exception("Object data must be a hfv data frame.", false);
  }
  
  // Check data frame for the correct variables, and pull out variables as vectors
  CharacterVector data_vars = data_.attr("names");
  int data_vars_num = data_vars.length();
  int data_rows = data_.nrows();
  int indices_length = indices_.length();
  int stagecol_length = stagecol_.length();
  
  if (historical) {
    if (check_age) {
      throw Rcpp::exception("Package lefko3 does not currently support historical age-by-stage analyses.", false);
    }
    
    if (indices_length < 3) {
      throw Rcpp::exception("Object indices must contain stage index variables for times t+1, t, and t-1 if historical = TRUE.", false);
    }
    
    if (stagecol_length < 3) {
      throw Rcpp::exception("Object stagecol must contain stage classifications for times t+1, t, and t-1 if historical = TRUE.", false);
    }
  } else {
    if (indices_length < 2) {
      throw Rcpp::exception("Object indices must contain stage index variables for times t+1 and t if historical = FALSE.", false);
    }
    
    if (stagecol_length < 2) {
      throw Rcpp::exception("Object stagecol must contain stage classifications for times t+1 and t if historical = FALSE.", false);
    }
  }
  
  // Loop to find column numbers of variables given as strings
  for (int i = 0; i < data_vars_num; i++) {
    if (!year2_num_yn) {
      if (stringcompare_hard(as<std::string>(data_vars(i)), year2_)) {
        year2_num_ = i;
        years_supplied = true;
      }
    }
    
    if (!agecol_num_yn) {
      if (stringcompare_hard(as<std::string>(data_vars(i)), agecol_)) {
        agecol_num_ = i;
        ages_supplied = true;
      }
    }
    
    if (!indices_num_yn) {
      if (stringcompare_hard(as<std::string>(data_vars(i)), as<std::string>(indices_(0)))) {
        stage3index_ = i;
        indices3_supplied = true;
      }
      
      if (stringcompare_hard(as<std::string>(data_vars(i)), as<std::string>(indices_(1)))) {
        stage2index_ = i;
        indices2_supplied = true;
      }
      
      if (historical) {
        if (stringcompare_hard(as<std::string>(data_vars(i)), as<std::string>(indices_(2)))) {
          stage1index_ = i;
          indices1_supplied = true;
        }
      }
    }
    
    if (!stagecol_num_yn) {
      if (stringcompare_hard(as<std::string>(data_vars(i)), as<std::string>(stagecol_(0)))) {
        stage3_ = i;
        stages3_supplied = true;
      }
      if (stringcompare_hard(as<std::string>(data_vars(i)), as<std::string>(stagecol_(1)))) {
        stage2_ = i;
        stages2_supplied = true;
      }
      
      if (historical) {
        if (stringcompare_hard(as<std::string>(data_vars(i)), as<std::string>(stagecol_(2)))) {
          stage1_ = i;
          stages1_supplied = true;
        }
      }
    }
  }
  
  // Check if variable column numbers actually fit within data frame provided
  if (year2_num_yn) {
    if (year2_num_ < 0 || year2_num_ >= data_vars_num) {
      throw Rcpp::exception("Variable year2 not defined within data frame as input.", false);
      
    } else {
      years_supplied = true;
    }
  }
  
  if (agecol_num_yn) {
    if (agecol_num_ < 0 || agecol_num_ >= data_vars_num) {
      throw Rcpp::exception("Variable agecol not defined within data frame as input.", false);
      
    } else {
      ages_supplied = true;
    }
  }
  
  if (indices_num_yn) {
    if (indices_num_(0) >= 0 && indices_num_(0) < data_vars_num) {
      stage3index_ = indices_num_(0);
      indices3_supplied = true;
    }
    
    if (indices_num_(1) >= 0 && indices_num_(1) < data_vars_num) {
      stage2index_ = indices_num_(1);
      indices2_supplied = true;
    }
    if (historical) {
      if (indices_num_(2) >= 0 && indices_num_(2) < data_vars_num) {
        stage1index_ = indices_num_(2);
        indices1_supplied = true;
      }
    }
  }
  
  if (stagecol_num_yn) {
    if (stagecol_num_(0) >= 0 && stagecol_num_(0) < data_vars_num) {
      stage3_ = stagecol_num_(0);
      stages3_supplied = true;
    }
    
    if (stagecol_num_(1) >= 0 && stagecol_num_(1) < data_vars_num) {
      stage2_ = stagecol_num_(1);
      stages2_supplied = true;
    }
    
    if (historical) {
      if (stagecol_num_(2) >= 0 && stagecol_num_(2) < data_vars_num) {
        stage1_ = stagecol_num_(2);
        stages1_supplied = true;
      }
    }
  }
  
  // New data frame variable assignments
  arma::ivec data_year2;
  if (years_supplied) data_year2 = as<arma::ivec>(data_[year2_num_]);
  
  if (data_year2.n_elem != data_rows) {
    throw Rcpp::exception("Object year2 does not contain a valid variable.", false);
  }
  
  IntegerVector data_agecol;
  
  if (check_age) {
    if (ages_supplied) {
      data_agecol = data_[agecol_];
      
      if (data_agecol.length() != data_rows) {
        throw Rcpp::exception("Object agecol does not contain a valid variable for age at time t.", false);
      }
    } else {
      throw Rcpp::exception("Object agecol must be provided if check_age = TRUE.", false);
    }
  }
  
  CharacterVector data_stage3;
  CharacterVector data_stage2;
  CharacterVector data_stage1;
  
  arma::ivec data_stage3index;
  arma::ivec data_stage2index;
  arma::ivec data_stage1index;
  
  if (check_stage) {
    if (stages3_supplied) {
      data_stage3 = as<CharacterVector>(data_[stage3_]);
      
      // Rcout << "\n Length of data_stage3: " << data_stage3.length() << "\n";
      // Rcout << "First entry in data_stage3: " << data_stage3(0) << "\n";
      
      if (data_stage3.length() != data_rows) {
        throw Rcpp::exception("Object stagecol does not contain a valid variable for stage at time t+1.", false);
      }
    }
    if (stages2_supplied) {
      data_stage2 = as<CharacterVector>(data_[stage2_]);
      
      // Rcout << "\n Length of data_stage2: " << data_stage2.length() << "\n";
      // Rcout << "First entry in data_stage2: " << data_stage2(0) << "\n";
      
      if (data_stage2.length() != data_rows) {
        throw Rcpp::exception("Object stagecol does not contain a valid variable for stage at time t.", false);
      }
    }
    
    if (indices3_supplied) {
      data_stage3index = as<arma::ivec>(data_[stage3index_]);
      
      // Rcout << "\n Length of data_stage3index: " << data_stage3index.n_elem << "\n";
      // Rcout << "First entry in data_stage3index: " << data_stage3index(0) << "\n";
      
      if (data_stage3index.n_elem != data_rows) {
        throw Rcpp::exception("Object indices does not contain a valid variable for stage at time t+1.", false);
      }
    }
    if (indices2_supplied) {
      data_stage2index = as<arma::ivec>(data_[stage2index_]);
      
      // Rcout << "First entry in data_stage2index: " << data_stage2index(0) << "\n";
      // Rcout << "\n Length of data_stage2index: " << data_stage2index.n_elem << "\n";
      
      if (data_stage2index.n_elem != data_rows) {
        throw Rcpp::exception("Object indices does not contain a valid variable for stage at time t.", false);
      }
    }
    
    if (!stages3_supplied && !indices3_supplied) {
      throw Rcpp::exception("Objects indices and/or stagecol must be provided if check_stage = TRUE.", false);
    }
    if (!stages2_supplied && !indices2_supplied) {
      throw Rcpp::exception("Objects indices and/or stagecol must be provided if check_stage = TRUE.", false);
    }
    
    if (historical) {
      if (stages1_supplied) {
        data_stage1 = data_[stage1_];
        
        if (data_stage1.length() != data_rows) {
          throw Rcpp::exception("Object stagecol does not contain a valid variable for stage at time t-1.", false);
        }
      }
      if (indices3_supplied) {
        data_stage1index = as<arma::ivec>(data_[stage1index_]);
        
        if (data_stage1index.n_elem != data_rows) {
          throw Rcpp::exception("Object indices does not contain a valid variable for stage at time t-1.", false);
        }
      }
      
      if (!stages1_supplied && !indices1_supplied) {
        throw Rcpp::exception("Objects indices and/or stagecol must be provided if check_stage = TRUE.", false);
      }
    }
  }
  
  // Now we develop vectors of all values
  arma::ivec all_year2 = unique(data_year2);
  int year2_num = all_year2.n_elem;
  int years_num = year2_num + 1;
  
  arma::ivec all_years (years_num);
  for (int i = 0; i < year2_num; i++) {
    all_years(i) = all_year2(i);
  }
  all_years(year2_num) = all_years(year2_num - 1) + 1;
  
  CharacterVector all_stage3;
  CharacterVector all_stage2;
  CharacterVector all_stage1;
  CharacterVector all_stages;
  CharacterVector all_stages_t1a;
  
  arma::ivec all_stage3index;
  arma::ivec all_stage2index;
  arma::ivec all_stage1index;
  arma::ivec all_stageindices;
  arma::ivec all_stageindices_t1a;
  
  int num_stages = 0;
  int num_stages_t1a = 0;
  
  if (check_stage) {
    if (indices3_supplied && indices2_supplied) {
      all_stage3index = unique(data_stage3index);
      all_stage2index = unique(data_stage2index);
      arma::ivec all_stage32indices = unique(join_cols(all_stage3index, all_stage2index));
      
      if (historical) {
        all_stage1index = unique(data_stage1index);
        all_stageindices = unique(join_cols(all_stage32indices, all_stage1index));
      } else {
        all_stageindices = all_stage32indices;
      }
      
      int s_length = all_stageindices.n_elem;
      // Rcout << "\n Length of all_stageindices: " << s_length << "\n";
      
      if (stages3_supplied && stages2_supplied) {
        CharacterVector temp_stages (s_length);
        for (int i = 0; i < s_length; i++) {
          arma::uvec found_stages3 = find(data_stage3index == all_stageindices(i));
          
          if (found_stages3.n_elem > 0) {
            temp_stages(i) = data_stage3(found_stages3(0));
          } else {
            arma::uvec found_stages2 = find(data_stage2index == all_stageindices(i));
            
            if (found_stages2.n_elem > 0) {
              temp_stages(i) = data_stage2(found_stages2(0));
              
            } else if (historical && stages1_supplied) {
              arma::uvec found_stages1 = find(data_stage1index == all_stageindices(i));
              
              if (found_stages1.n_elem > 0) {
                temp_stages(i) = data_stage1(found_stages1(0));
              }
            }
          }
        }
        all_stages = temp_stages;
        num_stages = all_stages.length();
        num_stages_t1a = all_stages.length();
      } else {
        num_stages = all_stageindices.n_elem;
        num_stages_t1a = all_stageindices.n_elem;
      }
    } else if (stages3_supplied && stages2_supplied) {
      all_stage3 = sort_unique(data_stage3);
      all_stage2 = sort_unique(data_stage2);
      
      CharacterVector all_stage32 = union_(all_stage3, all_stage2);
      
      if (historical) {
        all_stage1 = sort_unique(data_stage1);
        
        all_stages = union_(all_stage32, all_stage1);
        
      } else {
        all_stages = all_stage32;
      }
      
      num_stages = all_stages.length();
      num_stages_t1a = all_stages.length();
      IntegerVector temp_stageindices = seq(1, num_stages);
      all_stageindices = as<arma::ivec>(temp_stageindices);
    }
    
    if (remove_stage_.length() != 0 && !remove_stage_num_yn) {
      IntegerVector transfer_remove_stage_num (remove_stage_.length());
      IntegerVector transfer_remove_stage_num_t1a (remove_stage_.length());
      int internal_counter = 0;
      int internal_counter_t1a = 0;
      
      for (int i = 0; i < num_stages; i++) {
        for (int j = 0; j < remove_stage_.length(); j++) {
          if (all_stages(i) == remove_stage_(j)) {
            transfer_remove_stage_num(internal_counter) = i;
            remove_stage_num_yn = true;
            internal_counter++;
            
            if (historical) {
              if (t1_allow_.length() > 0) {
                
                bool found_stage = false;
                for (int k = 0; k < t1_allow_.length(); k++) {
                  if (all_stages(i) == t1_allow_(k)) {
                    found_stage = true;
                    // Rcout << "Found stage (1) " << all_stages(i) << "\n";
                    
                    transfer_remove_stage_num_t1a = shrink(transfer_remove_stage_num_t1a);
                  }
                }
                
                if (!found_stage) {
                  transfer_remove_stage_num_t1a(internal_counter_t1a) = i;
                  internal_counter_t1a++;
                }
              } else if (t1_allow_num_yn) {
                bool found_stage = false;
                for (int k = 0; k < t1_allow_num_.length(); k++) {
                  if (all_stageindices(i) == t1_allow_num_(k)) {
                    found_stage = true;
                    // Rcout << "Found stage (2) " << all_stages(i) << "\n";
                    
                    transfer_remove_stage_num_t1a = shrink(transfer_remove_stage_num_t1a);
                  }
                }
                
                if (!found_stage) {
                  transfer_remove_stage_num_t1a(internal_counter_t1a) = i;
                  internal_counter_t1a++;
                }
              }
            }
          }
        }
      }
      remove_stage_num_ = transfer_remove_stage_num;
      if (historical) {
        remove_stage_num_t1a = clone(transfer_remove_stage_num_t1a);
      } else {
        remove_stage_num_t1a = clone(transfer_remove_stage_num);
      }
      
      
      
      // Error checks
      // Rcout << "all_stages\n";
      // for (int i=0; i < all_stages.length(); i++) {
      //   Rcout << "stage " << i << " " << all_stages(i) << "\n";
      // }
      // 
      // Rcout << "remove_stage_num_\n";
      // for (int i=0; i < remove_stage_num_.length(); i++) {
      //   Rcout << "remove_stage_num_ " << i << " " << all_stages(remove_stage_num_(i)) << "\n";
      // }
      // 
      // Rcout << "remove_stage_num_\n";
      // for (int i=0; i < remove_stage_num_t1a.length(); i++) {
      //   Rcout << "remove_stage_num_t1a " << i << " " << all_stages(remove_stage_num_t1a(i)) << "\n";
      // }
      
      if (!remove_stage_num_yn) {
        Rf_warningcall(R_NilValue,
          "Stage(s) provided in option remove_stage could not be found and so will be ignored.");
      }
      
      IntegerVector unique_rsn = unique(remove_stage_num_);
      IntegerVector unique_rsn_t1a = unique(remove_stage_num_t1a);
      
      if (unique_rsn.length() < remove_stage_num_.length() || 
          unique_rsn_t1a.length() < remove_stage_num_t1a.length()) {
        Rf_warningcall(R_NilValue, "Some stages supplied in either option remove_stage or option t1_allow are duplicates and will be ignored.");
        
        remove_stage_num_ = clone(unique_rsn);
        unique_rsn_t1a = clone(unique_rsn_t1a);
      }
      
      // Redetermination of the number of stages to be used
      IntegerVector sorted_remst = remove_stage_num_.sort();
      IntegerVector rev_sorted_remst = rev(sorted_remst);
      
      CharacterVector rep_all_stages = clone(all_stages);
      CharacterVector orig_all_stages = clone(all_stages);
      IntegerVector rep_all_stageindices = as<IntegerVector>(wrap(all_stageindices));
      IntegerVector orig_rep_all_stageindices = clone(rep_all_stageindices);
      
      for (int i = 0; i < rev_sorted_remst.length(); i++) {
        rep_all_stages.erase(rev_sorted_remst(i));
        rep_all_stageindices.erase(rev_sorted_remst(i));
      }
      all_stages = rep_all_stages;
      all_stageindices = as<arma::ivec>(rep_all_stageindices);
      num_stages = all_stages.length();
      
      if (historical) {
        IntegerVector sorted_remst_t1a = remove_stage_num_t1a.sort();
        IntegerVector rev_sorted_remst_t1a = rev(sorted_remst_t1a);
        
        CharacterVector rep_all_stages_t1a = clone(orig_all_stages);
        IntegerVector rep_all_stageindices_t1a = clone(orig_rep_all_stageindices);
        
        for (int i = 0; i < rev_sorted_remst_t1a.length(); i++) {
          rep_all_stages_t1a.erase(rev_sorted_remst_t1a(i));
          rep_all_stageindices_t1a.erase(rev_sorted_remst_t1a(i));
        }
        all_stages_t1a = rep_all_stages_t1a;
        all_stageindices_t1a = as<arma::ivec>(rep_all_stageindices_t1a);
        num_stages_t1a = all_stages_t1a.length();
      }
    } else if (historical) {
      all_stages_t1a = all_stages;
      all_stageindices_t1a = all_stageindices;
      num_stages_t1a = all_stages.length();
    }
  }
  
  IntegerVector all_ages;
  int num_ages = 0;
  
  if (check_age) {
    IntegerVector all_age2 = sort_unique(data_agecol);
    int num_age2 = all_age2.length();
    num_ages = num_age2 + 1;
    
    IntegerVector all_ages_ (num_age2 + 1);
    
    for (int i = 0; i < num_age2; i++) {
      all_ages_(i) = all_age2(i);
    }
    all_ages_(num_age2) = all_ages_(num_age2 - 1) + 1;
    
    all_ages = all_ages_;
  }
  
  // New data frame variables
  int new_df_rows = 0;
  
  if (check_stage && !check_age) {
    if (!historical) {
      new_df_rows = years_num * num_stages;
      
    } else {
      new_df_rows = years_num * num_stages * num_stages_t1a;
    }
    
  } else if (!check_stage && check_age) {
    new_df_rows = years_num * num_ages;
    
  } else if (check_stage && check_age) {
    new_df_rows = years_num * num_stages * num_ages;
  }
  
  IntegerVector new_row_names (new_df_rows);
  CharacterVector new_row_id (new_df_rows);
  IntegerVector new_stageindex (new_df_rows);
  CharacterVector new_stage (new_df_rows);
  CharacterVector new_stage2 (new_df_rows);
  CharacterVector new_stage1 (new_df_rows);
  IntegerVector new_age (new_df_rows);
  IntegerVector new_year (new_df_rows);
  IntegerVector new_frequency (new_df_rows);
  NumericVector new_actualprop (new_df_rows);
  
  // Rcout << "\n new_df_rows: " << new_df_rows << "\n";
  
  if (check_stage && !check_age) {
    if (!historical) {
      IntegerVector first_years = rep(all_years(0), num_stages);
      CharacterVector first_stages2 = clone(all_stages);
      IntegerVector first_indices = as<IntegerVector>(wrap(all_stageindices));
      
      CharacterVector dud_stages (num_stages);
      for (int i = 0; i < num_stages; i++) {
        dud_stages(i) = "";
      }
      
      CharacterVector first_stages1 (num_stages);
      for (int i = 0; i < num_stages; i++) {
        first_stages1(i) = "";
      }
      
      for (int i = 1; i < years_num; i++) {
        IntegerVector next_years = rep(all_years(i), num_stages);
        first_years = concat_int(first_years, next_years);
        
        first_stages2 = concat_str(first_stages2, all_stages);
        first_stages1 = concat_str(first_stages1, dud_stages);
        
        first_indices = concat_int(first_indices, as<IntegerVector>(wrap(all_stageindices)));
      }
      new_year = first_years;
      
      new_stage2 = first_stages2;
      new_stage1 = first_stages1;
      new_stageindex = first_indices;
      
    } else {
      IntegerVector first_years = rep(all_years(0), num_stages * num_stages_t1a);
      CharacterVector first_stages2 = clone(all_stages);
      
      CharacterVector first_stages1_short (num_stages_t1a);
      CharacterVector first_stages1;
      
      for (int i = 0; i < num_stages_t1a; i++) {
        for (int j = 0; j < num_stages_t1a; j++) {
          first_stages1_short(j) = all_stages_t1a(i);
        }
        
        if (i > 0) {
          first_stages1 = concat_str(first_stages1, first_stages1_short);
          first_stages2 = concat_str(first_stages2, all_stages);
        } else {
          first_stages1 = clone(first_stages1_short);
        }
      }
      
      CharacterVector new_first_stages1 = clone(first_stages1);
      CharacterVector new_first_stages2 = clone(first_stages2);
      
      for (int i = 1; i < years_num; i++) {
        IntegerVector next_years = rep(all_years(i), num_stages * num_stages_t1a);
        first_years = concat_int(first_years, next_years);
        
        new_first_stages1 = concat_str(new_first_stages1, first_stages1);
        new_first_stages2 = concat_str(new_first_stages2, first_stages2);
      }
      
      new_year = first_years;
      
      new_stage2 = new_first_stages2;
      new_stage1 = new_first_stages1;
    }
    
  } else if (!check_stage && check_age) {
  
    IntegerVector first_years = rep(all_years(0), num_ages);
    IntegerVector first_ages = clone(all_ages);
    
    for (int i = 1; i < years_num; i++) {
      IntegerVector next_years = rep(all_years(i), num_ages);
      first_years = concat_int(first_years, next_years);
      
      first_ages = concat_int(first_ages, all_ages);
    }
    
    new_year = first_years;
    new_age = first_ages;
    
  } else if (check_stage && check_age) {
    
    CharacterVector first_stages2 = clone(all_stages);
    IntegerVector first_indices = as<IntegerVector>(wrap(all_stageindices));
    
    CharacterVector dud_stages (num_stages);
    for (int i = 0; i < num_stages; i++) {
      dud_stages(i) = "";
    }
    
    CharacterVector first_stages1 (num_stages);
    for (int i = 0; i < num_stages; i++) {
      first_stages1(i) = "";
    }
    
    IntegerVector first_ages = rep(all_ages(0), num_stages);
    for (int i = 1; i < num_ages; i++) {
      IntegerVector next_ages = rep(all_ages(i), num_stages);
      first_ages = concat_int(first_ages, next_ages);
      
      first_stages2 = concat_str(first_stages2, all_stages);
      first_stages1 = concat_str(first_stages1, dud_stages);
    }
    
    IntegerVector first_years = rep(all_years(0), num_stages * num_ages);
    IntegerVector cfirst_ages = clone(first_ages);
    CharacterVector cfirst_stages2 = clone(first_stages2);
    CharacterVector cfirst_stages1 = clone(first_stages1);
    
    for (int i = 1; i < years_num; i++) {
      IntegerVector next_years = rep(all_years(i), num_stages * num_ages);
      first_years = concat_int(first_years, next_years);
      
      first_ages = concat_int(first_ages, cfirst_ages);
      
      first_stages2 = concat_str(first_stages2, cfirst_stages2);
      first_stages1 = concat_str(first_stages1, cfirst_stages1);
      
      first_indices = concat_int(first_indices, as<IntegerVector>(wrap(all_stageindices)));
    }
    
    new_year = first_years;
    new_age = first_ages;
    new_stage2 = first_stages2;
    new_stage1 = first_stages1;
    new_stageindex = first_indices;
    
  }
  
  /*
  Rcout << "Pre-loop\n";
  Rcout << "\n Length of new_row_names: " << new_row_names.length() << "\n";
  Rcout << " First entry: " << new_row_names(0) << "\n";
  
  Rcout << "\n Length of new_row_id: " << new_row_id.length() << "\n";
  Rcout << " First entry: " << new_row_id(0) << "\n";
  
  Rcout << "\n Length of new_age: " << new_age.length() << "\n";
  Rcout << " First entry: " << new_age(0) << "\n";
  
  Rcout << "\n Length of new_stage: " << new_stage.length() << "\n";
  Rcout << " First entry: " << new_stage(0) << "\n";
  
  Rcout << "\n Length of new_stage2: " << new_stage2.length() << "\n";
  Rcout << " First entry: " << new_stage2(0) << "\n";
  
  Rcout << "\n Length of new_stage1: " << new_stage1.length() << "\n";
  Rcout << " First entry: " << new_stage1(0) << "\n";
  
  Rcout << "V";
  */
  
  // Main analysis loop
  IntegerVector pop_by_years (years_num);
  
  for (int i = 0; i < new_df_rows; i++) {
    new_row_names(i) = i + 1;
    
    new_row_id(i) = new_year(i);
    new_row_id(i) += " ";
    
    // Here we set the stage and row id designations
    if (check_stage && !check_age) {
      if (!historical) {
        new_stage(i) = new_stage2(i);
        
      } else {
        new_stage(i) = new_stage2(i);
        new_stage(i) += " ";
        new_stage(i) += new_stage1(i);
        
      }
      new_row_id(i) += new_stage(i);
      
    } else if (!check_stage && check_age) {
      new_row_id(i) += new_age(i);
      
    } else if (check_stage && check_age) {
      new_stage(i) = new_stage2(i);
      
      new_row_id(i) += new_stage(i);
      new_row_id(i) += " ";
      new_row_id(i) += new_age(i);
      
    }
    
    // Now we'll find the frequencies of age-stage-year combos
    arma::uvec year_guys = find(data_year2 == new_year(i));
    int year_guys_length = year_guys.n_elem;
    
    if (check_stage && !check_age) {
      if (year_guys_length > 0) {
        for (int j = 0; j < year_guys_length; j++) {
          if (stringcompare_hard(as<std::string>(data_stage2(year_guys(j))), as<std::string>(new_stage2(i)))) {
            if (!historical) {
              arma::uvec year_found = find(all_years == new_year(i));
              pop_by_years(year_found(0))++;
              
              new_frequency(i)++;
            } else {
            
              if (stringcompare_hard(as<std::string>(data_stage1(year_guys(j))), as<std::string>(new_stage1(i)))) {
                arma::uvec year_found = find(all_years == new_year(i));
                pop_by_years(year_found(0))++;
                
                new_frequency(i)++;
              }
            }
          }
        }
      } else {
        arma::uvec year_guys3 = find(data_year2 == new_year(i) - 1);
        int year_guys_length3 = year_guys3.n_elem;
        
        for (int j = 0; j < year_guys_length3; j++) {
          if (stringcompare_hard(as<std::string>(data_stage3(year_guys3(j))), as<std::string>(new_stage2(i)))) {
            if (!historical) {
              arma::uvec year_found = find(all_years == new_year(i));
              pop_by_years(year_found(0))++;
              
              new_frequency(i)++;
            } else {
              if (stringcompare_hard(as<std::string>(data_stage2(year_guys3(j))), as<std::string>(new_stage1(i)))) {
                arma::uvec year_found = find(all_years == new_year(i));
                pop_by_years(year_found(0))++;
                
                new_frequency(i)++;
              }
            }
          }
        }
      }
      
      
      
    } else if (check_stage && check_age) {
      if (year_guys_length > 0) {
        for (int j = 0; j < year_guys_length; j++) {
          if (stringcompare_hard(as<std::string>(data_stage2(year_guys(j))), as<std::string>(new_stage2(i)))) {
            if (data_agecol(year_guys(j)) == new_age(i)) {
              
              arma::uvec year_found = find(all_years == new_year(i));
              pop_by_years(year_found(0))++;
              
              new_frequency(i)++;
            }
          }
        }
      } else {
        arma::uvec year_guys3 = find(data_year2 == new_year(i) - 1);
        int year_guys_length3 = year_guys3.n_elem;
        
        for (int j = 0; j < year_guys_length3; j++) {
          if (stringcompare_hard(as<std::string>(data_stage3(year_guys3(j))), as<std::string>(new_stage2(i)))) {
            if (data_agecol(year_guys3(j)) == new_age(i) - 1) {
              
              arma::uvec year_found = find(all_years == new_year(i));
              pop_by_years(year_found(0))++;
              
              new_frequency(i)++;
            }
          }
        }
      }
    } else if (!check_stage && check_age) {
      if (year_guys_length > 0) {
        for (int j = 0; j < year_guys_length; j++) {
          if (data_agecol(year_guys(j)) == new_age(i)) {
            
            arma::uvec year_found = find(all_years == new_year(i));
            pop_by_years(year_found(0))++;
            
            new_frequency(i)++;
          }
        }
      } else {
        arma::uvec year_guys3 = find(data_year2 == new_year(i) - 1);
        int year_guys_length3 = year_guys3.n_elem;
        
        for (int j = 0; j < year_guys_length3; j++) {
          if (data_agecol(year_guys3(j)) == new_age(i) - 1) {
            
            arma::uvec year_found = find(all_years == new_year(i));
            pop_by_years(year_found(0))++;
            
            new_frequency(i)++;
          }
        }
      }
    }
  }
  
  /*
  Rcout << "Post-loop\n";
  Rcout << "\n Length of new_row_id: " << new_row_id.length() << "\n";
  Rcout << " First entry: " << new_row_id(0) << "\n";
  
  Rcout << "\n Length of new_stage: " << new_stage.length() << "\n";
  Rcout << " First entry: " << new_stage(0) << "\n";
  
  Rcout << "\n Length of new_age: " << new_age.length() << "\n";
  Rcout << " First entry: " << new_age(0) << "\n";
  
  Rcout << "\n Length of new_year: " << new_year.length() << "\n";
  Rcout << " First entry: " << new_year(0) << "\n";
  
  Rcout << "\n Length of new_frequency: " << new_frequency.length() << "\n";
  Rcout << " First entry: " << new_frequency(0) << "\n";
  
  Rcout << "\n Length of new_actualprop: " << new_actualprop.length() << "\n";
  Rcout << " First entry: " << new_actualprop(0) << "\n";
  */
  
  for (int i = 0; i < new_df_rows; i++) {
    arma::uvec year_found = find(all_years == new_year(i));
    
    if (pop_by_years(year_found(0)) > 0) {
      new_actualprop(i) = new_frequency(i) / static_cast<double>(pop_by_years(year_found(0)));
    } else {
      new_actualprop(i) = 0;
    }
  }
  
  
  // Structure the output data frame
  List output;
  
  if (check_stage && !check_age) {
    List raw_output (8);
    
    raw_output(0) = new_row_id;
    raw_output(1) = new_stageindex;
    raw_output(2) = new_stage;
    raw_output(3) = new_stage2;
    raw_output(4) = new_stage1;
    raw_output(5) = new_year;
    raw_output(6) = new_frequency;
    raw_output(7) = new_actualprop;
    
    CharacterVector output_names = {"rowid", "stageindex", "stage", "stage2",
      "stage1", "year2", "Freq", "actual_prop"};
    raw_output.attr("names") = output_names;
    output = raw_output;
    
  } else if (!check_stage && check_age) {
    List raw_output (5);
    
    raw_output(0) = new_row_id;
    raw_output(1) = new_age;
    raw_output(2) = new_year;
    raw_output(3) = new_frequency;
    raw_output(4) = new_actualprop;
    
    CharacterVector output_names = {"rowid", "age", "year2", "Freq",
      "actual_prop"};
    raw_output.attr("names") = output_names;
    output = raw_output;
    
  } else if (check_stage && check_age) {
    List raw_output (9);
    
    raw_output(0) = new_row_id;
    raw_output(1) = new_stageindex;
    raw_output(2) = new_stage;
    raw_output(3) = new_stage2;
    raw_output(4) = new_stage1;
    raw_output(5) = new_age;
    raw_output(6) = new_year;
    raw_output(7) = new_frequency;
    raw_output(8) = new_actualprop;
    
    CharacterVector output_names = {"rowid", "stageindex", "stage", "stage2",
      "stage1", "age", "year2", "Freq", "actual_prop"};
    raw_output.attr("names") = output_names;
    output = raw_output;
  
  }
  
  output.attr("row.names") = new_row_names;
  output.attr("class") = "data.frame";
  return output;
}

//' Check and Reorganize Density Input Table Into Usable Format
//' 
//' Function \code{density_reassess()} takes a density input table as supplied
//' by the \code{\link{density_input}()} function, and checks and rearranges it
//' into a single, complete density input table.
//' 
//' @name density_reassess
//' 
//' @param stageframe The correct stageframe, already modified by
//' \code{\link{.sf_reassess}()}.
//' @param dens_inp The density input data frame as is toward the end of
//' \code{\link{density_input}()}.
//' @param agestages The agestages element from the used \code{lefkoMat} object.
//' Only needed if an age-by-stage MPM will be used.
//' @param historical A logical value denoting whether MPM is historical.
//' Defaults to \code{FALSE}.
//' @param agebystage A logical value denoting whether MPM is age-by-stage.
//' Defaults to \code{FALSE}.
//' 
//' @return A corrected density input deta frame, usable in density-dependent
//' MPM creation.
//' 
//' @keywords internal
//' @noRd
Rcpp::DataFrame density_reassess(DataFrame stageframe, DataFrame dens_inp,
  Nullable<DataFrame> agestages, bool historical = false,
  bool agebystage = false) {
  
  StringVector stagevec = as<StringVector>(stageframe["stage"]);
  arma::ivec stageidvec = as<arma::ivec>(stageframe["stage_id"]);
  NumericVector minagevec = as<NumericVector>(stageframe["min_age"]);
  NumericVector maxagevec = as<NumericVector>(stageframe["max_age"]);
  arma::uvec repvec = as<arma::uvec>(stageframe["repstatus"]);
  arma::uvec obsvec = as<arma::uvec>(stageframe["obsstatus"]);
  arma::uvec propvec = as<arma::uvec>(stageframe["propstatus"]);
  arma::uvec immvec = as<arma::uvec>(stageframe["immstatus"]);
  arma::uvec matvec = as<arma::uvec>(stageframe["matstatus"]);
  arma::uvec indvec = as<arma::uvec>(stageframe["indataset"]);
  arma::ivec groupvec = as<arma::ivec>(stageframe["group"]);
  
  int no_stages = indvec.n_elem;
  arma::ivec alive(no_stages, fill::ones);

  // Identify all groups in the stageframe
  arma::ivec all_groups = unique(groupvec);
  int no_groups = all_groups.n_elem;
  StringVector group_text(no_groups);
  
  for (int i = 0; i < no_groups; i++) {
    group_text(i) = "group";
    group_text(i) += std::to_string(all_groups(i));
  }
  
  StringVector stage3_di = dens_inp["stage3"];
  StringVector stage2_di = dens_inp["stage2"];
  StringVector stage1_di = dens_inp["stage1"];
  IntegerVector age2_di = dens_inp["age2"];
  IntegerVector style_di = dens_inp["style"];
  IntegerVector time_delay_di = dens_inp["time_delay"];
  NumericVector alpha_di = dens_inp["alpha"];
  NumericVector beta_di = dens_inp["beta"];
  IntegerVector type_di = dens_inp["type"];
  IntegerVector type_t12_di = dens_inp["type_t12"];
  int di_rows = stage3_di.length();
  
  StringVector unique_stages = unique(stagevec);
  StringVector extra_terms = {"rep", "nrep", "immat", "mat", "prop", "npr", "all", "obs", "nobs"};
  
  int no_newstages = unique_stages.length();
  int no_extraterms = extra_terms.length();
  
  StringVector all_possible_stage_terms(no_newstages + no_extraterms + no_groups);
  for (int i = 0; i < no_newstages; i++) {
    all_possible_stage_terms(i) = unique_stages(i);
  }
  for (int i = 0; i < no_extraterms; i++) {
    all_possible_stage_terms(i + no_newstages) = extra_terms(i);
  }
  for (int i = 0; i < no_groups; i++) {
    all_possible_stage_terms(i + no_newstages + no_extraterms) = group_text(i);
  }
  
  DataFrame agestages_;
  arma::ivec agestages_stageid;
  StringVector agestages_stage;
  arma::ivec agestages_age;
  int agestages_rows = 0;
  int age_min = 0;
  
  if (agestages.isNotNull()) {
    if (is<NumericVector>(agestages)) {
      if (agebystage) {
        Rf_warningcall(R_NilValue, "Function density_input() requires an agestages object if input MPM is age-by-stage.");
        agebystage = false;
      }
      
    } else if (is<DataFrame>(agestages)) {
      agestages_ = as<DataFrame>(agestages);
      
      if (agestages_.length() == 1) {
        agebystage = false;
        
      } else {
        agestages_stage = as<StringVector>(agestages_["stage"]);
        agestages_stageid = as<arma::ivec>(agestages_["stage_id"]);
        agestages_age = as<arma::ivec>(agestages_["age"]);
        agestages_rows = agestages_stage.length();
        age_min = agestages_age.min();
        
        agebystage = true;
      }
    } else {
      throw Rcpp::exception("Object input as agestages is not recognized.", false);
    }
  }
  
  if (historical && agebystage) {
    throw Rcpp::exception("MPMs cannot be both historical and age-by-stage.", false);
  }
  
  // Check for good entries density input data frame
  for (int i = 0; i < stage3_di.length(); i++) {
    int s3di_count {0};
    int s2di_count {0};
    int s1di_count {0};
    
    for (int j = 0; j < all_possible_stage_terms.length(); j++) {
      std::string s3used = as<std::string>(stage3_di(i));
      std::string s2used = as<std::string>(stage2_di(i));
      std::string s1used = as<std::string>(stage1_di(i));
      
      for (int k = 0; k < s3used.size(); k++) {
        s3used[k] = tolower(s3used[k]);
      }
      for (int k = 0; k < s3used.size(); k++) {
        s2used[k] = tolower(s2used[k]);
      }
      for (int k = 0; k < s3used.size(); k++) {
        s1used[k] = tolower(s1used[k]);
      }
      
      if (stringcompare_hard(s3used, as<std::string>(all_possible_stage_terms(j)))) s3di_count++;
      if (stringcompare_hard(as<std::string>(stage3_di(i)), as<std::string>(all_possible_stage_terms(j)))) s3di_count++;
      
      if (stringcompare_hard(s2used, as<std::string>(all_possible_stage_terms(j)))) s2di_count++;
      if (stringcompare_hard(as<std::string>(stage2_di(i)), as<std::string>(all_possible_stage_terms(j)))) s2di_count++;
      
      if (stringcompare_hard(s3used, "notalive")) {
        throw Rcpp::exception("Stage NotAlive is not allowed.", false);
      }
      if (stringcompare_hard(s2used, "notalive")) {
        throw Rcpp::exception("Stage NotAlive is not allowed.", false);
      }
      
      if (historical) {
        if (stringcompare_hard(s1used, as<std::string>(all_possible_stage_terms(j)))) s1di_count++;
        if (stringcompare_hard(as<std::string>(stage1_di(i)), as<std::string>(all_possible_stage_terms(j)))) s1di_count++;
        
        if (stringcompare_hard(s1used, "notalive")) {
          throw Rcpp::exception("Stage NotAlive is not allowed.", false);
        }
      } 
    }
    
    if (s3di_count == 0) {
      throw Rcpp::exception("Stage names in density input frame (stage3) must match stageframe",
        false);
    }
    if (s2di_count == 0) {
      throw Rcpp::exception("Stage names in density input frame (stage2) must match stageframe",
        false);
    }
    if (historical) {
      if (s1di_count == 0) {
        throw Rcpp::exception("Stage names in density input frame (stage1) must match stageframe",
          false);
      }
    }
  }
    
  IntegerVector s1_calls (di_rows, 1);
  IntegerVector s2_calls (di_rows, 1);
  IntegerVector s3_calls (di_rows, 1);
  IntegerVector s3_planned (di_rows, 1);
  IntegerVector s2_planned (di_rows, 1);
  IntegerVector s1_planned (di_rows, 1);
  
  IntegerVector s123_calls (di_rows, 1);
  
  List age3_calls (di_rows);
  List age2_calls (di_rows);
  
  List stageid3_calls (di_rows);
  List stageid2_calls (di_rows);
  
  arma::uvec prop_stages = find(propvec);
  arma::uvec prop0_stages = find(propvec == 0);
  arma::uvec imm_stages = find(immvec);
  arma::uvec alive_stages = find(alive);
  arma::uvec mat_stages = find(matvec);
  arma::uvec rep_stages = find(repvec);
  arma::uvec rep0_stages = find(repvec == 0);
  arma::uvec mat_rep0_stages = intersect(mat_stages, rep0_stages);
  arma::uvec obs_stages = find(obsvec);
  arma::uvec obs0_stages = find(obsvec == 0);
  arma::uvec all_stages = find(alive); // 7 "all"
  int no_current_group {0};
  
  // Now we build the expanded and edited density input frame
  for (int i = 0; i < di_rows; i++) {
    
    std::string s3used = as<std::string>(stage3_di(i));
    std::string s2used = as<std::string>(stage2_di(i));
    std::string s1used = as<std::string>(stage1_di(i));
    
    for (int j = 0; j < s3used.size(); j++) {
      s3used[j] = tolower(s3used[j]);
    }
    for (int j = 0; j < s2used.size(); j++) {
      s2used[j] = tolower(s2used[j]);
    }
    for (int j = 0; j < s1used.size(); j++) {
      s1used[j] = tolower(s1used[j]);
    }
    
    // Time t+1
    if (stringcompare_hard(s3used, "prop")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < prop_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (prop_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (prop_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < prop_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (prop_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (prop_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
        s3_calls(i) = prop_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s3used, "npr")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < prop0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (prop0_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (prop0_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < prop0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (prop0_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (prop0_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
         s3_calls(i) = prop0_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s3used, "immat")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < imm_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (imm_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (imm_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < imm_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (imm_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (imm_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
        s3_calls(i) = imm_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s3used, "mat")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < mat_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (mat_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (mat_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < mat_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (mat_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (mat_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
         s3_calls(i) = mat_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s3used, "rep")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < rep_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (rep_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (rep_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < rep_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (rep_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (rep_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
        s3_calls(i) = rep_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s3used, "nrep")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < mat_rep0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (mat_rep0_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (mat_rep0_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < mat_rep0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (mat_rep0_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (mat_rep0_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
        s3_calls(i) = mat_rep0_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s3used, "obs")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < obs_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (obs_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (obs_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < obs_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (obs_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (obs_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
        s3_calls(i) = obs_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s3used, "nobs")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < obs0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (obs0_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (obs0_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < obs0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (obs0_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (obs0_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
        s3_calls(i) = obs0_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s3used, "all")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < all_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (all_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  found_stages++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  found_stages++;
                }
              }
            } else {
              if (all_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) found_stages++;
              }
            }
          }
        }
        
        s3_calls(i) = found_stages;
        IntegerVector age3_vec(found_stages);
        IntegerVector stageid3_vec(found_stages);
        int a3_counter = 0;
        
        for (int j = 0; j < all_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (type_di(i) == 1) {
              if (all_stages(j) == agestages_stageid(k) - 1) {
                if (IntegerVector::is_na(age2_di(i))) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                } else if (agestages_age(k) == (age2_di(i) + 1)) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            } else {
              if (all_stages(j) == agestages_stageid(k) - 1) {
                if (agestages_age(k) == age_min) {
                  age3_vec(a3_counter) = agestages_age(k);
                  stageid3_vec(a3_counter) = agestages_stageid(k);
                  a3_counter++;
                }
              }
            }
          }
        }
        age3_calls(i) = age3_vec;
        stageid3_calls(i) = stageid3_vec;
        
      } else {
        s3_calls(i) = all_stages.n_elem;
      }
      
    } else {
      for (int j = 0; j < no_groups; j++) {
        if (stage3_di(i) == group_text(j)) {
          arma::uvec current_group = find(groupvec == j);
          no_current_group = current_group.n_elem;
          
          if (agebystage) {
            int found_stages = 0;
            
            for (int j = 0; j < current_group.n_elem; j++) {
              for (int k = 0; k < agestages_rows; k++) {
                if (type_di(i) == 1) {
                  if (current_group(j) == agestages_stageid(k) - 1) {
                    if (IntegerVector::is_na(age2_di(i))) {
                      found_stages++;
                    } else if (agestages_age(k) == (age2_di(i) + 1)) {
                      found_stages++;
                    }
                  }
                } else {
                  if (current_group(j) == agestages_stageid(k) - 1) {
                    if (agestages_age(k) == age_min) found_stages++;
                  }
                }
              }
            }
            
            s3_calls(i) = found_stages;
            IntegerVector age3_vec(found_stages);
            IntegerVector stageid3_vec(found_stages);
            int a3_counter = 0;
            
            for (int j = 0; j < current_group.n_elem; j++) {
              for (int k = 0; k < agestages_rows; k++) {
                if (type_di(i) == 1) {
                  if (current_group(j) == agestages_stageid(k) - 1) {
                    if (IntegerVector::is_na(age2_di(i))) {
                      age3_vec(a3_counter) = agestages_age(k);
                      stageid3_vec(a3_counter) = agestages_stageid(k);
                      a3_counter++;
                    } else if (agestages_age(k) == (age2_di(i) + 1)) {
                      age3_vec(a3_counter) = agestages_age(k);
                      stageid3_vec(a3_counter) = agestages_stageid(k);
                      a3_counter++;
                    }
                  }
                } else {
                  if (current_group(j) == agestages_stageid(k) - 1) {
                    if (agestages_age(k) == age_min) {
                      age3_vec(a3_counter) = agestages_age(k);
                      stageid3_vec(a3_counter) = agestages_stageid(k);
                      a3_counter++;
                    }
                  }
                }
              }
            }
            age3_calls(i) = age3_vec;
            stageid3_calls(i) = stageid3_vec;
            
          } else {
            s3_calls(i) = no_current_group;
          }
        }
      }
    }
    if (s3_calls(i) == 0) s3_calls(i) = 1;
    
    // Time t
    if (stringcompare_hard(s2used, "prop")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < prop_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (prop_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < prop_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (prop_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = prop_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s2used, "npr")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < prop0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (prop0_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < prop0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (prop0_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = prop0_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s2used, "immat")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < imm_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (imm_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < imm_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (imm_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = imm_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s2used, "mat")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < mat_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (mat_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < mat_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (mat_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = mat_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s2used, "rep")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < rep_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (rep_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < rep_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (rep_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = rep_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s2used, "nrep")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < mat_rep0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (mat_rep0_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < mat_rep0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (mat_rep0_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = mat_rep0_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s2used, "obs")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < obs_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (obs_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < obs_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (obs_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = obs_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s2used, "nobs")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < obs0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (obs0_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < obs0_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (obs0_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = obs0_stages.n_elem;
      }
      
    } else if (stringcompare_hard(s2used, "all")) {
      if (agebystage) {
        int found_stages = 0;
        
        for (int j = 0; j < all_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (all_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                found_stages++;
              } else if (agestages_age(k) == age2_di(i)) {
                found_stages++;
              }
            }
          }
        }
        
        s2_calls(i) = found_stages;
        IntegerVector age2_vec(found_stages);
        IntegerVector stageid2_vec(found_stages);
        int a2_counter = 0;
        
        for (int j = 0; j < all_stages.n_elem; j++) {
          for (int k = 0; k < agestages_rows; k++) {
            if (all_stages(j) == agestages_stageid(k) - 1) {
              if (IntegerVector::is_na(age2_di(i))) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              } else if (agestages_age(k) == age2_di(i)) {
                age2_vec(a2_counter) = agestages_age(k);
                stageid2_vec(a2_counter) = agestages_stageid(k);
                a2_counter++;
              }
            }
          }
        }
        age2_calls(i) = age2_vec;
        stageid2_calls(i) = stageid2_vec;
        
      } else {
        s2_calls(i) = all_stages.n_elem;
      }
    } else {
      for (int j = 0; j < no_groups; j++) {
        if (stage2_di(i) == group_text(j)) {
          arma::uvec current_group = find(groupvec == j);
          no_current_group = current_group.n_elem;
          
          if (agebystage) {
            int found_stages = 0;
            
            for (int j = 0; j < current_group.n_elem; j++) {
              for (int k = 0; k < agestages_rows; k++) {
                if (current_group(j) == agestages_stageid(k) - 1) {
                  if (IntegerVector::is_na(age2_di(i))) {
                    found_stages++;
                  } else if (agestages_age(k) == age2_di(i)) {
                    found_stages++;
                  }
                }
              }
            }
            
            s2_calls(i) = found_stages;
            IntegerVector age2_vec(found_stages);
            IntegerVector stageid2_vec(found_stages);
            int a2_counter = 0;
            
            for (int j = 0; j < current_group.n_elem; j++) {
              for (int k = 0; k < agestages_rows; k++) {
                if (current_group(j) == agestages_stageid(k) - 1) {
                  if (IntegerVector::is_na(age2_di(i))) {
                    age2_vec(a2_counter) = agestages_age(k);
                    stageid2_vec(a2_counter) = agestages_stageid(k);
                    a2_counter++;
                  } else if (agestages_age(k) == age2_di(i)) {
                    age2_vec(a2_counter) = agestages_age(k);
                    stageid2_vec(a2_counter) = agestages_stageid(k);
                    a2_counter++;
                  }
                }
              }
            }
            age2_calls(i) = age2_vec;
            stageid2_calls(i) = stageid2_vec;
            
          } else {
            s2_calls(i) = no_current_group;
          }
        }
      }
    }
    if (s2_calls(i) == 0) s2_calls(i) = 1;
    
    // Time t-1
    if (stringcompare_hard(s1used, "prop")) {
      s1_calls(i) = prop_stages.n_elem;
    } else if (stringcompare_hard(s1used, "npr")) {
      s1_calls(i) = prop0_stages.n_elem;
    } else if (stringcompare_hard(s1used, "immat")) {
      s1_calls(i) = imm_stages.n_elem;
    } else if (stringcompare_hard(s1used, "mat")) {
      s1_calls(i) = mat_stages.n_elem;
    } else if (stringcompare_hard(s1used, "rep")) {
      s1_calls(i) = rep_stages.n_elem;
    } else if (stringcompare_hard(s1used, "nrep")) {
      s1_calls(i) = mat_rep0_stages.n_elem;
    } else if (stringcompare_hard(s1used, "obs")) {
      s1_calls(i) = obs_stages.n_elem;
    } else if (stringcompare_hard(s1used, "nobs")) {
      s1_calls(i) = obs0_stages.n_elem;
    } else if (stringcompare_hard(s1used, "all")) {
      s1_calls(i) = all_stages.n_elem;
    } else if (StringVector::is_na(stage1_di(i))) {
      s1_calls(i) = 1;
    } else {
      for (int j = 0; j < no_groups; j++) {
        if (stage1_di(i) == group_text(j)) {
          arma::uvec current_group = find(groupvec == j);
          no_current_group = current_group.n_elem;
          
          s1_calls(i) = no_current_group;
        }
      }
    }
    if (s1_calls(i) == 0) s1_calls(i) = 1;
    
    s123_calls(i) = s3_calls(i) * s2_calls(i) * s1_calls(i);
  }
  
  NumericVector basepoints(di_rows, 0.0);
  for (int i = 0; i < (di_rows - 1); i++) {
    basepoints(i+1) = basepoints(i) + s123_calls(i);
  }
  
  // New output data frame set-up
  int newdi_rows = sum(s123_calls);
  
  StringVector stage3_newdi(newdi_rows);
  StringVector stage2_newdi(newdi_rows);
  StringVector stage1_newdi(newdi_rows);
  IntegerVector age2_newdi(newdi_rows);
  IntegerVector style_newdi(newdi_rows);
  NumericVector alpha_newdi(newdi_rows);
  NumericVector beta_newdi(newdi_rows);
  IntegerVector time_delay_newdi(newdi_rows);
  IntegerVector type_newdi(newdi_rows);
  IntegerVector type_t12_newdi(newdi_rows);
  
  int overall_counter {0};
  // int group_check {0};
  
  int group_baseline3 {0};
  int group_baseline2 {0};
  int group_baseline1 {0};
  
  int group_ratchet3 {0};
  int group_ratchet2 {0};
  int group_ratchet1 {0};
  
  int prevl3 {0};
  int prevl2 {0};
  int prevl1 {0};
  
  for (int i = 0; i < di_rows; i++) {
    overall_counter = 0;
    
    int age3_counter = 0;
    int age2_counter = 0;
    
    std::string s3used = as<std::string>(stage3_di(i));
    std::string s2used = as<std::string>(stage2_di(i));
    std::string s1used = as<std::string>(stage1_di(i));
    
    for (int j = 0; j < s3used.size(); j++) {
      s3used[j] = tolower(s3used[j]);
    }
    for (int j = 0; j < s2used.size(); j++) {
      s2used[j] = tolower(s2used[j]);
    }
    for (int j = 0; j < s1used.size(); j++) {
      s1used[j] = tolower(s1used[j]);
    }
    
    for (int j = 0; j < s1_calls(i); j++) {
      for (int k = 0; k < s2_calls(i); k++) {
        for (int l = 0; l < s3_calls(i); l++) {
          
          // Time t+1
          if (stringcompare_hard(s3used, "prop")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(prop_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "npr")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(prop0_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "immat")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(imm_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "mat")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(mat_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "rep")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(rep_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "nrep")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(mat_rep0_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "obs")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(obs_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "nobs")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(obs0_stages(l));
            }
            
          } else if (stringcompare_hard(s3used, "all")) {
            if (agebystage) {
              IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
              IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
              int a3v_length = age3_vec.length();
              
              stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
              
              age3_counter++;
              if (age3_counter == a3v_length) age3_counter = 0;
            
            } else {
              stage3_newdi(basepoints(i) + overall_counter) = stagevec(all_stages(l));
            }
            
          } else {
            int group_check = 0;
            
            for (int m = 0; m < no_groups; m++) {
              if (stage3_di(i) == group_text(m)) {
                if (agebystage) {
                  group_check = 1;
                  
                  IntegerVector age3_vec = as<IntegerVector>(age3_calls(i));
                  IntegerVector stageid3_vec = as<IntegerVector>(stageid3_calls(i));
                  int a3v_length = age3_vec.length();
                  
                  stage3_newdi(basepoints(i) + overall_counter) = stagevec((stageid3_vec(age3_counter) - 1));
                  
                  age3_counter++;
                  if (age3_counter == a3v_length) age3_counter = 0;
                
                } else {
                  if (l == 0) group_ratchet3 = 0;
                  if (l != prevl3 && l != 0) group_ratchet3 += 1;
                  
                  group_check = 1;
                  arma::uvec current_group = find(groupvec == m);
                  int current_group_length = current_group.n_elem;
                  if (group_ratchet3 > (current_group_length - 1)) {
                    group_ratchet3 = 0;
                  }
                  
                  if (group_ratchet3 == 0) {
                    group_baseline3 = l;
                  }
                  
                  stage3_newdi(basepoints(i) + overall_counter) = 
                    stagevec(current_group(l - group_baseline3));
                  
                  prevl3 = l;
                  
                }
              }
            }
            if (group_check == 0) {
              stage3_newdi(basepoints(i) + overall_counter) = stage3_di(i);
            }
            
            group_check = 0;
          }
          
          // Set up of most core variables in output data frame
          age2_newdi(basepoints(i) + overall_counter) = age2_di(i);
          style_newdi(basepoints(i) + overall_counter) = style_di(i);
          alpha_newdi(basepoints(i) + overall_counter) = alpha_di(i);
          beta_newdi(basepoints(i) + overall_counter) = beta_di(i);
          time_delay_newdi(basepoints(i) + overall_counter) = time_delay_di(i);
          type_newdi(basepoints(i) + overall_counter) = type_di(i);
          type_t12_newdi(basepoints(i) + overall_counter) = type_t12_di(i);
          
          // Time t
          if (stringcompare_hard(s2used, "prop")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(prop_stages(k));
            }
          } else if (stringcompare_hard(s2used, "npr")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(prop0_stages(k));
            }
          } else if (stringcompare_hard(s2used, "immat")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(imm_stages(k));
            }
          } else if (stringcompare_hard(s2used, "mat")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(mat_stages(k));
            }
          } else if (stringcompare_hard(s2used, "rep")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(rep_stages(k));
            }
          } else if (stringcompare_hard(s2used, "nrep")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(mat_rep0_stages(k));
            }
          } else if (stringcompare_hard(s2used, "obs")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(obs_stages(k));
            }
          } else if (stringcompare_hard(s2used, "nobs")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(obs0_stages(k));
            }
          } else if (stringcompare_hard(s2used, "all")) {
            if (agebystage) {
              IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
              IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
              int a2v_length = age2_vec.length();
              
              stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
              age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
              
              age2_counter++;
              if (age2_counter == a2v_length) age2_counter = 0;
            
            } else {
              stage2_newdi(basepoints(i) + overall_counter) = stagevec(all_stages(k));
            }
          } else {
            int group_check = 0;
            
            for (int m = 0; m < no_groups; m++) {
              if (stage2_di(i) == group_text(m)) {
                if (agebystage) {
                  group_check = 1;
                  
                  IntegerVector age2_vec = as<IntegerVector>(age2_calls(i));
                  IntegerVector stageid2_vec = as<IntegerVector>(stageid2_calls(i));
                  int a2v_length = age2_vec.length();
                  
                  stage2_newdi(basepoints(i) + overall_counter) = stagevec((stageid2_vec(age2_counter) - 1));
                  age2_newdi(basepoints(i) + overall_counter) = age2_vec(age2_counter);
                  
                  age2_counter++;
                  if (age2_counter == a2v_length) age2_counter = 0;
                
                } else {
                  if (k == 0) group_ratchet2 = 0;
                  if (k != prevl2 && k != 0) group_ratchet2 += 1;
                  
                  group_check = 1;
                  arma::uvec current_group = find(groupvec == m);
                  int current_group_length = current_group.n_elem;
                  if (group_ratchet2 > (current_group_length - 1)) {
                    group_ratchet2 = 0;
                  }
                  
                  if (group_ratchet2 == 0) {
                    group_baseline2 = k;
                  }
                  
                  stage2_newdi(basepoints(i) + overall_counter) =
                    stagevec(current_group(k - group_baseline2));
                  
                  prevl2 = k;
                }
              }
            }
             
            if (group_check == 0) {
             stage2_newdi(basepoints(i) + overall_counter) = stage2_di(i);
            }
            
            group_check = 0;
          }
          
          // Time t-1
          if (stringcompare_hard(s1used, "prop")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(prop_stages(j));
          } else if (stringcompare_hard(s1used, "npr")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(prop0_stages(j));
          } else if (stringcompare_hard(s1used, "immat")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(imm_stages(j));
          } else if (stringcompare_hard(s1used, "mat")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(mat_stages(j));
          } else if (stringcompare_hard(s1used, "rep")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(rep_stages(j));
          } else if (stringcompare_hard(s1used, "nrep")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(mat_rep0_stages(j));
          } else if (stringcompare_hard(s1used, "obs")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(obs_stages(j));
          } else if (stringcompare_hard(s1used, "nobs")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(obs0_stages(j));
          } else if (stringcompare_hard(s1used, "all")) {
            stage1_newdi(basepoints(i) + overall_counter) = stagevec(all_stages(j));
          } else {
            int group_check = 0;
            
            for (int m = 0; m < no_groups; m++) {
              if (stage1_di(i) == group_text(m)) {
                if (j == 0) group_ratchet1 = 0;
                if (j != prevl1 && j != 0) group_ratchet1 += 1;
                
                group_check = 1;
                arma::uvec current_group = find(groupvec == m);
                int current_group_length = current_group.n_elem;
                if (group_ratchet1 > (current_group_length - 1)) {
                  group_ratchet1 = 0;
                }
                
                if (group_ratchet1 == 0) {
                  group_baseline1 = j;
                }
                
                stage1_newdi(basepoints(i) + overall_counter) =
                  stagevec(current_group(j - group_baseline1));
                
                prevl1 = j;
              }
            }
            
            if (group_check == 0) {
              stage1_newdi(basepoints(i) + overall_counter) = stage1_di(i);
            }
            
            group_check = 0;
          }
          
          overall_counter++;
        }
      }
    }
  }
  
  // Output final set-up
  Rcpp::List new_di(10);
  
  new_di(0) = stage3_newdi;
  new_di(1) = stage2_newdi;
  new_di(2) = stage1_newdi;
  new_di(3) = age2_newdi;
  new_di(4) = style_newdi;
  new_di(5) = alpha_newdi;
  new_di(6) = beta_newdi;
  new_di(7) = time_delay_newdi;
  new_di(8) = type_newdi;
  new_di(9) = type_t12_newdi;
  
  CharacterVector namevec = {"stage3", "stage2", "stage1", "age2", "style",
    "alpha", "beta", "time_delay", "type", "type_t12"};
  CharacterVector newclass = {"data.frame"};
  new_di.attr("names") = namevec;
  new_di.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, newdi_rows);
  new_di.attr("class") = newclass;
  
  return new_di;
}

//' Create a Data Frame of Density Dependence Relationships in Matrix Elements
//' 
//' Function \code{density_input()} provides all necessary data to incorporate
//' density dependence into a \code{lefkoMat} object, a list of matrices, or a
//' single matrix. Four forms of density dependence are allowed, including the
//' Ricker function, the Beverton-Holt function, the Usher function, and the
//' logistic function. In each case, density must have an effect with at least a
//' one time-step delay (see Notes). The resulting data frame provides a guide
//' for other \code{lefko3} functions to modify matrix elements by density.
//'
//' @name density_input
//' 
//' @param mpm The \code{lefkoMat} object that will be subject to density
//' dependent projection.
//' @param stage3 A vector showing the name or number of the stage in occasion
//' \emph{t}+1 in the transitions to be affected by density. Abbreviations for
//' groups of stages are also usable (see Notes).
//' @param stage2 A vector showing the name or number of the stage in occasion
//' \emph{t} in the transition to be affected by density. Abbreviations for
//' groups of stages are also usable (see Notes).
//' @param stage1 A vector showing the name or number of the stage in occasion
//' \emph{t}-1 in the transition to be affected by density. Only needed if a
//' historical MPM is used. Abbreviations for groups of stages are also usable
//' (see Notes).
//' @param age2 A vector showing the age of the stage in occasion \emph{t} in the
//' transition to be affected by density. Only needed if an age-by-stage MPM is
//' used.
//' @param style A vector coding for the style of density dependence on each
//' transition subject to density dependence. Options include \code{1},
//' \code{ricker}, \code{ric}, or \code{r} for the Ricker function; \code{2},
//' \code{beverton}, \code{bev}, and \code{b} for the Beverton-Holt function;
//' \code{3}, \code{usher}, \code{ush}, and \code{u} for the Usher function; and
//' \code{4}, \code{logistic}, \code{log}, and \code{l} for the logistic
//' function. If only a single code is provided, then all noted transitions are
//' assumed to be subject to this style of density dependence. Defaults to
//' \code{ricker}.
//' @param time_delay An integer vector indicating the number of occasions back
//' on which density dependence operates. Defaults to \code{1}, and may not equal
//' any integer less than 1. If a single number is input, then all noted
//' transitions are assumed to be subject to this time delay.  Defaults to
//' \code{1}.
//' @param alpha A vector indicating the numeric values to use as the
//' alpha term in the two parameter Ricker, Beverton-Holt, or Usher function, or
//' the value of the carrying capacity \emph{K} to use in the logistic equation
//' (see \code{Notes} section for more on this term). If a single number is
//' provided, then all noted transitions are assumed to be subject to this value
//' of alpha. Defaults to \code{1}.
//' @param beta A vector indicating the numeric values to use as the beta term in
//' the two parameter Ricker, Beverton-Holt, or Usher function. Used to indicate
//' whether to use \emph{K} as a hard limit in the logistic equation (see section
//' \code{Notes} below). If a single number is provided, then all noted
//' transitions are assumed to be subject to this value of \code{beta}. Defaults
//' to \code{1}.
//' @param type A vector denoting the kind of transition between occasions
//' \emph{t} and \emph{t}+1 to be replaced. This should be entered as \code{1},
//' \code{S}, or \code{s} for the replacement of a survival transition; or 
//' \code{2}, \code{F}, or \code{f} for the replacement of a fecundity
//' transition. If empty or not provided, then defaults to \code{1} for survival
//' transition.
//' @param type_t12 An optional vector denoting the kind of transition between
//' occasions \emph{t}-1 and \emph{t}. Only necessary if a historical MPM in
//' deVries format is desired. This should be entered as \code{1}, \code{S}, or
//' \code{s} for a survival transition; or \code{2}, \code{F}, or \code{f} for a
//' fecundity transitions. Defaults to \code{1} for survival transition, with
//' impacts only on the construction of deVries-format hMPMs.
//' 
//' @return A data frame of class \code{lefkoDens}. This object can be used as
//' input in function \code{\link{projection3}()}.
//' 
//' Variables in this object include the following:
//' \item{stage3}{Stage at occasion \emph{t}+1 in the transition to be replaced.}
//' \item{stage2}{Stage at occasion \emph{t} in the transition to be replaced.}
//' \item{stage1}{Stage at occasion \emph{t}-1 in the transition to be replaced,
//' if applicable.}
//' \item{age2}{Age at occasion \emph{t} in the transition to be replaced, if
//' applicable.}
//' \item{style}{Style of density dependence, coded as 1, 2, 3, or 4 for the
//' Ricker, Beverton-Holt, Usher, or logistic function, respectively.}
//' \item{time_delay}{The time delay on density dependence, in time steps.}
//' \item{alpha}{The value of alpha in the Ricker, Beverton-Holt, or Usher
//' function, or the value of carrying capacity, \emph{K}, in the logistic
//' function.}
//' \item{beta}{The value of beta in the Ricker, Beverton-Holt, or Usher
//' function.}
//' \item{type}{Designates whether the transition from occasion \emph{t} to
//' occasion \emph{t}+1 is a survival transition probability (1), or a fecundity
//' rate (2).}
//' \item{type_t12}{Designates whether the transition from occasion \emph{t}-1 to
//' occasion \emph{t} is a survival transition probability (1), a fecundity rate
//' (2).}
//' 
//' @section Notes:
//' This function provides inputs when density dependence is operationalized
//' directly on matrix elements. It can be used in both \code{projection3()} and
//' \code{f_projection3()}. Users wishing to modify vital rate functions by
//' density dependence functions for use in function-based projections with
//' function \code{f_projection3()} should use function \code{density_vr()} to
//' provide the correct inputs.
//' 
//' The parameters \code{alpha} and \code{beta} are applied according to the
//' two-parameter Ricker function, the two-parameter Beverton-Holt function, the
//' two-parameter Usher function, or the one-parameter logistic function.
//' Although the default is that a 1 time step delay is assumed, greater time
//' delays can be set through the \code{time_delay} option.
//' 
//' Entries in \code{stage3}, \code{stage2}, and \code{stage1} can include
//' abbreviations for groups of stages. Use \code{rep} if all reproductive stages
//' are to be used, \code{nrep} if all mature but non-reproductive stages are to
//' be used, \code{mat} if all mature stages are to be used, \code{immat} if all
//' immature stages are to be used, \code{prop} if all propagule stages are to be
//' used, \code{npr} if all non-propagule stages are to be used, \code{obs} if
//' all observable stages are to be used, \code{nobs} if all unobservable stages
//' are to be used, and leave empty or use \code{all} if all stages in stageframe
//' are to be used.
//' 
//' When using the logistic function, it is possible that the time delay used in
//' density dependent simulations will cause matrix elements to become negative.
//' To prevent this behavior, set the associated \code{beta} term to \code{1.0}.
//' Doing so will set \code{K} as the hard limit in the logistic equation,
//' essentially setting a minimum limit at \code{0} for all matrix elements
//' modified.
//' 
//' @seealso \code{\link{start_input}()}
//' @seealso \code{\link{projection3}()}
//' 
//' @examples
//' \donttest{
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
//' 
//' e3d <- density_input(ehrlen3mean, stage3 = c("Sd", "Sdl"),
//'   stage2 = c("rep", "rep"), stage1 = c("all", "all"), style = 1,
//'   time_delay = 1, alpha = 1, beta = 0, type = c(2, 2), type_t12 = c(1, 1))
//' 
//' lathproj <- projection3(ehrlen3, nreps = 5, stochastic = TRUE, substoch = 2,
//'   density = e3d)
//' }
//' 
//' @export density_input
// [[Rcpp::export(density_input)]]
DataFrame density_input(List mpm, RObject stage3, RObject stage2,
  Nullable<RObject> stage1 = R_NilValue, Nullable<RObject> age2 = R_NilValue,
  Nullable<RObject> style = R_NilValue, Nullable<RObject> time_delay = R_NilValue,
  Nullable<RObject> alpha = R_NilValue, Nullable<RObject> beta = R_NilValue,
  Nullable<RObject> type = R_NilValue, Nullable<RObject> type_t12 = R_NilValue) {
  
  bool historical = false;
  bool agebystage = false;
  
  // Check quality of mpm input
  StringVector mpm_class_vec  = mpm.attr("class");
  std::string mpm_class = as<std::string>(mpm_class_vec(0));
  
  CharacterVector mpm_elems = mpm.names();
  int mpm_name_check = 0;
  for (int i = 0; i < mpm_elems.length(); i++) {
    if (stringcompare_hard(as<std::string>(mpm_elems(i)), "ahstages")) mpm_name_check++;
    if (stringcompare_hard(as<std::string>(mpm_elems(i)), "hstages")) mpm_name_check++;
    if (stringcompare_hard(as<std::string>(mpm_elems(i)), "agestages")) mpm_name_check++;
  }
  
  if (!(mpm_name_check == 3 && stringcompare_hard(mpm_class, "lefkoMat"))) {
    throw Rcpp::exception("This function requires a lefkoMat object as input.", false);
  }
  
  DataFrame stageframe = as<DataFrame>(mpm["ahstages"]);
  DataFrame hstages = as<DataFrame>(mpm["hstages"]);
  DataFrame agestages = as<DataFrame>(mpm["agestages"]);
  
  if (hstages.length() > 1) {
    historical = true;
  }
  if (agestages.length() > 1) {
    agebystage = true;
  } 
  if (stageframe.length() == 1) {
    throw Rcpp::exception("Input lefkoMat object does not appear to have a stageframe, which should be element ahstages.",
      false);
  }
  
  CharacterVector ahstages_elems = stageframe.names();
  int ahs_name_check = 0;
  for (int i = 0; i < ahstages_elems.length(); i++) {
    if (stringcompare_hard(as<std::string>(ahstages_elems(i)), "stage_id")) ahs_name_check++;
    if (stringcompare_hard(as<std::string>(ahstages_elems(i)), "stage")) ahs_name_check++;
  }
  
  if (ahs_name_check != 2) {
    throw Rcpp::exception("Stageframe appears to be modified. Please make sure that element ahstages in the lefkoMat object includes both a stage column holding stage names and a stage_id column holding unique, stage identifying integers.", 
      false);
  }
  
  IntegerVector stage_id_sf = as<IntegerVector>(stageframe["stage_id"]);
  StringVector stage_sf = as<StringVector>(stageframe["stage"]);
  int no_stages = stage_sf.length();
  
  // Density input vector standardization
  StringVector stage3_names;
  StringVector stage2_names;
  StringVector stage1_names;
  
  if (is<StringVector>(stage3)) {
    StringVector stage3_names_ = as<StringVector>(stage3);
    stage3_names = stage3_names_;
    
  } else if (is<IntegerVector>(stage3)) {
    arma::ivec stage3_ids = as<arma::ivec>(stage3);
    int stage3_entries = stage3_ids.n_elem;
    
    arma::uvec bad_lows = find(stage3_ids < 1);
    arma::uvec bad_highs = find(stage3_ids > no_stages);
    
    if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
      throw Rcpp::exception("Stageframe contains invalid entries in stage_id column.", false);
    }
    
    StringVector new_stage3 (stage3_entries);
    for (int i = 0; i < stage3_entries; i++) {
      new_stage3(i) = stage_sf(stage3_ids(i) - 1);
    }
    
    stage3_names = new_stage3;
    
  } else {
    throw Rcpp::exception("Option stage3 is not a recognized input type.", false);
  }
  
  if (is<StringVector>(stage2)) {
    StringVector stage2_names_ = as<StringVector>(stage2);
    stage2_names = stage2_names_;
    
  } else if (is<IntegerVector>(stage2)) {
    arma::ivec stage2_ids = as<arma::ivec>(stage2);
    int stage2_entries = stage2_ids.n_elem;
    
    arma::uvec bad_lows = find(stage2_ids < 1);
    arma::uvec bad_highs = find(stage2_ids > no_stages);
    
    if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
      throw Rcpp::exception("Stageframe contains invalid entries in stage_id column.", false);
    }
    
    StringVector new_stage2 (stage2_entries);
    for (int i = 0; i < stage2_entries; i++) {
      new_stage2(i) = stage_sf(stage2_ids(i) - 1);
    }
    
    stage2_names = new_stage2;
    
  } else {
    throw Rcpp::exception("Option stage2 is not a recognized input type.", false);
  }
  
  if (stage1.isNotNull()) {
    if (is<IntegerVector>(stage1)) {
      arma::ivec stage1_ids = as<arma::ivec>(stage1);
      int stage1_entries = stage1_ids.n_elem;
      
      arma::uvec bad_lows = find(stage1_ids < 1);
      arma::uvec bad_highs = find(stage1_ids > no_stages);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Stageframe contains invalid entries in stage_id column.", false);
      }
      
      StringVector new_stage1 (stage1_entries);
      for (int i = 0; i < stage1_entries; i++) {
        new_stage1(i) = stage_sf(stage1_ids(i) - 1);
      }
      
      stage1_names = new_stage1;
      
    } else if (is<StringVector>(stage1)) {
      StringVector stage1_names_ = as<StringVector>(stage1);
      stage1_names = stage1_names_;
      
    } else {
      throw Rcpp::exception("Option stage1 is not a recognized input type.", false);
    }
    
    if (!historical) {
      throw Rcpp::exception("Do not include option stage1 for ahistorical MPMs.", false);
    }
    
  } else if (historical) {
    throw Rcpp::exception("Option stage1 is needed for historical MPMs.", false);
    
  } else {
    StringVector stage1_names_ = {NA_STRING};
    stage1_names = stage1_names_;
  }
  
  if (stage3_names.length() != stage2_names.length()) {
    throw Rcpp::exception("All transitions to modify require information at least for stage2 and stage3. These inputs must also be of equal length.", 
      false);
  }
  
  if (historical) {
    if (stage3_names.length() != stage1_names.length()) {
      throw Rcpp::exception("All historical transitions to modify require information for stage1, stage2, and stage3. These inputs must also be of equal length.", 
        false);
    }
  }
  
  IntegerVector age2_vec;
  IntegerVector style_vec;
  IntegerVector time_delay_vec;
  NumericVector alpha_vec;
  NumericVector beta_vec;
  IntegerVector type_vec;
  IntegerVector type_t12_vec;
  
  if (age2.isNotNull()) {
    if (!agebystage) {
      throw Rcpp::exception("Do not use the age2 option unless the MPM is age-by-stage.", false);
    }
    
    if (is<IntegerVector>(age2)) {
      IntegerVector age2_vec_ = as<IntegerVector>(age2);
      age2_vec = age2_vec_;
      
      arma::ivec age2_arma = as<arma::ivec>(age2_vec_);
      arma::uvec bad_lows = find(age2_arma < 0);
      
      if (bad_lows.n_elem > 0) {
        throw Rcpp::exception("Negative ages are not allowed.", false);
      }
    } else {
      throw Rcpp::exception("Only integers are allowed in option age2.", false);
    }
  } else {
    IntegerVector age2_vec_ = {NA_INTEGER};
    age2_vec = age2_vec_;
  }
  
  if (style.isNotNull()) {
    StringVector ricker_style = {"1", "ricker", "ricke", "rick", "ric", "ri", "r"};
    StringVector beverton_style = {"2", "beverton", "beverto", "bevert", "bever",
      "beve", "bev", "be", "b", "holt", "hol", "ho", "h"};
    StringVector usher_style = {"3", "usher", "ushe", "ush", "us", "u"};
    StringVector logistic_style = {"4", "logistic", "logisti", "logist", "logis",
      "logi", "log", "lo", "l"};
    
    if (is<IntegerVector>(style)) {
      IntegerVector style_ = as<IntegerVector>(style);
      style_vec = style_;
      
      arma::ivec style_arma = as<arma::ivec>(style_);
      
      arma::uvec bad_lows = find(style_arma < 1);
      arma::uvec bad_highs = find(style_arma > 4);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Invalid density dependence style entered.", false);
      }
      
    } else if (is<NumericVector>(style)) {
      IntegerVector style_ = as<IntegerVector>(style);
      style_vec = style_;
      
      arma::ivec style_arma = as<arma::ivec>(style_);
      
      arma::uvec bad_lows = find(style_arma < 1);
      arma::uvec bad_highs = find(style_arma > 4);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Invalid density dependence style entered.", false);
      }
      
    } else if (is<StringVector>(style)) {
      StringVector style_stringvec = as<StringVector>(style);
      int style_elems = style_stringvec.length();
      
      IntegerVector style_vec_ (style_elems, 0);
      
      for (int i = 0; i < style_elems; i++) {
        std::string ssv = as<std::string>(style_stringvec(i));
        for (int j = 0; j < ssv.size(); j++) {
          ssv[j] = tolower(ssv[j]);
        }
        
        for (int j = 0; j < ricker_style.length(); j++) {
          if (stringcompare_hard(ssv, as<std::string>(ricker_style(j)))) {
            style_vec_(i) = 1;
          }
        }
        
        for (int j = 0; j < beverton_style.length(); j++) {
          if (stringcompare_hard(ssv, as<std::string>(beverton_style(j)))) {
            style_vec_(i) = 2;
          }
        }
        
        for (int j = 0; j < usher_style.length(); j++) {
          if (stringcompare_hard(ssv, as<std::string>(usher_style(j)))) {
            style_vec_(i) = 3;
          }
        }
        
        for (int j = 0; j < logistic_style.length(); j++) {
          if (stringcompare_hard(ssv, as<std::string>(logistic_style(j)))) {
            style_vec_(i) = 1;
          }
        }
      }
      
      if (min(style_vec_) == 0) {
        throw Rcpp::exception("Some density dependence styles were not recognized.", false);
      }
      
      style_vec = style_vec_;
      
    } else {
      throw Rcpp::exception("Invalid entry for option style.", false);
      
    }
  } else {
    IntegerVector style_vec_ = {1};
    style_vec = style_vec_;
  }
  
  if (time_delay.isNotNull()) {
    if (is<IntegerVector>(time_delay)) {
      IntegerVector time_delay_vec_ = as<IntegerVector>(time_delay);
      time_delay_vec = time_delay_vec_;
      
      arma::ivec time_delay_arma = as<arma::ivec>(time_delay_vec_);
      arma::uvec bad_lows = find(time_delay_arma < 1);
      
      if (bad_lows.n_elem > 0) {
        throw Rcpp::exception("Time delays less than 1 time step are not allowed.", false);
      }
    } else if (is<NumericVector>(time_delay)) {
      IntegerVector time_delay_vec_ = as<IntegerVector>(time_delay);
      time_delay_vec = time_delay_vec_;
      
      arma::ivec time_delay_arma = as<arma::ivec>(time_delay_vec_);
      arma::uvec bad_lows = find(time_delay_arma < 1);
      
      if (bad_lows.n_elem > 0) {
        throw Rcpp::exception("Time delays less than 1 time step are not allowed.", false);
      }
    } else {
      throw Rcpp::exception("Only integers are allowed in option time_delay.", false);
    }
  } else {
    IntegerVector time_delay_vec_ = {1};
    time_delay_vec = time_delay_vec_;
  }
  
  if (alpha.isNotNull()) {
    if (is<NumericVector>(alpha)) {
      NumericVector alpha_vec_ = as<NumericVector>(alpha);
      alpha_vec = alpha_vec_;
      
    } else {
      throw Rcpp::exception("Only floating point decimals are allowed in option alpha.", false);
    }
    
  } else {
    NumericVector alpha_vec_ = {1};
    alpha_vec = alpha_vec_;
  }
  
  if (beta.isNotNull()) {
    if (is<NumericVector>(beta)) {
      NumericVector beta_vec_ = as<NumericVector>(beta);
      beta_vec = beta_vec_;
      
    } else {
      throw Rcpp::exception("Only floating point decimals are allowed in option beta.", false);
    }
    
  } else {
    NumericVector beta_vec_ = {1};
    beta_vec = beta_vec_;
  }
  
  if (type.isNotNull()) {
    if (is<IntegerVector>(type)) {
      IntegerVector type_vec_ = as<IntegerVector>(type);
      type_vec = type_vec_;
      
      arma::ivec type_arma = as<arma::ivec>(type_vec_);
      arma::uvec bad_lows = find(type_arma < 1);
      arma::uvec bad_highs = find(type_arma > 2);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Transition types may only be type 1 or 2.", false);
      }
    } else if (is<StringVector>(type)) {
      StringVector type_stringvec = as<StringVector>(type);
      int type_elems = type_stringvec.length();
      
      IntegerVector type_vec_ (type_elems, 0);
      
      for (int i = 0; i < type_elems; i++) {
        if (stringcompare_simple(as<std::string>(type_stringvec(i)), "s", true)) {
          type_vec_(i) = 1;
        } else if (stringcompare_simple(as<std::string>(type_stringvec(i)), "1", true)) {
          type_vec_(i) = 1;
        } else if (stringcompare_simple(as<std::string>(type_stringvec(i)), "f", true)) {
          type_vec_(i) = 2;
        } else if (stringcompare_simple(as<std::string>(type_stringvec(i)), "2", true)) {
          type_vec_(i) = 2;
        } else {
          throw Rcpp::exception("Invalid entry in option type.", false);
        }
      }
      type_vec = type_vec_;
      
    } else if (is<NumericVector>(type)) {
      IntegerVector type_vec_ = as<IntegerVector>(type);
      type_vec = type_vec_;
      
      arma::ivec type_arma = as<arma::ivec>(type_vec_);
      arma::uvec bad_lows = find(type_arma < 1);
      arma::uvec bad_highs = find(type_arma > 2);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Transition types may only be type 1 or 2.", false);
      }
      
    } else {
      throw Rcpp::exception("Only integers are allowed in option type.", false);
    }
    
  } else {
    IntegerVector type_vec_ = {1};
    type_vec = type_vec_;
  }
  
  if (type_t12.isNotNull()) {
    if (is<IntegerVector>(type_t12)) {
      IntegerVector type_t12_vec_ = as<IntegerVector>(type_t12);
      type_t12_vec = type_t12_vec_;
      
      arma::ivec type_t12_arma = as<arma::ivec>(type_t12_vec_);
      arma::uvec bad_lows = find(type_t12_arma < 1);
      arma::uvec bad_highs = find(type_t12_arma > 2);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Historical transition types may only be type 1 or 2.", false);
      }
    } else if (is<NumericVector>(type_t12)) {
      IntegerVector type_t12_vec_ = as<IntegerVector>(type_t12);
      type_t12_vec = type_t12_vec_;
      
      arma::ivec type_t12_arma = as<arma::ivec>(type_t12_vec_);
      arma::uvec bad_lows = find(type_t12_arma < 1);
      arma::uvec bad_highs = find(type_t12_arma > 2);
      
      if (bad_lows.n_elem > 0 || bad_highs.n_elem > 0) {
        throw Rcpp::exception("Historical transition types may only be type 1 or 2.", false);
      }
      
    } else if (is<StringVector>(type_t12)) {
      StringVector type_t12_stringvec = as<StringVector>(type_t12);
      int type_t12_elems = type_t12_stringvec.length();
      
      IntegerVector type_t12_vec_ (type_t12_elems, 0);
      
      for (int i = 0; i < type_t12_elems; i++) {
        if (stringcompare_simple(as<std::string>(type_t12_stringvec(i)), "s", true)) {
          type_t12_vec_(i) = 1;
        } else if (stringcompare_simple(as<std::string>(type_t12_stringvec(i)), "1", true)) {
          type_t12_vec_(i) = 1;
        } else if (stringcompare_simple(as<std::string>(type_t12_stringvec(i)), "f", true)) {
          type_t12_vec_(i) = 2;
        } else if (stringcompare_simple(as<std::string>(type_t12_stringvec(i)), "2", true)) {
          type_t12_vec_(i) = 2;
        } else {
          throw Rcpp::exception("Invalid entry in option type_t12.", false);
        }
      }
      type_t12_vec = type_t12_vec_;
      
    } else {
      throw Rcpp::exception("Only integers are allowed in option type_t12.", false);
    }
    
  } else {
    IntegerVector type_t12_vec_ = {1};
    type_t12_vec = type_t12_vec_;
  }
  
  
  StringVector new_stage3_names;
  StringVector new_stage2_names;
  StringVector new_stage1_names;
  IntegerVector new_age2_vec;
  IntegerVector new_style_vec;
  IntegerVector new_time_delay_vec;
  NumericVector new_alpha_vec;
  NumericVector new_beta_vec;
  IntegerVector new_type_vec;
  IntegerVector new_type_t12_vec;
  
  IntegerVector vec_lengths = {static_cast<int>(stage3_names.length()),
    static_cast<int>(stage2_names.length()), static_cast<int>(stage1_names.length()),
    static_cast<int>(age2_vec.length()), static_cast<int>(style_vec.length()),
    static_cast<int>(time_delay_vec.length()), static_cast<int>(alpha_vec.length()),
    static_cast<int>(beta_vec.length()), static_cast<int>(type_vec.length()),
    static_cast<int>(type_t12_vec.length())};
  int vec_max = max(vec_lengths);
  
  if (stage3_names.length() != vec_max) {
    if (stage3_names.length() == 1) {
      StringVector ns3n (vec_max, stage3_names(0));
      new_stage3_names = ns3n;
      
    } else {
      throw Rcpp::exception("Vector stage3 is not the correct length.", false);
    }
    
  } else {
    new_stage3_names = stage3_names;
  }
  
  if (stage2_names.length() != vec_max) {
    if (stage2_names.length() == 1) {
      StringVector ns2n (vec_max, stage2_names(0));
      new_stage2_names = ns2n;
      
    } else {
      throw Rcpp::exception("Vector stage2 is not the correct length.", false);
    }
    
  } else {
    new_stage2_names = stage2_names;
  }
  
  if (stage1_names.length() != vec_max) {
    if (stage1_names.length() == 1) {
      StringVector ns1n (vec_max, stage1_names(0));
      new_stage1_names = ns1n;
      
    } else {
      throw Rcpp::exception("Vector stage1 is not the correct length.", false);
    }
    
  } else {
    new_stage1_names = stage1_names;
  }
  
  if (age2_vec.length() != vec_max) {
    if (age2_vec.length() == 1) {
      IntegerVector na2v (vec_max, age2_vec(0));
      new_age2_vec = na2v;
      
    } else {
      throw Rcpp::exception("Vector age2 is not the correct length.", false);
    }
    
  } else {
    new_age2_vec = age2_vec;
  }
  
  if (style_vec.length() != vec_max) {
    if (style_vec.length() == 1) {
      IntegerVector ns2v (vec_max, style_vec(0));
      new_style_vec = ns2v;
      
    } else {
      throw Rcpp::exception("Vector style is not the correct length.", false);
    }
    
  } else {
    new_style_vec = style_vec;
  }
  
  if (time_delay_vec.length() != vec_max) {
    if (time_delay_vec.length() == 1) {
      IntegerVector ntd2v (vec_max, time_delay_vec(0));
      new_time_delay_vec = ntd2v;
      
    } else {
      throw Rcpp::exception("Vector time_delay is not the correct length.", false);
    }
    
  } else {
    new_time_delay_vec = time_delay_vec;
  }
  
  if (type_vec.length() != vec_max) {
    if (type_vec.length() == 1) {
      IntegerVector nt2v (vec_max, type_vec(0));
      new_type_vec = nt2v;
      
    } else {
      throw Rcpp::exception("Vector type is not the correct length.", false);
    }
    
  } else {
    new_type_vec = type_vec;
  }
  
  if (type_t12_vec.length() != vec_max) {
    if (type_t12_vec.length() == 1) {
      IntegerVector ntt122v (vec_max, type_t12_vec(0));
      new_type_t12_vec = ntt122v;
      
    } else {
      throw Rcpp::exception("Vector type_t12 is not the correct length.", false);
    }
    
  } else {
    new_type_t12_vec = type_t12_vec;
  }
  
  if (alpha_vec.length() != vec_max) {
    if (alpha_vec.length() == 1) {
      NumericVector nav (vec_max, alpha_vec(0));
      new_alpha_vec = nav;
      
    } else {
      throw Rcpp::exception("Vector alpha is not the correct length.", false);
    }
    
  } else {
    new_alpha_vec = alpha_vec;
  }
  
  if (beta_vec.length() != vec_max) {
    if (beta_vec.length() == 1) {
      NumericVector nbv (vec_max, beta_vec(0));
      new_beta_vec = nbv;
      
    } else {
      throw Rcpp::exception("Vector beta is not the correct length.", false);
    }
    
  } else {
    new_beta_vec = beta_vec;
  }
  
  DataFrame output_tab = DataFrame::create(_["stage3"] = new_stage3_names,
    _["stage2"] = new_stage2_names, _["stage1"] = new_stage1_names,
    _["age2"] = new_age2_vec, _["style"] = new_style_vec,
    _["time_delay"] = new_time_delay_vec, _["alpha"] = new_alpha_vec,
    _["beta"] = new_beta_vec, _["type"] = new_type_vec,
    _["type_t12"] = new_type_t12_vec);
  
  DataFrame output = density_reassess(stageframe, output_tab, agestages,
    historical, agebystage);
  
  StringVector needed_classes {"data.frame", "lefkoDens"};
  output.attr("class") = needed_classes;

  StringVector stage3_final = as<StringVector>(output["stage3"]);
  StringVector stage2_final = as<StringVector>(output["stage2"]);
  StringVector stage1_final = as<StringVector>(output["stage1"]);
  StringVector age2_final = as<StringVector>(output["age2"]);
  
  int final_rows = stage3_final.length();
  
  StringVector check_elems(final_rows);
  
  for (int i = 0; i < final_rows; i++) {
    check_elems(i) = stage3_final(i);
    check_elems(i) += " ";
    check_elems(i) += stage2_final(i);
    check_elems(i) += " ";
    if (!StringVector::is_na(stage1_final(i))) {
      check_elems(i) += stage1_final(i);
      check_elems(i) += " ";
    }
    if (!StringVector::is_na(age2_final(i))) check_elems(i) += age2_final(i);
  }
  
  StringVector unique_elems = unique(check_elems);
  if (unique_elems.length() != final_rows) {
    Rf_warningcall(R_NilValue, "Some transitions appear to be listed multiple times. This may cause errors in analysis.");
  }
  
  return output;
}

