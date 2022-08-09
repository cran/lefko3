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
