#include <RcppArmadillo.h>
#include <LefkoUtils.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace LefkoUtils;
using namespace LefkoMats;

//' Main Formula Creation for Function \code{modelsearch()}
//'
//' Function \code{praxis()} is the workhorse function used by function
//' \code{stovokor} to create individual vital rate model formulae, which are
//' then used as input in the global model calls used in function
//' \code{\link{modelsearch}()}.
//' 
//' @name praxis
//' 
//' @param surv A vector of strings indicating the names of the variables coding
//' survival.
//' @param obs A vector of strings indicating the names of the variables coding
//' observation status.
//' @param size A vector of strings indicating the names of the variables coding
//' primary size.
//' @param sizeb A vector of strings indicating the names of the variables
//' coding secondary size.
//' @param sizec A vector of strings indicating the names of the variables
//' coding tertiary size.
//' @param repst A vector of strings indicating the names of the variables
//' coding reproductive status.
//' @param fec A vector of strings indicating the names of the variables coding
//' fecundity.
//' @param matstat A vector of strings indicating the names of the variables
//' coding for maturity status.
//' @param historical A logical value indicating whether to create global models
//' with historical effects.
//' @param response An integer coding for the vital rate to be modeled. Codes
//' are as follows: 1) survival, 2) observation status, 3) primary size,
//' 4) secondary size, 5) tertiary size, 6) reproductive status, 7) fecundity,
//' and 8) maturity status.
//' @param suite A string indicating the scope of independent factors included
//' in the global models. Options include \code{"full"}, \code{"main"},
//' \code{"size"}, \code{"rep"}, and \code{"const"}.
//' @param approach A string indicating whether to use mixed model encoding 
//' (\code{"mixed"}) or GLM encoding (\code{"glm"}).
//' @param nojuvs A logical value indicating that juvenile rates should be
//' estimated (\code{FALSE}) or not (\code{TRUE}).
//' @param juvsize A logical value indicating whether to include size terms in
//' juvenile models.
//' @param indiv A string indicating the name of the variable coding individual
//' identity.
//' @param patch A string indicating the name of the variable coding patch
//' identity.
//' @param year A string indicating the name of the variable coding time
//' \emph{t}.
//' @param age A string indicating the name of the variable coding age.
//' @param densitycol The name of the density variable, or \code{"none"}.
//' @param indcova_raw A vector of strings indicating the names in times
//' \emph{t}+1, \emph{t}, and \emph{t}-1 of a specific individual covariate used
//' in the dataset.
//' @param indcovb_raw A vector of strings indicating the names in times
//' \emph{t}+1, \emph{t}, and \emph{t}-1 of a specific individual covariate used
//' in the dataset.
//' @param indcovc_raw A vector of strings indicating the names in times
//' \emph{t}+1, \emph{t}, and \emph{t}-1 of a specific individual covariate used
//' in the dataset.
//' @param sizebused A logical value denoting if secondary size variables are to
//' be used.
//' @param sizecused A logical value denoting if tertiary size variables are to
//' be used.
//' @param grouptest A logical value indicating whether to test for group
//' effect.
//' @param ageused A logical value indicating whether to test for age effect.
//' @param densityused A logical value indicating whether the density variable
//' is to be used.
//' @param indcovaused Logical value indicating whether individual covariate a
//' is used.
//' @param indcovbused Logical value indicating whether individual covariate b
//' is used.
//' @param indcovcused Logical value indicating whether individual covariate c
//' is used.
//' @param pasrand A logical value indicating whether to treat patch as a random
//' variable within mixed models.
//' @param yasrand A logical value indicating whether to treat year as a random
//' variable within mixed models.
//' @param iaasrand A logical value indicating whether to treat indcova as
//' random.
//' @param ibasrand A logical value indicating whether to treat indcovb as
//' random.
//' @param icasrand A logical value indicating whether to treat indcovc as
//' random.
//' @param iaasfac A logical value indicating whether to treat indcova as a
//' factor variable.
//' @param ibasfac A logical value indicating whether to treat indcovb as a
//' factor variable.
//' @param icasfac A logical value indicating whether to treat indcovc as a
//' factor variable.
//' @param fectime An integer indicating whether to use reproductive output in
//' time \emph{t} (2) or time \emph{t}+1 (3) as the response for fecundity.
//' @param repstcheck A boolean variable denoting whether reproductive status is
//' being tested as a response within the suite of vital rates being estimated.
//' 
//' @return A list with four elements. The first and second are both one-element
//' string vectors, with the first coding for the adult vital rate global model,
//' and the second coding the juvenile vital rate global model. The third and
//' fourth code the number of terms tested in each of these models,
//' respectively.
//'
//' @keywords internal
//' @noRd
List praxis(const StringVector& surv, const StringVector& obs,
  const StringVector& size, const StringVector& sizeb,
  const StringVector& sizec, const StringVector& repst, const StringVector& fec,
  const StringVector& matstat, bool historical, int response,
  const String& suite, const String& approach, const bool nojuvs,
  bool juvsize, const String& indiv, const String& patch,
  const String& year, const String& age, const String& densitycol,
  const StringVector indcova_raw, const StringVector indcovb_raw,
  const StringVector indcovc_raw, const bool sizebused, const bool sizecused,
  const bool grouptest, const bool ageused, const bool densityused,
  const bool indcovaused, const bool indcovbused, const bool indcovcused,
  bool pasrand, bool yasrand, bool iaasrand, bool ibasrand, bool icasrand,
  const bool iaasfac, const bool ibasfac, const bool icasfac, int fectime,
  bool repstcheck) {
  
  if (nojuvs) juvsize = false;
  
  String randomtackonp {""};
  String randomtackony {""};
  String randomtackoni {""};
  String randomtackonia {""};
  String randomtackonib {""};
  String randomtackonic {""};
  String randomtackon {""};
  
  String fixedtackong {""};
  String fixedtackonp {""};
  String fixedtackony {""};
  String fixedtackonia {""};
  String fixedtackonib {""};
  String fixedtackonic {""};
  String fixedtackon {""};
  
  String sizesuffix;
  String jsizesuffix;
  String fecsuffix;
  
  String fullmainmodel;
  String juvmainmodel;
  String adult_out_model;
  String juv_out_model;
  
  int covcount {0};
  if (indcova_raw(1) != "none" && indcovaused) covcount += 1;
  if (indcovb_raw(1) != "none" && indcovbused) covcount += 1;
  if (indcovc_raw(1) != "none" && indcovcused) covcount += 1;
  
  int total_terms {0};
  int fixedcovcounter {0};
  int randomcovcounter {0};
  
  StringVector indcova (3);
  StringVector indcovb (3);
  StringVector indcovc (3);
  
  if (indcovaused) {
    if (!iaasrand && iaasfac) {
      String major_pain1 = "as.factor(";
      major_pain1 += indcova_raw(1);
      major_pain1 += ")";
      
      String major_pain2 = "as.factor(";
      major_pain2 += indcova_raw(2);
      major_pain2 += ")";
      
      indcova(1) = major_pain1;
      indcova(2) = major_pain2;
    } else indcova = indcova_raw;
  }
  if (indcovbused) {
    if (!ibasrand && ibasfac) {
      String major_pain1 = "as.factor(";
      major_pain1 += indcovb_raw(1);
      major_pain1 += ")";
      
      String major_pain2 = "as.factor(";
      major_pain2 += indcovb_raw(2);
      major_pain2 += ")";
      
      indcovb(1) = major_pain1;
      indcovb(2) = major_pain2;
    } else indcovb = indcovb_raw;
  }
  if (indcovcused) {
    if (!icasrand && icasfac) {
      String major_pain1 = "as.factor(";
      major_pain1 += indcovc_raw(1);
      major_pain1 += ")";
      
      String major_pain2 = "as.factor(";
      major_pain2 += indcovc_raw(2);
      major_pain2 += ")";
      
      indcovc(1) = major_pain1;
      indcovc(2) = major_pain2;
    } else indcovc = indcovc_raw;
  }
  
  // Year, patch, indiv, and group - all possible factor variables
  // Year, patch, and indiv may be random
  if (approach != "mixed") {
    yasrand = false;
    pasrand = false;
  }
  
  if (grouptest) {
    fixedtackong += "as.factor(group2)";
    fixedcovcounter += 1;
    total_terms++;
    
    if (historical) {
      fixedcovcounter += 1;
      fixedtackong += " + as.factor(group1)";
      total_terms++;
    }
    fixedtackon += fixedtackong;
  }
  
  if (year!= "none") {
    if (yasrand) {
      randomtackony += "(1 | ";
      randomtackony += year;
      randomtackony += ")";
      randomcovcounter += 1;
      total_terms++;
      
      randomtackon += randomtackony;
      
    } else {
      fixedtackony += "as.factor(";
      fixedtackony += year;
      fixedtackony += ")";
      total_terms++;
      fixedcovcounter += 1;
      
      if (fixedcovcounter > 1) fixedtackon += " + ";
      fixedtackon += fixedtackony;
    }
  }
  
  if (patch!= "none") {
    if (pasrand) {
      randomtackonp += "(1 | ";
      randomtackonp += patch;
      randomtackonp += ")";
      randomcovcounter += 1;
      total_terms++;
      
      if (randomcovcounter > 1) randomtackon += " + ";
      randomtackon += randomtackonp;
      
    } else {
      fixedtackonp += "as.factor(";
      fixedtackonp += patch;
      fixedtackonp += ")";
      total_terms++;
      fixedcovcounter += 1;
      
      if (fixedcovcounter > 1) fixedtackon += " + ";
      fixedtackon += fixedtackonp;
    }
  }
  
  if (indiv != "none" && approach == "mixed") {
    randomtackoni += "(1 | ";
    randomtackoni += indiv;
    randomtackoni += ")";
    total_terms++;
    randomcovcounter += 1;
    
    if (randomcovcounter > 1) randomtackon += " + ";
    randomtackon += randomtackoni;
  }
  
  // Add individual covariates to tacked-on sections
  if (indcova(1) != "none" && indcovaused) {
    if (!iaasrand) {
      fixedtackonia += indcova(1);
      fixedcovcounter += 1;
      total_terms++;
      
      if (historical) {
        fixedtackonia += " + ";
        fixedtackonia += indcova(2);
        fixedcovcounter += 1;
        total_terms++;
      }
      
      if ((fixedcovcounter > 2 && historical) || (fixedcovcounter > 1 && !historical))
        fixedtackon += " + ";
      fixedtackon += fixedtackonia;
      
    } else {
      randomtackonia += "(1 | ";
      randomtackonia += indcova(1);
      randomtackonia += ")";
      randomcovcounter += 1;
      total_terms++;
      
      if (historical) {
        randomtackonia += " + (1 | ";
        randomtackonia += indcova(2);
        randomtackonia += ")";
        randomcovcounter += 1;
        total_terms++;
      }
      
      if ((randomcovcounter > 2 && historical) || (randomcovcounter > 1 && !historical))
        randomtackon += " + ";
      randomtackon += randomtackonia;
    }
  }
  
  if (indcovb(1) != "none" && indcovbused) {
    if (!ibasrand) {
      fixedtackonib += indcovb(1);
      fixedcovcounter += 1;
      total_terms++;
      
      if (historical) {
        fixedtackonib += " + ";
        fixedtackonib += indcovb(2);
        fixedcovcounter += 1;
        total_terms++;
      }
      
      if ((fixedcovcounter > 2 && historical) || (fixedcovcounter > 1 && !historical))
        fixedtackon += " + ";
      fixedtackon += fixedtackonib;
      
    } else {
      randomtackonib += "(1 | ";
      randomtackonib += indcovb(1);
      randomtackonib += ")";
      randomcovcounter += 1;
      total_terms++;
      
      if (historical) {
        randomtackonib += " + (1 | ";
        randomtackonib += indcovb(2);
        randomtackonib += ")";
        randomcovcounter += 1;
        total_terms++;
      }
      
      if ((randomcovcounter > 2 && historical) || (randomcovcounter > 1 && !historical))
        randomtackon += " + ";
      randomtackon += randomtackonib;
    }
  }
  
  if (indcovc(1) != "none" && indcovcused) {
    if (!icasrand) {
      fixedtackonic += indcovc(1);
      fixedcovcounter += 1;
      total_terms++;
      
      if (historical) {
        fixedtackonic += " + ";
        fixedtackonic += indcovc(2);
        total_terms++;
      }
      
      if ((fixedcovcounter > 2 && historical) || (fixedcovcounter > 1 && !historical))
        fixedtackon += " + ";
      fixedtackon += fixedtackonic;
      
    } else {
      randomtackonic += "(1 | ";
      randomtackonic += indcovc(1);
      randomtackonic += ")";
      randomcovcounter += 1;
      total_terms++;
      
      if (historical) {
        randomtackonic += " + (1 | ";
        randomtackonic += indcovc(2);
        randomtackonic += ")";
        randomcovcounter += 1;
        total_terms++;
      }
      
      if ((randomcovcounter > 2 && historical) || (randomcovcounter > 1 && !historical))
        randomtackon += " + ";
      randomtackon += randomtackonic;
    }
  }
  
  if (suite == "full" && !iaasrand && !ibasrand) { 
    if ((indcova(1) != "none" && indcovb(1) != "none") && (indcovaused && indcovbused)) {
      fixedtackonib += " + ";
      fixedtackonib += indcova(1);
      fixedtackonib += ":";
      fixedtackonib += indcovb(1);
      fixedcovcounter += 1;
      total_terms++;
      
      if (historical) {
        fixedtackonib += " + ";
        fixedtackonib += indcova(2);
        fixedtackonib += ":";
        fixedtackonib += indcovb(2);
        fixedcovcounter += 1;
        total_terms++;
        
        fixedtackonib += " + ";
        fixedtackonib += indcova(1);
        fixedtackonib += ":";
        fixedtackonib += indcovb(2);
        fixedcovcounter += 1;
        total_terms++;
        
        fixedtackonib += " + ";
        fixedtackonib += indcova(2);
        fixedtackonib += ":";
        fixedtackonib += indcovb(1);
        fixedcovcounter += 1;
        total_terms++;
      }
    }
  }
  if (suite == "full" && !iaasrand && !icasrand) {
    if ((indcova(1) != "none" && indcovc(1) != "none") && (indcovaused && indcovcused)) {
      fixedtackonic += " + ";
      fixedtackonic += indcova(1);
      fixedtackonic += ":";
      fixedtackonic += indcovc(1);
      fixedcovcounter += 1;
      total_terms++;
      
      if (historical) {
        fixedtackonic += " + ";
        fixedtackonic += indcova(2);
        fixedtackonic += ":";
        fixedtackonic += indcovc(2);
        fixedcovcounter += 1;
        total_terms++;
        
        fixedtackonic += " + ";
        fixedtackonic += indcova(1);
        fixedtackonic += ":";
        fixedtackonic += indcovc(2);
        fixedcovcounter += 1;
        total_terms++;
        
        fixedtackonic += " + ";
        fixedtackonic += indcova(2);
        fixedtackonic += ":";
        fixedtackonic += indcovc(1);
        fixedcovcounter += 1;
        total_terms++;
      }
    }
  }
  if (suite == "full" && !ibasrand && !icasrand) {
    if ((indcovb(1) != "none" && indcovc(1) != "none") && (indcovbused && indcovcused)) {
      fixedtackonic += " + ";
      fixedtackonic += indcovb(1);
      fixedtackonic += ":";
      fixedtackonic += indcovc(1);
      fixedcovcounter += 1;
      total_terms++;
      
      if (historical) {
        fixedtackonic += " + ";
        fixedtackonic += indcovb(2);
        fixedtackonic += ":";
        fixedtackonic += indcovc(2);
        fixedcovcounter += 1;
        total_terms++;
        
        fixedtackonic += " + ";
        fixedtackonic += indcovb(1);
        fixedtackonic += ":";
        fixedtackonic += indcovc(2);
        fixedcovcounter += 1;
        total_terms++;
        
        fixedtackonic += " + ";
        fixedtackonic += indcovb(2);
        fixedtackonic += ":";
        fixedtackonic += indcovc(1);
        fixedcovcounter += 1;
        total_terms++;
      }
    }
  }
  
  // Main model patterns
  int modelcounter {0};
  int juvmodelcounter {0};
  
  // First the juvenile model pattern
  if (!nojuvs) {
    juvmainmodel = " ~ ";
    
    if (suite != "const") {
      if (juvsize && suite != "repst") {
        juvmainmodel += size(1);
        juvmodelcounter = 1;
        
        if (sizebused) {
          if (juvmodelcounter > 0) juvmainmodel += " + ";
          juvmainmodel += sizeb(1);
          juvmodelcounter += 1;
          
          if (suite == "full") {
            juvmainmodel += " + ";
            juvmainmodel += size(1);
            juvmainmodel += ":";
            juvmainmodel += sizeb(1);
            juvmodelcounter += 1;
          }
        }
        if (sizecused) {
          if (juvmodelcounter > 0) juvmainmodel += " + ";
          juvmainmodel += sizec(1);
          juvmodelcounter += 1;
          
          if (suite == "full") {
            juvmainmodel += " + ";
            juvmainmodel += size(1);
            juvmainmodel += ":";
            juvmainmodel += sizec(1);
            juvmodelcounter += 1;
            
            if (sizebused) {
              juvmainmodel += " + ";
              juvmainmodel += sizeb(1);
              juvmainmodel += ":";
              juvmainmodel += sizec(1);
              juvmodelcounter += 1;
            }
          }
        }
        if (densityused) {
          if (juvmodelcounter > 0) juvmainmodel += " + ";
          juvmainmodel += densitycol;
          juvmodelcounter += 1;
          
          if (suite == "full") {
            juvmainmodel += " + ";
            juvmainmodel += size(1);
            juvmainmodel += ":";
            juvmainmodel += densitycol;
            juvmodelcounter += 1;
            
            if (sizebused) {
              juvmainmodel += " + ";
              juvmainmodel += sizeb(1);
              juvmainmodel += ":";
              juvmainmodel += densitycol;
              juvmodelcounter += 1;
            }
            if (sizecused) {
              juvmainmodel += " + ";
              juvmainmodel += sizec(1);
              juvmainmodel += ":";
              juvmainmodel += densitycol;
              juvmodelcounter += 1;
            }
          }
        }
      } else if (densityused) {
        if (juvmodelcounter > 0) juvmainmodel += " + ";
        juvmainmodel += densitycol;
        juvmodelcounter += 1;
      } else if (fixedcovcounter == 0) {
        juvmainmodel += "1";
        juvmodelcounter += 1;
      }
    } else  if (fixedcovcounter == 0) {
      juvmainmodel += "1";
      juvmodelcounter += 1;
    }
    
    if (juvmodelcounter == 0) {
      juvmainmodel += "1";
      juvmodelcounter += 1;
    }
    if (fixedcovcounter > 0) {
      if (juvmodelcounter > 0) juvmainmodel += " + ";
      juvmainmodel += fixedtackon;
    }
    if (randomcovcounter > 0) {
      if (juvmodelcounter > 0) juvmainmodel += " + ";
      juvmainmodel += randomtackon;
    }
    
    juvmodelcounter += total_terms;
  }
  
  // Adult model pattern
  fullmainmodel = " ~ ";
  
  if (age != "none") {
    fullmainmodel += age;
    modelcounter = 1;
  }
  
  if (densityused) {
    if (modelcounter > 0) fullmainmodel += " + ";
    fullmainmodel += densitycol;
    modelcounter += 1;
  }
  
  if (suite == "full" || suite == "main" || suite == "size" || suite == "repst") {
    if (suite != "repst") {
      if (modelcounter > 0) fullmainmodel += " + ";
      fullmainmodel += size(1);
      modelcounter += 1;
      
      if (suite == "full" && densityused) {
        fullmainmodel += " + ";
        fullmainmodel += size(1);
        fullmainmodel += ":";
        fullmainmodel += densitycol;
        modelcounter += 1;
      }
      
      if (historical) {
        fullmainmodel += " + ";
        fullmainmodel += size(2);
        modelcounter += 1;

        if (suite != "main") {
          fullmainmodel += " + ";
          fullmainmodel += size(1);
          fullmainmodel += ":";
          fullmainmodel += size(2);
          modelcounter += 1;
        }
        
        if (suite == "full" && densityused) {
          fullmainmodel += " + ";
          fullmainmodel += size(2);
          fullmainmodel += ":";
          fullmainmodel += densitycol;
          modelcounter += 1;
        }
      }
      
      if (sizebused) {
        if (modelcounter > 0) fullmainmodel += " + ";
        fullmainmodel += sizeb(1);
        modelcounter += 1;
        
        if (suite == "full" && densityused) {
          fullmainmodel += " + ";
          fullmainmodel += sizeb(1);
          fullmainmodel += ":";
          fullmainmodel += densitycol;
          modelcounter += 1;
        }
        
        if (historical) {
          fullmainmodel += " + ";
          fullmainmodel += sizeb(2);
          modelcounter += 1;

          if (suite != "main") {
            fullmainmodel += " + ";
            fullmainmodel += sizeb(1);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += size(1);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += sizeb(1);
            fullmainmodel += ":";
            fullmainmodel += size(2);
            modelcounter += 1;
          }
          
          if (suite == "full" && densityused) {
            fullmainmodel += " + ";
            fullmainmodel += sizeb(2);
            fullmainmodel += ":";
            fullmainmodel += densitycol;
            modelcounter += 1;
          }
        }
      }
      
      if (sizecused) {
        if (modelcounter > 0) fullmainmodel += " + ";
        fullmainmodel += sizec(1);
        modelcounter += 1;
        
        if (suite == "full" && densityused) {
          fullmainmodel += " + ";
          fullmainmodel += sizec(1);
          fullmainmodel += ":";
          fullmainmodel += densitycol;
          modelcounter += 1;
        }
        
        if (historical) {
          fullmainmodel += " + ";
          fullmainmodel += sizec(2);
          modelcounter += 1;
          
          if (suite != "main") {
            fullmainmodel += " + ";
            fullmainmodel += sizec(1);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += size(1);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += sizec(1);
            fullmainmodel += ":";
            fullmainmodel += size(2);
            modelcounter += 1;
            
            if (sizebused) {
              fullmainmodel += " + ";
              fullmainmodel += sizeb(1);
              fullmainmodel += ":";
              fullmainmodel += sizec(2);
              modelcounter += 1;
            
              fullmainmodel += " + ";
              fullmainmodel += sizec(1);
              fullmainmodel += ":";
              fullmainmodel += sizeb(2);
              modelcounter += 1;
            }
          }
          
          if (suite == "full" && densityused) {
            fullmainmodel += " + ";
            fullmainmodel += sizec(2);
            fullmainmodel += ":";
            fullmainmodel += densitycol;
            modelcounter += 1;
          }
        }
      }
    }
    
    if (suite != "size") {
      if (modelcounter > 0) fullmainmodel += " + ";
      fullmainmodel += repst(1);
      modelcounter += 1;
      
      if (suite == "full" && densityused) {
        fullmainmodel += " + ";
        fullmainmodel += repst(1);
        fullmainmodel += ":";
        fullmainmodel += densitycol;
        modelcounter += 1;
      }
      
      if (historical) {
        fullmainmodel += " + ";
        fullmainmodel += repst(2);
        modelcounter += 1;
        
        if (suite == "repst" || suite == "full") {
          fullmainmodel += " + ";
          fullmainmodel += repst(1);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          modelcounter += 1;
        }
        
        if (suite == "full" && densityused) {
          fullmainmodel += " + ";
          fullmainmodel += repst(2);
          fullmainmodel += ":";
          fullmainmodel += densitycol;
          modelcounter += 1;
        }
      }
    }
    
    if (suite == "full") {
      fullmainmodel += " + ";
      fullmainmodel += size(1);
      fullmainmodel += ":";
      fullmainmodel += repst(1);
      modelcounter += 1;
      
      if (sizebused) {
        fullmainmodel += " + ";
        fullmainmodel += sizeb(1);
        fullmainmodel += ":";
        fullmainmodel += repst(1);
        modelcounter += 1;
      }
      
      if (sizecused) {
        fullmainmodel += " + ";
        fullmainmodel += sizec(1);
        fullmainmodel += ":";
        fullmainmodel += repst(1);
        modelcounter += 1;
      }
      
      if (historical) {
        fullmainmodel += " + ";
        fullmainmodel += repst(2);
        fullmainmodel += ":";
        fullmainmodel += size(2);
        modelcounter += 1;
        
        fullmainmodel += " + ";
        fullmainmodel += size(1);
        fullmainmodel += ":";
        fullmainmodel += size(2);
        modelcounter += 1;
        
        fullmainmodel += " + ";
        fullmainmodel += repst(1);
        fullmainmodel += ":";
        fullmainmodel += repst(2);
        modelcounter += 1;
        
        fullmainmodel += " + ";
        fullmainmodel += size(1);
        fullmainmodel += ":";
        fullmainmodel += repst(2);
        modelcounter += 1;
        
        fullmainmodel += " + ";
        fullmainmodel += repst(1);
        fullmainmodel += ":";
        fullmainmodel += size(2);
        modelcounter += 1;
        
        if (sizebused) {
          fullmainmodel += " + ";
          fullmainmodel += sizeb(2);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          modelcounter += 1;
          
          fullmainmodel += " + ";
          fullmainmodel += sizeb(1);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          modelcounter += 1;
          
          fullmainmodel += " + ";
          fullmainmodel += repst(1);
          fullmainmodel += ":";
          fullmainmodel += sizeb(2);
          modelcounter += 1;
        }
        
        if (sizecused) {
          fullmainmodel += " + ";
          fullmainmodel += sizec(2);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          modelcounter += 1;
          
          fullmainmodel += " + ";
          fullmainmodel += sizec(1);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          modelcounter += 1;
          
          fullmainmodel += " + ";
          fullmainmodel += repst(1);
          fullmainmodel += ":";
          fullmainmodel += sizec(2);
          modelcounter += 1;
        }
      }
      
      if (age != "none" && ageused) {
        fullmainmodel += " + ";
        fullmainmodel += age;
        fullmainmodel += ":";
        fullmainmodel += size(1);
        modelcounter += 1;
        
        if (sizebused) {
          fullmainmodel += " + ";
          fullmainmodel += age;
          fullmainmodel += ":";
          fullmainmodel += sizeb(1);
          modelcounter += 1;
        }
        
        if (sizecused) {
          fullmainmodel += " + ";
          fullmainmodel += age;
          fullmainmodel += ":";
          fullmainmodel += sizec(1);
          modelcounter += 1;
        }
        
        fullmainmodel += " + ";
        fullmainmodel += age;
        fullmainmodel += ":";
        fullmainmodel += repst(1);
        modelcounter += 1;
        
        if (densityused) {
          fullmainmodel += " + ";
          fullmainmodel += age;
          fullmainmodel += ":";
          fullmainmodel += densitycol;
          modelcounter += 1;
        }
        
        if (historical) {
          fullmainmodel += " + ";
          fullmainmodel += age;
          fullmainmodel += ":";
          fullmainmodel += size(2);
          modelcounter += 1;
          
          if (sizebused) {
            fullmainmodel += " + ";
            fullmainmodel += age;
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            modelcounter += 1;
          }
          if (sizecused) {
            fullmainmodel += " + ";
            fullmainmodel += age;
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            modelcounter += 1;
          }
          
          fullmainmodel += " + ";
          fullmainmodel += age;
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          modelcounter += 1;
        }
      }
      
      if (indcova(1) != "none" && !iaasrand && indcovaused) {
        fullmainmodel += " + ";
        fullmainmodel += indcova(1);
        fullmainmodel += ":";
        fullmainmodel += size(1);
        modelcounter += 1;
        
        if (sizebused) {
          fullmainmodel += " + ";
          fullmainmodel += indcova(1);
          fullmainmodel += ":";
          fullmainmodel += sizeb(1);
          modelcounter += 1;
        }
        
        if (sizecused) {
          fullmainmodel += " + ";
          fullmainmodel += indcova(1);
          fullmainmodel += ":";
          fullmainmodel += sizec(1);
          modelcounter += 1;
        }
        
        fullmainmodel += " + ";
        fullmainmodel += indcova(1);
        fullmainmodel += ":";
        fullmainmodel += repst(1);
        modelcounter += 1;
        
        if (densityused) {
          fullmainmodel += " + ";
          fullmainmodel += indcova(1);
          fullmainmodel += ":";
          fullmainmodel += densitycol;
          modelcounter += 1;
        }
        
        if (historical && indcova(2) != "none") {
          fullmainmodel += " + ";
          fullmainmodel += indcova(2);
          fullmainmodel += ":";
          fullmainmodel += size(2);
          modelcounter += 1;
          
          if (sizebused) {
            fullmainmodel += " + ";
            fullmainmodel += indcova(2);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcova(1);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcova(2);
            fullmainmodel += ":";
            fullmainmodel += sizeb(1);
            modelcounter += 1;
          }
          
          if (sizecused) {
            fullmainmodel += " + ";
            fullmainmodel += indcova(2);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcova(1);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcova(2);
            fullmainmodel += ":";
            fullmainmodel += sizec(1);
            modelcounter += 1;
          }
        
          fullmainmodel += " + ";
          fullmainmodel += indcova(2);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          modelcounter += 1;
          
          if (densityused) {
            fullmainmodel += " + ";
            fullmainmodel += indcova(2);
            fullmainmodel += ":";
            fullmainmodel += densitycol;
            modelcounter += 1;
          }
        }
      }
      if (indcovb(1) != "none" && !ibasrand && indcovbused) {
        fullmainmodel += " + ";
        fullmainmodel += indcovb(1);
        fullmainmodel += ":";
        fullmainmodel += size(1);
        modelcounter += 1;
        
        if (sizebused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovb(1);
          fullmainmodel += ":";
          fullmainmodel += sizeb(1);
          modelcounter += 1;
        }
        
        if (sizecused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovb(1);
          fullmainmodel += ":";
          fullmainmodel += sizec(1);
          modelcounter += 1;
        }
        
        fullmainmodel += " + ";
        fullmainmodel += indcovb(1);
        fullmainmodel += ":";
        fullmainmodel += repst(1);
        modelcounter += 1;
        
        if (densityused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovb(1);
          fullmainmodel += ":";
          fullmainmodel += densitycol;
          modelcounter += 1;
        }
        
        if (historical && indcovb(2) != "none") {
          fullmainmodel += " + ";
          fullmainmodel += indcovb(2);
          fullmainmodel += ":";
          fullmainmodel += size(2);
          modelcounter += 1;
          
          if (sizebused) {
            fullmainmodel += " + ";
            fullmainmodel += indcovb(2);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcovb(1);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcovb(2);
            fullmainmodel += ":";
            fullmainmodel += sizeb(1);
            modelcounter += 1;
          }
          
          if (sizecused) {
            fullmainmodel += " + ";
            fullmainmodel += indcovb(2);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcovb(1);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcovb(2);
            fullmainmodel += ":";
            fullmainmodel += sizec(1);
            modelcounter += 1;
          }
        
          fullmainmodel += " + ";
          fullmainmodel += indcovb(2);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          modelcounter += 1;
          
          if (densityused) {
            fullmainmodel += " + ";
            fullmainmodel += indcovb(2);
            fullmainmodel += ":";
            fullmainmodel += densitycol;
            modelcounter += 1;
          }
        }
      }
      
      if (indcovc(1) != "none" && !icasrand && indcovcused) {
        fullmainmodel += " + ";
        fullmainmodel += indcovc(1);
        fullmainmodel += ":";
        fullmainmodel += size(1);
        modelcounter += 1;
        
        if (sizebused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovc(1);
          fullmainmodel += ":";
          fullmainmodel += sizeb(1);
          modelcounter += 1;
        }
        
        if (sizecused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovc(1);
          fullmainmodel += ":";
          fullmainmodel += sizec(1);
          modelcounter += 1;
        }
        
        fullmainmodel += " + ";
        fullmainmodel += indcovc(1);
        fullmainmodel += ":";
        fullmainmodel += repst(1);
        modelcounter += 1;
        
        if (densityused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovc(1);
          fullmainmodel += ":";
          fullmainmodel += densitycol;
          modelcounter += 1;
        }
        
        if (historical && indcovc(2) != "none") {
          fullmainmodel += " + ";
          fullmainmodel += indcovc(2);
          fullmainmodel += ":";
          fullmainmodel += size(2);
          modelcounter += 1;
          
          if (sizebused) {
            fullmainmodel += " + ";
            fullmainmodel += indcovc(2);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcovc(1);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcovc(2);
            fullmainmodel += ":";
            fullmainmodel += sizeb(1);
            modelcounter += 1;
          }
          
          if (sizecused) {
            fullmainmodel += " + ";
            fullmainmodel += indcovc(2);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcovc(1);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            modelcounter += 1;
            
            fullmainmodel += " + ";
            fullmainmodel += indcovc(2);
            fullmainmodel += ":";
            fullmainmodel += sizec(1);
            modelcounter += 1;
          }
        
          fullmainmodel += " + ";
          fullmainmodel += indcovc(2);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          modelcounter += 1;
          
          if (densityused) {
            fullmainmodel += " + ";
            fullmainmodel += indcovc(2);
            fullmainmodel += ":";
            fullmainmodel += densitycol;
            modelcounter += 1;
          }
        }
      }
    }
    
    if (fixedcovcounter > 0) {
      if (modelcounter > 0) fullmainmodel += " + ";
      fullmainmodel += fixedtackon;
    }
    
    if (randomcovcounter > 0) {
      if (modelcounter > 0) fullmainmodel += " + ";
      fullmainmodel += randomtackon;
    }
    
  } else if (suite == "rep") {
    if (age != "none" || covcount > 0) fullmainmodel += " + ";
    fullmainmodel += repst(1);
    modelcounter += 1;
    
    if (historical && repstcheck) {
      fullmainmodel += " + ";
      fullmainmodel += repst(2);
      modelcounter += 1;
      
      fullmainmodel += " + ";
      fullmainmodel += repst(1);
      fullmainmodel += ":";
      fullmainmodel += repst(2);
      modelcounter += 1;
    }
    
    if (fixedcovcounter > 0) {
      if (modelcounter > 0) fullmainmodel += " + ";
      fullmainmodel += fixedtackon;
    }
    
    if (randomcovcounter > 0) {
      if (modelcounter > 0) fullmainmodel += " + ";
      fullmainmodel += randomtackon;
    }
    
  } else if (suite == "cons") {
    if (fixedcovcounter == 0 && modelcounter == 0) fullmainmodel += "1";
    modelcounter += 1;
    
    if (fixedcovcounter > 0) {
      if (modelcounter > 0) fullmainmodel += " + ";
      fullmainmodel += fixedtackon;
    }
    
    if (randomcovcounter > 0) {
      if (modelcounter > 0) fullmainmodel += " + ";
      fullmainmodel += randomtackon;
    }
    
  } else {
    fullmainmodel = "none";
  }
  
  // Build the global models
  if (response == 1) {
    String surv0_first(surv[0]);
    adult_out_model = surv0_first;
    adult_out_model += fullmainmodel;
    
    if (!nojuvs) {
      juv_out_model = surv0_first;
      juv_out_model += juvmainmodel;
      
    } else {
      juv_out_model = "none";
    }
    
  } else if (response == 2) {
    String obs0_first(obs[0]);
    adult_out_model = obs0_first;
    adult_out_model += fullmainmodel;
    
    if (!nojuvs) {
      juv_out_model = obs0_first;
      juv_out_model += juvmainmodel;
    } else {
      juv_out_model = "none";
    }
  } else if (response == 3) {
    String size0_first(size[0]);
    adult_out_model = size0_first;
    adult_out_model += fullmainmodel;
    
    if (!nojuvs) {
      juv_out_model = size0_first;
      juv_out_model += juvmainmodel;
    } else {
      juv_out_model = "none";
    }
  } else if (response == 4) {
    String sizeb0_first(sizeb[0]);
    adult_out_model = sizeb0_first;
    adult_out_model += fullmainmodel;
    
    if (!nojuvs) {
      juv_out_model = sizeb0_first;
      juv_out_model += juvmainmodel;
    } else {
      juv_out_model = "none";
    }
  } else if (response == 5) {
    String sizec0_first(sizec[0]);
    adult_out_model = sizec0_first;
    adult_out_model += fullmainmodel;
    
    if (!nojuvs) {
      juv_out_model = sizec0_first;
      juv_out_model += juvmainmodel;
    } else {
      juv_out_model = "none";
    }
  } else if (response == 6) {
    String repst0_first(repst[0]);
    adult_out_model = repst0_first;
    adult_out_model += fullmainmodel;
    
    if (!nojuvs) {
      juv_out_model = repst0_first;
      juv_out_model += juvmainmodel;
    } else {
      juv_out_model = "none";
    }
  } else if (response == 7) {
    if (fectime == 3) {
      String fec0_first(fec[0]);
      adult_out_model = fec0_first;
    } else {
      String fec1_first(fec[1]);
      adult_out_model = fec1_first;
    }
    adult_out_model += fullmainmodel;
  } else  if (response == 8) {
    if (!nojuvs) {
      String matstat0_first(matstat[0]);
      juv_out_model = matstat0_first;
      juv_out_model += juvmainmodel;
    } else {
      juv_out_model = "none";
    }
  }
  modelcounter += total_terms;
  
  List output(4);
  output(0) = adult_out_model;
  output(1) = juv_out_model;
  output(2) = modelcounter;
  output(3) = juvmodelcounter;
  
  CharacterVector out_names = {"adult_model", "juv_model", "total_terms",
    "juv_total_terms"};
  
  if (fullmainmodel == "none") {
    output["adult_model"] = 1;
    
    if (!nojuvs) {
      output["juv_model"] = 1;
    }
  }
  if (nojuvs) {
    output["juv_model"] = 0;
  }
  
  return output;
}

//' Main Formula Creation for Function \code{modelsearch()}
//'
//' Function \code{stovokor()} creates the list of formulae to be used as input
//' in the global model calls used in function \code{\link{modelsearch}()}.
//' 
//' @name stovokor
//' 
//' @param surv A vector of strings indicating the names of the variables coding
//' survival.
//' @param obs A vector of strings indicating the names of the variables coding
//' observation status.
//' @param size A vector of strings indicating the names of the variables coding
//' primary size.
//' @param sizeb A vector of strings indicating the names of the variables
//' coding secondary size.
//' @param sizec A vector of strings indicating the names of the variables
//' coding tertiary size.
//' @param repst A vector of strings indicating the names of the variables
//' coding reproductive status.
//' @param fec A vector of strings indicating the names of the variables coding
//' fecundity.
//' @param matstat A vector of strings indicating the names of the variables
//' coding for maturity status.
//' @param vitalrates A vector of strings indicating which vital rates will be
//' estimated.
//' @param historical A logical value indicating whether to create global models
//' with historical effects.
//' @param suite A string indicating the scope of independent factors included
//' in the global models. Options include \code{"full"}, \code{"main"},
//' \code{"size"}, \code{"rep"}, and \code{"const"}.
//' @param approach A string indicating whether to use mixed model encoding 
//' (\code{"mixed"}) or GLM encoding (\code{"glm"}).
//' @param nojuvs A logical value indicating that juvenile rates should be
//' estimated (\code{FALSE}) or not (\code{TRUE}).
//' @param juvsize A logical value indicating whether to include size terms in
//' juvenile models.
//' @param indiv A string indicating the name of the variable coding individual
//' identity.
//' @param patch A string indicating the name of the variable coding patch
//' identity.
//' @param year A string indicating the name of the variable coding time
//' \emph{t}.
//' @param age A string indicating the name of the variable coding age.
//' @param densitycol The name of the density variable, or \code{"none"}.
//' @param indcova A vector of strings indicating the names in times \emph{t}+1,
//' \emph{t}, and \emph{t}-1 of a specific individual covariate used in the
//' dataset.
//' @param indcovb A vector of strings indicating the names in times \emph{t}+1,
//' \emph{t}, and \emph{t}-1 of a specific individual covariate used in the
//' dataset.
//' @param indcovc A vector of strings indicating the names in times \emph{t}+1,
//' \emph{t}, and \emph{t}-1 of a specific individual covariate used in the
//' dataset.
//' @param sizebused A logical vector indicating if secondary size variables are
//' to be tested in each of the 14 vital rate models.
//' @param sizecused A logical vector indicating if tertiary size variables are
//' to be tested in each of the 14 vital rate models.
//' @param grouptest A logical vector indicating if group is to be tested in
//' each of the 14 vital rate models.
//' @param ageused A logical vector indicating if age is to be tested in each of
//' the 14 vital rate models.
//' @param densityused A logical vector indicating if density is to be tested in
//' each of the 14 vital rate models.
//' @param indcovaused A logical vector indicating if individual covariate a is
//' to be tested in each of the 14 vital rate models.
//' @param indcovbused A logical vector indicating if individual covariate b is
//' to be tested in each of the 14 vital rate models.
//' @param indcovcused A logical vector indicating if individual covariate c is
//' to be tested in each of the 14 vital rate models.
//' @param pasrand A logical value indicating whether to treat patch as a random
//' variable within mixed models.
//' @param yasrand A logical value indicating whether to treat year as a random
//' variable within mixed models.
//' @param iaasrand A logical value indicating whether to treat indcova as
//' random.
//' @param ibasrand A logical value indicating whether to treat indcovb as
//' random.
//' @param icasrand A logical value indicating whether to treat indcovc as
//' random.
//' @param iaasfac A logical value indicating whether to treat indcova as a
//' factor variable.
//' @param ibasfac A logical value indicating whether to treat indcovb as a
//' factor variable.
//' @param icasfac A logical value indicating whether to treat indcovc as a
//' factor variable.
//' @param fectime An integer indicating whether to use reproductive output in
//' time \emph{t} (2) or time \emph{t}+1 (3) as the response for fecundity.
//' @param size_zero A boolean variable indicating whether the primary size
//' model is zero-inflated.
//' @param sizeb_zero A boolean variable indicating whether the secondary size
//' model is zero-inflated.
//' @param sizec_zero A boolean variable indicating whether the tertiary size
//' model is zero-inflated.
//' @param jsize_zero A boolean variable indicating whether the juvenile primary
//' size model is zero-inflated.
//' @param jsizeb_zero A boolean variable indicating whether the juvenile
//' secondary size model is zero-inflated.
//' @param jsizec_zero A boolean variable indicating whether the juvenile
//' tertiary size model is zero-inflated.
//' 
//' @return A list of four lists. The first list includes the 14 main global
//' models covering all 14 vital rate responses, followed by an associated
//' \code{paramnames} object, and a vector of 14 integers showing the number of
//' terms tested in each respective model. The next three lists repeat this
//' structure without a new paramnames object, but for a reduced set of models
//' (\code{alternate}), a glm-only version of the models (\code{glm.alternate}),
//' and a version without any individual covariates (\code{nocovs.alternate}).
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export(.stovokor)]]
List stovokor(const StringVector& surv, const StringVector& obs,
  const StringVector& size, const StringVector& sizeb,
  const StringVector& sizec, const StringVector& repst, const StringVector& fec,
  const StringVector& matstat, const StringVector& vitalrates, bool historical,
  StringVector& suite, const String& approach, const bool nojuvs,
  bool juvsize, const String& indiv, const String& patch,
  const String& year, const String& age, const String& densitycol,
  const StringVector& indcova, const StringVector& indcovb,
  const StringVector& indcovc, const bool sizebused, const bool sizecused,
  const LogicalVector& grouptest, const LogicalVector& ageused,
  const LogicalVector& densityused, const LogicalVector& indcovaused,
  const LogicalVector& indcovbused, const LogicalVector& indcovcused,
  const bool pasrand, const bool yasrand, const bool iaasrand,
  const bool ibasrand, const bool icasrand, const bool iaasfac,
  const bool ibasfac, const bool icasfac, const int fectime,
  const bool size_zero, const bool sizeb_zero, const bool sizec_zero,
  const bool jsize_zero, const bool jsizeb_zero, const bool jsizec_zero) {
  
  if (nojuvs) juvsize = false;
  
  int nvitalrates = vitalrates.length();
  bool survcheck {false};
  bool obscheck {false};
  bool sizecheck {false};
  bool repstcheck {false};
  bool feccheck  {false};
  
  const int sizel {static_cast<int>(size.length())};
  const int repstl {static_cast<int>(repst.length())};
  
  String fullsurvmodel {"none"};
  String fullobsmodel {"none"};
  String fullsizemodel {"none"};
  String fullsizebmodel {"none"};
  String fullsizecmodel {"none"};
  String fullrepstmodel {"none"};
  String fullfecmodel {"none"};
  
  String juvsurvmodel {"none"};
  String juvobsmodel {"none"};
  String juvsizemodel {"none"};
  String juvsizebmodel {"none"};
  String juvsizecmodel {"none"};
  String juvrepstmodel {"none"};
  String juvmatstmodel {"none"};
  
  String alt_fullsurvmodel {"none"};
  String alt_fullobsmodel {"none"};
  String alt_fullsizemodel {"none"};
  String alt_fullsizebmodel {"none"};
  String alt_fullsizecmodel {"none"};
  String alt_fullrepstmodel {"none"};
  String alt_fullfecmodel {"none"};
  
  String alt_juvsurvmodel {"none"};
  String alt_juvobsmodel {"none"};
  String alt_juvsizemodel {"none"};
  String alt_juvsizebmodel {"none"};
  String alt_juvsizecmodel {"none"};
  String alt_juvrepstmodel {"none"};
  String alt_juvmatstmodel {"none"};
  
  String nocovs_fullsurvmodel {"none"};
  String nocovs_fullobsmodel {"none"};
  String nocovs_fullsizemodel {"none"};
  String nocovs_fullsizebmodel {"none"};
  String nocovs_fullsizecmodel {"none"};
  String nocovs_fullrepstmodel {"none"};
  String nocovs_fullfecmodel {"none"};
  
  String nocovs_juvsurvmodel {"none"};
  String nocovs_juvobsmodel {"none"};
  String nocovs_juvsizemodel {"none"};
  String nocovs_juvsizebmodel {"none"};
  String nocovs_juvsizecmodel {"none"};
  String nocovs_juvrepstmodel {"none"};
  String nocovs_juvmatstmodel {"none"};
  
  String glm_fullsurvmodel {"none"};
  String glm_fullobsmodel {"none"};
  String glm_fullsizemodel {"none"};
  String glm_fullsizebmodel {"none"};
  String glm_fullsizecmodel {"none"};
  String glm_fullrepstmodel {"none"};
  String glm_fullfecmodel {"none"};
  
  String glm_juvsurvmodel {"none"};
  String glm_juvobsmodel {"none"};
  String glm_juvsizemodel {"none"};
  String glm_juvsizebmodel {"none"};
  String glm_juvsizecmodel {"none"};
  String glm_juvrepstmodel {"none"};
  String glm_juvmatstmodel {"none"};
  
  IntegerVector total_terms (14);
  IntegerVector alt_total_terms (14);
  IntegerVector glm_total_terms (14);
  IntegerVector nocovs_total_terms (14);
  
  // Determine which vital rates need global model formulae
  for (int i = 0; i < nvitalrates; i++) {
    if (vitalrates(i) == "surv") survcheck = 1;
    if (vitalrates(i) == "obs") obscheck = 1;
    if (vitalrates(i) == "size") sizecheck = 1;
    if (vitalrates(i) == "repst") repstcheck = 1;
    if (vitalrates(i) == "fec") feccheck = 1;
  }
  
  String suite_element;
  
  for (int i = 0; i < suite.length(); i++) {
    suite_element = String(suite(i));
    // Tests if inputs are appropriate for suite
    if (suite_element == "full" || suite_element == "main") {
      if (historical) {
        if (sizel != 3 && repstl != 3) {
          if (sizel == 2 && repstl == 2) {
            historical = false;
          } else if (repstl == 3) {
            suite(i) = "rep";
          } else if (sizel == 3) {
            suite(i) = "size";
          }
        }
      } else {
        if (sizel != 2 && repstl != 2) {
          if (sizel != 3 && repstl != 3) {
            suite(i) = "const";
          } else if (repstl > 1 && sizel < 2) {
            suite(i) = "rep";
          } else if (sizel > 1 && repstl < 2) {
            suite(i) = "size";
          }
        }
      }
    } else if (suite_element == "rep") {
      if (historical) {
        if (repstl != 3) {
          if (repstl == 2) {
            historical = false;
          } else {
            suite(i) = "const";
          }
        }
      } else {
        if (repstl != 2) {
          if (repstl != 3) {
            suite(i) = "const";
          }
        }
      }
    } else if (suite_element == "size") {
      if (historical) {
        if (sizel != 3) {
          if (sizel == 2) {
            historical = false;
          } else {
            suite(i) = "const";
          }
        }
      } else {
        if (sizel != 2) {
          if (sizel != 3) {
            suite(i) = "const";
          }
        }
      }
    }
  }
  
  StringVector alt_suite (14);
  
  for (int i = 0; i < 14; i++) {
    if (suite(i) == "full") alt_suite(i) = "main";
    if (suite(i) == "main") alt_suite(i) = "size";
    if (suite(i) == "size") alt_suite(i) = "const";
    if (suite(i) == "rep") alt_suite(i) = "const";
    if (suite(i) == "const") alt_suite(i) = "const";
  }
  
  // Build global models
  if (survcheck) {
    List surv_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec, matstat,
      historical, 1, String(suite(0)), approach, nojuvs, juvsize, indiv, patch,
      year, age, densitycol, indcova, indcovb, indcovc, sizebused, sizecused,
      grouptest(0), ageused(0), densityused(0), indcovaused(0), indcovbused(0),
      indcovcused(0), pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac,
      ibasfac, icasfac, fectime, repstcheck);
    
    String sp0_proxy(as<StringVector>(surv_prax[0]));
    fullsurvmodel = sp0_proxy;
    total_terms(0) = static_cast<int>(surv_prax(2));
    
    List alt_surv_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec, matstat,
      historical, 1, String(alt_suite(0)), approach, nojuvs, juvsize, indiv, patch,
      year, age, densitycol, indcova, indcovb, indcovc, sizebused, sizecused,
      grouptest(0), ageused(0), densityused(0), indcovaused(0), indcovbused(0),
      indcovcused(0), pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac,
      ibasfac, icasfac, fectime, repstcheck);
    
    String alt_sp0_proxy(as<StringVector>(alt_surv_prax[0]));
    alt_fullsurvmodel = alt_sp0_proxy;
    alt_total_terms(0) = static_cast<int>(alt_surv_prax(2));
    
    if (!nojuvs) {
      String sp1_proxy(as<StringVector>(surv_prax[1]));
      juvsurvmodel = sp1_proxy;
      total_terms(7) = static_cast<int>(surv_prax(3));
      
      String alt_sp1_proxy(as<StringVector>(alt_surv_prax[1]));
      alt_juvsurvmodel = sp1_proxy;
      alt_total_terms(7) = static_cast<int>(alt_surv_prax(3));
      
      List matstat_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 8, String(suite(13)), approach, nojuvs, juvsize,
        indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(13), ageused(13), densityused(13),
        indcovaused(13), indcovbused(13), indcovcused(13), pasrand, yasrand,
        iaasrand, ibasrand, icasrand, iaasfac, ibasfac, icasfac, fectime,
        repstcheck);
        
      String mt1_proxy(as<StringVector>(matstat_prax[1]));
      juvmatstmodel = mt1_proxy;
      total_terms(13) = static_cast<int>(matstat_prax(3));
      
      List alt_matstat_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 8, String(alt_suite(13)), approach, nojuvs, juvsize,
        indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(13), ageused(13), densityused(13),
        indcovaused(13), indcovbused(13), indcovcused(13), pasrand, yasrand,
        iaasrand, ibasrand, icasrand, iaasfac, ibasfac, icasfac, fectime,
        repstcheck);
        
      String alt_mt1_proxy(as<StringVector>(alt_matstat_prax[1]));
      alt_juvmatstmodel = alt_mt1_proxy;
      alt_total_terms(13) = static_cast<int>(alt_matstat_prax(3));
    }
    
    if (approach == "mixed") {
      List glm_surv_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 1, String(suite(0)), "glm", nojuvs, juvsize, indiv,
        patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
        sizecused, grouptest(0), ageused(0), densityused(0), indcovaused(0),
        indcovbused(0), indcovcused(0), false, false, false, false, false,
        iaasfac, ibasfac, icasfac, fectime, repstcheck);
      
      String glm_sp0_proxy(as<StringVector>(glm_surv_prax[0]));
      glm_fullsurvmodel = glm_sp0_proxy;
      glm_total_terms(0) = static_cast<int>(glm_surv_prax(2));
      
      if (!nojuvs) {
        String glm_sp1_proxy(as<StringVector>(glm_surv_prax[1]));
        glm_juvsurvmodel = glm_sp1_proxy;
        glm_total_terms(7) = static_cast<int>(glm_surv_prax(3));
        
        List glm_matstat_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
          matstat, historical, 8, String(suite(13)), "glm", nojuvs, juvsize,
          indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
          sizebused, sizecused, grouptest(13), ageused(13), densityused(13),
          indcovaused(13), indcovbused(13), indcovcused(13), false, false,
          false, false, false, iaasfac, ibasfac, icasfac, fectime, repstcheck);
          
        String glm_mt1_proxy(as<StringVector>(glm_matstat_prax[1]));
        glm_juvmatstmodel = glm_mt1_proxy;
        glm_total_terms(13) = static_cast<int>(glm_matstat_prax(3));
      }
    }
    
    int surv_check_indivs = static_cast<int>(indcovaused(0)) + 
      static_cast<int>(indcovbused(0)) + static_cast<int>(indcovcused(0));
    
    if (surv_check_indivs > 0) { 
      List nocovs_surv_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 1, String(suite(0)), approach, nojuvs, juvsize,
        indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(0), ageused(0), densityused(0), false,
        false, false, pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac,
        ibasfac, icasfac, fectime, repstcheck);
      
      String nocovs_sp0_proxy(as<StringVector>(nocovs_surv_prax[0]));
      nocovs_fullsurvmodel = nocovs_sp0_proxy;
      nocovs_total_terms(0) = static_cast<int>(nocovs_surv_prax(2));
      
       if (!nojuvs) {
        String nocovs_sp1_proxy(as<StringVector>(nocovs_surv_prax[1]));
        nocovs_juvsurvmodel = nocovs_sp1_proxy;
        nocovs_total_terms(7) = static_cast<int>(nocovs_surv_prax(3));
      }
    }
    
    int matstat_check_indivs = static_cast<int>(indcovaused(13)) + 
      static_cast<int>(indcovbused(13)) + static_cast<int>(indcovcused(13));
    
    if (matstat_check_indivs > 0 && !nojuvs) { 
      List nocovs_matstat_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 8, String(suite(13)), approach, nojuvs, juvsize,
        indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(13), ageused(13), densityused(13),
        false, false, false, pasrand, yasrand, iaasrand, ibasrand, icasrand,
        iaasfac, ibasfac, icasfac, fectime, repstcheck);
        
      String nocovs_mt1_proxy(as<StringVector>(nocovs_matstat_prax[1]));
      nocovs_juvmatstmodel = nocovs_mt1_proxy;
      nocovs_total_terms(13) = static_cast<int>(nocovs_matstat_prax(3));
    }
  }
  
  if (obscheck) {
    List obs_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec, matstat,
      historical, 2, String(suite(1)), approach, nojuvs, juvsize, indiv, patch,
      year, age, densitycol, indcova, indcovb, indcovc, sizebused, sizecused,
      grouptest(1), ageused(1), densityused(1), indcovaused(1), indcovbused(1),
      indcovcused(1), pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac,
      ibasfac, icasfac, fectime, repstcheck);
    
    String ob0_proxy(as<StringVector>(obs_prax[0]));
    fullobsmodel = ob0_proxy;
    total_terms(1) = static_cast<int>(obs_prax(2));
    
    List alt_obs_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
      matstat, historical, 2, String(alt_suite(1)), approach, nojuvs, juvsize,
      indiv, patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
      sizecused, grouptest(1), ageused(1), densityused(1), indcovaused(1),
      indcovbused(1), indcovcused(1), pasrand, yasrand, iaasrand, ibasrand,
      icasrand, iaasfac, ibasfac, icasfac, fectime, repstcheck);
    
    String alt_ob0_proxy(as<StringVector>(alt_obs_prax[0]));
    alt_fullobsmodel = alt_ob0_proxy;
    alt_total_terms(1) = static_cast<int>(alt_obs_prax(2));
    
    if (!nojuvs) {
      String ob1_proxy(as<StringVector>(obs_prax[1]));
      juvobsmodel = ob1_proxy;
      total_terms(8) = static_cast<int>(obs_prax(3));
      
      String alt_ob1_proxy(as<StringVector>(alt_obs_prax[1]));
      alt_juvobsmodel = alt_ob1_proxy;
      alt_total_terms(8) = static_cast<int>(alt_obs_prax(3));
    }
    
    if (approach == "mixed") { 
      List glm_obs_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 2, String(suite(1)), "glm", nojuvs, juvsize, indiv,
        patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
        sizecused, grouptest(1), ageused(1), densityused(1), indcovaused(1),
        indcovbused(1), indcovcused(1), false, false, false, false, false,
        iaasfac, ibasfac, icasfac, fectime, repstcheck);
      
      String glm_ob0_proxy(as<StringVector>(glm_obs_prax[0]));
      glm_fullobsmodel = glm_ob0_proxy;
      glm_total_terms(1) = static_cast<int>(glm_obs_prax(2));
      
      if (!nojuvs) {
        String glm_ob1_proxy(as<StringVector>(glm_obs_prax[1]));
        glm_juvobsmodel = glm_ob1_proxy;
        glm_total_terms(8) = static_cast<int>(glm_obs_prax(3));
      }
    }
    
    int obs_check_indivs = static_cast<int>(indcovaused(1)) + 
      static_cast<int>(indcovbused(1)) + static_cast<int>(indcovcused(1));
    
    if (obs_check_indivs > 0) { 
      List nocovs_obs_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 2, String(suite(1)), approach, nojuvs, juvsize,
        indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(1), ageused(1), densityused(1), false,
        false, false, pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac,
        ibasfac, icasfac, fectime, repstcheck);
      
      String nocovs_ob0_proxy(as<StringVector>(nocovs_obs_prax[0]));
      nocovs_fullobsmodel = nocovs_ob0_proxy;
      nocovs_total_terms(1) = static_cast<int>(nocovs_obs_prax(2));
      
       if (!nojuvs) {
        String nocovs_ob1_proxy(as<StringVector>(nocovs_obs_prax[1]));
        nocovs_juvobsmodel = nocovs_ob1_proxy;
        nocovs_total_terms(8) = static_cast<int>(nocovs_obs_prax(3));
      }
    }
  }
  
  if (sizecheck) {
    String current_suite(suite(2));
    
    List size_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec, matstat,
      historical, 3, current_suite, approach, nojuvs, juvsize, indiv, patch,
      year, age, densitycol, indcova, indcovb, indcovc, sizebused, sizecused,
      grouptest(2), ageused(2), densityused(2), indcovaused(2), indcovbused(2),
      indcovcused(2), pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac,
      ibasfac, icasfac, fectime, repstcheck);
    
    total_terms(2) = static_cast<int>(size_prax(2));
    
    List alt_size_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
      matstat, historical, 3, String(alt_suite(2)), approach, nojuvs, juvsize,
      indiv, patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
      sizecused, grouptest(2), ageused(2), densityused(2), indcovaused(2),
      indcovbused(2), indcovcused(2), pasrand, yasrand, iaasrand, ibasrand,
      icasrand, iaasfac, ibasfac, icasfac, fectime, repstcheck);
    
    alt_total_terms(2) = static_cast<int>(alt_size_prax(2));
    
    if (stringcompare_hard(current_suite, "full") && total_terms(2) > 18 && size_zero) {
      suite(2) = "main";
      alt_suite(2) = "size";
      
      List size_prax_main = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 3, String(suite(2)), approach, nojuvs, juvsize,
        indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(2), ageused(2), densityused(2),
        indcovaused(2), indcovbused(2), indcovcused(2), pasrand, yasrand,
        iaasrand, ibasrand, icasrand, iaasfac, ibasfac, icasfac, fectime,
        repstcheck);
      
      String sz0_proxy(as<StringVector>(size_prax_main[0]));
      fullsizemodel = sz0_proxy;
      total_terms(2) = static_cast<int>(size_prax_main(2));
      
      List alt_size_prax_main = praxis(surv, obs, size, sizeb, sizec, repst,
        fec, matstat, historical, 3, String(alt_suite(2)), approach, nojuvs,
        juvsize, indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(2), ageused(2), densityused(2),
        indcovaused(2), indcovbused(2), indcovcused(2), pasrand, yasrand,
        iaasrand, ibasrand, icasrand, iaasfac, ibasfac, icasfac, fectime,
        repstcheck);
      
      String alt_sz0_proxy(as<StringVector>(alt_size_prax_main[0]));
      alt_fullsizemodel = alt_sz0_proxy;
      alt_total_terms(2) = static_cast<int>(alt_size_prax_main(2));
      
      if (!nojuvs && jsize_zero) {
        String sz1_proxy(as<StringVector>(size_prax_main[1]));
        juvsizemodel = sz1_proxy;
        total_terms(9) = static_cast<int>(size_prax_main(3));
        
        String alt_sz1_proxy(as<StringVector>(alt_size_prax_main[1]));
        alt_juvsizemodel = alt_sz1_proxy;
        alt_total_terms(9) = static_cast<int>(alt_size_prax_main(3));
      }
    } else {
      String sz0_proxy(as<StringVector>(size_prax[0]));
      fullsizemodel = sz0_proxy;
      
      String alt_sz0_proxy(as<StringVector>(alt_size_prax[0]));
      alt_fullsizemodel = alt_sz0_proxy;
      
      if (!nojuvs) {
        String sz1_proxy(as<StringVector>(size_prax[1]));
        juvsizemodel = sz1_proxy;
        total_terms(9) = static_cast<int>(size_prax(3));
        
        String alt_sz1_proxy(as<StringVector>(alt_size_prax[1]));
        alt_juvsizemodel = alt_sz1_proxy;
        alt_total_terms(9) = static_cast<int>(alt_size_prax(3));
      }
    }
    
    if (approach == "mixed") {
      List glm_size_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 3, current_suite, "glm", nojuvs, juvsize, indiv,
        patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
        sizecused, grouptest(2), ageused(2), densityused(2), indcovaused(2),
        indcovbused(2), indcovcused(2), false, false, false, false, false,
        iaasfac, ibasfac, icasfac, fectime, repstcheck);
      
      glm_total_terms(2) = static_cast<int>(glm_size_prax(2));
      
      if (!nojuvs) {
        String glm_sz1_proxy(as<StringVector>(glm_size_prax[1]));
        glm_juvsizemodel = glm_sz1_proxy;
        glm_total_terms(9) = static_cast<int>(glm_size_prax(3));
      }
    }
    
    int size_check_indivs = static_cast<int>(indcovaused(2)) + 
      static_cast<int>(indcovbused(2)) + static_cast<int>(indcovcused(2));
    
    if (size_check_indivs > 0) { 
      List nocovs_size_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 3, String(suite(2)), approach, nojuvs, juvsize,
        indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(2), ageused(2), densityused(2), false,
        false, false, pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac,
        ibasfac, icasfac, fectime, repstcheck);
      
      String nocovs_sz0_proxy(as<StringVector>(nocovs_size_prax[0]));
      nocovs_fullsizemodel = nocovs_sz0_proxy;
      nocovs_total_terms(2) = static_cast<int>(nocovs_size_prax(2));
      
       if (!nojuvs) {
        String nocovs_sz1_proxy(as<StringVector>(nocovs_size_prax[1]));
        nocovs_juvsizemodel = nocovs_sz1_proxy;
        nocovs_total_terms(9) = static_cast<int>(nocovs_size_prax(3));
      }
    }
    
    if (sizebused) { 
      String current_suite(suite(3));
      
      List sizeb_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 4, current_suite, approach, nojuvs, juvsize, indiv,
        patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
        sizecused, grouptest(3), ageused(3), densityused(3), indcovaused(3),
        indcovbused(3), indcovcused(3), pasrand, yasrand, iaasrand, ibasrand,
        icasrand, iaasfac, ibasfac, icasfac, fectime, repstcheck);
      
      total_terms(3) = static_cast<int>(sizeb_prax(2));
      
      List alt_sizeb_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 4, String(alt_suite(3)), approach, nojuvs, juvsize,
        indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(3), ageused(3), densityused(3),
        indcovaused(3), indcovbused(3), indcovcused(3), pasrand, yasrand,
        iaasrand, ibasrand, icasrand, iaasfac, ibasfac, icasfac, fectime,
        repstcheck);
      
      alt_total_terms(3) = static_cast<int>(alt_sizeb_prax(2));
      
      if (stringcompare_hard(current_suite, "full") && total_terms(3) > 18 && sizeb_zero) { 
        suite(3) = "main";
        alt_suite(3) = "size";
        
        List sizeb_prax_main = praxis(surv, obs, size, sizeb, sizec, repst, fec,
          matstat, historical, 4, String(suite(3)), approach, nojuvs, juvsize,
          indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
          sizebused, sizecused, grouptest(3), ageused(3), densityused(3),
          indcovaused(3), indcovbused(3), indcovcused(3), pasrand, yasrand,
          iaasrand, ibasrand, icasrand, iaasfac, ibasfac, icasfac, fectime,
          repstcheck);
        
        String szb0_proxy(as<StringVector>(sizeb_prax_main[0]));
        fullsizebmodel = szb0_proxy;
        total_terms(3) = static_cast<int>(sizeb_prax_main(2));
        
        List alt_sizeb_prax_main = praxis(surv, obs, size, sizeb, sizec, repst,
          fec, matstat, historical, 4, String(alt_suite(3)), approach, nojuvs,
          juvsize, indiv, patch, year, age, densitycol, indcova, indcovb,
          indcovc, sizebused, sizecused, grouptest(3), ageused(3),
          densityused(3), indcovaused(3), indcovbused(3), indcovcused(3),
          pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac, ibasfac,
          icasfac, fectime, repstcheck);
        
        String alt_szb0_proxy(as<StringVector>(alt_sizeb_prax_main[0]));
        alt_fullsizebmodel = alt_szb0_proxy;
        alt_total_terms(3) = static_cast<int>(alt_sizeb_prax_main(2));
        
        if (!nojuvs && jsizeb_zero) {
          String szb1_proxy(as<StringVector>(sizeb_prax_main[1]));
          juvsizebmodel = szb1_proxy;
          total_terms(10) = static_cast<int>(sizeb_prax_main(3));
          
          String alt_szb1_proxy(as<StringVector>(alt_sizeb_prax_main[1]));
          alt_juvsizebmodel = alt_szb1_proxy;
          alt_total_terms(10) = static_cast<int>(alt_sizeb_prax_main(3));
        }
      } else {
        String szb0_proxy(as<StringVector>(sizeb_prax[0]));
        fullsizebmodel = szb0_proxy;
        
        String alt_szb0_proxy(as<StringVector>(alt_sizeb_prax[0]));
        alt_fullsizebmodel = alt_szb0_proxy;
        
        if (!nojuvs) {
          String szb1_proxy(as<StringVector>(sizeb_prax[1]));
          juvsizebmodel = szb1_proxy;
          total_terms(10) = static_cast<int>(sizeb_prax(3));
          
          String alt_szb1_proxy(as<StringVector>(alt_sizeb_prax[1]));
          alt_juvsizebmodel = alt_szb1_proxy;
          alt_total_terms(10) = static_cast<int>(alt_sizeb_prax(3));
        }
      }
      
      if (approach == "mixed") {
        List glm_sizeb_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
          matstat, historical, 4, current_suite, "glm", nojuvs, juvsize, indiv,
          patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
          sizecused, grouptest(3), ageused(3), densityused(3), indcovaused(3),
          indcovbused(3), indcovcused(3), false, false, false, false, false,
          iaasfac, ibasfac, icasfac, fectime, repstcheck);
        
        glm_total_terms(3) = static_cast<int>(glm_sizeb_prax(3));
        
        if (!nojuvs) {
          String glm_szb1_proxy(as<StringVector>(glm_sizeb_prax[1]));
          glm_juvsizebmodel = glm_szb1_proxy;
          glm_total_terms(10) = static_cast<int>(glm_sizeb_prax(3));
        }
      }
      
      int sizeb_check_indivs = static_cast<int>(indcovaused(3)) + 
        static_cast<int>(indcovbused(3)) + static_cast<int>(indcovcused(3));
      
      if (sizeb_check_indivs > 0) { 
        List nocovs_sizeb_prax = praxis(surv, obs, size, sizeb, sizec, repst,
          fec, matstat, historical, 4, String(suite(3)), approach, nojuvs,
          juvsize, indiv, patch, year, age, densitycol, indcova, indcovb,
          indcovc, sizebused, sizecused, grouptest(3), ageused(3),
          densityused(3), false, false, false, pasrand, yasrand, iaasrand,
          ibasrand, icasrand, iaasfac, ibasfac, icasfac, fectime, repstcheck);
        
        String nocovs_szb0_proxy(as<StringVector>(nocovs_sizeb_prax[0]));
        nocovs_fullsizebmodel = nocovs_szb0_proxy;
        nocovs_total_terms(3) = static_cast<int>(nocovs_sizeb_prax(2));
        
         if (!nojuvs) {
          String nocovs_szb1_proxy(as<StringVector>(nocovs_sizeb_prax[1]));
          nocovs_juvsizebmodel = nocovs_szb1_proxy;
          nocovs_total_terms(10) = static_cast<int>(nocovs_sizeb_prax(3));
        }
      }
    }
    
    if (sizecused) { 
      String current_suite(suite(4));
        
      List sizec_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 5, current_suite, approach, nojuvs, juvsize, indiv,
        patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
        sizecused, grouptest(4), ageused(4), densityused(4), indcovaused(4),
        indcovbused(4), indcovcused(4), pasrand, yasrand, iaasrand, ibasrand,
        icasrand, iaasfac, ibasfac, icasfac, fectime, repstcheck);
      
      total_terms(4) = static_cast<int>(sizec_prax(2));
      
      List alt_sizec_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 5, String(alt_suite(4)), approach, nojuvs, juvsize,
        indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(4), ageused(4), densityused(4),
        indcovaused(4), indcovbused(4), indcovcused(4), pasrand, yasrand,
        iaasrand, ibasrand, icasrand, iaasfac, ibasfac, icasfac, fectime,
        repstcheck);
      
      alt_total_terms(4) = static_cast<int>(alt_sizec_prax(2));
      
      if (stringcompare_hard(current_suite, "full") && total_terms(4) > 18 && sizec_zero) {
        suite(4) = "main";
        alt_suite(4) = "size";
      
        List sizec_prax_main = praxis(surv, obs, size, sizeb, sizec, repst, fec,
          matstat, historical, 5, String(suite(4)), approach, nojuvs, juvsize,
          indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
          sizebused, sizecused, grouptest(4), ageused(4), densityused(4),
          indcovaused(4), indcovbused(4), indcovcused(4), pasrand, yasrand,
          iaasrand, ibasrand, icasrand, iaasfac, ibasfac, icasfac, fectime,
          repstcheck);
        
        String szc0_proxy(as<StringVector>(sizec_prax_main[0]));
        fullsizecmodel = szc0_proxy;
        total_terms(4) = static_cast<int>(sizec_prax_main(2));
        
        List alt_sizec_prax_main = praxis(surv, obs, size, sizeb, sizec, repst,
          fec, matstat, historical, 5, String(alt_suite(4)), approach, nojuvs,
          juvsize, indiv, patch, year, age, densitycol, indcova, indcovb,
          indcovc, sizebused, sizecused, grouptest(4), ageused(4),
          densityused(4), indcovaused(4), indcovbused(4), indcovcused(4),
          pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac, ibasfac,
          icasfac, fectime, repstcheck);
        
        String alt_szc0_proxy(as<StringVector>(alt_sizec_prax_main[0]));
        alt_fullsizecmodel = alt_szc0_proxy;
        alt_total_terms(4) = static_cast<int>(alt_sizec_prax_main(2));
        
        if (!nojuvs && jsizec_zero) {
          String szc1_proxy(as<StringVector>(sizec_prax_main[1]));
          juvsizecmodel = szc1_proxy;
          total_terms(11) = static_cast<int>(sizec_prax_main(3));
          
          String alt_szc1_proxy(as<StringVector>(alt_sizec_prax_main[1]));
          alt_juvsizecmodel = alt_szc1_proxy;
          alt_total_terms(11) = static_cast<int>(alt_sizec_prax_main(3));
        }
      } else {
        String szc0_proxy(as<StringVector>(sizec_prax[0]));
        fullsizecmodel = szc0_proxy;
        
        String alt_szc0_proxy(as<StringVector>(alt_sizec_prax[0]));
        alt_fullsizecmodel = alt_szc0_proxy;
        
        if (!nojuvs) {
          String szc1_proxy(as<StringVector>(sizec_prax[1]));
          juvsizecmodel = szc1_proxy;
          total_terms(11) = static_cast<int>(sizec_prax(3));
          
          String alt_szc1_proxy(as<StringVector>(alt_sizec_prax[1]));
          alt_juvsizecmodel = alt_szc1_proxy;
          alt_total_terms(11) = static_cast<int>(alt_sizec_prax(3));
        }
      }
      
      if (approach == "mixed") {
        List glm_sizec_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
          matstat, historical, 5, current_suite, "glm", nojuvs, juvsize, indiv,
          patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
          sizecused, grouptest(4), ageused(4), densityused(4), indcovaused(4),
          indcovbused(4), indcovcused(4), false, false, false, false, false,
          iaasfac, ibasfac, icasfac, fectime, repstcheck);
        
        glm_total_terms(4) = static_cast<int>(glm_sizec_prax(4));
        
        if (!nojuvs) {
          String glm_szc1_proxy(as<StringVector>(glm_sizec_prax[1]));
          glm_juvsizecmodel = glm_szc1_proxy;
          glm_total_terms(11) = static_cast<int>(glm_sizec_prax(3));
        }
      }
      
      int sizec_check_indivs = static_cast<int>(indcovaused(4)) + 
        static_cast<int>(indcovbused(4)) + static_cast<int>(indcovcused(4));
      
      if (sizec_check_indivs > 0) { 
        List nocovs_sizec_prax = praxis(surv, obs, size, sizeb, sizec, repst,
          fec, matstat, historical, 5, String(suite(4)), approach, nojuvs,
          juvsize, indiv, patch, year, age, densitycol, indcova, indcovb,
          indcovc, sizebused, sizecused, grouptest(4), ageused(4),
          densityused(4), false, false, false, pasrand, yasrand, iaasrand,
          ibasrand, icasrand, iaasfac, ibasfac, icasfac, fectime, repstcheck);
        
        String nocovs_szc0_proxy(as<StringVector>(nocovs_sizec_prax[0]));
        nocovs_fullsizecmodel = nocovs_szc0_proxy;
        nocovs_total_terms(4) = static_cast<int>(nocovs_sizec_prax(2));
        
         if (!nojuvs) {
          String nocovs_szc1_proxy(as<StringVector>(nocovs_sizec_prax[1]));
          nocovs_juvsizecmodel = nocovs_szc1_proxy;
          nocovs_total_terms(11) = static_cast<int>(nocovs_sizec_prax(3));
        }
      }
    }
  }
  
  if (repstcheck) {
    List repst_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec, matstat,
      historical, 6, String(suite(5)), approach, nojuvs, juvsize, indiv, patch,
      year, age, densitycol, indcova, indcovb, indcovc, sizebused, sizecused,
      grouptest(5), ageused(5), densityused(5), indcovaused(5), indcovbused(5),
      indcovcused(5), pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac,
      ibasfac, icasfac, fectime, repstcheck);
    
    String rp0_proxy(as<StringVector>(repst_prax[0]));
    fullrepstmodel = rp0_proxy;
    total_terms(5) = static_cast<int>(repst_prax(2));
    
    List alt_repst_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
      matstat, historical, 6, String(alt_suite(5)), approach, nojuvs, juvsize,
      indiv, patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
      sizecused, grouptest(5), ageused(5), densityused(5), indcovaused(5),
      indcovbused(5), indcovcused(5), pasrand, yasrand, iaasrand, ibasrand,
      icasrand, iaasfac, ibasfac, icasfac, fectime, repstcheck);
    
    String alt_rp0_proxy(as<StringVector>(alt_repst_prax[0]));
    alt_fullrepstmodel = alt_rp0_proxy;
    alt_total_terms(5) = static_cast<int>(alt_repst_prax(2));
    
    if (!nojuvs) {
      String rp1_proxy(as<StringVector>(repst_prax[1]));
      juvrepstmodel = rp1_proxy;
      total_terms(12) = static_cast<int>(repst_prax(3));
      
      String alt_rp1_proxy(as<StringVector>(alt_repst_prax[1]));
      alt_juvrepstmodel = alt_rp1_proxy;
      alt_total_terms(12) = static_cast<int>(alt_repst_prax(3));
    }
    
    if (approach == "mixed") { 
      List glm_repst_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 6, String(suite(5)), "glm", nojuvs, juvsize, indiv,
        patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
        sizecused, grouptest(5), ageused(5), densityused(5), indcovaused(5),
        indcovbused(5), indcovcused(5), false, false, false, false, false,
        iaasfac, ibasfac, icasfac, fectime, repstcheck);
      
      String glm_rp0_proxy(as<StringVector>(glm_repst_prax[0]));
      glm_fullrepstmodel = glm_rp0_proxy;
      glm_total_terms(5) = static_cast<int>(glm_repst_prax(2));
      
      if (!nojuvs) {
        String glm_rp1_proxy(as<StringVector>(glm_repst_prax[1]));
        glm_juvrepstmodel = glm_rp1_proxy;
        glm_total_terms(12) = static_cast<int>(glm_repst_prax(3));
      }
    }
    
    int repst_check_indivs = static_cast<int>(indcovaused(5)) + 
      static_cast<int>(indcovbused(5)) + static_cast<int>(indcovcused(5));
    
    if (repst_check_indivs > 0) { 
      List nocovs_repst_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 6, String(suite(5)), approach, nojuvs, juvsize,
        indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(5), ageused(5), densityused(5), false,
        false, false, pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac,
        ibasfac, icasfac, fectime, repstcheck);
      
      String nocovs_rp0_proxy(as<StringVector>(nocovs_repst_prax[0]));
      nocovs_fullrepstmodel = nocovs_rp0_proxy;
      nocovs_total_terms(5) = static_cast<int>(nocovs_repst_prax(2));
      
       if (!nojuvs) {
        String nocovs_rp1_proxy(as<StringVector>(nocovs_repst_prax[1]));
        nocovs_juvrepstmodel = nocovs_rp1_proxy;
        nocovs_total_terms(12) = static_cast<int>(nocovs_repst_prax(3));
      }
    }
  }
  
  if (feccheck) {
    List fec_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec, matstat,
      historical, 7, String(suite(6)), approach, nojuvs, juvsize, indiv, patch,
      year, age, densitycol, indcova, indcovb, indcovc, sizebused, sizecused,
      grouptest(6), ageused(6), densityused(6), indcovaused(6), indcovbused(6),
      indcovcused(6), pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac,
      ibasfac, icasfac, fectime, repstcheck);
    
    String fc0_proxy(as<StringVector>(fec_prax[0]));
    fullfecmodel = fc0_proxy;
    total_terms(6) = static_cast<int>(fec_prax(2));
    
    List alt_fec_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
      matstat, historical, 7, String(alt_suite(6)), approach, nojuvs, juvsize,
      indiv, patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
      sizecused, grouptest(6), ageused(6), densityused(6), indcovaused(6),
      indcovbused(6), indcovcused(6), pasrand, yasrand, iaasrand, ibasrand,
      icasrand, iaasfac, ibasfac, icasfac, fectime, repstcheck);
    
    String alt_fc0_proxy(as<StringVector>(alt_fec_prax[0]));
    alt_fullfecmodel = alt_fc0_proxy;
    alt_total_terms(6) = static_cast<int>(alt_fec_prax(2));
    
    if (approach == "mixed") { 
      List glm_fec_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 7, String(suite(6)), "glm", nojuvs, juvsize, indiv,
        patch, year, age, densitycol, indcova, indcovb, indcovc, sizebused,
        sizecused, grouptest(6), ageused(6), densityused(6), indcovaused(6),
        indcovbused(6), indcovcused(6), false, false, false, false, false,
        iaasfac, ibasfac, icasfac, fectime, repstcheck);
      
      String glm_fc0_proxy(as<StringVector>(glm_fec_prax[0]));
      glm_fullfecmodel = glm_fc0_proxy;
      glm_total_terms(6) = static_cast<int>(glm_fec_prax(2));
    }
    
    int fec_check_indivs = static_cast<int>(indcovaused(6)) + 
      static_cast<int>(indcovbused(6)) + static_cast<int>(indcovcused(6));
    
    if (fec_check_indivs > 0) { 
      List nocovs_fec_prax = praxis(surv, obs, size, sizeb, sizec, repst, fec,
        matstat, historical, 7, String(suite(6)), approach, nojuvs, juvsize,
        indiv, patch, year, age, densitycol, indcova, indcovb, indcovc,
        sizebused, sizecused, grouptest(6), ageused(6), densityused(6), false,
        false, false, pasrand, yasrand, iaasrand, ibasrand, icasrand, iaasfac,
        ibasfac, icasfac, fectime, repstcheck);
      
      String nocovs_fc0_proxy(as<StringVector>(nocovs_fec_prax[0]));
      nocovs_fullfecmodel = nocovs_fc0_proxy;
      nocovs_total_terms(6) = static_cast<int>(nocovs_fec_prax(2));
    }
  }
  
  // New paramnames object
  Rcpp::DataFrame paramnames = paramnames_skeleton(false);
  int pm_varno = static_cast<int>(paramnames.nrows());
  StringVector mainparams = as<StringVector>(paramnames["mainparams"]);
  
  StringVector modelparams (31);
  
  for (int i = 0; i < pm_varno; i++) {
    modelparams(i) = "none"; // Default value
    
    if (stringcompare_hard(as<std::string>(mainparams(i)), "year2")) modelparams(i) = year;
    if (stringcompare_hard(as<std::string>(mainparams(i)), "individ")) modelparams(i) = indiv;
    if (stringcompare_hard(as<std::string>(mainparams(i)), "patch")) modelparams(i) = patch;
    if (stringcompare_hard(as<std::string>(mainparams(i)), "surv3")) modelparams(i) = surv(0);
    if (stringcompare_hard(as<std::string>(mainparams(i)), "obs3")) modelparams(i) = obs(0);
    if (stringcompare_hard(as<std::string>(mainparams(i)), "size3")) modelparams(i) = size(0);
    if (stringcompare_hard(as<std::string>(mainparams(i)), "repst3")) modelparams(i) = repst(0);
    if (fectime == 3 && stringcompare_hard(as<std::string>(mainparams(i)), "fec3")) modelparams(i) = fec(0);
    if (fectime == 2 && stringcompare_hard(as<std::string>(mainparams(i)), "fec2")) modelparams(i) = fec(1);
    if (stringcompare_hard(as<std::string>(mainparams(i)), "size2")) modelparams(i) = size(1);
    
    if (sizebused) {
      if (stringcompare_hard(as<std::string>(mainparams(i)), "sizeb3")) modelparams(i) = sizeb(0);
      if (stringcompare_hard(as<std::string>(mainparams(i)), "sizeb2")) modelparams(i) = sizeb(1);
      if (historical && stringcompare_hard(as<std::string>(mainparams(i)), "sizeb1")) modelparams(i) = sizeb(2);
      
    } else {
      if (stringcompare_hard(as<std::string>(mainparams(i)), "sizeb3")) modelparams(i) = "none";
      if (stringcompare_hard(as<std::string>(mainparams(i)), "sizeb2")) modelparams(i) = "none";
      if (historical && stringcompare_hard(as<std::string>(mainparams(i)), "sizeb1")) modelparams(i) = "none";
    }
    
    if (sizecused) {
      if (stringcompare_hard(as<std::string>(mainparams(i)), "sizec3")) modelparams(i) = sizec(0);
      if (stringcompare_hard(as<std::string>(mainparams(i)), "sizec2")) modelparams(i) = sizec(1);
      if (historical && stringcompare_hard(as<std::string>(mainparams(i)), "sizec1")) modelparams(i) = sizec(2);
      
    } else {
      if (stringcompare_hard(as<std::string>(mainparams(i)), "sizec3")) modelparams(i) = "none";
      if (stringcompare_hard(as<std::string>(mainparams(i)), "sizec2")) modelparams(i) = "none";
      if (historical && stringcompare_hard(as<std::string>(mainparams(i)), "sizec1")) modelparams(i) = "none";
    }
    
    if (stringcompare_hard(as<std::string>(mainparams(i)), "repst2")) modelparams(i) = repst(1);
    if (stringcompare_hard(as<std::string>(mainparams(i)), "matst3")) modelparams(i) = matstat(0);
    if (stringcompare_hard(as<std::string>(mainparams(i)), "matstat2")) modelparams(i) = matstat(1);
    if (stringcompare_hard(as<std::string>(mainparams(i)), "age")) modelparams(i) = age;
    if (densityused && stringcompare_hard(as<std::string>(mainparams(i)), "density")) modelparams(i) = densitycol;
    if (grouptest && stringcompare_hard(as<std::string>(mainparams(i)), "group2")) modelparams(i) = "group2";
    if (indcovaused && stringcompare_hard(as<std::string>(mainparams(i)), "indcova2")) modelparams(i) = indcova(1);
    if (indcovbused && stringcompare_hard(as<std::string>(mainparams(i)), "indcovb2")) modelparams(i) = indcovb(1);
    if (indcovcused && stringcompare_hard(as<std::string>(mainparams(i)), "indcovc2")) modelparams(i) = indcovc(1);
    
    // Further historical terms, if used
    if (historical) {
      if (stringcompare_hard(as<std::string>(mainparams(i)), "size1")) modelparams(i) = size(2);
      if (stringcompare_hard(as<std::string>(mainparams(i)), "repst1")) modelparams(i) = repst(2);
      if (stringcompare_hard(as<std::string>(mainparams(i)), "group1")) modelparams(i) = "group1";
      
      if (indcovaused && stringcompare_hard(as<std::string>(mainparams(i)), "indcova1")) modelparams(i) = indcova(2);
      if (indcovbused && stringcompare_hard(as<std::string>(mainparams(i)), "indcovb1")) modelparams(i) = indcovb(2);
      if (indcovcused && stringcompare_hard(as<std::string>(mainparams(i)), "indcovc1")) modelparams(i) = indcovc(2);
    }
  }
  paramnames["modelparams"] = modelparams;
  
  Rcpp::List main = List::create(Named("full.surv.model") = fullsurvmodel,
    _["full.obs.model"] = fullobsmodel, _["full.size.model"] = fullsizemodel,
    _["full.sizeb.model"] = fullsizebmodel, _["full.sizec.model"] = fullsizecmodel,
    _["full.repst.model"] = fullrepstmodel, _["full.fec.model"] = fullfecmodel,
    _["juv.surv.model"] = juvsurvmodel, _["juv.obs.model"] = juvobsmodel,
    _["juv.size.model"] = juvsizemodel, _["juv.sizeb.model"] = juvsizebmodel,
    _["juv.sizec.model"] = juvsizecmodel, _["juv.repst.model"] = juvrepstmodel,
    _["juv.matst.model"] = juvmatstmodel, _["paramnames"] = paramnames,
    _["total_terms"] = total_terms);
  
  Rcpp::List alternate = List::create(Named("full.surv.model") = alt_fullsurvmodel,
    _["full.obs.model"] = alt_fullobsmodel, _["full.size.model"] = alt_fullsizemodel,
    _["full.sizeb.model"] = alt_fullsizebmodel, _["full.sizec.model"] = alt_fullsizecmodel,
    _["full.repst.model"] = alt_fullrepstmodel, _["full.fec.model"] = alt_fullfecmodel,
    _["juv.surv.model"] = alt_juvsurvmodel, _["juv.obs.model"] = alt_juvobsmodel,
    _["juv.size.model"] = alt_juvsizemodel, _["juv.sizeb.model"] = alt_juvsizebmodel,
    _["juv.sizec.model"] = alt_juvsizecmodel, _["juv.repst.model"] = alt_juvrepstmodel,
    _["juv.matst.model"] = alt_juvmatstmodel, _["total_terms"] = alt_total_terms);
  
  Rcpp::List glm_alternate = List::create(Named("full.surv.model") = glm_fullsurvmodel,
    _["full.obs.model"] = glm_fullobsmodel, _["full.size.model"] = glm_fullsizemodel,
    _["full.sizeb.model"] = glm_fullsizebmodel, _["full.sizec.model"] = glm_fullsizecmodel,
    _["full.repst.model"] = glm_fullrepstmodel, _["full.fec.model"] = glm_fullfecmodel,
    _["juv.surv.model"] = glm_juvsurvmodel, _["juv.obs.model"] = glm_juvobsmodel,
    _["juv.size.model"] = glm_juvsizemodel, _["juv.sizeb.model"] = glm_juvsizebmodel,
    _["juv.sizec.model"] = glm_juvsizecmodel, _["juv.repst.model"] = glm_juvrepstmodel,
    _["juv.matst.model"] = glm_juvmatstmodel, _["total_terms"] = glm_total_terms);
  
  Rcpp::List nocovs_alternate = List::create(Named("full.surv.model") = nocovs_fullsurvmodel,
    _["full.obs.model"] = nocovs_fullobsmodel, _["full.size.model"] = nocovs_fullsizemodel,
    _["full.sizeb.model"] = nocovs_fullsizebmodel, _["full.sizec.model"] = nocovs_fullsizecmodel,
    _["full.repst.model"] = nocovs_fullrepstmodel, _["full.fec.model"] = nocovs_fullfecmodel,
    _["juv.surv.model"] = nocovs_juvsurvmodel, _["juv.obs.model"] = nocovs_juvobsmodel,
    _["juv.size.model"] = nocovs_juvsizemodel, _["juv.sizeb.model"] = nocovs_juvsizebmodel,
    _["juv.sizec.model"] = nocovs_juvsizecmodel, _["juv.repst.model"] = nocovs_juvrepstmodel,
    _["juv.matst.model"] = nocovs_juvmatstmodel, _["total_terms"] = nocovs_total_terms);
  
  List output = List::create(_["main"] = main, _["alternate"] = alternate,
    _["glm.alternate"] = glm_alternate, _["nocovs.alternate"] = nocovs_alternate);
  
  return output;
}

//' Creates a Skeleton Paramnames Object for Use in Function-based Modeling
//' 
//' Creates a simple skeleton \code{paramnames} object that can be entered as
//' input in functions \code{\link{flefko2}()}, \code{\link{flefko3}()}, and
//' \code{\link{aflefko2}()}.
//' 
//' @name create_pm
//' 
//' @param name_terms A logical value indicating whether to start each variable
//' name as \code{none} if \code{FALSE}, or as the default \code{modelparams}
//' name if \code{TRUE}. Defaults to \code{FALSE}.
//' 
//' @return A three column data frame, of which the first describes the
//' parameters in reasonably plain English, the second gives the name of the
//' parameter within the MPM generating functions, and the third is to be
//' edited with the names of the variables as they appear in the models.
//' 
//' @section Notes:
//' The third column in the resulting object should be edited with the names only
//' of those variables actually used in vital rate modeling. This
//' \code{paramnames} object should apply to all models used in a single MPM
//' building exercise. So, for example, if the models used include random terms,
//' then they should all have the same random terms. Fixed terms can vary,
//' however.
//' 
//' @examples 
//' our_pm <- create_pm()
//' our_pm
//' 
//' @export
// [[Rcpp::export(create_pm)]]
DataFrame create_pm(bool name_terms = false) {
  
  DataFrame output = paramnames_skeleton(name_terms);
    
  return output;
}

