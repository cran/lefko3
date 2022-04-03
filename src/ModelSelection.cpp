#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Main Formula Creation for Function \code{modelsearch()}
//'
//' Function \code{.stovokor()} creates formulae to be used as input in the
//' global model calls used in function \code{\link{modelsearch}()}.
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
//' @param age A string indicating the name of the variable coding age.
//' @param indcova A vector of strings indicating the names in times \emph{t}+1,
//' \emph{t}, and \emph{t}-1 of a specific individual covariate used in the
//' dataset.
//' @param indcovb A vector of strings indicating the names in times \emph{t}+1,
//' \emph{t}, and \emph{t}-1 of a specific individual covariate used in the
//' dataset.
//' @param indcovc A vector of strings indicating the names in times \emph{t}+1,
//' \emph{t}, and \emph{t}-1 of a specific individual covariate used in the
//' dataset.
//' @param indiv A string indicating the name of the variable coding individual
//' identity.
//' @param patch A string indicating the name of the variable coding patch
//' identity.
//' @param year A string indicating the name of the variable coding time
//' \emph{t}.
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
//' @param fectime An integer indicating whether to use reproductive output in
//' time \emph{t} (2) or time \emph{t}+1 (3) as the response for fecundity.
//' @param juvsize A logical value indicating whether to include size terms in
//' juvenile models.
//' @param sizebused A logical value denoting if secondary size variables are to
//' be used.
//' @param sizecused A logical value denoting if tertiary size variables are to
//' be used.
//' @param grouptest A logical value indicating whether to test for group
//' effect.
//' @param densitycol The name of the density variable, or \code{"none"}.
//' @param densityused A logical value indicating whether the density variable
//' is to be used.
//' @param indcovaused Logical value indicating whether individual covariate a
//' is used.
//' @param indcovbused Logical value indicating whether individual covariate b
//' is used.
//' @param indcovcused Logical value indicating whether individual covariate c
//' is used.
//' 
//' @return Vector of 9 strings, each a formula to be used as input in function.
//' \code{modelsearch()}.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export(.stovokor)]]
List stovokor(StringVector surv, StringVector obs, StringVector size,
  StringVector sizeb, StringVector sizec, StringVector repst, StringVector fec,
  StringVector matstat, StringVector vitalrates, bool historical, String suite,
  String approach, bool nojuvs, String age, StringVector indcova,
  StringVector indcovb, StringVector indcovc, String indiv, String patch,
  String year, bool pasrand, bool yasrand, bool iaasrand, bool ibasrand,
  bool icasrand, int fectime, bool juvsize, bool sizebused, bool sizecused,
  bool grouptest, String densitycol, bool densityused, bool indcovaused,
  bool indcovbused, bool indcovcused) {
  
  if (nojuvs) juvsize = false;
  
  int nvitalrates = vitalrates.length();
  bool survcheck = 0;
  bool obscheck = 0;
  bool sizecheck = 0;
  bool repstcheck = 0;
  bool feccheck = 0;
  
  int sizel = size.length();
  int repstl = repst.length();
  
  String fullsurvmodel = "none";
  String fullobsmodel = "none";
  String fullsizemodel = "none";
  String fullsizebmodel = "none";
  String fullsizecmodel = "none";
  String fullrepstmodel = "none";
  String fullfecmodel = "none";
  
  String juvsurvmodel = "none";
  String juvobsmodel = "none";
  String juvsizemodel = "none";
  String juvsizebmodel = "none";
  String juvsizecmodel = "none";
  String juvrepstmodel = "none";
  String juvmatstmodel = "none";
  
  String randomtackonp = "";
  String randomtackony = "";
  String randomtackoni = "";
  String randomtackonia = "";
  String randomtackonib = "";
  String randomtackonic = "";
  String randomtackon = "";
  
  String fixedtackong = "";
  String fixedtackonp = "";
  String fixedtackony = "";
  String fixedtackonia = "";
  String fixedtackonib = "";
  String fixedtackonic = "";
  String fixedtackon = "";
  
  String sizesuffix;
  String jsizesuffix;
  String fecsuffix;
  
  String fullmainmodel;
  String juvmainmodel;
  
  int covcount {0};
  if (indcova(1) != "none") covcount += 1;
  if (indcovb(1) != "none") covcount += 1;
  if (indcovc(1) != "none") covcount += 1;
  
  int modelcounter {0};
  int juvmodelcounter {0};
  int fixedcovcounter {0};
  int randomcovcounter {0};
  
  int indcova_ran_counter {0};
  int indcovb_ran_counter {0};
  int indcovc_ran_counter {0};
  int year_ran_counter {0};
  int patch_ran_counter {0};
  int indiv_ran_counter {0};
  
  int indcova_fix_counter {0};
  int indcovb_fix_counter {0};
  int indcovc_fix_counter {0};
  int year_fix_counter {0};
  int patch_fix_counter {0};
  int group_fix_counter {0};
  
  // This section determines which vital rates need global model formulae
  for (int i = 0; i < nvitalrates; i++) {
    if (vitalrates(i) == "surv") survcheck = 1;
    if (vitalrates(i) == "obs") obscheck = 1;
    if (vitalrates(i) == "size") sizecheck = 1;
    if (vitalrates(i) == "repst") repstcheck = 1;
    if (vitalrates(i) == "fec") feccheck = 1;
  }
  
  // This section tests to see if the inputs are appropriate for the suite
  if (suite == "full" || suite == "main") {
    if (historical) {
      if (sizel != 3 && repstl != 3) {
        if (sizel == 2 && repstl == 2) {
          historical = false;
        } else if (repstl == 3) {
          suite = "rep";
        } else if (sizel == 3) {
          suite = "size";
        }
      }
    } else {
      if (sizel != 2 && repstl != 2) {
        if (sizel != 3 && repstl != 3) {
          suite = "const";
        } else if (repstl > 1 && sizel < 2) {
          suite = "rep";
        } else if (sizel > 1 && repstl < 2) {
          suite = "size";
        }
      }
    }
  } else if (suite == "rep") {
    if (historical) {
      if (repstl != 3) {
        if (repstl == 2) {
          historical = false;
        } else {
          suite = "const";
        }
      }
    } else {
      if (repstl != 2) {
        if (repstl != 3) {
          suite = "const";
        }
      }
    }
  } else if (suite == "size") {
    if (historical) {
      if (sizel != 3) {
        if (sizel == 2) {
          historical = false;
        } else {
          suite = "const";
        }
      }
    } else {
      if (sizel != 2) {
        if (sizel != 3) {
          suite = "const";
        }
      }
    }
  }
  
  // Here we determine the nature of the potentially random variables
  if (approach != "mixed") {
    yasrand = false;
    pasrand = false;
  }
  
  if (grouptest) {
    fixedcovcounter += 1;
    fixedtackong += "as.factor(group2)";
    
    if (historical) {
      fixedcovcounter += 1;
      fixedtackong += "+ as.factor(group1)";
    }
    group_fix_counter += 1;
  }
  
  if (year!= "none") {
    if (yasrand) {
      randomtackony += "(1 | ";
      randomtackony += year;
      randomtackony += ")";
      
      year_ran_counter += 1;
      randomcovcounter += 1;
    } else {
      fixedtackony += "as.factor(";
      fixedtackony += year;
      fixedtackony += ")";
      
      year_fix_counter += 1;
      fixedcovcounter += 1;
    }
  }
  
  if (patch!= "none") {
    if (pasrand) {
      randomtackonp += "(1 | ";
      randomtackonp += patch;
      randomtackonp += ")";
      
      patch_ran_counter += 1;
      randomcovcounter += 1;
    } else {
      fixedtackonp += "as.factor(";
      fixedtackonp += patch;
      fixedtackonp += ")";
      
      patch_fix_counter += 1;
      fixedcovcounter += 1;
    }
  }
  
  if (indiv != "none" && approach == "mixed") {
    randomtackoni += "(1 | ";
    randomtackoni += indiv;
    randomtackoni += ")";
    
    indiv_ran_counter += 1;
    randomcovcounter += 1;
  }
  
  // Now we add the individual covariates to the tacked-on sections
  if (indcova(1) != "none") {
    if (!iaasrand) {
      fixedtackonia += indcova(1);
      if (historical) {
        fixedtackonia += " + ";
        fixedtackonia += indcova(2);
      }
      
      indcova_fix_counter += 1;
      fixedcovcounter += 1;
    } else {
      randomtackonia += "(1 | ";
      randomtackonia += indcova(1);
      randomtackonia += ")";
      
      if (historical) {
        randomtackonia += " + (1 | ";
        randomtackonia += indcova(2);
        randomtackonia += ") + ";
      }
      
      indcova_ran_counter += 1;
      randomcovcounter += 1;
    }
  }
  if (indcovb(1) != "none") {
    if (!ibasrand) {
      fixedtackonib += indcovb(1);
      
      if (historical) {
        fixedtackonib += " + ";
        fixedtackonib += indcovb(2);
      }
      
      indcovb_fix_counter += 1;
      fixedcovcounter += 1;
    } else {
      randomtackonib += "(1 | ";
      randomtackonib += indcovb(1);
      randomtackonib += ")";
      
      if (historical) {
        randomtackonib += " + (1 | ";
        randomtackonib += indcovb(2);
        randomtackonib += ")";
      }
      
      indcovb_ran_counter += 1;
      randomcovcounter += 1;
    }
  }
  if (indcovc(1) != "none") {
    if (!icasrand) {
      fixedtackonic += indcovc(1);
      if (historical) {
        fixedtackonic += " + ";
        fixedtackonic += indcovc(2);
      }
      
      indcovc_fix_counter += 1;
      fixedcovcounter += 1;
    } else {
      randomtackonic += "(1 | ";
      randomtackonic += indcovc(1);
      randomtackonic += ")";
      
      if (historical) {
        randomtackonic += " + (1 | ";
        randomtackonic += indcovc(2);
        randomtackonic += ")";
      }
      
      indcovc_ran_counter += 1;
      randomcovcounter += 1;
    }
  }
  if (suite == "full" && !iaasrand && !ibasrand) {
    if (indcova(1) != "none" && indcovb(1) != "none") {
      fixedtackonib += " + ";
      fixedtackonib += indcova(1);
      fixedtackonib += ":";
      fixedtackonib += indcovb(1);
      
      if (historical) {
        fixedtackonib += " + ";
        fixedtackonib += indcova(2);
        fixedtackonib += ":";
        fixedtackonib += indcovb(2);
        
        fixedtackonib += " + ";
        fixedtackonib += indcova(1);
        fixedtackonib += ":";
        fixedtackonib += indcovb(2);
        
        fixedtackonib += " + ";
        fixedtackonib += indcova(2);
        fixedtackonib += ":";
        fixedtackonib += indcovb(1);
      }
    }
  }
  if (suite == "full" && !iaasrand && !icasrand) {
    if (indcova(1) != "none" && indcovc(1) != "none") {
      fixedtackonic += " + ";
      fixedtackonic += indcova(1);
      fixedtackonic += ":";
      fixedtackonic += indcovc(1);
      
      if (historical) {
        fixedtackonic += " + ";
        fixedtackonic += indcova(2);
        fixedtackonic += ":";
        fixedtackonic += indcovc(2);
        
        fixedtackonic += " + ";
        fixedtackonic += indcova(1);
        fixedtackonic += ":";
        fixedtackonic += indcovc(2);
        
        fixedtackonic += " + ";
        fixedtackonic += indcova(2);
        fixedtackonic += ":";
        fixedtackonic += indcovc(1);
      }
    }
  }
  if (suite == "full" && !ibasrand && !icasrand) {
    if (indcovb(1) != "none" && indcovc(1) != "none") {
      fixedtackonic += " + ";
      fixedtackonic += indcovb(1);
      fixedtackonic += ":";
      fixedtackonic += indcovc(1);
      
      if (historical) {
        fixedtackonic += " + ";
        fixedtackonic += indcovb(2);
        fixedtackonic += ":";
        fixedtackonic += indcovc(2);
        
        fixedtackonic += " + ";
        fixedtackonic += indcovb(1);
        fixedtackonic += ":";
        fixedtackonic += indcovc(2);
        
        fixedtackonic += " + ";
        fixedtackonic += indcovb(2);
        fixedtackonic += ":";
        fixedtackonic += indcovc(1);
      }
    }
  }
  
  if (indcova_ran_counter > 0) {
    randomtackon += " + ";
    randomtackon += randomtackonia;
  }
  if (indcovb_ran_counter > 0) {
    randomtackon += " + ";
    randomtackon += randomtackonib;
  }
  if (indcovc_ran_counter > 0) {
    randomtackon += " + ";
    randomtackon += randomtackonic;
  }
  if (year_ran_counter > 0) {
    randomtackon += " + ";
    randomtackon += randomtackony;
  }
  if (patch_ran_counter > 0) {
    randomtackon += " + ";
    randomtackon += randomtackonp;
  }
  if (indiv_ran_counter > 0) {
    randomtackon += " + ";
    randomtackon += randomtackoni;
  }
  
  if (group_fix_counter > 0) {
    fixedtackon += " + ";
    fixedtackon += fixedtackong;
  }
  if (indcova_fix_counter > 0) {
    fixedtackon += " + ";
    fixedtackon += fixedtackonia;
  }
  if (indcovb_fix_counter > 0) {
    fixedtackon += " + ";
    fixedtackon += fixedtackonib;
  }
  if (indcovc_fix_counter > 0) {
    fixedtackon += " + ";
    fixedtackon += fixedtackonic;
  }
  if (year_fix_counter > 0) {
    fixedtackon += " + ";
    fixedtackon += fixedtackony;
  }
  if (patch_fix_counter > 0) {
    fixedtackon += " + ";
    fixedtackon += fixedtackonp;
  }
  
  // Main model patterns
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
            
            if (sizebused) {
              juvmainmodel += " + ";
              juvmainmodel += sizeb(1);
              juvmainmodel += ":";
              juvmainmodel += sizec(1);
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
            
            if (sizebused) {
              juvmainmodel += " + ";
              juvmainmodel += sizeb(1);
              juvmainmodel += ":";
              juvmainmodel += densitycol;
            }
            if (sizecused) {
              juvmainmodel += " + ";
              juvmainmodel += sizec(1);
              juvmainmodel += ":";
              juvmainmodel += densitycol;
            }
          }
        }
      } else if (densityused) {
        if (juvmodelcounter > 0) juvmainmodel += " + ";
        juvmainmodel += densitycol;
        juvmodelcounter += 1;
      } else if (fixedcovcounter == 0) {
        juvmainmodel += "1";
      }
    } else  if (fixedcovcounter == 0) {
      juvmainmodel += "1";
    }
    
    if (fixedcovcounter > 0) juvmainmodel += fixedtackon;
    if (randomcovcounter > 0) juvmainmodel += randomtackon;
  }
    
  // Now the adult model pattern
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
      }
      
      if (historical) {
        fullmainmodel += " + ";
        fullmainmodel += size(2);
        
        if (suite != "main") {
          fullmainmodel += " + ";
          fullmainmodel += size(1);
          fullmainmodel += ":";
          fullmainmodel += size(2);
        }
        
        if (suite == "full" && densityused) {
          fullmainmodel += " + ";
          fullmainmodel += size(2);
          fullmainmodel += ":";
          fullmainmodel += densitycol;
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
        }
        
        if (historical) {
          fullmainmodel += " + ";
          fullmainmodel += sizeb(2);
        
          if (suite != "main") {
            fullmainmodel += " + ";
            fullmainmodel += sizeb(1);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            
            fullmainmodel += " + ";
            fullmainmodel += size(1);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            
            fullmainmodel += " + ";
            fullmainmodel += sizeb(1);
            fullmainmodel += ":";
            fullmainmodel += size(2);
          }
          
          if (suite == "full" && densityused) {
            fullmainmodel += " + ";
            fullmainmodel += sizeb(2);
            fullmainmodel += ":";
            fullmainmodel += densitycol;
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
        }
        
        if (historical) {
          fullmainmodel += " + ";
          fullmainmodel += sizec(2);
        
          if (suite != "main") {
            fullmainmodel += " + ";
            fullmainmodel += sizec(1);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            
            fullmainmodel += " + ";
            fullmainmodel += size(1);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            
            fullmainmodel += " + ";
            fullmainmodel += sizec(1);
            fullmainmodel += ":";
            fullmainmodel += size(2);
            
            if (sizebused) {
              fullmainmodel += " + ";
              fullmainmodel += sizeb(1);
              fullmainmodel += ":";
              fullmainmodel += sizec(2);
              
              fullmainmodel += " + ";
              fullmainmodel += sizec(1);
              fullmainmodel += ":";
              fullmainmodel += sizeb(2);
            }
          }
          
          if (suite == "full" && densityused) {
            fullmainmodel += " + ";
            fullmainmodel += sizec(2);
            fullmainmodel += ":";
            fullmainmodel += densitycol;
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
      }
      
      if (historical) {
        fullmainmodel += " + ";
        fullmainmodel += repst(2);
        
        if (suite == "repst" || suite == "full") {
          fullmainmodel += " + ";
          fullmainmodel += repst(1);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
        }
        
        if (suite == "full" && densityused) {
          fullmainmodel += " + ";
          fullmainmodel += repst(2);
          fullmainmodel += ":";
          fullmainmodel += densitycol;
        }
      }
    }
    
    if (suite == "full") {
      fullmainmodel += " + ";
      fullmainmodel += size(1);
      fullmainmodel += ":";
      fullmainmodel += repst(1);
      
      if (sizebused) {
        fullmainmodel += " + ";
        fullmainmodel += sizeb(1);
        fullmainmodel += ":";
        fullmainmodel += repst(1);
      }
      
      if (sizecused) {
        fullmainmodel += " + ";
        fullmainmodel += sizec(1);
        fullmainmodel += ":";
        fullmainmodel += repst(1);
      }
      
      if (historical) {
        fullmainmodel += " + ";
        fullmainmodel += repst(2);
        fullmainmodel += ":";
        fullmainmodel += size(2);
        
        fullmainmodel += " + ";
        fullmainmodel += size(1);
        fullmainmodel += ":";
        fullmainmodel += size(2);
        
        fullmainmodel += " + ";
        fullmainmodel += repst(1);
        fullmainmodel += ":";
        fullmainmodel += repst(2);
        
        fullmainmodel += " + ";
        fullmainmodel += size(1);
        fullmainmodel += ":";
        fullmainmodel += repst(2);
        
        fullmainmodel += " + ";
        fullmainmodel += repst(1);
        fullmainmodel += ":";
        fullmainmodel += size(2);
        
        if (sizebused) {
          fullmainmodel += " + ";
          fullmainmodel += sizeb(2);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          
          fullmainmodel += " + ";
          fullmainmodel += sizeb(1);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          
          fullmainmodel += " + ";
          fullmainmodel += repst(1);
          fullmainmodel += ":";
          fullmainmodel += sizeb(2);
        }
        
        if (sizecused) {
          fullmainmodel += " + ";
          fullmainmodel += sizec(2);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          
          fullmainmodel += " + ";
          fullmainmodel += sizec(1);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          
          fullmainmodel += " + ";
          fullmainmodel += repst(1);
          fullmainmodel += ":";
          fullmainmodel += sizec(2);
        }
      }
      
      if (age != "none") {
        fullmainmodel += " + ";
        fullmainmodel += age;
        fullmainmodel += ":";
        fullmainmodel += size(1);
        
        if (sizebused) {
          fullmainmodel += " + ";
          fullmainmodel += age;
          fullmainmodel += ":";
          fullmainmodel += sizeb(1);
        }
        if (sizecused) {
          fullmainmodel += " + ";
          fullmainmodel += age;
          fullmainmodel += ":";
          fullmainmodel += sizec(1);
        }
        
        fullmainmodel += " + ";
        fullmainmodel += age;
        fullmainmodel += ":";
        fullmainmodel += repst(1);
        
        if (densityused) {
          fullmainmodel += " + ";
          fullmainmodel += age;
          fullmainmodel += ":";
          fullmainmodel += densitycol;
        }
        
        if (historical) {
          fullmainmodel += " + ";
          fullmainmodel += age;
          fullmainmodel += ":";
          fullmainmodel += size(2);
          
          if (sizebused) {
            fullmainmodel += " + ";
            fullmainmodel += age;
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
          }
          if (sizecused) {
            fullmainmodel += " + ";
            fullmainmodel += age;
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
          }
          
          fullmainmodel += " + ";
          fullmainmodel += age;
          fullmainmodel += ":";
          fullmainmodel += repst(2);
        }
      }
      
      if (indcova(1) != "none" && !iaasrand) {
        fullmainmodel += " + ";
        fullmainmodel += indcova(1);
        fullmainmodel += ":";
        fullmainmodel += size(1);
        
        if (sizebused) {
          fullmainmodel += " + ";
          fullmainmodel += indcova(1);
          fullmainmodel += ":";
          fullmainmodel += sizeb(1);
        }
        if (sizecused) {
          fullmainmodel += " + ";
          fullmainmodel += indcova(1);
          fullmainmodel += ":";
          fullmainmodel += sizec(1);
        }
        
        fullmainmodel += " + ";
        fullmainmodel += indcova(1);
        fullmainmodel += ":";
        fullmainmodel += repst(1);
        
        if (densityused) {
          fullmainmodel += " + ";
          fullmainmodel += indcova(1);
          fullmainmodel += ":";
          fullmainmodel += densitycol;
        }
        
        if (historical && indcova(2) != "none") {
          fullmainmodel += " + ";
          fullmainmodel += indcova(2);
          fullmainmodel += ":";
          fullmainmodel += size(2);
          
          if (sizebused) {
            fullmainmodel += " + ";
            fullmainmodel += indcova(2);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            
            fullmainmodel += " + ";
            fullmainmodel += indcova(1);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            
            fullmainmodel += " + ";
            fullmainmodel += indcova(2);
            fullmainmodel += ":";
            fullmainmodel += sizeb(1);
          }
          if (sizecused) {
            fullmainmodel += " + ";
            fullmainmodel += indcova(2);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            
            fullmainmodel += " + ";
            fullmainmodel += indcova(1);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            
            fullmainmodel += " + ";
            fullmainmodel += indcova(2);
            fullmainmodel += ":";
            fullmainmodel += sizec(1);
          }
        
          fullmainmodel += " + ";
          fullmainmodel += indcova(2);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          
          if (densityused) {
            fullmainmodel += " + ";
            fullmainmodel += indcova(2);
            fullmainmodel += ":";
            fullmainmodel += densitycol;
          }
        }
      }
      if (indcovb(1) != "none" && !ibasrand) {
        fullmainmodel += " + ";
        fullmainmodel += indcovb(1);
        fullmainmodel += ":";
        fullmainmodel += size(1);
        
        if (sizebused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovb(1);
          fullmainmodel += ":";
          fullmainmodel += sizeb(1);
        }
        if (sizecused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovb(1);
          fullmainmodel += ":";
          fullmainmodel += sizec(1);
        }
        
        fullmainmodel += " + ";
        fullmainmodel += indcovb(1);
        fullmainmodel += ":";
        fullmainmodel += repst(1);
        
        if (densityused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovb(1);
          fullmainmodel += ":";
          fullmainmodel += densitycol;
        }
        
        if (historical && indcovb(2) != "none") {
          fullmainmodel += " + ";
          fullmainmodel += indcovb(2);
          fullmainmodel += ":";
          fullmainmodel += size(2);
          
          if (sizebused) {
            fullmainmodel += " + ";
            fullmainmodel += indcovb(2);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            
            fullmainmodel += " + ";
            fullmainmodel += indcovb(1);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            
            fullmainmodel += " + ";
            fullmainmodel += indcovb(2);
            fullmainmodel += ":";
            fullmainmodel += sizeb(1);
          }
          if (sizecused) {
            fullmainmodel += " + ";
            fullmainmodel += indcovb(2);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            
            fullmainmodel += " + ";
            fullmainmodel += indcovb(1);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            
            fullmainmodel += " + ";
            fullmainmodel += indcovb(2);
            fullmainmodel += ":";
            fullmainmodel += sizec(1);
          }
        
          fullmainmodel += " + ";
          fullmainmodel += indcovb(2);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          
          if (densityused) {
            fullmainmodel += " + ";
            fullmainmodel += indcovb(2);
            fullmainmodel += ":";
            fullmainmodel += densitycol;
          }
        }
      }
      
      if (indcovc(1) != "none" && !icasrand) {
        fullmainmodel += " + ";
        fullmainmodel += indcovc(1);
        fullmainmodel += ":";
        fullmainmodel += size(1);
        
        if (sizebused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovc(1);
          fullmainmodel += ":";
          fullmainmodel += sizeb(1);
        }
        if (sizecused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovc(1);
          fullmainmodel += ":";
          fullmainmodel += sizec(1);
        }
        
        fullmainmodel += " + ";
        fullmainmodel += indcovc(1);
        fullmainmodel += ":";
        fullmainmodel += repst(1);
        
        if (densityused) {
          fullmainmodel += " + ";
          fullmainmodel += indcovc(1);
          fullmainmodel += ":";
          fullmainmodel += densitycol;
        }
        
        if (historical && indcovc(2) != "none") {
          fullmainmodel += " + ";
          fullmainmodel += indcovc(2);
          fullmainmodel += ":";
          fullmainmodel += size(2);
          
          if (sizebused) {
            fullmainmodel += " + ";
            fullmainmodel += indcovc(2);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            
            fullmainmodel += " + ";
            fullmainmodel += indcovc(1);
            fullmainmodel += ":";
            fullmainmodel += sizeb(2);
            
            fullmainmodel += " + ";
            fullmainmodel += indcovc(2);
            fullmainmodel += ":";
            fullmainmodel += sizeb(1);
          }
          if (sizecused) {
            fullmainmodel += " + ";
            fullmainmodel += indcovc(2);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            
            fullmainmodel += " + ";
            fullmainmodel += indcovc(1);
            fullmainmodel += ":";
            fullmainmodel += sizec(2);
            
            fullmainmodel += " + ";
            fullmainmodel += indcovc(2);
            fullmainmodel += ":";
            fullmainmodel += sizec(1);
          }
        
          fullmainmodel += " + ";
          fullmainmodel += indcovc(2);
          fullmainmodel += ":";
          fullmainmodel += repst(2);
          
          if (densityused) {
            fullmainmodel += " + ";
            fullmainmodel += indcovc(2);
            fullmainmodel += ":";
            fullmainmodel += densitycol;
          }
        }
      }
    }
    
    //if (modelcounter > 0 && fixedcovcounter > 0) fullmainmodel += " + ";
    
    if (fixedcovcounter > 0) fullmainmodel += fixedtackon;
    if (randomcovcounter > 0) fullmainmodel += randomtackon;
  } else if (suite == "rep") {
    if (age != "none" || covcount > 0) fullmainmodel += " + ";
    fullmainmodel += repst(1);
    
    if (historical && repstcheck) {
      fullmainmodel += " + ";
      fullmainmodel += repst(2);
      
      fullmainmodel += " + ";
      fullmainmodel += repst(1);
      fullmainmodel += ":";
      fullmainmodel += repst(2);
    }
    
    if (fixedcovcounter > 0) fullmainmodel += fixedtackon;
    if (randomcovcounter > 0) fullmainmodel += randomtackon;
  } else if (suite == "cons") {
    if (fixedcovcounter == 0 && modelcounter == 0) fullmainmodel += "1";
    
    //if (modelcounter > 0) fullmainmodel += " + "; // This added
    if (fixedcovcounter > 0) fullmainmodel += fixedtackon;
    
    if (randomcovcounter > 0) fullmainmodel += randomtackon;
    
  } else {
    fullmainmodel = "none";
  }
  
  // Now we will build the global models
  if (survcheck) {
    fullsurvmodel = surv(0);
    fullsurvmodel += fullmainmodel;
    
    if (!nojuvs) {
      juvsurvmodel = surv(0);
      juvsurvmodel += juvmainmodel;
      
      juvmatstmodel = matstat(0);
      juvmatstmodel += juvmainmodel;
    } else {
    juvsurvmodel = "none";
    juvmatstmodel = "none";
    }
    
  } else {
    fullsurvmodel = "none";
  }
  
  if (obscheck) {
    fullobsmodel = obs(0);
    fullobsmodel += fullmainmodel;
    
    if (!nojuvs) {
      juvobsmodel = obs(0);
      juvobsmodel += juvmainmodel;
    } else {
      juvobsmodel = "none";
    }
  } else {
    fullobsmodel = "none";
  }
  
  if (sizecheck) {
    fullsizemodel = size(0);
    fullsizemodel += fullmainmodel;
    
    if (!nojuvs) {
      juvsizemodel = size(0);
      juvsizemodel += juvmainmodel;
    } else {
      juvsizemodel = "none";
    }
    
    if (sizebused) {
      fullsizebmodel = sizeb(0);
      fullsizebmodel += fullmainmodel;
      
      if (!nojuvs) {
        juvsizebmodel = sizeb(0);
        juvsizebmodel += juvmainmodel;
      } else {
        juvsizebmodel = "none";
      }
    }
    if (sizecused) {
      fullsizecmodel = sizec(0);
      fullsizecmodel += fullmainmodel;
      
      if (!nojuvs) {
        juvsizecmodel = sizec(0);
        juvsizecmodel += juvmainmodel;
      } else {
        juvsizecmodel = "none";
      }
    }
  } else {
    fullsizecmodel = "none";
  }
  
  if (repstcheck) {
    fullrepstmodel = repst(0);
    fullrepstmodel += fullmainmodel;
    
    if (!nojuvs) {
      juvrepstmodel = repst(0);
      juvrepstmodel += juvmainmodel;
    } else {
      juvrepstmodel = "none";
    }
  } else {
    fullrepstmodel = "none";
  }
  
  if (feccheck) {
    if (fectime == 3) {
      fullfecmodel = fec(0);
    } else {
      fullfecmodel = fec(1);
    }
    fullfecmodel += fullmainmodel;
  } else {
    fullfecmodel = "none";
  }
  
  StringVector fullnames {"time t", "individual", "patch", "alive in time t+1",
    "observed in time t+1", "sizea in time t+1", "sizeb in time t+1",
    "sizec in time t+1", "reproductive status in time t+1",
    "fecundity in time t+1", "fecundity in time t", "sizea in time t",
    "sizea in time t-1", "sizeb in time t", "sizeb in time t-1", "sizec in time t", 
    "sizec in time t-1", "reproductive status in time t",
    "reproductive status in time t-1", "maturity status in time t+1",
    "maturity status in time t", "age in time t", "density in time t",
    "individual covariate a in time t", "individual covariate a in time t-1",
    "individual covariate b in time t", "individual covariate b in time t-1",
    "individual covariate c in time t", "individual covariate c in time t-1",
    "stage group in time t", "stage group in time t-1"};
  StringVector mainparams {"year2", "individ", "patch", "surv3", "obs3",
    "size3", "sizeb3", "sizec3", "repst3", "fec3", "fec2", "size2", "size1",
    "sizeb2", "sizeb1", "sizec2", "sizec1", "repst2", "repst1", "matst3",
    "matst2", "age", "density", "indcova2", "indcova1", "indcovb2", "indcovb1",
    "indcovc2", "indcovc1", "group2", "group1"};
  
  StringVector modelparams (31);
  modelparams(0) = year;
  modelparams(1) = indiv;
  modelparams(2) = patch;
  modelparams(3) = surv(0);
  modelparams(4) = obs(0);
  modelparams(5) = size(0);
  if (sizebused) {modelparams(6) = sizeb(0);} else {modelparams(6) = "none";}
  if (sizecused) {modelparams(7) = sizec(0);} else {modelparams(7) = "none";}
  modelparams(8) = repst(0);
  if (fectime == 3) {modelparams(9) = fec(0);} else {modelparams(9) = "none";}
  if (fectime == 2) {modelparams(10) = fec(1);} else {modelparams(10) = "none";}
  modelparams(11) = size(1);
  if (historical) {modelparams(12) = size(2);} else {modelparams(12) = "none";}
  if (sizebused) {modelparams(13) = sizeb(1);} else {modelparams(13) = "none";}
  if (sizebused && historical) {modelparams(14) = sizeb(2);} else {modelparams(14) = "none";}
  if (sizecused) {modelparams(15) = sizec(1);} else {modelparams(15) = "none";}
  if (sizecused && historical) {modelparams(16) = sizec(2);} else {modelparams(16) = "none";}
  modelparams(17) = repst(1);
  if (historical) {modelparams(18) = repst(2);} else {modelparams(18) = "none";}
  modelparams(19) = matstat(0);
  modelparams(20) = matstat(1);
  modelparams(21) = age;
  if (densityused) {modelparams(22) = densitycol;} else {modelparams(22) = "none";}
  if (indcovaused) {modelparams(23) = indcova(1);} else {modelparams(23) = "none";}
  if (indcovaused && historical) {modelparams(24) = indcova(2);} else {modelparams(24) = "none";}
  if (indcovbused) {modelparams(25) = indcovb(1);} else {modelparams(25) = "none";}
  if (indcovbused && historical) {modelparams(26) = indcovb(2);} else {modelparams(26) = "none";}
  if (indcovcused) {modelparams(27) = indcovc(1);} else {modelparams(27) = "none";}
  if (indcovcused && historical) {modelparams(28) = indcovc(2);} else {modelparams(28) = "none";}
  if (grouptest) {
    modelparams(29) = "group2";
    if (historical) {modelparams(30) = "group1";} else {modelparams(30) = "group1";}
  } else {
    modelparams(29) = "none";
    modelparams(30) = "none";
  }
  
  Rcpp::DataFrame paramnames = DataFrame::create(Named("parameter_names") = fullnames,
    _["mainparams"] = mainparams, _["modelparams"] = modelparams);
  
  Rcpp::List output = List::create(Named("full.surv.model") = fullsurvmodel,
    _["full.obs.model"] = fullobsmodel, _["full.size.model"] = fullsizemodel,
    _["full.sizeb.model"] = fullsizebmodel, _["full.sizec.model"] = fullsizecmodel,
    _["full.repst.model"] = fullrepstmodel, _["full.fec.model"] = fullfecmodel,
    _["juv.surv.model"] = juvsurvmodel, _["juv.obs.model"] = juvobsmodel,
    _["juv.size.model"] = juvsizemodel, _["juv.sizeb.model"] = juvsizebmodel,
    _["juv.sizec.model"] = juvsizecmodel, _["juv.repst.model"] = juvrepstmodel,
    _["juv.matst.model"] = juvmatstmodel, _["paramnames"] = paramnames);
  
  if (fullsurvmodel == "none") {
    output["full.surv.model"] = 1;
    
    if (!nojuvs) {
      output["juv.surv.model"] = 1;
      output["juv.matst.model"] = 1;
    }
  }
  if (nojuvs) {
    output["juv.surv.model"] = 0;
  }
  
  if (fullobsmodel == "none") {
    output["full.obs.model"] = 1;
    
    if (!nojuvs) {
      output["juv.obs.model"] = 1;
    }
  }
  if (nojuvs) {output["juv.obs.model"] = 0;}
  
  if (fullsizemodel == "none") {
    output["full.size.model"] = 1;
    
    if (!nojuvs) {
      output["juv.size.model"] = 1;
    }
  }
  if (nojuvs) {output["juv.size.model"] = 0;}
  
  if (fullsizebmodel == "none") {
    output["full.sizeb.model"] = 1;
    
    if (!nojuvs) {
      output["juv.sizeb.model"] = 1;
    }
  }
  if (nojuvs) {output["juv.sizeb.model"] = 0;}
  
  if (fullsizecmodel == "none") {
    output["full.sizec.model"] = 1;
    
    if (!nojuvs) {
      output["juv.sizec.model"] = 1;
    }
  }
  if (nojuvs) {output["juv.sizec.model"] = 0;}
  
  if (fullrepstmodel == "none") {
    output["full.repst.model"] = 1;
    if (!nojuvs) {output["juv.repst.model"] = 1;}
  }
  if (nojuvs) {output["juv.repst.model"] = 0;}
  
  if (fullfecmodel == "none") output["full.fec.model"] = 1;
  
  return output;
}

