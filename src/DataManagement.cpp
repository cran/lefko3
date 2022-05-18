#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Create Vertical Structure for Horizontal Data Frame Input
//' 
//' Function \code{pfj()} powers the R function \code{\link{verticalize3}()},
//' creating the vertical structure and rearranging the data in that shape.
//' 
//' @param data The horizontal data file.
//' @param stageframe The stageframe object describing the life history model.
//' This should be the full stageframe.
//' @param noyears The number of years or observation periods in the dataset.
//' @param firstyear The first year or time of observation.
//' @param popidcol Column number corresponding to the identity of the
//' population for each individual.
//' @param patchidcol Column number corresponding to the identity of the patch
//' for each individual.
//' @param individcol Column number corresponding to the identity of each 
//' individual.
//' @param blocksize The number of variables corresponding to each time step in 
//' the input dataset designated in \code{data}.
//' @param xcol Vector of column numbers corresponding to the x coordinate of
//' each individual in Cartesian space.
//' @param ycol Vector of column numbers corresponding to the y coordinate of
//' each individual in Cartesian space.
//' @param juvcol Vector of column numbers that marks individuals in immature
//' stages within the dataset.
//' @param sizeacol Vector of column numbers corresponding to the first or main
//' size variable associated with the first year or observation time in the
//' dataset.
//' @param sizebcol Vector of column numbers corresponding to the second size
//' variable associated with the first year or observation time in the dataset.
//' @param sizeccol Vector of column numbers corresponding to the third size
//' variable associated with the first year or observation time in the dataset.
//' @param repstracol Vector of column numbers corresponding to the main 
//' variable coding the production of reproductive structures associated with
//' the first year or observation period in the input dataset.
//' @param repstrbcol Vector of column numbers corresponding to a second
//' variable coding the production of reproductive structures associated with
//' the first year or observation period in the input dataset.
//' @param fecacol Vector of column numbers corresponding to the main variable
//' coding for fecundity associated with the first year or observation period in
//' the dataset.
//' @param fecbcol Vector of column numbers corresponding to a second variable
//' coding for fecundity associated with the first year or observation period in
//' the dataset.
//' @param indcovacol Vector of column numbers corresponding to an individual
//' covariate.
//' @param indcovbcol Vector of column numbers corresponding to an individual
//' covariate.
//' @param indcovccol Vector of column numbers corresponding to an individual
//' covariate.
//' @param aliveacol Vector of column numbers that details whether an individual
//' is alive at a given time.
//' @param deadacol Vector of column numbers that details whether an individual
//' is dead at a given time.
//' @param obsacol Vector of column numbers that details whether an individual
//' is in an observable stage at a given time.
//' @param nonobsacol Vector of column numbers that details whether an
//' individual is in an unobservable stage at a given time.
//' @param censorcol Vector of column numbers corresponding to the first entry
//' of a censor variable.
//' @param stagecol Vector of column numbers corresponding to the first entry of
//' a column designating stages.
//' @param repstrrel This is a scalar modifier for that makes the variable in
//' \code{repstrbcol} equivalent to \code{repstracol}.
//' @param fecrel This is a scalar modifier for that makes the variable in
//' \code{fecbcol} equivalent to \code{fecacol}.
//' @param NAas0 If TRUE, then all NA entries for size and fecundity variables
//' will be set to 0.
//' @param NRasRep If TRUE, then will treat non-reproductive but mature
//' individuals as reproductive during stage assignment.
//' @param RepasObs If TRUE, then will treat individuals with size 0 as observed
//' if and only if they are reproductive. Otherwise, all individuals with size 0
//' are treated as not observed.
//' @param NOasObs If TRUE, then will treat unobserved individuals as observed
//' during stage assignment.
//' @param stassign A logical value indicating whether to assign stages.
//' @param stszcol Integer describing which size variable or combination of size
//' variables to use in stage estimation.
//' @param censorkeep The value of the censoring variable identifying data
//' that should be included in analysis. Defaults to 0, but may take any numeric
//' value including NA.
//' @param censbool A logical variable determining whether NA denotes the value
//' of the censoring variable identifying data to keep. If used, then will set
//' all NAs to 0 and all other values to 1, treating 0 as the value to keep.
//' @param censrepeat A logical value indicating whether censor variable is a
//' single static column, or whether censor variables repeat across blocks.
//' @param coordsrepeat A logical value indicating whether coordinate variables
//' are single static columns, or whether they repeat across blocks.
//' @param quiet A logical value indicating whether to silense warnings.
//' 
//' @return The output is currently a 7 element list, where each element is a
//' data frame with the same number of rows.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.pfj)]]
Rcpp::List pfj(DataFrame data, DataFrame stageframe, int noyears, int firstyear,
  int popidcol, int patchidcol, int individcol, int blocksize, arma::ivec xcol,
  arma::ivec ycol, arma::ivec juvcol, arma::ivec sizeacol, arma::ivec sizebcol,
  arma::ivec sizeccol, arma::ivec repstracol, arma::ivec repstrbcol,
  arma::ivec fecacol, arma::ivec fecbcol, arma::ivec indcovacol,
  arma::ivec indcovbcol, arma::ivec indcovccol, arma::ivec aliveacol,
  arma::ivec deadacol, arma::ivec obsacol, arma::ivec nonobsacol,
  arma::ivec censorcol, arma::ivec stagecol, double repstrrel, double fecrel,
  bool NAas0, bool NRasRep, bool RepasObs, bool NOasObs, bool stassign,
  int stszcol, double censorkeep, bool censbool, bool censrepeat,
  bool coordsrepeat, bool quiet) {
  
  Rcpp::List output_longlist(81);
  int noindivs = data.nrows();
  int novars = data.length();
  
  // Vector length checks and standardizing vectors of data frame columns
  if (xcol.n_elem == 1) {
    xcol.resize(noyears);
    
    if (!coordsrepeat || blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        xcol(i) = xcol(0);
      }
    } else if (xcol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        xcol(i) = xcol(i-1) + blocksize;
        if (xcol(i) >= novars) {
          throw Rcpp::exception("Vector xcol is too small. Either noyears is too small or not enough data blocks have been input.");
        }
      }
    } else {
      xcol.fill(-1);
    }
  }
  
  if (ycol.n_elem == 1) {
    ycol.resize(noyears);
    
    if (!coordsrepeat || blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        ycol(i) = ycol(0);
      }
    } else if (ycol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        ycol(i) = ycol(i-1) + blocksize;
        if (ycol(i) >= novars) {
          throw Rcpp::exception("Vector ycol is too small. Either noyears is too small or not enough data blocks have been input.");
        }
      }
    } else {
      ycol.fill(-1);
    }
  }
  
  if (juvcol.n_elem == 1) {
    juvcol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        juvcol(i) = juvcol(0);
      }
    } else if (juvcol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        juvcol(i) = juvcol(i-1) + blocksize;
        if (juvcol(i) >= novars) {
          throw Rcpp::exception("Vector juvcol is too small. Either noyears is too small or not enough data blocks have been input.");
        }
      }
    } else {
      juvcol.fill(-1);
    }
  }
  
  if (sizeacol.n_elem == 1) {
    sizeacol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        sizeacol(i) = sizeacol(0);
      }
    } else if (sizeacol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        sizeacol(i) = sizeacol(i-1) + blocksize;
        if (sizeacol(i) >= novars) {
          throw Rcpp::exception("Vector sizeacol is too small. Either noyears is too small or not enough data blocks have been input.");
        }
      }
    } else {
      sizeacol.fill(-1);
    }
  }
  
  if (sizebcol.n_elem == 1) {
    sizebcol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        sizebcol(i) = sizebcol(0);
      }
    } else if (sizebcol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        sizebcol(i) = sizebcol(i-1) + blocksize;
        if (sizebcol(i) >= novars) {
          throw Rcpp::exception("Vector sizebcol is too small. Either noyears is too small or not enough data blocks have been input.");
        }
      }
    } else {
      sizebcol.fill(-1);
    }
  }
  
  if (sizeccol.n_elem == 1) {
    sizeccol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        sizeccol(i) = sizeccol(0);
      }
    } else if (sizeccol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        sizeccol(i) = sizeccol(i-1) + blocksize;
        if (sizeccol(i) >= novars) {
          throw Rcpp::exception("Vector sizeccol is too small. Either noyears is too small or not enough data blocks have been input.");
        }
      }
    } else {
      sizeccol.fill(-1);
    }
  }
  
  if (repstracol.n_elem == 1) {
    repstracol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        repstracol(i) = repstracol(0);
      }
    } else if (repstracol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        repstracol(i) = repstracol(i-1) + blocksize;
        if (repstracol(i) >= novars) {
          throw Rcpp::exception("Vector repstracol is too small. Either noyears is too small or not enough data blocks have been input.");
        }
      }
    } else {
      repstracol.fill(-1);
    }
  }
  
  if (repstrbcol.n_elem == 1) {
    repstrbcol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        repstrbcol(i) = repstrbcol(0);
      }
    } else if (repstrbcol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        repstrbcol(i) = repstrbcol(i-1) + blocksize;
        if (repstrbcol(i) >= novars) {
          throw Rcpp::exception("Vector repstrbcol is too small. Either noyears is too small or not enough data blocks have been input.");
        }
      }
    } else {
      repstrbcol.fill(-1);
    }
  }
  
  if (fecacol.n_elem == 1) {
    fecacol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        fecacol(i) = fecacol(0);
      }
    } else if (fecacol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        fecacol(i) = fecacol(i-1) + blocksize;
        if (fecacol(i) >= novars) {
          throw Rcpp::exception("Vector fecacol is too small. Either noyears is too small or not enough data blocks have been input.");
        }
      }
    } else {
      fecacol.fill(-1);
    }
  }
  
  if (fecbcol.n_elem == 1) {
    fecbcol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        fecbcol(i) = fecbcol(0);
      }
    } else if (fecbcol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        fecbcol(i) = fecbcol(i-1) + blocksize;
        if (fecbcol(i) >= novars) {
          throw Rcpp::exception("Vector fecbcol is too small. Either noyears is too small or not enough data blocks have been input.");
        }
      }
    } else {
      fecbcol.fill(-1);
    }
  }
  
  if (indcovacol.n_elem == 1) {
    indcovacol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        indcovacol(i) = indcovacol(0);
      }
    } else if (indcovacol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        indcovacol(i) = indcovacol(i-1) + blocksize;
        if (indcovacol(i) >= novars) {
          throw Rcpp::exception("Vector indcovacol is too small. Either noyears is too small or not enough data blocks have been input.");
        }
      }
    } else {
      indcovacol.fill(-1);
    }
  }
  
  if (indcovbcol.n_elem == 1) {
    indcovbcol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        indcovbcol(i) = indcovbcol(0);
      }
    } else if (indcovbcol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        indcovbcol(i) = indcovbcol(i-1) + blocksize;
        if (indcovbcol(i) >= novars) {
          throw Rcpp::exception("EVector indcovbcol is too small. ither noyears is too small or not enough data blocks have been input.");
        }
      }
    } else {
      indcovbcol.fill(-1);
    }
  }
  
  if (indcovccol.n_elem == 1) {
    indcovccol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        indcovccol(i) = indcovccol(0);
      }
    } else if (indcovccol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        indcovccol(i) = indcovccol(i-1) + blocksize;
        if (indcovccol(i) >= novars) {
          throw Rcpp::exception("Vector indcovccol is too small. Either noyears is too small or not enough data blocks have been input.");
        }
      }
    } else {
      indcovccol.fill(-1);
    }
  }
  
  if (aliveacol.n_elem == 1) {
    aliveacol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        aliveacol(i) = aliveacol(0);
      }
    } else if (aliveacol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        aliveacol(i) = aliveacol(i-1) + blocksize;
        if (aliveacol(i) >= novars) {
          throw Rcpp::exception("Vector aliveacol is too small. Either noyears is too small or not enough data blocks have been input.");
        }
      }
    } else {
      aliveacol.fill(-1);
    }
  }
  
  if (deadacol.n_elem == 1) {
    deadacol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        deadacol(i) = deadacol(0);
      }
    } else if (deadacol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        deadacol(i) = deadacol(i-1) + blocksize;
        if (deadacol(i) >= novars) {
          throw Rcpp::exception("Vector deadacol is too small. Either noyears is too small or not enough data blocks have been input.");
        }
      }
    } else {
      deadacol.fill(-1);
    }
  }
  
  if (obsacol.n_elem == 1) {
    obsacol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        obsacol(i) = obsacol(0);
      }
    } else if (obsacol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        obsacol(i) = obsacol(i-1) + blocksize;
        if (obsacol(i) >= novars) {
          throw Rcpp::exception("Vector obsacol is too small. Either noyears is too small or not enough data blocks have been input.");
        }
      }
    } else {
      obsacol.fill(-1);
    }
  }
  
  if (nonobsacol.n_elem == 1) {
    nonobsacol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        nonobsacol(i) = nonobsacol(0);
      }
    } else if (nonobsacol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        nonobsacol(i) = nonobsacol(i-1) + blocksize;
        if (nonobsacol(i) >= novars) {
          throw Rcpp::exception("Vector nonobsacol is too small. Either noyears is too small or not enough data blocks have been input.");
        }
      }
    } else {
      nonobsacol.fill(-1);
    }
  }
  
  if (censorcol.n_elem == 1) {
    censorcol.resize(noyears);
    
    if (blocksize == 0 || !censrepeat) {
      for (int i = 1; i < noyears; i++) {
        censorcol(i) = censorcol(0);
      }
    } else if (blocksize > 0 && censorcol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        censorcol(i) = censorcol(i-1) + blocksize;
        if (censorcol(i) >= novars) {
          throw Rcpp::exception("Vector censorcol is too small. Either noyears is too small or not enough data blocks have been input.");
        }
      }
    } else {
      censorcol.fill(-1);
    }
  }
  
  if (stagecol.n_elem == 1) {
    stagecol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        stagecol(i) = stagecol(0);
      }
    } else if (stagecol(0) != -1) {
      for (int i = 1; i < noyears; i++) {
        stagecol(i) = stagecol(i-1) + blocksize;
        if (stagecol(i) >= novars) {
          throw Rcpp::exception("Vector stagecol is too small. Either noyears is too small or not enough data blocks have been input.");
        }
      }
    } else {
      stagecol.fill(-1);
    }
  }
  
  // Variable initialization for most vectors
  double livcheck1 {0.0};
  double livcheck2 {0.0};
  double livcheck3 {0.0};
  
  double nonobssub1 {0.0};
  double nonobssub2 {0.0};
  double nonobssub3 {0.0};
  
  double stagesize1 {0.0};
  double stagesize2 {0.0};
  double stagesize3 {0.0};
  double stagesize1b {0.0};
  double stagesize2b {0.0};
  double stagesize3b {0.0};
  double stagesize1c {0.0};
  double stagesize2c {0.0};
  double stagesize3c {0.0};
  
  arma::uvec stagemini1;
  arma::uvec stagemaxi1;
  arma::uvec stagemini2;
  arma::uvec stagemaxi2;
  arma::uvec stagemini3;
  arma::uvec stagemaxi3;
  arma::uvec stageobs;
  arma::uvec stagerep;
  arma::uvec stagemat;
  arma::uvec cs1;
  arma::uvec cs2;
  arma::uvec cs3;
  arma::uvec cs4;
  int choicestage;
  
  Rcpp::StringVector sfname = as<StringVector>(stageframe["stage"]);
  Rcpp::NumericVector repstat = as<NumericVector>(stageframe["repstatus"]);
  Rcpp::NumericVector obsstat = as<NumericVector>(stageframe["obsstatus"]);
  Rcpp::NumericVector matstat = as<NumericVector>(stageframe["matstatus"]);
  arma::vec indataset = as<arma::vec>(stageframe["indataset"]);
  arma::vec sfszmin = as<arma::vec>(stageframe["sizebin_min"]);
  arma::vec sfszmax = as<arma::vec>(stageframe["sizebin_max"]);
  arma::vec sfszminb = as<arma::vec>(stageframe["sizebinb_min"]);
  arma::vec sfszmaxb = as<arma::vec>(stageframe["sizebinb_max"]);
  arma::vec sfszminc = as<arma::vec>(stageframe["sizebinc_min"]);
  arma::vec sfszmaxc = as<arma::vec>(stageframe["sizebinc_max"]);
  
  arma::vec repstatarma = as<arma::vec>(repstat);
  arma::vec obsstatarma = as<arma::vec>(obsstat);
  arma::vec matstatarma = as<arma::vec>(matstat);
  arma::vec sfszminarma = sfszmin;
  arma::vec sfszmaxarma = sfszmax;
  arma::vec sfszminarmab = sfszminb;
  arma::vec sfszmaxarmab = sfszmaxb;
  arma::vec sfszminarmac = sfszminc;
  arma::vec sfszmaxarmac = sfszmaxc;
  int stagenum = sfszmaxarma.n_elem;
  
  arma::uvec instages = find(indataset == 1);
  int instagenum = instages.n_elem;
  
  arma::uvec stageid (stagenum);
  arma::uvec instageid (instagenum);
  arma::vec insfszminarma (instagenum);
  arma::vec insfszmaxarma (instagenum);
  arma::vec inrepstatarma (instagenum);
  arma::vec inobsstatarma (instagenum);
  arma::vec inmatstatarma (instagenum);
  arma::vec insfszminarmab (instagenum);
  arma::vec insfszmaxarmab (instagenum);
  arma::vec insfszminarmac (instagenum);
  arma::vec insfszmaxarmac (instagenum);

  int inplace {0};
  for (int i = 0; i < stagenum; i++) {
    stageid(i) = i+1;
    
    if (indataset(i) == 1) {
      instageid(inplace) = i + 1;
      insfszminarma(inplace) = sfszminarma(i);
      insfszmaxarma(inplace) = sfszmaxarma(i);
      insfszminarmab(inplace) = sfszminarmab(i);
      insfszmaxarmab(inplace) = sfszmaxarmab(i);
      insfszminarmac(inplace) = sfszminarmac(i);
      insfszmaxarmac(inplace) = sfszmaxarmac(i);
      inrepstatarma(inplace) = repstatarma(i);
      inobsstatarma(inplace) = obsstatarma(i);
      inmatstatarma(inplace) = matstatarma(i);
      
      inplace++;
    }
  }
  
  Rcpp::StringVector popidx (noindivs);
  Rcpp::StringVector patchidx (noindivs);
  Rcpp::StringVector individx (noindivs);
  Rcpp::StringVector stage1x (noindivs);
  Rcpp::StringVector stage2x (noindivs);
  Rcpp::StringVector stage3x (noindivs);
  
  Rcpp::NumericVector xpos1x (noindivs);
  Rcpp::NumericVector ypos1x (noindivs);
  Rcpp::NumericVector xpos2x (noindivs);
  Rcpp::NumericVector ypos2x (noindivs);
  Rcpp::NumericVector xpos3x (noindivs);
  Rcpp::NumericVector ypos3x (noindivs);
  
  Rcpp::NumericVector sizea1x (noindivs);
  Rcpp::NumericVector sizea2x (noindivs);
  Rcpp::NumericVector sizea3x (noindivs);
  Rcpp::NumericVector repstra1x (noindivs);
  Rcpp::NumericVector repstra2x (noindivs);
  Rcpp::NumericVector repstra3x (noindivs);
  Rcpp::NumericVector feca1x (noindivs);
  Rcpp::NumericVector feca2x (noindivs);
  Rcpp::NumericVector feca3x (noindivs);
  Rcpp::NumericVector sizeb1x (noindivs);
  Rcpp::NumericVector sizeb2x (noindivs);
  Rcpp::NumericVector sizeb3x (noindivs);
  Rcpp::NumericVector repstrb1x (noindivs);
  Rcpp::NumericVector repstrb2x (noindivs);
  Rcpp::NumericVector repstrb3x (noindivs);
  Rcpp::NumericVector fecb1x (noindivs);
  Rcpp::NumericVector fecb2x (noindivs);
  Rcpp::NumericVector fecb3x (noindivs);
  Rcpp::NumericVector sizec1x (noindivs);
  Rcpp::NumericVector sizec2x (noindivs);
  Rcpp::NumericVector sizec3x (noindivs);
  
  Rcpp::NumericVector indcova1x (noindivs);
  Rcpp::NumericVector indcova2x (noindivs);
  Rcpp::NumericVector indcova3x (noindivs);
  Rcpp::NumericVector indcovb1x (noindivs);
  Rcpp::NumericVector indcovb2x (noindivs);
  Rcpp::NumericVector indcovb3x (noindivs);
  Rcpp::NumericVector indcovc1x (noindivs);
  Rcpp::NumericVector indcovc2x (noindivs);
  Rcpp::NumericVector indcovc3x (noindivs);
  
  Rcpp::NumericVector censor1x (noindivs);
  Rcpp::NumericVector censor2x (noindivs);
  Rcpp::NumericVector censor3x (noindivs);
  Rcpp::NumericVector alivegiven1x (noindivs);
  Rcpp::NumericVector alivegiven2x (noindivs);
  Rcpp::NumericVector alivegiven3x (noindivs);
  Rcpp::NumericVector deadgiven1x (noindivs);
  Rcpp::NumericVector deadgiven2x (noindivs);
  Rcpp::NumericVector deadgiven3x (noindivs);
  Rcpp::NumericVector obsgiven1x (noindivs);
  Rcpp::NumericVector obsgiven2x (noindivs);
  Rcpp::NumericVector obsgiven3x (noindivs);
  Rcpp::NumericVector nonobsgiven1x (noindivs);
  Rcpp::NumericVector nonobsgiven2x (noindivs); 
  Rcpp::NumericVector nonobsgiven3x (noindivs);
  Rcpp::NumericVector juvgiven1x (noindivs);
  Rcpp::NumericVector juvgiven2x (noindivs);
  Rcpp::NumericVector juvgiven3x (noindivs);
  nonobsgiven1x.fill(-1.0);
  nonobsgiven2x.fill(-1.0);
  nonobsgiven3x.fill(-1.0);
  
  Rcpp::NumericVector zerovec (noindivs, 0.0);
  Rcpp::NumericVector onevec (noindivs, 1.0);
  Rcpp::NumericVector negonevec (noindivs, -1.0);
  // zerovec.fill(0.0);
  // onevec.fill(1.0);
  // negonevec.fill(-1.0);
  
  int ndflength = noindivs * (noyears - 1);
  
  Rcpp::NumericVector rowid (ndflength);
  Rcpp::StringVector popid (ndflength);
  Rcpp::StringVector patchid (ndflength);
  Rcpp::StringVector individ (ndflength);
  Rcpp::NumericVector year2 (ndflength);
  Rcpp::NumericVector xpos1 (ndflength);
  Rcpp::NumericVector ypos1 (ndflength);
  Rcpp::NumericVector xpos2 (ndflength); 
  Rcpp::NumericVector ypos2 (ndflength);
  Rcpp::NumericVector xpos3 (ndflength);
  Rcpp::NumericVector ypos3 (ndflength);
  
  Rcpp::NumericVector sizea1 (ndflength);
  Rcpp::NumericVector sizea2 (ndflength);
  Rcpp::NumericVector sizea3 (ndflength);
  Rcpp::NumericVector sizea10 (ndflength);
  Rcpp::NumericVector sizea20 (ndflength);
  Rcpp::NumericVector sizea30 (ndflength);
  Rcpp::NumericVector sizeb1 (ndflength);
  Rcpp::NumericVector sizeb2 (ndflength); 
  Rcpp::NumericVector sizeb3 (ndflength);
  Rcpp::NumericVector sizeb10 (ndflength);
  Rcpp::NumericVector sizeb20 (ndflength); 
  Rcpp::NumericVector sizeb30 (ndflength);
  Rcpp::NumericVector sizec1 (ndflength);
  Rcpp::NumericVector sizec2 (ndflength);
  Rcpp::NumericVector sizec3 (ndflength);
  Rcpp::NumericVector sizec10 (ndflength);
  Rcpp::NumericVector sizec20 (ndflength);
  Rcpp::NumericVector sizec30 (ndflength);
  
  Rcpp::NumericVector repstra1 (ndflength);
  Rcpp::NumericVector repstra2 (ndflength);
  Rcpp::NumericVector repstra3 (ndflength);
  Rcpp::NumericVector repstra10 (ndflength);
  Rcpp::NumericVector repstra20 (ndflength);
  Rcpp::NumericVector repstra30 (ndflength);
  Rcpp::NumericVector repstrb1 (ndflength);
  Rcpp::NumericVector repstrb2 (ndflength);
  Rcpp::NumericVector repstrb3 (ndflength);
  Rcpp::NumericVector repstrb10 (ndflength);
  Rcpp::NumericVector repstrb20 (ndflength);
  Rcpp::NumericVector repstrb30 (ndflength);
  
  Rcpp::NumericVector feca1 (ndflength);
  Rcpp::NumericVector feca2 (ndflength);
  Rcpp::NumericVector feca3 (ndflength);
  Rcpp::NumericVector feca10 (ndflength);
  Rcpp::NumericVector feca20 (ndflength);
  Rcpp::NumericVector feca30 (ndflength);
  Rcpp::NumericVector fecb1 (ndflength);
  Rcpp::NumericVector fecb2 (ndflength);
  Rcpp::NumericVector fecb3 (ndflength);
  Rcpp::NumericVector fecb10 (ndflength);
  Rcpp::NumericVector fecb20 (ndflength);
  Rcpp::NumericVector fecb30 (ndflength);
  
  Rcpp::NumericVector indcova1 (ndflength);
  Rcpp::NumericVector indcova2 (ndflength);
  Rcpp::NumericVector indcova3 (ndflength);
  Rcpp::NumericVector indcovb1 (ndflength);
  Rcpp::NumericVector indcovb2 (ndflength);
  Rcpp::NumericVector indcovb3 (ndflength);
  Rcpp::NumericVector indcovc1 (ndflength);
  Rcpp::NumericVector indcovc2 (ndflength);
  Rcpp::NumericVector indcovc3 (ndflength);
  
  Rcpp::NumericVector censor1 (ndflength);
  Rcpp::NumericVector censor2 (ndflength);
  Rcpp::NumericVector censor3 (ndflength);
  Rcpp::NumericVector alivegiven1 (ndflength);
  Rcpp::NumericVector alivegiven2 (ndflength);
  Rcpp::NumericVector alivegiven3 (ndflength);
  Rcpp::NumericVector deadgiven1 (ndflength);
  Rcpp::NumericVector deadgiven2 (ndflength);
  Rcpp::NumericVector deadgiven3 (ndflength);
  Rcpp::NumericVector obsgiven1 (ndflength);
  Rcpp::NumericVector obsgiven2 (ndflength);
  Rcpp::NumericVector obsgiven3 (ndflength);
  Rcpp::NumericVector nonobsgiven1 (ndflength);
  Rcpp::NumericVector nonobsgiven2 (ndflength); 
  Rcpp::NumericVector nonobsgiven3 (ndflength);
  Rcpp::NumericVector juvgiven1 (ndflength);
  Rcpp::NumericVector juvgiven2 (ndflength);
  Rcpp::NumericVector juvgiven3 (ndflength);
  
  Rcpp::NumericVector addedsize1 (ndflength);
  Rcpp::NumericVector addedsize2 (ndflength);
  Rcpp::NumericVector addedsize3 (ndflength);
  Rcpp::NumericVector addedrepstr1 (ndflength);
  Rcpp::NumericVector addedrepstr2 (ndflength);
  Rcpp::NumericVector addedrepstr3 (ndflength);
  Rcpp::NumericVector addedfec1 (ndflength);
  Rcpp::NumericVector addedfec2 (ndflength);
  Rcpp::NumericVector addedfec3 (ndflength);
  
  Rcpp::NumericVector spryn1 (ndflength);
  Rcpp::NumericVector spryn2 (ndflength);
  Rcpp::NumericVector spryn3 (ndflength);
  Rcpp::NumericVector repyn1 (ndflength);
  Rcpp::NumericVector repyn2 (ndflength);
  Rcpp::NumericVector repyn3 (ndflength);
  Rcpp::NumericVector fecyn1 (ndflength);
  Rcpp::NumericVector fecyn2 (ndflength);
  Rcpp::NumericVector fecyn3 (ndflength);
  Rcpp::NumericVector matstat1 (ndflength);
  Rcpp::NumericVector matstat2 (ndflength);
  Rcpp::NumericVector matstat3 (ndflength);
  
  Rcpp::NumericVector alive1 (ndflength);
  Rcpp::NumericVector alive2 (ndflength);
  Rcpp::NumericVector alive3 (ndflength);
  
  Rcpp::StringVector stage1 (ndflength);
  Rcpp::StringVector stage2 (ndflength);
  Rcpp::StringVector stage3 (ndflength);
  Rcpp::NumericVector stage1num (ndflength);
  Rcpp::NumericVector stage2num (ndflength);
  Rcpp::NumericVector stage3num (ndflength);
  
  Rcpp::NumericVector firstseenx (noindivs, 0.0);
  Rcpp::NumericVector lastseenx (noindivs, 0.0);
  Rcpp::NumericVector firstseen (ndflength);
  Rcpp::NumericVector lastseen (ndflength);
  Rcpp::NumericVector obsage (ndflength);
  Rcpp::NumericVector obslifespan (ndflength);
  // firstseenx.fill(0.0);
  // lastseenx.fill(0.0);
  
  for (int i = 0; i < noindivs; i++) {
    // firstseenx[i] = 0.0;
    // lastseenx[i] = 0.0;
    
    Rcpp::NumericVector temp1 (noindivs, 0.0);
    Rcpp::NumericVector temp2 (noindivs, 0.0);
    Rcpp::NumericVector temp3 (noindivs, 0.0);
    Rcpp::NumericVector temp4 (noindivs, 0.0);
    Rcpp::NumericVector temp5 (noindivs, 0.0);
    
    for (int k = 0; k < noyears; k++) {
      if (sizeacol(k) != -1 ) {
        temp1 = clone(as<NumericVector>(data[sizeacol(k)]));
      } else {temp1 = clone(zerovec);}
      if (sizebcol(k) != -1 ) {
        temp2 = clone(as<NumericVector>(data[sizebcol(k)]));
      } else {temp2 = clone(zerovec);}
      if (sizeccol(k) != -1 ) {
        temp3 = clone(as<NumericVector>(data[sizeccol(k)]));
      } else {temp3 = clone(zerovec);}
      if (repstracol(k) != -1 ) {
        temp4 = clone(as<NumericVector>(data[repstracol(k)]));
      } else {temp4 = clone(zerovec);}
      if (repstrbcol(k) != -1 ) {
        temp5 = clone(as<NumericVector>(data[repstrbcol(k)]));
      } else {temp5 = clone(zerovec);}
      
      if (NumericVector::is_na(temp1(i))) {temp1(i) = 0.0;}
      if (NumericVector::is_na(temp2(i))) {temp2(i) = 0.0;}
      if (NumericVector::is_na(temp3(i))) {temp3(i) = 0.0;}
      if (NumericVector::is_na(temp4(i))) {temp4(i) = 0.0;}
      if (NumericVector::is_na(temp5(i))) {temp5(i) = 0.0;}
      
      if (temp1(i) > 0.0 && firstseenx(i) == 0.0) {
        firstseenx(i) = firstyear + static_cast<double>(k);
      } else if (temp2(i) > 0.0 && firstseenx(i) == 0.0) {
        firstseenx(i) = firstyear + static_cast<double>(k);
      } else if (temp3(i) > 0.0 && firstseenx(i) == 0.0) {
        firstseenx(i) = firstyear + static_cast<double>(k);
      } else if (temp4(i) > 0.0 && firstseenx(i) == 0.0) {
        firstseenx(i) = firstyear + static_cast<double>(k);
      } else if (temp5(i) > 0.0 && firstseenx(i) == 0.0) {
        firstseenx(i) = firstyear + static_cast<double>(k);
      }
      
      if (temp1(i) > 0 && firstseenx(i) != 0.0) {
        lastseenx(i) = firstyear + static_cast<double>(k);
      } else if (temp2(i) > 0 && firstseenx(i) != 0.0) {
        lastseenx(i) = firstyear + static_cast<double>(k);
      } else if (temp3(i) > 0 && firstseenx(i) != 0.0) {
       lastseenx(i) = firstyear + static_cast<double>(k);
      } else if (temp4(i) > 0 && firstseenx(i) != 0.0) {
        lastseenx(i) = firstyear + static_cast<double>(k);
      } else if (temp5(i) > 0 && firstseenx(i) != 0.0) {
        lastseenx(i) = firstyear + static_cast<double>(k);
      }
    }
  }
  
  // Set up main loop
  for (int j = 0; j < (noyears - 1); j++) {
    Rcpp::checkUserInterrupt();
    
    if (popidcol > -1) popidx = as<StringVector>(data[popidcol]);
    if (patchidcol > -1) patchidx = as<StringVector>(data[patchidcol]);
    if (individcol > -1) individx = as<StringVector>(data[individcol]);
    
    if (j == 0) { // This is implemented in the first year
      
      sizea1x = clone(zerovec);
      sizea2x = data[sizeacol(0)];
      sizea3x = data[sizeacol(1)];
      
      if (stagecol(0) > -1) {
        stage1x.fill("NotAlive");
        stage2x = data[stagecol(0)];
        stage3x = data[stagecol(1)];
      }
      
      if (xcol(0) > -1) {
        xpos1x = clone(zerovec);
        xpos2x = data[xcol(0)];
        xpos3x = data[xcol(1)];
      }
      
      if (ycol(0) > -1) {
        ypos1x = clone(zerovec);
        ypos2x = data[ycol(0)];
        ypos3x = data[ycol(1)];
      }
      
      if (repstracol(0) > -1) {
        repstra1x = clone(zerovec);
        repstra2x = data[repstracol(0)];
        repstra3x = data[repstracol(1)];
      }
      
      if (fecacol(0) > -1) {
        feca1x = clone(zerovec);
        feca2x = data[fecacol(0)];
        feca3x = data[fecacol(1)];
      }
      
      if (sizebcol(0) > -1) {
        sizeb1x = clone(zerovec);
        sizeb2x = data[sizebcol(0)];
        sizeb3x = data[sizebcol(1)];
      }
      
      if (repstrbcol(0) > -1) {
        repstrb1x = clone(zerovec);
        repstrb2x = data[repstrbcol(0)];
        repstrb3x = data[repstrbcol(1)];
      }
      
      if (fecbcol(0) > -1) {
        fecb1x = clone(zerovec);
        fecb2x = data[fecbcol(0)];
        fecb3x = data[fecbcol(1)];
      }
      
      if (sizeccol(0) > -1) {
        sizec1x = clone(zerovec);
        sizec2x = data[sizeccol(0)];
        sizec3x = data[sizeccol(1)];
      }
      
      if (indcovacol(0) > -1) {
        indcova1x = clone(zerovec);
        indcova2x = data[indcovacol(0)];
        indcova3x = data[indcovacol(1)];
      }
      
      if (indcovbcol(0) > -1) {
        indcovb1x = clone(zerovec);
        indcovb2x = data[indcovbcol(0)];
        indcovb3x = data[indcovbcol(1)];
      }
      
      if (indcovccol(0) > -1) {
        indcovc1x = clone(zerovec);
        indcovc2x = data[indcovccol(0)];
        indcovc3x = data[indcovccol(1)];
      }
      
      // This section decides how to deal with censor variables
      if (censorcol(0) > -1) {
        
        if (censrepeat) {
        
          if (censorkeep == 0) {
            censor1x = clone(onevec);
          } else {
            censor1x = clone(zerovec);
          }
          
          censor2x = data[censorcol(0)];
          censor3x = data[censorcol(1)];
          
        } else {
        
          if (censorkeep == 0) {
            censor1x = clone(onevec);
          } else {
            censor1x = clone(zerovec);
          }
          
          censor2x = data[censorcol(0)];
          censor3x = data[censorcol(0)];
        }
      }
      
      if (aliveacol(0) > -1) {
        alivegiven1x = clone(zerovec);
        alivegiven2x = data[aliveacol(0)];
        alivegiven3x = data[aliveacol(1)];
      }
      
      if (deadacol(0) > -1) {
        deadgiven1x = clone(zerovec);
        deadgiven2x = data[deadacol(0)];
        deadgiven3x = data[deadacol(1)];
      }
      
      if (obsacol(0) > -1) {
        obsgiven1x = clone(zerovec);
        obsgiven2x = data[obsacol(0)];
        obsgiven3x = data[obsacol(1)];
      }
      
      if (nonobsacol(0) > -1) {
        nonobsgiven1x = negonevec;
        nonobsgiven2x = data[nonobsacol(0)];
        nonobsgiven3x = data[nonobsacol(1)];
      }
      
      if (juvcol(0) > -1) {
        juvgiven1x = clone(zerovec);
        juvgiven2x = data[juvcol(0)];
        juvgiven3x = data[juvcol(1)];
      }
      
    } else { // This portion of the if statement establishes what gets done in years after the 1st year
      
      sizea1x = data[sizeacol(j - 1)];
      sizea2x = data[sizeacol(j)];
      sizea3x = data[sizeacol(j + 1)];
      
      if (stagecol(0) > -1) {
        stage1x = data[stagecol(j - 1)];
        stage2x = data[stagecol(j)];
        stage3x = data[stagecol(j + 1)];
      }
      
      if (xcol(0) > -1) {
        xpos1x = data[xcol(j - 1)];
        xpos2x = data[xcol(j)];
        xpos3x = data[xcol(j + 1)];
      }
      
      if (ycol(0) > -1) {
        ypos1x = data[ycol(j - 1)];
        ypos2x = data[ycol(j)];
        ypos3x = data[ycol(j + 1)];
      }
      
      if (repstracol(0) > -1) {
        repstra1x = data[repstracol(j - 1)];
        repstra2x = data[repstracol(j)];
        repstra3x = data[repstracol(j + 1)];
      }
      
      if (fecacol(0) > -1) {
        feca1x = data[fecacol(j - 1)];
        feca2x = data[fecacol(j)];
        feca3x = data[fecacol(j + 1)];
      }
      
      if (sizebcol(0) > -1) {
        sizeb1x = data[sizebcol(j - 1)];
        sizeb2x = data[sizebcol(j)];
        sizeb3x = data[sizebcol(j + 1)];
      }
      
      if (repstrbcol(0) > -1) {
        repstrb1x = data[repstrbcol(j - 1)];
        repstrb2x = data[repstrbcol(j)];
        repstrb3x = data[repstrbcol(j + 1)];
      }
      
      if (fecbcol(0) > -1) {
        fecb1x = data[fecbcol(j - 1)];
        fecb2x = data[fecbcol(j)];
        fecb3x = data[fecbcol(j + 1)];
      }
      
      if (sizeccol(0) > -1) {
        sizec1x = data[sizeccol(j - 1)];
        sizec2x = data[sizeccol(j)];
        sizec3x = data[sizeccol(j + 1)];
      }
      
      if (indcovacol(0) > -1) {
        indcova1x = data[indcovacol(j - 1)];
        indcova2x = data[indcovacol(j)];
        indcova3x = data[indcovacol(j + 1)];
      }
      
      if (indcovbcol(0) > -1) {
        indcovb1x = data[indcovbcol(j - 1)];
        indcovb2x = data[indcovbcol(j)];
        indcovb3x = data[indcovbcol(j + 1)];
      }
      
      if (indcovccol(0) > -1) {
        indcovc1x = data[indcovccol(j - 1)];
        indcovc2x = data[indcovccol(j)];
        indcovc3x = data[indcovccol(j + 1)];
      }
      
      // This section decides how to deal with censor variables
      if (censorcol(0) > -1) {
        if (censrepeat) {
          censor1x = data[censorcol(j - 1)];
          censor2x = data[censorcol(j)];
          censor3x = data[censorcol(j + 1)];
        } else {
          censor1x = data[censorcol(0)];
          censor2x = data[censorcol(0)];
          censor3x = data[censorcol(0)];
        }
      }
      
      if (aliveacol(0) > -1) {
        alivegiven1x = data[aliveacol(j - 1)];
        alivegiven2x = data[aliveacol(j)];
        alivegiven3x = data[aliveacol(j + 1)];
      }
      
      if (deadacol(0) > -1) {
        deadgiven1x = data[deadacol(j - 1)];
        deadgiven2x = data[deadacol(j)];
        deadgiven3x = data[deadacol(j + 1)];
      }
      
      if (obsacol(0) > -1) {
        obsgiven1x = data[obsacol(j - 1)];
        obsgiven2x = data[obsacol(j)];
        obsgiven3x = data[obsacol(j + 1)];
      }
      
      if (nonobsacol(0) > -1) {
        nonobsgiven1x = data[nonobsacol(j - 1)];
        nonobsgiven2x = data[nonobsacol(j)];
        nonobsgiven3x = data[nonobsacol(j + 1)];
      }
      
      if (juvcol(0) > -1) {
        juvgiven1x = data[juvcol(j - 1)];
        juvgiven2x = data[juvcol(j)];
        juvgiven3x = data[juvcol(j + 1)];
      }
    }
    
    for (int i = 0; i < noindivs; i++) {
      livcheck1 = 0.0;
      livcheck2 = 0.0;
      livcheck3 = 0.0;
      
      nonobssub1 = -1.0;
      nonobssub2 = -1.0;
      nonobssub3 = -1.0;
      
      rowid[(i + (j * noindivs))] = static_cast<double>(i) + 1.0;
      popid[(i + (j * noindivs))] = popidx[i];
      patchid[(i + (j * noindivs))] = patchidx[i];
      individ[(i + (j * noindivs))] = individx[i];
      year2[(i + (j * noindivs))] = static_cast<double>(firstyear) + static_cast<double>(j);
      xpos1[(i + (j * noindivs))] = xpos1x[i];
      ypos1[(i + (j * noindivs))] = ypos1x[i];
      xpos2[(i + (j * noindivs))] = xpos2x[i];
      ypos2[(i + (j * noindivs))] = ypos2x[i];
      xpos3[(i + (j * noindivs))] = xpos3x[i];
      ypos3[(i + (j * noindivs))] = ypos3x[i];
      
      sizea1[(i + (j * noindivs))] = sizea1x[i];
      sizea2[(i + (j * noindivs))] = sizea2x[i];
      sizea3[(i + (j * noindivs))] = sizea3x[i];
      sizea10[(i + (j * noindivs))] = sizea1x[i];
      sizea20[(i + (j * noindivs))] = sizea2x[i];
      sizea30[(i + (j * noindivs))] = sizea3x[i];
      sizeb1[(i + (j * noindivs))] = sizeb1x[i];
      sizeb2[(i + (j * noindivs))] = sizeb2x[i];
      sizeb3[(i + (j * noindivs))] = sizeb3x[i];
      sizeb10[(i + (j * noindivs))] = sizeb1x[i];
      sizeb20[(i + (j * noindivs))] = sizeb2x[i];
      sizeb30[(i + (j * noindivs))] = sizeb3x[i];
      sizec1[(i + (j * noindivs))] = sizec1x[i];
      sizec2[(i + (j * noindivs))] = sizec2x[i];
      sizec3[(i + (j * noindivs))] = sizec3x[i];
      sizec10[(i + (j * noindivs))] = sizec1x[i];
      sizec20[(i + (j * noindivs))] = sizec2x[i];
      sizec30[(i + (j * noindivs))] = sizec3x[i];
      
      repstra1[(i + (j * noindivs))] = repstra1x[i];
      repstra2[(i + (j * noindivs))] = repstra2x[i];
      repstra3[(i + (j * noindivs))] = repstra3x[i];
      repstra10[(i + (j * noindivs))] = repstra1x[i];
      repstra20[(i + (j * noindivs))] = repstra2x[i];
      repstra30[(i + (j * noindivs))] = repstra3x[i];
      repstrb1[(i + (j * noindivs))] = repstrb1x[i];
      repstrb2[(i + (j * noindivs))] = repstrb2x[i];
      repstrb3[(i + (j * noindivs))] = repstrb3x[i];
      repstrb10[(i + (j * noindivs))] = repstrb1x[i];
      repstrb20[(i + (j * noindivs))] = repstrb2x[i];
      repstrb30[(i + (j * noindivs))] = repstrb3x[i];
      
      feca1[(i + (j * noindivs))] = feca1x[i];
      feca2[(i + (j * noindivs))] = feca2x[i];
      feca3[(i + (j * noindivs))] = feca3x[i];
      feca10[(i + (j * noindivs))] = feca1x[i];
      feca20[(i + (j * noindivs))] = feca2x[i];
      feca30[(i + (j * noindivs))] = feca3x[i];
      fecb1[(i + (j * noindivs))] = fecb1x[i];
      fecb2[(i + (j * noindivs))] = fecb2x[i];
      fecb3[(i + (j * noindivs))] = fecb3x[i];
      fecb10[(i + (j * noindivs))] = fecb1x[i];
      fecb20[(i + (j * noindivs))] = fecb2x[i];
      fecb30[(i + (j * noindivs))] = fecb3x[i];
      
      indcova1[(i + (j * noindivs))] = indcova1x[i];
      indcova2[(i + (j * noindivs))] = indcova2x[i];
      indcova3[(i + (j * noindivs))] = indcova3x[i];
      indcovb1[(i + (j * noindivs))] = indcovb1x[i];
      indcovb2[(i + (j * noindivs))] = indcovb2x[i];
      indcovb3[(i + (j * noindivs))] = indcovb3x[i];
      indcovc1[(i + (j * noindivs))] = indcovc1x[i];
      indcovc2[(i + (j * noindivs))] = indcovc2x[i];
      indcovc3[(i + (j * noindivs))] = indcovc3x[i];
      
      if (stagecol(0) > -1) {
        stage1[(i + (j * noindivs))] = stage1x[i];
        stage2[(i + (j * noindivs))] = stage2x[i];
        stage3[(i + (j * noindivs))] = stage3x[i];
        
        stage1num[(i + (j * noindivs))] = 0.0; // Should probably redo this to create stage indices...
        stage2num[(i + (j * noindivs))] = 0.0;
        stage3num[(i + (j * noindivs))] = 0.0;
      }
      
      if (censbool) { // Here we develop the censoring variable
        
        if (NumericVector::is_na(censor1x[i])) {
          censor1[(i + (j * noindivs))] = 0.0;
        } else if (j != 0) {
          censor1[(i + (j * noindivs))] = 1.0;
        }
        
        if (NumericVector::is_na(censor2x[i])) {
          censor2[(i + (j * noindivs))] = 0.0;
        } else {
          censor2[(i + (j * noindivs))] = 1.0;
        }
        
        if (NumericVector::is_na(censor3x[i])) {
          censor3[(i + (j * noindivs))] = 0.0;
        } else {
          censor3[(i + (j * noindivs))] = 1.0;
        }
      } else {
        censor1[(i + (j * noindivs))] = censor1x[i];
        censor2[(i + (j * noindivs))] = censor2x[i];
        censor3[(i + (j * noindivs))] = censor3x[i];
      }
      
      alivegiven1[(i + (j * noindivs))] = alivegiven1x[i];
      alivegiven2[(i + (j * noindivs))] = alivegiven2x[i];
      alivegiven3[(i + (j * noindivs))] = alivegiven3x[i];
      deadgiven1[(i + (j * noindivs))] = deadgiven1x[i];
      deadgiven2[(i + (j * noindivs))] = deadgiven2x[i];
      deadgiven3[(i + (j * noindivs))] = deadgiven3x[i];
      obsgiven1[(i + (j * noindivs))] = obsgiven1x[i];
      obsgiven2[(i + (j * noindivs))] = obsgiven2x[i];
      obsgiven3[(i + (j * noindivs))] = obsgiven3x[i];
      nonobsgiven1[(i + (j * noindivs))] = nonobsgiven1x[i];
      nonobsgiven2[(i + (j * noindivs))] = nonobsgiven2x[i];
      nonobsgiven3[(i + (j * noindivs))] = nonobsgiven3x[i];
      juvgiven1[(i + (j * noindivs))] = juvgiven1x[i];
      juvgiven2[(i + (j * noindivs))] = juvgiven2x[i];
      juvgiven3[(i + (j * noindivs))] = juvgiven3x[i];
      
      if (NumericVector::is_na(sizea10[(i + (j * noindivs))])) sizea10[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(sizea20[(i + (j * noindivs))])) sizea20[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(sizea30[(i + (j * noindivs))])) sizea30[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(sizeb10[(i + (j * noindivs))])) sizeb10[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(sizeb20[(i + (j * noindivs))])) sizeb20[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(sizeb30[(i + (j * noindivs))])) sizeb30[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(sizec10[(i + (j * noindivs))])) sizec10[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(sizec20[(i + (j * noindivs))])) sizec20[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(sizec30[(i + (j * noindivs))])) sizec30[(i + (j * noindivs))] = 0.0;
      
      if (NumericVector::is_na(repstra10[(i + (j * noindivs))])) repstra10[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(repstra20[(i + (j * noindivs))])) repstra20[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(repstra30[(i + (j * noindivs))])) repstra30[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(repstrb10[(i + (j * noindivs))])) repstrb10[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(repstrb20[(i + (j * noindivs))])) repstrb20[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(repstrb30[(i + (j * noindivs))])) repstrb30[(i + (j * noindivs))] = 0.0;
      
      if (NumericVector::is_na(feca10[(i + (j * noindivs))])) feca10[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(feca20[(i + (j * noindivs))])) feca20[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(feca30[(i + (j * noindivs))])) feca30[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(fecb10[(i + (j * noindivs))])) fecb10[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(fecb20[(i + (j * noindivs))])) fecb20[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(fecb30[(i + (j * noindivs))])) fecb30[(i + (j * noindivs))] = 0.0;
      
      if (NumericVector::is_na(indcova1[(i + (j * noindivs))])) indcova1[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(indcova2[(i + (j * noindivs))])) indcova2[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(indcova3[(i + (j * noindivs))])) indcova3[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(indcovb1[(i + (j * noindivs))])) indcovb1[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(indcovb2[(i + (j * noindivs))])) indcovb2[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(indcovb3[(i + (j * noindivs))])) indcovb3[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(indcovc1[(i + (j * noindivs))])) indcovc1[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(indcovc2[(i + (j * noindivs))])) indcovc2[(i + (j * noindivs))] = 0.0;
      if (NumericVector::is_na(indcovc3[(i + (j * noindivs))])) indcovc3[(i + (j * noindivs))] = 0.0;
      
      if (NumericVector::is_na(juvgiven1[(i + (j * noindivs))])) {
        juvgiven1[(i + (j * noindivs))] = 0.0;
      } else if (juvgiven1[(i + (j * noindivs))] != 0.0) {
        juvgiven1[(i + (j * noindivs))] = 1.0;
      }
      
      if (NumericVector::is_na(juvgiven2[(i + (j * noindivs))])) {
        juvgiven2[(i + (j * noindivs))] = 0.0;
      } else if (juvgiven2[(i + (j * noindivs))] != 0.0) {
        juvgiven2[(i + (j * noindivs))] = 1.0;
      }
      
      if (NumericVector::is_na(juvgiven3[(i + (j * noindivs))])) {
        juvgiven3[(i + (j * noindivs))] = 0.0;
      } else if (juvgiven3[(i + (j * noindivs))] != 0.0) {
        juvgiven3[(i + (j * noindivs))] = 1.0;
      }
      
      matstat1[(i + (j * noindivs))] = 1.0 - juvgiven1[(i + (j * noindivs))];
      matstat2[(i + (j * noindivs))] = 1.0 - juvgiven2[(i + (j * noindivs))];
      matstat3[(i + (j * noindivs))] = 1.0 - juvgiven3[(i + (j * noindivs))];
      
      addedsize1[(i + (j * noindivs))] = sizea10[(i + (j * noindivs))] + 
        sizeb10[(i + (j * noindivs))] + sizec10[(i + (j * noindivs))];
      addedsize2[(i + (j * noindivs))] = sizea20[(i + (j * noindivs))] + 
        sizeb20[(i + (j * noindivs))] + sizec20[(i + (j * noindivs))];
      addedsize3[(i + (j * noindivs))] = sizea30[(i + (j * noindivs))] + 
        sizeb30[(i + (j * noindivs))] + sizec30[(i + (j * noindivs))];
      
      addedrepstr1[(i + (j * noindivs))] = repstra10[(i + (j * noindivs))] + 
        (repstrb10[(i + (j * noindivs))] * repstrrel);
      addedrepstr2[(i + (j * noindivs))] = repstra20[(i + (j * noindivs))] + 
        (repstrb20[(i + (j * noindivs))] * repstrrel);
      addedrepstr3[(i + (j * noindivs))] = repstra30[(i + (j * noindivs))] + 
        (repstrb30[(i + (j * noindivs))] * repstrrel);
      
      addedfec1[(i + (j * noindivs))] = feca10[(i + (j * noindivs))] + 
        (fecb10[(i + (j * noindivs))] * fecrel);
      addedfec2[(i + (j * noindivs))] = feca20[(i + (j * noindivs))] + 
        (fecb20[(i + (j * noindivs))] * fecrel);
      addedfec3[(i + (j * noindivs))] = feca30[(i + (j * noindivs))] + 
        (fecb30[(i + (j * noindivs))] * fecrel);
      
      if (NumericVector::is_na(nonobsgiven1[(i + (j * noindivs))])) {
        nonobssub1 = -1.0;
      } else if (nonobsgiven1[(i + (j * noindivs))] >= 0) {
        nonobsgiven1[(i + (j * noindivs))] = 1.0;
        nonobssub1 = 1.0;
      }
      
      if (NumericVector::is_na(nonobsgiven2[(i + (j * noindivs))])) {
        nonobssub2 = -1.0;
      } else if (nonobsgiven2[(i + (j * noindivs))] >= 0) {
        nonobsgiven2[(i + (j * noindivs))] = 1.0;
        nonobssub2 = 1.0;
      }
      
      if (NumericVector::is_na(nonobsgiven3[(i + (j * noindivs))])) {
        nonobssub3 = -1.0;
      } else if (nonobsgiven3[(i + (j * noindivs))] >= 0) {
        nonobsgiven3[(i + (j * noindivs))] = 1.0;
        nonobssub3 = 1.0;
      }
      
      // Observation status
      if (addedsize1[(i + (j * noindivs))] > 0 || obsgiven1[(i + (j * noindivs))] > 0 || nonobssub1 == 0.0) {
        spryn1[(i + (j * noindivs))] = 1.0;
      }
      if (addedsize2[(i + (j * noindivs))] > 0 || obsgiven2[(i + (j * noindivs))] > 0 || nonobssub2 == 0.0) {
        spryn2[(i + (j * noindivs))] = 1.0;
      } 
      if (addedsize3[(i + (j * noindivs))] > 0 || obsgiven3[(i + (j * noindivs))] > 0 || nonobssub3 == 0.0) {
        spryn3[(i + (j * noindivs))] = 1.0;
      }
      
      // Here is the correction if size 0 individuals that reproduce are actually observed
      if (RepasObs == 1 && addedsize1[(i + (j * noindivs))] == 0 &&
          addedrepstr1[(i + (j * noindivs))] != 0.0) {
        spryn1[(i + (j * noindivs))] = 1.0;
      }
      if (RepasObs == 1 && addedsize2[(i + (j * noindivs))] == 0 &&
          addedrepstr2[(i + (j * noindivs))] != 0.0) {
        spryn2[(i + (j * noindivs))] = 1.0;
      }
      if (RepasObs == 1 && addedsize3[(i + (j * noindivs))] == 0 &&
          addedrepstr3[(i + (j * noindivs))] != 0.0) {
        spryn3[(i + (j * noindivs))] = 1.0;
      }
      
      if (addedrepstr1[(i + (j * noindivs))] > 0) repyn1[(i + (j * noindivs))] = 1.0;
      if (addedrepstr2[(i + (j * noindivs))] > 0) repyn2[(i + (j * noindivs))] = 1.0;
      if (addedrepstr3[(i + (j * noindivs))] > 0) repyn3[(i + (j * noindivs))] = 1.0;
      
      if (addedfec1[(i + (j * noindivs))] > 0) fecyn1[(i + (j * noindivs))] = 1.0;
      if (addedfec2[(i + (j * noindivs))] > 0) fecyn2[(i + (j * noindivs))] = 1.0;
      if (addedfec3[(i + (j * noindivs))] > 0) fecyn3[(i + (j * noindivs))] = 1.0;
      
      if (nonobssub1 >= 0) {
        livcheck1 = addedsize1[(i + (j * noindivs))] + addedrepstr1[(i + (j * noindivs))] + 
            spryn1[(i + (j * noindivs))] + nonobssub1;
      } else {
        livcheck1 = addedsize1[(i + (j * noindivs))] + addedrepstr1[(i + (j * noindivs))] + 
            spryn1[(i + (j * noindivs))];
      }
      if (nonobssub2 >= 0) {
        livcheck2 = addedsize2[(i + (j * noindivs))] + addedrepstr2[(i + (j * noindivs))] + 
            spryn2[(i + (j * noindivs))] + nonobssub2;
      } else {
        livcheck2 = addedsize2[(i + (j * noindivs))] + addedrepstr2[(i + (j * noindivs))] + 
            spryn2[(i + (j * noindivs))];
      }
      if (nonobssub3 >= 0) {
        livcheck3 = addedsize3[(i + (j * noindivs))] + addedrepstr3[(i + (j * noindivs))] + 
            spryn3[(i + (j * noindivs))] + nonobssub3;
      } else {
        livcheck3 = addedsize3[(i + (j * noindivs))] + addedrepstr3[(i + (j * noindivs))] + 
            spryn3[(i + (j * noindivs))];
      }
      
      if (livcheck1 > 0) alive1[(i + (j * noindivs))] = 1.0;
      if (livcheck2 > 0) alive2[(i + (j * noindivs))] = 1.0;
      if (livcheck3 > 0) alive3[(i + (j * noindivs))] = 1.0;
      
      if (alivegiven1[(i + (j * noindivs))] > 0) alive1[(i + (j * noindivs))] = 1.0;
      if (alivegiven2[(i + (j * noindivs))] > 0) alive2[(i + (j * noindivs))] = 1.0;
      if (alivegiven3[(i + (j * noindivs))] > 0) alive3[(i + (j * noindivs))] = 1.0;
      
      if (deadgiven1[(i + (j * noindivs))] > 0) alive1[(i + (j * noindivs))] = 0.0;
      if (deadgiven2[(i + (j * noindivs))] > 0) alive2[(i + (j * noindivs))] = 0.0;
      if (deadgiven3[(i + (j * noindivs))] > 0) alive3[(i + (j * noindivs))] = 0.0;
      
      // Corrections to living status, lifespan, and age based on later sightings
      if (firstseenx[i] <= (j + firstyear) && lastseenx[i] >= (j + firstyear)) { // This loop needs to be altered to check all years for lastseen adjustments, once firstseen is determined
        firstseen[(i + (j * noindivs))] = firstseenx[i];
        lastseen[(i + (j * noindivs))] = lastseenx[i];
        
        obsage[(i + (j * noindivs))] = year2[(i + (j * noindivs))] - firstseen[(i + (j * noindivs))];
        obslifespan[(i + (j * noindivs))] = lastseen[(i + (j * noindivs))] - firstseen[(i + (j * noindivs))];
        
        alive2[(i + (j * noindivs))] = 1;
        
        if (lastseenx[i] >= (j + firstyear + 1)) {
          alive3[(i + (j * noindivs))] = 1.0;
        }
        if (firstseenx[i] <= (j + firstyear - 1)) {
          alive1[(i + (j * noindivs))] = 1.0;
        }
      }
      
      // Here we take care of stage assignments
      if (stassign && stagecol(0) == -1) {
        if (stszcol == 8) {
          stagesize1 = sizea10[(i + (j * noindivs))];
          stagesize2 = sizea20[(i + (j * noindivs))];
          stagesize3 = sizea30[(i + (j * noindivs))];
          stagesize1b = sizeb10[(i + (j * noindivs))];
          stagesize2b = sizeb20[(i + (j * noindivs))];
          stagesize3b = sizeb30[(i + (j * noindivs))];
          stagesize1c = sizec10[(i + (j * noindivs))];
          stagesize2c = sizec20[(i + (j * noindivs))];
          stagesize3c = sizec30[(i + (j * noindivs))];
          
          arma::uvec stagemini1a = find(insfszminarma < stagesize1);
          arma::uvec stagemaxi1a = find(insfszmaxarma >= stagesize1);
          arma::uvec stagemini1b = find(insfszminarmab < stagesize1b);
          arma::uvec stagemaxi1b = find(insfszmaxarmab >= stagesize1b);
          arma::uvec stagemini1c = find(insfszminarmac < stagesize1c);
          arma::uvec stagemaxi1c = find(insfszmaxarmac >= stagesize1c);
          arma::uvec stagemini2a = find(insfszminarma < stagesize2);
          arma::uvec stagemaxi2a = find(insfszmaxarma >= stagesize2);
          arma::uvec stagemini2b = find(insfszminarmab < stagesize2b);
          arma::uvec stagemaxi2b = find(insfszmaxarmab >= stagesize2b);
          arma::uvec stagemini2c = find(insfszminarmac < stagesize2c);
          arma::uvec stagemaxi2c = find(insfszmaxarmac >= stagesize2c);
          arma::uvec stagemini3a = find(insfszminarma < stagesize3);
          arma::uvec stagemaxi3a = find(insfszmaxarma >= stagesize3);
          arma::uvec stagemini3b = find(insfszminarmab < stagesize3b);
          arma::uvec stagemaxi3b = find(insfszmaxarmab >= stagesize3b);
          arma::uvec stagemini3c = find(insfszminarmac < stagesize3c);
          arma::uvec stagemaxi3c = find(insfszmaxarmac >= stagesize3c);
          
          arma::uvec stagemini1d = intersect(stagemini1a, stagemini1b);
          stagemini1 = intersect(stagemini1c, stagemini1d);
          arma::uvec stagemaxi1d = intersect(stagemaxi1a, stagemaxi1b);
          stagemaxi1 = intersect(stagemaxi1c, stagemaxi1d);
          arma::uvec stagemini2d = intersect(stagemini2a, stagemini2b);
          stagemini2 = intersect(stagemini2c, stagemini2d);
          arma::uvec stagemaxi2d = intersect(stagemaxi2a, stagemaxi2b);
          stagemaxi2 = intersect(stagemaxi2c, stagemaxi2d);
          arma::uvec stagemini3d = intersect(stagemini3a, stagemini3b);
          stagemini3 = intersect(stagemini3c, stagemini3d);
          arma::uvec stagemaxi3d = intersect(stagemaxi3a, stagemaxi3b);
          stagemaxi3 = intersect(stagemaxi3c, stagemaxi3d);
          
        } else if (stszcol == 7) {
          stagesize1 = sizeb10[(i + (j * noindivs))];
          stagesize2 = sizeb20[(i + (j * noindivs))];
          stagesize3 = sizeb30[(i + (j * noindivs))];
          stagesize1b = sizec10[(i + (j * noindivs))];
          stagesize2b = sizec20[(i + (j * noindivs))];
          stagesize3b = sizec30[(i + (j * noindivs))];
          
          arma::uvec stagemini1a = find(insfszminarma < stagesize1);
          arma::uvec stagemaxi1a = find(insfszmaxarma >= stagesize1);
          arma::uvec stagemini1b = find(insfszminarmab < stagesize1b);
          arma::uvec stagemaxi1b = find(insfszmaxarmab >= stagesize1b);
          arma::uvec stagemini2a = find(insfszminarma < stagesize2);
          arma::uvec stagemaxi2a = find(insfszmaxarma >= stagesize2);
          arma::uvec stagemini2b = find(insfszminarmab < stagesize2b);
          arma::uvec stagemaxi2b = find(insfszmaxarmab >= stagesize2b);
          arma::uvec stagemini3a = find(insfszminarma < stagesize3);
          arma::uvec stagemaxi3a = find(insfszmaxarma >= stagesize3);
          arma::uvec stagemini3b = find(insfszminarmab < stagesize3b);
          arma::uvec stagemaxi3b = find(insfszmaxarmab >= stagesize3b);
          
          stagemini1 = intersect(stagemini1a, stagemini1b);
          stagemaxi1 = intersect(stagemaxi1a, stagemaxi1b);
          stagemini2 = intersect(stagemini2a, stagemini2b);
          stagemaxi2 = intersect(stagemaxi2a, stagemaxi2b);
          stagemini3 = intersect(stagemini3a, stagemini3b);
          stagemaxi3 = intersect(stagemaxi3a, stagemaxi3b);
          
        } else if (stszcol == 6) {
          stagesize1 = sizea10[(i + (j * noindivs))];
          stagesize2 = sizea20[(i + (j * noindivs))];
          stagesize3 = sizea30[(i + (j * noindivs))];
          stagesize1b = sizec10[(i + (j * noindivs))];
          stagesize2b = sizec20[(i + (j * noindivs))];
          stagesize3b = sizec30[(i + (j * noindivs))];
          
          arma::uvec stagemini1a = find(insfszminarma < stagesize1);
          arma::uvec stagemaxi1a = find(insfszmaxarma >= stagesize1);
          arma::uvec stagemini1b = find(insfszminarmab < stagesize1b);
          arma::uvec stagemaxi1b = find(insfszmaxarmab >= stagesize1b);
          arma::uvec stagemini2a = find(insfszminarma < stagesize2);
          arma::uvec stagemaxi2a = find(insfszmaxarma >= stagesize2);
          arma::uvec stagemini2b = find(insfszminarmab < stagesize2b);
          arma::uvec stagemaxi2b = find(insfszmaxarmab >= stagesize2b);
          arma::uvec stagemini3a = find(insfszminarma < stagesize3);
          arma::uvec stagemaxi3a = find(insfszmaxarma >= stagesize3);
          arma::uvec stagemini3b = find(insfszminarmab < stagesize3b);
          arma::uvec stagemaxi3b = find(insfszmaxarmab >= stagesize3b);
          
          stagemini1 = intersect(stagemini1a, stagemini1b);
          stagemaxi1 = intersect(stagemaxi1a, stagemaxi1b);
          stagemini2 = intersect(stagemini2a, stagemini2b);
          stagemaxi2 = intersect(stagemaxi2a, stagemaxi2b);
          stagemini3 = intersect(stagemini3a, stagemini3b);
          stagemaxi3 = intersect(stagemaxi3a, stagemaxi3b);
          
        } else if (stszcol == 5) {
          stagesize1 = sizea10[(i + (j * noindivs))];
          stagesize2 = sizea20[(i + (j * noindivs))];
          stagesize3 = sizea30[(i + (j * noindivs))];
          stagesize1b = sizeb10[(i + (j * noindivs))];
          stagesize2b = sizeb20[(i + (j * noindivs))];
          stagesize3b = sizeb30[(i + (j * noindivs))];
          
          arma::uvec stagemini1a = find(insfszminarma < stagesize1);
          arma::uvec stagemaxi1a = find(insfszmaxarma >= stagesize1);
          arma::uvec stagemini1b = find(insfszminarmab < stagesize1b);
          arma::uvec stagemaxi1b = find(insfszmaxarmab >= stagesize1b);
          arma::uvec stagemini2a = find(insfszminarma < stagesize2);
          arma::uvec stagemaxi2a = find(insfszmaxarma >= stagesize2);
          arma::uvec stagemini2b = find(insfszminarmab < stagesize2b);
          arma::uvec stagemaxi2b = find(insfszmaxarmab >= stagesize2b);
          arma::uvec stagemini3a = find(insfszminarma < stagesize3);
          arma::uvec stagemaxi3a = find(insfszmaxarma >= stagesize3);
          arma::uvec stagemini3b = find(insfszminarmab < stagesize3b);
          arma::uvec stagemaxi3b = find(insfszmaxarmab >= stagesize3b);
          
          stagemini1 = intersect(stagemini1a, stagemini1b);
          stagemaxi1 = intersect(stagemaxi1a, stagemaxi1b);
          stagemini2 = intersect(stagemini2a, stagemini2b);
          stagemaxi2 = intersect(stagemaxi2a, stagemaxi2b);
          stagemini3 = intersect(stagemini3a, stagemini3b);
          stagemaxi3 = intersect(stagemaxi3a, stagemaxi3b);
          
        } else if (stszcol == 4) {
          stagesize1 = addedsize1[(i + (j * noindivs))];
          stagesize2 = addedsize2[(i + (j * noindivs))];
          stagesize3 = addedsize3[(i + (j * noindivs))];
          
          stagemini1 = find(insfszminarma < stagesize1);
          stagemaxi1 = find(insfszmaxarma >= stagesize1);
          stagemini2 = find(insfszminarma < stagesize2);
          stagemaxi2 = find(insfszmaxarma >= stagesize2);
          stagemini3 = find(insfszminarma < stagesize3);
          stagemaxi3 = find(insfszmaxarma >= stagesize3);
          
        } else if (stszcol == 3) {
          stagesize1 = sizec10[(i + (j * noindivs))];
          stagesize2 = sizec20[(i + (j * noindivs))];
          stagesize3 = sizec30[(i + (j * noindivs))];
          
          stagemini1 = find(insfszminarma < stagesize1);
          stagemaxi1 = find(insfszmaxarma >= stagesize1);
          stagemini2 = find(insfszminarma < stagesize2);
          stagemaxi2 = find(insfszmaxarma >= stagesize2);
          stagemini3 = find(insfszminarma < stagesize3);
          stagemaxi3 = find(insfszmaxarma >= stagesize3);
          
        } else if (stszcol == 2) {
          stagesize1 = sizeb10[(i + (j * noindivs))];
          stagesize2 = sizeb20[(i + (j * noindivs))];
          stagesize3 = sizeb30[(i + (j * noindivs))];
          
          stagemini1 = find(insfszminarma < stagesize1);
          stagemaxi1 = find(insfszmaxarma >= stagesize1);
          stagemini2 = find(insfszminarma < stagesize2);
          stagemaxi2 = find(insfszmaxarma >= stagesize2);
          stagemini3 = find(insfszminarma < stagesize3);
          stagemaxi3 = find(insfszmaxarma >= stagesize3);
          
        } else {
          stagesize1 = sizea10[(i + (j * noindivs))];
          stagesize2 = sizea20[(i + (j * noindivs))];
          stagesize3 = sizea30[(i + (j * noindivs))];
          
          stagemini1 = find(insfszminarma < stagesize1);
          stagemaxi1 = find(insfszmaxarma >= stagesize1);
          stagemini2 = find(insfszminarma < stagesize2);
          stagemaxi2 = find(insfszmaxarma >= stagesize2);
          stagemini3 = find(insfszminarma < stagesize3);
          stagemaxi3 = find(insfszmaxarma >= stagesize3);
        }
        
        // Stage 2
        stagerep = find(inrepstatarma == repyn2[(i + (j * noindivs))]);
        stagemat = find(inmatstatarma == matstat2[(i + (j * noindivs))]);
        stageobs = find(inobsstatarma == spryn2[(i + (j * noindivs))]);
        
        cs1 = intersect(stagemini2, stagemaxi2);
        
        if (NRasRep && NOasObs) {
          cs4 = cs1;
        } else if (NOasObs) {
          cs4 = intersect(stagerep, cs1);
        } else if (NRasRep) {
          cs2 = intersect(stageobs, stagemat);
          cs4 = intersect(cs1, cs2);
        } else {
          cs2 = intersect(stageobs, stagemat);
          cs3 = intersect(stagerep, cs2);
          cs4 = intersect(cs1, cs3);
        }
        
        if (cs4.n_elem > 0 && alive2[(i + (j * noindivs))] == 1) {
          choicestage = instageid(cs4[0]) - 1;
          stage2num[(i + (j * noindivs))] = static_cast<double>(choicestage) + 1.0;
          
          stage2[(i + (j * noindivs))] = sfname[choicestage];
        } else {
          stage2[(i + (j * noindivs))] = "NotAlive";
        }
        
        // Stage 1
        stagerep = find(inrepstatarma == repyn1[(i + (j * noindivs))]);
        stagemat = find(inmatstatarma == matstat1[(i + (j * noindivs))]);
        stageobs = find(inobsstatarma == spryn1[(i + (j * noindivs))]);
        
        cs1 = intersect(stagemini1, stagemaxi1);
        
        if (NRasRep && NOasObs) {
          cs4 = cs1;
        } else if (NOasObs) {
          cs4 = intersect(stagerep, cs1);
        } else if (NRasRep) {
          cs2 = intersect(stageobs, stagemat);
          cs4 = intersect(cs1, cs2);
        } else {
          cs2 = intersect(stageobs, stagemat);
          cs3 = intersect(stagerep, cs2);
          cs4 = intersect(cs1, cs3);
        }
        
        if (cs4.n_elem > 0 && alive1[(i + (j * noindivs))] == 1) {
          choicestage = instageid(cs4[0]) - 1;
          stage1num[(i + (j * noindivs))] = static_cast<double>(choicestage) + 1.0;
          
          stage1[(i + (j * noindivs))] = sfname[choicestage];
        } else {
          stage1[(i + (j * noindivs))] = "NotAlive";
          matstat1[(i + (j * noindivs))] = 0.0;
        }
        
        // Stage 3
        stagerep = find(inrepstatarma == repyn3[(i + (j * noindivs))]);
        stagemat = find(inmatstatarma == matstat3[(i + (j * noindivs))]);
        stageobs = find(inobsstatarma == spryn3[(i + (j * noindivs))]);
        
        cs1 = intersect(stagemini3, stagemaxi3);
        
        if (NRasRep && NOasObs) {
          cs4 = cs1;
        } else if (NOasObs) {
          cs4 = intersect(stagerep, cs1);
        } else if (NRasRep) {
          cs2 = intersect(stageobs, stagemat);
          cs4 = intersect(cs1, cs2);
        } else {
          cs2 = intersect(stageobs, stagemat);
          cs3 = intersect(stagerep, cs2);
          cs4 = intersect(cs1, cs3);
        }
        
        // Here we create exceptions based on stage assignment problems in time t+1
        if (cs4.n_elem == 1 && alive3[(i + (j * noindivs))] == 1) {
          choicestage = instageid(cs4[0]) - 1;
          stage3num[(i + (j * noindivs))] = static_cast<double>(choicestage) + 1.0;
          
          stage3[(i + (j * noindivs))] = sfname[choicestage];
          
        } else if (alive3[(i + (j * noindivs))] != 1) {
          stage3[(i + (j * noindivs))] = "NotAlive";
          
        } else if (cs4.n_elem == 0) {
          stage3[(i + (j * noindivs))] = "NoMatch";
          
          if (!quiet) Rf_warningcall(R_NilValue, 
              "Some stages occurring in the dataset do not match any characteristics in the input stageframe.");
          
        } else if (cs4.n_elem > 1) {
          if (!quiet) Rf_warningcall(R_NilValue, 
              "Some stages in the input stageframe appear to have the same description. Please make sure that all stages included in the stageframe are defined with unique sets of characteristics.");
          
        } else {
          throw Rcpp::exception("Stage assignment error.", false);
        }
      } // stassign if statement
    } // i loop
  } // j loop
  
  if (NAas0) {
    sizea1 = sizea10;
    sizea2 = sizea20;
    sizea3 = sizea30;
    
    sizeb1 = sizeb10;
    sizeb2 = sizeb20;
    sizeb3 = sizeb30;
    
    sizec1 = sizec10;
    sizec2 = sizec20;
    sizec3 = sizec30;
    
    repstra1 = repstra10;
    repstra2 = repstra20;
    repstra3 = repstra30;
    
    repstrb1 = repstrb10;
    repstrb2 = repstrb20;
    repstrb3 = repstrb30;
    
    feca1 = feca10;
    feca2 = feca20;
    feca3 = feca30;
    
    fecb1 = fecb10;
    fecb2 = fecb20;
    fecb3 = fecb30;
  }
  
  output_longlist(0) = rowid;
  output_longlist(1) = popid;
  output_longlist(2) = patchid;
  output_longlist(3) = individ;
  output_longlist(4) = year2;
  output_longlist(5) = firstseen;
  output_longlist(6) = lastseen;
  output_longlist(7) = obsage;
  output_longlist(8) = obslifespan;
  
  output_longlist(9) = xpos1;
  output_longlist(10) = ypos1;
  output_longlist(11) = sizea1;
  output_longlist(12) = sizeb1;
  output_longlist(13) = sizec1;
  output_longlist(14) = addedsize1;
  output_longlist(15) = repstra1;
  output_longlist(16) = repstrb1;
  output_longlist(17) = addedrepstr1;
  output_longlist(18) = feca1;
  output_longlist(19) = fecb1;
  output_longlist(20) = addedfec1;
  output_longlist(21) = indcova1;
  output_longlist(22) = indcovb1;
  output_longlist(23) = indcovc1;
  output_longlist(24) = censor1;
  output_longlist(25) = juvgiven1;
  
  output_longlist(26) = spryn1;
  output_longlist(27) = repyn1;
  output_longlist(28) = fecyn1;
  output_longlist(29) = matstat1;
  output_longlist(30) = alive1;
  output_longlist(31) = stage1;
  output_longlist(32) = stage1num;
  
  output_longlist(33) = xpos2;
  output_longlist(34) = ypos2;
  output_longlist(35) = sizea2;
  output_longlist(36) = sizeb2;
  output_longlist(37) = sizec2;
  output_longlist(38) = addedsize2;
  output_longlist(39) = repstra2;
  output_longlist(40) = repstrb2;
  output_longlist(41) = addedrepstr2;
  output_longlist(42) = feca2;
  output_longlist(43) = fecb2;
  output_longlist(44) = addedfec2;
  output_longlist(45) = indcova2;
  output_longlist(46) = indcovb2;
  output_longlist(47) = indcovc2;
  output_longlist(48) = censor2;
  output_longlist(49) = juvgiven2;
  
  output_longlist(50) = spryn2;
  output_longlist(51) = repyn2;
  output_longlist(52) = fecyn2;
  output_longlist(53) = matstat2;
  output_longlist(54) = alive2;
  output_longlist(55) = stage2;
  output_longlist(56) = stage2num;

  output_longlist(57) = xpos3;
  output_longlist(58) = ypos3;
  output_longlist(59) = sizea3;
  output_longlist(60) = sizeb3;
  output_longlist(61) = sizec3;
  output_longlist(62) = addedsize3;
  output_longlist(63) = repstra3;
  output_longlist(64) = repstrb3;
  output_longlist(65) = addedrepstr3;
  output_longlist(66) = feca3;
  output_longlist(67) = fecb3;
  output_longlist(68) = addedfec3;
  output_longlist(69) = indcova3;
  output_longlist(70) = indcovb3;
  output_longlist(71) = indcovc3;
  output_longlist(72) = censor3;
  output_longlist(73) = juvgiven3;
  
  output_longlist(74) = spryn3;
  output_longlist(75) = repyn3;
  output_longlist(76) = fecyn3;
  output_longlist(77) = matstat3;
  output_longlist(78) = alive3;
  output_longlist(79) = stage3;
  output_longlist(80) = stage3num;
  
  Rcpp::CharacterVector varnames {"rowid", "popid", "patchid", "individ",
    "year2", "firstseen", "lastseen", "obsage", "obslifespan","xpos1", "ypos1",
    "sizea1", "sizeb1", "sizec1", "size1added", "repstra1", "repstrb1",
    "repstr1added", "feca1", "fecb1", "fec1added", "indcova1", "indcovb1",
    "indcovc1", "censor1", "juvgiven1", "obsstatus1", "repstatus1",
    "fecstatus1", "matstatus1", "alive1", "stage1", "stage1index",
    "xpos2", "ypos2", "sizea2", "sizeb2", "sizec2", "size2added", "repstra2",
    "repstrb2", "repstr2added", "feca2", "fecb2", "fec2added", "indcova2",
    "indcovb2", "indcovc2", "censor2", "juvgiven2", "obsstatus2", "repstatus2",
    "fecstatus2", "matstatus2", "alive2", "stage2", "stage2index",
    "xpos3", "ypos3", "sizea3", "sizeb3", "sizec3", "size3added", "repstra3",
    "repstrb3", "repstr3added", "feca3", "fecb3", "fec3added", "indcova3",
    "indcovb3", "indcovc3", "censor3", "juvgiven3", "obsstatus3", "repstatus3",
    "fecstatus3", "matstatus3", "alive3", "stage3", "stage3index"};
  
  output_longlist.attr("names") = varnames;
  output_longlist.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, alive3.length());
  StringVector needed_classes {"data.frame", "hfvdata"};
  output_longlist.attr("class") = needed_classes; // data.frame
  
  return output_longlist;
}

//' Create Historical Vertical Structure for Ahistorical Vertical Data Frame
//' 
//' Function \code{jpf()} is the core kernel for function
//' \code{\link{historicalize3}()}, creating the historical, vertical structure
//' and rearranging the data in that shape.
//'
//' @param data The horizontal data file.
//' @param stageframe The stageframe object identifying the life history model
//' being operationalized. This should be the full stageframe.
//' @param popidcol Column number corresponding to the identity of the
//' population for each individual.
//' @param patchidcol Column number corresponding to the identity of the patch
//' for each individual.
//' @param individcol Column number corresponding to the identity of each 
//' individual.
//' @param year2col Column number of year or occasion in time \emph{t}.
//' @param year3col Column number of year or occasion in time \emph{t}+1.
//' @param xcol Column number corresponding to the x coordinate of each
//' individual in Cartesian space.
//' @param ycol Column number corresponding to the y coordinate of each
//' individual in Cartesian space.
//' @param juv2col Column number coding for status as a juvenile in time
//' \emph{t}.
//' @param juv3col Column number coding for status as a juvenile in time
//' \emph{t}+1.
//' @param sizea2col Column number corresponding to the primary size variable in
//' time \emph{t}.
//' @param sizea3col Column number corresponding to the primary size variable in
//' time \emph{t}+1.
//' @param sizeb2col Column number corresponding to the secondary size variable
//' in time \emph{t}.
//' @param sizeb3col Column number corresponding to the secondary size variable
//' in time \emph{t}+1.
//' @param sizec2col Column number corresponding to the tertiary size variable
//' in time \emph{t}.
//' @param sizec3col Column number corresponding to the tertiary size variable
//' in time \emph{t}+1.
//' @param repstra2col Column number corresponding to the main variable coding
//' for the production of reproductive structures, such as flowers, in time
//' \emph{t}.
//' @param repstra3col Column number corresponding to the main variable coding
//' for the production of reproductive structures, such as flowers, in time
//' \emph{t}+1.
//' @param repstrb2col Column number corresponding to a second variable coding
//' for the production of reproductive structures, such as flowers, in time
//' \emph{t}.
//' @param repstrb3col Column number corresponding to a second variable coding
//' for the production of reproductive structures, such as flowers, in time
//' \emph{t}+1.
//' @param feca2col Column number corresponding to the main variable coding for
//' fecundity in time \emph{t}.
//' @param feca3col Column number corresponding to the main variable coding for
//' fecundity in time \emph{t}+1.
//' @param fecb2col Column number corresponding to a second variable coding for
//' fecundity in time \emph{t}.
//' @param fecb3col Column number corresponding to a second variable coding for
//' fecundity in time \emph{t}+1.
//' @param indcova2col Column number corresponding to an individual covariate in
//' time \emph{t}.
//' @param indcova3col Column number corresponding to an individual covariate in
//' time \emph{t}+1.
//' @param indcovb2col Column number corresponding to an individual covariate in
//' time \emph{t}.
//' @param indcovb3col Column number corresponding to an individual covariate in
//' time \emph{t}+1.
//' @param indcovc2col Column number corresponding to an individual covariate in
//' time \emph{t}.
//' @param indcovc3col Column number corresponding to an individual covariate in
//' time \emph{t}+1.
//' @param alive2col Column number detailing whether an individual is alive in 
//' time \emph{t}.
//' @param alive3col Column number detailing whether an individual is alive in 
//' time \emph{t}+1.
//' @param dead2col Column number detailing whether an individual is dead in 
//' time \emph{t}.
//' @param dead3col Column number detailing whether an individual is dead in 
//' time \emph{t}+1.
//' @param obs2col Column number detailing whether an individual is in an
//' observable stage in time \emph{t}.
//' @param obs3col Column number detailing whether an individual is in an
//' observable stage in time \emph{t}+1.
//' @param nonobs2col Column number detailing whether an individual is in an
//' unobservable stage in time \emph{t}.
//' @param nonobs3col Column number detailing whether an individual is in an
//' unobservable stage in time \emph{t}+1.
//' @param repstrrel This is a scalar multiplier for that makes the variable in
//' \code{repstrb2col} equivalent to \code{repstra2col}.
//' @param fecrel This is a scalar multiplier for that makes the variable in
//' \code{fecb2col} equivalent to \code{feca2col}.
//' @param stage2col Column number corresponding to life history stage in time
//' \emph{t}.
//' @param stage3col Column number corresponding to life history stage in time
//' \emph{t}+1.
//' @param censorcol Column number corresponding to a censor variable within the
//' dataset.
//' @param NAas0 If \code{TRUE}, then all \code{NA} entries for size and
//' fecundity variables will be set to 0.
//' @param NRasRep If \code{TRUE}, then non-reproductive but mature individuals
//' will be treated as reproductive during stage assignment.
//' @param NOasObs If TRUE, then will treat unobserved individuals as observed
//' during stage assignment.
//' @param stassign A logical value indicating whether to assign stages.
//' @param stszcol Integer describing which size variable to use in stage 
//' estimation. Numbers 1 through 8 are possible.
//' @param censorkeep Numeric value of censor variable, denoting elements to
//' keep. If \code{NA} is to be used, then set this variable to \code{0} and set
//' \code{censbool = TRUE}.
//' @param censbool A logical variable determining whether \code{NA} denotes the
//' value of the censoring variable identifying data to keep.
//' @param quiet A logical value indicating whether to silense warnings.
//' 
//' @return The output is currently a list coerced into the data frame class,
//' and is read as a data frame by R. It is secondarily in class \code{hfvdata}.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.jpf)]]
Rcpp::List jpf(DataFrame data, DataFrame stageframe, int popidcol,
  int patchidcol, int individcol, int year2col, int year3col, int xcol,
  int ycol, int juv2col, int juv3col, int sizea2col, int sizea3col,
  int sizeb2col, int sizeb3col, int sizec2col, int sizec3col, int repstra2col,
  int repstra3col, int repstrb2col, int repstrb3col, int feca2col, int feca3col,
  int fecb2col, int fecb3col, int indcova2col, int indcova3col, int indcovb2col,
  int indcovb3col, int indcovc2col, int indcovc3col, int alive2col,
  int alive3col, int dead2col, int dead3col, int obs2col, int obs3col,
  int nonobs2col, int nonobs3col, double repstrrel, double fecrel,
  int stage2col, int stage3col, int censorcol, bool NAas0, bool NRasRep,
  bool NOasObs, bool stassign, int stszcol, double censorkeep, bool censbool,
  bool quiet) {
  
  Rcpp::List output_longlist(81);
  
  int norows = data.nrows();
  
  Rcpp::NumericVector ckcheck (1);
  ckcheck(0) = censorkeep;
  double crazycensor;
  Rcpp::NumericVector censfillvec (norows);
  if (NumericVector::is_na(ckcheck(0))) {
    crazycensor = 0;
  } else{
    crazycensor = censorkeep;
  }
  
  Rcpp::NumericVector zerovec (norows);
  Rcpp::NumericVector negonevec (norows);
  zerovec.fill(0);
  negonevec.fill(-1);
  
  Rcpp::StringVector individx(norows);
  individx = data[individcol];
  Rcpp::StringVector allindivs = unique(individx);
  int noindivs = allindivs.size(); // Total number of individuals in dataset
  
  Rcpp::IntegerVector year3x;
  Rcpp::IntegerVector year2x = as<IntegerVector>(data[year2col]);
  if (year3col != -1) {
    year3x = as<IntegerVector>(data[year3col]);
  } else {year3x = zerovec;}
  
  Rcpp::IntegerVector yearall2x = sort_unique(year2x);
  int firstyear = min(yearall2x);
  int noyears = yearall2x.size(); // Total number of observation periods
  
  int ndflength = noyears * noindivs; // Initial length of final hfv dataset
  int currentyear {0};
  int currentindiv {-1};
  int ndfindex {0};
  int prevyrindex {0};
  int nextyrindex {0};
  double livcheck2 {0.0};
  double livcheck3 {0.0};
  
  Rcpp::StringVector sfname = as<StringVector>(stageframe["stage"]);
  Rcpp::NumericVector repstat = as<NumericVector>(stageframe["repstatus"]);
  Rcpp::NumericVector obsstat = as<NumericVector>(stageframe["obsstatus"]);
  Rcpp::NumericVector matstat = as<NumericVector>(stageframe["matstatus"]);
  arma::vec indataset = as<arma::vec>(stageframe["indataset"]);
  arma::vec sfszmin = as<arma::vec>(stageframe["sizebin_min"]);
  arma::vec sfszmax = as<arma::vec>(stageframe["sizebin_max"]);
  arma::vec sfszminb = as<arma::vec>(stageframe["sizebinb_min"]);
  arma::vec sfszmaxb = as<arma::vec>(stageframe["sizebinb_max"]);
  arma::vec sfszminc = as<arma::vec>(stageframe["sizebinc_min"]);
  arma::vec sfszmaxc = as<arma::vec>(stageframe["sizebinc_max"]);
  
  arma::vec repstatarma = as<arma::vec>(repstat);
  arma::vec obsstatarma = as<arma::vec>(obsstat);
  arma::vec matstatarma = as<arma::vec>(matstat);
  arma::vec sfszminarma = sfszmin;
  arma::vec sfszmaxarma = sfszmax;
  arma::vec sfszminarmab = sfszminb;
  arma::vec sfszmaxarmab = sfszmaxb;
  arma::vec sfszminarmac = sfszminc;
  arma::vec sfszmaxarmac = sfszmaxc;
  int stagenum = sfszmaxarma.n_elem; // Total number of life stages in stageframe
  
  arma::uvec instages = find(indataset == 1); 
  int instagenum = instages.n_elem; // Total number of stages in dataset
  
  arma::uvec stageid(stagenum);
  arma::uvec instageid(instagenum);
  arma::vec insfszminarma(instagenum);
  arma::vec insfszmaxarma(instagenum);
  arma::vec inrepstatarma(instagenum);
  arma::vec inobsstatarma(instagenum);
  arma::vec inmatstatarma(instagenum);
  arma::vec insfszminarmab (instagenum);
  arma::vec insfszmaxarmab (instagenum);
  arma::vec insfszminarmac (instagenum);
  arma::vec insfszmaxarmac (instagenum);
  
  // This section creates vectors describing only life history stages in the dataset
  int inplace {0};
  for (int i = 0; i < stagenum; i++) {
    stageid(i) = i+1;
    
    if (indataset(i) == 1) {
      instageid(inplace) = i + 1;
      insfszminarma(inplace) = sfszminarma(i);
      insfszmaxarma(inplace) = sfszmaxarma(i);
      insfszminarmab(inplace) = sfszminarmab(i);
      insfszmaxarmab(inplace) = sfszmaxarmab(i);
      insfszminarmac(inplace) = sfszminarmac(i);
      insfszmaxarmac(inplace) = sfszmaxarmac(i);
      inrepstatarma(inplace) = repstatarma(i);
      inobsstatarma(inplace) = obsstatarma(i);
      inmatstatarma(inplace) = matstatarma(i);
      
      inplace++;
    }
  }
  
  // The following set of variable definitions is different from those used in pfj.
  // Here, these variables are defined by the row structure of the data, whereas in
  // pfj they correspond to the individuals in the data frame.
  Rcpp::StringVector popidx;
  Rcpp::StringVector patchidx;
  Rcpp::NumericVector xpos2x;
  Rcpp::NumericVector ypos2x;
  Rcpp::NumericVector sizea2x;
  Rcpp::NumericVector sizea3x;
  Rcpp::NumericVector repstra2x;
  Rcpp::NumericVector repstra3x;
  Rcpp::NumericVector feca2x;
  Rcpp::NumericVector feca3x;
  Rcpp::NumericVector sizeb2x;
  Rcpp::NumericVector sizeb3x;
  Rcpp::NumericVector repstrb2x;
  Rcpp::NumericVector repstrb3x;
  Rcpp::NumericVector fecb2x;
  Rcpp::NumericVector fecb3x;
  Rcpp::NumericVector sizec2x;
  Rcpp::NumericVector sizec3x;
  
  Rcpp::NumericVector indcova2x;
  Rcpp::NumericVector indcova3x;
  Rcpp::NumericVector indcovb2x;
  Rcpp::NumericVector indcovb3x;
  Rcpp::NumericVector indcovc2x;
  Rcpp::NumericVector indcovc3x;
  
  Rcpp::NumericVector censor2x;
  Rcpp::NumericVector alivegiven2x;
  Rcpp::NumericVector alivegiven3x;
  Rcpp::NumericVector deadgiven2x;
  Rcpp::NumericVector deadgiven3x;
  Rcpp::NumericVector obsgiven2x;
  Rcpp::NumericVector obsgiven3x;
  Rcpp::NumericVector nonobsgiven2x; 
  Rcpp::NumericVector nonobsgiven3x;
  Rcpp::NumericVector juvgiven2x;
  Rcpp::NumericVector juvgiven3x;
  nonobsgiven3x.fill(-1);
  
  Rcpp::NumericVector firstseenx (noindivs);
  Rcpp::NumericVector lastseenx (noindivs);
  firstseenx.fill(-1);
  lastseenx.fill(-1);
  
  Rcpp::NumericVector sizea20x (norows);
  Rcpp::NumericVector sizeb20x (norows);
  Rcpp::NumericVector sizec20x (norows);
  Rcpp::NumericVector repstra20x (norows);
  Rcpp::NumericVector repstrb20x (norows);
  Rcpp::NumericVector feca20x (norows);
  Rcpp::NumericVector feca30x (norows);
  Rcpp::NumericVector juvgiven20x (norows);
  
  Rcpp::NumericVector sizea30x (norows);
  Rcpp::NumericVector sizeb30x (norows);
  Rcpp::NumericVector sizec30x (norows);
  Rcpp::NumericVector repstra30x (norows);
  Rcpp::NumericVector repstrb30x (norows);
  Rcpp::NumericVector fecb20x (norows);
  Rcpp::NumericVector fecb30x (norows);
  Rcpp::NumericVector juvgiven30x (norows);
  
  Rcpp::StringVector stage2x;
  Rcpp::StringVector stage3x;
  
  // Assign values from the dataset
  if (popidcol != -1) {popidx = as<StringVector>(data[popidcol]);} else {popidx = as<StringVector>(clone(zerovec));}
  if (patchidcol != -1) {patchidx = as<StringVector>(data[patchidcol]);} else {patchidx = as<StringVector>(clone(zerovec));}
  if (censorcol != -1) {censor2x = as<NumericVector>(data[censorcol]);} else {censor2x = clone(zerovec);}
  
  if (sizea2col != -1) {sizea2x = as<NumericVector>(data[sizea2col]);} else {sizea2x = clone(zerovec);}
  if (sizeb2col != -1) {sizeb2x = as<NumericVector>(data[sizeb2col]);} else {sizeb2x = clone(zerovec);}
  if (sizec2col != -1) {sizec2x = as<NumericVector>(data[sizec2col]);} else {sizec2x = clone(zerovec);}
  if (repstra2col != -1) {repstra2x = as<NumericVector>(data[repstra2col]);} else {repstra2x = clone(zerovec);}
  if (repstrb2col != -1) {repstrb2x = as<NumericVector>(data[repstrb2col]);} else {repstrb2x = clone(zerovec);}
  if (feca2col != -1) {feca2x = as<NumericVector>(data[feca2col]);} else {feca2x = clone(zerovec);}
  if (fecb2col != -1) {fecb2x = as<NumericVector>(data[fecb2col]);} else {fecb2x = clone(zerovec);}
  
  if (indcova2col != -1) {indcova2x = as<NumericVector>(data[indcova2col]);} else {indcova2x = clone(zerovec);}
  if (indcova3col != -1) {indcova3x = as<NumericVector>(data[indcova3col]);} else {indcova3x = clone(zerovec);}
  if (indcovb2col != -1) {indcovb2x = as<NumericVector>(data[indcovb2col]);} else {indcovb2x = clone(zerovec);}
  if (indcovb3col != -1) {indcovb3x = as<NumericVector>(data[indcovb3col]);} else {indcovb3x = clone(zerovec);}
  if (indcovc2col != -1) {indcovc2x = as<NumericVector>(data[indcovc2col]);} else {indcovc2x = clone(zerovec);}
  if (indcovc3col != -1) {indcovc3x = as<NumericVector>(data[indcovc3col]);} else {indcovc3x = clone(zerovec);}
  
  if (juv2col != -1) {juvgiven2x = as<NumericVector>(data[juv2col]);} else {juvgiven2x = clone(zerovec);}
  if (obs2col != -1) {obsgiven2x = as<NumericVector>(data[obs2col]);} else {obsgiven2x = clone(negonevec);}
  if (nonobs2col != -1) {nonobsgiven2x = as<NumericVector>(data[nonobs2col]);} else {nonobsgiven2x = clone(negonevec);}
  if (alive2col != -1) {alivegiven2x = as<NumericVector>(data[alive2col]);} else {alivegiven2x = clone(negonevec);}
  if (dead2col != -1) {deadgiven2x = as<NumericVector>(data[dead2col]);} else {deadgiven2x = clone(negonevec);}
  if (xcol != -1) {xpos2x = as<NumericVector>(data[xcol]);} else {xpos2x = clone(negonevec);}
  if (ycol != -1) {ypos2x = as<NumericVector>(data[ycol]);} else {ypos2x = clone(negonevec);}
  
  if (sizea3col != -1) {sizea3x = as<NumericVector>(data[sizea3col]);} else {sizea3x = clone(zerovec);}
  if (sizeb3col != -1) {sizeb3x = as<NumericVector>(data[sizeb3col]);} else {sizeb3x = clone(zerovec);}
  if (sizec3col != -1) {sizec3x = as<NumericVector>(data[sizec3col]);} else {sizec3x = clone(zerovec);}
  if (repstra3col != -1) {repstra3x = as<NumericVector>(data[repstra3col]);} else {repstra3x = clone(zerovec);}
  if (repstrb3col != -1) {repstrb3x = as<NumericVector>(data[repstrb3col]);} else {repstrb3x = clone(zerovec);}
  if (feca3col != -1) {feca3x = as<NumericVector>(data[feca3col]);} else {feca3x = clone(zerovec);}
  if (fecb3col != -1) {fecb3x = as<NumericVector>(data[fecb3col]);} else {fecb3x = clone(zerovec);}
  if (juv3col != -1) {juvgiven3x = as<NumericVector>(data[juv3col]);} else {juvgiven3x = clone(zerovec);}
  if (obs3col != -1) {obsgiven3x = as<NumericVector>(data[obs3col]);} else {obsgiven3x = clone(negonevec);}
  if (nonobs3col != -1) {nonobsgiven3x = as<NumericVector>(data[nonobs3col]);} else {nonobsgiven3x = clone(negonevec);}
  if (alive3col != -1) {alivegiven3x = as<NumericVector>(data[alive3col]);} else {alivegiven3x = clone(negonevec);}
  if (dead3col != -1) {deadgiven3x = as<NumericVector>(data[dead3col]);} else {deadgiven3x = clone(negonevec);}
  
  if (stage2col != -1) {stage2x = as<StringVector>(data[stage2col]);} else {stage2x = as<StringVector>(clone(zerovec));}
  if (stage3col != -1) {stage3x = as<StringVector>(data[stage3col]);} else {stage3x = as<StringVector>(clone(zerovec));}
  
  Rcpp::NumericVector rowid (ndflength);
  Rcpp::StringVector popid (ndflength);
  Rcpp::StringVector patchid (ndflength);
  Rcpp::StringVector individ (ndflength);
  Rcpp::NumericVector censor2 (ndflength);
  Rcpp::NumericVector xpos2 (ndflength);
  Rcpp::NumericVector ypos2 (ndflength);
  Rcpp::NumericVector sizea2 (ndflength);
  Rcpp::NumericVector sizeb2 (ndflength);
  Rcpp::NumericVector sizec2 (ndflength);
  Rcpp::NumericVector sizeadded2 (ndflength);
  Rcpp::NumericVector repstra2 (ndflength);
  Rcpp::NumericVector repstrb2 (ndflength);
  Rcpp::NumericVector repstradded2 (ndflength);
  Rcpp::NumericVector feca2 (ndflength);
  Rcpp::NumericVector fecb2 (ndflength);
  Rcpp::NumericVector fecadded2 (ndflength);
  
  arma::ivec year2 (ndflength);
  year2.zeros();
  
  Rcpp::NumericVector indcova2 (ndflength);
  Rcpp::NumericVector indcovb2 (ndflength);
  Rcpp::NumericVector indcovc2 (ndflength);
  
  Rcpp::IntegerVector repstatus2 (ndflength);
  Rcpp::IntegerVector fecstatus2 (ndflength);
  Rcpp::IntegerVector obsstatus2 (ndflength);
  Rcpp::NumericVector juvgiven2 (ndflength);
  Rcpp::IntegerVector matstat2 (ndflength);
  
  censor2.fill(crazycensor);
  
  xpos2.fill(0.0);
  ypos2.fill(0.0);
  sizea2.fill(0.0);
  sizeb2.fill(0.0);
  sizec2.fill(0.0);
  sizeadded2.fill(0.0);
  repstra2.fill(0.0);
  repstrb2.fill(0.0);
  repstradded2.fill(0.0);
  feca2.fill(0.0);
  fecb2.fill(0.0);
  fecadded2.fill(0.0);
  indcova2.fill(0.0);
  indcovb2.fill(0.0);
  indcovc2.fill(0.0);
  repstatus2.fill(0);
  fecstatus2.fill(0);
  obsstatus2.fill(0);
  juvgiven2.fill(0.0);
  matstat2.fill(0);
  
  Rcpp::NumericVector censor3 (ndflength);
  Rcpp::NumericVector xpos3 (ndflength);
  Rcpp::NumericVector ypos3 (ndflength);
  Rcpp::NumericVector sizea3 (ndflength);
  Rcpp::NumericVector sizeb3 (ndflength);
  Rcpp::NumericVector sizec3 (ndflength);
  Rcpp::NumericVector sizeadded3 (ndflength);
  Rcpp::NumericVector repstra3 (ndflength);
  Rcpp::NumericVector repstrb3 (ndflength);
  Rcpp::NumericVector repstradded3 (ndflength);
  Rcpp::NumericVector feca3 (ndflength);
  Rcpp::NumericVector fecb3 (ndflength);
  Rcpp::NumericVector fecadded3 (ndflength);
  
  Rcpp::NumericVector indcova3 (ndflength);
  Rcpp::NumericVector indcovb3 (ndflength);
  Rcpp::NumericVector indcovc3 (ndflength);
  
  Rcpp::IntegerVector repstatus3 (ndflength);
  Rcpp::IntegerVector fecstatus3 (ndflength);
  Rcpp::IntegerVector obsstatus3 (ndflength);
  Rcpp::NumericVector juvgiven3 (ndflength);
  Rcpp::IntegerVector matstat3 (ndflength);
  
  censor3.fill(crazycensor);
  
  xpos3.fill(0.0);
  ypos3.fill(0.0);
  sizea3.fill(0.0);
  sizeb3.fill(0.0);
  sizec3.fill(0.0);
  sizeadded3.fill(0.0);
  repstra3.fill(0.0);
  repstrb3.fill(0.0);
  repstradded3.fill(0.0);
  feca3.fill(0.0);
  fecb3.fill(0.0);
  fecadded3.fill(0.0);
  indcova3.fill(0.0);
  indcovb3.fill(0.0);
  indcovc3.fill(0.0);
  repstatus3.fill(0);
  fecstatus3.fill(0);
  obsstatus3.fill(0);
  juvgiven3.fill(0.0);
  matstat3.fill(0);
  
  Rcpp::NumericVector censor1 (ndflength);
  Rcpp::NumericVector xpos1 (ndflength);
  Rcpp::NumericVector ypos1 (ndflength);
  Rcpp::NumericVector sizea1 (ndflength);
  Rcpp::NumericVector sizeb1 (ndflength);
  Rcpp::NumericVector sizec1 (ndflength);
  Rcpp::NumericVector sizeadded1 (ndflength);
  Rcpp::NumericVector repstra1 (ndflength);
  Rcpp::NumericVector repstrb1 (ndflength);
  Rcpp::NumericVector repstradded1 (ndflength);
  Rcpp::NumericVector feca1 (ndflength);
  Rcpp::NumericVector fecb1 (ndflength);
  Rcpp::NumericVector fecadded1 (ndflength);
  
  Rcpp::NumericVector indcova1 (ndflength);
  Rcpp::NumericVector indcovb1 (ndflength);
  Rcpp::NumericVector indcovc1 (ndflength);
  
  Rcpp::IntegerVector repstatus1 (ndflength);
  Rcpp::IntegerVector fecstatus1 (ndflength);
  Rcpp::IntegerVector obsstatus1 (ndflength);
  Rcpp::NumericVector juvgiven1 (ndflength);
  Rcpp::IntegerVector matstat1 (ndflength);
  
  censor1.fill(crazycensor);
  
  xpos1.fill(0.0);
  ypos1.fill(0.0);
  sizea1.fill(0.0);
  sizeb1.fill(0.0);
  sizec1.fill(0.0);
  sizeadded1.fill(0.0);
  repstra1.fill(0.0);
  repstrb1.fill(0.0);
  repstradded1.fill(0.0);
  feca1.fill(0.0);
  fecb1.fill(0.0);
  fecadded1.fill(0.0);
  indcova1.fill(0.0);
  indcovb1.fill(0.0);
  indcovc1.fill(0.0);
  repstatus1.fill(0);
  fecstatus1.fill(0);
  obsstatus1.fill(0);
  juvgiven1.fill(0.0);
  matstat1.fill(0);
  
  // This section introduces variables used to check whether the censor variable
  // has been checked and set at each step
  arma::uvec indivnum (ndflength);
  arma::uvec censor2check (ndflength);
  indivnum.zeros();
  censor2check.zeros();
  
  // Here we introduce some derived variables that require extra looping or
  // other control parameters
  arma::ivec firstseen (ndflength);
  arma::ivec lastseen (ndflength);
  arma::ivec obsage (ndflength);
  arma::ivec obslifespan (ndflength);
  Rcpp::NumericVector alive1 (ndflength);
  Rcpp::NumericVector alive2 (ndflength);
  Rcpp::NumericVector alive3 (ndflength);
  firstseen.fill(-1);
  lastseen.fill(-1);
  alive1.fill(0.0);
  alive2.fill(0.0);
  alive3.fill(0.0);
  
  Rcpp::StringVector stage1 (ndflength);
  Rcpp::StringVector stage2 (ndflength);
  Rcpp::StringVector stage3 (ndflength);
  Rcpp::NumericVector stage1num (ndflength);
  Rcpp::NumericVector stage2num (ndflength);
  Rcpp::NumericVector stage3num (ndflength);
  stage1num.fill(0.0);
  stage2num.fill(0.0);
  stage3num.fill(0.0);
  
  double stagesize1 {0.0};
  double stagesize2 {0.0};
  double stagesize3 {0.0};
  double stagesize1b {0.0};
  double stagesize2b {0.0};
  double stagesize3b {0.0};
  double stagesize1c {0.0};
  double stagesize2c {0.0};
  double stagesize3c {0.0};
  
  arma::uvec stagemini1;
  arma::uvec stagemaxi1;
  arma::uvec stagemini2;
  arma::uvec stagemaxi2;
  arma::uvec stagemini3;
  arma::uvec stagemaxi3;
  arma::uvec stageobs;
  arma::uvec stagerep;
  arma::uvec stagemat;
  arma::uvec cs1;
  arma::uvec cs2;
  arma::uvec cs3;
  arma::uvec cs4;
  int choicestage {0};
  
  // Main loop, which creates the main new dataset rows. Establishes state in
  // time t for all cases in which an individual is observed
  for (int i = 0; i < norows; i++) { // i corresponds to row in the old dataset
    if (i % 5 == 0) Rcpp::checkUserInterrupt();
    
    for (int j = 0; j < noyears; j++) {
      // This establishes a place marker for vectors corresponding to the current year
      if (year2x[i] == yearall2x[j]) currentyear = j;
    }
    
    // This establishes a place marker variable corresponding to the current individual
    currentindiv = -1;
    for (int k = 0; k < noindivs; k++) {
      if (individx[i] == allindivs[k]) {
        currentindiv = k;
        indivnum[i] = k;
      }
    }
    
    // This establishes the row in the new dataset being created
    ndfindex = (noyears * currentindiv) + currentyear;
    
    if (NumericVector::is_na(sizea2x[i])) {
      sizea20x[i] = 0.0;
      if (NAas0) {sizea2x[i] = 0.0;}
    } else {sizea20x[i] = sizea2x[i];}
    if (NumericVector::is_na(sizeb2x[i])) {
      sizeb20x[i] = 0.0;
      if (NAas0) {sizeb2x[i] = 0.0;}
    } else {sizeb20x[i] = sizeb2x[i];}
    if (NumericVector::is_na(sizec2x[i])) {
      sizec20x[i] = 0.0;
      if (NAas0) {sizec2x[i] = 0.0;}
    } else {sizec20x[i] = sizec2x[i];}
    if (NumericVector::is_na(repstra2x[i])) {
      repstra20x[i] = 0.0;
      if (NAas0) {repstra2x[i] = 0.0;}
    } else {repstra20x[i] = repstra2x[i];}
    if (NumericVector::is_na(repstrb2x[i])) {
      repstrb20x[i] = 0.0;
      if (NAas0) {repstrb2x[i] = 0.0;}
    } else {repstrb20x[i] = repstrb2x[i];}
    if (NumericVector::is_na(feca2x[i])) {
      feca20x[i] = 0.0;
      if (NAas0) {feca2x[i] = 0.0;}
    } else {feca20x[i] = feca2x[i];}
    if (NumericVector::is_na(fecb2x[i])) {
      fecb20x[i] = 0.0;
      if (NAas0) {fecb2x[i] = 0.0;}
    } else {fecb20x[i] = fecb2x[i];}
    
    if (NumericVector::is_na(juvgiven2x[i])) {
      juvgiven20x[i] = 0.0;
    } else if (juvgiven2x[i] != 0) {
      juvgiven20x[i] = 1.0;
    } else {juvgiven20x[i] = 0.0;}
    
    // Here we develop the censoring variable
    if (censbool && censorcol != -1) {
      // This section provides replacements in cases where NA designates data to keep
      if (NumericVector::is_na(censor2x[i])) {
        censor2[ndfindex] = 0.0;
        censor2check[ndfindex] = 1;
      } else {
        censor2[ndfindex] = 1.0;
        censor2check[ndfindex] = 1;
      }
    } else if (censorcol != -1) {
      if (censorkeep == 0 && NumericVector::is_na(censor2x[i])) {
        censor2[ndfindex] = 1.0;
        censor2check[ndfindex] = 1;
      } else if (censorkeep == 1 && NumericVector::is_na(censor2x[i])) {
        censor2[ndfindex] = 0.0;
        censor2check[ndfindex] = 1;
      } else {
        censor2[ndfindex] = censor2x[i];
        censor2check[ndfindex] = 1;
      }
    }
    
    rowid[ndfindex] = static_cast<double>(i);
    popid[ndfindex] = popidx[i];
    patchid[ndfindex] = patchidx[i];
    individ[ndfindex] = allindivs[currentindiv];
    year2[ndfindex] = yearall2x[currentyear];
    xpos2[ndfindex] = xpos2x[i];
    ypos2[ndfindex] = ypos2x[i];
    sizea2[ndfindex] = sizea2x[i];
    sizeb2[ndfindex] = sizeb2x[i];
    sizec2[ndfindex] = sizec2x[i];
    sizeadded2[ndfindex] = sizea20x[i] + sizeb20x[i] + sizec20x[i];
    
    repstra2[ndfindex] = repstra2x[i];
    repstrb2[ndfindex] = repstrb2x[i];
    repstradded2[ndfindex] = repstra20x[i] + (repstrb20x[i] * repstrrel);
    
    feca2[ndfindex] = feca2x[i];
    fecb2[ndfindex] = fecb2x[i];
    fecadded2[ndfindex] = feca20x[i] + (fecb20x[i] * fecrel);
    
    indcova2[ndfindex] = indcova2x[i];
    indcovb2[ndfindex] = indcovb2x[i];
    indcovc2[ndfindex] = indcovc2x[i];
    
    if (repstradded2[ndfindex] > 0) {
      repstatus2[ndfindex] = 1;} else {repstatus2[ndfindex] = 0;
    }
    if (NumericVector::is_na(obsgiven2x[i])) {
      obsstatus2[ndfindex] = 0;
    } else if (obsgiven2x[i] > 0) {
      obsstatus2[ndfindex] = 1;
    } else if (obsgiven2x[i] == -1 && (sizeadded2[ndfindex] + repstradded2[ndfindex]) > 0) {
      obsstatus2[ndfindex] = 1;
    } else {obsstatus2[ndfindex] = 0;}
    
    if (nonobsgiven2x[i] >= 0) {
      obsstatus2[ndfindex] = 0;
    } else if (nonobsgiven2x[i] == 0) {
      obsstatus2[ndfindex] = 1;
    }
    
    juvgiven2[ndfindex] = juvgiven20x[i];
    matstat2[ndfindex] = 1 - juvgiven2[ndfindex];
    
    if (alivegiven2x[i] > 0) {
      alive2[ndfindex] = 1.0;
    } else if (alivegiven2x[i] == 0) {
      alive2[ndfindex] = 0.0;
    }
    
    if (deadgiven2x[i] > 0) {
      alive2[ndfindex] = 0.0;
    } else if (deadgiven2x[i] == 0) {
      alive2[ndfindex] = 1.0;
    }
    
    livcheck2 = sizeadded2[ndfindex] + repstradded2[ndfindex] + obsstatus2[ndfindex];
    
    if (livcheck2 > 0 && firstseenx[currentindiv] == -1) {
      firstseenx[currentindiv] = currentyear + firstyear;
      lastseenx[currentindiv] = currentyear + firstyear;
      alive2[ndfindex] = 1.0;
    } else if (livcheck2 > 0) {
      lastseenx[currentindiv] = currentyear + firstyear;
      alive2[ndfindex] = 1.0;
    }
    
    if (alive2[ndfindex] == 1 && matstat2[ndfindex] == 1) {
      matstat3[ndfindex] = 1;
    }
    
    if (stage2col != -1 && alive2[ndfindex] == 1) {
      stage2[ndfindex] = stage2x[i];
    } else if (stassign) {stage2[ndfindex] = "NotAlive";}
    
    
    // Now we work on time t+1 for the last possible time t (which is technically
    // the second to last time), in cases where t+1 columns are provided
    if (currentyear == (noyears - 1)) {
      if (censbool && censorcol != -1) { // Here we develop the censoring variable for the last time
        if (NumericVector::is_na(censor2x[i])) {
          censor3[ndfindex] = 0.0;
        } else {
          censor3[ndfindex] = 1.0;
        }
      } else if (censorcol != -1) {     // New section through 1908
        if (censorkeep == 0 && NumericVector::is_na(censor2x[i])) {
          censor3[ndfindex] = 1.0;
        } else if (censorkeep == 1 && NumericVector::is_na(censor2x[i])) {
          censor3[ndfindex] = 0.0;
        } else {
          censor3[ndfindex] = censor2x[i];
        }
      }
      
      if (NumericVector::is_na(juvgiven3x[i])) {
        juvgiven30x[i] = 0.0;
      } else if (juvgiven3x[i] != 0) {
        juvgiven30x[i] = 1.0;
      } else {juvgiven30x[i] = 0.0;}
      
      if (sizea3col != -1) {
        if (NumericVector::is_na(sizea3x[i])) {
          sizea30x[i] = 0.0;
          if (NAas0) {sizea3x[i] = 0.0;}
        } else {sizea30x[i] = sizea3x[i];}
        sizea3[ndfindex] = sizea3x[i];
      }
      if (sizeb3col != -1) {
        if (NumericVector::is_na(sizeb3x[i])) {
          sizeb30x[i] = 0.0;
          if (NAas0) {sizeb3x[i] = 0.0;}
        } else {sizeb30x[i] = sizeb3x[i];}
        sizeb3[ndfindex] = sizeb3x[i];
      }
      if (sizec3col != -1) {
        if (NumericVector::is_na(sizec3x[i])) {
          sizec30x[i] = 0.0;
          if (NAas0) {sizec3x[i] = 0.0;}
        } else {sizec30x[i] = sizec3x[i];}
        sizec3[ndfindex] = sizec3x[i];
      }
      sizeadded3[ndfindex] = sizea30x[i] + sizeb30x[i] + sizec30x[i];
      
      if (repstra3col != -1) {
        if (NumericVector::is_na(repstra3x[i])) {
          repstra30x[i] = 0.0;
          if (NAas0) {repstra3x[i] = 0.0;}
        } else {repstra30x[i] = repstra3x[i];}
        repstra3[ndfindex] = repstra3x[i];
      }
      if (repstrb3col != -1) {
        if (NumericVector::is_na(repstrb3x[i])) {
          repstrb30x[i] = 0.0;
          if (NAas0) {repstrb3x[i] = 0.0;}
        } else {repstrb30x[i] = repstrb3x[i];}
        repstrb3[ndfindex] = repstrb3x[i];
      }
      repstradded3[ndfindex] = repstra30x[i] + (repstrb30x[i] * repstrrel);
      
      if (feca3col != -1) {
        if (NumericVector::is_na(feca3x[i])) {
          feca30x[i] = 0.0;
          if (NAas0) {feca3x[i] = 0.0;}
        } else {feca30x[i] = feca3x[i];}
        feca3[ndfindex] = feca3x[i];
      }
      if (fecb3col != -1) {
        if (NumericVector::is_na(fecb3x[i])) {
          fecb30x[i] = 0.0;
          if (NAas0) {feca3x[i] = 0.0;}
        } else {fecb30x[i] = fecb3x[i];}
        fecb3[ndfindex] = fecb3x[i];
      }
      fecadded3[ndfindex] = feca30x[i] + (fecb30x[i] * fecrel);
      if (fecadded3[ndfindex] > 0) {fecstatus3[ndfindex] = 1;}
      
      if (repstradded3[ndfindex] > 0) {repstatus3[ndfindex] = 1;} else {repstatus3[ndfindex] = 0;}
      
      if (indcova3col != -1) {
        indcova3[ndfindex] = indcova3x[i];
      }
      if (indcovb3col != -1) {
        indcovb3[ndfindex] = indcovb3x[i];
      }
      if (indcovc3col != -1) {
        indcovc3[ndfindex] = indcovc3x[i];
      }
      
      if (NumericVector::is_na(obsgiven3x[i])) {
        obsstatus3[ndfindex] = 0;
      } else if (obsgiven3x[i] > 0) {
        obsstatus3[ndfindex] = 1;
      } else if (obsgiven3x[i] == -1 && (sizeadded3[ndfindex] + repstradded3[ndfindex]) > 0) {
        obsstatus3[ndfindex] = 1;
      } else {obsstatus3[ndfindex] = 0;}
      
      if (nonobsgiven3x[i] >= 0) {
        obsstatus3[ndfindex] = 0;
      } else if (nonobsgiven3x[i] == 0) {
        obsstatus3[ndfindex] = 1;
      }
      
      juvgiven3[ndfindex] = juvgiven30x[i];
      matstat3[ndfindex] = 1 - juvgiven3[ndfindex];
      
      if (alivegiven3x[i] > 0) {
        alive3[ndfindex] = 1.0;
      } else if (alivegiven3x[i] == 0) {
        alive3[ndfindex] = 0.0;
      }
      if (deadgiven3x[i] > 0) {
        alive3[ndfindex] = 0.0;
      } else if (deadgiven3x[i] == 0) {
        alive3[ndfindex] = 1.0;
      }
      
      livcheck3 = sizeadded3[ndfindex] + repstradded3[ndfindex] + obsstatus3[ndfindex];
      
      if (firstseenx[currentindiv] == -1 && livcheck3 > 0) {
        firstseenx[currentindiv] = currentyear + firstyear + 1;
        lastseenx[currentindiv] = currentyear + firstyear + 1;
        alive3[ndfindex] = 1;
        if (juv3col == -1 && repstradded2[ndfindex] > 0) {
          matstat3[ndfindex] = 1;
        }
      } else if (livcheck3 > 0) {
        lastseenx[currentindiv] = currentyear + firstyear + 1;
        alive3[ndfindex] = 1;
        if (juv3col == -1 && repstradded2[ndfindex] > 0) {
          matstat3[ndfindex] = 1;
        }
      }
      
      if (stage3col != -1 && alive3[ndfindex] == 1) {
        stage3[ndfindex] = stage2x[i];
      } else if (stassign) {stage3[ndfindex] = "NotAlive";}
    } // End of currentyear if statement
  } // End of i loop
  
  
  // Now a loop that establishes most states in time t+1 and t-1,
  // and stages in all times
  for (int i = 0; i < ndflength; i++) { // i refers to rows in the final dataset
    if (i % 5 == 0) Rcpp::checkUserInterrupt();
    
    // This short section deals with correcting info for individuals that are
    // unobserved for long periods
    if (i > 0 && rowid[i] == 0) {
      if (year2[i-1] < lastseen[i-1] && (year2[i-1] + 1) < (firstyear + noyears)) {
        individ[i] = individ[i-1];
        popid[i] = popid[i-1];
        patchid[i] = patchid[i-1];
        year2[i] = year2[i-1] + 1;
        xpos2[i] = xpos2[i-1]; // This is new
        ypos2[i] = ypos2[i-1]; // This is new
        
        if (matstat2[i-1] == 1) {
          matstat2[i] = 1;
          
          if (year2[i+1] <= lastseen[i]) {
            matstat3[i] = 1;
          }
        }
      }
    }
    
    // Now the normal calculations
    currentindiv = -1;
    for (int k = 0; k < noindivs; k++) {
      if (individ[i] == allindivs[k]) currentindiv = k;
    }
    
    if (currentindiv != -1) { // This makes sure that we are only dealing with real individuals in the dataset
      if (year2[i] <= lastseenx[currentindiv] && year2[i] >= firstseenx[currentindiv] && year2[i] < (firstyear + noyears)) {
        firstseen[i] = firstseenx[currentindiv];
        lastseen[i] = lastseenx[currentindiv];
        
        obsage[i] = year2[i] - firstseen[i];
        obslifespan[i] = lastseen[i] - firstseen[i];
      }
      
      if (year2[i] >= firstseen[i] && year2[i] <= lastseen[i] && alive2[i] == 0) {
        alive2[i] = 1.0;
      } else if (alive2[i] == 0) {
        alive2[i] = 0.0;
      }
      
      if ((year2[i] + 1) >= firstseen[i] && (year2[i] + 1) <= lastseen[i] && alive3[i] == 0) {
        alive3[i] = 1.0;
      } else if (alive3[i] == 0) {
        alive3[i] = 0.0;
      }
      
      currentyear = year2[i] - firstyear;
      
      if (fecadded2[i] > 0) {
        fecstatus2[i] = 1;
      } else {
        fecstatus2[i] = 0;
      }
      
      if (stage2col != -1 && alive2[i] == 1) {
        if (obsstatus2[i] == 0) {
          stage2[i] = "NotObserved";
        }
      } else if (stage2col != -1 && alive2[i] == 0 && stassign) {
        stage2[i] = "NotAlive";
      }
      
      if (currentyear < (noyears - 1)) { 
        nextyrindex = (noyears * currentindiv) + (currentyear + 1);
        
        if (censor2[i] == censorkeep && alive2[i] == 1) {
          censor3[i] = censorkeep;
        } else {
          censor3[i] = censor2[nextyrindex];
        }
        
        xpos3[i] = xpos2[nextyrindex];
        ypos3[i] = ypos2[nextyrindex];
        
        sizea3[i] = sizea2[nextyrindex];
        sizeb3[i] = sizeb2[nextyrindex];
        sizec3[i] = sizec2[nextyrindex];
        sizeadded3[i] = sizeadded2[nextyrindex];
        
        repstra3[i] = repstra2[nextyrindex];
        repstrb3[i] = repstrb2[nextyrindex];
        repstradded3[i] = repstradded2[nextyrindex];
        
        feca3[i] = feca2[nextyrindex];
        fecb3[i] = fecb2[nextyrindex];
        fecadded3[i] = fecadded2[nextyrindex];
        
        indcova3[i] = indcova2[nextyrindex];
        indcovb3[i] = indcovb2[nextyrindex];
        indcovc3[i] = indcovc2[nextyrindex];
        
        if (fecadded3[i] > 0) {
          fecstatus3[i] = 1;
        } else {
          fecstatus3[i] = 0;
        }
        
        repstatus3[i] = repstatus2[nextyrindex];
        obsstatus3[i] = obsstatus2[nextyrindex];
        juvgiven3[i] = juvgiven2[nextyrindex];
        matstat3[i] = matstat2[nextyrindex];
        
        if (matstat2[i] == 1 && stage3[i] != "NotAlive") {
          matstat3[i] = 1;
        }
        
        if (stage2col != -1 && alive3[i] == 1 && stassign) {
          if (obsstatus3[i] == 0) {
            stage3[i] = "NotObserved";
          } else {stage3[i] = stage2[nextyrindex];}
        } else if (stage2col != -1 && alive3[i] == 0 && stassign) {
          stage3[i] = "NotAlive";
        }
      }
      
      if (currentyear > 0  && year2[i] < (firstyear + noyears)) {
        prevyrindex = (noyears * currentindiv) + (currentyear - 1);
        
        alive1[i] = alive2[prevyrindex];
        
        if (censor2(i) == censorkeep && alive1(i) == 1) {
          if (censbool) {
            censor1(i) = 0;
          } else {
            censor1[i] = censorkeep;
          }
        } else {
          censor1[i] = censor2[prevyrindex];
        }
        
        xpos1[i] = xpos2[prevyrindex];
        ypos1[i] = ypos2[prevyrindex];
        
        sizea1[i] = sizea2[prevyrindex];
        sizeb1[i] = sizeb2[prevyrindex];
        sizec1[i] = sizec2[prevyrindex];
        sizeadded1[i] = sizeadded2[prevyrindex];
        
        repstra1[i] = repstra2[prevyrindex];
        repstrb1[i] = repstrb2[prevyrindex];
        repstradded1[i] = repstradded2[prevyrindex];
        
        feca1[i] = feca2[prevyrindex];
        fecb1[i] = fecb2[prevyrindex];
        fecadded1[i] = fecadded2[prevyrindex];
        
        indcova1[i] = indcova2[prevyrindex];
        indcovb1[i] = indcovb2[prevyrindex];
        indcovc1[i] = indcovc2[prevyrindex];
        
        if (fecadded1[i] > 0) {
          fecstatus1[i] = 1;
        } else {
          fecstatus1[i] = 0;
        }
        
        repstatus1[i] = repstatus2[prevyrindex];
        obsstatus1[i] = obsstatus2[prevyrindex];
        juvgiven1[i] = juvgiven2[prevyrindex];
        matstat1[i] = matstat2[prevyrindex];
        
        if (stage2col != -1 && alive1[i] == 1 && stassign) {
          if (obsstatus1[i] == 0) {
            stage1[i] = "NotObserved";
          } else {
            stage1[i] = stage2[prevyrindex];
          }
        } else if (stage2col != -1 && alive1[i] == 0 && stassign) {
          stage1[i] = "NotAlive";
        }
      }
      
      // Coordinaete corrections
      if (xpos1[i] != xpos2[i] && xpos1[i] == 0.0) {
        xpos1[i] = xpos2[i];
      }
      if (xpos3[i] != xpos2[i] && xpos3[i] == 0.0) {
        xpos3[i] = xpos2[i];
      }
      if (ypos1[i] != ypos2[i] && ypos1[i] == 0.0) {
        ypos1[i] = ypos2[i];
      }
      if (ypos3[i] != ypos2[i] && ypos3[i] == 0.0) {
        ypos3[i] = ypos2[i];
      }
      
      // Stage assignments
      if (stassign && stage2col == -1) {
        if (stszcol == 8) {
          stagesize1 = sizea1[i];
          stagesize2 = sizea2[i];
          stagesize3 = sizea3[i];
          stagesize1b = sizeb1[i];
          stagesize2b = sizeb2[i];
          stagesize3b = sizeb3[i];
          stagesize1c = sizec1[i];
          stagesize2c = sizec2[i];
          stagesize3c = sizec3[i];
          
          arma::uvec stagemini1a = find(insfszminarma < stagesize1);
          arma::uvec stagemaxi1a = find(insfszmaxarma >= stagesize1);
          arma::uvec stagemini1b = find(insfszminarmab < stagesize1b);
          arma::uvec stagemaxi1b = find(insfszmaxarmab >= stagesize1b);
          arma::uvec stagemini1c = find(insfszminarmac < stagesize1c);
          arma::uvec stagemaxi1c = find(insfszmaxarmac >= stagesize1c);
          arma::uvec stagemini2a = find(insfszminarma < stagesize2);
          arma::uvec stagemaxi2a = find(insfszmaxarma >= stagesize2);
          arma::uvec stagemini2b = find(insfszminarmab < stagesize2b);
          arma::uvec stagemaxi2b = find(insfszmaxarmab >= stagesize2b);
          arma::uvec stagemini2c = find(insfszminarmac < stagesize2c);
          arma::uvec stagemaxi2c = find(insfszmaxarmac >= stagesize2c);
          arma::uvec stagemini3a = find(insfszminarma < stagesize3);
          arma::uvec stagemaxi3a = find(insfszmaxarma >= stagesize3);
          arma::uvec stagemini3b = find(insfszminarmab < stagesize3b);
          arma::uvec stagemaxi3b = find(insfszmaxarmab >= stagesize3b);
          arma::uvec stagemini3c = find(insfszminarmac < stagesize3c);
          arma::uvec stagemaxi3c = find(insfszmaxarmac >= stagesize3c);
          
          arma::uvec stagemini1d = intersect(stagemini1a, stagemini1b);
          stagemini1 = intersect(stagemini1c, stagemini1d);
          arma::uvec stagemaxi1d = intersect(stagemaxi1a, stagemaxi1b);
          stagemaxi1 = intersect(stagemaxi1c, stagemaxi1d);
          arma::uvec stagemini2d = intersect(stagemini2a, stagemini2b);
          stagemini2 = intersect(stagemini2c, stagemini2d);
          arma::uvec stagemaxi2d = intersect(stagemaxi2a, stagemaxi2b);
          stagemaxi2 = intersect(stagemaxi2c, stagemaxi2d);
          arma::uvec stagemini3d = intersect(stagemini3a, stagemini3b);
          stagemini3 = intersect(stagemini3c, stagemini3d);
          arma::uvec stagemaxi3d = intersect(stagemaxi3a, stagemaxi3b);
          stagemaxi3 = intersect(stagemaxi3c, stagemaxi3d);
          
        } else if (stszcol == 7) {
          stagesize1 = sizeb1[i];
          stagesize2 = sizeb2[i];
          stagesize3 = sizeb3[i];
          stagesize1b = sizec1[i];
          stagesize2b = sizec2[i];
          stagesize3b = sizec3[i];
          
          arma::uvec stagemini1a = find(insfszminarma < stagesize1);
          arma::uvec stagemaxi1a = find(insfszmaxarma >= stagesize1);
          arma::uvec stagemini1b = find(insfszminarmab < stagesize1b);
          arma::uvec stagemaxi1b = find(insfszmaxarmab >= stagesize1b);
          arma::uvec stagemini2a = find(insfszminarma < stagesize2);
          arma::uvec stagemaxi2a = find(insfszmaxarma >= stagesize2);
          arma::uvec stagemini2b = find(insfszminarmab < stagesize2b);
          arma::uvec stagemaxi2b = find(insfszmaxarmab >= stagesize2b);
          arma::uvec stagemini3a = find(insfszminarma < stagesize3);
          arma::uvec stagemaxi3a = find(insfszmaxarma >= stagesize3);
          arma::uvec stagemini3b = find(insfszminarmab < stagesize3b);
          arma::uvec stagemaxi3b = find(insfszmaxarmab >= stagesize3b);
          
          stagemini1 = intersect(stagemini1a, stagemini1b);
          stagemaxi1 = intersect(stagemaxi1a, stagemaxi1b);
          stagemini2 = intersect(stagemini2a, stagemini2b);
          stagemaxi2 = intersect(stagemaxi2a, stagemaxi2b);
          stagemini3 = intersect(stagemini3a, stagemini3b);
          stagemaxi3 = intersect(stagemaxi3a, stagemaxi3b);
          
        } else if (stszcol == 6) {
          stagesize1 = sizea1[i];
          stagesize2 = sizea2[i];
          stagesize3 = sizea3[i];
          stagesize1b = sizec1[i];
          stagesize2b = sizec2[i];
          stagesize3b = sizec3[i];
          
          arma::uvec stagemini1a = find(insfszminarma < stagesize1);
          arma::uvec stagemaxi1a = find(insfszmaxarma >= stagesize1);
          arma::uvec stagemini1b = find(insfszminarmab < stagesize1b);
          arma::uvec stagemaxi1b = find(insfszmaxarmab >= stagesize1b);
          arma::uvec stagemini2a = find(insfszminarma < stagesize2);
          arma::uvec stagemaxi2a = find(insfszmaxarma >= stagesize2);
          arma::uvec stagemini2b = find(insfszminarmab < stagesize2b);
          arma::uvec stagemaxi2b = find(insfszmaxarmab >= stagesize2b);
          arma::uvec stagemini3a = find(insfszminarma < stagesize3);
          arma::uvec stagemaxi3a = find(insfszmaxarma >= stagesize3);
          arma::uvec stagemini3b = find(insfszminarmab < stagesize3b);
          arma::uvec stagemaxi3b = find(insfszmaxarmab >= stagesize3b);
          
          stagemini1 = intersect(stagemini1a, stagemini1b);
          stagemaxi1 = intersect(stagemaxi1a, stagemaxi1b);
          stagemini2 = intersect(stagemini2a, stagemini2b);
          stagemaxi2 = intersect(stagemaxi2a, stagemaxi2b);
          stagemini3 = intersect(stagemini3a, stagemini3b);
          stagemaxi3 = intersect(stagemaxi3a, stagemaxi3b);
          
        } else if (stszcol == 5) {
          stagesize1 = sizea1[i];
          stagesize2 = sizea2[i];
          stagesize3 = sizea3[i];
          stagesize1b = sizeb1[i];
          stagesize2b = sizeb2[i];
          stagesize3b = sizeb3[i];
          
          arma::uvec stagemini1a = find(insfszminarma < stagesize1);
          arma::uvec stagemaxi1a = find(insfszmaxarma >= stagesize1);
          arma::uvec stagemini1b = find(insfszminarmab < stagesize1b);
          arma::uvec stagemaxi1b = find(insfszmaxarmab >= stagesize1b);
          arma::uvec stagemini2a = find(insfszminarma < stagesize2);
          arma::uvec stagemaxi2a = find(insfszmaxarma >= stagesize2);
          arma::uvec stagemini2b = find(insfszminarmab < stagesize2b);
          arma::uvec stagemaxi2b = find(insfszmaxarmab >= stagesize2b);
          arma::uvec stagemini3a = find(insfszminarma < stagesize3);
          arma::uvec stagemaxi3a = find(insfszmaxarma >= stagesize3);
          arma::uvec stagemini3b = find(insfszminarmab < stagesize3b);
          arma::uvec stagemaxi3b = find(insfszmaxarmab >= stagesize3b);
          
          stagemini1 = intersect(stagemini1a, stagemini1b);
          stagemaxi1 = intersect(stagemaxi1a, stagemaxi1b);
          stagemini2 = intersect(stagemini2a, stagemini2b);
          stagemaxi2 = intersect(stagemaxi2a, stagemaxi2b);
          stagemini3 = intersect(stagemini3a, stagemini3b);
          stagemaxi3 = intersect(stagemaxi3a, stagemaxi3b);
          
        } else if (stszcol == 4) {
          stagesize1 = sizeadded1[i];
          stagesize2 = sizeadded2[i];
          stagesize3 = sizeadded3[i];
          
          stagemini1 = find(insfszminarma < stagesize1);
          stagemaxi1 = find(insfszmaxarma >= stagesize1);
          stagemini2 = find(insfszminarma < stagesize2);
          stagemaxi2 = find(insfszmaxarma >= stagesize2);
          stagemini3 = find(insfszminarma < stagesize3);
          stagemaxi3 = find(insfszmaxarma >= stagesize3);
          
        } else if (stszcol == 3) {
          stagesize1 = sizec1[i];
          stagesize2 = sizec2[i];
          stagesize3 = sizec3[i];
          
          stagemini1 = find(insfszminarma < stagesize1);
          stagemaxi1 = find(insfszmaxarma >= stagesize1);
          stagemini2 = find(insfszminarma < stagesize2);
          stagemaxi2 = find(insfszmaxarma >= stagesize2);
          stagemini3 = find(insfszminarma < stagesize3);
          stagemaxi3 = find(insfszmaxarma >= stagesize3);
        } else if (stszcol == 2) {
          stagesize1 = sizeb1[i];
          stagesize2 = sizeb2[i];
          stagesize3 = sizeb3[i];
          
          stagemini1 = find(insfszminarma < stagesize1);
          stagemaxi1 = find(insfszmaxarma >= stagesize1);
          stagemini2 = find(insfszminarma < stagesize2);
          stagemaxi2 = find(insfszmaxarma >= stagesize2);
          stagemini3 = find(insfszminarma < stagesize3);
          stagemaxi3 = find(insfszmaxarma >= stagesize3);
          
        } else {
          stagesize1 = sizea1[i];
          stagesize2 = sizea2[i];
          stagesize3 = sizea3[i];
          
          stagemini1 = find(insfszminarma < stagesize1);
          stagemaxi1 = find(insfszmaxarma >= stagesize1);
          stagemini2 = find(insfszminarma < stagesize2);
          stagemaxi2 = find(insfszmaxarma >= stagesize2);
          stagemini3 = find(insfszminarma < stagesize3);
          stagemaxi3 = find(insfszmaxarma >= stagesize3);
        }
        
        // Stage 2
        stagerep = find(inrepstatarma == repstatus2[i]);
        stagemat = find(inmatstatarma == matstat2[i]);
        stageobs = find(inobsstatarma == obsstatus2[i]);
        
        cs1 = intersect(stagemini2, stagemaxi2);
        
        if (NRasRep && NOasObs) {
          cs4 = cs1;
        } else if (NOasObs) {
          cs4 = intersect(stagerep, cs1);
        } else if (NRasRep) {
          cs2 = intersect(stageobs, stagemat);
          cs4 = intersect(cs1, cs2);
        } else {
          cs2 = intersect(stageobs, stagemat);
          cs3 = intersect(stagerep, cs2);
          cs4 = intersect(cs1, cs3);
        }
        
        if (cs4.n_elem == 1 && alive2[i] == 1) {
          choicestage = instageid(cs4[0]) - 1;
          stage2num[i] = choicestage + 1;
          
          stage2[i] = sfname[choicestage];
        } else if (cs4.n_elem == 0 && alive2[i] == 1) {
          stage2[i] = "NoMatch";
          
          if (!quiet) Rf_warningcall(R_NilValue,
            "Some stages occurring in the dataset do not match any characteristics in the input stageframe.");
        } else if (alive2[i] != 1)  {
          stage2[i] = "NotAlive";
        } else if (cs4.n_elem > 1) {
          if (!quiet) Rf_warningcall(R_NilValue,
            "Some stages in the input stageframe appear to have the same description. Please make sure that all stages included in the stageframe are defined with unique sets of characteristics.");
          
        } else {
          throw Rcpp::exception("Stage assignment error.", false);
        }
        
        // Stage 1
        stagerep = find(inrepstatarma == repstatus1[i]);
        stagemat = find(inmatstatarma == matstat1[i]);
        stageobs = find(inobsstatarma == obsstatus1[i]);
        
        cs1 = intersect(stagemini1, stagemaxi1);
        
        if (NRasRep && NOasObs) {
          cs4 = cs1;
        } else if (NOasObs) {
          cs4 = intersect(stagerep, cs1);
        } else if (NRasRep) {
          cs2 = intersect(stageobs, stagemat);
          cs4 = intersect(cs1, cs2);
        } else {
          cs2 = intersect(stageobs, stagemat);
          cs3 = intersect(stagerep, cs2);
          cs4 = intersect(cs1, cs3);
        }
        
        if (cs4.n_elem == 1 && alive1[i] == 1) {
          choicestage = instageid(cs4[0]) - 1;
          stage1num[i] = choicestage + 1;
          
          stage1[i] = sfname[choicestage];
        } else if (cs4.n_elem == 0 && alive1[i] == 1) {
          stage1[i] = "NoMatch";
          
          if (!quiet) Rf_warningcall(R_NilValue,
            "Some stages occurring in the dataset do not match any characteristics in the input stageframe.");
        } else if (alive1[i] != 1) {
          stage1[i] = "NotAlive";
          matstat1[i] = 0;
        } else if (cs4.n_elem > 1) {
          if (!quiet) Rf_warningcall(R_NilValue,
            "Some stages in the input stageframe appear to have the same description. Please make sure that all stages included in the stageframe are defined with unique sets of characteristics.");
          
        } else {
          throw Rcpp::exception("Stage assignment error.", false);
        }
        
        // Stage 3
        stagerep = find(inrepstatarma == repstatus3[i]);
        stagemat = find(inmatstatarma == matstat3[i]);
        stageobs = find(inobsstatarma == obsstatus3[i]);
        
        cs1 = intersect(stagemini3, stagemaxi3);
        
        if (NRasRep && NOasObs) {
          cs4 = cs1;
        } else if (NOasObs) {
          cs4 = intersect(stagerep, cs1);
        } else if (NRasRep) {
          cs2 = intersect(stageobs, stagemat);
          cs4 = intersect(cs1, cs2);
        } else {
          cs2 = intersect(stageobs, stagemat);
          cs3 = intersect(stagerep, cs2);
          cs4 = intersect(cs1, cs3);
        }
        
        // Here we create exceptions based on stage assignment problems in time t+1
        if (cs4.n_elem == 1 && alive3[i] == 1) {
          choicestage = instageid(cs4[0]) - 1;
          stage3num[i] = choicestage + 1;
          
          stage3[i] = sfname[choicestage];
        } else if (alive3[i] != 1) {
          stage3[i] = "NotAlive";
        } else if (cs4.n_elem == 0) {
          stage3[i] = "NoMatch";
          
          if (!quiet) Rf_warningcall(R_NilValue,
            "Some stages occurring in the dataset do not match any characteristics in the input stageframe.");
        } else if (cs4.n_elem > 1) {
          if (!quiet) Rf_warningcall(R_NilValue,
            "Some stages in the input stageframe appear to have the same description. Please make sure that all stages included in the stageframe are defined with unique sets of characteristics.");
          
        } else {
          throw Rcpp::exception("Stage assignment error.", false);
        }
      } // stassign if statement
    } // currentindiv if statement
  } // i loop
  
  // Now let's check to see if censor variables have been fully and properly assigned - new section
  if (censorcol != -1) {
    arma::uvec censorzeros = find(censor2check == 0);
    
    for (int i = 0; i < noindivs; i++) {
      arma::uvec indivindices = sort(find(indivnum == i));
      
      arma::uvec indivzeros = intersect(indivindices, censorzeros);
      arma::ivec years_utilized = year2.elem(indivindices);
      
      int newcount = indivzeros.n_elem;
      
      if (newcount > 0) {
        for (int j = 0; j < newcount; j++) {
          if (censbool) {
            censor2(indivzeros(j)) = 0.0;
          } else {
            censor2(indivzeros(j)) = censorkeep;
          }
          
          arma::uvec currenttracking = find(indivindices == indivzeros(j));
          
          if (currenttracking.n_elem > 0) {
            
            int yearnow = year2(indivindices(currenttracking(0)));
            int yearprior = yearnow - 1;
            int yearnext = yearnow + 1;
            
            arma::uvec priorvec = find(years_utilized == yearprior);
            arma::uvec nextvec = find(years_utilized == yearnext);
            
            if (nextvec.n_elem > 0) {
              if (censbool) {
                censor1(nextvec(0)) = 0.0;
              } else {
                censor1(nextvec(0)) = censorkeep;
              }
            }
            
            if (priorvec.n_elem > 0) {
              if (censbool) {
                censor3(priorvec(0)) = 0.0;
              } else {
                censor3(priorvec(0)) = censorkeep;
              }
            } // priorvec if statement
          } // currenttracking if statement
        } // newcount for loop
      } // newcount if statement
    } // noindivs for loop
  } // end censor correction section
  
  // Output list creation
  output_longlist(0) = rowid;
  output_longlist(1) = popid;
  output_longlist(2) = patchid;
  output_longlist(3) = individ;
  output_longlist(4) = Rcpp::IntegerVector(year2.begin(), year2.end());
  output_longlist(5) = Rcpp::IntegerVector(firstseen.begin(), firstseen.end());
  output_longlist(6) = Rcpp::IntegerVector(lastseen.begin(), lastseen.end());
  output_longlist(7) = Rcpp::IntegerVector(obsage.begin(), obsage.end());
  output_longlist(8) = Rcpp::IntegerVector(obslifespan.begin(), obslifespan.end());
  
  output_longlist(9) = xpos1;
  output_longlist(10) = ypos1;
  output_longlist(11) = sizea1;
  output_longlist(12) = sizeb1;
  output_longlist(13) = sizec1;
  output_longlist(14) = sizeadded1;
  output_longlist(15) = repstra1;
  output_longlist(16) = repstrb1;
  output_longlist(17) = repstradded1;
  output_longlist(18) = feca1;
  output_longlist(19) = fecb1;
  output_longlist(20) = fecadded1;
  output_longlist(21) = indcova1;
  output_longlist(22) = indcovb1;
  output_longlist(23) = indcovc1;
  output_longlist(24) = censor1;
  output_longlist(25) = juvgiven1;
  
  output_longlist(26) = obsstatus1;
  output_longlist(27) = repstatus1;
  output_longlist(28) = fecstatus1;
  output_longlist(29) = matstat1;
  output_longlist(30) = alive1;
  output_longlist(31) = stage1;
  output_longlist(32) = stage1num;
  
  output_longlist(33) = xpos2;
  output_longlist(34) = ypos2;
  output_longlist(35) = sizea2;
  output_longlist(36) = sizeb2;
  output_longlist(37) = sizec2;
  output_longlist(38) = sizeadded2;
  output_longlist(39) = repstra2;
  output_longlist(40) = repstrb2;
  output_longlist(41) = repstradded2;
  output_longlist(42) = feca2;
  output_longlist(43) = fecb2;
  output_longlist(44) = fecadded2;
  output_longlist(45) = indcova2;
  output_longlist(46) = indcovb2;
  output_longlist(47) = indcovc2;
  output_longlist(48) = censor2;
  output_longlist(49) = juvgiven2;
  
  output_longlist(50) = obsstatus2;
  output_longlist(51) = repstatus2;
  output_longlist(52) = fecstatus2;
  output_longlist(53) = matstat2;
  output_longlist(54) = alive2;
  output_longlist(55) = stage2;
  output_longlist(56) = stage2num;

  output_longlist(57) = xpos3;
  output_longlist(58) = ypos3;
  output_longlist(59) = sizea3;
  output_longlist(60) = sizeb3;
  output_longlist(61) = sizec3;
  output_longlist(62) = sizeadded3;
  output_longlist(63) = repstra3;
  output_longlist(64) = repstrb3;
  output_longlist(65) = repstradded3;
  output_longlist(66) = feca3;
  output_longlist(67) = fecb3;
  output_longlist(68) = fecadded3;
  output_longlist(69) = indcova3;
  output_longlist(70) = indcovb3;
  output_longlist(71) = indcovc3;
  output_longlist(72) = censor3;
  output_longlist(73) = juvgiven3;
  
  output_longlist(74) = obsstatus3;
  output_longlist(75) = repstatus3;
  output_longlist(76) = fecstatus3;
  output_longlist(77) = matstat3;
  output_longlist(78) = alive3;
  output_longlist(79) = stage3;
  output_longlist(80) = stage3num;
  
  Rcpp::CharacterVector varnames {"rowid", "popid", "patchid", "individ",
    "year2", "firstseen", "lastseen", "obsage", "obslifespan","xpos1", "ypos1",
    "sizea1", "sizeb1", "sizec1", "size1added", "repstra1", "repstrb1",
    "repstr1added", "feca1", "fecb1", "fec1added", "indcova1", "indcovb1",
    "indcovc1", "censor1", "juvgiven1", "obsstatus1", "repstatus1",
    "fecstatus1", "matstatus1", "alive1", "stage1", "stage1index",
    "xpos2", "ypos2", "sizea2", "sizeb2", "sizec2", "size2added", "repstra2",
    "repstrb2", "repstr2added", "feca2", "fecb2", "fec2added", "indcova2",
    "indcovb2", "indcovc2", "censor2", "juvgiven2", "obsstatus2", "repstatus2",
    "fecstatus2", "matstatus2", "alive2", "stage2", "stage2index",
    "xpos3", "ypos3", "sizea3", "sizeb3", "sizec3", "size3added", "repstra3",
    "repstrb3", "repstr3added", "feca3", "fecb3", "fec3added", "indcova3",
    "indcovb3", "indcovc3", "censor3", "juvgiven3", "obsstatus3", "repstatus3",
    "fecstatus3", "matstatus3", "alive3", "stage3", "stage3index"};
  
  output_longlist.attr("names") = varnames;
  output_longlist.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, alive3.length());
  StringVector needed_classes {"data.frame", "hfvdata"};
  output_longlist.attr("class") = needed_classes; // data.frame
  
  return output_longlist;
}

//' Estimate Radial Density in Cartesian Space
//' 
//' Function \code{.density3()} estimates radial density on the basis of
//' Cartesian coordinates and spacing information supplied as input. It is used
//' internally by \code{\link{historicalize3}()} and
//' \code{\link{verticalize3}()}.
//' 
//' @param data Demographic dataset in historical vertical format.
//' @param xcol Number of column in \code{data} corresponding to x position.
//' @param ycol Number of column in \code{data} corresponding to y position.
//' @param yearcol Number of column in \code{data} corresponding to occasion
//' \emph{t}.
//' @param spacing Resolution of density estimation, as a scalar numeric.
//' 
//' @return This function returns a vector counting the number of individuals
//' within the specified distance of each individual in the historically
//' formatted vertical dataset.
//' 
//' @section Notes:
//' The process used to estimate density is one in which the distances between
//' all pairs of individuals are calculated via the Pythagorean theorem, and
//' then individual density equals the number of these individuals with
//' distances within the number input as \code{spacing}, respectively for each
//' individual.
//' 
//' This function assumes that all individuals are alive in time \emph{t}, and
//' so data should be filtered appropriately beforehand. Any rows with NA in X
//' or Y will not be counted, and density is estimated specific to time \emph{t}.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.density3)]]
Rcpp::NumericVector density3(Rcpp::DataFrame data, int xcol, int ycol, int yearcol,
  double spacing) {
  
  int data_size = data.length();
  int data_n = data.nrows();
  
  if (xcol < 0 || ycol < 0 || xcol > data_size || ycol > data_size) {
    throw Rcpp::exception("Input column numbers for X and/or Y are outside the range of the dataset", 
      false);
  }
  
  int xcol_true = xcol - 1;
  int ycol_true = ycol - 1;
  
  NumericVector Xdata = as<NumericVector>(data[xcol_true]);
  NumericVector Ydata = as<NumericVector>(data[ycol_true]);
  NumericVector yeardata = as<NumericVector>(data[(yearcol - 1)]);
  
  double ref_x {0.0};
  double ref_y {0.0};
  double est_a {0.0};
  double est_b {0.0};
  double est_c {0.0};
  int counted_n {0};
  
  arma::vec density(data_n);
  density.zeros();
  
  for (int i = 0; i < data_n; i++) {
    ref_x = Xdata(i);
    ref_y = Ydata(i);
    
    if (!NumericVector::is_na(Xdata(i)) && !NumericVector::is_na(Ydata(i))) {
      for (int j = 0; j < data_n; j++) {
        if (!NumericVector::is_na(Xdata(j)) && !NumericVector::is_na(Ydata(j)) && 
            yeardata(i) == yeardata(j)) {
          
          est_a = Xdata(j) - ref_x;
          est_b = Ydata(j) - ref_y;
          
          est_c = (est_a * est_a) + (est_b * est_b);
          est_c = sqrt(est_c);
          
          if (est_c < spacing) counted_n++;
        }
      }
    } else {
      counted_n = 1;
    }
    
    density(i) = counted_n - 1;
    counted_n = 0;
  }
  
  return Rcpp::NumericVector(density.begin(), density.end());
}

//' Create Element Index for Matrix Estimation
//' 
//' Function \code{.simplepizzle()} creates a data frame object used by function
//' \code{\link{hist_null}()} to provide an index for estimation of null
//' historical matrices from ahistorical MPM inputs.
//' 
//' @param StageFrame The stageframe object identifying the life history model
//' being operationalized.
//' @param format Integer indicating whether historical matrices should be in
//' (1) Ehrlen or (2) deVries format.
//' 
//' @return The output is composed of three elements:
//' \item{ahstages}{A new stageframe, which only differs from the input
//' stageframe in deVries format.}
//' \item{hstages}{A new historical stage-pair index for the new historical
//' matrices.}
//' \item{allstages}{A large data frame describing every element to be estimated
//' in the new historical matrices}.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.simplepizzle)]]
Rcpp::List simplepizzle(DataFrame StageFrame, int format) {
  
  arma::vec newstageid = as<arma::vec>(StageFrame["stage_id"]);
  StringVector origstageid = as<StringVector>(StageFrame["stage"]);
  arma::vec binsizectr = as<arma::vec>(StageFrame["sizebin_center"]);
  arma::vec repstatus = as<arma::vec>(StageFrame["repstatus"]);
  arma::vec obsstatus = as<arma::vec>(StageFrame["obsstatus"]);
  arma::vec immstatus = as<arma::vec>(StageFrame["immstatus"]);
  arma::vec matstatus = as<arma::vec>(StageFrame["matstatus"]);
  arma::vec indata = as<arma::vec>(StageFrame["indataset"]);
  arma::vec binsizewidth = as<arma::vec>(StageFrame["sizebin_width"]);
  arma::vec minage = as<arma::vec>(StageFrame["min_age"]);
  arma::vec maxage = as<arma::vec>(StageFrame["max_age"]);
  arma::vec group = as<arma::vec>(StageFrame["group"]);
  arma::vec entrystage = as<arma::vec>(StageFrame["entrystage"]);
  
  arma::vec binsizebctr = as<arma::vec>(StageFrame["sizebinb_center"]);
  arma::vec binsizecctr = as<arma::vec>(StageFrame["sizebinc_center"]);
  arma::vec binsizebwidth = as<arma::vec>(StageFrame["sizebinb_width"]);
  arma::vec binsizecwidth = as<arma::vec>(StageFrame["sizebinc_width"]);
  
  // This section determines the length of the matrix map data frame
  int nostages = newstageid.n_elem;
  int nostages_nounborn = nostages;
  int totallength {0};
  if (format == 2)  nostages = nostages + 1;
  
  arma::vec almostborn (nostages, fill::zeros);
  
  if (format == 2) {
    arma::vec oldorigsize = as<arma::vec>(StageFrame["original_size"]);
    arma::vec oldorigbsize = as<arma::vec>(StageFrame["original_size_b"]);
    arma::vec oldorigcsize = as<arma::vec>(StageFrame["original_size_c"]);
    arma::vec oldpropstatus = as<arma::vec>(StageFrame["propstatus"]);
    arma::vec oldbinhalfwidthraw = as<arma::vec>(StageFrame["binhalfwidth_raw"]);
    arma::vec oldsizebinmin = as<arma::vec>(StageFrame["sizebin_min"]);
    arma::vec oldsizebinmax = as<arma::vec>(StageFrame["sizebin_max"]);
    arma::vec oldbinhalfwidthbraw = as<arma::vec>(StageFrame["binhalfwidthb_raw"]);
    arma::vec oldsizebinbmin = as<arma::vec>(StageFrame["sizebinb_min"]);
    arma::vec oldsizebinbmax = as<arma::vec>(StageFrame["sizebinb_max"]);
    arma::vec oldbinhalfwidthcraw = as<arma::vec>(StageFrame["binhalfwidthc_raw"]);
    arma::vec oldsizebincmin = as<arma::vec>(StageFrame["sizebinc_min"]);
    arma::vec oldsizebincmax = as<arma::vec>(StageFrame["sizebinc_max"]);
    Rcpp::StringVector oldcomments = as<StringVector>(StageFrame["comments"]);
    arma::vec oldentrystage = as<arma::vec>(StageFrame["entrystage"]);
    
    nostages_nounborn = nostages - 1;
    totallength = (2 * nostages_nounborn * nostages_nounborn *
      nostages);
    
    arma::vec newstageidvec(nostages, fill::zeros);
    Rcpp::StringVector newstagevec(nostages);
    arma::vec neworigsizevec(nostages, fill::zeros);
    arma::vec neworigsizebvec(nostages, fill::zeros);
    arma::vec neworigsizecvec(nostages, fill::zeros);
    arma::vec newminagevec(nostages, fill::zeros);
    arma::vec newmaxagevec(nostages, fill::zeros);
    arma::vec newrepstatusvec(nostages, fill::zeros);
    arma::vec newobsstatusvec(nostages);
    arma::vec newpropstatusvec(nostages, fill::zeros);
    arma::vec newimmstatusvec(nostages, fill::zeros);
    arma::vec newmatstatusvec(nostages, fill::zeros);
    arma::vec newindatasetvec(nostages, fill::zeros);
    arma::vec newbinhalfwidthrawvec(nostages, fill::zeros);
    arma::vec newsizebinminvec(nostages, fill::zeros);
    arma::vec newsizebinmaxvec(nostages, fill::zeros);
    arma::vec newsizebincentervec(nostages, fill::zeros);
    arma::vec newsizebinwidthvec(nostages, fill::zeros);
    
    arma::vec newbinhalfwidthbrawvec(nostages, fill::zeros);
    arma::vec newsizebinbminvec(nostages, fill::zeros);
    arma::vec newsizebinbmaxvec(nostages, fill::zeros);
    arma::vec newsizebinbcentervec(nostages, fill::zeros);
    arma::vec newsizebinbwidthvec(nostages, fill::zeros);
    
    arma::vec newbinhalfwidthcrawvec(nostages, fill::zeros);
    arma::vec newsizebincminvec(nostages, fill::zeros);
    arma::vec newsizebincmaxvec(nostages, fill::zeros);
    arma::vec newsizebinccentervec(nostages, fill::zeros);
    arma::vec newsizebincwidthvec(nostages, fill::zeros);
    
    arma::vec newgroupvec(nostages, fill::zeros);
    Rcpp::StringVector newcomments(nostages);
    arma::vec newentrystage(nostages, fill::zeros);
    arma::vec newalmostborn(nostages, fill::zeros);
    
    for (int i = 0; i < nostages_nounborn; i++) {
      newstageidvec(i) = newstageid(i);
      newstagevec(i) = origstageid(i);
      neworigsizevec(i) = oldorigsize(i);
      neworigsizebvec(i) = oldorigbsize(i);
      neworigsizecvec(i) = oldorigcsize(i);
      newminagevec(i) = minage(i);
      newmaxagevec(i) = maxage(i);
      newrepstatusvec(i) = repstatus(i);
      newobsstatusvec(i) = obsstatus(i);
      newpropstatusvec(i) = oldpropstatus(i);
      newimmstatusvec(i) = immstatus(i);
      newmatstatusvec(i) = matstatus(i);
      newindatasetvec(i) = indata(i);
      
      newbinhalfwidthrawvec(i) = oldbinhalfwidthraw(i);
      newsizebinminvec(i) = oldsizebinmin(i);
      newsizebinmaxvec(i) = oldsizebinmax(i);
      newsizebincentervec(i) = binsizectr(i);
      newsizebinwidthvec(i) = binsizewidth(i);
      
      newbinhalfwidthbrawvec(i) = oldbinhalfwidthbraw(i);
      newsizebinbminvec(i) = oldsizebinbmin(i);
      newsizebinbmaxvec(i) = oldsizebinbmax(i);
      newsizebinbcentervec(i) = binsizebctr(i);
      newsizebinbwidthvec(i) = binsizebwidth(i);
      
      newbinhalfwidthcrawvec(i) = oldbinhalfwidthcraw(i);
      newsizebincminvec(i) = oldsizebincmin(i);
      newsizebincmaxvec(i) = oldsizebincmax(i);
      newsizebinccentervec(i) = binsizecctr(i);
      newsizebincwidthvec(i) = binsizecwidth(i);
      
      newgroupvec(i) = group(i);
      newcomments(i) = oldcomments(i);
      newentrystage(i) = oldentrystage(i);
    }
    
    newstageidvec(nostages_nounborn) = newstageid(nostages_nounborn - 1.0) + 1.0;
    newstagevec(nostages_nounborn) = "AlmostBorn";
    neworigsizevec(nostages_nounborn) = 0.0;
    neworigsizebvec(nostages_nounborn) = 0.0;
    neworigsizecvec(nostages_nounborn) = 0.0;
    newminagevec(nostages_nounborn) = NA_REAL;
    newmaxagevec(nostages_nounborn) = NA_REAL;
    newrepstatusvec(nostages_nounborn) = 0.0;
    newobsstatusvec(nostages_nounborn) = 1.0;
    newpropstatusvec(nostages_nounborn) = 0.0;
    newimmstatusvec(nostages_nounborn) = 1.0;
    newmatstatusvec(nostages_nounborn) = 0.0;
    newindatasetvec(nostages_nounborn) = 1.0;
    
    newbinhalfwidthrawvec(nostages_nounborn) = 0.0;
    newsizebinminvec(nostages_nounborn) = 0.0;
    newsizebinmaxvec(nostages_nounborn) = 0.0;
    newsizebincentervec(nostages_nounborn) = 0.0;
    newsizebinwidthvec(nostages_nounborn) = 0.0;
    
    newbinhalfwidthbrawvec(nostages_nounborn) = 0.0;
    newsizebinbminvec(nostages_nounborn) = 0.0;
    newsizebinbmaxvec(nostages_nounborn) = 0.0;
    newsizebinbcentervec(nostages_nounborn) = 0.0;
    newsizebinbwidthvec(nostages_nounborn) = 0.0;
    
    newbinhalfwidthcrawvec(nostages_nounborn) = 0.0;
    newsizebincminvec(nostages_nounborn) = 0.0;
    newsizebincmaxvec(nostages_nounborn) = 0.0;
    newsizebinccentervec(nostages_nounborn) = 0.0;
    newsizebincwidthvec(nostages_nounborn) = 0.0;
    
    newgroupvec(nostages_nounborn) = 0.0;
    newcomments(nostages_nounborn) = "Almost Born";
    newentrystage(nostages_nounborn) = 0.0;
    newalmostborn(nostages_nounborn) = 1.0;
    
    Rcpp::List new_stageframe(32);
    
    new_stageframe(0) = newstageidvec;
    new_stageframe(1) = newstagevec;
    new_stageframe(2) = neworigsizevec;
    new_stageframe(3) = neworigsizebvec;
    new_stageframe(4) = neworigsizecvec;
    new_stageframe(5) = newminagevec;
    new_stageframe(6) = newmaxagevec;
    new_stageframe(7) = newrepstatusvec;
    new_stageframe(8) = newobsstatusvec;
    new_stageframe(9) = newpropstatusvec;
    new_stageframe(10) = newimmstatusvec;
    new_stageframe(11) = newmatstatusvec;
    new_stageframe(12) = newindatasetvec;
    
    new_stageframe(13) = newbinhalfwidthrawvec;
    new_stageframe(14) = newsizebinminvec;
    new_stageframe(15) = newsizebinmaxvec;
    new_stageframe(16) = newsizebincentervec;
    new_stageframe(17) = newsizebinwidthvec;
    
    new_stageframe(18) = newbinhalfwidthbrawvec;
    new_stageframe(19) = newsizebinbminvec;
    new_stageframe(20) = newsizebinbmaxvec;
    new_stageframe(21) = newsizebinbcentervec;
    new_stageframe(22) = newsizebinbwidthvec;
    
    new_stageframe(23) = newbinhalfwidthcrawvec;
    new_stageframe(24) = newsizebincminvec;
    new_stageframe(25) = newsizebincmaxvec;
    new_stageframe(26) = newsizebinccentervec;
    new_stageframe(27) = newsizebincwidthvec;
    
    new_stageframe(28) = newgroupvec;
    new_stageframe(29) = newcomments;
    new_stageframe(30) = newentrystage;
    new_stageframe(31) = newalmostborn;
    
    CharacterVector sfnamevec = {"stage_id", "stage", "original_size", "size_b",
      "size_c", "min_age", "max_age", "repstatus", "obsstatus", "propstatus",
      "immstatus", "matstatus", "indataset", "binhalfwidth_raw", "sizebin_min",
      "sizebin_max", "sizebin_center", "sizebin_width", "binhalfwidthb_raw",
      "sizebinb_min", "sizebinb_max", "sizebinb_center", "sizebinb_width",
      "binhalfwidthc_raw", "sizebinc_min", "sizebinc_max", "sizebinc_center",
      "sizebinc_width", "group", "comments", "entrystage", "almostborn"};
    
    new_stageframe.attr("names") = sfnamevec;
    new_stageframe.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, newstageidvec.n_elem);
    new_stageframe.attr("class") = "data.frame";
    
    StageFrame = new_stageframe;
    
    newstageid = as<arma::vec>(new_stageframe["stage_id"]);
    origstageid = as<StringVector>(new_stageframe["stage"]);
    binsizectr = as<arma::vec>(new_stageframe["sizebin_center"]);
    repstatus = as<arma::vec>(new_stageframe["repstatus"]);
    obsstatus = as<arma::vec>(new_stageframe["obsstatus"]);
    immstatus = as<arma::vec>(new_stageframe["immstatus"]);
    matstatus = as<arma::vec>(new_stageframe["matstatus"]);
    indata = as<arma::vec>(new_stageframe["indataset"]);
    binsizewidth = as<arma::vec>(new_stageframe["sizebin_width"]);
    minage = as<arma::vec>(new_stageframe["min_age"]);
    maxage = as<arma::vec>(new_stageframe["max_age"]);
    group = as<arma::vec>(new_stageframe["group"]);
    
    binsizebctr = as<arma::vec>(new_stageframe["sizebinb_center"]);
    binsizecctr = as<arma::vec>(new_stageframe["sizebinc_center"]);
    binsizebwidth = as<arma::vec>(new_stageframe["sizebinb_width"]);
    binsizecwidth = as<arma::vec>(new_stageframe["sizebinc_width"]);
    almostborn = as<arma::vec>(new_stageframe["almostborn"]);
    entrystage = as<arma::vec>(new_stageframe["entrystage"]);
    
  } else {
    totallength = (nostages * nostages * nostages);
  }
  
  // New let's create the new hstages
  Rcpp::List hstages(4);
  
  int hstages_length = nostages * nostages_nounborn;
  
  Rcpp::IntegerVector stid2(hstages_length);
  Rcpp::IntegerVector stid1(hstages_length);
  Rcpp::StringVector st2(hstages_length);
  Rcpp::StringVector st1(hstages_length);
  
  for (int j = 0; j < nostages; j++) {
    for (int i = 0; i < nostages_nounborn; i++) {
      stid2((j * (nostages_nounborn)) + i) = i + 1;
      stid1((j * (nostages_nounborn)) + i) = j + 1;
      st2((j * (nostages_nounborn)) + i) = origstageid(i);
      st1((j * (nostages_nounborn)) + i) = origstageid(j);
    }
  }
  
  hstages(0) = stid2;
  hstages(1) = stid1;
  hstages(2) = st2;
  hstages(3) = st1;
  
  CharacterVector hsnamevec = {"stage_id_2", "stage_id_1", "stage_2", "stage_1"};
  hstages.attr("names") = hsnamevec;
  hstages.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, stid2.length());
  hstages.attr("class") = "data.frame";
  
  // Here we set up the vectors that will be put together into the matrix map
  // data frame
  arma::vec stage3(totallength, fill::zeros);
  arma::vec stage2n(totallength, fill::zeros);
  arma::vec stage2o(totallength, fill::zeros);
  arma::vec stage1(totallength, fill::zeros);
  
  arma::vec size3(totallength, fill::zeros);
  arma::vec size2n(totallength, fill::zeros);
  arma::vec size2o(totallength, fill::zeros);
  arma::vec size1(totallength, fill::zeros);
  
  arma::vec sizeb3(totallength, fill::zeros);
  arma::vec sizeb2n(totallength, fill::zeros);
  arma::vec sizeb2o(totallength, fill::zeros);
  arma::vec sizeb1(totallength, fill::zeros);
  
  arma::vec sizec3(totallength, fill::zeros);
  arma::vec sizec2n(totallength, fill::zeros);
  arma::vec sizec2o(totallength, fill::zeros);
  arma::vec sizec1(totallength, fill::zeros);
  
  arma::vec obs3(totallength, fill::zeros);
  arma::vec obs2n(totallength, fill::zeros);
  arma::vec obs2o(totallength, fill::zeros);
  arma::vec obs1(totallength, fill::zeros);
  
  arma::vec rep3(totallength, fill::zeros);
  arma::vec rep2n(totallength, fill::zeros);
  arma::vec rep2o(totallength, fill::zeros);
  arma::vec rep1(totallength, fill::zeros);
  
  arma::vec mat3(totallength, fill::zeros);
  arma::vec mat2n(totallength, fill::zeros);
  arma::vec mat2o(totallength, fill::zeros);
  arma::vec mat1(totallength, fill::zeros);
  
  arma::vec imm3(totallength, fill::zeros);
  arma::vec imm2n(totallength, fill::zeros);
  arma::vec imm2o(totallength, fill::zeros);
  arma::vec imm1(totallength, fill::zeros);
  
  arma::vec repentry3(totallength, fill::zeros);
  arma::vec repentry2o(totallength, fill::zeros);
  arma::vec almostborn2n(totallength, fill::zeros);
  arma::vec almostborn1(totallength, fill::zeros);
  
  arma::vec binwidth(totallength, fill::zeros);
  arma::vec binbwidth(totallength, fill::zeros);
  arma::vec bincwidth(totallength, fill::zeros);
  
  arma::vec indata3(totallength, fill::zeros);
  arma::vec indata2n(totallength, fill::zeros);
  arma::vec indata2o(totallength, fill::zeros);
  arma::vec indata1(totallength, fill::zeros);
  
  arma::vec minage3(totallength, fill::zeros);
  arma::vec minage2(totallength, fill::zeros);
  arma::vec maxage3(totallength, fill::zeros);
  arma::vec maxage2(totallength, fill::zeros);
  
  arma::vec grp3(totallength, fill::zeros);
  arma::vec grp2n(totallength, fill::zeros);
  arma::vec grp2o(totallength, fill::zeros);
  arma::vec grp1(totallength, fill::zeros);
  
  arma::vec actualage(totallength, fill::zeros);
  arma::vec index321u(totallength);
  arma::vec index321f(totallength);
  arma::vec index21(totallength);
  arma::vec indatalong(totallength, fill::zeros);
  arma::vec aliveequal(totallength);
  arma::vec included(totallength, fill::zeros);
  index321u.fill(-1);
  index321f.fill(-1);
  index21.fill(-1);
  aliveequal.fill(-1);
  
  long long int currentindex {0};
  
  // Now we cover the main data frame creation loops
  // When style = 0, this will create AllStages for the historical case
  if (format == 2) {
    for (int time1 = 0; time1 < nostages; time1++) {
      for (int time2o = 0; time2o < nostages_nounborn; time2o++) {
        for (int time2n = 0; time2n < nostages; time2n++) {
          for (int time3 = 0; time3 < nostages; time3++) {
            
            if (time3 != nostages_nounborn) {
              if (time2n == time2o || time2n == nostages_nounborn) {
                
                included(currentindex) = 1;
                
                stage3(currentindex) = newstageid(time3);
                stage2n(currentindex) = newstageid(time2n);
                stage2o(currentindex) = newstageid(time2o);
                stage1(currentindex) = newstageid(time1);
                
                size3(currentindex) = binsizectr(time3);
                size2n(currentindex) = binsizectr(time2n);
                size2o(currentindex) = binsizectr(time2o);
                size1(currentindex) = binsizectr(time1);
                
                sizeb3(currentindex) = binsizebctr(time3);
                sizeb2n(currentindex) = binsizebctr(time2n);
                sizeb2o(currentindex) = binsizebctr(time2o);
                sizeb1(currentindex) = binsizebctr(time1);
                
                if (NumericVector::is_na(sizeb3(currentindex))) sizeb3(currentindex) = 0;
                if (NumericVector::is_na(sizeb2n(currentindex))) sizeb2n(currentindex) = 0;
                if (NumericVector::is_na(sizeb2o(currentindex))) sizeb2o(currentindex) = 0;
                if (NumericVector::is_na(sizeb1(currentindex))) sizeb1(currentindex) = 0;
                
                sizec3(currentindex) = binsizecctr(time3);
                sizec2n(currentindex) = binsizecctr(time2n);
                sizec2o(currentindex) = binsizecctr(time2o);
                sizec1(currentindex) = binsizecctr(time1);
                
                if (NumericVector::is_na(sizec3(currentindex))) sizec3(currentindex) = 0;
                if (NumericVector::is_na(sizec2n(currentindex))) sizec2n(currentindex) = 0;
                if (NumericVector::is_na(sizec2o(currentindex))) sizec2o(currentindex) = 0;
                if (NumericVector::is_na(sizec1(currentindex))) sizec1(currentindex) = 0;
                
                obs3(currentindex) = obsstatus(time3);
                obs2n(currentindex) = obsstatus(time2n);
                obs2o(currentindex) = obsstatus(time2o);
                obs1(currentindex) = obsstatus(time1);
                
                rep3(currentindex) = repstatus(time3);
                rep2n(currentindex) = repstatus(time2n);
                rep2o(currentindex) = repstatus(time2o);
                rep1(currentindex) = repstatus(time1);
                
                mat3(currentindex) = matstatus(time3);
                mat2n(currentindex) = matstatus(time2n);
                mat2o(currentindex) = matstatus(time2o);
                mat1(currentindex) = matstatus(time1);
                
                imm3(currentindex) = immstatus(time3);
                imm2n(currentindex) = immstatus(time2n);
                imm2o(currentindex) = immstatus(time2o);
                imm1(currentindex) = immstatus(time1);
                
                repentry3(currentindex) = entrystage(time3);
                repentry2o(currentindex) = entrystage(time2o);
                almostborn2n(currentindex) = almostborn(time2n);
                almostborn1(currentindex) = almostborn(time1);
                
                indata3(currentindex) = indata(time3);
                indata2n(currentindex) = indata(time2n);
                indata2o(currentindex) = indata(time2o);
                indata1(currentindex) = indata(time1);
                
                binwidth(currentindex) = binsizewidth(time3);
                binbwidth(currentindex) = binsizebwidth(time3);
                bincwidth(currentindex) = binsizecwidth(time3);
                
                if (NumericVector::is_na(binbwidth(currentindex))) binbwidth(currentindex) = 0;
                if (NumericVector::is_na(bincwidth(currentindex))) bincwidth(currentindex) = 0;
                
                minage3(currentindex) = minage(time3);
                minage2(currentindex) = minage(time2o);
                maxage3(currentindex) = maxage(time3);
                maxage2(currentindex) = maxage(time2o);
                actualage(currentindex) = 0;
                
                grp3(currentindex) = group(time3);
                grp2n(currentindex) = group(time2n);
                grp2o(currentindex) = group(time2o);
                grp1(currentindex) = group(time1);
                
                aliveequal(currentindex) = (stage3(currentindex) - 1) + 
                  ((stage2n(currentindex) - 1) * nostages_nounborn) + 
                  ((stage2o(currentindex) - 1) * nostages * nostages_nounborn) + 
                  ((stage1(currentindex) - 1) * nostages_nounborn * nostages *
                    nostages_nounborn);
                
                index321u(currentindex) = (stage3(currentindex) - 1) + 
                  ((stage2n(currentindex) - 1) * nostages_nounborn) + 
                  ((stage2o(currentindex) - 1) * nostages * nostages_nounborn) + 
                  ((stage1(currentindex) - 1) * nostages_nounborn * nostages *
                    nostages_nounborn);
                index321f(currentindex) = index321u(currentindex);
                  
                index21(currentindex) = (stage3(currentindex) - 1) + 
                  ((stage2o(currentindex) - 1) * nostages_nounborn); // Used to be 2o and 1
                
                indatalong(currentindex) = indata3(currentindex) * indata2n(currentindex) * 
                  indata2o(currentindex) * indata1(currentindex);
                
                // Now a few corrections
                if ((almostborn1(currentindex) > 0.0 && repentry2o(currentindex) < 1.0)) {
                    index321u(currentindex) = -1.0;
                    index321f(currentindex) = -1.0;
                }
                if (almostborn2n(currentindex) > 0.0) {
                  if (repentry3(currentindex) < 1.0 || rep2o(currentindex) < 1.0) {
                    index321u(currentindex) = -1.0;
                    index321f(currentindex) = -1.0;
                  }
                }
                if (time2n == time2o) {
                  index321f(currentindex) = -1.0;
                }
                currentindex += 1;
              } // if (time2n == tim2o || time2n == nostages_nounborn) statement
            } // if (time3n != dead_stage) statement
          } // time3 loop
        } // time2n loop
      } // time2o loop
    } // time1 loop 
    
  } else if (format == 1) { // Historical MPM in Ehrlen format
    for (int time1 = 0; time1 < nostages; time1++) {
      for (int time2o = 0; time2o < nostages; time2o++) {
        for (int time3 = 0; time3 < nostages; time3++) {
          
          included(currentindex) = 1;
          
          stage3(currentindex) = newstageid(time3);
          stage2n(currentindex) = newstageid(time2o);
          stage2o(currentindex) = newstageid(time2o);
          stage1(currentindex) = newstageid(time1);
          
          size3(currentindex) = binsizectr(time3);
          size2n(currentindex) = binsizectr(time2o);
          size2o(currentindex) = binsizectr(time2o);
          size1(currentindex) = binsizectr(time1);
          
          sizeb3(currentindex) = binsizebctr(time3);
          sizeb2n(currentindex) = binsizebctr(time2o);
          sizeb2o(currentindex) = binsizebctr(time2o);
          sizeb1(currentindex) = binsizebctr(time1);
          
          if (NumericVector::is_na(sizeb3(currentindex))) sizeb3(currentindex) = 0;
          if (NumericVector::is_na(sizeb2n(currentindex))) sizeb2n(currentindex) = 0;
          if (NumericVector::is_na(sizeb2o(currentindex))) sizeb2o(currentindex) = 0;
          if (NumericVector::is_na(sizeb1(currentindex))) sizeb1(currentindex) = 0;
                
          sizec3(currentindex) = binsizecctr(time3);
          sizec2n(currentindex) = binsizecctr(time2o);
          sizec2o(currentindex) = binsizecctr(time2o);
          sizec1(currentindex) = binsizecctr(time1);
          
          if (NumericVector::is_na(sizec3(currentindex))) sizec3(currentindex) = 0;
          if (NumericVector::is_na(sizec2n(currentindex))) sizec2n(currentindex) = 0;
          if (NumericVector::is_na(sizec2o(currentindex))) sizec2o(currentindex) = 0;
          if (NumericVector::is_na(sizec1(currentindex))) sizec1(currentindex) = 0;
          
          obs3(currentindex) = obsstatus(time3);
          obs2n(currentindex) = obsstatus(time2o);
          obs2o(currentindex) = obsstatus(time2o);
          obs1(currentindex) = obsstatus(time1);
          
          rep3(currentindex) = repstatus(time3);
          rep2n(currentindex) = repstatus(time2o);
          rep2o(currentindex) = repstatus(time2o);
          rep1(currentindex) = repstatus(time1);
          
          mat3(currentindex) = matstatus(time3);
          mat2n(currentindex) = matstatus(time2o);
          mat2o(currentindex) = matstatus(time2o);
          mat1(currentindex) = matstatus(time1);
          
          imm3(currentindex) = immstatus(time3);
          imm2n(currentindex) = immstatus(time2o);
          imm2o(currentindex) = immstatus(time2o);
          imm1(currentindex) = immstatus(time1);
          
          repentry3(currentindex) = 0;
          repentry2o(currentindex) = entrystage(time2o);
          almostborn1(currentindex) = almostborn(time1);
          
          indata3(currentindex) = indata(time3);
          indata2n(currentindex) = indata(time2o);
          indata2o(currentindex) = indata(time2o);
          indata1(currentindex) = indata(time1);
          
          binwidth(currentindex) = binsizewidth(time3);
          binbwidth(currentindex) = binsizebwidth(time3);
          bincwidth(currentindex) = binsizecwidth(time3);
          
          if (NumericVector::is_na(binbwidth(currentindex))) binbwidth(currentindex) = 0;
          if (NumericVector::is_na(bincwidth(currentindex))) bincwidth(currentindex) = 0;
          
          minage3(currentindex) = minage(time3);
          minage2(currentindex) = minage(time2o);
          maxage3(currentindex) = maxage(time3);
          maxage2(currentindex) = maxage(time2o);
          actualage(currentindex) = 0;
          
          grp3(currentindex) = group(time3);
          grp2n(currentindex) = group(time2o);
          grp2o(currentindex) = group(time2o);
          grp1(currentindex) = group(time1);
          
          aliveequal(currentindex) = (stage3(currentindex) - 1) + ((stage2n(currentindex) - 1) * 
              (nostages - 1)) + ((stage2o(currentindex) - 1) * (nostages - 1) * (nostages - 1)) + 
            ((stage1(currentindex) - 1) * (nostages - 1) * (nostages - 1) * (nostages - 1));
            
          index321u(currentindex) = (stage3(currentindex) - 1) + 
            ((stage2n(currentindex) - 1) * nostages) + 
            ((stage2n(currentindex) - 1) * nostages * nostages) + 
            ((stage1(currentindex) - 1) * nostages * nostages * nostages);
          index321f(currentindex) = index321u(currentindex);
          index21(currentindex) = (stage3(currentindex) - 1) + ((stage2n(currentindex) - 1) * nostages);
          
          indatalong(currentindex) = indata3(currentindex) * indata2n(currentindex) * 
            indata2o(currentindex) * indata1(currentindex);
          
          currentindex += 1;
        } // time3 loop
      } // time2o loop
    } // time1 loop 
  }
  
  int stage3_length = stage3.n_elem;
  
  Rcpp::List output_longlist(63);
  
  output_longlist(0) = Rcpp::NumericVector(stage3.begin(), stage3.end());
  output_longlist(1) = Rcpp::NumericVector(stage2n.begin(), stage2n.end());
  output_longlist(2) = Rcpp::NumericVector(stage2o.begin(), stage2o.end());
  output_longlist(3) = Rcpp::NumericVector(stage1.begin(), stage1.end());
  output_longlist(4) = Rcpp::NumericVector(size3.begin(), size3.end());
  output_longlist(5) = Rcpp::NumericVector(size2n.begin(), size2n.end());
  output_longlist(6) = Rcpp::NumericVector(size2o.begin(), size2o.end());
  output_longlist(7) = Rcpp::NumericVector(size1.begin(), size1.end());
  output_longlist(8) = Rcpp::NumericVector(sizeb3.begin(), sizeb3.end());
  output_longlist(9) = Rcpp::NumericVector(sizeb2n.begin(), sizeb2n.end());
  
  output_longlist(10) = Rcpp::NumericVector(sizeb2o.begin(), sizeb2o.end());
  output_longlist(11) = Rcpp::NumericVector(sizeb1.begin(), sizeb1.end());
  output_longlist(12) = Rcpp::NumericVector(sizec3.begin(), sizec3.end());
  output_longlist(13) = Rcpp::NumericVector(sizec2n.begin(), sizec2n.end());
  output_longlist(14) = Rcpp::NumericVector(sizec2o.begin(), sizec2o.end());
  output_longlist(15) = Rcpp::NumericVector(sizec1.begin(), sizec1.end());
  output_longlist(16) = Rcpp::NumericVector(obs3.begin(), obs3.end());
  output_longlist(17) = Rcpp::NumericVector(obs2n.begin(), obs2n.end());
  output_longlist(18) = Rcpp::NumericVector(obs2o.begin(), obs2o.end());
  output_longlist(19) = Rcpp::NumericVector(obs1.begin(), obs1.end());
  
  output_longlist(20) = Rcpp::NumericVector(rep3.begin(), rep3.end());
  output_longlist(21) = Rcpp::NumericVector(rep2n.begin(), rep2n.end());
  output_longlist(22) = Rcpp::NumericVector(rep2o.begin(), rep2o.end());
  output_longlist(23) = Rcpp::NumericVector(rep1.begin(), rep1.end());
  output_longlist(24) = Rcpp::NumericVector(mat3.begin(), mat3.end());
  output_longlist(25) = Rcpp::NumericVector(mat2n.begin(), mat2n.end());
  output_longlist(26) = Rcpp::NumericVector(mat2o.begin(), mat2o.end());
  output_longlist(27) = Rcpp::NumericVector(mat1.begin(), mat1.end());
  output_longlist(28) = Rcpp::NumericVector(imm3.begin(), imm3.end());
  output_longlist(29) = Rcpp::NumericVector(imm2n.begin(), imm2n.end());
  
  output_longlist(30) = Rcpp::NumericVector(imm2o.begin(), imm2o.end());
  output_longlist(31) = Rcpp::NumericVector(imm1.begin(), imm1.end());
  output_longlist(32) = Rcpp::NumericVector(repentry3.begin(), repentry3.end());
  output_longlist(33) = Rcpp::NumericVector(indata3.begin(), indata3.end());
  output_longlist(34) = Rcpp::NumericVector(indata2n.begin(), indata2n.end());
  output_longlist(35) = Rcpp::NumericVector(indata2o.begin(), indata2o.end());
  output_longlist(36) = Rcpp::NumericVector(indata1.begin(), indata1.end());
  output_longlist(37) = Rcpp::NumericVector(binwidth.begin(), binwidth.end());
  output_longlist(38) = Rcpp::NumericVector(binbwidth.begin(), binbwidth.end());
  output_longlist(39) = Rcpp::NumericVector(bincwidth.begin(), bincwidth.end());
  
  output_longlist(40) = Rcpp::NumericVector(minage3.begin(), minage3.end());
  output_longlist(41) = Rcpp::NumericVector(minage2.begin(), minage2.end());
  output_longlist(42) = Rcpp::NumericVector(maxage3.begin(), maxage3.end());
  output_longlist(43) = Rcpp::NumericVector(maxage2.begin(), maxage2.end());
  output_longlist(44) = Rcpp::NumericVector(actualage.begin(), actualage.end());
  
  output_longlist(45) = Rcpp::NumericVector(grp3.begin(), grp3.end());
  output_longlist(46) = Rcpp::NumericVector(grp2n.begin(), grp2n.end());
  output_longlist(47) = Rcpp::NumericVector(grp2o.begin(), grp2o.end());
  output_longlist(48) = Rcpp::NumericVector(grp1.begin(), grp1.end());
  
  output_longlist(49) = Rcpp::NumericVector(indatalong.begin(), indatalong.end());
  output_longlist(50) = Rcpp::NumericVector(stage3_length, -1.0);
  output_longlist(51) = Rcpp::NumericVector(stage3_length, -1.0);
  output_longlist(52) = Rcpp::NumericVector(stage3_length, -1.0);
  output_longlist(53) = Rcpp::NumericVector(stage3_length, -1.0);
  output_longlist(54) = Rcpp::NumericVector(stage3_length, 1.0);
  output_longlist(55) = Rcpp::NumericVector(stage3_length, 1.0);
  
  output_longlist(56) = Rcpp::NumericVector(aliveequal.begin(), aliveequal.end());
  output_longlist(57) = Rcpp::NumericVector(index321u.begin(), index321u.end());
  output_longlist(58) = Rcpp::NumericVector(index21.begin(), index21.end());
  output_longlist(59) = Rcpp::NumericVector(repentry2o.begin(), repentry2o.end());
  output_longlist(60) = Rcpp::NumericVector(almostborn2n.begin(), almostborn2n.end());
  output_longlist(61) = Rcpp::NumericVector(almostborn1.begin(), almostborn1.end());
  output_longlist(62) = Rcpp::NumericVector(index321f.begin(), index321f.end());
  
  CharacterVector namevec = {"stage3", "stage2n", "stage2o", "stage1", "size3",
    "size2n", "size2o", "size1", "sizeb3", "sizeb2n", "sizeb2o", "sizeb1", 
    "sizec3", "sizec2n", "sizec2o", "sizec1", "obs3", "obs2n", "obs2o", "obs1",
    "rep3", "rep2n", "rep2o", "rep1", "mat3", "mat2n", "mat2o", "mat1", "imm3",
    "imm2n", "imm2o", "imm1", "repentry3", "indata3", "indata2n", "indata2o",
    "indata1", "binwidth", "binbwidth", "bincwidth", "minage3", "minage2",
    "maxage3", "maxage2", "actualage", "group3", "group2n", "group2o", "group1",
    "indata", "ovgiven_t", "ovest_t", "ovgiven_f", "ovest_f", "ovsurvmult",
    "ovfecmult", "aliveandequal", "index321u", "index21", "repentry2o",
    "almostborn2n", "almostborn1", "index321f"};
  output_longlist.attr("names") = namevec;
  output_longlist.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, stage3_length);
  output_longlist.attr("class") = "data.frame";
  
  Rcpp::List output = Rcpp::List::create(Named("ahstages") = StageFrame,
    _["hstages"] = hstages, _["allstages"] = output_longlist);
  return output;
}

