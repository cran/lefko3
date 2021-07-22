#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Re-index Projection Matrix On Basis of Overwrite Table
//' 
//' This function takes matrix indices provided by functions
//' \code{\link{rlefko3}()}, \code{\link{rlefko2}()}, \code{\link{flefko3}()},
//' \code{\link{flefko2}()}, and \code{\link{aflefko2}()} and updates them with
//' information provided in the overwrite table used as input in that function.
//' 
//' @param allst321 Vector containing the original element-by-element matrix
//' index.
//' @param idx321old Vector containing the indices of matrix elements to be
//' updated.
//' @param idx321new Vector containing the replacement matrix element indices.
//' @param convtype Vector denoting survival transition (1), fecundity (2), or
//' fecundity multiplier (3).
//' @param eststag3 Vector of new stages in time \emph{t}+1.
//' @param gvnrate Vector of replacement transition values.
//' @param multipl Vector of fecundity multipliers.
//' 
//' @return Vector of updated matrix indices and, where appropriate, replacement
//' matrix element values.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat ovreplace(arma::vec allst321, arma::vec idx321old,
  arma::vec idx321new, arma::vec convtype, arma::vec eststag3, 
  arma::vec gvnrate, arma::vec multipl) {
  
  int n = idx321new.n_elem;
  
  arma::mat replacements(allst321.n_elem, 5);
  replacements.fill(-1);
  
  for (int i = 0; i < n; i++) {
    arma::uvec correctplace = find(allst321 == idx321old[i]);
    
    int m = correctplace.n_elem; 
    
    for (int j = 0; j < m; j++) {
      if (convtype[i] == 1) {
        if (gvnrate[i] >= 0) {replacements(correctplace[j], 0) = gvnrate[i];}
        if (eststag3[i] != -1 && idx321new[i] >= 0) {replacements(correctplace[j], 1) = idx321new[i];}
      }
      
      if (convtype[i] == 2) {
        if (gvnrate[i] >= 0) {replacements(correctplace[j], 2) = gvnrate[i];}
        if (eststag3[i] != -1 && idx321new[i] >= 0) {replacements(correctplace[j], 3) = idx321new[i];}
      }
      
      if (convtype[i] == 3) {
        if (multipl[i] > 0) {replacements(correctplace[j], 4) = multipl[i];}
      }
    }
  }
  
  return replacements;
}

//' Create Vertical Structure for Horizontal Data Frame Input
//' 
//' Function \code{pfj()} powers the R function \code{\link{verticalize3}()},
//' creating the vertical structure and rearranging the data in that shape.
//' 
//' @param data The horizontal data file.
//' @param stageframe The stageframe object identifying the life history model
//' being operationalized. This should be the full stageframe.
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
//' @param stassign A logical value indicating whether to assign stages.
//' @param stszcol Column number describing which size variable to use in stage 
//' estimation.
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
//' 
//' @return The output is currently a 7 element list, where each element is a
//' data frame with the same number of rows.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List pfj(DataFrame data, DataFrame stageframe, int noyears, int firstyear,
  int popidcol, int patchidcol, int individcol, int blocksize, arma::ivec xcol,
  arma::ivec ycol, arma::ivec juvcol, arma::ivec sizeacol, arma::ivec sizebcol,
  arma::ivec sizeccol, arma::ivec repstracol, arma::ivec repstrbcol,
  arma::ivec fecacol, arma::ivec fecbcol, arma::ivec indcovacol,
  arma::ivec indcovbcol, arma::ivec indcovccol, arma::ivec aliveacol,
  arma::ivec deadacol, arma::ivec obsacol, arma::ivec nonobsacol,
  arma::ivec censorcol, arma::ivec stagecol, double repstrrel, double fecrel,
  bool NAas0, bool NRasRep, bool RepasObs, bool stassign, int stszcol,
  double censorkeep, bool censbool, bool censrepeat, bool coordsrepeat) {
  
  // Here we standardize variable vector lengths  
  if (xcol.n_elem == 1) {
    xcol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        xcol(i) = xcol(0);
      }
    } else {
      for (int i = 1; i < noyears; i++) {
        xcol(i) = xcol(i-1) + blocksize;
      }
    }
  }
  
  if (ycol.n_elem == 1) {
    ycol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        ycol(i) = ycol(0);
      }
    } else {
      for (int i = 1; i < noyears; i++) {
        ycol(i) = ycol(i-1) + blocksize;
      }
    }
  }
  
  if (juvcol.n_elem == 1) {
    juvcol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        juvcol(i) = juvcol(0);
      }
    } else {
      for (int i = 1; i < noyears; i++) {
        juvcol(i) = juvcol(i-1) + blocksize;
      }
    }
  }
  
  if (sizeacol.n_elem == 1) {
    sizeacol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        sizeacol(i) = sizeacol(0);
      }
    } else {
      for (int i = 1; i < noyears; i++) {
        sizeacol(i) = sizeacol(i-1) + blocksize;
      }
    }
  }
  
  if (sizebcol.n_elem == 1) {
    sizebcol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        sizebcol(i) = sizebcol(0);
      }
    } else {
      for (int i = 1; i < noyears; i++) {
        sizebcol(i) = sizebcol(i-1) + blocksize;
      }
    }
  }
  
  if (sizeccol.n_elem == 1) {
    sizeccol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        sizeccol(i) = sizeccol(0);
      }
    } else {
      for (int i = 1; i < noyears; i++) {
        sizeccol(i) = sizeccol(i-1) + blocksize;
      }
    }
  }
  
  if (repstracol.n_elem == 1) {
    repstracol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        repstracol(i) = repstracol(0);
      }
    } else {
      for (int i = 1; i < noyears; i++) {
        repstracol(i) = repstracol(i-1) + blocksize;
      }
    }
  }
  
  if (repstrbcol.n_elem == 1) {
    repstrbcol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        repstrbcol(i) = repstrbcol(0);
      }
    } else {
      for (int i = 1; i < noyears; i++) {
        repstrbcol(i) = repstrbcol(i-1) + blocksize;
      }
    }
  }
  
  if (fecacol.n_elem == 1) {
    fecacol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        fecacol(i) = fecacol(0);
      }
    } else {
      for (int i = 1; i < noyears; i++) {
        fecacol(i) = fecacol(i-1) + blocksize;
      }
    }
  }
  
  if (fecbcol.n_elem == 1) {
    fecbcol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        fecbcol(i) = fecbcol(0);
      }
    } else {
      for (int i = 1; i < noyears; i++) {
        fecbcol(i) = fecbcol(i-1) + blocksize;
      }
    }
  }
  
  if (indcovacol.n_elem == 1) {
    indcovacol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        indcovacol(i) = indcovacol(0);
      }
    } else {
      for (int i = 1; i < noyears; i++) {
        indcovacol(i) = indcovacol(i-1) + blocksize;
      }
    }
  }
  
  if (indcovbcol.n_elem == 1) {
    indcovbcol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        indcovbcol(i) = indcovbcol(0);
      }
    } else {
      for (int i = 1; i < noyears; i++) {
        indcovbcol(i) = indcovbcol(i-1) + blocksize;
      }
    }
  }
  
  if (indcovccol.n_elem == 1) {
    indcovccol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        indcovccol(i) = indcovccol(0);
      }
    } else {
      for (int i = 1; i < noyears; i++) {
        indcovccol(i) = indcovccol(i-1) + blocksize;
      }
    }
  }
  
  if (aliveacol.n_elem == 1) {
    aliveacol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        aliveacol(i) = aliveacol(0);
      }
    } else {
      for (int i = 1; i < noyears; i++) {
        aliveacol(i) = aliveacol(i-1) + blocksize;
      }
    }
  }
  
  if (deadacol.n_elem == 1) {
    deadacol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        deadacol(i) = deadacol(0);
      }
    } else {
      for (int i = 1; i < noyears; i++) {
        deadacol(i) = deadacol(i-1) + blocksize;
      }
    }
  }
  
  if (obsacol.n_elem == 1) {
    obsacol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        obsacol(i) = obsacol(0);
      }
    } else {
      for (int i = 1; i < noyears; i++) {
        obsacol(i) = obsacol(i-1) + blocksize;
      }
    }
  }
  
  if (nonobsacol.n_elem == 1) {
    nonobsacol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        nonobsacol(i) = nonobsacol(0);
      }
    } else {
      for (int i = 1; i < noyears; i++) {
        nonobsacol(i) = nonobsacol(i-1) + blocksize;
      }
    }
  }
  
  if (censorcol.n_elem == 1) {
    censorcol.resize(noyears);
    
    if (blocksize == 0 || !censrepeat) {
      for (int i = 1; i < noyears; i++) {
        censorcol(i) = censorcol(0);
      }
    } else if (blocksize > 0) {
      for (int i = 1; i < noyears; i++) {
        censorcol(i) = censorcol(i-1) + blocksize;
      }
    }
  }
  
  if (stagecol.n_elem == 1) {
    stagecol.resize(noyears);
    
    if (blocksize == 0) {
      for (int i = 1; i < noyears; i++) {
        stagecol(i) = stagecol(0);
      }
    } else {
      for (int i = 1; i < noyears; i++) {
        stagecol(i) = stagecol(i-1) + blocksize;
      }
    }
  }
  
  int noindivs = data.nrows();
  
  double livcheck1 {0};
  double livcheck2 {0};
  double livcheck3 {0};
  
  double nonobssub1 {0};
  double nonobssub2 {0};
  double nonobssub3 {0};
  
  double stagesize1 {0};
  double stagesize2 {0};
  double stagesize3 {0};
  
  arma::uvec stagemini;
  arma::uvec stagemaxi;
  arma::uvec stageobs;
  arma::uvec stagerep;
  arma::uvec stagemat;
  arma::uvec cs1;
  arma::uvec cs2;
  arma::uvec cs3;
  arma::uvec cs4;
  int choicestage;
  
  Rcpp::StringVector sfname = stageframe["stage"];
  Rcpp::NumericVector repstat = stageframe["repstatus"];
  Rcpp::NumericVector obsstat = stageframe["obsstatus"];
  Rcpp::NumericVector matstat = stageframe["matstatus"];
  arma::vec indataset = stageframe["indataset"];
  arma::vec sfszmin = stageframe["sizebin_min"];
  arma::vec sfszmax = stageframe["sizebin_max"];
  
  arma::vec repstatarma = repstat;
  arma::vec obsstatarma = obsstat;
  arma::vec matstatarma = matstat;
  arma::vec sfszminarma = sfszmin;
  arma::vec sfszmaxarma = sfszmax;
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
  
  int inplace {0};
  for (int i = 0; i < stagenum; i++) {
    stageid(i) = i+1;
    
    if (indataset(i) == 1) {
      instageid(inplace) = i + 1;
      insfszminarma(inplace) = sfszminarma(i);
      insfszmaxarma(inplace) = sfszmaxarma(i);
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
  nonobsgiven1x.fill(-1);
  nonobsgiven2x.fill(-1);
  nonobsgiven3x.fill(-1);
  
  Rcpp::NumericVector zerovec (noindivs);
  Rcpp::NumericVector onevec (noindivs);
  Rcpp::NumericVector negonevec (noindivs);
  zerovec.fill(0);
  onevec.fill(1);
  negonevec.fill(-1);
  
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
  
  Rcpp::NumericVector firstseenx (noindivs);
  Rcpp::NumericVector lastseenx (noindivs);
  Rcpp::NumericVector firstseen (ndflength);
  Rcpp::NumericVector lastseen (ndflength);
  Rcpp::NumericVector obsage (ndflength);
  Rcpp::NumericVector obslifespan (ndflength);
  firstseenx.fill(0);
  lastseenx.fill(0);
  
  Rcpp::NumericVector temp1 (noindivs);
  Rcpp::NumericVector temp2 (noindivs);
  Rcpp::NumericVector temp3 (noindivs);
  Rcpp::NumericVector temp4 (noindivs);
  Rcpp::NumericVector temp5 (noindivs);
  
  for (int i = 0; i < noindivs; i++) {
    firstseenx[i] = 0;
    lastseenx[i] = 0;
    
    for (int k = 0; k < noyears; k++) {
      temp1 = data[sizeacol(k)];
      if (sizebcol(k) != -1 ) {temp2 = data[sizebcol(k)];} else {temp2 = zerovec;}
      if (sizeccol(k) != -1 ) {temp3 = data[sizeccol(k)];} else {temp3 = zerovec;}
      if (repstracol(k) != -1 ) {temp4 = data[repstracol(k)];} else {temp4 = zerovec;}
      if (repstrbcol(k) != -1 ) {temp5 = data[repstrbcol(k)];} else {temp5 = zerovec;}
      
      if (NumericVector::is_na(temp1[i])) {temp1[i] = 0;}
      if (NumericVector::is_na(temp2[i])) {temp2[i] = 0;}
      if (NumericVector::is_na(temp3[i])) {temp3[i] = 0;}
      if (NumericVector::is_na(temp4[i])) {temp4[i] = 0;}
      if (NumericVector::is_na(temp5[i])) {temp5[i] = 0;}
      
      if (temp1[i] > 0 && firstseenx[i] == 0) {
        firstseenx[i] = firstyear + k;
      } else if (temp2[i] > 0 && firstseenx[i] == 0) {
        firstseenx[i] = firstyear + k;
      } else if (temp3[i] > 0 && firstseenx[i] == 0) {
        firstseenx[i] = firstyear + k;
      } else if (temp4[i] > 0 && firstseenx[i] == 0) {
        firstseenx[i] = firstyear + k;
      } else if (temp5[i] > 0 && firstseenx[i] == 0) {
        firstseenx[i] = firstyear + k;
      }
      
      if (temp1[i] > 0 && firstseenx[i] != 0) {
        lastseenx[i] = firstyear + k;
      } else if (temp2[i] > 0 && firstseenx[i] != 0) {
        lastseenx[i] = firstyear + k;
      } else if (temp3[i] > 0 && firstseenx[i] != 0) {
        lastseenx[i] = firstyear + k;
      } else if (temp4[i] > 0 && firstseenx[i] != 0) {
        lastseenx[i] = firstyear + k;
      } else if (temp5[i] > 0 && firstseenx[i] != 0) {
        lastseenx[i] = firstyear + k;
      }
    }
  }

  // Set up main loop
  for (int j = 0; j < (noyears - 1); j++) {
    if (popidcol > -1) popidx = data[popidcol];
    if (patchidcol > -1) patchidx = data[patchidcol];
    if (individcol > -1) individx = data[individcol];
    
    if (j == 0) { // This is implemented in the first year
      
      sizea1x = zerovec;
      sizea2x = data[sizeacol(0)];
      sizea3x = data[sizeacol(1)];
      
      if (stagecol(0) > -1) {
        stage1x.fill("NotAlive");
        stage2x = data[stagecol(0)];
        stage3x = data[stagecol(1)];
      }
      
      if (xcol(0) > -1) {
        xpos1x = zerovec;
        xpos2x = data[xcol(0)];
        xpos3x = data[xcol(1)];
      }
      
      if (ycol(0) > -1) {
        ypos1x = zerovec;
        ypos2x = data[ycol(0)];
        ypos3x = data[ycol(1)];
      }
      
      if (repstracol(0) > -1) {
        repstra1x = zerovec;
        repstra2x = data[repstracol(0)];
        repstra3x = data[repstracol(1)];
      }
      
      if (fecacol(0) > -1) {
        feca1x = zerovec;
        feca2x = data[fecacol(0)];
        feca3x = data[fecacol(1)];
      }
      
      if (sizebcol(0) > -1) {
        sizeb1x = zerovec;
        sizeb2x = data[sizebcol(0)];
        sizeb3x = data[sizebcol(1)];
      }
      
      if (repstrbcol(0) > -1) {
        repstrb1x = zerovec;
        repstrb2x = data[repstrbcol(0)];
        repstrb3x = data[repstrbcol(1)];
      }
      
      if (fecbcol(0) > -1) {
        fecb1x = zerovec;
        fecb2x = data[fecbcol(0)];
        fecb3x = data[fecbcol(1)];
      }
      
      if (sizeccol(0) > -1) {
        sizec1x = zerovec;
        sizec2x = data[sizeccol(0)];
        sizec3x = data[sizeccol(1)];
      }
      
      if (indcovacol(0) > -1) {
        indcova1x = zerovec;
        indcova2x = data[indcovacol(0)];
        indcova3x = data[indcovacol(1)];
      }
      
      if (indcovbcol(0) > -1) {
        indcovb1x = zerovec;
        indcovb2x = data[indcovbcol(0)];
        indcovb3x = data[indcovbcol(1)];
      }
      
      if (indcovccol(0) > -1) {
        indcovc1x = zerovec;
        indcovc2x = data[indcovccol(0)];
        indcovc3x = data[indcovccol(1)];
      }
      
      // This section decides how to deal with censor variables
      if (censorcol(0) > -1) {
        
        if (censrepeat) {
        
          if (censorkeep == 0) {
            censor1x = onevec;
          } else {
            censor1x = zerovec;
          }
          
          censor2x = data[censorcol(0)];
          censor3x = data[censorcol(1)];
          
        } else {
        
          if (censorkeep == 0) {
            censor1x = onevec;
          } else {
            censor1x = zerovec;
          }
          
          censor2x = data[censorcol(0)];
          censor3x = data[censorcol(0)];
        }
      }
      
      
      if (aliveacol(0) > -1) {
        alivegiven1x = zerovec;
        alivegiven2x = data[aliveacol(0)];
        alivegiven3x = data[aliveacol(1)];
      }
      
      if (deadacol(0) > -1) {
        deadgiven1x = zerovec;
        deadgiven2x = data[deadacol(0)];
        deadgiven3x = data[deadacol(1)];
      }
      
      if (obsacol(0) > -1) {
        obsgiven1x = zerovec;
        obsgiven2x = data[obsacol(0)];
        obsgiven3x = data[obsacol(1)];
      }
      
      if (nonobsacol(0) > -1) {
        nonobsgiven1x = negonevec;
        nonobsgiven2x = data[nonobsacol(0)];
        nonobsgiven3x = data[nonobsacol(1)];
      }
      
      if (juvcol(0) > -1) {
        juvgiven1x = zerovec;
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
      livcheck1 = 0;
      livcheck2 = 0;
      livcheck3 = 0;
      
      nonobssub1 = -1;
      nonobssub2 = -1;
      nonobssub3 = -1;
      
      rowid[(i + (j * noindivs))] = i + 1;
      popid[(i + (j * noindivs))] = popidx[i];
      patchid[(i + (j * noindivs))] = patchidx[i];
      individ[(i + (j * noindivs))] = individx[i];
      year2[(i + (j * noindivs))] = firstyear + j;
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
        
        stage1num[(i + (j * noindivs))] = 0;
        stage2num[(i + (j * noindivs))] = 0;
        stage3num[(i + (j * noindivs))] = 0;
      }
      
      
      if (censbool) { // Here we develop the censoring variable
        
        if (NumericVector::is_na(censor1x[i])) {
          censor1[(i + (j * noindivs))] = 0;
        } else if (j != 0) {
          censor1[(i + (j * noindivs))] = 1;
        }
        
        if (NumericVector::is_na(censor2x[i])) {
          censor2[(i + (j * noindivs))] = 0;
        } else {
          censor2[(i + (j * noindivs))] = 1;
        }
        
        if (NumericVector::is_na(censor3x[i])) {
          censor3[(i + (j * noindivs))] = 0;
        } else {
          censor3[(i + (j * noindivs))] = 1;
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
      
      if (NumericVector::is_na(sizea10[(i + (j * noindivs))])) sizea10[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(sizea20[(i + (j * noindivs))])) sizea20[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(sizea30[(i + (j * noindivs))])) sizea30[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(sizeb10[(i + (j * noindivs))])) sizeb10[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(sizeb20[(i + (j * noindivs))])) sizeb20[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(sizeb30[(i + (j * noindivs))])) sizeb30[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(sizec10[(i + (j * noindivs))])) sizec10[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(sizec20[(i + (j * noindivs))])) sizec20[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(sizec30[(i + (j * noindivs))])) sizec30[(i + (j * noindivs))] = 0;
      
      if (NumericVector::is_na(repstra10[(i + (j * noindivs))])) repstra10[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(repstra20[(i + (j * noindivs))])) repstra20[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(repstra30[(i + (j * noindivs))])) repstra30[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(repstrb10[(i + (j * noindivs))])) repstrb10[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(repstrb20[(i + (j * noindivs))])) repstrb20[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(repstrb30[(i + (j * noindivs))])) repstrb30[(i + (j * noindivs))] = 0;
      
      if (NumericVector::is_na(feca10[(i + (j * noindivs))])) feca10[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(feca20[(i + (j * noindivs))])) feca20[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(feca30[(i + (j * noindivs))])) feca30[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(fecb10[(i + (j * noindivs))])) fecb10[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(fecb20[(i + (j * noindivs))])) fecb20[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(fecb30[(i + (j * noindivs))])) fecb30[(i + (j * noindivs))] = 0;
      
      if (NumericVector::is_na(indcova1[(i + (j * noindivs))])) indcova1[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(indcova2[(i + (j * noindivs))])) indcova2[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(indcova3[(i + (j * noindivs))])) indcova3[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(indcovb1[(i + (j * noindivs))])) indcovb1[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(indcovb2[(i + (j * noindivs))])) indcovb2[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(indcovb3[(i + (j * noindivs))])) indcovb3[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(indcovc1[(i + (j * noindivs))])) indcovc1[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(indcovc2[(i + (j * noindivs))])) indcovc2[(i + (j * noindivs))] = 0;
      if (NumericVector::is_na(indcovc3[(i + (j * noindivs))])) indcovc3[(i + (j * noindivs))] = 0;
      
      if (NumericVector::is_na(juvgiven1[(i + (j * noindivs))])) {
        juvgiven1[(i + (j * noindivs))] = 0;
      } else if (juvgiven1[(i + (j * noindivs))] != 0) {
        juvgiven1[(i + (j * noindivs))] = 1;
      }
      
      if (NumericVector::is_na(juvgiven2[(i + (j * noindivs))])) {
        juvgiven2[(i + (j * noindivs))] = 0;
      } else if (juvgiven2[(i + (j * noindivs))] != 0) {
        juvgiven2[(i + (j * noindivs))] = 1;
      }
      
      if (NumericVector::is_na(juvgiven3[(i + (j * noindivs))])) {
        juvgiven3[(i + (j * noindivs))] = 0;
      } else if (juvgiven3[(i + (j * noindivs))] != 0) {
        juvgiven3[(i + (j * noindivs))] = 1;
      }
      
      matstat1[(i + (j * noindivs))] = 1 - juvgiven1[(i + (j * noindivs))];
      matstat2[(i + (j * noindivs))] = 1 - juvgiven2[(i + (j * noindivs))];
      matstat3[(i + (j * noindivs))] = 1 - juvgiven3[(i + (j * noindivs))];
      
      addedsize1[(i + (j * noindivs))] = sizea10[(i + (j * noindivs))] + sizeb10[(i + (j * noindivs))] + sizec10[(i + (j * noindivs))];
      addedsize2[(i + (j * noindivs))] = sizea20[(i + (j * noindivs))] + sizeb20[(i + (j * noindivs))] + sizec20[(i + (j * noindivs))];
      addedsize3[(i + (j * noindivs))] = sizea30[(i + (j * noindivs))] + sizeb30[(i + (j * noindivs))] + sizec30[(i + (j * noindivs))];
      
      addedrepstr1[(i + (j * noindivs))] = repstra10[(i + (j * noindivs))] + (repstrb10[(i + (j * noindivs))] * repstrrel);
      addedrepstr2[(i + (j * noindivs))] = repstra20[(i + (j * noindivs))] + (repstrb20[(i + (j * noindivs))] * repstrrel);
      addedrepstr3[(i + (j * noindivs))] = repstra30[(i + (j * noindivs))] + (repstrb30[(i + (j * noindivs))] * repstrrel);
      
      addedfec1[(i + (j * noindivs))] = feca10[(i + (j * noindivs))] + (fecb10[(i + (j * noindivs))] * fecrel);
      addedfec2[(i + (j * noindivs))] = feca20[(i + (j * noindivs))] + (fecb20[(i + (j * noindivs))] * fecrel);
      addedfec3[(i + (j * noindivs))] = feca30[(i + (j * noindivs))] + (fecb30[(i + (j * noindivs))] * fecrel);
      
      if (NumericVector::is_na(nonobsgiven1[(i + (j * noindivs))])) {
        nonobssub1 = -1;
      } else if (nonobsgiven1[(i + (j * noindivs))] >= 0) {
        nonobsgiven1[(i + (j * noindivs))] = 1;
        nonobssub1 = 1;
      }
      
      if (NumericVector::is_na(nonobsgiven2[(i + (j * noindivs))])) {
        nonobssub2 = -1;
      } else if (nonobsgiven2[(i + (j * noindivs))] >= 0) {
        nonobsgiven2[(i + (j * noindivs))] = 1;
        nonobssub2 = 1;
      }
      
      if (NumericVector::is_na(nonobsgiven3[(i + (j * noindivs))])) {
        nonobssub3 = -1;
      } else if (nonobsgiven3[(i + (j * noindivs))] >= 0) {
        nonobsgiven3[(i + (j * noindivs))] = 1;
        nonobssub3 = 1;
      }
      
      // Observation status
      if (addedsize1[(i + (j * noindivs))] > 0 || obsgiven1[(i + (j * noindivs))] > 0 || nonobssub1 == 0) {
        spryn1[(i + (j * noindivs))] = 1;
      }
      if (addedsize2[(i + (j * noindivs))] > 0 || obsgiven2[(i + (j * noindivs))] > 0 || nonobssub2 == 0) {
        spryn2[(i + (j * noindivs))] = 1;
      } 
      if (addedsize3[(i + (j * noindivs))] > 0 || obsgiven3[(i + (j * noindivs))] > 0 || nonobssub3 == 0) {
        spryn3[(i + (j * noindivs))] = 1;
      }
      
      // Here is thre correction is size 0 individuals that reproduce are actually observed
      if (RepasObs == 1 && addedsize1[(i + (j * noindivs))] == 0 && addedrepstr1[(i + (j * noindivs))] != 0) {
        spryn1[(i + (j * noindivs))] = 1;
      }
      if (RepasObs == 1 && addedsize2[(i + (j * noindivs))] == 0 && addedrepstr2[(i + (j * noindivs))] != 0) {
        spryn2[(i + (j * noindivs))] = 1;
      }
      if (RepasObs == 1 && addedsize3[(i + (j * noindivs))] == 0 && addedrepstr3[(i + (j * noindivs))] != 0) {
        spryn3[(i + (j * noindivs))] = 1;
      }
      
      if (addedrepstr1[(i + (j * noindivs))] > 0) repyn1[(i + (j * noindivs))] = 1;
      if (addedrepstr2[(i + (j * noindivs))] > 0) repyn2[(i + (j * noindivs))] = 1;
      if (addedrepstr3[(i + (j * noindivs))] > 0) repyn3[(i + (j * noindivs))] = 1;
      
      if (addedfec1[(i + (j * noindivs))] > 0) fecyn1[(i + (j * noindivs))] = 1;
      if (addedfec2[(i + (j * noindivs))] > 0) fecyn2[(i + (j * noindivs))] = 1;
      if (addedfec3[(i + (j * noindivs))] > 0) fecyn3[(i + (j * noindivs))] = 1;
      
      if (nonobssub1 >= 0) {
        livcheck1 = addedsize1[(i + (j * noindivs))] + addedrepstr1[(i + (j * noindivs))] + spryn1[(i + (j * noindivs))] + nonobssub1;
      } else {
        livcheck1 = addedsize1[(i + (j * noindivs))] + addedrepstr1[(i + (j * noindivs))] + spryn1[(i + (j * noindivs))];
      }
      if (nonobssub2 >= 0) {
        livcheck2 = addedsize2[(i + (j * noindivs))] + addedrepstr2[(i + (j * noindivs))] + spryn2[(i + (j * noindivs))] + nonobssub2;
      } else {
        livcheck2 = addedsize2[(i + (j * noindivs))] + addedrepstr2[(i + (j * noindivs))] + spryn2[(i + (j * noindivs))];
      }
      if (nonobssub3 >= 0) {
        livcheck3 = addedsize3[(i + (j * noindivs))] + addedrepstr3[(i + (j * noindivs))] + spryn3[(i + (j * noindivs))] + nonobssub3;
      } else {
        livcheck3 = addedsize3[(i + (j * noindivs))] + addedrepstr3[(i + (j * noindivs))] + spryn3[(i + (j * noindivs))];
      }
      
      if (livcheck1 > 0) alive1[(i + (j * noindivs))] = 1;
      if (livcheck2 > 0) alive2[(i + (j * noindivs))] = 1;
      if (livcheck3 > 0) alive3[(i + (j * noindivs))] = 1;
      
      if (alivegiven1[(i + (j * noindivs))] > 0) alive1[(i + (j * noindivs))] = 1;
      if (alivegiven2[(i + (j * noindivs))] > 0) alive2[(i + (j * noindivs))] = 1;
      if (alivegiven3[(i + (j * noindivs))] > 0) alive3[(i + (j * noindivs))] = 1;
      
      if (deadgiven1[(i + (j * noindivs))] > 0) alive1[(i + (j * noindivs))] = 0;
      if (deadgiven2[(i + (j * noindivs))] > 0) alive2[(i + (j * noindivs))] = 0;
      if (deadgiven3[(i + (j * noindivs))] > 0) alive3[(i + (j * noindivs))] = 0;
      
      // Corrections to living status, lifespan, and age based on later sightings
      if (firstseenx[i] <= (j + firstyear) && lastseenx[i] >= (j + firstyear)) { // This loop needs to be altered to check all years for lastseen adjustments, once firstseen is determined
        firstseen[(i + (j * noindivs))] = firstseenx[i];
        lastseen[(i + (j * noindivs))] = lastseenx[i];
        
        obsage[(i + (j * noindivs))] = year2[(i + (j * noindivs))] - firstseen[(i + (j * noindivs))];
        obslifespan[(i + (j * noindivs))] = lastseen[(i + (j * noindivs))] - firstseen[(i + (j * noindivs))];
        
        alive2[(i + (j * noindivs))] = 1;
        
        if (lastseenx[i] >= (j + firstyear + 1)) {alive3[(i + (j * noindivs))] = 1;}
        if (firstseenx[i] <= (j + firstyear - 1)) {alive1[(i + (j * noindivs))] = 1;}
      }
      
      // Here we take care of stage assignments
      if (stassign && stagecol(0) == -1) {
        
        if (stszcol == 4) {
          
          stagesize1 = addedsize1[(i + (j * noindivs))];
          stagesize2 = addedsize2[(i + (j * noindivs))];
          stagesize3 = addedsize3[(i + (j * noindivs))];
          
        } else if (stszcol == 3) {
          
          stagesize1 = sizec10[(i + (j * noindivs))];
          stagesize2 = sizec20[(i + (j * noindivs))];
          stagesize3 = sizec30[(i + (j * noindivs))];
          
        } else if (stszcol == 2) {
          
          stagesize1 = sizeb10[(i + (j * noindivs))];
          stagesize2 = sizeb20[(i + (j * noindivs))];
          stagesize3 = sizeb30[(i + (j * noindivs))];
          
        } else {
          
          stagesize1 = sizea10[(i + (j * noindivs))];
          stagesize2 = sizea20[(i + (j * noindivs))];
          stagesize3 = sizea30[(i + (j * noindivs))];
          
        }
        
        // Stage 2
        stagemini = find(insfszminarma < stagesize2);
        stagemaxi = find(insfszmaxarma >= stagesize2);
        stagerep = find(inrepstatarma == repyn2[(i + (j * noindivs))]);
        stagemat = find(inmatstatarma == matstat2[(i + (j * noindivs))]);
        stageobs = find(inobsstatarma == spryn2[(i + (j * noindivs))]);
        
        cs1 = intersect(stagemini, stagemaxi);
        cs2 = intersect(stageobs, stagemat);
        
        if (NRasRep) {
          cs4 = intersect(cs1, cs2);
        } else {
          cs3 = intersect(stagerep, cs2);
          cs4 = intersect(cs1, cs3);
        }
        
        if (cs4.n_elem > 0 && alive2[(i + (j * noindivs))] == 1) {
          choicestage = instageid(cs4[0]) - 1;
          stage2num[(i + (j * noindivs))] = choicestage + 1;
          
          stage2[(i + (j * noindivs))] = sfname[choicestage];
        } else {
          stage2[(i + (j * noindivs))] = "NotAlive";
        }
        
        // Stage 1
        stagemini = find(insfszminarma < stagesize1);
        stagemaxi = find(insfszmaxarma >= stagesize1);
        stagerep = find(inrepstatarma == repyn1[(i + (j * noindivs))]);
        stagemat = find(inmatstatarma == matstat1[(i + (j * noindivs))]);
        stageobs = find(inobsstatarma == spryn1[(i + (j * noindivs))]);
        
        cs1 = intersect(stagemini, stagemaxi);
        cs2 = intersect(stageobs, stagemat);
        
        if (NRasRep) {
          cs4 = intersect(cs1, cs2);
        } else {
          cs3 = intersect(stagerep, cs2);
          cs4 = intersect(cs1, cs3);
        }
        
        if (cs4.n_elem > 0 && alive1[(i + (j * noindivs))] == 1) {
          choicestage = instageid(cs4[0]) - 1;
          stage1num[(i + (j * noindivs))] = choicestage + 1;
          
          stage1[(i + (j * noindivs))] = sfname[choicestage];
        } else {
          stage1[(i + (j * noindivs))] = "NotAlive";
          matstat1[(i + (j * noindivs))] = 0;
        }
        
        // Stage 3
        stagemini = find(insfszminarma < stagesize3);
        stagemaxi = find(insfszmaxarma >= stagesize3);
        stagerep = find(inrepstatarma == repyn3[(i + (j * noindivs))]);
        stagemat = find(inmatstatarma == matstat3[(i + (j * noindivs))]);
        stageobs = find(inobsstatarma == spryn3[(i + (j * noindivs))]);
        
        cs1 = intersect(stagemini, stagemaxi);
        cs2 = intersect(stageobs, stagemat);
        
        if (NRasRep) {
          cs4 = intersect(cs1, cs2);
        } else {
          cs3 = intersect(stagerep, cs2);
          cs4 = intersect(cs1, cs3);
        }
        
        // Here we create exceptions based on stage assignment problems in time t+1
        if (cs4.n_elem == 1 && alive3[(i + (j * noindivs))] == 1) {
          
          choicestage = instageid(cs4[0]) - 1;
          stage3num[(i + (j * noindivs))] = choicestage + 1;
          
          stage3[(i + (j * noindivs))] = sfname[choicestage];
          
        } else if (alive3[(i + (j * noindivs))] != 1) {
          
          stage3[(i + (j * noindivs))] = "NotAlive";
          
        } else if (cs4.n_elem == 0) {
          
          stage3[(i + (j * noindivs))] = "NoMatch";
          
          Rf_warningcall(R_NilValue, "Some stages occurring in the dataset do not match any characteristics in the input stageframe.");
          
        } else if (cs4.n_elem > 1) {
          
          Rf_warningcall(R_NilValue, "Some stages in the input stageframe appear to have the same description. Please make sure that all stages included in the stageframe are defined with unique sets of characteristics.");
          
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
  
  Rcpp::DataFrame df0 = DataFrame::create(Named("rowid") = rowid,
    _["popid"] = popid, _["patchid"] = patchid, _["individ"] = individ,
    _["year2"] = year2, _["firstseen"] = firstseen, _["lastseen"] = lastseen,
    _["obsage"] = obsage, _["obslifespan"] = obslifespan);
  
  Rcpp::DataFrame df1 = DataFrame::create(Named("xpos1") = xpos1,
    _["ypos1"] = ypos1, _["sizea1"] = sizea1, _["sizeb1"] = sizeb1,
    _["sizec1"] = sizec1, _["sizeadded1"] = addedsize1,
    _["repstra1"] = repstra1, _["repstrb1"] = repstrb1,
    _["repstradded1"] = addedrepstr1, _["feca1"] = feca1, _["fecb1"] = fecb1,
    _["fecadded1"] = addedfec1, _["indcova1"] = indcova1,
    _["indcovb1"] = indcovb1, _["indcovc1"] = indcovc1, _["censor1"] = censor1,
    _["juvgiven1"] = juvgiven1);
  
  Rcpp::DataFrame df1a = DataFrame::create(Named("obsstatus1") = spryn1,
    _["repstatus1"] = repyn1, _["fecstatus1"] = fecyn1,
    _["matstatus1"] = matstat1, _["alive1"] = alive1, _["stage1"] = stage1,
    _["stage1index"] = stage1num);
  
  Rcpp::DataFrame df2 = DataFrame::create(Named("xpos2") = xpos2,
    _["ypos2"] = ypos2,  _["sizea2"] = sizea2, _["sizeb2"] = sizeb2,
    _["sizec2"] = sizec2, _["sizeadded2"] = addedsize2,
    _["repstra2"] = repstra2, _["repstrb2"] = repstrb2,
    _["repstradded2"] = addedrepstr2, _["feca2"] = feca2, _["fecb2"] = fecb2,
    _["fecadded2"] = addedfec2, _["indcova2"] = indcova2,
    _["indcovb2"] = indcovb2, _["indcovc2"] = indcovc2, _["censor2"] = censor2,
    _["juvgiven2"] = juvgiven2);
  
  Rcpp::DataFrame df2a = DataFrame::create(Named("obsstatus2") = spryn2,
    _["repstatus2"] = repyn2, _["fecstatus2"] = fecyn2,
    _["matstatus2"] = matstat2, _["alive2"] = alive2, _["stage2"] = stage2,
    _["stage2index"] = stage2num);
  
  Rcpp::DataFrame df3 = DataFrame::create(Named("xpos3") = xpos3,
    _["ypos3"] = ypos3, _["sizea3"] = sizea3, _["sizeb3"] = sizeb3,
    _["sizec3"] = sizec3, _["sizeadded3"] = addedsize3,
    _["repstra3"] = repstra3, _["repstrb3"] = repstrb3,
    _["repstradded3"] = addedrepstr3, _["feca3"] = feca3, _["fecb3"] = fecb3,
    _["fecadded3"] = addedfec3, _["indcova3"] = indcova3,
    _["indcovb3"] = indcovb3, _["indcovc3"] = indcovc3, _["censor3"] = censor3,
    _["juvgiven3"] = juvgiven3);
  
  Rcpp::DataFrame df3a = DataFrame::create(Named("obsstatus3") = spryn3,
    _["repstatus3"] = repyn3, _["fecstatus3"] = fecyn3,
    _["matstatus3"] = matstat3, _["alive3"] = alive3, _["stage3"] = stage3,
    _["stage3index"] = stage3num);
  
  
  return Rcpp::List::create(Named("a") = df0, Named("b") = df1,
    Named("c") = df1a, Named("d") = df2, Named("e") = df2a, Named("f") = df3,
    Named("g") = df3a);
}

//' Create Historical Vertical Structure for Ahistorical Vertical Data Frame
//' 
//' Function \code{jpf()} powers the R function \code{\link{historicalize3}()},
//' creating the historical, vertical structure and rearranging the data in that
//' shape.
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
//' @param year2col Column number of year or time step in time \emph{t}.
//' @param year3col Column number of year or time step in time \emph{t}+1.
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
//' the production of reproductive structures, such as flowers, in time \emph{t}.
//' @param repstra3col Column number corresponding to the main variable coding
//' the production of reproductive structures, such as flowers, in time
//' \emph{t}+1.
//' @param repstrb2col Column number corresponding to a second variable coding
//' the production of reproductive structures, such as flowers, in time \emph{t}.
//' @param repstrb3col Column number corresponding to a second variable coding
//' the production of reproductive structures, such as flowers, in time
//' \emph{t}+1.
//' @param feca2col Column number corresponding to the main variable coding for
//' fecundity in time \emph{t}.
//' @param feca3col Column number corresponding to the main variable coding for
//' fecundity in time \emph{t}+1.
//' @param fecb2col Column number corresponding to a second variable coding for
//' fecundity in time \emph{t}.
//' @param fecb3col Column number corresponding to a second variable coding for
//' fecundity in time \emph{t}+1.
//' @param indcova2col Column number corresponding to an individual covariate
//' in time \emph{t}.
//' @param indcova3col Column number corresponding to an individual covariate
//' in time \emph{t}+1.
//' @param indcovb2col Column number corresponding to an individual covariate
//' in time \emph{t}.
//' @param indcovb3col Column number corresponding to an individual covariate
//' in time \emph{t}+1.
//' @param indcovc2col Column number corresponding to an individual covariate
//' in time \emph{t}.
//' @param indcovc3col Column number corresponding to an individual covariate
//' in time \emph{t}+1.
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
//' unobservable stage in time \emph{t}.
//' @param nonobs2col Column number detailing whether an individual is in an
//' unobservable stage in time \emph{t}.
//' @param nonobs3col Column number detailing whether an individual is in an
//' unobservable stage in time \emph{t}+1.
//' @param repstrrel This is a scalar modifier for that makes the variable in
//' \code{repstrb2col} equivalent to \code{repstra2col}.
//' @param fecrel This is a scalar modifier for that makes the variable in
//' \code{fecb2col} equivalent to \code{feca2col}.
//' @param stage2col Column number corresponding to life history stage in time
//' \emph{t}.
//' @param stage3col Column number corresponding to life history stage in time
//' \emph{t}+1.
//' @param censorcol Column number corresponding to a censor variable within the
//' dataset.
//' @param NAas0 If TRUE, then all NA entries for size and fecundity variables
//' will be set to 0.
//' @param NRasRep If TRUE, then will treat non-reproductive but mature
//' individuals as reproductive during stage assignment.
//' @param stassign A logical value indicating whether to assign stages.
//' @param stszcol Column number describing which size variable to use in stage 
//' estimation.
//' @param censorkeep Numeric value of censor variable, denoting elements to
//' keep. If NA is to be used, then set this variable to 0 and set
//' \code{censbool = TRUE}.
//' @param censbool A logical variable determining whether NA denotes the value
//' of the censoring variable identifying data to keep.
//' 
//' @return The output is currently a 7 element list, where each element is a
//' data frame with the same number of rows.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
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
  bool stassign, int stszcol, double censorkeep, bool censbool) {
  
  int norows = data.nrows(); // The number of data points in the demographic dataset
  
  Rcpp::NumericVector ckcheck (1); // This section through 1378 is new
  ckcheck(0) = censorkeep;
  double crazycensor;
  Rcpp::NumericVector censfillvec (norows);
  if (NumericVector::is_na(ckcheck(0))) {
    crazycensor = 0;
  } else{
    crazycensor = censorkeep;
  }                               // End new section
  
  Rcpp::NumericVector zerovec (norows);
  Rcpp::NumericVector negonevec (norows);
  zerovec.fill(0);
  negonevec.fill(-1);
  
  Rcpp::StringVector individx(norows);
  individx = data[individcol];
  Rcpp::StringVector allindivs = unique(individx);
  int noindivs = allindivs.size(); // This is the total number of individuals in the dataset
  
  Rcpp::IntegerVector year2x(norows);
  Rcpp::IntegerVector year3x(norows);
  year2x = data[year2col];
  if (year3col != -1) {year3x = data[year3col];} else {year3x = zerovec;}
  Rcpp::IntegerVector yearall2x = sort_unique(year2x);
  int firstyear = min(yearall2x);
  int noyears = yearall2x.size(); // This is the total number of observation periods
  
  int ndflength = noyears * noindivs; // This is the initial length of the final vertical dataset with all year x indiv combos
  int currentyear {0};
  int currentindiv {-1};
  int ndfindex {0};
  int prevyrindex {0};
  int nextyrindex {0};
  double livcheck2 {0};
  double livcheck3 {0};
  
  Rcpp::StringVector sfname = stageframe["stage"]; // This section reads in the stageframe
  Rcpp::NumericVector repstat = stageframe["repstatus"];
  Rcpp::NumericVector obsstat = stageframe["obsstatus"];
  Rcpp::NumericVector matstat = stageframe["matstatus"];
  arma::vec indataset = stageframe["indataset"];
  arma::vec sfszmin = stageframe["sizebin_min"];
  arma::vec sfszmax = stageframe["sizebin_max"];
  
  arma::vec repstatarma = repstat;
  arma::vec obsstatarma = obsstat;
  arma::vec matstatarma = matstat;
  arma::vec sfszminarma = sfszmin;
  arma::vec sfszmaxarma = sfszmax;
  int stagenum = sfszmaxarma.n_elem; // This is the total number of life history stages in the stageframe
  
  arma::uvec instages = find(indataset == 1); 
  int instagenum = instages.n_elem; // This is the total number of life history stages in the demographic dataset
  
  arma::uvec stageid(stagenum);
  arma::uvec instageid(instagenum);
  arma::vec insfszminarma(instagenum);
  arma::vec insfszmaxarma(instagenum);
  arma::vec inrepstatarma(instagenum);
  arma::vec inobsstatarma(instagenum);
  arma::vec inmatstatarma(instagenum);
  
  int inplace {0}; // This section creates vectors describing only life history stages that actually occur in the dataset
  for (int i = 0; i < stagenum; i++) {
    stageid(i) = i+1;
    
    if (indataset(i) == 1) {
      instageid(inplace) = i + 1;
      insfszminarma(inplace) = sfszminarma(i);
      insfszmaxarma(inplace) = sfszmaxarma(i);
      inrepstatarma(inplace) = repstatarma(i);
      inobsstatarma(inplace) = obsstatarma(i);
      inmatstatarma(inplace) = matstatarma(i);
      
      inplace++;
    }
  }
  
  // The following set of variable definitions is different from those used in pfj.
  // Here, these variables are defined by the row structure of the data, whereas in
  // pfj they correspond to the individuals in the data frame.
  Rcpp::StringVector popidx (norows);
  Rcpp::StringVector patchidx (norows);
  Rcpp::NumericVector xpos2x (norows);
  Rcpp::NumericVector ypos2x (norows);
  Rcpp::NumericVector xpos3x (norows);
  Rcpp::NumericVector ypos3x (norows);
  Rcpp::NumericVector sizea2x (norows);
  Rcpp::NumericVector sizea3x (norows);
  Rcpp::NumericVector repstra2x (norows);
  Rcpp::NumericVector repstra3x (norows);
  Rcpp::NumericVector feca2x (norows);
  Rcpp::NumericVector feca3x (norows);
  Rcpp::NumericVector sizeb2x (norows);
  Rcpp::NumericVector sizeb3x (norows);
  Rcpp::NumericVector repstrb2x (norows);
  Rcpp::NumericVector repstrb3x (norows);
  Rcpp::NumericVector fecb2x (norows);
  Rcpp::NumericVector fecb3x (norows);
  Rcpp::NumericVector sizec2x (norows);
  Rcpp::NumericVector sizec3x (norows);
  
  Rcpp::NumericVector indcova2x (norows);
  Rcpp::NumericVector indcova3x (norows);
  Rcpp::NumericVector indcovb2x (norows);
  Rcpp::NumericVector indcovb3x (norows);
  Rcpp::NumericVector indcovc2x (norows);
  Rcpp::NumericVector indcovc3x (norows);
  
  Rcpp::NumericVector censor2x (norows);
  Rcpp::NumericVector censor3x (norows);
  Rcpp::NumericVector alivegiven2x (norows);
  Rcpp::NumericVector alivegiven3x (norows);
  Rcpp::NumericVector deadgiven2x (norows);
  Rcpp::NumericVector deadgiven3x (norows);
  Rcpp::NumericVector obsgiven2x (norows);
  Rcpp::NumericVector obsgiven3x (norows);
  Rcpp::NumericVector nonobsgiven2x (norows); 
  Rcpp::NumericVector nonobsgiven3x (norows);
  Rcpp::NumericVector juvgiven2x (norows);
  Rcpp::NumericVector juvgiven3x (norows);
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
  
  Rcpp::StringVector stage2x (norows);
  Rcpp::StringVector stage3x (norows);
  
  // Assign values from the dataset
  if (popidcol != -1) {popidx = data[popidcol];} else {popidx = zerovec;}
  if (patchidcol != -1) {patchidx = data[patchidcol];} else {patchidx = zerovec;}
  if (censorcol != -1) {censor2x = data[censorcol];} else {censor2x = zerovec;}
  
  sizea2x = data[sizea2col];
  if (sizeb2col != -1) {sizeb2x = data[sizeb2col];} else {sizeb2x = zerovec;}
  if (sizec2col != -1) {sizec2x = data[sizec2col];} else {sizec2x = zerovec;}
  if (repstra2col != -1) {repstra2x = data[repstra2col];} else {repstra2x = zerovec;}
  if (repstrb2col != -1) {repstrb2x = data[repstrb2col];} else {repstrb2x = zerovec;}
  if (feca2col != -1) {feca2x = data[feca2col];} else {feca2x = zerovec;}
  if (fecb2col != -1) {fecb2x = data[fecb2col];} else {fecb2x = zerovec;}
  
  if (indcova2col != -1) {indcova2x = data[indcova2col];} else {indcova2x = zerovec;}
  if (indcova3col != -1) {indcova3x = data[indcova3col];} else {indcova3x = zerovec;}
  if (indcovb2col != -1) {indcovb2x = data[indcovb2col];} else {indcovb2x = zerovec;}
  if (indcovb3col != -1) {indcovb3x = data[indcovb3col];} else {indcovb3x = zerovec;}
  if (indcovc2col != -1) {indcovc2x = data[indcovc2col];} else {indcovc2x = zerovec;}
  if (indcovc3col != -1) {indcovc3x = data[indcovc3col];} else {indcovc3x = zerovec;}
  
  if (juv2col != -1) {juvgiven2x = data[juv2col];} else {juvgiven2x = zerovec;}
  if (obs2col != -1) {obsgiven2x = data[obs2col];} else {obsgiven2x.fill(-1);} //Set to -1 to make nonobs vs dead designation easier
  if (nonobs2col != -1) {nonobsgiven2x = data[nonobs2col];} else {nonobsgiven2x.fill(-1);}
  if (alive2col != -1) {alivegiven2x = data[alive2col];} else {alivegiven2x.fill(-1);}
  if (dead2col != -1) {deadgiven2x = data[dead2col];} else {deadgiven2x.fill(-1);}
  if (xcol != -1) {xpos2x = data[xcol];} else {xpos2x.fill(-1);}
  if (ycol != -1) {ypos2x = data[ycol];} else {ypos2x.fill(-1);}
  
  if (sizea3col != -1) {sizea3x = data[sizea3col];} else {sizea3x = zerovec;}
  if (sizeb3col != -1) {sizeb3x = data[sizeb3col];} else {sizeb3x = zerovec;}
  if (sizec3col != -1) {sizec3x = data[sizec3col];} else {sizec3x = zerovec;}
  if (repstra3col != -1) {repstra3x = data[repstra3col];} else {repstra3x = zerovec;}
  if (repstrb3col != -1) {repstrb3x = data[repstrb3col];} else {repstrb3x = zerovec;}
  if (feca3col != -1) {feca3x = data[feca3col];} else {feca3x = zerovec;}
  if (fecb3col != -1) {fecb3x = data[fecb3col];} else {fecb3x = zerovec;}
  if (juv3col != -1) {juvgiven3x = data[juv3col];} else {juvgiven3x = zerovec;}
  if (obs3col != -1) {obsgiven3x = data[obs3col];} else {obsgiven3x.fill(-1);} //Set to -1 to make nonobs vs dead designation easier
  if (nonobs3col != -1) {nonobsgiven3x = data[nonobs3col];} else {nonobsgiven3x.fill(-1);}
  if (alive3col != -1) {alivegiven3x = data[alive3col];} else {alivegiven3x.fill(-1);}
  if (dead3col != -1) {deadgiven3x = data[dead3col];} else {deadgiven3x.fill(-1);}
  
  if (stage2col != -1) {stage2x = data[stage2col];}
  if (stage3col != -1) {stage3x = data[stage3col];}
  
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
  
  arma::ivec year2 (ndflength);   // valgrind mentioned this line, which used to be an Rcpp::IntegerVector
  year2.zeros();
  
  Rcpp::NumericVector indcova2 (ndflength);
  Rcpp::NumericVector indcovb2 (ndflength);
  Rcpp::NumericVector indcovc2 (ndflength);
  
  Rcpp::IntegerVector repstatus2 (ndflength);
  Rcpp::IntegerVector fecstatus2 (ndflength);
  Rcpp::IntegerVector obsstatus2 (ndflength);
  Rcpp::NumericVector juvgiven2 (ndflength);
  Rcpp::IntegerVector matstat2 (ndflength);
  
  censor2.fill(crazycensor); //This variable used to be filled with 0
  
  xpos2.fill(0);
  ypos2.fill(0);
  sizea2.fill(0);
  sizeb2.fill(0);
  sizec2.fill(0);
  sizeadded2.fill(0);
  repstra2.fill(0);
  repstrb2.fill(0);
  repstradded2.fill(0);
  feca2.fill(0);
  fecb2.fill(0);
  fecadded2.fill(0);
  indcova2.fill(0);
  indcovb2.fill(0);
  indcovc2.fill(0);
  repstatus2.fill(0);
  fecstatus2.fill(0);
  obsstatus2.fill(0);
  juvgiven2.fill(0);
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
  
  censor3.fill(crazycensor); // This used to be filled with 0
  
  xpos3.fill(0);
  ypos3.fill(0);
  sizea3.fill(0);
  sizeb3.fill(0);
  sizec3.fill(0);
  sizeadded3.fill(0);
  repstra3.fill(0);
  repstrb3.fill(0);
  repstradded3.fill(0);
  feca3.fill(0);
  fecb3.fill(0);
  fecadded3.fill(0);
  indcova3.fill(0);
  indcovb3.fill(0);
  indcovc3.fill(0);
  repstatus3.fill(0);
  fecstatus3.fill(0);
  obsstatus3.fill(0);
  juvgiven3.fill(0);
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
  
  censor1.fill(crazycensor); // This used to be filled with 0
  
  xpos1.fill(0);
  ypos1.fill(0);
  sizea1.fill(0);
  sizeb1.fill(0);
  sizec1.fill(0);
  sizeadded1.fill(0);
  repstra1.fill(0);
  repstrb1.fill(0);
  repstradded1.fill(0);
  feca1.fill(0);
  fecb1.fill(0);
  fecadded1.fill(0);
  indcova1.fill(0);
  indcovb1.fill(0);
  indcovc1.fill(0);
  repstatus1.fill(0);
  fecstatus1.fill(0);
  obsstatus1.fill(0);
  juvgiven1.fill(0);
  matstat1.fill(0);
  
  // This section introduces variables used to check whether the censor variable has been checked and set at each step
  arma::uvec indivnum (ndflength);          // New Line
  arma::uvec censor2check (ndflength);      // New line
  indivnum.zeros();                         // New line
  censor2check.zeros();                     // New line
  
  // Here we introduce some derived variables that require extra looping or other control parameters
  arma::ivec firstseen (ndflength);
  arma::ivec lastseen (ndflength);
  arma::ivec obsage (ndflength);
  arma::ivec obslifespan (ndflength);
  Rcpp::NumericVector alive1 (ndflength);
  Rcpp::NumericVector alive2 (ndflength);
  Rcpp::NumericVector alive3 (ndflength);
  firstseen.fill(-1);
  lastseen.fill(-1);
  alive1.fill(0);
  alive2.fill(0);
  alive3.fill(0);
  
  Rcpp::StringVector stage1 (ndflength);
  Rcpp::StringVector stage2 (ndflength);
  Rcpp::StringVector stage3 (ndflength);
  Rcpp::NumericVector stage1num (ndflength);
  Rcpp::NumericVector stage2num (ndflength);
  Rcpp::NumericVector stage3num (ndflength);
  stage1num.fill(0);
  stage2num.fill(0);
  stage3num.fill(0);
  
  double stagesize1 {0};
  double stagesize2 {0};
  double stagesize3 {0};
  
  arma::uvec stagemini;
  arma::uvec stagemaxi;
  arma::uvec stageobs;
  arma::uvec stagerep;
  arma::uvec stagemat;
  arma::uvec cs1;
  arma::uvec cs2;
  arma::uvec cs3;
  arma::uvec cs4;
  int choicestage {0};
  
  // Main loop, which creates the main new dataset rows. Establishes state in time t for all cases in which an individual is observed
  for (int i = 0; i < norows; i++) { // Variable i corresponds to row in the old dataset
    for (int j = 0; j < noyears; j++) { // This establishes a place marker for vectors corresponding to the current year
      if (year2x[i] == yearall2x[j]) currentyear = j;
    }
    
    currentindiv = -1;
    for (int k = 0; k < noindivs; k++) { // This establishes a place marker variable corresponding to the current individual
      if (individx[i] == allindivs[k]) {
        currentindiv = k;
        indivnum[i] = k;
      }
    }
    
    ndfindex = (noyears * currentindiv) + currentyear; // This establishes the row in the new dataset being created
    
    if (NumericVector::is_na(sizea2x[i])) {
      sizea20x[i] = 0;
      if (NAas0) {sizea2x[i] = 0;}
    } else {sizea20x[i] = sizea2x[i];}
    if (NumericVector::is_na(sizeb2x[i])) {
      sizeb20x[i] = 0;
      if (NAas0) {sizeb2x[i] = 0;}
    } else {sizeb20x[i] = sizeb2x[i];}
    if (NumericVector::is_na(sizec2x[i])) {
      sizec20x[i] = 0;
      if (NAas0) {sizec2x[i] = 0;}
    } else {sizec20x[i] = sizec2x[i];}
    if (NumericVector::is_na(repstra2x[i])) {
      repstra20x[i] = 0;
      if (NAas0) {repstra2x[i] = 0;}
    } else {repstra20x[i] = repstra2x[i];}
    if (NumericVector::is_na(repstrb2x[i])) {
      repstrb20x[i] = 0;
      if (NAas0) {repstrb2x[i] = 0;}
    } else {repstrb20x[i] = repstrb2x[i];}
    if (NumericVector::is_na(feca2x[i])) {
      feca20x[i] = 0;
      if (NAas0) {feca2x[i] = 0;}
    } else {feca20x[i] = feca2x[i];}
    if (NumericVector::is_na(fecb2x[i])) {
      fecb20x[i] = 0;
      if (NAas0) {fecb2x[i] = 0;}
    } else {fecb20x[i] = fecb2x[i];}
    
    if (NumericVector::is_na(juvgiven2x[i])) {
      juvgiven20x[i] = 0;
    } else if (juvgiven2x[i] != 0) {
      juvgiven20x[i] = 1;
    } else {juvgiven20x[i] = 0;}
    
    // Here we develop the censoring variable
    if (censbool && censorcol != -1) { // This provides a replacement in cases where NA designates data to keep
      if (NumericVector::is_na(censor2x[i])) {
        censor2[ndfindex] = 0;
        censor2check[ndfindex] = 1;
      } else {
        censor2[ndfindex] = 1;
        censor2check[ndfindex] = 1;
      }
    } else if (censorcol != -1) {
      if (censorkeep == 0 && NumericVector::is_na(censor2x[i])) {
        censor2[ndfindex] = 1;
        censor2check[ndfindex] = 1;
      } else if (censorkeep == 1 && NumericVector::is_na(censor2x[i])) {
        censor2[ndfindex] = 0;
        censor2check[ndfindex] = 1;
      } else {
        censor2[ndfindex] = censor2x[i];
        censor2check[ndfindex] = 1;
      }
    }
    
    rowid[ndfindex] = i;
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
    
    if (repstradded2[ndfindex] > 0) {repstatus2[ndfindex] = 1;} else {repstatus2[ndfindex] = 0;}
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
    
    if (alivegiven2x[i] > 0) {alive2[ndfindex] = 1;} else if (alivegiven2x[i] == 0) {alive2[ndfindex] = 0;}
    if (deadgiven2x[i] > 0) {alive2[ndfindex] = 0;} else if (deadgiven2x[i] == 0) {alive2[ndfindex] = 1;}
    
    livcheck2 = sizeadded2[ndfindex] + repstradded2[ndfindex] + obsstatus2[ndfindex];
    
    if (livcheck2 > 0 && firstseenx[currentindiv] == -1) {
      firstseenx[currentindiv] = currentyear + firstyear;
      lastseenx[currentindiv] = currentyear + firstyear;
      alive2[ndfindex] = 1;
    } else if (livcheck2 > 0) {
      lastseenx[currentindiv] = currentyear + firstyear;
      alive2[ndfindex] = 1;
    }
    
    if (alive2[ndfindex] == 1 && matstat2[ndfindex] == 1) {
      matstat3[ndfindex] = 1;
    }
    
    if (stage2col != -1 && alive2[ndfindex] == 1) {
      stage2[ndfindex] = stage2x[i];
    } else if (stassign) {stage2[ndfindex] = "NotAlive";}
    
    
    // Now we work on time t+1 for the last possible time t (which is technically the second to last time), in cases where t+1 columns are provided
    if (currentyear == (noyears - 1)) {
      if (censbool && censorcol != -1) { // Here we develop the censoring variable for the last time
        if (NumericVector::is_na(censor2x[i])) {
          censor3[ndfindex] = 0;
        } else {
          censor3[ndfindex] = 1;
        }
      } else if (censorcol != -1) {     // New section through 1908
        if (censorkeep == 0 && NumericVector::is_na(censor2x[i])) {
          censor3[ndfindex] = 1;
        } else if (censorkeep == 1 && NumericVector::is_na(censor2x[i])) {
          censor3[ndfindex] = 0;
        } else {
          censor3[ndfindex] = censor2x[i];
        }
      }
      
      if (NumericVector::is_na(juvgiven3x[i])) {
        juvgiven30x[i] = 0;
      } else if (juvgiven3x[i] != 0) {
        juvgiven30x[i] = 1;
      } else {juvgiven30x[i] = 0;}
      
      if (sizea3col != -1) {
        if (NumericVector::is_na(sizea3x[i])) {
          sizea30x[i] = 0;
          if (NAas0) {sizea3x[i] = 0;}
        } else {sizea30x[i] = sizea3x[i];}
        sizea3[ndfindex] = sizea3x[i];
      }
      if (sizeb3col != -1) {
        if (NumericVector::is_na(sizeb3x[i])) {
          sizeb30x[i] = 0;
          if (NAas0) {sizeb3x[i] = 0;}
        } else {sizeb30x[i] = sizeb3x[i];}
        sizeb3[ndfindex] = sizeb3x[i];
      }
      if (sizec3col != -1) {
        if (NumericVector::is_na(sizec3x[i])) {
          sizec30x[i] = 0;
          if (NAas0) {sizec3x[i] = 0;}
        } else {sizec30x[i] = sizec3x[i];}
        sizec3[ndfindex] = sizec3x[i];
      }
      sizeadded3[ndfindex] = sizea30x[i] + sizeb30x[i] + sizec30x[i];
      
      if (repstra3col != -1) {
        if (NumericVector::is_na(repstra3x[i])) {
          repstra30x[i] = 0;
          if (NAas0) {repstra3x[i] = 0;}
        } else {repstra30x[i] = repstra3x[i];}
        repstra3[ndfindex] = repstra3x[i];
      }
      if (repstrb3col != -1) {
        if (NumericVector::is_na(repstrb3x[i])) {
          repstrb30x[i] = 0;
          if (NAas0) {repstrb3x[i] = 0;}
        } else {repstrb30x[i] = repstrb3x[i];}
        repstrb3[ndfindex] = repstrb3x[i];
      }
      repstradded3[ndfindex] = repstra30x[i] + (repstrb30x[i] * repstrrel);
      
      if (feca3col != -1) {
        if (NumericVector::is_na(feca3x[i])) {
          feca30x[i] = 0;
          if (NAas0) {feca3x[i] = 0;}
        } else {feca30x[i] = feca3x[i];}
        feca3[ndfindex] = feca3x[i];
      }
      if (fecb3col != -1) {
        if (NumericVector::is_na(fecb3x[i])) {
          fecb30x[i] = 0;
          if (NAas0) {feca3x[i] = 0;}
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
      
      if (alivegiven3x[i] > 0) {alive3[ndfindex] = 1;} else if (alivegiven3x[i] == 0) {alive3[ndfindex] = 0;}
      if (deadgiven3x[i] > 0) {alive3[ndfindex] = 0;} else if (deadgiven3x[i] == 0) {alive3[ndfindex] = 1;}
      
      livcheck3 = sizeadded3[ndfindex] + repstradded3[ndfindex] + obsstatus3[ndfindex];
      
      if (firstseenx[currentindiv] == -1 && livcheck3 > 0) {
        firstseenx[currentindiv] = currentyear + firstyear + 1;
        lastseenx[currentindiv] = currentyear + firstyear + 1;
        alive3[ndfindex] = 1;
        if (juv3col == -1 && repstradded2[ndfindex] > 0) {matstat3[ndfindex] = 1;}
      } else if (livcheck3 > 0) {
        lastseenx[currentindiv] = currentyear + firstyear + 1;
        alive3[ndfindex] = 1;
        if (juv3col == -1 && repstradded2[ndfindex] > 0) {matstat3[ndfindex] = 1;}
      }
      
      if (stage3col != -1 && alive3[ndfindex] == 1) {
        stage3[ndfindex] = stage2x[i];
      } else if (stassign) {stage3[ndfindex] = "NotAlive";}
    } // End of currentyear if statement
  } // End of i loop
  
  
  // Now a loop that establishes most states in time t+1 and t-1, and stages in all times
  for (int i = 0; i < ndflength; i++) { // Here variable i refers to rows in the final dataset
    
    // This short section deals with correcting info for individuals that are unobserved for long periods
    if (i > 0 && rowid[i] == 0) {
      if (year2[i-1] < lastseen[i-1] && (year2[i-1] + 1) < (firstyear + noyears)) {
        individ[i] = individ[i-1];
        popid[i] = popid[i-1];
        patchid[i] = patchid[i-1];
        year2[i] = year2[i-1] + 1;
        
        if (matstat2[i-1] == 1) {
          matstat2[i] = 1;
          
          if (year2[i+1] <= lastseen[i]) {
            matstat3[i] = 1;
          }
        }
      }
    }
    // valgrind mentioned this line
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
        alive2[i] = 1;
      } else if (alive2[i] == 0) {
        alive2[i] = 0;
      }
      
      if ((year2[i] + 1) >= firstseen[i] && (year2[i] + 1) <= lastseen[i] && alive3[i] == 0) {
        alive3[i] = 1;
      } else if (alive3[i] == 0) {
        alive3[i] = 0;
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
        
        if (censor2[i] == censorkeep && alive2[i] == 1) {      // New section through 2095
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
        } else if (stage2col != -1 && alive3[i] == 0 && stassign) {stage3[i] = "NotAlive";}
      }
      
      if (currentyear > 0  && year2[i] < (firstyear + noyears)) {
        prevyrindex = (noyears * currentindiv) + (currentyear - 1);
        
        alive1[i] = alive2[prevyrindex];
        
        if (censor2(i) == censorkeep && alive1(i) == 1) {  // New section through 2153
          if (censbool) {
            censor1(i) = 0;
          } else {
            censor1[i] = censorkeep;
          }
        } else {
          censor1[i] = censor2[prevyrindex];
        }
        
        xpos1[i] = xpos2[prevyrindex];
        ypos1[i] = ypos3[prevyrindex];
        
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
          } else {stage1[i] = stage2[prevyrindex];}
        } else if (stage2col != -1 && alive1[i] == 0 && stassign) {stage1[i] = "NotAlive";}
      }
      
      // Stage assignments
      if (stassign && stage2col == -1) {
        
        if (stszcol == 4) {
          
          stagesize1 = sizeadded1[i];
          stagesize2 = sizeadded2[i];
          stagesize3 = sizeadded3[i];
          
        } else if (stszcol == 3) {
          
          stagesize1 = sizec1[i];
          stagesize2 = sizec2[i];
          stagesize3 = sizec3[i];
          
        } else if (stszcol == 2) {
          
          stagesize1 = sizeb1[i];
          stagesize2 = sizeb2[i];
          stagesize3 = sizeb3[i];
          
        } else {
          
          stagesize1 = sizea1[i];
          stagesize2 = sizea2[i];
          stagesize3 = sizea3[i];
        }
        
        // Stage 2
        stagemini = find(insfszminarma < stagesize2);
        stagemaxi = find(insfszmaxarma >= stagesize2);
        stagerep = find(inrepstatarma == repstatus2[i]);
        stagemat = find(inmatstatarma == matstat2[i]);
        stageobs = find(inobsstatarma == obsstatus2[i]);
        
        cs1 = intersect(stagemini, stagemaxi);
        cs2 = intersect(stageobs, stagemat);
        
        if (NRasRep) {
          cs4 = intersect(cs1, cs2);
        } else {
          cs3 = intersect(stagerep, cs2);
          cs4 = intersect(cs1, cs3);
        }
        
        if (cs4.n_elem > 0 && alive2[i] == 1) {
          choicestage = instageid(cs4[0]) - 1;
          stage2num[i] = choicestage + 1;
          
          stage2[i] = sfname[choicestage];
        } else {
          stage2[i] = "NotAlive";
        }
        
        // Stage 1
        stagemini = find(insfszminarma < stagesize1);
        stagemaxi = find(insfszmaxarma >= stagesize1);
        stagerep = find(inrepstatarma == repstatus1[i]);
        stagemat = find(inmatstatarma == matstat1[i]);
        stageobs = find(inobsstatarma == obsstatus1[i]);
        
        cs1 = intersect(stagemini, stagemaxi);
        cs2 = intersect(stageobs, stagemat);
        
        if (NRasRep) {
          cs4 = intersect(cs1, cs2);
        } else {
          cs3 = intersect(stagerep, cs2);
          cs4 = intersect(cs1, cs3);
        }
        
        if (cs4.n_elem > 0 && alive1[i] == 1) {
          choicestage = instageid(cs4[0]) - 1;
          stage1num[i] = choicestage + 1;
          
          stage1[i] = sfname[choicestage];
        } else {
          stage1[i] = "NotAlive";
          matstat1[i] = 0;
        }
        
        // Stage 3
        stagemini = find(insfszminarma < stagesize3);
        stagemaxi = find(insfszmaxarma >= stagesize3);
        stagerep = find(inrepstatarma == repstatus3[i]);
        stagemat = find(inmatstatarma == matstat3[i]);
        stageobs = find(inobsstatarma == obsstatus3[i]);
        
        cs1 = intersect(stagemini, stagemaxi);
        cs2 = intersect(stageobs, stagemat);
        
        if (NRasRep) {
          cs4 = intersect(cs1, cs2);
        } else {
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
          
          Rf_warningcall(R_NilValue, "Some stages occurring in the dataset do not match any characteristics in the input stageframe.");
          
        } else if (cs4.n_elem > 1) {
          
          Rf_warningcall(R_NilValue, "Some stages in the input stageframe appear to have the same description. Please make sure that all stages included in the stageframe are defined with unique sets of characteristics.");
          
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
            censor2(indivzeros(j)) = 0;
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
                censor1(nextvec(0)) = 0;
              } else {
                censor1(nextvec(0)) = censorkeep;
              }
            }
            
            if (priorvec.n_elem > 0) {
              if (censbool) {
                censor3(priorvec(0)) = 0;
              } else {
                censor3(priorvec(0)) = censorkeep;  // valgrind mentioned this line
              }
            } // priorvec if statement
          } // currenttracking if statement
        } // newcount for loop
      } // newcount if statement
    } // noindivs for loop
  } // end censor correction section
  
  Rcpp::DataFrame df0 = DataFrame::create(Named("rowid") = rowid,
    _["popid"] = popid, _["patchid"] = patchid, _["individ"] = individ,
    _["year2"] = year2, _["firstseen"] = firstseen, _["lastseen"] = lastseen,
    _["obsage"] = obsage, _["obslifespan"] = obslifespan);
  
  Rcpp::DataFrame df1 = DataFrame::create(Named("xpos1") = xpos1,
    _["ypos1"] = ypos1, _["sizea1"] = sizea1, _["sizeb1"] = sizeb1,
    _["sizec1"] = sizec1, _["sizeadded1"] = sizeadded1,
    _["repstra1"] = repstra1, _["repstrb1"] = repstrb1,
    _["repstradded1"] = repstradded1, _["feca1"] = feca1, _["fecb1"] = fecb1,
    _["fecaadded1"] = fecadded1, _["indcova1"] = indcova1,
    _["indcovb1"] = indcovb1, _["indcovc1"] = indcovc1, _["censor1"] = censor1,
    _["juvgiven1"] = juvgiven1);
  
  Rcpp::DataFrame df1a = DataFrame::create(Named("obsstatus1") = obsstatus1,
    _["repstatus1"] = repstatus1, _["fecstatus1"] = fecstatus1,
    _["matstatus1"] = matstat1, _["alive1"] = alive1, _["stage1"] = stage1,
    _["stage1index"] = stage1num);
  
  Rcpp::DataFrame df2 = DataFrame::create(Named("xpos2") = xpos2,
    _["ypos2"] = ypos2, _["sizea2"] = sizea2, _["sizeb2"] = sizeb2,
    _["sizec2"] = sizec2, _["sizeadded2"] = sizeadded2,
    _["repstra2"] = repstra2, _["repstrb2"] = repstrb2,
    _["repstradded2"] = repstradded2, _["feca2"] = feca2, _["fecb2"] = fecb2,
    _["fecaadded2"] = fecadded2, _["indcova2"] = indcova2,
    _["indcovb2"] = indcovb2, _["indcovc2"] = indcovc2, _["censor2"] = censor2,
    _["juvgiven2"] = juvgiven2);
  
  Rcpp::DataFrame df2a = DataFrame::create(Named("obsstatus2") = obsstatus2,
    _["repstatus2"] = repstatus2, _["fecstatus2"] = fecstatus2,
    _["matstatus2"] = matstat2, _["alive2"] = alive2, _["stage2"] = stage2,
    _["stage2index"] = stage2num);
  
  Rcpp::DataFrame df3 = DataFrame::create(Named("xpos3") = xpos3,
    _["ypos3"] = ypos3, _["sizea3"] = sizea3, _["sizeb3"] = sizeb3,
    _["sizec3"] = sizec3, _["sizeadded3"] = sizeadded3,
    _["repstra3"] = repstra3, _["repstrb3"] = repstrb3,
    _["repstradded3"] = repstradded3, _["feca3"] = feca3, _["fecb3"] = fecb3,
    _["fecaadded3"] = fecadded3, _["indcova3"] = indcova3,
    _["indcovb3"] = indcovb3, _["indcovc3"] = indcovc3, _["censor3"] = censor3,
    _["juvgiven3"] = juvgiven3);
  
  Rcpp::DataFrame df3a = DataFrame::create(Named("obsstatus3") = obsstatus3,
    _["repstatus3"] = repstatus3, _["fecstatus3"] = fecstatus3,
    _["matstatus3"] = matstat3, _["alive3"] = alive3, _["stage3"] = stage3,
    _["stage3index"] = stage3num);
  
  return Rcpp::List::create(Named("a") = df0, Named("b") = df1,
    Named("c") = df1a, Named("d") = df2, Named("e") = df2a, Named("f") = df3,
    Named("g") = df3a);
}

//' Create Element Index for Matrix Estimation
//' 
//' Function \code{theoldpizzle()} create a data frame object spread across
//' three 20-element long list objects that is used by \code{jerzeibalowski()}
//' to estimate function-derived matrices, and also provides the element map for
//' raw matrices.
//'
//' @param StageFrame The stageframe object identifying the life history model
//' being operationalized.
//' @param OverWrite The overwrite table used in analysis, as modified by 
//' \code{.overwrite_reassess}. Must be processed via \code{.overwrite_reassess}
//' rather than being a raw overwrite or supplement table.
//' @param repmatrix The reproductive matrix used in analysis.
//' @param finalage The final age to be used in analysis.
//' @param format Indicates whether historical matrices should be in (1) Ehrlen
//' or (2) deVries format.
//' @param style The style of analysis, where 0 is historical, 1 is ahistorical,
//' and 2 is age-by-stage.
//' @param cont Denotes whether age-by-stage matrix continues past the final age.
//' 
//' @return The output is currently a three element list, where each element is
//' a data frame with the same number of rows. These are combined into one data
//' frame by \code{\link{flefko3}()}, \code{\link{flefko2}()}, and
//' \code{\link{aflefko2}()}.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
Rcpp::List theoldpizzle(DataFrame StageFrame, DataFrame OverWrite,
  arma::mat repmatrix, int finalage, int format, int style, int cont) {
  
  StringVector ovstage3 = OverWrite["stage3"];
  StringVector ovstage2 = OverWrite["stage2"];
  StringVector ovstage1 = OverWrite["stage1"];
  StringVector oveststage3 = OverWrite["eststage3"];
  StringVector oveststage2 = OverWrite["eststage2"];
  StringVector oveststage1 = OverWrite["eststage1"];
  arma::vec ovgivenrate = OverWrite["givenrate"];
  arma::vec ovmultiplier = OverWrite["multiplier"];
  arma::vec ovconvtype = OverWrite["convtype"];
  arma::vec ovconvt12 = OverWrite["convtype_t12"];
  int ovrows = ovconvtype.n_elem;
  
  int totalages = finalage + 1;
  
  arma::vec ovindex3(ovrows * totalages);
  arma::vec ovindex2(ovrows * totalages);
  arma::vec ovindex1(ovrows * totalages);
  arma::vec ovnew3(ovrows * totalages);
  arma::vec ovnew2(ovrows * totalages);
  arma::vec ovnew1(ovrows * totalages);
  arma::vec ovindexold321(ovrows * totalages);
  arma::vec ovindexnew321(ovrows * totalages);
  arma::vec ovnewgivenrate(ovrows * totalages);
  arma::vec ovnewmultiplier(ovrows * totalages);
  arma::vec ovconvtypeage(ovrows * totalages);
  ovindex3.fill(-1);
  ovindex2.fill(-1);
  ovindex1.fill(-1);
  ovnew3.fill(-1);
  ovnew2.fill(-1);
  ovnew1.fill(-1);
  ovindexold321.fill(-1);
  ovindexnew321.fill(-1);
  ovnewgivenrate.fill(-1);
  ovnewmultiplier.zeros();
  ovconvtypeage.fill(-1);
  
  arma::vec newstageid = StageFrame["stage_id"];
  StringVector origstageid = StageFrame["stage"];
  arma::vec binsizectr = StageFrame["bin_size_ctr"];
  arma::vec repstatus = StageFrame["repstatus"];
  arma::vec obsstatus = StageFrame["obsstatus"];
  arma::vec immstatus = StageFrame["immstatus"];
  arma::vec matstatus = StageFrame["matstatus"];
  arma::vec indata = StageFrame["indataset"];
  arma::vec binsizewidth = StageFrame["bin_size_width"];
  arma::vec alive = StageFrame["alive"];
  arma::vec minage = StageFrame["min_age"];
  arma::vec maxage = StageFrame["max_age"];
  
  // This section determines the length of the matrix map data frame
  int nostages = newstageid.n_elem;
  int nostages_nodead = nostages - 1;
  int nostages_nounborn = nostages;
  int nostages_nodead_nounborn = nostages_nodead;
  int prior_stage = -1;
  arma::vec ovrepentry_prior(nostages);
  ovrepentry_prior.zeros();
  
  int totallength {0};
  
  if (style == 2) {
    totallength = (nostages * nostages * (finalage + 1) * (finalage + 1));
  } else if (style == 1) {
    totallength = (nostages * nostages_nodead);
  } else {
    if (format == 2) {
      nostages_nodead_nounborn = nostages - 2;
      prior_stage = nostages_nodead_nounborn;
      nostages_nounborn = nostages - 1;
      totallength = (2 * nostages_nodead_nounborn * nostages_nounborn * nostages_nounborn);
    } else {
      totallength = (nostages * (nostages_nodead * nostages_nodead));
    }
  }
  
  // This section sets up the repmatrix. First we will determine whether the
  // repmatrix has been entered in historical or ahistorical format, since this
  // does not necessarily match the MPM type
  int reprows = repmatrix.n_rows;

  int repmattype = 0;
  
  if (reprows == (nostages - 1) || reprows == (nostages - 2)) {
    repmattype = 1; // The repmatrix is ahistorical in dimensions
  } else if (reprows == ((nostages - 1) * (nostages - 1)) || reprows == ((nostages - 2) * (nostages - 2))) {
    repmattype = 2; // The repmatrix is historical in dimensions
  }
  
  // Here we set up the vectors that will be put together into the matrix map
  // data frame
  arma::vec stage3(totallength);
  arma::vec stage2n(totallength);
  arma::vec stage2o(totallength);
  arma::vec stage1(totallength);
  stage3.zeros();
  stage2n.zeros();
  stage2o.zeros();
  stage1.zeros();
  
  arma::vec size3(totallength);
  arma::vec size2n(totallength);
  arma::vec size2o(totallength);
  arma::vec size1(totallength);
  size3.zeros();
  size2n.zeros();
  size2o.zeros();
  size1.zeros();
  
  arma::vec obs3(totallength);
  arma::vec obs2n(totallength);
  arma::vec obs2o(totallength);
  arma::vec obs1(totallength);
  obs3.zeros();
  obs2n.zeros();
  obs2o.zeros();
  obs1.zeros();
  
  arma::vec rep3(totallength);
  arma::vec rep2n(totallength);
  arma::vec rep2o(totallength);
  arma::vec rep1(totallength);
  rep3.zeros();
  rep2n.zeros();
  rep2o.zeros();
  rep1.zeros();
  
  arma::vec mat3(totallength);
  arma::vec mat2n(totallength);
  arma::vec mat2o(totallength);
  arma::vec mat1(totallength);
  mat3.zeros();
  mat2n.zeros();
  mat2o.zeros();
  mat1.zeros();
  
  arma::vec imm3(totallength);
  arma::vec imm2n(totallength);
  arma::vec imm2o(totallength);
  arma::vec imm1(totallength);
  imm3.zeros();
  imm2n.zeros();
  imm2o.zeros();
  imm1.zeros();
  
  arma::vec repentry3(totallength);
  arma::vec binwidth(totallength);
  repentry3.zeros();
  binwidth.zeros();
  
  arma::vec indata3(totallength);
  arma::vec indata2n(totallength);
  arma::vec indata2o(totallength);
  arma::vec indata1(totallength);
  indata3.zeros();
  indata2n.zeros();
  indata2o.zeros();
  indata1.zeros();
  
  arma::vec minage3(totallength);
  arma::vec minage2(totallength);
  arma::vec maxage3(totallength);
  arma::vec maxage2(totallength);
  minage3.zeros();
  minage2.zeros();
  maxage3.zeros();
  maxage2.zeros();
  
  arma::vec actualage(totallength);
  arma::vec index321(totallength);
  arma::vec index21(totallength);
  arma::vec indatalong(totallength);
  arma::vec aliveequal(totallength);
  arma::vec included(totallength);
  actualage.zeros();
  index321.fill(-1);
  index21.fill(-1);
  indatalong.zeros();
  aliveequal.fill(-1);
  included.zeros();
  
  arma::mat asadditions(totallength, 5);
  arma::vec ovgivent(totallength);
  arma::vec ovestt(totallength);
  arma::vec ovgivenf(totallength);
  arma::vec ovestf(totallength);
  arma::vec ovrepentry(totallength);
  asadditions.zeros();
  ovgivent.fill(-1);
  ovestt.fill(-1);
  ovgivenf.fill(-1);
  ovestf.fill(-1);
  ovrepentry.zeros();
  
  int repm_elem {-1};
  double deadandnasty {0};
  // long long int lifeindex {0};
  long long int currentindex {0};

  // This step changes the stage names to stage numbers per the input stageframe for styles 0 and 1
  if (style < 2) {
    if (ovrows > 1 || ovconvtype(0) != -1) {
      for (int i = 0; i < ovrows; i++) { // Loop across overwrite rows
        for (int j = 0; j < nostages; j++) { // Loop across stageframe rows
          if (ovstage3(i) == origstageid(j)) {
            ovindex3(i) = newstageid(j);
          }
          
          if (ovstage2(i) == origstageid(j)) {
            ovindex2(i) = newstageid(j);
          }
          
          if (ovstage1(i) == origstageid(j)) {
            ovindex1(i) = newstageid(j);
          }
          
          if (oveststage3(i) == origstageid(j)) {
            ovnew3(i) = newstageid(j);
          }
          
          if (oveststage2(i) == origstageid(j)) {
            ovnew2(i) = newstageid(j);
          }
          
          if (oveststage1(i) == origstageid(j)) {
            ovnew1(i) = newstageid(j);
          }
        } // j for loop
      } // i for loop
    } // ovrows if statement
  } // style if statement
  
  // Now we cover the main data frame creation loops
  // When style = 0, this will create AllStages for the historical case
  if (style == 0 && format == 2) {
    
    if (ovrows > 1 || ovconvtype(0) != -1) {
      for (int i = 0; i < ovrows; i++) { // Loop across overwrite rows\
        
        // Here are the new versions
        if (ovconvtype(i) > 1) { // This catches all changes to fecundity and reproductive multipliers
          ovindexold321(i) = (ovindex3(i) - 1) + (prior_stage * nostages) + 
            ((ovindex2(i) - 1) * nostages * nostages) + 
            ((ovindex1(i) - 1) * nostages * nostages * nostages);
            
          ovindexnew321(i) = (ovnew3(i) - 1) + (prior_stage * nostages) + 
            ((ovnew2(i) - 1) * nostages * nostages) + 
            ((ovnew1(i) - 1) * nostages * nostages * nostages);
        } else if (ovconvt12(i) == 2) { // This catches all survival terms with historical reproduction events
          ovindexold321(i) = (ovindex3(i) - 1) + ((ovindex2(i) - 1) * nostages) + 
            ((ovindex2(i) - 1) * nostages * nostages) + 
            (prior_stage * nostages * nostages * nostages);
            
          ovindexnew321(i) = (ovnew3(i) - 1) + ((ovnew2(i) - 1) * nostages) + 
            ((ovnew2(i) - 1) * nostages * nostages) + 
            (prior_stage * nostages * nostages * nostages);
        } else { // This refers to full survival transitions
          ovindexold321(i) = (ovindex3(i) - 1) + ((ovindex2(i) - 1) * nostages) + 
            ((ovindex2(i) - 1) * nostages * nostages) + 
            ((ovindex1(i) - 1) * nostages * nostages * nostages);
            
          ovindexnew321(i) = (ovnew3(i) - 1) + ((ovnew2(i) - 1) * nostages) + 
            ((ovnew2(i) - 1) * nostages * nostages) + 
            ((ovnew1(i) - 1) * nostages * nostages * nostages);
        }
        if (ovindexold321(i) < 0) ovindexold321(i) = -1;
        if (ovindexnew321(i) < 0) ovindexnew321(i) = -1;
        
        if (!NumericVector::is_na(ovgivenrate(i))) {
          ovnewgivenrate(i) = ovgivenrate(i);
        }
        if (!NumericVector::is_na(ovmultiplier(i))) {
          ovnewmultiplier(i) = ovmultiplier(i);
        }
        
        if (ovconvtype(i) == 3) {
          for (int j = 0; j < nostages; j++) {
            if (origstageid(j) == ovstage3(i)) ovrepentry_prior(j) = 1;
          }
        }
      } // i for loop
    } // ovrows if statement
    
    for (int time1 = 0; time1 < nostages_nodead; time1++) {
      for (int time2o = 0; time2o < nostages_nodead_nounborn; time2o++) {
        for (int time2n = 0; time2n < nostages; time2n++) {
          for (int time3 = 0; time3 < nostages; time3++) {
            
            if (time3 != prior_stage) {
              if (time2n == time2o || time2n == prior_stage){
                
                included(currentindex) = 1;
                
                stage3(currentindex) = newstageid(time3);
                stage2n(currentindex) = newstageid(time2n);
                stage2o(currentindex) = newstageid(time2o);
                stage1(currentindex) = newstageid(time1);
                size3(currentindex) = binsizectr(time3);
                size2n(currentindex) = binsizectr(time2n);
                size2o(currentindex) = binsizectr(time2o);
                size1(currentindex) = binsizectr(time1);
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
                
                // This statement fills in the repentry info from the repmatrix
                if (time2n == prior_stage && time3 < prior_stage && time2o < prior_stage) {
                  if (repmattype == 1) { 
                    repm_elem = time3 + (time2o * nostages_nodead_nounborn);
                  } else if (repmattype == 2) {
                    repm_elem = time3 + (time2o * nostages_nodead_nounborn) + 
                      (time2o * nostages_nodead_nounborn * nostages_nodead_nounborn) +
                      (time1 * nostages_nodead_nounborn * nostages_nodead_nounborn * nostages_nodead_nounborn);
                  } else repm_elem = -1;
                  
                  if (repmatrix(repm_elem) > 0) {
                    repentry3(currentindex) = repmatrix(repm_elem);
                    if (repentry3(currentindex) == 0 && ovrepentry_prior(time3) != 0) {
                      repentry3(currentindex) = 1;
                    } 
                  }
                } else repentry3(currentindex) = 0;
                
                indata3(currentindex) = indata(time3);
                indata2n(currentindex) = indata(time2n);
                indata2o(currentindex) = indata(time2o);
                indata1(currentindex) = indata(time1);
                binwidth(currentindex) = binsizewidth(time3);
                minage3(currentindex) = minage(time3);
                minage2(currentindex) = minage(time2o);
                maxage3(currentindex) = maxage(time3);
                maxage2(currentindex) = maxage(time2o);
                actualage(currentindex) = 0;
                
                if (stage3(currentindex) == nostages || stage2n(currentindex) == nostages) {
                  deadandnasty = 1;
                } else if (stage2o(currentindex) == nostages || stage1(currentindex) == nostages) {
                  deadandnasty = 1;
                } else {
                  deadandnasty = 0;
                }
                
                if (deadandnasty == 0) {
                  // The next index variable gives the element in the final matrix
                  aliveequal(currentindex) = (stage3(currentindex) - 1) + 
                    ((stage2n(currentindex) - 1) * nostages_nodead_nounborn) + 
                    ((stage2o(currentindex) - 1) * nostages_nodead * nostages_nodead_nounborn) + 
                    ((stage1(currentindex) - 1) * nostages_nodead_nounborn * nostages_nodead * nostages_nodead_nounborn);
                  
                  // The next two index variables are used by ovreplace
                  index321(currentindex) = (stage3(currentindex) - 1) + 
                    ((stage2n(currentindex) - 1) * nostages) + 
                    ((stage2o(currentindex) - 1) * nostages * nostages) + 
                    ((stage1(currentindex) - 1) * nostages * nostages * nostages);
                    
                  index21(currentindex) = (stage2o(currentindex) - 1) + 
                    ((stage1(currentindex) - 1) * nostages);
                }
                
                indatalong(currentindex) = indata3(currentindex) * indata2n(currentindex) * 
                  indata2o(currentindex) * indata1(currentindex);
                
                currentindex += 1;
              } // if (time2n == tim2o || time2n == prior_stage) statement
            } // if (time3n != dead_stage) statement
          } // time3 loop
        } // time2n loop
      } // time2o loop
    } // time1 loop 
    
    if (ovrows > 1 || ovconvtype(0) != -1) {
      asadditions = ovreplace(index321, ovindexold321, ovindexnew321, ovconvtype, ovnew3, ovnewgivenrate, ovnewmultiplier);
      
      ovgivent = asadditions.col(0);
      ovestt = asadditions.col(1);
      ovgivenf = asadditions.col(2);
      ovestf = asadditions.col(3);
      
      ovrepentry = asadditions.col(4);
      
      arma::uvec workedupindex = find(ovrepentry > 0);
      int changedreps = workedupindex.n_elem;
      
      if (changedreps > 0) {
        for (int i = 0; i < changedreps; i++) {
          repentry3(workedupindex(i)) = ovrepentry(workedupindex(i));
        }
      }
    } // ovreplace if statement
  } else if (style == 0 && format == 1) { // Historical MPM in Ehrlen format
    
    if (ovrows > 1 || ovconvtype(0) != -1) {
      for (int i = 0; i < ovrows; i++) { // Loop across overwrite rows\
        
        // Here are the new versions
        ovindexold321(i) = (ovindex3(i) - 1) + ((ovindex2(i) - 1) * nostages_nodead_nounborn) + 
          ((ovindex2(i) - 1) * nostages_nodead_nounborn * nostages_nodead_nounborn) + 
          ((ovindex1(i) - 1) * nostages_nodead_nounborn * nostages_nodead_nounborn * nostages_nodead_nounborn);
          
        ovindexnew321(i) = (ovnew3(i) - 1) + ((ovnew2(i) - 1) * nostages_nodead) + 
          ((ovnew2(i) - 1) * nostages_nodead_nounborn * nostages_nodead_nounborn) + 
          ((ovnew1(i) - 1) * nostages_nodead_nounborn * nostages_nodead_nounborn * nostages_nodead_nounborn);
        
        if (ovindexold321(i) < 0) ovindexold321(i) = -1;
        if (ovindexnew321(i) < 0) ovindexnew321(i) = -1;
        
        if (!NumericVector::is_na(ovgivenrate(i))) {
          ovnewgivenrate(i) = ovgivenrate(i);
        }
        if (!NumericVector::is_na(ovmultiplier(i))) {
          ovnewmultiplier(i) = ovmultiplier(i);
        }
      } // i for loop
    } // ovrows if statement
    
    for (int time1 = 0; time1 < nostages_nodead; time1++) {
      for (int time2o = 0; time2o < nostages_nodead; time2o++) {
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
          
          //This next section determines repentry3 on the basis of the input repmatrix
          if (time3 < nostages_nodead_nounborn) {
            if (repmattype == 1) {
              repm_elem = time3 + (time2o * nostages_nodead_nounborn);
            } else if (repmattype == 2) {
              repm_elem = time3 + (time2o * nostages_nodead_nounborn) + 
                (time2o * nostages_nodead_nounborn * nostages_nodead_nounborn) +
                (time1 * nostages_nodead_nounborn * nostages_nodead_nounborn * nostages_nodead_nounborn);
            } else {
              repm_elem = -1;
            }
          }
          
          if(repm_elem > -1) {
            if (repmatrix(repm_elem) > 0) {
              repentry3(currentindex) = repmatrix(repm_elem);
            }
          }
          
          if (time3 < nostages_nodead_nounborn) {
            if (repmattype == 1) { // Ahistorical repmatrix
              repentry3(currentindex) = repmatrix((time3 + (nostages_nodead_nounborn * time2o)));
            } else if (repmattype == 2) {  // Historical repmatrix
              repentry3(currentindex) = repmatrix((time3 + (nostages_nodead_nounborn * time2o)) + 
                ((nostages_nodead_nounborn * nostages_nodead_nounborn * time2o)) +
                (nostages_nodead_nounborn * nostages_nodead_nounborn * nostages_nodead_nounborn * time1));
            }
          } else {
            repentry3(currentindex) = 0;
          }
          
          indata3(currentindex) = indata(time3);
          indata2n(currentindex) = indata(time2o);
          indata2o(currentindex) = indata(time2o);
          indata1(currentindex) = indata(time1);
          binwidth(currentindex) = binsizewidth(time3);
          minage3(currentindex) = minage(time3);
          minage2(currentindex) = minage(time2o);
          maxage3(currentindex) = maxage(time3);
          maxage2(currentindex) = maxage(time2o);
          actualage(currentindex) = 0;
          
          if (stage3(currentindex) == nostages || stage2n(currentindex) == nostages) {
            deadandnasty = 1;
          } else if (stage2o(currentindex) == nostages || stage1(currentindex) == nostages) {
            deadandnasty = 1;
          } else {
            deadandnasty = 0;
          }
          
          if (deadandnasty == 0) {
            aliveequal(currentindex) = (stage3(currentindex) - 1) + ((stage2n(currentindex) - 1) * (nostages - 1)) + 
              ((stage2o(currentindex) - 1) * (nostages - 1) * (nostages - 1)) + 
              ((stage1(currentindex) - 1) * (nostages - 1) * (nostages - 1) * (nostages - 1));
            
            index321(currentindex) = (stage3(currentindex) - 1) + 
              ((stage2n(currentindex) - 1) * nostages_nodead_nounborn) + 
              ((stage2n(currentindex) - 1) * nostages_nodead_nounborn * nostages_nodead_nounborn) + 
              ((stage1(currentindex) - 1) * nostages_nodead_nounborn * nostages_nodead_nounborn * 
                nostages_nodead_nounborn);
            index21(currentindex) = (stage2n(currentindex) - 1) + ((stage1(currentindex) - 1) * nostages);
          }
          
          indatalong(currentindex) = indata3(currentindex) * indata2n(currentindex) * 
            indata2o(currentindex) * indata1(currentindex);
          
          currentindex += 1;
        } // time3 loop
      } // time2o loop
    } // time1 loop 
    
    if (ovrows > 1 || ovconvtype(0) != -1) {
      asadditions = ovreplace(index321, ovindexold321, ovindexnew321, ovconvtype, ovnew3, ovnewgivenrate, ovnewmultiplier);
      
      ovgivent = asadditions.col(0);
      ovestt = asadditions.col(1);
      ovgivenf = asadditions.col(2);
      ovestf = asadditions.col(3);
      
      ovrepentry = asadditions.col(4);
      
      arma::uvec workedupindex = find(ovrepentry > 0);
      int changedreps = workedupindex.n_elem;
      
      if (changedreps > 0) {
        for (int i = 0; i < changedreps; i++) {
          repentry3(workedupindex(i)) = ovrepentry(workedupindex(i));
        }
      }
    } // ovreplace if statement
  } else if (style == 1) { // This will take care of the ahistorical case
    
    if (ovrows > 1 || ovconvtype(0) != -1) {
      for (int i = 0; i < ovrows; i++) { // Loop across overwrite rows
      
        ovindexold321(i) = (ovindex3(i) - 1) + ((ovindex2(i) - 1) * nostages);
        ovindexnew321(i) = (ovnew3(i) - 1) + ((ovnew2(i) - 1) * nostages);
        
        if (ovindexold321(i) < 0) ovindexold321(i) = -1;
        if (ovindexnew321(i) < 0) ovindexnew321(i) = -1;
        
        if (!NumericVector::is_na(ovgivenrate(i))) {
          ovnewgivenrate(i) = ovgivenrate(i);
        }
        if (!NumericVector::is_na(ovmultiplier(i))) {
          ovnewmultiplier(i) = ovmultiplier(i);
        }
      } // i for loop
    } // ovrows if statement
    
    for (int time2n = 0; time2n < nostages_nodead; time2n++) {
      for (int time3 = 0; time3 < nostages; time3++) {
        
        stage3(currentindex) = newstageid(time3);
        stage2n(currentindex) = newstageid(time2n);
        stage2o(currentindex) = newstageid(time2n);
        stage1(currentindex) = 0;
        size3(currentindex) = binsizectr(time3);
        size2n(currentindex) = binsizectr(time2n);
        size2o(currentindex) = binsizectr(time2n);
        size1(currentindex) = 0;
        obs3(currentindex) = obsstatus(time3);
        obs2n(currentindex) = obsstatus(time2n);
        obs2o(currentindex) = obsstatus(time2n);
        obs1(currentindex) = 0;
        rep3(currentindex) = repstatus(time3);
        rep2n(currentindex) = repstatus(time2n);
        rep2o(currentindex) = repstatus(time2n);
        rep1(currentindex) = 0;
        mat3(currentindex) = matstatus(time3);
        mat2n(currentindex) = matstatus(time2n);
        mat2o(currentindex) = matstatus(time2n);
        mat1(currentindex) = matstatus(time2n);
        imm3(currentindex) = immstatus(time3);
        imm2n(currentindex) = immstatus(time2n);
        imm2o(currentindex) = immstatus(time2n);
        imm1(currentindex) = immstatus(time2n);
        
        if (time3 < nostages_nodead) {
          repentry3(currentindex) = repmatrix((time3 + (nostages_nodead * time2n)));
        } else {
          repentry3(currentindex) = 0;
        }
        
        indata3(currentindex) = indata(time3);
        indata2n(currentindex) = indata(time2n);
        indata2o(currentindex) = indata(time2n);
        indata1(currentindex) = 1;
        binwidth(currentindex) = binsizewidth(time3);
        minage3(currentindex) = minage(time3);
        minage2(currentindex) = minage(time2n);
        maxage3(currentindex) = maxage(time3);
        maxage2(currentindex) = maxage(time2n);
        actualage(currentindex) = 0;
        
        if (stage3(currentindex) == nostages || stage2n(currentindex) == nostages) {
          deadandnasty = 1;
        } else {
          deadandnasty = 0;
        }
        
        if (deadandnasty == 0) {
          aliveequal(currentindex) = (stage3(currentindex) - 1) + ((stage2n(currentindex) - 1) * nostages_nodead);

          index321(currentindex) = (stage3(currentindex) - 1) + ((stage2n(currentindex) - 1) * nostages);
          index21(currentindex) = (stage2n(currentindex) - 1);
        }
        
        indatalong(currentindex) = indata3(currentindex) * indata2n(currentindex) * 
          indata2o(currentindex);
          
        currentindex += 1;
        
      } // time3 loop
    } // time2n loop
    
    if (ovrows > 1 || ovconvtype(0) != -1) {
      asadditions = ovreplace(index321, ovindexold321, ovindexnew321, ovconvtype, ovnew3, ovnewgivenrate, ovnewmultiplier);
      
      ovgivent = asadditions.col(0);
      ovestt = asadditions.col(1);
      ovgivenf = asadditions.col(2);
      ovestf = asadditions.col(3);
      
      ovrepentry = asadditions.col(4);
      
      arma::uvec workedupindex = find(ovrepentry > 0);
      int changedreps = workedupindex.n_elem;
      
      if (changedreps > 0) {
        for (int i = 0; i < changedreps; i++) {
          repentry3(workedupindex(i)) = ovrepentry(workedupindex(i));
        }
      }
    } // ovreplace if statement
  } else if (style == 2) { // This takes care of the stage x age case
    
    int age3 {0};
    
    for (int time3 = 0; time3 < nostages; time3++) {
      if (NumericVector::is_na(maxage(time3))) {
        maxage(time3) = finalage + cont;
      }
    }
    
    // This sets up the overwrite tables
    if (ovrows > 1 || ovconvtype(0) != -1) {
      // This first set of loops establishes a number of indices
      for (int age2 = 0; age2 < totalages; age2++) {
        for (int i = 0; i < ovrows; i++) { // Loop across overwrite rows
          for (int j = 0; j < nostages; j++) { // Loop across stageframe rows
            ovconvtypeage(i + (ovrows * age2)) = ovconvtype(i);
              
            if (age2 < (totalages - 1)) {
              if (ovconvtype(i) == 1) {
                age3 = age2 + 1;
              } else {
                age3 = 0;
              }
              
              if (ovstage3(i) == origstageid(j)) {
                ovindex3(i + (ovrows * age2)) = newstageid(j) - 1;
              }
              
              if (ovstage2(i) == origstageid(j)) {
                ovindex2(i + (ovrows * age2)) = newstageid(j) - 1;
              }
              
              if (oveststage3(i) == origstageid(j)) {
                ovnew3(i + (ovrows * age2)) = newstageid(j) - 1;
              }
              
              if (oveststage2(i) == origstageid(j)) {
                ovnew2(i + (ovrows * age2)) = newstageid(j) - 1;
              }
              
              if (ovindex3(i + (ovrows * age2)) != -1 && ovindex2(i + (ovrows * age2)) != -1) {
                ovindexold321(i + (ovrows * age2)) = ovindex3(i + (ovrows * age2)) + (age3 * nostages) + 
                  (ovindex2(i + (ovrows * age2)) * nostages * totalages) + (age2 * nostages * nostages * totalages);
              }
              
              if (ovnew3(i + (ovrows * age2)) != -1 && ovnew2(i + (ovrows * age2)) != -1) {
                ovindexnew321(i + (ovrows * age2)) = ovnew3(i + (ovrows * age2)) + (age3 * nostages) + 
                  (ovnew2(i + (ovrows * age2)) * nostages * totalages) + (age2 * nostages * nostages * totalages);
              }
              
              if (!NumericVector::is_na(ovgivenrate(i))) {
                ovnewgivenrate(i + (ovrows * age2)) = ovgivenrate(i);
              }
              if (!NumericVector::is_na(ovmultiplier(i))) {
                ovnewmultiplier(i + (ovrows * age2)) = ovmultiplier(i);
              }
            } else {
              if (ovconvtype(i) == 1) {
                age3 = age2;
              } else {
                age3 = 0;
              }
              
              if (ovstage3(i) == origstageid(j)) {
                ovindex3(i + (ovrows * age2)) = newstageid(j) - 1;
              }
              
              if (ovstage2(i) == origstageid(j)) {
                ovindex2(i + (ovrows * age2)) = newstageid(j) - 1;
              }
              
              if (oveststage3(i) == origstageid(j)) {
                ovnew3(i + (ovrows * age2)) = newstageid(j) - 1;
              }
              
              if (oveststage2(i) == origstageid(j)) {
                ovnew2(i + (ovrows * age2)) = newstageid(j) - 1;
              }
              
              if (ovindex3(i + (ovrows * age2)) != -1 && ovindex2(i + (ovrows * age2)) != -1) {
                ovindexold321(i + (ovrows * age2)) = ovindex3(i + (ovrows * age2)) + (age3 * nostages) + 
                  (ovindex2(i + (ovrows * age2)) * nostages * totalages) + (age2 * nostages * nostages * totalages);
              }
              
              if (ovnew3(i + (ovrows * age2)) != -1 && ovnew2(i + (ovrows * age2)) != -1) {
                ovindexnew321(i + (ovrows * age2)) = ovnew3(i + (ovrows * age2)) + (age3 * nostages) + 
                  (ovnew2(i + (ovrows * age2)) * nostages * totalages) + (age2 * nostages * nostages * totalages);
              }
              if (!NumericVector::is_na(ovgivenrate(i))) {
                ovnewgivenrate(i + (ovrows * age2)) = ovgivenrate(i);
              }
              if (!NumericVector::is_na(ovmultiplier(i))) {
                ovnewmultiplier(i + (ovrows * age2)) = ovmultiplier(i);
              }
            }
          } // j for loop
          
        if (ovindexold321(i) < 0) ovindexold321(i) = -1;
        if (ovindexnew321(i) < 0) ovindexnew321(i) = -1;
          
        } // i for loop
      }
    } // ovrows if statement
    
    for (int age2 = 0; age2 <= finalage; age2++) {
      if (age2 < finalage) { // This first loop takes care of transitions from one age to the next
        for (int time2n = 0; time2n < nostages; time2n++) {
          for (int time3 = 0; time3 < nostages; time3++) {
            
            // First survival
            age3 = age2 + 1;
            currentindex = time3 + (age3 * nostages) + 
              (time2n * nostages * totalages) + (age2 * nostages * nostages * totalages);
            
            stage3(currentindex) = newstageid(time3);
            stage2n(currentindex) = newstageid(time2n);
            stage2o(currentindex) = newstageid(time2n);
            stage1(currentindex) = 0;
            size3(currentindex) = binsizectr(time3);
            size2n(currentindex) = binsizectr(time2n);
            size2o(currentindex) = binsizectr(time2n);
            size1(currentindex) = 0;
            obs3(currentindex) = obsstatus(time3);
            obs2n(currentindex) = obsstatus(time2n);
            obs2o(currentindex) = obsstatus(time2n);
            obs1(currentindex) = 0;
            rep3(currentindex) = repstatus(time3);
            rep2n(currentindex) = repstatus(time2n);
            rep2o(currentindex) = repstatus(time2n);
            rep1(currentindex) = 0;
            mat3(currentindex) = matstatus(time3);
            mat2n(currentindex) = matstatus(time2n);
            mat2o(currentindex) = matstatus(time2n);
            mat1(currentindex) = 0;
            imm3(currentindex) = immstatus(time3);
            imm2n(currentindex) = immstatus(time2n);
            imm2o(currentindex) = immstatus(time2n);
            imm1(currentindex) = 0;
            
            repentry3(currentindex) = 0;
            
            indata3(currentindex) = indata(time3);
            indata2n(currentindex) = indata(time2n);
            indata2o(currentindex) = indata(time2n);
            indata1(currentindex) = 0;
            binwidth(currentindex) = binsizewidth(time3);
            minage3(currentindex) = minage(time3);
            minage2(currentindex) = minage(time2n);
            maxage3(currentindex) = maxage(time3);
            maxage2(currentindex) = maxage(time2n);
            actualage(currentindex) = age2;
            
            // The next indexer includes the following order: (1st # of age blocks) + (1st # of stage cols) +
            // (1st # of age rows) + stage in time 3
            index321(currentindex) = currentindex;
            indatalong(currentindex) = 1;
            
            // This section identifies elements with non-zero entries by their element number in the final matrix
            if (alive(time2n) == 1 && alive(time3) == 1) {
              if (age2 >= minage2(currentindex) && age2 < maxage3(currentindex)) { 
                
                // Survival transitions
                aliveequal(currentindex) = (age2 * (nostages - 1) * (nostages - 1) * totalages) + 
                  (time2n * (nostages - 1) * totalages) + (age3 * (nostages - 1)) + time3;
              }
            }
            
            if (time3 < nostages_nodead && time2n < nostages_nodead) {
              
              if (repmatrix((time3 + (nostages_nodead * time2n))) > 0) {
                
                // Now fecundity
                age3 = 0;
                currentindex = time3 + (age3 * nostages) + 
                  (time2n * nostages * totalages) + (age2 * nostages * nostages * totalages);
                
                stage3(currentindex) = newstageid(time3);
                stage2n(currentindex) = newstageid(time2n);
                stage2o(currentindex) = newstageid(time2n);
                stage1(currentindex) = 0;
                size3(currentindex) = binsizectr(time3);
                size2n(currentindex) = binsizectr(time2n);
                size2o(currentindex) = binsizectr(time2n);
                size1(currentindex) = 0;
                obs3(currentindex) = obsstatus(time3);
                obs2n(currentindex) = obsstatus(time2n);
                obs2o(currentindex) = obsstatus(time2n);
                obs1(currentindex) = 0;
                rep3(currentindex) = repstatus(time3);
                rep2n(currentindex) = repstatus(time2n);
                rep2o(currentindex) = repstatus(time2n);
                rep1(currentindex) = 0;
                mat3(currentindex) = matstatus(time3);
                mat2n(currentindex) = matstatus(time2n);
                mat2o(currentindex) = matstatus(time2n);
                mat1(currentindex) = 0;
                imm3(currentindex) = immstatus(time3);
                imm2n(currentindex) = immstatus(time2n);
                imm2o(currentindex) = immstatus(time2n);
                imm1(currentindex) = 0;
                
                if (rep2n(currentindex) > 0 && time3 < nostages_nodead && time2n < nostages_nodead) {
                  repentry3(currentindex) = repmatrix((time3 + (nostages_nodead * time2n)));
                } else repentry3(currentindex) = 0;
                
                indata3(currentindex) = indata(time3);
                indata2n(currentindex) = indata(time2n);
                indata2o(currentindex) = indata(time2n);
                indata1(currentindex) = 0;
                binwidth(currentindex) = binsizewidth(time3);
                minage3(currentindex) = minage(time3);
                minage2(currentindex) = minage(time2n);
                maxage3(currentindex) = maxage(time3);
                maxage2(currentindex) = maxage(time2n);
                actualage(currentindex) = age2;
                
                // The next indexer includes the following order: (1st # of age blocks) + (1st # of stage cols) +
                // (1st # of age rows) + stage in time 3
                index321(currentindex) = currentindex;
                indatalong(currentindex) = 1;
                
                // This section identifies elements with non-zero entries by their element number in the final matrix
                if (alive(time2n) == 1 && alive(time3) == 1) {
                  if (age2 >= minage2(currentindex) && age2 <= maxage2(currentindex)) { 
                    
                    // Fecundity transitions
                    aliveequal(currentindex) = (age2 * (nostages - 1) * (nostages - 1) * totalages) + 
                      (time2n * (nostages - 1) * totalages) + (age3 * (nostages - 1)) + time3;
                  }
                } // if statement leading to aliveequal assignment
              } // if statement yielding fecundity estimation
            } // if statement checking time3 and time2n
          } // time3 loop
        } // time2n loop
      } else if (cont == 1) { // This is the self-loop on the final age if the organism can live past the final age
        for (int time2n = 0; time2n < nostages; time2n++) {
          for (int time3 = 0; time3 < nostages; time3++) {
            
            // First survival
            age3 = age2;
            currentindex = time3 + (age3 * nostages) + 
              (time2n * nostages * totalages) + (age2 * nostages * nostages * totalages);
            
            stage3(currentindex) = newstageid(time3);
            stage2n(currentindex) = newstageid(time2n);
            stage2o(currentindex) = newstageid(time2n);
            stage1(currentindex) = 0;
            size3(currentindex) = binsizectr(time3);
            size2n(currentindex) = binsizectr(time2n);
            size2o(currentindex) = binsizectr(time2n);
            size1(currentindex) = 0;
            obs3(currentindex) = obsstatus(time3);
            obs2n(currentindex) = obsstatus(time2n);
            obs2o(currentindex) = obsstatus(time2n);
            obs1(currentindex) = 0;
            rep3(currentindex) = repstatus(time3);
            rep2n(currentindex) = repstatus(time2n);
            rep2o(currentindex) = repstatus(time2n);
            rep1(currentindex) = 0;
            mat3(currentindex) = matstatus(time3);
            mat2n(currentindex) = matstatus(time2n);
            mat2o(currentindex) = matstatus(time2n);
            mat1(currentindex) = 0;
            imm3(currentindex) = immstatus(time3);
            imm2n(currentindex) = immstatus(time2n);
            imm2o(currentindex) = immstatus(time2n);
            imm1(currentindex) = 0;
            
            repentry3(currentindex) = 0;
            
            indata3(currentindex) = indata(time3);
            indata2n(currentindex) = indata(time2n);
            indata2o(currentindex) = indata(time2n);
            indata1(currentindex) = 0;
            binwidth(currentindex) = binsizewidth(time3);
            minage3(currentindex) = minage(time3);
            minage2(currentindex) = minage(time2n);
            maxage3(currentindex) = maxage(time3);
            maxage2(currentindex) = maxage(time2n);
            actualage(currentindex) = age2;
            
            // The next indexer includes the following order: (1st # of age blocks) + (1st # of stage cols) +
            // (1st # of age rows) + stage in time 3
            index321(currentindex) = currentindex;
            indatalong(currentindex) = 1;
            
            // This section identifies elements with non-zero entries by their element number in the final matrix
            if (alive(time2n) == 1 && alive(time3) == 1) {
              if (age2 >= minage2(currentindex) && age2 < maxage3(currentindex)) { 
                
                // Survival transitions
                aliveequal(currentindex) = (age2 * (nostages - 1) * (nostages - 1) * totalages) + 
                  (time2n * (nostages - 1) * totalages) + (age3 * (nostages - 1)) + time3;
              }
            }
            
            if (time3 < nostages_nodead && time2n < nostages_nodead) {
            
              if (repmatrix((time3 + (nostages_nodead * time2n))) > 0) {
                
                // Now fecundity
                age3 = 0;
                currentindex = time3 + (age3 * nostages) + 
                  (time2n * nostages * totalages) + (age2 * nostages * nostages * totalages);
                
                stage3(currentindex) = newstageid(time3);
                stage2n(currentindex) = newstageid(time2n);
                stage2o(currentindex) = newstageid(time2n);
                stage1(currentindex) = 0;
                size3(currentindex) = binsizectr(time3);
                size2n(currentindex) = binsizectr(time2n);
                size2o(currentindex) = binsizectr(time2n);
                size1(currentindex) = 0;
                obs3(currentindex) = obsstatus(time3);
                obs2n(currentindex) = obsstatus(time2n);
                obs2o(currentindex) = obsstatus(time2n);
                obs1(currentindex) = 0;
                rep3(currentindex) = repstatus(time3);
                rep2n(currentindex) = repstatus(time2n);
                rep2o(currentindex) = repstatus(time2n);
                rep1(currentindex) = 0;
                mat3(currentindex) = matstatus(time3);
                mat2n(currentindex) = matstatus(time2n);
                mat2o(currentindex) = matstatus(time2n);
                mat1(currentindex) = 0;
                imm3(currentindex) = immstatus(time3);
                imm2n(currentindex) = immstatus(time2n);
                imm2o(currentindex) = immstatus(time2n);
                imm1(currentindex) = 0;
                
                if (rep2n(currentindex) == 1) {
                  repentry3(currentindex) = repmatrix((time3 + (nostages_nodead * time2n)));
                } else repentry3(currentindex) = 0;
                
                indata3(currentindex) = indata(time3);
                indata2n(currentindex) = indata(time2n);
                indata2o(currentindex) = indata(time2n);
                indata1(currentindex) = 0;
                binwidth(currentindex) = binsizewidth(time3);
                minage3(currentindex) = minage(time3);
                minage2(currentindex) = minage(time2n);
                maxage3(currentindex) = maxage(time3);
                maxage2(currentindex) = maxage(time2n);
                actualage(currentindex) = age2;
                
                // The next indexer includes the following order: (1st # of age blocks) + (1st # of stage cols) +
                // (1st # of age rows) + stage in time 3
                index321(currentindex) = currentindex;
                indatalong(currentindex) = 1;
                
                // This section identifies elements with non-zero entries by their element number in the final matrix
                if (alive(time2n) == 1 && alive(time3) == 1) {
                  if (age2 >= minage2(currentindex) && age2 <= maxage2(currentindex)) { 
                    
                    // Fecundity transitions
                    aliveequal(currentindex) = (age2 * (nostages - 1) * (nostages - 1) * totalages) + 
                      (time2n * (nostages - 1) * totalages) + (age3 * (nostages - 1)) + time3;
                  }
                } // if statement leading to aliveequal assignment
              } // if statement yielding fecundity estimation
            } // if statement checking time3 and time2n
          } // time3 loop
        } // time2n loop
      }// if-else statement
    } // age2 loop
    
    if (ovrows > 1 || ovconvtype(0) != -1) {
      asadditions = ovreplace(index321, ovindexold321, ovindexnew321, ovconvtypeage,
        ovnew3, ovnewgivenrate, ovnewmultiplier);
      
      ovgivent = asadditions.col(0);
      ovestt = asadditions.col(1);
      ovgivenf = asadditions.col(2);
      ovestf = asadditions.col(3);
      
      ovrepentry = asadditions.col(4);
      
      arma::uvec workedupindex = find(ovrepentry > 0);
      int changedreps = workedupindex.n_elem;
      
      if (changedreps > 0) {
        for (int i = 0; i < changedreps; i++) {
          repentry3(workedupindex(i)) = ovrepentry(workedupindex(i));
        }
      }
    } // ovreplace if statement
  } // Age by stage loop (style == 2)
  
  DataFrame df0 = DataFrame::create(Named("stage3") = stage3,
    _["stage2n"] = stage2n, _["stage2o"] = stage2o, _["stage1"] = stage1,
    _["size3"] = size3, _["size2n"] = size2n, _["size2o"] = size2o,
    _["size1"] = size1, _["obs3"] = obs3, _["obs2n"] = obs2n,
    _["obs2o"] = obs2o, _["obs1"] = obs1, _["rep3"] = rep3, _["rep2n"] = rep2n,
    _["rep2o"] = rep2o, _["rep1"] = rep1, _["mat3"] = mat3, _["mat2n"] = mat2n, 
    _["mat2o"] = mat2o, _["mat1"] = mat1);
  
  DataFrame df1 = DataFrame::create(Named("imm3") = imm3, _["imm2n"] = imm2n,
    _["imm2o"] = imm2o, _["imm1"] = imm1, _["repentry3"] = repentry3,
    _["indata3"] = indata3, _["indata2n"] = indata2n, _["indata2o"] = indata2o,
    _["indata1"] = indata1, _["binwidth"] = binwidth, _["minage3"] = minage3,
    _["minage2"] = minage2, _["maxage3"] = maxage3, _["maxage2"] = maxage2,
    _["actualage"] = actualage, _["indata"] = indatalong, _["ovgiven_t"] = ovgivent,
    _["ovest_t"] = ovestt, _["ovgiven_f"] = ovgivenf, _["ovest_f"] = ovestf);
  
  DataFrame df2 = DataFrame::create(Named("aliveandequal") = aliveequal, 
    _["index321"] = index321, _["index21"] = index21);
  
  return Rcpp::List::create(Named("a") = df0, Named("b") = df1, Named("c") = df2);
}

