#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Re-index Projection Matrix On Basis of Overwrite Table
//' 
//' This function takes matrix indices provided by functions \code{\link{rlefko3}()},
//' \code{\link{rlefko2}()}, \code{\link{flefko3}()}, and \code{\link{flefko2}()} and updates
//' them with information provided in the overwrite table used as input in that
//' function.
//' 
//' @param allst321 Vector containing the original element-by-element matrix index.
//' @param idx321old Vector containing the indices of matrix elements to be updated.
//' @param idx321new Vector containing the replacement matrix element indices.
//' @param convtype Vector denoting survival-transition (1) or fecundity (2).
//' @param eststag3 Vector of new stages in time \emph{t}+1.
//' @param gvnrate Vector of replacement transition values.
//' 
//' @return Vector of updated matrix indices and, where appropriate, replacement
//' matrix element values.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat ovreplace(arma::vec allst321, arma::vec idx321old, arma::vec idx321new, arma::vec convtype,
                    arma::vec eststag3, arma::vec gvnrate) {
  int n = idx321new.n_elem;
  
  arma::mat replacements(allst321.n_elem, 4);
  replacements.fill(-1);
  
  for (int i = 0; i < n; i++) {
    arma::uvec correctplace = find(allst321 == idx321old[i]);
    
    int m = correctplace.n_elem; 
    
    for (int j = 0; j < m; j++) {
      if (convtype[i] == 1) {
        if (gvnrate[i] != -1) {replacements(correctplace[j], 1) = gvnrate[i];}
        if (eststag3[i] != -1) {replacements(correctplace[j], 0) = idx321new[i];}
      }
      
      if (convtype[i] == 2) {
        if (gvnrate[i] != -1) {replacements(correctplace[j], 3) = gvnrate[i];}
        if (eststag3[i] != -1) {replacements(correctplace[j], 2) = idx321new[i];}
      }
    }
  }
  
  return replacements;
}

//' Make Horizontal Data Frame Vertical
//' 
//' Function \code{pfj()} powers the R function \code{\link{verticalize3}()}, creating
//' the vertical structure and rearranging the data in that shape.
//' 
//' @param data The horizontal data file.
//' @param stageframe The stageframe object identifying the life history model
//' being operationalized.
//' @param noyears The number of years or observation periods in the dataset.
//' @param firstyear The first year or time of observation.
//' @param popidcol Column number corresponding to the identity of the population 
//' for each individual.
//' @param patchidcol Column number corresponding to the identity of the patch for
//' each individual.
//' @param individcol Column number corresponding to the identity of each 
//' individual.
//' @param blocksize The number of variables corresponding to each time step in 
//' the input dataset designated in \code{data}.
//' @param xcol Column number corresponding to the x coordinate of each individual
//' in Cartesian space.
//' @param ycol Column number corresponding to the y coordinate of each individual
//' in Cartesian space.
//' @param juvcol Column number that marks individuals in immature stages within
//' the dataset.
//' @param size1col Column number corresponding to the first or main size variable 
//' associated with the first year or observation time in the dataset.
//' @param size2col Column number corresponding to the second size variable
//' associated with the first year or observation time in the dataset.
//' @param size3col Column number corresponding to the third size variable
//' associated with the first year or observation time in the dataset.
//' @param repstr1col Column number corresponding to the main variable coding the
//' production of reproductive structures associated with the first year or 
//' observation period in the input dataset.
//' @param repstr2col Column number corresponding to a secone variable coding the
//' production of reproductive structures associated with the first year or 
//' observation period in the input dataset.
//' @param fec1col Column number corresponding to the main variable coding for
//' fecundity associated with the first year or observation period in the dataset.
//' @param fec2col Column number corresponding to a second variable coding for
//' fecundity associated with the first year or observation period in the dataset.
//' @param alive1col Column number that details whether an individual is alive at
//' a given time.
//' @param dead1col Column number that details whether an individual is dead at
//' a given time.
//' @param obs1col Column number that details whether an individual is in an
//' observable stage at a given time.
//' @param nonobs1col Column number that details whether an individual is in an
//' unobservable stage at a given time.
//' @param censorcol Column number corresponding to the first entry of a censor 
//' variable.
//' @param stagecol Column number corresponding to the first entry of a column
//' designating stages.
//' @param repstrrel This is a scalar modifier for that makes the variable in
//' \code{repstr2col} equivalent to \code{repstr1col}.
//' @param fecrel This is a scalar modifier for that makes the variable in
//' \code{fec2col} equivalent to \code{fec1col}.
//' @param NAas0 If TRUE, then all NA entries for size and fecundity variables
//' will be set to 0.
//' @param NRasRep If TRUE, then will treat non-reproductive but mature
//' individuals as reproductive during stage assignment.
//' @param stassign A logical value indicating whether to assign stages.
//' @param stszcol Column number describing which size variable to use in stage 
//' estimation.
//' @param stagenum The number of stages in the dataset.
//' @param censorkeep The value of the censoring variable identifying data
//' that should be included in analysis. Defaults to 0, but may take any value
//' including NA.
//' @param censbool A logical variable determining whether NA denotes the value
//' of the censoring variable identifying data to keep.
//' 
//' @return The output is currently a 7 element list, where each element is a
//' data frame with the same number of rows.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List pfj(DataFrame data, DataFrame stageframe, int noyears, int firstyear, int popidcol, 
         int patchidcol, int individcol, int blocksize, int xcol, int ycol, int juvcol, 
         int size1col, int size2col, int size3col, int repstr1col, int repstr2col, 
         int fec1col, int fec2col, int alive1col, int dead1col, int obs1col, 
         int nonobs1col, int censorcol, int stagecol, double repstrrel, double fecrel,
         bool NAas0, bool NRasRep, bool stassign, int stszcol, int stagenum, bool censbool) {
  
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
  
  Rcpp::StringVector sfname(stagenum);
  Rcpp::NumericVector sfszmin(stagenum);
  Rcpp::NumericVector sfszmax(stagenum);
  Rcpp::NumericVector obsstat(stagenum);
  Rcpp::NumericVector repstat(stagenum);
  Rcpp::NumericVector matstat(stagenum);
  arma::vec sfszminarma(stagenum);
  arma::vec sfszmaxarma(stagenum);
  arma::vec obsstatarma(stagenum);
  arma::vec repstatarma(stagenum);
  arma::vec matstatarma(stagenum);
  
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
  Rcpp::NumericVector negonevec (noindivs);
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
  Rcpp::NumericVector fecb3  (ndflength);
  Rcpp::NumericVector fecb10 (ndflength);
  Rcpp::NumericVector fecb20 (ndflength);
  Rcpp::NumericVector fecb30  (ndflength);
  
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
  
  if (stassign) {
    sfname = stageframe[0];
    repstat = stageframe[2];
    obsstat = stageframe[3];
    matstat = stageframe[6];
    sfszmin = stageframe[9];
    sfszmax = stageframe[10];
    
    repstatarma = repstat;
    obsstatarma = obsstat;
    matstatarma = matstat;
    sfszminarma = sfszmin;
    sfszmaxarma = sfszmax;
  }
  
  for (int i = 0; i < noindivs; i++) {
    firstseenx[i] = 0;
    lastseenx[i] = 0;
    
    for (int k = 0; k < noyears; k++) {
      temp1 = data[size1col + (k * blocksize)];
      if (size2col != -1 ) {temp2 = data[size2col + (k * blocksize)];} else {temp2 = zerovec;}
      if (size3col != -1 ) {temp3 = data[size3col + (k * blocksize)];} else {temp3 = zerovec;}
      if (repstr1col != -1 ) {temp4 = data[repstr1col + (k * blocksize)];} else {temp4 = zerovec;}
      if (repstr2col != -1 ) {temp5 = data[repstr2col + (k * blocksize)];} else {temp5 = zerovec;}
    
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
    
    if (j == 0) {
      
      sizea1x = zerovec;
      sizea2x = data[size1col];
      sizea3x = data[size1col + blocksize];
      
      if (stagecol > -1) {
        stage1x.fill("NotAlive");
        stage2x = data[stagecol];
        stage3x = data[stagecol + blocksize];
      }
      
      if (xcol > -1) {
        xpos1x = zerovec;
        xpos2x = data[xcol];
        xpos3x = data[xcol + blocksize];
      }
      
      if (ycol > -1) {
        ypos1x = zerovec;
        ypos2x = data[ycol];
        ypos3x = data[ycol + blocksize];
      }
      
      if (repstr1col > -1) {
        repstra1x = zerovec;
        repstra2x = data[repstr1col];
        repstra3x = data[repstr1col + blocksize];
      }
      
      if (fec1col > -1) {
        feca1x = zerovec;
        feca2x = data[fec1col];
        feca3x = data[fec1col + blocksize];
      }
      
      if (size2col > -1) {
        sizeb1x = zerovec;
        sizeb2x = data[size2col];
        sizeb3x = data[size2col + blocksize];
      }
      
      if (repstr2col > -1) {
        repstrb1x = zerovec;
        repstrb2x = data[repstr2col];
        repstrb3x = data[repstr2col + blocksize];
      }
      
      if (fec2col > -1) {
        fecb1x = zerovec;
        fecb2x = data[fec2col];
        fecb3x = data[fec2col + blocksize];
      }
      
      if (size3col > -1) {
        sizec1x = zerovec;
        sizec2x = data[size3col];
        sizec3x = data[size3col + blocksize];
      }
      
      if (censorcol > -1) {
        censor1x = zerovec;
        censor2x = data[censorcol];
        censor3x = data[censorcol + blocksize];
      }
      
      if (alive1col > -1) {
        alivegiven1x = zerovec;
        alivegiven2x = data[alive1col];
        alivegiven3x = data[alive1col + blocksize];
      }
      
      if (dead1col > -1) {
        deadgiven1x = zerovec;
        deadgiven2x = data[dead1col];
        deadgiven3x = data[dead1col + blocksize];
      }
      
      if (obs1col > -1) {
        obsgiven1x = zerovec;
        obsgiven2x = data[obs1col];
        obsgiven3x = data[obs1col + blocksize];
      }
      
      if (nonobs1col > -1) {
        nonobsgiven1x = negonevec;
        nonobsgiven2x = data[nonobs1col];
        nonobsgiven3x = data[nonobs1col + blocksize];
      }
      
      if (juvcol > -1) {
        juvgiven1x = zerovec;
        juvgiven2x = data[juvcol];
        juvgiven3x = data[juvcol + blocksize];
      }
    } else if (j == (noyears - 1)) {
      
      sizea1x = data[size1col + ((j - 1) * blocksize)];
      sizea2x = data[size1col + (j * blocksize)];
      sizea3x = zerovec;
      
      if (stagecol > -1) {
        stage1x = data[stagecol + ((j - 1) * blocksize)];
        stage2x = data[stagecol + (j * blocksize)];
        stage3x.fill("Unknown");
      }
      
      if (xcol > -1) {
        xpos1x = data[xcol + ((j - 1) * blocksize)];
        xpos2x = data[xcol + (j * blocksize)];
        xpos3x = zerovec;
      }
      
      if (ycol > -1) {
        ypos1x = data[ycol + ((j - 1) * blocksize)];
        ypos2x = data[ycol + (j * blocksize)];
        ypos3x = zerovec;
      }
      
      if (repstr1col > -1) {
        repstra1x = data[repstr1col + ((j - 1) * blocksize)];
        repstra2x = data[repstr1col + (j * blocksize)];
        repstra3x = zerovec;
      }
      
      if (fec1col > -1) {
        feca1x = data[fec1col + ((j - 1) * blocksize)];
        feca2x = data[fec1col + (j * blocksize)];
        feca3x = zerovec;
      }
      
      if (size2col > -1) {
        sizeb1x = data[size2col + ((j - 1) * blocksize)];
        sizeb2x = data[size2col + (j * blocksize)];
        sizeb3x = zerovec;
      }
      
      if (repstr2col > -1) {
        repstrb1x = data[repstr2col + ((j - 1) * blocksize)];
        repstrb2x = data[repstr2col + (j * blocksize)];
        repstrb3x = zerovec;
      }
      
      if (fec2col > -1) {
        fecb1x = data[fec2col + ((j - 1) * blocksize)];
        fecb2x = data[fec2col + (j * blocksize)];
        fecb3x = zerovec;
      }
      
      if (size3col > -1) {
        sizec1x = data[size3col + ((j - 1) * blocksize)];
        sizec2x = data[size3col + (j * blocksize)];
        sizec3x = zerovec;
      }
      
      if (censorcol > -1) {
        censor1x = data[censorcol + ((j - 1) * blocksize)];
        censor2x = data[censorcol + (j * blocksize)];
        censor3x = zerovec;
      }
      
      if (alive1col > -1) {
        alivegiven1x = data[alive1col + ((j - 1) * blocksize)];
        alivegiven2x = data[alive1col + (j * blocksize)];
        alivegiven3x = zerovec;
      }
      
      if (dead1col > -1) {
        deadgiven1x = data[dead1col + ((j - 1) * blocksize)];
        deadgiven2x = data[dead1col + (j * blocksize)];
        deadgiven3x = zerovec;
      }
      
      if (obs1col > -1) {
        obsgiven1x = data[obs1col + ((j - 1) * blocksize)];
        obsgiven2x = data[obs1col + (j * blocksize)];
        obsgiven3x = zerovec;
      }
      
      if (nonobs1col > -1) {
        nonobsgiven1x = data[nonobs1col + ((j - 1) * blocksize)];
        nonobsgiven2x = data[nonobs1col + (j * blocksize)];
        nonobsgiven3x = negonevec;
      }
      
      if (juvcol > -1) {
        juvgiven1x = data[juvcol + ((j - 1) * blocksize)];
        juvgiven2x = data[juvcol + (j * blocksize)];
        juvgiven3x = zerovec;
      }
    } else {
      
      sizea1x = data[size1col + ((j - 1) * blocksize)];
      sizea2x = data[size1col + (j * blocksize)];
      sizea3x = data[size1col + ((j + 1) * blocksize)];
      
      if (stagecol > -1) {
        stage1x = data[stagecol + ((j - 1) * blocksize)];
        stage2x = data[stagecol + (j * blocksize)];
        stage3x = data[stagecol + ((j + 1) * blocksize)];
      }
      
      if (xcol > -1) {
        xpos1x = data[xcol + ((j - 1) * blocksize)];
        xpos2x = data[xcol + (j * blocksize)];
        xpos3x = data[xcol + ((j + 1) * blocksize)];
      }
      
      if (ycol > -1) {
        ypos1x = data[ycol + ((j - 1) * blocksize)];
        ypos2x = data[ycol + (j * blocksize)];
        ypos3x = data[ycol + ((j + 1) * blocksize)];
      }
      
      if (repstr1col > -1) {
        repstra1x = data[repstr1col + ((j - 1) * blocksize)];
        repstra2x = data[repstr1col + (j * blocksize)];
        repstra3x = data[repstr1col + ((j + 1) * blocksize)];
      }
      
      if (fec1col > -1) {
        feca1x = data[fec1col + ((j - 1) * blocksize)];
        feca2x = data[fec1col + (j * blocksize)];
        feca3x = data[fec1col + ((j + 1) * blocksize)];
      }
      
      if (size2col > -1) {
        sizeb1x = data[size2col + ((j - 1) * blocksize)];
        sizeb2x = data[size2col + (j * blocksize)];
        sizeb3x = data[size2col + ((j + 1) * blocksize)];
      }
      
      if (repstr2col > -1) {
        repstrb1x = data[repstr2col + ((j - 1) * blocksize)];
        repstrb2x = data[repstr2col + (j * blocksize)];
        repstrb3x = data[repstr2col + ((j + 1) * blocksize)];
      }
      
      if (fec2col > -1) {
        fecb1x = data[fec2col + ((j - 1) * blocksize)];
        fecb2x = data[fec2col + (j * blocksize)];
        fecb3x = data[fec2col + ((j + 1) * blocksize)];
      }
      
      if (size3col > -1) {
        sizec1x = data[size3col + ((j - 1) * blocksize)];
        sizec2x = data[size3col + (j * blocksize)];
        sizec3x = data[size3col + ((j + 1) * blocksize)];
      }
      
      if (censorcol > -1) {
        censor1x = data[censorcol + ((j - 1) * blocksize)];
        censor2x = data[censorcol + (j * blocksize)];
        censor3x = data[censorcol + ((j + 1) * blocksize)];
      }
      
      if (alive1col > -1) {
        alivegiven1x = data[alive1col + ((j - 1) * blocksize)];
        alivegiven2x = data[alive1col + (j * blocksize)];
        alivegiven3x = data[alive1col + ((j + 1) * blocksize)];
      }
      
      if (dead1col > -1) {
        deadgiven1x = data[dead1col + ((j - 1) * blocksize)];
        deadgiven2x = data[dead1col + (j * blocksize)];
        deadgiven3x = data[dead1col + ((j + 1) * blocksize)];
      }
      
      if (obs1col > -1) {
        obsgiven1x = data[obs1col + ((j - 1) * blocksize)];
        obsgiven2x = data[obs1col + (j * blocksize)];
        obsgiven3x = data[obs1col + ((j + 1) * blocksize)];
      }
      
      if (nonobs1col > -1) {
        nonobsgiven1x = data[nonobs1col + ((j - 1) * blocksize)];
        nonobsgiven2x = data[nonobs1col + (j * blocksize)];
        nonobsgiven3x = data[nonobs1col + ((j + 1) * blocksize)];
      }
      
      if (juvcol > -1) {
        juvgiven1x = data[juvcol + ((j - 1) * blocksize)];
        juvgiven2x = data[juvcol + (j * blocksize)];
        juvgiven3x = data[juvcol + ((j + 1) * blocksize)];
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
      
      if (stagecol > -1) {
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
      
      if (addedsize1[(i + (j * noindivs))] > 0 || obsgiven1[(i + (j * noindivs))] > 0 || nonobssub1 == 0) {
        spryn1[(i + (j * noindivs))] = 1;
      }
      if (addedsize2[(i + (j * noindivs))] > 0 || obsgiven2[(i + (j * noindivs))] > 0 || nonobssub2 == 0) {
        spryn2[(i + (j * noindivs))] = 1;
      } 
      if (addedsize3[(i + (j * noindivs))] > 0 || obsgiven3[(i + (j * noindivs))] > 0 || nonobssub3 == 0) {
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
      if (stassign && stagecol == -1) {
        
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
        stagemini = find(sfszminarma < stagesize2);
        stagemaxi = find(sfszmaxarma >= stagesize2);
        stagerep = find(repstatarma == repyn2[(i + (j * noindivs))]);
        stagemat = find(matstatarma == matstat2[(i + (j * noindivs))]);
        stageobs = find(obsstatarma == spryn2[(i + (j * noindivs))]);
        
        cs1 = intersect(stagemini, stagemaxi);
        cs2 = intersect(stageobs, stagemat);
        
        if (NRasRep) {
          cs4 = intersect(cs1, cs2);
        } else {
          cs3 = intersect(stagerep, cs2);
          cs4 = intersect(cs1, cs3);
        }
        
        if (cs4.n_elem > 0 && alive2[(i + (j * noindivs))] == 1) {
          choicestage = cs4[0];
          stage2num[(i + (j * noindivs))] = choicestage + 1;
          
          stage2[(i + (j * noindivs))] = sfname[choicestage];
        } else {
          stage2[(i + (j * noindivs))] = "NotAlive";
        }
        
        // Stage 1
        stagemini = find(sfszminarma < stagesize1);
        stagemaxi = find(sfszmaxarma >= stagesize1);
        stagerep = find(repstatarma == repyn1[(i + (j * noindivs))]);
        stagemat = find(matstatarma == matstat1[(i + (j * noindivs))]);
        stageobs = find(obsstatarma == spryn1[(i + (j * noindivs))]);
        
        cs1 = intersect(stagemini, stagemaxi);
        cs2 = intersect(stageobs, stagemat);
        
        if (NRasRep) {
          cs4 = intersect(cs1, cs2);
        } else {
          cs3 = intersect(stagerep, cs2);
          cs4 = intersect(cs1, cs3);
        }
        
        if (cs4.n_elem > 0 && alive1[(i + (j * noindivs))] == 1) {
          choicestage = cs4[0];
          stage1num[(i + (j * noindivs))] = choicestage + 1;
          
          stage1[(i + (j * noindivs))] = sfname[choicestage];
        } else {
          stage1[(i + (j * noindivs))] = "NotAlive";
          matstat1[(i + (j * noindivs))] = 0;
        }
        
        // Stage 3
        stagemini = find(sfszminarma < stagesize3);
        stagemaxi = find(sfszmaxarma >= stagesize3);
        stagerep = find(repstatarma == repyn3[(i + (j * noindivs))]);
        stagemat = find(matstatarma == matstat3[(i + (j * noindivs))]);
        stageobs = find(obsstatarma == spryn3[(i + (j * noindivs))]);
        
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
          
          choicestage = cs4[0];
          stage3num[(i + (j * noindivs))] = choicestage + 1;
          
          stage3[(i + (j * noindivs))] = sfname[choicestage];
          
        } else if (alive3[(i + (j * noindivs))] != 1) {
          
          stage3[(i + (j * noindivs))] = "NotAlive";

        } else if (cs4.n_elem == 0) {
          
          stage3[(i + (j * noindivs))] = "NoMatch";
          
          Rcpp::warning("Some stages occurring in the dataset do not match any characteristics in the input stageframe.");
          
        } else if (cs4.n_elem > 1) {
          
          Rcpp::warning("Some stages in the input stageframe appear to have the same description. Please make sure that all stages included in the stageframe are defined with unique sets of characteristics.");
          
        } else {
          
          Rcpp::stop("Stage assignment error.");
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
  
  Rcpp::DataFrame df0 = DataFrame::create(Named("rowid") = rowid, _["popid"] = popid, _["patchid"] = patchid, 
                                          _["individ"] = individ,  _["year2"] = year2, _["firstseen"] = firstseen,
                                          _["lastseen"] = lastseen, _["obsage"] = obsage, _["obslifespan"] = obslifespan);
  
  Rcpp::DataFrame df1 = DataFrame::create(Named("xpos1") = xpos1, _["ypos1"] = ypos1, _["sizea1"] = sizea1,
                                          _["sizeb1"] = sizeb1, _["sizec1"] = sizec1, _["addedsize1"] = addedsize1,
                                          _["repstra1"] = repstra1, _["repstrb1"] = repstrb1, _["addedrepstr1"] = addedrepstr1,
                                          _["feca1"] = feca1, _["fecb1"] = fecb1, _["addedfec1"] = addedfec1,
                                          _["censor1"] = censor1, _["juvgiven1"] = juvgiven1);
  
  Rcpp::DataFrame df1a = DataFrame::create(Named("obsstatus1") = spryn1, _["repstatus1"] = repyn1, _["fecstatus1"] = fecyn1,
                                           _["matstatus1"] = matstat1, _["alive1"] = alive1, _["stage1"] = stage1,
                                           _["stage1index"] = stage1num);
  
  Rcpp::DataFrame df2 = DataFrame::create(Named("xpos2") = xpos2,  _["ypos2"] = ypos2,  _["sizea2"] = sizea2,
                                          _["sizeb2"] = sizeb2, _["sizec2"] = sizec2, _["addedsize2"] = addedsize2,
                                          _["repstra2"] = repstra2, _["repstrb2"] = repstrb2, _["addedrepstr2"] = addedrepstr2,
                                          _["feca2"] = feca2, _["fecb2"] = fecb2, _["addedfec2"] = addedfec2,
                                          _["censor2"] = censor2, _["juvgiven2"] = juvgiven2);
  
  Rcpp::DataFrame df2a = DataFrame::create(Named("obsstatus2") = spryn2, _["repstatus2"] = repyn2, _["fecstatus2"] = fecyn2,
                                           _["matstatus2"] = matstat2, _["alive2"] = alive2, _["stage2"] = stage2,
                                           _["stage2index"] = stage2num);
  
  Rcpp::DataFrame df3 = DataFrame::create(Named("xpos3") = xpos3, _["ypos3"] = ypos3, _["sizea3"] = sizea3, 
                                          _["sizeb3"] = sizeb3, _["sizec3"] = sizec3, _["addedsize3"] = addedsize3,
                                          _["repstra3"] = repstra3, _["repstrb3"] = repstrb3, _["addedrepstr3"] = addedrepstr3,
                                          _["feca3"] = feca3, _["fecb3"] = fecb3, _["addedfec3"] = addedfec3,
                                          _["censor3"] = censor3, _["juvgiven3"] = juvgiven3);
  
  Rcpp::DataFrame df3a = DataFrame::create(Named("obsstatus3") = spryn3, _["repstatus3"] = repyn3, _["fecstatus3"] = fecyn3,
                                           _["matstatus3"] = matstat3, _["alive3"] = alive3, _["stage3"] = stage3,
                                           _["stage3index"] = stage3num);
  
  
  return Rcpp::List::create(Named("a") = df0, Named("b") = df1, Named("c") = df1a, Named("d") = df2, Named("e") = df2a,
                                  Named("f") = df3, Named("g") = df3a);
}

//' Make Vertical Data Frame Historical
//' 
//' Function \code{jpf()} powers the R function \code{\link{historicalize3}()}, creating
//' the historical, vertical structure and rearranging the data in that shape.
//'
//' @param data The horizontal data file.
//' @param stageframe The stageframe object identifying the life history model
//' being operationalized.
//' @param popidcol Column number corresponding to the identity of the population 
//' for each individual.
//' @param patchidcol Column number corresponding to the identity of the patch for
//' each individual.
//' @param individcol Column number corresponding to the identity of each 
//' individual.
//' @param year2col Column number of year or time step in time \emph{t}.
//' @param year3col Column number of year or time step in time \emph{t}+1.
//' @param xcol Column number corresponding to the x coordinate of each individual
//' in Cartesian space.
//' @param ycol Column number corresponding to the y coordinate of each individual
//' in Cartesian space.
//' @param juv2col Column number coding for status as a juvenile in time \emph{t}.
//' @param juv3col Column number coding for status as a juvenile in time \emph{t}+1.
//' @param sizea2col Column number corresponding to the primary size variable in
//' time \emph{t}.
//' @param sizea3col Column number corresponding to the primary size variable in
//' time \emph{t}+1.
//' @param sizeb2col Column number corresponding to the secondary size variable in
//' time \emph{t}.
//' @param sizeb3col Column number corresponding to the secondary size variable in
//' time \emph{t}+1.
//' @param sizec2col Column number corresponding to the tertiary size variable in
//' time \emph{t}.
//' @param sizec3col Column number corresponding to the tertiary size variable in
//' time \emph{t}+1.
//' @param repstra2col Column number corresponding to the main variable coding the
//' production of reproductive structures, such as flowers, in time \emph{t}.
//' @param repstra3col Column number corresponding to the main variable coding the
//' production of reproductive structures, such as flowers, in time \emph{t}+1.
//' @param repstrb2col Column number corresponding to a second variable coding the
//' production of reproductive structures, such as flowers, in time
//' @param repstrb3col Column number corresponding to a second variable coding the
//' production of reproductive structures, such as flowers, in time \emph{t}+1.
//' @param feca2col Column number corresponding to the main variable coding for
//' fecundity in time \emph{t}.
//' @param feca3col Column number corresponding to the main variable coding for
//' fecundity in time \emph{t}+1.
//' @param fecb2col Column number corresponding to a second variable coding for
//' fecundity in time \emph{t}.
//' @param fecb3col Column number corresponding to a second variable coding for
//' fecundity in time \emph{t}+1.
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
//' @param stage2col Column number corresponding to life history stage in time \emph{t}.
//' @param stage3col Column number corresponding to life history stage in time \emph{t}+1.
//' @param censorcol Column number corresponding to a censor variable within the 
//' dataset.
//' @param NAas0 If TRUE, then all NA entries for size and fecundity variables
//' will be set to 0.
//' @param NRasRep If TRUE, then will treat non-reproductive but mature
//' individuals as reproductive during stage assignment.
//' @param stassign A logical value indicating whether to assign stages.
//' @param stszcol Column number describing which size variable to use in stage 
//' estimation.
//' @param stagenum The number of stages in the dataset.
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
               int patchidcol, int individcol, int year2col, int year3col, int xcol, int ycol,
               int juv2col, int juv3col, int sizea2col, int sizea3col, int sizeb2col, int sizeb3col,
               int sizec2col, int sizec3col, int repstra2col, int repstra3col, int repstrb2col,
               int repstrb3col, int feca2col, int feca3col, int fecb2col, int fecb3col, int alive2col,
               int alive3col, int dead2col, int dead3col, int obs2col, int obs3col, int nonobs2col,
               int nonobs3col, double repstrrel, double fecrel, int stage2col, int stage3col, 
               int censorcol, bool NAas0, bool NRasRep, bool stassign, int stszcol, int stagenum,
               bool censbool) {
  
  int norows = data.nrows();
  
  Rcpp::NumericVector zerovec (norows);
  Rcpp::NumericVector negonevec (norows);
  negonevec.fill(-1);
  
  Rcpp::StringVector individx(norows);
  individx = data[individcol];
  Rcpp::StringVector allindivs = unique(individx);
  int noindivs = allindivs.size();
  
  Rcpp::IntegerVector year2x(norows);
  Rcpp::IntegerVector year3x(norows);
  year2x = data[year2col];
  if (year3col != -1) {year3x = data[year3col];} else {year3x = zerovec;}
  Rcpp::IntegerVector yearall2x = sort_unique(year2x);
  int firstyear = min(yearall2x);
  int noyears = yearall2x.size();
  
  int ndflength = noyears * noindivs;
  int currentyear {0};
  int prevyear {0};
  int nextyear {0};
  int currentindiv {0};
  int ndfindex {0};
  int prevyrindex {0};
  int nextyrindex {0};
  double livcheck2 {0};
  double livcheck3 {0};
  
  Rcpp::StringVector sfname(stagenum);
  Rcpp::NumericVector sfszmin(stagenum);
  Rcpp::NumericVector sfszmax(stagenum);
  Rcpp::NumericVector obsstat(stagenum);
  Rcpp::NumericVector repstat(stagenum);
  Rcpp::NumericVector matstat(stagenum);
  arma::vec sfszminarma(stagenum);
  arma::vec sfszmaxarma(stagenum);
  arma::vec obsstatarma(stagenum);
  arma::vec repstatarma(stagenum);
  arma::vec matstatarma(stagenum);
  
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
  Rcpp::IntegerVector year2 (ndflength);
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
  Rcpp::NumericVector repstatus2 (ndflength);
  Rcpp::NumericVector fecstatus2 (ndflength);
  Rcpp::NumericVector obsstatus2 (ndflength);
  Rcpp::NumericVector juvgiven2 (ndflength);
  Rcpp::NumericVector matstat2 (ndflength);
  
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
  Rcpp::NumericVector repstatus3 (ndflength);
  Rcpp::NumericVector fecstatus3 (ndflength);
  Rcpp::NumericVector obsstatus3 (ndflength);
  Rcpp::NumericVector juvgiven3 (ndflength);
  Rcpp::NumericVector matstat3 (ndflength);
  
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
  Rcpp::NumericVector repstatus1 (ndflength);
  Rcpp::NumericVector fecstatus1 (ndflength);
  Rcpp::NumericVector obsstatus1 (ndflength);
  Rcpp::NumericVector juvgiven1 (ndflength);
  Rcpp::NumericVector matstat1 (ndflength);
  
  Rcpp::NumericVector firstseen (ndflength);
  Rcpp::NumericVector lastseen (ndflength);
  Rcpp::NumericVector obsage (ndflength);
  Rcpp::NumericVector obslifespan (ndflength);
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
  
  if (stassign) {
    sfname = stageframe[0];
    repstat = stageframe[2];
    obsstat = stageframe[3];
    matstat = stageframe[6];
    sfszmin = stageframe[9];
    sfszmax = stageframe[10];
    
    repstatarma = repstat;
    obsstatarma = obsstat;
    matstatarma = matstat;
    sfszminarma = sfszmin;
    sfszmaxarma = sfszmax;
  }
  
  // Initialize main loop, which focuses on rows and establishes state in time t
  for (int i = 0; i < norows; i++) {
    for (int j = 0; j < noyears; j++) {
      if (year2x[i] == yearall2x[j]) currentyear = j;
    }
    nextyear = currentyear + 1;
    prevyear = currentyear - 1;
    
    for (int k = 0; k < noindivs; k++) {
      if (individx[i] == allindivs[k]) currentindiv = k;
    }
    
    ndfindex = (noyears * currentindiv) + currentyear;
    
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
    
    if (censbool && censorcol != -1) { // Here we develop the censoring variable
      if (NumericVector::is_na(censor2x[i])) {
        censor2[ndfindex] = 0;
      } else {
        censor2[ndfindex] = 1;
      }
    } else if (censorcol != -1) {
      censor2[ndfindex] = censor2x[i];
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
    } else {stage2[ndfindex] = "NotAlive";}
    
    //Now we work on time t+1 in cases where t+1 columns are provided
    if (currentyear == (noyears - 1)) {
      if (censbool && censorcol != -1) { // Here we develop the censoring variable for the last time
        if (NumericVector::is_na(censor2x[i])) {
          censor3[ndfindex] = 0;
        } else {
          censor3[ndfindex] = 1;
        }
      } else if (censorcol != -1) {
        censor3[ndfindex] = censor2x[i];
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
      } else {stage3[ndfindex] = "NotAlive";}
    }
  }
  
  // Now a loop that establishes most states in time t+1 and t-1, and stages in all times
  for (int i = 0; i < ndflength; i++) {
    
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
    
    // Now the normal calculations
    for (int k = 0; k < noindivs; k++) {
      if (individ[i] == allindivs[k]) currentindiv = k;
    }
    
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
    } else if (stage2col != -1 && alive2[i] == 0) {
      stage2[i] = "NotAlive";
    }
    
    if (currentyear < (noyears - 1)) { 
      nextyrindex = (noyears * currentindiv) + (currentyear + 1);
      
      censor3[i] = censor2[nextyrindex];
      xpos3[i] = xpos2[nextyrindex];
      ypos3[i] = ypos3[nextyrindex];
      
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
      
      if (stage2col != -1 && alive3[i] == 1) {
        if (obsstatus3[i] == 0) {
          stage3[i] = "NotObserved";
        } else {stage3[i] = stage2[nextyrindex];}
      } else if (stage2col != -1 && alive3[i] == 0) {stage3[i] = "NotAlive";}
    }
    
    if (currentyear > 0  && year2[i] < (firstyear + noyears)) {
      prevyrindex = (noyears * currentindiv) + (currentyear - 1);
      
      censor1[i] = censor2[prevyrindex];
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
      
      if (fecadded1[i] > 0) {
        fecstatus1[i] = 1;
      } else {
        fecstatus1[i] = 0;
      }
      
      repstatus1[i] = repstatus2[prevyrindex];
      obsstatus1[i] = obsstatus2[prevyrindex];
      juvgiven1[i] = juvgiven2[prevyrindex];
      matstat1[i] = matstat2[prevyrindex];
      
      alive1[i] = alive2[prevyrindex];
      
      if (stage2col != -1 && alive1[i] == 1) {
        if (obsstatus1[i] == 0) {
          stage1[i] = "NotObserved";
        } else {stage1[i] = stage2[prevyrindex];}
      } else if (stage2col != -1 && alive1[i] == 0) {stage1[i] = "NotAlive";}
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
      stagemini = find(sfszminarma < stagesize2);
      stagemaxi = find(sfszmaxarma >= stagesize2);
      stagerep = find(repstatarma == repstatus2[i]);
      stagemat = find(matstatarma == matstat2[i]);
      stageobs = find(obsstatarma == obsstatus2[i]);
      
      cs1 = intersect(stagemini, stagemaxi);
      cs2 = intersect(stageobs, stagemat);
      
      if (NRasRep) {
        cs4 = intersect(cs1, cs2);
      } else {
        cs3 = intersect(stagerep, cs2);
        cs4 = intersect(cs1, cs3);
      }
      
      if (cs4.n_elem > 0 && alive2[i] == 1) {
        choicestage = cs4[0];
        stage2num[i] = choicestage + 1;
        
        stage2[i] = sfname[choicestage];
      } else {
        stage2[i] = "NotAlive";
      }
      
      // Stage 1
      stagemini = find(sfszminarma < stagesize1);
      stagemaxi = find(sfszmaxarma >= stagesize1);
      stagerep = find(repstatarma == repstatus1[i]);
      stagemat = find(matstatarma == matstat1[i]);
      stageobs = find(obsstatarma == obsstatus1[i]);
      
      cs1 = intersect(stagemini, stagemaxi);
      cs2 = intersect(stageobs, stagemat);
      
      if (NRasRep) {
        cs4 = intersect(cs1, cs2);
      } else {
        cs3 = intersect(stagerep, cs2);
        cs4 = intersect(cs1, cs3);
      }
      
      if (cs4.n_elem > 0 && alive1[i] == 1) {
        choicestage = cs4[0];
        stage1num[i] = choicestage + 1;
        
        stage1[i] = sfname[choicestage];
      } else {
        stage1[i] = "NotAlive";
        matstat1[i] = 0;
      }
      
      // Stage 3
      stagemini = find(sfszminarma < stagesize3);
      stagemaxi = find(sfszmaxarma >= stagesize3);
      stagerep = find(repstatarma == repstatus3[i]);
      stagemat = find(matstatarma == matstat3[i]);
      stageobs = find(obsstatarma == obsstatus3[i]);
      
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
        
        choicestage = cs4[0];
        stage3num[i] = choicestage + 1;
        
        stage3[i] = sfname[choicestage];
        
      } else if (alive3[i] != 1) {
        
        stage3[i] = "NotAlive";

      } else if (cs4.n_elem == 0) {
        
        stage3[i] = "NoMatch";
        
        Rcpp::warning("Some stages occurring in the dataset do not match any characteristics in the input stageframe.");
        
      } else if (cs4.n_elem > 1) {
        
        Rcpp::warning("Some stages in the input stageframe appear to have the same description. Please make sure that all stages included in the stageframe are defined with unique sets of characteristics.");
        
      } else {
        
        Rcpp::stop("Stage assignment error.");
      }
    } // stassign if statement
  } // i loop
  
  Rcpp::DataFrame df0 = DataFrame::create(Named("rowid") = rowid, _["popid"] = popid, _["patchid"] = patchid, 
                                          _["individ"] = individ,  _["year2"] = year2, _["firstseen"] = firstseen,
                                          _["lastseen"] = lastseen, _["obsage"] = obsage, _["obslifespan"] = obslifespan);
  
  Rcpp::DataFrame df1 = DataFrame::create(Named("xpos1") = xpos1, _["ypos1"] = ypos1, _["sizea1"] = sizea1, 
                                          _["sizeb1"] = sizeb1, _["sizec1"] = sizec1, _["sizeadded1"] = sizeadded1, 
                                            _["repstra1"] = repstra1, _["repstrb1"] = repstrb1, _["repstradded1"] = repstradded1, 
                                              _["feca1"] = feca1, _["fecb1"] = fecb1, _["fecaadded1"] = fecadded1, 
                                                _["censor1"] = censor1, _["juvgiven1"] = juvgiven1);
  Rcpp::DataFrame df1a = DataFrame::create(Named("obsstatus1") = obsstatus1, _["repstatus1"] = repstatus1, 
                                           _["fecstatus1"] = fecstatus1, _["matstatus1"] = matstat1, _["alive1"] = alive1, 
                                             _["stage1"] = stage1, _["stage1index"] = stage1num);
  
  Rcpp::DataFrame df2 = DataFrame::create(Named("xpos2") = xpos2, _["ypos2"] = ypos2, _["sizea2"] = sizea2, 
                                          _["sizeb2"] = sizeb2, _["sizec2"] = sizec2, _["sizeadded2"] = sizeadded2, 
                                            _["repstra2"] = repstra2, _["repstrb2"] = repstrb2, _["repstradded2"] = repstradded2, 
                                              _["feca2"] = feca2, _["fecb2"] = fecb2, _["fecaadded2"] = fecadded2, 
                                                _["censor2"] = censor2, _["juvgiven2"] = juvgiven2);
  Rcpp::DataFrame df2a = DataFrame::create(Named("obsstatus2") = obsstatus2, _["repstatus2"] = repstatus2, 
                                           _["fecstatus2"] = fecstatus2, _["matstatus2"] = matstat2, _["alive2"] = alive2,
                                             _["stage2"] = stage2, _["stage2index"] = stage2num);
  
  Rcpp::DataFrame df3 = DataFrame::create(Named("xpos3") = xpos3, _["ypos3"] = ypos3, _["sizea3"] = sizea3,
                                          _["sizeb3"] = sizeb3, _["sizec3"] = sizec3, _["sizeadded3"] = sizeadded3, 
                                            _["repstra3"] = repstra3, _["repstrb3"] = repstrb3, _["repstradded3"] = repstradded3,  
                                              _["feca3"] = feca3, _["fecb3"] = fecb3, _["fecaadded3"] = fecadded3, 
                                                _["censor3"] = censor3, _["juvgiven3"] = juvgiven3);
  Rcpp::DataFrame df3a = DataFrame::create(Named("obsstatus3") = obsstatus3, _["repstatus3"] = repstatus3, 
                                           _["fecstatus3"] = fecstatus3, _["matstatus3"] = matstat3, _["alive3"] = alive3, 
                                             _["stage3"] = stage3, _["stage3index"] = stage3num);
  
  return Rcpp::List::create(Named("a") = df0, Named("b") = df1, Named("c") = df1a, Named("d") = df2, Named("e") = df2a,
                                  Named("f") = df3, Named("g") = df3a);
}

