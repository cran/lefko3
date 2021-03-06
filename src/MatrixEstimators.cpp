#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Estimate All Elements of Raw Historical Matrix
//' 
//' Function \code{specialpatrolgroup} swiftly calculates matrix transitions in
//' raw historical matrices, and serves as the core workhorse function behind
//' \code{\link{rlefko3}()}.
//' 
//' @param sge9l The Allstages data frame developed for rlefko3() covering stage
//' pairs across times \emph{t}+1, \emph{t} and \emph{t}-1. Generally termed
//' \code{stageexpansion9}.
//' @param sge3 The data frame covering all stages in times \emph{t} and
//' \emph{t}-1. Generally termed \code{stageexpansion3}.
//' @param MainData The demographic dataset modified to hold \code{usedfec}
//' columns.
//' @param StageFrame The full stageframe for the analysis.
//' @param repmatrix The modified repmatrix used in the course of computation.
//' This is used particularly when deVries-format hMPMs are desired.
//' @param format Indicates whether to output Ehrlen-format hMPMs (1) or
//' deVries-format hMPMs (2).
//' @param err_switch If set to 1, then will also output probsrates and
//' stage2fec.
//' 
//' @return List of three matrices, including the survival-transition (U)
//' matrix, the fecundity matrix (F), and the sum (A) matrix, with A first.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List specialpatrolgroup(DataFrame sge9l, DataFrame sge3, DataFrame MainData,
  DataFrame StageFrame, int format, int err_switch) {
  
  arma::vec sge9stage3 = sge9l["a.stage3"];
  arma::vec sge9fec32 = sge9l["b.repentry3"];
  arma::vec sge9rep2 = sge9l["a.rep2o"];
  arma::vec sge9indata32 = sge9l["b.indata"];
  arma::vec sge9ovgivent = sge9l["b.ovgiven_t"];
  arma::vec sge9ovgivenf = sge9l["b.ovgiven_f"];
  arma::vec sge9ovestt = sge9l["b.ovest_t"];
  arma::vec sge9ovestf = sge9l["b.ovest_f"];
  arma::vec sge9index321 = sge9l["c.index321"];
  arma::vec sge9index21 = sge9l["c.index21"]; // This is sge92index - not sure if this is needed
  arma::vec aliveandequal = sge9l["c.aliveandequal"];
  
  arma::vec sge3rep2 = sge3["rep2n"];
  arma::vec sge3fec32 = sge3["fec32n"];
  arma::vec sge3index21 = sge3["index21"];
  arma::vec sge3stage2n = sge3["stage2n"];
  arma::vec sge3stage3 = sge3["stage3"];
  
  arma::vec dataindex321 = MainData["index321"];
  arma::vec dataindex21 = MainData["pairindex21"];
  arma::vec dataalive3 = MainData["alive3"];
  arma::vec datausedfec2 = MainData["usedfec2"];
  arma::vec dataindex3 = MainData["index3"];
  arma::vec dataindex2 = MainData["index2"];
  arma::vec dataindex1 = MainData["index1"];

  arma::vec sfsizes = StageFrame["bin_size_ctr"];
  int nostages = sfsizes.n_elem;
  
  int n = dataindex321.n_elem;
  int no2stages = nostages - 1;
  int noelems = sge9index321.n_elem;
  
  if (format == 2) {
    no2stages = no2stages - 1;
  }
  
  int matrixdim = (nostages - 1) * no2stages;
  
  arma::vec probsrates0(noelems); // 1st vec = # indivs (3 trans)
  arma::vec probsrates1(noelems); // 2nd vec = total indivs for pair stage
  arma::vec probsrates2(noelems); // 3rd vec = total indivs alive for pair stage
  arma::vec probsrates3(noelems); // 4th vec = total fec for pair stage
  probsrates0.zeros();
  probsrates1.zeros();
  probsrates2.zeros();
  probsrates3.zeros();
  
  arma::mat stage2fec(sge3index21.n_elem, 3); //1st col = # indivs total, 2nd col = no indivs alive, 3rd col = sum fec
  stage2fec.zeros();
  
  // These next structures develop the prior stage
  arma::vec probsrates0p(noelems); // 1st vec = # indivs (3 trans)
  arma::vec probsrates1p(noelems); // 2nd vec = total indivs for pair stage
  arma::vec probsrates2p(noelems); // 3rd vec = total indivs alive for pair stage
  arma::vec probsrates3p(noelems); // 4th vec = total fec for pair stage
  probsrates0p.zeros();
  probsrates1p.zeros();
  probsrates2p.zeros();
  probsrates3p.zeros();
  
  arma::mat stage2fecp(sge3index21.n_elem, 3); //1st col = # indivs total, 2nd col = no indivs alive, 3rd col = sum fec
  stage2fecp.zeros();
  
  // The final matrices, though empty
  arma::mat tmatrix(matrixdim, matrixdim); // Main output U matrix
  arma::mat fmatrix(matrixdim, matrixdim); // Main output F matrix
  tmatrix.zeros();
  fmatrix.zeros();
  
  arma::uvec all_repentries = find(sge9fec32 > 0);
  arma::vec all_entry_stages = arma::unique(sge9stage3(all_repentries));
  int aes_count = all_entry_stages.n_elem;
  
  arma::mat dataindex321_prior(n, aes_count);
  dataindex321_prior.fill(-1);
  
  // This section creates an alternative index for use in fecundity calculations under deVries format
  if (format == 2) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < aes_count; j++) {
        dataindex321_prior(i, j) = (all_entry_stages(j) - 1) + ((nostages - 2) * nostages) + 
          ((dataindex2(i) - 1) * nostages * nostages) + 
          ((dataindex1(i) - 1) * nostages * nostages * nostages);
      }
    }
  }
  
  // This main loop counts individuals going through transitions and sums their
  // fecundities, and then adds that info to the 3-trans and 2-trans tables
  for (int i = 0; i < n; i++) { 
    arma::uvec choiceelement = find(sge9index321 == dataindex321(i)); // Added this now
    
    stage2fec((dataindex21(i)), 0) = stage2fec((dataindex21(i)), 0) + 1; // Yields sum of all individuals with particular transition
    
    if (choiceelement.n_elem > 0) {
      probsrates0(choiceelement(0)) = probsrates0(choiceelement(0)) + 1; // Yields sum of all individuals with particular transition
      
      if (dataalive3(i) > 0) {
        stage2fec((dataindex21(i)), 1) = stage2fec((dataindex21(i)), 1) + 1;
      }
      
      stage2fec((dataindex21(i)), 2) = stage2fec((dataindex21(i)), 2) + datausedfec2(i);
    }
    
    if (format == 2) {
      for (int j = 0; j < aes_count; j++) {
        arma::uvec choiceelementp = find(sge9index321 == dataindex321_prior(i, j)); // Added this now
        
        stage2fecp((dataindex21(i)), 0) = stage2fecp((dataindex21(i)), 0) + 1;
      
        if (choiceelementp.n_elem > 0) {
          probsrates0p(choiceelementp(0)) = probsrates0p(choiceelementp(0)) + 1;
          
          if (dataalive3(i) > 0) {
            stage2fecp((dataindex21(i)), 1) = stage2fecp((dataindex21(i)), 1) + 1;
          }
          
          stage2fecp((dataindex21(i)), 2) = stage2fecp((dataindex21(i)), 2) + datausedfec2(i);
        }
      }
    }
  }
  
  // The next bit puts together core data to be used to estimate matrix elements
  for (int i = 0; i < noelems; i++) {
    int baseindex21 = sge9index21(i);
    
    if (baseindex21 > -1) {
      arma::uvec coreelementsforchoice = find(sge3index21 == baseindex21);
      unsigned int thechosenone = coreelementsforchoice(0);
      
      probsrates1(i) = stage2fec(thechosenone, 0);
      probsrates2(i) = stage2fec(thechosenone, 1);
      probsrates3(i) = stage2fec(thechosenone, 2);
      
      if (format == 2) {
        arma::uvec coreelementsforchoicep = find(sge3index21 == baseindex21);
        unsigned int thechosenonep = coreelementsforchoicep(0);
        
        probsrates1p(i) = stage2fecp(thechosenonep, 0);
        probsrates2p(i) = stage2fecp(thechosenonep, 1);
        probsrates3p(i) = stage2fecp(thechosenonep, 2);
      }
    }
  }
  
  // Here we create the matrices
  for (int elem3 = 0; elem3 < noelems; elem3++) {
    
    if (aliveandequal(elem3) != -1) {
      
      tmatrix(aliveandequal(elem3)) = probsrates0(elem3) / probsrates1(elem3); // Survival
      
      // Fecundity
      if (format == 2) {
        fmatrix(aliveandequal(elem3)) = sge9fec32(elem3) * sge9rep2(elem3) * probsrates3p(elem3) / probsrates1p(elem3);
      } else {
        fmatrix(aliveandequal(elem3)) = sge9fec32(elem3) * sge9rep2(elem3) * probsrates3(elem3) / probsrates1(elem3);
      }
    }
  }
  
  // Now we will correct transitions and rates for given stuff
  arma::uvec ovgiventind = find(sge9ovgivent != -1);
  arma::uvec ovgivenfind = find(sge9ovgivenf != -1);
  int ovgtn = ovgiventind.n_elem;
  int ovgfn = ovgivenfind.n_elem;
  
  if (ovgtn > 0) {
    for (int i = 0; i < ovgtn; i++) {
      int matrixelement2 = aliveandequal(ovgiventind(i));
      
      tmatrix(matrixelement2) = sge9ovgivent(ovgiventind(i));
    }
  }
  
  if (ovgfn > 0) {
    for (int i = 0; i < ovgfn; i++) {
      int matrixelement2 = aliveandequal(ovgivenfind(i));
      
      fmatrix(matrixelement2) = sge9ovgivenf(ovgivenfind(i));
    }
  }
  
  // This section replaces transitions for proxy values as given in the overwrite table  
  arma::uvec ovesttind = find(sge9ovestt != -1);
  arma::uvec ovestfind = find(sge9ovestf != -1);
  int ovestn = ovesttind.n_elem;
  int ovesfn = ovestfind.n_elem;
  
  if (ovestn > 0) {
    for (int i = 0; i < ovestn; i++) {
      arma::uvec replacement = find(sge9index321 == sge9ovestt(ovesttind(i)));
      
      if (replacement.n_elem > 0) {
        tmatrix(aliveandequal(ovesttind(i))) = tmatrix(aliveandequal(replacement(0)));
      }
      
    }
  }
  
  if (ovesfn > 0) {
    for (int i = 0; i < ovesfn; i++) {
      arma::uvec replacement = find(sge9index321 == sge9ovestf(ovestfind(i)));
      
      if (replacement.n_elem > 0) {
        fmatrix(aliveandequal(ovestfind(i))) = fmatrix(aliveandequal(replacement(0)));
      }
    }
  }
  
  // The next bit changes NA to 0
  tmatrix(find_nonfinite(tmatrix)).zeros();
  fmatrix(find_nonfinite(fmatrix)).zeros();
  
  arma::mat amatrix = tmatrix + fmatrix; // Create the A matrix
  
  if (err_switch == 1) {
    arma::mat concatenated_crap = arma::join_horiz(sge9index321, aliveandequal);
    concatenated_crap = arma::join_horiz(concatenated_crap, probsrates0);
    concatenated_crap = arma::join_horiz(concatenated_crap, probsrates1);
    concatenated_crap = arma::join_horiz(concatenated_crap, probsrates2);
    concatenated_crap = arma::join_horiz(concatenated_crap, probsrates3);
    concatenated_crap = arma::join_horiz(concatenated_crap, sge9fec32);
    concatenated_crap = arma::join_horiz(concatenated_crap, probsrates0p);
    concatenated_crap = arma::join_horiz(concatenated_crap, probsrates1p);
    concatenated_crap = arma::join_horiz(concatenated_crap, probsrates2p);
    concatenated_crap = arma::join_horiz(concatenated_crap, probsrates3p);
    
    arma::mat s2f = arma::join_horiz(stage2fec, stage2fecp);

    return List::create(Named("A") = amatrix, _["U"] = tmatrix, _["F"] = fmatrix,
      _["concrp"] = concatenated_crap, _["s2f"] = s2f, _["dataprior"] = dataindex321_prior);
  } else {
    return List::create(Named("A") = amatrix, _["U"] = tmatrix, _["F"] = fmatrix);
  }
}

//' Estimate All Elements of Raw Ahistorical Population Projection Matrix
//' 
//' Function \code{normalpatrolgroup} swiftly calculates matrix transitions in
//' raw ahistorical matrices, and serves as the core workhorse function behind
//' \code{\link{rlefko2}()}.
//' 
//' @param sge3 The Allstages data frame developed for rlefko2() covering stage
//' pairs across times \emph{t}+1 and \emph{t}. Generally termed
//' \code{stageexpansion3}.
//' @param sge2 The data frame covering all stages in time \emph{t}. Generally
//' termed \code{stageexpansion2}.
//' @param MainData The demographic dataset modified to hold \code{usedfec}
//' columns.
//' @param StageFrame The full stageframe for the analysis.
//' 
//' @return List of three matrices, including the survival-transition (U)
//' matrix, the fecundity matrix (F), and the sum (A) matrix, with A first.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List normalpatrolgroup(DataFrame sge3, DataFrame sge2, DataFrame MainData,
  DataFrame StageFrame) {
  
  arma::vec sge3fec32 = sge3["b.repentry3"];
  arma::vec sge3rep2 = sge3["a.rep2n"];
  arma::vec sge3indata32 = sge3["b.indata"];
  arma::vec sge3ovgivent = sge3["b.ovgiven_t"];
  arma::vec sge3ovgivenf = sge3["b.ovgiven_f"];
  arma::vec sge3ovestt = sge3["b.ovest_t"];
  arma::vec sge3ovestf = sge3["b.ovest_f"];
  arma::vec sge3index32 = sge3["c.index321"];
  arma::vec sge3index2 = sge3["a.stage2n"];
  arma::vec aliveandequal = sge3["c.aliveandequal"];
  
  arma::vec sge2rep2 = sge2["rep2"];
  arma::vec sge2fec3 = sge2["fec3"];
  arma::vec sge2index2 = sge2["index2"];
  arma::vec sge2stage2 = sge2["stage2"];
  
  arma::vec dataindex32 = MainData["index32"];
  arma::vec dataindex2 = MainData["index2"];
  arma::vec dataalive3 = MainData["alive3"];
  arma::vec datausedfec2 = MainData["usedfec2"];
  
  arma::vec sfsizes = StageFrame["bin_size_ctr"];
  int nostages = sfsizes.n_elem;
  
  int n = dataindex32.n_elem;
  int no2stages = sge2index2.n_elem - 1; // The -1 removes the dead stage, which is still within sge2
  int noelems = sge3index32.n_elem;
  
  arma::vec probsrates0(noelems); // 1st vec = # indivs (3 trans)
  arma::vec probsrates1(noelems); // 2nd vec = total indivs for pair stage
  arma::vec probsrates2(noelems); // 3rd vec = total indivs alive for pair stage
  arma::vec probsrates3(noelems); // 4th vec = total fec for pair stage
  probsrates0.zeros();
  probsrates1.zeros();
  probsrates2.zeros();
  probsrates3.zeros();
  
  arma::mat stage2fec(no2stages, 3); //1st col = # indivs total, 2nd col = no indivs alive, 3rd col = sum fec
  stage2fec.zeros();
  
  arma::mat tmatrix((nostages-1), (nostages-1)); // Main output U matrix
  arma::mat fmatrix((nostages-1), (nostages-1)); // Main output F matrix
  tmatrix.zeros();
  fmatrix.zeros();
  
  // This main loop counts individuals going through transitions and sums their
  // fecundities, and then adds that info to the 3-trans and 2-trans tables
  for (int i = 0; i < n; i++) { 
    
    // The next line yields sum of all individuals with particular transition
    probsrates0(dataindex32(i)) = probsrates0(dataindex32(i)) + 1; 
    
    // The next line yields sum of all individuals with particular transition
    stage2fec((dataindex2(i)), 0) = stage2fec((dataindex2(i)), 0) + 1; 
    if (dataalive3(i) > 0) {
      stage2fec((dataindex2(i)), 1) = stage2fec((dataindex2(i)), 1) + 1;
    }
    
    stage2fec((dataindex2(i)), 2) = stage2fec((dataindex2(i)), 2) + datausedfec2(i);
    
  }
  
  // This next loop populates vectors of individuals according to stage in time t
  for (int i = 0; i < no2stages; i++) {
    unsigned int foradding = ((sge2stage2(i) - 1) * nostages);
    
    for (int j = 0; j < nostages; j++) {
      unsigned int entry = foradding + j;
      
      probsrates1(entry) = stage2fec(i, 0);
      probsrates2(entry) = stage2fec(i, 1);
      probsrates3(entry) = stage2fec(i, 2);
    }
  }
  
  // Here we populate the main U and F matrices
  for (int elem3 = 0; elem3 < noelems; elem3++) {
    
    if (aliveandequal(elem3) != -1) {
      
      // The next lines DO leave NaNs in the matrices when 0 individuals are summed through in probsrates1
      tmatrix(aliveandequal(elem3)) = probsrates0(elem3) / probsrates1(elem3);
      fmatrix(aliveandequal(elem3)) = sge3fec32(elem3) * sge3rep2(elem3) * probsrates3(elem3) / probsrates1(elem3);
    }
  }
  
  // This section corrects for transitions given in the overwrite table
  arma::uvec ovgiventind = find(sge3ovgivent != -1);
  arma::uvec ovgivenfind = find(sge3ovgivenf != -1);
  int ovgtn = ovgiventind.n_elem;
  int ovgfn = ovgivenfind.n_elem;
  
  if (ovgtn > 0) {
    for (int i = 0; i < ovgtn; i++) {
      int matrixelement2 = aliveandequal(ovgiventind(i));
      
      tmatrix(matrixelement2) = sge3ovgivent(ovgiventind(i));
    }
  }
  
  if (ovgfn > 0) {
    for (int i = 0; i < ovgfn; i++) {
      int matrixelement2 = aliveandequal(ovgivenfind(i));
      
      fmatrix(matrixelement2) = sge3ovgivenf(ovgivenfind(i));
    }
  }
  
  // This section replaces transitions with proxies as given in the overwrite table
  arma::uvec ovesttind = find(sge3ovestt != -1);
  arma::uvec ovestfind = find(sge3ovestf != -1);
  int ovestn = ovesttind.n_elem;
  int ovesfn = ovestfind.n_elem;
  
  if (ovestn > 0) {
    for (int i = 0; i < ovestn; i++) {
      arma::uvec replacement = find(sge3index32 == sge3ovestt(ovesttind(i)));
      
      tmatrix(aliveandequal(ovesttind(i))) = tmatrix(aliveandequal(replacement(0)));
    }
  }
  
  if (ovesfn > 0) {
    for (int i = 0; i < ovesfn; i++) {
      arma::uvec replacement = find(sge3index32 == sge3ovestf(ovestfind(i)));
      
      fmatrix(aliveandequal(ovestfind(i))) = fmatrix(aliveandequal(replacement(0)));
    }
  }
  
  // The next bit changes NAs to 0
  tmatrix(find_nonfinite(tmatrix)).zeros();
  fmatrix(find_nonfinite(fmatrix)).zeros();
  
  arma::mat amatrix = tmatrix + fmatrix; // Create the A matrix
  
  return List::create(Named("A") = amatrix, _["U"] = tmatrix, _["F"] = fmatrix);
}

//' Estimate All Elements of Function-based Population Projection Matrix
//' 
//' Function \code{jerzeibalowski} swiftly calculates matrix elements in
//' function-based population projection matrices. Used in
//' \code{\link{flefko3}()}, \code{\link{flefko2}()}, and
//' \code{\link{aflefko2}()}.
//' 
//' @param ppy A data frame with one row, showing the population, patch, and
//' year.
//' @param AllStages A large data frame giving all required inputs for vital
//' rate estimation other than the vital rate model coefficients themselves.
//' Contains a row for each ultimate matrix element.
//' @param stageframe The modified stageframe used in matrix calculations.
//' @param matrixformat An integer representing the style of matrix to develop.
//' Options include Ehrlen-format hMPM (1), deVries-format hMPM (2), ahMPM (3),
//' and age-by-stage MPM (4).
//' @param survproxy List of coefficients estimated in model of survival.
//' @param obsproxy List of coefficients estimated in model of observation.
//' @param sizeproxy List of coefficients estimated in model of size.
//' @param repstproxy List of coefficients estimated in model of reproductive 
//' status.
//' @param fecproxy List of coefficients estimated in model of fecundity.
//' @param jsurvproxy List of coefficients estimated in model of juvenile
//' survival.
//' @param jobsproxy List of coefficients estimated in model of juvenile
//' observation.
//' @param jsizeproxy List of coefficients estimated in model of juvenile size.
//' @param jrepstproxy List of coefficients estimated in model of juvenile
//' reproductive status.
//' @param inda A numeric value equal to the value of individual covariate a to
//' be used in analysis.
//' @param indb A numeric value equal to the value of individual covariate b to
//' be used in analysis.
//' @param indc A numeric value equal to the value of individual covariate c to
//' be used in analysis.
//' @param survdev Scalar value to be added to the y-intercept of the linear
//' model of survival probability.
//' @param obsdev Scalar value to be added to the y-intercept of the linear
//' model of observation probability.
//' @param sizedev Scalar value to be added to the y-intercept of the linear
//' model of size transition.
//' @param repstdev Scalar value to be added to the y-intercept of the linear
//' model of reproduction probability.
//' @param fecdev Scalar value to be added to the y-intercept of the linear
//' model of fecundity.
//' @param jsurvdev Scalar value to be added to the y-intercept of the linear
//' model of juvenile survival probability.
//' @param jobsdev Scalar value to be added to the y-intercept of the linear
//' model of juvenile observation probability.
//' @param jsizedev Scalar value to be added to the y-intercept of the linear
//' model of juvenile size transition.
//' @param jrepstdev Scalar value to be added to the y-intercept of the linear
//' model of juvenile reproduction probability.
//' @param fecmod A scalar multiplier for fecundity.
//' @param summedvars Summed variance-covariance terms in Poisson size
//' distribution.
//' @param sigma Standard deviation of Gaussian size distribution.
//' @param jsummedvars Summed variance-covariance terms in Poisson juvenile size
//' distribution.
//' @param jsigma Standard deviation of Gaussian juvenile size distribution.
//' @param maxsize The maximum size to be used in element estimation.
//' @param finalage The final age to be included in age-by-stage MPM estimation.
//' @param sizedist Designates whether size is Gaussian (2), Poisson (0), or
//' negative binomial (1) distributed.
//' @param fecdist Designates whether fecundity is Gaussian (2), Poisson (0), or
//' negative binomial (1) distributed.
//' @param negfec Logical value denoting whether to change negative estimated
//` fecundity to 0.
//' 
//' @return A list of 3 matrices, including the main MPM (A), the survival-
//' transition matrix (U), anf a fecundity matrix (F). With tweaking, can also
//' produce a 4 column matrix showing survival probability, observation
//' probability, reproduction probability, and size transition probability, for
//' each element of the final MPM.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List jerzeibalowski(DataFrame ppy, DataFrame AllStages, DataFrame stageframe,
  int matrixformat, List survproxy, List obsproxy, List sizeproxy,
  List repstproxy, List fecproxy, List jsurvproxy, List jobsproxy,
  List jsizeproxy, List jrepstproxy, double inda, double indb, double indc,
  double survdev, double obsdev, double sizedev, double repstdev, double fecdev,
  double jsurvdev, double jobsdev, double jsizedev, double jrepstdev,
  double fecmod, double summedvars, double sigma, double jsummedvars,
  double jsigma, double maxsize, unsigned int finalage, int sizedist, int fecdist,
  bool negfec) {
  
  // The DataFrame AllStages introduces variables used in size and fecundity calculations. This DataFrame
  // is broken up into long vectors composed of input sizes and related variables for these calculations. 
  // The "model" Lists bring in the vital rate models, and include random coefficients
  // where needed. We also have a number of extra variables, that include such info as whether to use
  // the Poisson, negative binomial, and Gaussian for size and fecundity calculations. If either sizedist
  // or fecdist equals 0, then the Poisson is used. If either equals 1, then the negative binomial is 
  // used. If 2, then the Gaussian.
  
  // In the code below, the function decides on an appropriate single loop routine based on which coefficient
  // vectors have length greater than 1. Only a single routine will run because of the series of else statements
  
  // Determines the size of the matrix
  StringVector stagenames = stageframe["stage"];
  int nostages = stagenames.length();
  unsigned long matrixdim {0};
  
  int nostages_counter = nostages;
  for (int i=0; i < nostages_counter; i++) {
    if (stagenames(i) == "AlmostBorn") nostages -= 1;  
    if (stagenames(i) == "Dead") nostages -= 1;
  }
  
  if (matrixformat == 1) { // Ehrlen-format hMPM
    matrixdim = nostages * nostages;
  } else if (matrixformat == 2) { // deVries-format hMPM
    matrixdim = nostages * (nostages + 1);
  } else if (matrixformat == 3) { // ahMPM
    matrixdim = nostages;
  } else if (matrixformat == 4) { // age-by-stage MPM
    matrixdim = nostages * (finalage + 1);
  }
  
  bool sizezero = FALSE;
  bool feczero = FALSE;
  bool jsizezero = FALSE;
  double sizezerosum {0};
  double feczerosum {0};
  double jsizezerosum {0};
  
  arma::vec survcoefs = survproxy["coefficients"];
  arma::vec obscoefs = obsproxy["coefficients"];
  arma::vec sizecoefs = sizeproxy["coefficients"];
  arma::vec repstcoefs = repstproxy["coefficients"];
  arma::vec feccoefs = fecproxy["coefficients"];
  arma::vec jsurvcoefs = jsurvproxy["coefficients"];
  arma::vec jobscoefs = jobsproxy["coefficients"];
  arma::vec jsizecoefs = jsizeproxy["coefficients"];
  arma::vec jrepstcoefs = jrepstproxy["coefficients"];
  
  int sizetrunc = sizeproxy["trunc"];
  int jsizetrunc = jsizeproxy["trunc"];
  
  int survl = survcoefs.n_elem;
  int obsl = obscoefs.n_elem;
  int sizel = sizecoefs.n_elem;
  int repstl = repstcoefs.n_elem;
  int fecl = feccoefs.n_elem;
  int jsurvl = jsurvcoefs.n_elem;
  int jobsl = jobscoefs.n_elem;
  int jsizel = jsizecoefs.n_elem;
  int jrepstl = jrepstcoefs.n_elem;
  
  for (int i = 46; i < 92; i++) {
    if (sizel > 1) {
      sizezerosum += sizecoefs(i);
    } else sizezerosum = 0;
    
    if (fecl > 1) {
      feczerosum += feccoefs(i);
    } else feczerosum = 0;
    
    if (jsizel > 1) {
      jsizezerosum += jsizecoefs(i);
    } else jsizezerosum = 0;
  }
  
  arma::vec survyear = survproxy["years"];
  arma::vec obsyear = obsproxy["years"];
  arma::vec sizeyear = sizeproxy["years"];
  arma::vec repstyear = repstproxy["years"];
  arma::vec fecyear = fecproxy["years"];
  arma::vec jsurvyear = jsurvproxy["years"];
  arma::vec jobsyear = jobsproxy["years"];
  arma::vec jsizeyear = jsizeproxy["years"];
  arma::vec jrepstyear = jrepstproxy["years"];
  
  arma::vec survpatch = survproxy["patches"];
  arma::vec obspatch = obsproxy["patches"];
  arma::vec sizepatch = sizeproxy["patches"];
  arma::vec repstpatch = repstproxy["patches"];
  arma::vec fecpatch = fecproxy["patches"];
  arma::vec jsurvpatch = jsurvproxy["patches"];
  arma::vec jobspatch = jobsproxy["patches"];
  arma::vec jsizepatch = jsizeproxy["patches"];
  arma::vec jrepstpatch = jrepstproxy["patches"];
  
  if (NumericVector::is_na(survpatch(0))) {survpatch(0) = 0;}
  if (NumericVector::is_na(obspatch(0))) {obspatch(0) = 0;}
  if (NumericVector::is_na(sizepatch(0))) {sizepatch(0) = 0;}
  if (NumericVector::is_na(repstpatch(0))) {repstpatch(0) = 0;}
  if (NumericVector::is_na(fecpatch(0))) {fecpatch(0) = 0;}
  if (NumericVector::is_na(jsurvpatch(0))) {jsurvpatch(0) = 0;}
  if (NumericVector::is_na(jobspatch(0))) {jobspatch(0) = 0;}
  if (NumericVector::is_na(jsizepatch(0))) {jsizepatch(0) = 0;}
  if (NumericVector::is_na(jrepstpatch(0))) {jrepstpatch(0) = 0;}
  
  arma::vec sizeyearzi = sizeproxy["zeroyear"];
  arma::vec fecyearzi = fecproxy["zeroyear"];
  arma::vec jsizeyearzi = jsizeproxy["zeroyear"];
  
  if (!NumericVector::is_na(sizeyearzi(0)) && sizezerosum != 0) sizezero = TRUE;
  if (!NumericVector::is_na(fecyearzi(0)) && feczerosum != 0) feczero = TRUE;
  if (!NumericVector::is_na(jsizeyearzi(0)) && jsizezerosum != 0) jsizezero = TRUE;
  
  if (NumericVector::is_na(sizeyearzi(0))) {sizeyearzi(0) = 0;}
  if (NumericVector::is_na(fecyearzi(0))) {fecyearzi(0) = 0;}
  if (NumericVector::is_na(jsizeyearzi(0))) {jsizeyearzi(0) = 0;}
  
  arma::vec sizepatchzi = sizeproxy["zeropatch"];
  arma::vec fecpatchzi = fecproxy["zeropatch"];
  arma::vec jsizepatchzi = jsizeproxy["zeropatch"];
  
  if (NumericVector::is_na(sizepatchzi(0))) {sizepatchzi(0) = 0;}
  if (NumericVector::is_na(fecpatchzi(0))) {fecpatchzi(0) = 0;}
  if (NumericVector::is_na(jsizepatchzi(0))) {jsizepatchzi(0) = 0;}
  
  double sizesigma = sizeproxy["sigma"];
  double fecsigma = fecproxy["sigma"];
  double jsizesigma = jsizeproxy["sigma"];
  
  if (NumericVector::is_na(sizesigma)) {
    if (sizedist == 1) {
      sizesigma = 1;
    } else {
      sizesigma = 0;
    }
  }
  if (NumericVector::is_na(fecsigma)) {
    if (fecdist == 1) {
      fecsigma = 1;
    } else {
      fecsigma = 0;
    }
  }
  if (NumericVector::is_na(jsizesigma)) {
    if (sizedist == 1) {
      jsizesigma = 1;
    } else {
      jsizesigma = 0;
    }
  }
  
  arma::vec stage3 = AllStages["a.stage3"];
  arma::vec stage2n = AllStages["a.stage2n"];
  arma::vec stage2o = AllStages["a.stage2o"];
  Rcpp::NumericVector sz3 = AllStages["a.size3"];
  Rcpp::NumericVector sz2n = AllStages["a.size2n"];
  Rcpp::NumericVector sz2o = AllStages["a.size2o"];
  Rcpp::NumericVector sz1 = AllStages["a.size1"];
  Rcpp::NumericVector ob3 = AllStages["a.obs3"];
  Rcpp::NumericVector fl3 = AllStages["a.rep3"];
  Rcpp::NumericVector fl2n = AllStages["a.rep2n"];
  Rcpp::NumericVector fl2o = AllStages["a.rep2o"];
  Rcpp::NumericVector fl1 = AllStages["a.rep1"];
  Rcpp::NumericVector mat3 = AllStages["a.mat3"];
  Rcpp::NumericVector mat2n = AllStages["a.mat2n"];
  Rcpp::NumericVector mat2o = AllStages["a.mat2o"];
  Rcpp::NumericVector mat1 = AllStages["a.mat1"];
  Rcpp::NumericVector immat2n = AllStages["b.imm2n"];
  Rcpp::NumericVector immat2o = AllStages["b.imm2o"];
  Rcpp::NumericVector immat1 = AllStages["b.imm1"];
  Rcpp::NumericVector indata2n = AllStages["b.indata2n"];
  Rcpp::NumericVector indata2o = AllStages["b.indata2o"];
  Rcpp::NumericVector repentry = AllStages["b.repentry3"];
  Rcpp::NumericVector binwidth3 = AllStages["b.binwidth"];
  Rcpp::NumericVector actualage2 = AllStages["b.actualage"];
  Rcpp::NumericVector indata = AllStages["b.indata"];
  arma::vec ovestt = AllStages["b.ovest_t"];
  Rcpp::NumericVector ovgivent = AllStages["b.ovgiven_t"];
  arma::vec ovestf = AllStages["b.ovest_f"];
  Rcpp::NumericVector ovgivenf = AllStages["b.ovgiven_f"];
  arma::uvec index321 = AllStages["c.index321"];
  Rcpp::NumericVector aliveandequal = AllStages["c.aliveandequal"];
  
  int n = stage3.n_elem;
  
  arma::uvec years = ppy["yearorder"];
  int yearnumber = years(0) - 1;
  
  arma::uvec patches = ppy["patchorder"];
  int patchnumber = patches(0) - 1;
  
  arma::uvec replacetvec = find(ovestt != -1);
  arma::uvec replacefvec = find(ovestf != -1);
  int replacementst = replacetvec.n_elem;
  int replacementsf = replacefvec.n_elem;
  int repindex {0};
  int properindex {0};
  int proxyindex {0};
  
  arma::mat out(n, 4);  // 0 matrix with n rows & 4 columns: 0 surv, 1 obs, 2 repst, 3 size, >3 are test variables
  arma::mat survtransmat(matrixdim, matrixdim);
  arma::mat fectransmat(matrixdim, matrixdim);
  
  out.zeros();
  survtransmat.zeros();
  fectransmat.zeros();
  
  double mu {0};
  double lambda {0};
  double lambda_preexp {0};
  double mu_preexp {0};
  
  // The following loop runs through each line of AllStages, and so runs through each estimable element in the matrix
  for(int i = 0; i < n; i++) {
    unsigned int k = aliveandequal(i);
    
    if (ovgivent(i) == -1 && indata(i) == 1 && stage2n(i) == stage2o(i)) {
      if ((mat2n(i) == 1 && mat3(i) == 1) || (mat2o(i) == 1 && mat3(i) == 1)) {
        // Adult survival transitions
        
        arma::vec preout(4);
        
        if (survl > 1) {
          
          preout(0) = (survcoefs(0) + (survcoefs(1) * fl1(i)) + (survcoefs(2) * fl2n(i)) +
            (survcoefs(3) * sz1(i)) + (survcoefs(4) * sz2o(i)) + (survcoefs(5) * fl2n(i) * fl1(i)) + 
            (survcoefs(6) * sz2o(i) * sz1(i)) + (survcoefs(7) * sz1(i) * fl1(i)) +
            (survcoefs(8) * sz2o(i) * fl2n(i)) + (survcoefs(9) * sz2o(i) * fl1(i)) + 
            (survcoefs(10) * sz1(i) * fl2n(i)) + (survcoefs(11) * actualage2(i)) + 
            (survcoefs(12) * actualage2(i) * sz1(i)) + (survcoefs(13) * actualage2(i) * sz2o(i)) + 
            (survcoefs(14) * actualage2(i) * fl1(i)) + (survcoefs(15) * actualage2(i) * fl2n(i)) + 
            (survcoefs(16) * inda) + (survcoefs(17) * indb) + (survcoefs(18) * indc) + 
            (survcoefs(19) * inda) + (survcoefs(20) * indb) + (survcoefs(21) * indc) + 
            (survcoefs(22) * inda * sz2o(i)) + (survcoefs(23) * indb * sz2o(i)) + 
            (survcoefs(24) * indc * sz2o(i)) + (survcoefs(25) * inda * fl2n(i)) + 
            (survcoefs(26) * indb * fl2n(i)) + (survcoefs(27) * indc * fl2n(i)) + 
            (survcoefs(28) * inda * sz1(i)) + (survcoefs(29) * indb * sz1(i)) +
            (survcoefs(30) * indc * sz1(i)) + (survcoefs(31) * inda * fl1(i)) + 
            (survcoefs(32) * indb * fl1(i)) + (survcoefs(33) * indc * fl1(i)) + 
            (survcoefs(34) * inda * indb) + (survcoefs(35) * inda * indc) + 
            (survcoefs(36) * indb * indc) + (survcoefs(37) * inda * indb) + 
            (survcoefs(38) * inda * indc) + (survcoefs(39) * indb * indc) + 
            (survcoefs(40) * inda * indb) + (survcoefs(41) * inda * indb) +
            (survcoefs(42) * inda * indc) + (survcoefs(43) * inda * indc) + 
            (survcoefs(44) * indb * indc) + (survcoefs(45) * indb * indc) + 
            survpatch(patchnumber) + survyear(yearnumber) + survdev);
          
          if (preout(0) > 700) preout(0) = 700; // This catches numbers too high to be dealt with properly
          out(i, 0) = exp(preout(0)) / (1 + exp(preout(0)));
          
        } else {
          out(i, 0) = survcoefs(0);
        }
        
        if (obsl > 1) {
          
          preout(1) = (obscoefs(0) + (obscoefs(1) * fl1(i)) + (obscoefs(2) * fl2n(i)) + (obscoefs(3) * sz1(i)) + 
            (obscoefs(4) * sz2o(i)) + (obscoefs(5) * fl2n(i) * fl1(i)) + (obscoefs(6) * sz2o(i) * sz1(i)) +
            (obscoefs(7) * sz1(i) * fl1(i)) + (obscoefs(8) * sz2o(i) * fl2n(i)) + (obscoefs(9) * sz2o(i) * fl1(i)) + 
            (obscoefs(10) * sz1(i) * fl2n(i)) + (obscoefs(11) * actualage2(i)) + (obscoefs(12) * actualage2(i) * sz1(i)) + 
            (obscoefs(13) * actualage2(i) * sz2o(i)) + (obscoefs(14) * actualage2(i) * fl1(i)) + 
            (obscoefs(15) * actualage2(i) * fl2n(i)) + (obscoefs(16) * inda) + (obscoefs(17) * indb) + 
            (obscoefs(18) * indc) + (obscoefs(19) * inda) + (obscoefs(20) * indb) + (obscoefs(21) * indc) + 
            (obscoefs(22) * inda * sz2o(i)) + (obscoefs(23) * indb * sz2o(i)) + (obscoefs(24) * indc * sz2o(i)) + 
            (obscoefs(25) * inda * fl2n(i)) + (obscoefs(26) * indb * fl2n(i)) + (obscoefs(27) * indc * fl2n(i)) + 
            (obscoefs(28) * inda * sz1(i)) + (obscoefs(29) * indb * sz1(i)) + (obscoefs(30) * indc * sz1(i)) + 
            (obscoefs(31) * inda * fl1(i)) + (obscoefs(32) * indb * fl1(i)) + (obscoefs(33) * indc * fl1(i)) + 
            (obscoefs(34) * inda * indb) + (obscoefs(35) * inda * indc) + (obscoefs(36) * indb * indc) + 
            (obscoefs(37) * inda * indb) + (obscoefs(38) * inda * indc) + (obscoefs(39) * indb * indc) + 
            (obscoefs(40) * inda * indb) + (obscoefs(41) * inda * indb) + (obscoefs(42) * inda * indc) + 
            (obscoefs(43) * inda * indc) + (obscoefs(44) * indb * indc) + (obscoefs(45) * indb * indc) + 
            obspatch(patchnumber) + obsyear(yearnumber) + obsdev);
          
          if (preout(1) > 700) preout(1) = 700; // This catches numbers too high to be dealt with properly
          out(i, 1) = exp(preout(1)) / (1 + exp(preout(1)));
          
        } else {
          out(i, 1) = obscoefs(0);
        }
        
        if (ob3(i) == 1 || obsl == 1) {
          
          if (sizel > 1) {
            if (sizedist == 0) {
              // Poisson size distribution
              
              if (sizezero == TRUE && sz3(i) == 0) {
                lambda_preexp = (sizecoefs(46) + (sizecoefs(47) * fl1(i)) + 
                  (sizecoefs(48) * fl2n(i)) + (sizecoefs(49) * sz1(i)) +
                  (sizecoefs(50) * sz2o(i)) + (sizecoefs(51) * fl2n(i) * fl1(i)) +
                  (sizecoefs(52) * sz2o(i) * sz1(i)) + (sizecoefs(53) * sz1(i) * fl1(i)) + 
                  (sizecoefs(54) * sz2o(i) * fl2n(i)) + (sizecoefs(55) * sz2o(i) * fl1(i)) + 
                  (sizecoefs(56) * sz1(i) * fl2n(i)) + (sizecoefs(57) * actualage2(i)) + 
                  (sizecoefs(58) * actualage2(i) * sz1(i)) + (sizecoefs(59) * actualage2(i) * sz2o(i)) + 
                  (sizecoefs(60) * actualage2(i) * fl1(i)) + (sizecoefs(61) * actualage2(i) * fl2n(i)) + 
                  (sizecoefs(62) * inda) + (sizecoefs(63) * indb) + (sizecoefs(64) * indc) + 
                  (sizecoefs(65) * inda) + (sizecoefs(66) * indb) + (sizecoefs(67) * indc) + 
                  (sizecoefs(68) * inda * sz2o(i)) + (sizecoefs(69) * indb * sz2o(i)) + 
                  (sizecoefs(70) * indc * sz2o(i)) + (sizecoefs(71) * inda * fl2n(i)) + 
                  (sizecoefs(72) * indb * fl2n(i)) + (sizecoefs(73) * indc * fl2n(i)) + 
                  (sizecoefs(74) * inda * sz1(i)) + (sizecoefs(75) * indb * sz1(i)) + 
                  (sizecoefs(76) * indc * sz1(i)) + (sizecoefs(77) * inda * fl1(i)) + 
                  (sizecoefs(78) * indb * fl1(i)) + (sizecoefs(79) * indc * fl1(i)) + 
                  (sizecoefs(80) * inda * indb) + (sizecoefs(81) * inda * indc) +
                  (sizecoefs(82) * indb * indc) + (sizecoefs(83) * inda * indb) + 
                  (sizecoefs(84) * inda * indc) + (sizecoefs(85) * indb * indc) + 
                  (sizecoefs(86) * inda * indb) + (sizecoefs(87) * inda * indb) + 
                  (sizecoefs(88) * inda * indc) + (sizecoefs(89) * inda * indc) +
                  (sizecoefs(90) * indb * indc) + (sizecoefs(91) * indb * indc) + 
                  sizepatchzi(patchnumber) + sizeyearzi(yearnumber) + sizedev + (summedvars / 2));
                
                if (lambda_preexp > 700) {
                  lambda = exp(static_cast<double>(700));
                } else lambda = exp(lambda_preexp);
                out(i, 3) = (lambda) / (1 + (lambda));
                
              } else {
                double sizefac = sz3(i) * tgamma(sz3(i));
                lambda_preexp = (sizecoefs(0) + (sizecoefs(1) * fl1(i)) + 
                  (sizecoefs(2) * fl2n(i)) + (sizecoefs(3) * sz1(i)) +
                  (sizecoefs(4) * sz2o(i)) + (sizecoefs(5) * fl2n(i) * fl1(i)) + 
                  (sizecoefs(6) * sz2o(i) * sz1(i)) + (sizecoefs(7) * sz1(i) * fl1(i)) + 
                  (sizecoefs(8) * sz2o(i) * fl2n(i)) + (sizecoefs(9) * sz2o(i) * fl1(i)) + 
                  (sizecoefs(10) * sz1(i) * fl2n(i)) + (sizecoefs(11) * actualage2(i)) + 
                  (sizecoefs(12) * actualage2(i) * sz1(i)) + (sizecoefs(13) * actualage2(i) * sz2o(i)) + 
                  (sizecoefs(14) * actualage2(i) * fl1(i)) + (sizecoefs(15) * actualage2(i) * fl2n(i)) + 
                  (sizecoefs(16) * inda) + (sizecoefs(17) * indb) + (sizecoefs(18) * indc) + 
                  (sizecoefs(19) * inda) + (sizecoefs(20) * indb) + (sizecoefs(21) * indc) + 
                  (sizecoefs(22) * inda * sz2o(i)) + (sizecoefs(23) * indb * sz2o(i)) + 
                  (sizecoefs(24) * indc * sz2o(i)) + (sizecoefs(25) * inda * fl2n(i)) + 
                  (sizecoefs(26) * indb * fl2n(i)) + (sizecoefs(27) * indc * fl2n(i)) + 
                  (sizecoefs(28) * inda * sz1(i)) + (sizecoefs(29) * indb * sz1(i)) + 
                  (sizecoefs(30) * indc * sz1(i)) + (sizecoefs(31) * inda * fl1(i)) + 
                  (sizecoefs(32) * indb * fl1(i)) + (sizecoefs(33) * indc * fl1(i)) + 
                  (sizecoefs(34) * inda * indb) + (sizecoefs(35) * inda * indc) + 
                  (sizecoefs(36) * indb * indc) + (sizecoefs(37) * inda * indb) + 
                  (sizecoefs(38) * inda * indc) + (sizecoefs(39) * indb * indc) + 
                  (sizecoefs(40) * inda * indb) + (sizecoefs(41) * inda * indb) + 
                  (sizecoefs(42) * inda * indc) + (sizecoefs(43) * inda * indc) + 
                  (sizecoefs(44) * indb * indc) + (sizecoefs(45) * indb * indc) + 
                  sizepatch(patchnumber) + sizeyear(yearnumber) + sizedev + 
                  (summedvars / 2));
                
                if (lambda_preexp > 500) {
                  lambda = exp(static_cast<double>(500));
                } else lambda = exp(lambda_preexp);
                
                if (sizetrunc == 1) {
                  out(i, 3) = ((pow(lambda, sz3(i)) * exp(-1 * lambda)) / sizefac) / (1 - (exp(-1 * lambda)));
                } else {
                  out(i, 3) = ((pow(lambda, sz3(i)) * exp(-1 * lambda)) / sizefac);
                }
              }
              
            } else if (sizedist == 1) {
              // Negative binomial size distribution
              
              if (sizezero == TRUE && sz3(i) == 0) {
                mu_preexp = (sizecoefs(46) + (sizecoefs(47) * fl1(i)) + 
                  (sizecoefs(48) * fl2n(i)) + (sizecoefs(49) * sz1(i)) +
                  (sizecoefs(50) * sz2o(i)) + (sizecoefs(51) * fl2n(i) * fl1(i)) + 
                  (sizecoefs(52) * sz2o(i) * sz1(i)) + (sizecoefs(53) * sz1(i) * fl1(i)) + 
                  (sizecoefs(54) * sz2o(i) * fl2n(i)) + (sizecoefs(55) * sz2o(i) * fl1(i)) + 
                  (sizecoefs(56) * sz1(i) * fl2n(i)) + (sizecoefs(57) * actualage2(i)) + 
                  (sizecoefs(58) * actualage2(i) * sz1(i)) + (sizecoefs(59) * actualage2(i) * sz2o(i)) + 
                  (sizecoefs(60) * actualage2(i) * fl1(i)) + (sizecoefs(61) * actualage2(i) * fl2n(i)) + 
                  (sizecoefs(62) * inda) + (sizecoefs(63) * indb) + (sizecoefs(64) * indc) + 
                  (sizecoefs(65) * inda) + (sizecoefs(66) * indb) + (sizecoefs(67) * indc) + 
                  (sizecoefs(68) * inda * sz2o(i)) + (sizecoefs(69) * indb * sz2o(i)) + 
                  (sizecoefs(70) * indc * sz2o(i)) + (sizecoefs(71) * inda * fl2n(i)) + 
                  (sizecoefs(72) * indb * fl2n(i)) + (sizecoefs(73) * indc * fl2n(i)) + 
                  (sizecoefs(74) * inda * sz1(i)) + (sizecoefs(75) * indb * sz1(i)) + 
                  (sizecoefs(76) * indc * sz1(i)) + (sizecoefs(77) * inda * fl1(i)) + 
                  (sizecoefs(78) * indb * fl1(i)) + (sizecoefs(79) * indc * fl1(i)) + 
                  (sizecoefs(80) * inda * indb) + (sizecoefs(81) * inda * indc) +
                  (sizecoefs(82) * indb * indc) + (sizecoefs(83) * inda * indb) + 
                  (sizecoefs(84) * inda * indc) + (sizecoefs(85) * indb * indc) + 
                  (sizecoefs(86) * inda * indb) + (sizecoefs(87) * inda * indb) +
                  (sizecoefs(88) * inda * indc) + (sizecoefs(89) * inda * indc) + 
                  (sizecoefs(90) * indb * indc) + (sizecoefs(91) * indb * indc) + 
                  sizepatchzi(patchnumber) + sizeyearzi(yearnumber) + sizedev + 
                  (summedvars / 2));
                
                if (mu_preexp > 700) {
                  mu = exp(static_cast<double>(700));
                } else mu = exp(mu_preexp);
                
                out(i, 3) = (mu) / (1 + (mu));
              } else {
                mu_preexp = (sizecoefs(0) + (sizecoefs(1) * fl1(i)) +
                  (sizecoefs(2) * fl2n(i)) + (sizecoefs(3) * sz1(i)) +
                  (sizecoefs(4) * sz2o(i)) + (sizecoefs(5) * fl2n(i) * fl1(i)) +
                  (sizecoefs(6) * sz2o(i) * sz1(i)) + (sizecoefs(7) * sz1(i) * fl1(i)) + 
                  (sizecoefs(8) * sz2o(i) * fl2n(i)) + (sizecoefs(9) * sz2o(i) * fl1(i)) + 
                  (sizecoefs(10) * sz1(i) * fl2n(i)) + (sizecoefs(11) * actualage2(i)) + 
                  (sizecoefs(12) * actualage2(i) * sz1(i)) + (sizecoefs(13) * actualage2(i) * sz2o(i)) + 
                  (sizecoefs(14) * actualage2(i) * fl1(i)) + (sizecoefs(15) * actualage2(i) * fl2n(i)) + 
                  (sizecoefs(16) * inda) + (sizecoefs(17) * indb) + (sizecoefs(18) * indc) + 
                  (sizecoefs(19) * inda) + (sizecoefs(20) * indb) + (sizecoefs(21) * indc) + 
                  (sizecoefs(22) * inda * sz2o(i)) + (sizecoefs(23) * indb * sz2o(i)) + 
                  (sizecoefs(24) * indc * sz2o(i)) + (sizecoefs(25) * inda * fl2n(i)) + 
                  (sizecoefs(26) * indb * fl2n(i)) + (sizecoefs(27) * indc * fl2n(i)) + 
                  (sizecoefs(28) * inda * sz1(i)) + (sizecoefs(29) * indb * sz1(i)) + 
                  (sizecoefs(30) * indc * sz1(i)) + (sizecoefs(31) * inda * fl1(i)) + 
                  (sizecoefs(32) * indb * fl1(i)) + (sizecoefs(33) * indc * fl1(i)) + 
                  (sizecoefs(34) * inda * indb) + (sizecoefs(35) * inda * indc) + 
                  (sizecoefs(36) * indb * indc) + (sizecoefs(37) * inda * indb) + 
                  (sizecoefs(38) * inda * indc) + (sizecoefs(39) * indb * indc) + 
                  (sizecoefs(40) * inda * indb) + (sizecoefs(41) * inda * indb) +
                  (sizecoefs(42) * inda * indc) + (sizecoefs(43) * inda * indc) + 
                  (sizecoefs(44) * indb * indc) + (sizecoefs(45) * indb * indc) + 
                  sizepatch(patchnumber) + sizeyear(yearnumber) + sizedev);
                
                if (mu_preexp > 500) {
                  mu = exp(static_cast<double>(500));
                } else mu = exp(mu_preexp);
                
                double theta = sizesigma;
                if (theta > 1000000000) {
                  theta = 1000000000;
                }
                double alpha = 1 / theta;
                
                int y = static_cast<int>(sz3(i));
                
                double log_leftie = 0;
                for (int j = 0; j < y; j++) {
                  log_leftie = log(static_cast<double>(j) + theta) - log(static_cast<double>(j) + 1) + log_leftie;
                }
                
                double log_amu = log(alpha) + log(mu);
                double log_mid = -1 * theta * log(1 + (alpha * mu));
                
                double log_rightie = sz3(i) * (log_amu - log(1 + (alpha * mu)));
                
                double raw_prob = log_leftie + log_mid + log_rightie;
                
                if (sizetrunc == 1) {
                  double zero_raw_prob = log_mid;
                  
                  out(i, 3) = exp(raw_prob) / (1 - exp(zero_raw_prob));
                } else {
                  out(i, 3) = exp(raw_prob);
                }
                
                // out(i, 4) = log_leftie;
                // out(i, 5) = log_mid;
                // out(i, 6) = log_rightie;
              }
              
            } else if (sizedist == 2) {
              // Gaussian size distribution
              
              double sigma2 = sigma * sigma;
              preout(3) = (sizecoefs(0) + (sizecoefs(1) * fl1(i)) + (sizecoefs(2) * fl2n(i)) + 
                (sizecoefs(3) * sz1(i)) + (sizecoefs(4) * sz2o(i)) + (sizecoefs(5) * fl2n(i) * fl1(i)) + 
                (sizecoefs(6) * sz2o(i) * sz1(i)) + (sizecoefs(7) * sz1(i) * fl1(i)) + 
                (sizecoefs(8) * sz2o(i) * fl2n(i)) + (sizecoefs(9) * sz2o(i) * fl1(i)) + 
                (sizecoefs(10) * sz1(i) * fl2n(i)) + (sizecoefs(11) * actualage2(i)) + 
                (sizecoefs(12) * actualage2(i) * sz1(i)) + (sizecoefs(13) * actualage2(i) * sz2o(i)) + 
                (sizecoefs(14) * actualage2(i) * fl1(i)) + (sizecoefs(15) * actualage2(i) * fl2n(i)) + 
                (sizecoefs(16) * inda) + (sizecoefs(17) * indb) + (sizecoefs(18) * indc) + 
                (sizecoefs(19) * inda) + (sizecoefs(20) * indb) + (sizecoefs(21) * indc) + 
                (sizecoefs(22) * inda * sz2o(i)) + (sizecoefs(23) * indb * sz2o(i)) + 
                (sizecoefs(24) * indc * sz2o(i)) + (sizecoefs(25) * inda * fl2n(i)) + 
                (sizecoefs(26) * indb * fl2n(i)) + (sizecoefs(27) * indc * fl2n(i)) + 
                (sizecoefs(28) * inda * sz1(i)) + (sizecoefs(29) * indb * sz1(i)) + 
                (sizecoefs(30) * indc * sz1(i)) + (sizecoefs(31) * inda * fl1(i)) + 
                (sizecoefs(32) * indb * fl1(i)) + (sizecoefs(33) * indc * fl1(i)) + 
                (sizecoefs(34) * inda * indb) + (sizecoefs(35) * inda * indc) + 
                (sizecoefs(36) * indb * indc) + (sizecoefs(37) * inda * indb) + 
                (sizecoefs(38) * inda * indc) + (sizecoefs(39) * indb * indc) + 
                (sizecoefs(40) * inda * indb) + (sizecoefs(41) * inda * indb) + 
                (sizecoefs(42) * inda * indc) + (sizecoefs(43) * inda * indc) + 
                (sizecoefs(44) * indb * indc) + (sizecoefs(45) * indb * indc) + 
                sizepatch(patchnumber) + sizeyear(yearnumber) + sizedev);
              
              out(i, 3) = (exp(-1 * (pow((sz3(i) - preout(3)), 2) / (2*sigma2))) / 
                ((pow((2 * M_PI), 0.5)) * sigma)) * binwidth3(i);
              
            } else {
              out(i, 3) = sizedist - 3;
            }
          }
          
          if (repstl > 1) {
            
            preout(2) = (repstcoefs(0) + (repstcoefs(1) * fl1(i)) + (repstcoefs(2) * fl2n(i)) + 
              (repstcoefs(3) * sz1(i)) + (repstcoefs(4) * sz2o(i)) + (repstcoefs(5) * fl2n(i) * fl1(i)) + 
              (repstcoefs(6) * sz2o(i) * sz1(i)) + (repstcoefs(7) * sz1(i) * fl1(i)) + 
              (repstcoefs(8) * sz2o(i) * fl2n(i)) + (repstcoefs(9) * sz2o(i) * fl1(i)) + 
              (repstcoefs(10) * sz1(i) * fl2n(i)) + (repstcoefs(11) * actualage2(i)) + 
              (repstcoefs(12) * actualage2(i) * sz1(i)) + (repstcoefs(13) * actualage2(i) * sz2o(i)) + 
              (repstcoefs(14) * actualage2(i) * fl1(i)) + (repstcoefs(15) * actualage2(i) * fl2n(i)) + 
              (repstcoefs(16) * inda) + (repstcoefs(17) * indb) + (repstcoefs(18) * indc) + 
              (repstcoefs(19) * inda) + (repstcoefs(20) * indb) + (repstcoefs(21) * indc) + 
              (repstcoefs(22) * inda * sz2o(i)) + (repstcoefs(23) * indb * sz2o(i)) + 
              (repstcoefs(24) * indc * sz2o(i)) + (repstcoefs(25) * inda * fl2n(i)) + 
              (repstcoefs(26) * indb * fl2n(i)) + (repstcoefs(27) * indc * fl2n(i)) + 
              (repstcoefs(28) * inda * sz1(i)) + (repstcoefs(29) * indb * sz1(i)) + 
              (repstcoefs(30) * indc * sz1(i)) + (repstcoefs(31) * inda * fl1(i)) + 
              (repstcoefs(32) * indb * fl1(i)) + (repstcoefs(33) * indc * fl1(i)) + 
              (repstcoefs(34) * inda * indb) + (repstcoefs(35) * inda * indc) +
              (repstcoefs(36) * indb * indc) + (repstcoefs(37) * inda * indb) + 
              (repstcoefs(38) * inda * indc) + (repstcoefs(39) * indb * indc) + 
              (repstcoefs(40) * inda * indb) + (repstcoefs(41) * inda * indb) +
              (repstcoefs(42) * inda * indc) + (repstcoefs(43) * inda * indc) + 
              (repstcoefs(44) * indb * indc) + (repstcoefs(45) * indb * indc) + 
              repstpatch(patchnumber) + repstyear(yearnumber) + repstdev);
            
            if (preout(2) > 700) preout(2) = 700;
            
            out(i, 2) = exp(preout(2)) / (1 + exp(preout(2)));
            
            if (fl3(i) == 0) {
              out(i, 2) = 1 - out(i, 2);
            }
            
          } else {
            if (fl3(i) == 0) {
              out(i, 2) = 1 - repstcoefs(0);
            } else if (fl3(i) == 1) {
              out(i, 2) = repstcoefs(0);
            } else {
              out(i, 2) = 0;
            }
          }
          
        } else {
          out(i, 1) = 1 - out(i, 1);
          out(i, 3) = 1;
          out(i, 2) = 1;
        }
        
        survtransmat(k) = out(i, 0) * out(i, 1) * out(i, 2) * out(i, 3);
        
        
      } else if (immat2n(i) == 1 && immat1(i) == 1 && jsurvl > 0) {
        // Juvenile to adult transitions
        
        arma::vec preout(4);
        
        if (jsurvl > 1) {
          preout(0) = (jsurvcoefs(0) + jsurvpatch(patchnumber) + jsurvyear(yearnumber) + jsurvdev);
          
          if (preout(0) > 700) preout(0) = 700;
          
          out(i, 0) = exp(preout(0)) / (1 + exp(preout(0)));
        } else {
          out(i, 0) = jsurvcoefs(0);
        }
        
        if (jobsl > 1) {
          preout(1) = (jobscoefs(0) + jobspatch(patchnumber) + jobsyear(yearnumber) + jobsdev);
          
          if (preout(1) > 700) preout(1) = 700;
          
          out(i, 1) = exp(preout(1)) / (1 + exp(preout(1)));
          
        } else {
          out(i, 1) = jobscoefs(0);
        }
        
        if (ob3(i) == 1 || jobsl == 1) {
          if (jsizel > 1) {
            if (sizedist == 0) {
              // Poisson size distribution
              
              if (jsizezero == TRUE && sz3(i) == 0) {
                lambda_preexp = (jsizecoefs(46) + jsizepatchzi(patchnumber) + 
                  jsizeyearzi(yearnumber) + jsizedev + (jsummedvars / 2));
                
                if (lambda_preexp > 700) {
                  lambda = exp(static_cast<double>(700));
                } else lambda = exp(lambda_preexp);
                
                out(i, 3) = (lambda) / (1 + (lambda));
                
              } else {
                double sizefac = sz3(i) * tgamma(sz3(i));
                lambda_preexp = (jsizecoefs(0) + jsizepatch(patchnumber) + 
                  jsizeyear(yearnumber) + jsizedev + (jsummedvars / 2));
                
                if (lambda_preexp > 500) {
                  lambda = exp(static_cast<double>(500));
                } else lambda = exp(lambda_preexp);
                
                if (jsizetrunc == 1) {
                  out(i, 3) = ((pow(lambda, sz3(i)) * exp(-1 * lambda)) / sizefac) / (1 - (exp(-1 * lambda)));
                } else {
                  out(i, 3) = ((pow(lambda, sz3(i)) * exp(-1 * lambda)) / sizefac);
                }
              }
            } else if (sizedist == 1) {
              // Negative binomial size distribution
              
              if (jsizezero == TRUE && sz3(i) == 0) {
                mu_preexp = (jsizecoefs(46) + jsizepatchzi(patchnumber) + jsizeyearzi(yearnumber) + jsizedev);
                
                if (mu_preexp > 700) {
                  mu = exp(static_cast<double>(700));
                } else mu = exp(mu_preexp);
                
                out(i, 3) = (mu) / (1 + (mu));
              } else {
                mu = (jsizecoefs(0) + jsizepatch(patchnumber) + jsizeyear(yearnumber) + jsizedev);
                
                if (mu_preexp > 500) {
                  mu = exp(static_cast<double>(500));
                } else mu = exp(mu_preexp);
                
                double theta = jsizesigma;
                if (theta > 1000000000) {
                  theta = 1000000000;
                }
                double alpha = 1 / theta;
                
                int y = static_cast<int>(sz3(i));
                
                double log_leftie = 0;
                for (int j = 0; j < y; j++) {
                  log_leftie = log(static_cast<double>(j) + theta) - log(static_cast<double>(j) + 1) + log_leftie;
                }
                
                double log_amu = log(alpha) + log(mu);
                double log_mid = -1 * theta * log(1 + (alpha * mu));
                
                double log_rightie = sz3(i) * (log_amu - log(1 + (alpha * mu)));
                
                double raw_prob = log_leftie + log_mid + log_rightie;
                
                if (jsizetrunc == 1) {
                  double zero_raw_prob = log_mid;
                  
                  out(i, 3) = exp(raw_prob) / (1 - exp(zero_raw_prob));
                } else {
                  out(i, 3) = exp(raw_prob);
                }
              }
            } else if (sizedist == 2) {
              // Gaussian size distribution
              
              double sigma2 = jsigma * jsigma;
              preout(3) = (jsizecoefs(0) + jsizepatch(patchnumber) + jsizeyear(yearnumber) + jsizedev);
              
              out(i, 3) = (exp(-1 * (pow((sz3(i) - preout(3)), 2) / (2*sigma2))) / 
                ((pow((2 * M_PI), 0.5)) * jsigma)) * binwidth3(i);
              
            } else {
              out(i, 3) = sizedist - 3;
            }
          }
          
          if (jrepstl > 1) {
            
            preout(2) = (jrepstcoefs(0) + jrepstpatch(patchnumber) + jrepstyear(yearnumber) + jrepstdev);
            
            if (preout(2) > 700) preout(2) = 700;
            
            out(i, 2) = exp(preout(2)) / (1 + exp(preout(2)));
            
            if (fl3(i) == 0) {
              out(i, 2) = 1 - out(i, 2);
            }
            
          } else {
            if (fl3(i) == 0) {
              out(i, 2) = 1 - jrepstcoefs(0);
            } else if (fl3(i) == 1) {
              out(i, 2) = jrepstcoefs(0);
            } else {
              out(i, 2) = 0;
            }
          }
          
        } else {
          out(i, 1) = 1 - out(i, 1);
          out(i, 3) = 1;
          out(i, 2) = 1;
        }
        
        survtransmat(k) = out(i, 0) * out(i, 1) * out(i, 2) * out(i, 3);
        
      }
    } else if (ovgivent(i) != -1) {
      // All other transitions
      
      survtransmat(k) = ovgivent(i);
    }
    
    // This next block calculates fecundity
    if (indata2n(i) == 1 && fecl > 0) {
      if (fl2o(i) == 1 && ovgivenf(i) == -1) {
        
        double preoutx;
        
        if (matrixformat != 2) {
          if (fecdist < 3 && fecl > 1) {
            if (feczero == TRUE && sz3(i) == 0) {
              
              preoutx = (feccoefs(46) + (feccoefs(47) * fl1(i)) + (feccoefs(48) * fl2n(i)) +
                (feccoefs(49) * sz1(i)) + (feccoefs(50) * sz2o(i)) + 
                (feccoefs(51) * fl2n(i) * fl1(i)) + (feccoefs(52) * sz2o(i) * sz1(i)) +
                (feccoefs(53) * sz1(i) * fl1(i)) + (feccoefs(54) * sz2o(i) * fl2n(i)) + 
                (feccoefs(55) * sz2o(i) * fl1(i)) + (feccoefs(56) * sz1(i) * fl2n(i)) + 
                (feccoefs(57) * actualage2(i)) + (feccoefs(58) * actualage2(i) * sz1(i)) + 
                (feccoefs(59) * actualage2(i) * sz2o(i)) + (feccoefs(60) * actualage2(i) * fl1(i)) + 
                (feccoefs(61) * actualage2(i) * fl2n(i)) + (feccoefs(62) * inda) +
                (feccoefs(63) * indb) + (feccoefs(64) * indc) + (feccoefs(65) * inda) +
                (feccoefs(66) * indb) + (feccoefs(67) * indc) + (feccoefs(68) * inda * sz2o(i)) + 
                (feccoefs(69) * indb * sz2o(i)) + (feccoefs(70) * indc * sz2o(i)) + 
                (feccoefs(71) * inda * fl2n(i)) + (feccoefs(72) * indb * fl2n(i)) + 
                (feccoefs(73) * indc * fl2n(i)) + (feccoefs(74) * inda * sz1(i)) + 
                (feccoefs(75) * indb * sz1(i)) + (feccoefs(76) * indc * sz1(i)) + 
                (feccoefs(77) * inda * fl1(i)) + (feccoefs(78) * indb * fl1(i)) + 
                (feccoefs(79) * indc * fl1(i)) + (feccoefs(80) * inda * indb) + 
                (feccoefs(81) * inda * indc) + (feccoefs(82) * indb * indc) + 
                (feccoefs(83) * inda * indb) + (feccoefs(84) * inda * indc) + 
                (feccoefs(85) * indb * indc) + (feccoefs(86) * inda * indb) + 
                (feccoefs(87) * inda * indb) + (feccoefs(88) * inda * indc) + 
                (feccoefs(89) * inda * indc) + (feccoefs(90) * indb * indc) + 
                (feccoefs(91) * indb * indc) + fecpatchzi(patchnumber) + 
                fecyearzi(yearnumber) + fecdev);
                
            } else {
              
              preoutx = (feccoefs(0) + (feccoefs(1) * fl1(i)) + (feccoefs(2) * fl2n(i)) + 
                (feccoefs(3) * sz1(i)) + (feccoefs(4) * sz2o(i)) + 
                (feccoefs(5) * fl2n(i) * fl1(i)) + (feccoefs(6) * sz2o(i) * sz1(i)) +
                (feccoefs(7) * sz1(i) * fl1(i)) + (feccoefs(8) * sz2o(i) * fl2n(i)) +
                (feccoefs(9) * sz2o(i) * fl1(i)) + (feccoefs(10) * sz1(i) * fl2n(i)) +
                (feccoefs(11) * actualage2(i)) + (feccoefs(12) * actualage2(i) * sz1(i)) + 
                (feccoefs(13) * actualage2(i) * sz2o(i)) + (feccoefs(14) * actualage2(i) * fl1(i)) + 
                (feccoefs(15) * actualage2(i) * fl2n(i)) + (feccoefs(16) * inda) + (feccoefs(17) * indb) + 
                (feccoefs(18) * indc) + (feccoefs(19) * inda) + (feccoefs(20) * indb) +
                (feccoefs(21) * indc) + (feccoefs(22) * inda * sz2o(i)) +
                (feccoefs(23) * indb * sz2o(i)) + (feccoefs(24) * indc * sz2o(i)) + 
                (feccoefs(25) * inda * fl2n(i)) + (feccoefs(26) * indb * fl2n(i)) +
                (feccoefs(27) * indc * fl2n(i)) + (feccoefs(28) * inda * sz1(i)) + 
                (feccoefs(29) * indb * sz1(i)) + (feccoefs(30) * indc * sz1(i)) + 
                (feccoefs(31) * inda * fl1(i)) + (feccoefs(32) * indb * fl1(i)) + 
                (feccoefs(33) * indc * fl1(i)) + (feccoefs(34) * inda * indb) + 
                (feccoefs(35) * inda * indc) + (feccoefs(36) * indb * indc) + 
                (feccoefs(37) * inda * indb) + (feccoefs(38) * inda * indc) +
                (feccoefs(39) * indb * indc) + (feccoefs(40) * inda * indb) + 
                (feccoefs(41) * inda * indb) + (feccoefs(42) * inda * indc) + 
                (feccoefs(43) * inda * indc) + (feccoefs(44) * indb * indc) + 
                (feccoefs(45) * indb * indc) + fecpatch(patchnumber) + 
                fecyear(yearnumber) + fecdev);
            }
            
            if (fecdist != 2) {
              // Poisson and negative binomial fecundity
              
              if (feczero == TRUE && sz3(i) == 0) {
                
                if (preoutx > 700) preoutx = 700;
                
                fectransmat(k) = (exp(preoutx) / (1+exp(preoutx))) * fecmod * repentry(i);
                
              } else {
              
                if (preoutx > 700) preoutx = 700;
                
                fectransmat(k) = exp(preoutx) * fecmod * repentry(i);
              }
            } else {
              // Gaussian fecundity
              fectransmat(k) = preoutx * fecmod * repentry(i);
              
              if (negfec && fectransmat(k) < 0) {
                fectransmat(k) = 0;
              }
            }
            
          } else if (fecl > 1) {
            // All others with estimated models
            
            fectransmat(k) = fecdist - 3;
          } else {
            fectransmat(k) = feccoefs(0);
          }
        } else if (stage2n(i) == (nostages+1)) {
          // This propagates fecundity in deVries-formatted hMPMs
          if (fecdist < 3 && fecl > 1) {
            if (feczero == TRUE && sz3(i) == 0) {
              
              preoutx = (feccoefs(46) + (feccoefs(47) * fl1(i)) + (feccoefs(48) * fl2o(i)) +
                (feccoefs(49) * sz1(i)) + (feccoefs(50) * sz2o(i)) + 
                (feccoefs(51) * fl2o(i) * fl1(i)) + (feccoefs(52) * sz2o(i) * sz1(i)) +
                (feccoefs(53) * sz1(i) * fl1(i)) + (feccoefs(54) * sz2o(i) * fl2o(i)) + 
                (feccoefs(55) * sz2o(i) * fl1(i)) + (feccoefs(56) * sz1(i) * fl2o(i)) + 
                (feccoefs(57) * actualage2(i)) + (feccoefs(58) * actualage2(i) * sz1(i)) + 
                (feccoefs(59) * actualage2(i) * sz2o(i)) + (feccoefs(60) * actualage2(i) * fl1(i)) + 
                (feccoefs(61) * actualage2(i) * fl2o(i)) + (feccoefs(62) * inda) +
                (feccoefs(63) * indb) + (feccoefs(64) * indc) + (feccoefs(65) * inda) +
                (feccoefs(66) * indb) + (feccoefs(67) * indc) + (feccoefs(68) * inda * sz2o(i)) + 
                (feccoefs(69) * indb * sz2o(i)) + (feccoefs(70) * indc * sz2o(i)) + 
                (feccoefs(71) * inda * fl2o(i)) + (feccoefs(72) * indb * fl2o(i)) + 
                (feccoefs(73) * indc * fl2o(i)) + (feccoefs(74) * inda * sz1(i)) + 
                (feccoefs(75) * indb * sz1(i)) + (feccoefs(76) * indc * sz1(i)) + 
                (feccoefs(77) * inda * fl1(i)) + (feccoefs(78) * indb * fl1(i)) + 
                (feccoefs(79) * indc * fl1(i)) + (feccoefs(80) * inda * indb) + 
                (feccoefs(81) * inda * indc) + (feccoefs(82) * indb * indc) + 
                (feccoefs(83) * inda * indb) + (feccoefs(84) * inda * indc) + 
                (feccoefs(85) * indb * indc) + (feccoefs(86) * inda * indb) + 
                (feccoefs(87) * inda * indb) + (feccoefs(88) * inda * indc) + 
                (feccoefs(89) * inda * indc) + (feccoefs(90) * indb * indc) + 
                (feccoefs(91) * indb * indc) + fecpatchzi(patchnumber) + 
                fecyearzi(yearnumber) + fecdev);
                
            } else {
              
              preoutx = (feccoefs(0) + (feccoefs(1) * fl1(i)) + (feccoefs(2) * fl2o(i)) + 
                (feccoefs(3) * sz1(i)) + (feccoefs(4) * sz2o(i)) + 
                (feccoefs(5) * fl2o(i) * fl1(i)) + (feccoefs(6) * sz2o(i) * sz1(i)) +
                (feccoefs(7) * sz1(i) * fl1(i)) + (feccoefs(8) * sz2o(i) * fl2o(i)) +
                (feccoefs(9) * sz2o(i) * fl1(i)) + (feccoefs(10) * sz1(i) * fl2o(i)) +
                (feccoefs(11) * actualage2(i)) + (feccoefs(12) * actualage2(i) * sz1(i)) + 
                (feccoefs(13) * actualage2(i) * sz2o(i)) + (feccoefs(14) * actualage2(i) * fl1(i)) + 
                (feccoefs(15) * actualage2(i) * fl2o(i)) + (feccoefs(16) * inda) + (feccoefs(17) * indb) + 
                (feccoefs(18) * indc) + (feccoefs(19) * inda) + (feccoefs(20) * indb) +
                (feccoefs(21) * indc) + (feccoefs(22) * inda * sz2o(i)) +
                (feccoefs(23) * indb * sz2o(i)) + (feccoefs(24) * indc * sz2o(i)) + 
                (feccoefs(25) * inda * fl2o(i)) + (feccoefs(26) * indb * fl2o(i)) +
                (feccoefs(27) * indc * fl2o(i)) + (feccoefs(28) * inda * sz1(i)) + 
                (feccoefs(29) * indb * sz1(i)) + (feccoefs(30) * indc * sz1(i)) + 
                (feccoefs(31) * inda * fl1(i)) + (feccoefs(32) * indb * fl1(i)) + 
                (feccoefs(33) * indc * fl1(i)) + (feccoefs(34) * inda * indb) + 
                (feccoefs(35) * inda * indc) + (feccoefs(36) * indb * indc) + 
                (feccoefs(37) * inda * indb) + (feccoefs(38) * inda * indc) +
                (feccoefs(39) * indb * indc) + (feccoefs(40) * inda * indb) + 
                (feccoefs(41) * inda * indb) + (feccoefs(42) * inda * indc) + 
                (feccoefs(43) * inda * indc) + (feccoefs(44) * indb * indc) + 
                (feccoefs(45) * indb * indc) + fecpatch(patchnumber) + 
                fecyear(yearnumber) + fecdev);
            }
            
            if (fecdist != 2) {
              // Poisson and negative binomial fecundity
              
              if (feczero == TRUE && sz3(i) == 0) {
                
                if (preoutx > 700) preoutx = 700;
                
                fectransmat(k) = (exp(preoutx) / (1+exp(preoutx))) * fecmod * repentry(i);
                
              } else {
              
                if (preoutx > 700) preoutx = 700;
                
                fectransmat(k) = exp(preoutx) * fecmod * repentry(i);
              }
            } else {
              // Gaussian fecundity
              fectransmat(k) = preoutx * fecmod * repentry(i);
              
              if (negfec && fectransmat(k) < 0) {
                fectransmat(k) = 0;
              }
            }
            
          } else if (fecl > 1) {
            // All others with estimated models
            
            fectransmat(k) = fecdist - 3;
          } else {
            fectransmat(k) = feccoefs(0);
          }
        }

      } else if (ovgivenf(i) != -1 ) {
        fectransmat(k) = ovgivenf(i);
      }
    } else if (ovgivenf(i) != -1 ) {
      fectransmat(k) = ovgivenf(i);
    }
  }
  
  if (replacementst > 0) {
    for (int i = 0; i < replacementst; i++) {
      repindex = replacetvec(i); // AllStages index
      properindex = aliveandequal(repindex);
      arma::uvec rightindex = find(index321 == ovestt(repindex));
      proxyindex = aliveandequal(rightindex(0));
      survtransmat(properindex) = survtransmat(proxyindex);
    }
  }
  
  if (replacementsf > 0) {
    for (int i = 0; i < replacementsf; i++) {
      repindex = replacefvec(i); // AllStages index
      properindex = aliveandequal(repindex);
      arma::uvec rightindex = find(index321 == ovestf(repindex));
      proxyindex = aliveandequal(rightindex(0));
      fectransmat(properindex) = fectransmat(proxyindex);
    }
  }
  
  arma::mat amatrix = survtransmat + fectransmat;
  
  return List::create(Named("A") = amatrix, _["U"] = survtransmat,
    _["F"] = fectransmat, _["out"] = out);
}

