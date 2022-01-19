// [[Rcpp::depends(RcppArmadillo)]]
#define BOOST_DISABLE_ASSERTS

#include <RcppArmadillo.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>

using namespace Rcpp;
using namespace arma;

//' Estimate All Elements of Raw Historical Matrix
//' 
//' Function \code{.specialpatrolgroup()} swiftly calculates matrix transitions
//' in raw historical matrices, and serves as the core workhorse function behind
//' \code{\link{rlefko3}()}.
//' 
//' @param sge9l The Allstages data frame developed for \code{rlefko3()}
//' covering stage pairs across times \emph{t}+1, \emph{t} and \emph{t}-1.
//' Generally termed \code{stageexpansion9}.
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
// [[Rcpp::export(.specialpatrolgroup)]]
List specialpatrolgroup(DataFrame sge9l, DataFrame sge3, DataFrame MainData,
  DataFrame StageFrame, int format, int err_switch) {
  
  arma::vec sge9stage3 = sge9l["stage3"];
  arma::vec sge9fec32 = sge9l["repentry3"];
  arma::vec sge9rep2 = sge9l["rep2o"];
  arma::vec sge9indata32 = sge9l["indata"];
  arma::vec sge9ovgivent = sge9l["ovgiven_t"];
  arma::vec sge9ovgivenf = sge9l["ovgiven_f"];
  arma::vec sge9ovestt = sge9l["ovest_t"];
  arma::vec sge9ovestf = sge9l["ovest_f"];
  arma::vec sge9ovsurvmult = sge9l["ovsurvmult"];
  arma::vec sge9ovfecmult = sge9l["ovfecmult"];
  arma::vec sge9index321 = sge9l["index321"];
  arma::vec sge9index21 = sge9l["index21"]; // This is sge92index - not sure if this is needed
  arma::vec aliveandequal = sge9l["aliveandequal"];
  
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

  arma::vec sfsizes = StageFrame["sizebin_center"];
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
  
  arma::mat stage2fec(sge3index21.n_elem, 3, fill::zeros); //1st col = # indivs total, 2nd col = no indivs alive, 3rd col = sum fec
  
  // These next structures develop the prior stage
  arma::vec probsrates0p(noelems, fill::zeros); // 1st vec = # indivs (3 trans)
  arma::vec probsrates1p(noelems, fill::zeros); // 2nd vec = total indivs for pair stage
  arma::vec probsrates2p(noelems, fill::zeros); // 3rd vec = total indivs alive for pair stage
  arma::vec probsrates3p(noelems, fill::zeros); // 4th vec = total fec for pair stage
  
  arma::mat stage2fecp(sge3index21.n_elem, 3, fill::zeros); //1st col = # indivs total, 2nd col = no indivs alive, 3rd col = sum fec
  
  // The final matrices, though empty
  arma::mat tmatrix(matrixdim, matrixdim, fill::zeros); // Main output U matrix
  arma::mat fmatrix(matrixdim, matrixdim, fill::zeros); // Main output F matrix
  
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
        arma::uvec choiceelementp = find(sge9index321 == dataindex321_prior(i, j));
        
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
      if (sge9ovsurvmult(elem3) < 0) sge9ovsurvmult(elem3) = 1.0;
      
      tmatrix(aliveandequal(elem3)) = probsrates0(elem3)* sge9ovsurvmult(elem3) /
        probsrates1(elem3); // Survival
      
      // Fecundity
      if (sge9ovfecmult(elem3) < 0) sge9ovfecmult(elem3) = 1.0;
      if (format == 2) {
        fmatrix(aliveandequal(elem3)) = sge9fec32(elem3) * sge9rep2(elem3) *
          probsrates3p(elem3) * sge9ovfecmult(elem3) / probsrates1p(elem3);
      } else {
        fmatrix(aliveandequal(elem3)) = sge9fec32(elem3) * sge9rep2(elem3) *
          probsrates3(elem3) * sge9ovfecmult(elem3) / probsrates1(elem3);
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
  
  // The next bit changes NAs to 0
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
//' Function \code{.normalpatrolgroup()} swiftly calculates matrix transitions
//' in raw ahistorical matrices, and serves as the core workhorse function
//' behind \code{\link{rlefko2}()}.
//' 
//' @param sge3 The Allstages data frame developed for \code{rlefko2()} covering
//' stage pairs across times \emph{t}+1 and \emph{t}. Generally termed
//' \code{stageexpansion3}.
//' @param sge2 The data frame covering all stages in time \emph{t}. Generally
//' termed \code{stageexpansion2}.
//' @param MainData The demographic dataset modified to hold \code{usedfec} and
//' \code{usedstage} columns.
//' @param StageFrame The full stageframe for the analysis.
//' 
//' @return List of three matrices, including the survival-transition (U)
//' matrix, the fecundity matrix (F), and the sum (A) matrix, with A first.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.normalpatrolgroup)]]
List normalpatrolgroup(DataFrame sge3, DataFrame sge2, DataFrame MainData,
  DataFrame StageFrame) {
  
  arma::vec sge3fec32 = sge3["repentry3"];
  arma::vec sge3rep2 = sge3["rep2n"];
  arma::vec sge3indata32 = sge3["indata"];
  arma::vec sge3ovgivent = sge3["ovgiven_t"];
  arma::vec sge3ovgivenf = sge3["ovgiven_f"];
  arma::vec sge3ovestt = sge3["ovest_t"];
  arma::vec sge3ovestf = sge3["ovest_f"];
  arma::vec sge3ovsurvmult = sge3["ovsurvmult"];
  arma::vec sge3ovfecmult = sge3["ovfecmult"];
  arma::vec sge3index32 = sge3["index321"];
  arma::vec sge3index2 = sge3["stage2n"];
  arma::vec aliveandequal = sge3["aliveandequal"];
  
  arma::vec sge2rep2 = sge2["rep2"];
  arma::vec sge2fec3 = sge2["fec3"];
  arma::vec sge2index2 = sge2["index2"];
  arma::vec sge2stage2 = sge2["stage2"];
  
  arma::vec dataindex32 = MainData["index32"];
  arma::vec dataindex2 = MainData["index2"];
  arma::vec dataalive3 = MainData["alive3"];
  arma::vec datausedfec2 = MainData["usedfec2"];
  
  arma::vec sfsizes = StageFrame["sizebin_center"];
  int nostages = sfsizes.n_elem;
  
  int n = dataindex32.n_elem;
  int no2stages = sge2index2.n_elem - 1; // The -1 removes the dead stage, which is still within sge2
  int noelems = sge3index32.n_elem;
  
  arma::vec probsrates0(noelems, fill::zeros); // 1st vec = # indivs (3 trans)
  arma::vec probsrates1(noelems, fill::zeros); // 2nd vec = total indivs for pair stage
  arma::vec probsrates2(noelems, fill::zeros); // 3rd vec = total indivs alive for pair stage
  arma::vec probsrates3(noelems, fill::zeros); // 4th vec = total fec for pair stage
  
  arma::mat stage2fec(no2stages, 3, fill::zeros); //1st col = # indivs total, 2nd col = no indivs alive, 3rd col = sum fec
  
  arma::mat tmatrix((nostages-1), (nostages-1), fill::zeros); // Main output U matrix
  arma::mat fmatrix((nostages-1), (nostages-1), fill::zeros); // Main output F matrix
  
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
      if (sge3ovsurvmult(elem3) < 0) sge3ovsurvmult(elem3) = 1.0;
      tmatrix(aliveandequal(elem3)) = probsrates0(elem3) * sge3ovsurvmult(elem3) / 
        probsrates1(elem3);
        
      if (sge3ovfecmult(elem3) < 0) sge3ovfecmult(elem3) = 1.0;
      fmatrix(aliveandequal(elem3)) = sge3fec32(elem3) * sge3rep2(elem3) * 
        probsrates3(elem3) * sge3ovfecmult(elem3) / probsrates1(elem3);
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

//' Estimate All Elements of Raw Ahistorical Population Projection Matrix
//' 
//' Function \code{.minorpatrolgroup()} swiftly calculates matrix transitions
//' in raw Leslie MPMs, and is used internally in \code{\link{rleslie}()}.
//' 
//' @param MainData The demographic dataset modified internally to have needed
//' variables for living status, reproduction status, and fecundity.
//' @param StageFrame The full stageframe for the analysis.
//' @param fectime An integer coding to estimate fecundity using time \emph{t}
//' (\code{2}) or time \emph{t}+1 \code{(3)}.
//' @param cont Should a self-loop transition be estimated for the final age.
//' @param lastage An integer coding for the last age to use in matrix
//' construction.
//' 
//' @return List of three matrices, including the survival-transition (U)
//' matrix, the fecundity matrix (F), and the sum (A) matrix, with A first.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.minorpatrolgroup)]]
Rcpp::List minorpatrolgroup(DataFrame MainData, DataFrame StageFrame,
  int fectime, bool cont, double fec_mod) {
  
  arma::ivec data_age = MainData["usedobsage"];
  arma::vec data_alive2 = MainData["usedalive2"];
  arma::vec data_alive3 = MainData["usedalive3"];
  arma::vec data_usedfec = MainData["usedfec2"];
  arma::vec data_usedfec3 = MainData["usedfec3"];
  arma::vec data_usedrepst2 = MainData["usedrepst2"];
  arma::vec data_usedrepst3 = MainData["usedrepst3"];

  if (fectime == 3) {
    data_usedfec = data_usedfec3;
  }
  
  IntegerVector sf_minage = StageFrame["min_age"];
  IntegerVector sf_maxage = StageFrame["max_age"];
  IntegerVector sf_repstatus = StageFrame["repstatus"];
  int noages = sf_minage.length();
  
  arma::vec probsrates0(noages, fill::zeros); // 1st vec = total indivs in t
  arma::vec probsrates1(noages, fill::zeros); // 2nd vec = total indivs alive in t+1
  arma::vec probsrates2(noages, fill::zeros); // 3rd vec = total fec
  
  arma::mat tmatrix(noages, noages, fill::zeros); // Main output U matrix
  arma::mat fmatrix(noages, noages, fill::zeros); // Main output F matrix
  
  // This main loop counts individuals going through transitions and calculates
  // survival and fecundity
  arma::uvec data_allalive = find(data_alive2);
  int survsum {0};
  double fecsum {0};
  
  for (int i = 0; i < noages; i++) { 
    
    arma::uvec data_indices = find(data_age == sf_minage(i));
    arma::uvec aget_alive = intersect(data_allalive, data_indices);
    int num_aget_alive = aget_alive.n_elem;
    
    if (num_aget_alive > 0) {
      for (int j = 0; j < num_aget_alive; j++) {
        if (data_alive3(aget_alive(j)) > 0) survsum++;
        fecsum = fecsum + data_usedfec(aget_alive(j));
      }
    }
    
    probsrates0(i) = num_aget_alive;
    if (num_aget_alive > 0) {
      probsrates1(i) = static_cast<double>(survsum) / static_cast<double>(num_aget_alive);
      if (sf_repstatus(i) > 0) probsrates2(i) = fec_mod * fecsum / static_cast<double>(num_aget_alive);
    } else {
      probsrates1(i) = 0;
      probsrates2(i) = 0;
    }
    
    if (i < (noages - 1)) tmatrix(i+1, i) = probsrates1(i);
    fmatrix(0, i) = probsrates2(i);
    
    survsum = 0;
    fecsum = 0;
  }
  
  tmatrix(find_nonfinite(tmatrix)).zeros();
  fmatrix(find_nonfinite(fmatrix)).zeros();
  
  arma::mat amatrix = tmatrix + fmatrix; // Create the A matrix
  
  return List::create(Named("A") = amatrix, _["U"] = tmatrix, _["F"] = fmatrix);
}

//' Creates Matrices of Year and Patch Terms in Models
//' 
//' Function \code{.revelations()} creates a matrix holding either the year or
//' patch coefficients from all vital rate models. This reduces memory load in
//' functions \code{\link{.jerzeibalowski}()}, which may be important in some
//' systems or compilers.
//' 
//' @param survproxy The proxy vital rate model covering survival from the main
//' matrix estimator function.
//' @param obsproxy The proxy vital rate model covering observation status from
//' the main matrix estimator function.
//' @param sizeproxy The proxy vital rate model covering primary size from the
//' main matrix estimator function.
//' @param sizebproxy The proxy vital rate model covering secondary size from
//' the main matrix estimator function.
//' @param sizecproxy The proxy vital rate model covering tertiary size from the
//' main matrix estimator function.
//' @param repstproxy The proxy vital rate model covering reproductive status
//' from the main matrix estimator function.
//' @param fecproxy The proxy vital rate model covering fecundity from the main
//' matrix estimator function.
//' @param jsurvproxy The proxy vital rate model covering juvenile survival from
//' the main matrix estimator function.
//' @param jobsproxy The proxy vital rate model covering juvenile observation
//' status from the main matrix estimator function.
//' @param jsizeproxy The proxy vital rate model covering juvenile primary size
//' from the main matrix estimator function.
//' @param jsizebproxy The proxy vital rate model covering juvenile secondary
//' size from the main matrix estimator function.
//' @param jsizecproxy The proxy vital rate model covering juvenile tertiary
//' size from the main matrix estimator function.
//' @param jrepstproxy The proxy vital rate model covering juvenile reproductive
//' status from the main matrix estimator function.
//' @param jmatstproxy The proxy vital rate model covering juvenile probability
//' of becoming mature from the main matrix estimator function.
//' @param mat_switch An integer coding for year (\code{1}) or patch (\code{2}).
//' 
//' @return A matrix with 14 columns corresponding to the number of vital rates
//' and number of columns equal to the number of year or patches.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.revelations)]]
NumericMatrix revelations(List survproxy, List obsproxy, List sizeproxy,
  List sizebproxy, List sizecproxy, List repstproxy, List fecproxy,
  List jsurvproxy, List jobsproxy, List jsizeproxy, List jsizebproxy,
  List jsizecproxy, List jrepstproxy,List jmatstproxy, int mat_switch) {
  
  NumericMatrix final_mat;
  
  if (mat_switch == 1) {
    Rcpp::DataFrame survyear_df(survproxy["years"]);
    Rcpp::DataFrame obsyear_df(obsproxy["years"]);
    Rcpp::DataFrame sizeyear_df(sizeproxy["years"]);
    Rcpp::DataFrame sizebyear_df(sizebproxy["years"]);
    Rcpp::DataFrame sizecyear_df(sizecproxy["years"]);
    Rcpp::DataFrame repstyear_df(repstproxy["years"]);
    Rcpp::DataFrame fecyear_df(fecproxy["years"]);
    Rcpp::DataFrame jsurvyear_df(jsurvproxy["years"]);
    Rcpp::DataFrame jobsyear_df(jobsproxy["years"]);
    Rcpp::DataFrame jsizeyear_df(jsizeproxy["years"]);
    Rcpp::DataFrame jsizebyear_df(jsizebproxy["years"]);
    Rcpp::DataFrame jsizecyear_df(jsizecproxy["years"]);
    Rcpp::DataFrame jrepstyear_df(jrepstproxy["years"]);
    Rcpp::DataFrame jmatstyear_df(jmatstproxy["years"]);
    
    NumericVector survyear = survyear_df[0];
    NumericVector obsyear = obsyear_df[0];
    NumericVector sizeyear = sizeyear_df[0];
    NumericVector sizebyear = sizebyear_df[0];
    NumericVector sizecyear = sizecyear_df[0];
    NumericVector repstyear = repstyear_df[0];
    NumericVector fecyear = fecyear_df[0];
    NumericVector jsurvyear = jsurvyear_df[0];
    NumericVector jobsyear = jobsyear_df[0];
    NumericVector jsizeyear = jsizeyear_df[0];
    NumericVector jsizebyear = jsizebyear_df[0];
    NumericVector jsizecyear = jsizecyear_df[0];
    NumericVector jrepstyear = jrepstyear_df[0];
    NumericVector jmatstyear = jmatstyear_df[0];
    
    int matrows = survyear.length();
    
    NumericMatrix year_mat(matrows, 14);
    year_mat(_, 0) = survyear;
    year_mat(_, 1) = obsyear;
    year_mat(_, 2) = sizeyear;
    year_mat(_, 3) = sizebyear;
    year_mat(_, 4) = sizecyear;
    year_mat(_, 5) = repstyear;
    year_mat(_, 6) = fecyear;
    year_mat(_, 7) = jsurvyear;
    year_mat(_, 8) = jobsyear;
    year_mat(_, 9) = jsizeyear;
    year_mat(_, 10) = jsizebyear;
    year_mat(_, 11) = jsizecyear;
    year_mat(_, 12) = jrepstyear;
    year_mat(_, 13) = jmatstyear;
    
    final_mat = year_mat;
    
  } else if (mat_switch == 2) {
    
    Rcpp::DataFrame survpatch_df(survproxy["patches"]);
    Rcpp::DataFrame obspatch_df(obsproxy["patches"]);
    Rcpp::DataFrame sizepatch_df(sizeproxy["patches"]);
    Rcpp::DataFrame sizebpatch_df(sizebproxy["patches"]);
    Rcpp::DataFrame sizecpatch_df(sizecproxy["patches"]);
    Rcpp::DataFrame repstpatch_df(repstproxy["patches"]);
    Rcpp::DataFrame fecpatch_df(fecproxy["patches"]);
    Rcpp::DataFrame jsurvpatch_df(jsurvproxy["patches"]);
    Rcpp::DataFrame jobspatch_df(jobsproxy["patches"]);
    Rcpp::DataFrame jsizepatch_df(jsizeproxy["patches"]);
    Rcpp::DataFrame jsizebpatch_df(jsizebproxy["patches"]);
    Rcpp::DataFrame jsizecpatch_df(jsizecproxy["patches"]);
    Rcpp::DataFrame jrepstpatch_df(jrepstproxy["patches"]);
    Rcpp::DataFrame jmatstpatch_df(jmatstproxy["patches"]);
    
    NumericVector survpatch = survpatch_df[0];
    NumericVector obspatch = obspatch_df[0];
    NumericVector sizepatch = sizepatch_df[0];
    NumericVector sizebpatch = sizebpatch_df[0];
    NumericVector sizecpatch = sizecpatch_df[0];
    NumericVector repstpatch = repstpatch_df[0];
    NumericVector fecpatch = fecpatch_df[0];
    NumericVector jsurvpatch = jsurvpatch_df[0];
    NumericVector jobspatch = jobspatch_df[0];
    NumericVector jsizepatch = jsizepatch_df[0];
    NumericVector jsizebpatch = jsizebpatch_df[0];
    NumericVector jsizecpatch = jsizecpatch_df[0];
    NumericVector jrepstpatch = jrepstpatch_df[0];
    NumericVector jmatstpatch = jmatstpatch_df[0];
    
    int matrows = survpatch.length();
    
    NumericMatrix patch_mat(matrows, 14);
    patch_mat(_, 0) = survpatch;
    patch_mat(_, 1) = obspatch;
    patch_mat(_, 2) = sizepatch;
    patch_mat(_, 3) = sizebpatch;
    patch_mat(_, 4) = sizecpatch;
    patch_mat(_, 5) = repstpatch;
    patch_mat(_, 6) = fecpatch;
    patch_mat(_, 7) = jsurvpatch;
    patch_mat(_, 8) = jobspatch;
    patch_mat(_, 9) = jsizepatch;
    patch_mat(_, 10) = jsizebpatch;
    patch_mat(_, 11) = jsizecpatch;
    patch_mat(_, 12) = jrepstpatch;
    patch_mat(_, 13) = jmatstpatch;
    
    final_mat = patch_mat;
  }
  
  return final_mat;
}

//' Creates a Summation of Most Terms Needed in Vital Rate Calculation
//' 
//' Function \code{.rimeotam()} provides the majority of the work in creating
//' the linear model sum to be used in vital rate estimation in the MPM. Works
//' specifically with functions \code{\link{.jerzeibalowski}()} and
//' \code{\link{.motherbalowski}()}.
//' 
//' @param maincoefs The coefficients portion of the vital rate model proxy.
//' @param fl1_i Reproductive status in time \emph{t}*-1.
//' @param fl2n_i Reproductive status in time \emph{t}.
//' @param sz1_i Primary size in time \emph{t}-1.
//' @param sz2o_i Primary size in time \emph{t}.
//' @param szb1_i Secondary size in time \emph{t}-1.
//' @param szb2o_i Secondary size in time \emph{t}.
//' @param szc1_i Tertiary size in time \emph{t}-1.
//' @param szc2o_i Tertiary size in time \emph{t}.
//' @param aage2_i Used age in time \emph{t}.
//' @param inda_1 Value of numeric individual covariate a in time \emph{t}-1.
//' @param inda_2 Value of numeric individual covariate a in time \emph{t}.
//' @param indb_1 Value of numeric individual covariate b in time \emph{t}-1.
//' @param indb_2 Value of numeric individual covariate b in time \emph{t}.
//' @param indc_1 Value of numeric individual covariate c in time \emph{t}-1.
//' @param indc_2 Value of numeric individual covariate c in time \emph{t}.
//' @param used_dens Density value used.
//' @param zi A logical value indicating whether model coefficients refer to the
//' zero inflation portion of a model.
//' 
//' @return A single numeric value giving the sum of the products of the linear
//' coefficients and the used status values.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.rimeotam)]]
double rimeotam(NumericVector maincoefs, double fl1_i, double fl2n_i, double sz1_i,
  double sz2o_i, double szb1_i, double szb2o_i, double szc1_i, double szc2o_i,
  double aage2_i, double inda_1, double inda_2, double indb_1, double indb_2,
  double indc_1, double indc_2, double used_dens, bool zi) {
  
  int add1 {0};
  int add2 {0};
  
  if (zi) {
    add1 = 46;
    add2 = 100;
  }
  
  double parti = maincoefs(0 + add1) + (maincoefs(1 + add1) * fl1_i) + (maincoefs(2 + add1) * fl2n_i) +
    (maincoefs(3 + add1) * sz1_i) + (maincoefs(4 + add1) * sz2o_i) + (maincoefs(5 + add1) * fl2n_i * fl1_i) + 
    (maincoefs(6 + add1) * sz2o_i * sz1_i) + (maincoefs(7 + add1) * sz1_i * fl1_i) +
    (maincoefs(8 + add1) * sz2o_i * fl2n_i) + (maincoefs(9 + add1) * sz2o_i * fl1_i) + 
    (maincoefs(10 + add1) * sz1_i * fl2n_i) + (maincoefs(11 + add1) * aage2_i) + 
    (maincoefs(12 + add1) * aage2_i * sz1_i) + (maincoefs(13 + add1) * aage2_i * sz2o_i) + 
    (maincoefs(14 + add1) * aage2_i * fl1_i) + (maincoefs(15 + add1) * aage2_i * fl2n_i) + 
    (maincoefs(16 + add1) * inda_2) + (maincoefs(17 + add1) * indb_2) + (maincoefs(18 + add1) * indc_2) + 
    (maincoefs(19 + add1) * inda_1) + (maincoefs(20 + add1) * indb_1) + (maincoefs(21 + add1) * indc_1) + 
    (maincoefs(22 + add1) * inda_2 * sz2o_i) + (maincoefs(23 + add1) * indb_2 * sz2o_i) + 
    (maincoefs(24 + add1) * indc_2 * sz2o_i) + (maincoefs(25 + add1) * inda_2 * fl2n_i) + 
    (maincoefs(26 + add1) * indb_2 * fl2n_i) + (maincoefs(27 + add1) * indc_2 * fl2n_i) + 
    (maincoefs(28 + add1) * inda_1 * sz1_i) + (maincoefs(29 + add1) * indb_1 * sz1_i) +
    (maincoefs(30 + add1) * indc_1 * sz1_i) + (maincoefs(31 + add1) * inda_1 * fl1_i) + 
    (maincoefs(32 + add1) * indb_1 * fl1_i) + (maincoefs(33 + add1) * indc_1 * fl1_i) + 
    (maincoefs(34 + add1) * inda_2 * indb_2) + (maincoefs(35 + add1) * inda_2 * indc_2) + 
    (maincoefs(36 + add1) * indb_2 * indc_2) + (maincoefs(37 + add1) * inda_1 * indb_1) + 
    (maincoefs(38 + add1) * inda_1 * indc_1) + (maincoefs(39 + add1) * indb_1 * indc_1) + 
    (maincoefs(40 + add1) * inda_2 * indb_1) + (maincoefs(41 + add1) * inda_1 * indb_2) +
    (maincoefs(42 + add1) * inda_2 * indc_1) + (maincoefs(43 + add1) * inda_1 * indc_2) + 
    (maincoefs(44 + add1) * indb_2 * indc_1) + (maincoefs(45 + add1) * indb_1 * indc_2);
    
  double partii = (maincoefs(100 + add2) * szb2o_i) + (maincoefs(101 + add2) * szb1_i) + 
    (maincoefs(102 + add2) * szc2o_i) + (maincoefs(103 + add2) * szc1_i) + 
    (maincoefs(104 + add2) * used_dens) + (maincoefs(105 + add2) * szb1_i * szb2o_i);
    
  double partiii = (maincoefs(106 + add2) * szc1_i * szc2o_i) + (maincoefs(107 + add2) * sz1_i * szb1_i) + 
    (maincoefs(108 + add2) * sz1_i * szc1_i) + (maincoefs(109 + add2) * szb1_i * szc1_i) + 
    (maincoefs(110 + add2) * sz2o_i * szb2o_i) + (maincoefs(111 + add2) * sz2o_i * szc2o_i) + 
    (maincoefs(112 + add2) * szb2o_i * szc2o_i) + (maincoefs(113 + add2) * sz1_i * szb2o_i) + 
    (maincoefs(114 + add2) * sz1_i * szc2o_i) + (maincoefs(115 + add2) * szb1_i * szc2o_i) + 
    (maincoefs(116 + add2) * sz2o_i * szb1_i) + (maincoefs(117 + add2) * sz2o_i * szc1_i) + 
    (maincoefs(118 + add2) * szb2o_i * szc1_i) + (maincoefs(119 + add2) * sz2o_i * used_dens) + 
    (maincoefs(120 + add2) * szb2o_i * used_dens) + (maincoefs(121 + add2) * szc2o_i * used_dens) + 
    (maincoefs(122 + add2) * sz1_i * used_dens) + (maincoefs(123 + add2) * szb1_i * used_dens) + 
    (maincoefs(124 + add2) * szc1_i * used_dens) + (maincoefs(125 + add2) * fl2n_i * used_dens) + 
    (maincoefs(126 + add2) * fl1_i * used_dens) + (maincoefs(127 + add2) * szb2o_i * fl2n_i) + 
    (maincoefs(128 + add2) * szc2o_i * fl2n_i) + 0 + (maincoefs(130 + add2) * szb1_i * fl1_i) + 
    (maincoefs(131 + add2) * szb2o_i * fl1_i) + (maincoefs(132 + add2) * szb1_i * fl2n_i) + 
    (maincoefs(133 + add2) * szc1_i * fl1_i) + (maincoefs(134 + add2) * szc2o_i * fl1_i) + 
    (maincoefs(135 + add2) * szc1_i * fl2n_i) + (maincoefs(136 + add2) * szb2o_i * aage2_i) + 
    (maincoefs(137 + add2) * szc2o_i * aage2_i) + (maincoefs(138 + add2) * used_dens * aage2_i) + 
    (maincoefs(139 + add2) * szb1_i * aage2_i) + (maincoefs(140 + add2) * szc1_i * aage2_i);
    
  double partiv = (maincoefs(141 + add2) * inda_2 * szb2o_i) + (maincoefs(142 + add2) * inda_2 * szc2o_i) + 
    (maincoefs(143 + add2) * inda_2 * used_dens) + (maincoefs(144 + add2) * inda_1 * szb1_i) + 
    (maincoefs(145 + add2) * inda_1 * szc1_i) + (maincoefs(146 + add2) * inda_1 * szb2o_i) + 
    (maincoefs(147 + add2) * inda_1 * szc2o_i) + (maincoefs(148 + add2) * inda_2 * szb1_i) + 
    (maincoefs(149 + add2) * inda_2 * szc1_i) + (maincoefs(150 + add2) * inda_1 * used_dens);
    
  double partv = (maincoefs(151 + add2) * indb_2 * szb2o_i) + (maincoefs(152 + add2) * indb_2 * szc2o_i) + 
    (maincoefs(153 + add2) * indb_2 * used_dens) + (maincoefs(154 + add2) * indb_1 * szb1_i) + 
    (maincoefs(155 + add2) * indb_1 * szc1_i) + (maincoefs(156 + add2) * indb_1 * szb2o_i) + 
    (maincoefs(157 + add2) * indb_1 * szc2o_i) + (maincoefs(158 + add2) * indb_2 * szb1_i) + 
    (maincoefs(159 + add2) * indb_2 * szc1_i) + (maincoefs(160 + add2) * indb_1 * used_dens);
    
  double partvi = (maincoefs(161 + add2) * indc_2 * szb2o_i) + (maincoefs(162 + add2) * indc_2 * szc2o_i) + 
    (maincoefs(163 + add2) * indc_2 * used_dens) + (maincoefs(164 + add2) * indc_1 * szb1_i) + 
    (maincoefs(165 + add2) * indc_1 * szc1_i) + (maincoefs(166 + add2) * indc_1 * szb2o_i) + 
    (maincoefs(167 + add2) * indc_1 * szc2o_i) + (maincoefs(168 + add2) * indc_2 * szb1_i) + 
    (maincoefs(169 + add2) * indc_2 * szc1_i) + (maincoefs(170 + add2) * indc_1 * used_dens);
    
  double partvii = (maincoefs(171 + add2) * inda_2 * sz1_i) + (maincoefs(172 + add2) * indb_2 * sz1_i) + 
    (maincoefs(173 + add2) * indc_2 * sz1_i) + (maincoefs(174 + add2) * inda_1 * sz2o_i) + 
    (maincoefs(175 + add2) * indb_1 * sz2o_i) + (maincoefs(176 + add2) * indc_1 * sz2o_i) + 
    (maincoefs(177 + add2) * inda_2 * fl1_i) + (maincoefs(178 + add2) * indb_2 * fl1_i) + 
    (maincoefs(179 + add2) * indc_2 * fl1_i) + (maincoefs(180 + add2) * inda_1 * fl2n_i) + 
    (maincoefs(181 + add2) * indb_1 * fl2n_i) + (maincoefs(182 + add2) * indc_1 * fl2n_i);
  
  double albatross = parti + partii + partiii + partiv + partv + partvi + partvii;
  
  return albatross;
}

//' Counts Numbers of Elements in Each Random Individual Covariate Portion of
//' Model
//' 
//' Function \code{.foi_counter()} counts the number of elements in each random
//' individual covariate and returns that as a vector.
//' 
//' @param modelproxy A list holding the contents of a model processed with
//' function \code{\link{.modelextract}()}
//' @param zi A logical value indicating whether to focus on the zero-inflation
//' parameters.
//' 
//' @return A 6 element vector holding the numbers of elements in each random
//' individual covariate in a model (either the cont portion or the zi portion).
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.foi_counter)]]
arma::ivec foi_counter(List modelproxy, bool zi) {
  
  arma::ivec return_vec(6, fill::zeros);
  
  if (!zi) {
    Rcpp::DataFrame modelinda2r_df(modelproxy["indcova2s"]);
    Rcpp::DataFrame modelinda1r_df(modelproxy["indcova1s"]);
    Rcpp::DataFrame modelindb2r_df(modelproxy["indcovb2s"]);
    Rcpp::DataFrame modelindb1r_df(modelproxy["indcovb1s"]);
    Rcpp::DataFrame modelindc2r_df(modelproxy["indcovc2s"]);
    Rcpp::DataFrame modelindc1r_df(modelproxy["indcovc1s"]);
    
    arma::vec modelinda2r = modelinda2r_df[0];
    arma::vec modelinda1r = modelinda1r_df[0];
    arma::vec modelindb2r = modelindb2r_df[0];
    arma::vec modelindb1r = modelindb1r_df[0];
    arma::vec modelindc2r = modelindc2r_df[0];
    arma::vec modelindc1r = modelindc1r_df[0];
    
    int v1_l = modelinda2r.n_elem;
    int v2_l = modelinda1r.n_elem;
    int v3_l = modelindb2r.n_elem;
    int v4_l = modelindb1r.n_elem;
    int v5_l = modelindc2r.n_elem;
    int v6_l = modelindc1r.n_elem;
    
    return_vec = {v1_l, v2_l, v3_l, v4_l, v5_l, v6_l};
  } else {
    Rcpp::DataFrame modelinda2r_df(modelproxy["zeroindcova2s"]);
    Rcpp::DataFrame modelinda1r_df(modelproxy["zeroindcova1s"]);
    Rcpp::DataFrame modelindb2r_df(modelproxy["zeroindcovb2s"]);
    Rcpp::DataFrame modelindb1r_df(modelproxy["zeroindcovb1s"]);
    Rcpp::DataFrame modelindc2r_df(modelproxy["zeroindcovc2s"]);
    Rcpp::DataFrame modelindc1r_df(modelproxy["zeroindcovc1s"]);
  
    arma::vec modelinda2r = modelinda2r_df[0];
    arma::vec modelinda1r = modelinda1r_df[0];
    arma::vec modelindb2r = modelindb2r_df[0];
    arma::vec modelindb1r = modelindb1r_df[0];
    arma::vec modelindc2r = modelindc2r_df[0];
    arma::vec modelindc1r = modelindc1r_df[0];
    
    int v1_l = modelinda2r.n_elem;
    int v2_l = modelinda1r.n_elem;
    int v3_l = modelindb2r.n_elem;
    int v4_l = modelindb1r.n_elem;
    int v5_l = modelindc2r.n_elem;
    int v6_l = modelindc1r.n_elem;
    
    return_vec = {v1_l, v2_l, v3_l, v4_l, v5_l, v6_l};
  }
  
  return return_vec;
}

//' Create Vector of Random Individual Covariate Terms
//' 
//' Function \code{.flightoficarus()} creates vectors of random covariate
//' terms.
//' 
//' @param modelproxy A model proxy list extracted with function
//' \code{\link{.modelextract}()}.
//' 
//' @return A vector of numeric values for random categorical terms. The order
//' is: 1) cov a time 2, 2) cov a time 1, 3) cov b time 2, 4) cov b time 1,
//' 5) cov c time 2, and 6) cov c time 1. Rows may vary, but must be the same
//' length for each model.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.flightoficarus)]]
NumericVector flightoficarus(List modelproxy) {
  Rcpp::DataFrame modelinda2r_df(modelproxy["indcova2s"]);
  Rcpp::DataFrame modelinda1r_df(modelproxy["indcova1s"]);
  Rcpp::DataFrame modelindb2r_df(modelproxy["indcovb2s"]);
  Rcpp::DataFrame modelindb1r_df(modelproxy["indcovb1s"]);
  Rcpp::DataFrame modelindc2r_df(modelproxy["indcovc2s"]);
  Rcpp::DataFrame modelindc1r_df(modelproxy["indcovc1s"]);
  
  NumericVector modelinda2r = modelinda2r_df[0];
  NumericVector modelinda1r = modelinda1r_df[0];
  NumericVector modelindb2r = modelindb2r_df[0];
  NumericVector modelindb1r = modelindb1r_df[0];
  NumericVector modelindc2r = modelindc2r_df[0];
  NumericVector modelindc1r = modelindc1r_df[0];
  
  int v1_l = modelinda2r.length();
  int v2_l = modelinda1r.length();
  int v3_l = modelindb2r.length();
  int v4_l = modelindb1r.length();
  int v5_l = modelindc2r.length();
  int v6_l = modelindc1r.length();
  int vec_length = v1_l + v2_l + v3_l + v4_l + v5_l + v6_l;
  
  NumericVector final_vec(vec_length);
  int all_counter {0};
  
  for (int i = 0; i < v1_l; i++) {
    final_vec(all_counter) = modelinda2r(i);
    all_counter++;
  }
  for (int i = 0; i < v2_l; i++) {
    final_vec(all_counter) = modelinda1r(i);
    all_counter++;
  }
  for (int i = 0; i < v3_l; i++) {
    final_vec(all_counter) = modelindb2r(i);
    all_counter++;
  }
  for (int i = 0; i < v4_l; i++) {
    final_vec(all_counter) = modelindb1r(i);
    all_counter++;
  }
  for (int i = 0; i < v5_l; i++) {
    final_vec(all_counter) = modelindc2r(i);
    all_counter++;
  }
  for (int i = 0; i < v6_l; i++) {
    final_vec(all_counter) = modelindc1r(i);
    all_counter++;
  }
  
  return final_vec;
}

//' Create Concatenated Vector of Random Individual Covariate Term Names
//' 
//' Function \code{.bootson()} creates a concatenated string vector holding all
//' covariate term names.
//' 
//' @param modelproxy A model proxy list extracted with function
//' \code{\link{.modelextract}()}.
//' 
//' @return A vector holding all covariate name terms. The order is: 1) cov a
//' time 2, 2) cov a time 1, 3) cov b time 2, 4) cov b time 1, 5) cov c time 2,
//' and 6) cov c time 1. Note that the element order is the same as in function
//' \code{\link{.flightoficarus}()}.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.bootson)]]
StringVector bootson(List modelproxy) {
  Rcpp::DataFrame modelinda2r_df(modelproxy["indcova2s"]);
  Rcpp::DataFrame modelinda1r_df(modelproxy["indcova1s"]);
  Rcpp::DataFrame modelindb2r_df(modelproxy["indcovb2s"]);
  Rcpp::DataFrame modelindb1r_df(modelproxy["indcovb1s"]);
  Rcpp::DataFrame modelindc2r_df(modelproxy["indcovc2s"]);
  Rcpp::DataFrame modelindc1r_df(modelproxy["indcovc1s"]);
  
  StringVector modelinda2r_rownames = modelinda2r_df.attr("row.names");
  StringVector modelinda1r_rownames = modelinda1r_df.attr("row.names");
  StringVector modelindb2r_rownames = modelindb2r_df.attr("row.names");
  StringVector modelindb1r_rownames = modelindb1r_df.attr("row.names");
  StringVector modelindc2r_rownames = modelindc2r_df.attr("row.names");
  StringVector modelindc1r_rownames = modelindc1r_df.attr("row.names");
  
  int v1_l = modelinda2r_rownames.length();
  int v2_l = modelinda1r_rownames.length();
  int v3_l = modelindb2r_rownames.length();
  int v4_l = modelindb1r_rownames.length();
  int v5_l = modelindc2r_rownames.length();
  int v6_l = modelindc1r_rownames.length();
  int vec_length = v1_l + v2_l + v3_l + v4_l + v5_l + v6_l;

  StringVector final_vec(vec_length);
  int all_counter {0};
  
  for (int i = 0; i < v1_l; i++) {
    final_vec(all_counter) = modelinda2r_rownames(i);
    all_counter++;
  }
  for (int i = 0; i < v2_l; i++) {
    final_vec(all_counter) = modelinda1r_rownames(i);
    all_counter++;
  }
  for (int i = 0; i < v3_l; i++) {
    final_vec(all_counter) = modelindb2r_rownames(i);
    all_counter++;
  }
  for (int i = 0; i < v4_l; i++) {
    final_vec(all_counter) = modelindb1r_rownames(i);
    all_counter++;
  }
  for (int i = 0; i < v5_l; i++) {
    final_vec(all_counter) = modelindc2r_rownames(i);
    all_counter++;
  }
  for (int i = 0; i < v6_l; i++) {
    final_vec(all_counter) = modelindc1r_rownames(i);
    all_counter++;
  }
  
  return final_vec;
}

//' Create Vector of Random Individual Covariate Terms for Zero-Inflated Models
//' 
//' Function \code{.zero_flightoficarus()} creates vectors of random covariate
//' terms from the binomial portion of a zero-inflated model.
//' 
//' @param modelproxy A model proxy list extracted with function
//' \code{\link{.modelextract}()}.
//' 
//' @return A vector of numeric values for random categorical terms. The order
//' is: 1) cov a time 2, 2) cov a time 1, 3) cov b time 2, 4) cov b time 1,
//' 5) cov c time 2, and 6) cov c time 1. Rows may vary, but must be the same
//' length for each model.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.zero_flightoficarus)]]
NumericVector zero_flightoficarus(List modelproxy) {
  Rcpp::DataFrame modelinda2r_df(modelproxy["zeroindcova2s"]);
  Rcpp::DataFrame modelinda1r_df(modelproxy["zeroindcova1s"]);
  Rcpp::DataFrame modelindb2r_df(modelproxy["zeroindcovb2s"]);
  Rcpp::DataFrame modelindb1r_df(modelproxy["zeroindcovb1s"]);
  Rcpp::DataFrame modelindc2r_df(modelproxy["zeroindcovc2s"]);
  Rcpp::DataFrame modelindc1r_df(modelproxy["zeroindcovc1s"]);
  
  NumericVector modelinda2r = modelinda2r_df[0];
  NumericVector modelinda1r = modelinda1r_df[0];
  NumericVector modelindb2r = modelindb2r_df[0];
  NumericVector modelindb1r = modelindb1r_df[0];
  NumericVector modelindc2r = modelindc2r_df[0];
  NumericVector modelindc1r = modelindc1r_df[0];
  
  int v1_l = modelinda2r.length();
  int v2_l = modelinda1r.length();
  int v3_l = modelindb2r.length();
  int v4_l = modelindb1r.length();
  int v5_l = modelindc2r.length();
  int v6_l = modelindc1r.length();
  int vec_length = v1_l + v2_l + v3_l + v4_l + v5_l + v6_l;
  
  NumericVector final_vec(vec_length);
  int all_counter {0};
  
  for (int i = 0; i < v1_l; i++) {
    final_vec(all_counter) = modelinda2r(i);
    all_counter++;
  }
  for (int i = 0; i < v2_l; i++) {
    final_vec(all_counter) = modelinda1r(i);
    all_counter++;
  }
  for (int i = 0; i < v3_l; i++) {
    final_vec(all_counter) = modelindb2r(i);
    all_counter++;
  }
  for (int i = 0; i < v4_l; i++) {
    final_vec(all_counter) = modelindb1r(i);
    all_counter++;
  }
  for (int i = 0; i < v5_l; i++) {
    final_vec(all_counter) = modelindc2r(i);
    all_counter++;
  }
  for (int i = 0; i < v6_l; i++) {
    final_vec(all_counter) = modelindc1r(i);
    all_counter++;
  }
  
  return final_vec;
}

//' Create Concatenated Vector of Random Individual Covariate Term Names from
//' a Zero-Inflated Model
//' 
//' Function \code{.zero_bootson()} creates a concatenated string vector holding
//' all covariate term names from the binomial portion of a zero-inflated model.
//' 
//' @param modelproxy A model proxy list extracted with function
//' \code{\link{.modelextract}()}.
//' 
//' @return A vector holding all covariate name terms. The order is: 1) cov a
//' time 2, 2) cov a time 1, 3) cov b time 2, 4) cov b time 1, 5) cov c time 2,
//' and 6) cov c time 1. Note that the element order is the same as in function
//' \code{\link{.zero_flightoficarus}()}.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.zero_bootson)]]
StringVector zero_bootson(List modelproxy) {
  Rcpp::DataFrame modelinda2r_df(modelproxy["zeroindcova2s"]);
  Rcpp::DataFrame modelinda1r_df(modelproxy["zeroindcova1s"]);
  Rcpp::DataFrame modelindb2r_df(modelproxy["zeroindcovb2s"]);
  Rcpp::DataFrame modelindb1r_df(modelproxy["zeroindcovb1s"]);
  Rcpp::DataFrame modelindc2r_df(modelproxy["zeroindcovc2s"]);
  Rcpp::DataFrame modelindc1r_df(modelproxy["zeroindcovc1s"]);
  
  StringVector modelinda2r_rownames = modelinda2r_df.attr("row.names");
  StringVector modelinda1r_rownames = modelinda1r_df.attr("row.names");
  StringVector modelindb2r_rownames = modelindb2r_df.attr("row.names");
  StringVector modelindb1r_rownames = modelindb1r_df.attr("row.names");
  StringVector modelindc2r_rownames = modelindc2r_df.attr("row.names");
  StringVector modelindc1r_rownames = modelindc1r_df.attr("row.names");
  
  int v1_l = modelinda2r_rownames.length();
  int v2_l = modelinda1r_rownames.length();
  int v3_l = modelindb2r_rownames.length();
  int v4_l = modelindb1r_rownames.length();
  int v5_l = modelindc2r_rownames.length();
  int v6_l = modelindc1r_rownames.length();
  int vec_length = v1_l + v2_l + v3_l + v4_l + v5_l + v6_l;

  StringVector final_vec(vec_length);
  int all_counter {0};
  
  for (int i = 0; i < v1_l; i++) {
    final_vec(all_counter) = modelinda2r_rownames(i);
    all_counter++;
  }
  for (int i = 0; i < v2_l; i++) {
    final_vec(all_counter) = modelinda1r_rownames(i);
    all_counter++;
  }
  for (int i = 0; i < v3_l; i++) {
    final_vec(all_counter) = modelindb2r_rownames(i);
    all_counter++;
  }
  for (int i = 0; i < v4_l; i++) {
    final_vec(all_counter) = modelindb1r_rownames(i);
    all_counter++;
  }
  for (int i = 0; i < v5_l; i++) {
    final_vec(all_counter) = modelindc2r_rownames(i);
    all_counter++;
  }
  for (int i = 0; i < v6_l; i++) {
    final_vec(all_counter) = modelindc1r_rownames(i);
    all_counter++;
  }
  
  return final_vec;
}

//' Create Index of Element Numbers for Random Individual Covariate Terms
//' 
//' Function \code{.foi_index()} creates a matrix indexing the end points of
//' each random individual covariate in the utilized vectors.
//' 
//' @param surv_proxy Adult survival model proxy.
//' @param obs_proxy Adult observation status model proxy.
//' @param size_proxy Adult primary size model proxy.
//' @param sizeb_proxy Adult secondary size model proxy.
//' @param sizec_proxy Adult tertiary size model proxy.
//' @param repst_proxy Adult reproductive status model proxy.
//' @param fec_proxy Adult fecundity model proxy.
//' @param jsurv_proxy Juvenile survival model proxy.
//' @param jobs_proxy Juvenile observation status model proxy.
//' @param jsize_proxy Juvenile primary size model proxy.
//' @param jsizeb_proxy Juvenile secondary size model proxy.
//' @param jsizec_proxy Juvenile tertiary size model proxy.
//' @param jrepst_proxy Juvenile reproductive status model proxy.
//' @param jmatst_proxy Juvenile maturity status model proxy.
//' 
//' @return An integer matrix with 6 rows and 20 columns. The columns contain
//' the number of elements in each random individual covariate term, with the
//' row order being: 1) cov a t2, 2) cov a t1, 3) cov b t2, 4) cov b t1,
//' 5) cov c t2, and 6) cov c t1.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.foi_index)]]
arma::imat foi_index(List surv_proxy, List obs_proxy, List size_proxy, 
  List sizeb_proxy, List sizec_proxy, List repst_proxy, List fec_proxy,
  List jsurv_proxy, List jobs_proxy, List jsize_proxy, List jsizeb_proxy,
  List jsizec_proxy, List jrepst_proxy, List jmatst_proxy) {
  
  arma::ivec surv_fc = foi_counter(surv_proxy, false);
  arma::ivec obs_fc = foi_counter(obs_proxy, false);
  arma::ivec size_fc = foi_counter(size_proxy, false);
  arma::ivec sizeb_fc = foi_counter(sizeb_proxy, false);
  arma::ivec sizec_fc = foi_counter(sizec_proxy, false);
  arma::ivec repst_fc = foi_counter(repst_proxy, false);
  arma::ivec fec_fc = foi_counter(fec_proxy, false);
  arma::ivec jsurv_fc = foi_counter(jsurv_proxy, false);
  arma::ivec jobs_fc = foi_counter(jobs_proxy, false);
  arma::ivec jsize_fc = foi_counter(jsize_proxy, false);
  arma::ivec jsizeb_fc = foi_counter(jsizeb_proxy, false);
  arma::ivec jsizec_fc = foi_counter(jsizec_proxy, false);
  arma::ivec jrepst_fc = foi_counter(jrepst_proxy, false);
  arma::ivec jmatst_fc = foi_counter(jmatst_proxy, false);
  arma::ivec size_fc_zi = foi_counter(size_proxy, true);
  arma::ivec sizeb_fc_zi = foi_counter(sizeb_proxy, true);
  arma::ivec sizec_fc_zi = foi_counter(sizec_proxy, true);
  arma::ivec fec_fc_zi = foi_counter(fec_proxy, true);
  arma::ivec jsize_fc_zi = foi_counter(jsize_proxy, true);
  arma::ivec jsizeb_fc_zi = foi_counter(jsizeb_proxy, true);
  arma::ivec jsizec_fc_zi = foi_counter(jsizec_proxy, true);
  
  arma::imat final_mat(6, 21, fill::zeros);
  
  for (int i = 0; i < 6; i++) {
    final_mat(i, 0) = surv_fc(i);
    final_mat(i, 1) = obs_fc(i);
    final_mat(i, 2) = size_fc(i);
    final_mat(i, 3) = sizeb_fc(i);
    final_mat(i, 4) = sizec_fc(i);
    final_mat(i, 5) = repst_fc(i);
    final_mat(i, 6) = fec_fc(i);
    final_mat(i, 7) = jsurv_fc(i);
    final_mat(i, 8) = jobs_fc(i);
    final_mat(i, 9) = jsize_fc(i);
    final_mat(i, 10) = jsizeb_fc(i);
    final_mat(i, 11) = jsizec_fc(i);
    final_mat(i, 12) = jrepst_fc(i);
    final_mat(i, 13) = size_fc_zi(i);
    final_mat(i, 14) = sizeb_fc_zi(i);
    final_mat(i, 15) = sizec_fc_zi(i);
    final_mat(i, 16) = fec_fc_zi(i);
    final_mat(i, 17) = jsize_fc_zi(i);
    final_mat(i, 18) = jsizeb_fc_zi(i);
    final_mat(i, 19) = jsizec_fc_zi(i);
    final_mat(i, 20) = jrepst_fc(i);
  }
  
  return final_mat;
}

//' Estimate Value for Vital Rate Based on Inputs
//' 
//' Function \code{.preouterator()} calculates the value of the vital rate
//' called for by the function.
//' 
//' @param modelproxy A model_proxy object derived from function
//' \code{\link(.modelextract)()}.
//' @param maincoefs The coefficients portion of the vital rate model proxy.
//' @param randindex An integer matrix indexing all random covariates for all
//' vital rates.
//' @param dev_terms A numeric vector containing the deviations to the linear
//' models input by the user. The order is: survival, observation status, size,
//' size_b, size_c, reproductive status, fecundity, juvenile survival, juvenile
//' observation status, juvenile size, juvenile size_b, juvenile size_c,
//' and juvenile reproductive status.
//' @param svsigmas A vector of sigma and summedvar terms from vital rate
//' models, in the order of: summedvars, sigma, summedvarsb, sigmab,
//' summedvarsc, sigmac, jsummedvars, jsigma, jsummedvarsb, jsigmab,
//' jsummedvarsc, and jsigmac. Summedvar terms are summed variance-covariance
//' terms in Poisson and negative binomial size distributions, and sigma terms
//' are standard deviations in the Gaussian size distribution.
//' @param vitalyear A matrix with year coefficients for all vital rates.
//' @param vitalpatch A matrix with patch coefficients for all vital rates.
//' @param chosen_r2inda A string identifying random covariate a in time t.
//' @param chosen_r1inda A string identifying random covariate a in time t-1.
//' @param chosen_r2indb A string identifying random covariate b in time t.
//' @param chosen_r1indb A string identifying random covariate b in time t-1.
//' @param chosen_r2indc A string identifying random covariate c in time t.
//' @param chosen_r1indc A string identifying random covariate c in time t-1.
//' @param status_terms A NumericVector containing, in order: fl1_i, fl2n_i,
//' sz1_i, sz2o_i, szb1_i, szb2o_i, szc1_i, szc2o_i, aage2_i, inda_1, inda_2,
//' indb_1, indb_2, indc_1, indc_2, used_dens, sz3_i, szb3_i, szc3_i,
//' binwidth3_i, binbwidth3_i, and bincwidth3_i.
//' @param modelgroups2 A vector of group slope coefficients for time t.
//' @param modelgroups1 A vector of group slope coefficients for time t-1.
//' @param modelgroups2zi A vector of zero-inflation model group slope
//' coefficients for time t.
//' @param modelgroups1zi A vector of zero-inflation model group slope
//' coefficients for time t-1.
//' @param modelyearzi A vector of zero-inflation model time slope coefficients.
//' @param modelpatchzi A vector of zero-inflation model patch slope coefficients.
//' @param modelind A vector of individual covariate slope coefficients.
//' @param modelind_rownames A string vector with the names of the individual
//' covariate coefficients.
//' @param modelindzi A vector of individual covariate slope coefficients.
//' @param modelind_rownames_zi A string vector with the names of the individual
//' covariate coefficients.
//' @param zi A logical value indicating whether model coefficients refer to the
//' zero inflation portion of a model.
//' @param grp2o_i Stage group number in time t.
//' @param grp1_i Stage group number in time t-1.
//' @param patchnumber An integer index for pop-patch.
//' @param yearnumber An integer index for monitoring occasion in time t.
//' @param vitaldist A parameter specifying the distribution of the vital rate.
//' Current options are: Poisson (0), negative binomial (1), Gaussian (2),
//' Gamma (3), and binomial (4).
//' @param vitalrate An integer specifying the vital rate. 1 = surv, 2 = obs,
//' 3 = size, 4 = sizeb, 5 = sizec, 6 = repst, 7 = fec, 8 = jsurv, 9 = jobs,
//' 10 = jsize, 11 = jsizeb, 12 = jsizec, 13 = jrepst, 14 = jmatst.
//' @param exp_tol A numeric value indicating the maximum limit for the
//' \code{exp()} function to be used in vital rate calculations. Defaults to
//' \code{700.0}.
//' @param theta_tol A numeric value indicating a maximum value for theta in
//' negative binomial probability density estimation. Defaults to
//' \code{100000000.0}.
//' @param ipm_cdf A logical value indicating whether to use the cumulative
//' density function to estimate size transitions in continuous distributions
//' (\code{true}), or the midpoint method (\code{false}).
//' @param matrixformat An integer representing the style of matrix to develop.
//' Options include Ehrlen-format hMPM (1), deVries-format hMPM (2), ahMPM (3),
//' and age-by-stage MPM (4).
//' @param fecmod A scalar multiplier for fecundity.
//' @param repentry_i Rep entry value for time t+1.
//' @param negfec A logical value denoting whether to change negative estimated
//' fecundity to 0.
//' @param stage2n_i Numeric index of stage in time t.
//' @param nostages The total number of stages in the stageframe.
//' @param modeltrunc An integer coding for zero-truncation status.
//' @param modelsigma A double numeric holding the standard deviation of the
//' parameter distribution.
//' 
//' @return A class double numeric value for the vital rate being estimated.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.preouterator)]]
double preouterator(List modelproxy, NumericVector maincoefs, arma::imat randindex,
  NumericVector dev_terms, NumericVector svsigmas, NumericMatrix vitalyear,
  NumericMatrix vitalpatch, String chosen_r2inda, String chosen_r1inda,
  String chosen_r2indb, String chosen_r1indb, String chosen_r2indc,
  String chosen_r1indc, NumericVector status_terms, NumericVector modelgroups2,
  NumericVector modelgroups1, NumericVector modelgroups2zi,
  NumericVector modelgroups1zi, NumericVector modelyearzi,
  NumericVector modelpatchzi, NumericVector modelind,
  StringVector modelind_rownames, NumericVector modelindzi,
  StringVector modelind_rownames_zi, bool zi, double grp2o_i,
  double grp1_i, int patchnumber, int yearnumber, int vitaldist, int vitalrate,
  double exp_tol, double theta_tol, bool ipm_cdf, int matrixformat,
  double fecmod, double repentry_i, bool negfec, double stage2n_i,
  int nostages, int modeltrunc, double modelsigma) {
  
  double preout {0.0};
  double all_out {0.0};
  
  int placeholder = vitalrate - 1;
  int placeholder_zi = placeholder + 11;
  int vitaltype {0}; // Binomial vital rates
  if (vitalrate == 3 || vitalrate == 4 || vitalrate == 5) {
    vitaltype = 1; // Size
  } else if (vitalrate == 10 || vitalrate == 11 || vitalrate == 12) {
    vitaltype = 1; // Juv size
  } else if (vitalrate == 7) {
    vitaltype = 2; // Fecundity
  }
  
  // This section occurs in all vital rates
  double mainsum = rimeotam(maincoefs, status_terms(0), status_terms(1),
    status_terms(2), status_terms(3), status_terms(4), status_terms(5),
    status_terms(6), status_terms(7), status_terms(8), status_terms(9),
    status_terms(10), status_terms(11), status_terms(12), status_terms(13),
    status_terms(14), status_terms(15), zi);
  
  bool zi_processing = false;
  
  if (vitaltype == 1) {
    if (vitalrate == 3 || vitalrate == 10) {
      if (zi && status_terms(16) == 0.0) zi_processing = true;
    } else if (vitalrate == 4 || vitalrate == 11) {
      if (zi && status_terms(17) == 0.0) zi_processing = true;
    } else if (vitalrate == 5 || vitalrate == 12) {
      if (zi && status_terms(18) == 0.0) zi_processing = true;
    } 
  } else if (vitaltype == 2) {
    if (zi && status_terms(16) == 0 && status_terms(17) == 0 &&
      status_terms(18) == 0 && vitaldist < 2) zi_processing = true;  
  }
  
  if (!zi_processing) {
    // Random covariate processing
    double chosen_randcova2 {0.0};
    if (chosen_r2inda != "none") {
      for (int indcount = 0; indcount < randindex(0, placeholder); indcount++) {
        if (chosen_r2inda == modelind_rownames(indcount)) {
          chosen_randcova2 = modelind(indcount);
        }
      }
    }
    double chosen_randcova1 {0.0};
    if (chosen_r1inda != "none") {
      int delectable_sum = randindex(0, placeholder);
      for (int indcount = 0; indcount < randindex(1, placeholder); indcount++) {
        if (chosen_r1inda == modelind_rownames(indcount + delectable_sum)) {
          chosen_randcova1 = modelind(indcount + delectable_sum);
        }
      }
    }
    double chosen_randcovb2 {0.0};
    if (chosen_r2indb != "none") {
      int delectable_sum = randindex(0, placeholder) + randindex(1, placeholder);
      for (int indcount = 0; indcount < randindex(2, placeholder); indcount++) {
        if (chosen_r2indb == modelind_rownames(indcount + delectable_sum)) {
          chosen_randcovb2 = modelind(indcount + delectable_sum);
        }
      }
    }
    double chosen_randcovb1 {0.0};
    if (chosen_r1indb != "none") {
      int delectable_sum = randindex(0, placeholder) + randindex(1, placeholder) +
        randindex(2, placeholder);
      for (int indcount = 0; indcount < randindex(3, placeholder); indcount++) {
        if (chosen_r1indb == modelind_rownames(indcount + delectable_sum)) {
          chosen_randcovb1 = modelind(indcount + delectable_sum);
        }
      }
    }
    double chosen_randcovc2 {0.0};
    if (chosen_r2indc != "none") {
      int delectable_sum = randindex(0, placeholder) + randindex(1, placeholder) +
        randindex(2, placeholder) + randindex(3, placeholder);
      for (int indcount = 0; indcount < randindex(4, placeholder); indcount++) {
        if (chosen_r2indc == modelind_rownames(indcount + delectable_sum)) {
          chosen_randcovc2 = modelind(indcount + delectable_sum);
        }
      }
    }
    double chosen_randcovc1 {0.0};
    if (chosen_r1indc != "none") {
      int delectable_sum = randindex(0, placeholder) + randindex(1, placeholder) +
        randindex(2, placeholder) + randindex(3, placeholder) + randindex(4, placeholder);
      for (int indcount = 0; indcount < randindex(5, placeholder); indcount++) {
        if (chosen_r1indc == modelind_rownames(indcount + delectable_sum)) {
          chosen_randcovc1 = modelind(indcount + delectable_sum);
        }
      }
    }
    
    preout = (mainsum + chosen_randcova2 + chosen_randcova1 +
      chosen_randcovb2 + chosen_randcovb1 + chosen_randcovc2 +
      chosen_randcovc1 + modelgroups2(grp2o_i) + modelgroups1(grp1_i) + 
      vitalpatch(patchnumber, placeholder) + vitalyear(yearnumber, placeholder) +
      dev_terms(placeholder));
      
    if (preout > exp_tol && vitaldist < 2) preout = exp_tol;
  } else {
    // Only for size and fec
    double chosen_randcova2zi {0.0};
    if (chosen_r2inda != "none") {
      for (int indcount = 0; indcount < randindex(0, placeholder_zi); indcount++) {
        if (chosen_r2inda == modelind_rownames_zi(indcount)) {
          chosen_randcova2zi = modelindzi(indcount);
        }
      }
    }
    double chosen_randcova1zi {0.0};
    if (chosen_r1inda != "none") {
      int delectable_sum = randindex(0, placeholder_zi);
      for (int indcount = 0; indcount < randindex(1, placeholder_zi); indcount++) {
        if (chosen_r1inda == modelind_rownames_zi(indcount + delectable_sum)) {
          chosen_randcova1zi = modelindzi(indcount + delectable_sum);
        }
      }
    }
    double chosen_randcovb2zi {0.0};
    if (chosen_r2indb != "none") {
      int delectable_sum = randindex(0, placeholder_zi) + randindex(1, placeholder_zi);
      for (int indcount = 0; indcount < randindex(2, placeholder_zi); indcount++) {
        if (chosen_r2indb == modelind_rownames_zi(indcount + delectable_sum)) {
          chosen_randcovb2zi = modelindzi(indcount + delectable_sum);
        }
      }
    }
    double chosen_randcovb1zi {0.0};
    if (chosen_r1indb != "none") {
      int delectable_sum = randindex(0, placeholder_zi) + randindex(1, placeholder_zi) +
        randindex(2, placeholder_zi);
      for (int indcount = 0; indcount < randindex(3, placeholder_zi); indcount++) {
        if (chosen_r1indb == modelind_rownames_zi(indcount + delectable_sum)) {
          chosen_randcovb1zi = modelindzi(indcount + delectable_sum);
        }
      }
    }
    double chosen_randcovc2zi {0.0};
    if (chosen_r2indc != "none") {
      int delectable_sum = randindex(0, placeholder_zi) + randindex(1, placeholder_zi) +
        randindex(2, placeholder_zi) + randindex(3, placeholder_zi);
      for (int indcount = 0; indcount < randindex(4, placeholder_zi); indcount++) {
        if (chosen_r2indc == modelind_rownames_zi(indcount + delectable_sum)) {
          chosen_randcovc2zi = modelindzi(indcount + delectable_sum);
        }
      }
    }
    double chosen_randcovc1zi {0.0};
    if (chosen_r1indc != "none") {
      int delectable_sum = randindex(0, placeholder_zi) + randindex(1, placeholder_zi) +
        randindex(2, placeholder_zi) + randindex(3, placeholder_zi) + randindex(4, placeholder_zi);
      for (int indcount = 0; indcount < randindex(5, placeholder_zi); indcount++) {
        if (chosen_r1indc == modelind_rownames_zi(indcount + delectable_sum)) {
          chosen_randcovc1zi = modelindzi(indcount + delectable_sum);
        }
      }
    }
    
    preout = (mainsum + chosen_randcova2zi + chosen_randcova1zi +
      chosen_randcovb2zi + chosen_randcovb1zi + chosen_randcovc2zi +
      chosen_randcovc1zi + modelgroups2zi(grp2o_i) + modelgroups1zi(grp1_i) + 
      modelpatchzi(patchnumber) + modelyearzi(yearnumber) +
      dev_terms(placeholder));
  }
  
  if (vitaltype == 0) {
    if (preout > exp_tol) preout = exp_tol;
      
    double pre_exp = exp(preout);
    all_out = pre_exp / (1.0 + pre_exp);
    
    // Rcout << "Binomial: preout: " << preout << " pre_exp: " << pre_exp <<
    //   " all_out: " << all_out << "\n";
  } else if (vitaltype == 1) {
    
    double Used_size3 = status_terms(16);
    double Used_binwidth3 = status_terms(19);
    double Used_svsigmas0 = svsigmas(0);
    double Used_svsigmas1 = svsigmas(1);
    
    if (vitalrate == 4) {
      Used_size3 = status_terms(17);
      Used_binwidth3 = status_terms(20);
      
      Used_svsigmas0 = svsigmas(2);
      Used_svsigmas1 = svsigmas(3);
    } else if (vitalrate == 5) {
      Used_size3 = status_terms(18);
      Used_binwidth3 = status_terms(21);
      
      Used_svsigmas0 = svsigmas(4);
      Used_svsigmas1 = svsigmas(5);
    } else if (vitalrate == 10) {
      Used_svsigmas0 = svsigmas(6);
      Used_svsigmas1 = svsigmas(7);
    } else if (vitalrate == 11) {
      Used_size3 = status_terms(17);
      Used_binwidth3 = status_terms(20);
      
      Used_svsigmas0 = svsigmas(8);
      Used_svsigmas1 = svsigmas(9);
    } else if (vitalrate == 12) {
      Used_size3 = status_terms(18);
      Used_binwidth3 = status_terms(21);
      
      Used_svsigmas0 = svsigmas(10);
      Used_svsigmas1 = svsigmas(11);
    }
    
    if (zi_processing) {
      preout = preout + (Used_svsigmas0 / 2.0);
      
      if (preout > exp_tol) preout = exp_tol;
      
      double pre_exp = exp(preout);
      all_out = pre_exp / (1.0 + pre_exp);
      
      // Rcout << "ZI Binomial: preout: " << preout << " pre_exp: " << pre_exp <<
      //   " all_out: " << all_out << "\n";
      
    } else {
      if (vitaldist == 0) {
        // Poisson distribution
        
        preout = preout + (Used_svsigmas0 / 2.0);
        
        if (preout > exp_tol) preout = exp_tol;
        double lambda = exp(preout);
        
        double upper_boundary = (Used_size3 + (Used_binwidth3 / 2));
        double upper_boundary_int = floor(upper_boundary);
        
        double lower_boundary = (Used_size3 - (Used_binwidth3 / 2));
        double lower_boundary_int = floor(lower_boundary);
        
        if (ipm_cdf) {
          if (lower_boundary_int < 0.0) lower_boundary_int = 0.0;
          
          double sizefac {1.0};
          if (upper_boundary_int > 0.0) {
            sizefac = upper_boundary_int * tgamma(upper_boundary_int);
          }
          double main_out = boost::math::tgamma((upper_boundary_int + 1), lambda) / sizefac;
          
          if (upper_boundary_int > lower_boundary_int) {
            double sizefac_low {1.0};
            if (lower_boundary_int > 0.0) {
              sizefac_low = lower_boundary_int * tgamma(lower_boundary_int);
            }
            all_out = main_out - boost::math::tgamma((lower_boundary_int + 1), lambda) / sizefac_low;
          } else {
            all_out = main_out;
          }
          
          // Rcout << "Poisson cdf: upper_boundary_int: " << upper_boundary_int << 
          //   " lower_boundary_int: " << lower_boundary_int << " lambda: " << lambda << 
          //   " all_out: " << all_out << "\n";
          
        } else {
          int y = static_cast<int>(upper_boundary_int);
          int y0 = static_cast<int>(lower_boundary_int);
          if (y0 < -1) y0 = -1;
          
          double current_prob {0.0};
          
          for (int summed_size = (y0 + 1); summed_size <= y; summed_size++) {
            double sizefac {1.0};
            if (Used_size3 > 0.0) {
              sizefac = Used_size3 * tgamma(Used_size3);
            }
            
            double den_corr {1.0};
            if (modeltrunc == 1) den_corr = (1.0 - (exp(-1 * lambda)));
            
            current_prob += ((pow(lambda, Used_size3) * exp(-1.0 * lambda)) / sizefac) / den_corr;
          }
          all_out = current_prob;
          
          // Rcout << "Poisson mid: upper_boundary_int: " << upper_boundary_int <<
          //   " lower_boundary_int: " << lower_boundary_int << " lambda: " << lambda << 
          //   " current_prob: " << current_prob << "\n";
        }
        
      } else if (vitaldist == 1) {
        // Negative binomial
        
        double mu = exp(preout);
        
        double theta = modelproxy["sigma"];
        if (NumericVector::is_na(theta)) theta = 1.0;
        if (theta > theta_tol) theta = theta_tol;
        double alpha = 1.0 / theta;
        
        double upper_boundary = (Used_size3 + (Used_binwidth3 / 2));
        double upper_boundary_int = floor(upper_boundary);
        int y = static_cast<int>(upper_boundary_int);
        
        double lower_boundary = (Used_size3 - (Used_binwidth3 / 2));
        double lower_boundary_int = floor(lower_boundary);
        int y0 = static_cast<int>(lower_boundary_int);
        if (y0 < -1) y0 = -1;
        
        double log_amu = log(alpha) + log(mu);
        double log_mid = -1.0 * theta * log(1.0 + (alpha * mu));
        double den_corr {1.0};
        if (modeltrunc == 1) den_corr = 1.0 - exp(log_mid);
        
        double current_prob {0.0};
        
        for (int summed_size = (y0 + 1); summed_size <= y; summed_size++) {
          double log_leftie = 0.0;
          for (int j = 0; j < summed_size; j++) {
            log_leftie = log(static_cast<double>(j) + theta) - log(static_cast<double>(j) + 1.0) + log_leftie;
          }
          double log_rightie = static_cast<double>(summed_size) * (log_amu - log(1.0 + (alpha * mu)));
          
          double raw_prob = log_leftie + log_mid + log_rightie;
          
          current_prob += exp(raw_prob) / den_corr;
        }
        all_out = current_prob;
        
        // Rcout << "Negbin: y: " << y << " y0: " << y0 << " alpha: " << alpha <<
        //   " mu: " << mu << " current_prob: " << current_prob << "\n";
        
      } else if (vitaldist == 2) {
        // Gaussian size distribution, assuming midpoint
        
        if (ipm_cdf) {
          double lower_size = Used_size3 - (0.5 * Used_binwidth3);
          double upper_size = Used_size3 + (0.5 * Used_binwidth3);
          
          double lower_prob = normcdf(lower_size, preout, Used_svsigmas1);
          double upper_prob = normcdf(upper_size, preout, Used_svsigmas1);
          
          all_out = upper_prob - lower_prob;
          
          // Rcout << "Gaussian cdf: upper_size: " << upper_size << " lower_size: " <<
          //   lower_size << " Used_size3: " << Used_size3 << " Used_binwidth3: " <<
          //   Used_binwidth3 << " preout: " << preout << " Used_svsigmas1: " <<
          //   Used_svsigmas1 << " upper_prob: " << upper_prob << " lower_prob: " <<
          //   lower_prob << " all_out: " << all_out << "\n";
        } else {
          double sigma2 = Used_svsigmas1 * Used_svsigmas1;
          
          all_out = (exp(-1 * (pow((Used_size3 - preout), 2) / (2.0 * sigma2))) / 
            ((pow((2 * M_PI), 0.5)) * Used_svsigmas1));
          all_out = all_out * Used_binwidth3; // This is the midpoint integration
          
          // Rcout << "Gaussian mid: Used_size3: " << Used_size3 << " Used_binwidth3: " <<
          //   Used_binwidth3 << " Used_svsigmas1: " << Used_svsigmas1 << " preout: " <<
          //   preout << " all_out: " << all_out << "\n";
        }
      } else if (vitaldist == 3) {
        // Gamma size distribution, assuming midpoint
        
        double E_y = 1 / preout;
        double sigma2 = Used_svsigmas1 * Used_svsigmas1;
        double alpha = 1.0 / sigma2;
        double beta = (alpha / E_y);
        
        if (ipm_cdf) {
          double lower_size = Used_size3 - (0.5 * Used_binwidth3);
          double upper_size = Used_size3 + (0.5 * Used_binwidth3);
          
          double lower_prob = boost::math::gamma_p(alpha, (beta * lower_size));
          double upper_prob = boost::math::gamma_p(alpha, (beta * upper_size));
          
          all_out = upper_prob - lower_prob;
          
          // Rcout << "Gamma cdf: upper_size: " << upper_size << " lower_size: " <<
          //   lower_size << " alpha: " << alpha << " beta: " << beta << " upper_prob: " <<
          //   upper_prob << " lower_prob: " << lower_prob << " all_out: " << all_out << "\n";
        } else {
          
          all_out = pow(beta, alpha) * (1.0 / tgamma(alpha)) * 
            pow(Used_size3, (alpha - 1.0)) * exp(-1.0 * beta * Used_size3);
          all_out = all_out * Used_binwidth3; // This is the midpoint integration
          
          // Rcout << "Gamma mid: Used_size3: " << Used_size3 << " Used_binwidth3: " <<
          //   Used_binwidth3 << " alpha: " << alpha << " beta: " << beta <<
          //   " all_out: " << all_out << "\n";
        }
      }
    }
  } else if (vitaltype == 2) {
    if (matrixformat != 2 || stage2n_i != static_cast<double>(nostages+1)) {
      if (vitaldist == 0 || vitaldist == 1) {
        // Poisson and negative binomial fecundity
        if (preout > exp_tol) preout = exp_tol;
        
        if (zi_processing) {
          
          all_out = (exp(preout) / (1.0 + exp(preout))) * fecmod * repentry_i;
          
        } else {
          
          all_out = exp(preout) * fecmod * repentry_i;
        }
      } else if (vitaldist == 2) {
        // Gaussian fecundity
        all_out = preout * fecmod * repentry_i;
        
        if (negfec && all_out < 0.0) all_out = 0.0;
        
      } else if (vitaldist == 3) {
        // Gamma fecundity
        all_out = (1.0 / preout) * fecmod * repentry_i;
      } else {
        all_out = maincoefs(0);
      }
    } else if (stage2n_i == static_cast<double>(nostages+1)) {
      // This propagates fecundity in deVries-formatted hMPMs
      if (vitaldist == 0 || vitaldist == 1) {
        // Poisson and negative binomial fecundity
        
        if (preout > exp_tol) preout = exp_tol;
            
        if (zi_processing) {
          
          all_out = (exp(preout) / (1.0 + exp(preout))) * fecmod * repentry_i;
          
        } else {
          
            all_out = exp(preout) * fecmod * repentry_i;
          
        }
      } else if (vitaldist == 2) {
        // Gaussian fecundity
        all_out = preout * fecmod * repentry_i;
        
        if (negfec && all_out < 0.0) {
          all_out = 0.0;
        }
      } else if (vitaldist == 3) {
        // Gamma fecundity
        all_out = (1.0 / preout) * fecmod * repentry_i;
      }
    } else {
      all_out = maincoefs(0);
    }
  }
  
  return(all_out);
}

//' Estimate All Elements of Function-based Population Projection Matrix
//' 
//' Function \code{.jerzeibalowski()} swiftly calculates matrix elements in
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
//' @param jmatstproxy List of coefficients estimated in model of juvenile
//' maturity probability.
//' @param f2_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{a} at
//' each time \emph{t} to be used in analysis.
//' @param f1_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{a} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{b} at
//' each time \emph{t} to be used in analysis.
//' @param f1_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{b} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{c} at
//' each time \emph{t} to be used in analysis.
//' @param f1_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{c} at
//' each time \emph{t}-1 to be used in analysis.
//' @param r2_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{a} at each time \emph{t} to be used in analysis.
//' @param r1_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{a} at each time \emph{t}-1 to be used in analysis.
//' @param r2_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{b} at each time \emph{t} to be used in analysis.
//' @param r1_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{b} at each time \emph{t}-1 to be used in analysis.
//' @param r2_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{c} at each time \emph{t} to be used in analysis.
//' @param r1_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{c} at each time \emph{t}-1 to be used in analysis.
//' @param dev_terms A numeric vector containing the deviations to the linear
//' models input by the user. The order is: survival, observation status, size,
//' size_b, size_c, reproductive status, fecundity, juvenile survival, juvenile
//' observation status, juvenile size, juvenile size_b, juvenile size_c,
//' juvenile reproductive status, and juvenile maturity status.
//' @param dens A numeric value equal to the density to be used in calculations.
//' @param fecmod A scalar multiplier for fecundity.
//' @param svsigmas A vector of sigma and summedvar terms from vital rate
//' models, in the order of: summedvars, sigma, summedvarsb, sigmab,
//' summedvarsc, sigmac, jsummedvars, jsigma, jsummedvarsb, jsigmab,
//' jsummedvarsc, and jsigmac. Summedvar terms are summed variance-covariance
//' terms in Poisson and negative binomial size distributions, and sigma terms
//' are standard deviations in the Gaussian size distribution.
//' @param maxsize The maximum primary size to be used in element estimation.
//' @param maxsizeb The maximum secondary size to be used in element estimation.
//' @param maxsizec The maximum tertiary size to be used in element estimation.
//' @param finalage The final age to be included in age-by-stage MPM estimation.
//' @param sizedist Designates whether primary size is Gamma (3), Gaussian (2),
//' Poisson (0), or negative binomial (1) distributed.
//' @param sizebdist Designates whether secondary size is Gamma (3),
//' Gaussian (2), Poisson (0), or negative binomial (1) distributed.
//' @param sizecdist Designates whether tertiary size is Gamma (3),
//' Gaussian (2), Poisson (0), or negative binomial (1) distributed.
//' @param fecdist Designates whether fecundity is Gamma (3), Gaussian (2),
//' Poisson (0), or negative binomial (1) distributed.
//' @param negfec A logical value denoting whether to change negative estimated
//' fecundity to 0.
//' @param exp_tol A numeric value indicating the maximum limit for the
//' \code{exp()} function to be used in vital rate calculations. Defaults to
//' \code{700.0}.
//' @param theta_tol A numeric value indicating a maximum value for theta in
//' negative binomial probability density estimation. Defaults to
//' \code{100000000.0}.
//' @param ipm_method A string indicating which method should be used to
//' estimate size transitions in cases with continuous distributions. Options
//' include \code{"midpoint"}, which uses the midpoint method, and \code{"cdf"},
//' which uses the cumulative density function.
//' 
//' @return A list with 4 elements. The first 3 elements are matrices, including
//' the main MPM (A), the survival-transition matrix (U), and a fecundity matrix
//' (F). The last element is a 6 column matrix showing survival probability,
//' observation probability, reproduction probability, sizea transition
//' probability, sizeb transition probability, and sizec transition probability
//' for each element of the final MPM.
//' 
//' @section Notes:
//' 
//' The DataFrame AllStages introduces variables used in size and fecundity
//' calculations. This DataFrame is broken up into long vectors composed of
//' input sizes and related variables for these calculations. The "model" Lists
//' bring in the vital rate models, and include random coefficients where
//' needed. We also have a number of extra variables, that include such info as
//' whether to use the Poisson, negative binomial, Gamma, or Gaussian
//' distributions for size and fecundity calculations. If \code{sizedist},
//' \code{sizebdist}, \code{sizecdist}, or \code{fecdist} equals 0, 1, 2, or 3,
//' then the Poisson, negative binomial, Gaussian, or Gamma is used,
//' respectively.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.jerzeibalowski)]]
List jerzeibalowski(DataFrame ppy, DataFrame AllStages, DataFrame stageframe,
  int matrixformat, List survproxy, List obsproxy, List sizeproxy,
  List sizebproxy, List sizecproxy, List repstproxy, List fecproxy,
  List jsurvproxy, List jobsproxy, List jsizeproxy, List jsizebproxy,
  List jsizecproxy, List jrepstproxy, List jmatstproxy, NumericVector f2_inda,
  NumericVector f1_inda, NumericVector f2_indb, NumericVector f1_indb,
  NumericVector f2_indc, NumericVector f1_indc, StringVector r2_inda,
  StringVector r1_inda, StringVector r2_indb, StringVector r1_indb,
  StringVector r2_indc, StringVector r1_indc, NumericVector dev_terms,
  double dens, double fecmod, NumericVector svsigmas, double maxsize,
  double maxsizeb, double maxsizec, unsigned int finalage, int sizedist,
  int sizebdist, int sizecdist, int fecdist, bool negfec,
  double exp_tol = 700.0, double theta_tol = 100000000.0,
  String ipm_method = "cdf") {
  
  bool ipm_cdf = true;
  if (ipm_method == "midpoint") ipm_cdf = false;
  
  // Determines the size of the matrix
  StringVector stagenames = stageframe["stage"];
  int nostages = stagenames.length();
  unsigned long matrixdim {0};
  
  int nostages_counter = nostages;
  for (int i = 0; i < nostages_counter; i++) {
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
  
  // Proxy model imports and settings
  bool sizezero = false;
  bool sizebzero = false;
  bool sizeczero = false;
  bool feczero = false;
  bool jsizezero = false;
  bool jsizebzero = false;
  bool jsizeczero = false;
  
  NumericVector survcoefs = survproxy["coefficients"];
  NumericVector obscoefs = obsproxy["coefficients"];
  NumericVector sizecoefs = sizeproxy["coefficients"];
  NumericVector sizebcoefs = sizebproxy["coefficients"];
  NumericVector sizeccoefs = sizecproxy["coefficients"];
  NumericVector repstcoefs = repstproxy["coefficients"];
  NumericVector feccoefs = fecproxy["coefficients"];
  NumericVector jsurvcoefs = jsurvproxy["coefficients"];
  NumericVector jobscoefs = jobsproxy["coefficients"];
  NumericVector jsizecoefs = jsizeproxy["coefficients"];
  NumericVector jsizebcoefs = jsizebproxy["coefficients"];
  NumericVector jsizeccoefs = jsizecproxy["coefficients"];
  NumericVector jrepstcoefs = jrepstproxy["coefficients"];
  NumericVector jmatstcoefs = jmatstproxy["coefficients"];
  
  int survl = survcoefs.length();
  int obsl = obscoefs.length();
  int sizel = sizecoefs.length();
  int sizebl = sizebcoefs.length();
  int sizecl = sizeccoefs.length();
  int repstl = repstcoefs.length();
  int fecl = feccoefs.length();
  int jsurvl = jsurvcoefs.length();
  int jobsl = jobscoefs.length();
  int jsizel = jsizecoefs.length();
  int jsizebl = jsizebcoefs.length();
  int jsizecl = jsizeccoefs.length();
  int jrepstl = jrepstcoefs.length();
  int jmatstl = jmatstcoefs.length();
  
  int sizetrunc = sizeproxy["trunc"];
  int sizebtrunc = sizebproxy["trunc"];
  int sizectrunc = sizecproxy["trunc"];
  int fectrunc = fecproxy["trunc"];
  int jsizetrunc = jsizeproxy["trunc"];
  int jsizebtrunc = jsizebproxy["trunc"];
  int jsizectrunc = jsizecproxy["trunc"];
  
  double sizesigma = sizeproxy["sigma"];
  double sizebsigma = sizebproxy["sigma"];
  double sizecsigma = sizecproxy["sigma"];
  double fecsigma = fecproxy["sigma"];
  double jsizesigma = jsizeproxy["sigma"];
  double jsizebsigma = jsizebproxy["sigma"];
  double jsizecsigma = jsizecproxy["sigma"];
  
  if (NumericVector::is_na(sizesigma)) {
    if (sizedist == 1) {
      sizesigma = 1.0;
    } else {
      sizesigma = 0.0;
    }
  }
  if (NumericVector::is_na(sizebsigma)) {
    if (sizebdist == 1) {
      sizebsigma = 1.0;
    } else {
      sizebsigma = 0.0;
    }
  }
  if (NumericVector::is_na(sizecsigma)) {
    if (sizecdist == 1) {
      sizecsigma = 1.0;
    } else {
      sizecsigma = 0.0;
    }
  }
  if (NumericVector::is_na(fecsigma)) {
    if (fecdist == 1) {
      fecsigma = 1.0;
    } else {
      fecsigma = 0.0;
    }
  }
  if (NumericVector::is_na(jsizesigma)) {
    if (sizedist == 1) {
      jsizesigma = 1.0;
    } else {
      jsizesigma = 0.0;
    }
  }
  if (NumericVector::is_na(jsizebsigma)) {
    if (sizebdist == 1) {
      jsizebsigma = 1.0;
    } else {
      jsizebsigma = 0.0;
    }
  }
  if (NumericVector::is_na(jsizecsigma)) {
    if (sizecdist == 1) {
      jsizecsigma = 1.0;
    } else {
      jsizecsigma = 0.0;
    }
  }
  
  NumericMatrix vital_year = revelations(survproxy, obsproxy, sizeproxy, sizebproxy,
    sizecproxy, repstproxy, fecproxy, jsurvproxy, jobsproxy, jsizeproxy,
    jsizebproxy, jsizecproxy, jrepstproxy, jmatstproxy, 1);
  
  NumericMatrix vital_patch = revelations(survproxy, obsproxy, sizeproxy,
    sizebproxy, sizecproxy, repstproxy, fecproxy, jsurvproxy, jobsproxy,
    jsizeproxy, jsizebproxy, jsizecproxy, jrepstproxy, jmatstproxy, 2);
  
  arma::imat rand_index = foi_index(survproxy, obsproxy, sizeproxy, sizebproxy,
    sizecproxy, repstproxy, fecproxy, jsurvproxy, jobsproxy, jsizeproxy,
    jsizebproxy, jsizecproxy, jrepstproxy, jmatstproxy);
  
  Rcpp::DataFrame sizeyearzi_df(sizeproxy["zeroyear"]);
  Rcpp::DataFrame sizebyearzi_df(sizebproxy["zeroyear"]);
  Rcpp::DataFrame sizecyearzi_df(sizecproxy["zeroyear"]);
  Rcpp::DataFrame fecyearzi_df(fecproxy["zeroyear"]);
  Rcpp::DataFrame jsizeyearzi_df(jsizeproxy["zeroyear"]);
  Rcpp::DataFrame jsizebyearzi_df(jsizebproxy["zeroyear"]);
  Rcpp::DataFrame jsizecyearzi_df(jsizecproxy["zeroyear"]);
  
  NumericVector sizeyearzi = sizeyearzi_df[0];
  NumericVector sizebyearzi = sizebyearzi_df[0];
  NumericVector sizecyearzi = sizecyearzi_df[0];
  NumericVector fecyearzi = fecyearzi_df[0];
  NumericVector jsizeyearzi = jsizeyearzi_df[0];
  NumericVector jsizebyearzi = jsizebyearzi_df[0];
  NumericVector jsizecyearzi = jsizecyearzi_df[0];
  
  NumericVector dud_yearzi(sizeyearzi.length());
  
  NumericVector unisyzi = unique(sizeyearzi);
  NumericVector unisyzbi = unique(sizebyearzi);
  NumericVector unisyzci = unique(sizecyearzi);
  NumericVector unijsyzi = unique(jsizeyearzi);
  NumericVector unijsyzbi = unique(jsizebyearzi);
  NumericVector unijsyzci = unique(jsizecyearzi);
  NumericVector unifeci = unique(fecyearzi);
  
  if (unisyzi.length() > 1 || unisyzi(0) != 0) sizezero = true;
  if (unisyzbi.length() > 1 || unisyzbi(0) != 0) sizebzero = true;
  if (unisyzci.length() > 1 || unisyzci(0) != 0) sizeczero = true;
  if (unifeci.length() > 1 || unifeci(0) != 0) feczero = true;
  if (unijsyzi.length() > 1 || unijsyzi(0) != 0) jsizezero = true;
  if (unijsyzbi.length() > 1 || unijsyzbi(0) != 0) jsizebzero = true;
  if (unijsyzci.length() > 1 || unijsyzci(0) != 0) jsizeczero = true;
  
  Rcpp::DataFrame sizepatchzi_df(sizeproxy["zeropatch"]);
  Rcpp::DataFrame sizebpatchzi_df(sizebproxy["zeropatch"]);
  Rcpp::DataFrame sizecpatchzi_df(sizecproxy["zeropatch"]);
  Rcpp::DataFrame fecpatchzi_df(fecproxy["zeropatch"]);
  Rcpp::DataFrame jsizepatchzi_df(jsizeproxy["zeropatch"]);
  Rcpp::DataFrame jsizebpatchzi_df(jsizebproxy["zeropatch"]);
  Rcpp::DataFrame jsizecpatchzi_df(jsizecproxy["zeropatch"]);
  
  NumericVector sizepatchzi = sizepatchzi_df[0];
  NumericVector sizebpatchzi = sizebpatchzi_df[0];
  NumericVector sizecpatchzi = sizecpatchzi_df[0];
  NumericVector fecpatchzi = fecpatchzi_df[0];
  NumericVector jsizepatchzi = jsizepatchzi_df[0];
  NumericVector jsizebpatchzi = jsizebpatchzi_df[0];
  NumericVector jsizecpatchzi = jsizecpatchzi_df[0];
  
  NumericVector dud_patchzi(sizepatchzi.length());
  
  Rcpp::DataFrame survgroups2_df(survproxy["groups2"]);
  Rcpp::DataFrame obsgroups2_df(obsproxy["groups2"]);
  Rcpp::DataFrame sizegroups2_df(sizeproxy["groups2"]);
  Rcpp::DataFrame sizebgroups2_df(sizebproxy["groups2"]);
  Rcpp::DataFrame sizecgroups2_df(sizecproxy["groups2"]);
  Rcpp::DataFrame repstgroups2_df(repstproxy["groups2"]);
  Rcpp::DataFrame fecgroups2_df(fecproxy["groups2"]);
  Rcpp::DataFrame jsurvgroups2_df(jsurvproxy["groups2"]);
  Rcpp::DataFrame jobsgroups2_df(jobsproxy["groups2"]);
  Rcpp::DataFrame jsizegroups2_df(jsizeproxy["groups2"]);
  Rcpp::DataFrame jsizebgroups2_df(jsizebproxy["groups2"]);
  Rcpp::DataFrame jsizecgroups2_df(jsizecproxy["groups2"]);
  Rcpp::DataFrame jrepstgroups2_df(jrepstproxy["groups2"]);
  Rcpp::DataFrame jmatstgroups2_df(jmatstproxy["groups2"]);
  
  NumericVector survgroups2 = survgroups2_df[0];
  NumericVector obsgroups2 = obsgroups2_df[0];
  NumericVector sizegroups2 = sizegroups2_df[0];
  NumericVector sizebgroups2 = sizebgroups2_df[0];
  NumericVector sizecgroups2 = sizecgroups2_df[0];
  NumericVector repstgroups2 = repstgroups2_df[0];
  NumericVector fecgroups2 = fecgroups2_df[0];
  NumericVector jsurvgroups2 = jsurvgroups2_df[0];
  NumericVector jobsgroups2 = jobsgroups2_df[0];
  NumericVector jsizegroups2 = jsizegroups2_df[0];
  NumericVector jsizebgroups2 = jsizebgroups2_df[0];
  NumericVector jsizecgroups2 = jsizecgroups2_df[0];
  NumericVector jrepstgroups2 = jrepstgroups2_df[0];
  NumericVector jmatstgroups2 = jmatstgroups2_df[0];
  
  Rcpp::DataFrame survgroups1_df(survproxy["groups1"]);
  Rcpp::DataFrame obsgroups1_df(obsproxy["groups1"]);
  Rcpp::DataFrame sizegroups1_df(sizeproxy["groups1"]);
  Rcpp::DataFrame sizebgroups1_df(sizebproxy["groups1"]);
  Rcpp::DataFrame sizecgroups1_df(sizecproxy["groups1"]);
  Rcpp::DataFrame repstgroups1_df(repstproxy["groups1"]);
  Rcpp::DataFrame fecgroups1_df(fecproxy["groups1"]);
  Rcpp::DataFrame jsurvgroups1_df(jsurvproxy["groups1"]);
  Rcpp::DataFrame jobsgroups1_df(jobsproxy["groups1"]);
  Rcpp::DataFrame jsizegroups1_df(jsizeproxy["groups1"]);
  Rcpp::DataFrame jsizebgroups1_df(jsizebproxy["groups1"]);
  Rcpp::DataFrame jsizecgroups1_df(jsizecproxy["groups1"]);
  Rcpp::DataFrame jrepstgroups1_df(jrepstproxy["groups1"]);
  Rcpp::DataFrame jmatstgroups1_df(jmatstproxy["groups1"]);
  
  NumericVector survgroups1 = survgroups1_df[0];
  NumericVector obsgroups1 = obsgroups1_df[0];
  NumericVector sizegroups1 = sizegroups1_df[0];
  NumericVector sizebgroups1 = sizebgroups1_df[0];
  NumericVector sizecgroups1 = sizecgroups1_df[0];
  NumericVector repstgroups1 = repstgroups1_df[0];
  NumericVector fecgroups1 = fecgroups1_df[0];
  NumericVector jsurvgroups1 = jsurvgroups1_df[0];
  NumericVector jobsgroups1 = jobsgroups1_df[0];
  NumericVector jsizegroups1 = jsizegroups1_df[0];
  NumericVector jsizebgroups1 = jsizebgroups1_df[0];
  NumericVector jsizecgroups1 = jsizecgroups1_df[0];
  NumericVector jrepstgroups1 = jrepstgroups1_df[0];
  NumericVector jmatstgroups1 = jmatstgroups1_df[0];
  
  Rcpp::DataFrame sizegroups2zi_df(sizeproxy["zerogroups2"]);
  Rcpp::DataFrame sizebgroups2zi_df(sizebproxy["zerogroups2"]);
  Rcpp::DataFrame sizecgroups2zi_df(sizecproxy["zerogroups2"]);
  Rcpp::DataFrame fecgroups2zi_df(fecproxy["zerogroups2"]);
  Rcpp::DataFrame jsizegroups2zi_df(jsizeproxy["zerogroups2"]);
  Rcpp::DataFrame jsizebgroups2zi_df(jsizebproxy["zerogroups2"]);
  Rcpp::DataFrame jsizecgroups2zi_df(jsizecproxy["zerogroups2"]);
  
  NumericVector sizegroups2zi = sizegroups2zi_df[0];
  NumericVector sizebgroups2zi = sizebgroups2zi_df[0];
  NumericVector sizecgroups2zi = sizecgroups2zi_df[0];
  NumericVector fecgroups2zi = fecgroups2zi_df[0];
  NumericVector jsizegroups2zi = jsizegroups2zi_df[0];
  NumericVector jsizebgroups2zi = jsizebgroups2zi_df[0];
  NumericVector jsizecgroups2zi = jsizecgroups2zi_df[0];
  
  NumericVector dud_groups2zi(jsizecyearzi.length());
  
  Rcpp::DataFrame sizegroups1zi_df(sizeproxy["zerogroups1"]);
  Rcpp::DataFrame sizebgroups1zi_df(sizebproxy["zerogroups1"]);
  Rcpp::DataFrame sizecgroups1zi_df(sizecproxy["zerogroups1"]);
  Rcpp::DataFrame fecgroups1zi_df(fecproxy["zerogroups1"]);
  Rcpp::DataFrame jsizegroups1zi_df(jsizeproxy["zerogroups1"]);
  Rcpp::DataFrame jsizebgroups1zi_df(jsizebproxy["zerogroups1"]);
  Rcpp::DataFrame jsizecgroups1zi_df(jsizecproxy["zerogroups1"]);
  
  NumericVector sizegroups1zi = sizegroups1zi_df[0];
  NumericVector sizebgroups1zi = sizebgroups1zi_df[0];
  NumericVector sizecgroups1zi = sizecgroups1zi_df[0];
  NumericVector fecgroups1zi = fecgroups1zi_df[0];
  NumericVector jsizegroups1zi = jsizegroups1zi_df[0];
  NumericVector jsizebgroups1zi = jsizebgroups1zi_df[0];
  NumericVector jsizecgroups1zi = jsizecgroups1zi_df[0];
  
  NumericVector dud_groups1zi(jsizecyearzi.length());
  
  NumericVector survind = flightoficarus(survproxy);
  NumericVector obsind = flightoficarus(obsproxy);
  NumericVector sizeind = flightoficarus(sizeproxy);
  NumericVector sizebind = flightoficarus(sizebproxy);
  NumericVector sizecind = flightoficarus(sizecproxy);
  NumericVector repstind = flightoficarus(repstproxy);
  NumericVector fecind = flightoficarus(fecproxy);
  NumericVector jsurvind = flightoficarus(jsurvproxy);
  NumericVector jobsind = flightoficarus(jobsproxy);
  NumericVector jsizeind = flightoficarus(jsizeproxy);
  NumericVector jsizebind = flightoficarus(jsizebproxy);
  NumericVector jsizecind = flightoficarus(jsizecproxy);
  NumericVector jrepstind = flightoficarus(jrepstproxy);
  NumericVector jmatstind = flightoficarus(jmatstproxy);
  
  NumericVector sizeindzi = zero_flightoficarus(sizeproxy);
  NumericVector sizebindzi = zero_flightoficarus(sizebproxy);
  NumericVector sizecindzi = zero_flightoficarus(sizecproxy);
  NumericVector fecindzi = zero_flightoficarus(fecproxy);
  NumericVector jsizeindzi = zero_flightoficarus(jsizeproxy);
  NumericVector jsizebindzi = zero_flightoficarus(jsizebproxy);
  NumericVector jsizecindzi = zero_flightoficarus(jsizecproxy);
  
  StringVector survind_rownames = bootson(survproxy);
  StringVector obsind_rownames = bootson(obsproxy);
  StringVector sizeind_rownames = bootson(sizeproxy);
  StringVector sizebind_rownames = bootson(sizebproxy);
  StringVector sizecind_rownames = bootson(sizecproxy);
  StringVector repstind_rownames = bootson(repstproxy);
  StringVector fecind_rownames = bootson(fecproxy);
  StringVector jsurvind_rownames = bootson(jsurvproxy);
  StringVector jobsind_rownames = bootson(jobsproxy);
  StringVector jsizeind_rownames = bootson(jsizeproxy);
  StringVector jsizebind_rownames = bootson(jsizebproxy);
  StringVector jsizecind_rownames = bootson(jsizecproxy);
  StringVector jrepstind_rownames = bootson(jrepstproxy);
  StringVector jmatstind_rownames = bootson(jmatstproxy);
  
  StringVector sizeind_rownames_zi = zero_bootson(sizeproxy);
  StringVector sizebind_rownames_zi = zero_bootson(sizebproxy);
  StringVector sizecind_rownames_zi = zero_bootson(sizecproxy);
  StringVector fecind_rownames_zi = zero_bootson(fecproxy);
  StringVector jsizeind_rownames_zi = zero_bootson(jsizeproxy);
  StringVector jsizebind_rownames_zi = zero_bootson(jsizebproxy);
  StringVector jsizecind_rownames_zi = zero_bootson(jsizecproxy);
  
  // AllStages import and settings
  Rcpp::NumericVector stage3_num = AllStages["stage3"];
  Rcpp::NumericVector stage2n_num = AllStages["stage2n"];
  Rcpp::NumericVector stage2o_num = AllStages["stage2o"];
  arma::vec stage3 = as<arma::vec>(stage3_num);
  arma::vec stage2n = as<arma::vec>(stage2n_num);
  arma::vec stage2o = as<arma::vec>(stage2o_num);
  
  Rcpp::NumericVector sz3 = AllStages["size3"];
  Rcpp::NumericVector sz2n = AllStages["size2n"];
  Rcpp::NumericVector sz2o = AllStages["size2o"];
  Rcpp::NumericVector sz1 = AllStages["size1"];
  Rcpp::NumericVector szb3 = AllStages["sizeb3"];
  Rcpp::NumericVector szb2n = AllStages["sizeb2n"];
  Rcpp::NumericVector szb2o = AllStages["sizeb2o"];
  Rcpp::NumericVector szb1 = AllStages["sizeb1"];
  Rcpp::NumericVector szc3 = AllStages["sizec3"];
  Rcpp::NumericVector szc2n = AllStages["sizec2n"];
  Rcpp::NumericVector szc2o = AllStages["sizec2o"];
  Rcpp::NumericVector szc1 = AllStages["sizec1"];
  Rcpp::NumericVector ob3 = AllStages["obs3"];
  Rcpp::NumericVector fl3 = AllStages["rep3"];
  Rcpp::NumericVector fl2n = AllStages["rep2n"];
  Rcpp::NumericVector fl2o = AllStages["rep2o"];
  Rcpp::NumericVector fl1 = AllStages["rep1"];
  Rcpp::NumericVector mat3 = AllStages["mat3"];
  Rcpp::NumericVector mat2n = AllStages["mat2n"];
  Rcpp::NumericVector mat2o = AllStages["mat2o"];
  Rcpp::NumericVector mat1 = AllStages["mat1"];
  Rcpp::NumericVector immat2n = AllStages["imm2n"];
  Rcpp::NumericVector immat2o = AllStages["imm2o"];
  Rcpp::NumericVector immat1 = AllStages["imm1"];
  
  Rcpp::NumericVector repentry = AllStages["repentry3"];
  Rcpp::NumericVector indata2n = AllStages["indata2n"];
  Rcpp::NumericVector indata2o = AllStages["indata2o"];
  Rcpp::NumericVector binwidth3 = AllStages["binwidth"];
  Rcpp::NumericVector binbwidth3 = AllStages["binbwidth"];
  Rcpp::NumericVector bincwidth3 = AllStages["bincwidth"];
  Rcpp::NumericVector actualage2 = AllStages["actualage"];
  
  Rcpp::NumericVector grp3 = AllStages["group3"];
  Rcpp::NumericVector grp2n = AllStages["group2n"];
  Rcpp::NumericVector grp2o = AllStages["group2o"];
  Rcpp::NumericVector grp1 = AllStages["group1"];
  
  Rcpp::NumericVector ovestt_num = AllStages["ovest_t"];
  arma::vec ovestt = as<arma::vec>(ovestt_num);
  
  Rcpp::NumericVector ovestf_num = AllStages["ovest_f"];
  arma::vec ovestf = as<arma::vec>(ovestf_num);
  
  Rcpp::NumericVector indata = AllStages["indata"];
  Rcpp::NumericVector ovgivent = AllStages["ovgiven_t"];
  Rcpp::NumericVector ovgivenf = AllStages["ovgiven_f"];
  
  Rcpp::NumericVector ovsurvmult = AllStages["ovsurvmult"];
  Rcpp::NumericVector ovfecmult = AllStages["ovfecmult"];
  
  Rcpp::IntegerVector index321_int = AllStages["index321"];
  arma::uvec index321 = as<arma::uvec>(index321_int);
  
  Rcpp::NumericVector aliveandequal = AllStages["aliveandequal"];
  
  int n = stage3.n_elem;
  
  arma::uvec replacetvec = find(ovestt != -1);
  arma::uvec replacefvec = find(ovestf != -1);
  int replacementst = replacetvec.n_elem;
  int replacementsf = replacefvec.n_elem;
  int repindex {0};
  int properindex {0};
  int proxyindex {0};
  
  // listofyears import and settings
  arma::uvec years = ppy["yearorder"];
  int yearnumber = years(0) - 1;
  
  arma::uvec patches = ppy["patchorder"];
  int patchnumber = patches(0) - 1;
  
  // Determination of choices of fixed and random individual covariates
  double inda1 = f1_inda(yearnumber);
  double indb1 = f1_indb(yearnumber);
  double indc1 = f1_indc(yearnumber);
  double inda2 = f2_inda(yearnumber);
  double indb2 = f2_indb(yearnumber);
  double indc2 = f2_indc(yearnumber);
  
  String chosen_r2inda = r2_inda(yearnumber);
  String chosen_r1inda = r1_inda(yearnumber);
  String chosen_r2indb = r2_indb(yearnumber);
  String chosen_r1indb = r1_indb(yearnumber);
  String chosen_r2indc = r2_indc(yearnumber);
  String chosen_r1indc = r1_indc(yearnumber);
  
  // The output matrix to collect conditional probabilities
  // Matrix out is 0 matrix with n rows & 6 columns: 0 surv, 1 obs, 2 repst,
  // 3 size, 4 size_b, 5 size_c, 6 matst, >6 are test variables
  arma::mat out(n, 7);
  arma::mat survtransmat(matrixdim, matrixdim);
  arma::mat fectransmat(matrixdim, matrixdim);
  
  out.zeros();
  survtransmat.zeros();
  fectransmat.zeros();
  
  // The following loop runs through each line of AllStages, and so runs through
  // each estimable element in the matrix
  for(int i = 0; i < n; i++) {
    unsigned int k = aliveandequal(i);
    
    out(i, 6) = 1.0; // Initialization of maturity status probability for typical case
    
    Rcpp::NumericVector statusterms = {fl1(i), fl2n(i), sz1(i), sz2o(i),
      szb1(i), szb2o(i), szc1(i), szc2o(i), actualage2(i), inda1, inda2, indb1,
      indb2, indc1, indc2, dens, sz3(i), szb3(i), szc3(i), binwidth3(i),
      binbwidth3(i), bincwidth3(i)};
    
    if (ovgivent(i) == -1 && indata(i) == 1 && stage2n(i) == stage2o(i)) {
      if ((mat2n(i) == 1 && mat3(i) == 1) || (mat2o(i) == 1 && mat3(i) == 1)) {
        
        // Adult survival transitions
        if (survl > 1) {
          out(i, 0) = preouterator(survproxy, survcoefs, rand_index, dev_terms,
            svsigmas, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
            chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
            statusterms, survgroups2, survgroups1, dud_groups2zi, dud_groups1zi,
            dud_yearzi, dud_patchzi, survind, survind_rownames, sizeindzi,
            sizeind_rownames_zi, false, grp2o(i), grp1(i), patchnumber,
            yearnumber, 4, 1, exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod,
            repentry(i), negfec, stage2n(i), nostages, 0, 0.0);
          
        } else {
          out(i, 0) = survcoefs(0);
        }
        
        if (obsl > 1) {
          out(i, 1) = preouterator(obsproxy, obscoefs, rand_index, dev_terms,
            svsigmas, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
            chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
            statusterms, obsgroups2, obsgroups1, dud_groups2zi, dud_groups1zi,
            dud_yearzi, dud_patchzi, obsind, obsind_rownames, sizeindzi,
            sizeind_rownames_zi, false, grp2o(i), grp1(i), patchnumber,
            yearnumber, 4, 2, exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod,
            repentry(i), negfec, stage2n(i), nostages, 0, 0.0);
          
        } else {
          out(i, 1) = obscoefs(0);
        }
        
        if (ob3(i) == 1 || obsl == 1) {
          
          if (sizel > 1) {
            bool used_sizezero = false;
            if (sizezero && sz3(i) == 0) used_sizezero = sizezero;
            
            out(i, 3) = preouterator(sizeproxy, sizecoefs, rand_index,
              dev_terms, svsigmas, vital_year, vital_patch, chosen_r2inda,
              chosen_r1inda, chosen_r2indb, chosen_r1indb, chosen_r2indc,
              chosen_r1indc, statusterms, sizegroups2, sizegroups1,
              sizegroups2zi, sizegroups1zi, sizeyearzi, sizepatchzi, sizeind,
              sizeind_rownames, sizeindzi, sizeind_rownames_zi, used_sizezero,
              grp2o(i), grp1(i), patchnumber, yearnumber, sizedist, 3, exp_tol,
              theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
              stage2n(i), nostages, sizetrunc, sizesigma);
              
          } else {
            out(i, 3) = 1.0;
          }
          
          if (sizebl > 1) {
            bool used_sizebzero = false;
            if (sizebzero && szb3(i) == 0) used_sizebzero = sizebzero;
            
            out(i, 4) = preouterator(sizebproxy, sizebcoefs, rand_index,
              dev_terms, svsigmas, vital_year, vital_patch, chosen_r2inda,
              chosen_r1inda, chosen_r2indb, chosen_r1indb, chosen_r2indc,
              chosen_r1indc, statusterms, sizebgroups2, sizebgroups1,
              sizebgroups2zi, sizebgroups1zi, sizebyearzi, sizebpatchzi,
              sizebind, sizebind_rownames, sizebindzi, sizebind_rownames_zi,
              used_sizebzero, grp2o(i), grp1(i), patchnumber, yearnumber, sizebdist,
              4, exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i),
              negfec, stage2n(i), nostages, sizebtrunc, sizebsigma);
          } else {
            out(i, 4) = 1.0;
          }
          
          if (sizecl > 1) {
            bool used_sizeczero = false;
            if (sizeczero && szc3(i) == 0) used_sizeczero = sizeczero;
            
            out(i, 5) = preouterator(sizecproxy, sizeccoefs, rand_index,
              dev_terms, svsigmas, vital_year, vital_patch, chosen_r2inda,
              chosen_r1inda, chosen_r2indb, chosen_r1indb, chosen_r2indc,
              chosen_r1indc, statusterms, sizecgroups2, sizecgroups1,
              sizecgroups2zi, sizecgroups1zi, sizecyearzi, sizecpatchzi,
              sizecind, sizecind_rownames, sizecindzi, sizecind_rownames_zi,
              used_sizeczero, grp2o(i), grp1(i), patchnumber, yearnumber, sizecdist,
              5, exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i),
              negfec, stage2n(i), nostages, sizectrunc, sizecsigma);
          } else {
            out(i, 5) = 1.0;
          }
          
          if (repstl > 1) {
            out(i, 2) = preouterator(repstproxy, repstcoefs, rand_index,
              dev_terms, svsigmas, vital_year, vital_patch, chosen_r2inda,
              chosen_r1inda, chosen_r2indb, chosen_r1indb, chosen_r2indc,
              chosen_r1indc, statusterms, repstgroups2, repstgroups1,
              dud_groups2zi, dud_groups1zi, dud_yearzi, dud_patchzi, repstind,
              repstind_rownames, sizeindzi, sizeind_rownames_zi, false,
              grp2o(i), grp1(i), patchnumber, yearnumber, 4, 6, exp_tol,
              theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
              stage2n(i), nostages, 0, 0.0);
              
            if (fl3(i) == 0) {
              out(i, 2) = 1.0 - out(i, 2);
            }
          } else {
            if (fl3(i) == 0) {
              out(i, 2) = 1.0 - repstcoefs(0);
            } else if (fl3(i) == 1) {
              out(i, 2) = repstcoefs(0);
            } else {
              out(i, 2) = 0.0;
            }
          }
        } else {
          out(i, 1) = 1.0 - out(i, 1);
          out(i, 2) = 1.0;
          out(i, 3) = 1.0;
          out(i, 4) = 1.0;
          out(i, 5) = 1.0;
          out(i, 6) = 1.0;
        }
        
        survtransmat(k) = out(i, 0) * out(i, 1) * out(i, 2) * out(i, 3) *
          out(i, 4) * out(i, 5) * out(i, 6);
        
      } else if (immat2n(i) == 1 && immat1(i) == 1 && jsurvl > 0) {
        // Juvenile to adult transitions
        
        if (jmatstl > 1) {
          double mat_predicted = preouterator(jmatstproxy, jmatstcoefs,
            rand_index, dev_terms, svsigmas, vital_year, vital_patch,
            chosen_r2inda, chosen_r1inda, chosen_r2indb, chosen_r1indb,
            chosen_r2indc, chosen_r1indc, statusterms, jmatstgroups2,
            jmatstgroups1, dud_groups2zi, dud_groups1zi, dud_yearzi,
            dud_patchzi, jmatstind, jmatstind_rownames, jsizeindzi,
            jsizeind_rownames_zi, false, grp2o(i), grp1(i), patchnumber,
            yearnumber, 4, 14, exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod,
            repentry(i), negfec, stage2n(i), nostages, 0, 0.0);
          
          if (mat3(i) > 0.5) {
            out(i, 6) = mat_predicted;
          } else {
            out(i, 6) = 1 - mat_predicted;
          }
        } else {
          if (mat3(i) > 0.5) {
            out(i, 6) = 1;
          } else {
            out(i, 6) = 0;
          }
        }
        
        if (jsurvl > 1) {
          out(i, 0) = preouterator(jsurvproxy, jsurvcoefs, rand_index,
            dev_terms, svsigmas, vital_year, vital_patch, chosen_r2inda,
            chosen_r1inda, chosen_r2indb, chosen_r1indb, chosen_r2indc,
            chosen_r1indc, statusterms, jsurvgroups2, jsurvgroups1,
            dud_groups2zi, dud_groups1zi, dud_yearzi, dud_patchzi, jsurvind,
            jsurvind_rownames, jsizeindzi, jsizeind_rownames_zi, false,
            grp2o(i), grp1(i), patchnumber, yearnumber, 4, 8, exp_tol,
            theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
            stage2n(i), nostages, 0, 0.0);
        } else {
          out(i, 0) = jsurvcoefs(0);
        }
        
        if (jobsl > 1) {
          out(i, 1) = preouterator(jobsproxy, jobscoefs, rand_index, dev_terms,
            svsigmas, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
            chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
            statusterms, jobsgroups2, jobsgroups1, dud_groups2zi, dud_groups1zi,
            dud_yearzi, dud_patchzi, jobsind, jobsind_rownames, jsizeindzi,
            jsizeind_rownames_zi, false, grp2o(i), grp1(i), patchnumber,
            yearnumber, 4, 9, exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod,
            repentry(i), negfec, stage2n(i), nostages, 0, 0.0);
        } else {
          out(i, 1) = jobscoefs(0);
        }
        
        if (ob3(i) == 1 || jobsl == 1) {
          if (jsizel > 1) {
            out(i, 3) = preouterator(jsizeproxy, jsizecoefs, rand_index,
              dev_terms, svsigmas, vital_year, vital_patch, chosen_r2inda,
              chosen_r1inda, chosen_r2indb, chosen_r1indb, chosen_r2indc,
              chosen_r1indc, statusterms, jsizegroups2, jsizegroups1,
              jsizegroups2zi, jsizegroups1zi, jsizeyearzi, jsizepatchzi,
              jsizeind, jsizeind_rownames, jsizeindzi, jsizeind_rownames_zi,
              jsizezero, grp2o(i), grp1(i), patchnumber, yearnumber, sizedist,
              10, exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod,
              repentry(i), negfec, stage2n(i), nostages, jsizetrunc,
              jsizesigma);
          } else {
            out(i, 3) = 1.0;
          }
          
          if (jsizebl > 1) {
            out(i, 4) = preouterator(jsizebproxy, jsizebcoefs, rand_index,
              dev_terms, svsigmas, vital_year, vital_patch, chosen_r2inda,
              chosen_r1inda, chosen_r2indb, chosen_r1indb, chosen_r2indc,
              chosen_r1indc, statusterms, jsizebgroups2, jsizebgroups1,
              jsizebgroups2zi, jsizebgroups1zi, jsizebyearzi, jsizebpatchzi,
              jsizebind, jsizebind_rownames, jsizebindzi, jsizebind_rownames_zi,
              jsizebzero, grp2o(i), grp1(i), patchnumber, yearnumber, sizebdist,
              11, exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod,
              repentry(i), negfec, stage2n(i), nostages, jsizebtrunc,
              jsizebsigma);
          } else {
            out(i, 4) = 1.0;
          }
          
          if (jsizecl > 1) {
            out(i, 5) = preouterator(jsizecproxy, jsizeccoefs, rand_index,
              dev_terms, svsigmas, vital_year, vital_patch, chosen_r2inda,
              chosen_r1inda, chosen_r2indb, chosen_r1indb, chosen_r2indc,
              chosen_r1indc, statusterms, jsizecgroups2, jsizecgroups1,
              jsizecgroups2zi, jsizecgroups1zi, jsizecyearzi, jsizecpatchzi,
              jsizecind, jsizecind_rownames, jsizecindzi, jsizecind_rownames_zi,
              jsizeczero, grp2o(i), grp1(i), patchnumber, yearnumber, sizecdist,
              12, exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod,
              repentry(i), negfec, stage2n(i), nostages, jsizectrunc,
              jsizecsigma);
          } else {
            out(i, 5) = 1.0;
          }
          
          if (jrepstl > 1) {
            out(i, 2) = preouterator(jrepstproxy, jrepstcoefs, rand_index,
              dev_terms, svsigmas, vital_year, vital_patch, chosen_r2inda,
              chosen_r1inda, chosen_r2indb, chosen_r1indb, chosen_r2indc,
              chosen_r1indc, statusterms, jrepstgroups2, jrepstgroups1,
              dud_groups2zi, dud_groups1zi, dud_yearzi, dud_patchzi, jrepstind,
              jrepstind_rownames, jsizeindzi, jsizeind_rownames_zi, false,
              grp2o(i), grp1(i), patchnumber, yearnumber, 4, 13, exp_tol,
              theta_tol, ipm_cdf, matrixformat, fecmod, repentry(i), negfec,
              stage2n(i), nostages, 0, 0.0);
              
            if (fl3(i) == 0) {
              out(i, 2) = 1.0 - out(i, 2);
            }
          } else {
            if (fl3(i) == 0) {
              out(i, 2) = 1.0 - jrepstcoefs(0);
            } else if (fl3(i) == 1) {
              out(i, 2) = jrepstcoefs(0);
            } else {
              out(i, 2) = 0.0;
            }
          }
        } else {
          out(i, 1) = 1.0 - out(i, 1);
          out(i, 2) = 1.0;
          out(i, 3) = 1.0;
          out(i, 4) = 1.0;
          out(i, 5) = 1.0;
          out(i, 6) = 1.0;
        }
        
        survtransmat(k) = out(i, 0) * out(i, 1) * out(i, 2) * out(i, 3) *
          out(i, 4) * out(i, 5) * out(i, 6);
      }
    } else if (ovgivent(i) != -1) {
      // All other transitions
      
      survtransmat(k) = ovgivent(i);
    }
    
    // This next block calculates fecundity
    if (indata2n(i) == 1 && fecl > 0) {
      if (fl2o(i) == 1 && ovgivenf(i) == -1) {
        
        fectransmat(k) = preouterator(fecproxy, feccoefs, rand_index, dev_terms,
          svsigmas, vital_year, vital_patch, chosen_r2inda, chosen_r1inda,
          chosen_r2indb, chosen_r1indb, chosen_r2indc, chosen_r1indc,
          statusterms, fecgroups2, fecgroups1, fecgroups2zi, fecgroups1zi,
          fecyearzi, fecpatchzi, fecind, fecind_rownames, fecindzi,
          fecind_rownames_zi, feczero, grp2o(i), grp1(i), patchnumber, yearnumber,
          fecdist, 7, exp_tol, theta_tol, ipm_cdf, matrixformat, fecmod,
          repentry(i), negfec, stage2n(i), nostages, fectrunc, fecsigma);
        
      } else if (ovgivenf(i) != -1 ) {
        fectransmat(k) = ovgivenf(i);
      }
    } else if (ovgivenf(i) != -1 ) {
      fectransmat(k) = ovgivenf(i);
    }
  }
  
  double ov_mult {0};
  if (replacementst > 0) {
    for (int i = 0; i < replacementst; i++) {
      
      repindex = replacetvec(i); // AllStages index
      properindex = aliveandequal(repindex);
      arma::uvec rightindex = find(index321 == ovestt(repindex));
      
      if (rightindex.n_elem > 0) {
        proxyindex = aliveandequal(rightindex(0));
        
        ov_mult = ovsurvmult(i);
        if (ov_mult < 0.0) ov_mult = 1.0;
        
        survtransmat(properindex) = survtransmat(proxyindex) * ov_mult;
      }
    }
  }
  
  if (replacementsf > 0) {
    for (int i = 0; i < replacementsf; i++) {
      
      repindex = replacefvec(i); // AllStages index
      properindex = aliveandequal(repindex);
      arma::uvec rightindex = find(index321 == ovestf(repindex));
      
      if (rightindex.n_elem > 0) {
        proxyindex = aliveandequal(rightindex(0));
        
        ov_mult = ovfecmult(i);
        if (ov_mult < 0.0) ov_mult = 1.0;
        
        fectransmat(properindex) = fectransmat(proxyindex) * ov_mult;
      }
    }
  }
  
  arma::mat amatrix = survtransmat + fectransmat;
  
  return List::create(Named("A") = amatrix, _["U"] = survtransmat,
    _["F"] = fectransmat, _["out"] = out);
}

//' Create Historically Structured Version of ahMPM
//' 
//' Function \code{.thefifthhousemate()} takes an ahistorical MPM as input, and
//' uses the \code{allstages} index to create a historically structured version
//' of it.
//' 
//' @param mpm The original ahMPM, supplied as a \code{lefkoMat} object.
//' @param allstages The index dataframe developed by
//' \code{\link{.simplepizzle}()}.
//' @param stageframe The ahistorical stageframe supplied by
//' \code{\link{.simplepizzle}()}.
//' @param format Integer indicating whether historical matrices should be in
//' (1) Ehrlen or (2) deVries format.
//' 
//' @return This will return a list of lists. The first list is composed of all
//' new \code{A} matrices. The second list is composed of all new \code{U}
//' matrices. The third list is composed of all new \code{F} matrices.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.thefifthhousemate)]]
Rcpp::List thefifthhousemate (List mpm, DataFrame allstages,
  DataFrame stageframe, int format) {
  Rcpp::List old_Umats = mpm["U"];
  Rcpp::List old_Fmats = mpm["F"];
  
  Rcpp::IntegerVector stageid = stageframe["stage_id"];
  int nostages = stageid.length();
  int nocols = nostages * nostages;
  
  if (format == 2) nocols = (nostages -1) * nostages;
  
  Rcpp::IntegerVector old_index = allstages["index21"];
  Rcpp::IntegerVector new_index = allstages["index321"];
  
  int num_mats = old_Umats.length();
  int index_elems = new_index.length();
  
  Rcpp::List new_Umats(num_mats);
  Rcpp::List new_Fmats(num_mats);
  Rcpp::List new_Amats(num_mats);
  
  arma::mat new_U(nocols, nocols, fill::zeros);
  arma::mat new_F(nocols, nocols, fill::zeros);
  arma::mat new_A(nocols, nocols, fill::zeros);
  arma::mat old_U(nostages, nostages, fill::zeros);
  arma::mat old_F(nostages, nostages, fill::zeros);
  
  for (int i = 0; i < num_mats; i++) {
    new_U.zeros();
    new_F.zeros();
    new_A.zeros();
    old_U.zeros();
    old_F.zeros();
    
    old_U = as<arma::mat>(old_Umats(i));
    old_F = as<arma::mat>(old_Fmats(i));
    
    for (int j = 0; j < index_elems; j++) {
      new_U(new_index(j)) = old_U(old_index(j));
      new_F(new_index(j)) = old_F(old_index(j));
    }
    
    new_A = new_U + new_F;
    new_Umats(i) = new_U;
    new_Fmats(i) = new_F;
    new_Amats(i) = new_A;
  }
  
  Rcpp::List output = List::create(Named("A") = new_Amats, _["U"] = new_Umats,
    _["F"] = new_Fmats);
  return output;
}

//' Creates Matrices of Year and Patch Terms in Leslie Models
//' 
//' Function \code{.revelations_leslie()} creates a matrix holding either the
//' year or patch coefficients from Leslie vital rate models. This reduces
//' memory load in function \code{\link{.motherbalowski}()}.
//' 
//' @param survproxy The proxy vital rate model covering survival from the main
//' matrix estimator function.
//' @param fecproxy The proxy vital rate model covering fecundity from the main
//' matrix estimator function.
//' 
//' @return A matrix with 2 columns corresponding to the number of vital rates
//' and number of columns equal to the number of year or patches.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.revelations_leslie)]]
arma::mat revelations_leslie(List survproxy, List fecproxy, int mat_switch) {
  
  arma::mat final_mat;
  
  if (mat_switch == 1) {
    Rcpp::DataFrame survyear_df(survproxy["years"]);
    Rcpp::DataFrame fecyear_df(fecproxy["years"]);

    arma::vec survyear = survyear_df[0];
    arma::vec fecyear = fecyear_df[0];

    int matrows = survyear.n_elem;
    
    arma::mat year_mat(matrows, 2, fill::zeros);
    year_mat.col(0) = survyear;
    year_mat.col(1) = fecyear;
    
    final_mat = year_mat;
    
  } else if (mat_switch == 2) {
    
    Rcpp::DataFrame survpatch_df(survproxy["patches"]);
    Rcpp::DataFrame fecpatch_df(fecproxy["patches"]);

    arma::vec survpatch = survpatch_df[0];
    arma::vec fecpatch = fecpatch_df[0];

    int matrows = survpatch.n_elem;
    
    arma::mat patch_mat(matrows, 2, fill::zeros);
    patch_mat.col(0) = survpatch;
    patch_mat.col(1) = fecpatch;

    final_mat = patch_mat;
  }
  
  return final_mat;
}

//' Create Index of Element Numbers for Random Individual Covariate Terms in
//' Leslie Models
//' 
//' Function \code{.foi_index_leslie()} creates a matrix indexing the end points
//' of each random individual covariate in the utilized vectors. Used in
//' function \code{\link{.motherbalowski}()}.
//' 
//' @param surv_proxy Adult survival model proxy.
//' @param fec_proxy Adult fecundity model proxy.
//' 
//' @return An integer matrix with 6 rows and 3 columns. The columns contain the
//' number of elements in each random individual covariate term, with the row
//' order being: 1) cov a t2, 2) cov a t1, 3) cov b t2, 4) cov b t1,
//' 5) cov c t2, and 6) cov c t1.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.foi_index_leslie)]]
arma::imat foi_index_leslie(List surv_proxy, List fec_proxy) {
  
  arma::ivec surv_fc = foi_counter(surv_proxy, false);
  arma::ivec fec_fc = foi_counter(fec_proxy, false);
  arma::ivec fec_fc_zi = foi_counter(fec_proxy, true);
  
  arma::imat final_mat(6, 3, fill::zeros);
  
  for (int i = 0; i < 6; i++) {
    final_mat(i, 0) = surv_fc(i);
    final_mat(i, 1) = fec_fc(i);
    final_mat(i, 2) = fec_fc_zi(i);
  }
  
  return final_mat;
}

//' Estimate All Elements of Function-based Population Projection Matrix
//' 
//' Function \code{.motherbalowski()} swiftly calculates matrix elements in
//' function-based Leslie population projection matrices. Used in
//' \code{\link{fleslie}()}.
//' 
//' @param ppy A data frame with one row, showing the population, patch, and
//' year.
//' @param ageframe The modified stageframe used in matrix calculations.
//' @param survproxy List of coefficients estimated in model of survival.
//' @param fecproxy List of coefficients estimated in model of fecundity.
//' @param f2_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{a} at
//' each time \emph{t} to be used in analysis.
//' @param f1_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{a} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{b} at
//' each time \emph{t} to be used in analysis.
//' @param f1_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{b} at
//' each time \emph{t}-1 to be used in analysis.
//' @param f2_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{c} at
//' each time \emph{t} to be used in analysis.
//' @param f1_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual covariate \code{c} at
//' each time \emph{t}-1 to be used in analysis.
//' @param r2_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{a} at each time \emph{t} to be used in analysis.
//' @param r1_inda A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{a} at each time \emph{t}-1 to be used in analysis.
//' @param r2_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{b} at each time \emph{t} to be used in analysis.
//' @param r1_indb A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{b} at each time \emph{t}-1 to be used in analysis.
//' @param r2_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{c} at each time \emph{t} to be used in analysis.
//' @param r1_indc A numeric vector of length equal to the number of years,
//' holding values equal to the mean value of individual random covariate
//' \code{c} at each time \emph{t}-1 to be used in analysis.
//' @param surv_dev A numeric value indicating the deviation to the linear
//' model of survival input by the user.
//' @param fec_dev A numeric value indicating the deviation to the linear
//' model of fecundity input by the user.
//' @param dens A numeric value equal to the density to be used in calculations.
//' @param fecmod A scalar multiplier for fecundity.
//' @param finalage The final age to be included in Leslie MPM estimation.
//' @param fecdist Designates whether fecundity is Gamma (3), Gaussian (2),
//' Poisson (0), or negative binomial (1) distributed.
//' @param negfec A logical value denoting whether to change negative estimated
//' fecundity to 0.
//' @param exp_tol A numeric value indicating the maximum limit for the
//' \code{exp()} function to be used in vital rate calculations. Defaults to
//' \code{700.0}.
//' @param theta_tol A numeric value indicating a maximum value for theta in
//' negative binomial probability density estimation. Defaults to
//' \code{100000000.0}.
//' 
//' @return A list of 3 matrices, including the main MPM (A), the survival-
//' transition matrix (U), and a fecundity matrix (F).
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.motherbalowski)]]
List motherbalowski(DataFrame ppy, DataFrame ageframe, List survproxy,
  List fecproxy, NumericVector f2_inda, NumericVector f1_inda,
  NumericVector f2_indb, NumericVector f1_indb, NumericVector f2_indc,
  NumericVector f1_indc, StringVector r2_inda, StringVector r1_inda,
  StringVector r2_indb, StringVector r1_indb, StringVector r2_indc,
  StringVector r1_indc, double surv_dev, double fec_dev, double dens,
  double fecmod, unsigned int finalage, int fecdist, bool negfec,
  double exp_tol = 700.0, double theta_tol = 100000000.0) {
  
  // Determines the size of the matrix
  StringVector sf_agenames = ageframe["stage"];
  IntegerVector sf_minage = ageframe["min_age"];
  IntegerVector sf_maxage = ageframe["max_age"];
  IntegerVector sf_repstatus = ageframe["repstatus"];
  int noages = sf_minage.length();
  
  bool cont = false;
  if (sf_maxage(noages - 1) == NA_INTEGER) {
    cont = true;
  }
  
  // Proxy model imports and settings
  bool feczero = false;

  NumericVector survcoefs = survproxy["coefficients"];
  NumericVector feccoefs = fecproxy["coefficients"];

  int survl = survcoefs.length();
  int fecl = feccoefs.length();

  double fecsigma = fecproxy["sigma"];

  if (NumericVector::is_na(fecsigma)) {
    if (fecdist == 1) {
      fecsigma = 1.0;
    } else {
      fecsigma = 0.0;
    }
  }

  arma::mat vital_year = revelations_leslie(survproxy, fecproxy, 1);
  arma::mat vital_patch = revelations_leslie(survproxy, fecproxy, 2);
  
  Rcpp::DataFrame fecyearzi_df(fecproxy["zeroyear"]);
  Rcpp::DataFrame fecpatchzi_df(fecproxy["zeropatch"]);
  Rcpp::DataFrame survgroups2_df(survproxy["groups2"]);
  Rcpp::DataFrame fecgroups2_df(fecproxy["groups2"]);
  Rcpp::DataFrame survgroups1_df(survproxy["groups1"]);
  Rcpp::DataFrame fecgroups1_df(fecproxy["groups1"]);
  Rcpp::DataFrame fecgroups2zi_df(fecproxy["zerogroups2"]);
  Rcpp::DataFrame fecgroups1zi_df(fecproxy["zerogroups1"]);
  
  arma::vec fecyearzi = fecyearzi_df[0];
  arma::vec fecpatchzi = fecpatchzi_df[0];
  arma::vec survgroups2 = survgroups2_df[0];
  arma::vec fecgroups2 = fecgroups2_df[0];
  arma::vec survgroups1 = survgroups1_df[0];
  arma::vec fecgroups1 = fecgroups1_df[0];
  arma::vec fecgroups2zi = fecgroups2zi_df[0];
  arma::vec fecgroups1zi = fecgroups1zi_df[0];
  
  if (fecyearzi.n_elem > 1 || fecyearzi(0) != 0) feczero = true;
  
  arma::vec survind = flightoficarus(survproxy);
  arma::vec fecind = flightoficarus(fecproxy);
  arma::vec fecindzi = zero_flightoficarus(fecproxy);
  
  arma::imat rand_index = foi_index_leslie(survproxy, fecproxy);
  
  StringVector survind_rownames = bootson(survproxy);
  StringVector fecind_rownames = bootson(fecproxy);
  StringVector fecind_rownames_zi = zero_bootson(fecproxy);
  
  // listofyears import and settings
  arma::uvec years = ppy["yearorder"];
  int yearnumber = years(0) - 1;
  
  arma::uvec patches = ppy["patchorder"];
  int patchnumber = patches(0) - 1;
  
  // Determination of choices of fixed and random individual covariates
  double inda1 = f1_inda(yearnumber);
  double indb1 = f1_indb(yearnumber);
  double indc1 = f1_indc(yearnumber);
  double inda2 = f2_inda(yearnumber);
  double indb2 = f2_indb(yearnumber);
  double indc2 = f2_indc(yearnumber);
  
  String chosen_r2inda = r2_inda(yearnumber);
  String chosen_r1inda = r1_inda(yearnumber);
  String chosen_r2indb = r2_indb(yearnumber);
  String chosen_r1indb = r1_indb(yearnumber);
  String chosen_r2indc = r2_indc(yearnumber);
  String chosen_r1indc = r1_indc(yearnumber);
  
  // The output matrices
  arma::mat survtransmat(noages, noages, fill::zeros);
  arma::mat fectransmat(noages, noages, fill::zeros);
  
  // The following loop runs through each age, and so runs through
  // each estimable element in the matrix
  for(int i = 0; i < noages; i++) {
    // Adult survival transitions
    
    double preout {0.0};
    
    if (survl > 1) {
      
      double chosen_randcova2 {0.0};
      if (chosen_r2inda != "none") {
        for (int indcount = 0; indcount < rand_index(0, 0); indcount++) {
          if (chosen_r2inda == survind_rownames(indcount)) {
            chosen_randcova2 = survind(indcount);
          }
        }
      }
      double chosen_randcova1 {0.0};
      if (chosen_r1inda != "none") {
        int delectable_sum = rand_index(0, 0);
        for (int indcount = 0; indcount < rand_index(1, 0); indcount++) {
          if (chosen_r1inda == survind_rownames(indcount + delectable_sum)) {
            chosen_randcova1 = survind(indcount + delectable_sum);
          }
        }
      }
      double chosen_randcovb2 {0.0};
      if (chosen_r2indb != "none") {
        int delectable_sum = rand_index(0, 0) + rand_index(1, 0);
        for (int indcount = 0; indcount < rand_index(2, 0); indcount++) {
          if (chosen_r2indb == survind_rownames(indcount + delectable_sum)) {
            chosen_randcovb2 = survind(indcount + delectable_sum);
          }
        }
      }
      double chosen_randcovb1 {0.0};
      if (chosen_r1indb != "none") {
        int delectable_sum = rand_index(0, 0) + rand_index(1, 0) + rand_index(2, 0);
        for (int indcount = 0; indcount < rand_index(3, 0); indcount++) {
          if (chosen_r1indb == survind_rownames(indcount + delectable_sum)) {
            chosen_randcovb1 = survind(indcount + delectable_sum);
          }
        }
      }
      double chosen_randcovc2 {0.0};
      if (chosen_r2indc != "none") {
        int delectable_sum = rand_index(0, 0) + rand_index(1, 0) + rand_index(2, 0) +
          rand_index(3, 0);
        for (int indcount = 0; indcount < rand_index(4, 0); indcount++) {
          if (chosen_r2indc == survind_rownames(indcount + delectable_sum)) {
            chosen_randcovc2 = survind(indcount + delectable_sum);
          }
        }
      }
      double chosen_randcovc1 {0.0};
      if (chosen_r1indc != "none") {
        int delectable_sum = rand_index(0, 0) + rand_index(1, 0) + rand_index(2, 0) +
          rand_index(3, 0) + rand_index(4, 0);
        for (int indcount = 0; indcount < rand_index(5, 0); indcount++) {
          if (chosen_r1indc == survind_rownames(indcount + delectable_sum)) {
            chosen_randcovc1 = survind(indcount + delectable_sum);
          }
        }
      }
      
      double mainsum = rimeotam(survcoefs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, static_cast<double>(i), inda1, inda2, indb1, indb2, indc1, indc2,
        dens, false);
      
      preout = (mainsum + chosen_randcova2 + chosen_randcova1 +
        chosen_randcovb2 + chosen_randcovb1 + chosen_randcovc2 +
        chosen_randcovc1 + survgroups2(0) + survgroups1(0) + 
        vital_patch(patchnumber, 0) + vital_year(yearnumber, 0) + surv_dev);

      if (preout > exp_tol) preout = exp_tol; // This catches numbers too high to be dealt with properly
      if (i < (noages - 1)) {
        survtransmat(i+1, i) = exp(preout) / (1.0 + exp(preout));
      } else {
        if (cont) {
          survtransmat(i, i) = exp(preout) / (1.0 + exp(preout));
        }
      }
    } else {
      if (i < (noages - 1)) {
        survtransmat(i+1, i) = survcoefs(0);
      } else {
        survtransmat(i, i) = survcoefs(0);
      }
    }
    
    // This next block calculates fecundity
    if (fecl > 0) {
      if (sf_repstatus(i) == 1) {
        
        double chosen_randcova2 {0.0};
        if (chosen_r2inda != "none") {
          for (int indcount = 0; indcount < rand_index(0, 6); indcount++) {
            if (chosen_r2inda == fecind_rownames(indcount)) {
              chosen_randcova2 = fecind(indcount);
            }
          }
        }
        double chosen_randcova1 {0.0};
        if (chosen_r1inda != "none") {
          int delectable_sum = rand_index(0, 6);
          for (int indcount = 0; indcount < rand_index(1, 6); indcount++) {
            if (chosen_r1inda == fecind_rownames(indcount + delectable_sum)) {
              chosen_randcova1 = fecind(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovb2 {0.0};
        if (chosen_r2indb != "none") {
          int delectable_sum = rand_index(0, 6) + rand_index(1, 6);
          for (int indcount = 0; indcount < rand_index(2, 6); indcount++) {
            if (chosen_r2indb == fecind_rownames(indcount + delectable_sum)) {
              chosen_randcovb2 = fecind(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovb1 {0.0};
        if (chosen_r1indb != "none") {
          int delectable_sum = rand_index(0, 6) + rand_index(1, 6) + rand_index(2, 6);
          for (int indcount = 0; indcount < rand_index(3, 6); indcount++) {
            if (chosen_r1indb == fecind_rownames(indcount + delectable_sum)) {
              chosen_randcovb1 = fecind(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovc2 {0.0};
        if (chosen_r2indc != "none") {
          int delectable_sum = rand_index(0, 6) + rand_index(1, 6) + rand_index(2, 6) +
            rand_index(3, 6);
          for (int indcount = 0; indcount < rand_index(4, 6); indcount++) {
            if (chosen_r2indc == fecind_rownames(indcount + delectable_sum)) {
              chosen_randcovc2 = fecind(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovc1 {0.0};
        if (chosen_r1indc != "none") {
          int delectable_sum = rand_index(0, 6) + rand_index(1, 6) + rand_index(2, 6) +
            rand_index(3, 6) + rand_index(4, 6);
          for (int indcount = 0; indcount < rand_index(5, 6); indcount++) {
            if (chosen_r1indc == fecind_rownames(indcount + delectable_sum)) {
              chosen_randcovc1 = fecind(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcova2zi {0.0};
        if (chosen_r2inda != "none") {
          for (int indcount = 0; indcount < rand_index(0, 16); indcount++) {
            if (chosen_r2inda == fecind_rownames_zi(indcount)) {
              chosen_randcova2zi = fecindzi(indcount);
            }
          }
        }
        double chosen_randcova1zi {0.0};
        if (chosen_r1inda != "none") {
          int delectable_sum = rand_index(0, 16);
          for (int indcount = 0; indcount < rand_index(1, 16); indcount++) {
            if (chosen_r1inda == fecind_rownames_zi(indcount + delectable_sum)) {
              chosen_randcova1zi = fecindzi(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovb2zi {0.0};
        if (chosen_r2indb != "none") {
          int delectable_sum = rand_index(0, 16) + rand_index(1, 16);
          for (int indcount = 0; indcount < rand_index(2, 16); indcount++) {
            if (chosen_r2indb == fecind_rownames_zi(indcount + delectable_sum)) {
              chosen_randcovb2zi = fecindzi(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovb1zi {0.0};
        if (chosen_r1indb != "none") {
          int delectable_sum = rand_index(0, 16) + rand_index(1, 16) + rand_index(2, 16);
          for (int indcount = 0; indcount < rand_index(3, 16); indcount++) {
            if (chosen_r1indb == fecind_rownames_zi(indcount + delectable_sum)) {
              chosen_randcovb1zi = fecindzi(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovc2zi {0.0};
        if (chosen_r2indc != "none") {
          int delectable_sum = rand_index(0, 16) + rand_index(1, 16) + rand_index(2, 16) +
            rand_index(3, 16);
          for (int indcount = 0; indcount < rand_index(4, 16); indcount++) {
            if (chosen_r2indc == fecind_rownames_zi(indcount + delectable_sum)) {
              chosen_randcovc2zi = fecindzi(indcount + delectable_sum);
            }
          }
        }
        double chosen_randcovc1zi {0.0};
        if (chosen_r1indc != "none") {
          int delectable_sum = rand_index(0, 16) + rand_index(1, 16) + rand_index(2, 16) +
            rand_index(3, 16) + rand_index(4, 16);
          for (int indcount = 0; indcount < rand_index(5, 16); indcount++) {
            if (chosen_r1indc == fecind_rownames_zi(indcount + delectable_sum)) {
              chosen_randcovc1zi = fecindzi(indcount + delectable_sum);
            }
          }
        }
            
        double preoutx {0.0};
        
        if (fecdist < 4 && fecl > 1) {
          if (feczero) {
            
            double mainsum = rimeotam(feccoefs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, static_cast<double>(i), inda1, inda2, indb1, indb2, indc1,
              indc2, dens, true);
            
            preoutx = (mainsum + chosen_randcova2zi + chosen_randcova1zi +
              chosen_randcovb2zi + chosen_randcovb1zi + chosen_randcovc2zi +
              chosen_randcovc1zi + fecgroups2zi(0) + fecgroups1zi(0) + 
              fecpatchzi(patchnumber) + fecyearzi(yearnumber) + fec_dev);
            
          } else {
            
            double mainsum = rimeotam(feccoefs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, static_cast<double>(i), inda1, inda2, indb1, indb2, indc1,
              indc2, dens, false);
            
            preoutx = (mainsum + chosen_randcova2 + chosen_randcova1 +
              chosen_randcovb2 + chosen_randcovb1 + chosen_randcovc2 +
              chosen_randcovc1 + fecgroups2(0) + fecgroups1(0) + 
              vital_patch(patchnumber, 1) + vital_year(yearnumber, 1) + fec_dev);
          }
          
          if (fecdist == 0 || fecdist == 1) {
            // Poisson and negative binomial fecundity
            
            if (feczero) {
              if (preoutx > exp_tol) preoutx = exp_tol;
              
              fectransmat(0, i) = (exp(preoutx) / (1.0 + exp(preoutx))) * fecmod;
            } else {
              if (preoutx > exp_tol) preoutx = exp_tol;
              
              fectransmat(0, i) = exp(preoutx) * fecmod;
            }
          } else if (fecdist == 2) {
            // Gaussian fecundity
            fectransmat(0, i) = preoutx * fecmod;
            
            if (negfec && fectransmat(0, i) < 0.0) {
              fectransmat(0, i) = 0.0;
            }
          } else if (fecdist == 3) {
            // Gamma fecundity
            fectransmat(0, i) = (1.0 / preoutx) * fecmod;
          }
          
        } else if (fecl > 1) {
          
          // All others with estimated models
          fectransmat(0, i) = 0.0;
        } else {
          fectransmat(0, i) = feccoefs(0);
        }
      }
    }
  }
  
  arma::mat amatrix = survtransmat + fectransmat;
  
  return List::create(Named("A") = amatrix, _["U"] = survtransmat,
    _["F"] = fectransmat);
}

