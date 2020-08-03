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

//' Estimate All Elements in Raw Historical Matrix
//' 
//' This function swiftly calculates matrix transitions in raw historical matrices,
//' and serves as the core workhorse function behind \code{\link{rlefko3}()}.
//' 
//' @param sge9l A 7 column matrix containing information on fecundity potential,
//' reproductive status, presence in main dataset, supplied given rates (survival
//' and fecundity), and estimated proxy rates (survival and fecundity), 
//' respectively, for all combinations of stage pairs at times \emph{t}+1 and \emph{t}, 
//' and times \emph{t} and \emph{t}-1.
//' @param sge3 A 2 column matrix containing reproductive status and fecundity 
//' potential of stage pairs.
//' @param maindata A 2 column matrix of raw data denoting status as alive and
//' offspring produced.
//' @param sge93index Full matrix index vector denoting each element with respect to
//' From stage pair (times \emph{t} and \emph{t}-1) and To stage pair (times \emph{t}+1 and
//' \emph{t}).
//' @param sge92index Vector denoting From stage pair index in full matrix.
//' @param sge32index Vector indicating overall stage pair index (baseline From
//' pair for survival estimation).
//' @param sge33 Vector with stage at time \emph{t}+1 in \code{sge32index}.
//' @param sge32 Vector with stage at time \emph{t} in \code{sge32index}.
//' @param data3221 Vector of stage-pair combination indices in raw dataset.
//' @param data21 Vector of stage-pair indices in raw dataset, corresponding to
//' \code{sge32index} and used in survival estimation.
//' @param nostages The number of ahistorical stages.
//' 
//' @return Matrix composed of the survival-transitions (U) matrix as the first 
//' column and the fecundity matrix (F) as the second column.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat specialpatrolgroup(arma::mat sge9l, arma::mat sge3, arma::mat maindata, arma::uvec sge93index,
                             arma::uvec sge92index, arma::uvec sge32index, arma::uvec sge33, arma::uvec sge32,
                             arma::uvec data3221, arma::uvec data21, int nostages) {
  
  int n = maindata.n_rows;
  int no2stages = sge3.n_rows;
  int noelems = sge9l.n_rows;
  
  arma::mat probsrates(noelems, 4); // 1st col = # indivs (3 trans), 2nd col = total indivs for pair stage,
  // 3rd col = total indivs alive for pair stage, 4th col = total fec for pair stage,
  probsrates.zeros();
  
  arma::mat stage2fec(no2stages, 3); //1st col = # indivs total, 2nd col = no indivs alive, 3rd col = sum fec
  stage2fec.zeros();
  
  for (int i = 0; i < n; i++) { // This main loop counts individuals going through transitions and sums their
    // fecundities, and then adds that info to the 3-trans and 2-trans tables
    
    probsrates((data3221(i)), 0) = probsrates((data3221(i)), 0) + 1; // Yields sum of all individuals with particuar transition
    
    stage2fec((data21(i)), 0) = stage2fec((data21(i)), 0) + 1; // Yields sum of all individuals with particuar transition
    if (maindata(i, 0) > 0) {
      stage2fec((data21(i)), 1) = stage2fec((data21(i)), 1) + 1;
    }
    
    stage2fec((data21(i)), 2) = stage2fec((data21(i)), 2) + maindata(i,1);
    
  }
  
  for (int i = 0; i < no2stages; i++) {
    unsigned int foradding = ((sge33(i) - 1) * nostages) + ((sge33(i) - 1) * nostages * nostages) + 
      ((sge32(i) - 1) * nostages * nostages * nostages);
    
    for (int j = 0; j < nostages; j++) {
      unsigned int entry = foradding + j;
      
      probsrates(entry, 1) = stage2fec(i, 0);
      probsrates(entry, 2) = stage2fec(i, 1);
      probsrates(entry, 3) = stage2fec(i, 2);
    }
  }
  
  arma::mat out(noelems, 2); // Main output matrix
  
  out.col(0) = probsrates.col(0) / probsrates.col(1); // Survival
  out.col(1) = sge9l.col(0) % sge9l.col(1) % probsrates.col(3) / probsrates.col(1); // Fecundity
  
  out.elem(find_nonfinite(out)).zeros();
  
  // Now we will correct transitions and rates for given stuff
  
  double ovgiventsum = arma::accu(sge9l.col(3)) / noelems;
  double ovgivenfsum = arma::accu(sge9l.col(4)) / noelems;
  
  double ovesttsum = arma::accu(sge9l.col(5)) / noelems;
  double ovestfsum = arma::accu(sge9l.col(6)) / noelems;
  
  if ((ovgiventsum + ovgivenfsum + ovesttsum + ovestfsum) > -4) {
    
    arma::uvec ovgiventind = find(sge9l.col(3) != -1);
    arma::uvec ovgivenfind = find(sge9l.col(4) != -1);
    arma::uvec ovesttind = find(sge9l.col(5) != -1);
    arma::uvec ovestfind = find(sge9l.col(6) != -1);
    
    int ovgtn = ovgiventind.n_elem;
    int ovgfn = ovgivenfind.n_elem;
    int ovestn = ovesttind.n_elem;
    int ovesfn = ovestfind.n_elem;
    
    if (ovgtn > 0) {
      for (int i = 0; i < ovgtn; i++) {
        out(ovgiventind(i), 0) = sge9l(ovgiventind(i), 3);
      }
    }
    
    if (ovgfn > 0) {
      for (int i = 0; i < ovgfn; i++) {
        out(ovgivenfind(i), 0) = sge9l(ovgivenfind(i), 4);
      }
    }
    
    if (ovestn > 0) {
      for (int i = 0; i < ovestn; i++) {
        unsigned int replacement = sge9l(ovesttind(i), 5);
        
        out(ovesttind(i), 0) = out(replacement, 0);
      }
    }
    
    if (ovesfn > 0) {
      for (int i = 0; i < ovesfn; i++) {
        unsigned int replacement = sge9l(ovestfind(i), 6);
        
        out(ovestfind(i), 0) = out(replacement, 0);
      }
    }
  }
  
  return out;
}

//' Estimate All Elements in Raw Ahistorical Population Projection Matrices
//' 
//' This function swiftly calculates matrix transitions in raw ahistorical matrices,
//' and serves as the core workhorse function behind \code{\link{rlefko2}()}.
//' 
//' @param sge9l A 7 column matrix containing information on fecundity potential,
//' reproductive status, presence in main dataset, supplied given rates (survival
//' and fecundity), and estimated proxy rates (survival and fecundity), 
//' respectively, for all combinations of stages at times \emph{t}+1 and \emph{t}.
//' @param sge3 A 2 column matrix containing reproductive status and fecundity 
//' potential of ahistorical stages.
//' @param maindata A 2 column matrix of raw data denoting status as alive and
//' offspring produced.
//' @param sge93index Full matrix index vector denoting each element with respect to
//' stage pairs (times \emph{t}+1 and \emph{t}).
//' @param sge92index Vector denoting From stage in full matrix.
//' @param sge32index Vector indicating overall From stage (baseline From for
//' survival estimation).
//' @param sge33 Vector with stage at time \emph{t} in \code{sge92index}.
//' @param sge32 Vector with stage at time \emph{t} in \code{sge32index}.
//' @param data3221 Vector of stage-pair indices in raw dataset.
//' @param data21 Vector of stage indices in raw dataset, corresponding to
//' \code{sge32index} and used in survival estimation.
//' @param nostages The number of ahistorical stages.
//' 
//' @return Matrix composed of the survival-transitions (U) matrix as the first 
//' column and the fecundity matrix (F) as the second column.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat normalpatrolgroup(arma::mat sge9l, arma::mat sge3, arma::mat maindata, arma::uvec sge93index,
                            arma::uvec sge92index, arma::uvec sge32index, arma::uvec sge32,
                            arma::uvec data3221, arma::uvec data21, int nostages) {
  
  int n = maindata.n_rows;
  int no2stages = sge3.n_rows;
  int noelems = sge9l.n_rows;
  
  arma::mat probsrates(noelems, 4); // 1st col = # indivs (3 trans), 2nd col = total indivs for pair stage,
  // 3rd col = total indivs alive for pair stage, 4th col = total fec for pair stage,
  probsrates.zeros();
  
  arma::mat stage2fec(no2stages, 3); //1st col = # indivs total, 2nd col = no indivs alive, 3rd col = sum fec
  stage2fec.zeros();
  
  for (int i = 0; i < n; i++) { // This main loop counts individuals going through transitions and sums their
    // fecundities, and then adds that info to the 3-trans and 2-trans tables
    
    probsrates((data3221(i)), 0) = probsrates((data3221(i)), 0) + 1; // Yields sum of all individuals with particuar transition
    
    stage2fec((data21(i)), 0) = stage2fec((data21(i)), 0) + 1; // Yields sum of all individuals with particuar transition
    if (maindata(i, 0) > 0) {
      stage2fec((data21(i)), 1) = stage2fec((data21(i)), 1) + 1;
    }
    
    stage2fec((data21(i)), 2) = stage2fec((data21(i)), 2) + maindata(i,1);
    
  }
  
  for (int i = 0; i < no2stages; i++) {
    unsigned int foradding = ((sge32(i) - 1) * nostages);
    
    for (int j = 0; j < nostages; j++) {
      unsigned int entry = foradding + j;
      
      probsrates(entry, 1) = stage2fec(i, 0);
      probsrates(entry, 2) = stage2fec(i, 1);
      probsrates(entry, 3) = stage2fec(i, 2);
    }
  }
  
  arma::mat out(noelems, 2); // Main output matrix
  
  out.col(0) = probsrates.col(0) / probsrates.col(1);
  out.col(1) = sge9l.col(0) % sge9l.col(1) % probsrates.col(3) / probsrates.col(1);
  
  out.elem(find_nonfinite(out)).zeros();
  
  // Now we will correct transitions and rates for given stuff
  
  double ovgiventsum = arma::accu(sge9l.col(3)) / noelems;
  double ovgivenfsum = arma::accu(sge9l.col(4)) / noelems;
  
  double ovesttsum = arma::accu(sge9l.col(5)) / noelems;
  double ovestfsum = arma::accu(sge9l.col(6)) / noelems;
  
  if ((ovgiventsum + ovgivenfsum + ovesttsum + ovestfsum) > -4) {
    
    arma::uvec ovgiventind = find(sge9l.col(3) != -1);
    arma::uvec ovgivenfind = find(sge9l.col(4) != -1);
    arma::uvec ovesttind = find(sge9l.col(5) != -1);
    arma::uvec ovestfind = find(sge9l.col(6) != -1);
    
    int ovgtn = ovgiventind.n_elem;
    int ovgfn = ovgivenfind.n_elem;
    int ovestn = ovesttind.n_elem;
    int ovesfn = ovestfind.n_elem;
    
    if (ovgtn > 0) {
      for (int i = 0; i < ovgtn; i++) {
        out(ovgiventind(i), 0) = sge9l(ovgiventind(i), 3);
      }
    }
    
    if (ovgfn > 0) {
      for (int i = 0; i < ovgfn; i++) {
        out(ovgivenfind(i), 0) = sge9l(ovgivenfind(i), 4);
      }
    }
    
    if (ovestn > 0) {
      for (int i = 0; i < ovestn; i++) {
        unsigned int replacement = sge9l(ovesttind(i), 5);
        
        out(ovesttind(i), 0) = out(replacement, 0);
      }
    }
    
    if (ovesfn > 0) {
      for (int i = 0; i < ovesfn; i++) {
        unsigned int replacement = sge9l(ovestfind(i), 6);
        
        out(ovestfind(i), 0) = out(replacement, 0);
      }
    }
  }
  
  return out;
}

//' Creates Index of Estimable Population Projection Matrix Elements
//' 
//' This function identifies which matrix elements in a projection matrix are
//' logically possible to estimate. In both historical and ahistorical matrices,
//' this will effectively remove the Dead stage from the final matrix output. In
//' addition, in historical matrices, this will identify matrix elements 
//' corresponding to From and To stage-pair combinations in which stages at time \emph{t}
//' are equal. Used in \code{\link{rlefko3}()}, \code{\link{rlefko2}()}, \code{\link{flefko3}()}, 
//' and \code{\link{flefko2}()}.
//' 
//' @param mainindex Should only living elements (0) be identified, or should stages
//' equal at time \emph{t} be identified (1)?
//' @param allstages A 4 column matrix identifying stage at time \emph{t}+1, stage at
//' time {t} (if historical, then in To stage-pair), stage at time \emph{t} (if
//' historical, then in From stage pair), and stage at time \emph{t}-1 if historical
//' or stage at time \emph{t} if ahistorical.
//' 
//' @return Vector of estimable matrix element indices.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::uvec hoffmannofstuttgart(int mainindex, arma::mat allstages) {
  if (mainindex == 0) {
    
    arma::vec stages3 = allstages.col(0);
    double death = stages3.max();
    
    arma::uvec aliverows3 = find(stages3 != death);
    arma::uvec aliverows2n = find(allstages.col(1) != death);
    arma::uvec aliverows2o = find(allstages.col(2) != death);
    arma::uvec aliverows1 = find(allstages.col(3) != death);
    
    arma::uvec aliverows32o = intersect(aliverows3, aliverows2o);
    arma::uvec aliverows2n1 = intersect(aliverows2n, aliverows1);
    
    arma::uvec aliverows = intersect(aliverows32o, aliverows2n1);
    
    return aliverows;
    
  } else if (mainindex == 1) {
    
    arma::vec stageratios = allstages.col(1) / allstages.col(2);
    
    arma::uvec equal2 = find(stageratios == 1);

    return equal2;
    
  } else {
    
    arma::vec stages3 = allstages.col(0);
    arma::vec stages2 = allstages.col(1);
    double death = stages3.max();
    
    int endbit = stages3.n_elem;
    
    arma::uvec aliverows(endbit);
    
    aliverows.zeros();
    
    for (int i = 0; i < endbit; i++) {
      if (stages3(i) != death && stages2(i) != death) aliverows(i) = 1;
    }
    
    return aliverows;
  }
}

//' Estimate All Elements in Function-based Population Projection Matrices
//' 
//' This function swiftly calculates matrix elements in function-based population
//' projection matrices. Used in \code{\link{flefko3}()} and \code{\link{flefko2}()}.
//' 
//' @param survcoefs Vector of coefficients estimated in model of survival.
//' @param obscoefs Vector of coefficients estimated in model of observation.
//' @param sizecoefs Vector of coefficients estimated in model of size.
//' @param repstcoefs Vector of coefficients estimated in model of reproductive 
//' status.
//' @param feccoefs Vector of coefficients estimated in model of fecundity.
//' @param jsurvcoefs Vector of coefficients estimated in model of juvenile
//' survival.
//' @param jobscoefs Vector of coefficients estimated in model of juvenile
//' observation.
//' @param jsizecoefs Vector of coefficients estimated in model of juvenile size.
//' @param jrepstcoefs Vector of coefficients estimated in model of juvenile
//' reproductive status.
//' @param stage3 Vector of stage index at time \emph{t}+1 across projection matrix.
//' @param stage2n Vector of stage index at time \emph{t} across projection matrix.
//' @param stage1 Vector of stage index at time \emph{t}-1 across projection matrix.
//' @param sz3 Vector of size at time \emph{t}+1 across projection matrix.
//' @param sz2n Vector of size at time \emph{t} across projection matrix.
//' @param sz1 Vector of size at time \emph{t}-1 across projection matrix.
//' @param fl3 Vector of reproductive status at time \emph{t}+1 across projection 
//' matrix.
//' @param fl2n Vector of reproductive status at time \emph{t} across projection 
//' matrix.
//' @param fl1 Vector of reproductive status at time \emph{t}-1 across projection 
//' matrix.
//' @param ob3 Vector of observation status at time \emph{t}+1 across projection 
//' matrix.
//' @param ob2n Vector of observation status at time \emph{t} across projection 
//' matrix.
//' @param ob1 Vector of observation status at time \emph{t}-1 across projection 
//' matrix.
//' @param immat3 Vector of immaturity status at time \emph{t}+1 across projection 
//' matrix.
//' @param immat2n Vector of immaturity status at time \emph{t} across projection 
//' matrix.
//' @param indata2 Vector indicating whether stage time time \emph{t} in the To
//' stage pair is estimable.
//' @param indata Vector indicating whether the combination of stages at times
//' \emph{t}+1 and \emph{t} in the To stage, and time \emph{t} in the From stage,
//' are estimable.
//' @param aliveandequal Vector indicating whether the transition indicated by the
//' respective From and To stage-pairs is estimable.
//' @param repentry Vector identifying whether stage at time {t}+1 is a reproductive
//' entry stage, and whether stage at time \emph{t} is potentially reproductive.
//' @param ovgivent Vector of given survival-transition probabilities for 
//' replacement.
//' @param ovgivenf Vector of given fecundity rates for replacement.
//' @param binwidth3 Width of size class bins associated with stage in time \emph{t}+1.
//' @param numofsizes4 Number of elements in main matrix.
//' @param fecmod A scalar multiplier for fecundity.
//' @param summedvars Summed variance-covariance terms in Poisson size distribution.
//' @param sigma Standard deviation of Gaussian size distribution.
//' @param jsummedvars Summed variance-covariance terms in Poisson juvenile size
//' distribution.
//' @param jsigma Standard deviation of Gaussian juvenile size distribution.
//' @param maxsize The maximum size to be used in element estimation.
//' @param sizedist Designates whether size is Gaussian (2), Poisson (0), or
//' negative binomial (1) distributed.
//' @param fecdist Designates whether fecundity is Gaussian (2), Poisson (0), or
//' negative binomial (1) distributed.
//' 
//' @return A 6 column matrix showing survival probability, observation probability,
//' reproduction probability, size transition probability, overall 
//' survival-transition probability, and fecundity rate, for each element of the
//' final population projection matrix.
//' 
//' @keywords internal
//' @noRd
//[[Rcpp::export]]
arma::mat jerzeibalowski(arma::vec survcoefs, arma::vec obscoefs, arma::vec sizecoefs, 
                         arma::vec repstcoefs, arma::vec feccoefs, arma::vec jsurvcoefs,
                         arma::vec jobscoefs, arma::vec jsizecoefs, arma::vec jrepstcoefs,
                         arma::vec stage3, arma::vec stage2n, arma::vec stage1, 
                         arma::vec sz3, arma::vec sz2n, arma::vec sz1, arma::vec fl3, 
                         arma::vec fl2n, arma::vec fl1, arma::vec ob3, arma::vec ob2n, 
                         arma::vec ob1, arma::vec immat3, arma::vec immat2n, arma::vec indata2,
                         arma::vec indata, arma::uvec aliveandequal, arma::vec repentry,
                         arma::vec ovgivent, arma::vec ovgivenf, arma::vec binwidth3,
                         unsigned long numofsizes4, double fecmod, double summedvars, double sigma,
                         double jsummedvars, double jsigma, double maxsize, int sizedist, int fecdist) {
  
  // Lines 1-3 introduce coefficient vectors from modeling, and include random coefficients where 
  // needed. Lines 3-8 of the input above are typically long vectors composed of input sizes and the 
  // such for calculations. Lines 8-11 introduces variables used in size and fecundity calculations.
  // The first three provide needed input for the Poisson, negative binomial, and Gaussian 
  // distributions, and include variables are designed to tell C++ what distributions are used to  
  // model size and fecundity. If either equals 0, then the Poisson is used. If either equals 1, 
  // then the negative binomial is used. If 2, then the Gaussian. Please note the numofsizes4 
  // variable, which should be the number of stages NOT INCLUDING DEATH raised to the power of 4.
  
  // In the code below, the function decides on an appropriate single loop routine based on which coefficient
  // vectors have length greater than 1. Only a single routine will run because of the series of else statements
  
  int n = stage3.n_elem;
  
  int survl = survcoefs.n_elem;
  int obsl = obscoefs.n_elem;
  int repstl = repstcoefs.n_elem;
  int fecl = feccoefs.n_elem;
  int jsurvl = jsurvcoefs.n_elem;
  int jobsl = jobscoefs.n_elem;
  int jrepstl = jrepstcoefs.n_elem;
  
  arma::mat out(numofsizes4, 6);  // Should be a 0 matrix with numofsizes4 rows and 6 columns: 
                                  // 0 is surv, 1 is obs, 2 is repst, 3 is size, 
                                  // 4 is survival-transition, and 5 is fec
  
  out.zeros();
  
  for(int i = 0; i < n; i++) {
    unsigned int k = aliveandequal(i);
    
    if (ovgivent(i) == -1 && indata(i) == 1) {
      if (immat2n(i) == 0 && immat3(i) == 0) {
        // Adult survival transitions
        
        arma::vec preout(4);
        
        if (survl > 1) {
          preout(0) = (survcoefs(0) + survcoefs(16) + (survcoefs(4) * sz2n(i)) + (survcoefs(3) * sz1(i)) +
            (survcoefs(2) * fl2n(i)) + (survcoefs(1) * fl1(i)) + (survcoefs(6) * sz2n(i) * sz1(i)) +
            (survcoefs(5) * fl2n(i) * fl1(i)) + (survcoefs(7) * sz1(i) * fl1(i)) +
            (survcoefs(8) * sz2n(i) * fl2n(i)) + (survcoefs(10) * sz1(i) * fl2n(i)) +
            (survcoefs(9) * sz2n(i) * fl1(i)) + survcoefs(17) + survcoefs(18));
          
          out(k, 0) = exp(preout(0)) / (1 + exp(preout(0)));
          
        } else {
          out(k, 0) = survcoefs(0);
        }
        
        if (obsl > 1) {
          preout(1) = (obscoefs(0) + obscoefs(16) + (obscoefs(4) * sz2n(i)) + (obscoefs(3) * sz1(i)) + 
            (obscoefs(2) * fl2n(i)) + (obscoefs(1) * fl1(i)) + (obscoefs(6) * sz2n(i) * sz1(i)) +
            (obscoefs(5) * fl2n(i) * fl1(i)) + (obscoefs(7) * sz1(i) * fl1(i)) +
            (obscoefs(8) * sz2n(i) * fl2n(i)) + (obscoefs(10) * sz1(i) * fl2n(i)) + 
            (obscoefs(9) * sz2n(i) * fl1(i)) + obscoefs(17) + obscoefs(18));
          
          out(k, 1) = exp(preout(1)) / (1 + exp(preout(1)));
          
        } else {
          out(k, 1) = obscoefs(0);
        }
        
        if (ob3(i) == 1) {
          if (sizedist == 0) {
            // Poisson size distribution
            
            double sizefac = sz3(i) * tgamma(sz3(i));
            double lambda = exp(sizecoefs(0) + sizecoefs(16) + (sizecoefs(4) * sz2n(i)) + (sizecoefs(3) * sz1(i)) + 
                                (sizecoefs(2) * fl2n(i)) + (sizecoefs(1) * fl1(i)) + (sizecoefs(6) * sz2n(i) * sz1(i)) +
                                (sizecoefs(5) * fl2n(i) * fl1(i)) + (sizecoefs(7) * sz1(i) * fl1(i)) +
                                (sizecoefs(8) * sz2n(i) * fl2n(i)) + (sizecoefs(10) * sz1(i) * fl2n(i)) + 
                                (sizecoefs(9) * sz2n(i) * fl1(i)) + sizecoefs(17) + sizecoefs(18) + (summedvars / 2));
            
            out(k, 3) = ((pow(lambda, sz3(i)) * exp(-1 * lambda)) / sizefac);
          } else if (sizedist == 1) {
            // Negative binomial size distribution
            
            double mu = exp(sizecoefs(0) + sizecoefs(16) + (sizecoefs(4) * sz2n(i)) + (sizecoefs(3) * sz1(i)) + 
                            (sizecoefs(2) * fl2n(i)) + (sizecoefs(1) * fl1(i)) + (sizecoefs(6) * sz2n(i) * sz1(i)) +
                            (sizecoefs(5) * fl2n(i) * fl1(i)) + (sizecoefs(7) * sz1(i) * fl1(i)) +
                            (sizecoefs(8) * sz2n(i) * fl2n(i)) + (sizecoefs(10) * sz1(i) * fl2n(i)) + 
                            (sizecoefs(9) * sz2n(i) * fl1(i)) + sizecoefs(17) + sizecoefs(18));
            double y = sz3(i);
            double r = maxsize - y;
            
            double p = mu / (mu + r);
            
            double binoc = tgamma(maxsize) / ((y*tgamma(y)) * tgamma(r));
            double midterm = pow((1-p), r);
            double rightterm = pow(p, y);
            
            out(k, 3) = (binoc * midterm * rightterm);
            
          } else if (sizedist == 2) {
            // Gaussian size distribution
            
            double sigma2 = sigma * sigma;
            preout(3) = (sizecoefs(0) + sizecoefs(16) + (sizecoefs(4) * sz2n(i)) + (sizecoefs(3) * sz1(i)) + 
              (sizecoefs(2) * fl2n(i)) + (sizecoefs(1) * fl1(i)) + (sizecoefs(6) * sz2n(i) * sz1(i)) +
              (sizecoefs(5) * fl2n(i) * fl1(i)) + (sizecoefs(7) * sz1(i) * fl1(i)) +
              (sizecoefs(8) * sz2n(i) * fl2n(i)) + (sizecoefs(10) * sz1(i) * fl2n(i)) + 
              (sizecoefs(9) * sz2n(i) * fl1(i)) + sizecoefs(17) + sizecoefs(18));
            
            out(k, 3) = (exp(-1 * (pow((sz3(i) - preout(3)), 2) / (2*sigma2))) / ((pow((2 * M_PI), 0.5)) * sigma)) * binwidth3(i);
            
          } else {
            out(k, 3) = sizedist - 3;
          }
          
          if (repstl > 1) {
            
            preout(2) = (repstcoefs(0) + repstcoefs(16) + (repstcoefs(4) * sz2n(i)) + (repstcoefs(3) * sz1(i)) + 
              (repstcoefs(2) * fl2n(i)) + (repstcoefs(1) * fl1(i)) + (repstcoefs(6) * sz2n(i) * sz1(i)) +
              (repstcoefs(5) * fl2n(i) * fl1(i)) + (repstcoefs(7) * sz1(i) * fl1(i)) +
              (repstcoefs(8) * sz2n(i) * fl2n(i)) + (repstcoefs(10) * sz1(i) * fl2n(i)) + 
              (repstcoefs(9) * sz2n(i) * fl1(i)) + repstcoefs(17) + repstcoefs(18));
            
            out(k, 2) = exp(preout(2)) / (1 + exp(preout(2)));
            
            if (fl3(i) == 0) {
              out(k, 2) = 1 - out(k, 2);
            }
            
          } else {
            if (fl3(i) == 0) {
              out(k, 2) = 1 - repstcoefs(0);
            } else if (fl3(i) == 1) {
              out(k, 2) = repstcoefs(0);
            } else {
              out(k, 2) = 0;
            }
          }
          
        } else {
          out(k, 1) = 1 - out(k, 1);
          out(k, 3) = 1;
          out(k, 2) = 1;
        }
        
        out(k, 4) = out(k, 0) * out(k, 1) * out(k, 2) * out(k, 3);
        
        
      } else if (immat2n(i) == 1 && immat3(i) == 0 && jsurvl > 0) {
        // Juvenile to adult transitions
        
        arma::vec preout(4);
        
        if (jsurvl > 1) {
          preout(0) = (jsurvcoefs(0) + jsurvcoefs(1) + jsurvcoefs(2) + jsurvcoefs(3));
          
          out(k, 0) = exp(preout(0)) / (1 + exp(preout(0)));
        } else {
          out(k, 0) = jsurvcoefs(0);
        }
        
        if (jobsl > 1) {
          preout(1) = (jobscoefs(0) + jobscoefs(1) + jobscoefs(2) + jobscoefs(3));
          
          out(k, 1) = exp(preout(1)) / (1 + exp(preout(1)));
          
        } else {
          out(k, 1) = jobscoefs(0);
        }
        
        if (ob3(i) == 1) {
          if (sizedist == 0) {
            // Poisson size distribution
            
            double sizefac = sz3(i) * tgamma(sz3(i));
            double lambda = exp(jsizecoefs(0) + jsizecoefs(1) + (jsummedvars / 2));
            
            out(k, 3) = ((pow(lambda, sz3(i)) * exp(-1 * lambda)) / sizefac);
            
          } else if (sizedist == 1) {
            // Negative binomial size distribution
            
            double mu = exp(jsizecoefs(0));
            double y = sz3(i);
            double r = maxsize - y;
            
            double p = mu / (mu + r);
            
            double binoc = tgamma(maxsize) / ((y*tgamma(y)) * tgamma(r));
            double midterm = pow((1-p), r);
            double rightterm = pow(p, y);
            
            out(k, 3) = (binoc * midterm * rightterm);
            
          } else if (sizedist == 2) {
            // Gaussian size distribution
            
            double sigma2 = jsigma * jsigma;
            preout(3) = (jsizecoefs(0) + jsizecoefs(1) + jsizecoefs(2) + jsizecoefs(3));
            
            out(k, 3) = (exp(-1 * (pow((sz3(i) - preout(3)), 2) / (2*sigma2))) / ((pow((2 * M_PI), 0.5)) * jsigma)) * binwidth3(i);
            
          } else {
            out(k, 3) = sizedist - 3;
          }
          
          if (jrepstl > 1) {
            
            preout(2) = (jrepstcoefs(0) + jrepstcoefs(1) + jrepstcoefs(2) + jrepstcoefs(3));
            
            out(k, 2) = exp(preout(2)) / (1 + exp(preout(2)));
            
            if (fl3(i) == 0) {
              out(k, 2) = 1 - out(k, 2);
            }
            
          } else {
            if (fl3(i) == 0) {
              out(k, 2) = 1 - jrepstcoefs(0);
            } else if (fl3(i) == 1) {
              out(k, 2) = jrepstcoefs(0);
            } else {
              out(k, 2) = 0;
            }
          }
          
        } else {
          out(k, 1) = 1 - out(k, 1);
          out(k, 3) = 1;
          out(k, 2) = 1;
        }
        
        out(k, 4) = out(k, 0) * out(k, 1) * out(k, 2) * out(k, 3);
        
      }
    } else if (ovgivent(i) != -1) {
      // All other transitions
      
      out(k, 4) = ovgivent(i);
    }
    
    if (indata2(i) == 1 && fecl > 0) {
      // Fecundity
      if (fl2n(i) == 1 && ovgivenf(i) == -1) {
        
        double preoutx;
        
        if (fecdist < 3) {
          preoutx = (feccoefs(0) + feccoefs(16) + (feccoefs(4) * sz2n(i)) + (feccoefs(3) * sz1(i)) + 
            (feccoefs(2) * fl2n(i)) + (feccoefs(1) * fl1(i)) + (feccoefs(6) * sz2n(i) * sz1(i)) +
            (feccoefs(5) * fl2n(i) * fl1(i)) + (feccoefs(7) * sz1(i) * fl1(i)) +
            (feccoefs(8) * sz2n(i) * fl2n(i)) + (feccoefs(10) * sz1(i) * fl2n(i)) + 
            (feccoefs(9) * sz2n(i) * fl1(i)) + feccoefs(17) + feccoefs(18));
          
          if (fecdist != 2) {
            // Poisson and negative binomial fecundity
            
            out(k, 5) = exp(preoutx) * fecmod * repentry(i);
          } else {
            // Gaussian fecundity
            out(k, 5) = preoutx * fecmod * repentry(i);
          }
          
        } else {
          // All others
          
          out(k, 5) = fecdist - 3;
        }
        
      } else if (ovgivenf(i) != -1 ) {
        out(k, 5) = ovgivenf(i);
      }
    } else if (ovgivenf(i) != -1 ) {
      out(k, 5) = ovgivenf(i);
    }
  }
  
  return out;
}

//' Estimates Mean Population Projection Matrix, Using Summed U and F Matrices
//' 
//' This function estimates the mean population projection matrices, treating the
//' mean as arithmetic across space but either arithmetic or geometric across time.
//' It differs from \code{\link{turbogeodiesel}()} in that it estimates the \code{A} matrix
//' as a sum of the associated \code{U} and \code{F} matrices. Used to power the
//' \code{\link{lmean}()} function.
//' 
//' @param geom Should the mean across time be geometric (1) or arithmetic (0)?
//' @param sparse Should 0s be ignored when some matrices include non-zero entries
//' in common elements?
//' @param numofpops Number of populations to be analyzed.
//' @param numofpatches Number of patches to be analyzed, where this number should
//' include a patch total across all populations.
//' @param numofyears Number of time steps to be analyzed.
//' @param loy2c Matrix denoting the population, patch, and time step designation
//' of each matrix.
//' @param Umats A matrix with all U matrices turned into columns.
//' @param Fmats A matrix with all F matrices turned into columns.
//' 
//' @return A matrix with 3n columns, where n is the sum of the number of patches and
//' populations. Each pop/patch has its own set of three columns denoting survival,
//' fecundity, and the overall sum of the previous two columns.
//' 
//' @keywords internal
//' @noRd
//[[Rcpp::export]]
arma::mat geodiesel(int geom, int sparse, int numofpops, int numofpatches, int numofyears, 
                    arma::mat loy2c, arma::mat Umats, arma::mat Fmats) {
  
  
  int elements = Umats.n_rows;
  int matrices = Umats.n_cols - 1;
  
  arma::vec indexmatU = Umats.col(matrices);
  arma::vec indexmatF = Fmats.col(matrices);
  
  arma::mat out(elements, ((numofpatches*3) + (numofpops * 3)));
  out.zeros();
  
  
  if (geom == 1) {
    if (sparse == 0) {
      // Geometric without sparse
      
      for (int i = 0; i < elements; i++) {
        
        if (indexmatU(i) != 0 || indexmatF(i) != 0) {
          
          for (int j = 0; j < matrices; j++) {
            //This section creates the sum of logs for each non-zero positive element
            
            if (Umats(i, j) != 0 && out(i, loy2c(j, 1)) != 1000) {
              out(i, loy2c(j, 1)) = out(i, loy2c(j, 1)) + (log(Umats(i, j)) / loy2c(j, 4));
            } else {
              out(i, loy2c(j, 1)) = 1000;
            }
            
            if (Fmats(i, j) != 0 && out(i, (numofpatches + loy2c(j, 1))) != 1000) {
              out(i, (numofpatches + loy2c(j, 1))) = out(i, (numofpatches + loy2c(j, 1))) + (log(Fmats(i, j)) / loy2c(j, 4));
            } else {
              out(i, (numofpatches + loy2c(j, 1))) = 1000;
            }
          }
          
          for (int k = 0; k < numofpatches; k++) {
            // This section takes e^logtotal for the U and F elements and creates the A elements,
            // and then creates the grand means
            
            // U elements
            if (out(i, k) != 1000) {
              out(i, k) = exp(out(i, k));
            } else {
              out(i, k) = 0;
            }
            
            // F elements
            if (out(i, (k + numofpatches)) != 1000) {
              out(i, (k + numofpatches)) = exp(out(i, (k + numofpatches)));
            } else {
              out(i, (k + numofpatches)) = 0;
            }
            
            // A elements
            out(i, (k + 2*numofpatches)) = out(i, k) + out(i, (k + numofpatches));
            
            // Grand mean elements
            out(i, (3*numofpatches)) = out(i, (3*numofpatches)) + (out(i, k) / numofpatches);
            out(i, (3*numofpatches + 1)) = out(i, (3*numofpatches + 1)) + (out(i, (numofpatches + k)) / numofpatches);
            out(i, (3*numofpatches + 2)) = out(i, (3*numofpatches)) + out(i, (3*numofpatches + 1));
          }
        }
      }
    } else {
      // Geometric with sparse
      arma::uvec totalvec(numofpatches*2);
      
      for (int i = 0; i < elements; i++) {
        
        if (indexmatU(i) != 0 || indexmatF(i) != 0) {
          
          for (int l = 0; l < numofpatches; l++) {
            //This section creates a vector to keep track of the dividend in mean calculations
            
            int properindex = numofyears * l;
            totalvec(l) = loy2c(properindex, 4);
            totalvec(l + numofpatches) = loy2c(properindex, 4);
          }
          
          for (int j = 0; j < matrices; j++) {
            //This section creates the sum of logs for each non-zero positive element
            
            if (Umats(i, j) != 0) {
              out(i, loy2c(j, 1)) = out(i, loy2c(j, 1)) + log(Umats(i, j));
            } else {
              totalvec(loy2c(j, 1)) = totalvec(loy2c(j, 1)) - 1;
            }
            
            if (Fmats(i, j) != 0) {
              out(i, (numofpatches + loy2c(j, 1))) = out(i, (numofpatches + loy2c(j, 1))) + log(Fmats(i, j));
            } else {
              totalvec(numofpatches + loy2c(j, 1)) = totalvec(numofpatches + loy2c(j, 1)) - 1;
            }
          }
          
          for (int k = 0; k < numofpatches; k++) {
            // This section takes e^logtotal for the U and F elements and creates the A elements,
            // and then creates the grand means
            
            // U elements
            if (totalvec(k) != 0) {
              out(i, k) = exp((out(i, k) / totalvec(k)));
            } else {
              out(i, k) = 0;
            }
            
            // F elements
            if (totalvec(k + numofpatches) != 0) {
              out(i, (k + numofpatches)) = exp((out(i, (k + numofpatches)) / totalvec(k + numofpatches)));
            } else {
              out(i, (k + numofpatches)) = 0;
            }
            
            // A elements
            out(i, (k + 2*numofpatches)) = out(i, k) + out(i, (k + numofpatches));
            
            // Grand mean elements
            out(i, (3*numofpatches)) = out(i, (3*numofpatches)) + (out(i, k) / numofpatches);
            out(i, (3*numofpatches + 1)) = out(i, (3*numofpatches + 1)) + (out(i, (numofpatches + k)) / numofpatches);
            out(i, (3*numofpatches + 2)) = out(i, (3*numofpatches)) + out(i, (3*numofpatches + 1));
          }
        }
      }
    }
  } else {
    if (sparse == 0) {
      // Arithmetic without sparse
      
      for (int i = 0; i < elements; i++) {
        
        if (indexmatU(i) != 0 || indexmatF(i) != 0) {
          
          for (int j = 0; j < matrices; j++) {
            //This section creates the means for U and F elements
            
            out(i, loy2c(j, 1)) = out(i, loy2c(j, 1)) + ((Umats(i, j)) / loy2c(j, 4));
            
            out(i, (numofpatches + loy2c(j, 1))) = out(i, (numofpatches + loy2c(j, 1))) + ((Fmats(i, j)) / loy2c(j, 4));
          }
          
          for (int k = 0; k < numofpatches; k++) {
            // This section creates the A elements,
            // and then creates the grand means
            
            // A elements
            out(i, (k + 2*numofpatches)) = out(i, k) + out(i, (k + numofpatches));
            
            // Grand mean elements
            out(i, (3*numofpatches)) = out(i, (3*numofpatches)) + (out(i, k) / numofpatches);
            out(i, (3*numofpatches + 1)) = out(i, (3*numofpatches + 1)) + (out(i, (numofpatches + k)) / numofpatches);
            out(i, (3*numofpatches + 2)) = out(i, (3*numofpatches)) + out(i, (3*numofpatches + 1));
          }
        }
      }
    } else {
      // Arithmetic with sparse
      
      arma::uvec totalvec(numofpatches*2);
      
      for (int i = 0; i < elements; i++) {
        
        if (indexmatU(i) != 0 || indexmatF(i) != 0) {
          
          for (int l = 0; l < numofpatches; l++) {
            //This section creates a vector to keep track of the dividend in mean calculations
            
            int properindex = numofyears * l;
            totalvec(l) = loy2c(properindex, 4);
            totalvec(l + numofpatches) = loy2c(properindex, 4);
          }
          
          for (int j = 0; j < matrices; j++) {
            //This section creates the sum of logs for each non-zero positive element
            
            if (Umats(i, j) != 0) {
              out(i, loy2c(j, 1)) = out(i, loy2c(j, 1)) + (Umats(i, j));
            } else {
              totalvec(loy2c(j, 1)) = totalvec(loy2c(j, 1)) - 1;
            }
            
            if (Fmats(i, j) != 0) {
              out(i, (numofpatches + loy2c(j, 1))) = out(i, (numofpatches + loy2c(j, 1))) + (Fmats(i, j));
            } else {
              totalvec(numofpatches + loy2c(j, 1)) = totalvec(numofpatches + loy2c(j, 1)) - 1;
            }
          }
          
          for (int k = 0; k < numofpatches; k++) {
            // This section takes sum for the U and F elements and creates the A elements,
            // and then creates the grand means
            
            // U elements
            if (totalvec(k) != 0) {
              out(i, k) = (out(i, k) / totalvec(k));
            } else {
              out(i, k) = 0;
            }
            
            // F elements
            if (totalvec(k + numofpatches) != 0) {
              out(i, (k + numofpatches)) = (out(i, (k + numofpatches)) / totalvec(k + numofpatches));
            } else {
              out(i, (k + numofpatches)) = 0;
            }
            
            // A elements
            out(i, (k + 2*numofpatches)) = out(i, k) + out(i, (k + numofpatches));
            
            // Grand mean elements
            out(i, (3*numofpatches)) = out(i, (3*numofpatches)) + (out(i, k) / numofpatches);
            out(i, (3*numofpatches + 1)) = out(i, (3*numofpatches + 1)) + (out(i, (numofpatches + k)) / numofpatches);
            out(i, (3*numofpatches + 2)) = out(i, (3*numofpatches)) + out(i, (3*numofpatches + 1));
          }
        }
      }
    }
  }
  
  return out;
}

//' Estimates Mean Population Projection Matrix
//' 
//' This function estimates the mean population projection matrices, treating the
//' mean as arithmetic across space but either arithmetic or geometric across time.
//' It differs from \code{\link{geodiesel}()} in that it estimates the \code{A} matrix
//' indepenently of the associated \code{U} and \code{F} matrices. Used to power the
//' \code{\link{lmean}()} function.
//' 
//' @param geom Should the mean across time be geometric (1) or arithmetic (0)?
//' @param sparse Should 0s be ignored when some matrices include non-zero entries
//' in common elements?
//' @param numofpops Number of populations to be analyzed.
//' @param numofpatches Number of patches to be analyzed, where this number should
//' include a patch total across all populations.
//' @param numofyears Number of time steps to be analyzed.
//' @param loy2c Matrix denoting the population, patch, and time step designation
//' of each matrix.
//' @param Umats A matrix with all U matrices turned into columns.
//' @param Fmats A matrix with all F matrices turned into columns.
//' @param Amats A matrix with all A matrices turned into columns.
//' 
//' @return A matrix with 3n columns, where n is the sum of the number of patches and
//' populations. Each pop/patch has its own set of three columns denoting survival,
//' fecundity, and the overall projection.
//' 
//' @keywords internal
//' @noRd
//[[Rcpp::export]]
arma::mat turbogeodiesel(int geom, int sparse, int numofpops, int numofpatches, int numofyears, 
                         arma::mat loy2c, arma::mat Umats, arma::mat Fmats, arma::mat Amats) {
  
  
  int elements = Umats.n_rows;
  int matrices = Umats.n_cols - 1;
  
  arma::vec indexmatU = Umats.col(matrices);
  arma::vec indexmatF = Fmats.col(matrices);
  
  arma::mat out(elements, ((numofpatches*3) + (numofpops * 3)));
  out.zeros();
  
  
  if (geom == 1) {
    if (sparse == 0) {
      // Geometric without sparse
      
      for (int i = 0; i < elements; i++) {
        
        if (indexmatU(i) != 0 || indexmatF(i) != 0) {
          
          for (int j = 0; j < matrices; j++) {
            //This section creates the sum of logs for each non-zero positive element
            
            if (Umats(i, j) != 0 && out(i, loy2c(j, 1)) != 1000) {
              out(i, loy2c(j, 1)) = out(i, loy2c(j, 1)) + (log(Umats(i, j)) / loy2c(j, 4));
            } else {
              out(i, loy2c(j, 1)) = 1000;
            }
            
            if (Fmats(i, j) != 0 && out(i, (numofpatches + loy2c(j, 1))) != 1000) {
              out(i, (numofpatches + loy2c(j, 1))) = out(i, (numofpatches + loy2c(j, 1))) + (log(Fmats(i, j)) / loy2c(j, 4));
            } else {
              out(i, (numofpatches + loy2c(j, 1))) = 1000;
            }
            
            if (Amats(i, j) != 0 && out(i, (2 * numofpatches + loy2c(j, 1))) != 1000) {
              out(i, (2 * numofpatches + loy2c(j, 1))) = out(i, (2 * numofpatches + loy2c(j, 1))) + (log(Amats(i, j)) / loy2c(j, 4));
            } else {
              out(i, (2 * numofpatches + loy2c(j, 1))) = 1000;
            }
          }
          
          for (int k = 0; k < numofpatches; k++) {
            // This section takes e^logtotal for the U and F elements and creates the A elements,
            // and then creates the grand means
            
            // U elements
            if (out(i, k) != 1000) {
              out(i, k) = exp(out(i, k));
            } else {
              out(i, k) = 0;
            }
            
            // F elements
            if (out(i, (k + numofpatches)) != 1000) {
              out(i, (k + numofpatches)) = exp(out(i, (k + numofpatches)));
            } else {
              out(i, (k + numofpatches)) = 0;
            }
            
            // A elements
            if (out(i, (k + 2 * numofpatches)) != 1000) {
              out(i, (k + 2 * numofpatches)) = exp(out(i, (k + 2 * numofpatches)));
            } else {
              out(i, (k + 2 * numofpatches)) = 0;
            }
            
            // Grand mean elements
            out(i, (3*numofpatches)) = out(i, (3*numofpatches)) + (out(i, k) / numofpatches);
            out(i, (3*numofpatches + 1)) = out(i, (3*numofpatches + 1)) + (out(i, (numofpatches + k)) / numofpatches);
            out(i, (3*numofpatches + 2)) = out(i, (3*numofpatches)) + out(i, (3*numofpatches + 1));
          }
        }
      }
    } else {
      // Geometric with sparse
      arma::uvec totalvec(numofpatches*3);
      
      for (int i = 0; i < elements; i++) {
        
        if (indexmatU(i) != 0 || indexmatF(i) != 0) {
          
          for (int l = 0; l < numofpatches; l++) {
            //This section creates a vector to keep track of the dividend in mean calculations
            
            int properindex = numofyears * l;
            totalvec(l) = loy2c(properindex, 4);
            totalvec(l + numofpatches) = loy2c(properindex, 4);
            totalvec(l + 2 * numofpatches) = loy2c(properindex, 4);
          }
          
          for (int j = 0; j < matrices; j++) {
            //This section creates the sum of logs for each non-zero positive element
            
            if (Umats(i, j) != 0) {
              out(i, loy2c(j, 1)) = out(i, loy2c(j, 1)) + log(Umats(i, j));
            } else {
              totalvec(loy2c(j, 1)) = totalvec(loy2c(j, 1)) - 1;
            }
            
            if (Fmats(i, j) != 0) {
              out(i, (numofpatches + loy2c(j, 1))) = out(i, (numofpatches + loy2c(j, 1))) + log(Fmats(i, j));
            } else {
              totalvec(numofpatches + loy2c(j, 1)) = totalvec(numofpatches + loy2c(j, 1)) - 1;
            }
            
            if (Amats(i, j) != 0) {
              out(i, (2 * numofpatches + loy2c(j, 1))) = out(i, (2 * numofpatches + loy2c(j, 1))) + log(Amats(i, j));
            } else {
              totalvec(2 * numofpatches + loy2c(j, 1)) = totalvec(2 * numofpatches + loy2c(j, 1)) - 1;
            }
          }
          
          for (int k = 0; k < numofpatches; k++) {
            // This section takes e^logtotal for the U and F elements and creates the A elements,
            // and then creates the grand means
            
            // U elements
            if (totalvec(k) != 0) {
              out(i, k) = exp((out(i, k) / totalvec(k)));
            } else {
              out(i, k) = 0;
            }
            
            // F elements
            if (totalvec(k + numofpatches) != 0) {
              out(i, (k + numofpatches)) = exp((out(i, (k + numofpatches)) / totalvec(k + numofpatches)));
            } else {
              out(i, (k + numofpatches)) = 0;
            }
            
            // A elements
            if (totalvec(k + 2 * numofpatches) != 0) {
              out(i, (k + 2 * numofpatches)) = exp((out(i, (k + 2 * numofpatches)) / totalvec(k + 2 * numofpatches)));
            } else {
              out(i, (k + 2 * numofpatches)) = 0;
            }
            
            // Grand mean elements
            out(i, (3*numofpatches)) = out(i, (3*numofpatches)) + (out(i, k) / numofpatches);
            out(i, (3*numofpatches + 1)) = out(i, (3*numofpatches + 1)) + (out(i, (numofpatches + k)) / numofpatches);
            out(i, (3*numofpatches + 2)) = out(i, (3*numofpatches)) + out(i, (3*numofpatches + 1));
          }
        }
      }
    }
  } else {
    if (sparse == 0) {
      // Arithmetic without sparse
      
      for (int i = 0; i < elements; i++) {
        
        if (indexmatU(i) != 0 || indexmatF(i) != 0) {
          
          for (int j = 0; j < matrices; j++) {
            //This section creates the means for U, F, and A elements
            
            out(i, loy2c(j, 1)) = out(i, loy2c(j, 1)) + ((Umats(i, j)) / loy2c(j, 4));
            
            out(i, (numofpatches + loy2c(j, 1))) = out(i, (numofpatches + loy2c(j, 1))) + ((Fmats(i, j)) / loy2c(j, 4));
            
            out(i, (2 * numofpatches + loy2c(j, 1))) = out(i, (2 * numofpatches + loy2c(j, 1))) + ((Amats(i, j)) / loy2c(j, 4));
          }
          
          for (int k = 0; k < numofpatches; k++) {
            // This section creates the grand means
            
            out(i, (3*numofpatches)) = out(i, (3*numofpatches)) + (out(i, k) / numofpatches);
            out(i, (3*numofpatches + 1)) = out(i, (3*numofpatches + 1)) + (out(i, (numofpatches + k)) / numofpatches);
            out(i, (3*numofpatches + 2)) = out(i, (3*numofpatches)) + out(i, (3*numofpatches + 1));
          }
        }
      }
    } else {
      // Arithmetic with sparse
      
      arma::uvec totalvec(numofpatches*3);
      
      for (int i = 0; i < elements; i++) {
        
        if (indexmatU(i) != 0 || indexmatF(i) != 0) {
          
          for (int l = 0; l < numofpatches; l++) {
            //This section creates a vector to keep track of the dividend in mean calculations
            
            int properindex = numofyears * l;
            totalvec(l) = loy2c(properindex, 4);
            totalvec(l + numofpatches) = loy2c(properindex, 4);
            totalvec(l + 2 * numofpatches) = loy2c(properindex, 4);
          }
          
          for (int j = 0; j < matrices; j++) {
            //This section creates the sum of logs for each non-zero positive element
            
            if (Umats(i, j) != 0) {
              out(i, loy2c(j, 1)) = out(i, loy2c(j, 1)) + (Umats(i, j));
            } else {
              totalvec(loy2c(j, 1)) = totalvec(loy2c(j, 1)) - 1;
            }
            
            if (Fmats(i, j) != 0) {
              out(i, (numofpatches + loy2c(j, 1))) = out(i, (numofpatches + loy2c(j, 1))) + (Fmats(i, j));
            } else {
              totalvec(numofpatches + loy2c(j, 1)) = totalvec(numofpatches + loy2c(j, 1)) - 1;
            }
            
            if (Amats(i, j) != 0) {
              out(i, (2 * numofpatches + loy2c(j, 1))) = out(i, (2 * numofpatches + loy2c(j, 1))) + (Amats(i, j));
            } else {
              totalvec(2 * numofpatches + loy2c(j, 1)) = totalvec(2 * numofpatches + loy2c(j, 1)) - 1;
            }
          }
          
          for (int k = 0; k < numofpatches; k++) {
            // This section takes sum for the U and F elements and creates the A elements,
            // and then creates the grand means
            
            // U elements
            if (totalvec(k) != 0) {
              out(i, k) = (out(i, k) / totalvec(k));
            } else {
              out(i, k) = 0;
            }
            
            // F elements
            if (totalvec(k + numofpatches) != 0) {
              out(i, (k + numofpatches)) = (out(i, (k + numofpatches)) / totalvec(k + numofpatches));
            } else {
              out(i, (k + numofpatches)) = 0;
            }
            
            // A elements
            if (totalvec(k + 2 * numofpatches) != 0) {
              out(i, (k + 2 * numofpatches)) = (out(i, (k + 2 * numofpatches)) / totalvec(k + 2 * numofpatches));
            } else {
              out(i, (k + 2 * numofpatches)) = 0;
            }
            
            // Grand mean elements
            out(i, (3*numofpatches)) = out(i, (3*numofpatches)) + (out(i, k) / numofpatches);
            out(i, (3*numofpatches + 1)) = out(i, (3*numofpatches + 1)) + (out(i, (numofpatches + k)) / numofpatches);
            out(i, (3*numofpatches + 2)) = out(i, (3*numofpatches)) + out(i, (3*numofpatches + 1));
          }
        }
      }
    }
  }
  
  return out;
}


