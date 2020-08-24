#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

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
//' @param ppy A data frame with one row, showing the population, patch, and year.
//' @param AllStages A large data frame giving all required inputs for vital rate
//' estimation other than the vital rate model coefficients themselves. Contains
//' a row for each ultimate matrix element.
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
//' @param survdev Scalar value to be added to the y-intercept of the linear model
//' of survival probability.
//' @param obsdev Scalar value to be added to the y-intercept of the linear model
//' of observation probability.
//' @param sizedev Scalar value to be added to the y-intercept of the linear model
//' of size transition.
//' @param repstdev Scalar value to be added to the y-intercept of the linear model
//' of reproduction probability.
//' @param fecdev Scalar value to be added to the y-intercept of the linear model
//' of fecundity.
//' @param jsurvdev Scalar value to be added to the y-intercept of the linear model
//' of juvenile survival probability.
//' @param jobsdev Scalar value to be added to the y-intercept of the linear model
//' of juvenile observation probability.
//' @param jsizedev Scalar value to be added to the y-intercept of the linear model
//' of juvenile size transition.
//' @param jrepstdev Scalar value to be added to the y-intercept of the linear model
//' of juvenile reproduction probability.
//' @param numofsizes4 Number of elements in main matrix.
//' @param matrixdim Number of rows (and columns) in the final matrix.
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
//' @param negfec Logical value denoting whether to change negative estimated
//` fecundity to 0.
//' 
//' @return A list of 3 matrices, including the main MPM (A), the survival-transition
//' matrix (U), anf a fecundity matrix (F). With tweaking, can also produce a 4 column 
//' matrix showing survival probability, observation probability, reproduction
//' probability, and size transition probability, for each element of the final MPM.
//' 
//' @keywords internal
//' @noRd
//[[Rcpp::export]]
List jerzeibalowski(DataFrame ppy, DataFrame AllStages, List survproxy, List obsproxy,
                    List sizeproxy, List repstproxy, List fecproxy, List jsurvproxy,
                    List jobsproxy, List jsizeproxy, List jrepstproxy, double survdev,
                    double obsdev, double sizedev, double repstdev, double fecdev,
                    double jsurvdev, double jobsdev, double jsizedev, double jrepstdev,
                    unsigned long numofsizes4, unsigned long matrixdim, double fecmod, 
                    double summedvars, double sigma, double jsummedvars, double jsigma, 
                    double maxsize, int sizedist, int fecdist, bool negfec) {
  
  
  // The DataFrame introduces variables used in size and fecundity calculations. This DataFrame is
  // broken up into long vectors composed of input sizes and related variables for these calculations. 
  // The "model" Lists bring in the vital rate models, and include random coefficients
  // where needed. We also have a number of extra variables, that include such info as whether to use
  // the Poisson, negative binomial, and Gaussian for size and fecundity calculations. If either sizedist
  // or fecdist equals 0, then the Poisson is used. If either equals 1, then the negative binomial is 
  // used. If 2, then the Gaussian. Please note the numofsizes4 variable, which should be the number of
  // stages NOT INCLUDING DEATH raised to the power of 4.
  
  // In the code below, the function decides on an appropriate single loop routine based on which coefficient
  // vectors have length greater than 1. Only a single routine will run because of the series of else statements
  
  arma::vec survcoefs = survproxy["coefficients"];
  arma::vec obscoefs = obsproxy["coefficients"];
  arma::vec sizecoefs = sizeproxy["coefficients"];
  arma::vec repstcoefs = repstproxy["coefficients"];
  arma::vec feccoefs = fecproxy["coefficients"];
  arma::vec jsurvcoefs = jsurvproxy["coefficients"];
  arma::vec jobscoefs = jobsproxy["coefficients"];
  arma::vec jsizecoefs = jsizeproxy["coefficients"];
  arma::vec jrepstcoefs = jrepstproxy["coefficients"];
  
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
  
  arma::vec survyear = survproxy["years"];
  arma::vec obsyear = obsproxy["years"];
  arma::vec sizeyear = sizeproxy["years"];
  arma::vec repstyear = repstproxy["years"];
  arma::vec fecyear = fecproxy["years"];
  arma::vec jsurvyear = jsurvproxy["years"];
  arma::vec jobsyear = jobsproxy["years"];
  arma::vec jsizeyear = jsizeproxy["years"];
  arma::vec jrepstyear = jrepstproxy["years"];
  
  arma::vec stage3 = AllStages[0];
  Rcpp::NumericVector stage2n = AllStages[1];
  Rcpp::NumericVector stage1 = AllStages[3];
  Rcpp::NumericVector sz3 = AllStages[4];
  Rcpp::NumericVector sz2n = AllStages[5];
  Rcpp::NumericVector sz1 = AllStages[7];
  Rcpp::NumericVector ob3 = AllStages[8];
  Rcpp::NumericVector ob2n = AllStages[9];
  Rcpp::NumericVector ob1 = AllStages[11];
  Rcpp::NumericVector fl3 = AllStages[12];
  Rcpp::NumericVector fl2n = AllStages[13];
  Rcpp::NumericVector fl1 = AllStages[15];
  Rcpp::NumericVector mat3 = AllStages[16];
  Rcpp::NumericVector mat2n = AllStages[17];
  Rcpp::NumericVector mat1 = AllStages[19];
  Rcpp::NumericVector immat3 = AllStages[20];
  Rcpp::NumericVector immat2n = AllStages[21];
  Rcpp::NumericVector immat1 = AllStages[23];
  Rcpp::NumericVector indata2 = AllStages[26];
  Rcpp::NumericVector repentry = AllStages[24];
  Rcpp::NumericVector binwidth3 = AllStages[29];
  Rcpp::NumericVector minage3 = AllStages[30];
  Rcpp::NumericVector minage2 = AllStages[31];
  Rcpp::NumericVector maxage3 = AllStages[32];
  Rcpp::NumericVector maxage2 = AllStages[33];
  Rcpp::NumericVector actualage2 = AllStages[34];
  arma::uvec index321 = AllStages[35];
  arma::vec ovestt = AllStages[36];
  Rcpp::NumericVector ovgivent = AllStages[37];
  arma::vec ovestf = AllStages[38];
  Rcpp::NumericVector ovgivenf = AllStages[39];
  Rcpp::NumericVector indata = AllStages[40];
  Rcpp::NumericVector aliveandequal = AllStages[41];
  
  int n = stage3.n_elem;
  
  arma::uvec years = ppy[5];
  int yearnumber = years(0) - 1;
  
  arma::uvec patches = ppy[4];
  int patchnumber = patches(0) - 1;
  
  int survl = survcoefs.n_elem;
  int obsl = obscoefs.n_elem;
  int sizel = sizecoefs.n_elem;
  int repstl = repstcoefs.n_elem;
  int fecl = feccoefs.n_elem;
  int jsurvl = jsurvcoefs.n_elem;
  int jobsl = jobscoefs.n_elem;
  int jsizel = jsizecoefs.n_elem;
  int jrepstl = jrepstcoefs.n_elem;
  
  arma::uvec replacetvec = find(ovestt != -1);
  arma::uvec replacefvec = find(ovestf != -1);
  int replacementst = replacetvec.n_elem;
  int replacementsf = replacefvec.n_elem;
  int repindex {0};
  int properindex {0};
  int proxyindex {0};
  
  arma::mat out(numofsizes4, 4);  // 0 matrix with numofsizes4 rows & 4 columns: 0 surv, 1 obs, 2 repst, 3 size
  arma::mat survtransmat(matrixdim, matrixdim);
  arma::mat fectransmat(matrixdim, matrixdim);
  
  out.zeros();
  survtransmat.zeros();
  fectransmat.zeros();
  
  for(int i = 0; i < n; i++) {
    unsigned int k = aliveandequal(i);
    
    if (ovgivent(i) == -1 && indata(i) == 1) {
      if (mat2n(i) == 1 && mat3(i) == 1) {
        // Adult survival transitions
        
        arma::vec preout(4);
        
        if (survl > 1) {
          
          preout(0) = (survcoefs(0) + (survcoefs(4) * sz2n(i)) + (survcoefs(3) * sz1(i)) +
            (survcoefs(2) * fl2n(i)) + (survcoefs(1) * fl1(i)) + (survcoefs(6) * sz2n(i) * sz1(i)) +
            (survcoefs(5) * fl2n(i) * fl1(i)) + (survcoefs(7) * sz1(i) * fl1(i)) +
            (survcoefs(8) * sz2n(i) * fl2n(i)) + (survcoefs(10) * sz1(i) * fl2n(i)) +
            (survcoefs(9) * sz2n(i) * fl1(i)) + survpatch(patchnumber) + survyear(yearnumber) + survdev);
          
          out(k, 0) = exp(preout(0)) / (1 + exp(preout(0)));
          
        } else {
          out(k, 0) = survcoefs(0);
        }
        
        if (obsl > 1) {
          
          preout(1) = (obscoefs(0) + (obscoefs(4) * sz2n(i)) + (obscoefs(3) * sz1(i)) + 
            (obscoefs(2) * fl2n(i)) + (obscoefs(1) * fl1(i)) + (obscoefs(6) * sz2n(i) * sz1(i)) +
            (obscoefs(5) * fl2n(i) * fl1(i)) + (obscoefs(7) * sz1(i) * fl1(i)) +
            (obscoefs(8) * sz2n(i) * fl2n(i)) + (obscoefs(10) * sz1(i) * fl2n(i)) + 
            (obscoefs(9) * sz2n(i) * fl1(i)) + obspatch(patchnumber) + obsyear(yearnumber) + obsdev);
          
          out(k, 1) = exp(preout(1)) / (1 + exp(preout(1)));
          
        } else {
          out(k, 1) = obscoefs(0);
        }
        
        if (ob3(i) == 1) {
          
          if (sizel > 1) {
            if (sizedist == 0) {
              // Poisson size distribution
              
              double sizefac = sz3(i) * tgamma(sz3(i));
              double lambda = exp(sizecoefs(0) + (sizecoefs(4) * sz2n(i)) + (sizecoefs(3) * sz1(i)) + 
                                  (sizecoefs(2) * fl2n(i)) + (sizecoefs(1) * fl1(i)) + (sizecoefs(6) * sz2n(i) * sz1(i)) +
                                  (sizecoefs(5) * fl2n(i) * fl1(i)) + (sizecoefs(7) * sz1(i) * fl1(i)) +
                                  (sizecoefs(8) * sz2n(i) * fl2n(i)) + (sizecoefs(10) * sz1(i) * fl2n(i)) + 
                                  (sizecoefs(9) * sz2n(i) * fl1(i)) + sizepatch(patchnumber) + sizeyear(yearnumber) + 
                                  sizedev + (summedvars / 2));
              
              out(k, 3) = ((pow(lambda, sz3(i)) * exp(-1 * lambda)) / sizefac);
              
            } else if (sizedist == 1) {
              // Negative binomial size distribution
              
              double mu = exp(sizecoefs(0) + (sizecoefs(4) * sz2n(i)) + (sizecoefs(3) * sz1(i)) + 
                              (sizecoefs(2) * fl2n(i)) + (sizecoefs(1) * fl1(i)) + (sizecoefs(6) * sz2n(i) * sz1(i)) +
                              (sizecoefs(5) * fl2n(i) * fl1(i)) + (sizecoefs(7) * sz1(i) * fl1(i)) +
                              (sizecoefs(8) * sz2n(i) * fl2n(i)) + (sizecoefs(10) * sz1(i) * fl2n(i)) + 
                              (sizecoefs(9) * sz2n(i) * fl1(i)) + sizepatch(patchnumber) + sizeyear(yearnumber) + sizedev);
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
              preout(3) = (sizecoefs(0) + (sizecoefs(4) * sz2n(i)) + (sizecoefs(3) * sz1(i)) + 
                (sizecoefs(2) * fl2n(i)) + (sizecoefs(1) * fl1(i)) + (sizecoefs(6) * sz2n(i) * sz1(i)) +
                (sizecoefs(5) * fl2n(i) * fl1(i)) + (sizecoefs(7) * sz1(i) * fl1(i)) +
                (sizecoefs(8) * sz2n(i) * fl2n(i)) + (sizecoefs(10) * sz1(i) * fl2n(i)) + 
                (sizecoefs(9) * sz2n(i) * fl1(i)) + sizepatch(patchnumber) + sizeyear(yearnumber) + sizedev);
              
              out(k, 3) = (exp(-1 * (pow((sz3(i) - preout(3)), 2) / (2*sigma2))) / ((pow((2 * M_PI), 0.5)) * sigma)) * binwidth3(i);
              
            } else {
              out(k, 3) = sizedist - 3;
            }
          }
          
          if (repstl > 1) {
            
            preout(2) = (repstcoefs(0) + (repstcoefs(4) * sz2n(i)) + (repstcoefs(3) * sz1(i)) + 
              (repstcoefs(2) * fl2n(i)) + (repstcoefs(1) * fl1(i)) + (repstcoefs(6) * sz2n(i) * sz1(i)) +
              (repstcoefs(5) * fl2n(i) * fl1(i)) + (repstcoefs(7) * sz1(i) * fl1(i)) +
              (repstcoefs(8) * sz2n(i) * fl2n(i)) + (repstcoefs(10) * sz1(i) * fl2n(i)) + 
              (repstcoefs(9) * sz2n(i) * fl1(i))  + repstpatch(patchnumber) + repstyear(yearnumber) + repstdev);
            
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
        
        survtransmat(k) = out(k, 0) * out(k, 1) * out(k, 2) * out(k, 3);
        
        
      } else if (immat2n(i) == 1 && immat1(i) == 1 && jsurvl > 0) {
        // Juvenile to adult transitions
        
        arma::vec preout(4);
        
        if (jsurvl > 1) {
          preout(0) = (jsurvcoefs(0) + jsurvpatch(patchnumber) + jsurvyear(yearnumber) + jsurvdev);
          
          out(k, 0) = exp(preout(0)) / (1 + exp(preout(0)));
        } else {
          out(k, 0) = jsurvcoefs(0);
        }
        
        if (jobsl > 1) {
          preout(1) = (jobscoefs(0) + jobspatch(patchnumber) + jobsyear(yearnumber) + jobsdev);
          
          out(k, 1) = exp(preout(1)) / (1 + exp(preout(1)));
          
        } else {
          out(k, 1) = jobscoefs(0);
        }
        
        if (ob3(i) == 1) {
          if (jsizel > 1) {
            if (sizedist == 0) {
              // Poisson size distribution
              
              double sizefac = sz3(i) * tgamma(sz3(i));
              double lambda = exp(jsizecoefs(0) + jsizepatch(patchnumber) + jsizeyear(yearnumber) + jsizedev + (jsummedvars / 2));
              
              out(k, 3) = ((pow(lambda, sz3(i)) * exp(-1 * lambda)) / sizefac);
              
            } else if (sizedist == 1) {
              // Negative binomial size distribution
              
              double mu = exp(jsizecoefs(0) + jsizepatch(patchnumber) + jsizeyear(yearnumber) + jsizedev);
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
              preout(3) = (jsizecoefs(0) + jsizepatch(patchnumber) + jsizeyear(yearnumber) + jsizedev);
              
              out(k, 3) = (exp(-1 * (pow((sz3(i) - preout(3)), 2) / (2*sigma2))) / ((pow((2 * M_PI), 0.5)) * jsigma)) * binwidth3(i);
              
            } else {
              out(k, 3) = sizedist - 3;
            }
          }
          
          if (jrepstl > 1) {
            
            preout(2) = (jrepstcoefs(0) + jrepstpatch(patchnumber) + jrepstyear(yearnumber) + jrepstdev);
            
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
        
        survtransmat(k) = out(k, 0) * out(k, 1) * out(k, 2) * out(k, 3);
        
      }
    } else if (ovgivent(i) != -1) {
      // All other transitions
      
      survtransmat(k) = ovgivent(i);
    }
    
    if (indata2(i) == 1 && fecl > 0) {
      // Fecundity
      if (fl2n(i) == 1 && ovgivenf(i) == -1) {
        
        double preoutx;
        
        if (fecdist < 3) {
          preoutx = (feccoefs(0) + (feccoefs(4) * sz2n(i)) + (feccoefs(3) * sz1(i)) + 
            (feccoefs(2) * fl2n(i)) + (feccoefs(1) * fl1(i)) + (feccoefs(6) * sz2n(i) * sz1(i)) +
            (feccoefs(5) * fl2n(i) * fl1(i)) + (feccoefs(7) * sz1(i) * fl1(i)) +
            (feccoefs(8) * sz2n(i) * fl2n(i)) + (feccoefs(10) * sz1(i) * fl2n(i)) + 
            (feccoefs(9) * sz2n(i) * fl1(i)) + fecpatch(patchnumber) + fecyear(yearnumber) + fecdev);
          
          if (fecdist != 2) {
            // Poisson and negative binomial fecundity
            
            fectransmat(k) = exp(preoutx) * fecmod * repentry(i);
          } else {
            // Gaussian fecundity
            fectransmat(k) = preoutx * fecmod * repentry(i);
            if (negfec && fectransmat(k) < 0) {
              fectransmat(k) = 0;
            }
          }
          
        } else {
          // All others
          
          fectransmat(k) = fecdist - 3;
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
  
  return List::create(Named("A") = amatrix, _["U"] = survtransmat, _["F"] = fectransmat);
}


