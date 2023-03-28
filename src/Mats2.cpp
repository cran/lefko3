#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <LefkoUtils.h>

using namespace Rcpp;
using namespace arma;
using namespace LefkoUtils;
using namespace LefkoMats;

//' Create Skeleton Stageframe
//' 
//' Function \code{sf_skeleton()} creates a skeleton \code{stageframe} object.
//' 
//' @name sf_skeleton
//' 
//' @param stages The number of stages, as an integer.
//' @param standard A logical value indicating whether to create a standard
//' \code{stageframe} object (\code{TRUE}, the default), or a reassessed
//' \code{stageframe} object as created by function \code{mpm_create()}
//' (\code{FALSE}).
//' 
//' @return A dat frame of class \code{stageframe}.
//' 
//' @export
// [[Rcpp::export(sf_skeleton)]]
Rcpp::DataFrame sf_skeleton(int stages, bool standard = true) { 
  
  bool reassessed = false;
  if (!standard) reassessed = true;
  
  DataFrame new_sf = LefkoMats::sf_core(stages, reassessed, false);
  
  return new_sf;
}

//' Create Element Index for Matrix Estimation
//' 
//' Function \code{simplepizzle()} creates a data frame object used by function
//' \code{\link{hist_null}()} to provide an index for estimation of null
//' historical matrices from ahistorical MPM inputs.
//' 
//' @name simplepizzle
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
//' in the new historical matrices.}
//' 
//' @keywords internal
//' @noRd
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
  int nostages = static_cast<int>(newstageid.n_elem);
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
    new_stageframe.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER,
      static_cast<int>(newstageidvec.n_elem));
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
  hstages.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER,
    stid2.length());
  hstages.attr("class") = "data.frame";
  
  // Set up vectors to put together into matrix map data frame
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
  
  // Main data frame creation loops
  // If style = 0, this will create AllStages for the historical case
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
  
  int stage3_length = static_cast<int>(stage3.n_elem);
  
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

//' Create Historically Structured Version of ahMPM
//' 
//' Function \code{thefifthhousemate()} takes an ahistorical MPM as input, and
//' uses the \code{allstages} index to create a historically structured version
//' of it.
//' 
//' @name thefifthhousemate
//' 
//' @param mpm The original ahMPM, supplied as a \code{lefkoMat} object.
//' @param allstages The index dataframe named \code{allstages}, in the third
//' element of output developed by \code{simplepizzle()}.
//' @param hstages The index dataframe named \code{hstages}, in the second
//' element of output developed by \code{simplepizzle()}.
//' @param stageframe The ahistorical stageframe supplied by
//' \code{simplepizzle()}.
//' @param format Integer indicating whether historical matrices should be in
//' (1) Ehrlen or (2) deVries format.
//' 
//' @return This will return a list of lists. The first list is composed of all
//' new \code{A} matrices. The second list is composed of all new \code{U}
//' matrices. The third list is composed of all new \code{F} matrices.
//' 
//' @keywords internal
//' @noRd
Rcpp::List thefifthhousemate(List mpm, DataFrame allstages, DataFrame hstages,
  DataFrame stageframe, int format) {
  Rcpp::List old_Umats = as<List>(mpm["U"]);
  Rcpp::List old_Fmats = as<List>(mpm["F"]);
  
  Rcpp::IntegerVector stageid = stageframe["stage_id"];
  int nostages = stageid.length();
  int nocols = nostages * nostages;
  
  if (format == 2) nocols = (nostages -1) * nostages;
  
  Rcpp::IntegerVector old_index = allstages["index21"];
  Rcpp::IntegerVector new_indexu = allstages["index321u"];
  Rcpp::IntegerVector new_indexf = allstages["index321f"];
  arma::ivec stage2o_all = as<arma::ivec>(allstages["stage2o"]);
  arma::ivec stage1_all = as<arma::ivec>(allstages["stage1"]);
  arma::ivec hst_s_id_1 = as<arma::ivec>(hstages["stage_id_1"]);
  arma::ivec hst_s_id_2 = as<arma::ivec>(hstages["stage_id_2"]);
  
  int num_mats = old_Umats.length();
  int index_elems = new_indexu.length();
  
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
      if (new_indexu(j) > -1.0) {
        new_U(new_indexu(j)) = old_U(old_index(j));
      }
      if (new_indexf(j) > -1.0) {
        new_F(new_indexf(j)) = old_F(old_index(j));
      }
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

//' Create Historical MPMs Assuming No Influence of Individual History
//' 
//' Function \code{hist_null()} uses ahistorical MPMs to create the equivalent
//' MPMs in the structure of historical MPMs. These MPMs have the same
//' dimensions and stage structure of hMPMs but assume no influence of
//' individual history, and so can be compared to actual hMPMs.
//' 
//' @name hist_null
//' 
//' @param mpm An ahistorical MPM of class \code{lefkoMat}.
//' @param format An integer stipulating whether historical matrices should be
//' produced in Ehrlen format (\code{1}) or deVries format (\code{2}).
//' @param err_check A logical value indicating whether to output the main index
//' data frames used to sort elements in the matrices.
//' 
//' @return An object of class \code{lefkoMat}, with the same list structure as
//' the input object, but with \code{A}, \code{U}, and \code{F} elements
//' replaced with lists of historically-structured matrices, and with element
//' \code{hstages} changed from \code{NA} to an index of stage pairs
//' corresponding to the rows and columns of the new matrices. If
//' \code{err_check = TRUE}, then a list of three data frames showing the values
//' used to determine matrix element index values is also exported.
//' 
//' @section Notes:
//' This function does not currently identify biologically impossible
//' transitions. Ahistorical transition values are placed in all theoretically
//' possible historical transitions.
//' 
//' @examples
//' sizevector <- c(1, 1, 2, 3)
//' stagevector <- c("Sdl", "Veg", "SmFlo", "LFlo")
//' repvector <- c(0, 0, 1, 1)
//' obsvector <- c(1, 1, 1, 1)
//' matvector <- c(0, 1, 1, 1)
//' immvector <- c(1, 0, 0, 0)
//' propvector <- c(0, 0, 0, 0)
//' indataset <- c(1, 1, 1, 1)
//' binvec <- c(0.5, 0.5, 0.5, 0.5)
//' 
//' anthframe <- sf_create(sizes = sizevector, stagenames = stagevector,
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
//'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
//'   propstatus = propvector)
//' 
//' # POPN C 2003-2004
//' XC3 <- matrix(c(0, 0, 1.74, 1.74,
//' 0.208333333, 0, 0, 0.057142857,
//' 0.041666667, 0.076923077, 0, 0,
//' 0.083333333, 0.076923077, 0.066666667, 0.028571429), 4, 4, byrow = TRUE)
//' 
//' # 2004-2005
//' XC4 <- matrix(c(0, 0, 0.3, 0.6,
//' 0.32183908, 0.142857143, 0, 0,
//' 0.16091954, 0.285714286, 0, 0,
//' 0.252873563, 0.285714286, 0.5, 0.6), 4, 4, byrow = TRUE)
//' 
//' mats_list <- list(XC3, XC4)
//' yr_ord <- c(1, 2)
//' pch_ord <- c(1, 1)
//' 
//' anth_lefkoMat <- create_lM(mats_list, anthframe, hstages = NA, historical = FALSE,
//'   poporder = 1, patchorder = pch_ord, yearorder = yr_ord)
//'   
//' nullmodel1 <- hist_null(anth_lefkoMat, 1) # Ehrlen format
//' nullmodel2 <- hist_null(anth_lefkoMat, 2) # deVries format
//' 
//' @export hist_null
// [[Rcpp::export(hist_null)]]
Rcpp::List hist_null (RObject mpm, int format = 1, bool err_check = false) {
  
  List mpm_;
  DataFrame mpm_ahstages;
  DataFrame mpm_agestages;
  DataFrame mpm_labels;
  IntegerVector mpm_data_qc;
  DataFrame mpm_model_qc;
  
  if (is<List>(mpm)) {
    List mpm_thru(mpm);
    mpm_ = mpm_thru;
    
    StringVector mpm_class_vec = mpm_.attr("class");
    std::string mpm_class = as<std::string>(mpm_class_vec[0]);
    if (!LefkoUtils::stringcompare_hard(mpm_class, "lefkoMat")) {
      throw Rcpp::exception("Object mpm is of unrecognized class.", false);
    }
    
    DataFrame mpm_ahstages_ = as<DataFrame>(mpm_["ahstages"]);
    DataFrame mpm_hstages_ = as<DataFrame>(mpm_["hstages"]);
    DataFrame mpm_agestages_ = as<DataFrame>(mpm_["agestages"]);
    DataFrame mpm_labels_ = as<DataFrame>(mpm_["labels"]);
    
    if (mpm_ahstages.length() == 1) {
      throw Rcpp::exception("Object mpm is of unrecognized class.", false);
    }
    mpm_ahstages = mpm_ahstages_;
    
    if (mpm_hstages_.length() != 1) {
      throw Rcpp::exception("Input MPM must be ahistorical.", false);
    }

    if (mpm_agestages_.length() != 1) {
      throw Rcpp::exception("Input MPM must be ahistorical, and cannot be age-by-stage.",
        false);
    }
    mpm_agestages = mpm_agestages_;
    
    if (mpm_labels_.length() == 1) {
      throw Rcpp::exception("Object mpm is of unrecognized class.", false);
    }
    mpm_labels = mpm_labels_;
    
  } else {
    throw Rcpp::exception("Object mpm is of unrecognized class.", false);
  }
  
  int tot_dims = mpm_.length();
  int datqc_position = -1;
  int modqc_position = -1;
  
  StringVector mpm_components = mpm_.names();
  for (int i = 0; i < tot_dims; i++) {
    if (stringcompare_hard(as<std::string>(mpm_components(i)), "modelqc")) {
      modqc_position = i;
    }
    if (stringcompare_hard(as<std::string>(mpm_components(i)), "dataqc")) {
      datqc_position = i;
    }
  }
  
  if (datqc_position != -1) {
    IntegerVector mpm_data_qc_ = as<IntegerVector>(mpm_["dataqc"]);
    mpm_data_qc = mpm_data_qc_;
  }
  if (modqc_position != -1) {
    DataFrame mpm_model_qc_ = as<DataFrame>(mpm_["modelqc"]);
    mpm_model_qc = mpm_model_qc_;
  }
  
  List allstages = simplepizzle(mpm_ahstages, format);
  
  DataFrame allindices = as<DataFrame>(allstages["allstages"]);
  DataFrame all_hstages = as<DataFrame>(allstages["hstages"]);
  DataFrame all_ahstages = as<DataFrame>(allstages["ahstages"]);
  
  List redone_mpms = thefifthhousemate(mpm_, allindices, all_hstages,
    all_ahstages, format);
  
  List A_mpms = as<List>(redone_mpms["A"]);
  List U_mpms = as<List>(redone_mpms["U"]);
  List F_mpms = as<List>(redone_mpms["F"]);
  
  int tot_U_elems {0};
  int tot_F_elems {0};
  int tot_mats = A_mpms.length();
  
  for (int i = 0; i < tot_mats; i++) {
    arma::mat current_Umat = as<arma::mat>(U_mpms(i));
    arma::mat current_Fmat = as<arma::mat>(F_mpms(i));
    
    arma::uvec found_U_uvec = find(current_Umat);
    arma::uvec found_F_uvec = find(current_Fmat);
    int found_U = static_cast<int>(found_U_uvec.n_elem);
    int found_F = static_cast<int>(found_F_uvec.n_elem);
    
    tot_U_elems += found_U;
    tot_F_elems += found_F;
  }
  
  IntegerVector new_mqc = {tot_U_elems, tot_F_elems, tot_mats};
  
  if (err_check) tot_dims++;
  
  List final_output(tot_dims);
  CharacterVector listnames (tot_dims);
  
  final_output(0) = A_mpms;
  listnames(0) = "A";
  final_output(1) = U_mpms;
  listnames(1) = "U";
  final_output(2) = F_mpms;
  listnames(2) = "F";
  final_output(3) = all_hstages;
  listnames(3) = "hstages";
  final_output(4) = mpm_agestages;
  listnames(4) = "agestages";
  final_output(5) = all_ahstages;
  listnames(5) = "ahstages";
  final_output(6) = mpm_labels;
  listnames(6) = "labels";
  final_output(7) = new_mqc;
  listnames(7) = "matrixqc";
  
  int next_dude = 8;
  if (modqc_position > -1) {
    final_output(next_dude) = mpm_model_qc;
    listnames(next_dude) = "modelqc";
    next_dude++;
  }
  if (datqc_position > -1) {
    final_output(next_dude) = mpm_data_qc;
    listnames(next_dude) = "dataqc";
  }
  
  if (err_check) {
    final_output(tot_dims - 1) = allstages;
    listnames(tot_dims - 1) = "err_check";
  }
  
  final_output.attr("names") = listnames;
  StringVector needed_classes {"lefkoMat"};
  final_output.attr("class") = needed_classes;
  
  return final_output;
}

//' Estimate Mean Projection Matrices
//' 
//' Function \code{lmean()} estimates mean projection matrices as element-wise
//' arithmetic means. It produces \code{lefkoMat} objects if provided with them,
//' or single matrices in a simple one-element list if provided a list of
//' matrices.
//' 
//' @name lmean
//' 
//' @param mats A \code{lefkoMat} object, or a list of square matrices of equal
//' dimension.
//' @param matsout A string identifying which means to estimate. Option
//' \code{"pop"} indicates population-level only, \code{"patch"} indicates
//' patch-level only, and \code{"all"} indicates that both patch- and
//' population-level means should be estimated. Defaults to \code{"all"}.
//' @param force_sparse A logical value identifying whether to output the mean
//' matrices in sparse format, if input as standard matrices.
//' 
//' @return Yields a \code{lefkoMat} object with the following characteristics:
//' 
//' \item{A}{A list of full mean projection matrices in order of sorted
//' populations, patches, and years. These are typically estimated as the sums
//' of the associated mean \code{U} and \code{F} matrices. All matrices output
//' in either the \code{matrix} class, or the \code{dgCMatrix} class.}
//' \item{U}{A list of mean survival-transition matrices sorted as in \code{A}.
//' All matrices output in the \code{matrix} class.}
//' \item{F}{A list of mean fecundity matrices sorted as in \code{A}. All
//' matrices output in the \code{matrix} class.}
//' \item{hstages}{A data frame showing the pairing of ahistorical stages used
//' to create historical stage pairs. Given if the MPM is historical.}
//' \item{ahstages}{A data frame detailing the characteristics of associated
//' ahistorical stages.}
//' \item{labels}{A data frame detailing the order of population, patch, and
//' year of each mean matrix. If \code{pop}, \code{patch}, or \code{year2} are
//' \code{NA} in the original \code{labels} set, then these will be re-labeled
//' as \code{A}, \code{1}, or \code{1}, respectively.}
//' \item{matrixqc}{A short vector describing the number of non-zero elements in
//' \code{U} and \code{F} mean matrices, and the number of annual matrices.}
//' \item{modelqc}{This is the \code{qc} portion of the \code{modelsuite} input.
//' Only output from \code{lefkoMat} objects resulting from function-based
//' estimation.}
//' \item{dataqc}{A vector showing the numbers of individuals and rows in the
//' vertical dataset used as input. Only output from \code{lefkoMat} objects
//' resulting from raw matrix estimation.}
//' 
//' @examples
//' data(cypdata)
//' 
//' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
//' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
//'   "XLg")
//' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
//' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
//' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
//' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
//' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
//' 
//' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
//'   propstatus = propvector, immstatus = immvector, indataset = indataset,
//'   binhalfwidth = binvec)
//' 
//' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
//'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
//'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
//'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
//'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
//'   NRasRep = TRUE)
//' 
//' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
//'     "XSm", "Sm", "SD", "P1"),
//'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
//'     "rep"),
//'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
//'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
//'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
//'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
//'   stageframe = cypframe_raw, historical = FALSE)
//' 
//' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
//'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
//'   size = c("size3added", "size2added"), supplement = cypsupp2r,
//'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
//' 
//' cyp2mean <- lmean(cypmatrix2r)
//' 
//' @export
// [[Rcpp::export(lmean)]]
Rcpp::List lmean(RObject mats, Nullable<String> matsout = R_NilValue,
  bool force_sparse = false) {
  
  StringVector mats_class;
  StringVector mats_elements;
  String matsout_;
  
  List u_mats;
  List f_mats;
  List a_mats;
  
  DataFrame listofyears;
  DataFrame ahstages;
  DataFrame agestages;
  DataFrame hstages;
  DataFrame modelqc;
  IntegerVector dataqc;
  
  bool historical {false};
  bool poponly {false};
  bool patchonly {false};
  bool dataqc_found {false};
  bool modelqc_found {false};
  bool mat_input {true};
  
  List gd_output;
  
  if (matsout.isNotNull()) {
    String matsout__(matsout);
    matsout_ = matsout__;
    StringVector matsout_choices = {"pop", "patch", "all"};
    int found_one {0};
    for (int i = 0; i < 3; i++) {
      if (stringcompare_simple(matsout_, String(matsout_choices(i)))) {
        found_one++;
      }
    }
    if (found_one == 0 && matsout_ == "") {
      matsout_ = "all";
    } else if (found_one != 1) { 
      throw Rcpp::exception("Argument matsout must equal 'all', 'pop', or 'patch'.",
        false);
    }
  } else {
    matsout_ = "all";
  }
  
  if (stringcompare_simple(matsout_, "all")) { 
    poponly = true;
    patchonly = true;
    
  } else if (stringcompare_simple(matsout_, "pat")) {
    patchonly = true;
    
  } else if (stringcompare_simple(matsout_, "pop")) { 
    poponly = true;
    
  }
  
  if (is<NumericMatrix>(mats)) {
    throw Rcpp::exception("Matrix averaging cannot be performed on a single matrix.",
      false);
    
  } else if (is<List>(mats)) {
    List mats_ = as<List>(mats);
    if (mats_.hasAttribute("class")) {
      mats_class = as<StringVector>(mats_.attr("class"));
    }
    if (mats_.hasAttribute("names")) {
      mats_elements = as<StringVector>(mats_.attr("names"));
    }
    
    int mats_class_length = mats_class.length();
    bool found_it {false};
    for (int i = 0; i < mats_class_length; i++) {
      if (stringcompare_hard(String(mats_class(i)), "lefkoMat")) {
        found_it = true;
      }
    }
    
    if (found_it) {
      // lefkoMat input
      StringVector lefkoMat_required = {"A", "U", "F", "ahstages", "labels"};
      int found_element {0};
      
      for (int i = 0; i < lefkoMat_required.length(); i++) {
        for (int j = 0; j < mats_elements.length(); j++) {
          if (stringcompare_hard(String(lefkoMat_required(i)), String(mats_elements(j)))) {
            found_element++;
          }
        }
      }
      if (found_element < 5) {
        throw Rcpp::exception("Argument mats does not appear to be a lefkoMat object.",
          false);
      }
      
      u_mats = as<List>(mats_["U"]);
      f_mats = as<List>(mats_["F"]);
      ahstages = as<DataFrame>(mats_["ahstages"]);
      agestages = as<DataFrame>(mats_["agestages"]);
      hstages = as<DataFrame>(mats_["hstages"]);
      
      CharacterVector hstages_names = hstages.names();
      
      if (hstages_names.length() > 1) historical = true;
      
      for (int i = 0; i < mats_elements.length(); i++) {
        if (stringcompare_hard(String(mats_elements(i)), "dataqc")) {
          dataqc = as<IntegerVector>(mats_["dataqc"]);
          dataqc_found = true;
        }
        
        if (stringcompare_hard(String(mats_elements(i)), "modelqc")) {
          modelqc = as<DataFrame>(mats_["modelqc"]);
          modelqc_found = true;
        }
      }
      
      // Quality control for listofyears
      DataFrame labels = as<DataFrame>(mats_["labels"]);
      listofyears = loy_inator(labels, true);
      
      if (is<S4>(u_mats(0))) mat_input = false;
      
      // Matrix averaging
      if (historical) {
        gd_output = turbogeodiesel(listofyears, u_mats, f_mats, hstages, 
          agestages, ahstages, patchonly, poponly, mat_input, 0);
        
      } else {
        gd_output = geodiesel(listofyears, u_mats, f_mats, agestages, ahstages,
          patchonly, poponly, mat_input, 0);
      }
    
      gd_output["ahstages"] = ahstages;
      gd_output["hstages"] = hstages;
      gd_output["agestages"] = agestages;
      if (dataqc_found) gd_output["dataqc"] = dataqc;
      if (modelqc_found) gd_output["modelqc"] = modelqc;
      
      StringVector gd_class = {"lefkoMat"};
      gd_output.attr("class") = gd_class;
      
    } else {
      // List of NumericMatrix input
      int mats_length = mats_.length();
      
      if (is<S4>(mats_(0))) mat_input = false;
      
      arma::mat core_mat;
      arma::sp_mat core_mat_sp;
      int mat_rows {0};
      int mat_cols {0};
      
      for (int i = 0; i < mats_length; i++) {
        RObject test_bit = as<RObject>(mats_[i]);
        
        if (!is<NumericMatrix>(test_bit) && !is<S4>(test_bit)) {
          String eat_my_shorts = "Argument mats must be either a lefkoMat object ";
          String eat_my_shorts1 = "or a list of numeric matrices.";
          eat_my_shorts += eat_my_shorts1;
          
          throw Rcpp::exception(eat_my_shorts.get_cstring(), false);
        }
        
        if (i == 0) {
          if (mat_input) {
            core_mat = as<arma::mat>(mats_[0]);
            
            mat_rows = core_mat.n_rows;
            mat_cols = core_mat.n_cols;
          } else {
            core_mat_sp = as<arma::sp_mat>(mats_[0]);
            
            mat_rows = core_mat_sp.n_rows;
            mat_cols = core_mat_sp.n_cols;
          }
          
          if (mat_rows != mat_cols) {
            throw Rcpp::exception("Input matrices must be square.", false);
          }
          
          if (mat_input) {
            core_mat = core_mat / mats_length;
          } else {
            core_mat_sp = core_mat_sp / mats_length;
          }
        } else {
          if (mat_input) {
            arma::mat next_mat = as<arma::mat>(mats_[i]);
            
            int next_rows = next_mat.n_rows;
            int next_cols = next_mat.n_cols;
            
            if (next_rows != mat_rows || next_cols != mat_cols) {
              throw Rcpp::exception("All input matrices must have the same dimensions.",
                false);
            }
            
            core_mat = core_mat + (next_mat / mats_length);
          } else {
            arma::sp_mat next_mat_sp = as<arma::sp_mat>(mats_[i]);
            
            int next_rows = next_mat_sp.n_rows;
            int next_cols = next_mat_sp.n_cols;
            
            if (next_rows != mat_rows || next_cols != mat_cols) {
              throw Rcpp::exception("All input matrices must have the same dimensions.",
                false);
            }
            
            core_mat_sp = core_mat_sp + (next_mat_sp / mats_length);
          }
        }
        
        if (mat_input && force_sparse) {
          arma::sp_mat core_mat_sp_(core_mat);
          core_mat_sp = core_mat_sp_;
        }
        
        if (mat_input && !force_sparse) {
          List A = List::create(_["A"] = core_mat);
          gd_output = A;
        } else {
          List A = List::create(_["A"] = core_mat_sp);
          gd_output = A;
        }
      }
    }
  }
  
  return gd_output;
}

