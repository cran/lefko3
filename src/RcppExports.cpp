// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// ovreplace
arma::mat ovreplace(arma::vec allst321, arma::vec idx321old, arma::vec idx321new, arma::vec convtype, arma::vec eststag3, arma::vec gvnrate);
RcppExport SEXP _lefko3_ovreplace(SEXP allst321SEXP, SEXP idx321oldSEXP, SEXP idx321newSEXP, SEXP convtypeSEXP, SEXP eststag3SEXP, SEXP gvnrateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type allst321(allst321SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type idx321old(idx321oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type idx321new(idx321newSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type convtype(convtypeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eststag3(eststag3SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gvnrate(gvnrateSEXP);
    rcpp_result_gen = Rcpp::wrap(ovreplace(allst321, idx321old, idx321new, convtype, eststag3, gvnrate));
    return rcpp_result_gen;
END_RCPP
}
// specialpatrolgroup
arma::mat specialpatrolgroup(arma::mat sge9l, arma::mat sge3, arma::mat maindata, arma::uvec sge93index, arma::uvec sge92index, arma::uvec sge32index, arma::uvec sge33, arma::uvec sge32, arma::uvec data3221, arma::uvec data21, int nostages);
RcppExport SEXP _lefko3_specialpatrolgroup(SEXP sge9lSEXP, SEXP sge3SEXP, SEXP maindataSEXP, SEXP sge93indexSEXP, SEXP sge92indexSEXP, SEXP sge32indexSEXP, SEXP sge33SEXP, SEXP sge32SEXP, SEXP data3221SEXP, SEXP data21SEXP, SEXP nostagesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sge9l(sge9lSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sge3(sge3SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type maindata(maindataSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type sge93index(sge93indexSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type sge92index(sge92indexSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type sge32index(sge32indexSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type sge33(sge33SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type sge32(sge32SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type data3221(data3221SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type data21(data21SEXP);
    Rcpp::traits::input_parameter< int >::type nostages(nostagesSEXP);
    rcpp_result_gen = Rcpp::wrap(specialpatrolgroup(sge9l, sge3, maindata, sge93index, sge92index, sge32index, sge33, sge32, data3221, data21, nostages));
    return rcpp_result_gen;
END_RCPP
}
// normalpatrolgroup
arma::mat normalpatrolgroup(arma::mat sge9l, arma::mat sge3, arma::mat maindata, arma::uvec sge93index, arma::uvec sge92index, arma::uvec sge32index, arma::uvec sge32, arma::uvec data3221, arma::uvec data21, int nostages);
RcppExport SEXP _lefko3_normalpatrolgroup(SEXP sge9lSEXP, SEXP sge3SEXP, SEXP maindataSEXP, SEXP sge93indexSEXP, SEXP sge92indexSEXP, SEXP sge32indexSEXP, SEXP sge32SEXP, SEXP data3221SEXP, SEXP data21SEXP, SEXP nostagesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type sge9l(sge9lSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sge3(sge3SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type maindata(maindataSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type sge93index(sge93indexSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type sge92index(sge92indexSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type sge32index(sge32indexSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type sge32(sge32SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type data3221(data3221SEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type data21(data21SEXP);
    Rcpp::traits::input_parameter< int >::type nostages(nostagesSEXP);
    rcpp_result_gen = Rcpp::wrap(normalpatrolgroup(sge9l, sge3, maindata, sge93index, sge92index, sge32index, sge32, data3221, data21, nostages));
    return rcpp_result_gen;
END_RCPP
}
// hoffmannofstuttgart
arma::uvec hoffmannofstuttgart(int mainindex, arma::mat allstages);
RcppExport SEXP _lefko3_hoffmannofstuttgart(SEXP mainindexSEXP, SEXP allstagesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type mainindex(mainindexSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type allstages(allstagesSEXP);
    rcpp_result_gen = Rcpp::wrap(hoffmannofstuttgart(mainindex, allstages));
    return rcpp_result_gen;
END_RCPP
}
// jerzeibalowski
arma::mat jerzeibalowski(arma::vec survcoefs, arma::vec obscoefs, arma::vec sizecoefs, arma::vec repstcoefs, arma::vec feccoefs, arma::vec jsurvcoefs, arma::vec jobscoefs, arma::vec jsizecoefs, arma::vec jrepstcoefs, arma::vec stage3, arma::vec stage2n, arma::vec stage1, arma::vec sz3, arma::vec sz2n, arma::vec sz1, arma::vec fl3, arma::vec fl2n, arma::vec fl1, arma::vec ob3, arma::vec ob2n, arma::vec ob1, arma::vec immat3, arma::vec immat2n, arma::vec indata2, arma::vec indata, arma::uvec aliveandequal, arma::vec repentry, arma::vec ovgivent, arma::vec ovgivenf, arma::vec binwidth3, unsigned long numofsizes4, double fecmod, double summedvars, double sigma, double jsummedvars, double jsigma, double maxsize, int sizedist, int fecdist);
RcppExport SEXP _lefko3_jerzeibalowski(SEXP survcoefsSEXP, SEXP obscoefsSEXP, SEXP sizecoefsSEXP, SEXP repstcoefsSEXP, SEXP feccoefsSEXP, SEXP jsurvcoefsSEXP, SEXP jobscoefsSEXP, SEXP jsizecoefsSEXP, SEXP jrepstcoefsSEXP, SEXP stage3SEXP, SEXP stage2nSEXP, SEXP stage1SEXP, SEXP sz3SEXP, SEXP sz2nSEXP, SEXP sz1SEXP, SEXP fl3SEXP, SEXP fl2nSEXP, SEXP fl1SEXP, SEXP ob3SEXP, SEXP ob2nSEXP, SEXP ob1SEXP, SEXP immat3SEXP, SEXP immat2nSEXP, SEXP indata2SEXP, SEXP indataSEXP, SEXP aliveandequalSEXP, SEXP repentrySEXP, SEXP ovgiventSEXP, SEXP ovgivenfSEXP, SEXP binwidth3SEXP, SEXP numofsizes4SEXP, SEXP fecmodSEXP, SEXP summedvarsSEXP, SEXP sigmaSEXP, SEXP jsummedvarsSEXP, SEXP jsigmaSEXP, SEXP maxsizeSEXP, SEXP sizedistSEXP, SEXP fecdistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type survcoefs(survcoefsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type obscoefs(obscoefsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sizecoefs(sizecoefsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type repstcoefs(repstcoefsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type feccoefs(feccoefsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type jsurvcoefs(jsurvcoefsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type jobscoefs(jobscoefsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type jsizecoefs(jsizecoefsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type jrepstcoefs(jrepstcoefsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type stage3(stage3SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type stage2n(stage2nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type stage1(stage1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sz3(sz3SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sz2n(sz2nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sz1(sz1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type fl3(fl3SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type fl2n(fl2nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type fl1(fl1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ob3(ob3SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ob2n(ob2nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ob1(ob1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type immat3(immat3SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type immat2n(immat2nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type indata2(indata2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type indata(indataSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type aliveandequal(aliveandequalSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type repentry(repentrySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ovgivent(ovgiventSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ovgivenf(ovgivenfSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type binwidth3(binwidth3SEXP);
    Rcpp::traits::input_parameter< unsigned long >::type numofsizes4(numofsizes4SEXP);
    Rcpp::traits::input_parameter< double >::type fecmod(fecmodSEXP);
    Rcpp::traits::input_parameter< double >::type summedvars(summedvarsSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type jsummedvars(jsummedvarsSEXP);
    Rcpp::traits::input_parameter< double >::type jsigma(jsigmaSEXP);
    Rcpp::traits::input_parameter< double >::type maxsize(maxsizeSEXP);
    Rcpp::traits::input_parameter< int >::type sizedist(sizedistSEXP);
    Rcpp::traits::input_parameter< int >::type fecdist(fecdistSEXP);
    rcpp_result_gen = Rcpp::wrap(jerzeibalowski(survcoefs, obscoefs, sizecoefs, repstcoefs, feccoefs, jsurvcoefs, jobscoefs, jsizecoefs, jrepstcoefs, stage3, stage2n, stage1, sz3, sz2n, sz1, fl3, fl2n, fl1, ob3, ob2n, ob1, immat3, immat2n, indata2, indata, aliveandequal, repentry, ovgivent, ovgivenf, binwidth3, numofsizes4, fecmod, summedvars, sigma, jsummedvars, jsigma, maxsize, sizedist, fecdist));
    return rcpp_result_gen;
END_RCPP
}
// geodiesel
arma::mat geodiesel(int geom, int sparse, int numofpops, int numofpatches, int numofyears, arma::mat loy2c, arma::mat Umats, arma::mat Fmats);
RcppExport SEXP _lefko3_geodiesel(SEXP geomSEXP, SEXP sparseSEXP, SEXP numofpopsSEXP, SEXP numofpatchesSEXP, SEXP numofyearsSEXP, SEXP loy2cSEXP, SEXP UmatsSEXP, SEXP FmatsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type geom(geomSEXP);
    Rcpp::traits::input_parameter< int >::type sparse(sparseSEXP);
    Rcpp::traits::input_parameter< int >::type numofpops(numofpopsSEXP);
    Rcpp::traits::input_parameter< int >::type numofpatches(numofpatchesSEXP);
    Rcpp::traits::input_parameter< int >::type numofyears(numofyearsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type loy2c(loy2cSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Umats(UmatsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Fmats(FmatsSEXP);
    rcpp_result_gen = Rcpp::wrap(geodiesel(geom, sparse, numofpops, numofpatches, numofyears, loy2c, Umats, Fmats));
    return rcpp_result_gen;
END_RCPP
}
// turbogeodiesel
arma::mat turbogeodiesel(int geom, int sparse, int numofpops, int numofpatches, int numofyears, arma::mat loy2c, arma::mat Umats, arma::mat Fmats, arma::mat Amats);
RcppExport SEXP _lefko3_turbogeodiesel(SEXP geomSEXP, SEXP sparseSEXP, SEXP numofpopsSEXP, SEXP numofpatchesSEXP, SEXP numofyearsSEXP, SEXP loy2cSEXP, SEXP UmatsSEXP, SEXP FmatsSEXP, SEXP AmatsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type geom(geomSEXP);
    Rcpp::traits::input_parameter< int >::type sparse(sparseSEXP);
    Rcpp::traits::input_parameter< int >::type numofpops(numofpopsSEXP);
    Rcpp::traits::input_parameter< int >::type numofpatches(numofpatchesSEXP);
    Rcpp::traits::input_parameter< int >::type numofyears(numofyearsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type loy2c(loy2cSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Umats(UmatsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Fmats(FmatsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Amats(AmatsSEXP);
    rcpp_result_gen = Rcpp::wrap(turbogeodiesel(geom, sparse, numofpops, numofpatches, numofyears, loy2c, Umats, Fmats, Amats));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_lefko3_ovreplace", (DL_FUNC) &_lefko3_ovreplace, 6},
    {"_lefko3_specialpatrolgroup", (DL_FUNC) &_lefko3_specialpatrolgroup, 11},
    {"_lefko3_normalpatrolgroup", (DL_FUNC) &_lefko3_normalpatrolgroup, 10},
    {"_lefko3_hoffmannofstuttgart", (DL_FUNC) &_lefko3_hoffmannofstuttgart, 2},
    {"_lefko3_jerzeibalowski", (DL_FUNC) &_lefko3_jerzeibalowski, 39},
    {"_lefko3_geodiesel", (DL_FUNC) &_lefko3_geodiesel, 8},
    {"_lefko3_turbogeodiesel", (DL_FUNC) &_lefko3_turbogeodiesel, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_lefko3(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
