#include <RcppArmadillo.h>
// #include <string.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Append NumericVector to the End of Another NumericVector
//' 
//' This function appends one NumericVector fully to another.
//' 
//' @param A Any NumericVector.
//' @param B Any other NumericVector.
//' 
//' @return Returns a new NumericVector with elements of vector A followed by
//' elements of vector B.
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export(.concat_dbl)]]
NumericVector concat_dbl(NumericVector x, NumericVector y) {
  
  std::vector<double> xconv = as<std::vector<double> >(x);
  std::vector<double> yconv = as<std::vector<double> >(y);
  std::vector<double> z(xconv.size() + yconv.size());
  
  std::copy(xconv.begin(), xconv.end(), z.begin());
  std::copy(yconv.begin(), yconv.end(), z.begin() + xconv.size());
  
  NumericVector zconv(z.begin(), z.end());
  
  return(zconv);
}

//' Append IntegerVector to the End of Another IntegerVector
//' 
//' Returns a new IntegerVector with elements of vector A followed by
//' elements of vector B.
//' 
//' @param A Any IntegerVector.
//' @param B Any other IntegerVector.
//' 
//' @return Returns a new IntegerVector with elements of vector A followed by
//' elements of vector B.
//' 
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export(.concat_int)]]
IntegerVector concat_int(IntegerVector x, IntegerVector y) {
  
  std::vector<long long int> xconv = as<std::vector<long long int> >(x);
  std::vector<long long int> yconv = as<std::vector<long long int> >(y);
  std::vector<long long int> z(x.size() + y.size());
  
  std::copy(xconv.begin(), xconv.end(), z.begin());
  std::copy(yconv.begin(), yconv.end(), z.begin() + xconv.size());
  
  IntegerVector zconv(z.begin(), z.end());
  
  return(zconv);
}

//' Append StringVector to the End of Another StringVector
//' 
//' Returns a new StringVector with elements of vector A followed by
//' elements of vector B.
//' 
//' @param A Any StringVector.
//' @param B Any other StringVector.
//' 
//' @return Returns a new StringVector with elements of vector A followed by
//' elements of vector B.
//' 
//'
//' @keywords internal
//' @noRd
// [[Rcpp::export(.concat_str)]]
StringVector concat_str(StringVector x, StringVector y) {
  
  std::vector<std::string> xconv = as<std::vector<std::string> >(x);
  std::vector<std::string> yconv = as<std::vector<std::string> >(y);
  std::vector<std::string> z(x.size() + y.size());
  
  std::copy(x.begin(), x.end(), z.begin());
  std::copy(y.begin(), y.end(), z.begin() + x.size());
  
  StringVector zconv(z.begin(), z.end());
  
  return(zconv);
}

//' Create Stageframe for Population Matrix Projection Analysis
//' 
//' Function \code{sf_create()} returns a data frame describing each ahistorical
//' life history stage in the life history model. This data frame can be used as 
//' input into MPM creation functions including \code{\link{flefko3}()}, 
//' \code{\link{flefko2}()}, \code{\link{aflefko2}()}, \code{\link{rlefko3}()},
//' and \code{\link{rlefko2}()}, in which it determines how each stage is
//' treated during matrix estimation.
//' 
//' @param sizes A numeric vector of the typical or representative size of each
//' life history stage. If making function-based MPMs, then this should be a
//' vector composed of the midpoints of each size bin. If denoting the boundary
//' of an automated size classification group, then should denote the absolute
//' minimum size of that group, or the absolute size of that group (see Notes).
//' @param stagenames A vector of stage names, in the same order as elements in
//' sizes. Can also be set to \code{ipm} for automated size classification (see
//' Notes section).
//' @param sizesb An optional numeric vector for a second size metric for each
//' life history stage. Only to be used if stages are defined by at least two
//' size metrics in all cases. Same issues apply as in \code{sizes}.
//' @param sizesc An optional numeric vector for a third size metric for each
//' life history stage. Only to be used if stages are defined by at least three
//' size metrics in all cases. Same issues apply as in \code{sizes}.
//' @param repstatus A vector denoting the binomial reproductive status of each
//' life history stage. Defaults to 1.
//' @param obsstatus A vector denoting the binomial observation status of each
//' life history stage. Defaults to 1, but may be changed for unobservable 
//' stages.
//' @param propstatus A vector denoting whether each life history stage is a 
//' propagule. Such stages are generally only used in fecundity estimation. 
//' Defaults to 0.
//' @param matstatus A vector denoting whether each stage is mature. Must be
//' composed of binomial values if given. Defaults to 1 for all stages defined 
//' in \code{sizes}.
//' @param immstatus A vector denoting whether each stage is immature. Must be
//' composed of binomial values if given. Defaults to the complement of vector
//' \code{matstatus}.
//' @param minage An optional vector denoting the minimum age at which a stage
//' can occur. Only used in age x stage matrix development. Defaults to NA.
//' @param maxage An optional vector denoting the maximum age at which a stage
//' should occur. Only used in age x stage matrix development. Defaults to NA.
//' @param indataset A vector designating which stages are found within the 
//' dataset. While \code{\link{rlefko2}()} and \code{\link{rlefko3}()} can use
//' all stages in the input dataset, \code{\link{flefko3}()} and
//' \code{\link{flefko2}()} can only handle size-classified stages with
//' non-overlapping combinations of size and status variables. Stages that do
//' not actually exist within the dataset should be marked as 0 in this vector.
//' @param binhalfwidth A numeric vector giving the half-width of size bins.
//' Required to classify individuals appropriately within size classes.
//' Defaults to 0.5 for all sizes.
//' @param binhalfwidthb A numeric vector giving the half-width of size bins
//' used for the optional second size metric. Required to classify individuals
//' appropriately with two or three size classes. Defaults to 0.5 for all sizes.
//' @param binhalfwidthc A numeric vector giving the half-width of size bins
//' used for the optional third size metric. Required to classify individuals
//' appropriately with three size classes. Defaults to 0.5 for all sizes.
//' @param group An integer vector providing information on each respective
//' stage's size classification group. If used, then function-based MPM creation
//' functions \code{\link{flefko2}()}, \code{\link{flefko3}()}, and
//' \code{\link{aflefko2}()} will estimate transitions only within these groups
//' and for allowed cross-group transitions noted within the supplement table.
//' Defaults to 0.
//' @param comments An optional vector of text entries holding useful text
//' descriptions of all stages.
//' @param roundsize This parameter sets the precision of size classification,
//' and equals the number of digits used in rounding sizes. Defaults to 5.
//' @param roundsizeb This parameter sets the precision of size classification
//' in the optional second size metric, and equals the number of digits used in
//' rounding sizes. Defaults to 5.
//' @param roundsizec This parameter sets the precision of size classification
//' in the optional third size metric, and equals the number of digits used in
//' rounding sizes. Defaults to 5.
//' @param ipmbins An integer giving the number of size bins to create using the
//' primary size classification variable. This number is in addition to any
//' stages that are not size classified. Defaults to 100, and numbers greater
//' than this yield a warning about the loss of statistical power and increasing
//' chance of matrix over-parameterization resulting from increasing numbers of
//' stages.
//' @param ipmbinsb An optional integer giving the number of size bins to create
//' using the secondary size classification variable. This number is in addition
//' to any stages that are not size classified, as well as in addition to any
//' automated size classification using the primary and tertiary size variables.
//' Defaults to NA, and must be set to a positive integer for automated size
//' classification to progress.
//' @param ipmbinsc An optional integer giving the number of size bins to create
//' using the tertiary size classification variable. This number is in addition
//' to any stages that are not size classified, as well as in addition to any
//' automated size classification using the primary and secondary size
//' variables. Defaults to NA, and must be set to a positive integer for
//' automated size classification to progress.
//' 
//' @return A data frame of class \code{stageframe}, which includes information
//' on the stage name, size, reproductive status, observation status, propagule 
//' status, immaturity status, maturity status, presence within the core dataset, 
//' stage group classification, raw bin half-width, and the minimum, 
//' center, and maximum of each size bin, as well as its width. If minimum and
//' maximum ages were specified, then these are also included. Also includes an 
//' empty string variable that can be used to describe stages meaningfully. This
//' object can be used as the \code{stageframe} input for \code{\link{flefko3}()} 
//' \code{\link{flefko2}()}, \code{\link{rlefko3}()}, and \code{\link{rlefko2}()}.
//' 
//' Variables in this data frame include the following:
//' \item{stage}{The unique names of the stages to be analyzed.}
//' \item{size}{The typical or representative size at which each stage occurs.}
//' \item{size_b}{Size at which each stage occurs in terms of a second size
//' variable, if one exists.}
//' \item{size_c}{Size at which each stage occurs in terms of a third size
//' variable, if one exists.}
//' \item{min_age}{The minimum age at which the stage may occur.}
//' \item{max_age}{The maximum age at which the stage may occur.}
//' \item{repstatus}{A binomial variable showing whether each stage is
//' reproductive.}
//' \item{obsstatus}{A binomial variable showing whether each stage is
//' observable.}
//' \item{propstatus}{A binomial variable showing whether each stage is a
//' propagule.}
//' \item{immstatus}{A binomial variable showing whether each stage can occur as
//' immature.}
//' \item{matstatus}{A binomial variable showing whether each stage occurs in
//' maturity.}
//' \item{indataset}{A binomial variable describing whether each stage occurs in
//' the input dataset.}
//' \item{binhalfwidth_raw}{The half-width of the size bin, as input.}
//' \item{sizebin_min}{The minimum size at which the stage may occur.}
//' \item{sizebin_max}{The maximum size at which the stage may occur.}
//' \item{sizebin_center}{The midpoint of the size bin at which the stage may
//' occur.}
//' \item{sizebin_width}{The width of the size bin corresponding to the stage.}
//' \item{binhalfwidthb_raw}{The half-width of the size bin of a second size
//' variable, as input.}
//' \item{sizebinb_min}{The minimum size at which the stage may occur.}
//' \item{sizebinb_max}{The maximum size at which the stage may occur.}
//' \item{sizebinb_center}{The midpoint of the size bin at which the stage may
//' occur, in terms of a second size variable.}
//' \item{sizebinb_width}{The width of the size bin corresponding to the stage,
//' in terms of a second size variable.}
//' \item{binhalfwidthc_raw}{The half-width of the size bin of a third size
//' variable, as input.}
//' \item{sizebinc_min}{The minimum size at which the stage may occur, in terms
//' of a third size variable.}
//' \item{sizebinc_max}{The maximum size at which the stage may occur, in terms
//' of a third size variable.}
//' \item{sizebinc_center}{The midpoint of the size bin at which the stage may
//' occur, in terms of a third size variable.}
//' \item{sizebinc_width}{The width of the size bin corresponding to the stage,
//' in terms of a third size variable.}
//' \item{group}{An integer denoting the size classification group that the
//' stage falls within.}
//' \item{comments}{A text field for stage descriptions.}
//' 
//' @section Notes:
//' If an IPM or function-based matrix with automated size classification is
//' desired, then two stages that occur within the dataset and represent the
//' lower and upper size limits of the IPM must be marked with \code{ipm} in
//' the stagenames vector. These stages should have all characteristics other
//' than size equal, and the size input for whichever size will be classified
//' automatically must include the minimum in one stage and the maximum in the
//' other. The actual characteristics of the first stage encountered in the
//' inputs will be used as the template for the creation of these sizes. Note
//' that \code{ipm} refers to size classification with the primary size
//' variable. To automate size classification with the secondary size variable,
//' use \code{ipmb}, and to automate size classification with the tertiary size
//' variable, use \code{ipmc}. To nest automated size classifications, use 
//' \code{ipmab} for the primary and secondary size variables, \code{ipmac} for
//' the primary and tertiary size variables, \code{ipmbc} for the secondary and
//' tertiary size variables, and \code{ipmabc} for all three size variables.
//' The primary size variable can also be set with \code{ipma}.
//' 
//' If two or more groups of stages, each with its own characteristics, are to
//' be developed for an IPM or function-based MPM, then an even number of stages
//' with two stages marking the minimum and maximum size of each group should be
//' marked with the same code as given above, with all other characteristics
//' equal within each group.
//' 
//' Stage classification groups set with the \code{group} variable create zones
//' within function-based matrices in which survival transitions are estimated.
//' These groups should not be set if transitions are possible between all
//' stages regardless of group. To denote specific transitions as estimable
//' between stage groups, use the \code{\link{supplemental}()} function.
//' 
//' @examples
//' # Lathyrus example
//' data(lathyrus)
//' 
//' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
//' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
//' repvector <- c(0, 0, 0, 0, 0, 1, 0)
//' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
//' matvector <- c(0, 0, 1, 1, 1, 1, 1)
//' immvector <- c(1, 1, 0, 0, 0, 0, 0)
//' propvector <- c(1, 0, 0, 0, 0, 0, 0)
//' indataset <- c(0, 1, 1, 1, 1, 1, 1)
//' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
//' 
//' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector,
//'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
//'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
//'   propstatus = propvector)
//' 
//' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
//'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9,
//'   juvcol = "Seedling1988", sizeacol = "Volume88", repstracol = "FCODE88",
//'   fecacol = "Intactseed88", deadacol = "Dead1988",
//'   nonobsacol = "Dormant1988", stageassign = lathframe, stagesize = "sizea",
//'   censorcol = "Missing1988", censorkeep = NA, censor = TRUE)
//' 
//' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl", "mat"),
//'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep", "Sdl"),
//'   stage1 = c("Sd", "rep", "Sd", "rep", "npr", "npr", "Sd"),
//'   eststage3 = c(NA, NA, NA, NA, NA, NA, "mat"),
//'   eststage2 = c(NA, NA, NA, NA, NA, NA, "Sdl"),
//'   eststage1 = c(NA, NA, NA, NA, NA, NA, "NotAlive"),
//'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
//'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054, NA),
//'   type = c(1, 1, 1, 1, 3, 3, 1), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
//'   stageframe = lathframe, historical = TRUE)
//' 
//' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = "all", 
//'   stages = c("stage3", "stage2", "stage1"), supplement = lathsupp3,
//'   yearcol = "year2", indivcol = "individ")
//' 
//' ehrlen3mean <- lmean(ehrlen3)
//' ehrlen3mean$A[[1]]
//' 
//' # Cypripedium example
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
//' cyp2mean
//' 
//' @export sf_create
// [[Rcpp::export]]
Rcpp::List sf_create (NumericVector sizes, Nullable<StringVector> stagenames = R_NilValue,
  Nullable<NumericVector> sizesb = R_NilValue, Nullable<NumericVector> sizesc = R_NilValue,
  Nullable<IntegerVector> repstatus = R_NilValue, Nullable<IntegerVector> obsstatus = R_NilValue,
  Nullable<IntegerVector> propstatus = R_NilValue, Nullable<IntegerVector> matstatus = R_NilValue,
  Nullable<IntegerVector> immstatus = R_NilValue, Nullable<NumericVector> minage = R_NilValue,
  Nullable<NumericVector> maxage = R_NilValue, Nullable<IntegerVector> indataset = R_NilValue,
  Nullable<NumericVector> binhalfwidth = R_NilValue, Nullable<NumericVector> binhalfwidthb = R_NilValue,
  Nullable<NumericVector> binhalfwidthc = R_NilValue, Nullable<IntegerVector> group = R_NilValue,
  Nullable<StringVector> comments = R_NilValue, int roundsize = 5, int roundsizeb = 5,
  int roundsizec = 5, int ipmbins = 100, int ipmbinsb = NA_INTEGER, int ipmbinsc = NA_INTEGER) {
  
  Rcpp::List output_longlist(29);
  Rcpp::CharacterVector varnames {"stage", "size", "size_b", "size_c", "min_age", "max_age",
    "repstatus", "obsstatus", "propstatus", "immstatus", "matstatus", "indataset",
    "binhalfwidth_raw", "sizebin_min", "sizebin_max", "sizebin_center", "sizebin_width",
    "binhalfwidthb_raw", "sizebinb_min", "sizebinb_max", "sizebinb_center", "sizebinb_width",
    "binhalfwidthc_raw", "sizebinc_min", "sizebinc_max", "sizebinc_center", "sizebinc_width",
    "group", "comments"};
  
  int matsize = sizes.size(); // Core vector size
  int used_sizes {1};
  arma::uvec ipm_calls_a (matsize, fill::zeros);
  arma::uvec ipm_calls_b (matsize, fill::zeros);
  arma::uvec ipm_calls_c (matsize, fill::zeros);
  arma::uvec ipm_calls_ab (matsize, fill::zeros);
  arma::uvec ipm_calls_bc (matsize, fill::zeros);
  arma::uvec ipm_calls_ac (matsize, fill::zeros);
  arma::uvec ipm_calls_abc (matsize, fill::zeros);
  
  StringVector stagenames_true (matsize, NA_STRING);
  NumericVector sizesb_true (matsize, NA_REAL);
  NumericVector sizesc_true (matsize, NA_REAL);
  IntegerVector repstatus_true (matsize, 1);
  IntegerVector obsstatus_true (matsize, 1);
  IntegerVector propstatus_true (matsize, 0);
  IntegerVector matstatus_true (matsize, 1);
  IntegerVector immstatus_true (matsize, 0);
  IntegerVector indataset_true (matsize, 1);
  NumericVector minage_true (matsize, NA_REAL);
  NumericVector maxage_true (matsize, NA_REAL);
  NumericVector binhalfwidth_true (matsize, 0.5);
  NumericVector binhalfwidthb_true (matsize, NA_REAL);
  NumericVector binhalfwidthc_true (matsize, NA_REAL);
  NumericVector sizebin_min (matsize, NA_REAL);
  NumericVector sizebin_max (matsize, NA_REAL);
  NumericVector sizebin_center (matsize, NA_REAL);
  NumericVector sizebin_width (matsize, NA_REAL);
  NumericVector sizebinb_min (matsize, NA_REAL);
  NumericVector sizebinb_max (matsize, NA_REAL);
  NumericVector sizebinb_center (matsize, NA_REAL);
  NumericVector sizebinb_width (matsize, NA_REAL);
  NumericVector sizebinc_min (matsize, NA_REAL);
  NumericVector sizebinc_max (matsize, NA_REAL);
  NumericVector sizebinc_center (matsize, NA_REAL);
  NumericVector sizebinc_width (matsize, NA_REAL);
  IntegerVector group_true (matsize, 0);
  StringVector comments_true (matsize, "No description");
  
  bool no_age = true;
  
  if (stagenames.isNotNull()) {
    Rcpp::StringVector stagenames_thru(stagenames);
    
    if (stagenames_thru.length() == matsize) {
      stagenames_true = stagenames_thru;
      
      StringVector st_t(stagenames_thru.size());
      
      std::transform(stagenames_thru.begin(), stagenames_thru.end(), st_t.begin(), 
        make_string_transformer(tolower));
      
      for (int i = 0; i < stagenames_thru.size(); i++) {
        String check_elem = trimws(st_t(i));
        
        if (check_elem == "ipm" || check_elem == "ipma" || check_elem == "ipm_a" ||
          check_elem == "ipm1" || check_elem == "ipm_1") {
          ipm_calls_a(i) = 1;
        } else if (check_elem == "ipmb" || check_elem == "ipm_b" || check_elem == "ipm2" ||
          check_elem == "ipm_2") {
          ipm_calls_b(i) = 1;
        } else if (check_elem == "ipmc" || check_elem == "ipm_c" || check_elem == "ipm3" ||
          check_elem == "ipm_3") {
          ipm_calls_c(i) = 1;
        } else if (check_elem == "ipmab" || check_elem == "ipm_ab" || check_elem == "ipm12" ||
          check_elem == "ipm_12") {
          ipm_calls_ab(i) = 1;
        } else if (check_elem == "ipmac" || check_elem == "ipm_ac" || check_elem == "ipm13" ||
          check_elem == "ipm_13") {
          ipm_calls_ac(i) = 1;
        } else if (check_elem == "ipmbc" || check_elem == "ipm_bc" || check_elem == "ipm23" ||
          check_elem == "ipm_23") {
          ipm_calls_bc(i) = 1;
        } else if (check_elem == "ipmabc" || check_elem == "ipm_abc" || check_elem == "ipm123" ||
          check_elem == "ipm_123") {
          ipm_calls_abc(i) = 1;
        }
      }
    } else if (stagenames_thru.length() == 1) {
      
      if (!stagenames_thru.is_na(0)) {
        for (int i = 0; i < matsize; i++) {
          Rcpp::String part1(stagenames_thru(0));
          part1 += (static_cast<char>(i + 1));
          
          stagenames_true(i) = part1;
        }
      } else {
        for (int i = 0; i < matsize; i++) {
          Rcpp::String part1("Stage");
          part1 += (static_cast<char>(i + 1));
          
          stagenames_true(i) = part1;
        }
      }
    } else {
      throw Rcpp::exception("Vector stagenames should be the same length as vector sizes.", false);
    }
  } else {
    for (int i = 0; i < matsize; i++) {
      Rcpp::String part1("Stage");
      part1 += (static_cast<char>(i + 1));
      
      stagenames_true(i) = part1;
    }
  }
  
  if (sizesb.isNotNull()) {
    Rcpp::NumericVector sizesb_thru(sizesb);
    used_sizes++;
    
    if (sizesb_thru.length() == matsize) {
      sizesb_true = sizesb_thru;
    } else if (sizesb_thru.length() == 1) {
      NumericVector try_size (matsize, sizesb_thru(0));
      sizesb_true = try_size;
    } else {
      throw Rcpp::exception("Vector sizesb should be the same length as vector sizes.", false);
    }
  }
  if (sizesc.isNotNull()) {
    Rcpp::NumericVector sizesc_thru(sizesc);
    used_sizes++;
    
    if (sizesc_thru.length() == matsize) {
      sizesc_true = sizesc_thru;
    } else if (sizesc_thru.length() == 1) {
      NumericVector try_size (matsize, sizesc_thru(0));
      sizesc_true = try_size;
    } else {
      throw Rcpp::exception("Vector sizesc should be the same length as vector sizes.", false);
    }
    
    if (used_sizes != 3) {
      throw Rcpp::exception("Vector sizesc should only be set if vector sizesb is also set.", false);
    }
  }
  if (minage.isNotNull()) {
    Rcpp::NumericVector minage_thru(minage);
    
    if (minage_thru.length() == matsize) {
      minage_true = minage_thru;
    } else if (minage_thru.length() == 1) {
      NumericVector try_mna (matsize, minage_thru(0));
      minage_true = try_mna;
    } else {
      throw Rcpp::exception("Vector minage should be the same length as vector sizes.", false);
    }
    
    no_age = false;
  }
  if (maxage.isNotNull()) {
    Rcpp::NumericVector maxage_thru(maxage);
    
    if (maxage_thru.length() == matsize) {
      maxage_true = maxage_thru;
    } else if (maxage_thru.length() == 1) {
      NumericVector try_mxa (matsize, maxage_thru(0));
      maxage_true = try_mxa;
    } else {
      throw Rcpp::exception("Vector maxage should be the same length as vector sizes.", false);
    }
    if (no_age) {
      throw Rcpp::exception("Vector minage is required if vector maxage is provided.", false);
    }
  } else if (!no_age) {
    throw Rcpp::exception("Vector maxage is required if vector minage is provided.", false);
  }
  
  if (repstatus.isNotNull()) {
    Rcpp::IntegerVector repstatus_thru(repstatus);
    
    if (repstatus_thru.length() == matsize) {
      repstatus_true = repstatus_thru;
    } else if (repstatus_thru.length() == 1) {
      IntegerVector try_rep (matsize, repstatus_thru(0));
      repstatus_true = try_rep;
    } else {
      throw Rcpp::exception("Vector repstatus should be the same length as vector sizes.", false);
    }
    
    if (max(repstatus_true) > 1 || min(repstatus_true) < 0) {
      throw Rcpp::exception("Vector repstatus should be composed only of 0s and 1s.", false);
    }
  }
  if (obsstatus.isNotNull()) {
    Rcpp::IntegerVector obsstatus_thru(obsstatus);
    
    if (obsstatus_thru.length() == matsize) {
      obsstatus_true = obsstatus_thru;
    } else if (obsstatus_thru.length() == 1) {
      IntegerVector try_obs (matsize, obsstatus_thru(0));
      obsstatus_true = try_obs;
    } else {
      throw Rcpp::exception("Vector obsstatus should be the same length as vector sizes.", false);
    }
    
    if (max(obsstatus_true) > 1 || min(obsstatus_true) < 0) {
      throw Rcpp::exception("Vector obsstatus should be composed only of 0s and 1s.", false);
    }
  }
  if (propstatus.isNotNull()) {
    Rcpp::IntegerVector propstatus_thru(propstatus);
    
    if (propstatus_thru.length() == matsize) {
      propstatus_true = propstatus_thru;
    } else if (propstatus_thru.length() == 1) {
      IntegerVector try_prop (matsize, propstatus_thru(0));
      propstatus_true = try_prop;
    } else {
      throw Rcpp::exception("Vector propstatus should be the same length as vector sizes.", false);
    }
    
    if (max(propstatus_true) > 1 || min(propstatus_true) < 0) {
      throw Rcpp::exception("Vector propstatus should be composed only of 0s and 1s.", false);
    }
  }
  if (matstatus.isNotNull()) {
    Rcpp::IntegerVector matstatus_thru(matstatus);
    
    if (matstatus_thru.length() == matsize) {
      matstatus_true = matstatus_thru;
    } else if (matstatus_thru.length() == 1) {
      IntegerVector try_mat (matsize, matstatus_thru(0));
      matstatus_true = try_mat;
    } else {
      throw Rcpp::exception("Vector matstatus should be the same length as vector sizes.", false);
    }
    
    if (max(matstatus_true) > 1 || min(matstatus_true) < 0) {
      throw Rcpp::exception("Vector matstatus should be composed only of 0s and 1s.", false);
    }
  }
  if (immstatus.isNotNull()) {
    Rcpp::IntegerVector immstatus_thru(immstatus);
    
    if (immstatus_thru.length() == matsize) {
      immstatus_true = immstatus_thru;
    } else if (immstatus_thru.length() == 1) {
      IntegerVector try_imm (matsize, immstatus_thru(0));
      immstatus_true = try_imm;
    } else {
      throw Rcpp::exception("Vector immstatus should be the same length as vector sizes.", false);
    }
    
    if (max(immstatus_true) > 1 || min(immstatus_true) < 0) {
      throw Rcpp::exception("Vector immstatus should be composed only of 0s and 1s.", false);
    }
  } else {
    IntegerVector trialones(matsize, 1);
    immstatus_true = trialones - matstatus_true;
  }
  if (indataset.isNotNull()) {
    Rcpp::IntegerVector indataset_thru(indataset);
    
    if (indataset_thru.length() == matsize) {
      indataset_true = indataset_thru;
    } else if (indataset_thru.length() == 1) {
      IntegerVector try_ind (matsize, indataset_thru(0));
      indataset_true = try_ind;
    } else {
      throw Rcpp::exception("Vector indataset should be the same length as vector sizes.", false);
    }
    
    if (max(indataset_true) > 1 || min(indataset_true) < 0) {
      throw Rcpp::exception("Vector indataset should be composed only of 0s and 1s.", false);
    }
  }
  
  if (binhalfwidth.isNotNull()) {
    Rcpp::NumericVector binhalfwidth_thru(binhalfwidth);
    
    if (binhalfwidth_thru.length() == matsize) {
      binhalfwidth_true = binhalfwidth_thru;
    } else if (binhalfwidth_thru.length() == 1) {
      NumericVector try_hwa (matsize, binhalfwidth_thru(0));
      binhalfwidth_true = try_hwa;
    } else {
      throw Rcpp::exception("Vector binhalfwidth should be the same length as vector sizes.", false);
    }
  }
  if (binhalfwidthb.isNotNull()) {
    if (used_sizes == 1) {
      throw Rcpp::exception("Vector binhalfwidthb should only be used if multiple size variables are being used for classification.", false);
    }
    
    Rcpp::NumericVector binhalfwidthb_thru(binhalfwidthb);
    
    if (binhalfwidthb_thru.length() == matsize) {
      binhalfwidthb_true = binhalfwidthb_thru;
    } else if (binhalfwidthb_thru.length() == 1) {
      NumericVector try_hwb (matsize, binhalfwidthb_thru(0));
      binhalfwidthb_true = try_hwb;
    } else {
      throw Rcpp::exception("Vector binhalfwidthb should be the same length as vector sizes.", false);
    }
  } else if (used_sizes > 1) {
    for (int i = 0; i < binhalfwidthb_true.length(); i++) {
      binhalfwidthb_true(i) = 0.5;
    }
  }
  if (binhalfwidthc.isNotNull()) {
    if (used_sizes < 3) {
      throw Rcpp::exception("Vector binhalfwidthc should only be used if three size variables are being used for classification.", false);
    }
    
    Rcpp::NumericVector binhalfwidthc_thru(binhalfwidthc);
    
    if (binhalfwidthc_thru.length() == matsize) {
      binhalfwidthc_true = binhalfwidthc_thru;
    } else if (binhalfwidthc_thru.length() == 1) {
      NumericVector try_hwc (matsize, binhalfwidthc_thru(0));
      binhalfwidthc_true = try_hwc;
    } else {
      throw Rcpp::exception("Vector binhalfwidthc should be the same length as vector sizes.", false);
    }
  } else if (used_sizes > 2) {
    for (int i = 0; i < binhalfwidthc_true.length(); i++) {
      binhalfwidthc_true(i) = 0.5;
    }
  }
  
  if (group.isNotNull()) {
    Rcpp::IntegerVector group_thru(group);
    
    if (group_thru.length() == matsize) {
      group_true = group_thru;
    } else if (group_thru.length() == 1) {
      IntegerVector try_grp (matsize, group_thru(0));
      group_true = try_grp;
    } else {
      throw Rcpp::exception("Vector group should be the same length as vector sizes.", false);
    }
    
    if (min(group_true) < 0) {
      throw Rcpp::exception("Please use only positive integers for group designations.", false);
    }
  }
  
  if (comments.isNotNull()) {
    Rcpp::StringVector comments_thru(comments);
    
    if (comments_thru.length() == matsize) {
      comments_true = comments_thru;
    } else if (comments_thru.length() == 1) {
      StringVector try_com (matsize, comments_thru(0));
      comments_true = try_com;
    } else {
      throw Rcpp::exception("Comments vector should be the same length as vector sizes.", false);
    }
  }
  
  
  if (ipmbins < 2) {
    throw Rcpp::exception("Please enter a valid integer greater than 1 for ipmbins option.", false);
  } else if (ipmbins > 100) {
    Rf_warningcall(R_NilValue, "High ipmbin numbers may lead to dramatic decreases in statistical power and overparameterized matrices.");
  }
  
  // Now the automated size classification processing
  int ipm_a = sum(ipm_calls_a);
  int ipm_b = sum(ipm_calls_b);
  int ipm_c = sum(ipm_calls_c);
  int ipm_ab = sum(ipm_calls_ab);
  int ipm_bc = sum(ipm_calls_bc);
  int ipm_ac = sum(ipm_calls_ac);
  int ipm_abc = sum(ipm_calls_abc);
  
  // This vector will act as an unique marker for most characteristics
  arma::ivec pair_check (matsize, fill::zeros); 
  arma::uvec entries_to_delete (1);
  int no_entries_to_delete {0};
  
  for (int i = 0; i < matsize; i++) {
    pair_check(i) = repstatus_true(i) * 1000000 + obsstatus_true(i) * 100000 + 
      propstatus_true(i) * 10000 + immstatus_true(i) * 1000 + matstatus_true(i) * 100 +
      indataset_true(i) * 10 + group_true(i);
  }

  // Automated size classification with sizea
  if (ipm_a > 0) {
    if (IntegerVector::is_na(ipmbins)) {
      throw Rcpp::exception("The number of bins for automated size classification of the primary size variable has not been set.", false);
    }
    
    if (ipm_a % 2 != 0) {
      throw Rcpp::exception("The ipm designation must specify both the start size and the end size, requiring an even number of calls. Calls for automated size classification must be matched and not overlap.", false);
    }
    
    arma::uvec check_elems = find(ipm_calls_a); // This vector points out which rows have ipm designations
    
    arma::ivec called_pair_check = pair_check.elem(check_elems); // 
    int called_pair_check_length = called_pair_check.n_elem;
    
    for (int i = 0; i < called_pair_check_length; i++) {
      int trial_1 = called_pair_check(i);
      int match_count {0};
      int main_index_1 {i};
      int go_ahead {0};
      
      for (int j = 0; j < called_pair_check_length; j++) {
        if (trial_1 == called_pair_check(j) && i != j) {
          match_count++;
          main_index_1 = j;
          
          if (j > i) {
            go_ahead = 1;
          } else {
            go_ahead = 0;
          }
        }
      }
      
      if (go_ahead == 1) {
        if (match_count > 1) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. More than 2 stages with the same characteristics marked ipm cannot be handled.", false);
        } else if (match_count == 0) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. Single stages with unique characteristics marked ipm cannot be handled.", false);
        }
        
        // Now we work out the minimum and maximum sizes
        double minsize {0};
        double maxsize {0};
        
        if (sizes(check_elems(i)) > sizes(check_elems(main_index_1))) {
          maxsize = sizes(check_elems(i));
          minsize = sizes(check_elems(main_index_1));
        } else if (sizes(check_elems(i)) < sizes(check_elems(main_index_1))) {
          minsize = sizes(check_elems(i));
          maxsize = sizes(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        
        // Keep track of the entries to delete from the original vectors
        if (no_entries_to_delete == 0) {
          entries_to_delete.resize(2);
          
          entries_to_delete(0) = check_elems(i);
          entries_to_delete(1) = check_elems(main_index_1);
        } else {
          int entries_del_old_size = entries_to_delete.size();
          entries_to_delete.resize(entries_del_old_size + 2);
          
          entries_to_delete(entries_del_old_size) = check_elems(i);
          entries_to_delete(entries_del_old_size + 1) = check_elems(main_index_1);
        }
        
        no_entries_to_delete++;
        no_entries_to_delete++;
        
        // Now we will create the new stages
        double full_range = maxsize - minsize;
        double standard_increment = full_range / static_cast<double>(ipmbins);
        double standard_midpoint = standard_increment / 2;
        
        // New vectors to append
        NumericVector newsizes (ipmbins, 0.0);
        StringVector newstagenames (ipmbins, "");
        NumericVector newsizesb (ipmbins, 0.0);
        NumericVector newsizesc (ipmbins, 0.0);
        NumericVector newminage (ipmbins, 0.0);
        NumericVector newmaxage (ipmbins, 0.0);
        IntegerVector newrepstatus (ipmbins, 0);
        IntegerVector newobsstatus (ipmbins, 0);
        IntegerVector newpropstatus (ipmbins, 0);
        IntegerVector newmatstatus (ipmbins, 0);
        IntegerVector newimmstatus (ipmbins, 0);
        IntegerVector newindataset (ipmbins, 0);
        NumericVector newbinhalfwidth (ipmbins, 0.0);
        NumericVector newbinhalfwidthb (ipmbins, 0.0);
        NumericVector newbinhalfwidthc (ipmbins, 0.0);
        IntegerVector newgroup (ipmbins, 0);
        StringVector newcomments (ipmbins, "");
        
        for (int j = 0; j < ipmbins; j++) {
          newsizes(j) = minsize + (j*standard_increment) + standard_midpoint;
          newsizesb(j) = sizesb_true(check_elems(i));
          newsizesc(j) = sizesc_true(check_elems(i));
          newminage(j) = minage_true(check_elems(i));
          newmaxage(j) = maxage_true(check_elems(i));
          newrepstatus(j) = repstatus_true(check_elems(i));
          newobsstatus(j) = obsstatus_true(check_elems(i));
          newpropstatus(j) = propstatus_true(check_elems(i));
          newmatstatus(j) = matstatus_true(check_elems(i));
          newimmstatus(j) = immstatus_true(check_elems(i));
          newindataset(j) = indataset_true(check_elems(i));
          newbinhalfwidth(j) = standard_midpoint;
          newbinhalfwidthb(j) = binhalfwidthb_true(check_elems(i));
          newbinhalfwidthc(j) = binhalfwidthc_true(check_elems(i));
          newgroup(j) = group_true(check_elems(i));
          newcomments(j) = comments_true(check_elems(i));
          
          std::string sizenums = std::to_string(newsizes(j));
          newstagenames(j) = "sza_" + sizenums.substr(0, 6);
          newstagenames(j) += "_";
          newstagenames(j) += std::to_string(newgroup(j));
        }
        
        sizes = concat_dbl(sizes, newsizes);
        sizesb_true = concat_dbl(sizesb_true, newsizesb);
        sizesc_true = concat_dbl(sizesc_true, newsizesc);
        minage_true = concat_dbl(minage_true, newminage);
        maxage_true = concat_dbl(maxage_true, newmaxage);
        repstatus_true = concat_int(repstatus_true, newrepstatus);
        obsstatus_true = concat_int(obsstatus_true, newobsstatus);
        propstatus_true = concat_int(propstatus_true, newpropstatus);
        matstatus_true = concat_int(matstatus_true, newmatstatus);
        immstatus_true = concat_int(immstatus_true, newimmstatus);
        indataset_true = concat_int(indataset_true, newindataset);
        binhalfwidth_true = concat_dbl(binhalfwidth_true, newbinhalfwidth);
        binhalfwidthb_true = concat_dbl(binhalfwidthb_true, newbinhalfwidthb);
        binhalfwidthc_true = concat_dbl(binhalfwidthc_true, newbinhalfwidthc);
        stagenames_true = concat_str(stagenames_true, newstagenames);
        group_true = concat_int(group_true, newgroup);
        comments_true = concat_str(comments_true, newcomments);
      }
    }
  }
  
  // IPM classification with sizeb
  if (ipm_b > 0) {
    if (IntegerVector::is_na(ipmbinsb)) {
      throw Rcpp::exception("The number of bins for automated size classification of the secondary size variable has not been set.", false);
    }
    
    if (ipm_b % 2 != 0) {
      throw Rcpp::exception("The ipm designation must specify both the start size and the end size, requiring an even number of calls. Calls for automated size classification must be matched and not overlap.", false);
    }
    
    arma::uvec check_elems = find(ipm_calls_b); // This vector points out which rows have ipm designations
    
    arma::ivec called_pair_check = pair_check.elem(check_elems); // 
    int called_pair_check_length = called_pair_check.n_elem;
    
    for (int i = 0; i < called_pair_check_length; i++) {
      int trial_1 = called_pair_check(i);
      int match_count {0};
      int main_index_1 {i};
      int go_ahead {0};
      
      for (int j = 0; j < called_pair_check_length; j++) {
        if (trial_1 == called_pair_check(j) && i != j) {
          match_count++;
          main_index_1 = j;
          
          if (j > i) {
            go_ahead = 1;
          } else {
            go_ahead = 0;
          }
        }
      }

      if (go_ahead == 1) {
        
        if (match_count > 1) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. More than 2 stages with the same characteristics marked ipm cannot be handled.", false);
        } else if (match_count == 0) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. Single stages with unique characteristics marked ipm cannot be handled.", false);
        }
        
        // Now we work out the minimum and maximum sizes
        double minsize {0};
        double maxsize {0};
        
        if (sizesb_true(check_elems(i)) > sizesb_true(check_elems(main_index_1))) {
          maxsize = sizesb_true(check_elems(i));
          minsize = sizesb_true(check_elems(main_index_1));
        } else if (sizesb_true(check_elems(i)) < sizesb_true(check_elems(main_index_1))) {
          minsize = sizesb_true(check_elems(i));
          maxsize = sizesb_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        
        // Keep track of the entries to delete from the original vectors
        if (no_entries_to_delete == 0) {
          entries_to_delete.resize(2);
          
          entries_to_delete(0) = check_elems(i);
          entries_to_delete(1) = check_elems(main_index_1);
        } else {
          int entries_del_old_size = entries_to_delete.size();
          entries_to_delete.resize(entries_del_old_size + 2);
          
          entries_to_delete(entries_del_old_size) = check_elems(i);
          entries_to_delete(entries_del_old_size + 1) = check_elems(main_index_1);
        }
        
        no_entries_to_delete++;
        no_entries_to_delete++;

        // Now we will create the new stages
        double full_range = maxsize - minsize;
        double standard_increment = full_range / static_cast<double>(ipmbins);
        double standard_midpoint = standard_increment / 2;
        
        // New vectors to append
        NumericVector newsizes (ipmbinsb, 0.0);
        StringVector newstagenames (ipmbinsb, "");
        NumericVector newsizesb (ipmbinsb, 0.0);
        NumericVector newsizesc (ipmbinsb, 0.0);
        NumericVector newminage (ipmbinsb, 0.0);
        NumericVector newmaxage (ipmbinsb, 0.0);
        IntegerVector newrepstatus (ipmbinsb, 0);
        IntegerVector newobsstatus (ipmbinsb, 0);
        IntegerVector newpropstatus (ipmbinsb, 0);
        IntegerVector newmatstatus (ipmbinsb, 0);
        IntegerVector newimmstatus (ipmbinsb, 0);
        IntegerVector newindataset (ipmbinsb, 0);
        NumericVector newbinhalfwidth (ipmbinsb, 0.0);
        NumericVector newbinhalfwidthb (ipmbinsb, 0.0);
        NumericVector newbinhalfwidthc (ipmbinsb, 0.0);
        IntegerVector newgroup (ipmbinsb, 0);
        StringVector newcomments (ipmbinsb, "");
        
        for (int j = 0; j < ipmbinsb; j++) {
          newsizesb(j) = minsize + (j*standard_increment) + standard_midpoint;
          newsizes(j) = sizes(check_elems(i));
          newsizesc(j) = sizesc_true(check_elems(i));
          newminage(j) = minage_true(check_elems(i));
          newmaxage(j) = maxage_true(check_elems(i));
          newrepstatus(j) = repstatus_true(check_elems(i));
          newobsstatus(j) = obsstatus_true(check_elems(i));
          newpropstatus(j) = propstatus_true(check_elems(i));
          newmatstatus(j) = matstatus_true(check_elems(i));
          newimmstatus(j) = immstatus_true(check_elems(i));
          newindataset(j) = indataset_true(check_elems(i));
          newbinhalfwidth(j) = standard_midpoint;
          newbinhalfwidthb(j) = binhalfwidthb_true(check_elems(i));
          newbinhalfwidthc(j) = binhalfwidthc_true(check_elems(i));
          newgroup(j) = group_true(check_elems(i));
          newcomments(j) = comments_true(check_elems(i));
          
          std::string sizenums = std::to_string(newsizesb(j));
          newstagenames(j) = "szb_" + sizenums.substr(0, 6);
          newstagenames(j) += "_";
          newstagenames(j) += std::to_string(newgroup(j));
        }
        
        sizes = concat_dbl(sizes, newsizes);
        sizesb_true = concat_dbl(sizesb_true, newsizesb);
        sizesc_true = concat_dbl(sizesc_true, newsizesc);
        minage_true = concat_dbl(minage_true, newminage);
        maxage_true = concat_dbl(maxage_true, newmaxage);
        repstatus_true = concat_int(repstatus_true, newrepstatus);
        obsstatus_true = concat_int(obsstatus_true, newobsstatus);
        propstatus_true = concat_int(propstatus_true, newpropstatus);
        matstatus_true = concat_int(matstatus_true, newmatstatus);
        immstatus_true = concat_int(immstatus_true, newimmstatus);
        indataset_true = concat_int(indataset_true, newindataset);
        binhalfwidth_true = concat_dbl(binhalfwidth_true, newbinhalfwidth);
        binhalfwidthb_true = concat_dbl(binhalfwidthb_true, newbinhalfwidthb);
        binhalfwidthc_true = concat_dbl(binhalfwidthc_true, newbinhalfwidthc);
        stagenames_true = concat_str(stagenames_true, newstagenames);
        group_true = concat_int(group_true, newgroup);
        comments_true = concat_str(comments_true, newcomments);
      }
    }
  }
  
  // IPM classification with sizec
  if (ipm_c > 0) {
    if (IntegerVector::is_na(ipmbinsc)) {
      throw Rcpp::exception("The number of bins for automated size classification of the tertiary size variable has not been set.", false);
    }
    
    if (ipm_c % 2 != 0) {
      throw Rcpp::exception("The ipm designation must specify both the start size and the end size, requiring an even number of calls. Calls for automated size classification must be matched and not overlap.", false);
    }
    
    arma::uvec check_elems = find(ipm_calls_c); // This vector points out which rows have ipm designations
    
    arma::ivec called_pair_check = pair_check.elem(check_elems); // 
    int called_pair_check_length = called_pair_check.n_elem;
    
    for (int i = 0; i < called_pair_check_length; i++) {
      int trial_1 = called_pair_check(i);
      int match_count {0};
      int main_index_1 {i};
      int go_ahead {0};
      
      for (int j = 0; j < called_pair_check_length; j++) {
        if (trial_1 == called_pair_check(j) && i != j) {
          match_count++;
          main_index_1 = j;
          
          if (j > i) {
            go_ahead = 1;
          } else {
            go_ahead = 0;
          }
        }
      }
      
      if (go_ahead == 1) {
        
        if (match_count > 1) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. More than 2 stages with the same characteristics marked ipm cannot be handled.", false);
        } else if (match_count == 0) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. Single stages with unique characteristics marked ipm cannot be handled.", false);
        }
        
        // Now we work out the minimum and maximum sizes
        double minsize {0};
        double maxsize {0};
        
        if (sizesc_true(check_elems(i)) > sizesc_true(check_elems(main_index_1))) {
          maxsize = sizesc_true(check_elems(i));
          minsize = sizesc_true(check_elems(main_index_1));
        } else if (sizesc_true(check_elems(i)) < sizesc_true(check_elems(main_index_1))) {
          minsize = sizesc_true(check_elems(i));
          maxsize = sizesc_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        
        // Keep track of the entries to delete from the original vectors
        if (no_entries_to_delete == 0) {
          entries_to_delete.resize(2);
          
          entries_to_delete(0) = check_elems(i);
          entries_to_delete(1) = check_elems(main_index_1);
        } else {
          int entries_del_old_size = entries_to_delete.size();
          entries_to_delete.resize(entries_del_old_size + 2);
          
          entries_to_delete(entries_del_old_size) = check_elems(i);
          entries_to_delete(entries_del_old_size + 1) = check_elems(main_index_1);
        }
        
        no_entries_to_delete++;
        no_entries_to_delete++;
        
        // Now we will create the new stages
        double full_range = maxsize - minsize;
        double standard_increment = full_range / static_cast<double>(ipmbins);
        double standard_midpoint = standard_increment / 2;
        
        // New vectors to append
        NumericVector newsizes (ipmbinsc, 0.0);
        StringVector newstagenames (ipmbinsc, "");
        NumericVector newsizesb (ipmbinsc, 0.0);
        NumericVector newsizesc (ipmbinsc, 0.0);
        NumericVector newminage (ipmbinsc, 0.0);
        NumericVector newmaxage (ipmbinsc, 0.0);
        IntegerVector newrepstatus (ipmbinsc, 0);
        IntegerVector newobsstatus (ipmbinsc, 0);
        IntegerVector newpropstatus (ipmbinsc, 0);
        IntegerVector newmatstatus (ipmbinsc, 0);
        IntegerVector newimmstatus (ipmbinsc, 0);
        IntegerVector newindataset (ipmbinsc, 0);
        NumericVector newbinhalfwidth (ipmbinsc, 0.0);
        NumericVector newbinhalfwidthb (ipmbinsc, 0.0);
        NumericVector newbinhalfwidthc (ipmbinsc, 0.0);
        IntegerVector newgroup (ipmbinsc, 0);
        StringVector newcomments (ipmbinsc, "");
        
        for (int j = 0; j < ipmbinsc; j++) {
          newsizesc(j) = minsize + (j*standard_increment) + standard_midpoint;
          newsizesb(j) = sizesb_true(check_elems(i));
          newsizes(j) = sizes(check_elems(i));
          newminage(j) = minage_true(check_elems(i));
          newmaxage(j) = maxage_true(check_elems(i));
          newrepstatus(j) = repstatus_true(check_elems(i));
          newobsstatus(j) = obsstatus_true(check_elems(i));
          newpropstatus(j) = propstatus_true(check_elems(i));
          newmatstatus(j) = matstatus_true(check_elems(i));
          newimmstatus(j) = immstatus_true(check_elems(i));
          newindataset(j) = indataset_true(check_elems(i));
          newbinhalfwidth(j) = standard_midpoint;
          newbinhalfwidthb(j) = binhalfwidthb_true(check_elems(i));
          newbinhalfwidthc(j) = binhalfwidthc_true(check_elems(i));
          newgroup(j) = group_true(check_elems(i));
          newcomments(j) = comments_true(check_elems(i));
          
          std::string sizenums = std::to_string(newsizesc(j));
          newstagenames(j) = "szc_" + sizenums.substr(0, 6);
          newstagenames(j) += "_";
          newstagenames(j) += std::to_string(newgroup(j));
        }
        
        sizes = concat_dbl(sizes, newsizes);
        sizesb_true = concat_dbl(sizesb_true, newsizesb);
        sizesc_true = concat_dbl(sizesc_true, newsizesc);
        minage_true = concat_dbl(minage_true, newminage);
        maxage_true = concat_dbl(maxage_true, newmaxage);
        repstatus_true = concat_int(repstatus_true, newrepstatus);
        obsstatus_true = concat_int(obsstatus_true, newobsstatus);
        propstatus_true = concat_int(propstatus_true, newpropstatus);
        matstatus_true = concat_int(matstatus_true, newmatstatus);
        immstatus_true = concat_int(immstatus_true, newimmstatus);
        indataset_true = concat_int(indataset_true, newindataset);
        binhalfwidth_true = concat_dbl(binhalfwidth_true, newbinhalfwidth);
        binhalfwidthb_true = concat_dbl(binhalfwidthb_true, newbinhalfwidthb);
        binhalfwidthc_true = concat_dbl(binhalfwidthc_true, newbinhalfwidthc);
        stagenames_true = concat_str(stagenames_true, newstagenames);
        group_true = concat_int(group_true, newgroup);
        comments_true = concat_str(comments_true, newcomments);
      }
    }
  }
  
  // Automated size classification with sizea & sizeb together
  if (ipm_ab > 0) {
    if (IntegerVector::is_na(ipmbins)) {
      throw Rcpp::exception("The number of bins for automated size classification of the primary size variable has not been set.", false);
    }
    
    if (IntegerVector::is_na(ipmbinsb)) {
      throw Rcpp::exception("The number of bins for automated size classification of the primary size variable has not been set.", false);
    }
    
    if (ipm_ab % 2 != 0) {
      throw Rcpp::exception("The ipm designation must specify both the start size and the end size, requiring an even number of calls. Calls for automated size classification must be matched and not overlap.", false);
    }
    
    arma::uvec check_elems = find(ipm_calls_ab); // This vector points out which rows have ipm designations
    
    arma::ivec called_pair_check = pair_check.elem(check_elems); // 
    int called_pair_check_length = called_pair_check.n_elem;
    
    for (int i = 0; i < called_pair_check_length; i++) {
      int trial_1 = called_pair_check(i);
      int match_count {0};
      int main_index_1 {i};
      int go_ahead {0};
      
      for (int j = 0; j < called_pair_check_length; j++) {
        if (trial_1 == called_pair_check(j) && i != j) {
          match_count++;
          main_index_1 = j;
          
          if (j > i) {
            go_ahead = 1;
          } else {
            go_ahead = 0;
          }
        }
      }
      
      if (go_ahead == 1) {
        
        if (match_count > 1) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. More than 2 stages with the same characteristics marked ipm cannot be handled.", false);
        } else if (match_count == 0) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. Single stages with unique characteristics marked ipm cannot be handled.", false);
        }
        
        // Now we work out the minimum and maximum sizes
        double minsize_a {0};
        double maxsize_a {0};
        double minsize_b {0};
        double maxsize_b {0};
        
        if (sizes(check_elems(i)) > sizes(check_elems(main_index_1))) {
          maxsize_a = sizes(check_elems(i));
          minsize_a = sizes(check_elems(main_index_1));
        } else if (sizes(check_elems(i)) < sizes(check_elems(main_index_1))) {
          minsize_a = sizes(check_elems(i));
          maxsize_a = sizes(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        if (sizesb_true(check_elems(i)) > sizesb_true(check_elems(main_index_1))) {
          maxsize_b = sizesb_true(check_elems(i));
          minsize_b = sizesb_true(check_elems(main_index_1));
        } else if (sizesb_true(check_elems(i)) < sizesb_true(check_elems(main_index_1))) {
          minsize_b = sizesb_true(check_elems(i));
          maxsize_b = sizesb_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        
        // Keep track of the entries to delete from the original vectors
        if (no_entries_to_delete == 0) {
          entries_to_delete.resize(2);
          
          entries_to_delete(0) = check_elems(i);
          entries_to_delete(1) = check_elems(main_index_1);
        } else {
          int entries_del_old_size = entries_to_delete.size();
          entries_to_delete.resize(entries_del_old_size + 2);
          
          entries_to_delete(entries_del_old_size) = check_elems(i);
          entries_to_delete(entries_del_old_size + 1) = check_elems(main_index_1);
        }
        
        no_entries_to_delete++;
        no_entries_to_delete++;
        
        // Now we will create the new stages
        double full_range_a = maxsize_a - minsize_a;
        double standard_increment_a = full_range_a / static_cast<double>(ipmbins);
        double standard_midpoint_a = standard_increment_a / 2;
        
        double full_range_b = maxsize_b - minsize_b;
        double standard_increment_b = full_range_b / static_cast<double>(ipmbinsb);
        double standard_midpoint_b = standard_increment_b / 2;
        
        // New vectors to append
        NumericVector newsizes ((ipmbins * ipmbinsb), 0.0);
        StringVector newstagenames ((ipmbins * ipmbinsb), "");
        NumericVector newsizesb ((ipmbins * ipmbinsb), 0.0);
        NumericVector newsizesc ((ipmbins * ipmbinsb), 0.0);
        NumericVector newminage ((ipmbins * ipmbinsb), 0.0);
        NumericVector newmaxage ((ipmbins * ipmbinsb), 0.0);
        IntegerVector newrepstatus ((ipmbins * ipmbinsb), 0);
        IntegerVector newobsstatus ((ipmbins * ipmbinsb), 0);
        IntegerVector newpropstatus ((ipmbins * ipmbinsb), 0);
        IntegerVector newmatstatus ((ipmbins * ipmbinsb), 0);
        IntegerVector newimmstatus ((ipmbins * ipmbinsb), 0);
        IntegerVector newindataset ((ipmbins * ipmbinsb), 0);
        NumericVector newbinhalfwidth ((ipmbins * ipmbinsb), 0.0);
        NumericVector newbinhalfwidthb ((ipmbins * ipmbinsb), 0.0);
        NumericVector newbinhalfwidthc ((ipmbins * ipmbinsb), 0.0);
        IntegerVector newgroup ((ipmbins * ipmbinsb), 0);
        StringVector newcomments ((ipmbins * ipmbinsb), "");
        
        for (int j = 0; j < ipmbins; j++) {
          for (int k = 0; k < ipmbinsb; k++) {
            newsizes((j * ipmbinsb) + k) = minsize_a + (j*standard_increment_a) + standard_midpoint_a;
            newsizesb((j * ipmbinsb) + k) = minsize_b + (k*standard_increment_b) + standard_midpoint_b;;
            newsizesc((j * ipmbinsb) + k) = sizesc_true(check_elems(i));
            newminage((j * ipmbinsb) + k) = minage_true(check_elems(i));
            newmaxage((j * ipmbinsb) + k) = maxage_true(check_elems(i));
            newrepstatus((j * ipmbinsb) + k) = repstatus_true(check_elems(i));
            newobsstatus((j * ipmbinsb) + k) = obsstatus_true(check_elems(i));
            newpropstatus((j * ipmbinsb) + k) = propstatus_true(check_elems(i));
            newmatstatus((j * ipmbinsb) + k) = matstatus_true(check_elems(i));
            newimmstatus((j * ipmbinsb) + k) = immstatus_true(check_elems(i));
            newindataset((j * ipmbinsb) + k) = indataset_true(check_elems(i));
            newbinhalfwidth((j * ipmbinsb) + k) = standard_midpoint_a;
            newbinhalfwidthb((j * ipmbinsb) + k) = standard_midpoint_b;
            newbinhalfwidthc((j * ipmbinsb) + k) = binhalfwidthc_true(check_elems(i));
            newgroup((j * ipmbinsb) + k) = group_true(check_elems(i));
            newcomments((j * ipmbinsb) + k) = comments_true(check_elems(i));
            
            std::string sizenums1 = std::to_string(newsizes((j * ipmbinsb) + k));
            std::string sizenums2 = std::to_string(newsizesb((j * ipmbinsb) + k));
            newstagenames((j * ipmbinsb) + k) = "sza_" + sizenums1.substr(0, 5);
            newstagenames((j * ipmbinsb) + k) += "_szb_" + sizenums2.substr(0, 5);
            newstagenames((j * ipmbinsb) + k) += "_";
            newstagenames((j * ipmbinsb) + k) += std::to_string(newgroup((j * ipmbinsb) + k));
          }
        }
        
        sizes = concat_dbl(sizes, newsizes);
        sizesb_true = concat_dbl(sizesb_true, newsizesb);
        sizesc_true = concat_dbl(sizesc_true, newsizesc);
        minage_true = concat_dbl(minage_true, newminage);
        maxage_true = concat_dbl(maxage_true, newmaxage);
        repstatus_true = concat_int(repstatus_true, newrepstatus);
        obsstatus_true = concat_int(obsstatus_true, newobsstatus);
        propstatus_true = concat_int(propstatus_true, newpropstatus);
        matstatus_true = concat_int(matstatus_true, newmatstatus);
        immstatus_true = concat_int(immstatus_true, newimmstatus);
        indataset_true = concat_int(indataset_true, newindataset);
        binhalfwidth_true = concat_dbl(binhalfwidth_true, newbinhalfwidth);
        binhalfwidthb_true = concat_dbl(binhalfwidthb_true, newbinhalfwidthb);
        binhalfwidthc_true = concat_dbl(binhalfwidthc_true, newbinhalfwidthc);
        stagenames_true = concat_str(stagenames_true, newstagenames);
        group_true = concat_int(group_true, newgroup);
        comments_true = concat_str(comments_true, newcomments);
      }
    }
  }
  
  // Automated size classification with sizea & sizec together
  if (ipm_ac > 0) {
    if (IntegerVector::is_na(ipmbins)) {
      throw Rcpp::exception("The number of bins for automated size classification of the primary size variable has not been set.", false);
    }
    
    if (IntegerVector::is_na(ipmbinsc)) {
      throw Rcpp::exception("The number of bins for automated size classification of the primary size variable has not been set.", false);
    }
    
    if (ipm_ac % 2 != 0) {
      throw Rcpp::exception("The ipm designation must specify both the start size and the end size, requiring an even number of calls. Calls for automated size classification must be matched and not overlap.", false);
    }
    
    arma::uvec check_elems = find(ipm_calls_ac); // This vector points out which rows have ipm designations
    
    arma::ivec called_pair_check = pair_check.elem(check_elems); // 
    int called_pair_check_length = called_pair_check.n_elem;
    
    for (int i = 0; i < called_pair_check_length; i++) {
      int trial_1 = called_pair_check(i);
      int match_count {0};
      int main_index_1 {i};
      int go_ahead {0};
      
      for (int j = 0; j < called_pair_check_length; j++) {
        if (trial_1 == called_pair_check(j) && i != j) {
          match_count++;
          main_index_1 = j;
          
          if (j > i) {
            go_ahead = 1;
          } else {
            go_ahead = 0;
          }
        }
      }
      
      if (go_ahead == 1) {
        if (match_count > 1) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. More than 2 stages with the same characteristics marked ipm cannot be handled.", false);
        } else if (match_count == 0) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. Single stages with unique characteristics marked ipm cannot be handled.", false);
        }
        
        // Now we work out the minimum and maximum sizes
        double minsize_a {0};
        double maxsize_a {0};
        double minsize_c {0};
        double maxsize_c {0};
        
        if (sizes(check_elems(i)) > sizes(check_elems(main_index_1))) {
          maxsize_a = sizes(check_elems(i));
          minsize_a = sizes(check_elems(main_index_1));
        } else if (sizes(check_elems(i)) < sizes(check_elems(main_index_1))) {
          minsize_a = sizes(check_elems(i));
          maxsize_a = sizes(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        if (sizesc_true(check_elems(i)) > sizesc_true(check_elems(main_index_1))) {
          maxsize_c = sizesc_true(check_elems(i));
          minsize_c = sizesc_true(check_elems(main_index_1));
        } else if (sizesc_true(check_elems(i)) < sizesc_true(check_elems(main_index_1))) {
          minsize_c = sizesc_true(check_elems(i));
          maxsize_c = sizesc_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        
        // Keep track of the entries to delete from the original vectors
        if (no_entries_to_delete == 0) {
          entries_to_delete.resize(2);
          
          entries_to_delete(0) = check_elems(i);
          entries_to_delete(1) = check_elems(main_index_1);
        } else {
          int entries_del_old_size = entries_to_delete.size();
          entries_to_delete.resize(entries_del_old_size + 2);
          
          entries_to_delete(entries_del_old_size) = check_elems(i);
          entries_to_delete(entries_del_old_size + 1) = check_elems(main_index_1);
        }
        
        no_entries_to_delete++;
        no_entries_to_delete++;
        
        // Now we will create the new stages
        double full_range_a = maxsize_a - minsize_a;
        double standard_increment_a = full_range_a / static_cast<double>(ipmbins);
        double standard_midpoint_a = standard_increment_a / 2;
        
        double full_range_c = maxsize_c - minsize_c;
        double standard_increment_c = full_range_c / static_cast<double>(ipmbinsc);
        double standard_midpoint_c = standard_increment_c / 2;
        
        // New vectors to append
        NumericVector newsizes ((ipmbins * ipmbinsc), 0.0);
        StringVector newstagenames ((ipmbins * ipmbinsc), "");
        NumericVector newsizesb ((ipmbins * ipmbinsc), 0.0);
        NumericVector newsizesc ((ipmbins * ipmbinsc), 0.0);
        NumericVector newminage ((ipmbins * ipmbinsc), 0.0);
        NumericVector newmaxage ((ipmbins * ipmbinsc), 0.0);
        IntegerVector newrepstatus ((ipmbins * ipmbinsc), 0);
        IntegerVector newobsstatus ((ipmbins * ipmbinsc), 0);
        IntegerVector newpropstatus ((ipmbins * ipmbinsc), 0);
        IntegerVector newmatstatus ((ipmbins * ipmbinsc), 0);
        IntegerVector newimmstatus ((ipmbins * ipmbinsc), 0);
        IntegerVector newindataset ((ipmbins * ipmbinsc), 0);
        NumericVector newbinhalfwidth ((ipmbins * ipmbinsc), 0.0);
        NumericVector newbinhalfwidthb ((ipmbins * ipmbinsc), 0.0);
        NumericVector newbinhalfwidthc ((ipmbins * ipmbinsc), 0.0);
        IntegerVector newgroup ((ipmbins * ipmbinsc), 0);
        StringVector newcomments ((ipmbins * ipmbinsc), "");
        
        for (int j = 0; j < ipmbins; j++) {
          for (int k = 0; k < ipmbinsc; k++) {
            newsizes((j * ipmbinsc) + k) = minsize_a + (j*standard_increment_a) + standard_midpoint_a;
            newsizesc((j * ipmbinsc) + k) = minsize_c + (k*standard_increment_c) + standard_midpoint_c;;
            newsizesb((j * ipmbinsc) + k) = sizesb_true(check_elems(i));
            newminage((j * ipmbinsc) + k) = minage_true(check_elems(i));
            newmaxage((j * ipmbinsc) + k) = maxage_true(check_elems(i));
            newrepstatus((j * ipmbinsc) + k) = repstatus_true(check_elems(i));
            newobsstatus((j * ipmbinsc) + k) = obsstatus_true(check_elems(i));
            newpropstatus((j * ipmbinsc) + k) = propstatus_true(check_elems(i));
            newmatstatus((j * ipmbinsc) + k) = matstatus_true(check_elems(i));
            newimmstatus((j * ipmbinsc) + k) = immstatus_true(check_elems(i));
            newindataset((j * ipmbinsc) + k) = indataset_true(check_elems(i));
            newbinhalfwidth((j * ipmbinsc) + k) = standard_midpoint_a;
            newbinhalfwidthb((j * ipmbinsc) + k) = binhalfwidthb_true(check_elems(i));
            newbinhalfwidthc((j * ipmbinsc) + k) = standard_midpoint_c;
            newgroup((j * ipmbinsc) + k) = group_true(check_elems(i));
            newcomments((j * ipmbinsc) + k) = comments_true(check_elems(i));
            
            std::string sizenums1 = std::to_string(newsizes((j * ipmbinsc) + k));
            std::string sizenums2 = std::to_string(newsizesc((j * ipmbinsc) + k));
            newstagenames((j * ipmbinsc) + k) = "sza_" + sizenums1.substr(0, 5);
            newstagenames((j * ipmbinsc) + k) += "_szc_" + sizenums2.substr(0, 5);
            newstagenames((j * ipmbinsc) + k) += "_";
            newstagenames((j * ipmbinsc) + k) += std::to_string(newgroup((j * ipmbinsc)));
          }
        }
        
        sizes = concat_dbl(sizes, newsizes);
        sizesb_true = concat_dbl(sizesb_true, newsizesb);
        sizesc_true = concat_dbl(sizesc_true, newsizesc);
        minage_true = concat_dbl(minage_true, newminage);
        maxage_true = concat_dbl(maxage_true, newmaxage);
        repstatus_true = concat_int(repstatus_true, newrepstatus);
        obsstatus_true = concat_int(obsstatus_true, newobsstatus);
        propstatus_true = concat_int(propstatus_true, newpropstatus);
        matstatus_true = concat_int(matstatus_true, newmatstatus);
        immstatus_true = concat_int(immstatus_true, newimmstatus);
        indataset_true = concat_int(indataset_true, newindataset);
        binhalfwidth_true = concat_dbl(binhalfwidth_true, newbinhalfwidth);
        binhalfwidthb_true = concat_dbl(binhalfwidthb_true, newbinhalfwidthb);
        binhalfwidthc_true = concat_dbl(binhalfwidthc_true, newbinhalfwidthc);
        stagenames_true = concat_str(stagenames_true, newstagenames);
        group_true = concat_int(group_true, newgroup);
        comments_true = concat_str(comments_true, newcomments);
      }
    }
  }
  
  // Automated size classification with sizeb & sizec together
  if (ipm_bc > 0) {
    if (IntegerVector::is_na(ipmbinsb)) {
      throw Rcpp::exception("The number of bins for automated size classification of the primary size variable has not been set.", false);
    }
    
    if (IntegerVector::is_na(ipmbinsc)) {
      throw Rcpp::exception("The number of bins for automated size classification of the primary size variable has not been set.", false);
    }
    
    if (ipm_bc % 2 != 0) {
      throw Rcpp::exception("The ipm designation must specify both the start size and the end size, requiring an even number of calls. Calls for automated size classification must be matched and not overlap.", false);
    }
    
    arma::uvec check_elems = find(ipm_calls_bc); // This vector points out which rows have ipm designations
    
    arma::ivec called_pair_check = pair_check.elem(check_elems); // 
    int called_pair_check_length = called_pair_check.n_elem;
    
    for (int i = 0; i < called_pair_check_length; i++) {
      int trial_1 = called_pair_check(i);
      int match_count {0};
      int main_index_1 {i};
      int go_ahead {0};
      
      for (int j = 0; j < called_pair_check_length; j++) {
        if (trial_1 == called_pair_check(j) && i != j) {
          match_count++;
          main_index_1 = j;
          
          if (j > i) {
            go_ahead = 1;
          } else {
            go_ahead = 0;
          }
        }
      }
      
      if (go_ahead == 1) {
        if (match_count > 1) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. More than 2 stages with the same characteristics marked ipm cannot be handled.", false);
        } else if (match_count == 0) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. Single stages with unique characteristics marked ipm cannot be handled.", false);
        }
        
        // Now we work out the minimum and maximum sizes
        double minsize_b {0};
        double maxsize_b {0};
        double minsize_c {0};
        double maxsize_c {0};
        
        if (sizesb_true(check_elems(i)) > sizesb_true(check_elems(main_index_1))) {
          maxsize_b = sizesb_true(check_elems(i));
          minsize_b = sizesb_true(check_elems(main_index_1));
        } else if (sizesb_true(check_elems(i)) < sizesb_true(check_elems(main_index_1))) {
          minsize_b = sizesb_true(check_elems(i));
          maxsize_b = sizesb_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        if (sizesc_true(check_elems(i)) > sizesc_true(check_elems(main_index_1))) {
          maxsize_c = sizesc_true(check_elems(i));
          minsize_c = sizesc_true(check_elems(main_index_1));
        } else if (sizesc_true(check_elems(i)) < sizesc_true(check_elems(main_index_1))) {
          minsize_c = sizesc_true(check_elems(i));
          maxsize_c = sizesc_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        
        // Keep track of the entries to delete from the original vectors
        if (no_entries_to_delete == 0) {
          entries_to_delete.resize(2);
          
          entries_to_delete(0) = check_elems(i);
          entries_to_delete(1) = check_elems(main_index_1);
        } else {
          int entries_del_old_size = entries_to_delete.size();
          entries_to_delete.resize(entries_del_old_size + 2);
          
          entries_to_delete(entries_del_old_size) = check_elems(i);
          entries_to_delete(entries_del_old_size + 1) = check_elems(main_index_1);
        }
        
        no_entries_to_delete++;
        no_entries_to_delete++;
        
        // Now we will create the new stages
        double full_range_b = maxsize_b - minsize_b;
        double standard_increment_b = full_range_b / static_cast<double>(ipmbinsb);
        double standard_midpoint_b = standard_increment_b / 2;
        
        double full_range_c = maxsize_c - minsize_c;
        double standard_increment_c = full_range_c / static_cast<double>(ipmbinsc);
        double standard_midpoint_c = standard_increment_c / 2;
        
        // New vectors to append
        NumericVector newsizes ((ipmbinsb * ipmbinsc), 0.0);
        StringVector newstagenames ((ipmbinsb * ipmbinsc), "");
        NumericVector newsizesb ((ipmbinsb * ipmbinsc), 0.0);
        NumericVector newsizesc ((ipmbinsb * ipmbinsc), 0.0);
        NumericVector newminage ((ipmbinsb * ipmbinsc), 0.0);
        NumericVector newmaxage ((ipmbinsb * ipmbinsc), 0.0);
        IntegerVector newrepstatus ((ipmbinsb * ipmbinsc), 0);
        IntegerVector newobsstatus ((ipmbinsb * ipmbinsc), 0);
        IntegerVector newpropstatus ((ipmbinsb * ipmbinsc), 0);
        IntegerVector newmatstatus ((ipmbinsb * ipmbinsc), 0);
        IntegerVector newimmstatus ((ipmbinsb * ipmbinsc), 0);
        IntegerVector newindataset ((ipmbinsb * ipmbinsc), 0);
        NumericVector newbinhalfwidth ((ipmbinsb * ipmbinsc), 0.0);
        NumericVector newbinhalfwidthb ((ipmbinsb * ipmbinsc), 0.0);
        NumericVector newbinhalfwidthc ((ipmbinsb * ipmbinsc), 0.0);
        IntegerVector newgroup ((ipmbinsb * ipmbinsc), 0);
        StringVector newcomments ((ipmbinsb * ipmbinsc), "");
        
        for (int j = 0; j < ipmbinsb; j++) {
          for (int k = 0; k < ipmbinsc; k++) {
            newsizesb((j * ipmbinsc) + k) = minsize_b + (j*standard_increment_b) + standard_midpoint_b;
            newsizesc((j * ipmbinsc) + k) = minsize_c + (k*standard_increment_c) + standard_midpoint_c;;
            newsizes((j * ipmbinsc) + k) = sizes(check_elems(i));
            newminage((j * ipmbinsc) + k) = minage_true(check_elems(i));
            newmaxage((j * ipmbinsc) + k) = maxage_true(check_elems(i));
            newrepstatus((j * ipmbinsc) + k) = repstatus_true(check_elems(i));
            newobsstatus((j * ipmbinsc) + k) = obsstatus_true(check_elems(i));
            newpropstatus((j * ipmbinsc) + k) = propstatus_true(check_elems(i));
            newmatstatus((j * ipmbinsc) + k) = matstatus_true(check_elems(i));
            newimmstatus((j * ipmbinsc) + k) = immstatus_true(check_elems(i));
            newindataset((j * ipmbinsc) + k) = indataset_true(check_elems(i));
            newbinhalfwidth((j * ipmbinsc) + k) = binhalfwidth_true(check_elems(i));
            newbinhalfwidthb((j * ipmbinsc) + k) = standard_midpoint_b;
            newbinhalfwidthc((j * ipmbinsc) + k) = standard_midpoint_c;
            newgroup((j * ipmbinsc) + k) = group_true(check_elems(i));
            newcomments((j * ipmbinsc) + k) = comments_true(check_elems(i));
            
            std::string sizenums1 = std::to_string(newsizesb((j * ipmbinsc) + k));
            std::string sizenums2 = std::to_string(newsizesc((j * ipmbinsc) + k));
            newstagenames((j * ipmbinsc) + k) = "szb_" + sizenums1.substr(0, 5);
            newstagenames((j * ipmbinsc) + k) += "_szc_" + sizenums2.substr(0, 5);
            newstagenames((j * ipmbinsc) + k) += "_";
            newstagenames((j * ipmbinsc) + k) += std::to_string(newgroup((j * ipmbinsc) + k));
          }
        }
        
        sizes = concat_dbl(sizes, newsizes);
        sizesb_true = concat_dbl(sizesb_true, newsizesb);
        sizesc_true = concat_dbl(sizesc_true, newsizesc);
        minage_true = concat_dbl(minage_true, newminage);
        maxage_true = concat_dbl(maxage_true, newmaxage);
        repstatus_true = concat_int(repstatus_true, newrepstatus);
        obsstatus_true = concat_int(obsstatus_true, newobsstatus);
        propstatus_true = concat_int(propstatus_true, newpropstatus);
        matstatus_true = concat_int(matstatus_true, newmatstatus);
        immstatus_true = concat_int(immstatus_true, newimmstatus);
        indataset_true = concat_int(indataset_true, newindataset);
        binhalfwidth_true = concat_dbl(binhalfwidth_true, newbinhalfwidth);
        binhalfwidthb_true = concat_dbl(binhalfwidthb_true, newbinhalfwidthb);
        binhalfwidthc_true = concat_dbl(binhalfwidthc_true, newbinhalfwidthc);
        stagenames_true = concat_str(stagenames_true, newstagenames);
        group_true = concat_int(group_true, newgroup);
        comments_true = concat_str(comments_true, newcomments);
      }
    }
  }
  
  // Automated size classification with sizea, sizeb, & sizec together
  if (ipm_abc > 0) {
    if (IntegerVector::is_na(ipmbins)) {
      throw Rcpp::exception("The number of bins for automated size classification of the primary size variable has not been set.", false);
    }
    
    if (IntegerVector::is_na(ipmbinsb)) {
      throw Rcpp::exception("The number of bins for automated size classification of the primary size variable has not been set.", false);
    }
    
    if (IntegerVector::is_na(ipmbinsc)) {
      throw Rcpp::exception("The number of bins for automated size classification of the primary size variable has not been set.", false);
    }
    
    if (ipm_abc % 2 != 0) {
      throw Rcpp::exception("The ipm designation must specify both the start size and the end size, requiring an even number of calls. Calls for automated size classification must be matched and not overlap.", false);
    }
    
    arma::uvec check_elems = find(ipm_calls_abc); // This vector points out which rows have ipm designations
    
    arma::ivec called_pair_check = pair_check.elem(check_elems); // 
    int called_pair_check_length = called_pair_check.n_elem;
    
    for (int i = 0; i < called_pair_check_length; i++) {
      int trial_1 = called_pair_check(i);
      int match_count {0};
      int main_index_1 {i};
      int go_ahead {0};
      
      for (int j = 0; j < called_pair_check_length; j++) {
        if (trial_1 == called_pair_check(j) && i != j) {
          match_count++;
          main_index_1 = j;
          
          if (j > i) {
            go_ahead = 1;
          } else {
            go_ahead = 0;
          }
        }
      }
      
      if (go_ahead == 1) {
        if (match_count > 1) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. More than 2 stages with the same characteristics marked ipm cannot be handled.", false);
        } else if (match_count == 0) {
          throw Rcpp::exception("Stages marked ipm must have equal characteristics in pairs only, corresponding to size minimum and maximum. Single stages with unique characteristics marked ipm cannot be handled.", false);
        }
        
        // Now we work out the minimum and maximum sizes
        double minsize_a {0};
        double maxsize_a {0};
        double minsize_b {0};
        double maxsize_b {0};
        double minsize_c {0};
        double maxsize_c {0};
        
        if (sizes(check_elems(i)) > sizes(check_elems(main_index_1))) {
          maxsize_a = sizes(check_elems(i));
          minsize_a = sizes(check_elems(main_index_1));
        } else if (sizes(check_elems(i)) < sizes(check_elems(main_index_1))) {
          minsize_a = sizes(check_elems(i));
          maxsize_a = sizes(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        if (sizesb_true(check_elems(i)) > sizesb_true(check_elems(main_index_1))) {
          maxsize_b = sizesb_true(check_elems(i));
          minsize_b = sizesb_true(check_elems(main_index_1));
        } else if (sizesb_true(check_elems(i)) < sizesb_true(check_elems(main_index_1))) {
          minsize_b = sizesb_true(check_elems(i));
          maxsize_b = sizesb_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        if (sizesc_true(check_elems(i)) > sizesc_true(check_elems(main_index_1))) {
          maxsize_c = sizesc_true(check_elems(i));
          minsize_c = sizesc_true(check_elems(main_index_1));
        } else if (sizesc_true(check_elems(i)) < sizesc_true(check_elems(main_index_1))) {
          minsize_c = sizesc_true(check_elems(i));
          maxsize_c = sizesc_true(check_elems(main_index_1));
        } else {
          throw Rcpp::exception("Size values used in IPM classification must differ.", false);
        }
        
        // Keep track of the entries to delete from the original vectors
        if (no_entries_to_delete == 0) {
          entries_to_delete.resize(2);
          
          entries_to_delete(0) = check_elems(i);
          entries_to_delete(1) = check_elems(main_index_1);
        } else {
          int entries_del_old_size = entries_to_delete.size();
          entries_to_delete.resize(entries_del_old_size + 2);
          
          entries_to_delete(entries_del_old_size) = check_elems(i);
          entries_to_delete(entries_del_old_size + 1) = check_elems(main_index_1);
        }
        
        no_entries_to_delete++;
        no_entries_to_delete++;
        
        // Now we will create the new stages
        double full_range_a = maxsize_a - minsize_a;
        double standard_increment_a = full_range_a / static_cast<double>(ipmbins);
        double standard_midpoint_a = standard_increment_a / 2;
        
        double full_range_b = maxsize_b - minsize_b;
        double standard_increment_b = full_range_b / static_cast<double>(ipmbinsb);
        double standard_midpoint_b = standard_increment_b / 2;
        
        double full_range_c = maxsize_c - minsize_c;
        double standard_increment_c = full_range_c / static_cast<double>(ipmbinsc);
        double standard_midpoint_c = standard_increment_c / 2;
        
        // New vectors to append
        NumericVector newsizes ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        StringVector newstagenames ((ipmbins * ipmbinsb * ipmbinsc), "");
        NumericVector newsizesb ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        NumericVector newsizesc ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        NumericVector newminage ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        NumericVector newmaxage ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        IntegerVector newrepstatus ((ipmbins * ipmbinsb * ipmbinsc), 0);
        IntegerVector newobsstatus ((ipmbins * ipmbinsb * ipmbinsc), 0);
        IntegerVector newpropstatus ((ipmbins * ipmbinsb * ipmbinsc), 0);
        IntegerVector newmatstatus ((ipmbins * ipmbinsb * ipmbinsc), 0);
        IntegerVector newimmstatus ((ipmbins * ipmbinsb * ipmbinsc), 0);
        IntegerVector newindataset ((ipmbins * ipmbinsb * ipmbinsc), 0);
        NumericVector newbinhalfwidth ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        NumericVector newbinhalfwidthb ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        NumericVector newbinhalfwidthc ((ipmbins * ipmbinsb * ipmbinsc), 0.0);
        IntegerVector newgroup ((ipmbins * ipmbinsb * ipmbinsc), 0);
        StringVector newcomments ((ipmbins * ipmbinsb * ipmbinsc), "");
        
        for (int j = 0; j < ipmbins; j++) {
          for (int k = 0; k < ipmbinsb; k++) {
            for (int l = 0; l < ipmbinsc; l++) {
              newsizes((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = minsize_a + 
                (j*standard_increment_a) + standard_midpoint_a;
              newsizesb((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = minsize_b + 
                (k*standard_increment_b) + standard_midpoint_b;
              newsizesc((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = minsize_c + 
                (l*standard_increment_c) + standard_midpoint_c;
              newminage((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = minage_true(check_elems(i));
              newmaxage((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = maxage_true(check_elems(i));
              newrepstatus((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = repstatus_true(check_elems(i));
              newobsstatus((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = obsstatus_true(check_elems(i));
              newpropstatus((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = propstatus_true(check_elems(i));
              newmatstatus((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = matstatus_true(check_elems(i));
              newimmstatus((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = immstatus_true(check_elems(i));
              newindataset((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = indataset_true(check_elems(i));
              newbinhalfwidth((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = standard_midpoint_a;
              newbinhalfwidthb((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = standard_midpoint_b;
              newbinhalfwidthc((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = standard_midpoint_c;
              newgroup((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = group_true(check_elems(i));
              newcomments((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = comments_true(check_elems(i));
              
              std::string sizenums1 = std::to_string(newsizes((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l));
              std::string sizenums2 = std::to_string(newsizesb((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l));
              std::string sizenums3 = std::to_string(newsizesc((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l));
              newstagenames((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) = "sza_" + sizenums1.substr(0, 4);
              newstagenames((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) += "_szb_" + sizenums2.substr(0, 4);
              newstagenames((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) += "_szc_" + sizenums3.substr(0, 4);
              newstagenames((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) += "_";
              newstagenames((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l) += 
                std::to_string(newgroup((j * ipmbinsc * ipmbinsb) + (k * ipmbinsc) + l));
            }
          }
        }
        
        sizes = concat_dbl(sizes, newsizes);
        sizesb_true = concat_dbl(sizesb_true, newsizesb);
        sizesc_true = concat_dbl(sizesc_true, newsizesc);
        minage_true = concat_dbl(minage_true, newminage);
        maxage_true = concat_dbl(maxage_true, newmaxage);
        repstatus_true = concat_int(repstatus_true, newrepstatus);
        obsstatus_true = concat_int(obsstatus_true, newobsstatus);
        propstatus_true = concat_int(propstatus_true, newpropstatus);
        matstatus_true = concat_int(matstatus_true, newmatstatus);
        immstatus_true = concat_int(immstatus_true, newimmstatus);
        indataset_true = concat_int(indataset_true, newindataset);
        binhalfwidth_true = concat_dbl(binhalfwidth_true, newbinhalfwidth);
        binhalfwidthb_true = concat_dbl(binhalfwidthb_true, newbinhalfwidthb);
        binhalfwidthc_true = concat_dbl(binhalfwidthc_true, newbinhalfwidthc);
        stagenames_true = concat_str(stagenames_true, newstagenames);
        group_true = concat_int(group_true, newgroup);
        comments_true = concat_str(comments_true, newcomments);
      }
    }
  }
  
  // Delete unneeded elements - run for loop backwards to delete elements properly
  arma::uvec unique_delete = unique(entries_to_delete);
  int no_unique_delete = unique_delete.n_elem;
  
  if (no_entries_to_delete > 0) {
    for (int i = 0; i < no_unique_delete; i++) {
      sizes.erase(unique_delete(no_unique_delete - (i + 1)));
      sizesb_true.erase(unique_delete(no_unique_delete - (i + 1)));
      sizesc_true.erase(unique_delete(no_unique_delete - (i + 1)));
      stagenames_true.erase(unique_delete(no_unique_delete - (i + 1)));
      
      minage_true.erase(unique_delete(no_unique_delete - (i + 1)));
      maxage_true.erase(unique_delete(no_unique_delete - (i + 1)));
      
      repstatus_true.erase(unique_delete(no_unique_delete - (i + 1)));
      obsstatus_true.erase(unique_delete(no_unique_delete - (i + 1)));
      propstatus_true.erase(unique_delete(no_unique_delete - (i + 1)));
      matstatus_true.erase(unique_delete(no_unique_delete - (i + 1)));
      immstatus_true.erase(unique_delete(no_unique_delete - (i + 1)));
      indataset_true.erase(unique_delete(no_unique_delete - (i + 1)));
      
      binhalfwidth_true.erase(unique_delete(no_unique_delete - (i + 1)));
      binhalfwidthb_true.erase(unique_delete(no_unique_delete - (i + 1)));
      binhalfwidthc_true.erase(unique_delete(no_unique_delete - (i + 1)));
      
      group_true.erase(unique_delete(no_unique_delete - (i + 1)));
      
      comments_true.erase(unique_delete(no_unique_delete - (i + 1)));
    }
  }
  
  // Post-IPM calculations
  matsize = sizes.size(); // Redefined length
  int elems_to_grow = matsize - sizebin_min.length();
  NumericVector zeros_to_append (elems_to_grow, NA_REAL);
  
  sizebin_min = concat_dbl(sizebin_min, zeros_to_append);
  sizebin_max = concat_dbl(sizebin_max, zeros_to_append);
  sizebin_center = concat_dbl(sizebin_center, zeros_to_append);
  sizebin_width = concat_dbl(sizebin_width, zeros_to_append);
  
  sizebinb_min = concat_dbl(sizebinb_min, zeros_to_append);
  sizebinb_max = concat_dbl(sizebinb_max, zeros_to_append);
  sizebinb_center = concat_dbl(sizebinb_center, zeros_to_append);
  sizebinb_width = concat_dbl(sizebinb_width, zeros_to_append);
  
  sizebinc_min = concat_dbl(sizebinc_min, zeros_to_append);
  sizebinc_max = concat_dbl(sizebinc_max, zeros_to_append);
  sizebinc_center = concat_dbl(sizebinc_center, zeros_to_append);
  sizebinc_width = concat_dbl(sizebinc_width, zeros_to_append);
  
  for (int i = 0; i < matsize; i++) {
    sizebin_min(i) = sizes(i) - binhalfwidth_true(i);
    sizebin_max(i) = sizes(i) + binhalfwidth_true(i);
    sizebin_center(i) = sizebin_min(i) + ((sizebin_max(i) - sizebin_min(i))/ 2);
    sizebin_width(i) = sizebin_max(i) - sizebin_min(i);
    
    if (used_sizes > 1) {
      sizebinb_min(i) = sizesb_true(i) - binhalfwidthb_true(i);
      sizebinb_max(i) = sizesb_true(i) + binhalfwidthb_true(i);
      sizebinb_center(i) = sizebinb_min(i) + ((sizebinb_max(i) - sizebinb_min(i))/ 2);
      sizebinb_width(i) = sizebinb_max(i) - sizebinb_min(i);
    }
    if (used_sizes > 2) {
      sizebinc_min(i) = sizesc_true(i) - binhalfwidthc_true(i);
      sizebinc_max(i) = sizesc_true(i) + binhalfwidthc_true(i);
      sizebinc_center(i) = sizebinc_min(i) + ((sizebinc_max(i) - sizebinc_min(i))/ 2);
      sizebinc_width(i) = sizebinc_max(i) - sizebinc_min(i);
    }
  }
  sizebin_min = round(sizebin_min, roundsize);
  sizebin_max = round(sizebin_max, roundsize);
  sizebin_center = round(sizebin_center, roundsize);
  sizebin_width = round(sizebin_width, roundsize);
  
  if (used_sizes > 1) {
    sizebinb_min = round(sizebinb_min, roundsizeb);
    sizebinb_max = round(sizebinb_max, roundsizeb);
    sizebinb_center = round(sizebinb_center, roundsizeb);
    sizebinb_width = round(sizebinb_width, roundsizeb);
  }
  if (used_sizes > 2) {
    sizebinc_min = round(sizebinc_min, roundsizec);
    sizebinc_max = round(sizebinc_max, roundsizec);
    sizebinc_center = round(sizebinc_center, roundsizec);
    sizebinc_width = round(sizebinc_width, roundsizec);
  }
  
  output_longlist(0) = stagenames_true;
  output_longlist(1) = sizes;
  output_longlist(2) = sizesb_true;
  output_longlist(3) = sizesc_true;
  
  output_longlist(4) = minage_true;
  output_longlist(5) = maxage_true;
  output_longlist(6) = repstatus_true;
  output_longlist(7) = obsstatus_true;
  output_longlist(8) = propstatus_true;
  output_longlist(9) = immstatus_true;
  output_longlist(10) = matstatus_true;
  
  output_longlist(11) = indataset_true;
  
  output_longlist(12) = binhalfwidth_true;
  output_longlist(13) = sizebin_min;
  output_longlist(14) = sizebin_max;
  output_longlist(15) = sizebin_center;
  output_longlist(16) = sizebin_width;
  
  output_longlist(17) = binhalfwidthb_true;
  output_longlist(18) = sizebinb_min;
  output_longlist(19) = sizebinb_max;
  output_longlist(20) = sizebinb_center;
  output_longlist(21) = sizebinb_width;
  
  output_longlist(22) = binhalfwidthc_true;
  output_longlist(23) = sizebinc_min;
  output_longlist(24) = sizebinc_max;
  output_longlist(25) = sizebinc_center;
  output_longlist(26) = sizebinc_width;
  
  output_longlist(27) = group_true;
  output_longlist(28) = comments_true;
  
  output_longlist.attr("names") = varnames;
  output_longlist.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, matsize);
  StringVector needed_classes {"data.frame", "stageframe"};
  output_longlist.attr("class") = needed_classes; // data.frame
  
  return output_longlist;
}

//' Standardize Stageframe For MPM Analysis
//' 
//' Function \code{.sf_reassess()} takes a stageframe as input, and uses
//' information supplied there and through the supplement, reproduction and
//' overwrite tables to rearrange this into a format usable by the matrix
//' creation functions, \code{\link{flefko3}()}, \code{\link{flefko2}()},
//' \code{\link{aflefko2}()}, \code{\link{rlefko3}()}, and
//' \code{\link{rlefko2}()}.
//' 
//' @param stageframe The original stageframe.
//' @param supplement The original supplemental data input
//' (class \code{lefkoSD}). Can also equal NA.
//' @param overwrite An overwrite table.
//' @param repmatrix The original reproduction matrix. Can also equal NA or 0.
//' @param agemat A logical value indicating whether MPM is age-by-stage.
//' @param historical A logical value indicating whether MPM is historical.
//' @param format An integer indicating whether matrices will be in Ehrlen format
//' (if set to 1), or deVries format (if set to 2). Setting to deVries format
//' adds one extra stage to account for the prior status of newborns.
//' 
//' @return This function returns a list with a modified stageframe usable in MPM
//' construction, an associated reproduction matrix, and a general supplement
//' table that takes over the input supplement and overwrite tables. Note that
//' if a supplement is provided and a repmatrix is not, or if repmatrix is set
//' to 0, then it will be assumed that a repmatrix should not be used.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.sf_reassess)]]
Rcpp::List sf_reassess(DataFrame stageframe, Nullable<DataFrame> supplement,
  Nullable<DataFrame> overwrite, Nullable<NumericMatrix> repmatrix,
  bool agemat = false, bool historical = false, int format = 1) {
  
  bool skiprepmat = false;
  bool supp_provided = false;
  bool repm_provided = false;
  bool over_provided = false;
  
  Rcpp::DataFrame supplement_true;
  arma::mat repmatrix_true;
  int newsupp_rows {0};
  
  StringVector stagevec = stageframe["stage"];
  NumericVector origsizevec = stageframe["size"];
  NumericVector origsizebvec = stageframe["size_b"];
  NumericVector origsizecvec = stageframe["size_c"];
  NumericVector minagevec = stageframe["min_age"];
  NumericVector maxagevec = stageframe["max_age"];
  arma::uvec repvec = stageframe["repstatus"];
  arma::uvec obsvec = stageframe["obsstatus"];
  arma::uvec propvec = stageframe["propstatus"];
  arma::uvec immvec = stageframe["immstatus"];
  arma::uvec matvec = stageframe["matstatus"];
  arma::uvec indvec = stageframe["indataset"];
  NumericVector binvec = stageframe["binhalfwidth_raw"];
  NumericVector binbvec = stageframe["binhalfwidthb_raw"];
  NumericVector bincvec = stageframe["binhalfwidthc_raw"];
  NumericVector sizeminvec = stageframe["sizebin_min"];
  NumericVector sizemaxvec = stageframe["sizebin_max"];
  NumericVector sizectrvec = stageframe["sizebin_center"];
  NumericVector sizewidthvec = stageframe["sizebin_width"];
  NumericVector sizebminvec = stageframe["sizebinb_min"];
  NumericVector sizebmaxvec = stageframe["sizebinb_max"];
  NumericVector sizebctrvec = stageframe["sizebinb_center"];
  NumericVector sizebwidthvec = stageframe["sizebinb_width"];
  NumericVector sizecminvec = stageframe["sizebinc_min"];
  NumericVector sizecmaxvec = stageframe["sizebinc_max"];
  NumericVector sizecctrvec = stageframe["sizebinc_center"];
  NumericVector sizecwidthvec = stageframe["sizebinc_width"];
  arma::ivec groupvec = stageframe["group"];
  StringVector comvec = stageframe["comments"];
  
  StringVector stage3_supp;
  StringVector stage2_supp;
  StringVector stage1_supp;
  StringVector eststage3_supp;
  StringVector eststage2_supp;
  StringVector eststage1_supp;
  NumericVector givenrate_supp;
  NumericVector multiplier_supp;
  IntegerVector convtype_supp;
  IntegerVector convtype_t12_supp;
  int supp_rows {0};
  
  int stageframe_length = repvec.n_elem;
  arma::uvec repentryvec(stageframe_length, fill::zeros);
  
  // Identify all groups in the stageframe
  arma::ivec all_groups = unique(groupvec);
  int no_groups = all_groups.n_elem;
  StringVector group_text(no_groups);
  
  for (int i = 0; i < no_groups; i++) {
    group_text(i) = "group";
    group_text(i) += std::to_string(all_groups(i));
  }
  
  if (supplement.isNotNull()) {
    supp_provided = true;
    
    Rcpp::DataFrame supplement_thru(supplement);
    supplement_true = supplement_thru;
    
    stage3_supp = supplement_true["stage3"];
    stage2_supp = supplement_true["stage2"];
    stage1_supp = supplement_true["stage1"];
    eststage3_supp = supplement_true["eststage3"];
    eststage2_supp = supplement_true["eststage2"];
    eststage1_supp = supplement_true["eststage1"];
    givenrate_supp = supplement_true["givenrate"];
    multiplier_supp = supplement_true["multiplier"]; // Not in overwrite()
    convtype_supp = supplement_true["convtype"];
    convtype_t12_supp = supplement_true["convtype_t12"];
    supp_rows = stage3_supp.length();
  }
  
  if (repmatrix.isNotNull()) {
    repm_provided = true;
    
    NumericMatrix repmatrix_thru(repmatrix);
    repmatrix_true = as<arma::mat>(repmatrix_thru);
    
    arma::rowvec rep_sums = sum(repmatrix_true, 0);
    arma::vec rep_rowsums = sum(repmatrix_true, 1);
    double repmat_sum = sum(rep_sums);
    
    if (rep_sums.n_elem != stageframe_length) {
      throw Rcpp::exception("Object repmatrix must have the rows and columns equal to the number of rows in the stageframe.", false);
    }
    
    if (repmat_sum > 0) {
      for (int i = 0; i < rep_sums.n_elem; i++) {
        if (rep_sums(i) > 0) {
          repvec(i) = 1;
        } else {
          repvec(i) = 0;
        }
        
        if (rep_rowsums(i) > 0) {
          repentryvec(i) = 1;
        } else {
          repentryvec(i) = 0;
        }
      }
    } else {
      repentryvec = propvec + immvec;
      int rev_checksum = sum(repentryvec);
      
      if (rev_checksum == 0) {
        repentryvec(0) = 1;
        Rf_warning("No information on reproductive entry stages provided. Assuming the first stage is the entry stage into the life cycle.");
      }
      
      arma::mat token_mat(stageframe_length, stageframe_length, fill::zeros);
        
      arma::uvec repentry_calls = find(repentryvec);
      arma::uvec rep_calls = find(repvec);
      
      if (repentry_calls.n_elem > 0 && rep_calls.n_elem > 0) {
        for (int i = 0; i < repentry_calls.n_elem; i++) {
          for (int j = 0; j < rep_calls.n_elem; j++) {
            token_mat(repentry_calls(i), rep_calls(j)) = 1;
          }
        }
      }
      repmatrix_true = token_mat;
    }
  }  else if (supp_provided) {
    arma::ivec cv_supp_arma = as<arma::ivec>(convtype_supp);
    arma::uvec mult_elems = find(cv_supp_arma == 3);
    
    int fec_mults = mult_elems.n_elem;
    // arma::uvec needed_reprods(fec_mults, fill::zeros);
    for (int i = 0; i < fec_mults; i++) {
      for (int j = 0; j < stageframe_length; j++) {
        if (stage3_supp(mult_elems(i)) == stagevec(j)) {
          repentryvec(j) = 1;
        }
      }
      
      if (stage3_supp(mult_elems(i)) == "prop") {
        arma::uvec propvec_1elems = find(propvec);
        
        int propvec_no1s = propvec_1elems.n_elem;
        for (int j = 0; j < propvec_no1s; j++) {
          repentryvec(propvec_1elems(j)) = 1;
        }
      }
      
      if (stage3_supp(mult_elems(i)) == "immat") {
        arma::uvec immvec_1elems = find(immvec);
        
        int immvec_no1s = immvec_1elems.n_elem;
        for (int j = 0; j < immvec_no1s; j++) {
          repentryvec(immvec_1elems(j)) = 1;
        }
      }
      
      for (int j = 0; j < no_groups; j++) {
        if (stage3_supp(mult_elems(i)) == group_text(j)) {
          arma::uvec group_elems = find(groupvec == all_groups(j));
          
          int group_noelems = group_elems.n_elem;
          for (int k = 0; k < group_noelems; k++) {
            repentryvec(group_elems(k)) = 1;
          }
        }
      }
    }
    
    if (agemat) {
      arma::ivec needed_reprods(stageframe_length, fill::zeros);
      arma::vec needed_mults(stageframe_length, fill::zeros);
      
      for (int i = 0; i < fec_mults; i++) {
        for (int j = 0; j < stageframe_length; j++) {
          if (stage2_supp(mult_elems(i)) == stagevec(j)) {
            needed_reprods(j) = 1;
            needed_mults(j) = multiplier_supp(mult_elems(i));
          }
          
          if (stage2_supp(mult_elems(i)) == "rep") {
            arma::uvec repvec_1elems = find(repvec);
            
            int repvec_no1s = repvec_1elems.n_elem;
            for (int k = 0; k < repvec_no1s; k++) {
              needed_reprods(repvec_1elems(k)) = 1;
              needed_mults(repvec_1elems(k)) = multiplier_supp(mult_elems(i));
            }
          }
          
          for (int k = 0; k < no_groups; k++) {
            if (stage2_supp(mult_elems(i)) == group_text(k)) {
              arma::uvec group_elems = find(groupvec == all_groups(k));
              
              int group_noelems = group_elems.n_elem;
              for (int l = 0; l < group_noelems; l++) {
                needed_reprods(group_elems(l)) = 1;
                needed_mults(group_elems(l)) = multiplier_supp(mult_elems(i));
              }
            }
          }
        }
      }
      
      arma::mat token_mat(stageframe_length, stageframe_length, fill::zeros);
      
      arma::uvec repentry_calls = find(repentryvec);
      arma::uvec neededrep_calls = find(needed_reprods);
      
      if (repentry_calls.n_elem > 0 && neededrep_calls.n_elem > 0) {
        for (int i = 0; i < repentry_calls.n_elem; i++) {
          for (int j = 0; j < neededrep_calls.n_elem; j++) {
            token_mat(repentry_calls(i), neededrep_calls(j)) = needed_mults(neededrep_calls(j));
          }
        }
      }
      repmatrix_true = token_mat;
      skiprepmat = false;
    
    } else {
      arma::mat token_mat(stageframe_length, stageframe_length, fill::zeros);
      repmatrix_true = token_mat;
      skiprepmat = true;
    }
  } else {
    for (int i = 0; i < stageframe_length; i++) {
      if (propvec(i) > 0) {
        repentryvec(i) = 1;
      } else if (immvec(i) > 0) {
        repentryvec(i) = 1;
      }
    }
    
    int rev_checksum = sum(repentryvec);
    if(rev_checksum == 0) {
      repentryvec(0) = 1;
      Rf_warning("No information on reproductive entry stages provided. Assuming the first stage is the entry stage into the life cycle.");
    }
    
    arma::mat token_mat(stageframe_length, stageframe_length, fill::zeros);
      
    arma::uvec repentry_calls = find(repentryvec);
    arma::uvec rep_calls = find(repvec);
    
    if (repentry_calls.n_elem > 0 && rep_calls.n_elem > 0) {
      for (int i = 0; i < repentry_calls.n_elem; i++) {
        for (int j = 0; j < rep_calls.n_elem; j++) {
          token_mat(repentry_calls(i), rep_calls(j)) = 1;
        }
      }
    }
    repmatrix_true = token_mat;
  }
  
  // Now we reorder the stageframe
  arma::uvec prop_stages = find(propvec);
  arma::uvec prop0_stages = find(propvec == 0);
  arma::uvec imm_stages = find(immvec);
  arma::uvec imm0_stages = find(immvec == 0);
  arma::uvec p0_im_stages = intersect(prop0_stages, imm_stages);
  
  arma::uvec mat_stages = find(matvec);
  arma::uvec mat_imm0_stages = intersect(mat_stages, imm0_stages);
  arma::uvec rep_stages = find(repvec);
  arma::uvec rep0_stages = find(repvec == 0);
  arma::uvec rep0_mat_imm0_stages = intersect(mat_imm0_stages, rep0_stages);
  arma::uvec mat_rep0_stages = intersect(mat_stages, rep0_stages);
  arma::uvec neworder(stageframe_length, fill::zeros);
  arma::uvec tracked_elems(stageframe_length, fill::ones);
  
  int no_prop_stages = prop_stages.n_elem;
  int no_p0_im_stages = p0_im_stages.n_elem;
  int no_r0_im_mt_stages = rep0_mat_imm0_stages.n_elem;
  
  int counter {0};
  for (int j = 0; j < no_prop_stages; j++) {
    neworder(counter) = prop_stages(j);
    tracked_elems(prop_stages(j)) = 0;
    counter++;
  }
  for (int j = 0; j < no_p0_im_stages; j++) {
    neworder(counter) = p0_im_stages(j);
    tracked_elems(p0_im_stages(j)) = 0;
    counter++;
  }
  for (int j = 0; j < no_r0_im_mt_stages; j++) {
    neworder(counter) = rep0_mat_imm0_stages(j);
    tracked_elems(rep0_mat_imm0_stages(j)) = 0;
    counter++;
  }
  
  arma::uvec remaining_elems = find(tracked_elems);
  int no_remaining_elems = remaining_elems.n_elem;
  for (int j = 0; j < no_remaining_elems; j++) {
    neworder(counter) = remaining_elems(j);
    counter++;
  }
  
  StringVector newstagevec(stageframe_length + format);
  NumericVector neworigsizevec(stageframe_length + format);
  NumericVector neworigsizebvec(stageframe_length + format);
  NumericVector neworigsizecvec(stageframe_length + format);
  NumericVector newminagevec(stageframe_length + format);
  NumericVector newmaxagevec(stageframe_length + format);
  arma::uvec newrepvec((stageframe_length + format), fill::zeros);
  arma::uvec newobsvec((stageframe_length + format), fill::zeros);
  arma::uvec newpropvec((stageframe_length + format), fill::zeros);
  arma::uvec newimmvec((stageframe_length + format), fill::zeros);
  arma::uvec newmatvec((stageframe_length + format), fill::zeros);
  arma::uvec newrepentryvec((stageframe_length + format), fill::zeros);
  arma::uvec newindvec((stageframe_length + format), fill::zeros);
  NumericVector newbinvec((stageframe_length + format));
  NumericVector newbinbvec((stageframe_length + format));
  NumericVector newbincvec((stageframe_length + format));
  NumericVector newsizeminvec((stageframe_length + format));
  NumericVector newsizemaxvec((stageframe_length + format));
  NumericVector newsizectrvec((stageframe_length + format));
  NumericVector newsizewidthvec((stageframe_length + format));
  NumericVector newsizebminvec((stageframe_length + format));
  NumericVector newsizebmaxvec((stageframe_length + format));
  NumericVector newsizebctrvec((stageframe_length + format));
  NumericVector newsizebwidthvec((stageframe_length + format));
  NumericVector newsizecminvec((stageframe_length + format));
  NumericVector newsizecmaxvec((stageframe_length + format));
  NumericVector newsizecctrvec((stageframe_length + format));
  NumericVector newsizecwidthvec((stageframe_length + format));
  arma::ivec newgroupvec((stageframe_length + format), fill::zeros);
  StringVector newcomvec((stageframe_length + format));
  arma::uvec newalive((stageframe_length + format), fill::zeros);
  
  arma::mat repmat1(stageframe_length, stageframe_length);
  arma::mat repmat2(stageframe_length, stageframe_length);
  
  for (int i = 0; i < stageframe_length; i++) {
    newstagevec(i) = stagevec(neworder(i));
    neworigsizevec(i) = origsizevec(neworder(i));
    neworigsizebvec(i) = origsizebvec(neworder(i));
    neworigsizecvec(i) = origsizecvec(neworder(i));
    newminagevec(i) = minagevec(neworder(i));
    newmaxagevec(i) = maxagevec(neworder(i));
    newrepvec(i) = repvec(neworder(i));
    newobsvec(i) = obsvec(neworder(i));
    newpropvec(i) = propvec(neworder(i));
    newimmvec(i) = immvec(neworder(i));
    newmatvec(i) = matvec(neworder(i));
    newrepentryvec(i) = repentryvec(neworder(i));
    newindvec(i) = indvec(neworder(i));
    newbinvec(i) = binvec(neworder(i));
    newbinbvec(i) = binbvec(neworder(i));
    newbincvec(i) = bincvec(neworder(i));
    newsizeminvec(i) = sizeminvec(neworder(i));
    newsizemaxvec(i) = sizemaxvec(neworder(i));
    newsizectrvec(i) = sizectrvec(neworder(i));
    newsizewidthvec(i) = sizewidthvec(neworder(i));
    newsizebminvec(i) = sizebminvec(neworder(i));
    newsizebmaxvec(i) = sizebmaxvec(neworder(i));
    newsizebctrvec(i) = sizebctrvec(neworder(i));
    newsizebwidthvec(i) = sizebwidthvec(neworder(i));
    newsizecminvec(i) = sizecminvec(neworder(i));
    newsizecmaxvec(i) = sizecmaxvec(neworder(i));
    newsizecctrvec(i) = sizecctrvec(neworder(i));
    newsizecwidthvec(i) = sizecwidthvec(neworder(i));
    newgroupvec(i) = groupvec(neworder(i));
    newcomvec(i) = comvec(neworder(i));
    newalive(i) = 1;
    
    repmat1.col(i) = repmatrix_true.col(neworder(i));
  }
  
  for (int i = 0; i < stageframe_length; i++) {
    repmat2.row(i) = repmat1.row(neworder(i));
  }
  
  if (format == 2) {
    newstagevec(stageframe_length) = "AlmostBorn";
    neworigsizevec(stageframe_length) = 0;
    neworigsizebvec(stageframe_length) = 0;
    neworigsizecvec(stageframe_length) = 0;
    newrepvec(stageframe_length) = 0;
    newobsvec(stageframe_length) = 1;
    newpropvec(stageframe_length) = 0;
    newimmvec(stageframe_length) = 1;
    newmatvec(stageframe_length) = 0;
    newrepentryvec(stageframe_length) = 0;
    newindvec(stageframe_length) = 1;
    newbinvec(stageframe_length) = 0;
    newbinbvec(stageframe_length) = 0;
    newbincvec(stageframe_length) = 0;
    newsizeminvec(stageframe_length) = 0;
    newsizemaxvec(stageframe_length) = 0;
    newsizectrvec(stageframe_length) = 0;
    newsizewidthvec(stageframe_length) = 0;
    newsizebminvec(stageframe_length) = 0;
    newsizebmaxvec(stageframe_length) = 0;
    newsizebctrvec(stageframe_length) = 0;
    newsizebwidthvec(stageframe_length) = 0;
    newsizecminvec(stageframe_length) = 0;
    newsizecmaxvec(stageframe_length) = 0;
    newsizecctrvec(stageframe_length) = 0;
    newsizecwidthvec(stageframe_length) = 0;
    newgroupvec(stageframe_length) = no_groups + 1;
    newcomvec(stageframe_length) = "Almost born (t-1)";
    newalive(stageframe_length) = 0;
    
    if (agemat) {
      newminagevec(stageframe_length) = 0;
      newmaxagevec(stageframe_length) = 0;
    } else {
      newminagevec(stageframe_length) = NA_REAL;
      newmaxagevec(stageframe_length) = NA_REAL;
    }
  }
  
  newstagevec(stageframe_length + (format - 1)) = "Dead";
  neworigsizevec(stageframe_length + (format - 1)) = 0;
  neworigsizebvec(stageframe_length + (format - 1)) = 0;
  neworigsizecvec(stageframe_length + (format - 1)) = 0;
  newrepvec(stageframe_length + (format - 1)) = 0;
  newobsvec(stageframe_length + (format - 1)) = 1;
  newpropvec(stageframe_length + (format - 1)) = 0;
  newimmvec(stageframe_length + (format - 1)) = 0;
  newmatvec(stageframe_length + (format - 1)) = 1;
  newrepentryvec(stageframe_length + (format - 1)) = 0;
  newindvec(stageframe_length + (format - 1)) = 1;
  newbinvec(stageframe_length + (format - 1)) = 0;
  newbinbvec(stageframe_length + (format - 1)) = 0;
  newbincvec(stageframe_length + (format - 1)) = 0;
  newsizeminvec(stageframe_length + (format - 1)) = 0;
  newsizemaxvec(stageframe_length + (format - 1)) = 0;
  newsizectrvec(stageframe_length + (format - 1)) = 0;
  newsizewidthvec(stageframe_length + (format - 1)) = 0;
  newsizebminvec(stageframe_length + (format - 1)) = 0;
  newsizebmaxvec(stageframe_length + (format - 1)) = 0;
  newsizebctrvec(stageframe_length + (format - 1)) = 0;
  newsizebwidthvec(stageframe_length + (format - 1)) = 0;
  newsizecminvec(stageframe_length + (format - 1)) = 0;
  newsizecmaxvec(stageframe_length + (format - 1)) = 0;
  newsizecctrvec(stageframe_length + (format - 1)) = 0;
  newsizecwidthvec(stageframe_length + (format - 1)) = 0;
  newgroupvec(stageframe_length + (format - 1)) = no_groups + 1;
  newcomvec(stageframe_length + (format - 1)) = "Dead";
  newalive(stageframe_length + (format - 1)) = 0;
  
  if (agemat) {
    newminagevec(stageframe_length + (format - 1)) = 0;
    newmaxagevec(stageframe_length + (format - 1)) = 0;
  } else {
    newminagevec(stageframe_length + (format - 1)) = NA_REAL;
    newmaxagevec(stageframe_length + (format - 1)) = NA_REAL;
  }
  
  arma::ivec stage_id;
  if (format == 1) {
    stage_id = linspace<arma::ivec>(1, (stageframe_length + 1), (stageframe_length + 1));
  } else {
    stage_id = linspace<arma::ivec>(1, (stageframe_length + 2), (stageframe_length + 2));
  }
  
  StringVector newstagevec_check = unique(newstagevec);
  if (newstagevec_check.length() < newstagevec.length()) {
    throw Rcpp::exception("All stage names must be unique.", false);
  }
  
  Rcpp::List newstageframe(32);
  
  newstageframe(0) = Rcpp::IntegerVector(stage_id.begin(), stage_id.end());
  newstageframe(1) = newstagevec;
  newstageframe(2) = neworigsizevec;
  newstageframe(3) = neworigsizebvec;
  newstageframe(4) = neworigsizecvec;
  newstageframe(5) = newminagevec;
  newstageframe(6) = newmaxagevec;
  newstageframe(7) = Rcpp::IntegerVector(newrepvec.begin(), newrepvec.end());
  newstageframe(8) = Rcpp::IntegerVector(newobsvec.begin(), newobsvec.end());
  newstageframe(9) = Rcpp::IntegerVector(newpropvec.begin(), newpropvec.end());
  newstageframe(10) = Rcpp::IntegerVector(newimmvec.begin(), newimmvec.end());
  newstageframe(11) = Rcpp::IntegerVector(newmatvec.begin(), newmatvec.end());
  newstageframe(12) = Rcpp::IntegerVector(newrepentryvec.begin(), newrepentryvec.end());
  newstageframe(13) = Rcpp::IntegerVector(newindvec.begin(), newindvec.end());
  newstageframe(14) = newbinvec;
  newstageframe(15) = newsizeminvec;
  newstageframe(16) = newsizemaxvec;
  newstageframe(17) = newsizectrvec;
  newstageframe(18) = newsizewidthvec;
  newstageframe(19) = newbinbvec;
  newstageframe(20) = newsizebminvec;
  newstageframe(21) = newsizebmaxvec;
  newstageframe(22) = newsizebctrvec;
  newstageframe(23) = newsizebwidthvec;
  newstageframe(24) = newbincvec;
  newstageframe(25) = newsizecminvec;
  newstageframe(26) = newsizecmaxvec;
  newstageframe(27) = newsizecctrvec;
  newstageframe(28) = newsizecwidthvec;
  newstageframe(29) = Rcpp::IntegerVector(newgroupvec.begin(), newgroupvec.end());
  newstageframe(30) = newcomvec;
  newstageframe(31) = Rcpp::IntegerVector(newalive.begin(), newalive.end());
  
  CharacterVector namevec = {"stage_id", "stage", "original_size", "original_size_b",
    "original_size_c", "min_age", "max_age", "repstatus", "obsstatus", "propstatus",
    "immstatus", "matstatus", "entrystage", "indataset", "binhalfwidth_raw", "sizebin_min",
    "sizebin_max", "sizebin_center", "sizebin_width", "binhalfwidthb_raw", "sizebinb_min",
    "sizebinb_max", "sizebinb_center", "sizebinb_width", "binhalfwidthc_raw", "sizebinc_min",
    "sizebinc_max", "sizebinc_center", "sizebinc_width", "group", "comments", "alive"};
  CharacterVector newclasses = {"data.frame", "stageframe"};
  newstageframe.attr("names") = namevec;
  newstageframe.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, (stageframe_length + format));
  newstageframe.attr("class") = newclasses;
  
  // From here on, we cover the new version of .overwrite_reassess
  if (overwrite.isNotNull()) {
    if (supp_provided) {
      throw Rcpp::exception("Please input either a supplement table or an overwrite table, not both.", false);
    }
    
    Rcpp::DataFrame supplement_thru(supplement);
    supplement_true = supplement_thru;
    
    stage3_supp = supplement_true["stage3"];
    stage2_supp = supplement_true["stage2"];
    stage1_supp = supplement_true["stage1"];
    eststage3_supp = supplement_true["eststage3"];
    eststage2_supp = supplement_true["eststage2"];
    eststage1_supp = supplement_true["eststage1"];
    givenrate_supp = supplement_true["givenrate"];
    convtype_supp = supplement_true["convtype"];
    convtype_t12_supp = supplement_true["convtype_t12"];
    
    supp_rows = givenrate_supp.length();
    multiplier_supp = Rcpp::NumericVector::create(1.0, supp_rows); // Not in overwrite()
    
    over_provided = true;
  }
  
  Rcpp::List newsupplement(10);
  
  if (over_provided || supp_provided) {
    StringVector unique_stages = unique(newstagevec);
    StringVector extra_terms = {"rep", "nrep", "immat", "mat", "prop", "npr", "all"};
    
    int no_newstages = unique_stages.length();
    int no_extraterms = extra_terms.length();
    
    StringVector all_possible_stage_terms(no_newstages + no_extraterms + no_groups);
    for (int i = 0; i < no_newstages; i++) {
      all_possible_stage_terms(i) = unique_stages(i);
    }
    for (int i = 0; i < no_extraterms; i++) {
      all_possible_stage_terms(i + no_newstages) = extra_terms(i);
    }
    for (int i = 0; i < no_groups; i++) {
      all_possible_stage_terms(i + no_newstages + no_extraterms) = group_text(i);
    }
    
    // Check for good entries in the supplement or overwrite table
    for (int i = 0; i < stage3_supp.length(); i++) {
      int s3supp_count {0};
      int s2supp_count {0};
      int s1supp_count {0};
      
      bool ests3_used = false;
      bool ests2_used = false;
      bool ests1_used = false;
      
      int ests3supp_count {0};
      int ests2supp_count {0};
      int ests1supp_count {0};
      
      for (int j = 0; j < all_possible_stage_terms.length(); j++) {
        if (stage3_supp(i) == all_possible_stage_terms(j)) s3supp_count++;
        if (stage2_supp(i) == all_possible_stage_terms(j)) s2supp_count++;
        
        if (!StringVector::is_na(eststage3_supp(i))) {
          ests3_used = true;
          if (eststage3_supp(i) == all_possible_stage_terms(j)) ests3supp_count++;
        }
        if (!StringVector::is_na(eststage2_supp(i))) {
          ests2_used = true;
          if (eststage2_supp(i) == all_possible_stage_terms(j)) ests2supp_count++;
        }
        
        if (historical) {
          if (stage1_supp(i) == all_possible_stage_terms(j)) s1supp_count++;
          
          if (!StringVector::is_na(eststage1_supp(i))) {
            ests1_used = true;
            if (eststage1_supp(i) == all_possible_stage_terms(j)) ests1supp_count++;
          }
        } 
      }
      
      if (s3supp_count == 0) {
        throw Rcpp::exception("Stage names in supplement or overwrite table (stage3) must match stageframe",
          false);
      }
      if (s2supp_count == 0) {
        throw Rcpp::exception("Stage names in supplement or overwrite table (stage2) must match stageframe",
          false);
      }
      if (ests3_used) {
        if (s3supp_count == 0) {
          throw Rcpp::exception("Stage names in supplement or overwrite table (eststage3) must match stageframe",
            false);
        }
      }
      if (ests2_used) {
        if (s2supp_count == 0) {
          throw Rcpp::exception("Stage names in supplement or overwrite table (eststage2) must match stageframe",
            false);
        }
      }
      if (historical) {
        if (s1supp_count == 0) {
          throw Rcpp::exception("Stage names in supplement or overwrite table (stage1) must match stageframe",
            false);
        }
        if (ests1_used) {
          if (s1supp_count == 0) {
            throw Rcpp::exception("Stage names in supplement or overwrite table (eststage1) must match stageframe",
              false);
          }
        }
      }
    }
    
    IntegerVector s1_calls (supp_rows, 1);
    IntegerVector s2_calls (supp_rows, 1);
    IntegerVector s3_calls (supp_rows, 1);
    IntegerVector ests1_calls (supp_rows, 1);
    IntegerVector ests2_calls (supp_rows, 1);
    IntegerVector ests3_calls (supp_rows, 1);
    IntegerVector s3_planned (supp_rows, 1);
    IntegerVector s2_planned (supp_rows, 1);
    IntegerVector s1_planned (supp_rows, 1);
      
    IntegerVector s123_calls (supp_rows, 1);
    
    // This section creates the new indices for the edited supplement/overwrite table
    arma::uvec newprop_stages = find(newpropvec);
    arma::uvec newprop0_stages = find(newpropvec == 0);
    arma::uvec newimm_stages = find(newimmvec);
    arma::uvec newalive_stages = find(newalive);
    arma::uvec newmat_stages1 = find(newmatvec);
    arma::uvec newmat_stages = intersect(newalive_stages, newmat_stages1);
    arma::uvec newrep_stages = find(newrepvec);
    arma::uvec newrep0_stages = find(newrepvec == 0);
    arma::uvec newmat_rep0_stages = intersect(newmat_stages, newrep0_stages);
    arma::uvec all_stages = find(newalive); // 7 "all"
    int no_current_group {0};
    
    // Now we build the expanded and edited supplement table
    for (int i = 0; i < supp_rows; i++) {
      if (stage3_supp(i) == "prop") {
        s3_calls(i) = newprop_stages.n_elem;
      } else if (stage3_supp(i) == "npr") {
        s3_calls(i) = newprop0_stages.n_elem;
      } else if (stage3_supp(i) == "immat") {
        s3_calls(i) = newimm_stages.n_elem;
      } else if (stage3_supp(i) == "mat") {
        s3_calls(i) = newmat_stages.n_elem;
      } else if (stage3_supp(i) == "rep") {
        s3_calls(i) = newrep_stages.n_elem;
      } else if (stage3_supp(i) == "nrep") {
        s3_calls(i) = newmat_rep0_stages.n_elem;
      } else if (stage3_supp(i) == "all") {
        s3_calls(i) = all_stages.n_elem;
      } else {
        for (int j = 0; j < no_groups; j++) {
          if (stage3_supp(i) == group_text(j)) {
            arma::uvec current_group = find(newgroupvec == j);
            no_current_group = current_group.n_elem;
            
            s3_calls(i) = no_current_group;
          }
        }
      }
      if (s3_calls(i) == 0) s3_calls(i) = 1;
      
      if (eststage3_supp(i) == "prop") {
        ests3_calls(i) = newprop_stages.n_elem;
      } else if (eststage3_supp(i) == "npr") {
        ests3_calls(i) = newprop0_stages.n_elem;
      } else if (eststage3_supp(i) == "immat") {
        ests3_calls(i) = newimm_stages.n_elem;
      } else if (eststage3_supp(i) == "mat") {
        ests3_calls(i) = newmat_stages.n_elem;
      } else if (eststage3_supp(i) == "rep") {
        ests3_calls(i) = newrep_stages.n_elem;
      } else if (eststage3_supp(i) == "nrep") {
        ests3_calls(i) = newmat_rep0_stages.n_elem;
      } else if (eststage3_supp(i) == "all") {
        ests3_calls(i) = all_stages.n_elem;
      } else {
        for (int j = 0; j < no_groups; j++) {
          if (eststage3_supp(i) == group_text(j)) {
            arma::uvec current_group = find(newgroupvec == j);
            no_current_group = current_group.n_elem;
            
            ests3_calls(i) = no_current_group;
          }
        }
      }
      if (ests3_calls(i) == 0) ests3_calls(i) = 1;
      
      if (stage2_supp(i) == "prop") {
        s2_calls(i) = newprop_stages.n_elem;
      } else if (stage2_supp(i) == "npr") {
        s2_calls(i) = newprop0_stages.n_elem;
      } else if (stage2_supp(i) == "immat") {
        s2_calls(i) = newimm_stages.n_elem;
      } else if (stage2_supp(i) == "mat") {
        s2_calls(i) = newmat_stages.n_elem;
      } else if (stage2_supp(i) == "rep") {
        s2_calls(i) = newrep_stages.n_elem;
      } else if (stage2_supp(i) == "nrep") {
        s2_calls(i) = newmat_rep0_stages.n_elem;
      } else if (stage2_supp(i) == "all") {
        s2_calls(i) = all_stages.n_elem;
      } else {
        for (int j = 0; j < no_groups; j++) {
          if (stage2_supp(i) == group_text(j)) {
            arma::uvec current_group = find(newgroupvec == j);
            no_current_group = current_group.n_elem;
            
            s2_calls(i) = no_current_group;
          }
        }
      }
      if (s2_calls(i) == 0) s2_calls(i) = 1;
      
      if (eststage2_supp(i) == "prop") {
        ests2_calls(i) = newprop_stages.n_elem;
      } else if (eststage2_supp(i) == "npr") {
        ests2_calls(i) = newprop0_stages.n_elem;
      } else if (eststage2_supp(i) == "immat") {
        ests2_calls(i) = newimm_stages.n_elem;
      } else if (eststage2_supp(i) == "mat") {
        ests2_calls(i) = newmat_stages.n_elem;
      } else if (eststage2_supp(i) == "rep") {
        ests2_calls(i) = newrep_stages.n_elem;
      } else if (eststage2_supp(i) == "nrep") {
        ests2_calls(i) = newmat_rep0_stages.n_elem;
      } else if (eststage2_supp(i) == "all") {
        ests2_calls(i) = all_stages.n_elem;
      } else {
        for (int j = 0; j < no_groups; j++) {
          if (eststage2_supp(i) == group_text(j)) {
            arma::uvec current_group = find(newgroupvec == j);
            no_current_group = current_group.n_elem;
            
            ests2_calls(i) = no_current_group;
          }
        }
      }
      if (ests2_calls(i) == 0) ests2_calls(i) = 1;
      
      if (stage1_supp(i) == "prop") {
        s1_calls(i) = newprop_stages.n_elem;
      } else if (stage1_supp(i) == "npr") {
        s1_calls(i) = newprop0_stages.n_elem;
      } else if (stage1_supp(i) == "immat") {
        s1_calls(i) = newimm_stages.n_elem;
      } else if (stage1_supp(i) == "mat") {
        s1_calls(i) = newmat_stages.n_elem;
      } else if (stage1_supp(i) == "rep") {
        s1_calls(i) = newrep_stages.n_elem;
      } else if (stage1_supp(i) == "nrep") {
        s1_calls(i) = newmat_rep0_stages.n_elem;
      } else if (stage1_supp(i) == "all") {
        s1_calls(i) = all_stages.n_elem;
      } else if (StringVector::is_na(stage1_supp(i))) {
        s1_calls(i) = 1;
      } else {
        for (int j = 0; j < no_groups; j++) {
          if (stage1_supp(i) == group_text(j)) {
            arma::uvec current_group = find(newgroupvec == j);
            no_current_group = current_group.n_elem;
            
            s1_calls(i) = no_current_group;
          }
        }
      }
      if (s1_calls(i) == 0) s1_calls(i) = 1;
      
      if (eststage1_supp(i) == "prop") {
        ests1_calls(i) = newprop_stages.n_elem;
      } else if (eststage1_supp(i) == "npr") {
        ests1_calls(i) = newprop0_stages.n_elem;
      } else if (eststage1_supp(i) == "immat") {
        ests1_calls(i) = newimm_stages.n_elem;
      } else if (eststage1_supp(i) == "mat") {
        ests1_calls(i) = newmat_stages.n_elem;
      } else if (eststage1_supp(i) == "rep") {
        ests1_calls(i) = newrep_stages.n_elem;
      } else if (eststage1_supp(i) == "nrep") {
        ests1_calls(i) = newmat_rep0_stages.n_elem;
      } else if (eststage1_supp(i) == "all") {
        ests1_calls(i) = all_stages.n_elem;
      } else if (StringVector::is_na(eststage1_supp(i))) {
        ests1_calls(i) = 1;
      } else {
        for (int j = 0; j < no_groups; j++) {
          if (eststage1_supp(i) == group_text(j)) {
            arma::uvec current_group = find(newgroupvec == j);
            no_current_group = current_group.n_elem;
            
            ests1_calls(i) = no_current_group;
          }
        }
      }
      if (ests1_calls(i) == 0) ests1_calls(i) = 1;
      
      if (!StringVector::is_na(eststage3_supp(i))) {
        if (eststage3_supp(i) != stage3_supp(i)) {
          if (s3_calls(i) == 1 && ests3_calls(i) > 1) {
            s3_planned(i) = ests3_calls(i);
          } else if (s3_calls(i) > 1 && ests3_calls(i) > 1) {
            throw Rcpp::exception("If stage group shorthand is used to designate both a transition and a proxy, then the shorthand group must be the same in both cases.", false);
          }
        } else {
          s3_planned(i) = s3_calls(i);
        }
      } else {
        s3_planned(i) = s3_calls(i);
      }
      
      if (!StringVector::is_na(eststage2_supp(i))) {
        if (eststage2_supp(i) != stage2_supp(i)) {
          if (s2_calls(i) == 1 && ests2_calls(i) > 1) {
            s2_planned(i) = ests2_calls(i);
          } else if (s2_calls(i) > 1 && ests2_calls(i) > 1) {
            throw Rcpp::exception("If stage group shorthand is used to designate both a transition and a proxy, then the shorthand group must be the same in both cases.", false);
          }
        } else {
          s2_planned(i) = s2_calls(i);
        }
      } else {
        s2_planned(i) = s2_calls(i);
      }
      
      if (!StringVector::is_na(eststage1_supp(i))) {
        if (historical && eststage1_supp(i) != stage1_supp(i)) {
          if (s1_calls(i) == 1 && ests1_calls(i) > 1) {
            s1_planned(i) = ests1_calls(i);
          } else if (s1_calls(i) > 1 && ests1_calls(i) > 1) {
            throw Rcpp::exception("If stage group shorthand is used to designate both a transition and a proxy, then the shorthand group must be the same in both cases.", false);
          }
        } else if (historical) {
          s1_planned(i) = s1_calls(i);
        } else if (!historical) {
          s1_planned(i) = 1;
        }
      } else {
        s1_planned(i) = s1_calls(i);
      }
      
      s123_calls(i) = s3_planned(i) * s2_planned(i) * s1_planned(i);
    }
    
    NumericVector basepoints(supp_rows, 0.0);
    for (int i = 0; i < (supp_rows - 1); i++) {
      basepoints(i+1) = basepoints(i) + s123_calls(i);
    }
    
    newsupp_rows = sum(s123_calls);
    
    StringVector stage3_newsupp(newsupp_rows);
    StringVector stage2_newsupp(newsupp_rows);
    StringVector stage1_newsupp(newsupp_rows);
    StringVector eststage3_newsupp(newsupp_rows);
    StringVector eststage2_newsupp(newsupp_rows);
    StringVector eststage1_newsupp(newsupp_rows);
    NumericVector givenrate_newsupp(newsupp_rows);
    IntegerVector convtype_newsupp(newsupp_rows);
    IntegerVector convtype_t12_newsupp(newsupp_rows);
    NumericVector multiplier_newsupp(newsupp_rows);
    
    int overall_counter {0};
    int group_check {0};
    
    int group_baseline3 {0};
    int group_baseline2 {0};
    int group_baseline1 {0};
    int group_baselinee3 {0};
    int group_baselinee2 {0};
    int group_baselinee1 {0};
    
    int group_ratchet3 {0};
    int group_ratchet2 {0};
    int group_ratchet1 {0};
    int group_ratchete3 {0};
    int group_ratchete2 {0};
    int group_ratchete1 {0};
    
    int prevl3 {0};
    int prevl2 {0};
    int prevl1 {0};
    int prevle3 {0};
    int prevle2 {0};
    int prevle1 {0};
    
    for (int i = 0; i < supp_rows; i++) {
      overall_counter = 0;
      for (int j = 0; j < s1_planned(i); j++) {
        for (int k = 0; k < s2_planned(i); k++) {
          for (int l = 0; l < s3_planned(i); l++) {
            if (stage3_supp(i) == "prop") {
              stage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop_stages(l));
            } else if (stage3_supp(i) == "npr") {
              stage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop0_stages(l));
            } else if (stage3_supp(i) == "immat") {
              stage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newimm_stages(l));
            } else if (stage3_supp(i) == "mat") {
              stage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_stages(l));
            } else if (stage3_supp(i) == "rep") {
              stage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newrep_stages(l));
            } else if (stage3_supp(i) == "nrep") {
              stage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_rep0_stages(l));
            } else if (stage3_supp(i) == "all") {
              stage3_newsupp(basepoints(i) + overall_counter) = newstagevec(all_stages(l));
            } else {
              for (int m = 0; m < no_groups; m++) {
                if (stage3_supp(i) == group_text(m)) {
                  if (l == 0) group_ratchet3 = 0;
                  if (l != prevl3 && l != 0) group_ratchet3 += 1;
                  
                  group_check = 1;
                  arma::uvec current_group = find(newgroupvec == m);
                  int current_group_length = current_group.n_elem;
                  if (group_ratchet3 > (current_group_length - 1)) {
                    group_ratchet3 = 0;
                  }
                  
                  if (group_ratchet3 == 0) {
                    group_baseline3 = l;
                  }
                  
                  stage3_newsupp(basepoints(i) + overall_counter) = 
                    newstagevec(current_group(l - group_baseline3));
                  
                  prevl3 = l;
                }
              }
              
              if (group_check == 0) {
               stage3_newsupp(basepoints(i) + overall_counter) = stage3_supp(i);
              }
              
              group_check = 0;
            }
            givenrate_newsupp(basepoints(i) + overall_counter) = givenrate_supp(i);
            multiplier_newsupp(basepoints(i) + overall_counter) = multiplier_supp(i);
            convtype_newsupp(basepoints(i) + overall_counter) = convtype_supp(i);
            convtype_t12_newsupp(basepoints(i) + overall_counter) = convtype_t12_supp(i);
            
            if (eststage3_supp(i) == "prop") {
              eststage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop_stages(l));
            } else if (eststage3_supp(i) == "npr") {
              eststage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop0_stages(l));
            } else if (eststage3_supp(i) == "immat") {
              eststage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newimm_stages(l));
            } else if (eststage3_supp(i) == "mat") {
              eststage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_stages(l));
            } else if (eststage3_supp(i) == "rep") {
              eststage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newrep_stages(l));
            } else if (eststage3_supp(i) == "nrep") {
              eststage3_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_rep0_stages(l));
            } else if (eststage3_supp(i) == "all") {
              eststage3_newsupp(basepoints(i) + overall_counter) = newstagevec(all_stages(l));
            } else {
              for (int m = 0; m < no_groups; m++) {
                if (eststage3_supp(i) == group_text(m)) {
                  if (l == 0) group_ratchete3 = 0;
                  if (l != prevle3 && l != 0) group_ratchete3 += 1;
                  
                  group_check = 1;
                  arma::uvec current_group = find(newgroupvec == m);
                  int current_group_length = current_group.n_elem;
                  if (group_ratchete3 > (current_group_length - 1)) {
                    group_ratchete3 = 0;
                  }
                  
                  if (group_ratchete3 == 0) {
                    group_baselinee3 = l;
                  }
                  
                  eststage3_newsupp(basepoints(i) + overall_counter) = 
                    newstagevec(current_group(l - group_baselinee3));
                  
                  prevle3 = l;
                }
              }
              
              if (group_check == 0) {
                eststage3_newsupp(basepoints(i) + overall_counter) = eststage3_supp(i);
              }
              
              group_check = 0;
            }
            
            if (stage2_supp(i) == "prop") {
              stage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop_stages(k));
            } else if (stage2_supp(i) == "npr") {
              stage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop0_stages(k));
            } else if (stage2_supp(i) == "immat") {
              stage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newimm_stages(k));
            } else if (stage2_supp(i) == "mat") {
              stage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_stages(k));
            } else if (stage2_supp(i) == "rep") {
              stage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newrep_stages(k));
            } else if (stage2_supp(i) == "nrep") {
              stage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_rep0_stages(k));
            } else if (stage2_supp(i) == "all") {
              stage2_newsupp(basepoints(i) + overall_counter) = newstagevec(all_stages(k));
            } else {
              for (int m = 0; m < no_groups; m++) {
                if (stage2_supp(i) == group_text(m)) {
                  if (k == 0) group_ratchet2 = 0;
                  if (k != prevl2 && k != 0) group_ratchet2 += 1;
                  
                  group_check = 1;
                  arma::uvec current_group = find(newgroupvec == m);
                  int current_group_length = current_group.n_elem;
                  if (group_ratchet2 > (current_group_length - 1)) {
                    group_ratchet2 = 0;
                  }
                  
                  if (group_ratchet2 == 0) {
                    group_baseline2 = k;
                  }
                  
                  stage2_newsupp(basepoints(i) + overall_counter) =
                    newstagevec(current_group(k - group_baseline2));
                  
                  prevl2 = k;
                }
              }
              
              if (group_check == 0) {
               stage2_newsupp(basepoints(i) + overall_counter) = stage2_supp(i);
              }
              
              group_check = 0;
            }
            
            if (eststage2_supp(i) == "prop") {
              eststage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop_stages(k));
            } else if (eststage2_supp(i) == "npr") {
              eststage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop0_stages(k));
            } else if (eststage2_supp(i) == "immat") {
              eststage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newimm_stages(k));
            } else if (eststage2_supp(i) == "mat") {
              eststage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_stages(k));
            } else if (eststage2_supp(i) == "rep") {
              eststage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newrep_stages(k));
            } else if (eststage2_supp(i) == "nrep") {
              eststage2_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_rep0_stages(k));
            } else if (eststage2_supp(i) == "all") {
              eststage2_newsupp(basepoints(i) + overall_counter) = newstagevec(all_stages(k));
            } else {
              for (int m = 0; m < no_groups; m++) {
                if (eststage2_supp(i) == group_text(m)) {
                  if (k == 0) group_ratchete2 = 0;
                  if (k != prevle2 && k != 0) group_ratchete2 += 1;
                  
                  group_check = 1;
                  arma::uvec current_group = find(newgroupvec == m);
                  int current_group_length = current_group.n_elem;
                  if (group_ratchete2 > (current_group_length - 1)) {
                    group_ratchete2 = 0;
                  }
                  
                  if (group_ratchete2 == 0) {
                    group_baselinee2 = k;
                  }
                  
                  eststage2_newsupp(basepoints(i) + overall_counter) =
                    newstagevec(current_group(k - group_baselinee2));
                  
                  prevle2 = k;
                }
              }
              
              if (group_check == 0) {
                eststage2_newsupp(basepoints(i) + overall_counter) = eststage2_supp(i);
              }
              
              group_check = 0;
            }
            
            if (stage1_supp(i) == "prop") {
              stage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop_stages(j));
            } else if (stage1_supp(i) == "npr") {
              stage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop0_stages(j));
            } else if (stage1_supp(i) == "immat") {
              stage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newimm_stages(j));
            } else if (stage1_supp(i) == "mat") {
              stage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_stages(j));
            } else if (stage1_supp(i) == "rep") {
              stage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newrep_stages(j));
            } else if (stage1_supp(i) == "nrep") {
              stage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_rep0_stages(j));
            } else if (stage1_supp(i) == "all") {
              stage1_newsupp(basepoints(i) + overall_counter) = newstagevec(all_stages(j));
            } else {
              for (int m = 0; m < no_groups; m++) {
                if (stage1_supp(i) == group_text(m)) {
                  if (j == 0) group_ratchet1 = 0;
                  if (j != prevl1 && j != 0) group_ratchet1 += 1;
                  
                  group_check = 1;
                  arma::uvec current_group = find(newgroupvec == m);
                  int current_group_length = current_group.n_elem;
                  if (group_ratchet1 > (current_group_length - 1)) {
                    group_ratchet1 = 0;
                  }
                  
                  if (group_ratchet1 == 0) {
                    group_baseline1 = j;
                  }
                  
                  stage1_newsupp(basepoints(i) + overall_counter) =
                    newstagevec(current_group(j - group_baseline1));
                  
                  prevl1 = j;
                }
              }
              
              if (group_check == 0) {
               stage1_newsupp(basepoints(i) + overall_counter) = stage1_supp(i);
              }
              
              group_check = 0;
            }
            
            if (eststage1_supp(i) == "prop") {
              eststage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop_stages(j));
            } else if (eststage1_supp(i) == "npr") {
              eststage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newprop0_stages(j));
            } else if (eststage1_supp(i) == "immat") {
              eststage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newimm_stages(j));
            } else if (eststage1_supp(i) == "mat") {
              eststage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_stages(j));
            } else if (eststage1_supp(i) == "rep") {
              eststage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newrep_stages(j));
            } else if (eststage1_supp(i) == "nrep") {
              eststage1_newsupp(basepoints(i) + overall_counter) = newstagevec(newmat_rep0_stages(j));
            } else if (eststage1_supp(i) == "all") {
              eststage1_newsupp(basepoints(i) + overall_counter) = newstagevec(all_stages(j));
            } else {
              for (int m = 0; m < no_groups; m++) {
                if (eststage1_supp(i) == group_text(m)) {
                  if (j == 0) group_ratchete1 = 0;
                  if (j != prevle1 && j != 0) group_ratchete1 += 1;
                  
                  group_check = 1;
                  arma::uvec current_group = find(newgroupvec == m);
                  int current_group_length = current_group.n_elem;
                  if (group_ratchete1 > (current_group_length - 1)) {
                    group_ratchete1 = 0;
                  }
                  
                  if (group_ratchete1 == 0) {
                    group_baselinee1 = j;
                  }
                  
                  eststage1_newsupp(basepoints(i) + overall_counter) =
                    newstagevec(current_group(j - group_baselinee1));
                  
                  prevle1 = j;
                }
              }
              
              if (group_check == 0) {
                eststage1_newsupp(basepoints(i) + overall_counter) = eststage1_supp(i);
              }
              
              group_check = 0;
            }
            
            overall_counter++;
          }
        }
      }
    }
    
    newsupplement(0) = stage3_newsupp;
    newsupplement(1) = stage2_newsupp;
    newsupplement(2) = stage1_newsupp;
    newsupplement(3) = eststage3_newsupp;
    newsupplement(4) = eststage2_newsupp;
    newsupplement(5) = eststage1_newsupp;
    newsupplement(6) = givenrate_newsupp;
    newsupplement(7) = multiplier_newsupp;
    newsupplement(8) = convtype_newsupp;
    newsupplement(9) = convtype_t12_newsupp;
  } else {
    newsupplement(0) = NULL;
    newsupplement(1) = NULL;
    newsupplement(2) = NULL;
    newsupplement(3) = NULL;
    newsupplement(4) = NULL;
    newsupplement(5) = NULL;
    newsupplement(6) = NULL;
    newsupplement(7) = NULL;
    newsupplement(8) = NULL;
    newsupplement(9) = NULL;
  }
  
  CharacterVector su_namevec = {"stage3", "stage2", "stage1", "eststage3", "eststage2",
    "eststage1", "givenrate", "multiplier", "convtype", "convtype_t12"};
  CharacterVector su_newclasses = {"data.frame", "lefkoSD"};
  newsupplement.attr("names") = su_namevec;
  newsupplement.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, newsupp_rows);
  newsupplement.attr("class") = su_newclasses;
  
  return Rcpp::List::create(_["stageframe"] = newstageframe, _["repmatrix"] = repmat2,
    _["ovtable"] = newsupplement);
}

