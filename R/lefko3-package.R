#' @title Historical and Ahistorical Population Projection Matrix Analysis
#' 
#' @description This package creates population matrix projection models (MPMs)
#' for use in population ecological analyses. Its specialty is the estimation
#' of historical MPMs, which are 2-dimensional matrices comprising 3 monitoring
#' occasions (2 time steps or periods) of demographic information. The package
#' constructs both function-based and raw MPMs for both standard ahistorical
#' (i.e. 2 occasions, 1 period) and historical analyses, and can also produce
#' age-by-stage MPMs and IPMs. It also includes powerful functions to
#' standardize demographic datasets.
#' 
#' @details The lefko3 package provides six categories of functions:
#' 
#' 1. Data transformation and handling functions
#' 
#' 2. Functions determining population characteristics from vertical data
#' 
#' 3. Model building and selection
#' 
#' 4. Matrix / integral projection model creation functions
#' 
#' 5. Population dynamics analysis functions
#' 
#' 6. Functions describing, summarizing, or visualizing MPMs and derived
#' structures
#' 
#' @details lefko3 also includes example datasets complete with sample code.
#' 
#' @docType package
#' @author Richard P. Shefferson <cdorm@g.ecc.u-tokyo.ac.jp>
#' @author Johan Ehrlén
#' @references Shefferson, R.P., J. Ehrlen, and S. Kurokawa. 2021. 
#' \emph{lefko3}: analyzing individual history through size-classified matrix 
#' population models. \emph{Methods in Ecology and Evolution} 12(2): 378-382.
#' @import Rcpp
#' @importFrom glmmTMB fixef glmmTMB nbinom2 ranef truncated_nbinom2 truncated_poisson
#' @importFrom lme4 fixef glmer lmer ranef VarCorr
#' @importFrom MASS glm.nb
#' @importFrom MuMIn dredge
#' @importFrom pscl zeroinfl
#' @importFrom Rcpp evalCpp
#' @importFrom SparseM as.matrix.csr image
#' @importFrom stats getCall glm lm na.action na.fail na.omit rnorm sd setNames xtabs
#' @importFrom stats as.formula median pchisq poisson var
#' @importFrom VGAM posnegbinomial pospoisson vglm
#' @useDynLib lefko3
#' @name lefko3
NULL
