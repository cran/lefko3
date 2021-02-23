#' @title Historical and Ahistorical Population Projection Matrix Analysis
#' 
#' @description This package creates population matrix projection models (MPMs)
#' for use in population ecological analyses. Its specialty is the estimation
#' of historical MPMs, which are 2-dimensional matrices comprising 3 time steps
#' of demographic information. The package constructs both function-based and
#' raw MPMs for both standard ahistorical (i.e. 2 time interval) and historical
#' analyses.
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
#' @importFrom Rcpp evalCpp
#' @importFrom stats na.action na.fail na.omit lm glm getCall xtabs sd rnorm setNames
#' @importFrom stats as.formula median pchisq var poisson
#' @importFrom pscl zeroinfl
#' @importFrom MuMIn dredge
#' @importFrom MASS glm.nb
#' @importFrom lme4 lmer glmer fixef ranef VarCorr
#' @importFrom glmmTMB glmmTMB nbinom2 fixef ranef
#' @useDynLib lefko3
#' @name lefko3
NULL
