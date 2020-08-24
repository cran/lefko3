#' Create Function-based Historical Population Projection Matrices
#' 
#' \code{flefko3()} returns a list of historical population projection 
#' matrices corresponding to the patches and years given, as well as the 
#' associated component transition and fecundity matrices, data frames 
#' detailing the characteristics of the ahistorical stages used and the
#' historical stage pairs created, and a data frame characterizing the patch 
#' and year combinations corresponding to these matrices. Note that, unlike 
#' \code{\link{rlefko3}()}, this function currently does not currently distinguish
#' populations.
#'
#' @param year A variable corresponding to year or observation time, or a set
#' of such values, given in values associated with the year term used in linear 
#' model development. Can also equal \code{all}, in which case matrices will
#' be estimated for all years. Defaults to \code{all}.
#' @param patch A variable designating which patches or subpopulations will
#' have matrices estimated. Should be set to specific patch names, or to 
#' \code{all} if all patches should have matrices estimated. Defaults to \code{all}.
#' @param stageframe A stageframe object that includes information on the size,
#' observation status, propagule status, immaturity status, and maturity status
#' of each ahistorical stage. Should also incorporate bin widths if size is
#' continuous.
#' @param repmatrix A matrix composed mostly of 0s, with non-zero values for
#' each potentially new individual (row) born to each reproductive stage
#' (column). Non-zero entries correspond to multipliers for fecundity, with 1
#' equaling full fecundity.
#' @param data The original historical demographic data frame used to 
#' estimate vital rates (class \code{hfvdata}). The original data frame is 
#' required in order to initialize years and patches properly.
#' @param modelsuite An optional suite of models of class \code{lefkoMod}.
#' If given, then \code{surv_model}, \code{obs_model}, \code{size_model}, \code{repst_model}, 
#' \code{fec_model}, \code{paramnames}, \code{yearcol}, and \code{patchcol} are not required. 
#' Note that the modeling exercise used to develop the models must have 
#' tested the impact of time \emph{t}-1 for this to work properly.
#' @param surv_model A linear model predicting survival probability. This can 
#' be a model of class \code{glm} or \code{glmer}, and requires a predicted binomial 
#' variable under a logit link. If given, then will overwrite any survival 
#' probability model given in \code{modelsuite}. This model must have been developed 
#' in a modeling exercise testing the impacts of times \emph{t} and \emph{t}-1.
#' @param obs_model A linear model predicting sprouting or observation
#' probability. This can be a model of class \code{glm} or \code{glmer}, and requires a 
#' predicted binomial variable under a logit link. If given, then 
#' will overwrite any observation probability model given in \code{modelsuite}.
#' This model must have been developed in a modeling exercise testing the impacts 
#' of times \emph{t} and \emph{t}-1.
#' @param size_model A linear model predicting size. This can be a model of class
#' \code{glm} or \code{glmer}, both of which require a predicted poisson variable under a 
#' log link, or a model of class \code{lm} or \code{lmer}, in which a Gaussian response is 
#' assumed. If given, then will overwrite any size model given in \code{modelsuite}.  
#' This model must have been developed in a modeling exercise testing the impacts 
#' of times \emph{t} and \emph{t}-1.
#' @param repst_model A linear model predicting reproduction probability. This 
#' can be a model of class \code{glm} or \code{glmer}, and requires a predicted binomial 
#' variable under a logit link. If given, then will overwrite any reproduction 
#' probability model given in \code{modelsuite}.  This model must have been developed 
#' in a modeling exercise testing the impacts of times \emph{t} and \emph{t}-1.
#' @param fec_model A linear model predicting fecundity. This can be a model
#' of class \code{glm} or \code{glmer}, and requires a predicted poisson variable under a 
#' log link. If given, then will overwrite any fecundity model given in 
#' \code{modelsuite}. This model must have been developed in a modeling exercise 
#' testing the impacts of times \emph{t} and \emph{t}-1.
#' @param jsurv_model A linear model predicting juvenile survival probability.
#' This can be a model of class \code{glm} or \code{glmer}, and requires a predicted 
#' binomial variable under a logit link. If given, then will overwrite any 
#' juvenile survival probability model given in \code{modelsuite}. This model must 
#' have been developed in a modeling exercise testing the impacts of times \emph{t} 
#' and \emph{t}-1.
#' @param jobs_model A linear model predicting juvenile sprouting or observation
#' probability. This can be a model of class \code{glm} or \code{glmer}, and requires a 
#' predicted binomial variable under a logit link. If given, then 
#' will overwrite any juvenile observation probability model given in \code{modelsuite}.
#' This model must have been developed in a modeling exercise testing the impacts 
#' of times \emph{t} and \emph{t}-1.
#' @param jsize_model A linear model predicting juvenile size. This can be a model
#' of class \code{glm} or \code{glmer}, both of which require a predicted poisson variable 
#' under a log link, or a model of class \code{lm} or \code{lmer}, in which a Gaussian 
#' response is assumed. If given, then will overwrite any juvenile size model 
#' given in \code{modelsuite}. This model must have been developed in a modeling 
#' exercise testing the impacts of times \emph{t} and \emph{t}-1.
#' @param jrepst_model A linear model predicting reproduction probability of a 
#' mature individual that was immature in the previous year. This can be a model 
#' of class \code{glm} or \code{glmer}, and requires a predicted binomial variable under a 
#' logit link. If given, then will overwrite any reproduction probability model 
#' given in \code{modelsuite}. This model must have been developed in a modeling 
#' exercise testing the impacts of times \emph{t} and \emph{t}-1.
#' @param paramnames A dataframe with two columns, the first showing the general
#' model terms that will be used in matrix creation, and the second showing the
#' equivalent terms used in modeling. Only required if `modelsuite` is not 
#' supplied.
#' @param surv_dev A numeric value to be added to the y-intercept in the linear
#' model for survival probability.
#' @param obs_dev A numeric value to be added to the y-intercept in the linear
#' model for observation probability.
#' @param size_dev A numeric value to be added to the y-intercept in the linear
#' model for size.
#' @param repst_dev A numeric value to be added to the y-intercept in the linear
#' model for probability of reproduction.
#' @param fec_dev A numeric value to be added to the y-intercept in the linear
#' model for fecundity.
#' @param jsurv_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile survival probability.
#' @param jobs_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile observation probability.
#' @param jsize_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile size.
#' @param jrepst_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile reproduction probability.
#' @param repmod A scalar multiplier of fecundity. Defaults to 1.
#' @param overwrite A data frame developed with the \code{\link{overwrite}()} function,
#' describing transitions to be overwritten either with given values or 
#' with other estimated transitions.
#' @param yearcol The variable name or column number corresponding to year 
#' in time \emph{t} in the dataset. Not needed if a \code{modelsuite} is supplied.
#' @param patchcol The variable name or column number corresponding to patch in 
#' the dataset. Not needed if a \code{modelsuite} is supplied.
#' @param year.as.random A logical term indicating whether coefficients for
#' missing patches within vital rate models should be estimated as random
#' intercepts. Defaults to FALSE, in which case missing time step coefficients 
#' are set to 0.
#' @param patch.as.random A logical term indicating whether coefficients for
#' missing patches within vital rate models should be estimated as random
#' intercepts. Defaults to FALSE, in which case missing patch coefficients 
#' are set to 0.
#' @param randomseed A numeric value used as a seed to generate random 
#' estimates for missing time step and patch coefficients, if either
#' \code{year.as.random} or \code{patch.as.random} is set to TRUE.
#' @param negfec A logical value denoting whether fecundity values estimated
#' to be negative should be reset to 0. Defaults to FALSE.
#' @param reduce A logical value denoting whether to remove historical stages
#' associated solely with 0 transitions. These are only removed in cases where
#' the associated row and column sums in ALL matrices estimated equal 0. 
#' Defaults to FALSE.
#'
#' @return If all inputs are properly formatted, then this function will return
#' either an object of class \code{lefkoMat}. Output includes:
#' 
#' \item{A}{A list of full projection matrices in order of sorted patches
#' and years.}
#' \item{U}{A list of survival-transition matrices sorted as in \code{A}.}
#' \item{F}{A list of fecundity matrices sorted as in \code{A}.}
#' \item{hstages}{A data frame matrix showing the pairing of ahistorical 
#' stages used to create historical stage pairs.}
#' \item{ahstages}{A data frame detailing the characteristics of associated
#' ahistorical stages.}
#' \item{labels}{A data frame giving the patch and year of each matrix in
#' order. In \code{flefko3()}, only one population may be analyzed at once, and 
#' so \code{pop = NA}.}
#' \item{matrixqc}{A short vector describing the number of non-zero elements
#' in \code{U} and \code{F} matrices, and the number of annual matrices.}
#' \item{modelqc}{This is the \code{qc} portion of the \code{modelsuite} input.}
#'
#' @examples
#' \donttest{
#' data(lathyrus)
#' 
#' sizevector <- c(0, 4.6, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9)
#' stagevector <- c("Sd", "Sdl", "Dorm", "Sz1nr", "Sz2nr", "Sz3nr", "Sz4nr", "Sz5nr",
#'                  "Sz6nr", "Sz7nr", "Sz8nr", "Sz9nr", "Sz1r", "Sz2r", "Sz3r", "Sz4r",
#'                  "Sz5r", "Sz6r", "Sz7r", "Sz8r", "Sz9r")
#' repvector <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 4.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
#'             0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
#' 
#' lathframeln <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector,
#'                          obsstatus = obsvector, matstatus = matvector, immstatus = immvector,
#'                          indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvertln <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT",
#'                            individcol = "GENET", blocksize = 9, juvcol = "Seedling1988",
#'                            size1col = "lnVol88", repstr1col = "FCODE88",
#'                            fec1col = "Intactseed88", dead1col = "Dead1988",
#'                            nonobs1col = "Dormant1988", stageassign = lathframeln,
#'                            stagesize = "sizea", censorcol = "Missing1988",
#'                            censorkeep = NA, NAas0 = TRUE, censor = TRUE)
#' 
#' lathvertln$feca2 <- round(lathvertln$feca2)
#' lathvertln$feca1 <- round(lathvertln$feca1)
#' lathvertln$feca3 <- round(lathvertln$feca3)
#' 
#' lathrepmln <- matrix(0, 21, 21)
#' lathrepmln[1, c(13:21)] <- 0.345
#' lathrepmln[2, c(13:21)] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"),
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' lathmodelsln3 <- modelsearch(lathvertln, historical = TRUE, approach = "lme4", suite = "main", 
#'                              vitalrates = c("surv", "obs", "size", "repst", "fec"), 
#'                              juvestimate = "Sdl",bestfit = "AICc&k", sizedist = "gaussian", 
#'                              fecdist = "poisson", indiv = "individ", patch = "patchid", 
#'                              year = "year2",year.as.random = TRUE, patch.as.random = TRUE,
#'                              show.model.tables = TRUE, quiet = TRUE)
#' 
#' lathmat3ln <- flefko3(year = "all", patch = "all", stageframe = lathframeln, 
#'                       modelsuite = lathmodelsln3, data = lathvertln, 
#'                       repmatrix = lathrepmln, overwrite = lathover3,
#'                       patchcol = "patchid", yearcol = "year2", year.as.random = FALSE,
#'                       patch.as.random = FALSE, reduce = FALSE)
#' summary(lathmat3ln)
#' }
#' 
#' @export
flefko3 <- function(year = "all", patch = "all", stageframe, repmatrix = NA, data = NA, 
                    modelsuite = NA, surv_model = NA, obs_model = NA, size_model = NA, 
                    repst_model = NA, fec_model = NA, jsurv_model = NA, jobs_model = NA, 
                    jsize_model = NA, jrepst_model = NA, paramnames = NA, surv_dev = 0, 
                    obs_dev = 0, size_dev = 0, repst_dev = 0, fec_dev = 0, jsurv_dev = 0, 
                    jobs_dev = 0, jsize_dev = 0, jrepst_dev = 0, repmod = 1, overwrite = NA, 
                    yearcol = NA, patchcol = NA, year.as.random = FALSE, patch.as.random = FALSE, 
                    randomseed = 0, negfec = FALSE, reduce = FALSE) {
  
  if (all(is.na(modelsuite)) & all(is.na(paramnames))) {
    warning("Function may not work properly without a dataframe of model parameters or equivalents supplied either through the modelsuite option or through the paramnames input parameter.")
  } else if (!all(is.na(modelsuite))) {
    paramnames <- modelsuite$paramnames
    yearcol <- paramnames$modelparams[which(paramnames$mainparams == "year2")]
    patchcol <- paramnames$modelparams[which(paramnames$mainparams == "patch")]
  }
  
  if (all(is.na(data))) {stop("Need original vertical dataset to set proper limits on year and patch.", call. = FALSE)}
  if (!any(class(data) == "data.frame")) {stop("Need original vertical dataset used in modeling to proceed.", call. = FALSE)}
  
  if (is.character(yearcol)) {
    choicevar <- which(names(data) == yearcol);
    mainyears <- sort(unique(data[,choicevar]))
  } else if (is.numeric(yearcol)) {
    mainyears <- sort(unique(data[, yearcol]));
  } else {
    stop("Need appropriate year column designation.", call. = FALSE)
  }
  
  if (any(is.character(year))) {
    if (is.element(tolower(year), "all")) {
      year <- mainyears
    } else {
      stop("Year designation not recognized.", call. = FALSE)
    }
  }
  
  if (length(year) == 0 | all(is.na(year) == TRUE) | any(is.na(year))) {
    stop("This function cannot proceed without being given a specific year, or a suite of years.", call. = FALSE)
  }
  
  if (is.na(patch) & !is.na(patchcol)) {
    warning("Matrix creation may not proceed properly without input in the patch option when using a modelsuite in which patch is designated.")
  }
  
  if (is.character(patchcol)) {
    choicevar <- which(names(data) == patchcol);
    mainpatches <- sort(unique(as.character(data[,choicevar])))
  } else if (is.numeric(patchcol)) {
    mainpatches <- sort(unique(as.character(data[, patchcol])));
  } else {
    mainpatches <- NA
  }
  
  if (any(is.character(patch))) {
    if (is.element(tolower(patch), "all")) {
      patch <- mainpatches
    } else if (!is.element(patch, mainpatches)) {
      stop("Patch designation not recognized.", call. = FALSE)
    }
  }
  
  if (all(is.na(repmatrix))) {
    repmatrix <- matrix(0, dim(stageframe)[1], dim(stageframe)[1])
    repstages <- which(stageframe$repstatus == 1)
    repmatrix[1, repstages] <- 1
  }
  
  if (any(!suppressWarnings(!is.na(as.numeric(as.character(stageframe$bin_size_ctr)))))) {
    stop("Function flefko3() requires size to be numeric rather than categorical.", call. = FALSE)
  }
  
  melchett <- .sf_reassess(stageframe, repmatrix, overwrite)
  stageframe <- melchett[[1]]
  repmatrix <- melchett[[2]]
  
  stageframe$new_stage_id <- as.numeric(stageframe$new_stage_id)
  stageframe$original_size <- as.numeric(stageframe$original_size)
  stageframe$bin_size_min <- as.numeric(stageframe$bin_size_min)
  stageframe$bin_size_ctr <- as.numeric(stageframe$bin_size_ctr)
  stageframe$bin_size_max <- as.numeric(stageframe$bin_size_max)
  stageframe$bin_size_width <- as.numeric(stageframe$bin_size_width)
  stageframe$bin_raw_halfwidth <- as.numeric(stageframe$bin_raw_halfwidth)
  stageframe$repstatus <- as.numeric(stageframe$repstatus)
  stageframe$obsstatus <- as.numeric(stageframe$obsstatus)
  stageframe$propstatus <- as.numeric(stageframe$propstatus)
  stageframe$immstatus <- as.numeric(stageframe$immstatus)
  stageframe$matstatus <- as.numeric(stageframe$matstatus)
  stageframe$indataset <- as.numeric(stageframe$indataset)
  stageframe$alive <- as.numeric(stageframe$alive)
  
  stageframe$stageno <- c(1:dim(stageframe)[1])
  rownames(stageframe) <- stageframe$stageno
  stageframe$fullstage <- apply(as.matrix(c(1:dim(stageframe)[1])), 1, function(X) {
    paste(stageframe$bin_size_ctr[X], stageframe$repstatus[X])
  })
  
  instages <- length(stageframe$new_stage_id)
  
  allstages.size <- expand.grid(size3 = stageframe$bin_size_ctr, size2n = stageframe$bin_size_ctr, 
                                size2o = stageframe$bin_size_ctr, size1 = stageframe$bin_size_ctr)
  allstages.obs <- expand.grid(obs3 = stageframe$obsstatus, obs2n = stageframe$obsstatus, 
                               obs2o = stageframe$obsstatus, obs1 = stageframe$obsstatus)
  allstages.rep <- expand.grid(rep3 = stageframe$repstatus, rep2n = stageframe$repstatus, 
                               rep2o = stageframe$repstatus, rep1 = stageframe$repstatus)
  allstages.mat <- expand.grid(mat3 = stageframe$matstatus, mat2n = stageframe$matstatus,
                               mat2o = stageframe$matstatus, mat1 = stageframe$matstatus)
  allstages.imm <- expand.grid(imm3 = stageframe$immstatus, imm2n = stageframe$immstatus,
                               imm2o = stageframe$immstatus, imm1 = stageframe$immstatus)
  allstages.re32 <- as.vector(apply(as.matrix(1:dim(cbind(rbind(repmatrix, 0), 0))[1]), 
                                    1, function(X) {as.vector(cbind(rbind(repmatrix, 0), 0)[,X])}))
  allstages.re33 <- expand.grid(allstages.re32, allstages.re32)
  allstages.ind <- expand.grid(indata3 = stageframe$indataset, indata2n = stageframe$indataset,
                               indata2o = stageframe$indataset, indata1 = stageframe$indataset)
  allstages.stages <- expand.grid(stage3 = stageframe$new_stage_id, stage2n = stageframe$new_stage_id,
                                  stage2o = stageframe$new_stage_id, stage1 = stageframe$new_stage_id)
  allstages.bins <- rep(stageframe$bin_size_width, (instages^3))
  allstages <- cbind.data.frame(allstages.stages, allstages.size, allstages.obs, allstages.rep, 
                                allstages.mat, allstages.imm, allstages.re33[,1], allstages.ind,
                                allstages.bins)
  names(allstages) <- c("stage3", "stage2n", "stage2o", "stage1", "size3", "size2n", "size2o", "size1", 
                        "obs3", "obs2n", "obs2o", "obs1", "rep3", "rep2n", "rep2o", "rep1", 
                        "mat3", "mat2n", "mat2o", "mat1", "imm3", "imm2n", "imm2o", "imm1",
                        "repentry3", "indata3", "indata2n", "indata2o", "indata1", "binwidth")
  
  allstages$minage3 <- 0
  allstages$minage2 <- 0
  allstages$maxage3 <- 0
  allstages$maxage2 <- 0
  allstages$actualage2 <- 0
  
  allstages$index321 <- allstages$stage3 + (allstages$stage2n * instages) + (allstages$stage1 * instages * instages)
  
  if (!all(is.na(overwrite))) {
    
    overwrite <- .overwrite_reassess(overwrite, stageframe, historical = TRUE)
    if (dim(overwrite)[1] == 0) overwrite <- NA
  }
  
  if (!all(is.na(overwrite))) {
    
    overwrite$index3 <- apply(as.matrix(c(1:dim(overwrite)[1])), 1, function(X) {
      stageframe$stageno[which(stageframe$orig_stage_id == overwrite$stage3[X])]})
    overwrite$index2 <- apply(as.matrix(c(1:dim(overwrite)[1])), 1, function(X) {
      stageframe$stageno[which(stageframe$orig_stage_id == overwrite$stage2[X])]})
    overwrite$index1 <- apply(as.matrix(c(1:dim(overwrite)[1])), 1, function(X) {
      stageframe$stageno[which(stageframe$orig_stage_id == overwrite$stage1[X])]})
    overwrite$new3 <- apply(as.matrix(c(1:dim(overwrite)[1])), 1, function(X) {
      if(!is.na(overwrite$eststage3[X])) {
        stageframe$stageno[which(stageframe$orig_stage_id == overwrite$eststage3[X])]
      } else {
        overwrite$index3[X]
      }
    })
    
    overwrite$new2 <- apply(as.matrix(c(1:dim(overwrite)[1])), 1, function(X) {
      if(!is.na(overwrite$eststage2[X])) {
        stageframe$stageno[which(stageframe$orig_stage_id == overwrite$eststage2[X])]
      } else {
        overwrite$index2[X]
      }
    })
    
    overwrite$new1 <- apply(as.matrix(c(1:dim(overwrite)[1])), 1, function(X) {
      if(!is.na(overwrite$eststage1[X])) {
        stageframe$stageno[which(stageframe$orig_stage_id == overwrite$eststage1[X])]
      } else {
        overwrite$index1[X]
      }
    })
    
    overwrite$indexold321 <- overwrite$index3 + (overwrite$index2 * instages) + (overwrite$index1 * instages * instages)
    overwrite$indexnew321 <- overwrite$new3 + (overwrite$new2 * instages) + (overwrite$new1 * instages * instages)
    overwrite$givenrate[which(is.na(overwrite$givenrate))] <- 0
    
    overwrite$new3[which(is.na(overwrite$eststage3))] <- -1
    overwrite$givenrate[which(is.na(overwrite$givenrate))] <- -1
    overwrite$givenrate[which(overwrite$givenrate == 0)] <- -1
    
    allstageadditions <- ovreplace(allstages$index321, overwrite$indexold321, overwrite$indexnew321, 
                                   overwrite$convtype, overwrite$new3, overwrite$givenrate)
    
    allstages$ovest_t <- allstageadditions[,1]
    allstages$ovgiven_t <- allstageadditions[,2]
    allstages$ovest_f <- allstageadditions[,3]
    allstages$ovgiven_f <- allstageadditions[,4]
    
  } else {
    
    allstages$ovgiven_t <- -1
    allstages$ovest_t <- -1
    allstages$ovgiven_f <- -1
    allstages$ovest_f <- -1
    
  }
  
  if (class(modelsuite) == "lefkoMod") {
    if(is.na(surv_model)) {surv_model <- modelsuite$survival_model}
    if(is.na(obs_model)) {obs_model <- modelsuite$observation_model}
    if(is.na(size_model)) {size_model <- modelsuite$size_model}
    if(is.na(repst_model)) {repst_model <- modelsuite$repstatus_model}
    if(is.na(fec_model)) {fec_model <- modelsuite$fecundity_model}
    
    if(is.na(jsurv_model)) {jsurv_model <- modelsuite$juv_survival_model}
    if(is.na(jobs_model)) {jobs_model <- modelsuite$juv_observation_model}
    if(is.na(jsize_model)) {jsize_model <- modelsuite$juv_size_model}
    if(is.na(jrepst_model)) {jrepst_model <- modelsuite$juv_reproduction_model}
  }
  
  set.seed(randomseed)
  
  surv_proxy <- .modelextract(surv_model, paramnames, mainyears, mainpatches)
  obs_proxy <- .modelextract(obs_model, paramnames, mainyears, mainpatches)
  
  sigma <- 0
  rvarssummed <- 0
  sizedist <- 1
  
  size_proxy <- .modelextract(size_model, paramnames, mainyears, mainpatches)
  
  if (size_proxy$family == "poisson") {
    sizedist <- 0
    if (!all(is.na(size_proxy$variances))) {
      rvarssummed <- sum(size_proxy$variances[,"vcov"])
    } else {
      rvarssummed <- 0
    }
  } else if (size_proxy$family == "gaussian") {
    sizedist <- 2
    sigma <- size_proxy$sigma
  } else {
    sizedist <- 1
  }
  
  repst_proxy <- .modelextract(repst_model, paramnames, mainyears, mainpatches)
  fec_proxy <- .modelextract(fec_model, paramnames, mainyears, mainpatches)
  
  if (fec_proxy$family == "poisson") {
    fecdist <- 0
  } else if (fec_proxy$family == "gaussian") {
    fecdist <- 2
  } else {
    fecdist <- 1
  }
  
  jsurv_proxy <- .modelextract(jsurv_model, paramnames, mainyears, mainpatches)
  jobs_proxy <- .modelextract(jobs_model, paramnames, mainyears, mainpatches)
  
  jsigma <- 0
  jrvarssummed <- 0
  
  jsize_proxy <- .modelextract(jsize_model, paramnames, mainyears, mainpatches)
  
  if (jsize_proxy$family == "poisson") {
    if (!all(is.na(jsize_proxy$variances))) {
      jrvarssummed <- sum(jsize_proxy$variances[,"vcov"])
    } else {
      jrvarssummed <- 0
    }
  } else if (jsize_proxy$family == "gaussian") {
    jsigma <- jsize_proxy$sigma
  }
  
  jrepst_proxy <- .modelextract(jrepst_model, paramnames, mainyears, mainpatches)
  
  if (!all(is.na(patch))) {
    listofyears <- apply(as.matrix(patch), 1, function(X) {
      output <- cbind.data.frame(NA, X, as.matrix(year));
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    listofyears <- do.call(rbind.data.frame, listofyears)
    listofyears$poporder <- NA
    listofyears$patchorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainpatches == listofyears$patch[X])})
    listofyears$yearorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainyears == listofyears$year2[X])})
    
  } else {
    
    listofyears <- cbind.data.frame(NA, NA, as.matrix(year))
    names(listofyears) <- c("pop", "patch", "year2")
    
    listofyears$poporder <- NA
    listofyears$patchorder <- 1
    listofyears$yearorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainyears == listofyears$year2[X])})
  }
  
  maxsize <- max(c(allstages$size3, allstages$size2n, allstages$size2o, allstages$size1), na.rm = TRUE)
  allstages$indata <- allstages$indata3 * allstages$indata2n * allstages$indata2o * allstages$indata1
  aliverows <- hoffmannofstuttgart(0, as.matrix(allstages))
  
  allstages <- allstages[(aliverows + 1),]
  
  aliveandequal <- hoffmannofstuttgart(1, as.matrix(allstages))
  
  total.matrix.dim <- (length(stageframe$bin_size_ctr) - 1)^2
  patch.elements <- total.matrix.dim^2
  
  yearlist <- split(listofyears, seq(nrow(listofyears)))
  
  allstages <- allstages[(aliveandequal + 1),]
  allstages$aliveandequal <- as.vector(aliveandequal)
  allstages$r_aliveandequal <- as.vector(aliveandequal) + 1
  
  madsexmadrigal <- lapply(yearlist, jerzeibalowski, allstages, surv_proxy, obs_proxy, size_proxy, repst_proxy,
                           fec_proxy, jsurv_proxy, jobs_proxy, jsize_proxy, jrepst_proxy, surv_dev, obs_dev, 
                           size_dev, repst_dev, fec_dev, jsurv_dev, jobs_dev, jsize_dev, jrepst_dev, patch.elements,
                           total.matrix.dim, repmod, rvarssummed, sigma, jrvarssummed, jsigma, maxsize, sizedist, 
                           fecdist, negfec)
  
  a_list <- lapply(madsexmadrigal, function(X) {X$A})
  u_list <- lapply(madsexmadrigal, function(X) {X$U})
  f_list <- lapply(madsexmadrigal, function(X) {X$F})
  
  ahstages <- stageframe[1:(dim(stageframe)[1] - 1),]
  
  pairings1 <- expand.grid(stcod3 = ahstages$new_stage_id, stcod2 = ahstages$new_stage_id)
  pairings2 <- expand.grid(stage3 = ahstages$orig_stage_id, stage2 = ahstages$orig_stage_id)
  hstages <- cbind.data.frame(pairings2, pairings1)
  
  qcoutput1 <- NA
  qcoutput2 <- NA
  
  totalutransitions <- sum(unlist(lapply(u_list, function(X) {length(which(X != 0))})))
  totalftransitions <- sum(unlist(lapply(f_list, function(X) {length(which(X != 0))})))
  totalmatrices <- length(u_list)
  
  qcoutput1 <- c(totalutransitions, totalftransitions, totalmatrices)
  
  if (is.element("qc", names(modelsuite))) {qcoutput2 <- modelsuite$qc}
  
  if (reduce == TRUE) {
    drops <- .reducer3(a_list, u_list, f_list, hstages)
    
    a_list <- drops$A
    u_list <- drops$U
    f_list <- drops$F
    hstages <- drops$hstages
  }
  
  output <- list(A = a_list, U = u_list, F = f_list, hstages = hstages, ahstages = ahstages, 
                 labels = listofyears[,c(1:3)], matrixqc = qcoutput1, modelqc = qcoutput2)
  
  class(output) <- "lefkoMat"
  
  return(output)
}

#' Reduce Matrix Dimensions By Eliminating Empty Stages
#' 
#' \code{.reducer3()} identifies empty stages in a set of historical matrices and
#' removes them from all matrices. It is used within \code{\link{flefko3}()} and 
#' \code{\link{rlefko3}()}.
#' 
#' @param A List of population projection matrices, from a \code{lefkoMat} object.
#' @param U List of surviva-transition matrices corresponding to \code{A}.
#' @param F List of fecundity matrices corresponding to \code{A}.
#' @param hstages Data frame giving the names and identities of historical stage 
#' pairs used to create matrices.
#' 
#' @return Returns a list of reduced \code{A}, \code{U}, and \code{F} matrices, plus the reduced
#' \code{hstages} object.
#' 
#' @keywords internal
#' @noRd
.reducer3 <- function(A, U, F, hstages) {
  stagepatterns <- lapply(A, function(X) {
    matrix.sums <- colSums(X) + rowSums(X)
    return(matrix.sums)
  })
  
  used.stages.mat <- do.call("rbind", stagepatterns)
  used.stages.ovr <- colSums(used.stages.mat)
  keep.stages <- which(used.stages.ovr > 0)
  
  Ared <- lapply(A, function(X) {
    return(X[keep.stages, keep.stages])
  })
  
  Ured <- lapply(U, function(X) {
    return(X[keep.stages, keep.stages])
  })
  
  Fred <- lapply(F, function(X) {
    return(X[keep.stages, keep.stages])
  })
  
  hstred <- hstages[keep.stages,]
  
  return(list(A = Ared, U = Ured, F = Fred, hstages = hstred))
}

#' Create Function-based Ahistorical Population Projection Matrices
#'
#' \code{flefko2()} returns a list of ahistorical population projection 
#' matrices corresponding to the patches and years given, as well as the 
#' associated component transition and fecundity matrices, a data frame 
#' detailing the characteristics of ahistorical stages, and a data frame
#' characterizing the patch and year combinations corresponding to these 
#' matrices. Note that, unlike \code{\link{rlefko2}()}, this function currently does
#' not currently distinguish populations.
#'
#' @param year A variable corresponding to year or observation time, or a set
#' of such values, given in values associated with the year term used in linear 
#' model development. Can also equal \code{all}, in which case matrices will
#' be estimated for all years. Defaults to \code{all}.
#' @param patch A variable designating which patches or subpopulations will
#' have matrices estimated. Should be set to specific patch names, or to 
#' `"all"` if all patches should have matrices estimated. Defaults to \code{all}.
#' @param stageframe A stageframe object that includes information on the size,
#' observation status, propagule status, immaturity status, and maturity status
#' of each ahistorical stage. Should also incorporate bin widths if size is
#' continuous.
#' @param repmatrix A matrix composed mostly of 0s, with non-zero values for
#' each potentially new individual (row) born to each reproductive stage
#' (column). Non-zero entries correspond to multipliers for fecundity, with 1
#' equaling full fecundity.
#' @param modelsuite An optional suite of models of class \code{lefkoMod}.
#' If given, then \code{surv_model}, \code{obs_model}, \code{size_model}, \code{repst_model}, 
#' \code{fec_model}, \code{paramnames}, \code{yearcol}, and \code{patchcol} are not required. 
#' Note that the modeling exercise used to develop the models must have 
#' tested only the impact of time \emph{t} for this to work properly.
#' @param surv_model A linear model predicting survival probability. This can 
#' be a model of class \code{glm} or \code{glmer}, and requires a predicted binomial 
#' variable under a logit link. If given, then will overwrite any survival 
#' probability model given in \code{modelsuite}. This model must have been developed 
#' in a modeling exercise testing only the impacts of time \emph{t}.
#' @param obs_model A linear model predicting sprouting or observation
#' probability. This can be a model of class \code{glm} or \code{glmer}, and requires a 
#' predicted binomial variable under a logit link. If given, then 
#' will overwrite any observation probability model given in \code{modelsuite}.
#' This model must have been developed in a modeling exercise testing only
#' the impacts of time \emph{t}.
#' @param size_model A linear model predicting size. This can be a model of class
#' \code{glm} or \code{glmer}, both of which require a predicted poisson variable under a 
#' log link, or a model of class \code{lm} or \code{lmer}, in which a Gaussian response is 
#' assumed. If given, then will overwrite any size model given in \code{modelsuite}.  
#' This model must have been developed in a modeling exercise testing only
#' the impacts of time \emph{t}.
#' @param repst_model A linear model predicting reproduction probability. This 
#' can be a model of class \code{glm} or \code{glmer}, and requires a predicted binomial 
#' variable under a logit link. If given, then will overwrite any reproduction 
#' probability model given in \code{modelsuite}.  This model must have been developed 
#' in a modeling exercise testing only the impacts of time \emph{t}.
#' @param fec_model A linear model predicting fecundity. This can be a model
#' of class \code{glm} or \code{glmer}, and requires a predicted poisson variable under a 
#' log link. If given, then will overwrite any fecundity model given in 
#' \code{modelsuite}. This model must have been developed in a modeling exercise 
#' testing only the impacts of time \emph{t}.
#' @param jsurv_model A linear model predicting juvenile survival probability.
#' This can be a model of class \code{glm} or \code{glmer}, and requires a predicted 
#' binomial variable under a logit link. If given, then will overwrite any 
#' juvenile survival probability model given in \code{modelsuite}. This model must 
#' have been developed in a modeling exercise testing only the impacts of time \emph{t}.
#' @param jobs_model A linear model predicting juvenile sprouting or observation
#' probability. This can be a model of class \code{glm} or \code{glmer}, and requires a 
#' predicted binomial variable under a logit link. If given, then 
#' will overwrite any juvenile observation probability model given in \code{modelsuite}.
#' This model must have been developed in a modeling exercise testing only the 
#' impacts of time \emph{t}.
#' @param jsize_model A linear model predicting juvenile size. This can be a model
#' of class \code{glm} or \code{glmer}, both of which require a predicted poisson variable 
#' under a log link, or a model of class \code{lm} or \code{lmer}, in which a Gaussian 
#' response is assumed. If given, then will overwrite any juvenile size model 
#' given in \code{modelsuite}. This model must have been developed in a modeling 
#' exercise testing only the impacts of time \emph{t}.
#' @param jrepst_model A linear model predicting reproduction probability of a 
#' mature individual that was immature in the previous year. This can be a model 
#' of class \code{glm} or \code{glmer}, and requires a predicted binomial variable under a 
#' logit link. If given, then will overwrite any reproduction probability model 
#' given in \code{modelsuite}. This model must have been developed in a modeling 
#' exercise testing only the impacts of time \emph{t}.
#' @param paramnames A dataframe with two columns, the first showing the general
#' model terms that will be used in matrix creation, and the second showing the
#' equivalent terms used in modeling. Only required if \code{modelsuite} is not 
#' supplied.
#' @param surv_dev A numeric value to be added to the y-intercept in the linear
#' model for survival probability.
#' @param obs_dev A numeric value to be added to the y-intercept in the linear
#' model for observation probability.
#' @param size_dev A numeric value to be added to the y-intercept in the linear
#' model for size.
#' @param repst_dev A numeric value to be added to the y-intercept in the linear
#' model for probability of reproduction.
#' @param fec_dev A numeric value to be added to the y-intercept in the linear
#' model for fecundity.
#' @param jsurv_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile survival probability.
#' @param jobs_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile observation probability.
#' @param jsize_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile size.
#' @param jrepst_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile reproduction probability.
#' @param repmod A scalar multiplier of fecundity. Defaults to 1.
#' @param overwrite A data frame developed with the `overwrite` function,
#' describing transitions to be overwritten either with given values or with 
#' other estimated transitions.
#' @param data The original historical demographic data frame used to 
#' estimate vital rates (class \code{hfvdata}). The original data frame is 
#' required in order to initialize years and patches properly.
#' @param yearcol The variable name or column number corresponding to year 
#' in time \emph{t} in the dataset. Not needed if a \code{modelsuite} is supplied.
#' @param patchcol The variable name or column number corresponding to patch in 
#' the dataset. Not needed if a \code{modelsuite} is supplied.
#' @param year.as.random A logical term indicating whether coefficients for
#' missing patches within vital rate models should be estimated as random
#' intercepts. Defaults to FALSE, in which case missing time step coefficients 
#' are set to 0.
#' @param patch.as.random A logical term indicating whether coefficients for
#' missing patches within vital rate models should be estimated as random
#' intercepts. Defaults to FALSE, in which case missing patch coefficients 
#' are set to 0.
#' @param randomseed A numeric value used as a seed to generate random 
#' estimates for missing time step and patch coefficients, if either
#' `year.as.random` or `patch.as.random` is set to TRUE.
#' @param negfec A logical value denoting whether fecundity values estimated
#' to be negative should be reset to 0. Defaults to FALSE.
#' @param reduce A logical value denoting whether to remove ahistorical stages
#' associated solely with 0 transitions. These are only removed in cases where
#' the associated row and column sums in ALL matrices estimated equal 0. 
#' Defaults to FALSE.
#'
#' @return If all inputs are properly formatted, then this function will return
#' either an object of class \code{lefkoMat}. Output includes:
#'
#' \item{A}{A list of full projection matrices in order of sorted patches
#' and years.}
#' \item{U}{A list of survival-transition matrices sorted as in \code{A}.}
#' \item{F}{A list of fecundity matrices sorted as in \code{A}.}
#' \item{hstages}{Null for ahistorical matrices.}
#' \item{ahstages}{A data frame detailing the characteristics of associated
#' ahistorical stages.}
#' \item{labels}{A data frame giving the patch and year of each matrix in
#' order. In \code{flefko2()}, only one population may be analyzed at once, and 
#' so `pop = NA`.}
#' \item{matrixqc}{A short vector describing the number of non-zero elements
#' in \code{U} and \code{F} matrices, and the number of annual matrices.}
#' \item{modelqc}{This is the \code{qc} portion of the modelsuite input.}
#' 
#' Please note that this function will yield incorrect estimates if the models
#' utilized incorporate state in time \emph{t}-1. Only use models developed testing
#' ahistorical effects.
#'
#' @examples
#' \donttest{
#' data(lathyrus)
#' 
#' sizevector <- c(0, 4.6, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9)
#' stagevector <- c("Sd", "Sdl", "Dorm", "Sz1nr", "Sz2nr", "Sz3nr", "Sz4nr", "Sz5nr",
#'                  "Sz6nr", "Sz7nr", "Sz8nr", "Sz9nr", "Sz1r", "Sz2r", "Sz3r", "Sz4r",
#'                  "Sz5r", "Sz6r", "Sz7r", "Sz8r", "Sz9r")
#' repvector <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 4.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
#'             0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
#' 
#' 
#' lathframeln <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector,
#'                          obsstatus = obsvector, matstatus = matvector, immstatus = immvector,
#'                          indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvertln <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT",
#'                            individcol = "GENET", blocksize = 9, juvcol = "Seedling1988",
#'                            size1col = "lnVol88", repstr1col = "FCODE88",
#'                            fec1col = "Intactseed88", dead1col = "Dead1988",
#'                            nonobs1col = "Dormant1988", stageassign = lathframeln,
#'                            stagesize = "sizea", censorcol = "Missing1988",
#'                            censorkeep = NA, NAas0 = TRUE, censor = TRUE)
#' 
#' 
#' lathvertln$feca2 <- round(lathvertln$feca2)
#' lathvertln$feca1 <- round(lathvertln$feca1)
#' lathvertln$feca3 <- round(lathvertln$feca3)
#' 
#' lathrepmln <- matrix(0, 21, 21)
#' lathrepmln[1, c(13:21)] <- 0.345
#' lathrepmln[2, c(13:21)] <- 0.054
#' 
#' lathover2 <- overwrite(stage3 = c("Sd", "Sdl"), stage2 = c("Sd", "Sd"),
#'                        givenrate = c(0.345, 0.054))
#' 
#' 
#' lathmodelsln2 <- modelsearch(lathvertln, historical = FALSE, approach = "lme4", suite = "main",
#'                              vitalrates = c("surv", "obs", "size", "repst", "fec"), 
#'                              juvestimate = "Sdl", bestfit = "AICc&k", sizedist = "gaussian", 
#'                              fecdist = "poisson", indiv = "individ", patch = "patchid", 
#'                              year = "year2", year.as.random = TRUE, patch.as.random = TRUE,
#'                              show.model.tables = TRUE)
#'                              
#' lathmat2ln <- flefko2(year = "all", patch = "all", stageframe = lathframeln, 
#'                       modelsuite = lathmodelsln2, data = lathvertln, 
#'                       repmatrix = lathrepmln, overwrite = lathover2,
#'                       patchcol = "patchid", yearcol = "year2",
#'                       year.as.random = FALSE, patch.as.random = FALSE, 
#'                       reduce = FALSE)
#' 
#' summary(lathmat2ln)
#' }
#' 
#' @export
flefko2 <- function(year = "all", patch = "all", stageframe, repmatrix = NA, data = NA, 
                    modelsuite = NA, surv_model = NA, obs_model = NA, size_model = NA, 
                    repst_model = NA, fec_model = NA, jsurv_model = NA, jobs_model = NA, 
                    jsize_model = NA, jrepst_model = NA, paramnames = NA, surv_dev = 0, 
                    obs_dev = 0, size_dev = 0, repst_dev = 0, fec_dev = 0, jsurv_dev = 0, 
                    jobs_dev = 0, jsize_dev = 0, jrepst_dev = 0, repmod = 1, overwrite = NA, 
                    yearcol = NA, patchcol = NA, year.as.random = FALSE, patch.as.random = FALSE, 
                    randomseed = 0, negfec = FALSE, reduce = FALSE) {
  
  if (all(is.na(modelsuite)) & all(is.na(paramnames))) {
    warning("Function may not work properly without a dataframe of model parameters or equivalents supplied either through the modelsuite option or through the paramnames input parameter.")
  } else if (!all(is.na(modelsuite))) {
    paramnames <- modelsuite$paramnames
    yearcol <- paramnames$modelparams[which(paramnames$mainparams == "year2")]
    patchcol <- paramnames$modelparams[which(paramnames$mainparams == "patch")]
  }
  
  if (all(is.na(data))) {stop("Need original vertical dataset to set proper limits on year and patch.", call. = FALSE)}
  if (!any(class(data) == "data.frame")) {stop("Need original vertical dataset used in modeling to proceed.", call. = FALSE)}
  
  if (is.character(yearcol)) {
    choicevar <- which(names(data) == yearcol);
    mainyears <- sort(unique(data[,choicevar]))
  } else if (is.numeric(yearcol)) {
    mainyears <- sort(unique(data[, yearcol]));
  } else {
    stop("Need appropriate year column designation.", call. = FALSE)
  }
  
  if (any(is.character(year))) {
    if (is.element(tolower(year), "all")) {
      year <- mainyears
    } else {
      stop("Year designation not recognized.", call. = FALSE)
    }
  }
  
  if (length(year) == 0 | all(is.na(year) == TRUE) | any(is.na(year))) {
    stop("This function cannot proceed without being given a specific year, or a suite of years.", call. = FALSE)
  }
  
  if (is.na(patch) & !is.na(patchcol)) {
    warning("Matrix creation may not proceed properly without input in the patch option when using a modelsuite in which patch is designated.")
  }
  
  if (is.character(patchcol)) {
    choicevar <- which(names(data) == patchcol);
    mainpatches <- sort(unique(as.character(data[,choicevar])))
  } else if (is.numeric(patchcol)) {
    mainpatches <- sort(unique(as.character(data[, patchcol])));
  } else {
    mainpatches <- NA
  }
  
  if (any(is.character(patch))) {
    if (is.element(tolower(patch), "all")) {
      patch <- mainpatches
    } else if (!is.element(patch, mainpatches)) {
      stop("Patch designation not recognized.", call. = FALSE)
    }
  }
  
  if (all(is.na(repmatrix))) {    #This bit needs to be redone to check for fecundity values in the dataset and use them to determine the proper reproductive stages
    repmatrix <- matrix(0, dim(stageframe)[1], dim(stageframe)[1])
    repstages <- which(stageframe$repstatus == 1)
    repmatrix[1, repstages] <- 1
  }
  
  if (any(!suppressWarnings(!is.na(as.numeric(as.character(stageframe$size)))))) {
    stop("Function flefko2() requires size to be numeric rather than categorical.")
  }
  
  melchett <- .sf_reassess(stageframe, repmatrix, overwrite)
  stageframe <- melchett[[1]]
  repmatrix <- melchett[[2]]
  
  stageframe$new_stage_id <- as.numeric(stageframe$new_stage_id)
  stageframe$original_size <- as.numeric(stageframe$original_size)
  stageframe$bin_size_min <- as.numeric(stageframe$bin_size_min)
  stageframe$bin_size_ctr <- as.numeric(stageframe$bin_size_ctr)
  stageframe$bin_size_max <- as.numeric(stageframe$bin_size_max)
  stageframe$bin_size_width <- as.numeric(stageframe$bin_size_width)
  stageframe$bin_raw_halfwidth <- as.numeric(stageframe$bin_raw_halfwidth)
  stageframe$repstatus <- as.numeric(stageframe$repstatus)
  stageframe$obsstatus <- as.numeric(stageframe$obsstatus)
  stageframe$propstatus <- as.numeric(stageframe$propstatus)
  stageframe$immstatus <- as.numeric(stageframe$immstatus)
  stageframe$matstatus <- as.numeric(stageframe$matstatus)
  stageframe$indataset <- as.numeric(stageframe$indataset)
  stageframe$alive <- as.numeric(stageframe$alive)
  
  stageframe$stageno <- c(1:dim(stageframe)[1])
  rownames(stageframe) <- stageframe$stageno
  stageframe$fullstage <- apply(as.matrix(c(1:dim(stageframe)[1])), 1, function(X) {
    paste(stageframe$bin_size_ctr[X], stageframe$repstatus[X])
  })
  
  instages <- length(stageframe$new_stage_id)
  
  #This next portion creates a main reference table for the C++-based matrix populator functions, which are made to work on both flefko2() and flefko3()
  allstages.size <- expand.grid(size3 = stageframe$bin_size_ctr, size2n = stageframe$bin_size_ctr)
  allstages.size <- cbind.data.frame(allstages.size, allstages.size[,2], 0)
  allstages.obs <- expand.grid(obs3 = stageframe$obsstatus, obs2n = stageframe$obsstatus)
  allstages.obs <- cbind.data.frame(allstages.obs, allstages.obs[,2], 0)
  allstages.rep <- expand.grid(rep3 = stageframe$repstatus, rep2n = stageframe$repstatus)
  allstages.rep <- cbind.data.frame(allstages.rep, allstages.rep[,2], 0)
  allstages.mat <- expand.grid(mat3 = stageframe$matstatus, mat2n = stageframe$matstatus)
  allstages.mat <- cbind.data.frame(allstages.mat, allstages.mat[,2], allstages.mat[,2])
  allstages.imm <- expand.grid(imm3 = stageframe$immstatus, imm2n = stageframe$immstatus)
  allstages.imm <- cbind.data.frame(allstages.imm, allstages.imm[,2], allstages.imm[,2])
  allstages.re3 <- cbind(c(matrix(rbind(cbind(repmatrix, 0), 0)))) #Here we take repmatrix, add a 0 row and column for death, and vectorize
  allstages.ind <- expand.grid(indata3 = stageframe$indataset, indata2n = stageframe$indataset)
  allstages.ind <- cbind.data.frame(allstages.ind, allstages.ind[,2], 1)
  allstages.stages <- expand.grid(stage3 = stageframe$new_stage_id, stage2n = stageframe$new_stage_id)
  allstages.stages <- cbind.data.frame(allstages.stages, allstages.stages[,2], 0)
  allstages.bins <- rep(stageframe$bin_size_width, instages)
  allstages <- cbind.data.frame(allstages.stages, allstages.size, allstages.obs, allstages.rep, 
                                allstages.mat, allstages.imm, allstages.re3, allstages.ind, 
                                allstages.bins)
  names(allstages) <- c("stage3", "stage2n", "stage2o", "stage1", "size3", "size2n", "size2o", "size1", 
                        "obs3", "obs2n", "obs2o", "obs1", "rep3", "rep2n", "rep2o", "rep1", 
                        "mat3", "mat2n", "mat2o", "mat1", "imm3", "imm2n", "imm2o", "imm1",
                        "repentry3", "indata3", "indata2n", "indata2o", "indata1", "binwidth")
  
  allstages$minage3 <- 0
  allstages$minage2 <- 0
  allstages$maxage3 <- 0
  allstages$maxage2 <- 0
  allstages$actualage2 <- 0
  
  allstages$index32 <- allstages$stage3 + (allstages$stage2n * instages)
  
  if (!all(is.na(overwrite))) {
    
    overwrite <- .overwrite_reassess(overwrite, stageframe, historical = FALSE)
    if (dim(overwrite)[1] == 0) overwrite <- NA
  }
  
  if (!all(is.na(overwrite))) {
    
    overwrite$index3 <- apply(as.matrix(c(1:dim(overwrite)[1])), 1, function(X) {
      stageframe$stageno[which(stageframe$orig_stage_id == overwrite$stage3[X])]})
    overwrite$index2 <- apply(as.matrix(c(1:dim(overwrite)[1])), 1, function(X) {
      stageframe$stageno[which(stageframe$orig_stage_id == overwrite$stage2[X])]})
    overwrite$new3 <- apply(as.matrix(c(1:dim(overwrite)[1])), 1, function(X) {
      if(!is.na(overwrite$eststage3[X])) {
        stageframe$stageno[which(stageframe$orig_stage_id == overwrite$eststage3[X])]
      } else {overwrite$index3[X]}
    })
    overwrite$new2 <- apply(as.matrix(c(1:dim(overwrite)[1])), 1, function(X) {
      if(!is.na(overwrite$eststage2[X])) {
        stageframe$stageno[which(stageframe$orig_stage_id == overwrite$eststage2[X])]
      } else {overwrite$index2[X]}
    })
    overwrite$indexold32 <- overwrite$index3 + (overwrite$index2 * instages)
    overwrite$indexnew32 <- overwrite$new3 + (overwrite$new2 * instages)
    overwrite$givenrate[which(is.na(overwrite$givenrate))] <- 0
    
    overwrite$new3[which(is.na(overwrite$eststage3))] <- -1
    overwrite$givenrate[which(is.na(overwrite$givenrate))] <- -1
    overwrite$givenrate[which(overwrite$givenrate == 0)] <- -1
    
    allstageadditions <- ovreplace(allstages$index32, overwrite$indexold32, overwrite$indexnew32, 
                                   overwrite$convtype, overwrite$new3, overwrite$givenrate)
    
    allstages$ovest_t <- allstageadditions[,1]
    allstages$ovgiven_t <- allstageadditions[,2]
    allstages$ovest_f <- allstageadditions[,3]
    allstages$ovgiven_f <- allstageadditions[,4]
  } else {
    
    allstages$ovgiven_t <- -1
    allstages$ovest_t <- -1
    allstages$ovgiven_f <- -1
    allstages$ovest_f <- -1
  }
  
  if (class(modelsuite) == "lefkoMod") {
    if(is.na(surv_model)) {surv_model <- modelsuite$survival_model}
    if(is.na(obs_model)) {obs_model <- modelsuite$observation_model}
    if(is.na(size_model)) {size_model <- modelsuite$size_model}
    if(is.na(repst_model)) {repst_model <- modelsuite$repstatus_model}
    if(is.na(fec_model)) {fec_model <- modelsuite$fecundity_model}
    
    if(is.na(jsurv_model)) {jsurv_model <- modelsuite$juv_survival_model}
    if(is.na(jobs_model)) {jobs_model <- modelsuite$juv_observation_model}
    if(is.na(jsize_model)) {jsize_model <- modelsuite$juv_size_model}
    if(is.na(jrepst_model)) {jrepst_model <- modelsuite$juv_reproduction_model}
  }
  
  set.seed(randomseed)
  
  surv_proxy <- .modelextract(surv_model, paramnames, mainyears, mainpatches)
  obs_proxy <- .modelextract(obs_model, paramnames, mainyears, mainpatches)
  
  sigma <- 0
  rvarssummed <- 0
  sizedist <- 1
  
  size_proxy <- .modelextract(size_model, paramnames, mainyears, mainpatches)
  
  if (size_proxy$family == "poisson") {
    sizedist <- 0
    if (!all(is.na(size_proxy$variances))) {
      rvarssummed <- sum(size_proxy$variances[,"vcov"])
    } else {
      rvarssummed <- 0
    }
  } else if (size_proxy$family == "gaussian") {
    sizedist <- 2
    sigma <- size_proxy$sigma
  } else {
    sizedist <- 1
  }
  
  repst_proxy <- .modelextract(repst_model, paramnames, mainyears, mainpatches)
  fec_proxy <- .modelextract(fec_model, paramnames, mainyears, mainpatches)
  
  if (fec_proxy$family == "poisson") {
    fecdist <- 0
  } else if (fec_proxy$family == "gaussian") {
    fecdist <- 2
  } else {
    fecdist <- 1
  }
  
  jsurv_proxy <- .modelextract(jsurv_model, paramnames, mainyears, mainpatches)
  jobs_proxy <- .modelextract(jobs_model, paramnames, mainyears, mainpatches)
  
  jsigma <- 0
  jrvarssummed <- 0
  
  jsize_proxy <- .modelextract(jsize_model, paramnames, mainyears, mainpatches)
  
  if (jsize_proxy$family == "poisson") {
    if (!all(is.na(jsize_proxy$variances))) {
      jrvarssummed <- sum(jsize_proxy$variances[,"vcov"])
    } else {
      jrvarssummed <- 0
    }
  } else if (jsize_proxy$family == "gaussian") {
    jsigma <- jsize_proxy$sigma
  }
  
  jrepst_proxy <- .modelextract(jrepst_model, paramnames, mainyears, mainpatches)
  
  if (!all(is.na(patch))) {
    listofyears <- apply(as.matrix(patch), 1, function(X) {
      output <- cbind.data.frame(NA, X, as.matrix(year));
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    listofyears <- do.call(rbind.data.frame, listofyears)
    listofyears$poporder <- NA
    listofyears$patchorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainpatches == listofyears$patch[X])})
    listofyears$yearorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainyears == listofyears$year2[X])})
    
  } else {
    listofyears <- cbind.data.frame(NA, NA, as.matrix(year))
    names(listofyears) <- c("pop", "patch", "year2")
    
    listofyears$poporder <- NA
    listofyears$patchorder <- 1
    listofyears$yearorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainyears == listofyears$year2[X])})
  }
  
  maxsize <- max(c(allstages$size3, allstages$size2n, allstages$size2o), na.rm = TRUE)
  
  allstages$indata <- allstages$indata3 * allstages$indata2n * allstages$indata2o
  aliverows <- hoffmannofstuttgart(0, as.matrix(allstages))
  allstages <- allstages[(aliverows + 1),]
  
  allstages$aliveandequal <- c(0:(length(aliverows) - 1))
  allstages$r_aliveandequal <- allstages$aliveandequal + 1
  
  total.matrix.dim <- length(stageframe$bin_size_ctr) - 1
  patch.elements <- total.matrix.dim^2
  
  yearlist <- split(listofyears, seq(nrow(listofyears)))
  
  madsexmadrigal <- lapply(yearlist, jerzeibalowski, allstages, surv_proxy, obs_proxy, size_proxy, repst_proxy,
                           fec_proxy, jsurv_proxy, jobs_proxy, jsize_proxy, jrepst_proxy, surv_dev, obs_dev, 
                           size_dev, repst_dev, fec_dev, jsurv_dev, jobs_dev, jsize_dev, jrepst_dev, patch.elements,
                           total.matrix.dim, repmod, rvarssummed, sigma, jrvarssummed, jsigma, maxsize, sizedist, 
                           fecdist, negfec)
  
  a_list <- lapply(madsexmadrigal, function(X) {X$A})
  u_list <- lapply(madsexmadrigal, function(X) {X$U})
  f_list <- lapply(madsexmadrigal, function(X) {X$F})
  
  ahstages <- stageframe[1:(dim(stageframe)[1] - 1),]
  
  pairings1 <- expand.grid(stcod3 = stageframe$new_stage_id, stcod2 = stageframe$new_stage_id)
  pairings2 <- expand.grid(stage3 = stageframe$orig_stage_id, stage2 = stageframe$orig_stage_id)
  hstages <- cbind.data.frame(pairings2, pairings1)
  
  qcoutput1 <- NA
  qcoutput2 <- NA
  
  totalutransitions <- sum(unlist(lapply(u_list, function(X) {length(which(X != 0))})))
  totalftransitions <- sum(unlist(lapply(f_list, function(X) {length(which(X != 0))})))
  totalmatrices <- length(u_list)
  
  qcoutput1 <- c(totalutransitions, totalftransitions, totalmatrices)
  
  if (is.element("qc", names(modelsuite))) {qcoutput2 <- modelsuite$qc}
  
  if (reduce == TRUE) {
    drops <- .reducer2(a_list, u_list, f_list, ahstages)
    
    a_list <- drops$A
    u_list <- drops$U
    f_list <- drops$F
    ahstages <- drops$ahstages
  }
  
  output <- list(A = a_list, U = u_list, F = f_list, hstages = NA, ahstages = ahstages, 
                 labels = listofyears[,c(1:3)], matrixqc = qcoutput1, modelqc = qcoutput2)
  class(output) <- "lefkoMat"
  
  return(output)
}

#' Reduce Matrix Dimensions By Eliminating Empty Stages
#' 
#' \code{.reducer2()} identifies empty stages in a set of ahistorical matrices and
#' removes themfrom all matrices. It is used within \code{\link{flefko2}()} and 
#' \code{\link{rlefko2}()}.
#' 
#' @param A List of population projection matrices, from a \code{lefkoMat} object.
#' @param U List of surviva-transition matrices corresponding to \code{A}.
#' @param F List of fecundity matrices corresponding to \code{A}.
#' @param ahstages Data frame giving the names and identities of ahistorical 
#' stages used to create matrices.
#' 
#' @return Returns a list of reduced \code{A}, \code{U}, and \code{F} matrices, plus the reduced
#' \code{ahstages} object.
#' 
#' @keywords internal
#' @noRd
.reducer2 <- function(A, U, F, ahstages) {
  stagepatterns <- lapply(A, function(X) {
    matrix.sums <- colSums(X) + rowSums(X)
    return(matrix.sums)
  })
  
  used.stages.mat <- do.call("rbind", stagepatterns)
  used.stages.ovr <- colSums(used.stages.mat)
  keep.stages <- which(used.stages.ovr > 0)
  
  Ared <- lapply(A, function(X) {
    return(X[keep.stages, keep.stages])
  })
  
  Ured <- lapply(U, function(X) {
    return(X[keep.stages, keep.stages])
  })
  
  Fred <- lapply(F, function(X) {
    return(X[keep.stages, keep.stages])
  })
  
  ahstred <- ahstages[keep.stages,]
  
  return(list(A = Ared, U = Ured, F = Fred, ahstages = ahstred))
}

#' Create Raw Historical Population Projection Matrices
#' 
#' \code{rlefko3()} returns a list of raw historical population projection 
#' matrices, as well as the associated component transition and fecundity
#' matrices, data frames describing the ahistorical stages used and the pairing
#' of ahistorical stages used to create historical paired stages, and a data
#' frame describing the population, patch, and year associated with each matrix.
#' 
#' @param data A vertical demographic data frame, with variables corresponding 
#' to the naming conventions in \code{\link{verticalize3}()}.
#' @param stageframe A stageframe object that includes information on the size,
#' observation status, propagule status, immaturity status, and maturity status
#' of each ahistorical stage. Should also incorporate bin widths if size is
#' continuous.
#' @param year A variable corresponding to year or observation time, or a set
#' of such values, given in values associated with the year term used in linear 
#' model development. Can also equal \code{all}, in which case matrices will
#' be estimated for all years. Defaults to \code{all}.
#' @param pop A variable designating which populations will have matrices 
#' estimated. Should be set to specific population names, or to \code{all} if all 
#' populations should have matrices estimated.
#' @param patch A variable designating which patches or subpopulations will
#' have matrices estimated. Should be set to specific patch names, or to 
#' \code{all} if all patches should have matrices estimated.
#' @param censor If TRUE, then data will be removed according to the variable
#' set in \code{censorcol}, such that only data with censor values equal to 1
#' will remain. Defaults to FALSE.
#' @param stages An optional but important vector denoting the names of the
#' variables within the main vertical dataset coding for the names of the stages
#' of each individual in times \emph{t}+1, \emph{t}, and \emph{t}-1, respectively. The names of
#' the stages in these variables should match those used in the \code{stageframe}
#' input for this analysis exactly. If left blank, then \code{rlefko3()} will attempt
#' to infer stages by matching values of \code{alive}, \code{size}, \code{repst}, and \code{matst}
#' to characteristics noted in the associated \code{stageframe}.
#' @param alive A vector of names of binomial variables corresponding to status 
#' as alive (1) or dead (0) in times \emph{t}+1, \emph{t}, and \emph{t}-1, respectively.
#' @param size A vector of names of variables coding size in times \emph{t}+1, \emph{t},
#' and \emph{t}-1, respectively. Defaults to \code{c("sizea3", "sizea2", "sizea1")}.
#' @param repst A vector of names of variables coding reproductive status in 
#' times \emph{t}+1, \emph{t}, and \emph{t}-1, respectively. Defaults to \code{c("repstatus3", 
#' "repstatus2", "repstatus1")}.
#' @param matst A vector of names of variables coding maturity status in 
#' times \emph{t}+1, \emph{t}, and \emph{t}-1, respectively. Defaults to \code{c("matstatus3", 
#' "matstatus2", "matstatus1")}. Must be supplied if \code{stages} is not provided.
#' @param fec A vector of names of variables coding fecundity in times \emph{t}+1,
#' \emph{t}, and \emph{t}-1, respectively. Defaults to \code{c("feca3", "feca2", "feca1")}.
#' @param repmatrix A matrix composed mostly of 0s, with non-zero values for
#' each potentially new individual (row) born to each reproductive stage
#' (column). Non-zero entries correspond to multipliers for fecundity, with 1
#' equaling full fecundity.
#' @param overwrite A data frame developed with the \code{\link{overwrite}()} function,
#' describing transitions to be overwritten either with given values or with 
#' other estimated transitions.
#' @param yearcol The variable name or column number corresponding to year 
#' in time \emph{t} in the dataset.
#' @param popcol The variable name or column number corresponding to the
#' identity of the population.
#' @param patchcol The variable name or column number corresponding to patch in 
#' the dataset.
#' @param indivcol The variable name or column number coding individual 
#' identity.
#' @param censorcol The variable name or column number denoting the
#' censor status. Only needed if \code{censor = TRUE}.
#' @param reduce A logical value denoting whether to remove historical stages
#' associated with only zero transitions. These are removed only if all row
#' and column sums in ALL matrices estimated equal 0. Defaults to FALSE.
#'
#' @return If all inputs are properly formatted, then this function will return
#' either an object of class \code{lefkoMat}. Output includes:
#'
#' \item{A}{A list of full projection matrices in order of sorted populations,
#' patches, and years.}
#' \item{U}{A list of survival-transition matrices sorted as in \code{A}.}
#' \item{F}{A list of fecundity matrices sorted as in \code{A}.}
#' \item{hstages}{A data frame matrix showing the pairing of ahistorical 
#' stages used to create historical stage pairs.}
#' \item{ahstages}{A data frame detailing the characteristics of associated
#' ahistorical stages.}
#' \item{labels}{A data frame giving the population, patch, and year of each 
#' matrix in order.}
#' \item{matrixqc}{A short vector describing the number of non-zero elements
#' in \code{U} and \code{F} matrices, and the number of annual matrices.}
#' \item{dataqc}{A vector showing the numbers of individuals and rows in the
#' vertical dataset used as input.}
#'
#' @examples
#' data(cypdata)
#' 
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", 
#'                  "Sm", "Md", "Lg", "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
#'                           repstatus = repvector, obsstatus = obsvector,
#'                           matstatus = matvector, propstatus = propvector,
#'                           immstatus = immvector, indataset = indataset,
#'                           binhalfwidth = binvec)
#' 
#' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
#'                           patchidcol = "patch", individcol = "plantid",
#'                           blocksize = 4, size1col = "Inf2.04", size2col = "Inf.04",
#'                           size3col = "Veg.04", repstr1col = "Inf.04",
#'                           repstr2col = "Inf2.04", fec1col = "Pod.04",
#'                           stageassign = cypframe_raw, stagesize = "sizeadded",
#'                           NAas0 = TRUE, NRasRep = TRUE)
#' 
#' rep_cyp_raw <- matrix(0, 11, 11)
#' rep_cyp_raw[1:2,7:11] <- 0.5
#' 
#' cypover3r <- overwrite(stage3 = c("SD", "SD", "P1", "P1", "P2", "P3", "D", "XSm", "Sm",
#'                        "SL", "SL", "SL"), stage2 = c("SD", "SD", "SD", "SD", "P1", "P2",
#'                        "P3", "P3", "P3", "P3", "SL", "SL"), stage1 = c("SD", "rep", "SD",
#'                        "rep", "SD", "P1", "P2", "P2", "P2", "P2", "P3", "SL"),
#'                        eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA, NA),
#'                        eststage2 = c(NA, NA, NA, NA, NA, NA, "D", "D", "D", NA, NA, NA),
#'                        eststage1 = c(NA, NA, NA, NA, NA, NA, "D", "D", "D", NA, NA, NA),
#'                        givenrate = c(0.1, 0.1, 0.2, 0.2, 0.2, 0.2, NA, NA, NA, 0.25, 0.4, 0.4),
#'                        type = c("S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S"))
#' 
#' cypmatrix3r <- rlefko3(data = cypraw_v1, stageframe = cypframe_raw, year = "all",
#'                        patch = "all", stages = c("stage3", "stage2", "stage1"),
#'                        size = c("size3added", "size2added", "size1added"),
#'                        repmatrix = rep_cyp_raw, overwrite = cypover3r,
#'                        yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' 
#' summary(cypmatrix3r)
#' 
#' @export
rlefko3 <- function(data, stageframe, year = "all", pop = NA, patch = NA, censor = FALSE, 
                    stages = NA, alive = c("alive3", "alive2", "alive1"), 
                    size = c("sizea3", "sizea2", "sizea1"), repst = c("repstatus3", "repstatus2", "repstatus1"), 
                    matst = c("matstatus3", "matstatus2", "matstatus1"), fec = c("feca3", "feca2", "feca1"), 
                    repmatrix = NA, overwrite = NA, yearcol = NA, popcol = NA, patchcol = NA, indivcol = NA, 
                    censorcol = NA, reduce = FALSE) {
  
  tocensor <- indataset <- alive2 <- popused <- patchused <- yearused <- NULL
  
  if (all(is.na(data))) {
    stop("Need original vertical dataset to proceed.")
  }
  
  if (all(is.na(stages))) {
    if ((length(alive) != 3)) {
      stop("This function requires stage informationfor each of times t+1, t, and t-1. In the absence of stage columns in the dataset, it requires three variables for living/dead status, size, reproductive status, and maturity status, for each of times t+1, t, and t-1.")
    }
    if ((length(size) != 3)) {
      stop("This function requires stage informationfor each of times t+1, t, and t-1. In the absence of stage columns in the dataset, it requires three variables for living/dead status, size, reproductive status, and maturity status, for each of times t+1, t, and t-1.")
    }
    if (!all(is.na(repst))) {
      if ((length(repst) != 3)) {
        stop("This function requires stage informationfor each of times t+1, t, and t-1. In the absence of stage columns in the dataset, it requires three variables for living/dead status, size, reproductive status, and maturity status, for each of times t+1, t, and t-1.")
      }
    }   
    if (!all(is.na(matst))) {
      if ((length(matst) != 3)) {
        stop("This function requires stage informationfor each of times t+1, t, and t-1. In the absence of stage columns in the dataset, it requires three variables for living/dead status, size, reproductive status, and maturity status, for each of times t+1, t, and t-1.")
      }
    }   
  } else if (length(stages) != 3) {
    stop("This function requires stage informationfor each of times t+1, t, and t-1.")
  }
  
  if ((length(fec) != 3)) {
    stop("This function requires two variables for fecundity, for each of times t+1, t, and t-1.")
  }
  
  if (any(is.character(year)) & any(class(data) == "hfvdata")) {
    if (is.element(tolower(year), "all")) {
      if (is.character(yearcol)) {
        choicevar <- which(names(data) == yearcol)
        year <- sort(unique(data[,choicevar]))[-1] #Here we remove the 1st year because it has no time t-1 in the dataset
        yearcol <- choicevar
      } else if (all(is.numeric(year))) {
        year <- sort(unique(data[,yearcol]))[-1]
      } else {stop("Cannot understand year designation.")}
    }
  }
  
  if (length(year) == 0 | all(is.na(year) == TRUE)) {
    stop("This function cannot proceed without being given a specific year, or a suite of years.")
  }
  
  if (!all(is.element(year, sort(unique(data[,yearcol]))))) {
    stop("Dataset does not contain one or more of the requested years.")
  }
  
  if (all(is.na(repmatrix))) {    #This bit needs to be redone to check for fecundity values in the dataset and use them to determine the proper reproductive stages
    repmatrix <- matrix(0, dim(stageframe)[1], dim(stageframe)[1])
    repstages <- which(stageframe$repstatus == 1)
    repmatrix[1, repstages] <- 1
  }
  
  if (censor == TRUE) {
    if(all(is.na(censorcol)) == TRUE) {
      stop("Cannot censor the data without a proper censor variable.")
    }
    
    for (i in c(1:length(censorcol))) {
      if (is.character(censorcol)) {
        data$tocensor <- data[,which(names(data) == censorcol[i])]
      } else {
        data$tocensor <- data[,censorcol[i]]
      }
      data <- subset(data, tocensor == 1)
    }   
  }
  
  if (!all(is.na(pop)) & !all(is.na(patch))) {
    if (is.na(popcol) | is.na(patchcol)) {stop("Need population and patch designation variables to proceed.")}
    
    if (is.element(tolower(pop), "all")) {
      if (is.character(popcol)) {popcol <- which(names(data) == popcol)}
      
      pops <- unique(data[,popcol])
    } else {pops <- pop}
    
    if (is.element(tolower(patch), "all")) {
      if (is.character(patchcol)) {patchcol <- which(names(data) == patchcol)}
      
      listofpatches <- apply(as.matrix(pops), 1, function(X) {
        patchfolly <- subset(data, popcol == X);
        output <- cbind.data.frame(X, sort(unique(patchfolly[,patchcol])));
        names(output) <- c("pop", "patch");
        return(output);
      })
      
      if (length(listofpatches) > 1) {
        listofpatches <- do.call(rbind.data.frame, listofpatches)
      }
    } else {listofpatches <- expand.grid(pop = pops, patch = patch)}
    
    listofyears <- apply(as.matrix(listofpatches), 1, function(X) {
      checkyrdata <- subset(data, popcol = X[1]);
      checkyrdata <- subset(checkyrdata, patchcol = X[2])
      output <- cbind.data.frame(X[1], X[2], sort(unique(checkyrdata[,yearcol]))[-1]);
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
  } else if (all(is.na(pop)) & !all(is.na(patch))) {
    if (is.na(patchcol)) {stop("Need patch designation variable to proceed.")}
    
    if (is.element(tolower(patch), "all")) {
      if (is.character(patchcol)) {patchcol <- which(names(data) == patchcol)}
      
      patches <- sort(unique(data[,patchcol]))
    } else {patches <- patch}
    
    listofyears <- apply(as.matrix(patches), 1, function(X) {
      checkyrdata <- subset(data, patchcol = X);
      output <- cbind.data.frame(NA, X, sort(unique(checkyrdata[,yearcol]))[-1]);
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    if (length(listofyears) > 1) {
      listofyears <- do.call(rbind.data.frame, listofyears)
    }
  } else if (!all(is.na(pop)) & all(is.na(patch))) {
    if (is.na(popcol)) {stop("Need population designation variable to proceed.")}
    
    if (is.element(tolower(pop), "all")) {
      if (is.character(popcol)) {popcol <- which(names(data) == popcol)}
      
      pops <- unique(data[,popcol])
    } else {pops <- pop}
    
    listofyears <- apply(as.matrix(pops), 1, function(X) {
      checkyrdata <- subset(data, popcol = X);
      output <- cbind.data.frame(X, NA, sort(unique(checkyrdata[,yearcol]))[-1]);
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    if (length(listofyears) > 1) {
      listofyears <- do.call(rbind.data.frame, listofyears)
    }
  } else if (all(is.na(pop)) & all(is.na(patch))) {
    listofyears <- cbind.data.frame(NA, NA, year)
    names(listofyears) <- c("pop", "patch", "year2")
  }
  
  melchett <- .sf_reassess(stageframe, repmatrix, overwrite)
  stageframe <- melchett[[1]]
  repmatrix <- melchett[[2]]
  
  stageframe$new_stage_id <- as.numeric(stageframe$new_stage_id)
  stageframe$original_size <- as.numeric(stageframe$original_size)
  stageframe$bin_size_min <- as.numeric(stageframe$bin_size_min)
  stageframe$bin_size_ctr <- as.numeric(stageframe$bin_size_ctr)
  stageframe$bin_size_max <- as.numeric(stageframe$bin_size_max)
  stageframe$bin_size_width <- as.numeric(stageframe$bin_size_width)
  stageframe$bin_raw_halfwidth <- as.numeric(stageframe$bin_raw_halfwidth)
  stageframe$repstatus <- as.numeric(stageframe$repstatus)
  stageframe$obsstatus <- as.numeric(stageframe$obsstatus)
  stageframe$propstatus <- as.numeric(stageframe$propstatus)
  stageframe$immstatus <- as.numeric(stageframe$immstatus)
  stageframe$matstatus <- as.numeric(stageframe$matstatus)
  stageframe$indataset <- as.numeric(stageframe$indataset)
  stageframe$alive <- as.numeric(stageframe$alive)
  
  stageframe$stageno <- c(1:dim(stageframe)[1])
  rownames(stageframe) <- stageframe$stageno
  
  data$alive1 <- data[,which(names(data) == alive[3])]
  data$alive2 <- data[,which(names(data) == alive[2])]
  data$alive3 <- data[,which(names(data) == alive[1])]
  
  instageframe <- subset(stageframe, indataset == 1)
  
  if (all(is.na(stages))) {
    if (length(size) > 1) {
      data$usedsize1 <- data[,which(names(data) == size[3])]
      data$usedsize2 <- data[,which(names(data) == size[2])]
      data$usedsize3 <- data[,which(names(data) == size[1])]
    } else {
      warning("Without stage columns, lefko3 MPM estimation functions generally require size variables. Failure to include size variables may lead to odd results.")
    }
    if (length(repst) > 1) {
      data$usedrepst1 <- data[,which(names(data) == repst[3])]
      data$usedrepst2 <- data[,which(names(data) == repst[2])]
      data$usedrepst3 <- data[,which(names(data) == repst[1])]
    } 
    if (length(matst) > 1) {
      data$usedmatstatus1 <- data[,which(names(data) == matst[3])]
      data$usedmatstatus2 <- data[,which(names(data) == matst[2])]
      data$usedmatstatus3 <- data[,which(names(data) == matst[1])]
    } 
    
    data$usedstage1 <- apply(as.matrix(c(1:dim(data)[1])), 1, function(X) {
      if (is.na(data$usedsize1[X])) {
        data$usedsize1[X] <- 0
      }
      mainstages <- intersect(which(instageframe$bin_size_min < data$usedsize1[X]), 
                              which(instageframe$bin_size_max >= data$usedsize1[X]))
      jmstages <- which(instageframe$immstatus == (1 - data$usedmatstatus1[X]))
      obsstages <- which(instageframe$obsstatus == data$obsstatus1[X])
      repstages <- which(instageframe$repstatus == data$repstatus1[X])
      alivestage1 <- which(instageframe$alive == data$alive1[X])
      
      choicestage <- intersect(intersect(intersect(mainstages, jmstages), intersect(obsstages, repstages)), alivestage1)
      
      if (length(choicestage) == 0) choicestage <- which(instageframe$new_stage_id == max(instageframe$new_stage_id))

      if (length(choicestage) == 0) {
        stop("Stage characteristics mismatch dataset. Consider using the stages option, particularly if the vertical file was created with NRasRep = TRUE in verticalize3() or historicalize3().")
      }
      
      return(instageframe$orig_stage_id[choicestage])
    })
    
    data$usedstage2 <- apply(as.matrix(c(1:dim(data)[1])), 1, function(X) {
      if (is.na(data$usedsize2[X])) {
        data$usedsize2[X] <- 0
      }
      mainstages <- intersect(which(instageframe$bin_size_min < data$usedsize2[X]), 
                              which(instageframe$bin_size_max >= data$usedsize2[X]))
      jmstages <- which(instageframe$immstatus == (1 - data$usedmatstatus2[X]))
      obsstages <- which(instageframe$obsstatus == data$obsstatus2[X])
      repstages <- which(instageframe$repstatus == data$repstatus2[X])
      
      choicestage <- intersect(intersect(mainstages, jmstages), intersect(obsstages, repstages))
      
      if (length(choicestage) == 0) {
        stop("Stage characteristics mismatch dataset. Consider using the stages option, particularly if the vertical file was created with NRasRep = TRUE in verticalize3() or historicalize3().")
      }
      
      return(instageframe$orig_stage_id[choicestage])
    })
    
    data$usedstage3 <- apply(as.matrix(c(1:dim(data)[1])), 1, function(X) {
      if (is.na(data$usedsize3[X])) {
        data$usedsize3[X] <- 0
      }
      mainstages <- intersect(which(instageframe$bin_size_min < data$usedsize3[X]), 
                              which(instageframe$bin_size_max >= data$usedsize3[X]))
      jmstages <- which(instageframe$immstatus == (1 - data$usedmatstatus3[X]))
      obsstages <- which(instageframe$obsstatus == data$obsstatus3[X])
      repstages <- which(instageframe$repstatus == data$repstatus3[X])
      alivestage3 <- which(instageframe$alive == data$alive3[X])
      
      choicestage <- intersect(intersect(intersect(mainstages, jmstages), intersect(obsstages, repstages)), alivestage3)
      
      if (length(choicestage) == 0) choicestage <- which(instageframe$new_stage_id == max(instageframe$new_stage_id))
      
      return(instageframe$orig_stage_id[choicestage])
    })
    
  } else if (length(stages) > 1) {
    if (is.numeric(stages[2])) {
      data$usedstage1 <- data[, stages[3]]
      data$usedstage1[which(data$usedstage1 == "NotAlive")] <- "Dead"
      data$usedstage2 <- data[, stages[2]]
      data$usedstage3 <- data[, stages[1]]
      data$usedstage3[which(data$usedstage3 == "NotAlive")] <- "Dead"
    } else {
      data$usedstage1 <- data[,which(names(data) == stages[3])]
      data$usedstage1[which(data$usedstage1 == "NotAlive")] <- "Dead"
      data$usedstage2 <- data[,which(names(data) == stages[2])]
      data$usedstage3 <- data[,which(names(data) == stages[1])]
      data$usedstage3[which(data$usedstage3 == "NotAlive")] <- "Dead"
    }
    stages.used <- sort(unique(c(data$usedstage2, data$usedstage3)))
    
    if (length(setdiff(stages.used, stageframe$orig_stage_id)) > 0 & !is.element("NotAlive", stages.used)) {
      stop("Some stages in dataset do not match those detailed in the input stageframe.", .call = FALSE)
    }
  }
  
  if (length(fec) > 1) {
    data$usedfec1 <- data[,which(names(data) == fec[3])]
    data$usedfec2 <- data[,which(names(data) == fec[2])]
    data$usedfec3 <- data[,which(names(data) == fec[1])]
    
    data$usedfec1[which(is.na(data$usedfec1))] <- 0
    data$usedfec2[which(is.na(data$usedfec2))] <- 0
    data$usedfec3[which(is.na(data$usedfec3))] <- 0
  } else {
    warning("Lefko3 MPM estimation functions generally require fecundity variables. Failure to include fecundity variables leads to matrices composed only of survival transitions.")
  } 
  
  stageexpansion3a <- cbind.data.frame(expand.grid(size2o = stageframe$bin_size_ctr, size1 = stageframe$bin_size_ctr), 
                                       expand.grid(rep2o = stageframe$repstatus, rep1 = stageframe$repstatus),
                                       expand.grid(indata2o = stageframe$indataset, indata1 = stageframe$indataset),
                                       expand.grid(stage2o = stageframe$stageno, stage1 = stageframe$stageno),
                                       fec2o1 = c(cbind(rbind(repmatrix, 0), 0)))
  
  stageexpansion3a$indata2o1 <- stageexpansion3a$indata2o * stageexpansion3a$indata1
  
  stageexpansion3b <- cbind.data.frame(expand.grid(size3 = stageframe$bin_size_ctr, size2n = stageframe$bin_size_ctr), 
                                       expand.grid(rep3 = stageframe$repstatus, rep2n = stageframe$repstatus),
                                       expand.grid(indata3 = stageframe$indataset, indata2n = stageframe$indataset),
                                       expand.grid(stage3 = stageframe$stageno, stage2n = stageframe$stageno),
                                       fec32n = c(cbind(rbind(repmatrix, 0), 0)))
  
  stageexpansion3b$indata32n <- stageexpansion3b$indata3 * stageexpansion3b$indata2n
  
  instages <- length(stageframe$new_stage_id)
  
  stageexpansion3b$pairindex <- apply(as.matrix(c(1:dim(stageexpansion3b)[1])), 1, function(X) {
    (stageexpansion3b$stage3[X] - 1) + ((stageexpansion3b$stage2n[X] - 1) * instages)
  })
  
  stageexpansion9 <- cbind.data.frame(expand.grid(index3 = stageexpansion3b$stage3, index2o = stageexpansion3a$stage2o), 
                                      expand.grid(index2n = stageexpansion3b$stage2n, index1 = stageexpansion3a$stage1), 
                                      expand.grid(size3 = stageexpansion3b$size3, size2o = stageexpansion3a$size2o),  
                                      expand.grid(size2n = stageexpansion3b$size2n, size1 = stageexpansion3a$size1), 
                                      expand.grid(rep3 = stageexpansion3b$rep3, rep2o = stageexpansion3a$rep2o), 
                                      expand.grid(rep2n = stageexpansion3b$rep2n, rep1 = stageexpansion3a$rep1), 
                                      expand.grid(indata32n = stageexpansion3b$indata32n, indata2o1 = 
                                                    stageexpansion3a$indata2o1), 
                                      expand.grid(fec32n = stageexpansion3b$fec32n, fec2o1 = stageexpansion3a$fec2o1))
  
  stageexpansion9$index3221 <- (stageexpansion9$index3 - 1) + ((stageexpansion9$index2n - 1) * instages) + 
    ((stageexpansion9$index2o - 1) * instages * instages) + ((stageexpansion9$index1 - 1) * instages * instages * instages)
  stageexpansion9$indata3221 <- stageexpansion9$indata32n * stageexpansion9$indata2o1
  stageexpansion9$indata3221[which(stageexpansion9$index2n / stageexpansion9$index2o != 1)] <- 0
  
  stageexpansion9$stagepair2n1 <- (stageexpansion9$index2n - 1) + ((stageexpansion9$index1 - 1) * instages)
  
  stageexpansion3 <- stageexpansion3b
  names(stageexpansion3) <- c("size3", "size2", "rep3", "rep2", "indata3", "indata2", "stage3", "stage2",
                              "fec32", "indata32", "pairindex")
  stageexpansion3$stcod3 <- apply(as.matrix(c(1:dim(stageexpansion3)[1])), 1, function(X) {
    stageframe$orig_stage_id[which(stageframe$stageno == stageexpansion3$stage3[X])]
  })
  stageexpansion3$stcod2 <- apply(as.matrix(c(1:dim(stageexpansion3)[1])), 1, function(X) {
    stageframe$orig_stage_id[which(stageframe$stageno == stageexpansion3$stage2[X])]
  })
  
  stageexpansion9$ovest_t <- 0
  stageexpansion9$ovgiven_t <- 0
  stageexpansion9$ovest_f <- 0
  stageexpansion9$ovgiven_f <- 0
  
  if (!all(is.na(overwrite))) {
    
    overwrite <- .overwrite_reassess(overwrite, stageframe, historical = TRUE)
    if (dim(overwrite)[1] == 0) overwrite <- NA
  }
  
  if (!all(is.na(overwrite))) {
    
    overwrite$index3 <- apply(as.matrix(c(1:dim(overwrite)[1])), 1, function(X) {
      stageframe$stageno[which(stageframe$orig_stage_id == overwrite$stage3[X])]})
    overwrite$index2 <- apply(as.matrix(c(1:dim(overwrite)[1])), 1, function(X) {
      stageframe$stageno[which(stageframe$orig_stage_id == overwrite$stage2[X])]})
    overwrite$index1 <- apply(as.matrix(c(1:dim(overwrite)[1])), 1, function(X) {
      stageframe$stageno[which(stageframe$orig_stage_id == overwrite$stage1[X])]})
    overwrite$new3 <- apply(as.matrix(c(1:dim(overwrite)[1])), 1, function(X) {
      if(!is.na(overwrite$eststage3[X])) {
        stageframe$stageno[which(stageframe$orig_stage_id == overwrite$eststage3[X])]
      } else {overwrite$index3[X]}
    })
    overwrite$new2 <- apply(as.matrix(c(1:dim(overwrite)[1])), 1, function(X) {
      if(!is.na(overwrite$eststage2[X])) {
        stageframe$stageno[which(stageframe$orig_stage_id == overwrite$eststage2[X])]
      } else {overwrite$index2[X]}
    })
    overwrite$new1 <- apply(as.matrix(c(1:dim(overwrite)[1])), 1, function(X) {
      if(!is.na(overwrite$eststage1[X])) {
        stageframe$stageno[which(stageframe$orig_stage_id == overwrite$eststage1[X])]
      } else {overwrite$index1[X]}
    })
    
    overwrite$indexold3221 <- (overwrite$index3 - 1) + ((overwrite$index2 - 1) * instages) + 
      ((overwrite$index2 - 1) * instages * instages) + ((overwrite$index1 - 1)  * instages * instages * instages)
    overwrite$indexnew3221 <- (overwrite$new3 - 1) + ((overwrite$new2 - 1) * instages) + 
      ((overwrite$new2 - 1) * instages * instages) + ((overwrite$new1 - 1)  * instages * instages * instages)
    overwrite$givenrate[which(is.na(overwrite$givenrate))] <- 0
    
    overwrite$new3[which(is.na(overwrite$eststage3))] <- -1
    overwrite$givenrate[which(is.na(overwrite$givenrate))] <- -1
    overwrite$givenrate[which(overwrite$givenrate == 0)] <- -1
    
    allstageadditions <- ovreplace(stageexpansion9$index3221, overwrite$indexold3221, overwrite$indexnew3221, 
                                   overwrite$convtype, overwrite$new3, overwrite$givenrate)
    
    stageexpansion9$ovest_t <- allstageadditions[,1]
    stageexpansion9$ovgiven_t <- allstageadditions[,2]
    stageexpansion9$ovest_f <- allstageadditions[,3]
    stageexpansion9$ovgiven_f <- allstageadditions[,4]
  }
  
  data <- subset(data, alive2 == 1)
  
  data$index1 <- apply(as.matrix(data$usedstage1), 1, function(X) {
    instageframe$stageno[which(instageframe$orig_stage_id == X)]
    
  })
  data$index2 <- apply(as.matrix(data$usedstage2), 1, function(X) {
    instageframe$stageno[which(instageframe$orig_stage_id == X)]
  })
  data$index3 <- apply(as.matrix(data$usedstage3), 1, function(X) {
    instageframe$stageno[which(instageframe$orig_stage_id == X)]
  })
  data$index3221 <- apply(as.matrix(c(1:length(data$usedstage1))), 1, function(X) {
    ((data$index3[X] - 1) + ((data$index2[X] - 1) * instages) + ((data$index2[X] - 1) * instages * instages) + 
       ((data$index1[X] - 1)  * instages * instages * instages))
  })
  data$pairindex21 <- apply(as.matrix(c(1:dim(data)[1])), 1, function(X) {
    (data$index2[X] - 1) + ((data$index1[X] - 1) * instages)})
  
  data$usedfec2[which(is.na(data$usedfec2))] <- 0
  
  if(is.element(0, unique(data$index1))) {
    warning("Data (stage at time t-1) contains some stages not identified in stageframe. Note that all stage characteristics must match, including reproductive status.")
  }
  if(is.element(0, unique(data$index2))) {
    warning("Data (stage at time t) contains some stages not identified in stageframe. Note that all stage characteristics must match, including reproductive status.")
  }
  if(is.element(0, unique(data$index3))) {
    warning("Data (stage at time t+1) contains some stages not identified in stageframe. Note that all stage characteristics must match, including reproductive status.")
  }
  
  aliverows <- hoffmannofstuttgart(0, as.matrix(stageexpansion9[, c("index3", "index2n", "index2o", "index1")]))
  
  madsexmadrigal <- apply(listofyears, 1, function(X) {
    passed_data <- data
    if (!is.na(X[1])) {
      passed_data$popused <- passed_data[,popcol];
      passed_data <- subset(passed_data, popused == X[1]);
    }
    if (!is.na(X[2])) {
      passed_data$patchused <- passed_data[,patchcol];
      passed_data <- subset(passed_data, patchused == X[2]);
    }
    if (!is.na(X[3])) {
      passed_data$yearused <- passed_data[,yearcol];
      passed_data <- subset(passed_data, yearused == X[3]);
    }
    
    .rlefko3_core(data = passed_data, stageframe = stageframe, stageexpansion9 = stageexpansion9, 
                  aliverows = aliverows, stageexpansion3 = stageexpansion3)
  })
  
  a_list <- lapply(madsexmadrigal, function(X) {X$A})
  u_list <- lapply(madsexmadrigal, function(X) {X$U})
  f_list <- lapply(madsexmadrigal, function(X) {X$F})
  
  qcoutput1 <- NA
  qcoutput2 <- NA
  
  totalutransitions <- sum(unlist(lapply(u_list, function(X) {length(which(X != 0))})))
  totalftransitions <- sum(unlist(lapply(f_list, function(X) {length(which(X != 0))})))
  totalmatrices <- length(u_list)
  
  qcoutput1 <- c(totalutransitions, totalftransitions, totalmatrices)
  
  indivs <- NA
  if (!all(is.na(indivcol))) {
    if (all(is.character(indivcol))) {indivcol <- which(names(data) == indivcol)[1]}
    indivs <- length(unique(data[,indivcol]))
  }
  qcoutput2 <- c(indivs, dim(data)[1])
  
  morebitstolose <- unique(c(which(stageexpansion3$stage3 == dim(stageframe)[1]), which(stageexpansion3$stage2 == dim(stageframe)[1])))
  stageexpansion3 <- stageexpansion3[-morebitstolose,]
  
  hstages <- stageexpansion3[,c("stcod3", "stcod2", "stage3", "stage2")]
  
  if (reduce == TRUE) {
    drops <- .reducer3(a_list, u_list, f_list, hstages)
    
    a_list <- drops$A
    u_list <- drops$U
    f_list <- drops$F
    hstages <- drops$hstages
  }
  
  output <- list(A = a_list, U = u_list, F = f_list, hstages = hstages, 
                 ahstages = stageframe[1:(dim(stageframe)[1] - 1),], labels = listofyears,
                 matrixqc = qcoutput1, dataqc = qcoutput2)
  class(output) <- "lefkoMat"
  
  return(output)
}

#' Core Wrapper Powering Historical Raw Matrix Estimation
#' 
#' \code{.rlefko3_core()} pulls together all required data and matrix conditions
#' and feeds them into\code{\link{specialpatrolgroup}()}, an Rcpp function designed
#' to quickly estimate historical raw matrices.
#' 
#' @param data Original demographic dataset showing the states and fates of
#' individuals.
#' @param stageframe Original stageframe used in analysis.
#' @param stageexpansion9 Massive data frame created by \code{\link{flefko3}()} including
#' general life history and statistical characteristics identifying every element
#' in the projection matrix.
#' @param aliverows Vector identifying elements in the projection matrix that are
#' both alive and logically possible.
#' @param stageexpansion3 Data frame characterizing all stage pairs used in
#' historical matrix estimation.
#' 
#' @return Returns a list containing one set of matrices corresponding to a 
#' specific population, patch, and time step combination.
#' 
#' @keywords internal
#' @noRd
.rlefko3_core <- function(data, stageframe, stageexpansion9 = NA, aliverows = NA, 
                          stageexpansion3 = NA) {
  
  if (all(is.na(stageexpansion9)) | all(is.na(stageexpansion3))) {
    stop("Error processing stageframe data.")
  }
  
  lifedeathandsex <- specialpatrolgroup(as.matrix(stageexpansion9[, c("fec32n", "rep2n", "indata3221", 
                                                                      "ovgiven_t", "ovgiven_f", "ovest_t", "ovest_f")]), 
                                        as.matrix(stageexpansion3[, c("rep2", "fec32")]), 
                                        as.matrix(data[, c("alive3", "usedfec2")]), stageexpansion9$index3221, 
                                        stageexpansion9$stagepair2n1, stageexpansion3$pairindex, 
                                        stageexpansion3$stage3, stageexpansion3$stage2, data$index3221, 
                                        data$pairindex21, max(stageframe$stageno))
  
  lifeanddeath <- lifedeathandsex[(aliverows + 1),]
  stageexpansion3 <- stageexpansion3[-(unique(c(which(stageexpansion3$stage3 == dim(stageframe)[1]), which(stageexpansion3$stage2 == dim(stageframe)[1])))),]
  
  total.matrix.dim <- dim(stageexpansion3)[1] 
  
  matrix.u <- matrix(lifeanddeath[,1], nrow = total.matrix.dim, ncol = total.matrix.dim)
  
  matrix.f <- matrix(lifeanddeath[,2], nrow = total.matrix.dim, ncol = total.matrix.dim)
  
  matrix.a <- matrix.u + matrix.f
  
  new.matrix <- list(A = matrix.a, U = matrix.u, F = matrix.f)
  
  return(new.matrix)
}

#' Create Raw Ahistorical Population Projection Matrices
#'
#' \code{rlefko2()} returns an ahistorical population projection matrix or, if
#' multiple years are given, then a list of such matrices, as well as the  
#' associated component transition and fecundity matrices, a data frame showing 
#' the pairing ofsingle ahistorical stages used to create historical paired  
#' stages, a dataframe detailing the associated ahistorical stages, and a 
#' vector of years corresponding to the matrices (time \emph{t}).
#' 
#' @param data A vertical demographic data frame, with variables corresponding 
#' to the naming conventions in \code{\link{verticalize3}()}.
#' @param stageframe A stageframe object that includes information on the size,
#' observation status, propagule status, immaturity status, and maturity status
#' of each ahistorical stage. Should also incorporate bin widths if size is
#' continuous.
#' @param year A variable corresponding to year or observation time, or a set
#' of such values, given in values associated with the year term used in linear 
#' model development. Can also equal \code{all}, in which case matrices will
#' be estimated for all years. Defaults to \code{all}.
#' @param pop A variable designating which populations will have matrices 
#' estimated. Should be set to specific population names, or to \code{all} if all 
#' populations should have matrices estimated.
#' @param patch A variable designating which patches or subpopulations will
#' have matrices estimated. Should be set to specific patch names, or to 
#' \code{all} if all patches should have matrices estimated.
#' @param censor If TRUE, then data will be removed according to the variable
#' set in \code{censorcol}, such that only data with censor values equal to 1
#' will remain. Defaults to FALSE.
#' @param stages An optional but important vector denoting the names of the
#' variables within the main vertical dataset coding for the names of the stages
#' of each individual in times \emph{t}+1 and \emph{t}, respectively. The names of
#' the stages in these variables should match those used in the \code{stageframe}
#' input for this analysis exactly. If left blank, then \code{rlefko2()} will attempt
#' to infer stages by matching values of \code{alive}, \code{size}, \code{repst}, and \code{matst}
#' to characteristics noted in the associated \code{stageframe}.
#' @param alive A vector of names of binomial variables corresponding to status 
#' as alive (1) or dead (0) in times \emph{t}+1, \emph{t}, and \emph{t}-1, respectively.
#' @param size A vector of names of variables coding size in times \emph{t}+1 and \emph{t},
#' respectively. Defaults to \code{c("sizea3", "sizea2")}.
#' @param repst A vector of names of variables coding reproductive status in 
#' times \emph{t}+1 and \emph{t}, respectively. Defaults to \code{c("repstatus3", 
#' "repstatus2")}.
#' @param matst A vector of names of variables coding maturity status in 
#' times \emph{t}+1 and \emph{t}, respectively. Defaults to \code{c("matstatus3", 
#' "matstatus2")}. Must be supplied if \code{stages} is not provided.
#' @param fec A vector of names of variables coding fecundity in times \emph{t}+1
#' and \emph{t}, respectively. Defaults to \code{c("feca3", "feca2")}.
#' @param repmatrix A matrix composed mostly of 0s, with non-zero values for
#' each potentially new individual (row) born to each reproductive stage
#' (column). Non-zero entries correspond to multipliers for fecundity, with 1
#' equaling full fecundity.
#' @param overwrite A data frame developed with the \code{\link{overwrite}()} function,
#' describing transitions to be overwritten either with given values or with 
#' other estimated transitions.
#' @param yearcol The variable name or column number corresponding to year 
#' in time \emph{t} in the dataset.
#' @param popcol The variable name or column number corresponding to the
#' identity of the population.
#' @param patchcol The variable name or column number corresponding to patch in 
#' the dataset.
#' @param indivcol The variable name or column number coding individual 
#' identity.
#' @param censorcol The variable name or column number denoting the
#' censor status. Only needed if \code{censor = TRUE}.
#' @param reduce A logical value denoting whether to remove historical stages
#' associated with only zero transitions. These are removed only if all row
#' and column sums in ALL matrices estimated equal 0. Defaults to FALSE.
#' 
#' @return If all inputs are properly formatted, then this function will return
#' an object of class \code{lefkoMat}. This includes:
#' 
#' \item{A}{A list of full projection matrices in order of sorted populations,
#' patches, and years.}
#' \item{U}{A list of survival-transition matrices sorted as in \code{A}.}
#' \item{F}{A list of fecundity matrices sorted as in \code{A}.}
#' \item{hstages}{Null for ahistorical matrices.}
#' \item{ahstages}{A data frame detailing the characteristics of associated
#' ahistorical stages.}
#' \item{labels}{A data frame giving the population, patch, and year of each 
#' matrix in order.}
#' \item{matrixqc}{A short vector describing the number of non-zero elements
#' in \code{U} and \code{F} matrices, and the number of annual matrices.}
#' \item{dataqc}{A vector showing the numbers of individuals and rows in the
#' vertical dataset used as input.}
#'
#' @examples
#' data(cypdata)
#' 
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", 
#'                  "Sm", "Md", "Lg", "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
#'                           repstatus = repvector, obsstatus = obsvector,
#'                           matstatus = matvector, propstatus = propvector,
#'                           immstatus = immvector, indataset = indataset,
#'                           binhalfwidth = binvec)
#' 
#' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
#'                           patchidcol = "patch", individcol = "plantid",
#'                           blocksize = 4, size1col = "Inf2.04", size2col = "Inf.04",
#'                           size3col = "Veg.04", repstr1col = "Inf.04",
#'                           repstr2col = "Inf2.04", fec1col = "Pod.04",
#'                           stageassign = cypframe_raw, stagesize = "sizeadded",
#'                           NAas0 = TRUE, NRasRep = TRUE)
#' 
#' rep_cyp_raw <- matrix(0, 11, 11)
#' rep_cyp_raw[1:2,7:11] <- 0.5
#' 
#' cypover2r <- overwrite(stage3 = c("SD", "P1", "P2", "P3", "D", "XSm", "Sm", "SL", "SL"),
#'                        stage2 = c("SD", "SD", "P1", "P2", "P3", "P3", "P3", "P3", "SL"),
#'                        eststage3 = c(NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'                        eststage2 = c(NA, NA, NA, NA, "D", "D", "D", NA, NA),
#'                        givenrate = c(0.1, 0.2, 0.2, 0.2, NA, NA, NA, 0.25, 0.4),
#'                        type = c("S", "S", "S", "S", "S", "S", "S", "S", "S"))
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, year = "all",
#'                        patch = "all", stages = c("stage3", "stage2", "stage1"),
#'                        size = c("size3added", "size2added"),
#'                        repmatrix = rep_cyp_raw, overwrite = cypover2r,
#'                        yearcol = "year2", patchcol = "patchid",
#'                        indivcol = "individ")
#' cypmatrix2r$A[[1]]

#' @export
rlefko2 <- function(data, stageframe, year = "all", pop = NA, patch = NA, censor = FALSE, 
                    stages = NA, alive = c("alive3", "alive2"), size = c("sizea3", "sizea2"), 
                    repst = c("repstatus3", "repstatus2"), matst = c("matstatus3", "matstatus2"),
                    fec = c("feca3", "feca2"), repmatrix = NA, overwrite = NA, yearcol = NA, popcol = NA, 
                    patchcol = NA, indivcol = NA, censorcol = NA, reduce = FALSE) {
  
  tocensor <- indataset <- alive2 <- popused <- patchused <- yearused <- NULL
  
  if (all(is.na(data))) {
    stop("Need original vertical dataset to proceed.")
  }
  
  if (all(is.na(stages))) {
    if (!(length(alive) > 1)) {
      stop("This function requires stage information. In the absence of stage columns in the dataset, it requires two variables for living/dead status, size, reproductive status, and maturity status, for each of times t+1 and t.")
    }
    if (!(length(size) > 1)) {
      stop("This function requires stage information. In the absence of stage columns in the dataset, it requires two variables for living/dead status, size, reproductive status, and maturity status, for each of times t+1 and t.")
    }
    if (!all(is.na(repst))) {
      if (!(length(repst) > 1)) {
        stop("This function requires stage information. In the absence of stage columns in the dataset, it requires two variables for living/dead status, size, reproductive status, and maturity status, for each of times t+1 and t.")
      }
    }   
    if (!all(is.na(matst))) {
      if (!(length(matst) > 1)) {
        stop("This function requires stage information. In the absence of stage columns in the dataset, it requires two variables for living/dead status, size, reproductive status, and maturity status, for each of times t+1 and t.")
      }
    }   
  }
  
  if (!(length(fec) > 1)) {
    stop("This function requires two variables for fecundity, for each of times t+1 and t.")
  }
  
  if (any(is.character(year)) & any(class(data) == "hfvdata")) {
    if (is.element(tolower(year), "all")) {
      if (is.character(yearcol)) {choicevar <- which(names(data) == yearcol);
      year <- sort(unique(data[,choicevar]))
      yearcol <- choicevar
      } else if (all(is.numeric(year))) {
        year <- sort(unique(data[,yearcol]))
      } else {stop("Cannot understand year designation.")}
    }
  }
  
  if (length(year) == 0 | all(is.na(year) == TRUE)) {
    stop("This function cannot proceed without being given a specific year, or a suite of years.")
  }
  
  if (!all(is.element(year, sort(unique(data[,yearcol]))))) {
    stop("Dataset does not contain one or more of the requested years.")
  }
  
  if (all(is.na(repmatrix))) {    
    repmatrix <- matrix(0, dim(stageframe)[1], dim(stageframe)[1])
    repstages <- which(stageframe$repstatus == 1)
    repmatrix[1, repstages] <- 1
  }
  
  if (censor == TRUE) {
    if(all(is.na(censorcol)) == TRUE) {stop("Cannot censor the data without a proper censor variable.")}
    
    for (i in c(1:length(censorcol))) {
      if (is.character(censorcol)) {
        data$tocensor <- data[,which(names(data) == censorcol[i])]} else {data$tocensor <- data[,censorcol[i]]
        }
      data <- subset(data, tocensor == 1)
    }   
  }
  
  if (!all(is.na(pop)) & !all(is.na(patch))) {
    if (is.na(popcol) | is.na(patchcol)) {stop("Need population and patch designation variables to proceed.")}
    
    if (is.element(tolower(pop), "all")) {
      if (is.character(popcol)) {popcol <- which(names(data) == popcol)}
      
      pops <- unique(data[,popcol])
    } else {pops <- pop}
    
    if (is.element(tolower(patch), "all")) {
      if (is.character(patchcol)) {patchcol <- which(names(data) == patchcol)}
      
      listofpatches <- apply(as.matrix(pops), 1, function(X) {
        patchfolly <- subset(data, popcol == X);
        output <- cbind.data.frame(X, unique(patchfolly[,yearcol]));
        names(output) <- c("pop", "patch");
        return(output);
      })
      
      if (length(listofpatches) > 1) {
        listofpatches <- do.call(rbind.data.frame, listofpatches)
      }
    } else {listofpatches <- expand.grid(pop = pops, patch = patch)}
    
    listofyears <- apply(as.matrix(listofpatches), 1, function(X) {
      checkyrdata <- subset(data, popcol = X[1]);
      checkyrdata <- subset(checkyrdata, patchcol = X[2])
      output <- cbind.data.frame(X[1], X[2], unique(checkyrdata[,yearcol]));
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
  } else if (all(is.na(pop)) & !all(is.na(patch))) {
    if (is.na(patchcol)) {stop("Need patch designation variable to proceed.")}
    
    if (is.element(tolower(patch), "all")) {
      if (is.character(patchcol)) {patchcol <- which(names(data) == patchcol)}
      
      patches <- unique(data[,patchcol])
    } else {patches <- patch}
    
    listofyears <- apply(as.matrix(patches), 1, function(X) {
      checkyrdata <- subset(data, patchcol = X);
      output <- cbind.data.frame(NA, X, unique(checkyrdata[,yearcol]));
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    if (length(listofyears) > 1) {
      listofyears <- do.call(rbind.data.frame, listofyears)
    }
  } else if (!all(is.na(pop)) & all(is.na(patch))) {
    if (is.na(popcol)) {stop("Need population designation variable to proceed.")}
    
    if (is.element(tolower(pop), "all")) {
      if (is.character(popcol)) {popcol <- which(names(data) == popcol)}
      
      pops <- unique(data[,popcol])
    } else {pops <- pop}
    
    listofyears <- apply(as.matrix(pops), 1, function(X) {
      checkyrdata <- subset(data, popcol = X);
      output <- cbind.data.frame(X, NA, unique(checkyrdata[,yearcol]));
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    if (length(listofyears) > 1) {
      listofyears <- do.call(rbind.data.frame, listofyears)
    }
  } else if (all(is.na(pop)) & all(is.na(patch))) {
    listofyears <- cbind.data.frame(NA, NA, year)
    names(listofyears) <- c("pop", "patch", "year2")
  }
  
  melchett <- .sf_reassess(stageframe, repmatrix, overwrite)
  stageframe <- melchett[[1]]
  repmatrix <- melchett[[2]]
  
  stageframe$new_stage_id <- as.numeric(stageframe$new_stage_id)
  stageframe$original_size <- as.numeric(stageframe$original_size)
  stageframe$bin_size_min <- as.numeric(stageframe$bin_size_min)
  stageframe$bin_size_ctr <- as.numeric(stageframe$bin_size_ctr)
  stageframe$bin_size_max <- as.numeric(stageframe$bin_size_max)
  stageframe$bin_size_width <- as.numeric(stageframe$bin_size_width)
  stageframe$bin_raw_halfwidth <- as.numeric(stageframe$bin_raw_halfwidth)
  stageframe$repstatus <- as.numeric(stageframe$repstatus)
  stageframe$obsstatus <- as.numeric(stageframe$obsstatus)
  stageframe$propstatus <- as.numeric(stageframe$propstatus)
  stageframe$immstatus <- as.numeric(stageframe$immstatus)
  stageframe$matstatus <- as.numeric(stageframe$matstatus)
  stageframe$indataset <- as.numeric(stageframe$indataset)
  stageframe$alive <- as.numeric(stageframe$alive)
  
  stageframe$stageno <- c(1:dim(stageframe)[1])
  rownames(stageframe) <- stageframe$stageno
  
  data$alive2 <- data[,which(names(data) == alive[2])]
  data$alive3 <- data[,which(names(data) == alive[1])]
  
  instageframe <- subset(stageframe, indataset == 1)
  
  if (all(is.na(stages))) {
    if (length(size) > 1) {
      data$usedsize2 <- data[,which(names(data) == size[2])]
      data$usedsize3 <- data[,which(names(data) == size[1])]
    } else {
      warning("Without stage columns, lefko3 MPM estimation functions generally require size variables. Failure to include size variables may lead to odd results.")
    }
    if (length(repst) > 1) {
      data$usedrepst2 <- data[,which(names(data) == repst[2])]
      data$usedrepst3 <- data[,which(names(data) == repst[1])]
    } 
    if (length(matst) > 1) {
      data$usedmatstatus2 <- data[,which(names(data) == matst[2])]
      data$usedmatstatus3 <- data[,which(names(data) == matst[1])]
    } 
    
    data$usedstage2 <- apply(as.matrix(c(1:dim(data)[1])), 1, function(X) {
      if (is.na(data$usedsize2[X])) {
        data$usedsize2[X] <- 0
      }
      mainstages <- intersect(which(instageframe$bin_size_min < data$usedsize2[X]), 
                              which(instageframe$bin_size_max >= data$usedsize2[X]))
      jmstages <- which(instageframe$immstatus == (1 - data$usedmatstatus2[X]))
      obsstages <- which(instageframe$obsstatus == data$obsstatus2[X])
      repstages <- which(instageframe$repstatus == data$repstatus2[X])
      
      choicestage <- intersect(intersect(mainstages, jmstages), intersect(obsstages, repstages))

      if (length(choicestage) == 0) {
        stop("Stage characteristics mismatch dataset. Consider using the stages option, particularly if the vertical file was created with NRasRep = TRUE in verticalize3() or historicalize3().")
      }
      
      return(instageframe$orig_stage_id[choicestage])
    })
    
    data$usedstage3 <- apply(as.matrix(c(1:dim(data)[1])), 1, function(X) {
      if (is.na(data$usedsize3[X])) {
        data$usedsize3[X] <- 0
      }
      mainstages <- intersect(which(instageframe$bin_size_min < data$usedsize3[X]), 
                              which(instageframe$bin_size_max >= data$usedsize3[X]))
      jmstages <- which(instageframe$immstatus == (1 - data$usedmatstatus3[X]))
      obsstages <- which(instageframe$obsstatus == data$obsstatus3[X])
      repstages <- which(instageframe$repstatus == data$repstatus3[X])
      alivestage3 <- which(instageframe$alive == data$alive3[X])
      
      choicestage <- intersect(intersect(intersect(mainstages, jmstages), intersect(obsstages, repstages)), alivestage3)

      if (length(choicestage) == 0) {
        stop("Stage characteristics mismatch dataset. Consider using the stages option, particularly if the vertical file was created with NRasRep = TRUE in verticalize3() or historicalize3().")
      }
      
      if (length(choicestage) == 0) choicestage <- which(instageframe$new_stage_id == max(instageframe$new_stage_id))
      
      return(instageframe$orig_stage_id[choicestage])
    })
    
  } else if (length(stages) > 1) {
    if (is.numeric(stages[2])) {
      data$usedstage2 <- data[, stages[2]]
      data$usedstage3 <- data[, stages[1]]
      data$usedstage3[which(data$usedstage3 == "NotAlive")] <- "Dead"
    } else {
      data$usedstage2 <- data[,which(names(data) == stages[2])]
      data$usedstage3 <- data[,which(names(data) == stages[1])]
      data$usedstage3[which(data$usedstage3 == "NotAlive")] <- "Dead"
    }
    stages.used <- sort(unique(c(data$usedstage2, data$usedstage3)))
    
    if (length(setdiff(stages.used, stageframe$orig_stage_id)) > 0) {
      stop("Some stages in dataset do not match those detailed in the input stageframe.", .call = FALSE)
    }
  }
  
  if (length(fec) > 1) {
    data$usedfec2 <- data[,which(names(data) == fec[2])]
    data$usedfec3 <- data[,which(names(data) == fec[1])]
    
    data$usedfec2[which(is.na(data$usedfec2))] <- 0
    data$usedfec3[which(is.na(data$usedfec3))] <- 0
  } else {
    warning("Lefko3 MPM estimation functions generally require fecundity variables. Failure to include fecundity variables leads to matrices composed only of survival transitions.")
  } 
  
  stageexpansion3 <- cbind.data.frame(expand.grid(size3 = stageframe$bin_size_ctr, size2 = stageframe$bin_size_ctr), 
                                      expand.grid(rep3 = stageframe$repstatus, rep2 = stageframe$repstatus),
                                      expand.grid(indata3 = stageframe$indataset, indata2 = stageframe$indataset),
                                      expand.grid(stage3 = stageframe$stageno, stage2 = stageframe$stageno),
                                      fec32 = c(cbind(rbind(repmatrix, 0), 0)))
  
  instages <- length(stageframe$new_stage_id)
  
  stageexpansion3$indata32 <- stageexpansion3$indata3 * stageexpansion3$indata2
  stageexpansion3$index32 <- (stageexpansion3$stage3 - 1) + ((stageexpansion3$stage2 - 1) * instages)
  stageexpansion3$index2 <- stageexpansion3$stage2 - 1
  
  stageexpansion3$ovest_t <- 0
  stageexpansion3$ovgiven_t <- 0
  stageexpansion3$ovest_f <- 0
  stageexpansion3$ovgiven_f <- 0
  
  if (!all(is.na(overwrite))) {
    
    overwrite <- .overwrite_reassess(overwrite, stageframe, historical = FALSE)
    if (dim(overwrite)[1] == 0) overwrite <- NA
  }
  
  if (!all(is.na(overwrite))) {
    
    overwrite$index3 <- apply(as.matrix(c(1:dim(overwrite)[1])), 1, function(X) {
      stageframe$stageno[which(stageframe$orig_stage_id == overwrite$stage3[X])]})
    overwrite$index2 <- apply(as.matrix(c(1:dim(overwrite)[1])), 1, function(X) {
      stageframe$stageno[which(stageframe$orig_stage_id == overwrite$stage2[X])]})
    overwrite$new3 <- apply(as.matrix(c(1:dim(overwrite)[1])), 1, function(X) {
      if(!is.na(overwrite$eststage3[X])) {
        stageframe$stageno[which(stageframe$orig_stage_id == overwrite$eststage3[X])]
      } else {
        overwrite$index3[X]
      }
    })
    overwrite$new2 <- apply(as.matrix(c(1:dim(overwrite)[1])), 1, function(X) {
      if(!is.na(overwrite$eststage2[X])) {
        stageframe$stageno[which(stageframe$orig_stage_id == overwrite$eststage2[X])]
      } else {
        overwrite$index2[2]
      }
    })
    overwrite$indexold32 <- (overwrite$index3 - 1) + ((overwrite$index2 - 1) * instages)
    overwrite$indexnew32 <- (overwrite$new3 - 1) + ((overwrite$new2 - 1) * instages)
    overwrite$givenrate[which(is.na(overwrite$givenrate))] <- 0
    
    overwrite$new3[which(is.na(overwrite$eststage3))] <- -1
    overwrite$givenrate[which(is.na(overwrite$givenrate))] <- -1
    overwrite$givenrate[which(overwrite$givenrate == 0)] <- -1
    
    allstageadditions <- ovreplace(stageexpansion3$index32, overwrite$indexold32, overwrite$indexnew32, 
                                   overwrite$convtype, overwrite$new3, overwrite$givenrate)
    
    stageexpansion3$ovest_t <- allstageadditions[,1]
    stageexpansion3$ovgiven_t <- allstageadditions[,2]
    stageexpansion3$ovest_f <- allstageadditions[,3]
    stageexpansion3$ovgiven_f <- allstageadditions[,4]
    
  }
  
  stageexpansion2 <- cbind.data.frame(stage2 = as.numeric(stageframe$stageno), size2 = as.numeric(stageframe$bin_size_ctr), 
                                      rep2 = as.numeric(stageframe$repstatus), indata2 = as.numeric(stageframe$indataset),
                                      index2 = (as.numeric(stageframe$stageno) - 1), fec3 = c(rowSums(repmatrix), 0))
  stageexpansion2$fec3[which(stageexpansion2$fec3 > 0)] <- 1
  
  data <- subset(data, alive2 == 1)
  
  data$index2 <- apply(as.matrix(data$usedstage2), 1, function(X) {
    instageframe$stageno[which(instageframe$orig_stage_id == X)] - 1
  })
  data$index2[which(is.na(data$index2))] <- 0
  data$index3 <- apply(as.matrix(data$usedstage3), 1, function(X) {
    instageframe$stageno[which(instageframe$orig_stage_id == X)] - 1
  })
  data$index3[which(is.na(data$index3))] <- 0
  data$index32 <- apply(as.matrix(c(1:length(data$usedstage2))), 1, function(X) {
    (data$index3[X] + (data$index2[X] * instages))
  })
  
  if(is.element(0, unique(data$index2))) {
    warning("Data (stage at time t) contains some stages not identified in stageframe. Note that all stage characteristics must match, including reproductive status.")
  }
  if(is.element(0, unique(data$index3))) {
    warning("Data (stage at time t+1) contains some stages not identified in stageframe. Note that all stage characteristics must match, including reproductive status.")
  }
  
  aliverows <- hoffmannofstuttgart(0, as.matrix(stageexpansion3[, c("stage3", "stage2", "stage2", "stage2")]))
  
  madsexmadrigal <- apply(listofyears, 1, function(X) {
    passed_data <- data
    if (!is.na(X[1])) {
      passed_data$popused <- passed_data[,popcol];
      passed_data <- subset(passed_data, popused == X[1]);
    }
    if (!is.na(X[2])) {
      passed_data$patchused <- passed_data[,patchcol];
      passed_data <- subset(passed_data, patchused == X[2]);
    }
    if (!is.na(X[3])) {
      passed_data$yearused <- passed_data[,yearcol];
      passed_data <- subset(passed_data, yearused == X[3]);
    }
    .rlefko2_core(data = passed_data, stageframe = stageframe, stageexpansion3 = stageexpansion3, 
                  aliverows = aliverows, stageexpansion2 = stageexpansion2)
  })
  
  a_list <- lapply(madsexmadrigal, function(X) {X$A})
  u_list <- lapply(madsexmadrigal, function(X) {X$U})
  f_list <- lapply(madsexmadrigal, function(X) {X$F})
  
  qcoutput1 <- NA
  qcoutput2 <- NA
  
  totalutransitions <- sum(unlist(lapply(u_list, function(X) {length(which(X != 0))})))
  totalftransitions <- sum(unlist(lapply(f_list, function(X) {length(which(X != 0))})))
  totalmatrices <- length(u_list)
  
  qcoutput1 <- c(totalutransitions, totalftransitions, totalmatrices)
  
  indivs <- NA
  if (!all(is.na(indivcol))) {
    if (all(is.character(indivcol))) {indivcol <- which(names(data) == indivcol)[1]}
    indivs <- length(unique(data[,indivcol]))
  }
  qcoutput2 <- c(indivs, dim(data)[1])
  
  ahstages <-  stageframe[1:(dim(stageframe)[1] - 1),]
  
  if (reduce == TRUE) {
    drops <- .reducer2(a_list, u_list, f_list, ahstages)
    
    a_list <- drops$A
    u_list <- drops$U
    f_list <- drops$F
    ahstages <- drops$ahstages
  }
  
  output <- list(A = a_list, U = u_list, F = f_list, hstages = NULL, ahstages = ahstages, 
                 labels = listofyears, matrixqc = qcoutput1, dataqc = qcoutput2)
  class(output) <- "lefkoMat"
  
  return(output)
}

#' Core Wrapper Powering Ahistorical Raw Matrix Estimation
#' 
#' \code{.rlefko2_core()} pulls together all required data and matrix conditions
#' and feeds them into\code{\link{normalpatrolgroup}()}, an Rcpp function designed
#' to quickly estimate ahistorical raw matrices.
#' 
#' @param data Original demographic dataset showing the states and fates of
#' individuals.
#' @param stageframe Original stageframe used in analysis.
#' @param stageexpansion3 Data frame created by \code{\link{flefko2}()} including
#' general life history and statistical characteristics identifying every element
#' in the projection matrix.
#' @param aliverows Vector identifying elements in the projection matrix that are
#' both alive and logically possible.
#' @param stageexpansion2 Data frame characterizing all stages used in
#' ahistorical matrix estimation.
#' 
#' @return Returns a list containing one set of matrices corresponding to a 
#' specific population, patch, and time step combination.
#' 
#' @keywords internal
#' @noRd
.rlefko2_core <- function(data, stageframe, stageexpansion3 = NA, aliverows = NA, stageexpansion2 = NA) {
  
  if (all(is.na(stageexpansion3)) | all(is.na(stageexpansion2))) {
    stop("Error processing stageframe data.")
  }
  
  lifedeathandsex <- normalpatrolgroup(as.matrix(stageexpansion3[, c("fec32", "rep2", "indata32",
                                                                     "ovgiven_t", "ovgiven_f", "ovest_t", "ovest_f")]), 
                                       as.matrix(stageexpansion2[, c("rep2", "fec3")]), 
                                       as.matrix(data[, c("alive3", "usedfec2")]), stageexpansion3$index32, 
                                       stageexpansion3$index2, stageexpansion2$index2, 
                                       stageexpansion2$stage2, data$index32, 
                                       data$index2, max(stageframe$stageno))
  
  lifeanddeath <- lifedeathandsex[(aliverows + 1),]
  stageexpansion2 <- stageexpansion2[-(unique(which(stageexpansion2$stage2 == dim(stageframe)[1]))),]
  
  total.matrix.dim <- dim(stageexpansion2)[1] 
  
  matrix.u <- matrix(lifeanddeath[,1], nrow = total.matrix.dim, ncol = total.matrix.dim)
  
  matrix.f <- matrix(lifeanddeath[,2], nrow = total.matrix.dim, ncol = total.matrix.dim)
  
  matrix.a <- matrix.u + matrix.f
  
  new.matrix <- list(A = matrix.a, U = matrix.u, F = matrix.f)
  
  return(new.matrix)
}

#' Summary of Class "lefkoMat"
#'
#' A function to simplify the viewing of basic information describing the matrices
#' produced through functions \code{\link{flefko3}()}, \code{\link{flefko2}()}, \code{\link{rlefko3}()},
#' and \code{\link{rlefko2}()}.
#' 
#' @param object An object of class \code{lefkoMat}.
#' @param ... Other parameters.
#' 
#' @return A summary of the object, showing the number of each type of matrix, the
#' number of annual matrices, the number of estimated (non-zero) elements across
#' all matrices and per matrix, the number of unique transitions in the dataset,
#' and the number of individuals.
#' 
#' @examples
#' data(cypdata)
#' 
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", 
#'                  "Sm", "Md", "Lg", "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
#'                           repstatus = repvector, obsstatus = obsvector,
#'                           matstatus = matvector, propstatus = propvector,
#'                           immstatus = immvector, indataset = indataset,
#'                           binhalfwidth = binvec)
#' 
#' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
#'                           patchidcol = "patch", individcol = "plantid",
#'                           blocksize = 4, size1col = "Inf2.04", size2col = "Inf.04",
#'                           size3col = "Veg.04", repstr1col = "Inf.04",
#'                           repstr2col = "Inf2.04", fec1col = "Pod.04",
#'                           stageassign = cypframe_raw, stagesize = "sizeadded",
#'                           NAas0 = TRUE, NRasRep = TRUE)
#' 
#' rep_cyp_raw <- matrix(0, 11, 11)
#' rep_cyp_raw[1:2,7:11] <- 0.5
#' 
#' cypover2r <- overwrite(stage3 = c("SD", "P1", "P2", "P3", "D", "XSm", "Sm", "SL", "SL"),
#'                        stage2 = c("SD", "SD", "P1", "P2", "P3", "P3", "P3", "P3", "SL"),
#'                        eststage3 = c(NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'                        eststage2 = c(NA, NA, NA, NA, "D", "D", "D", NA, NA),
#'                        givenrate = c(0.1, 0.2, 0.2, 0.2, NA, NA, NA, 0.25, 0.4),
#'                        type = c("S", "S", "S", "S", "S", "S", "S", "S", "S"))
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, year = "all",
#'                        patch = "all", stages = c("stage3", "stage2", "stage1"),
#'                        size = c("size3added", "size2added"),
#'                        repmatrix = rep_cyp_raw, overwrite = cypover2r,
#'                        yearcol = "year2", patchcol = "patchid",
#'                        indivcol = "individ")
#' summary(cypmatrix2r)
#' 
#' @export
summary.lefkoMat <- function(object, ...) {
  
  matrices <- object
  
  matdim <- dim(matrices$A[[1]])[1]
  
  mqca <- matrices$matrixqc[1]
  mqcb <- matrices$matrixqc[2]
  mqcc <- matrices$matrixqc[3]
  
  writeLines(paste0("\nThis lefkoMat object contains ", mqcc, " matrices."))
  writeLines(paste0("\nEach matrix is a square matrix with ", matdim, " rows and columns, and a total of ", matdim*matdim, " elements."))
  writeLines(paste0("A total of ", mqca, " survival transitions were estimated, with ", 
                    mqca / mqcc, " per matrix."))
  writeLines(paste0("A total of ", mqcb, " fecundity transitions were estimated, with ", 
                    mqcb / mqcc, " per matrix."))
  
  if (is.element("dataqc", names(matrices))) {
    dqca <- matrices$dataqc[1]
    dqcb <- matrices$dataqc[2]
    
    writeLines(paste0("\nThe dataset contains a total of ", dqca, " unique individuals and ", dqcb, " unique transitions."))
  }
  
  if (is.element("modelqc", names(matrices))) {
    moqc12 <- matrices$modelqc[1,2]
    moqc22 <- matrices$modelqc[2,2]
    moqc32 <- matrices$modelqc[3,2]
    moqc42 <- matrices$modelqc[4,2]
    moqc52 <- matrices$modelqc[5,2]
    moqc62 <- matrices$modelqc[6,2]
    moqc72 <- matrices$modelqc[7,2]
    moqc82 <- matrices$modelqc[8,2]
    moqc92 <- matrices$modelqc[9,2]
    
    moqc13 <- matrices$modelqc[1,3]
    moqc23 <- matrices$modelqc[2,3]
    moqc33 <- matrices$modelqc[3,3]
    moqc43 <- matrices$modelqc[4,3]
    moqc53 <- matrices$modelqc[5,3]
    moqc63 <- matrices$modelqc[6,3]
    moqc73 <- matrices$modelqc[7,3]
    moqc83 <- matrices$modelqc[8,3]
    moqc93 <- matrices$modelqc[9,3]
    
    writeLines("\nVital rate modeling quality control:\n")
    writeLines(paste0("Survival estimated with ", moqc12, " individuals and ", moqc13, " individual transitions."))
    writeLines(paste0("Observation estimated with ", moqc22, " individuals and ", moqc23, " individual transitions."))
    writeLines(paste0("Size estimated with ", moqc32, " individuals and ", moqc33, " individual transitions."))
    writeLines(paste0("Reproductive status estimated with ", moqc42, " individuals and ", moqc43, " individual transitions."))
    writeLines(paste0("Fecundity estimated with ", moqc52, " individuals and ", moqc53, " individual transitions."))
    writeLines(paste0("Juvenile survival estimated with ", moqc62, " individuals and ", moqc63, " individual transitions."))
    writeLines(paste0("Juvenile observation estimated with ", moqc72, " individuals and ", moqc73, " individual transitions."))
    writeLines(paste0("Juvenile size estimated with ", moqc82, " individuals and ", moqc83, " individual transitions."))
    writeLines(paste0("Juvenile reproductive status estimated with ", moqc92, " individuals and ", moqc93, " individual transitions."))
  }
  return()
}

