#' Develop Best-fit Vital Rate Estimation Models For Matrix Development
#' 
#' \code{modelsearch()} returns both a table of vital rate estimating models and a
#' best-fit model for each major vital rate estimated. The final output can be 
#' used as input in other functions within this package.
#' 
#' @param data The vertical dataset to be used for analysis. The dataset should 
#' ideally be of class \code{hfvdata}, but will also work with data frames 
#' formatted similarly to the output format provided by functions 
#' \code{\link{verticalize3}()} or \code{\link{historicalize3}()}, as long as all 
#' needed variables are properly designated.
#' @param historical A logical variable denoting whether to assess the effects
#' of state in time \emph{t}-1 in addition to state in time \emph{t}. Defaults 
#' to TRUE.
#' @param approach Designates the approach to be taken for model building. The 
#' default is \code{lme4}, which uses the mixed model approach utilized in 
#' package 'lme4'. Other options include \code{glm}, which uses the \code{lm}, \code{glm},
#' and related functions in base R.
#' @param suite This describes the global model for each vital rate estimation
#' and has the following possible values: \code{full}, includes main effects and 
#' all two-way interactions of size and reproductive status; \code{main}, 
#' includes main effects only of size and reproductive status; \code{size}, 
#' includes only size (also interactions between size in historical model); 
#' \code{rep}, includes only reproductive status (also interactions between 
#' status in historical model); \code{cons}, all vital rates estimated only as 
#' y-intercepts. If \code{approach = "glm"} and \code{year.as.random = FALSE}, then 
#' year is also included as a fixed effect, and, in the case of \code{full}, 
#' included in two-way interactions. Defaults to \code{size}.
#' @param vitalrates A vector describing which vital rates will be estimated
#' via linear modeling, with the following options: \code{surv}, survival 
#' probability; \code{obs}, observation probability; \code{size}, overall size; \code{repst}, 
#' probability of reproducing; and \code{fec}, amount of reproduction (overall 
#' fecundity). Defaults to \code{c("surv", "size", "fec")}.
#' @param juvestimate An optional variable denoting the stage name of the 
#' juvenile stage in the vertical dataset. If not NA, and \code{stage} is also 
#' given (see below), then vital rates listed in \code{vitalrates} other than 
#' \code{fec} will also be estimated from the juvenile stage to all adult stages. 
#' Defaults to NA, in which case juvenile vital rates are not estimated.
#' @param juvsize A logical variable denoting whether size should be used as
#' a term in models involving transition from the juvenile stage. Defaults 
#' to FALSE, and is only used if \code{juvestimate} does not equal NA.
#' @param bestfit A variable indicating the model selection criterion for the 
#' choice of best-fit model. The default is \code{AICc&k}, which chooses the 
#' best-fit model as the model with the lowest AICc or, if not the same 
#' model, then the model that has the lowest degrees of freedom among models 
#' with delta AICc <= 2.0. Alternatively, \code{AICc} may be chosen, in which case
#' the best-fit model is simply the model with the lowest AICc value.
#' @param sizedist The probability distribution used to model size. Options 
#' include \code{gaussian} for the Normal distribution (default), \code{poisson} 
#' for the Poisson distribution, and \code{negbin} for the negative binomial
#' distribution.
#' @param fecdist The probability distribution used to model fecundity. Options
#' include \code{gaussian} for the Normal distribution (default), \code{poisson} 
#' for the Poisson distribution, and \code{negbin} for the negative binomial 
#' distribution.
#' @param fectime A variable indicating which year of fecundity to use as 
#' the response term in fecundity models. Options include \code{2}, which refers 
#' to time \emph{t}, and \code{3}, which refers to time \emph{t}+1. Defaults to \code{2}.
#' @param censor A vector denoting the names of censoring variables in the 
#' dataset, in order from time \emph{t}+1, followed by time \emph{t}, 
#' and lastly followed by time \emph{t}-1. Defaults to NA.
#' @param indiv A variable indicating the variable name coding individual 
#' identity. Defaults to \code{individ}.
#' @param patch A variable indicating the variable name coding for patch, 
#' where patches are defined as permanent subgroups within the study
#' population. Defaults to NA.
#' @param year A variable indicating the variable coding for  observation 
#' time in time \emph{t}. Defaults to \code{year2}.
#' @param surv A vector indicating the variable names coding for status as 
#' alive or dead in times \emph{t}+1, \emph{t}, and \emph{t}-1, respectively. 
#' Defaults to \code{c("alive3", "alive2", "alive1")}.
#' @param obs A vector indicating the variable names coding for observation
#' status in times \emph{t}+1, \emph{t}, and \emph{t}-1, respectively. 
#' Defaults to \code{c("obsstatus3", "obsstatus2", "obsstatus1")}.
#' @param size A vector indicating the variable names coding for size in 
#' times \emph{t}+1, \emph{t}, and \emph{t}-1, respectively. Defaults to 
#' \code{c("sizea3", "sizea2", "sizea1")}.
#' @param repst A vector indicating the variable names coding for 
#' reproductive status in times \emph{t}+1, \emph{t}, and \emph{t}-1, 
#' respectively. Defaults to \code{c("repstatus3", "repstatus2", "repstatus1")}.
#' @param fec A vector indicating the variable names coding for fecundity
#' in times \emph{t}+1, \emph{t}, and \emph{t}-1, respectively. Defaults to 
#' \code{c("feca3", "feca2", "feca1")}.
#' @param stage A vector indicating the variables coding for stage in times 
#' \emph{t}+1, \emph{t}, and \emph{t}-1. Defaults to 
#' \code{c("stage3", "stage2", "stage1")}.
#' @param age Designates the name of the variable corresponding to age in the
#' vertical dataset. Defaults to NA, in which case age is not included in linear
#' models. Should only be used if building age x stage matrices.
#' @param year.as.random If set to TRUE and \code{approach = "lme4"}, then \code{year} is 
#' included as a random factor. If set to FALSE, then \code{year} is included as a 
#' fixed factor. All other combinations of logical value and \code{approach} lead to 
#' \code{year} not being included in modeling.
#' @param patch.as.random If set to TRUE and \code{approach = "lme4"}, then \code{patch}
#' is included as a random factor. If set to FALSE and \code{approach = "glm"}, then 
#' \code{patch} is included as a fixed factor. All other combinations of logical 
#' value and \code{approach} lead to \code{patch} not being included in modeling. 
#' @param show.model.tables If set to TRUE, then includes full modeling tables
#' in the output. Defaults to TRUE.
#' @param quiet If set to TRUE, then model building and selection will proceed
#' without warnings and diagnostic messages being issued. Note that this will
#' not affect warnings and messages generated during global model development.
#' Defaults to FALSE.
#' 
#' @return This function yields an object of class \code{lefkoMod}, which is a 
#' list in which the first 8 elements are the best-fit models for survival, 
#' observation status, size, reproductive status, fecundity, juvenile survival, 
#' juvenile observation, and juvenile size, respectively, followed by eight 
#' elements corresponding to the model tables for each of these vital rates, in
#' order, followed by a single character element denoting the criterion used for
#' model selection, as follows:
#' 
#' \item{survival_model}{Best-fit model of the binomial probability of survival
#' from time \emph{t} to time \emph{t}+1. Defaults to 1.}
#' \item{observation_model}{Best-fit model of the binomial probability of 
#' observation in time \emph{t}+1 given survival to that time. Defaults to 1.}
#' \item{size_model}{Best-fit model of size in time \emph{t}+1 given 
#' survival to and observation in that time. Defaults to 1.}
#' \item{repstatus_model}{Best-fit model of the binomial probability of
#' reproduction in time \emph{t}+1, given survival to and observation in that
#' time. Defaults to 1.}
#' \item{fecundity_model}{Best-fit model of fecundity in time \emph{t}+1 given 
#' survival to, and observation and reproduction in that time. Defaults to 1.}
#' \item{juv_survival_model}{Best-fit model of the binomial probability of 
#' survival from time \emph{t} to time \emph{t}+1 of an immature individual. 
#' Defaults to 1.}
#' \item{juv_observation_model}{Best-fit model of the binomial probability of 
#' observation in time \emph{t}+1 given survival to that time of 
#' an immature individual. Defaults to 1.}
#' \item{juv_size_model}{Best-fit model of size in time \emph{t}+1 given 
#' survival to and observation in that time of an immature individual. Defaults 
#' to 1.}
#' \item{juv_reproduction_model}{Best-fit model of the binomial probability 
#' of reproduction in time \emph{t}+1, given survival to and observation in 
#' that time of an individual that was immature in time \emph{t}. This model is 
#' technically not a model of reproduction probability for individuals that are 
#' immature, rather reproduction probability here is given for individuals that 
#' are mature in time \emph{t}+1 but were immature in time \emph{t}. Defaults 
#' to 1.}
#' \item{survival_table}{Full dredge model table of survival probability.}
#' \item{observation_table}{Full dredge model table of observation
#' probability.}
#' \item{size_table}{Full dredge model table of size.}
#' \item{repstatus_table}{Full dredge model table of reproduction probability.}
#' \item{fecundity_table}{Full dredge model table of fecundity.}
#' \item{juv_survival_table}{Full dredge model table of immature survival 
#' probability.}
#' \item{juv_observation_table}{Full dredge model table of immature 
#' observation probability.}
#' \item{juv_size_table}{Full dredge model table of immature size.}
#' \item{juv_reproduction_table}{Full dredge model table of immature 
#' reproduction probability.}
#' \item{criterion}{Vharacter variable denoting the criterion used to 
#' determine the best-fit model.}
#' \item{qc}{Data frame with three variables: 1) Name of vital rate, 2) number
#' of individuals used to model that vital rate, and 3) number of individual
#' transitions used to model that vital rate.}
#' 
#' The mechanics governing model building are fairly robust to errors and 
#' exceptions. The function attempts to build global models, and simplifies 
#' models automatically via several steps should model building fail. Model 
#' selection proceeds via the \link[MuMIn]{dredge} function in package 'MuMIn', and 
#' defaults to the global model should that fail.
#' 
#' This function is set to run in a verbose fashion, so that any errors and 
#' warnings developed during model building, model analysis, and model selection 
#' can be seen and, if necessary, dealt with. Interpretations of errors during 
#' global model analysis may be found in documentation in base R for functions 
#' \code{lm} and \code{glm} used in analysis of models without random terms, and packages 
#' 'lme4' and 'glmmTMB' for mixed models (see \link[lme4]{glmer} and \link[glmmTMB]{glmmTMB}, respectively). 
#' Package 'MuMIn' is used for model dredging (see \link[MuMIn]{dredge}), and errors 
#' and warnings during dredging can be interpreted using the documentation for 
#' that package. The \code{quiet = TRUE} option can be used to silence dredging 
#' warnings, but users should note that automated model selection can be viewed 
#' as mindless, and so great care should be taken to ensure that the models run
#' are sensical, and that issues of model quality are properly dealt with.
#' 
#' Note that care must be taken to build models that test the impacts of state 
#' in time \emph{t}-1 for historical models, and that do not test these impacts 
#' for ahistorical models. Ahistorical matrix modeling particularly will yield
#' transition estimates biased low if historical terms from models are ignored.
#' This can be at the start of modeling by setting \code{historical = FALSE} for
#' the ahistorical case, and \code{historical = TRUE} for the historical case.
#' 
#' Model building and selection may fail if NAs exist within variables used in
#' modeling. If NAs represent 0 entries, then please change all NAs to 0, as 
#' with the \code{NAas0 = TRUE} option in function \code{\link{verticalize3}()}.
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
#' lathmodelsln2 <- modelsearch(lathvertln, historical = FALSE, approach = "lme4", suite = "main",
#'                              vitalrates = c("surv", "obs", "size", "repst", "fec"), 
#'                              juvestimate = "Sdl", bestfit = "AICc&k", sizedist = "gaussian", 
#'                              fecdist = "poisson", indiv = "individ", patch = "patchid", 
#'                              year = "year2", year.as.random = TRUE, patch.as.random = TRUE,
#'                              show.model.tables = TRUE)
#' 
#' lathmodelsln2
#' }
#' 
#' @export
modelsearch <- function(data, historical = TRUE, approach = "lme4", suite = "size", 
                        vitalrates = c("surv", "size", "fec"), juvestimate = NA,
                        juvsize = FALSE, bestfit = "AICc&k", sizedist = "gaussian", fecdist = "gaussian", 
                        fectime = 2, censor = NA, indiv = "individ", patch = NA, year = "year2", 
                        surv = c("alive3", "alive2", "alive1"), obs = c("obsstatus3", "obsstatus2", "obsstatus1"), 
                        size = c("sizea3", "sizea2", "sizea1"), repst = c("repstatus3", "repstatus2", "repstatus1"), 
                        fec = c("feca3", "feca2", "feca1"), stage = c("stage3", "stage2", "stage1"), 
                        age = NA, year.as.random = TRUE, patch.as.random = TRUE, show.model.tables = TRUE,
                        quiet = FALSE) {
  
  old <- options() #This function requires changes to options(na.action) in order for the lme4::dredge routines to work properly
  on.exit(options(old)) #This will reset options() to user originals when the function exits
  
  censor1 <- censor2 <- censor3 <- NULL
  
  #This first section deals with errors and exceptions in input options, and also sets up some important variables
  if (all(class(data) != "hfvdata")) {warning("This function was made to work with standardized historically-formatted vertical datasets, as provided by the verticalize() and historicalize() functions. Failure to format the input data properly and designate needed variables appropriately may result in nonsensical output.")}
  
  if (!requireNamespace("MuMIn", quietly = TRUE)) {stop("Package MuMIn needed for this function to work. Please install it.", call. = FALSE)}
  if (!requireNamespace("stringr", quietly = TRUE)) {stop("Package stringr needed for this function to work. Please install it.", call. = FALSE)}
  if (approach == "lme4" & !requireNamespace("lme4", quietly = TRUE)) {
    if (sizedist == "negbin" & !requireNamespace("glmmTMB", quietly = TRUE)) {stop("Package glmmTMB needed to develop mixed size models with a negative binomial distribution.")}
    if (fecdist == "negbin" & !requireNamespace("glmmTMB", quietly = TRUE)) {stop("Package glmmTMB needed to develop mixed fecundity models with a negative binomial distribution.")}
    stop("Package lme4 needed for this function to work. Please install it.", call. = FALSE)
  }
  distoptions <- c("gaussian", "poisson", "negbin")
  packoptions <- c("lme4", "glm")
  
  if (!is.element(approach, packoptions)) {stop("Please enter a valid package option, currently either lme4, glm, or nlme.", call. = FALSE)}
  if (!is.element(sizedist, distoptions)) {stop("Please enter a valid assumed size distribution, currently either gaussian, poisson, or negbin.", call. = FALSE)}
  if (!is.element(fecdist, distoptions)) {stop("Please enter a valid assumed fecundity distribution, currently either gaussian, poisson, or negbin.", call. = FALSE)}
  
  if (length(censor) > 3) {
    stop("Censor variables should be included either as 1 variable per row in the historical data file (1 variable in the dataset), or as 1 variable per timestep within each historical data file (2 or 3 variables in the dataset). No more than 3 variables are allowed, and if more than one are supplied, then they are assumed to be in order of time t+1, time t, and time t-1, respectively.", call. = FALSE)
  }
  if (length(indiv) > 1) {stop("Only one individual identification variable is allowed.", call. = FALSE)}
  if (length(year) > 1) {stop("Only one time step variable is allowed, and it must refer to time step t.", call. = FALSE)}
  if (length(patch) > 1) {stop("Only one patch variable is allowed, and it must refer to time step t.", call. = FALSE)}
  
  if (is.element("surv", vitalrates)) {
    if (length(surv) > 3 | length(surv) == 1) {stop("This function requires 2 (if ahistorical) or 3 (if historical) survival variables from the dataset as input parameters.", call. = FALSE)}
  }
  if (is.element("obs", vitalrates)) {
    if (length(obs) > 3 | length(obs) == 1) {stop("This function requires 2 (if ahistorical) or 3 (if historical) observation variables from the dataset as input parameters.", call. = FALSE)}
  }
  if (is.element("size", vitalrates)) {
    if (length(size) > 3 | length(size) == 1) {stop("This function requires 2 (if ahistorical) or 3 (if historical) size variables from the dataset as input parameters.", call. = FALSE)}
  }
  if (is.element("repst", vitalrates)) {
    if (length(repst) > 3 | length(repst) == 1) {stop("This function requires 2 (if ahistorical) or 3 (if historical) reproductive status variables from the dataset as input parameters.", call. = FALSE)}
  }
  if (is.element("fec", vitalrates)) {
    if (length(fec) > 3 | length(fec) == 1) {stop("This function requires 2 (if ahistorical) or 3 (if historical) fecundity variables from the dataset as input parameters.", call. = FALSE)}
  }
  
  if (fectime != 2 & fectime != 3) {
    stop("The fectime option must equal either 2 or 3, depending on whether the fecundity response term is for time t or time t+1, respectively.", call. = FALSE)
  }
  
  if (!is.na(age)) {
    if (length(which(names(data) == age)) == 0) {
      stop("Age must equal either NA or the exact name of the variable denoting age in the dataset.", call. = FALSE)
    } else {
      agecol <- which(names(data) == age)
    }
  }
  
  if (historical == FALSE) { #This portion eliminates time t-1 from analysis in the ahistorical case
    if (length(surv) > 2) {surv <- surv[1:2]}
    if (length(obs) > 2) {obs <- obs[1:2]}
    if (length(size) > 2) {size <- size[1:2]}
    if (length(repst) > 2) {repst <- repst[1:2]}
    if (length(fec) > 2) {fec <- fec[1:2]}
  }
  
  #The next section creates the text lines needed for the main model calls, based on function input
  noyears <- length(surv)
  
  full.surv.model <- 1
  full.obs.model <- 1
  full.size.model <- 1
  full.repst.model <- 1
  full.fec.model <- 1
  
  tack.on.indiv <- NA
  tack.on.patch <- NA
  tack.on.year <- NA
  
  if (suite == "full") {
    if (noyears == 3) {
      main.body <- "sz1 + sz2 + fl1 + fl2 + sz1:sz2 + fl1:fl2 + sz1:fl1 + sz2:fl2 + sz1:fl2 + sz2:fl1"
    } else if (noyears == 2) {
      main.body <- "sz2 + fl2 + sz2:fl2"
    } else {
      stop("The surv vector must contain two entries for an ahistorical modelsearch and three entries for a historical modelsearch.", call. = FALSE)
    }
    
    if (!is.na(age)) {
      if (noyears == 3) {
        age.body <- " + ag + ag:sz1 + ag:sz2 + ag:fl1 + ag:fl2"
      } else if (noyears == 2) {
        age.body <- " + ag + ag:sz2 + ag:fl2"
      }
      main.body <- paste0(main.body, age.body)
    }
  } else if (suite == "main") {
    if (noyears == 3) {
      main.body <- "sz1 + sz2 + fl1 + fl2"
    } else if (noyears == 2) {
      main.body <- "sz2 + fl2"
    } else {
      stop("The surv vector must contain two entries for an ahistorical modelsearch and three entries for a historical modelsearch.", call. = FALSE)
    }
    
    if (!is.na(age)) {
      age.body <- " + ag"
      main.body <- paste0(main.body, age.body)
    }
  } else if (suite == "size") {
    if (noyears == 3) {
      main.body <- "sz1 + sz2 + sz1:sz2"
    } else if (noyears == 2) {
      main.body <- "sz2"
    } else {
      stop("The surv vector must contain two entries for an ahistorical modelsearch and three entries for a historical modelsearch.", call. = FALSE)
    }
    
    if (!is.na(age)) {
      if (noyears == 3) {
        age.body <- " + ag + ag:sz1 + ag:sz2"
      } else if (noyears == 2) {
        age.body <- " + ag + ag:sz2"
      }
      main.body <- paste0(main.body, age.body)
    }
  } else if (suite == "rep") {
    if (noyears == 3) {full.surv.model <- "al3 ~ fl1 + fl2 + fl1:fl2"}
    if (noyears == 3) {
      main.body <- "fl1 + fl2 + fl1:fl2"
    } else if (noyears == 2) {
      main.body <- "fl2"
    } else {
      stop("The surv vector must contain two entries for an ahistorical modelsearch and three entries for a historical modelsearch.", call. = FALSE)
    }
    
    if (!is.na(age)) {
      if (noyears == 3) {
        age.body <- " + ag + ag:fl1 + ag:fl2"
      } else if (noyears == 2) {
        age.body <- " + ag + ag:fl2"
      }
      main.body <- paste0(main.body, age.body)
    }
  } else if (suite == "cons") {
    main.body <- "1"
  } else {
    stop("Suite option not recognized.", .call = FALSE)
  }
  
  #Here we link up the response and independent terms
  if (length(surv) > 1) {full.surv.model <- paste0("al3 ~ ", main.body)}
  if (length(obs) > 1) {full.obs.model <- paste0("ob3 ~ ", main.body)}
  if (length(size) > 1) {full.size.model <- paste0("sz3 ~ ", main.body)}
  if (length(repst) > 1) {full.repst.model <- paste0("fl3 ~ ", main.body)}
  if (length(fec) > 1) {
    if (fectime == 2) {
      full.fec.model <- paste0("rp2 ~ ", main.body)
    } else if (fectime == 3) {
      full.fec.model <- paste0("rp3 ~ ", main.body)
    }
  }
  
  juv.surv.model <- 1
  juv.obs.model <- 1
  juv.size.model <- 1
  juv.repst.model <- 1
  
  if (!is.na(juvestimate)) {
    if (!any(is.element(stage, names(data)))) {stop("Names of stage variables do not match dataset.", call. = FALSE)}
    
    stage3col <- which(names(data) == stage[1])
    stage2col <- which(names(data) == stage[2])
    if (length(stage) == 3) {stage1col <- which(names(data) == stage[3])} else {stage1col <- 0}
    
    if (!is.element(juvestimate, unique(data[,stage2col]))) {stop("Stage input in juvestimate not recognized in stage2 variable in dataset", call. = FALSE)}
    
    if (full.surv.model != 1 & juvsize == FALSE) {
      juv.surv.model <- "al3 ~ 1"
    } else if (full.surv.model != 1 & juvsize == TRUE) {
      juv.surv.model <- "al3 ~ sz2"
    }
    
    if (full.obs.model != 1) {
      juv.obs.model <- "ob3 ~ 1"
    } else if (full.obs.model != 1 & juvsize == TRUE) {
      juv.obs.model <- "ob3 ~ sz2"
    }
    
    if (full.size.model != 1 & juvsize == FALSE) {
      juv.size.model <- "sz3 ~ 1"
    } else if (full.size.model != 1 & juvsize == TRUE) {
      juv.size.model <- "sz3 ~ sz2"
    }
    
    if (full.repst.model != 1) {
      juv.repst.model <- "fl3 ~ 1"
    } else if (full.repst.model != 1 & juvsize == TRUE) {
      juv.repst.model <- "fl3 ~ sz2"
    }
    
  }
  
  if (!is.element("surv", vitalrates)) {
    full.surv.model <- 1
    juv.surv.model <- 1
  }
  
  if (!is.element("obs", vitalrates)) {
    full.obs.model <- 1
    juv.obs.model <- 1
  }
  
  if (!is.element("size", vitalrates)) {
    full.size.model <- 1
    juv.size.model <- 1
  }
  
  if (!is.element("repst", vitalrates)) {
    full.repst.model <- 1
    juv.repst.model <- 1
  }
  
  if (!is.element("fec", vitalrates)) {
    full.fec.model <- 1
  }
  
  if (approach == "lme4") {
    if (!all(is.na(indiv))) {
      tack.on.indiv <- " + (1 | ind)";
      tack.on.indiv <- gsub("ind", indiv, tack.on.indiv);
      
      if (full.surv.model != 1) {full.surv.model <- paste0(full.surv.model, tack.on.indiv)}
      if (full.obs.model != 1) {full.obs.model <- paste0(full.obs.model, tack.on.indiv)}
      if (full.size.model != 1) {full.size.model <- paste0(full.size.model, tack.on.indiv)}
      if (full.repst.model != 1) {full.repst.model <- paste0(full.repst.model, tack.on.indiv)}
      if (full.fec.model != 1) {full.fec.model <- paste0(full.fec.model, tack.on.indiv)}
      
      if (juv.surv.model != 1) {juv.surv.model <- paste0(juv.surv.model, tack.on.indiv)}
      if (juv.obs.model != 1) {juv.obs.model <- paste0(juv.obs.model, tack.on.indiv)}
      if (juv.size.model != 1) {juv.size.model <- paste0(juv.size.model, tack.on.indiv)}
      if (juv.repst.model != 1) {juv.repst.model <- paste0(juv.repst.model, tack.on.indiv)}
    };
    
    if (year.as.random == TRUE) {
      tack.on.year <- " + (1 | yr)";
      
      if (full.surv.model != 1) {full.surv.model <- paste0(full.surv.model, tack.on.year)}
      if (full.obs.model != 1) {full.obs.model <- paste0(full.obs.model, tack.on.year)}
      if (full.size.model != 1) {full.size.model <- paste0(full.size.model, tack.on.year)}
      if (full.repst.model != 1) {full.repst.model <- paste0(full.repst.model, tack.on.year)}
      if (full.fec.model != 1) {full.fec.model <- paste0(full.fec.model, tack.on.year)}
      
      if (juv.surv.model != 1) {juv.surv.model <- paste0(juv.surv.model, tack.on.year)}
      if (juv.obs.model != 1) {juv.obs.model <- paste0(juv.obs.model, tack.on.year)}
      if (juv.size.model != 1) {juv.size.model <- paste0(juv.size.model, tack.on.year)}
      if (juv.repst.model != 1) {juv.repst.model <- paste0(juv.repst.model, tack.on.year)}
    };
    
    if (!all(is.na(patch))) {
      if (patch.as.random == TRUE) {
        tack.on.patch <- " + (1 | pat)";
      } else {
        tack.on.patch <- " + as.factor(pat)";
      }
      tack.on.patch <- gsub("pat", patch, tack.on.patch);
      
      if (full.surv.model != 1) {full.surv.model <- paste0(full.surv.model, tack.on.patch)}
      if (full.obs.model != 1) {full.obs.model <- paste0(full.obs.model, tack.on.patch)}
      if (full.size.model != 1) {full.size.model <- paste0(full.size.model, tack.on.patch)}
      if (full.repst.model != 1) {full.repst.model <- paste0(full.repst.model, tack.on.patch)}
      if (full.fec.model != 1) {full.fec.model <- paste0(full.fec.model, tack.on.patch)}
      
      if (juv.surv.model != 1) {juv.surv.model <- paste0(juv.surv.model, tack.on.patch)}
      if (juv.obs.model != 1) {juv.obs.model <- paste0(juv.obs.model, tack.on.patch)}
      if (juv.size.model != 1) {juv.size.model <- paste0(juv.size.model, tack.on.patch)}
      if (juv.repst.model != 1) {juv.repst.model <- paste0(juv.repst.model, tack.on.patch)}
    }
  } else if (approach == "glm") {
    if (year.as.random == FALSE) {
      tack.on.year <- " + as.factor(yr)";
      if (suite == "full") {
        tack.on.year <- paste0(tack.on.year, " + sz2:as.factor(yr) + fl2:as.factor(yr)");
        if (length(surv) == 3) {tack.on.year <- paste0(tack.on.year, " + sz1:as.factor(yr) + fl1:as.factor(yr)")}
        if (full.surv.model != 1) {full.surv.model <- paste0(full.surv.model, tack.on.year)}
        if (full.obs.model != 1) {full.obs.model <- paste0(full.obs.model, tack.on.year)}
        if (full.size.model != 1) {full.size.model <- paste0(full.size.model, tack.on.year)}
        if (full.repst.model != 1) {full.repst.model <- paste0(full.repst.model, tack.on.year)}
        if (full.fec.model != 1) {full.fec.model <- paste0(full.fec.model, tack.on.year)}
        
        if (juv.surv.model != 1) {juv.surv.model <- paste0(juv.surv.model, tack.on.year)}
        if (juv.obs.model != 1) {juv.obs.model <- paste0(juv.obs.model, tack.on.year)}
        if (juv.size.model != 1) {juv.size.model <- paste0(juv.size.model, tack.on.year)}
        if (juv.repst.model != 1) {juv.repst.model <- paste0(juv.repst.model, tack.on.year)}
      }
      if (suite != "full") {
        if (full.surv.model != 1) {full.surv.model <- paste0(full.surv.model, tack.on.year)}
        if (full.obs.model != 1) {full.obs.model <- paste0(full.obs.model, tack.on.year)}
        if (full.size.model != 1) {full.size.model <- paste0(full.size.model, tack.on.year)}
        if (full.repst.model != 1) {full.repst.model <- paste0(full.repst.model, tack.on.year)}
        if (full.fec.model != 1) {full.fec.model <- paste0(full.fec.model, tack.on.year)}
        
        if (juv.surv.model != 1) {juv.surv.model <- paste0(juv.surv.model, tack.on.year)}
        if (juv.obs.model != 1) {juv.obs.model <- paste0(juv.obs.model, tack.on.year)}
        if (juv.size.model != 1) {juv.size.model <- paste0(juv.size.model, tack.on.year)}
        if (juv.repst.model != 1) {juv.repst.model <- paste0(juv.repst.model, tack.on.year)}
      }
    }
    if (!all(is.na(patch)) & patch.as.random == FALSE) {
      tack.on.patch <- " + as.factor(pat)";
      tack.on.patch <- gsub("pat", patch, tack.on.patch);
      
      if (full.surv.model != 1) {full.surv.model <- paste0(full.surv.model, tack.on.patch)}
      if (full.obs.model != 1) {full.obs.model <- paste0(full.obs.model, tack.on.patch)}
      if (full.size.model != 1) {full.size.model <- paste0(full.size.model, tack.on.patch)}
      if (full.repst.model != 1) {full.repst.model <- paste0(full.repst.model, tack.on.patch)}
      if (full.fec.model != 1) {full.fec.model <- paste0(full.fec.model, tack.on.patch)}
      
      if (juv.surv.model != 1) {juv.surv.model <- paste0(juv.surv.model, tack.on.patch)}
      if (juv.obs.model != 1) {juv.obs.model <- paste0(juv.obs.model, tack.on.patch)}
      if (juv.size.model != 1) {juv.size.model <- paste0(juv.size.model, tack.on.patch)}
      if (juv.repst.model != 1) {juv.repst.model <- paste0(juv.repst.model, tack.on.patch)}
    }
  } else {
    stop("Approach option not recognized.", .call = FALSE)
  }
  
  #Now we create a dataframe for use to interpret linear models
  mainparams <- c("year2", "individ", "patch", "surv3", "obs3", "size3", "repst3", "fec3", "fec2", 
                  "size2", "size1", "repst2", "repst1", "age")
  if (historical == TRUE) {
    modelparams <- c(year, indiv, patch, surv[1], obs[1], size[1], repst[1], fec[1], fec[2], size[2], 
                     size[3], repst[2], repst[3], age)
  } else {
    modelparams <- c(year, indiv, patch, surv[1], obs[1], size[1], repst[1], fec[1], fec[2], size[2], 
                     NA, repst[2], NA, age)
  }
  paramnames <- cbind.data.frame(mainparams, modelparams)
  paramnames$mainparams <- as.character(paramnames$mainparams)
  paramnames$modelparams <- as.character(paramnames$modelparams)
  
  #We now replace factor codes for the appropriate factors from the function input
  if (full.surv.model != 1) {
    full.surv.model <- gsub("yr", year, full.surv.model)
    full.surv.model <- gsub("al3", surv[1], full.surv.model)
    full.surv.model <- gsub("sz2", size[2], full.surv.model)
    if (length(size) == 3) {
      full.surv.model <- gsub("sz1", size[3], full.surv.model)
    }
    full.surv.model <- gsub("fl2", repst[2], full.surv.model)
    if (length(repst) == 3) {
      full.surv.model <- gsub("fl1", repst[3], full.surv.model)
    }
    if (!is.na(age)) {
      full.surv.model <- gsub("ag", age, full.surv.model)
    }
  }
  
  if (full.obs.model != 1) {
    full.obs.model <- gsub("yr", year, full.obs.model)
    full.obs.model <- gsub("ob3", obs[1], full.obs.model)
    full.obs.model <- gsub("sz2", size[2], full.obs.model)
    if (length(size) == 3) {
      full.obs.model <- gsub("sz1", size[3], full.obs.model)
    }
    full.obs.model <- gsub("fl2", repst[2], full.obs.model)
    if (length(repst) == 3) {
      full.obs.model <- gsub("fl1", repst[3], full.obs.model)
    }
    if (!is.na(age)) {
      full.obs.model <- gsub("ag", age, full.obs.model)
    }
  }
  
  if (full.size.model != 1) {
    full.size.model <- gsub("yr", year, full.size.model)
    full.size.model <- gsub("sz3", size[1], full.size.model)
    full.size.model <- gsub("sz2", size[2], full.size.model)
    if (length(size) == 3) {
      full.size.model <- gsub("sz1", size[3], full.size.model)
    }
    full.size.model <- gsub("fl2", repst[2], full.size.model)
    if (length(repst) == 3) {
      full.size.model <- gsub("fl1", repst[3], full.size.model)
    }
    if (!is.na(age)) {
      full.size.model <- gsub("ag", age, full.size.model)
    }
  }
  
  if (full.repst.model != 1) {
    full.repst.model <- gsub("yr", year, full.repst.model)
    full.repst.model <- gsub("fl3", repst[1], full.repst.model)
    full.repst.model <- gsub("sz2", size[2], full.repst.model)
    if (length(size) == 3) {
      full.repst.model <- gsub("sz1", size[3], full.repst.model)
    }
    full.repst.model <- gsub("fl2", repst[2], full.repst.model)
    if (length(repst) == 3) {
      full.repst.model <- gsub("fl1", repst[3], full.repst.model)
    }
    if (!is.na(age)) {
      full.repst.model <- gsub("ag", age, full.repst.model)
    }
  }
  
  if (full.fec.model != 1) {
    full.fec.model <- gsub("yr", year, full.fec.model)
    if (fectime == 2) {
      full.fec.model <- gsub("rp2", fec[2], full.fec.model)
    } else if (fectime == 3) {
      full.fec.model <- gsub("rp3", fec[1], full.fec.model)
    }
    full.fec.model <- gsub("sz2", size[2], full.fec.model)
    if (length(size) == 3) {
      full.fec.model <- gsub("sz1", size[3], full.fec.model)
    }
    full.fec.model <- gsub("fl2", repst[2], full.fec.model)
    if (length(repst) == 3) {
      full.fec.model <- gsub("fl1", repst[3], full.fec.model)
    }
    if (!is.na(age)) {
      full.fec.model <- gsub("ag", age, full.fec.model)
    }
  }
  
  if (juv.surv.model != 1) {
    juv.surv.model <- gsub("yr", year, juv.surv.model)
    juv.surv.model <- gsub("al3", surv[1], juv.surv.model)
    if (juvsize == TRUE) {
      juv.surv.model <- gsub("sz2", size[2], juv.surv.model)
    }
  }
  
  if (juv.obs.model != 1) {
    juv.obs.model <- gsub("yr", year, juv.obs.model)
    juv.obs.model <- gsub("ob3", obs[1], juv.obs.model)
    if (juvsize == TRUE) {
      juv.obs.model <- gsub("sz2", size[2], juv.obs.model)
    }
  }
  
  if (juv.size.model != 1) {
    juv.size.model <- gsub("yr", year, juv.size.model)
    juv.size.model <- gsub("sz3", size[1], juv.size.model)
    if (juvsize == TRUE) {
      juv.size.model <- gsub("sz2", size[2], juv.size.model)
    }
  }
  
  if (juv.repst.model != 1) {
    juv.repst.model <- gsub("yr", year, juv.repst.model)
    juv.repst.model <- gsub("fl3", repst[1], juv.repst.model)
    if (juvsize == TRUE) {
      juv.repst.model <- gsub("sz2", size[2], juv.repst.model)
    }
  }
  
  #Now we need to create the input datasets
  if (!all(is.na(censor))) {
    if (length(censor) == 1) {
      data$censor2 <- data[, which(names(data) == censor[1])]
    } else {
      data$censor3 <- data[, which(names(data) == censor[1])]
      data$censor2 <- data[, which(names(data) == censor[2])]
      if (length(censor) > 2) {data$censor1 <- data[, which(names(data) == censor[3])]}
    }
  } else {
    data$censor2 <- 1
  }
  
  data <- subset(data, censor2 == 1)
  
  if (!all(is.na(censor))) {
    if (length(censor) > 1) {
      data <- subset(data, censor3 == 1)
      if (length(censor) > 2) {
        data <- subset(data, censor1 == 1)
      }
    }
  }
  
  if (!is.na(juvestimate)) {
    juvindivs <- which(data[,stage2col] == juvestimate)
    adultindivs <- setdiff(c(1:length(data[,stage2col])), juvindivs)
    
    juvsurv.data <- subset(data, data[,stage2col] == juvestimate & data[,which(names(data) == surv[2])] == 1)
    
    if (suite == "full" | suite == "main" | suite == "size") {
      if (any(is.na(juvsurv.data[, which(names(juvsurv.data) == size[2])]))) {
        warning("NAs in size variables may cause model selection to fail.")
      }
      
      if (historical == TRUE) {
        if (any(is.na(juvsurv.data[, which(names(juvsurv.data) == size[3])]))) {
          warning("NAs in size variables may cause model selection to fail.")
        }
      }
    }
    
    if (suite == "full" | suite == "main" | suite == "rep") {
      if (any(is.na(juvsurv.data[, which(names(juvsurv.data) == repst[2])]))) {
        warning("NAs in reproductive status variables may cause model selection to fail.")
      }
      
      if (historical == TRUE) {
        if (any(is.na(juvsurv.data[, which(names(juvsurv.data) == repst[3])]))) {
          warning("NAs in reproductive status variables may cause model selection to fail.")
        }
      }
    }
    
    if (is.element(0, juvsurv.data$matstatus3)) {
      warning("Modelsearch() assumes that all juveniles either die or transition to maturity within 1 year. Some individuals in this dataset appear to live longer as juveniles than assumptions allow.")
    }
    
    juvobs.data <- subset(juvsurv.data, juvsurv.data[, which(names(juvsurv.data) == surv[1])] == 1)
    if (full.obs.model != 1) {
      juvsize.data <- subset(juvobs.data, juvobs.data[, which(names(juvobs.data) == obs[1])] == 1)
      juvsize.data <- juvsize.data[which(!is.na(juvsize.data[, which(names(juvsize.data) == size[1])])),]
    } else {
      juvsize.data <- juvobs.data
      juvsize.data <- juvsize.data[which(!is.na(juvsize.data[, which(names(juvsize.data) == size[1])])),]
    }
    juvrepst.data <- juvsize.data
    
    if (dim(juvsurv.data)[1] < 100) {
      warning("Juvenile dataset is very small, and some models may fail given the size.")
    }
    
    data <- data[adultindivs,] #This line resets the main dataset to adults only
  }
  
  surv.data <- subset(data, data[,which(names(data) == surv[2])] == 1)
  
  if (suite == "full" | suite == "main" | suite == "size") {
    if (any(is.na(surv.data[, which(names(surv.data) == size[2])]))) {
      warning("NAs in size variables may cause model selection to fail.")
    }
    
    if (historical == TRUE) {
      if (any(is.na(surv.data[, which(names(surv.data) == size[3])]))) {
        warning("NAs in size variables may cause model selection to fail.")
      }
    }
  }
  
  if (suite == "full" | suite == "main" | suite == "rep") {
    if (any(is.na(surv.data[, which(names(surv.data) == repst[2])]))) {
      warning("NAs in reproductive status variables may cause model selection to fail.")
    }
    
    if (historical == TRUE) {
      if (any(is.na(surv.data[, which(names(surv.data) == repst[3])]))) {
        warning("NAs in reproductive status variables may cause model selection to fail.")
      }
    }
  }
  
  if (dim(surv.data)[1] < 100) {
    warning("Dataset is very small, and some models may fail given the size.")
  }
  
  if(any(!suppressWarnings(!is.na(as.numeric(as.character(surv.data[, which(names(surv.data) == size[1])])))))) {
    warning("Modelsearch(), flefko3(), and flefko2() are made to work with numeric size variables. Use of categorical variables may result in errors and unexpected behavior.")
  }
  
  obs.data <- subset(surv.data, surv.data[, which(names(surv.data) == surv[1])] == 1)
  if (full.obs.model != 1) {
    size.data <- subset(obs.data, obs.data[, which(names(obs.data) == obs[1])] == 1)
    size.data <- size.data[which(!is.na(size.data[, which(names(size.data) == size[1])])),]
  } else {
    size.data <- obs.data
    size.data <- size.data[which(!is.na(size.data[, which(names(size.data) == size[1])])),]
  }
  repst.data <- size.data
  if (full.repst.model != 1) {
    fec.data <- subset(surv.data, surv.data[, which(names(repst.data) == repst[2])] == 1)
    if (fectime == 2) {
      fec.data <- fec.data[which(!is.na(fec.data[, which(names(fec.data) == fec[2])])),]
    } else if (fectime == 3) {
      fec.data <- fec.data[which(!is.na(fec.data[, which(names(fec.data) == fec[1])])),]
    }
  } else {
    fec.data <- surv.data
    if (fectime == 2) {
      fec.data <- fec.data[which(!is.na(fec.data[, which(names(fec.data) == fec[2])])),]
    } else if (fectime == 3) {
      fec.data <- fec.data[which(!is.na(fec.data[, which(names(fec.data) == fec[1])])),]
    }
  }
  
  #Now we check for exceptions to size in the dataset
  if (sizedist == "poisson") {
    if (is.element(0, size.data[, which(names(size.data) == size[1])])) {
      stop("Cannot assume a Poisson distribution for size in time t+1 because that variable includes 0s.", call. = FALSE)
    }
    
    if (any(size.data[, which(names(size.data) == size[1])] != round(size.data[, which(names(size.data) == size[1])]))) {
      stop("Size variables must be composed only of integers greater than 0 for the Poisson distribution to be used.")
    }
    
    if (!is.na(juvestimate)) {
      if (is.element(0, juvsize.data[, which(names(juvsize.data) == size[1])])) {
        stop("Cannot assume a Poisson distribution for size in time t+1 because that variable includes 0s.", call. = FALSE)
      }
      
      if (any(juvsize.data[, which(names(juvsize.data) == size[1])] != round(size.data[, which(names(juvsize.data) == size[1])]))) {
        stop("Size variables must be composed only of integers greater than 0 for the Poisson distribution to be used.")
      }
    }
  } else if (sizedist == "negbin") {
    if (is.element(0, size.data[, which(names(size.data) == size[1])])) {
      stop("Cannot assume a negative binomial distribution for size in time t+1 because that variable includes 0s.", call. = FALSE)
    }
    
    if (any(size.data[, which(names(size.data) == size[1])] != round(size.data[, which(names(size.data) == size[1])]))) {
      stop("Size variables must be composed only of integers greater than 0 for the negative binomial distribution to be used.")
    }
    
    if (!is.na(juvestimate)) {
      if (is.element(0, juvsize.data[, which(names(juvsize.data) == size[1])])) {
        stop("Cannot assume a negative binomial distribution for size in time t+1 because that variable includes 0s.", call. = FALSE)
      }
      
      if (any(juvsize.data[, which(names(juvsize.data) == size[1])] != round(juvsize.data[, which(names(juvsize.data) == size[1])]))) {
        stop("Size variables must be composed only of integers greater than 0 for the negative binomial distribution to be used.")
      }
    }
  }
  
  if (fecdist == "poisson") {
    if (fectime == 2) {
      usedfec <- which(names(fec.data) == fec[2])
    } else if (fectime == 3) {
      usedfec <- which(names(fec.data) == fec[1])
    }
    if (is.element(0, fec.data[, usedfec])) {
      warning("WARNING: Fecundity in time t cannot be Poisson-distributed and include 0s. Will develop fecundity models excluding all 0s. Consider adding a reproductive status variable to absorb 0 values.\n")
      fec.data <- subset(fec.data, fec.data[, usedfec] > 0)
    }
    
    if (any(fec.data[, usedfec] != round(fec.data[, usedfec]))) {
      stop("Fecundity variables must be composed only of integers greater than 0 for the Poisson distribution to be used.")
    }
  } else if (fecdist == "negbin") {
    if (is.element(0, fec.data[, usedfec])) {
      warning("WARNING: Fecundity in time t cannot conform to the negative binomial distribution and include 0s. Will develop fecundity models excluding all 0s. Consider adding a reproductive status variable to absorb 0 values.\n")
      fec.data <- subset(fec.data, fec.data[, usedfec] > 0)
    }
    
    if (any(fec.data[, usedfec] != round(fec.data[, usedfec]))) {
      stop("Fecundity variables must be composed only of integers greater than 0 for the negative binomial distribution to be used.")
    }
  }
  
  #Now we run the modeling exercises
  if (full.surv.model == 1) {surv.global.model <- 1}
  if (full.obs.model == 1) {obs.global.model <- 1}
  if (full.size.model == 1) {size.global.model <- 1}
  if (full.repst.model == 1) {repst.global.model <- 1}
  if (full.fec.model == 1) {fec.global.model <- 1}
  
  if (juv.surv.model == 1) {juv.surv.global.model <- 1}
  if (juv.obs.model == 1) {juv.obs.global.model <- 1}
  if (juv.size.model == 1) {juv.size.global.model <- 1}
  if (juv.repst.model == 1) {juv.repst.global.model <- 1}
  
  surv.table <- NA
  obs.table <- NA
  size.table <- NA
  repst.table <- NA
  fec.table <- NA
  
  juvsurv.table <- NA
  juvobs.table <- NA
  juvsize.table <- NA
  juvrepst.table <- NA
  
  surv.bf <- NA
  obs.bf <- NA
  size.bf <- NA
  repst.bf <- NA
  fec.bf <- NA
  
  juvsurv.bf <- NA
  juvobs.bf <- NA
  juvsize.bf <- NA
  juvrepst.bf <- NA
  
  #A few more corrections to the model structure, used in running the global models
  correction.sz1sz2 <- gsub("sz1", size[3], " + sz1:sz2", fixed = TRUE)
  correction.sz1sz2 <- gsub("sz2", size[2], correction.sz1sz2, fixed = TRUE)
  correction.fl1fl2 <- gsub("fl1", repst[3], " + fl1:fl2", fixed = TRUE)
  correction.fl1fl2 <- gsub("fl2", repst[2], correction.fl1fl2, fixed = TRUE)
  correction.sz1fl1 <- gsub("sz1", size[3], " + sz1:fl1", fixed = TRUE)
  correction.sz1fl1 <- gsub("fl1", repst[3], correction.sz1fl1, fixed = TRUE)
  correction.sz2fl2 <- gsub("sz2", size[2], " + sz2:fl2", fixed = TRUE)
  correction.sz2fl2 <- gsub("fl2", repst[2], correction.sz2fl2, fixed = TRUE)
  correction.sz1fl2 <- gsub("sz1", size[3], " + sz1:fl2", fixed = TRUE)
  correction.sz1fl2 <- gsub("fl2", repst[2], correction.sz1fl2, fixed = TRUE)
  correction.sz2fl1 <- gsub("sz2", size[2], " + sz2:fl1", fixed = TRUE)
  correction.sz2fl1 <- gsub("fl1", repst[3], correction.sz2fl1, fixed = TRUE)
  
  correction.indiv <- gsub("individ", indiv, tack.on.indiv, fixed = TRUE)
  correction.year <- gsub("yr", year, tack.on.year, fixed = TRUE)
  correction.patch <- tack.on.patch #This is already correctly formatted
  
  #Here we run the global models
  if (approach == "lme4") {
    if (full.surv.model != 1) {
      if (is.element(0, surv.data$alive3) & is.element(1, surv.data$alive3)) {
        writeLines("\nDeveloping global model of survival probability...\n"); 
        surv.global.model <- try(lme4::glmer(formula = full.surv.model, data = surv.data, family = "binomial"),
                                 silent = TRUE)
        
        if (class(surv.global.model) == "try-error") {
          nox.surv.model <- full.surv.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.surv.model <- gsub(correction.sz1sz2, "", nox.surv.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.surv.model <- gsub(correction.fl1fl2, "", nox.surv.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.surv.model <- gsub(correction.sz1fl1, "", nox.surv.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.surv.model <- gsub(correction.sz2fl2, "", nox.surv.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.surv.model <- gsub(correction.sz1fl2, "", nox.surv.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.surv.model <- gsub(correction.sz2fl1, "", nox.surv.model, fixed = TRUE)
          }
          
          if (nox.surv.model != full.surv.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            surv.global.model <- try(lme4::glmer(formula = nox.surv.model, data = surv.data, family = "binomial"),
                                     silent = TRUE)
          }
          
          nopat.surv.model <- nox.surv.model
          if (!is.na(correction.patch)) {
            nopat.surv.model <- gsub(correction.patch, "", nopat.surv.model, fixed = TRUE)
          }
          
          if (class(surv.global.model) == "try-error") {
            
            if (nox.surv.model != nopat.surv.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              surv.global.model <- try(lme4::glmer(formula = nopat.surv.model, data = surv.data, family = "binomial"),
                                       silent = TRUE)
            }
          }
          
          noyr.surv.model <- nopat.surv.model
          if (!is.na(correction.year)) {
            noyr.surv.model <- gsub(correction.year, "", noyr.surv.model, fixed = TRUE)
          }
          
          if (class(surv.global.model) == "try-error") {
            
            if (noyr.surv.model != nopat.surv.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              surv.global.model <- try(lme4::glmer(formula = noyr.surv.model, data = surv.data, family = "binomial"),
                                       silent = TRUE)
            }
          }
          
          noind.surv.model <- noyr.surv.model
          if (!is.na(correction.indiv)) {
            noind.surv.model <- gsub(correction.indiv, "", noind.surv.model, fixed = TRUE)
          }
          
          if (class(surv.global.model) == "try-error") {
            
            if (noind.surv.model != noyr.surv.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              surv.global.model <- try(lme4::glmer(formula = noind.surv.model, data = surv.data, family = "binomial"),
                                       silent = TRUE)
            }
          }
          
          if (class(surv.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for survival.")
          }
        }
      } else if (!is.element(0, surv.data$alive3)) {
        writeLines("\nSurvival response is constant so will not model it.")
        full.surv.model <- 1
        surv.global.model <- 1
      } else {
        writeLines("\nSurvival response is constant so will not model it.")
        full.surv.model <- 1
        surv.global.model <- 0
      }
    }
    
    if (full.obs.model != 1) {
      if (is.element(0, obs.data$obsstatus3) & is.element(1, obs.data$obsstatus3)) {
        writeLines("\nDeveloping global model of observation probability...\n"); 
        obs.global.model <- try(lme4::glmer(formula = full.obs.model, data = obs.data, family = "binomial"),
                                silent = TRUE)
        
        if (class(obs.global.model) == "try-error") {
          nox.obs.model <- full.obs.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.obs.model <- gsub(correction.sz1sz2, "", nox.obs.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.obs.model <- gsub(correction.fl1fl2, "", nox.obs.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.obs.model <- gsub(correction.sz1fl1, "", nox.obs.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.obs.model <- gsub(correction.sz2fl2, "", nox.obs.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.obs.model <- gsub(correction.sz1fl2, "", nox.obs.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.obs.model <- gsub(correction.sz2fl1, "", nox.obs.model, fixed = TRUE)
          }
          
          if (nox.obs.model != full.obs.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            obs.global.model <- try(lme4::glmer(formula = nox.obs.model, data = obs.data, family = "binomial"),
                                    silent = TRUE)
          }
          
          nopat.obs.model <- nox.obs.model
          if (!is.na(correction.patch)) {
            nopat.obs.model <- gsub(correction.patch, "", nopat.obs.model, fixed = TRUE)
          }
          
          if (class(obs.global.model) == "try-error") {
            
            if (nox.obs.model != nopat.obs.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              obs.global.model <- try(lme4::glmer(formula = nopat.obs.model, data = obs.data, family = "binomial"),
                                      silent = TRUE)
            }
          }
          
          noyr.obs.model <- nopat.obs.model
          if (!is.na(correction.year)) {
            noyr.obs.model <- gsub(correction.year, "", noyr.obs.model, fixed = TRUE)
          }
          
          if (class(obs.global.model) == "try-error") {
            
            if (noyr.obs.model != nopat.obs.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              obs.global.model <- try(lme4::glmer(formula = noyr.obs.model, data = obs.data, family = "binomial"),
                                      silent = TRUE)
            }
          }
          
          noind.obs.model <- noyr.obs.model
          if (!is.na(correction.indiv)) {
            noind.obs.model <- gsub(correction.indiv, "", noyr.obs.model, fixed = TRUE)
          }
          
          if (class(obs.global.model) == "try-error") {
            
            if (noind.obs.model != noyr.obs.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              obs.global.model <- try(lme4::glmer(formula = noind.obs.model, data = obs.data, family = "binomial"),
                                      silent = TRUE)
            }
          }
          
          if (class(obs.global.model) == "try-error") {
            writeLines("Could not properly estimate a global model for observation status.")
          }
        }
      } else if (!is.element(0, obs.data$obsstatus3)) {
        writeLines("\nObservation response is constant so will not model it.")
        full.obs.model <- 1
        obs.global.model <- 1
      } else {
        writeLines("\nObservation response is constant so will not model it.")
        full.obs.model <- 1
        obs.global.model <- 0
      }
    }
    
    if (full.size.model != 1) {
      if (sizedist == "gaussian") {
        writeLines("\nDeveloping global model of size (Gaussian)...\n");
        size.global.model <- try(lme4::lmer(formula = full.size.model, data = size.data),
                                 silent = TRUE)
        
        if (class(size.global.model) == "try-error") {
          nox.size.model <- full.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.size.model <- gsub(correction.sz1sz2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.size.model <- gsub(correction.fl1fl2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.size.model <- gsub(correction.sz1fl1, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.size.model <- gsub(correction.sz2fl2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.size.model <- gsub(correction.sz1fl2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.size.model <- gsub(correction.sz2fl1, "", nox.size.model, fixed = TRUE)
          }
          
          if (nox.size.model != full.size.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            size.global.model <- try(lme4::lmer(formula = nox.size.model, data = size.data), silent = TRUE)
          }
          
          nopat.size.model <- nox.size.model
          if (!is.na(correction.patch)) {
            nopat.size.model <- gsub(correction.patch, "", nopat.size.model, fixed = TRUE)
          }
          
          if (class(size.global.model) == "try-error") {
            
            if (nox.size.model != nopat.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              size.global.model <- try(lme4::lmer(formula = nopat.size.model, data = size.data),
                                       silent = TRUE)
            }
          }
          
          noyr.size.model <- nopat.size.model
          if (!is.na(correction.year)) {
            noyr.size.model <- gsub(correction.year, "", noyr.size.model, fixed = TRUE)
          }
          
          if (class(size.global.model) == "try-error") {
            
            if (noyr.size.model != nopat.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              size.global.model <- try(lme4::lmer(formula = noyr.size.model, data = size.data),
                                       silent = TRUE)
            }
          }
          
          noind.size.model <- noyr.size.model
          if (!is.na(correction.indiv)) {
            noind.size.model <- gsub(correction.indiv, "", noind.size.model, fixed = TRUE)
          }
          
          if (class(size.global.model) == "try-error") {
            
            if (noind.size.model != noyr.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              size.global.model <- try(lme4::lmer(formula = noind.size.model, data = size.data),
                                       silent = TRUE)
            }
          }
          
          if (class(size.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for size.")
          }
        }
      } else if (sizedist == "poisson") {
        writeLines("\nDeveloping global model of size (Poisson)...\n");
        size.global.model <- try(lme4::glmer(formula = full.size.model, data = size.data, family = "poisson"),
                                 silent = TRUE)
        
        if (class(size.global.model) == "try-error") {
          nox.size.model <- full.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.size.model <- gsub(correction.sz1sz2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.size.model <- gsub(correction.fl1fl2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.size.model <- gsub(correction.sz1fl1, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.size.model <- gsub(correction.sz2fl2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.size.model <- gsub(correction.sz1fl2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.size.model <- gsub(correction.sz2fl1, "", nox.size.model, fixed = TRUE)
          }
          
          if (nox.size.model != full.size.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            size.global.model <- try(lme4::glmer(formula = nox.size.model, data = size.data, family = "poisson"),
                                     silent = TRUE)
          }
          
          nopat.size.model <- nox.size.model
          if (!is.na(correction.patch)) {
            nopat.size.model <- gsub(correction.patch, "", nopat.size.model, fixed = TRUE)
          }
          
          if (class(size.global.model) == "try-error") {
            
            if (nox.size.model != nopat.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              size.global.model <- try(lme4::glmer(formula = nopat.size.model, data = size.data, family = "poisson"),
                                       silent = TRUE)
            }
          }
          
          noyr.size.model <- nopat.size.model
          if (!is.na(correction.year)) {
            noyr.size.model <- gsub(correction.year, "", noyr.size.model, fixed = TRUE)
          }
          
          if (class(size.global.model) == "try-error") {
            
            if (noyr.size.model != nopat.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              size.global.model <- try(lme4::glmer(formula = noyr.size.model, data = size.data, family = "poisson"),
                                       silent = TRUE)
            }
          }
          
          noind.size.model <- noyr.size.model
          if (!is.na(correction.indiv)) {
            noind.size.model <- gsub(correction.indiv, "", noind.size.model, fixed = TRUE)
          }
          
          if (class(size.global.model) == "try-error") {
            
            if (noind.size.model != noyr.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              size.global.model <- try(lme4::glmer(formula = noind.size.model, data = size.data, family = "poisson"),
                                       silent = TRUE)
            }
          }
          
          if (class(size.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for size.")
          }
        }
      } else if (sizedist == "negbin") {
        writeLines("\nDeveloping global model of size (negative binomial)...\n");
        size.global.model <- try(glmmTMB::glmmTMB(formula = as.formula(full.size.model), data = size.data,
                                                  ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
        
        if (class(size.global.model) == "try-error") {
          nox.size.model <- full.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.size.model <- gsub(correction.sz1sz2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.size.model <- gsub(correction.fl1fl2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.size.model <- gsub(correction.sz1fl1, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.size.model <- gsub(correction.sz2fl2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.size.model <- gsub(correction.sz1fl2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.size.model <- gsub(correction.sz2fl1, "", nox.size.model, fixed = TRUE)
          }
          
          if (nox.size.model != full.size.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            size.global.model <- try(glmmTMB::glmmTMB(formula = as.formula(nox.size.model), data = size.data,
                                                      ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
          }
          
          nopat.size.model <- nox.size.model
          if (!is.na(correction.patch)) {
            nopat.size.model <- gsub(correction.patch, "", nopat.size.model, fixed = TRUE)
          }
          
          if (class(size.global.model) == "try-error") {
            
            if (nox.size.model != nopat.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              size.global.model <- try(glmmTMB::glmmTMB(formula = as.formula(nopat.size.model), data = size.data,
                                                        ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
            }
          }
          
          noyr.size.model <- nopat.size.model
          if (!is.na(correction.year)) {
            noyr.size.model <- gsub(correction.year, "", noyr.size.model, fixed = TRUE)
          }
          
          if (class(size.global.model) == "try-error") {
            
            if (noyr.size.model != nopat.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              size.global.model <- try(glmmTMB::glmmTMB(formula = as.formula(noyr.size.model), data = size.data,
                                                        ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
            }
          }
          
          noind.size.model <- noyr.size.model
          if (!is.na(correction.indiv)) {
            noind.size.model <- gsub(correction.indiv, "", noind.size.model, fixed = TRUE)
          }
          
          if (class(size.global.model) == "try-error") {
            
            if (noind.size.model != noyr.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              size.global.model <- try(glmmTMB::glmmTMB(formula = as.formula(noind.size.model), data = size.data,
                                                        ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
            }
          }
          
          if (class(size.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for size.")
          }
        }
      }
    }
    
    if (full.repst.model != 1) {
      if (is.element(0, repst.data$repstatus3) & is.element(1, repst.data$repstatus3)) {
        writeLines("\nDeveloping global model of the probability of reproduction...\n"); 
        repst.global.model <- try(lme4::glmer(formula = full.repst.model, data = repst.data, family = "binomial"),
                                  silent = TRUE)
        
        if (class(repst.global.model) == "try-error") {
          nox.repst.model <- full.repst.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.repst.model <- gsub(correction.sz1sz2, "", nox.repst.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.repst.model <- gsub(correction.fl1fl2, "", nox.repst.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.repst.model <- gsub(correction.sz1fl1, "", nox.repst.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.repst.model <- gsub(correction.sz2fl2, "", nox.repst.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.repst.model <- gsub(correction.sz1fl2, "", nox.repst.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.repst.model <- gsub(correction.sz2fl1, "", nox.repst.model, fixed = TRUE)
          }
          
          if (nox.repst.model != full.repst.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            repst.global.model <- try(lme4::glmer(formula = nox.repst.model, data = repst.data, family = "binomial"),
                                      silent = TRUE)
          }
          
          nopat.repst.model <- nox.repst.model
          if (!is.na(correction.patch)) {
            nopat.repst.model <- gsub(correction.patch, "", nopat.repst.model, fixed = TRUE)
          }
          
          if (class(repst.global.model) == "try-error") {
            
            if (nox.repst.model != nopat.repst.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              repst.global.model <- try(lme4::glmer(formula = nopat.repst.model, data = repst.data, family = "binomial"),
                                        silent = TRUE)
            }
          }
          
          noyr.repst.model <- nopat.repst.model
          if (!is.na(correction.year)) {
            noyr.repst.model <- gsub(correction.year, "", noyr.repst.model, fixed = TRUE)
          }
          
          if (class(repst.global.model) == "try-error") {
            
            if (noyr.repst.model != nopat.repst.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              repst.global.model <- try(lme4::glmer(formula = noyr.repst.model, data = repst.data, family = "binomial"),
                                        silent = TRUE)
            }
          }
          
          noind.repst.model <- noyr.repst.model
          if (!is.na(correction.indiv)) {
            noind.repst.model <- gsub(correction.indiv, "", noind.repst.model, fixed = TRUE)
          }
          
          if (class(repst.global.model) == "try-error") {
            
            if (noind.repst.model != noyr.repst.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              repst.global.model <- try(lme4::glmer(formula = noind.repst.model, data = repst.data, family = "binomial"),
                                        silent = TRUE)
            }
          }
          
          if (class(repst.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for reproduction status.")
          }
        }
      } else if (!is.element(0, repst.data$repstatus3)) {
        writeLines("\nReproductive status response is constant so will not model it.")
        full.repst.model <- 1
        repst.global.model <- 1
      } else {
        writeLines("\nReproductive status response is constant so will not model it.")
        full.repst.model <- 1
        repst.global.model <- 0
      }
    }
    
    if (full.fec.model != 1) {
      if (fecdist == "gaussian") {
        writeLines("\nDeveloping global model of fecundity (Gaussian)...\n");
        fec.global.model <- try(lme4::lmer(formula = full.fec.model, data = fec.data), silent = TRUE)
        
        if (class(fec.global.model) == "try-error") {
          nox.fec.model <- full.fec.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.fec.model <- gsub(correction.sz1sz2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.fec.model <- gsub(correction.fl1fl2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.fec.model <- gsub(correction.sz1fl1, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.fec.model <- gsub(correction.sz2fl2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.fec.model <- gsub(correction.sz1fl2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.fec.model <- gsub(correction.sz2fl1, "", nox.fec.model, fixed = TRUE)
          }
          
          if (nox.fec.model != full.fec.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            fec.global.model <- try(lme4::lmer(formula = nox.fec.model, data = fec.data),
                                    silent = TRUE)
          }
          
          nopat.fec.model <- nox.fec.model
          if (!is.na(correction.patch)) {
            nopat.fec.model <- gsub(correction.patch, "", nopat.fec.model, fixed = TRUE)
          }
          
          if (class(fec.global.model) == "try-error") {
            
            if (nox.fec.model != nopat.fec.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              fec.global.model <- try(lme4::lmer(formula = nopat.fec.model, data = fec.data),
                                      silent = TRUE)
            }
          }
          
          noyr.fec.model <- nopat.fec.model
          if (!is.na(correction.patch)) {
            noyr.fec.model <- gsub(correction.year, "", noyr.fec.model, fixed = TRUE)
          }
          
          if (class(fec.global.model) == "try-error") {
            
            if (noyr.fec.model != nopat.fec.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              fec.global.model <- try(lme4::lmer(formula = noyr.fec.model, data = fec.data),
                                      silent = TRUE)
            }
          }
          
          noind.fec.model <- noyr.fec.model
          if (!is.na(correction.indiv)) {
            noind.fec.model <- gsub(correction.indiv, "", noind.fec.model, fixed = TRUE)
          }
          
          if (class(fec.global.model) == "try-error") {
            
            if (noind.fec.model != noyr.fec.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              fec.global.model <- try(lme4::lmer(formula = noind.fec.model, data = fec.data),
                                      silent = TRUE)
            }
          }
          
          if (class(fec.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for fecundity.")
          }
        }
      } else if (fecdist == "poisson") {
        writeLines("\nDeveloping global model of fecundity (Poisson)...\n");
        fec.global.model <- try(lme4::glmer(formula = full.fec.model, data = fec.data, family = "poisson"),
                                silent = TRUE)
        
        if (class(fec.global.model) == "try-error") {
          nox.fec.model <- full.fec.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.fec.model <- gsub(correction.sz1sz2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.fec.model <- gsub(correction.fl1fl2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.fec.model <- gsub(correction.sz1fl1, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.fec.model <- gsub(correction.sz2fl2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.fec.model <- gsub(correction.sz1fl2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.fec.model <- gsub(correction.sz2fl1, "", nox.fec.model, fixed = TRUE)
          }
          
          if (nox.fec.model != full.fec.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            fec.global.model <- try(lme4::glmer(formula = nox.fec.model, data = fec.data, family = "poisson"),
                                    silent = TRUE)
          }
          
          nopat.fec.model <- nox.fec.model
          if (!is.na(correction.patch)) {
            nopat.fec.model <- gsub(correction.patch, "", nopat.fec.model, fixed = TRUE)
          }
          
          if (class(fec.global.model) == "try-error") {
            
            if (nox.fec.model != nopat.fec.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              fec.global.model <- try(lme4::glmer(formula = nopat.fec.model, data = fec.data, family = "poisson"),
                                      silent = TRUE)
            }
          }
          
          noyr.fec.model <- nopat.fec.model
          if (!is.na(correction.year)) {
            noyr.fec.model <- gsub(correction.year, "", noyr.fec.model, fixed = TRUE)
          }
          
          if (class(fec.global.model) == "try-error") {
            
            if (noyr.fec.model != nopat.fec.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              fec.global.model <- try(lme4::glmer(formula = noyr.fec.model, data = fec.data, family = "poisson"),
                                      silent = TRUE)
            }
          }
          
          noind.fec.model <- noyr.fec.model
          if (!is.na(correction.indiv)) {
            noind.fec.model <- gsub(correction.indiv, "", noind.fec.model, fixed = TRUE)
          }
          
          if (class(fec.global.model) == "try-error") {
            
            if (noind.fec.model != noyr.fec.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              fec.global.model <- try(lme4::glmer(formula = noind.fec.model, data = fec.data, family = "poisson"),
                                      silent = TRUE)
            }
          }
          
          if (class(fec.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for fecundity.")
          }
        }
      } else if (fecdist == "negbin") {
        writeLines("\nDeveloping global model of fecundity (negative binomial)...\n");
        fec.global.model <- try(glmmTMB::glmmTMB(formula = as.formula(full.fec.model), data = fec.data,
                                                 ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
        
        if (class(fec.global.model) == "try-error") {
          nox.fec.model <- full.fec.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.fec.model <- gsub(correction.sz1sz2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.fec.model <- gsub(correction.fl1fl2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.fec.model <- gsub(correction.sz1fl1, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.fec.model <- gsub(correction.sz2fl2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.fec.model <- gsub(correction.sz1fl2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.fec.model <- gsub(correction.sz2fl1, "", nox.fec.model, fixed = TRUE)
          }
          
          if (nox.fec.model != full.fec.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            fec.global.model <- try(glmmTMB::glmmTMB(formula = as.formula(nox.fec.model), data = fec.data,
                                                     ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
          }
          
          nopat.fec.model <- nox.fec.model
          if (!is.na(correction.patch)) {
            nopat.fec.model <- gsub(correction.patch, "", nopat.fec.model, fixed = TRUE)
          }
          
          if (class(fec.global.model) == "try-error") {
            
            if (nox.fec.model != nopat.fec.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              fec.global.model <- try(glmmTMB::glmmTMB(formula = as.formula(nopat.fec.model), data = fec.data,
                                                       ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
            }
          }
          
          noyr.fec.model <- nopat.fec.model
          if (!is.na(correction.year)) {
            noyr.fec.model <- gsub(correction.year, "", noyr.fec.model, fixed = TRUE)
          }
          
          if (class(fec.global.model) == "try-error") {
            
            if (noyr.fec.model != nopat.fec.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              fec.global.model <- try(glmmTMB::glmmTMB(formula = as.formula(noyr.fec.model), data = fec.data,
                                                       ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
            }
          }
          
          noind.fec.model <- noyr.fec.model
          if (!is.na(correction.indiv)) {
            noind.fec.model <- gsub(correction.indiv, "", noind.fec.model, fixed = TRUE)
          }
          
          if (class(fec.global.model) == "try-error") {
            
            if (noind.fec.model != noyr.fec.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              fec.global.model <- try(glmmTMB::glmmTMB(formula = as.formula(noind.fec.model), data = fec.data,
                                                       ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
            }
          }
          
          if (class(fec.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for fecundity.")
          }
        }
      }
    }
    
    if (juv.surv.model != 1) {
      if (is.element(0, juvsurv.data$alive3) & is.element(1, juvsurv.data$alive3)) {
        writeLines("\nDeveloping global model of juvenile survival probability...\n"); 
        juv.surv.global.model <- try(lme4::glmer(formula = juv.surv.model, data = juvsurv.data, 
                                                 family = "binomial"), silent = TRUE)
        
        if (class(juv.surv.global.model) == "try-error") {
          nox.juv.surv.model <- juv.surv.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.surv.model <- gsub(correction.sz1sz2, "", nox.juv.surv.model, fixed = TRUE)
          }
          
          if (nox.juv.surv.model != juv.surv.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            juv.surv.global.model <- try(lme4::glmer(formula = nox.juv.surv.model, data = juvsurv.data, 
                                                     family = "binomial"), silent = TRUE)
          }
          
          nopat.juv.surv.model <- nox.juv.surv.model
          if (!is.na(correction.patch)) {
            nopat.juv.surv.model <- gsub(correction.patch, "", nopat.juv.surv.model, fixed = TRUE)
          }
          
          if (class(juv.surv.global.model) == "try-error") {
            
            if (nox.juv.surv.model != nopat.juv.surv.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              juv.surv.global.model <- try(lme4::glmer(formula = nopat.juv.surv.model, data = juvsurv.data, 
                                                       family = "binomial"), silent = TRUE)
            }
          }
          
          noyr.juv.surv.model <- nopat.juv.surv.model
          if (!is.na(correction.year)) {
            noyr.juv.surv.model <- gsub(correction.year, "", noyr.juv.surv.model, fixed = TRUE)
          }
          
          if (class(juv.surv.global.model) == "try-error") {
            
            if (noyr.juv.surv.model != nopat.juv.surv.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              juv.surv.global.model <- try(lme4::glmer(formula = noyr.juv.surv.model, data = juvsurv.data, 
                                                       family = "binomial"), silent = TRUE)
            }
          }
          
          noind.juv.surv.model <- noyr.juv.surv.model
          if (!is.na(correction.indiv)) {
            noind.juv.surv.model <- gsub(correction.indiv, "", noind.juv.surv.model, fixed = TRUE)
          }
          
          if (class(juv.surv.global.model) == "try-error") {
            
            if (noind.juv.surv.model != noyr.juv.surv.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              juv.surv.global.model <- try(lme4::glmer(formula = noind.juv.surv.model, data = juvsurv.data, 
                                                       family = "binomial"), silent = TRUE)
            }
          }
          
          if (class(juv.surv.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for juvenile survival.")
          }
        }
      } else if (!is.element(0, juvsurv.data$alive3)) {
        writeLines("\nJuvenile survival response is constant so will not model it.")
        juv.surv.model <- 1
        juv.surv.global.model <- 1
      } else {
        writeLines("\nJuvenile survival response is constant so will not model it.")
        juv.surv.model <- 1
        juv.surv.global.model <- 0
      }
    }
    
    if (juv.obs.model != 1) {
      if (is.element(0, juvobs.data$obsstatus3) & is.element(1, juvobs.data$obsstatus3)) {
        writeLines("\nDeveloping global model of juvenile observation probability...\n"); 
        juv.obs.global.model <- try(lme4::glmer(formula = juv.obs.model, data = juvobs.data, 
                                                family = "binomial"), silent = TRUE)
        
        if (class(juv.obs.global.model) == "try-error") {
          nox.juv.obs.model <- juv.obs.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.obs.model <- gsub(correction.sz1sz2, "", nox.juv.obs.model, fixed = TRUE)
          }
          
          if (nox.juv.obs.model != juv.obs.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            juv.obs.global.model <- try(lme4::glmer(formula = nox.juv.obs.model, data = juvobs.data, 
                                                    family = "binomial"), silent = TRUE)
          }
          
          nopat.juv.obs.model <- nox.juv.obs.model
          if (!is.na(correction.patch)) {
            nopat.juv.obs.model <- gsub(correction.patch, "", nopat.juv.obs.model, fixed = TRUE)
          }
          
          if (class(juv.obs.global.model) == "try-error") {
            
            if (nox.juv.obs.model != nopat.juv.obs.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              juv.obs.global.model <- try(lme4::glmer(formula = nopat.juv.obs.model, data = juvobs.data, 
                                                      family = "binomial"), silent = TRUE)
            }
          }
          
          noyr.juv.obs.model <- nopat.juv.obs.model
          if (!is.na(correction.year)) {
            noyr.juv.obs.model <- gsub(correction.year, "", noyr.juv.obs.model, fixed = TRUE)
          }
          
          if (class(juv.obs.global.model) == "try-error") {
            
            if (noyr.juv.obs.model != nopat.juv.obs.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              juv.obs.global.model <- try(lme4::glmer(formula = noyr.juv.obs.model, data = juvobs.data, 
                                                      family = "binomial"), silent = TRUE)
            }
          }
          
          noind.juv.obs.model <- noyr.juv.obs.model
          if (!is.na(correction.indiv)) {
            noind.juv.obs.model <- gsub(correction.indiv, "", noind.juv.obs.model, fixed = TRUE)
          }
          
          if (class(juv.obs.global.model) == "try-error") {
            
            if (noind.juv.obs.model != noyr.juv.obs.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              juv.obs.global.model <- try(lme4::glmer(formula = noind.juv.obs.model, data = juvobs.data, 
                                                      family = "binomial"), silent = TRUE)
            }
          }
          
          if (class(juv.obs.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for juvenile observation status.")
          }
        }
      } else if (!is.element(0, juvobs.data$obsstatus3)) {
        writeLines("\nJuvenile observation response is constant so will not model it.")
        juv.obs.model <- 1
        juv.obs.global.model <- 1
      } else {
        writeLines("\nJuvenile observation response is constant so will not model it.")
        juv.obs.model <- 1
        juv.obs.global.model <- 0
      }
    }
    
    if (juv.size.model != 1) {
      if (sizedist == "gaussian") {
        writeLines("\nDeveloping global model of juvenile size (Gaussian)...\n");
        juv.size.global.model <- try(lme4::lmer(formula = juv.size.model, data = juvsize.data),
                                     silent = TRUE)
        
        if (class(juv.size.global.model) == "try-error") {
          nox.juv.size.model <- juv.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.size.model <- gsub(correction.sz1sz2, "", nox.juv.size.model, fixed = TRUE)
          }
          
          if (nox.juv.size.model != juv.size.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            juv.size.global.model <- try(lme4::lmer(formula = nox.juv.size.model, data = juvsize.data), 
                                         silent = TRUE)
          }
          
          nopat.juv.size.model <- nox.juv.size.model
          if (!is.na(correction.patch)) {
            nopat.juv.size.model <- gsub(correction.patch, "", nopat.juv.size.model, fixed = TRUE)
          }
          
          if (class(juv.size.global.model) == "try-error") {
            
            if (nox.juv.size.model != nopat.juv.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              juv.size.global.model <- try(lme4::lmer(formula = nopat.juv.size.model, data = juvsize.data),
                                           silent = TRUE)
            }
          }
          
          noyr.juv.size.model <- nopat.juv.size.model
          if (!is.na(correction.year)) {
            noyr.juv.size.model <- gsub(correction.year, "", noyr.juv.size.model, fixed = TRUE)
          }
          
          if (class(juv.size.global.model) == "try-error") {
            
            if (noyr.juv.size.model != nopat.juv.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              juv.size.global.model <- try(lme4::lmer(formula = noyr.juv.size.model, data = juvsize.data),
                                           silent = TRUE)
            }
          }
          
          noind.juv.size.model <- noyr.juv.size.model
          if (!is.na(correction.indiv)) {
            noind.juv.size.model <- gsub(correction.indiv, "", noind.juv.size.model, fixed = TRUE)
          }
          
          if (class(juv.size.global.model) == "try-error") {
            
            if (noind.juv.size.model != noyr.juv.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              juv.size.global.model <- try(lme4::lmer(formula = noind.juv.size.model, data = juvsize.data),
                                           silent = TRUE)
            }
          }
          
          if (class(juv.size.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for juvenile size.")
          }
        }
      } else if (sizedist == "poisson") {
        writeLines("\nDeveloping global model of juvenile size (Poisson)...\n");
        juv.size.global.model <- try(lme4::glmer(formula = juv.size.model, data = juvsize.data, 
                                                 family = "poisson"), silent = TRUE)
        
        if (class(juv.size.global.model) == "try-error") {
          nox.juv.size.model <- juv.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.size.model <- gsub(correction.sz1sz2, "", nox.juv.size.model, fixed = TRUE)
          }
          
          if (nox.juv.size.model != juv.size.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            juv.size.global.model <- try(lme4::glmer(formula = nox.juv.size.model, data = juvsize.data, 
                                                     family = "poisson"), silent = TRUE)
          }
          
          nopat.juv.size.model <- nox.juv.size.model
          if (!is.na(correction.patch)) {
            nopat.juv.size.model <- gsub(correction.patch, "", nopat.juv.size.model, fixed = TRUE)
          }
          
          if (class(juv.size.global.model) == "try-error") {
            
            if (nox.juv.size.model != nopat.juv.size.model) {
              writeLines("Global model estimation difficulties. Attempting a global model without a patch term.")
              juv.size.global.model <- try(lme4::glmer(formula = nopat.juv.size.model, data = juvsize.data, 
                                                       family = "poisson"), silent = TRUE)
            }
          }
          
          noyr.juv.size.model <- nopat.juv.size.model
          if (!is.na(correction.year)) {
            noyr.juv.size.model <- gsub(correction.year, "", noyr.juv.size.model, fixed = TRUE)
          }
          
          if (class(juv.size.global.model) == "try-error") {
            
            if (noyr.juv.size.model != nopat.juv.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              juv.size.global.model <- try(lme4::glmer(formula = noyr.juv.size.model, data = juvsize.data, 
                                                       family = "poisson"), silent = TRUE)
            }
          }
          
          noind.juv.size.model <- noyr.juv.size.model
          if (!is.na(correction.indiv)) {
            noind.juv.size.model <- gsub(correction.indiv, "", noind.juv.size.model, fixed = TRUE)
          }
          
          if (class(juv.size.global.model) == "try-error") {
            
            if (noind.juv.size.model != noyr.juv.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              juv.size.global.model <- try(lme4::glmer(formula = noind.juv.size.model, data = juvsize.data, 
                                                       family = "poisson"), silent = TRUE)
            }
          }
          
          if (class(juv.size.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for juvenile size.")
          }
        }
      } else if (sizedist == "negbin") {
        writeLines("\nDeveloping global model of juvenile size (negative binomial)...\n");
        juv.size.global.model <- try(glmmTMB::glmmTMB(formula = as.formula(juv.size.model), data = juvsize.data, 
                                                      ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
        
        if (class(juv.size.global.model) == "try-error") {
          nox.juv.size.model <- juv.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.size.model <- gsub(correction.sz1sz2, "", nox.juv.size.model, fixed = TRUE)
          }
          
          if (nox.juv.size.model != juv.size.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            juv.size.global.model <- try(glmmTMB::glmmTMB(formula = as.formula(nox.juv.size.model), 
                                                          data = juvsize.data, ziformula=~0, family = glmmTMB::nbinom2), 
                                         silent = TRUE)
          }
          
          nopat.juv.size.model <- nox.juv.size.model
          if (!is.na(correction.patch)) {
            nopat.juv.size.model <- gsub(correction.patch, "", nopat.juv.size.model, fixed = TRUE)
          }
          
          if (class(juv.size.global.model) == "try-error") {
            
            if (nox.juv.size.model != nopat.juv.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              juv.size.global.model <- try(glmmTMB::glmmTMB(formula = as.formula(nopat.juv.size.model), 
                                                            data = juvsize.data, ziformula=~0, family = glmmTMB::nbinom2), 
                                           silent = TRUE)
            }
          }
          
          noyr.juv.size.model <- nopat.juv.size.model
          if (!is.na(correction.year)) {
            noyr.juv.size.model <- gsub(correction.year, "", noyr.juv.size.model, fixed = TRUE)
          }
          
          if (class(juv.size.global.model) == "try-error") {
            
            if (noyr.juv.size.model != nopat.juv.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              juv.size.global.model <- try(glmmTMB::glmmTMB(formula = as.formula(noyr.juv.size.model), 
                                                            data = juvsize.data, ziformula=~0, family = glmmTMB::nbinom2), 
                                           silent = TRUE)
            }
          }
          
          noind.juv.size.model <- noyr.juv.size.model
          if (!is.na(correction.indiv)) {
            noind.juv.size.model <- gsub(correction.indiv, "", noind.juv.size.model, fixed = TRUE)
          }
          
          if (class(juv.size.global.model) == "try-error") {
            
            if (noind.juv.size.model != noyr.juv.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              juv.size.global.model <- try(glmmTMB::glmmTMB(formula = as.formula(noind.juv.size.model), 
                                                            data = juvsize.data, ziformula=~0, family = glmmTMB::nbinom2), 
                                           silent = TRUE)
            }
          }
          
          if (class(juv.size.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for juvenile size.")
          }
        }
      }
    }
    
    if (juv.repst.model != 1) {
      if (is.element(0, juvrepst.data$repstatus3) & is.element(1, juvrepst.data$repstatus3)) {
        writeLines("\nDeveloping global model of juvenile reproduction probability...\n"); 
        juv.repst.global.model <- try(lme4::glmer(formula = juv.repst.model, data = juvrepst.data, 
                                                  family = "binomial"), silent = TRUE)
        
        if (class(juv.repst.global.model) == "try-error") {
          nox.juv.repst.model <- juv.repst.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.repst.model <- gsub(correction.sz1sz2, "", nox.juv.repst.model, fixed = TRUE)
          }
          
          if (nox.juv.repst.model != juv.repst.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            juv.repst.global.model <- try(lme4::glmer(formula = nox.juv.repst.model, data = juvrepst.data, 
                                                      family = "binomial"), silent = TRUE)
          }
          
          nopat.juv.repst.model <- nox.juv.repst.model
          if (!is.na(correction.patch)) {
            nopat.juv.repst.model <- gsub(correction.patch, "", nopat.juv.repst.model, fixed = TRUE)
          }
          
          if (class(juv.repst.global.model) == "try-error") {
            
            if (nox.juv.repst.model != nopat.juv.repst.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              juv.repst.global.model <- try(lme4::glmer(formula = nopat.juv.repst.model, data = juvrepst.data,
                                                        family = "binomial"), silent = TRUE)
            }
          }
          
          noyr.juv.repst.model <- nopat.juv.repst.model
          if (!is.na(correction.year)) {
            noyr.juv.repst.model <- gsub(correction.year, "", noyr.juv.repst.model, fixed = TRUE)
          }
          
          if (class(juv.repst.global.model) == "try-error") {
            
            if (noyr.juv.repst.model != nopat.juv.repst.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              juv.repst.global.model <- try(lme4::glmer(formula = noyr.juv.repst.model, data = juvrepst.data,
                                                        family = "binomial"), silent = TRUE)
            }
          }
          
          noind.juv.repst.model <- noyr.juv.repst.model
          if (!is.na(correction.indiv)) {
            noind.juv.repst.model <- gsub(correction.indiv, "", noind.juv.repst.model, fixed = TRUE)
          }
          
          if (class(juv.repst.global.model) == "try-error") {
            
            if (noind.juv.repst.model != noyr.juv.repst.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              juv.repst.global.model <- try(lme4::glmer(formula = noind.juv.repst.model, data = juvrepst.data,
                                                        family = "binomial"), silent = TRUE)
            }
          }
          
          if (class(juv.repst.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for juvenile reproduction status.")
          }
        }
      } else if (!is.element(0, juvrepst.data$repstatus3)) {
        writeLines("\nJuvenile reproduction response is constant so will not model it.")
        juv.repst.model <- 1
        juv.repst.global.model <- 1
      } else {
        writeLines("\nJuvenile reproduction response is constant so will not model it.")
        juv.repst.model <- 1
        juv.repst.global.model <- 0
      }
    }
    
    writeLines("\nAll global models developed.\n")
  }
  
  if (approach == "glm") {
    if (full.surv.model != 1) {
      if (is.element(0, surv.data$alive3) & is.element(1, surv.data$alive3)) {
        writeLines("\nDeveloping global model of survival probability...\n"); 
        surv.global.model <- try(glm(formula = full.surv.model, data = surv.data, family = "binomial"),
                                 silent = TRUE)
        
        if (class(surv.global.model) == "try-error") {
          nox.surv.model <- full.surv.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.surv.model <- gsub(correction.sz1sz2, "", nox.surv.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.surv.model <- gsub(correction.fl1fl2, "", nox.surv.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.surv.model <- gsub(correction.sz1fl1, "", nox.surv.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.surv.model <- gsub(correction.sz2fl2, "", nox.surv.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.surv.model <- gsub(correction.sz1fl2, "", nox.surv.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.surv.model <- gsub(correction.sz2fl1, "", nox.surv.model, fixed = TRUE)
          }
          
          if (nox.surv.model != full.surv.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            surv.global.model <- try(glm(formula = nox.surv.model, data = surv.data, family = "binomial"),
                                     silent = TRUE)
          }
          
          nopat.surv.model <- nox.surv.model
          if (!is.na(correction.patch)) {
            nopat.surv.model <- gsub(correction.patch, "", nopat.surv.model, fixed = TRUE)
          }
          
          if (class(surv.global.model) == "try-error") {
            
            if (nox.surv.model != nopat.surv.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              surv.global.model <- try(glm(formula = nopat.surv.model, data = surv.data, family = "binomial"),
                                       silent = TRUE)
            }
          }
          
          noyr.surv.model <- nopat.surv.model
          if (!is.na(correction.year)) {
            noyr.surv.model <- gsub(correction.year, "", noyr.surv.model, fixed = TRUE)
          }
          
          if (class(surv.global.model) == "try-error") {
            
            if (noyr.surv.model != nopat.surv.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              surv.global.model <- try(glm(formula = noyr.surv.model, data = surv.data, family = "binomial"),
                                       silent = TRUE)
            }
          }
          
          noind.surv.model <- noyr.surv.model
          if (!is.na(correction.indiv)) {
            noind.surv.model <- gsub(correction.indiv, "", noind.surv.model, fixed = TRUE)
          }
          
          if (class(surv.global.model) == "try-error") {
            
            if (noind.surv.model != noyr.surv.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              surv.global.model <- try(glm(formula = noind.surv.model, data = surv.data, family = "binomial"),
                                       silent = TRUE)
            }
          }
          
          if (class(surv.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for survival.")
          }
        }
      } else if (!is.element(0, surv.data$alive3)) {
        writeLines("\nSurvival response is constant so will not model it.")
        full.surv.model <- 1
        surv.global.model <- 1
      } else {
        writeLines("\nSurvival response is constant so will not model it.")
        full.surv.model <- 1
        surv.global.model <- 0
      }
    }
    
    if (full.obs.model != 1) {
      if (is.element(0, obs.data$obsstatus3) & is.element(1, obs.data$obsstatus3)) {
        writeLines("\nDeveloping global model of observation probability...\n"); 
        obs.global.model <- try(glm(formula = full.obs.model, data = obs.data, family = "binomial"),
                                silent = TRUE)
        
        if (class(obs.global.model) == "try-error") {
          nox.obs.model <- full.obs.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.obs.model <- gsub(correction.sz1sz2, "", nox.obs.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.obs.model <- gsub(correction.fl1fl2, "", nox.obs.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.obs.model <- gsub(correction.sz1fl1, "", nox.obs.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.obs.model <- gsub(correction.sz2fl2, "", nox.obs.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.obs.model <- gsub(correction.sz1fl2, "", nox.obs.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.obs.model <- gsub(correction.sz2fl1, "", nox.obs.model, fixed = TRUE)
          }
          
          if (nox.obs.model != full.obs.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            obs.global.model <- try(glm(formula = nox.obs.model, data = obs.data, family = "binomial"),
                                    silent = TRUE)
          }
          
          nopat.obs.model <- nox.obs.model
          if (!is.na(correction.patch)) {
            nopat.obs.model <- gsub(correction.patch, "", nopat.obs.model, fixed = TRUE)
          }
          
          if (class(obs.global.model) == "try-error") {
            
            if (nox.obs.model != nopat.obs.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              obs.global.model <- try(glm(formula = nopat.obs.model, data = obs.data, family = "binomial"),
                                      silent = TRUE)
            }
          }
          
          noyr.obs.model <- nopat.obs.model
          if (!is.na(correction.year)) {
            noyr.obs.model <- gsub(correction.year, "", noyr.obs.model, fixed = TRUE)
          }
          
          if (class(obs.global.model) == "try-error") {
            
            if (noyr.obs.model != nopat.obs.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              obs.global.model <- try(glm(formula = noyr.obs.model, data = obs.data, family = "binomial"),
                                      silent = TRUE)
            }
          }
          
          noind.obs.model <- noyr.obs.model
          if (!is.na(correction.indiv)) {
            noind.obs.model <- gsub(correction.indiv, "", noind.obs.model, fixed = TRUE)
          }
          
          if (class(obs.global.model) == "try-error") {
            
            if (noind.obs.model != noyr.obs.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              obs.global.model <- try(glm(formula = noind.obs.model, data = obs.data, family = "binomial"),
                                      silent = TRUE)
            }
          }
          
          if (class(obs.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for observation status.")
          }
        }
      } else if (!is.element(0, obs.data$obsstatus3)) {
        writeLines("\nObservation response is constant so will not model it.")
        full.obs.model <- 1
        obs.global.model <- 1
      } else {
        writeLines("\nObservation response is constant so will not model it.")
        full.obs.model <- 1
        obs.global.model <- 0
      }
    }
    
    if (full.size.model != 1) {
      if (sizedist == "gaussian") {
        writeLines("\nDeveloping global model of size (Gaussian)...\n");
        size.global.model <- try(lm(formula = full.size.model, data = size.data), silent = TRUE)
        
        if (class(size.global.model) == "try-error") {
          nox.size.model <- full.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.size.model <- gsub(correction.sz1sz2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.size.model <- gsub(correction.fl1fl2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.size.model <- gsub(correction.sz1fl1, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.size.model <- gsub(correction.sz2fl2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.size.model <- gsub(correction.sz1fl2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.size.model <- gsub(correction.sz2fl1, "", nox.size.model, fixed = TRUE)
          }
          
          if (nox.size.model != full.size.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            size.global.model <- try(lm(formula = nox.size.model, data = size.data), silent = TRUE)
          }
          
          nopat.size.model <- nox.size.model
          if (!is.na(correction.patch)) {
            nopat.size.model <- gsub(correction.patch, "", nopat.size.model, fixed = TRUE)
          }
          
          if (class(size.global.model) == "try-error") {
            
            if (nox.size.model != nopat.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              size.global.model <- try(lm(formula = nopat.size.model, data = size.data), silent = TRUE)
            }
          }
          
          noyr.size.model <- nopat.size.model
          if (!is.na(correction.year)) {
            noyr.size.model <- gsub(correction.year, "", noyr.size.model, fixed = TRUE)
          }
          
          if (class(size.global.model) == "try-error") {
            
            if (noyr.size.model != nopat.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              size.global.model <- try(lm(formula = noyr.size.model, data = size.data), silent = TRUE)
            }
          }
          
          noind.size.model <- noyr.size.model
          if (!is.na(correction.indiv)) {
            noind.size.model <- gsub(correction.indiv, "", noind.size.model, fixed = TRUE)
          }
          
          if (class(size.global.model) == "try-error") {
            
            if (noind.size.model != noyr.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              size.global.model <- try(lm(formula = noind.size.model, data = size.data), silent = TRUE)
            }
          }
          
          if (class(size.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for size.")
          }
        }
      } else if (sizedist == "poisson") {
        writeLines("\nDeveloping global model of size (Poisson)...\n");
        size.global.model <- try(glm(formula = full.size.model, data = size.data, family = "poisson"),
                                 silent = TRUE)
        
        if (class(size.global.model) == "try-error") {
          nox.size.model <- full.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.size.model <- gsub(correction.sz1sz2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.size.model <- gsub(correction.fl1fl2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.size.model <- gsub(correction.sz1fl1, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.size.model <- gsub(correction.sz2fl2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.size.model <- gsub(correction.sz1fl2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.size.model <- gsub(correction.sz2fl1, "", nox.size.model, fixed = TRUE)
          }
          
          if (nox.size.model != full.size.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            size.global.model <- try(glm(formula = nox.size.model, data = size.data, family = "poisson"),
                                     silent = TRUE)
          }
          
          nopat.size.model <- nox.size.model
          if (!is.na(correction.patch)) {
            nopat.size.model <- gsub(correction.patch, "", nopat.size.model, fixed = TRUE)
          }
          
          if (class(size.global.model) == "try-error") {
            
            if (nox.size.model != nopat.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              size.global.model <- try(glm(formula = nopat.size.model, data = size.data, family = "poisson"),
                                       silent = TRUE)
            }
          }
          
          noyr.size.model <- nopat.size.model
          if (!is.na(correction.year)) {
            noyr.size.model <- gsub(correction.year, "", noyr.size.model, fixed = TRUE)
          }
          
          if (class(size.global.model) == "try-error") {
            
            if (noyr.size.model != nopat.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              size.global.model <- try(glm(formula = noyr.size.model, data = size.data, family = "poisson"),
                                       silent = TRUE)
            }
          }
          
          noind.size.model <- noyr.size.model
          if (!is.na(correction.indiv)) {
            noind.size.model <- gsub(correction.indiv, "", noind.size.model, fixed = TRUE)
          }
          
          if (class(size.global.model) == "try-error") {
            
            if (noind.size.model != noyr.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              size.global.model <- try(glm(formula = noind.size.model, data = size.data, family = "poisson"),
                                       silent = TRUE)
            }
          }
          
          if (class(size.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for size.")
          }
        }
      } else if (sizedist == "negbin") {
        writeLines("\nDeveloping global model of size (negative binomial)...\n");
        size.global.model <- try(glm.nb(formula = full.size.model, data = size.data), silent = TRUE)
        
        if (class(size.global.model) == "try-error") {
          nox.size.model <- full.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.size.model <- gsub(correction.sz1sz2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.size.model <- gsub(correction.fl1fl2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.size.model <- gsub(correction.sz1fl1, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.size.model <- gsub(correction.sz2fl2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.size.model <- gsub(correction.sz1fl2, "", nox.size.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.size.model <- gsub(correction.sz2fl1, "", nox.size.model, fixed = TRUE)
          }
          
          if (nox.size.model != full.size.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            size.global.model <- try(glm.nb(formula = nox.size.model, data = size.data),
                                     silent = TRUE)
          }
          
          nopat.size.model <- nox.size.model
          if (!is.na(correction.patch)) {
            nopat.size.model <- gsub(correction.patch, "", nopat.size.model, fixed = TRUE)
          }
          
          if (class(size.global.model) == "try-error") {
            
            if (nox.size.model != nopat.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              size.global.model <- try(glm.nb(formula = nopat.size.model, data = size.data),
                                       silent = TRUE)
            }
          }
          
          noyr.size.model <- nopat.size.model
          if (!is.na(correction.year)) {
            noyr.size.model <- gsub(correction.year, "", noyr.size.model, fixed = TRUE)
          }
          
          if (class(size.global.model) == "try-error") {
            
            if (noyr.size.model != nopat.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              size.global.model <- try(glm.nb(formula = noyr.size.model, data = size.data),
                                       silent = TRUE)
            }
          }
          
          noind.size.model <- noyr.size.model
          if (!is.na(correction.indiv)) {
            noind.size.model <- gsub(correction.indiv, "", noind.size.model, fixed = TRUE)
          }
          
          if (class(size.global.model) == "try-error") {
            
            if (noind.size.model != noyr.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              size.global.model <- try(glm.nb(formula = noind.size.model, data = size.data),
                                       silent = TRUE)
            }
          }
          
          if (class(size.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for size.")
          }
        }
      }
    }
    
    if (full.repst.model != 1) {
      if (is.element(0, repst.data$repstatus3) & is.element(1, repst.data$repstatus3)) {
        writeLines("\nDeveloping global model of the probability of reproduction...\n"); 
        repst.global.model <- try(glm(formula = full.repst.model, data = repst.data, family = "binomial"),
                                  silent = TRUE)
        
        if (class(repst.global.model) == "try-error") {
          nox.repst.model <- full.repst.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.repst.model <- gsub(correction.sz1sz2, "", nox.repst.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.repst.model <- gsub(correction.fl1fl2, "", nox.repst.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.repst.model <- gsub(correction.sz1fl1, "", nox.repst.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.repst.model <- gsub(correction.sz2fl2, "", nox.repst.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.repst.model <- gsub(correction.sz1fl2, "", nox.repst.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.repst.model <- gsub(correction.sz2fl1, "", nox.repst.model, fixed = TRUE)
          }
          
          if (nox.repst.model != full.repst.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            repst.global.model <- try(glm(formula = nox.repst.model, data = repst.data, family = "binomial"),
                                      silent = TRUE)
          }
          
          nopat.repst.model <- nox.repst.model
          if (!is.na(correction.patch)) {
            nopat.repst.model <- gsub(correction.patch, "", nopat.repst.model, fixed = TRUE)
          }
          
          if (class(repst.global.model) == "try-error") {
            
            if (nox.repst.model != nopat.repst.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              repst.global.model <- try(glm(formula = nopat.repst.model, data = repst.data, family = "binomial"),
                                        silent = TRUE)
            }
          }
          
          noyr.repst.model <- nopat.repst.model
          if (!is.na(correction.year)) {
            noyr.repst.model <- gsub(correction.year, "", noyr.repst.model, fixed = TRUE)
          }
          
          if (class(repst.global.model) == "try-error") {
            
            if (noyr.repst.model != nopat.repst.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              repst.global.model <- try(glm(formula = noyr.repst.model, data = repst.data, family = "binomial"),
                                        silent = TRUE)
            }
          }
          
          noind.repst.model <- noyr.repst.model
          if (!is.na(correction.indiv)) {
            noind.repst.model <- gsub(correction.indiv, "", noind.repst.model, fixed = TRUE)
          }
          
          if (class(repst.global.model) == "try-error") {
            
            if (noind.repst.model != noyr.repst.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              repst.global.model <- try(glm(formula = noind.repst.model, data = repst.data, family = "binomial"),
                                        silent = TRUE)
            }
          }
          
          if (class(repst.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for reproduction status.")
          }
        }
      } else if (!is.element(0, repst.data$repstatus3)) {
        writeLines("\nReproductive status response is constant so will not model it.")
        full.repst.model <- 1
        repst.global.model <- 1
      } else {
        writeLines("\nReproductive status response is constant so will not model it.")
        full.repst.model <- 1
        repst.global.model <- 0
      }
    }
    
    if (full.fec.model != 1) {
      if (fecdist == "gaussian") {
        writeLines("\nDeveloping global model of fecundity (Gaussian)...\n");
        fec.global.model <- try(lm(formula = full.fec.model, data = fec.data), silent = TRUE)
        
        if (class(fec.global.model) == "try-error") {
          nox.fec.model <- full.fec.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.fec.model <- gsub(correction.sz1sz2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.fec.model <- gsub(correction.fl1fl2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.fec.model <- gsub(correction.sz1fl1, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.fec.model <- gsub(correction.sz2fl2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.fec.model <- gsub(correction.sz1fl2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.fec.model <- gsub(correction.sz2fl1, "", nox.fec.model, fixed = TRUE)
          }
          
          if (nox.fec.model != full.fec.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            fec.global.model <- try(lm(formula = nox.fec.model, data = fec.data), silent = TRUE)
          }
          
          nopat.fec.model <- nox.fec.model
          if (!is.na(correction.patch)) {
            nopat.fec.model <- gsub(correction.patch, "", nopat.fec.model, fixed = TRUE)
          }
          
          if (class(fec.global.model) == "try-error") {
            
            if (nox.fec.model != nopat.fec.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              fec.global.model <- try(lm(formula = nopat.fec.model, data = fec.data), silent = TRUE)
            }
          }
          
          noyr.fec.model <- nopat.fec.model
          if (!is.na(correction.year)) {
            noyr.fec.model <- gsub(correction.year, "", noyr.fec.model, fixed = TRUE)
          }
          
          if (class(fec.global.model) == "try-error") {
            
            if (noyr.fec.model != nopat.fec.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              fec.global.model <- try(lm(formula = noyr.fec.model, data = fec.data), silent = TRUE)
            }
          }
          
          noind.fec.model <- noyr.fec.model
          if (!is.na(correction.indiv)) {
            noind.fec.model <- gsub(correction.indiv, "", noind.fec.model, fixed = TRUE)
          }
          
          if (class(fec.global.model) == "try-error") {
            
            if (noind.fec.model != noyr.fec.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              fec.global.model <- try(lm(formula = noind.fec.model, data = fec.data), silent = TRUE)
            }
          }
          
          if (class(fec.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for fecundity.")
          }
        }
      } else if (fecdist == "poisson") {
        writeLines("\nDeveloping global model of fecundity (Poisson)...\n");
        fec.global.model <- try(glm(formula = full.fec.model, data = fec.data, family = "poisson"),
                                silent = TRUE)
        
        if (class(fec.global.model) == "try-error") {
          nox.fec.model <- full.fec.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.fec.model <- gsub(correction.sz1sz2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.fec.model <- gsub(correction.fl1fl2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.fec.model <- gsub(correction.sz1fl1, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.fec.model <- gsub(correction.sz2fl2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.fec.model <- gsub(correction.sz1fl2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.fec.model <- gsub(correction.sz2fl1, "", nox.fec.model, fixed = TRUE)
          }
          
          if (nox.fec.model != full.fec.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            fec.global.model <- try(glm(formula = nox.fec.model, data = fec.data, family = "poisson"),
                                    silent = TRUE)
          }
          
          nopat.fec.model <- nox.fec.model
          if (!is.na(correction.patch)) {
            nopat.fec.model <- gsub(correction.patch, "", nopat.fec.model, fixed = TRUE)
          }
          
          if (class(fec.global.model) == "try-error") {
            
            if (nox.fec.model != nopat.fec.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              fec.global.model <- try(glm(formula = nopat.fec.model, data = fec.data, family = "poisson"),
                                      silent = TRUE)
            }
          }
          
          noyr.fec.model <- nopat.fec.model
          if (!is.na(correction.year)) {
            noyr.fec.model <- gsub(correction.year, "", noyr.fec.model, fixed = TRUE)
          }
          
          if (class(fec.global.model) == "try-error") {
            
            if (noyr.fec.model != nopat.fec.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              fec.global.model <- try(glm(formula = noyr.fec.model, data = fec.data, family = "poisson"),
                                      silent = TRUE)
            }
          }
          
          noind.fec.model <- noyr.fec.model
          if (!is.na(correction.indiv)) {
            noind.fec.model <- gsub(correction.indiv, "", noind.fec.model, fixed = TRUE)
          }
          
          if (class(fec.global.model) == "try-error") {
            
            if (noind.fec.model != noyr.fec.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              fec.global.model <- try(glm(formula = noind.fec.model, data = fec.data, family = "poisson"),
                                      silent = TRUE)
            }
          }
          
          if (class(fec.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for fecundity.")
          }
        }
      } else if (fecdist == "negbin") {
        writeLines("\nDeveloping global model of fecundity (negative binomial)...\n");
        fec.global.model <- try(glm.nb(formula = full.fec.model, data = fec.data), silent = TRUE)
        
        if (class(fec.global.model) == "try-error") {
          nox.fec.model <- full.fec.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.fec.model <- gsub(correction.sz1sz2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.fl1fl2)) {
            nox.fec.model <- gsub(correction.fl1fl2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl1)) {
            nox.fec.model <- gsub(correction.sz1fl1, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl2)) {
            nox.fec.model <- gsub(correction.sz2fl2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz1fl2)) {
            nox.fec.model <- gsub(correction.sz1fl2, "", nox.fec.model, fixed = TRUE)
          }
          if (!is.na(correction.sz2fl1)) {
            nox.fec.model <- gsub(correction.sz2fl1, "", nox.fec.model, fixed = TRUE)
          }
          
          if (nox.fec.model != full.fec.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            fec.global.model <- try(glm.nb(formula = nox.fec.model, data = fec.data),
                                    silent = TRUE)
          }
          
          nopat.fec.model <- nox.fec.model
          if (!is.na(correction.patch)) {
            nopat.fec.model <- gsub(correction.patch, "", nopat.fec.model, fixed = TRUE)
          }
          
          if (class(fec.global.model) == "try-error") {
            
            if (nox.fec.model != nopat.fec.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              fec.global.model <- try(glm.nb(formula = nopat.fec.model, data = fec.data),
                                      silent = TRUE)
            }
          }
          
          noyr.fec.model <- nopat.fec.model
          if (!is.na(correction.year)) {
            noyr.fec.model <- gsub(correction.year, "", noyr.fec.model, fixed = TRUE)
          }
          
          if (class(fec.global.model) == "try-error") {
            
            if (noyr.fec.model != nopat.fec.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              fec.global.model <- try(glm.nb(formula = noyr.fec.model, data = fec.data),
                                      silent = TRUE)
            }
          }
          
          noind.fec.model <- noyr.fec.model
          if (!is.na(correction.indiv)) {
            noind.fec.model <- gsub(correction.indiv, "", noind.fec.model, fixed = TRUE)
          }
          
          if (class(fec.global.model) == "try-error") {
            
            if (noind.fec.model != noyr.fec.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              fec.global.model <- try(glm.nb(formula = noind.fec.model, data = fec.data),
                                      silent = TRUE)
            }
          }
          
          if (class(fec.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for fecundity.")
          }
        }
      }
    }
    
    if (juv.surv.model != 1) {
      if (is.element(0, juvsurv.data$alive3) & is.element(1, juvsurv.data$alive3)) {
        writeLines("\nDeveloping global model of juvenile survival probability...\n"); 
        juv.surv.global.model <- try(glm(formula = juv.surv.model, data = juvsurv.data, family = "binomial"),
                                     silent = TRUE)
        
        if (class(juv.surv.global.model) == "try-error") {
          nox.juv.surv.model <- juv.surv.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.surv.model <- gsub(correction.sz1sz2, "", nox.juv.surv.model, fixed = TRUE)
          }
          
          if (nox.juv.surv.model != juv.surv.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            juv.surv.global.model <- try(glm(formula = nox.juv.surv.model, data = juvsurv.data, family = "binomial"),
                                         silent = TRUE)
          }
          
          nopat.juv.surv.model <- nox.juv.surv.model
          if (!is.na(correction.patch)) {
            nopat.juv.surv.model <- gsub(correction.patch, "", nopat.juv.surv.model, fixed = TRUE)
          }
          
          if (class(juv.surv.global.model) == "try-error") {
            
            if (nox.juv.surv.model != nopat.juv.surv.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              juv.surv.global.model <- try(glm(formula = nopat.juv.surv.model, data = juvsurv.data, family = "binomial"),
                                           silent = TRUE)
            }
          }
          
          noyr.juv.surv.model <- nopat.juv.surv.model
          if (!is.na(correction.year)) {
            noyr.juv.surv.model <- gsub(correction.year, "", noyr.juv.surv.model, fixed = TRUE)
          }
          
          if (class(juv.surv.global.model) == "try-error") {
            
            if (noyr.juv.surv.model != nopat.juv.surv.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              juv.surv.global.model <- try(glm(formula = noyr.juv.surv.model, data = juvsurv.data, family = "binomial"),
                                           silent = TRUE)
            }
          }
          
          noind.juv.surv.model <- noyr.juv.surv.model
          if (!is.na(correction.indiv)) {
            noind.juv.surv.model <- gsub(correction.indiv, "", noind.juv.surv.model, fixed = TRUE)
          }
          
          if (class(juv.surv.global.model) == "try-error") {
            
            if (noind.juv.surv.model != noyr.juv.surv.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              juv.surv.global.model <- try(glm(formula = noind.juv.surv.model, data = juvsurv.data, family = "binomial"),
                                           silent = TRUE)
            }
          }
          
          if (class(juv.surv.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for juvenile survival.")
          }
        }
      } else if (!is.element(0, juvsurv.data$alive3)) {
        writeLines("\nJuvenile survival response is constant so will not model it.")
        juv.surv.model <- 1
        juv.surv.global.model <- 1
      } else {
        writeLines("\nJuvenile survival response is constant so will not model it.")
        juv.surv.model <- 1
        juv.surv.global.model <- 0
      }
    }
    
    if (juv.obs.model != 1) {
      if (is.element(0, juvobs.data$obsstatus3) & is.element(1, juvobs.data$obsstatus3)) {
        writeLines("\nDeveloping global model of juvenile observation probability...\n"); 
        juv.obs.global.model <- try(glm(formula = juv.obs.model, data = juvobs.data, family = "binomial"),
                                    silent = TRUE)
        
        if (class(juv.obs.global.model) == "try-error") {
          nox.juv.obs.model <- juv.obs.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.obs.model <- gsub(correction.sz1sz2, "", nox.juv.obs.model, fixed = TRUE)
          }
          
          if (nox.juv.obs.model != juv.obs.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            juv.obs.global.model <- try(glm(formula = nox.juv.obs.model, data = juvobs.data, family = "binomial"),
                                        silent = TRUE)
          }
          
          nopat.juv.obs.model <- nox.juv.obs.model
          if (!is.na(correction.patch)) {
            nopat.juv.obs.model <- gsub(correction.patch, "", nopat.juv.obs.model, fixed = TRUE)
          }
          
          if (class(juv.obs.global.model) == "try-error") {
            
            if (nox.juv.obs.model != nopat.juv.obs.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              juv.obs.global.model <- try(glm(formula = nopat.juv.obs.model, data = juvobs.data, family = "binomial"),
                                          silent = TRUE)
            }
          }
          
          noyr.juv.obs.model <- nopat.juv.obs.model
          if (!is.na(correction.year)) {
            noyr.juv.obs.model <- gsub(correction.year, "", noyr.juv.obs.model, fixed = TRUE)
          }
          
          if (class(juv.obs.global.model) == "try-error") {
            
            if (noyr.juv.obs.model != nopat.juv.obs.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              juv.obs.global.model <- try(glm(formula = noyr.juv.obs.model, data = juvobs.data, family = "binomial"),
                                          silent = TRUE)
            }
          }
          
          noind.juv.obs.model <- noyr.juv.obs.model
          if (!is.na(correction.indiv)) {
            noind.juv.obs.model <- gsub(correction.indiv, "", noind.juv.obs.model, fixed = TRUE)
          }
          
          if (class(juv.obs.global.model) == "try-error") {
            
            if (noind.juv.obs.model != noyr.juv.obs.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              juv.obs.global.model <- try(glm(formula = noind.juv.obs.model, data = juvobs.data, family = "binomial"),
                                          silent = TRUE)
            }
          }
          
          if (class(juv.obs.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for juvenile observation status.")
          }
        }
      } else if (!is.element(0, juvobs.data$obsstatus3)) {
        writeLines("\nJuvenile observation response is constant so will not model it.")
        juv.obs.model <- 1
        juv.obs.global.model <- 1
      } else {
        writeLines("\nJuvenile observation response is constant so will not model it.")
        juv.obs.model <- 1
        juv.obs.global.model <- 0}
    }
    
    if (juv.size.model != 1) {
      if (sizedist == "gaussian") {
        writeLines("\nDeveloping global model of size (Gaussian)...\n");
        juv.size.global.model <- try(lm(formula = juv.size.model, data = juvsize.data), silent = TRUE)
        
        if (class(juv.size.global.model) == "try-error") {
          nox.juv.size.model <- juv.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.size.model <- gsub(correction.sz1sz2, "", nox.juv.size.model, fixed = TRUE)
          }
          
          if (nox.juv.size.model != juv.size.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            juv.size.global.model <- try(lm(formula = nox.juv.size.model, data = juvsize.data),
                                         silent = TRUE)
          }
          
          nopat.juv.size.model <- nox.juv.size.model
          if (!is.na(correction.patch)) {
            nopat.juv.size.model <- gsub(correction.patch, "", nopat.juv.size.model, fixed = TRUE)
          }
          
          if (class(juv.size.global.model) == "try-error") {
            
            if (nox.juv.size.model != nopat.juv.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              juv.size.global.model <- try(lm(formula = nopat.juv.size.model, data = juvsize.data),
                                           silent = TRUE)
            }
          }
          
          noyr.juv.size.model <- nopat.juv.size.model
          if (!is.na(correction.year)) {
            noyr.juv.size.model <- gsub(correction.year, "", noyr.juv.size.model, fixed = TRUE)
          }
          
          if (class(juv.size.global.model) == "try-error") {
            
            if (noyr.juv.size.model != nopat.juv.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              juv.size.global.model <- try(lm(formula = noyr.juv.size.model, data = juvsize.data),
                                           silent = TRUE)
            }
          }
          
          noind.juv.size.model <- noyr.juv.size.model
          if (!is.na(correction.indiv)) {
            noind.juv.size.model <- gsub(correction.indiv, "", noind.juv.size.model, fixed = TRUE)
          }
          
          if (class(juv.size.global.model) == "try-error") {
            
            if (noind.juv.size.model != noyr.juv.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              juv.size.global.model <- try(lm(formula = noind.juv.size.model, data = juvsize.data),
                                           silent = TRUE)
            }
          }
          
          if (class(juv.size.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for juvenile size.")
          }
        }
      } else if (sizedist == "poisson") {
        writeLines("\nDeveloping global model of size (Poisson)...\n");
        juv.size.global.model <- try(glm(formula = juv.size.model, data = juvsize.data, family = "poisson"),
                                     silent = TRUE)
        
        if (class(juv.size.global.model) == "try-error") {
          nox.juv.size.model <- juv.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.size.model <- gsub(correction.sz1sz2, "", nox.juv.size.model, fixed = TRUE)
          }
          
          if (nox.juv.size.model != juv.size.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            juv.size.global.model <- try(glm(formula = nox.juv.size.model, data = juvsize.data, family = "poisson"),
                                         silent = TRUE)
          }
          
          nopat.juv.size.model <- nox.juv.size.model
          if (!is.na(correction.patch)) {
            nopat.juv.size.model <- gsub(correction.patch, "", nopat.juv.size.model, fixed = TRUE)
          }
          
          if (class(juv.size.global.model) == "try-error") {
            
            if (nox.juv.size.model != nopat.juv.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              juv.size.global.model <- try(glm(formula = nopat.juv.size.model, data = juvsize.data, family = "poisson"),
                                           silent = TRUE)
            }
          }
          
          noyr.juv.size.model <- nopat.juv.size.model
          if (!is.na(correction.year)) {
            noyr.juv.size.model <- gsub(correction.year, "", noyr.juv.size.model, fixed = TRUE)
          }
          
          if (class(juv.size.global.model) == "try-error") {
            
            if (noyr.juv.size.model != nopat.juv.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              juv.size.global.model <- try(glm(formula = noyr.juv.size.model, data = juvsize.data, family = "poisson"),
                                           silent = TRUE)
            }
          }
          
          noind.juv.size.model <- noyr.juv.size.model
          if (!is.na(correction.indiv)) {
            noind.juv.size.model <- gsub(correction.indiv, "", noind.juv.size.model, fixed = TRUE)
          }
          
          if (class(juv.size.global.model) == "try-error") {
            
            if (noind.juv.size.model != noyr.juv.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              juv.size.global.model <- try(glm(formula = noind.juv.size.model, data = juvsize.data, family = "poisson"),
                                           silent = TRUE)
            }
          }
          
          if (class(juv.size.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for juvenile size.")
          }
        }
      } else if (sizedist == "negbin") {
        writeLines("\nDeveloping global model of size (negative binomial)...\n");
        juv.size.global.model <- try(glm.nb(formula = juv.size.model, data = juvsize.data),
                                     silent = TRUE)
        
        if (class(juv.size.global.model) == "try-error") {
          nox.juv.size.model <- juv.size.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.size.model <- gsub(correction.sz1sz2, "", nox.juv.size.model, fixed = TRUE)
          }
          
          if (nox.juv.size.model != juv.size.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            juv.size.global.model <- try(glm.nb(formula = nox.juv.size.model, data = juvsize.data),
                                         silent = TRUE)
          }
          
          nopat.juv.size.model <- nox.juv.size.model
          if (!is.na(correction.patch)) {
            nopat.juv.size.model <- gsub(correction.patch, "", nopat.juv.size.model, fixed = TRUE)
          }
          
          if (class(juv.size.global.model) == "try-error") {
            
            if (nox.juv.size.model != nopat.juv.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              juv.size.global.model <- try(glm.nb(formula = nopat.juv.size.model, data = juvsize.data),
                                           silent = TRUE)
            }
          }
          
          noyr.juv.size.model <- nopat.juv.size.model
          if (!is.na(correction.year)) {
            noyr.juv.size.model <- gsub(correction.year, "", noyr.juv.size.model, fixed = TRUE)
          }
          
          if (class(juv.size.global.model) == "try-error") {
            
            if (noyr.juv.size.model != nopat.juv.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              juv.size.global.model <- try(glm.nb(formula = noyr.juv.size.model, data = juvsize.data),
                                           silent = TRUE)
            }
          }
          
          noind.juv.size.model <- noyr.juv.size.model
          if (!is.na(correction.indiv)) {
            noind.juv.size.model <- gsub(correction.indiv, "", noind.juv.size.model, fixed = TRUE)
          }
          
          if (class(juv.size.global.model) == "try-error") {
            
            if (noind.juv.size.model != noyr.juv.size.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              juv.size.global.model <- try(glm.nb(formula = noind.juv.size.model, data = juvsize.data),
                                           silent = TRUE)
            }
          }
          
          if (class(juv.size.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for juvenile size.")
          }
        }
      }
    }
    
    if (juv.repst.model != 1) {
      if (is.element(0, juvrepst.data$repstatus3) & is.element(1, juvrepst.data$repstatus3)) {
        writeLines("\nDeveloping global model of juvenile reproduction probability...\n"); 
        juv.repst.global.model <- try(glm(formula = juv.repst.model, data = juvrepst.data, family = "binomial"),
                                      silent = TRUE)
        
        if (class(juv.repst.global.model) == "try-error") {
          nox.juv.repst.model <- juv.repst.model
          
          if (!is.na(correction.sz1sz2)) {
            nox.juv.repst.model <- gsub(correction.sz1sz2, "", nox.juv.repst.model, fixed = TRUE)
          }
          
          if (nox.juv.repst.model != juv.repst.model) {
            writeLines("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
            juv.repst.global.model <- try(glm(formula = nox.juv.repst.model, data = juvrepst.data, family = "binomial"),
                                          silent = TRUE)
          }
          
          nopat.juv.repst.model <- nox.juv.repst.model
          if (!is.na(correction.patch)) {
            nopat.juv.repst.model <- gsub(correction.patch, "", nopat.juv.repst.model, fixed = TRUE)
          }
          
          if (class(juv.repst.global.model) == "try-error") {
            
            if (nox.juv.repst.model != nopat.juv.repst.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
              juv.repst.global.model <- try(glm(formula = nopat.juv.repst.model, data = juvrepst.data, family = "binomial"),
                                            silent = TRUE)
            }
          }
          
          noyr.juv.repst.model <- nopat.juv.repst.model
          if (!is.na(correction.year)) {
            noyr.juv.repst.model <- gsub(correction.year, "", noyr.juv.repst.model, fixed = TRUE)
          }
          
          if (class(juv.repst.global.model) == "try-error") {
            
            if (noyr.juv.repst.model != nopat.juv.repst.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
              juv.repst.global.model <- try(glm(formula = noyr.juv.repst.model, data = juvrepst.data, family = "binomial"),
                                            silent = TRUE)
            }
          }
          
          noind.juv.repst.model <- noyr.juv.repst.model
          if (!is.na(correction.indiv)) {
            noind.juv.repst.model <- gsub(correction.indiv, "", noind.juv.repst.model, fixed = TRUE)
          }
          
          if (class(juv.repst.global.model) == "try-error") {
            
            if (noind.juv.repst.model != noyr.juv.repst.model) {
              writeLines("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
              juv.repst.global.model <- try(glm(formula = noind.juv.repst.model, data = juvrepst.data, family = "binomial"),
                                            silent = TRUE)
            }
          }
          
          if (class(juv.repst.global.model) == "try-error") {
            writeLines("\nCould not properly estimate a global model for juvenile observation status.")
          }
        }
      } else if (!is.element(0, juvrepst.data$repstatus3)) {
        writeLines("\nJuvenile observation response is constant so will not model it.")
        juv.repst.model <- 1
        juv.repst.global.model <- 1
      } else {
        writeLines("\nJuvenile observation response is constant so will not model it.")
        juv.repst.model <- 1
        juv.repst.global.model <- 0}
    }
    
    writeLines("All global models developed.\n")
  }
  
  #This is the section where we dredge the models
  used.criterion <- gsub("&k", "", bestfit)
  
  if (suite != "cons") {options(na.action = "na.fail");
    if (full.surv.model != 1) {
      options(na.action = "na.fail")
      writeLines("\nCommencing dredge of survival probability...\n")
      if (!quiet) {
        surv.table <- try(dredge(surv.global.model, rank = used.criterion), silent = TRUE)
      } else {
        surv.table <- suppressWarnings(suppressMessages(try(dredge(surv.global.model, rank = used.criterion), silent = TRUE)))
      }
      if (any(class(surv.table) == "try-error")) {warning("Dredge of survival probability failed.")}
    }
    
    if (full.obs.model != 1) {
      options(na.action = "na.fail")
      writeLines("\nCommencing dredge of observation probability...\n")
      if (!quiet) {
        obs.table <- try(dredge(obs.global.model, rank = used.criterion), silent = TRUE)
      } else {
        obs.table <- suppressWarnings(suppressMessages(try(dredge(obs.global.model, rank = used.criterion), silent = TRUE)))
      }
      if (any(class(obs.table) == "try-error")) {warning("Dredge of observation probability failed.")}
    }
    
    if (full.size.model != 1) {
      options(na.action = "na.fail")
      writeLines("\nCommencing dredge of size...\n")
      if (!quiet) {
        size.table <- try(dredge(size.global.model, rank = used.criterion), silent = TRUE)
      } else {
        size.table <- suppressWarnings(suppressMessages(try(dredge(size.global.model, rank = used.criterion), silent = TRUE)))
      }
      if (any(class(size.table) == "try-error")) {warning("Dredge of size response failed.")}
    }
    
    if (full.repst.model != 1) {
      options(na.action = "na.fail")
      writeLines("\nCommencing dredge of reproduction probability...\n")
      if (!quiet) {
        repst.table <- try(dredge(repst.global.model, rank = used.criterion), silent = TRUE)
      } else {
        repst.table <- suppressWarnings(suppressMessages(try(dredge(repst.global.model, rank = used.criterion), silent = TRUE)))
      }
      if (any(class(repst.table) == "try-error")) {warning("Dredge of reproductive status probability failed.")}
    }
    
    if (full.fec.model != 1) {
      options(na.action = "na.fail")
      writeLines("\nCommencing dredge of fecundity...\n")
      if (!quiet) {
        fec.table <- try(dredge(fec.global.model, rank = used.criterion), silent = TRUE)
      } else {
        fec.table <- suppressWarnings(suppressMessages(try(dredge(fec.global.model, rank = used.criterion), silent = TRUE)))
      }
      if (any(class(fec.table) == "try-error")) {warning("Dredge of fecundity response failed.")}
    }
    
    if (juv.surv.model != 1) {
      options(na.action = "na.fail")
      writeLines("\nCommencing dredge of juvenile survival probability...\n")
      if (!quiet) {
        juvsurv.table <- try(dredge(juv.surv.global.model, rank = used.criterion), silent = TRUE)
      } else {
        juvsurv.table <- suppressWarnings(suppressMessages(try(dredge(juv.surv.global.model, rank = used.criterion), silent = TRUE)))
      }
      if (any(class(juvsurv.table) == "try-error")) {warning("Dredge of juvenilesurvival probability failed.")}
    }
    
    if (juv.obs.model != 1) {
      options(na.action = "na.fail")
      writeLines("\nCommencing dredge of juvenile observation probability...\n")
      if (!quiet) {
        juvobs.table <- try(dredge(juv.obs.global.model, rank = used.criterion), silent = TRUE)
      } else {
        juvobs.table <- suppressWarnings(suppressMessages(try(dredge(juv.obs.global.model, rank = used.criterion), silent = TRUE)))
      }
      if (any(class(juvobs.table) == "try-error")) {warning("Dredge of juvenile observation probability failed.")}
    }
    
    if (juv.size.model != 1) {
      options(na.action = "na.fail")
      writeLines("\nCommencing dredge of juvenile size...\n")
      if (!quiet) {
        juvsize.table <- try(dredge(juv.size.global.model, rank = used.criterion), silent = TRUE)
      } else {
        juvsize.table <- suppressWarnings(suppressMessages(try(dredge(juv.size.global.model, rank = used.criterion), silent = TRUE)))
      }
      if (any(class(juvsize.table) == "try-error")) {warning("Dredge of juvenile size response failed.")}
    }
    
    if (juv.repst.model != 1) {
      options(na.action = "na.fail")
      writeLines("\nCommencing dredge of juvenile reproduction status...\n")
      if (!quiet) {
        juvrepst.table <- try(dredge(juv.repst.global.model, rank = used.criterion), silent = TRUE)
      } else {
        juvrepst.table <- suppressWarnings(suppressMessages(try(dredge(juv.repst.global.model, rank = used.criterion), silent = TRUE)))
      }
      if (any(class(juvrepst.table) == "try-error")) {warning("Dredge of juvenile reproduction status failed.")}
    }
    
    writeLines("\nFinished dredging all vital rates.\n")
  }

  if (stringr::str_detect(bestfit, "&k")) {
    if (any(class(surv.table) == "model.selection")) {
      relevant.surv.models <- which(surv.table$delta <= 2)
      df.surv.models <- which(surv.table$df[relevant.surv.models] == min(surv.table$df[relevant.surv.models]))
      if (length(df.surv.models) > 1) {
        df.surv.models <- min(df.surv.models)
      }
      if (quiet) {
        surv.bf <- suppressWarnings(suppressMessages(eval(getCall(surv.table, min(df.surv.models)))))
      } else {
        surv.bf <- eval(getCall(surv.table, min(df.surv.models)))
      }
      surv.ind <- length(unique(surv.data[, which(names(surv.data) == indiv)]))
      surv.trans <- dim(surv.data)[1]
    } else if (any(class(surv.table) == "try-error")) {
      surv.bf <- surv.global.model
      surv.ind <- length(unique(surv.data[, which(names(surv.data) == indiv)]))
      surv.trans <- dim(surv.data)[1]
    } else {
      surv.bf <- surv.global.model
      surv.ind <- 0
      surv.trans <- 0
    }
    
    if (any(class(obs.table) == "model.selection")) {
      relevant.obs.models <- which(obs.table$delta <= 2)
      df.obs.models <- which(obs.table$df[relevant.obs.models] == min(obs.table$df[relevant.obs.models]))
      if (length(df.obs.models) > 1) {
        df.obs.models <- min(df.obs.models)
      }
      if (quiet) {
        obs.bf <- suppressWarnings(suppressMessages(eval(getCall(obs.table, min(df.obs.models)))))
      } else {
        obs.bf <- eval(getCall(obs.table, min(df.obs.models)))
      }
      obs.ind <- length(unique(obs.data[, which(names(obs.data) == indiv)]))
      obs.trans <- dim(obs.data)[1]
    } else if (any(class(obs.table) == "try-error")) {
      obs.bf <- obs.global.model
      obs.ind <- length(unique(obs.data[, which(names(obs.data) == indiv)]))
      obs.trans <- dim(obs.data)[1]
    } else {
      obs.bf <- obs.global.model
      obs.ind <- 0
      obs.trans <- 0
    }
    
    if (any(class(size.table) == "model.selection")) {
      relevant.size.models <- which(size.table$delta <= 2)
      df.size.models <- which(size.table$df[relevant.size.models] == min(size.table$df[relevant.size.models]))
      if (length(df.size.models) > 1) {
        df.size.models <- min(df.size.models)
      }
      if (quiet) {
        size.bf <- suppressWarnings(suppressMessages(eval(getCall(size.table, min(df.size.models)))))
      } else {
        size.bf <- eval(getCall(size.table, min(df.size.models)))
      }
      size.ind <- length(unique(size.data[, which(names(size.data) == indiv)]))
      size.trans <- dim(size.data)[1]
    } else if (any(class(size.table) == "try-error")) {
      size.bf <- size.global.model
      size.ind <- length(unique(size.data[, which(names(size.data) == indiv)]))
      size.trans <- dim(size.data)[1]
    } else {
      size.bf <- size.global.model
      size.ind <- 0
      size.trans <- 0
    }
    
    if (any(class(repst.table) == "model.selection")) {
      relevant.repst.models <- which(repst.table$delta <= 2)
      df.repst.models <- which(repst.table$df[relevant.repst.models] == min(repst.table$df[relevant.repst.models]))
      if (length(df.repst.models) > 1) {
        df.repst.models <- min(df.repst.models)
      }
      if (quiet) {
        repst.bf <- suppressWarnings(suppressMessages(eval(getCall(repst.table, min(df.repst.models)))))
      } else {
        repst.bf <- eval(getCall(repst.table, min(df.repst.models)))
      }
      repst.ind <- length(unique(repst.data[, which(names(repst.data) == indiv)]))
      repst.trans <- dim(repst.data)[1]
    } else if (any(class(repst.table) == "try-error")) {
      repst.bf <- repst.global.model
      repst.ind <- length(unique(repst.data[, which(names(repst.data) == indiv)]))
      repst.trans <- dim(repst.data)[1]
    } else {
      repst.bf <- repst.global.model
      repst.ind <- 0
      repst.trans <- 0
    }
    
    if (any(class(fec.table) == "model.selection")) {
      relevant.fec.models <- which(fec.table$delta <= 2)
      df.fec.models <- which(fec.table$df[relevant.fec.models] == min(fec.table$df[relevant.fec.models]))
      if (length(df.fec.models) > 1) {
        df.fec.models <- min(df.fec.models)
      }
      if (quiet) {
        fec.bf <- suppressWarnings(suppressMessages(eval(getCall(fec.table, min(df.fec.models)))))
      } else {
        fec.bf <- eval(getCall(fec.table, min(df.fec.models)))
      }
      fec.ind <- length(unique(fec.data[, which(names(fec.data) == indiv)]))
      fec.trans <- dim(fec.data)[1]
    } else if (any(class(fec.table) == "try-error")) {
      fec.bf <- fec.global.model
      fec.ind <- length(unique(fec.data[, which(names(fec.data) == indiv)]))
      fec.trans <- dim(fec.data)[1]
    } else {
      fec.bf <- fec.global.model
      fec.ind <- 0
      fec.trans <- 0
    }
    
    if (any(class(juvsurv.table) == "model.selection")) {
      relevant.juv.surv.models <- which(juvsurv.table$delta <= 2)
      df.juv.surv.models <- which(juvsurv.table$df[relevant.juv.surv.models] == min(juvsurv.table$df[relevant.juv.surv.models]))
      if (length(df.juv.surv.models) > 1) {
        df.juv.surv.models <- min(df.juv.surv.models)
      }
      if (quiet) {
        juvsurv.bf <- suppressWarnings(suppressMessages(eval(getCall(juvsurv.table, min(df.juv.surv.models)))))
      } else {
        juvsurv.bf <- eval(getCall(juvsurv.table, min(df.juv.surv.models)))
      }
      juvsurv.ind <- length(unique(juvsurv.data[, which(names(juvsurv.data) == indiv)]))
      juvsurv.trans <- dim(juvsurv.data)[1]
    } else if (any(class(juvsurv.table) == "try-error")) {
      juvsurv.bf <- juv.surv.global.model
      juvsurv.ind <- length(unique(juvsurv.data[, which(names(juvsurv.data) == indiv)]))
      juvsurv.trans <- dim(juvsurv.data)[1]
    } else {
      juvsurv.bf <- juv.surv.global.model
      juvsurv.ind <- 0
      juvsurv.trans <- 0
    }
    
    if (any(class(juvobs.table) == "model.selection")) {
      relevant.juv.obs.models <- which(juvobs.table$delta <= 2)
      df.juv.obs.models <- which(juvobs.table$df[relevant.juv.obs.models] == min(juvobs.table$df[relevant.juv.obs.models]))
      if (length(df.juv.obs.models) > 1) {
        df.juv.obs.models <- min(df.juv.obs.models)
      }
      if (quiet) {
        juvobs.bf <- suppressWarnings(suppressMessages(eval(getCall(juvobs.table, min(df.juv.obs.models)))))
      } else {
        juvobs.bf <- eval(getCall(juvobs.table, min(df.juv.obs.models)))
      }
      juvobs.ind <- length(unique(juvobs.data[, which(names(juvobs.data) == indiv)]))
      juvobs.trans <- dim(juvobs.data)[1]
    } else if (any(class(juvobs.table) == "try-error")) {
      juvobs.bf <- juv.obs.global.model
      juvobs.ind <- length(unique(juvobs.data[, which(names(juvobs.data) == indiv)]))
      juvobs.trans <- dim(juvobs.data)[1]
    } else {
      juvobs.bf <- juv.obs.global.model
      juvobs.ind <- 0
      juvobs.trans <- 0
    }
    
    if (any(class(juvsize.table) == "model.selection")) {
      relevant.juv.size.models <- which(juvsize.table$delta <= 2)
      df.juv.size.models <- which(juvsize.table$df[relevant.juv.size.models] == min(juvsize.table$df[relevant.juv.size.models]))
      if (length(df.juv.size.models) > 1) {
        df.juv.size.models <- min(df.juv.size.models)
      }
      if (quiet) {
        juvsize.bf <- suppressWarnings(suppressMessages(eval(getCall(juvsize.table, min(df.juv.size.models)))))
      } else {
        juvsize.bf <- eval(getCall(juvsize.table, min(df.juv.size.models)))
      }
      juvsize.ind <- length(unique(juvsize.data[, which(names(juvsize.data) == indiv)]))
      juvsize.trans <- dim(juvsize.data)[1]
    } else if (any(class(juvsize.table) == "try-error")) {
      juvsize.bf <- juv.size.global.model
      juvsize.ind <- length(unique(juvsize.data[, which(names(juvsize.data) == indiv)]))
      juvsize.trans <- dim(juvsize.data)[1]
    } else {
      juvsize.bf <- juv.size.global.model
      juvsize.ind <- 0
      juvsize.trans <- 0
    }
    
    if (any(class(juvrepst.table) == "model.selection")) {
      relevant.juv.repst.models <- which(juvrepst.table$delta <= 2)
      df.juv.repst.models <- which(juvrepst.table$df[relevant.juv.repst.models] == min(juvrepst.table$df[relevant.juv.repst.models]))
      if (length(df.juv.repst.models) > 1) {
        df.juv.repst.models <- min(df.juv.repst.models)
      }
      if (quiet) {
        juvrepst.bf <- suppressWarnings(suppressMessages(eval(getCall(juvrepst.table, min(df.juv.repst.models)))))
      } else {
        juvrepst.bf <- eval(getCall(juvrepst.table, min(df.juv.repst.models)))
      }
      juvrepst.ind <- length(unique(juvrepst.data[, which(names(juvrepst.data) == indiv)]))
      juvrepst.trans <- dim(juvrepst.data)[1]
    } else if (any(class(juvrepst.table) == "try-error")) {
      juvrepst.bf <- juv.repst.global.model
      juvrepst.ind <- length(unique(juvrepst.data[, which(names(juvrepst.data) == indiv)]))
      juvrepst.trans <- dim(juvrepst.data)[1]
    } else {
      juvrepst.bf <- juv.repst.global.model
      juvrepst.ind <- 0
      juvrepst.trans <- 0
    }
  }
  
  if (!stringr::str_detect(bestfit, "&k")) {
    if (any(class(surv.table) == "model.selection")) {
      if (quiet) {
        surv.bf <- suppressWarnings(suppressMessages(eval(getCall(surv.table, 1))))
      } else {
        surv.bf <- eval(getCall(surv.table, 1))
      }
      surv.ind <- length(unique(surv.data[, which(names(surv.data) == indiv)]))
      surv.trans <- dim(surv.data)[1]
    } else if (any(class(surv.table) == "try-error")) {
      surv.bf <- surv.global.model
      surv.ind <- length(unique(surv.data[, which(names(surv.data) == indiv)]))
      surv.trans <- dim(surv.data)[1]
    } else {
      surv.bf <- surv.global.model
      surv.ind <- 0
      surv.trans <- 0
    }
    
    if (any(class(obs.table) == "model.selection")) {
      if (quiet) {
        obs.bf <- suppressWarnings(suppressMessages(eval(getCall(obs.table, 1))))
      } else {
        obs.bf <- eval(getCall(obs.table, 1))
      }
      obs.ind <- length(unique(obs.data[, which(names(obs.data) == indiv)]))
      obs.trans <- dim(obs.data)[1]
    } else if (any(class(obs.table) == "try-error")) {
      obs.bf <- obs.global.model
      obs.ind <- length(unique(obs.data[, which(names(obs.data) == indiv)]))
      obs.trans <- dim(obs.data)[1]
    } else {
      obs.bf <- obs.global.model
      obs.ind <- 0
      obs.trans <- 0
    }
    
    if (any(class(size.table) == "model.selection")) {
      if (quiet) {
        size.bf <- suppressWarnings(suppressMessages(eval(getCall(size.table, 1))))
      } else {
        size.bf <- eval(getCall(size.table, 1))
      }
      size.ind <- length(unique(size.data[, which(names(size.data) == indiv)]))
      size.trans <- dim(size.data)[1]
    } else if (any(class(size.table) == "try-error")) {
      size.bf <- size.global.model
      size.ind <- length(unique(size.data[, which(names(size.data) == indiv)]))
      size.trans <- dim(size.data)[1]
    } else {
      size.bf <- size.global.model
      size.ind <- 0
      size.trans <- 0
    }
    
    if (any(class(repst.table) == "model.selection")) {
      if (quiet) {
        repst.bf <- suppressWarnings(suppressMessages(eval(getCall(repst.table, 1))))
      } else {
        repst.bf <- eval(getCall(repst.table, 1))
      }
      repst.ind <- length(unique(repst.data[, which(names(repst.data) == indiv)]))
      repst.trans <- dim(repst.data)[1]
    } else if (any(class(repst.table) == "try-error")) {
      repst.bf <- repst.global.model
      repst.ind <- length(unique(repst.data[, which(names(repst.data) == indiv)]))
      repst.trans <- dim(repst.data)[1]
    } else {
      repst.bf <- repst.global.model
      repst.ind <- 0
      repst.trans <- 0
    }
    
    if (any(class(fec.table) == "model.selection")) {
      if (quiet) {
        fec.bf <- suppressWarnings(suppressMessages(eval(getCall(fec.table, 1))))
      } else {
        fec.bf <- eval(getCall(fec.table, 1))
      }
      fec.ind <- length(unique(fec.data[, which(names(fec.data) == indiv)]))
      fec.trans <- dim(fec.data)[1]
    } else if (any(class(fec.table) == "try-error")) {
      fec.bf <- fec.global.model
      fec.ind <- length(unique(fec.data[, which(names(fec.data) == indiv)]))
      fec.trans <- dim(fec.data)[1]
    } else {
      fec.bf <- fec.global.model
      fec.ind <- 0
      fec.trans <- 0
    }
    
    if (any(class(juvsurv.table) == "model.selection")) {
      if (quiet) {
        juvsurv.bf <- suppressWarnings(suppressMessages(eval(getCall(juvsurv.table, 1))))
      } else {
        juvsurv.bf <- eval(getCall(juvsurv.table, 1))
      }
      juvsurv.ind <- length(unique(juvsurv.data[, which(names(juvsurv.data) == indiv)]))
      juvsurv.trans <- dim(juvsurv.data)[1]
    } else if (any(class(juvsurv.table) == "try-error")) {
      juvsurv.bf <- juv.surv.global.model
      juvsurv.ind <- length(unique(juvsurv.data[, which(names(juvsurv.data) == indiv)]))
      juvsurv.trans <- dim(juvsurv.data)[1]
    } else {
      juvsurv.bf <- juv.surv.global.model
      juvsurv.ind <- 0
      juvsurv.trans <- 0
    }
    
    if (any(class(juvobs.table) == "model.selection")) {
      if (quiet) {
        juvobs.bf <- suppressWarnings(suppressMessages(eval(getCall(juvobs.table, 1))))
      } else {
        juvobs.bf <- eval(getCall(juvobs.table, 1))
      }
      juvobs.ind <- length(unique(juvobs.data[, which(names(juvobs.data) == indiv)]))
      juvobs.trans <- dim(juvobs.data)[1]
    } else if (any(class(juvobs.table) == "try-error")) {
      juvobs.bf <- juv.obs.global.model
      juvobs.ind <- length(unique(juvobs.data[, which(names(juvobs.data) == indiv)]))
      juvobs.trans <- dim(juvobs.data)[1]
    } else {
      juvobs.bf <- juv.obs.global.model
      juvobs.ind <- 0
      juvobs.trans <- 0
    }
    
    if (any(class(juvsize.table) == "model.selection")) {
      if (quiet) {
        juvsize.bf <- suppressWarnings(suppressMessages(eval(getCall(juvsize.table, 1))))
      } else {
        juvsize.bf <- eval(getCall(juvsize.table, 1))
      }
      juvsize.ind <- length(unique(juvsize.data[, which(names(juvsize.data) == indiv)]))
      juvsize.trans <- dim(juvsize.data)[1]
    } else if (any(class(juvsize.table) == "try-error")) {
      juvsize.bf <- juv.size.global.model
      juvsize.ind <- length(unique(juvsize.data[, which(names(juvsize.data) == indiv)]))
      juvsize.trans <- dim(juvsize.data)[1]
    } else {
      juvsize.bf <- juv.size.global.model
      juvsize.ind <- 0
      juvsize.trans <- 0
    }
    
    if (any(class(juvrepst.table) == "model.selection")) {
      if (quiet) {
        juvrepst.bf <- suppressWarnings(suppressMessages(eval(getCall(juvrepst.table, 1))))
      } else {
        juvrepst.bf <- eval(getCall(juvrepst.table, 1))
      }
      juvrepst.ind <- length(unique(juvrepst.data[, which(names(juvrepst.data) == indiv)]))
      juvrepst.trans <- dim(juvrepst.data)[1]
    } else if (any(class(juvrepst.table) == "try-error")) {
      juvrepst.bf <- juv.repst.global.model
      juvrepst.ind <- length(unique(juvrepst.data[, which(names(juvrepst.data) == indiv)]))
      juvrepst.trans <- dim(juvrepst.data)[1]
    } else {
      juvrepst.bf <- juv.repst.global.model
      juvrepst.ind <- 0
      juvrepst.trans <- 0
    }
  }
  
  qcoutput <- cbind.data.frame(c("survival", "observation", "size", "reproduction", "fecundity", 
                                 "juvenile_survival", "juvnile_observation", "juvenile_size", 
                                 "juvenile_reproduction"), c(surv.ind, obs.ind, size.ind, repst.ind, 
                                 fec.ind, juvsurv.ind, juvobs.ind, juvsize.ind, juvrepst.ind), 
                                 c(surv.trans, obs.trans, size.trans, repst.trans, fec.trans, 
                                 juvsurv.trans, juvobs.trans, juvsize.trans, juvrepst.trans))
  names(qcoutput) <- c("vital_rate", "individuals", "transitions")
  
  if (show.model.tables == FALSE) {
    surv.table <- NA
    obs.table <- NA
    size.table <- NA
    repst.table <- NA
    fec.table <- NA
    
    juvsurv.table <- NA
    juvobs.table <- NA
    juvsize.table <- NA
    juvrepst.table <- NA
  }
  
  #Now we develop the final output, creating a new S3 class to do it
  full.output <- list(survival_model = surv.bf, observation_model = obs.bf, size_model = size.bf, 
                      repstatus_model = repst.bf, fecundity_model = fec.bf, 
                      juv_survival_model = juvsurv.bf, juv_observation_model = juvobs.bf, 
                      juv_size_model = juvsize.bf, juv_reproduction_model = juvrepst.bf,
                      survival_table = surv.table, observation_table = obs.table, 
                      size_table = size.table, repstatus_table = repst.table,
                      fecundity_table = fec.table, juv_survival_table = juvsurv.table, 
                      juv_observation_table = juvobs.table, juv_size_table = juvsize.table,
                      juv_reproduction_table = juvrepst.table, paramnames = paramnames, 
                      criterion = bestfit, qc = qcoutput)
  class(full.output) <- "lefkoMod"
  
  return(full.output)
}

#' Summary of Class "lefkoMod"
#' 
#' A function to simplify the immediate viewable output for an R object of class
#' \code{lefkoMod}. This function shows the best-fit models, summarizes the numbers
#' of models in the model tables, shows the criterion used to determine the
#' best-fit models, and provides some basic quality control information about
#' the models.
#' 
#' @param object An R object of class \code{lefkoMod} resulting from 
#' \code{\link{modelsearch}()}.
#' @param ... Other parameters.
#' 
#' @return A summary of the object, showing the best-fit models for all vital
#' rates, with constants of 0 or 1 used for unestimated models. Thie is
#' followed by a summary of the number of models tested per vital rate, and
#' a table showing the names of the parameters used to model vital rates and
#' represent tested factors. At the end is a section describing the number
#' of individuals and individual transitions used to estimate each vital
#' rate.
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
#' lathmodelsln2 <- modelsearch(lathvertln, historical = FALSE, approach = "lme4", suite = "main",
#'                              vitalrates = c("surv", "obs", "size", "repst", "fec"), 
#'                              juvestimate = "Sdl", bestfit = "AICc&k", sizedist = "gaussian", 
#'                              fecdist = "poisson", indiv = "individ", patch = "patchid", 
#'                              year = "year2", year.as.random = TRUE, patch.as.random = TRUE,
#'                              show.model.tables = TRUE)
#' 
#' summary(lathmodelsln2)
#' }
#' 
#' @export
summary.lefkoMod <- function(object, ...) {
  
  modelsuite <- object

  totalmodels <- length(which(c(is.numeric(modelsuite$survival_model), is.numeric(modelsuite$observation_model), 
                                is.numeric(modelsuite$size_model), is.numeric(modelsuite$repstatus_model), 
                                is.numeric(modelsuite$fecundity_model), is.numeric(modelsuite$juv_survival_model),
                                is.numeric(modelsuite$juv_observation_model), is.numeric(modelsuite$juv_size_model),
                                is.numeric(modelsuite$juv_reproduction_model)) == FALSE))
  
  writeLines(paste0("This LefkoMod object includes ", totalmodels, " linear models."))
  writeLines(paste0("Best-fit model criterion used: ", modelsuite$criterion))
  writeLines("\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  writeLines("Survival model:")
  print(modelsuite$survival_model)
  writeLines("\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  writeLines("\nObservation model:")
  print(modelsuite$observation_model)
  writeLines("\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  writeLines("\nSize model:")
  print(modelsuite$size_model)
  writeLines("\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  writeLines("\nReproductive status model:")
  print(modelsuite$repstatus_model)
  writeLines("\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  writeLines("\nFecundity model:")
  print(modelsuite$fecundity_model)
  writeLines("\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  writeLines("Juvenile survival model:")
  print(modelsuite$juv_survival_model)
  writeLines("\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  writeLines("\nJuvenile observation model:")
  print(modelsuite$juv_observation_model)
  writeLines("\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  writeLines("\nJuvenile size model:")
  print(modelsuite$juv_size_model)
  writeLines("\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  writeLines("\nJuvenile reproduction model:")
  print(modelsuite$juv_reproduction_model)
  
  writeLines("\n\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  if (!is.logical(modelsuite$survival_model) & is.element("model.selection", class(modelsuite$survival_table))) {
    writeLines(paste0("\nNumber of models in survival table:", dim(modelsuite$survival_table)[1]))
  } else if (!is.logical(modelsuite$survival_model)) {
    writeLines("\nNumber of models in survival table: 1")
  } else {
    writeLines("\nSurvival table not estimated")
  }
  
  if (!is.logical(modelsuite$observation_model) & is.element("model.selection", class(modelsuite$observation_table))) {
    writeLines(paste0("\nNumber of models in observation table:", dim(modelsuite$observation_table)[1]))
  } else if (!is.logical(modelsuite$observation_model)) {
    writeLines("\nNumber of models in observation table: 1")
  } else {
    writeLines("\nObservation table not estimated")
  }
  
  if (!is.logical(modelsuite$size_model) & is.element("model.selection", class(modelsuite$size_table))) {
    writeLines(paste0("\nNumber of models in size table:", dim(modelsuite$size_table)[1]))
  } else if (!is.logical(modelsuite$size_model)) {
    writeLines("\nNumber of models in size table: 1")
  } else {
    writeLines("\nSize table not estimated")
  }
  
  if (!is.logical(modelsuite$repstatus_model) & is.element("model.selection", class(modelsuite$repstatus_table))) {
    writeLines(paste0("\nNumber of models in reproduction status table:", dim(modelsuite$repstatus_table)[1]))
  } else if (!is.logical(modelsuite$repstatus_model)) {
    writeLines("\nNumber of models in reproduction status table: 1")
  } else {
    writeLines("\nReproduction status table not estimated")
  }
  
  if (!is.logical(modelsuite$fecundity_model) & is.element("model.selection", class(modelsuite$fecundity_table))) {
    writeLines(paste0("\nNumber of models in fecundity table:", dim(modelsuite$fecundity_table)[1]))
  } else if (!is.logical(modelsuite$fecundity_model)) {
    writeLines("\nNumber of models in fecundity table: 1")
  } else {
    writeLines("\nFecundity table not estimated")
  }
  
  if (!is.logical(modelsuite$juv_survival_model) & is.element("model.selection", class(modelsuite$juv_survival_table))) {
    writeLines(paste0("\nNumber of models in juvenile survival table:", dim(modelsuite$juv_survival_table)[1]))
  } else if (!is.logical(modelsuite$juv_survival_model)) {
    writeLines("\nNumber of models in juvenile survival table: 1")
  } else {
    writeLines("\nJuvenile survival table not estimated")
  }
  
  if (!is.logical(modelsuite$juv_observation_model) & is.element("model.selection", class(modelsuite$juv_observation_table))) {
    writeLines(paste0("\nNumber of models in juvenile observation table:", dim(modelsuite$juv_observation_table)[1]))
  } else if (!is.logical(modelsuite$juv_observation_model)) {
    writeLines("\nNumber of models in juvenile observation table: 1")
  } else {
    writeLines("\nJuvenile observation table not estimated")
  }
  
  if (!is.logical(modelsuite$juv_size_model) & is.element("model.selection", class(modelsuite$juv_size_table))) {
    writeLines(paste0("\nNumber of models in juvenile size table:", dim(modelsuite$juv_size_table)[1]))
  } else if (!is.logical(modelsuite$juv_size_model)) {
    writeLines("\nNumber of models in juvenile size table: 1")
  } else {
    writeLines("\nJuvenile size table not estimated")
  }
  
  if (!is.logical(modelsuite$juv_reproduction_model) & is.element("model.selection", class(modelsuite$juv_reproduction_table))) {
    writeLines(paste0("\nNumber of models in juvenile reproduction table:", dim(modelsuite$juv_reproduction_table)[1]))
  } else if (!is.logical(modelsuite$juv_reproduction_model)) {
    writeLines("\nNumber of models in juvenile reproduction table: 1")
  } else {
    writeLines("\nJuvenile reproduction table not estimated")
  }
  
  
  writeLines("\n\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  writeLines("\nGeneral model parameter names (column 1), and specific names used in these models (column 2):")
  print.data.frame(modelsuite$paramnames[c(1,2)])
  
  writeLines("\n\n\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500")
  writeLines("\nQuality control:\n")
  if (modelsuite$qc[1,2] > 0) {
    writeLines(paste0("Survival estimated with ", modelsuite$qc[1,2], " individuals and ", modelsuite$qc[1,3], " individual transitions."))
  } else {
    writeLines("Survival not estimated.")
  }
  
  if (modelsuite$qc[2,2] > 0) {
    writeLines(paste0("Observation estimated with ", modelsuite$qc[2,2], " individuals and ", modelsuite$qc[2,3], " individual transitions."))
  } else {
    writeLines("Observation probability not estimated.")
  }
  
  if (modelsuite$qc[3,2] > 0) {
    writeLines(paste0("Size estimated with ", modelsuite$qc[3,2], " individuals and ", modelsuite$qc[3,3], " individual transitions."))
  } else {
    writeLines("Size transition not estimated.")
  }
  
  if (modelsuite$qc[4,2] > 0) {
    writeLines(paste0("Reproductive status estimated with ", modelsuite$qc[4,2], " individuals and ", modelsuite$qc[4,3], " individual transitions."))
  } else {
    writeLines("Reproduction probability not estimated.")
  }
  
  if (modelsuite$qc[5,2] > 0) {
    writeLines(paste0("Fecundity estimated with ", modelsuite$qc[5,2], " individuals and ", modelsuite$qc[5,3], " individual transitions."))
  } else {
    writeLines("Fecundity not estimated.")
  }
  
  if (modelsuite$qc[6,2] > 0) {
    writeLines(paste0("Juvenile survival estimated with ", modelsuite$qc[6,2], " individuals and ", modelsuite$qc[6,3], " individual transitions."))
  } else {
    writeLines("Juvenile survival not estimated.")
  }
  
  if (modelsuite$qc[7,2] > 0) {
    writeLines(paste0("Juvenile observation estimated with ", modelsuite$qc[7,2], " individuals and ", modelsuite$qc[7,3], " individual transitions."))
  } else {
    writeLines("Juvenile observation probability not estimated.")
  }
  
  if (modelsuite$qc[8,2] > 0) {
    writeLines(paste0("Juvenile size estimated with ", modelsuite$qc[8,2], " individuals and ", modelsuite$qc[8,3], " individual transitions."))
  } else {
    writeLines("Juvenile size transition not estimated.")
  }
  
  if (modelsuite$qc[9,2] > 0) {
    writeLines(paste0("Juvenile reproduction estimated with ", modelsuite$qc[9,2], " individuals and ", modelsuite$qc[9,3], " individual transitions."))
  } else {
    writeLines("Juvenile reproduction probability not estimated.")
  }
  
  return()
}

#' Extract Required Coefficient Values From Vital Rate Models
#' 
#' \code{.modelextract()} is an S3 generic function designed to extract the
#' estimated coefficient values from linear models created to estimate vital rates
#' in \code{'lefko3'}. Used to supply coefficients to \code{\link{flefko3}()} and 
#' \code{\link{flefko2}()}.
#' 
#' @param model Model estimated through functions supported by \code{'lefko3'}.
#' @param paramnames Data frame giving the names of standard coefficients required
#' by matrix creation functions.
#' 
#' @return This function typically returns a list with elements corresponding to
#' different classes of coefficients.
#' 
#' @seealso \code{\link{.modelextract.glmmTMB}()}
#' @seealso \code{\link{.modelextract.merMod}()}
#' @seealso \code{\link{.modelextract.lmerMod}()}
#' @seealso \code{\link{.modelextract.glm}()}
#' @seealso \code{\link{.modelextract.lm}()}
#' @seealso \code{\link{.modelextract.numeric}()}
#' @seealso \code{\link{.modelextract.logical}()}
#' 
#' @keywords internal
#' @noRd
.modelextract <- function(model, ...) UseMethod(".modelextract")

#' Extract Required Coefficient Values From glmmTMB-estimated Vital Rate Models
#' 
#' \code{.modelextract.glmmTMB()} extracts coefficient values from linear models 
#' estimated through the \code{\link[glmmTMB]{glmmTMB}()} function, to estimate
#' vital rates in \code{'lefko3'}. Used to supply coefficients to \code{\link{flefko3}()} and 
#' \code{\link{flefko2}()}.
#' 
#' @param model Model estimated through \code{\link[glmmTMB]{glmmTMB}()}.
#' @param paramnames Data frame giving the names of standard coefficients required
#' by matrix creation functions.
#' 
#' @return This function returns a list with the following elements:
#' \item{coefficients}{Vector of fixed effect coefficients.}
#' \item{years}{Vector of time step coefficients, typically random effect.}
#' \item{patches}{Vector of patch coefficients, typically random effect.}
#' \item{variances}{Residual variance terms associated with model.}
#' \item{family}{Distribution of response term.}
#' \item{sigma}{Standard deviation of response.}
#' 
#' @seealso \code{\link{.modelextract}()}
#' @seealso \code{\link{.modelextract.merMod}()}
#' @seealso \code{\link{.modelextract.lmerMod}()}
#' @seealso \code{\link{.modelextract.glm}()}
#' @seealso \code{\link{.modelextract.lm}()}
#' @seealso \code{\link{.modelextract.numeric}()}
#' @seealso \code{\link{.modelextract.logical}()}
#' 
#' @keywords internal
#' @noRd
.modelextract.glmmTMB <- function(model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random) {
  
  size2.coef <- 0
  size1.coef <- 0
  flw2.coef <- 0
  flw1.coef <- 0
  
  size1.size2.coef <- 0
  flw1.flw2.coef <- 0
  size1.flw1.coef <- 0
  size2.flw2.coef <- 0
  size1.flw2.coef <- 0
  size2.flw1.coef <- 0
  
  age.coef <- 0
  age.size1.coef <- 0
  age.size2.coef <- 0
  age.flw1.coef <- 0
  age.flw2.coef <- 0
  
  yintercept <- glmmTMB::fixef(model)[["cond"]]["(Intercept)"]
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "repst2")), 2]])) {flw2.coef <- glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "repst2")), 2]]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "repst1")), 2]])) {flw1.coef <- glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "repst1")), 2]]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])])) {flw1.flw2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])])) {flw1.flw2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "size2")), 2]])) {size2.coef <- glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "size2")), 2]]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "size1")), 2]])) {size1.coef <- glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "size1")), 2]]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])])) {size1.size2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])])) {size1.size2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])])) {size1.flw1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])])) {size2.flw2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])])) {size1.flw1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])])) {size2.flw2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])])) {size2.flw1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])])) {size1.flw2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])])) {size2.flw1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])])) {size1.flw2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "age")), 2]])) {age.coef <- glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "age")), 2]]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])])) {age.size1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])])) {age.size1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])])) {age.size2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])])) {age.size2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])])) {age.flw1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])])) {age.flw1.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])]}
  
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])])) {age.flw2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])]}
  if (!is.na(glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])])) {age.flw2.coef <- glmmTMB::fixef(model)[["cond"]][paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])]}
  
  year.coefs <- glmmTMB::ranef(model)[["cond"]][[paramnames[(which(paramnames$mainparams == "year2")), 2]]]
  
  if (!all(is.na(glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "patch")), 2]]))) {
    patch.coefs <- glmmTMB::fixef(model)[["cond"]][paramnames[(which(paramnames$mainparams == "patch")), 2]]
  } else if (!all(is.na(glmmTMB::ranef(model)[["cond"]][[paramnames[(which(paramnames$mainparams == "patch")), 2]]]))) {
    patch.coefs <- glmmTMB::ranef(model)[["cond"]][[paramnames[(which(paramnames$mainparams == "patch")), 2]]]
  } else {
    patch.coefs <- 0
  }
  
  coef.vec <- c(yintercept, flw1.coef, flw2.coef, size1.coef, size2.coef, flw1.flw2.coef, size1.size2.coef, 
                size1.flw1.coef, size2.flw2.coef, size2.flw1.coef, size1.flw2.coef, age.coef, age.size1.coef, 
                age.size2.coef, age.flw1.coef, age.flw2.coef)
  
  rvars <- NA
  
  coef.list <- list(coefficients = coef.vec, years = year.coefs, patches = patch.coefs,
                    variances = rvars, family = stats::family(model)$family, sigma = NA)
  
  #Here we correct for missing years
  yeardiff <- setdiff(mainyears, rownames(coef.list$years))
  
  if (length(yeardiff) > 0) {
    if (year.as.random == TRUE) {
      newdevs <- rnorm(length(yeardiff), mean = 0, sd = sd(coef.list$years[,1]))
    } else {
      newdevs <- rep(0, (length(yeardiff)))
    }
    coef.list$years <- rbind.data.frame(setNames(as.data.frame(as.matrix(newdevs, ncol = 1)), names(coef.list$years)), coef.list$years)
    yearlabels <- c(yeardiff, rownames(coef.list$years)[(length(yeardiff) + 1):(length(coef.list$years[,1]))])
    coef.list$years <- as.data.frame(coef.list$years[order(yearlabels),])
    rownames(coef.list$years) <- sort(yearlabels)
    
  }
  
  #This next part addresses missing patches
  if (class(coef.list$patches) == "data.frame") {
    patchdiff <- setdiff(mainpatches, rownames(coef.list$patches))
    
    if (length(patchdiff) > 0) {
      if (patch.as.random == TRUE) {
        newdevs <- rnorm(length(patchdiff), mean = 0, sd = sd(coef.list$patches[,1]))
      } else {
        newdevs <- rep(0, (length(patchdiff)))
      }
      coef.list$patches <- rbind.data.frame(setNames(as.data.frame(as.matrix(newdevs, ncol = 1)), names(coef.list$patches)), coef.list$patches)
      patchlabels <- c(patchdiff, rownames(coef.list$patches)[(length(patchdiff) + 1):(length(coef.list$patches[,1]))])
      coef.list$patches <- as.data.frame(coef.list$patches[order(patchlabels),])
      rownames(coef.list$patches) <- sort(patchlabels)
      
    }
  }
  
  coef.list$years <- coef.list$years[,1]
  coef.list$patches <- coef.list$patches[,1]
  
  return(coef.list)
}

#' Extract Required Coefficient Values From glmer-estimated Vital Rate Models
#' 
#' \code{.modelextract.merMod()} extracts coefficient values from linear models 
#' estimated through the \code{\link[lme4]{glmer}()} function, to estimate
#' vital rates in \code{'lefko3'}. Used to supply coefficients to \code{\link{flefko3}()} and 
#' \code{\link{flefko2}()}.
#' 
#' @param model Model estimated through \code{\link[lme4]{glmer}()}.
#' @param paramnames Data frame giving the names of standard coefficients required
#' by matrix creation functions.
#' 
#' @return This function returns a list with the following elements:
#' \item{coefficients}{Vector of fixed effect coefficients.}
#' \item{years}{Vector of time step coefficients, typically random effect.}
#' \item{patches}{Vector of patch coefficients, typically random effect.}
#' \item{variances}{Residual variance terms associated with model.}
#' \item{family}{Distribution of response term.}
#' \item{sigma}{Standard deviation of response.}
#' 
#' @seealso \code{\link{.modelextract}()}
#' @seealso \code{\link{.modelextract.glmmTMB}()}
#' @seealso \code{\link{.modelextract.lmerMod}()}
#' @seealso \code{\link{.modelextract.glm}()}
#' @seealso \code{\link{.modelextract.lm}()}
#' @seealso \code{\link{.modelextract.numeric}()}
#' @seealso \code{\link{.modelextract.logical}()}
#' 
#' @keywords internal
#' @noRd
.modelextract.merMod <- function(model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random) {
  
  size2.coef <- 0
  size1.coef <- 0
  flw2.coef <- 0
  flw1.coef <- 0
  
  size1.size2.coef <- 0
  flw1.flw2.coef <- 0
  size1.flw1.coef <- 0
  size2.flw2.coef <- 0
  size1.flw2.coef <- 0
  size2.flw1.coef <- 0
  
  age.coef <- 0
  age.size1.coef <- 0
  age.size2.coef <- 0
  age.flw1.coef <- 0
  age.flw2.coef <- 0
  
  yintercept <- lme4::fixef(model)["(Intercept)"]
  
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "repst2")), 2]])) {flw2.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "repst2")), 2]]}
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "repst1")), 2]])) {flw1.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "repst1")), 2]]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])])) {flw1.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])])) {flw1.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])]}
  
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "size2")), 2]])) {size2.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "size2")), 2]]}
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "size1")), 2]])) {size1.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "size1")), 2]]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])])) {size1.size2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])])) {size1.size2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])])) {size1.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])])) {size2.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])])) {size1.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])])) {size2.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])])) {size2.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])])) {size1.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])])) {size2.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])])) {size1.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])]}
  
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "age")), 2]])) {age.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "age")), 2]]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])])) {age.size1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])])) {age.size1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])])) {age.size2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])])) {age.size2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])])) {age.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])])) {age.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])])) {age.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])])) {age.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])]}
  
  year.coefs <- lme4::ranef(model)[[paramnames[(which(paramnames$mainparams == "year2")), 2]]]
  
  if (!all(is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "patch")), 2]]))) {
    patch.coefs <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "patch")), 2]]
  } else if (!all(is.na(lme4::ranef(model)[[paramnames[(which(paramnames$mainparams == "patch")), 2]]]))) {
    patch.coefs <- lme4::ranef(model)[[paramnames[(which(paramnames$mainparams == "patch")), 2]]]
  } else {
    patch.coefs <- NA
  }
  
  coef.vec <- c(yintercept, flw1.coef, flw2.coef, size1.coef, size2.coef, flw1.flw2.coef, size1.size2.coef, 
                size1.flw1.coef, size2.flw2.coef, size2.flw1.coef, size1.flw2.coef, age.coef, age.size1.coef, 
                age.size2.coef, age.flw1.coef, age.flw2.coef)
  
  rvars <- as.data.frame(lme4::VarCorr(model))
  
  coef.list <- list(coefficients = coef.vec, years = year.coefs, patches = patch.coefs, 
                    variances = rvars, family = stats::family(model)$family, sigma = NA)
  
  #Here we correct for missing years
  yeardiff <- setdiff(mainyears, rownames(coef.list$years))
  
  if (length(yeardiff) > 0) {
    if (year.as.random == TRUE) {
      newdevs <- rnorm(length(yeardiff), mean = 0, sd = sd(coef.list$years[,1]))
    } else {
      newdevs <- rep(0, (length(yeardiff)))
    }
    coef.list$years <- rbind.data.frame(setNames(as.data.frame(as.matrix(newdevs, ncol = 1)), names(coef.list$years)), coef.list$years)
    yearlabels <- c(yeardiff, rownames(coef.list$years)[(length(yeardiff) + 1):(length(coef.list$years[,1]))])
    coef.list$years <- as.data.frame(coef.list$years[order(yearlabels),])
    rownames(coef.list$years) <- sort(yearlabels)
  }
  
  #This next part addresses missing patches
  if (class(coef.list$patches) == "data.frame") {
    patchdiff <- setdiff(mainpatches, rownames(coef.list$patches))
    
    if (length(patchdiff) > 0) {
      if (patch.as.random == TRUE) {
        newdevs <- rnorm(length(patchdiff), mean = 0, sd = sd(coef.list$patches[,1]))
      } else {
        newdevs <- rep(0, (length(patchdiff)))
      }
      coef.list$patches <- rbind.data.frame(setNames(as.data.frame(as.matrix(newdevs, ncol = 1)), names(coef.list$patches)), coef.list$patches)
      patchlabels <- c(patchdiff, rownames(coef.list$patches)[(length(patchdiff) + 1):(length(coef.list$patches[,1]))])
      coef.list$patches <- as.data.frame(coef.list$patches[order(patchlabels),])
      rownames(coef.list$patches) <- sort(patchlabels)
    }
  }
  
  if (class(coef.list$years) == "data.frame") {coef.list$years <- coef.list$years[,1]}
  if (class(coef.list$patches ) == "data.frame") {coef.list$patches <- coef.list$patches[,1]}
  
  return(coef.list)
}

#' Extract Required Coefficient Values From lmer-estimated Vital Rate Models
#' 
#' \code{.modelextract.lmerMod()} extracts coefficient values from linear models 
#' estimated through the \code{\link[lme4]{lmer}()} function, to estimate
#' vital rates in \code{'lefko3'}. Used to supply coefficients to \code{\link{flefko3}()} and 
#' \code{\link{flefko2}()}.
#' 
#' @param model Model estimated through \code{\link[lme4]{lmer}()}.
#' @param paramnames Data frame giving the names of standard coefficients required
#' by matrix creation functions.
#' 
#' @return This function returns a list with the following elements:
#' \item{coefficients}{Vector of fixed effect coefficients.}
#' \item{years}{Vector of time step coefficients, typically random effect.}
#' \item{patches}{Vector of patch coefficients, typically random effect.}
#' \item{variances}{Residual variance terms associated with model.}
#' \item{family}{Distribution of response term.}
#' \item{sigma}{Standard deviation of response.}
#' 
#' @seealso \code{\link{.modelextract}()}
#' @seealso \code{\link{.modelextract.glmmTMB}()}
#' @seealso \code{\link{.modelextract.merMod}()}
#' @seealso \code{\link{.modelextract.glm}()}
#' @seealso \code{\link{.modelextract.lm}()}
#' @seealso \code{\link{.modelextract.numeric}()}
#' @seealso \code{\link{.modelextract.logical}()}
#' 
#' @keywords internal
#' @noRd
.modelextract.lmerMod <- function(model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random) {
  
  size2.coef <- 0
  size1.coef <- 0
  flw2.coef <- 0
  flw1.coef <- 0
  
  size1.size2.coef <- 0
  flw1.flw2.coef <- 0
  size1.flw1.coef <- 0
  size2.flw2.coef <- 0
  size1.flw2.coef <- 0
  size2.flw1.coef <- 0
  
  age.coef <- 0
  age.size1.coef <- 0
  age.size2.coef <- 0
  age.flw1.coef <- 0
  age.flw2.coef <- 0
  
  yintercept <- lme4::fixef(model)["(Intercept)"]
  
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "repst2")), 2]])) {flw2.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "repst2")), 2]]}
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "repst1")), 2]])) {flw1.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "repst1")), 2]]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])])) {flw1.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])])) {flw1.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])]}
  
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "size2")), 2]])) {size2.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "size2")), 2]]}
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "size1")), 2]])) {size1.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "size1")), 2]]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])])) {size1.size2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])])) {size1.size2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])])) {size1.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])])) {size2.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])])) {size1.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])])) {size2.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])])) {size2.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])])) {size1.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])])) {size2.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])])) {size1.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])]}
  
  if (!is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "age")), 2]])) {age.coef <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "age")), 2]]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])])) {age.size1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])])) {age.size1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])])) {age.size2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])])) {age.size2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])])) {age.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])])) {age.flw1.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])]}
  
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])])) {age.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])]}
  if (!is.na(lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])])) {age.flw2.coef <- lme4::fixef(model)[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])]}
  
  year.coefs <- lme4::ranef(model)[[paramnames[(which(paramnames$mainparams == "year2")), 2]]]
  
  if (!all(is.na(lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "patch")), 2]]))) {
    patch.coefs <- lme4::fixef(model)[paramnames[(which(paramnames$mainparams == "patch")), 2]]
  } else if (!all(is.na(lme4::ranef(model)[[paramnames[(which(paramnames$mainparams == "patch")), 2]]]))) {
    patch.coefs <- lme4::ranef(model)[[paramnames[(which(paramnames$mainparams == "patch")), 2]]]
  } else {
    patch.coefs <- NA
  }
  
  sigmax <- attr(summary(model)$varcor, "sc")
  
  coef.vec <- c(yintercept, flw1.coef, flw2.coef, size1.coef, size2.coef, flw1.flw2.coef, size1.size2.coef, 
                size1.flw1.coef, size2.flw2.coef, size2.flw1.coef, size1.flw2.coef, age.coef, age.size1.coef, 
                age.size2.coef, age.flw1.coef, age.flw2.coef)
  
  rvars <- as.data.frame(lme4::VarCorr(model))
  
  coef.list <- list(coefficients = coef.vec, years = year.coefs, patches = patch.coefs, 
                    variances = rvars, family = "gaussian", sigma = sigmax)
  
  #Here we correct for missing years
  yeardiff <- setdiff(mainyears, rownames(coef.list$years))
  
  if (length(yeardiff) > 0) {
    if (year.as.random == TRUE) {
      newdevs <- rnorm(length(yeardiff), mean = 0, sd = sd(coef.list$years[,1]))
    } else {
      newdevs <- rep(0, (length(yeardiff)))
    }
    coef.list$years <- rbind.data.frame(setNames(as.data.frame(as.matrix(newdevs, ncol = 1)), names(coef.list$years)), coef.list$years)
    yearlabels <- c(yeardiff, rownames(coef.list$years)[(length(yeardiff) + 1):(length(coef.list$years[,1]))])
    coef.list$years <- as.data.frame(coef.list$years[order(yearlabels),])
    rownames(coef.list$years) <- sort(yearlabels)
  }
  
  #This next part addresses missing patches
  if (class(coef.list$patches) == "data.frame") {
    patchdiff <- setdiff(mainpatches, rownames(coef.list$patches))
    
    if (length(patchdiff) > 0) {
      if (patch.as.random == TRUE) {
        newdevs <- rnorm(length(patchdiff), mean = 0, sd = sd(coef.list$patches[,1]))
      } else {
        newdevs <- rep(0, (length(patchdiff)))
      }
      coef.list$patches <- rbind.data.frame(setNames(as.data.frame(as.matrix(newdevs, ncol = 1)), names(coef.list$patches)), coef.list$patches)
      patchlabels <- c(patchdiff, rownames(coef.list$patches)[(length(patchdiff) + 1):(length(coef.list$patches[,1]))])
      coef.list$patches <- as.data.frame(coef.list$patches[order(patchlabels),])
      rownames(coef.list$patches) <- sort(patchlabels)
    }
  }
  
  if (class(coef.list$years) == "data.frame") {coef.list$years <- coef.list$years[,1]}
  if (class(coef.list$patches ) == "data.frame") {coef.list$patches <- coef.list$patches[,1]}
  
  return(coef.list)
}

#' Extract Required Coefficient Values From glm-estimated Vital Rate Models
#' 
#' \code{.modelextract.glm()} extracts coefficient values from linear models 
#' estimated through the \code{\link[stats]{glm}()} function, to estimate
#' vital rates in \code{'lefko3'}. Used to supply coefficients to \code{\link{flefko3}()} and 
#' \code{\link{flefko2}()}.
#' 
#' @param model Model estimated through \code{\link[stats]{glm}()}.
#' @param paramnames Data frame giving the names of standard coefficients required
#' by matrix creation functions.
#' 
#' @return This function returns a list with the following elements:
#' \item{coefficients}{Vector of fixed effect coefficients.}
#' \item{years}{Vector of time step coefficients.}
#' \item{patches}{Vector of patch coefficients.}
#' \item{variances}{Residual variance terms associated with model.}
#' \item{family}{Distribution of response term.}
#' \item{sigma}{Standard deviation of response.}
#' 
#' @seealso \code{\link{.modelextract}()}
#' @seealso \code{\link{.modelextract.glmmTMB}()}
#' @seealso \code{\link{.modelextract.merMod}()}
#' @seealso \code{\link{.modelextract.lmerMod}()}
#' @seealso \code{\link{.modelextract.lm}()}
#' @seealso \code{\link{.modelextract.numeric}()}
#' @seealso \code{\link{.modelextract.logical}()}
#' 
#' @keywords internal
#' @noRd
.modelextract.glm <- function(model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random) {
  
  family.coef <- NA
  
  if (model$family["family"] == "binomial") {family.coef <- "binomial"}
  if (model$family["family"] == "poisson") {family.coef <- "poisson"}
  if (is.element("negbin", class(model))) {family.coef <- "negbin"}
  
  if (is.na(family.coef)) {stop("Model distribution not recognized.")}
  
  size2.coef <- 0
  size1.coef <- 0
  flw2.coef <- 0
  flw1.coef <- 0
  
  size1.size2.coef <- 0
  flw1.flw2.coef <- 0
  size1.flw1.coef <- 0
  size2.flw2.coef <- 0
  size1.flw2.coef <- 0
  size2.flw1.coef <- 0
  
  age.coef <- 0
  age.size1.coef <- 0
  age.size2.coef <- 0
  age.flw1.coef <- 0
  age.flw2.coef <- 0
  
  yintercept <- model$coefficients["(Intercept)"]
  
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "repst2")), 2]])) {flw2.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "repst2")), 2]]}
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "repst1")), 2]])) {flw1.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "repst1")), 2]]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])])) {flw1.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])])) {flw1.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])]}
  
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "size2")), 2]])) {size2.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "size2")), 2]]}
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "size1")), 2]])) {size1.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "size1")), 2]]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])])) {size1.size2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])])) {size1.size2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])])) {size1.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])])) {size2.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])])) {size1.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])])) {size2.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])])) {size1.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])])) {size2.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])])) {size1.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])])) {size2.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])]}
  
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "age")), 2]])) {age.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "age")), 2]]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])])) {age.size1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])])) {age.size1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])])) {age.size2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])])) {age.size2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])])) {age.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])])) {age.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])])) {age.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])])) {age.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])]}
  
  year.coef <- 0
  year.name.l.vec <- apply(as.matrix(names(model$coefficients)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "year2")), 2], X)})
  year.trace <- rep(0, length(year.name.l.vec))
  year.trace[which(year.name.l.vec == TRUE)] <- 1
  
  if (sum(year.trace) > 0) {
    year.vec <- rep(0, (length(which(year.name.l.vec == TRUE)) + 1))
    year.vec[2:length(year.vec)] <- model$coefficients[which(year.name.l.vec == TRUE)]
    year.name.vec <- c(paramnames[(which(paramnames$mainparams == "year2")), 2], names(model$coefficients)[which(apply(as.matrix(names(model$coefficients)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "year2")), 2], X)}) == TRUE)])
    year.name.vec <- apply(as.matrix(year.name.vec), 1, function(X) {gsub(paramnames[(which(paramnames$mainparams == "year2")), 2], "", X)})
    names(year.vec) <- year.name.vec
  } else {
    year.vec <- 0
    warning("Year was not included in the modeling. Consider setting year.as.random to FALSE.")
  }
  
  patch.vec <- NA
  patch.name.l.vec <- apply(as.matrix(names(model$coefficients)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "patch")), 2], X)})
  patch.trace <- rep(0, length(patch.name.l.vec))
  patch.trace[which(patch.name.l.vec == TRUE)] <- 1
  
  if (sum(patch.trace) > 0) {
    patch.vec <- rep(0, (length(which(patch.name.l.vec == TRUE)) + 1))
    patch.vec[2:length(patch.vec)] <- model$coefficients[which(patch.name.l.vec == TRUE)]
    patch.name.vec <- c(paramnames[(which(paramnames$mainparams == "patch")), 2], names(model$coefficients)[which(apply(as.matrix(names(model$coefficients)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "patch")), 2], X)}) == TRUE)])
    patch.name.vec <- apply(as.matrix(patch.name.vec), 1, function(X) {gsub(paramnames[(which(paramnames$mainparams == "patch")), 2], "", X)})
    patch.name.vec <- apply(as.matrix(patch.name.vec), 1, function(X) {gsub("as.factor\\()", "", X)})
    names(patch.vec) <- patch.name.vec
  }
  
  coef.vec <- c(yintercept, flw1.coef, flw2.coef, size1.coef, size2.coef, flw1.flw2.coef, size1.size2.coef, 
                size1.flw1.coef, size2.flw2.coef, size2.flw1.coef, size1.flw2.coef, age.coef, age.size1.coef, 
                age.size2.coef, age.flw1.coef, age.flw2.coef)
  
  coef.list <- list(coefficients = coef.vec, years = year.vec, patches = patch.vec,
                    variances = NA, family = family.coef, sigma = NA)
  
  #Here we correct for missing years
  yeardiff <- setdiff(mainyears, rownames(coef.list$years))
  
  if (length(yeardiff) > 0) {
    newdevs <- rep(0, (length(yeardiff)))
    coef.list$years <- rbind.data.frame(setNames(as.data.frame(as.matrix(newdevs, ncol = 1)), names(coef.list$years)), coef.list$years)
    yearlabels <- c(yeardiff, rownames(coef.list$years)[(length(yeardiff) + 1):(length(coef.list$years[,1]))])
    coef.list$years <- as.data.frame(coef.list$years[order(yearlabels),])
    rownames(coef.list$years) <- sort(yearlabels)
  }
  
  #This next part addresses missing patches
  if (class(coef.list$patches) == "data.frame") {
    patchdiff <- setdiff(mainpatches, rownames(coef.list$patches))
    
    if (length(patchdiff) > 0) {
      newdevs <- rep(0, (length(patchdiff)))
      coef.list$patches <- rbind.data.frame(setNames(as.data.frame(as.matrix(newdevs, ncol = 1)), names(coef.list$patches)), coef.list$patches)
      patchlabels <- c(patchdiff, rownames(coef.list$patches)[(length(patchdiff) + 1):(length(coef.list$patches[,1]))])
      coef.list$patches <- as.data.frame(coef.list$patches[order(patchlabels),])
      rownames(coef.list$patches) <- sort(patchlabels)
    }
  }
  
  if (class(coef.list$years) == "data.frame") {coef.list$years <- coef.list$years[,1]}
  if (class(coef.list$patches ) == "data.frame") {coef.list$patches <- coef.list$patches[,1]}
  
  return(coef.list)
}

#' Extract Required Coefficient Values From lm-estimated Vital Rate Models
#' 
#' \code{.modelextract.lm()} extracts coefficient values from linear models 
#' estimated through the \code{\link[stats]{lm}()} function, to estimate
#' vital rates in \code{'lefko3'}. Used to supply coefficients to \code{\link{flefko3}()} and 
#' \code{\link{flefko2}()}.
#' 
#' @param model Model estimated through \code{\link[stats]{lm}()}.
#' @param paramnames Data frame giving the names of standard coefficients required
#' by matrix creation functions.
#' 
#' @return This function returns a list with the following elements:
#' \item{coefficients}{Vector of fixed effect coefficients.}
#' \item{years}{Vector of time step coefficients.}
#' \item{patches}{Vector of patch coefficients.}
#' \item{variances}{Residual variance terms associated with model.}
#' \item{family}{Distribution of response term.}
#' \item{sigma}{Standard deviation of response.}
#' 
#' @seealso \code{\link{.modelextract}()}
#' @seealso \code{\link{.modelextract.glmmTMB}()}
#' @seealso \code{\link{.modelextract.merMod}()}
#' @seealso \code{\link{.modelextract.lmerMod}()}
#' @seealso \code{\link{.modelextract.glm}()}
#' @seealso \code{\link{.modelextract.numeric}()}
#' @seealso \code{\link{.modelextract.logical}()}
#' 
#' @keywords internal
#' @noRd
.modelextract.lm <- function(model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random) {
  
  size2.coef <- 0
  size1.coef <- 0
  flw2.coef <- 0
  flw1.coef <- 0
  
  size1.size2.coef <- 0
  flw1.flw2.coef <- 0
  size1.flw1.coef <- 0
  size2.flw2.coef <- 0
  size1.flw2.coef <- 0
  size2.flw1.coef <- 0
  
  age.coef <- 0
  age.size1.coef <- 0
  age.size2.coef <- 0
  age.flw1.coef <- 0
  age.flw2.coef <- 0
  
  yintercept <- model$coefficients["(Intercept)"]
  
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "repst2")), 2]])) {flw2.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "repst2")), 2]]}
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "repst1")), 2]])) {flw1.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "repst1")), 2]]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])])) {flw1.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])])) {flw1.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])]}
  
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "size2")), 2]])) {size2.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "size2")), 2]]}
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "size1")), 2]])) {size1.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "size1")), 2]]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])])) {size1.size2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])])) {size1.size2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])])) {size1.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])])) {size2.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])])) {size1.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])])) {size2.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])])) {size1.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])])) {size2.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])])) {size1.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])])) {size2.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])]}
  
  if (!is.na(model$coefficients[paramnames[(which(paramnames$mainparams == "age")), 2]])) {age.coef <- model$coefficients[paramnames[(which(paramnames$mainparams == "age")), 2]]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])])) {age.size1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "size1")), 2])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])])) {age.size1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size1")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])])) {age.size2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "size2")), 2])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])])) {age.size2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "size2")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])])) {age.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "repst1")), 2])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])])) {age.flw1.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst1")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])]}
  
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])])) {age.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "age")), 2], ":", paramnames[(which(paramnames$mainparams == "repst2")), 2])]}
  if (!is.na(model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])])) {age.flw2.coef <- model$coefficients[paste0(paramnames[(which(paramnames$mainparams == "repst2")), 2], ":", paramnames[(which(paramnames$mainparams == "age")), 2])]}
  
  year.coef <- 0
  year.name.l.vec <- apply(as.matrix(names(model$coefficients)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "year2")), 2], X)})
  year.trace <- rep(0, length(year.name.l.vec))
  year.trace[which(year.name.l.vec == TRUE)] <- 1
  
  if (sum(year.trace) > 0) {
    year.vec <- rep(0, (length(which(year.name.l.vec == TRUE)) + 1))
    year.vec[2:length(year.vec)] <- model$coefficients[which(year.name.l.vec == TRUE)]
    year.name.vec <- c(paramnames[(which(paramnames$mainparams == "year2")), 2], names(model$coefficients)[which(apply(as.matrix(names(model$coefficients)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "year2")), 2], X)}) == TRUE)])
    year.name.vec <- apply(as.matrix(year.name.vec), 1, function(X) {gsub(paramnames[(which(paramnames$mainparams == "year2")), 2], "", X)})
    names(year.vec) <- year.name.vec
  } else {
    year.vec <- 0
    warning("Year was not included in the modeling. Consider setting year.as.random to FALSE.")
  }
  
  patch.vec <- NA
  patch.name.l.vec <- apply(as.matrix(names(model$coefficients)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "patch")), 2], X)})
  patch.trace <- rep(0, length(patch.name.l.vec))
  patch.trace[which(patch.name.l.vec == TRUE)] <- 1
  
  if (sum(patch.trace) > 0) {
    patch.vec <- rep(0, (length(which(patch.name.l.vec == TRUE)) + 1))
    patch.vec[2:length(patch.vec)] <- model$coefficients[which(patch.name.l.vec == TRUE)]
    patch.name.vec <- c(paramnames[(which(paramnames$mainparams == "patch")), 2], names(model$coefficients)[which(apply(as.matrix(names(model$coefficients)), 1, function(X) {grepl(paramnames[(which(paramnames$mainparams == "patch")), 2], X)}) == TRUE)])
    patch.name.vec <- apply(as.matrix(patch.name.vec), 1, function(X) {gsub(paramnames[(which(paramnames$mainparams == "patch")), 2], "", X)})
    patch.name.vec <- apply(as.matrix(patch.name.vec), 1, function(X) {gsub("as.factor\\()", "", X)})
    names(patch.vec) <- patch.name.vec
  }
  
  sigmax <- stats::sigma(model)
  
  coef.vec <- c(yintercept, flw1.coef, flw2.coef, size1.coef, size2.coef, flw1.flw2.coef, size1.size2.coef, 
                size1.flw1.coef, size2.flw2.coef, size2.flw1.coef, size1.flw2.coef, age.coef, age.size1.coef, 
                age.size2.coef, age.flw1.coef, age.flw2.coef)
  
  coef.list <- list(coefficients = coef.vec, years = year.vec, patches = patch.vec,
                    variances = NA, family = "gaussian", sigma = sigmax)
  
  #Here we correct for missing years
  yeardiff <- setdiff(mainyears, rownames(coef.list$years))
  
  if (length(yeardiff) > 0) {
    newdevs <- rep(0, (length(yeardiff)))
    coef.list$years <- rbind.data.frame(setNames(as.data.frame(as.matrix(newdevs, ncol = 1)), names(coef.list$years)), coef.list$years)
    yearlabels <- c(yeardiff, rownames(coef.list$years)[(length(yeardiff) + 1):(length(coef.list$years[,1]))])
    coef.list$years <- as.data.frame(coef.list$years[order(yearlabels),])
    rownames(coef.list$years) <- sort(yearlabels)
  }
  
  #This next part addresses missing patches
  if (class(coef.list$patches) == "data.frame") {
    patchdiff <- setdiff(mainpatches, rownames(coef.list$patches))
    
    if (length(patchdiff) > 0) {
      newdevs <- rep(0, (length(patchdiff)))
      coef.list$patches <- rbind.data.frame(setNames(as.data.frame(as.matrix(newdevs, ncol = 1)), names(coef.list$patches)), coef.list$patches)
      patchlabels <- c(patchdiff, rownames(coef.list$patches)[(length(patchdiff) + 1):(length(coef.list$patches[,1]))])
      coef.list$patches <- as.data.frame(coef.list$patches[order(patchlabels),])
      rownames(coef.list$patches) <- sort(patchlabels)
    }
  }
  
  if (class(coef.list$years) == "data.frame") {coef.list$years <- coef.list$years[,1]}
  if (class(coef.list$patches ) == "data.frame") {coef.list$patches <- coef.list$patches[,1]}
  
  return(coef.list)
}

#' Return Simple List When Model Is Scalar Numeric
#' 
#' \code{.modelextract.numeric()} returns NA when a vital rate model is simply a
#' scalar numeric value. Used to supply coefficients to \code{\link{flefko3}()} and 
#' \code{\link{flefko2}()}.
#' 
#' @param model Model estimated through functions supported by \code{'lefko3'}.
#' @param paramnames Data frame giving the names of standard coefficients required
#' by matrix creation functions.
#' 
#' @return This function returns a list with single values for matrix estimation.
#' 
#' @seealso \code{\link{.modelextract}()}
#' @seealso \code{\link{.modelextract.glmmTMB}()}
#' @seealso \code{\link{.modelextract.merMod}()}
#' @seealso \code{\link{.modelextract.lmerMod}()}
#' @seealso \code{\link{.modelextract.glm}()}
#' @seealso \code{\link{.modelextract.lm}()}
#' @seealso \code{\link{.modelextract.logical}()}
#' 
#' @keywords internal
#' @noRd
.modelextract.numeric <- function(model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random) {
  coef.list <- list(coefficients = c(model), years = c(0), patches = c(0), 
                    variances = as.data.frame(0), family = 1, sigma = 1)
  
  return(coef.list)
}

#' Return Simple List When Model Is Logical
#' 
#' \code{.modelextract.logical()} returns NA when a vital rate model is simply a
#' logical value. Used to supply coefficients to \code{\link{flefko3}()} and 
#' \code{\link{flefko2}()}.
#' 
#' @param model Model estimated through functions supported by \code{'lefko3'}.
#' @param paramnames Data frame giving the names of standard coefficients required
#' by matrix creation functions.
#' 
#' @return This function returns NA.
#' 
#' @seealso \code{\link{.modelextract}()}
#' @seealso \code{\link{.modelextract.glmmTMB}()}
#' @seealso \code{\link{.modelextract.merMod}()}
#' @seealso \code{\link{.modelextract.lmerMod}()}
#' @seealso \code{\link{.modelextract.glm}()}
#' @seealso \code{\link{.modelextract.lm}()}
#' @seealso \code{\link{.modelextract.numeric}()}
#' 
#' @keywords internal
#' @noRd
.modelextract.logical <- function(model, paramnames, mainyears, mainpatches, year.as.random, patch.as.random) {
  coef.list <- list(coefficients = c(0), years = c(0), patches = c(0), 
                    variances = as.data.frame(0), family = 1, sigma = 1)
  
  return(coef.list)
}
