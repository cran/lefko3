#' Develop Best-fit Vital Rate Estimation Models for MPM Development
#' 
#' Function \code{modelsearch()} runs exhaustive model building and selection
#' for each vital rate needed to estimate a function-based MPM or IPM. It
#' returns best-fit models for each vital rate, model table showing all models
#' tested, and model quality control data. The final output can be used as input
#' in other functions within this package.
#' 
#' @param data The vertical dataset to be used for analysis. This dataset should 
#' be of class \code{hfvdata}, but can also be a data frame formatted similarly
#' to the output format provided by functions \code{\link{verticalize3}()} or
#' \code{\link{historicalize3}()}, as long as all needed variables are properly
#' designated.
#' @param stageframe The stageframe characterizing the life history model used.
#' Optional unless \code{test.group = TRUE}, in which case it is required.
#' Defaults to \code{NULL}.
#' @param historical A logical variable denoting whether to assess the effects
#' of state in occasion \emph{t}-1, in addition to state in occasion \emph{t}.
#' Defaults to \code{TRUE}.
#' @param approach The statistical approach to be taken for model building. The 
#' default is \code{"mixed"}, which uses the mixed model approach utilized in 
#' packages \code{lme4} and \code{glmmTMB}. Other options include \code{"glm"},
#' which uses generalized linear modeling assuming that all factors are fixed.
#' @param suite This describes the global model for each vital rate estimation,
#' and has the following possible values: \code{full}, includes main effects and
#' all two-way interactions of size and reproductive status; \code{main},
#' includes main effects only of size and reproductive status; \code{size},
#' includes only size (also interactions between size in historical model);
#' \code{rep}, includes only reproductive status (also interactions between
#' status in historical model); and \code{cons}, all vital rates estimated only
#' as y-intercepts. If \code{approach = "glm"} and
#' \code{year.as.random = FALSE}, then year is also included as a fixed effect,
#' and, in the case of \code{full}, included in two-way interactions. Defaults
#' to \code{size}.
#' @param bestfit A variable indicating the model selection criterion for the
#' choice of best-fit model. The default is \code{AICc&k}, which chooses the 
#' best-fit model as the model with the lowest AICc or, if not the same model,
#' then the model that has the lowest degrees of freedom among models with
#' \eqn{\Delta AICc <= 2.0}. Alternatively, \code{AICc} may be chosen, in which
#' case the best-fit model is simply the model with the lowest AICc value.
#' @param vitalrates A vector describing which vital rates will be estimated via
#' linear modeling, with the following options: \code{surv}, survival
#' probability; \code{obs}, observation probability; \code{size}, overall size;
#' \code{repst}, probability of reproducing; and \code{fec}, amount of
#' reproduction (overall fecundity). Defaults to \code{c("surv", "size", "fec")}.
#' @param surv A vector indicating the variable names coding for status as alive
#' or dead in occasions \emph{t}+1, \emph{t}, and \emph{t}-1, respectively.
#' Defaults to \code{c("alive3", "alive2", "alive1")}.
#' @param obs A vector indicating the variable names coding for observation
#' status in occasions \emph{t}+1, \emph{t}, and \emph{t}-1, respectively.
#' Defaults to \code{c("obsstatus3", "obsstatus2", "obsstatus1")}.
#' @param size A vector indicating the variable names coding for the primary
#' size variable on occasions \emph{t}+1, \emph{t}, and \emph{t}-1,
#' respectively. Defaults to \code{c("sizea3", "sizea2", "sizea1")}.
#' @param sizeb A vector indicating the variable names coding for the secondary
#' size variable on occasions \emph{t}+1, \emph{t}, and \emph{t}-1,
#' respectively. Defaults to \code{c(NA, NA, NA)}, in which case \code{sizeb} is
#' not used.
#' @param sizec A vector indicating the variable names coding for the tertiary
#' size variable on occasions \emph{t}+1, \emph{t}, and \emph{t}-1,
#' respectively. Defaults to \code{c(NA, NA, NA)}, in which case \code{sizec} is
#' not used.
#' @param repst A vector indicating the variable names coding for reproductive
#' status in occasions \emph{t}+1, \emph{t}, and \emph{t}-1, respectively.
#' Defaults to \code{c("repstatus3", "repstatus2", "repstatus1")}.
#' @param fec A vector indicating the variable names coding for fecundity in
#' occasions \emph{t}+1, \emph{t}, and \emph{t}-1, respectively. Defaults to
#' \code{c("feca3", "feca2", "feca1")}.
#' @param stage A vector indicating the variable names coding for stage in
#' occasions \emph{t}+1, \emph{t}, and \emph{t}-1. Defaults to
#' \code{c("stage3", "stage2", "stage1")}.
#' @param indiv A text value indicating the variable name coding individual
#' identity. Defaults to \code{"individ"}.
#' @param patch A text value indicating the variable name coding for patch,
#' where patches are defined as permanent subgroups within the study population.
#' Defaults to \code{NA}.
#' @param year A text value indicating the variable coding for observation
#' occasion \emph{t}. Defaults to \code{year2}.
#' @param density A text value indicating the name of the variable coding for
#' spatial density, should the user wish to test spatial density as a fixed
#' factor affecting vital rates. Defaults to \code{NA}.
#' @param sizedist The probability distribution used to model primary size.
#' Options include \code{"gaussian"} for the Normal distribution (default),
#' \code{"poisson"} for the Poisson distribution, \code{"negbin"} for the
#' negative binomial distribution (quadratic parameterization), and 
#' \code{"gamma"} for the Gamma distribution.
#' @param sizebdist The probability distribution used to model secondary size.
#' Options include \code{"gaussian"} for the Normal distribution,
#' \code{"poisson"} for the Poisson distribution, \code{"negbin"} for the
#' negative binomial distribution (quadratic parameterization), and
#' \code{"gamma"} for the Gamma distribution. Defaults to \code{NA}.
#' @param sizecdist The probability distribution used to model tertiary size.
#' Options include \code{"gaussian"} for the Normal distribution,
#' \code{"poisson"} for the Poisson distribution, \code{"negbin"} for the
#' negative binomial distribution (quadratic parameterization), and
#' \code{"gamma"} for the Gamma distribution. Defaults to \code{NA}.
#' @param fecdist The probability distribution used to model fecundity. Options
#' include \code{"gaussian"} for the Normal distribution (default),
#' \code{"poisson"} for the Poisson distribution, \code{"negbin"} for the
#' negative binomial distribution (quadratic parameterization), and
#' \code{"gamma"} for the Gamma distribution.
#' @param size.zero A logical variable indicating whether the primary size
#' distribution should be zero-inflated. Only applies to Poisson and negative
#' binomial distributions. Defaults to \code{FALSE}.
#' @param sizeb.zero A logical variable indicating whether the secondary size
#' distribution should be zero-inflated. Only applies to Poisson and negative
#' binomial distributions. Defaults to \code{FALSE}.
#' @param sizec.zero  A logical variable indicating whether the tertiary size
#' distribution should be zero-inflated. Only applies to Poisson and negative
#' binomial distributions. Defaults to \code{FALSE}.
#' @param size.trunc A logical variable indicating whether the primary size
#' distribution should be zero-truncated. Only applies to Poisson and negative
#' binomial distributions. Defaults to \code{FALSE}. Cannot be \code{TRUE} if
#' \code{size.zero = TRUE}.
#' @param sizeb.trunc A logical variable indicating whether the secondary size
#' distribution should be zero-truncated. Only applies to Poisson and negative
#' binomial distributions. Defaults to \code{FALSE}. Cannot be \code{TRUE} if
#' \code{sizeb.zero = TRUE}.
#' @param sizec.trunc A logical variable indicating whether the tertiary size
#' distribution should be zero-truncated. Only applies to Poisson and negative
#' binomial distributions. Defaults to \code{FALSE}. Cannot be \code{TRUE} if
#' \code{sizec.zero = TRUE}.
#' @param fec.zero A logical variable indicating whether the fecundity
#' distribution should be zero-inflated. Only applies to Poisson and negative
#' binomial distributions. Defaults to \code{FALSE}.
#' @param fec.trunc A logical variable indicating whether the fecundity
#' distribution should be zero-truncated. Only applies to the Poisson and
#' negative binomial distributions. Defaults to \code{FALSE}. Cannot be
#' \code{TRUE} if \code{fec.zero = TRUE}.
#' @param patch.as.random If set to \code{TRUE} and \code{approach = "mixed"},
#' then \code{patch} is included as a random factor. If set to \code{FALSE} and
#' \code{approach = "glm"}, then \code{patch} is included as a fixed factor. All
#' other combinations of logical value and \code{approach} lead to \code{patch}
#' not being included in modeling. Defaults to \code{TRUE}.
#' @param year.as.random If set to \code{TRUE} and \code{approach = "mixed"},
#' then \code{year} is included as a random factor. If set to \code{FALSE}, then
#' \code{year} is included as a fixed factor. All other combinations of logical
#' value and \code{approach} lead to \code{year} not being included in modeling.
#' Defaults to \code{TRUE}.
#' @param juvestimate An optional variable denoting the stage name of the
#' juvenile stage in the vertical dataset. If not \code{NA}, and \code{stage} is
#' also given (see below), then vital rates listed in \code{vitalrates} other
#' than \code{fec} will also be estimated from the juvenile stage to all adult
#' stages. Defaults to \code{NA}, in which case juvenile vital rates are not
#' estimated.
#' @param juvsize A logical variable denoting whether size should be used as a
#' term in models involving transition from the juvenile stage. Defaults to
#' \code{FALSE}, and is only used if \code{juvestimate} does not equal
#' \code{NA}.
#' @param jsize.zero A logical variable indicating whether the primary size
#' distribution of juveniles should be zero-inflated. Only applies to Poisson
#' and negative binomial distributions. Defaults to \code{FALSE}.
#' @param jsizeb.zero A logical variable indicating whether the secondary size
#' distribution of juveniles should be zero-inflated. Only applies to Poisson
#' and negative binomial distributions. Defaults to \code{FALSE}.
#' @param jsizec.zero A logical variable indicating whether the tertiary size
#' distribution of juveniles should be zero-inflated. Only applies to Poisson
#' and negative binomial distributions. Defaults to \code{FALSE}.
#' @param jsize.trunc A logical variable indicating whether the primary size
#' distribution in juveniles should be zero-truncated. Defaults to \code{FALSE}.
#' Cannot be \code{TRUE} if \code{jsize.zero = TRUE}.
#' @param jsizeb.trunc A logical variable indicating whether the secondary size
#' distribution in juveniles should be zero-truncated. Defaults to \code{FALSE}.
#' Cannot be \code{TRUE} if \code{jsizeb.zero = TRUE}.
#' @param jsizec.trunc A logical variable indicating whether the tertiary size
#' distribution in juveniles should be zero-truncated. Defaults to \code{FALSE}.
#' Cannot be \code{TRUE} if \code{jsizec.zero = TRUE}.
#' @param fectime A variable indicating which year of fecundity to use as the
#' response term in fecundity models. Options include \code{2}, which refers to
#' occasion \emph{t}, and \code{3}, which refers to occasion \emph{t}+1.
#' Defaults to \code{2}.
#' @param censor A vector denoting the names of censoring variables in the
#' dataset, in order from occasion \emph{t}+1, followed by occasion \emph{t},
#' and lastly followed by occasion \emph{t}-1. Defaults to \code{NA}.
#' @param age Designates the name of the variable corresponding to age in the
#' vertical dataset. Defaults to \code{NA}, in which case age is not included in
#' linear models. Should only be used if building Leslie or age x stage
#' matrices.
#' @param indcova Vector designating the names in occasions \emph{t}+1,
#' \emph{t}, and \emph{t}-1 of an individual covariate. Defaults to \code{NA}.
#' @param indcovb Vector designating the names in occasions \emph{t}+1,
#' \emph{t}, and \emph{t}-1 of a second individual covariate. Defaults to
#' \code{NA}.
#' @param indcovc Vector designating the names in occasions \emph{t}+1,
#' \emph{t}, and \emph{t}-1 of a third individual covariate. Defaults to
#' \code{NA}.
#' @param random.indcova A logical value indicating whether \code{indcova}
#' should be treated as a random categorical factor, rather than as a fixed
#' factor. Defaults to \code{FALSE}.
#' @param random.indcovb A logical value indicating whether \code{indcovb}
#' should be treated as a random categorical factor, rather than as a fixed
#' factor. Defaults to \code{FALSE}.
#' @param random.indcovc A logical value indicating whether \code{indcovc}
#' should be treated as a random categorical factor, rather than as a fixed
#' factor. Defaults to \code{FALSE}.
#' @param test.group A logical value indicating whether to include the
#' \code{group} variable from the input \code{stageframe} as a fixed categorical
#' variable in linear models. Defaults to \code{FALSE}.
#' @param show.model.tables If set to TRUE, then includes full modeling tables
#' in the output. Defaults to \code{TRUE}.
#' @param global.only If set to TRUE, then only global models will be built and
#' evaluated. Defaults to \code{FALSE}.
#' @param quiet If set to TRUE, then model building and selection will proceed
#' with most warnings and diagnostic messages silenced. Defaults to
#' \code{FALSE}.
#' 
#' @return This function yields an object of class \code{lefkoMod}, which is a
#' list in which the first 13 elements are the best-fit models for survival,
#' observation status, primary size, secondary size, tertiary size,
#' reproductive status, fecundity, juvenile survival, juvenile observation,
#' juvenile primary size, juvenile secondary size, juvenile tertiary size, and
#' juvenile transition to reproduction, respectively, followed by 13 elements
#' corresponding to the model tables for each of these vital rates, in order,
#' followed by a data frame showing the order and names of variables used in
#' modeling, followed by a single character element denoting the criterion used
#' for model selection, and ending on a data frame with quality control data:
#' 
#' \item{survival_model}{Best-fit model of the binomial probability of survival
#' from occasion \emph{t} to occasion \emph{t}+1. Defaults to \code{1}.}
#' \item{observation_model}{Best-fit model of the binomial probability of 
#' observation in occasion \emph{t}+1 given survival to that occasion. Defaults
#' to \code{1}.}
#' \item{size_model}{Best-fit model of the primary size metric on occasion
#' \emph{t}+1 given survival to and observation in that occasion. Defaults to
#' \code{1}.}
#' \item{sizeb_model}{Best-fit model of the secondary size metric on occasion
#' \emph{t}+1 given survival to and observation in that occasion. Defaults to
#' \code{1}.}
#' \item{sizec_model}{Best-fit model of the tertiary size metric on occasion
#' \emph{t}+1 given survival to and observation in that occasion. Defaults to
#' \code{1}.}
#' \item{repstatus_model}{Best-fit model of the binomial probability of
#' reproduction in occasion \emph{t}+1, given survival to and observation in
#' that occasion. Defaults to \code{1}.}
#' \item{fecundity_model}{Best-fit model of fecundity in occasion \emph{t}+1
#' given survival to, and observation and reproduction in that occasion.
#' Defaults to \code{1}.}
#' \item{juv_survival_model}{Best-fit model of the binomial probability of
#' survival from occasion \emph{t} to occasion \emph{t}+1 of an immature
#' individual. Defaults to \code{1}.}
#' \item{juv_observation_model}{Best-fit model of the binomial probability of 
#' observation in occasion \emph{t}+1 given survival to that occasion of an
#' immature individual. Defaults to \code{1}.}
#' \item{juv_size_model}{Best-fit model of the primary size metric on occasion
#' \emph{t}+1 given survival to and observation in that occasion of an immature
#' individual. Defaults to \code{1}.}
#' \item{juv_sizeb_model}{Best-fit model of the secondary size metric on
#' occasion \emph{t}+1 given survival to and observation in that occasion of an
#' immature individual. Defaults to \code{1}.}
#' \item{juv_sizec_model}{Best-fit model of the tertiary size metric on occasion
#' \emph{t}+1 given survival to and observation in that occasion of an immature
#' individual. Defaults to \code{1}.}
#' \item{juv_reproduction_model}{Best-fit model of the binomial probability of
#' reproduction in occasion \emph{t}+1, given survival to and observation in
#' that occasion of an individual that was immature in occasion \emph{t}. This
#' model is technically not a model of reproduction probability for individuals
#' that are immature, rather reproduction probability here is given for
#' individuals that are mature in occasion \emph{t}+1 but immature in occasion
#' \emph{t}. Defaults to \code{1}.}
#' \item{survival_table}{Full dredge model table of survival probability.}
#' \item{observation_table}{Full dredge model table of observation probability.}
#' \item{size_table}{Full dredge model table of the primary size variable.}
#' \item{sizeb_table}{Full dredge model table of the secondary size variable.}
#' \item{sizec_table}{Full dredge model table of the tertiary size variable.}
#' \item{repstatus_table}{Full dredge model table of reproduction probability.}
#' \item{fecundity_table}{Full dredge model table of fecundity.}
#' \item{juv_survival_table}{Full dredge model table of immature survival 
#' probability.}
#' \item{juv_observation_table}{Full dredge model table of immature observation
#' probability.}
#' \item{juv_size_table}{Full dredge model table of primary size in immature
#' individuals.}
#' \item{juv_sizeb_table}{Full dredge model table of secondary size in immature
#' individuals.}
#' \item{juv_sizec_table}{Full dredge model table of tertiary size in immature
#' individuals.}
#' \item{juv_reproduction_table}{Full dredge model table of immature
#' reproduction probability.}
#' \item{criterion}{Character variable denoting the criterion used to determine
#' the best-fit model.}
#' \item{qc}{Data frame with four variables: 1) Name of vital rate, 2) number
#' of individuals used to model that vital rate, 3) number of individual
#' transitions used to model that vital rate, and 4) accuracy of model expressed
#' as percent of predicted responses equal to actual responses (only in binomial
#' models).}
#' 
#' @section Notes:
#' The mechanics governing model building are fairly robust to errors and
#' exceptions. The function attempts to build global models, and simplifies
#' models automatically should model building fail. Model building proceeds
#' through the functions \code{\link[stats]{lm}()} (GLM with Gaussian response),
#' \code{\link[stats]{glm}()} (GLM with Poisson, Gamma, or binomial response),
#' \code{\link[MASS]{glm.nb}()} (GLM with negative binomial response),
#' \code{\link[pscl]{zeroinfl}()} (GLM with zero-inflated Poisson or negative
#' binomial response), \code{\link[VGAM]{vglm}()} (GLM with zero-truncated
#' Poisson or negative binomial response), \code{\link[lme4]{lmer}()} (mixed
#' model with Gaussian response), \code{\link[lme4]{glmer}()} (mixed model with
#' binomial, Poisson, or Gamma response), and \code{\link[glmmTMB]{glmmTMB}()}
#' (mixed model with negative binomial, or zero-truncated or zero-inflated
#' Poisson or negative binomial response). See documentation related to these
#' functions for further information. Any response term that is invariable in
#' the dataset will lead to a best-fit model for that response represented by a
#' single constant value.
#' 
#' Exhaustive model building and selection proceeds via the
#' \code{\link[MuMIn]{dredge}()} function in package \code{MuMIn}. This function
#' is verbose, so that any errors and warnings developed during model building,
#' model analysis, and model selection can be found and dealt with.
#' Interpretations of errors during global model analysis may be found in
#' documentation for the functions and packages mentioned. Package \code{MuMIn}
#' is used for model dredging (see \link[MuMIn]{dredge}()), and errors and
#' warnings during dredging can be interpreted using the documentation for that
#' package. Errors occurring during dredging lead to the adoption of the global
#' model as the best-fit, and the user should view all logged errors and
#' warnings to determine the best way to proceed. The \code{quiet = TRUE} option
#' can be used to silence dredge warnings, but users should note that automated
#' model selection can be viewed as a black box, and so care should be taken to
#' ensure that the models run make biological sense, and that model quality is
#' prioritized.
#' 
#' Exhaustive model selection through dredging works best with larger datasets
#' and fewer tested parameters. Setting \code{suite = "full"} may initiate a
#' dredge that takes a dramatically long time, particularly if the model is
#' historical, individual covariates are used, or a zero-inflated distribution
#' is assumed. In such cases, the number of models built and tested will run at
#' least in the millions. Small datasets will also increase the error associated
#' with these tests, leading to adoption of simpler models overall.
#'
#' Care must be taken to build models that test the impacts of state in occasion
#' \emph{t}-1 for historical models, and that do not test these impacts for
#' ahistorical models. Ahistorical matrix modeling particularly will yield
#' biased transition estimates if historical terms from models are ignored. This
#' can be dealt with at the start of modeling by setting 
#' \code{historical = FALSE} for the ahistorical case, and 
#' \code{historical = TRUE} for the historical case.
#' 
#' This function handles generalized linear models (GLMs) under zero-inflated
#' distributions using the \code{\link[pscl]{zeroinfl}()} function, and zero-
#' truncated distributions using the \code{\link[VGAM]{vglm}()} function. Model
#' dredging may fail with these functions, leading to the global model being
#' accepted as the best-fit model. However, model dredges of mixed models work
#' for all distributions. We encourage the use of mixed models in all cases.
#' 
#' The negative binomial and truncated negative binomial distributions use the
#' quadratic structure emphasized in Hardin and Hilbe (2018, 4th Edition of
#' Generalized Linear Models and Extensions). The truncated negative binomial
#' distribution may fail to predict size probabilities correctly when dispersion
#' is near that expected of the Poisson distribution. To prevent this problem,
#' we have integrated a cap on the overdispersion parameter. However, when using
#' this distribution, please check the matrix column sums to make sure that they
#' do not predict survival greater than 1.0. If they do, then please use either
#' the negative binomial distribution or the zero-truncated Poisson
#' distribution.
#' 
#' If density dependence is explored through function \code{modelsearch()},
#' then the interpretation of density is not the full population size but rather
#' the spatial density term included in the dataset.
#' 
#' @examples
#' \donttest{
#' # Lathyrus example
#' data(lathyrus)
#' 
#' sizevector <- c(0, 4.6, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8,
#'   9)
#' stagevector <- c("Sd", "Sdl", "Dorm", "Sz1nr", "Sz2nr", "Sz3nr", "Sz4nr",
#'   "Sz5nr", "Sz6nr", "Sz7nr", "Sz8nr", "Sz9nr", "Sz1r", "Sz2r", "Sz3r", 
#'   "Sz4r", "Sz5r", "Sz6r", "Sz7r", "Sz8r", "Sz9r")
#' repvector <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
#'   0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 4.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 
#'   0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
#' 
#' lathframeln <- sf_create(sizes = sizevector, stagenames = stagevector, 
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector, 
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec, 
#'   propstatus = propvector)
#' 
#' lathvertln <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
#'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9, 
#'   juvcol = "Seedling1988", sizeacol = "lnVol88", repstracol = "Intactseed88",
#'   fecacol = "Intactseed88", deadacol = "Dead1988", 
#'   nonobsacol = "Dormant1988", stageassign = lathframeln, stagesize = "sizea",
#'   censorcol = "Missing1988", censorkeep = NA, NAas0 = TRUE, censor = TRUE)
#' 
#' lathvertln$feca2 <- round(lathvertln$feca2)
#' lathvertln$feca1 <- round(lathvertln$feca1)
#' lathvertln$feca3 <- round(lathvertln$feca3)
#' 
#' lathmodelsln3 <- modelsearch(lathvertln, historical = TRUE, 
#'   approach = "mixed", suite = "main", 
#'   vitalrates = c("surv", "obs", "size", "repst", "fec"), juvestimate = "Sdl",
#'   bestfit = "AICc&k", sizedist = "gaussian", fecdist = "poisson", 
#'   indiv = "individ", patch = "patchid", year = "year2",year.as.random = TRUE,
#'   patch.as.random = TRUE, show.model.tables = TRUE, quiet = TRUE)
#' 
#' # Here we use supplemental() to provide overwrite and reproductive info
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "mat", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "Sd", "Sdl", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "Sd", "rep", "Sd", "mat", "mat"),
#'   eststage3 = c(NA, NA, NA, NA, "mat", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, "Sdl", NA, NA),
#'   eststage1 = c(NA, NA, NA, NA, "Sdl", NA, NA),
#'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 1, 1, 3, 3), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
#'   stageframe = lathframeln, historical = TRUE)
#' 
#' lathmat3ln <- flefko3(year = "all", patch = "all", stageframe = lathframeln, 
#'   modelsuite = lathmodelsln3, data = lathvertln, supplement = lathsupp3, 
#'   patchcol = "patchid", yearcol = "year2", year.as.random = TRUE,
#'   patch.as.random = TRUE, reduce = FALSE)
#' 
#' summary(lathmat3ln)
#' }
#' 
#' @export
modelsearch <- function(data, stageframe = NULL, historical = TRUE,
  approach = "mixed", suite = "size", bestfit = "AICc&k",
  vitalrates = c("surv", "size", "fec"), surv = c("alive3", "alive2", "alive1"),
  obs = c("obsstatus3", "obsstatus2", "obsstatus1"),
  size = c("sizea3", "sizea2", "sizea1"), sizeb = c(NA, NA, NA), 
  sizec = c(NA, NA, NA), repst = c("repstatus3", "repstatus2", "repstatus1"),
  fec = c("feca3", "feca2", "feca1"), stage = c("stage3", "stage2", "stage1"),
  indiv = "individ", patch = NA, year = "year2", density = NA,
  sizedist = "gaussian", sizebdist = NA, sizecdist = NA, fecdist = "gaussian",
  size.zero = FALSE, sizeb.zero = FALSE, sizec.zero = FALSE, size.trunc = FALSE,
  sizeb.trunc = FALSE, sizec.trunc = FALSE, fec.zero = FALSE, fec.trunc = FALSE,
  patch.as.random = TRUE, year.as.random = TRUE, juvestimate = NA,
  juvsize = FALSE, jsize.zero = FALSE, jsizeb.zero = FALSE, jsizec.zero = FALSE,
  jsize.trunc = FALSE, jsizeb.trunc = FALSE, jsizec.trunc = FALSE, fectime = 2,
  censor = NA, age = NA, indcova = NA, indcovb = NA, indcovc = NA,
  random.indcova = FALSE, random.indcovb = FALSE, random.indcovc = FALSE,
  test.group = FALSE, show.model.tables = TRUE, global.only = FALSE,
  quiet = FALSE) {
  
  censor1 <- censor2 <- censor3 <- surv.data <- obs.data <- size.data <- NULL
  repst.data <- fec.data <- juvsurv.data <- juvobs.data <- NULL
  juvsize.data <- juvrepst.data <- usedfec <- NULL
  
  sizeb_used <- FALSE
  sizec_used <- FALSE
  density_used <- FALSE
  indcova_used <- FALSE
  indcovb_used <- FALSE
  indcovc_used <- FALSE
  
  #Input testing, input standardization, and exception handling
  if (all(class(data) != "hfvdata")) {
    warning("This function was made to work with standardized historically
      formatted vertical datasets, as provided by the verticalize() and
      historicalize() functions. Failure to format the input data properly and
      designate needed variables appropriately may result in nonsensical output.",
      call. = FALSE)
  }
  
  if (test.group) {
    if (is.null(stageframe)) {
      stop("Cannot test groups without inclusion of appropriate stageframe.",
        call. = FALSE)
    } else if (!is.element("stageframe", class(stageframe))) {
      stop("Cannot test groups without inclusion of appropriate stageframe.",
        call. = FALSE)
    } else {
      all_groups <- unique(stageframe$group)
      
      if (length(all_groups) > 1) {
        data$group2 <- apply(as.matrix(data$stage2), 1, function(X){
          found_group <- which(stageframe$stage == X)
          if (length(found_group) != 1) {
            if (!is.element(tolower(X), c("notalive", "dead", "almostborn"))) {
              warning("Some group calls appear to lack a positive ID. Please check
                group identifications.", call. = FALSE)
            }
            group_call <- 0
          } else {
            group_call <- stageframe$group[found_group]
          }
          return(group_call)
        })
        data$group3 <- apply(as.matrix(data$stage3), 1, function(X){
          found_group <- which(stageframe$stage == X)
          if (length(found_group) != 1) {
            if (!is.element(tolower(X), c("notalive", "dead", "almostborn"))) {
              warning("Some group calls appear to lack a positive ID. Please check
                group identifications.", call. = FALSE)
            }
            group_call <- 0
          } else {
            group_call <- stageframe$group[found_group]
          }
          return(group_call)
        })
        data$group1 <- apply(as.matrix(data$stage1), 1, function(X){
          found_group <- which(stageframe$stage == X)
          if (length(found_group) != 1) {
            if (!is.element(tolower(X), c("notalive", "dead", "almostborn"))) {
              warning("Some group calls appear to lack a positive ID. Please check
                group identifications.", call. = FALSE)
            }
            group_call <- 0
          } else {
            group_call <- stageframe$group[found_group]
          }
          return(group_call)
        })
      } else {
        test.group <- FALSE
        warning("Only one stage group found, so will not test group.", call. = FALSE)
      }
    }
  }
  
  if (!requireNamespace("MuMIn", quietly = TRUE)) {
    stop("Package MuMIn is required. Please install it.",
      call. = FALSE)
  }
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package stringr needed for this function to work. Please install it.",
      call. = FALSE)
  }
  
  if (size.zero & size.trunc) {
    stop("Size distribution cannot be both zero-inflated and zero-truncated.
      Please set size.zero, size.trunc, or both to FALSE.", call. = FALSE)
  }
  if (sizeb.zero & sizeb.trunc) {
    stop("Sizeb distribution cannot be both zero-inflated and zero-truncated.
      Please set size.zero, size.trunc, or both to FALSE.", call. = FALSE)
  }
  if (sizec.zero & sizec.trunc) {
    stop("Sizec distribution cannot be both zero-inflated and zero-truncated.
      Please set size.zero, size.trunc, or both to FALSE.", call. = FALSE)
  }
  if (fec.zero & fec.trunc) {
    stop("Fecundity distribution cannot be both zero-inflated and zero-truncated.
      Please set fec.zero, fec.trunc, or both to FALSE.", call. = FALSE)
  }
  if (jsize.zero & jsize.trunc) {
    stop("Juvenile primary size distribution cannot be both zero-inflated and zero-truncated.
      Please set jsize.zero, jsize.trunc, or both to FALSE.", call. = FALSE)
  }
  if (jsizeb.zero & jsizeb.trunc) {
    stop("Juvenile secondary size distribution cannot be both zero-inflated and zero-truncated.
      Please set jsize.zero, jsize.trunc, or both to FALSE.", call. = FALSE)
  }
  if (jsizec.zero & jsizec.trunc) {
    stop("Juvenile tertiary size distribution cannot be both zero-inflated and zero-truncated.
      Please set jsize.zero, jsize.trunc, or both to FALSE.", call. = FALSE)
  }
  
  approach <- tolower(approach)
  sizedist <- tolower(sizedist)
  sizebdist <- tolower(sizebdist)
  sizecdist <- tolower(sizecdist)
  fecdist <- tolower(fecdist)
  
  if (is.element(approach, c("lme4", "glmmtmb", "mixed"))) {approach <- "mixed"}
  
  if (approach == "mixed" & !requireNamespace("lme4", quietly = TRUE)) {
    stop("Package lme4 is required for mixed models. Please install it.", call. = FALSE)
  } else if (approach == "mixed" & !requireNamespace("glmmTMB", quietly = TRUE)) {
    if (sizedist == "negbin" | sizebdist == "negbin" | sizecdist == "negbin" | fecdist == "negbin") {
      stop("Package glmmTMB needed to develop mixed size or fecundity models
        with a negative binomial distribution.", call. = FALSE)
    }
  } else if (sizedist != "gaussian" | any(sizebdist != "gaussian", na.rm = TRUE) | 
      any(sizecdist != "gaussian", na.rm = TRUE)) {
    if (!requireNamespace("glmmTMB", quietly = TRUE)) {
      if (size.trunc | sizeb.trunc | sizec.trunc | size.zero | sizeb.zero | sizec.zero) {
        stop("Package glmmTMB needed to develop mixed size models with
          zero-truncated or zero-inflated distributions.", call. = FALSE)
      }
    }
  } else if (fecdist != "gaussian") {
    if (fec.trunc & !requireNamespace("glmmTMB", quietly = TRUE)) {
      stop("Package glmmTMB needed to develop mixed fecundity models with
        zero-truncated distribution.", call. = FALSE)}
    if (fec.zero & !requireNamespace("glmmTMB", quietly = TRUE)) {
      stop("Package glmmTMB needed to develop mixed fecundity models with
        zero-inflated distribution.", call. = FALSE)}
  }
  
  if (approach == "glm") {
    if (sizedist != "gaussian" | any(sizebdist != "gaussian", na.rm = TRUE) |
        any(sizecdist != "gaussian", na.rm = TRUE)) {
      if (size.trunc | sizeb.trunc | sizec.trunc) {
        if (!requireNamespace("VGAM", quietly = TRUE)) {
        stop("Package VGAM needed to develop non-Gaussian size GLMs with 
          zero-truncated distribution.", call. = FALSE)
        }
      } else if (size.zero | sizeb.zero | sizec.zero) {
        if (!requireNamespace("pscl", quietly = TRUE)) {
        stop("Package pscl needed to develop non-Gaussian size GLMs with
          zero-inflated distribution.", call. = FALSE)
        }
      }
    }
    if (fecdist != "gaussian") {
      if (fec.trunc & !requireNamespace("VGAM", quietly = TRUE)) {
        stop("Package VGAM needed to develop non-Gaussian fecundity GLMs with
          zero-truncated distribution.", call. = FALSE)}
      if (fec.zero & !requireNamespace("pscl", quietly = TRUE)) {
        stop("Package pscl needed to develop non-Gaussian fecundity GLMs with
          zero-inflated distribution.", call. = FALSE)}
    }
  }
  
  distoptions <- c("gaussian", "poisson", "negbin", "gamma")
  packoptions <- c("mixed", "glm") #The mixed option now handles all mixed models
  
  if (!is.element(approach, packoptions)) {
    stop("Please enter a valid approach, currently either 'mixed' or 'glm'.", 
      call. = FALSE)}
  if (!is.element(sizedist, distoptions)) {
    stop("Please enter a valid assumed size distribution, currently limited to
      gaussian, poisson, negbin, and gamma.", call. = FALSE)}
  if (!is.element(sizebdist, c(NA, distoptions))) {
    stop("Please enter a valid assumed sizeb distribution, currently limited to
      gaussian, poisson, negbin, and gamma.", call. = FALSE)}
  if (!is.element(sizecdist, c(NA, distoptions))) {
    stop("Please enter a valid assumed sizec distribution, currently limited to
      gaussian, poisson, negbin, and gamma.", call. = FALSE)}
  if (!is.element(fecdist, distoptions)) {
    stop("Please enter a valid assumed fecundity distribution, currently limited
      to gaussian, poisson, negbin, and gamma.", call. = FALSE)}
  
  if (length(censor) > 3) {
    stop("Censor variables should be included either as 1 variable per row in
      the historical data file (1 variable in the dataset), or as 1 variable for
      each of occasions t+1, t, and, if historical, t-1 (2 or 3 variables in the
      dataset). No more than 3 variables are allowed. If more than one are
      supplied, then they are assumed to be in order of occasion t+1,
      occasion t, and occasion t-1, respectively.", call. = FALSE)
  }
  if (length(indiv) > 1) {
    stop("Only one individual identification variable is allowed.",
      call. = FALSE)
  }
  if (length(year) > 1) {
    stop("Only one observation occasion variable is allowed, and it must refer
      to occasion t.", call. = FALSE)
  }
  if (length(patch) > 1) {
    stop("Only one patch variable is allowed.", call. = FALSE)
  }
  if (length(density) > 1) {
    stop("Only one density variable is allowed.", call. = FALSE)
  }
  
  if (is.element("surv", vitalrates)) {
    if (length(surv) > 3 | length(surv) == 1) {
      stop("This function requires 2 (if ahistorical) or 3 (if historical)
        survival variables from the dataset as input parameters.",
        call. = FALSE)}
    if (all(is.numeric(surv))) {
      if (!any(surv < 1) & !any(surv > length(names(data)))) {
        surv <- names(data)[surv]
      } else {
        stop("Survival variables do not match those in the dataset.",
          call. = FALSE)
      }
    }
    if (any(!is.element(surv, names(data)))) {
      stop("Survival variables need to match those used in the dataset.", call. = FALSE)
    }
  }
  if (is.element("obs", vitalrates)) {
    if (length(obs) > 3 | length(obs) == 1) {
      stop("This function requires 2 (if ahistorical) or 3 (if historical)
        observation variables from the dataset as input parameters.",
        call. = FALSE)}
    if (all(is.numeric(obs))) {
      if (!any(obs < 1) & !any(obs > length(names(data)))) {
        obs <- names(data)[obs]
      } else {
        stop("Observation variables do not match those in the dataset.",
          call. = FALSE)
      }
    }
    if (any(!is.element(obs, names(data)))) {
      stop("Observation status variables need to match those used in the dataset.",
        call. = FALSE)
    }
  }
  if (is.element("size", vitalrates)) {
    if (length(size) > 3 | length(size) == 1) {
      stop("This function requires 2 (if ahistorical) or 3 (if historical) size
        variables from the dataset as input parameters.", call. = FALSE)}
    if (all(is.numeric(size))) {
      if (!any(size < 1) & !any(size > length(names(data)))) {
        size <- names(data)[size]
      } else {
        stop("Size variables do not match those in the dataset.",
          call. = FALSE)
      }
    }
    if (any(!is.element(size, names(data)))) {
      stop("Size variables need to match those used in the dataset.", call. = FALSE)
    }
  }
  if (all(!is.na(sizeb))) {
    if (is.na(sizebdist)) {
      stop("Need valid choice of distribution for secondary size.", call. = FALSE)
    }
    if (length(sizeb) > 3 | length(sizeb) == 1) {
      stop("This function requires 2 (if ahistorical) or 3 (if historical)
        secondary size variables from the dataset as input parameters.",
        call. = FALSE)}
    if (all(is.numeric(sizeb))) {
      if (!any(sizeb < 1) & !any(sizeb > length(names(data)))) {
        sizeb <- names(data)[sizeb]
      } else {
        stop("Secondary variables do not match those in the dataset.",
          call. = FALSE)
      }
    }
    if (any(!is.element(sizeb, names(data)))) {
      stop("Secondary size variables need to match those used in the dataset.",
        call. = FALSE)
    }
    sizeb_used <- 1
  }
  if (all(!is.na(sizec))) {
    if (is.na(sizecdist)) {
      stop("Need valid choice of distribution for tertiary size.", call. = FALSE)
    }
    if (length(sizec) > 3 | length(sizec) == 1) {
      stop("This function requires 2 (if ahistorical) or 3 (if historical)
        tertiary size variables from the dataset as input parameters.",
        call. = FALSE)}
    if (all(is.numeric(sizec))) {
      if (!any(sizec < 1) & !any(sizec > length(names(data)))) {
        sizec <- names(data)[sizec]
      } else {
        stop("Tertiary size variables do not match those in the dataset.",
          call. = FALSE)
      }
    }
    if (any(!is.element(sizec, names(data)))) {
      stop("Secondary size variables need to match those used in the dataset.",
        call. = FALSE)
    }
    sizec_used <- 1
  }
  if (is.element("repst", vitalrates)) {
    if (length(repst) > 3 | length(repst) == 1) {
      stop("This function requires 2 (if ahistorical) or 3 (if historical)
        reproductive status variables from the dataset as input parameters.",
        call. = FALSE)}
    if (all(is.numeric(repst))) {
      if (!any(repst < 1) & !any(repst > length(names(data)))) {
        repst <- names(data)[repst]
      } else {
        stop("Reproductive status variables do not match those in the dataset.",
          call. = FALSE)
      }
    }
    if (any(!is.element(repst, names(data)))) {
      stop("Reproductive status variables need to match those used in the dataset.",
        call. = FALSE)
    }
  }
  if (is.element("fec", vitalrates)) {
    if (length(fec) > 3 | length(fec) == 1) {
      stop("This function requires 2 (if ahistorical) or 3 (if historical)
        fecundity variables from the dataset as input parameters.",
        call. = FALSE)}
    if (all(is.numeric(fec))) {
      if (!any(fec < 1) & !any(fec > length(names(data)))) {
        fec <- names(data)[fec]
      } else {
        stop("Fecundity variables do not match those in the dataset.",
          call. = FALSE)
      }
    }
    if (any(!is.element(fec, names(data)))) {
      stop("Fecundity variables need to match those used in the dataset.",
        call. = FALSE)
    }
  }
  
  if (fectime != 2 & fectime != 3) {
    stop("The fectime option must equal either 2 or 3, depending on whether the
      fecundity response term is for occasion t or occasion t+1, respectively.
      The default is 2, corresponding to occasion t.", call. = FALSE)
  }
  
  if (!is.na(age)) {
    if (length(which(names(data) == age)) == 0) {
      stop("Variable age must either equal the exact name of the variable
           denoting age in the dataset, or be set to NA.", call. = FALSE)
    } else {
      agecol <- which(names(data) == age)
    }
  } else {age <- "none"}
  
  if (any(!is.na(indcova))) {
    if (length(indcova) > 3) {
      warning("Vector indcova holds the exact names of an individual covariate across
        occasions t+1, t, and t-1. Only the first three elements will be used.",
        call. = FALSE)
      indcova <- indcova[1:3]
      indcova_used <- TRUE
    } else if (length(indcova) == 1) {
      warning("Vector indcova requires the names of an individual covariate across
        occasions t+1, t, and, if historical, t-1. Only 1 variable name was supplied,
        so this individual covariate will not be used.", call. = FALSE)
      indcova <- c("none", "none", "none")
    } else {
      if (length(which(names(data) == indcova[2])) == 0 && indcova[2] != "none") {
        stop("Vector indcova must either equal either the exact names of an
          individual covariate across occasions t+1, t, and t-1, or be set to NA.",
          call. = FALSE)
      } else {
        indcova2col <- which(names(data) == indcova[2])
        indcova_used <- TRUE
      }
      
      if (length(indcova) == 3) {
        if (length(which(names(data) == indcova[3])) == 0 && indcova[3] != "none") {
          stop("Vector indcova must either equal either the exact names of an
            individual covariate across occasions t+1, t, and t-1, or be set to NA.",
            call. = FALSE)
        } else {
          indcova1col <- which(names(data) == indcova[3])
          indcova_used <- TRUE
        }
      } else {
        indcova[3] <- "none"
      }
    }
  } else {indcova <- c("none", "none", "none")}
  
  if (any(!is.na(indcovb))) {
    if (length(indcovb) > 3) {
      warning("Vector indcovb holds the exact names of an individual covariate
        across occasions t+1, t, and t-1. Only the first three elements will be used.",
        call. = FALSE)
      indcovb <- indcovb[1:3]
      indcovb_used <- TRUE
    } else if (length(indcovb) == 1) {
      warning("Vector indcovb requires the names of an individual covariate across
        occasions t+1, t, and, if historical, t-1. Only 1 variable name was supplied,
        so this individual covariate will not be used.", call. = FALSE)
      indcovb <- c("none", "none", "none")
    } else {
      if (length(which(names(data) == indcovb[2])) == 0 && indcovb[2] != "none") {
        stop("Vector indcovb must either equal either the exact names of an
          individual covariate across occasions t+1, t, and t-1, or be set to NA.",
          call. = FALSE)
      } else {
        indcovb2col <- which(names(data) == indcovb[2])
        indcovb_used <- TRUE
      }
      
      if (length(indcovb) == 3) {
        if (length(which(names(data) == indcovb[3])) == 0 && indcovb[3] != "none") {
          stop("Vector indcovb must either equal either the exact names of an
            individual covariate across occasions t+1, t, and t-1, or be set to NA.",
            call. = FALSE)
        } else {
          indcovb1col <- which(names(data) == indcovb[3])
          indcovb_used <- TRUE
        }
      } else {
        indcovb[3] <- "none"
      }
    }
  } else {indcovb <- c("none", "none", "none")}
  
  if (any(!is.na(indcovc))) {
    if (length(indcovc) > 3) {
      warning("Vector indcovc holds the exact names of an individual covariate
        across occasions t+1, t, and t-1. Only the first three elements will be used.",
        call. = FALSE)
      indcovc <- indcovc[1:3]
      indcovc_used <- TRUE
    } else if (length(indcovc) == 1) {
      warning("Vector indcovc requires the names of an individual covariate across
        occasions t+1, t, and, if historical, t-1. Only 1 variable name was supplied,
        so this individual covariate will not be used.", call. = FALSE)
      indcovc <- c("none", "none", "none")
    } else {
      if (length(which(names(data) == indcovc[2])) == 0 && indcovc[2] != "none") {
        stop("Vector indcovc must either equal either the exact names of an
          individual covariate across occasions t+1, t, and t-1, or be set to NA.",
          call. = FALSE)
      } else {
        indcovc2col <- which(names(data) == indcovc[2])
        indcovc_used <- TRUE
      }
      
      if (length(indcovc) == 3) {
        if (length(which(names(data) == indcovc[3])) == 0 && indcovc[3] != "none") {
          stop("Vector indcovc must either equal either the exact names of an
            individual covariate across occasions t+1, t, and t-1, or be set to NA.",
            call. = FALSE)
        } else {
          indcovc1col <- which(names(data) == indcovc[3])
          indcovc_used <- TRUE
        }
      } else {
        indcovc[3] <- "none"
      }
    }
  } else {indcovc <- c("none", "none", "none")}
  
  if (!is.na(indiv)) {
    if (!is.numeric(indiv) & length(which(names(data) == indiv)) == 0) {
      stop("Variable indiv must either equal the exact name of the variable denoting
        individual identity in the dataset, or be set to NA.", call. = FALSE)
    } else if (is.character(indiv)){
      indivcol <- which(names(data) == indiv)
      
      if (any(is.na(data[,indivcol])) & approach == "mixed") {
        warning("NAs in individual ID variable may cause unexpected behavior in mixed
          model building. Please rename all individuals with unique names, avoiding NAs.",
          call. = FALSE)
      }
    } else if (is.numeric(indiv)) {
      if (!any(indiv < 1) & !any(indiv > length(names(data)))) {
        indivcol <- names(data)[indiv]
      }
    } else {
      stop("Unable to interpret indiv variable.", call. = FALSE)
    }
  } else {indiv <- "none"}
  
  if (!is.na(patch)) {
    if (!is.numeric(patch) & length(which(names(data) == patch)) == 0) {
      stop("Variable patch must either equal the exact name of the variable denoting
        patch identity in the dataset, or be set to NA.", call. = FALSE)
    } else if (is.character(patch)) {
      patchcol <- which(names(data) == patch)
    } else if (is.numeric(patch)) {
      if (!any(patch < 1) & !any(patch > length(names(data)))) {
        patchcol <- names(data)[patch]
      }
    } else {
      stop("Unable to interpret patch variable.", call. = FALSE)
    }
  } else {patch <- "none"}
  
  if (!is.na(year)) {
    if (!is.numeric(year) & length(which(names(data) == year)) == 0) {
      stop("Variable year must either equal the exact name of the variable denoting
        occasion t in the dataset, or be set to NA.", call. = FALSE)
    } else if (is.character(year)) {
      yearcol <- which(names(data) == year)
    } else if (is.numeric(year)) {
      if (!any(year < 1) & !any(year > length(names(data)))) {
        yearcol <- names(data)[year]
      }
    } else {
      stop("Unable to interpret year variable.", call. = FALSE)
    }
  } else {year <- "none"}
  
  if (!is.na(density)) {
    if (!is.numeric(density) & length(which(names(data) == density)) == 0) {
      stop("Variable density must either equal the exact name of the variable
        denoting spatial density in occasion t in the dataset, or be set to NA.",
        call. = FALSE)
    } else if (is.character(density)) {
      if (is.element(density, names(data))) {
        density_used <- TRUE
      } else {
        stop("Unable to interpret density variable.", call. = FALSE)
      }
    } else if (is.numeric(density)) {
      if (!any(density < 1) & !any(density > length(names(data)))) {
        density <- names(data)[density]
        density_used <- TRUE
      }
    } else {
      stop("Unable to interpret density variable.", call. = FALSE)
    }
  } else {density <- "none"}
  
  # Here we test the dataset for appropriate stage names
  if (!is.na(juvestimate)) {
    if (!any(is.element(stage, names(data)))) {
      stop("Names of stage variables do not match the dataset.", call. = FALSE)}
    
    stage3col <- which(names(data) == stage[1])
    stage2col <- which(names(data) == stage[2])
    if (length(stage) == 3) {stage1col <- which(names(data) == stage[3])} else {stage1col <- 0}
    
    if (!is.element(juvestimate, unique(data[,stage2col]))) {
      stop("The stage declared as juvenile via juvestimate is not recognized within
        the stage2 variable in the dataset.", call. = FALSE)
    }
  }
  
  if (!is.logical(juvsize)) {
    stop("Option juvsize must be set to either TRUE or FALSE.", call. = FALSE)
  }
  
  #Now we check whether the best-fit criterion is appropriate
  if (!is.element(tolower(bestfit), c("aicc", "aicc&k"))) {
    stop("Bestfit must equal either 'AICc' or 'AICc&k', with the latter as the default.",
      call. = FALSE)
  }
  #This variable will be used once dredging is done to determine best-fit models
  used.criterion <- gsub("&k", "", bestfit)
  
  if (approach != "mixed") {
    if (random.indcova | random.indcovb | random.indcovc) {
      warning("Random covariates can only be included in mixed models. Setting
        random.indcova, random.indcovb, and random.indcovc to FALSE.",
        call. = FALSE)
      random.indcova <- FALSE
      random.indcovb <- FALSE
      random.indcovc <- FALSE
    }
  }
  
  #The next section creates the text lines needed for the main model calls, based on function input
  formulae <- .stovokor(surv, obs, size, sizeb, sizec, repst, fec, vitalrates,
    historical, suite, approach, is.na(juvestimate), age, indcova, indcovb,
    indcovc, indiv, patch, year, patch.as.random, year.as.random,
    random.indcova, random.indcovb, random.indcovc, fectime, juvsize,
    sizeb_used, sizec_used, test.group, density, density_used, indcova_used,
    indcovb_used, indcovc_used)
  
  if (suite == "full") {
    alt.formulae <- .stovokor(surv, obs, size, sizeb, sizec, repst, fec, vitalrates,
      historical, suite = "main", approach, is.na(juvestimate), age, indcova,
      indcovb, indcovc, indiv, patch, year, patch.as.random, year.as.random,
      random.indcova, random.indcovb, random.indcovc, fectime, juvsize,
      sizeb_used, sizec_used, test.group, density, density_used, indcova_used,
      indcovb_used, indcovc_used)
  } else {
    alt.formulae <- list(full.surv.model = NA, full.obs.model = NA, full.size.model = NA,
      full.sizeb.model = NA, full.sizec.model = NA, full.repst.model = NA, full.fec.model = NA,
      juv.surv.model = NA, juv.obs.model = NA, juv.size.model = NA, just.sizeb.model = NA,
      juve.sizec.model = NA, juve.repst.model = NA)
  }
  
  #Now we need to create the input datasets
  if (!all(is.na(censor))) {
    if (length(censor) == 1) {
      data$censor2 <- data[, which(names(data) == censor[1])]
    } else {
      data$censor3 <- data[, which(names(data) == censor[1])]
      data$censor2 <- data[, which(names(data) == censor[2])]
      if (length(censor) > 2) {
        data$censor1 <- data[, which(names(data) == censor[3])]
      }
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
    juvsurv.ind <- length(unique(juvsurv.data[, which(names(juvsurv.data) == indiv)]))
    juvsurv.trans <- dim(juvsurv.data)[1]
    
    if (suite == "full" | suite == "main" | suite == "size") {
      if (any(is.na(juvsurv.data[, which(names(juvsurv.data) == size[2])]))) {
        warning("NAs in size variables may cause model selection to fail.",
          call. = FALSE)
      }
      
      if (historical == TRUE) {
        if (any(is.na(juvsurv.data[, which(names(juvsurv.data) == size[3])]))) {
          warning("NAs in size variables may cause model selection to fail.",
            call. = FALSE)
        }
      }
    }
    
    if (suite == "full" | suite == "main" | suite == "rep") {
      if (any(is.na(juvsurv.data[, which(names(juvsurv.data) == repst[2])]))) {
        warning("NAs in reproductive status variables may cause model selection to fail.",
          call. = FALSE)
      }
      
      if (historical == TRUE) {
        if (any(is.na(juvsurv.data[, which(names(juvsurv.data) == repst[3])]))) {
          warning("NAs in reproductive status variables may cause model selection to fail.",
            call. = FALSE)
        }
      }
    }
    
    if (is.element(0, juvsurv.data$matstatus3)) {
      warning("Function modelsearch() assumes that all juveniles either die or transition
        to maturity within 1 year. Some individuals in this dataset appear to live
        longer as juveniles than assumptions allow.", call. = FALSE)
    }
    
    juvobs.data <- subset(juvsurv.data, juvsurv.data[, which(names(juvsurv.data) == surv[1])] == 1)
    juvobs.ind <- length(unique(juvobs.data[, which(names(juvobs.data) == indiv)]))
    juvobs.trans <- dim(juvobs.data)[1]
    
    if (formulae$full.obs.model != 1) {
      juvsize.data <- subset(juvobs.data, juvobs.data[, which(names(juvobs.data) == obs[1])] == 1)
      juvsize.data <- juvsize.data[which(!is.na(juvsize.data[, which(names(juvsize.data) == size[1])])),]
      
      if (sizeb_used) {
        juvsizeb.data <- subset(juvobs.data, juvobs.data[, which(names(juvobs.data) == obs[1])] == 1)
        juvsizeb.data <- juvsizeb.data[which(!is.na(juvsizeb.data[, which(names(juvsizeb.data) == sizeb[1])])),]
      }
      
      if (sizec_used) {
        juvsizec.data <- subset(juvobs.data, juvobs.data[, which(names(juvobs.data) == obs[1])] == 1)
        juvsizec.data <- juvsizec.data[which(!is.na(juvsizec.data[, which(names(juvsizec.data) == sizec[1])])),]
      }
    } else {
      juvsize.data <- juvobs.data
      juvsize.data <- juvsize.data[which(!is.na(juvsize.data[, which(names(juvsize.data) == size[1])])),]
      
      if (sizeb_used) {
        juvsizeb.data <- juvobs.data
        juvsizeb.data <- juvsizeb.data[which(!is.na(juvsizeb.data[, which(names(juvsizeb.data) == sizeb[1])])),]
      }
      
      if (sizec_used) {
        juvsizec.data <- juvobs.data
        juvsizec.data <- juvsizec.data[which(!is.na(juvsizec.data[, which(names(juvsizec.data) == sizec[1])])),]
      }
    }
    juvsize.ind <- length(unique(juvsize.data[, which(names(juvsize.data) == indiv)]))
    juvsize.trans <- dim(juvsize.data)[1]
    
    if (sizeb_used) {
      juvsizeb.ind <- length(unique(juvsizeb.data[, which(names(juvsizeb.data) == indiv)]))
      juvsizeb.trans <- dim(juvsizeb.data)[1]
    }
    if (sizec_used) {
      juvsizec.ind <- length(unique(juvsizec.data[, which(names(juvsizec.data) == indiv)]))
      juvsizec.trans <- dim(juvsizec.data)[1]
    }
    
    juvrepst.data <- juvsize.data
    juvrepst.ind <- length(unique(juvrepst.data[, which(names(juvrepst.data) == indiv)]))
    juvrepst.trans <- dim(juvrepst.data)[1]
    
    data <- data[adultindivs,] #This line resets the main dataset to adults only
  }
  
  surv.data <- subset(data, data[,which(names(data) == surv[2])] == 1)
  surv.ind <- length(unique(surv.data[, which(names(surv.data) == indiv)]))
  surv.trans <- dim(surv.data)[1]
  
  if(any(!suppressWarnings(!is.na(as.numeric(as.character(surv.data[, which(names(surv.data) == size[1])])))))) {
    warning("Modelsearch(), flefko3(), flefko2(), and aflefko2() are made to work
      with numeric size variables. Use of categorical variables may result in
      errors and unexpected behavior.", call. = FALSE)
  }
  if (suite == "full" | suite == "main" | suite == "size") {
    if (any(is.na(surv.data[, which(names(surv.data) == size[2])]))) {
      warning("NAs in size variables may cause model selection to fail.",
        call. = FALSE)
    }
    if (sizeb_used) {
      if (any(is.na(surv.data[, which(names(surv.data) == sizeb[2])]))) {
        warning("NAs in size variables may cause model selection to fail.",
          call. = FALSE)
      }
    }
    if (sizec_used) {
      if (any(is.na(surv.data[, which(names(surv.data) == sizec[2])]))) {
        warning("NAs in size variables may cause model selection to fail.",
          call. = FALSE)
      }
    }
    
    if (historical == TRUE) {
      if (any(is.na(surv.data[, which(names(surv.data) == size[3])]))) {
        warning("NAs in size variables may cause model selection to fail.",
          call. = FALSE)
      }
      if (sizeb_used) {
        if (any(is.na(surv.data[, which(names(surv.data) == sizeb[3])]))) {
          warning("NAs in size variables may cause model selection to fail.",
            call. = FALSE)
        }
      }
      if (sizec_used) {
        if (any(is.na(surv.data[, which(names(surv.data) == sizec[3])]))) {
          warning("NAs in size variables may cause model selection to fail.",
            call. = FALSE)
        }
      }
    }
  }
  if (suite == "full" | suite == "main" | suite == "rep") {
    if (any(is.na(surv.data[, which(names(surv.data) == repst[2])]))) {
      warning("NAs in reproductive status variables may cause model selection to fail.",
        call. = FALSE)
    }
    
    if (historical == TRUE) {
      if (any(is.na(surv.data[, which(names(surv.data) == repst[3])]))) {
        warning("NAs in reproductive status variables may cause model selection to fail.",
          call. = FALSE)
      }
    }
  }
  
  surv.uns <- unique(surv.data[,which(names(surv.data) == surv[1])])
  if (length(surv.uns) == 1) {
    warning("Survival to occasion t+1 appears to be constant, and so will be set to a constant.",
      call. = FALSE)
    formulae$full.surv.model <- surv.uns[1]
    alt.formulae$full.surv.model <- surv.uns[1]
  }
  
  obs.data <- subset(surv.data, surv.data[, which(names(surv.data) == surv[1])] == 1)
  obs.ind <- length(unique(obs.data[, which(names(obs.data) == indiv)]))
  obs.trans <- dim(obs.data)[1]
  
  if (formulae$full.obs.model != 1) {
    obs.uns <- unique(obs.data[,which(names(obs.data) == obs[1])])
    if (length(obs.uns) == 1) {
      warning("Observation in occasion t+1 appears to be constant, and so will be
        set to a constant.", call. = FALSE)
    formulae$full.obs.model <- obs.uns[1]
    alt.formulae$full.obs.model <- obs.uns[1]
    }
  }

  if (formulae$full.obs.model != 1) {
    size.data <- subset(obs.data, obs.data[, which(names(obs.data) == obs[1])] == 1)
    size.data <- size.data[which(!is.na(size.data[, which(names(size.data) == size[1])])),]
    
    if (sizeb_used) {
      sizeb.data <- subset(obs.data, obs.data[, which(names(obs.data) == obs[1])] == 1)
      sizeb.data <- sizeb.data[which(!is.na(sizeb.data[, which(names(sizeb.data) == sizeb[1])])),]
    }
    if (sizec_used) {
      sizec.data <- subset(obs.data, obs.data[, which(names(obs.data) == obs[1])] == 1)
      sizec.data <- sizec.data[which(!is.na(sizec.data[, which(names(sizec.data) == sizec[1])])),]
    }
  } else {
    size.data <- obs.data
    size.data <- size.data[which(!is.na(size.data[, which(names(size.data) == size[1])])),]
    
    if (sizeb_used) {
      sizeb.data <- obs.data
      sizeb.data <- sizeb.data[which(!is.na(sizeb.data[, which(names(sizeb.data) == sizeb[1])])),]
    }
    if (sizec_used) {
      sizec.data <- obs.data
      sizec.data <- sizec.data[which(!is.na(sizec.data[, which(names(sizec.data) == sizec[1])])),]
    }
  }
  size.ind <- length(unique(size.data[, which(names(size.data) == indiv)]))
  size.trans <- dim(size.data)[1]
  
  if (sizeb_used) {
    sizeb.ind <- length(unique(sizeb.data[, which(names(sizeb.data) == indiv)]))
    sizeb.trans <- dim(sizeb.data)[1]
  }
  if (sizec_used) {
    sizec.ind <- length(unique(sizec.data[, which(names(sizec.data) == indiv)]))
    sizec.trans <- dim(sizec.data)[1]
  }
  
  if (formulae$full.size.model != 1) {
    size.uns <- unique(size.data[,which(names(size.data) == size[1])])
    if (length(size.uns) == 1) {
      warning("Size in occasion t+1 appears to be constant, and so will be set to a constant.",
        call. = FALSE)
    formulae$full.size.model <- size.uns[1]
    alt.formulae$full.size.model <- size.uns[1]
    }
  }
  
  if (formulae$full.sizeb.model != 1 & sizeb_used) {
    sizeb.uns <- unique(sizeb.data[,which(names(sizeb.data) == sizeb[1])])
    if (length(sizeb.uns) == 1) {
      warning("Sizeb in occasion t+1 appears to be constant, and so will be set to a constant.",
        call. = FALSE)
    formulae$full.sizeb.model <- sizeb.uns[1]
    alt.formulae$full.sizeb.model <- sizeb.uns[1]
    }
  }
  
  if (formulae$full.sizec.model != 1 & sizec_used) {
    sizec.uns <- unique(sizec.data[,which(names(sizec.data) == sizec[1])])
    if (length(sizec.uns) == 1) {
      warning("Sizec in occasion t+1 appears to be constant, and so will be set to a constant.",
        call. = FALSE)
    formulae$full.sizec.model <- sizec.uns[1]
    alt.formulae$full.sizec.model <- sizec.uns[1]
    }
  }
  
  repst.data <- size.data
  repst.ind <- length(unique(repst.data[, which(names(repst.data) == indiv)]))
  repst.trans <- dim(repst.data)[1]
  
  if (formulae$full.repst.model != 1) {
    repst.uns <- unique(repst.data[,which(names(repst.data) == repst[1])])
    if (length(repst.uns) == 1) {
      warning("Reproductive status in occasion t+1 appears to be constant, and so will
        be set to a constant.", call. = FALSE)
      formulae$full.repst.model <- repst.uns[1]
      alt.formulae$full.repst.model <- repst.uns[1]
    }
  }
  
  if (formulae$full.repst.model != 1) {
    fec.data <- subset(surv.data, surv.data[, which(names(repst.data) == repst[2])] == 1)
    if (fectime == 2) {
      fec.data <- fec.data[which(!is.na(fec.data[, which(names(fec.data) == fec[2])])),]
      fec.data <- fec.data[which(fec.data[, which(names(fec.data) == repst[2])] == 1),]
    } else if (fectime == 3) {
      fec.data <- fec.data[which(!is.na(fec.data[, which(names(fec.data) == fec[1])])),]
      fec.data <- fec.data[which(fec.data[, which(names(fec.data) == repst[1])] == 1),]
    }
  } else {
    fec.data <- surv.data
    if (fectime == 2) {
      fec.data <- fec.data[which(!is.na(fec.data[, which(names(fec.data) == fec[2])])),]
    } else if (fectime == 3) {
      fec.data <- fec.data[which(!is.na(fec.data[, which(names(fec.data) == fec[1])])),]
    }
  }
  fec.ind <- length(unique(fec.data[, which(names(fec.data) == indiv)]))
  fec.trans <- dim(fec.data)[1]
  
  if (formulae$full.fec.model != 1) {
    if (fectime == 2) {
      fec.uns <- unique(fec.data[,which(names(fec.data) == fec[2])])
    } else {
      fec.uns <- unique(fec.data[,which(names(fec.data) == fec[1])])
    }
    
    if (length(fec.uns) == 1) {
      warning("Fecundity appears to be constant, and so will be set to a constant.",
        call. = FALSE)
      formulae$full.fec.model <- fec.uns[1]
      alt.formulae$full.fec.model <- fec.uns[1]
    }
  }
  
  #Now we check for exceptions to size and fecundity in the dataset
  if (!is.numeric(formulae$full.size.model)) {
    if (sizedist == "poisson") {
      if (any(size.data[, which(names(size.data) == size[1])] != 
          round(size.data[, which(names(size.data) == size[1])]))) {
        stop("Size variables must be composed only of integers for the Poisson
          distribution to be used.", call. = FALSE)
      }
      
      if (!is.na(juvestimate)) {
        if (any(juvsize.data[, which(names(juvsize.data) == size[1])] != 
            round(juvsize.data[, which(names(juvsize.data) == size[1])]))) {
          stop("Size variables must be composed only of integers for the Poisson
            distribution to be used.", call. = FALSE)
        }
      }
    } else if (sizedist == "negbin") {
      if (any(size.data[, which(names(size.data) == size[1])] != 
          round(size.data[, which(names(size.data) == size[1])]))) {
        stop("Size variables must be composed only of integers for the negative
          binomial distribution to be used.", call. = FALSE)
      }
      
      if (!is.na(juvestimate)) {
        if (any(juvsize.data[, which(names(juvsize.data) == size[1])] != 
            round(juvsize.data[, which(names(juvsize.data) == size[1])]))) {
          stop("Size variables must be composed only of integers for the negative
            binomial distribution to be used.", call. = FALSE)
        }
      }
    }
  } else if (sizedist == "gamma") {
    if (any(size.data[, which(names(size.data) == size[1])] < 0, na.rm = TRUE)) {
      stop("Size variables must be non-negative for the gamma distribution to be used.",
        call. = FALSE)
    }
    
    if (!is.na(juvestimate)) {
      if (any(juvsize.data[, which(names(juvsize.data) == size[1])]  < 0, na.rm = TRUE)) {
        stop("Size variables must be non-negative for the gamma distribution to be used.",
          call. = FALSE)
      }
    }
  }
    
  if (!is.numeric(formulae$full.sizeb.model) & sizeb_used) {
    if (sizebdist == "poisson") {
      if (any(sizeb.data[, which(names(sizeb.data) == sizeb[1])] != 
          round(sizeb.data[, which(names(sizeb.data) == sizeb[1])]))) {
        stop("Secondary size variables must be composed only of integers for the Poisson
          distribution to be used.", call. = FALSE)
      }
      
      if (!is.na(juvestimate)) {
        if (any(juvsizeb.data[, which(names(juvsizeb.data) == sizeb[1])] != 
            round(juvsizeb.data[, which(names(juvsizeb.data) == sizeb[1])]))) {
          stop("Secondary size variables must be composed only of integers for the Poisson
            distribution to be used.", call. = FALSE)
        }
      }
    } else if (sizebdist == "negbin") {
      if (any(sizeb.data[, which(names(sizeb.data) == sizeb[1])] != 
          round(sizeb.data[, which(names(sizeb.data) == sizeb[1])]))) {
        stop("Secondary size variables must be composed only of integers for the negative
          binomial distribution to be used.", call. = FALSE)
      }
      
      if (!is.na(juvestimate)) {
        if (any(juvsizeb.data[, which(names(juvsizeb.data) == sizeb[1])] != 
            round(juvsizeb.data[, which(names(juvsizeb.data) == sizeb[1])]))) {
          stop("Secondary size variables must be composed only of integers for the negative
            binomial distribution to be used.", call. = FALSE)
        }
      }
    }
  } else if (sizedist == "gamma" & sizeb_used) {
    if (any(sizeb.data[, which(names(sizeb.data) == sizeb[1])] < 0, na.rm = TRUE)) {
      stop("Secondary size variables must be non-negative for the gamma distribution to be used.",
        call. = FALSE)
    }
    
    if (!is.na(juvestimate)) {
      if (any(juvsizeb.data[, which(names(juvsizeb.data) == sizeb[1])]  < 0, na.rm = TRUE)) {
        stop("Secondary size variables must be non-negative for the gamma distribution to be used.",
          call. = FALSE)
      }
    }
  }
  
  if (!is.numeric(formulae$full.sizec.model) & sizec_used) {
    if (sizecdist == "poisson") {
      if (any(sizec.data[, which(names(sizec.data) == sizec[1])] != 
          round(sizec.data[, which(names(sizec.data) == sizec[1])]))) {
        stop("Tertiary size variables must be composed only of integers for the Poisson
          distribution to be used.", call. = FALSE)
      }
      
      if (!is.na(juvestimate)) {
        if (any(juvsizec.data[, which(names(juvsizec.data) == sizec[1])] != 
            round(juvsizec.data[, which(names(juvsizec.data) == sizec[1])]))) {
          stop("Tertiary size variables must be composed only of integers for the Poisson
            distribution to be used.", call. = FALSE)
        }
      }
    } else if (sizecdist == "negbin") {
      if (any(sizec.data[, which(names(sizec.data) == sizec[1])] != 
          round(sizec.data[, which(names(sizec.data) == sizec[1])]))) {
        stop("Tertiary size variables must be composed only of integers for the negative
          binomial distribution to be used.", call. = FALSE)
      }
      
      if (!is.na(juvestimate)) {
        if (any(juvsizec.data[, which(names(juvsizec.data) == sizec[1])] != 
            round(juvsizec.data[, which(names(juvsizec.data) == sizec[1])]))) {
          stop("Tertiary size variables must be composed only of integers for the negative
            binomial distribution to be used.", call. = FALSE)
        }
      }
    }
  } else if (sizedist == "gamma" & sizec_used) {
    if (any(sizec.data[, which(names(sizec.data) == sizec[1])] < 0, na.rm = TRUE)) {
      stop("Tertiary size variables must be non-negative for the gamma distribution to be used.",
        call. = FALSE)
    }
    
    if (!is.na(juvestimate)) {
      if (any(juvsizec.data[, which(names(juvsizec.data) == sizec[1])]  < 0, na.rm = TRUE)) {
        stop("Tertiary size variables must be non-negative for the gamma distribution to be used.",
          call. = FALSE)
      }
    }
  }
  
  if (is.element("fec", vitalrates)) {
    if (fectime == 2) {
      usedfec <- which(names(fec.data) == fec[2])
    } else if (fectime == 3) {
      usedfec <- which(names(fec.data) == fec[1])
    }
  }
  
  if (fecdist == "poisson" & !is.numeric(formulae$full.fec.model)) {
    
    if (any(fec.data[, usedfec] != round(fec.data[, usedfec]))) {
      stop("Fecundity variables must be composed only of integers for the Poisson
        distribution to be used.", call. = FALSE)
    }
  } else if (fecdist == "negbin" & !is.numeric(formulae$full.fec.model)) {
    if (any(fec.data[, usedfec] != round(fec.data[, usedfec]))) {
      stop("Fecundity variables must be composed only of integers for the negative
        binomial distribution to be used.", call. = FALSE)
    }
  } else if (fecdist == "gamma" & any(fec.data[, usedfec] < 0, na.rm = TRUE)) {
    stop("Fecundity variables must be non-negative for the gamma distribution to be used.",
      call. = FALSE)
  }
  
  #Now we run the modeling exercises
  if (is.numeric(formulae$full.surv.model)) {surv.global.model <- formulae$full.surv.model}
  if (is.numeric(formulae$full.obs.model)) {obs.global.model <- formulae$full.obs.model}
  if (is.numeric(formulae$full.size.model)) {size.global.model <- formulae$full.size.model}
  if (is.numeric(formulae$full.sizeb.model)) {sizeb.global.model <- formulae$full.sizeb.model}
  if (is.numeric(formulae$full.sizec.model)) {sizec.global.model <- formulae$full.sizec.model}
  if (is.numeric(formulae$full.repst.model)) {repst.global.model <- formulae$full.repst.model}
  if (is.numeric(formulae$full.fec.model)) {fec.global.model <- formulae$full.fec.model}
  
  if (is.numeric(formulae$juv.surv.model)) {juv.surv.global.model <- formulae$juv.surv.model}
  if (is.numeric(formulae$juv.obs.model)) {juv.obs.global.model <- formulae$juv.obs.model}
  if (is.numeric(formulae$juv.size.model)) {juv.size.global.model <- formulae$juv.size.model}
  if (is.numeric(formulae$juv.sizeb.model)) {juv.sizeb.global.model <- formulae$juv.sizeb.model}
  if (is.numeric(formulae$juv.sizec.model)) {juv.sizec.global.model <- formulae$juv.sizec.model}
  if (is.numeric(formulae$juv.repst.model)) {juv.repst.global.model <- formulae$juv.repst.model}
  
  surv.table <- NA
  obs.table <- NA
  size.table <- NA
  sizeb.table <- NA
  sizec.table <- NA
  repst.table <- NA
  fec.table <- NA
  
  juvsurv.table <- NA
  juvobs.table <- NA
  juvsize.table <- NA
  juvsizeb.table <- NA
  juvsizec.table <- NA
  juvrepst.table <- NA
  
  surv.bf <- NA
  obs.bf <- NA
  size.bf <- NA
  sizeb.bf <- NA
  sizec.bf <- NA
  repst.bf <- NA
  fec.bf <- NA
  
  juvsurv.bf <- NA
  juvobs.bf <- NA
  juvsize.bf <- NA
  juvsizeb.bf <- NA
  juvsizec.bf <- NA
  juvrepst.bf <- NA
  
  #A few more corrections to the model structure, used in running the global models
  correction.indiv <- gsub("individ", indiv, " + (1 | individ)", fixed = TRUE)
  if(approach == "mixed" & year.as.random) {
    correction.year <- gsub("yr", year, " + (1 | yr)", fixed = TRUE)
  } else {
    correction.year <- gsub("yr", year, " + yr", fixed = TRUE)
  }
  if (approach == "mixed" & patch.as.random) {
    correction.patch <- gsub("patch", patch, " + (1 | patch)", fixed = TRUE)
  } else {
    correction.patch <- gsub("patch", patch, " + patch", fixed = TRUE)
  }
  
  #Here we run the global models
  if (!is.numeric(formulae$full.surv.model) & nchar(formulae$full.surv.model) > 1) {
    
    if (is.element(0, surv.data$alive3) & is.element(1, surv.data$alive3)) {
      
      surv.global.list <- .headmaster_ritual(vrate = 1, approach = approach,
        dist = "binom", zero = FALSE, truncz = FALSE, quiet = quiet,
        usedformula = formulae$full.surv.model, subdata = surv.data,
        vind = surv.ind, vtrans = surv.trans, suite = suite,
        global.only = global.only, criterion = used.criterion,
        bestfit = bestfit, correction.patch, correction.year, correction.indiv,
        altformula = alt.formulae$full.surv.model)
      
      surv.global.model <- surv.global.list$model
      surv.ind <- surv.global.list$ind
      surv.trans <- surv.global.list$trans
      surv.table <- surv.global.list$table
      surv.bf <- surv.global.list$bf.model
      
      surv.accuracy <- .accu_predict(surv.bf, surv.data, surv[1])
      
    } else if (!is.element(0, surv.data$alive3)) {
      if (!quiet) {message("\nSurvival response is constant so will not model it.")}
      formulae$full.surv.model <- 1
      surv.global.model <- 1
      surv.bf <- 1
      surv.ind <- 0
      surv.trans <- 0
      surv.accuracy <- NA
    } else {
      if (!quiet) {message("\nSurvival response is constant so will not model it.")}
      formulae$full.surv.model <- 1
      surv.global.model <- 0
      surv.bf <- 0
      surv.ind <- 0
      surv.trans <- 0
      surv.accuracy <- NA
    }
  } else {
    surv.global.model <- 1
    surv.bf <- 1
    surv.ind <- 0
    surv.trans <- 0
    surv.accuracy <- NA
  }
  
  if (!is.numeric(formulae$full.obs.model) & nchar(formulae$full.obs.model) > 1) {
    if (is.element(0, obs.data$obsstatus3) & is.element(1, obs.data$obsstatus3)) {
        
      obs.global.list <- .headmaster_ritual(vrate = 2, approach = approach, 
        dist = "binom", zero = FALSE, truncz = FALSE, quiet = quiet, 
        usedformula = formulae$full.obs.model, subdata = obs.data,
        vind = obs.ind, vtrans = obs.trans, suite = suite,
        global.only = global.only, criterion = used.criterion,
        bestfit = bestfit, correction.patch, correction.year, correction.indiv,
        altformula = alt.formulae$full.obs.model)
      
      obs.global.model <- obs.global.list$model
      obs.ind <- obs.global.list$ind
      obs.trans <- obs.global.list$trans
      obs.table <- obs.global.list$table
      obs.bf <- obs.global.list$bf.model
      
      obs.accuracy <- .accu_predict(obs.bf, obs.data, obs[1])
      
    } else if (!is.element(0, obs.data$obsstatus3)) {
      if (!quiet) {message("\nObservation response is constant so will not model it.")}
      formulae$full.obs.model <- 1
      obs.global.model <- 1
      obs.bf <- 1
      obs.ind <- 0
      obs.trans <- 0
      obs.accuracy <- NA
    } else {
      if (!quiet) {message("\nObservation response is constant so will not model it.")}
      formulae$full.obs.model <- 1
      obs.global.model <- 0
      obs.bf <- 0
      obs.ind <- 0
      obs.trans <- 0
      obs.accuracy <- NA
    }
  } else {
    obs.global.model <- 1
    obs.bf <- 1
    obs.ind <- 0
    obs.trans <- 0
    obs.accuracy <- NA
  }
  
  if (!is.numeric(formulae$full.size.model) & nchar(formulae$full.size.model) > 1) {
    
    size.global.list <- .headmaster_ritual(vrate = 3, approach = approach, 
      dist = sizedist, zero = size.zero, truncz = size.trunc, quiet = quiet, 
      usedformula = formulae$full.size.model, subdata = size.data,
      vind = size.ind, vtrans = size.trans, suite = suite,
      global.only = global.only, criterion = used.criterion, bestfit = bestfit,
      correction.patch, correction.year, correction.indiv,
      altformula = alt.formulae$full.size.model)
    
    size.global.model <- size.global.list$model
    size.ind <- size.global.list$ind
    size.trans <- size.global.list$trans
    size.table <- size.global.list$table
    size.bf <- size.global.list$bf.model
    
    size.accuracy <- NA
    
  } else {
    size.global.model <- 1
    size.bf <- 1
    size.ind <- 0
    size.trans <- 0
    size.accuracy <- NA
  }
  
  if (!is.numeric(formulae$full.sizeb.model) & nchar(formulae$full.sizeb.model) > 1) {
    
    sizeb.global.list <- .headmaster_ritual(vrate = 10, approach = approach, 
      dist = sizebdist, zero = sizeb.zero, truncz = sizeb.trunc, quiet = quiet, 
      usedformula = formulae$full.sizeb.model, subdata = sizeb.data,
      vind = sizeb.ind, vtrans = sizeb.trans, suite = suite,
      global.only = global.only, criterion = used.criterion, bestfit = bestfit,
      correction.patch, correction.year, correction.indiv,
      altformula = alt.formulae$full.sizeb.model)
    
    sizeb.global.model <- sizeb.global.list$model
    sizeb.ind <- sizeb.global.list$ind
    sizeb.trans <- sizeb.global.list$trans
    sizeb.table <- sizeb.global.list$table
    sizeb.bf <- sizeb.global.list$bf.model
    
    sizeb.accuracy <- NA
    
  } else {
    sizeb.global.model <- 1
    sizeb.bf <- 1
    sizeb.ind <- 0
    sizeb.trans <- 0
    sizeb.accuracy <- NA
  }
  
  if (!is.numeric(formulae$full.sizec.model) & nchar(formulae$full.sizec.model) > 1) {
    
    sizec.global.list <- .headmaster_ritual(vrate = 11, approach = approach, 
      dist = sizecdist, zero = sizec.zero, truncz = sizec.trunc, quiet = quiet, 
      usedformula = formulae$full.sizec.model, subdata = sizec.data,
      vind = sizec.ind, vtrans = sizec.trans, suite = suite,
      global.only = global.only, criterion = used.criterion, bestfit = bestfit,
      correction.patch, correction.year, correction.indiv,
      altformula = alt.formulae$full.sizec.model)
    
    sizec.global.model <- sizec.global.list$model
    sizec.ind <- sizec.global.list$ind
    sizec.trans <- sizec.global.list$trans
    sizec.table <- sizec.global.list$table
    sizec.bf <- sizec.global.list$bf.model
    
    sizec.accuracy <- NA
    
  } else {
    sizec.global.model <- 1
    sizec.bf <- 1
    sizec.ind <- 0
    sizec.trans <- 0
    sizec.accuracy <- NA
  }
  
  if (!is.numeric(formulae$full.repst.model) & nchar(formulae$full.repst.model) > 1) {
    if (is.element(0, repst.data$repstatus3) & is.element(1, repst.data$repstatus3)) {
      repst.global.list <- .headmaster_ritual(vrate = 4, approach = approach, 
        dist = "binom", zero = FALSE, truncz = FALSE, quiet = quiet, 
        usedformula = formulae$full.repst.model, subdata = repst.data,
        vind = repst.ind, vtrans = repst.trans, suite = suite,
        global.only = global.only, criterion = used.criterion,
        bestfit = bestfit, correction.patch, correction.year, correction.indiv,
        altformula = alt.formulae$full.repst.model)
      
      repst.global.model <- repst.global.list$model
      repst.ind <- repst.global.list$ind
      repst.trans <- repst.global.list$trans
      repst.table <- repst.global.list$table
      repst.bf <- repst.global.list$bf.model
      
      repst.accuracy <- .accu_predict(repst.bf, repst.data, repst[1])
      
    } else if (!is.element(0, repst.data$repstatus3)) {
      if (!quiet) {message("\nReproductive status response is constant so will not model it.")}
      formulae$full.repst.model <- 1
      repst.global.model <- 1
      repst.bf <- 1
      repst.ind <- 0
      repst.trans <- 0
      repst.accuracy <- NA
    } else {
      if (!quiet) {message("\nReproductive status response is constant so will not model it.")}
      formulae$full.repst.model <- 1
      repst.global.model <- 0
      repst.bf <- 0
      repst.ind <- 0
      repst.trans <- 0
      repst.accuracy <- NA
    }
  } else {
    repst.global.model <- 1
    repst.bf <- 1
    repst.ind <- 0
    repst.trans <- 0
    repst.accuracy <- NA
  }
  
  if (!is.numeric(formulae$full.fec.model) & nchar(formulae$full.fec.model) > 1) {
    
    fec.global.list <- .headmaster_ritual(vrate = 5, approach = approach, 
      dist = fecdist, zero = fec.zero, truncz = fec.trunc, quiet = quiet, 
      usedformula = formulae$full.fec.model, subdata = fec.data,
      vind = fec.ind, vtrans = fec.trans, suite = suite,
      global.only = global.only, criterion = used.criterion, bestfit = bestfit,
      correction.patch, correction.year, correction.indiv,
      altformula = alt.formulae$full.fec.model)
    
    fec.global.model <- fec.global.list$model
    fec.ind <- fec.global.list$ind
    fec.trans <- fec.global.list$trans
    fec.table <- fec.global.list$table
    fec.bf <- fec.global.list$bf.model
    fec.accuracy <- NA
    
  } else {
    fec.global.model <- 1
    fec.bf <- 1
    fec.ind <- 0
    fec.trans <- 0
    fec.accuracy <- NA
  }
  
  if (!is.numeric(formulae$juv.surv.model) & nchar(formulae$juv.surv.model) > 1) {
    
    if (is.element(0, juvsurv.data$alive3) & is.element(1, juvsurv.data$alive3)) {
      
      juv.surv.global.list <- .headmaster_ritual(vrate = 6, approach = approach,
        dist = "binom", zero = FALSE, truncz = FALSE, quiet = quiet,
        usedformula = formulae$juv.surv.model, subdata = juvsurv.data,
        vind = juvsurv.ind, vtrans = juvsurv.trans, suite = suite,
        global.only = global.only, criterion = used.criterion,
        bestfit = bestfit, correction.patch, correction.year, correction.indiv,
        altformula = alt.formulae$juv.surv.model)
      
      juv.surv.global.model <- juv.surv.global.list$model
      juvsurv.ind <- juv.surv.global.list$ind
      juvsurv.trans <- juv.surv.global.list$trans
      juvsurv.table <- juv.surv.global.list$table
      juvsurv.bf <- juv.surv.global.list$bf.model
      
      juvsurv.accuracy <- .accu_predict(juvsurv.bf, juvsurv.data, surv[1])
        
    } else if (!is.element(0, juvsurv.data$alive3)) {
      if (!quiet) {message("\nJuvenile survival response is constant so will not model it.")}
      formulae$juv.surv.model <- 1
      juv.surv.global.model <- 1
      juvsurv.bf <- 1
      juvsurv.ind <- 0
      juvsurv.trans <- 0
      juvsurv.accuracy <- NA
    } else {
      if (!quiet) {message("\nJuvenile survival response is constant so will not model it.")}
      formulae$juv.surv.model <- 1
      juv.surv.global.model <- 0
      juvsurv.bf <- 0
      juvsurv.ind <- 0
      juvsurv.trans <- 0
      juvsurv.accuracy <- NA
    }
  } else {
    juv.surv.global.model <- 1
    juvsurv.bf <- 1
    juvsurv.ind <- 0
    juvsurv.trans <- 0
    juvsurv.accuracy <- NA
  }
  
  if (!is.numeric(formulae$juv.obs.model) & nchar(formulae$juv.obs.model) > 1) {
    
    if (is.element(0, juvobs.data$obsstatus3) & is.element(1, juvobs.data$obsstatus3)) {
      
      juv.obs.global.list <- .headmaster_ritual(vrate = 7, approach = approach,
        dist = "binom", zero = FALSE, truncz = FALSE, quiet = quiet,
        usedformula = formulae$juv.obs.model, subdata = juvobs.data,
        vind = juvobs.ind, vtrans = juvobs.trans, suite = suite,
        global.only = global.only, criterion = used.criterion,
        bestfit = bestfit, correction.patch, correction.year, correction.indiv,
        altformula = alt.formulae$juv.obs.model)
      
      juv.obs.global.model <- juv.obs.global.list$model
      juvobs.ind <- juv.obs.global.list$ind
      juvobs.trans <- juv.obs.global.list$trans
      juvobs.table <- juv.obs.global.list$table
      juvobs.bf <- juv.obs.global.list$bf.model
      
      juvobs.accuracy <- .accu_predict(juvobs.bf, juvobs.data, obs[1])
        
    } else if (!is.element(0, juvobs.data$obsstatus3)) {
      if (!quiet) {message("\nJuvenile observation response is constant so will not model it.")}
      formulae$juv.obs.model <- 1
      juv.obs.global.model <- 1
      juvobs.bf <- 1
      juvobs.ind <- 0
      juvobs.trans <- 0
      juvobs.accuracy <- NA
    } else {
      if (!quiet) {message("\nJuvenile observation response is constant so will not model it.")}
      formulae$juv.obs.model <- 1
      juv.obs.global.model <- 0
      juvobs.bf <- 0
      juvobs.ind <- 0
      juvobs.trans <- 0
      juvobs.accuracy <- NA
    }
  } else {
    juv.obs.global.model <- 1
    juvobs.bf <- 1
    juvobs.ind <- 0
    juvobs.trans <- 0
    juvobs.accuracy <- NA
  }
  
  if (!is.numeric(formulae$juv.size.model) & nchar(formulae$juv.size.model) > 1) {
    
    juv.size.global.list <- .headmaster_ritual(vrate = 8, approach = approach, 
      dist = sizedist, zero = jsize.zero, truncz = jsize.trunc, quiet = quiet, 
      usedformula = formulae$juv.size.model, subdata = juvsize.data,
      vind = juvsize.ind, vtrans = juvsize.trans, suite = suite,
      global.only = global.only, criterion = used.criterion, bestfit = bestfit,
      correction.patch, correction.year, correction.indiv,
      altformula = alt.formulae$juv.size.model)
    
    juv.size.global.model <- juv.size.global.list$model
    juvsize.ind <- juv.size.global.list$ind
    juvsize.trans <- juv.size.global.list$trans
    juvsize.table <- juv.size.global.list$table
    juvsize.bf <- juv.size.global.list$bf.model
    juvsize.accuracy <- NA
    
  } else {
    juv.size.global.model <- 1
    juvsize.bf <- 1
    juvsize.ind <- 0
    juvsize.trans <- 0
    juvsize.accuracy <- NA
  }
  
  if (!is.numeric(formulae$juv.sizeb.model) & nchar(formulae$juv.sizeb.model) > 1) {
    
    juv.sizeb.global.list <- .headmaster_ritual(vrate = 12, approach = approach, 
      dist = sizebdist, zero = jsizeb.zero, truncz = jsizeb.trunc, quiet = quiet, 
      usedformula = formulae$juv.sizeb.model, subdata = juvsizeb.data,
      vind = juvsizeb.ind, vtrans = juvsizeb.trans, suite = suite,
      global.only = global.only, criterion = used.criterion, bestfit = bestfit,
      correction.patch, correction.year, correction.indiv,
      altformula = alt.formulae$juv.sizeb.model)
    
    juv.sizeb.global.model <- juv.sizeb.global.list$model
    juvsizeb.ind <- juv.sizeb.global.list$ind
    juvsizeb.trans <- juv.sizeb.global.list$trans
    juvsizeb.table <- juv.sizeb.global.list$table
    juvsizeb.bf <- juv.sizeb.global.list$bf.model
    juvsizeb.accuracy <- NA
    
  } else {
    juv.sizeb.global.model <- 1
    juvsizeb.bf <- 1
    juvsizeb.ind <- 0
    juvsizeb.trans <- 0
    juvsizeb.accuracy <- NA
  }
  
  if (!is.numeric(formulae$juv.sizec.model) & nchar(formulae$juv.sizec.model) > 1) {
    
    juv.sizec.global.list <- .headmaster_ritual(vrate = 13, approach = approach, 
      dist = sizecdist, zero = jsizec.zero, truncz = jsizec.trunc, quiet = quiet, 
      usedformula = formulae$juv.sizec.model, subdata = juvsizec.data,
      vind = juvsizec.ind, vtrans = juvsizec.trans, suite = suite,
      global.only = global.only, criterion = used.criterion, bestfit = bestfit,
      correction.patch, correction.year, correction.indiv,
      altformula = alt.formulae$juv.sizec.model)
    
    juv.sizec.global.model <- juv.sizec.global.list$model
    juvsizec.ind <- juv.sizec.global.list$ind
    juvsizec.trans <- juv.sizec.global.list$trans
    juvsizec.table <- juv.sizec.global.list$table
    juvsizec.bf <- juv.sizec.global.list$bf.model
    juvsizec.accuracy <- NA
    
  } else {
    juv.sizec.global.model <- 1
    juvsizec.bf <- 1
    juvsizec.ind <- 0
    juvsizec.trans <- 0
    juvsizec.accuracy <- NA
  }
  
  if (!is.numeric(formulae$juv.repst.model) & nchar(formulae$juv.repst.model) > 1) {
    if (is.element(0, juvrepst.data$repstatus3) & is.element(1, juvrepst.data$repstatus3)) {
      juv.repst.global.list <- .headmaster_ritual(vrate = 9, approach = approach, 
        dist = "binom", zero = FALSE, truncz = FALSE, quiet = quiet, 
        usedformula = formulae$juv.repst.model, subdata = juvrepst.data,
        vind = juvrepst.ind, vtrans = juvrepst.trans, suite = suite,
        global.only = global.only, criterion = used.criterion,
        bestfit = bestfit, correction.patch, correction.year, correction.indiv,
        altformula = alt.formulae$juv.repst.model)
      
      juv.repst.global.model <- juv.repst.global.list$model
      juvrepst.ind <- juv.repst.global.list$ind
      juvrepst.trans <- juv.repst.global.list$trans
      juvrepst.table <- juv.repst.global.list$table
      juvrepst.bf <- juv.repst.global.list$bf.model
      
      juvrepst.accuracy <- .accu_predict(juvrepst.bf, juvrepst.data, repst[1])
      
    } else if (!is.element(0, juvrepst.data$repstatus3)) {
      if (!quiet) {message("\nJuvenile reproductive status response is constant so will not model it.")}
      formulae$juv.repst.model <- 1
      juv.repst.global.model <- 1
      juvrepst.bf <- 1
      juvrepst.ind <- 0
      juvrepst.trans <- 0
      juvrepst.accuracy <- NA
    } else {
      if (!quiet) {message("\nJuvenile reproductive status response is constant so will not model it.")}
      formulae$juv.repst.model <- 1
      juv.repst.global.model <- 0
      juvrepst.bf <- 0
      juvrepst.ind <- 0
      juvrepst.trans <- 0
      juvrepst.accuracy <- NA
    }
  } else {
    juv.repst.global.model <- 1
    juvrepst.bf <- 1
    juvrepst.ind <- 0
    juvrepst.trans <- 0
    juvrepst.accuracy <- NA
  }
  
  if (!quiet & global.only == FALSE) {message("\nFinished selecting best-fit models.\n")}
  
  qcoutput <- cbind.data.frame(
    c("survival", "observation", "size", "sizeb", "sizec", "reproduction", "fecundity",
      "juvenile_survival", "juvnile_observation", "juvenile_size", "juvenile_sizeb",
      "juvenile_sizec", "juvenile_reproduction"),
    c(surv.ind, obs.ind, size.ind, sizeb.ind, sizec.ind, repst.ind, fec.ind,
      juvsurv.ind, juvobs.ind, juvsize.ind, juvsizeb.ind, juvsizec.ind, juvrepst.ind),
    c(surv.trans, obs.trans, size.trans, sizeb.trans, sizec.trans, repst.trans, fec.trans,
      juvsurv.trans, juvobs.trans, juvsize.trans, juvsizeb.trans, juvsizec.trans, juvrepst.trans),
    c(surv.accuracy, obs.accuracy, size.accuracy, sizeb.accuracy, sizec.accuracy,
      repst.accuracy, fec.accuracy, juvsurv.accuracy, juvobs.accuracy, juvsize.accuracy, 
      juvsizeb.accuracy, juvsizec.accuracy, juvrepst.accuracy)
  )
  names(qcoutput) <- c("vital_rate", "individuals", "transitions", "accuracy")
  
  if (show.model.tables == FALSE) {
    surv.table <- NA
    obs.table <- NA
    size.table <- NA
    sizeb.table <- NA
    sizec.table <- NA
    repst.table <- NA
    fec.table <- NA
    
    juvsurv.table <- NA
    juvobs.table <- NA
    juvsize.table <- NA
    juvsizeb.table <- NA
    juvsizec.table <- NA
    juvrepst.table <- NA
  }
  
  if (global.only) {bestfit <- "global model only"}
  
  #Now we develop the final output, creating a new S3 class to do it
  full.output <- list(survival_model = surv.bf, observation_model = obs.bf,
    size_model = size.bf, sizeb_model = sizeb.bf, sizec_model = sizec.bf,
    repstatus_model = repst.bf, fecundity_model = fec.bf, 
    juv_survival_model = juvsurv.bf, juv_observation_model = juvobs.bf, 
    juv_size_model = juvsize.bf, juv_sizeb_model = juvsizeb.bf, 
    juv_sizec_model = juvsizec.bf, juv_reproduction_model = juvrepst.bf,
    survival_table = surv.table, observation_table = obs.table, 
    size_table = size.table, sizeb_table = sizeb.table,
    sizec_table = sizec.table, repstatus_table = repst.table,
    fecundity_table = fec.table, juv_survival_table = juvsurv.table, 
    juv_observation_table = juvobs.table, juv_size_table = juvsize.table,
    juv_sizeb_table = juvsizeb.table, juv_sizec_table = juvsizec.table,
    juv_reproduction_table = juvrepst.table, paramnames = formulae$paramnames, 
    criterion = bestfit, qc = qcoutput)
  class(full.output) <- "lefkoMod"
  
  return(full.output)
}

#' Full Model Selection for Single Vital Rate
#'
#' Function \code{.headmaster_ritual()} creates global model, dredges, and finds
#' the best-fit model given a core vital rate, dataset, and starting model
#' formula, and then passes it back to function \code{\link{modelsearch}()}.
#'
#' @param vrate An integer value indicating which model to build, as follows:
#' \code{1}: Adult survival probability, \code{2}: Adult observation
#' probability, \code{3}: Adult primary size, \code{4}: Adult reproduction
#' probability, \code{5}: Adult fecundity rate, \code{6}: Juvenile survival
#' probability, \code{7}: Juvenile observation probability, \code{9}: Juvenile
#' reproduction probability, \code{10}: Adult secondary size, \code{11}: Adult
#' tertiary size, \code{12}: Juvenile secondary size, and \code{13}: Juvenile
#' tertiary size.
#' @param approach A text variable indicating whether to use a mixed model
#' approach. If \code{"mixed"}, then assumes a mixed model approach. If
#' \code{"glm"}, will assume a GLM approach. Defaults to \code{"mixed"}.
#' @param dist Distribution of response. Used mostly in size and fecundity
#' models, which use \code{"gaussian"}, \code{"poisson"}, \code{"negbin"}, or
#' \code{"gamma"} distributions. Can also be set to \code{"binom"} for
#' probabilities.
#' @param zero A logical value indicating whether to use a zero-inflated
#' distribution. Only needed for non-Gaussian size or fecundity response models.
#' Defaults to \code{FALSE}.
#' @param truncz A logical value indicating whether to use a zero-truncated
#' distribution. Only needed for non-Gaussian size or fecundity response models.
#' Defaults to \code{FALSE}.
#' @param quiet A logical value indicating whether to issue diagnostic messages
#' while running. Defaults to \code{FALSE}.
#' @param usedformula The specific formula to use in modeling, entered as text
#' in list format.
#' @param subdata The data frame to use in modeling.
#' @param vind This is the individual number term previously developed for QC
#' output.
#' @param vtrans This is the transition number term previously developed for QC
#' output.
#' @param suite This describes the global model for each vital rate estimation
#' and has the following possible values: \code{full}, includes main effects and
#' all two-way interactions of size and reproductive status; \code{main},
#' includes main effects only of size and reproductive status; \code{size},
#' includes only size (also interactions between size in historical model);
#' \code{rep}, includes only reproductive status (also interactions between
#' status in historical model); \code{cons}, all vital rates estimated only as
#' y-intercepts. If \code{approach = "glm"} and \code{year.as.random = FALSE},
#' then year is also included as a fixed effect, and, in the case of
#' \code{full}, included in two-way interactions. Defaults to \code{size}.
#' @param global.only If set to TRUE, then only global models will be built and
#' evaluated. Defaults to FALSE.
#' @param criterion This refers to the model selection criterion used to rank
#' models when running the \code{\link[MuMIn]{dredge}()} function.
#' @param bestfit A variable indicating the model selection criterion for the
#' choice of best-fit model. The default is \code{AICc&k}, which chooses the 
#' best-fit model as the model with the lowest AICc or, if not the same model,
#' then the model that has the lowest degrees of freedom among models with
#' \eqn{\Delta AICc <= 2.0}. Alternatively, \code{AICc} may be chosen, in which
#' case the best-fit model is simply the model with the lowest AICc value.
#' @param correction.patch Text term denoting patch.
#' @param correction.year Text term denoting year.
#' @param correction.indiv Text term denoting individual identity.
#' @param altformula A list giving an alternative set of formulae to work with,
#' provided that \code{suite = "full"}. In other cases, a list of \code{NA}
#' values.
#'
#' @return Function \code{ms_binom()} outputs a list containing a global model,
#' the number of individuals and transitions used in modeling, the best-fit
#' model, and a dredge model table.
#' 
#' @keywords internal
#' @noRd
.headmaster_ritual <- function(vrate, approach = "mixed", dist = NA,
  zero = FALSE, truncz = FALSE, quiet = FALSE, usedformula, subdata, vind,
  vtrans, suite, global.only = FALSE, criterion = "AICc", bestfit = "AICc&k",
  correction.patch, correction.year, correction.indiv, altformula) {
  
  old <- options() #This function requires changes to options(na.action) in order for the lme4::dredge routines to work properly
  on.exit(options(old)) #This will reset options() to user originals when the function exits
  
  if (!is.na(dist)) {
    if (!is.element(dist, c("gaussian", "poisson", "negbin", "binom", "gamma"))) {
      stop("Response distribution not recognized.", call. = FALSE)
    }
  }
  
  if (vrate == 1) {
    if (!quiet) {message("\nDeveloping global model of survival probability...\n");}
    binom.model = TRUE
  } else if (vrate == 2) {
    if (!quiet) {message("\nDeveloping global model of observation probability...\n");}
    binom.model = TRUE
  } else if (vrate == 3) {
    if (!quiet) {message("\nDeveloping global model of size...\n");}
    binom.model = FALSE
  } else if (vrate == 4) {
    if (!quiet) {message("\nDeveloping global model of reproduction probability...\n");}
    binom.model = TRUE
  } else if (vrate == 5) {
    if (!quiet) {message("\nDeveloping global model of fecundity...\n");}
    binom.model = FALSE
  } else if (vrate == 6) {
    if (!quiet) {message("\nDeveloping global model of juvenile survival probability...\n");}
    binom.model = TRUE
  } else if (vrate == 7) {
    if (!quiet) {message("\nDeveloping global model of juvenile observation probability...\n");}
    binom.model = TRUE
  } else if (vrate == 8) {
    if (!quiet) {message("\nDeveloping global model of juvenile size...\n");}
    binom.model = FALSE
  } else if (vrate == 9) {
    if (!quiet) {message("\nDeveloping global model of juvenile reproduction probability...\n");}
    binom.model = TRUE
  } else if (vrate == 10) {
    if (!quiet) {message("\nDeveloping global model of secondary size...\n");}
    binom.model = FALSE
  } else if (vrate == 11) {
    if (!quiet) {message("\nDeveloping global model of tertiary size...\n");}
    binom.model = FALSE
  } else if (vrate == 12) {
    if (!quiet) {message("\nDeveloping global model of juvenile secondary size...\n");}
    binom.model = FALSE
  } else if (vrate == 13) {
    if (!quiet) {message("\nDeveloping global model of juvenile tertiary size...\n");}
    binom.model = FALSE
  } else {
    stop("Parameter to model not recognized.", call. = FALSE)
  }
  
  global.model <- .levindurosier(usedformula, subdata, approach, binom.model,
    dist, truncz, zero)
  
  if (any(class(global.model) == "try-error")) {
    if (!is.na(altformula)) {
      nox.model <- altformula
    } else {
      nox.model <- usedformula
    }
    
    if (nox.model != usedformula) {
      if (!quiet) {
        message("\nInitial global model estimation failed. Attempting a global model without interaction terms.\n")
      }
      global.model <- .levindurosier(nox.model, subdata, approach, binom.model,
        dist, truncz, zero)
    }
    
    if (any(class(global.model) == "try-error")) {
      nopat.model <- nox.model
      
      if (!is.na(correction.patch)) {
        nopat.model <- gsub(correction.patch, "", nopat.model, fixed = TRUE)
      }
      
      if (nox.model != nopat.model) {
        if (!quiet) {
          message("\nGlobal model estimation difficulties. Attempting a global model without a patch term.")
        }
        global.model <- .levindurosier(nopat.model, subdata, approach, binom.model,
          dist, truncz, zero)
      }
    }
    
    if (any(class(global.model) == "try-error")) {
      noyr.model <- nopat.model
      
      if (!is.na(correction.year)) {
        noyr.model <- gsub(correction.year, "", noyr.model, fixed = TRUE)
      }
      
      if (noyr.model != nopat.model) {
        if (!quiet) {
          message("\nGlobal model estimation difficulties. Attempting a global model without a year term.")
        }
        global.model <- .levindurosier(noyr.model, subdata, approach, binom.model,
          dist, truncz, zero)
      }
    }
    
    if (any(class(global.model) == "try-error")) {
      noind.model <- noyr.model
      
      if (!is.na(correction.indiv)) {
        noind.model <- gsub(correction.indiv, "", noind.model, fixed = TRUE)
      }
      
      if (noind.model != noyr.model) {
        if (!quiet) {
          message("\nGlobal model estimation difficulties. Attempting a global model without an individual identity term.")
        }
        global.model <- .levindurosier(noind.model, subdata, approach, binom.model,
          dist, truncz, zero)
      }
    }
    
    if (any(class(global.model) == "try-error")) {
      if (!quiet) {
        if (vrate == 1) {
          message("\nCould not properly estimate a global model for survival probability.")
        } else if (vrate == 2) {
          message("\nCould not properly estimate a global model for observation probability.")
        } else if (vrate == 3) {
          message("\nCould not properly estimate a global model for primary size.")
        } else if (vrate == 4) {
          message("\nCould not properly estimate a global model for reproduction probability.")
        } else if (vrate == 5) {
          message("\nCould not properly estimate a global model for fecundity.")
        } else if (vrate == 6) {
          message("\nCould not properly estimate a global model for juvenile survival probability.")
        } else if (vrate == 7) {
          message("\nCould not properly estimate a global model for juvenile observation probability.")
        } else if (vrate == 8) {
          message("\nCould not properly estimate a global model for juvenile size.")
        } else if (vrate == 9) {
          message("\nCould not properly estimate a global model for juvenile reproduction probability.")
        } else if (vrate == 10) {
          message("\nCould not properly estimate a global model for secondary primary size.")
        } else if (vrate == 11) {
          message("\nCould not properly estimate a global model for tertiary size.")
        } else if (vrate == 12) {
          message("\nCould not properly estimate a global model for juvenile secondary size.")
        } else if (vrate == 13) {
          message("\nCould not properly estimate a global model for juvenile tertiary size.")
        }
      }
      
      global.model <- 1
      vind <- 0
      vtrans <- 0
    }
  }
  
  if (vrate == 1) {
    if (!quiet) {message("\nGlobal model of survival probability developed. Proceeding with model dredge...\n");}
  } else if (vrate == 2) {
    if (!quiet) {message("\nGlobal model of observation probability developed. Proceeding with model dredge...\n");}
  } else if (vrate == 3) {
    if (!quiet) {message("\nGlobal model of primary size developed. Proceeding with model dredge...\n");}
  } else if (vrate == 4) {
    if (!quiet) {message("\nGlobal model of reproduction probability developed. Proceeding with model dredge...\n");}
  } else if (vrate == 5) {
    if (!quiet) {message("\nGlobal model of fecundity developed. Proceeding with model dredge...\n");}
  } else if (vrate == 6) {
    if (!quiet) {message("\nGlobal model of juvenile survival probability developed. Proceeding with model dredge...\n");}
  } else if (vrate == 7) {
    if (!quiet) {message("\nGlobal model of juvenile observation probability developed. Proceeding with model dredge...\n");}
  } else if (vrate == 8) {
    if (!quiet) {message("\nGlobal model of juvenile primary size developed. Proceeding with model dredge...\n");}
  } else if (vrate == 9) {
    if (!quiet) {message("\nGlobal model of juvenile reproduction probability developed. Proceeding with model dredge...\n");}
  } else if (vrate == 10) {
    if (!quiet) {message("\nGlobal model of secondary size developed. Proceeding with model dredge...\n");}
  } else if (vrate == 11) {
    if (!quiet) {message("\nGlobal model of tertiary size developed. Proceeding with model dredge...\n");}
  } else if (vrate == 12) {
    if (!quiet) {message("\nGlobal model of juvenile secondary size developed. Proceeding with model dredge...\n");}
  } else if (vrate == 13) {
    if (!quiet) {message("\nGlobal model of juvenile tertiary size developed. Proceeding with model dredge...\n");}
  }
  
  model.table <- NA
  
  #This is the section where we dredge the models
  if (suite != "cons" & global.only == FALSE) {

    if (usedformula != 1) {
      options(na.action = "na.fail")
      if (!quiet) {
        model.table <- try(MuMIn::dredge(global.model, rank = criterion), silent = TRUE)
      } else {
        model.table <- suppressWarnings(suppressMessages(try(MuMIn::dredge(global.model, 
                rank = criterion), silent = TRUE)))
      }
      
      if (any(class(model.table) == "try-error")) {
        if (vrate == 1) {
          message("\nDredge of survival probability failed.\n")
        } else if (vrate == 2) {
          message("\nDredge of observation probability failed.\n")
        } else if (vrate == 3) {
          message("\nDredge of primary size failed.\n")
        } else if (vrate == 4) {
          message("\nDredge of reproduction probability failed.\n")
        } else if (vrate == 5) {
          message("\nDredge of fecundity failed.\n")
        } else if (vrate == 6) {
          message("\nDredge of juvenile survival probability failed.\n")
        } else if (vrate == 7) {
          message("\nDredge of juvenile observation probability failed.\n")
        } else if (vrate == 8) {
          message("\nDredge of juvenile primary size failed.\n")
        } else if (vrate == 9) {
          message("\nDredge of juvenile reproduction probability failed.\n")
        } else if (vrate == 10) {
          message("\nDredge of secondary size failed.\n")
        } else if (vrate == 11) {
          message("\nDredge of tertiary size failed.\n")
        } else if (vrate == 12) {
          message("\nDredge of juvenile secondary size failed.\n")
        } else if (vrate == 13) {
          message("\nDredge of juvenile tertiary size failed.\n")
        }
      }
    }
  }
  
  #Here we extract the best-fit model
  if (stringr::str_detect(bestfit, "&k")) {
    if (any(class(model.table) == "model.selection")) {
      if (!quiet) {
        if (vrate == 1) {
          message("\nExtracting best-fit model of survival probability.\n")
        } else if (vrate == 2) {
          message("\nExtracting best-fit model of observation probability.\n")
        } else if (vrate == 3) {
          message("\nExtracting best-fit model of primary size.\n")
        } else if (vrate == 4) {
          message("\nExtracting best-fit model of reproduction probability.\n")
        } else if (vrate == 5) {
          message("\nExtracting best-fit model of fecundity.\n")
        } else if (vrate == 6) {
          message("\nExtracting best-fit model of juvenile survival probability.\n")
        } else if (vrate == 7) {
          message("\nExtracting best-fit model of juvenile observation probability.\n")
        } else if (vrate == 8) {
          message("\nExtracting best-fit model of juvenile primary size.\n")
        } else if (vrate == 9) {
          message("\nExtracting best-fit model of juvenile reproduction probability.\n")
        } else if (vrate == 10) {
          message("\nExtracting best-fit model of secondary size.\n")
        } else if (vrate == 11) {
          message("\nExtracting best-fit model of tertiary size.\n")
        } else if (vrate == 12) {
          message("\nExtracting best-fit model of juvenile secondary size.\n")
        } else if (vrate == 13) {
          message("\nExtracting best-fit model of juvenile tertiary size.\n")
        }
      }
      
      relevant.models <- which(model.table$delta <= 2)
      df.models <- which(model.table$df[relevant.models] == min(model.table$df[relevant.models]))
      if (length(df.models) > 1) {
        df.models <- min(df.models)
      }
      if (quiet) {
        model.bf <- suppressWarnings(suppressMessages(eval(stats::getCall(model.table, min(df.models)))))
      } else {
        model.bf <- eval(stats::getCall(model.table, min(df.models)))
      }
    } else {
      model.bf <- global.model
    }
    
  } else if (!stringr::str_detect(bestfit, "&k")) {
    if (any(class(model.table) == "model.selection")) {
      if (!quiet) {
        message("\nExtracting best-fit model.\n")
      }
      
      if (quiet) {
        model.bf <- suppressWarnings(suppressMessages(eval(stats::getCall(model.table, 1))))
      } else {
        model.bf <- eval(stats::getCall(model.table, 1))
      }
    } else {
      model.bf <- global.model
    }
  }
  
  output <- list(model = global.model, ind = vind, trans = vtrans,
    bf.model = model.bf, table = model.table)
  
  return(output)
}

#' Core Global Model Builder for .headmaster_ritual()
#' 
#' A function that gets used repeatedly in \code{\link{.headmaster_ritual}()} to
#' build the global model to be dredged in function \code{modelsearch()}.
#' 
#' @param usedformula The formula to be used in the linear modeling call.
#' @param subdata The data subset to be used in the linear modeling call.
#' @param approach Statistical approach, currently either "mixed" or "glm".
#' @param binom.model A logical value indicating whether to fit a binomial
#' distribution.
#' @param dist A string indicating whether the response is "gaussian",
#' "poisson", "negbin", or "gamma", if it is not binomial.
#' @param truncz A logical value indicating whether to use a zero-truncated
#' distribution.
#' @param zero A logical value indicating whether to use a zero-inflated
#' distribution.
#' 
#' @return This function returns a fit linear model of class generated by the
#' appropriate linear modeling function.
#' 
#' @keywords internal
#' @noRd
.levindurosier <- function(usedformula, subdata, approach, binom.model, dist,
  truncz, zero) {
  if (approach == "mixed" & !is.na(dist)) {
    if (binom.model) {
      global.model <- try(lme4::glmer(formula = stats::as.formula(usedformula), 
          data = subdata, family = "binomial"), silent = TRUE)
    } else {
      if (dist == "gaussian") {
        global.model <- try(lme4::lmer(formula = stats::as.formula(usedformula),
            data = subdata), silent = TRUE)
      } else if (dist == "gamma") {
          global.model <- try(lme4::glmer(formula = stats::as.formula(usedformula),
              data = subdata, family = "Gamma"), silent = TRUE)
      } else if (!truncz) {
        if (dist == "poisson" & !zero) {
          global.model <- try(lme4::glmer(formula = stats::as.formula(usedformula),
              data = subdata, family = "poisson"), silent = TRUE)
        } else if (dist == "poisson" & zero) {
          global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(usedformula),
              data = subdata, ziformula=~., family = "poisson"), silent = TRUE)
        } else if (dist == "negbin" & !zero) {
          global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(usedformula),
              data = subdata, ziformula=~0, family = glmmTMB::nbinom2), silent = TRUE)
        } else if (dist == "negbin" & zero) {
          global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(usedformula),
              data = subdata, ziformula=~., family = glmmTMB::nbinom2), silent = TRUE)
        }
      } else if (truncz) {
        if (dist == "poisson" & !zero) {
          global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(usedformula),
              data = subdata, family = glmmTMB::truncated_poisson), silent = TRUE)
        } else if (dist == "negbin" & !zero) {
          global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(usedformula),
              data = subdata, ziformula=~0, family = glmmTMB::truncated_nbinom2), silent = TRUE)
        }
      }
    }
  } else if (approach == "glm" & !is.na(dist)) {
    if (binom.model) {
      global.model <- try(stats::glm(formula = stats::as.formula(usedformula),
          data = subdata, family = "binomial"), silent = TRUE)
    } else {
      if (dist == "gaussian") {
        global.model <- try(stats::lm(formula = stats::as.formula(usedformula),
            data = subdata), silent = TRUE)
      } else if (dist == "gamma") {
        global.model <- try(stats::glm(formula = stats::as.formula(usedformula),
            data = subdata, family = "Gamma"), silent = TRUE)
      }else if (!truncz) {
        if (dist == "poisson" & !zero) {
          global.model <- try(stats::glm(formula = stats::as.formula(usedformula),
              data = subdata, family = "poisson"), silent = TRUE)
        } else if (dist == "poisson" & zero) {
          global.model <- try(pscl::zeroinfl(formula = stats::as.formula(usedformula),
              data = subdata, dist = "poisson"), silent = TRUE)
        } else if (dist == "negbin" & !zero) {
          global.model <- try(MASS::glm.nb(formula = stats::as.formula(usedformula),
              data = subdata), silent = TRUE)
        } else if (dist == "negbin" & zero) {
          global.model <- try(pscl::zeroinfl(formula = stats::as.formula(usedformula), 
              data = subdata, dist = "negbin"), silent = TRUE)
        }
      } else if (truncz) {
        usedformula <- gsub(" + 1", "", usedformula, fixed = TRUE)
        
        if (dist == "poisson" & !zero) {
          global.model <- try(VGAM::vglm(formula = stats::as.formula(usedformula),
              data = subdata, family = VGAM::pospoisson()), silent = TRUE)
        } else if (dist == "negbin" & !zero) {
          global.model <- try(VGAM::vglm(formula = stats::as.formula(usedformula),
              data = subdata, family = VGAM::posnegbinomial()), silent = TRUE)
        }
      }
    }
  } else {
    stop("Modeling approach not recognized.", call. = FALSE)
  }
}

#' Estimate Accuracy of Binomial Model
#' 
#' Function \code{.accu_predict} estimates the accuracy of vital rate models. It
#' currently only handles this for binomial models, and performs it as a
#' comparison of the actual and predicted responses from the respective model.
#' 
#' @param bestfitmodel The best-fit model to be passed.
#' @param givendata The dataset to be used for accuracy testing.
#' @param param The name of the response parameter to be used in accuracy
#' testing.
#' 
#' @return A single numeric value giving the proportion of responses accurately
#' predicted.
#' 
#' @keywords internal
#' @noRd
.accu_predict <- function(bestfitmodel, givendata, param) {
  pred_vec <- NULL
  
  accuracy <- NA 
  
  if (!is.numeric(bestfitmodel)) {
    pred_vec <- stats::predict(bestfitmodel, newdata = givendata,
      type = "response")
    pred_vec <- round(pred_vec)
    
    test_vec <- givendata[,which(names(givendata) == param)]
    
    results_vec <- pred_vec - test_vec
    
    accuracy <- length(which(results_vec == 0)) / length(results_vec)
  }
  
  return(accuracy)
}

#' Summary of Class "lefkoMod"
#' 
#' A function to summarize objects of class \code{lefkoMod}. This function shows
#' the best-fit models, summarizes the numbers of models in the model tables,
#' shows the criterion used to determine the best-fit models, and provides some
#' basic quality control information.
#' 
#' @param object An R object of class \code{lefkoMod} resulting from
#' \code{\link{modelsearch}()}.
#' @param ... Other parameters.
#' 
#' @return A summary of the object, showing the best-fit models for all vital
#' rates, with constants of 0 or 1 used for unestimated models. This is followed
#' by a summary of the number of models tested per vital rate, and a table
#' showing the names of the parameters used to model vital rates and represent
#' tested factors. At the end is a section describing the numbers of individuals
#' and of individual transitions used to estimate each vital rate best-fit
#' model, along with the accuracy of each binomial model.
#' 
#' @examples
#' \donttest{
#' # Lathyrus example
#' data(lathyrus)
#' 
#' sizevector <- c(0, 4.6, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8,
#'   9)
#' stagevector <- c("Sd", "Sdl", "Dorm", "Sz1nr", "Sz2nr", "Sz3nr", "Sz4nr", 
#'   "Sz5nr", "Sz6nr", "Sz7nr", "Sz8nr", "Sz9nr", "Sz1r", "Sz2r", "Sz3r", 
#'   "Sz4r", "Sz5r", "Sz6r", "Sz7r", "Sz8r", "Sz9r")
#' repvector <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#'   0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 4.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
#'   0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
#' 
#' lathframeln <- sf_create(sizes = sizevector, stagenames = stagevector, 
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
#'   propstatus = propvector)
#' 
#' lathvertln <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, 
#'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9,
#'   juvcol = "Seedling1988", sizeacol = "lnVol88", repstracol = "Intactseed88",
#'   fecacol = "Intactseed88", deadacol = "Dead1988", 
#'   nonobsacol = "Dormant1988", stageassign = lathframeln,
#'   stagesize = "sizea", censorcol = "Missing1988", censorkeep = NA,
#'   NAas0 = TRUE, censor = TRUE)
#' 
#' lathvertln$feca2 <- round(lathvertln$feca2)
#' lathvertln$feca1 <- round(lathvertln$feca1)
#' lathvertln$feca3 <- round(lathvertln$feca3)
#' 
#' lathmodelsln2 <- modelsearch(lathvertln, historical = FALSE, 
#'   approach = "mixed", suite = "main",
#'   vitalrates = c("surv", "obs", "size", "repst", "fec"), juvestimate = "Sdl",
#'   bestfit = "AICc&k", sizedist = "gaussian", fecdist = "poisson",
#'   indiv = "individ", patch = "patchid", year = "year2",
#'   year.as.random = TRUE, patch.as.random = TRUE, show.model.tables = TRUE,
#'   quiet = TRUE)
#' 
#' summary(lathmodelsln2)
#' }
#' 
#' @export
summary.lefkoMod <- function(object, ...) {
  
  modelsuite <- object

  totalmodels <- length(which(c(is.numeric(modelsuite$survival_model),
    is.numeric(modelsuite$observation_model), is.numeric(modelsuite$size_model),
    is.numeric(modelsuite$sizeb_model), is.numeric(modelsuite$sizec_model),
    is.numeric(modelsuite$repstatus_model), is.numeric(modelsuite$fecundity_model),
    is.numeric(modelsuite$juv_survival_model), is.numeric(modelsuite$juv_observation_model),
    is.numeric(modelsuite$juv_size_model), is.numeric(modelsuite$juv_sizeb_model),
    is.numeric(modelsuite$juv_sizec_model), is.numeric(modelsuite$juv_reproduction_model)) == FALSE))
  
  writeLines(paste0("This LefkoMod object includes ", totalmodels, " linear models."))
  writeLines(paste0("Best-fit model criterion used: ", modelsuite$criterion))
  writeLines(paste0("\n\n\n"))
  writeLines("Survival model:")
  print(modelsuite$survival_model)
  writeLines(paste0("\n\n"))
  writeLines("\nObservation model:")
  print(modelsuite$observation_model)
  writeLines(paste0("\n\n"))
  writeLines("\nSize model:")
  print(modelsuite$size_model)
  writeLines(paste0("\n\n"))
  writeLines("\nSecondary size model:")
  print(modelsuite$sizeb_model)
  writeLines(paste0("\n\n"))
  writeLines("\nTertiary size model:")
  print(modelsuite$sizec_model)
  writeLines(paste0("\n\n"))
  writeLines("\nReproductive status model:")
  print(modelsuite$repstatus_model)
  writeLines(paste0("\n\n"))
  writeLines("\nFecundity model:")
  print(modelsuite$fecundity_model)
  writeLines(paste0("\n\n\n"))
  writeLines("Juvenile survival model:")
  print(modelsuite$juv_survival_model)
  writeLines(paste0("\n\n"))
  writeLines("\nJuvenile observation model:")
  print(modelsuite$juv_observation_model)
  writeLines(paste0("\n\n"))
  writeLines("\nJuvenile size model:")
  print(modelsuite$juv_size_model)
  writeLines(paste0("\n\n"))
  writeLines("\nJuvenile secondary size model:")
  print(modelsuite$juv_sizeb_model)
  writeLines(paste0("\n\n"))
  writeLines("\nJuvenile tertiary size model:")
  print(modelsuite$juv_sizec_model)
  writeLines(paste0("\n\n"))
  writeLines("\nJuvenile reproduction model:")
  print(modelsuite$juv_reproduction_model)
  
  writeLines(paste0("\n\n\n"))
  if (!is.logical(modelsuite$survival_model) & 
      is.element("model.selection", class(modelsuite$survival_table))) {
    writeLines(paste0("\nNumber of models in survival table: ", 
        dim(modelsuite$survival_table)[1]))
  } else if (!is.logical(modelsuite$survival_model)) {
    writeLines("\nNumber of models in survival table: 1")
  } else {
    writeLines("\nSurvival table not estimated")
  }
  
  if (!is.logical(modelsuite$observation_model) & 
      is.element("model.selection", class(modelsuite$observation_table))) {
    writeLines(paste0("\nNumber of models in observation table: ", 
        dim(modelsuite$observation_table)[1]))
  } else if (!is.logical(modelsuite$observation_model)) {
    writeLines("\nNumber of models in observation table: 1")
  } else {
    writeLines("\nObservation table not estimated")
  }
  
  if (!is.logical(modelsuite$size_model) & 
      is.element("model.selection", class(modelsuite$size_table))) {
    writeLines(paste0("\nNumber of models in size table: ", 
        dim(modelsuite$size_table)[1]))
  } else if (!is.logical(modelsuite$size_model)) {
    writeLines("\nNumber of models in size table: 1")
  } else {
    writeLines("\nSize table not estimated")
  }
  
  if (!is.logical(modelsuite$sizeb_model) & 
      is.element("model.selection", class(modelsuite$sizeb_table))) {
    writeLines(paste0("\nNumber of models in secondary size table: ", 
        dim(modelsuite$sizeb_table)[1]))
  } else if (!is.logical(modelsuite$sizeb_model)) {
    writeLines("\nNumber of models in secondary size table: 1")
  } else {
    writeLines("\nSecondary size table not estimated")
  }
  
  if (!is.logical(modelsuite$sizec_model) & 
      is.element("model.selection", class(modelsuite$sizec_table))) {
    writeLines(paste0("\nNumber of models in tertiary size table: ", 
        dim(modelsuite$sizec_table)[1]))
  } else if (!is.logical(modelsuite$sizec_model)) {
    writeLines("\nNumber of models in tertiary size table: 1")
  } else {
    writeLines("\nTertiary size table not estimated")
  }
  
  if (!is.logical(modelsuite$repstatus_model) & 
      is.element("model.selection", class(modelsuite$repstatus_table))) {
    writeLines(paste0("\nNumber of models in reproduction status table: ", 
        dim(modelsuite$repstatus_table)[1]))
  } else if (!is.logical(modelsuite$repstatus_model)) {
    writeLines("\nNumber of models in reproduction status table: 1")
  } else {
    writeLines("\nReproduction status table not estimated")
  }
  
  if (!is.logical(modelsuite$fecundity_model) & 
      is.element("model.selection", class(modelsuite$fecundity_table))) {
    writeLines(paste0("\nNumber of models in fecundity table: ", 
        dim(modelsuite$fecundity_table)[1]))
  } else if (!is.logical(modelsuite$fecundity_model)) {
    writeLines("\nNumber of models in fecundity table: 1")
  } else {
    writeLines("\nFecundity table not estimated")
  }
  
  if (!is.logical(modelsuite$juv_survival_model) & 
      is.element("model.selection", class(modelsuite$juv_survival_table))) {
    writeLines(paste0("\nNumber of models in juvenile survival table: ", 
        dim(modelsuite$juv_survival_table)[1]))
  } else if (!is.logical(modelsuite$juv_survival_model)) {
    writeLines("\nNumber of models in juvenile survival table: 1")
  } else {
    writeLines("\nJuvenile survival table not estimated")
  }
  
  if (!is.logical(modelsuite$juv_observation_model) & 
      is.element("model.selection", class(modelsuite$juv_observation_table))) {
    writeLines(paste0("\nNumber of models in juvenile observation table: ", 
        dim(modelsuite$juv_observation_table)[1]))
  } else if (!is.logical(modelsuite$juv_observation_model)) {
    writeLines("\nNumber of models in juvenile observation table: 1")
  } else {
    writeLines("\nJuvenile observation table not estimated")
  }
  
  if (!is.logical(modelsuite$juv_size_model) & 
      is.element("model.selection", class(modelsuite$juv_size_table))) {
    writeLines(paste0("\nNumber of models in juvenile size table: ", 
        dim(modelsuite$juv_size_table)[1]))
  } else if (!is.logical(modelsuite$juv_size_model)) {
    writeLines("\nNumber of models in juvenile size table: 1")
  } else {
    writeLines("\nJuvenile size table not estimated")
  }
  
  if (!is.logical(modelsuite$juv_sizeb_model) & 
      is.element("model.selection", class(modelsuite$juv_sizeb_table))) {
    writeLines(paste0("\nNumber of models in juvenile secondary size table: ", 
        dim(modelsuite$juv_sizeb_table)[1]))
  } else if (!is.logical(modelsuite$juv_sizeb_model)) {
    writeLines("\nNumber of models in juvenile secondary size table: 1")
  } else {
    writeLines("\nJuvenile secondary size table not estimated")
  }
  
  if (!is.logical(modelsuite$juv_sizec_model) & 
      is.element("model.selection", class(modelsuite$juv_sizec_table))) {
    writeLines(paste0("\nNumber of models in juvenile tertiary size table: ", 
        dim(modelsuite$juv_sizec_table)[1]))
  } else if (!is.logical(modelsuite$juv_sizec_model)) {
    writeLines("\nNumber of models in juvenile tertiary size table: 1")
  } else {
    writeLines("\nJuvenile tertiary size table not estimated")
  }
  
  if (!is.logical(modelsuite$juv_reproduction_model) & 
      is.element("model.selection", class(modelsuite$juv_reproduction_table))) {
    writeLines(paste0("\nNumber of models in juvenile reproduction table: ", 
        dim(modelsuite$juv_reproduction_table)[1]))
  } else if (!is.logical(modelsuite$juv_reproduction_model)) {
    writeLines("\nNumber of models in juvenile reproduction table: 1")
  } else {
    writeLines("\nJuvenile reproduction table not estimated")
  }
  
  writeLines(paste0("\n\n\n"))
  writeLines("\nGeneral model parameter names (column 1), and ")
  writeLines("specific names used in these models (column 2): ")
  print.data.frame(modelsuite$paramnames[c(1,2)])
  
  writeLines(paste0("\n\n\n"))
  writeLines("\nQuality control:\n")
  if (modelsuite$qc[1,2] > 0) {
    writeLines(paste0("Survival estimated with ", modelsuite$qc[1,2], 
        " individuals and ", modelsuite$qc[1,3], " individual transitions."))
    writeLines(paste0("Survival accuracy is ", round(modelsuite$qc[1,4], 3), "."))
  } else {
    writeLines("Survival not estimated.")
  }
  
  if (modelsuite$qc[2,2] > 0) {
    writeLines(paste0("Observation estimated with ", modelsuite$qc[2,2], 
        " individuals and ", modelsuite$qc[2,3], " individual transitions."))
    writeLines(paste0("Observation accuracy is ", round(modelsuite$qc[2,4], 3), "."))
  } else {
    writeLines("Observation probability not estimated.")
  }
  
  if (modelsuite$qc[3,2] > 0) {
    writeLines(paste0("Size estimated with ", modelsuite$qc[3,2], 
        " individuals and ", modelsuite$qc[3,3], " individual transitions."))
  } else {
    writeLines("Size transition not estimated.")
  }
  
  if (modelsuite$qc[4,2] > 0) {
    writeLines(paste0("Secondary size estimated with ", modelsuite$qc[4,2], 
        " individuals and ", modelsuite$qc[4,3], " individual transitions."))
  } else {
    writeLines("Secondary size transition not estimated.")
  }
  
  if (modelsuite$qc[5,2] > 0) {
    writeLines(paste0("tertiry ize estimated with ", modelsuite$qc[5,2], 
        " individuals and ", modelsuite$qc[5,3], " individual transitions."))
  } else {
    writeLines("Tertiary size transition not estimated.")
  }
  
  if (modelsuite$qc[6,2] > 0) {
    writeLines(paste0("Reproductive status estimated with ", modelsuite$qc[6,2], 
        " individuals and ", modelsuite$qc[6,3], " individual transitions."))
   writeLines(paste0("Reproductive status accuracy is ", round(modelsuite$qc[6,4], 3), "."))
   } else {
    writeLines("Reproduction probability not estimated.")
  }
  
  if (modelsuite$qc[7,2] > 0) {
    writeLines(paste0("Fecundity estimated with ", modelsuite$qc[7,2], 
        " individuals and ", modelsuite$qc[7,3], " individual transitions."))
  } else {
    writeLines("Fecundity not estimated.")
  }
  
  if (modelsuite$qc[8,2] > 0) {
    writeLines(paste0("Juvenile survival estimated with ", modelsuite$qc[8,2], 
        " individuals and ", modelsuite$qc[8,3], " individual transitions."))
   writeLines(paste0("Juvenile survival accuracy is ", round(modelsuite$qc[8,4], 3), "."))
   } else {
    writeLines("Juvenile survival not estimated.")
  }
  
  if (modelsuite$qc[9,2] > 0) {
    writeLines(paste0("Juvenile observation estimated with ", modelsuite$qc[9,2], 
        " individuals and ", modelsuite$qc[9,3], " individual transitions."))
   writeLines(paste0("Juvenile observation accuracy is ", round(modelsuite$qc[9,4], 3), "."))
  } else {
    writeLines("Juvenile observation probability not estimated.")
  }
  
  if (modelsuite$qc[10,2] > 0) {
    writeLines(paste0("Juvenile size estimated with ", modelsuite$qc[10,2], 
        " individuals and ", modelsuite$qc[10,3], " individual transitions."))
  } else {
    writeLines("Juvenile size transition not estimated.")
  }
  
  if (modelsuite$qc[11,2] > 0) {
    writeLines(paste0("Juvenile secondary size estimated with ", modelsuite$qc[11,2],
        " individuals and ", modelsuite$qc[11,3], " individual transitions."))
  } else {
    writeLines("Juvenile secondary size transition not estimated.")
  }
  
  if (modelsuite$qc[12,2] > 0) {
    writeLines(paste0("Juvenile tertiary size estimated with ", modelsuite$qc[12,2],
        " individuals and ", modelsuite$qc[12,3], " individual transitions."))
  } else {
    writeLines("Juvenile tertiary size transition not estimated.")
  }
  
  if (modelsuite$qc[13,2] > 0) {
    writeLines(paste0("Juvenile reproduction estimated with ", modelsuite$qc[13,2], " individuals and ", modelsuite$qc[13,3], " individual transitions."))
   writeLines(paste0("Juvenile reproductive status accuracy is ", round(modelsuite$qc[13,4], 3), "."))
  } else {
    writeLines("Juvenile reproduction probability not estimated.")
  }
}

#' Extracts Year, Patch, and Other Terms from lefkoMod Objects
#' 
#' Function \code{.mystic_rhythms()} extracts year and patch terms from
#' \code{lefkoMod} objects.
#' 
#' @param model The actual \code{lefkoMod} object.
#' @param fix_zi A logical value indicating whether to use the conditional or
#' zero-inflated portion of the model output. Defaults to \code{TRUE}, in which
#' case the conditional is assumed.
#' @param model_code The integer index of the model.
#' @param random_switch A logical value indicating whether to treat the term as
#' random. Defaults to \code{FALSE}.
#' @param yp_switch An integer indicating whether to extract year or patch
#' terms, or something else.
#' @param paramnames The \code{paramnames} portion of the \code{lefkoMod}
#' object.
#' @param mainnames A vector of all times, patches, or other objects to be
#' extracted.
#' @param chosen_term A string giving the name of the term to be extracted from
#' the \code{paramnames} portion of the \code{lefkoMod} object. Defaults to
#' \code{NULL}.
#' 
#' @return A data frame showing the values of year or patch terms.
#' 
#' @keywords internal
#' @noRd
.mystic_rhythms <- function(model, fix_zi = TRUE, model_code,
  random_switch = FALSE, yp_switch, paramnames, mainnames,
  chosen_term = NULL) {
  
  yp_coefs <- yp_names <- yp_trace <- NULL
  
  if (yp_switch == 1) {
    extracted_term <- "year2"
  } else if (yp_switch == 2) {
    extracted_term <- "patch"
  } else if (yp_switch == 3) {
    extracted_term <- chosen_term
  } else {
    stop("Error in year/patch term extraction.", call. = FALSE)
  }
  
  if (fix_zi) {
    model_focus <- "cond"
  } else {
    model_focus <- "zi"
  }
  
  if (length(mainnames) == 1) {
    if (paramnames[(which(paramnames$mainparams == extracted_term)), "modelparams"] == "none") {
      yp_coefs <- as.data.frame(c(0))
      names(yp_coefs) <- extracted_term
      
      return(yp_coefs)
    }
  }
   
  if (random_switch) {
    if (model_code == 1) {
      check_model <- glmmTMB::fixef(model)[[model_focus]]
      
      if (length(check_model) > 0) {
        yp_coefs <- glmmTMB::ranef(model)[[model_focus]][paramnames[(which(paramnames$mainparams == extracted_term)), "modelparams"]][[1]]
        yp_names <- names(glmmTMB::ranef(model)[[model_focus]][paramnames[(which(paramnames$mainparams == extracted_term)), "modelparams"]])
      }
    } else if (model_code < 4 & fix_zi) {
      yp_coefs <- lme4::ranef(model)[paramnames[(which(paramnames$mainparams == extracted_term)), "modelparams"]][[1]]
      yp_names <- names(lme4::ranef(model)[paramnames[(which(paramnames$mainparams == extracted_term)), "modelparams"]][[1]])
    } else if (model_code < 4 & !fix_zi) {
      yp_coefs <- as.data.frame(rep(0, length(mainnames)))
      rownames(yp_coefs) <- mainnames
      names(yp_coefs) <- paramnames[(which(paramnames$mainparams == extracted_term)), "modelparams"]
    } else if (model_code > 7) {
      if (!all(is.null(mainnames))) {
        yp_coefs <- as.data.frame(rep(0, length(mainnames)))
        rownames(yp_coefs) <- mainnames
      } else if (!all(is.na(mainnames))) {
        yp_coefs <- as.data.frame(rep(0, length(mainnames)))
        rownames(yp_coefs) <- mainnames
      } else {
        yp_coefs <- as.data.frame(c(0))
      }
      names(yp_coefs) <- paramnames[(which(paramnames$mainparams == extracted_term)), "modelparams"]
      yp_names <- rownames(yp_coefs)
    } else {
      stop("Random terms are only allowed in glmmTMB and lme4 objects.", call. = FALSE)
    }
  } else {
    if (model_code == 1) {
      check_model <- glmmTMB::fixef(model)[[model_focus]]
      
      if (length(check_model) > 0) {
        yp_trace <- apply(as.matrix(names(glmmTMB::fixef(model)[[model_focus]])), 1, function(X) {
          grepl(paramnames[(which(paramnames$mainparams == extracted_term)), "modelparams"], X)
        })
        yp_coefs <- glmmTMB::fixef(model)[[model_focus]][yp_trace]
        yp_names <- names(glmmTMB::fixef(model)[[model_focus]][yp_trace])
      }
    } else if (model_code < 4) {
      if (fix_zi) {
        yp_trace <- apply(as.matrix(names(lme4::fixef(model))), 1, function(X) {
          grepl(paramnames[(which(paramnames$mainparams == extracted_term)), "modelparams"], X)
        })
        yp_coefs <- lme4::fixef(model)[yp_trace]
        yp_names <- names(lme4::fixef(model)[yp_trace])     
      } else {
        yp_coefs <- as.data.frame(rep(0, length(mainnames)))
        rownames(yp_coefs) <- mainnames
        names(yp_coefs) <- paramnames[(which(paramnames$mainparams == extracted_term)), "modelparams"]
        yp_names <- rownames(yp_coefs)
      }
    } else if (model_code == 4) {
      if (fix_zi) {
        yp_trace <- apply(as.matrix(names(model$coefficients$count)), 1, function(X) {
          grepl(paramnames[(which(paramnames$mainparams == extracted_term)), "modelparams"], X)
        })
        yp_coefs <- model$coefficients$count[yp_trace]
        yp_names <- names(model$coefficients$count[yp_trace])
      } else {
        yp_trace <- apply(as.matrix(names(model$coefficients$zero)), 1, function(X) {
          grepl(paramnames[(which(paramnames$mainparams == extracted_term)), "modelparams"], X)
        })
        yp_coefs <- model$coefficients$zero[yp_trace]
        yp_names <- names(model$coefficients$zero[yp_trace])
      }
    } else if (model_code == 5) {
      if (fix_zi) {
        yp_trace <- apply(as.matrix(names(model@coefficients)), 1, function(X) {
          grepl(paramnames[(which(paramnames$mainparams == extracted_term)), "modelparams"], X)
        })
        yp_coefs <- model@coefficients[yp_trace]
        yp_names <- names(model@coefficients[yp_trace])
      } else {
        yp_coefs <- as.data.frame(rep(0, length(mainnames)))
        rownames(yp_coefs) <- mainnames
        names(yp_coefs) <- paramnames[(which(paramnames$mainparams == extracted_term)), "modelparams"]
        yp_names <- rownames(yp_coefs)
      }
    } else if (model_code < 8) {
      if (fix_zi) {
        yp_trace <- apply(as.matrix(names(model$coefficients)), 1, function(X) {
          grepl(paramnames[(which(paramnames$mainparams == extracted_term)), "modelparams"], X)
        })
        yp_coefs <- model$coefficients[yp_trace]
        yp_names <- names(model$coefficients[yp_trace])
      } else {
        yp_coefs <- as.data.frame(rep(0, length(mainnames)))
        rownames(yp_coefs) <- mainnames
        names(yp_coefs) <- paramnames[(which(paramnames$mainparams == extracted_term)), "modelparams"]
        yp_names <- rownames(yp_coefs)
      }
    } else if (model_code < 10) {
      yp_coefs <- as.data.frame(rep(0, length(mainnames)))
      names(yp_coefs) <- paramnames[(which(paramnames$mainparams == extracted_term)), "modelparams"]
      rownames(yp_coefs) <- mainnames
      yp_names <- rownames(yp_coefs)
      if (names(yp_coefs) == "none") names(yp_coefs) <- extracted_term
    } else {
      stop("Error in fixed year/patch/other extraction.", call. = FALSE)
    }
  }
  
  if (!all(is.na(yp_coefs)) & !random_switch) {
    yp_names <- as.vector(apply(as.matrix(yp_names), 1, function(X) {
      a <- gsub(paramnames[(which(paramnames$mainparams == extracted_term)), "modelparams"], "", X)
      b <- gsub("as.factor\\(\\)", "", a)
      
      return(b)
    }))
    if (!is.data.frame(yp_coefs)) {
      names(yp_coefs) <- yp_names
      yp_coefs <- as.data.frame(yp_coefs)
      names(yp_coefs) <- paramnames[(which(paramnames$mainparams == extracted_term)), "modelparams"]
    }
  } else if (all(is.na(yp_coefs)) | all(is.null(yp_coefs)))  {
    yp_coefs <- as.data.frame(rep(0, length(mainnames)))
    rownames(yp_coefs) <- mainnames
    names(yp_coefs) <- paramnames[(which(paramnames$mainparams == extracted_term)), "modelparams"]
    if (names(yp_coefs) == "none") names(yp_coefs) <- extracted_term
  }
  
  return(yp_coefs)
}

#' Main Fixed Effect Extraction Function
#' 
#' Function \code{.handinglove()} extracts fixed factor coefficients from input
#' models.
#' 
#' @param fixed_model The fixed, conditional, or zero-inflated fixed portion of
#' the model output.
#' @param model_code An integer signifying which class of linear model to
#' assume. Use: \code{1} for \code{glmmTMB}, \code{2} for \code{lmerMod},
#' \code{3} for \code{merMod}, \code{4} for \code{zeroinfl}, \code{5} for
#' \code{vglm}.
#' @param paramnames The \code{paramnames} portion of the \code{lefkoMod} object
#' from which the model is derived.
#' @param factor1 The name of the primary factor to extract.
#' @param factor2 If extracting an interaction, then the name of the second
#' factor involved in the interaction.
#' @param interaction A logical value indicating whether to test for an
#' interaction between \code{factor1} and \code{factor2}.
#' @param err_check If \code{TRUE}, then leaves missing factor values as
#' \code{NA}.
#' 
#' @return A single numeric value corresponding to the factor coefficient, or a
#' \code{0} value otherwise. If \code{err_check = TRUE}, then any conditions
#' leading to \code{NA} values return \code{NA} values.
#' 
#' @keywords internal
#' @noRd
.handinglove <- function(fixed_model, model_code, paramnames, factor1,
  factor2 = NULL, interaction = FALSE, err_check = FALSE) {
  
  first_check <- second_check <- NULL
  
  if (model_code < 8) {
    if (!interaction) {
      if (factor1 != "(Intercept)") {
        first_check <- fixed_model[paramnames[(which(paramnames$mainparams == factor1)), "modelparams"]]
      } else {
        if (model_code == 5) {
          first_check <- fixed_model["(Intercept):1"]
        } else {
          first_check <- fixed_model["(Intercept)"]
        }
      }
      
      
    } else {
      first_check <- fixed_model[paste0(paramnames[(which(paramnames$mainparams == factor1)),
          "modelparams"], ":", paramnames[(which(paramnames$mainparams == factor2)), "modelparams"])]
      second_check <- fixed_model[paste0(paramnames[(which(paramnames$mainparams == factor2)),
          "modelparams"], ":", paramnames[(which(paramnames$mainparams == factor1)), "modelparams"])]
    }
  } else {
    if (model_code > 7 & factor1 == "(Intercept)") {
      first_check <- fixed_model
    } else {
      first_check <- NA
      second_check <- NA
    }
  }
  
  if (is.numeric(first_check) & length(first_check) > 0) {
    if (!is.na(first_check)) {
      return(first_check)
    } else if (is.numeric(second_check) & length(second_check) > 0) {
      if(!is.na(second_check)) {
        return(second_check)
      } else {
        return(0)
      }
    } else {
      return(0)
    }
  } else if (is.numeric(second_check) & length(second_check) > 0) {
    return(second_check)
  } else if (err_check) {
    return(NA)
  } else {
    return(0)
  }
}

#' Imputes Missing Year and Patch Values
#' 
#' Function \code{.moosewhistle()} takes a vector of year, patch, or other
#' values, finds which terms are missing, and imputes them either with \code{0}
#' values or with random values.
#' 
#' @param coefvec The vector of year, patch, or other coefficients
#' @param random_switch A logical value indicating whether term is assumed to be
#' a random factor.
#' @param mainnames A vector of all term names in the correct order.
#' 
#' @return An updated version of \code{coefvec} with missing values imputed and
#' inserted in the correct order.
#' 
#' @keywords internal
#' @noRd
.moosewhistle <- function(coefvec, random_switch, mainnames) {
  
  if (!all(is.null(coefvec)) & !all(is.na(mainnames))) {
    namesdiff <- setdiff(mainnames, rownames(coefvec))
    
    if (length(namesdiff) > 0) {
      if (random_switch) {
        newdevs <- rnorm(length(namesdiff), mean = 0, sd = sd(coefvec[,1]))
      } else {
        newdevs <- rep(0, (length(namesdiff)))
      }
      updated_df <- as.data.frame(as.matrix(newdevs, ncol = 1))
      names(updated_df) <- names(coefvec)
      rownames(updated_df) <- namesdiff
      
      updated_df <- rbind.data.frame(updated_df, coefvec)
      updated_df <- updated_df[with(updated_df, order(rownames(updated_df))),, drop = FALSE]
      
      updated_vec <- updated_df
    } else {
      updated_vec <- coefvec
    }
  } else if (!all(is.na(mainnames))) {
    updated_vec <- as.data.frame(rep(0, length(mainnames)))
    names(updated_vec) <- names(coefvec)
  } else {
    updated_vec <- as.data.frame(c(0))
    names(updated_vec) <- names(coefvec)
  }
  
  return(updated_vec)
}

#' Extract Coefficients From glmmTMB-estimated Linear Vital Rate Models
#' 
#' Function \code{.modelextract()} extracts coefficient values from linear
#' models estimated through the \code{\link[glmmTMB]{glmmTMB}()} function, to
#' estimate vital rates in \code{'lefko3'}. Used to supply coefficients to
#' \code{\link{flefko3}()}, \code{\link{flefko2}()}, and
#' \code{\link{aflefko2}()}.
#' 
#' @param model A linear model estimated through one of the methods used in
#' function \code{\link{modelsearch}()}.
#' @param paramnames Data frame giving the names of standard coefficients
#' required by matrix creation functions.
#' @param mainyears A vector of the names of the monitoring occasions.
#' @param mainpatches A vector of the names of the patches.
#' @param maingroups A vector of the names of all stage groups.
#' @param mainindcova A vector denoting values of individual covariate \code{a}.
#' @param mainindcovb A vector denoting values of individual covariate \code{b}.
#' @param mainindcovc A vector denoting values of individual covariate \code{c}.
#' @param year.as.random A logical value indicating whether to treat years as
#' random factors.
#' @param patch.as.random A logical value indicating whether to treat patches as
#' random factors.
#' @param inda.as.random A logical value indicating whether to treat individual
#' covariate \code{a} as random and categorical. Defaults to \code{FALSE}.
#' @param indb.as.random A logical value indicating whether to treat individual
#' covariate \code{b} as random and categorical. Defaults to \code{FALSE}.
#' @param indc.as.random A logical value indicating whether to treat individual
#' covariate \code{c} as random and categorical. Defaults to \code{FALSE}.
#' @param err_check A logical value indicating whether to allow \code{NA} values
#' to pass into the model extracted output. Otherwise, all \code{NA} values are
#' converted to \code{0} values.
#' 
#' @return This function returns a list with the following elements:
#' \item{coefficients}{Vector of fixed effect coefficients.}
#' \item{years}{Vector of occasion coefficients, typically random.}
#' \item{patches}{Vector of patch coefficients, typically random.}
#' \item{groups2}{Vector of group coefficients for time t.}
#' \item{groups1}{Vector of group coefficients for time t-1.}
#' \item{indcova2}{Vector of individual covariate \code{a} values for time t.}
#' \item{indcova1}{Vector of individual covariate \code{a} values for time t-1.}
#' \item{indcovb2}{Vector of individual covariate \code{b} values for time t.}
#' \item{indcovb1}{Vector of individual covariate \code{b} values for time t-1.}
#' \item{indcovc2}{Vector of individual covariate \code{c} values for time t.}
#' \item{indcovc1}{Vector of individual covariate \code{c} values for time t-1.}
#' \item{zeroyear}{Vector of zero-inflated occasion coefficients, typically
#' random.}
#' \item{zeropatch}{Vector of zero-inflated patch coefficients, typically
#' random.}
#' \item{zerogroups2}{Vector of zero-inflated group coefficients for time t.}
#' \item{zerogroups1}{Vector of zero-inflated group coefficients for time t-1.}
#' \item{zeroindcova2}{Vector of zero-inflated individual covariate \code{a}
#' values for time t.}
#' \item{zeroindcova1}{Vector of zero-inflated individual covariate \code{a}
#' values for time t-1.}
#' \item{zeroindcovb2}{Vector of zero-inflated individual covariate \code{b}
#' values for time t.}
#' \item{zeroindcovb1}{Vector of zero-inflated individual covariate \code{b}
#' values for time t-1.}
#' \item{zeroindcovc2}{Vector of zero-inflated individual covariate \code{c}
#' values for time t.}
#' \item{zeroindcovc1}{Vector of zero-inflated individual covariate \code{c}
#' values for time t-1.}
#' \item{variances}{Residual variance terms associated with model.}
#' \item{family}{Distribution of response term.}
#' \item{sigma}{Standard deviation or disperson parameter of response.}
#' \item{trunc}{A binomial value indicating whether the model is from a
#' truncated distribution (1) or not (0).}
#' 
#' @keywords internal
#' @noRd
.modelextract <- function(model, paramnames, mainyears, mainpatches, maingroups,
  mainindcova, mainindcovb, mainindcovc, year.as.random, patch.as.random,
  inda.as.random = FALSE, indb.as.random = FALSE, indc.as.random = FALSE,
  err_check = FALSE) {
  
  year.zi <- patch.zi <- year.zi.names <- patch.zi.names <- NULL
  
  model_type <- viewedfamily <- NULL
  rvars <- NA
  sigma <- NA
  trunc0 <- 0
  
  if (any(class(model) == "glmmTMB")) {
    model_type <- 1
    
    fix_model <- glmmTMB::fixef(model)[["cond"]]
    zi_model <- glmmTMB::fixef(model)[["zi"]]
    
    viewedfamily <- stats::family(model)$family
    
    if (viewedfamily == "truncated_nbinom2") {
      viewedfamily <- "nbinom2"
      trunc0 <- 1
    } else if (viewedfamily == "truncated_poisson") {
      viewedfamily <- "poisson"
      trunc0 <- 1
    }
    sigma <-  glmmTMB::sigma(model)
    rvars <- glmmTMB::VarCorr(model)
    
  } else if (any(class(model) == "lmerMod")) {
    model_type <- 2
    
    fix_model <- lme4::fixef(model)
    zi_model <- NA
    rvars <- as.data.frame(lme4::VarCorr(model))
    
    viewedfamily <- "gaussian"
    sigma <- attr(summary(model)$varcor, "sc")
    
  } else if (any(class(model) == "merMod") | any(class(model) == "glmerMod")) {
    model_type <- 3
    
    fix_model <- lme4::fixef(model)
    zi_model <- NA
    rvars <- as.data.frame(lme4::VarCorr(model))
    sigma <- stats::sigma(model)
    
    viewedfamily <- stats::family(model)$family
    
  } else if (any(class(model) == "zeroinfl")) {
    model_type <- 4
    
    fix_model <- model$coefficients$count
    zi_model <- model$coefficients$zero
    
    if (model$dist == "binomial") {
      viewedfamily <- "binomial"
      sigma <- 1
    }
    if (model$dist == "poisson") {
      viewedfamily <- "poisson"
      sigma <- 1
    }
    if (model$dist == "negbin") {
      viewedfamily <- "negbin"
      sigma <- model$theta
    }
    
  } else if (any(class(model) == "vglm")) {
    model_type <- 5
    
    fix_model <- model@coefficients
    zi_model <- NA
    
    if (model@family@vfamily == "posnegbinomial") {
      viewedfamily <- "negbin"
      sigma <- model@dispersion
    } else {
      sigma <- 1
    }
    
    if (model@family@vfamily == "pospoisson") {
      viewedfamily <- "poisson"
    }
    
  } else if (any(class(model) == "glm")) {
    model_type <- 7
    
    fix_model <- model$coefficients
    zi_model <- NA
    
    if (is.element("gaussian", class(model))) {
      viewedfamily <- "gaussian"
      sigma <- stats::sigma(model)
    }
    if (is.element("negbin", class(model))) {
      viewedfamily <- "negbin"
      sigma <- model$theta
    }
    if (model$family["family"] == "Gamma") {
      viewedfamily <- "gamma"
      sigma <- stats::sigma(model)
    }
    if (model$family["family"] == "binomial") {viewedfamily <- "binomial"}
    if (model$family["family"] == "poisson") {viewedfamily <- "poisson"}
    
  } else if (any(class(model) == "lm")) {
    model_type <- 6
    
    fix_model <- model$coefficients
    zi_model <- NA
    
    viewedfamily <- "gaussian"
    sigma <- stats::sigma(model)
    
  } else if (class(model) == "numeric") {
    model_type <- 8
    
    fix_model <- model
    zi_model <- NA
    rvars <- as.data.frame(0)
    
    viewedfamily <- "constant"
    family <- 1
    sigma <- 1
    
  } else if (class(model) == "logical") {
    model_type <- 9
    
    fix_model <- NA
    zi_model <- NA
    rvars <- as.data.frame(0)
    
    viewedfamily <- "constant"
    family <- 1
    sigma <- 1
    
  }
    
  year.coefs <- .mystic_rhythms(model, fix_zi = TRUE, model_code = model_type,
    random_switch = year.as.random, yp_switch = 1, paramnames = paramnames,
    mainnames = mainyears)
  year.zi <- .mystic_rhythms(model, fix_zi = FALSE, model_code = model_type,
    random_switch = year.as.random, yp_switch = 1, paramnames = paramnames,
    mainnames = mainyears)
  patch.coefs <- .mystic_rhythms(model, fix_zi = TRUE, model_code = model_type,
    random_switch = patch.as.random, yp_switch = 2, paramnames = paramnames,
    mainnames = mainpatches)
  patch.zi <- .mystic_rhythms(model, fix_zi = FALSE, model_code = model_type,
    random_switch = patch.as.random, yp_switch = 2, paramnames = paramnames,
    mainnames = mainpatches)
    
  if (model_type == 1 | model_type == 4) {
    year.zi <- .moosewhistle(year.zi, year.as.random, mainyears)
    patch.zi <- .moosewhistle(patch.zi, patch.as.random, mainpatches)
  } 
  year.coefs <- .moosewhistle(year.coefs, year.as.random, mainyears)
  patch.coefs <- .moosewhistle(patch.coefs, patch.as.random, mainpatches)
  
  group2.coefs <- .mystic_rhythms(model, fix_zi = TRUE, model_code = model_type,
    random_switch = FALSE, yp_switch = 3, paramnames = paramnames,
    mainnames = maingroups, chosen_term = "group2")
  group1.coefs <- .mystic_rhythms(model, fix_zi = TRUE, model_code = model_type,
    random_switch = FALSE, yp_switch = 3, paramnames = paramnames,
    mainnames = maingroups, chosen_term = "group1")
  group2.coefs <- .moosewhistle(group2.coefs, FALSE, maingroups)
  group1.coefs <- .moosewhistle(group1.coefs, FALSE, maingroups)
  
  group2.zi <- .mystic_rhythms(model, fix_zi = FALSE, model_code = model_type,
    random_switch = FALSE, yp_switch = 3, paramnames = paramnames,
    mainnames = maingroups, chosen_term = "group2")
  group1.zi <- .mystic_rhythms(model, fix_zi = FALSE, model_code = model_type,
    random_switch = FALSE, yp_switch = 3, paramnames = paramnames,
    mainnames = maingroups, chosen_term = "group1")
  group2.zi <- .moosewhistle(group2.zi, FALSE, maingroups)
  group1.zi <- .moosewhistle(group1.zi, FALSE, maingroups)
  
  indcova1.coefs <- .mystic_rhythms(model, fix_zi = TRUE, model_code = model_type,
    random_switch = inda.as.random, yp_switch = 3, paramnames = paramnames,
    mainnames = mainindcova, chosen_term = "indcova1")
  indcovb1.coefs <- .mystic_rhythms(model, fix_zi = TRUE, model_code = model_type,
    random_switch = indb.as.random, yp_switch = 3, paramnames = paramnames,
    mainnames = mainindcovb, chosen_term = "indcovb1")
  indcovc1.coefs <- .mystic_rhythms(model, fix_zi = TRUE, model_code = model_type,
    random_switch = indc.as.random, yp_switch = 3, paramnames = paramnames,
    mainnames = mainindcovc, chosen_term = "indcovc1")
  indcova2.coefs <- .mystic_rhythms(model, fix_zi = TRUE, model_code = model_type,
    random_switch = inda.as.random, yp_switch = 3, paramnames = paramnames,
    mainnames = mainindcova, chosen_term = "indcova2")
  indcovb2.coefs <- .mystic_rhythms(model, fix_zi = TRUE, model_code = model_type,
    random_switch = indb.as.random, yp_switch = 3, paramnames = paramnames,
    mainnames = mainindcovb, chosen_term = "indcovb2")
  indcovc2.coefs <- .mystic_rhythms(model, fix_zi = TRUE, model_code = model_type,
    random_switch = indc.as.random, yp_switch = 3, paramnames = paramnames,
    mainnames = mainindcovc, chosen_term = "indcovc2")
  
  indcova1.zi <- .mystic_rhythms(model, fix_zi = FALSE, model_code = model_type,
    random_switch = inda.as.random, yp_switch = 3, paramnames = paramnames,
    mainnames = mainindcova, chosen_term = "indcova1")
  indcovb1.zi <- .mystic_rhythms(model, fix_zi = FALSE, model_code = model_type,
    random_switch = indb.as.random, yp_switch = 3, paramnames = paramnames,
    mainnames = mainindcovb, chosen_term = "indcovb1")
  indcovc1.zi <- .mystic_rhythms(model, fix_zi = FALSE, model_code = model_type,
    random_switch = indc.as.random, yp_switch = 3, paramnames = paramnames,
    mainnames = mainindcovc, chosen_term = "indcovc1")
  indcova2.zi <- .mystic_rhythms(model, fix_zi = FALSE, model_code = model_type,
    random_switch = inda.as.random, yp_switch = 3, paramnames = paramnames,
    mainnames = mainindcova, chosen_term = "indcova2")
  indcovb2.zi <- .mystic_rhythms(model, fix_zi = FALSE, model_code = model_type,
    random_switch = indb.as.random, yp_switch = 3, paramnames = paramnames,
    mainnames = mainindcovb, chosen_term = "indcovb2")
  indcovc2.zi <- .mystic_rhythms(model, fix_zi = FALSE, model_code = model_type,
    random_switch = indc.as.random, yp_switch = 3, paramnames = paramnames,
    mainnames = mainindcovc, chosen_term = "indcovc2")
  
  if (is.null(viewedfamily)) {
    stop("Model distribution not recognized.", call. = FALSE)
  }
  
  yintercept.coef <-.handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "(Intercept)", factor2 = NULL, interaction = FALSE, err_check = err_check)
  flw2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "repst2", factor2 = NULL, interaction = FALSE, err_check = err_check)
  flw1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "repst1", factor2 = NULL, interaction = FALSE, err_check = err_check)
  flw1.flw2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "repst1", factor2 = "repst2", interaction = TRUE, err_check = err_check)

  size2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size2", factor2 = NULL, interaction = FALSE, err_check = err_check)
  size1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size1", factor2 = NULL, interaction = FALSE, err_check = err_check)
  size1.size2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size1", factor2 = "size2", interaction = TRUE, err_check = err_check)
  
  size1.flw1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size1", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  size2.flw2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size2", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  
  size2.flw1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size2", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  size1.flw2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size1", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  
  age.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "age", factor2 = NULL, interaction = FALSE, err_check = err_check)
  age.size1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size1", factor2 = "age", interaction = TRUE, err_check = err_check)
  age.size2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size2", factor2 = "age", interaction = TRUE, err_check = err_check)
  age.flw1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "repst1", factor2 = "age", interaction = TRUE, err_check = err_check)
  age.flw2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "repst2", factor2 = "age", interaction = TRUE, err_check = err_check)
  
  #Here are the old individual covariate coefficients
  inda2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = NULL, interaction = FALSE, err_check = err_check)
  indb2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb2", factor2 = NULL, interaction = FALSE, err_check = err_check)
  indc2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = NULL, interaction = FALSE, err_check = err_check)
  inda1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = NULL, interaction = FALSE, err_check = err_check)
  indb1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb1", factor2 = NULL, interaction = FALSE, err_check = err_check)
  indc1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = NULL, interaction = FALSE, err_check = err_check)
  
  inda2.size2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "size2", interaction = TRUE, err_check = err_check)
  indb2.size2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb2", factor2 = "size2", interaction = TRUE, err_check = err_check)
  indc2.size2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = "size2", interaction = TRUE, err_check = err_check)
  inda2.flw2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  indb2.flw2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb2", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  indc2.flw2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  
  inda1.size1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "size1", interaction = TRUE, err_check = err_check)
  indb1.size1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb1", factor2 = "size1", interaction = TRUE, err_check = err_check)
  indc1.size1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "size1", interaction = TRUE, err_check = err_check)
  inda1.flw1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  indb1.flw1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb1", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  indc1.flw1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  
  inda2.indb2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "indcovb2", interaction = TRUE, err_check = err_check)
  inda2.indc2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "indcovc2", interaction = TRUE, err_check = err_check)
  indb2.indc2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = "indcovb2", interaction = TRUE, err_check = err_check)
  inda1.indb1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "indcovb1", interaction = TRUE, err_check = err_check)
  inda1.indc1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "indcovc1", interaction = TRUE, err_check = err_check)
  indb1.indc1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "indcovb1", interaction = TRUE, err_check = err_check)
  
  inda2.indb1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "indcovb1", interaction = TRUE, err_check = err_check)
  inda1.indb2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "indcovb2", interaction = TRUE, err_check = err_check)
  inda2.indc1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "indcova2", interaction = TRUE, err_check = err_check)
  inda1.indc2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "indcovc2", interaction = TRUE, err_check = err_check)
  indb2.indc1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "indcovb2", interaction = TRUE, err_check = err_check)
  indb1.indc2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = "indcovb1", interaction = TRUE, err_check = err_check)
  
  #New conditional
  inda2.size1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "size1", interaction = TRUE, err_check = err_check)
  indb2.size1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb2", factor2 = "size1", interaction = TRUE, err_check = err_check)
  indc2.size1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = "size1", interaction = TRUE, err_check = err_check)
  inda2.flw1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  indb2.flw1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb2", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  indc2.flw1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  
  inda1.size2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "size2", interaction = TRUE, err_check = err_check)
  indb1.size2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb1", factor2 = "size2", interaction = TRUE, err_check = err_check)
  indc1.size2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "size2", interaction = TRUE, err_check = err_check)
  inda1.flw2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  indb1.flw2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb1", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  indc1.flw2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  
  sizeb2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizeb2", factor2 = NULL, interaction = FALSE, err_check = err_check)
  sizeb1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizeb1", factor2 = NULL, interaction = FALSE, err_check = err_check)
  sizec2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec2", factor2 = NULL, interaction = FALSE, err_check = err_check)
  sizec1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec1", factor2 = NULL, interaction = FALSE, err_check = err_check)
  
  dens.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "density", factor2 = NULL, interaction = FALSE, err_check = err_check)
  
  sizeb1.sizeb2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizeb1", factor2 = "sizeb2", interaction = TRUE, err_check = err_check)
  sizec1.sizec2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec1", factor2 = "sizec2", interaction = TRUE, err_check = err_check)
  sizea1.sizeb1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size1", factor2 = "sizeb1", interaction = TRUE, err_check = err_check)
  sizea1.sizec1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size1", factor2 = "sizec1", interaction = TRUE, err_check = err_check)
  sizeb1.sizec1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec1", factor2 = "sizeb1", interaction = TRUE, err_check = err_check)
  sizea2.sizeb2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size2", factor2 = "sizeb2", interaction = TRUE, err_check = err_check)
  sizea2.sizec2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size2", factor2 = "sizec2", interaction = TRUE, err_check = err_check)
  sizeb2.sizec2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec2", factor2 = "sizeb2", interaction = TRUE, err_check = err_check)
  sizea1.sizeb2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size1", factor2 = "sizeb2", interaction = TRUE, err_check = err_check)
  sizea1.sizec2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size1", factor2 = "sizec2", interaction = TRUE, err_check = err_check)
  sizeb1.sizec2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec2", factor2 = "sizeb1", interaction = TRUE, err_check = err_check)
  sizea2.sizeb1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size2", factor2 = "sizeb1", interaction = TRUE, err_check = err_check)
  sizea2.sizec1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size2", factor2 = "sizec1", interaction = TRUE, err_check = err_check)
  sizeb2.sizec1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec1", factor2 = "sizeb2", interaction = TRUE, err_check = err_check)
  
  sizea2.dens.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size2", factor2 = "density", interaction = TRUE, err_check = err_check)
  sizea1.dens.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size1", factor2 = "density", interaction = TRUE, err_check = err_check)
  sizeb2.dens.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizeb2", factor2 = "density", interaction = TRUE, err_check = err_check)
  sizeb1.dens.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizeb1", factor2 = "density", interaction = TRUE, err_check = err_check)
  sizec2.dens.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec2", factor2 = "density", interaction = TRUE, err_check = err_check)
  sizec1.dens.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec1", factor2 = "density", interaction = TRUE, err_check = err_check)
  flw1.dens.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "repst1", factor2 = "density", interaction = TRUE, err_check = err_check)
  flw2.dens.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "repst2", factor2 = "density", interaction = TRUE, err_check = err_check)
  
  sizeb2.flw2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizeb2", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  sizeb1.flw2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizeb1", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  sizec2.flw2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec2", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  sizec1.flw2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec1", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  sizeb2.flw1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizeb2", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  sizeb1.flw1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizeb1", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  sizec2.flw1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec2", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  sizec1.flw1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec1", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  
  sizeb2.age.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizeb2", factor2 = "age", interaction = TRUE, err_check = err_check)
  sizeb1.age.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizeb1", factor2 = "age", interaction = TRUE, err_check = err_check)
  sizec2.age.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec2", factor2 = "age", interaction = TRUE, err_check = err_check)
  sizec1.age.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec1", factor2 = "age", interaction = TRUE, err_check = err_check)
  dens.age.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "density", factor2 = "age", interaction = TRUE, err_check = err_check)
  
  inda2.sizeb2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "sizeb2", interaction = TRUE, err_check = err_check)
  inda2.sizec2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "sizec2", interaction = TRUE, err_check = err_check)
  inda2.sizeb1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "sizeb1", interaction = TRUE, err_check = err_check)
  inda2.sizec1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "sizec1", interaction = TRUE, err_check = err_check)
  inda1.sizeb2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "sizeb2", interaction = TRUE, err_check = err_check)
  inda1.sizec2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "sizec2", interaction = TRUE, err_check = err_check)
  inda1.sizeb1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "sizeb1", interaction = TRUE, err_check = err_check)
  inda1.sizec1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "sizec1", interaction = TRUE, err_check = err_check)
  inda2.dens.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "density", interaction = TRUE, err_check = err_check)
  inda1.dens.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "density", interaction = TRUE, err_check = err_check)
  
  indb2.sizeb2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb2", factor2 = "sizeb2", interaction = TRUE, err_check = err_check)
  indb2.sizec2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb2", factor2 = "sizec2", interaction = TRUE, err_check = err_check)
  indb2.sizeb1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb2", factor2 = "sizeb1", interaction = TRUE, err_check = err_check)
  indb2.sizec1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb2", factor2 = "sizec1", interaction = TRUE, err_check = err_check)
  indb1.sizeb2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb1", factor2 = "sizeb2", interaction = TRUE, err_check = err_check)
  indb1.sizec2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb1", factor2 = "sizec2", interaction = TRUE, err_check = err_check)
  indb1.sizeb1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb1", factor2 = "sizeb1", interaction = TRUE, err_check = err_check)
  indb1.sizec1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb1", factor2 = "sizec1", interaction = TRUE, err_check = err_check)
  indb2.dens.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb2", factor2 = "density", interaction = TRUE, err_check = err_check)
  indb1.dens.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb1", factor2 = "density", interaction = TRUE, err_check = err_check)
  
  indc2.sizeb2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = "sizeb2", interaction = TRUE, err_check = err_check)
  indc2.sizec2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = "sizec2", interaction = TRUE, err_check = err_check)
  indc2.sizeb1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = "sizeb1", interaction = TRUE, err_check = err_check)
  indc2.sizec1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = "sizec1", interaction = TRUE, err_check = err_check)
  indc1.sizeb2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "sizeb2", interaction = TRUE, err_check = err_check)
  indc1.sizec2.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "sizec2", interaction = TRUE, err_check = err_check)
  indc1.sizeb1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "sizeb1", interaction = TRUE, err_check = err_check)
  indc1.sizec1.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "sizec1", interaction = TRUE, err_check = err_check)
  indc2.dens.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = "density", interaction = TRUE, err_check = err_check)
  indc1.dens.coef <- .handinglove(fix_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "density", interaction = TRUE, err_check = err_check)
  
  #These are old zero-inflated coefficients
  yintercept.zi <-.handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "(Intercept)", factor2 = NULL, interaction = FALSE, err_check = err_check)
  flw2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "repst2", factor2 = NULL, interaction = FALSE, err_check = err_check)
  flw1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "repst1", factor2 = NULL, interaction = FALSE, err_check = err_check)
  flw1.flw2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "repst1", factor2 = "repst2", interaction = TRUE, err_check = err_check)

  size2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size2", factor2 = NULL, interaction = FALSE, err_check = err_check)
  size1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size1", factor2 = NULL, interaction = FALSE, err_check = err_check)
  size1.size2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size1", factor2 = "size2", interaction = TRUE, err_check = err_check)
  
  size1.flw1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size1", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  size2.flw2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size2", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  
  size2.flw1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size2", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  size1.flw2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size1", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  
  age.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "age", factor2 = NULL, interaction = FALSE, err_check = err_check)
  age.size1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size1", factor2 = "age", interaction = TRUE, err_check = err_check)
  age.size2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size2", factor2 = "age", interaction = TRUE, err_check = err_check)
  age.flw1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "repst1", factor2 = "age", interaction = TRUE, err_check = err_check)
  age.flw2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "repst2", factor2 = "age", interaction = TRUE, err_check = err_check)
  
  #Here are the individual covariate bits
  inda2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = NULL, interaction = FALSE, err_check = err_check)
  indb2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb2", factor2 = NULL, interaction = FALSE, err_check = err_check)
  indc2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = NULL, interaction = FALSE, err_check = err_check)
  inda1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = NULL, interaction = FALSE, err_check = err_check)
  indb1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb1", factor2 = NULL, interaction = FALSE, err_check = err_check)
  indc1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = NULL, interaction = FALSE, err_check = err_check)
  
  inda2.size2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "size2", interaction = TRUE, err_check = err_check)
  indb2.size2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb2", factor2 = "size2", interaction = TRUE, err_check = err_check)
  indc2.size2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = "size2", interaction = TRUE, err_check = err_check)
  inda2.flw2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  indb2.flw2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb2", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  indc2.flw2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  
  inda1.size1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "size1", interaction = TRUE, err_check = err_check)
  indb1.size1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb1", factor2 = "size1", interaction = TRUE, err_check = err_check)
  indc1.size1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "size1", interaction = TRUE, err_check = err_check)
  inda1.flw1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  indb1.flw1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb1", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  indc1.flw1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  
  inda2.indb2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "indcovb2", interaction = TRUE, err_check = err_check)
  inda2.indc2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "indcovc2", interaction = TRUE, err_check = err_check)
  indb2.indc2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = "indcovb2", interaction = TRUE, err_check = err_check)
  inda1.indb1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "indcovb1", interaction = TRUE, err_check = err_check)
  inda1.indc1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "indcovc1", interaction = TRUE, err_check = err_check)
  indb1.indc1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "indcovb1", interaction = TRUE, err_check = err_check)
  
  inda2.indb1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "indcovb1", interaction = TRUE, err_check = err_check)
  inda1.indb2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "indcovb2", interaction = TRUE, err_check = err_check)
  inda2.indc1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "indcova2", interaction = TRUE, err_check = err_check)
  inda1.indc2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "indcovc2", interaction = TRUE, err_check = err_check)
  indb2.indc1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "indcovb2", interaction = TRUE, err_check = err_check)
  indb1.indc2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = "indcovb1", interaction = TRUE, err_check = err_check)
  
  #New
  inda2.size1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "size1", interaction = TRUE, err_check = err_check)
  indb2.size1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb2", factor2 = "size1", interaction = TRUE, err_check = err_check)
  indc2.size1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = "size1", interaction = TRUE, err_check = err_check)
  inda2.flw1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  indb2.flw1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb2", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  indc2.flw1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  
  inda1.size2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "size2", interaction = TRUE, err_check = err_check)
  indb1.size2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb1", factor2 = "size2", interaction = TRUE, err_check = err_check)
  indc1.size2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "size2", interaction = TRUE, err_check = err_check)
  inda1.flw2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  indb1.flw2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb1", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  indc1.flw2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  
  sizeb2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizeb2", factor2 = NULL, interaction = FALSE, err_check = err_check)
  sizeb1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizeb1", factor2 = NULL, interaction = FALSE, err_check = err_check)
  sizec2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec2", factor2 = NULL, interaction = FALSE, err_check = err_check)
  sizec1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec1", factor2 = NULL, interaction = FALSE, err_check = err_check)
  
  dens.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "density", factor2 = NULL, interaction = FALSE, err_check = err_check)
  
  sizeb1.sizeb2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizeb1", factor2 = "sizeb2", interaction = TRUE, err_check = err_check)
  sizec1.sizec2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec1", factor2 = "sizec2", interaction = TRUE, err_check = err_check)
  sizea1.sizeb1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size1", factor2 = "sizeb1", interaction = TRUE, err_check = err_check)
  sizea1.sizec1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size1", factor2 = "sizec1", interaction = TRUE, err_check = err_check)
  sizeb1.sizec1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec1", factor2 = "sizeb1", interaction = TRUE, err_check = err_check)
  sizea2.sizeb2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size2", factor2 = "sizeb2", interaction = TRUE, err_check = err_check)
  sizea2.sizec2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size2", factor2 = "sizec2", interaction = TRUE, err_check = err_check)
  sizeb2.sizec2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec2", factor2 = "sizeb2", interaction = TRUE, err_check = err_check)
  sizea1.sizeb2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size1", factor2 = "sizeb2", interaction = TRUE, err_check = err_check)
  sizea1.sizec2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size1", factor2 = "sizec2", interaction = TRUE, err_check = err_check)
  sizeb1.sizec2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec2", factor2 = "sizeb1", interaction = TRUE, err_check = err_check)
  sizea2.sizeb1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size2", factor2 = "sizeb1", interaction = TRUE, err_check = err_check)
  sizea2.sizec1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size2", factor2 = "sizec1", interaction = TRUE, err_check = err_check)
  sizeb2.sizec1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec1", factor2 = "sizeb2", interaction = TRUE, err_check = err_check)
  
  sizea2.dens.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size2", factor2 = "density", interaction = TRUE, err_check = err_check)
  sizea1.dens.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "size1", factor2 = "density", interaction = TRUE, err_check = err_check)
  sizeb2.dens.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizeb2", factor2 = "density", interaction = TRUE, err_check = err_check)
  sizeb1.dens.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizeb1", factor2 = "density", interaction = TRUE, err_check = err_check)
  sizec2.dens.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec2", factor2 = "density", interaction = TRUE, err_check = err_check)
  sizec1.dens.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec1", factor2 = "density", interaction = TRUE, err_check = err_check)
  flw1.dens.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "repst1", factor2 = "density", interaction = TRUE, err_check = err_check)
  flw2.dens.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "repst2", factor2 = "density", interaction = TRUE, err_check = err_check)
  
  sizeb2.flw2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizeb2", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  sizeb1.flw2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizeb1", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  sizec2.flw2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec2", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  sizec1.flw2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec1", factor2 = "repst2", interaction = TRUE, err_check = err_check)
  sizeb2.flw1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizeb2", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  sizeb1.flw1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizeb1", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  sizec2.flw1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec2", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  sizec1.flw1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec1", factor2 = "repst1", interaction = TRUE, err_check = err_check)
  
  sizeb2.age.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizeb2", factor2 = "age", interaction = TRUE, err_check = err_check)
  sizeb1.age.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizeb1", factor2 = "age", interaction = TRUE, err_check = err_check)
  sizec2.age.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec2", factor2 = "age", interaction = TRUE, err_check = err_check)
  sizec1.age.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "sizec1", factor2 = "age", interaction = TRUE, err_check = err_check)
  dens.age.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "density", factor2 = "age", interaction = TRUE, err_check = err_check)
  
  inda2.sizeb2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "sizeb2", interaction = TRUE, err_check = err_check)
  inda2.sizec2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "sizec2", interaction = TRUE, err_check = err_check)
  inda2.sizeb1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "sizeb1", interaction = TRUE, err_check = err_check)
  inda2.sizec1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "sizec1", interaction = TRUE, err_check = err_check)
  inda1.sizeb2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "sizeb2", interaction = TRUE, err_check = err_check)
  inda1.sizec2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "sizec2", interaction = TRUE, err_check = err_check)
  inda1.sizeb1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "sizeb1", interaction = TRUE, err_check = err_check)
  inda1.sizec1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "sizec1", interaction = TRUE, err_check = err_check)
  inda2.dens.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova2", factor2 = "density", interaction = TRUE, err_check = err_check)
  inda1.dens.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcova1", factor2 = "density", interaction = TRUE, err_check = err_check)
  
  indb2.sizeb2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb2", factor2 = "sizeb2", interaction = TRUE, err_check = err_check)
  indb2.sizec2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb2", factor2 = "sizec2", interaction = TRUE, err_check = err_check)
  indb2.sizeb1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb2", factor2 = "sizeb1", interaction = TRUE, err_check = err_check)
  indb2.sizec1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb2", factor2 = "sizec1", interaction = TRUE, err_check = err_check)
  indb1.sizeb2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb1", factor2 = "sizeb2", interaction = TRUE, err_check = err_check)
  indb1.sizec2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb1", factor2 = "sizec2", interaction = TRUE, err_check = err_check)
  indb1.sizeb1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb1", factor2 = "sizeb1", interaction = TRUE, err_check = err_check)
  indb1.sizec1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb1", factor2 = "sizec1", interaction = TRUE, err_check = err_check)
  indb2.dens.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb2", factor2 = "density", interaction = TRUE, err_check = err_check)
  indb1.dens.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovb1", factor2 = "density", interaction = TRUE, err_check = err_check)
  
  indc2.sizeb2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = "sizeb2", interaction = TRUE, err_check = err_check)
  indc2.sizec2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = "sizec2", interaction = TRUE, err_check = err_check)
  indc2.sizeb1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = "sizeb1", interaction = TRUE, err_check = err_check)
  indc2.sizec1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = "sizec1", interaction = TRUE, err_check = err_check)
  indc1.sizeb2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "sizeb2", interaction = TRUE, err_check = err_check)
  indc1.sizec2.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "sizec2", interaction = TRUE, err_check = err_check)
  indc1.sizeb1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "sizeb1", interaction = TRUE, err_check = err_check)
  indc1.sizec1.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "sizec1", interaction = TRUE, err_check = err_check)
  indc2.dens.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc2", factor2 = "density", interaction = TRUE, err_check = err_check)
  indc1.dens.zi <- .handinglove(zi_model, model_code = model_type, paramnames = paramnames,
    factor1 = "indcovc1", factor2 = "density", interaction = TRUE, err_check = err_check)
  
  
  # Elements 1-92 are the original coefficients in lefko3 version 3.
  # Elements 1-46 cover the conditional. Elements 47-92 cover the zi model.
  # Elements 93-100 are 0s added to bring the vector length to 100.
  # Elements 101-283 are new with version 4.0.0.
  # Elements 101-183 are the new conditional model coefficients.
  # Elements 201-283 are the new zi model coefficients.
  # Elements 184-200 are 0s added to create the right vector length.
  if (model_type > 7) {
    coef.vec <- yintercept.coef
  } else {
    coef.vec <- c(yintercept.coef, flw1.coef, flw2.coef, size1.coef, size2.coef,
      flw1.flw2.coef, size1.size2.coef, size1.flw1.coef, size2.flw2.coef,
      size2.flw1.coef, size1.flw2.coef, age.coef, age.size1.coef, age.size2.coef,
      age.flw1.coef, age.flw2.coef, inda2.coef, indb2.coef, indc2.coef, inda1.coef,
      indb1.coef, indc1.coef, inda2.size2.coef, indb2.size2.coef, indc2.size2.coef,
      inda2.flw2.coef, indb2.flw2.coef, indc2.flw2.coef, inda1.size1.coef,
      indb1.size1.coef, indc1.size1.coef, inda1.flw1.coef, indb1.flw1.coef,
      indc1.flw1.coef, inda2.indb2.coef, inda2.indc2.coef, indb2.indc2.coef,
      inda1.indb1.coef, inda1.indc1.coef, indb1.indc1.coef, inda2.indb1.coef, 
      inda1.indb2.coef, inda2.indc1.coef, inda1.indc2.coef, indb2.indc1.coef,
      indb1.indc2.coef, 
      
      yintercept.zi, flw1.zi, flw2.zi, size1.zi, size2.zi, flw1.flw2.zi, size1.size2.zi, 
      size1.flw1.zi, size2.flw2.zi, size2.flw1.zi, size1.flw2.zi, age.zi, age.size1.zi, 
      age.size2.zi, age.flw1.zi, age.flw2.zi, inda2.zi, indb2.zi, indc2.zi, inda1.zi,
      indb1.zi, indc1.zi, inda2.size2.zi, indb2.size2.zi, indc2.size2.zi, inda2.flw2.zi,
      indb2.flw2.zi, indc2.flw2.zi, inda1.size1.zi, indb1.size1.zi, indc1.size1.zi,
      inda1.flw1.zi, indb1.flw1.zi, indc1.flw1.zi, inda2.indb2.zi, inda2.indc2.zi, 
      indb2.indc2.zi, inda1.indb1.zi, inda1.indc1.zi, indb1.indc1.zi, inda2.indb1.zi, 
      inda1.indb2.zi, inda2.indc1.zi, inda1.indc2.zi, indb2.indc1.zi, indb1.indc2.zi,
      
      0, 0, 0, 0, 0, 0, 0, 0,
      
      sizeb2.coef, sizeb1.coef, sizec2.coef, sizec1.coef, dens.coef, sizeb1.sizeb2.coef,
      sizec1.sizec2.coef, sizea1.sizeb1.coef, sizea1.sizec1.coef, sizeb1.sizec1.coef,
      sizea2.sizeb2.coef, sizea2.sizec2.coef, sizeb2.sizec2.coef, sizea1.sizeb2.coef,
      sizea1.sizec2.coef, sizeb1.sizec2.coef, sizea2.sizeb1.coef, sizea2.sizec1.coef,
      sizeb2.sizec1.coef, sizea2.dens.coef, sizeb2.dens.coef, sizec2.dens.coef,
      sizea1.dens.coef, sizeb1.dens.coef, sizec1.dens.coef, flw2.dens.coef,
      flw1.dens.coef, sizeb2.flw2.coef, sizec2.flw2.coef, 0, sizeb1.flw1.coef,
      sizeb2.flw1.coef, sizeb1.flw2.coef, sizec1.flw1.coef, sizec2.flw1.coef,
      sizec1.flw2.coef, sizeb2.age.coef, sizec2.age.coef, dens.age.coef, sizeb1.age.coef,
      sizec1.age.coef,
      
      inda2.sizeb2.coef, inda2.sizec2.coef, inda2.dens.coef, inda1.sizeb1.coef,
      inda1.sizec1.coef, inda1.sizeb2.coef, inda1.sizec2.coef, inda2.sizeb1.coef,
      inda2.sizec1.coef, inda1.dens.coef,
      
      indb2.sizeb2.coef, indb2.sizec2.coef, indb2.dens.coef, indb1.sizeb1.coef,
      indb1.sizec1.coef, indb1.sizeb2.coef, indb1.sizec2.coef, indb2.sizeb1.coef,
      indb2.sizec1.coef, indb1.dens.coef,
      
      indc2.sizeb2.coef, indc2.sizec2.coef, indc2.dens.coef, indc1.sizeb1.coef,
      indc1.sizec1.coef, indc1.sizeb2.coef, indc1.sizec2.coef, indc2.sizeb1.coef,
      indc2.sizec1.coef, indc1.dens.coef,
      
      inda2.size1.coef, indb2.size1.coef, indc2.size1.coef, inda1.size2.coef,
      indb1.size2.coef, indc1.size2.coef, inda2.flw1.coef, indb2.flw1.coef,
      indc2.flw1.coef, inda1.flw2.coef, indb1.flw2.coef, indc1.flw2.coef,
      
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      
      sizeb2.zi, sizeb1.zi, sizec2.zi, sizec1.zi, dens.zi, sizeb1.sizeb2.zi,
      sizec1.sizec2.zi, sizea1.sizeb1.zi, sizea1.sizec1.zi, sizeb1.sizec1.zi,
      sizea2.sizeb2.zi, sizea2.sizec2.zi, sizeb2.sizec2.zi, sizea1.sizeb2.zi,
      sizea1.sizec2.zi, sizeb1.sizec2.zi, sizea2.sizeb1.zi, sizea2.sizec1.zi,
      sizeb2.sizec1.zi, sizea2.dens.zi, sizeb2.dens.zi, sizec2.dens.zi,
      sizea1.dens.zi, sizeb1.dens.zi, sizec1.dens.zi, flw2.dens.zi,
      flw1.dens.zi, sizeb2.flw2.zi, sizec2.flw2.zi, 0, sizeb1.flw1.zi,
      sizeb2.flw1.zi, sizeb1.flw2.zi, sizec1.flw1.zi, sizec2.flw1.zi,
      sizec1.flw2.zi, sizeb2.age.zi, sizec2.age.zi, dens.age.zi, sizeb1.age.zi,
      sizec1.age.zi,
      
      inda2.sizeb2.zi, inda2.sizec2.zi, inda2.dens.zi, inda1.sizeb1.zi,
      inda1.sizec1.zi, inda1.sizeb2.zi, inda1.sizec2.zi, inda2.sizeb1.zi,
      inda2.sizec1.zi, inda1.dens.zi,
      
      indb2.sizeb2.zi, indb2.sizec2.zi, indb2.dens.zi, indb1.sizeb1.zi,
      indb1.sizec1.zi, indb1.sizeb2.zi, indb1.sizec2.zi, indb2.sizeb1.zi,
      indb2.sizec1.zi, indb1.dens.zi,
      
      indc2.sizeb2.zi, indc2.sizec2.zi, indc2.dens.zi, indc1.sizeb1.zi,
      indc1.sizec1.zi, indc1.sizeb2.zi, indc1.sizec2.zi, indc2.sizeb1.zi,
      indc2.sizec1.zi, indc1.dens.zi,
      
      inda2.size1.zi, indb2.size1.zi, indc2.size1.zi, inda1.size2.zi,
      indb1.size2.zi, indc1.size2.zi, inda2.flw1.zi, indb2.flw1.zi,
      indc2.flw1.zi, inda1.flw2.zi, indb1.flw2.zi, indc1.flw2.zi)
  }
  
  coef.list <- list(coefficients = coef.vec, years = year.coefs,
    patches = patch.coefs, groups2 = group2.coefs, groups1 = group1.coefs,
    indcova2s = indcova2.coefs, indcova1s = indcova1.coefs,
    indcovb2s = indcovb2.coefs, indcovb1s = indcovb1.coefs,
    indcovc2s = indcovc2.coefs, indcovc1s = indcovc1.coefs, zeroyear = year.zi,
    zeropatch = patch.zi, zerogroups2 = group2.zi, zerogroups1 = group1.zi,
    zeroindcova2s = indcova2.zi, zeroindcova1s = indcova1.zi,
    zeroindcovb2s = indcovb2.zi, zeroindcovb1s = indcovb1.zi,
    zeroindcovc2s = indcovc2.zi, zeroindcovc1s = indcovc1.zi, variances = rvars,
    family = viewedfamily, sigma = sigma, trunc = trunc0)
  
  return(coef.list)
}

