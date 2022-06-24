#' Develop Best-fit Vital Rate Estimation Models for MPM Development
#' 
#' Function \code{modelsearch()} runs exhaustive model building and selection
#' for each vital rate needed to estimate a function-based MPM or IPM. It
#' returns best-fit models for each vital rate, model table showing all models
#' tested, and model quality control data. The final output can be used as input
#' in other functions within this package.
#' 
#' @name modelsearch
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
#' status in historical model); \code{age}, all vital rates estimated with age
#' and y-intercepts only; \code{cons}, all vital rates estimated only as
#' y-intercepts. If \code{approach = "glm"} and \code{year.as.random = FALSE},
#' then year is also included as a fixed effect, and, in the case of
#' \code{full}, included in two-way interactions. Defaults to \code{size}.
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
#' reproduction (overall fecundity). May also be set to
#' \code{vitalrates = "leslie"}, which is equivalent to setting
#' \code{c("surv", "fec")} for a Leslie MPM. This choice also determines how
#' internal data subsetting for vital rate model estimation will work. Defaults
#' to \code{c("surv", "size", "fec")}.
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
#' @param matstat A vector indicating the variable names coding for maturity
#' status in occasions \emph{t}+1, \emph{t}, and \emph{t}-1. Defaults to
#' \code{c("matstatus3", "matstatus2", "matstatus1")}.
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
#' @param age Designates the name of the variable corresponding to age in time
#' \emph{t} in the vertical dataset. Defaults to \code{NA}, in which case age
#' is not included in linear models. Should only be used if building Leslie or
#' age x stage matrices.
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
#' @param accuracy A logical value indicating whether to test accuracy of
#' models. See \code{Notes} section for details on hoew accuracy is assessed.
#' Defaults to \code{TRUE}.
#' @param quiet If set to TRUE, then model building and selection will proceed
#' with most warnings and diagnostic messages silenced. Defaults to
#' \code{FALSE}.
#' 
#' @return This function yields an object of class \code{lefkoMod}, which is a
#' list in which the first 14 elements are the best-fit models for survival,
#' observation status, primary size, secondary size, tertiary size,
#' reproductive status, fecundity, juvenile survival, juvenile observation,
#' juvenile primary size, juvenile secondary size, juvenile tertiary size,
#' juvenile transition to reproduction, and juvenile transition to maturity,
#' respectively. This is followed by 14 elements corresponding to the model
#' tables for each of these vital rates, in order, followed by a data frame
#' showing the order and names of variables used in modeling, followed by a
#' single character element denoting the criterion used for model selection, and
#' ending on a data frame with quality control data:
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
#' \item{juv_maturity_model}{Best-fit model of the binomial probability of
#' becoming mature in occasion \emph{t}+1, given survival to that occasion of an
#' individual that was immature in occasion \emph{t}. Defaults to \code{1}.}
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
#' \item{juv_maturity_table}{Full dredge model table of the probability of
#' an immature individual transitioning to maturity.}
#' \item{criterion}{Character variable denoting the criterion used to determine
#' the best-fit model.}
#' \item{qc}{Data frame with five variables: 1) Name of vital rate, 2) number
#' of individuals used to model that vital rate, 3) number of individual
#' transitions used to model that vital rate, 4) parameter distribution used to
#' model the vital rats, and 5) accuracy of model, given as detailed in Notes
#' section.}
#' 
#' @section Notes:
#' Setting \code{suite = "cons"} prevents the inclusion of size and reproductive
#' status as fixed, independent factors in modeling. However, it does not
#' prevent any other terms from being included. Density, age, individual
#' covariates, individual identity, patch, and year may all be included.
#' 
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
#' When \code{modelsearch()} is called, it first builds global models for all
#' vital rates and runs them. If a global model fails, then the function
#' proceeds by dropping any two-way interactions and trying again. If this
#' fails, then function will continue to drop key terms (typically random terms
#' or individual covariates) until something rune. The last attempt if running a
#' mixed set of models is to try a glm version of the original failed model, and
#' use that as a global model if it runs properly. Finally, if all attempts
#' fail, then the function returns a \code{1}.
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
#' Accuracy of vital rate models is calculated differently depending on vital
#' rate and assumed distribution. For all vital rates assuming a binomial
#' distribution, including survival, observation status, reproductive status,
#' and juvenile version of these, accuracy is calculated as the percent of
#' predicted responses equal to actual responses. In all other models, accuracy
#' is actually the conditional R-wquared using package \code{MuMIn}'s
#' \code{\link[MuMIn]{r.squaredGLMM}()} function, estimated via the delta
#' method. When this method fails, \code{modelsearch()} calculates McFadden's
#' pseudo R-squared. If this fails, then \code{NA} is returned.
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
#' Users building vital rate models for Leslie matrices must set
#' \code{vitalrates = c("surv", "fec")} or \code{vitalrates = "leslie"} rather
#' than the default, because only survival and fecundity should be estimated in
#' these cases. Also, the \code{suite} setting can be set to either \code{age}
#' or \code{cons}, as the results will be exactly the same.
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
#'   reduce = FALSE)
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
  matstat = c("matstatus3", "matstatus2", "matstatus1"),
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
  accuracy = TRUE, quiet = FALSE) {
  
  censor1 <- censor2 <- censor3 <- surv.data <- obs.data <- size.data <- NULL
  repst.data <- fec.data <- juvsurv.data <- juvobs.data <- NULL
  juvsize.data <- juvrepst.data <- usedfec <- NULL
  patchcol <- yearcol <- extra_factors <- 0
  
  sizeb_used <- sizec_used <- density_used <- indcova_used <- FALSE
  indcovb_used <- indcovc_used <- FALSE
  
  total_vars <- length(names(data))
  
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
        extra_factors <- extra_factors + 1;
        
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
  
  # Here we will use text matching to identify the linear modeling approach and distributions
  vitalrates <- tolower(vitalrates)
  approach <- tolower(approach)
  suite <- tolower(suite)
  bestfit <- tolower(bestfit)
  sizedist <- tolower(sizedist)
  sizebdist <- tolower(sizebdist)
  sizecdist <- tolower(sizecdist)
  fecdist <- tolower(fecdist)
  
  appr_length <- length(grep("mix", approach)) + length(grep("lme", approach)) +
    length(grep("tmb", approach))
  if (appr_length > 0) {
    approach <- "mixed"
  }
  
  if(length(grep("fu", suite)) > 0) {
    suite <- "full"
  } else if(length(grep("fl", suite)) > 0) {
    suite <- "full"
  } else if(length(grep("ma", suite)) > 0) {
    suite <- "main"
  } else if(length(grep("mn", suite)) > 0) {
    suite <- "main"
  } else if(length(grep("si", suite)) > 0) {
    suite <- "size"
  } else if(length(grep("sz", suite)) > 0) {
    suite <- "size"
  } else if(length(grep("re", suite)) > 0) {
    suite <- "rep"
  } else if(length(grep("rp", suite)) > 0) {
    suite <- "rep"
  } else if(length(grep("co", suite)) > 0) {
    suite <- "cons"
  } else if(length(grep("ag", suite)) > 0) {
    if (is.na(age)) {
      stop("Age variable required for age-based and age-by-stage MPMs.", call. = FALSE)
    }
    suite <- "cons"
  }
  
  if (length(grep("&k", bestfit)) > 0) {
    bestfit <- "aicc&k"
  }
  
  if (length(grep("gaus", sizedist)) > 0) {
    sizedist <- "gaussian"
  } else if (length(grep("gam", sizedist)) > 0) {
    sizedist <- "gamma"
  } else if (length(grep("pois", sizedist)) > 0) {
    sizedist <- "poisson"
  } else if (length(grep("neg", sizedist)) > 0) {
    sizedist <- "negbin"
  }
  
  if (length(grep("gaus", sizebdist)) > 0) {
    sizebdist <- "gaussian"
  } else if (length(grep("gam", sizebdist)) > 0) {
    sizebdist <- "gamma"
  } else if (length(grep("pois", sizebdist)) > 0) {
    sizebdist <- "poisson"
  } else if (length(grep("neg", sizebdist)) > 0) {
    sizebdist <- "negbin"
  }
  
  if (length(grep("gaus", sizecdist)) > 0) {
    sizecdist <- "gaussian"
  } else if (length(grep("gam", sizecdist)) > 0) {
    sizecdist <- "gamma"
  } else if (length(grep("pois", sizecdist)) > 0) {
    sizecdist <- "poisson"
  } else if (length(grep("neg", sizecdist)) > 0) {
    sizecdist <- "negbin"
  }
  
  if (length(grep("gaus", fecdist)) > 0) {
    fecdist <- "gaussian"
  } else if (length(grep("gam", fecdist)) > 0) {
    fecdist <- "gamma"
  } else if (length(grep("pois", fecdist)) > 0) {
    fecdist <- "poisson"
  } else if (length(grep("neg", fecdist)) > 0) {
    fecdist <- "negbin"
  }
  
  if (length(grep("su", vitalrates)) > 0) {
    vitalrates[grep("su", vitalrates)] <- "surv"
  }
  if (length(grep("sr", vitalrates)) > 0) {
    vitalrates[grep("sr", vitalrates)] <- "surv"
  }
  if (length(grep("ob", vitalrates)) > 0) {
    vitalrates[grep("ob", vitalrates)] <- "obs"
  }
  if (length(grep("si", vitalrates)) > 0) {
    vitalrates[grep("si", vitalrates)] <- "size"
  }
  if (length(grep("sz", vitalrates)) > 0) {
    vitalrates[grep("sz", vitalrates)] <- "size"
  }
  if (length(grep("re", vitalrates)) > 0) {
    vitalrates[grep("re", vitalrates)] <- "repst"
  }
  if (length(grep("rp", vitalrates)) > 0) {
    vitalrates[grep("rp", vitalrates)] <- "repst"
  }
  if (length(grep("fe", vitalrates)) > 0) {
    vitalrates[grep("fe", vitalrates)] <- "fec"
  }
  if (length(grep("fc", vitalrates)) > 0) {
    vitalrates[grep("fc", vitalrates)] <- "fec"
  }
  if (length(vitalrates) == 1 & length(grep("lesl", vitalrates)) > 0) {
    vitalrates <- c("surv", "fec")
  }
  
  if (approach == "mixed" & !requireNamespace("lme4", quietly = TRUE)) {
    stop("Package lme4 is required for mixed models. Please install it.", call. = FALSE)
  } else if (approach == "mixed" & !requireNamespace("glmmTMB", quietly = TRUE)) {
    if (any(sizedist == "negbin", na.rm = TRUE) | any(sizebdist == "negbin", na.rm = TRUE) |
      any(sizecdist == "negbin", na.rm = TRUE) | any(fecdist == "negbin", na.rm = TRUE)) {
      
      stop("Package glmmTMB needed to develop mixed size or fecundity models
        with a negative binomial distribution.", call. = FALSE)
    }
  } else if (any(sizedist != "gaussian", na.rm = TRUE) | any(sizebdist != "gaussian", na.rm = TRUE) | 
      any(sizecdist != "gaussian", na.rm = TRUE)) {
    if (!requireNamespace("glmmTMB", quietly = TRUE)) {
      if (size.trunc | sizeb.trunc | sizec.trunc | size.zero | sizeb.zero | sizec.zero) {
        stop("Package glmmTMB needed to develop mixed size models with
          zero-truncated or zero-inflated distributions.", call. = FALSE)
      }
    }
  } else if (any(fecdist != "gaussian", na.rm = TRUE)) {
    if (fec.trunc & !requireNamespace("glmmTMB", quietly = TRUE)) {
      stop("Package glmmTMB needed to develop mixed fecundity models with
        zero-truncated distribution.", call. = FALSE)}
    if (fec.zero & !requireNamespace("glmmTMB", quietly = TRUE)) {
      stop("Package glmmTMB needed to develop mixed fecundity models with
        zero-inflated distribution.", call. = FALSE)}
  }
  
  if (approach == "glm") {
    if (any(sizedist != "gaussian", na.rm = TRUE) | any(sizebdist != "gaussian", na.rm = TRUE) |
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
    if (any(fecdist != "gaussian", na.rm = TRUE)) {
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
    stop("Please enter a valid primary size distribution, currently limited to
      'gaussian', 'poisson', 'negbin', and 'gamma'.", call. = FALSE)}
  if (!is.element(sizebdist, c(NA, distoptions))) {
    stop("Please enter a valid secondary size distribution, currently limited to
      'gaussian', 'poisson', 'negbin', and 'gamma'.", call. = FALSE)}
  if (!is.element(sizecdist, c(NA, distoptions))) {
    stop("Please enter a valid tertiary size distribution, currently limited to
      'gaussian', 'poisson', 'negbin', and 'gamma'.", call. = FALSE)}
  if (!is.element(fecdist, distoptions)) {
    stop("Please enter a valid fecundity distribution, currently limited to
      'gaussian', 'poisson', 'negbin', and 'gamma'.", call. = FALSE)}
  
  if (length(censor) > 3) {
    stop("Censor variables should be included either as 1 variable per row in the
      historical data frame (1 variable in the dataset), or as 1 variable for each of
      occasions t+1, t, and, if historical, t-1 (2 or 3 variables in the data frame).
      No more than 3 variables are allowed. If more than one are supplied, then they
      are assumed to be in order of occasion t+1, t, and t-1, respectively.",
      call. = FALSE)
  }
  if (length(indiv) > 1) {
    stop("Only one individual identification variable is allowed.",
      call. = FALSE)
  }
  if (length(year) > 1) {
    stop("Only one time variable is allowed, and it must refer to time t.", call. = FALSE)
  }
  if (length(patch) > 1) {
    stop("Only one patch variable is allowed.", call. = FALSE)
  }
  if (length(density) > 1) {
    stop("Only one density variable is allowed.", call. = FALSE)
  }
  if (length(age) > 1) {
    stop("Only one age variable is allowed.", call. = FALSE)
  }
  
  if (is.element("surv", vitalrates)) {
    if (length(surv) > 3 | length(surv) == 1) {
      stop("This function requires 2 (if ahistorical) or 3 (if historical)
        survival variables as input parameters.", call. = FALSE)}
    if (all(is.numeric(surv))) {
      if (any(surv < 1) | any(surv > total_vars)) {
        stop("Survival variables do not match data frame.",
          call. = FALSE)
      } else {
        surv <- names(data)[surv]
      }
    }
    if (any(!is.element(surv, names(data)))) {
      stop("Survival variables must match data frame.", call. = FALSE)
    }
  }
  
  if (is.element("obs", vitalrates)) {
    if (length(obs) > 3 | length(obs) == 1) {
      stop("This function requires 2 (if ahistorical) or 3 (if historical)
        observation variables as input parameters.",
        call. = FALSE)}
    if (all(is.numeric(obs))) {
      if (any(obs < 1) | any(obs > total_vars)) {
        stop("Observation variables do not match data frame.",
          call. = FALSE)
      } else {
       obs <- names(data)[obs]
      }
    }
    if (any(!is.element(obs, names(data)))) {
      stop("Observation status variables must match data frame.",
        call. = FALSE)
    }
  }
  
  if (is.element("size", vitalrates)) {
    if (length(size) > 3 | length(size) == 1) {
      stop("This function requires 2 (if ahistorical) or 3 (if historical) size
        variables as input parameters.", call. = FALSE)}
    if (all(is.numeric(size))) {
      if (any(size < 1) | any(size > total_vars)) {
        stop("Size variables do not match data frame.",
          call. = FALSE)
      } else {
        size <- names(data)[size]
      }
    }
    if (any(!is.element(size, names(data)))) {
      stop("Size variables must match data frame.", call. = FALSE)
    }
    if (all(!is.na(sizeb))) {
      if (is.na(sizebdist)) {
        stop("Need valid choice of distribution for secondary size.", call. = FALSE)
      }
      if (length(sizeb) > 3 | length(sizeb) == 1) {
        stop("This function requires 2 (if ahistorical) or 3 (if historical)
          secondary size variables as input parameters.",
          call. = FALSE)}
      if (all(is.numeric(sizeb))) {
        if (any(sizeb < 1) | any(sizeb > total_vars)) {
          stop("Secondary size variables do not match data frame.",
            call. = FALSE)
        } else {
          sizeb <- names(data)[sizeb]
        }
      }
      if (any(!is.element(sizeb, names(data)))) {
        stop("Secondary size variables must match data frame.",
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
          tertiary size variables as input parameters.",
          call. = FALSE)}
      if (all(is.numeric(sizec))) {
        if (any(sizec < 1) | any(sizec > total_vars)) {
          stop("Tertiary size variables do not match data frame.",
            call. = FALSE)
        } else {
          sizec <- names(data)[sizec]
        }
      }
      if (any(!is.element(sizec, names(data)))) {
        stop("Tertiary size variables must match data frame.",
          call. = FALSE)
      }
      sizec_used <- 1
    }
  }
  
  if (is.element("repst", vitalrates)) {
    if (length(repst) > 3 | length(repst) == 1) {
      stop("This function requires 2 (if ahistorical) or 3 (if historical)
        reproductive status variables as input parameters.",
        call. = FALSE)}
    if (all(is.numeric(repst))) {
      if (any(repst < 1) | any(repst > total_vars)) {
        stop("Reproductive status variables do not match data frame.",
          call. = FALSE)
      } else {
        repst <- names(data)[repst]
      }
    }
    if (any(!is.element(repst, names(data)))) {
      stop("Reproductive status variables must match data frame.",
        call. = FALSE)
    }
  }
  
  if (is.element("fec", vitalrates)) {
    if (length(fec) > 3 | length(fec) == 1) {
      stop("This function requires 2 (if ahistorical) or 3 (if historical)
        fecundity variables as input parameters.",
        call. = FALSE)}
    if (all(is.numeric(fec))) {
      if (any(fec < 1) | any(fec > total_vars)) {
        stop("Fecundity variables do not match data frame.",
          call. = FALSE)
      } else {
        fec <- names(data)[fec]
      }
    }
    if (any(!is.element(fec, names(data)))) {
      stop("Fecundity variables must match data frame.",
        call. = FALSE)
    }
  }
  
  if (fectime != 2 & fectime != 3) {
    stop("fectime must equal 2 or 3, depending on whether fecundity occurs in
      time t or t+1, respectively (the default is 2).", call. = FALSE)
  }
  
  if (all(!is.na(age))) {
    if (all(is.numeric(age))) {
      if (age > 0 & age <= total_vars) {
        agecol <- age
        age <- names(data)[agecol]
      }
    } else if (length(which(names(data) == age)) == 0) {
      stop("Variable age must either equal the name of the variable denoting age,
        or be set to NA.", call. = FALSE)
    } else {
      agecol <- which(names(data) == age)
    }
    extra_factors <- extra_factors + 1
  } else {age <- "none"}
  
  if (any(!is.na(indcova))) {
    if (length(indcova) > 3) {
      warning("Vector indcova holds the exact names of an individual covariate across
        times t+1, t, and t-1. Only the first three elements will be used.",
        call. = FALSE)
      indcova <- indcova[1:3]
      indcova_used <- TRUE
    } else if (length(indcova) == 1) {
      warning("Vector indcova requires the names of an individual covariate across
        times t+1, t, and, if historical, t-1. Only 1 variable name was supplied,
        so this individual covariate will not be used.", call. = FALSE)
      indcova <- c("none", "none", "none")
    }
    
    if (all(is.numeric(indcova), na.rm = TRUE)) {
      if (all(indcova > 0, na.rm = TRUE) & all(indcova <= total_vars)) {
        indcova2col <- indcova[2]
        if (!is.na(indcova[3])) indcova1col <- indcova[3]
        indcova_used <- TRUE
      }
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
  
  if (indcova_used & !random.indcova) {
    extra_factors <- extra_factors + 1
  }
  
  if (any(!is.na(indcovb))) {
    if (length(indcovb) > 3) {
      warning("Vector indcovb holds the exact names of an individual covariate
        across times t+1, t, and t-1. Only the first three elements will be used.",
        call. = FALSE)
      indcovb <- indcovb[1:3]
      indcovb_used <- TRUE
    } else if (length(indcovb) == 1) {
      warning("Vector indcovb requires the names of an individual covariate across
        times t+1, t, and, if historical, t-1. Only 1 variable name was supplied,
        so this individual covariate will not be used.", call. = FALSE)
      indcovb <- c("none", "none", "none")
    }
    
    if (all(is.numeric(indcovb), na.rm = TRUE)) {
      if (all(indcovb > 0, na.rm = TRUE) & all(indcovb <= total_vars)) {
        indcovb2col <- indcovb[2]
        if (!is.na(indcovb[3])) indcovb1col <- indcovb[3]
        indcovb_used <- TRUE
      }
    } else {
      if (length(which(names(data) == indcovb[2])) == 0 && indcovb[2] != "none") {
        stop("Vector indcovb must either equal either the exact names of an
          individual covariate across times t+1, t, and t-1, or be set to NA.",
          call. = FALSE)
      } else {
        indcovb2col <- which(names(data) == indcovb[2])
        indcovb_used <- TRUE
      }
      
      if (length(indcovb) == 3) {
        if (length(which(names(data) == indcovb[3])) == 0 && indcovb[3] != "none") {
          stop("Vector indcovb must either equal either the exact names of an
            individual covariate across times t+1, t, and t-1, or be set to NA.",
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
  
  if (indcovb_used & !random.indcovb) {
    extra_factors <- extra_factors + 1
  }
  
  if (any(!is.na(indcovc))) {
    if (length(indcovc) > 3) {
      warning("Vector indcovc holds the exact names of an individual covariate
        across times t+1, t, and t-1. Only the first three elements will be used.",
        call. = FALSE)
      indcovc <- indcovc[1:3]
      indcovc_used <- TRUE
    } else if (length(indcovc) == 1) {
      warning("Vector indcovc requires the names of an individual covariate across
        times t+1, t, and, if historical, t-1. Only 1 variable name was supplied,
        so this individual covariate will not be used.", call. = FALSE)
      indcovc <- c("none", "none", "none")
    } 
    
    if (all(is.numeric(indcovc), na.rm = TRUE)) {
      if (all(indcovc > 0, na.rm = TRUE) & all(indcovc <= total_vars)) {
        indcovc2col <- indcovc[2]
        if (!is.na(indcovc[3])) indcovc1col <- indcovc[3]
        indcovc_used <- TRUE
      }
    } else {
      if (length(which(names(data) == indcovc[2])) == 0 && indcovc[2] != "none") {
        stop("Vector indcovc must either equal either the exact names of an
          individual covariate across times t+1, t, and t-1, or be set to NA.",
          call. = FALSE)
      } else {
        indcovc2col <- which(names(data) == indcovc[2])
        indcovc_used <- TRUE
      }
      
      if (length(indcovc) == 3) {
        if (length(which(names(data) == indcovc[3])) == 0 && indcovc[3] != "none") {
          stop("Vector indcovc must either equal either the exact names of an
            individual covariate across times t+1, t, and t-1, or be set to NA.",
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
  
  if (indcovc_used & !random.indcovc) {
    extra_factors <- extra_factors + 1
  }
  
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
      if (any(indiv < 1) | any(indiv > total_vars)) {
        stop("Unable to interpret indiv variable.", call. = FALSE)
      } else {
        indivcol <- indiv
        indiv <- names(data)[indivcol]
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
      if (any(patch < 1) | any(patch > total_vars)) {
        stop("Unable to interpret patch variable.", call. = FALSE)
      } else {
        patchcol <- patch  # Used to be names(data)[patch]
        patch <- names(data)[patchcol] 
      }
    } else {
      stop("Unable to interpret patch variable.", call. = FALSE)
    }
  } else {patch <- "none"}
  
  if (patchcol > 0 & !patch.as.random) {
    extra_factors <- extra_factors + 1
  }
  
  if (!is.na(year)) {
    if (!is.numeric(year) & length(which(names(data) == year)) == 0) {
      stop("Variable year must either equal the exact name of the variable denoting
        occasion t in the dataset, or be set to NA.", call. = FALSE)
    } else if (is.character(year)) {
      yearcol <- which(names(data) == year)
    } else if (is.numeric(year)) {
      if (any(year < 1) | any(year > total_vars)) {
        stop("Unable to interpret year variable.", call. = FALSE)
      } else {
        yearcol <- year  # Used to be names(data)[year]
        year <- names(data)[yearcol] 
      }
    } else {
      stop("Unable to interpret year variable.", call. = FALSE)
    }
  } else {year <- "none"}
  
  if (yearcol > 0 & !year.as.random) {
    extra_factors <- extra_factors + 1
  }
  
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
      if (any(density < 1) | any(density > total_vars)) {
        stop("Unable to interpret density variable.", call. = FALSE)
      } else {
        density <- names(data)[density]
        density_used <- TRUE
      }
    } else {
      stop("Unable to interpret density variable.", call. = FALSE)
    }
  } else {density <- "none"}
  
  if (density_used) {
    extra_factors <- extra_factors + 1
  }
  
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
    if (length(matstat) > 3 | length(matstat) == 1) {
      stop("This function requires 2 (if ahistorical) or 3 (if historical) maturity status
        variables as input parameters if juvenile parameters are to be estimated.",
        call. = FALSE)}
    if (all(is.numeric(matstat))) {
      if (any(matstat < 1) | any(matstat > total_vars)) {
        stop("Maturity status variables do not match data frame.",
          call. = FALSE)
      } else {
        matstat <- names(data)[matstat]
      }
    }
    if (any(!is.element(matstat, names(data)))) {
      stop("Maturity status variables must match data frame.",
        call. = FALSE)
    }
  }
  
  if (!is.logical(juvsize)) {
    stop("Option juvsize must be set to either TRUE or FALSE.", call. = FALSE)
  }
  
  #Now we check whether the best-fit criterion is appropriate
  if (!is.element(bestfit, c("aicc", "aicc&k"))) {
    stop("bestfit must equal either 'AICc' or 'AICc&k', with the latter as the default.",
      call. = FALSE)
  }
  #This variable will be used once dredging is done to determine best-fit models
  used.criterion <- "AICc"
  
  if (approach != "mixed") {
    if (random.indcova | random.indcovb | random.indcovc) {
      warning("Random covariates can only be included in mixed models. Setting
        random.indcova, random.indcovb, and random.indcovc to FALSE.",
        call. = FALSE)
      random.indcova <- FALSE
      random.indcovb <- FALSE
      random.indcovc <- FALSE
    }
    
    if (patch.as.random) {
      warning("Patch can only be random in mixed models. Setting patch.as.random to FALSE.",
        call. = FALSE)
      patch.as.random = FALSE
    }
    
    if (year.as.random) {
      warning("Year can only be random in mixed models. Setting year.as.random to FALSE.",
        call. = FALSE)
      year.as.random = FALSE
    }
  }
  
  #The next section creates the text lines needed for the main model calls, based on function input
  formulae <- .stovokor(surv, obs, size, sizeb, sizec, repst, fec, matstat,
    vitalrates, historical, suite, approach, is.na(juvestimate), age, indcova,
    indcovb, indcovc, indiv, patch, year, patch.as.random, year.as.random,
    random.indcova, random.indcovb, random.indcovc, fectime, juvsize,
    sizeb_used, sizec_used, test.group, density, density_used, indcova_used,
    indcovb_used, indcovc_used)
  
  if (suite == "full") {
    alt.formulae <- .stovokor(surv, obs, size, sizeb, sizec, repst, fec,
      matstat, vitalrates, historical, suite = "main", approach,
      is.na(juvestimate), age, indcova, indcovb, indcovc, indiv, patch, year,
      patch.as.random, year.as.random, random.indcova, random.indcovb,
      random.indcovc, fectime, juvsize, sizeb_used, sizec_used, test.group,
      density, density_used, indcova_used, indcovb_used, indcovc_used)
  } else {
    alt.formulae <- list(full.surv.model = NA, full.obs.model = NA, full.size.model = NA,
      full.sizeb.model = NA, full.sizec.model = NA, full.repst.model = NA, full.fec.model = NA,
      juv.surv.model = NA, juv.obs.model = NA, juv.size.model = NA, just.sizeb.model = NA,
      juve.sizec.model = NA, juve.repst.model = NA, juv.matst.model = NA)
  }
  
  alt.nocovs.formulae <- .stovokor(surv, obs, size, sizeb, sizec, repst, fec,
    matstat, vitalrates, historical, suite, approach, is.na(juvestimate), age,
    indcova, indcovb, indcovc, indiv, patch, year, patch.as.random,
    year.as.random, FALSE, FALSE, FALSE, fectime, juvsize,
    sizeb_used, sizec_used, test.group, density, density_used, FALSE,
    FALSE, FALSE)
  
  if (approach == "mixed") {
    alt.glm.formulae <- .stovokor(surv, obs, size, sizeb, sizec, repst, fec,
      matstat, vitalrates, historical, suite, "glm", is.na(juvestimate), age,
      indcova, indcovb, indcovc, indiv, patch, year, FALSE, FALSE, FALSE, FALSE,
      FALSE, fectime, juvsize, sizeb_used, sizec_used, test.group, density,
      density_used, indcova_used, indcovb_used, indcovc_used)
  } else {
    alt.glm.formulae <- formulae
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
    if (fectime == 2) {
      fec.data <- subset(surv.data, surv.data[, which(names(repst.data) == repst[2])] == 1)
      fec.data <- fec.data[which(!is.na(fec.data[, which(names(fec.data) == fec[2])])),]
    } else if (fectime == 3) {
      fec.data <- subset(surv.data, surv.data[, which(names(repst.data) == repst[1])] == 1)
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
    if (sizedist == "poisson" | sizedist == "negbin") {
      if (any(size.data[, which(names(size.data) == size[1])] != 
          round(size.data[, which(names(size.data) == size[1])]))) {
        stop("Size variables must be composed only of integers for the Poisson or
          negative binomial distributions to be used.", call. = FALSE)
      }
      
      if (!is.na(juvestimate)) {
        if (any(juvsize.data[, which(names(juvsize.data) == size[1])] != 
            round(juvsize.data[, which(names(juvsize.data) == size[1])]))) {
          stop("Size variables must be composed only of integers for the Poisson or
            negative binomial distributions to be used.", call. = FALSE)
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
    if (sizebdist == "poisson" | sizebdist == "negbin") {
      if (any(sizeb.data[, which(names(sizeb.data) == sizeb[1])] != 
          round(sizeb.data[, which(names(sizeb.data) == sizeb[1])]))) {
        stop("Secondary size variables must be composed only of integers for the Poisson
          or negative binomial distributions to be used.", call. = FALSE)
      }
      
      if (!is.na(juvestimate)) {
        if (any(juvsizeb.data[, which(names(juvsizeb.data) == sizeb[1])] != 
            round(juvsizeb.data[, which(names(juvsizeb.data) == sizeb[1])]))) {
          stop("Secondary size variables must be composed only of integers for the Poisson
            or negative binomial distributions to be used.", call. = FALSE)
        }
      }
    }
  } else if (sizebdist == "gamma" & sizeb_used) {
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
    if (sizecdist == "poisson" | sizecdist == "negbin") {
      if (any(sizec.data[, which(names(sizec.data) == sizec[1])] != 
          round(sizec.data[, which(names(sizec.data) == sizec[1])]))) {
        stop("Tertiary size variables must be composed only of integers for the Poisson
          or negative binomial distributions to be used.", call. = FALSE)
      }
      
      if (!is.na(juvestimate)) {
        if (any(juvsizec.data[, which(names(juvsizec.data) == sizec[1])] != 
            round(juvsizec.data[, which(names(juvsizec.data) == sizec[1])]))) {
          stop("Tertiary size variables must be composed only of integers for the Poisson
            or negative binomial distributions to be used.", call. = FALSE)
        }
      }
    }
  } else if (sizecdist == "gamma" & sizec_used) {
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
  
  if (is.element(fecdist, c("poisson", "negbin")) & !is.numeric(formulae$full.fec.model)) {
    if (any(fec.data[, usedfec] != round(fec.data[, usedfec]))) {
      stop("Fecundity variables must be composed only of integers for the Poisson
        or negative binomial distributions to be used.", call. = FALSE)
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
  if (is.numeric(formulae$juv.matst.model)) {juv.matst.global.model <- formulae$juv.matst.model}
  if (is.character(formulae$juv.matst.model)) {
    if (formulae$juv.matst.model == "none") {
      juv.matst.global.model <- 1
    }
  }
  
  surv.table <- obs.table <- size.table <- sizeb.table <- sizec.table <- repst.table <- NA
  juvsurv.table <- juvobs.table <- juvsize.table <- juvsizeb.table <- juvsizec.table <- NA
  fec.table <- juvrepst.table <- juvmatst.table <- NA
  
  surv.bf <- obs.bf <- size.bf <- sizeb.bf <- sizec.bf <- repst.bf <- fec.bf <- NA
  juvsurv.bf <- juvobs.bf <- juvsize.bf <- juvsizeb.bf <- juvsizec.bf <- juvrepst.bf <- juvmatst.bf <- NA
  
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
    
    chosen_var <- which(names(surv.data) == surv[1])
    if (is.element(0, surv.data[, chosen_var]) & is.element(1, surv.data[, chosen_var])) {
      
      surv.global.list <- .headmaster_ritual(vrate = 1, approach = approach,
        dist = "binom", zero = FALSE, truncz = FALSE, quiet = quiet,
        usedformula = formulae$full.surv.model, subdata = surv.data,
        vind = surv.ind, vtrans = surv.trans, suite = suite,
        global.only = global.only, criterion = used.criterion,
        bestfit = bestfit, correction.patch, correction.year, correction.indiv,
        alt_formula = alt.formulae$full.surv.model,
        alt_nocovsformula = alt.nocovs.formulae$full.surv.model,
        alt_glmformula = alt.glm.formulae$full.surv.model, extra_fac = extra_factors)
      
      surv.global.model <- surv.global.list$model
      surv.ind <- surv.global.list$ind
      surv.trans <- surv.global.list$trans
      surv.table <- surv.global.list$table
      surv.bf <- surv.global.list$bf.model
      
      surv.accuracy <- .accu_predict(bestfitmodel = surv.bf,
        subdata = surv.data, param = surv[1], quiet = quiet, check = accuracy)
      
    } else if (!is.element(0, surv.data[, chosen_var])) {
      if (!quiet) {message("\nSurvival response is constant so will not model it.")}
      formulae$full.surv.model <- 1
      surv.global.model <- 1
      surv.bf <- 1
      surv.ind <- 0
      surv.trans <- 0
      surv.accuracy <- 1
    } else {
      if (!quiet) {message("\nSurvival response is constant so will not model it.")}
      formulae$full.surv.model <- 1
      surv.global.model <- 0
      surv.bf <- 0
      surv.ind <- 0
      surv.trans <- 0
      surv.accuracy <- 1
    }
  } else {
    surv.global.model <- 1
    surv.bf <- 1
    surv.ind <- 0
    surv.trans <- 0
    surv.accuracy <- NA
  }
  
  chosen_var <- which(names(obs.data) == obs[1])
  if (!is.numeric(formulae$full.obs.model) & nchar(formulae$full.obs.model) > 1) {
    if (is.element(0, obs.data[, chosen_var]) & is.element(1, obs.data[, chosen_var])) {
        
      obs.global.list <- .headmaster_ritual(vrate = 2, approach = approach, 
        dist = "binom", zero = FALSE, truncz = FALSE, quiet = quiet, 
        usedformula = formulae$full.obs.model, subdata = obs.data,
        vind = obs.ind, vtrans = obs.trans, suite = suite,
        global.only = global.only, criterion = used.criterion,
        bestfit = bestfit, correction.patch, correction.year, correction.indiv,
        alt_formula = alt.formulae$full.obs.model,
        alt_nocovsformula = alt.nocovs.formulae$full.obs.model,
        alt_glmformula = alt.glm.formulae$full.obs.model, extra_fac = extra_factors)
      
      obs.global.model <- obs.global.list$model
      obs.ind <- obs.global.list$ind
      obs.trans <- obs.global.list$trans
      obs.table <- obs.global.list$table
      obs.bf <- obs.global.list$bf.model
      
      obs.accuracy <- .accu_predict(bestfitmodel = obs.bf, subdata = obs.data,
        param = obs[1], quiet = quiet, check = accuracy)
      
    } else if (!is.element(0, obs.data[, chosen_var])) {
      if (!quiet) {message("\nObservation response is constant so will not model it.")}
      formulae$full.obs.model <- 1
      obs.global.model <- 1
      obs.bf <- 1
      obs.ind <- 0
      obs.trans <- 0
      obs.accuracy <- 1
    } else {
      if (!quiet) {message("\nObservation response is constant so will not model it.")}
      formulae$full.obs.model <- 1
      obs.global.model <- 0
      obs.bf <- 0
      obs.ind <- 0
      obs.trans <- 0
      obs.accuracy <- 1
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
      alt_formula = alt.formulae$full.size.model,
      alt_nocovsformula = alt.nocovs.formulae$full.size.model,
      alt_glmformula = alt.glm.formulae$full.size.model, extra_fac = extra_factors,
      null_model = TRUE)
    
    size.global.model <- size.global.list$model
    size.ind <- size.global.list$ind
    size.trans <- size.global.list$trans
    size.table <- size.global.list$table
    size.bf <- size.global.list$bf.model
    size.null <- size.global.list$null.model
    
    size.accuracy <- .accu_predict(bestfitmodel = size.bf,
      subdata = size.data, param = size[1], style = 2,
      nullmodel = size.null, quiet = quiet, check = accuracy)
    
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
      alt_formula = alt.formulae$full.sizeb.model,
      alt_nocovsformula = alt.nocovs.formulae$full.sizeb.model,
      alt_glmformula = alt.glm.formulae$full.sizeb.model, extra_fac = extra_factors,
      null_model = TRUE)
    
    sizeb.global.model <- sizeb.global.list$model
    sizeb.ind <- sizeb.global.list$ind
    sizeb.trans <- sizeb.global.list$trans
    sizeb.table <- sizeb.global.list$table
    sizeb.bf <- sizeb.global.list$bf.model
    sizeb.null <- sizeb.global.list$null.model
    
    sizeb.accuracy <- .accu_predict(bestfitmodel = sizeb.bf,
      subdata = sizeb.data, param = sizeb[1], style = 2,
      nullmodel = sizeb.null, quiet = quiet, check = accuracy)
    
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
      alt_formula = alt.formulae$full.sizec.model,
      alt_nocovsformula = alt.nocovs.formulae$full.sizec.model,
      alt_glmformula = alt.glm.formulae$full.sizec.model, extra_fac = extra_factors,
      null_model = TRUE)
    
    sizec.global.model <- sizec.global.list$model
    sizec.ind <- sizec.global.list$ind
    sizec.trans <- sizec.global.list$trans
    sizec.table <- sizec.global.list$table
    sizec.bf <- sizec.global.list$bf.model
    sizec.null <- sizec.global.list$null.model
    
    sizec.accuracy <- .accu_predict(bestfitmodel = sizec.bf,
      subdata = sizec.data, param = sizec[1], style = 2,
      nullmodel = sizec.null, quiet = quiet, check = accuracy)
    
  } else {
    sizec.global.model <- 1
    sizec.bf <- 1
    sizec.ind <- 0
    sizec.trans <- 0
    sizec.accuracy <- NA
  }
  
  if (!is.numeric(formulae$full.repst.model) & nchar(formulae$full.repst.model) > 1) {
    chosen_var <- which(names(repst.data) == repst[1])
    if (is.element(0, repst.data[, chosen_var]) & is.element(1, repst.data[, chosen_var])) {
      repst.global.list <- .headmaster_ritual(vrate = 4, approach = approach, 
        dist = "binom", zero = FALSE, truncz = FALSE, quiet = quiet, 
        usedformula = formulae$full.repst.model, subdata = repst.data,
        vind = repst.ind, vtrans = repst.trans, suite = suite,
        global.only = global.only, criterion = used.criterion,
        bestfit = bestfit, correction.patch, correction.year, correction.indiv,
        alt_formula = alt.formulae$full.repst.model,
        alt_nocovsformula = alt.nocovs.formulae$full.repst.model,
        alt_glmformula = alt.glm.formulae$full.repst.model, extra_fac = extra_factors)
      
      repst.global.model <- repst.global.list$model
      repst.ind <- repst.global.list$ind
      repst.trans <- repst.global.list$trans
      repst.table <- repst.global.list$table
      repst.bf <- repst.global.list$bf.model
      
      repst.accuracy <- .accu_predict(bestfitmodel = repst.bf,
        subdata = repst.data, param = repst[1], quiet = quiet, check = accuracy)
      
    } else if (!is.element(0, repst.data[, chosen_var])) {
      if (!quiet) {message("\nReproductive status response is constant so will not model it.")}
      formulae$full.repst.model <- 1
      repst.global.model <- 1
      repst.bf <- 1
      repst.ind <- 0
      repst.trans <- 0
      repst.accuracy <- 1
    } else {
      if (!quiet) {message("\nReproductive status response is constant so will not model it.")}
      formulae$full.repst.model <- 1
      repst.global.model <- 0
      repst.bf <- 0
      repst.ind <- 0
      repst.trans <- 0
      repst.accuracy <- 1
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
      alt_formula = alt.formulae$full.fec.model,
      alt_nocovsformula = alt.nocovs.formulae$full.fec.model,
      alt_glmformula = alt.glm.formulae$full.fec.model, extra_fac = extra_factors,
      null_model = TRUE)
    
    fec.global.model <- fec.global.list$model
    fec.ind <- fec.global.list$ind
    fec.trans <- fec.global.list$trans
    fec.table <- fec.global.list$table
    fec.bf <- fec.global.list$bf.model
    fec.null <- fec.global.list$null.model
    
    fec.accuracy <- .accu_predict(bestfitmodel = fec.bf, subdata = fec.data,
      param = names(fec.data)[usedfec], style = 2, nullmodel = fec.null,
      quiet = quiet, check = accuracy)
    
  } else {
    fec.global.model <- 1
    fec.bf <- 1
    fec.ind <- 0
    fec.trans <- 0
    fec.accuracy <- NA
  }
  
  if (!is.numeric(formulae$juv.surv.model) & nchar(formulae$juv.surv.model) > 1) {
    chosen_var <- which(names(juvsurv.data) == surv[1])
    if (is.element(0, juvsurv.data[, chosen_var]) & is.element(1, juvsurv.data[, chosen_var])) {
      
      juv.surv.global.list <- .headmaster_ritual(vrate = 6, approach = approach,
        dist = "binom", zero = FALSE, truncz = FALSE, quiet = quiet,
        usedformula = formulae$juv.surv.model, subdata = juvsurv.data,
        vind = juvsurv.ind, vtrans = juvsurv.trans, suite = suite,
        global.only = global.only, criterion = used.criterion,
        bestfit = bestfit, correction.patch, correction.year, correction.indiv,
        alt_formula = alt.formulae$juv.surv.model,
        alt_nocovsformula = alt.nocovs.formulae$juv.surv.model,
        alt_glmformula = alt.glm.formulae$juv.surv.model, extra_fac = extra_factors)
      
      juv.surv.global.model <- juv.surv.global.list$model
      juvsurv.ind <- juv.surv.global.list$ind
      juvsurv.trans <- juv.surv.global.list$trans
      juvsurv.table <- juv.surv.global.list$table
      juvsurv.bf <- juv.surv.global.list$bf.model
      
      juvsurv.accuracy <- .accu_predict(bestfitmodel = juvsurv.bf,
        subdata = juvsurv.data, param = surv[1], quiet = quiet,
        check = accuracy)
        
    } else if (!is.element(0, juvsurv.data[, chosen_var])) {
      if (!quiet) {message("\nJuvenile survival response is constant so will not model it.")}
      formulae$juv.surv.model <- 1
      juv.surv.global.model <- 1
      juvsurv.bf <- 1
      juvsurv.ind <- 0
      juvsurv.trans <- 0
      juvsurv.accuracy <- 1
    } else {
      if (!quiet) {message("\nJuvenile survival response is constant so will not model it.")}
      formulae$juv.surv.model <- 1
      juv.surv.global.model <- 0
      juvsurv.bf <- 0
      juvsurv.ind <- 0
      juvsurv.trans <- 0
      juvsurv.accuracy <- 1
    }
    
    chosen_var <- which(names(juvobs.data) == matstat[1])
    if (is.element(0, juvobs.data[, chosen_var]) & is.element(1, juvobs.data[, chosen_var])) {
      
      juv.matst.global.list <- .headmaster_ritual(vrate = 14, approach = approach,
        dist = "binom", zero = FALSE, truncz = FALSE, quiet = quiet,
        usedformula = formulae$juv.matst.model, subdata = juvobs.data,
        vind = juvobs.ind, vtrans = juvobs.trans, suite = suite,
        global.only = global.only, criterion = used.criterion,
        bestfit = bestfit, correction.patch, correction.year, correction.indiv,
        alt_formula = alt.formulae$juv.matst.model,
        alt_nocovsformula = alt.nocovs.formulae$juv.matst.model,
        alt_glmformula = alt.glm.formulae$juv.matst.model, extra_fac = extra_factors)
      
      juv.matst.global.model <- juv.matst.global.list$model
      juvmatst.ind <- juv.matst.global.list$ind
      juvmatst.trans <- juv.matst.global.list$trans
      juvmatst.table <- juv.matst.global.list$table
      juvmatst.bf <- juv.matst.global.list$bf.model
      
      juvmatst.accuracy <- .accu_predict(bestfitmodel = juvmatst.bf,
        subdata = juvobs.data, param = matstat[1], quiet = quiet,
        check = accuracy)
        
    } else if (!is.element(0, juvobs.data[, chosen_var])) {
      if (!quiet) {message("\nJuvenile maturity status response is constant so will not model it.")}
      formulae$juv.matst.model <- 1
      juv.matst.global.model <- 1
      juvmatst.bf <- 1
      juvmatst.ind <- 0
      juvmatst.trans <- 0
      juvmatst.accuracy <- 1
    } else {
      if (!quiet) {message("\nJuvenile maturity status response is constant so will not model it.")}
      formulae$juv.matst.model <- 1
      juv.matst.global.model <- 0
      juvmatst.bf <- 0
      juvmatst.ind <- 0
      juvmatst.trans <- 0
      juvmatst.accuracy <- 1
    }
  } else {
    juv.surv.global.model <- 1
    juvsurv.bf <- 1
    juvsurv.ind <- 0
    juvsurv.trans <- 0
    juvsurv.accuracy <- NA
    
    juv.matst.global.model <- 1
    juvmatst.bf <- 1
    juvmatst.ind <- 0
    juvmatst.trans <- 0
    juvmatst.accuracy <- NA
  }
  
  if (!is.numeric(formulae$juv.obs.model) & nchar(formulae$juv.obs.model) > 1) {
    chosen_var <- which(names(juvobs.data) == obs[1])
    if (is.element(0, juvobs.data[, chosen_var]) & is.element(1, juvobs.data[, chosen_var])) {
      
      juv.obs.global.list <- .headmaster_ritual(vrate = 7, approach = approach,
        dist = "binom", zero = FALSE, truncz = FALSE, quiet = quiet,
        usedformula = formulae$juv.obs.model, subdata = juvobs.data,
        vind = juvobs.ind, vtrans = juvobs.trans, suite = suite,
        global.only = global.only, criterion = used.criterion,
        bestfit = bestfit, correction.patch, correction.year, correction.indiv,
        alt_formula = alt.formulae$juv.obs.model,
        alt_nocovsformula = alt.nocovs.formulae$juv.obs.model,
        alt_glmformula = alt.glm.formulae$juv.obs.model, extra_fac = extra_factors)
      
      juv.obs.global.model <- juv.obs.global.list$model
      juvobs.ind <- juv.obs.global.list$ind
      juvobs.trans <- juv.obs.global.list$trans
      juvobs.table <- juv.obs.global.list$table
      juvobs.bf <- juv.obs.global.list$bf.model
      
      juvobs.accuracy <- .accu_predict(bestfitmodel = juvobs.bf,
        subdata = juvobs.data, param = obs[1], quiet = quiet, check = accuracy)
        
    } else if (!is.element(0, juvobs.data[, chosen_var])) {
      if (!quiet) {message("\nJuvenile observation response is constant so will not model it.")}
      formulae$juv.obs.model <- 1
      juv.obs.global.model <- 1
      juvobs.bf <- 1
      juvobs.ind <- 0
      juvobs.trans <- 0
      juvobs.accuracy <- 1
    } else {
      if (!quiet) {message("\nJuvenile observation response is constant so will not model it.")}
      formulae$juv.obs.model <- 1
      juv.obs.global.model <- 0
      juvobs.bf <- 0
      juvobs.ind <- 0
      juvobs.trans <- 0
      juvobs.accuracy <- 1
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
      alt_formula = alt.formulae$juv.size.model,
      alt_nocovsformula = alt.nocovs.formulae$juv.size.model,
      alt_glmformula = alt.glm.formulae$juv.size.model, extra_fac = extra_factors,
      null_model = TRUE)
    
    juv.size.global.model <- juv.size.global.list$model
    juvsize.ind <- juv.size.global.list$ind
    juvsize.trans <- juv.size.global.list$trans
    juvsize.table <- juv.size.global.list$table
    juvsize.bf <- juv.size.global.list$bf.model
    juvsize.null <- juv.size.global.list$null.model
    
    juvsize.accuracy <- .accu_predict(bestfitmodel = juvsize.bf,
      subdata = juvsize.data, param = size[1], style = 2,
      nullmodel = juvsize.null, quiet = quiet, check = accuracy)
    
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
      alt_formula = alt.formulae$juv.sizeb.model,
      alt_nocovsformula = alt.nocovs.formulae$juv.sizeb.model,
      alt_glmformula = alt.glm.formulae$juv.sizeb.model, extra_fac = extra_factors,
      null_model = TRUE)
    
    juv.sizeb.global.model <- juv.sizeb.global.list$model
    juvsizeb.ind <- juv.sizeb.global.list$ind
    juvsizeb.trans <- juv.sizeb.global.list$trans
    juvsizeb.table <- juv.sizeb.global.list$table
    juvsizeb.bf <- juv.sizeb.global.list$bf.model
    juvsizeb.null <- juv.sizeb.global.list$null.model
    
    juvsizeb.accuracy <- .accu_predict(bestfitmodel = juvsizeb.bf,
      subdata = juvsizeb.data, param = sizeb[1], style = 2,
      nullmodel = juvsizeb.null, quiet = quiet, check = accuracy)
    
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
      alt_formula = alt.formulae$juv.sizec.model,
      alt_nocovsformula = alt.nocovs.formulae$juv.sizec.model,
      alt_glmformula = alt.glm.formulae$juv.sizec.model, extra_fac = extra_factors,
      null_model = TRUE)
    
    juv.sizec.global.model <- juv.sizec.global.list$model
    juvsizec.ind <- juv.sizec.global.list$ind
    juvsizec.trans <- juv.sizec.global.list$trans
    juvsizec.table <- juv.sizec.global.list$table
    juvsizec.bf <- juv.sizec.global.list$bf.model
    juvsizec.null <- juv.sizec.global.list$null.model
    
    juvsizec.accuracy <- .accu_predict(bestfitmodel = juvsizec.bf,
      subdata = juvsizec.data, param = sizec[1], style = 2,
      nullmodel = juvsizec.null, quiet = quiet, check = accuracy)
    
  } else {
    juv.sizec.global.model <- 1
    juvsizec.bf <- 1
    juvsizec.ind <- 0
    juvsizec.trans <- 0
    juvsizec.accuracy <- NA
  }
  
  if (!is.numeric(formulae$juv.repst.model) & nchar(formulae$juv.repst.model) > 1) {
    chosen_var <- which(names(juvrepst.data) == repst[1])
    if (is.element(0, juvrepst.data[, chosen_var]) & is.element(1, juvrepst.data[, chosen_var])) {
      juv.repst.global.list <- .headmaster_ritual(vrate = 9, approach = approach, 
        dist = "binom", zero = FALSE, truncz = FALSE, quiet = quiet, 
        usedformula = formulae$juv.repst.model, subdata = juvrepst.data,
        vind = juvrepst.ind, vtrans = juvrepst.trans, suite = suite,
        global.only = global.only, criterion = used.criterion,
        bestfit = bestfit, correction.patch, correction.year, correction.indiv,
        alt_formula = alt.formulae$juv.repst.model,
        alt_nocovsformula = alt.nocovs.formulae$juv.repst.model,
        alt_glmformula = alt.glm.formulae$juv.repst.model, extra_fac = extra_factors)
      
      juv.repst.global.model <- juv.repst.global.list$model
      juvrepst.ind <- juv.repst.global.list$ind
      juvrepst.trans <- juv.repst.global.list$trans
      juvrepst.table <- juv.repst.global.list$table
      juvrepst.bf <- juv.repst.global.list$bf.model
      
      juvrepst.accuracy <- .accu_predict(bestfitmodel = juvrepst.bf,
        subdata = juvrepst.data, param = repst[1], quiet = quiet,
        check = accuracy)
      
    } else if (!is.element(0, juvrepst.data[, chosen_var])) {
      if (!quiet) {message("\nJuvenile reproductive status response is constant so will not model it.")}
      formulae$juv.repst.model <- 1
      juv.repst.global.model <- 1
      juvrepst.bf <- 1
      juvrepst.ind <- 0
      juvrepst.trans <- 0
      juvrepst.accuracy <- 1
    } else {
      if (!quiet) {message("\nJuvenile reproductive status response is constant so will not model it.")}
      formulae$juv.repst.model <- 1
      juv.repst.global.model <- 0
      juvrepst.bf <- 0
      juvrepst.ind <- 0
      juvrepst.trans <- 0
      juvrepst.accuracy <- 1
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
      "juvenile_sizec", "juvenile_reproduction", "juvenile_maturity"),
    c(surv.ind, obs.ind, size.ind, sizeb.ind, sizec.ind, repst.ind, fec.ind, juvsurv.ind,
      juvobs.ind, juvsize.ind, juvsizeb.ind, juvsizec.ind, juvrepst.ind, juvmatst.ind),
    c(surv.trans, obs.trans, size.trans, sizeb.trans, sizec.trans, repst.trans, fec.trans,
      juvsurv.trans, juvobs.trans, juvsize.trans, juvsizeb.trans, juvsizec.trans,
      juvrepst.trans, juvmatst.trans),
    c("binomial", "binomial", sizedist, sizebdist, sizecdist, "binomial", fecdist,
      "binomial", "binomial", sizedist, sizebdist, sizecdist, "binomial", "binomial"),
    c(surv.accuracy, obs.accuracy, size.accuracy, sizeb.accuracy, sizec.accuracy,
      repst.accuracy, fec.accuracy, juvsurv.accuracy, juvobs.accuracy, juvsize.accuracy, 
      juvsizeb.accuracy, juvsizec.accuracy, juvrepst.accuracy, juvmatst.accuracy)
  )
  names(qcoutput) <- c("vital_rate", "individuals", "transitions",
    "distribution", "accuracy")
  
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
    juvmatst.table <- NA
  }
  
  if (global.only) {bestfit <- "global model only"}
  
  #Now we develop the final output, creating a new S3 class to do it
  full.output <- list(survival_model = surv.bf, observation_model = obs.bf,
    size_model = size.bf, sizeb_model = sizeb.bf, sizec_model = sizec.bf,
    repstatus_model = repst.bf, fecundity_model = fec.bf, 
    juv_survival_model = juvsurv.bf, juv_observation_model = juvobs.bf, 
    juv_size_model = juvsize.bf, juv_sizeb_model = juvsizeb.bf, 
    juv_sizec_model = juvsizec.bf, juv_reproduction_model = juvrepst.bf,
    juv_maturity_model = juvmatst.bf, survival_table = surv.table,
    observation_table = obs.table, size_table = size.table,
    sizeb_table = sizeb.table, sizec_table = sizec.table,
    repstatus_table = repst.table, fecundity_table = fec.table,
    juv_survival_table = juvsurv.table, juv_observation_table = juvobs.table,
    juv_size_table = juvsize.table, juv_sizeb_table = juvsizeb.table,
    juv_sizec_table = juvsizec.table, juv_reproduction_table = juvrepst.table,
    juv_maturity_table = juvmatst.table, paramnames = formulae$paramnames, 
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
#' @name .headmaster_ritual
#' 
#' @param vrate An integer value indicating which model to build, as follows:
#' \code{1}: Adult survival probability, \code{2}: Adult observation
#' probability, \code{3}: Adult primary size, \code{4}: Adult reproduction
#' probability, \code{5}: Adult fecundity rate, \code{6}: Juvenile survival
#' probability, \code{7}: Juvenile observation probability, \code{9}: Juvenile
#' reproduction probability, \code{10}: Adult secondary size, \code{11}: Adult
#' tertiary size, \code{12}: Juvenile secondary size, \code{13}: Juvenile
#' tertiary size, and \code{14}: Juvenile maturity status.
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
#' @param alt_formula A string indicating an alternative formula to work with,
#' provided that \code{suite = "full"}. In other cases, \code{NA}.
#' @param alt_nocovsformula An alternative formula without individual
#' covariates.
#' @param alt_glmformula An alternative formula assuming all fixed factors and
#' \code{approach = "glm"}.
#' @param extra_fac An integer giving the number of extra fixed factors to
#' include in the model. Primarily used to determine whether to dredge models
#' if \code{suite = "cons"}.
#' @param null_model A logical value indicating whether to extract the null
#' (y-intercept only) model.
#'
#' @return Function \code{ms_binom()} outputs a list containing a global model,
#' the number of individuals and transitions used in modeling, the best-fit
#' model, the null model (y-intercept only), and a dredge model table.
#' 
#' @keywords internal
#' @noRd
.headmaster_ritual <- function(vrate, approach = "mixed", dist = NA,
  zero = FALSE, truncz = FALSE, quiet = FALSE, usedformula, subdata, vind,
  vtrans, suite, global.only = FALSE, criterion = "AICc", bestfit = "AICc&k",
  correction.patch, correction.year, correction.indiv, alt_formula,
  alt_nocovsformula, alt_glmformula, extra_fac, null_model = FALSE) {
  
  old <- options()
  on.exit(options(old))
  
  model.null <- null.model.num <- NA
  override <- cutswitch <- FALSE
  
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
  } else if (vrate == 14) {
    if (!quiet) {message("\nDeveloping global model of juvenile maturity status...\n");}
    binom.model = TRUE
  } else {
    stop("Vital rate model not recognized.", call. = FALSE)
  }
  
  global.model <- .levindurosier(usedformula, subdata, approach, binom.model,
    dist, truncz, zero, quiet)
  
  if (any(class(global.model) == "try-error")) {
    if (!is.na(alt_formula)) {
      nox.model <- alt_formula
    } else {
      nox.model <- usedformula
    }
    
    if (nox.model != usedformula) {
      if (!quiet) {
        message("\nInitial global model estimation failed.
          Attempting a global model without interaction terms.\n")
      }
      global.model <- .levindurosier(nox.model, subdata, approach, binom.model,
        dist, truncz, zero, quiet)
    }
    
    if (any(class(global.model) == "try-error")) {
      nopat.model <- nox.model
      
      if (!is.na(correction.patch)) {
        nopat.model <- gsub(correction.patch, "", nopat.model, fixed = TRUE)
      }
      
      if (nox.model != nopat.model) {
        if (!quiet) {
          message("\nGlobal model estimation difficulties.
            Attempting a global model without a patch term.")
        }
        global.model <- .levindurosier(nopat.model, subdata, approach, binom.model,
          dist, truncz, zero, quiet)
      }
    }
    
    if (any(class(global.model) == "try-error")) {
      noyr.model <- nopat.model
      
      if (!is.na(correction.year)) {
        noyr.model <- gsub(correction.year, "", noyr.model, fixed = TRUE)
      }
      
      if (noyr.model != nopat.model) {
        if (!quiet) {
          message("\nGlobal model estimation difficulties.
            Attempting a global model without a year term.")
        }
        global.model <- .levindurosier(noyr.model, subdata, approach, binom.model,
          dist, truncz, zero, quiet)
      }
    }
    
    if (any(class(global.model) == "try-error")) {
      nocovs.model <- alt_nocovsformula
      
      if (!quiet) {
        message("\nGlobal model estimation difficulties.
          Attempting a global model without individual covariates.")
      }
      global.model <- .levindurosier(nocovs.model, subdata, approach, binom.model,
        dist, truncz, zero, quiet)
    }
    
    if (any(class(global.model) == "try-error")) {
      nocovpat.model <- nocovs.model
      
      if (!is.na(correction.patch)) {
        nocovpat.model <- gsub(correction.patch, "", nocovpat.model, fixed = TRUE)
      }
      
      if (nocovpat.model != nocovs.model) {
        if (!quiet) {
          message("\nGlobal model estimation difficulties.
            Attempting a global model without individual covariates and a patch term.")
        }
        global.model <- .levindurosier(nocovpat.model, subdata, approach, binom.model,
          dist, truncz, zero, quiet)
      }
    }
    
    if (any(class(global.model) == "try-error")) {
      nocovyr.model <- nocovpat.model
      
      if (!is.na(correction.year)) {
        nocovyr.model <- gsub(correction.year, "", nocovyr.model, fixed = TRUE)
      }
      
      if (nocovyr.model != nocovpat.model) {
        if (!quiet) {
          message("\nGlobal model estimation difficulties.
            Attempting a global model without individual covariates and patch and year terms.")
        }
        global.model <- .levindurosier(nocovyr.model, subdata, approach, binom.model,
          dist, truncz, zero, quiet)
      }
    }
    
    if (any(class(global.model) == "try-error")) {
      noind.model <- usedformula
      
      if (!is.na(correction.indiv)) {
        noind.model <- gsub(correction.indiv, "", noind.model, fixed = TRUE)
      }
      
      if (noind.model != usedformula) {
        if (!quiet) {
          message("\nGlobal model estimation difficulties.
            Attempting a global model without an individual identity term.")
        }
        global.model <- .levindurosier(noind.model, subdata, approach, binom.model,
          dist, truncz, zero, quiet)
      }
    }
    
    # The no random term version
    if (any(class(global.model) == "try-error") & approach == "mixed") {
      noran.model <- alt_glmformula
      
      if (!quiet) {
        message("\nGlobal model estimation difficulties.
          Attempting a global GLM model without random terms.")
      }
      global.model <- .levindurosier(noind.model, subdata, approach = "glm",
        binom.model, dist, truncz, zero, quiet)
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
        } else if (vrate == 14) {
          message("\nCould not properly estimate a global model for juvenile maturity status.")
        }
      }
      
      global.model <- 1
      vind <- 0
      vtrans <- 0
      model.bf <- global.model
      model.null <- global.model
      model.table <- NA
      
      cutswitch <- TRUE
    }
  }
  
  if (!cutswitch) {
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
    } else if (vrate == 14) {
      if (!quiet) {message("\nGlobal model of juvenile maturity status developed. Proceeding with model dredge...\n");}
    }
    
    model.table <- NA
    
    if (suite != "cons") override <- TRUE
    if (extra_fac > 0) override <- TRUE
    
    #This is the section where we dredge the models
    if (override & !global.only & !any(class(global.model) == "vglm")) {
  
      if (usedformula != 1) {
        options(na.action = "na.fail")
        if (!quiet) {
          model.table <- try(MuMIn::dredge(global.model, rank = criterion), silent = FALSE)
        } else {
          model.table <- suppressWarnings(suppressMessages(try(MuMIn::dredge(global.model, 
                  rank = criterion), silent = TRUE)))
        }
        null.model.num <- which(model.table$df == min(model.table$df))[1]
        
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
          } else if (vrate == 14) {
            message("\nDredge of juvenile maturity status failed.\n")
          }
        }
      }
    }
    
    #Here we extract the best-fit model
    if (length(grep("&k", bestfit)) > 0) {
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
          } else if (vrate == 14) {
            message("\nExtracting best-fit model of juvenile maturity status.\n")
          }
        }
        
        relevant.models <- which(model.table$delta <= 2)
        df.models <- which(model.table$df[relevant.models] == min(model.table$df[relevant.models]))
        if (length(df.models) > 1) {
          df.models <- min(df.models)
        }
        if (quiet) {
          model.bf <- suppressWarnings(suppressMessages(eval(stats::getCall(model.table, min(df.models)))))
          if (null_model) {
            model.null <- suppressWarnings(suppressMessages(eval(stats::getCall(model.table, null.model.num))))
          }
        } else {
          model.bf <- eval(stats::getCall(model.table, min(df.models)))
          if (null_model) {
            model.null <- eval(stats::getCall(model.table, null.model.num))
          }
        }
      } else {
        model.bf <- global.model
      }
      
    } else {
      if (any(class(model.table) == "model.selection")) {
        if (!quiet) {
          message("\nExtracting best-fit model.\n")
        }
        
        if (quiet) {
          model.bf <- suppressWarnings(suppressMessages(eval(stats::getCall(model.table, 1))))
          if (null_model) {
            model.null <- suppressWarnings(suppressMessages(eval(stats::getCall(model.table, null.model.num))))
          }
        } else {
          model.bf <- eval(stats::getCall(model.table, 1))
          if (null_model) {
            model.null <- eval(stats::getCall(model.table, null.model.num))
          }
        }
      } else {
        model.bf <- global.model
      }
    }
  }
  output <- list(model = global.model, ind = vind, trans = vtrans,
    bf.model = model.bf, null.model = model.null, table = model.table)
  
  return(output)
}

#' Core Global Model Builder for .headmaster_ritual()
#' 
#' A function that gets used repeatedly in \code{\link{.headmaster_ritual}()} to
#' build the global model to be dredged in function \code{\link{modelsearch}()}.
#' 
#' @name .levindurosier

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
#' @param quiet A logical value indicating whether warning messages should be
#' suppressed (\code{TRUE}) or displayed (\code{FALSE}).
#' 
#' @return This function returns a fit linear model of class generated by the
#' appropriate linear modeling function.
#' 
#' @keywords internal
#' @noRd
.levindurosier <- function(usedformula, subdata, approach, binom.model, dist,
  truncz, zero, quiet = FALSE) {
  
  if (!quiet) {
    if (approach == "mixed" & !is.na(dist)) {
      if (binom.model) {
        global.model <- try(lme4::glmer(formula = stats::as.formula(usedformula), 
            data = subdata, family = "binomial"), silent = quiet)
      } else {
        if (dist == "gaussian") {
          global.model <- try(lme4::lmer(formula = stats::as.formula(usedformula),
              data = subdata), silent = quiet)
        } else if (dist == "gamma") {
            global.model <- try(lme4::glmer(formula = stats::as.formula(usedformula),
                data = subdata, family = "Gamma"), silent = quiet)
        } else if (!truncz) {
          if (dist == "poisson" & !zero) {
            global.model <- try(lme4::glmer(formula = stats::as.formula(usedformula),
                data = subdata, family = "poisson"), silent = quiet)
          } else if (dist == "poisson" & zero) {
            global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(usedformula),
                data = subdata, ziformula=~., family = "poisson"), silent = quiet)
          } else if (dist == "negbin" & !zero) {
            global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(usedformula),
                data = subdata, ziformula=~0, family = glmmTMB::nbinom2), silent = quiet)
          } else if (dist == "negbin" & zero) {
            global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(usedformula),
                data = subdata, ziformula=~., family = glmmTMB::nbinom2), silent = quiet)
          }
        } else if (truncz) {
          if (dist == "poisson" & !zero) {
            global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(usedformula),
                data = subdata, family = glmmTMB::truncated_poisson), silent = quiet)
          } else if (dist == "negbin" & !zero) {
            global.model <- try(glmmTMB::glmmTMB(formula = stats::as.formula(usedformula),
                data = subdata, ziformula=~0, family = glmmTMB::truncated_nbinom2), silent = quiet)
          }
        }
      }
    } else if (approach == "glm" & !is.na(dist)) {
      if (binom.model) {
        global.model <- try(stats::glm(formula = stats::as.formula(usedformula),
            data = subdata, family = "binomial"), silent = quiet)
      } else {
        if (dist == "gaussian") {
          global.model <- try(stats::lm(formula = stats::as.formula(usedformula),
              data = subdata), silent = quiet)
        } else if (dist == "gamma") {
          global.model <- try(stats::glm(formula = stats::as.formula(usedformula),
              data = subdata, family = "Gamma"), silent = quiet)
        } else if (!truncz) {
          if (dist == "poisson" & !zero) {
            global.model <- try(stats::glm(formula = stats::as.formula(usedformula),
                data = subdata, family = "poisson"), silent = quiet)
          } else if (dist == "poisson" & zero) {
            global.model <- try(pscl::zeroinfl(formula = stats::as.formula(usedformula),
                data = subdata, dist = "poisson"), silent = quiet)
          } else if (dist == "negbin" & !zero) {
            global.model <- try(MASS::glm.nb(formula = stats::as.formula(usedformula),
                data = subdata), silent = quiet)
          } else if (dist == "negbin" & zero) {
            global.model <- try(pscl::zeroinfl(formula = stats::as.formula(usedformula), 
                data = subdata, dist = "negbin"), silent = quiet)
          }
        } else if (truncz) {
          usedformula <- gsub(" + 1", "", usedformula, fixed = TRUE)
          
          if (dist == "poisson" & !zero) {
            global.model <- try(VGAM::vglm(formula = stats::as.formula(usedformula),
                data = subdata, family = VGAM::pospoisson()), silent = quiet)
          } else if (dist == "negbin" & !zero) {
            global.model <- try(VGAM::vglm(formula = stats::as.formula(usedformula),
                data = subdata, family = VGAM::posnegbinomial()), silent = quiet)
          }
        }
      }
    } else {
      stop("Modeling approach not recognized.", call. = FALSE)
    }
  } else {
    if (approach == "mixed" & !is.na(dist)) {
      if (binom.model) {
        global.model <- suppressWarnings(suppressMessages(try(lme4::glmer(formula =
            stats::as.formula(usedformula), data = subdata, family = "binomial"),
            silent = quiet)))
      } else {
        if (dist == "gaussian") {
          global.model <- suppressWarnings(suppressMessages(try(lme4::lmer(formula =
              stats::as.formula(usedformula), data = subdata), silent = quiet)))
        } else if (dist == "gamma") {
            global.model <- suppressWarnings(suppressMessages(try(lme4::glmer(formula =
                stats::as.formula(usedformula), data = subdata, family = "Gamma"),
                silent = quiet)))
        } else if (!truncz) {
          if (dist == "poisson" & !zero) {
            global.model <- suppressWarnings(suppressMessages(try(lme4::glmer(formula =
                stats::as.formula(usedformula), data = subdata, family = "poisson"),
                silent = quiet)))
          } else if (dist == "poisson" & zero) {
            global.model <- suppressWarnings(suppressMessages(try(glmmTMB::glmmTMB(formula = 
                stats::as.formula(usedformula), data = subdata, ziformula=~., family = "poisson"),
                silent = quiet)))
          } else if (dist == "negbin" & !zero) {
            global.model <- suppressWarnings(suppressMessages(try(glmmTMB::glmmTMB(formula = 
                stats::as.formula(usedformula), data = subdata, ziformula=~0,
                family = glmmTMB::nbinom2), silent = quiet)))
          } else if (dist == "negbin" & zero) {
            global.model <- suppressWarnings(suppressMessages(try(glmmTMB::glmmTMB(formula = 
                stats::as.formula(usedformula), data = subdata, ziformula=~.,
                family = glmmTMB::nbinom2), silent = quiet)))
          }
        } else if (truncz) {
          if (dist == "poisson" & !zero) {
            global.model <- suppressWarnings(suppressMessages(try(glmmTMB::glmmTMB(formula = 
                stats::as.formula(usedformula), data = subdata,
                family = glmmTMB::truncated_poisson), silent = quiet)))
          } else if (dist == "negbin" & !zero) {
            global.model <- suppressWarnings(suppressMessages(try(glmmTMB::glmmTMB(formula = 
                stats::as.formula(usedformula), data = subdata, ziformula=~0,
                family = glmmTMB::truncated_nbinom2), silent = quiet)))
          }
        }
      }
    } else if (approach == "glm" & !is.na(dist)) {
      if (binom.model) {
        global.model <- suppressWarnings(suppressMessages(try(stats::glm(formula =
            stats::as.formula(usedformula), data = subdata, family = "binomial"),
            silent = quiet)))
      } else {
        if (dist == "gaussian") {
          global.model <- suppressWarnings(suppressMessages(try(stats::lm(formula =
              stats::as.formula(usedformula), data = subdata), silent = quiet)))
        } else if (dist == "gamma") {
          global.model <- suppressWarnings(suppressMessages(try(stats::glm(formula =
              stats::as.formula(usedformula), data = subdata, family = "Gamma"),
              silent = quiet)))
        } else if (!truncz) {
          if (dist == "poisson" & !zero) {
            global.model <- suppressWarnings(suppressMessages(try(stats::glm(formula =
                stats::as.formula(usedformula), data = subdata, family = "poisson"),
                silent = quiet)))
          } else if (dist == "poisson" & zero) {
            global.model <- suppressWarnings(suppressMessages(try(pscl::zeroinfl(formula = 
                stats::as.formula(usedformula), data = subdata, dist = "poisson"),
                silent = quiet)))
          } else if (dist == "negbin" & !zero) {
            global.model <- suppressWarnings(suppressMessages(try(MASS::glm.nb(formula = 
                stats::as.formula(usedformula), data = subdata), silent = quiet)))
          } else if (dist == "negbin" & zero) {
            global.model <- suppressWarnings(suppressMessages(try(pscl::zeroinfl(formula = 
                stats::as.formula(usedformula), data = subdata, dist = "negbin"),
                silent = quiet)))
          }
        } else if (truncz) {
          usedformula <- gsub(" + 1", "", usedformula, fixed = TRUE)
          
          if (dist == "poisson" & !zero) {
            global.model <- suppressWarnings(suppressMessages(try(VGAM::vglm(formula =
                stats::as.formula(usedformula), data = subdata, family = VGAM::pospoisson()),
                silent = quiet)))
          } else if (dist == "negbin" & !zero) {
            global.model <- suppressWarnings(suppressMessages(try(VGAM::vglm(formula =
                stats::as.formula(usedformula), data = subdata, family = VGAM::posnegbinomial()),
                silent = quiet)))
          }
        }
      }
    } else {
      stop("Modeling approach not recognized.", call. = FALSE)
    }
  }
  return(global.model)
}

#' Estimate Accuracy of Binomial Model
#' 
#' Function \code{.accu_predict} estimates the accuracy of vital rate models.
#' Accuracy of vital rate models is calculated differently depending on vital
#' rate and assumed distribution. For all vital rates assuming a binomial
#' distribution, including survival, observation status, reproductive status,
#' and juvenile version of these, accuracy is calculated as the percent of
#' predicted responses equal to actual responses. In all other models, accuracy
#' is actually the conditional R\textsuperscript{2} using package \code{MuMIn}'s
#' \code{\link[MuMIn]{r.squaredGLMM}()} function, estimated via the delta
#' method. When this method fails, \code{modelsearch()} calculates McFadden's
#' pseudo-R\textsuperscript{2}. If this fails, then \code{NA} is returned.
#' #' currently only handles this for binomial models, and performs it as a
#' comparison of the actual and predicted responses from the respective model.
#' 
#' @name .accu_predict
#' 
#' @param bestfitmodel The best-fit model to be passed.
#' @param subdata The dataset to be used for accuracy testing.
#' @param param The name of the response parameter to be used in accuracy
#' testing.
#' @param style An integer indicating whether to calculate binomial accuracy (1)
#' or package \code{MuMIn}'s conditional R2 (2). Defaults to \code{1}.
#' @param nullmodel A null (y-intercept only) model from the model table. Only
#' used to calculate McFadden's pseudo-R2. Defaults to \code{NA}.
#' @param quiet A logical value indicating whether to allow warning messages to
#' be displayed (\code{FALSE}) or suppressed (\code{TRUE}).
#' Defaults to \code{FALSE}.
#' @param check A logical value indicating whether to test accuracy. Defaults to
#' \code{TRUE}.
#' 
#' @return A single numeric value giving the proportion of responses accurately
#' predicted.
#' 
#' @keywords internal
#' @noRd
.accu_predict <- function(bestfitmodel, subdata = NA, param, style = 1,
  nullmodel = NA, quiet = FALSE, check = TRUE) {
  
  pred_vec <- NULL
  accuracy <- NA 
  
  if (check & !any(class(bestfitmodel) == "vglm")) {
    if (!is.numeric(bestfitmodel) & !is.element("NULL", class(bestfitmodel))) {
      if (style == 1) {
        
        pred_vec <- stats::predict(bestfitmodel, newdata = subdata,
          type = "response")
        pred_vec <- round(pred_vec)
        
        test_vec <- subdata[,which(names(subdata) == param)]
        
        results_vec <- pred_vec - test_vec
        
        accuracy <- length(which(results_vec == 0)) / length(results_vec)
        
      } else if (style == 2) {
        
        if (quiet) {
          accuracy.table <- suppressWarnings(suppressMessages(try(MuMIn::r.squaredGLMM(bestfitmodel),
              silent = quiet)))
        } else {
          accuracy.table <- try(MuMIn::r.squaredGLMM(bestfitmodel), silent = quiet)
        }
        
        if (any(class(accuracy.table) == "try-error")) {
          if (all(is.numeric(bestfitmodel))) {
            all_dat_values <- unique(subdata[,which(names(subdata) == param)])
            if (length(all_dat_values) == 1) {
              accuracy <- 1
            }
          } else {
            bf_log <- logLik(bestfitmodel)[1]
            null_log <- logLik(nullmodel)[1]
          
            accuracy <- 1 - (bf_log / null_log)
          }
          
          if (is.numeric(accuracy)) {
            if (!is.na(accuracy)) {
              if (accuracy < 0 | accuracy > 1) accuracy <- NA
            }
          } else {accuracy <- NA}
          
        } else if (!suppressWarnings(all(is.na(nullmodel)))) {
          if (is.element("delta", rownames(accuracy.table))) {
            accuracy <- accuracy.table["delta","R2c"]
          } else {
            accuracy <- accuracy.table[1,"R2c"]
          }
        }
      }
    }
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
#' @name summary.lefkoMod
#' 
#' @param object An R object of class \code{lefkoMod} resulting from
#' \code{\link{modelsearch}()}.
#' @param ... Other parameters currently not utilized.
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
    is.numeric(modelsuite$juv_sizec_model), is.numeric(modelsuite$juv_reproduction_model),
    is.numeric(modelsuite$juv_maturity_model)) == FALSE))
  
  writeLines(paste0("This LefkoMod object includes ", totalmodels, " linear models."))
  writeLines(paste0("Best-fit model criterion used: ", modelsuite$criterion))
  writeLines(paste0("\n\n"))
  writeLines("Survival model:")
  print(modelsuite$survival_model)
  writeLines(paste0("\n"))
  writeLines("\nObservation model:")
  print(modelsuite$observation_model)
  writeLines(paste0("\n"))
  writeLines("\nSize model:")
  print(modelsuite$size_model)
  writeLines(paste0("\n"))
  writeLines("\nSecondary size model:")
  print(modelsuite$sizeb_model)
  writeLines(paste0("\n"))
  writeLines("\nTertiary size model:")
  print(modelsuite$sizec_model)
  writeLines(paste0("\n"))
  writeLines("\nReproductive status model:")
  print(modelsuite$repstatus_model)
  writeLines(paste0("\n"))
  writeLines("\nFecundity model:")
  print(modelsuite$fecundity_model)
  writeLines(paste0("\n"))
  writeLines("Juvenile survival model:")
  print(modelsuite$juv_survival_model)
  writeLines(paste0("\n"))
  writeLines("\nJuvenile observation model:")
  print(modelsuite$juv_observation_model)
  writeLines(paste0("\n"))
  writeLines("\nJuvenile size model:")
  print(modelsuite$juv_size_model)
  writeLines(paste0("\n"))
  writeLines("\nJuvenile secondary size model:")
  print(modelsuite$juv_sizeb_model)
  writeLines(paste0("\n"))
  writeLines("\nJuvenile tertiary size model:")
  print(modelsuite$juv_sizec_model)
  writeLines(paste0("\n"))
  writeLines("\nJuvenile reproduction model:")
  print(modelsuite$juv_reproduction_model)
  writeLines(paste0("\n"))
  writeLines("\nJuvenile maturity model:")
  print(modelsuite$juv_maturity_model)
  
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
  
  if (!is.logical(modelsuite$juv_maturity_model) & 
      is.element("model.selection", class(modelsuite$juv_maturity_table))) {
    writeLines(paste0("\nNumber of models in juvenile maturity table: ", 
        dim(modelsuite$juv_maturity_table)[1]))
  } else if (!is.logical(modelsuite$juv_maturity_model)) {
    writeLines("\nNumber of models in juvenile maturity table: 1")
  } else {
    writeLines("\nJuvenile maturity table not estimated")
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
    writeLines(paste0("Survival accuracy is ", round(modelsuite$qc[1,"accuracy"], 3), "."))
  } else {
    writeLines("Survival not estimated.")
  }
  
  if (modelsuite$qc[2,2] > 0) {
    writeLines(paste0("Observation estimated with ", modelsuite$qc[2,2], 
        " individuals and ", modelsuite$qc[2,3], " individual transitions."))
    writeLines(paste0("Observation accuracy is ", round(modelsuite$qc[2,"accuracy"], 3), "."))
  } else {
    writeLines("Observation probability not estimated.")
  }
  
  if (modelsuite$qc[3,2] > 0) {
    writeLines(paste0("Primary size estimated with ", modelsuite$qc[3,2], 
        " individuals and ", modelsuite$qc[3,3], " individual transitions."))
    writeLines(paste0("Primary size pseudo R-squared is ",
      round(modelsuite$qc[3,"accuracy"], 3), "."))
  } else {
    writeLines("Primary size transition not estimated.")
  }
  
  if (modelsuite$qc[4,2] > 0) {
    writeLines(paste0("Secondary size estimated with ", modelsuite$qc[4,2], 
        " individuals and ", modelsuite$qc[4,3], " individual transitions."))
    writeLines(paste0("Secondary size pseudo R-squared is ",
      round(modelsuite$qc[4,"accuracy"], 3), "."))
  } else {
    writeLines("Secondary size transition not estimated.")
  }
  
  if (modelsuite$qc[5,2] > 0) {
    writeLines(paste0("Tertiary size estimated with ", modelsuite$qc[5,2], 
        " individuals and ", modelsuite$qc[5,3], " individual transitions."))
    writeLines(paste0("Tertiary size pseudo R-squared is ",
      round(modelsuite$qc[5,"accuracy"], 3), "."))
  } else {
    writeLines("Tertiary size transition not estimated.")
  }
  
  if (modelsuite$qc[6,2] > 0) {
    writeLines(paste0("Reproductive status estimated with ", modelsuite$qc[6,2], 
        " individuals and ", modelsuite$qc[6,3], " individual transitions."))
   writeLines(paste0("Reproductive status accuracy is ", round(modelsuite$qc[6,"accuracy"], 3), "."))
   } else {
    writeLines("Reproduction probability not estimated.")
  }
  
  if (modelsuite$qc[7,2] > 0) {
    writeLines(paste0("Fecundity estimated with ", modelsuite$qc[7,2], 
        " individuals and ", modelsuite$qc[7,3], " individual transitions."))
    writeLines(paste0("Fecundity pseudo R-squared is ",
        round(modelsuite$qc[7,"accuracy"], 3), "."))
  } else {
    writeLines("Fecundity not estimated.")
  }
  
  if (modelsuite$qc[8,2] > 0) {
    writeLines(paste0("Juvenile survival estimated with ", modelsuite$qc[8,2], 
        " individuals and ", modelsuite$qc[8,3], " individual transitions."))
   writeLines(paste0("Juvenile survival accuracy is ", round(modelsuite$qc[8,"accuracy"], 3), "."))
   } else {
    writeLines("Juvenile survival not estimated.")
  }
  
  if (modelsuite$qc[9,2] > 0) {
    writeLines(paste0("Juvenile observation estimated with ", modelsuite$qc[9,2], 
        " individuals and ", modelsuite$qc[9,3], " individual transitions."))
   writeLines(paste0("Juvenile observation accuracy is ", round(modelsuite$qc[9,"accuracy"], 3), "."))
  } else {
    writeLines("Juvenile observation probability not estimated.")
  }
  
  if (modelsuite$qc[10,2] > 0) {
    writeLines(paste0("Juvenile primary size estimated with ", modelsuite$qc[10,2], 
        " individuals and ", modelsuite$qc[10,3], " individual transitions."))
    writeLines(paste0("Juvenile primary size pseudo R-squared is ",
      round(modelsuite$qc[10,"accuracy"], 3), "."))
  } else {
    writeLines("Juvenile primary size transition not estimated.")
  }
  
  if (modelsuite$qc[11,2] > 0) {
    writeLines(paste0("Juvenile secondary size estimated with ", modelsuite$qc[11,2],
        " individuals and ", modelsuite$qc[11,3], " individual transitions."))
    writeLines(paste0("Juvenile secondary size pseudo R-squared is ",
      round(modelsuite$qc[11,"accuracy"], 3), "."))
  } else {
    writeLines("Juvenile secondary size transition not estimated.")
  }
  
  if (modelsuite$qc[12,2] > 0) {
    writeLines(paste0("Juvenile tertiary size estimated with ", modelsuite$qc[12,2],
        " individuals and ", modelsuite$qc[12,3], " individual transitions."))
    writeLines(paste0("Juvenile tertiary size pseudo R-squared is ",
      round(modelsuite$qc[12,"accuracy"], 3), "."))
  } else {
    writeLines("Juvenile tertiary size transition not estimated.")
  }
  
  if (modelsuite$qc[13,2] > 0) {
    writeLines(paste0("Juvenile reproduction estimated with ", modelsuite$qc[13,2],
        " individuals and ", modelsuite$qc[13,3], " individual transitions."))
    writeLines(paste0("Juvenile reproductive status accuracy is ",
        round(modelsuite$qc[13,"accuracy"], 3), "."))
  } else {
    writeLines("Juvenile reproduction probability not estimated.")
  }
  
  if (modelsuite$qc[14,2] > 0) {
    writeLines(paste0("Juvenile maturity estimated with ", modelsuite$qc[14,2],
        " individuals and ", modelsuite$qc[14,3], " individual transitions."))
    writeLines(paste0("Juvenile maturity status accuracy is ",
        round(modelsuite$qc[14,"accuracy"], 3), "."))
  } else {
    writeLines("Juvenile maturity probability not estimated.")
  }
}

#' Import Vital Rate Model Factor Values for Function-based MPM Development
#' 
#' Function \code{vrm_import()} builds a skeleton list holding data frames and
#' vectors that can be used to import coefficient values for the factors of the
#' vital rate models used to build function-based MPMs or run function-based
#' projections.
#' 
#' @name vrm_import
#' 
#' @param years A numeric vector of the years or times at time \code{t} to be
#' modeled.
#' @param patches A string or numeric vector of the patch names to be modeled.
#' @param groups An integer vector of stage groups to be modeled. Defaults to a
#' vector with a single element with value \code{0}.
#' @param interactions A logical value indicating whether to include two-way
#' interactions between main effects (\code{TRUE}), or only main effects
#' (\code{FALSE}). Defaults to \code{FALSE}.
#' @param zi A logical value indicating whether to include coefficients for the
#' binomial components of zero-inflation models. Defaults to \code{FALSE}.
#' @param cat.indcova If individual covariate a is categorical, then this term
#' should equal a string vector of the names of the categories. Defaults to
#' \code{NULL}, in which case individual covariate a is either not used or is
#' numeric.
#' @param cat.indcovb If individual covariate b is categorical, then this term
#' should equal a string vector of the names of the categories. Defaults to
#' \code{NULL}, in which case individual covariate b is either not used or is
#' numeric.
#' @param cat.indcovc If individual covariate c is categorical, then this term
#' should equal a string vector of the names of the categories. Defaults to
#' \code{NULL}, in which case individual covariate c is either not used or is
#' numeric.
#' @param dist.sizea A string value giving the distribution of the variable
#' coding primary size. Can equal \code{"none"}, \code{"gamma"},
#' \code{"gaussian"}, \code{"poisson"}, \code{"negbin"}, or \code{"constant"}.
#' Defaults to \code{"gaussian"}.
#' @param dist.sizeb A string value giving the distribution of the variable
#' coding secondary size. Can equal \code{"none"}, \code{"gamma"},
#' \code{"gaussian"}, \code{"poisson"}, \code{"negbin"}, or \code{"constant"}.
#' Defaults to \code{"constant"}.
#' @param dist.sizec A string value giving the distribution of the variable
#' coding tertiary size. Can equal \code{"none"}, \code{"gamma"},
#' \code{"gaussian"}, \code{"poisson"}, \code{"negbin"}, or \code{"constant"}.
#' Defaults to \code{"constant"}.
#' @param dist.fec A string value giving the distribution of the variable
#' coding fecundity. Can equal \code{"none"}, \code{"gamma"},
#' \code{"gaussian"}, \code{"poisson"}, or \code{"negbin"}. Defaults to
#' \code{"gaussian"}.
#' @param trunc.sizea A logical value indicating whether the distribution of the
#' primary size variable should be zero-truncated. Defaults to \code{FALSE}.
#' Currently only works with the Poisson and negative binomial distributions.
#' @param trunc.sizeb A logical value indicating whether the distribution of the
#' secondary size variable should be zero-truncated. Defaults to \code{FALSE}.
#' Currently only works with the Poisson and negative binomial distributions.
#' @param trunc.sizec A logical value indicating whether the distribution of the
#' tertiary size variable should be zero-truncated. Defaults to \code{FALSE}.
#' Currently only works with the Poisson and negative binomial distributions.
#' @param trunc.fec A logical value indicating whether the distribution of the
#' fecundity variable should be zero-truncated. Defaults to \code{FALSE}.
#' Currently only works with the Poisson and negative binomial distributions.
#' @param use.juv A logical value indicating whether to utilize juvenile vital
#' rates. If \code{FALSE}, then all juvenile vital rates will be set to
#' \code{constant} distributions. Defaults to \code{FALSE}.
#' 
#' @return A list of class \code{vrm_input}, with up to 13 elements including:
#' \item{vrm_frame}{A data frame holding the main slope coefficients for the
#' linear vital rate models.}
#' \item{year_frame}{A data frame holding the main slope coefficients for the
#' year at time t terms in the linear vital rate models.}
#' \item{patch_frame}{A data frame holding the main slope coefficients for the
#' patch terms in the linear vital rate models.}
#' \item{group2_frame}{A data frame holding the main slope coefficients for the
#' stage group terms in time \emph{t} in the linear vital rate models.}
#' \item{group1_frame}{A data frame holding the main slope coefficients for the
#' stage group terms in time \emph{t}-1 in the linear vital rate models.}
#' \item{dist_frame}{A data frame giving the distributions of all variables,
#' including primary, secondary, and tertiary size, and fecundity. Some
#' variables begin as \code{constant}.}
#' \item{indcova2_frame}{A data frame holding the main slope coefficients for
#' the categorical individual covariate a terms in time \emph{t} in the linear
#' vital rate models.}
#' \item{indcova1_frame}{A data frame holding the main slope coefficients for
#' the categorical individual covariate a terms in time \emph{t}-1 in the linear
#' vital rate models.}
#' \item{indcovb2_frame}{A data frame holding the main slope coefficients for
#' the categorical individual covariate b terms in time \emph{t} in the linear
#' vital rate models.}
#' \item{indcovb1_frame}{A data frame holding the main slope coefficients for
#' the categorical individual covariate b terms in time \emph{t}-1 in the linear
#' vital rate models.}
#' \item{indcovc2_frame}{A data frame holding the main slope coefficients for
#' the categorical individual covariate c terms in time \emph{t} in the linear
#' vital rate models.}
#' \item{indcovc1_frame}{A data frame holding the main slope coefficients for
#' the categorical individual covariate c terms in time \emph{t}-1 in the linear
#' vital rate models.}
#' \item{st_frame}{A data frame holding values of sigma or theta for use in
#' Gaussian or negative binomial response terms, respectively.}
#' 
#' The first element, called \code{vrm_frame}, is a data frame with the
#' following 18 variables:
#' \item{main_effect_1}{The main effect for which coefficients are to be
#' entered.}
#' \item{main_1_defined}{A more natural explanation of \code{main_effect_1}.}
#' \item{main_effect_2}{If given, then indicates another effect in a two-way
#' interaction with \code{main_effect_1}.}
#' \item{main_2_defined}{A more natural explanation of \code{main_effect_2}.}
#' \item{surv}{A vector of coefficients for the factors in the model of adult
#' survival.}
#' \item{obs}{A vector of coefficients for the factors in the model of adult
#' observation status.}
#' \item{sizea}{A vector of coefficients for the factors in the model of adult
#' primary size.}
#' \item{sizeb}{A vector of coefficients for the factors in the model of adult
#' secondary size.}
#' \item{sizec}{A vector of coefficients for the factors in the model of adult
#' tertiary size.}
#' \item{repst}{A vector of coefficients for the factors in the model of adult
#' reproductive status.}
#' \item{fec}{A vector of coefficients for the factors in the model of adult
#' fecundity.}
#' \item{jsurv}{A vector of coefficients for the factors in the model of
#' juvenile survival.}
#' \item{jobs}{A vector of coefficients for the factors in the model of juvenile
#' observation status.}
#' \item{jsize}{A vector of coefficients for the factors in the model of
#' juvenile primary size.}
#' \item{jsizeb}{A vector of coefficients for the factors in the model of
#' juvenile secondary size.}
#' \item{jsizec}{A vector of coefficients for the factors in the model of
#' juvenile tertiary size.}
#' \item{jrepst}{A vector of coefficients for the factors in the model of
#' juvenile reproductive status, for individuals maturing in the current time
#' step.}
#' \item{jmat}{A vector of coefficients for the factors in the model of maturity
#' status, for individuals capable of maturing at the current time step.}
#' \item{sizea_zi}{A vector of coefficients for the factors in the binomial
#' component of the zero-inflated model of adult primary size, if zero-inflated
#' models are being used.}
#' \item{sizeb_zi}{A vector of coefficients for the factors in the binomial
#' component of the zero-inflated model of adult secondary size, if
#' zero-inflated models are being used.}
#' \item{sizec_zi}{A vector of coefficients for the factors in the binomial
#' component of the zero-inflated model of adult tertiary size, if zero-inflated
#' models are being used.}
#' \item{fec_zi}{A vector of coefficients for the factors in the binomial
#' component of the zero-inflated model of fecundity, if zero-inflated models
#' are being used.}
#' \item{jsizea_zi}{A vector of coefficients for the factors in the binomial
#' component of the zero-inflated model of juvenile primary size, if
#' zero-inflated models are being used.}
#' \item{jsizeb_zi}{A vector of coefficients for the factors in the binomial
#' component of the zero-inflated model of juvenile secondary size, if
#' zero-inflated models are being used.}
#' \item{jsizec_zi}{A vector of coefficients for the factors in the binomial
#' component of the zero-inflated model of juvenile tertiary size, if
#' zero-inflated models are being used.}
#' 
#' @section Notes:
#' All coefficients across all data frames are initially set to \code{0}. After
#' using this function to create the skeleton list, all relevant coefficient
#' values should be set to non-zero values equal to the respective slope from
#' the appropriate linear model, and any vital rate model to be used should
#' have its distribution set to \code{"binom"}, \code{"gaussian"},
#' \code{"gamma"}, \code{"poisson"}, or \code{"negbin"}. Unused vital rates
#' should be set to \code{"constant"}, and the first element of the correspoding
#' column in \code{vrm_frame} (corresponding to the y-intercept) should be set
#' to the constant value to utilize (generally \code{1}). If no values are
#' manually edited, then function-based MPM generator functions will not be able
#' to generate valid MPMs.
#' 
#' Users should never change the labels or the order of terms in the data frames
#' and vectors produced via this function, nor should they ever changes the
#' names of core list elements in the \code{vrm_input} object. Doing so will
#' result either in fatal errors or erroneous matrix calculations.
#' 
#' Using the \code{vrm_input} approach to building function-based MPMs requires
#' careful attention to the stageframe. Although no hfv data frame needs to be
#' entered, stages for which vital rates are to be estimated via linear models
#' parameterized with coefficients provided via function \code{vrm_import()}
#' should be marked as occurring within the dataset, while stages for which
#' the provided coefficients should not be used should be marked as not
#' occurring within the dataset.
#' 
#' Coefficients added to zero-inflation models can only be added to primary
#' size, secondary size, tertiary size, fecundity, and the juvenile versions of
#' primary, secondary, and tertiary size. Care must be taken to include zero-
#' inflated coefficients only for variables without size-truncated
#' distributions. Adding such terms will result in fatal errors during matrix
#' creation.
#' 
#' @examples
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 0, 1, 7100)
#' stagevector <- c("Sd", "Sdl", "Dorm", "ipm", "ipm")
#' repvector <- c(0, 0, 0, 1, 1)
#' obsvector <- c(0, 1, 0, 1, 1)
#' matvector <- c(0, 0, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1)
#' binvec <- c(0, 100, 0.5, 1, 1)
#' comments <- c("Dormant seed", "Seedling", "Dormant", "ipm adult stage",
#'   "ipm adult stage")
#' lathframeipm <- sf_create(sizes = sizevector, stagenames = stagevector, 
#'   repstatus = repvector, obsstatus = obsvector, propstatus = propvector,
#'   immstatus = immvector, matstatus = matvector, comments = comments,
#'   indataset = indataset, binhalfwidth = binvec, ipmbins = 100, roundsize = 3)
#' 
#' lathsupp2 <- supplemental(stage3 = c("Sd", "Sdl", "Sd", "Sdl"),
#'   stage2 = c("Sd", "Sd", "rep", "rep"),
#'   givenrate = c(0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 3, 3), stageframe = lathframeipm, historical = FALSE)
#' 
#' lath_vrm <- vrm_import(years = c(1988:1990), zi = TRUE, dist.fec = "negbin",
#'   use.juv = TRUE)
#' 
#' lath_vrm$vrm_frame$surv[1] <- 2.32571
#' lath_vrm$vrm_frame$surv[2] <- 0.00109
#' lath_vrm$vrm_frame$obs[1] <- 2.230
#' lath_vrm$vrm_frame$sizea[1] <- 164.0695
#' lath_vrm$vrm_frame$sizea[2] <- 0.6211
#' lath_vrm$vrm_frame$fec[1] <- 1.517
#' lath_vrm$vrm_frame$fec_zi[1] <- 6.252765
#' lath_vrm$vrm_frame$fec_zi[2] <- -0.007313
#' lath_vrm$vrm_frame$jsurv[1] <- 1.03
#' lath_vrm$vrm_frame$jobs[1] <- 10.390
#' lath_vrm$vrm_frame$jsizea[1] <- 3.0559
#' lath_vrm$vrm_frame$jsizea[2] <- 0.8482
#' 
#' lath_vrm$year_frame$fec[c(1:3)] <- c(-0.41749627, 0.51421684, -0.07964038)
#' lath_vrm$year_frame$fec_zi[c(1:3)] <- c(3.741475e-07, -7.804715e-08,
#'   -2.533755e-07)
#' lath_vrm$year_frame$sizea[c(1:3)] <- c(96.3244, -240.8036, 144.4792)
#' lath_vrm$year_frame$jobs[c(1:3)] <- c(-0.7459843, 0.6118826, -0.9468618)
#' lath_vrm$year_frame$jsizea[c(1:3)] <- c(0.5937962, 1.4551236, -2.0489198)
#' 
#' lath_vrm$dist_frame$dist[2] <- "binom"
#' lath_vrm$dist_frame$dist[9] <- "binom"
#' 
#' lath_vrm$st_frame[3] <- 503.6167
#' lath_vrm$st_frame[7] <- 0.2342114
#' lath_vrm$st_frame[10] <- 5.831
#' 
#' lath_vrm$vrm_frame$sizeb[1] <- 1
#' lath_vrm$vrm_frame$sizec[1] <- 1
#' lath_vrm$vrm_frame$repst[1] <- 1
#' lath_vrm$vrm_frame$jsizeb[1] <- 1
#' lath_vrm$vrm_frame$jsizec[1] <- 1
#' lath_vrm$vrm_frame$jrepst[1] <- 1
#' lath_vrm$vrm_frame$jmatst[1] <- 1
#' 
#' lathmat2_importipm <- flefko2(stageframe = lathframeipm,
#'   modelsuite = lath_vrm, supplement = lathsupp2, reduce = FALSE)
#' 
#' summary(lathmat2_importipm)
#' 
#' @export
vrm_import <- function(years = NULL, patches = c(1), groups = c(0),
  interactions = FALSE, zi = FALSE, cat.indcova = NULL, cat.indcovb = NULL,
  cat.indcovc = NULL, dist.sizea = "gaussian", dist.sizeb = "constant",
  dist.sizec = "constant", dist.fec = "gaussian", trunc.sizea = FALSE,
  trunc.sizeb = FALSE, trunc.sizec = FALSE, trunc.fec = FALSE,
  use.juv = FALSE) {
  
  dist.sizea <- tolower(dist.sizea)
  dist.sizeb <- tolower(dist.sizeb)
  dist.sizec <- tolower(dist.sizec)
  dist.fec <- tolower(dist.fec)
  
  dist.possibilities <- c("none", "gamma", "gaussian", "poisson", "negbin",
    "constant")
  
  if (length(dist.sizea) != 1 | !is.element(dist.sizea[1], dist.possibilities)) {
    stop("Option dist.sizea must be set to a single value. Options include:
      none, gamma, gaussian, poisson, negbin, or constant.", call. = FALSE)
  }
  if (length(dist.sizeb) != 1 | !is.element(dist.sizeb[1], dist.possibilities)) {
    stop("Option dist.sizeb must be set to a single value. Options include:
      none, gamma, gaussian, poisson, negbin, or constant.", call. = FALSE)
  }
  if (length(dist.sizec) != 1 | !is.element(dist.sizec[1], dist.possibilities)) {
    stop("Option dist.sizec must be set to a single value. Options include:
      none, gamma, gaussian, poisson, negbin, or constant.", call. = FALSE)
  }
  if (length(dist.fec) != 1 | !is.element(dist.fec[1], dist.possibilities)) {
    stop("Option dist.fec must be set to a single value. Options include:
      none, gamma, gaussian, poisson, negbin, or constant.", call. = FALSE)
  }
  
  if (trunc.sizea) {
    if (dist.sizea == "poisson") {
      dist.sizea <- "truncated_poisson"
    } else if (dist.sizea == "negbin") {
      dist.sizea <- "truncated_negbin"
    } else {
      stop("Zero truncation can currently be incorporated only in Poisson and
          negative binomial models.", call. = FALSE)
    }
  }
  if (trunc.sizeb) {
    if (dist.sizeb == "poisson") {
      dist.sizeb <- "truncated_poisson"
    } else if (dist.sizeb == "negbin") {
      dist.sizeb <- "truncated_negbin"
    } else {
      stop("Zero truncation can currently be incorporated only in Poisson and
          negative binomial models.", call. = FALSE)
    }
  }
  if (trunc.sizec) {
    if (dist.sizec == "poisson") {
      dist.sizec <- "truncated_poisson"
    } else if (dist.sizec == "negbin") {
      dist.sizec <- "truncated_negbin"
    } else {
      stop("Zero truncation can currently be incorporated only in Poisson and
          negative binomial models.", call. = FALSE)
    }
  }
  if (trunc.fec) {
    if (dist.fec == "poisson") {
      dist.fec <- "truncated_poisson"
    } else if (dist.fec == "negbin") {
      dist.fec <- "truncated_negbin"
    } else {
      stop("Zero truncation can currently be incorporated only in Poisson and
          negative binomial models.", call. = FALSE)
    }
  }
  
  if (is.null(years)) stop("Option years must equal a numeric vector.", call. = FALSE);
  if (!is.numeric(years)) stop("Option years must equal a numeric vector.", call. = FALSE);
  
  num_years <- length(years)
  num_patches <- length(patches)
  num_groups <- length(groups)
  
  main_effects <- c("intercept", "size2", "size1", "sizeb2", "sizeb1", "sizec2",
    "sizec1", "repst2", "repst1", "age", "density", "indcova2", "indcova1",
    "indcovb2", "indcovb1", "indcovc2", "indcovc1")
  main_defined <- c("y-intercept", "sizea in time t", "sizea in time t-1",
    "sizeb in time t", "sizeb in time t-1", "sizec in time t", "sizec in time t-1",
    "reproductive status in time t", "reproductive status in time t-1",
    "age in time t", "density in time t", "individual covariate a in time t",
    "individual covariate a in time t-1", "individual covariate b in time t",
    "individual covariate b in time t-1", "individual covariate c in time t",
    "individual covariate c in time t-1")
  
  main1 <- main_effects
  main1_def <- main_defined
  
  basic_length <- 17
  if (interactions) {
    main2 <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
    main2_def <- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "")
    
    basic_length <- basic_length + 111
    
    extended_terms_1 <- c(main_effects[9], main_effects[3], main_effects[3],
      main_effects[2], main_effects[2], main_effects[3], main_effects[10],
      main_effects[10], main_effects[10], main_effects[10], main_effects[12],
      main_effects[14], main_effects[16], main_effects[12], main_effects[14],
      main_effects[16], main_effects[13], main_effects[15], main_effects[17],
      main_effects[13], main_effects[15], main_effects[17], main_effects[12],
      main_effects[12], main_effects[14], main_effects[13], main_effects[13],
      main_effects[15], main_effects[12], main_effects[13], main_effects[12],
      main_effects[13], main_effects[14], main_effects[15], main_effects[4],
      main_effects[6], main_effects[3], main_effects[3], main_effects[5],
      main_effects[2], main_effects[2], main_effects[4], main_effects[3],
      main_effects[3], main_effects[5], main_effects[2], main_effects[2],
      main_effects[4], main_effects[11], main_effects[11], main_effects[11],
      main_effects[11], main_effects[11], main_effects[11], main_effects[11],
      main_effects[11], main_effects[4], main_effects[6], main_effects[5],
      main_effects[4], main_effects[5], main_effects[7], main_effects[6],
      main_effects[7], main_effects[4], main_effects[6], main_effects[11],
      main_effects[5], main_effects[7], main_effects[12], main_effects[12],
      main_effects[12], main_effects[13], main_effects[13], main_effects[13],
      main_effects[13], main_effects[12], main_effects[12], main_effects[13],
      main_effects[14], main_effects[14], main_effects[14], main_effects[15],
      main_effects[15], main_effects[15], main_effects[15], main_effects[14],
      main_effects[14], main_effects[15], main_effects[16], main_effects[16],
      main_effects[16], main_effects[17], main_effects[17], main_effects[17],
      main_effects[17], main_effects[16], main_effects[16], main_effects[17],
      main_effects[12], main_effects[14], main_effects[16], main_effects[13],
      main_effects[15], main_effects[17], main_effects[12], main_effects[14],
      main_effects[16], main_effects[13], main_effects[15], main_effects[17])
    
    extended_terms_2 <- c(main_effects[8], main_effects[2], main_effects[9],
      main_effects[8], main_effects[9], main_effects[8], main_effects[3],
      main_effects[2], main_effects[9], main_effects[8], main_effects[2],
      main_effects[2], main_effects[2], main_effects[8], main_effects[8],
      main_effects[8], main_effects[3], main_effects[3], main_effects[3],
      main_effects[9], main_effects[9], main_effects[9], main_effects[14],
      main_effects[16], main_effects[16], main_effects[15], main_effects[17],
      main_effects[17], main_effects[15], main_effects[14], main_effects[17],
      main_effects[16], main_effects[17], main_effects[16], main_effects[5],
      main_effects[7], main_effects[5], main_effects[7], main_effects[7],
      main_effects[4], main_effects[6], main_effects[6], main_effects[4],
      main_effects[6], main_effects[6], main_effects[5], main_effects[7],
      main_effects[7], main_effects[2], main_effects[4], main_effects[6],
      main_effects[3], main_effects[5], main_effects[7], main_effects[8],
      main_effects[9], main_effects[8], main_effects[8], main_effects[9],
      main_effects[9], main_effects[8], main_effects[9], main_effects[9],
      main_effects[8], main_effects[10], main_effects[10], main_effects[10],
      main_effects[10], main_effects[10], main_effects[4], main_effects[6],
      main_effects[11], main_effects[5], main_effects[7], main_effects[4],
      main_effects[6], main_effects[5], main_effects[7], main_effects[11],
      main_effects[4], main_effects[6], main_effects[11], main_effects[5],
      main_effects[7], main_effects[4], main_effects[6], main_effects[5],
      main_effects[7], main_effects[11], main_effects[4], main_effects[6],
      main_effects[11], main_effects[5], main_effects[7], main_effects[4],
      main_effects[6], main_effects[5], main_effects[7], main_effects[11],
      main_effects[3], main_effects[3], main_effects[3], main_effects[2],
      main_effects[2], main_effects[2], main_effects[9], main_effects[9],
      main_effects[9], main_effects[8], main_effects[8], main_effects[8])
    
    ext_terms_1_def <- c(main_defined[9], main_defined[3], main_defined[3],
      main_defined[2], main_defined[2], main_defined[3], main_defined[10],
      main_defined[10], main_defined[10], main_defined[10], main_defined[12],
      main_defined[14], main_defined[16], main_defined[12], main_defined[14],
      main_defined[16], main_defined[13], main_defined[15], main_defined[17],
      main_defined[13], main_defined[15], main_defined[17], main_defined[12],
      main_defined[12], main_defined[14], main_defined[13], main_defined[13],
      main_defined[15], main_defined[12], main_defined[13], main_defined[12],
      main_defined[13], main_defined[14], main_defined[15], main_defined[4],
      main_defined[6], main_defined[3], main_defined[3], main_defined[5],
      main_defined[2], main_defined[2], main_defined[4], main_defined[3],
      main_defined[3], main_defined[5], main_defined[2], main_defined[2],
      main_defined[4], main_defined[11], main_defined[11], main_defined[11],
      main_defined[11], main_defined[11], main_defined[11], main_defined[11],
      main_defined[11], main_defined[4], main_defined[6], main_defined[5],
      main_defined[4], main_defined[5], main_defined[7], main_defined[6],
      main_defined[7], main_defined[4], main_defined[6], main_defined[11],
      main_defined[5], main_defined[7], main_defined[12], main_defined[12],
      main_defined[12], main_defined[13], main_defined[13], main_defined[13],
      main_defined[13], main_defined[12], main_defined[12], main_defined[13],
      main_defined[14], main_defined[14], main_defined[14], main_defined[15],
      main_defined[15], main_defined[15], main_defined[15], main_defined[14],
      main_defined[14], main_defined[15], main_defined[16], main_defined[16],
      main_defined[16], main_defined[17], main_defined[17], main_defined[17],
      main_defined[17], main_defined[16], main_defined[16], main_defined[17],
      main_defined[12], main_defined[14], main_defined[16], main_defined[13],
      main_defined[15], main_defined[17], main_defined[12], main_defined[14],
      main_defined[16], main_defined[13], main_defined[15], main_defined[17])
    
    ext_terms_2_def <- c(main_defined[8], main_defined[2], main_defined[9],
      main_defined[8], main_defined[9], main_defined[8], main_defined[3],
      main_defined[2], main_defined[9], main_defined[8], main_defined[2],
      main_defined[2], main_defined[2], main_defined[8], main_defined[8],
      main_defined[8], main_defined[3], main_defined[3], main_defined[3],
      main_defined[9], main_defined[9], main_defined[9], main_defined[14],
      main_defined[16], main_defined[16], main_defined[15], main_defined[17],
      main_defined[17], main_defined[15], main_defined[14], main_defined[17],
      main_defined[16], main_defined[17], main_defined[16], main_defined[5],
      main_defined[7], main_defined[5], main_defined[7], main_defined[7],
      main_defined[4], main_defined[6], main_defined[6], main_defined[4],
      main_defined[6], main_defined[6], main_defined[5], main_defined[7],
      main_defined[7], main_defined[2], main_defined[4], main_defined[6],
      main_defined[3], main_defined[5], main_defined[7], main_defined[8],
      main_defined[9], main_defined[8], main_defined[8], main_defined[9],
      main_defined[9], main_defined[8], main_defined[9], main_defined[9],
      main_defined[8], main_defined[10], main_defined[10], main_defined[10],
      main_defined[10], main_defined[10], main_defined[4], main_defined[6],
      main_defined[11], main_defined[5], main_defined[7], main_defined[4],
      main_defined[6], main_defined[5], main_defined[7], main_defined[11],
      main_defined[4], main_defined[6], main_defined[11], main_defined[5],
      main_defined[7], main_defined[4], main_defined[6], main_defined[5],
      main_defined[7], main_defined[11], main_defined[4], main_defined[6],
      main_defined[11], main_defined[5], main_defined[7], main_defined[4],
      main_defined[6], main_defined[5], main_defined[7], main_defined[11],
      main_defined[3], main_defined[3], main_defined[3], main_defined[2],
      main_defined[2], main_defined[2], main_defined[9], main_defined[9],
      main_defined[9], main_defined[8], main_defined[8], main_defined[8])
    
    main1 <- c(main1, extended_terms_1)
    main2 <- c(main2, extended_terms_2)
    main1_def <- c(main1_def, ext_terms_1_def)
    main2_def <- c(main2_def, ext_terms_2_def)
  }
  basic_double_vec <- rep(0.0, basic_length)
  year_double_vec <- rep(0.0, num_years)
  patch_double_vec <- rep(0.0, num_patches)
  group_double_vec <- rep(0.0, num_groups)
  
  out_data <- if (!interactions) {
    cbind.data.frame(main_effect_1 = main1, main_1_defined = main1_def,
      surv = basic_double_vec, obs = basic_double_vec, sizea = basic_double_vec,
      sizeb = basic_double_vec, sizec = basic_double_vec, repst = basic_double_vec,
      fec = basic_double_vec, jsurv = basic_double_vec, jobs = basic_double_vec,
      jsizea = basic_double_vec, jsizeb = basic_double_vec, jsizec = basic_double_vec,
      jrepst = basic_double_vec, jmatst = basic_double_vec)
  } else {
    cbind.data.frame(main_effect_1 = main1, main_1_defined = main1_def,
      main_effect_2 = main2, main_2_defined = main2_def, surv = basic_double_vec,
      obs = basic_double_vec, sizea = basic_double_vec, sizeb = basic_double_vec,
      sizec = basic_double_vec, repst = basic_double_vec, fec = basic_double_vec,
      jsurv = basic_double_vec, jobs = basic_double_vec, jsizea = basic_double_vec,
      jsizeb = basic_double_vec, jsizec = basic_double_vec,
      jrepst = basic_double_vec, jmatst = basic_double_vec)
  }
  
  term_names <- c("surv", "obs", "sizea", "sizeb", "sizec", "repst", "fec",
    "jsurv", "jobs", "jsizea", "jsizeb", "jsizec", "jrepst", "jmatst")
  
  year_frame <- cbind.data.frame(years, year_double_vec, year_double_vec,
    year_double_vec, year_double_vec, year_double_vec, year_double_vec,
    year_double_vec, year_double_vec, year_double_vec, year_double_vec,
    year_double_vec, year_double_vec, year_double_vec, year_double_vec)
  patch_frame <- cbind.data.frame(patches, patch_double_vec, patch_double_vec,
    patch_double_vec, patch_double_vec, patch_double_vec, patch_double_vec,
    patch_double_vec, patch_double_vec, patch_double_vec, patch_double_vec,
    patch_double_vec, patch_double_vec, patch_double_vec, patch_double_vec)
  group2_frame <- cbind.data.frame(groups, group_double_vec, group_double_vec,
    group_double_vec, group_double_vec, group_double_vec, group_double_vec,
    group_double_vec, group_double_vec, group_double_vec, group_double_vec,
    group_double_vec, group_double_vec, group_double_vec, group_double_vec)
  
  names(year_frame)[2:15] <- term_names
  names(patch_frame)[2:15] <- term_names
  names(group2_frame)[2:15] <- term_names
  
  if (zi) {
    out_data <- cbind.data.frame(out_data, sizea_zi = basic_double_vec,
      sizeb_zi = basic_double_vec, sizec_zi = basic_double_vec,
      fec_zi = basic_double_vec, jsizea_zi = basic_double_vec,
      jsizeb_zi = basic_double_vec, jsizec_zi = basic_double_vec)
    year_frame <- cbind.data.frame(year_frame, sizea_zi = year_double_vec,
      sizeb_zi = year_double_vec, sizec_zi = year_double_vec,
      fec_zi = year_double_vec, jsizea_zi = year_double_vec,
      jsizeb_zi = year_double_vec, jsizec_zi = year_double_vec)
    patch_frame <- cbind.data.frame(patch_frame, sizea_zi = patch_double_vec,
      sizeb_zi = patch_double_vec, sizec_zi = patch_double_vec,
      fec_zi = patch_double_vec, jsizea_zi = patch_double_vec,
      jsizeb_zi = patch_double_vec, jsizec_zi = patch_double_vec)
    group2_frame <- cbind.data.frame(group2_frame, sizea_zi = group_double_vec,
      sizeb_zi = group_double_vec, sizec_zi = group_double_vec,
      fec_zi = group_double_vec, jsizea_zi = group_double_vec,
      jsizeb_zi = group_double_vec, jsizec_zi = group_double_vec)
  }
  
  dist_frame <- cbind.data.frame(response = term_names,
    dist = c("binom", "constant", dist.sizea, dist.sizeb, dist.sizec, "constant",
      dist.fec, "binom", "constant", dist.sizea, dist.sizeb, dist.sizec,
      "constant", "constant"))
  if (!use.juv) {
    dist_frame$dist[8] = "constant"
    dist_frame$dist[10] = "constant"
    dist_frame$dist[11] = "constant"
    dist_frame$dist[12] = "constant"
  }
  st_frame <- rep(1.0, 14)
  names(st_frame) <- term_names
  
  out_list <- list(vrm_frame = out_data, year_frame = year_frame,
    patch_frame = patch_frame, group2_frame = group2_frame,
    group1_frame = group2_frame, dist_frame = dist_frame, st_frame = st_frame)
  
  if (!is.null(cat.indcova)) {
    ica_length <- length(cat.indcova)
    ica_double_vec <- rep(0.0, ica_length)
    ica_frame <- cbind.data.frame(cat.indcova, ica_double_vec, ica_double_vec,
      ica_double_vec, ica_double_vec, ica_double_vec, ica_double_vec,
      ica_double_vec, ica_double_vec, ica_double_vec, ica_double_vec,
      ica_double_vec, ica_double_vec, ica_double_vec, ica_double_vec)
    names(ica_frame)[1] <- "indcova"
    names(ica_frame)[2:15] <- term_names
    
    if (zi) {
      ica_frame <- cbind.data.frame(ica_frame, sizea_zi = ica_double_vec,
      sizeb_zi = ica_double_vec, sizec_zi = ica_double_vec,
      fec_zi = ica_double_vec, jsizea_zi = ica_double_vec,
      jsizeb_zi = ica_double_vec, jsizec_zi = ica_double_vec)
    }
    
    out_list$indcova2_frame <- ica_frame
    out_list$indcova1_frame <- ica_frame
  }
  
  if (!is.null(cat.indcovb)) {
    icb_length <- length(cat.indcovb)
    icb_double_vec <- rep(0.0, icb_length)
    icb_frame <- cbind.data.frame(cat.indcovb, icb_double_vec, icb_double_vec,
      icb_double_vec, icb_double_vec, icb_double_vec, icb_double_vec,
      icb_double_vec, icb_double_vec, icb_double_vec, icb_double_vec,
      icb_double_vec, icb_double_vec, icb_double_vec, icb_double_vec)
    names(icb_frame)[1] <- "indcovb"
    names(icb_frame)[2:15] <- term_names
    
    if (zi) {
      icb_frame <- cbind.data.frame(icb_frame, sizea_zi = icb_double_vec,
      sizeb_zi = icb_double_vec, sizec_zi = icb_double_vec,
      fec_zi = icb_double_vec, jsizea_zi = icb_double_vec,
      jsizeb_zi = icb_double_vec, jsizec_zi = icb_double_vec)
    }
    
    out_list$indcovb2_frame <- icb_frame
    out_list$indcovb1_frame <- icb_frame
  }
  
  if (!is.null(cat.indcovc)) {
    icc_length <- length(cat.indcovc)
    icc_double_vec <- rep(0.0, icc_length)
    icc_frame <- cbind.data.frame(cat.indcovc, icc_double_vec, icc_double_vec,
      icc_double_vec, icc_double_vec, icc_double_vec, icc_double_vec,
      icc_double_vec, icc_double_vec, icc_double_vec, icc_double_vec,
      icc_double_vec, icc_double_vec, icc_double_vec, icc_double_vec)
    names(icc_frame)[1] <- "indcovc"
    names(icc_frame)[2:15] <- term_names
    
    if (zi) {
      icc_frame <- cbind.data.frame(icc_frame, sizea_zi = icc_double_vec,
      sizeb_zi = icc_double_vec, sizec_zi = icc_double_vec,
      fec_zi = icc_double_vec, jsizea_zi = icc_double_vec,
      jsizeb_zi = icc_double_vec, jsizec_zi = icc_double_vec)
    }
    
    out_list$indcovc2_frame <- icc_frame
    out_list$indcovc1_frame <- icc_frame
  }
  
  class(out_list) <- "vrm_input"
  
  return(out_list)
}

