#' Create Function-based Historical Matrix Projection Model
#' 
#' Function \code{flefko3()} returns function-based historical MPMs
#' corresponding to the patches and occasion times given, including the
#' associated component transition and fecundity matrices, data frames detailing
#' the characteristics of the ahistorical stages used and historical stage pairs
#' created, and a data frame characterizing the patch and occasion time
#' combinations corresponding to these matrices.
#' 
#' @param year A variable corresponding to the observation occasion, or a set
#' of such values, given in values associated with the year term used in linear 
#' model development. Defaults to \code{"all"}, in which case matrices will be
#' estimated for all occasion times.
#' @param patch A variable designating which patches or subpopulations will have
#' matrices estimated. Defaults to \code{"all"}, but can also be set to specific
#' patch names.
#' @param stageframe A stageframe object that includes information on the size,
#' observation status, propagule status, immaturity status, and maturity status
#' of each ahistorical stage. Should also incorporate bin widths if size is
#' continuous.
#' @param supplement An optional data frame of class \code{lefkoSD} that
#' provides supplemental data that should be incorporated into the MPM. Three
#' kinds of data may be integrated this way: transitions to be estimated via the
#' use of proxy transitions, transition overwrites from the literature or
#' supplemental studies, and transition multipliers for survival and fecundity.
#' This data frame should be produced using the \code{\link{supplemental}()}
#' function. Can be used in place of or in addition to an overwrite table (see 
#' \code{overwrite} below) and a reproduction matrix (see \code{repmatrix}
#' below).
#' @param repmatrix An optional reproduction matrix. This matrix is composed
#' mostly of 0s, with non-zero entries acting as element identifiers and
#' multipliers for fecundity (with 1 equaling full fecundity). If left blank,
#' and no \code{supplement} is provided, then \code{flefko3()} will assume that
#' all stages marked as reproductive produce offspring at 1x that of estimated
#' fecundity, and that offspring production will yield the first stage noted as
#' propagule or immature.  To prevent this behavior, input just \code{0}, which
#' will result in fecundity being estimated only for transitions noted in
#' \code{supplement} above. May be the dimensions of either a historical or an
#' ahistorical matrix. If the latter, then all stages will be used in occasion
#' \emph{t}-1 for each suggested ahistorical transition.
#' @param overwrite An optional data frame developed with the
#' \code{\link{overwrite}()} function describing transitions to be overwritten
#' either with given values or with other estimated transitions. Note that this
#' function supplements overwrite data provided in \code{supplement}.
#' @param data The historical vertical demographic data frame used to estimate
#' vital rates (class \code{hfvdata}), which is required to initialize times and
#' patches properly.
#' @param modelsuite An optional \code{lefkoMod} object holding the vital rate
#' models. If given, then \code{surv_model}, \code{obs_model}, 
#' \code{size_model}, \code{sizeb_model}, \code{sizec_model},
#' \code{repst_model}, \code{fec_model}, \code{jsurv_model}, \code{jobs_model},
#' \code{jsize_model}, \code{jsizeb_model}, \code{jsizec_model},
#' \code{jrepst_model}, \code{jmatst_model}, \code{paramnames}, \code{yearcol},
#' and \code{patchcol} are not required. One or more of these  models should
#' include size or reproductive status in occasion \emph{t}-1.
#' @param surv_model A linear model predicting survival probability. This can 
#' be a model of class \code{glm} or \code{glmer}, and requires a predicted
#' binomial variable under a logit link. If given, then will overwrite any
#' survival probability model given in \code{modelsuite}. This model must have
#' been developed in a modeling exercise testing the impacts of occasions
#' \emph{t} and \emph{t}-1.
#' @param obs_model A linear model predicting sprouting or observation
#' probability. This can be a model of class \code{glm} or \code{glmer}, and
#' requires a predicted binomial variable under a logit link. If given, then 
#' will overwrite any observation probability model given in \code{modelsuite}.
#' This model must have been developed in a modeling exercise testing the
#' impacts of occasions \emph{t} and \emph{t}-1.
#' @param size_model A linear model predicting primary size. This can be a model
#' of class \code{glm}, \code{glmer}, \code{glmmTMB}, \code{zeroinfl},
#' \code{vglm}, \code{lm}, or \code{lmer}. If given, then will overwrite any
#' primary size model given in \code{modelsuite}. This model must have been
#' developed in a modeling exercise testing the impacts of occasions \emph{t}
#' and \emph{t}-1.
#' @param sizeb_model A linear model predicting secondary size. This can be a
#' modelvof class \code{glm}, \code{glmer}, \code{glmmTMB}, \code{zeroinfl},
#' \code{vglm}, \code{lm}, or \code{lmer}. If given, then will overwrite any
#' secondary size model given in \code{modelsuite}. This model must have been
#' developed in a modeling exercise testing the impacts of occasions \emph{t}
#' and \emph{t}-1.
#' @param sizec_model A linear model predicting tertiary size. This can be a
#' modelvof class \code{glm}, \code{glmer}, \code{glmmTMB}, \code{zeroinfl},
#' \code{vglm}, \code{lm}, or \code{lmer}. If given, then will overwrite any
#' tertiary size model given in \code{modelsuite}. This model must have been
#' developed in a modeling exercise testing the impacts of occasions \emph{t}
#' and \emph{t}-1.
#' @param repst_model A linear model predicting reproduction probability. This 
#' can be a model of class \code{glm} or \code{glmer}, and requires a predicted
#' binomial variable under a logit link. If given, then will overwrite any
#' reproduction probability model given in \code{modelsuite}. This model must
#' have been developed in a modeling exercise testing the impacts of occasions
#' \emph{t} and \emph{t}-1.
#' @param fec_model A linear model predicting fecundity. This can be a model of
#' class \code{glm}, \code{glmer}, \code{glmmTMB}, \code{zeroinfl},
#' \code{vglm}, \code{lm}, or \code{lmer}. If given, then will overwrite any 
#' fecundity model given in \code{modelsuite}. This model must have been
#' developed in a modeling exercise testing the impacts of occasions \emph{t}
#' and \emph{t}-1.
#' @param jsurv_model A linear model predicting juvenile survival probability.
#' This can be a model of class \code{glm} or \code{glmer}, and requires a
#' predicted binomial variable under a logit link. If given, then will overwrite
#' any juvenile survival probability model given in \code{modelsuite}. This
#' model must have been developed in a modeling exercise testing the impacts of
#' occasions \emph{t} and \emph{t}-1.
#' @param jobs_model A linear model predicting juvenile sprouting or observation
#' probability. This can be a model of class \code{glm} or \code{glmer}, and
#' requires a predicted binomial variable under a logit link. If given, then
#' will overwrite any juvenile observation probability model given in 
#' \code{modelsuite}. This model must have been developed in a modeling exercise
#' testing the impacts of occasions \emph{t} and \emph{t}-1.
#' @param jsize_model A linear model predicting juvenile primary size. This
#' can be a model of class \code{glm}, \code{glmer}, \code{glmmTMB},
#' \code{zeroinfl}, \code{vglm}, \code{lm}, or \code{lmer}. If given, then will
#' overwrite any juvenile primary size model given in \code{modelsuite}. This
#' model must have been developed in a modeling exercise testing the impacts of 
#' occasions \emph{t} and \emph{t}-1.
#' @param jsizeb_model A linear model predicting juvenile secondary size. This
#' can be a model of class \code{glm}, \code{glmer}, \code{glmmTMB},
#' \code{zeroinfl}, \code{vglm}, \code{lm}, or \code{lmer}. If given, then will
#' overwrite any juvenile secondary size model given in \code{modelsuite}. This
#' model must have been developed in a modeling exercise testing the impacts of 
#' occasions \emph{t} and \emph{t}-1.
#' @param jsizec_model A linear model predicting juvenile tertiary size. This
#' can be a model of class \code{glm}, \code{glmer}, \code{glmmTMB},
#' \code{zeroinfl}, \code{vglm}, \code{lm}, or \code{lmer}. If given, then will
#' overwrite any juvenile tertiary size model given in \code{modelsuite}. This
#' model must have been developed in a modeling exercise testing the impacts of 
#' occasions \emph{t} and \emph{t}-1.
#' @param jrepst_model A linear model predicting reproduction probability of a 
#' mature individual that was immature in time \emph{t}. This can be a model
#' of class \code{glm} or \code{glmer}, and requires a predicted binomial
#' variable under a logit link. If given, then will overwrite any reproduction
#' probability model given in \code{modelsuite}. This model must have been
#' developed in a modeling exercise testing the impacts of occasions \emph{t}
#' and \emph{t}-1.
#' @param jmatst_model A linear model predicting maturity probability of an 
#' individual that was immature in time \emph{t}. This can be a model of class
#' \code{glm} or \code{glmer}, and requires a predicted binomial variable under
#' a logit link. If given, then will overwrite any maturity probability model
#' given in \code{modelsuite}. This model must have been developed in a modeling
#' exercise testing the impacts of occasions \emph{t} and \emph{t}-1.
#' @param paramnames A dataframe with three columns, the first describing all
#' terms used in linear modeling, the second (must be called \code{mainparams}),
#' showing the general model terms that will be used in matrix creation (users
#' should use \code{\link{modelsearch}()} at least once to see the proper
#' names to be used in this column), and the third showing the equivalent terms
#' used in modeling (must be named \code{modelparams}). Only required if
#' \code{modelsuite} is not supplied.
#' @param inda Can be a single value to use for individual covariate \code{a}
#' in all matrices, a pair of values to use for times \emph{t} and \emph{t}-1 in
#' historical matrices, or a vector of such values corresponding to each
#' occasion in option \code{year}. Defaults to \code{NULL}.
#' @param indb Can be a single value to use for individual covariate \code{b}
#' in all matrices, a pair of values to use for times \emph{t} and \emph{t}-1 in
#' historical matrices, or a vector of such values corresponding to each
#' occasion in option \code{year}. Defaults to \code{NULL}.
#' @param indc Can be a single value to use for individual covariate \code{c}
#' in all matrices, a pair of values to use for times \emph{t} and \emph{t}-1 in
#' historical matrices, or a vector of such values corresponding to each
#' occasion in option \code{year}. Defaults to \code{NULL}.
#' @param surv_dev A numeric value to be added to the y-intercept in the linear
#' model for survival probability.
#' @param obs_dev A numeric value to be added to the y-intercept in the linear
#' model for observation probability.
#' @param size_dev A numeric value to be added to the y-intercept in the linear
#' model for primary size.
#' @param sizeb_dev A numeric value to be added to the y-intercept in the linear
#' model for secondary size.
#' @param sizec_dev A numeric value to be added to the y-intercept in the linear
#' model for tertiary size.
#' @param repst_dev A numeric value to be added to the y-intercept in the linear
#' model for probability of reproduction.
#' @param fec_dev A numeric value to be added to the y-intercept in the linear
#' model for fecundity.
#' @param jsurv_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile survival probability.
#' @param jobs_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile observation probability.
#' @param jsize_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile primary size.
#' @param jsizeb_dev A numeric value to be added to the y-intercept in the
#' linear model for juvenile secondary size.
#' @param jsizec_dev A numeric value to be added to the y-intercept in the
#' linear model for juvenile tertiary size.
#' @param jrepst_dev A numeric value to be added to the y-intercept in the
#' linear model for juvenile reproduction probability.
#' @param jmatst_dev A numeric value to be added to the y-intercept in the
#' linear model for juvenile maturity probability.
#' @param density A numeric value indicating density value to use to propagate
#' matrices. Only needed if density is an explanatory term used in linear
#' models. Defaults to \code{NA}.
#' @param repmod A scalar multiplier of fecundity. Defaults to \code{1}.
#' @param yearcol The variable name or column number corresponding to year 
#' in occasion \emph{t} in the dataset. Not needed if \code{modelsuite} is
#' supplied.
#' @param patchcol The variable name or column number corresponding to patch in 
#' the dataset. Not needed if \code{modelsuite} is supplied.
#' @param year.as.random A logical term indicating whether coefficients for
#' missing occasions within vital rate models should be estimated as random
#' intercepts. Defaults to \code{FALSE}, in which case missing monitoring
#' occasion coefficients are set to \code{0}.
#' @param patch.as.random A logical term indicating whether coefficients for
#' missing patches within vital rate models should be estimated as random
#' intercepts. Defaults to \code{FALSE}, in which case missing patch
#' coefficients are set to \code{0}.
#' @param random.inda A logical value denoting whether to treat individual
#' covariate \code{a} as a random, categorical variable. Otherwise is treated as
#' a fixed, numeric variable. Defaults to \code{FALSE}.
#' @param random.indb A logical value denoting whether to treat individual
#' covariate \code{b} as a random, categorical variable. Otherwise is treated as
#' a fixed, numeric variable. Defaults to \code{FALSE}.
#' @param random.indc A logical value denoting whether to treat individual
#' covariate \code{c} as a random, categorical variable. Otherwise is treated as
#' a fixed, numeric variable. Defaults to \code{FALSE}.
#' @param randomseed A numeric value used as a seed to generate random estimates
#' for missing occasion and patch coefficients, if either \code{year.as.random}
#' or \code{patch.as.random} is set to \code{TRUE}. Defaults to 
#' \code{\link{set.seed}()} default.
#' @param negfec A logical value denoting whether fecundity values estimated to
#' be negative should be reset to \code{0}. Defaults to \code{FALSE}.
#' @param format A string indicating whether to estimate matrices in
#' \code{ehrlen} format or \code{deVries} format. The latter adds one extra
#' prior stage to account for the prior state of newborns. Defaults to
#' \code{ehrlen} format.
#' @param ipm_method A string indicating what method to use to estimate size
#' transition probabilities if size is treated as continuous. Options include
#' \code{"midpoint"}, which utilizes the midpoint method, and \code{"CDF"},
#' which uses the cumulative distribution function. Defaults to \code{"CDF"}.
#' @param reduce A logical value denoting whether to remove historical stages
#' associated solely with 0 transitions. These are only removed in cases where
#' the associated row and column sums in ALL matrices estimated equal 0. 
#' Defaults to \code{FALSE}.
#' @param err_check A logical value indicating whether to append matrices of
#' vital rate probabilities associated with each matrix to the \code{lefkoMat}
#' object generated. These matrices are developed internally and can be used for
#' error checking. Defaults to \code{FALSE}.
#' @param exp_tol A numeric value used to indicate a maximum value to set
#' exponents to in the core kernel to prevent numerical overflow. Defaults to
#' \code{700}.
#' @param theta_tol A numeric value used to indicate a maximum value to theta as
#' used in the negative binomial probability density kernel. Defaults to
#' \code{100000000}, but can be reset to other values during error checking.
#'
#' @return If all inputs are properly formatted, then this function will return
#' an object of class \code{lefkoMat}, which is a list that holds the matrix
#' projection model and all of its metadata. Its structure is a list with the
#' following elements:
#' 
#' \item{A}{A list of full projection matrices in order of sorted patches and
#' occasion times. All matrices output in R's \code{matrix} class.}
#' \item{U}{A list of survival transition matrices sorted as in \code{A}. All 
#' matrices output in R's \code{matrix} class.}
#' \item{F}{A list of fecundity matrices sorted as in \code{A}. All matrices 
#' output in R's \code{matrix} class.}
#' \item{hstages}{A data frame matrix showing the pairing of ahistorical stages
#' used to create historical stage pairs.}
#' \item{agestages}{A data frame showing age-stage pairs. In this function, it
#' is set to \code{NA}. Only used in output to function \code{aflefko2}().}
#' \item{ahstages}{A data frame detailing the characteristics of associated
#' ahistorical stages, in the form of a modified stageframe that includes
#' status as an entry stage through reproduction.}
#' \item{labels}{A data frame giving the population, patch, and year of each
#' matrix in order. In \code{flefko3()}, only one population may be analyzed at
#' once, and so \code{pop = NA}.}
#' \item{matrixqc}{A short vector describing the number of non-zero elements in
#' \code{U} and \code{F} matrices, and the number of annual matrices.}
#' \item{modelqc}{This is the \code{qc} portion of the \code{modelsuite} input.}
#' \item{prob_out}{An optional element only added if \code{err_check = TRUE}.
#' This is a list of vital rate probability matrices, with 6 columns in the
#' order of survival, observation probability, reproduction probability, primary
#' size transition probability, secondary size transition probability, and
#' tertiary size transition probability.}
#'
#' @section Notes:
#' Unlike \code{\link{rlefko2}()}, \code{\link{rlefko3}()}, and
#' \code{\link{rleslie}()}, this function does not currently distinguish
#' populations.
#' 
#' The default behavior of this function is to estimate fecundity with regards
#' to transitions specified via associated fecundity multipliers in either
#' \code{supplement} or \code{repmatrix}. If both of these fields are left
#' empty, then fecundity will be estimated at full for all transitions leading
#' from reproductive stages to immature and propagule stages. However, if a
#' \code{supplement} is provided and a \code{repmatrix} is not, or if
#' \code{repmatrix} is set to \code{0}, then only fecundity transitions noted in
#' the \code{supplement} will be set to non-zero values. To use the default
#' behavior of setting all reproductive stages to reproduce at full fecundity
#' into immature and propagule stages, but also incorporate given or proxy
#' survival transitions, input those given and proxy transitions through the
#' \code{overwrite} option.
#' 
#' The reproduction matrix (field \code{repmatrix}) may be supplied as either
#' historical or ahistorical. If provided as ahistorical, then \code{flefko3()}
#' will assume that all historical transitions involving stages noted for
#' occasions \emph{t} and \emph{t}+1 should be set to the respective fecundity
#' multipliers noted.
#' 
#' Users may at times wish to estimate MPMs using a dataset incorporating
#' multiple patches or subpopulations, but without discriminating between those
#' patches or subpopulations. Should the aim of analysis be a general MPM that
#' does not distinguish these patches or subpopulations, the \code{patchcol}
#' variable should be set to \code{NA}, which is the default.
#'
#' Input options including multiple variable names must be entered in the order
#' of variables in occasion \emph{t}+1, \emph{t}, and \emph{t}-1. Rearranging
#' the order will lead to erroneous calculations, and will may lead to fatal
#' errors.
#' 
#' Care should be taken to match the random status of year and patch to the
#' states of those variables within the modelsuite. If they do not match, then
#' they will be treated as zeroes in vital rate estimation.
#' 
#' The \code{ipm_method} function gives the option of using two different means
#' of estimating the probability of size transition. The midpoint method
#' (\code{"midpoint"}) refers to the method in which the probability is
#' estimated by first estimating the probability associated with transition from
#' the exact size at the midpoint of the size class using the corresponding
#' probability density function, and then multiplying that value by the bin
#' width of the size class. Doak et al. 2021 (Ecological Monographs) noted that
#' this method can produce biased results, with total size transitions
#' associated with a specific size not totaling to 1.0 and even specific size
#' transition probabilities capable of being estimated at values greater than
#' 1.0. The alternative and default method, \code{"CDF"}, uses the corresponding
#' cumulative density function to estimate the probability of size transition as
#' the cumulative probability of size transition at the greater limit of the
#' size class minus the cumulative probability of size transition at the lower
#' limit of the size class. The latter method avoids this bias. Note, however,
#' that both methods are exact and unbiased for negative binomial and Poisson
#' distributions.
#' 
#' Under the Gaussian and gamma size distributions, the number of estimated
#' parameters may differ between the two \code{ipm_method} settings. Because
#' the midpoint method has a tendency to incorporate upward bias in the
#' estimation of size transition probabilities, it is more likely to yield non-
#' zero values when the true probability is extremely close to 0. This will
#' result in the \code{summary.lefkoMat} function yielding higher numbers of
#' estimated parameters than the \code{ipm_method = "CDF"} yields in some cases.
#' 
#' Using the \code{err_check} option will produce a matrix of 6 columns, each
#' characterizing a different vital rate. The product of each row yields an
#' element in the associated \code{$U} matrix. The number and order of elements
#' in each column of this matrix matches the associated matrix in column vector
#' format. Use of this option is generally for the purposes of debugging code.
#'`
#' @seealso \code{\link{flefko2}()}
#' @seealso \code{\link{aflefko2}()}
#' @seealso \code{\link{fleslie}()}
#' @seealso \code{\link{rlefko3}()}
#' @seealso \code{\link{rlefko2}()}
#' @seealso \code{\link{rleslie}()}
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
#' 
#' #Cypripedium example using three size metrics for classification
#' data(cypdata)
#' sizevector.f <- c(0, 0, 0, 0, 0, 0, seq(1, 12, by = 1), seq(0, 9, by = 1),
#'   seq(0, 8, by = 1), seq(0, 7, by = 1), seq(0, 6, by = 1), seq(0, 5, by = 1),
#'   seq(0, 4, by = 1), seq(0, 3, by = 1), 0, 1, 2, 0, 1, 0, 
#'   0, 0, 1, 0)
#' sizebvector.f <- c(0, 0, 0, 0, 0, 0, rep(0, 12), rep(1, 10), rep(2, 9),
#'   rep(3, 8), rep(4, 7), rep(5, 6), rep(6, 5), rep(7, 4), rep(8, 3), 9, 9, 10, 
#'   0, 1, 1, 2)
#' sizecvector.f <- c(0, 0, 0, 0, 0, 0, rep(0, 12), rep(0, 10), rep(0, 9),
#'   rep(0, 8), rep(0, 7), rep(0, 6), rep(0, 5), rep(0, 4), 0, 0, 0, 0, 0, 0, 
#'   1, 1, 1, 1)
#' stagevector.f <- c("DS", "P1", "P2", "P3", "Sdl", "Dorm", "V1 I0 D0",
#'   "V2 I0 D0", "V3 I0 D0", "V4 I0 D0", "V5 I0 D0", "V6 I0 D0", "V7 I0 D0",
#'   "V8 I0 D0", "V9 I0 D0", "V10 I0 D0", "V11 I0 D0", "V12 I0 D0", "V0 I1 D0",
#'   "V1 I1 D0", "V2 I1 D0", "V3 I1 D0", "V4 I1 D0", "V5 I1 D0", "V6 I1 D0",
#'   "V7 I1 D0", "V8 I1 D0", "V9 I1 D0", "V0 I2 D0", "V1 I2 D0", "V2 I2 D0",
#'   "V3 I2 D0", "V4 I2 D0", "V5 I2 D0", "V6 I2 D0", "V7 I2 D0", "V8 I2 D0",
#'   "V0 I3 D0", "V1 I3 D0", "V2 I3 D0", "V3 I3 D0", "V4 I3 D0", "V5 I3 D0",
#'   "V6 I3 D0", "V7 I3 D0", "V0 I4 D0", "V1 I4 D0", "V2 I4 D0", "V3 I4 D0",
#'   "V4 I4 D0", "V5 I4 D0", "V6 I4 D0", "V0 I5 D0", "V1 I5 D0", "V2 I5 D0",
#'   "V3 I5 D0", "V4 I5 D0", "V5 I5 D0", "V0 I6 D0", "V1 I6 D0", "V2 I6 D0",
#'   "V3 I6 D0", "V4 I6 D0", "V0 I7 D0", "V1 I7 D0", "V2 I7 D0", "V3 I7 D0",
#'   "V0 I8 D0", "V1 I8 D0", "V2 I8 D0", "V0 I9 D0", "V1 I9 D0", "V0 I10 D0",
#'   "V0 I0 D1", "V0 I1 D1", "V1 I1 D1", "V0 I2 D1")
#' repvector.f <- c(0, 0, 0, 0, 0, rep(0, 13), rep(1, 59))
#' obsvector.f <- c(0, 0, 0, 0, 0, 0, rep(1, 71))
#' matvector.f <- c(0, 0, 0, 0, 0, rep(1, 72))
#' immvector.f <- c(0, 1, 1, 1, 1, rep(0, 72))
#' propvector.f <- c(1, rep(0, 76))
#' indataset.f <- c(0, 0, 0, 0, 0, rep(1, 72))
#' binvec.f <- c(0, 0, 0, 0, 0, rep(0.5, 72))
#' binbvec.f <- c(0, 0, 0, 0, 0, rep(0.5, 72))
#' bincvec.f <- c(0, 0, 0, 0, 0, rep(0.5, 72))
#' 
#' vertframe.f <- sf_create(sizes = sizevector.f, sizesb = sizebvector.f,
#'   sizesc = sizecvector.f, stagenames = stagevector.f, repstatus = repvector.f,
#'   obsstatus = obsvector.f, propstatus = propvector.f, immstatus = immvector.f,
#'   matstatus = matvector.f, indataset = indataset.f, binhalfwidth = binvec.f,
#'   binhalfwidthb = binbvec.f, binhalfwidthc = bincvec.f)
#' 
#' vert.data.f <- verticalize3(cypdata, noyears = 6, firstyear = 2004,
#'   individcol = "plantid", blocksize = 4, sizeacol = "Veg.04",
#'   sizebcol = "Inf.04", sizeccol = "Inf2.04", repstracol = "Inf.04",
#'   repstrbcol = "Inf2.04", fecacol = "Pod.04", censorcol = "censor",
#'   censorkeep = 1, censorRepeat = FALSE, stageassign = vertframe.f,
#'   stagesize = "sizeabc", NAas0 = TRUE, censor = FALSE)
#' 
#' vertmodels3f <- modelsearch(vert.data.f, historical = TRUE, suite = "main",
#'   sizeb = c("sizeb3", "sizeb2", "sizeb1"), sizec = c("sizec3", "sizec2", "sizec1"),
#'   approach = "glm", vitalrates = c("surv", "obs", "size", "repst", "fec"),
#'   sizedist = "negbin", sizebdist = "poisson", sizecdist = "poisson",
#'   fecdist = "poisson", patch.as.random = TRUE, year.as.random = TRUE)
#' 
#' vertsupp3f <- supplemental(stage3 = c("DS", "P1", "DS", "P1", "P2", "P2", "P3",
#'     "Sdl", "Sdl", "Sdl", "Dorm", "V1 I0 D0", "V2 I0 D0", "V3 I0 D0", "Dorm",
#'     "V1 I0 D0", "V2 I0 D0", "V3 I0 D0", "DS", "P1"),
#'   stage2 = c("DS", "DS", "DS", "DS", "P1", "P1", "P2", "P3", "Sdl", "Sdl", "Sdl",
#'     "Sdl", "Sdl", "Sdl", "Sdl", "Sdl", "Sdl", "Sdl", "rep", "rep"),
#'   stage1 = c("DS", "DS", "rep", "rep", "DS", "rep", "P1", "P2", "P3", "Sdl",
#'     "Sdl", "Sdl", "Sdl", "Sdl", "P3", "P3", "P3", "P3", "mat", "mat"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "Dorm", "V1 I0 D0",
#'     "V2 I0 D0", "V3 I0 D0", "Dorm", "V1 I0 D0", "V2 I0 D0", "V3 I0 D0", NA, NA), 
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "V1 I0 D0", "V1 I0 D0",
#'     "V1 I0 D0", "V1 I0 D0", "V1 I0 D0", "V1 I0 D0", "V1 I0 D0", "V1 I0 D0", NA, NA),
#'   eststage1 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "V1 I0 D0", "V1 I0 D0",
#'     "V1 I0 D0", "V1 I0 D0", "V1 I0 D0", "V1 I0 D0", "V1 I0 D0", "V1 I0 D0", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.10, 0.20, 0.20, 0.20, 0.20, 0.25, 0.40, 0.40, NA,
#'     NA, NA, NA, NA, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
#'     NA, NA, 0.5 * 5000, 0.5 * 5000),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   stageframe = vertframe.f, historical = TRUE)
#' 
#' vert.mats.f3 <- flefko3(stageframe = vertframe.f, supplement = vertsupp3f, 
#'   data = vert.data.f, modelsuite = vertmodels3f)
#' summary(vert.mats.f3)
#' }
#' 
#' @export
flefko3 <- function(year = "all", patch = "all", stageframe, supplement = NULL,
  repmatrix = NULL, overwrite = NULL, data = NA, modelsuite = NA,
  surv_model = NA, obs_model = NA, size_model = NA, sizeb_model = NA,
  sizec_model = NA, repst_model = NA, fec_model = NA, jsurv_model = NA,
  jobs_model = NA, jsize_model = NA, jsizeb_model = NA, jsizec_model = NA,
  jrepst_model = NA, jmatst_model = NA, paramnames = NA, inda = NULL, indb = NULL, indc = NULL,
  surv_dev = 0, obs_dev = 0, size_dev = 0, sizeb_dev = 0, sizec_dev = 0,
  repst_dev = 0, fec_dev = 0, jsurv_dev = 0, jobs_dev = 0, jsize_dev = 0,
  jsizeb_dev = 0, jsizec_dev = 0, jrepst_dev = 0, jmatst_dev = 0, density = NA, repmod = 1,
  yearcol = NA, patchcol = NA, year.as.random = FALSE, patch.as.random = FALSE,
  random.inda = FALSE, random.indb = FALSE, random.indc = FALSE, 
  randomseed = NA, negfec = FALSE, format = "ehrlen", ipm_method = "CDF",
  reduce = FALSE, err_check = FALSE, exp_tol = 700, theta_tol = 100000000) {
  
  indanames <- indbnames <- indcnames <- NULL
  
  if (tolower(format) == "ehrlen") {
    format_int <- 1
  } else if (tolower(format) == "devries") {
    format_int <- 2
  } else {
    stop("The format parameter must be set to either 'ehrlen' or 'deVries'.",
      call. = FALSE)
  }
  
  ipm_method <- tolower(ipm_method)
  if (length(grep("mi", ipm_method)) > 0) {
    ipm_method <- "midpoint"
  } else if (length(grep("cd", ipm_method)) > 0) {
    ipm_method <- "cdf"
  } else {
    stop("Option ipm_method not recognized.", call. = FALSE)
  }
  
  if (all(is.na(modelsuite)) & all(is.na(paramnames))) {
    warning("Function may not work properly without a dataframe of model parameters or
      equivalents supplied either through the modelsuite option or through the paramnames
      input parameter.", call. = FALSE)
  } else if (!all(is.na(modelsuite))) {
    paramnames <- modelsuite$paramnames
    yearcol <- paramnames$modelparams[which(paramnames$mainparams == "year2")]
    patchcol <- paramnames$modelparams[which(paramnames$mainparams == "patch")]
  }
  
  if (all(is.na(data))) {
    stop("Need original vertical dataset to set proper limits on year and patch.", 
      call. = FALSE)
  }
  if (!any(class(data) == "data.frame")) {
    stop("Need original vertical dataset used in modeling to proceed.",
      call. = FALSE)
  }
  if (!any(class(data) == "hfvdata")) {
    warning("Dataset used as input is not of class hfvdata. Will assume that the
      dataset has been formatted equivalently.", call. = FALSE)
  }
  
  stageframe_vars <- c("stage", "size", "size_b", "size_c", "min_age", "max_age",
    "repstatus", "obsstatus", "propstatus", "immstatus", "matstatus", "indataset",
    "binhalfwidth_raw", "sizebin_min", "sizebin_max", "sizebin_center",
    "sizebin_width", "binhalfwidthb_raw", "sizebinb_min", "sizebinb_max",
    "sizebinb_center", "sizebinb_width", "binhalfwidthc_raw", "sizebinc_min",
    "sizebinc_max", "sizebinc_center", "sizebinc_width", "group", "comments")
  if (any(!is.element(names(stageframe), stageframe_vars))) {
    stop("Please use properly formatted stageframe as input.", call. = FALSE)
  }
  
  if (is.character(yearcol)) {
    choicevar <- which(names(data) == yearcol);
    mainyears <- sort(unique(data[,choicevar]))
  } else if (is.numeric(yearcol)) {
    mainyears <- sort(unique(data[, yearcol]));
  } else {
    stop("Need appropriate year column designation.", call. = FALSE)
  }
  
  if (any(is.character(year))) {
    if (is.element("all", tolower(year))) {
      year <- mainyears
    } else {
      stop("Year designation not recognized.", call. = FALSE)
    }
  }
  
  if (length(year) == 0 | all(is.na(year) == TRUE) | any(is.na(year))) {
    stop("This function cannot proceed without a specific occasion, or a suite of
      occasions, designated via the year option. NA entries are not allowed.",
      call. = FALSE)
  }
  
  if (all(is.na(patch)) & !is.na(patchcol)) {
    warning("Matrix creation may not proceed properly without input in the patch
      option when using a modelsuite in which patch is designated.",
      call. = FALSE)
  }
  
  if (is.character(patchcol) & patchcol != "none") {
    choicevar <- which(names(data) == patchcol);
    mainpatches <- sort(unique(as.character(data[,choicevar])))
  } else if (is.numeric(patchcol)) {
    mainpatches <- sort(unique(as.character(data[, patchcol])));
  } else {
    mainpatches <- NA
  }
  
  if (any(is.character(patch))) {
    if (is.element("all", tolower(patch))) {
      patch <- mainpatches
    } else if (!all(is.element(patch, mainpatches))) {
      stop("Patch designation not recognized.", call. = FALSE)
    }
  }
  
  if (!is.null(inda)) {
    if (!is.numeric(inda) & !random.inda) {
      stop("Individual covariate vector a must be numeric if not set to random.",
        call. = FALSE)
    }
    
    if (!is.element(length(inda), c(1, 2, length(year)))) {
      stop("Individual covariate vector a must be empty, or include 1, 2, or as
        many elements as occasions to be modeled.", call. = FALSE)
    }
    
    if (random.inda) {
      indacol <- paramnames$modelparams[which(paramnames$mainparams == "indcova2")]
      if (indacol == "none") {
        stop("Individual covariate a not recognized in the modelsuite", call. = FALSE)
      }
      
      indacol <- which(names(data) == indacol)
      
      if (length(indacol) > 0) {
        indanames <- sort(unique(data[, indacol]))
      } else {
        stop("Individual covariate a not recognized in the modelsuite", call. = FALSE)
      }
      
      if (any(!is.element(inda, indanames))) {
        stop("Entered value for individual covariate a does not exist in the data.",
          call. = FALSE)
      }
      
      if (length(inda) == 1) {
        r1.inda <- rep(as.character(inda), length(year))
        r2.inda <- rep(as.character(inda), length(year))
      } else if (length(inda) == 2 & length(year) != 2) {
        r1.inda <- rep(as.character(inda[1]), length(year))
        r2.inda <- rep(as.character(inda[2]), length(year))
      } else if (length(inda) == length(year)) {
        r2.inda <- as.character(inda)
        r1.inda <- c("none", r2.inda[1:(length(inda) - 1)])
      }
      
      f1.inda <- rep(0, length(year))
      f2.inda <- rep(0, length(year))
      
    } else {
      indanames <- c(0)
      
      if (length(inda) == 1) {
        f1.inda <- rep(inda, length(year))
        f2.inda <- rep(inda, length(year))
      } else if (length(inda) == 2 & length(year) != 2) {
        f1.inda <- rep(inda[1], length(year))
        f2.inda <- rep(inda[2], length(year))
      } else if (length(inda) == length(year)) {
        f2.inda <- inda
        f1.inda <- c(0, f2.inda[1:(length(inda) - 1)])
      }
      r2.inda <- rep("none", length(year))
      r1.inda <- rep("none", length(year))
    }
  } else {
    indanames <- c(0)
    
    f1.inda <- rep(0, length(year))
    f2.inda <- rep(0, length(year))
    r2.inda <- rep("none", length(year))
    r1.inda <- rep("none", length(year))
  }
  
  if (!is.null(indb)) {
    if (!is.numeric(indb) & !random.indb) {
      stop("Individual covariate vector b must be numeric if not set to random.",
        call. = FALSE)
    }
    
    if (!is.element(length(indb), c(1, 2, length(year)))) {
      stop("Individual covariate vector b must be empty, or include 1, 2, or as
        many elements as occasions to be modeled.", call. = FALSE)
    }
    
    if (random.indb) {
      indbcol <- paramnames$modelparams[which(paramnames$mainparams == "indcovb2")]
      if (indbcol == "none") {
        stop("Individual covariate b not recognized in the modelsuite", call. = FALSE)
      }
      
      indbcol <- which(names(data) == indbcol)
      
      if (length(indbcol) > 0) {
        indbnames <- sort(unique(data[, indbcol]))
      } else {
        stop("Individual covariate b not recognized in the modelsuite", call. = FALSE)
      }
      
      if (any(!is.element(indb, indbnames))) {
        stop("Entered value for individual covariate b does not exist in the data.",
          call. = FALSE)
      }
      
      if (length(indb) == 1) {
        r1.indb <- rep(as.character(indb), length(year))
        r2.indb <- rep(as.character(indb), length(year))
      } else if (length(indb) == 2 & length(year) != 2) {
        r1.indb <- rep(as.character(indb[1]), length(year))
        r2.indb <- rep(as.character(indb[2]), length(year))
      } else if (length(indb) == length(year)) {
        r2.indb <- as.character(indb)
        r1.indb <- c("none", r2.indb[1:(length(indb) - 1)])
      }
      
      f1.indb <- rep(0, length(year))
      f2.indb <- rep(0, length(year))
      
    } else {
      indbnames <- c(0)
      
      if (length(indb) == 1) {
        f1.indb <- rep(indb, length(year))
        f2.indb <- rep(indb, length(year))
      } else if (length(indb) == 2 & length(year) != 2) {
        f1.indb <- rep(indb[1], length(year))
        f2.indb <- rep(indb[2], length(year))
      } else if (length(indb) == length(year)) {
        f2.indb <- indb
        f1.indb <- c(0, f2.indb[1:(length(indb) - 1)])
      }
      r2.indb <- rep("none", length(year))
      r1.indb <- rep("none", length(year))
    }
  } else {
    indbnames <- c(0)
    
    f1.indb <- rep(0, length(year))
    f2.indb <- rep(0, length(year))
    r2.indb <- rep("none", length(year))
    r1.indb <- rep("none", length(year))
  }
  
  if (!is.null(indc)) {
    if (!is.numeric(indc) & !random.indc) {
      stop("Individual covariate vector c must be numeric if not set to random.",
        call. = FALSE)
    }
    
    if (!is.element(length(indc), c(1, 2, length(year)))) {
      stop("Individual covariate vector c must be empty, or include 1, 2, or as
        many elements as occasions to be modeled.", call. = FALSE)
    }
    
    if (random.indc) {
      indccol <- paramnames$modelparams[which(paramnames$mainparams == "indcovc2")]
      if (indccol == "none") {
        stop("Individual covariate c not recognized in the modelsuite", call. = FALSE)
      }
      
      indccol <- which(names(data) == indccol)
      
      if (length(indccol) > 0) {
        indcnames <- sort(unique(data[, indccol]))
      } else {
        stop("Individual covariate c not recognized in the modelsuite", call. = FALSE)
      }
      
      if (any(!is.element(indc, indcnames))) {
        stop("Entered value for individual covariate c does not exist in the data.",
          call. = FALSE)
      }
      
      if (length(indc) == 1) {
        r1.indc <- rep(as.character(indc), length(year))
        r2.indc <- rep(as.character(indc), length(year))
      } else if (length(indc) == 2 & length(year) != 2) {
        r1.indc <- rep(as.character(indc[1]), length(year))
        r2.indc <- rep(as.character(indc[2]), length(year))
      } else if (length(indc) == length(year)) {
        r2.indc <- as.character(indc)
        r1.indc <- c("none", r2.indc[1:(length(indc) - 1)])
      }
      
      f1.indc <- rep(0, length(year))
      f2.indc <- rep(0, length(year))
      
    } else {
      indcnames <- c(0)
      
      if (length(indc) == 1) {
        f1.indc <- rep(indc, length(year))
        f2.indc <- rep(indc, length(year))
      } else if (length(indc) == 2 & length(year) != 2) {
        f1.indc <- rep(indc[1], length(year))
        f2.indc <- rep(indc[2], length(year))
      } else if (length(indc) == length(year)) {
        f2.indc <- indc
        f1.indc <- c(0, f2.indc[1:(length(indc) - 1)])
      }
      r2.indc <- rep("none", length(year))
      r1.indc <- rep("none", length(year))
    }
  } else {
    indcnames <- c(0)
    
    f1.indc <- rep(0, length(year))
    f2.indc <- rep(0, length(year))
    r2.indc <- rep("none", length(year))
    r1.indc <- rep("none", length(year))
  }
  
  maingroups <- sort(unique(stageframe$group))
  
  if (!all(is.na(density))) {
    if (!all(is.numeric(density))) {
      stop("Density value must be numeric.", call. = FALSE)
    }
    
    if (any(is.na(density))) {
      density[which(is.na(density))] <- 0
    }
  } else {
    density <- 0
  }
  
  if (all(is.na(repmatrix)) & all(is.na(supplement))) {
    warning("Neither supplemental data nor a reproduction matrix have been supplied.
      All fecundity transitions will be inferred from the stageframe.",
      call. = FALSE)
  } else if (all(is.na(repmatrix)) & any(class(supplement) == "lefkoSD")) {
    checkconv <- supplement$convtype
    
    if (!is.element(3, checkconv)) {
      warning("Supplemental data does not include fecundity information, and a reproduction
        matrix has not been supplied. All fecundity transitions will be inferred from the
        stageframe.", call. = FALSE)
    }
  }
  
  stagenum_init <- dim(stageframe)[1]
  if (!all(is.null(repmatrix))) {
    if (any(class(repmatrix) == "matrix")) {
      if (dim(repmatrix)[1] != stagenum_init & dim(repmatrix)[1] != stagenum_init^2) {
        stop("The repmatrix provided must be a square matrix with dimensions
          equal to the number of stages in the stageframe, or the square thereof.",
          call. = FALSE)
      }
      
      if (dim(repmatrix)[2] != stagenum_init & dim(repmatrix)[2] != stagenum_init^2) {
        stop("The repmatrix provided must be a square matrix with dimensions
          equal to the number of stages in the stageframe, or the square thereof.",
          call. = FALSE)
      }
    }
  }
  
  if (any(!suppressWarnings(!is.na(as.numeric(as.character(stageframe$sizebin_center)))))) {
    stop("Function flefko3() requires size to be numeric rather than categorical.",
      call. = FALSE)
  }
  
  melchett <- .sf_reassess(stageframe, supplement, overwrite, repmatrix,
    agemat = FALSE, historical = TRUE, format = format_int)
  stageframe <- melchett$stageframe
  repmatrix <- melchett$repmatrix
  ovtable <- melchett$ovtable
  
  if (!all(is.na(overwrite)) | !all(is.na(supplement))) {
    if(any(duplicated(ovtable[,1:3]))) {
      stop("Multiple entries with different values for the same stage transition are not allowed
        in the supplemental or overwrite table. If modifying a historical table to perform an
        ahistorical analysis, then this may be due to different given rates of substitutions
        caused by dropping stage at occasion t-1. Please eliminate duplicate transitions.",
        call. = FALSE)
    }
  }
  
  # Next the data frame carrying all raw values and element indices for matrix element estimation
  allstages <- .theoldpizzle(stageframe, ovtable, repmatrix, finalage = 0,
    format = format_int, style = 0, cont = 0)
  
  maxsize <- max(c(allstages$size3, allstages$size2n, allstages$size2o), na.rm = TRUE)
  maxsizeb <- max(c(allstages$sizeb3, allstages$sizeb2n, allstages$sizeb2o), na.rm = TRUE)
  maxsizec <- max(c(allstages$sizec3, allstages$sizec2n, allstages$sizec2o), na.rm = TRUE)
  
  allstages <- allstages[(which(allstages$index321 != -1)),]
  
  # Now we work up the models
  if (class(modelsuite) == "lefkoMod") {
    if(is.na(surv_model)) {surv_model <- modelsuite$survival_model}
    if(is.na(obs_model)) {obs_model <- modelsuite$observation_model}
    if(is.na(size_model)) {size_model <- modelsuite$size_model}
    if(is.na(sizeb_model)) {sizeb_model <- modelsuite$sizeb_model}
    if(is.na(sizec_model)) {sizec_model <- modelsuite$sizec_model}
    if(is.na(repst_model)) {repst_model <- modelsuite$repstatus_model}
    if(is.na(fec_model)) {fec_model <- modelsuite$fecundity_model}
    
    if(is.na(jsurv_model)) {jsurv_model <- modelsuite$juv_survival_model}
    if(is.na(jobs_model)) {jobs_model <- modelsuite$juv_observation_model}
    if(is.na(jsize_model)) {jsize_model <- modelsuite$juv_size_model}
    if(is.na(jsizeb_model)) {jsizeb_model <- modelsuite$juv_sizeb_model}
    if(is.na(jsizec_model)) {jsizec_model <- modelsuite$juv_sizec_model}
    if(is.na(jrepst_model)) {jrepst_model <- modelsuite$juv_reproduction_model}
    if(is.na(jmatst_model)) {jmatst_model <- modelsuite$juv_maturity_model}
  }
  
  if (is.na(randomseed)) {
    set.seed(NULL)
  } else {
    set.seed(randomseed)
  }
  
  surv_proxy <- .modelextract(surv_model, paramnames, mainyears, mainpatches, 
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  obs_proxy <- .modelextract(obs_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  sigma <- 0
  rvarssummed <- 0
  sizedist <- 1
  
  size_proxy <- .modelextract(size_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  if (size_proxy$family == "poisson") {
    sizedist <- 0
    if (!all(is.na(size_proxy$variances))) {
      if (is.element("vcov", names(size_proxy$variances))) {
        rvarssummed <- sum(size_proxy$variances[,"vcov"])
      } else if (is.element("cond", names(size_proxy$variances))) {
        rvarssummed <- sum(unlist(lapply(c(1:length(size_proxy$variances[["cond"]])), function(X) {
          stddev.est <- attr(size_proxy$variances[["cond"]][[X]], "stddev")
          return(stddev.est * stddev.est)
        })))
      }
    } else {
      rvarssummed <- 0
    }
  } else if (size_proxy$family == "gaussian") {
    sizedist <- 2
    sigma <- size_proxy$sigma
  } else if (size_proxy$family == "gamma") {
    sizedist <- 3
    
  } else {
    sizedist <- 1
  }
  
  sigmab <- 0
  rvarssummedb <- 0
  sizebdist <- 1
  
  sizeb_proxy <- .modelextract(sizeb_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  if (sizeb_proxy$family == "poisson") {
    sizebdist <- 0
    if (!all(is.na(sizeb_proxy$variances))) {
      if (is.element("vcov", names(sizeb_proxy$variances))) {
        rvarssummedb <- sum(sizeb_proxy$variances[,"vcov"])
      } else if (is.element("cond", names(sizeb_proxy$variances))) {
        rvarssummedb <- sum(unlist(lapply(c(1:length(sizeb_proxy$variances[["cond"]])), function(X) {
          stddev.est <- attr(sizeb_proxy$variances[["cond"]][[X]], "stddev")
          return(stddev.est * stddev.est)
        })))
      }
    } else {
      rvarssummedb <- 0
    }
  } else if (sizeb_proxy$family == "gaussian") {
    sizebdist <- 2
    sigmab <- sizeb_proxy$sigma
  } else if (sizeb_proxy$family == "gamma") {
    sizebdist <- 3
    
  } else {
    sizebdist <- 1
  }
  
  sigmac <- 0
  rvarssummedc <- 0
  sizecdist <- 1
  
  sizec_proxy <- .modelextract(sizec_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  if (sizec_proxy$family == "poisson") {
    sizecdist <- 0
    if (!all(is.na(sizec_proxy$variances))) {
      if (is.element("vcov", names(sizec_proxy$variances))) {
        rvarssummedc <- sum(sizec_proxy$variances[,"vcov"])
      } else if (is.element("cond", names(sizec_proxy$variances))) {
        rvarssummedc <- sum(unlist(lapply(c(1:length(sizec_proxy$variances[["cond"]])), function(X) {
          stddev.est <- attr(sizec_proxy$variances[["cond"]][[X]], "stddev")
          return(stddev.est * stddev.est)
        })))
      }
    } else {
      rvarssummedc <- 0
    }
  } else if (sizec_proxy$family == "gaussian") {
    sizecdist <- 2
    sigmac <- sizec_proxy$sigma
  } else if (sizec_proxy$family == "gamma") {
    sizecdist <- 3
    
  } else {
    sizecdist <- 1
  }
  
  repst_proxy <- .modelextract(repst_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  fec_proxy <- .modelextract(fec_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  if (fec_proxy$family == "poisson") {
    fecdist <- 0
  } else if (fec_proxy$family == "gaussian") {
    fecdist <- 2
  } else if (fec_proxy$family == "gamma") {
    fecdist <- 3
  } else {
    fecdist <- 1
  }
  
  jsurv_proxy <- .modelextract(jsurv_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  jobs_proxy <- .modelextract(jobs_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  jsigma <- 0
  jrvarssummed <- 0
  
  jsize_proxy <- .modelextract(jsize_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  if (jsize_proxy$family == "poisson") {
    if (!all(is.na(jsize_proxy$variances))) {
      if (is.element("vcov", names(jsize_proxy$variances))) {
        jrvarssummed <- sum(jsize_proxy$variances[,"vcov"])
      } else if (is.element("cond", names(jsize_proxy$variances))) {
        jrvarssummed <- sum(unlist(lapply(c(1:length(jsize_proxy$variances[["cond"]])), function(X) {
          stddev.est <- attr(jsize_proxy$variances[["cond"]][[X]], "stddev")
          return(stddev.est * stddev.est)
        })))
      }
    } else {
      jrvarssummed <- 0
    }
  } else if (jsize_proxy$family == "gaussian") {
    jsigma <- jsize_proxy$sigma
  }
  
  jsigmab <- 0
  jrvarssummedb <- 0
  
  jsizeb_proxy <- .modelextract(jsizeb_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  if (jsizeb_proxy$family == "poisson") {
    if (!all(is.na(jsizeb_proxy$variances))) {
      if (is.element("vcov", names(jsizeb_proxy$variances))) {
        jrvarssummedb <- sum(jsizeb_proxy$variances[,"vcov"])
      } else if (is.element("cond", names(jsizeb_proxy$variances))) {
        jrvarssummedb <- sum(unlist(lapply(c(1:length(jsizeb_proxy$variances[["cond"]])), function(X) {
          stddev.est <- attr(jsizeb_proxy$variances[["cond"]][[X]], "stddev")
          return(stddev.est * stddev.est)
        })))
      }
    } else {
      jrvarssummedb <- 0
    }
  } else if (jsizeb_proxy$family == "gaussian") {
    jsigmab <- jsizeb_proxy$sigma
  }
  
  jsigmac <- 0
  jrvarssummedc <- 0
  
  jsizec_proxy <- .modelextract(jsizec_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  if (jsizec_proxy$family == "poisson") {
    if (!all(is.na(jsizec_proxy$variances))) {
      if (is.element("vcov", names(jsizec_proxy$variances))) {
        jrvarssummedc <- sum(jsizec_proxy$variances[,"vcov"])
      } else if (is.element("cond", names(jsizec_proxy$variances))) {
        jrvarssummedc <- sum(unlist(lapply(c(1:length(jsizec_proxy$variances[["cond"]])), function(X) {
          stddev.est <- attr(jsizec_proxy$variances[["cond"]][[X]], "stddev")
          return(stddev.est * stddev.est)
        })))
      }
    } else {
      jrvarssummedc <- 0
    }
  } else if (jsizec_proxy$family == "gaussian") {
    jsigmac <- jsizec_proxy$sigma
  }
  
  jrepst_proxy <- .modelextract(jrepst_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  jmatst_proxy <- .modelextract(jmatst_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  # This creates a list of pop, patch, and year in order of matrix
  if (!all(is.na(patch))) {
    listofyears <- apply(as.matrix(patch), 1, function(X) {
      output <- cbind.data.frame("1", X, as.matrix(year), stringsAsFactors = FALSE);
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    listofyears <- do.call(rbind.data.frame, listofyears)
    listofyears$poporder <- 1
    listofyears$patchorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainpatches == listofyears$patch[X])})
    listofyears$yearorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainyears == listofyears$year2[X])})
    
  } else {
    
    listofyears <- cbind.data.frame("1", "1", as.matrix(year), stringsAsFactors = FALSE)
    names(listofyears) <- c("pop", "patch", "year2")
    
    listofyears$poporder <- 1
    listofyears$patchorder <- 1
    listofyears$yearorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainyears == listofyears$year2[X])})
  }
  
  # A few extra tidbits required for the core matrix estimator to work
  yearlist <- split(listofyears, seq(nrow(listofyears)))
  
  # The next call runs the core matrix estimator function and creates all matrices
  madsexmadrigal <- lapply(yearlist, .jerzeibalowski, allstages, stageframe, 
    format_int, surv_proxy, obs_proxy, size_proxy, sizeb_proxy, sizec_proxy, 
    repst_proxy, fec_proxy, jsurv_proxy, jobs_proxy, jsize_proxy, jsizeb_proxy, 
    jsizec_proxy, jrepst_proxy, jmatst_proxy, f2.inda, f1.inda, f2.indb,
    f1.indb, f2.indc, f1.indc, r2.inda, r1.inda, r2.indb, r1.indb, r2.indc,
    r1.indc, c(surv_dev, obs_dev, size_dev, sizeb_dev, sizec_dev, repst_dev, fec_dev,
      jsurv_dev, jobs_dev, jsize_dev, jsizeb_dev, jsizec_dev, jrepst_dev, jmatst_dev),
    density, repmod, c(rvarssummed, sigma, rvarssummedb, sigmab, rvarssummedc,
      sigmac, jrvarssummed, jsigma, jrvarssummedb, jsigmab, jrvarssummedc, jsigmac),
    maxsize, maxsizeb, maxsizec, 0, sizedist, sizebdist, sizecdist, fecdist,
    negfec, exp_tol, theta_tol, ipm_method)
  
  a_list <- lapply(madsexmadrigal, function(X) {X$A})
  u_list <- lapply(madsexmadrigal, function(X) {X$U})
  f_list <- lapply(madsexmadrigal, function(X) {X$F})
  
  if (err_check) {
    out_list <- lapply(madsexmadrigal, function(X) {
        aaa <- X$out
        colnames(aaa) <- c("surv", "obs", "repst", "sizea", "sizeb", "sizec", "matst")
        
        return(aaa)
      }
    )
  }
  
  ahstages <- stageframe[1:(dim(stageframe)[1] - 1),]
  
  pairings1 <- expand.grid(stage_id_2 = stageframe$stage_id[1:(dim(stageframe)[1] - format_int)], 
    stage_id_1 = stageframe$stage_id[1:(dim(stageframe)[1] - 1)])
  pairings2 <- expand.grid(stage_2 = stageframe$stage[1:(dim(stageframe)[1] - format_int)], 
    stage_1 = stageframe$stage[1:(dim(stageframe)[1] - 1)])
  hstages <- cbind.data.frame(pairings1, pairings2)
  
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
  
  rownames(hstages) <- c(1:dim(hstages)[1])
  
  if (!err_check) {
    output <- list(A = a_list, U = u_list, F = f_list, hstages = hstages,
      agestages = NA, ahstages = ahstages, labels = listofyears[,c(1:3)],
      matrixqc = qcoutput1, modelqc = qcoutput2)
  } else {
    output <- list(A = a_list, U = u_list, F = f_list, hstages = hstages,
      agestages = NA, ahstages = ahstages, labels = listofyears[,c(1:3)],
      matrixqc = qcoutput1, modelqc = qcoutput2, prob_out = out_list)
  }
  
  class(output) <- "lefkoMat"
  
  return(output)
}

#' Reduce Matrix Dimensions By Eliminating Empty Stages
#' 
#' \code{.reducer3()} identifies empty stages in a set of historical matrices
#' and removes them from all matrices. It is used within \code{\link{flefko3}()}
#' and \code{\link{rlefko3}()}.
#' 
#' @param A List of population projection matrices, from a \code{lefkoMat}
#' object.
#' @param U List of surviva-transition matrices corresponding to \code{A}.
#' @param F List of fecundity matrices corresponding to \code{A}.
#' @param hstages Data frame giving the names and identities of historical stage
#' pairs used to create matrices.
#' 
#' @return Returns a list of reduced \code{A}, \code{U}, and \code{F} matrices,
#' plus the reduced \code{hstages} object.
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

#' Create Function-based Ahistorical Matrix Projection Model
#'
#' Function \code{flefko2()} returns ahistorical MPMs corresponding to the
#' patches and occasion times given, including the associated component
#' transition and fecundity matrices, a data frame detailing the characteristics
#' of the ahistorical stages used, and a data frame characterizing the patch and
#' occasion time combinations corresponding to these matrices.
#'
#' @param year A variable corresponding to observation occasion, or a set of
#' such values, given in values associated with the year term used in linear
#' model development. Defaults to \code{"all"}, in which case matrices will be
#' estimated for all occasion times.
#' @param patch A variable designating which patches or subpopulations will have
#' matrices estimated. Defaults to \code{"all"}, but can also be set to specific
#' patch names.
#' @param stageframe A stageframe object that includes information on the size,
#' observation status, propagule status, immaturity status, and maturity status
#' of each ahistorical stage. Should also incorporate bin widths if size is
#' continuous.
#' @param supplement An optional data frame of class \code{lefkoSD} that
#' provides supplemental data that should be incorporated into the MPM. Three
#' kinds of data may be integrated this way: transitions to be estimated via the
#' use of proxy transitions, transition overwrites from the literature or
#' supplemental studies, and transition multipliers for survival and fecundity.
#' This data frame should be produced using the \code{\link{supplemental}()}
#' function. Can be used in place of or in addition to an overwrite table (see 
#' \code{overwrite} below) and a reproduction matrix (see \code{repmatrix}
#' below).
#' @param repmatrix An optional reproduction matrix. This matrix is composed
#' mostly of 0s, with non-zero entries acting as element identifiers and
#' multipliers for fecundity (with 1 equaling full fecundity). If left blank,
#' and no \code{supplement} is provided, then \code{flefko2()} will assume that
#' all stages marked as reproductive produce offspring at 1x that of estimated
#' fecundity, and that offspring production will yield the first stage noted as
#' propagule or immature.  To prevent this behavior, input just \code{0}, which
#' will result in fecundity being estimated only for transitions noted in
#' \code{supplement} above. Must be the dimensions of an ahistorical matrix.
#' @param overwrite An optional data frame developed with the
#' \code{\link{overwrite}()} function describing transitions to be overwritten
#' either with given values or with other estimated transitions. Note that this
#' function supplements overwrite data provided in \code{supplement}.
#' @param data The historical vertical demographic data frame used to estimate
#' vital rates (class \code{hfvdata}). The original data frame is required in
#' order to initialize times and patches properly.
#' @param modelsuite An optional \code{lefkoMod} object holding the vital rate
#' models. If given, then \code{surv_model}, \code{obs_model}, 
#' \code{size_model}, \code{sizeb_model}, \code{sizec_model},
#' \code{repst_model}, \code{fec_model}, \code{jsurv_model}, \code{jobs_model},
#' \code{jsize_model}, \code{jsizeb_model}, \code{jsizec_model},
#' \code{jrepst_model}, \code{jmatst_model}, \code{paramnames}, \code{yearcol},
#' and \code{patchcol} are not required. No models should include size or
#' reproductive status in occasion \emph{t}-1.
#' @param surv_model A linear model predicting survival probability. This can
#' be a model of class \code{glm} or \code{glmer}, and requires a predicted
#' binomial variable under a logit link. If given, then will overwrite any
#' survival probability model given in \code{modelsuite}. This model must have
#' been developed in a modeling exercise testing only the impacts of occasion 
#' \emph{t}.
#' @param obs_model A linear model predicting sprouting or observation
#' probability. This can be a model of class \code{glm} or \code{glmer}, and
#' requires a predicted binomial variable under a logit link. If given, then
#' will overwrite any observation probability model given in \code{modelsuite}.
#' This model must have been developed in a modeling exercise testing only the
#' impacts of occasion \emph{t}.
#' @param size_model A linear model predicting primary size. This can be a model
#' of class \code{glm}, \code{glmer}, \code{glmmTMB}, \code{zeroinfl},
#' \code{vglm}, \code{lm}, or \code{lmer}. If given, then will overwrite any
#' primary size model given in \code{modelsuite}. This model must have been
#' developed in a modeling exercise testing only the impacts of occasion
#' \emph{t}.
#' @param sizeb_model A linear model predicting secondary size. This can be a
#' modelvof class \code{glm}, \code{glmer}, \code{glmmTMB}, \code{zeroinfl},
#' \code{vglm}, \code{lm}, or \code{lmer}. If given, then will overwrite any
#' secondary size model given in \code{modelsuite}. This model must have been
#' developed in a modeling exercise testing only the impacts of occasion
#' \emph{t}.
#' @param sizec_model A linear model predicting tertiary size. This can be a
#' modelvof class \code{glm}, \code{glmer}, \code{glmmTMB}, \code{zeroinfl},
#' \code{vglm}, \code{lm}, or \code{lmer}. If given, then will overwrite any
#' tertiary size model given in \code{modelsuite}. This model must have been
#' developed in a modeling exercise testing only the impacts of occasion
#' \emph{t}.
#' @param repst_model A linear model predicting reproduction probability. This
#' can be a model of class \code{glm} or \code{glmer}, and requires a predicted
#' binomial variable under a logit link. If given, then will overwrite any
#' reproduction probability model given in \code{modelsuite}. This model must
#' have been developed in a modeling exercise testing only the impacts of
#' occasion \emph{t}.
#' @param fec_model A linear model predicting fecundity. This can be a model of
#' class \code{glm}, \code{glmer}, \code{glmmTMB}, \code{zeroinfl},
#' \code{vglm}, \code{lm}, or \code{lmer}. If given, then will overwrite any 
#' fecundity model given in \code{modelsuite}. This model must have been
#' developed in a modeling exercise testing only the impacts of occasion
#' \emph{t}.
#' @param jsurv_model A linear model predicting juvenile survival probability.
#' This can be a model of class \code{glm} or \code{glmer}, and requires a
#' predicted binomial variable under a logit link. If given, then will overwrite
#' any juvenile survival probability model given in \code{modelsuite}. This
#' model must have been developed in a modeling exercise testing only the
#' impacts of occasion \emph{t}.
#' @param jobs_model A linear model predicting juvenile sprouting or observation
#' probability. This can be a model of class \code{glm} or \code{glmer}, and
#' requires a predicted binomial variable under a logit link. If given, then
#' will overwrite any juvenile observation probability model given in 
#' \code{modelsuite}. This model must have been developed in a modeling exercise
#' testing only the impacts of occasion \emph{t}.
#' @param jsize_model A linear model predicting juvenile primary size. This
#' can be a model of class \code{glm}, \code{glmer}, \code{glmmTMB},
#' \code{zeroinfl}, \code{vglm}, \code{lm}, or \code{lmer}. If given, then will
#' overwrite any juvenile primary size model given in \code{modelsuite}. This
#' model must have been developed in a modeling exercise testing only the
#' impacts of occasion \emph{t}.
#' @param jsizeb_model A linear model predicting juvenile secondary size. This
#' can be a model of class \code{glm}, \code{glmer}, \code{glmmTMB},
#' \code{zeroinfl}, \code{vglm}, \code{lm}, or \code{lmer}. If given, then will
#' overwrite any juvenile secondary size model given in \code{modelsuite}. This
#' model must have been developed in a modeling exercise testing only the
#' impacts of occasion \emph{t}.
#' @param jsizec_model A linear model predicting juvenile tertiary size. This
#' can be a model of class \code{glm}, \code{glmer}, \code{glmmTMB},
#' \code{zeroinfl}, \code{vglm}, \code{lm}, or \code{lmer}. If given, then will
#' overwrite any juvenile tertiary size model given in \code{modelsuite}. This
#' model must have been developed in a modeling exercise testing only the
#' impacts of occasion \emph{t}.
#' @param jrepst_model A linear model predicting reproduction probability of a 
#' mature individual that was immature in time \emph{t}. This can be a model
#' of class \code{glm} or \code{glmer}, and requires a predicted binomial
#' variable under a logit link. If given, then will overwrite any reproduction
#' probability model given in \code{modelsuite}. This model must have been
#' developed in a modeling exercise testing only the impacts of occasion
#' \emph{t}.
#' @param jmatst_model A linear model predicting maturity probability of an 
#' individual that was immature in time \emph{t}. This can be a model of class
#' \code{glm} or \code{glmer}, and requires a predicted binomial variable under
#' a logit link. If given, then will overwrite any maturity probability model
#' given in \code{modelsuite}. This model must have been developed in a modeling
#' exercise testing only the impacts of occasion \emph{t}.
#' @param paramnames A dataframe with three columns, the first describing all
#' terms used in linear modeling, the second (must be called \code{mainparams}),
#' showing the general model terms that will be used in matrix creation (users
#' should use \code{\link{modelsearch}()} at least once to see the proper
#' names to be used in this column), and the third showing the equivalent terms
#' used in modeling (must be named \code{modelparams}). Only required if
#' \code{modelsuite} is not supplied.
#' @param inda Can be a single value to use for individual covariate \code{a}
#' in all matrices, a pair of values to use for times \emph{t} and \emph{t}-1 in
#' historical matrices, or a vector of such values corresponding to each
#' occasion in option \code{year}. Defaults to \code{NULL}.
#' @param indb Can be a single value to use for individual covariate \code{b}
#' in all matrices, a pair of values to use for times \emph{t} and \emph{t}-1 in
#' historical matrices, or a vector of such values corresponding to each
#' occasion in option \code{year}. Defaults to \code{NULL}.
#' @param indc Can be a single value to use for individual covariate \code{c}
#' in all matrices, a pair of values to use for times \emph{t} and \emph{t}-1 in
#' historical matrices, or a vector of such values corresponding to each
#' occasion in option \code{year}. Defaults to \code{NULL}.
#' @param surv_dev A numeric value to be added to the y-intercept in the linear
#' model for survival probability.
#' @param obs_dev A numeric value to be added to the y-intercept in the linear
#' model for observation probability.
#' @param size_dev A numeric value to be added to the y-intercept in the linear
#' model for primary size.
#' @param sizeb_dev A numeric value to be added to the y-intercept in the linear
#' model for secondary size.
#' @param sizec_dev A numeric value to be added to the y-intercept in the linear
#' model for tertiary size.
#' @param repst_dev A numeric value to be added to the y-intercept in the linear
#' model for probability of reproduction.
#' @param fec_dev A numeric value to be added to the y-intercept in the linear
#' model for fecundity.
#' @param jsurv_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile survival probability.
#' @param jobs_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile observation probability.
#' @param jsize_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile primary size.
#' @param jsizeb_dev A numeric value to be added to the y-intercept in the
#' linear model for juvenile secondary size.
#' @param jsizec_dev A numeric value to be added to the y-intercept in the
#' linear model for juvenile tertiary size.
#' @param jrepst_dev A numeric value to be added to the y-intercept in the
#' linear model for juvenile reproduction probability.
#' @param jmatst_dev A numeric value to be added to the y-intercept in the
#' linear model for juvenile maturity probability.
#' @param density A numeric value indicating density value to use to propagate
#' matrices. Only needed if density is an explanatory term used in linear
#' models. Defaults to \code{NA}.
#' @param repmod A scalar multiplier of fecundity. Defaults to \code{1}.
#' @param yearcol The variable name or column number corresponding to year in
#' occasion \emph{t} in the dataset. Not needed if a \code{modelsuite} is
#' supplied.
#' @param patchcol The variable name or column number corresponding to patch in
#' the dataset. Not needed if a \code{modelsuite} is supplied.
#' @param year.as.random A logical term indicating whether coefficients for
#' missing occasions within vital rate models should be estimated as random
#' intercepts. Defaults to \code{FALSE}, in which case missing monitoring
#' occasion coefficients are set to \code{0}.
#' @param patch.as.random A logical term indicating whether coefficients for
#' missing patches within vital rate models should be estimated as random
#' intercepts. Defaults to \code{FALSE}, in which case missing patch
#' coefficients are set to \code{0}.
#' @param random.inda A logical value denoting whether to treat individual
#' covariate \code{a} as a random, categorical variable. Otherwise is treated as
#' a fixed, numeric variable. Defaults to \code{FALSE}.
#' @param random.indb A logical value denoting whether to treat individual
#' covariate \code{b} as a random, categorical variable. Otherwise is treated as
#' a fixed, numeric variable. Defaults to \code{FALSE}.
#' @param random.indc A logical value denoting whether to treat individual
#' covariate \code{c} as a random, categorical variable. Otherwise is treated as
#' a fixed, numeric variable. Defaults to \code{FALSE}.
#' @param randomseed A numeric value used as a seed to generate random estimates
#' for missing occasion and patch coefficients, if either \code{year.as.random}
#' or \code{patch.as.random} is set to \code{TRUE}. Defaults to
#' \code{\link{set.seed}()} default.
#' @param negfec A logical value denoting whether fecundity values estimated to
#' be negative should be reset to \code{0}. Defaults to \code{FALSE}.
#' @param ipm_method A string indicating what method to use to estimate size
#' transition probabilities if size is treated as continuous. Options include
#' \code{"midpoint"}, which utilizes the midpoint method, and \code{"CDF"},
#' which uses the cumulative distribution function. Defaults to \code{"CDF"}.
#' @param reduce A logical value denoting whether to remove ahistorical stages
#' associated solely with 0 transitions. These are only removed in cases where
#' the associated row and column sums in ALL matrices estimated equal 0.
#' Defaults to \code{FALSE}.
#' @param err_check A logical value indicating whether to append matrices of
#' vital rate probabilities associated with each matrix to the \code{lefkoMat}
#' object generated. These matrices are developed internally and can be used for
#' error checking. Defaults to \code{FALSE}.
#' @param exp_tol A numeric value used to indicate a maximum value to set
#' exponents to in the core kernel to prevent numerical overflow. Defaults to
#' \code{700}.
#' @param theta_tol A numeric value used to indicate a maximum value to theta as
#' used in the negative binomial probability density kernel. Defaults to
#' \code{100000000}, but can be reset to other values during error checking.
#'
#' @return If all inputs are properly formatted, then this function will return
#' an object of class \code{lefkoMat}, which is a list that holds the matrix
#' projection model and all of its metadata. Its structure is a list with the
#' following elements:
#'
#' \item{A}{A list of full projection matrices in order of sorted patches and
#' occasion times. All matrices output in R's \code{matrix} class.}
#' \item{U}{A list of survival transition matrices sorted as in \code{A}. All 
#' matrices output in R's \code{matrix} class.}
#' \item{F}{A list of fecundity matrices sorted as in \code{A}. All matrices 
#' output in R's \code{matrix} class.}
#' \item{hstages}{A data frame matrix showing the pairing of ahistorical stages
#' used to create historical stage pairs. Set to \code{NA} for ahistorical
#' matrices.}
#' \item{agestages}{A data frame showing age-stage pairs. In this function, it
#' is set to \code{NA}. Only used in output to function \code{aflefko2}().}
#' \item{ahstages}{A data frame detailing the characteristics of associated
#' ahistorical stages, in the form of a modified stageframe that includes
#' status as an entry stage through reproduction.}
#' \item{labels}{A data frame giving the population, patch, and year of each
#' matrix in order. In \code{flefko2()}, only one population may be analyzed at
#' once, and so \code{pop = NA}.}
#' \item{matrixqc}{A short vector describing the number of non-zero elements in
#' \code{U} and \code{F} matrices, and the number of matrices.}
#' \item{modelqc}{This is the \code{qc} portion of the modelsuite input.}
#' \item{prob_out}{An optional element only added if \code{err_check = TRUE}.
#' This is a list of vital rate probability matrices, with 6 columns in the
#' order of survival, observation probability, reproduction probability, primary
#' size transition probability, secondary size transition probability, and
#' tertiary size transition probability.}
#' 
#' @section Notes:
#' Unlike \code{\link{rlefko2}()}, \code{\link{rlefko3}()}, and
#' \code{\link{rleslie}()}, this function does not currently distinguish
#' populations.
#' 
#' This function will yield incorrect estimates if the models utilized
#' incorporate state in occasion \emph{t}-1. Only use models developed testing
#' for ahistorical effects.
#' 
#' The default behavior of this function is to estimate fecundity with regards
#' to transitions specified via associated fecundity multipliers in either
#' \code{supplement} or \code{repmatrix}. If both of these fields are left
#' empty, then fecundity will be estimated at full for all transitions leading
#' from reproductive stages to immature and propagule stages. However, if a
#' \code{supplement} is provided and a \code{repmatrix} is not, or if
#' \code{repmatrix} is set to \code{0}, then only fecundity transitions noted in
#' the \code{supplement} will be set to non-zero values. To use the default
#' behavior of setting all reproductive stages to reproduce at full fecundity
#' into immature and propagule stages, but also incorporate given or proxy
#' survival transitions, input those given and proxy transitions through the
#' \code{overwrite} option.
#' 
#' The reproduction matrix (field \code{repmatrix}) may only be supplied as
#' ahistorical. If provided as historical, then \code{flefko2()} will fail and
#' produce an error.
#' 
#' Users may at times wish to estimate MPMs using a dataset incorporating
#' multiple patches or subpopulations, but without discriminating between those
#' patches or subpopulations. Should the aim of analysis be a general MPM that
#' does not distinguish these patches or subpopulations, the \code{patchcol}
#' variable should be set to \code{NA}, which is the default.
#'
#' Input options including multiple variable names must be entered in the order
#' of variables in occasion \emph{t}+1 and \emph{t}. Rearranging the order will
#' lead to erroneous calculations, and may lead to fatal errors.
#'
#' Care should be taken to match the random status of year and patch to the
#' states of those variables within the modelsuite. If they do not match, then
#' they will be treated as zeroes in vital rate estimation.
#' 
#' The \code{ipm_method} function gives the option of using two different means
#' of estimating the probability of size transition. The midpoint method
#' (\code{"midpoint"}) refers to the method in which the probability is
#' estimated by first estimating the probability associated with transition from
#' the exact size at the midpoint of the size class using the corresponding
#' probability density function, and then multiplying that value by the bin
#' width of the size class. Doak et al. 2021 (Ecological Monographs) noted that
#' this method can produce biased results, with total size transitions
#' associated with a specific size not totaling to 1.0 and even specific size
#' transition probabilities capable of being estimated at values greater than
#' 1.0. The alternative and default method, \code{"CDF"}, uses the corresponding
#' cumulative density function to estimate the probability of size transition as
#' the cumulative probability of size transition at the greater limit of the
#' size class minus the cumulative probability of size transition at the lower
#' limit of the size class. The latter method avoids this bias. Note, however,
#' that both methods are exact and unbiased for negative binomial and Poisson
#' distributions.
#' 
#' Under the Gaussian and gamma size distributions, the number of estimated
#' parameters may differ between the two \code{ipm_method} settings. Because
#' the midpoint method has a tendency to incorporate upward bias in the
#' estimation of size transition probabilities, it is more likely to yield non-
#' zero values when the true probability is extremely close to 0. This will
#' result in the \code{summary.lefkoMat} function yielding higher numbers of
#' estimated parameters than the \code{ipm_method = "CDF"} yields in some cases.
#' 
#' Using the \code{err_check} option will produce a matrix of 6 columns, each
#' characterizing a different vital rate. The product of each row yields an
#' element in the associated \code{$U} matrix. The number and order of elements
#' in each column of this matrix matches the associated matrix in column vector
#' format. Use of this option is generally for the purposes of debugging code.
#'
#' @seealso \code{\link{flefko3}()}
#' @seealso \code{\link{aflefko2}()}
#' @seealso \code{\link{fleslie}()}
#' @seealso \code{\link{rlefko3}()}
#' @seealso \code{\link{rlefko2}()}
#' @seealso \code{\link{rleslie}()}
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
#' # Here we use supplemental to provide overwrite and reproductive info
#' lathsupp2 <- supplemental(stage3 = c("Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "rep", "rep"),
#'   givenrate = c(0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 3, 3), stageframe = lathframeln, historical = FALSE)
#' 
#' lathmat2ln <- flefko2(year = "all", patch = "all", stageframe = lathframeln, 
#'   modelsuite = lathmodelsln2, data = lathvertln, supplement = lathsupp2,
#'   patchcol = "patchid", yearcol = "year2", year.as.random = FALSE,
#'   patch.as.random = FALSE, reduce = FALSE)
#' 
#' summary(lathmat2ln)
#' 
#' #Cypripedium example using three size metrics for classification
#' rm(list=ls(all=TRUE))
#' 
#' data(cypdata)
#' sizevector.f <- c(0, 0, 0, 0, 0, 0, seq(1, 12, by = 1), seq(0, 9, by = 1),
#'   seq(0, 8, by = 1), seq(0, 7, by = 1), seq(0, 6, by = 1), seq(0, 5, by = 1),
#'   seq(0, 4, by = 1), seq(0, 3, by = 1), 0, 1, 2, 0, 1, 0, 
#'   0, 0, 1, 0)
#' sizebvector.f <- c(0, 0, 0, 0, 0, 0, rep(0, 12), rep(1, 10), rep(2, 9),
#'   rep(3, 8), rep(4, 7), rep(5, 6), rep(6, 5), rep(7, 4), rep(8, 3), 9, 9, 10, 
#'   0, 1, 1, 2)
#' sizecvector.f <- c(0, 0, 0, 0, 0, 0, rep(0, 12), rep(0, 10), rep(0, 9),
#'   rep(0, 8), rep(0, 7), rep(0, 6), rep(0, 5), rep(0, 4), 0, 0, 0, 0, 0, 0, 
#'   1, 1, 1, 1)
#' stagevector.f <- c("DS", "P1", "P2", "P3", "Sdl", "Dorm", "V1 I0 D0",
#'   "V2 I0 D0", "V3 I0 D0", "V4 I0 D0", "V5 I0 D0", "V6 I0 D0", "V7 I0 D0",
#'   "V8 I0 D0", "V9 I0 D0", "V10 I0 D0", "V11 I0 D0", "V12 I0 D0", "V0 I1 D0",
#'   "V1 I1 D0", "V2 I1 D0", "V3 I1 D0", "V4 I1 D0", "V5 I1 D0", "V6 I1 D0",
#'   "V7 I1 D0", "V8 I1 D0", "V9 I1 D0", "V0 I2 D0", "V1 I2 D0", "V2 I2 D0",
#'   "V3 I2 D0", "V4 I2 D0", "V5 I2 D0", "V6 I2 D0", "V7 I2 D0", "V8 I2 D0",
#'   "V0 I3 D0", "V1 I3 D0", "V2 I3 D0", "V3 I3 D0", "V4 I3 D0", "V5 I3 D0",
#'   "V6 I3 D0", "V7 I3 D0", "V0 I4 D0", "V1 I4 D0", "V2 I4 D0", "V3 I4 D0",
#'   "V4 I4 D0", "V5 I4 D0", "V6 I4 D0", "V0 I5 D0", "V1 I5 D0", "V2 I5 D0",
#'   "V3 I5 D0", "V4 I5 D0", "V5 I5 D0", "V0 I6 D0", "V1 I6 D0", "V2 I6 D0",
#'   "V3 I6 D0", "V4 I6 D0", "V0 I7 D0", "V1 I7 D0", "V2 I7 D0", "V3 I7 D0",
#'   "V0 I8 D0", "V1 I8 D0", "V2 I8 D0", "V0 I9 D0", "V1 I9 D0", "V0 I10 D0",
#'   "V0 I0 D1", "V0 I1 D1", "V1 I1 D1", "V0 I2 D1")
#' repvector.f <- c(0, 0, 0, 0, 0, rep(0, 13), rep(1, 59))
#' obsvector.f <- c(0, 0, 0, 0, 0, 0, rep(1, 71))
#' matvector.f <- c(0, 0, 0, 0, 0, rep(1, 72))
#' immvector.f <- c(0, 1, 1, 1, 1, rep(0, 72))
#' propvector.f <- c(1, rep(0, 76))
#' indataset.f <- c(0, 0, 0, 0, 0, rep(1, 72))
#' binvec.f <- c(0, 0, 0, 0, 0, rep(0.5, 72))
#' binbvec.f <- c(0, 0, 0, 0, 0, rep(0.5, 72))
#' bincvec.f <- c(0, 0, 0, 0, 0, rep(0.5, 72))
#' 
#' vertframe.f <- sf_create(sizes = sizevector.f, sizesb = sizebvector.f,
#'   sizesc = sizecvector.f, stagenames = stagevector.f, repstatus = repvector.f,
#'   obsstatus = obsvector.f, propstatus = propvector.f, immstatus = immvector.f,
#'   matstatus = matvector.f, indataset = indataset.f, binhalfwidth = binvec.f,
#'   binhalfwidthb = binbvec.f, binhalfwidthc = bincvec.f)
#' 
#' vert.data.f <- verticalize3(cypdata, noyears = 6, firstyear = 2004,
#'   individcol = "plantid", blocksize = 4, sizeacol = "Veg.04",
#'   sizebcol = "Inf.04", sizeccol = "Inf2.04", repstracol = "Inf.04",
#'   repstrbcol = "Inf2.04", fecacol = "Pod.04", censorcol = "censor",
#'   censorkeep = 1, censorRepeat = FALSE, stageassign = vertframe.f,
#'   stagesize = "sizeabc", NAas0 = TRUE, censor = FALSE)
#' 
#' vertmodels2f <- modelsearch(vert.data.f, historical = FALSE, suite = "main", 
#'   sizeb = c("sizeb3", "sizeb2", "sizeb1"), sizec = c("sizec3", "sizec2", "sizec1"),
#'   approach = "glm", vitalrates = c("surv", "obs", "size", "repst", "fec"),
#'   sizedist = "negbin", sizebdist = "poisson", sizecdist = "poisson",
#'   fecdist = "poisson", patch.as.random = TRUE, year.as.random = TRUE)
#' 
#' vertsupp2f <- supplemental(stage3 = c("DS", "P1", "P2", "P3", "Sdl", "Sdl",
#'     "Dorm", "V1 I0 D0", "V2 I0 D0", "V3 I0 D0", "DS", "P1"),
#'   stage2 = c("DS", "DS", "P1", "P2", "P3", "Sdl", "Sdl", "Sdl", "Sdl", "Sdl",
#'     "rep", "rep"), 
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "Dorm", "V1 I0 D0", "V2 I0 D0",
#'     "V3 I0 D0", NA, NA), 
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "V1 I0 D0", "V1 I0 D0", "V1 I0 D0",
#'     "V1 I0 D0", NA, NA), 
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, 0.40, NA, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5 * 5000, 0.5 * 5000),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3), stageframe = vertframe.f,
#'   historical = FALSE)
#' 
#' vert.mats.f2 <- flefko2(stageframe = vertframe.f, supplement = vertsupp2f, 
#'   data = vert.data.f, modelsuite = vertmodels2f)
#' summary(vert.mats.f2)
#' }
#' 
#' @export
flefko2 <- function(year = "all", patch = "all", stageframe, supplement = NULL,
  repmatrix = NULL, overwrite = NULL, data = NA, modelsuite = NA,
  surv_model = NA, obs_model = NA, size_model = NA, sizeb_model = NA,
  sizec_model = NA, repst_model = NA, fec_model = NA, jsurv_model = NA,
  jobs_model = NA, jsize_model = NA, jsizeb_model = NA, jsizec_model = NA,
  jrepst_model = NA, jmatst_model = NA, paramnames = NA, inda = NULL, indb = NULL, indc = NULL,
  surv_dev = 0, obs_dev = 0, size_dev = 0, sizeb_dev = 0, sizec_dev = 0,
  repst_dev = 0, fec_dev = 0, jsurv_dev = 0, jobs_dev = 0, jsize_dev = 0,
  jsizeb_dev = 0, jsizec_dev = 0, jrepst_dev = 0, jmatst_dev = 0, density = NA, repmod = 1,
  yearcol = NA, patchcol = NA, year.as.random = FALSE, patch.as.random = FALSE,
  random.inda = FALSE, random.indb = FALSE, random.indc = FALSE,
  randomseed = NA, negfec = FALSE, ipm_method = "CDF", reduce = FALSE,
  err_check = FALSE, exp_tol = 700, theta_tol = 100000000) {
  
  indanames <- indbnames <- indcnames <- NULL
  
  if (all(is.na(modelsuite)) & all(is.na(paramnames))) {
    warning("Function may not work properly without a dataframe of model parameters or
      equivalents supplied either through the modelsuite option or through the paramnames
      input parameter.", call. = FALSE)
  } else if (!all(is.na(modelsuite))) {
    paramnames <- modelsuite$paramnames
    yearcol <- paramnames$modelparams[which(paramnames$mainparams == "year2")]
    patchcol <- paramnames$modelparams[which(paramnames$mainparams == "patch")]
  }
  
  ipm_method <- tolower(ipm_method)
  if (length(grep("mi", ipm_method)) > 0) {
    ipm_method <- "midpoint"
  } else if (length(grep("cd", ipm_method)) > 0) {
    ipm_method <- "cdf"
  } else {
    stop("Option ipm_method not recognized.", call. = FALSE)
  }
  
  if (all(is.na(data))) {
    stop("Need original vertical dataset to set proper limits on year and patch.",
      call. = FALSE)
  }
  if (!any(class(data) == "data.frame")) {
    stop("Need original vertical dataset used in modeling to proceed.",
      call. = FALSE)
  }
  if (!any(class(data) == "hfvdata")) {
    warning("Dataset used as input is not of class hfvdata. Will assume that the
      dataset has been formatted equivalently.", call. = FALSE)
  }
  
  stageframe_vars <- c("stage", "size", "size_b", "size_c", "min_age", "max_age",
    "repstatus", "obsstatus", "propstatus", "immstatus", "matstatus", "indataset",
    "binhalfwidth_raw", "sizebin_min", "sizebin_max", "sizebin_center",
    "sizebin_width", "binhalfwidthb_raw", "sizebinb_min", "sizebinb_max",
    "sizebinb_center", "sizebinb_width", "binhalfwidthc_raw", "sizebinc_min",
    "sizebinc_max", "sizebinc_center", "sizebinc_width", "group", "comments")
  if (any(!is.element(names(stageframe), stageframe_vars))) {
    stop("Please use properly formatted stageframe as input.", call. = FALSE)
  }
  
  if (is.character(yearcol)) {
    choicevar <- which(names(data) == yearcol);
    mainyears <- sort(unique(data[,choicevar]))
  } else if (is.numeric(yearcol)) {
    mainyears <- sort(unique(data[, yearcol]));
  } else {
    stop("Need appropriate year column designation.", call. = FALSE)
  }
  
  if (any(is.character(year))) {
    if (is.element("all", tolower(year))) {
      year <- mainyears
    } else {
      stop("Year designation not recognized.", call. = FALSE)
    }
  }
  
  if (length(year) == 0 | all(is.na(year) == TRUE) | any(is.na(year))) {
    stop("This function cannot proceed without a specific occasion, or a suite of
      occasions, designated via the year option. NA entries are not allowed.",
      call. = FALSE)
  }
  
  if (all(is.na(patch)) & !is.na(patchcol)) {
    warning("Matrix creation may not proceed properly without input in the patch
      option when using a modelsuite in which patch is designated.",
      call. = FALSE)
  }
  
  if (is.character(patchcol) & patchcol != "none") {
    choicevar <- which(names(data) == patchcol);
    mainpatches <- sort(unique(as.character(data[,choicevar])))
  } else if (is.numeric(patchcol)) {
    mainpatches <- sort(unique(as.character(data[, patchcol])));
  } else {
    mainpatches <- NA
  }
  
  if (any(is.character(patch))) {
    if (is.element("all", tolower(patch))) {
      patch <- mainpatches
    } else if (!is.element(patch, mainpatches)) {
      stop("Patch designation not recognized.", call. = FALSE)
    }
  }
  
  if (!is.null(inda)) {
    if (!is.numeric(inda) & !random.inda) {
      stop("Individual covariate vector a must be numeric if not set to random.",
        call. = FALSE)
    }
    
    if (!is.element(length(inda), c(1, 2, length(year)))) {
      stop("Individual covariate vector a must be empty, or include 1, 2, or as
        many elements as occasions to be modeled.", call. = FALSE)
    }
    
    if (random.inda) {
      indacol <- paramnames$modelparams[which(paramnames$mainparams == "indcova2")]
      if (indacol == "none") {
        stop("Individual covariate a not recognized in the modelsuite", call. = FALSE)
      }
      
      indacol <- which(names(data) == indacol)
      
      if (length(indacol) > 0) {
        indanames <- sort(unique(data[, indacol]))
      } else {
        stop("Individual covariate a not recognized in the modelsuite", call. = FALSE)
      }
      
      if (any(!is.element(inda, indanames))) {
        stop("Entered value for individual covariate a does not exist in the data.",
          call. = FALSE)
      }
      
      if (length(inda) == 1) {
        r1.inda <- rep(as.character(inda), length(year))
        r2.inda <- rep(as.character(inda), length(year))
      } else if (length(inda) == 2 & length(year) != 2) {
        r1.inda <- rep(as.character(inda[1]), length(year))
        r2.inda <- rep(as.character(inda[2]), length(year))
      } else if (length(inda) == length(year)) {
        r2.inda <- as.character(inda)
        r1.inda <- c("none", r2.inda[1:(length(inda) - 1)])
      }
      
      f1.inda <- rep(0, length(year))
      f2.inda <- rep(0, length(year))
      
    } else {
      indanames <- c(0)
      
      if (length(inda) == 1) {
        f1.inda <- rep(inda, length(year))
        f2.inda <- rep(inda, length(year))
      } else if (length(inda) == 2 & length(year) != 2) {
        f1.inda <- rep(inda[1], length(year))
        f2.inda <- rep(inda[2], length(year))
      } else if (length(inda) == length(year)) {
        f2.inda <- inda
        f1.inda <- c(0, f2.inda[1:(length(inda) - 1)])
      }
      r2.inda <- rep("none", length(year))
      r1.inda <- rep("none", length(year))
    }
  } else {
    indanames <- c(0)
    
    f1.inda <- rep(0, length(year))
    f2.inda <- rep(0, length(year))
    r2.inda <- rep("none", length(year))
    r1.inda <- rep("none", length(year))
  }
  
  if (!is.null(indb)) {
    if (!is.numeric(indb) & !random.indb) {
      stop("Individual covariate vector b must be numeric if not set to random.",
        call. = FALSE)
    }
    
    if (!is.element(length(indb), c(1, 2, length(year)))) {
      stop("Individual covariate vector b must be empty, or include 1, 2, or as
        many elements as occasions to be modeled.", call. = FALSE)
    }
    
    if (random.indb) {
      indbcol <- paramnames$modelparams[which(paramnames$mainparams == "indcovb2")]
      if (indbcol == "none") {
        stop("Individual covariate b not recognized in the modelsuite", call. = FALSE)
      }
      
      indbcol <- which(names(data) == indbcol)
      
      if (length(indbcol) > 0) {
        indbnames <- sort(unique(data[, indbcol]))
      } else {
        stop("Individual covariate b not recognized in the modelsuite", call. = FALSE)
      }
      
      if (any(!is.element(indb, indbnames))) {
        stop("Entered value for individual covariate b does not exist in the data.",
          call. = FALSE)
      }
      
      if (length(indb) == 1) {
        r1.indb <- rep(as.character(indb), length(year))
        r2.indb <- rep(as.character(indb), length(year))
      } else if (length(indb) == 2 & length(year) != 2) {
        r1.indb <- rep(as.character(indb[1]), length(year))
        r2.indb <- rep(as.character(indb[2]), length(year))
      } else if (length(indb) == length(year)) {
        r2.indb <- as.character(indb)
        r1.indb <- c("none", r2.indb[1:(length(indb) - 1)])
      }
      
      f1.indb <- rep(0, length(year))
      f2.indb <- rep(0, length(year))
      
    } else {
      indbnames <- c(0)
      
      if (length(indb) == 1) {
        f1.indb <- rep(indb, length(year))
        f2.indb <- rep(indb, length(year))
      } else if (length(indb) == 2 & length(year) != 2) {
        f1.indb <- rep(indb[1], length(year))
        f2.indb <- rep(indb[2], length(year))
      } else if (length(indb) == length(year)) {
        f2.indb <- indb
        f1.indb <- c(0, f2.indb[1:(length(indb) - 1)])
      }
      r2.indb <- rep("none", length(year))
      r1.indb <- rep("none", length(year))
    }
  } else {
    indbnames <- c(0)
    
    f1.indb <- rep(0, length(year))
    f2.indb <- rep(0, length(year))
    r2.indb <- rep("none", length(year))
    r1.indb <- rep("none", length(year))
  }
  
  if (!is.null(indc)) {
    if (!is.numeric(indc) & !random.indc) {
      stop("Individual covariate vector c must be numeric if not set to random.",
        call. = FALSE)
    }
    
    if (!is.element(length(indc), c(1, 2, length(year)))) {
      stop("Individual covariate vector c must be empty, or include 1, 2, or as
        many elements as occasions to be modeled.", call. = FALSE)
    }
    
    if (random.indc) {
      indccol <- paramnames$modelparams[which(paramnames$mainparams == "indcovc2")]
      if (indccol == "none") {
        stop("Individual covariate c not recognized in the modelsuite", call. = FALSE)
      }
      
      indccol <- which(names(data) == indccol)
      
      if (length(indccol) > 0) {
        indcnames <- sort(unique(data[, indccol]))
      } else {
        stop("Individual covariate c not recognized in the modelsuite", call. = FALSE)
      }
      
      if (any(!is.element(indc, indcnames))) {
        stop("Entered value for individual covariate c does not exist in the data.",
          call. = FALSE)
      }
      
      if (length(indc) == 1) {
        r1.indc <- rep(as.character(indc), length(year))
        r2.indc <- rep(as.character(indc), length(year))
      } else if (length(indc) == 2 & length(year) != 2) {
        r1.indc <- rep(as.character(indc[1]), length(year))
        r2.indc <- rep(as.character(indc[2]), length(year))
      } else if (length(indc) == length(year)) {
        r2.indc <- as.character(indc)
        r1.indc <- c("none", r2.indc[1:(length(indc) - 1)])
      }
      
      f1.indc <- rep(0, length(year))
      f2.indc <- rep(0, length(year))
      
    } else {
      indcnames <- c(0)
      
      if (length(indc) == 1) {
        f1.indc <- rep(indc, length(year))
        f2.indc <- rep(indc, length(year))
      } else if (length(indc) == 2 & length(year) != 2) {
        f1.indc <- rep(indc[1], length(year))
        f2.indc <- rep(indc[2], length(year))
      } else if (length(indc) == length(year)) {
        f2.indc <- indc
        f1.indc <- c(0, f2.indc[1:(length(indc) - 1)])
      }
      r2.indc <- rep("none", length(year))
      r1.indc <- rep("none", length(year))
    }
  } else {
    indcnames <- c(0)
    
    f1.indc <- rep(0, length(year))
    f2.indc <- rep(0, length(year))
    r2.indc <- rep("none", length(year))
    r1.indc <- rep("none", length(year))
  }
  
  maingroups <- sort(unique(stageframe$group))
  
  if (!all(is.na(density))) {
    if (!all(is.numeric(density))) {
      stop("Density value must be numeric.", call. = FALSE)
    }
    
    if (any(is.na(density))) {
      density[which(is.na(density))] <- 0
    }
  } else {
    density <- 0
  }
  
  if (all(is.na(repmatrix)) & all(is.na(supplement))) {
    warning("Neither supplemental data nor a reproduction matrix have been supplied.
      All fecundity transitions will be inferred from the stageframe.",
      call. = FALSE)
  } else if (all(is.na(repmatrix)) & any(class(supplement) == "lefkoSD")) {
    checkconv <- supplement$convtype
    
    if (!is.element(3, checkconv)) {
      warning("Supplemental data does not include fecundity information, and a reproduction
        matrix has not been supplied. All fecundity transitions will be inferred from the
        stageframe.", call. = FALSE)
    }
  }
  
  stagenum_init <- dim(stageframe)[1]
  if (!all(is.null(repmatrix))) {
    if (any(class(repmatrix) == "matrix")) {
      if (dim(repmatrix)[1] != stagenum_init | dim(repmatrix)[2] != stagenum_init) {
        stop("The repmatrix provided must be a square matrix with dimensions
          equal to the number of stages in the stageframe.", call. = FALSE)
      }
    }
  }
  
  if (any(!suppressWarnings(!is.na(as.numeric(as.character(stageframe$size)))))) {
    stop("Function flefko2() requires size to be numeric rather than categorical.",
      call. = FALSE)
  }
  
  melchett <- .sf_reassess(stageframe, supplement, overwrite, repmatrix,
    agemat = FALSE, historical = FALSE, format = 1)
  stageframe <- melchett$stageframe
  repmatrix <- melchett$repmatrix
  ovtable <- melchett$ovtable
  
  if (!all(is.na(overwrite)) | !all(is.na(supplement))) {
    if(any(duplicated(ovtable[,1:3]))) {
      stop("Multiple entries with different values for the same stage transition are not allowed
        in the supplemental or overwrite table. If modifying a historical table to perform an
        ahistorical analysis, then this may be due to different given rates of substitutions
        caused by dropping stage at occasion t-1. Please eliminate duplicate transitions.",
        call. = FALSE)
    }
  }
  
  # Next the data frame for the C++-based matrix populator functions
  allstages <- .theoldpizzle(stageframe, ovtable, repmatrix, finalage = 0,
    format = 1, style = 1, cont = 0)
  
  maxsize <- max(c(allstages$size3, allstages$size2n, allstages$size2o), na.rm = TRUE)
  maxsizeb <- max(c(allstages$sizeb3, allstages$sizeb2n, allstages$sizeb2o), na.rm = TRUE)
  maxsizec <- max(c(allstages$sizec3, allstages$sizec2n, allstages$sizec2o), na.rm = TRUE)
  
  allstages <- allstages[(which(allstages$index321 != -1)),]
  
  # Now we will work up the vital rate models
  if (class(modelsuite) == "lefkoMod") {
    if(is.na(surv_model)) {surv_model <- modelsuite$survival_model}
    if(is.na(obs_model)) {obs_model <- modelsuite$observation_model}
    if(is.na(size_model)) {size_model <- modelsuite$size_model}
    if(is.na(sizeb_model)) {sizeb_model <- modelsuite$sizeb_model}
    if(is.na(sizec_model)) {sizec_model <- modelsuite$sizec_model}
    if(is.na(repst_model)) {repst_model <- modelsuite$repstatus_model}
    if(is.na(fec_model)) {fec_model <- modelsuite$fecundity_model}
    
    if(is.na(jsurv_model)) {jsurv_model <- modelsuite$juv_survival_model}
    if(is.na(jobs_model)) {jobs_model <- modelsuite$juv_observation_model}
    if(is.na(jsize_model)) {jsize_model <- modelsuite$juv_size_model}
    if(is.na(jsizeb_model)) {jsizeb_model <- modelsuite$juv_sizeb_model}
    if(is.na(jsizec_model)) {jsizec_model <- modelsuite$juv_sizec_model}
    if(is.na(jrepst_model)) {jrepst_model <- modelsuite$juv_reproduction_model}
    if(is.na(jmatst_model)) {jmatst_model <- modelsuite$juv_maturity_model}
  }
  
  if (is.na(randomseed)) {
    set.seed(NULL)
  } else {
    set.seed(randomseed)
  }
  
  surv_proxy <- .modelextract(surv_model, paramnames, mainyears, mainpatches, 
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  obs_proxy <- .modelextract(obs_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  sigma <- 0
  rvarssummed <- 0
  sizedist <- 1
  
  size_proxy <- .modelextract(size_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  if (size_proxy$family == "poisson") {
    sizedist <- 0
    if (!all(is.na(size_proxy$variances))) {
      if (is.element("vcov", names(size_proxy$variances))) {
        rvarssummed <- sum(size_proxy$variances[,"vcov"])
      } else if (is.element("cond", names(size_proxy$variances))) {
        rvarssummed <- sum(unlist(lapply(c(1:length(size_proxy$variances[["cond"]])), function(X) {
          stddev.est <- attr(size_proxy$variances[["cond"]][[X]], "stddev")
          return(stddev.est * stddev.est)
        })))
      }
    } else {
      rvarssummed <- 0
    }
  } else if (size_proxy$family == "gaussian") {
    sizedist <- 2
    sigma <- size_proxy$sigma
  } else if (size_proxy$family == "gamma") {
    sizedist <- 3
    
  } else {
    sizedist <- 1
  }
  
  sigmab <- 0
  rvarssummedb <- 0
  sizebdist <- 1
  
  sizeb_proxy <- .modelextract(sizeb_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  if (sizeb_proxy$family == "poisson") {
    sizebdist <- 0
    if (!all(is.na(sizeb_proxy$variances))) {
      if (is.element("vcov", names(sizeb_proxy$variances))) {
        rvarssummedb <- sum(sizeb_proxy$variances[,"vcov"])
      } else if (is.element("cond", names(sizeb_proxy$variances))) {
        rvarssummedb <- sum(unlist(lapply(c(1:length(sizeb_proxy$variances[["cond"]])), function(X) {
          stddev.est <- attr(sizeb_proxy$variances[["cond"]][[X]], "stddev")
          return(stddev.est * stddev.est)
        })))
      }
    } else {
      rvarssummedb <- 0
    }
  } else if (sizeb_proxy$family == "gaussian") {
    sizebdist <- 2
    sigmab <- sizeb_proxy$sigma
  } else if (sizeb_proxy$family == "gamma") {
    sizebdist <- 3
    
  } else {
    sizebdist <- 1
  }
  
  sigmac <- 0
  rvarssummedc <- 0
  sizecdist <- 1
  
  sizec_proxy <- .modelextract(sizec_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  if (sizec_proxy$family == "poisson") {
    sizecdist <- 0
    if (!all(is.na(sizec_proxy$variances))) {
      if (is.element("vcov", names(sizec_proxy$variances))) {
        rvarssummedc <- sum(sizec_proxy$variances[,"vcov"])
      } else if (is.element("cond", names(sizec_proxy$variances))) {
        rvarssummedc <- sum(unlist(lapply(c(1:length(sizec_proxy$variances[["cond"]])), function(X) {
          stddev.est <- attr(sizec_proxy$variances[["cond"]][[X]], "stddev")
          return(stddev.est * stddev.est)
        })))
      }
    } else {
      rvarssummedc <- 0
    }
  } else if (sizec_proxy$family == "gaussian") {
    sizecdist <- 2
    sigmac <- sizec_proxy$sigma
  } else if (sizec_proxy$family == "gamma") {
    sizecdist <- 3
    
  } else {
    sizecdist <- 1
  }
  
  repst_proxy <- .modelextract(repst_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  fec_proxy <- .modelextract(fec_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  if (fec_proxy$family == "poisson") {
    fecdist <- 0
  } else if (fec_proxy$family == "gaussian") {
    fecdist <- 2
  } else if (fec_proxy$family == "gamma") {
    fecdist <- 3
  } else {
    fecdist <- 1
  }
  
  jsurv_proxy <- .modelextract(jsurv_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  jobs_proxy <- .modelextract(jobs_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  jsigma <- 0
  jrvarssummed <- 0
  
  jsize_proxy <- .modelextract(jsize_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  if (jsize_proxy$family == "poisson") {
    if (!all(is.na(jsize_proxy$variances))) {
      if (is.element("vcov", names(jsize_proxy$variances))) {
        jrvarssummed <- sum(jsize_proxy$variances[,"vcov"])
      } else if (is.element("cond", names(jsize_proxy$variances))) {
        jrvarssummed <- sum(unlist(lapply(c(1:length(jsize_proxy$variances[["cond"]])), function(X) {
          stddev.est <- attr(jsize_proxy$variances[["cond"]][[X]], "stddev")
          return(stddev.est * stddev.est)
        })))
      }
    } else {
      jrvarssummed <- 0
    }
  } else if (jsize_proxy$family == "gaussian") {
    jsigma <- jsize_proxy$sigma
  }
  
  jsigmab <- 0
  jrvarssummedb <- 0
  
  jsizeb_proxy <- .modelextract(jsizeb_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  if (jsizeb_proxy$family == "poisson") {
    if (!all(is.na(jsizeb_proxy$variances))) {
      if (is.element("vcov", names(jsizeb_proxy$variances))) {
        jrvarssummedb <- sum(jsizeb_proxy$variances[,"vcov"])
      } else if (is.element("cond", names(jsizeb_proxy$variances))) {
        jrvarssummedb <- sum(unlist(lapply(c(1:length(jsizeb_proxy$variances[["cond"]])), function(X) {
          stddev.est <- attr(jsizeb_proxy$variances[["cond"]][[X]], "stddev")
          return(stddev.est * stddev.est)
        })))
      }
    } else {
      jrvarssummedb <- 0
    }
  } else if (jsizeb_proxy$family == "gaussian") {
    jsigmab <- jsizeb_proxy$sigma
  }
  
  jsigmac <- 0
  jrvarssummedc <- 0
  
  jsizec_proxy <- .modelextract(jsizec_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  if (jsizec_proxy$family == "poisson") {
    if (!all(is.na(jsizec_proxy$variances))) {
      if (is.element("vcov", names(jsizec_proxy$variances))) {
        jrvarssummedc <- sum(jsizec_proxy$variances[,"vcov"])
      } else if (is.element("cond", names(jsizec_proxy$variances))) {
        jrvarssummedc <- sum(unlist(lapply(c(1:length(jsizec_proxy$variances[["cond"]])), function(X) {
          stddev.est <- attr(jsizec_proxy$variances[["cond"]][[X]], "stddev")
          return(stddev.est * stddev.est)
        })))
      }
    } else {
      jrvarssummedc <- 0
    }
  } else if (jsizec_proxy$family == "gaussian") {
    jsigmac <- jsizec_proxy$sigma
  }
  
  jrepst_proxy <- .modelextract(jrepst_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  jmatst_proxy <- .modelextract(jmatst_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  # Next we create a list of pops, patches, and years in order of matrix
  if (!all(is.na(patch))) {
    listofyears <- apply(as.matrix(patch), 1, function(X) {
      output <- cbind.data.frame("1", X, as.matrix(year), stringsAsFactors = FALSE);
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    listofyears <- do.call(rbind.data.frame, listofyears)
    listofyears$poporder <- 1
    listofyears$patchorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainpatches == listofyears$patch[X])})
    listofyears$yearorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainyears == listofyears$year2[X])})
    
  } else {
    listofyears <- cbind.data.frame("1", "1", as.matrix(year), stringsAsFactors = FALSE)
    names(listofyears) <- c("pop", "patch", "year2")
    
    listofyears$poporder <- 1
    listofyears$patchorder <- 1
    listofyears$yearorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainyears == listofyears$year2[X])})
  }
  
  # A few extra tidbits required for the core matrix estimator to work
  yearlist <- split(listofyears, seq(nrow(listofyears)))
  
  # The next line calls the core matrix estimator function
  madsexmadrigal <- lapply(yearlist, .jerzeibalowski, allstages, stageframe, 3,
    surv_proxy, obs_proxy, size_proxy, sizeb_proxy, sizec_proxy, repst_proxy,
    fec_proxy, jsurv_proxy, jobs_proxy, jsize_proxy, jsizeb_proxy, jsizec_proxy,
    jrepst_proxy, jmatst_proxy, f2.inda, f1.inda, f2.indb, f1.indb, f2.indc,
    f1.indc, r2.inda, r1.inda, r2.indb, r1.indb, r2.indc, r1.indc, 
    c(surv_dev, obs_dev, size_dev, sizeb_dev, sizec_dev, repst_dev, fec_dev,
      jsurv_dev, jobs_dev, jsize_dev, jsizeb_dev, jsizec_dev, jrepst_dev, jmatst_dev),
    density, repmod, c(rvarssummed, sigma, rvarssummedb, sigmab, rvarssummedc,
      sigmac, jrvarssummed, jsigma, jrvarssummedb, jsigmab, jrvarssummedc, jsigmac),
    maxsize, maxsizeb, maxsizec, 0, sizedist, sizebdist, sizecdist, fecdist,
    negfec, exp_tol, theta_tol, ipm_method)
  
  a_list <- lapply(madsexmadrigal, function(X) {X$A})
  u_list <- lapply(madsexmadrigal, function(X) {X$U})
  f_list <- lapply(madsexmadrigal, function(X) {X$F})
  
  if (err_check) {
    out_list <- lapply(madsexmadrigal, function(X) {
        aaa <- X$out
        colnames(aaa) <- c("surv", "obs", "repst", "sizea", "sizeb", "sizec", "matst")
        
        return(aaa)
      }
    )
  }
  
  ahstages <- stageframe[1:(dim(stageframe)[1] - 1),]
  
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
  
  if (!err_check) {
    output <- list(A = a_list, U = u_list, F = f_list, hstages = NA,
      agestages = NA, ahstages = ahstages, labels = listofyears[,c(1:3)],
      matrixqc = qcoutput1, modelqc = qcoutput2)
  } else {
    output <- list(A = a_list, U = u_list, F = f_list, hstages = NA,
      agestages = NA, ahstages = ahstages, labels = listofyears[,c(1:3)],
      matrixqc = qcoutput1, modelqc = qcoutput2, prob_out = out_list)
  }
  
  class(output) <- "lefkoMat"
  
  return(output)
}

#' Reduce Matrix Dimensions By Eliminating Empty Stages
#' 
#' \code{.reducer2()} identifies empty stages in a set of ahistorical matrices
#' and removes them from all matrices. It also removes the associated rows in
#' the associated \code{ahstages} or \code{agestages} object. It is used within
#' \code{\link{flefko2}()}, \code{\link{aflefko2}()}, and
#' \code{\link{rlefko2}()}.
#' 
#' @param A List of population projection matrices, from a \code{lefkoMat}
#' object.
#' @param U List of surviva-transition matrices corresponding to \code{A}.
#' @param F List of fecundity matrices corresponding to \code{A}.
#' @param ahstages Data frame giving the names and identities of ahistorical 
#' stages used to create matrices.
#' 
#' @return Returns a list of reduced \code{A}, \code{U}, and \code{F} matrices,
#' plus the reduced \code{ahstages} object. Note that this can also work on
#' \code{agestages}, if passed instead of \code{ahstages}.
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

#' Create Raw Historical Matrix Projection Model
#' 
#' Function \code{rlefko3()} returns raw historical MPMs, including the
#' associated component transition and fecundity matrices, data frames
#' describing the ahistorical stages used and the historical paired stages, and
#' a data frame describing the population, patch, and occasion time associated
#' with each matrix.
#' 
#' @param data A vertical demographic data frame, with variables corresponding
#' to the naming conventions in \code{\link{verticalize3}()}.
#' @param stageframe A stageframe object that includes information on the size,
#' observation status, propagule status, immaturity status, and maturity status
#' of each ahistorical stage. Should also incorporate bin widths if size is
#' continuous.
#' @param year A variable corresponding to observation occasion, or a set of
#' such values, given in values associated with the \code{year} term used in
#' vital rate model development. Can also equal \code{"all"}, in which case
#' matrices will be estimated for all occasions. Defaults to \code{"all"}.
#' @param pop A variable designating which populations will have matrices
#' estimated. Should be set to specific population names, or to \code{"all"} if
#' all populations should have matrices estimated.
#' @param patch A variable designating which patches or subpopulations will have
#' matrices estimated. Should be set to specific patch names, or to \code{"all"}
#' if matrices should be estimated for all patches. Defaults to \code{"all"}.
#' @param censor If \code{TRUE}, then data will be removed according to the
#' variable set in \code{censorcol}, such that only data with censor values
#' equal to \code{censorkeep} will remain. Defaults to \code{FALSE}.
#' @param stages An optional vector denoting the names of the variables within
#' the main vertical dataset coding for the stages of each individual in
#' occasions \emph{t}+1, \emph{t}, and \emph{t}-1. The names of stages in these
#' variables should match those used in the \code{stageframe} exactly. If left
#' blank, then \code{rlefko3()} will attempt to infer stages by matching values
#' of \code{alive}, \code{size}, \code{repst}, and \code{matst} to
#' characteristics noted in the associated \code{stageframe}.
#' @param alive A vector of names of binomial variables corresponding to status
#' as alive (\code{1}) or dead (\code{0}) in occasions \emph{t}+1, \emph{t}, and
#' \emph{t}-1, respectively.
#' @param size A vector of names of variables coding the primary size variable
#' in occasions \emph{t}+1, \emph{t}, and \emph{t}-1, respectively. Defaults to 
#' \code{c("sizea3", "sizea2", "sizea1")}.
#' @param sizeb A vector of names of variables coding the secondary size
#' variable in occasions \emph{t}+1, \emph{t}, and \emph{t}-1, respectively.
#' Defaults to \code{c(NA, NA, NA)}.
#' @param sizec A vector of names of variables coding the tertiary size
#' variable in occasions \emph{t}+1, \emph{t}, and \emph{t}-1, respectively.
#' Defaults to \code{c(NA, NA, NA)}.
#' @param repst A vector of names of variables coding reproductive status in
#' occasions \emph{t}+1, \emph{t}, and \emph{t}-1, respectively. Defaults to 
#' \code{c("repstatus3", "repstatus2", "repstatus1")}. Must be supplied if
#' \code{stages} is not provided.
#' @param matst A vector of names of variables coding maturity status in
#' occasions \emph{t}+1, \emph{t}, and \emph{t}-1, respectively. Defaults to
#' \code{c("matstatus3", "matstatus2", "matstatus1")}. Must be supplied if
#' \code{stages} is not provided.
#' @param fec A vector of names of variables coding fecundity in occasions
#' \emph{t}+1, \emph{t}, and \emph{t}-1, respectively. Defaults to
#' \code{c("feca3", "feca2", "feca1")}.
#' @param supplement An optional data frame of class \code{lefkoSD} that
#' provides supplemental data that should be incorporated into the MPM. Three
#' kinds of data may be integrated this way: transitions to be estimated via the
#' use of proxy transitions, transition overwrites from the literature or
#' supplemental studies, and transition multipliers for fecundity. This data
#' frame should be produced using the \code{\link{supplemental}()} function.
#' Should be used in place of or in addition to an overwrite table (see 
#' \code{overwrite} below) and a reproduction matrix (see \code{repmatrix}
#' below).
#' @param repmatrix An optional reproduction matrix. This matrix is composed
#' mostly of 0s, with non-zero entries acting as element identifiers and
#' multipliers for fecundity (with 1 equaling full fecundity). If left blank,
#' and no \code{supplement} is provided, then \code{rlefko3()} will assume that
#' all stages marked as reproductive produce offspring at 1x that of estimated
#' fecundity, and that offspring production will yield the first stage noted as
#' propagule or immature. To prevent this behavior, input just \code{0}, which
#' will result in fecundity being estimated only for transitions noted in
#' \code{supplement} above. May be the dimensions of either a historical or an
#' ahistorical matrix. If the latter, then all stages will be used in occasion
#' \emph{t}-1 for each suggested ahistorical transition.
#' @param overwrite An optional data frame developed with the
#' \code{\link{overwrite}()} function describing transitions to be overwritten
#' either with given values or with other estimated transitions. Note that this
#' function supplements overwrite data provided in \code{supplement}.
#' @param yearcol The variable name or column number corresponding to occasion
#' \emph{t} in the dataset.
#' @param popcol The variable name or column number corresponding to the
#' identity of the population.
#' @param patchcol The variable name or column number corresponding to patch in 
#' the dataset.
#' @param indivcol The variable name or column number coding individual
#' identity.
#' @param censorcol The variable name or column number denoting the censor
#' status. Only needed if \code{censor = TRUE}.
#' @param censorkeep The value of the censor variable denoting data elements to
#' keep. Defaults to \code{0}.
#' @param format A string indicating whether to estimate matrices in
#' \code{ehrlen} format or \code{deVries} format. The latter adds one unborn
#' prior stage to account for the prior state of newborns. Defaults to
#' \code{ehrlen} format.
#' @param reduce A logical value denoting whether to remove historical stages
#' associated exclusively with zero transitions. These are removed only if the
#' respective row and column sums in ALL matrices estimated equal 0. Defaults to
#' \code{FALSE}.
#' @param err_check A logical value indicating whether to append extra
#' information used in matrix calculation within the output list. Used for
#' development debugging purposes. Defaults to \code{FALSE}.
#'
#' @return If all inputs are properly formatted, then this function will return
#' an object of class \code{lefkoMat}, which is a list that holds the matrix
#' projection model and all of its metadata. Its structure is a list with the
#' following elements:
#'
#' \item{A}{A list of full projection matrices in order of sorted populations,
#' patches, and occasions. All matrices output in the \code{matrix} class.}
#' \item{U}{A list of survival transition matrices sorted as in \code{A}. All 
#' matrices output in the \code{matrix} class.}
#' \item{F}{A list of fecundity matrices sorted as in \code{A}. All matrices 
#' output in the \code{matrix} class.}
#' \item{hstages}{A data frame matrix showing the pairing of ahistorical stages
#' used to create historical stage pairs.}
#' \item{agestages}{A data frame showing age-stage pairs. In this function, it
#' is set to NA. Only used in output to function \code{aflefko2}().}
#' \item{ahstages}{A data frame detailing the characteristics of associated
#' ahistorical stages, in the form of a modified stageframe that includes
#' status as an entry stage through reproduction.}
#' \item{labels}{A data frame giving the population, patch, and year of each 
#' matrix in order.}
#' \item{matrixqc}{A short vector describing the number of non-zero elements in
#' \code{U} and \code{F} matrices, and the number of annual matrices.}
#' \item{dataqc}{A vector showing the numbers of individuals and rows in the
#' vertical dataset used as input.}
#'
#' @section Notes:
#' The default behavior of this function is to estimate fecundity with regards
#' to transitions specified via associated fecundity multipliers in either
#' \code{supplement} or \code{repmatrix}. If both of these fields are left
#' empty, then fecundity will be estimated at full for all transitions leading
#' from reproductive stages to immature and propagule stages. However, if a
#' \code{supplement} is provided and a \code{repmatrix} is not, or if
#' \code{repmatrix} is set to 0, then only fecundity transitions noted in the
#' supplement will be set to non-zero values. To use the default behavior of
#' setting all reproductive stages to reproduce at full fecundity into immature
#' and propagule stages but incorporate given or proxy survival transitions,
#' input those given and proxy transitions through the \code{overwrite} option.
#' 
#' The reproduction matrix (field \code{repmatrix}) may be supplied as either
#' historical or ahistorical. If provided as ahistorical, then \code{flefko3()}
#' will assume that all historical transitions involving stages noted for
#' occasions \emph{t} and \emph{t}+1 should be set to the respective fecundity
#' multipliers noted.
#' 
#' Users may at times wish to estimate MPMs using a dataset incorporating
#' multiple patches or subpopulations. Should the aim of analysis be a general
#' MPM that does not distinguish these patches or subpopulations, the
#' \code{patchcol} variable should be left to \code{NA}, which is the default.
#' Otherwise the variable identifying patch needs to be named.
#'
#' Input options including multiple variable names must be entered in the order
#' of variables in occasion \emph{t}+1, \emph{t}, and \emph{t}-1. Rearranging
#' the order WILL lead to erroneous calculations, and may lead to
#' fatal errors.
#'
#' Although this function is capable of assigning stages given an input
#' stageframe, it lacks the power of \code{\link{verticalize3}()} and
#' \code{\link{historicalize3}()} in this regard. Users are strongly
#' encouraged to use the latter two functions for stage assignment.
#' 
#' @seealso \code{\link{flefko3}()}
#' @seealso \code{\link{flefko2}()}
#' @seealso \code{\link{aflefko2}()}
#' @seealso \code{\link{fleslie}()}
#' @seealso \code{\link{rlefko2}()}
#' @seealso \code{\link{rleslie}()}
#' 
#' @examples
#' # Lathyrus example
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, 
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector, 
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec, 
#'   propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, 
#'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9, 
#'   juvcol = "Seedling1988", sizeacol = "Volume88", repstracol = "FCODE88", 
#'   fecacol = "Intactseed88", deadacol = "Dead1988", nonobsacol = "Dormant1988", 
#'   stageassign = lathframe, stagesize = "sizea", censorcol = "Missing1988", 
#'   censorkeep = NA, censor = TRUE)
#' 
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl", "mat"),
#'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep", "Sdl"),
#'   stage1 = c("Sd", "rep", "Sd", "rep", "npr", "npr", "Sd"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "mat"),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "Sdl"),
#'   eststage1 = c(NA, NA, NA, NA, NA, NA, "NotAlive"),
#'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054, NA),
#'   type = c(1, 1, 1, 1, 3, 3, 1), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
#'   stageframe = lathframe, historical = TRUE)
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = "all", 
#'   stages = c("stage3", "stage2", "stage1"), supplement = lathsupp3,
#'   yearcol = "year2", indivcol = "individ")
#' 
#' summary(ehrlen3)
#' 
#' # Cypripedium example
#' rm(list=ls(all=TRUE))
#' data(cypdata)
#' 
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
#'   "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   propstatus = propvector, immstatus = immvector, indataset = indataset,
#'   binhalfwidth = binvec)
#' 
#' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
#'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
#'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
#'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
#'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
#'   NRasRep = TRUE)
#' 
#' cypsupp3r <- supplemental(stage3 = c("SD", "SD", "P1", "P1", "P2", "P3", "SL",
#'     "D", "XSm", "Sm", "D", "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL",
#'     "SL", "SL", "rep", "rep"),
#'   stage1 = c("SD", "rep", "SD", "rep", "SD", "P1", "P2", "P3", "P3", "P3",
#'     "SL", "SL", "SL", "mat", "mat"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", "D", "XSm", "Sm",
#'     NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm", "XSm",
#'     "XSm", NA, NA),
#'   eststage1 = c(NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm", "XSm",
#'     "XSm", NA, NA),
#'   givenrate = c(0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.25, NA, NA, NA, NA, NA, NA,
#'     NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   type_t12 = c(1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
#'   stageframe = cypframe_raw, historical = TRUE)
#' 
#' cypmatrix3r <- rlefko3(data = cypraw_v1, stageframe = cypframe_raw,
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added", "size1added"),
#'   supplement = cypsupp3r, yearcol = "year2", patchcol = "patchid",
#'   indivcol = "individ")
#' 
#' summary(cypmatrix3r)
#' 
#' @export
rlefko3 <- function(data, stageframe, year = "all", pop = NA, patch = NA,
  censor = FALSE, stages = NA, alive = c("alive3", "alive2", "alive1"),
  size = c("sizea3", "sizea2", "sizea1"), sizeb = c(NA, NA, NA),
  sizec = c(NA, NA, NA), repst = c("repstatus3", "repstatus2", "repstatus1"),
  matst = c("matstatus3", "matstatus2", "matstatus1"),
  fec = c("feca3", "feca2", "feca1"), supplement = NULL, repmatrix = NULL,
  overwrite = NULL, yearcol = NA, popcol = NA, patchcol = NA, indivcol = NA,
  censorcol = NA, censorkeep = 0, format = "ehrlen", reduce = FALSE,
  err_check = FALSE) {
  
  instageframe <- tocensor <- indataset <- alive2 <- popused <- patchused <- yearused <- NULL
  
  sizeb_used <- 0
  sizec_used <- 0
  
  if (tolower(format) == "ehrlen") {
    format_int <- 1
  } else if (tolower(format) == "devries") {
    format_int <- 2
  } else {
    stop("The format parameter must be set to either 'ehrlen' or 'deVries'.",
      call. = FALSE)
  }
  
  if (all(is.na(data))) {
    stop("Need original vertical dataset to proceed.", call. = FALSE)
  }
  
  if (!any(class(data) == "data.frame")) {
    stop("Need original vertical dataset to proceed. This dataset must be in
      historical vertical format.", call. = FALSE)
  }
  
  if (!any(class(data) == "hfvdata")) {
    warning("Dataset used as input is not of class hfvdata. Will assume that the
      dataset has been formatted equivalently.", call. = FALSE)
  }
  
  stageframe_vars <- c("stage", "size", "size_b", "size_c", "min_age", "max_age",
    "repstatus", "obsstatus", "propstatus", "immstatus", "matstatus", "indataset",
    "binhalfwidth_raw", "sizebin_min", "sizebin_max", "sizebin_center",
    "sizebin_width", "binhalfwidthb_raw", "sizebinb_min", "sizebinb_max",
    "sizebinb_center", "sizebinb_width", "binhalfwidthc_raw", "sizebinc_min",
    "sizebinc_max", "sizebinc_center", "sizebinc_width", "group", "comments")
  if (any(!is.element(names(stageframe), stageframe_vars))) {
    stop("Please use properly formatted stageframe as input.", call. = FALSE)
  }
  
  if (all(is.na(stages))) {
    if ((length(alive) != 3)) {
      stop("This function requires stage information for each of occasions t+1,
        t, and t-1. In the absence of stage columns in the dataset, it requires
        the input of data for living/dead status, size, reproductive status, and
        maturity status, for each of occasions t+1, t, and t-1.", call. = FALSE)
    }
    if ((length(size) != 3)) {
      stop("This function requires stage information for each of occasions t+1,
        t, and t-1. In the absence of stage columns in the dataset, it requires
        the input of data for living/dead status, size, reproductive status, and
        maturity status, for each of occasions t+1, t, and t-1.", call. = FALSE)
    }
    if (!all(is.na(repst))) {
      if ((length(repst) != 3)) {
        stop("This function requires stage information for each of occasions t+1,
          t, and t-1. In the absence of stage columns in the dataset, it requires
          the input of data for living/dead status, size, reproductive status, and
          maturity status, for each of occasions t+1, t, and t-1.", call. = FALSE)
      }
    }   
    if (!all(is.na(matst))) {
      if ((length(matst) != 3)) {
        stop("This function requires stage information for each of occasions t+1,
          t, and t-1. In the absence of stage columns in the dataset, it requires
          the input of data for living/dead status, size, reproductive status, and
          maturity status, for each of occasions t+1, t, and t-1.", call. = FALSE)
      }
    }   
  } else if (length(stages) != 3) {
    stop("This function requires stage information for each of occasions t+1, t,
      and t-1.", call. = FALSE)
  }
  
  if ((length(fec) != 3)) {
    stop("This function requires three variables for fecundity, for each of
      occasions t+1, t, and t-1.", call. = FALSE)
  }
  
  if (is.character(yearcol)) {
    choicevar <- which(names(data) == yearcol);
    mainyears <- sort(unique(data[,choicevar]))[-1] # Occasion 1 is unusable, so removed
  } else if (is.numeric(yearcol)) {
    mainyears <- sort(unique(data[, yearcol]))[-1]
  } else {
    stop("Need appropriate year column designation.", call. = FALSE)
  }
  
  if (any(is.character(year))) {
    if (is.element("all", tolower(year))) {
      year <- mainyears
    } else {
      stop("Year designation not recognized.", call. = FALSE)
    }
  }
  
  if (length(year) == 0 | all(is.na(year) == TRUE) | any(is.na(year))) {
    stop("This function cannot proceed without being given a specific year, or a
      suite of years. NA entries are not allowed.", call. = FALSE)
  }
  
  if (!all(is.element(year, mainyears))) {
    stop("Dataset does not contain one or more of the requested years. Note that
      matrices cannot be made for the first year in a historical dataset.",
      call. = FALSE)
  }
  
  if (censor == TRUE) {
    if(all(is.na(censorcol)) == TRUE) {
      stop("Cannot censor the data without a proper censor variable.", call. = FALSE)
    }
    
    if (all(is.character(censorcol))) {
      if (!all(is.element(censorcol, names(data)))) {
        stop("Censor variable names input for censorcol do not match any
          variable names in the dataset.", call. = FALSE)
      }
    }
    
    censorcolsonly <- data[,censorcol]
    sleeplessnights <- apply(as.matrix(c(1:dim(censorcolsonly)[1])), 1, function(X) {
      crazyvec <- if(is.element(censorkeep, censorcolsonly[X,])) {
        return(X);
      } else {
        return(NA);
      }
    })
    sleeplessnights <- sleeplessnights[!is.na(sleeplessnights)]
    
    data <- data[sleeplessnights,]
  }
  
  if (!all(is.na(pop)) & !all(is.na(patch))) {
    if (is.na(popcol) | is.na(patchcol)) {
      stop("Need population and patch designation variables to proceed.", 
        call. = FALSE)
    }
    
    if (is.element("all", tolower(pop))) {
      if (is.character(popcol)) {popcol <- which(names(data) == popcol)}
      
      pops <- unique(data[,popcol])
    } else {pops <- pop}
    
    if (is.element("all", tolower(patch))) {
      if (is.character(patchcol)) {patchcol <- which(names(data) == patchcol)}
      
      listofpatches <- apply(as.matrix(pops), 1, function(X) {
        patchfolly <- subset(data, popcol == X);
        output <- cbind.data.frame(X, sort(unique(patchfolly[,patchcol])),
          stringsAsFactors = FALSE);
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
      output <- cbind.data.frame(X[1], X[2], sort(unique(checkyrdata[,yearcol]))[-1],
        stringsAsFactors = FALSE);
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
  } else if (all(is.na(pop)) & !all(is.na(patch))) {
    if (!is.na(popcol)) {
      popcol <- NA
    }
    
    if (is.na(patchcol)) {
      stop("Need patch designation variable to proceed.", call. = FALSE)
    }
    
    if (is.element("all", tolower(patch))) {
      if (is.character(patchcol)) {patchcol <- which(names(data) == patchcol)}
      
      patches <- sort(unique(data[,patchcol]))
    } else {patches <- patch}
    
    listofyears <- apply(as.matrix(patches), 1, function(X) {
      checkyrdata <- subset(data, patchcol = X);
      output <- cbind.data.frame("1", X, sort(unique(checkyrdata[,yearcol]))[-1],
        stringsAsFactors = FALSE);
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    if (length(listofyears) > 1) {
      listofyears <- do.call(rbind.data.frame, listofyears)
    } else {
      listofyears <- listofyears[[1]]
    }
  } else if (!all(is.na(pop)) & all(is.na(patch))) {
    if (is.na(popcol)) {
      stop("Need population designation variable to proceed.", call. = FALSE)
    }
    
    if (!is.na(patchcol)) {
      patchcol <- NA
    }
    
    if (is.element("all", tolower(pop))) {
      if (is.character(popcol)) {popcol <- which(names(data) == popcol)}
      
      pops <- unique(data[,popcol])
    } else {pops <- pop}
    
    listofyears <- apply(as.matrix(pops), 1, function(X) {
      checkyrdata <- subset(data, popcol = X);
      output <- cbind.data.frame(X, "1", sort(unique(checkyrdata[,yearcol]))[-1],
        stringsAsFactors = FALSE);
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    if (length(listofyears) > 1) {
      listofyears <- do.call(rbind.data.frame, listofyears)
    } else {
      listofyears <- listofyears[[1]]
    }
  } else if (all(is.na(pop)) & all(is.na(patch))) {
    if (!is.na(popcol)) {
      popcol <- NA
    }
    if (!is.na(patchcol)) {
      patchcol <- NA
    }
    
    listofyears <- cbind.data.frame("1", "1", year, stringsAsFactors = FALSE)
    names(listofyears) <- c("pop", "patch", "year2")
  }
  
  identifiedyearrows <- which(is.element(listofyears$year2, year))
  if (length(identifiedyearrows) == 0) {
    stop("Cannot recognize input year(s)", call. = FALSE)
  } else {
    listofyears <- listofyears[identifiedyearrows,]
  }
  yearlist <- split(listofyears, seq(nrow(listofyears)))
  
  stagenum_init <- dim(stageframe)[1]
  if (!all(is.null(repmatrix))) {
    if (any(class(repmatrix) == "matrix")) {
      if (dim(repmatrix)[1] != stagenum_init & dim(repmatrix)[1] != stagenum_init^2) {
        stop("The repmatrix must be a square matrix with dimensions equal to the
          number of stages in the stageframe, or the square thereof.",
          call. = FALSE)
      }
      
      if (dim(repmatrix)[2] != stagenum_init & dim(repmatrix)[2] != stagenum_init^2) {
        stop("The repmatrix must be a square matrix with dimensions equal to the
          number of stages in the stageframe, or the square thereof.",
          call. = FALSE)
      }
    }
  }
  
  melchett <- .sf_reassess(stageframe, supplement, overwrite, repmatrix,
    agemat = FALSE, historical = TRUE, format = format_int)
  stageframe <- melchett$stageframe
  repmatrix <- melchett$repmatrix
  ovtable <- melchett$ovtable
  
  if (!all(is.na(overwrite)) | !all(is.na(supplement))) {
    
    if(any(duplicated(ovtable[,1:3]))) {
      stop("Multiple entries with different values for the same stage transition
        are not allowed in the supplemental or overwrite table. If performing
        a historical analysis, then please remember that all stages must be
        clearly defined in all three times, including time t-1. Please eliminate
        duplicate transitions.", call. = FALSE)
    }
  }
  
  data$alive1 <- data[,which(names(data) == alive[3])]
  data$alive2 <- data[,which(names(data) == alive[2])]
  data$alive3 <- data[,which(names(data) == alive[1])]
  
  if (all(is.na(stages))) {
    if (length(size) > 2) {
      size1met <- which(names(data) == size[3])
      size2met <- which(names(data) == size[2])
      size3met <- which(names(data) == size[1])
      
      if (length(size1met) == 1 & length(size2met) == 1 & length(size3met) == 1) {
        data$usedsize1 <- data[,size1met]
        data$usedsize2 <- data[,size2met]
        data$usedsize3 <- data[,size3met]
      } else {
        stop("Entered size variable names do not strictly correspond to single
          variables in the dataset.", call. = FALSE)
      }
      
      if (!all(is.na(sizeb)) & length(sizeb) > 2) {
        size1met <- which(names(data) == sizeb[3])
        size2met <- which(names(data) == sizeb[2])
        size3met <- which(names(data) == sizeb[1])
        
        if (length(size1met) == 1 & length(size2met) == 1 & length(size3met) == 1) {
          sizeb_used <- 1
          data$usedsize1b <- data[,size1met]
          data$usedsize2b <- data[,size2met]
          data$usedsize3b <- data[,size3met]
        } else {
          stop("Entered sizeb variable names do not strictly correspond to single
          variables in the dataset, or fewer than 3 variable names have been
          entered.", call. = FALSE)
        }
      }
      
      if (!all(is.na(sizec)) & length(sizec) > 2) {
        size1met <- which(names(data) == sizec[3])
        size2met <- which(names(data) == sizec[2])
        size3met <- which(names(data) == sizec[1])
        
        if (length(size1met) == 1 & length(size2met) == 1 & length(size3met) == 1) {
          sizec_used <- 1
          data$usedsize1c <- data[,size1met]
          data$usedsize2c <- data[,size2met]
          data$usedsize3c <- data[,size3met]
        } else {
          stop("Entered sizec variable names do not strictly correspond to single
          variables in the dataset, or fewer than 3 variable names have been
          entered.", call. = FALSE)
        }
      }
    } else {
      warning("Without stage columns, rlefko3() requires size variables as input.
        Failure to include size variables may lead to odd results.", call. = FALSE)
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
      mainstages <- intersect(which(stageframe$bin_size_min < data$usedsize1[X]), 
        which(stageframe$bin_size_max >= data$usedsize1[X]))
      if (sizeb_used == 1) {
        if (is.na(data$usedsize1b[X])) {
          data$usedsize1b[X] <- 0
        }
        mainstages <- intersect(which(instageframe$sizebinb_min < data$usedsize1b[X]), mainstages)
        mainstages <- intersect(which(instageframe$sizebinb_max >= data$usedsize1b[X]), mainstages)
      }
      if (sizec_used == 1) {
        if (is.na(data$usedsize1c[X])) {
          data$usedsize1c[X] <- 0
        }
        mainstages <- intersect(which(instageframe$sizebinc_min < data$usedsize1c[X]), mainstages)
        mainstages <- intersect(which(instageframe$sizebinc_max >= data$usedsize1c[X]), mainstages)
      }
      
      jmstages <- which(stageframe$immstatus == (1 - data$usedmatstatus1[X]))
      obsstages <- which(stageframe$obsstatus == data$obsstatus1[X])
      repstages <- which(stageframe$repstatus == data$repstatus1[X])
      alivestage1 <- which(stageframe$alive == data$alive1[X])
      
      choicestage <- intersect(intersect(intersect(mainstages, jmstages),
          intersect(obsstages, repstages)), alivestage1)
      
      if (length(choicestage) == 0) {
        choicestage <- which(stageframe$stage_id == max(stageframe$stage_id))
        if (data$alive1[X] != 0) {
          stop("Stage characteristics mismatch dataset. Consider using the stages
            option, particularly if the vertical file was created with NRasRep = TRUE
            in verticalize3() or historicalize3().", call. = FALSE)
        }
      }
      
      return(as.character(stageframe$stage[choicestage]))
    })
    
    data$usedstage2 <- apply(as.matrix(c(1:dim(data)[1])), 1, function(X) {
      if (is.na(data$usedsize2[X])) {
        data$usedsize2[X] <- 0
      }
      mainstages <- intersect(which(stageframe$bin_size_min < data$usedsize2[X]), 
        which(stageframe$bin_size_max >= data$usedsize2[X]))
      if (sizeb_used == 1) {
        if (is.na(data$usedsize2b[X])) {
          data$usedsize2b[X] <- 0
        }
        mainstages <- intersect(which(instageframe$sizebinb_min < data$usedsize2b[X]), mainstages)
        mainstages <- intersect(which(instageframe$sizebinb_max >= data$usedsize2b[X]), mainstages)
      }
      if (sizec_used == 1) {
        if (is.na(data$usedsize2c[X])) {
          data$usedsize2c[X] <- 0
        }
        mainstages <- intersect(which(instageframe$sizebinc_min < data$usedsize2c[X]), mainstages)
        mainstages <- intersect(which(instageframe$sizebinc_max >= data$usedsize2c[X]), mainstages)
      }
      jmstages <- which(stageframe$immstatus == (1 - data$usedmatstatus2[X]))
      obsstages <- which(stageframe$obsstatus == data$obsstatus2[X])
      repstages <- which(stageframe$repstatus == data$repstatus2[X])
      
      choicestage <- intersect(intersect(mainstages, jmstages),
          intersect(obsstages, repstages))
      
      if (length(choicestage) == 0) {
        choicestage <- which(stageframe$stage_id == max(stageframe$stage_id))
        if (data$alive2[X] != 0) {
          stop("Stage characteristics mismatch dataset. Consider using the stages
            option, particularly if the vertical file was created with NRasRep = TRUE
            in verticalize3() or historicalize3().", call. = FALSE)
        }
      }
      
      return(as.character(stageframe$stage[choicestage]))
    })
    
    data$usedstage3 <- apply(as.matrix(c(1:dim(data)[1])), 1, function(X) {
      if (is.na(data$usedsize3[X])) {
        data$usedsize3[X] <- 0
      }
      mainstages <- intersect(which(stageframe$bin_size_min < data$usedsize3[X]), 
        which(stageframe$bin_size_max >= data$usedsize3[X]))
      if (sizeb_used == 1) {
        if (is.na(data$usedsize3b[X])) {
          data$usedsize3b[X] <- 0
        }
        mainstages <- intersect(which(instageframe$sizebinb_min < data$usedsize3b[X]), mainstages)
        mainstages <- intersect(which(instageframe$sizebinb_max >= data$usedsize3b[X]), mainstages)
      }
      if (sizec_used == 1) {
        if (is.na(data$usedsize3c[X])) {
          data$usedsize3c[X] <- 0
        }
        mainstages <- intersect(which(instageframe$sizebinc_min < data$usedsize3c[X]), mainstages)
        mainstages <- intersect(which(instageframe$sizebinc_max >= data$usedsize3c[X]), mainstages)
      }
      jmstages <- which(stageframe$immstatus == (1 - data$usedmatstatus3[X]))
      obsstages <- which(stageframe$obsstatus == data$obsstatus3[X])
      repstages <- which(stageframe$repstatus == data$repstatus3[X])
      alivestage3 <- which(stageframe$alive == data$alive3[X])
      
      choicestage <- intersect(intersect(intersect(mainstages, jmstages),
        intersect(obsstages, repstages)), alivestage3)
      
      if (length(choicestage) == 0) {
        choicestage <- which(instageframe$stage_id == max(instageframe$stage_id))
        if (data$alive3[X] != 0) {
          stop("Stage characteristics mismatch dataset. Consider using the stages
            option, particularly if the vertical file was created with NRasRep = TRUE
            in verticalize3() or historicalize3().", call. = FALSE)
        }
      }
      
      return(as.character(stageframe$stage[choicestage]))
    })
  } else if (length(stages) > 2) {
    if (is.numeric(stages[2])) {
      data$usedstage1 <- data[, stages[3]]
      data$usedstage2 <- data[, stages[2]]
      data$usedstage3 <- data[, stages[1]]
      
      data$usedstage1 <- as.character(data$usedstage1)
      data$usedstage2 <- as.character(data$usedstage2)
      data$usedstage3 <- as.character(data$usedstage3)
      
      if (is.element("NotAlive", unique(data$usedstage1))) {
        data$usedstage1[which(data$usedstage1 == "NotAlive")] <- "Dead"
      }
      if (is.element("NotAlive", unique(data$usedstage3))) {
        data$usedstage3[which(data$usedstage3 == "NotAlive")] <- "Dead"
      }
    } else {
      data$usedstage1 <- data[,which(names(data) == stages[3])]
      data$usedstage2 <- data[,which(names(data) == stages[2])]
      data$usedstage3 <- data[,which(names(data) == stages[1])]
      
      data$usedstage1 <- as.character(data$usedstage1)
      data$usedstage2 <- as.character(data$usedstage2)
      data$usedstage3 <- as.character(data$usedstage3)
      
      if (is.element("NotAlive", unique(data$usedstage1))) {
        data$usedstage1[which(data$usedstage1 == "NotAlive")] <- "Dead"
      }
      if (is.element("NotAlive", unique(data$usedstage3))) {
        data$usedstage3[which(data$usedstage3 == "NotAlive")] <- "Dead"
      }
    }
    stages.used <- sort(unique(c(data$usedstage2, data$usedstage3)))
    
    if (length(setdiff(stages.used, stageframe$stage)) > 0 & !is.element("NotAlive", stages.used)) {
      stop("Some stages in dataset do not match those detailed in the input stageframe.",
        call. = FALSE)
    }
  }
  
  if (length(fec) > 2) {
    data$usedfec1 <- data[,which(names(data) == fec[3])]
    data$usedfec2 <- data[,which(names(data) == fec[2])]
    data$usedfec3 <- data[,which(names(data) == fec[1])]
    
    data$usedfec1[which(is.na(data$usedfec1))] <- 0
    data$usedfec2[which(is.na(data$usedfec2))] <- 0
    data$usedfec3[which(is.na(data$usedfec3))] <- 0
  } else {
    warning("Function rlefko3() requires 3 fecundity variables, for times t+1, t,
      and t-1. Failure to include fecundity variables leads to matrices composed
      only of survival transitions.", call. = FALSE)
  } 
  
  # Here we search for NotAlive entries in the dataset, and then alter the dataset
  # based on eststage entries in the ovtable. Finally, we remove the NotAlive
  # entries from the ovtable table
  flubbleindices <- which(tolower(ovtable$eststage1) == "notalive")
  if (length(flubbleindices) > 0) {
    flubble <- ovtable[flubbleindices,]
    for (i in c(1:dim(flubble)[1])) {
      datamatch_t1 <- which(tolower(data$stage1) == "notalive")
      datamatch_t2 <- which(data$usedstage2 == flubble$eststage2[i])
      datamatch_t3 <- which(data$usedstage3 == flubble$eststage3[i])
      
      finalshowdown <- intersect(intersect(datamatch_t1, datamatch_t2), datamatch_t3)
      
      if (length(finalshowdown) > 0) {
        data$usedstage1[finalshowdown] <- flubble$stage1[i]
        data$usedstage2[finalshowdown] <- flubble$stage2[i]
        data$usedstage3[finalshowdown] <- flubble$stage3[i]
      }
    }
    ovtable <- ovtable[-flubbleindices,]
  }
  
  # This section creates stageexpansion9, which is a data frame that holds values for stage transitions from paired stages
  # in occasions t and t-1 to paired stages in occasions t and t+1
  stageexpansion9 <- .theoldpizzle(stageframe, ovtable, repmatrix, finalage = 0,
    format = format_int, style = 0, cont = 0)
  
  #Here we reformat the repmatrix if deVries format is chosen
  harpoon <- if (format_int == 2) {
    prior_reps <- colSums(repmatrix)
    prior_reps[which(prior_reps > 0)] <- 1
    
    cbind(rbind(repmatrix, prior_reps, 0), 0, 0)
  } else {
    harpoon <- cbind(rbind(repmatrix, 0), 0)
  }
  
  # Stageexpansion3 is a dataframe created to hold values for paired stages in occasions t and t-1 only
  stageexpansion3 <- cbind.data.frame(expand.grid(size3 = stageframe$sizebin_center, 
    size2n = stageframe$sizebin_center), expand.grid(sizeb3 = stageframe$sizebinb_center, 
    sizeb2n = stageframe$sizebinb_center), expand.grid(sizec3 = stageframe$sizebinc_center, 
    sizec2n = stageframe$sizebinc_center), expand.grid(rep3 = stageframe$repstatus, 
    rep2n = stageframe$repstatus), expand.grid(indata3 = stageframe$indataset, 
    indata2n = stageframe$indataset), expand.grid(stage3 = stageframe$stage_id,
    stage2n = stageframe$stage_id), fec32n = c(harpoon))
  
  stageexpansion3$indata32n <- stageexpansion3$indata3 * stageexpansion3$indata2n
  
  instages <- length(stageframe$stage_id) #Total number of stages, including the dead stage
  stageexpansion3$index21 <- apply(as.matrix(c(1:dim(stageexpansion3)[1])), 1, function(X) {
    ((stageexpansion3$stage3[X] - 1) + ((stageexpansion3$stage2n[X] - 1) * instages))
  })
  stageexpansion3$stcod3 <- apply(as.matrix(c(1:dim(stageexpansion3)[1])), 1, function(X) {
    stageframe$stage[which(stageframe$stage_id == stageexpansion3$stage3[X])]
  })
  stageexpansion3$stcod2 <- apply(as.matrix(c(1:dim(stageexpansion3)[1])), 1, function(X) {
    stageframe$stage[which(stageframe$stage_id == stageexpansion3$stage2n[X])]
  })
  
  # Now we will add a number of indices to the dataset
  data <- subset(data, alive2 == 1)
  
  data$index1 <- apply(as.matrix(data$usedstage1), 1, function(X) {
    stageframe$stage_id[which(stageframe$stage == X)]
  })
  data$index2 <- apply(as.matrix(data$usedstage2), 1, function(X) {
    stageframe$stage_id[which(stageframe$stage == X)]
  })
  data$index3 <- apply(as.matrix(data$usedstage3), 1, function(X) {
    stageframe$stage_id[which(stageframe$stage == X)]
  })
  
  data$index321 <- apply(as.matrix(c(1:length(data$usedstage1))), 1, function(X) {
    if (format_int == 2) {
      ((data$index3[X] - 1) + ((data$index2[X] - 1) * instages) +
        ((data$index2[X] - 1) * instages * instages) + 
        ((data$index1[X] - 1) * instages * instages * instages))
    } else {
      ((data$index3[X] - 1) + ((data$index2[X] - 1) * (instages - 1)) + 
        ((data$index2[X] - 1) * (instages - 1) * (instages - 1)) + 
        ((data$index1[X] - 1) * (instages - 1) * (instages - 1) * (instages - 1)))
    }
  })
  data$pairindex21 <- apply(as.matrix(c(1:dim(data)[1])), 1, function(X) {
    return(((data$index2[X] - 1) + ((data$index1[X] - 1) * instages)))
  })
  data$usedfec2[which(is.na(data$usedfec2))] <- 0
  
  if(is.element(0, unique(data$index1))) {
    warning("Data (stage at occasion t-1) contains stages not identified in stageframe.
      All stage characteristics must match, including reproductive status.",
      call. = FALSE)
  }
  if(is.element(0, unique(data$index2))) {
    warning("Data (stage at occasion t) contains stages not identified in stageframe.
      All stage characteristics must match, including reproductive status.",
      call. = FALSE)
  }
  if(is.element(0, unique(data$index3))) {
    warning("Data (stage at occasion t+1) contains stages not identified in stageframe.
      All stage characteristics must match, including reproductive status.",
      call. = FALSE)
  }
  
  madsexmadrigal <- lapply(yearlist, function(X) {
    passed_data <- data
    if (!is.na(X$pop[1]) & !is.na(popcol)) {
      passed_data$popused <- passed_data[,popcol];
      passed_data <- subset(passed_data, popused == X$pop[1]);
    }
    if (!is.na(X$patch[1]) & !is.na(patchcol)) {
      passed_data$patchused <- passed_data[,patchcol];
      passed_data <- subset(passed_data, patchused == X$patch[1]);
    }
    if (!is.na(X$year2[1])) {
      passed_data$yearused <- passed_data[,yearcol];
      passed_data <- subset(passed_data, yearused == X$year2[1]);
    }
    if (err_check) { err_push <- 1} else {err_push <- 0}
    
    .specialpatrolgroup(sge9l = stageexpansion9, sge3 = stageexpansion3,
      MainData = passed_data, StageFrame = stageframe, format = format_int, 
        err_switch = err_push)
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
  } else {
    warning("Individual identity variable not provided, affecting quality control output.",
      call. = FALSE)
  }
  qcoutput2 <- c(indivs, dim(data)[1])
  
  morebitstolose <- unique(c(which(stageexpansion3$stcod3 == "Dead"),
    which(stageexpansion3$stcod2 == "Dead"),
    which(stageexpansion3$stcod3 == "AlmostBorn")))
  stageexpansion3 <- stageexpansion3[-morebitstolose,]
  
  hstages <- stageexpansion3[,c("stage3", "stage2n", "stcod3", "stcod2")]
  names(hstages) = c("stage_id_2", "stage_id_1", "stage_2", "stage_1")
  
  if (reduce == TRUE) {
    drops <- .reducer3(a_list, u_list, f_list, hstages)
    
    a_list <- drops$A
    u_list <- drops$U
    f_list <- drops$F
    hstages <- drops$hstages
  }
  
  rownames(hstages) <- c(1:dim(hstages)[1])
  
  if (!err_check) {
    output <- list(A = a_list, U = u_list, F = f_list, hstages = hstages, 
      agestages = NA, ahstages = stageframe[1:(dim(stageframe)[1] - 1),],
      labels = listofyears, matrixqc = qcoutput1, dataqc = qcoutput2)
  } else {
    err1_concrp <- lapply(madsexmadrigal, function(X) {X$concrp})
    err2_s2f <- lapply(madsexmadrigal, function(X) {X$s2f})
    err3_dpr <- lapply(madsexmadrigal, function(X) {X$dataprior})
    
    output <- list(A = a_list, U = u_list, F = f_list, hstages = hstages, 
      agestages = NA, ahstages = stageframe[1:(dim(stageframe)[1] - 1),],
      labels = listofyears, matrixqc = qcoutput1, dataqc = qcoutput2,
      err1_concrp = err1_concrp, err2_s2f = err2_s2f, err3_dpr = err3_dpr)
  }
  
  class(output) <- "lefkoMat"
  
  return(output)
}

#' Create Raw Ahistorical Matrix Projection Model
#'
#' Function \code{rlefko2()} returns raw ahistorical MPMs, including the
#' associated component transition and fecundity matrices, a data frame
#' describing the ahistorical stages used, and a data frame describing the
#' population, patch, and occasion time associated with each matrix.
#' 
#' @param data A vertical demographic data frame, with variables corresponding 
#' to the naming conventions in \code{\link{verticalize3}()}.
#' @param stageframe A stageframe object that includes information on the size,
#' observation status, propagule status, immaturity status, and maturity status
#' of each ahistorical stage. Should also incorporate bin widths if size is
#' continuous.
#' @param year A variable corresponding to observation occasion, or a set
#' of such values, given in values associated with the \code{year} term used in
#' vital rate model development. Can also equal \code{"all"}, in which case
#' matrices will be estimated for all occasion times. Defaults to \code{"all"}.
#' @param pop A variable designating which populations will have matrices
#' estimated. Should be set to specific population names, or to \code{"all"} if
#' all populations should have matrices estimated.
#' @param patch A variable designating which patches or subpopulations will have
#' matrices estimated. Should be set to specific patch names, or to \code{"all"}
#' if matrices should be estimated for all patches. Defaults to \code{"all"}.
#' @param censor If \code{TRUE}, then data will be removed according to the
#' variable set in \code{censorcol}, such that only data with censor values
#' equal to \code{censorkeep} will remain. Defaults to \code{FALSE}.
#' @param stages An optional vector denoting the names of the variables within
#' the main vertical dataset coding for the stages of each individual in
#' occasions \emph{t}+1, \emph{t}, and \emph{t}-1. The names of stages in these
#' variables should match those used in the \code{stageframe} exactly. If left
#' blank, then \code{rlefko3()} will attempt to infer stages by matching values
#' of \code{alive}, \code{size}, \code{repst}, and \code{matst} to
#' characteristics noted in the associated \code{stageframe}.
#' @param alive A vector of names of binomial variables corresponding to status
#' as alive (\code{1}) or dead (\code{0}) in occasions \emph{t}+1 ans \emph{t},
#' respectively.
#' @param size A vector of names of variables coding the primary size variable
#' in occasions \emph{t}+1 and \emph{t}, respectively. Defaults to
#' \code{c("sizea3", "sizea2")}.
#' @param sizeb A vector of names of variables coding the secondary size
#' variable in occasions \emph{t}+1 and \emph{t}, respectively. Defaults to
#' \code{c(NA, NA)}.
#' @param sizec A vector of names of variables coding the tertiary size
#' variable in occasions \emph{t}+1 and \emph{t}, respectively. Defaults to
#' \code{c(NA, NA)}.
#' @param repst A vector of names of variables coding reproductive status in
#' occasions \emph{t}+1 and \emph{t}, respectively. Defaults to 
#' \code{c("repstatus3", "repstatus2")}. Must be supplied if \code{stages} is
#' not provided.
#' @param matst A vector of names of variables coding maturity status in
#' occasions \emph{t}+1 and \emph{t}, respectively. Defaults to
#' \code{c("matstatus3", "matstatus2")}. Must be supplied if \code{stages} is
#' not provided.
#' @param fec A vector of names of variables coding fecundity in occasions
#' \emph{t}+1 and \emph{t}, respectively. Defaults to \code{c("feca3", "feca2")}.
#' @param supplement An optional data frame of class \code{lefkoSD} that
#' provides supplemental data that should be incorporated into the MPM. Three
#' kinds of data may be integrated this way: transitions to be estimated via the
#' use of proxy transitions, transition overwrites from the literature or
#' supplemental studies, and transition multipliers for fecundity. This data
#' frame should be produced using the \code{\link{supplemental}()} function.
#' Should be used in place of or in addition to an overwrite table (see 
#' \code{overwrite} below) and a reproduction matrix (see \code{repmatrix}
#' below).
#' @param repmatrix An optional reproduction matrix. This matrix is composed
#' mostly of 0s, with non-zero entries acting as element identifiers and
#' multipliers for fecundity (with 1 equaling full fecundity). If left blank,
#' and no \code{supplement} is provided, then \code{rlefko2()} will assume that
#' all stages marked as reproductive produce offspring at 1x that of estimated
#' fecundity, and that offspring production will yield the first stage noted as
#' propagule or immature. To prevent this behavior, input just \code{0}, which
#' will result in fecundity being estimated only for transitions noted in
#' \code{supplement} above. Must be the dimensions of an ahistorical matrix.
#' @param overwrite An optional data frame developed with the
#' \code{\link{overwrite}()} function describing transitions to be overwritten
#' either with given values or with other estimated transitions. Note that this
#' function supplements overwrite data provided in \code{supplement}.
#' @param yearcol The variable name or column number corresponding to occasion 
#' \emph{t} in the dataset.
#' @param popcol The variable name or column number corresponding to the
#' identity of the population.
#' @param patchcol The variable name or column number corresponding to patch in
#' the dataset.
#' @param indivcol The variable name or column number coding individual
#' identity.
#' @param censorcol The variable name or column number denoting the censor
#' status. Only needed if \code{censor = TRUE}.
#' @param censorkeep The value of the censor variable denoting data elements to
#' keep. Defaults to \code{0}.
#' @param reduce A logical value denoting whether to remove historical stages
#' associated with only zero transitions. These are removed only if the
#' respective row and column sums in ALL matrices estimated equal 0. Defaults to
#' \code{FALSE}.
#' 
#' @return If all inputs are properly formatted, then this function will return
#' an object of class \code{lefkoMat}, which is a list that holds the matrix
#' projection model and all of its metadata. Its structure is a list with the
#' following elements:
#' 
#' \item{A}{A list of full projection matrices in order of sorted populations,
#' patches, and occasions. All matrices output in the \code{matrix} class.}
#' \item{U}{A list of survival transition matrices sorted as in \code{A}. All 
#' matrices output in the \code{matrix} class.}
#' \item{F}{A list of fecundity matrices sorted as in \code{A}. All matrices 
#' output in the \code{matrix} class.}
#' \item{hstages}{A data frame matrix showing the pairing of ahistorical stages
#' used to create historical stage pairs. Set to NA for ahistorical matrices.}
#' \item{agestages}{A data frame showing age-stage pairs. In this function, it
#' is set to NA. Only used in output to function \code{aflefko2}().}
#' \item{ahstages}{A data frame detailing the characteristics of associated
#' ahistorical stages, in the form of a modified stageframe that includes
#' status as an entry stage through reproduction.}
#' \item{labels}{A data frame giving the population, patch, and year of each 
#' matrix in order.}
#' \item{matrixqc}{A short vector describing the number of non-zero elements
#' in \code{U} and \code{F} matrices, and the number of annual matrices.}
#' \item{dataqc}{A vector showing the numbers of individuals and rows in the
#' vertical dataset used as input.}
#'
#' @section Notes:
#' The default behavior of this function is to estimate fecundity with regards
#' to transitions specified via associated fecundity multipliers in either
#' \code{supplement} or \code{repmatrix}. If both of these fields are left
#' empty, then fecundity will be estimated at full for all transitions leading
#' from reproductive stages to immature and propagule stages. However, if a
#' \code{supplement} is provided and a \code{repmatrix} is not, or if
#' \code{repmatrix} is set to 0, then only fecundity transitions noted in the
#' supplement will be set to non-zero values. To use the default behavior of
#' setting all reproductive stages to reproduce at full fecundity into immature
#' and propagule stages but also incorporate given or proxy survival
#' transitions, input those given and proxy transitions through the
#' \code{overwrite} options.
#' 
#' The reproduction matrix (field \code{repmatrix}) may only be supplied as
#' ahistorical. If provided as historical, then \code{rlefko2()} will fail and
#' produce an error.
#' 
#' Users may at times wish to estimate MPMs using a dataset incorporating
#' multiple patches or subpopulations. Should the aim of analysis be a general
#' MPM that does not distinguish these patches or subpopulations, the
#' \code{patchcol} variable should be left to \code{NA}, which is the default.
#' Otherwise the variable identifying patch needs to be named.
#'
#' Input options including multiple variable names must be entered in the order
#' of variables in occasion \emph{t}+1 and \emph{t}. Rearranging the order WILL
#' lead to erroneous calculations, and may lead to fatal errors.
#' 
#' Although this function is capable of assigning stages given an input
#' stageframe, it lacks the power of \code{\link{verticalize3}()} and
#' \code{\link{historicalize3}()} in this regard. Users are strongly
#' encouraged to use the latter two functions for stage assignment.
#' 
#' @seealso \code{\link{flefko3}()}
#' @seealso \code{\link{flefko2}()}
#' @seealso \code{\link{aflefko2}()}
#' @seealso \code{\link{fleslie}()}
#' @seealso \code{\link{rlefko3}()}
#' @seealso \code{\link{rleslie}()}
#' 
#' @examples
#' # Lathyrus example
#' data(lathyrus)
#' 
#' sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
#' stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
#' repvector <- c(0, 0, 0, 0, 0, 1, 0)
#' obsvector <- c(0, 1, 1, 1, 1, 1, 0)
#' matvector <- c(0, 0, 1, 1, 1, 1, 1)
#' immvector <- c(1, 1, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)
#' 
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, 
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector, 
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec, 
#'   propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, 
#'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9, 
#'   juvcol = "Seedling1988", sizeacol = "Volume88", repstracol = "FCODE88", 
#'   fecacol = "Intactseed88", deadacol = "Dead1988", nonobsacol = "Dormant1988", 
#'   stageassign = lathframe, stagesize = "sizea", censorcol = "Missing1988", 
#'   censorkeep = NA, censor = TRUE)
#' 
#' lathsupp2 <- supplemental(stage3 = c("Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "rep", "rep"),
#'   givenrate = c(0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 3, 3), stageframe = lathframe, historical = FALSE)
#' 
#' ehrlen2 <- rlefko2(data = lathvert, stageframe = lathframe, year = "all", 
#'   stages = c("stage3", "stage2"), supplement = lathsupp2, yearcol = "year2",
#'   indivcol = "individ")
#' 
#' summary(ehrlen2)
#' 
#' # Cypripedium example
#' rm(list=ls(all=TRUE))
#' data(cypdata)
#' 
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
#'   "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   propstatus = propvector, immstatus = immvector, indataset = indataset,
#'   binhalfwidth = binvec)
#' 
#' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
#'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
#'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
#'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
#'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
#'   NRasRep = TRUE)
#' 
#' # Here we use supplemental() to provide overwrite and reproductive info
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   stageframe = cypframe_raw, historical = FALSE)
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' 
#' cypmatrix2r$A[[1]]
#' 
#' @export
rlefko2 <- function(data, stageframe, year = "all", pop = NA, patch = NA,
  censor = FALSE, stages = NA, alive = c("alive3", "alive2"),
  size = c("sizea3", "sizea2"), sizeb = c(NA, NA), sizec = c(NA, NA),
  repst = c("repstatus3", "repstatus2"),
  matst = c("matstatus3", "matstatus2"), fec = c("feca3", "feca2"),
  supplement = NULL, repmatrix = NULL, overwrite = NULL, yearcol = NA,
  popcol = NA, patchcol = NA, indivcol = NA, censorcol = NA, censorkeep = 0,
  reduce = FALSE) {
  
  tocensor <- indataset <- alive2 <- popused <- patchused <- yearused <- NULL
  
  sizeb_used <- 0
  sizec_used <- 0
  
  if (all(is.na(data))) {
    stop("Need original vertical dataset to proceed.", call. = FALSE)
  }
  
  if (!any(class(data) == "data.frame")) {
    stop("Need original vertical dataset to proceed. This dataset must be in
      historical vertical format.", call. = FALSE)
  }
  
  if (!any(class(data) == "hfvdata")) {
    warning("Dataset used as input is not of class hfvdata. Will assume that the
      dataset has been formatted equivalently.", call. = FALSE)
  }
  
  stageframe_vars <- c("stage", "size", "size_b", "size_c", "min_age", "max_age",
    "repstatus", "obsstatus", "propstatus", "immstatus", "matstatus", "indataset",
    "binhalfwidth_raw", "sizebin_min", "sizebin_max", "sizebin_center",
    "sizebin_width", "binhalfwidthb_raw", "sizebinb_min", "sizebinb_max",
    "sizebinb_center", "sizebinb_width", "binhalfwidthc_raw", "sizebinc_min",
    "sizebinc_max", "sizebinc_center", "sizebinc_width", "group", "comments")
  if (any(!is.element(names(stageframe), stageframe_vars))) {
    stop("Please use properly formatted stageframe as input.", call. = FALSE)
  }
  
  if (all(is.na(stages))) {
    if (!(length(alive) > 1)) {
      stop("This function requires stage information for each of occasions t+1
        and t. In the absence of stage columns in the dataset, it requires two
        variables for living/dead status, size, reproductive status, and
        maturity status, for each of occasions t+1 and t.", call. = FALSE)
    }
    if (!(length(size) > 1)) {
      stop("This function requires stage information for each of occasions t+1
        and t. In the absence of stage columns in the dataset, it requires two
        variables for living/dead status, size, reproductive status, and
        maturity status, for each of occasions t+1 and t.", call. = FALSE)
    }
    if (!all(is.na(repst))) {
      if (!(length(repst) > 1)) {
        stop("This function requires stage information for each of occasions t+1
          and t. In the absence of stage columns in the dataset, it requires two
          variables for living/dead status, size, reproductive status, and
          maturity status, for each of occasions t+1 and t.", call. = FALSE)
      }
    }   
    if (!all(is.na(matst))) {
      if (!(length(matst) > 1)) {
        stop("This function requires stage information for each of occasions t+1
          and t. In the absence of stage columns in the dataset, it requires two
          variables for living/dead status, size, reproductive status, and
          maturity status, for each of occasions t+1 and t.", call. = FALSE)
      }
    }   
  }
  
  if (!(length(fec) > 1)) {
    stop("This function requires two variables for fecundity, for each of
      occasions t+1 and t.", call. = FALSE)
  }
  
  if (is.character(yearcol)) {
    choicevar <- which(names(data) == yearcol);
    mainyears <- sort(unique(data[,choicevar]))
  } else if (is.numeric(yearcol)) {
    mainyears <- sort(unique(data[, yearcol]))
  } else {
    stop("Need appropriate year column designation.", call. = FALSE)
  }
  
  if (any(is.character(year))) {
    if (is.element("all", tolower(year))) {
      year <- mainyears
    } else {
      stop("Year designation not recognized.", call. = FALSE)
    }
  }
  
  if (length(year) == 0 | all(is.na(year) == TRUE) | any(is.na(year))) {
    stop("This function cannot proceed without being given a specific year, or a
      suite of years. NA entries are not allowed.", call. = FALSE)
  }
  
  if (!all(is.element(year, mainyears))) {
    stop("Dataset does not contain one or more of the requested years. Note that
      matrices cannot be made for the first year in a historical dataset.",
      call. = FALSE)
  }
  
  if (censor == TRUE) {
    if(all(is.na(censorcol)) == TRUE) {
      stop("Cannot censor the data without a proper censor variable.",
        call. = FALSE)
    }
    
    if (all(is.character(censorcol))) {
      if (!all(is.element(censorcol, names(data)))) {
        stop("Censor variable names input for censorcol do not match any
          variable names in the dataset.", call. = FALSE)
      }
    }
    
    censorcolsonly <- data[,censorcol]
    sleeplessnights <- apply(as.matrix(c(1:dim(censorcolsonly)[1])), 1, function(X) {
      crazyvec <- if(is.element(censorkeep, censorcolsonly[X,])) {
        return(X);
      } else {
        return(NA);
      }
    })
    sleeplessnights <- sleeplessnights[!is.na(sleeplessnights)]
    
    data <- data[sleeplessnights,]
  }
  
  if (!all(is.na(pop)) & !all(is.na(patch))) {
    if (is.na(popcol) | is.na(patchcol)) {
      stop("Need population and patch designation variables to proceed.",
        call. = FALSE)
    }
    
    if (is.element("all", tolower(pop))) {
      if (is.character(popcol)) {popcol <- which(names(data) == popcol)}
      
      pops <- unique(data[,popcol])
    } else {pops <- pop}
    
    if (is.element("all", tolower(patch))) {
      if (is.character(patchcol)) {patchcol <- which(names(data) == patchcol)}
      
      listofpatches <- apply(as.matrix(pops), 1, function(X) {
        patchfolly <- subset(data, popcol == X);
        output <- cbind.data.frame(X, unique(patchfolly[,yearcol]), 
          stringsAsFactors = FALSE);
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
      output <- cbind.data.frame(X[1], X[2], unique(checkyrdata[,yearcol]),
        stringsAsFactors = FALSE);
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
  } else if (all(is.na(pop)) & !all(is.na(patch))) {
    if (!is.na(popcol)) {
      popcol <- NA
    }
    
    if (is.na(patchcol)) {
      stop("Need patch designation variable to proceed.", call. = FALSE)
    }
    
    if (is.element("all", tolower(patch))) {
      if (is.character(patchcol)) {patchcol <- which(names(data) == patchcol)}
      
      patches <- unique(data[,patchcol])
    } else {patches <- patch}
    
    listofyears <- apply(as.matrix(patches), 1, function(X) {
      checkyrdata <- subset(data, patchcol = X);
      output <- cbind.data.frame("1", X, unique(checkyrdata[,yearcol]),
        stringsAsFactors = FALSE);
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    if (length(listofyears) > 1) {
      listofyears <- do.call(rbind.data.frame, listofyears)
    } else {
      listofyears <- listofyears[[1]]
    }
  } else if (!all(is.na(pop)) & all(is.na(patch))) {
    if (is.na(popcol)) {
      stop("Need population designation variable to proceed.", call. = FALSE)
    }
    
    if (!is.na(patchcol)) {
      patchcol <- NA
    }
    
    if (is.element("all", tolower(pop))) {
      if (is.character(popcol)) {popcol <- which(names(data) == popcol)}
      
      pops <- unique(data[,popcol])
    } else {pops <- pop}
    
    listofyears <- apply(as.matrix(pops), 1, function(X) {
      checkyrdata <- subset(data, popcol = X);
      output <- cbind.data.frame(X, "1", unique(checkyrdata[,yearcol]),
        stringsAsFactors = FALSE);
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    if (length(listofyears) > 1) {
      listofyears <- do.call(rbind.data.frame, listofyears)
    } else {
      listofyears <- listofyears[[1]]
    }
  } else if (all(is.na(pop)) & all(is.na(patch))) {
    if (!is.na(popcol)) {
      popcol <- NA
    }
    if (!is.na(patchcol)) {
      patchcol <- NA
    }
    
    listofyears <- cbind.data.frame("1", "1", year, stringsAsFactors = FALSE)
    names(listofyears) <- c("pop", "patch", "year2")
  }
  
  identifiedyearrows <- which(is.element(listofyears$year2, year))
  if (length(identifiedyearrows) == 0) {
    stop("Cannot recognize input year(s)", call. = FALSE)
  } else {
    listofyears <- listofyears[identifiedyearrows,]
  }
  yearlist <- split(listofyears, seq(nrow(listofyears)))
  
  stagenum_init <- dim(stageframe)[1]
  if (!all(is.null(repmatrix))) {
    if (any(class(repmatrix) == "matrix")) {
      if (dim(repmatrix)[1] != stagenum_init) {
        stop("The repmatrix must be a square matrix with dimensions equal to the
          number of stages in the stageframe.", call. = FALSE)
      }
      
      if (dim(repmatrix)[2] != stagenum_init) {
        stop("The repmatrix must be a square matrix with dimensions equal to the
          number of stages in the stageframe.", call. = FALSE)
      }
    }
  }
  
  melchett <- .sf_reassess(stageframe, supplement, overwrite, repmatrix,
    agemat = FALSE, historical = FALSE, format = 1)
  stageframe <- melchett$stageframe
  repmatrix <- melchett$repmatrix
  ovtable <- melchett$ovtable
  
  if (!all(is.na(overwrite)) | !all(is.na(supplement))) {
    
    if(any(duplicated(ovtable[,1:3]))) {
      stop("Multiple entries with different values for the same stage transition
        are not allowed in the supplemental or overwrite table. If modifying a
        historical table to perform an ahistorical analysis, then this may be
        due to different given rates of substitutions caused by dropping stage
        at occasion t-1. Please eliminate duplicate transitions.",
        call. = FALSE)
    }
  }
  
  data$alive2 <- data[,which(names(data) == alive[2])]
  data$alive3 <- data[,which(names(data) == alive[1])]
  
  instageframe <- subset(stageframe, indataset == 1)
  instages <- dim(stageframe)[1] #This is actually the total number of stages, including the dead stage
  
  if (all(is.na(stages))) {
    if (length(size) > 1) {
      size2met <- which(names(data) == size[2])
      size3met <- which(names(data) == size[1])
      
      if (length(size2met) == 1 & length(size3met) == 1) {
        data$usedsize2 <- data[,size2met]
        data$usedsize3 <- data[,size3met]
      } else {
        stop("Entered size variable names do not strictly correspond to single
          variables in the dataset.", call. = FALSE)
      }
      
      if (!all(is.na(sizeb)) & length(sizeb) > 1) {
        size2met <- which(names(data) == sizeb[2])
        size3met <- which(names(data) == sizeb[1])
        
        if (length(size2met) == 1 & length(size3met) == 1) {
          sizeb_used <- 1
          data$usedsize2b <- data[,size2met]
          data$usedsize3b <- data[,size3met]
        } else {
          stop("Entered sizeb variable names do not strictly correspond to single
          variables in the dataset, or fewer than 2 variable names have been
          entered.", call. = FALSE)
        }
      }
      
      if (!all(is.na(sizec)) & length(sizec) > 1) {
        size2met <- which(names(data) == sizec[2])
        size3met <- which(names(data) == sizec[1])
        
        if (length(size2met) == 1 & length(size3met) == 1) {
          sizec_used <- 1
          data$usedsize2c <- data[,size2met]
          data$usedsize3c <- data[,size3met]
        } else {
          stop("Entered sizec variable names do not strictly correspond to single
          variables in the dataset, or fewer than 2 variable names have been
          entered.", call. = FALSE)
        }
      }
    } else {
      warning("Without stage columns, rlefko3() requires size variables as input.
        Failure to include size variables may lead to odd results.", call. = FALSE)
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
      mainstages <- intersect(which(instageframe$sizebin_min < data$usedsize2[X]), 
        which(instageframe$sizebin_max >= data$usedsize2[X]))
      if (sizeb_used == 1) {
        if (is.na(data$usedsize2b[X])) {
          data$usedsize2b[X] <- 0
        }
        mainstages <- intersect(which(instageframe$sizebinb_min < data$usedsize2b[X]), mainstages)
        mainstages <- intersect(which(instageframe$sizebinb_max >= data$usedsize2b[X]), mainstages)
      }
      if (sizec_used == 1) {
        if (is.na(data$usedsize2c[X])) {
          data$usedsize2c[X] <- 0
        }
        mainstages <- intersect(which(instageframe$sizebinc_min < data$usedsize2c[X]), mainstages)
        mainstages <- intersect(which(instageframe$sizebinc_max >= data$usedsize2c[X]), mainstages)
      }
      jmstages <- which(instageframe$immstatus == (1 - data$usedmatstatus2[X]))
      obsstages <- which(instageframe$obsstatus == data$obsstatus2[X])
      repstages <- which(instageframe$repstatus == data$repstatus2[X])
        
      choicestage <- intersect(intersect(mainstages, jmstages),
        intersect(obsstages, repstages))
      
      if (length(choicestage) == 0) {
        choicestage <- which(stageframe$stage_id == max(stageframe$stage_id))
        if (data$alive2[X] != 0) {
          stop("Stage characteristics mismatch dataset. Consider using the stages
            option, particularly if the vertical file was created with NRasRep = TRUE
            in verticalize3() or historicalize3().", call. = FALSE)
        }
      }
      
      return(as.character(instageframe$stage[choicestage]))
    })
    
    data$usedstage3 <- apply(as.matrix(c(1:dim(data)[1])), 1, function(X) {
      if (is.na(data$usedsize3[X])) {
        data$usedsize3[X] <- 0
      }
      mainstages <- intersect(which(instageframe$sizebin_min < data$usedsize3[X]), 
        which(instageframe$sizebin_max >= data$usedsize3[X]))
      if (sizeb_used == 1) {
        if (is.na(data$usedsize3b[X])) {
          data$usedsize3b[X] <- 0
        }
        mainstages <- intersect(which(instageframe$sizebinb_min < data$usedsize3b[X]), mainstages)
        mainstages <- intersect(which(instageframe$sizebinb_max >= data$usedsize3b[X]), mainstages)
      }
      if (sizec_used == 1) {
        if (is.na(data$usedsize3c[X])) {
          data$usedsize3c[X] <- 0
        }
        mainstages <- intersect(which(instageframe$sizebinc_min < data$usedsize3c[X]), mainstages)
        mainstages <- intersect(which(instageframe$sizebinc_max >= data$usedsize3c[X]), mainstages)
      }
      jmstages <- which(instageframe$immstatus == (1 - data$usedmatstatus3[X]))
      obsstages <- which(instageframe$obsstatus == data$obsstatus3[X])
      repstages <- which(instageframe$repstatus == data$repstatus3[X])
      alivestage3 <- which(instageframe$alive == data$alive3[X])
      
      choicestage <- intersect(intersect(intersect(mainstages, jmstages),
        intersect(obsstages, repstages)), alivestage3)
      
      if (length(choicestage) == 0) {
        choicestage <- which(instageframe$stage_id == max(instageframe$stage_id))
        if (data$alive3[X] != 0) {
          stop("Stage characteristics mismatch dataset. Consider using the stages
            option, particularly if the vertical file was created with NRasRep = TRUE
            in verticalize3() or historicalize3().", call. = FALSE)
        }
      }
      
      return(as.character(instageframe$stage[choicestage]))
    })
  } else if (length(stages) > 1) {
    if (is.numeric(stages[2])) {
      data$usedstage2 <- data[, stages[2]]
      data$usedstage3 <- data[, stages[1]]
      
      data$usedstage2 <- as.character(data$usedstage2)
      data$usedstage3 <- as.character(data$usedstage3)
      
      if (is.element("NotAlive", unique(data$usedstage3))) {
        data$usedstage3[which(data$usedstage3 == "NotAlive")] <- "Dead"
      }
    } else {
      data$usedstage2 <- data[,which(names(data) == stages[2])]
      data$usedstage3 <- data[,which(names(data) == stages[1])]
      
      data$usedstage2 <- as.character(data$usedstage2)
      data$usedstage3 <- as.character(data$usedstage3)
      
      if (is.element("NotAlive", unique(data$usedstage3))) {
        data$usedstage3[which(data$usedstage3 == "NotAlive")] <- "Dead"
      }
    }
    stages.used <- sort(unique(c(data$usedstage2, data$usedstage3)))
    
    if (length(setdiff(stages.used, stageframe$stage)) > 0) {
      stop("Some stages in dataset do not match those detailed in the input stageframe.",
        call. = FALSE)
    }
  }
  
  if (length(fec) > 1) {
    data$usedfec2 <- data[,which(names(data) == fec[2])]
    data$usedfec3 <- data[,which(names(data) == fec[1])]
    
    data$usedfec2[which(is.na(data$usedfec2))] <- 0
    data$usedfec3[which(is.na(data$usedfec3))] <- 0
  } else {
    warning("Function rlefko2() requires 2 fecundity variables, for times t+1 and t. 
      Failure to include fecundity variables leads to matrices composed only of 
      survival transitions.", call. = FALSE)
  } 
  
  # This section creates stageexpansion3, a data frame with stage transition values from occasion t to t+1
  stageexpansion3 <- .theoldpizzle(stageframe, ovtable, repmatrix, finalage = 0, 
    format = 1, style = 1, cont = 0)
  
  # Stageexpansion2 is a dataframe created to hold values for stages in occasion t only
  stageexpansion2 <- cbind.data.frame(stage2 = as.numeric(stageframe$stage_id),
    size2 = as.numeric(stageframe$sizebin_center), sizeb2 = as.numeric(stageframe$sizebinb_center),
    sizec2 = as.numeric(stageframe$sizebinc_center), rep2 = as.numeric(stageframe$repstatus),
    indata2 = as.numeric(stageframe$indataset), index2 = (as.numeric(stageframe$stage_id) - 1),
    fec3 = c(rowSums(repmatrix), 0))
  stageexpansion2$fec3[which(stageexpansion2$fec3 > 0)] <- 1
  
  data <- subset(data, alive2 == 1)
  
  data$index2 <- apply(as.matrix(data$usedstage2), 1, function(X) {
    instageframe$stage_id[which(instageframe$stage == X)] - 1
  })
  data$index2[which(is.na(data$index2))] <- 0
  data$index3 <- apply(as.matrix(data$usedstage3), 1, function(X) {
    instageframe$stage_id[which(instageframe$stage == X)] - 1
  })
  data$index3[which(is.na(data$index3))] <- 0
  data$index32 <- apply(as.matrix(c(1:length(data$usedstage2))), 1, function(X) {
    (data$index3[X] + (data$index2[X] * instages))
  })
  
  if(is.element(0, unique(data$index2))) {
    warning("Data (stage at occasion t) contains stages not identified in stageframe.
      All stage characteristics must match, including reproductive status.", 
      call. = FALSE)
  }
  if(is.element(0, unique(data$index3))) {
    warning("Data (stage at occasion t+1) contains stages not identified in stageframe.
      All stage characteristics must match, including reproductive status.",
      call. = FALSE)
  }
  
  # This section runs the core matrix estimator
  madsexmadrigal <- lapply(yearlist, function(X) {
    passed_data <- data
    if (!is.na(X$pop[1]) & !is.na(popcol)) {
      passed_data$popused <- passed_data[,popcol];
      passed_data <- subset(passed_data, popused == X$pop[1]);
    }
    if (!is.na(X$patch[1]) & !is.na(patchcol)) {
      passed_data$patchused <- passed_data[,patchcol];
      passed_data <- subset(passed_data, patchused == X$patch[1]);
    }
    if (!is.na(X$year2[1])) {
      passed_data$yearused <- passed_data[,yearcol];
      passed_data <- subset(passed_data, yearused == X$year2[1]);
    }
    .normalpatrolgroup(sge3 = stageexpansion3, sge2 = stageexpansion2, 
      MainData = passed_data, StageFrame = stageframe)
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
  } else {
    warning("Individual identity variable not provided, affecting quality control output.",
      call. = FALSE)
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
  
  output <- list(A = a_list, U = u_list, F = f_list, hstages = NA,
    agestages = NA, ahstages = ahstages, labels = listofyears,
    matrixqc = qcoutput1, dataqc = qcoutput2)
  
  class(output) <- "lefkoMat"
  
  return(output)
}

#' Create Function-based Ahistorical Age x Stage Matrix Projection Model
#'
#' Function \code{aflefko2()} returns ahistorical age x stage MPMs corresponding
#' to the patches and occasion times given, including the associated component
#' transition and fecundity matrices, data frames detailing the characteristics
#' of ahistorical stages and the exact age-stage combinations corresponding to
#' rows and columns in estimated matrices, and a data frame characterizing the
#' patch and occasion time combinations corresponding to these matrices.
#'
#' @param year A variable corresponding to observation occasion, or a set
#' of such values, given in values associated with the year term used in linear 
#' model development. Defaults to \code{"all"}, in which case matrices will be
#' estimated for all occasions.
#' @param patch A variable designating which patches or subpopulations will have
#' matrices estimated. Defaults to \code{"all"}, but can also be set to specific
#' patch names.
#' @param stageframe A stageframe object that includes information on the size,
#' observation status, propagule status, immaturity status, and maturity status
#' of each ahistorical stage. Should also incorporate bin widths if size is
#' continuous.
#' @param supplement An optional data frame of class \code{lefkoSD} that
#' provides supplemental data that should be incorporated into the MPM. Three
#' kinds of data may be integrated this way: transitions to be estimated via the
#' use of proxy transitions, transition overwrites from the literature or
#' supplemental studies, and transition multipliers for survival and fecundity.
#' This data frame should be produced using the \code{\link{supplemental}()}
#' function. Can be used in place of or in addition to an overwrite table (see 
#' \code{overwrite} below) and a reproduction matrix (see \code{repmatrix}
#' below).
#' @param repmatrix An optional reproduction matrix. This matrix is composed
#' mostly of 0s, with non-zero entries acting as element identifiers and
#' multipliers for fecundity (with 1 equaling full fecundity). If left blank,
#' and no \code{supplement} is provided, then \code{aflefko2()} will assume that
#' all stages marked as reproductive produce offspring at 1x that of estimated
#' fecundity, and that offspring production will yield the first stage noted as
#' propagule or immature.  To prevent this behavior, input just \code{0}, which
#' will result in fecundity being estimated only for transitions noted in
#' \code{supplement} above. Must be the dimensions of an ahistorical stage-based
#' matrix.
#' @param overwrite An optional data frame developed with the
#' \code{\link{overwrite}()} function describing transitions to be overwritten
#' either with given values or with other estimated transitions. Note that this
#' function supplements overwrite data provided in \code{supplement}.
#' @param data The historical vertical demographic data frame used to estimate
#' vital rates (class \code{hfvdata}). The original data frame is required in
#' order to initialize occasions and patches properly.
#' @param modelsuite An optional \code{lefkoMod} object holding the vital rate
#' models. If given, then \code{surv_model}, \code{obs_model}, 
#' \code{size_model}, \code{sizeb_model}, \code{sizec_model},
#' \code{repst_model}, \code{fec_model}, \code{jsurv_model}, \code{jobs_model},
#' \code{jsize_model}, \code{jsizeb_model}, \code{jsizec_model},
#' \code{jrepst_model}, \code{jmatst_model}, \code{paramnames}, \code{yearcol},
#' and \code{patchcol} are not required. No models should include size or
#' reproductive status in occasion \emph{t}-1.
#' @param surv_model A linear model predicting survival probability. This can be
#' a model of class \code{glm} or \code{glmer}, and requires a predicted
#' binomial variable under a logit link. If given, then will overwrite any
#' survival probability model given in \code{modelsuite}. This model must have
#' been developed in a modeling exercise testing only the impacts of occasion
#' \emph{t}.
#' @param obs_model A linear model predicting sprouting or observation
#' probability. This can be a model of class \code{glm} or \code{glmer}, and
#' requires a predicted binomial variable under a logit link. If given, then
#' will overwrite any observation probability model given in \code{modelsuite}.
#' This model must have been developed in a modeling exercise testing only the
#' impacts of occasion \emph{t}.
#' @param size_model A linear model predicting primary size. This can be a model
#' of class \code{glm}, \code{glmer}, \code{glmmTMB}, \code{zeroinfl},
#' \code{vglm}, \code{lm}, or \code{lmer}. If given, then will overwrite any
#' primary size model given in \code{modelsuite}. This model must have been
#' developed in a modeling exercise testing only the impacts of occasion
#' \emph{t}.
#' @param sizeb_model A linear model predicting secondary size. This can be a
#' modelvof class \code{glm}, \code{glmer}, \code{glmmTMB}, \code{zeroinfl},
#' \code{vglm}, \code{lm}, or \code{lmer}. If given, then will overwrite any
#' secondary size model given in \code{modelsuite}. This model must have been
#' developed in a modeling exercise testing only the impacts of occasion
#' \emph{t}.
#' @param sizec_model A linear model predicting tertiary size. This can be a
#' modelvof class \code{glm}, \code{glmer}, \code{glmmTMB}, \code{zeroinfl},
#' \code{vglm}, \code{lm}, or \code{lmer}. If given, then will overwrite any
#' tertiary size model given in \code{modelsuite}. This model must have been
#' developed in a modeling exercise testing only the impacts of occasion
#' \emph{t}.
#' @param repst_model A linear model predicting reproduction probability. This
#' can be a model of class \code{glm} or \code{glmer}, and requires a predicted
#' binomial variable under a logit link. If given, then will overwrite any
#' reproduction probability model given in \code{modelsuite}. This model must
#' have been developed in a modeling exercise testing only the impacts of
#' occasion \emph{t}.
#' @param fec_model A linear model predicting fecundity. This can be a model of
#' class \code{glm}, \code{glmer}, \code{glmmTMB}, \code{zeroinfl},
#' \code{vglm}, \code{lm}, or \code{lmer}. If given, then will overwrite any 
#' fecundity model given in \code{modelsuite}. This model must have been
#' developed in a modeling exercise testing only the impacts of occasion
#' \emph{t}.
#' @param jsurv_model A linear model predicting juvenile survival probability.
#' This can be a model of class \code{glm} or \code{glmer}, and requires a
#' predicted binomial variable under a logit link. If given, then will overwrite
#' any juvenile survival probability model given in \code{modelsuite}. This
#' model must have been developed in a modeling exercise testing only the
#' impacts of occasion \emph{t}.
#' @param jobs_model A linear model predicting juvenile sprouting or observation
#' probability. This can be a model of class \code{glm} or \code{glmer}, and
#' requires a predicted binomial variable under a logit link. If given, then
#' will overwrite any juvenile observation probability model given in
#' \code{modelsuite}. This model must have been developed in a modeling exercise
#' testing only the impacts of occasion \emph{t}.
#' @param jsize_model A linear model predicting juvenile primary size. This
#' can be a model of class \code{glm}, \code{glmer}, \code{glmmTMB},
#' \code{zeroinfl}, \code{vglm}, \code{lm}, or \code{lmer}. If given, then will
#' overwrite any juvenile primary size model given in \code{modelsuite}. This
#' model must have been developed in a modeling exercise testing only the
#' impacts of occasion \emph{t}.
#' @param jsizeb_model A linear model predicting juvenile secondary size. This
#' can be a model of class \code{glm}, \code{glmer}, \code{glmmTMB},
#' \code{zeroinfl}, \code{vglm}, \code{lm}, or \code{lmer}. If given, then will
#' overwrite any juvenile secondary size model given in \code{modelsuite}. This
#' model must have been developed in a modeling exercise testing only the
#' impacts of occasion \emph{t}.
#' @param jsizec_model A linear model predicting juvenile tertiary size. This
#' can be a model of class \code{glm}, \code{glmer}, \code{glmmTMB},
#' \code{zeroinfl}, \code{vglm}, \code{lm}, or \code{lmer}. If given, then will
#' overwrite any juvenile tertiary size model given in \code{modelsuite}. This
#' model must have been developed in a modeling exercise testing only the
#' impacts of occasion \emph{t}.
#' @param jrepst_model A linear model predicting reproduction probability of an 
#' individual that was immature in time \emph{t}. This can be a model of class
#' \code{glm} or \code{glmer}, and requires a predicted binomial variable under
#' a logit link. If given, then will overwrite any reproduction probability
#' model given in \code{modelsuite}. This model must have been developed in a
#' modeling exercise testing only the impacts of occasion \emph{t}.
#' @param jmatst_model A linear model predicting maturity probability of an 
#' individual that was immature in time \emph{t}. This can be a model of class
#' \code{glm} or \code{glmer}, and requires a predicted binomial variable under
#' a logit link. If given, then will overwrite any reproduction probability
#' model given in \code{modelsuite}. This model must have been developed in a
#' modeling exercise testing only the impacts of occasion \emph{t}.
#' @param paramnames A dataframe with three columns, the first describing all
#' terms used in linear modeling, the second (must be called \code{mainparams}),
#' showing the general model terms that will be used in matrix creation (users
#' should use \code{\link{modelsearch}()} at least once to see the proper
#' names to be used in this column), and the third showing the equivalent terms
#' used in modeling (must be named \code{modelparams}). Only required if
#' \code{modelsuite} is not supplied.
#' @param inda Can be a single value to use for individual covariate \code{a}
#' in all matrices, a pair of values to use for times \emph{t} and \emph{t}-1 in
#' historical matrices, or a vector of such values corresponding to each
#' occasion in option \code{year}. Defaults to \code{NULL}.
#' @param indb Can be a single value to use for individual covariate \code{b}
#' in all matrices, a pair of values to use for times \emph{t} and \emph{t}-1 in
#' historical matrices, or a vector of such values corresponding to each
#' occasion in option \code{year}. Defaults to \code{NULL}.
#' @param indc Can be a single value to use for individual covariate \code{c}
#' in all matrices, a pair of values to use for times \emph{t} and \emph{t}-1 in
#' historical matrices, or a vector of such values corresponding to each
#' occasion in option \code{year}. Defaults to \code{NULL}.
#' @param surv_dev A numeric value to be added to the y-intercept in the linear
#' model for survival probability.
#' @param obs_dev A numeric value to be added to the y-intercept in the linear
#' model for observation probability.
#' @param size_dev A numeric value to be added to the y-intercept in the linear
#' model for primary size.
#' @param sizeb_dev A numeric value to be added to the y-intercept in the linear
#' model for secondary size.
#' @param sizec_dev A numeric value to be added to the y-intercept in the linear
#' model for tertiary size.
#' @param repst_dev A numeric value to be added to the y-intercept in the linear
#' model for probability of reproduction.
#' @param fec_dev A numeric value to be added to the y-intercept in the linear
#' model for fecundity.
#' @param jsurv_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile survival probability.
#' @param jobs_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile observation probability.
#' @param jsize_dev A numeric value to be added to the y-intercept in the linear
#' model for juvenile primary size.
#' @param jsizeb_dev A numeric value to be added to the y-intercept in the
#' linear model for juvenile secondary size.
#' @param jsizec_dev A numeric value to be added to the y-intercept in the
#' linear model for juvenile tertiary size.
#' @param jrepst_dev A numeric value to be added to the y-intercept in the
#' linear model for juvenile reproduction probability.
#' @param jmatst_dev A numeric value to be added to the y-intercept in the
#' linear model for juvenile maturity probability.
#' @param density A numeric value indicating density value to use to propagate
#' matrices. Only needed if density is an explanatory term used in linear
#' models. Defaults to \code{NA}.
#' @param repmod A scalar multiplier of fecundity. Defaults to \code{1}.
#' @param yearcol The variable name or column number corresponding to year in
#' occasion \emph{t} in the dataset. Not needed if a \code{modelsuite} is
#' supplied.
#' @param patchcol The variable name or column number corresponding to patch in 
#' the dataset. Not needed if a \code{modelsuite} is supplied.
#' @param year.as.random A logical term indicating whether coefficients for
#' missing occasions within vital rate models should be estimated as random
#' intercepts. Defaults to \code{FALSE}, in which case missing monitoring
#' occasion coefficients are set to \code{0}.
#' @param patch.as.random A logical term indicating whether coefficients for
#' missing patches within vital rate models should be estimated as random
#' intercepts. Defaults to \code{FALSE}, in which case missing patch
#' coefficients are set to \code{0}.
#' @param random.inda A logical value denoting whether to treat individual
#' covariate \code{a} as a random, categorical variable. Otherwise is treated as
#' a fixed, numeric variable. Defaults to \code{FALSE}.
#' @param random.indb A logical value denoting whether to treat individual
#' covariate \code{b} as a random, categorical variable. Otherwise is treated as
#' a fixed, numeric variable. Defaults to \code{FALSE}.
#' @param random.indc A logical value denoting whether to treat individual
#' covariate \code{c} as a random, categorical variable. Otherwise is treated as
#' a fixed, numeric variable. Defaults to \code{FALSE}.
#' @param final_age The final age to model in the matrix, where the first age
#' will be age 0.
#' @param continue A logical value designating whether to allow continued
#' survival of individuals past the final age noted in the stageframe, using the 
#' demographic characteristics of the final age. Defaults to \code{TRUE}.
#' @param randomseed A numeric value used as a seed to generate random estimates
#' for missing occasion and patch coefficients, if either \code{year.as.random}
#' or \code{patch.as.random} is set to \code{TRUE}. Defaults to
#' \code{\link{set.seed}()} default.
#' @param negfec A logical value denoting whether fecundity values estimated to
#' be negative should be reset to \code{0}. Defaults to \code{FALSE}.
#' @param ipm_method A string indicating what method to use to estimate size
#' transition probabilities if size is treated as continuous. Options include
#' \code{"midpoint"}, which utilizes the midpoint method, and \code{"CDF"},
#' which uses the cumulative distribution function. Defaults to \code{"CDF"}.
#' @param reduce A logical value denoting whether to remove ahistorical stages
#' associated solely with 0 transitions. These are only removed in cases where
#' the associated row and column sums in ALL matrices estimated equal 0.
#' Defaults to \code{FALSE}.
#' @param err_check A logical value indicating whether to append matrices of
#' vital rate probabilities associated with each matrix to the \code{lefkoMat}
#' object generated. These matrices are developed internally and can be used for
#' error checking. Defaults to \code{FALSE}.
#' @param exp_tol A numeric value used to indicate a maximum value to set
#' exponents to in the core kernel to prevent numerical overflow. Defaults to
#' \code{700}.
#' @param theta_tol A numeric value used to indicate a maximum value to theta as
#' used in the negative binomial probability density kernel. Defaults to
#' \code{100000000}, but can be reset to other values during error checking.
#'
#' @return If all inputs are properly formatted, then this function will return
#' an object of class \code{lefkoMat}, which is a list that holds the matrix
#' projection model and all of its metadata. Its structure is a list with the
#' following elements:
#'
#' \item{A}{A list of full projection matrices in order of sorted patches and
#' occasions. All matrices output in R's \code{matrix} class.}
#' \item{U}{A list of survival transition matrices sorted as in \code{A}. All 
#' matrices output in R's \code{matrix} class.}
#' \item{F}{A list of fecundity matrices sorted as in \code{A}. All matrices 
#' output in R's \code{matrix} class.}
#' \item{hstages}{A data frame matrix showing the pairing of ahistorical stages
#' used to create historical stage pairs. Set to \code{NA} for age-by-stage
#' MPMs.}
#' \item{agestages}{A data frame showing the stage number and stage name
#' corresponding to \code{ahstages}, as well as the associated age, of each
#' row in each age-by-stage matrix.}
#' \item{ahstages}{A data frame detailing the characteristics of associated
#' ahistorical stages, in the form of a modified stageframe that includes
#' status as an entry stage through reproduction.}
#' \item{labels}{A data frame giving the patch and year of each matrix in order.
#' In \code{aflefko2()}, only one population may be analyzed at once, and so
#' \code{pop = NA}}
#' \item{matrixqc}{A short vector describing the number of non-zero elements
#' in \code{U} and \code{F} matrices, and the number of annual matrices.}
#' \item{modelqc}{This is the \code{qc} portion of the modelsuite input.}
#' \item{prob_out}{An optional element only added if \code{err_check = TRUE}.
#' This is a list of vital rate probability matrices, with 6 columns in the
#' order of survival, observation probability, reproduction probability, primary
#' size transition probability, secondary size transition probability, and
#' tertiary size transition probability.}
#' 
#' @section Notes:
#' Unlike \code{\link{rlefko2}()} \code{\link{rlefko3}()}, and
#' \code{\link{rleslie}()}, this function does not currently distinguish
#' populations.
#' 
#' This function will yield incorrect estimates if the models utilized
#' incorporate state in occasion \emph{t}-1. Only use models developed testing
#' for ahistorical effects.
#' 
#' The default behavior of this function is to estimate fecundity with regards
#' to transitions specified via associated fecundity multipliers in either
#' \code{supplement} or \code{repmatrix}. If both of these fields are left
#' empty, then fecundity will be estimated at full for all transitions leading
#' from reproductive stages to immature and propagule stages. However, if a
#' \code{supplement} is provided and a \code{repmatrix} is not, or if
#' \code{repmatrix} is set to 0, then only fecundity transitions noted in the
#' supplement will be set to non-zero values. To use the default behavior of
#' setting all reproductive stages to reproduce at full fecundity into immature
#' and propagule stages but also incorporate given or proxy
#' survival transitions, input those given and proxy transitions through the
#' \code{overwrite} option.
#' 
#' The reproduction matrix (field \code{repmatrix}) may only be supplied as
#' ahistorical. If provided as historical, then \code{aflefko2()} will fail and
#' produce an error.
#' 
#' Stageframes used in this function should include ages for minimum and maximum
#' age for each stage. \code{NA}s are treated as \code{0}s in minimum age, and
#' as \code{final_age} for maximum age.
#' 
#' Users may at times wish to estimate MPMs using a dataset incorporating
#' multiple patches or subpopulations, but without discriminating between those
#' patches or subpopulations. Should the aim of analysis be a general MPM that
#' does not distinguish these patches or subpopulations, the \code{patchcol}
#' variable should be set to \code{NA}, which is the default.
#'
#' Input options including multiple variable names must be entered in the order
#' of variables in occasion \emph{t}+1 and \emph{t}. Rearranging the order will
#' lead to erroneous calculations, and may lead to fatal errors.
#'
#' Care should be taken to match the random status of year and patch to the
#' states of those variables within the modelsuite. If they do not match, then
#' they will be treated as zeroes in vital rate estimation.
#' 
#' The \code{ipm_method} function gives the option of using two different means
#' of estimating the probability of size transition. The midpoint method
#' (\code{"midpoint"}) refers to the method in which the probability is
#' estimated by first estimating the probability associated with transition from
#' the exact size at the midpoint of the size class using the corresponding
#' probability density function, and then multiplying that value by the bin
#' width of the size class. Doak et al. 2021 (Ecological Monographs) noted that
#' this method can produce biased results, with total size transitions
#' associated with a specific size not totaling to 1.0 and even specific size
#' transition probabilities capable of being estimated at values greater than
#' 1.0. The alternative and default method, \code{"CDF"}, uses the corresponding
#' cumulative density function to estimate the probability of size transition as
#' the cumulative probability of size transition at the greater limit of the
#' size class minus the cumulative probability of size transition at the lower
#' limit of the size class. The latter method avoids this bias. Note, however,
#' that both methods are exact and unbiased for the Poisson and negative
#' binomial distributions.
#' 
#' Under the Gaussian and gamma size distributions, the number of estimated
#' parameters may differ between the two \code{ipm_method} settings. Because
#' the midpoint method has a tendency to incorporate upward bias in the
#' estimation of size transition probabilities, it is more likely to yield non-
#' zero values when the true probability is extremely close to 0. This will
#' result in the \code{summary.lefkoMat} function yielding higher numbers of
#' estimated parameters than the \code{ipm_method = "CDF"} yields in some cases.
#' 
#' Using the \code{err_check} option will produce a matrix of 6 columns, each
#' characterizing a different vital rate. The product of each row yields an
#' element in the associated \code{$U} matrix. The number and order of elements
#' in each column of this matrix matches the associated matrix in column vector
#' format. Use of this option is generally for the purposes of debugging code.
#' 
#' @seealso \code{\link{flefko3}()}
#' @seealso \code{\link{flefko2}()}
#' @seealso \code{\link{fleslie}()}
#' @seealso \code{\link{rlefko3}()}
#' @seealso \code{\link{rlefko2}()}
#' @seealso \code{\link{rleslie}()}
#' 
#' @examples
#' \donttest{
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
#' minima <- c(0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
#' maxima <- c(NA, 1, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
#'   NA, NA, NA, NA, NA)
#' binvec <- c(0, 4.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
#'   0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
#' 
#' lathframeln <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
#'   propstatus = propvector, minage = minima, maxage = maxima)
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
#'   indiv = "individ", patch = "patchid", year = "year2", age = "obsage",
#'   year.as.random = TRUE, patch.as.random = TRUE, show.model.tables = TRUE,
#'   quiet = TRUE)
#' 
#' # Here we use supplemental() to provide overwrite and reproductive info
#' lathsupp2 <- supplemental(stage3 = c("Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "rep", "rep"),
#'   givenrate = c(0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 3, 3), stageframe = lathframeln, historical = FALSE)
#' 
#' lathmat2age <- aflefko2(year = "all", patch = "all", 
#'   stageframe = lathframeln, modelsuite = lathmodelsln2, data = lathvertln,
#'   supplement = lathsupp2, patchcol = "patchid",
#'   yearcol = "year2", year.as.random = FALSE, patch.as.random = FALSE,
#'   final_age = 2, continue = TRUE, reduce = FALSE)
#' 
#' summary(lathmat2age)
#' 
#' }
#' @export
aflefko2 <- function(year = "all", patch = "all", stageframe, supplement = NULL,
  repmatrix = NULL, overwrite = NULL, data = NA, modelsuite = NA,
  surv_model = NA, obs_model = NA, size_model = NA, sizeb_model = NA,
  sizec_model = NA, repst_model = NA, fec_model = NA, jsurv_model = NA,
  jobs_model = NA, jsize_model = NA, jsizeb_model = NA, jsizec_model = NA,
  jrepst_model = NA, jmatst_model = NA, paramnames = NA, inda = NULL, indb = NULL, indc = NULL,
  surv_dev = 0, obs_dev = 0, size_dev = 0, sizeb_dev = 0, sizec_dev = 0,
  repst_dev = 0, fec_dev = 0, jsurv_dev = 0, jobs_dev = 0, jsize_dev = 0,
  jsizeb_dev = 0, jsizec_dev = 0, jrepst_dev = 0, jmatst_dev = 0, density = NA, repmod = 1,
  yearcol = NA, patchcol = NA, year.as.random = FALSE, patch.as.random = FALSE,
  random.inda = FALSE, random.indb = FALSE, random.indc = FALSE, final_age = 10,
  continue = TRUE, randomseed = NA, negfec = FALSE, ipm_method = "CDF",
  reduce = FALSE, err_check = FALSE, exp_tol = 700, theta_tol = 100000000) {
  
  indanames <- indbnames <- indcnames <- NULL
  
  if (all(is.na(modelsuite)) & all(is.na(paramnames))) {
    warning("Function may not work properly without a dataframe of model parameters or
      equivalents supplied either through the modelsuite option or through the paramnames
      input parameter.", call. = FALSE)
  } else if (!all(is.na(modelsuite))) {
    paramnames <- modelsuite$paramnames
    yearcol <- paramnames$modelparams[which(paramnames$mainparams == "year2")]
    patchcol <- paramnames$modelparams[which(paramnames$mainparams == "patch")]
  }
  
  ipm_method <- tolower(ipm_method)
  if (length(grep("mi", ipm_method)) > 0) {
    ipm_method <- "midpoint"
  } else if (length(grep("cd", ipm_method)) > 0) {
    ipm_method <- "cdf"
  } else {
    stop("Option ipm_method not recognized.", call. = FALSE)
  }
  
  if (all(is.na(data))) {
    stop("Need original vertical dataset to set proper limits on year and patch.", 
      call. = FALSE)
  }
  if (!any(class(data) == "data.frame")) {
    stop("Need original vertical dataset used in modeling to proceed.",
      call. = FALSE)
  }
  if (!any(class(data) == "hfvdata")) {
    warning("Dataset used as input is not of class hfvdata. Will assume that the
      dataset has been formatted equivalently.", call. = FALSE)
  }
  
  stageframe_vars <- c("stage", "size", "size_b", "size_c", "min_age", "max_age",
    "repstatus", "obsstatus", "propstatus", "immstatus", "matstatus", "indataset",
    "binhalfwidth_raw", "sizebin_min", "sizebin_max", "sizebin_center",
    "sizebin_width", "binhalfwidthb_raw", "sizebinb_min", "sizebinb_max",
    "sizebinb_center", "sizebinb_width", "binhalfwidthc_raw", "sizebinc_min",
    "sizebinc_max", "sizebinc_center", "sizebinc_width", "group", "comments")
  if (any(!is.element(names(stageframe), stageframe_vars))) {
    stop("Please use properly formatted stageframe as input.", call. = FALSE)
  }
  
  if (is.character(yearcol)) {
    choicevar <- which(names(data) == yearcol);
    mainyears <- sort(unique(data[,choicevar]))
  } else if (is.numeric(yearcol)) {
    mainyears <- sort(unique(data[, yearcol]));
  } else {
    stop("Need appropriate year column designation.", call. = FALSE)
  }
  
  if (any(is.character(year))) {
    if (is.element("all", tolower(year))) {
      year <- mainyears
    } else {
      stop("Year designation not recognized.", call. = FALSE)
    }
  }
  
  if (length(year) == 0 | all(is.na(year) == TRUE) | any(is.na(year))) {
    stop("This function cannot proceed without a specific occasion, or a suite of
      occasions, designated via the year option. NA entries are not allowed.",
      call. = FALSE)
  }
  
  if (all(is.na(patch)) & !is.na(patchcol)) {
    warning("Matrix creation may not proceed properly without input in the patch
      option when using a modelsuite in which patch is designated.",
      call. = FALSE)
  }
  
  if (is.character(patchcol) & patchcol != "none") {
    choicevar <- which(names(data) == patchcol);
    mainpatches <- sort(unique(as.character(data[,choicevar])))
  } else if (is.numeric(patchcol)) {
    mainpatches <- sort(unique(as.character(data[, patchcol])));
  } else {
    mainpatches <- NA
  }
  
  if (any(is.character(patch))) {
    if (is.element("all", tolower(patch))) {
      patch <- mainpatches
    } else if (!is.element(patch, mainpatches)) {
      stop("Patch designation not recognized.", call. = FALSE)
    }
  }
  
  if (!is.null(inda)) {
    if (!is.numeric(inda) & !random.inda) {
      stop("Individual covariate vector a must be numeric if not set to random.",
        call. = FALSE)
    }
    
    if (!is.element(length(inda), c(1, 2, length(year)))) {
      stop("Individual covariate vector a must be empty, or include 1, 2, or as
        many elements as occasions to be modeled.", call. = FALSE)
    }
    
    if (random.inda) {
      indacol <- paramnames$modelparams[which(paramnames$mainparams == "indcova2")]
      if (indacol == "none") {
        stop("Individual covariate a not recognized in the modelsuite", call. = FALSE)
      }
      
      indacol <- which(names(data) == indacol)
      
      if (length(indacol) > 0) {
        indanames <- sort(unique(data[, indacol]))
      } else {
        stop("Individual covariate a not recognized in the modelsuite", call. = FALSE)
      }
      
      if (any(!is.element(inda, indanames))) {
        stop("Entered value for individual covariate a does not exist in the data.",
          call. = FALSE)
      }
      
      if (length(inda) == 1) {
        r1.inda <- rep(as.character(inda), length(year))
        r2.inda <- rep(as.character(inda), length(year))
      } else if (length(inda) == 2 & length(year) != 2) {
        r1.inda <- rep(as.character(inda[1]), length(year))
        r2.inda <- rep(as.character(inda[2]), length(year))
      } else if (length(inda) == length(year)) {
        r2.inda <- as.character(inda)
        r1.inda <- c("none", r2.inda[1:(length(inda) - 1)])
      }
      
      f1.inda <- rep(0, length(year))
      f2.inda <- rep(0, length(year))
      
    } else {
      indanames <- c(0)
      
      if (length(inda) == 1) {
        f1.inda <- rep(inda, length(year))
        f2.inda <- rep(inda, length(year))
      } else if (length(inda) == 2 & length(year) != 2) {
        f1.inda <- rep(inda[1], length(year))
        f2.inda <- rep(inda[2], length(year))
      } else if (length(inda) == length(year)) {
        f2.inda <- inda
        f1.inda <- c(0, f2.inda[1:(length(inda) - 1)])
      }
      r2.inda <- rep("none", length(year))
      r1.inda <- rep("none", length(year))
    }
  } else {
    indanames <- c(0)
    
    f1.inda <- rep(0, length(year))
    f2.inda <- rep(0, length(year))
    r2.inda <- rep("none", length(year))
    r1.inda <- rep("none", length(year))
  }
  
  if (!is.null(indb)) {
    if (!is.numeric(indb) & !random.indb) {
      stop("Individual covariate vector b must be numeric if not set to random.",
        call. = FALSE)
    }
    
    if (!is.element(length(indb), c(1, 2, length(year)))) {
      stop("Individual covariate vector b must be empty, or include 1, 2, or as
        many elements as occasions to be modeled.", call. = FALSE)
    }
    
    if (random.indb) {
      indbcol <- paramnames$modelparams[which(paramnames$mainparams == "indcovb2")]
      if (indbcol == "none") {
        stop("Individual covariate b not recognized in the modelsuite", call. = FALSE)
      }
      
      indbcol <- which(names(data) == indbcol)
      
      if (length(indbcol) > 0) {
        indbnames <- sort(unique(data[, indbcol]))
      } else {
        stop("Individual covariate b not recognized in the modelsuite", call. = FALSE)
      }
      
      if (any(!is.element(indb, indbnames))) {
        stop("Entered value for individual covariate b does not exist in the data.",
          call. = FALSE)
      }
      
      if (length(indb) == 1) {
        r1.indb <- rep(as.character(indb), length(year))
        r2.indb <- rep(as.character(indb), length(year))
      } else if (length(indb) == 2 & length(year) != 2) {
        r1.indb <- rep(as.character(indb[1]), length(year))
        r2.indb <- rep(as.character(indb[2]), length(year))
      } else if (length(indb) == length(year)) {
        r2.indb <- as.character(indb)
        r1.indb <- c("none", r2.indb[1:(length(indb) - 1)])
      }
      
      f1.indb <- rep(0, length(year))
      f2.indb <- rep(0, length(year))
      
    } else {
      indbnames <- c(0)
      
      if (length(indb) == 1) {
        f1.indb <- rep(indb, length(year))
        f2.indb <- rep(indb, length(year))
      } else if (length(indb) == 2 & length(year) != 2) {
        f1.indb <- rep(indb[1], length(year))
        f2.indb <- rep(indb[2], length(year))
      } else if (length(indb) == length(year)) {
        f2.indb <- indb
        f1.indb <- c(0, f2.indb[1:(length(indb) - 1)])
      }
      r2.indb <- rep("none", length(year))
      r1.indb <- rep("none", length(year))
    }
  } else {
    indbnames <- c(0)
    
    f1.indb <- rep(0, length(year))
    f2.indb <- rep(0, length(year))
    r2.indb <- rep("none", length(year))
    r1.indb <- rep("none", length(year))
  }
  
  if (!is.null(indc)) {
    if (!is.numeric(indc) & !random.indc) {
      stop("Individual covariate vector c must be numeric if not set to random.",
        call. = FALSE)
    }
    
    if (!is.element(length(indc), c(1, 2, length(year)))) {
      stop("Individual covariate vector c must be empty, or include 1, 2, or as
        many elements as occasions to be modeled.", call. = FALSE)
    }
    
    if (random.indc) {
      indccol <- paramnames$modelparams[which(paramnames$mainparams == "indcovc2")]
      if (indccol == "none") {
        stop("Individual covariate c not recognized in the modelsuite", call. = FALSE)
      }
      
      indccol <- which(names(data) == indccol)
      
      if (length(indccol) > 0) {
        indcnames <- sort(unique(data[, indccol]))
      } else {
        stop("Individual covariate c not recognized in the modelsuite", call. = FALSE)
      }
      
      if (any(!is.element(indc, indcnames))) {
        stop("Entered value for individual covariate c does not exist in the data.",
          call. = FALSE)
      }
      
      if (length(indc) == 1) {
        r1.indc <- rep(as.character(indc), length(year))
        r2.indc <- rep(as.character(indc), length(year))
      } else if (length(indc) == 2 & length(year) != 2) {
        r1.indc <- rep(as.character(indc[1]), length(year))
        r2.indc <- rep(as.character(indc[2]), length(year))
      } else if (length(indc) == length(year)) {
        r2.indc <- as.character(indc)
        r1.indc <- c("none", r2.indc[1:(length(indc) - 1)])
      }
      
      f1.indc <- rep(0, length(year))
      f2.indc <- rep(0, length(year))
      
    } else {
      indcnames <- c(0)
      
      if (length(indc) == 1) {
        f1.indc <- rep(indc, length(year))
        f2.indc <- rep(indc, length(year))
      } else if (length(indc) == 2 & length(year) != 2) {
        f1.indc <- rep(indc[1], length(year))
        f2.indc <- rep(indc[2], length(year))
      } else if (length(indc) == length(year)) {
        f2.indc <- indc
        f1.indc <- c(0, f2.indc[1:(length(indc) - 1)])
      }
      r2.indc <- rep("none", length(year))
      r1.indc <- rep("none", length(year))
    }
  } else {
    indcnames <- c(0)
    
    f1.indc <- rep(0, length(year))
    f2.indc <- rep(0, length(year))
    r2.indc <- rep("none", length(year))
    r1.indc <- rep("none", length(year))
  }
  
  maingroups <- sort(unique(stageframe$group))
  
  if (!all(is.na(density))) {
    if (!all(is.numeric(density))) {
      stop("Density value must be numeric.", call. = FALSE)
    }
    
    if (any(is.na(density))) {
      density[which(is.na(density))] <- 0
    }
  } else {
    density <- 0
  }
  
  if (all(is.na(repmatrix)) & all(is.na(supplement))) {
    warning("Neither supplemental data nor a reproduction matrix have been supplied. 
      All fecundity transitions will be inferred from the stageframe.",
      call. = FALSE)
  } else if (all(is.na(repmatrix)) & any(class(supplement) == "lefkoSD")) {
    checkconv <- supplement$convtype
    
    if (!is.element(3, checkconv)) {
      warning("Supplemental data does not include fecundity information, and a reproduction
        matrix has not been supplied. All fecundity transitions will be inferred from the
        stageframe.", call. = FALSE)
    }
  }
  
  stagenum_init <- dim(stageframe)[1]
  if (!all(is.na(repmatrix))) {
    if (any(class(repmatrix) == "matrix")) {
      if (dim(repmatrix)[1] != stagenum_init | dim(repmatrix)[2] != stagenum_init) {
        stop("The repmatrix provided must be a square matrix with dimensions
          equal to the number of stages in the stageframe.", call. = FALSE)
      }
    }
  }
  
  if (any(!suppressWarnings(!is.na(as.numeric(as.character(stageframe$bin_size_ctr)))))) {
    stop("Function aflefko2() requires size to be numeric rather than categorical.", 
      call. = FALSE)
  }
  
  melchett <- .sf_reassess(stageframe, supplement, overwrite, repmatrix,
    agemat = TRUE, historical = FALSE, format = 1)
  stageframe <- melchett$stageframe
  repmatrix <- melchett$repmatrix
  ovtable <- melchett$ovtable
  
  if (!all(is.na(overwrite)) | !all(is.na(supplement))) {
    if(any(duplicated(ovtable[,1:3]))) {
      stop("Multiple entries with different values for the same stage transition are not allowed
        in the supplemental or overwrite table. If modifying a historical table to perform an
        ahistorical analysis, then this may be due to different given rates of substitutions
        caused by dropping stage at occasion t-1. Please eliminate duplicate transitions.",
        call. = FALSE)
    }
  }
  
  # Next the data frame for the C++-based matrix populator functions
  allstages <- .theoldpizzle(stageframe, ovtable, repmatrix, finalage = final_age, 
    format = 1, style = 2, cont = continue)
  
  maxsize <- max(c(allstages$size3, allstages$size2n, allstages$size2o), na.rm = TRUE)
  maxsizeb <- max(c(allstages$sizeb3, allstages$sizeb2n, allstages$sizeb2o), na.rm = TRUE)
  maxsizec <- max(c(allstages$sizec3, allstages$sizec2n, allstages$sizec2o), na.rm = TRUE)
  
  allstages <- allstages[(which(allstages$aliveandequal != -1)),]
  
  if (class(modelsuite) == "lefkoMod") {
    if(is.na(surv_model)) {surv_model <- modelsuite$survival_model}
    if(is.na(obs_model)) {obs_model <- modelsuite$observation_model}
    if(is.na(size_model)) {size_model <- modelsuite$size_model}
    if(is.na(sizeb_model)) {sizeb_model <- modelsuite$sizeb_model}
    if(is.na(sizec_model)) {sizec_model <- modelsuite$sizec_model}
    if(is.na(repst_model)) {repst_model <- modelsuite$repstatus_model}
    if(is.na(fec_model)) {fec_model <- modelsuite$fecundity_model}
    
    if(is.na(jsurv_model)) {jsurv_model <- modelsuite$juv_survival_model}
    if(is.na(jobs_model)) {jobs_model <- modelsuite$juv_observation_model}
    if(is.na(jsize_model)) {jsize_model <- modelsuite$juv_size_model}
    if(is.na(jsizeb_model)) {jsizeb_model <- modelsuite$juv_sizeb_model}
    if(is.na(jsizec_model)) {jsizec_model <- modelsuite$juv_sizec_model}
    if(is.na(jrepst_model)) {jrepst_model <- modelsuite$juv_reproduction_model}
    if(is.na(jmatst_model)) {jmatst_model <- modelsuite$juv_maturity_model}
  }
  
  if (is.na(randomseed)) {
    set.seed(NULL)
  } else {
    set.seed(randomseed)
  }
  
  surv_proxy <- .modelextract(surv_model, paramnames, mainyears, mainpatches, 
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  obs_proxy <- .modelextract(obs_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  sigma <- 0
  rvarssummed <- 0
  sizedist <- 1
  
  size_proxy <- .modelextract(size_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  if (size_proxy$family == "poisson") {
    sizedist <- 0
    if (!all(is.na(size_proxy$variances))) {
      if (is.element("vcov", names(size_proxy$variances))) {
        rvarssummed <- sum(size_proxy$variances[,"vcov"])
      } else if (is.element("cond", names(size_proxy$variances))) {
        rvarssummed <- sum(unlist(lapply(c(1:length(size_proxy$variances[["cond"]])), function(X) {
          stddev.est <- attr(size_proxy$variances[["cond"]][[X]], "stddev")
          return(stddev.est * stddev.est)
        })))
      }
    } else {
      rvarssummed <- 0
    }
  } else if (size_proxy$family == "gaussian") {
    sizedist <- 2
    sigma <- size_proxy$sigma
  } else if (size_proxy$family == "gamma") {
    sizedist <- 3
    
  } else {
    sizedist <- 1
  }
  
  sigmab <- 0
  rvarssummedb <- 0
  sizebdist <- 1
  
  sizeb_proxy <- .modelextract(sizeb_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  if (sizeb_proxy$family == "poisson") {
    sizebdist <- 0
    if (!all(is.na(sizeb_proxy$variances))) {
      if (is.element("vcov", names(sizeb_proxy$variances))) {
        rvarssummedb <- sum(sizeb_proxy$variances[,"vcov"])
      } else if (is.element("cond", names(sizeb_proxy$variances))) {
        rvarssummedb <- sum(unlist(lapply(c(1:length(sizeb_proxy$variances[["cond"]])), function(X) {
          stddev.est <- attr(sizeb_proxy$variances[["cond"]][[X]], "stddev")
          return(stddev.est * stddev.est)
        })))
      }
    } else {
      rvarssummedb <- 0
    }
  } else if (sizeb_proxy$family == "gaussian") {
    sizebdist <- 2
    sigmab <- sizeb_proxy$sigma
  } else if (sizeb_proxy$family == "gamma") {
    sizebdist <- 3
    
  } else {
    sizebdist <- 1
  }
  
  sigmac <- 0
  rvarssummedc <- 0
  sizecdist <- 1
  
  sizec_proxy <- .modelextract(sizec_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  if (sizec_proxy$family == "poisson") {
    sizecdist <- 0
    if (!all(is.na(sizec_proxy$variances))) {
      if (is.element("vcov", names(sizec_proxy$variances))) {
        rvarssummedc <- sum(sizec_proxy$variances[,"vcov"])
      } else if (is.element("cond", names(sizec_proxy$variances))) {
        rvarssummedc <- sum(unlist(lapply(c(1:length(sizec_proxy$variances[["cond"]])), function(X) {
          stddev.est <- attr(sizec_proxy$variances[["cond"]][[X]], "stddev")
          return(stddev.est * stddev.est)
        })))
      }
    } else {
      rvarssummedc <- 0
    }
  } else if (sizec_proxy$family == "gaussian") {
    sizecdist <- 2
    sigmac <- sizec_proxy$sigma
  } else if (sizec_proxy$family == "gamma") {
    sizecdist <- 3
    
  } else {
    sizecdist <- 1
  }
  
  repst_proxy <- .modelextract(repst_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  fec_proxy <- .modelextract(fec_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  if (fec_proxy$family == "poisson") {
    fecdist <- 0
  } else if (fec_proxy$family == "gaussian") {
    fecdist <- 2
  } else if (fec_proxy$family == "gamma") {
    fecdist <- 3
  } else {
    fecdist <- 1
  }
  
  jsurv_proxy <- .modelextract(jsurv_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  jobs_proxy <- .modelextract(jobs_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  jsigma <- 0
  jrvarssummed <- 0
  
  jsize_proxy <- .modelextract(jsize_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  if (jsize_proxy$family == "poisson") {
    if (!all(is.na(jsize_proxy$variances))) {
      if (is.element("vcov", names(jsize_proxy$variances))) {
        jrvarssummed <- sum(jsize_proxy$variances[,"vcov"])
      } else if (is.element("cond", names(jsize_proxy$variances))) {
        jrvarssummed <- sum(unlist(lapply(c(1:length(jsize_proxy$variances[["cond"]])), function(X) {
          stddev.est <- attr(jsize_proxy$variances[["cond"]][[X]], "stddev")
          return(stddev.est * stddev.est)
        })))
      }
    } else {
      jrvarssummed <- 0
    }
  } else if (jsize_proxy$family == "gaussian") {
    jsigma <- jsize_proxy$sigma
  }
  
  jsigmab <- 0
  jrvarssummedb <- 0
  
  jsizeb_proxy <- .modelextract(jsizeb_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  if (jsizeb_proxy$family == "poisson") {
    if (!all(is.na(jsizeb_proxy$variances))) {
      if (is.element("vcov", names(jsizeb_proxy$variances))) {
        jrvarssummedb <- sum(jsizeb_proxy$variances[,"vcov"])
      } else if (is.element("cond", names(jsizeb_proxy$variances))) {
        jrvarssummedb <- sum(unlist(lapply(c(1:length(jsizeb_proxy$variances[["cond"]])), function(X) {
          stddev.est <- attr(jsizeb_proxy$variances[["cond"]][[X]], "stddev")
          return(stddev.est * stddev.est)
        })))
      }
    } else {
      jrvarssummedb <- 0
    }
  } else if (jsizeb_proxy$family == "gaussian") {
    jsigmab <- jsizeb_proxy$sigma
  }
  
  jsigmac <- 0
  jrvarssummedc <- 0
  
  jsizec_proxy <- .modelextract(jsizec_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  if (jsizec_proxy$family == "poisson") {
    if (!all(is.na(jsizec_proxy$variances))) {
      if (is.element("vcov", names(jsizec_proxy$variances))) {
        jrvarssummedc <- sum(jsizec_proxy$variances[,"vcov"])
      } else if (is.element("cond", names(jsizec_proxy$variances))) {
        jrvarssummedc <- sum(unlist(lapply(c(1:length(jsizec_proxy$variances[["cond"]])), function(X) {
          stddev.est <- attr(jsizec_proxy$variances[["cond"]][[X]], "stddev")
          return(stddev.est * stddev.est)
        })))
      }
    } else {
      jrvarssummedc <- 0
    }
  } else if (jsizec_proxy$family == "gaussian") {
    jsigmac <- jsizec_proxy$sigma
  }
  
  jrepst_proxy <- .modelextract(jrepst_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  jmatst_proxy <- .modelextract(jmatst_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  # This creates a list of pop, patch, and year in order of matrix
  if (!all(is.na(patch))) {
    listofyears <- apply(as.matrix(patch), 1, function(X) {
      output <- cbind.data.frame("1", X, as.matrix(year), stringsAsFactors = FALSE);
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    listofyears <- do.call(rbind.data.frame, listofyears)
    listofyears$poporder <- 1
    listofyears$patchorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainpatches == listofyears$patch[X])})
    listofyears$yearorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainyears == listofyears$year2[X])})
    
  } else {
    
    listofyears <- cbind.data.frame("1", "1", as.matrix(year), stringsAsFactors = FALSE)
    names(listofyears) <- c("pop", "patch", "year2")
    
    listofyears$poporder <- 1
    listofyears$patchorder <- 1
    listofyears$yearorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainyears == listofyears$year2[X])})
  }
  
  yearlist <- split(listofyears, seq(nrow(listofyears)))
  
  madsexmadrigal <- lapply(yearlist, .jerzeibalowski, allstages, stageframe, 4,
    surv_proxy, obs_proxy, size_proxy, sizeb_proxy, sizec_proxy, repst_proxy,
    fec_proxy, jsurv_proxy, jobs_proxy, jsize_proxy, jsizeb_proxy, jsizec_proxy,
    jrepst_proxy, jmatst_proxy, f2.inda, f1.inda, f2.indb, f1.indb, f2.indc,
    f1.indc, r2.inda, r1.inda, r2.indb, r1.indb, r2.indc, r1.indc,
    c(surv_dev, obs_dev, size_dev, sizeb_dev, sizec_dev, repst_dev, fec_dev,
      jsurv_dev, jobs_dev, jsize_dev, jsizeb_dev, jsizec_dev, jrepst_dev, jmatst_dev),
    density, repmod, c(rvarssummed, sigma, rvarssummedb, sigmab, rvarssummedc,
      sigmac, jrvarssummed, jsigma, jrvarssummedb, jsigmab, jrvarssummedc, jsigmac),
    maxsize, maxsizeb, maxsizec, final_age, sizedist, sizebdist, sizecdist,
    fecdist, negfec, exp_tol, theta_tol, ipm_method)
  
  a_list <- lapply(madsexmadrigal, function(X) {X$A})
  u_list <- lapply(madsexmadrigal, function(X) {X$U})
  f_list <- lapply(madsexmadrigal, function(X) {X$F})
  
  if (err_check) {
    out_list <- lapply(madsexmadrigal, function(X) {
        aaa <- X$out
        colnames(aaa) <- c("surv", "obs", "repst", "sizea", "sizeb", "sizec", "matst")
        
        return(aaa)
      }
    )
  }
  
  ahstages <- stageframe[1:(dim(stageframe)[1] - 1),]
  
  agestages3 <- ahstages[rep(seq_len(nrow(ahstages)), (final_age + 1)), c(1,2)]
  agestages2 <- rep(c(0:final_age), each = nrow(ahstages))
  agestages <- cbind.data.frame(agestages3, agestages2)
  names(agestages) <- c("stage_id", "stage", "age")
  
  qcoutput1 <- NA
  qcoutput2 <- NA
  
  totalutransitions <- sum(unlist(lapply(u_list, function(X) {length(which(X != 0))})))
  totalftransitions <- sum(unlist(lapply(f_list, function(X) {length(which(X != 0))})))
  totalmatrices <- length(u_list)
  
  qcoutput1 <- c(totalutransitions, totalftransitions, totalmatrices)
  
  if (is.element("qc", names(modelsuite))) {qcoutput2 <- modelsuite$qc}
  
  if (reduce == TRUE) {
    drops <- .reducer2(a_list, u_list, f_list, agestages)
    
    a_list <- drops$A
    u_list <- drops$U
    f_list <- drops$F
    agestages <- drops$ahstages
  }
  
  if (!err_check) {
    output <- list(A = a_list, U = u_list, F = f_list, hstages = NA,
      agestages = agestages, ahstages = ahstages, labels = listofyears[,c(1:3)],
      matrixqc = qcoutput1, modelqc = qcoutput2)
  } else {
    output <- list(A = a_list, U = u_list, F = f_list, hstages = NA,
      agestages = agestages, ahstages = ahstages, labels = listofyears[,c(1:3)],
      matrixqc = qcoutput1, modelqc = qcoutput2, prob_out = out_list)
  }
  
  class(output) <- "lefkoMat"
  
  return(output)
}

#' Create Function-based Age-based (Leslie) Matrix Projection Model
#'
#' Function \code{fleslie()} returns age-based (Leslie) MPMs corresponding to
#' the patches and occasion times given, including the associated component
#' transition and fecundity matrices, data frames detailing the characteristics
#' of the exact ages corresponding to rows and columns in estimated matrices,
#' and a data frame characterizing the patch and occasion time combinations
#' corresponding to these matrices.
#'
#' @param year A variable corresponding to observation occasion, or a set
#' of such values, given in values associated with the year term used in linear 
#' model development. Defaults to \code{"all"}, in which case matrices will be
#' estimated for all occasions.
#' @param patch A variable designating which patches or subpopulations will have
#' matrices estimated. Defaults to \code{"all"}, but can also be set to specific
#' patch names.
#' @param data The historical vertical demographic data frame used to estimate
#' vital rates (class \code{hfvdata}). The original data frame is required in
#' order to initialize occasions and patches properly.
#' @param modelsuite An optional \code{lefkoMod} object holding the vital rate
#' models. If given, then \code{surv_model}, \code{fec_model},
#' \code{paramnames}, \code{yearcol}, and \code{patchcol} are not required. No
#' models should include size or reproductive status in occasion \emph{t}-1.
#' @param surv_model A linear model predicting survival probability. This can be
#' a model of class \code{glm} or \code{glmer}, and requires a predicted
#' binomial variable under a logit link. If given, then will overwrite any
#' survival probability model given in \code{modelsuite}. This model must have
#' been developed in a modeling exercise testing only the impacts of occasion
#' \emph{t}.
#' @param fec_model A linear model predicting fecundity. This can be a model of
#' class \code{glm}, \code{glmer}, \code{glmmTMB}, \code{zeroinfl},
#' \code{vglm}, \code{lm}, or \code{lmer}. If given, then will overwrite any 
#' fecundity model given in \code{modelsuite}. This model must have been
#' developed in a modeling exercise testing only the impacts of occasion
#' \emph{t}.
#' @param paramnames A dataframe with three columns, the first describing all
#' terms used in linear modeling, the second (must be called \code{mainparams}),
#' showing the general model terms that will be used in matrix creation (users
#' should use \code{\link{modelsearch}()} at least once to see the proper
#' names to be used in this column), and the third showing the equivalent terms
#' used in modeling (must be named \code{modelparams}). Only required if
#' \code{modelsuite} is not supplied.
#' @param start_age The age from which to start the matrix. Defaults to
#' \code{NA}, in which case the earliest age in the dataset is used.
#' @param last_age The final age to use in the matrix. Defaults to \code{NA}, in
#' which case the highest age in the dataset is used.
#' @param fecage_min The minimum age at which reproduction is possible. Defaults
#' to \code{NA}, which is interpreted to mean that fecundity should be assessed
#' starting in the first described age.
#' @param fecage_max The maximum age at which reproduction is possible. Defaults
#' to \code{NA}, which is interpreted to mean that fecundity should be assessed
#' until the final observed age.
#' @param continue A logical value designating whether to allow continued
#' survival of individuals past the final age noted in the stageframe, using the 
#' demographic characteristics of the final age. Defaults to \code{TRUE}.
#' @param inda Can be a single value to use for individual covariate \code{a}
#' in all matrices, a pair of values to use for times \emph{t} and \emph{t}-1 in
#' historical matrices, or a vector of such values corresponding to each
#' occasion in option \code{year}. Defaults to \code{NULL}.
#' @param indb Can be a single value to use for individual covariate \code{b}
#' in all matrices, a pair of values to use for times \emph{t} and \emph{t}-1 in
#' historical matrices, or a vector of such values corresponding to each
#' occasion in option \code{year}. Defaults to \code{NULL}.
#' @param indc Can be a single value to use for individual covariate \code{c}
#' in all matrices, a pair of values to use for times \emph{t} and \emph{t}-1 in
#' historical matrices, or a vector of such values corresponding to each
#' occasion in option \code{year}. Defaults to \code{NULL}.
#' @param surv_dev A numeric value to be added to the y-intercept in the linear
#' model for survival probability.
#' @param fec_dev A numeric value to be added to the y-intercept in the linear
#' model for fecundity.
#' @param density A numeric value indicating density value to use to propagate
#' matrices. Only needed if density is an explanatory term used in linear
#' models. Defaults to \code{NA}.
#' @param repmod A scalar multiplier of fecundity. Defaults to \code{1}.
#' @param yearcol The variable name or column number corresponding to year in
#' occasion \emph{t} in the dataset. Not needed if a \code{modelsuite} is
#' supplied.
#' @param patchcol The variable name or column number corresponding to patch in 
#' the dataset. Not needed if a \code{modelsuite} is supplied.
#' @param agecol The name or column number of the variable coding for age in
#' \code{data}. Defaults to \code{"obsage"}.
#' @param year.as.random A logical term indicating whether coefficients for
#' missing occasions within vital rate models should be estimated as random
#' intercepts. Defaults to \code{FALSE}, in which case missing monitoring
#' occasion coefficients are set to \code{0}.
#' @param patch.as.random A logical term indicating whether coefficients for
#' missing patches within vital rate models should be estimated as random
#' intercepts. Defaults to \code{FALSE}, in which case missing patch
#' coefficients are set to \code{0}.
#' @param random.inda A logical value denoting whether to treat individual
#' covariate \code{a} as a random, categorical variable. Otherwise is treated as
#' a fixed, numeric variable. Defaults to \code{FALSE}.
#' @param random.indb A logical value denoting whether to treat individual
#' covariate \code{b} as a random, categorical variable. Otherwise is treated as
#' a fixed, numeric variable. Defaults to \code{FALSE}.
#' @param random.indc A logical value denoting whether to treat individual
#' covariate \code{c} as a random, categorical variable. Otherwise is treated as
#' a fixed, numeric variable. Defaults to \code{FALSE}.
#' @param randomseed A numeric value used as a seed to generate random estimates
#' for missing occasion and patch coefficients, if either \code{year.as.random}
#' or \code{patch.as.random} is set to \code{TRUE}. Defaults to
#' \code{\link{set.seed}()} default.
#' @param negfec A logical value denoting whether fecundity values estimated to
#' be negative should be reset to \code{0}. Defaults to \code{FALSE}.
#' @param exp_tol A numeric value used to indicate a maximum value to set
#' exponents to in the core kernel to prevent numerical overflow. Defaults to
#' \code{700}.
#' @param theta_tol A numeric value used to indicate a maximum value to theta as
#' used in the negative binomial probability density kernel. Defaults to
#' \code{100000000}, but can be reset to other values during error checking.
#'
#' @return If all inputs are properly formatted, then this function will return
#' an object of class \code{lefkoMat}, which is a list that holds the matrix
#' projection model and all of its metadata. Its structure is a list with the
#' following elements:
#'
#' \item{A}{A list of full projection matrices in order of sorted patches and
#' occasions. All matrices output in R's \code{matrix} class.}
#' \item{U}{A list of survival transition matrices sorted as in \code{A}. All 
#' matrices output in R's \code{matrix} class.}
#' \item{F}{A list of fecundity matrices sorted as in \code{A}. All matrices 
#' output in R's \code{matrix} class.}
#' \item{hstages}{Set to \code{NA} for Leslie MPMs.}
#' \item{agestages}{Set to \code{NA} for Leslie MPMs.}
#' \item{ahstages}{A data frame detailing the characteristics of associated
#' ages, in the form of a modified stageframe including reproduction status.}
#' \item{labels}{A data frame giving the patch and year of each matrix in order.
#' In \code{aflefko2()}, only one population may be analyzed at once, and so
#' \code{pop = NA}}
#' \item{matrixqc}{A short vector describing the number of non-zero elements
#' in \code{U} and \code{F} matrices, and the number of annual matrices.}
#' \item{modelqc}{This is the \code{qc} portion of the modelsuite input.}
#' 
#' @section Notes:
#' Unlike \code{\link{rlefko2}()}, \code{\link{rlefko3}()}, and
#' \code{\link{rleslie}()}, this function does not currently distinguish
#' populations.
#' 
#' This function will yield incorrect estimates if the models utilized
#' incorporate state in occasion \emph{t}-1, or any size or reproductive status
#' terms.
#' 
#' Users may at times wish to estimate MPMs using a dataset incorporating
#' multiple patches or subpopulations, but without discriminating between those
#' patches or subpopulations. Should the aim of analysis be a general MPM that
#' does not distinguish these patches or subpopulations, the \code{patchcol}
#' variable should be set to \code{NA}, which is the default.
#'
#' Input options including multiple variable names must be entered in the order
#' of variables in occasion \emph{t}+1 and \emph{t}. Rearranging the order will
#' lead to erroneous calculations, and may lead to fatal errors.
#'
#' Care should be taken to match the random status of year and patch to the
#' states of those variables within the modelsuite. If they do not match, then
#' they will be treated as zeroes in vital rate estimation.
#' 
#' @seealso \code{\link{flefko3}()}
#' @seealso \code{\link{flefko2}()}
#' @seealso \code{\link{aflefko2}()}
#' @seealso \code{\link{rlefko3}()}
#' @seealso \code{\link{rlefko2}()}
#' @seealso \code{\link{rleslie}()}
#' 
#' @examples
#' \donttest{
#' data(lathyrus)
#' 
#' lathvert_base <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
#'   patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9,
#'   sizeacol = "Volume88", repstracol = "FCODE88", fecacol = "Intactseed88",
#'   deadacol = "Dead1988", censorcol = "Missing1988", censorkeep = NA,
#'   censor = TRUE, NRasRep = TRUE, NOasObs = TRUE)
#' 
#' lathvert_base$feca3 <- round(lathvert_base$feca3)
#' lathvert_base$feca2 <- round(lathvert_base$feca2)
#' lathvert_base$feca1 <- round(lathvert_base$feca1)
#' 
#' lathvert_age <- subset(lathvert_base, firstseen > 1988)
#' 
#' lathmodels2_age <- modelsearch(lathvert_age, historical = FALSE,
#'   approach = "mixed", suite = "cons", bestfit = "AICc&k", age = "obsage",
#'   vitalrates = c("surv", "fec"), fecdist = "poisson", indiv = "individ",
#'   year = "year2", year.as.random = TRUE, show.model.tables = TRUE,
#'   quiet = TRUE)
#' 
#' lathmat2fleslie <- fleslie(year = "all", data = lathvert_age,
#'   modelsuite = lathmodels2_age, yearcol = "year2", fecage_min = 1)
#' 
#' summary(lathmat2fleslie)
#' }
#' @export
fleslie <- function(year = "all", patch = "all", data = NA, modelsuite = NA,
  surv_model = NA, fec_model = NA, paramnames = NA, start_age = NA,
  last_age = NA, fecage_min = NA, fecage_max = NA, continue = TRUE, inda = NULL,
  indb = NULL, indc = NULL, surv_dev = 0, fec_dev = 0, density = NA, repmod = 1,
  yearcol = 0, patchcol = NA, agecol = "obsage", year.as.random = FALSE,
  patch.as.random = FALSE, random.inda = FALSE, random.indb = FALSE,
  random.indc = FALSE, randomseed = NA, negfec = FALSE, exp_tol = 700,
  theta_tol = 100000000) {
  
  indanames <- indbnames <- indcnames <- NULL
  
  if (all(is.na(modelsuite)) & all(is.na(paramnames))) {
    warning("Function may not work properly without a dataframe of model 
      parameters or equivalents supplied either through the modelsuite option or 
      through the paramnames input parameter.", call. = FALSE)
  } else if (!all(is.na(modelsuite))) {
    paramnames <- modelsuite$paramnames
    yearcol <- paramnames$modelparams[which(paramnames$mainparams == "year2")]
    patchcol <- paramnames$modelparams[which(paramnames$mainparams == "patch")]
  }
  
  if (all(is.na(data))) {
    stop("Need original vertical dataset to set proper limits on year and patch.", 
      call. = FALSE)
  }
  if (!any(class(data) == "data.frame")) {
    stop("Need original vertical dataset used in modeling to proceed.",
      call. = FALSE)
  }
  if (!any(class(data) == "hfvdata")) {
    warning("Dataset used as input is not of class hfvdata. Will assume that the
      dataset has been formatted equivalently.", call. = FALSE)
  }
  
  # Now we will check which ages are used
  if (is.character(agecol)) {
    choicevar <- which(names(data) == agecol);
    data$usedobsage <- data[,choicevar]
    mainages <- sort(unique(data[,choicevar]))
  } else if (is.numeric(agecol)) {
    mainages <- sort(unique(data[, agecol]))
    data$usedobsage <- data[,agecol]
  } else {
    stop("Need appropriate age column designation.", call. = FALSE)
  }
  age_limit <- max(mainages) + 1
  
  if (is.na(start_age)) {start_age <- min(mainages, na.rm = TRUE)}
  if (is.na(last_age)) {last_age <- max(mainages, na.rm = TRUE) + 1}
  if (is.na(fecage_min)) {fecage_min <- start_age}
  if (is.na(fecage_max)) {fecage_max <- last_age}
  
  start_age <- as.integer(start_age)
  last_age <- as.integer(last_age)
  fecage_min <- as.integer(fecage_min)
  fecage_max <- as.integer(fecage_max)
  
  if (start_age > age_limit || last_age > age_limit) {
    warning("Entered start_age or last_age is beyond what is found in the dataset.",
      call. = FALSE)
  }
  if (fecage_min > age_limit || fecage_max > age_limit) {
    warning("Entered fecage_min or fecage_max is beyond what is found in the dataset.",
      call. = FALSE)
  }
  
  if (last_age < (start_age + 1)) {
    stop("Please include at least 2 ages, and set last_age to be greater than start_age.",
      call. = FALSE)
  }
  if (fecage_max < fecage_min) {
    stop("Please set fecage_max to be greater than or equal to fecage_min.",
      call. = FALSE)
  }
  
  #Now we'll check time and patch variables
  if (is.character(yearcol)) {
    choicevar <- which(names(data) == yearcol);
    mainyears <- sort(unique(data[,choicevar]))
  } else if (is.numeric(yearcol)) {
    mainyears <- sort(unique(data[, yearcol]));
  } else {
    stop("Need appropriate year column designation.", call. = FALSE)
  }
  
  if (any(is.character(year))) {
    if (is.element("all", tolower(year))) {
      year <- mainyears
    } else {
      stop("Year designation not recognized.", call. = FALSE)
    }
  }
  
  if (length(year) == 0 | all(is.na(year) == TRUE) | any(is.na(year))) {
    stop("This function cannot proceed without a specific occasion, or a suite of
      occasions, designated via the year option. NA entries are not allowed.",
      call. = FALSE)
  }
  
  if (all(is.na(patch)) & !is.na(patchcol)) {
    warning("Matrix creation may not proceed properly without input in the patch
      option when using a modelsuite in which patch is designated.",
      call. = FALSE)
  }
  
  if (is.character(patchcol) & patchcol != "none") {
    choicevar <- which(names(data) == patchcol);
    mainpatches <- sort(unique(as.character(data[,choicevar])))
  } else if (is.numeric(patchcol)) {
    mainpatches <- sort(unique(as.character(data[, patchcol])));
  } else {
    mainpatches <- NA
  }
  
  if (any(is.character(patch))) {
    if (is.element("all", tolower(patch))) {
      patch <- mainpatches
    } else if (!all(is.element(patch, mainpatches))) {
      stop("Patch designation not recognized.", call. = FALSE)
    }
  }
  
  if (!is.null(inda)) {
    if (!is.numeric(inda) & !random.inda) {
      stop("Individual covariate vector a must be numeric if not set to random.",
        call. = FALSE)
    }
    
    if (!is.element(length(inda), c(1, 2, length(year)))) {
      stop("Individual covariate vector a must be empty, or include 1, 2, or as
        many elements as occasions to be modeled.", call. = FALSE)
    }
    
    if (random.inda) {
      indacol <- paramnames$modelparams[which(paramnames$mainparams == "indcova2")]
      if (indacol == "none") {
        stop("Individual covariate a not recognized in the modelsuite", call. = FALSE)
      }
      
      indacol <- which(names(data) == indacol)
      
      if (length(indacol) > 0) {
        indanames <- sort(unique(data[, indacol]))
      } else {
        stop("Individual covariate a not recognized in the modelsuite", call. = FALSE)
      }
      
      if (any(!is.element(inda, indanames))) {
        stop("Entered value for individual covariate a does not exist in the data.",
          call. = FALSE)
      }
      
      if (length(inda) == 1) {
        r1.inda <- rep(as.character(inda), length(year))
        r2.inda <- rep(as.character(inda), length(year))
      } else if (length(inda) == 2 & length(year) != 2) {
        r1.inda <- rep(as.character(inda[1]), length(year))
        r2.inda <- rep(as.character(inda[2]), length(year))
      } else if (length(inda) == length(year)) {
        r2.inda <- as.character(inda)
        r1.inda <- c("none", r2.inda[1:(length(inda) - 1)])
      }
      
      f1.inda <- rep(0, length(year))
      f2.inda <- rep(0, length(year))
      
    } else {
      indanames <- c(0)
      
      if (length(inda) == 1) {
        f1.inda <- rep(inda, length(year))
        f2.inda <- rep(inda, length(year))
      } else if (length(inda) == 2 & length(year) != 2) {
        f1.inda <- rep(inda[1], length(year))
        f2.inda <- rep(inda[2], length(year))
      } else if (length(inda) == length(year)) {
        f2.inda <- inda
        f1.inda <- c(0, f2.inda[1:(length(inda) - 1)])
      }
      r2.inda <- rep("none", length(year))
      r1.inda <- rep("none", length(year))
    }
  } else {
    indanames <- c(0)
    
    f1.inda <- rep(0, length(year))
    f2.inda <- rep(0, length(year))
    r2.inda <- rep("none", length(year))
    r1.inda <- rep("none", length(year))
  }
  
  if (!is.null(indb)) {
    if (!is.numeric(indb) & !random.indb) {
      stop("Individual covariate vector b must be numeric if not set to random.",
        call. = FALSE)
    }
    
    if (!is.element(length(indb), c(1, 2, length(year)))) {
      stop("Individual covariate vector b must be empty, or include 1, 2, or as
        many elements as occasions to be modeled.", call. = FALSE)
    }
    
    if (random.indb) {
      indbcol <- paramnames$modelparams[which(paramnames$mainparams == "indcovb2")]
      if (indbcol == "none") {
        stop("Individual covariate b not recognized in the modelsuite", call. = FALSE)
      }
      
      indbcol <- which(names(data) == indbcol)
      
      if (length(indbcol) > 0) {
        indbnames <- sort(unique(data[, indbcol]))
      } else {
        stop("Individual covariate b not recognized in the modelsuite", call. = FALSE)
      }
      
      if (any(!is.element(indb, indbnames))) {
        stop("Entered value for individual covariate b does not exist in the data.",
          call. = FALSE)
      }
      
      if (length(indb) == 1) {
        r1.indb <- rep(as.character(indb), length(year))
        r2.indb <- rep(as.character(indb), length(year))
      } else if (length(indb) == 2 & length(year) != 2) {
        r1.indb <- rep(as.character(indb[1]), length(year))
        r2.indb <- rep(as.character(indb[2]), length(year))
      } else if (length(indb) == length(year)) {
        r2.indb <- as.character(indb)
        r1.indb <- c("none", r2.indb[1:(length(indb) - 1)])
      }
      
      f1.indb <- rep(0, length(year))
      f2.indb <- rep(0, length(year))
      
    } else {
      indbnames <- c(0)
      
      if (length(indb) == 1) {
        f1.indb <- rep(indb, length(year))
        f2.indb <- rep(indb, length(year))
      } else if (length(indb) == 2 & length(year) != 2) {
        f1.indb <- rep(indb[1], length(year))
        f2.indb <- rep(indb[2], length(year))
      } else if (length(indb) == length(year)) {
        f2.indb <- indb
        f1.indb <- c(0, f2.indb[1:(length(indb) - 1)])
      }
      r2.indb <- rep("none", length(year))
      r1.indb <- rep("none", length(year))
    }
  } else {
    indbnames <- c(0)
    
    f1.indb <- rep(0, length(year))
    f2.indb <- rep(0, length(year))
    r2.indb <- rep("none", length(year))
    r1.indb <- rep("none", length(year))
  }
  
  if (!is.null(indc)) {
    if (!is.numeric(indc) & !random.indc) {
      stop("Individual covariate vector c must be numeric if not set to random.",
        call. = FALSE)
    }
    
    if (!is.element(length(indc), c(1, 2, length(year)))) {
      stop("Individual covariate vector c must be empty, or include 1, 2, or as
        many elements as occasions to be modeled.", call. = FALSE)
    }
    
    if (random.indc) {
      indccol <- paramnames$modelparams[which(paramnames$mainparams == "indcovc2")]
      if (indccol == "none") {
        stop("Individual covariate c not recognized in the modelsuite", call. = FALSE)
      }
      
      indccol <- which(names(data) == indccol)
      
      if (length(indccol) > 0) {
        indcnames <- sort(unique(data[, indccol]))
      } else {
        stop("Individual covariate c not recognized in the modelsuite", call. = FALSE)
      }
      
      if (any(!is.element(indc, indcnames))) {
        stop("Entered value for individual covariate c does not exist in the data.",
          call. = FALSE)
      }
      
      if (length(indc) == 1) {
        r1.indc <- rep(as.character(indc), length(year))
        r2.indc <- rep(as.character(indc), length(year))
      } else if (length(indc) == 2 & length(year) != 2) {
        r1.indc <- rep(as.character(indc[1]), length(year))
        r2.indc <- rep(as.character(indc[2]), length(year))
      } else if (length(indc) == length(year)) {
        r2.indc <- as.character(indc)
        r1.indc <- c("none", r2.indc[1:(length(indc) - 1)])
      }
      
      f1.indc <- rep(0, length(year))
      f2.indc <- rep(0, length(year))
      
    } else {
      indcnames <- c(0)
      
      if (length(indc) == 1) {
        f1.indc <- rep(indc, length(year))
        f2.indc <- rep(indc, length(year))
      } else if (length(indc) == 2 & length(year) != 2) {
        f1.indc <- rep(indc[1], length(year))
        f2.indc <- rep(indc[2], length(year))
      } else if (length(indc) == length(year)) {
        f2.indc <- indc
        f1.indc <- c(0, f2.indc[1:(length(indc) - 1)])
      }
      r2.indc <- rep("none", length(year))
      r1.indc <- rep("none", length(year))
    }
  } else {
    indcnames <- c(0)
    
    f1.indc <- rep(0, length(year))
    f2.indc <- rep(0, length(year))
    r2.indc <- rep("none", length(year))
    r1.indc <- rep("none", length(year))
  }
  
  if (!all(is.na(density))) {
    if (!all(is.numeric(density))) {
      stop("Density value must be numeric.", call. = FALSE)
    }
    
    if (any(is.na(density))) {
      density[which(is.na(density))] <- 0
    }
  } else {
    density <- 0
  }
  
  ahages <- .sf_leslie(min_age = start_age, max_age = last_age,
    min_fecage = fecage_min, max_fecage = fecage_max, cont = continue)
  ahages$min_age <- as.integer(ahages$min_age)
  ahages$max_age <- as.integer(ahages$max_age)
  maingroups <- 0
  
  #Now we will work out the vital rate models
  if (class(modelsuite) == "lefkoMod") {
    if(is.na(surv_model)) {surv_model <- modelsuite$survival_model}
    if(is.na(fec_model)) {fec_model <- modelsuite$fecundity_model}
  }
  
  if (is.na(randomseed)) {
    set.seed(NULL)
  } else {
    set.seed(randomseed)
  }
  
  surv_proxy <- .modelextract(surv_model, paramnames, mainyears, mainpatches, 
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  fec_proxy <- .modelextract(fec_model, paramnames, mainyears, mainpatches,
    maingroups, indanames, indbnames, indcnames, year.as.random,
    patch.as.random, random.inda, random.indb, random.indc, err_check = FALSE)
  
  if (fec_proxy$family == "poisson") {
    fecdist <- 0
  } else if (fec_proxy$family == "gaussian") {
    fecdist <- 2
  } else if (fec_proxy$family == "gamma") {
    fecdist <- 3
  } else {
    fecdist <- 1
  }
   
  # Next we create a list of pops, patches, and years in order of matrix
  if (!all(is.na(patch))) {
    listofyears <- apply(as.matrix(patch), 1, function(X) {
      output <- cbind.data.frame("1", X, as.matrix(year), stringsAsFactors = FALSE);
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    listofyears <- do.call(rbind.data.frame, listofyears)
    listofyears$poporder <- 1
    listofyears$patchorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainpatches == listofyears$patch[X])})
    listofyears$yearorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainyears == listofyears$year2[X])})
    
  } else {
    listofyears <- cbind.data.frame("1", "1", as.matrix(year), stringsAsFactors = FALSE)
    names(listofyears) <- c("pop", "patch", "year2")
    
    listofyears$poporder <- 1
    listofyears$patchorder <- 1
    listofyears$yearorder <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {which(mainyears == listofyears$year2[X])})
  }
  
  # A few extra tidbits required for the core matrix estimator to work
  yearlist <- split(listofyears, seq(nrow(listofyears)))
  
  madsexmadrigal <- lapply(yearlist, .motherbalowski, ahages, surv_proxy,
    fec_proxy, f2.inda, f1.inda, f2.indb, f1.indb, f2.indc, f1.indc, 
    r2.inda, r1.inda, r2.indb, r1.indb, r2.indc, r1.indc, surv_dev, fec_dev,
    density, repmod, last_age, fecdist, negfec, exp_tol, theta_tol)
  
  a_list <- lapply(madsexmadrigal, function(X) {X$A})
  u_list <- lapply(madsexmadrigal, function(X) {X$U})
  f_list <- lapply(madsexmadrigal, function(X) {X$F})
  
  qcoutput1 <- NA
  qcoutput2 <- NA
  
  totalutransitions <- sum(unlist(lapply(u_list, function(X) {length(which(X != 0))})))
  totalftransitions <- sum(unlist(lapply(f_list, function(X) {length(which(X != 0))})))
  totalmatrices <- length(u_list)
  
  qcoutput1 <- c(totalutransitions, totalftransitions, totalmatrices)
  
  if (is.element("qc", names(modelsuite))) {qcoutput2 <- modelsuite$qc}
  
  output <- list(A = a_list, U = u_list, F = f_list, hstages = NA,
    agestages = NA, ahstages = ahages, labels = listofyears[,c(1:3)],
    matrixqc = qcoutput1, modelqc = qcoutput2)
  
  class(output) <- "lefkoMat"
  
  return(output)
}

#' Create Raw Leslie (Age-based) Matrix Projection Model
#'
#' Function \code{rleslie()} returns raw Leslie MPMs, including the
#' associated component transition and fecundity matrices, a data frame
#' describing the ages used, and a data frame describing the population, patch,
#' and occasion time associated with each matrix.
#' 
#' @param data A vertical demographic data frame, with variables corresponding 
#' to the naming conventions in \code{\link{verticalize3}()}.
#' @param start_age The age from which to start the matrix. Defaults to
#' \code{NA}, in which case the earliest age in the dataset is used.
#' @param last_age The final age to use in the matrix. Defaults to \code{NA}, in
#' which case the highest age in the dataset is used.
#' @param continue A logical value designating whether to allow continued
#' survival of individuals past the final age noted in the stageframe, using the 
#' demographic characteristics of the final age. Defaults to \code{TRUE}.
#' @param fecage_min The minimum age at which reproduction is possible. Defaults
#' to \code{NA}, which is interpreted to mean that fecundity should be assessed
#' starting in the first described age.
#' @param fecage_max The maximum age at which reproduction is possible. Defaults
#' to \code{NA}, which is interpreted to mean that fecundity should be assessed
#' until the final observed age.
#' @param alive A vector of names of binomial variables corresponding to status
#' as alive (\code{1}) or dead (\code{0}) in occasions \emph{t}+1 ans \emph{t},
#' respectively.
#' @param repst A vector of names of variables coding reproductive status in
#' occasions \emph{t}+1 and \emph{t}, respectively. Defaults to 
#' \code{c("repstatus3", "repstatus2")}.
#' @param fec A vector of names of variables coding fecundity in occasions
#' \emph{t}+1 and \emph{t}, respectively. Defaults to \code{c("feca3", "feca2")}.
#' @param agecol The name or column number of the variable coding for age in
#' \code{data}. Defaults to \code{"obsage"}.
#' @param year A variable corresponding to observation occasion, or a set
#' of such values, given in values associated with the \code{year} term used in
#' vital rate model development. Can also equal \code{"all"}, in which case
#' matrices will be estimated for all occasion times. Defaults to \code{"all"}.
#' @param pop A variable designating which populations will have matrices
#' estimated. Should be set to specific population names, or to \code{"all"} if
#' all populations should have matrices estimated.
#' @param patch A variable designating which patches or subpopulations will have
#' matrices estimated. Should be set to specific patch names, or to \code{"all"}
#' if matrices should be estimated for all patches. Defaults to \code{"all"}.
#' @param yearcol The variable name or column number corresponding to occasion 
#' \emph{t} in the dataset.
#' @param popcol The variable name or column number corresponding to the
#' identity of the population.
#' @param patchcol The variable name or column number corresponding to patch in
#' the dataset.
#' @param indivcol The variable name or column number coding individual
#' identity.
#' @param censor If \code{TRUE}, then data will be removed according to the
#' variable set in \code{censorcol}, such that only data with censor values
#' equal to \code{censorkeep} will remain. Defaults to \code{FALSE}.
#' @param censorcol The variable name or column number denoting the censor
#' status. Only needed if \code{censor = TRUE}.
#' @param censorkeep The value of the censor variable denoting data elements to
#' keep. Defaults to \code{0}.
#' @param fectime An integer indicating whether to estimate fecundity using
#' the variable given for \code{fec} in time \emph{t} (\code{2}) or time
#' \emph{t}+1 (\code{3}).
#' @param fecmod A scalar multiplier for fecundity. Defaults to \code{1.0}.
#' 
#' @return If all inputs are properly formatted, then this function will return
#' an object of class \code{lefkoMat}, which is a list that holds the matrix
#' projection model and all of its metadata. Its structure is a list with the
#' following elements:
#' 
#' \item{A}{A list of full projection matrices in order of sorted populations,
#' patches, and occasions. All matrices output in the \code{matrix} class.}
#' \item{U}{A list of survival transition matrices sorted as in \code{A}. All 
#' matrices output in the \code{matrix} class.}
#' \item{F}{A list of fecundity matrices sorted as in \code{A}. All matrices 
#' output in the \code{matrix} class.}
#' \item{hstages}{A data frame matrix showing the pairing of ahistorical stages
#' used to create historical stage pairs. Set to NA for ahistorical matrices.}
#' \item{agestages}{A data frame showing age-stage pairs. In this function, it
#' is set to NA. Only used in output to function \code{aflefko2}().}
#' \item{ahstages}{A data frame detailing the characteristics of associated
#' ahistorical stages, in the form of a modified stageframe that includes
#' status as an entry stage through reproduction.}
#' \item{labels}{A data frame giving the population, patch, and year of each 
#' matrix in order.}
#' \item{matrixqc}{A short vector describing the number of non-zero elements
#' in \code{U} and \code{F} matrices, and the number of annual matrices.}
#' \item{dataqc}{A vector showing the numbers of individuals and rows in the
#' vertical dataset used as input.}
#'
#' @section Notes:
#' In order to accomodate survival to time \emph{t}+1 in the final year of a
#' study, the maximum age assessed if no input \code{last_age} is provided is
#' one time step past the final described age.
#' 
#' Users may at times wish to estimate MPMs using a dataset incorporating
#' multiple patches or subpopulations. Should the aim of analysis be a general
#' MPM that does not distinguish these patches or subpopulations, the
#' \code{patchcol} variable should be left to \code{NA}, which is the default.
#' Otherwise the variable identifying patch needs to be named.
#'
#' Input options including multiple variable names must be entered in the order
#' of variables in occasion \emph{t}+1 and \emph{t}. Rearranging the order WILL
#' lead to erroneous calculations, and may lead to fatal errors.
#' 
#' @seealso \code{\link{flefko3}()}
#' @seealso \code{\link{flefko2}()}
#' @seealso \code{\link{aflefko2}()}
#' @seealso \code{\link{fleslie}()}
#' @seealso \code{\link{rlefko3}()}
#' @seealso \code{\link{rlefko2}()}
#' 
#' @examples
#' # Cypripedium example
#' data(cypdata)
#' 
#' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
#'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
#'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
#'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
#'   NAas0 = TRUE, NRasRep = TRUE)
#' 
#' cyp_rl <- rleslie(data = cypraw_v1, start_age = 0, last_age = 4, continue = TRUE,
#'   fecage_min = 3, fecage_max = 4, year = "all", pop = NA, patch = "all",
#'   yearcol = "year2", popcol = NA, patchcol = "patchid")
#' cyp_rl
#' 
#' @export
rleslie <- function(data, start_age = NA, last_age = NA, continue = TRUE,
  fecage_min = NA, fecage_max = NA, alive = c("alive3", "alive2", "alive1"),
  repst = c("repstatus3", "repstatus2", "repstatus1"),
  fec = c("feca3", "feca2", "feca1"), agecol = "obsage", year = "all", pop = NA,
  patch = NA, yearcol = NA, popcol = NA, patchcol = NA, indivcol = NA,
  censor = FALSE, censorcol = NA, censorkeep = 0, fectime = 2, fecmod = 1.0) {
  
  popused <- patchused <- yearused <- NULL
  
  if (all(is.na(data))) {
    stop("Need original vertical dataset to proceed.", call. = FALSE)
  }
  
  if (!any(class(data) == "data.frame")) {
    stop("Need original vertical dataset to proceed. This dataset must be in
      historical vertical format.", call. = FALSE)
  }
  
  if (!any(class(data) == "hfvdata")) {
    warning("Dataset used as input is not of class hfvdata. Will assume that the
      dataset has been formatted equivalently.", call. = FALSE)
  }
  
  if (is.character(agecol)) {
    choicevar <- which(names(data) == agecol);
    data$usedobsage <- data[,choicevar]
    mainages <- sort(unique(data[,choicevar]))
  } else if (is.numeric(agecol)) {
    mainages <- sort(unique(data[, agecol]))
    data$usedobsage <- data[,agecol]
  } else {
    stop("Need appropriate age column designation.", call. = FALSE)
  }
  age_limit <- max(mainages) + 1
  
  if (is.na(start_age)) {start_age <- min(mainages, na.rm = TRUE)}
  if (is.na(last_age)) {last_age <- max(mainages, na.rm = TRUE) + 1}
  if (is.na(fecage_min)) {fecage_min <- start_age}
  if (is.na(fecage_max)) {fecage_max <- last_age}
  
  start_age <- as.integer(start_age)
  last_age <- as.integer(last_age)
  fecage_min <- as.integer(fecage_min)
  fecage_max <- as.integer(fecage_max)
  
  if (start_age > age_limit || last_age > age_limit) {
    stop("Entered start_age or last_age is beyond what is found in the dataset.",
      call. = FALSE)
  }
  if (fecage_min > age_limit || fecage_max > age_limit) {
    stop("Entered fecage_min or fecage_max is beyond what is found in the dataset.",
      call. = FALSE)
  }
  
  if (last_age < (start_age + 1)) {
    stop("Please include at least 2 ages, and set last_age to be greater than start_age.",
      call. = FALSE)
  }
  if (fecage_max < fecage_min) {
    stop("Please set fecage_max to be greater than or equal to fecage_min.",
      call. = FALSE)
  }
  
  if (!is.element(fectime, c(2,3))) {
    stop("Option fectime must equal either 2 or 3.", call. = FALSE)
  }
  
  if (censor == TRUE) {
    if(all(is.na(censorcol)) == TRUE) {
      stop("Cannot censor the data without a proper censor variable.", call. = FALSE)
    }
    
    if (all(is.character(censorcol))) {
      if (!all(is.element(censorcol, names(data)))) {
        stop("Censor variable names input for censorcol do not match any
          variable names in the dataset.", call. = FALSE)
      }
    }
    
    censorcolsonly <- data[,censorcol]
    sleeplessnights <- apply(as.matrix(c(1:dim(censorcolsonly)[1])), 1, function(X) {
      crazyvec <- if(is.element(censorkeep, censorcolsonly[X,])) {
        return(X);
      } else {
        return(NA);
      }
    })
    sleeplessnights <- sleeplessnights[!is.na(sleeplessnights)]
    
    data <- data[sleeplessnights,]
  }
  
  if (is.character(yearcol)) {
    choicevar <- which(names(data) == yearcol);
    mainyears <- sort(unique(data[,choicevar]))[-1] # Occasion 1 is unusable, so removed
  } else if (is.numeric(yearcol)) {
    mainyears <- sort(unique(data[, yearcol]))[-1]
  } else {
    stop("Need appropriate year column designation.", call. = FALSE)
  }
  
  if (any(is.character(year))) {
    if (is.element("all", tolower(year))) {
      year <- mainyears
    } else {
      stop("Year designation not recognized.", call. = FALSE)
    }
  }
  
  if (length(year) == 0 | all(is.na(year) == TRUE) | any(is.na(year))) {
    stop("This function cannot proceed without being given a specific year, or a
      suite of years. NA entries are not allowed.", call. = FALSE)
  }
  
  if (!all(is.element(year, mainyears))) {
    stop("Dataset does not contain one or more of the requested years. Note that
      matrices cannot be made for the first year in a historical dataset.",
      call. = FALSE)
  }
  
  if (!all(is.na(pop)) & !all(is.na(patch))) {
    if (is.na(popcol) | is.na(patchcol)) {
      stop("Need population and patch designation variables to proceed.", 
        call. = FALSE)
    }
    
    if (is.element("all", tolower(pop))) {
      if (is.character(popcol)) {popcol <- which(names(data) == popcol)}
      
      pops <- unique(data[,popcol])
    } else {pops <- pop}
    
    if (is.element("all", tolower(patch))) {
      if (is.character(patchcol)) {patchcol <- which(names(data) == patchcol)}
      
      listofpatches <- apply(as.matrix(pops), 1, function(X) {
        patchfolly <- subset(data, popcol == X);
        output <- cbind.data.frame(X, sort(unique(patchfolly[,patchcol])),
          stringsAsFactors = FALSE);
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
      output <- cbind.data.frame(X[1], X[2], sort(unique(checkyrdata[,yearcol]))[-1],
        stringsAsFactors = FALSE);
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
  } else if (all(is.na(pop)) & !all(is.na(patch))) {
    if (!is.na(popcol)) {
      popcol <- NA
    }
    
    if (is.na(patchcol)) {
      stop("Need patch designation variable to proceed.", call. = FALSE)
    }
    
    if (is.element("all", tolower(patch))) {
      if (is.character(patchcol)) {patchcol <- which(names(data) == patchcol)}
      
      patches <- sort(unique(data[,patchcol]))
    } else {patches <- patch}
    
    listofyears <- apply(as.matrix(patches), 1, function(X) {
      checkyrdata <- subset(data, patchcol = X);
      output <- cbind.data.frame("1", X, sort(unique(checkyrdata[,yearcol]))[-1],
        stringsAsFactors = FALSE);
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    if (length(listofyears) > 1) {
      listofyears <- do.call(rbind.data.frame, listofyears)
    } else {
      listofyears <- listofyears[[1]]
    }
  } else if (!all(is.na(pop)) & all(is.na(patch))) {
    if (is.na(popcol)) {
      stop("Need population designation variable to proceed.", call. = FALSE)
    }
    
    if (!is.na(patchcol)) {
      patchcol <- NA
    }
    
    if (is.element("all", tolower(pop))) {
      if (is.character(popcol)) {popcol <- which(names(data) == popcol)}
      
      pops <- unique(data[,popcol])
    } else {pops <- pop}
    
    listofyears <- apply(as.matrix(pops), 1, function(X) {
      checkyrdata <- subset(data, popcol = X);
      output <- cbind.data.frame(X, "1", sort(unique(checkyrdata[,yearcol]))[-1],
        stringsAsFactors = FALSE);
      names(output) <- c("pop", "patch", "year2");
      return(output)
    })
    
    if (length(listofyears) > 1) {
      listofyears <- do.call(rbind.data.frame, listofyears)
    } else {
      listofyears <- listofyears[[1]]
    }
  } else if (all(is.na(pop)) & all(is.na(patch))) {
    if (!is.na(popcol)) {
      popcol <- NA
    }
    if (!is.na(patchcol)) {
      patchcol <- NA
    }
    
    listofyears <- cbind.data.frame("1", "1", year, stringsAsFactors = FALSE)
    names(listofyears) <- c("pop", "patch", "year2")
  }
  
  identifiedyearrows <- which(is.element(listofyears$year2, year))
  if (length(identifiedyearrows) == 0) {
    stop("Cannot recognize input year(s)", call. = FALSE)
  } else {
    listofyears <- listofyears[identifiedyearrows,]
  }
  yearlist <- split(listofyears, seq(nrow(listofyears)))
  
  if (length(alive) > 1) {
    if (length(alive) > 2) {
      data$usedalive1 <- data[,which(names(data) == alive[3])]
    }
    data$usedalive2 <- data[,which(names(data) == alive[2])]
    data$usedalive3 <- data[,which(names(data) == alive[1])]
  } else {
    stop("Function rleslie cannot proceed without the names of at least the 2
      variables coding for status as alive or dead in times t+1 and t.")
  }
  
  if (length(fec) > 1) {
    if (length(fec) > 2) {
      data$usedfec1 <- data[,which(names(data) == fec[3])]
      data$usedfec1[which(is.na(data$usedfec1))] <- 0
    }
    data$usedfec2 <- data[,which(names(data) == fec[2])]
    data$usedfec2[which(is.na(data$usedfec2))] <- 0
    
    data$usedfec3 <- data[,which(names(data) == fec[1])]
    data$usedfec3[which(is.na(data$usedfec3))] <- 0
  } else {
    warning("Function rleslie() requires at least 2 fecundity variables, for
      times t+1 and t (time t-1 may also be input). Failure to include fecundity
      variables leads to matrices composed only of survival transitions.",
      call. = FALSE)
  } 
  
  if (length(repst) > 1) {
    if (length(repst) > 2) {
      data$usedrepst1 <- data[,which(names(data) == repst[3])]
    }
    data$usedrepst2 <- data[,which(names(data) == repst[2])]
    data$usedrepst3 <- data[,which(names(data) == repst[1])]
  } else {
    stop("Function rleslie() cannot proceed without the names of at least the 2 
      variables coding for reproductive status in times t+1 and t.",
      call. = FALSE)
  }
  
  ahages <- .sf_leslie(min_age = start_age, max_age = last_age,
    min_fecage = fecage_min, max_fecage = fecage_max, cont = continue)
  ahages$min_age <- as.integer(ahages$min_age)
  ahages$max_age <- as.integer(ahages$max_age)
  
  madsexmadrigal <- lapply(yearlist, function(X) {
    passed_data <- data
    if (!is.na(X$pop[1]) & !is.na(popcol)) {
      passed_data$popused <- passed_data[,popcol];
      passed_data <- subset(passed_data, popused == X$pop[1]);
    }
    if (!is.na(X$patch[1]) & !is.na(patchcol)) {
      passed_data$patchused <- passed_data[,patchcol];
      passed_data <- subset(passed_data, patchused == X$patch[1]);
    }
    if (!is.na(X$year2[1])) {
      passed_data$yearused <- passed_data[,yearcol];
      passed_data <- subset(passed_data, yearused == X$year2[1]);
    }
    .minorpatrolgroup(MainData = passed_data, StageFrame = ahages,
      fectime = fectime, cont = continue, fec_mod = fecmod)
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
  } else {
    warning("Individual identity variable not provided, affecting quality control output.",
      call. = FALSE)
  }
  qcoutput2 <- c(indivs, dim(data)[1])
  
  output <- list(A = a_list, U = u_list, F = f_list, hstages = NA,
    agestages = NA, ahstages = ahages, labels = listofyears,
    matrixqc = qcoutput1, dataqc = qcoutput2)
  
  class(output) <- "lefkoMat"
  
  return(output)
}

#' Summary of Class "lefkoMat"
#'
#' A function to simplify the viewing of basic information describing the
#' matrices produced through functions \code{\link{flefko3}()},
#' \code{\link{flefko2}()}, \code{\link{rlefko3}()}, \code{\link{rlefko2}()},
#' \code{\link{aflefko2}()}, \code{\link{rleslie}()}, and
#' \code{\link{fleslie}()}.
#' 
#' @param object An object of class \code{lefkoMat}.
#' @param colsums A logical value indicating whether column sums should be shown
#' for U matrices, allowing users to check stage survival probabilities.
#' Defaults to TRUE.
#' @param ... Other parameters.
#' 
#' @return A summary of the object, showing the number of each type of matrix,
#' the number of annual matrices, the number of estimated (non-zero) elements
#' across all matrices and per matrix, the number of unique transitions in the
#' dataset, the number of individuals, and summaries of the column sums of the
#' survival-transition matrices. This function will also yield warnings if any
#' survival-transition matrices include elements outside of the interval [0,1],
#' if any fecundity matrices contain negative elements, and if any matrices
#' include NA values.
#' 
#' @section Notes:
#' Under the Gaussian and gamma size distributions, the number of estimated
#' parameters may differ between the two \code{ipm_method} settings. Because
#' the midpoint method has a tendency to incorporate upward bias in the
#' estimation of size transition probabilities, it is more likely to yield non-
#' zero values when the true probability is extremely close to 0. This will
#' result in the \code{summary.lefkoMat} function yielding higher numbers of
#' estimated parameters than the \code{ipm_method = "CDF"} yields in some cases.
#' 
#' @examples
#' data(cypdata)
#' 
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
#'   "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   propstatus = propvector, immstatus = immvector, indataset = indataset,
#'   binhalfwidth = binvec)
#' 
#' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
#'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
#'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
#'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
#'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
#'   NRasRep = TRUE)
#' 
#' # Here we use supplemental() to provide overwrite and reproductive info
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
#'   stageframe = cypframe_raw, historical = FALSE)
#' 
#' cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added"), supplement = cypsupp2r,
#'   yearcol = "year2", patchcol = "patchid", indivcol = "individ")
#' 
#' summary(cypmatrix2r)
#' 
#' @export
summary.lefkoMat <- function(object, colsums = TRUE, ...) {
  
  matrices <- object
  
  matdim <- dim(matrices$A[[1]])[1]
  
  mqca <- matrices$matrixqc[1]
  mqcb <- matrices$matrixqc[2]
  mqcc <- matrices$matrixqc[3]
  
  totalpops <- length(unique(matrices$labels$pop))
  totalpatches <- length(unique(matrices$labels$patch))
  totalyears <- length(unique(matrices$labels$year2))
  
  if (!all(is.na(matrices$hstages))) {
    histmark <- "historical"
  } else {
    histmark <- "ahistorical"
  }
  
  if (mqcc == 1) {
    writeLines(paste0("\nThis ", histmark, " lefkoMat object contains ", mqcc,
        " matrix."))
  } else {
    writeLines(paste0("\nThis ", histmark, " lefkoMat object contains ", mqcc,
        " matrices."))
  }
  writeLines(paste0("\nEach matrix is square with ", matdim,
    " rows and columns, and a total of ", matdim*matdim, " elements."))
  
  mqac <- mqca / mqcc
  if (mqac != floor(mqac)) mqac <- round(mqac, digits = 3)
  
  if (!all(is.na(mqac))) {
    mqbc <- mqcb / mqcc
    if (mqbc != floor(mqbc)) mqbc <- round(mqbc, digits = 3)
    
    writeLines(paste0("A total of ", mqca, " survival transitions were estimated, with ", 
        mqac, " per matrix."))
    writeLines(paste0("A total of ", mqcb, " fecundity transitions were estimated, with ", 
        mqbc, " per matrix."))
  } else {
    writeLines(paste0("A total of ", mqca, " transitions were estimated, with ", 
        mqac, " per matrix. Positions of survival vs fecundity transitions are not known."))
  }
  
  grammar_pops <- " populations, "
  grammar_patches <- " patches, and "
  grammar_years <- " time steps."
  if (totalpops == 1) grammar_pops <- " population, "
  if (totalpatches == 1) grammar_patches <- " patch, and "
  if (totalyears == 1) grammar_years <- " time step."
  
  writeLines(paste0("This lefkoMat object covers ", totalpops, grammar_pops,
      totalpatches, grammar_patches, totalyears, grammar_years))
  
  if (is.element("dataqc", names(matrices))) {
    dqca <- matrices$dataqc[1]
    dqcb <- matrices$dataqc[2]
    
    if (!is.na(dqca) & !is.na(dqcb)) {
      writeLines(paste0("\nThe dataset contains a total of ", dqca, " unique individuals and ", dqcb, " unique transitions."))
    } else if (!is.na(dqca)) {
      writeLines(paste0("\nThe dataset contains a total of ", dqca, " unique individuals. Number of unique transitions not known."))
    } else {
      writeLines(paste0("\nThe dataset contains a total of ", dqcb, " unique transitions. Number of unique individuals not known."))
    }
  }
  
  if (is.element("modelqc", names(matrices))) {
    if (class(matrices$modelqc) == "data.frame") {
      moqc12 <- matrices$modelqc[1,2]
      moqc22 <- matrices$modelqc[2,2]
      moqc32 <- matrices$modelqc[3,2]
      moqc42 <- matrices$modelqc[4,2]
      moqc52 <- matrices$modelqc[5,2]
      moqc62 <- matrices$modelqc[6,2]
      moqc72 <- matrices$modelqc[7,2]
      moqc82 <- matrices$modelqc[8,2]
      moqc92 <- matrices$modelqc[9,2]
      moqc102 <- matrices$modelqc[10,2]
      moqc112 <- matrices$modelqc[11,2]
      moqc122 <- matrices$modelqc[12,2]
      moqc132 <- matrices$modelqc[13,2]
      moqc142 <- matrices$modelqc[14,2]
      
      moqc13 <- matrices$modelqc[1,3]
      moqc23 <- matrices$modelqc[2,3]
      moqc33 <- matrices$modelqc[3,3]
      moqc43 <- matrices$modelqc[4,3]
      moqc53 <- matrices$modelqc[5,3]
      moqc63 <- matrices$modelqc[6,3]
      moqc73 <- matrices$modelqc[7,3]
      moqc83 <- matrices$modelqc[8,3]
      moqc93 <- matrices$modelqc[9,3]
      moqc103 <- matrices$modelqc[10,3]
      moqc113 <- matrices$modelqc[11,3]
      moqc123 <- matrices$modelqc[12,3]
      moqc133 <- matrices$modelqc[13,3]
      moqc143 <- matrices$modelqc[14,3]
      
      writeLines("\nVital rate modeling quality control:\n")
      
      if (moqc12 > 0) {
        writeLines(paste0("Survival estimated with ", moqc12, " individuals and ", moqc13, " individual transitions."))
      } else {
        writeLines("Survival not estimated.")
      }
      
      if (moqc22 > 0) {
        writeLines(paste0("Observation estimated with ", moqc22, " individuals and ", moqc23, " individual transitions."))
      } else {
        writeLines("Observation probability not estimated.")
      }
      
      if (moqc32 > 0) {
        writeLines(paste0("Primary size estimated with ", moqc32, " individuals and ", moqc33, " individual transitions."))
      } else {
        writeLines("Primary size transition not estimated.")
      }
      
      if (moqc42 > 0) {
        writeLines(paste0("Secondary size estimated with ", moqc42, " individuals and ", moqc43, " individual transitions."))
      } else {
        writeLines("Secondary size transition not estimated.")
      }
      
      if (moqc52 > 0) {
        writeLines(paste0("Tertiary size estimated with ", moqc52, " individuals and ", moqc53, " individual transitions."))
      } else {
        writeLines("Tertiary size transition not estimated.")
      }
      
      if (moqc62 > 0) {
        writeLines(paste0("Reproductive status estimated with ", moqc62, " individuals and ", moqc63, " individual transitions."))
      } else {
        writeLines("Reproduction probability not estimated.")
      }
      
      if (moqc72 > 0) {
        writeLines(paste0("Fecundity estimated with ", moqc72, " individuals and ", moqc73, " individual transitions."))
      } else {
        writeLines("Fecundity not estimated.")
      }
      
      if (moqc82 > 0) {
        writeLines(paste0("Juvenile survival estimated with ", moqc82, " individuals and ", moqc83, " individual transitions."))
      } else {
        writeLines("Juvenile survival not estimated.")
      }
      
      if (moqc92 > 0) {
        writeLines(paste0("Juvenile observation estimated with ", moqc92, " individuals and ", moqc93, " individual transitions."))
      } else {
        writeLines("Juvenile observation probability not estimated.")
      }
      
      if (moqc102 > 0) {
        writeLines(paste0("Juvenile primary size estimated with ", moqc102, " individuals and ", moqc103, " individual transitions."))
      } else {
        writeLines("Juvenile primary size transition not estimated.")
      }
      
      if (moqc112 > 0) {
        writeLines(paste0("Juvenile secondary size estimated with ", moqc112, " individuals and ", moqc113, " individual transitions."))
      } else {
        writeLines("Juvenile secondary size transition not estimated.")
      }
      
      if (moqc122 > 0) {
        writeLines(paste0("Juvenile tertiary size estimated with ", moqc122, " individuals and ", moqc123, " individual transitions."))
      } else {
        writeLines("Juvenile tertiary size transition not estimated.")
      }
      
      if (moqc132 > 0) {
        writeLines(paste0("Juvenile reproduction estimated with ", moqc132, " individuals and ", moqc133, " individual transitions."))
      } else {
        writeLines("Juvenile reproduction probability not estimated.")
      }
      
      if (moqc142 > 0) {
        writeLines(paste0("Juvenile maturity transition estimated with ", moqc142, " individuals and ", moqc143, " individual transitions."))
      } else {
        writeLines("Juvenile maturity transition probability not estimated.")
      }
    }
  }
  
  dethonthetoilet <- apply(as.matrix(c(1:length(matrices$U))), 1, function(X) {
      summary(colSums(matrices$U[[X]]))
    }
  )
  
  sexinthelavatory <- apply(as.matrix(c(1:length(matrices$U))), 1, function(X) {
      summary(colSums(matrices$F[[X]]))
    }
  )
  
  dethintheurinal <- apply(as.matrix(c(1:length(matrices$U))), 1, function(X) {
      any(is.na(matrices$A[[X]]))
    }
  )
  
  if (colsums) {
    writeLines("\nSurvival probability sum check (each matrix represented by column in order):")
    print(dethonthetoilet, digits = 3)
  }
  
  if (max(dethonthetoilet) > 1) {
    warning("Some matrices include stages with survival probability greater than 1.0.", call. = FALSE)
  }
  
  if (min(dethonthetoilet) < 0) {
    warning("Some matrices include stages with survival probability less than 0.0.", call. = FALSE)
  }
  
  if (min(sexinthelavatory) < 0) {
    warning("Some matrices include stages with fecundity less than 0.0.", call. = FALSE)
  }
  
  if (any(dethintheurinal)) {
    warning("Some matrices include NA values.", call. = FALSE)
  }
}

