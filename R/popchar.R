#' Create Overwrite Table for MPM Development
#' 
#' \code{overwrite()} returns a data frame describing which particular
#' transitions within an ahistorical or historical projection matrix to
#' overwrite with either given rates and probabilities, or other estimated
#' transitions.
#'
#' @param stage3 The name of the stage in occasion \emph{t}+1 in the transition
#' to be replaced. Abbreviations for groups of stages are also allowed
#' (see Notes).
#' @param stage2 The name of the stage in occasion \emph{t} in the transition to
#' be replaced. Abbreviations for groups of stages are also allowed (see Notes).
#' @param stage1 The name of the stage in occasion \emph{t}-1 in the transition
#' to be replaced. Only needed if a historical matrix is to be produced.
#' Abbreviations for groups of stages are also allowed (see Notes).
#' @param eststage3 The name of the stage to replace \code{stage3}. Only needed
#' if a transition will be replaced by another estimated transition.
#' @param eststage2 The name of the stage to replace \code{stage2}. Only needed
#' if a transition will be replaced by another estimated transition.
#' @param eststage1 The name of the stage to replace \code{stage1}. Only needed
#' if a transition will be replaced by another estimated transition, and the
#' matrix to be estimated is historical.
#' @param givenrate A fixed rate or probability to replace for the transition
#' described by \code{stage3}, \code{stage2}, and \code{stage1}.
#' @param type A vector denoting the kind of transition between occasions
#' \emph{t} and \emph{t}+1 to be replaced. This should be entered as \code{1},
#' \code{S}, or \code{s} for the replacement of a survival transition; or
#' \code{2}, \code{F}, or \code{f} for the replacement of a fecundity
#' transition. If empty or not provided, then defaults to \code{1} for survival
#' transition.
#' @param type_t12 An optional vector denoting the kind of transition between
#' occasions \emph{t}-1 and \emph{t}. Only necessary if a historical MPM in
#' deVries format is desired. This should be entered as \code{1}, \code{S}, or
#' \code{s} for a survival transition; or \code{2}, \code{F}, or \code{f} for a
#' fecundity transitions. Defaults to \code{1} for survival transition, with
#' impacts only on the construction of deVries-format hMPMs.
#'
#' @return A data frame that puts the above vectors together and can be used as
#' input in \code{\link{flefko3}()}, \code{\link{flefko2}()},
#' \code{\link{rlefko3}()},\code{\link{rlefko2}()}, and
#' \code{\link{aflefko2}()}.
#' 
#' Variables in this data frame include the following:
#' \item{stage3}{Stage at occasion \emph{t}+1 in the transition to be replaced.}
#' \item{stage2}{Stage at occasion \emph{t} in the transition to be replaced.}
#' \item{stage1}{Stage at occasion \emph{t}-1 in the transition to be replaced.}
#' \item{eststage3}{Stage at occasion \emph{t}+1 in the transition to replace
#' the transition designated by \code{stage3}, \code{stage2}, and \code{stage1}.}
#' \item{eststage2}{Stage at occasion \emph{t} in the transition to replace the
#' transition designated by \code{stage3}, \code{stage2}, and \code{stage1}.}
#' \item{eststage1}{Stage at occasion \emph{t}-1 in the transition to replace
#' the transition designated by \code{stage3}, \code{stage2}, and \code{stage1}.}
#' \item{givenrate}{A constant to be used as the value of the transition.}
#' \item{convtype}{Designates whether the transition from occasion \emph{t} to
#' occasion \emph{t}+1 is a survival-transition probability (1) or a fecundity
#' rate (2).}
#' \item{convtype_t12}{Designates whether the transition from occasion
#' \emph{t}-1 to occasion \emph{t} is a survival transition probability (1), a
#' fecundity rate (2).}
#' 
#' @section Notes:
#' Entries in \code{stage3}, \code{stage2}, and \code{stage1} can include
#' abbreviations for groups of stages. Use \code{rep} if all reproductive stages
#' are to be used, \code{nrep} if all mature but non-reproductive stages are to
#' be used, \code{mat} if all mature stages are to be used, \code{immat} if all
#' immature stages are to be used, \code{prop} if all propagule stages are to be
#' used, \code{npr} if all non-propagule stages are to be used, and leave empty
#' or use \code{all} if all stages in stageframe are to be used.
#' 
#' @examples
#' cypover2r <- overwrite(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
#'     "XSm", "Sm"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL"),
#'   eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm"),
#'   eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm"),
#'   givenrate = c(0.1, 0.2, 0.2, 0.2, 0.25, NA, NA, NA),
#'   type = c("S", "S", "S", "S", "S", "S", "S", "S"))
#' 
#' cypover2r
#' 
#' cypover3r <- overwrite(stage3 = c("SD", "SD", "P1", "P1", "P2", "P3", "SL", 
#'     "D", "XSm", "Sm", "D", "XSm", "Sm"),
#'   stage2 = c("SD", "SD", "SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL",
#'     "SL", "SL", "SL"),
#'   stage1 = c("SD", "rep", "SD", "rep", "SD", "P1", "P2", "P3", "P3", "P3",
#'     "SL", "SL", "SL"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", "D", "XSm",
#'     "Sm"),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm",
#'     "XSm", "XSm"),
#'   eststage1 = c(NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm",
#'     "XSm", "XSm"),
#'   givenrate = c(0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.25, NA, NA, NA, NA, NA, NA),
#'   type = c("S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S"))
#' 
#' cypover3r
#' 
#' @export
overwrite <- function(stage3, stage2, stage1 = NA, eststage3 = NA, 
  eststage2 = NA, eststage1 = NA, givenrate = NA, type = NA, type_t12 = NA) {
  
  if (length(stage3) != length(stage2)) {
    stop("All transitions to overwrite require information at least for stage2
      and stage3. These inputs must also be of equal length.", call. = FALSE)
  }
  
  fulllength <- max(length(stage3), length(stage2), length(stage1), length(eststage3), 
    length(eststage2), length(eststage1), length(givenrate), length(type))
  
  if (length(stage3) != fulllength) {
    stop("Please provide all input vectors in the same order.", call. = FALSE)
  }
  
  if (length(stage1) < fulllength) {
    missinglength <- fulllength - length(stage1)
    stage1 <- as.character(append(stage1, rep(NA, missinglength)))
  }
  if (length(eststage3) < fulllength) {
    missinglength <- fulllength - length(eststage3)
    eststage3 <- as.character(append(eststage3, rep(NA, missinglength)))
  }
  if (length(eststage2) < fulllength) {    
    missinglength <- fulllength - length(eststage2)
    eststage2 <- as.character(append(eststage2, rep(NA, missinglength)))
  }
  if (length(eststage1) < fulllength) {
    missinglength <- fulllength - length(eststage1)
    eststage1 <- as.character(append(eststage1, rep(NA, missinglength)))
  }
  if (length(givenrate) < fulllength) {
    missinglength <- fulllength - length(givenrate)
    givenrate <- as.numeric(append(givenrate, rep(NA, missinglength)))
  }
  if (length(type) < fulllength) {
    missinglength <- fulllength - length(type)
    type <- as.character(append(type, rep(NA, missinglength)))
  }
  if (length(type_t12) < fulllength) {
    missinglength <- fulllength - length(type_t12)
    type_t12 <- as.character(append(type_t12, rep(NA, missinglength)))
  }
  
  if(!all(is.na(type))) {
    convtype <- rep(1, length(type))
    convtype[which(type == "F")] <- 2
    convtype[which(type == "f")] <- 2
  } else {
    convtype <- rep(1, length(stage3))
  }
  if(!all(is.na(type_t12))) {
    convtype_t12 <- rep(1, length(type_t12))
    convtype_t12[which(type_t12 == "F")] <- 2
    convtype_t12[which(type_t12 == "f")] <- 2
    convtype_t12[which(type_t12 == "2")] <- 2
  } else {
    convtype_t12 <- rep(1, length(stage3))
  }
  
  fullpack <- cbind.data.frame(stage3, stage2, stage1, eststage3, eststage2,
    eststage1, givenrate, convtype, convtype_t12, stringsAsFactors = FALSE)
  
  return(fullpack)
}

#' Create a Data Frame of Supplemental Data for MPM Development
#' 
#' Function \code{supplemental()} provides all necessary supplemental data for
#' matrix estimation, particularly bringing together data on proxy rates, data
#' to overwrite existing rates, identified reproductive transitions complete,
#' and fecundity multipliers.
#'
#' @param stage3 The name of the stage in occasion \emph{t}+1 in the transition
#' to be replaced. Abbreviations for groups of stages are also usable (see
#' Notes).
#' @param stage2 The name of the stage in occasion \emph{t} in the transition to
#' be replaced. Abbreviations for groups of stages are also usable (see Notes).
#' @param stage1 The name of the stage in occasion \emph{t}-1 in the transition
#' to be replaced. Only needed if a historical matrix is to be produced.
#' Abbreviations for groups of stages are also usable (see Notes).
#' @param eststage3 The name of the stage to replace \code{stage3} in a proxy
#' transition. Only needed if a transition will be replaced by another estimated
#' transition.
#' @param eststage2 The name of the stage to replace \code{stage2} in a proxy
#' transition. Only needed if a transition will be replaced by another estimated
#' transition.
#' @param eststage1 The name of the stage to replace \code{stage1} in a proxy
#' historical transition. Only needed if a transition will be replaced by
#' another estimated transition, and the matrix to be estimated is historical.
#' Stage \code{NotAlive} is also possible for raw hMPMs as a means of handling
#' the prior stage for individuals entering the population in occasion \emph{t}.
#' @param givenrate A fixed rate or probability to replace for the transition
#' described by \code{stage3}, \code{stage2}, and \code{stage1}.
#' @param multiplier A vector of numeric multipliers for fecundity or for proxy
#' transitions. Defaults to \code{1}.
#' @param type A vector denoting the kind of transition between occasions
#' \emph{t} and \emph{t}+1 to be replaced. This should be entered as \code{1},
#' \code{S}, or \code{s} for the replacement of a survival transition; \code{2},
#' \code{F}, or \code{f} for the replacement of a fecundity transition; or
#' \code{3}, \code{R}, or \code{r} for a fecundity multiplier. If empty or not
#' provided, then defaults to \code{1} for survival transition.
#' @param type_t12 An optional vector denoting the kind of transition between
#' occasions \emph{t}-1 and \emph{t}. Only necessary if a historical MPM in
#' deVries format is desired. This should be entered as \code{1}, \code{S}, or
#' \code{s} for a survival transition; or \code{2}, \code{F}, or \code{f} for a
#' fecundity transitions. Defaults to \code{1} for survival transition, with
#' impacts only on the construction of deVries-format hMPMs.
#' @param stageframe The stageframe being used to produce the MPMs in the study.
#' @param historical A logical value indicating whether the MPMs intended will
#' be historical or ahistorical. Defaults to TRUE.
#'
#' @return A data frame of class \code{lefkoSD}. This object can be used as
#' input in \code{\link{flefko3}()}, \code{\link{flefko2}()}, 
#' \code{\link{rlefko3}()}, \code{\link{rlefko2}()}, and 
#' \code{\link{aflefko2}()}.
#' 
#' Variables in this object include the following:
#' \item{stage3}{Stage at occasion \emph{t}+1 in the transition to be replaced.}
#' \item{stage2}{Stage at occasion \emph{t} in the transition to be replaced.}
#' \item{stage1}{Stage at occasion \emph{t}-1 in the transition to be replaced.}
#' \item{eststage3}{Stage at occasion \emph{t}+1 in the transition to replace
#' the transition designated by \code{stage3}, \code{stage2}, and 
#' \code{stage1}.}
#' \item{eststage2}{Stage at occasion \emph{t} in the transition to replace the
#' transition designated by \code{stage3}, \code{stage2}, and \code{stage1}.}
#' \item{eststage1}{Stage at occasion \emph{t}-1 in the transition to replace
#' the transition designated by \code{stage3}, \code{stage2}, and 
#' \code{stage1}.}
#' \item{givenrate}{A constant to be used as the value of the transition.}
#' \item{multiplier}{A multiplier for proxy transitions or for fecundity.}
#' \item{convtype}{Designates whether the transition from occasion \emph{t} to
#' occasion \emph{t}+1 is a survival transition probability (1), a fecundity
#' rate (2), or a fecundity multiplier (3).}
#' \item{convtype_t12}{Designates whether the transition from occasion
#' \emph{t}-1 to occasion \emph{t} is a survival transition probability (1), a
#' fecundity rate (2).}
#' 
#' @section Notes:
#' Negative values are not allowed in \code{givenrate} and \code{multiplier}
#' input.
#' 
#' Fecundity multiplier data supplied via the \code{supplemental()} function
#' acts in the same way as non-zero entries supplied via a reproductive matrix,
#' but gets priority in all matrix creations. Thus, in cases where fecundity
#' multipliers are provided for the same function via the reproductive matrix
#' and function \code{supplemental()}, the latter is used.
#' 
#' Entries in \code{stage3}, \code{stage2}, and \code{stage1} can include
#' abbreviations for groups of stages. Use \code{rep} if all reproductive stages
#' are to be used, \code{nrep} if all mature but non-reproductive stages are to
#' be used, \code{mat} if all mature stages are to be used, \code{immat} if all
#' immature stages are to be used, \code{prop} if all propagule stages are to be
#' used, \code{npr} if all non-propagule stages are to be used, and leave empty
#' or use \code{all} if all stages in stageframe are to be used. Also use
#' \code{groupX} to denote all stages in group X (e.g. \code{group1} will use
#' all stages in the respective stageframe's group 1).
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
#'   fecacol = "Intactseed88", deadacol = "Dead1988",
#'   nonobsacol = "Dormant1988", stageassign = lathframe, stagesize = "sizea",
#'   censorcol = "Missing1988", censorkeep = NA, censor = TRUE)
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
#' ehrlen3mean <- lmean(ehrlen3)
#' ehrlen3mean$A[[1]]
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
#' cyp2mean <- lmean(cypmatrix2r)
#' cyp2mean
#' 
#' @export
supplemental <- function(stage3, stage2, stage1 = NA, eststage3 = NA,
  eststage2 = NA, eststage1 = NA, givenrate = NA, multiplier = 1, type = NA,
  type_t12 = NA, stageframe, historical = TRUE) {
  
  if (all(class(stageframe) != "stageframe")) {
    stop("A regular stageframe, as output from the sf_create() function, is
      required for function supplemental().", call. = FALSE)
  }
  
  if (!is.element("stage", names(stageframe))) {
    stop("Stageframe appears to be modified. Please make sure that a stage
      column exists holding stage names.", call. = FALSE)
  }
  
  if (length(stage3) != length(stage2)) {
    stop("All transitions to overwrite require information at least for stage2
      and stage3. These inputs must also be of equal length.", call. = FALSE)
  }
  
  fulllength <- max(length(stage3), length(stage2), length(stage1),
    length(eststage3), length(eststage2), length(eststage1), length(givenrate),
    length(type))
  
  if (length(stage3) != fulllength) {
    stop("Please provide all input vectors in the same order.", call. = FALSE)
  }
  
  if (length(stage1) < fulllength) {
    missinglength <- fulllength - length(stage1)
    stage1 <- as.character(append(stage1, rep(NA, missinglength)))
  }
  if (length(eststage3) < fulllength) {
    missinglength <- fulllength - length(eststage3)
    eststage3 <- as.character(append(eststage3, rep(NA, missinglength)))
  }
  if (length(eststage2) < fulllength) {    
    missinglength <- fulllength - length(eststage2)
    eststage2 <- as.character(append(eststage2, rep(NA, missinglength)))
  }
  if (length(eststage1) < fulllength) {
    missinglength <- fulllength - length(eststage1)
    eststage1 <- as.character(append(eststage1, rep(NA, missinglength)))
  }
  if (length(givenrate) < fulllength) {
    missinglength <- fulllength - length(givenrate)
    givenrate <- as.numeric(append(givenrate, rep(NA, missinglength)))
  }
  if (length(multiplier) < fulllength) {
    missinglength <- fulllength - length(multiplier)
    multiplier <- as.numeric(append(multiplier, rep(1, missinglength)))
  }
  if (any(is.na(multiplier))) {
    multNAs <- which(is.na(multiplier))
    multiplier[multNAs] <- 1
  }
  if (length(type) < fulllength) {
    missinglength <- fulllength - length(type)
    type <- as.character(append(type, rep(NA, missinglength)))
  }
  if (length(type_t12) < fulllength) {
    missinglength <- fulllength - length(type_t12)
    type_t12 <- as.character(append(type_t12, rep(NA, missinglength)))
  }
  
  ltype <- tolower(type)
  typeall <- unique(ltype)
  if (!all(is.element(typeall, c(NA, "1", "2", "3", "f", "r", "s")))) {
    stop("Variable type must include only 1, 2, 3, s, r, and f. All other
      entries are not allowed.", call. = FALSE)
  }
  ltype_t12 <- tolower(type_t12)
  typeall_t12 <- unique(ltype_t12)
  if (!all(is.element(typeall_t12, c(NA, "1", "2", "f", "s")))) {
    stop("Variable type_t12 must include only 1, 2, s, and f. All other entries
      are not allowed.", call. = FALSE)
  }
  
  convtype <- rep(1, length(type))
  convtype[which(ltype == "2")] <- 2
  convtype[which(ltype == "3")] <- 3
  convtype[which(ltype == "f")] <- 2
  convtype[which(ltype == "r")] <- 3
  
  convtype_t12 <- rep(1, length(type_t12))
  convtype_t12[which(ltype_t12 == "2")] <- 2
  convtype_t12[which(ltype_t12 == "f")] <- 2

  all.stages.sf <- stageframe$stage
  
  all.stages.inp <- unique(c(stage3, stage2, stage1, eststage3, eststage2, eststage1))
  
  mismatches <- !is.element(all.stages.inp, c(all.stages.sf, NA))
  
  if (length(which(mismatches)) > 0) {
    extrastuff <- tolower(all.stages.inp[which(mismatches)])
    
    unique_groups <- unique(stageframe$group)
    group_labels <- apply(as.matrix(c(1:length(unique_groups))), 1, function(X) {
      return(paste0("group", X))
    })
    
    wildcard_list <- c("all", "rep", "nrep", "mat", "immat", "prop", "npr",
      "notalive", group_labels)
    unaccountedfor <- extrastuff[which(!is.element(extrastuff, wildcard_list))]
    
    if (!all(is.element(extrastuff, wildcard_list))) {
      stop(paste("The following stage names input in supplemental() do not match
        the stageframe:", paste(unaccountedfor, collapse = ' ')), call. = FALSE)
    }
  }
  
  if (is.element("notalive", tolower(c(stage1, stage2, stage3, eststage2, eststage3)))) {
    stop("Stage NotAlive is only allowed in the input for eststage1.",
      call. = FALSE)
  }
  
  if (any(multiplier[which(!is.na(givenrate))] != 1)) {
    warning("Multipliers assigned to given rates will be ignored in MPM creation.",
      call. = FALSE)
  }
  
  if (any(givenrate < 0, na.rm = TRUE)) {
    stop("Given rates cannot be negative.", call. = FALSE)
  }
  
  if (any(multiplier < 0, na.rm = TRUE)) {
    stop("Multipliers cannot be negative.", call. = FALSE)
  }
  
  output <- cbind.data.frame(stage3, stage2, stage1, eststage3, eststage2,
    eststage1, givenrate, multiplier, convtype, convtype_t12,
    stringsAsFactors = FALSE)
  
  class(output) <- append(class(output), "lefkoSD")
  
  #Final check
  all12s <- which(convtype != 3)
  all3s <- which(convtype == 3)
  
  if (length(all12s) > 0) {
    givenests <- which(!is.na(eststage3))
    givengivens <- which(!is.na(givenrate))
    givenmults <- which(multiplier != 1)
    
    givens <- unique(union(givenests, union(givengivens, givenmults)))
    if (length(givens) < length(all12s)) {
      stop("Some given rates or proxy transitions do not appear to be given.",
        call. = FALSE)
    }
  }
  
  return(output)
}

#' Test Overdispersion and Zero Inflation in Size and Fecundity Distributions
#' 
#' Function \code{sf_distrib} takes a historically formatted vertical data as
#' input and tests whether size and fecundity data are dispersed according to a
#' Poisson distribution (where mean = variance), and whether the number of 0s
#' exceeds expectations.
#'
#' @param data A historical vertical data file, which is a data frame of class
#' \code{hfvdata}.
#' @param sizea A vector holding the name or column number of the variables
#' corresponding to primary size in occasions *t*+1 and *t*. Input only if
#' \code{sizea} is to be tested.
#' @param sizeb A vector holding the name or column number of the variables
#' corresponding to secondary size in occasions *t*+1 and *t*. Input only if
#' \code{sizeb} is to be tested.
#' @param sizec A vector holding the name or column number of the variables
#' corresponding to tertiary size in occasions *t*+1 and *t*. Input only if
#' \code{sizec} is to be tested.
#' @param obs3 The name or column number of the variable corresponding to
#' observation status in occasion *t+1*. This should be used if observation
#' status will be used as a vital rate to absorb states of size = 0.
#' @param fec A vector holding the names or column numbers of the variables
#' corresponding to in occasions *t*+1 and *t*. Input only if \code{fec} is to
#' be tested.
#' @param repst A vector holding the names or column numbers of the variables
#' corresponding to reproductive status in occasions *t*+1 and *t*. If not
#' provided, then fecundity will be tested without subsetting to only
#' reproductive individuals.
#' @param zisizea A logical value indicating whether to conduct a test of zero
#' inflation in primary size. Defaults to \code{TRUE}.
#' @param zisizeb A logical value indicating whether to conduct a test of zero
#' inflation in secondary size. Defaults to \code{TRUE}.
#' @param zisizec A logical value indicating whether to conduct a test of zero
#' inflation in tertiary size. Defaults to \code{TRUE}.
#' @param zifec A logical value indicating whether to conduct a test of zero
#' inflation in fecundity. Defaults to TRUE.
#' @param fectime An integer indicating whether to treat fecundity as occurring
#' in time *t* (\code{2}) or time *t*+1 (\code{3}). Defaults to \code{2}.
#' @param show.size A logical value indicating whether to show the output for
#' tests of size. Defaults to \code{TRUE}.
#' @param show.fec A logical value indicating whether to show the output for
#' tests of fecundity. Defaults to \code{TRUE}.
#'
#' @return Produces text describing the degree and significance of difference
#' from expected dispersion, and the degree and significance of zero inflation.
#' The tests are chi-squared score tests based on the expectations of 
#' mean = variance, and 0s as abundant as predicted by the value of lambda
#' estimated from the dataset. See van der Broek (1995) for more details.
#' 
#' @section Notes:
#' This function subsets the data in the same way as \code{\link{modelsearch}()}
#' before testing underlying distributions, making the output much more
#' appropriate than a simple analysis of size and fecundity variables in
#' \code{data}.
#' 
#' The specific test used for overdispersion is a chi-squared test of the
#' dispersion parameter estimated using a generalized linear model predicting
#' the response given size in occasion *t*, under a quasi-Poisson distribution.
#' 
#' The specific test used for zero-inflation is the chi-squared test presented
#' in van der Broek (1995).
#' 
#' @examples
#' # Lathyrux example
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
#' # The following will only test fecundity, since size is Gaussian.
#' # Zero-inflation will not be assessed in this example, since 0 values in
#' # fecundity have been excluded in the life history model.
#' 
#' sf_distrib(lathvertln, sizea = c("sizea3", "sizea2"), fec = c("feca3", "feca2"),
#'   repst = c("repstatus3", "repstatus2"), zifec = FALSE)
#' 
#' # Cypripedium example
#' rm(list=ls(all=TRUE))
#' 
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
#' sf_distrib(cypraw_v1, sizea = c("size3added", "size2added"),
#'   fec = c("feca3", "feca2"), repst = c("repstatus3", "repstatus2"),
#'   zisizea = TRUE)
#' 
#' @export
sf_distrib <- function(data, sizea = NA, sizeb = NA, sizec = NA, obs3 = NA,
  fec = NA, repst = NA, zisizea = TRUE, zisizeb = TRUE, zisizec = TRUE,
  zifec = TRUE, fectime = 2, show.size = TRUE, show.fec = TRUE) {
  
  alive3 <- NULL
  
  if (!any(class(data) == "hfvdata")) {
    stop("Function sf_distrib requires an object of class hfvdata as input.",
      call. = FALSE)
  }
  
  if (all(is.na(c(sizea, sizeb, sizec)))) {
    stop("Function sf_distrib requires a size variable, even if only fecundity is
      to be tested. Please designate at least one such variable", call. = FALSE)
  }
  
  if (!all(is.na(fec)) & all(is.na(repst))) {
    warning("Fecundity will be assessed without subsetting to only reproductive
      individuals. If incorrect, then please provide reproductive status
      variables (repst) in time t+1 and t.", call. = FALSE)
  }
  
  sdata <- subset(data, alive3 == 1)
  
  if (length(obs3) > 1) {
    stop("Please input only one column number or name in obs3.", call. = FALSE)
  }
  if (!all(is.na(obs3))) {
    if (is.numeric(obs3)) {
      if (!is.element(obs3, c(1:length(names(data))))) {
        stop("Column number given in obs is not recognized.", call. = FALSE)
      } else {
        sdata <- sdata[which(sdata[,obs3] == 1),]
      }
    } else if (is.character(obs3)) {
      obs3low <- tolower(obs3)
      datanames <- tolower(names(sdata))
      
      obs3proxy <- grep(obs3low, datanames, fixed = TRUE)
      
      if (length(obs3proxy) == 0) {
        stop("Name of obs3 variable does not match any variable in dataset.",
          call. = FALSE)
      } else if (length(obs3proxy) > 1) {
        stop("Obs3 variable name appears to match several variable names in
          dataset.", call. = FALSE)
      }
      
      sdata <- sdata[which(sdata[,obs3proxy] == 1),]
    }
  }
  
  sizeatest <- .knightswhosaynee(sdata, sizea, 1, zisizea, sizea, sizeb, sizec,
    repst, fectime, show.size)
  sizebtest <- .knightswhosaynee(sdata, sizeb, 2, zisizeb, sizea, sizeb, sizec,
    repst, fectime, show.size)
  sizectest <- .knightswhosaynee(sdata, sizec, 3, zisizec, sizea, sizeb, sizec,
    repst, fectime, show.size)
  
  fectest <- .knightswhosaynee(sdata, fec, 4, zifec, sizea, sizeb, sizec, repst,
    fectime, show.fec)
}

#' Internal Test of Size and Fecundity For Overdispersion and Zero-Inflation
#' 
#' This internal function automates the testing of size and fecundity variables
#' for overdispersion and zero-inflation, It is used within
#' \code{\link{sf_distrib}()}.
#' 
#' @param data_used The modified data frame to be used, typically \code{sdata}.
#' @param variable_vec A vector giving the names or numbers of the size or
#' fecundity variables to be tested in occasions \emph{t}+1 and \emph{t},
#' respectively.
#' @param term_used An integer denoting which variable is being worked with. If
#' \code{1}, then \code{sizea}; if \code{2}, then \code{sizeb}, if \code{3},
#' then \code{sizec}, and if \code{4}, then \code{fec}.
#' @param zi_state A logical value indicating whether to test for
#' zero-inflation.
#' @param size_a The vector of variable names or column numbers denoting
#' primary size (\code{sizea}). Used in fecundity assessment.
#' @param size_b The vector of variable names or column numbers denoting
#' secondary size (\code{sizeb}). Used in fecundity assessment.
#' @param size_c The vector of variable names or column numbers denoting
#' tertiary size (\code{sizec}). Used in fecundity assessment.
#' @param repst The vector of variable names or column numbers denoting
#' reproductive status (\code{repst}). Used in fecundity assessment.
#' @param fectime An integer denoting whether fecundity is assessed in time
#' *t*+1 (\code{3}) or time *t* (\code{2}). Used in fecundity assessment.
#' @param show_var A logical variable indicating whether to show the results of
#' tests for the particular variable in question.
#' 
#' @return This function produces text in the console giving the results of the
#' tests of overdispersion and zero inflation. No specific object is returned.
#' 
#' @section Notes:
#' This function will not test for overdispersion and zero inflation in
#' non-integer variables.
#' 
#' @keywords internal
#' @noRd
.knightswhosaynee <- function(data_used, variable_vec, term_used, zi_state,
  size_a, size_b, size_c, repst, fectime, show_var) {
  
  var_used <- FALSE
  
  full_term <- c("sizea", "sizeb", "sizec", "fec")
  full_inenglish <- c("Primary size", "Secondary size", "Tertiary size", "Fecundity")
  full_inenglish_small <- tolower(full_inenglish)
  
  if (!all(is.na(variable_vec))) {
    if (length(variable_vec) < 2) {
      
      stop(paste0("Function sf_distrib requires ", full_term[term_used], " in occasions t and t+1
        to function. Please designate these variables within a 2-element vector as input for",
        full_term[term_used]), call. = FALSE)
    } else if (length(variable_vec) > 2) {
      if (any(is.na(variable_vec[c(1,2)]))) {
        stop(paste0("NAs cannot be included within the first two elements of ", full_term[term_used], "."),
          call. = FALSE)
      }
      warning(paste0("Only the first two entries will be used in option ", full_term[term_used],
        " corresponding to occasions t+1 and t, respectively."), call. = FALSE)
    }
    if (all(is.numeric(variable_vec))) {
      if (any(!is.element(variable_vec, c(1:length(names(data_used)))))) {
        stop(paste0("Column numbers given in ", full_term[term_used]," are not recognized."), call. = FALSE)
      }
    } else if (all(is.character(variable_vec))) {
      if (any(!is.element(variable_vec, names(data_used)))) {
        stop(paste0("Column names given in ", full_term[term_used]," are not recognized."), call. = FALSE)
      }
    }
    
    var_used <- TRUE
    if (term_used < 4) {
      var3data <- data_used[, variable_vec[1]]
      var2data <- data_used[, variable_vec[2]]
    } else {
      if (!all(is.na(size_a)) & length(size_a) == 2) {
        var2 <- size_a
      } else if (!all(is.na(size_b)) & length(size_b) == 2) {
        var2 <- size_b
      } else if (!all(is.na(size_c)) & length(size_c) == 2) {
        var2 <- size_c
      }
      
      if (!all(is.na(repst))) {
        if (length(repst) < 2) {
          stop("Fecundity cannot be tested without repst holding the names of the reproductive status
            variables for times t+1 and t.", call. = FALSE)
        }
        if (fectime == 3) {
          repst_data <- data_used[which(data_used[,repst[1]] == 1),]
          var3data <- repst_data[, variable_vec[1]]
          var2data <- repst_data[, var2[1]]
        } else if (fectime == 2) {
          repst_data <- data_used[which(data_used[,repst[2]] == 1),]
          var3data <- repst_data[, variable_vec[2]]
          var2data <- repst_data[, var2[2]]
        } else {
          stop("Timing of fecundity not recognized.", call. = FALSE)
        }
        
      } else {
        if (fectime == 3) {
          var3data <- data_used[, variable_vec[1]]
          var2data <- data_used[, var2[1]]
        } else {
          var3data <- data_used[, variable_vec[2]]
          var2data <- data_used[, var2[2]]
        }
      }
    }
  
    #Here is the test of overdispersion
    jvmean <- mean(var3data)
    jvvar <- stats::var(var3data)
    
    v_pmodel <- stats::glm(var3data ~ var2data)
    v_qpmodel <- stats::glm(var3data ~ var2data)
    v_disp <- summary(v_qpmodel)$dispersion
    v_df <- summary(v_pmodel)$df.residual
    
    jvodchip <- stats::pchisq(v_disp * v_df, v_df, lower = FALSE)
    
    jvintcheck <- var3data%%1
    if (length(which(jvintcheck != 0)) > 0) {
      writeLines(paste0("Non-integer values detected, so will not test for overdispersion and zero-inflation in ",
          full_term[term_used]))
      show_var <- FALSE
      zi_state <- FALSE
    }
    
    if (show_var) {
      writeLines(paste0("Mean ", full_term[term_used]," is ", signif(jvmean, digits = 4)))
      writeLines(paste0("\nThe variance in ", full_term[term_used]," is ", signif(jvvar, digits = 4)))
      writeLines("\nThe probability of this dispersion level by chance assuming that")
      writeLines(paste0("the true mean ", full_term[term_used]," = variance in ", full_term[term_used], ","))
      writeLines(paste0("and an alternative hypothesis of overdispersion, is ", signif(jvodchip, digits = 4)))
      
      if (jvodchip <= 0.05 & jvvar > jvmean) {
        writeLines(paste0("\n", full_inenglish[term_used]," is significantly overdispersed."))
      } else if (jvodchip <= 0.05 & jvvar < jvmean) {
        writeLines(paste0("\n", full_inenglish[term_used]," is significantly underdispersed."))
      } else {
        writeLines(paste0("\nDispersion level in ", full_inenglish_small[term_used]," matches expectation."))
      }
    }
    
    #Here is the test of zero inflation
    if (zi_state) {
      
      v0est <- exp(-jvmean) #Estimated lambda
      v0n0 <- sum(var3data == 0) #Actual no of zeroes
      
      v0exp <- length(var3data) * v0est #Expected no of zeroes
      
      jvdbs <- (v0n0 - v0exp)^2 / (v0exp * (1 - v0est) - length(var3data) * jvmean * (v0est^2))
      jvzichip <- stats::pchisq(jvdbs, df = 1, lower.tail = FALSE)
      
      if (v0n0 < v0exp & jvzichip < 0.50) { #Correction for lower than expected numbers of 0s
        jvzichip <- 1 - jvzichip
      }
      
      if (show_var) {
        writeLines(paste0("\nMean lambda in ", full_term[term_used]," is ", signif(v0est, digits = 4)))
        writeLines(paste0("The actual number of 0s in ", full_term[term_used]," is ", v0n0))
        writeLines(paste0("The expected number of 0s in ", full_term[term_used]," under the null hypothesis is ", signif(v0exp, digits = 4)))
        writeLines(paste0("The probability of this deviation in 0s from expectation by chance is ", signif(jvzichip, digits = 4)))
        
        if (jvzichip <= 0.05 & v0n0 > v0exp) {
          writeLines(paste0("\n", full_inenglish[term_used]," is significantly zero-inflated.\n"))
        } else {
          writeLines(paste0("\n", full_inenglish[term_used]," is not significantly zero-inflated."))
          
          if (v0n0 == 0) {
            writeLines(paste0(full_inenglish[term_used]," does not appear to include 0s, suggesting
              that a zero-truncated distribution may be warranted."))
          }
          writeLines("\n")
        }
        
        writeLines("\n\n")
      }
    }
  }
}

#' Create a Data Frame of Elements Subject to Density Dependence
#' 
#' Function \code{density_input()} provides all necessary data to incorporate
#' density dependence into a \code{lefkoMat} object, a list of matrices, or a
#' single matrix. Three forms of density dependence are allowed, including the
#' Ricker function, the Beverton-Holt function, the Usher function, and the
#' logistic function. In each case, density must have an effect with at least a
#' one time-step delay (see Notes).
#'
#' @param mpm The \code{lefkoMat} object that will be subject to density
#' dependent projection.
#' @param stage3 A vector showing the name or number of the stage in occasion
#' \emph{t}+1 in the transitions to be affected by density. Abbreviations for
#' groups of stages are also usable (see Notes).
#' @param stage2 A vector showing the name or number of the stage in occasion
#' \emph{t} in the transition to be affected by density. Abbreviations for
#' groups of stages are also usable (see Notes).
#' @param stage1 A vector showing the name or number of the stage in occasion
#' \emph{t}-1 in the transition to be affected by density. Only needed if a
#' historical MPM is used. Abbreviations for groups of stages are also usable
#' (see Notes).
#' @param age2 A vector showing the age of the stage in occasion \emph{t} in the
#' transition to be affected by density. Only needed if an age-by-stage MPM is
#' used.
#' @param style A vector coding for the style of density dependence on each
#' transition subject to density dependence. Options include \code{1},
#' \code{ricker}, \code{ric}, or \code{r} for the Ricker function; \code{2},
#' \code{beverton}, \code{bev}, and \code{b} for the Beverton-Holt function;
#' \code{3}, \code{usher}, \code{ush}, and \code{u} for the Usher function; and
#' \code{4}, \code{logistic}, \code{log}, and \code{l} for the logistic
#' function. If only a single code is provided, then all noted transitions are
#' assumed to be subject to this style of density dependence. Defaults to
#' \code{ricker}.
#' @param time_delay A vector indicating the number of occasions back on which
#' density dependence operates. Defaults to 1, and may not equal any number less
#' than 1. If a single number is input, then all noted transitions are assumed
#' to be subject to this time delay.
#' @param alpha A vector indicating the numeric values to use as the
#' alpha term in the two parameter Ricker, Beverton-Holt, or Usher function, or
#' the value of the carrying capacity \emph{K} to use in the logistic equation
#' (see \code{Notes} section for more on this term). If a single number is
#' provided, then all noted transitions are assumed to be subject to this value
#' of alpha.
#' @param beta A vector indicating the numeric values to use as the beta term in
#' the two parameter Ricker, Beverton-Holt, or Usher function. Not used in the
#' logistic equation. If a single number is provided, then all noted transitions
#' are assumed to be subject to this value of beta.
#' @param type A vector denoting the kind of transition between occasions
#' \emph{t} and \emph{t}+1 to be replaced. This should be entered as \code{1},
#' \code{S}, or \code{s} for the replacement of a survival transition; or 
#' \code{2}, \code{F}, or \code{f} for the replacement of a fecundity
#' transition. If empty or not provided, then defaults to \code{1} for survival
#' transition.
#' @param type_t12 An optional vector denoting the kind of transition between
#' occasions \emph{t}-1 and \emph{t}. Only necessary if a historical MPM in
#' deVries format is desired. This should be entered as \code{1}, \code{S}, or
#' \code{s} for a survival transition; or \code{2}, \code{F}, or \code{f} for a
#' fecundity transitions. Defaults to \code{1} for survival transition, with
#' impacts only on the construction of deVries-format hMPMs.
#' 
#' @return A data frame of class \code{lefkoDens}. This object can be used as
#' input in function \code{\link{projection3}()}.
#' 
#' Variables in this object include the following:
#' \item{stage3}{Stage at occasion \emph{t}+1 in the transition to be replaced.}
#' \item{stage2}{Stage at occasion \emph{t} in the transition to be replaced.}
#' \item{stage1}{Stage at occasion \emph{t}-1 in the transition to be replaced,
#' if applicable.}
#' \item{age2}{Age at occasion \emph{t} in the transition to be replaced, if
#' applicable.}
#' \item{style}{Style of density dependence, coded as 1, 2, 3, or 4 for the
#' Ricker, Beverton-Holt, Usher, or logistic function, respectively.}
#' \item{time_delay}{The time delay on density dependence, in time steps.}
#' \item{alpha}{The value of alpha in the Ricker, Beverton-Holt, or Usher
#' function, or the value of carrying capacity, \emph{K}, in the logistic
#' function.}
#' \item{beta}{The value of beta in the Ricker, Beverton-Holt, or Usher
#' function.}
#' \item{type}{Designates whether the transition from occasion \emph{t} to
#' occasion \emph{t}+1 is a survival transition probability (1), or a fecundity
#' rate (2).}
#' \item{type_t12}{Designates whether the transition from occasion \emph{t}-1 to
#' occasion \emph{t} is a survival transition probability (1), a fecundity rate
#' (2).}
#' 
#' @section Notes:
#' The parameters \code{alpha} and \code{beta} are applied according to the
#' two-parameter Ricker function, the two-parameter Beverton-Holt function, the
#' two-parameter Usher function, or the one-parameter logistic function.
#' Although the default is that a 1 time step delay is assumed), greater time
#' delays can be set through the \code{time_delay} option.
#' 
#' Entries in \code{stage3}, \code{stage2}, and \code{stage1} can include
#' abbreviations for groups of stages. Use \code{rep} if all reproductive stages
#' are to be used, \code{nrep} if all mature but non-reproductive stages are to
#' be used, \code{mat} if all mature stages are to be used, \code{immat} if all
#' immature stages are to be used, \code{prop} if all propagule stages are to be
#' used, \code{npr} if all non-propagule stages are to be used, and leave empty
#' or use \code{all} if all stages in stageframe are to be used.
#' 
#' @seealso \code{\link{start_input}()}
#' @seealso \code{\link{projection3}()}
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
#'   fecacol = "Intactseed88", deadacol = "Dead1988",
#'   nonobsacol = "Dormant1988", stageassign = lathframe, stagesize = "sizea",
#'   censorcol = "Missing1988", censorkeep = NA, censor = TRUE)
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
#' ehrlen3mean <- lmean(ehrlen3)
#' 
#' e3d <- density_input(ehrlen3mean, stage3 = c("Sd", "Sdl"),
#'   stage2 = c("rep", "rep"), stage1 = c("all", "all"), style = 1,
#'   time_delay = 1, alpha = 1, beta = 0, type = c(2, 2), type_t12 = c(1, 1))
#' 
#' @export
density_input <- function(mpm, stage3, stage2, stage1 = NA, age2 = NA,
  style = 1, time_delay = 1, alpha = NA, beta = NA, type = NA, type_t12 = NA) {
  
  lstyle <- ltype <- ltype_t12 <- NULL
  
  if (!is.element("lefkoMat", class(mpm))) {
    stop("This function requires a lefkoMat object as input.", call. = FALSE)
  }
  
  stageframe <- mpm$ahstages
  
  if (all(is.na(mpm$hstages))) {
    historical <- FALSE
  } else {
    historical <- TRUE
  }
  if (all(is.na(mpm$agestages))) {
    agebystage <- FALSE
  } else {
    agebystage <- TRUE
  }
  
  if (agebystage & all(is.na(age2))) {
    stop("Density inputs for age-by-stage MPMs require the ages in time t of
      all transitions subject to density.", call. = FALSE)
  } else if (!agebystage & !all(is.na(age2))) {
    stop("Please use the age2 option only for age-by-stage MPMs.", call. = FALSE)
  }
  
  if (all(is.na(stageframe))) {
    stop("The input lefkoMat object does not appear to have a stageframe,
      which should be included as element ahstages.", call. = FALSE)
  }
  
  if (!is.element("stage", names(stageframe))) {
    stop("Stageframe appears to be modified. Please make sure that a stage
      column exists holding stage names.", call. = FALSE)
  }
  
  if (length(stage3) != length(stage2)) {
    stop("All transitions to modify require information at least for stage2 and
      stage3. These inputs must also be of equal length.", call. = FALSE)
  }
  
  if (historical & all(is.na(stage1))) {
    stop("Historical projection analysis requires that stage in time t-1 be
      designated for all density dependent transitions.", call. = FALSE)
  } else if (!historical & !all(is.na(stage1))) {
    stop("Ahistorical projection analysis cannot include designated stages in
      time t-1.", call. = FALSE)
  }
  if (historical & any(is.na(type_t12))) {
    stop("Historical projection analysis requires the kind of transition
      occurring between times t-1 and t to be described in option type_t12.",
      call. = FALSE)
  } else if (!historical & !all(is.na(type_t12))) {
    stop("Ahistorical projection analysis cannot include historical transitions.
      Please leave option type_t12 empty.", call. = FALSE)
  }
  
  full_length <- max(length(stage3), length(stage2), length(stage1),
    length(age2), length(style), length(time_delay), length(alpha),
    length(beta), length(type), length(type_t12))
  
  if (length(stage1) == 1 & full_length > 1) {
    stage1 <- rep(stage1, full_length)
  }
  
  if (length(age2) == 1 & full_length > 1) {
    age2 <- rep(age2, full_length)
  }
  
  if (length(style) == 1 & full_length > 1) {
    style <- rep(style, full_length)
  }
  
  if (length(time_delay) == 1 & full_length > 1) {
    time_delay <- rep(time_delay, full_length)
  }
  
  if (length(alpha) == 1 & full_length > 1) {
    alpha <- rep(alpha, full_length)
  }
  
  if (length(beta) == 1 & full_length > 1) {
    beta <- rep(beta, full_length)
  }
  
  if (length(type) == 1 & full_length > 1) {
    type <- rep(type, full_length)
  }
  
  if (length(type_t12) == 1 & full_length > 1) {
    type_t12 <- rep(type_t12, full_length)
  }
  
  if (length(stage3) != full_length) {
    stop("Please provide all input vectors in the same order.", call. = FALSE)
  }
  
  if (all(is.character(style))) {
    style <- tolower(style)
    
    ricker_style <- c("1", "ricker", "ricke", "rick", "ric", "ri", "r")
    beverton_style <- c("2", "beverton", "beverto", "bevert", "bever", "beve",
      "bev", "be", "b", "holt", "hol", "ho", "h")
    usher_style <- c("3", "usher", "ushe", "ush", "us", "u")
    logistic_style <- c("4", "logistic", "logisti", "logist", "logis", "logi",
      "log", "lo", "l")
    
    unknown_style <- which(!is.element(style, c(ricker_style, beverton_style, logistic_style)))
    if (length(unknown_style) > 0) {
      stop(paste0("Unknown style code used: ", style[unknown_style], ". Cannot
        process."), call. = FALSE)
    }
    
    lstyle <- rep(1, full_length)
    lstyle[which(is.element(style, beverton_style))] <- 2
    lstyle[which(is.element(style, usher_style))] <- 3
    lstyle[which(is.element(style, logistic_style))] <- 4
    
    style <- lstyle
  } else if (all(is.numeric(style))) {
    if (any(style > 4) | any(style < 1)) {
      stop("Unknown style code used. Please only use numbers 1, 2, 3, or 4.",
        call. = FALSE)
    }
  } else {
    stop("Input style codes do not conform to accepted inputs.", call. = FALSE)
  }
  
  if (!all(as.integer(time_delay) == time_delay)) {
    stop("Input for time_delay must be integer.", call. = FALSE)
  } else if (any(time_delay < 1)) {
    stop("Input for time_delay must be an integer no smaller than 1.",
      call. = FALSE)
  }
  
  if (any(!is.numeric(alpha)) & !all(is.na(alpha))) {
    stop("Option alpha must be either NA or a numeric value.", call. = FALSE)
  }
  if (any(!is.numeric(beta)) & !all(is.na(beta))) {
    stop("Option beta must be either NA or a numeric value.", call. = FALSE)
  }
  
  if (any(is.character(type))) {
    type <- tolower(type)
    
    if (!all(is.element(type, c(NA, "1", "2", "f", "s")))) {
      stop("Variable type must include only 1, 2, s, and f. All other entries
        are not allowed.", call. = FALSE)
    }
    ltype <- rep(1, full_length)
    ltype[which(is.element(type, c("2", "f")))] <- 2
    
    type <- ltype
  } else if (all(is.numeric(type))) {
    if (!all(is.element(type, c(NA, 1, 2)))) {
      stop("Variable type must include only 1, 2, s, and f. All other entries
        are not allowed.", call. = FALSE)
    }
  }
  if (any(is.character(type_t12))) {
    type_t12 <- tolower(type_t12)
    
    if (!all(is.element(type_t12, c(NA, "1", "2", "f", "s")))) {
      stop("Variable type_t12 must include only 1, 2, s, and f. All other
        entries are not allowed.", call. = FALSE)
    }
    ltype_t12 <- rep(1, full_length)
    ltype_t12[which(is.element(type_t12, c("2", "f")))] <- 2
    
    type_t12 <- ltype_t12
  } else if (all(is.numeric(type_t12))) {
    if (!all(is.element(type_t12, c(NA, 1, 2)))) {
      stop("Variable type_t12 must include only 1, 2, s, and f. All other
        entries are not allowed.", call. = FALSE)
    }
  }
  
  if (length(stage1) < full_length) {
    stop("Vector stage1 must be of the same length as stage 2 and stage3.",
      call. = FALSE)
  }
  if (length(age2) < full_length) {
    stop("Vector age2 must be of the same length as stage 2 and stage 3.",
      call. = FALSE)
  }
  if (length(alpha) < full_length) {
    stop("Vector alpha must be of the same length as stage 2 and stage3.",
      call. = FALSE)
  }
  if (length(beta) < full_length) {    
    stop("Vector beta must be of the same length as stage 2 and stage3.",
      call. = FALSE)
  }
  if (length(style) < full_length) {
    stop("Vector style must be of the same length as stage 2 and stage3.",
      call. = FALSE)
  }
  if (length(time_delay) < full_length) {
    stop("Vector time_delay must be of the same length as stage 2 and stage3.",
      call. = FALSE)
  }
  if (length(type) < full_length) {
    stop("Vector type must be of the same length as stage 2 and stage3.",
      call. = FALSE)
  }
  if (length(type_t12) < full_length) {
    stop("Vector type_t12 must be of the same length as stage 2 and stage3.",
      call. = FALSE)
  }
  
  all.stages.inp <- unique(c(stage3, stage2, stage1))
  
  mismatches <- !is.element(all.stages.inp, c(stageframe$stage, NA))
  
  if (length(which(mismatches)) > 0) {
    extrastuff <- tolower(all.stages.inp[which(mismatches)])
    
    unaccountedfor <- extrastuff[which(!is.element(extrastuff, 
          c("all", "rep", "nrep", "mat", "immat", "prop", "npr", "notalive")))]
    
    if (!all(is.element(extrastuff, c("all", "rep", "nrep", "mat", "immat",
          "prop", "npr", "notalive")))) {
      stop(paste("The following stage names input in supplemental() do not match
          the stageframe:", paste(unaccountedfor, collapse = ' ')),
        call. = FALSE)
    }
  }
  
  if (is.element("notalive", tolower(c(stage1, stage2, stage3)))) {
    stop("Stage NotAlive is not allowed.", call. = FALSE)
  }
  
  output_tab <- cbind.data.frame(stage3, stage2, stage1, age2, style,
    time_delay, alpha, beta, type, type_t12, stringsAsFactors = FALSE)
  
  output <- .density_reassess(stageframe, mpm$agestages, output_tab, historical,
    agebystage)
  
  if (!historical & !agebystage) {
    out_check <- unique(output[,c("stage3", "stage2")])
  } else if (historical & !agebystage) {
    out_check <- unique(output[,c("stage3", "stage2", "stage1")])
  } else if (!historical & agebystage) {
    out_check <- unique(output[,c("stage3", "stage2", "age2")])
  }
  if (dim(out_check)[1] < dim(output)[1]) {
    warning("Some transitions appear to be listed multiple times. This may cause
      errors in analysis.", call. = FALSE)
  }
  
  class(output) <- append(class(output), "lefkoDens")
  
  return(output)
}

#' Check and Reorganize Density Input Table Into Usable Format
#' 
#' Function \code{.density_reassess()} takes a density input table as supplied
#' by the \code{\link{density_input}()} function, and checks and rearranges it
#' into a single, complete density input table.
#' 
#' @param stageframe The correct stageframe, already modified by
#' \code{\link{.sf_reassess}()}.
#' @param agestages The agestages element from the used \code{lefkoMat} object.
#' @param dens_inp The density input data frame as is toward the end of
#' \code{\link{density_input}()}.
#' @param historical A logical value denoting whether MPM is historical.
#' @param agebystage A logical value denoting whether MPM is age-by-stage.
#' 
#' @return A corrected overwrite table, usable in MPM creation.
#' 
#' @keywords internal
#' @noRd
.density_reassess <- function(stageframe, agestages, dens_inp, historical, agebystage) {
  
  if (!all(is.na(dens_inp))) {
    shrubbery <- unique(dens_inp)
    
    if(any(duplicated(shrubbery[,1:4]))) {
      stop("Multiple entries with different values for the same stage transition
        are not allowed in the density input table. If modifying a historical
        table to perform an ahistorical analysis, then this may be due to
        different given rates of substitutions caused by dropping stage at
        occasion t-1. Please eliminate duplicate transitions.", call. = FALSE)
    }
    
  } else {
    stop("No recognized density data was input.", call. = FALSE)
  }
  
  if (!all(is.na(shrubbery))) {
    #First we make sure that the data is in the right format  
    shrubbery$stage3 <- as.character(shrubbery$stage3)
    shrubbery$stage2 <- as.character(shrubbery$stage2)
    shrubbery$stage1 <- as.character(shrubbery$stage1)
    shrubbery$age2 <- as.character(shrubbery$age2)
    
    #Stage at occasion t-1
    reassessed <- apply(as.matrix(c(1:dim(shrubbery)[1])), 1, function(X) {
      checkna2vec <- c(shrubbery[X, "stage3"], shrubbery[X, "stage2"])
      
      if (!all(!is.na(checkna2vec))) {
        stop("All entries for stage2 and stage3 in density input table must
          refer to possible life history stages and cannot include NAs.",
          call. = FALSE)
      }
      
      if (!is.na(shrubbery[X, "stage1"])) {
        if (is.element(shrubbery[X, "stage1"], stageframe$stage)) {
          return(shrubbery[X,])
        } else if (is.element(shrubbery[X, "stage1"], as.character(stageframe$stage_id))) {
          shrubbery.small <- cbind.data.frame(
            stage3 = shrubbery[X, "stage3"], stage2 = shrubbery[X, "stage1"],
            stage1 = stageframe$stage[which(as.character(stageframe$stage_id) == shrubbery[X, "stage1"])],
            age2 = shrubbery[X, "age2"], style = shrubbery[X, "style"],
            alpha = shrubbery[X, "alpha"], beta = shrubbery[X, "beta"],
            time_delay = shrubbery[X, "time_delay"], type = shrubbery[X, "type"],
            type_t12 = shrubbery[X, "type_t12"],
            stringsAsFactors = FALSE)
          
          return(shrubbery.small)
          
        } else if (shrubbery[X, "stage1"] == "rep") {
          shrubbery.small <- cbind.data.frame(stage3 = shrubbery[X, "stage3"], 
            stage2 = shrubbery[X, "stage2"], stage1 = stageframe$stage[which(stageframe$repstatus == 1)], 
            age2 = shrubbery[X, "age2"], style = shrubbery[X, "style"], 
            alpha = shrubbery[X, "alpha"], beta = shrubbery[X, "beta"],
            time_delay = shrubbery[X, "time_delay"], type = shrubbery[X, "type"],
            type_t12 = shrubbery[X, "type_t12"], stringsAsFactors = FALSE)
          
          return(shrubbery.small)
          
        } else if (shrubbery[X, "stage1"] == "nrep") {
          shrubbery.small <- cbind.data.frame(stage3 = shrubbery[X, "stage3"], 
            stage2 = shrubbery[X, "stage2"], 
            stage1 = stageframe$stage[intersect(which(stageframe$repstatus == 0),
                which(stageframe$matstatus == 1))], age2 = shrubbery[X, "age2"],
            style = shrubbery[X, "style"], alpha = shrubbery[X, "alpha"],
            beta = shrubbery[X, "beta"], time_delay = shrubbery[X, "time_delay"],
            type = shrubbery[X, "type"], type_t12 = shrubbery[X, "type_t12"],
            stringsAsFactors = FALSE)
          
          return(shrubbery.small)
          
        } else if (shrubbery[X, "stage1"] == "immat") {
          shrubbery.small <- cbind.data.frame(stage3 = shrubbery[X, "stage3"], 
            stage2 = shrubbery[X, "stage2"], stage1 = stageframe$stage[which(stageframe$immstatus == 1)], 
            age2 = shrubbery[X, "age2"], style = shrubbery[X, "style"],
            alpha = shrubbery[X, "alpha"], beta = shrubbery[X, "beta"],
            time_delay = shrubbery[X, "time_delay"], type = shrubbery[X, "type"],
            type_t12 = shrubbery[X, "type_t12"], stringsAsFactors = FALSE)
          
          return(shrubbery.small)
          
        } else if (shrubbery[X, "stage1"] == "mat") {
          shrubbery.small <- cbind.data.frame(stage3 = shrubbery[X, "stage3"], 
            stage2 = shrubbery[X, "stage2"], stage1 = stageframe$stage[which(stageframe$matstatus == 1)], 
            age2 = shrubbery[X, "age2"], style = shrubbery[X, "style"],
            alpha = shrubbery[X, "alpha"], beta = shrubbery[X, "beta"],
            time_delay = shrubbery[X, "time_delay"], type = shrubbery[X, "type"],
            type_t12 = shrubbery[X, "type_t12"], stringsAsFactors = FALSE)
          
          return(shrubbery.small)
          
        } else if (shrubbery[X, "stage1"] == "prop") {
          shrubbery.small <- cbind.data.frame(stage3 = shrubbery[X, "stage3"], 
            stage2 = shrubbery[X, "stage2"], stage1 = stageframe$stage[which(stageframe$propstatus == 1)], 
            age2 = shrubbery[X, "age2"], style = shrubbery[X, "style"],
            alpha = shrubbery[X, "alpha"], beta = shrubbery[X, "beta"],
            time_delay = shrubbery[X, "time_delay"], type = shrubbery[X, "type"],
            type_t12 = shrubbery[X, "type_t12"], stringsAsFactors = FALSE)
          
          return(shrubbery.small)
          
        } else if (shrubbery[X, "stage1"] == "npr") {
          shrubbery.small <- cbind.data.frame(stage3 = shrubbery[X, "stage3"], 
            stage2 = shrubbery[X, "stage2"], stage1 = stageframe$stage[which(stageframe$propstatus == 0)], 
            age2 = shrubbery[X, "age2"], style = shrubbery[X, "style"],
            alpha = shrubbery[X, "alpha"], beta = shrubbery[X, "beta"],
            time_delay = shrubbery[X, "time_delay"], type = shrubbery[X, "type"],
            type_t12 = shrubbery[X, "type_t12"], stringsAsFactors = FALSE)
          
          return(shrubbery.small)
          
        } else if (shrubbery[X, "stage1"] == "all") {
          shrubbery.small <- cbind.data.frame(stage3 = shrubbery[X, "stage3"], 
            stage2 = shrubbery[X, "stage2"], stage1 = stageframe$stage, 
            age2 = shrubbery[X, "age2"], style = shrubbery[X, "style"],
            alpha = shrubbery[X, "alpha"], beta = shrubbery[X, "beta"],
            time_delay = shrubbery[X, "time_delay"], type = shrubbery[X, "type"],
            type_t12 = shrubbery[X, "type_t12"], stringsAsFactors = FALSE)
          
          return(shrubbery.small)
        }
      } else {
        return(shrubbery[X,])
      }
    })
    shrubbery <- do.call(rbind.data.frame, reassessed)
    
    #Stage at occasion t
    reassessed <- apply(as.matrix(c(1:dim(shrubbery)[1])), 1, function(X) {
      if (!is.na(shrubbery[X, "stage2"])) {
        if (is.element(shrubbery[X, "stage2"], stageframe$stage)) {
          return(shrubbery[X,])
        } else if (is.element(shrubbery[X, "stage2"], as.character(stageframe$stage_id))) {
          shrubbery.small <- cbind.data.frame(
            stage3 = shrubbery[X, "stage3"], 
            stage2 = stageframe$stage[which(as.character(stageframe$stage_id) == shrubbery[X, "stage2"])],
            stage1 = shrubbery[X, "stage1"], age2 = shrubbery[X, "age2"],
            style = shrubbery[X, "style"], alpha = shrubbery[X, "alpha"],
            beta = shrubbery[X, "beta"], time_delay = shrubbery[X, "time_delay"],
            type = shrubbery[X, "type"], type_t12 = shrubbery[X, "type_t12"],
            stringsAsFactors = FALSE)
          
          return(shrubbery.small)
          
        } else if (shrubbery[X, "stage2"] == "rep") {
          shrubbery.small <- cbind.data.frame(stage3 = shrubbery[X, "stage3"], 
            stage2 = stageframe$stage[which(stageframe$repstatus == 1)],
            stage1 = shrubbery[X, "stage1"], age2 = shrubbery[X, "age2"],
            style = shrubbery[X, "style"], alpha = shrubbery[X, "alpha"],
            beta = shrubbery[X, "beta"], time_delay = shrubbery[X, "time_delay"],
            type = shrubbery[X, "type"], type_t12 = shrubbery[X, "type_t12"],
            stringsAsFactors = FALSE)
          
          return(shrubbery.small)
          
        } else if (shrubbery[X, "stage2"] == "nrep") {
          shrubbery.small <- cbind.data.frame(stage3 = shrubbery[X, "stage3"], 
            stage2 = stageframe$stage[intersect(which(stageframe$repstatus == 0),
                which(stageframe$matstatus == 1))], stage1 = shrubbery[X, "stage1"], 
            age2 = shrubbery[X, "age2"], style = shrubbery[X, "style"],
            alpha = shrubbery[X, "alpha"], beta = shrubbery[X, "beta"],
            time_delay = shrubbery[X, "time_delay"], type = shrubbery[X, "type"],
            type_t12 = shrubbery[X, "type_t12"], stringsAsFactors = FALSE)
          
          return(shrubbery.small)
          
        } else if (shrubbery[X, "stage2"] == "immat") {
          shrubbery.small <- cbind.data.frame(stage3 = shrubbery[X, "stage3"], 
            stage2 = stageframe$stage[which(stageframe$immstatus == 1)],
            stage1 = shrubbery[X, "stage1"], age2 = shrubbery[X, "age2"],
            style = shrubbery[X, "style"], alpha = shrubbery[X, "alpha"],
            beta = shrubbery[X, "beta"], time_delay = shrubbery[X, "time_delay"],
            type = shrubbery[X, "type"], type_t12 = shrubbery[X, "type_t12"],
            stringsAsFactors = FALSE)
          
          return(shrubbery.small)
          
        } else if (shrubbery[X, "stage2"] == "mat") {
          shrubbery.small <- cbind.data.frame(stage3 = shrubbery[X, "stage3"], 
            stage2 = stageframe$stage[which(stageframe$matstatus == 1)],
            stage1 = shrubbery[X, "stage1"], age2 = shrubbery[X, "age2"],
            style = shrubbery[X, "style"], alpha = shrubbery[X, "alpha"],
            beta = shrubbery[X, "beta"], time_delay = shrubbery[X, "time_delay"],
            type = shrubbery[X, "type"], type_t12 = shrubbery[X, "type_t12"],
            stringsAsFactors = FALSE)
          
          return(shrubbery.small)
          
        } else if (shrubbery[X, "stage2"] == "prop") {
          shrubbery.small <- cbind.data.frame(stage3 = shrubbery[X, "stage3"], 
            stage2 = stageframe$stage[which(stageframe$propstatus == 1)],
            stage1 = shrubbery[X, "stage1"], age2 = shrubbery[X, "age2"],
            style = shrubbery[X, "style"], alpha = shrubbery[X, "alpha"],
            beta = shrubbery[X, "beta"], time_delay = shrubbery[X, "time_delay"],
            type = shrubbery[X, "type"], type_t12 = shrubbery[X, "type_t12"],
            stringsAsFactors = FALSE)
          
          return(shrubbery.small)
          
        } else if (shrubbery[X, "stage2"] == "npr") {
          shrubbery.small <- cbind.data.frame(stage3 = shrubbery[X, "stage3"], 
            stage2 = stageframe$stage[which(stageframe$propstatus == 0)],
            stage1 = shrubbery[X, "stage1"], age2 = shrubbery[X, "age2"],
            style = shrubbery[X, "style"], alpha = shrubbery[X, "alpha"],
            beta = shrubbery[X, "beta"], time_delay = shrubbery[X, "time_delay"],
            type = shrubbery[X, "type"], type_t12 = shrubbery[X, "type_t12"],
            stringsAsFactors = FALSE)
          
          return(shrubbery.small)
          
        } else if (shrubbery[X, "stage2"] == "all") {
          shrubbery.small <- cbind.data.frame(stage3 = shrubbery[X, "stage3"], 
            stage2 = stageframe$stage, stage1 = shrubbery[X, "stage1"], 
            age2 = shrubbery[X, "age2"], style = shrubbery[X, "style"],
            alpha = shrubbery[X, "alpha"], beta = shrubbery[X, "beta"],
            time_delay = shrubbery[X, "time_delay"], type = shrubbery[X, "type"],
            type_t12 = shrubbery[X, "type_t12"], stringsAsFactors = FALSE)
          
          return(shrubbery.small)
        }
      } else {
        return(shrubbery[X,])
      }
    })
    shrubbery <- do.call(rbind.data.frame, reassessed)
    
    #Stage at occasion t+1
    reassessed <- apply(as.matrix(c(1:dim(shrubbery)[1])), 1, function(X) {
      if (!is.na(shrubbery[X, "stage3"])) {
        if (is.element(shrubbery[X, "stage3"], stageframe$stage)) {
          return(shrubbery[X,])
        } else if (is.element(shrubbery[X, "stage3"], as.character(stageframe$stage_id))) {
          shrubbery.small <- cbind.data.frame(
            stage3 = stageframe$stage[which(as.character(stageframe$stage_id) == shrubbery[X, "stage3"])], 
            stage2 = shrubbery[X, "stage2"], stage1 = shrubbery[X, "stage1"], 
            age2 = shrubbery[X, "age2"], style = shrubbery[X, "style"],
            alpha = shrubbery[X, "alpha"], beta = shrubbery[X, "beta"],
            time_delay = shrubbery[X, "time_delay"], type = shrubbery[X, "type"],
            type_t12 = shrubbery[X, "type_t12"], stringsAsFactors = FALSE)
          
          return(shrubbery.small)
          
        } else if (shrubbery[X, "stage3"] == "rep") {
          shrubbery.small <- cbind.data.frame(stage3 = stageframe$stage[which(stageframe$repstatus == 1)], 
            stage2 = shrubbery[X, "stage2"], stage1 = shrubbery[X, "stage1"], 
            age2 = shrubbery[X, "age2"], style = shrubbery[X, "style"],
            alpha = shrubbery[X, "alpha"], beta = shrubbery[X, "beta"],
            time_delay = shrubbery[X, "time_delay"], type = shrubbery[X, "type"],
            type_t12 = shrubbery[X, "type_t12"], stringsAsFactors = FALSE)
          
          return(shrubbery.small)
          
        } else if (shrubbery[X, "stage3"] == "nrep") {
          shrubbery.small <- cbind.data.frame(stage3 = stageframe$stage[intersect(which(stageframe$repstatus == 0),
                which(stageframe$matstatus == 1))], 
            stage2 = shrubbery[X, "stage2"], stage1 = shrubbery[X, "stage1"], 
            age2 = shrubbery[X, "age2"], style = shrubbery[X, "style"],
            alpha = shrubbery[X, "alpha"], beta = shrubbery[X, "beta"],
            time_delay = shrubbery[X, "time_delay"], type = shrubbery[X, "type"],
            type_t12 = shrubbery[X, "type_t12"], stringsAsFactors = FALSE)
          
          return(shrubbery.small)
          
        } else if (shrubbery[X, "stage3"] == "immat") {
          shrubbery.small <- cbind.data.frame(stage3 = stageframe$stage[which(stageframe$immstatus == 1)], 
            stage2 = shrubbery[X, "stage2"], stage1 = shrubbery[X, "stage1"], 
            age2 = shrubbery[X, "age2"], style = shrubbery[X, "style"],
            alpha = shrubbery[X, "alpha"], beta = shrubbery[X, "beta"],
            time_delay = shrubbery[X, "time_delay"], type = shrubbery[X, "type"],
            type_t12 = shrubbery[X, "type_t12"], stringsAsFactors = FALSE)
          
          return(shrubbery.small)
          
        } else if (shrubbery[X, "stage3"] == "mat") {
          shrubbery.small <- cbind.data.frame(stage3 = stageframe$stage[which(stageframe$matstatus == 1)], 
            stage2 = shrubbery[X, "stage2"], stage1 = shrubbery[X, "stage1"], 
            age2 = shrubbery[X, "age2"], style = shrubbery[X, "style"],
            alpha = shrubbery[X, "alpha"], beta = shrubbery[X, "beta"],
            time_delay = shrubbery[X, "time_delay"], type = shrubbery[X, "type"],
            type_t12 = shrubbery[X, "type_t12"], stringsAsFactors = FALSE)
          
          return(shrubbery.small)
          
        } else if (shrubbery[X, "stage3"] == "prop") {
          shrubbery.small <- cbind.data.frame(stage3 = stageframe$stage[which(stageframe$propstatus == 1)], 
            stage2 = shrubbery[X, "stage2"], stage1 = shrubbery[X, "stage1"], 
            age2 = shrubbery[X, "age2"], style = shrubbery[X, "style"],
            alpha = shrubbery[X, "alpha"], beta = shrubbery[X, "beta"],
            time_delay = shrubbery[X, "time_delay"], type = shrubbery[X, "type"],
            type_t12 = shrubbery[X, "type_t12"], stringsAsFactors = FALSE)
          
          return(shrubbery.small)
          
        } else if (shrubbery[X, "stage3"] == "npr") {
          shrubbery.small <- cbind.data.frame(stage3 = stageframe$stage[which(stageframe$propstatus == 0)], 
            stage2 = shrubbery[X, "stage2"], stage1 = shrubbery[X, "stage1"], 
            age2 = shrubbery[X, "age2"], style = shrubbery[X, "style"],
            alpha = shrubbery[X, "alpha"], beta = shrubbery[X, "beta"],
            time_delay = shrubbery[X, "time_delay"], type = shrubbery[X, "type"],
            type_t12 = shrubbery[X, "type_t12"], stringsAsFactors = FALSE)
          
          return(shrubbery.small)
          
        } else if (shrubbery[X, "stage3"] == "all") {
          shrubbery.small <- cbind.data.frame(stage3 = stageframe$stage, 
            stage2 = shrubbery[X, "stage2"], stage1 = shrubbery[X, "stage1"], 
            age2 = shrubbery[X, "age2"], style = shrubbery[X, "style"],
            alpha = shrubbery[X, "alpha"], beta = shrubbery[X, "beta"],
            time_delay = shrubbery[X, "time_delay"], type = shrubbery[X, "type"],
            type_t12 = shrubbery[X, "type_t12"], stringsAsFactors = FALSE)
          
          return(shrubbery.small)
          
        }
      } else {
        return(shrubbery[X,])
      }
    })
    shrubbery <- do.call(rbind.data.frame, reassessed)
    
    if (agebystage) {
      #Age at occasion t
      reassessed <- apply(as.matrix(c(1:dim(shrubbery)[1])), 1, function(X) {
        if (!is.na(shrubbery[X, "age2"])) {
          if (is.element(shrubbery[X, "age2"], agestages$age)) {
            return(shrubbery[X,])
          } else if (shrubbery[X, "age2"] == "all") {
            identified_ages <- agestages$age[which(agestages$stage == shrubbery[X, "stage2"])]
            
            if (length(identified_ages > 0)) {
              shrubbery.small <- cbind.data.frame(stage3 = shrubbery[X, "stage3"], 
                stage2 = shrubbery[X, "stage2"], stage1 = shrubbery[X, "stage1"], 
                age2 = identified_ages, style = shrubbery[X, "style"],
                alpha = shrubbery[X, "alpha"], beta = shrubbery[X, "beta"],
                time_delay = shrubbery[X, "time_delay"], type = shrubbery[X, "type"],
                type_t12 = shrubbery[X, "type_t12"], stringsAsFactors = FALSE)
            } else {
              shrubbery.small <- shrubbery
            }
            return(shrubbery.small)
          }
        } else {
          return(shrubbery[X,])
        }
      })
      shrubbery <- do.call(rbind.data.frame, reassessed)
    }
    
    #Now a bit of a check to remove entries that are not allowed
    stufftoremove <- unique(c(which(shrubbery$stage1 == "Dead"),
      which(shrubbery$stage2 == "Dead"), which(shrubbery$stage3 == "Dead")))
    
    if (length(stufftoremove) > 0) {
      if (stufftoremove[1] > 0) {
        shrubbery <- shrubbery[-stufftoremove,]
      }
    }
  }
  
  return(shrubbery)
}

#' Create a Starting Vector for Population Projection
#' 
#' Function \code{start_input()} creates a data frame summarizing the non-zero
#' elements of the start vector for use in population projection analysis via
#' function \code{\link{projection3}()}.
#'
#' @param mpm The lefkoMat object to be used in projection analysis.
#' @param stage2 A vector showing the name or number of a stage in occasion
#' \emph{t} that should be set to a positive number of individuals in the start
#' vector. Abbreviations for groups of stages are also usable (see Notes).
#' This input is required and has no default input.
#' @param stage1 A vector showing the name or number of a stage in occasion
#' \emph{t}-1 that should be set to a positive number of individuals in the
#' start vector. Abbreviations for groups of stages are also usable (see Notes).
#' This is only used for historical MPMs, since the rows of hMPMs correspond to
#' stage-pairs in times \emph{t} and \emph{t}-1 together. Only required for
#' historical MPMs, and will result in errors if otherwise used.
#' @param age2 A vector showing the age of each respective stage in occasion
#' \emph{t} that should be set to a positie number of individuals in the start
#' vector. Only used for age-by-stage MPMs. Defaults to NA.
#' @param value A vector showing the values, in order, of the number of
#' individuals set for the stage or stage-pair in question. Defaults to 1.
#' 
#' @return A list of class \code{lefkoStart}, with 4 objects, which can be used
#' as input in function \code{\link{projection3}()}. The last three include the
#' \code{ahstages}, \code{hstages}, and \code{agestages} objects from the
#' \code{lefkoMat} object supplied in \code{mpm}. The first element in the list
#' is a data frame with the following variables:
#' 
#' \item{stage2}{Stage at occasion \emph{t}.}
#' \item{stage_id_2}{The stage number associated with \code{stage2}.}
#' \item{stage1}{Stage at occasion \emph{t}-1, if historical. Otherwise NA.}
#' \item{stage_id_1}{The stage number associated with \code{stage1}.}
#' \item{age2}{The age of individuals in \code{stage2} and, if applicable,
#' \code{stage1}. Only used in age-by-stage MPMs.}
#' \item{row_num}{A number indicating the respective starting vector element.}
#' \item{value}{Number of individuals in corresponding stage or stage-pair.}
#' 
#' @section Notes:
#' Entries in \code{stage2}, and \code{stage1} can include abbreviations for
#' groups of stages. Use \code{rep} if all reproductive stages are to be used,
#' \code{nrep} if all mature but non-reproductive stages are to be used,
#' \code{mat} if all mature stages are to be used, \code{immat} if all immature
#' stages are to be used, \code{prop} if all propagule stages are to be used,
#' \code{npr} if all non-propagule stages are to be used, and leave empty or use
#' \code{all} if all stages in stageframe are to be used.
#' 
#' @seealso \code{\link{density_input}()}
#' @seealso \code{\link{projection3}()}
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
#'   fecacol = "Intactseed88", deadacol = "Dead1988",
#'   nonobsacol = "Dormant1988", stageassign = lathframe, stagesize = "sizea",
#'   censorcol = "Missing1988", censorkeep = NA, censor = TRUE)
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
#' ehrlen3mean <- lmean(ehrlen3)
#' 
#' e3m_sv <- start_input(ehrlen3mean, stage2 = "Sd", stage1 = "Sd", value = 1000)
#' 
#' @export
start_input <- function(mpm, stage2, stage1 = NA, age2 = NA, value = 1) {
  
  mpmrows <- stage2_id <- stage1_id <- start_vec <- NULL
  
  if (all(class(mpm) != "lefkoMat")) {
    stop("A regular lefkoMat object is required as input.", call. = FALSE)
  }
  
  if (!is.element("stage", names(mpm$ahstages))) {
    stop("Stageframe appears to be modified. Please make sure that a stage
      column exists holding stage names.", call. = FALSE)
  }
  
  if (all(is.na(mpm$hstages)) | all(is.null(mpm$hstages))) {
    historical <- FALSE
  } else {
    historical <- TRUE
  }
  if (all(is.na(mpm$agestages)) | all(is.null(mpm$agestages))) {
    agebystage <- FALSE
  } else {
    agebystage <- TRUE
  }
  
  if (historical & all(is.na(stage1))) {
    stop("Historical projection analysis requires that stage in time t-1 be
      designated for all stage pairs.", call. = FALSE)
  } else if (!historical & !all(is.na(stage1))) {
    stop("Ahistorical projection analysis cannot include designated stages in
      time t-1.", call. = FALSE)
  }
  
  full_length <- max(length(stage2), length(stage1), length(age2), length(value))
  
  if (length(value) == 1 & full_length > 1) {
    value <- rep(value, full_length)
  }
  
  if (length(stage1) == 1 & full_length > 1) {
    stage1 <- rep(stage1, full_length)
  }
  
  if (length(age2) == 1 & full_length > 1) {
    age2 <- rep(age2, full_length)
  }
  
  if (length(stage2) != full_length) {
    stop("Option stage2 is required for all stages or stage-pairs to set to
      non-zero values.", call. = FALSE)
  }
  
  if (all(is.character(stage2))) {
    
    unknown_stage2 <- which(!is.element(tolower(stage2), c(tolower(mpm$ahstages$stage),
        c("all", "rep", "nrep", "mat", "immat", "prop", "npr"))))
    if (length(unknown_stage2) > 0) {
      stop(paste0("Unknown stage designations used in stage2: ",
        stage2[unknown_stage2]), call. = FALSE)
    }
    
    reassessed <- apply(as.matrix(c(1:length(stage2))), 1, function(X) {
      if (!is.na(stage2[X])) {
        if (is.element(stage2[X], mpm$ahstages$stage)) {
          shrubbery.small <- cbind.data.frame(stage2 = stage2[X], stage1 = stage1[X],
            age2 = age2[X], value = value[X], stringsAsFactors = FALSE)
          return(shrubbery.small)
        } else if (is.element(stage2[X], as.character(mpm$ahstages$stage_id))) {
          shrubbery.small <- cbind.data.frame(stage2 = as.numeric(stage2[X]),
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "rep") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$repstatus == 1)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "all") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage,
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "nrep") {
          shrubbery.small <- cbind.data.frame(
            stage2 = mpm$ahstages$stage[intersect(which(mpm$ahstages$repstatus == 0),
                which(mpm$ahstages$matstatus == 1))],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) =="mat") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$matstatus == 1)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "immat") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$immstatus == 1)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "prop") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$propstatus == 1)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        } else if (tolower(stage2[X]) == "npr") {
          shrubbery.small <- cbind.data.frame(stage2 = mpm$ahstages$stage[which(mpm$ahstages$propstatus == 0)],
            stage1 = stage1[X], age2 = age2[X], value = value[X],
            stringsAsFactors = FALSE)
            
          return(shrubbery.small)
        }
      }
    })
    
    shrubbery <- do.call(rbind.data.frame, reassessed)
  } else if (all(is.numeric(stage2)) & !any(is.na(stage2))) {
    stage2_id <- stage2
    
    if (any(stage2_id > max(mpm$ahstages$stage_id)) | any(stage2_id < min(mpm$ahstages$stage_id))) {
      stop("Unknown stage2 codes used.", call. = FALSE)
    }
    
    stage2 <- apply(as.matrix(stage2_id), 1, function(X) {
      return(mpm$ahstages$stage[X])
    })
    
    shrubbery <- cbind.data.frame(stage2 = stage2, stage1 = stage1, age2 = age2,
      value = value, stringsAsFactors = FALSE)
  } else {
    stop("Input stage2 codes do not conform to accepted inputs.", call. = FALSE)
  }
  
  if (historical) {
    if (all(is.character(shrubbery$stage1)) & !all(is.na(shrubbery$stage1))) {
      
      unknown_stage1 <- which(!is.element(tolower(stage1), c(tolower(mpm$ahstages$stage),
          c("all", "rep", "nrep", "mat", "immat", "prop", "npr", "almostborn"))))
      if (length(unknown_stage1) > 0) {
      stop(paste0("Unknown stage designations used in stage1: ",
        stage1[unknown_stage1]), call. = FALSE)
      }
      
      reassessed <- apply(as.matrix(c(1:length(shrubbery$stage2))), 1, function(X) {
        if (!is.na(shrubbery$stage1[X])) {
          if (is.element(shrubbery$stage1[X], mpm$ahstages$stage)) {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = shrubbery$stage1[X], age2 = shrubbery$age2[X],
              value = shrubbery$value[X], stringsAsFactors = FALSE)
            return(shrubbery.small)
          } else if (is.element(stage1[X], as.character(mpm$ahstages$stage_id))) {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = as.numeric(shrubbery$stage1[X]), age2 = age2[X],
              value = value[X], stringsAsFactors = FALSE)
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "rep") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$repstatus == 1)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "all") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage, age2 = shrubbery$age2[X],
              value = shrubbery$value[X], stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "nrep") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[intersect(which(mpm$ahstages$repstatus == 0),
                  which(mpm$ahstages$matstatus == 1))],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) =="mat") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$matstatus == 1)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "immat") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$immstatus == 1)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "prop") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$propstatus == 1)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          } else if (tolower(shrubbery$stage1[X]) == "npr") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = mpm$ahstages$stage[which(mpm$ahstages$propstatus == 0)],
              age2 = shrubbery$age2[X], value = shrubbery$value[X],
              stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          }
        }
      })
      
      shrubbery <- do.call(rbind.data.frame, reassessed)
      
    } else if (all(is.numeric(shrubbery$stage1)) & !any(is.na(shrubbery$stage1))) {
      stage1_id <- shrubbery$stage1
      
      if (any(stage1_id > max(mpm$ahstages$stage_id)) | any(stage1_id < min(mpm$ahstages$stage_id))) {
        stop("Unknown stage1 codes used.", call. = FALSE)
      }
      
      stage1 <- apply(as.matrix(stage1_id), 1, function(X) {
        return(mpm$ahstages$stage[X])
      })
      shrubbery <- cbind.data.frame(stage2 = shrubbery$stage2, stage1 = stage1,
        age2 = shrubbery$age2, value = shrubbery$value, stringsAsFactors = FALSE)
    } else {
      stop("Input stage1 codes do not conform to accepted inputs.", call. = FALSE)
    }
  }
  
  if (agebystage) {
    if (all(is.character(shrubbery$age2)) & !all(is.na(shrubbery$age2))) {
      
      unknown_age2 <- which(!is.element(tolower(age2), c(tolower(mpm$agestages$age),
          "all")))
      if (length(unknown_age2) > 0) {
        stop(paste0("Unknown age designations used in age2: ",
          stage1[unknown_age2]), call. = FALSE)
      }
      
      reassessed <- apply(as.matrix(c(1:length(shrubbery$stage2))), 1, function(X) {
        if (!is.na(shrubbery$age2[X])) {
          common_ages <- unique(mpm$agestages$age[which(mpm$agestages$stage == shrubbery$stage2[X])])
          
          if (is.element(shrubbery$age2[X], as.character(common_ages))) {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = shrubbery$stage1[X], age2 = as.numeric(shrubbery$age2[X]),
              value = shrubbery$value[X], stringsAsFactors = FALSE)
            return(shrubbery.small)
          } else if (tolower(shrubbery$age2[X]) == "all") {
            shrubbery.small <- cbind.data.frame(stage2 = shrubbery$stage2[X],
              stage1 = shrubbery$stage1[X], age2 = common_ages,
              value = shrubbery$value[X], stringsAsFactors = FALSE)
              
            return(shrubbery.small)
          }
        }
      }) 
      
      shrubbery <- do.call(rbind.data.frame, reassessed)
    } else if (all(is.numeric(shrubbery$age2)) & !any(is.na(shrubbery$age2))) {
      if (any(shrubbery$age2 > max(mpm$agestages$age)) | any(shrubbery$age2 < min(mpm$agestages$age))) {
        stop("Unknown age2 values used.", call. = FALSE)
      }
      
      shrubbery <- cbind.data.frame(stage2 = shrubbery$stage2, stage1 = shrubbery$stage1,
        age2 = shrubbery$age2, value = shrubbery$value, stringsAsFactors = FALSE)
    } else {
      stop("Input stage1 codes do not conform to accepted inputs.", call. = FALSE)
    }
  }
  
  if (!all(is.numeric(value))) {
    stop("Object value must be composed only of valid numbers.", call. = FALSE)
  }
  
  shrubbery$stage2_id <- apply(as.matrix(shrubbery$stage2), 1, function(X) {
    return(mpm$ahstages$stage_id[which(mpm$ahstages$stage == X)])
  })
  shrubbery$stage1_id <- apply(as.matrix(shrubbery$stage1), 1, function(X) {
    possible_option <- mpm$ahstages$stage_id[which(mpm$ahstages$stage == X)]
    if (length(possible_option) > 0) return(possible_option) else return(NA)
  })
  
  full_length <- dim(shrubbery)[1]
  
  if (!historical & !agebystage) {
    if (dim(mpm$A[[1]])[1] != dim(mpm$ahstages)[1]) {
      stop("This ahistorical mpm includes matrices with dimensions that do not
        match expectation.", call. = FALSE)
    }
    
    start_vec <- shrubbery$stage2_id
    
  } else if (agebystage & !historical) {
    if (dim(mpm$A[[1]])[1] != dim(mpm$agestages)[1]) {
      stop("This age-by-stage mpm includes matrices with dimensions that do not
        match expectation.", call. = FALSE)
    }
    
    if (any(is.na(shrubbery$age2)) | any(!is.numeric(shrubbery$age2))) {
      stop("Option age2 must include only numbers for age-by-stage MPMs.",
        call. = FALSE)
    }
    
    if (any(shrubbery$age2 < min(mpm$agestages$age)) | any(shrubbery$age2 > max(mpm$agestages$age))) {
      stop("Option age2 can only take ages shown in element agestages within the input MPM.",
        call. = FALSE)
    }
    
    start_vec <- apply(as.matrix(c(1:full_length)), 1, function(X) {
      vec2 <- which(mpm$agestages$stage_id == shrubbery$stage2_id[X])
      vec1 <- which(mpm$agestages$age == shrubbery$age2[X])
      
      return(intersect(vec2, vec1)[1])
    })
    
  } else if (historical & !agebystage) {
    if (dim(mpm$A[[1]])[1] != dim(mpm$hstages)[1]) {
      stop("This historical mpm includes matrices with dimensions that do not
        match expectation.", call. = FALSE)
    }
    
    start_vec <- apply(as.matrix(c(1:full_length)), 1, function(X) {
      vec2 <- which(mpm$hstages$stage_id_2 == shrubbery$stage2_id[X])
      vec1 <- which(mpm$hstages$stage_id_1 == shrubbery$stage1_id[X])
      
      return(intersect(vec2, vec1)[1])
    })
    
  } else {
    stop("Format of mpm not recognized.", call. = FALSE)
  }
  
  output_tab <- cbind.data.frame(shrubbery$stage2, shrubbery$stage2_id,
    shrubbery$stage1, shrubbery$stage1_id, shrubbery$age2, start_vec,
    shrubbery$value, stringsAsFactors = FALSE)
  
  names(output_tab) <- c("stage2", "stage_id_2", "stage1", "stage_id_1", "age2",
    "row_num", "value")
  
  if (!historical & !agebystage) {
    out_check <- unique(output_tab[,c("stage2", "stage_id_2")])
  } else if (historical & !agebystage) {
    out_check <- unique(output_tab[,c("stage2", "stage_id_2", "stage1", "stage_id_1")])
  } else if (!historical & agebystage) {
    out_check <- unique(output_tab[,c("stage2", "stage_id_2", "age2")])
  }
  if (dim(out_check)[1] < dim(output_tab)[1]) {
    warning("Some stages, stage-pairs, or age-stages appear to be listed
      multiple times. This may cause errors in analysis.", call. = FALSE)
  }

  class(output_tab) <- append(class(output_tab), "lefkoSV")
  
  return(output_tab)
}

#' Calculate Actual Stage or Stage-Pair Distributions
#' 
#' Function \code{actualstage3()} shows the frequencies and proportions of
#' each stage or stage pair in each year.
#' 
#' @param data A demographic dataset in hfv format.
#' @param historical A logical value indicating whether the stage structure
#' should be ahistorical (\code{FALSE}) or historical (\code{TRUE}). Defaults to
#' \code{FALSE}.
#' @param year2 A string value indicating the name of the variable coding for
#' monitoring occasion at time \emph{t}.
#' @param indices A vector of three strings, indicating the stage indices for
#' times \emph{t}+1, \emph{t}, and \emph{t}-1, respectively, in \code{data}.
#' Defaults to \code{c("stage3index", "stage2index", "stage1index")}.
#' @param stagecol A vector of three strings, indicating the stage name columns
#' for times \emph{t}+1, \emph{t}, and \emph{t}-1, respectively, in \code{data}.
#' Defaults to \code{stagecol = c("stage3", "stage2", "stage1")}.
#' 
#' @return A data frame with five variables:
#' \item{rowid}{A string identifier term, equal to the monitoring occasion in
#' time \emph{t} and the stage index.}
#' \item{stageindex}{The stageframe index of the stage.}
#' \item{stage}{The name of each stage, or \code{NA}.}
#' \item{year2}{Monitoring occasion in time \emph{t}.}
#' \item{frequency}{The number of individuals in the respective stage and time.}
#' \item{actual_prop}{The proportion of individuals alive in time \emph{t} in
#' the respective stage.}
#' 
#' @section Notes:
#' This function produces frequencies and proportions of stages in hfv formatted
#' data using stage index variables rather than stage name variables, and so
#' requires the former. The latter is only required if the user wants to know
#' the associated stage names.
#' 
#' Frequencies and proportions will be calculated for all times, including the
#' last time, which is generally found in the \code{stage3} columns of the last
#' \code{year2} entry in object \code{data}. The default is to treat the
#' \code{year2} entry for that time as \code{max(year2) + 1}.
#' 
#' Note that no stageframe is required for this function to operate. Stage
#' names and their order are inferred directly from the object \code{data}.
#' 
#' @examples
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 3, 6, 11, 19.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
#'   "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1.5, 1.5, 3.5, 5)
#' comments <- c("Dormant seed", "1st yr protocorm", "2nd yr protocorm",
#'   "3rd yr protocorm", "Seedling", "Dormant adult",
#'   "Extra small adult (1 shoot)", "Small adult (2-4 shoots)",
#'   "Medium adult (5-7 shoots)", "Large adult (8-14 shoots)",
#'   "Extra large adult (>14 shoots)")
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector, 
#'   repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
#'   propstatus = propvector, immstatus = immvector, indataset = indataset, 
#'   binhalfwidth = binvec, comments = comments)
#' 
#' cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004, 
#'   patchidcol = "patch", individcol = "plantid", blocksize = 4,
#'   sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
#'   repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
#'   stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
#'   NRasRep = TRUE)
#' 
#' all_stage_props <- actualstage3(cypraw_v1)
#' all_stage_props
#' 
#' @export
actualstage3 <- function(data, historical = FALSE, year2 = "year2",
  indices = c("stage3index", "stage2index", "stage1index"),
  stagecol = c("stage3", "stage2", "stage1")) {
  
  aaa.data <- ordered_stages <- ordered_indices <- NULL
  stagenames <- FALSE
  
  if (length(indices) < 2) {
    stop("Object indices must contain the names of 3 variables corresponding to stage index
      in times t+1, t, and t-1, respectively", call. = FALSE)
  }
  if (!all(is.element(indices, names(data)))) {
    stop("Object indices must contain the names of 3 variables corresponding to stage index
      in times t+1, t, and t-1, respectively", call. = FALSE)
  }
  
  if (all(is.element(stagecol, names(data)))) {
    stagenames <- TRUE
  }
  
  if (length(year2) != 1) {
    stop("Object year2 must equal the name of the variable denoting monitoring occasion in time t",
      call. = FALSE)
  }
  if (!is.element(year2, names(data))) {
    stop("Object year2 must equal the name of the variable denoting monitoring occasion in time t",
      call. = FALSE)
  }
  stage3index <- indices[1]
  stage2index <- indices[2]
  stage1index <- indices[3]
  
  names(data)[which(names(data) == stage3index)] <- "stage3index"
  names(data)[which(names(data) == stage2index)] <- "stage2index"
  names(data)[which(names(data) == stage1index)] <- "stage1index"
  names(data)[which(names(data) == year2)] <- "year2"
  
  if (stagenames) {
    stage3name <- stagecol[1]
    stage2name <- stagecol[2]
    stage1name <- stagecol[3]
    names(data)[which(names(data) == stage3name)] <- "stage3"
    names(data)[which(names(data) == stage2name)] <- "stage2"
    names(data)[which(names(data) == stage1name)] <- "stage1"
  }
  
  data <- data[, c("year2", "stage1", "stage2", "stage3", "stage1index", "stage2index", "stage3index")]
  all_years <- sort(unique(data$year2), decreasing = TRUE)
  bits_to_tack_on <- data[which(data$year2 == all_years[1]),]
  bits_to_tack_on$stage1index <- bits_to_tack_on$stage2index
  bits_to_tack_on$stage2index <- bits_to_tack_on$stage3index
  bits_to_tack_on$stage1 <- bits_to_tack_on$stage2
  bits_to_tack_on$stage2 <- bits_to_tack_on$stage3
  bits_to_tack_on$year2 <- all_years[1] + 1
  data <- rbind.data.frame(data, bits_to_tack_on)
  
  if (!historical) {
    if (stagenames) {
      data$stages <- data$stage2
    }
    data$stageindex <- data$stage2index
    
    aaa.data <- as.data.frame(xtabs(~ stage2index + year2, data), stringsAsFactors = FALSE)
    names(aaa.data)[which(names(aaa.data) == "stage2index")] <- "stageindex"
    
    ordered_indices <- sort(unique(as.numeric(aaa.data$stageindex)))
  } else {
    data$stageindex <- apply(as.matrix(c(1:dim(data)[1])), 1, function(X) {
      paste(data$stage1index[X], data$stage2index[X])
    })
    
    if (stagenames) {
      data$stages <- apply(as.matrix(c(1:dim(data)[1])), 1, function(X) {
        paste(data$stage1[X], data$stage2[X])
      })
    }
    
    aaa.data <- as.data.frame(xtabs(~ stageindex + year2, data), stringsAsFactors = FALSE)
    
    ordered_indices <- sort(unique(aaa.data$stageindex))
  }
  
  ordered_stages <- apply(as.matrix(ordered_indices), 1, function(X) {
    acmecanning <- data$stages[which(data$stageindex == X)[1]]
    if (length(acmecanning) < 1) {
      acmecanning <- NA
    }
    return(acmecanning)
  })
  
  totalr <- as.data.frame(xtabs(~ year2, data), stringsAsFactors = FALSE)
  
  aaa.data$actual_prop <- apply(as.matrix(c(1:dim(aaa.data)[1])), 1, function(X) {
    a <- aaa.data$Freq[X] / totalr$Freq[which(totalr$year2 == aaa.data$year2[X])]
    return(a)
  })
  
  aaa.data$rowid <- apply(as.matrix(c(1:dim(aaa.data)[1])), 1, function(X) {
    paste(aaa.data$year2[X], aaa.data$stageindex[X])
  })
  
  if (stagenames) {
    aaa.data$stage <- apply(as.matrix(aaa.data$stageindex), 1, function(X) {
      return(ordered_stages[which(ordered_indices == X)])
    })
  } else {
    aaa.data$stage <- NA
  }
  
  aaa.data <- aaa.data[, c("rowid", "stageindex", "stage", "year2", "Freq",
    "actual_prop")]
  names(aaa.data)[which(names(aaa.data) == "Freq")] <- "frequency"
  
  return(aaa.data)
}

