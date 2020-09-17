#' Estimate Mean Projection Matrices
#' 
#' \code{lmean()} estimates mean projection matrices. The function differs from a
#' typical element-by-element mean matrix estimator through two options. First,
#' it allows mean matrix elements to be estimated as geometric means across time. 
#' Spatial means are always developed as arithmetic means. Second, it allows 
#' element means to be estimated ignoring 0s in cases where some elements are 
#' not zero. The function takes \code{lefkoMat} objects as input, and returns the 
#' same class of object.
#' 
#' @param mats A \code{lefkoMat} object holding population projection matrices.
#' @param time A variable designating whether element means be computed as 
#' geometric (\code{geometric} or \code{g}) or arithmetic (\code{arithmetic} or \code{a}) across 
#' time. Please note that using the \code{geometric} option does not yield a 
#' geometric mean matrix - it yields a matrix of geometric mean elements, which 
#' is theoretically different. Defaults to \code{arithmetic}. 
#' @param sparse If TRUE, then all 0s will be ignored in elements that include
#' other numbers across matrices. Only elements that equal 0 in all matrices
#' (structural zeroes) will be exempt. Defaults to FALSE.
#' @param AasSum If TRUE, then lefkoMat matrix means are estimated as means of
#' U and F matrices, and then summed to estimate A matrices. If FALSE, then
#' mean A matrices are estimated as means of A matrices within lefkoMat objects.
#' Defaults to TRUE.
#' 
#' @return Yields a \code{lefkoMat} object with the following characteristics:
#' 
#' \item{A}{A list of full mean projection matrices in order of sorted
#' populations,patches, and years. These are typically estimated as the sums
#' of associated mean \code{U} and \code{F} matrices.}
#' \item{U}{A list of mean survival-transition matrices sorted as in \code{A}.}
#' \item{F}{A list of mean fecundity matrices sorted as in \code{A}.}
#' \item{hstages}{A data frame matrix showing the pairing of ahistorical 
#' stages used to create historical stage pairs.}
#' \item{ahstages}{A data frame detailing the characteristics of associated
#' ahistorical stages.}
#' \item{labels}{A data frame giving the population, patch, and year of each 
#' mean matrix in order. If \code{pop}, \code{patch}, or \code{year2} are all NA in the
#' original \code{labels} set, then these will be re-labeled as \code{A}, \code{1}, or \code{1},
#' respectively.}
#' \item{matrixqc}{A short vector describing the number of non-zero elements
#' in \code{U} and \code{F} mean matrices, and the number of annual matrices.}
#' \item{modelqc}{The \code{qc} portion of the modelsuite input, if used.}
#' \item{dataqc}{A vector showing the numbers of individuals and rows in the
#' vertical dataset used as input.}
#' 
#' @examples
#' \donttest{
#' data(lathyrus)
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
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          size1col = "Volume88", repstr1col = "FCODE88", 
#'                          fec1col = "Intactseed88", dead1col = "Dead1988", 
#'                          nonobs1col = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' ehrlen3mean$A[[1]]
#' }
#' 
#' @export
lmean <- function(mats, time = "arithmetic", sparse = FALSE, AasSum = TRUE) {
  
  if (class(mats) != "lefkoMat") {
    stop("An object of class lefkoMat is required as input.")
  }
  
  time <- tolower(time)
  
  time.possible <- c("arithmetic", "a", "geometric", "g")
  
  if (!is.element(time, time.possible)) {
    stop("Invalid expression given for time. Please choose 'g' or 'geometric' for the geometric mean, or 'a' or 'arithmetic' for the arithmetic mean.")
  }
  
  if (time == "a") time <- "arithmetic"
  if (time == "g") time <- "geometric"
  
  #First we create an index based on the input lefkoMat object
  listofyears <- mats$labels
  
  if (all(is.na(listofyears$pop))) {
    listofyears$pop <- 1
  }
  
  if (all(is.na(listofyears$patch))) {
    listofyears$patch <- 1
  }
  
  if (all(is.na(listofyears$year2))) {
    listofyears$year2 <- 1
  }
  
  listofyears$poppatch <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {
    paste(listofyears$pop[X], listofyears$patch[X])
  })
  
  listofyears$popc <- apply(as.matrix(listofyears$pop), 1, function(X) {which(unique(listofyears$pop) == X)-1})
  listofyears$poppatchc <- apply(as.matrix(listofyears$poppatch), 1, function(X) {which(unique(listofyears$poppatch) == X)-1})
  listofyears$year2c <- apply(as.matrix(listofyears$year2), 1, function(X) {which(unique(listofyears$year2) == X)-1})
  
  listofyears$patchesinpop <- apply(as.matrix(c(1:length(listofyears$poppatchc))), 1, function(X) {length(unique(listofyears$poppatchc[which(listofyears$popc == listofyears$popc[X])]))})
  listofyears$yearsinpatch <- apply(as.matrix(c(1:length(listofyears$year2c))), 1, function(X) {length(unique(listofyears$year2c[which(listofyears$poppatchc == listofyears$poppatchc[X])]))})
  
  #  loy2c <- as.matrix(listofyears[,5:9])
  
  numofpops <- length(unique(listofyears$popc))
  numofpatches <- length(unique(listofyears$poppatchc))
  numofyears <- length(unique(listofyears$year2c))
  
  allmatricesU <- do.call("cbind", lapply(mats$U, as.vector))
  allmatricesU <- cbind(allmatricesU, rowSums(allmatricesU))
  
  allmatricesF <- do.call("cbind", lapply(mats$F, as.vector))
  allmatricesF <- cbind(allmatricesF, rowSums(allmatricesF))
  
  if (!AasSum) {
    allmatricesA <- do.call("cbind", lapply(mats$A, as.vector))
    allmatricesA <- cbind(allmatricesA, rowSums(allmatricesA))
  }
  
  if (time == "geometric") pushtime <- 1 else pushtime <- 0
  if (sparse == TRUE) pushsparse <- 1 else pushsparse <- 0
  
  if (!AasSum) {
    UFmats <- turbogeodiesel(listofyears, allmatricesU, allmatricesF, allmatricesA, pushtime, 
                             pushsparse, numofpops, numofpatches, numofyears)
  } else {
    UFmats <- geodiesel(listofyears, allmatricesU, allmatricesF, pushtime, pushsparse, numofpops, 
                        numofpatches, numofyears)
  }
  
  matsdim <- dim(mats$A[[1]])[1]
  
  meanUmats <- lapply(c(1:numofpatches, ((numofpatches*3) + 1):((numofpatches*3) + numofpops)), function(i) {
    matrix(UFmats[,i], ncol = matsdim, nrow = matsdim)}
  )
  
  meanFmats <- lapply(c((numofpatches + 1):(2 * numofpatches), ((numofpatches*3) + numofpops + 1):((numofpatches*3) + (2 * numofpops))), function(i) {
    matrix(UFmats[,i], ncol = matsdim, nrow = matsdim)}
  )
  
  meanAmats <- lapply(c((2 * numofpatches + 1):(3 * numofpatches), ((numofpatches*3) + (2 * numofpops) + 1):((numofpatches*3) + (3 * numofpops))), function(i) {
    matrix(UFmats[,i], ncol = matsdim, nrow = matsdim)}
  )
  
  listofyears$groupid <- listofyears$poppatch
  listofyears$sorter1 <- apply(as.matrix(c(1:dim(listofyears)[1])), 1, function(X) {
    paste(listofyears$pop[X], listofyears$patch[X], listofyears$year2[X])
  })
  
  poppatch.combos <- sort(unique(listofyears$poppatch))
  pop.combos <- sort(unique(listofyears$pop))
  sortorder1 <- apply(as.matrix(poppatch.combos), 1, function(X) {which(listofyears$poppatch == X)[1]})
  sortorder2 <- apply(as.matrix(pop.combos), 1, function(X) {which(listofyears$pop == X)[1]})
  
  finallabels <- rbind.data.frame(listofyears[sortorder1,c(1:2)], listofyears[sortorder2,c(1:2)])
  finallabels$patch[c(length(sortorder1) + 1):(length(sortorder1) + length(sortorder2))] <- "All"
  
  #Matrix QC
  qcoutput1 <- NA
  
  totalutransitions <- sum(unlist(lapply(meanUmats, function(X) {length(which(X != 0))})))
  totalftransitions <- sum(unlist(lapply(meanFmats, function(X) {length(which(X != 0))})))
  totalmatrices <- length(meanUmats)
  
  qcoutput1 <- c(totalutransitions, totalftransitions, totalmatrices)
  
  output <- NA
  
  if (is.element("modelqc", names(mats))) {
    output <- list(A = meanAmats, U = meanUmats, F = meanFmats, hstages = mats$hstages, 
                   ahstages = mats$ahstages, labels = finallabels, matrixqc = qcoutput1, 
                   modelqc = mats$modelqc)
  } else if (is.element("dataqc", names(mats))) {
    output <- list(A = meanAmats, U = meanUmats, F = meanFmats, hstages = mats$hstages, 
                   ahstages = mats$ahstages, labels = finallabels, matrixqc = qcoutput1, 
                   dataqc = mats$dataqc)
  }
  
  class(output) <- "lefkoMat"
  
  return(output)
}

#' Dominant Eigenvalue and Deterministic Population Growth Rate Estimation
#' 
#' \code{lambda3()} is a generic function that returns the dominant eigenvalue
#' of a matrix, and set of dominant eigenvalues of a set of matrices. Unlike the
#' 'popbio' package's \link[popbio]{lambda}() function, it is particularly
#' made to handle very large and sparse matrices supplied as \code{lefkoMat} 
#' objects or as individual matrices. This function can handle large and sparse 
#' matrices, and so can be used with large historical matrices, IPMs, 
#' age x stage matrices, as well as smaller ahistorical matrices.
#' 
#' @param mats A lefkoMat object, or a population projection matrix, for which
#' the dominant eigenvalue is desired.
#' 
#' @return The value returned depends on the class of the \code{mats} argument.
#' 
#' @seealso \code{\link{lambda3.lefkoMat}()}
#' @seealso \code{\link{lambda3.matrix}()}
#' 
#' @examples
#' \donttest{
#' data(lathyrus)
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
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          size1col = "Volume88", repstr1col = "FCODE88", 
#'                          fec1col = "Intactseed88", dead1col = "Dead1988", 
#'                          nonobs1col = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' lambda3(ehrlen3mean)
#' }
#' 
#' @export
lambda3 <- function(mats) UseMethod("lambda3")

#' Estimate Deterministic Population Growth Rate for a lefkoMat Object
#' 
#' \code{lambda3.lefkoMat()} returns the dominant eigenvalues of projection
#' matrices supplied as \code{lefkoMat} objects. This function can handle large 
#' and sparsematrices, and so can be used with large historical matrices, IPMs, 
#' age x stage matrices, as well as smaller ahistorical matrices.
#' 
#' @param mats An object of class \code{lefkoMat}.
#' 
#' @return This function returns the dominant eigenvalue of each \code{$A} matrix
#' in the lefkoMat object input. This is given as the largest real part of all 
#' eigenvalues estimated via the \code{\link[RSpectra]{eigs}()} function in package 
#' 'RSpectra'. The output includes a data frame showing the population, patch,
#' and lambda estimate for each \code{$A} matrix within the object. Row names
#' correspond to the number of the matrix within the \code{$A} element of the 
#' \code{lefkoMat} object.
#' 
#' @seealso \code{\link{lambda3}()}
#' @seealso \code{\link{lambda3.matrix}()}
#' 
#' @examples
#' \donttest{
#' data(lathyrus)
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
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          size1col = "Volume88", repstr1col = "FCODE88", 
#'                          fec1col = "Intactseed88", dead1col = "Dead1988", 
#'                          nonobs1col = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' lambda3(ehrlen3mean)
#' }
#' 
#' @export
lambda3.lefkoMat <- function(mats) {
  baldrick <- if (any(class(mats$A) == "matrix")) {
    
    lambda3.matrix(mats$A)
    
  } else if (class(mats$A) == "list") {
    
    unlist(lapply(mats$A, lambda3.matrix))
    
  } else {
    
    stop("Input not recognized.")
    
  }
  
  output <- cbind.data.frame(mats$labels, baldrick)
  rownames(output) <- c(1:length(baldrick))
  
  names(output)[length(names(output))] <- "lambda"
  
  return(output)
}

#' Estimate Deterministic Population Growth Rate of a Projection Matrix
#' 
#' \code{lambda3.matrix()} returns the dominant eigenvalues of a single
#' projection matrix. This function can handle large and sparse matrices, 
#' and so can be used with large historical matrices, IPMs, age x stage 
#' matrices, as well as smaller ahistorical matrices.
#' 
#' @param mats A population projection matrix of class \code{matrix}.
#'
#' @return This function returns the dominant eigenvalue of the matrix. This
#' is given as the largest real part of all eigenvalues estimated via the 
#' \code{\link[RSpectra]{eigs}()} function in package 'RSpectra'.
#' 
#' @seealso \code{\link{lambda3}()}
#' @seealso \code{\link{lambda3.lefkoMat}()}
#' 
#' @examples
#' \donttest{
#' data(lathyrus)
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
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          size1col = "Volume88", repstr1col = "FCODE88", 
#'                          fec1col = "Intactseed88", dead1col = "Dead1988", 
#'                          nonobs1col = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' lambda3(ehrlen3mean$A[[1]])
#' }
#' 
#' @export
lambda3.matrix <- function(mats)
{
  decomp <- RSpectra::eigs(mats, k = 1, which = "LR")
  
  lambda <- Re(decomp$values[1]);
  
  return(lambda)
}

#' Stable Stage Distribution Estimation
#' 
#' \code{stablestage3()} is a generic function that returns the stable stage 
#' distribution for a population projection matrix or set of matrices. Unlike 
#' the 'popbio' package's \code{\link[popbio]{stable.stage}()} function, it is 
#' particularly made to handle very large and sparse matrices supplied as 
#' \code{lefkoMat} objects or as individual matrices. This function can
#' handle large and sparse matrices, and so can be used with large historical
#' matrices, IPMs, age x stage matrices, as well as smaller ahistorical 
#' matrices.
#' 
#' @param mats A lefkoMat object, or population projection matrix, for which
#' the stable stage distribution is desired.
#' 
#' @return The value returned depends on the class of the \code{mats} argument.
#' 
#' @seealso \code{\link{stablestage3.lefkoMat}()}
#' @seealso \code{\link{stablestage3.matrix}()}
#' 
#' @examples
#' \donttest{
#' data(lathyrus)
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
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          size1col = "Volume88", repstr1col = "FCODE88", 
#'                          fec1col = "Intactseed88", dead1col = "Dead1988", 
#'                          nonobs1col = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' stablestage3(ehrlen3mean)
#' }
#' 
#' @export
stablestage3 <- function(mats) UseMethod("stablestage3")

#' Estimate Stable Stage Distribution for a lefkoMat Object
#' 
#' \code{stablestage3.lefkoMat()} returns the stable stage distributions for all
#' \code{$A} matrices in an object of class \code{lefkoMat}. This function can 
#' handle large and sparse matrices, and so can be used with large historical 
#' matrices, IPMs, age x stage matrices, as well as smaller ahistorical 
#' matrices.
#' 
#' @param mats An object of class \code{lefkoMat}.
#' 
#' @return This function returns the stable stage distributions corresponding
#' to the matrices in a \code{lefkoMat} object. The stable stage distribution is 
#' given as the right eigenvector associated with largest real part of all 
#' eigenvalues estimated via the \code{\link[RSpectra]{eigs}()} function in 
#' package 'RSpectra' divided by the sum of the associated right eigenvector. 
#' 
#' The actual output depends on whether the \code{lefkoMat} object used as input
#' is ahistorical or historical. If the former, then a single data frame is
#' output. This data frame includes the number of the matrix within the \code{$A} 
#' element of the input \code{lefkoMat} object, followed by the original stage id, 
#' the new stage id (numeric and assigned through \code{\link{sf_create}()}), the 
#' original given size, and the estimated proportion of the stable stage 
#' distribution within a variable called \code{ss_prop}.
#' 
#' If a historical matrix is used as input, then two data frames are output
#' into a list object. The \code{$hist} element contains a data frame where the 
#' stable stage distribution is given in terms of across-year stage pairs.
#' The structure includes the matrix number, the original and new 
#' designations for stages in times \emph{t} and \emph{t}-1, respectively, followed
#' by the estimated proportion of the stable stage distribution for that matrix. 
#' The \code{$ahist} element contains the stable stage distribution in stages
#' as given in the original stageframe. It includes a data frame with the new 
#' stage designation from the associated stageframe for the \code{lefkoMat}
#' object, the corresponding matrix, and the stable stage distribution estimated 
#' as the sum of distribution elements from \code{$hist} corresponding to the 
#' equivalent stage in time \emph{t}, irrespective of stage in time \emph{t}-1.
#' 
#' @seealso \code{\link{stablestage3}()}
#' @seealso \code{\link{stablestage3.matrix}()}
#' 
#' @examples
#' \donttest{
#' data(lathyrus)
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
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          size1col = "Volume88", repstr1col = "FCODE88", 
#'                          fec1col = "Intactseed88", dead1col = "Dead1988", 
#'                          nonobs1col = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' stablestage3(ehrlen3mean)
#' }
#' 
#' @export
stablestage3.lefkoMat <- function(mats) {
  baldrick <- if (any(class(mats$A) == "matrix")) {
    
    stablestage3.matrix(mats$A)
    
  } else if (class(mats$A) == "list") {
    
    unlist(lapply(mats$A, stablestage3.matrix))
    
  } else {
    
    stop("Input not recognized.")
    
  }
  
  if (class(mats$A) == "list") {
    multiplier <- length(mats$A)
  } else multiplier <- 1
  
  if (all(is.na(mats$hstages))) {
    labels <- mats$ahstages[,1:3]
    
    modlabels <- cbind.data.frame(as.matrix(rep(c(1:multiplier), each = dim(labels)[1])), do.call("rbind.data.frame", apply(as.matrix(c(1:multiplier)), 1, function(X) return(labels))))
    
    output <- cbind.data.frame(modlabels, baldrick)
    names(output) <- c("matrix", "new_stage_id", "orig_stage_id", "original_size", "ss_prop")
    rownames(output) <- c(1:dim(output)[1])
  } else {
    labels <- mats$hstages
    
    modlabels <- cbind.data.frame(as.matrix(rep(c(1:multiplier), each = dim(labels)[1])), do.call("rbind.data.frame", apply(as.matrix(c(1:multiplier)), 1, function(X) return(labels))))
    
    outputh <- cbind.data.frame(modlabels, baldrick)
    names(outputh) <- c("matrix", "orig_stage_id_2", "orig_stage_id_1", "new_stage_id_2", "new_stage_id_1", "ss_prop")
    rownames(outputh) <- c(1:dim(outputh)[1])
    
    ahlabels <- mats$ahstages[,"new_stage_id"]
    ss2 <- c(apply(as.matrix(c(1:multiplier)), 1, function(X) {
      rightset <- subset(outputh, matrix == X)
      apply(as.matrix(ahlabels), 1, function(Y) {
        sum(rightset$ss_prop[which(rightset$new_stage_id_2 == Y)])
      })
    }))
    outputah <- cbind.data.frame(rep(c(1:multiplier), each = length(ahlabels)), rep(ahlabels, multiplier), ss2)
    names(outputah) <- c("matrix", "new_stage_id", "ss_prop")
    rownames(outputah) <- c(1:dim(outputah)[1])
    
    output <-list(hist = outputh, ahist = outputah)
  }
  
  return(output)
}

#' Estimate Stable Stage Distribution for a Population Projection Matrix
#' 
#' \code{stablestage3.matrix()} returns the stable stage distribution for a 
#' population projection matrix. This function can handle large and sparse 
#' matrices, and so can be used with large historical matrices, IPMs, 
#' age x stage matrices, as well as smaller ahistorical matrices.
#' 
#' @param mats A population projection matrix of class \code{matrix}.
#' 
#' @return This function returns the stable stage distribution corresponding to
#' the input matrix. The stable stage distribution is given as the right 
#' eigenvector associated with largest real part of the eigenvalues estimated 
#' for the matrix via the \code{\link[RSpectra]{eigs}()} function in package 
#' 'RSpectra', divided by the sum of the associated right eigenvector. 
#' 
#' @seealso \code{\link{stablestage3}()}
#' @seealso \code{\link{stablestage3.lefkoMat}()}
#' 
#' @examples
#' \donttest{
#' data(lathyrus)
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
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          size1col = "Volume88", repstr1col = "FCODE88", 
#'                          fec1col = "Intactseed88", dead1col = "Dead1988", 
#'                          nonobs1col = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ")
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' stablestage3(ehrlen3mean$A[[1]])
#' }
#' 
#' @export
stablestage3.matrix <- function(mats)
{
  decomp <- RSpectra::eigs(mats, k = 1, which = "LR")
  
  w <- zapsmall(Re(decomp$vectors[,1]))
  wcorr <- w / sum(w)
  
  return(wcorr)
}

#' Reproductive Value Estimation
#' 
#' \code{repvalue3()} is a generic function that estimates returns the 
#' reproductive values of stages in a population projection matrix or a set of
#' matrices. The specifics of estimation vary with the class of input object.
#' However, unlike the 'popbio' package's \code{\link[popbio]{reproductive.value}()}
#' function, it is particularly made to handle very large and sparse matrices
#' supplied as \code{lefkoMat} objects or as individual matrices. This function
#' can handle large and sparse matrices, and so can be used with large
#' historical matrices, IPMs, age x stage matrices, as well as smaller
#' ahistorical matrices.
#' 
#' @param mats A lefkoMat object, or population projection matrix.
#' 
#' @return The value returned depends on the class of the \code{mats} argument.
#' 
#' @seealso \code{\link{repvalue3.lefkoMat}()}
#' @seealso \code{\link{repvalue3.matrix}()}
#' 
#' @examples
#' \donttest{
#' data(lathyrus)
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
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          size1col = "Volume88", repstr1col = "FCODE88", 
#'                          fec1col = "Intactseed88", dead1col = "Dead1988", 
#'                          nonobs1col = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ", reduce = TRUE)
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' repvalue3(ehrlen3mean)
#' }
#' 
#' @export
repvalue3 <- function(mats) UseMethod("repvalue3")

#' Estimate Reproductive Value for a lefkoMat Object
#' 
#' \code{repvalue3.lefkoMat()} returns the reproductive values for stages in a
#' set of population projection matrices. This function can handle large and 
#' sparse matrices, and so can be used with large historical matrices, IPMs, 
#' age x stage matrices, as well as smaller ahistorical matrices.
#' 
#' @param mats An object of class \code{lefkoMat} object.
#' 
#' @return This function returns the reproductive values for stages of a
#' population projection matrix. The nature of the output depends on whether the
#' lefkoMat object used as input is ahistorical or historical. In both cases,
#' raw reproductive values are estimated as the left eigenvector associated with
#' the largest real part of the dominant eigenvalue estimated via the \code{\link[RSpectra]{eigs}()}
#' function in package 'RSpectra', divided either by the first non-zero element
#' of the left eigenvector.
#' 
#' If an ahistorical matrix set is used as input, then the output is a data
#' frame that includes the number of the matrix within the `$A` element of the 
#' input \code{lefkoMat} object, followed by the original stage id, the new
#' stage id (numeric and assigned through \code{\link{sf_create}())}, the original 
#' given size, and the reproductive value estimate within a variable called
#' \code{repvalue}. 
#' 
#' If a historical matrix set is used as input, then a list with two elements is
#' output. The first element is a data frame showing the reproductive values
#' given in terms of across-year stage pairs, as estimated in the procedure
#' described above. The second element is another data frame showing the
#' reproductive values of the basic stages in the associated stageframe. The 
#' reproductive values in this second data frame are estimated via the approach
#' developed in Ehrlen (2000), in which each ahistorical stage's reproductive 
#' value is the average of the RVs summed by stage at time \emph{t} weighted by 
#' the proportion of that stage pair within the historical stable stage 
#' distribution associated with the matrix.
#' 
#' @seealso \code{\link{repvalue3}()}
#' @seealso \code{\link{repvalue3.matrix}()}
#' 
#' @examples
#' \donttest{
#' data(lathyrus)
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
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          size1col = "Volume88", repstr1col = "FCODE88", 
#'                          fec1col = "Intactseed88", dead1col = "Dead1988", 
#'                          nonobs1col = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ", reduce = TRUE)
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' repvalue3(ehrlen3mean)
#' }
#' 
#' @export
repvalue3.lefkoMat <- function(mats) {
  
  if (any(class(mats$A) == "matrix")) {
    
    decomp <- RSpectra::eigs(t(mats$A), k = 1, which = "LR")
    
    baldrick <- zapsmall(Re(decomp$vectors[,1]))
    
  } else if (class(mats$A) == "list") {
    
    baldrick <- unlist(lapply(mats$A, function(X) {
      
      decomp <- RSpectra::eigs(t(X), k = 1, which = "LR")
      
      v <- zapsmall(Re(decomp$vectors[,1]))
      
      return(v)
      
    }))
    
  } else {
    
    stop("Input not recognized.")
    
  }
  
  if (all(is.na(mats$hstages))) {
    labels <- mats$ahstages[,1:3]
  } else {
    labels <- mats$hstages
  }
  
  if (class(mats$A) == "list") {
    multiplier <- length(mats$A)
  } else multiplier <- 1
  
  if (all(is.na(mats$hstages))) {
    labels <- mats$ahstages[,1:3]
    
    modlabels <- cbind.data.frame(as.matrix(rep(c(1:multiplier), each = dim(labels)[1])), do.call("rbind.data.frame", apply(as.matrix(c(1:multiplier)), 1, function(X) return(labels))))
    
    output <- cbind.data.frame(modlabels, baldrick)
    names(output) <- c("matrix", "new_stage_id", "orig_stage_id", "original_size", "left_vector")
    rownames(output) <- c(1:dim(output)[1])
    
    output$rep_value <- apply(as.matrix(c(1:dim(output)[1])), 1, function(X) {
      matsub <- subset(output, matrix == output$matrix[X])
      entrystage <- min(which(abs(matsub$left_vector) > 0))
      return(output$left_vector[X] / matsub$left_vector[entrystage])
    })
    
  } else {
    ss3 <- stablestage3(mats)
    rahist <- ss3$ahist
    rhist <-ss3$hist
    rhist$ss3sum <- apply(as.matrix(c(1:dim(rhist)[1])), 1, function(X) {
      rahist$ss_prop[intersect(which(rahist$new_stage_id == rhist$new_stage_id_2[X]), 
                               which(rahist$matrix == rhist$matrix[X]))]
    })
    rhist$sscorr <- rhist$ss_prop / rhist$ss3sum
    rhist$sscorr[which(is.na(rhist$sscorr))] <- 0
    rhist$rv3raw <- baldrick
    
    rhist$rep_value_unc <- apply(as.matrix(c(1:dim(rhist)[1])), 1, function(X) {
      matsub <- subset(rhist, matrix == rhist$matrix[X])
      entrystage <- min(which(abs(matsub$rv3raw) > 0))
      return(rhist$rv3raw[X] / matsub$rv3raw[entrystage])
    })
    
    rhist$rv3 <- rhist$sscorr * rhist$rep_value_unc
    names(rhist)[which(names(rhist) == "rv3raw")] <- "left_vector"
    names(rhist)[which(names(rhist) == "rv3")] <- "rep_value"
    outputh <- rhist[,c(1,2,3,4,5,9,11)]
    
    ahlabels <- mats$ahstages[,"new_stage_id"]
    rv2 <- Re(c(apply(as.matrix(c(1:multiplier)), 1, function(X) {
      rightset <- subset(outputh, matrix == X)
      apply(as.matrix(ahlabels), 1, function(Y) {
        sum(rightset$rep_value[which(rightset$new_stage_id_2 == Y)])
      })
    })))
    outputah <- cbind.data.frame(rep(c(1:multiplier), each = length(ahlabels)), rep(ahlabels, multiplier), rv2)
    names(outputah) <- c("matrix", "new_stage_id", "rep_value_unc")
    rownames(outputah) <- c(1:dim(outputah)[1])
    
    outputah$rep_value <- apply(as.matrix(c(1:dim(outputah)[1])), 1, function(X) {
      matsub <- subset(outputah, matrix == outputah$matrix[X])
      entrystage <- min(which(abs(matsub$rep_value_unc) > 0))
      return(outputah$rep_value_unc[X] / matsub$rep_value_unc[entrystage])
    })
    outputah <- outputah[,c(1, 2, 4)]
    
    output <-list(hist = outputh, ahist = outputah)
    
  }
  
  return(output)
}

#' Estimate Reproductive Value for a Population Projection Matrix
#' 
#' \code{repvalue3.matrix()} returns the reproductive values for stages in a 
#' population projection matrix. The function assumes that the matrix is
#' ahistorical and provides standard reproductive values, meaning that
#' the overall reproductive values of basic life history stages in a historical
#' matrix are not provided (the \code{\link{repvalue3.lefkoMat}()} function estimates
#' these on the basis of stage description information provided in the 
#' \code{lefkoMat} object used as input in that function). This function can
#' handle large and sparse matrices, and so can be used with large historical
#' matrices, IPMs, age x stage matrices, as well as smaller ahistorical 
#' matrices.
#' 
#' @param mats A population projection matrix.
#' 
#' @return This function returns a vector data frame characterizing the 
#' reproductive values for stages of a population projection matrix. This is 
#' given as the left eigenvector associated with largest real part of the
#' dominant eigenvalue estimated via the \code{\link[RSpectra]{eigs}()} function in package 
#' 'RSpectra', divided by the first non-zero element of the left eigenvector. 
#' 
#' @seealso \code{\link{repvalue3}()}
#' @seealso \code{\link{repvalue3.lefkoMat}()}
#' 
#' @examples
#' \donttest{
#' data(lathyrus)
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
#' lathframe <- sf_create(sizes = sizevector, stagenames = stagevector, repstatus = repvector, 
#'                        obsstatus = obsvector, matstatus = matvector, immstatus = immvector, 
#'                        indataset = indataset, binhalfwidth = binvec, propstatus = propvector)
#' 
#' lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, patchidcol = "SUBPLOT", 
#'                          individcol = "GENET", blocksize = 9, juvcol = "Seedling1988", 
#'                          size1col = "Volume88", repstr1col = "FCODE88", 
#'                          fec1col = "Intactseed88", dead1col = "Dead1988", 
#'                          nonobs1col = "Dormant1988", stageassign = lathframe, 
#'                          stagesize = "sizea", censorcol = "Missing1988", 
#'                          censorkeep = NA, censor = TRUE)
#' 
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathover3 <- overwrite(stage3 = c("Sd", "Sd", "Sdl"), stage2 = c("Sd", "Sd", "Sd"), 
#'                        stage1 = c("Sd", "rep", "rep"), givenrate = c(0.345, 0.054))
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = c(1989, 1990), 
#'                    stages = c("stage3", "stage2", "stage1"), repmatrix = lathrepm, 
#'                    overwrite = lathover3, yearcol = "year2", 
#'                    indivcol = "individ", reduce = TRUE)
#' 
#' ehrlen3mean <- lmean(ehrlen3)
#' repvalue3(ehrlen3mean$A[[1]])
#' }
#' 
#' @export
repvalue3.matrix <- function(mats)
{
  decomp <- RSpectra::eigs(t(mats), k = 1, which = "LR")
  
  v <- zapsmall(Re(decomp$vectors[,1]))
  
  vcorr <- v / v[min(which(v != 0))]
  
  return(vcorr)
}


