#' Estimate Mean Projection Matrices
#' 
#' \code{lmean()} estimates mean projection matrices as element-wise arithmetic
#' means.
#' 
#' @param mats A \code{lefkoMat} object.
#' @param matsout A string identifying which means to estimate. Option
#' \code{"pop"} indicates population-level only, \code{"patch"} indicates
#' patch-level only, and \code{"all"} indicates that both patch- and
#' population-level means should be estimated. Defaults to \code{"all"}.
#' 
#' @return Yields a \code{lefkoMat} object with the following characteristics:
#' 
#' \item{A}{A list of full mean projection matrices in order of sorted
#' populations, patches, and years. These are typically estimated as the sums of
#' the associated mean \code{U} and \code{F} matrices. All matrices output in
#' the \code{matrix} class.}
#' \item{U}{A list of mean survival-transition matrices sorted as in \code{A}.
#' All matrices output in the \code{matrix} class.}
#' \item{F}{A list of mean fecundity matrices sorted as in \code{A}. All
#' matrices output in the \code{matrix} class.}
#' \item{hstages}{A data frame showing the pairing of ahistorical stages used to
#' create historical stage pairs. Given if the MPM is historical.}
#' \item{ahstages}{A data frame detailing the characteristics of associated
#' ahistorical stages.}
#' \item{labels}{A data frame detailing the order of population, patch, and year 
#' of each mean matrix. If \code{pop}, \code{patch}, or \code{year2} are NA in
#' the original \code{labels} set, then these will be re-labeled as \code{A},
#' \code{1}, or \code{1}, respectively.}
#' \item{matrixqc}{A short vector describing the number of non-zero elements in
#' \code{U} and \code{F} mean matrices, and the number of annual matrices.}
#' \item{modelqc}{This is the \code{qc} portion of the \code{modelsuite} input.
#' Only output from \code{lefkoMat} objects resulting from function-based
#' estimation.}
#' \item{dataqc}{A vector showing the numbers of individuals and rows in the
#' vertical dataset used as input. Only output from \code{lefkoMat} objects
#' resulting from raw matrix estimation.}
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
lmean <- function(mats, matsout = "all") {
  
  agestages <- NA
  
  if (class(mats) != "lefkoMat") {
    stop("An object of class lefkoMat is required as input.")
  }
  
  if (is.element("agestages", names(mats))) {
    agestages <- mats$agestages
  }
  
  matsout.possible <- c("all", "pop", "patch")
  
  if (!is.element(tolower(matsout), matsout.possible)) {
    stop("The matsout option must take a value of all, pop, or patch.")
  }
  
  #First we create an index based on the input lefkoMat object
  listofyears <- mats$labels
  
  if (all(is.na(listofyears$pop))) {
    listofyears$pop <- "1"
  }
  
  if (all(is.na(listofyears$patch))) {
    listofyears$patch <- "1"
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
  
  numofpops <- length(unique(listofyears$popc))
  numofpatches <- length(unique(listofyears$poppatchc))
  numofyears <- length(unique(listofyears$year2c))
  
  listofyears$poppatchc <- as.numeric(listofyears$poppatchc)
  
  if (matsout == "all") {
    poponly <- 1;
    patchonly <- 1;
  } else if (matsout == "patch") {
    poponly <- 0;
    patchonly <- 1;
  } else {
    poponly <- 1;
    patchonly <- 0;
  }
  
  if (length(unique(listofyears$poppatchc)) == 1) {
    poponly <- 0
    patchonly <- 1
  }
  
  if (!all(is.na(mats$hstages))) {
    output <- .turbogeodiesel(listofyears, mats$U, mats$F, mats$hstages, 
      agestages, mats$ahstages, patchonly, poponly)
  } else {
    output <- .geodiesel(listofyears, mats$U, mats$F, agestages, mats$ahstages,
      patchonly, poponly)
    output$hstages <- NA
  }
  
  if (is.element("dataqc", names(mats))) {
    output$dataqc <- mats$dataqc
  }
  if (is.element("modelqc", names(mats))) {
    output$dataqc <- mats$modelqc
  }
  
  class(output) <- "lefkoMat"
  
  return(output)
}

#' Estimate Dominant Eigenvalue and Deterministic Population Growth Rate
#' 
#' \code{lambda3()} is a generic function that returns the dominant eigenvalue
#' of a matrix, and set of dominant eigenvalues of a set of matrices. It can
#' handle very large and sparse matrices supplied as \code{lefkoMat} objects or
#' as individual matrices, and can be used with large historical matrices, IPMs, 
#' age x stage matrices, as well as smaller ahistorical matrices.
#' 
#' @param mats A lefkoMat object, or a single projection matrix, for which the
#' dominant eigenvalue is desired.
#' @param ... Other parameters.
#' 
#' @return The value returned depends on the class of the \code{mats} argument.
#' 
#' @seealso \code{\link{lambda3.lefkoMat}()}
#' @seealso \code{\link{lambda3.matrix}()}
#' @seealso \code{\link{slambda3}()}
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
#' lambda3(ehrlen3mean)
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
#' lambda3(cypmatrix2r)
#' 
#' @export
lambda3 <- function(mats, ...) UseMethod("lambda3")

#' Estimate Deterministic Population Growth Rates of lefkoMat Matrices
#' 
#' \code{lambda3.lefkoMat()} returns the dominant eigenvalues of all projection
#' matrices supplied within \code{lefkoMat} objects. This function can handle
#' large and sparse matrices, and so can be used with large historical matrices,
#' IPMs, age x stage matrices, as well as ahistorical matrices.
#' 
#' @param mats An object of class \code{lefkoMat}.
#' @param sparse A text string indicating whether to use sparse matrix encoding
#' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
#' \code{"auto"}.
#' @param ... Other parameters.
#' 
#' @return This function returns the dominant eigenvalue of each \code{$A}
#' matrix in \code{mats}. The output includes a data frame showing the
#' population, patch, and lambda estimate for each \code{A} matrix. Row names
#' correspond to the order of the matrix within the \code{A} element of
#' \code{mats}.
#' 
#' @section Notes:
#' The \code{sparse} option allows the function to utilize underlying methods
#' of either dense or sparse matrix manipulation in order to speed up processing
#' time and prevent memory shortages. Under the \code{auto} setting, the
#' function will determine whether the matrix is sparse and act accordingly.
#' For extremely large, sparse matrices, the user may simply set
#' \code{sparse = "yes"} to save time further and force the use of sparse format
#' in calculations.
#' 
#' @seealso \code{\link{lambda3}()}
#' @seealso \code{\link{lambda3.matrix}()}
#' @seealso \code{\link{slambda3}()}
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
#' lambda3(ehrlen3mean)
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
#' lambda3(cypmatrix2r)
#' 
#' @export
lambda3.lefkoMat <- function(mats, sparse = "auto", ...) {
  
  if (is.element(tolower(sparse), c("y", "yes", "yea", "yeah", "t", "true"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(sparse), c("n", "no", "non", "nah", "f", "false"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats$A[[1]])
    dense_elements <- length(which(mats$A[[1]] != 0))
    
    if ((dense_elements / elements_total) < 0.5) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  
  baldrick <- if (any(class(mats$A) == "matrix")) {
    .lambda3matrix(mats$A, sparsemethod)

  } else if (class(mats$A) == "list") {
    unlist(lapply(mats$A, .lambda3matrix, sparsemethod))

  } else {
    stop("Input not recognized.")
  }
  
  if (length(mats$labels$pop) == 2 & length(baldrick == 1)) {
    output <- cbind.data.frame(pop = 1, patch = 1, baldrick)
  } else {
    output <- cbind.data.frame(mats$labels, baldrick)
    rownames(output) <- c(1:length(baldrick))
  }
  
  names(output)[length(names(output))] <- "lambda"
  
  return(output)
}

#' Estimate Deterministic Population Growth Rate of Single Projection Matrix
#' 
#' \code{lambda3.matrix()} returns the dominant eigenvalue of a single
#' projection matrix. This function can handle large and sparse matrices, so can
#' be used with large historical matrices, IPMs, age x stage matrices, as well
#' as ahistorical matrices.
#' 
#' @param mats A population projection matrix of class \code{matrix}.
#' @param sparse A text string indicating whether to use sparse matrix encoding
#' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
#' \code{"auto"}.
#' @param ... Other parameters.
#' 
#' @return This function returns the dominant eigenvalue of the matrix.
#' 
#' @section Notes:
#' The \code{sparse} option allows the function to utilize underlying methods
#' of either dense or sparse matrix manipulation in order to speed up processing
#' time and prevent memory shortages. Under the \code{auto} setting, the
#' function will determine whether the matrix is sparse and act accordingly.
#' For extremely large sparse matrices, the user may simply set
#' \code{sparse = "yes"} to save time further and force the use of sparse format
#' in calculations.
#' 
#' @seealso \code{\link{lambda3}()}
#' @seealso \code{\link{lambda3.lefkoMat}()}
#' @seealso \code{\link{slambda3}()}
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
#' lambda3(ehrlen3mean$A[[1]])
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
#' lambda3(cypmatrix2r$A[[1]])
#' 
#' @export
lambda3.matrix <- function(mats, sparse = "auto", ...)
{
  if (is.element(tolower(sparse), c("y", "yes", "yea", "yeah", "t", "true"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(sparse), c("n", "no", "non", "nah", "f", "false"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats)
    dense_elements <- length(which(mats != 0))
    
    if ((dense_elements / elements_total) < 0.5) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  
  lambda <- .lambda3matrix(mats, sparsemethod);

  return(lambda)
}

#' Estimate Stable Stage Distribution
#' 
#' \code{stablestage3()} is a generic function that returns the stable stage 
#' distribution for a population projection matrix or set of matrices. This
#' function is made to handle very large and sparse matrices supplied as 
#' \code{lefkoMat} objects or as individual matrices, and can be used with large
#' historical matrices, IPMs, age x stage matrices, as well as ahistorical
#' matrices.
#' 
#' @param mats A lefkoMat object, or population projection matrix, for which the
#' stable stage distribution is desired.
#' @param ... Other parameters.
#' 
#' @return The value returned depends on the class of the \code{mats} argument.
#' See related functions for details.
#' 
#' @seealso \code{\link{stablestage3.lefkoMat}()}
#' @seealso \code{\link{stablestage3.matrix}()}
#' 
#' @examples
#' # Lathyrus deterministic example
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
#' stablestage3(ehrlen3mean)
#' 
#' # Cypripedium stochastic example
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
#' stablestage3(cypmatrix2r, stochastic = TRUE)
#' 
#' @export
stablestage3 <- function(mats, ...) UseMethod("stablestage3")

#' Estimate Stable Stage Distribution of Matrices in lefkoMat Object
#' 
#' \code{stablestage3.lefkoMat()} returns the deterministic stable stage
#' distributions of all \code{A} matrices in an object of class \code{lefkoMat},
#' as well as the long-run projected mean stage distribution in stochastic
#' analysis. This function can handle large and sparse matrices, and so can be
#' used with large historical matrices, IPMs, age x stage matrices, as well as
#' ahistorical matrices.
#' 
#' @param mats An object of class \code{lefkoMat}.
#' @param stochastic A logical value indicating whether to use deterministic
#' (\code{FALSE}) or stochastic (\code{TRUE}) analysis. Defaults to
#' \code{FALSE}.
#' @param times An integer variable indicating number of occasions to project if
#' using stochastic analysis. Defaults to 10000.
#' @param tweights An optional vector indicating the probability weighting to
#' use for each matrix in stochastic simulations. If not given, then defaults to
#' equal weighting.
#' @param seed A number to use as a random number seed.
#' @param sparse A text string indicating whether to use sparse matrix encoding
#' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
#' \code{"auto"}.
#' @param ... Other parameters.
#' 
#' @return This function returns the stable stage distributions (and long-run
#' mean stage distributions in stochastic analysis) corresponding to the
#' matrices in a \code{lefkoMat} object.
#' 
#' The output depends on whether the \code{lefkoMat} object used as input is
#' ahistorical or historical, and whether the analysis is deterministic or
#' stochastic. If deterministic and ahistorical, then a single data frame is
#' output, which includes the number of the matrix within the \code{A} element
#' of the input \code{lefkoMat} object, followed by the stage id (numeric and
#' assigned through \code{\link{sf_create}()}), the stage name, and the
#' estimated proportion of the stable stage distribution (\code{ss_prop}). If
#' stochastic and ahistorical, then a single data frame is output starting with
#' the number of the population-patch (\code{matrix_set}), a string
#' concatenating the names of the population and the patch (\code{poppatch}),
#' the assigned stage id number (\code{stage_id}), and the stage name
#' (\code{stage}), and the long-run average stage distribution (\code{ss_prop}).
#' 
#' If a historical matrix is used as input, then two data frames are output
#' into a list object. The \code{hist} element describes the historical
#' stage-pair distribution, while the \code{ahist} element describes the stage
#' distribution. If deterministic, then \code{hist} contains a data frame
#' including the matrix number (\code{matrix}), the numeric stage designations for
#' stages in occasions \emph{t} and \emph{t}-1, (\code{stage_id_2} and
#' \code{stage_id_1}, respectively), followed by the respective stage names (
#' \code{stage_2} and \code{stage_1}), and ending with the estimated stable
#' stage-pair distribution. The associated \code{ahist} element is as before. If
#' stochastic, then the \code{hist} element contains a single data frame with
#' the number of the population-patch (\code{matrix_set}), a string
#' concatenating the names of the population and the patch (\code{poppatch}),
#' the assigned stage id numbers in times \emph{t} and \emph{t}-1 (
#' \code{stage_id_2} and \code{stage_id_2}, respectively), and the associated
#' stage names (\code{stage_2} and \code{stage_1}, respectively), and the
#' long-run average stage distribution (\code{ss_prop}). The associated
#' \code{ahist} element is as before in the ahistorical, stochastic case.
#'
#' In addition to the data frames noted above, stochastic analysis will result
#' in the additional output of a list of matrices containing the actual
#' projected stage distributions across all projected occasions, in the order of
#' population-patch combinations in the \code{lefkoMat} input.
#'
#' @section Notes:
#' In stochastic analysis, the projected mean distribution is the arithmetic
#' mean across the final 1000 projected occasions if the simulation is at least
#' 2000 projected occasions long. If between 500 and 2000 projected occasions
#' long, then only the final 200 are used, and if fewer than 500 occasions are
#' used, then all are used. Note that because stage distributions in stochastic
#' simulations can change greatly in the initial portion of the run, we
#' encourage a minimum of 2000 projected occasions per simulation, with 10000
#' preferred.
#' 
#' @seealso \code{\link{stablestage3}()}
#' @seealso \code{\link{stablestage3.matrix}()}
#' 
#' @examples
#' # Lathyrus deterministic example
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
#' stablestage3(ehrlen3mean)
#' 
#' # Cypripedium stochastic example
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
#' stablestage3(cypmatrix2r, stochastic = TRUE)
#' 
#' @export
stablestage3.lefkoMat <- function(mats, stochastic = FALSE, times = 10000,
  tweights = NA, seed = NA, sparse = "auto", ...) {
  
  matrix_set <- NULL
  
  if (is.element(tolower(sparse), c("y", "yes", "yea", "yeah", "t", "true"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(sparse), c("n", "no", "non", "nah", "f", "false"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats$A[[1]])
    dense_elements <- length(which(mats$A[[1]] != 0))
    
    if ((dense_elements / elements_total) < 0.5) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  
  if (!stochastic) {
    baldrick <- if (any(class(mats$A) == "matrix")) {
      .ss3matrix(mats$A, sparsemethod)
      
    } else if (class(mats$A) == "list") {
      unlist(lapply(mats$A, .ss3matrix, sparsemethod))
      
    } else {
      stop("Input not recognized.")
    }
  } else {
    
    if (!is.na(seed)) {
      set.seed(seed)
    }
    
    mats$labels$poppatch <- paste(mats$labels$pop, mats$labels$patch)
    used_poppatches <- as.list(unique(mats$labels$poppatch))
    
    #Here we get the full stage distribution series for all occasions, as a list
    princegeorge <- lapply(used_poppatches, function(X) {
      used_slots <- which(mats$labels$poppatch == X)
      
      if (length(used_slots) < 2) {
        warning("Only 1 annual matrix found for some population-patch combinations. Stochastic analysis requires multiple annual matrices per population-patch combination.", call. = FALSE)
      }
      
      if (!is.na(tweights)) {
        if (length(tweights) != length(used_slots)) {
          if (length(tweights) == length(mats$A)) {
            used_weights <- tweights[used_slots] / sum(tweights[used_slots])
          } else {
            stop("Option tweights must be either NA, or a numeric vector equal to the number of years supplied or matrices supplied.",
              call. = FALSE)
          }
        } else {
          used_weights <- tweights / sum(tweights)
        }
      } else {
        used_weights <- rep(1, length(used_slots))
        used_weights <- used_weights / sum(used_weights)
      }
      
      theprophecy <- sample(used_slots, times, replace = TRUE, prob = used_weights) - 1
      starter <- .ss3matrix(mats$A[[used_slots[1]]], sparsemethod)
      
      theseventhmatrix <- .proj3(starter, mats$A, theprophecy, 1, 0, 0)
      
      ssonly <- theseventhmatrix[((dim(mats$A[[1]])[1]) + 1):(2 *(dim(mats$A[[1]])[1])),]
      
      return(ssonly)
    })
    
    # Now we create the mean distributions
    baldrick <- unlist(
      lapply(princegeorge, function(X) {
        if (times > 2000) {
          usedX <- X[,(times - 999):(times)]
        } else if (times > 500) {
          usedX <- X[,(times-199):(times)]
        } else {
          usedX <- X
        }
        apply(usedX, 1, mean)
      })
    )
  }
  
  # The final bits sort everything and clean it up, and create the ahistorical
  # version if a historical input was used
  if (class(mats$A) == "list") {
    if (!stochastic) {
      multiplier <- length(mats$A)
    } else {
      multiplier <- length(used_poppatches)
    }
  } else multiplier <- 1
  
  if (all(is.na(mats$hstages))) {
    labels_orig <- mats$ahstages[,1:2]
    mat_dims <- dim(mats$A[[1]])[1]
    if (dim(labels_orig)[1] == mat_dims) {
      labels <- labels_orig
    } else {
      newmult <- mat_dims / dim(labels_orig)[1]
      if (mat_dims %% dim(labels_orig)[1] != 0) {
        stop("Matrices do not appear to be ahistorical, historical, or age x stage. Cannot proceed. Please make sure that matrix dimensions match stage descriptions.", call. = FALSE)
      }
      age_bit <- c(apply(as.matrix(c(0:(newmult-1))), 1, rep, dim(labels_orig)[1]))
      core_labels <- cbind.data.frame(age_bit, do.call("rbind.data.frame", replicate(newmult, labels_orig, simplify = FALSE)))
      core_labels$agestage_id <- c(1:length(core_labels$stage_id))
      core_labels$agestage <- apply(as.matrix(core_labels$agestage_id), 1, function(X) {
        paste(core_labels$age_bit[X], core_labels$stage[X])
      })
      labels <- core_labels
      names(labels) <- c("age", "stage_id", "stage", "agestage_id", "agestage")
    }
    
    if (!stochastic) {
      modlabels <- cbind.data.frame(as.matrix(rep(c(1:multiplier), each = dim(labels)[1])), 
        do.call("rbind.data.frame", apply(as.matrix(c(1:multiplier)), 1, function(X) return(labels))))
      output <- cbind.data.frame(modlabels, baldrick)
      if (is.element("age", names(labels))) {
        names(output) <- c("matrix", "age", "stage_id", "stage", "agestage_id", "agestage", "ss_prop")
      } else {
        names(output) <- c("matrix", "stage_id", "stage", "ss_prop")
      }
    } else {
      modlabels <- cbind.data.frame(as.matrix(rep(unlist(used_poppatches), each = dim(labels)[1])), 
        do.call("rbind.data.frame", apply(as.matrix(c(1:multiplier)), 1, function(X) return(labels))))
      unipoppatches <- unique(modlabels[,1])
      
      modnums <- apply(as.matrix(modlabels[,1]), 1, function(X) {
        return(which(unipoppatches == X))
      })
      
      output <- cbind.data.frame(modnums, modlabels, baldrick)
      if (is.element("age", names(labels))) {
        names(output) <- c("matrix_set", "poppatch", "age", "stage_id", "stage",
          "agestage_id", "agestage", "ss_prop")
      } else {
        names(output) <- c("matrix_set", "poppatch", "stage_id", "stage",
          "ss_prop")
      }
    }
    rownames(output) <- c(1:dim(output)[1])
  } else {
    labels <- mats$hstages
    
    if (!stochastic) {
      modlabels <- cbind.data.frame(as.matrix(rep(c(1:multiplier), each = dim(labels)[1])), 
        do.call("rbind.data.frame", apply(as.matrix(c(1:multiplier)), 1, function(X) return(labels))))
      
      outputh <- cbind.data.frame(modlabels, baldrick)
      
      names(outputh) <- c("matrix", "stage_id_2", "stage_id_1", "stage_2",
        "stage_1", "ss_prop")
    } else {
      modlabels <- cbind.data.frame(as.matrix(rep(unlist(used_poppatches), each = dim(labels)[1])), 
        do.call("rbind.data.frame", apply(as.matrix(c(1:multiplier)), 1, function(X) return(labels))))
      unipoppatches <- unique(modlabels[,1])
      
      modnums <- apply(as.matrix(modlabels[,1]), 1, function(X) {
        return(which(unipoppatches == X))
      })
      
      outputh <- cbind.data.frame(modnums, modlabels, baldrick)
      names(outputh) <- c("matrix_set", "poppatch", "stage_id_2", "stage_id_1",
        "stage_2", "stage_1", "ss_prop")
    }
    
    rownames(outputh) <- c(1:dim(outputh)[1])
    
    ahlabels <- mats$ahstages[,c("stage_id", "stage")]
    ss2 <- c(apply(as.matrix(c(1:multiplier)), 1, function(X) {
      rightset <- if (!stochastic) {
        subset(outputh, matrix == X)
      } else {
        subset(outputh, matrix_set == X)
      }
      apply(as.matrix(ahlabels[,1]), 1, function(Y) {
        sum(rightset$ss_prop[which(rightset$stage_id_2 == Y)])
      })
    }))
    
    if (!stochastic) {
      outputah <- cbind.data.frame(rep(c(1:multiplier), each = length(ahlabels[,1])), 
        rep(ahlabels[,1], multiplier), rep(ahlabels[,2], multiplier), ss2)
      names(outputah) <- c("matrix", "stage_id", "stage", "ss_prop")
    } else {
      modlabels <- cbind.data.frame(as.matrix(rep(unlist(used_poppatches), each = dim(ahlabels)[1])), 
        do.call("rbind.data.frame", apply(as.matrix(c(1:multiplier)), 1, function(X) return(ahlabels))))
      modnums <- apply(as.matrix(modlabels[,1]), 1, function(X) {
        return(which(unipoppatches == X))
      })
      
      outputah <- cbind.data.frame(modnums, modlabels, ss2)
      names(outputah) <- c("matrix_set", "poppatch", "stage_id", "stage", "ss_prop")
    }
    
    rownames(outputah) <- c(1:dim(outputah)[1])
    
    if (!stochastic) {
      output <- list(hist = outputh, ahist = outputah)
    } else {
      output <- list(hist = outputh, ahist = outputah, projections = princegeorge)
    }
  }
  
  return(output)
}

#' Estimate Stable Stage Distribution of a Single Population Projection Matrix
#' 
#' \code{stablestage3.matrix()} returns the stable stage distribution for a 
#' population projection matrix. This function can handle large and sparse
#' matrices, and so can be used with large historical matrices, IPMs, 
#' age x stage matrices, as well as smaller ahistorical matrices.
#' 
#' @param mats A population projection matrix of class \code{matrix}.
#' @param sparse A text string indicating whether to use sparse matrix encoding
#' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
#' \code{"auto"}.
#' @param ... Other parameters.
#' 
#' @return This function returns the stable stage distribution corresponding to
#' the input matrix.
#' 
#' @seealso \code{\link{stablestage3}()}
#' @seealso \code{\link{stablestage3.lefkoMat}()}
#' 
#' @examples
#' #Lathyrus example
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
#' stablestage3(ehrlen3mean$A[[1]])
#' 
#' # Cypripedium stochastic example
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
#' stablestage3(cypmatrix2r$A[[1]])
#' 
#' @export
stablestage3.matrix <- function(mats, sparse = "auto", ...)
{
  if (is.element(tolower(sparse), c("y", "yes", "yea", "yeah", "t", "true"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(sparse), c("n", "no", "non", "nah", "f", "false"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats)
    dense_elements <- length(which(mats != 0))
    
    if ((dense_elements / elements_total) < 0.5) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  
  wcorr <- .ss3matrix(mats, sparsemethod)
  
  return(wcorr)
}

#' Estimate Reproductive Value
#' 
#' \code{repvalue3()} is a generic function that estimates returns the
#' reproductive values of stages in a population projection matrix or a set of
#' matrices. The specifics of estimation vary with the class of input object.
#' This function is made to handle very large and sparse matrices supplied as
#' \code{lefkoMat} objects or as individual matrices, and can be used with large
#' historical matrices, IPMs, age x stage matrices, as well as ahistorical
#' matrices.
#' 
#' @param mats A lefkoMat object, or population projection matrix.
#' @param ... Other parameters.
#' 
#' @return The value returned depends on the class of the \code{mats} argument.
#' See related functions for details.
#' 
#' @seealso \code{\link{repvalue3.lefkoMat}()}
#' @seealso \code{\link{repvalue3.matrix}()}
#' 
#' @examples
#' # Lathyrus deterministic example
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
#' repvalue3(ehrlen3mean)
#' 
#' # Cypripedium stochastic example
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
#' repvalue3(cypmatrix2r, stochastic = TRUE)
#' 
#' @export
repvalue3 <- function(mats, ...) UseMethod("repvalue3")

#' Estimate Reproductive Value Vectors of Matrices in a lefkoMat Object
#' 
#' \code{repvalue3.lefkoMat()} returns the reproductive values for stages in a
#' set of population projection matrices provided as a \code{lefkoMat} object.
#' This function can handle large and sparse matrices, and so can be used with
#' large historical matrices, IPMs, age x stage matrices, as well as ahistorical
#' matrices.
#' 
#' @param mats An object of class \code{lefkoMat} object.
#' @param stochastic A logical value indicating whether to use deterministic
#' (\code{FALSE}) or stochastic (\code{TRUE}) analysis. Defaults to
#' \code{FALSE}.
#' @param times An integer variable indicating number of occasions to project if
#' using stochastic analysis. Defaults to 10000.
#' @param tweights An optional vector indicating the probability weighting to
#' use for each matrix in stochastic simulations. If not given, then defaults to
#' equal weighting.
#' @param seed A number to use as a random number seed.
#' @param sparse A text string indicating whether to use sparse matrix encoding
#' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
#' \code{"auto"}.
#' @param ... Other parameters.
#' 
#' @return This function returns the asymptotic reproductive value vectors if
#' deterministic analysis is chosen, and long-run mean reproductive value
#' vectors if stochastic analysis is chosen.
#' 
#' The output depends on whether the \code{lefkoMat} object used as input is
#' ahistorical or historical, and whether the analysis is deterministic or
#' stochastic. If deterministic and ahistorical, then a single data frame is
#' output, which includes the number of the matrix within the \code{A} element
#' of the input \code{lefkoMat} object, followed by the stage id (numeric and
#' assigned through \code{\link{sf_create}()}), the stage name, and the
#' estimated proportion of the reproductive value vector (\code{rep_value}). If
#' stochastic and ahistorical, then a single data frame is output starting with
#' the number of the population-patch (\code{matrix_set}), a string
#' concatenating the names of the population and the patch (\code{poppatch}),
#' the assigned stage id number (\code{stage_id}), and the stage name
#' (\code{stage}), and the long-run mean reproductive value vector
#' (\code{rep_value}).
#' 
#' If a historical matrix is used as input, then two data frames are output
#' into a list object. The \code{hist} element describes the historical
#' stage-pair reproductive values, while the \code{ahist} element describes the
#' stage reproductive values. If deterministic, then \code{hist} contains a data
#' frame including the matrix number (\code{matrix}), the numeric stage
#' designations for stages in occasions \emph{t} and \emph{t}-1,
#' (\code{stage_id_2} and \code{stage_id_1}, respectively), followed by the
#' respective stage names (\code{stage_2} and \code{stage_1}), and ending with
#' the estimated reproductive values (\code{rep_value}). The associated
#' \code{ahist} element is as before. If stochastic, then the \code{hist}
#' element contains a single data frame with the number of the population-patch
#' (\code{matrix_set}), a string concatenating the names of the population and
#' the patch (\code{poppatch}), the assigned stage id numbers in times \emph{t}
#' and \emph{t}-1 (\code{stage_id_2} and \code{stage_id_2}, respectively), and
#' the associated stage names (\code{stage_2} and \code{stage_1}, respectively),
#' and the long-run mean reproductive values (\code{rep_value}). The associated
#' \code{ahist} element is as before in the ahistorical, stochastic case.
#'
#' In addition to the data frames noted above, stochastic analysis will result
#' in the additional output of a list of matrices containing the actual
#' projected reproductive value vectors across all projected occasions, in the
#' order of population-patch combinations in the \code{lefkoMat} input.
#'
#' @section Notes:
#' In stochastic analysis, the projected mean reproductive value vector is the
#' arithmetic mean across the final projected 1000 occasions if the simulation
#' is at least 2000 projected occasions long. If between 500 and 2000 projected
#' occasions long, then only the final 200 are used, and if fewer than 500
#' occasions are used, then all are used. Note that because reproductive values
#' in stochastic simulations can change greatly in the initial portion of the
#' run, we encourage a minimum 2000 projected occasions per simulation, with
#' 10000 preferred.
#' 
#' @seealso \code{\link{repvalue3}()}
#' @seealso \code{\link{repvalue3.matrix}()}
#' 
#' @examples
#' # Lathyrus deterministic example
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
#' repvalue3(ehrlen3mean)
#' 
#' # Cypripedium stochastic example
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
#' repvalue3(cypmatrix2r, stochastic = TRUE)
#' 
#' @export
repvalue3.lefkoMat <- function(mats, stochastic = FALSE, times = 10000,
  tweights = NA, seed = NA, sparse = "auto", ...) {
  
  matrix_set <- poppatch <- NULL
  
  if (is.element(tolower(sparse), c("y", "yes", "yea", "yeah", "t", "true"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(sparse), c("n", "no", "non", "nah", "f", "false"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats$A[[1]])
    dense_elements <- length(which(mats$A[[1]] != 0))
    
    if ((dense_elements / elements_total) < 0.5) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  
  if (!stochastic) {
    baldrick <- if (any(class(mats$A) == "matrix")) {
      
      almost_final <- .rv3matrix(mats$A, sparsemethod)
      
    } else if (class(mats$A) == "list") {
      final <- unlist(lapply(mats$A, function(X) {
        almost_final <- .rv3matrix(X, sparsemethod)
        return(almost_final/almost_final[which(almost_final == min(almost_final[which(almost_final > 0)])[1])])
      }))
    } else {
      stop("Input not recognized.")
    }
  } else {
    
    if (!is.na(seed)) {
      set.seed(seed)
    }
    
    mats$labels$poppatch <- paste(mats$labels$pop, mats$labels$patch)
    used_poppatches <- as.list(unique(mats$labels$poppatch))
    
    # Here we get the full stage distribution and reproductive value vector series for all occasions, as a list
    # Stage distributions are the top half the matrix, and reproductive value is at the bottom
    princegeorge <- lapply(used_poppatches, function(X) {
      used_slots <- which(mats$labels$poppatch == X)
      
      if (length(used_slots) < 2) {
        warning("Only 1 annual matrix found for some population-patch combinations. Stochastic analysis requires multiple annual matrices per population-patch combination.", call. = FALSE)
      }
      
      if (!is.na(tweights)) {
        if (length(tweights) != length(used_slots)) {
          if (length(tweights) == length(mats$A)) {
            used_weights <- tweights[used_slots] / sum(tweights[used_slots])
          } else {
            stop("Option tweights must be either NA, or a numeric vector equal to the number of years supplied or matrices supplied.",
              call. = FALSE)
          }
        } else {
          used_weights <- tweights / sum(tweights)
        }
      } else {
        used_weights <- rep(1, length(used_slots))
        used_weights <- used_weights / sum(used_weights)
      }
      
      theprophecy <- sample(used_slots, times, replace = TRUE, prob = used_weights) - 1
      starter <- .ss3matrix(mats$A[[used_slots[1]]], sparsemethod)
      
      theseventhmatrix <- .proj3(starter, mats$A, theprophecy, 1, 0, 0)
      almostall <- theseventhmatrix[((dim(mats$A[[1]])[1]) + 1):(3 * (dim(mats$A[[1]])[1])),]
      
      return(almostall)
    })
    
    # Now we create the mean distributions
    baldrick <- unlist(
      lapply(princegeorge, function(X) {
        if (times > 2000) {
          usedX1 <- X[1:(dim(mats$A[[1]])[1]), (times - 998):(times+1)]
          usedX2 <- X[((dim(mats$A[[1]])[1]) + 1):(2*(dim(mats$A[[1]])[1])),1:1000]
        } else if (times > 500) {
          usedX1 <- X[1:(dim(mats$A[[1]])[1]), (times - 198):(times+1)]
          usedX2 <- X[((dim(mats$A[[1]])[1]) + 1):(2*(dim(mats$A[[1]])[1])), 1:200]
        } else {
          usedX1 <- X[1:(dim(mats$A[[1]])[1]),]
          usedX2 <- X[((dim(mats$A[[1]])[1]) + 1):(2*(dim(mats$A[[1]])[1])),]
        }
        meanX1 <- apply(usedX1, 1, mean)
        meanX2 <- apply(usedX2, 1, mean)
        meanX2 <- zapsmall(meanX2)
        meanX2 <- meanX2 / meanX2[(which(meanX2 > 0)[1])]
        meanX <- c(meanX1, meanX2)
        
        return(meanX)
      })
    )
    
    princegeorge <- lapply(princegeorge, function(X) {
      return(X[((dim(mats$A[[1]])[1]) + 1):(2*(dim(mats$A[[1]])[1])),])
    })
  }
  
  if (class(mats$A) == "list") {
    if (!stochastic) {
      multiplier <- length(mats$A)
    } else {
      multiplier <- length(used_poppatches)
    }
  } else multiplier <- 1
  
  if (all(is.na(mats$hstages))) {
    labels_orig <- mats$ahstages[,1:2]
    mat_dims <- dim(mats$A[[1]])[1]
    if (dim(labels_orig)[1] == mat_dims) {
      labels <- labels_orig
    } else {
      newmult <- mat_dims / dim(labels_orig)[1]
      if (mat_dims %% dim(labels_orig)[1] != 0) {
        stop("Matrices do not appear to be ahistorical, historical, or age x stage. Cannot proceed. Please make sure that matrix dimensions match stage descriptions.", call. = FALSE)
      }
      age_bit <- c(apply(as.matrix(c(0:(newmult-1))), 1, rep, dim(labels_orig)[1]))
      core_labels <- cbind.data.frame(age_bit, do.call("rbind.data.frame", replicate(newmult, labels_orig, simplify = FALSE)))
      core_labels$agestage_id <- c(1:length(core_labels$stage_id))
      core_labels$agestage <- apply(as.matrix(core_labels$agestage_id), 1, function(X) {
        paste(core_labels$age_bit[X], core_labels$stage[X])
      })
      labels <- core_labels
      names(labels) <- c("age", "stage_id", "stage", "agestage_id", "agestage")
    }
    
    if (!stochastic) {
      modlabels <- cbind.data.frame(as.matrix(rep(c(1:multiplier), each = dim(labels)[1])), 
        do.call("rbind.data.frame", apply(as.matrix(c(1:multiplier)), 1, function(X) return(labels))))
      output <- cbind.data.frame(modlabels, baldrick)
      if (is.element("age", names(labels))) {
        names(output) <- c("matrix", "age", "stage_id", "stage", "agestage_id", "agestage", "rep_value")
      } else {
        names(output) <- c("matrix", "stage_id", "stage", "rep_value")
      }
    } else {
      modlabels <- cbind.data.frame(as.matrix(rep(unlist(used_poppatches), each = dim(labels)[1])), 
        do.call("rbind.data.frame", apply(as.matrix(c(1:multiplier)), 1, function(X) return(labels))))
      unipoppatches <- unique(modlabels[,1])
      
      modnums <- apply(as.matrix(modlabels[,1]), 1, function(X) {
        return(which(unipoppatches == X))
      })
      output <- cbind.data.frame(modnums, modlabels, baldrick[(mat_dims+1):(2*mat_dims)])
      if (is.element("age", names(labels))) {
        names(output) <- c("matrix_set", "poppatch", "age", "stage_id", "stage", "agestage_id", "agestage", "rep_value")
      } else {
        names(output) <- c("matrix_set", "poppatch", "stage_id", "stage", "rep_value")
      }
    }
    
    rownames(output) <- c(1:dim(output)[1])
    
  } else {
    
    # This section translates historical results to ahistorical and then cleans everything up
    if (!stochastic) {
      ss3 <- stablestage3.lefkoMat(mats, sparse = sparse) #stablestage3.lefkoMat
      rahist <- ss3$ahist
      rhist <-ss3$hist
      rhist$ss3sum <- apply(as.matrix(c(1:dim(rhist)[1])), 1, function(X) {
        rahist$ss_prop[intersect(which(rahist$stage_id == rhist$stage_id_2[X]), 
          which(rahist$matrix == rhist$matrix[X]))]
      })
      rhist$sscorr <- rhist$ss_prop / rhist$ss3sum
      rhist$sscorr[which(is.na(rhist$sscorr))] <- 0
      rhist$rep_value <- baldrick
      
      rhist$rv3raw <- rhist$sscorr * rhist$rep_value
      
      outputh <- rhist[,c("matrix", "stage_id_2", "stage_id_1", "stage_2", "stage_1", "rep_value")]
      
      ahlabels <- mats$ahstages[,c("stage_id", "stage")]
      rv2 <- Re(c(apply(as.matrix(c(1:multiplier)), 1, function(X) {
        rightset <- subset(rhist, matrix == X)
        apply(as.matrix(ahlabels[,1]), 1, function(Y) {
          sum(rightset$rv3raw[which(rightset$stage_id_2 == Y)])
        })
      })))
      outputah <- cbind.data.frame(rep(c(1:multiplier), each = length(ahlabels[,1])),
        rep(ahlabels[,1], multiplier), rep(ahlabels[,2], multiplier), rv2)
      names(outputah) <- c("matrix", "stage_id", "stage", "rep_value_unc")
      
      outputah$rep_value <- apply(as.matrix(c(1:dim(outputah)[1])), 1, function(X) {
        matsub <- subset(outputah, matrix == outputah$matrix[X])
        entrystage <- min(which(abs(matsub$rep_value_unc) > 0))
        return(outputah$rep_value_unc[X] / matsub$rep_value_unc[entrystage])
      })
      
      outputah <- outputah[,c("matrix", "stage_id", "stage", "rep_value")]
      rownames(outputah) <- c(1:dim(outputah)[1])
      
    } else {
      
      ss_unlisted <- apply(as.matrix(c(1:multiplier)), 1, function(X) {
        return(baldrick[((2 * (X - 1) * (dim(mats$A[[1]])[1])) + 1): (((2 * (X - 1)) + 1) * (dim(mats$A[[1]])[1]))])
      })
      rv_unlisted <- apply(as.matrix(c(1:multiplier)), 1, function(X) {
        return(baldrick[((((2 * (X - 1)) + 1) * (dim(mats$A[[1]])[1])) + 1): (X * 2 * (dim(mats$A[[1]])[1]))])
      })
      
      ss_sums <- apply(as.matrix(c(1:multiplier)), 1, function(X) {
        sum_vecs <- apply(as.matrix(mats$hstages$stage_id_2), 1, function(Y) {
          sum(ss_unlisted[(which(mats$hstages$stage_id_2 == Y)) + ((X - 1) * dim(mats$hstages)[1])])
        })
        return(sum_vecs)
      })
      
      ahistsize <- length(mats$ahstages$stage_id)
      histsize <- length(mats$hstages$stage_id_2)
      
      poppatchcol <- rep(unlist(used_poppatches), each = ahistsize)
      unipoppatches <- unique(poppatchcol)
      poppatchnum <- apply(as.matrix(poppatchcol), 1, function(X) {
        return(which(unipoppatches == X))
      })
      
      rahist <- cbind.data.frame(matrix_set = poppatchnum, poppatch = poppatchcol, 
        stage_id = rep(mats$ahstages$stage_id, times = multiplier),
        stage = rep(mats$ahstages$stage, times = multiplier))
      
      poppatchcol <- rep(unlist(used_poppatches), each = histsize)
      poppatchnum <- apply(as.matrix(poppatchcol), 1, function(X) {
        return(which(unipoppatches == X))
      })
      
      rhist <- cbind.data.frame(matrix_set = poppatchnum, poppatch = poppatchcol,
        stage_id_2 = rep(mats$hstages$stage_id_2, times = multiplier),
        stage_id_1 = rep(mats$hstages$stage_id_1, times = multiplier),
        stage_2 = rep(mats$hstages$stage_2, times = multiplier),
        stage_1 = rep(mats$hstages$stage_1, times = multiplier))
      
      rhist$ss_prop <- as.vector(ss_unlisted)
      rhist$ss3sum <- as.vector(ss_sums)
      rhist$sscorr <- rhist$ss_prop / rhist$ss3sum
      rhist$sscorr[which(is.na(rhist$sscorr))] <- 0
      rhist$rep_value <- as.vector(rv_unlisted)
      
      rhist$rv3raw <- rhist$sscorr * rhist$rep_value
      
      outputh <- rhist[,c("matrix_set", "poppatch", "stage_id_2", "stage_id_1",
        "stage_2", "stage_1", "rep_value")]
      
      ahlabels <- mats$ahstages[,c("stage_id", "stage")]
      rv2 <- Re(c(apply(as.matrix(unlist(used_poppatches)), 1, function(X) {
        rightset <- subset(rhist, poppatch == X)
        apply(as.matrix(ahlabels[,1]), 1, function(Y) {
          sum(rightset$rv3raw[which(rightset$stage_id_2 == Y)])
        })
      })))
      
      poppatchcol <- rep(unlist(used_poppatches), each = ahistsize)
      
      outputah <- cbind.data.frame(rep(c(1:multiplier), each = length(ahlabels[,1])),
        poppatchcol, rep(ahlabels[,1], multiplier), rep(ahlabels[,2], multiplier), rv2)
      names(outputah) <- c("matrix_set", "poppatch", "stage_id", "stage",
        "rep_value_unc")
      
      outputah$rep_value <- apply(as.matrix(c(1:dim(outputah)[1])), 1, function(X) {
        matsub <- subset(outputah, poppatch == outputah$poppatch[X])
        
        if (!all(matsub$rep_value_unc == 0)) {
          entrystage <- min(which(abs(matsub$rep_value_unc) > 0))
          return(outputah$rep_value_unc[X] / matsub$rep_value_unc[entrystage])
        } else return(0)
      })
      
      outputah <- outputah[,c("matrix_set", "poppatch", "stage_id", "stage", "rep_value")]
      rownames(outputah) <- c(1:dim(outputah)[1])
    }
    
    if (!stochastic) {
      output <-list(hist = outputh, ahist = outputah)
    } else {
      output <-list(hist = outputh, ahist = outputah, projections = princegeorge)
    }
  }
  return(output)
}

#' Estimate Reproductive Value Vector for a Single Population Projection Matrix
#' 
#' \code{repvalue3.matrix()} returns the reproductive values for stages in a 
#' population projection matrix. The function makes no assumptions about whether
#' the matrix is ahistorical and simply provides standard reproductive values
#' corresponding to each row, meaning that the overall reproductive values of
#' basic life history stages in a historical matrix are not provided (the 
#' \code{\link{repvalue3.lefkoMat}()} function estimates these on the basis of
#' stage description information provided in the \code{lefkoMat} object used as
#' input in that function). This function can handle large and sparse matrices,
#' and so can be used with large historical matrices, IPMs, age x stage
#' matrices, as well as smaller ahistorical matrices.
#' 
#' @param mats A population projection matrix.
#' @param sparse A text string indicating whether to use sparse matrix encoding
#' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
#' \code{"auto"}.
#' @param ... Other parameters.
#' 
#' @return This function returns a vector data frame characterizing the 
#' reproductive values for stages of a population projection matrix. This is 
#' given as the left eigenvector associated with largest real part of the
#' dominant eigenvalue, divided by the first non-zero element of the left 
#' eigenvector.
#' 
#' @seealso \code{\link{repvalue3}()}
#' @seealso \code{\link{repvalue3.lefkoMat}()}
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
#' repvalue3(ehrlen3mean$A[[1]])
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
#' repvalue3(cypmatrix2r$A[[1]])
#' 
#' @export
repvalue3.matrix <- function(mats, sparse = "auto", ...)
{
  if (is.element(tolower(sparse), c("y", "yes", "yea", "yeah", "t", "true"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(sparse), c("n", "no", "non", "nah", "f", "false"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats)
    dense_elements <- length(which(mats != 0))
    
    if ((dense_elements / elements_total) < 0.5) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  
  v <- .rv3matrix(mats, sparsemethod)

  return(v)
}

#' Estimate Sensitivity of Population Growth Rate to Matrix Elements
#' 
#' \code{sensitivity3()} is a generic function that returns the sensitivity of
#' the population growth rate to the elements of the matrices in a matrix
#' population model. Currently, this function estimates both deterministic and
#' stochastic sensitivities, where the growth rate is \eqn{\lambda} in the
#' former case and the log of the stochastic \eqn{\lambda} in the latter case.
#' This function is made to handle very large and sparse matrices supplied as
#' \code{lefkoMat} objects, as lists of matrices, and as individual matrices.
#' 
#' @param mats A lefkoMat object, or population projection matrix, for which
#' the stable stage distribution is desired.
#' @param ... Other parameters
#' 
#' @return The value returned depends on the class of the \code{mats} argument.
#' 
#' @seealso \code{\link{sensitivity3.lefkoMat}()}
#' @seealso \code{\link{sensitivity3.matrix}()}
#' @seealso \code{\link{sensitivity3.list}()}
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
#' sensitivity3(ehrlen3mean)
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
#' sensitivity3(cypmatrix2r)
#' 
#' @export
sensitivity3 <- function(mats, ...) UseMethod("sensitivity3")

#' Estimate Sensitivity of Population Growth Rate of a lefkoMat Object
#' 
#' \code{sensitivity3.lefkoMat()} returns the sensitivities of population growth
#' rate to elements of all \code{$A} matrices in an object of class
#' \code{lefkoMat}. If deterministic, then \eqn{\lambda} is taken as the
#' population growth rate. If stochastic, then the log of stochastic
#' \eqn{\lambda}, or the log stochastic growth rate, is taken as the population
#' growth rate. This function can handle large and sparse matrices, and so can
#' be used with large historical matrices, IPMs, age x stage matrices, as well
#' as smaller ahistorical matrices.
#' 
#' @param mats An object of class \code{lefkoMat}.
#' @param stochastic A logical value determining whether to conduct a
#' deterministic (FALSE) or stochastic (TRUE) sensitivity analysis. Defaults to
#' FALSE.
#' @param steps The number of occasions to project forward in stochastic
#' simulation. Defaults to 10,000.
#' @param time_weights Numeric vector denoting the probabilistic weightings of
#' annual matrices. Defaults to equal weighting among occasions.
#' @param sparse A text string indicating whether to use sparse matrix encoding
#' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
#' \code{"auto"}.
#' @param append_mats A logical value indicating whether to include the original
#' A, U, and F matrices in the output \code{lefkoSens} object.
#' @param ... Other parameters.
#' 
#' @return This function returns an object of class \code{lefkoSens}, which is a
#' list of 8 elements. The first, \code{h_sensmats}, is a list of historical
#' sensitivity matrices (NULL if an ahMPM is used as input). The second,
#' \code{ah_elasmats}, is a list of either ahistorical sensitivity matrices if
#' an ahMPM is used as input, or, if an hMPM is used as input, then the result
#' is a list of ahistorical matrices based on the equivalent historical
#' dependencies assumed in the input historical matrices. The third element,
#' \code{h_stages}, is a data frame showing historical stage pairs (NULL if
#' ahMPM used as input). The fourth element, \code{agestages}, show the order of
#' age-stage combinations, if age-by-stage MPMs have been supplied. The fifth
#' element, \code{ah_stages}, is a data frame showing the order of ahistorical
#' stages. The last 3 elements are the A, U, and F portions of the input.
#' 
#' @section Notes:
#' Deterministic sensitivities are estimated as eqn. 9.14 in Caswell (2001,
#' Matrix Population Models). Stochastic sensitivities are estimated as eqn.
#' 14.97 in Caswell (2001). Note that stochastic sensitivities are of the log of
#' the stochastic \eqn{\lambda}.
#'
#' @seealso \code{\link{sensitivity3}()}
#' @seealso \code{\link{sensitivity3.matrix}()}
#' @seealso \code{\link{sensitivity3.list}()}
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
#' sensitivity3(ehrlen3, stochastic = TRUE)
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
#' sensitivity3(cypmatrix2r)
#' 
#' @export
sensitivity3.lefkoMat <- function(mats, stochastic = FALSE, steps = 10000,
  time_weights = NA, sparse = "auto", append_mats = FALSE, ...) {
  
  if (is.element(tolower(sparse), c("y", "yes", "yea", "yeah", "t", "true"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(sparse), c("n", "no", "non", "nah", "f", "false"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats$A[[1]])
    dense_elements <- length(which(mats$A[[1]] != 0))
    
    if ((dense_elements / elements_total) < 0.5) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  
  if (stochastic & length(mats$A) < 2) {
    stop("Stochastic sensitivity estimation cannot be completed with fewer than 2 annual matrices.", call. = FALSE)
  }
  
  if (!is.element("agestages", names(mats))) {
    mats$agestages <- NA
  }
  
  if (!stochastic) {
    # Deterministic sensitivity analysis
    
    baldrick <- if (any(class(mats$A) == "matrix")) {
      .sens3matrix(mats$A, sparsemethod)
      
    } else if (class(mats$A) == "list") {
      if (all(is.na(mats$hstages))) {
        lapply(mats$A, .sens3matrix, sparsemethod)
      } else {
        lapply(mats$A, .sens3hlefko, mats$ahstages, mats$hstages)
      }
    } else {
      stop("Input not recognized.")
    }
    
    if (all(is.na(mats$hstages))) {
      
      ahlabels <- mats$ahstages
      
      output <- list(h_sensmats = NULL, ah_sensmats = baldrick, h_stages = NULL, 
        ah_stages = ahlabels)
      
    } else {
      
      he_list <- lapply(baldrick, function(X) {X$h_smat})
      ahe_list <- lapply(baldrick, function(X) {X$ah_smat})
      
      hlabels <- mats$hstages
      ahlabels <- mats$ahstages
      
      output <- list(h_sensmats = he_list, ah_sensmats = ahe_list,
        h_stages = hlabels, ah_stages = ahlabels)
    }
  } else {
    # Stochastic sensitivity analysis
    
    if(!any(is.na(time_weights))) {
      returned_cubes <- .stoch_senselas(mats, times = steps, style = 1,
        tweights = time_weights) 
    } else {
      returned_cubes <- .stoch_senselas(mats, times = steps, style = 1) 
    }
    
    main_cube <- returned_cubes[[1]]
    ah_cube <- returned_cubes[[2]]
    
    returned_main <- lapply(as.list(c(1:dim(main_cube)[3])), function(X) {
        return(main_cube[,,X])
      }
    )
    
    if (!all(is.na(mats$hstages))) {
      returned_ah <- lapply(as.list(c(1:dim(ah_cube)[3])), function(X) {
          return(ah_cube[,,X])
        }
      )
      output <- list(h_sensmats = returned_main, ah_sensmats = returned_ah,
        h_stages = mats$hstages, agestages = mats$agestages,
        ah_stages = mats$ahstages)
    } else {
      output <- list(h_sensmats = NULL, ah_sensmats = returned_main,
        h_stages = mats$hstages, agestages = mats$agestages,
        ah_stages = mats$ahstages)
    }
  }
  
  if (append_mats) {
    output$A <- mats$A
    output$U <- mats$U
    output$F <- mats$F
  }
  
  class(output) <- "lefkoSens"
  
  return(output)
}

#' Estimate Sensitivity of Population Growth Rate of a Single Matrix
#' 
#' \code{sensitivity3.matrix()} returns the sensitivities of \eqn{\lambda} to
#' elements of a single matrix. Because this handles only one matrix, the
#' sensitivities are inherently deterministic and based on the dominant eigen
#' value as the best metric of the population growth rate. This function can
#' handle large and sparse matrices, and so can be used with large historical
#' matrices, IPMs, age x stage matrices, as well as smaller ahistorical
#' matrices.
#' 
#' @param mats An object of class \code{matrix}.
#' @param sparse A text string indicating whether to use sparse matrix encoding
#' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
#' \code{"auto"}.
#' @param ... Other parameters.
#' 
#' @return This function returns a single deterministic sensitivity matrix.
#' 
#' @seealso \code{\link{sensitivity3}()}
#' @seealso \code{\link{sensitivity3.lefkoMat}()}
#' @seealso \code{\link{sensitivity3.list}()}
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
#' sensitivity3(ehrlen3mean$A[[1]])
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
#' sensitivity3(cypmatrix2r$A[[1]])
#' 
#' @export
sensitivity3.matrix <- function(mats, sparse = "auto", ...)
{
  if (is.element(tolower(sparse), c("y", "yes", "yea", "yeah", "t", "true"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(sparse), c("n", "no", "non", "nah", "f", "false"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats)
    dense_elements <- length(which(mats != 0))
    
    if ((dense_elements / elements_total) < 0.5) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  wcorr <- .sens3matrix(mats, sparsemethod)
  
  return(wcorr)
}

#' Estimate Sensitivity of Population Growth Rate of a List of Matrices
#' 
#' \code{sensitivity3.list()} returns the sensitivities of population growth
#' rate to elements of matrices supplied in a list. The sensitivity analysis can
#' be deterministic or stochastic, but if the latter then at least two A
#' matrices must be included in the list. This function can handle large and
#' sparse matrices, and so can be used with large historical matrices, IPMs,
#' age x stage matrices, as well as smaller ahistorical matrices.
#' 
#' @param mats An object of class \code{matrix}.
#' @param stochastic A logical value determining whether to conduct a
#' deterministic (FALSE) or stochastic (TRUE) sensitivity analysis. Defaults to
#' FALSE.
#' @param steps The number of occasions to project forward in stochastic
#' simulation. Defaults to 10,000.
#' @param time_weights Numeric vector denoting the probabilistic weightings of
#' annual matrices. Defaults to equal weighting among occasions.
#' @param historical A logical value indicating whether matrices are historical.
#' Defaults to FALSE.
#' @param sparse A text string indicating whether to use sparse matrix encoding
#' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
#' \code{"auto"}.
#' @param append_mats A logical value indicating whether to include the original
#' matrices input as object \code{mats} in the output \code{lefkoSense} object.
#' Defaults to FALSE.
#' @param ... Other parameters.
#' 
#' @return This function returns an object of class \code{lefkoSens}, which is a
#' list of 8 elements. The first, \code{h_sensmats}, is a list of historical
#' sensitivity matrices (NULL if an ahMPM is used as input). The second,
#' \code{ah_elasmats}, is a list of ahistorical sensitivity matrices if an ahMPM
#' is used as input (NULL if an hMPM is used as input). The third element,
#' \code{h_stages}, the fourth element, \code{agestages}, and the fifth element,
#' \code{ah_stages}, are NULL. The last 3 elements include the original A
#' matrices supplied (as the \code{A} element), followed by NULLs for the U and
#' F elements.
#' 
#' @section Notes:
#' Deterministic sensitivities are estimated as eqn. 9.14 in Caswell (2001,
#' Matrix Population Models). Stochastic sensitivities are estimated as eqn.
#' 14.97 in Caswell (2001). Note that stochastic sensitivities are with regard
#' to the log of the stochastic \eqn{\lambda}.
#'
#' Currently, this function does not estimate equivalent ahistorical stochastic
#' sensitivities for input historical matrices, due to the lack of guidance
#' input on the order of stages (such guidance is provided within
#' \code{lefkoMat} objects).
#'
#' @seealso \code{\link{sensitivity3}()}
#' @seealso \code{\link{sensitivity3.lefkoMat}()}
#' @seealso \code{\link{sensitivity3.matrix}()}
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
#' sensitivity3(ehrlen3$A)
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
#' sensitivity3(cypmatrix2r$A)
#' 
#' @export
sensitivity3.list <- function(mats, stochastic = FALSE, steps = 10000,
  time_weights = NA, historical = FALSE, sparse = "auto", append_mats = FALSE,
  ...) {
  
  if (is.element(tolower(sparse), c("y", "yes", "yea", "yeah", "t", "true"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(sparse), c("n", "no", "non", "nah", "f", "false"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats[[1]])
    dense_elements <- length(which(mats[[1]] != 0))
    
    if ((dense_elements / elements_total) < 0.5) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  
  if(length(setdiff(unlist(lapply(mats, class)), c("matrix", "array"))) > 0) {
    stop("Input list must be composed only of numeric matrices.", call. = FALSE)
  }
  if (stochastic & length(mats) < 2) {
    stop("Stochastic sensitivity estimation cannot be completed with fewer than 2 annual matrices.", call. = FALSE)
  }
  
  if (!stochastic) {
    # Deterministic sensitivity analysis
    baldrick <- lapply(mats, .sens3matrix, sparsemethod)
    
    if (historical) {
      output <- list(h_sensmats = baldrick, ah_sensmats = NULL, h_stages = NULL,
        agestages = NULL, ah_stages = NULL, A = mats, U = NULL, F = NULL)
    } else {
      output <- list(h_sensmats = NULL, ah_sensmats = baldrick, h_stages = NULL,
        agestages = NULL, ah_stages = NULL, A = mats, U = NULL, F = NULL)
    }
    
  } else {
    # Stochastic sensitivity analysis
    
    if(!any(is.na(time_weights))) {
      returned_cube <- .stoch_senselas(mats, times = steps, style = 1,
        tweights = time_weights)[[1]]
    } else {
      returned_cube <- .stoch_senselas(mats, times = steps, style = 1)[[1]]
    }
    
    returned_list <- lapply(as.list(c(1:dim(returned_cube)[3])), function(X) {
        return(returned_cube[,,X])
      }
    )
    
    if (historical) {
      output <- list(h_sensmats = returned_list, ah_sensmats = NULL,
        h_stages = NULL, agestages = NULL, ah_stages = NULL)
    } else {
      output <- list(h_sensmats = NULL, ah_sensmats = returned_list,
        h_stages = NULL, agestages = NULL, ah_stages = NULL)
    }
  }
  
  if (append_mats) {
    output$A <- mats
  }
  
  class(output) <- "lefkoSens"
  
  return(output)
}

#' Estimate Elasticity of Population Growth Rate to Matrix Elements
#' 
#' \code{elasticity3()} is a generic function that returns the elasticity of
#' the population growth rate to the elements of the matrices in a matrix
#' population model. Currently, this function estimates both deterministic and
#' stochastic elasticities, where the growth rate is \eqn{\lambda} in the former
#' case and the log of the stochastic \eqn{\lambda} in the latter case. This
#' function is made to handle very large and sparse matrices supplied as
#' \code{lefkoMat} objects, as lists of matrices, and as individual matrices.
#' 
#' @param mats A lefkoMat object, or population projection matrix, for which
#' the stable stage distribution is desired.
#' @param ... Other parameters.
#' 
#' @return The value returned depends on the class of the \code{mats} argument.
#' 
#' @seealso \code{\link{elasticity3.lefkoMat}()}
#' @seealso \code{\link{elasticity3.matrix}()}
#' @seealso \code{\link{elasticity3.list}()}
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
#' elasticity3(ehrlen3mean)
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
#' elasticity3(cypmatrix2r)
#' 
#' @export
elasticity3 <- function(mats, ...) UseMethod("elasticity3")

#' Estimate Elasticity of Population Growth Rate of a lefkoMat Object
#' 
#' \code{elasticity3.lefkoMat()} returns the elasticities of population growth
#' rate to elements of all \code{$A} matrices in an object of class
#' \code{lefkoMat}. If deterministic, then \eqn{\lambda} is taken as the
#' population growth rate. If stochastic, then stochastic \eqn{\lambda}, or
#' the stochastic growth rate, is taken as the population growth rate. This
#' function can handle large and sparse matrices, and so can be used with large
#' historical matrices, IPMs, age x stage matrices, as well as smaller
#' ahistorical matrices.
#' 
#' @param mats An object of class \code{lefkoMat}.
#' @param stochastic A logical value determining whether to conduct a
#' deterministic (FALSE) or stochastic (TRUE) elasticity analysis. Defaults to
#' FALSE.
#' @param steps The number of occasions to project forward in stochastic
#' simulation. Defaults to 10,000.
#' @param time_weights Numeric vector denoting the probabilistic weightings of
#' annual matrices. Defaults to equal weighting among occasions.
#' @param sparse A text string indicating whether to use sparse matrix encoding
#' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
#' \code{"auto"}.
#' @param append_mats A logical value indicating whether to include the original
#' A, U, and F matrices in the output \code{lefkoElas} object.
#' @param ... Other parameters.
#' 
#' @return This function returns an object of class \code{lefkoElas}, which is a
#' list with 8 elements. The first, \code{h_elasmats}, is a list of historical
#' elasticity matrices (NULL if an ahMPM is used as input). The second,
#' \code{ah_elasmats}, is a list of either ahistorical elasticity matrices if an
#' ahMPM is used as input, or, if an hMPM is used as input, then the result is a
#' list of elasticity matrices in which historical elasticities have been summed
#' by the stage in occasions \emph{t} and \emph{t}+1 to produce
#' historically-corrected elasticity matrices, which are equivalent in dimension
#' to ahistorical elasticity matrices but reflect the effects of stage in
#' occasion \emph{t}-1. The third element, \code{h_stages}, is a data frame
#' showing historical stage pairs (NULL if ahMPM used as input). The fourth
#' element, \code{agestages}, shows age-stage combinations in the order used in
#' age-by-stage MPMs, if suppled. The fifth element, \code{ah_stages}, is a data
#' frame showing the order of ahistorical stages. The last 3 elements are the A,
#' U, and F portions of the input.
#' 
#' @section Notes:
#' Deterministic elasticities are estimated as eqn. 9.72 in Caswell (2001,
#' Matrix Population Models). Stochastic elasticities are estimated as eqn.
#' 14.99 in Caswell (2001). Note that stochastic elasticities are of the
#' stochastic \eqn{\lambda}, while stochastic sensitivities are with regard to
#' the log of the stochastic \eqn{\lambda}.
#' 
#' @seealso \code{\link{elasticity3}()}
#' @seealso \code{\link{elasticity3.matrix}()}
#' @seealso \code{\link{elasticity3.list}()}
#' @seealso \code{\link{summary.lefkoElas}()}
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
#' elasticity3(ehrlen3, stochastic = TRUE)
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
#' elasticity3(cypmatrix2r)
#' 
#' @export
elasticity3.lefkoMat <- function(mats, stochastic = FALSE, steps = 10000,
  time_weights = NA, sparse = "auto", append_mats = FALSE, ...) {
  
  if (is.element(tolower(sparse), c("y", "yes", "yea", "yeah", "t", "true"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(sparse), c("n", "no", "non", "nah", "f", "false"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats$A[[1]])
    dense_elements <- length(which(mats$A[[1]] != 0))
    
    if ((dense_elements / elements_total) < 0.5) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  
  if (!is.element("agestages", names(mats))) {
    mats$agestages <- NULL
  }
  
  if (!stochastic) {
    # Deterministic elasticity analysis
    
    baldrick <- if (any(class(mats$A) == "matrix")) {
      .elas3matrix(mats$A, sparsemethod)
    } else if (class(mats$A) == "list") {
      
      if (all(is.na(mats$hstages))) {
        lapply(mats$A, .elas3matrix, sparsemethod)
      } else {
        lapply(mats$A, .elas3hlefko, mats$ahstages, mats$hstages)
      }
    } else {
      stop("Input not recognized.")
    }
    
    if (class(mats$A) == "list") {
      multiplier <- length(mats$A)
    } else multiplier <- 1
    
    if (all(is.na(mats$hstages))) {
      
      ahlabels <- mats$ahstages #Originally only the first two columns
      
      output <- list(h_elasmats = NULL, ah_elasmats = baldrick, h_stages = NULL,
        agestages = mats$agestages, ah_stages = ahlabels)
    } else {
      
      he_list <- lapply(baldrick, function(X) {X$h_emat})
      ahe_list <- lapply(baldrick, function(X) {X$ah_emat})
      
      hlabels <- mats$hstages
      ahlabels <- mats$ahstages #Originally only the first two columns
      
      output <- list(h_elasmats = he_list, ah_elasmats = ahe_list, 
        h_stages = hlabels, agestages = mats$agestages, ah_stages = ahlabels)
    }
  } else {
    # Stochastic elasticity analysis
    if(!any(is.na(time_weights))) {
      returned_cubes <- .stoch_senselas(mats, times = steps, style = 2,
        tweights = time_weights) 
    } else {
      returned_cubes <- .stoch_senselas(mats, times = steps, style = 2) 
    }
    
    main_cube <- returned_cubes[[1]]
    ah_cube <- returned_cubes[[2]]
    
    returned_list <- lapply(as.list(c(1:dim(main_cube)[3])), function(X) {
      return(main_cube[,,X])
      }
    )
    
    if (!all(is.na(mats$hstages))) {
      ah_list <- lapply(as.list(c(1:dim(ah_cube)[3])), function(X) {
      return(ah_cube[,,X])
      }
    )
      
      output <- list(h_elasmats = returned_list, ah_elasmats = ah_list,
        h_stages = mats$hstages, agestages = mats$agestages, 
        ah_stages = mats$ahstages)
    } else {
      output <- list(h_elasmats = NULL, ah_elasmats = returned_list,
        h_stages = mats$hstages, agestages = mats$agestages, 
        ah_stages = mats$ahstages)
    }
  }
  
  if (append_mats) {
    output$A <- mats$A
    output$U <- mats$U
    output$F <- mats$F
  }
  
  class(output) <- "lefkoElas"
  
  return(output)
}

#' Estimate Elasticity of Population Growth Rate of a Single Matrix
#' 
#' \code{elasticity3.matrix()} returns the elasticities of lambda to elements
#' of a single matrix.  Because this handles only one matrix, the elasticities
#' are inherently deterministic and based on the dominant eigen value as the
#' best metric of the population growth rate. This function can handle large and
#' sparse matrices, and so can be used with large historical matrices, IPMs,
#' age x stage matrices, as well as smaller ahistorical matrices.
#' 
#' @param mats An object of class \code{matrix}.
#' @param sparse A text string indicating whether to use sparse matrix encoding
#' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
#' \code{"auto"}.
#' @param ... Other parameters.
#' 
#' @return This function returns a single elasticity matrix.
#' 
#' @seealso \code{\link{elasticity3}()}
#' @seealso \code{\link{elasticity3.lefkoMat}()}
#' @seealso \code{\link{elasticity3.list}()}
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
#' elasticity3(ehrlen3mean$A[[1]])
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
#' elasticity3(cypmatrix2r$A[[1]])
#' 
#' @export
elasticity3.matrix <- function(mats, sparse = "auto", ...)
{
  if (is.element(tolower(sparse), c("y", "yes", "yea", "yeah", "t", "true"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(sparse), c("n", "no", "non", "nah", "f", "false"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats)
    dense_elements <- length(which(mats != 0))
    
    if ((dense_elements / elements_total) < 0.5) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  wcorr <- .elas3matrix(mats, sparsemethod)
  
  return(wcorr)
}

#' Estimate Elasticity of Population Growth Rate of a List of Matrices
#' 
#' \code{elasticity3.list()} returns the elasticities of lambda to elements
#' of a single matrix. This function can handle large and sparse matrices, and 
#' so can be used with large historical matrices, IPMs, age x stage matrices,
#' as well as smaller ahistorical matrices.
#' 
#' @param mats A list of objects of class \code{matrix}.
#' @param stochastic A logical value determining whether to conduct a
#' deterministic (FALSE) or stochastic (TRUE) elasticity analysis. Defaults to
#' FALSE.
#' @param steps The number of occasions to project forward in stochastic
#' simulation. Defaults to 10,000.
#' @param time_weights Numeric vector denoting the probabilistic weightings of
#' annual matrices. Defaults to equal weighting among occasions.
#' @param historical A logical value denoting whether the input matrices are
#' historical. Defaults to FALSE.
#' @param sparse A text string indicating whether to use sparse matrix encoding
#' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
#' \code{"auto"}.
#' @param append_mats A logical value indicating whether to include the original
#' matrices input as object \code{mats} in the output \code{lefkoElas} object.
#' @param ... Other parameters.
#' 
#' @return This function returns an object of class \code{lefkoElas}, which is a
#' list with 8 elements. The first, \code{h_elasmats}, is a list of historical
#' elasticity matrices, though in the standard list case it returns a NULL
#' value. The second, \code{ah_elasmats}, is a list of ahistorical elasticity
#' matrices. The third element, \code{h_stages}, the fourth element,
#' \code{agestages}, and the fifth element, \code{ah_stages}, are set to NULL.
#' The last 3 elements are the original A matrices in element A, followed by
#' NULL values for the U and F elements.
#' 
#' @section Notes:
#' Deterministic elasticities are estimated as eqn. 9.72 in Caswell (2001,
#' Matrix Population Models). Stochastic elasticities are estimated as eqn.
#' 14.99 in Caswell (2001). Note that stochastic elasticities are of stochastic
#' \eqn{\lambda}, while stochastic sensitivities are with regard to the log of
#' the stochastic \eqn{\lambda}.
#' 
#' @seealso \code{\link{elasticity3}()}
#' @seealso \code{\link{elasticity3.lefkoMat}()}
#' @seealso \code{\link{elasticity3.matrix}()}
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
#' elasticity3(ehrlen3$A, stochastic = TRUE)
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
#' elasticity3(cypmatrix2r$A)
#' 
#' @export
elasticity3.list <- function(mats, stochastic = FALSE, steps = 10000,
  time_weights = NA, historical = FALSE, sparse = "auto", append_mats = FALSE,
  ...) {
  
  if (is.element(tolower(sparse), c("y", "yes", "yea", "yeah", "t", "true"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(sparse), c("n", "no", "non", "nah", "f", "false"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats[[1]])
    dense_elements <- length(which(mats[[1]] != 0))
    
    if ((dense_elements / elements_total) < 0.5) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  
  if(length(setdiff(unlist(lapply(mats, class)), c("matrix", "array"))) > 0) {
    stop("Input list must be composed only of numeric matrices.", call. = FALSE)
  }
  if (stochastic & length(mats) < 2) {
    stop("Stochastic elasticity estimation cannot be completed with fewer than 2 annual matrices.", call. = FALSE)
  }
  
  if (!stochastic) {
    # Deterministic elasticity analysis
    
    allelems <- length(mats[[1]])
    allnzs <- length(which(mats[[1]] > 0))
    nzprop <- allnzs / allelems
    
    baldrick <- lapply(mats, .elas3matrix, sparsemethod)
    
    if (historical) {
      output <- list(h_elasmats = baldrick, ah_elasmats = NULL, h_stages = NULL,
        ah_stages = NULL)
    } else {
      output <- list(h_elasmats = NULL, ah_elasmats = baldrick, h_stages = NULL,
        ah_stages = NULL)
    }
  } else {
    # Stochastic elasticity analysis
    
    if(!any(is.na(time_weights))) {
      returned_cube <- .stoch_senselas(mats, times = steps, style = 2,
        tweights = time_weights)[[1]]
    } else {
      returned_cube <- .stoch_senselas(mats, times = steps, style = 2)[[1]]
    }
    
    returned_list <- lapply(as.list(c(1:dim(returned_cube)[3])), function(X) {
      return(returned_cube[,,X])
      }
    )
    
    if (historical) {
      output <- list(h_elasmats = returned_list, ah_elasmats = NULL,
        h_stages = NULL, ah_stages = NULL)
    } else {
      output <- list(h_elasmats = NULL, ah_elasmats = returned_list,
        h_stages = NULL, ah_stages = NULL)
    }
  }
  
  if (append_mats) {
    output$A <- mats
  }
  
  class(output) <- "lefkoElas"
  
  return(output)
}

#' Conduct a Life Table Response Experiment
#' 
#' \code{ltre3()} is a generic function that returns life table response
#' experiment (LTRE) or stochastic LTRE matrices for the input projection
#' matrices.
#' 
#' @param mats A lefkoMat object, population projection matrix, or list of
#' population projection matrices.
#' @param refmats A reference lefkoMat object, or matrix, for use as the
#' control. If missing, then is set to the same object as \code{mats}.
#' @param ... Other parameters.
#' 
#' @return The value returned depends on the class of the \code{mats} argument.
#' 
#' @section Notes:
#' Deterministic LTRE is one-way, fixed, and based on the sensitivities of the
#' matrix midway between each input matrix and the reference matrix, per Caswell
#' (2001, Matrix Population Models, Sinauer Associates, MA, USA). Stochastic
#' LTRE is per Davison et al. (2010, doi: 10.1111/j.1365-2745.2009.01611.x).
#' 
#' @seealso \code{\link{ltre3.lefkoMat}()}
#' @seealso \code{\link{summary.lefkoLTRE}()}
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
#' ltre3(ehrlen3)
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
#' ltre3(cypmatrix2r)
#' 
#' @export
ltre3 <- function(mats, refmats, ...) UseMethod("ltre3")

#' Conduct a Life Table Response Experiment of a lefkoMat Object
#' 
#' \code{ltre3.lefkoMat()} returns a set of matrices of one-way LTRE (life table
#' response experiment) or stochastic LTRE matrices contributions.
#' 
#' @param mats An object of class \code{lefkoMat}.
#' @param refmats A reference lefkoMat object, or matrix, for use as the
#' control. If missing, then is set to the same object as \code{mats}.
#' @param ref A numeric value indicating which matrix or matrices in
#' \code{refmats} to use as the control. The numbers used must correspond to the
#' number of the matrices in the \code{labels} element of the associated
#' \code{lefkoMat} object. The default setting, NA, uses all entries in
#' \code{refmats}.
#' @param stochastic A logical value determining whether to conduct a
#' deterministic (FALSE) or stochastic (TRUE) elasticity analysis. Defaults to
#' FALSE.
#' @param steps The number of occasions to project forward in stochastic
#' simulation. Defaults to 10,000.
#' @param burnin The number of initial steps to ignore in stochastic projection
#' when calculating stochastic elasticities. Must be smaller than \code{steps}.
#' Defaults to 3000.
#' @param time_weights Numeric vector denoting the probabilistic weightings of
#' all matrices. Defaults to equal weighting among matrices.
#' @param sparse A string indicating whether to use sparse matrix encoding
#' (\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
#' \code{"auto"}.
#' @param rseed Optional numeric value corresponding to the random seed for
#' stochastic simulation.
#' @param append_mats A logical value denoting whether to include the original
#' A, U, and F matrices in the returned \code{lefkoLTRE} object. Defaults to
#' FALSE.
#' @param ... Other parameters.
#' 
#' @return This function returns an object of class \code{lefkoLTRE}. This
#' includes a list of LTRE matrices as object \code{ltre_det} if a deterministic
#' LTRE is called for, or a list of mean-value LTRE matrices as object
#' \code{ltre_mean} and a list of SD-value LTRE matrices as object
#' \code{ltre_sd} if a stochastic LTRE is called for. This is followed by the
#' stageframe as object \code{ahstages}, the order of historical stages as
#' object \code{hstages}, the age-by-stage order as object \code{agestages}, the
#' order of matrices as object \code{labels}, and, if requested, the original A,
#' U, and F matrices.
#' 
#' @section Notes:
#' Deterministic LTRE is one-way, fixed, and based on the sensitivities of the
#' matrix midway between each input matrix and the reference matrix, per Caswell
#' (2001, Matrix Population Models, Sinauer Associates, MA, USA). Stochastic
#' LTRE is simulated per Davison et al. (2010) Journal of Ecology 98:255-267
#' (doi: 10.1111/j.1365-2745.2009.01611.x).
#'
#' Default behavior for stochastic LTRE uses the full population provided in
#' \code{mats} as the reference if no \code{refmats} and \code{ref} is provided.
#' If no \code{refmats} is provided but \code{ref} is, then the matrices noted
#' in \code{ref} are used as the reference matrix set. Year and patch order is
#' utilized from object \code{mats}, but not from object \code{refmats}, in
#' which each matrix is assumed to represent a different year from one
#' population. This function cannot currently handle multiple populations within
#' the same \code{mats} object (although such analysis is possible if these
#' populations are designated as patches instead).
#' 
#' @seealso \code{\link{ltre3}()}
#' @seealso \code{\link{summary.lefkoLTRE}()}
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
#' ltre3(ehrlen3, stochastic = TRUE)
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
#' ltre3(cypmatrix2r)
#' 
#' @export
ltre3.lefkoMat <- function(mats, refmats = NA, ref = NA, stochastic = FALSE,
  steps = 10000, burnin = 3000, time_weights = NA, sparse = "auto", rseed = NA,
  append_mats = FALSE, ...) {
  
  if (!all(is.na(rseed))) set.seed(rseed);
  
  if (is.element(tolower(sparse), c("y", "yes", "yea", "yeah", "t", "true"))) {
    sparsemethod <- 1
  } else if (is.element(tolower(sparse), c("n", "no", "non", "nah", "f", "false"))) {
    sparsemethod <- 0
  } else {
    elements_total <- length(mats$A[[1]])
    dense_elements <- length(which(mats$A[[1]] != 0))
    
    if ((dense_elements / elements_total) < 0.5) {
      sparsemethod <- 1
    } else sparsemethod <- 0
  }
  
  if (!is.element("agestages", names(mats))) {
    mats$agestages <- NULL
  }
  
  if (all(is.na(refmats))) {
    warning("Matrices input as mats will also be used as reference matrices.", call. = FALSE)
  } else if (is.element("lefkoMat", class(refmats))) {
    refmats <- refmats$A
  } else if (is.element("matrix", class(refmats))) {
    refmats <- list(refmats)
  } else {
    stop("Object refmats not recognized. Use only objects of class lefkoMat, list, or matrix.", call. = FALSE)
  }
  
  if (all(is.numeric(ref))) {
    if (length(ref) > 1) {
      if (!stochastic) {
        message("This function currently conducts deterministic LTREs against only single matrices, or mean matrices. Will use the mean.")
      }
    }
  } else if (all(is.na(ref))) {
    message("Using all refmats matrices as reference matrices.")
    
    if (!all(is.na(refmats))) {
      ref <- c(1:length(refmats))
    } else {
      ref <- c(1:length(mats$A))
    }
  }
  
  if (!stochastic & length(ref) > 1) {
    meanout <- 1
  } else if (stochastic & length(ref) == 1) {
    stop("Stochastic LTRE requires multiple matrices in object refmats.", call. = FALSE)
  } else if (stochastic & length(mats$A) == 1) {
    stop("Stochastic LTRE requires multiple matrices in object mats.", call. = FALSE)
  } else {
    meanout <- 0;
  }
  
  if (all(is.na(refmats))) {
    if (any(ref > length(mats$A)) | any(ref < 1)) {
      stop("Numbers of matrices provided in object ref must be integers between 1 and the total number of matrices in object mats (or refmats, if given).", call. = FALSE)
    }
  } else {
    
    if (dim(mats$A[[1]])[1] != dim(refmats[[1]])[1] | dim(mats$A[[1]])[2] != dim(refmats[[1]])[2]) {
      stop("Matrices used in objects mats and refmats must be of equal dimension.", call. = FALSE)
    }
  }
  
  ref <- ref - 1;
  
  if (!stochastic) {
    # Deterministic LTRE analysis
    
    if (all(is.na(refmats))) {
      baldrick <- .ltre3matrix(mats$A, refnum = ref, mean = meanout,
        sparse = sparsemethod)
    } else {
      baldrick <- .ltre3matrix(mats$A, refnum = ref, refmats_ = refmats,
        mean = meanout, sparse = sparsemethod)
    }
    
    returned_list <- lapply(as.list(c(1:dim(baldrick)[3])), function(X) {
      return(baldrick[,,X])
    })
    output <- list(ltre_det = returned_list, ahstages = mats$ahstages,
      agestages = mats$agestages, hstages = mats$hstages, labels = mats$labels)
    
  } else {
    # Stochastic LTRE analysis
    
    if (burnin >= steps) {
      stop("Option burnin must be smaller than option steps.", call. = FALSE)
    }
    
    #LOY table development
    listofyears <- mats$labels
    
    if (all(is.na(listofyears$pop))) {
      listofyears$pop <- "1"
    }
    
    if (all(is.na(listofyears$patch))) {
      listofyears$patch <- "1"
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
    
    numofpops <- length(unique(listofyears$popc))
    numofpatches <- length(unique(listofyears$poppatchc))
    numofyears <- length(unique(listofyears$year2c))
    
    listofyears$poppatchc <- as.numeric(listofyears$poppatchc)
    
    if (all(is.na(refmats))) {
      if (all(is.na(time_weights))) {
        baldrick <- .sltre3matrix(mats$A, loy = listofyears, refnum = ref,
          steps = steps, burnin = burnin, sparse = sparsemethod)
      } else {
        baldrick <- .sltre3matrix(mats$A, loy = listofyears, refnum = ref,
          tweights_ = time_weights, steps = steps, burnin = burnin,
          sparse = sparsemethod)
      }
    } else {
      if (all(is.na(time_weights))) {
        baldrick <- .sltre3matrix(mats$A, loy = listofyears, refnum = ref,
          refmats_ = refmats, steps = steps, burnin = burnin,
          sparse = sparsemethod)
      } else {
        baldrick <- .sltre3matrix(mats$A, loy = listofyears, refnum = ref,
          refmats_ = refmats, tweights_ = time_weights, steps = steps,
          burnin = burnin, sparse = sparsemethod)
      }
    }
    
    alabels <- mats$labels[,c("pop", "patch")]
    alabels <- unique(alabels)
    
    output <- list(ltre_mean = baldrick$cont_mean, ltre_sd = baldrick$cont_sd,
      ahstages = mats$ahstages, agestages = mats$agestages,
      hstages = mats$hstages, labels = alabels)
  }
  
  if (append_mats) {
    output$A <- mats$A
    output$U <- mats$U
    output$F <- mats$F
  }
  
  rownames(output$labels) <- c(1:length(output$labels$patch))
  
  class(output) <- "lefkoLTRE"
  
  return(output)
}

#' Summarize lefkoElas Objects
#' 
#' Function \code{summary.lefkoElas()} summarizes \code{lefkoElas} objects.
#' Particularly, it breaks down elasticity values by the kind of ahistorical
#' and, if applicable, historical transition.
#'
#' @param object A \code{lefkoElas} object.
#' @param ... Other parameters currently not utilized.
#' 
#' @return A list composed of 2 data frames. The first, \code{hist}, is a data
#' frame showing the summed elasticities for all 16 kinds of historical
#' transition per matrix, with each column corresponding to each elasticity
#' matrix in order. The second, \code{ahist}, is a data frame showing the
#' summed elasticities for all 4 kinds of ahistorical transition per matrix,
#' with each column corresponding to each elasticity matrix in order.
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
#' lathsupp2 <- supplemental(stage3 = c("Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "rep", "rep"),
#'   givenrate = c(0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 3, 3), stageframe = lathframe, historical = FALSE)
#'   
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = "all", 
#'   stages = c("stage3", "stage2", "stage1"), supplement = lathsupp3,
#'   yearcol = "year2", indivcol = "individ")
#' 
#' ehrlen2 <- rlefko2(data = lathvert, stageframe = lathframe, year = "all",
#'   stages = c("stage3", "stage2"), supplement = lathsupp2,
#'   yearcol = "year2", indivcol = "individ")
#' 
#' ehrlen3elas <- elasticity3(ehrlen3)
#' ehrlen2elas <- elasticity3(ehrlen2)
#' 
#' summary(ehrlen3elas)
#' summary(ehrlen2elas)
#' 
#' @export
summary.lefkoElas <- function(object, ...) {
  elasmats <- object
  
  num_h_mats <- length(elasmats$h_elasmats)
  num_ah_mats <- length(elasmats$ah_elasmats)
  
  used_emats <- if(num_h_mats == 0) elasmats$ah_elasmats else elasmats$h_elasmats
  
  used_iterations <- if(num_h_mats > 0) num_h_mats else num_ah_mats
  
  if (num_h_mats == 0) {
    if (!all(is.null(elasmats$agestages))) {
      if (!all(is.na(elasmats$agestages))) {
        if (is.element("stage_id", names(elasmats$agestages))) {
          new_ahstages_list <- apply(as.matrix(c(1:length(elasmats$agestages$stage_id))), 
            1, function(X) {
              return(elasmats$ah_stages[which(elasmats$ah_stages$stage_id == elasmats$agestages$stage_id[X]),])
            })
          new_ahstages <- do.call("rbind.data.frame", new_ahstages_list)
          indices <- .bambi2(new_ahstages)
        } else {
          indices <- .bambi2(elasmats$ah_stages)
        }
      } else {
        indices <- .bambi2(elasmats$ah_stages)
      }
    } else {
      indices <- .bambi2(elasmats$ah_stages)
    }
    
  } else {
    indices <- .bambi3(elasmats$ah_stages, elasmats$h_stages)
  }
  
  for (i in c(1:used_iterations)) {
    trialguy <- .demolition3(used_emats[[i]], indices)
    
    if (i == 1) {
      if (num_h_mats > 0) hist <- trialguy$hist[,c(1,2)]
      ahist <- trialguy$ahist[,c(1,2)]
      
      if (num_h_mats > 0) names(hist)[2] <- "matrix1"
      names(ahist)[2] <- "matrix1"
    } else {
      if (num_h_mats > 0) hist <- cbind.data.frame(hist, trialguy$hist[,2])
      ahist <- cbind.data.frame(ahist, trialguy$ahist[,2])
      if (num_h_mats > 0) names(hist)[(i+1)] <- paste0("matrix", i)
      names(ahist)[(i+1)] <- paste0("matrix", i)
    }
  }
  
  if (num_h_mats > 0) {
    output <- list(hist = hist, ahist = ahist)
  } else {
    output <- list(hist = NULL, ahist = ahist)
  }
  
  return (output)
}

#' Summarize lefkoLTRE Objects
#' 
#' Function \code{summary.lefkoLTRE()} summarizes \code{lefkoLTRE} objects.
#' Particularly, it breaks down LTRE contributions by the kind of ahistorical
#' and, if applicable, historical transition.
#'
#' @param object A \code{lefkoLTRE} object.
#' @param ... Other parameters currently not utilized.
#' 
#' @return A list composed of 2 (if deterministic) or 4 (if stochastic) data
#' frames. If deterministic, then \code{hist_det} is a data
#' frame showing the summed LTRE contributions for all 16 kinds of historical
#' transition per matrix, with each column corresponding to each A matrix in
#' order, followed by all summed positive and all summed negative contributions.
#' Object \code{ahist_det} is a data frame showing the summed LTRE
#' contributions for all 4 kinds of ahistorical transition per matrix, with
#' order as before, followed by summed positive and summed negative
#' contributions. If stochastic, then \code{hist_mean} and \code{hist_sd} are
#' the summed LTRE contributions for the mean vital rates and variability in
#' vital rates, respectively, according to all 16 historical transition types,
#' followed by summed positive and negative contributions, and \code{ahist_mean}
#' and \code{ahist_sd} are the equivalent ahistorical versions.
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
#' lathsupp2 <- supplemental(stage3 = c("Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "rep", "rep"),
#'   givenrate = c(0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 3, 3), stageframe = lathframe, historical = FALSE)
#'   
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = "all", 
#'   stages = c("stage3", "stage2", "stage1"), supplement = lathsupp3,
#'   yearcol = "year2", indivcol = "individ")
#' 
#' ehrlen2 <- rlefko2(data = lathvert, stageframe = lathframe, year = "all",
#'   stages = c("stage3", "stage2"), supplement = lathsupp2,
#'   yearcol = "year2", indivcol = "individ")
#' 
#' ehrlen3ltre <- ltre3(ehrlen3)
#' summary(ehrlen3ltre)
#' 
#' @export
summary.lefkoLTRE <- function(object, ...) {
  
  ltremats <- object
  
  if (is.element("ltre_det", names(object))) {
    ltretype <- 1 # Deterministic LTRE
  } else if (is.element("ltre_sd", names(object))) {
    ltretype <- 2
  } else {
    stop("Input object is of unrecognized type. Use only lefkoLTRE obects with this function.", call. = FALSE)
  }
  
  numstages <- if (ltretype == 1) {
    dim(object$ltre_det[[1]])[1]
  } else {
    dim(object$ltre_sd[[1]])[1]
  }
  
  if(!all(is.na(object$hstages))) {
    if(numstages == dim(object$hstages)[1]) {
      historical <- TRUE
    } else {
      historical <- FALSE
    }
  } else {
    historical <- FALSE
  }
  
  if (!historical) {
    if (!all(is.null(object$agestages))) {
      if (!all(is.na(object$agestages))) {
        if (is.element("stage_id", names(object$agestages))) {
          new_ahstages_list <- apply(as.matrix(c(1:length(object$agestages$stage_id))), 
            1, function(X) {
              return(object$ahstages[which(object$ahstages$stage_id == object$agestages$stage_id[X]),])
            })
          new_ahstages <- do.call("rbind.data.frame", new_ahstages_list)
          indices <- .bambi2(new_ahstages)
        } else {
          indices <- .bambi2(object$ahstages)
        }
      } else {
        indices <- .bambi2(object$ahstages)
      }
    } else {
      indices <- .bambi2(object$ahstages)
    }
  } else {
    indices <- .bambi3(object$ahstages, object$hstages)
  }
  
  used_iterations <- if(ltretype == 1) {
    length(object$ltre_det)
  } else {
    length(object$ltre_sd)
  }
  
  for (i in c(1:used_iterations)) {
    trialguy1 <- if (ltretype == 1) {
      .demolition3(object$ltre_det[[i]], indices)
    } else {
      .demolition3(object$ltre_mean[[i]], indices)
    }
    
    if (ltretype == 2) {
      trialguy2 <- .demolition3(object$ltre_sd[[i]], indices)
    }
    
    if (i == 1) {
      hist1 <- trialguy1$hist
      ahist1 <- trialguy1$ahist
      if (historical) names(hist1)[2] <- "matrix1"
      if (historical) names(hist1)[3] <- "matrix1_pos"
      if (historical) names(hist1)[4] <- "matrix1_neg"
      names(ahist1)[2] <- "matrix1"
      names(ahist1)[3] <- "matrix1_pos"
      names(ahist1)[4] <- "matrix1_neg"
      
      if (ltretype == 2) {
        hist2 <- trialguy2$hist
        ahist2 <- trialguy2$ahist
        if (historical) names(hist2)[2] <- "matrix1"
        if (historical) names(hist2)[3] <- "matrix1_pos"
        if (historical) names(hist2)[4] <- "matrix1_neg"
        names(ahist2)[2] <- "matrix1"
        names(ahist2)[3] <- "matrix1_pos"
        names(ahist2)[4] <- "matrix1_neg"
      }
    } else {
      if (historical) hist1 <- cbind.data.frame(hist1, trialguy1$hist[,c(2:4)])
      ahist1 <- cbind.data.frame(ahist1, trialguy1$ahist[,c(2:4)])
      if (historical) names(hist1)[(((i-1)*3)+2)] <- paste0("matrix", i)
      if (historical) names(hist1)[(((i-1)*3)+3)] <- paste0("matrix", i, "_pos")
      if (historical) names(hist1)[(((i-1)*3)+4)] <- paste0("matrix", i, "_neg")
      names(ahist1)[(((i-1)*3)+2)] <- paste0("matrix", i)
      names(ahist1)[(((i-1)*3)+3)] <- paste0("matrix", i, "_pos")
      names(ahist1)[(((i-1)*3)+4)] <- paste0("matrix", i, "_neg")
      
      if (ltretype == 2) {
        if (historical) hist2 <- cbind.data.frame(hist2, trialguy2$hist[,c(2:4)])
        ahist2 <- cbind.data.frame(ahist2, trialguy2$ahist[,c(2:4)])
        if (historical) names(hist2)[(((i-1)*3)+2)] <- paste0("matrix", i)
        if (historical) names(hist2)[(((i-1)*3)+3)] <- paste0("matrix", i, "_pos")
        if (historical) names(hist2)[(((i-1)*3)+4)] <- paste0("matrix", i, "_neg")
        names(ahist2)[(((i-1)*3)+2)] <- paste0("matrix", i)
        names(ahist2)[(((i-1)*3)+3)] <- paste0("matrix", i, "_pos")
        names(ahist2)[(((i-1)*3)+4)] <- paste0("matrix", i, "_neg")
      }
    }
  }
  
  output <- if (ltretype == 1) {
    list(hist_det = hist1, ahist_det = ahist1)
  } else {
    list(hist_mean = hist1, hist_sd = hist2, ahist_mean = ahist1, ahist_sd = ahist2)
  }
  
  return (output)
}

#' Summarize lefkoProj Objects
#' 
#' Function \code{summary.lefkoProj()} summarizes \code{lefkoProj} objects.
#' Particularly, it breaks down the data frames provided in the 
#' \code{projection} element in ways meaningful for those running simulations.
#'
#' @param object A \code{lefkoProj} object.
#' @param threshold A threshold population size to be searched for in
#' projections. Defaults to 1.
#' @param milepost A numeric vector indicating at which points in the projection
#' to assess detailed results. Can be input as integer values, in which case
#' each number must be between 1 and the total number of occasions projected in
#' each projection, or decimals between 0 and 1, which would then be translated
#' into the corresponding projection steps of the total. Defaults to
#' \code{c(0, 0.25, 0.50, 0.75, 1.00)}.
#' @param ... Other parameters currently not utilized.
#' 
#' @return Apart from a statement of the results, we have the following item in
#' the output:
#' \item{milepost_sums}{A data frame showing the number of replicates at each
#' of the milepost times that is above the threshold population/patch size.}
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
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "Sd", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 1, 3, 3), type_t12 = c(1, 2, 1, 2, 1, 1),
#'   stageframe = lathframe, historical = TRUE)
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe,
#'   year = c(1989, 1990), stages = c("stage3", "stage2", "stage1"),
#'   repmatrix = lathrepm, supplement = lathsupp3, yearcol = "year2",
#'   indivcol = "individ")
#' 
#' lathproj <- projection3(ehrlen3, nreps = 5, stochastic = TRUE)
#' summary(lathproj)
#' 
#' # Cypripedium example
#' rm(list = ls(all=TRUE))
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
#'   supplement = cypsupp3r, yearcol = "year2", 
#'   patchcol = "patchid", indivcol = "individ")
#' 
#' cypstoch <- projection3(cypmatrix3r, nreps = 5, stochastic = TRUE)
#' summary(cypstoch)
#' 
#' @export
summary.lefkoProj <- function(object, threshold = 1,
  milepost = c(0, 0.25, 0.50, 0.75, 1.00), ...) {
  
  poppatches <- length(object$pop_size)
  
  nreps <- object$control[1]
  times <- object$control[2]
  
  if (any(milepost < 0)) {
    stop("Option milepost may not take negative values.", call. = FALSE)
  }
  if (any(milepost > times)) {
    stop("Option milepost may not take values higher than the number of actual number of projected occasions.", call. = FALSE)
  }
  
  if (all(milepost >= 0) & all(milepost <= 1)) {
    milepost <- floor(milepost * times) + 1
  } else if (any(milepost == 0)) {
    milepost <- milepost + 1
  }
  
  if (nreps > 1) {
    milepost_sums <- apply(as.matrix(c(1:poppatches)), 1, function (X) {
      phew <- apply(as.matrix(object$pop_size[[X]][,milepost]), 2, function(Y) {
        return(length(which(as.vector(Y) >= threshold)))
      })
      return(phew)
    })
  } else {
     milepost_sums <- apply(as.matrix(c(1:poppatches)), 1, function (X) {
      phew <- apply(as.matrix(object$pop_size[[X]][,milepost]), 1, function(Y) {
        return(length(which(as.vector(Y) >= threshold)))
      })
      return(phew)
    })
  }
  
  if (is.element("matrix", class(milepost_sums))) {
    rownames(milepost_sums) <- milepost
    
    col_labels <- apply(object$labels, 1, function(X) {
      paste(X[1], X[2])
    })
    colnames(milepost_sums) <- col_labels
  } else {
    rownames(milepost_sums) <- milepost
  }
  
  writeLines(paste0("\nThe input lefkoProj object covers ", poppatches,
    " population-patches."), con = stdout())
  writeLines(paste0("It includes ", times, " projected steps per replicate and ",
    nreps, " replicates."), con = stdout())
  writeLines(paste0("The number of replicates with population size above the threshold size of ", threshold,
    " is as in"), con = stdout())
  writeLines(paste0("the following matrix, with pop-patches given by column and milepost times given by row: \n"),
    con = stdout())
  
  return (milepost_sums)
}

#' Plot Projection Simulations
#' 
#' Function \code{plot.lefkoProj()} produces plots of \code{lefkoProj} objects.
#' Acts as a convenient wrapper for the \code{plot.default()} function.
#' 
#' @param x A \code{lefkoProj} object.
#' @param variable The focus variable of the plot to produce. Defaults to
#' \code{"popsize"}, which produces line plots of the \code{popsize} element in
#' object \code{x}.
#' @param style A string denoting ther kind of plot to produce. Currently
#' limited to \code{"timeseries"}, which shows \code{variable} against time on
#' the x axis. Other choices include \code{"statespace"}, which plots
#' \code{variable} at one time on the x axis against the same variable in the
#' next time on the y axis.
#' @param repl The replicate to plot. Defaults to \code{"all"}, in which case
#' all replicates are plotted.
#' @param patch The patch to plot, as labeled in the \code{labels} element in
#' object \code{x}. Defaults to \code{"pop"}, in which case only the final
#' population-level projection is plotted. Can also be set to \code{"all"}, in
#' which case projections for all patches and population in the \code{labels}
#' element are plotted.
#' @param auto_ylim A logical value indicating whether the maximum of the y axis
#' should be determined automatically. Defaults to \code{TRUE}, but reverts to
#' \code{FALSE} if any setting for \code{ylim} is given.
#' @param auto_col A logical value indicating whether to shift the color of
#' lines associated with each patch automatically. Defaults to \code{TRUE}, but
#' reverts to \code{FALSE} if any setting for \code{col} is given.
#' @param auto_lty A logical value indicating whether to shift the line type
#' associated with each replicate automatically. Defaults to \code{TRUE}, but
#' reverts to \code{FALSE} if any setting for \code{lty} is given.
#' @param ... Other parameters used by functions \code{plot.default()} and
#' \code{lines()}.
#' 
#' @return A plot of the results of a \code{link{projection3}()} run.
#' 
#' @section Notes:
#' Output plots are currently limited to time series and state space plots of
#' population size.
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
#' lathrepm <- matrix(0, 7, 7)
#' lathrepm[1, 6] <- 0.345
#' lathrepm[2, 6] <- 0.054
#' 
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "Sd", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 1, 3, 3), type_t12 = c(1, 2, 1, 2, 1, 1),
#'   stageframe = lathframe, historical = TRUE)
#' 
#' ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe,
#'   year = c(1989, 1990), stages = c("stage3", "stage2", "stage1"),
#'   repmatrix = lathrepm, supplement = lathsupp3, yearcol = "year2",
#'   indivcol = "individ")
#' 
#' lathproj <- projection3(ehrlen3, nreps = 5, stochastic = TRUE)
#' plot(lathproj)
#' 
#' @export
plot.lefkoProj <- function(x, variable = "popsize", style = "time",
  repl = "all", patch = "pop", auto_ylim = TRUE, auto_col = TRUE,
  auto_lty = TRUE, ...) {
  
  further_args <- list(...)
  
  if (length(further_args) == 0) further_args <- list()
  
  if (!is.element("type", names(further_args))) {
    further_args$type <- "l"
  }
  if (!is.element("ylab", names(further_args))) {
    further_args$ylab <- "Population size"
  }
  if (!is.element("xlab", names(further_args))) {
    further_args$xlab <- "Time"
  }
  if (is.element("col", names(further_args))) {
    auto_col <- FALSE
  }
  if (is.element("lty", names(further_args))) {
    auto_lty <- FALSE
  }
  if (is.element("ylim", names(further_args))) {
    auto_ylim <- FALSE
  }
  basal_args <- further_args
  
  actual_patches <- c(1:length(x$labels$patch))
  actual_replicates <- c(1:x$control[1])
  
  if (length(grep("pop", variable)) > 0 | length(grep("size", variable)) > 0) {
    variable <- "popsize"
  }
  
  if (length(grep("stat", style)) > 0) {
    style <- "statespace"
  } else if (length(grep("tim", style)) > 0) {
    style <- "timeseries"
  }
  
  if (all(is.na(patch))) {
    patch <- actual_patches
  } else if (is.character(patch)) {
    if (any(grep("al", patch))) {
      patch <- actual_patches
    } else if (length(patch) == 1) {
      if (grep("po", patch)) {
        if (length(actual_patches) == 1) {
          patch <- actual_patches
        } else {
          patch <- which(x$labels$patch == 0)
        }
      }
    } else {
      patch <- as.numeric(patch)
      
      if (any(is.na(patch))) {
        stop("Setting patch not understood. Please do not combine text with numbers.",
          call. = FALSE)
      }
      if (!all(is.element(patch, actual_patches))) {
        stop("Setting patch not understood. Please check the labels element of
          your input lefkoProj object.", call. = FALSE)
      }
    }
  } else if (is.numeric(patch)) {
    if (any(!is.element(patch, actual_patches))) {
      stop("Setting patch not understood. Please check the labels element of
          your input lefkoProj object.", call. = FALSE)
    }
  }
  
  if (is.character(repl)) {
    if (any(grep("al", repl))) {
      repl <- actual_replicates
    } else {
      stop("Setting repl not understood.", call. = FALSE)
    }
  } else {
    if (any(!is.element(repl, actual_replicates))) {
      stop("Setting repl not understood.", call. = FALSE)
    }
  }
  
  if (variable != "popsize") {
    stop("Function plot.lefkoProj() does not current handle plots of elements
      other than `popsize`.", call. = FALSE)
  }
  
  used_col <- 1
  
  for (i in patch) {
    used_lty <- 1
    
    if (auto_ylim) {
      further_args$ylim <- c(0, max(x$pop_size[[i]], na.rm = TRUE))
    }
    if (auto_col) {
      further_args$col <- used_col
      basal_args$col <- used_col
    }
    if (auto_lty) {
      further_args$lty <- used_lty
      basal_args$lty <- used_lty
    }
      
    if (style == "timeseries") {
      c_xy <- xy.coords(x = c(1:length(x$pop_size[[i]][1,])), y = x$pop_size[[i]][1,])
      further_args$x <- c_xy
      
    } else if (style == "statespace") {
      c_xy <- xy.coords(x = x$pop_size[[i]][1,1:(dim(x$pop_size[[i]])[2] - 1)],
        y = x$pop_size[[i]][1,2:dim(x$pop_size[[i]])[2]])
      further_args$x <- c_xy
    }
    
    do.call("plot.default", further_args)
    
    if (length(repl) > 1) {
      for (j in c(2:x$control[1])) {
        if (style == "timeseries") {
          basal_args$x <- c(1:length(x$pop_size[[i]][1,]))
          basal_args$y <- x$pop_size[[i]][1,]
          
        } else if (style == "statespace") {
          basal_args$y <- x$pop_size[[i]][j,2:dim(x$pop_size[[i]])[2]]
          basal_args$x <- x$pop_size[[i]][j,1:(dim(x$pop_size[[i]])[2] - 1)]
        
        }
        do.call("lines", basal_args)
      }
      used_lty <- used_lty + 1;
    }
    used_col <- used_col + 1;
    if (used_col > length(palette())) used_col <- 1;
  }
}

