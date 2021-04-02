#' Estimate Mean Projection Matrices
#' 
#' \code{lmean()} estimates mean projection matrices as element-wise arithmetic
#' means.
#' 
#' @param mats A \code{lefkoMat} object.
#' @param matsout A string identifying which means to estimate. Option "pop"
#' indicates population-level only, "patch" indicates patch-level only, and
#' "all" indicates that both patch- and population-level means should be
#' estimated.
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
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 3, 3), stageframe = lathframe, historical = TRUE)
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, 0.40, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
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
  
  if (class(mats) != "lefkoMat") {
    stop("An object of class lefkoMat is required as input.")
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
    output <- turbogeodiesel(listofyears, mats$U, mats$F, mats$ahstages, mats$hstages, patchonly, poponly)
  } else {
    output <- geodiesel(listofyears, mats$U, mats$F, mats$ahstages, patchonly, poponly)
    output$hstages <- NA
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
#' 
#' @return The value returned depends on the class of the \code{mats} argument.
#' See related functions for more details.
#' 
#' @seealso \code{\link{lambda3.lefkoMat}()}
#' @seealso \code{\link{lambda3.matrix}()}
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
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 3, 3), stageframe = lathframe, historical = TRUE)
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, 0.40, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
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
lambda3 <- function(mats) UseMethod("lambda3")

#' Estimate Deterministic Population Growth Rates of lefkoMat Matrices
#' 
#' \code{lambda3.lefkoMat()} returns the dominant eigenvalues of all projection
#' matrices supplied within \code{lefkoMat} objects. This function can handle
#' large and sparse matrices, and so can be used with large historical matrices,
#' IPMs, age x stage matrices, as well as smaller ahistorical matrices.
#' 
#' @param mats An object of class \code{lefkoMat}.
#' 
#' @return This function returns the dominant eigenvalue of each \code{$A}
#' matrix in \code{mats}. For square matrices with fewer than 400 rows, this is
#' given as the largest real part of all eigenvalues estimated via the
#' \code{eig_gen}() function in the C++ Armadillo library. For larger matrices,
#' the function assumes that matrices are sparse and uses \code{eigs_gen}()
#' instead (the function handles all conversions automatically). The output
#' includes a data frame showing the population, patch, and lambda estimate for
#' each \code{$A} matrix. Row names correspond to the order of the matrix within
#' the \code{$A} element of \code{mats}.
#' 
#' @seealso \code{\link{lambda3}()}
#' @seealso \code{\link{lambda3.matrix}()}
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
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 3, 3), stageframe = lathframe, historical = TRUE)
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, 0.40, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
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
lambda3.lefkoMat <- function(mats) {
  baldrick <- if (any(class(mats$A) == "matrix")) {
    
    if (dim(mats$A)[1] > 400) {
      lambda3matrixsp(mats$A)
    } else {
      lambda3matrix(mats$A)
    }
    
  } else if (class(mats$A) == "list") {
    
    if (dim(mats$A[[1]])[1] > 400) {
      unlist(lapply(mats$A, lambda3matrixsp))
    } else {
      unlist(lapply(mats$A, lambda3matrix))
    }
    
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
#' as smaller ahistorical matrices.
#' 
#' @param mats A population projection matrix of class \code{matrix}.
#'
#' @return This function returns the dominant eigenvalue of the matrix. For
#' square matrices with fewer than 400 rows, this is given as the largest real
#' part of all eigenvalues estimated via the \code{eig_gen}() function in
#' package the C++ Armadillo library. For larger matrices, the matrix is assumed
#' to be sparse and \code{eigs_gen}() is used instead.
#' 
#' @seealso \code{\link{lambda3}()}
#' @seealso \code{\link{lambda3.lefkoMat}()}
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
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 3, 3), stageframe = lathframe, historical = TRUE)
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, 0.40, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
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
lambda3.matrix <- function(mats)
{
  if (dim(mats)[1] > 400) {
    lambda <- lambda3matrixsp(mats);
  } else {
    lambda <- lambda3matrix(mats);
  }
  
  return(lambda)
}

#' Estimate Stable Stage Distribution
#' 
#' \code{stablestage3()} is a generic function that returns the stable stage 
#' distribution for a population projection matrix or set of matrices. This
#' function is made to handle very large and sparse matrices supplied as 
#' \code{lefkoMat} objects or as individual matrices, and can be used with large
#' historical matrices, IPMs, age x stage matrices, as well as smaller
#' ahistorical matrices.
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
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 3, 3), stageframe = lathframe, historical = TRUE)
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, 0.40, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
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
#' distributions of all \code{$A} matrices in an object of class
#' \code{lefkoMat}, as well as the long-run projected mean stage distribution in
#' stochastic analysis. This function can handle large and sparse matrices, and
#' so can be used with large historical matrices, IPMs, age x stage matrices, as
#' well as smaller ahistorical matrices.
#' 
#' @param mats An object of class \code{lefkoMat}.
#' @param stochastic A logical value indicating whether to use deterministic
#' (\code{FALSE}) or stochastic (\code{TRUE}) analysis. Defaults to
#' \code{FALSE}.
#' @param times An integer variable indicating number of times to project if
#' using stochastic analysis. Defaults to 10000.
#' @param tweights An optional vector indicating the probability weighting to
#' use for each matrix in stochastic simulations. If not given, then defaults to
#' equal weighting.
#' @param seed A number to use as a random number seed.
#' @param ... Other parameters.
#' 
#' @return This function returns the stable stage distributions (and long-run
#' mean stage distributions in stochastic analysis) corresponding to the
#' matrices in a \code{lefkoMat} object.
#' 
#' The output depends on whether the \code{lefkoMat} object used as input is
#' ahistorical or historical, and whether the analysis is deterministic or
#' stochastic. If ahistorical, then a single data frame is output, which
#' includes the number of the matrix within the \code{$A} element of the input
#' \code{lefkoMat} object, followed by the stage id (numeric and assigned
#' through \code{\link{sf_create}()}), the stage name, and the estimated
#' proportion of the stable stage distribution (\code{ss_prop}).
#' 
#' If a historical matrix is used as input, then two data frames are output
#' into a list object. The \code{$hist} element contains a data frame in which 
#' the stable stage distribution is given in terms of across-year stage pairs.
#' The structure includes the matrix number, the numeric stage designations for
#' stages in times \emph{t} and \emph{t}-1, respectively, followed by the
#' respective stage names, and ending with the estimated proportion of the
#' stable stage distribution for that stage within its matrix (\code{ss_prop}).
#' The \code{$ahist} element contains the stable stage distribution in stages
#' as given in the original stageframe. It includes a data frame with the matrix 
#' of origin, the numeric stage designation, stage name, and the stable stage
#' distribution estimated as the sum of distribution elements from \code{$hist}
#' corresponding to the equivalent stage in time \emph{t}, irrespective of stage
#' in time \emph{t}-1.
#'
#' In addition to the data frames noted above, stochastic analysis will result
#' in the additional output of a list of matrices containing the actual
#' projected stage distributions across all projected times, in the order of
#' population-patch combinations in the \code{lefkoMat} input.
#'
#' @section Notes:
#' For square matrices with fewer than 400 rows, the stable stage distribution
#' is given as the right eigenvector associated with largest real part of all
#' eigenvalues estimated via the \code{eig_gen}() function in the C++ Armadillo
#' library divided by the sum of the associated right eigenvector. For larger
#' matrices, the function assumes that the matrix is sparse and conducts a
#' similar calculation but using the \code{eigs_gen}() for sparse matrix eigen
#' analysis.
#' 
#' In stochastic analysis, the projected mean distribution is the arithmetic
#' mean across the final 1000 projected times if the simulation is at least 2000
#' projected times long. If between 500 and 2000 projected times long, then only
#' the final 200 are used, and if fewer than 500 times are used, then all are
#' used. Note that because stage distributions in stochastic simulations can
#' change greatly in the initial portion of the run, we encourage a minimum of
#' 2000 projected times per simulation, with 10,000 preferred.
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
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 3, 3), stageframe = lathframe, historical = TRUE)
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, 0.40, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
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
  tweights = NA, seed = NA, ...) {
  
  if (!stochastic) {
    baldrick <- if (any(class(mats$A) == "matrix")) {
      
      if (dim(mats$A)[1] > 400) {
        ss3matrixsp(mats$A)
      } else {
        ss3matrix(mats$A)
      }
      
    } else if (class(mats$A) == "list") {
      
      if (dim(mats$A[[1]])[1] > 400) {
        unlist(lapply(mats$A, ss3matrixsp))
      } else {
        unlist(lapply(mats$A, ss3matrix))
      }
      
    } else {
      
      stop("Input not recognized.")
    }
  } else {
    
    if (!is.na(seed)) {
      set.seed(seed)
    }
    
    mats$labels$poppatch <- paste(mats$labels$pop, mats$labels$patch)
    used_poppatches <- as.list(unique(mats$labels$poppatch))
    
    #Here we get the full stage distribution series for all times, as a list
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
      starter <- if (dim(mats$A[[1]])[1] > 400) {
        ss3matrixsp(mats$A[[used_slots[1]]])
      } else {
        ss3matrix(mats$A[[used_slots[1]]])
      }
      
      theseventhmatrix <- proj3(starter, mats$A, theprophecy, 1, 0, 0)
      
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
      output <- cbind.data.frame(modlabels, baldrick)
      if (is.element("age", names(labels))) {
        names(output) <- c("poppatch", "age", "stage_id", "stage", "agestage_id", "agestage", "ss_prop")
      } else {
        names(output) <- c("poppatch", "stage_id", "stage", "ss_prop")
      }
    }
    rownames(output) <- c(1:dim(output)[1])
  } else {
    labels <- mats$hstages
    
    modlabels <- cbind.data.frame(as.matrix(rep(c(1:multiplier), each = dim(labels)[1])), 
      do.call("rbind.data.frame", apply(as.matrix(c(1:multiplier)), 1, function(X) return(labels))))
    
    outputh <- cbind.data.frame(modlabels, baldrick)
    names(outputh) <- c("matrix", "stage_id_2", "stage_id_1", "stage_2", "stage_1", "ss_prop")
    rownames(outputh) <- c(1:dim(outputh)[1])
    
    ahlabels <- mats$ahstages[,c("stage_id", "stage")]
    ss2 <- c(apply(as.matrix(c(1:multiplier)), 1, function(X) {
      rightset <- subset(outputh, matrix == X)
      apply(as.matrix(ahlabels[,1]), 1, function(Y) {
        sum(rightset$ss_prop[which(rightset$stage_id_2 == Y)])
      })
    }))
    outputah <- cbind.data.frame(rep(c(1:multiplier), each = length(ahlabels[,1])), 
      rep(ahlabels[,1], multiplier), rep(ahlabels[,2], multiplier), ss2)
    names(outputah) <- c("matrix", "stage_id", "stage", "ss_prop")
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
#' @param ... Other parameters.
#' 
#' @return This function returns the stable stage distribution corresponding to
#' the input matrix. For square matrices with fewer than 400 rows, the stable
#' stage distribution is given as the right eigenvector associated with largest
#' real part of the eigenvalues estimated for the matrix via the
#' \code{eig_gen}() function in the C++ Armadillo library, divided by the sum of
#' the associated right eigenvector. For larger matrices, the matrix is assumed
#' to be sparse and the calculation is conducted similarly but using
#' \code{eig_gens}() instead.
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
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 3, 3), stageframe = lathframe, historical = TRUE)
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, 0.40, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
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
stablestage3.matrix <- function(mats, ...)
{
  if (dim(mats)[1] > 400) {
    wcorr <- ss3matrixsp(mats)
  } else {
    wcorr <- ss3matrix(mats)
  }
  
  return(wcorr)
}

#' Estimate Reproductive Value
#' 
#' \code{repvalue3()} is a generic function that estimates returns the
#' reproductive values of stages in a population projection matrix or a set of
#' matrices. The specifics of estimation vary with the class of input object.
#' This function is made to handle very large and sparse matrices supplied as
#' \code{lefkoMat} objects or as individual matrices, and can be used with large
#' historical matrices, IPMs, age x stage matrices, as well as smaller
#' ahistorical matrices.
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
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 3, 3), stageframe = lathframe, historical = TRUE)
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, 0.40, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
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
#' large historical matrices, IPMs, age x stage matrices, as well as smaller
#' ahistorical matrices.
#' 
#' @param mats An object of class \code{lefkoMat} object.
#' @param stochastic A logical value indicating whether to use deterministic
#' (\code{FALSE}) or stochastic (\code{TRUE}) analysis. Defaults to
#' \code{FALSE}.
#' @param times An integer variable indicating number of times to project if
#' using stochastic analysis. Defaults to 10000.
#' @param tweights An optional vector indicating the probability weighting to
#' use for each matrix in stochastic simulations. If not given, then defaults to
#' equal weighting.
#' @param seed A number to use as a random number seed.
#' @param ... Other parameters.
#' 
#' @return This function returns the asymptotic reproductive value vectors if
#' deterministic analysis chosen, and long-run mean reproductive value vectors
#' if stochastic analysis was chosen.
#' 
#' The output depends on whether the \code{lefkoMat} object used as input is
#' ahistorical or historical, and whether the analysis is deterministic or
#' stochastic. If ahistorical, then a single data frame is output, which
#' includes the number of the matrix within the \code{$A} element of the input
#' \code{lefkoMat} object, followed by the stage id (numeric and assigned
#' through \code{\link{sf_create}()}), the stage name, and the estimated
#' reproductive value (\code{rep_value}). Reproductive values are scaled by the
#' first non-zero value.
#' 
#' If a historical matrix is used as input, then two data frames are output
#' into a list object. The \code{$hist} element contains a data frame in which 
#' the stable stage distribution is given in terms of across-year stage pairs.
#' The structure includes the matrix number, the numeric stage designations for
#' stages in times \emph{t} and \emph{t}-1, respectively, followed by the
#' respective stage names, and ending with the estimated reproductive value for
#' that stage within its matrix (\code{rep_value}). The \code{$ahist} element is
#' a data frame showing the reproductive values of the basic stages in the
#' associated stageframe. The reproductive values in this second data frame are
#' estimated via the approach developed in Ehrlen (2000), in which each
#' ahistorical stage's reproductive value is the average of the RVs summed by
#' stage at time \emph{t} weighted by the proportion of that stage pair within
#' the historical stable stage distribution associated with the matrix. Both
#' historical and ahistorical reproductive values are scaled to the first non-
#' zero reproductive value in each case.
#'
#' In addition to the data frames noted above, stochastic analysis will result
#' in the additional output of a list of matrices containing the actual
#' projected reproductive value vectors across all projected times, in the order
#' of population-patch combinations in the \code{lefkoMat} input.
#'
#' @section Notes:
#' For square matrices with fewer than 400 rows, the reproductive value vector
#' is given as the right eigenvector associated with largest real part of all
#' eigenvalues estimated via the \code{eig_gen}() function in the C++ Armadillo
#' library divided by the sum of the associated right eigenvector. For larger
#' matrices, the function assumes that the matrix is sparse and conducts a
#' similar calculation but using the \code{eigs_gen}() for sparse matrix eigen
#' analysis.
#' 
#' In stochastic analysis, the projected mean reproductive value vector is the
#' arithmetic mean across the final projected 1000 times if the simulation is at
#' least 2000 projected times long. If between 500 and 2000 projected times
#' long, then only the final 200 are used, and if fewer than 500 times are used,
#' then all are used. Note that because reproductive values in stochastic
#' simulations can change greatly in the initial portion of the run, we
#' encourage a minimum 2000 projected times per simulation, with 10,000
#' preferred.
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
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 3, 3), stageframe = lathframe, historical = TRUE)
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, 0.40, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
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
  tweights = NA, seed = NA, ...) {
  
  poppatch <- NULL
  
  if (!stochastic) {
    baldrick <- if (any(class(mats$A) == "matrix")) {
      
      if (dim(mats$A)[1] > 400) {
        almost_final <- rv3matrixsp(mats$A)
      } else {
        almost_final <- rv3matrix(mats$A)
      }
      
    } else if (class(mats$A) == "list") {
      
      final <- if (dim(mats$A[[1]])[1] > 400) {
        unlist(lapply(mats$A, function(X) {
            almost_final <- rv3matrixsp(X)
            return(almost_final/almost_final[which(almost_final == (almost_final[which(almost_final > 0)])[1])])
            }
          )
        )
      } else {
        unlist(lapply(mats$A, function(X) {
            almost_final <- rv3matrix(X)
            return(almost_final/almost_final[which(almost_final == (almost_final[which(almost_final > 0)])[1])])
            }
          )
        )
      }
      
    } else {
      
      stop("Input not recognized.")
    }
  } else {
    
    if (!is.na(seed)) {
      set.seed(seed)
    }
    
    mats$labels$poppatch <- paste(mats$labels$pop, mats$labels$patch)
    used_poppatches <- as.list(unique(mats$labels$poppatch))
    
    # Here we get the full stage distribution and reproductive value vector series for all times, as a list
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
      starter <- if (dim(mats$A[[1]])[1] > 400) {
        ss3matrixsp(mats$A[[used_slots[1]]])
      } else {
        ss3matrix(mats$A[[used_slots[1]]])
      }
      
      theseventhmatrix <- proj3(starter, mats$A, theprophecy, 1, 0, 0)
      
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
      output <- cbind.data.frame(modlabels, baldrick[(mat_dims+1):(2*mat_dims)])
      if (is.element("age", names(labels))) {
        names(output) <- c("poppatch", "age", "stage_id", "stage", "agestage_id", "agestage", "rep_value")
      } else {
        names(output) <- c("poppatch", "stage_id", "stage", "rep_value")
      }
    }
    
    rownames(output) <- c(1:dim(output)[1])
    
  } else {
    
    # This section translates historical results to ahistorical and then cleans everything up
    if (!stochastic) {
      ss3 <- stablestage3.lefkoMat(mats) #stablestage3.lefkoMat
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
      
      outputh <- rhist[,c("matrix", "stage_2", "stage_1", "stage_id_2", "stage_id_1", "rep_value")]
      
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
          sum(ss_unlisted[(which(mats$hstages$stage_id_2 == Y)) + ((X - 1) * dim(mats$hstages)[1]), 1])
        })
        return(sum_vecs)
      })
      
      ahistsize <- length(mats$ahstages$stage_id)
      histsize <- length(mats$hstages$stage_id_2)
      rahist <- cbind.data.frame(poppatch = (rep(unlist(used_poppatches), each = ahistsize, times = multiplier)), 
        stage_id = rep(mats$ahstages$stage_id, times = multiplier), stage = rep(mats$ahstages$stage, times = multiplier))
      rhist <- cbind.data.frame(poppatch = (rep(unlist(used_poppatches), each = histsize, times = multiplier)),
        stage_id_2 = rep(mats$hstages$stage_id_2, times = multiplier),
        stage_id_1 = rep(mats$hstages$stage_id_1, times = multiplier),
        stage_2 = rep(mats$hstages$stage_2, times = multiplier),
        stage_1 = rep(mats$hstages$stage_1, times = multiplier))
      
      rhist$ss_prop <- ss_unlisted
      rhist$ss3sum <- ss_sums
      rhist$sscorr <- rhist$ss_prop / rhist$ss3sum
      rhist$sscorr[which(is.na(rhist$sscorr))] <- 0
      rhist$rep_value <- rv_unlisted
      
      rhist$rv3raw <- rhist$sscorr * rhist$rep_value
      
      outputh <- rhist[,c("poppatch", "stage_2", "stage_1", "stage_id_2", "stage_id_1", "rep_value")]
      
      ahlabels <- mats$ahstages[,c("stage_id", "stage")]
      rv2 <- Re(c(apply(as.matrix(unlist(used_poppatches)), 1, function(X) {
        rightset <- subset(rhist, poppatch == X)
        apply(as.matrix(ahlabels[,1]), 1, function(Y) {
          sum(rightset$rv3raw[which(rightset$stage_id_2 == Y)])
        })
      })))
      outputah <- cbind.data.frame(rep(c(1:multiplier), each = length(ahlabels[,1])),
        rep(ahlabels[,1], multiplier), rep(ahlabels[,2], multiplier), rv2)
      names(outputah) <- c("poppatch", "stage_id", "stage", "rep_value_unc")
      
      outputah$rep_value <- apply(as.matrix(c(1:dim(outputah)[1])), 1, function(X) {
        matsub <- subset(outputah, poppatch == outputah$poppatch[X])
        entrystage <- min(which(abs(matsub$rep_value_unc) > 0))
        return(outputah$rep_value_unc[X] / matsub$rep_value_unc[entrystage])
      })
      
      outputah <- outputah[,c("poppatch", "stage_id", "stage", "rep_value")]
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
#' @param ... Other parameters.
#' 
#' @return This function returns a vector data frame characterizing the 
#' reproductive values for stages of a population projection matrix. This is 
#' given as the left eigenvector associated with largest real part of the
#' dominant eigenvalue, divided by the first non-zero element of the left 
#' eigenvector. Eigen analysis is handled by the \code{eig_gen}() function in
#' the C++ Armadillo library for square matrices with fewer than 400 rows, and
#' using the \code{eigs_gen}() function for larger matrices.
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
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 3, 3), stageframe = lathframe, historical = TRUE)
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, 0.40, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
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
repvalue3.matrix <- function(mats, ...)
{
  if (dim(mats)[1] > 400) {
    v <- rv3matrixsp(mats)
  } else {
    v <- rv3matrix(mats)
  }
  
  return(v)
}

#' Estimate Sensitivity of Population Growth Rate to Matrix Elements
#' 
#' \code{sensitivity3()} is a generic function that returns the sensitivity of
#' the population growth rate to the elements of the matrices in a matrix
#' population model. Currently, this function estimates both deterministic and
#' stochastic sensitivities, where the growth rate is \eqn{\lambda} in the former
#' case and the log of the stochastic \eqn{\lambda} in the latter case. This
#' function is made to handle very large and sparse matrices supplied as
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
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 3, 3), stageframe = lathframe, historical = TRUE)
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, 0.40, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
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
#' \code{lefkoMat}. If deterministic, then \eqn{\lambda} is taken as the population
#' growth rate. If stochastic, then the log of stochastic \eqn{\lambda}, or
#' the log stochastic growth rate, is taken as the population growth rate. This
#' function can handle large and sparse matrices, and so can be used with large
#' historical matrices, IPMs, age x stage matrices, as well as smaller
#' ahistorical matrices.
#' 
#' @param mats An object of class \code{lefkoMat}.
#' @param stochastic A logical value determining whether to conduct a
#' deterministic (FALSE) or stochastic (TRUE) sensitivity analysis. Defaults to
#' FALSE.
#' @param steps The number of times to project forward in stochastic simulation.
#' Defaults to 10,000.
#' @param time_weights Numeric vector denoting the probabilistic weightings of
#' annual matrices. Defaults to equal weighting among times.
#' @param ... Other parameters.
#' 
#' @return This function returns an object of class \code{lefkoSens}, which is a
#' list of 7 elements. The first, h_sensmats, is a list of historical sensitivity
#' matrices (NULL if an ahMPM is used as input). The second, ah_elasmats, is a
#' list of either ahistorical sensitivity matrices if an ahMPM is used as input,
#' or, if an hMPM is used as input, then the result is a list of ahistorical
#' matrices based on the equivalent historical dependencies assumed in the
#' input historical matrices. The third element, h_stages, is a data frame
#' showing historical stage pairs (NULL if ahMPM used as input). The fourth
#' element, ah_stages, is a data frame showing the order of ahistorical stages.
#' The last 3 elements are the A, U, and F portions of the input.
#' 
#' @section Notes:
#' Deterministic sensitivities are estimated as eqn. 9.14 in Caswell (2001,
#' Matrix Population Models). Stochastic sensitivities are estimated as eqn.
#' 14.97 in Caswell (2001). Note that stochastic sensitivities are of the log of
#' the stochastic \eqn{\lambda}.
#'
#' Currently, this function does not estimate equivalent ahistorical stochastic
#' sensitivities for input historical matrices.
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
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 3, 3), stageframe = lathframe, historical = TRUE)
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, 0.40, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
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
  time_weights = NA, ...) {
  
  if (stochastic & length(mats$A) < 2) {
    stop("Stochastic sensitivity estimation cannot be completed with fewer than 2 annual matrices.", call. = FALSE)
  }
  
  if (!stochastic) {
    # Deterministic sensitivity analysis
    
    baldrick <- if (any(class(mats$A) == "matrix")) {
      
      if (dim(mats$A)[1] > 400 & ((length(mats$A[which(mats$A[[1]] > 0)]) / length(mats$A)) <= 0.5)) {
        sens3matrixsp(mats$A)
      } else {
        sens3matrix(mats$A)
      }
      
    } else if (class(mats$A) == "list") {
      
      if (all(is.na(mats$hstages))) {
        
        if (dim(mats$A[[1]])[1] > 400 & ((length(mats$A[[1]][which(mats$A[[1]] > 0)]) / length(mats$A[[1]])) <= 0.5)) {
          lapply(mats$A, sens3matrixsp)
        } else {
          lapply(mats$A, sens3matrix)
        }
        
      } else {
        
        lapply(mats$A, sens3hlefko, mats$ahstages, mats$hstages)
      }
    } else {
      
      stop("Input not recognized.")
    }
    
    if (all(is.na(mats$hstages))) {
      
      ahlabels <- mats$ahstages
      
      output <- list(h_sensmats = NULL, ah_sensmats = baldrick, h_stages = NULL, 
        ah_stages = ahlabels, A = mats$A, U = mats$U, F = mats$F)
      
    } else {
      
      he_list <- lapply(baldrick, function(X) {X$h_smat})
      ahe_list <- lapply(baldrick, function(X) {X$ah_smat})
      
      hlabels <- mats$hstages
      
      ahlabels <- mats$ahstages
      
      output <- list(h_sensmats = he_list, ah_sensmats = ahe_list,
        h_stages = hlabels, ah_stages = ahlabels, A = mats$A, U = mats$U,
        F = mats$F)
    }
  } else {
    # Stochastic sensitivity analysis
    
    if(!any(is.na(time_weights))) {
      returned_cube <- stoch_senselas(mats, times = steps, style = 1,
        tweights = time_weights) 
    } else {
      returned_cube <- stoch_senselas(mats, times = steps, style = 1) 
    }
    
    returned_list <- lapply(as.list(c(1:dim(returned_cube)[3])), function(X) {
        return(returned_cube[,,X])
      }
    )
    
    if (!all(is.na(mats$hstages))) {
      output <- list(h_sensmats = returned_list, ah_sensmats = NULL,
        h_stages = mats$hstages, ah_stages = mats$ahstages, A = mats$A,
        U = mats$U, F = mats$F)
    } else {
      output <- list(h_sensmats = NULL, ah_sensmats = returned_list,
        h_stages = mats$hstages, ah_stages = mats$ahstages, A = mats$A,
        U = mats$U, F = mats$F)
    }
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
#' @param ... Other parameters.
#' 
#' @return This function returns a single sensitivity matrix.
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
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 3, 3), stageframe = lathframe, historical = TRUE)
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, 0.40, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
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
sensitivity3.matrix <- function(mats, ...)
{
  if (dim(mats)[1] > 400 & ((length(mats[which(mats > 0)]) / length(mats)) <= 0.5)) {
    wcorr <- sens3matrixsp(mats)
  } else {
    wcorr <- sens3matrix(mats)
  }
  
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
#' @param steps The number of times to project forward in stochastic simulation.
#' Defaults to 10,000.
#' @param time_weights Numeric vector denoting the probabilistic weightings of
#' annual matrices. Defaults to equal weighting among times.
#' @param ... Other parameters.
#' 
#' @return This function returns an object of class \code{lefkoSens}, which is a
#' list of 7 elements. The first, h_sensmats, is a list of historical
#' sensitivity matrices (NULL if an ahMPM is used as input). The second,
#' ah_elasmats, is a list of ahistorical sensitivity matrices if an ahMPM is
#' used as input (NULL if an hMPM is used as input). The third element,
#' h_stages, and the fourth element, ah_stages, are NULL. The last 3 elements
#' include the original A matrices supplied (as the A element), followed by
#' NULLs for the U and F elements.
#' 
#' @section Notes:
#' Deterministic sensitivities are estimated as eqn. 9.14 in Caswell (2001,
#' Matrix Population Models). Stochastic sensitivities are estimated as eqn.
#' 14.97 in Caswell (2001). Note that stochastic sensitivities are with regard
#' to the log of the stochastic \eqn{\lambda}.
#'
#' Determination of whether the input MPM is ahistorical or historical is made
#' on the basis of the dimensionality of input matrices - fewer than 400 rows
#' suggests an ahistorical MPM, while more indicates a historical MPM. This
#' designation does not impact the underlying algorithms in any way, so mistaken
#' calls can be fixed by moving the resulting sensitivity matrices to the
#' correct element of the output list.
#' 
#' Currently, this function does not estimate equivalent ahistorical stochastic
#' sensitivities for input historical matrices.
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
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 3, 3), stageframe = lathframe, historical = TRUE)
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, 0.40, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
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
  time_weights = NA, ...) {
  
  if(length(setdiff(unlist(lapply(mats, class)), c("matrix", "array"))) > 0) {
    stop("Input list must be composed only of numeric matrices.", call. = FALSE)
  }
  if (stochastic & length(mats) < 2) {
    stop("Stochastic sensitivity estimation cannot be completed with fewer than 2 annual matrices.", call. = FALSE)
  }
  
  if (!stochastic) {
    # Deterministic sensitivity analysis
    
    if (dim(mats[[1]])[1] > 400 & ((length(mats[[1]][which(mats[[1]] > 0)]) / length(mats[[1]])) <= 0.5)) {
      message("Matrices have more than 400 rows and columns, so assuming they are historical.")
      
      baldrick <- lapply(mats, sens3matrixsp)
      
      output <- list(h_sensmats = baldrick, ah_sensmats = NULL, h_stages = NULL,
        ah_stages = NULL, A = mats, U = NULL, F = NULL)
    } else {
      message("Matrices have fewer than 400 rows and columns, so assuming they are ahistorical.")
      
      baldrick <- lapply(mats, sens3matrix)
      
      output <- list(h_sensmats = NULL, ah_sensmats = baldrick, h_stages = NULL,
        ah_stages = NULL, A = mats, U = NULL, F = NULL)
    }
    
  } else {
    # Stochastic sensitivity analysis
    
    if(!any(is.na(time_weights))) {
      returned_cube <- stoch_senselas(mats, times = steps, style = 1,
        tweights = time_weights) 
    } else {
      returned_cube <- stoch_senselas(mats, times = steps, style = 1) 
    }
    
    returned_list <- lapply(as.list(c(1:dim(returned_cube)[3])), function(X) {
        return(returned_cube[,,X])
      }
    )
    
    if (dim(mats[[1]])[1] > 400) {
      message("Matrices have more than 400 rows and columns, so assuming they are historical.")
      
      output <- list(h_sensmats = returned_list, ah_sensmats = NULL,
        h_stages = NULL, ah_stages = NULL, A = mats, U = NULL, F = NULL)
    } else {
      message("Matrices have fewer than 400 rows and columns, so assuming they are ahistorical.")
      
      output <- list(h_sensmats = NULL, ah_sensmats = returned_list,
        h_stages = NULL, ah_stages = NULL, A = mats, U = NULL, F = NULL)
    }
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
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 3, 3), stageframe = lathframe, historical = TRUE)
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, 0.40, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
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
#' @param steps The number of times to project forward in stochastic simulation.
#' Defaults to 10,000.
#' @param time_weights Numeric vector denoting the probabilistic weightings of
#' annual matrices. Defaults to equal weighting among times.
#' @param ... Other parameters.
#' 
#' @return This function returns an object of class \code{lefkoElas}, which is a
#' list with 7 elements. The first, h_elasmats, is a list of historical
#' elasticity matrices (NULL if an ahMPM is used as input). The second,
#' ah_elasmats, is a list of either ahistorical elasticity matrices if an ahMPM
#' is used as input, or, if an hMPM is used as input, then the result is a list
#' of elasticity matrices in which historical elasticities have been summed by
#' the stage in times \emph{t} and \emph{t}+1 to produce historically-corrected
#' elasticity matrices, which are equivalent in dimension to ahistorical
#' elasticitymatrices but reflect the effects of stage in time t-1. The third
#' element, h_stages, is a data frame showing historical stage pairs (NULL if
#' ahMPM used as input). The fourth element, ah_stages, is a data frame showing
#' the order of ahistorical stages. The last 3 elements are the A, U, and F
#' portions of the input.
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
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 3, 3), stageframe = lathframe, historical = TRUE)
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, 0.40, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
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
  time_weights = NA, ...) {
  
  if (!stochastic) {
    # Deterministic elasticity analysis
    
    baldrick <- if (any(class(mats$A) == "matrix")) {
      
      if (dim(mats$A)[1] > 400 & ((length(mats$A[which(mats$A[[1]] > 0)]) / length(mats$A)) <= 0.5)) {
        elas3matrixsp(mats$A)
      } else {
        elas3matrix(mats$A)
      }
      
    } else if (class(mats$A) == "list") {
      
      if (all(is.na(mats$hstages))) {
        
        if (dim(mats$A[[1]])[1] > 400 & ((length(mats$A[[1]][which(mats$A[[1]] > 0)]) / length(mats$A[[1]])) <= 0.5)) {
          lapply(mats$A, elas3matrixsp)
        } else {
          lapply(mats$A, elas3matrix)
        }
        
      } else {
        
        lapply(mats$A, elas3hlefko, mats$ahstages, mats$hstages)
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
        ah_stages = ahlabels, A = mats$A, U = mats$U, F = mats$F)
    } else {
      
      he_list <- lapply(baldrick, function(X) {X$h_emat})
      ahe_list <- lapply(baldrick, function(X) {X$ah_emat})
      
      hlabels <- mats$hstages
      
      ahlabels <- mats$ahstages #Originally only the first two columns
      
      output <- list(h_elasmats = he_list, ah_elasmats = ahe_list, 
        h_stages = hlabels, ah_stages = ahlabels, A = mats$A, U = mats$U,
        F = mats$F)
    }
  } else {
    # Stochastic elasticity analysis
    
    if(!any(is.na(time_weights))) {
      returned_cube <- stoch_senselas(mats, times = steps, style = 2,
        tweights = time_weights) 
    } else {
      returned_cube <- stoch_senselas(mats, times = steps, style = 2) 
    }
    
    returned_list <- lapply(as.list(c(1:dim(returned_cube)[3])), function(X) {
      return(returned_cube[,,X])
      }
    )
    
    if (!all(is.na(mats$hstages))) {
      ah_list <- lapply(returned_list, shopliftersunite, mats$hstages, mats$ahstages)
      
      output <- list(h_elasmats = returned_list, ah_elasmats = ah_list,
        h_stages = mats$hstages, ah_stages = mats$ahstages, A = mats$A,
        U = mats$U, F = mats$F)
    } else {
      output <- list(h_elasmats = NULL, ah_elasmats = returned_list,
        h_stages = mats$hstages, ah_stages = mats$ahstages, A = mats$A,
        U = mats$U, F = mats$F)
    }
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
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 3, 3), stageframe = lathframe, historical = TRUE)
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, 0.40, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
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
elasticity3.matrix <- function(mats, ...)
{
  if (dim(mats)[1] > 400 & ((length(mats[which(mats > 0)]) / length(mats)) <= 0.5)) {
    wcorr <- elas3matrixsp(mats)
  } else {
    wcorr <- elas3matrix(mats)
  }
  
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
#' @param steps The number of times to project forward in stochastic simulation.
#' Defaults to 10,000.
#' @param time_weights Numeric vector denoting the probabilistic weightings of
#' annual matrices. Defaults to equal weighting among times.
#' @param ... Other parameters.
#' 
#' @return This function returns an object of class \code{lefkoElas}, which is a
#' list with 7 elements. The first, h_elasmats, is a list of historical
#' elasticity matrices, though in the standard list case it returns a NULL
#' value. The second, ah_elasmats, is a list of either ahistorical elasticity
#' matrices. Since the list used as input contains no stage information, all
#' matrices are assumed to be ahistorical and elasticities are provided via this
#' element. The third element, h_stages, and the fourth element, ah_stages, are
#' NULL. The last 3 elements are the original A matrices in element A, followed
#' by NULL values for the U and F elements.
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
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 3, 3), stageframe = lathframe, historical = TRUE)
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
#' cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'     "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL", "rep",
#'     "rep"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, 0.40, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
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
  time_weights = NA, ...) {
  
  if(length(setdiff(unlist(lapply(mats, class)), c("matrix", "array"))) > 0) {
    stop("Input list must be composed only of numeric matrices.", call. = FALSE)
  }
  if (stochastic & length(mats) < 2) {
    stop("Stochastic elasticity estimation cannot be completed with fewer than 2 annual matrices.", call. = FALSE)
  }
  
  if (!stochastic) {
    # Deterministic elasticity analysis
    
    if (dim(mats[[1]])[1] > 400 & ((length(mats[[1]][which(mats$A[[1]] > 0)]) / length(mats[[1]])) <= 0.5)) {
      baldrick <- lapply(mats, elas3matrixsp)
      
      output <- list(h_elasmats = NULL, ah_elasmats = baldrick, h_stages = NULL,
        ah_stages = NULL, A = mats, U = NULL, F = NULL)
    } else {
      baldrick <- lapply(mats, elas3matrix)
      
      output <- list(h_elasmats = NULL, ah_elasmats = baldrick, h_stages = NULL,
        ah_stages = NULL, A = mats, U = NULL, F = NULL)
    }
  } else {
    # Stochastic elasticity analysis
    
    if(!any(is.na(time_weights))) {
      returned_cube <- stoch_senselas(mats, times = steps, style = 2,
        tweights = time_weights) 
    } else {
      returned_cube <- stoch_senselas(mats, times = steps, style = 2) 
    }
    
    returned_list <- lapply(as.list(c(1:dim(returned_cube)[3])), function(X) {
      return(returned_cube[,,X])
      }
    )
    output <- list(h_elasmats = NULL, ah_elasmats = returned_list,
      h_stages = NULL, ah_stages = NULL, A = mats, U = NULL, F = NULL)
  }
  class(output) <- "lefkoElas"
  
  return(output)
}

#' Summarize lefkoElas Objects
#' 
#' Function \code{summary.lefkoElas()} summarizes \code{lefkoElas} objects.
#' Particularly, it breaks down elasticity values by the kind of ahistorical
#' and, if applicable, historical transition.
#'
#' @param object A \code{lefkoElas} object.
#' @param ... Other parameters.
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
#' lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sd", "Sdl"), 
#'   stage2 = c("Sd", "Sd", "Sd", "rep", "rep"),
#'   stage1 = c("Sd", "rep", "rep", "all", "all"), 
#'   givenrate = c(0.345, 0.345, 0.054, NA, NA),
#'   multiplier = c(NA, NA, NA, 0.345, 0.054),
#'   type = c(1, 1, 1, 3, 3), stageframe = lathframe, historical = TRUE)
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
    indices <- bambi2(elasmats$ah_stages)
  } else {
    indices <- bambi3(elasmats$ah_stages, elasmats$h_stages)
  }
  
  for (i in c(1:used_iterations)) {
    trialguy <- demolition3(used_emats[[i]], elasmats$A[[i]], elasmats$F[[i]], indices)
    
    if (i == 1) {
      hist <- trialguy$hist
      ahist <- trialguy$ahist
      if (num_h_mats > 0) names(hist)[2] <- "matrix1"
      names(ahist)[2] <- "matrix1"
    } else {
      if (num_h_mats > 0) hist <- cbind.data.frame(hist, trialguy$hist[,2])
      ahist <- cbind.data.frame(ahist, trialguy$ahist[,2])
      if (num_h_mats > 0) names(hist)[(i+1)] <- paste0("matrix", i)
      names(ahist)[(i+1)] <- paste0("matrix", i)
    }
  }
  
  output <- list(hist = hist, ahist = ahist)
  
  return (output)
}

