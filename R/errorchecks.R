#' Summary of Class "lefkoCondMat"
#' 
#' This function provides basic information summarizing the characteristics of
#' conditional matrices derived from a \code{lefkoMat} object.
#' 
#' @param object An object of class \code{lefkoMat}.
#' @param ... Other parameters.
#' 
#' @return A summary of the object, showing the number of historical matrices,
#' as well as the number of conditional matrices nested within each historical
#' matrix.
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
#' cypsupp3r <- supplemental(stage3 = c("SD", "P1", "SD", "P1", "P2", "P2",
#'   "P3", "SL", "SL", "SL", "D", "XSm", "Sm", "D", "XSm", "Sm", "SD", "P1"),
#'   stage2 = c("SD", "SD", "SD", "SD", "P1", "P1", "P2", "P3", "SL", "SL",
#'     "SL", "SL", "SL", "SL", "SL", "SL", "rep", "rep"),
#'   stage1 = c("SD", "rep", "SD", "rep", "SD", "rep", "P1", "P2", "P3", "SL",
#'     "SL", "SL", "SL", "P3", "P3", "P3", "mat", "mat"),
#'   eststage3 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm",
#'     "D", "XSm", "Sm", NA, NA),
#'   eststage2 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm",
#'     "XSm", "XSm", "XSm", NA, NA),
#'   eststage1 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm",
#'     "XSm", "XSm", "XSm", NA, NA),
#'   givenrate = c(0.10, 0.20, 0.10, 0.20, 0.20, 0.20, 0.20, 0.25, 0.40, 0.40,
#'     NA, NA, NA, NA, NA, NA, NA, NA),
#'   multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
#'     NA, 0.5, 0.5),
#'   type =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1, 1, 1, 3, 3),
#'   stageframe = cypframe_raw, historical = TRUE)
#' 
#' cypmatrix3r <- rlefko3(data = cypraw_v1, stageframe = cypframe_raw, 
#'   year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
#'   size = c("size3added", "size2added", "size1added"), 
#'   supplement = cypsupp3r, yearcol = "year2", patchcol = "patchid",
#'   indivcol = "individ")
#' 
#' cypcondmats <- cond_hmpm(cypmatrix3r)
#' 
#' summary(cypcondmats)
#' 
#' @export
summary.lefkoCondMat <- function(object, ...) {
  
  histmatrices <- object$Acond
  condmatrices <- histmatrices[[1]]
  firstcondmat <- condmatrices[[1]]
  
  numhistmats <- length(histmatrices)
  prevstages <- length(condmatrices)
  matdim <- dim(firstcondmat)
  
  writeLines(paste0("\nThis lefkoCondMat object contains ", prevstages, " conditional matrices per historical matrix, covering ", numhistmats, " historical matrices."))
  writeLines(paste0("Each conditional matrix is a square matrix with ", matdim[1], " rows and columns, and a total of ", matdim[1]*matdim[1], " elements."))
  writeLines(paste0("\nThe order of conditional matrices corresponding to stage in time t-1 is:\n", paste(object$ahstages$stage, collapse = " ")))
  writeLines("\nThe order of historical matrices is: \n")
  print.data.frame(object$labels)
  
  writeLines("\nThe order of conditional matrices matches the stage column in object $ahstages.")
  writeLines("The order of historical matrices follows that shown in object $labels.")
  
  return()
}
