% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/popdyn.R
\name{sensitivity3.list}
\alias{sensitivity3.list}
\title{Estimate Sensitivity of Population Growth Rate of a List of Matrices}
\usage{
\method{sensitivity3}{list}(
  mats,
  stochastic = FALSE,
  times = 10000,
  tweights = NA,
  historical = FALSE,
  seed = NA,
  sparse = "auto",
  append_mats = FALSE,
  ...
)
}
\arguments{
\item{mats}{An object of class \code{matrix}.}

\item{stochastic}{A logical value determining whether to conduct a
deterministic (FALSE) or stochastic (TRUE) sensitivity analysis. Defaults to
FALSE.}

\item{times}{The number of occasions to project forward in stochastic
simulation. Defaults to 10,000.}

\item{tweights}{An optional numeric vector or matrix denoting the
probabilities of choosing each matrix in a stochastic projection. If a matrix
is input, then a first-order Markovian environment is assumed, in which the
probability of choosing a specific annual matrix depends on which annual
matrix is currently chosen. If a vector is input, then the choice of annual
matrix is assumed to be independent of the current matrix. Defaults to equal
weighting among matrices.}

\item{historical}{A logical value indicating whether matrices are historical.
Defaults to \code{FALSE}.}

\item{seed}{A number to use as a random number seed in stochastic projection.}

\item{sparse}{A text string indicating whether to use sparse matrix encoding
(\code{"yes"}) or dense matrix encoding (\code{"no"}). Defaults to
\code{"auto"}, in which case sparse matrix encoding is used with square
matrices with at least 50 rows and no more than 50\% of elements with values
greater than zero.}

\item{append_mats}{A logical value indicating whether to include the original
matrices input as object \code{mats} in the output \code{lefkoSense} object.
Defaults to FALSE.}

\item{...}{Other parameters.}
}
\value{
This function returns an object of class \code{lefkoSens}, which is a
list of 8 elements. The first, \code{h_sensmats}, is a list of historical
sensitivity matrices (\code{NULL} if an ahMPM is used as input). The second,
\code{ah_elasmats}, is a list of ahistorical sensitivity matrices if an ahMPM
is used as input (\code{NULL} if an hMPM is used as input). The third
element, \code{hstages}, the fourth element, \code{agestages}, and the fifth
element, \code{ahstages}, are \code{NULL}. The last 3 elements include the
original A matrices supplied (as the \code{A} element), followed by
\code{NULL}s for the U and F elements.
}
\description{
\code{sensitivity3.list()} returns the sensitivities of population growth
rate to elements of matrices supplied in a list. The sensitivity analysis can
be deterministic or stochastic, but if the latter then at least two A
matrices must be included in the list. This function can handle large and
sparse matrices, and so can be used with large historical matrices, IPMs,
age x stage matrices, as well as smaller ahistorical matrices.
}
\section{Notes}{

All sensitivity matrix outputs from this function are in standard matrix
format.

Deterministic sensitivities are estimated as eqn. 9.14 in Caswell (2001,
Matrix Population Models). Stochastic sensitivities are estimated as eqn.
14.97 in Caswell (2001). Note that stochastic sensitivities are with regard
to the log of the stochastic \eqn{\lambda}.

Currently, this function does not estimate equivalent ahistorical stochastic
sensitivities for input historical matrices, due to the lack of guidance
input on the order of stages (guidance is provided within \code{lefkoMat}
objects).

Speed can sometimes be increased by shifting from automatic sparse matrix
determination to forced dense or sparse matrix projection. This will most
likely occur when matrices have between 30 and 300 rows and columns.
Defaults work best when matrices are very small and dense, or very large and
sparse.

The \code{time_weights} and \code{steps} arguments are now deprecated.
Instead, please use the \code{tweights} and \code{times} arguments.
}

\examples{
# Lathyrus example
data(lathyrus)

sizevector <- c(0, 100, 13, 127, 3730, 3800, 0)
stagevector <- c("Sd", "Sdl", "VSm", "Sm", "VLa", "Flo", "Dorm")
repvector <- c(0, 0, 0, 0, 0, 1, 0)
obsvector <- c(0, 1, 1, 1, 1, 1, 0)
matvector <- c(0, 0, 1, 1, 1, 1, 1)
immvector <- c(1, 1, 0, 0, 0, 0, 0)
propvector <- c(1, 0, 0, 0, 0, 0, 0)
indataset <- c(0, 1, 1, 1, 1, 1, 1)
binvec <- c(0, 100, 11, 103, 3500, 3800, 0.5)

lathframe <- sf_create(sizes = sizevector, stagenames = stagevector,
  repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
  immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
  propstatus = propvector)

lathvert <- verticalize3(lathyrus, noyears = 4, firstyear = 1988,
  patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9,
  juvcol = "Seedling1988", sizeacol = "Volume88", repstracol = "FCODE88",
  fecacol = "Intactseed88", deadacol = "Dead1988",
  nonobsacol = "Dormant1988", stageassign = lathframe, stagesize = "sizea",
  censorcol = "Missing1988", censorkeep = NA, censor = TRUE)

lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl", "mat"),
  stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep", "Sdl"),
  stage1 = c("Sd", "rep", "Sd", "rep", "npr", "npr", "Sd"),
  eststage3 = c(NA, NA, NA, NA, NA, NA, "mat"),
  eststage2 = c(NA, NA, NA, NA, NA, NA, "Sdl"),
  eststage1 = c(NA, NA, NA, NA, NA, NA, "NotAlive"),
  givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA, NA),
  multiplier = c(NA, NA, NA, NA, 0.345, 0.054, NA),
  type = c(1, 1, 1, 1, 3, 3, 1), type_t12 = c(1, 2, 1, 2, 1, 1, 1),
  stageframe = lathframe, historical = TRUE)

ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe, year = "all", 
  stages = c("stage3", "stage2", "stage1"), supplement = lathsupp3,
  yearcol = "year2", indivcol = "individ")

sensitivity3(ehrlen3$A)

# Cypripedium example
data(cypdata)

sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg",
  "XLg")
repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)

cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector,
  repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
  propstatus = propvector, immstatus = immvector, indataset = indataset,
  binhalfwidth = binvec)

cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
  patchidcol = "patch", individcol = "plantid", blocksize = 4,
  sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
  repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
  stageassign = cypframe_raw, stagesize = "sizeadded", NAas0 = TRUE,
  NRasRep = TRUE)

cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "D", 
    "XSm", "Sm", "SD", "P1"),
  stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep",
    "rep"),
  eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
  eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
  givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, NA, NA, NA, NA, NA),
  multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
  type =c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
  stageframe = cypframe_raw, historical = FALSE)

cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
  year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
  size = c("size3added", "size2added"), supplement = cypsupp2r,
  yearcol = "year2", patchcol = "patchid", indivcol = "individ")

sensitivity3(cypmatrix2r$A)

}
\seealso{
\code{\link{sensitivity3}()}

\code{\link{sensitivity3.lefkoMat}()}

\code{\link{sensitivity3.matrix}()}

\code{\link{sensitivity3.dgCMatrix}()}
}
