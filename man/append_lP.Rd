% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{append_lP}
\alias{append_lP}
\title{Append Projections Into New lefkoProj Object}
\usage{
append_lP(proj1 = NULL, proj2 = NULL)
}
\arguments{
\item{proj1}{A \code{lefkoProj} object.}

\item{proj2}{A second \code{lefkoProj} object, based on the same stageframe
as \code{proj1}.}
}
\value{
A list of class \code{lefkoProj}, which always includes the first
three elements of the following, and also includes the remaining elements
below when a \code{lefkoMat} object is used as input:
\item{projection}{A list of lists of matrices showing the total number of
individuals per stage per occasion. The first list corresponds to each
pop-patch followed by each population (this top-level list is a single
element in \code{f_projection3()}). The inner list corresponds to
replicates within each pop-patch or population.}
\item{stage_dist}{A list of lists of the actual stage distribution in each
occasion in each replicate in each pop-patch or population.}
\item{rep_value}{A list of lists of the actual reproductive value in each
occasion in each replicate in each pop-patch or population.}
\item{pop_size}{A list of matrices showing the total population size in each
occasion per replicate (row within data frame) per pop-patch or population
(list element). \code{NA} values will result if projections with different
numbers of time steps are appended.}
\item{labels}{A data frame showing the order of populations and patches in
item \code{projection}.}
\item{ahstages}{The original stageframe used in the study.}
\item{hstages}{A data frame showing the order of historical stage pairs.}
\item{agestages}{A data frame showing the order of age-stage pairs.}
\item{labels}{A short data frame indicating the population (always \code{1}),
and patch (either the numeric index of the single chosen patch, or \code{1}
in all other cases). Any pop-patches having the same designation across the
two input projections will be appended together.}
\item{control}{A data frame showing the number of replicates and time steps
corresponding to each set of projections, where each set corresponds to a
pop-patch within the labels object of each input projection.}
\item{density}{The data frame input under the density option. Only provided
if input by the user for at least one of the two projections. Output as a
nested list corresponding to each pop-patch - replicate.}
\item{density_vr}{The data frame input under the density_vr option. Only
provided if input by the user for at least one of the two projections.
Output as a nested list corresponding to each pop-patch - replicate.}
}
\description{
Function \code{append_lP()} combines two population projections. It takes
two \code{lefkoProj} objects and appends them into a new \code{lefkoProj}
object.
}
\section{Notes}{

\code{lefkoProj} objects resulting from previous appends can also be
appended.
}

\examples{
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
  stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "rep", "rep"),
  eststage3 = c(NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
  eststage2 = c(NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
  givenrate = c(0.1, 0.2, 0.2, 0.2, 0.25, NA, NA, NA, NA, NA),
  multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
  type = c(1, 1, 1, 1, 1, 1, 1, 1, 3, 3), stageframe = cypframe_raw,
  historical = FALSE)

cypmatrix2r_AB <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
  year = "all", patch = c("A", "B"), stages = c("stage3", "stage2"),
  size = c("size3added", "size2added"), supplement = cypsupp2r,
  yearcol = "year2",  patchcol = "patchid", indivcol = "individ")

cypmatrix2r_AC <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
  year = "all", patch = c("A", "C"), stages = c("stage3", "stage2"),
  size = c("size3added", "size2added"), supplement = cypsupp2r,
  yearcol = "year2",  patchcol = "patchid", indivcol = "individ")

cypproj1 <- projection3(cypmatrix2r_AB, nreps = 5, times = 15,
  stochastic = TRUE)
cypproj2 <- projection3(cypmatrix2r_AC, nreps = 10, times = 20,
  stochastic = TRUE)
cypproj3 <- append_lP(cypproj1, cypproj2)

}
\seealso{
\code{\link{projection3}()}
}
