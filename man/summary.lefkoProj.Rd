% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/popdyn.R
\name{summary.lefkoProj}
\alias{summary.lefkoProj}
\title{Summarize lefkoProj Objects}
\usage{
\method{summary}{lefkoProj}(
  object,
  threshold = 1,
  inf_alive = TRUE,
  milepost = c(0, 0.25, 0.5, 0.75, 1),
  ext_time = FALSE,
  ...
)
}
\arguments{
\item{object}{A \code{lefkoProj} object.}

\item{threshold}{A threshold population size to be searched for in
projections. Defaults to 1.}

\item{inf_alive}{A logical value indicating whether to treat infinitely
large population size as indicating that the population is still extant.
If \code{FALSE}, then the population is considered extinct. Defaults to
\code{TRUE}.}

\item{milepost}{A numeric vector indicating at which points in the projection
to assess detailed results. Can be input as integer values, in which case
each number must be between 1 and the total number of occasions projected in
each projection, or decimals between 0 and 1, which would then be translated
into the corresponding projection steps of the total. Defaults to
\code{c(0, 0.25, 0.50, 0.75, 1.00)}.}

\item{ext_time}{A logical value indicating whether to output extinction times
per population-patch. Defaults to \code{FALSE}.}

\item{...}{Other parameters currently not utilized.}
}
\value{
Apart from a statement of the results, this function outputs a list
with the following elements:
\item{milepost_sums}{A data frame showing the number of replicates at each
of the milepost times that is above the threshold population/patch size.}
\item{extinction_times}{A dataframe showing the numbers of replicates going
extinct (\code{ext_reps}) and mean extinction time (\code{ext_time}) per
population-patch. If \code{ext_time = FALSE}, then only outputs \code{NA}.}
}
\description{
Function \code{summary.lefkoProj()} summarizes \code{lefkoProj} objects.
Particularly, it breaks down the data frames provided in the 
\code{projection} element in ways meaningful for those running simulations.
}
\section{Notes}{

The \code{inf_alive} and \code{ext_time} options both assess whether
replicates have reached a value of \code{NaN} or \code{Inf}. If
\code{inf_alive = TRUE} or \code{ext_time = TRUE} and one of these values is
found, then the replicate is counted in the \code{milepost_sums} object if
the last numeric value in the replicate is above the \code{threshold} value,
and is counted as extant and not extinct if the last numeric value in the
replicate is above the extinction threshold of a single individual.

Extinction time is calculated on the basis of whether the replicate ever
falls below a single individual. A replicate with a positive population size
below 0.0 that manages to rise above 1.0 individual is still considered to
have gone extinct the first time it crossed below 1.0.

If the input \code{lefkoProj} object is a mixture of two or more other
\code{lefkoProj} objects, then mileposts will be given relative to the
maximum number of time steps noted.
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

lathrepm <- matrix(0, 7, 7)
lathrepm[1, 6] <- 0.345
lathrepm[2, 6] <- 0.054

lathsupp3 <- supplemental(stage3 = c("Sd", "Sd", "Sdl", "Sdl", "Sd", "Sdl"), 
  stage2 = c("Sd", "Sd", "Sd", "Sd", "rep", "rep"),
  stage1 = c("Sd", "rep", "Sd", "rep", "all", "all"), 
  givenrate = c(0.345, 0.345, 0.054, 0.054, NA, NA),
  multiplier = c(NA, NA, NA, NA, 0.345, 0.054),
  type = c(1, 1, 1, 1, 3, 3), type_t12 = c(1, 2, 1, 2, 1, 1),
  stageframe = lathframe, historical = TRUE)

ehrlen3 <- rlefko3(data = lathvert, stageframe = lathframe,
  year = c(1989, 1990), stages = c("stage3", "stage2", "stage1"),
  repmatrix = lathrepm, supplement = lathsupp3, yearcol = "year2",
  indivcol = "individ")

lathproj <- projection3(ehrlen3, nreps = 5, stochastic = TRUE)
summary(lathproj)

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

cypsupp3r <- supplemental(stage3 = c("SD", "SD", "P1", "P1", "P2", "P3", "SL",
    "D", "XSm", "Sm", "D", "XSm", "Sm", "mat", "mat", "mat", "SD", "P1"),
  stage2 = c("SD", "SD", "SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL",
    "SL", "SL", "D", "XSm", "Sm", "rep", "rep"),
  stage1 = c("SD", "rep", "SD", "rep", "SD", "P1", "P2", "P3", "P3", "P3",
    "SL", "SL", "SL", "SL", "SL", "SL", "mat", "mat"),
  eststage3 = c(NA, NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", "D", "XSm", "Sm",
    "mat", "mat", "mat", NA, NA),
  eststage2 = c(NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm", "XSm",
    "XSm", "D", "XSm", "Sm", NA, NA),
  eststage1 = c(NA, NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm", "XSm",
    "XSm", "XSm", "XSm", "XSm", NA, NA),
  givenrate = c(0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.25, NA, NA, NA, NA, NA, NA,
    NA, NA, NA, NA, NA),
  multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
    NA, 0.5, 0.5),
  type = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
  type_t12 = c(1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
  stageframe = cypframe_raw, historical = TRUE)

cypmatrix3r <- rlefko3(data = cypraw_v1, stageframe = cypframe_raw, 
  year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
  size = c("size3added", "size2added", "size1added"), 
  supplement = cypsupp3r, yearcol = "year2", 
  patchcol = "patchid", indivcol = "individ")

cypstoch <- projection3(cypmatrix3r, nreps = 5, stochastic = TRUE)
summary(cypstoch, ext_time = TRUE)

}
