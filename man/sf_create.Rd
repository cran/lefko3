% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/popchar.R
\name{sf_create}
\alias{sf_create}
\title{Create Stageframe for Population Matrix Projection Analysis}
\usage{
sf_create(
  sizes,
  stagenames = NA,
  repstatus = 1,
  obsstatus = 1,
  propstatus = NA,
  immstatus = NA,
  matstatus = 1,
  minage = NA,
  maxage = NA,
  indataset = NA,
  binhalfwidth = 0.5,
  comments = NA,
  ipmbins = 100,
  roundsize = 5
)
}
\arguments{
\item{sizes}{A numeric vector of the typical or representative size of each
life history stage.}

\item{stagenames}{An optional vector of stage names, in the same order as
elements in sizes. If an IPM or function-based matrix with many stages is
desired, then two stages that occur within the dataset and represent the 
lower and upper size limits of the IPM must be marked as \code{ipm} in this 
vector. These stages must be mature stages, and should have all 
characteristics other than size equal. If two or more groups of stages, each
with its own characteristics, are to be developed for an IPM, then an even
number of stages with two stages marking the minimum and maximum size of
each group should be marked, with all other characteristics equal within
each group.}

\item{repstatus}{A vector denoting the binomial reproductive status of each
life history stage. Defaults to 1.}

\item{obsstatus}{A vector denoting the binomial observation status of each
life history stage. Defaults to 1, but may be changed for unobservable 
stages.}

\item{propstatus}{A vector denoting whether each life history stage is a 
propagule. Such stages are generally only used in fecundity estimation. 
Defaults to NA.}

\item{immstatus}{A vector denoting whether each stage is immature. Must be
composed of binomial values if given. Defaults to NA.}

\item{matstatus}{A vector denoting whether each stage is mature. Must be
composed of binomial values if given. Defaults to 1 for all stages defined 
in \code{sizes}.}

\item{minage}{An optional vector denoting the minimum age at which a stage
can occur. Only used in age x stage matrix development. Defaults to NA.}

\item{maxage}{An optional vector denoting the maximum age at which a stage
should occur. Only used in age x stage matrix development. Defaults to NA.}

\item{indataset}{A vector designating which stages are found within the 
dataset. While \code{\link{rlefko2}()} and \code{\link{rlefko3}()} can use
all stages in the input dataset, \code{\link{flefko3}()} and
\code{\link{flefko2}()} can only handle size-classified stages with
non-overlapping combinations of size and reproductive status, plus one
immature stage. Stages that do not actually exist within the dataset should
be marked as 0 in this vector.}

\item{binhalfwidth}{A numeric vector giving the half-width of size bins.
Required to classify individuals appropriately within size classes.
Defaults to 0.5 for all sizes.}

\item{comments}{An optional vector of text entries holding useful text
descriptions of all stages.}

\item{ipmbins}{If an IPM is desired, then this parameter sets the number of
stages to create for that IPM. This number is in addition to any stages
that are not size-classified. Defaults to 100, and numbers greater than this
yield a warning about the loss of statistical power and increasing chance of
matrix over-parameterization resulting from increasing numbers of stages.}

\item{roundsize}{This parameter sets the precision of size classification,
and equals the number of digits used in rounding sizes. Defaults to 5.}
}
\value{
A data frame of class \code{stageframe}, which includes information
on the stage name, size, reproductive status, observation status, propagule 
status, immaturity status, maturity status, presence within the core dataset, 
counts of similarly sized stages, raw bin half-width, and the minimum, 
center, and maximum of each size bin, as well as its width. If minimum and
maximum ages were specified, then these are also included. Also includes an 
empty string variable that can be used to describe stages meaningfully. This
object can be used as the \code{stageframe} input for \code{\link{flefko3}()} 
\code{\link{flefko2}()}, \code{\link{rlefko3}()}, and \code{\link{rlefko2}()}.

Variables in this data frame include the following:
\item{stage}{The unique names of the stages to be analyzed.}
\item{size}{The typical or representative size at which each stage occurs.}
\item{repstatus}{A binomial variable showing whether each stage is
reproductive.}
\item{obsstatus}{A binomial variable showing whether each stage is
observable.}
\item{propstatus}{A binomial variable showing whether each stage is a
propagule.}
\item{immstatus}{A binomial variable showing whether each stage can occur as
immature.}
\item{matstatus}{A binomial variable showing whether each stage occurs in
maturity.}
\item{indataset}{A binomial variable describing whether each stage occurs in
the input dataset.}
\item{binhalfwidth_raw}{The half-width of the size bin, as input.}
\item{min_age}{The minimum age at which the stage may occur.}
\item{max_age}{The maximum age at which the stage may occur.}
\item{sizebin_min}{The minimum size at which the stage may occur.}
\item{sizebin_max}{The maximum size at which the stage may occur.}
\item{sizebin_center}{The centroid of the size bin at which the stage may
occur.}
\item{sizebin_width}{The width of the size bin corresponding to the stage.}
\item{comments}{A text field for stage descriptions.}
}
\description{
\code{sf_create()} returns a data frame describing each ahistorical life
history stage in the life history model. This data frame can be used as 
input into MPM creation functions such as \code{\link{flefko3}()}, where it
determines how each stage is treated during matrix estimation.
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

ehrlen3mean <- lmean(ehrlen3)
ehrlen3mean$A[[1]]

# Cypripedium example
rm(list=ls(all=TRUE))

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

cypsupp2r <- supplemental(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
    "XSm", "Sm", "SD", "P1"),
  stage2 = c("SD", "SD", "P1", "P2", "P3", "SL", "SL", "SL", "SL", "rep",
    "rep"),
  eststage3 = c(NA, NA, NA, NA, NA, NA, "D", "XSm", "Sm", NA, NA),
  eststage2 = c(NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", NA, NA),
  givenrate = c(0.10, 0.20, 0.20, 0.20, 0.25, 0.40, NA, NA, NA, NA, NA),
  multiplier = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.5, 0.5),
  type =c(1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3),
  stageframe = cypframe_raw, historical = FALSE)

cypmatrix2r <- rlefko2(data = cypraw_v1, stageframe = cypframe_raw, 
  year = "all", patch = "all", stages = c("stage3", "stage2", "stage1"),
  size = c("size3added", "size2added"), supplement = cypsupp2r,
  yearcol = "year2", patchcol = "patchid", indivcol = "individ")
                       
cyp2mean <- lmean(cypmatrix2r)
cyp2mean

}
