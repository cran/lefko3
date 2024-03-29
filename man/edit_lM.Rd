% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{edit_lM}
\alias{edit_lM}
\title{Edit an MPM based on Supplemental Data}
\usage{
edit_lM(
  mpm,
  pop = NULL,
  patch = NULL,
  year2 = NULL,
  stage3 = NULL,
  stage2 = NULL,
  stage1 = NULL,
  age2 = NULL,
  eststage3 = NULL,
  eststage2 = NULL,
  eststage1 = NULL,
  estage2 = NULL,
  givenrate = NULL,
  multiplier = NULL,
  type = NULL,
  type_t12 = NULL
)
}
\arguments{
\item{mpm}{The \code{lefkoMat} object to be edited.}

\item{pop}{A string vector denoting the populations to be edited. Defaults
to \code{NULL}, in which case all populations are edited.}

\item{patch}{A string vector denoting the patches to be edited. Defaults
to \code{NULL}, in which case all patches are edited.}

\item{year2}{A string vector denoting the years to be edited. Defaults
to \code{NULL}, in which case all years are edited.}

\item{stage3}{The name of the stage in occasion \emph{t}+1 in the transition
to be replaced. Abbreviations for groups of stages are also usable (see
\code{Notes}). Required in all stage-based and age-by-stage MPMs.}

\item{stage2}{The name of the stage in occasion \emph{t} in the transition
to be replaced. Abbreviations for groups of stages are also usable (see
\code{Notes}). Required in all stage-based and age-by-stage MPMs.}

\item{stage1}{The name of the stage in occasion \emph{t}-1 in the transition
to be replaced. Only needed if a historical matrix is to be produced.
Abbreviations for groups of stages are also usable (see \code{Notes}).
Required for historical stage-based MPMs.}

\item{age2}{An integer vector of the ages in occasion \emph{t} to use in
transitions to be changed or replaced. Required for all age- and
age-by-stage MPMs.}

\item{eststage3}{The name of the stage to replace \code{stage3} in a proxy
transition. Only needed if a transition will be replaced by another
estimated transition, and only in stage-based and age-by-stage MPMs.}

\item{eststage2}{The name of the stage to replace \code{stage2} in a proxy
transition. Only needed if a transition will be replaced by another
estimated transition, and only in stage-based and age-by-stage MPMs.}

\item{eststage1}{The name of the stage to replace \code{stage1} in a proxy
historical transition. Only needed if a transition will be replaced by
another estimated transition, and the matrix to be estimated is historical
and stage-based. Stage \code{NotAlive} is also possible for raw hMPMs as a
means of handling the prior stage for individuals entering the population in
occasion \emph{t}.}

\item{estage2}{The age at time \emph{t} to replace \code{age2} in a proxy
transition. Only needed if a transition will be replaced by another
estimated transition, and only in age-based and age-by-stage MPMs.}

\item{givenrate}{A fixed rate or probability to replace for the transition
described by \code{stage3}, \code{stage2}, and \code{stage1}.}

\item{multiplier}{A vector of numeric multipliers for fecundity or for proxy
transitions. Defaults to \code{1}.}

\item{type}{A vector denoting the kind of transition between occasions
\emph{t} and \emph{t}+1 to be replaced. This should be entered as \code{1},
\code{S}, or \code{s} for the replacement of a survival transition;
\code{2}, \code{F}, or \code{f} for the replacement of a fecundity
transition; or \code{3}, \code{R}, or \code{r} for a fecundity multiplier.
If empty or not provided, then defaults to \code{1} for survival transition.}

\item{type_t12}{An optional vector denoting the kind of transition between
occasions \emph{t}-1 and \emph{t}. Only necessary if a historical MPM in
deVries format is desired. This should be entered as \code{1}, \code{S}, or
\code{s} for a survival transition; or \code{2}, \code{F}, or \code{f} for a
fecundity transitions. Defaults to \code{1} for survival transition, with
impacts only on the construction of deVries-format hMPMs.}
}
\value{
An edited copy of the original MPM is returned, also as a
\code{lefkoMat} object.
}
\description{
Function \code{edit_lM()} edits existing \code{lefkoMat} objects with
external data supplied by the user. The effects are similar to function
\code{\link{supplemental}()}, though function \code{edit_lM()} allows
individuals matrices within \code{lefkoMat} objects to be edited after
creation, while \code{\link{supplemental}()} provides external data that
modifies all matrices within a \code{lefkoMat} object.
}
\section{Notes}{

Entries in \code{stage3}, \code{stage2}, and \code{stage1} can include
abbreviations for groups of stages. Use \code{rep} if all reproductive
stages are to be used, \code{nrep} if all mature but non-reproductive stages
are to be used, \code{mat} if all mature stages are to be used, \code{immat}
if all immature stages are to be used, \code{prop} if all propagule stages
are to be used, \code{npr} if all non-propagule stages are to be used,
\code{obs} if all observable stages are to be used, \code{nobs} if all
unobservable stages are to be used, and leave empty or use \code{all} if all
stages in stageframe are to be used. Also use \code{groupX} to denote all
stages in group X (e.g. \code{group1} will use all stages in the respective
stageframe's group 1).
}

\examples{
data(cypdata)

cypraw_v1 <- verticalize3(data = cypdata, noyears = 6, firstyear = 2004,
  patchidcol = "patch", individcol = "plantid", blocksize = 4,
  sizeacol = "Inf2.04", sizebcol = "Inf.04", sizeccol = "Veg.04",
  repstracol = "Inf.04", repstrbcol = "Inf2.04", fecacol = "Pod.04",
  age_offset = 3, NAas0 = TRUE, NRasRep = TRUE)

cyp_rl <- rleslie(data = cypraw_v1, start_age = 0, last_age = 6, continue = TRUE,
  fecage_min = 3, year = "all", pop = NA, patch = "all", yearcol = "year2",
  patchcol = "patchid", indivcol = "individ")

ddd1 <- edit_lM(cyp_rl, age2 = c(0, 1, 2, 3, 4, 5, 6),
  givenrate = c(0.25, 0.25, 0.4, 0.4, NA, NA, NA),
  multiplier = c(NA, NA, NA, NA, 2000, 2000, 2000),
  type = c(1, 1, 1, 1, 3, 3, 3))
  
ddd1 <- edit_lM(ddd1, age2 = 6, multiplier = 1.5, type = 3, patch = "B",
  year2 = "2005")

}
\seealso{
\code{\link{supplemental}()}
}
