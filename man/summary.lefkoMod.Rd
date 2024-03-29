% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelselection.R
\name{summary.lefkoMod}
\alias{summary.lefkoMod}
\title{Summary of Class "lefkoMod"}
\usage{
\method{summary}{lefkoMod}(object, ...)
}
\arguments{
\item{object}{An R object of class \code{lefkoMod} resulting from
\code{\link{modelsearch}()}.}

\item{...}{Other parameters currently not utilized.}
}
\value{
A summary of the object, showing the best-fit models for all vital
rates, with constants of 0 or 1 used for unestimated models. This is followed
by a summary of the number of models tested per vital rate, and a table
showing the names of the parameters used to model vital rates and represent
tested factors. At the end is a section describing the numbers of individuals
and of individual transitions used to estimate each vital rate best-fit
model, along with the accuracy of each binomial model.
}
\description{
A function to summarize objects of class \code{lefkoMod}. This function shows
the best-fit models, summarizes the numbers of models in the model tables,
shows the criterion used to determine the best-fit models, and provides some
basic quality control information.
}
\examples{
\donttest{
# Lathyrus example
data(lathyrus)

sizevector <- c(0, 4.6, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8,
  9)
stagevector <- c("Sd", "Sdl", "Dorm", "Sz1nr", "Sz2nr", "Sz3nr", "Sz4nr", 
  "Sz5nr", "Sz6nr", "Sz7nr", "Sz8nr", "Sz9nr", "Sz1r", "Sz2r", "Sz3r", 
  "Sz4r", "Sz5r", "Sz6r", "Sz7r", "Sz8r", "Sz9r")
repvector <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
obsvector <- c(0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
matvector <- c(0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
immvector <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0)
indataset <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
binvec <- c(0, 4.6, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)

lathframeln <- sf_create(sizes = sizevector, stagenames = stagevector, 
  repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
  immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
  propstatus = propvector)

lathvertln <- verticalize3(lathyrus, noyears = 4, firstyear = 1988, 
  patchidcol = "SUBPLOT", individcol = "GENET", blocksize = 9,
  juvcol = "Seedling1988", sizeacol = "lnVol88", repstracol = "Intactseed88",
  fecacol = "Intactseed88", deadacol = "Dead1988", 
  nonobsacol = "Dormant1988", stageassign = lathframeln,
  stagesize = "sizea", censorcol = "Missing1988", censorkeep = NA,
  NAas0 = TRUE, censor = TRUE)

lathvertln$feca2 <- round(lathvertln$feca2)
lathvertln$feca1 <- round(lathvertln$feca1)
lathvertln$feca3 <- round(lathvertln$feca3)

lathmodelsln2 <- modelsearch(lathvertln, historical = FALSE, 
  approach = "mixed", suite = "main",
  vitalrates = c("surv", "obs", "size", "repst", "fec"), juvestimate = "Sdl",
  bestfit = "AICc&k", sizedist = "gaussian", fecdist = "poisson",
  indiv = "individ", patch = "patchid", year = "year2",
  year.as.random = TRUE, patch.as.random = TRUE, show.model.tables = TRUE,
  quiet = "partial")

summary(lathmodelsln2)
}

}
