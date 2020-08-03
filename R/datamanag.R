#' Create Historical Vertical Data Frames From Horizontal Data Frames
#'
#' \code{verticalize3()} returns a vertically formatted demographic data frame 
#' organized to create historical projection matrices, given a horizontally
#' formatted input data frame.
#'
#' @param data The horizontal data file.
#' @param noyears The number of years or observation periods in the dataset.
#' @param firstyear The first year or time of observation.
#' @param popidcol A variable name or column number corresponding to the identity
#' of the population for each individual.
#' @param patchidcol A variable name or column number corresponding to the
#' identity of the patch for each individual, if patches have been designated
#' within populations.
#' @param individcol A variable name or column number corresponding to the
#' identity of each individual.
#' @param blocksize The number of variables corresponding to each time step in 
#' the input dataset designated in \code{data}.
#' @param xcol A variable name or column number corresponding to the x 
#' coordinate of each individual in Cartesian space.
#' @param ycol A variable name or column number corresponding to the y 
#' coordinate of each individual in Cartesian space.
#' @param juvcol A variable name or column number that marks individuals in
#' immature stages within the dataset. The \code{verticalize3()} function assumes 
#' that immature individuals are identified in this variable marked with a 
#' number equal to or greater than 1, and that mature individuals are marked 
#' as 0 or NA.
#' @param size1col A variable name or column number corresponding to the
#' size entry associated with the first year or observation time in the dataset.
#' @param size2col A second variable name or column number corresponding to the
#' size entry associated with the first year or observation time in the dataset.
#' @param size3col A third variable name or column number corresponding to the
#' size entry associated with the first year or observation time in the dataset.
#' @param repstr1col A variable name or column number corresponding to the
#' production of reproductive structures, such as flowers, associated with the 
#' first year or observation period in the input dataset. This can be binomial or 
#' count data, and is used to in analysis of the probability of reproduction.
#' @param repstr2col A second variable name or column number corresponding to
#' the production of reproductive structures, such as flowers, associated with
#' the first year or observation period in the input dataset. This can be  
#' binomial or count data, and is used to in analysis of the probability of
#' reproduction.
#' @param fec1col A variable name or column number denoting fecundity associated
#' with the first year or observation time in the input dataset. This may
#' represent egg counts, fruit counts, seed production, etc.
#' @param fec2col A second variable name or column number denoting fecundity
#' associated with the first year or observation time in the input dataset. This 
#' may represent egg counts, fruit counts, seed production, etc.
#' @param alive1col A variable name or column number that provides information
#' an whether an individual is alive at a given time. If used, living status
#' must be designated as binomial (living = 1, dead = 0).
#' @param dead1col A variable name or column number that provides information
#' an whether an individual is alive at a given time. If used, dead status
#' must be designated as binomial (dead = 1, living = 0).
#' @param obs1col A variable name or column number providing information on
#' whether an individual is in an observable stage at a given time. If used,
#' observation status must be designated as binomial (observed = 1, not 
#' observed = 0).
#' @param nonobs1col A variable name or column number providing information on
#' whether an individual is in an unobservable stage at a given time. If used,
#' observation status must be designated as binomial (observed = 0, not 
#' observed = 1).
#' @param censorcol A variable name or column number corresponding to 
#' the first entry of a censor variable, used to distinguish between entries to 
#' use and entries not to use, or to designate entries with special issues that 
#' require further attention. If used, this should be associated with the first 
#' year or observation time, and all other years or times must also have censor 
#' columns.
#' @param repstrrel This is a scalar modifier for that makes the variable in
#' \code{repstr2col} equivalent to \code{repstr1col}. This can be useful if two 
#' reproductive status variables have related but unequal units, for example
#' if \code{repstr1col} refers to one-flowered stems while \code{repstr2col} refers to
#' two-flowered stems.
#' @param fecrel This is a scalar modifier for that makes the variable in
#' \code{fec2col} equivalent to \code{fec1col}. This can be useful if two fecundity
#' variables have related but unequal units.
#' @param stageassign The stageframe object identifying the life history model
#' being operationalized.
#' @param stagesize A variable name or column number describing which size
#' variable to use in stage estimation. Defaults to NA, and can also take
#' \code{sizea}, \code{sizeb}, \code{sizec}, or \code{sizeadded}, depending on which size variable 
#' is chosen.
#' @param censorkeep The value of the censoring variable identifying data
#' that should be included in analysis. Defaults to 1, but may take any value
#' including NA.
#' @param censor A logical variable determining whether the output data
#' should be censored using the variable defined in \code{censorcol}. Defaults
#' to FALSE.
#' @param spacing The spacing at which density should be estimated, if density
#' estimation is desired and x and y coordinates are supplied. Given in the
#' same units as those used in the x and y coordinates given in \code{xcol} and
#' \code{ycol}. Defaults to NA.
#' @param NAas0 If TRUE, then all NA entries for size and fecundity variables
#' will be set to 0. This can help increase the sample size analyzed by
#' \code{\link{modelsearch}()}, but should only be used when it is clear that this
#' substitution is biologically realistic. Defaults to FALSE.
#' @param NRasRep If TRUE, then will treat non-reproductive but mature
#' individuals as reproductive during stage assignment. This can be useful,
#' for example, when a matrix is desired without separation of reproductive
#' and non-reproductive but mature stages of the same size. Only used if
#' \code{stageassign} is set to a stageframe. Defaults to FALSE.
#' @param reduce A logical variable determining whether invariant state
#' variables should be removed from the output dataset. For example, if
#' all living individuals are always observable, then variables identifying
#' observation status are invariant and can be removed. Defaults to TRUE.
#' 
#' @return If all inputs are properly formatted, then this function will output
#' a historical vertical data frame, meaning that the output data frame will have
#' three consecutive years of size and reproductive data per individual per row.
#' This data frame is in standard format for all functions used in \code{lefko3},
#' and so can be used without further modification.
#' 
#' Variables in this data frame include the following:
#' \item{rowid}{Unique identifier for the row of the data frame.}
#' \item{popid}{Unique identifier for the population, if given.}
#' \item{patchid}{Unique identifier for patch within population, if
#' given.}
#' \item{individ}{Unique identifier for the individual.}
#' \item{year2}{Year or time step at time \emph{t}.}
#' \item{xpos,ypos}{X and Y position in Cartesian space, if given.}
#' \item{sizea1,sizea2,sizea3}{Main size measurement in times \emph{t}-1, 
#' \emph{t}, and \emph{t}+1, respectively.}
#' \item{sizeb1,sizeb2,sizeb3}{Secondary size measurement in times 
#' \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{sizec1,sizec2,sizec3}{Tertiary measurement in times \emph{t}-1, 
#' \emph{t}, and \emph{t}+1, respectively.}
#' \item{censor1,censor2,censor3}{Censor state values in times \emph{t}-1,
#' \emph{t}, and \emph{t}+1.}
#' \item{repstra1,repstrb1,repstrc1}{Main, secondary, and tertiary numbers of 
#' reproductive structures in time \emph{t}-1.}
#' \item{repstra2,repstrb2,repstrc2}{Main, secondary, and tertiary numbers of 
#' reproductive structures in time \emph{t}.}
#' \item{repstra3,repstrb3,repstrc3}{Main, secondary, and tertiary numbers of 
#' reproductive structures in time \emph{t}+1.}
#' \item{feca1,fecb1}{Main and secondary numbers of offspring in time 
#' \emph{t}-1.}
#' \item{feca2,fecb2}{Main and secondary numbers of offspring in time 
#' \emph{t}.}
#' \item{feca3,fecb3}{Main and secondary numbers of offspring in time 
#' \emph{t}+1.}
#' \item{size1added,size2added,size3added}{Sum of primary, secondary, and
#' tertiary size measurements in timea \emph{t}-1, \emph{t}, and \emph{t}+1, 
#' respectively.}
#' \item{repstr1added,repstr2added,repstr3added}{Sum of primary, secondary,
#' and tertiary reproductive structures in times \emph{t}-1, \emph{t}, and 
#' \emph{t}+1, respectively.}
#' \item{fec1added,fec2added,fec3added}{Sum of primary and secondary 
#' fecundity in times \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{obsstatus1,obsstatus2,obsstatus3}{Binomial observation state in times
#' \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{repstatus1,repstatus2,repstatus3}{Binomial reproductive state in 
#' times \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{fecstatus1,fecstatus2,fecstatus3}{Binomial offspring production
#' state in times \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{firstseen}{Year or time step of first observation.}
#' \item{lastseen}{Year or time step of last observation.}
#' \item{xcorr,ycorr}{Overall x and y coordinates of individual in
#' Cartesian space.}
#' \item{alive1,alive2,alive3}{Binomial state as alive in times \emph{t}-1, 
#' \emph{t}, and \emph{t}+1, respectively.}
#' \item{obsage}{Observed age in time \emph{t}, assuming first observation
#' corresponds to age = 0.}
#' \item{obslifespan}{Observed lifespan, given as `lastseen - firstseen + 1`.}
#' \item{matstatus1,matstatus2,matstatus3}{Binomial state as mature.}
#' \item{density}{Density of individuals per unit designated in `spacing`.
#' Only given if spacing is not NA.}
#' 
#' @examples
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
#' summary(lathvert)
#' 
#' @export
verticalize3 <- function(data, noyears, firstyear, popidcol = 0, patchidcol = 0, individcol= 0, 
                         blocksize, xcol = 0, ycol = 0, juvcol = 0, size1col, size2col = 0, 
                         size3col = 0, repstr1col = 0, repstr2col = 0, fec1col = 0, fec2col = 0, 
                         alive1col = 0, dead1col = 0, obs1col = 0, nonobs1col = 0, censorcol = 0, 
                         repstrrel = 1, fecrel = 1, stageassign = NA, stagesize = NA, 
                         censorkeep = 1, censor = FALSE, spacing = NA, NAas0 = FALSE, 
                         NRasRep = FALSE, reduce = TRUE) {
  
  rowid <- alive2 <- indataset <- censor1 <- censor2 <- censor3 <- NULL
  
  popid <- NA
  patchid <- NA
  individ <- NA
  
  #This first section tests the input for valid entries
  if (is.character(popidcol)) {
    if (is.element(popidcol, names(data))) {
      true.popidcol <- which(names(data) == popidcol)
      popidcol <- true.popidcol
    } else {stop("Please enter popidcol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(patchidcol)) {
    if (is.element(patchidcol, names(data))) {
      true.patchidcol <- which(names(data) == patchidcol)
      patchidcol <- true.patchidcol
    } else {stop("Please enter patchidcol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(individcol)) {
    if (is.element(individcol, names(data))) {
      true.individcol <- which(names(data) == individcol)
      individcol <- true.individcol
    } else {stop("Please enter individcol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(xcol)) {
    if (is.element(xcol, names(data))) {
      true.xcol <- which(names(data) == xcol)
      xcol <- true.xcol
    } else {stop("Please enter xcol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(ycol)) {
    if (is.element(ycol, names(data))) {
      true.ycol <- which(names(data) == ycol)
      ycol <- true.ycol
    } else {stop("Please enter ycol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(juvcol)) {
    if (is.element(juvcol, names(data))) {
      true.juvcol <- which(names(data) == juvcol)
      juvcol <- true.juvcol
    } else {stop("Please enter juvcol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(size1col)) {
    if (is.element(size1col, names(data))) {
      true.size1col <- which(names(data) == size1col)
      size1col <- true.size1col
    } else {stop("Please enter size1col exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(size2col)) {
    if (is.element(size2col, names(data))) {
      true.size2col <- which(names(data) == size2col)
      size2col <- true.size2col
    } else {stop("Please enter size2col exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(size3col)) {
    if (is.element(size3col, names(data))) {
      true.size3col <- which(names(data) == size3col)
      size3col <- true.size3col
    } else {stop("Please enter size3col exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(repstr1col)) {
    if (is.element(repstr1col, names(data))) {
      true.repstr1col <- which(names(data) == repstr1col)
      repstr1col <- true.repstr1col
    } else {stop("Please enter repstr1col exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(repstr2col)) {
    if (is.element(repstr2col, names(data))) {
      true.repstr2col <- which(names(data) == repstr2col)
      repstr2col <- true.repstr2col
    } else {stop("Please enter repstr2col exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(fec1col)) {
    if (is.element(fec1col, names(data))) {
      true.fec1col <- which(names(data) == fec1col)
      fec1col <- true.fec1col
    } else {stop("Please enter fec1col exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(fec2col)) {
    if (is.element(fec2col, names(data))) {
      true.fec2col <- which(names(data) == fec2col)
      fec2col <- true.fec2col
    } else {stop("Please enter fec2col exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(censorcol)) {
    if (is.element(censorcol, names(data))) {
      true.censorcol <- which(names(data) == censorcol)
      censorcol <- true.censorcol
    } else {stop("Please enter censorcol exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(alive1col)) {
    if (is.element(alive1col, names(data))) {
      true.alive1col <- which(names(data) == alive1col)
      alive1col <- true.alive1col
    } else {stop("Please enter alive1col exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(dead1col)) {
    if (is.element(dead1col, names(data))) {
      true.dead1col <- which(names(data) == dead1col)
      dead1col <- true.dead1col
    } else {stop("Please enter dead1col exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(obs1col)) {
    if (is.element(obs1col, names(data))) {
      true.obs1col <- which(names(data) == obs1col)
      obs1col <- true.obs1col
    } else {stop("Please enter obs1col exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if (is.character(nonobs1col)) {
    if (is.element(nonobs1col, names(data))) {
      true.nonobs1col <- which(names(data) == nonobs1col)
      nonobs1col <- true.nonobs1col
    } else {stop("Please enter nonobs1col exactly as it appears in the dataset.", call. = FALSE)}
  }
  
  if(!all(is.na(stageassign))) {
    if(!is.element("stageframe", class(stageassign))) {
      stop("The stageassign option can only take NA or a stageframe object as input.", call. = FALSE)
    }
    if(length(intersect(stagesize, c("sizea", "sizeb", "sizec", "sizeadded"))) == 0) {
      stop("The stagesize option must equal NA, 'sizea', 'sizeb', 'sizec', or 'sizeadded'. No other values are permitted.", call. = FALSE)
    }
  }
  
  if (!is.na(spacing)) {
    if (xcol == 0 | ycol == 0) {
      stop("Density estimation cannot proceed without valid x and y coordinates.")
    }
    
    if (is.character(spacing)) {
      stop("The spacing option requires either a number, or defaults to NA.")
    }
  }
  
  #Now we come to the core part of the function that takes apart the dataset and makes it vertical
  if (popidcol > 0) {popid <- data[,popidcol];}
  if (patchidcol > 0) {patchid <- data[,patchidcol];}
  if (individcol > 0) {individ <- data[,individcol];}
  
  if (censor) {
    if (!is.element(censorkeep, data[,censorcol])) {
      fullcenvec <- as.vector(apply(as.matrix(c(1:noyears)), 1, function(X) {c(data[,(censorcol + (X - 1) * blocksize)])}))
      
      if (!is.element(censorkeep, fullcenvec)) {
        stop("Please enter a valid value for censorkeep. This value should occur in the censor variable within the dataset.", call. = FALSE)
      }
    }
  }
  
  noindivs <- dim(data)[1]
  
  popdatalist <- lapply(c(1:(noyears-1)), function(j) {
    currentyear <- firstyear + (j-1)
    
    if (j == 1)
    {
      yearblock <- cbind.data.frame(c(1:noindivs), popid, patchid, individ, currentyear)
      
      if (xcol > 0) {yearblock <- cbind.data.frame(yearblock, NA, NA, data[,xcol], data[,ycol], data[,xcol+blocksize], data[,ycol+blocksize])}
      if (xcol == 0) {yearblock <- cbind.data.frame(yearblock, NA, NA, NA, NA, NA, NA)}
      
      yearblock <- cbind.data.frame(yearblock, NA, data[,size1col], data[,(size1col+blocksize)], NA, data[,repstr1col], data[,(repstr1col+blocksize)])
      
      if (fec1col > 0) {
        yearblock <- cbind.data.frame(yearblock, NA, data[,fec1col], data[,(fec1col+blocksize)])
      } else {yearblock <- cbind.data.frame(yearblock, NA, NA, NA)}
      
      if (size2col > 0) {
        yearblock <- cbind.data.frame(yearblock, NA, data[,size2col], data[,(size2col+blocksize)])
      } else {yearblock <- cbind.data.frame(yearblock, NA, NA, NA)}
      
      if (repstr2col > 0) {
        yearblock <- cbind.data.frame(yearblock, NA, data[,repstr2col], data[,(repstr2col+blocksize)])
      } else {yearblock <- cbind.data.frame(yearblock, NA, NA, NA)}
      
      if (fec2col > 0) {
        yearblock <- cbind.data.frame(yearblock, NA, data[,fec2col], data[,(fec2col+blocksize)])
      } else {yearblock <- cbind.data.frame(yearblock, NA, NA, NA)}
      
      if (size3col > 0) {
        yearblock <- cbind.data.frame(yearblock, NA, data[,size3col], data[,(size3col+blocksize)])
      } else {yearblock <- cbind.data.frame(yearblock, NA, NA, NA)}
      
      if (censorcol > 0) {
        yearblock <- cbind.data.frame(yearblock, NA, data[,censorcol], data[,(censorcol+blocksize)])
      } else {yearblock <- cbind.data.frame(yearblock, NA, 1, 1)}
      
      if (alive1col > 0) {
        yearblock <- cbind.data.frame(yearblock, NA, data[,alive1col], data[,(alive1col+blocksize)])
      } else {yearblock <- cbind.data.frame(yearblock, NA, NA, NA)}
      
      if (dead1col > 0) {
        yearblock <- cbind.data.frame(yearblock, NA, data[,dead1col], data[,(dead1col+blocksize)])
      } else {yearblock <- cbind.data.frame(yearblock, NA, NA, NA)}
      
      if (obs1col > 0) {
        yearblock <- cbind.data.frame(yearblock, NA, data[,obs1col], data[,(obs1col+blocksize)])
      } else {yearblock <- cbind.data.frame(yearblock, NA, NA, NA)}
      
      if (nonobs1col > 0) {
        yearblock <- cbind.data.frame(yearblock, NA, data[,nonobs1col], data[,(nonobs1col+blocksize)])
      } else {yearblock <- cbind.data.frame(yearblock, NA, NA, NA)}
      
      if (juvcol > 0) {
        yearblock <- cbind.data.frame(yearblock, NA, data[,juvcol], data[,(juvcol+blocksize)])
      } else {yearblock <- cbind.data.frame(yearblock, NA, NA, NA)}
      
      names(yearblock) <- c("rowid", "popid", "patchid", "individ", "year2", "xpos1", "ypos1", "xpos2", 
                            "ypos2", "xpos3", "ypos3", "sizea1", "sizea2", "sizea3", "repstra1", 
                            "repstra2", "repstra3", "feca1", "feca2", "feca3", "sizeb1", "sizeb2", 
                            "sizeb3", "repstrb1", "repstrb2", "repstrb3", "fecb1", "fecb2", "fecb3", 
                            "sizec1", "sizec2", "sizec3", "censor1", "censor2", "censor3", "alivegiven1", 
                            "alivegiven2", "alivegiven3", "deadgiven1", "deadgiven2", "deadgiven3",
                            "obsgiven1", "obsgiven2", "obsgiven3", "nonobsgiven1", "nonobsgiven2", 
                            "nonobsgiven3", "juvgiven1", "juvgiven2", "juvgiven3")
    }
    
    if (j > 1)
    {
      yearblock <- cbind.data.frame(c(1:noindivs), popid, patchid, individ, currentyear)
      
      if (xcol > 0) {
        yearblock <- cbind.data.frame(yearblock, data[,xcol+(blocksize*(j-2))], data[,ycol+(blocksize*(j-2))], data[,xcol+(blocksize*(j-1))], data[,ycol+(blocksize*(j-1))], data[,xcol+blocksize*j], data[,ycol+blocksize*j])
      } else {yearblock <- cbind.data.frame(yearblock, NA, NA, NA, NA, NA, NA)}
      
      yearblock <- cbind.data.frame(yearblock, data[,(size1col+(blocksize*(j-2)))], data[,(size1col+(blocksize*(j-1)))], data[,(size1col+(blocksize*j))], data[,(repstr1col+(blocksize*(j-2)))], data[,(repstr1col+(blocksize*(j-1)))], data[,(repstr1col+(blocksize*j))])
      
      if (fec1col > 0) {
        yearblock <- cbind.data.frame(yearblock, data[,(fec1col+(blocksize*(j-2)))], data[,(fec1col+(blocksize*(j-1)))], data[,(fec1col+(blocksize*j))])
      } else {yearblock <- cbind.data.frame(yearblock, NA, NA, NA)}
      
      if (size2col > 0) {
        yearblock <- cbind.data.frame(yearblock, data[,(size2col+(blocksize*(j-2)))], data[,(size2col+(blocksize*(j-1)))], data[,(size2col+(blocksize*j))])
      } else {yearblock <- cbind.data.frame(yearblock, NA, NA, NA)}
      
      if (repstr2col > 0) {
        yearblock <- cbind.data.frame(yearblock, data[,(repstr2col+(blocksize*(j-2)))], data[,(repstr2col+(blocksize*(j-1)))], data[,(repstr2col+(blocksize*j))])
      } else {yearblock <- cbind.data.frame(yearblock, NA, NA, NA)}
      
      if (fec2col > 0) {
        yearblock <- cbind.data.frame(yearblock, data[,(fec2col+(blocksize*(j-2)))], data[,(fec2col+(blocksize*(j-1)))], data[,(fec2col+(blocksize*j))])
      } else {yearblock <- cbind.data.frame(yearblock, NA, NA, NA)}
      
      if (size3col > 0) {
        yearblock <- cbind.data.frame(yearblock, data[,(size3col+(blocksize*(j-2)))], data[,(size3col+(blocksize*(j-1)))], data[,(size3col+(blocksize*j))])
      } else {yearblock <- cbind.data.frame(yearblock, NA, NA, NA)}
      
      if (censorcol > 0) {
        yearblock <- cbind.data.frame(yearblock, data[,(censorcol+(blocksize*(j-2)))], data[,(censorcol+(blocksize*(j-1)))], data[,(censorcol+(blocksize*j))])
      } else {yearblock <- cbind.data.frame(yearblock, 1, 1, 1)}
      
      if (alive1col > 0) {
        yearblock <- cbind.data.frame(yearblock, data[,(alive1col+(blocksize*(j-2)))], data[,(alive1col+(blocksize*(j-1)))], data[,(alive1col+(blocksize*j))])
      } else {yearblock <- cbind.data.frame(yearblock, NA, NA, NA)}
      
      if (dead1col > 0) {
        yearblock <- cbind.data.frame(yearblock, data[,(dead1col+(blocksize*(j-2)))], data[,(dead1col+(blocksize*(j-1)))], data[,(dead1col+(blocksize*j))])
      } else {yearblock <- cbind.data.frame(yearblock, NA, NA, NA)}
      
      if (obs1col > 0) {
        yearblock <- cbind.data.frame(yearblock, data[,(obs1col+(blocksize*(j-2)))], data[,(obs1col+(blocksize*(j-1)))], data[,(obs1col+(blocksize*j))])
      } else {yearblock <- cbind.data.frame(yearblock, NA, NA, NA)}
      
      if (nonobs1col > 0) {
        yearblock <- cbind.data.frame(yearblock, data[,(nonobs1col+(blocksize*(j-2)))], data[,(nonobs1col+(blocksize*(j-1)))], data[,(nonobs1col+(blocksize*j))])
      } else {yearblock <- cbind.data.frame(yearblock, NA, NA, NA)}
      
      if (juvcol > 0) {
        yearblock <- cbind.data.frame(yearblock, data[,(juvcol+(blocksize*(j-2)))], data[,(juvcol+(blocksize*(j-1)))], data[,(juvcol+(blocksize*j))])
      } else {yearblock <- cbind.data.frame(yearblock, NA, NA, NA)}
      
      
      names(yearblock) <- c("rowid", "popid", "patchid", "individ", "year2", "xpos1", "ypos1", "xpos2", 
                            "ypos2", "xpos3", "ypos3", "sizea1", "sizea2", "sizea3", "repstra1", 
                            "repstra2", "repstra3", "feca1", "feca2", "feca3", "sizeb1", "sizeb2", 
                            "sizeb3", "repstrb1", "repstrb2", "repstrb3", "fecb1", "fecb2", "fecb3", 
                            "sizec1", "sizec2", "sizec3", "censor1", "censor2", "censor3", "alivegiven1", 
                            "alivegiven2", "alivegiven3", "deadgiven1", "deadgiven2", "deadgiven3",
                            "obsgiven1", "obsgiven2", "obsgiven3", "nonobsgiven1", "nonobsgiven2", 
                            "nonobsgiven3", "juvgiven1", "juvgiven2", "juvgiven3")
    }
    return(yearblock)
  })
  
  popdata <- do.call("rbind.data.frame", popdatalist)
  
  popdata$juvgiven1[which(is.na(popdata$juvgiven1))] <- 0
  popdata$juvgiven2[which(is.na(popdata$juvgiven2))] <- 0
  popdata$juvgiven3[which(is.na(popdata$juvgiven3))] <- 0
  popdata$juvgiven1[which(popdata$juvgiven1 > 0)] <- 1
  popdata$juvgiven2[which(popdata$juvgiven2 > 0)] <- 1
  popdata$juvgiven3[which(popdata$juvgiven3 > 0)] <- 1
  
  popdata$addedsize1 <- apply(as.matrix(c(1:dim(popdata)[1])), 1, function(X) {sum(popdata$sizea1[X], popdata$sizeb1[X], popdata$sizec1[X], na.rm = TRUE)})
  popdata$addedsize2 <- apply(as.matrix(c(1:dim(popdata)[1])), 1, function(X) {sum(popdata$sizea2[X], popdata$sizeb2[X], popdata$sizec2[X], na.rm = TRUE)})
  popdata$addedsize3 <- apply(as.matrix(c(1:dim(popdata)[1])), 1, function(X) {sum(popdata$sizea3[X], popdata$sizeb3[X], popdata$sizec3[X], na.rm = TRUE)})
  
  popdata$addedflower1 <- apply(as.matrix(c(1:dim(popdata)[1])), 1, function(X) {sum(popdata$repstra1[X], popdata$repstrb1[X] * repstrrel, na.rm = TRUE)})
  popdata$addedflower2 <- apply(as.matrix(c(1:dim(popdata)[1])), 1, function(X) {sum(popdata$repstra2[X], popdata$repstrb2[X] * repstrrel, na.rm = TRUE)})
  popdata$addedflower3 <- apply(as.matrix(c(1:dim(popdata)[1])), 1, function(X) {sum(popdata$repstra3[X], popdata$repstrb3[X] * repstrrel, na.rm = TRUE)})
  
  popdata$addedfruit1 <- apply(as.matrix(c(1:dim(popdata)[1])), 1, function(X) {sum(popdata$feca1[X], popdata$fecb1[X] * fecrel, na.rm = TRUE)})
  popdata$addedfruit2 <- apply(as.matrix(c(1:dim(popdata)[1])), 1, function(X) {sum(popdata$feca2[X], popdata$fecb2[X] * fecrel, na.rm = TRUE)})
  popdata$addedfruit3 <- apply(as.matrix(c(1:dim(popdata)[1])), 1, function(X) {sum(popdata$feca3[X], popdata$fecb3[X] * fecrel, na.rm = TRUE)})
  
  popdata$spryn1 <- 0
  popdata$spryn2 <- 0
  popdata$spryn3 <- 0
  
  popdata$spryn1[(which(popdata$sizea1 > 0))] <- 1
  popdata$spryn2[(which(popdata$sizea2 > 0))] <- 1
  popdata$spryn3[(which(popdata$sizea3 > 0))] <- 1
  popdata$spryn1[(which(popdata$sizeb1 > 0))] <- 1
  popdata$spryn2[(which(popdata$sizeb2 > 0))] <- 1
  popdata$spryn3[(which(popdata$sizeb3 > 0))] <- 1
  popdata$spryn1[(which(popdata$sizec1 > 0))] <- 1
  popdata$spryn2[(which(popdata$sizec2 > 0))] <- 1
  popdata$spryn3[(which(popdata$sizec3 > 0))] <- 1
  
  popdata$spryn1[(which(popdata$obsgiven1 == 1))] <- 1
  popdata$spryn2[(which(popdata$obsgiven2 == 1))] <- 1
  popdata$spryn3[(which(popdata$obsgiven3 == 1))] <- 1
  popdata$spryn1[(which(popdata$nonobsgiven1 == 1))] <- 0
  popdata$spryn2[(which(popdata$nonobsgiven2 == 1))] <- 0
  popdata$spryn3[(which(popdata$nonobsgiven3 == 1))] <- 0
  
  popdata$flryn1 <- 0
  popdata$flryn2 <- 0
  popdata$flryn3 <- 0
  
  popdata$flryn1[(which(popdata$repstra1 > 0))] <- 1
  popdata$flryn2[(which(popdata$repstra2 > 0))] <- 1
  popdata$flryn3[(which(popdata$repstra3 > 0))] <- 1
  popdata$flryn1[(which(popdata$repstrb1 > 0))] <- 1
  popdata$flryn2[(which(popdata$repstrb2 > 0))] <- 1
  popdata$flryn3[(which(popdata$repstrb3 > 0))] <- 1
  
  popdata$fecyn1 <- 0
  popdata$fecyn2 <- 0
  popdata$fecyn3 <- 0
  
  popdata$fecyn1[(which(popdata$feca1 > 0))] <- 1
  popdata$fecyn2[(which(popdata$feca2 > 0))] <- 1
  popdata$fecyn3[(which(popdata$feca3 > 0))] <- 1
  popdata$fecyn1[(which(popdata$fecb1 > 0))] <- 1
  popdata$fecyn2[(which(popdata$fecb2 > 0))] <- 1
  popdata$fecyn3[(which(popdata$fecb3 > 0))] <- 1
  
  #The dataset called popdata is already a vertical dataset, but we need to add some 
  #more variables to it, including age and corrected spatial bearings
  indivdata <- lapply(c(1:noindivs), function(i) {
    rowdata <- subset(popdata, rowid == i)
    
    yearsused <- c(rowdata$year2, (firstyear + noyears - 1))
    
    sproutpattern <- c(rowdata$spryn2, rowdata$spryn3[(which(rowdata$year2 == (firstyear + noyears - 2)))]) #This line can be used to add a CMR resighting history variable to the output
    
    sproutyears <- yearsused[which(sproutpattern == 1)]
    unsproutyears <- yearsused[which(sproutpattern == 0)]
    
    if (length(sproutyears) == 0) {warning(paste("Individual identified with no observation at rowID", i))}
    
    firstseen <- min(sproutyears)
    lastseen <- max(sproutyears)
    
    xcorr <- rowdata$xpos2
    ycorr <- rowdata$ypos2
    
    sproutxmean <- median(c(rowdata$xpos1[which(sproutpattern == 1)], rowdata$xpos2[which(sproutpattern == 1)], 
                            rowdata$xpos3[which(sproutpattern == 1)]), na.rm = TRUE)
    unsproutxmean <- median(c(rowdata$xpos1[which(sproutpattern == 0)], rowdata$xpos2[which(sproutpattern == 0)], 
                              rowdata$xpos3[which(sproutpattern == 0)]), na.rm = TRUE)
    sproutymean <- median(c(rowdata$ypos1[which(sproutpattern == 1)], rowdata$ypos2[which(sproutpattern == 1)], 
                            rowdata$ypos3[which(sproutpattern == 1)]), na.rm = TRUE)
    unsproutymean <- median(c(rowdata$ypos1[which(sproutpattern == 0)], rowdata$ypos2[which(sproutpattern == 0)], 
                              rowdata$ypos3[which(sproutpattern == 0)]), na.rm = TRUE)
    
    if (xcol > 0) {
      if(is.na(sproutxmean)) {sproutxmean <- 0}
      if(is.na(unsproutxmean)) {unsproutxmean <- 0}
      if (sproutxmean != unsproutxmean) {xcorr <- rep(sproutxmean, length(xcorr))}
    }
    
    if (ycol > 0) {
      if(is.na(sproutymean)) {sproutymean <- 0}
      if(is.na(unsproutymean)) {unsproutymean <- 0}
      if (sproutymean != unsproutymean) {ycorr <- rep(sproutymean, length(ycorr))}
    }
    
    finalrowdata <- cbind.data.frame(rowdata, firstseen, lastseen, xcorr, ycorr)
    names(finalrowdata) <- c("rowid", "popid", "patchid", "individ", "year2", "xpos1", "ypos1", 
                             "xpos2", "ypos2", "xpos3", "ypos3", "sizea1", "sizea2", "sizea3", 
                             "repstra1", "repstra2", "repstra3", "feca1", "feca2", "feca3", 
                             "sizeb1", "sizeb2", "sizeb3", "repstrb1", "repstrb2", "repstrb3", 
                             "fecb1", "fecb2", "fecb3", "sizec1", "sizec2", "sizec3", "censor1", 
                             "censor2", "censor3", "alivegiven1", "alivegiven2", "alivegiven3", 
                             "deadgiven1", "deadgiven2", "deadgiven3", "obsgiven1", "obsgiven2", 
                             "obsgiven3", "nonobsgiven1", "nonobsgiven2", "nonobsgiven3", 
                             "juvgiven1", "juvgiven2", "juvgiven3", "size1added", "size2added", 
                             "size3added", "repstr1added", "repstr2added", "repstr3added", 
                             "fec1added", "fec2added", "fec3added", "obsstatus1", "obsstatus2", 
                             "obsstatus3", "repstatus1", "repstatus2", "repstatus3", "fecstatus1", 
                             "fecstatus2", "fecstatus3", "firstseen", "lastseen", "xcorr", "ycorr")
    return(finalrowdata)
  })
  popdatanew <- do.call("rbind.data.frame", indivdata)
  rownames(popdatanew) <- c(1:dim(popdatanew)[1])
  
  popdatanew$alive1 <- 0
  popdatanew$alive2 <- 0
  popdatanew$alive3 <- 0
  
  popdatanew$alive1[(which((((popdatanew$year2 - 1) - popdatanew$firstseen) >= 0) & ((popdatanew$lastseen - (popdatanew$year2 - 1)) >= 0)))] <- 1
  popdatanew$alive2[(which(((popdatanew$year2 - popdatanew$firstseen) >= 0) & ((popdatanew$lastseen - popdatanew$year2) >= 0)))] <- 1
  popdatanew$alive3[(which((((popdatanew$year2 + 1) - popdatanew$firstseen) >= 0) & ((popdatanew$lastseen - (popdatanew$year2 + 1)) >= 0)))] <- 1
  
  popdatanew$alive1[which(popdatanew$alivegiven1 == 1)] <- 1
  popdatanew$alive2[which(popdatanew$alivegiven2 == 1)] <- 1
  popdatanew$alive3[which(popdatanew$alivegiven3 == 1)] <- 1
  popdatanew$alive1[which(popdatanew$deadgiven1 == 1)] <- 0
  popdatanew$alive2[which(popdatanew$deadgiven2 == 1)] <- 0
  popdatanew$alive3[which(popdatanew$deadgiven3 == 1)] <- 0
  
  popdatanew$obsage <- popdatanew$year2 - popdatanew$firstseen
  popdatanew$obslifespan <- popdatanew$lastseen - popdatanew$firstseen
  
  popdatareal <- subset(popdatanew, subset = (alive2 == 1))
  
  popdatareal$sizea1c <- popdatareal$sizea1
  popdatareal$sizea2c <- popdatareal$sizea2
  popdatareal$sizea3c <- popdatareal$sizea3
  popdatareal$sizeb1c <- popdatareal$sizeb1
  popdatareal$sizeb2c <- popdatareal$sizeb2
  popdatareal$sizeb3c <- popdatareal$sizeb3
  popdatareal$sizec1c <- popdatareal$sizec1
  popdatareal$sizec2c <- popdatareal$sizec2
  popdatareal$sizec3c <- popdatareal$sizec3
  
  popdatareal$repstra1c <- popdatareal$repstra1
  popdatareal$repstra2c <- popdatareal$repstra2
  popdatareal$repstra3c <- popdatareal$repstra3
  popdatareal$repstrb1c <- popdatareal$repstrb1
  popdatareal$repstrb2c <- popdatareal$repstrb2
  popdatareal$repstrb3c <- popdatareal$repstrb3
  
  popdatareal$feca1c <- popdatareal$feca1
  popdatareal$feca2c <- popdatareal$feca2
  popdatareal$feca3c <- popdatareal$feca3
  popdatareal$fecb1c <- popdatareal$fecb1
  popdatareal$fecb2c <- popdatareal$fecb2
  popdatareal$fecb3c <- popdatareal$fecb3
  
  popdatareal$sizea1c[which(is.na(popdatareal$sizea1))] <- 0
  popdatareal$sizea2c[which(is.na(popdatareal$sizea2))] <- 0
  popdatareal$sizea3c[which(is.na(popdatareal$sizea3))] <- 0
  
  popdatareal$repstra1c[which(is.na(popdatareal$repstra1))] <- 0
  popdatareal$repstra2c[which(is.na(popdatareal$repstra2))] <- 0
  popdatareal$repstra3c[which(is.na(popdatareal$repstra3))] <- 0
  
  popdatareal$sizeb1c[which(is.na(popdatareal$sizeb1))] <- 0
  popdatareal$sizeb2c[which(is.na(popdatareal$sizeb2))] <- 0
  popdatareal$sizeb3c[which(is.na(popdatareal$sizeb3))] <- 0
  
  popdatareal$repstrb1c[which(is.na(popdatareal$repstrb1))] <- 0
  popdatareal$repstrb2c[which(is.na(popdatareal$repstrb2))] <- 0
  popdatareal$repstrb3c[which(is.na(popdatareal$repstrb3))] <- 0
  
  popdatareal$sizec1c[which(is.na(popdatareal$sizec1))] <- 0
  popdatareal$sizec2c[which(is.na(popdatareal$sizec2))] <- 0
  popdatareal$sizec3c[which(is.na(popdatareal$sizec3))] <- 0
  
  popdatareal$feca1c[which(is.na(popdatareal$feca1))] <- 0
  popdatareal$feca2c[which(is.na(popdatareal$feca2))] <- 0
  popdatareal$feca3c[which(is.na(popdatareal$feca3))] <- 0
  
  popdatareal$fecb1c[which(is.na(popdatareal$fecb1))] <- 0
  popdatareal$fecb2c[which(is.na(popdatareal$fecb2))] <- 0
  popdatareal$fecb3c[which(is.na(popdatareal$fecb3))] <- 0
  
  if(!all(is.na(stageassign))) {
    if (stagesize == "sizeadded") {
      stagesizecol1 <- which(names(popdatareal) == "size1added")
      stagesizecol2 <- which(names(popdatareal) == "size2added")
      stagesizecol3 <- which(names(popdatareal) == "size3added")
    } else if (stagesize == "sizec") {
      stagesizecol1 <- which(names(popdatareal) == "sizec1c")
      stagesizecol2 <- which(names(popdatareal) == "sizec2c")
      stagesizecol3 <- which(names(popdatareal) == "sizec3c")
    } else if (stagesize == "sizeb") {
      stagesizecol1 <- which(names(popdatareal) == "sizeb1c")
      stagesizecol2 <- which(names(popdatareal) == "sizeb2c")
      stagesizecol3 <- which(names(popdatareal) == "sizeb3c")
    } else {
      stagesizecol1 <- which(names(popdatareal) == "sizea1c")
      stagesizecol2 <- which(names(popdatareal) == "sizea2c")
      stagesizecol3 <- which(names(popdatareal) == "sizea3c")
    }
    
    ltdframe <- subset(stageassign, indataset == 1)
    ltdframe$stagenames <- as.character(ltdframe$stagenames)
    
    popdatareal$stage1 <- apply(as.matrix(c(1:dim(popdatareal)[1])), 1, function(X) {
      if (popdatareal$alive1[X] == 1) {
        if (is.na(popdatareal[X, stagesizecol1])) {
          popdatareal[X, stagesizecol1] <- 0
        }
        mainstages <- intersect(which(ltdframe$sizebin_min < popdatareal[X, stagesizecol1]), 
                                which(ltdframe$sizebin_max >= popdatareal[X, stagesizecol1]))
        jmstages <- which(ltdframe$immstatus == popdatareal$juvgiven1[X])
        obsstages <- which(ltdframe$obsstatus == popdatareal$obsstatus1[X])
        repstages <- which(ltdframe$repstatus == popdatareal$repstatus1[X])
        
        if (!NRasRep) {
          choicestage <- intersect(intersect(mainstages, jmstages), intersect(obsstages, repstages))
        } else {
          choicestage <- intersect(intersect(mainstages, jmstages), obsstages)
        }
        
        if (all(is.na(choicestage))) {
          stop("Some stages occurring in the dataset do not match any characteristics in the input stageframe.", 
               .call = FALSE)
        } else if (length(choicestage) > 1) {
          stop("Some stages in the input stageframe appear to have the same description. Please make sure that all stages included in the stageframe are defined with unique sets of characteristics.", .call = FALSE)
        }
        
        return(ltdframe$stagenames[choicestage])
      } else return("NotAlive")
    })
    
    popdatareal$stage2 <- apply(as.matrix(c(1:dim(popdatareal)[1])), 1, function(X) {
      if (is.na(popdatareal[X, stagesizecol2])) {
        popdatareal[X, stagesizecol2] <- 0
      }
      mainstages <- intersect(which(ltdframe$sizebin_min < popdatareal[X, stagesizecol2]), 
                              which(ltdframe$sizebin_max >= popdatareal[X, stagesizecol2]))
      jmstages <- which(ltdframe$immstatus == popdatareal$juvgiven2[X])
      obsstages <- which(ltdframe$obsstatus == popdatareal$obsstatus2[X])
      repstages <- which(ltdframe$repstatus == popdatareal$repstatus2[X])
      
      if (!NRasRep) {
        choicestage <- intersect(intersect(mainstages, jmstages), intersect(obsstages, repstages))
      } else {
        choicestage <- intersect(intersect(mainstages, jmstages), obsstages)
      }
      
      if (all(is.na(choicestage))) {
        stop("Some stages occurring in the dataset do not match any characteristics in the input stageframe.", 
             .call = FALSE)
      } else if (length(choicestage) > 1) {
        stop("Some stages in the input stageframe appear to have the same description. Please make sure that all stages included in the stageframe are defined with unique sets of characteristics.", .call = FALSE)
      }
      
      return(ltdframe$stagenames[choicestage])
    })
    
    popdatareal$stage3 <- apply(as.matrix(c(1:dim(popdatareal)[1])), 1, function(X) {
      if (popdatareal$alive3[X] == 1) {
        if (is.na(popdatareal[X, stagesizecol3])) {
          popdatareal[X, stagesizecol3] <- 0
        }
        mainstages <- intersect(which(ltdframe$sizebin_min < popdatareal[X, stagesizecol3]), 
                                which(ltdframe$sizebin_max >= popdatareal[X, stagesizecol3]))
        jmstages <- which(ltdframe$immstatus == popdatareal$juvgiven3[X])
        obsstages <- which(ltdframe$obsstatus == popdatareal$obsstatus3[X])
        repstages <- which(ltdframe$repstatus == popdatareal$repstatus3[X])
        
        if (!NRasRep) {
          choicestage <- intersect(intersect(mainstages, jmstages), intersect(obsstages, repstages))
        } else {
          choicestage <- intersect(intersect(mainstages, jmstages), obsstages)
        }
        
        if (all(is.na(choicestage))) {
          stop("Some stages occurring in the dataset do not match any characteristics in the input stageframe.", 
               .call = FALSE)
        } else if (length(choicestage) > 1) {
          stop("Some stages in the input stageframe appear to have the same description. Please make sure that all stages included in the stageframe are defined with unique sets of characteristics.", .call = FALSE)
        }
        
        return(ltdframe$stagenames[choicestage])
      } else return("Dead")
    })
  }
  
  popdatareal$matstatus1 <- 1 - popdatareal$juvgiven1
  
  popdatareal$matstatus2 <- 1 - popdatareal$juvgiven2
  
  popdatareal$matstatus3 <- 1 - popdatareal$juvgiven3
  
  if (any(!is.na(popdatareal$popid))) {popdatareal$popid <- as.factor(popdatareal$popid)}
  
  if (any(!is.na(popdatareal$patchid))) {popdatareal$patchid <- as.factor(popdatareal$patchid)}
  
  if (censor) {
    if (!is.na(censorkeep)) {
      popdatareal <- subset(popdatareal, censor1 == censorkeep)
      popdatareal <- subset(popdatareal, censor2 == censorkeep)
      popdatareal <- subset(popdatareal, censor3 == censorkeep)
    } else {
      popdatareal <- subset(popdatareal, is.na(censor1))
      popdatareal <- subset(popdatareal, is.na(censor2))
      popdatareal <- subset(popdatareal, is.na(censor3))
    }
  }
  
  if (NAas0) {
    popdatareal$sizea1 <- popdatareal$sizea1c
    popdatareal$sizea2 <- popdatareal$sizea2c
    popdatareal$sizea3 <- popdatareal$sizea3c
    
    popdatareal$sizeb1 <- popdatareal$sizeb1c
    popdatareal$sizeb2 <- popdatareal$sizeb2c
    popdatareal$sizeb3 <- popdatareal$sizeb3c
    
    popdatareal$sizec1 <- popdatareal$sizec1c
    popdatareal$sizec2 <- popdatareal$sizec2c
    popdatareal$sizec3 <- popdatareal$sizec3c
    
    popdatareal$repstra1 <- popdatareal$repstra1c
    popdatareal$repstra2 <- popdatareal$repstra2c
    popdatareal$repstra3 <- popdatareal$repstra3c
    
    popdatareal$repstrb1 <- popdatareal$repstrb1c
    popdatareal$repstrb2 <- popdatareal$repstrb2c
    popdatareal$repstrb3 <- popdatareal$repstrb3c
    
    popdatareal$feca1 <- popdatareal$feca1c
    popdatareal$feca2 <- popdatareal$feca2c
    popdatareal$feca3 <- popdatareal$feca3c
    
    popdatareal$fecb1 <- popdatareal$fecb1c
    popdatareal$fecb2 <- popdatareal$fecb2c
    popdatareal$fecb3 <- popdatareal$fecb3c
  }
  
  popdatareal <- popdatareal[,-c(which(names(popdatareal) == "sizea1c"), which(names(popdatareal) == "sizea2c"), 
                                 which(names(popdatareal) == "sizea3c"), which(names(popdatareal) == "sizeb1c"), 
                                 which(names(popdatareal) == "sizeb2c"), which(names(popdatareal) == "sizeb3c"), 
                                 which(names(popdatareal) == "sizec1c"), which(names(popdatareal) == "sizec2c"), 
                                 which(names(popdatareal) == "sizec3c"), which(names(popdatareal) == "repstra1c"),
                                 which(names(popdatareal) == "repstra2c"), which(names(popdatareal) == "repstra3c"), 
                                 which(names(popdatareal) == "repstrb1c"), which(names(popdatareal) == "repstrb2c"), 
                                 which(names(popdatareal) == "repstrb3c"), which(names(popdatareal) == "feca1c"), 
                                 which(names(popdatareal) == "feca2c"), which(names(popdatareal) == "feca3c"), 
                                 which(names(popdatareal) == "fecb1c"), which(names(popdatareal) == "fecb2c"), 
                                 which(names(popdatareal) == "fecb3c"))]
  
  #A little clean-up
  popdatareal <- popdatareal[,-c(which(names(popdatareal) == "alivegiven1"))]
  popdatareal <- popdatareal[,-c(which(names(popdatareal) == "alivegiven2"))]
  popdatareal <- popdatareal[,-c(which(names(popdatareal) == "alivegiven3"))]
  popdatareal <- popdatareal[,-c(which(names(popdatareal) == "deadgiven1"))]
  popdatareal <- popdatareal[,-c(which(names(popdatareal) == "deadgiven2"))]
  popdatareal <- popdatareal[,-c(which(names(popdatareal) == "deadgiven3"))]
  popdatareal <- popdatareal[,-c(which(names(popdatareal) == "obsgiven1"))]
  popdatareal <- popdatareal[,-c(which(names(popdatareal) == "obsgiven2"))]
  popdatareal <- popdatareal[,-c(which(names(popdatareal) == "obsgiven3"))]
  popdatareal <- popdatareal[,-c(which(names(popdatareal) == "nonobsgiven1"))]
  popdatareal <- popdatareal[,-c(which(names(popdatareal) == "nonobsgiven2"))]
  popdatareal <- popdatareal[,-c(which(names(popdatareal) == "nonobsgiven3"))]
  popdatareal <- popdatareal[,-c(which(names(popdatareal) == "juvgiven1"))]
  popdatareal <- popdatareal[,-c(which(names(popdatareal) == "juvgiven2"))]
  popdatareal <- popdatareal[,-c(which(names(popdatareal) == "juvgiven3"))]
  
  if (!is.na(spacing)) {
    popdatareal$density <- .density3(popdatareal, which(names(popdatareal) == "xpos2"),
                                     which(names(popdatareal) == "ypos2"),
                                     which(names(popdatareal) == "year2"), spacing)
  }
  
  
  if (reduce) {
    if (all(is.na(popdatareal$xpos1)) | length(unique(popdatareal$xpos1)) == 1) {popdatareal <- popdatareal[,-c(which(names(popdatareal) =="xpos1"))]}
    if (all(is.na(popdatareal$ypos1)) | length(unique(popdatareal$ypos1)) == 1) {popdatareal <- popdatareal[,-c(which(names(popdatareal) =="ypos1"))]}
    if (all(is.na(popdatareal$xpos2)) | length(unique(popdatareal$xpos2)) == 1) {popdatareal <- popdatareal[,-c(which(names(popdatareal) =="xpos2"))]}
    if (all(is.na(popdatareal$ypos2)) | length(unique(popdatareal$ypos2)) == 1) {popdatareal <- popdatareal[,-c(which(names(popdatareal) =="ypos2"))]}
    if (all(is.na(popdatareal$xpos3)) | length(unique(popdatareal$xpos3)) == 1) {popdatareal <- popdatareal[,-c(which(names(popdatareal) =="xpos3"))]}
    if (all(is.na(popdatareal$ypos3)) | length(unique(popdatareal$ypos3)) == 1) {popdatareal <- popdatareal[,-c(which(names(popdatareal) =="ypos3"))]}
    if (all(is.na(popdatareal$xcorr)) | length(unique(popdatareal$xcorr)) == 1) {popdatareal <- popdatareal[,-c(which(names(popdatareal) =="xcorr"))]}
    if (all(is.na(popdatareal$ycorr)) | length(unique(popdatareal$ycorr)) == 1) {popdatareal <- popdatareal[,-c(which(names(popdatareal) =="ycorr"))]}
    
    if (!is.na(censorkeep)) {
      if (censorcol > 0 & censor) {
        if (all(popdatareal$censor1 == popdatareal$censor1[1])) {
          popdatareal <- popdatareal[,-c(which(names(popdatareal) =="censor1"))]
        }
      }
      if (censorcol > 0 & censor) {
        if (all(popdatareal$censor2 == popdatareal$censor2[1])) {
          popdatareal <- popdatareal[,-c(which(names(popdatareal) =="censor2"))]
        }
      }
      if (censorcol > 0 & censor) {
        if (all(popdatareal$censor3 == popdatareal$censor3[1])) {
          popdatareal <- popdatareal[,-c(which(names(popdatareal) =="censor3"))]
        }
      }
    } else {
      if (censorcol > 0 & censor) {
        if (all(is.na(popdatareal$censor1))) {
          popdatareal <- popdatareal[,-c(which(names(popdatareal) =="censor1"))]
        }
      }
      if (censorcol > 0 & censor) {
        if (all(is.na(popdatareal$censor2))) {
          popdatareal <- popdatareal[,-c(which(names(popdatareal) =="censor2"))]
        }
      }
      if (censorcol > 0 & censor) {
        if (all(is.na(popdatareal$censor3))) {
          popdatareal <- popdatareal[,-c(which(names(popdatareal) =="censor3"))]
        }
      }
    }
    
    if (all(is.na(popdatareal$sizea1)) | length(unique(popdatareal$sizea1)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="sizea1"))]
    }
    if (all(is.na(popdatareal$sizea2)) | length(unique(popdatareal$sizea2)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="sizea2"))]
    }
    if (all(is.na(popdatareal$sizea3)) | length(unique(popdatareal$sizea3)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="sizea3"))]
    }
    
    if (all(is.na(popdatareal$sizeb1)) | length(unique(popdatareal$sizeb1)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="sizeb1"))]
    }
    if (all(is.na(popdatareal$sizeb2)) | length(unique(popdatareal$sizeb2)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="sizeb2"))]
    }
    if (all(is.na(popdatareal$sizeb3)) | length(unique(popdatareal$sizeb3)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="sizeb3"))]
    }
    
    if (all(is.na(popdatareal$sizec1)) | length(unique(popdatareal$sizec1)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="sizec1"))]
    }
    if (all(is.na(popdatareal$sizec2)) | length(unique(popdatareal$sizec2)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="sizec2"))]
    }
    if (all(is.na(popdatareal$sizec3)) | length(unique(popdatareal$sizec3)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="sizec3"))]
    }
    
    if (isTRUE(all.equal(popdatareal$size1added, popdatareal$sizea1))) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "size1added"))]
    } else if (all(is.na(popdatareal$size1added)) | length(unique(popdatareal$size1added)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "size1added"))]
    }
    if (isTRUE(all.equal(popdatareal$size2added, popdatareal$sizea2))) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "size2added"))]
    } else if (all(is.na(popdatareal$size2added)) | length(unique(popdatareal$size2added)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "size2added"))]
    }
    if (isTRUE(all.equal(popdatareal$size3added, popdatareal$sizea3))) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "size3added"))]
    } else if (all(is.na(popdatareal$size3added)) | length(unique(popdatareal$size3added)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "size3added"))]
    }
    
    if (all(is.na(popdatareal$repstra1)) | length(unique(popdatareal$repstra1)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "repstra1"))]
    }
    if (all(is.na(popdatareal$repstra2)) | length(unique(popdatareal$repstra2)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "repstra2"))]
    }
    if (all(is.na(popdatareal$repstra3)) | length(unique(popdatareal$repstra3)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "repstra3"))]
    }
    
    if (all(is.na(popdatareal$repstrb1)) | length(unique(popdatareal$repstrb1)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "repstrb1"))]
    }
    if (all(is.na(popdatareal$repstrb2)) | length(unique(popdatareal$repstrb2)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "repstrb2"))]
    }
    if (all(is.na(popdatareal$repstrb3)) | length(unique(popdatareal$repstrb3)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) == "repstrb3"))]
    }
    
    if (isTRUE(all.equal(popdatareal$repstr1added, popdatareal$repstr1a))) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="repstr1added"))]
    } else if (all(is.na(popdatareal$repstr1added)) | length(unique(popdatareal$repstr1added)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="repstr1added"))]
    }
    if (isTRUE(all.equal(popdatareal$repstr2added, popdatareal$repstr2a))) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="repstr2added"))]
    } else if (all(is.na(popdatareal$repstr2added)) | length(unique(popdatareal$repstr2added)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="repstr2added"))]
    }
    if (isTRUE(all.equal(popdatareal$repstr3added, popdatareal$repstr3a))) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="repstr3added"))]
    } else if (all(is.na(popdatareal$repstr3added)) | length(unique(popdatareal$repstr3added)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="repstr3added"))]
    }
    
    if (all(is.na(popdatareal$feca1)) | length(unique(popdatareal$feca1)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="feca1"))]
    }
    if (all(is.na(popdatareal$feca2)) | length(unique(popdatareal$feca2)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="feca2"))]
    }
    if (all(is.na(popdatareal$feca3)) | length(unique(popdatareal$feca3)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="feca3"))]
    }
    
    if (all(is.na(popdatareal$fecb1)) | length(unique(popdatareal$fecb1)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fecb1"))]
    }
    if (all(is.na(popdatareal$fecb2)) | length(unique(popdatareal$fecb2)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fecb2"))]
    }
    if (all(is.na(popdatareal$fecb3)) | length(unique(popdatareal$fecb3)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fecb3"))]
    }
    
    if (isTRUE(all.equal(popdatareal$fec1added, popdatareal$feca1))) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fec1added"))]
    } else if (all(is.na(popdatareal$fec1added)) | length(unique(popdatareal$fec1added)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fec1added"))]
    }
    if (isTRUE(all.equal(popdatareal$fec2added, popdatareal$feca2))) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fec2added"))]
    } else if (all(is.na(popdatareal$fec2added)) | length(unique(popdatareal$fec2added)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fec2added"))]
    }
    if (isTRUE(all.equal(popdatareal$fec3added, popdatareal$feca3))) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fec3added"))]
    } else if (all(is.na(popdatareal$fec3added)) | length(unique(popdatareal$fec3added)) == 1) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fec3added"))]
    }
    
    if (all(popdatareal$obsstatus1 == popdatareal$obsstatus1[1])) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="obsstatus1"))]
    }
    if (all(popdatareal$obsstatus2 == popdatareal$obsstatus2[1])) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="obsstatus2"))]
    }
    if (all(popdatareal$obsstatus3 == popdatareal$obsstatus3[1])) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="obsstatus3"))]
    }
    
    if (all(popdatareal$repstatus1 == popdatareal$repstatus1[1])) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="repstatus1"))]
    }
    if (all(popdatareal$repstatus2 == popdatareal$repstatus2[1])) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="repstatus2"))]
    }
    if (all(popdatareal$repstatus3 == popdatareal$repstatus3[1])) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="repstatus3"))]
    }
    
    if (all(popdatareal$fecstatus1 == popdatareal$fecstatus1[1])) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fecstatus1"))]
    }
    if (all(popdatareal$fecstatus2 == popdatareal$fecstatus2[1])) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fecstatus2"))]
    }
    if (all(popdatareal$fecstatus3 == popdatareal$fecstatus3[1])) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="fecstatus3"))]
    }
    
    if (!censor) {
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="censor1"))]
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="censor2"))]
      popdatareal <- popdatareal[,-c(which(names(popdatareal) =="censor3"))]
    }
    
  }
  
  class(popdatareal) <- append(class(popdatareal), "hfvdata")
  
  return(popdatareal)
}

#' Create Historical Vertical Data Frame From Ahistorical Vertical Data Frame
#' 
#' \code{historicalize3()} returns a vertically formatted demographic data frame
#' organized to create historical projection matrices, given a vertically but
#' ahistorically formatted data frame. This data frame is in standard 'lefko3'
#' format and can be used in all functions in the package.
#'
#' @param data The horizontal data file.
#' @param popidcol A variable name or column number corresponding to the identity
#' of the population for each individual.
#' @param patchidcol A variable name or column number corresponding to the
#' identity of the patch for each individual, if patches have been designated
#' within populations.
#' @param individcol A variable name or column number corresponding to the
#' identity of each individual.
#' @param year2col A variable name or column number corresponding to the year
#' or time step in time \emph{t}.
#' @param year3col A variable name or column number corresponding to the year
#' or time step in time \emph{t}+1.
#' @param xcol A variable name or column number corresponding to the x 
#' coordinate of each individual in Cartesian space.
#' @param ycol A variable name or column number corresponding to the y 
#' coordinate of each individual in Cartesian space.
#' @param sizea2col A variable name or column number corresponding to
#' the primary size entry in time \emph{t}.
#' @param sizea3col A variable name or column number corresponding to
#' the primary size entry in time \emph{t}+1.
#' @param sizeb2col A variable name or column number corresponding to
#' the secondary size entry in time \emph{t}.
#' @param sizeb3col A variable name or column number corresponding to
#' the secondary size entry in time \emph{t}+1.
#' @param sizec2col A variable name or column number corresponding to
#' the tertiary size entry in time \emph{t}.
#' @param sizec3col A variable name or column number corresponding to
#' the tertiary size entry in time \emph{t}+1.
#' @param repstra2col A variable name or column number corresponding to the
#' production of reproductive structures, such as flowers, in time \emph{t}. 
#' This can be binomial or count data, and is used to in analysis of the 
#' probability of reproduction.
#' @param repstra3col A variable name or column number corresponding to the
#' production of reproductive structures, such as flowers, in time \emph{t}+1. 
#' This can be binomial or count data, and is used to in analysis of the 
#' probability of reproduction.
#' @param repstrb2col A second variable name or column number corresponding 
#' to the production of reproductive structures, such as flowers, in time
#' \emph{t}. This can be binomial or count data.
#' @param repstrb3col A second variable name or column number corresponding 
#' to the production of reproductive structures, such as flowers, in time
#' \emph{t}+1. This can be binomial or count data.
#' @param feca2col A variable name or column number corresponding to fecundity
#' in time \emph{t}. This may represent egg counts, fruit counts, seed 
#' production, etc.
#' @param feca3col A variable name or column number corresponding to fecundity
#' in time \emph{t}+1. This may represent egg counts, fruit counts, seed 
#' production, etc.
#' @param fecb2col A second variable name or column number corresponding to 
#' fecundity in time \emph{t}. This may represent egg counts, fruit counts, 
#' seed production, etc.
#' @param fecb3col A second variable name or column number corresponding to 
#' fecundity in time \emph{t}+1. This may represent egg counts, fruit counts, 
#' seed production, etc.
#' @param alive2col A variable name or column number that provides information
#' an whether an individual is alive in time \emph{t}. If used, living status
#' must be designated as binomial (living = 1, dead = 0).
#' @param alive3col A variable name or column number that provides information
#' an whether an individual is alive in time \emph{t}+1. If used, living status
#' must be designated as binomial (living = 1, dead = 0).
#' @param dead2col A variable name or column number that provides information
#' an whether an individual is dead in time \emph{t}. If used, dead status
#' must be designated as binomial (living = 0, dead = 1).
#' @param dead3col A variable name or column number that provides information
#' an whether an individual is dead in time \emph{t}+1. If used, dead status
#' must be designated as binomial (living = 0, dead = 1).
#' @param obs2col A variable name or column number providing information on
#' whether an individual is in an observable stage in time \emph{t}. If used,
#' observation status must be designated as binomial (observed = 1, not 
#' observed = 0).
#' @param obs3col A variable name or column number providing information on
#' whether an individual is in an observable stage in time \emph{t}+1. If
#' used, observation status must be designated as binomial (observed = 1, not 
#' observed = 0).
#' @param nonobs2col A variable name or column number providing information on
#' whether an individual is in an unobservable stage in time \emph{t}. If
#' used, observation status must be designated as binomial (observed = 0, not 
#' observed = 1).
#' @param nonobs3col A variable name or column number providing information on
#' whether an individual is in an unobservable stage in time \emph{t}+1. If
#' used, observation status must be designated as binomial (observed = 0, not 
#' observed = 1).
#' @param repstrrel This is a scalar modifier for that makes the variable in
#' \code{repstrb2col} equivalent to \code{repstra2col}. This can be useful if two 
#' reproductive status variables have related but unequal units, for example
#' if \code{repstrb2col} refers to one-flowered stems while \code{repstra2col} refers to
#' two-flowered stems.
#' @param fecrel This is a scalar modifier for that makes the variable in
#' \code{fecb2col} equivalent to \code{feca2col}. This can be useful if two fecundity
#' variables have related but unequal units.
#' @param stage2col A variable name or column number corresponding to
#' life history stage in time \emph{t}.
#' @param stage3col A variable name or column number corresponding to
#' life history stage in time \emph{t}+1.
#' @param juv2col A variable name or column number that marks individuals in
#' immature stages in time \emph{t}. The \code{historicalize3()} function assumes 
#' that immature individuals are identified in this variable marked with a 
#' number equal to or greater than 1, and that mature individuals are marked 
#' as 0 or NA.
#' @param juv3col A variable name or column number that marks individuals in
#' immature stages in time \emph{t}+1. The \code{historicalize3()} function assumes 
#' that immature individuals are identified in this variable marked with a 
#' number equal to or greater than 1, and that mature individuals are marked 
#' as 0 or NA.
#' @param stageassign The stageframe object identifying the life history model
#' being operationalized.
#' @param stagesize A variable name or column number describing which size
#' variable to use in stage estimation. Defaults to NA, and can also take
#' \code{sizea}, \code{sizeb}, \code{sizec}, or \code{sizeadded}, depending on which size variable 
#' is chosen.
#' @param censorcol A variable name or column number corresponding to a censor
#' variable within the dataset, used to distinguish between entries to use and
#' those to discard from analysis, or to designate entries with special issues 
#' that require further attention.
#' @param censorkeep The value of the censoring variable identifying data
#' that should be included in analysis. Defaults to 1, but may take any value
#' including NA.
#' @param censor A logical variable determining whether the output data
#' should be censored using the variable defined in \code{censorcol}. Defaults
#' to FALSE.
#' @param spacing The spacing at which density should be estimated, if density
#' estimation is desired and x and y coordinates are supplied. Given in the
#' same units as those used in the x and y coordinates given in \code{xcol} and
#' \code{ycol}. Defaults to NA.
#' @param NAas0 If TRUE, then all NA entries for size and fecundity variables
#' will be set to 0. This can help increase the sample size analyzed by
#' \code{\link{modelsearch}()}, but should only be used when it is clear that this
#' substitution is biologically realistic. Defaults to FALSE.
#' @param NRasRep If TRUE, then will treat non-reproductive but mature
#' individuals as reproductive during stage assignment. This can be useful,
#' for example, when a matrix is desired without separation of reproductive
#' and non-reproductive but mature stages of the same size. Only used if
#' \code{stageassign} is set to a stageframe. Defaults to FALSE.
#' @param reduce A logical variable determining whether invariant state
#' variables should be removed from the output dataset. For example, if
#' all living individuals are always observable, then variables identifying
#' observation status are invariant and can be removed. Defaults to TRUE.
#'
#' @return If all inputs are properly formatted, then this function will output
#' a historical vertical data frame (class \code{hfvdata}), meaning that the
#' output data frame will have three consecutive years of size and reproductive 
#' data per individual per row.
#' \item{popid}{Unique identifier for the population, if given.}
#' \item{patchid}{Unique identifier for the patch within the population, if
#' given.}
#' \item{year1, year2, year3}{Year or time step at times \emph{t}-1, \emph{t}, 
#' and \emph{t}+1.}
#' \item{xpos,ypos}{X and Y position in Cartesian space, if given.}
#' \item{sizea1,sizea2,sizea3}{Main size measurement in times \emph{t}-1, 
#' \emph{t}, and \emph{t}+1, respectively.}
#' \item{sizeb1,sizeb2,sizeb3}{Secondary size measurement in times \emph{t}-1, 
#' \emph{t}, and \emph{t}+1, respectively.}
#' \item{sizec1,sizec2,sizec3}{Tertiary size measurement in times \emph{t}-1, 
#' \emph{t}, and \emph{t}+1, respectively.}
#' \item{censor1,censor2,censor3}{Censor state values in times \emph{t}-1,
#' \emph{t}, and \emph{t}+1.}
#' \item{repstra1,repstrb1,repstrc1}{Main, secondary, and tertiary numbers of 
#' reproductive structures in time \emph{t}-1.}
#' \item{repstra2,repstrb2,repstrc2}{Main, secondary, and tertiary numbers of 
#' reproductive structures in time \emph{t}.}
#' \item{repstra3,repstrb3,repstrc3}{Main, secondary, and tertiary numbers of 
#' reproductive structures in time \emph{t}+1.}
#' \item{feca1,fecb1}{Main, secondary, and tertiary numbers of offspring in 
#' time \emph{t}-1.}
#' \item{feca2,fecb2}{Main, secondary, and tertiary numbers of offspring in 
#' time \emph{t}.}
#' \item{feca3,fecb3}{Main, secondary, and tertiary numbers of offspring in 
#' time \emph{t}+1.}
#' \item{size1added,size2added,size3added}{Sum of primary, secondary, and
#' tertiary size measurements in times \emph{t}-1, \emph{t}, and \emph{t}+1, 
#' respectively.}
#' \item{repstr1added,repstr2added,repstr3added}{Sum of primary, secondary, and
#' tertiary reproductive structures in times \emph{t}-1, \emph{t}, and 
#' \emph{t}+1, respectively.}
#' \item{fec1added,fec2added,fec3added}{Sum of primary, secondary, and
#' tertiary offspring numbers in times \emph{t}-1, \emph{t}, and \emph{t}+1, 
#' respectively.}
#' \item{obsstatus1,obsstatus2,obsstatus3}{Binomial observation state in 
#' times \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{repstatus1,repstatus2,repstatus3}{Binomial reproductive state in 
#' times \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{fecstatus1,fecstatus2,fecstatus3}{Binomial offspring production
#' state in times \emph{t}-1, \emph{t}, and \emph{t}+1, respectively.}
#' \item{firstseen}{Year or time step of first observation.}
#' \item{lastseen}{Year or time step of last observation.}
#' \item{xcorr,ycorr}{Overall x and y coordinates of individual in
#' Cartesian space.}
#' \item{alive1,alive2,alive3}{Binomial state as alive in times \emph{t}-1, 
#' \emph{t}, and \emph{t}+1, respectively.}
#' \item{obsage}{Observed age in time \emph{t}, assuming first observation
#' corresponds to age = 0.}
#' \item{obslifespan}{Observed lifespan, given as `lastseen - firstseen + 1`.}
#' \item{individ}{Unique identifier for the individual.}
#' \item{density}{Density of individuals per unit designated in `spacing`.
#' Only given if spacing is not NA.}
#' 
#' @examples
#' data(cypvert)
#' 
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg", "XLg")
#' repvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' obsvector <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)
#' matvector <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' immvector <- c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
#' propvector <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' indataset <- c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' binvec <- c(0, 0, 0, 0, 0, 0.5, 0.5, 1, 1, 2.5, 7)
#' 
#' cypframe_raw <- sf_create(sizes = sizevector, stagenames = stagevector, 
#'                           repstatus = repvector, obsstatus = obsvector,
#'                           matstatus = matvector, propstatus = propvector,
#'                           immstatus = immvector, indataset = indataset,
#'                           binhalfwidth = binvec)
#' 
#' cypframe_raw
#' 
#' cypraw_v2 <- historicalize3(data = cypvert, patchidcol = "patch", individcol = "plantid",
#'                             year2col = "year2", sizea2col = "Inf2.2", sizea3col = "Inf2.3",
#'                             sizeb2col = "Inf.2", sizeb3col = "Inf.3", sizec2col = "Veg.2",
#'                             sizec3col = "Veg.3", repstra2col = "Inf2.2", repstra3col = "Inf2.3",
#'                             repstrb2col = "Inf.2", repstrb3col = "Inf.3", feca2col = "Pod.2",
#'                             feca3col = "Pod.3", repstrrel = 2, stageassign = cypframe_raw,
#'                             stagesize = "sizeadded", censorcol = "censor", censor = FALSE,
#'                             NAas0 = TRUE, NRasRep = TRUE, reduce = TRUE)
#' summary(cypraw_v2)

#' @export
historicalize3 <- function(data, popidcol = 0, patchidcol = 0, individcol, year2col = 0, 
                           year3col = 0, xcol = 0, ycol = 0, sizea2col = 0, sizea3col = 0, 
                           sizeb2col = 0, sizeb3col = 0, sizec2col = 0, sizec3col = 0,
                           repstra2col = 0, repstra3col = 0, repstrb2col = 0, 
                           repstrb3col = 0, feca2col = 0, feca3col = 0,  fecb2col = 0, 
                           fecb3col = 0, alive2col = 0, alive3col = 0, dead2col = 0, 
                           dead3col = 0, obs2col = 0, obs3col = 0, nonobs2col = 0, 
                           nonobs3col = 0, repstrrel = 1, fecrel = 1, stage2col = 0, 
                           stage3col = 0, juv2col = 0, juv3col = 0, stageassign = NA, 
                           stagesize = NA, censorcol = 0, censorkeep = 1, 
                           censor = FALSE, spacing = NA, NAas0 = FALSE, NRasRep = FALSE, 
                           reduce = TRUE) {
  
  alive2 <- indataset <- censor1 <- censor2 <- censor3 <- NULL
  
  if (is.na(individcol)) {
    stop("Individual ID variable is required.", .call = FALSE)
  }
  
  if (is.na(year2col) & is.na(year3col)) {
    stop("Variable identifying either year2 (time t) or year3 (time t+1) is required.", .call = FALSE)
  }
  
  if (is.character(popidcol)) {
    if (is.element(popidcol, names(data))) {
      true.popidcol <- which(names(data) == popidcol)
      popidcol <- true.popidcol
    } else {
      stop("Please enter popidcol exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(patchidcol)) {
    if (is.element(patchidcol, names(data))) {
      true.patchidcol <- which(names(data) == patchidcol)
      patchidcol <- true.patchidcol
    } else {
      stop("Please enter patchidcol exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(individcol)) {
    if (is.element(individcol, names(data))) {
      true.individcol <- which(names(data) == individcol)
      individcol <- true.individcol
    } else {
      stop("Please enter individcol exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (year2col != 0) {
    if (is.character(year2col)) {
      if (is.element(year2col, names(data))) {
        true.year2col <- which(names(data) == year2col)
        year2col <- true.year2col
      } else {
        stop("Please enter year2col exactly as it appears in the dataset.", call. = FALSE)
      }
    }
  } 
  
  if (year3col != 0) {
    if (is.character(year3col)) {
      if (is.element(year3col, names(data))) {
        true.year3col <- which(names(data) == year3col)
        year3col <- true.year3col
      } else {
        stop("Please enter year3col exactly as it appears in the dataset.", call. = FALSE)
      }
    }
  }
  
  if (year2col != 0 & year3col == 0) {
    data$year3 <- data[,year2col] + 1
    year3col <- which(names(data) == "year3")
  } else if (year2col == 0 & year3col != 0) {
    data$year2 <- data[,year3col] - 1
    year2col <- which(names(data) == "year2")
  }
  
  if (is.character(xcol)) {
    if (is.element(xcol, names(data))) {
      true.xcol <- which(names(data) == xcol)
      xcol <- true.xcol
    } else {
      stop("Please enter xcol exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(ycol)) {
    if (is.element(ycol, names(data))) {
      true.ycol <- which(names(data) == ycol)
      ycol <- true.ycol
    } else {
      stop("Please enter ycol exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(sizea2col)) {
    if (is.element(sizea2col, names(data))) {
      true.sizea2col <- which(names(data) == sizea2col)
      sizea2col <- true.sizea2col
    } else {
      stop("Please enter sizea2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(sizea3col)) {
    if (is.element(sizea3col, names(data))) {
      true.sizea3col <- which(names(data) == sizea3col)
      sizea3col <- true.sizea3col
    } else {
      stop("Please enter sizea3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(sizeb2col)) {
    if (is.element(sizeb2col, names(data))) {
      true.sizeb2col <- which(names(data) == sizeb2col)
      sizeb2col <- true.sizeb2col
    } else {
      stop("Please enter sizeb2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(sizeb3col)) {
    if (is.element(sizeb3col, names(data))) {
      true.sizeb3col <- which(names(data) == sizeb3col)
      sizeb3col <- true.sizeb3col
    } else {
      stop("Please enter sizeb3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(sizec2col)) {
    if (is.element(sizec2col, names(data))) {
      true.sizec2col <- which(names(data) == sizec2col)
      sizec2col <- true.sizec2col
    } else {
      stop("Please enter sizec2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(sizec3col)) {
    if (is.element(sizec3col, names(data))) {
      true.sizec3col <- which(names(data) == sizec3col)
      sizec3col <- true.sizec3col
    } else {
      stop("Please enter sizec3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(repstra2col)) {
    if (is.element(repstra2col, names(data))) {
      true.repstra2col <- which(names(data) == repstra2col)
      repstra2col <- true.repstra2col
    } else {
      stop("Please enter repstra2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(repstra3col)) {
    if (is.element(repstra3col, names(data))) {
      true.repstra3col <- which(names(data) == repstra3col)
      repstra3col <- true.repstra3col
    } else {
      stop("Please enter repstra3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(repstrb2col)) {
    if (is.element(repstrb2col, names(data))) {
      true.repstrb2col <- which(names(data) == repstrb2col)
      repstrb2col <- true.repstrb2col
    } else {
      stop("Please enter repstrb2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(repstrb3col)) {
    if (is.element(repstrb3col, names(data))) {
      true.repstrb3col <- which(names(data) == repstrb3col)
      repstrb3col <- true.repstrb3col
    } else {
      stop("Please enter repstrb3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(feca2col)) {
    if (is.element(feca2col, names(data))) {
      true.feca2col <- which(names(data) == feca2col)
      feca2col <- true.feca2col
    } else {
      stop("Please enter feca2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(feca3col)) {
    if (is.element(feca3col, names(data))) {
      true.feca3col <- which(names(data) == feca3col)
      feca3col <- true.feca3col
    } else {
      stop("Please enter feca3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(fecb2col)) {
    if (is.element(fecb2col, names(data))) {
      true.fecb2col <- which(names(data) == fecb2col)
      fecb2col <- true.fecb2col
    } else {
      stop("Please enter fecb2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(fecb3col)) {
    if (is.element(fecb3col, names(data))) {
      true.fecb3col <- which(names(data) == fecb3col)
      fecb3col <- true.fecb3col
    } else {
      stop("Please enter fecb3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(alive2col)) {
    if (is.element(alive2col, names(data))) {
      true.alive2col <- which(names(data) == alive2col)
      alive2col <- true.alive2col
    } else {
      stop("Please enter alive2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(alive3col)) {
    if (is.element(alive3col, names(data))) {
      true.alive3col <- which(names(data) == alive3col)
      alive3col <- true.alive3col
    } else {
      stop("Please enter alive3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(dead2col)) {
    if (is.element(dead2col, names(data))) {
      true.dead2col <- which(names(data) == dead2col)
      dead2col <- true.dead2col
    } else {
      stop("Please enter dead2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(dead3col)) {
    if (is.element(dead3col, names(data))) {
      true.dead3col <- which(names(data) == dead3col)
      dead3col <- true.dead3col
    } else {
      stop("Please enter dead3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(obs2col)) {
    if (is.element(obs2col, names(data))) {
      true.obs2col <- which(names(data) == obs2col)
      obs2col <- true.obs2col
    } else {
      stop("Please enter obs2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(obs3col)) {
    if (is.element(obs3col, names(data))) {
      true.obs3col <- which(names(data) == obs3col)
      obs3col <- true.obs3col
    } else {
      stop("Please enter obs3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(nonobs2col)) {
    if (is.element(nonobs2col, names(data))) {
      true.nonobs2col <- which(names(data) == nonobs2col)
      nonobs2col <- true.nonobs2col
    } else {
      stop("Please enter nonobs2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(nonobs3col)) {
    if (is.element(nonobs3col, names(data))) {
      true.nonobs3col <- which(names(data) == nonobs3col)
      nonobs3col <- true.nonobs3col
    } else {
      stop("Please enter nonobs3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(stage2col)) {
    if (is.element(stage2col, names(data))) {
      true.stage2col <- which(names(data) == stage2col)
      stage2col <- true.stage2col
    } else {
      stop("Please enter stage2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(stage3col)) {
    if (is.element(stage3col, names(data))) {
      true.stage3col <- which(names(data) == stage3col)
      stage3col <- true.stage3col
    } else {
      stop("Please enter stage3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(juv2col)) {
    if (is.element(juv2col, names(data))) {
      true.juv2col <- which(names(data) == juv2col)
      juv2col <- true.juv2col
    } else {
      stop("Please enter juv2col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(juv3col)) {
    if (is.element(juv3col, names(data))) {
      true.juv3col <- which(names(data) == juv3col)
      juv3col <- true.juv3col
    } else {
      stop("Please enter juv3col exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (is.character(censorcol)) {
    if (is.element(censorcol, names(data))) {
      true.censorcol <- which(names(data) == censorcol)
      censorcol <- true.censorcol
    } else {
      stop("Please enter censorcol exactly as it appears in the dataset.", call. = FALSE)
    }
  }
  
  if (!is.na(spacing)) {
    if (xcol == 0 | ycol == 0) {
      stop("Density estimation cannot proceed without valid x and y coordinates.")
    }
    
    if (is.character(spacing)) {
      stop("The spacing option requires either a number, or defaults to NA.")
    }
  }
  
  id.vec <- sort(unique(data[, individcol]))
  
  norows <- dim(data)[1]
  data$rowid <- c(1:norows)
  
  full.set <- do.call(rbind, apply(as.matrix(id.vec), 1, function(X) {
    indiv.set <- .core_gather_ahvtohv(subset(data, data[,individcol] == X), popidcol, patchidcol, 
                                      year2col, year3col, xcol, ycol, sizea2col, sizea3col, 
                                      sizeb2col, sizeb3col, sizec2col, sizec3col, repstra2col, 
                                      repstra3col, repstrb2col, repstrb3col, feca2col, feca3col, 
                                      fecb2col, fecb3col, alive2col, alive3col, dead2col, dead3col,
                                      obs2col, obs3col, nonobs2col, nonobs3col, repstrrel, fecrel,
                                      stage2col, stage3col, juv2col, juv3col, censorcol)
    if (!all(is.na(indiv.set))) {
      indiv.set <- cbind.data.frame(indiv.set, X)
    } else {
      indiv.set <- rep(NA, 71)
    }
    return(indiv.set)
  }))
  
  names(full.set) <- c("popid", "patchid", "rowid",
                       
                       "year1",  "xpos1", "ypos1", "sizea1", "sizeb1", "sizec1", 
                       "repstra1", "repstrb1", "feca1", "fecb1", "alivegiven1", "deadgiven1", "obsgiven1", 
                       "nonobsgiven1", "stage1", "juvgiven1", "censor1", 
                       
                       "year2", "xpos2", "ypos2", "sizea2", "sizeb2", "sizec2", "repstra2", "repstrb2", 
                       "feca2", "fecb2", "alivegiven2", "deadgiven2", "obsgiven2", "nonobsgiven2", "stage2", 
                       "juvgiven2", "censor2", 
                       
                       "year3", "xpos3", "ypos3", "sizea3", "sizeb3", "sizec3", "repstra3", "repstrb3", 
                       "feca3", "fecb3", "alivegiven3", "deadgiven3", "obsgiven3", "nonobsgiven3", "stage3", 
                       "juvgiven3", "censor3", 
                       
                       "firstseen", "lastseen", "alive1", "alive2", "alive3", "obsstatus1", 
                       "obsstatus2", "obsstatus3", "matstatus1", "matstatus2", "matstatus3", "repstatus1",
                       "repstatus2", "repstatus3", "fecstatus1", "fecstatus2", "fecstatus3", "individ")
  
  rownames(full.set) <- c(1:dim(full.set)[1])
  
  if (any(is.na(full.set$alive1))) {
    full.set$alive1[(which(is.na(full.set$alive1)))] <- 0
  }
  
  if (any(is.na(full.set$alive2))) {
    full.set$alive2[(which(is.na(full.set$alive2)))] <- 0
  }
  
  if (any(is.na(full.set$alive3))) {
    full.set$alive3[(which(is.na(full.set$alive3)))] <- 0
  }
  
  full.set <- subset(full.set, alive2 == 1)
  
  full.set$obsage <- full.set$year2 - full.set$firstseen
  full.set$obslifespan <- full.set$lastseen - full.set$firstseen
  
  full.set$sizea1c <- full.set$sizea1
  full.set$sizea2c <- full.set$sizea2
  full.set$sizea3c <- full.set$sizea3
  full.set$sizeb1c <- full.set$sizeb1
  full.set$sizeb2c <- full.set$sizeb2
  full.set$sizeb3c <- full.set$sizeb3
  full.set$sizec1c <- full.set$sizec1
  full.set$sizec2c <- full.set$sizec2
  full.set$sizec3c <- full.set$sizec3
  
  full.set$repstra1c <- full.set$repstra1
  full.set$repstra2c <- full.set$repstra2
  full.set$repstra3c <- full.set$repstra3
  full.set$repstrb1c <- full.set$repstrb1
  full.set$repstrb2c <- full.set$repstrb2
  full.set$repstrb3c <- full.set$repstrb3
  
  full.set$feca1c <- full.set$feca1
  full.set$feca2c <- full.set$feca2
  full.set$feca3c <- full.set$feca3
  full.set$fecb1c <- full.set$fecb1
  full.set$fecb2c <- full.set$fecb2
  full.set$fecb3c <- full.set$fecb3
  
  full.set$sizea1c[which(is.na(full.set$sizea1))] <- 0
  full.set$sizea2c[which(is.na(full.set$sizea2))] <- 0
  full.set$sizea3c[which(is.na(full.set$sizea3))] <- 0
  
  full.set$repstra1c[which(is.na(full.set$repstra1))] <- 0
  full.set$repstra2c[which(is.na(full.set$repstra2))] <- 0
  full.set$repstra3c[which(is.na(full.set$repstra3))] <- 0
  
  full.set$sizeb1c[which(is.na(full.set$sizeb1))] <- 0
  full.set$sizeb2c[which(is.na(full.set$sizeb2))] <- 0
  full.set$sizeb3c[which(is.na(full.set$sizeb3))] <- 0
  
  full.set$repstrb1c[which(is.na(full.set$repstrb1))] <- 0
  full.set$repstrb2c[which(is.na(full.set$repstrb2))] <- 0
  full.set$repstrb3c[which(is.na(full.set$repstrb3))] <- 0
  
  full.set$sizec1c[which(is.na(full.set$sizec1))] <- 0
  full.set$sizec2c[which(is.na(full.set$sizec2))] <- 0
  full.set$sizec3c[which(is.na(full.set$sizec3))] <- 0
  
  full.set$feca1c[which(is.na(full.set$feca1))] <- 0
  full.set$feca2c[which(is.na(full.set$feca2))] <- 0
  full.set$feca3c[which(is.na(full.set$feca3))] <- 0
  
  full.set$fecb1c[which(is.na(full.set$fecb1))] <- 0
  full.set$fecb2c[which(is.na(full.set$fecb2))] <- 0
  full.set$fecb3c[which(is.na(full.set$fecb3))] <- 0
  
  full.set$juvgiven1[which(is.na(full.set$juvgiven1))] <- 0
  full.set$juvgiven2[which(is.na(full.set$juvgiven2))] <- 0
  full.set$juvgiven3[which(is.na(full.set$juvgiven3))] <- 0
  
  full.set$size1added <- full.set$sizea1c + full.set$sizeb1c + full.set$sizec1c
  full.set$size2added <- full.set$sizea2c + full.set$sizeb2c + full.set$sizec2c
  full.set$size3added <- full.set$sizea3c + full.set$sizeb3c + full.set$sizec3c
  
  full.set$repstr1added <- full.set$repstra1c + full.set$repstrb1c
  full.set$repstr2added <- full.set$repstra2c + full.set$repstrb2c
  full.set$repstr3added <- full.set$repstra3c + full.set$repstrb3c
  
  full.set$fec1added <- full.set$feca1c + full.set$fecb1c
  full.set$fec2added <- full.set$feca2c + full.set$fecb2c
  full.set$fec3added <- full.set$feca3c + full.set$fecb3c
  
  if(!all(is.na(stageassign))) {
    if (stagesize == "sizeadded") {
      stagesizecol1 <- which(names(full.set) == "size1added")
      stagesizecol2 <- which(names(full.set) == "size2added")
      stagesizecol3 <- which(names(full.set) == "size3added")
    } else if (stagesize == "sizec") {
      stagesizecol1 <- which(names(full.set) == "sizec1c")
      stagesizecol2 <- which(names(full.set) == "sizec2c")
      stagesizecol3 <- which(names(full.set) == "sizec3c")
    } else if (stagesize == "sizeb") {
      stagesizecol1 <- which(names(full.set) == "sizeb1c")
      stagesizecol2 <- which(names(full.set) == "sizeb2c")
      stagesizecol3 <- which(names(full.set) == "sizeb3c")
    } else {
      stagesizecol1 <- which(names(full.set) == "sizea1c")
      stagesizecol2 <- which(names(full.set) == "sizea2c")
      stagesizecol3 <- which(names(full.set) == "sizea3c")
    }
    
    ltdframe <- subset(stageassign, indataset == 1)
    ltdframe$stagenames <- as.character(ltdframe$stagenames)
    
    full.set$stage1 <- apply(as.matrix(c(1:dim(full.set)[1])), 1, function(X) {
      if (full.set$alive1[X] == 1) {
        if (is.na(full.set[X, stagesizecol1])) {
          full.set[X, stagesizecol1] <- 0
        }
        mainstages <- intersect(which(ltdframe$sizebin_min < full.set[X, stagesizecol1]), 
                                which(ltdframe$sizebin_max >= full.set[X, stagesizecol1]))
        jmstages <- which(ltdframe$immstatus == full.set$juvgiven1[X])
        obsstages <- which(ltdframe$obsstatus == full.set$obsstatus1[X])
        repstages <- which(ltdframe$repstatus == full.set$repstatus1[X])
        
        if (!NRasRep) {
          choicestage <- intersect(intersect(mainstages, jmstages), intersect(obsstages, repstages))
        } else {
          choicestage <- intersect(intersect(mainstages, jmstages), obsstages)
        }
        
        if (all(is.na(choicestage))) {
          stop("Some stages occurring in the dataset do not match any characteristics in the input stageframe.", 
               .call = FALSE)
        } else if (length(choicestage) > 1) {
          stop("Some stages in the input stageframe appear to have the same description. Please make sure that all stages included in the stageframe are defined with unique sets of characteristics.", .call = FALSE)
        }
        
        return(ltdframe$stagenames[choicestage])
      } else return("NotAlive")
    })
    
    full.set$stage2 <- apply(as.matrix(c(1:dim(full.set)[1])), 1, function(X) {
      if (is.na(full.set[X, stagesizecol2])) {
        full.set[X, stagesizecol2] <- 0
      }
      mainstages <- intersect(which(ltdframe$sizebin_min < full.set[X, stagesizecol2]), 
                              which(ltdframe$sizebin_max >= full.set[X, stagesizecol2]))
      jmstages <- which(ltdframe$immstatus == full.set$juvgiven2[X])
      obsstages <- which(ltdframe$obsstatus == full.set$obsstatus2[X])
      repstages <- which(ltdframe$repstatus == full.set$repstatus2[X])
      
      if (!NRasRep) {
        choicestage <- intersect(intersect(mainstages, jmstages), intersect(obsstages, repstages))
      } else {
        choicestage <- intersect(intersect(mainstages, jmstages), obsstages)
      }
      
      if (all(is.na(choicestage))) {
        stop("Some stages occurring in the dataset do not match any characteristics in the input stageframe.", 
             .call = FALSE)
      } else if (length(choicestage) > 1) {
        stop("Some stages in the input stageframe appear to have the same description. Please make sure that all stages included in the stageframe are defined with unique sets of characteristics.", .call = FALSE)
      }
      
      return(ltdframe$stagenames[choicestage])
    })
    
    full.set$stage3 <- apply(as.matrix(c(1:dim(full.set)[1])), 1, function(X) {
      if (full.set$alive3[X] == 1) {
        if (is.na(full.set[X, stagesizecol3])) {
          full.set[X, stagesizecol3] <- 0
        }
        mainstages <- intersect(which(ltdframe$sizebin_min < full.set[X, stagesizecol3]), 
                                which(ltdframe$sizebin_max >= full.set[X, stagesizecol3]))
        jmstages <- which(ltdframe$immstatus == full.set$juvgiven3[X])
        obsstages <- which(ltdframe$obsstatus == full.set$obsstatus3[X])
        repstages <- which(ltdframe$repstatus == full.set$repstatus3[X])
        
        if (!NRasRep) {
          choicestage <- intersect(intersect(mainstages, jmstages), intersect(obsstages, repstages))
        } else {
          choicestage <- intersect(intersect(mainstages, jmstages), obsstages)
        }
        
        if (all(is.na(choicestage))) {
          stop("Some stages occurring in the dataset do not match any characteristics in the input stageframe.", 
               .call = FALSE)
        } else if (length(choicestage) > 1) {
          stop("Some stages in the input stageframe appear to have the same description. Please make sure that all stages included in the stageframe are defined with unique sets of characteristics.", .call = FALSE)
        }
        
        return(ltdframe$stagenames[choicestage])
      } else return("Dead")
    })
  }
  
  if (censor) {
    full.set <- subset(full.set, censor1 != 0 & censor2 != 0)
    full.set <- subset(full.set, censor3 != 0)
  }
  
  if (NAas0) {
    full.set$sizea1 <- full.set$sizea1c
    full.set$sizea2 <- full.set$sizea2c
    full.set$sizea3 <- full.set$sizea3c
    
    full.set$sizeb1 <- full.set$sizeb1c
    full.set$sizeb2 <- full.set$sizeb2c
    full.set$sizeb3 <- full.set$sizeb3c
    
    full.set$sizec1 <- full.set$sizec1c
    full.set$sizec2 <- full.set$sizec2c
    full.set$sizec3 <- full.set$sizec3c
    
    full.set$repstra1 <- full.set$repstra1c
    full.set$repstra2 <- full.set$repstra2c
    full.set$repstra3 <- full.set$repstra3c
    
    full.set$repstrb1 <- full.set$repstrb1c
    full.set$repstrb2 <- full.set$repstrb2c
    full.set$repstrb3 <- full.set$repstrb3c
    
    full.set$feca1 <- full.set$feca1c
    full.set$feca2 <- full.set$feca2c
    full.set$feca3 <- full.set$feca3c
    
    full.set$fecb1 <- full.set$fecb1c
    full.set$fecb2 <- full.set$fecb2c
    full.set$fecb3 <- full.set$fecb3c
  }
  
  full.set <- full.set[,-c(which(names(full.set) == "sizea1c"), which(names(full.set) == "sizea2c"), 
                           which(names(full.set) == "sizea3c"), which(names(full.set) == "sizeb1c"), 
                           which(names(full.set) == "sizeb2c"), which(names(full.set) == "sizeb3c"), 
                           which(names(full.set) == "sizec1c"), which(names(full.set) == "sizec2c"), 
                           which(names(full.set) == "sizec3c"), which(names(full.set) == "repstra1c"),
                           which(names(full.set) == "repstra2c"), which(names(full.set) == "repstra3c"), 
                           which(names(full.set) == "repstrb1c"), which(names(full.set) == "repstrb2c"), 
                           which(names(full.set) == "repstrb3c"), which(names(full.set) == "feca1c"), 
                           which(names(full.set) == "feca2c"), which(names(full.set) == "feca3c"), 
                           which(names(full.set) == "fecb1c"), which(names(full.set) == "fecb2c"), 
                           which(names(full.set) == "fecb3c"))]
  
  full.set <- full.set[,-c(which(names(full.set) == "alivegiven1"))]
  full.set <- full.set[,-c(which(names(full.set) == "alivegiven2"))]
  full.set <- full.set[,-c(which(names(full.set) == "alivegiven3"))]
  full.set <- full.set[,-c(which(names(full.set) == "deadgiven1"))]
  full.set <- full.set[,-c(which(names(full.set) == "deadgiven2"))]
  full.set <- full.set[,-c(which(names(full.set) == "deadgiven3"))]
  full.set <- full.set[,-c(which(names(full.set) == "obsgiven1"))]
  full.set <- full.set[,-c(which(names(full.set) == "obsgiven2"))]
  full.set <- full.set[,-c(which(names(full.set) == "obsgiven3"))]
  full.set <- full.set[,-c(which(names(full.set) == "nonobsgiven1"))]
  full.set <- full.set[,-c(which(names(full.set) == "nonobsgiven2"))]
  full.set <- full.set[,-c(which(names(full.set) == "nonobsgiven3"))]
  full.set <- full.set[,-c(which(names(full.set) == "juvgiven1"))]
  full.set <- full.set[,-c(which(names(full.set) == "juvgiven2"))]
  full.set <- full.set[,-c(which(names(full.set) == "juvgiven3"))]
  
  full.set <- full.set[,-c(which(names(full.set) =="year1"))]
  full.set <- full.set[,-c(which(names(full.set) =="year3"))]
  
  if (!is.na(spacing)) {
    full.set$density <- .density3(full.set, which(names(full.set) == "xpos2"), 
                                  which(names(full.set) == "ypos2"), 
                                  which(names(full.set) == "year2"), spacing)
  }
  
  if (reduce) {
    if (all(is.na(full.set$xpos1)) | length(unique(full.set$xpos1)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) == "xpos1"))]
    }
    if (all(is.na(full.set$ypos1)) | length(unique(full.set$ypos1)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) == "ypos1"))]
    }
    
    if (all(is.na(full.set$xpos2)) | length(unique(full.set$xpos2)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) == "xpos2"))]
    }
    if (all(is.na(full.set$ypos2)) | length(unique(full.set$ypos2)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) == "ypos2"))]
    }
    if (all(is.na(full.set$xpos3)) | length(unique(full.set$xpos3)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) == "xpos3"))]
    }
    if (all(is.na(full.set$ypos3)) | length(unique(full.set$ypos3)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) == "ypos3"))]
    }
    
    if (!censor) {
      full.set <- full.set[,-c(which(names(full.set) =="censor1"))]
      full.set <- full.set[,-c(which(names(full.set) =="censor2"))]
      full.set <- full.set[,-c(which(names(full.set) =="censor3"))]
    }
    
    if (all(is.na(full.set$sizea1)) | length(unique(full.set$sizea1)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) =="sizea1"))]
    }
    if (all(is.na(full.set$sizea2)) | length(unique(full.set$sizea2)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) =="sizea2"))]
    }
    if (all(is.na(full.set$sizea3)) | length(unique(full.set$sizea3)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) =="sizea3"))]
    }
    
    if (all(is.na(full.set$sizeb1)) | length(unique(full.set$sizeb1)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) =="sizeb1"))]
    }
    if (all(is.na(full.set$sizeb2)) | length(unique(full.set$sizeb2)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) =="sizeb2"))]
    }
    if (all(is.na(full.set$sizeb3)) | length(unique(full.set$sizeb3)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) =="sizeb3"))]
    }
    
    if (all(is.na(full.set$sizec1)) | length(unique(full.set$sizec1)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) =="sizec1"))]
    }
    if (all(is.na(full.set$sizec2)) | length(unique(full.set$sizec2)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) =="sizec2"))]
    }
    if (all(is.na(full.set$sizec3)) | length(unique(full.set$sizec3)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) =="sizec3"))]
    }
    
    if (isTRUE(all.equal(full.set$size1added, full.set$sizea1))) {
      full.set <- full.set[,-c(which(names(full.set) == "size1added"))]
    } else if (all(is.na(full.set$size1added)) | length(unique(full.set$size1added)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) == "size1added"))]
    }
    if (isTRUE(all.equal(full.set$size2added, full.set$sizea2))) {
      full.set <- full.set[,-c(which(names(full.set) == "size2added"))]
    } else if (all(is.na(full.set$size2added)) | length(unique(full.set$size2added)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) == "size2added"))]
    }
    if (isTRUE(all.equal(full.set$size3added, full.set$sizea3))) {
      full.set <- full.set[,-c(which(names(full.set) == "size3added"))]
    } else if (all(is.na(full.set$size3added)) | length(unique(full.set$size3added)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) == "size3added"))]
    }
    
    if (all(is.na(full.set$repstra1)) | length(unique(full.set$repstra1)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) == "repstra1"))]
    }
    if (all(is.na(full.set$repstra2)) | length(unique(full.set$repstra2)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) == "repstra2"))]
    }
    if (all(is.na(full.set$repstra3)) | length(unique(full.set$repstra3)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) == "repstra3"))]
    }
    
    if (all(is.na(full.set$repstrb1)) | length(unique(full.set$repstrb1)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) == "repstrb1"))]
    }
    if (all(is.na(full.set$repstrb2)) | length(unique(full.set$repstrb2)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) == "repstrb2"))]
    }
    if (all(is.na(full.set$repstrb3)) | length(unique(full.set$repstrb3)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) == "repstrb3"))]
    }
    
    if (isTRUE(all.equal(full.set$repstr1added, full.set$repstr1a))) {
      full.set <- full.set[,-c(which(names(full.set) =="repstr1added"))]
    } else if (all(is.na(full.set$repstr1added)) | length(unique(full.set$repstr1added)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) =="repstr1added"))]
    }
    if (isTRUE(all.equal(full.set$repstr2added, full.set$repstr2a))) {
      full.set <- full.set[,-c(which(names(full.set) =="repstr2added"))]
    } else if (all(is.na(full.set$repstr2added)) | length(unique(full.set$repstr2added)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) =="repstr2added"))]
    }
    if (isTRUE(all.equal(full.set$repstr3added, full.set$repstr3a))) {
      full.set <- full.set[,-c(which(names(full.set) =="repstr3added"))]
    } else if (all(is.na(full.set$repstr3added)) | length(unique(full.set$repstr3added)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) =="repstr3added"))]
    }
    
    if (all(is.na(full.set$feca1)) | length(unique(full.set$feca1)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) =="feca1"))]
    }
    if (all(is.na(full.set$feca2)) | length(unique(full.set$feca2)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) =="feca2"))]
    }
    if (all(is.na(full.set$feca3)) | length(unique(full.set$feca3)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) =="feca3"))]
    }
    
    if (all(is.na(full.set$fecb1)) | length(unique(full.set$fecb1)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) =="fecb1"))]
    }
    if (all(is.na(full.set$fecb2)) | length(unique(full.set$fecb2)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) =="fecb2"))]
    }
    if (all(is.na(full.set$fecb3)) | length(unique(full.set$fecb3)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) =="fecb3"))]
    }
    
    if (isTRUE(all.equal(full.set$fec1added, full.set$feca1))) {
      full.set <- full.set[,-c(which(names(full.set) =="fec1added"))]
    } else if (all(is.na(full.set$fec1added)) | length(unique(full.set$fec1added)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) =="fec1added"))]
    }
    if (isTRUE(all.equal(full.set$fec2added, full.set$feca2))) {
      full.set <- full.set[,-c(which(names(full.set) =="fec2added"))]
    } else if (all(is.na(full.set$fec2added)) | length(unique(full.set$fec2added)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) =="fec2added"))]
    }
    if (isTRUE(all.equal(full.set$fec3added, full.set$feca3))) {
      full.set <- full.set[,-c(which(names(full.set) =="fec3added"))]
    } else if (all(is.na(full.set$fec3added)) | length(unique(full.set$fec3added)) == 1) {
      full.set <- full.set[,-c(which(names(full.set) =="fec3added"))]
    }
    
    if (all(full.set$obsstatus1 == full.set$obsstatus1[1])) {
      full.set <- full.set[,-c(which(names(full.set) =="obsstatus1"))]
    }
    if (all(full.set$obsstatus2 == full.set$obsstatus2[1])) {
      full.set <- full.set[,-c(which(names(full.set) =="obsstatus2"))]
    }
    if (all(full.set$obsstatus3 == full.set$obsstatus3[1])) {
      full.set <- full.set[,-c(which(names(full.set) =="obsstatus3"))]
    }
    
    if (all(full.set$repstatus1 == full.set$repstatus1[1])) {
      full.set <- full.set[,-c(which(names(full.set) =="repstatus1"))]
    }
    if (all(full.set$repstatus2 == full.set$repstatus2[1])) {
      full.set <- full.set[,-c(which(names(full.set) =="repstatus2"))]
    }
    if (all(full.set$repstatus3 == full.set$repstatus3[1])) {
      full.set <- full.set[,-c(which(names(full.set) =="repstatus3"))]
    }
    
    if (all(full.set$fecstatus1 == full.set$fecstatus1[1])) {
      full.set <- full.set[,-c(which(names(full.set) =="fecstatus1"))]
    }
    if (all(full.set$fecstatus2 == full.set$fecstatus2[1])) {
      full.set <- full.set[,-c(which(names(full.set) =="fecstatus2"))]
    }
    if (all(full.set$fecstatus3 == full.set$fecstatus3[1])) {
      full.set <- full.set[,-c(which(names(full.set) =="fecstatus3"))]
    }
  }
  
  class(full.set) <- append(class(full.set), "hfvdata")
  
  return(full.set)
}

#' @export
.core_gather_ahvtohv <- function(subdata, popidcol, patchidcol, year2col, year3col, xcol = 0, ycol = 0,
                                 sizea2col = 0, sizea3col = 0, sizeb2col = 0, sizeb3col = 0, 
                                 sizec2col = 0, sizec3col = 0, repstra2col = 0, repstra3col = 0, 
                                 repstrb2col = 0, repstrb3col = 0, feca2col = 0, feca3col = 0, 
                                 fecb2col = 0, fecb3col = 0, alive2col = 0, alive3col = 0, 
                                 dead2col = 0, dead3col = 0, obs2col = 0, obs3col = 0, nonobs2col = 0,
                                 nonobs3col = 0, repstrrel = 0, fecrel = 0, stage2col = 0,
                                 stage3col = 0, juv2col = 0, juv3col = 0, censorcol = 0) {
  
  year.vec <- sort(na.omit(unique(union(subdata[,year2col], subdata[,year3col]))))
  
  size.vecorator <- apply(as.matrix(year.vec), 1, function(X) {
    if (!is.element(X, subdata[, year2col])) { return(0)}
    base.size <- 0
    if (sizea2col != 0) {
      if (!is.na(subdata[(which(subdata[,year2col] == X)), sizea2col])) {
        base.size <- base.size + subdata[(which(subdata[,year2col] == X)), sizea2col]
      }
    }
    if (sizeb2col != 0) {
      if (!is.na(subdata[(which(subdata[,year2col] == X)), sizeb2col])) {
        base.size <- base.size + subdata[(which(subdata[,year2col] == X)), sizeb2col]
      }
    }
    if (sizec2col != 0) {
      if (!is.na(subdata[(which(subdata[,year2col] == X)), sizec2col])) {
        base.size <- base.size + subdata[(which(subdata[,year2col] == X)), sizec2col]
      }
    }
    
    return(base.size)
  })
  
  if (all(size.vecorator == 0)) {
    return()
  }
  
  if (min(year.vec) < min(year.vec[which(size.vecorator > 0)])) {
    year.vec <- year.vec[-(which(year.vec < min(year.vec[which(size.vecorator > 0)])))]
  }
  
  firstyr <- min(year.vec)
  lastyr <- max(year.vec)
  
  subdata <- subdata[order(subdata$year2), ]
  
  pop <- NA
  patch <- NA
  if(popidcol > 0) {pop <- subdata[1,popidcol]}
  if(patchidcol > 0) {patch <- subdata[1,patchidcol]}
  
  if (xcol == 0) {
    xcol <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (ycol == 0) {
    ycol <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (sizea2col == 0) {
    sizea2col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (sizea3col == 0) {
    sizea3col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (sizeb2col == 0) {
    sizeb2col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (sizeb3col == 0) {
    sizeb3col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (sizec2col == 0) {
    sizec2col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (sizec3col == 0) {
    sizec3col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (repstra2col == 0) {
    repstra2col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (repstra3col == 0) {
    repstra3col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (repstrb2col == 0) {
    repstrb2col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (repstrb3col == 0) {
    repstrb3col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (feca2col == 0) {
    feca2col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (feca3col == 0) {
    feca3col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (fecb2col == 0) {
    fecb2col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (fecb3col == 0) {
    fecb3col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (alive2col == 0) {
    alive2col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (alive3col == 0) {
    alive3col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (dead2col == 0) {
    dead2col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (dead3col == 0) {
    dead3col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (obs2col == 0) {
    obs2col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (obs3col == 0) {
    obs3col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (nonobs2col == 0) {
    nonobs2col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (nonobs3col == 0) {
    nonobs3col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (stage2col == 0) {
    stage2col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (stage3col == 0) {
    stage3col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (juv2col == 0) {
    juv2col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (juv3col == 0) {
    juv3col <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, NA)
  }
  if (censorcol == 0) {
    censorcol <- (dim(subdata)[2] + 1)
    subdata <- cbind.data.frame(subdata, 1)
  }
  
  row.col <- which(names(subdata) == "rowid")
  
  for (i in c(firstyr:lastyr)) {
    
    counter <- i - (firstyr - 1)
    
    if (counter == 2) {start.vec <- row.vec}
    
    if (i == firstyr) {
      year.t1 <- NA; xcol.t1 <- NA; ycol.t1 <- NA; sizea.t1 <- NA; sizeb.t1 <- NA; 
      sizec.t1 <- NA; repstra.t1 <- NA; repstrb.t1 <- NA; feca.t1 <- NA; fecb.t1 <- NA; 
      alivegvn.t1 <- NA; deadgvn.t1 <- NA; obsgvn.t1 <- NA; nonobsgvn.t1 <- NA; 
      stage.t1 <- NA; juvgvn.t1 <- NA; censor1 <- NA; 
      
      year.t2 <- i; rowid <- subdata[counter, row.col]; 
      xcol.t2 <- subdata[counter, xcol]; ycol.t2 <- subdata[counter, ycol]; 
      sizea.t2 <- subdata[counter, sizea2col]; sizeb.t2 <- subdata[counter, sizeb2col];
      sizec.t2 <- subdata[counter, sizec2col]; repstra.t2 <- subdata[counter, repstra2col];
      repstrb.t2 <- subdata[counter, repstrb2col]; feca.t2 <- subdata[counter, feca2col];
      fecb.t2 <- subdata[counter, fecb2col]; alivegvn.t2 <- subdata[counter, alive2col]; 
      deadgvn.t2 <- subdata[counter, dead2col]; obsgvn.t2 <- subdata[counter, obs2col]; 
      nonobsgvn.t2 <- subdata[counter, nonobs2col]; stage.t2 <- subdata[counter, stage2col]; 
      juvgvn.t2 <- subdata[counter, juv2col]; censor2 <- subdata[counter, censorcol]
    }
    
    if (i > firstyr) {
      year.t1 <- i - 1; xcol.t1 <- subdata[(counter - 1), xcol]; ycol.t1 <- subdata[(counter - 1), ycol]; 
      sizea.t1 <- subdata[(counter - 1), sizea2col]; sizeb.t1 <- subdata[(counter - 1), sizeb2col]; 
      sizec.t1 <- subdata[(counter - 1), sizec2col]; repstra.t1 <- subdata[(counter - 1), repstra2col]; 
      repstrb.t1 <- subdata[(counter - 1), repstrb2col]; feca.t1 <- subdata[(counter - 1), feca2col]; 
      fecb.t1 <- subdata[(counter - 1), fecb2col]; alivegvn.t1 <- subdata[(counter - 1), alive2col]; 
      deadgvn.t1 <- subdata[(counter - 1), dead2col]; obsgvn.t1 <- subdata[(counter - 1), obs2col];
      nonobsgvn.t1 <- subdata[(counter - 1), nonobs2col]; stage.t1 <- subdata[(counter - 1), stage2col]; 
      juvgvn.t1 <- subdata[(counter - 1), juv2col]; censor1 <- subdata[(counter - 1), censorcol];
      
      year.t2 <- i; rowid <- subdata[(counter - 1), row.col]; 
      xcol.t2 <- subdata[(counter - 1), xcol]; ycol.t2 <- subdata[(counter - 1), ycol]; 
      sizea.t2 <- subdata[(counter - 1), sizea3col]; sizeb.t2 <- subdata[(counter - 1), sizeb3col]; 
      sizec.t2 <- subdata[(counter - 1), sizec3col]; repstra.t2 <- subdata[(counter - 1), repstra3col]; 
      repstrb.t2 <- subdata[(counter - 1), repstrb3col]; feca.t2 <- subdata[(counter - 1), feca3col]; 
      fecb.t2 <- subdata[(counter - 1), fecb3col]; alivegvn.t2 <- subdata[(counter - 1), alive3col]; 
      deadgvn.t2 <- subdata[(counter - 1), dead3col]; obsgvn.t2 <- subdata[(counter - 1), obs3col];
      nonobsgvn.t2 <- subdata[(counter - 1), nonobs3col]; stage.t2 <- subdata[(counter - 1), stage3col]; 
      juvgvn.t2 <- subdata[(counter - 1), juv3col]; censor2 <- subdata[(counter - 1), censorcol];
    }
    
    if (i == lastyr) {
      year.t3 <- NA; xcol.t3 <- NA; ycol.t3 <- NA; sizea.t3 <- NA; sizeb.t3 <- NA; 
      sizec.t3 <- NA; repstra.t3 <- NA; repstrb.t3 <- NA; feca.t3 <- NA; fecb.t3 <- NA; 
      alivegvn.t3 <- NA; deadgvn.t3 <- NA; obsgvn.t3 <- NA; nonobsgvn.t3 <- NA; 
      stage.t3 <- NA; juvgvn.t3 <- NA; censor3 <- NA
    }
    if (i < lastyr) {
      year.t3 <- i + 1; xcol.t3 <- subdata[counter, xcol]; ycol.t3 <- subdata[counter, ycol]; 
      sizea.t3 <- subdata[counter, sizea3col]; sizeb.t3 <- subdata[counter, sizeb3col]; 
      sizec.t3 <- subdata[counter, sizec3col]; repstra.t3 <- subdata[counter, repstra3col]; 
      repstrb.t3 <- subdata[counter, repstrb3col]; feca.t3 <- subdata[counter, feca3col]; 
      fecb.t3 <- subdata[counter, fecb3col]; alivegvn.t3 <- subdata[counter, alive3col]; 
      deadgvn.t3 <- subdata[counter, dead3col]; obsgvn.t3 <- subdata[counter, obs3col]; 
      nonobsgvn.t3 <- subdata[counter, nonobs3col]; stage.t3 <- subdata[counter, stage3col]; 
      juvgvn.t3 <- subdata[counter, juv3col]; censor3 <- subdata[counter, censorcol]
    }
    
    row.vec <- cbind.data.frame(pop, patch, rowid, year.t1, xcol.t1, ycol.t1, sizea.t1, sizeb.t1, 
                                sizec.t1, repstra.t1, repstrb.t1, feca.t1, fecb.t1, 
                                alivegvn.t1, deadgvn.t1, obsgvn.t1, nonobsgvn.t1, stage.t1, 
                                juvgvn.t1, censor1, 
                                
                                year.t2, xcol.t2, ycol.t2, sizea.t2, sizeb.t2, sizec.t2, 
                                repstra.t2, repstrb.t2, feca.t2, fecb.t2, alivegvn.t2, 
                                deadgvn.t2, obsgvn.t2, nonobsgvn.t2, stage.t2, 
                                juvgvn.t2, censor2, 
                                
                                year.t3, xcol.t3, ycol.t3, sizea.t3, sizeb.t3, sizec.t3, 
                                repstra.t3, repstrb.t3, feca.t3, fecb.t3, alivegvn.t3,
                                deadgvn.t3, obsgvn.t3, nonobsgvn.t3, stage.t3, 
                                juvgvn.t3, censor3, firstyr, lastyr)
    
    if (counter == 2) {indiv.data <- rbind.data.frame(start.vec, row.vec)}
    if (counter > 2) {indiv.data <- rbind.data.frame(indiv.data, row.vec)}
  }
  
  names(indiv.data) <- c("popid", "patchid", "rowid", "year1", "xcol1", "ycol1", "sizea1", "sizeb1", 
                         "sizec1", "repstra1", "repstrb1", "feca1", "fecb1", "alivegvn1",
                         "deadgvn1", "obsgvn1", "nonobsgvn1", "stage1", "juvgvn1", "censor1",
                         
                         "year2", "xcol2", "ycol2", "sizea2", "sizeb2", "sizec2", "repstra2", 
                         "repstrb2", "feca2", "fecb2", "alivegvn2", "deadgvn2", "obsgvn2",
                         "nonobsgvn2", "stage2", "juvgvn2", "censor2", 
                         
                         "year3", "xcol3", "ycol3", "sizea3", "sizeb3", "sizec3", "repstra3", 
                         "repstrb3", "feca3", "fecb3", "alivegvn3", "deadgvn3", "obsgvn3",
                         "nonobsgvn3", "stage3", "juvgvn3", "censor3", "firstseen", "lastseen")
  
  indiv.data <- indiv.data[(which(!is.na(indiv.data$year3))),]
  
  indiv.data$alive1 <- 0
  indiv.data$alive2 <- 0
  indiv.data$alive3 <- 0
  
  indiv.data$alive1[unique(intersect(which(indiv.data$year1 >= indiv.data$firstseen), which(indiv.data$year1 <= indiv.data$lastseen)))] <- 1
  indiv.data$alive2[unique(intersect(which(indiv.data$year2 >= indiv.data$firstseen), which(indiv.data$year2 <= indiv.data$lastseen)))] <- 1
  indiv.data$alive3[unique(intersect(which(indiv.data$year3 >= indiv.data$firstseen), which(indiv.data$year3 <= indiv.data$lastseen)))] <- 1
  
  if (alive2col != 0) {
    indiv.data$alive1[which(indiv.data$alivegvn1 == 1)] <- 1
    indiv.data$alive2[which(indiv.data$alivegvn2 == 1)] <- 1
  }
  if (alive3col != 0) {
    indiv.data$alive2[which(indiv.data$alivegvn2 == 1)] <- 1
    indiv.data$alive3[which(indiv.data$alivegvn3 == 1)] <- 1
  }
  if (dead2col != 0) {
    indiv.data$alive1[which(indiv.data$deadgvn1 == 1)] <- 0
    indiv.data$alive2[which(indiv.data$deadgvn2 == 1)] <- 0
  }
  if (dead3col != 0) {
    indiv.data$alive2[which(indiv.data$deadgvn2 == 1)] <- 0
    indiv.data$alive3[which(indiv.data$deadgvn3 == 1)] <- 0
  }
  
  indiv.data$obsstatus1 <- 0
  indiv.data$obsstatus2 <- 0
  indiv.data$obsstatus3 <- 0
  
  indiv.data$obsstatus1[which(indiv.data$sizea1 > 0)] <- 1
  indiv.data$obsstatus2[which(indiv.data$sizea2 > 0)] <- 1
  indiv.data$obsstatus3[which(indiv.data$sizea3 > 0)] <- 1
  
  indiv.data$obsstatus1[which(indiv.data$sizeb1 > 0)] <- 1
  indiv.data$obsstatus2[which(indiv.data$sizeb2 > 0)] <- 1
  indiv.data$obsstatus3[which(indiv.data$sizeb3 > 0)] <- 1
  
  indiv.data$obsstatus1[which(indiv.data$sizec1 > 0)] <- 1
  indiv.data$obsstatus2[which(indiv.data$sizec2 > 0)] <- 1
  indiv.data$obsstatus3[which(indiv.data$sizec3 > 0)] <- 1
  
  indiv.data$obsstatus1[which(indiv.data$repstra1 > 0)] <- 1
  indiv.data$obsstatus2[which(indiv.data$repstra2 > 0)] <- 1
  indiv.data$obsstatus3[which(indiv.data$repstra3 > 0)] <- 1
  
  indiv.data$obsstatus1[which(indiv.data$repstrb1 > 0)] <- 1
  indiv.data$obsstatus2[which(indiv.data$repstrb2 > 0)] <- 1
  indiv.data$obsstatus3[which(indiv.data$repstrb3 > 0)] <- 1
  
  if (obs2col != 0) {
    indiv.data$obsstatus1[which(indiv.data$obsgvn1 == 1)] <- 1
    indiv.data$obsstatus2[which(indiv.data$obsgvn2 == 1)] <- 1
  }
  if (obs3col != 0) {
    indiv.data$obsstatus2[which(indiv.data$obsgvn2 == 1)] <- 1
    indiv.data$obsstatus3[which(indiv.data$obsgvn3 == 1)] <- 1
  }
  if (nonobs2col != 0) {
    indiv.data$obsstatus1[which(indiv.data$nonobsgvn1 == 1)] <- 0
    indiv.data$obsstatus2[which(indiv.data$nonobsgvn2 == 1)] <- 0
  }
  if (nonobs3col != 0) {
    indiv.data$obsstatus2[which(indiv.data$nonobsgvn2 == 1)] <- 0
    indiv.data$obsstatus3[which(indiv.data$nonobsgvn3 == 1)] <- 0
  }
  
  indiv.data$matstatus1 <- 1
  indiv.data$matstatus2 <- 1
  indiv.data$matstatus3 <- 1
  
  if (juv2col != 0) {
    indiv.data$matstatus1[which(indiv.data$juvgvn1 == 1)] <- 0
    indiv.data$matstatus2[which(indiv.data$juvgvn2 == 1)] <- 0
  }
  if (juv3col != 0) {
    indiv.data$matstatus2[which(indiv.data$juvgvn2 == 1)] <- 0
    indiv.data$matstatus3[which(indiv.data$juvgvn3 == 1)] <- 0
  }
  
  indiv.data$repstatus1 <- 0
  indiv.data$repstatus2 <- 0
  indiv.data$repstatus3 <- 0
  
  indiv.data$repstatus1[which(indiv.data$repstra1 > 0)] <- 1
  indiv.data$repstatus2[which(indiv.data$repstra2 > 0)] <- 1
  indiv.data$repstatus3[which(indiv.data$repstra3 > 0)] <- 1
  
  indiv.data$repstatus1[which(indiv.data$repstrb1 > 0)] <- 1
  indiv.data$repstatus2[which(indiv.data$repstrb2 > 0)] <- 1
  indiv.data$repstatus3[which(indiv.data$repstrb3 > 0)] <- 1
  
  indiv.data$fecstatus1 <- 0
  indiv.data$fecstatus2 <- 0
  indiv.data$fecstatus3 <- 0
  
  indiv.data$fecstatus1[which(indiv.data$feca1 > 0)] <- 1
  indiv.data$fecstatus2[which(indiv.data$feca2 > 0)] <- 1
  indiv.data$fecstatus3[which(indiv.data$feca3 > 0)] <- 1
  
  indiv.data$fecstatus1[which(indiv.data$fecb1 > 0)] <- 1
  indiv.data$fecstatus2[which(indiv.data$fecb2 > 0)] <- 1
  indiv.data$fecstatus3[which(indiv.data$fecb3 > 0)] <- 1
  
  return(indiv.data)
}

#' Estimate Density on Basis of Cartesian Coordinates
#' 
#' \code{.density3()} estimates density on the basis of Cartesian coordinates and
#' spacing information supplied as input. It is used internally by 
#' \code{\link{historicalize3}()} and \code{\link{verticalize3}()}.
#' 
#' @param data Demographic dataset in historical vertical format.
#' @param xcol Number of column in \code{data} corresponding to x position.
#' @param ycol Number of column in \code{data} corresponding to y position.
#' @param yearcol Number of column in \code{data} corresponding to time step.
#' @param spacing Resolution of density estimation, as a scalar numeric.
#' 
#' @return This function returns the original data frame supplied as \code{data} but
#' with a new variable added to the end of the data frame showing the local density
#' of the individual.
#' 
#' @keywords internal
#' @noRd
.density3 <- function(data, xcol, ycol, yearcol, spacing) {
  data$Xgen <- floor(data[,xcol] / spacing)
  data$Ygen <- floor(data[,ycol] / spacing)
  
  data$grid.id <- paste(data$Xgen, data$Ygen, sep = ",")
  
  grid.densities <- xtabs(paste("~ grid.id + ", yearcol), data = data)
  
  data$density.est <- apply(as.matrix(c(1:dim(data)[1])), 1, function(X) {
    grid.densities[as.character(data$grid.id[X]), as.character(data[X, yearcol])]
  })
  
  data$density.est[which(is.na(data$Xgen))] <- NA
  
  return(data)
}
