#' Create a Stageframe for Population Matrix Projection Analysis
#' 
#' \code{sf_create()} returns a data frame describing each ahistorical life
#' history stage. This data frame can be used as input into matrix-estimating 
#' functions such as \code{\link{flefko3}()}, where it instructs the function 
#' as to the treatment of each stage during matrix creation.
#'
#' @param sizes A numeric vector of the typical or representative size of each
#' life history stage.
#' @param stagenames An optional vector of stage names, in the same order as
#' elements in sizes. If an IPM or function-based matrix with many stages is
#' desired, then two stages that occur within the dataset and represent the 
#' lower and upper size limits of the IPM must be marked as \code{ipm} in this 
#' vector. Note that these stages must be mature stages, and should have all 
#' characteristics besides size equal. If two or more groups of stages, each
#' with its own characteristics, are to be developed for IPM or function-based
#' matrix development, then an even number of stages with two stages marking
#' the minimum and maximum size of each group should be marked, with all other
#' characteristics equal within each group.
#' @param repstatus A vector denoting the binomial reproductive status of each
#' life history stage. Defaults to 1.
#' @param obsstatus A vector denoting the binomial observation status of each
#' life history stage. Defaults to 1, but may be changed for unobservable 
#' stages.
#' @param propstatus A vector denoting whether each life history stage is a 
#' propagule. Such stages are generally only used in fecundity estimation. 
#' Defaults to NA.
#' @param immstatus A vector denoting whether each stage is immature. Must be
#' composed of binomial values if given. Defaults to NA.
#' @param matstatus A vector denoting whether each stage is mature. Must be
#' composed of binomial values if given. Defaults to 1 for all stages defined 
#' in \code{sizes}.
#' @param minage An optional vector denoting the minimum age at which a stage
#' can begin. Only used in age x stage matrix development. Defaults to NA.
#' @param maxage An optional vector denoting the maximum age at which a stage
#' should end. Only used in age x stage matrix development. Defaults to NA.
#' @param indataset A vector designating which stages are found within the 
#' dataset. While \code{\link{rlefko2}()} and \code{\link{rlefko3}()} can use all stages in
#' the input dataset, \code{\link{flefko3}()} and \code{\link{flefko2}()} can only handle
#' size-classified stages with non-overlapping combinations of size and 
#' reproductive status, plus one immature stage. In the latter scenario, other 
#' stages should be excluded by marking them as 0 in this vector.
#' @param binhalfwidth A numeric vector giving the half-width of size bins. Used
#' and required to classify individuals appropriately within size classes.
#' Defaults to 0.5 for all sizes.
#' @param ipmbins If an IPM is desired, then this parameter sets the number of
#' stages to create for that IPM. This number is in addition to any stages
#' that are not size-classified. Defaults to 100, and numbers greater than this
#' yield a warning about the loss of statistical power and increasing chance of
#' matrix over-parameterization resulting from increasing numbers of stages.
#' @param roundsize This parameter sets the precision of size classification,
#' and equals the number of digits used in rounding sizes. Defaults to 5.
#'
#' @return A stageframe object, which is a data frame that includes information
#' on the stage name, size, reproductive status, observation status, propagule 
#' status, immaturity status, maturity status, presence within the core dataset, 
#' counts of similarly sized stages, raw bin half-width, and the minimum, 
#' center, and maximum of each size bin, as well as its width. If minimum and
#' maximum ages were specified, then these are also included. Also includes an 
#' empty string variable that can be used to describe stages meaningfully. This
#' object can be used as the \code{stageframe} input for \code{\link{flefko3}()} 
#' \code{\link{flefko2}()}, \code{\link{rlefko3}()}, and \code{\link{rlefko2}()}.
#' 
#' @return Variables in this data frame include the following:
#' \item{stagenames}{The unique names of the stages to be analyzed.}
#' \item{size}{The typical or representative size at which the stage occurs.}
#' \item{repstatus}{A binomial variable designating whether the stage occurs
#' is reproductive.}
#' \item{obsstatus}{A binomial variable designating whether the stage occurs
#' is observable.}
#' \item{propstatus}{A binomial variable designating whether the stage occurs
#' as a propagule.}
#' \item{immstatus}{A binomial variable designating whether the stage occurs
#' in immaturity.}
#' \item{matstatus}{A binomial variable designating whether the stage occurs
#' in maturity.}
#' \item{indataset}{A binomial variable describing whether the stage occurs
#' in the input dataset.}
#' \item{binhalfwidth_raw}{The half-width input for the size bin corresponding
#' to the stage.}
#' \item{sizebin_min}{The minimum size at which the stage may occur.}
#' \item{sizebin_max}{The maximum size at which the stage may occur.}
#' \item{sizebin_center}{The centroid of the size bin at which the stage may
#' occur.}
#' \item{sizebin_width}{The width of the size bin corresponding to the stage.}
#' \item{min_age}{The minimum age at which the stage may occur.}
#' \item{max_age}{The maximum age at which the stage may occur.}
#' \item{comments}{A text field to hold verbal stage descriptions.}
#'  
#' @examples
#' sizevector <- c(0, 0, 0, 0, 0, 0, 1, 2.5, 4.5, 8, 17.5)
#' stagevector <- c("SD", "P1", "P2", "P3", "SL", "D", "XSm", "Sm", "Md", "Lg", 
#'                  "XLg")
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
#' cypframe_raw$comments[(cypframe_raw$stagenames == "SD")] <- "Dormant seed"
#' cypframe_raw$comments[(cypframe_raw$stagenames == "P1")] <- "1st yr protocorm"
#' cypframe_raw$comments[(cypframe_raw$stagenames == "P2")] <- "2nd yr protocorm"
#' cypframe_raw$comments[(cypframe_raw$stagenames == "P3")] <- "3rd yr protocorm"
#' cypframe_raw$comments[(cypframe_raw$stagenames == "SL")] <- "Seedling"
#' cypframe_raw$comments[(cypframe_raw$stagenames == "D")] <- "Dormant adult"
#' cypframe_raw$comments[(cypframe_raw$stagenames == "XSm")] <- "Extra small adult (1 shoot)"
#' cypframe_raw$comments[(cypframe_raw$stagenames == "Sm")] <- "Small adult (2-3 shoots)"
#' cypframe_raw$comments[(cypframe_raw$stagenames == "Md")] <- "Medium adult (4-5 shoots)"
#' cypframe_raw$comments[(cypframe_raw$stagenames == "Lg")] <- "Large adult (6-10 shoots)"
#' cypframe_raw$comments[(cypframe_raw$stagenames == "XLg")] <- "Extra large adult (>10 shoots)"
#' 
#' cypframe_raw
#' @export
sf_create <- function(sizes, stagenames = NA, repstatus = 1, obsstatus = 1, propstatus = NA, 
                      immstatus = NA, matstatus = 1, minage = NA, maxage = NA, indataset = NA, 
                      binhalfwidth = 0.5, ipmbins = 100, roundsize = 5) {
  
  #Initially we standardize the length of option vectors
  matsize <- length(sizes)
  
  if (is.na(stagenames[1]) & length(stagenames) == 1) {
    stagenames <- seq(1, matsize)
  }
  
  if (repstatus[1] == 1 & length(repstatus) == 1) {
    repstatus <- rep(1, matsize)
  } else if (repstatus[1] == 0 & length(repstatus) == 1) {
    repstatus <- rep(0, matsize)
  } else if (length(repstatus) != matsize) {
    stop("Input vectors are not equal in length.")
  }
  
  if (obsstatus[1] == 1 & length(obsstatus) == 1) {
    obsstatus <- rep(1, matsize)
  } else if (obsstatus[1] == 0 & length(obsstatus) == 1) {
    obsstatus <- rep(0, matsize)
  } else if (length(obsstatus) != matsize) {
    stop("Input vectors are not equal in length.")
  }
  
  if (matstatus[1] == 1 & length(matstatus) == 1) {
    matstatus <- rep(1, matsize)
  } else if (matstatus[1] == 0 & length(matstatus) == 1) {
    matstatus <- rep(0, matsize)
  } else if (length(matstatus) != matsize) {
    stop("Input vectors are not equal in length.")
  }
  
  if (is.na(immstatus[1]) & length(immstatus) == 1) {
    immstatus <- rep(0, matsize)
  }
  
  if (is.na(propstatus[1]) & length(propstatus) == 1) {
    propstatus <- rep(0, matsize)
  }
  
  if (is.na(minage[1]) & length(minage) == 1) {
    no_age <- TRUE
  } else {
    if (length(minage) != matsize | length(maxage) != matsize) {
      stop("Input vectors are not equal in length.")
    }
    
    no_age <- FALSE
  }
  
  if (is.na(indataset[1]) & length(indataset) == 1) {
    indataset <- seq(1, matsize)
  }
  
  if (is.numeric(binhalfwidth[1]) & length(binhalfwidth) == 1) {
    truebinvec <- rep(binhalfwidth[1], matsize)
    binhalfwidth <- truebinvec
  }
  
  #Here we check for illegal inputs and stop where necessary
  if (length(stagenames) != matsize) {
    stop("Stagename option must be either equal NA or 'ipm', or be a vector of the same length as the sizes vector.", call. = FALSE)
  }
  
  if (length(repstatus) != matsize) {
    stop("Repstatus option must either equal 0 or 1, or be a vector of 1's and 0's of the same length as the sizes vector.", call. = FALSE)
  }
  
  if (length(obsstatus) != matsize) {
    stop("Obsstatus option must either equal 0 or 1, or be a vector of 1's and 0's of the same length as the sizes vector.", call. = FALSE)
  }
  
  if (length(matstatus) != matsize) {
    stop("Matstatus option must either equal 0 or 1, or be a vector of 1's and 0's of the same length as the sizes vector.", call. = FALSE)
  }
  
  if (length(immstatus) != matsize) {
    stop("Immstatus option must either equal NA or be a vector of the same length as the sizes vector.", call. = FALSE)
  }
  
  if (length(propstatus) != matsize) {
    stop("Propstatus option must either equal NA or be a vector of the same length as the sizes vector.", call. = FALSE)
  }
  
  if (length(indataset) != matsize) {
    stop("Indataset option must either equal NA or be a vector of 1's and 0's of the same length as the sizes vector.", call. = FALSE)
  }
  
  if (length(binhalfwidth) != matsize) {
    stop("Binhalfwidth option must either equal a single number or be a numeric vector of the same length as the sizes vector.", call. = FALSE)
  }
  
  if (any(!is.numeric(sizes))) {
    warning("Size should be numeric if matrices will be created with flefko2() or flefko3().")
  }
  
  if (ipmbins != as.integer(ipmbins) | ipmbins < 2) {
    stop("Please enter a valid integer greater than 1 for ipmbins option.", .call = FALSE)
  } else if (ipmbins > 100) {
    warning("High ipmbin numbers may lead to dramatic decreases in statistical power and overparameterized matrices.")
  }
  
  if (no_age == TRUE) {
    sfmat <- cbind.data.frame(stagenames = as.character(stagenames), size = sizes, repstatus = repstatus, 
                              obsstatus = obsstatus, propstatus = propstatus, immstatus = immstatus, 
                              matstatus = matstatus, indataset = indataset, binhalfwidth_raw = binhalfwidth, 
                              min_age = NA, max_age = NA, stringsAsFactors = FALSE)
  } else {
    sfmat <- cbind.data.frame(stagenames = as.character(stagenames), size = sizes, repstatus = repstatus, 
                              obsstatus = obsstatus, propstatus = propstatus, immstatus = immstatus, 
                              matstatus = matstatus, indataset = indataset, binhalfwidth_raw = binhalfwidth,
                              min_age = minage, max_age = maxage, stringsAsFactors = FALSE)
  }
  
  if (any(tolower(stagenames) == "ipm")) {
    if ((length(which(tolower(stagenames) == "ipm")) %% 2) != 0) {
      stop("Pairs of stages must be marked as 'ipm' in stagenames, corresponding to the upper and lower size limits for ipm size class and bin calculation for each pair.", .call = FALSE)
    }
    
    if (any(matstatus[which(tolower(stagenames) == "ipm")] != 1)) {
      stop("Size classes used as bases for function-based stages must be mature classes. Please reassign all such class to mature status.", .call = FALSE)
    }
    
    if ((length(unique(sizes[(which(tolower(stagenames) == "ipm"))])) %% 2) != 0) {
      stop("Stages marked 'ipm' must differ in size within pairs.", .call = FALSE)
    }
    
    if (any(indataset[(which(tolower(stagenames) == "ipm"))] == 0)) {
      stop("All stages used to develop function-based stages must exist within the vertical dataset.", .call = FALSE)
    }
    
    ipmbase <- which(tolower(stagenames) == "ipm")
    ipmmarked <- sfmat[ipmbase,]
    ipmmarked$ipmindex <- ipmmarked$propstatus + (ipmmarked$immstatus * 10) + (ipmmarked$matstatus * 100) + 
      (ipmmarked$obsstatus * 1000) + (ipmmarked$repstatus * 10000)
    notin <- which(tolower(stagenames) != "ipm")
    
    ipmgrouplist <- split(ipmmarked, ipmmarked$ipmindex)
    
    divisions <- length(ipmbase) / 2
    
    if (divisions == 1) {
      currentbins <- ipmbins
      
    } else {
      spreadbase <- sapply(ipmgrouplist, function(X) {
        X$size[which(X$size == max(X$size))] - X$size[which(X$size == min(X$size))]
      })
      
      spread <- spreadbase / sum(spreadbase)
      currentbins <- round(spread * ipmbins)
    }
    
    addedon_list <- apply(as.matrix(c(1:length(currentbins))), 1, function(X) {
      .ipmerator(ipmgrouplist[[X]], sfmat, no_age, notin, currentbins[X], roundsize)
    })
    
    sfmat <- sfmat[-ipmbase,]
    
    addedon <- do.call("rbind.data.frame", addedon_list)
    sfmat <- rbind.data.frame(sfmat, addedon)
    
  }
  
  matsize <- length(sfmat$size)
  
  sfmat$sizebin_min <- apply(as.matrix(c(1:matsize)), 1, function(X) {
    return(round((sfmat$size[X] - sfmat$binhalfwidth[X]), digits = roundsize))
  })
  
  sfmat$sizebin_max <- apply(as.matrix(c(1:matsize)), 1, function(X) {
    return(round((sfmat$size[X] + sfmat$binhalfwidth[X]), digits = roundsize))
  })
  
  sfmat$sizebin_center <- apply(as.matrix(c(1:matsize)), 1, function(X) {
    return(round(sfmat$size[X], digits = roundsize))
  })
  
  sfmat$sizebin_width <- sfmat$sizebin_max - sfmat$sizebin_min
  
  sfmat$comments <- "Add text stage description if desired"
  
  class(sfmat) <- append(class(sfmat), "stageframe")
  
  return(sfmat)
}

#' Rewrite Stageframe To Reflect IPM Stages
#' 
#' \code{.ipmerator()} searches through the supplied stageframe and rearranges the
#' information if an IPM is desired. This allows the original input in the 
#' stageframe to be developed easily in shorthand, which this function then parses
#' into a long format for analysis.
#' 
#' @param ipmdata Stageframe rows corresponding to the IPM sections only.
#' @param maindata The full stageframe.
#' @param notin Stageframe rows corresonponding to non-IPM sections only.
#' @param currentbins Number of bins to use in IPM stage development.
#' @param roundsize Resolution at which to round size within size bins.
#' 
#' @return A rearranged stageframe lengthened to include all IPM stages.
#' 
#' @keywords internal
#' @noRd
.ipmerator <- function(ipmdata, maindata, no_age, notin, currentbins, roundsize) {
  
  isx <- c(which(ipmdata$size == min(ipmdata$size)), which(ipmdata$size == max(ipmdata$size)))
  loipmstage <- isx[1]
  hiipmstage <- isx[2]
  
  if (ipmdata$repstatus[loipmstage] != ipmdata$repstatus[hiipmstage] | ipmdata$obsstatus[loipmstage] != ipmdata$obsstatus[hiipmstage]) {
    stop("All input characteristics of stages used for IPM stage development must be equal.", .call = FALSE)
  }
  
  if (ipmdata$matstatus[loipmstage] != ipmdata$matstatus[hiipmstage] | ipmdata$immstatus[loipmstage] != ipmdata$immstatus[hiipmstage]) {
    stop("All input characteristics of stages used for IPM stage development must be equal.", .call = FALSE)
  }
  
  if (ipmdata$propstatus[loipmstage] != ipmdata$propstatus[hiipmstage]) {
    stop("All input characteristics of stages used for IPM stage development must be equal.", .call = FALSE)
  }
  
  if (no_age == FALSE) {
    lominagetest <- ipmdata$min_age[loipmstage]
    lomaxagetest <- ipmdata$max_age[loipmstage]
    himinagetest <- ipmdata$min_age[hiipmstage]
    himaxagetest <- ipmdata$max_age[hiipmstage]
    
    if (is.na(lominagetest)) {lominagetest <- 0}
    if (is.na(lomaxagetest)) {lomaxagetest <- 1001}
    if (is.na(himinagetest)) {himinagetest <- 0}
    if (is.na(himaxagetest)) {himaxagetest <- 1001}
    
    if (lominagetest != himinagetest | lomaxagetest != himaxagetest) {
      stop("All input characteristics of stages used for IPM stage development must be equal.", .call = FALSE)
    }
  }
  
  if (length(which(ipmdata$matstatus[notin] == 1)) > 0) {
    testsizes <- maindata$size[which(maindata$matstatus[notin] == 1)]
    defurc <- seq(from = ipmdata$size[loipmstage], to = ipmdata$size[hiipmstage], length.out = currentbins)
    ipmhalftest <- (defurc[2] - defurc[1]) / 2
    
    if ((ipmdata$size[loipmstage] - ipmhalftest) < max(testsizes)) {
      ipmdata$size[loipmstage] <- max(testsizes) + ipmhalftest
    }
  }
  
  extrasizes <- seq(from = ipmdata$size[loipmstage], to = ipmdata$size[hiipmstage], length.out = currentbins)
  
  extrastagenames <- apply(as.matrix(extrasizes), 1, function(X) {
    paste0("sz", round(X, roundsize), " rp", ipmdata$repstatus[loipmstage], " mt", 
           ipmdata$matstatus[loipmstage], " ob", ipmdata$obsstatus[loipmstage])
  })
  ipmbinhalfwidth <- (extrasizes[2] - extrasizes[1]) / 2
  
  newstuff <- cbind.data.frame(stagenames = extrastagenames, size = extrasizes, 
                               repstatus = rep(ipmdata$repstatus[loipmstage], length.out = currentbins),
                               obsstatus = rep(ipmdata$obsstatus[loipmstage], length.out = currentbins), 
                               propstatus = rep(ipmdata$propstatus[loipmstage], length.out = currentbins), 
                               immstatus = rep(ipmdata$immstatus[loipmstage], length.out = currentbins), 
                               matstatus = rep(ipmdata$matstatus[loipmstage], length.out = currentbins), 
                               indataset = rep(ipmdata$indataset[loipmstage], length.out = currentbins), 
                               binhalfwidth_raw = rep(ipmbinhalfwidth, length.out = currentbins))
  
  if (!no_age) {
    newstuff <- cbind.data.frame(newstuff, min_age = rep(ipmdata$min_age[loipmstage], length.out = currentbins),
                                 max_age = rep(ipmdata$max_age[loipmstage], length.out = currentbins))
  } else {
    newstuff <- cbind.data.frame(newstuff, min_age = NA, max_age = NA)
  }
  
  return(newstuff)
}

#' Standardize Stageframe Input For Analysis
#' 
#' \code{.sfreassess()} takes a stageframe as input, and used information supplied 
#' there and through the reproduction and overwrite matrices to rearrange this
#' into a format usable by the matrix creation functions, \code{\link{flefko3}()}, \code{\link{flefko2}()},
#' \code{\link{rlefko3}()}, and \code{\link{rlefko2}()}.
#' 
#' @param stageframe The original stageframe.
#' @param repmatrix The original reproduction matrix.
#' @param overwrite The original overwrite table, as supplied by the \code{\link{overwrite}()}
#' function.
#' 
#' @return This functon returns a list with a modified stageframe usable in
#' projection matrix construction, and an associated reproduction matrix.
#' 
#' @keywords internal
#' @noRd
.sf_reassess <- function(stageframe, repmatrix, overwrite) {
  if (any(grepl("stage", names(stageframe)))) {
    stage.column <- min(which(grepl("stage", names(stageframe))))
    stage.vec <- stageframe[,stage.column]
    stage.vec <- as.character(stage.vec)
  } else {stage.vec <- c(1:dim(stageframe)[1])}
  
  if (any(grepl("size", names(stageframe)))) {
    orig.size.column <- min(which(grepl("size", names(stageframe))))
    size.min <- which(grepl("bin_min", names(stageframe)))
    size.max <- which(grepl("bin_max", names(stageframe)))
    size.width <- which(grepl("bin_width", names(stageframe)))
    
    if (any(grepl("bin_ctr", names(stageframe)))) {
      size.ctr <- which(grepl("bin_ctr", names(stageframe)))
    } else if (any(grepl("bin_cen", names(stageframe)))) {
      size.ctr <- which(grepl("bin_cen", names(stageframe)))
    }
    
  } else {
    orig.size.column <- 1
    size.min <- 1
    size.ctr <- 1
    size.max <- 1
    size.width <- 1
    
    warning("Not certain which column in the input stageframe provides size information so using the first.");
  }
  orig.size.vec <- stageframe[,orig.size.column]
  size.min.vec <- stageframe[,size.min]
  size.ctr.vec <- stageframe[,size.ctr]
  size.max.vec <- stageframe[,size.max]
  size.width.vec <- stageframe[,size.width]
  
  if (any(grepl("rep", names(stageframe)))) {
    rep.column <- min(which(grepl("rep", names(stageframe)))) 
    rep.vec <- stageframe[,rep.column]
  }
  
  if (any(grepl("spr", names(stageframe)))) {
    obs.column <- min(which(grepl("spr", names(stageframe))))
    obs.vec <- stageframe[,obs.column]
  } else if (any(grepl("obs", names(stageframe)))) {
    obs.column <- min(which(grepl("obs", names(stageframe))))
    obs.vec <- stageframe[,obs.column]
  } else {
    obs.vec <- rep(1, length(orig.size.vec))
  }
  
  if (any(grepl("prop", names(stageframe)))) {
    prop.column <- min(which(grepl("prop", names(stageframe))))
    prop.vec <- stageframe[,prop.column]
  } else {
    prop.vec <- rep(0, length(orig.size.vec))
  }
  
  if (any(grepl("imm", names(stageframe)))) {
    imm.column <- min(which(grepl("imm", names(stageframe))))
    imm.vec <- stageframe[,imm.column]
  } else {
    imm.vec <- c(1, rep(0, (length(orig.size.vec) - 1)))
  }
  
  if (any(grepl("mat", names(stageframe)))) {
    mat.column <- min(which(grepl("mat", names(stageframe))))
    mat.vec <- stageframe[,mat.column]
  } else {
    mat.vec <- c(0, rep(1, (length(orig.size.vec) - 1)))
  }
  
  if (any(grepl("indata", names(stageframe)))) {
    ind.column <- min(which(grepl("indata", names(stageframe))))
    ind.vec <- stageframe[,ind.column]
  } else {
    ind.vec <- c(1, rep(1, (length(orig.size.vec) - 1)))
  }
  
  if (any(grepl("binhalf", names(stageframe)))) {
    bin.column <- min(which(grepl("binhalf", names(stageframe))))
    bin.vec <- stageframe[,bin.column]
  } else {
    bin.vec <- rep(1, length(orig.size.vec))
  }
  
  if (any(grepl("min_ag", names(stageframe)))) {
    minage.column <- min(which(grepl("min_ag", names(stageframe))))
    minage.vec <- stageframe[,minage.column]
  } else if (any(grepl("minag", names(stageframe)))) {
    minage.column <- min(which(grepl("minag", names(stageframe))))
    minage.vec <- stageframe[,minage.column]
  } else {
    minage.vec <- rep(NA, length(orig.size.vec))
  }
  
  if (any(grepl("max_ag", names(stageframe)))) {
    maxage.column <- min(which(grepl("max_ag", names(stageframe))))
    maxage.vec <- stageframe[,maxage.column]
  } else if (any(grepl("maxag", names(stageframe)))) {
    maxage.column <- min(which(grepl("maxag", names(stageframe))))
    maxage.vec <- stageframe[,maxage.column]
  } else {
    maxage.vec <- rep(NA, length(orig.size.vec))
  }
  
  if (any(grepl("comment", names(stageframe)))) {
    com.column <- min(which(grepl("comment", names(stageframe))))
    com.vec <- stageframe[,com.column]
  } else {
    com.vec <- rep("no description", length(orig.size.vec))
  }
  
  if (any(!is.na(minage.vec)) | any(!is.na(maxage.vec))) {
    age <- TRUE
  } else {
    age <- FALSE
  }
  
  if (is.matrix(repmatrix)) {
    rep.entry.stages <- apply(repmatrix, 1, sum)
    rep.entry.stages[which(rep.entry.stages > 0)] <- 1
    rep.col <- apply(repmatrix, 2, sum)
    rep.col[which(rep.col > 0)] <- 1
    
    stageframe$rep.yn <- rep.col
    rep.vec <- rep.col
    
  } else if (!is.matrix(repmatrix)) {
    if (sum(prop.vec + imm.vec) > 1) {
      rep.entry.stages <- prop.vec + imm.vec; rep.entry.stages[which(rep.entry.stages > 0)] <- 1
    } else {
      rep.entry.stages <- c(1, rep(0, (length(orig.size.vec) - 1)))
      warning("No information on reproductive entry stages provided. Assuming the first stage is the entry stage into the life cycle.")
    }
    
    if (!exists("rep.vec")) {
      rep.vec <- mat.vec
      warning("No information on reproductive stages given. Assuming all mature stages are reproductive.")
    }
    repmatrix <- matrix(0, ncol = length(mat.vec), nrow = length(mat.vec))
    repmatrix[which(rep.entry.stages > 0),which(mat.vec == 1)] <- 1
  }
  
  if (sum(prop.vec) > 0 & sum(imm.vec) > 0) {
    orig.stage.vec.r <- c(stage.vec[which(prop.vec == 1)], stage.vec[which(imm.vec == 1 & prop.vec == 0)], 
                          stage.vec[which(mat.vec == 1 & imm.vec == 0)])
    orig.size.vec.r <- c(orig.size.vec[which(prop.vec == 1)], orig.size.vec[which(imm.vec == 1 & prop.vec == 0)], 
                         orig.size.vec[which(mat.vec == 1 & imm.vec == 0)])
    rep.vec.r <- c(rep.vec[which(prop.vec == 1)], rep.vec[which(imm.vec == 1 & prop.vec == 0)], 
                   rep.vec[which(mat.vec == 1 & imm.vec == 0)])
    obs.vec.r <- c(obs.vec[which(prop.vec == 1)], obs.vec[which(imm.vec == 1 & prop.vec == 0)], 
                   obs.vec[which(mat.vec == 1 & imm.vec == 0)])
    prop.vec.r <- c(prop.vec[which(prop.vec == 1)], prop.vec[which(imm.vec == 1 & prop.vec == 0)], 
                    prop.vec[which(mat.vec == 1 & imm.vec == 0)])
    imm.vec.r <- c(imm.vec[which(prop.vec == 1)], imm.vec[which(imm.vec == 1 & prop.vec == 0)], 
                   imm.vec[which(mat.vec == 1 & imm.vec == 0)])
    mat.vec.r <- c(mat.vec[which(prop.vec == 1)], mat.vec[which(imm.vec == 1 & prop.vec == 0)], 
                   mat.vec[which(mat.vec == 1 & imm.vec == 0)])
    ind.vec.r <- c(ind.vec[which(prop.vec == 1)], ind.vec[which(imm.vec == 1 & prop.vec == 0)], 
                   ind.vec[which(mat.vec == 1 & imm.vec == 0)])
    bin.vec.r <- c(bin.vec[which(prop.vec == 1)], bin.vec[which(imm.vec == 1 & prop.vec == 0)], 
                   bin.vec[which(mat.vec == 1 & imm.vec == 0)])
    minage.vec.r <- c(minage.vec[which(prop.vec == 1)], minage.vec[which(imm.vec == 1 & prop.vec == 0)], 
                      minage.vec[which(mat.vec == 1 & imm.vec == 0)])
    maxage.vec.r <- c(maxage.vec[which(prop.vec == 1)], maxage.vec[which(imm.vec == 1 & prop.vec == 0)], 
                      maxage.vec[which(mat.vec == 1 & imm.vec == 0)])
    size.min.vec.r <- c(size.min.vec[which(prop.vec == 1)], size.min.vec[which(imm.vec == 1 & prop.vec == 0)], 
                        size.min.vec[which(mat.vec == 1 & imm.vec == 0)])
    size.max.vec.r <- c(size.max.vec[which(prop.vec == 1)], size.max.vec[which(imm.vec == 1 & prop.vec == 0)], 
                        size.max.vec[which(mat.vec == 1 & imm.vec == 0)])
    size.ctr.vec.r <- c(size.ctr.vec[which(prop.vec == 1)], size.ctr.vec[which(imm.vec == 1 & prop.vec == 0)], 
                        size.ctr.vec[which(mat.vec == 1 & imm.vec == 0)])
    size.width.vec.r <- c(size.width.vec[which(prop.vec == 1)], size.width.vec[which(imm.vec == 1 & prop.vec == 0)], 
                          size.width.vec[which(mat.vec == 1 & imm.vec == 0)])
    com.vec.r <- c(com.vec[which(prop.vec == 1)], com.vec[which(imm.vec == 1 & prop.vec == 0)], 
                   com.vec[which(mat.vec == 1 & imm.vec == 0)])
    
    repmatrix.r.1 <- repmatrix[c(which(prop.vec == 1), which(imm.vec == 1 & prop.vec == 0), 
                                 which(mat.vec == 1 & imm.vec == 0)),]
    repmatrix.r <- repmatrix.r.1[,c(which(prop.vec == 1), which(imm.vec == 1 & prop.vec == 0), 
                                    which(mat.vec == 1 & imm.vec == 0))]
    
  } else if (sum(prop.vec) == 0 & sum(imm.vec) > 0) {
    orig.stage.vec.r <- c(stage.vec[which(imm.vec == 1)], stage.vec[which(mat.vec == 1 & imm.vec == 0)])
    orig.size.vec.r <- c(orig.size.vec[which(imm.vec == 1)], orig.size.vec[which(mat.vec == 1 & imm.vec == 0)])
    rep.vec.r <- c(rep.vec[which(imm.vec == 1)], rep.vec[which(mat.vec == 1 & imm.vec == 0)])
    obs.vec.r <- c(obs.vec[which(imm.vec == 1)], obs.vec[which(mat.vec == 1 & imm.vec == 0)])
    prop.vec.r <- c(prop.vec[which(imm.vec == 1)], prop.vec[which(mat.vec == 1 & imm.vec == 0)])
    imm.vec.r <- c(imm.vec[which(imm.vec == 1)], imm.vec[which(mat.vec == 1 & imm.vec == 0)])
    mat.vec.r <- c(mat.vec[which(imm.vec == 1)], mat.vec[which(mat.vec == 1 & imm.vec == 0)])
    ind.vec.r <- c(ind.vec[which(imm.vec == 1)], ind.vec[which(mat.vec == 1 & imm.vec == 0)])
    bin.vec.r <- c(bin.vec[which(imm.vec == 1)], bin.vec[which(mat.vec == 1 & imm.vec == 0)])
    minage.vec.r <- c(minage.vec[which(imm.vec == 1)], minage.vec[which(mat.vec == 1 & imm.vec == 0)])
    maxage.vec.r <- c(maxage.vec[which(imm.vec == 1)], maxage.vec[which(mat.vec == 1 & imm.vec == 0)])
    size.min.vec.r <- c(size.min.vec[which(imm.vec == 1)], size.min.vec[which(mat.vec == 1 & imm.vec == 0)])
    size.max.vec.r <- c(size.max.vec[which(imm.vec == 1)], size.max.vec[which(mat.vec == 1 & imm.vec == 0)])
    size.ctr.vec.r <- c(size.ctr.vec[which(imm.vec == 1)], size.ctr.vec[which(mat.vec == 1 & imm.vec == 0)])
    size.width.vec.r <- c(size.width.vec[which(imm.vec == 1)], size.width.vec[which(mat.vec == 1 & imm.vec == 0)])
    com.vec.r <- c(com.vec[which(imm.vec == 1)], com.vec[which(mat.vec == 1 & imm.vec == 0)])
    
    repmatrix.r.1 <- repmatrix[c(which(imm.vec == 1), which(mat.vec == 1 & imm.vec == 0)),]
    repmatrix.r <- repmatrix.r.1[,c(which(imm.vec == 1), which(mat.vec == 1 & imm.vec == 0))]
    
  } else if (sum(prop.vec > 0) & sum(imm.vec == 0)) {
    orig.stage.vec.r <- c(stage.vec[which(prop.vec == 1)], stage.vec[which(mat.vec == 1 & prop.vec == 0)])
    orig.size.vec.r <- c(orig.size.vec[which(prop.vec == 1)], orig.size.vec[which(mat.vec == 1 & prop.vec == 0)])
    rep.vec.r <- c(rep.vec[which(prop.vec == 1)], rep.vec[which(mat.vec == 1 & prop.vec == 0)])
    obs.vec.r <- c(obs.vec[which(prop.vec == 1)], obs.vec[which(mat.vec == 1 & prop.vec == 0)])
    prop.vec.r <- c(prop.vec[which(prop.vec == 1)], prop.vec[which(mat.vec == 1 & prop.vec == 0)])
    imm.vec.r <- c(imm.vec[which(prop.vec == 1)], imm.vec[which(mat.vec == 1 & prop.vec == 0)])
    mat.vec.r <- c(mat.vec[which(prop.vec == 1)], mat.vec[which(mat.vec == 1 & prop.vec == 0)])
    ind.vec.r <- c(ind.vec[which(prop.vec == 1)], ind.vec[which(mat.vec == 1 & prop.vec == 0)])
    bin.vec.r <- c(bin.vec[which(prop.vec == 1)], bin.vec[which(mat.vec == 1 & prop.vec == 0)])
    minage.vec.r <- c(minage.vec[which(prop.vec == 1)], minage.vec[which(mat.vec == 1 & prop.vec == 0)])
    maxage.vec.r <- c(maxage.vec[which(prop.vec == 1)], maxage.vec[which(mat.vec == 1 & prop.vec == 0)])
    size.min.vec.r <- c(size.min.vec[which(prop.vec == 1)], size.min.vec[which(mat.vec == 1 & prop.vec == 0)])
    size.max.vec.r <- c(size.max.vec[which(prop.vec == 1)], size.max.vec[which(mat.vec == 1 & prop.vec == 0)])
    size.ctr.vec.r <- c(size.ctr.vec[which(prop.vec == 1)], size.ctr.vec[which(mat.vec == 1 & prop.vec == 0)])
    size.width.vec.r <- c(size.width.vec[which(prop.vec == 1)], size.width.vec[which(mat.vec == 1 & prop.vec == 0)])
    com.vec.r <- c(com.vec[which(prop.vec == 1)], com.vec[which(mat.vec == 1 & prop.vec == 0)])
    
    repmatrix.r.1 <- repmatrix[c(which(prop.vec == 1), which(mat.vec == 1 & prop.vec == 0)),]
    repmatrix.r <- repmatrix.r.1[,c(which(prop.vec == 1), which(mat.vec == 1 & prop.vec == 0))]
    
  } else if (sum(prop.vec == 0) & sum(imm.vec == 0)) {
    orig.stage.vec.r <- stage.vec
    orig.size.vec.r <- orig.size.vec
    rep.vec.r <- rep.vec
    obs.vec.r <- obs.vec
    prop.vec.r <- prop.vec
    imm.vec.r <- imm.vec
    mat.vec.r <- mat.vec
    ind.vec.r <- ind.vec
    bin.vec.r <- bin.vec
    minage.vec.r <- minage.vec
    maxage.vec.r <- maxage.vec
    size.min.vec.r <- size.min.vec
    size.max.vec.r <- size.max.vec
    size.ctr.vec.r <- size.ctr.vec
    size.width.vec.r <- size.width.vec
    com.vec.r <- com.vec
    
    repmatrix.r <- repmatrix
  }
  
  stageframe.reassessed <- cbind.data.frame(new_stage_id = as.numeric(c(1:length(orig.stage.vec.r))), 
                                            orig_stage_id = as.character(orig.stage.vec.r), 
                                            original_size = as.numeric(orig.size.vec.r), 
                                            bin_size_ctr = as.numeric(size.ctr.vec.r), 
                                            bin_size_min = as.numeric(size.min.vec.r), 
                                            bin_size_max = as.numeric(size.max.vec.r), 
                                            repstatus = as.numeric(rep.vec.r), 
                                            obsstatus = as.numeric(obs.vec.r), 
                                            propstatus = as.numeric(prop.vec.r), 
                                            immstatus = as.numeric(imm.vec.r), 
                                            matstatus = as.numeric(mat.vec.r), 
                                            indataset = as.numeric(ind.vec.r), 
                                            bin_size_width = as.numeric(size.width.vec.r), 
                                            bin_raw_halfwidth = as.numeric(bin.vec.r))
  stageframe.reassessed$alive <- 1
  
  stageframe.reassessed <- rbind.data.frame(stageframe.reassessed, c((length(orig.stage.vec.r) + 1), "Dead", 0, 
                                                                     0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0))
  
  if (age) {
    stageframe.reassessed <- cbind.data.frame(stageframe.reassessed, min_age = c(as.numeric(minage.vec.r), 0),
                                              max_age = c(as.numeric(maxage.vec.r), 0))
  } else {
    stageframe.reassessed <- cbind.data.frame(stageframe.reassessed, min_age = NA,
                                              max_age = NA)
  }
  
  com.vec.r <- c(com.vec.r, "Dead")
  
  stageframe.reassessed$comments <- com.vec.r
  
  if(any(duplicated(stageframe.reassessed$orig_stage_id) == TRUE)) {
    stop("All stage names provided in stageframe must be unique.")
  }
  
  if(!all(is.na(overwrite))) {
    if (length(setdiff(unique(na.omit(overwrite$stage3)), stageframe.reassessed$orig_stage_id)) > 0) {
      stop("Stage names in overwrite table (stage3) must match stages in stageframe.")
    }
    
    if (length(setdiff(unique(na.omit(overwrite$stage2)), stageframe.reassessed$orig_stage_id)) > 0) {
      stop("Stage names in overwrite table (stage2) must match stages in stageframe.")
    }
    
    if (length(setdiff(unique(na.omit(overwrite$stage1)), c(stageframe.reassessed$orig_stage_id, "rep"))) > 0) {
      stop("Stage names in overwrite table (stage1) must match stages in stageframe.")
    }
    
    if (length(setdiff(unique(na.omit(overwrite$eststage3)), stageframe.reassessed$orig_stage_id)) > 0) {
      stop("Stage names in overwrite table (eststage3) must match stages in stageframe.")
    }
    
    if (length(setdiff(unique(na.omit(overwrite$eststage2)), stageframe.reassessed$orig_stage_id)) > 0) {
      stop("Stage names in overwrite table (eststage2) must match stages in stageframe.")
    }
    
    if (length(setdiff(unique(na.omit(overwrite$eststage1)), stageframe.reassessed$orig_stage_id)) > 0) {
      stop("Stage names in overwrite table (eststage1) must match stages in stageframe.")
    }
  }
  
  if (!is.element("stageframe", class(stageframe.reassessed))) {
    class(stageframe.reassessed) <- append(class(stageframe.reassessed), "stageframe")
  }
  
  return(list(stageframe.reassessed, repmatrix.r))
}

#' Create an Overwrite Table for Population Matrix Development
#' 
#' \code{overwrite()} returns a data frame describing which particular transitions
#' within an ahistorical or historical projection matrix to overwrite with either
#' given rates and probabilities, or other estimated transitions.
#'
#' @param stage3 The name of the stage in time \emph{t}+1 in the transition to be 
#' replaced.
#' @param stage2 The name of the stage in time \emph{t} in the transition to be 
#' replaced.
#' @param stage1 The name of the stage in time \emph{t}-1 in the transition to be
#' replaced. Only needed if a historical matrix is to be produced. Use \code{rep}
#' if all reproductive stages are to be used, and leave empty or use \code{all}
#' if all stages in stageframe are to be used.
#' @param eststage3 The name of the stage to replace \code{stage3}. Only needed if
#' a transition will be replaced by another estimated transition.
#' @param eststage2 The name of the stage to replace \code{stage2}. Only needed if
#' a transition will be replaced by another estimated transition.
#' @param eststage1 The name of the stage to replace \code{stage1}. Only needed if
#' a transition will be replaced by another estimated transition, and the matrix
#' to be estimated is a historical matrix.
#' @param givenrate A fixed rate or probability to replace for the transition
#' described by \code{stage3}, \code{stage2}, and \code{stage1}.
#' @param type A vector denoting the kind of transition that will be replaced.
#' Should be entered as 1, S, or s for survival, or 2, F, or f for fecundity.
#' Defaults to 1, for survival transition.
#'
#' @return A data frame that puts the above vectors together and can be used as
#' input in \code{\link{flefko3}}, \code{\link{flefko2}}, \code{\link{rlefko3}},
#' and \code{\link{rlefko2}}.
#' 
#' Variables in this data frame include the following:
#' \item{stage3}{Stage at time \emph{t}+1 in the transition to be replaced.}
#' \item{stage2}{Stage at time \emph{t} in the transition to be replaced.}
#' \item{stage1}{Stage at time \emph{t}-1 in the transition to be replaced.}
#' \item{eststage3}{Stage at time \emph{t}+1 in the transition to replace the
#' transition designated by \code{stage3}, \code{stage2}, and \code{stage1}.}
#' \item{eststage2}{Stage at time \emph{t} in the transition to replace the
#' transition designated by \code{stage3}, \code{stage2}, and \code{stage1}.}
#' \item{eststage1}{Stage at time \emph{t}-1 in the transition to replace the
#' transition designated by \code{stage3}, \code{stage2}, and \code{stage1}.}
#' \item{givenrate}{A constant to be used as the value of the transition.}
#' \item{convtype}{Designates whether the transition is a survival-transition
#' probability (1) or a fecundity rate (2).}
#' 
#' @examples
#' cypover2r <- overwrite(stage3 = c("SD", "P1", "P2", "P3", "SL", "SL", "D", 
#'                        "XSm", "Sm"), stage2 = c("SD", "SD", "P1", "P2", "P3", 
#'                        "SL", "SL", "SL", "SL"), eststage3 = c(NA, NA, NA, NA, 
#'                        NA, NA, "D", "XSm", "Sm"), eststage2 = c(NA, NA, NA, NA, 
#'                        NA, NA, "XSm", "XSm", "XSm"), givenrate = c(0.1, 0.2, 
#'                        0.2, 0.2, 0.25, 0.4, NA, NA, NA), type = c("S", "S", "S",
#'                        "S", "S", "S", "S", "S", "S"))
#' 
#' cypover2r
#' 
#' cypover3r <- overwrite(stage3 = c("SD", "SD", "P1", "P1", "P2", "P3", "SL", 
#'                        "SL", "SL", "D", "XSm", "Sm", "D", "XSm", "Sm"), 
#'                        stage2 = c("SD", "SD", "SD", "SD", "P1", "P2", "P3", 
#'                        "SL", "SL", "SL", "SL", "SL", "SL", "SL", "SL"),
#'                        stage1 = c("SD", "rep", "SD", "rep", "SD", "P1", "P2", 
#'                        "P3", "SL", "P3", "P3", "P3", "SL", "SL", "SL"),
#'                        eststage3 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, "D", 
#'                        "XSm", "Sm", "D", "XSm", "Sm"), eststage2 = c(NA, NA, NA, 
#'                        NA, NA, NA, NA, NA, NA, "XSm", "XSm", "XSm", "XSm", "XSm",
#'                        "XSm"), eststage1 = c(NA, NA, NA, NA, NA, NA, NA, NA, NA,
#'                        "XSm", "XSm", "XSm", "XSm", "XSm", "XSm"), 
#'                        givenrate = c(0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.25, 0.4,
#'                        0.4, NA, NA, NA, NA, NA, NA), type = c("S", "S", "S", 
#'                        "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", "S", 
#'                        "S"))
#' 
#' cypover3r
#' 
#' @export
overwrite <- function(stage3, stage2, stage1 = NA, eststage3 = NA, eststage2 = NA, 
                       eststage1 = NA, givenrate = NA, type = NA) {
  if (length(stage3) != length(stage2)) {
    stop("All transitions to overwrite require information at least for stage2 and stage3. These inputs must also be of equal length.",
         .call = FALSE)
  }
  
  fulllength <- max(length(stage3), length(stage2), length(stage1), length(eststage3), 
                    length(eststage2), length(eststage1), length(givenrate), length(type))
  
  if (length(stage3) != fulllength) {
    stop("Please provide all input vectors in the same order.", .call = FALSE)
  }
  
  if (length(stage1) < fulllength) {
    missinglength <- fulllength - length(stage1)
    stage1 <- append(stage1, rep(NA, missinglength))
  }
  if (length(eststage3) < fulllength) {
    missinglength <- fulllength - length(eststage3)
    eststage3 <- append(eststage3, rep(NA, missinglength))
  }
  if (length(eststage2) < fulllength) {    
    missinglength <- fulllength - length(eststage2)
    eststage2 <- append(eststage2, rep(NA, missinglength))
  }
  if (length(eststage1) < fulllength) {
    missinglength <- fulllength - length(eststage1)
    eststage1 <- append(eststage1, rep(NA, missinglength))
  }
  if (length(givenrate) < fulllength) {
    missinglength <- fulllength - length(givenrate)
    givenrate <- append(givenrate, rep(NA, missinglength))
  }
  if (length(type) < fulllength) {
    missinglength <- fulllength - length(type)
    type <- append(type, rep(NA, missinglength))
  }
  
  if(!all(is.na(type))) {
    convtype <- rep(1, length(type))
    convtype[which(type == "F")] <- 2
    convtype[which(type == "f")] <- 2
  } else {
    convtype <- rep(1, length(stage3))
  }
  
  fullpack <- cbind.data.frame(stage3, stage2, stage1, eststage3, eststage2, eststage1, givenrate, convtype)

  return(fullpack)
}

#' Check and Reorganize Overwrite Table Into Usable Format
#' 
#' \code{.overwrite_reassess()} takes an overwrite table as supplied by the
#' \code{\link{overwrite}()} function, and checks and rearranges it to make it
#' usable by the main matrix creation functions, \code{\link{flefko3}()}, \code{\link{flefko2}()},
#' \code{\link{rlefko3}()}, and \code{\link{rlefko2}()}.
#' 
#' @param overwritetable The original overwrite table created by the \code{\link{overwrite}()}
#' function,
#' @param stageframe The correct stageframe, already modified by \code{\link{.sf_reassess}()}.
#' @param historical Logical value denoting whether matrix to create is historical.
#' 
#' @return A corrected overwrite table, usable in matrix creation.
#' 
#' @keywords internal
#' @noRd
.overwrite_reassess <- function(overwritetable, stageframe, historical) {
  
  if (!all(is.na(overwritetable))) {
    #First we make sure that the data is in the right format  
    overwritetable$stage3 <- as.character(overwritetable$stage3)
    overwritetable$stage2 <- as.character(overwritetable$stage2)
    overwritetable$stage1 <- as.character(overwritetable$stage1)
    overwritetable$eststage3 <- as.character(overwritetable$eststage3)
    overwritetable$eststage2 <- as.character(overwritetable$eststage2)
    overwritetable$eststage1 <- as.character(overwritetable$eststage1)
    
    #Now we will take some codes and replace them with actual stages, and check for some stop conditions
    if (historical == TRUE) {
      reassessed <- apply(as.matrix(c(1:dim(overwritetable)[1])), 1, function(X) {
        checkna2vec <- c(overwritetable[X, "stage3"], overwritetable[X, "stage2"])
        
        if (!all(!is.na(checkna2vec))) {
          stop("All entries for stage2 and stage3 in overwrite table must refer to possible life history stages and cannot include NAs.")
        }
        
        checknaestvec <- c(overwritetable[X, "eststage3"], overwritetable[X, "eststage2"], overwritetable[X, "eststage1"])
        
        if (all(is.na(checknaestvec)) | all(!is.na(checknaestvec))) {
          if (!is.na(overwritetable[X, "stage1"])) {
            if (is.element(overwritetable[X, "stage1"], stageframe$orig_stage_id)) {
              return(overwritetable[X,])
            } else if (overwritetable[X, "stage1"] == "rep") {
              shrubbery <- cbind.data.frame(stage3 = overwritetable[X, "stage3"], 
                                            stage2 = overwritetable[X, "stage2"], 
                                            stage1 = stageframe$orig_stage_id[which(stageframe$repstatus == 1)], 
                                            eststage3 = overwritetable[X, "eststage3"], 
                                            eststage2 = overwritetable[X, "eststage2"], 
                                            eststage1 = overwritetable[X, "eststage1"],
                                            givenrate = overwritetable[X, "givenrate"], 
                                            convtype = overwritetable[X, "convtype"])
              
              return(shrubbery)
            } else if (overwritetable[X, "stage1"] == "all") {
              shrubbery <- cbind.data.frame(stage3 = overwritetable[X, "stage3"], 
                                            stage2 = overwritetable[X, "stage2"], 
                                            stage1 = stageframe$orig_stage_id, 
                                            eststage3 = overwritetable[X, "eststage3"],
                                            eststage2 = overwritetable[X, "eststage2"], 
                                            eststage1 = overwritetable[X, "eststage1"],
                                            givenrate = overwritetable[X, "givenrate"], 
                                            convtype = overwritetable[X, "convtype"])
              
              return(shrubbery)
            } else {
              stop("Please use official stage names, NAs, or 'rep' in stage1.")
            }
            
          } else if (is.na(overwritetable[X, "stage1"])) {
            shrubbery <- cbind.data.frame(stage3 = overwritetable[X, "stage3"], 
                                          stage2 = overwritetable[X, "stage2"], 
                                          stage1 = stageframe$orig_stage_id, 
                                          eststage3 = overwritetable[X, "eststage3"],
                                          eststage2 = overwritetable[X, "eststage2"], 
                                          eststage1 = overwritetable[X, "eststage1"],
                                          givenrate = overwritetable[X, "givenrate"], 
                                          convtype = overwritetable[X, "convtype"])
            
            return(shrubbery)
          } else {
            return(overwritetable[X,])
          }
        } else {
          stop("If setting transitions equal to other estimated transitions, then eststage3, eststage2, and eststage1 must refer to possible life history stages and cannot include NAs.")
        }
      })
      
      reassessed <- do.call(rbind.data.frame, reassessed)
      
      checkcode <- apply(as.matrix(c(1:dim(reassessed)[1])), 1, function(X) {
        return(paste(reassessed$stage3[X], reassessed$stage2[X], reassessed$stage1[X]))
      })
      
    } else if (historical == FALSE) {
      reassessed <- overwritetable
      reassessed$stage1 <- NA
      
      if (length(reassessed$stage3) == 0) {
        stop("No overwrite transitions listed are for ahistorical matrix use. Overwrite table will not be used.")
      }
      
      checkcode <- apply(as.matrix(c(1:dim(reassessed)[1])), 1, function(X) {
        return(paste(reassessed$stage3[X], reassessed$stage2[X]))
      })
      
    } else {
      stop("Not clear if matrix approach is historical or ahistorical")
    }
    
    if (length(unique(checkcode)) < length(checkcode)) {
      dups <- checkcode[(which(duplicated(checkcode)))]
      
      duptest <- apply(as.matrix(c(1:length(dups))), 1, function(X) {
        gvrvec <- reassessed$givenrate[which(checkcode == dups[X])]
        
        if (length(unique(gvrvec)) == 1) {thisguy = TRUE} else {thisguy = FALSE}
        
        return(thisguy)
      })
      
      if (any(!duptest)) stop("Multiple entries with different values for the same stage transition are not allowed in the overwrite table. If performing an ahistorical analysis, then this may be due to different given rates of substitutions caused by dropping stage at time t-1, if historical transitions are used as input.")
    }
  } else {
    reassessed <- data.frame(stage3 = NA, stage2 = NA, stage1 = NA, eststage3 = NA, eststage2 = NA, 
                             eststage1 = NA, givenrate = NA, convtype = -1)
  }
  
  return(reassessed)
}

