## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----Ch9.0--------------------------------------------------------------------
rm(list=ls(all=TRUE))
library(lefko3)

## ----Ch9.1--------------------------------------------------------------------
sizevector <- c(1, 1, 2, 3) # These sizes are not from the original paper
stagevector <- c("Sdl", "Veg", "SmFlo", "LFlo")
repvector <- c(0, 0, 1, 1)
obsvector <- c(1, 1, 1, 1)
matvector <- c(0, 1, 1, 1)
immvector <- c(1, 0, 0, 0)
propvector <- c(1, 0, 0, 0)
indataset <- c(1, 1, 1, 1)
binvec <- c(0.5, 0.5, 0.5, 0.5)
comments <- c("Seedling", "Vegetative adult", "Small flowering",
  "Large flowering")

anthframe <- sf_create(sizes = sizevector, stagenames = stagevector,
  repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
  immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
  propstatus = propvector, comments = comments)
anthframe

## ----Ch9.2--------------------------------------------------------------------
# POPN C 2003-2004
XC3 <- matrix(c(0, 0, 1.74, 1.74,
  0.208333333, 0, 0, 0.057142857,
  0.041666667, 0.076923077, 0, 0,
  0.083333333, 0.076923077, 0.066666667, 0.028571429), 4, 4, byrow = TRUE)
XC3

## ----Ch9.3--------------------------------------------------------------------
# POPN C 2004-2005
XC4 <- matrix(c(0, 0, 0.3, 0.6,
0.32183908, 0.142857143, 0, 0,
0.16091954, 0.285714286, 0, 0,
0.252873563, 0.285714286, 0.5, 0.6), 4, 4, byrow = TRUE)

# POPN C 2005-2006
XC5 <- matrix(c(0, 0, 0.50625, 0.675,
0, 0, 0, 0.035714286,
0.1, 0.068965517, 0.0625, 0.107142857,
0.3, 0.137931034, 0, 0.071428571), 4, 4, byrow = TRUE)

# POPN E 2003-2004
XE3 <- matrix(c(0, 0, 2.44, 6.569230769,
0.196428571, 0, 0, 0,
0.125, 0.5, 0, 0,
0.160714286, 0.5, 0.133333333, 0.076923077), 4, 4, byrow = TRUE)

# POPN E 2004-2005
XE4 <- matrix(c(0, 0, 0.45, 0.646153846,
0.06557377, 0.090909091, 0.125, 0,
0.032786885, 0, 0.125, 0.076923077,
0.049180328, 0, 0.125, 0.230769231), 4, 4, byrow = TRUE)

# POPN E 2005-2006
XE5 <- matrix(c(0, 0, 2.85, 3.99,
0.083333333, 0, 0, 0,
0, 0, 0, 0,
0.416666667, 0.1, 0, 0.1), 4, 4, byrow = TRUE)

# POPN F 2003-2004
XF3 <- matrix(c(0, 0, 1.815, 7.058333333,
0.075949367, 0, 0.05, 0.083333333,
0.139240506, 0, 0, 0.25,
0.075949367, 0, 0, 0.083333333), 4, 4, byrow = TRUE)

# POPN F 2004-2005
XF4 <- matrix(c(0, 0, 1.233333333, 7.4,
0.223880597, 0, 0.111111111, 0.142857143,
0.134328358, 0.272727273, 0.166666667, 0.142857143,
0.119402985, 0.363636364, 0.055555556, 0.142857143), 4, 4, byrow = TRUE)

# POPN F 2005-2006
XF5 <- matrix(c(0, 0, 1.06, 3.372727273,
0.073170732, 0.025, 0.033333333, 0,
0.036585366, 0.15, 0.1, 0.136363636,
0.06097561, 0.225, 0.166666667, 0.272727273), 4, 4, byrow = TRUE)

# POPN G 2003-2004
XG3 <- matrix(c(0, 0, 0.245454545, 2.1,
0, 0, 0.045454545, 0,
0.125, 0, 0.090909091, 0,
0.125, 0, 0.090909091, 0.333333333), 4, 4, byrow = TRUE)

# POPN G 2004-2005
XG4 <- matrix(c(0, 0, 1.1, 1.54,
0.111111111, 0, 0, 0,
0, 0, 0, 0,
0.111111111, 0, 0, 0), 4, 4, byrow = TRUE)

# POPN G 2005-2006
XG5 <- matrix(c(0, 0, 0, 1.5,
0, 0, 0, 0,
0.090909091, 0, 0, 0,
0.545454545, 0.5, 0, 0.5), 4, 4, byrow = TRUE)

# POPN H (EXCLDED FROM ANALYSIS B/C OF UNREALISTIC ELASTICITIES)
XH3 <- matrix(c(0, 0, 0.1125, 1.05,
0.2, 0, 0, 0,
0, 0.5, 0, 0,
0.2, 0.5, 0, 0), 4, 4, byrow = TRUE)

XH4 <- matrix(c(0, 0, 0, 0,
0, 0, 0.5, 0,
0.8, 0.5, 0.25, 0.25,
0.2, 0, 0, 0.75), 4, 4, byrow = TRUE)

XH5 <- matrix(c(0, 0, 0.2, 1.05,
0, 0, 0, 0,
0.001, 0.001, 0.333333333, 0, #ELEMENTS (3,1),(4,1),(3,2) REPLACED W NONZERO
0.001, 0, 0, 0), 4, 4, byrow = TRUE)

# POPN L 2003-2004
XL3 <- matrix(c(0, 0, 1.785365854, 1.856521739,
0.128571429, 0, 0, 0.010869565,
0.028571429, 0, 0, 0,
0.014285714, 0, 0, 0.02173913), 4, 4, byrow = TRUE)

# POPN L 2004-2005
XL4 <- matrix(c(0, 0, 14.25, 16.625,
0.131443299, 0.057142857, 0, 0.25,
0.144329897, 0, 0, 0,
0.092783505, 0.2, 0, 0.25), 4, 4, byrow = TRUE)

# POPN L 2005-2006
XL5 <- matrix(c(0, 0, 0.594642857, 1.765909091,
0, 0, 0.017857143, 0,
0.021052632, 0.018518519, 0.035714286, 0.045454545,
0.021052632, 0.018518519, 0.035714286, 0.068181818), 4, 4, byrow = TRUE)

# POPN O 2003-2004
XO3 <- matrix(c(0, 0, 11.5, 2.775862069,
0.6, 0.285714286, 0.333333333, 0.24137931,
0.04, 0.142857143, 0, 0,
0.16, 0.285714286, 0, 0.172413793), 4, 4, byrow = TRUE)

# POPN O 2004-2005
XO4 <- matrix(c(0, 0, 3.78, 1.225,
0.28358209, 0.171052632, 0, 0.166666667,
0.084577114, 0.026315789, 0, 0.055555556,
0.139303483, 0.447368421, 0, 0.305555556), 4, 4, byrow = TRUE)

# POPN O 2005-2006
XO5 <- matrix(c(0, 0, 1.542857143, 1.035616438,
0.126984127, 0.105263158, 0.047619048, 0.054794521,
0.095238095, 0.157894737, 0.19047619, 0.082191781,
0.111111111, 0.223684211, 0, 0.356164384), 4, 4, byrow = TRUE)

# POPN Q 2003-2004
XQ3 <- matrix(c(0, 0, 0.15, 0.175,
0, 0, 0, 0,
0, 0, 0, 0,
1, 0, 0, 0), 4, 4, byrow = TRUE)

# POPN Q 2004-2005
XQ4 <- matrix(c(0, 0, 0, 0.25,
0, 0, 0, 0,
0, 0, 0, 0,
1, 0.666666667, 0, 1), 4, 4, byrow = TRUE)

# POPN Q 2005-2006
XQ5 <- matrix(c(0, 0, 0, 1.428571429,
0, 0, 0, 0.142857143,
0.25, 0, 0, 0,
0.25, 0, 0, 0.571428571), 4, 4, byrow = TRUE)

# POPN R 2003-2004
XR3 <- matrix(c(0, 0, 0.7, 0.6125,
0.25, 0, 0, 0.125,
0, 0, 0, 0,
0.25, 0.166666667, 0, 0.25), 4, 4, byrow = TRUE)

# POPN R 2004-2005
XR4 <- matrix(c(0, 0, 0, 0.6,
0.285714286, 0, 0, 0,
0.285714286, 0.333333333, 0, 0,
0.285714286, 0.333333333, 0, 1), 4, 4, byrow = TRUE)

# POPN R 2005-2006
XR5 <- matrix(c(0, 0, 0.7, 0.6125,
0, 0, 0, 0,
0, 0, 0, 0,
0.333333333, 0, 0.333333333, 0.625), 4, 4, byrow = TRUE)

## ----Ch9.4--------------------------------------------------------------------
mats_list <- list(XC3, XC4, XC5, XE3, XE4, XE5, XF3, XF4, XF5, XG3, XG4, XG5,
  XH3, XH4, XH5, XL3, XL4, XL5, XO3, XO4, XO5, XQ3, XQ4, XQ5, XR3, XR4, XR5)
#mats_list

## ----Ch9.5--------------------------------------------------------------------
pch_ord <- c("C", "C", "C", "E", "E", "E", "F", "F", "F", "G", "G", "G", "H",
  "H", "H", "L", "L", "L", "O", "O", "O", "Q", "Q", "Q", "R", "R", "R")
yr_ord <- c(2003, 2004, 2005, 2003, 2004, 2005, 2003, 2004, 2005, 2003, 2004,
  2005, 2003, 2004, 2005, 2003, 2004, 2005, 2003, 2004, 2005, 2003, 2004, 2005,
  2003, 2004, 2005)

anth_lM1 <- create_lM(mats_list, anthframe, historical = FALSE, poporder = 1,
  patchorder = pch_ord, yearorder = yr_ord)
#anth_lM1

## ----Ch9.6--------------------------------------------------------------------
summary(anth_lM1)

## ----Ch9.7--------------------------------------------------------------------
# POPN S 2003-2004
XS3 <- matrix(c(0, 0, 2.1, 0.816666667,
0.166666667, 0, 0, 0,
0, 0, 0, 0,
0, 0, 0, 0.166666667), 4, 4, byrow = TRUE)

# POPN S 2004-2005
XS4 <- matrix(c(0, 0, 0, 7,
0.333333333, 0.5, 0, 0,
0, 0, 0, 0,
0.333333333, 0, 0, 1), 4, 4, byrow = TRUE)

# POPN S 2005-2006
XS5 <- matrix(c(0, 0, 0, 1.4,
0, 0, 0, 0,
0, 0, 0, 0.2,
0.111111111, 0.75, 0, 0.2), 4, 4, byrow = TRUE)

anth_lM2 <- add_lM(anth_lM1, Amats = list(XS3, XS4, XS5), UFdecomp = TRUE,
  patch = "S", year = c(2003, 2004, 2005))

summary(anth_lM2)

## ----Ch9.8--------------------------------------------------------------------
anth_lM2$labels

## ----Ch9.9--------------------------------------------------------------------
anth_lM3 <- subset_lM(anth_lM2, patch = "S")
anth_lM4 <- subset_lM(anth_lM2, year = 2004)

summary(anth_lM3)
summary(anth_lM4)

## ----Ch9.10-------------------------------------------------------------------
anth_finallM <- delete_lM(anth_lM2, patch = "H")
summary(anth_finallM)

## ----Ch9.11-------------------------------------------------------------------
anth_finallM$labels

## ----Ch9.12, fig.cap = "Figure 9.2. Deterministic vs. stochastic lambda"------
anth_lmean <- lmean(anth_finallM)

lambda2 <- lambda3(anth_finallM)
lambda2m <- lambda3(anth_lmean)
set.seed(42)
sl2 <- slambda3(anth_finallM) #Stochastic growth rate
sl2$expa <- exp(sl2$a)

plot(lambda ~ year2, data = subset(lambda2, patch == "C"), ylim = c(0, 2.5),xlab = "Year",
  ylab = expression(lambda), type = "l", col = "gray", lty= 2, lwd = 2, bty = "n")
lines(lambda ~ year2, data = subset(lambda2, patch == "E"), col = "gray", lty= 2, lwd = 2)
lines(lambda ~ year2, data = subset(lambda2, patch == "F"), col = "gray", lty= 2, lwd = 2)
lines(lambda ~ year2, data = subset(lambda2, patch == "G"), col = "gray", lty= 2, lwd = 2)
lines(lambda ~ year2, data = subset(lambda2, patch == "L"), col = "gray", lty= 2, lwd = 2)
lines(lambda ~ year2, data = subset(lambda2, patch == "O"), col = "gray", lty= 2, lwd = 2)
lines(lambda ~ year2, data = subset(lambda2, patch == "Q"), col = "gray", lty= 2, lwd = 2)
lines(lambda ~ year2, data = subset(lambda2, patch == "R"), col = "gray", lty= 2, lwd = 2)
lines(lambda ~ year2, data = subset(lambda2, patch == "S"), col = "gray", lty= 2, lwd = 2)
abline(a = lambda2m$lambda[1], b = 0, lty = 1, lwd = 4, col = "orangered")
abline(a = sl2$expa[1], b = 0, lty = 1, lwd = 4, col = "darkred")
legend("topleft", c("det annual", "det mean", "stochastic"), lty = c(2, 1, 1),
  col = c("gray", "orangered", "darkred"), lwd = c(2, 4, 4), bty = "n")


## ----Ch9.13-------------------------------------------------------------------
trialltre_det <- ltre3(anth_finallM, refmats = NA, stochastic = FALSE,
  steps = 10000, time_weights = NA, sparse = "auto")
#trialltre_det

## ----Ch9.14-------------------------------------------------------------------
trialltre_sto <- ltre3(anth_finallM, refmats = NA, stochastic = TRUE,
  steps = 10000, time_weights = NA, sparse = "auto")
#trialltre_sto

## ----Ch9.15-------------------------------------------------------------------
writeLines("The greatest (i.e most positive) deterministic LTRE contribution: ")
max(trialltre_det$ltre_det[[1]])
writeLines("The greatest deterministic LTRE contribution is associated with element: ")
which(trialltre_det$ltre_det[[1]] == max(trialltre_det$ltre_det[[1]]))
writeLines("The lowest (i.e. most negative) deterministic LTRE contribution: ")
min(trialltre_det$ltre_det[[1]])
writeLines("The lowest deterministic LTRE contribution is associated with element: ")
which(trialltre_det$ltre_det[[1]] == min(trialltre_det$ltre_det[[1]]))

writeLines("\nThe greatest (i.e most positive) stochastic mean LTRE contribution: ")
max(trialltre_sto$ltre_mean[[1]])
writeLines("The greatest stochastic mean LTRE contribution is associated with element: ")
which(trialltre_sto$ltre_mean[[1]] == max(trialltre_sto$ltre_mean[[1]]))
writeLines("The lowest (i.e. most negative) stochastic mean LTRE contribution: ")
min(trialltre_sto$ltre_mean[[1]])
writeLines("The lowest stochastic mean LTRE contribution is associated with element: ")
which(trialltre_sto$ltre_mean[[1]] == min(trialltre_sto$ltre_mean[[1]]))

writeLines("\nThe greatest (i.e most positive) stochastic SD LTRE contribution: ")
max(trialltre_sto$ltre_sd[[1]])
writeLines("The greatest stochastic SD LTRE contribution is associated with element: ")
which(trialltre_sto$ltre_sd[[1]] == max(trialltre_sto$ltre_sd[[1]]))
writeLines("The lowest (i.e. most negative) stochastic SD LTRE contribution: ")
min(trialltre_sto$ltre_sd[[1]])
writeLines("The lowest stochastic SD LTRE contribution is associated with element: ")
which(trialltre_sto$ltre_sd[[1]] == min(trialltre_sto$ltre_sd[[1]]))

writeLines("\nTotal positive contribution of shifts in deterministic LTRE contributions: ")
sum(trialltre_det$ltre_det[[1]][which(trialltre_det$ltre_det[[1]] > 0)])
writeLines("Total negative contribution of shifts in deterministic LTRE contributions: ")
sum(trialltre_det$ltre_det[[1]][which(trialltre_det$ltre_det[[1]] < 0)])
writeLines("\nTotal positive contribution of shifts in stochastic mean LTRE contributions: ")
sum(trialltre_sto$ltre_mean[[1]][which(trialltre_sto$ltre_mean[[1]] > 0)])
writeLines("Total negative contribution of shifts in stochastic mean LTRE contributions: ")
sum(trialltre_sto$ltre_mean[[1]][which(trialltre_sto$ltre_mean[[1]] < 0)])
writeLines("\nTotal positive contribution of shifts in stochastic SD LTRE contributions: ")
sum(trialltre_sto$ltre_sd[[1]][which(trialltre_sto$ltre_sd[[1]] > 0)])
writeLines("Total negative contribution of shifts in stochastic SD LTRE contributions: ")
sum(trialltre_sto$ltre_sd[[1]][which(trialltre_sto$ltre_sd[[1]] < 0)])

## ----Ch9.16, fig.cap = "Figure 9.3. LTRE contributions by stage in population C"----
ltre_pos <- trialltre_det$ltre_det[[1]]
ltre_neg <- trialltre_det$ltre_det[[1]]
ltre_pos[which(ltre_pos < 0)] <- 0
ltre_neg[which(ltre_neg > 0)] <- 0

sltre_meanpos <- trialltre_sto$ltre_mean[[1]]
sltre_meanneg <- trialltre_sto$ltre_mean[[1]]
sltre_meanpos[which(sltre_meanpos < 0)] <- 0
sltre_meanneg[which(sltre_meanneg > 0)] <- 0

sltre_sdpos <- trialltre_sto$ltre_sd[[1]]
sltre_sdneg <- trialltre_sto$ltre_sd[[1]]
sltre_sdpos[which(sltre_sdpos < 0)] <- 0
sltre_sdneg[which(sltre_sdneg > 0)] <- 0

ltresums_pos <- cbind(colSums(ltre_pos), colSums(sltre_meanpos), colSums(sltre_sdpos))
ltresums_neg <- cbind(colSums(ltre_neg), colSums(sltre_meanneg), colSums(sltre_sdneg))

ltre_as_names <- trialltre_det$ahstages$stage

barplot(t(ltresums_pos), beside = T, col = c("black", "grey", "red"),
  ylim = c(-0.40, 0.10))
barplot(t(ltresums_neg), beside = T, col = c("black", "grey", "red"), add = TRUE)
abline(0, 0, lty= 3)
text(cex=1, y = -0.45, x = seq(from = 2, to = 3.98*length(ltre_as_names),
    by = 4), ltre_as_names, xpd=TRUE, srt=45)
legend("bottomleft", c("deterministic", "stochastic mean", "stochastic SD"),
  col = c("black", "grey", "red"), pch = 15, bty = "n")

## ----Ch9.17, fig.cap = "Figure 9.4. LTRE contributions by transition type in population C"----
det_ltre_summary <- summary(trialltre_det)
sto_ltre_summary <- summary(trialltre_sto)

ltresums_tpos <- cbind(det_ltre_summary$ahist_det$matrix1_pos,
  sto_ltre_summary$ahist_mean$matrix1_pos,
  sto_ltre_summary$ahist_sd$matrix1_pos)
ltresums_tneg <- cbind(det_ltre_summary$ahist_det$matrix1_neg,
  sto_ltre_summary$ahist_mean$matrix1_neg,
  sto_ltre_summary$ahist_sd$matrix1_neg)

barplot(t(ltresums_tpos), beside = T, col = c("black", "grey", "red"),
  ylim = c(-0.50, 0.10))
barplot(t(ltresums_tneg), beside = T, col = c("black", "grey", "red"),
  add = TRUE)
abline(0, 0, lty = 3)
text(cex=0.85, y = -0.54, x = seq(from = 2, to = 3.98*length(det_ltre_summary$ahist_det$category),
    by = 4), det_ltre_summary$ahist_det$category, xpd=TRUE, srt=45)
legend("bottomleft", c("deterministic", "stochastic mean", "stochastic SD"),
  col = c("black", "grey", "red"), pch = 15, bty = "n")

## ----Ch9.18-------------------------------------------------------------------
sltre_meanpos_toplot <- sto_ltre_summary$ahist_mean[,seq(from = 3, by = 3, length.out = 9)]
sltre_meanneg_toplot <- sto_ltre_summary$ahist_mean[,seq(from = 4, by = 3, length.out = 9)]
sltre_sdpos_toplot <- sto_ltre_summary$ahist_sd[,seq(from = 3, by = 3, length.out = 9)]
sltre_sdneg_toplot <- sto_ltre_summary$ahist_sd[,seq(from = 4, by = 3, length.out = 9)]

colnames(sltre_meanpos_toplot) <- trialltre_sto$labels$patch
colnames(sltre_meanneg_toplot) <- trialltre_sto$labels$patch
colnames(sltre_sdpos_toplot) <- trialltre_sto$labels$patch
colnames(sltre_sdneg_toplot) <- trialltre_sto$labels$patch

rownames(sltre_meanpos_toplot) <- c("Stasis", "Growth", "Shrinkage", "Fecundity")
rownames(sltre_meanneg_toplot) <- c("Stasis", "Growth", "Shrinkage", "Fecundity")
rownames(sltre_sdpos_toplot) <- c("Stasis", "Growth", "Shrinkage", "Fecundity")
rownames(sltre_sdneg_toplot) <- c("Stasis", "Growth", "Shrinkage", "Fecundity")

sltre_meanpos_toplot <- sltre_meanpos_toplot[c(4, 3, 1, 2),]
sltre_meanneg_toplot <- sltre_meanneg_toplot[c(4, 3, 1, 2),]
sltre_sdpos_toplot <- sltre_sdpos_toplot[c(4, 3, 1, 2),]
sltre_sdneg_toplot <- sltre_sdneg_toplot[c(4, 3, 1, 2),]

## ----Ch9.19, fig.cap = "Figure 9.5. Mean element sLTRE contributions by transition type"----
barplot(as.matrix(sltre_meanpos_toplot), main = "sLTRE impacts of shifts in mean", 
  xlab = "Population", ylab = "Contribution to stochastic growth rate",
  ylim = c(-0.8, 0.6), col = c("darkblue", "cyan", "yellow", "darkred"),
  legend.text = rownames(sltre_meanpos_toplot), beside = FALSE, args.legend = c(x = "topleft", bty = "n"))
barplot(as.matrix(sltre_meanneg_toplot), main = "sLTRE impacts of shifts in mean", 
  add = TRUE, col = c("darkblue", "cyan", "yellow", "darkred"),
  beside = FALSE)
abline(0, 0, lty = 2)

## ----Ch9.20, fig.cap = "Figure 9.6. Element SD sLTRE contributions by transition type"----
barplot(as.matrix(sltre_sdpos_toplot), main = "sLTRE impacts of shifts in SD", 
  xlab = "Population", ylab = "Contribution to stochastic growth rate",
  ylim = c(-0.06, 0.13), col = c("darkblue", "cyan", "yellow", "darkred"),
  legend.text = rownames(sltre_sdpos_toplot), beside = FALSE, args.legend = c(x = "topleft", bty = "n"))
barplot(as.matrix(sltre_sdneg_toplot), main = "sLTRE impacts of shifts in mean", 
  add = TRUE, col = c("darkblue", "cyan", "yellow", "darkred"),
  beside = FALSE)
abline(0, 0, lty = 2)

