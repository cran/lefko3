## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(lefko3)

## ----Ch9.1--------------------------------------------------------------------
sizevector <- c(1, 1, 2, 3) # These sizes are not from the original paper
stagevector <- c("Sdl", "Veg", "SmFlo", "LFlo")
repvector <- c(0, 0, 1, 1)
obsvector <- c(1, 1, 1, 1)
matvector <- c(0, 1, 1, 1)
immvector <- c(1, 0, 0, 0)
propvector <- c(0, 0, 0, 0)
indataset <- c(1, 1, 1, 1)
binvec <- c(0.5, 0.5, 0.5, 0.5)
comments <- c("Seedling", "Vegetative adult", "Small flowering",
  "Large flowering")

anthframe <- sf_create(sizes = sizevector, stagenames = stagevector,
  repstatus = repvector, obsstatus = obsvector, matstatus = matvector,
  immstatus = immvector, indataset = indataset, binhalfwidth = binvec,
  propstatus = propvector, comments = comments)
#anthframe

## ----Ch9.2--------------------------------------------------------------------
XC3 <- matrix(c(0, 0, 1.74, 1.74,     # POPN C 2003-2004
  0.208333333, 0, 0, 0.057142857,
  0.041666667, 0.076923077, 0, 0,
  0.083333333, 0.076923077, 0.066666667, 0.028571429), 4, 4, byrow = TRUE)
XC3

## ----Ch9.3--------------------------------------------------------------------
XC4 <- matrix(c(0, 0, 0.3, 0.6,     # POPN C 2004-2005
  0.32183908, 0.142857143, 0, 0,
  0.16091954, 0.285714286, 0, 0,
  0.252873563, 0.285714286, 0.5, 0.6), 4, 4, byrow = TRUE)
XC5 <- matrix(c(0, 0, 0.50625, 0.675,     # POPN C 2005-2006
  0, 0, 0, 0.035714286,
  0.1, 0.068965517, 0.0625, 0.107142857,
  0.3, 0.137931034, 0, 0.071428571), 4, 4, byrow = TRUE)

XE3 <- matrix(c(0, 0, 2.44, 6.569230769,     # POPN E 2003-2004
  0.196428571, 0, 0, 0,
  0.125, 0.5, 0, 0,
  0.160714286, 0.5, 0.133333333, 0.076923077), 4, 4, byrow = TRUE)
XE4 <- matrix(c(0, 0, 0.45, 0.646153846,     # POPN E 2004-2005
  0.06557377, 0.090909091, 0.125, 0,
  0.032786885, 0, 0.125, 0.076923077,
  0.049180328, 0, 0.125, 0.230769231), 4, 4, byrow = TRUE)
XE5 <- matrix(c(0, 0, 2.85, 3.99,     # POPN E 2005-2006
  0.083333333, 0, 0, 0,
  0, 0, 0, 0,
  0.416666667, 0.1, 0, 0.1), 4, 4, byrow = TRUE)

XF3 <- matrix(c(0, 0, 1.815, 7.058333333,     # POPN F 2003-2004
  0.075949367, 0, 0.05, 0.083333333,
  0.139240506, 0, 0, 0.25,
  0.075949367, 0, 0, 0.083333333), 4, 4, byrow = TRUE)
XF4 <- matrix(c(0, 0, 1.233333333, 7.4,     # POPN F 2004-2005
  0.223880597, 0, 0.111111111, 0.142857143,
  0.134328358, 0.272727273, 0.166666667, 0.142857143,
  0.119402985, 0.363636364, 0.055555556, 0.142857143), 4, 4, byrow = TRUE)
XF5 <- matrix(c(0, 0, 1.06, 3.372727273,     # POPN F 2005-2006
  0.073170732, 0.025, 0.033333333, 0,
  0.036585366, 0.15, 0.1, 0.136363636,
  0.06097561, 0.225, 0.166666667, 0.272727273), 4, 4, byrow = TRUE)

XG3 <- matrix(c(0, 0, 0.245454545, 2.1,     # POPN G 2003-2004
  0, 0, 0.045454545, 0,
  0.125, 0, 0.090909091, 0,
  0.125, 0, 0.090909091, 0.333333333), 4, 4, byrow = TRUE)
XG4 <- matrix(c(0, 0, 1.1, 1.54,     # POPN G 2004-2005
  0.111111111, 0, 0, 0,
  0, 0, 0, 0,
  0.111111111, 0, 0, 0), 4, 4, byrow = TRUE)
XG5 <- matrix(c(0, 0, 0, 1.5,     # POPN G 2005-2006
  0, 0, 0, 0,
  0.090909091, 0, 0, 0,
  0.545454545, 0.5, 0, 0.5), 4, 4, byrow = TRUE)

XL3 <- matrix(c(0, 0, 1.785365854, 1.856521739,     # POPN L 2003-2004
  0.128571429, 0, 0, 0.010869565,
  0.028571429, 0, 0, 0,
  0.014285714, 0, 0, 0.02173913), 4, 4, byrow = TRUE)
XL4 <- matrix(c(0, 0, 14.25, 16.625,     # POPN L 2004-2005
  0.131443299, 0.057142857, 0, 0.25,
  0.144329897, 0, 0, 0,
  0.092783505, 0.2, 0, 0.25), 4, 4, byrow = TRUE)
XL5 <- matrix(c(0, 0, 0.594642857, 1.765909091,     # POPN L 2005-2006
  0, 0, 0.017857143, 0,
  0.021052632, 0.018518519, 0.035714286, 0.045454545,
  0.021052632, 0.018518519, 0.035714286, 0.068181818), 4, 4, byrow = TRUE)

XO3 <- matrix(c(0, 0, 11.5, 2.775862069,     # POPN O 2003-2004
  0.6, 0.285714286, 0.333333333, 0.24137931,
  0.04, 0.142857143, 0, 0,
  0.16, 0.285714286, 0, 0.172413793), 4, 4, byrow = TRUE)
XO4 <- matrix(c(0, 0, 3.78, 1.225,     # POPN O 2004-2005
  0.28358209, 0.171052632, 0, 0.166666667,
  0.084577114, 0.026315789, 0, 0.055555556,
  0.139303483, 0.447368421, 0, 0.305555556), 4, 4, byrow = TRUE)
XO5 <- matrix(c(0, 0, 1.542857143, 1.035616438,     # POPN O 2005-2006
  0.126984127, 0.105263158, 0.047619048, 0.054794521,
  0.095238095, 0.157894737, 0.19047619, 0.082191781,
  0.111111111, 0.223684211, 0, 0.356164384), 4, 4, byrow = TRUE)

XQ3 <- matrix(c(0, 0, 0.15, 0.175,     # POPN Q 2003-2004
  0, 0, 0, 0,
  0, 0, 0, 0,
  1, 0, 0, 0), 4, 4, byrow = TRUE)
XQ4 <- matrix(c(0, 0, 0, 0.25,     # POPN Q 2004-2005
  0, 0, 0, 0,
  0, 0, 0, 0,
  1, 0.666666667, 0, 1), 4, 4, byrow = TRUE)
XQ5 <- matrix(c(0, 0, 0, 1.428571429,     # POPN Q 2005-2006
  0, 0, 0, 0.142857143,
  0.25, 0, 0, 0,
  0.25, 0, 0, 0.571428571), 4, 4, byrow = TRUE)

XR3 <- matrix(c(0, 0, 0.7, 0.6125,     # POPN R 2003-2004
  0.25, 0, 0, 0.125,
  0, 0, 0, 0,
  0.25, 0.166666667, 0, 0.25), 4, 4, byrow = TRUE)
XR4 <- matrix(c(0, 0, 0, 0.6,     # POPN R 2004-2005
  0.285714286, 0, 0, 0,
  0.285714286, 0.333333333, 0, 0,
  0.285714286, 0.333333333, 0, 1), 4, 4, byrow = TRUE)
XR5 <- matrix(c(0, 0, 0.7, 0.6125,     # POPN R 2005-2006
  0, 0, 0, 0,
  0, 0, 0, 0,
  0.333333333, 0, 0.333333333, 0.625), 4, 4, byrow = TRUE)

XS3 <- matrix(c(0, 0, 2.1, 0.816666667,     # POPN S 2003-2004
  0.166666667, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0.166666667), 4, 4, byrow = TRUE)
XS4 <- matrix(c(0, 0, 0, 7,     # POPN S 2004-2005
  0.333333333, 0.5, 0, 0,
  0, 0, 0, 0,
  0.333333333, 0, 0, 1), 4, 4, byrow = TRUE)
XS5 <- matrix(c(0, 0, 0, 1.4,     # POPN S 2005-2006
  0, 0, 0, 0,
  0, 0, 0, 0.2,
  0.111111111, 0.75, 0, 0.2), 4, 4, byrow = TRUE)

## ----Ch9.4--------------------------------------------------------------------
mats_list <- list(XC3, XC4, XC5, XE3, XE4, XE5, XF3, XF4, XF5, XG3, XG4, XG5,
  XL3, XL4, XL5, XO3, XO4, XO5, XQ3, XQ4, XQ5, XR3, XR4, XR5, XS3, XS4, XS5)

anth_lefkoMat <- create_lM(mats_list, anthframe, hstages = NA,
  historical = FALSE, poporder = 1,
  patchorder = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7,
    8, 8, 8, 9, 9, 9),
  yearorder = c(2003, 2004, 2005, 2003, 2004, 2005, 2003, 2004, 2005, 2003,
    2004, 2005, 2003, 2004, 2005, 2003, 2004, 2005, 2003, 2004, 2005, 2003,
    2004, 2005, 2003, 2004, 2005))
#anth_lefkoMat

## ----Ch9.5--------------------------------------------------------------------
#summary(anth_lefkoMat)

## ----Ch9.6, fig.cap = "Figure 9.2. Deterministic vs. stochastic lambda"-------
anth_lmean <- lmean(anth_lefkoMat)

lambda2 <- lambda3(anth_lefkoMat)
lambda2m <- lambda3(anth_lmean)
set.seed(42)
sl2 <- slambda3(anth_lefkoMat) #Stochastic growth rate
sl2$expa <- exp(sl2$a)

plot(lambda ~ year2, data = subset(lambda2, patch == 1), ylim = c(0, 2.5),xlab = "Year",
  ylab = expression(lambda), type = "l", col = "gray", lty= 2, lwd = 2, bty = "n")
lines(lambda ~ year2, data = subset(lambda2, patch == 2), col = "gray", lty= 2, lwd = 2)
lines(lambda ~ year2, data = subset(lambda2, patch == 3), col = "gray", lty= 2, lwd = 2)
lines(lambda ~ year2, data = subset(lambda2, patch == 4), col = "gray", lty= 2, lwd = 2)
lines(lambda ~ year2, data = subset(lambda2, patch == 5), col = "gray", lty= 2, lwd = 2)
lines(lambda ~ year2, data = subset(lambda2, patch == 6), col = "gray", lty= 2, lwd = 2)
lines(lambda ~ year2, data = subset(lambda2, patch == 7), col = "gray", lty= 2, lwd = 2)
lines(lambda ~ year2, data = subset(lambda2, patch == 8), col = "gray", lty= 2, lwd = 2)
lines(lambda ~ year2, data = subset(lambda2, patch == 9), col = "gray", lty= 2, lwd = 2)
abline(a = lambda2m$lambda[1], b = 0, lty = 1, lwd = 4, col = "orangered")
abline(a = sl2$expa[1], b = 0, lty = 1, lwd = 4, col = "darkred")
legend("topleft", c("det annual", "det mean", "stochastic"), lty = c(2, 1, 1),
  col = c("gray", "orangered", "darkred"), lwd = c(2, 4, 4), bty = "n")

## ----Ch9.7--------------------------------------------------------------------
trialltre_det <- ltre3(anth_lmean, sparse = "auto")
#trialltre_det

## ----Ch9.8--------------------------------------------------------------------
trialltre_sto <- ltre3(anth_lefkoMat, stochastic = TRUE, times = 10000,
  tweights = NA, sparse = "auto", seed = 42)
#trialltre_sto

## ----Ch9.9--------------------------------------------------------------------
# Highest (i.e most positive) deterministic LTRE contribution:
max(trialltre_det$cont_mean[[1]])
# Highest deterministic LTRE contribution is associated with element:
which(trialltre_det$cont_mean[[1]] == max(trialltre_det$cont_mean[[1]]))
# Lowest (i.e. most negative) deterministic LTRE contribution:
min(trialltre_det$cont_mean[[1]])
# Lowest deterministic LTRE contribution is associated with element:
which(trialltre_det$cont_mean[[1]] == min(trialltre_det$cont_mean[[1]]))

# Highest stochastic mean LTRE contribution:
max(trialltre_sto$cont_mean[[1]])
# Highest stochastic mean LTRE contribution is associated with element:
which(trialltre_sto$cont_mean[[1]] == max(trialltre_sto$cont_mean[[1]]))
# Lowest stochastic mean LTRE contribution:
min(trialltre_sto$cont_mean[[1]])
# Lowest stochastic mean LTRE contribution is associated with element:
which(trialltre_sto$cont_mean[[1]] == min(trialltre_sto$cont_mean[[1]]))

# Highest stochastic SD LTRE contribution:
max(trialltre_sto$cont_sd[[1]])
# Highest stochastic SD LTRE contribution is associated with element:
which(trialltre_sto$cont_sd[[1]] == max(trialltre_sto$cont_sd[[1]]))
# Lowest stochastic SD LTRE contribution:
min(trialltre_sto$cont_sd[[1]])
# Lowest stochastic SD LTRE contribution is associated with element:
which(trialltre_sto$cont_sd[[1]] == min(trialltre_sto$cont_sd[[1]]))

# Total positive deterministic LTRE contributions:
sum(trialltre_det$cont_mean[[1]][which(trialltre_det$cont_mean[[1]] > 0)])
# Total negative deterministic LTRE contributions:
sum(trialltre_det$cont_mean[[1]][which(trialltre_det$cont_mean[[1]] < 0)])
# Total positive stochastic mean LTRE contributions:
sum(trialltre_sto$cont_mean[[1]][which(trialltre_sto$cont_mean[[1]] > 0)])
# Total negative stochastic mean LTRE contributions:
sum(trialltre_sto$cont_mean[[1]][which(trialltre_sto$cont_mean[[1]] < 0)])
# Total positive stochastic SD LTRE contributions:
sum(trialltre_sto$cont_sd[[1]][which(trialltre_sto$cont_sd[[1]] > 0)])
# Total negative stochastic SD LTRE contributions:
sum(trialltre_sto$cont_sd[[1]][which(trialltre_sto$cont_sd[[1]] < 0)])

## ----Ch9.10, fig.cap = "Figure 9.3. LTRE contributions by stage"--------------
ltre_pos <- trialltre_det$cont_mean[[1]]
ltre_neg <- trialltre_det$cont_mean[[1]]
ltre_pos[which(ltre_pos < 0)] <- 0
ltre_neg[which(ltre_neg > 0)] <- 0

sltre_meanpos <- trialltre_sto$cont_mean[[1]]
sltre_meanneg <- trialltre_sto$cont_mean[[1]]
sltre_meanpos[which(sltre_meanpos < 0)] <- 0
sltre_meanneg[which(sltre_meanneg > 0)] <- 0

sltre_sdpos <- trialltre_sto$cont_sd[[1]]
sltre_sdneg <- trialltre_sto$cont_sd[[1]]
sltre_sdpos[which(sltre_sdpos < 0)] <- 0
sltre_sdneg[which(sltre_sdneg > 0)] <- 0

ltresums_pos <- cbind(colSums(ltre_pos), colSums(sltre_meanpos), colSums(sltre_sdpos))
ltresums_neg <- cbind(colSums(ltre_neg), colSums(sltre_meanneg), colSums(sltre_sdneg))

ltre_as_names <- trialltre_det$ahstages$stage

barplot(t(ltresums_pos), beside = T, col = c("black", "grey", "red"),
  ylim = c(-0.50, 0.10))
barplot(t(ltresums_neg), beside = T, col = c("black", "grey", "red"), add = TRUE)
abline(0, 0, lty= 3)
text(cex=1, y = -0.57, x = seq(from = 2, to = 3.98*length(ltre_as_names),
    by = 4), ltre_as_names, xpd=TRUE, srt=45)
legend("bottomleft", c("deterministic", "stochastic mean", "stochastic SD"),
  col = c("black", "grey", "red"), pch = 15, bty = "n")

## ----Ch9.11, fig.cap = "Figure 9.4. LTRE contributions by transition type"----
det_ltre_summary <- summary(trialltre_det)
sto_ltre_summary <- summary(trialltre_sto)

ltresums_tpos <- cbind(det_ltre_summary$ahist_mean$matrix1_pos,
  sto_ltre_summary$ahist_mean$matrix1_pos,
  sto_ltre_summary$ahist_sd$matrix1_pos)
ltresums_tneg <- cbind(det_ltre_summary$ahist_mean$matrix1_neg,
  sto_ltre_summary$ahist_mean$matrix1_neg,
  sto_ltre_summary$ahist_sd$matrix1_neg)

barplot(t(ltresums_tpos), beside = T, col = c("black", "grey", "red"),
  ylim = c(-0.55, 0.10))
barplot(t(ltresums_tneg), beside = T, col = c("black", "grey", "red"),
  add = TRUE)
abline(0, 0, lty = 3)
text(cex=0.85, y = -0.64, x = seq(from = 2,
    to = 3.98*length(det_ltre_summary$ahist_mean$category), by = 4),
    det_ltre_summary$ahist_mean$category, xpd=TRUE, srt=45)
legend("bottomleft", c("deterministic", "stochastic mean", "stochastic SD"),
  col = c("black", "grey", "red"), pch = 15, bty = "n")

