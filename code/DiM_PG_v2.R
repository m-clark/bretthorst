################################################################################
# ON THE DIFFERENCE OF MEANS
################################################################################

# original Mathematica code by Phil Gregory
# http://www.phas.ubc.ca/~gregory/gregory.html
# Ch. 9: Bayesian analysis of two independent samples
# Introduction
# This is a Mathematrica implementation of the probability calculations discussed in the book in the section entitled,     
# "Bayesian Comparison of Two Samples?". 
#
# based on a paper from G.L. Bretthorst (1993) "on the difference of means"
# http://bayes.wustl.edu/glb/diff.pdf
#
# see also Mathematica code by UM Studer (1996 + 1998) on the same problem (paper)
# taken code from there to allow calculations based only on summary statistics
# and not on raw data (see also Bretthorst, 1993, for an example taken from Jaynes)

# R code by Leo G?rtler 2017
# first = 2017-04-19
# last = 2017-27-04

# notes:
# *- introduce logs to integral calculations, but probably that won't help...
# *- very small numbers are slightly different from Mathematica -> e.g. e-230


################################################################################
# pre-calculations taken from UMS 1998
ums2pg <- function(inputvalues)
  #convert values from UMS implementation to PG implementation
{
  #mapping
  #n
  #NN <- Ni+Nii
  #
  #dbar
  #DD <- (Ni * Di + Nii * Dii) / NN
  #
  #d1squbar
  #Dsi <- (Ni-1) / Ni * si^2 + Di^2
  #
  #d2squbar
  #Dsii <- (Nii-1) / Nii * sii^2 + Dii^2
  #
  #dsqubar
  #DsD <- (Ni * Dsi + Nii * Dsii) / NN
  #
  #dstd
  #ss <- sqrt(NN * (DsD - DD^2) / (NN-1))
  
  d1bar <- inputvalues[["Di"]]
  d2bar <- inputvalues[["Dii"]]
  d1std <- inputvalues[["si"]]
  d2std <- inputvalues[["sii"]]
  n1 <- inputvalues[["Ni"]]
  n2 <- inputvalues[["Nii"]]
  low <- inputvalues[["L"]]
  high <- inputvalues[["H"]]
  sigma.low <- inputvalues[["sL"]]
  sigma.high <- inputvalues[["sH"]]
  snames <- inputvalues[["snames"]]
  ndelta <- inputvalues[["ndelta"]]
  nr <- inputvalues[["nr"]]
  
  ###########################
  # manual update of the data
  n <- n1 + n2
  dbar <- (d1bar * n1 + d2bar * n2) / n#(Ni+Nii)
  # Compute the mean square average of the 1st data set
  d1squbar <- (n1 - 1) / n1 * d1std ^ 2 + d1bar ^ 2
  # Compute the mean square average of the 2nd data set
  d2squbar <- (n2 - 1) / n2 * d2std ^ 2 + d2bar ^ 2
  # Compute the mean square average of the combined data set
  dsqubar <- (n1 * d1squbar + n2 * d2squbar) / n
  # Compute the standard deviation of the combined data set
  dstd <- sqrt(n * (dsqubar - dbar ^ 2) / (n - 1))
  ##########
  
  d1bar
  d2bar
  dbar
  d1squbar
  d2squbar
  d1std
  d2std
  dstd
  
  res <- list(
    snames = snames,
    n1 = n1,
    n2 = n2,
    n = n,
    low = low,
    high = high,
    sigma.low = sigma.low,
    sigma.high = sigma.high,
    d1bar = d1bar,
    d2bar = d2bar,
    dbar = dbar,
    d1std = d1std,
    d2std = d2std,
    dstd = dstd,
    d1squbar = d1squbar,
    d2squbar = d2squbar,
    dsqubar = dsqubar,
    ndelta = ndelta,
    nr = nr
  )
  
  return(res)
}
################################################################################
# call according to UMS scheme
#inputvalues <- list(snames=c("Jaynes.1","Jaynes.2"), si=6.48, Ni=4, sii=7.48, Nii=9, Di=50, Dii=42, L=34, H=58, sL=3, sH=10, ndelta=1000, nr=1000)
#inputvalues
#ums2pg(inputvalues)
#dim.res <- DiM.pg(invtyp="ums", inputvalues, print.res=TRUE)
#dim.res
#plot.DiM(DiM.res=dim.res, filling=TRUE)


################################################################################
DiM.pg <- function(invtyp=NULL, inputvalues=NULL, print.res=FALSE, graph=TRUE, dig=4)
{
#Bretthorst DiM
#adapted from Mathematica code by Phil Gregory

####################
# START INPUT VALUES
if(invtyp=="ums")
{
 d1 <- NULL
 d2 <- NULL
 
 ums2pg.res <- ums2pg(inputvalues)
 ums2pg.res
 snames <- ums2pg.res[["snames"]]

 n1 <- ums2pg.res[["n1"]]
 n2 <- ums2pg.res[["n2"]]
 n <- ums2pg.res[["n"]]

 low <- ums2pg.res[["low"]]
 high <- ums2pg.res[["high"]]
 sigma.low <- ums2pg.res[["sigma.low"]]
 sigma.high <- ums2pg.res[["sigma.high"]]

 d1bar <- ums2pg.res[["d1bar"]]
 d2bar <- ums2pg.res[["d2bar"]]
 dbar <- ums2pg.res[["dbar"]]
 
 d1std <- ums2pg.res[["d1std"]]
 d2std <- ums2pg.res[["d2std"]]
 dstd <- ums2pg.res[["dstd"]]
 
 d1squbar <- ums2pg.res[["d1squbar"]]
 d2squbar <- ums2pg.res[["d2squbar"]]
 dsqubar <- ums2pg.res[["dsqubar"]]

 ndelta <- inputvalues[["ndelta"]]
 nr <- inputvalues[["nr"]]
  
} else if(invtyp=="pg") {

################################################################################
# input values PG scheme

 snames <- inputvalues[["snames"]]
 d1 <- inputvalues[["d1"]]
 d2 <- inputvalues[["d2"]]
 ndelta <- inputvalues[["ndelta"]]
 nr <- inputvalues[["nr"]]
 high <- inputvalues[["high"]]
 low <- inputvalues[["low"]]
 sigma.low <- inputvalues[["sigma.low"]]
 sigma.high <- inputvalues[["sigma.high"]]

# sample sizes
 n1 <- length(d1)
 n2 <- length(d2)
 n <- n1 + n2

 # Compute the average of the 1st data set
 d1bar <- sum(d1)/n1 #mean
 # Compute the average of the 2nd data set
 d2bar <- sum(d2)/n2
 # Compute the average of the combined data set
 dbar <- sum(d1,d2)/n
 # Compute the mean square average of the 1st data set
 d1squbar <- sum(d1 * d1)/n1
 # Compute the mean square average of the 2nd data set
 d2squbar <- sum(d2 * d2)/n2
 # Compute the mean square average of the combined data set
 dsqubar <- (sum(d1 * d1) + sum(d2 * d2))/n
 # Neater way of computing the above using a vector dot product
 # (d1 %*% d1 + d2 %*% d2)/n
 # Compute the standard deviation of the 1st data set
 d1std <- sd(d1)
 # Compute the standard deviation of the 2nd data set
 d2std <- sd(d2)
 # Compute the standard deviation of the combined data set
 dstd <- sd(c(d1,d2))

 d1
 d2

} else {
 stop(paste("\nNo valid input type for input values given.\n\n",sep=""))
}

# Set prior limits (assumed the same for each data set) on mean (low,high),
# and prior limits (assumed the same for each data set) on the
# standard deviation (sigmalow, sigmahigh).
# range mean
Rc <- high - low
# range sd
Rsigma <- sigma.high / sigma.low


#summary input values regardless scheme
snames



n1
n2
n

low
high
sigma.low
sigma.high

d1bar
d2bar
dbar

d1std
d2std
dstd

d1squbar
d2squbar
dsqubar

Rc
Rsigma

ndelta
nr

# END INPUT VALUES
##################


################################################################################
# pCSk
# Compute pCSk = p(C,S|D_1,D_2,I) * p(D_1,D_2|I)
# According to the formulas given in text (see Appendix C  entitled: "Difference in Two Samples").
z <- n * (dsqubar - dbar^2)
uL <- sqrt(n/2) * (low - dbar)
uH <- sqrt(n/2) * (high - dbar)

#????????TODO
#sign <- -Sign[uH/uL]
(uH/uL) > 0

# error function
errf <- function(ERR) return( 2*pnorm(ERR*sqrt(2))-1 )

# formula to calculate area below the curve / integral
fnc1 <- function(sigma1)
{
 ( (2*pi)^(-n/2)*sqrt(pi/(2*n)) ) / ( 4*Rc*log(Rsigma) ) * sigma1^(-n) * exp(-z/(2*sigma1^2)) * ( errf(ERR=uH/sigma1) - errf(ERR=uL/sigma1) )
}
# integrate
pCSk <- integrate(fnc1, lower=sigma.low, upper=sigma.high)$value
pCSk


################################################################################
# pCSbark
# Compute pCSbark = p(C,Sbar|D_1,D_2,I) * p(D_1,D_2|I)
# Mathematica does not allow variable definitions containing a character with an overhead bar so we will use Sbar to represent S-with-a-bar 
# >>> pCS-with-a-bark becomes pCSbark
u1A <- function(A) return( n1/2 * (d1squbar - 2 * A * d1bar + A^2) )
u2A <- function(A) return( n2/2 * (d2squbar - 2 * A * d2bar + A^2) )
# formula to calculate area below the curve / integral
fnc2 <- function(A)
{
 (2*pi)^(-n/2)*gamma(n1/2)*gamma(n2/2) / ( 16*Rc*log(Rsigma)^2 ) *
 u1A(A)^(-n1/2) * u2A(A)^(-n2/2) *
 ( pgamma(u1A(A)/sigma.high^2,n1/2) - pgamma(u1A(A)/sigma.low^2,n1/2) ) *
 ( pgamma(u2A(A)/sigma.high^2,n2/2) - pgamma(u2A(A)/sigma.low^2,n2/2) )
}
# integrate
pCSbark <- integrate(fnc2, lower=low, upper=high)$value
pCSbark


################################################################################
# pCbarSk
# Compute pCbarSk = p(Cbar,S|D_1,D_2,I) * p(D_1,D_2|I)
# Mathematica does not allow variable definitions containing a character with an overhead bar so we will use Cbar to represent C-with-a-bar. 
# >>> pC-with-a-bar-Sk becomes pCbarSk
z1 <- n1 * (d1squbar - d1bar^2)
u1H <- sqrt(n1/2) * (high - d1bar)
u1L <- sqrt(n1/2) * (low - d1bar)
z2 <- n2 * (d2squbar - d2bar^2)
u2H <- sqrt(n2/2) * (high - d2bar)
u2L <- sqrt(n2/2) * (low - d2bar)

(u1H/u1L) > 0
(u2H/u2L) > 0
#??????
#sign1 <- -Sign[u1H/ u1L]
#sign2 <- -Sign[u2H/ u2L]

# formula to calculate area below the curve / integral
fnc3 <- function(sigma1)
{
 (2*pi)^(-n/2) * pi / ( 8 * Rc^2 * log(Rsigma) * sqrt(n1*n2) ) * sigma1^(-n+1) *
 exp(-(z1+z2)/(2*sigma1^2)) *
 (errf(ERR=u1H/sigma1) - errf(ERR=u1L/sigma1)) *
 (errf(ERR=u2H/sigma1) - errf(ERR=u2L/sigma1))
}
# integrate
pCbarSk <- integrate(fnc3, lower=sigma.low, upper=sigma.high)$value
pCbarSk


################################################################################
# pCbarSbark
# Compute pCbarSbark = p(Cbar,Sbar|D_1,D_2,I) * p(D_1,D_2|I)
# Mathematica does not allow variable definitions containing a character with an overhead bar so we will use CBar/ Sbar to represent
# C-with-a-bar and S-with-a-bar.
# >>> pC-with-a-barS-with-a-bark becomes pCbarSbark

(u1H/u1L) > 0
(u2H/u2L) > 0
#??????
#sign1 <- -Sign[u1H/ u1L]
#sign2 <- -Sign[u2H/ u2L]

# formula to calculate area below the curve / integral
fnc4.A <- function(sigma1)
{
 sigma1^(-n1) * exp(-z1/(2*sigma1^2)) *
 (errf(ERR=u1H/sigma1) - errf(ERR=u1L/sigma1))
}
# integrate
pCbarSbark.A <- integrate(fnc4.A, lower=sigma.low, upper=sigma.high)$value
pCbarSbark.A

# formula to calculate area below the curve / integral
fnc4.B <- function(sigma2)
{
 sigma2^(-n2) * exp(-z2/(2*sigma2^2)) *
 (errf(ERR=u2H/sigma2) - errf(ERR=u2L/sigma2))
}
# integrate
pCbarSbark.B <- integrate(fnc4.B, lower=sigma.low, upper=sigma.high)$value
pCbarSbark.B

pCbarSbark <- ( (2*pi)^(-n/2) * pi * pCbarSbark.A * pCbarSbark.B ) / ( 8 * Rc^2 * log(Rsigma)^2 * sqrt(n1*n2) )
pCbarSbark


################################################################################
# total probability
# Compute p(D_1,D_2|I) and the normalized probabilities

# tot = pCSk + pCSbark + pCbarSk + pCbarSbark = p(D_1,D_2|I)^2
tot <- pCSk + pCSbark + pCbarSk + pCbarSbark
# pCS = probability that the means &  standard  deviations are the same
pCS <- pCSk / tot
# pCSbar = probability that the means are the same &  standard  deviations are different
pCSbar <- pCSbark / tot
# pCbarS = probability that the means are different &  standard  deviations are the same
pCbarS <- pCbarSk / tot
# pCbarSbar = probability that the means  &  standard  deviations are different
pCbarSbar <- pCbarSbark / tot
# pC = probability that the means are the same independent of whether the standard deviations are the same or different
pC <- pCS + pCSbar
# pCbar = probability that the means are different  independent of whether the standard deviations are the same or different
pCbar <- pCbarS + pCbarSbar
# Odds in favor of different means
oddCbarC <- pCbar / pC
# pS = probability that the standard deviations are the same  independent of whether the means are the same or different
pS <- pCS + pCbarS
# pSbar = probability that the standard deviations are different  independent of whether the means are the same or different
pSbar <- pCSbar + pCbarSbar
# Odds in favor of different standard deviations
oddSbarS <- pSbar / pS
# Odds in favor of different means and/or standard deviations
odddiff <- (1 - pCS) / pCS


################################################################################
# create output tables

# table/ dataframe with input values
desc.df <- data.frame("No." = c(n1,n2,n1+n2),
                       "Standard Deviation" = c(d1std,d2std,dstd),
                       "Variance" = c(d1std^2,d2std^2,dstd^2),
                       "Mean" = c(d1bar,d2bar,dbar),
                       "Data set" = c(snames,"combined"),
                       check.names = FALSE)

# table/ dataframe with prior values/ information
prior.df <- data.frame("Numerical Example" = c("Prior Mean lower bound",
                                               "Prior Mean upper bound",
                                               "Prior Standard Deviation lower bound",
                                               "Prior Standard Deviation upper bound",
                                               "Number of steps for plotting p(delta | D_1, D_2, I)",
                                               "Number of steps for plotting p(r | D_1, D_2, I)"),
                       "Value" = c(low,high,sigma.low,sigma.high,ndelta,nr),
                       check.names = FALSE)

# table/ dataframe with resulting probabilities
posterior.df <- data.frame("Hypothesis" = c(
                                           "C, S       = Same Mean, same Standard Deviation",
                                           "Cbar, S    = Different Means, same Standard Deviation",
                                           "C, Sbar    = Same Mean, different Standard Deviations",
                                           "Cbar, Sbar = Different Means, different Standard Deviations",
                                           "C          = The Means are the same",
                                           "Cbar       = The Means are different",                            
                                           "S          = The Standard Deviations are the same",
                                           "Sbar       = The Standard Deviations are different",                            
                                           "C, S       = Same Means and Standard Deviations",
                                           "Cbar, Sbar = One or Both are different"),                           
                          "Probability" = c(pCS,pCbarS,pCSbar,pCbarSbar,pC,pCbar,pS,pSbar,pCS,1-pCS),
            check.names = FALSE)

# table/ dataframe with odds ratio results
oddsratio.res.df <- data.frame("Hypothesis" = c("The odds ratio in favour of a difference (means)",
                                                "The odds ratio in favour of a difference (standard deviations)",
                                                "The odds ratio in favour of a difference (means and standard deviations)",
                                                "The odds ratio in favour of the same (means)",
                                                "The odds ratio in favour of the same (standard deviations)",
                                                "The odds ratio in favour of the same (means and standard deviations)"),
                               "Odds Ratio" = c(oddCbarC,oddSbarS,odddiff,
                                                1/oddCbarC,1/oddSbarS,1/odddiff),
                               check.names = FALSE)
# print tables
if(print.res)
{
 cat("\n#####################################")
 cat("\n\n  Output on the difference in means")
 cat("\n\n\n")
 print(desc.df, right=FALSE, row.names=FALSE)
 cat("\n")
 print(prior.df, right=FALSE, row.names=FALSE)
 cat("\n")
 print(posterior.df, digits=dig, right=FALSE, row.names=FALSE)
 cat("\n")
 print(oddsratio.res.df, digits=dig+1, right=FALSE, row.names=FALSE)
 cat("\n#####################################")
 cat("\n")
}

 res <- list(snames=snames,
             d1=d1,
             d2=d2,
             n1=n1,
             n2=n2,
             n=n,
             low=low,
             high=high,
             sigma.low=sigma.low,
             sigma.high=sigma.high,
             d1bar=d1bar,
             d2bar=d2bar,
             dbar=dbar,
             d1std=d1std,
             d2std=d2std,
             dstd=dstd,
             d1squbar=d1squbar,
             d2squbar=d2squbar,
             dsqubar=dsqubar,
             Rc=Rc,
             Rsigma=Rsigma,
             nr=nr,
             ndelta=ndelta,
             descriptive=desc.df,
             prior=prior.df,
             posterior=posterior.df,
             OR=oddsratio.res.df)
 attr(res,"invtyp") <- invtyp

return(res)
}
# end Bretthorst adapted from Mathematica code by Phil Gregory
################################################################################
#call
#dim.res <- DiM.pg(invtyp="pg", inputvalues, print.res=TRUE)
#DiM.pg(invtyp="ums", inputvalues, print.res=TRUE)
#
#plot.DiM(DiM.res=dim.res, filling=TRUE)


################################################################################
plot.DiM <- function(DiM.res, cols=NULL, fac=1.05, lwd=1.9, alphav=0.2,
                     filling=FALSE, by1=TRUE, colb="red",
                     dig=2, scaleL=30, scaleH=4)
{

# error function
errf <- function(ERR) return( 2*pnorm(ERR*sqrt(2))-1 )

 if(is.null(cols)) cols <- c("red","darkgreen","magenta","orange","lightgreen","violet")
 cols.rgb <- col2rgb(cols[4:6])/255
 cols
 cols.rgb
  
##map values
 low <- DiM.res[["low"]]
 high <- DiM.res[["high"]]
 ndelta <- DiM.res[["ndelta"]]
 n1 <- DiM.res[["n1"]]
 n2 <- DiM.res[["n2"]]
 n <- DiM.res[["n"]]
 d1bar <- DiM.res[["d1bar"]]
 d2bar <- DiM.res[["d2bar"]]
 Rc <- DiM.res[["Rc"]]
 Rsigma <- DiM.res[["Rsigma"]]
 dsqubar <- DiM.res[["dsqubar"]]  
 dbar <- DiM.res[["dbar"]]
 sigma.high <- DiM.res[["sigma.high"]]
 sigma.low <- DiM.res[["sigma.low"]]
 d1squbar <- DiM.res[["d1squbar"]]
 d2squbar <- DiM.res[["d2squbar"]]
 posterior <- DiM.res[["posterior"]]
 d1std <- DiM.res[["d1std"]]
 d2std <- DiM.res[["d2std"]]
 nr <- DiM.res[["nr"]]

  
################################################################################
# p(delta|S,D_1,D_2,I) problem
# plot The difference in the means (delta)

# This section computes and plots the posterior probability density
# for the difference in the two means assuming the standard deviations are the same.

# create sequence for values to calculate the posterior probability density(ies)
delta.low <- low - high
delta.high <- high - low
delta.delta <- (delta.high - delta.low) / (ndelta - 1)
delta.sek <- seq(delta.low,delta.high,delta.delta)

beta.low <- 2 * low
beta.high <- 2 * high

delta.N <- (n1 - n2) / n

b <- (n1 * d1bar - n2 * d2bar) / (2 * n)

V <- function(delta, BETA) {
  (n / 2) * (
    dsqubar - (2 * delta * b) - (BETA * dbar) + (BETA ^ 2 / 4) + (delta ^ 2 /
                                                                    4) + (delta * BETA * delta.N / 2)
  )
}

fndbV <- function(BETA) {
  gamma(n / 2) / (8 * Rc ^ 2 * log(Rsigma)) * V(delta = delta, BETA = BETA) ^
    (-n / 2) *
    (pgamma(V(delta = delta, BETA = BETA) / sigma.high ^ 2, n / 2) - pgamma(V(delta =
                                                                                delta, BETA = BETA) / sigma.low ^ 2, n / 2))
}

# integrate for each element of delta.sek
# do it later with sapply
pdel.SD1D2I <- vector(mode="numeric",length=length(delta.sek))
for(i in 1:length(delta.sek))
{
# print(i)
 delta <- delta.sek[i]
 pdel.SD1D2I[i] <- integrate(fndbV, lower=beta.low, upper=beta.high)$value
}
#pdel.SD1D2I

pdelA <- pdel.SD1D2I / (delta.delta * sum(pdel.SD1D2I))
#pdelA
# create dataframe with input and results
pdelta.df <- data.frame(delta.sek, pdel.SD1D2I, pdelA)
#pdelta.df

# not necessary but in Mathematica
# xdel <- delta.sek #vector
# xdA <- pdelta.df[,c("delta.sek","pdelA")] #table/ data.frame with 2 cols

# plot p(delta|S,D_1,D_2,I)
 x.lims <- c(min(delta.sek),max(delta.sek))
 y.lims <- c(min(pdelA),max(pdelA)*fac)
 if(by1) par(ask=TRUE)
 par(mar=c(5,6,5,5))
 par("cex.axis" = 0.8)
 plot(
   x.lims[2],
   y.lims[2],
   col = "white",
   main = "",
   xlab = "",
   ylab = expression(paste("p(", delta, " | S, D" ["1"] , ", D" ["2"] , ", I)", sep =
                             "")),
   bty = "l",
   pre.plot = grid(),
   xlim = x.lims,
   ylim = y.lims
 )
 #rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col="grey95", border=NA)
 #grid(col="white", lwd=1, lty=1)
 #grid()
 if (filling)
   polygon(
     c(delta.sek[1], delta.sek, delta.sek[length(delta.sek)]),
     c(0, pdelA, 0),
     col = rgb(cols.rgb[1, 1], cols.rgb[2, 1], cols.rgb[3, 1], 0.5),
     border = NA
   )
 
 points(delta.sek, pdelA, col=cols[1], type="l", lty=1, lwd=lwd)
 axis(2)
 axis(1)
 mtext(expression(paste(delta," Mean Difference",sep="")), 3, line=2, cex=1.5)
 mtext("assuming the standard deviations are the same", 3, line=0.8)
 mtext(expression(paste(delta," = ", mu["1"]," ? ",mu["2"])), 1, line=3, cex=1.5, col=colb)


################################################################################
# p(delta|Sbar,D_1,D_2,I) problem
# This section computes and plots the posterior probability density
# for the difference in the two means assuming the standard deviations are different.

#redundant, not necessary
delta.low <- low - high
delta.high <- high - low
delta.delta <- (delta.high - delta.low) / (ndelta - 1)
beta.low <- 2 * low
beta.high <- 2 * high

w1 <- function(delta, BETA)
{
 (n1/2) * (d1squbar - (BETA + delta) * d1bar + ((BETA + delta)^2/4) )
}

w2 <- function(delta, BETA)
{
 (n2/2) * (d2squbar - (BETA - delta) * d2bar + ((BETA - delta)^2/4) )
} 

fndbVbar <- function(BETA)
{
 gamma(n1/2) * gamma(n2/2) / (16 * Rc^2 * log(Rsigma)^2) *
 w1(delta=delta,BETA=BETA)^(-n1/2) * w2(delta=delta,BETA=BETA)^(-n2/2) *
 ( pgamma(w1(delta=delta,BETA=BETA)/sigma.high^2,n1/2) - pgamma(w1(delta=delta,BETA=BETA)/sigma.low^2,n1/2) ) *
 ( pgamma(w2(delta=delta,BETA=BETA)/sigma.high^2,n2/2) - pgamma(w2(delta=delta,BETA=BETA)/sigma.low^2,n2/2) )
}

# integrate for each element of delta.sek
# do it later with sapply
pdel.SbarD1D2I <- vector(mode="numeric",length=length(delta.sek))
for(i in 1:length(delta.sek))
{
# print(i)
 delta <- delta.sek[i]
 pdel.SbarD1D2I[i] <- integrate(fndbVbar, lower=beta.low, upper=beta.high)$value
}
#pdel.SbarD1D2I

pdelB <- pdel.SbarD1D2I / (delta.delta * sum(pdel.SbarD1D2I))
#pdelB
# create dataframe with input and results
pdelta.df <- data.frame(pdelta.df, pdel.SbarD1D2I, pdelB)
#pdelta.df

# not necessary but in Mathematica
# xdel <- delta.sek #vector
# xdB <- pdelta.df[,c("delta.sek","pdelB")] #table/ data.frame with 2 cols

# plot p(delta|Sbar,D_1,D_2,I)
 x.lims <- c(min(delta.sek),max(delta.sek))
 y.lims <- c(min(pdelB),max(pdelB)*fac)
 if(by1) par(ask=TRUE)
 par(mar=c(5,6,5,5))
 par("cex.axis" = 0.8)
 plot(x.lims[2], y.lims[2],
     col="white",
     main="",
     xlab="",
     ylab=expression(paste("p(",delta," | ",bar(S),", D" ["1"] ,", D" ["2"] ,", I)",sep="")),
     pre.plot=grid(),
     bty="l",
     xlim=x.lims, ylim=y.lims)
 #rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col="grey95", border=NA)
 #grid(col="white", lwd=1, lty=1)
 if(filling) polygon(c(delta.sek[1],delta.sek,delta.sek[length(delta.sek)]), c(0,pdelB,0), col=rgb(cols.rgb[1,2],cols.rgb[2,2],cols.rgb[3,2],0.5), border=NA)
 points(delta.sek, pdelB, col=cols[2], type="l", lty=1, lwd=lwd)
 axis(2)
 axis(1)
 mtext(expression(paste(delta," Mean Difference",sep="")), 3, line=2, cex=1.5)
 mtext("assuming the standard deviations are different", 3, line=0.8)
 mtext(expression(paste(delta," = ", mu["1"]," ? ",mu["2"])), 1, line=3, cex=1.5, col=colb)


################################################################################
# p(delta|D_1,D_2,I) problem
# This section computes and plots the posterior probability density
# for the difference in the two means independent of whether the
# standard deviations are the same. It is a weighted average of
# p(delta|S,D_1,D_2,I) and p(delta|Sbar,D_1,D_2,I) where the weights
# are pS and pSbar, the probability that the standard deviations are
# the same or different, respectively.

# The code also computes the 95% credible region and the 99.5% credible region.
# The latter is used to choose suitable boundaries for plotting purposes so
# only the region with significant probability is displayed.

# to remember - s.a.:
# p(delta|S,D_1,D_2,I)
# xdel <- delta.sek                         #vector
# xda <- pdelta.df[,c("delta.sek","pdelA")] #table/ data.frame with 2 cols
#
# p(delta|Sbar,D_1,D_2,I)
# xdel <- delta.sek                         #vector
# xdB <- pdelta.df[,c("delta.sek","pdelB")] #table/ data.frame with 2 cols

# taken from results
pS <- posterior[7,2]
pSbar <- posterior[8,2]

pdelave <- (pS * pdelA + pSbar * pdelB)
#pdelave
pdelave <- pdelave / sum(pdelave)
#pdelave
rxd <- data.frame(pdelave,delta.sek) # delta.sek = xdel
delta.peak <- rxd[which(pdelave == sort(pdelave,dec=TRUE)[1]),"delta.sek"] #MA: sort[rxd][[ndelta,2]]
delta.mean <- pdelave %*% delta.sek # MA: pdelave.xdel
sortrxd <- rxd[with(rxd, order(pdelave,delta.sek)),] #MA: sort(rxd)

delta.min <- delta.high
delta.max <- delta.low
AREA.del <- 0
ndelta.i <- ndelta

#rxd
delta.peak
delta.mean
delta.min
delta.max
AREA.del
ndelta.i

i <- 0
while(AREA.del <= 0.95)
{
#  print(i)
#  print(AREA.del)
  AREA.del <- AREA.del + sortrxd[ndelta.i,"pdelave"] #MA: sortrxd[[ndelta.i,1]]
#  print(AREA.del)
  delta.crit <- sortrxd[ndelta.i,"delta.sek"]
#  print(delta.crit)
  if(delta.crit > delta.max) delta.max <- delta.crit #MA:sortrxd[[ndelta.i,2]] ,Null
  if(delta.crit < delta.min) delta.min <- delta.crit #MA:sortrxd[[ndelta.i,2]] ,Null
#  print(delta.min)
#  print(delta.max)
#  cat(paste("\n\n"))
  ndelta.i <- ndelta.i - 1
  i <- i +1
}

delta.bmin <- delta.min
delta.bmax <- delta.max

j <- 0
while( (AREA.del <= 0.995) && (ndelta.i > 0) )
{
#  print(j)
#  print(AREA.del)
  AREA.del <- AREA.del + sortrxd[ndelta.i,"pdelave"] #MA: sortrxd[[ndelta.i,1]]
#  print(AREA.del)
  delta.bcrit <- sortrxd[ndelta.i,"delta.sek"]
#  print(delta.bcrit)
  if(delta.bcrit > delta.bmax) delta.bmax <- delta.bcrit #MA:sortrxd[[ndelta.i,2]] ,Null
  if(delta.bcrit < delta.bmin) delta.bmin <- delta.bcrit #MA:sortrxd[[ndelta.i,1]] ,Null
#  print(delta.bmin)
#  print(delta.bmax)
#  cat(paste("\n\n"))
  ndelta.i <- ndelta.i - 1
  j <- j +1
}

delta.bmin
delta.bmax

xdA <- pdelta.df[,c("delta.sek","pdelA")]
xdB <- pdelta.df[,c("delta.sek","pdelB")]
#xdA
#xdB

#extract and sort according to delta.sek
#extract: sortrxd[(ndelta.i+1):ndelta,]
temp <- sortrxd[(ndelta.i+1):ndelta,]
#temp
pdelwta <- temp[with(temp, order(delta.sek,pdelave)),]
#pdelwta
#MA: extract + sort according to delta.sek #pdelwta <- sort[RotateLeft[Take[sortrxd,{nδi+1,nδ}],{0,1}]]
#tr <- t(pdelwta)
#pdelwt <- t( [{tr[[1]],tr[[2]] / delta.delta}
#pdelx <- tr[[1]]
pdelwt <- data.frame(pdelwta, pdelwt = pdelwta[,"pdelave"]/delta.delta)
#pdelwt
pdelx <- pdelwt[,"delta.sek"]
#pdelx

#MA: search intersection (R: a %in% b) of xdA (xdB) versus pdelx and use
#that for delta.sek & pdelA (pdelB) to create xdAsub (xdBsub)
xdAsub <- xdA[xdA[,"delta.sek"] %in% pdelx,]
xdBsub <- xdB[xdB[,"delta.sek"] %in% pdelx,]
#xdAsub
#xdBsub

#see whether values are ok
ndelta
ndelta.i
AREA.del
delta.peak
delta.mean
delta.min
delta.max

# Here we plot p(delta|v,D_1,D_2,I), p(delta|vbar,D_1,D_2,I) and
# the weighted average of the two. See plot legend for details.
# The peak delta, average delta and upper and lower boundaries
# of the 95% credible region are displayed along the top of the plot.
pdelwt.4plot <- pdelwt[,c("delta.sek","pdelwt")]
#pdelwt.4plot

# plot p(delta|v,D_1,D_2,I), p(delta|vbar,D_1,D_2,I)
y.lims <- c(min(c(pdelwt.4plot[,"pdelwt"],xdAsub[,"pdelA"],xdBsub[,"pdelB"])),
            max(c(pdelwt.4plot[,"pdelwt"],xdAsub[,"pdelA"],xdBsub[,"pdelB"]))*fac)
x.lims <- c(min(pdelwt.4plot[,"delta.sek"]),max(pdelwt.4plot[,"delta.sek"]))
#if(ylim.max * fac > 1) ylim.max <- 1 else ylim.max <- ylim.max * fac
 if(by1) par(ask=TRUE)
 par(mar=c(5,6,5,5))
 par(oma=c(2,1,1,1))
 par("cex.axis" = 0.8)
 plot(x.lims[2], y.lims[2], col="white",
     main="",
     xlab="",
     ylab="Probability Density",
     pre.plot=grid(),
     bty="l", xlim=x.lims, ylim=y.lims)
 #rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col="grey95", border=NA)
 #grid(col="white", lwd=1, lty=1)
 if(filling) polygon(c(xdAsub[1,"delta.sek"],xdAsub[,"delta.sek"],xdAsub[dim(xdAsub)[1],"delta.sek"]), c(0,xdAsub[,"pdelA"],0), col=rgb(cols.rgb[1,1],cols.rgb[2,1],cols.rgb[3,1],alphav), border=NA)
 points(xdAsub, type="l", lty=2, col=cols[1], lwd=lwd)
 if(filling) polygon(c(xdBsub[1,"delta.sek"],xdBsub[,"delta.sek"],xdBsub[dim(xdBsub)[1],"delta.sek"]), c(0,xdBsub[,"pdelB"],0), col=rgb(cols.rgb[1,2],cols.rgb[2,2],cols.rgb[3,2],alphav), border=NA)
 points(xdBsub, type="l", lty=3, col=cols[2], lwd=lwd)
 if(filling) polygon(pdelwt.4plot[,"delta.sek"], pdelwt.4plot[,"pdelwt"], col=rgb(cols.rgb[1,3],cols.rgb[2,3],cols.rgb[3,3],alphav), border=NA)
 points(pdelwt.4plot[,"delta.sek"], pdelwt.4plot[,"pdelwt"], type="l", lty=1, col=cols[3], lwd=lwd)
 axis(2)
 axis(1)
 mtext(expression(paste(delta," Mean Difference",sep="")), 3, line=3, cex=1.5)
 mtext("independent of the standard deviations being different/ the same", 3, line=1.5)
 mtext(expression(paste(delta," = ", mu["1"]," ? ",mu["2"])), 1, line=3, cex=1.5, col=colb)
# add a nice legend with information 
 mtext(eval(substitute(expression(paste("mode (peak) ",delta," = ",delta.peak, "  |  mean ",delta," = ",delta.mean)),
       list(delta.peak=round(delta.peak,dig), delta.mean=round(delta.mean,dig) ))) , 4, line=1, cex=0.9, col="black")
 mtext(eval(substitute(expression(paste("95% lower ",delta," = ", delta.min, "  |  95% upper ",delta," = ", delta.max)),
       list(delta.min=round(delta.min,dig), delta.max=round(delta.max,dig) ))) , 4, line=2, cex=0.9, col="black")
 par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
 plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n")
# add a nice legend with information
 legend("bottom", legend=c(expression(paste("p(",delta," | S, D" ["1"] ,", D" ["2"] ,", I)")),
                           expression(paste("p(",delta," | ",bar(S),", D" ["1"] ,", D" ["2"] ,", I)")),
                           expression(paste("p(",delta," | D" ["1"] ,", D" ["2"] ,", I)"))),
                           xpd=TRUE, horiz=TRUE, inset=c(0,0), y.intersp=2.4,
                           col=cols[1:3], lty=c(2,3,1), lwd=lwd, bty="n", cex=0.9)                         


################################################################################
# p(r|C,D_1,D_2,I) problem
# This section computes and plots the posterior probability density
# for the ratio of the standard deviations assuming the means are the same.
# The upper and lower bounds on the plot can be varied by the user
# by setting scaleL and scaleH factors.

#use as input values to allow to change them:
#scaleL <- 30
#scaleH <- 4

r.low <- d1std / (scaleL * d2std)
r.high <- scaleH * d1std / d2std
r.delta <- (r.high - r.low) / (nr - 1)
AL <- low
AH <- high

u <- function(r) ( (n1 * d1squbar) / r^2 ) + n2 * d2squbar
v <- function(r) ( (n1 * d1bar) / r^2 ) + n2 * d2bar
w <- function(r) n1 / r^2 + n2
z <- function(r) u(r=r) - v(r=r)^2 / w(r=r)
XH <- function(sigma) sqrt(w(r=r) / (2 * sigma^2)) * (AH - v(r=r)/w(r=r))
XL <- function(sigma) sqrt(w(r=r) / (2 * sigma^2)) * (AL - v(r=r)/w(r=r))

fnrS <- function(sigma)
{
 (2*pi)^(-n/2) * sqrt(pi/ (8*w(r=r))) / (Rc * log(Rsigma)^2) *
 r^(-(n1+1)) * sigma^(-n) *
 exp(-z(r=r)/ (2*sigma^2)) * ( errf(ERR=XH(sigma=sigma)) - errf(ERR=XL(sigma=sigma)) )
}

r.sek <- seq(r.low,r.high,r.delta)
#r.sek #MA: xr <- r.sek
# integrate for each element of delta.sek
# do it later with sapply
pr.CD1D2I <- vector(mode="numeric",length=length(r.sek))
for(i in 1:length(r.sek))
{
# print(i)
 r <- r.sek[i]
 pr.CD1D2I[i] <- integrate(fnrS, lower=sigma.low, upper=sigma.high)$value
}
#pr.CD1D2I
#note by LG:
#first value different to MA:
#MA: 3.15108e-230
#R: 3.206374e-230
#and consequently prB is different, too:
#MA: 1.0686e-210
#R: 1.087353e-210
#all other values are 100% identical...
prA <- pr.CD1D2I / (r.delta * sum(pr.CD1D2I))
#prA
# create dataframe with input and results
pr.df <- data.frame(r.sek, pr.CD1D2I, prA)
#pr.df

# not necessary from Mathematica
# xdel <- 
# xda <- t()

# plot p(r|C,D_1,D_2,I)
 x.lims <- c(min(r.sek),max(r.sek))
 y.lims <- c(min(prA),max(prA)*fac)
 if(by1) par(ask=TRUE)
 par(mar=c(5,6,5,5))
 plot(x.lims[2],y.lims[2], type="l", lty=1, col="white",
     main="",
     xlab="",
     ylab=expression(paste("p(r | C, D" ["1"] ,", D" ["2"] ,", I)",sep="")),
     pre.plot=grid(),
     bty="l", xlim=x.lims, ylim=y.lims)
 #rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col="grey95", border=NA)
 #grid(col="white", lwd=1, lty=1)
 if(filling) polygon(c(r.sek[1],r.sek,r.sek[length(r.sek)]), c(0,prA,0), col=rgb(cols.rgb[1,1],cols.rgb[2,1],cols.rgb[3,1],alphav), border=NA)
 points(r.sek, prA, col=cols[1], type="l", lty=1, lwd=lwd)
 axis(2)
 axis(1)
 mtext(expression(paste(sigma," Standard Deviation Ratio (r)",sep="")), 3, line=2, cex=1.5)
 mtext("assuming the means are the same", 3, line=0.8)
 mtext(expression(paste("r = ",sigma[1],"/",sigma[2])), 1, line=3, cex=1.5, col=colb)


################################################################################
# p(r|Cbar,D_1,D_2,I) problem
# This section computes and plots the posterior probability density
# for the ratio of the standard deviations assuming the means are the different. 

# redundant, not necessary
#scaleL <- 30
#scaleH <- 4
r.low <- d1std / (scaleL * d2std)
r.high <- scaleH * d1std / d2std
r.delta <- (r.high - r.low) / (nr - 1)
AL <- low
AH <- high

z1 <- n1 * (d1squbar - d1bar^2)
z2 <- n2 * (d2squbar - d2bar^2)

X1H <- function(r, sigma) sqrt(n1 / (2 * r^2 * sigma^2)) * (AH - d1bar)
X1L <- function(r, sigma) sqrt(n1 / (2 * r^2 * sigma^2)) * (AL - d1bar)
X2H <- function(sigma) sqrt(n2 / (2 * sigma^2)) * (AH - d2bar)
X2L <- function(sigma) sqrt(n2 / (2 * sigma^2)) * (AL - d2bar)

fnrSb <- function(sigma)
{
 ( (2*pi)^(-n/2) * pi ) / ( 4 * Rc^2 * log(Rsigma)^2 * sqrt(n1*n2) ) *
 r^(-n1) * sigma^(1-n) *
 exp( -(z1 / (2 * r^2 * sigma^2)) - (z2 / (2*sigma^2)) ) *
 ( errf(ERR=X1H(r=r, sigma=sigma)) - errf(ERR=X1L(r=r, sigma=sigma)) ) *
 ( errf(ERR=X2H(sigma=sigma)) - errf(ERR=X2L(sigma=sigma)) )
}

# integrate for each element of delta.sek
# do it later with sapply
pr.CbarD1D2I <- vector(mode="numeric",length=length(r.sek))
for(i in 1:length(r.sek))
{
# print(i)
 r <- r.sek[i]
 pr.CbarD1D2I[i] <- integrate(fnrSb, lower=sigma.low, upper=sigma.high)$value
}
#pr.CbarD1D2I
#note by LG:
#first value different to MA:
#MA: 4.06446e-230
#R: 4.135770e-230
#and consequently prB is different, too:
#MA: 4.31101e-211
#R: 4.386639e-211
#all other values are 100% identical...
prB <- pr.CbarD1D2I / (r.delta * sum(pr.CbarD1D2I))
#prB
# create dataframe with input and results
pr.df <- data.frame(pr.df, pr.CbarD1D2I, prB)
#pr.df

# not necessary from Mathematica
# xdel <- 
# xda <- t()

# plot p(r|Cbar,D_1,D_2,I)
 x.lims <- c(min(r.sek),max(r.sek))
 y.lims <- c(min(prB),max(prB)*fac)
 if(by1) par(ask=TRUE)
 par(mar=c(5,6,5,5))
 plot(x.lims[2], y.lims[2], col="white",
     main="",
     xlab="",
     ylab=expression(paste("p(r | ",bar(C),", D" ["1"] ,", D" ["2"] ,", I)",sep="")),
     pre.plot=grid(),
     bty="l", xlim=x.lims, ylim=y.lims)
 #rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col="grey95", border=NA)
 #grid(col="white", lwd=1, lty=1)
 if(filling) polygon(c(r.sek[1],r.sek,r.sek[length(r.sek)]), c(0,prB,0), col=rgb(cols.rgb[1,2],cols.rgb[2,2],cols.rgb[3,2],alphav), border=NA)
 points(r.sek, prB, col=cols[2], type="l", lty=1, lwd=lwd)
 axis(2)
 axis(1)      
 mtext(expression(paste(sigma," Standard Deviation Ratio (r)",sep="")), 3, line=2, cex=1.5)
 mtext("assuming the means are different", 3, line=0.8)
 mtext(expression(paste("r = ",sigma[1],"/",sigma[2])), 1, line=3, cex=1.5, col=colb)


################################################################################
# p(r|D_1,D_2,I)
# This section computes and plots the posterior probability density
# for the ratio of the standard deviations independent of whether
# the means are the same.
# It is a weighted average of p(r|C,D1,D2,I) and p(r|Cbar,D1,D2,I)
# where the weights are pC and pCbar, the probability that the means
# are the same or different, respectively.
# The code also computes the 95% credible region and the 99.9% credible region.
# The latter is used to choose suitable boundaries for plotting purposes
# so only the region with significant probability is displayed.

#note by LG:
#first (smalles values) are slightly different due to following errors from
#above -> prA and prB, not pC and pCbar -> probably due to precision on the
#not-log scale (ie. not using logs to integrate)

pC <- posterior[5,2]
pCbar <- posterior[6,2]
pC
pCbar

prave <- (pC * prA + pCbar * prB)
#prave
prave <- prave / sum(prave)
#prave

rxr <- data.frame(prave,r.sek) # r.sek = xr
r.peak <- rxr[which(prave == sort(prave,dec=TRUE)[1]),"r.sek"] #MA: sort[rxr][[nr,2]]
r.mean <- prave %*% r.sek # MA: pdelave.xr
sortrxr <- rxr[with(rxr, order(prave,delta.sek)),] #MA: sort(rxr)

r.min <- r.high
r.max <- r.low
AREA.r <- 0
nr.i <- nr

#rxr
r.peak
r.mean
r.min
r.max
AREA.r
nr.i

i <- 0
while(AREA.r <= 0.95)
{
#  print(i)
#  print(AREA.r)
  AREA.r <- AREA.r + sortrxr[nr.i,"prave"] #MA: sortrxr[[nr.i,1]]
#  print(AREA.r)
  r.crit <- sortrxr[nr.i,"r.sek"]
#  print(r.crit)
  if(r.crit > r.max) r.max <- r.crit #MA:sortrxr[[nr.i,2]] ,Null
  if(r.crit < r.min) r.min <- r.crit #MA:sortrxr[[nr.i,2]] ,Null
#  print(r.min)
#  print(r.max)
#  cat(paste("\n\n"))
  nr.i <- nr.i - 1
  i <- i +1
}  
r.min
r.max

r.bmin <- r.min
r.bmax <- r.max

j <- 0
while( (AREA.r <= 0.999) && (nr.i > 0) )
{
#  print(j)
#  print(AREA.r)
  AREA.r <- AREA.r + sortrxr[nr.i,"prave"] #MA: sortrxr[[nr.i,1]]
#  print(AREA.r)
  r.bcrit <- sortrxr[nr.i,"r.sek"]
#  print(r.bcrit)
  if( r.bcrit > r.bmax) r.bmax <- r.bcrit #MA:sortrxr[[nr.i,2]] ,Null
  if( r.bcrit < r.bmin) r.bmin <- r.bcrit #MA:sortrxr[[nr.i,1]] ,Null
#  print(r.bmin)
#  print(r.bmax)
#  cat(paste("\n\n"))
  nr.i <- nr.i - 1
  j <- j +1
}
r.bmin
r.bmax

xrA <- pr.df[,c("r.sek","prA")]
xrB <- pr.df[,c("r.sek","prB")]
#xrA
#xrB

#extract and sort according to delta.sek
#extract: sortrxd[(ndelta.i+1):ndelta,]
temp <- sortrxr[(nr.i+1):nr,]
#temp
prwta <- temp[with(temp, order(r.sek,prave)),]
#prwta
#MA: extract + sort according to r.sek #prwta <- sort[RotateLeft[Take[sortrxr,{nδi+1,nδ}],{0,1}]]
#tr <- t(prwta)
#prwt <- t( [{tr[[1]],tr[[2]] / r.delta}
#prx <- tr[[1]]
prwt <- data.frame(prwta, prwt = prwta[,"prave"]/r.delta)
#prwt
prx <- prwt[,"r.sek"]
#prx

#MA: search intersection (R: a %in% b) of xrA (xrB) versus prx and use
#that for r.sek & prA (prB) to create xrAsub (xrBsub)
xrAsub <- xrA[xrA[,"r.sek"] %in% prx,]
xrBsub <- xrB[xrB[,"r.sek"] %in% prx,]
#xrAsub
#xrBsub

#see whether values are ok
nr
nr.i
AREA.r
r.peak
r.mean
r.min
r.max

# Here we plot  p(r|s,D1,D2,I), p(r|sbar,D1,D2,I) and the weighted average
# of the two. See plot legend for details.
# The peak r, average r and upper and lower boundaries of the 95% credible
# region are displayed along the top of the plot.
prwt.4plot <- prwt[,c("r.sek","prwt")]
#prwt.4plot

# plot p(r|C,D_1,D_2,I), p(r|Cbar,D_1,D_2,I)
y.lims <- c(min(c(prwt.4plot[,"prwt"],xrAsub[,"prA"],xrBsub[,"prB"])),
            max(c(prwt.4plot[,"prwt"],xrAsub[,"prA"],xrBsub[,"prB"]))*fac)
x.lims <- c(min(prwt.4plot[,"r.sek"]), max(c(prwt.4plot[,"r.sek"])))
#if(ylim.max * fac > 1) ylim.max <- 1 else ylim.max <- ylim.max * fac
 if(by1) par(ask=TRUE)
 par(mar=c(5,6,5,5))
 par(oma=c(2,1,1,1))
 par("cex.axis" = 0.8)
 plot(x.lims[2], y.lims[2], col="white",
      main="",
      xlab="",
      ylab="Probability Density",
      pre.plot=grid(),
      bty="l",
      xlim=x.lims, ylim=y.lims)
 #rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col="grey95", border=NA)
 #grid(col="white", lwd=1, lty=1)
 if(filling) polygon(c(xrAsub[1,"r.sek"],xrAsub[,"r.sek"],xrAsub[dim(xrAsub)[1],"r.sek"]), c(0,xrAsub[,"prA"],0), col=rgb(cols.rgb[1,1],cols.rgb[2,1],cols.rgb[3,1],alphav), border=NA)
 points(xrAsub, type="l", lty=2, col=cols[1], lwd=lwd)
 if(filling) polygon(c(xrBsub[1,"r.sek"],xrBsub[,"r.sek"],xrBsub[dim(xrBsub)[1],"r.sek"]), c(0,xrBsub[,"prB"],0), col=rgb(cols.rgb[1,2],cols.rgb[2,2],cols.rgb[3,2],alphav), border=NA)
 points(xrBsub, type="l", lty=3, col=cols[2], lwd=lwd)
 if(filling) polygon(c(prwt.4plot[1,"r.sek"],prwt.4plot[,"r.sek"],prwt.4plot[dim(prwt.4plot)[1],"r.sek"]), c(0,prwt.4plot[,"prwt"],0), col=rgb(cols.rgb[1,3],cols.rgb[2,3],cols.rgb[3,3],alphav), border=NA)
 points(prwt.4plot[,"r.sek"], prwt.4plot[,"prwt"], type="l", lty=1, col=cols[3], lwd=lwd)
 axis(2)
 axis(1)
 mtext(expression(paste("r Standard Deviation Ratio",sep="")), 3, line=3, cex=1.5)
 mtext("independent of the means being different/ the same", 3, line=1.5)
 mtext(expression(paste("r = ",sigma[1],"/",sigma[2])), 1, line=3, cex=1.5, col=colb)

# add a nice legend with information 
 mtext(eval(substitute(expression(paste("mode (peak) r = ",r.peak, "  |  mean r = ",r.mean)),
       list(r.peak=round(r.peak,dig), r.mean=round(r.mean,dig) ))) , 4, line=1, cex=0.9, col="black")
 mtext(eval(substitute(expression(paste("95% lower r = ", r.min, "  |  95% upper r = ", r.max)),
       list(r.min=round(r.min,dig), r.max=round(r.max,dig) ))) , 4, line=2, cex=0.9, col="black")
 par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
 plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
# add a nice legend with information
# add a nice legend to explain curves and associated probability densities
 legend("bottom", legend=c(expression(paste("p(r | C, D" ["1"] ,", D" ["2"] ,", I)")),
                           expression(paste("p(r | ",bar(C),", D" ["1"] ,", D" ["2"] ,", I)")),
                           expression(paste("p(r | D" ["1"] ,", D" ["2"] ,", I)"))),
                           xpd=TRUE, horiz=TRUE, inset=c(0,0), y.intersp=2.4,
                           col=cols[1:3], lty=c(2,3,1), lwd=lwd, bty="n", cex=0.9)
 
}


################################################################################
# print results / tables
PGprint <- function(res, dig=4)
{
 
 cat("\n ####################################################################################\n")
 cat("\n  Output on the difference in means")
 cat("\n  -----------------------------------------------------\n")
 cat("  Paper by G.L. Bretthorst (1993)\n")
 cat("  R Implementation of a Mathematica script (P. Gregory)\n")
 cat("\n\n")
 print(res$descriptive, right=FALSE, row.names=FALSE)
 cat("\n")
 print(res$prior, right=FALSE, row.names=FALSE)
 cat("\n")
 print(res$posterior, digits=dig, right=FALSE, row.names=FALSE)
 cat("\n")
 print(res$OR, digits=dig+1, right=FALSE, row.names=FALSE)
 cat("\n ####################################################################################\n")
 cat("\n")

}
#call:
#PGprint(dim.res)
 
 
################################################################################
# END OF 'ON THE DIFFERENCE OF MEANS'
################################################################################
