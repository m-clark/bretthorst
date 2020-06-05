################################################################################
# On the difference in means
# paper: G.L. Bretthorst "On the difference of means" (1993)
# http://bayes.wustl.edu/glb/diff.pdf
#
# R code based on Mathematica code by Urban Studer (90's, Zürich/ CH)
# R code by Leo Gürtler
# first: 12-06-05
# latest: 21-06-06, 20-04-17

################################################################################
# success rates and integration bounds based on successes and total N
SucRatesIntBounds <- function(Si, Ni, Sii, Nii, smin, snames=c("sample1","sample2"))
{
# --------------------------------------------------
# Success rates and integration bounds in the case
# of the (conservative) Bayes-Laplace prior
# --------------------------------------------------
# necessary variables:
# N, S(uccesses) for sample 1 (i) and 2 (ii)

################################################################################
# defintion of constants
 Di <- (Si + 1) / (Ni + 2)
 si <- sqrt(Di * (1 - Di) / (Ni + 3))
 Dii <- (Sii + 1) / (Nii + 2)
 sii <- sqrt(Dii * (1 - Dii) / (Nii + 3))

################################################################################
# calculation
 Nmax <- max(Ni, Nii)
 Nmin <- min(Ni, Nii)
 sL <- floor(1000 * sqrt((Nmax+1) / (Nmax+3)) / (Nmax + 2)) / 1000
 sH <- ceiling(1000 / (2 * sqrt(Nmin + 3))) / 1000
 L <- floor(100 * (smin + 1) / (Nmax + 2)) / 100
 H <- 1 - L

 res <- list(Si=Si, Ni=Ni, Sii=Sii, Nii=Nii, smin=smin,
             Di=Di, si=si, Dii=Dii, sii=sii,
             L=L, H=H, sL=sL, sH=sH,
             snames=snames)
 attr(res,"typ") <- c("SRIB")
return(res)
}
#call
#res.SIB <- SucRatesIntBounds(Si=11, Ni=15, Sii=10, Nii=16, smin=0, snames=c("voluntary","non-voluntary"))


################################################################################
# --------------------------------------------------
# ON THE DIFFERENCE IN MEANS
# --------------------------------------------------
DiffinMeans <- function(inval=NULL, out=FALSE, smin=0)
{
################################################################################
# necessary variables in the function
#{ NN, DD, Dsi, Dsii, DsD, ss,
#  dd, lownum, upnum, low, up, psv, PSV,
#  zz, lowinum, upinum, lowiinum, upiinum, psbarv, PSbarV,
#  psvbar, PSVbar, psbarvbar, PSbarVbar, cc,
#  sv, sbarv, svbar, sbarvbar,
#  samemeans, diffmeans, samevars, diffvars, diffsets },

################################################################################
# extract values for calculations from inputv object (class "UMS")

  Di <- inval[["Di"]]
  Dii <- inval[["Dii"]]
  si <- inval[["si"]]
  sii <- inval[["sii"]]
  Ni <- inval[["Ni"]]
  Nii <- inval[["Nii"]]
  L <- inval[["L"]]
  H <- inval[["H"]]
  sL <- inval[["sL"]]
  sH <- inval[["sH"]]
  snames <- inval[["snames"]]
  smin <- inval[["smin"]]
  
 Di
 Dii
 si
 sii
 Ni
 Nii
 L
 H
 sL
 sH
 snames

################################################################################
# defintion of constants
 NN <- Ni+Nii
 DD <- (Ni * Di + Nii * Dii) / NN
 Dsi <- (Ni-1) / Ni * si^2 + Di^2
 Dsii <- (Nii-1) / Nii * sii^2 + Dii^2
 DsD <- (Ni * Dsi + Nii * Dsii) / NN
 ss <- sqrt(NN * (DsD - DD^2) / (NN-1))

 NN
 DD
 Dsi
 Dsii
 DsD
 ss

################################################################################
# p_sv
 dd <- NN * (DsD - DD^2)
 lownum <- NN * (L-DD)^2
 upnum  <- NN * (H-DD)^2

#low = lownum / (2*s^2)
#up  = upnum / (2*s^2)
 integpsv <- function(s)
 { 1 / (s^NN) * exp(-dd / (2 * s^2)) *
   ( pgamma(upnum/(2*s^2),1/2)*gamma(1/2) + pgamma(lownum/(2*s^2),1/2)*gamma(1/2) )
 }
 psv <- integrate(integpsv, lower=sL, upper=sH)$value
 psv
#???	(* + if (L-DD) < 0 *)
 PSV <- psv / sqrt(2*NN)
 PSV

################################################################################
# p_sbarv
 zz <- Ni * (Dsi - Di^2) + Nii * (Dsii - Dii^2)
 lowinum <- Ni * (L - Di)^2
 upinum <- Ni * (H - Di)^2
 lowiinum <- Nii * (L - Dii)^2
 upiinum <- Nii * (H - Dii)^2

 integpsbarv <- function(s)
 {
  1/(s^(NN-1)) * exp(-zz/(2*s^2)) *
  ( pgamma(upinum/(2*s^2),1/2)*gamma(1/2) + pgamma(lowinum/(2*s^2),1/2)*gamma(1/2) ) *
  ( pgamma(upiinum/(2*s^2),1/2)*gamma(1/2) + pgamma(lowiinum/(2*s^2),1/2)* gamma(1/2) )
 }
 psbarv <- integrate(integpsbarv, lower=sL, upper=sH)$value
 psbarv
#??? (* + if (L-DD) < 0 *)
 PSbarV <- psbarv / (2*(H-L) * sqrt(Ni*Nii))
 PSbarV

################################################################################
# p_svbar
# note Mathematica: gamma[a,z0,z1] = gamma[a,z1] - gamma[a,z0]
 UiA <- function(A) Ni*(Dsi-2*Di*A+A^2)/2
 UiiA <- function(A) Nii*(Dsii-2*Dii*A+A^2)/2
 integpsvbar <- function(A)
 {
# UiA <- Ni*(Dsi-2*Di*A+A^2)/2
# UiiA <- Nii*(Dsii-2*Dii*A+A^2)/2
  1/UiA(A)^(Ni/2) *
  1/UiiA(A)^(Nii/2) *
  ( pgamma(UiA(A)/(sH^2),Ni/2)*gamma(Ni/2) - pgamma(UiA(A)/(sL^2),Ni/2)*gamma(Ni/2) ) *
  ( pgamma(UiiA(A)/(sH^2),Nii/2)*gamma(Nii/2) - pgamma(UiiA(A)/(sL^2),Nii/2)*gamma(Nii/2) )
 }
 psvbar <- integrate(integpsvbar, lower=L, upper=H)$value
 psvbar
 PSVbar <- psvbar / (4*log(sH/sL))
 PSVbar

################################################################################
# p_sbarvbar
 UiA <- function(A) Ni*(Dsi-2*Di*A+A^2)/2
 UiiB <- function(B) Nii*(Dsii-2*Dii*B+B^2)/2
 integpsbarvbar1 <- function(A)
 {
  1/(UiA(A)^(Ni/2)) * ( pgamma(UiA(A)/(sH^2),Ni/2)*gamma(Ni/2) - pgamma(UiA(A)/(sL^2),Ni/2)*gamma(Ni/2) )
 }

 integpsbarvbar2 <- function(B)
 {
  1/(UiiB(B)^(Nii/2)) * ( pgamma(UiiB(B)/(sH^2),Nii/2)*gamma(Nii/2) - pgamma(UiiB(B)/(sL^2),Nii/2)*gamma(Nii/2) )
 }
 psbarvbar <- (integrate(integpsbarvbar1, lower=L, upper=H)$value *
               integrate(integpsbarvbar2, lower=L, upper=H)$value)
 psbarvbar
 PSbarVbar <- psbarvbar / (4*(H-L) * log(sH/sL))
 PSbarVbar

################################################################################
# calculate total probability (denominator of Bayes' Theorem)
# sum of all four hypotheses (probabilities)
 cc <- 1 / (PSV + PSbarV + PSVbar + PSbarVbar)
 cc

################################################################################
# compile input values, constants, and results (probabilities)
 sv <- 100 * cc * PSV
 sbarv <- 100 * cc * PSbarV
 svbar <- 100 * cc * PSVbar
 sbarvbar <- 100 * cc * PSbarVbar

 samemeans <- sv + svbar
 diffmeans <- sbarv + sbarvbar
 samevars <- sv + sbarv
 diffvars <- svbar + sbarvbar
 diffsets <- svbar + sbarv + sbarvbar

################################################################################
# calculate various odds ratios
 OR.diffmeans <- diffmeans/samemeans
 OR.samemeans <- samemeans/diffmeans
 OR.diffvars <- diffvars/samevars
 OR.samevars <- samevars/diffvars
 OR.diffsets <- diffsets/sv
 OR.samesets <- sv/diffsets

#compile results
#inputs
 constants <- data.frame(sample1=snames[1],sample2=snames[2],Di,Dii,Ni,si,Nii,sii,NN,DD,ss,L,H,sL,sH,Dsi,Dsii,DsD,dd,lownum,upnum,lowinum,upinum,lowiinum,upiinum,smin)
#unnormalized + total prob
 posteriorprobs.un <- data.frame(cc,PSV,PSbarV,PSVbar,PSbarVbar,psv,psbarv,psvbar,psbarvbar)
#normalized
 posteriorprobs.no <- data.frame(sv,sbarv,svbar,sbarvbar,samemeans,diffmeans,samevars,diffvars,samesets=100-diffsets,diffsets)
#odds ratios
 ORs <- data.frame(OR.diffmeans,OR.samemeans,OR.diffvars,OR.samevars,OR.diffsets,OR.samesets)

 constants
 posteriorprobs.un
 posteriorprobs.no
 ORs

# resulting table/ dataframe with input values
 input.df <- data.frame("No." = c(Ni,Nii,NN),
                        "Standard Deviation" = c(si,sii,ss),
                        "Mean" = c(Di,Dii,DD),
                        "Data set" = c(snames,"combined"),
                        check.names = FALSE)

# table/ dataframe with prior values/ information
 prior.df <- data.frame("Numerical Example" = c("Prior Mean lower bound",
                                                "Prior Mean upper bound",
                                                "Prior Standard Deviation lower bound",
                                                "Prior Standard Deviation upper bound"),
                        "Value" = c(L,H,sL,sH),
                        check.names = FALSE)

# table/ dataframe with resulting probabilities
 prob.res.df <- data.frame("Hypothesis (Probabilities)" = c(
                                            "Same Mean,       same Standard Deviation",
                                            "Different Means, same Standard Deviation",
                                            "Same Mean,       different Standard Deviations",
                                            "Different Means, different Standard Deviations",
                                            "The Means are the same",
                                            "The Means are different",
                                            "The Standard Deviations are the same",
                                            "The Standard Deviations are different",
                                            "The Data Sets are the same",
                                            "The Data Sets are different"),
                           "abbrev" = c("C&S","Cbar,S","C,Sbar","Cbar&Sbar","C","Cbar","S","Sbar","C&S","Cbar|Sbar"),
						   "p(H|D1,D2,I)" = c(sv,sbarv,svbar,sbarvbar,samemeans,diffmeans,samevars,diffvars,100-diffsets,diffsets),
                           check.names = FALSE)

# table/ dataframe with odds ratio results
 OR.res.df <- data.frame("Hypothesis in favor of ..." = c(
                                          "A difference in Means",
                                          "The same Means",
                                          "A difference in Standard Deviations",
                                          "The same Standard Deviations",
                                          "A difference in the Sets (different Means and/ or Standard Deviations)",
                                          "The same Sets (same Means and/ or Standard Deviations)"),
                         "Odds Ratio (OR)" = c(OR.diffmeans,OR.samemeans,
                                          OR.diffvars,OR.samevars,
                                          OR.diffsets,OR.samesets),
                         check.names = FALSE)
# print tables
 if(out)
 {
  cat(paste("\nShort output of 'the difference in means':\n\n"))
  print(input.df, right=FALSE, row.names=FALSE)
  cat(paste("\n"))
  print(prior.df, right=FALSE, row.names=FALSE)
  cat(paste("\n"))
  print(prob.res.df, right=FALSE, row.names=FALSE)
  cat(paste("\n"))
  print(OR.res.df, right=FALSE, row.names=FALSE)
  cat(paste("\n"))
 }

 res <- list(consts = constants,
             pp.un = posteriorprobs.un,
             pp.no = posteriorprobs.no,
             OR = ORs,
             iv.df = input.df,
             prior.df = prior.df,
             prob.df = prob.res.df,
             OR.df = OR.res.df)


return(res)
}
#call
#DiM.results <- DiffinMeans(inval=res.SIB, out=FALSE)


################################################################################
# print nice output
UMSprint <- function(results, dig=4)
{

 cat("\n##############################\n")
 cat("### ON THE DIFFERENCE IN MEANS\n")
 cat("### G.L. Bretthorst (1993)\n")

 consts <- results[["consts"]]
 Ni <- consts$Ni
 Nii <- consts$Nii
 Di <- consts$Di
 Dii <- consts$Dii
 si <- consts$si
 sii <- consts$sii
 NN <- consts$NN
 DD <- consts$DD
 ss <- consts$ss
 smin <- consts$smin
 L <- consts$L
 H <- consts$H
 sL <- consts$sL
 sH <- consts$sH
 
 pp.un <- results[["pp.un"]]
 PSV <-  pp.un$PSV
 PSbarV <-  pp.un$PSbarV
 PSVbar <- pp.un$PSVbar
 PSbVbar <-  pp.un$PSbVbar
 PSbarVbar <-  pp.un$PSbarVbar
 cc <- pp.un$cc
 
 pp.no <- results[["pp.no"]]
 sv <- pp.no$sv
 sbarv <- pp.no$sbarv
 svbar <- pp.no$svbar
 sbarvbar <- pp.no$sbarvbar
 samemeans <- pp.no$samemeans
 diffmeans <- pp.no$diffmeans
 samevars <- pp.no$samevars
 diffvars <- pp.no$diffvars
 diffsets <- pp.no$diffsets
 
 # descriptive statistics
 cat("\n-------------------------------- Descriptive Statistics --------------------\n\n")
 cat(paste("N_1 = ",Ni ," :\t\t\tMean_1 ± SD_1\t\t= ", round(Di,digits=dig)," ± ",round(si,digits=dig), "\n", sep=""))
 cat(paste("N_2 = ",Nii," :\t\t\tMean_2 ± SD_2\t\t= ", round(Dii,digits=dig)," ± ",round(sii,digits=dig), "\n", sep=""))
 cat(paste("N_total = N_1 + N_2  = ", NN ," :\tMean_comb ± SD_comb\t= ", round(DD,digits=dig)," ± ",round(ss,digits=dig), "\n", sep=""))
 cat("\n")
 cat(paste("Bounds on the Mean (smin = ",smin,"):\t\tMean_L = ",L, ",\tMean_H = ",H,"\n",sep=""))
 cat(paste("Bounds on the Standard Deviation:\t  SD_L = ",sL,",\t  SD_H = ",sH,"\n",sep=""))
 if(L < DD)
   {
    cat(paste("\nMean_L - Mean_comb < 0 = ", (L<DD), "\t(-> '+'-sign between Gamma-fcts o.k.)","\n", sep=""))
   } else cat(paste("\nMean_L - Mean_comb < 0 = ", (L<DD), "\t(-> '+'-sign between Gamma-fcts false!)","\n", sep=""))

# results
 cat("\n-------------------------------- Results -----------------------------------\n\n")
 cat(paste("p(sv | D_1, D_2, I)\t\t= const. ", signif(PSV,digits=dig), "\n", sep=""))
 cat(paste("p(sbarv | D_1, D_2, I)\t\t= const. ", signif(PSbarV,digits=dig), "\n", sep=""))
 cat(paste("p(svbar | D_1, D_2, I)\t\t= const. ", signif(PSVbar,digits=dig), "\n", sep=""))
 cat(paste("p(sbarvbar | D_1, D_2, I)\t= const. ", signif(PSbarVbar,digits=dig),"\n\n", sep=""))
 cat(paste("where\t\tconst.\t= ", signif(1/(4*(H-L) * log(sH/sL) * (2*pi)^(NN/2)),digits=dig),
 " / p(D_1,D_2|,I)\n\t\t\t= ",signif(cc,digits=dig), "\n",sep=""))

# model probabilities
 cat(paste("\n--------------- Model --------------------------------- Probability --------\n\n", sep=""))
 cat(paste("sv:\t\tSame Mean,      Same Variance:\t\t", round(sv,digits=dig), "\n", sep=""))
 cat(paste("sbarv:\t\tDifferent Mean, Same Variance:\t\t", round(sbarv,digits=dig), "\n", sep=""))
 cat(paste("svbar:\t\tSame Mean,      Different Variance:\t", round(svbar,digits=dig), "\n", sep=""))
 cat(paste("sbarvbar:\tDifferent Mean, Different Variance:\t", round(sbarvbar,digits=dig), "\n", sep=""))

# odds ratios
 cat("\n------------------------------ Odds Ratios ---------------------------------\n\n")
 cat(paste("The probability the means are the same is:  ", round(samemeans,digits=dig), "\n", sep=""))
 cat(paste("The probability the means are different is: ", round(diffmeans,digits=dig), "\n\n", sep=""))
 if(diffmeans > samemeans)
 {
  cat(paste("The odds ratio is ", round(diffmeans / samemeans,digits=dig),
             " to 1 in favor of different means.\n", sep=""))
 } else
  cat(paste("The odds ratio is ", round(samemeans / diffmeans,digits=dig),
    	       " to 1 in favor of the same means.\n", sep=""))
 cat("\n")

 cat(paste("The probability the variances are the same is:  ", round(samevars,digits=dig),
           "\nThe probability the variances are different is: ", round(diffvars,digits=dig), "\n\n", sep=""))
 if(diffvars > samevars)
 {
  cat(paste("The odds ratio is ", round(diffvars / samevars,digits=dig),
             " to 1 in favor of different variances\n", sep=""))
 } else
  cat(paste("The odds ratio is ", round(samevars / diffvars,digits=dig),
	       	  " to 1 in favor of the same variances\n", sep=""))
 cat("\n")

 cat(paste("The probability the data sets are the same is:  ", round(sv,digits=dig), "\n", sep=""))
 cat(paste("The probability the data sets are different is: ", round(diffsets,digits=dig), "\n\n", sep=""))
 if(diffsets > sv)
 {
  cat(paste("The odds ratio is ", round(diffsets / sv,digits=dig),
            " to 1 in favor of different means and/ or variances.\n",
            sep=""))
 } else
  cat(paste("The odds ratio is ", round(sv / diffsets,digits=dig),
      	 	  " to 1 in favor of the same means and variances.\n",sep=""))
 cat("\n-------------------------------- End ---------------------------------------\n\n")
 cat("\n")

}
################################################################################
#call
#UMSprint(results=DiM.results)
#
#call with writing output to file
#sink("UMSprint.txt")
#UMSprint(results=DiM.results)
#sink()


################################################################################
# plot graphical comparison
# only qualitative Bretthorst
UMSplot <- function(inval, dig=4, pdfout=FALSE, fname="UMSplot.pdf", loga=TRUE, fac=1.1)
{

  Si <- inval[["Si"]] 
  Sii <- inval[["Sii"]]
  Di <- inval[["Di"]]
  Dii <- inval[["Dii"]]
  si <- inval[["si"]]
  sii <- inval[["sii"]]
  Ni <- inval[["Ni"]]
  Nii <- inval[["Nii"]]
  L <- inval[["L"]]
  H <- inval[["H"]]
  sL <- inval[["sL"]]
  sH <- inval[["sH"]]
  snames <- inval[["snames"]]
  smin <- inval[["smin"]]

# determining plot range
 sigma <- min(si, sii)
 plr <- round(10 / (sqrt(2 * pi) * sigma) + 5) / 10
 sek <- seq(from=0, to=1, by=0.01)

# check whether to use the log() to plot (factorials!) = default
 if(loga)
 {
  consti.l <- lfactorial(Ni + 1) - ( lfactorial(Si) + lfactorial(Ni - Si) )
  constii.l <- lfactorial(Nii + 1) - ( lfactorial(Sii) + lfactorial(Nii - Sii) )
  pBLi.l <- function(x) { x^Si * (1 - x)^(Ni - Si) }
  probs.i <- exp(consti.l + log(pBLi.l(sek)))
  pBLii.l <- function(x) { x^Sii * (1 - x)^(Nii - Sii) }
  probs.ii <- exp(constii.l + log(pBLii.l(sek)))
  probs.i
  probs.ii
 } else {
# defining functions for later plotting
  consti <- factorial(Ni + 1) / ( factorial(Si) * factorial(Ni - Si) )
  pBLi <- function(x) { consti * x^Si * (1 - x)^(Ni - Si) }
  constii <- factorial(Nii + 1) / ( factorial(Sii) * factorial(Nii - Sii) )
  pBLii <- function(x) { constii * x^Sii * (1 - x)^(Nii - Sii) }
 # calculating probability densities
  probs.i <- pBLi(sek)
  probs.ii <- pBLii(sek)
 }

# calculate vertical range
 ylim.max <- max(c(probs.i,probs.ii)) * fac

### plot graphical comparison p(H[delta]|S_1/N_1,I) versus p(H[delta]|S_2/N_2,I)
 if(pdfout) pdf(fname,width=9,height=6,paper="A4r")
 par(mar=c(5,6,5,5))
 plot(sek, probs.i,
      xlab="", ylab="probabilities",
      main="",
      type="l", lty=1, lwd=1.75, col="red", bty="l",
      ylim=c(0,ylim.max))
 points(sek, probs.ii, type="l", lty="dashed", lwd=1.75, col="blue")
 mtext(expression(paste(delta," Mean Difference",sep="")), 3, line=2, cex=1.5)
 mtext(eval(substitute(expression(paste("qualitative comparison of ",Si,"/",Ni," (red) vs. ",Sii,"/",Nii," (blue)")),
       list(Si=Si,Ni=Ni,Sii=Sii,Nii=Nii))), 3, line=0.8)
 mtext(expression(paste(delta," = ", mu["1"]," – ",mu["2"])), 1, line=3, cex=1.5)
 mtext(
 eval(substitute(expression(paste(bar(x)["1"] ," = ",Di," | ",
                                  bar(x)["2"] ," = ",Dii," | ",
                                  sigma["1"] ," = ",si, " | ",
                                  sigma["2"] ," = ",sii, "")),list(Dii=round(Dii,dig),Di=round(Di,dig),si=round(si,dig),sii=round(sii,dig)))),
       4,line=1, cex=0.9)
                            
# add a nice legend to explain curves and associated probability densities
 legend("topleft", legend=c(expression(paste("p(H[",x["1"],"] | S" ["1"] ,"/ N" ["1"],", I)")),
                            expression(paste("p(H[",x["2"],"] | S" ["2"] ,"/ N" ["2"],", I)"))),
                            text.col=c("red","blue"), bty="n", cex=1, y.intersp=1.4,
                            title="probability densities", title.col="black", title.adj=1.1)
 if(pdfout) dev.off()
# alternative: add another nice legend with information
#legend("topright", legend=c(eval(substitute(expression(paste(bar(delta)["1"] ," = ",Di," (mean)")),list(Di=round(Di,dig)))),
#                            eval(substitute(expression(paste(bar(delta)["2"] ," = ",Dii," (mean)")),list(Dii=round(Dii,dig)))),
#                            eval(substitute(expression(paste(sigma["1"] ," = ",si, " (SD)")),list(si=round(si,dig)))),
#                            eval(substitute(expression(paste(sigma["2"] ," = ",sii, " (SD)")),list(sii=round(sii,dig))))),
#                            text.col=c("black"),
#                            bty="n", cex=1,
#                            y.intersp=1.4,
#                            title="values", title.col="black")

}
################################################################################
#call
#UMSplot(inputv=res.SIB)
#
#call with pdf output
#UMSplot(inputv=res.SIB,pdfout=TRUE,fname="UMSplot_1998_p48.pdf")

# --------------------------------------------------
# End: ON THE DIFFERENCE IN MEANS
# --------------------------------------------------

