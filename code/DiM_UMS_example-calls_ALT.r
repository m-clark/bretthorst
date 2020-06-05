rm(list=ls())
source("DiM_UMS.r")

################################################################################
# On the difference in means
# paper: G.L. Bretthorst "On the difference of means" (1993)
# http://bayes.wustl.edu/glb/diff.pdf
#
# R code based on Mathematica code by Urban Studer (90's, Zürich/ CH)


################################################################################
#works
res.SIB.NRFmf <- SucRatesIntBounds(Si=20, Ni=47, Sii=13, Nii=28, smin=0, snames=c("male","female"))
res.SIB.NRFmf
DiM.res.NRFmf <- DiffinMeans(inval=res.SIB.NRFmf, out=FALSE)
str(DiM.res.NRFmf)
UMSprint(results=DiM.res.NRFmf)
DiM.res.NRFmf$OR.df
DiM.res.NRFmf$prob.df
UMSplot(inval=res.SIB.NRFmf,pdfout=FALSE)


################################################################################
#EXAMPLE FOR FAILURE
# does not work with these values
res.SIB.NRFtotal <- SucRatesIntBounds(Si=(20+13), Ni=(47+28), Sii=(338 %/% 4), Nii=338, smin=0, snames=c("male","female"))
res.SIB.NRFtotal
DiM.res.NRF <- DiffinMeans(inval=res.SIB.NRFtotal, out=FALSE)
str(DiM.res.NRF)
UMSprint(results=DiM.res.NRFtotal)
#WATCH THE PLOT - BECOMES OVIOUS (ie. the result), but should work nevertheless
UMSplot(inval=res.SIB.NRFtotal,pdfout=FALSE)


################################################################################
# example values for success rates and integration bounds
# Studer 1998, p.47
res.SIB <- SucRatesIntBounds(Si=11, Ni=15, Sii=10, Nii=16, smin=0, snames=c("voluntary","non-voluntary"))
res.SIB
DiM.res.SIB <- DiffinMeans(inval=res.SIB, out=FALSE)
str(DiM.res.SIB)
UMSprint(results=DiM.res.SIB)
UMSplot(inval=res.SIB,pdfout=FALSE)


################################################################################
# example values for success rates and integration bounds
# Studer 1998, p.48
res.SIB <- SucRatesIntBounds(Si=20, Ni=31, Sii=17, Nii=27, smin=0, snames=c("voluntary","non-voluntary"))
res.SIB
DiM.res.SIB <- DiffinMeans(inval=res.SIB, out=FALSE)
str(DiM.res.SIB)
UMSprint(results=DiM.res.SIB)
UMSplot(inval=res.SIB,pdfout=FALSE)


################################################################################
# Gregory, 200x, p.xxx
inval <- list(Si=NULL,
     Ni = 12,
     Sii=NULL,
     Nii = 8,
     smin = 0,
     Di = 10.3167,
     si = 2.1771,
     Dii = 8.5875,
     sii = 1.28,
     L = 7,
     H = 12,
     sL = 1,
     sH = 4,
     snames = c("river.1","river.2")
    )
DiM.PG <- DiffinMeans(inval=inval, out=FALSE)
DiM.PG$prob.df
DiM.PG$OR.df

# Studer 1998, p.47
inval <- list(Si=NULL,
     Ni = 15,
     Sii=NULL,
     Nii = 16,
     smin = 0,
     Di = 0.70588,
     si = 0.1074,
     Dii = 0.61111,
     sii = 0.11184,
     L = 0.05,
     H = 0.95,
     sL = 0.052,
     sH = 0.118,
     snames = c("voluntary","non-voluntary")
    )
DiM.UMS <- DiffinMeans(inval=inval, out=FALSE)
DiM.UMS$prob.df
DiM.UMS$OR.df

# Bretthorst, 1993, p.189 (from Jaynes, 1976 + 1983)
inval <- list(Si=NULL,
     Ni = 4,
     Sii=NULL,
     Nii = 9,
     smin = 0,
     Di = 50,
     si = 6.48,
     Dii = 42,
     sii = 7.48,
     L = 34,
     H = 58,
     sL = 3,
     sH = 10,
     snames = c("Jaynes.1","Jaynes.2")
    )
DiM.ETJ <- DiffinMeans(inval=inval, out=FALSE)
DiM.ETJ$prob.df
DiM.ETJ$OR.df
################################################################################

######################################################
#GENERAL call and preparation
# Studer 1998, p.47
# Bretthorst, 1993, p.189 (from Jaynes, 1976 + 1983)	
inval <- list(
     Si=NULL,		#UMS specific -> successes group1
     Ni = 4,		#N group1
     Sii=NULL,		#UMS specific -> successes group2
     Nii = 9,		#N group2
     smin = 0,		#UMS specific -> bounds on the mean -> only for SucRatesIntBounds()
     Di = 50,		#mean group1
     si = 6.48,		#sd group1
     Dii = 42,		#mean group2
     sii = 7.48,	#sd group2
     L = 34,		#mean lower bound
     H = 58,		#mean upper bound
     sL = 3,		#variance lower bound
     sH = 10,		#variance upper bound
     snames = c("Jaynes.1","Jaynes.2")
    )
######################################################
