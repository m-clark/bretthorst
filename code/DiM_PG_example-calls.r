rm(list=ls())
source("DiM_PG.r")

################################################################################
# Phil Gregory
inputvalues <- list(snames = c("riverB.1","riverB.2"),
# sample 1
 d1 = c(13.2,13.8,8.7,9,8.6,9.9,14.2,9.7,10.7,8.3,8.5,9.2),
# sample 2
 d2 = c(8.9,9.1,8.3,6,7.7,9.9,9.9,8.9),

# Input priors and no. of steps in evaluation of p(r|D_1,D_2,I) & p(delta|D_1,D_2,I)
# ndelta = number of steps in delta parameter (mean difference)
 ndelta = 1000, #100
# nr = number of steps in r parameter (ratio of the standard deviations)
 nr = 1000, # 100

# Set prior limits (assumed the same for each data set) on mean (low,high),
# and prior limits (assumed the same for each data set) on the
# standard deviation (sigmalow, sigmahigh).
# upper mean
 high = 12,
# lower mean
 low = 7,
# upper sd
 sigma.high = 4,
# lower sd
 sigma.low = 1)

 # call according to Phil Gregory scheme
inputvalues
dim.res <- DiM.pg(invtyp="pg", inputvalues, print.res=TRUE)
str(dim.res)
plot.DiM(DiM.res=dim.res, filling=TRUE)


################################################################################
# INPUT VALUES
# BASED ON THE METHOD OF U M STUDER
# 
# Phil Gregory
# import original values
#
#> print(input.df, right=FALSE, row.names=FALSE)
# No. Standard Deviation Mean     Data set
# 12  2.177085           10.31667 riverB.1
#  8  1.279997            8.58750 riverB.2
# 20  2.025593            9.62500 combined
#> print(prior.df, right=FALSE, row.names=FALSE)
# Numerical Example                                   Value
# Prior mean lower bound                                7  
# Prior mean upper bound                               12  
# Prior standard deviation lower bound                  1  
# Prior standard deviation upper bound                  4  
# Number of steps for plotting p(delta | D_1, D_2, I) 100  
# Number of steps for plotting p(r | D_1, D_2, I)     100  
#> print(prob.res.df, digits=dig, right=FALSE, row.names=FALSE)
# Hypothesis                                                 Probability
# C,S       = same mean, same standard deviation             0.09758    
# Cbar,S    = different means, same standard deviation       0.28915    
# C,Sbar    = same mean, different standard deviations       0.14429    
# Cbar,Sbar = different means, different standard deviations 0.46898    
# C         = means are the same                             0.24187    
# Cbar      = means are different                            0.75813    
# S         = standard deviations are the same               0.38673    
# Sbar      = standard deviations are different              0.61327    
# C,S       = same means and standard deviations             0.09758    
# Cbar,Sbar = one or both are different                      0.90242    
#> print(oddsratio.res.df, digits=dig+1, right=FALSE, row.names=FALSE)
# Hypothesis                                                               Odds Ratio
# The odds ratio in favour of a difference (means)                         3.1345    
# The odds ratio in favour of a difference (standard deviations)           1.5858    
# The odds ratio in favour of a difference (means and standard deviations) 9.2481    

umsvalues <- list(snames=c("riverB.1","riverB.2"),
                    si=2.1771, Ni=12,
                    sii=1.28, Nii=8,
                    Di=10.3167, Dii=8.5875,
                    L=7, H=12, sL=1, sH=4,
                    ndelta=1000, nr=1000)
umsvalues
ums2pg(umsvalues)
pgbyums.res <- DiM.pg(invtyp="ums", inputvalues=umsvalues, print.res=TRUE)
pgbyums.res
plot.DiM(DiM.res=pgbyums.res, filling=TRUE)


################################################################################
# INPUT VALUES
# BASED ON THE METHOD OF U M STUDER
#
# UM Studer (1998, p.47)
# 
#> print(input.df, right=FALSE, row.names=FALSE)
# No. Standard Deviation Mean      Data set     
# 15  0.1074000          0.7058800 voluntary    
# 16  0.1118400          0.6111100 non-voluntary
# 31  0.1181302          0.6569665 combined     
#> print(prior.df, right=FALSE, row.names=FALSE)
# Numerical Example                                   Value  
# Prior mean lower bound                                0.050
# Prior mean upper bound                                0.950
# Prior standard deviation lower bound                  0.052
# Prior standard deviation upper bound                  0.118
# Number of steps for plotting p(delta | D_1, D_2, I) 100.000
# Number of steps for plotting p(r | D_1, D_2, I)     100.000
#> print(prob.res.df, digits=dig, right=FALSE, row.names=FALSE)
# Hypothesis                                                 Probability
# C,S       = same mean, same standard deviation             0.20275    
# Cbar,S    = different means, same standard deviation       0.49970    
# C,Sbar    = same mean, different standard deviations       0.07608    
# Cbar,Sbar = different means, different standard deviations 0.22147    
# C         = means are the same                             0.27883    
# Cbar      = means are different                            0.72117    
# S         = standard deviations are the same               0.70245    
# Sbar      = standard deviations are different              0.29755    
# C,S       = same means and standard deviations             0.20275    
# Cbar,Sbar = one or both are different                      0.79725    
#> print(oddsratio.res.df, digits=dig+1, right=FALSE, row.names=FALSE)
# Hypothesis                                                               Odds Ratio
# The odds ratio in favour of a difference (means)                         2.58639   
# The odds ratio in favour of a difference (standard deviations)           0.42360   
# The odds ratio in favour of a difference (means and standard deviations) 3.93224   
# The odds ratio in favour of the same (means)                             0.38664   
# The odds ratio in favour of the same (standard deviations)               2.36074   
# The odds ratio in favour of the same (means and standard deviations)     0.25431

# UM Studer 1998, p.47
# call according to UMS scheme
umsvalues <- list(snames=c("voluntary","non-voluntary"),
                    si=0.1074, Ni=15,
                    sii=0.11184, Nii=16,
                    Di=0.70588, Dii=0.61111,
                    L=0.05, H=0.95, sL=0.052, sH=0.118,
                    ndelta=1000, nr=1000)
umsvalues
ums2pg(umsvalues)
ums1.res <- DiM.pg(invtyp="ums", inputvalues=umsvalues, print.res=TRUE)
ums1.res
plot.DiM(DiM.res=ums1.res, filling=TRUE)


################################################################################
# INPUT VALUES
# BASED ON THE METHOD OF U M STUDER
# 
# GL Bretthorst 1993, p.189) after Jaynes 1976 + 1983
#
# No. Standard Deviation Mean     Data set
# 4  6.480000           50.00000 Jaynes.1
#  9  7.480000           42.00000 Jaynes.2
# 13  7.909937           44.46154 combined
#> print(prior.df, right=FALSE, row.names=FALSE)
# Numerical Example                                   Value
# Prior mean lower bound                               34  
# Prior mean upper bound                               58  
# Prior standard deviation lower bound                  3  
# Prior standard deviation upper bound                 10  
# Number of steps for plotting p(delta | D_1, D_2, I) 100  
# Number of steps for plotting p(r | D_1, D_2, I)     100  
#> print(prob.res.df, digits=dig, right=FALSE, row.names=FALSE)
# Hypothesis                                                 Probability
# C,S       = same mean, same standard deviation             0.1677     
# Cbar,S    = different means, same standard deviation       0.4157     
# C,Sbar    = same mean, different standard deviations       0.1069     
# Cbar,Sbar = different means, different standard deviations 0.3097     
# C         = means are the same                             0.2746     
# Cbar      = means are different                            0.7254     
# S         = standard deviations are the same               0.5834     
# Sbar      = standard deviations are different              0.4166     
# C,S       = same means and standard deviations             0.1677     
# Cbar,Sbar = one or both are different                      0.8323     
#> print(oddsratio.res.df, digits=dig+1, right=FALSE, row.names=FALSE)
# Hypothesis                                                               Odds Ratio
# The odds ratio in favour of a difference (means)                         2.64118   
# The odds ratio in favour of a difference (standard deviations)           0.71414   
# The odds ratio in favour of a difference (means and standard deviations) 4.96299   
# The odds ratio in favour of the same (means)                             0.37862   
# The odds ratio in favour of the same (standard deviations)               1.40028   
# The odds ratio in favour of the same (means and standard deviations)     0.20149   

# call according to UMS scheme
inputvalues <- list(snames=c("Jaynes.1","Jaynes.2"),
                    si=6.48, Ni=4, sii=7.48, Nii=9,
                    Di=50, Dii=42,
                    L=34, H=58, sL=3, sH=10,
                    ndelta=1000, nr=1000)
inputvalues
ums2pg(inputvalues)
jaynes.res <- DiM.pg(invtyp="ums", inputvalues, print.res=TRUE)
jaynes.res
plot.DiM(DiM.res=jaynes.res, filling=TRUE)


