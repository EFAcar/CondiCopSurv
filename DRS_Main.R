################################################################################
# Analysis of Diabetic Retinopathy Data 
################################################################################
# Geerdens, Acar and Janssen (2017)
# Conditional copula models for right-censored clustered event time data
################################################################################

# Data source:
######################################################
# http://www.mayo.edu/research/documents/diabeteshtml/doc-10027460


# Load the R library: CondiCopSurv
######################################################

# To _install_ package, set working directory to where file below is located and run the following command:
# install.packages(pkgs = "CondiCopSurv_1.0.tar.gz",
#                  repos = NULL, type = "source")


library(CondiCopSurv)
## library(survival)
## library(VineCopula)

VecBiCopPar2Tau = Vectorize(BiCopPar2Tau)
VecBiCopTau2Par = Vectorize(BiCopTau2Par)


# Load the data 
######################################################

diaret = read.table(file="DiaRet.txt", header = F)
names(diaret) = c("subject_id","laser_type", "treated_eye", 
                  "age_at_onset", "diabetes_type", 
                  "risk_treated", "status_treated", "time_treated", 
                  "risk_untreated", "status_untreated", "time_untreated")


X  <- diaret$age_at_onset
Y1 <- diaret$time_treated
Y2 <- diaret$time_untreated
status1 <- diaret$status_treated
status2 <- diaret$status_untreated
  
work.data <- data.frame(Y1, Y2, status1, status2, X)
n <- nrow(work.data)
n   # 197


# Set fixed points and pilot h values
######################################################

x <- seq(1, 58, by=1)
nx <- length(x)  # 58
pilot.val <- c(3, 5, 6, 8, 12, 17, 23, 31, 42, 57)


##################################################################
# CLAYTON COPULA
##################################################################

# See "DRS_Weibull.R" for details of parametric conditional marginal models
## source("DRS_Weibull.R")
#
# See "DRS_Beran.R" for details of nonparametric conditional marginal models
## source("DRS_Beran.R")
#
# See "DRS_Bandwidth.R" for bandwidth selection under different families
## source("DRS_Bandwidth.R")
#
# See "DRS_Initials.R" for determinations of initial values used in estimation 
## source("DRS_Initials.R")


# Conditional Copula Model with Weibull margins 
##################################################################


fit.PMar.Clayton = CondiCopSurv( Y1=Y1, Y2=Y2, status1=status1, status2=status2, 
                                 X=X, x=x, family=3,
                                 mar.method1="Weibull", mar.method2="Weibull", 
                                 cop.band=42, eta.init=eta.init.C1, update.eta = F,
                                 return.options = c(par=T, test=T,band=T, pval=F))


# Conditional Copula Model with Beran margins 
##################################################################

fit.NPMar.Clayton = CondiCopSurv( Y1=Y1, Y2=Y2, status1=status1, status2=status2, 
                                  X=X, x=x, family=3,
                                  mar.method1="Beran", mar.method2="Beran", 
                                  mar.band1=17, mar.band2=17, 
                                  cop.band=42, eta.init= eta.init.C2, update.eta = F,
                                  return.options = c(par=T, test=T,band=T, pval=F))

Tau.hat_PMar.Clayton =  VecBiCopPar2Tau(family=3, fit.PMar.Clayton$par.alter.x)
Tau.hat_NPMar.Clayton = VecBiCopPar2Tau(family=3, fit.NPMar.Clayton$par.alter.x)


# See "DRS_BootCI.R" for bootstrap confidence intervals (run in parallel on cluster!)
## source("DRS_BootCI.R")
#
# See "DRS_Test.R" for p-value calculation based on bootstrap procedure (run in parallel on cluster!)
## source("DRS_Test.R")
#



# Estimation Results under Clayton Copula:

postscript("DRS_Clayton.eps", height=5.5,width=12.5,
           horizontal = FALSE, onefile = FALSE, paper = "special", pointsize=12)

ind = !is.na(match(x,seq(1,57,2)))
par( mfrow=c(1,2), pin=c(1,1), mar = c(4, 4, 1.5, 0.5) + 0.1)  

pcol= "black" ; npcol="black"

plot(x,  Tau.hat_PMar.Clayton, type="l", col=pcol, ylim=c(0,1), lty=2, lwd=3.5, xlab= expression(x), ylab=expression(widehat(tau)(x)))
#lines(fixed.pts, mean.tau.hat_PMar, col=pcol, lty=2, lwd=2)
lines(x[ind], low.tau.hat_PMar[ind], col=pcol, lty=3, lwd=2.5)
lines(x[ind], high.tau.hat_PMar[ind], col=pcol, lty=3, lwd=2.5)
segments(x0=1, x1=58, y0=Tau.constant_PMar.Clayton, y1=Tau.constant_PMar.Clayton, col="gray65", lwd=2.5, lty=1)
lines(x[ind], Tau.linear_PMar.Clayton[ind], col="gray65", lty=1, lwd=3)

plot(x,  Tau.hat_NPMar.Clayton.x, type="l", col=npcol, ylim=c(0,1), lty=2, lwd=3.5, xlab= expression(x), ylab=expression(widetilde(tau)(x)))
#lines(fixed.pts, mean.tau.hat_NPMar, col=npcol, lty=2, lwd=2)
lines(x[ind], low.tau.hat_NPMar[ind], col=npcol, lty=3, lwd=2.5)
lines(x[ind], high.tau.hat_NPMar[ind], col=npcol, lty=3, lwd=2.5)
segments(x0=1, x1=58, y0=Tau.constant_NPMar.Clayton, y1=Tau.constant_NPMar.Clayton, col="gray65", lwd=2.5, lty=1)
lines(x[ind], Tau.linear_NPMar.Clayton[ind], col="gray65", lty=1, lwd=3)

dev.off()
quartz.options(reset = TRUE)

