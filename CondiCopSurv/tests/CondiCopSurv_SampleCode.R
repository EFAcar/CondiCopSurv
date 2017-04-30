###############################################
# Sample Code for the CondiCopSurv Package
###############################################

# An Illustration of the Analysis with Simulated Data

library(CondiCopSurv)
## library(survival)
## library(VineCopula)

VecBiCopPar2Tau = Vectorize(BiCopPar2Tau)
VecBiCopTau2Par = Vectorize(BiCopTau2Par)

# Clayton Copula - Concave Model - Sample Size n=250 - Moderate Censoring

family = 3
n=250
censoring = "moderate"

# Specification of the Kendall's Tau function:

tau.fnc = function(t){ -0.1*(t-3)^2 +0.7}

# Specification of the marginal parameters under the Weibull model

lambda1 = 0.5 ; rho1 = 1.5 ; b1 = 0.8
lambda2 = 0.5;  rho2 = 1.5 ; b2 = 0.8

true.mar.par = c(rho1, rho2, lambda1, lambda2, b1, b2)

# Specification of the parameters for the censoring time distribution

if(censoring == "low"){ lambdac = 1.5;  rhoc = 1.5}  # Low Censoring: about 20 % 
if(censoring == "moderate"){ lambdac = 1.5;  rhoc = 0.5}  # Moderate Censoring: about 50 %

# Specification of the pilot bandwidth values and fixed points

pilot.val = c(0.3, 0.5, 0.8, 1.2, 2, 3)
x = seq(2, 5, by=0.1)



############################################################################
#  1. Data Generation
############################################################################

# Part I: Generate data from Conditional Copula

X = runif(n,2,5)
true.tau = tau.fnc(X)
true.theta = VecBiCopTau2Par(family,true.tau)

if(family==3){true.eta = log(true.theta)}
if(family==4){true.eta = log(true.theta-1)}
if(family==5){true.eta = true.theta}

Udata = t(mapply(BiCopSim, 1, family=family, par=true.theta))


# Part II:  Obtain event times under the Weibull model (reparametrization is required)

U.temp1 = Udata[,1]
U.temp2 = Udata[,2]

shape1 = rho1 
shape2 = rho2
scale1 = (lambda1 * exp(b1*X))^(-1/rho1)
scale2 = (lambda2 * exp(b2*X))^(-1/rho2)

  shapec = rhoc
  scalec = (lambdac)^(-1/rhoc)
  
  ET1 = qweibull(U.temp1, shape = shape1, scale = scale1, lower.tail= F)
  ET2 = qweibull(U.temp2, shape = shape2, scale = scale2, lower.tail= F)
  CenT = rweibull(n, shape = shapec, scale = scalec)
  
  Y1 = pmin(ET1, CenT)
  Y2 = pmin(ET2, CenT)
  status1 = 1*(ET1 <= CenT)
  status2 = 1*(ET2 <= CenT)
  

############################################################################
#  2. Obtain Estimation and Test Results - Parametric Margins
############################################################################

# Fit Marginal Models
  
system.time({
    fit.PMar1 = CondiMarFit(Y1, status1, X=X, mar.method ="Weibull")
    fit.PMar2 = CondiMarFit(Y2, status2, X=X, mar.method ="Weibull")
})


# Estimation at a fixed point for a given bandwidth value
  
system.time({
    fit.PMar.Est.x = CondiCopSurv(Y1,  Y2,  status1, status2,  X=X, x=x[15], family=family, 
                                  mar.method1 ="Weibull", mar.method2="Weibull", cop.band = 2, 
                                  return.options = c(par=TRUE, test = FALSE, band=FALSE, pval=FALSE))
})
  
  
# Estimation at grid points for a given bandwidth value

system.time({
  fit.PMar.Est.grid = CondiCopSurv(Y1,  Y2,  status1, status2,  X=X, x=x, family=family, 
                               mar.method1 ="Weibull", mar.method2="Weibull", cop.band = 2, 
                               return.options = c(par=TRUE, test = FALSE, band=FALSE, pval=FALSE))
})



# Bandwidth Selection and Estimation 

system.time({
  fit.PMar.Est = CondiCopSurv(Y1,  Y2,  status1, status2,  X=X, family=family, 
                              mar.method1 ="Weibull", mar.method2="Weibull", cop.pilot = pilot.val, 
                              return.options = c(par=TRUE, test = FALSE, band=TRUE, pval=FALSE))
})



# Bandwidth Selection, Estimation and Test Statistic 

system.time({
  fit.PMar.TEst = CondiCopSurv(Y1,  Y2,  status1, status2,  X=X, family=family, 
                              mar.method1 ="Weibull", mar.method2="Weibull", cop.pilot = pilot.val, 
                              return.options = c(par=TRUE, test = TRUE, band=TRUE, pval=FALSE))
})



# Bandwidth Selection, Estimation, Testing and bootstrap-based p-value calculation (for B=1)

system.time({
  fit.PMar.TEst = CondiCopSurv(Y1,  Y2,  status1, status2,  X=X, family=family, 
                               mar.method1 ="Weibull", mar.method2="Weibull", cop.pilot = pilot.val,
                               return.options = c(par=TRUE, test = TRUE, band=FALSE, pval=TRUE), 
                               nBoot=1, # change for p-value calculation
                               boot.options = c(update.marpar = TRUE, update.copband = FALSE, update.marband=FALSE))
})






############################################################################
#  3. Obtain Estimation and Test Results - NonParametric Margins
############################################################################


# Fit Marginal Models

system.time({
  fit.NPMar1 = CondiMarFit(Y1, status1, X=X, mar.method ="Beran", pilot=pilot.val)
  fit.NPMar2 = CondiMarFit(Y2, status2, X=X, mar.method ="Beran", pilot=pilot.val)
})



# Estimation at a fixed point for a given bandwidth value

system.time({
  fit.NPMar.Est.x = CondiCopSurv(Y1,  Y2,  status1, status2,  X=X, x=x[15], family=family, 
                                mar.method1 = "Beran", mar.method2 = "Beran", 
                                mar.band1 = 2, mar.band2 = 2, cop.band = 2, 
                                return.options = c(par=TRUE, test = FALSE, band=FALSE, pval=FALSE))
})


# Estimation at grid points for a given bandwidth value

system.time({
  fit.NPMar.Est.grid = CondiCopSurv(Y1,  Y2,  status1, status2,  X=X, x=x, family=family, 
                                 mar.method1 = "Beran", mar.method2 = "Beran", 
                                 mar.band1 = 2, mar.band2 = 2, cop.band = 2, 
                                 return.options = c(par=TRUE, test = FALSE, band=FALSE, pval=FALSE))
})

# Bandwidth Selection and Estimation 

system.time({
  fit.NPMar.Est = CondiCopSurv(Y1,  Y2,  status1, status2,  X=X, x=x, family=family, 
                              mar.method1 = "Beran", mar.method2 = "Beran", 
                              mar.band1 = 2, mar.band2 = 2, cop.pilot = pilot.val, 
                              return.options = c(par=TRUE, test = FALSE, band=TRUE, pval=FALSE))
})


# Bandwidth Selection, Estimation, Testing

system.time({
  fit.NPMar.TEst = CondiCopSurv(Y1,  Y2,  status1, status2,  X=X, x=x, family=family, 
                                mar.method1 ="Beran", mar.method2="Beran", 
                                mar.pilot = pilot.val, cop.pilot = pilot.val, 
                                return.options = c(par=TRUE, test = TRUE, band=TRUE, pval=FALSE))
})


# Bandwidth Selection, Estimation, Testing and bootstrap-based p-value calculation (B=1)

system.time({
  fit.NPMar.TEst = CondiCopSurv(Y1,  Y2,  status1, status2,  X=X, x=x, family=family, 
                                mar.method1 ="Beran", mar.method2="Beran", 
                                mar.pilot = pilot.val, cop.pilot = pilot.val, 
                                return.options = c(par=TRUE, test = TRUE, band=TRUE, pval=TRUE), 
                                nBoot=1, # change for p-value calculation
                                boot.options = c(update.marpar = FALSE, update.copband = FALSE, update.marband=FALSE))
  
})




############################################################################
#  4. A Graphical Display of the Estimation Results in Kendall's tau scale
############################################################################

tau.PMar.Est = VecBiCopPar2Tau(fit.PMar.Est.grid$par.alter.x, family=family)
tau.NPMar.Est = VecBiCopPar2Tau(fit.NPMar.Est.grid$par.alter.x, family=family)

plot(x, tau.fnc(x), xlim=c(2,5), ylim = c(0,1),  xlab = expression(x), ylab = expression(tau(x)), type="l", lwd=2)
lines(x, tau.PMar.Est, col=2, lwd=2)
lines(x, tau.NPMar.Est, col=3, lwd=2)
legend("topright", paste(c("True Model", "Weibull", "Beran")),  cex=0.6, fill=c(1,2,3))




