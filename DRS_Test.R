################################################################################
# Analysis of Diabetic Retinopathy Data 
################################################################################
# Geerdens, Acar and Janssen (2017)
# Conditional copula models for right-censored clustered event time data
################################################################################


# LRT: Constant versus Linear Calibration Model 
#######################################################

# Conditional Copula Model with Weibull Margins

null_PMar.Clayton = optimize(CondiCopSurv.LLik ,interval=c(-10,5), family=3, p=0, 
                             u1=U.data1[,1], u2=U.data1[,2], status1=U.data1[,3], status2=U.data1[,4], 
                             X=X, control.finite=TRUE, negative=TRUE)

linear_PMar.Clayton = optim(par=c(null_PMar.Clayton$minimum, 0), fn=CondiCopSurv.LLik, family=3, p=1, 
                            u1=U.data1[,1], u2=U.data1[,2], status1=U.data1[,3], status2=U.data1[,4], 
                            X=X, control.finite=TRUE, negative=TRUE)

1-pchisq(2*(-linear_PMar.Clayton$value + null_PMar.Clayton$objective), df=1)
# 0.111

# Conditional Copula Model with Beran Margins

null_NPMar.Clayton   = optimize(CondiCopSurv.LLik ,interval=c(-10,5), family=3, p=0, 
                                u1=U.data2[,1], u2=U.data2[,2], status1=U.data2[,3], status2=U.data2[,4], 
                                X=X, control.finite=TRUE, negative=TRUE)

linear_NPMar.Clayton = optim(par= c(null_NPMar.Clayton$minimum, 0), fn=CondiCopSurv.LLik, family=3, p=1, 
                             u1=U.data2[,1], u2=U.data2[,2], status1=U.data2[,3], status2=U.data2[,4], 
                             X=X, control.finite=TRUE, negative=TRUE)

1-pchisq(2*(-linear_NPMar.Clayton$value + null_NPMar.Clayton$objective), df=1)
# 0.102


Tau.constant_PMar.Clayton = BiCopPar2Tau(family=3, exp(null_PMar.Clayton$minimum) )
Tau.constant_NPMar.Clayton = BiCopPar2Tau(family=3, exp(null_NPMar.Clayton$minimum) )

Tau.linear_PMar.Clayton = VecBiCopPar2Tau(family=3, exp(linear_PMar.Clayton$par[1]+linear_PMar.Clayton$par[2]*x))
Tau.linear_NPMar.Clayton = VecBiCopPar2Tau(family=3, exp(linear_NPMar.Clayton$par[1]+linear_NPMar.Clayton$par[2]*x))


# GLRT and Bootstrap P-Value Calculation (run in parallel on cluster)
######################################################################

# Use the CondiCopSurv to perform estimation and testing along with p-value calculation.  
# For the latter specify "pval=T" in return.options.


# Conditional Copula Model with Weibull margins 

cat("Running the code for B=1. Run it in parallel to obtain bootstrap p-values.\n")

fit.PMar.Clayton = CondiCopSurv( Y1=Y1, Y2=Y2, status1=status1, status2=status2, 
                                 X=X, x=x, family=3,
                                 mar.method1="Weibull", mar.method2="Weibull", 
                                 cop.band=42, eta.init=eta.init.C1, update.eta = F,
                                 nBoot=1, # change to 1000
                                 return.options = c(par=T, test=T,band=T, pval=T),
                                 boot.options = c(update.marpar = T, update.copband = F, update.marband=F))


# Conditional Copula Model with Beran margins 


fit.NPMar.Clayton = CondiCopSurv( Y1=Y1, Y2=Y2, status1=status1, status2=status2, 
                                  X=X, x=x, family=3,
                                  mar.method1="Beran", mar.method2="Beran", 
                                  mar.band1=17, mar.band2=17, 
                                  cop.band=42, eta.init= eta.init.C2, update.eta = F,
                                  nBoot=1,  # change to 1000
                                  return.options = c(par=T, test=T,band=T, pval=T), 
                                  boot.options = c(update.marpar = T, update.copband = F, update.marband=F))


pval.Par.Clayton  =  fit.PMar.Clayton$pval
pval.NPar.Clayton =  fit.NPMar.Clayton$pval

# p-values: 0.164 and 0.168

