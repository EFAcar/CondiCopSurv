################################################################################
# Analysis of Diabetic Retinopathy Data 
################################################################################
# Geerdens, Acar and Janssen (2017)
# Conditional copula models for right-censored clustered event time data
################################################################################


# Bootstrap Confidence Intervals (run in parallel on cluster)
################################################################

cat("Running the code for B=1. Run it in parallel to obtain bootstrap confidence intervals.\n")

B=1  # change to 1000
conf = 0.90

BootCI.PMar = matrix(NA,nrow=B ,ncol=nx)
BootCI.NPMar = matrix(NA,nrow=B ,ncol=nx)

for(b in 1:B){
  
  time.0 = Sys.time()
  
  ind = sample(1:n,n,replace=TRUE)
  boot.data = work.data[ind,]
  
  par.hat_PMar= CondiCopSurv(Y1=boot.data[,1], Y2=boot.data[,2], status1=boot.data[,3], status2=boot.data[,4], 
                             X=boot.data[,5], x=x, family=3, 
                             mar.method1 = "Weibull", mar.method2 = "Weibull", 
                             cop.band=42, update.eta = F,
                             return.options = c(par=T, test = F, band=F,  pval=F))$par.alter.x 
  
  
  par.hat_NPMar= CondiCopSurv.default(Y1=boot.data[,1], Y2=boot.data[,2], status1=boot.data[,3], status2=boot.data[,4], 
                              X=boot.data[,5], x=x, family=3, 
                              mar.method1 = "Beran", mar.method2 = "Beran", 
                              mar.band1=17, mar.band2=17, 
                              cop.band=42,
                              return.options = c(par=T, test=F,band=F, pval=F))$par.alter.x 
  
  
  # Tau calculation returns error if any parameter estimate >100. Therefore, replace these runs:
  cond = all(par.hat_PMar <100 & par.hat_PMar > 10e-5 ) &  all(par.hat_NPMar <100 & par.hat_NPMar > 10e-5 )
  
  if(cond){
    BootCI.PMar[b,]  = VecBiCopPar2Tau(family=3, par.hat_PMar)
    BootCI.NPMar[b,] = VecBiCopPar2Tau(family=3, par.hat_NPMar)
  }else{b=b-1}

  print(b)
  print(Sys.time()-time.0)
}

# Mean values over B=1000 Bootstrap samples:
mean.tau.hat_PMar = apply(BootCI.PMar, 2, mean, na.rm=T)
mean.tau.hat_NPMar = apply(BootCI.NPMar, 2, mean, na.rm=T)

# Lower bound for CI
low.tau.hat_PMar = apply(BootCI.PMar, 2, quantile, prob=0.5*(1-conf), na.rm=T)
low.tau.hat_NPMar = apply(BootCI.NPMar, 2, quantile, prob=0.5*(1-conf),na.rm=T)

# Upper bound for CI
high.tau.hat_PMar = apply(BootCI.PMar, 2, quantile, prob=1-0.5*(1-conf), na.rm=T)
high.tau.hat_NPMar = apply(BootCI.NPMar, 2, quantile, prob=1-0.5*(1-conf),na.rm=T)
