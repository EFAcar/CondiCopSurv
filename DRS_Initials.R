################################################################################
# Analysis of Diabetic Retinopathy Data 
################################################################################
# Geerdens, Acar and Janssen (2017)
# Conditional copula models for right-censored clustered event time data
################################################################################


# Initial Values for Local Likelihood Estimation
#####################################################

eta.Clayton.Par <- eta.Clayton.NPar <- rep(NA,nx)
h.PMar <- h.NPMar <-42

# Grid Search to Find Initial Values:

grid.C = expand.grid(seq(-0.5,0.5,by=0.1),seq(-0.5,0.5,by=0.1))
init.1 = Find.init(grid.C, u1=U.data1[,1], u2=U.data1[,2], status1=U.data1[,3], status2=U.data1[,4], 
                   X=X, x=x[1], family = 3, p=1, band=h.PMar)
init.2 = Find.init(grid.C, u1=U.data2[,1], u2=U.data2[,2], status1=U.data2[,3], status2=U.data2[,4], 
                   X=X, x=x[1], family = 3, p=1, band=h.NPMar)

eta.init.C1 <- eta.init.C2 <- matrix(NA, ncol=2, nrow=nx)
  
for(i in 1:nx){
  
  fit1.C = CondiCopSurv.LocFit( u1=U.data1[,1], u2=U.data1[,2], status1=U.data1[,3], status2=U.data1[,4], 
                                X=X, x=x[i], family = 3, band=h.PMar, eta.init=init.1)
  eta.Clayton.Par[i]  = fit1.C$eta[1]
  init.1 = fit1.C$eta

  
  fit2.C = CondiCopSurv.LocFit( u1=U.data2[,1], u2=U.data2[,2], status1=U.data2[,3], status2=U.data2[,4], 
                                X=X, x=x[i], family = 3, band=h.NPMar, eta.init=init.2)
  eta.Clayton.NPar[i]  = fit2.C$eta[1]
  init.2 = fit2.C$eta

  
  # Very few subjects with diabetes onset age >50
  # Hence carefully choose the initial values
  
  if(x[i] >50){
    grid.C = expand.grid(seq(1.58,1.78,by=0.01),seq(1,3,by=0.1))
    init.1 = Find.init( grid.C, u1=U.data1[,1], u2=U.data1[,2], status1=U.data1[,3], status2=U.data1[,4], 
                        X=X, x=x[i], family = 3, p=1, band=h.PMar)
    init.2 = Find.init(grid.C, u1=U.data2[,1], u2=U.data2[,2], status1=U.data2[,3], status2=U.data2[,4], 
                       X=X, x=x[i], family = 3, p=1,  band=h.NPMar)
  }
  eta.init.C1[i,] = init.1
  eta.init.C2[i,] = init.2
}




