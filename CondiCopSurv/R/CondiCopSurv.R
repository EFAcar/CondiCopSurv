#' @title Conditional Copula Model Fitting and Testing for Right Censored Survival Data
#' @description This function performs the local likelihood estimation of the calibration function in conditional copula models for right censored survival data given a covariate. The function also returns an approximate p-value for the Generalized Likelihood Ratio (GLR) test to assess the significance of the covariate effect.
#'
#' @param Y1,Y2 vectors of event times
#' @param status1,status2 vectors of censoring indicators taking values 0 (censored) or 1 (uncensored)
#' @param X vector of covariate values
#' @param x vector of covariate values at which the estimation is performed
#' @param family an integer defining the bivariate copula family:
#'
#' 3 = Clayton copula
#'
#' 4 = Gumbel copula
#'
#' 5 = Frank copula
#' @param mar.method1,mar.method2 method to obtain conditional marginal survival fuctions for each margin: "Weibull" or "Beran"
#' @param mar.par1,mar.par2 marginal parameters of the Weibull model (optional) 
#' @param mar.band1,mar.band2 bandwidth value used in the Beran estimator (optional)
#' @param cop.band bandwidth value used in the conditional copula estimation (optional)
#' @param mar.pilot,cop.pilot user-specified pilot bandwidth values (optional)
#' @param nband.mar,nband.cop desired length of bandwidth values if pilot values are not specified
#' @param band.method method used in bandwidth selection: three cross-validation based selectors 
#' @param eta.init initial value of calibration parameters for optimization
#' @param p degree of local polynomial in the calibration model: p=0 for local constant, p=1 for local linear
#' @param Kern Kernel function
#' @param optim.method method used in optimization, current options are "optim" and "nlminb"
#' @param control.finite logical indicator to exclude nonfinite loglikelihood contributions (default value = TRUE)
#' @param update.eta logical indicator to update the initial value of the calibration parameter for optimization
#' @param return.options logical vector specifying which results to return 
#' @param nBoot number of bootstrap samples in the test procedure 
#' @param boot.options logical vector specifying bootstrap control options
#' @export
CondiCopSurv = function(Y1, Y2, status1, status2, X, x, family, 
                        mar.method1 = c("Weibull", "Beran"), 
                        mar.method2 = c("Weibull", "Beran"), 
                        mar.par1, mar.par2, 
                        mar.band1, mar.band2, cop.band,
                        mar.pilot = NULL, cop.pilot = NULL,
                        nband.mar = 10, nband.cop = 10,
                        band.method = c("CV_Y","CV_T","CV_use"), 
                        eta.init, p = 1, Kern = KernEpa, 
                        optim.method = "nlminb", control.finite = TRUE, 
                        update.eta = FALSE, position,
                        return.options = c(par=TRUE, test = TRUE, band=TRUE, pval=TRUE), 
                        nBoot = 100, 
                        boot.options = c(update.marpar = TRUE, update.copband = FALSE, update.marband=FALSE), ...){


  mar.method1 = match.arg(mar.method1)
  mar.method2 = match.arg(mar.method2)
  band.method = match.arg(band.method)
  
  # Results for original data:
  
  data.res  =  CondiCopSurv.default (Y1=Y1, Y2=Y2, status1=status1, status2=status2, X=X, x=x, 
                                     family=family, mar.method1 = mar.method1, mar.method2 = mar.method2, 
                                     mar.par1 = mar.par1, mar.par2 = mar.par2, 
                                     mar.band1 = mar.band1, mar.band2 = mar.band2, cop.band = cop.band,
                                     mar.pilot = mar.pilot, cop.pilot = cop.pilot, 
                                     nband.mar = nband.mar, nband.cop = nband.cop, 
                                     band.method = band.method,
                                     eta.init = eta.init, p = p, Kern = Kern, optim.method = optim.method, 
                                     control.finite = control.finite, update.eta = update.eta, 
                                     position=position,
                                     return.options =  return.options, ...)
  
  if(!return.options["pval"]){
    return(data.res)
  }
  
  if(return.options["pval"]){
    
    cop.par.sim <-  data.res$par.null
    mar.par.sim <-  c(data.res$par.mar1, data.res$par.mar2)[c(1,4,2,5,3,6)]
    mar.band.sim <- c(data.res$band.mar1, data.res$band.mar2)
    
    if(!boot.options["update.marpar"]){ 
      boot.mar.par1 <- mar.par.sim[c(1,3,5)]
      boot.mar.par2 <- mar.par.sim[c(2,4,6)]
    }else{
      boot.mar.par1 <- NULL 
      boot.mar.par2 <- NULL
    }
    
    if(!boot.options["update.copband"]){ 
      boot.cop.band <- data.res$band.cop
    }else{ boot.cop.band <- NULL}
    
    
    if(!boot.options["update.marband"]){ 
      boot.mar.band1 <- mar.band.sim[1]
      boot.mar.band2 <- mar.band.sim[2]
    }else{
      boot.mar.band1 <- NULL
      boot.mar.band2 <- NULL
    }
    
    
    boot.test  <-  rep(NA, nBoot)
    boot.band.cop <-  rep(NA, nBoot)
    boot.band.mar1 <-  rep(NA, nBoot)
    boot.band.mar2 <-  rep(NA, nBoot)
    
    if(mar.method1!=mar.method2){
      stop("Hybrid-methods for margins are not yet implemented for the bootstrap samples.")}

    b=1
    while(b <= nBoot){
      
      if(mar.method1=="Weibull" && mar.method2=="Weibull"){ 
        
        boot.data  <-  BootSim.P(Y1=Y1, Y2=Y2, status1=status1, status2=status2, X=X, family=family, 
                                 cop.par= cop.par.sim, mar.par = mar.par.sim)
        
        boot.res <- CondiCopSurv.default(Y1=boot.data[,1], Y2=boot.data[,2], 
                                         status1=boot.data[,3], status2=boot.data[,4], 
                                         X=boot.data[,5], x=x, family=family, 
                                         mar.method1 = mar.method1, mar.method2 = mar.method2, 
                                         mar.par1=boot.mar.par1 , mar.par2= boot.mar.par2 , 
                                         cop.band=boot.cop.band, cop.pilot = cop.pilot, 
                                         nband.cop = nband.cop, eta.init = eta.init, p = p, 
                                         Kern = Kern, optim.method = optim.method, 
                                         control.finite = control.finite, update.eta = update.eta,
                                         position=position,...)
        boot.test[b]=boot.res$test.stat   
        boot.band.cop[b]=boot.res$band.cop
      }
      
      if(mar.method1=="Beran" && mar.method2=="Beran"){ 
        
        boot.data = BootSim.NP(Y1=Y1, Y2=Y2, status1=status1, status2=status2, X=X, family=family, 
                               cop.par = cop.par.sim, mar.band = c(mar.band.sim[1], mar.band.sim[2]))
        
        boot.res = CondiCopSurv.default (Y1=boot.data[,1], Y2=boot.data[,2], 
                                         status1=boot.data[,3], status2=boot.data[,4], 
                                         X=boot.data[,5], x=x, family=family, 
                                         mar.method1 = mar.method1, mar.method2 = mar.method2, 
                                         mar.band1= boot.mar.band1, mar.band2=  boot.mar.band2,
                                         cop.band=boot.cop.band,
                                         mar.pilot = mar.pilot, cop.pilot = cop.pilot, 
                                         nband.mar = nband.mar, nband.cop = nband.cop, 
                                         band.method = band.method, 
                                         eta.init = eta.init, p = p, Kern = Kern, 
                                         optim.method = optim.method, 
                                         control.finite = control.finite, update.eta = update.eta,
                                         position=position,...)
        boot.test[b]=boot.res$test.stat   
        boot.band.cop[b]=boot.res$band.cop
        boot.band.mar1[b]=boot.res$band.mar1
        boot.band.mar2[b]=boot.res$band.mar2
      }
      
      if(is.finite(boot.test[b]) & boot.test[b] >0 ){
        if(b==1 | b-trunc(b/10)*10==0) {cat("Completed bootstrap samples =", b, "\n")}
        b = b+1
      }else{b = b}
    }
    
    pval =  sum(boot.test >= data.res$test.stat )/nBoot
    
    
    # Return results:
    
    return(c(data.res, list(pval=pval, boot.test = boot.test, boot.band.cop=boot.band.cop,boot.band.mar1=boot.band.mar1,boot.band.mar2=boot.band.mar2)))
  }}

