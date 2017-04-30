#' @title Fitting a Conditional Copula to Right Censored Survival Data 
#' @description This function performs the local likelihood estimation of the calibration function in conditional copula models for right censored survival data when the conditional margins are estimated using the Weibull model.
#'
#' @param u1,u2 numeric vectors of equal length taking values in [0,1]
#' @param status1,status2 vectors of censoring indicators taking values 0 (censored) or 1 (uncensored)
#' @param X vector of covariate values
#' @param x vector of covariate values at which the estimation is performed
#' @param eta.init initial value(s) of calibration parameter for optimization (optional)
#' @param family an integer defining the bivariate copula family:
#'
#' 3 = Clayton copula
#'
#' 4 = Gumbel copula
#'
#' 5 = Frank copula
#' @param p degree of local polynomial in the calibration model: p=0 for local constant, p=1 for local linear
#' @param band bandwidth value used in the conditional copula estimation
#' @param Kern Kernel function
#' @param optim.method method used in optimization, current options are "optim" and "nlminb"
#' @param control.finite logical indicator to exclude nonfinite loglikelihood contributions (default value = FALSE)
#' @param update.eta logical indicator to update the initial value of the calibration parameter for optimization
#' @param position integer value specifying the starting location for the function estimation 
#' @export
CondiCopSurv.Fit = function(u1, u2, status1, status2, X, x, eta.init=NULL, 
                            family, p=1, band, Kern=KernEpa, optim.method = "nlminb", 
                            control.finite=TRUE, update.eta=FALSE, position){
  
  if(!(family %in% c(3, 4, 5)))
    stop("Copula family not yet implemented.")
  
  if(!(p %in% c(0, 1)))
    stop(paste("Local polynomial order not yet implemented."))
  
  data = data.frame(u1, u2, status1, status2, X)
  n = nrow(data)
  
  if(missing(x)){x=X}
  
  if(!missing(eta.init) & nrow(matrix(eta.init, ncol=(p+1)))>1){
    if(nrow(eta.init) != length(x)) 
      stop("Dimension of initial values does not match to the length of x.")
  }
  
  
  if(missing(eta.init)){
    null.fit = optimize(CondiCopSurv.LLik ,interval=Null.init(family), family=family, p=0, 
                        u1=u1, u2=u2, status1=status1,status2=status2, X=X, 
                        control.finite=control.finite, negative=TRUE)
    if(p==0){eta.init = null.fit$minimum }
    if(p==1){eta.init = c(null.fit$minimum, 0)
    eta.init = optim(par=eta.init, fn=CondiCopSurv.LLik, family=family, p=1, 
                     u1=u1, u2=u2, status1=status1,status2=status2, X=X, 
                     control.finite=control.finite, negative=TRUE)$par
    }
  }
  
  

  if(length(x)==1){
    eta.estim = CondiCopSurv.LocFit(u1=u1, u2=u2, status1=status1, status2=status2, X = X, x=x, 
                                    eta.init=eta.init, family=family, p=p, band=band, Kern=Kern, 
                                    optim.method = optim.method, control.finite=control.finite)$eta[1]
  }
  if(length(x)>1){
    if(update.eta==FALSE){
      if(nrow(matrix(eta.init, ncol=2))==1){
        easy.fnc = function(i){CondiCopSurv.LocFit(u1=u1, u2=u2, status1=status1, status2=status2, X = X, x=x[i], 
                                                   eta.init=eta.init, family=family, p=p, band=band, Kern=Kern, 
                                                   optim.method = optim.method, control.finite=control.finite)$eta[1]
          }
      } 
      if(nrow(matrix(eta.init, ncol=2)) > 1){
        easy.fnc = function(i){CondiCopSurv.LocFit(u1=u1, u2=u2, status1=status1, status2=status2, X = X, x=x[i], 
                                                   eta.init=eta.init[i,], family=family,p=p, band=band, Kern=Kern, 
                                                   optim.method = optim.method, control.finite=control.finite)$eta[1]
          }
      } 
      eta.estim = try(sapply(1:length(x), easy.fnc), silent=TRUE)
   
       }else{
      
      eta.estim = rep(NA, length(x))
      if(missing(position)) {start.val = ceiling(median(1:length(x)))}
      if(!missing(position)){start.val = position}
      
      x.new <- sort(x, index.return=T)$x
      x.ind <- sort(x, index.return=T)$ix
      if(nrow(matrix(eta.init, ncol=2)) > 1){eta.init = eta.init[start.val,]}
      
      # sequential fitting with updated parameters
      
      fit0 = CondiCopSurv.LocFit(u1, u2, status1, status2, X, x=x.new[start.val], eta.init=eta.init, 
                                 family = family, p=p, band=band, Kern=Kern, optim.method = optim.method, 
                                 control.finite=control.finite)$eta
      eta.estim[start.val] = fit0[1]
      fit=fit0
      
      for(i in (start.val-1):1){
        try_fit = try(CondiCopSurv.LocFit(u1, u2, status1, status2, X, x=x.new[i], eta.init=fit, family = family, p=p, band=band, Kern=Kern, optim.method = optim.method, control.finite=control.finite),TRUE)
        fail_ind=is(try_fit,"try-error")
        diff.fit = ifelse(sum(!is.na(eta.estim))==1, fit0[2],mean(diff(eta.estim), na.rm=T))
        stp=1
        
        
        while((fail_ind | abs(try_fit$eta[1]-fit[1]) > 1.5*diff.fit) & stp < 10){
          try_init = c(fit[1], diff.fit + runif(1, -0.1,0.1))
          try_fit = try(CondiCopSurv.LocFit(u1, u2, status1, status2, X, x=x.new[i], eta.init=try_init, family = family, p=p, band=band, Kern=Kern, optim.method = optim.method, control.finite=control.finite),TRUE)
          fail_ind=is(try_fit,"try-error")
          stp=stp+1
        }
        if(!fail_ind){
          fit = try_fit$eta
          eta.estim[i] = fit[1]
        }
      }
      
      fit = fit0
      for(i in (start.val+1):length(x)){
        try_fit = try(CondiCopSurv.LocFit(u1, u2, status1, status2, X, x=x.new[i], eta.init=fit, family = family, p=p, band=band, Kern=Kern, optim.method = optim.method, control.finite=control.finite),TRUE)
        fail_ind = is(try_fit,"try-error")
        diff.fit = ifelse(sum(!is.na(eta.estim[-i]))==1, fit0[2],mean(diff(eta.estim[-i]), na.rm=T))
        stp=1
        
        
        while((fail_ind | abs(try_fit$eta[1]-fit[1]) > 1.5*diff.fit) & stp < 10){
          try_init = c(fit[1], diff.fit + runif(1, -0.1,0.1))
          try_fit = try(CondiCopSurv.LocFit(u1, u2, status1, status2, X, x=x.new[i], eta.init=try_init, family = family, p=p, band=band, Kern=Kern, optim.method = optim.method, control.finite=control.finite),TRUE)
          fail_ind=is(try_fit,"try-error")
          stp=stp+1
        }
        if(!fail_ind){
          fit = try_fit$eta
          eta.estim[i] = fit[1]
        }
      }
      reverse.index = match(x,x.new)
      eta.estim = eta.estim[reverse.index]
    }
  }
  
  if(family==3){ cop.par = exp(eta.estim) }
  if(family==4){ cop.par = exp(eta.estim)+1 }
  if(family==5){ cop.par = eta.estim }
  
  if(length(x)==length(X) && x==X){
    alter_ll  = rep(NA,length(x))
    for(j in 1:length(X)){
      alter_ll[j]= CondiCopSurv.LLik(eta=eta.estim[j],family=family, p=0, u1=u1[j], u2=u2[j], status1=status1[j],status2=status2[j], X=X[j], control.finite=control.finite, negative=FALSE)
    }
    alter_LL = sum(alter_ll)
    
    return(list(eta = eta.estim, cop.par = cop.par, LLik = alter_LL))
    
  }else{
    return(list(eta = eta.estim, cop.par = cop.par))
  }
}