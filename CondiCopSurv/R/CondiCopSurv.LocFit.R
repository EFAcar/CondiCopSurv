#' @title Local Likelihood Estimation of the Dependence Parameter for Right Censored Survival Data
#' @description This function performs the local likelihood estimation of the calibration function in conditional copula models for right censored survival data
#'
#' @param u1,u2 numeric vectors of equal length taking values in [0,1]
#' @param status1,status2 vectors of censoring indicators taking values 0 (censored) or 1 (uncensored)
#' @param X vector of covariate values
#' @param x fixed numeric value in the covariate range
#' @param eta.init initial value of calibration parameters for optimization
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
#' @param optim.method method used in optimization; options are "optim" and "nlminb"
#' @param control.finite logical indicator to exclude nonfinite loglikelihood contributions (default value = FALSE)
#' @keywords internal
CondiCopSurv.LocFit = function(u1, u2, status1, status2, X, x, eta.init, family, p=1, 
                               band, Kern=KernEpa, optim.method = "nlminb", control.finite=TRUE){
  
  if(any(u1 < 0) || any(u1 >1) || any(u2 < 0) || any(u2 >1))
    stop("Provide uniform data for the conditional marginal survival functions.")
  
  if(!(family %in% c(3, 4, 5)))
    stop("Copula family not yet implemented.")
  
  if(!(p %in% c(0, 1)))
    stop(paste("Local polynomial order not yet implemented."))
  
  if(missing(eta.init)){
    null.fit = optimize(CondiCopSurv.LLik ,interval=Null.init(family), family=family, p=0, u1=u1, u2=u2, status1=status1,status2=status2, X=X, control.finite=control.finite, negative=TRUE)
    if(p==0){eta.init = null.fit$minimum }
    if(p==1){eta.init = c(null.fit$minimum, 0)
    eta.init = optim(par=eta.init, fn=CondiCopSurv.LLik, family=family, p=1, u1=u1, u2=u2, status1=status1,status2=status2, X=X, control.finite=control.finite, negative=TRUE)$par}
  }
  
  if(length(eta.init)!=p+1)
    stop(paste("Dimension of the calibration parameter should be", p+1))
  
  
  obj.fnc = function(eta){
    CondiCopSurv.LocLLik(eta, family=family, p=p, u1=u1, u2=u2, status1=status1 ,status2=status2, X=X, x=x, band=band, Kern=Kern, control.finite=control.finite)
  }
  
  if(optim.method=="optim"){
    opt.res = optim(eta.init,obj.fnc, method="L-BFGS-B")
    return(list(eta=opt.res$par, LocLLik =-opt.res$value, convergence = opt.res$convergence))
  }
  
  if(optim.method=="nlminb"){
    opt.res = nlminb(eta.init, obj.fnc)
    return(list(eta=opt.res$par, LocLLik =-opt.res$objective, convergence = opt.res$convergence))
  }
}