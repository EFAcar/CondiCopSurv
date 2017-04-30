#' @title Conditional Copula Model Fitting and Testing for Right Censored Survival Data (Default)
#' @description Default function that performs estimation and testing in conditional copula models for right censored survival data.
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
#' @param cop.band bandwidth value used in the conditional copula estimation
#' @param mar.pilot,cop.pilot user-specified pilot bandwidth values (optional)
#' @param nband.mar,nband.cop desired length of bandwidth values if pilot values are not specified
#' @param band.method method used in bandwidth selection: three cross-validation based selectors 
#' @param eta.init initial value(s) of calibration parameter for optimization (optional)
#' @param p degree of local polynomial in the calibration model: p=0 for local constant, p=1 for local linear
#' @param Kern Kernel function
#' @param optim.method method used in optimization, current options are "optim" and "nlminb"
#' @param control.finite logical indicator to exclude nonfinite loglikelihood contributions (default value = TRUE)
#' @param update.eta logical indicator to update the initial value of the calibration parameter for optimization
#' @param position integer value specifying the starting location for the function estimation 
#' @param return.options logical vector specifying which results to return 
CondiCopSurv.default = function(Y1, Y2, status1, status2, X, x, family, 
                                mar.method1 = c("Weibull","Beran"), mar.method2 = c("Weibull","Beran"), 
                                mar.par1, mar.par2, 
                                mar.band1, mar.band2, cop.band,
                                mar.pilot = NULL,cop.pilot = NULL, 
                                nband.mar = 10, nband.cop = 10, 
                                band.method = c("CV_Y","CV_T","CV_use"), 
                                eta.init, p = 1,  Kern = KernEpa, 
                                optim.method = "nlminb", control.finite = TRUE, 
                                update.eta = FALSE, position,
                                return.options = c(par=TRUE, test=TRUE, band=TRUE),...){

  if(missing(mar.par1)) mar.par1 <- NULL
  if(missing(mar.par2)) mar.par2 <- NULL
  if(missing(mar.band1)) mar.band1 <- NULL
  if(missing(mar.band2)) mar.band2 <- NULL
  
  band.method = match.arg(band.method)
  
  if(!(family %in% c(3, 4, 5)))
    stop("Copula family not yet implemented.")
  
  if(!(p %in% c(0, 1)))
    stop("Local polynomial order not yet implemented.")
  
  if(sum(return.options)==0)
    stop("There is nothing to return.")
  
  # Marginal Models 
  
  if(!is.null(mar.par1)){fit1 <- CondiMarFit(Y1, status1, X, mar.method=mar.method1, par=mar.par1, ... )}        #known parametric margins
  if(!is.null(mar.band1)){fit1 <- CondiMarFit(Y1, status1, X, mar.method=mar.method1, band=mar.band1, ... )}	   #nonparametric margins with known bandwidth
  
  if(!is.null(mar.par2)){fit2 <- CondiMarFit(Y2, status2, X, mar.method=mar.method2, par=mar.par2, ... )}        #known parametric margins
  if(!is.null(mar.band2)){fit2 <- CondiMarFit(Y2, status2, X, mar.method=mar.method2, band=mar.band2, ... )}	   #nonparametric margins with known bandwidth
  
  if(is.null(mar.par1) & mar.method1=="Weibull"){fit1 <- CondiMarFit(Y1, status1, X, mar.method=mar.method1, ... )}	 #parametric margins (estimated)
  if(is.null(mar.par2) & mar.method2=="Weibull"){fit2 <- CondiMarFit(Y2, status2, X, mar.method=mar.method2, ... )}	 #parametric margins (estimated)
  
  if(is.null(mar.band1) & mar.method1=="Beran"){
    fit1 <- CondiMarFit(Y1, status1, X, mar.method=mar.method1, pilot=mar.pilot, band.method=band.method, ... )}  #nonparametric margins (bandwidth selected)
  if(is.null(mar.band2) & mar.method2=="Beran"){
    fit2 <- CondiMarFit(Y2, status2, X, mar.method=mar.method2, pilot=mar.pilot, band.method=band.method, ... )}  #nonparametric margins (bandwidth selected)
  
  u1 = fit1$u
  u2 = fit2$u
  
  mar.par1 = fit1$par
  mar.par2 = fit2$par
  
  mar.band1 = fit1$band
  mar.band2 = fit2$band
  
  if(any(u1 < 0) || any(u1 >1) || any(u2 < 0) || any(u2 >1))
    stop("Provide uniform data for the conditional marginal survival functions.")
  
  data = data.frame(u1, u2, status1, status2,X)
  n = nrow(data)
  
  # Global Estimation
  
  fit.null <- optimize(CondiCopSurv.LLik,interval=Null.init(family), family=family, p=0, 
                       u1=u1, u2=u2, status1=status1,status2=status2, X=X, 
                       control.finite=control.finite, negative=TRUE)
  eta.null = fit.null$minimum
  null_LL = -fit.null$objective
  
  if(family==3){  par.null = exp(eta.null) }
  if(family==4){  par.null = exp(eta.null )+1 }
  if(family==5){  par.null = eta.null }
  
  
  # Local Estimation
  
  
    if(p==0){eta.init.0 = eta.null}
    if(p==1){eta.init.0 = c(eta.null, 0)
    eta.init.0 <- optim(par=eta.init.0, fn=CondiCopSurv.LLik, family=family, p=1, u1=u1, u2=u2, 
                      status1=status1,status2=status2, X=X,  
                      control.finite=control.finite, negative=TRUE)$par
    }

  
  if(!missing(cop.band)){choose.band=NULL}
  
  if(missing(cop.band) || is.null(cop.band)){
    choose.band = BandSelect(pilot = cop.pilot, nband = nband.cop, u1=u1, u2=u2, 
                             status1=status1, status2=status2, X=X, 
                             eta.init = eta.init.0, family=family, p=p, Kern=Kern,  
                             optim.method = optim.method, control.finite= control.finite)
    cop.band = choose.band$h
  }
  
  if(!missing(x)){
    
    if(missing(eta.init)){eta.init = eta.init.0}
    
    nx = length(x)
    fit.x <-  CondiCopSurv.Fit( u1=u1, u2=u2, status1=status1, status2=status2,
                                X=X, x=x, eta.init = eta.init, 
                                family=family, p=p, band= cop.band, Kern=Kern, 
                                optim.method = optim.method, control.finite=control.finite, 
                                update.eta = update.eta, position=position)
    eta.x = fit.x$eta
    par.x = fit.x$cop.par
  }
  
  # Test Statistic
  
  if(return.options["test"]){
    
    if(!missing(x) & nrow(matrix(eta.init, ncol=(p+1))) == length(x)){
      eta.init.X = matrix(NA,nrow=length(X), ncol=2 )
      eta.init.X[,1] = approx(x,eta.init[,1], xout=X)$y
      eta.init.X[,2] = approx(x,eta.init[,2], xout=X)$y
    }else{eta.init.X=NULL}
    
    
    fit.alter =  CondiCopSurv.Fit(u1=u1, u2=u2, status1=status1, status2=status2, X=X, 
                                  eta.init = eta.init.X, family=family, p=p, 
                                  band= cop.band, Kern=Kern, optim.method = optim.method, 
                                  control.finite=control.finite, 
                                  update.eta = update.eta, position=position)
    eta.alter = fit.alter$eta
    par.alter = fit.alter$cop.par
    alter_LL  = fit.alter$LLik
    
    lambda_test = alter_LL - null_LL
  }
  
  
  # Return results:
  
  if(return.options["par"]){
    par.list = list( null = par.null,mar1 = mar.par1, mar2= mar.par2)
    if(!missing(x)){ 
      par.list = c(par.list, list(alter.x = par.x))}
    if(return.options["test"]){ 
      par.list = c(par.list, list(alter = par.alter) )}
  }
  
  test.list = NULL
  band.list = NULL
  
  if(return.options["test"]){test.list = list(stat =lambda_test , null_LL = null_LL,  alter_LL = alter_LL )}
  if(return.options["band"]){
    band.list = c(cop=cop.band , mar1 = mar.band1, mar2=mar.band2)}
  
  return(c(par = par.list[return.options["par"]], test = test.list[return.options["test"]] ,  
           band = band.list[return.options["band"]] ) )
  
}