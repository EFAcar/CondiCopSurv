#' @title Bandwidth Selection for the Local Likelihood Estimation of the Dependence Parameter for Right Censored Survival Data
#' @description This function performs bandwidth selection for the local likelihood estimation of the calibration function in conditional copula models for right censored survival data given the conditional margins.
#'
#' @param pilot user-specified pilot bandwidth values (optional)
#' @param nband desired length of bandwidth values if pilot values are not specified
#' @param u1,u2 numeric vectors of equal length taking values in [0,1]
#' @param status1,status2 vectors of censoring indicators taking values 0 (censored) or 1 (uncensored)
#' @param X vector of covariate values
#' @param eta.init initial value of calibration parameters for optimization
#' @param family an integer defining the bivariate copula family:
#'
#' 3 = Clayton copula
#'
#' 4 = Gumbel copula
#'
#' 5 = Frank copula
#' @param p degree of local polynomial in the calibration model: p=0 for local constant, p=1 for local linear
#' @param Kern Kernel function
#' @param optim.method method used in optimization, current options are "optim" and "nlminb"
#' @param control.finite logical indicator to exclude nonfinite loglikelihood contributions (default value = FALSE)
#' @param subset.index a vector of indices of data points used in bandwidth selection (optional)
#' @export
BandSelect = function(pilot = NULL, nband = 10, u1, u2, status1, status2, X, 
                      eta.init, family, p=1, Kern=KernEpa, 
                      optim.method = "nlminb", control.finite=TRUE, subset.index = 1:n){

  if(missing(pilot)){
    nband <- nband+2 ; diffs <- diff(sort(X),1)
    h.min <- max(diffs); h.max <- max(X)-min(X)
    log.seq <- seq(from=log(h.min),to=log(h.max), length.out =nband)
    pilot.band <- round(exp(log.seq),2)
    pilot.band <- pilot.band[-(1:2)]
    nband <- length(pilot.band)
  }else{
    pilot.band <- pilot
    nband <- length(pilot.band)
  }

  if(any(u1 < 0) || any(u1 >1) || any(u2 < 0) || any(u2 >1))
    stop("Provide uniform data for the conditional marginal survival functions.")

  if(!(family %in% c(3, 4, 5)))
    stop("Copula family not yet implemented.")

  if(!(p %in% c(0, 1)))
    stop(paste("Local polynomial order not yet implemented."))

  if(missing(eta.init)){
    null.fit <- optimize(CondiCopSurv.LLik ,interval= Null.init(family), family=family, p=0, u1=u1, u2=u2, status1=status1,status2=status2, X=X, control.finite=control.finite, negative=TRUE)
    if(p==0){eta.init <- null.fit$minimum }
    if(p==1){eta.init <- c(null.fit$minimum, 0)
    eta.init <- optim(par=eta.init, fn=CondiCopSurv.LLik, family=family, p=1, u1=u1, u2=u2, status1=status1,status2=status2, X=X, control.finite=control.finite, negative=TRUE)$par}
  }

  if(length(eta.init)!=p+1)
    stop(paste("Dimension of the calibration parameter should be", p+1))

  n <- length(X)
  cvlogdens <- matrix(NA, nrow=n, ncol=nband)

  if(missing(subset.index)){
    index  = 1:n}else{index = subset.index}
  for(i in subset.index){
    easy.fnc <- function(h){CondiCopSurv.LocFit(u1=u1[-i], u2=u2[-i], status1=status1[-i], status2=status2[-i], X = X[-i], x=X[i], eta.init=eta.init, family=family,
                                             p=p, band=h, Kern=Kern, optim.method = optim.method, control.finite=control.finite)$eta[1]}
    eta.estim <- try(sapply(pilot.band, easy.fnc), silent=TRUE)

    if(mode(eta.estim)=="numeric" && all(is.na(eta.estim)==FALSE)){
      eta.estim_v <- matrix(eta.estim, ncol=nband)
      for(j in 1:nband){
      cvlogdens[i,j] <- CondiCopSurv.LLik(  eta= eta.estim_v[1,j],family=family, p=0, u1=u1[i], u2=u2[i], status1=status1[i], status2=status2[i], control.finite=control.finite, negative=FALSE)
      }
    }else{cvlogdens[i,] <- rep(NA,nband) }
  }

  use.ind <- which(apply((apply(cvlogdens, 2, is.finite)), 1, prod)==1)
  cvLIK <- colSums(as.matrix(cvlogdens[use.ind,]))

  if(!all(is.finite(cvLIK)) & !all(is.na(cvLIK))){
    pilot.band <- pilot.band[is.finite(cvLIK) & !is.na(cvLIK)]
    cvLIK <- cvLIK[is.finite(cvLIK) & !is.na(cvLIK)]
  }

  opt.band <-  pilot.band[which(cvLIK==max(cvLIK))]
  return(list(h=opt.band, h.pilot = pilot.band, cvLIK = cvLIK, length = length(use.ind)))

}


