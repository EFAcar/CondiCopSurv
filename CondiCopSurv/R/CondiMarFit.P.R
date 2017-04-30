#' @title Fitting the Conditional Survival Functions
#'
#' @description This function fits conditional marginal survival functions for right censored survival data.
#'
#' @param Y vector of event times
#' @param status vector of censoring indicators taking values 0 (censored) or 1 (uncensored)
#' @param X vector of covariate values
#' @param mar.par marginal parameters of the Weibull model: (optional) estimated if missing. 
#' @keywords internal
CondiMarFit.P = function(Y, status, X, par){
  if(missing(par)){
    fit <- survival::survreg(Surv(Y, status)~X, dist="weibull")
    unname(coefficients(fit))
    hat.rho <- 1/(fit$scale)
    hat.lambda <- exp(-fit$coefficients[[1]]/fit$scale)
    hat.b <-  -fit$coefficients[[2]]/fit$scale
    par <- c(hat.rho, hat.lambda, hat.b)
  }
  u <- exp(-par[2]* Y^par[1] * exp(par[3]*X))
  return(list(u=u, par = par))
}