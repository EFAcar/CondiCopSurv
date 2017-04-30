#' @title Fitting the Conditional Survival Functions
#'
#' @description This function fits conditional marginal survival functions for right censored survival data.
#'
#' @param Y vector of event times
#' @param status vector of censoring indicators taking values 0 (censored) or 1 (uncensored)
#' @param X vector of covariate values
#' @param mar.method method to obtain conditional marginal survival fuctions: "Weibull" or "Beran" 
#' @export
CondiMarFit = function(Y, status, X, mar.method=c("Weibull", "Beran"),...){
  mar.method <- match.arg(mar.method)
  if(mar.method == "Weibull"){ return(CondiMarFit.P(Y=Y,  status=status, X=X, ...))}
  if(mar.method == "Beran"){ return(CondiMarFit.NP(Y=Y, status=status, X=X, ...))}
} 