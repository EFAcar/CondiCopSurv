#' @title Conditional Copula Local Loglikelihood Function
#'
#' @description This function gives the objective function for the local likelihood estimation in conditional copula models for right censored survival data
#'
#' @param family an integer defining the bivariate copula family:
#'
#' 3 = Clayton copula
#'
#' 4 = Gumbel copula
#'
#' 5 = Frank copula
#' @param eta calibration parameter
#' @param p degree of local polynomial in the calibration model: p=0 for local constant, p=1 for local linear
#' @param u1,u2 numeric vectors of equal length taking values in [0,1]
#' @param status1,status2 vectors of censoring indicators taking values 0 (censored) or 1 (uncensored)
#' @param X vector of covariate values
#' @param x fixed numeric value in the covariate range
#' @param band bandwidth value used in the conditional copula estimation
#' @param Kern Kernel function
#' @param control.finite logical indicator to exclude nonfinite loglikelihood contributions (default value = FALSE)
#' @param negative logical indicator to return the loglikelihood or negative loglikelihood function (default value = TRUE)
#' @keywords internal
CondiCopSurv.LocLLik <- function( eta, family, p=1, u1, u2, status1, status2, X, x, 
                                  band, Kern=KernEpa, control.finite=FALSE, negative=TRUE){

  if(!(family %in% c(3, 4, 5)))
    stop("Copula family not yet implemented.")

  if(!(p %in% c(0, 1)))
    stop(paste("Local polynomial order not yet implemented."))

  if(length(eta)!=p+1)
    stop(paste("Dimension of the calibration parameter should be", p+1))

  if(family==3){
    copula_model <- function(eta, p, u1, u2, X, x){
      if(p==0){theta <- exp(eta[1])}
      if(p==1){theta <- exp(eta[1]+eta[2]*(X-x))}

      cdf <- ((u1^(-theta)) + (u2^(-theta))-1)^(-1/theta)
      partial1 <- u1^(-(theta+1)) * (u1^(-theta)+u2^(-theta)-1)^(-(1+1/theta))
      partial2 <- u2^(-(theta+1)) * (u1^(-theta)+u2^(-theta)-1)^(-(1+1/theta))
      pdf <- (theta+1)* u1^(-(theta+1)) * u2^(-(theta+1))* (u1^(-theta)+u2^(-theta)-1)^(-(2+1/theta))

      model <- list(cdf=cdf, partial1 = partial1, partial2=partial2,  pdf=pdf)
      return(model)
    }}


  if(family==4){
    copula_model <- function(eta, p, u1, u2, X, x){
      if(p==0){theta <- exp(eta[1])+1}
      if(p==1){theta <- exp(eta[1]+eta[2]*(X-x))+1}

      cdf <- exp(-((-log(u1))^(theta)+(-log(u2))^(theta))^(1/(theta)))
      partial1 <- (1/u1) * (-log(u1))^((theta)-1) * ((-log(u1))^(theta)+(-log(u2))^(theta))^(1/(theta)-1) * exp(-((-log(u1))^(theta)+(-log(u2))^(theta))^(1/(theta)))
      partial2 <- (1/u2) * (-log(u2))^((theta)-1) * ((-log(u1))^(theta)+(-log(u2))^(theta))^(1/(theta)-1) * exp(-((-log(u1))^(theta)+(-log(u2))^(theta))^(1/(theta)))
      pdf <- (1/u1) *(1/u2) * exp(-((-log(u1))^(theta)+(-log(u2))^(theta))^(1/(theta))) * (log(u1)*log(u2))^((theta)-1) * ((-log(u1))^(theta)+(-log(u2))^(theta))^(2*(1/(theta)-1)) * (1+ ((theta)-1)* ((-log(u1))^(theta)+(-log(u2))^(theta))^(-1/(theta)))

      model <- list(cdf=cdf, partial1 = partial1, partial2=partial2,  pdf=pdf)
      return(model)
    }}


  if(family==5){
    copula_model <- function(eta, p, u1, u2, X, x){
      if(p==0){theta <- eta[1]}
      if(p==1){theta <- eta[1]+eta[2]*(X-x)}

      cdf <- (-1/theta)*log(1+(((exp(-theta*u1)-1)*(exp(-theta*u2)-1))/(exp(-theta)-1)))
      partial1 <- exp(-(theta)*u1) * (exp(-(theta)*u2)-1) / ( (exp(-(theta))-1) + (exp(-(theta)*u1)-1) * (exp(-(theta)*u2)-1))
      partial2 <- exp(-(theta)*u2) * (exp(-(theta)*u1)-1) / ( (exp(-(theta))-1) + (exp(-(theta)*u1)-1) * (exp(-(theta)*u2)-1))
      pdf <- ((theta) * (1-exp(-(theta)))* exp(-(theta)*(u1+u2)))/ ((1-exp(-(theta))) - (1-exp(-(theta)*u1)) * (1-exp(-(theta)*u2)))^2

      model <- list(cdf=cdf, partial1 = partial1, partial2=partial2,  pdf=pdf)
      return(model)
    }}

  W <- (1/band)*Kern((X-x)/band)
  ind <- which(W>0)

  W <- W[ind] ; X <- X[ind]
  u1 <- u1[ind]; u2 <- u2[ind]
  status1 <- status1[ind]; status2 <- status2[ind]

  model <- copula_model(eta, p, u1, u2, X, x)

  term1 <- model$cdf
  term2 <- model$partial1
  term3 <- model$partial2
  term4 <- model$pdf

  delta1 <- (1-status1)*(1-status2)
  delta2 <- status1*(1-status2)
  delta3 <- (1-status1)*status2
  delta4 <- status1*status2

  new.term1 <- new.term2 <- new.term3 <- new.term4 <- rep(0,length(u1))
  if(length(which(delta1==1)) >0) {new.term1[delta1==1] <- log(term1[delta1==1])}
  if(length(which(delta2==1)) >0) {new.term2[delta2==1] <- log(term2[delta2==1])}
  if(length(which(delta3==1)) >0) {new.term3[delta3==1] <- log(term3[delta3==1])}
  if(length(which(delta4==1)) >0) {new.term4[delta4==1] <- log(term4[delta4==1])}

  get.terms <- (new.term1 + new.term2+ new.term3 + new.term4)*W
  if(control.finite==TRUE){ loglik <- sum(get.terms[is.finite(get.terms)]) }
  if(control.finite==FALSE){loglik <- sum(get.terms) }

  if(negative==TRUE){return(-loglik)
  }else{return(loglik)}
}

