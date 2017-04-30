#' @title Inverse of the Beran Estimator
#'
#' @description This function calculates the inverse of the Beran estimator required to generate bootstrap samples for the generalized likelihood ratio test.
#'
#' @param u specific value in [0,1] at which the inverse of the conditional survival function is calculated
#' @param x specific covariate value given which the inverse of the conditional survival function is calculated
#' @param Y vector of event times
#' @param status vector of censoring indicators taking values 0 (censored) or 1 (uncensored)
#' @param X vector of covariate values
#' @param band bandwidth parameter
#' @param Kern Kernel function
#' @keywords internal
BeranInv = function(u, x,  Y, status,  X, band, Kern=KernEpa){

n <- length(Y)

Shat <- sapply(1:n, function(j) 1-Beran(t=Y[j], x=x, Y=Y, status=status, X=X,  band=band, Kern=Kern))

if(u==1){invBeran <- 0}
else{ if (u < min(Shat)){invBeran <- Inf}
 else{invBeran <- min(Y[Shat <= u])}}

return(invBeran)
}

