#' @title The Beran Estimator
#'
#' @description This function defines the Beran estimator to estimate conditional marginal survival functions for right censored event time data given a covariate.
#'
#' @param t specific time value at which the conditional survival function is calculated
#' @param x specific covariate value given which the conditional survival function is calculated
#' @param Y vector of event times
#' @param status vector of censoring indicators taking values 0 (censored) or 1 (uncensored)
#' @param X vector of covariate values
#' @param band bandwidth parameter
#' @param Kern Kernel function
#' @keywords internal
Beran = function(t, x, Y, status, X, band, Kern=KernEpa){
  if(t < min(Y)){term3 = 1}
  else{ if(t>=max(Y)){term3=0}
    else{
      W <-  Kern((X-x)/band) / sum(Kern((X-x)/band))
      keep <- which(Y <= t)

      Y.keep <- Y[keep]
      W.keep <- W[keep]
      status.keep <- status[keep]
      n.keep <- length(Y.keep)

      term1 <- rep(NA,nrow=n.keep)
      term2 <- rep(NA,nrow=n.keep)

      for(i in 1:n.keep){
        term1[i] <- sum(as.numeric(Y>=Y.keep[i])*W)

        if(W.keep[i]!=0){term2[i] <- (1- W.keep[i]/term1[i])^status.keep[i]}
        else{term2[i] <- 1}
      }

      term3 <- prod(term2)
    }}
  return(1-term3)
}
