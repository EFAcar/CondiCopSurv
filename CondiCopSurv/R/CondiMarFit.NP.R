#' @title Fitting the Conditional Survival Functions
#'
#' @description This function fits conditional marginal survival functions for right censored survival data.
#'
#' @param Y vector of event times
#' @param status vector of censoring indicators taking values 0 (censored) or 1 (uncensored)
#' @param X vector of covariate values
#' @param band bandwidth value used in the Beran estimator (optional)
#' @param pilot user-specified pilot bandwidth values (optional)
#' @param nband desired length of bandwidth values if pilot values are not specified
#' @param band.method method used in bandwidth selection: three cross-validation based selectors 
#' @keywords internal
CondiMarFit.NP = function(Y, status, X, band, Kern=KernEpa, 
                          pilot=NULL, nband=10, 
                          band.method= c("CV_Y","CV_T","CV_use"), ...){
  
  band.method <- match.arg(band.method)
  
  if(missing(band)){
    if(is.null(pilot)){
      nband <- nband +2 ; diffs <- diff(sort(X),1)
      h.min <- max(diffs); h.max <- max(X)-min(X)
      log.seq <- seq(from=log(h.min),to=log(h.max), length.out =nband)
      pilot <- round(exp(log.seq),2)
      pilot <- pilot[-c(1:2)]
      nband <- length(pilot)
    }else{
      pilot <-pilot
      nband <-length(pilot)
    }   
    
    
    # Based on http://smj.sagepub.com.uml.idm.oclc.org/content/7/4/329.full.pdf+html
    
    if(band.method=="CV_Y"){
      dev.sum = rep(NA, nband)
      for(k in 1:nband){
        dev <- rep(NA, length(Y))
        for(i in 1:length(Y)){
          dev[i] <- sum(sapply(1:length(Y), function(j){
            (1*(Y[i] <= Y[j]) - Beran(t=Y[j],x=X[i], Y[-i], status[-i], X[-i], 
                                      band=pilot[k],Kern=Kern))^2}))
        }
        dev.sum[k] <- sum(dev)
      }
      band <- pilot[which.min(dev.sum)]
    }
    
    if(band.method=="CV_T"){
      dev.sum <- rep(NA, nband)
      for(k in 1:nband){
        dev <- rep(NA, length(Y))
        for(i in 1:length(Y)){
          dev[i] <- sum(sapply(1:length(Y), function(j){
            status[i]*status[j]*(1*(Y[i] <= Y[j]) - Beran(t=Y[j],x=X[i], Y[-i], status[-i], X[-i], 
                                                          band=pilot[k],Kern=Kern))^2}))
        }
        dev.sum[k] <- sum(dev)
      }
      band <- pilot[which.min(dev.sum)]
    }
    
    if(band.method=="CV_use"){
      CV_ind=matrix(0,nrow=length(Y),ncol=length(Y))
      for(i in 1:length(Y)){
        for(j in 1:length(Y)){
          if(status[i]==1 & status[j]==1){CV_ind[i,j] <- 1}
          if(status[i]==1 & status[j]==0 & Y[i]<=Y[j]){CV_ind[i,j] <- 1}
          if(status[i]==0 & status[j]==1 & Y[i]>=Y[j]){CV_ind[i,j] <- 1}
        }
        CV_ind[i,i] <- 1
      }
      
      dev.sum <- rep(NA, nband)
      for(k in 1:nband){
        dev <- rep(NA, length(Y))
        for(i in 1:length(Y)){
          dev[i] <- sum(sapply(1:length(Y), function(j){
            CV_ind[i,j]*(1*(Y[i] <= Y[j]) - Beran(t=Y[j],x=X[i], Y[-i], status[-i], X[-i], 
                                                  band=pilot[k],Kern=Kern))^2}))
        }
        dev.sum[k] <- sum(dev)
      }
      band <- pilot[which.min(dev.sum)]
    }
  }
  
  u <- sapply(1:length(Y), function(j) 1-Beran(t=Y[j], x=X[j], Y, status, X, band=band, Kern=KernEpa))
  return(list(u=u, band=band))
}

