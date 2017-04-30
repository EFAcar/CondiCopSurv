#' @title Simulating Event Time Data for GLRT Bootstrap Samples with Nonparametric Margins
#' @description This function simulates right censored event time data for boostrap samples when conditional margins are estimated nonparametrically using the Beran's estimator.
#'
#' @param Y1,Y2 vectors of event times
#' @param status1,status2 vectors of censoring indicators
#' @param X vector of covariate values
#' @param family an integer defining the bivariate copula family:
#'
#' 3 = Clayton copula
#'
#' 4 = Gumbel copula
#'
#' 5 = Frank copula
#' @param cop.par copula parameter
#' @param mar.band vector of bandwidth values for the estimation of conditional margins
#' @param Kern Kernel function
#' @importFrom VineCopula BiCopSim
#' @keywords internal
BootSim.NP = function(Y1, Y2, status1, status2, X, family, cop.par, mar.band, Kern=KernEpa){
  
  if(!(family %in% c(3, 4, 5)))
    stop("Copula family not yet implemented.")
  
  if(length(mar.band) != 2)
    stop("Provide a two-dimensional bandwidth vector for the conditional margins.")
  
  n = length(X)
  Udata = CDVine::BiCopSim(n, family, cop.par)
  U.temp1 = Udata[,1]; U.temp2 = Udata[,2]
  
  ET1 = sapply(1:n, function(i) BeranInv(u=U.temp1[i], x=X[i], Y=Y1, status=status1, X=X, band=mar.band[1], Kern=Kern))
  ET2 = sapply(1:n, function(i) BeranInv(u=U.temp2[i], x=X[i], Y=Y2, status=status2, X=X, band=mar.band[2], Kern=Kern))
  
  CenT = BootSim.Censor(Y1, Y2, status1, status2)
  
  Y1 = pmin(ET1, CenT)
  Y2 = pmin(ET2, CenT)
  status1 = 1*(ET1 <= CenT)
  status2 = 1*(ET2 <= CenT)
  X = X
  
  Boot.data = data.frame(Y1,Y2, status1, status2, X)
  return(Boot.data)
}
