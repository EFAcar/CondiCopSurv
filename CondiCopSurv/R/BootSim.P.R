#' @title Simulating Event Time Data for GLRT Bootstrap Samples with Parametric Margins
#' @description This function simulates right censored event time data for boostrap samples when conditional margins are estimated parametrically using a Weibull model.
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
#' @param mar.par vector of marginal parameters under the Weibull model. The entries are of the form: (rho1, rho2, lambda1, lambda2, b1, b2)
#' @importFrom VineCopula BiCopSim
#' @keywords internal
BootSim.P = function(Y1, Y2, status1, status2, X, family, cop.par, mar.par){
  
  if(!(family %in% c(3, 4, 5)))
    stop("Copula family not yet implemented.")
  
  if(length(mar.par) != 6)
    stop("Provide marginal parameters (rho1, rho2, lambda1, lambda2, b1, b2) for the conditional marginal survival functions under the Weibull model")
  
  n = length(X)
  Udata = CDVine::BiCopSim(n, family, cop.par)
  U.temp1 = Udata[,1]; U.temp2 = Udata[,2]
  
  rho1 = mar.par[1]    ; rho2 = mar.par[2]
  lambda1 = mar.par[3] ; lambda2 = mar.par[4]
  b1 = mar.par[5]      ; b2 = mar.par[6]
  
  shape1 = rho1; shape2 = rho2;
  scale1 = (lambda1 * exp(b1*X))^(-1/rho1)
  scale2 = (lambda2 * exp(b2*X))^(-1/rho2)
  
  ET1 = qweibull(U.temp1, shape = shape1, scale = scale1, lower.tail= F)
  ET2 = qweibull(U.temp2, shape = shape2, scale = scale2, lower.tail= F)
  
  CenT = BootSim.Censor(Y1, Y2, status1, status2)
  
  Y1 = pmin(ET1, CenT)
  Y2 = pmin(ET2, CenT)
  status1 = 1*(ET1 <= CenT)
  status2 = 1*(ET2 <= CenT)
  X = X
  
  Boot.data = data.frame(Y1,Y2, status1, status2, X)
  return(Boot.data)
}