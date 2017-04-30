#' @title Initial Values for Optimization
#'
#' @description These functions determine initial parameter values for the (local) likelihood estimation.
#' @keywords internal
Null.init = function(family){
  if(family==3 || family ==4){ intv  <-  c(-10,5)}
  if(family==5){intv  <-  c(-10,15)}				
  return(intv)
}
Find.init  =  function(grid.val, u1, u2, status1, status2, X, x, family, p, band, Kern=KernEpa){
  easy.fnc  =  function(eta){
    CondiCopSurv.LocLLik(eta, family=family, p=p, u1=u1, u2=u2, status1=status1, status2=status2, X=X, x=x, band=band, Kern=Kern, control.finite=TRUE, negative=FALSE)
  }
  locllik  <-  try(apply(grid.val,1, easy.fnc), silent=T)
  return(as.numeric(grid.val[which.max(locllik),]))
}


