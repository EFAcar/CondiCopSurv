#' Epanechnikov Kernel
#'
#' This function defines the Epanechnikov kernel used in the local likelihood estimation.
#'
#' @param t numeric vector
#' @keywords internal
KernEpa <-function(t){
  if(abs(t)<1){return(3/4*(1-t^2))}
  else {return(0)}  }
KernEpa=Vectorize(KernEpa)
