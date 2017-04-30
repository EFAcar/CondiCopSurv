#' Conditional Copula Models for Right Censored Event Time Data
#'
#' @name CondiCopSurv-package
#' @docType package
#' @examples
#' #--- An Illustration of the Analysis with Simulated Data  -----------------------------
#' 
#' 
#' library(CondiCopSurv)
#' VecBiCopPar2Tau = Vectorize(BiCopPar2Tau)
#' VecBiCopTau2Par = Vectorize(BiCopTau2Par)
#' 
#' myfamily = 3
#' n=250
#'
#' # marginal parameters under the Weibull model
#' lambda1 = 0.5 ; rho1 = 1.5 ; b1 = 0.8
#' lambda2 = 0.5;  rho2 = 1.5 ; b2 = 0.8
#'
#' # parameters for the censoring time distribution
#' lambdac = 1.5;  rhoc = 0.5         # Moderate Censoring: about 50 %
#'
#' # Kendall's Tau function
#' tau.fnc = function(t){ -0.1*(t-3)^2 +0.7}
#'
#' #--- Data Generation -------------------------------------------------------
#'
#' # Generate data from Conditional Copula
#' X = runif(n,2,5)
#' true.tau = tau.fnc(X)
#' true.theta = VecBiCopTau2Par(myfamily,true.tau)
#' true.eta = log(true.theta)
#' Udata = t(mapply(BiCopSim, 1, family=myfamily, par=true.theta))
#'
#'
#' # Obtain event times under the Weibull model (reparametrization is required)
#' shape1 = rho1; shape2 = rho2; shapec = rhoc
#' scale1 = (lambda1 * exp(b1*X))^(-1/rho1)
#' scale2 = (lambda2 * exp(b2*X))^(-1/rho2)
#' scalec = (lambdac)^(-1/rhoc)
#'
#' U.temp1 = Udata[,1]
#' U.temp2 = Udata[,2]
#'
#' ET1 = qweibull(U.temp1, shape = shape1, scale = scale1, lower.tail= F)
#' ET2 = qweibull(U.temp2, shape = shape2, scale = scale2, lower.tail= F)
#' CenT = rweibull(n, shape = shapec, scale = scalec)
#'
#' Y1 = pmin(ET1, CenT)
#' Y2 = pmin(ET2, CenT)
#' status1 = 1*(ET1 <= CenT)
#' status2 = 1*(ET2 <= CenT)
#' Obs.data = data.frame(Y1,Y2, status1, status2, X)
#'
#'
#' #--- Estimation ------------------------------------------------------------
#'
#' # Parametric Conditional Margins
#' fit1 = survreg(Surv(Y1, status1)~X, dist="weibull")
#' fit2 = survreg(Surv(Y2, status2)~X, dist="weibull")
#'
#' hat.rho1 = 1/(fit1$scale) ; hat.rho2 = 1/(fit2$scale)
#' hat.lambda1 = exp(-fit1$coefficients[1]/fit1$scale)
#' hat.lambda2 = exp(-fit2$coefficients[1]/fit2$scale)
#' hat.b1 = -fit1$coefficients[2]/fit1$scale
#' hat.b2 = -fit2$coefficients[2]/fit2$scale
#'
#' PMar.est = c(hat.rho1, hat.rho2, hat.lambda1,
#'   hat.lambda2, hat.b1, hat.b2)
#'
#' P.U1  = exp(-hat.lambda1* Y1^hat.rho1 * exp(hat.b1*X))
#' P.U2  = exp(-hat.lambda2* Y2^hat.rho2 * exp(hat.b2*X))
#'
#' U.data1 = data.frame(P.U1, P.U2, status1, status2, X)
#' hopt.1 = 1
#'
#'
#' # NonParametric Conditional Margins
#'
#' ### Joint Bandwidth Selection for NonParametric Conditional Margins and Conditional Copula
#' ## pilot.val = c(0.3, 0.5, 0.8, 1.2, 2, 3)
#' pilot.val = 2
#' cross.band = expand.grid(pilot.val,pilot.val,pilot.val)
#' system.time({
#'   band.NPMar = BandSelect.NPMar(pilot.mar = pilot.val,
#'     pilot.cop= pilot.val, Y1=Y1, Y2=Y2,
#'     status1=status1, status2=status2, X=X,
#'     family=myfamily, p=1, Kern=KernEpa)
#'   hopt.2 = band.NPMar$h
#'   U.data2 = band.NPMar$data
#' })
#'
#' # Local Likelihood Estimation
#' x = seq(2, 5, by=0.1)
#' nx = length(x)
#' fit.PMar <- CondiCopSurvFit.PMar(Y1=Y1, Y2=Y2,
#'                                  status1=status1, status2=status2,
#'                                  X=X, x=x, family=myfamily, p=1,
#'                                  band.cop=hopt.1)
#' fit.NPMar <- CondiCopSurvFit.NPMar(Y1=Y1, Y2=Y2,
#'                                    status1=status1, status2=status2,
#'                                    X=X, x=x, family=myfamily, p=1,
#'                                    band.cop=hopt.2[3],
#'                                    band.mar=hopt.2[1:2])
#'
#' # Plot Parametric and NonParametric Estimates
#' par(mar = c(4,4.2,.1,.1)+.1)
#' plot(x, BiCopPar2Tau(myfamily, fit.PMar$cop.par), ylim=c(0,1),
#'      type="b", col=4,
#'      xlab = expression(x), ylab = expression(tau(x)))
#' points(x, BiCopPar2Tau(myfamily,fit.NPMar$cop.par), type="b", col=3)
#' legend("topright", legend = c("Parametric", "Non-Parametric"),
#'        pch = 21, col = c(4,3), title = "Margins")
NULL