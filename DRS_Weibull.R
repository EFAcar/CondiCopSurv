################################################################################
# Analysis of Diabetic Retinopathy Data 
################################################################################
# Geerdens, Acar and Janssen (2017)
# Conditional copula models for right-censored clustered event time data
################################################################################

# load("DRS.Rdata")

################################################################################
# Parametric Estimation of the Conditional Margins under the Weibull model
################################################################################

# require(survival)

fit1 = survreg(Surv(Y1, status1)~X, data=work.data, dist="weibull")
fit2 = survreg(Surv(Y2, status2)~X, data=work.data, dist="weibull")

# Parameter Estimates 

hat.rho1 = 1/(fit1$scale) # 0.788
hat.rho2 = 1/(fit2$scale) # 0.830
hat.lambda1 = exp(-fit1$coefficients[1]/fit1$scale) # 0.0214
hat.lambda2 = exp(-fit2$coefficients[1]/fit2$scale) # 0.0219
hat.b1 = -fit1$coefficients[2]/fit1$scale  # -0.015
hat.b2 = -fit2$coefficients[2]/fit2$scale  #  0.014

# Standard Errors (Margin 1)

varscale = fit1$var[3,3]*fit1$scale^2 
covscale = fit1$var[1,3]*fit1$scale
var.hat.lambda1 = exp(-2*fit1$coef[1]/fit1$scale)*
  (fit1$var[1,1]/fit1$scale^2+fit1$coef[1]^2*varscale/fit1$scale^4 - 
     2*fit1$coef[1]*covscale/fit1$scale^3)
var.hat.rho1 = varscale/fit1$scale^4
covscale = fit2$var[2,3]*fit2$scale
var.hat.b1=diag(fit1$var)[2]/fit1$scale^2+ fit1$coef[2]^2*varscale/fit1$scale^4-
  2*fit1$coef[2]*covscale/fit1$scale^3

sd.Mar1 = rbind(rho=c(hat.rho1,sqrt(var.hat.rho1)), 
                lambda=c(hat.lambda1,sqrt(var.hat.lambda1)), 
                beta=c(hat.b1,sqrt(var.hat.b1)))
colnames(sd.Mar1)=c('parameter estimate','standard error')
print(round(sd.Mar1,3))


# Standard Errors (Margin 2)

varscale=fit2$var[3,3]*fit2$scale^2
covscale=fit2$var[1,3]*fit2$scale
var.hat.lambda2 = exp(-2*fit2$coef[1]/fit2$scale)*
  (fit2$var[1,1]/fit2$scale^2+fit2$coef[1]^2*varscale/fit2$scale^4-
     2*fit2$coef[1]*covscale/fit2$scale^3)
var.hat.rho2 = varscale/fit2$scale^4
covscale = fit2$var[2,3]*fit2$scale
var.hat.b2 = diag(fit2$var)[2]/fit2$scale^2+ fit2$coef[2]^2*varscale/fit2$scale^4-
  2*fit2$coef[2]*covscale/fit2$scale^3

sd.Mar2 = rbind(rho=c(hat.rho2,sqrt(var.hat.rho2)), 
                lambda=c(hat.lambda2,sqrt(var.hat.lambda2)), 
                beta=c(hat.b2,sqrt(var.hat.b2)))
colnames(sd.Mar2)=c('parameter estimate','standard error')
print(round(sd.Mar2,3))


# Copula Data - Parametric Margins

PMar.est = c(hat.rho1, hat.rho2, hat.lambda1 , hat.lambda2, hat.b1, hat.b2)
P.U1  = exp(-hat.lambda1* work.data$Y1^hat.rho1 * exp(hat.b1*work.data$X))
P.U2  = exp(-hat.lambda2* work.data$Y2^hat.rho2 * exp(hat.b2*work.data$X))
U.data1 = data.frame(P.U1, P.U2, work.data$status1, work.data$status2)



