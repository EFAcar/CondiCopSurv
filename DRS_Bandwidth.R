################################################################################
# Analysis of Diabetic Retinopathy Data 
################################################################################
# Geerdens, Acar and Janssen (2017)
# Conditional copula models for right-censored clustered event time data
################################################################################

# load("DRS.Rdata")

################################################################################
# Bandwidth Selection 
################################################################################


# Bandwidth Selection for Parametric Margins
#############################################
band.Clayton.Par = BandSelect(pilot=pilot.val, u1=U.data1[,1], u2=U.data1[,2], 
                              status1=status1, status2=status2, X=X, family=3)

# check the bandwidth values yielding the highest three CV-loglikelihood:

cbind(band.Clayton.Par$h.pilot[sort(band.Clayton.Par$cvLIK, index.return=T, decreasing=T)$ix], 
      sort(band.Clayton.Par$cvLIK, decreasing=T))

plot(band.Clayton.Par$h.pilot[-c(1:3)],band.Clayton.Par$cvLIK[-c(1:3)], 
     xlab="h", ylab="CV-logL", col=2 , type="l")






# Bandwidth Selection for NonParametric Margins
#################################################

band.Clayton.NPar = BandSelect(pilot= pilot.val, u1=U.data2[,1], u2=U.data2[,2], 
                               status1=status1, status2=status2, X=X, family=3)


# check the bandwidth values yielding the highest three CV-loglikelihood:

cbind(band.Clayton.NPar$h.pilot[sort(band.Clayton.NPar$cvLIK, index.return=T, decreasing=T)$ix], 
      sort(band.Clayton.NPar$cvLIK, decreasing=T))

plot(band.Clayton.NPar$h.pilot[-c(1:3)],band.Clayton.NPar$cvLIK[-c(1:3)], 
     xlab="h", ylab="CV-logL", col=2 , type="l")


##################################################################
### Notes: 
##################################################################
## A small bandwidth choice may lead to highly unstable estimates.
## In such cases, we suggest choosing the second best, or anoother reasonable, bandwidth value. 

band.Clayton.Par$h <- 42;  band.Clayton.NPar$h <- 42
#band.Frank.Par$h  <- 23;  band.Frank.NPar$h <-  23
#band.Gumbel.Par$h <- 31;  band.Gumbel.NPar$h <- 42

