################################################################################
# Analysis of Diabetic Retinopathy Data 
################################################################################
# Geerdens, Acar and Janssen (2017)
# Conditional copula models for right-censored clustered event time data
################################################################################

# load("DRS.Rdata")
# library(CondiCopSurv)


################################################################################
# Nonparametric Estimation of the Conditional Margins using the Beran's Estimator
################################################################################

# Bandwidth Selection

h1.result = CondiMarFit(Y= Y1, status = status1, X=X, mar.method="Beran",
                        Kern=KernEpa, pilot=pilot.val, band.method= "CV_use")
h1 = h1.result$band  # 17


h2.result = CondiMarFit.NP(Y= Y2, status = status2, X=X, mar.method="Beran",
                           Kern=KernEpa, pilot=pilot.val, band.method= "CV_use") 
h2 = h2.result$band  # 17


# Copula Data - Nonparametric Margins

NP.U1= h1.result$u
NP.U2= h2.result$u
U.data2 = data.frame(NP.U1, NP.U2, work.data$status1, work.data$status2)




