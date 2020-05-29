\name{SSCox-package}
\alias{SSCox-package}
\alias{SSCox}
\docType{package}
\title{
Sparse-smoothed Cox proportional hazard model for genomic profiles.
}
\description{
Functions to model Sparse-smoothed Cox proportional hazard model for genomic profiles. Moreover, the packages contains related plotting functions.
}
\details{
\tabular{ll}{
Package: \tab SSCox\cr
Type: \tab Package\cr
Version: \tab 1.0 \cr
Date: \tab 2016-07-03\cr
License: \tab GPL (>=2)\cr
}
library(SSCox)
}
\author{
Khaled Alqahtani

Maintainer: Khaled Alqahtani <mmkma@leeds.ac.uk>

}
\references{
Khaled Alqahtani, Henry M.wood, Charles Taylor, and Arief Gusnanto (2016) Sparse-smoothed Cox proportional hazard model for genomic profiles.
}

\keyword{ package }
\seealso{
 \code{\link{pre.SSCox}}, \code{\link{SSCox}}, \code{\link{summary.SSCox}}, \code{\link{plot.brandom}}, \code{\link{plot.Coxsnell}}, \code{\link{survpredict.SSCox}}, \code{\link{cv.SSCox}}, \code{\link{plot.cv.SSCox }}, and \code{\link{plot.path.SSCox}}
}
\examples{
	library(SSCox)
	# No real data example yet.
# Examples using simulated data
set.seed(124)

#clinical data
Age=sample(50:75,80,replace=TRUE)
Sex=c(rep(1,35), rep(0,45))
Survival= rexp(80,3)
Status=rbinom(80,1,0.8)

# Simulated Copy number alteations profile 
CNA <- matrix(rnorm(80*200),80)
CNAh20=matrix(nrow=length(Survival),ncol=20)
CNAh20[which(Survival>median(Survival)),1:10]=rnorm(400,6,1)
CNAh20[which(Survival<=median(Survival)),1:10]=rnorm(400,2,1)
CNAh20[which(Survival>median(Survival)),11:20]=rnorm(400,2,1)
CNAh20[which(Survival<=median(Survival)),11:20]=rnorm(400,6,1)
CNAg80=matrix(rnorm(6400,4,1),nrow=80,ncol=80)
CNA=cbind(CNAh20,CNAg80)
attr(CNA,"Chr") <- c(paste("chr",rep(1:19,rep(5,19)),sep=""),rep("chrX",5))
attr(CNA,"Pos") <- c(rep(1:5,19),1:5)

 # Sparse-smothed Cox PH
fit <- SSCox(Survival, Status,~Age+Sex,CNA,theta=1e-3,alphan=0.3,alphac=0.2)

#summary of Sparse-smothed Cox PH
summary.SSCox(fit)

# Plot random effects estimates
plot.brandom(fit, type="h")

# Plot Cox-Snell residuals
plot.Coxsnell(fit)

# (Simulated) New data for prediction
new.cna <- matrix(rnorm(100),1)
new.clin <- data.frame(Age=65, Sex=1) 

# Prediction on the new data
survpredict.SSCox(fit,new.clin=new.clin,new.cna=new.cna)

###

# Finding an optimal cross-validated partial likelihood
cv <- cv.SScox(Survival,Status,~Age+Sex,CNA,alphac.vec=0.2,alphan.vec=0.3,steps=10)

# plotting cross-validated partial likelihood curve
plot.cv.SSCox(cv)

# plotting the path 
plot.path.SSCox(cv)
}
