\name{survpredict.SSCox}
\alias{survpredict.SSCox}
\title{
Prediction of Survival probability based on new data. Also, It produce the survival prediction plot 
}
\description{
Prediction of Survival probability based on new data. Also, It produce the survival prediction plot 
}
\usage{
survpredict.SSCox(x, new.clin = 1, new.cna = NULL, plot = T)
}

\arguments{
  \item{x}{
Output from \code{\link{SSCox}}
}
  \item{new.clin}{
New clinical data (in data frame or matrix), whose number (and order) of columns must match the number of columns of the matrix  in the fitted model.
}
  \item{new.cna}{
New CNA profiles (in data frame or matrix), whose number (and order) of columns must match the number of columns of the matrix  in the fitted model.
}
  \item{plot}{
Should a survival plot be produced
}
}

\value{
\item{surv.est}{A predicted Survival probability}
}
\references{
Khaled Alqahtani, Henry M.wood, Charles Taylor, and Arief Gusnanto (2016) Sparse-smoothed Cox proportional hazard model for genomic profiles.}
\author{
Khaled Alqahtani
}
\examples{
	# No real data example yet.
# Examples using simulated data
set.seed(124)

#clinical data
Age=sample(50:75,80,replace=TRUE)
Sex=c(rep(1,35), rep(0,45))
Survival= rexp(80,3)
Status=rbinom(80,1,0.8)

# Simulated Copy number alterations profile 
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

# (Simulated) New data for prediction
new.cna <- matrix(rnorm(100),1)
new.clin <- data.frame(Age=65, Sex=1) 

# Prediction on the new data
survpredict.SSCox(fit,new.clin=new.clin,new.cna=new.cna)
}
\keyword{ prediction}

