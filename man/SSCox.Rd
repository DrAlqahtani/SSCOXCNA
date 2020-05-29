\name{SSCox}
\alias{SSCox}
\title{
Sparse-smothed Cox proportional hazard model on clinical data and copy number alteration profiles.

}
\description{
	Sparse-smothed Cox proportional hazard model on clinical data and copy number alteration profiles. 
}
\usage{
SSCox(time,event,formula,CNA,init.fix=matrix(0,ncol=1,nrow=dim(model.matrix(formula))[2]-1),init.random=matrix(0,ncol=1,nrow=ncol(CNA)),theta,alphan,alphac,diff=2,eps=1e-6,maxiter=500)}
\arguments{
  \item{time}{
 Survival time.
}
  \item{event}{
The status indicator, normally 0=alive, 1=dead.
}
  \item{formula}{
A formula object of patients clinical data (fix effects) after a ~ operator.
}
  \item{CNA}{
Matrix of copy number alterations. No missing value is allowed.
}
  \item{init.fix}{
vector of init.fixal values of the iteration for patients' clinical variables (Fixed effects' estimate)  . Default init.fixa value is zero for all variables.
}
  \item{init.random}{
  vector of init.fixal values of the iteration for copy number alterations variables (Random effects' estimate)  . Default init.random value is zero for all variables.	
}
  \item{theta}{
  Value of \code{thetas} to be evaluated
}

\item{alphan ,alphac}{The fixed wights of the tuning parameters for Normal and  Cauchy penalization; where the wight of Laplace is (alphal= 1-(alphac+alphan)). Each must be  a single positive numbers where the sum of each sets of wights equal to 1.}

\item{diff}{The order differences of random effects of Cauchy distribution. The default is \code{diff =2}, resulting in a second-order differences of random effects. }

  \item{eps}{
A termination value for the convergence criterion, by default set to be 1e-6
}
  \item{maxiter}{
The maximum number of iterations; by default set to be 500
}
}
\details{
	formula for clinical data, lets say x1,x2,x3, should be written as ~x1+x2+x3
}
\value{
	A list containing
\item{beta}{Estimates of the fixed effects (from clinical data).}
\item{b}{Estimates of the random effects (from CNA profiles).}
\item{nvar}{The total number of nonzero estimates of random effects .}
\item{ll}{Object of class "numeric". Returns the log partial likelihood of the fitted model.}
\item{iterations}{Object of class "numeric". Number of iterations used in the fitting process.}
\item{converged}{Object of class "logical". Whether the fitting process was judged to be converged.}
\item{varfix}{Variance matrix of the estimates of the fixed effects}
\item{H0}{The baseline cumulative hazard rate.}
\item{time}{The observed survival time.}
\item{status}{The status indicator, normally 0=alive, 1=dead (need to be used to check the residuals).}
\item{X}{The design matrix of fixed effects (need to be used to check the residuals).}
\item{Z}{The design matrix of random effects (need to be used to check the residuals).}
\item{Chr}{Chromosomes.}
\item{pos}{position of each chromosomes.}
\item{thetavalue}{Evaluated \code{theta}.}

}
\references{
	Khaled Alqahtani, Henry M.wood, Charles Taylor, and Arief Gusnanto (2016) Sparse-smoothed Cox proportional hazard model for genomic profiles.
}
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
}

\keyword{ Sparse-smothed Cox proportional hazard model  }

