\name{cv.SScox}
\alias{cv.SScox}
\alias{lambdamaxrgeomn}
\alias{logPLfixr}


\title{Cross-validation in Sparse-smothed Cox PH model.}

\description{Cross-validating Sparse-smothed Cox PH, using penalized likelihood cross-validation.}

\usage{

cv.SScox (time,event,formula,CNA, fold = 5,alphac.vec,alphan.vec,diff=2,steps,maxiter=400,eps=1e-6,print=TRUE ) 


}

\arguments{

 \item{time}{Survival time}
  \item{event}{The status indicator, normally 0=alive, 1=dead.}
  \item{formula}{a formula object of patients clinical data (fix effects) after a ~ operator.}
  \item{CNA}{Matrix of copy number alterations. No missing value is allowed.}
\item{fold}{The fold for cross-validation. May be supplied as a single number (between 3 and n) . Although \code{fold} can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. The default is \code{fold = 5}, resulting in  (5-fold) cross-validation. }

\item{alphac.vec, alphan.vec}{The fixed wights of the tuning parameters for Normal and  Cauchy penalization; where the wight of Laplace is (alphal.vec= 1-(alphac.vec+ alphan.vec)). Each must be either a single positive numbers or a vector where the sum of each sets of wights equal to 1.}
\item{diff}{the order differences of random effects of Cauchy distribution. The default is \code{diff =2}, resulting in a second-order differences of random effects. }

\item{steps}{The maximum number of steps between \code{mintheta} and \code{maxtheta} at which the cross-validated likelihood is to be calculated. \code{mintheta} and  \code{maxtheta} is the output of \code{lambdamaxrgeomn} which depends on the wight of Laplace.}

\item{maxiter}{The maximum number of iterations allowed in each fitting of the model. Set by default at 400 .}

\item{eps}{The convergence criterion.}
\item{print}{If print is TRUE, each steps of cross validation are printed.}

  }

\value{ An object of class \code{"cv.SSCox"} is returned, which is a
  list of the cross-validation fit.
	\item{cvmat}{Matrix of cross-validated partial likelihood.}
	\item{opt}{The optimal tuning parameters which maximize cross-validated partial likelihood.}
	\item{sd}{Standard errors for the cross-validated partial likelihood curve at each tuning parameter.}
	\item{range.theta}{The values of \code{theta} used in the fits.}
	\item{betar}{The estimation of random effects on the full data at each tuning parameter to be used in \code{plot.path.SSCox}.}
	
}

\references{
	Khaled Alqahtani, Henry M.wood, Charles Taylor, and Arief Gusnanto (2016) Sparse-smoothed Cox proportional hazard model for genomic profiles.
}
\author{
Khaled Alqahtani
}

\seealso{
	\code{\link{plot.cv.SSCox}}, and \code{\link{plot.path.SSCox}}
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

# Finding an optimal cross-validated partial likelihood

cv <- cv.SScox(Survival,Status,~Age+Sex,CNA,alphac.vec=0.2,alphan.vec=0.3,steps=10)

}