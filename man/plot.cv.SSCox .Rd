\name{plot.cv.SSCox}
\alias{plot.cv.SSCox }
\alias{error.bars}
\title{Plotting the cross-validated likelihood versus log of theta}

\description{Plotting the cross-validated likelihood versus  log of theta resulting from \code{\link{cv.SSCox}}.}

\usage{plot.cv.SSCox(x,...) }

\arguments{
\item{x}{This object can be created using the \code{\link{cv.SSCox}} function using different values of theta.}

\item{...}{Any other arguments will be forwarded to the plot function.}
}

\value{
A cross-validated likelihood plot is created.
}
\references{
	Khaled Alqahtani, Henry M.wood, Charles Taylor, and Arief Gusnanto (2016) Sparse-smoothed Cox proportional hazard model for genomic profiles
}
\author{
Khaled Alqahtani
}

\seealso{
	\code{\link{SSCox}}, and \code{\link{cv.SSCox}}
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

#Finding an optimal cross-validated partial likelihood

cv <- cv.SScox(Survival,Status,~Age+Sex,CNA,alphac.vec=0.2,alphan.vec=0.3,steps=10)

# plotting cross-validated partial likelihood curve
plot.cv.SSCox(cv)

}