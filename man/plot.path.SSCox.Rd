\name{plot.path.SSCox}
\alias{plot.path.SSCox}

\title{Plotting the Sparse-smothed Cox PH path.}

\description{Plotting a Sparse-smothed Cox PH  fitted resulting from \code{\link{cv.SSCox}}.}

\usage{plotpath(x, labelsize = 0.6, ...) }

\arguments{
\item{x}{This object can be created using the \code{\link{cv.SSCox}} function using different values of \code{theta}.}

\item{labelsize}{Sets the size of the variable labels in the plot. Set to zero for no variable labels.}
\item{...}{Any other arguments will be forwarded to the plot function.}
}

\value{
A plot is created
}
\references{
	Khaled Alqahtani, Henry M.wood, Charles Taylor, and Arief Gusnanto (2016) Sparse-smoothed Cox proportional hazard model for genomic profiles
}
\author{
Khaled Alqahtani
}

\seealso{
	\code{\link{SSCox}},\code{\link{cv.SSCox}}
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

#Plotting the Sparse-smothed Cox PH path
plot.path.SSCox(cv)

}