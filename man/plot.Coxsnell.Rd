\name{plot.Coxsnell}
\title{
 plot of Cox-Snell residuals.
}
\description{
	Plot of  Cox-snell residuals}
\usage{
plot.Coxsnell(x, xlab = NULL, ylab = NULL, ...)
}
\arguments{
  \item{x}{
Output of the Sparse-smothed Cox proportional hazard model \code{SSCox}.
}
  \item{xlab}{
Horizontal axis label.
}
  \item{ylab}{
Vertical axis label.
}
  \item{\dots}{
Further arguments to be passed on to \code{plot}.
}
}

\value{
A plot is created.
}
\references{
	Khaled Alqahtani, Henry M.wood, Charles Taylor, and Arief Gusnanto (2016) Sparse-smoothed Cox proportional hazard model for genomic profiles
}
\author{
Khaled Alqahtani
}

\seealso{
	\code{\link{SSCox}}
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

# Plot Cox-Snell residuals
plot.Coxsnell(fit)

}

\keyword{ plot }


