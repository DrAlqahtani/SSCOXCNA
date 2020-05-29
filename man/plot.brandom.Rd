\name{plot.brandom}
\alias{plot.brandom}
\alias{plot.loc}
\title{
 plot random effects estimates from \code{SSCox} output
}
\description{
	Plot random effects estimates from \code{SSCox} function. \code{plot.loc} is the main function to create the plot in \code{plot.brandom}.
}
\usage{
plot.brandom(x, Chromosome = NULL, xlab = "Genomic Regions", ylab = "", cex.axis = 0.8, ...) 
plot.loc(x,Chr=NULL, xlab="Genomic Regions", ylab="",cex.axis=0.8,  ...)
}
\arguments{
  \item{x}{
Output of the Sparse-smothed Cox proportional hazard model \code{SSCox}
.}
  \item{Chromosome}{
  	A vector of characters describing the chromosomes to be selected for plotting, for example \code{Chromosome=c("chr1", "chr2").}
}
  \item{xlab}{
Horizontal axis label.
}
  \item{ylab}{
Vertical axis label.
}
  \item{cex.axis}{
Character expansion for axes.
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


# Plot random effects estimates
plot.brandom(fit, type="h")

}

\keyword{ plot }


