\name{pre.SSCox}
\alias{pre.SSCox}
\alias{check.NA}
\title{
Function to prepare copy number alteration profiles for Sparse-smothed Cox PH model.
}
\description{
The function checks for distribution of missing values in the genomic regions, and optionally impute them. The imputation (if called) will take place only on genomic regions that are missing in a small number of observations (sequences).
}
\usage{
pre.SSCox(a, l = 1, impute = F, plot.it = F, exclude = c("chrX", "chrY", "chrM"))
}
\arguments{
  \item{a}{
A data on CNA profiles as output from function \code{\link{read.norm}.}
}
  \item{l}{
(When \code{impute} is TRUE) Maximum number of missing values per genomic region allowed to be imputed. Keep it low (for example 1 or 2).
}
  \item{impute}{
Should imputation be performed?
}
  \item{plot.it}{
Should a plot be produced?
}
  \item{exclude}{
A vector of characters, indicating chromosomes that are excluded for further analysis. The default is chromosomes X, Y, and M.
}
}
\details{
The function will report the distribution of missing values in the data. If \code{impute=TRUE}, then imputation will be performed on genomic regions (columns of CNA) if the number of missing values in each region/window is maximum \code{l}.
}
\value{
A matrix of CNA profiles, where some chromosomes (if applicable) are excluded and the missing values have been removed.
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
# Generating simulated data
set.seed(1234)
CNA <- matrix(rnorm(1000), 10)
CNA[sample(1:10,1),sample(1:100,10)] <- NA
attr(CNA,"Chr") <- c(paste("chr",rep(1:20,rep(4,20)),sep=""),rep("chrX",20))
attr(CNA,"Pos") <- c(rep(1:4,20),1:20)
#-- Finished generating simulated data


CNA.imputed  <- pre.SSCox(CNA, impute=T)
attr(CNA, "Chr")
attr(CNA.imputed, "Chr")

}
\keyword{missing values }
\keyword{imputation}
