survpredict.SSCox <-
function(x, new.clin=1, new.cna=NULL, plot = T )
{
lbeta <- length(x$beta)
lb <- length(x$b)
new.clin <- as.matrix(new.clin)
new.cna <- as.matrix(new.cna)
if(nrow(new.clin)!=nrow(new.cna)) stop("The number of observations in the clinical data is not the same as in the CNA profiles.")
if(ncol(new.clin)!=lbeta) stop("Number of columns in the new clinical data is not the same as the one in the fitted results.")
if(ncol(new.cna)!=lb) stop("Number of columns in the new CNA profiles is not the same as the one in the fitted results.")
S0=exp(-x$H0)
surv.est=S0^(exp(new.clin%*%x$beta+new.cna%*%x$b))
 if (plot == T) 
 plot(sort(x$time),surv.est,type="S",ylab="Survival  probability",xlab="Survival time")
return(surv.est=surv.est)
}
