summary.SSCox<-
function(x)
{
beta <- x$beta
lb <- length(beta)
sd.beta <- sqrt(diag(x$varfix))
t.beta <- beta/sd.beta
pval.t <- pnorm(-abs(t.beta))*2
names.beta <- names(beta)
if(is.null(names.beta)) names.beta <- paste("X",1:length(beta),sep="")
table.result = data.frame(Term=names.beta,
Estimate=round(beta,4), Std.error=round(sd.beta,4),
t.value=round(t.beta,4), p.value= format(pval.t, digits=4),row.names=1)
print(table.result)
}
