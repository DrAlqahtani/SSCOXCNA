plot.Coxsnell <-
function(x,xlab=NULL, ylab=NULL, ...)
{
# Function to plot AIC as a function of log theta
if(is.null(xlab)) xlab <-"Cox-Snell Residuals"
if(is.null(ylab)) ylab <- "Cumulative hazard"
require(survival)
Cumhazard <- x$H0
S0=exp(-Cumhazard)

Surr=list()
for(i in 1:length(x$time)){
Surr[[i]]=S0^(exp((t(x$X[i,])%*%x$beta)+(t(x$Z[i,])%*%x$b)))
}
Sort=sort(x$time)
act=c()
for(i in 1:length(x$time)){
	act[i]=which(Sort==x$time[i])
}
 surrind=c()
 for(i in 1:length(x$time)){
surrind[i]=Surr[[i]][act[i]]
}
coxsnelresed=-log(surrind)
fitres=survfit(coxph(Surv(coxsnelresed,x$status)~1,method='breslow'),type='aalen')

plot(fitres$time,-log(fitres$surv),type='s',las=1,xlab=xlab, 	ylab=ylab,lwd=2)
abline(0,1,col="gray70",lwd=3,lty=2)

}

