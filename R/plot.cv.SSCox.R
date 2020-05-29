	plot.cv.SSCox = function(x,...){
	
	error.bars <-
function(x, upper, lower, width = 0.02, ...)
{
	xlim <- range(x)
	barw <- diff(xlim) * width
	segments(x, upper, x, lower, ...)
	segments(x - barw, upper, x + barw, upper, ...)
	segments(x - barw, lower, x + barw, lower, ...)
	range(upper, lower)
}

  cvobject = x
  xlab = expression(log(theta))
  ylab = "CV partial Log Likelihood"
  cvUp <- cvobject$cvmat + cvobject$sd
  cvDn <- cvobject$cvmat - cvobject$sd
 plot.argument = list(x=log(cvobject$range.theta), y = cvobject$cvmat, xlab = xlab, ylab = ylab, ylim = range(cvUp, cvDn))
  do.call("plot", plot.argument)
  error.bars(log(cvobject$range.theta), cvUp, cvDn, width = 0.01, col ="darkgrey")
  points(log(cvobject$range.theta),cvobject$cvmat,pch=20,col="red")
}
