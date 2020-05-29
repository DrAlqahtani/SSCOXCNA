	plot.path.SSCox=function(x,labelsize = 0.6,....){
	
	 betas <- x$betar[,-1]
	 rownames(betas)=1:dim(x$betar)[1]
	 remove <- apply(betas, 1, function(bet) all(bet == 0) )
	 if (all(remove)) stop("all coefficients are zero for all values of lambda in this object")
	 theta=x$range.theta
	  # Adjust the margins to make sure the labels fit
  labwidth <- ifelse(labelsize > 0, max(strwidth(rownames(betas[!remove,]),"inches",labelsize)), 0)
  margins <- par("mai")
  par("mai" = c(margins[1:3], max(margins[4], labwidth*1.4)))

 # Plot
  matplot((1-(x$opt[3]+x$opt[4]))/sqrt(theta), t(betas[!remove,,drop=FALSE]), type ="l", ylab = "coefficient", xlab = expression(w[l]/sqrt(theta)), col=rainbow(sum(!remove)), xlim = rev(range((1-(x$opt[3]+x$opt[4]))/sqrt(theta))))
  
   if (labelsize > 0 && !is.null(rownames(betas))) {
    take <- which(!remove)
    for (i in 1:sum(!remove)) {
      j <- take[i]
      axis(4, at = betas[j,ncol(betas)], labels = rownames(betas)[j],
      las=1,cex.axis=labelsize, col.axis=rainbow(sum(!remove))[i], lty = (i-1) %% 5 + 1, col = rainbow(sum(!remove))[i])
    }
  }
# Reset the margins
  par("mai"=margins)

  return(invisible(NULL))
}
