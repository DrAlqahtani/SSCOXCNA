cv.SScox <-
function(time,event,formula,CNA, fold = 5,alphac.vec,alphan.vec,diff=2,steps,maxiter=400,eps=1e-6,print=TRUE ){
	#if we need include rho we can put it in input rho.vec=0

disx=model.matrix(formula)
x=disx[,-1] # removing the intersebt
disy=model.matrix(~  time+event )
time=disy[,2] # survival time
delta=disy[,3] # delta "cencoring"
data=cbind(time,delta,x)
data=data[order(data[,"time"]),]
delta=data[,"delta"]
z=data[,-c(1,2)]
rho.vec=0
CNAA=cbind(CNA,time)
CNAAF=CNAA[order(CNAA[,"time"]),]
CNAF=CNAAF[,-(ncol(CNAAF))]
if(length(z)==nrow(CNAF)){z=as.matrix(z)}
	
lambdamaxrgeomn=function(CNAF,delta,alphal,bn=matrix(0,nrow=ncol(CNAF)),m){
	
grs=function(bn){
			CNAF=CNAF
			delta=delta
daR=c() # find the dominater formula for each patient 
for(i in 1:nrow(CNAF)){
daR[i]=exp((CNAF[i,]%*%bn))
}
darR=c() #find the dominater formula for all in the risk 
for(j in 1:nrow(CNAF)){
	darR[j]=sum(daR[j:nrow(CNAF)])
	}
	RuR=list() #find the nominater formula for each patient 
for(u in 1:nrow(CNAF)) {
RuR[[u]]=CNAF[u,]*exp((CNAF[u,]%*%bn))
}
RUR=list() # find the nominater formula for all in the risk
for(h in 1:nrow(CNAF)){
	RUR[[h]]=Reduce('+', RuR[h:nrow(CNAF)])
}
	DDR=list() # find the score function 
	for(k in 1:nrow(CNAF)){
DDR[[k]]=delta[k]*(CNAF[k,]-(matrix(unlist(RUR[k]))/darR[k]))
}

scoreRnopenalty=Reduce('+', DDR)
scoreR=scoreRnopenalty
UBR=matrix(scoreR)
return(UBR)}

#lambdamax=grs(bn)/alphal
#lambdamaxsqrt=(lambdamax)^2
#thetamin=1/max(abs(lambdamaxsqrt))
thetamin=((alphal)^2)/((max(abs(grs(bn))))^2)
thetamax=thetamin/(0.01)^2
range.theta=c()
for(j in 0:m){
	range.theta[j+1]=thetamin*((thetamax/thetamin)^(j/m))
}
#if(alphal!=1){thetaminsqrt=thetaminsqrt2}

return(list(thetamaxmin=c(thetamin,thetamax),range.theta=range.theta))}

	
	logPLfixr<- function(time,event,formula,CNA,b,beta){
		disx=model.matrix(formula)
x=disx[,-1] # removing the intersebt
disy=model.matrix(~  time+event )
time=disy[,2] # survival time
delta=disy[,3] # delta "cencoring"
data=cbind(time,delta,x)
data=data[order(data[,"time"]),]
delta=data[,"delta"]
z=data[,-c(1,2)]
CNAA=cbind(CNA,time)
CNAAF=CNAA[order(CNAA[,"time"]),]
CNAF=CNAAF[,-(ncol(CNAAF))]
if(length(z)==nrow(CNAF)){z=as.matrix(z)}

		y=data[,"time"];d=delta;x.fix=z;x=CNAF;b=b;beta=beta
		
                n<-dim(x)[1]
                eta<-c(x.fix%*%beta+x%*%b)
                expeta<-exp(eta)

                ## traingular matrix Mi
                Mi<-round(lower.tri(matrix(1, n, n),diag=TRUE))[,d==1]
                
                like1<-sum(d*eta)
                like2<-sum(log(t(Mi)%*%expeta))
                          
                logPL<-like1-like2 ## partial likelihood
  
  return(logPL) 
  }

	
	n<-dim(CNA)[1]; p<-dim(CNA)[2];m=steps
   
    set.seed(12345) 
    foldi <-split(sample(c(1:72,74:n)), rep(1:fold, length = n-1))
   cvmat <-array(0,c(length(rho.vec),length(alphac.vec),length(alphan.vec),m+1))
   cvmatmat=matrix(nrow=m+1,ncol=fold)
   
    	outcv=matrix(0,ncol=1,nrow=ncol(CNAF))
    	thet.vecall=NULL
    	alphanall=NULL
    	alphacall=NULL
    	for (j in 1:length(rho.vec)){
    		for (k in 1:length(alphac.vec)){
    			for (l in 1:length(alphan.vec)){
    		if (print==TRUE){cat(paste( "alphal =",(1-(alphac.vec[k]+alphan.vec[l])), ", alphac =",alphac.vec[k],", alphan =",alphan.vec[l], "\n")) }
    		 rho<-rho.vec[j];alphac<-alphac.vec[k];alphan<-alphan.vec[l]
    		 if( 1-(alphac+alphan)==0){thet.vec=exp(seq(from=-13,to=-4,length=(m-round(m/3,0))))}
    		if( 1-(alphac+alphan)!=0){thet.vec=lambdamaxrgeomn(CNAF,delta,alphal=(1-(alphac+alphan)),m=m)$range.theta}
    		thet.vecall=cbind(thet.vecall,thet.vec)
    		alphanall=cbind(alphanall,alphan)
    		alphacall=cbind(alphacall,alphac)
    		for (i in 1:(m+1)) {
    			thet<-thet.vec[i]
    			cat(paste( "Theta =", thet.vec[i],"\n")) 
    			objectall <- SSCox(time,event,formula,CNA,init.random=outcv[,i],theta=thet,alphan=alphan,alphac=alphac,diff=diff,maxiter=400)
    			outcv = cbind(outcv,objectall$b)
    		for (ff in 1:fold) {
                omit <- sort(foldi[[ff]])
                
                y.train<-time[-omit]; d.train<-event[-omit];x.train<-CNA[-omit,] ;z.train=x[-omit,]
  object <- SSCox(y.train,d.train,~z.train,x.train,init.random=objectall$b,theta=thet,alphan=alphan,alphac=alphac,maxiter=400); b<-c(object$b);beta<-c(object$beta)

 
 ## partial log-likelihood 
                y.new<-time[omit]; d.new<-event[omit]; x.new<-as.matrix(CNA[omit,]);z.new<-x[omit,]
                 
               
  
  lik<-logPLfixr(y.new,d.new,~z.new,x.new,b,beta)
  cvmatmat[i,ff]= lik
  cvmat[j,k,l,i]<-cvmat[j,k,l,i]+lik
   }
   }
   }
   }
   }
   
   opt.idx<- which(cvmat==max(cvmat), arr.ind=T)
    
    rho.opt <- rho.vec[opt.idx[1]]
    alphac.opt <- alphac.vec[opt.idx[2]]
   alphan.opt <- alphan.vec[opt.idx[3]]
   thet.opt <- as.numeric(thet.vecall[ ,which(alphacall==alphac.opt&alphanall==alphan.opt)][opt.idx[4]])
    cat(paste("\n Optimal parameters for partial likelihood criterion : thet = ", thet.opt, ", ", sep = ""))
    cat(paste("rho = ", rho.opt, "\n", sep = ""))
        cat(paste("alphac = ", alphac.opt, "\n", sep = ""))
    cat(paste("alphan = ", alphan.opt, "\n", sep = ""))


opt<-c(thet.opt,rho.opt,alphac.opt,alphan.opt)
# lasso.idx<-which(cvmat==max(cvmat[,1,]),arr.ind=T)
# thet.opt <- thet.vec[lasso.idx[1]]
 # rho.opt <- rho.vec[lasso.idx[2]]
 # L1.opt<-c(thet.opt,rho.opt)
   sd=apply(cvmatmat,1,sd)/sqrt(fold)
   # sd2=apply(cvmatmat,1,sd)/sqrt(fold)
   range.theta=thet.vec
   cv <- list(cvmat = cvmat, opt = opt, sd=sd ,range.theta= range.theta,betar=outcv)
   invisible(cv)
   }
