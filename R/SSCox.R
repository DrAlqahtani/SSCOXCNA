SSCox <-
function(time,event,formula,CNA,init.fix=matrix(0,ncol=1,nrow=dim(model.matrix(formula))[2]-1),init.random=matrix(0,ncol=1,nrow=ncol(CNA)),theta,alphan,alphac,diff=2,eps=1e-6,maxiter=500){
NewatenF <- function (time,event,formula,CNA,init.fix=matrix(0,ncol=1,nrow=dim(model.matrix(formula))[2]-1),init.random=matrix(0,ncol=1,nrow=ncol(CNA)),eps,maxiter){
# Use the Newaten-Raphson algorithm to get the MLE

#input are :
#time=the survival time; event=vector of either cencoring or not;
#formula= ~ the explantory variables which need to be in the model
#init.fix= init.fixal value for the estimate, by defult set to be a vector of zeros.
#eps=a value which diside when should the newaten eteration stop 
#maxiter=the maximum number of eteratin;by defult set to be 30
     	
#output are:
# return at each iteration parameters estimates
# the inverse of fisher information (covariance variance matrix of the estimation )

## Design the data and the order disain matrix  	
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

## score function 
UBF=function(z,CNAF,b=init.fix,bn=init.random,delta){
	da=c() # find the dominater formula for each patient 
for(i in 1:nrow(z)){
da[i]=exp((z[i,]%*%b)+(CNAF[i,]%*%bn))
}
dar=c() #find the dominater formula for all in the risk 
for(j in 1:nrow(z)){
	dar[j]=sum(da[j:nrow(z)])
	}
	Ru=list() #find the nominater formula for each patient 
for(u in 1:nrow(z)) {
Ru[[u]]=z[u,]*exp((z[u,]%*%b)+(CNAF[u,]%*%bn))
}
RU=list() # find the nominater formula for all in the risk
for(h in 1:nrow(z)){
	RU[[h]]=Reduce('+', Ru[h:nrow(z)])
}
	DD=list() # find the score function 
	for(k in 1:nrow(z)){
DD[[k]]=delta[k]*(z[k,]-(matrix(unlist(RU[k]))/dar[k]))
}
score=Reduce('+', DD)
UBF=score

return(UBF)}

### fisher information

IBF=function(z,CNAF,b=init.fix,bn=init.random,delta){
library(Matrix)
n=nrow(z)
p=ncol(z)
 I.all=matrix(0,nrow=p,ncol=p)
 for(i in 1:n){
 	I=matrix(0,nrow=p,ncol=p)
ris=matrix(c(rep(0,i-1),rep(1,n-(i-1))),ncol=1) #to find the patient at risk 
w = t(ris)%*%exp(z%*%b+CNAF%*%bn) # scaler 
zbar=Diagonal(n,ris)%*%z
left=w*Diagonal(n,exp(z%*%b+CNAF%*%bn))
right=(exp(z%*%b+CNAF%*%bn))%*%t((exp(z%*%b+CNAF%*%bn)))
I=(delta[i]/(w^2))[1,1]*t(zbar)%*%(left-right)%*%zbar
I.all=I.all+I
}
as.matrix(I.all)}


## Newaten itteration    	
beta=init.fix
betar=init.random
out=init.fix
count=0
continue=T
  while(continue){
		count=count+1
		cat(paste("Count_fix = ",count, "\n", sep = ""))
		beta.old=beta
		betar.old=betar
		dlt=UBF(z,CNAF,beta.old,betar.old,delta)
		dlt2=IBF(z,CNAF,beta.old,betar.old,delta)
		beta=beta.old+(solve(dlt2)%*%(dlt))
		beta=as.matrix(beta)
		out = cbind(out,beta)
		continue = 	(abs(beta-beta.old)>rep(eps,length(beta))) && (count<= maxiter)
		}
		if(count>maxiter){
			warning("Maximum number of itteration reached")
			return(t(out))
		}
		
return(list(coef=t(out),fisher=dlt2,sdrif=-dlt2,cov=solve(dlt2),countf=count)) # return the covereance variance matrix and the itteration with the coeff
	}

# random estimate with convex optmization 
NewtenmixRsqrt=function (time,event,formula,CNA,init.fix=matrix(0,ncol=1,nrow=dim(model.matrix(formula))[2]-1),init.random=matrix(0,ncol=1,nrow=ncol(CNA)),theta,alphan,alphac,diff,eps,maxiter){

.solve <- function(a,b) {
  out <- try(qr.coef(qr(a, LAPACK=TRUE), b))
  if (is(out, "try-error")) stop("Matrix inversion failed. Please increase lambda1 and/or lambda2", call. = FALSE)
  return(out)
}

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
dtimes=time[delta==1]
Riskset=outer(time,dtimes,">=")	
lambda=1/theta
#theta=diag(rep(THET,ncol(CNAF)),ncol=ncol(CNAF),nrow=ncol(CNAF))
theta1=diag(theta,ncol=ncol(CNAF),nrow=ncol(CNAF))
#for(i in 2:ncol(CNAF)){theta[i-1,i]=THET*rho}
#for(i in 2:ncol(CNAF)){theta[i,i-1]=THET*rho}
intial_theta=theta1


if (diff==1){
library(Matrix)
nmat=ncol(CNAF)
bmat=matrix(c(2,-1),nmat,2,byrow=T)
B=bandSparse(nmat,k=c(0,1),diag=bmat,symm=TRUE)
ss=as.matrix(B)
#fix the first few row and coloms
ss[1,1]=1
ss[1,2]=-1
#fixthe lastfew row
ss[nmat,nmat]=1}

if(diff==2){
library(Matrix)
nmat=ncol(CNAF)
bmat=matrix(c(6,-4,1),nmat,3,byrow=T)
B=bandSparse(nmat,k=c(0,1,2),diag=bmat,symm=TRUE)
ss=as.matrix(B)
#fix the first few row and coloms
ss[1,1]=1
ss[1,2]=-2
ss[1,3]=1
ss[2,1]=-2
ss[2,2]=5
#fixthe lastfew row
ss[nmat,nmat]=1
ss[nmat,nmat-1]=-2
 ss[nmat-1,nmat]=-2
 ss[nmat-1,nmat-1]=5}

intial_theta_inv=(1/theta)*ss

#intial for fixed
b=init.fix
brs <- function(bn) {
	nzb <- (bn != 0)
	status=delta
	activeCNAF <- CNAF[,nzb, drop=FALSE]
	linpred <- drop(activeCNAF %*% bn[nzb]+z%*%b)
	 ws=drop(exp( linpred))
	   breslows=drop(1/(ws%*%Riskset))
	breslow=drop(Riskset%*%breslows)
		   #martingle resduals
	 residuals=status-(breslow*ws)
	   #the log liklihood
	   loglik=sum(log(breslows))+sum( linpred[delta==1])
	Pij <- outer(ws, breslows) 
	Pij[!Riskset] <- 0
	W <- list(P = Pij, diagW = breslow * ws)        # construct: W = diag(diagW) - P %*% t(P)
	absb=c()
 bsq=c()
 for(i in 1:ncol(CNAF)){
	absb[i]=abs(bn[i])
	 bsq[i]=(bn[i])^2
 }
 penlty=(sqrt(lambda)*(1-(alphan+alphac))*sum(absb))+(0.5*(lambda*alphan)*sum(bsq))+((alphac)*((ncol(CNAF)+1)/2)*(log(1+(t(bn)%*%intial_theta_inv%*%bn))))
 return(list(ll=loglik,penalty=penlty,residuals = residuals,W=W))
		 }
		 
		 grs=function(bn){
			CNAF=CNAF
			delta=delta
daR=c() # find the dominater formula for each patient 
for(i in 1:nrow(CNAF)){
daR[i]=exp((CNAF[i,]%*%bn)+(z[i,]%*%b))
}
darR=c() #find the dominater formula for all in the risk 
for(j in 1:nrow(CNAF)){
	darR[j]=sum(daR[j:nrow(CNAF)])
	}
	RuR=list() #find the nominater formula for each patient 
for(u in 1:nrow(CNAF)) {
RuR[[u]]=CNAF[u,]*exp((CNAF[u,]%*%bn)+(z[u,]%*%b))
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
#fdnor=(solve(intial_theta)%*%bn)
fdnor=(1/intial_theta[1,1])*bn
fdcaush=(((ncol(CNAF)+1)*intial_theta_inv%*%bn)/as.numeric((1+(t(bn)%*%intial_theta_inv%*%bn))))
scoreR=scoreRnopenalty-((alphan*fdnor)+(alphac*fdcaush))
scorr=c()
for(i in 1:ncol(CNAF)){
if(sign(bn[i])!=0){
scorr[i]=scoreR[i]-((1 - (alphan + alphac))*sqrt(lambda)*sign(bn[i]))
}
if(sign(bn[i])==0 && abs(scoreR[i])> ((1-(alphan+alphac))*sqrt(lambda))){
	scorr[i]=scoreR[i]-((1 - (alphan + alphac))*sqrt(lambda)*sign(scoreR[i]))
}
if(sign(bn[i])==0 && abs(scoreR[i])<= ((1-(alphan+alphac))*sqrt(lambda))){
	scorr[i]=0
}
}
UBR=matrix(scorr)

return(UBR)}

# second directed derivetive 
IBRs=function(bn){
	 brfit=brs(bn)
	 dir=grs(bn)
	nzb <- (bn != 0)
	 newb=(dir!=0)
	 active= nzb| newb
	 CNAFdir <- drop(CNAF[,active,drop=F] %*%dir[active])
	 if (is.list(brfit$W)) {
	curve <- (sum(CNAFdir * CNAFdir * brfit$W$diagW) - drop(crossprod(crossprod(brfit$W$P, CNAFdir))))
	}else if (length(brfit$W) > 1) {
            curve <- sum(CNAFdir * CNAFdir * brfit$W)
            }else {
            curve <- sum(CNAFdir *CNAFdir)
            }
 caus1=ncol(CNAF)+1
caus2=as.numeric(1+(t(bn)%*%intial_theta_inv%*%bn))
caus3=intial_theta_inv
caus4=2*(ncol(CNAF)+1)
caus5=(intial_theta_inv%*%bn)%*%(t(bn)%*%intial_theta_inv)
causd=as.numeric((1+(t(bn)%*%intial_theta_inv%*%bn))^2)
sdcaus=((caus1*caus2*caus3)-(caus4*caus5))/causd

I.allF=curve+sum((alphan/theta)*dir*dir)+(alphac)*((t(dir)%*%sdcaus%*%dir))
return(list(I.allF=I.allF,sdcaus=sdcaus))}
		 
		 betar=init.random
		 beta=init.fix
		 out=init.random
		  m <- length(betar)
  n <- nrow(CNAF)
 LL <- -Inf
  penalty <- penalty1 <- Inf
  active <- !logical(m)
  nvar <- m
count=0
NR=FALSE
#continue=T
 finished <- FALSE
  newfit <- TRUE
  retain <- 0.05
  cumsteps <- 0
   while(!finished ){
  	 nzb <- (betar!= 0)
  	 dlt=grs(betar)
  	oldLL <- LL
  	   oldpenalty <- penalty
  	   LL <-  brs(betar)$ll
  	   penalty<-brs(betar)$penalty
  	   finishedLL <- (2 * abs(LL - oldLL) / (2 * abs(LL - penalty) + 0.1) < eps)
		  finishedpen <- (2 * abs(penalty - oldpenalty) / (2 * abs(LL - penalty) + 0.1) < eps)
cumsteps <- 0
        newb= (dlt!=0)
  	 oldactive=active
		 active= nzb| newb
       oldnvar <- nvar
       nvar=sum(active)
       finishednvar <- !any(xor(active, oldactive))
        finished <- (finishedLL && finishedpen && finishednvar&NR) || (all(dlt == 0)) || (count == maxiter)	
         count=count+1
         cat(paste("count_random = ",count, "\n", sep = ""))
 dlt22=IBRs(betar)
		dlt2=dlt22$I.allF
		if(dlt2==0){betar=matrix(rep(0,dim(dlt)[1]))}
		else{
			topt=sum(dlt*dlt)/dlt2 
		 tedge <- numeric(ncol(CNAF))
		  tedge[active] <- -betar[active] / dlt[active]
		   tedge[tedge <= 0] <- 2 * topt
		    mintedge <- min(tedge)
   if(count>10){
		    if (is.list(brs(betar)$W)){
		    	activeCNAF=CNAF[,active,drop=FALSE]
		    	Xdw<- 	activeCNAF*matrix(sqrt( brs(betar)$W$diagW),nrow( activeCNAF),ncol(activeCNAF))
		    	hessian <- -crossprod(Xdw)+crossprod(crossprod(brs(betar)$W$P,activeCNAF))
		    	if((alphan+alphac)!=0){hessian=hessian-(diag(rep(alphan*lambda,sum(active)),nrow=sum(active),ncol=sum(active))+(alphac*(dlt22$sdcaus)[active,active]))}
		    }
		    NRbeta=betar[active]-drop(.solve(hessian,dlt[active]))
		   NR= as.numeric(topt)<mintedge & all(sign(NRbeta)==sign(betar[active]) )
		    #cat ("NR",NR,"\n")
		    if(topt<mintedge & all(sign(NRbeta)==sign(betar[active]) )){
		    		betar[active]=NRbeta
		    }
		    }
		    
		     if (topt>=mintedge){
		    	betar=betar+(mintedge*(dlt))
		    }
		
		if(count>10){
			if(topt<mintedge & !all(sign(NRbeta)==sign(betar[active]))){
		    	 	betar=betar+(as.numeric(topt)*(dlt))
		    }
         }
         
         if(count<=10){
		    	if(topt<mintedge){
		    		betar=betar+(as.numeric(topt)*(dlt))
		    		}
		    		}
		    		}
		#betar=betar+(min(topt,mintedge)*dlt)
  	 betar=as.matrix(betar)
		out = cbind(out,betar)
	}
		if(count==maxiter){
			return(list(betar=t(out)[dim(t(out))[1],],numnonzero=nvar,loglik=LL,penalty=penalty,iterations=count,converged=finished))
		}
return(list(betar=t(out)[dim(t(out))[1],],numnonzero=nvar,loglik=LL,penalty=penalty,iterations=count,converged=finished))
		}

theta=theta
cat(paste("Theta = ", theta, "\n", sep = ""))
KKF=list()
KKR=list()
rand=list()
fix=list()
fisher=list()
fisherf=list()
second=list()
cov=list()
covr=list()
KKF[[1]]=NewatenF(time,event,formula,CNA,eps=eps,maxiter=maxiter)
KKR[[1]]=NewtenmixRsqrt(time,event,formula,CNA,init.fix=KKF[[1]]$coef[dim(KKF[[1]]$coef)[1],],theta=theta,alphac=alphac,alphan=alphan,diff=diff,maxiter=maxiter,eps=eps)
rand[[1]]=KKR[[1]]$betar
fix[[1]]=KKF[[1]]$coef[dim(KKF[[1]]$coef)[1],]
fisherf[[1]]=KKF[[1]]$fisher
 cov[[1]]=KKF[[1]]$cov
 nvar=KKR[[1]]$numnonzero
 ll=KKR[[1]]$loglik
iterations=KKR[[1]]$iterations
 converged=KKR[[1]]$converged
beta=fix[[length(fix)]]
b=rand[[length(rand)]]	
fisherfix=fisherf[[length(fisherf)]]
covfix=cov[[length(cov)]]
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

n=nrow(CNAF) 
haz_base=c()
for(j in 1:n){
 		haz_base[j]=exp(t(beta)%*%z[j,]+(t(b)%*%CNAF[j,])) 
 		}
 		haz_base1=c()
    for(kj in 1:n){
    	haz_base1[kj]=delta[kj]/sum(haz_base[kj:n])
    }
    Cum_haz=c()
    for(f in 1:n){
    Cum_haz[f]=sum(haz_base1[1:f])
    }
 H_0=Cum_haz

if (ncol(disx)!=2){
j=seq(1:nrow(CNAF))
Risk=matrix(nrow=nrow(CNAF),ncol=1)
Risk[j,]=exp((z[j,]%*%beta)+(CNAF[j,]%*%b))

#Riskno=matrix(nrow=nrow(CNAF),ncol=1)
#Riskno[j,]=exp((x[j,]%*%beta)+(CNA[j,]%*%b))
# loob to find PPL

H=c()
for (i in 1: nrow(CNAF)) {
	H1=delta[i]*(z[i,]%*%beta+CNAF[i,]%*%b)
	H2=delta[i]*log(sum(Risk[i:nrow(CNAF),]))
	H[i]=H1-H2
}
}
cat(paste("# of non-zero coefficients= ", nvar, "\n", sep = ""))
return(list(beta=beta,b=b,nvar=nvar,ll=ll,iterations=iterations,converged=converged,varfix=covfix,H0=H_0,time=time,status=delta,X=z,Z=CNAF,thetavalue=theta,Chr=attr(CNA,"Chr"),Pos=attr(CNA,"Pos")))
}