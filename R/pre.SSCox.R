pre.SSCox <-
function(a, l=1, impute=F, plot.it=F, exclude=c("chrX", "chrY", "chrM"))
{
#Function to preprocess the CNA data from the files
attr.C <- attr(a,"Chr")
attr.P <- attr(a,"Pos")
if(!is.null(attr.C)){
 unique.C <- unique(attr.C)
   if(any(exclude%in%unique.C)){
     unselect <- attr.C%in%exclude
     a <- a[,!unselect]
     attr.C <- attr.C[!unselect]
     attr(a,"Chr") <- attr.C
     if(!is.null(attr.P)){
          attr.P <- attr.P[!unselect]
          attr(a,"Pos") <- attr.P
          }# end of if null
   } else {
     warning("No chromosome is excluded.")
    }
} else {
 warning("No chromosome is excluded.")
}

a2 <- check.NA(a, l=l, impute=impute, plot.it=plot.it)

invisible(return(a2))
}

check.NA <-
function(a,l=1, impute=F, plot.it=F){
#Function to check missing values and (optionally) impute them
attr.Chr <- attr(a,"Chr")
attr.Pos <- attr(a,"Pos")
da <- dim(a)
d <- ifelse(da[1] > da[2], 1, 2)
nc <- min(da)
nr <- max(da)
temp <- apply(a,d,function(b) sum(is.na(b)))# no of missing values per variable
temp.pct <- round(temp*100/nc,2) # in percentage missing
seq.missing.vec <- seq(0,100,by=10)
missing <- vector()
for(j in 1:length(seq.missing.vec)){
   missing[j] <- sum(temp.pct>=seq.missing.vec[j])
   }
missing.pct <- round(missing*100/nr,2)
which.l <- which(temp<=l & temp>0)
cat("----- Check NA function -------------------------\n")
 for(j in 1:length(missing)){
 cat("Missing at least",seq.missing.vec[j]," pct =", missing[j],"(",missing.pct[j],"percent out of",nr,")\n")
 }
cat("Number of variables with <= ",l," observations missing is ", length(which.l),"\n")
cat("or ",round(length(which.l)*100/nr,2),"percent\n")
cat("----- Check NA function -------------------------\n")

if(plot.it){
 plot(seq.missing.vec, missing.pct, type="l", xlab="Missing obs. (percent)",
  ylab="Number of variables (percent)")
}

if(impute){
  if(d==1){
    for(k in which.l){
      a[k,is.na(a[k,])] <- mean(a[k,], na.rm=T)
    } # end for k
  } else {
    for(k in which.l){
      a[is.na(a[,k]),k] <- mean(a[,k], na.rm=T)
    } # end for k
  } # end if else d=1
cat("Finished imputation for variables with up to ",l,"missing observations\n")
}#if impute

temp2 <- apply(a,d,function(b) sum(is.na(b)))# no of missing values per variable
selected <- which(temp2==0)
non.selected <- which(temp2>0)
 if(d==1){
   a2 <- a[selected,]
   } else {
   a2 <- a[,selected]
   } # end if d==1

if(!is.null(attr(a,"position"))){
   attr(a2,"position") <- attr(a,"position")
   }# end if null attr
attr(a2,"selected.id") <- selected
attr(a2,"non.selected.id") <- non.selected
attr(a2,"Chr") <- attr.Chr
attr(a2,"Pos") <- attr.Pos
return(invisible(a2))
}

