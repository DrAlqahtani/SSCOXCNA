plot.brandom <-
function(x, Chromosome=NULL, xlab="Genomic Regions", ylab="", cex.axis=0.8, ...)
{
#
plot.loc <-
function(a,Chr=NULL, xlab="Genomic Regions", ylab="",
cex.axis=0.8,  ...){
# Function to plot the mean of data matrix
# arief at maths.leeds.ac.uk
selected = c("chr1",  "chr2",  "chr3",  "chr4",  "chr5",  "chr6",  "chr7",  "chr8",  "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22") 

if(is.null(Chr)) Chr <- selected

if(is.null(colnames(a))){
   n.a <- names(a)
   } else {
   n.a <- colnames(a)
   }

a <- a[n.a%in%Chr]  
ch <- n.a[n.a%in%Chr]
id <- c(1:c(length(a)))
chro <- unique(ch)
mid.marks <- c()
end.marks <- c()
for(j in chro){
        temp <- id[ch==j]
        mid.marks <- c(mid.marks,floor(mean(range(temp))))
        end.marks <- c(end.marks, max(temp))
      }

   plot(id,a, axes=F, xlab=xlab, ylab=ylab, ...)
   box()
   if(length(end.marks)>1) abline(v=end.marks, lwd=0.5, col="grey70")
   abline(h=seq(0,20,by=0.5),lwd=0.6,col="grey80")
   axis(2, las=1, cex.axis=cex.axis)
   axis(1, at=mid.marks, labels=chro, cex.axis=cex.axis)
}
	
brand <- x$b
names(brand) <- x$Chr
plot.loc(brand,Chr=Chromosome, xlab=xlab, ylab=ylab, cex.axis=cex.axis, ...)
}
