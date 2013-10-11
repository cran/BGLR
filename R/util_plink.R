
#This function will read a bed file (binary file for genotypes in plink format)
read_bed=function(bed_file,bim_file,fam_file,na.strings=c("0","-9"),verbose=FALSE)
{
	#Extended map file (this gives the number of snps)
	bim=read.table(bim_file, comment.char="", as.is=TRUE, na.strings=na.strings)
	snps=as.character(bim[,2])

	#First 6 columns of ped file (this gives the number of individuals)
	fam=read.table(fam_file, comment.char="", as.is=TRUE, na.strings=na.strings)

	n=nrow(fam)
	p=length(snps)

	out=rep(0,n*p)

        if(verbose)
        {
           verbose=1
        }else 
        {   
           verbose=0
        }

	out=.C("read_bed",as.character(bed_file),as.integer(n),as.integer(p),as.integer(out),as.integer(verbose))[[4]]
        return(list(n=n,p=p,x=out))
}

#This function will read a ped file
#FIXME: It assumes that the missing value is 0
read_ped=function(ped_file)
{
	out=.Call("read_ped",ped_file)
	return(out)
}

#This function will write a bed file (binary file for genotypes in plink format)
write_bed=function(x,n,p,bed_file)
{
   	#Check inputs

   	if(!is.numeric(x)) stop("x should be an integer and numeric\n");
   	if(min(x)<0) stop("Supported codes are 0,1,2,3\n");
   	if(max(x)>3) stop("Supported codes are 0,1,2,3\n");
   	if(n<=0) stop("n should be bigger than 0");
   	if(p<=0) stop("p should be bigger than 0");
   	if(length(x)!=n*p) stop("length of x is not equal to n*p");
          
	#Function call
	.C("write_bed",as.character(bed_file), as.integer(n), as.integer(p), as.integer(x)) 
}
