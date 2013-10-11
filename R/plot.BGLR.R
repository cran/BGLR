plot.BGLR=function(x,...)
{
   ### Check that object is compatible
  if (!inherits(x, "BGLR"))
     stop("This function only works for objects of class `BGLR'");

  if(sum(is.na(x$y)))
  {
    index=is.na(x$y)
    plot(x$y[!index],x$yHat[!index],main="Training")
    
  }else{
    plot(x$y,x$yHat,main="Training")
  }
}
