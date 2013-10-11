predict.BGLR=function(object,newdata = NULL,...)
{
   ### Check that object is compatible
  if (!inherits(object, "BGLR"))
     stop("This function only works for objects of class `BGLR'");
}


