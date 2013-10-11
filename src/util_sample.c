#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rconfig.h>
#include <R_ext/Lapack.h>


/*
 * This is a generic function to sample betas in various models, including 
 * Bayesian LASSO, BayesA, Bayesian Ridge Regression, BayesCpi, etc.
 
 * For example, in the Bayesian LASSO, we wish to draw samples from the full 
 * conditional distribution of each of the elements in the vector bL. The full conditional 
 * distribution is normal with mean and variance equal to the solution (inverse of the coefficient of the left hand side)
 * of the following equation (See suplementary materials in de los Campos et al., 2009 for details),
   
    (1/varE x_j' x_j + 1/(varE tau_j^2)) bL_j = 1/varE x_j' e
 
    or equivalently, 
    
    mean= (1/varE x_j' e)/ (1/varE x_j' x_j + 1/(varE tau_j^2))
    variance= 1/ (1/varE x_j' x_j + 1/(varE tau_j^2))
    
    xj= the jth column of the incidence matrix
    
 *The notation in the routine is as follows:
 
 n: Number of rows in X
 pL: Number of columns in X
 XL: the matrix X stacked by columns
 XL2: vector with x_j' x_j, j=1,...,p
 bL: vector of regression coefficients
 e: vector with residuals, e=y-yHat, yHat= predicted values
 varBj: vector with variances, 
	For Bayesian LASSO, varBj=tau_j^2 * varE, j=1,...,p
	For Ridge regression, varBj=varB, j=1,...,p, varB is the variance common to all betas.
	For BayesA, varBj=varB_j, j=1,...,p
	For BayesCpi, varBj=varB, j=1,...,p, varB is the variance common to all betas
	
 varE: residual variance
 minAbsBeta: in some cases values of betas near to zero can lead to numerical problems in BL, 
             so, instead of using this tiny number we assingn them minAbsBeta
 
 */


SEXP sample_beta(SEXP n, SEXP pL, SEXP XL, SEXP xL2, SEXP bL, SEXP e, SEXP varBj, SEXP varE, SEXP minAbsBeta)
{
    double *xj;
    double *pXL, *pxL2, *pbL, *pe, *pvarBj;
    double b;
    int inc=1;
    double rhs,c,sigma2e, smallBeta;
    int j,i, rows, cols;

    SEXP list;	

    GetRNGstate();
	
    rows=INTEGER_VALUE(n);
    cols=INTEGER_VALUE(pL);
    sigma2e=NUMERIC_VALUE(varE);
    smallBeta=NUMERIC_VALUE(minAbsBeta);
	
    PROTECT(XL=AS_NUMERIC(XL));
    pXL=NUMERIC_POINTER(XL);

    PROTECT(xL2=AS_NUMERIC(xL2));
    pxL2=NUMERIC_POINTER(xL2);

    PROTECT(bL=AS_NUMERIC(bL));
    pbL=NUMERIC_POINTER(bL);

    PROTECT(e=AS_NUMERIC(e));
    pe=NUMERIC_POINTER(e);

    PROTECT(varBj=AS_NUMERIC(varBj));
    pvarBj=NUMERIC_POINTER(varBj);

    for(j=0; j<cols;j++)
    {
          xj=pXL+j*rows;
          b=pbL[j];
          //F77_NAME(daxpy)(&rows, &b,xj,&inc, pe, &inc);
          rhs=F77_NAME(ddot)(&rows,xj,&inc,pe,&inc)/sigma2e;
          rhs+=pxL2[j]*b/sigma2e;
  	  c=pxL2[j]/sigma2e + 1.0/pvarBj[j];
	  pbL[j]=rhs/c + sqrt(1.0/c)*norm_rand();

          b-=pbL[j];
          //b=-pbL[j];        

          F77_NAME(daxpy)(&rows, &b,xj,&inc, pe,&inc);
          
          if(fabs(pbL[j])<smallBeta)
          {
             pbL[j]=smallBeta;
          }
      }
        
      // Creating a list with 2 vector elements:
      PROTECT(list = allocVector(VECSXP, 2));
      // attaching bL vector to list:
      SET_VECTOR_ELT(list, 0, bL);
      // attaching e vector to list:
      SET_VECTOR_ELT(list, 1, e);

      PutRNGstate();

      UNPROTECT(6);

      return(list);
}

/*
  Routine for sampling coefficients for BayesCpi and BayesB
*/

SEXP sample_beta3(SEXP n, SEXP p, SEXP X, SEXP x2, SEXP b, SEXP d, SEXP error, SEXP varBj, SEXP varE, SEXP minAbsBeta, SEXP probInside)
{
  int i,j,rows,cols;
  double sigma2e, probIn, logOdds,tmp,betaj;
  double logOddsPrior;
  double rhs,c;
  double RSS, Xe, RSS_in, RSS_out;
  double *pX, *perror, *pb, *px2,*pvarBj;
  double *xj;
  int inc=1;
  double c1;
  int *pd;
  int change;
  SEXP list;

  cols=INTEGER_VALUE(p);
  rows=INTEGER_VALUE(n);
  sigma2e=NUMERIC_VALUE(varE);
  probIn=NUMERIC_VALUE(probInside);
  logOddsPrior=log(probIn/(1-probIn));

  PROTECT(X=AS_NUMERIC(X));
  pX=NUMERIC_POINTER(X);

  PROTECT(x2=AS_NUMERIC(x2));
  px2=NUMERIC_POINTER(x2);


  PROTECT(d=AS_INTEGER(d));
  pd=INTEGER_POINTER(d);

  PROTECT(b=AS_NUMERIC(b));
  pb=NUMERIC_POINTER(b);

  PROTECT(error=AS_NUMERIC(error));
  perror=NUMERIC_POINTER(error);

  PROTECT(varBj=AS_NUMERIC(varBj));
  pvarBj=NUMERIC_POINTER(varBj);

  GetRNGstate();

  c1=0.5/sigma2e;

  RSS=F77_NAME(ddot)(&rows,perror,&inc,perror,&inc);
  //RSS=F77_NAME(dnrm2)(&rows,perror,&inc); RSS=RSS*RSS;

  for(j=0; j<cols; j++)
  {
     xj=pX+j*rows;
     Xe=F77_NAME(ddot)(&rows,perror,&inc,xj,&inc);

     if(pd[j])
     {
        RSS_in=RSS;
        RSS_out = RSS_in +  pb[j]*pb[j]*px2[j] + 2*pb[j]*Xe;
     }else{
        RSS_out=RSS;
        RSS_in = RSS_out - pb[j]*pb[j]*px2[j] - 2*pb[j]*Xe;
     }

     logOdds=logOddsPrior+c1*(RSS_out-RSS_in);
     tmp=1.0/(1.0+exp(-logOdds));

     change=pd[j];

     pd[j]=unif_rand()<tmp ? 1 : 0;

     //Update residuals
     if(change!=pd[j])
     {
        if(pd[j]>change)
        {
                betaj=-pb[j];
                F77_NAME(daxpy)(&rows, &betaj,xj,&inc, perror, &inc);
                Xe=F77_NAME(ddot)(&rows,perror,&inc,xj,&inc);
        }else{
                betaj=pb[j];
                F77_NAME(daxpy)(&rows, &betaj,xj,&inc, perror, &inc);
        }
        RSS=F77_NAME(ddot)(&rows,perror,&inc,perror,&inc);
        //RSS=F77_NAME(dnrm2)(&rows,perror,&inc); RSS=RSS*RSS;
     }

     //Sample the coefficients
     if(pd[j]==0)
     {
        //Sample from the prior
        pb[j]=sqrt(pvarBj[j])*norm_rand();
     }else{
	   //Sampling from the conditional
           rhs=(px2[j]*pb[j] + Xe)/sigma2e;
           c=px2[j]/sigma2e + 1.0/pvarBj[j];
           tmp=rhs/c + sqrt(1.0/c)*norm_rand();

           betaj=pb[j]-tmp;
           F77_NAME(daxpy)(&rows, &betaj,xj,&inc, perror, &inc);
           RSS=F77_NAME(ddot)(&rows,perror,&inc,perror,&inc);
           //RSS=F77_NAME(dnrm2)(&rows,perror,&inc); RSS=RSS*RSS;
           pb[j]=tmp;
     }

  }

  // Creating a list with 3 vector elements
  PROTECT(list = allocVector(VECSXP,3));

  // attaching b vector to list
  SET_VECTOR_ELT(list, 0,d);

  // attaching error vector to list
  SET_VECTOR_ELT(list, 1, error);

  // attaching b to the list
  SET_VECTOR_ELT(list,2,b);

  PutRNGstate();

  UNPROTECT(7);

  return(list);

}


