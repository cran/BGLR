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


SEXP sample_beta_sparse(SEXP n, SEXP pL, SEXP XL, SEXP xL2, SEXP bL, SEXP e, SEXP varBj, SEXP varE, SEXP minAbsBeta)
{
//    double *xj;
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

/* this is the inner representation of a compressed column major sparse-matrix  
   which is implemented in R (package Matrix) as the class "dgCMatrix" */	
    int *InnerIndices = INTEGER(GET_SLOT(XL, install("i"))),
        *OuterStarts = INTEGER(GET_SLOT(XL, install("p"))),
        *dim = INTEGER(GET_SLOT(XL, install("Dim")));
    double *x = REAL(GET_SLOT(XL, install("x")));

    int row_index,this_column_length;


    pxL2=REAL(xL2);
    pbL=REAL(bL);
    pe=REAL(e);
    pvarBj=REAL(varBj);

    for(j=0; j<cols;j++)
    {

// get the non-zero elements of this column
          this_column_length = OuterStarts[j+1]-OuterStarts[j];

          b=pbL[j];

// sparse dotproduct
          rhs = 0;
          for(int rit=0;rit<this_column_length;rit++) {
          
            row_index = OuterStarts[j]+rit;
            rhs += x[row_index] * pe[InnerIndices[row_index]];
  
          }

          rhs = rhs/sigma2e;
          rhs+=pxL2[j]*b/sigma2e;
  	  c=pxL2[j]/sigma2e + 1.0/pvarBj[j];
	  pbL[j]=rhs/c + sqrt(1.0/c)*norm_rand();

          b-=pbL[j];  

// sparse vector-vector addition
          for(int rit=0;rit<this_column_length;rit++) {
          
            row_index = OuterStarts[j]+rit;
            pe[InnerIndices[row_index]] += x[row_index] * b;
  
          }

          
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

      UNPROTECT(1);

      return(list);
}

/*
  Routine for sampling coefficients for BayesCpi and BayesB
*/

SEXP sample_beta3_sparse(SEXP n, SEXP p, SEXP X, SEXP x2, SEXP b, SEXP d, SEXP error, SEXP varBj, SEXP varE, SEXP minAbsBeta, SEXP probInside)
{
  int i,j,rows,cols;
  double sigma2e, probIn, logOdds,tmp,betaj;
  double logOddsPrior;
  double rhs,c;
  double RSS, Xe, RSS_in, RSS_out;
  double *pX, *perror, *pb, *px2,*pvarBj;
  int inc=1;
  double c1;
  int *pd;
  int change;
  SEXP list;

/* this is the inner representation of a "dgCMatrix"  */	
  int *InnerIndices = INTEGER(GET_SLOT(X, install("i"))),
      *OuterStarts = INTEGER(GET_SLOT(X, install("p"))),
      *dim = INTEGER(GET_SLOT(X, install("Dim")));
  double *x = REAL(GET_SLOT(X, install("x")));

  int row_index,this_column_length;

  double RSS_nonzero;

  cols=INTEGER_VALUE(p);
  rows=INTEGER_VALUE(n);
  sigma2e=NUMERIC_VALUE(varE);
  probIn=NUMERIC_VALUE(probInside);
  logOddsPrior=log(probIn/(1-probIn));

  px2=REAL(x2);
  pd=INTEGER(d);
  pb=REAL(b);
  perror=REAL(error);
  pvarBj=REAL(varBj);

  GetRNGstate();

  c1=0.5/sigma2e;

  RSS=F77_NAME(ddot)(&rows,perror,&inc,perror,&inc);

  for(j=0; j<cols; j++)
  {

          this_column_length = OuterStarts[j+1]-OuterStarts[j];

// sparse dotproduct
// get RSS for the nonzeroes
          Xe = 0;
          RSS_nonzero = 0;
          for(int rit=0;rit<this_column_length;rit++) {
          
            row_index = OuterStarts[j]+rit;
            Xe += x[row_index] * perror[InnerIndices[row_index]];
            RSS_nonzero += perror[InnerIndices[row_index]] * perror[InnerIndices[row_index]];
  
          }



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

// Update RSS only for the elements in perror that change
        RSS -= RSS_nonzero;
        RSS_nonzero = 0;

        if(pd[j]>change)
        {

           betaj=-pb[j];
           Xe = 0;

// Here we update perror, RSS_nonzero and Xe in the same loop

          for(int rit=0;rit<this_column_length;rit++) {
          
            row_index = OuterStarts[j]+rit;
            perror[InnerIndices[row_index]] += x[row_index] * betaj;
            RSS_nonzero += perror[InnerIndices[row_index]] * perror[InnerIndices[row_index]];
            Xe += x[row_index] * perror[InnerIndices[row_index]];
  
          }


        }else{
                betaj=pb[j];

// sparse vector-vector addition
          for(int rit=0;rit<this_column_length;rit++) {
          
            row_index = OuterStarts[j]+rit;
            perror[InnerIndices[row_index]] += x[row_index] * betaj;
            RSS_nonzero += perror[InnerIndices[row_index]] * perror[InnerIndices[row_index]];
  
          }


        }

// Update RSS
           RSS += RSS_nonzero;

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

// update RSS only for the nonzeros
          RSS -= RSS_nonzero;
          RSS_nonzero = 0;
          for(int rit=0;rit<this_column_length;rit++) {
          
            row_index = OuterStarts[j]+rit;
            perror[InnerIndices[row_index]] += x[row_index] * betaj;
            RSS_nonzero += perror[InnerIndices[row_index]] * perror[InnerIndices[row_index]];
  
          }

// Update RSS
           RSS += RSS_nonzero;

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

  UNPROTECT(1);

  return(list);

}

